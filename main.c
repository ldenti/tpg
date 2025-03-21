#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khash.h"
#include "kseq.h"
#include "kvec.h"

#include "graph.h"
#include "misc.h"
#include "path.h"

KSTREAM_INIT(gzFile, gzread, 65536)
KHASH_MAP_INIT_STR(sti, int)
KHASH_SET_INIT_INT64(set)

typedef struct {
  char *chrom;
  int s;
  int e;
} region_t;

int main(int argc, char *argv[]) {
  double rt0 = realtime(), rt = rt0;

  char *gfa_fn = argv[1];
  char *bed_fn = argv[2]; /* XXX: assuming sorted bed */

  /* Parse BED file */
  char line[256];
  int i;
  char *p, *q;
  int regions_c = 64;
  int regions_n = 0;
  region_t *regions = malloc(regions_c * sizeof(region_t)); /* Regions */

  khiter_t k;
  int hret;
  khash_t(sti) *paths = kh_init(sti); /* Prefix sum of paths */

  FILE *file = fopen(bed_fn, "r");
  if (file == NULL) {
    fprintf(stderr, "Unable to open file\n");
    exit(EXIT_FAILURE);
  }
  while (fgets(line, sizeof(line), file)) {
    if (regions_n == regions_c) {
      regions = realloc(regions, 2 * regions_c * sizeof(region_t));
      regions_c *= 2;
    }
    for (i = 0, p = q = line;; ++p) {
      if (*p == '\t' || *p == 0) {
        *p = 0;
        if (i == 0) {
          regions[regions_n].chrom = malloc(p - q + 1);
          strncpy(regions[regions_n].chrom, q, p - q);
          regions[regions_n].chrom[p - q] = '\0';
          k = kh_put(sti, paths, regions[regions_n].chrom, &hret);
          if (hret == 1)
            kh_value(paths, k) = regions_n;
        } else if (i == 1) {
          regions[regions_n].s = atoi(q);
        } else if (i == 2) {
          regions[regions_n].e = atoi(q);
          break;
        } else {
          fprintf(stderr, "Error while reading BED file\n");
          exit(EXIT_FAILURE);
        }
        ++i;
        q = p + 1;
      }
    }
    ++regions_n;
  }
  fclose(file);

  /* Load graph (segments) */
  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph, 1);
  fprintf(stderr, "loaded %d vertices (%ld total size) in %.3f secs\n",
          graph->vertices->n, graph->vertices->text_n, realtime() - rt);

  /* Exract source/sink for each region */
  rt = realtime();
  region_t *vertices =
      malloc(regions_n * sizeof(region_t)); /* Regions as vertex pair */
  int vertices_n = 0;
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0) {
    exit(1);
  }
  kstream_t *ks = ks_init(fp);
  path_t *path = ph_init();
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, graph, path);
      else
        gfa_parse_W(s.s, graph, path);
      k = kh_get(sti, paths, path->idx);
      if (k != kh_end(paths)) {
        int p = kh_value(paths, k);
        int sp, ep;
        while (p < regions_n && strcmp(path->idx, regions[p].chrom) == 0) {
          sp = regions[p].s;
          ep = regions[p].e;
          int offset = -1;
          uint source = -1, sink = -1;
          int v;
          for (v = 0; v < kv_size(path->vertices); ++v) {
            int iv = ph_getv(path, v); /* graph space id */
            offset += sgms_get_l(graph->vertices, iv);
            if (source == -1 && offset >= sp)
              source = iv;
            if (source != -1 && sink == -1 && offset >= ep) {
              sink = iv;
              break;
            }
          }
          if (source != -1 && sink != -1 /* v < kv_size(path->vertices) */) {
            vertices[vertices_n] =
                (region_t){NULL, source, sink}; /* TODO: store ID */
            ++vertices_n;
            fprintf(stderr, "hit: %s %d %d\n", regions[p].chrom, sp, ep);
          }
          ++p;
        }
      }
      ph_clear(path);
    }
  }
  ks_destroy(ks);
  gzclose(fp);
  fprintf(stderr, "found %d/%d sources/sinks in %.3f secs\n", vertices_n,
          regions_n, realtime() - rt);
  /* */

  /* Reiterate over paths to get subhaplotypes */
  rt = realtime();
  khash_t(set) *VV = kh_init(set);
  khash_t(set) *EE = kh_init(set);

  int subhaps_n = 0;
  int subhaps_c = 128;
  uint32_t **subhaps = malloc(subhaps_c * sizeof(uint32_t *));

  fp = gzopen(gfa_fn, "r");
  ks = ks_init(fp);
  ph_clear(path);
  khiter_t k1, k2;
  int n1, n2, p1, p2;
  int source, sink;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, graph, path);
      else
        gfa_parse_W(s.s, graph, path);

      for (int vv = 0; vv < vertices_n; ++vv) {
        source = vertices[vv].s;
        sink = vertices[vv].e;
        k1 = kh_get(im, path->data, source << 4);
        k2 = kh_get(im, path->data, sink << 4);
        if (k1 != kh_end(path->data) && k2 != kh_end(path->data)) {
          n1 = kh_value(path->data, k1);
          k1 = kh_get(im, path->data, (source << 4) | n1);
          p1 = kh_value(path->data, k1);

          n2 = kh_value(path->data, k2);
          k2 = kh_get(im, path->data, (sink << 4) | n2);
          p2 = kh_value(path->data, k2);

          if (p1 > p2)
            continue;

          /* Store vertices and edges as sets */
          for (int i = p1; i <= p2; ++i) {
            kh_put(set, VV, kv_A(path->vertices, i) >> 1, &hret);
            if (i < p2) {
              kh_put(set, EE,
                     ((uint64_t)(kv_A(path->vertices, i) >> 1) << 32) |
                         (kv_A(path->vertices, i + 1) >> 1),
                     &hret);
            }
          }
          /* Copy subhaplotypes */
          if (subhaps_n == subhaps_c) {
            subhaps = realloc(subhaps, subhaps_c * 2 * sizeof(uint32_t *));
            subhaps_c *= 2;
          }
          subhaps[subhaps_n] =
              malloc((p2 - p1 + 2) * sizeof(uint)); /* +1 for size */
          subhaps[subhaps_n][0] = p2 - p1 + 1;
          memcpy(subhaps[subhaps_n] + 1, path->vertices.a + p1,
                 subhaps[subhaps_n][0] * sizeof(uint32_t));
          ++subhaps_n;
        }
      }
      ph_clear(path);
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  fprintf(stderr, "extracted subhaplotypes in %.3f secs\n", realtime() - rt);

  /* */

  /* Print GFA */
  uint64_t v1, v2, ee;
  int label_c = 128;
  char *label = malloc(label_c);
  for (k = kh_begin(VV); k != kh_end(VV); ++k) {
    if (kh_exist(VV, k)) {
      v1 = kh_key(VV, k);
      sgms_get(graph->vertices, v1, &label, &label_c);
      printf("S\t%d\t%s\n", graph->vertices->vnames[v1], label);
    }
  }
  free(label);

  for (k = kh_begin(EE); k != kh_end(EE); ++k) {
    if (kh_exist(EE, k)) {
      ee = kh_key(EE, k);
      v1 = ee >> 32;
      v2 = (uint32_t)ee;
      printf("L\t%d\t+\t%d\t+\t0M\n", graph->vertices->vnames[v1],
             graph->vertices->vnames[v2]);
    }
  }

  for (int sh = 0; sh < subhaps_n; ++sh) {
    /* TODO: dump meaningful path name as W line */
    printf("P\tSH.%d.%d\t", sh, subhaps[sh][0]);
    /* we are sure to have at least one vertex */
    printf("%d%c", graph->vertices->vnames[subhaps[sh][1] >> 1],
           (subhaps[sh][1] & 1) == 1 ? '+' : '-');
    for (int i = 2; i < subhaps[sh][0] + 1; ++i) {
      printf(",%d%c", graph->vertices->vnames[subhaps[sh][i] >> 1],
             (subhaps[sh][i] & 1) == 1 ? '+' : '-');
    }
    printf("\t*\n");
  }
  /* */

  kh_destroy(set, VV);
  kh_destroy(set, EE);
  for (int i = 0; i < subhaps_n; ++i)
    free(subhaps[i]);
  free(subhaps);
  ph_destroy(path);
  destroy_graph(graph);
  kh_destroy(sti, paths);
  for (int i = 0; i < regions_n; ++i)
    free(regions[i].chrom);
  free(regions);
  free(vertices);

  fprintf(stderr, "done in in %.3f secs\n", realtime() - rt0);

  return 0;
}
