#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"

#include "khash.h"
#include "kvec.h"
#include "path.h"

#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 65536)

typedef struct {
  char *chrom;
  int s;
  int e;
} region_t;

int parse_region(char *region, char *chr, int *s, int *e) {
  int i;
  char *p, *q;
  for (i = 0, p = q = region;; ++p) {
    if (*p == ':' || *p == '-' || *p == 0) {
      int c = *p;
      *p = 0;
      if (i == 0 && c == ':') {
        strncpy(chr, q, p - q);
        chr[p - q] = '\0';
      } else if (i == 1 && c == '-') {
        *s = atoi(q);
      } else if (i == 2 && c == 0) {
        *e = atoi(q);
        return 0;
      } else {
        return 1;
      }
      ++i;
      q = p + 1;
    }
  }
}

int main(int argc, char *argv[]) {
  double rt0 = realtime(), rt = rt0;

  char *gfa_fn = argv[1];
  char *region = argv[2];
  /* char *bed_fn = argv[2]; */

  /* /\* Parsing BED file *\/ */
  /* char line[256]; */
  /* int i; */
  /* char *p, *q; */
  /* int regions_c = 64; */
  /* int regions_n = 0; */
  /* region_t *regions = malloc(regions_c * sizeof(region_t)); */
  /* FILE *file = fopen(bed_fn, "r"); */
  /* if (file == NULL) { */
  /*   fprintf(stderr, "Unable to open file\n"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  /* while (fgets(line, sizeof(line), file)) { */
  /*   if (regions_n == regions_c) { */
  /*     regions = realloc(regions, 2 * regions_c * sizeof(region_t)); */
  /*     regions_c *= 2; */
  /*   } */
  /*   for (i = 0, p = q = line;; ++p) { */
  /*     if (*p == '\t' || *p == 0) { */
  /*       *p = 0; */
  /*       if (i == 0) { */
  /*         regions[regions_n].chrom = malloc(p - q + 1); */
  /*         strncpy(regions[regions_n].chrom, q, p - q); */
  /*         regions[regions_n].chrom[p - q] = '\0'; */
  /*       } else if (i == 1) { */
  /*         regions[regions_n].s = atoi(q); */
  /*       } else if (i == 2) { */
  /*         regions[regions_n].e = atoi(q); */
  /*         break; */
  /*       } else { */
  /*         fprintf(stderr, "Error while reading BED file\n"); */
  /*         exit(EXIT_FAILURE); */
  /*       } */
  /*       ++i; */
  /*       q = p + 1; */
  /*     } */
  /*   } */
  /*   ++regions_n; */
  /* } */
  /* fclose(file); */

  /* /\* XXX: assuming sorted bed *\/ */
  /* char **paths = malloc(regions_n * sizeof(char *)); */
  /* int paths_n = 0; */
  /* for (int i = 0; i < regions_n; ++i) { */
  /*   if (i == 0 || strcmp(paths[i - 1], regions[i].chrom) != 0) { */
  /*     paths[paths_n] = malloc(strlen(regions[i].chrom) + 1); */
  /*     strcpy(paths[paths_n++], regions[i].chrom); */
  /*   } */
  /* } */
  /* for (int i = 0; i < paths_n; ++i) { */
  /*   printf("%s ", paths[i]); */
  /* } */
  /* printf("\n"); */

  graph_t *graph = init_graph(gfa_fn);
  load_vertices(graph, 1);
  fprintf(stderr, "loaded %d vertices (%ld total size) in %.3f secs\n",
          graph->vertices->n, graph->vertices->text_n, realtime() - rt);

  char chrom[32];
  int sp, ep;
  parse_region(region, chrom, &sp, &ep);

  /* Get reference path */
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(gfa_fn, "r");
  if (fp == 0) {
    exit(1);
  }
  kstream_t *ks = ks_init(fp);
  path_t *path = init_path();
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, graph, path);
      else
        gfa_parse_W(s.s, graph, path);

      if (strcmp(path->idx, chrom) == 0) {
        break;
      }
      clear_path(path);
    }
  }
  ks_destroy(ks);
  gzclose(fp);
  /* */

  if (kv_size(path->vertices) == 0) {
    fprintf(stderr, "Path %s not found\n", chrom);
    exit(EXIT_FAILURE);
  }

  /* Get source/sink vertices by iterating over the path */
  rt = realtime();
  int offset = -1;
  uint source = -1, sink = -1;
  int v;
  for (v = 0; v < kv_size(path->vertices); ++v) {
    int iv = p_get_v(path, v); /* graph space id */
    offset += l2b_get_l(graph->vertices, iv);
    if (source == -1 && offset >= sp)
      source = iv;
    if (source != -1 && sink == -1 && offset >= ep) {
      sink = iv;
      break;
    }
  }
  if (source == -1 && sink == -1 /* v == kv_size(path->vertices) */) {
    fprintf(stderr, "Something bad happened\n");
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "found source (%d) and sink (%d) in %.3f secs\n", source,
          sink, realtime() - rt);

  /* */
  ks = ks_init(fp);
  clear_path(path);
  fp = gzopen(gfa_fn, "r");

  khiter_t k1, k2;
  int n1, n2, p1, p2;
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      if (s.s[0] == 'P')
        gfa_parse_P(s.s, graph, path);
      else
        gfa_parse_W(s.s, graph, path);

      k1 = kh_get(im, path->data, source << 4);
      k2 = kh_get(im, path->data, sink << 4);
      if (k1 != kh_end(path->data) && k2 != kh_end(path->data)) {
        n1 = kh_value(path->data, k1);
        k1 = kh_get(im, path->data, (source << 4) | n1);
        p1 = kh_value(path->data, k1);

        n2 = kh_value(path->data, k2);
        k2 = kh_get(im, path->data, (sink << 4) | n2);
        p2 = kh_value(path->data, k2);

        printf("%d %d\n", p1, p2);
        for (int i = p1; i <= p2; ++i) {
          printf("%d,", kv_A(path->vertices, i) >> 1);
        }
        printf("\n");
      }
      clear_path(path);
    }
  }
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  /* */

  destroy_path(path);
  destroy_graph(graph);

  /* for (int i = 0; i < regions_n; ++i) */
  /*   free(regions[i].chrom); */
  /* free(regions); */

  fprintf(stderr, "done in in %.3f secs\n", realtime() - rt0);

  return 0;
}
