#include "graph.h"

KSTREAM_INIT(gzFile, gzread, 65536)

graph_t *init_graph(char *fn) {
  graph_t *g = malloc(sizeof(graph_t));
  g->fn = fn;

  g->vertices = sgms_init();
  g->v_map = kh_init(im);

  return g;
}

void destroy_graph(graph_t *g) {
  sgms_destroy(g->vertices);
  kh_destroy(im, g->v_map);
  free(g);
}

int load_vertices(graph_t *g, int wseq) {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(g->fn, "r");
  if (fp == 0)
    return 0;
  kstream_t *ks = ks_init(fp);
  int hret;
  khint_t k;
  seg_t *seg = init_seg();
  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'S') {
      gfa_parse_S(s.s, seg, wseq);

      k = kh_put(im, g->v_map, seg->idx, &hret);
      kh_value(g->v_map, k) = g->vertices->n;

      sgms_add(g->vertices, seg->idx, seg->seq, seg->l);
    }
  }
  sgms_add(g->vertices, -1, "\0", 1);

  destroy_seg(seg);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);

  return 0;
}

int get_iidx(graph_t *g, int v) {
  return kh_value(g->v_map, kh_get(im, g->v_map, v));
}

/** === GFA PARSING ======= **/
void gfa_parse_S(char *s, seg_t *ret, int wseq) {
  int i;       // , is_ok = 0;
  char *p, *q; // *seg = 0, *seq = 0, *rest = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        ret->idx = atoi(q);
        // strcpy(ret->idx, q);
      } else if (i == 1) {
        ret->l = p - q;
        if (wseq == 1) {
          if (ret->l >= ret->c) {
            char *temp = realloc(ret->seq, (ret->l * 2) * sizeof(char));
            if (temp == NULL) {
              free(ret->seq);
              fprintf(stderr,
                      "Error while reallocating memory for segment %d\n",
                      ret->idx);
              exit(2);
            } else {
              ret->seq = temp;
            }
            ret->c = ret->l * 2;
          }
          strcpy(ret->seq, q);
        }
        // is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0)
        break;
    }
  }
  // if (!is_ok) { // something is missing
}

void gfa_parse_P(char *s, graph_t *g, path_t *path) {
  int i, v;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i == 0) {
        /* TODO: check for duplicates */
        strcpy(path->idx, q); /* TODO: reallocation */
      } else if (i == 1) {
        char strand = *(p - 1);
        qq = q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == ',') {
            int c = *qq;
            *qq = 0;
            strand = *(qq - 1);
            *(qq - 1) = 0;
            v = get_iidx(g, atoi(q));
            ph_addv(path, v, strand == '+');
            q = qq + 1;
            if (c == 0)
              break;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}

void gfa_parse_W(char *s, graph_t *g, path_t *path) {
  int i, v;
  char *p, *q, *qq;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      *p = 0;
      if (i < 2) {
        *p = '#';
        ++i;
        continue;
      } else if (i == 2) {
        strcpy(path->idx, q); /* TODO: reallocation */
      } else if (i == 5) {
        char strand = *q;
        ++q;
        for (qq = q;; ++qq) {
          if (*qq == 0 || *qq == '>' || *qq == '<') {
            int c = *qq;
            *qq = 0;
            v = get_iidx(g, atoi(q));
            ph_addv(path, v, strand == '>');
            q = qq + 1;
            if (c == 0)
              break;
            strand = c;
          }
        }
        break;
      }
      ++i, q = p + 1;
    }
  }
}
