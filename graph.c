#include "graph.h"

KSTREAM_INIT(gzFile, gzread, 65536)

graph_t *init_graph(char *fn) {
  graph_t *g = malloc(sizeof(graph_t));
  g->fn = fn;

  /* vertices */
  g->vertices = l2b_init();
  g->v_map = kh_init(im);

  g->cp = 512;
  /* g->paths = ph_init(); */

  return g;
}

void destroy_graph(graph_t *g) {
  l2b_destroy(g->vertices);
  kh_destroy(im, g->v_map);
  /* ph_destroy(g->paths); */
  free(g);
}

void print_size(graph_t *g) {
  /* float size = 0.0; */
  /* int n = 0; */

  /* /\* labels *\/ */
  /* size += g->vertices->text_n / 1024 / 1024; */
  /* /\* vertices ids *\/ */
  /* size += g->vertices->n * 4 / 1024 / 1024; */
  /* /\* offsets *\/ */
  /* size += g->vertices->n * 8 / 1024 / 1024; */
  /* /\* total elements *\/ */
  /* n = g->vertices->n; */
  /* printf("V %fMB %d\n", size, n); */

  /* size = 0.0; */
  /* n = 0; */

  /* printf("%fMB\n", (float)g->paths->np * 16 / 1024 / 1024); */
  /* printf("%fMB\n", (float)kv_size(g->paths->vertices) * 4 / 1024 / 1024); */
  /* /\* for (int i = 0; i < g->paths->np - 1; ++i) { *\/ */
  /* /\*   size += g->paths->data[i]->n_buckets * 2 * 8 / 1024 / 1024; *\/ */
  /* /\*   n += g->paths->data[i]->size; *\/ */
  /* /\* } *\/ */
  /* /\* printf("P %fMB %d\n", size, n); *\/ */
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
      l2b_add(g->vertices, seg->idx, seg->seq, seg->l);
      k = kh_put(im, g->v_map, seg->idx, &hret);
      kh_value(g->v_map, k) =
          g->vertices->n - 1; /* -1 since we already added current vertex */
    }
  }
  l2b_add(g->vertices, -1, "\0", 1);
  destroy_seg(seg);
  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
}

int get_iidx(graph_t *g, int v) {
  return kh_value(g->v_map, kh_get(im, g->v_map, v));
}

/* XXX: refactor this function and move parts to path */
/* void compress_path(graph_t *g, paths_t *paths) { */
/*   int pi = paths->np - 1; */
/*   int pl = ph_get_l(g->paths, pi); */

/*   /\* Get sorted array of "used" vertices *\/ */
/*   int *keys = malloc(pl * sizeof(int)); */
/*   int nv = 0; */
/*   khiter_t k; */
/*   for (k = kh_begin(g->data); k != kh_end(g->data); ++k) */
/*     if (kh_exist(g->data, k)) */
/*       keys[nv++] = kh_key(g->data, k); */
/*   ks_introsort(int, nv, keys); */

/*   /\* Build run-length encoded bit vector for path *\/ */
/*   paths->ropes[pi] = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN); */
/*   int p = 0; */
/*   if (keys[0] != 0) { */
/*     /\* printf("Inserting %d 0s\n", keys[0]); *\/ */
/*     rope_insert_run(paths->ropes[pi], 0, 0, keys[0], NULL); */
/*     p = keys[0]; */
/*   } */
/*   int rl = 1; */
/*   for (int i = 1; i < nv; ++i) { */
/*     if (keys[i] == keys[i - 1] + 1) { */
/*       ++rl; */
/*     } else { */
/*       /\* printf("Inserting %d 1s\n", rl); *\/ */
/*       rope_insert_run(paths->ropes[pi], p, 1, rl, NULL); */
/*       p += rl; */
/*       /\* printf("Inserting %d 0s\n", keys[i] - p); *\/ */
/*       rope_insert_run(paths->ropes[pi], p, 0, keys[i] - p, NULL); */
/*       p = keys[i]; */
/*       rl = 1; */
/*     } */
/*   } */
/*   /\* printf("Inserting %d 1s\n", rl); *\/ */
/*   rope_insert_run(paths->ropes[pi], p, 1, rl, NULL); */

/*   /\* Get offsets for each vertex using its repetitions. Replace values in
 * data */
/*    * hash table *\/ */
/*   int v = 0; */
/*   int ofx = 0; */
/*   int tmp = 0; */
/*   for (int i = 0; i < nv; ++i) { */
/*     v = keys[i]; */
/*     k = kh_get(im, g->data, v); */
/*     tmp = kh_value(g->data, k); */
/*     kh_value(g->data, k) = ofx; */
/*     ofx += tmp; */
/*   } */
/*   /\* printf("Ofx: "); *\/ */
/*   /\* for (int i = 0; i < nv; ++i) { *\/ */
/*   /\*   v = keys[i]; *\/ */
/*   /\*   k = kh_get(im, g->data, v); *\/ */
/*   /\*   printf("%d ", kh_value(g->data, k)); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */

/*   /\* Computing ordering *\/ */
/*   for (int i = 0; i < pl; ++i) { */
/*     v = ph_get_v(paths, pi, i) >> 1; */
/*     k = kh_get(im, g->data, v); */
/*     p = kh_value(g->data, v); */
/*     kv_push(int, paths->orders, i); */
/*     ++kh_value(g->data, v); */
/*   } */
/*   /\* printf("Ordering: "); *\/ */
/*   /\* for (int i = 0; i < pl; ++i) { *\/ */
/*   /\*   printf("%d ", paths->orders[pi][i]); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */

/*   /\* printf("Ofx: "); *\/ */
/*   /\* for (int i = 0; i < nv; ++i) { *\/ */
/*   /\*   v = keys[i]; *\/ */
/*   /\*   k = kh_get(im, g->data, v); *\/ */
/*   /\*   printf("%d ", kh_value(g->data, k)); *\/ */
/*   /\* } *\/ */
/*   /\* printf("\n"); *\/ */

/*   /\* Reset offsets and store paths->offsets *\/ */
/*   v = keys[0]; */
/*   kv_push(int, paths->offsets, 0); */
/*   for (int i = 1; i < nv; ++i) { */
/*     k = kh_get(im, g->data, v); */
/*     v = keys[i]; */
/*     kv_push(int, paths->offsets, kh_value(g->data, k)); */
/*   } */

/*   /\* Reset offsets and store in keys array. Then delta compress into */
/*    * paths->offsets *\/ */
/*   /\* v = keys[0]; *\/ */
/*   /\* keys[0] = 0; *\/ */
/*   /\* for (int i = 1; i < nv; ++i) { *\/ */
/*   /\*   k = kh_get(im, g->data, v); *\/ */
/*   /\*   v = keys[i]; *\/ */
/*   /\*   keys[i] = kh_value(g->data, k); *\/ */
/*   /\* } *\/ */
/*   /\* store delta compressed offsets *\/ */
/*   /\* paths->offsets[pi] = malloc(p4nbound32(nv)); *\/ */
/*   /\* p4nd1enc32((uint32_t *)keys, nv, paths->offsets[pi]); *\/ */

/*   /\* printf("%ld %d\n", kv_size(path->vertices), reps_n); *\/ */
/* } */

int load_paths(graph_t *g) {
  kstring_t s = {0, 0, 0};
  int dret;
  gzFile fp = gzopen(g->fn, "r");
  if (fp == 0) {
    exit(1);
  }
  kstream_t *ks = ks_init(fp);

  /* g->data = kh_init(im); */

  while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
    if (s.s[0] == 'P' || s.s[0] == 'W') {
      /* current path: g->paths->np */
      /* if (g->paths->np == 1000) */
      /*   break; */
      /* if (s.s[0] == 'P') */
      /*   gfa_parse_P(s.s, g); */
      /* else */
      /*   gfa_parse_W(s.s, g); */
      /* compress_path(g, g->paths); */
      /* kh_clear(im, g->data); */
    }
  }
  /* ph_close(g->paths); */
  /* kh_destroy(im, g->data); */

  free(s.s);
  ks_destroy(ks);
  gzclose(fp);
  return 0;
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
            p_add_v(path, v, strand == '+');
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
            p_add_v(path, v, strand == '>');
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
