#include "path.h"
#include "khash.h"

KSORT_INIT_GENERIC(int)

path_t *init_path() {
  path_t *path = malloc(1 * sizeof(path_t));
  path->idx = malloc(128 * sizeof(char));

  kv_init(path->vertices);
  path->data = kh_init(im);

  return path;
}

void clear_path(path_t *path) {
  path->vertices.n = 0;
  kh_clear(im, path->data);
}

void destroy_path(path_t *path) {
  free(path->idx);
  kv_destroy(path->vertices);
  kh_destroy(im, path->data);
  free(path);
}

void p_add_v(path_t *path, int v, int strand) {
  kv_push(uint, path->vertices, (v << 1) | strand);

  khiter_t k;
  int hret;

  /* This if we want "custom" data structure */
  /* k = kh_put(im, path->data, v, &hret); */
  /* if (ret != 0) */
  /*   /\* new element *\/ */
  /*   kh_value(path->data, k) = 0; */
  /* ++kh_value(path->data, k); */

  /* This if we want just an hash table */
  k = kh_put(im, path->data, (v << 4), &hret);
  assert(hret < 2);
  if (hret == 1) {
    kh_value(path->data, k) = 0;
  }
  ++kh_value(path->data, k);
  int ev = p_encode(v, kh_value(path->data, k));

  k = kh_put(im, path->data, ev, &hret);
  assert(hret == 1);
  kh_value(path->data, k) =
      kv_size(path->vertices) -
      1; /* -1 since we already added the current vertex */
}

int p_get_v(path_t *path, int p) { return kv_A(path->vertices, p) >> 1; }

/* TODO if we want custom data structure */
/* int p_build_ds(path_t *path) { */
/*   int plen = kv_size(path->vertices); */

/*   /\* Get sorted array of "used" vertices *\/ */
/*   int *keys = malloc(plen * sizeof(int)); */
/*   int nv = 0; */
/*   khiter_t k; */
/*   for (k = kh_begin(path->data); k != kh_end(path->data); ++k) */
/*     if (kh_exist(path->data, k)) */
/*       keys[nv++] = kh_key(path->data, k); */
/*   ks_introsort(int, nv, keys); */
/* /\* Build run-length encoded bit vector for path *\/ */
/*   paths->ropes[pi] = rope_init(ROPE_DEF_MAX_NODES, ROPE_DEF_BLOCK_LEN);
 */
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

/*   /\* Get offsets for each vertex using its repetitions. Replace
values in
 */
/* data */
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

/*   /\* Reset offsets and store in keys array. Then delta compress into
 */
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

/*   return 0; */
/* } */
