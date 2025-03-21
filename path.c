#include "path.h"

path_t *ph_init() {
  path_t *path = malloc(1 * sizeof(path_t));
  path->idx = malloc(128 * sizeof(char));

  kv_init(path->vertices);
  path->data = kh_init(im);

  return path;
}

void ph_clear(path_t *path) {
  path->vertices.n = 0;
  kh_clear(im, path->data);
}

void ph_destroy(path_t *path) {
  free(path->idx);
  kv_destroy(path->vertices);
  kh_destroy(im, path->data);
  free(path);
}

void ph_addv(path_t *path, int v, int strand) {
  /* Update hash table */
  khiter_t k;
  int hret;
  k = kh_put(im, path->data, (v << 4), &hret);
  assert(hret < 2);
  if (hret == 1) {
    kh_value(path->data, k) = 0;
  }
  ++kh_value(path->data, k);
  int ev = p_encode(v, kh_value(path->data, k));

  k = kh_put(im, path->data, ev, &hret);
  assert(hret == 1);
  kh_value(path->data, k) = kv_size(path->vertices);

  /* Add vertex to array */
  kv_push(uint, path->vertices, (v << 1) | strand);
}

int ph_getv(path_t *path, int p) { return kv_A(path->vertices, p) >> 1; }
