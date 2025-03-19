#include "labels.h"

labels_t *l2b_init() {
  labels_t *l2b = malloc(sizeof(labels_t));
  l2b->text_c = 4096;
  l2b->text = malloc(l2b->text_c * sizeof(char));
  l2b->text_n = 0;

  l2b->ofx_c = 256;
  l2b->n = 0;
  l2b->ofx = malloc(l2b->ofx_c * sizeof(uint64_t));
  l2b->vnames = malloc(l2b->ofx_c * sizeof(int));

  return l2b;
}

void l2b_destroy(labels_t *l2b) {
  free(l2b->text);
  free(l2b->ofx);
  free(l2b->vnames);
  free(l2b);
}

int l2b_add(labels_t *l2b, int idx, char *label, int l) {
  if (l2b->text_n + l >= l2b->text_c) {
    l2b->text =
        realloc(l2b->text, (uint64_t)(l2b->text_n + l) * 2 * sizeof(char));
    l2b->text_c = (uint64_t)(l2b->text_n + l) * 2;
  }
  strncpy(l2b->text + l2b->text_n, label, l);

  if (l2b->n == l2b->ofx_c) {
    l2b->vnames = realloc(l2b->vnames, 2 * l2b->ofx_c * sizeof(int));
    l2b->ofx = realloc(l2b->ofx, 2 * l2b->ofx_c * sizeof(uint64_t));
    l2b->ofx_c *= 2;
  }

  l2b->vnames[l2b->n] = idx;
  l2b->ofx[l2b->n] = l2b->text_n;
  ++l2b->n;
  l2b->text_n += l;

  return 0;
}

int l2b_get_l(labels_t *l2b, int i) { return l2b->ofx[i + 1] - l2b->ofx[i]; }

void l2b_get(labels_t *l2b, int i, char *label, int *c) {
  int l = l2b_get_l(l2b, i);
  if (l + 1 >= *c) {
    label = realloc(label, 2 * l);
    *c = 2 * l;
  }
  strncpy(label, l2b->text + l2b->ofx[i], l + 1);
  label[l] = '\0';
}
