#ifndef PS_LAB_H
#define PS_LAB_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *text;      /* concatenation of labels */
  uint64_t text_c; /* text capacity */
  uint64_t text_n; /* text actual length */
  uint64_t *ofx;   /* offsets along text */
  int *vnames;     /* vertex identifiers (as in GFA) */
  int ofx_c;       /* offsets capacity */
  int n;           /* number of vertices, aka number of offsets */
} labels_t;

labels_t *l2b_init();
void l2b_destroy(labels_t *l2b);
int l2b_add(labels_t *l2b, int idx, char *label, int l);
int l2b_get_l(labels_t *l2b, int i);
void l2b_get(labels_t *l2b, int i, char *label, int *c);

#endif
