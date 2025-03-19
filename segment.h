#ifndef PS_SEG_H
#define PS_SEG_H

#include <stdlib.h>

typedef struct {
  int idx;   /* identifier (as in gfa) - assuming integer */
  char *seq; /* sequence */
  int l;     /* length */
  int c;     /* capacity */
} seg_t;

/* Initialize a segment */
seg_t *init_seg();

/* Destroy the segment */
void destroy_seg(seg_t *seg);

#endif
