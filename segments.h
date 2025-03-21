#ifndef PS_SEGS_H
#define PS_SEGS_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Single segment (S-line) */
typedef struct {
  int idx;   /* identifier (as in gfa) - assuming integer */
  char *seq; /* sequence */
  int l;     /* length */
  int c;     /* capacity */
} seg_t;

/* Concatenation of multiple segments */
typedef struct {
  char *text;      /* concatenation of labels */
  uint64_t text_c; /* text capacity */
  uint64_t text_n; /* text actual length */
  uint64_t *ofx;   /* offsets along text */
  int *vnames;     /* vertex identifiers (as in GFA) */
  int ofx_c;       /* offsets capacity */
  int n;           /* number of vertices, aka number of offsets */
} segments_t;

/* Initialize a segment */
seg_t *init_seg();

/* Destroy the segment */
void destroy_seg(seg_t *seg);

/* Initialize segments */
segments_t *sgms_init();

/* Destroy segments */
void sgms_destroy(segments_t *sgms);

/* Push a new segment */
int sgms_add(segments_t *sgms, int idx, char *label, int l);

/* Get length of i-th segment */
int sgms_get_l(segments_t *sgms, int i);

/* Get label of i-th segment (+ reallocate it if needed) */
void sgms_get(segments_t *sgms, int i, char **label, int *c);

#endif
