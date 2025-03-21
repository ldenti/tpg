#ifndef PS_PATH_H
#define PS_PATH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khash.h"
#include "kvec.h"

KHASH_MAP_INIT_INT(im, int)

static inline int p_encode(int v, int n) { return (v << 4) | (n & 0xF); }

typedef struct {
  char *idx;             /* path identifier (as in gfa) */
  kvec_t(uint) vertices; /* identifiers of the vertices in path order:  31
                          * bits for id in graph space, 1 bit for strand */
  khash_t(im) * data;    /* vertex membership and number of occurrences */
} path_t;

/* Initialize a path for a graph with nv vertices */
path_t *ph_init();

/* Clear a path without deallocating almost anything */
void ph_clear(path_t *path);

/* Destroy a path */
void ph_destroy(path_t *path);

/* Add vertex v to path, reallocating it if needed */
void ph_addv(path_t *path, int v, int strand);

/* Get p-th vertex along the path */
int ph_getv(path_t *path, int p);

#endif
