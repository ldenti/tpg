#ifndef PS_PATH_H
#define PS_PATH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "khash.h"
#include "ksort.h"
#include "kvec.h"
#include "rope.h"

KHASH_MAP_INIT_INT(im, int)

static inline int p_encode(int v, int n) { return (v << 4) | (n & 0xF); }

typedef struct {
  char *idx;             /* path identifier (as in gfa) */
  kvec_t(uint) vertices; /* identifiers of the vertices in path order:  31
                          * bits for id in graph space, 1 bit for strand */

  /* just in construction */
  khash_t(im) * data; /* vertex membership and number of occurrences */
  /* --- */

  rope_t *rope;        /* Membership run-length encoded bit vector */
  kvec_t(int) orders;  /* Path ordering. Size: V */
  kvec_t(int) offsets; /* Offsets along orders array. Size: V */
} path_t;

/* Initialize a path for a graph with nv vertices */
path_t *init_path();

void clear_path(path_t *path);

int p_build_ds(path_t *path);

/* Destroy a path */
void destroy_path(path_t *path);

/* Add vertex v to path, reallocating it if needed */
void p_add_v(path_t *path, int v, int strand);

/* Get p-th vertex along the path */
int p_get_v(path_t *path, int p);

#endif
