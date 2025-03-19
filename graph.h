#ifndef PS_GRAPH_H
#define PS_GRAPH_H

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "kseq.h"
#include "ksort.h"
#include "kvec.h"

#include "labels.h"
#include "misc.h"
#include "path.h"
#include "segment.h"

/* #include "ic.h" */

/* KHASH_MAP_INIT_INT(im, int) */

/* Some assumptions I made:
 * - vertex identifiers in GFA are integers
 * - graph has less than INT_MAX nodes
 */

typedef struct {
  char *fn;            /* gfa file name */
  labels_t *vertices;  /* labels in graph vertex order */
  khash_t(im) * v_map; /* mapping between gfa idx and internal idx (position on
                          vertices) */
  path_t *paths;

  int np;
  int cp;
} graph_t;

/* Initialize a graph */
graph_t *init_graph(char *fn);

/* Destroy the graph */
void destroy_graph(graph_t *g);

/* Load all vertices from gfa */
int load_vertices(graph_t *g, int wseq);

/* Load all paths from gfa */
int load_paths(graph_t *g);

/* Convert vertex from GFA space to graph space */
int get_iidx(graph_t *g, int v);

/* GFA reading utilities */
void gfa_parse_S(char *s, seg_t *ret, int wseq);
void gfa_parse_P(char *s, graph_t *g, path_t *p);
void gfa_parse_W(char *s, graph_t *g, path_t *p);

#endif
