#ifndef PS_GRAPH_HPP
#define PS_GRAPH_HPP

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "khash.h"
#include "kseq.h"

#include "misc.hpp"
#include "path.hpp"
#include "segments.hpp"

#include "gbwt/gbwt.h"

/* Some assumptions I made:
 * - vertex identifiers in GFA are integers
 * - graph has less than INT_MAX nodes
 * - we do not need to store all paths here
 */

typedef struct {
  char *fn;             /* gfa file name */
  segments_t *vertices; /* labels in graph vertex order */
  khash_t(im) * v_map;  /* mapping between gfa idx and internal idx (position on
                           vertices) */
  /* We do not need to load all paths */
  /* path_t *paths; */
} graph_t;

/* Initialize a graph */
graph_t *init_graph(char *fn);

/* Destroy the graph */
void destroy_graph(graph_t *g);

/* Load all vertices from gfa */
int load_vertices(graph_t *g, int wseq);

/* Load all paths from gfa */
/* int load_paths(graph_t *g); */

/* Convert vertex from GFA space to graph space */
int get_iidx(graph_t *g, int v);

/* GFA reading utilities */
void gfa_parse_S(char *s, seg_t *ret, int wseq);
void gfa_parse_P(char *s, graph_t *g, path_t *p);
void gfa_parse_W(char *s, graph_t *g, path_t *p);

#endif
