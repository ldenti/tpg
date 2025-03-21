#include "segments.h"

/* === Single segment ==============================================*/

seg_t *init_seg() {
  seg_t *seg = malloc(1 * sizeof(seg_t));
  seg->l = 0;
  // seg->idx = (char *)malloc(1024);
  seg->seq = NULL; /* malloc(4096 * sizeof(char)); */
  seg->c = 0;

  return seg;
}

void destroy_seg(seg_t *seg) {
  if (seg->seq != NULL)
    free(seg->seq);
  free(seg);
}

/* ================================================================ */

/* === Multiple segments ========================================== */

segments_t *sgms_init() {
  segments_t *sgms = malloc(sizeof(segments_t));
  sgms->text_c = 4096;
  sgms->text = malloc(sgms->text_c * sizeof(char));
  sgms->text_n = 0;

  sgms->ofx_c = 256;
  sgms->n = 0;
  sgms->ofx = malloc(sgms->ofx_c * sizeof(uint64_t));
  sgms->vnames = malloc(sgms->ofx_c * sizeof(int));

  return sgms;
}

void sgms_destroy(segments_t *sgms) {
  free(sgms->text);
  free(sgms->ofx);
  free(sgms->vnames);
  free(sgms);
}

int sgms_add(segments_t *sgms, int idx, char *label, int l) {
  if (sgms->text_n + l >= sgms->text_c) {
    sgms->text =
        realloc(sgms->text, (uint64_t)(sgms->text_n + l) * 2 * sizeof(char));
    sgms->text_c = (uint64_t)(sgms->text_n + l) * 2;
  }
  strncpy(sgms->text + sgms->text_n, label, l);

  if (sgms->n == sgms->ofx_c) {
    sgms->vnames = realloc(sgms->vnames, 2 * sgms->ofx_c * sizeof(int));
    sgms->ofx = realloc(sgms->ofx, 2 * sgms->ofx_c * sizeof(uint64_t));
    sgms->ofx_c *= 2;
  }

  sgms->vnames[sgms->n] = idx;
  sgms->ofx[sgms->n] = sgms->text_n;
  ++sgms->n;
  sgms->text_n += l;

  return 0;
}

int sgms_get_l(segments_t *sgms, int i) {
  return sgms->ofx[i + 1] - sgms->ofx[i];
}

void sgms_get(segments_t *sgms, int i, char **label, int *c) {
  int l = sgms_get_l(sgms, i);
  if (l + 1 >= *c) {
    *label = realloc(*label, 2 * l);
    *c = 2 * l;
  }
  strncpy(*label, sgms->text + sgms->ofx[i], l);
  *(*label + l) = '\0';
}

/* ================================================================ */
