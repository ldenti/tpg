#include <getopt.h>
#include <stdlib.h>

#include "graph.hpp"
#include "misc.hpp"
#include "usage.hpp"

int main_build(int argc, char *argv[]) {
  int c;
  char *gbz_fn;
  bool force = false;

  static struct option long_options[] = {{"force", no_argument, NULL, 'f'},
                                         {"help", no_argument, NULL, 'h'},
                                         {NULL, 0, NULL, 0}};
  while ((c = getopt_long(argc, argv, "fh", long_options, NULL)) != -1) {
    switch (c) {
    case 'f':
      force = true;
      break;
    case 'h':
      fprintf(stderr, BUILD_USAGE_MESSAGE);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 1) {
    fprintf(stderr, BUILD_USAGE_MESSAGE);
    exit(EXIT_FAILURE);
  }
  gbz_fn = argv[optind++];

  double rt = realtime();
  Graph graph(gbz_fn, true);
  graph.load();
  fprintf(stderr, "[I::%s] loaded graph in in %.3f secs\n", __func__,
          realtime() - rt);
  if (!graph.has_horders || force) {
    rt = realtime();
    graph.build_horders();
    fprintf(stderr, "[I::%s] built horders in in %.3f secs\n", __func__,
            realtime() - rt);
    rt = realtime();
    graph.serialize_horders();
    fprintf(stderr, "[I::%s] built horders in in %.3f secs\n", __func__,
            realtime() - rt);
  } else {
    fprintf(stderr,
            "[I::%s] horders already on disk. Not recreating (use -f)\n",
            __func__);
  }

  return 0;
}
