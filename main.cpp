#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "misc.hpp"
#include "usage.hpp"

int main_build(int argc, char *argv[]);
int main_extract(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  double rt = realtime();

  if (argc < 2) {
    fprintf(stderr, "tpg [build|extract] -h\n");
    exit(EXIT_FAILURE);
  }

  if (strcmp(argv[1], "build") == 0)
    main_build(argc - 1, argv + 1);
  else if (strcmp(argv[1], "extract") == 0)
    main_extract(argc - 1, argv + 1);
  else {
    fprintf(stderr, "%s\n", MAIN_USAGE);
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "[I::%s] done in in %.3f secs\n", __func__, realtime() - rt);

  return 0;
}
