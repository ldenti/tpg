#include <filesystem>
#include <getopt.h>
#include <stdlib.h>
#include <string>

#include "gbwt/fast_locate.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

#include "misc.hpp"
#include "usage.hpp"

int main_build(int argc, char *argv[]) {
  double rt;
  int c;
  bool force = false;
  bool verbose = false;
  std::string gbz_fn;

  static struct option long_options[] = {{"force", no_argument, NULL, 'f'},
                                         {"verbose", no_argument, NULL, 'v'},
                                         {"help", no_argument, NULL, 'h'},
                                         {NULL, 0, NULL, 0}};
  while ((c = getopt_long(argc, argv, "fh", long_options, NULL)) != -1) {
    switch (c) {
    case 'f':
      force = true;
      break;
    case 'v':
      verbose = true;
      break;
    case 'h':
      fprintf(stderr, "%s", BUILD_USAGE_MESSAGE);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 1) {
    fprintf(stderr, "%s", BUILD_USAGE_MESSAGE);
    exit(EXIT_FAILURE);
  }
  gbz_fn = argv[optind++];

  // Build the R-index
  rt = realtime();
  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  fprintf(stderr,
          "[I::%s] restored graph (%lu paths, vertices in [%lld,%lld]) in %.3f "
          "secs\n",
          __func__, gbz.index.metadata.path_names.size(),
          gbz.graph.min_node_id(), gbz.graph.max_node_id(), realtime() - rt);
  if (verbose)
    gbwt::printStatistics(gbz.index, gbz_fn, std::cerr);

  std::string ri_fn = gbz_fn + ".ri";
  if (!std::filesystem::exists(ri_fn) || force) {
    rt = realtime();
    gbwt::FastLocate fl;
    fl = gbwt::FastLocate(gbz.index);
    fprintf(stderr, "[I::%s] built R-index in %.3f sec\n", __func__,
            realtime() - rt);

    rt = realtime();
    std::ofstream out;
    out.open(ri_fn, std::ofstream::out);
    fl.serialize(out);
    out.close();
    fprintf(stderr, "[M::%s] stored R-index in %.3f sec\n", __func__,
            realtime() - rt);
  }

  return 0;
}
