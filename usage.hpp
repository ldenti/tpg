#ifndef TPG_USAGE_HPP
#define TPG_USAGE_HPP

static const char *const VERSION = "tpg, v0.0.1";

static const char *const MAIN_USAGE = "Usage: tpg [build|extract] -h\n";

static const char *const BUILD_USAGE_MESSAGE =
    "Usage: tpg build [options] <graph.gbz>\n"
    "Options:\n"
    "        -f,--force      force construction\n"
    "        -v,--verbose    print more messages to stderr\n"
    "        -h,--help       display this help and exit\n"
    "\n";

static const char *const EXTRACT_USAGE_MESSAGE =
    "Usage: tpg extract [options] <graph.gbz> <regions.bed>\n"
    "Options:\n"
    "        -r,--rp         use this as reference path (default: CHM13)\n"
    "        -f,--fasta      output FASTA (default: GFA)\n"
    "        -v,--verbose    print more messages to stderr\n"
    "        -h,--help       display this help and exit\n"
    "\n";

#endif
