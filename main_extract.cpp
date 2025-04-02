#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <getopt.h>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "graph.hpp"
#include "misc.hpp"
#include "usage.hpp"

typedef struct {
  std::string chrom;
  uint s;
  uint e;
} region_t;

int main_extract(int argc, char *argv[]) {
  double rt;
  int c;
  std::string fn;
  std::string bed_fn;
  bool is_gbz = true;
  bool verbose = false;
  std::string ref_name = "CHM13";

  static struct option long_options[] = {
      {"gfa", no_argument, NULL, 'g'},
      {"refpath", required_argument, NULL, 'r'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      {NULL, 0, NULL, 0}};
  while ((c = getopt_long(argc, argv, "gr:vh", long_options, NULL)) != -1) {
    switch (c) {
    case 'g':
      is_gbz = false;
      break;
    case 'r':
      ref_name = optarg;
      break;
    case 'v':
      verbose = true;
      break;
    case 'h':
      fprintf(stderr, "%s", EXTRACT_USAGE_MESSAGE);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 2) {
    fprintf(stderr, "%s", EXTRACT_USAGE_MESSAGE);
    exit(EXIT_FAILURE);
  }
  fn = argv[optind++];
  bed_fn = argv[optind++];

  // Load graph
  rt = realtime();
  Graph graph(fn, is_gbz);
  graph.load();
  if (!graph.has_horders) {
    fprintf(stderr, "[E::%s] graph has no horders. Halting\n", __func__);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr,
          "[I::%s] loaded graph (%lu paths, vertices in [%lld,%lld]) in %.3f "
          "secs\n",
          __func__, graph.gbz.index.metadata.path_names.size(),
          graph.gbz.graph.min_node_id(), graph.gbz.graph.max_node_id(),
          realtime() - rt);

  // Parse BED file
  rt = realtime();
  char line[256];
  int i;
  int regions_n = 0;
  char *p, *q;

  std::vector<region_t> regions;
  std::map<std::string, int> paths; // Prefix sum of paths

  FILE *file = fopen(bed_fn.c_str(), "r");
  if (file == NULL) {
    fprintf(stderr, "[E::%s] Unable to open file\n", __func__);
    exit(EXIT_FAILURE);
  }
  while (fgets(line, sizeof(line), file)) {
    for (i = 0, p = q = line;; ++p) {
      if (*p == '\t' || *p == 0) {
        *p = 0;
        if (i == 0) {
          regions.push_back({std::string(q, p - q), 0, 0});
          if (paths.find(regions[regions_n].chrom) == paths.end()) {
            paths[regions[regions_n].chrom] = regions_n;
          }
        } else if (i == 1) {
          regions[regions_n].s = atoi(q);
        } else if (i == 2) {
          regions[regions_n].e = atoi(q);
          break;
        } else {
          std::cerr << "Error while reading BED file" << std::endl;
          exit(EXIT_FAILURE);
        }
        ++i;
        q = p + 1;
      }
    }
    ++regions_n;
  }
  fclose(file);
  fprintf(stderr, "[I::%s] parsed BED file in %.3f secs\n", __func__,
          realtime() - rt);

  // Exract source/sink for each region
  rt = realtime();
  std::vector<region_t> ends; // Sink/source pairs for "each" region

  for (size_t i = 0; i < graph.gbz.index.metadata.path_names.size(); ++i) {
    std::string sample_name = graph.gbz.index.metadata.fullPath(i).sample_name;
    std::string contig_name = graph.gbz.index.metadata.fullPath(i).contig_name;
    // size_t offset = graph.gbz.index.metadata.fullPath(i).offset;
    // XXX: assuming offset == 0 since we are on a reference path

    if (sample_name.compare(ref_name) != 0)
      // we are not on the reference path
      continue;
    if (paths.find(contig_name) == paths.end()) {
      // we have no regions on this chromosome
      continue;
    }

    gbwt::size_type path_idx = gbwt::Path::encode(i, 0);
    // XXX: assuming reference path to be on + strand. Is this enough?
    gbwt::vector_type path = graph.gbz.index.extract(path_idx);

    int p = paths[contig_name]; // index of first region on this contig, using
                                // prefix sum

    if (verbose)
      fprintf(stderr, "[I::%s] traversing %s:%s\n", __func__,
              sample_name.c_str(), contig_name.c_str());

    int sp, ep;
    while (p < regions_n && contig_name.compare(regions[p].chrom) == 0) {
      sp = regions[p].s;
      ep = regions[p].e;
      int offset = 0;
      uint source = -1, sink = -1;

      for (const unsigned int &v : path) {
        // v is the encoded vertex (with strand bit)
        // vv is the gbwtgraph identifier
        int vv = gbwt::Node::id(v);
        assert(graph.gbz.graph.has_node(vv));
        gbwtgraph::handle_t vh = graph.gbz.graph.node_to_handle(v);

        offset += graph.gbz.graph.get_length(vh);
        if (source == (uint)-1 && offset >= sp) {
          source = vv;
        }
        if (source != (uint)-1 && sink == (uint)-1 && offset >= ep) {
          sink = vv;
          break;
        }
      }

      if (source != (uint)-1 && sink != (uint)-1) {
        ends.push_back(
            {contig_name + ":" + std::to_string(sp) + "-" + std::to_string(ep),
             source, sink});
        if (verbose)
          fprintf(stderr, "[I::%s] hit: %s %s (%d) %s (%d) --- %s\n", __func__,
                  regions[p].chrom.c_str(), graph.get_gfa_idx(source).c_str(),
                  sp, graph.get_gfa_idx(sink).c_str(), ep,
                  ends.back().chrom.c_str());
      }
      ++p;
    }
  }
  fprintf(stderr, "[I::%s] found %lu/%lu regions in %.3f secs\n", __func__,
          ends.size(), regions.size(), realtime() - rt);

  // Reiterate over paths to get subhaplotypes
  rt = realtime();

  std::set<std::pair<std::string, std::pair<uint32_t, uint32_t>>>
      VV; // we cannot encode node identifiers or handles since vertices longer
          // than 1024 are broken down. Only way I found to store each vertex
          // only once (one element per "real" vertex) is by using this pair
  // all other identifiers in edges and subhaps are encoded (has strand bit
  // since we need it in GFA)
  std::set<uint64_t> EE;
  std::vector<std::vector<uint32_t>> subhaps;
  std::vector<std::string> subhaps_names;

  for (const region_t &ss : ends) {
    uint source = ss.s;
    uint sink = ss.e;

    if (verbose)
      fprintf(stderr, "=== %s>%s ===\n", graph.get_gfa_idx(source).c_str(),
              graph.get_gfa_idx(sink).c_str());
    // sources
    std::map<gbwt::size_type, std::vector<gbwt::size_type>> source_hits =
        graph.locate(source);
    if (source_hits.size() == 0)
      continue;
    std::vector<gbwt::size_type> source_hits_paths;
    for (const auto &h : source_hits) {
      source_hits_paths.push_back(h.first);
    }
    // sinks
    std::map<gbwt::size_type, std::vector<gbwt::size_type>> sink_hits =
        graph.locate(sink);
    if (sink_hits.size() == 0)
      continue;
    std::vector<gbwt::size_type> sink_hits_paths;
    for (const auto &h : sink_hits) {
      sink_hits_paths.push_back(h.first);
    }
    // intersection
    std::vector<gbwt::size_type> intersection;
    std::sort(source_hits_paths.begin(), source_hits_paths.end());
    std::sort(sink_hits_paths.begin(), sink_hits_paths.end());
    std::set_intersection(source_hits_paths.begin(), source_hits_paths.end(),
                          sink_hits_paths.begin(), sink_hits_paths.end(),
                          std::back_inserter(intersection));
    // we iterate over paths passing through both source and sink

    for (const gbwt::size_type &p : intersection) {
      // path_idx is encoding, with strand. we need the real id (pp) to query
      // the gbwtgraph
      int pp = gbwt::Path::id(p);
      std::string sample_name =
          graph.gbz.index.metadata.fullPath(pp).sample_name;
      std::string contig_name =
          graph.gbz.index.metadata.fullPath(pp).contig_name;
      size_t haplotype = graph.gbz.index.metadata.fullPath(pp).haplotype;

      // XXX: I think checking one strand is enough, so I take +
      gbwt::size_type path_idx =
          gbwt::Path::encode(pp, 0); // this could be p (if p is not on -)
      gbwt::vector_type path = graph.gbz.index.extract(path_idx);
      assert(!path.empty());
      if (path.empty()) {
        fprintf(stderr, "!!! %s %s\n", sample_name.c_str(),
                contig_name.c_str());
        continue;
      }

      int p1, p2;
      p1 = *min_element(source_hits[p].begin(), source_hits[p].end());
      p2 = *max_element(sink_hits[p].begin(), sink_hits[p].end());

      if (p1 > p2) {
        // Some paths may be inverted. We just invert the two positions
        std::swap(p1, p2);
        // if (sample_name.compare("HG03579") == 0) {
        //   for (const auto &x : path)
        //     std::cerr << graph.get_gfa_idx(gbwt::Node::id(x)) << " ";
        //   std::cerr << std::endl;
        //   for (const auto &x : source_hits[p])
        //     std::cerr << x << " ";
        //   std::cerr << std::endl;
        //   for (const auto &x : sink_hits[p])
        //     std::cerr << x << " ";
        //   std::cerr << std::endl;
        // }
      }
      if (verbose)
        fprintf(stderr, "[W::%s] extracting from %s#%s\n", __func__,
                sample_name.c_str(), contig_name.c_str());

      gbwtgraph::subpath_type subpath =
          gbwtgraph::get_subpath(path, p1, p2 + 1);
      VV.insert(graph.gbz.graph.get_segment(
          graph.gbz.graph.node_to_handle(subpath.first[0])));
      subhaps_names.push_back(sample_name + "#" + std::to_string(haplotype) +
                              "#" + ss.chrom);
      subhaps.push_back({});
      subhaps.back().push_back(subpath.first[0]);
      for (uint j = 1; j < subpath.second; ++j) {
        VV.insert(graph.gbz.graph.get_segment(
            graph.gbz.graph.node_to_handle(subpath.first[j])));
        EE.insert(((uint64_t)subpath.first[j - 1] << 32) |
                  ((uint64_t)subpath.first[j]));
        subhaps.back().push_back(subpath.first[j]);
      }
    }
  }
  fprintf(stderr,
          "[I::%s] extracted %lu subhaplotypes (%lu vertices) in %.3f secs\n",
          __func__, subhaps.size(), VV.size(), realtime() - rt);

  // Print GFA
  rt = realtime();
  for (const std::pair<std::string, std::pair<uint32_t, uint32_t>> &v : VV) {
    printf("S\t%s\t", v.first.c_str());
    for (size_t i = v.second.first; i < v.second.second; ++i) {
      printf(
          "%s",
          graph.gbz.graph.get_sequence(graph.gbz.graph.get_handle(i)).c_str());
    }
    printf("\n");
  }

  for (const auto &ee : EE) {
    uint32_t v1 = ee >> 32;
    uint32_t v2 = (uint32_t)ee;

    printf("L\t%s\t%c\t%s\t%c\t0M\n",
           graph.get_gfa_idx(gbwt::Node::id(v1)).c_str(),
           gbwt::Node::is_reverse(v1) ? '-' : '+',
           graph.get_gfa_idx(gbwt::Node::id(v2)).c_str(),
           gbwt::Node::is_reverse(v2) ? '-' : '+');
  }

  for (size_t i = 0; i < subhaps.size(); ++i) {
    printf("P\t%s\t%s%c", subhaps_names[i].c_str(),
           graph.get_gfa_idx(gbwt::Node::id(subhaps[i][0])).c_str(),
           gbwt::Node::is_reverse(subhaps[i][0]) ? '-' : '+');
    // this works since we are sure to have at least one vertex
    for (uint j = 1; j < subhaps[i].size(); ++j)
      printf(",%s%c", graph.get_gfa_idx(gbwt::Node::id(subhaps[i][j])).c_str(),
             gbwt::Node::is_reverse(subhaps[i][j]) ? '-' : '+');
    printf("\t*\n");
  }

  return 0;
}
