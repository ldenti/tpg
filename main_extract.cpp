#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "gbwt/fast_locate.h"
#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

#include "misc.hpp"
#include "usage.hpp"

typedef struct {
  std::string chrom;
  uint s;
  uint e;
  std::string name;
} region_t; // BED line

typedef struct {
  std::string chrom;
  uint v1;
  uint offset1;
  uint v2;
  uint offset2;
  std::string name;
} hit_t; // BED line

typedef struct {
  std::string chrom;      // chromosome
  std::string name;       // name, i.e., same as BED region
  gbwt::vector_type path; // path
  uint offset1, offset2;  // offsets on first and last vertex
} subhap_t;

std::string get_gfa_idx(const gbwtgraph::GBZ gbz, gbwt::node_type v) {
  return gbz.graph.get_segment_name(gbz.graph.get_handle(v));
}

std::string get_gfa_seq(const gbwtgraph::GBZ gbz, gbwt::node_type v) {
  return gbz.graph.get_sequence(gbz.graph.get_handle(v));
}

int main_extract(int argc, char *argv[]) {
  double rt;
  int c;
  std::string gbz_fn;
  std::string bed_fn;
  std::string ref_name = "CHM13";
  bool fasta = false;
  bool verbose = false;

  static struct option long_options[] = {
      {"refpath", required_argument, NULL, 'r'},
      {"fasta", required_argument, NULL, 'f'},
      {"verbose", no_argument, NULL, 'v'},
      {"help", no_argument, NULL, 'h'},
      {NULL, 0, NULL, 0}};
  while ((c = getopt_long(argc, argv, "r:fvh", long_options, NULL)) != -1) {
    switch (c) {
    case 'r':
      ref_name = optarg;
      break;
    case 'f':
      fasta = true;
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
  gbz_fn = argv[optind++];
  bed_fn = argv[optind++];

  // === Load graph ===
  rt = realtime();
  gbwtgraph::GBZ gbz;
  sdsl::simple_sds::load_from(gbz, gbz_fn);
  fprintf(stderr,
          "[I::%s] restored graph (%lu paths, vertices in [%lld,%lld]) in %.3f "
          "secs\n",
          __func__, gbz.index.metadata.path_names.size(),
          gbz.graph.min_node_id(), gbz.graph.max_node_id(), realtime() - rt);

  std::string ri_fn = gbz_fn + ".ri";
  if (!std::filesystem::exists(ri_fn)) {
    fprintf(stderr, "[E::%s] cannot find %s\n", __func__, ri_fn.c_str());
    exit(EXIT_FAILURE);
  }
  gbwt::FastLocate fl;
  std::ifstream in;
  in.open(ri_fn, std::ifstream::in);
  fl.load(in);
  fl.setGBWT(gbz.index);
  in.close();
  fprintf(stderr, "[M::%s] restored R-index in %.3f sec\n", __func__,
          realtime() - rt);

  // === Parse BED file ===
  rt = realtime();
  char line[256];
  int i;
  int regions_n = 0;
  char *p, *q;

  std::vector<region_t> regions;
  std::map<std::string, int> paths; // Prefix sum of paths

  FILE *file = fopen(bed_fn.c_str(), "r");
  if (file == NULL) {
    fprintf(stderr, "[E::%s] Unable to open BED file %s\n", __func__,
            bed_fn.c_str());
    exit(EXIT_FAILURE);
  }
  while (fgets(line, sizeof(line), file)) {
    for (i = 0, p = q = line;; ++p) {
      if (*p == '\n' || *p == '\t' || *p == 0) {
        *p = 0;
        if (i == 0) {
          regions.push_back({std::string(q, p - q), 0, 0, ""});
          if (paths.find(regions[regions_n].chrom) == paths.end()) {
            paths[regions[regions_n].chrom] = regions_n;
          }
        } else if (i == 1) {
          regions[regions_n].s = atoi(q);
        } else if (i == 2) {
          regions[regions_n].e = atoi(q);
        } else if (i == 3) {
          regions[regions_n].name = q;
          break;
        } else {
          fprintf(stderr, "[E::%s] Error while reading BED file\n", __func__);
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

  // === Exract source/sink for each region ===
  rt = realtime();
  std::vector<hit_t> ends; // Sink/source pairs for "each" region
  for (size_t i = 0; i < gbz.index.metadata.path_names.size(); ++i) {
    std::string sample_name = gbz.index.metadata.fullPath(i).sample_name;
    std::string contig_name = gbz.index.metadata.fullPath(i).contig_name;
    size_t offset = gbz.index.metadata.fullPath(i).offset;

    if (sample_name.compare(ref_name) != 0)
      // we are not on the reference path
      continue;
    if (paths.find(contig_name) == paths.end()) {
      // we have no regions on this chromosome
      continue;
    }

    // assuming offset == 0 since we are on a reference path
    if (offset > 0) {
      fprintf(stderr,
              "[E::%s] path %s is not end-to-end (offset: %ld). Please, select "
              "correct reference path using `-r` option\n",
              __func__, sample_name.c_str(), offset);
      exit(EXIT_FAILURE);
    }

    gbwt::size_type path_idx = gbwt::Path::encode(i, 0);
    // XXX: assuming reference path to be on + strand. Is this enough?
    gbwt::vector_type path = gbz.index.extract(path_idx);

    int p = paths[contig_name]; // index of first region on this contig, using
                                // prefix sum

    if (verbose)
      fprintf(stderr, "[I::%s] traversing %s:%s\n", __func__,
              sample_name.c_str(), contig_name.c_str());

    int sp, ep;
    std::string name;
    while (p < regions_n && contig_name.compare(regions[p].chrom) == 0) {
      sp = regions[p].s;
      ep = regions[p].e;
      name = regions[p].name;
      int offset = 0;
      uint source = -1, sink = -1;
      uint source_offset = -1, sink_offset = -1;

      for (const unsigned int &v : path) {
        // v is the encoded vertex (with strand bit)
        // vv is the gbwtgraph identifier
        int vv = gbwt::Node::id(v);
        assert(gbz.graph.has_node(vv));
        gbwtgraph::handle_t vh = gbz.graph.node_to_handle(v);

        offset += gbz.graph.get_length(vh);
        if (source == (uint)-1 && offset >= sp) {
          source = vv;
          source_offset = offset - sp;
        }
        if (source != (uint)-1 && sink == (uint)-1 && offset >= ep) {
          sink = vv;
          sink_offset = offset - ep;
          break;
        }
      }

      if (source != (uint)-1 && sink != (uint)-1) {
        ends.push_back(
            {contig_name + ":" + std::to_string(sp) + "-" + std::to_string(ep),
             source, source_offset, sink, sink_offset, name});
        if (verbose)
          fprintf(stderr, "[I::%s] hit: %s %s (%d) %s (%d) --- %s\n ", __func__,
                  regions[p].chrom.c_str(), get_gfa_idx(gbz, source).c_str(),
                  sp, get_gfa_idx(gbz, sink).c_str(), ep,
                  ends.back().chrom.c_str());
      }
      ++p;
    }
  }
  fprintf(stderr, "[I::%s] found %lu/%lu regions in %.3f secs\n", __func__,
          ends.size(), regions.size(), realtime() - rt);

  // Reiterate over paths to get subhaplotypes
  rt = realtime();

  std::vector<subhap_t> subhaps;
  for (const hit_t &ss : ends) {
    std::string chrom = ss.chrom;
    uint source = ss.v1;
    uint source_offset = ss.offset1;
    uint sink = ss.v2;
    uint sink_offset = ss.offset2;
    std::string name = ss.name;

    if (verbose)
      fprintf(stderr, "=== %s>%s ===\n", get_gfa_idx(gbz, source).c_str(),
              get_gfa_idx(gbz, sink).c_str());

    for (int strand1 = 0; strand1 < 2; ++strand1) {
      for (int strand2 = 0; strand2 < 2; ++strand2) {
        uint e_source = gbwt::Node::encode(source, strand1);
        uint e_sink = gbwt::Node::encode(sink, strand2);
        std::vector<gbwt::size_type> pvisits1 = fl.decompressSA(e_source);
        std::vector<gbwt::size_type> pvisits2 = fl.decompressSA(e_sink);

        for (uint i = 0; i < pvisits1.size(); ++i) {
          for (uint j = 0; j < pvisits2.size(); ++j) {
            gbwt::size_type pidx1 = fl.seqId(pvisits1[i]);
            gbwt::size_type pidx2 = fl.seqId(pvisits2[j]);
            if (pidx1 != pidx2)
              continue;
            gbwt::size_type pp = gbwt::Path::id(pidx1);
            std::string sample_name =
                gbz.index.metadata.fullPath(pp).sample_name;
            std::string contig_name =
                gbz.index.metadata.fullPath(pp).contig_name;
            int haplotype = gbz.index.metadata.fullPath(pp).haplotype;
            // if (sample_name.compare("_gbwt_ref") != 0)
            //   continue;
            if (fl.seqOffset(pvisits1[i]) < fl.seqOffset(pvisits2[j]))
              continue;

            subhap_t subhap{chrom,
                            name + "#" + sample_name + "#" +
                                std::to_string(haplotype),
                            {},
                            source_offset,
                            sink_offset};
            gbwt::edge_type position = std::make_pair(e_source, i);
            //   if (!gbz.index.contains(position)) {
            //     vv1 = gbwt::Node::encode(v1, 1);
            //     position = std::make_pair(vv1, k);
            //     assert(gbz.index.contains(position));
            //   }

            while (position.first != e_sink) {
              // position.first is encoded with strand
              subhap.path.push_back(position.first);
              position = gbz.index.LF(position);
            }
            subhap.path.push_back(position.first);

            subhaps.push_back(subhap);
          }
        }
      }
    }
  }
  fprintf(stderr, "[I::%s] extracted %lu subhaplotypes in %.3f secs\n",
          __func__, subhaps.size(), realtime() - rt);

  rt = realtime();
  if (fasta) {
    // TODO: cut prefix and suffix
    for (const subhap_t &subhap : subhaps) {
      std::cout << ">" << subhap.name << " " << subhap.chrom << std::endl;
      uint i = 0;
      for (const auto &v : subhap.path) {
        std::string seq = get_gfa_seq(gbz, v >> 1);
        if (i == 0) {
          seq = seq.substr(seq.size() - subhap.offset1 - 1, seq.size());
        } else if (i == subhap.path.size() - 1) {
          seq = seq.substr(0, seq.size() - subhap.offset2);
        }
        std::cout << seq;
        ++i;
      }
      std::cout << std::endl;
    }
  } else {
    // get segments and links

    /* we cannot encode node identifiers or handles since vertices longer
     * than 1024 are broken down. Only way I found to store each vertex
     * only once (one element per "real" vertex) is by using this pair
     * all other identifiers in edges and subhaps are encoded (has strand bit
     * since we need it in GFA) */
    std::set<std::pair<std::string, std::pair<uint32_t, uint32_t>>> VV;
    std::set<uint64_t> EE;

    for (const subhap_t &subhap : subhaps) {
      VV.insert(
          gbz.graph.get_segment(gbz.graph.node_to_handle(subhap.path[0])));

      for (uint j = 1; j < subhap.path.size(); ++j) {
        VV.insert(
            gbz.graph.get_segment(gbz.graph.node_to_handle(subhap.path[j])));
        EE.insert(((uint64_t)subhap.path[j - 1] << 32) |
                  ((uint64_t)subhap.path[j]));
      }
    }

    // Print GFA
    for (const std::pair<std::string, std::pair<uint32_t, uint32_t>> &v : VV) {
      std::cout << "S\t" << v.first << "\t";
      for (size_t i = v.second.first; i < v.second.second; ++i) {
        std::cout << gbz.graph.get_sequence(gbz.graph.get_handle(i));
      }
      std::cout << std::endl;
    }

    for (const auto &ee : EE) {
      uint32_t v1 = ee >> 32;
      uint32_t v2 = (uint32_t)ee;
      std::cout << "L\t" << get_gfa_idx(gbz, gbwt::Node::id(v1)) << "\t"
                << (gbwt::Node::is_reverse(v1) ? '-' : '+') << "\t"
                << get_gfa_idx(gbz, gbwt::Node::id(v2)) << "\t"
                << (gbwt::Node::is_reverse(v2) ? '-' : '+') << "\t" << "0M"
                << std::endl;
    }

    for (const subhap_t &subhap : subhaps) {
      std::cout << "P" << "\t" << subhap.name << "#" << subhap.chrom << "\t"
                << get_gfa_idx(gbz, gbwt::Node::id(subhap.path[0]))
                << (gbwt::Node::is_reverse(subhap.path[0]) ? '-' : '+');
      // this works since we are sure to have at least one vertex
      for (uint j = 1; j < subhap.path.size(); ++j)
        std::cout << "," << get_gfa_idx(gbz, gbwt::Node::id(subhap.path[j]))
                  << (gbwt::Node::is_reverse(subhap.path[j]) ? '-' : '+');
      std::cout << "\t*" << std::endl;
    }
  }
  fprintf(stderr, "[I::%s] dumped %s in %.3f secs\n", __func__,
          fasta ? "FASTA" : "GFA", realtime() - rt);

  return 0;
}
