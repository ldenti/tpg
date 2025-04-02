#ifndef TPG_GRAPH_HPP
#define TPG_GRAPH_HPP

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gbwt/gbwt.h"
#include "gbwtgraph/gbz.h"
#include "sdsl/simple_sds.hpp"

class Graph {
public:
  std::string fp;
  bool is_gbz;

  // gbz-based
  gbwtgraph::GBZ gbz;
  bool has_horders;
  sdsl::int_vector<0> sampled_horders;

public:
  Graph(const std::string &fp, bool is_gfa);
  int load();

  int build_horders();
  int serialize_horders() const;

  std::map<gbwt::size_type, std::vector<gbwt::size_type>>
  locate(gbwt::node_type v) const;
  gbwt::size_type get_horder(gbwt::size_type p, gbwt::node_type v,
                             bool first) const;
  size_t get_nsamples() const;
  std::string get_gfa_idx(gbwt::node_type v) const;
};

#endif
