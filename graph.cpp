#include "graph.hpp"

Graph::Graph(const std::string &fp, bool gbz) : fp(fp) { is_gbz = gbz; }

int Graph::load() {
  sdsl::simple_sds::load_from(gbz, fp);

  try {
    sdsl::simple_sds::load_from(sampled_horders, fp + ".ho");
    has_horders = true;
  } catch (const sdsl::simple_sds::CannotOpenFile &) {
    has_horders = false;
  };
  return 0;
}

// TODO: how to report the strand of the vertex? do we need this info?

std::map<gbwt::size_type, std::vector<gbwt::size_type>>
Graph::locate(gbwt::node_type v) const {
  std::map<gbwt::size_type, std::vector<gbwt::size_type>> result;

  v = gbwt::Node::encode(v, 0);
  gbwt::SearchState searchstate = gbz.index.find(v);
  for (gbwt::size_type offset = searchstate.range.first;
       offset <= searchstate.range.second; ++offset) {
    gbwt::edge_type edge(searchstate.node, offset);

    // Move to next samples position
    int jumps = 0;
    while (gbz.index.tryLocate(edge) == gbwt::invalid_sequence()) {
      edge = gbz.index.LF(edge);
      ++jumps;
    }

    gbwt::size_type record_start =
        gbz.index.da_samples.start(gbz.index.toComp(edge.first));
    auto iter = gbz.index.da_samples.sampled_offsets.predecessor(record_start +
                                                                 edge.second);
    assert(iter->second == record_start + edge.second);
    // iter->first is the position of the sampled value

    uint sampled_ph = sampled_horders[iter->first];
    uint ph;
    ph = (sampled_ph & 1) ? (sampled_ph >> 1) - jumps
                          : (sampled_ph >> 1) + jumps;

    result[gbz.index.da_samples.array[iter->first]].push_back(ph);
  }

  // v = gbwt::Node::encode(v, 1);
  // searchstate = gbz.index.find(v);
  // for (gbwt::size_type offset = searchstate.range.first;
  //      offset <= searchstate.range.second; ++offset) {
  //   gbwt::edge_type edge(searchstate.node, offset);

  //   // Move to next samples position
  //   int jumps = 0;
  //   while (gbz.index.tryLocate(edge) == gbwt::invalid_sequence()) {
  //     edge = gbz.index.LF(edge);
  //     ++jumps;
  //   }

  //   gbwt::size_type record_start =
  //       gbz.index.da_samples.start(gbz.index.toComp(edge.first));
  //   auto iter = gbz.index.da_samples.sampled_offsets.predecessor(record_start
  //   +
  //                                                                edge.second);
  //   assert(iter->second == record_start + edge.second);
  //   // iter->first is the position of the sampled value

  //   uint sampled_ph = sampled_horders[iter->first];
  //   uint ph;
  //   ph = (sampled_ph & 1) ? (sampled_ph >> 1) - jumps
  //                         : (sampled_ph >> 1) + jumps;

  //   result[gbz.index.da_samples.array[iter->first]].push_back(ph);
  // }

  return result;
}

int Graph::serialize_horders() const {
  std::ofstream out;
  out.open(fp + ".ho", std::ofstream::out | std::ofstream::app);
  sampled_horders.simple_sds_serialize(out);
  out.close();

  return 0;
}

int Graph::build_horders() {
  sampled_horders =
      sdsl::int_vector<0>(get_nsamples(), 0, 32); // FIXME: is 16 enough?

  // std::string sample_name, contig_name;
  for (size_t i = 0; i < gbz.index.metadata.path_names.size(); ++i) {
    // sample_name = gbz.index.metadata.fullPath(i).sample_name;
    // contig_name = gbz.index.metadata.fullPath(i).contig_name;

    // TODO: can we get the - inserting position from the + position?
    for (int is_reverse = 0; is_reverse < 2; ++is_reverse) {
      gbwt::size_type path_idx = gbwt::Path::encode(i, is_reverse);
      gbwt::vector_type path = gbz.index.extract(path_idx);

      gbwt::node_type v = path[0];
      gbwt::SearchState ss = gbz.index.find(v);
      gbwt::size_type offset;
      for (offset = ss.range.first; offset <= ss.range.second; ++offset) {
        gbwt::size_type locate = gbz.index.locate(ss.node, offset);
        if (path_idx == locate)
          // we found the record offset we need
          break;
      }
      assert(offset <= ss.range.second);

      int n = 0;
      gbwt::edge_type edge(ss.node, offset);
      while (edge.first != gbwt::ENDMARKER) {
        if (gbz.index.tryLocate(edge) != gbwt::invalid_sequence()) {
          // position is sampled, so we get the index where to insert this
          // sample
          gbwt::size_type record_start =
              gbz.index.da_samples.start(gbz.index.toComp(edge.first));
          auto iter = gbz.index.da_samples.sampled_offsets.predecessor(
              record_start + edge.second);
          assert(iter->second == record_start + edge.second);
          // iter->first is the position to store
          assert(iter->first < get_nsamples());

          sampled_horders[iter->first] =
              !is_reverse ? ((n << 1) | 1) : ((path.size() - n - 1) << 1);
          // last bit stores if we are on + or - strand so that we later know if
          // we have to add or subtract the jumps. Strand here is
          // std::cerr << ">" << gbwt::Node::id(edge.first) << "/" << edge.first
          //           << ":" << edge.second << " " << iter->first << "="
          //           << sampled_horders[iter->first] << (is_reverse ? "+" :
          //           "-")
          //           << std::endl;
        }

        edge = gbz.index.LF(edge);
        ++n;
      }
    }
  }
  // std::cerr << gbz.index.da_samples.array << std::endl;
  // std::cerr << sampled_horders << std::endl;

  return 0;
}

size_t Graph::get_nsamples() const { return gbz.index.da_samples.size(); }

gbwt::size_type Graph::get_horder(gbwt::size_type p, gbwt::node_type v,
                                  bool first) const {
  bool hit = false;
  gbwt::size_type res = first ? -1 : 0;

  gbwt::SearchState ss = gbz.index.find(v);
  for (gbwt::size_type offset = ss.range.first; offset <= ss.range.second;
       ++offset) {
    gbwt::edge_type edge(ss.node, offset);

    // Move to next samples position
    int jumps = 0;
    while (gbz.index.tryLocate(edge) == gbwt::invalid_sequence()) {
      edge = gbz.index.LF(edge);
      ++jumps;
    }
    gbwt::size_type record_start =
        gbz.index.da_samples.start(gbz.index.toComp(edge.first));
    auto iter = gbz.index.da_samples.sampled_offsets.predecessor(record_start +
                                                                 edge.second);
    if (gbz.index.da_samples.array[iter->first] != p)
      continue;

    assert(iter->second == record_start + edge.second);
    // iter->first is the position of the sampled value

    uint sampled_ph = sampled_horders[iter->first];
    gbwt::size_type ph = (sampled_ph & 1) ? (sampled_ph >> 1) - jumps
                                          : (sampled_ph >> 1) + jumps;
    if (first)
      res = res < ph ? res : ph;
    else
      res = res < ph ? ph : res;
    hit = true;
  }
  return hit ? res : -1;
}

std::string Graph::get_gfa_idx(gbwt::node_type v) const {
  return gbz.graph.get_segment_name(gbz.graph.get_handle(v));
}
