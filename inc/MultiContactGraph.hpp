#ifndef GFASE_MULTICONTACTGRAPH_HPP
#define GFASE_MULTICONTACTGRAPH_HPP

#include "IncrementalIdMap.hpp"
#include "BubbleGraph.hpp"          // TODO: Only need `contact_map_t` definition, consider breaking it out

#include "handlegraph/handle_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "sparsepp/spp.h"
#include "Filesystem.hpp"
#include "edge.hpp"

using ghc::filesystem::path;
using spp::sparse_hash_map;
using handlegraph::handle_t;
using handlegraph::edge_t;
using handlegraph::nid_t;
using bdsg::HashGraph;


#include <functional>
#include <utility>
#include <thread>
#include <ostream>
#include <array>
#include <set>
#include <map>
#include <set>

using std::function;
using std::atomic;
using std::mutex;
using std::pair;
using std::ostream;
using std::array;
using std::map;
using std::set;


namespace gfase{


using alt_component_t = pair <set<int32_t>, set<int32_t> >;


class NonBipartiteEdgeException : public std::runtime_error{
public:
    alt_component_t component;
    int32_t a;
    int32_t b;
    string message;

    const char* what() const noexcept override;

    NonBipartiteEdgeException(alt_component_t& c, int32_t a, int32_t b);
};


class MultiNode {
public:
    // To find adjacent nodes
    set<int32_t> neighbors;

    // Total reads located on this node > mapQ threshold, regardless of pair mapQ
    int64_t coverage;

    // Sequence length of this node
    int32_t length;

    // To find linked/opposing node in a bubble
    set<int32_t> alts;

    // Which set does this node belong to
    int8_t partition;

    MultiNode(int8_t partition);
    bool has_alt() const;
};


class MultiContactGraph {
    // Edge map which will only store edges in sorted order {min(a,b), max(a,b)}
    // This is maintained in duplicate with node-level adjacency because there needs to be a 1x fast iteration method
    // to compute total consistency scores
    sparse_hash_map<pair<int32_t,int32_t>, int32_t> edge_weights;
    sparse_hash_map<int32_t,MultiNode> nodes;

    static const array<string,3> colors;

    // No safety checks built in, only should be called when it's known that the nodes exist and the edge does not.
    void insert_edge(int32_t a, int32_t b, int32_t weight);

public:
    // Constructors
    MultiContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map);
    MultiContactGraph();

    // Editing
    void remove_edge(int32_t a, int32_t b);
    void try_insert_edge(int32_t a, int32_t b);
    void try_insert_edge(int32_t a, int32_t b, int32_t weight);
    void increment_edge_weight(int32_t a, int32_t b, int32_t value);
    void increment_coverage(int32_t a, int64_t value);
    void set_node_coverage(int32_t a, int64_t value);
    void set_node_length(int32_t a, int32_t length);
    void insert_node(int32_t id);
    void insert_node(int32_t id, int8_t partition);
    void try_insert_node(int32_t id);
    void try_insert_node(int32_t id, int8_t partition);
    void remove_node(int32_t id);
    void set_partition(int32_t id, int8_t partition);
    void add_alt(int32_t a, int32_t b);
    size_t edge_count(int32_t id);

    // Iterating and accessing
    void for_each_node_neighbor(int32_t id, const function<void(int32_t id_other, const MultiNode& n)>& f) const;
    void for_each_node(const function<void(int32_t id, const MultiNode& n)>& f) const;
    void for_each_node(const function<void(int32_t id)>& f) const;
    void for_each_edge(const function<void(const pair<int32_t,int32_t>, int32_t weight)>& f) const;
    void for_each_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f);
    void for_each_double_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f);
    void get_node_ids(vector<int32_t>& ids);
    bool has_alt(int32_t id) const;
    bool has_node(int32_t id) const;
    bool has_edge(int32_t a, int32_t b) const;
    int64_t get_node_coverage(int32_t id) const;
    int32_t get_node_length(int32_t id) const;
    int32_t get_edge_weight(int32_t a, int32_t b) const;
    void get_alt_component(int32_t id, bool validate, alt_component_t& component) const;
    int8_t get_partition(int32_t id) const;

    // Optimization
    double get_score(const MultiNode& a, const MultiNode& b, int32_t weight) const;
    double compute_total_consistency_score() const;
    double compute_consistency_score(int32_t id) const;
    void get_partitions(vector <pair <int32_t,int8_t> >& partitions) const;
    void set_partitions(const vector <pair <int32_t,int8_t> >& partitions);
    void randomize_partitions();

    // Misc
    void write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map) const;
    void write_node_data(path output_path, IncrementalIdMap<string>& id_map) const;
    void assert_component_is_valid(const alt_component_t& component) const;
    void validate_alts();
    size_t size();
    void resize();
};

ostream& operator<<(ostream& o, const gfase::MultiNode& n);

void random_multicontact_phase_search(
        MultiContactGraph contact_graph,
        const vector<int32_t>& ids,
        vector <pair <int32_t,int8_t> >& best_partitions,
        atomic<double>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations);

}

#endif //GFASE_MULTICONTACTGRAPH_HPP
