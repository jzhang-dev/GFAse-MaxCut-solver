#include "ContactGraph.hpp"

#include <ostream>

using std::ofstream;
using std::ostream;
using std::cerr;
using std::min;
using std::max;


namespace gfase{


Node::Node(int8_t partition):
        neighbors(),
        partition(partition),
        alt(-1)
{}


bool Node::has_alt() const{
    return alt >= 0;
}


const array<string,3> ContactGraph::colors = {"Cornflower Blue", "Plum", "Tomato"};


ContactGraph::ContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map) {
    string s;

    for (const auto& [a,sub_map]: contact_map){
        try_insert_node(a);

        for (const auto& [b,count]: sub_map){
            try_insert_node(b);
            insert_edge(a,b,count);
        }
    }
}


pair<int32_t,int32_t> edge(int32_t a, int32_t b){
    return {min(a,b), max(a,b)};
}


void ContactGraph::insert_edge(int32_t a, int32_t b, int32_t weight){
    edge_weights.emplace(edge(a,b),weight);
    nodes.at(a).neighbors.emplace(b);
    nodes.at(b).neighbors.emplace(a);
}


void ContactGraph::try_insert_edge(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,0);
    }
}


void ContactGraph::add_alt(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add alt with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto& node_a = nodes.at(a);
    auto& node_b = nodes.at(b);

    // Remove any edge that may exist between a and b
    remove_edge(a,b);

    node_a.alt = b;
    node_b.alt = a;

    // Chose arbitrary partition assignment
    node_a.partition = -1;
    node_b.partition = 1;
}


void ContactGraph::try_insert_edge(int32_t a, int32_t b, int32_t weight){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,weight);
    }
}


void ContactGraph::increment_edge_weight(int32_t a, int32_t b, int32_t value){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        result->second += value;
    }
}


void ContactGraph::remove_edge(int32_t a, int32_t b){
    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        edge_weights.erase(result);
        nodes.at(a).neighbors.erase(b);
        nodes.at(b).neighbors.erase(a);
    }
}


void ContactGraph::for_each_node_neighbor(int32_t id, const function<void(int32_t id_other, const Node& n)>& f) const{
    for (auto& id_other: nodes.at(id).neighbors){
        f(id_other, nodes.at(id_other));
    }
}


void ContactGraph::for_each_node(const function<void(int32_t id, const Node& n)>& f) const{
    for (auto& [id,node]: nodes){
        f(id, nodes.at(id));
    }
}


void ContactGraph::for_each_edge(const function<void(const pair<int32_t,int32_t> edge, int32_t weight)>& f) const{
    for (auto& [e, weight]: edge_weights){
        f(e, weight);
    }
}


void ContactGraph::insert_node(int32_t id, int8_t partition){
    if (partition < -1 or partition > 1){
        throw runtime_error("ERROR: can't assign partition index outside of {-1,0,1}");
    }

    nodes.emplace(id, partition);
}


void ContactGraph::insert_node(int32_t id){
    nodes.emplace(id, 0);
}


void ContactGraph::try_insert_node(int32_t id){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, 0);
    }
}


void ContactGraph::try_insert_node(int32_t id, int8_t partition){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, partition);
    }
}


void ContactGraph::set_partition(int32_t id, int8_t partition) {
    auto& node = nodes.at(id);
    node.partition = partition;

    // If this node is linked to an alt, alt must be maintained in an opposite state
    if (node.has_alt()){
        if (partition == 0) {
            throw runtime_error("ERROR: cannot set 0 partition for bubble: " + to_string(id));
        }

        nodes.at(node.alt).partition *= -1;
    }
}


void ContactGraph::remove_node(int32_t id){
    vector <pair <int32_t, int32_t> > to_be_removed;

    for_each_node_neighbor(id, [&](int32_t id_other, const Node& n){
        to_be_removed.emplace_back(edge(id, id_other));
    });

    for (auto& e: to_be_removed){
        remove_edge(e.first, e.second);
    }

    // Make sure there is no dangling reference to this node in its alt
    auto& n = nodes.at(id);
    if (n.has_alt()) {
        auto alt_id = n.alt;
        nodes.at(alt_id).alt = -1;
    }

    nodes.erase(id);
}


size_t ContactGraph::edge_count(int32_t id){
    return nodes.at(id).neighbors.size();
}


size_t ContactGraph::size(){
    return nodes.size();
}


int64_t ContactGraph::get_score(const Node& a, const Node& b, int32_t weight) const{
    int64_t score = 0;

    auto p_a = a.partition;
    auto p_b = b.partition;

    if (p_a != 0 and p_b != 0) {
        score = p_a * p_b * weight;
    }
    else{
        score = weight/2;
    }

    return score;
}


int64_t ContactGraph::compute_consistency_score(int32_t id) const{
    int64_t score = 0;

    auto n = nodes.at(id);

    for_each_node_neighbor(id, [&](int32_t id_other, const Node& n_other){
        score += get_score(n, n_other, edge_weights.at(edge(id, id_other)));
//        cerr << '\t' << id << "<->" << id_other << ' ' << int(n.partition) << 'x' << int(n_other.partition) << 'x' << edge_weights.at(edge(id, id_other)) << ' ' << score << '\n';
    });

    if (n.has_alt()){
        auto& n_alt = nodes.at(n.alt);
        for_each_node_neighbor(n.alt, [&](int32_t id_other, const Node& n_other){
            score += get_score(n_alt, n_other, edge_weights.at(edge(n.alt, id_other)));
        });
    }

    return score;
}


int64_t ContactGraph::compute_total_consistency_score() const{
    int64_t score = 0;

    for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto& a = nodes.at(edge.first);
        auto& b = nodes.at(edge.second);
        score += get_score(a, b, weight);
    });

    return score;
}


void ContactGraph::get_partitions(vector <pair <int32_t,int8_t> >& partitions) const{
    partitions.clear();
    partitions.resize(nodes.size());

    size_t i=0;
    for (auto& [n, node]: nodes){
        partitions[i] = {n,node.partition};
        i++;
    }
}


void ContactGraph::randomize_partitions(){
    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,2);

    for (auto& [id,node]: nodes){
        int8_t p;
        if (node.has_alt()){
            // Only allow {1,-1} for known bubbles
            p = int8_t((uniform_distribution(rng) % 2));

            if (p == 0){
                p = -1;
            }

            node.partition = p;
            nodes.at(node.alt).partition = -p;
        }
        else{
            // Allow {1,0,-1}
            p = int8_t((uniform_distribution(rng) % 3) - 1);

            node.partition = p;
        }
    }
}


void ContactGraph::set_partitions(const vector <pair <int32_t,int8_t> >& partitions){
    for (auto& [n, p]: partitions){
        set_partition(n, p);
    }
}


void ContactGraph::get_node_ids(vector<int32_t>& ids){
    ids.resize(nodes.size());

    size_t i = 0;

    for (auto& [id,node]: nodes){
        ids[i] = id;
        i++;
    }
}


bool ContactGraph::has_alt(int32_t id) const{
    return nodes.at(id).has_alt();
}


bool ContactGraph::has_node(int32_t id) const{
    return nodes.find(id) != nodes.end();
}


ostream& operator<<(ostream& o, const Node& n){
    o << '\t' << "partition: " << int(n.partition) << '\n';
    o << '\t' << "neighbors: ";

    for (auto& id: n.neighbors){
        o << id << ' ';
    }

    return o;
}


void ContactGraph::write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map){
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Phase" << ',' << "Color" << '\n';

    for (auto& [id,node]: nodes){
        auto name = id_map.get_name(id);
        file << name << ',' << int(node.partition) << ',' << colors[node.partition+1] << '\n';
    }
}


void random_phase_search(
        ContactGraph contact_graph,
        const vector<int32_t>& ids,
        vector <pair <int32_t,int8_t> >& best_partitions,
        atomic<int64_t>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(contact_graph.size()-1));

    int64_t total_score;

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i<((contact_graph.size()/30) + 1); i++) {
            auto r = ids[uniform_distribution(rng)];

            int8_t p;
            if (contact_graph.has_alt(r)){
                // Only allow {1,-1}
                p = int8_t((uniform_distribution(rng) % 2));

                if (p == 0){
                    p = -1;
                }
            }
            else{
                // Allow {1,0,-1}
                p = int8_t((uniform_distribution(rng) % 3) - 1);
            }

            contact_graph.set_partition(r, p);
        }

        for (size_t i=0; i<contact_graph.size(); i++) {
            auto n = ids[uniform_distribution(rng)];

            int64_t max_score = std::numeric_limits<int64_t>::min();
            int8_t p_max = 0;

            if (contact_graph.edge_count(n) == 0){
                continue;
            }

            bool has_alt = contact_graph.has_alt(n);

            for (int8_t p=-1; p<=1; p++) {
                // If the node has an "alt" it can't be made neutral
                if (has_alt and p==0){
                    continue;
                }

                contact_graph.set_partition(n, p);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = p;
                }

//                cerr << n << ' ' << int(has_alt) << ' ' << int(p) << ' ' << score << ' ' << max_score << '\n';
            }

            contact_graph.set_partition(n, p_max);
        }

        total_score = contact_graph.compute_total_consistency_score();

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            contact_graph.get_partitions(best_partitions);
        }
        else {
            contact_graph.set_partitions(best_partitions);
        }

        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << std::flush;
//        for (auto& [n,p]: best_partitions){
//            cerr << '(' << n << ',' << int(p) << ") ";
//        }
        cerr << '\n';

        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }

    contact_graph.set_partitions(best_partitions);
}



}