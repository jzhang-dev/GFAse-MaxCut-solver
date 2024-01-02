#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "optimize.hpp"
#include "CLI11.hpp"

using gfase::NonBipartiteEdgeException;
using gfase::IncrementalIdMap;
using gfase::MultiContactGraph;
using gfase::alt_component_t;
using ghc::filesystem::path;
using CLI::App;

#include <iostream>

using std::exception;
using std::cerr;

int main(int argc, char* argv[]){

    path id_path;
    path graph_path;
    path output_dir;
    size_t n_threads = 1;
    size_t core_iterations = 200;
    size_t sample_size = 30;
    size_t n_rounds = 2;
    

    CLI::App app{"App description"};
    app.add_option(
        "-i,--id_path",
        id_path,
        "")
        ->required();
    app.add_option(
        "-g,--graph_path",
        graph_path,
        "")
        ->required();
    app.add_option(
        "-o,--output_dir",
        output_dir,
        "")
        ->required();
    app.add_option(
            "-c,--core_iterations",
            core_iterations,
            "(Default = "+ to_string(core_iterations) + ")\tNumber of iterations to use for each shallow convergence in the sampling process. The final phasing round uses 3*core_iterations.");

    app.add_option(
            "-s,--sample_size",
            sample_size,
            "(Default = "+ to_string(sample_size) + ")\tHow many shallowly converged phase states to sample from. This is also the maximum usable concurrency (n_threads) for this stage of the pipeline.");

    app.add_option(
            "-r,--n_rounds",
            n_rounds,
            "(Default = " + to_string(n_rounds) + ")\tHow many rounds to sample and merge.");
    app.add_option(
            "-t,--threads",
            n_threads,
            "(Default = " + to_string(n_threads) + ")\tMaximum number of threads to use.");
    CLI11_PARSE(app, argc, argv);

    cerr << "Load ID map" << '\n';
    IncrementalIdMap<string> id_map(id_path);

    cerr << "Load graph" << '\n';
    MultiContactGraph contact_graph(graph_path, id_map);

    cerr << "Infer alts from Shasta names" << '\n';
    contact_graph.get_alts_from_shasta_names(id_map);

    cerr << "Remove nodes that don't have any involvement in bubbles" << '\n';
    vector<int32_t> to_be_deleted;
    contact_graph.for_each_node([&](int32_t id){
        if (not contact_graph.has_alt(id)){
            to_be_deleted.emplace_back(id);
        }
    });

    for (auto& id: to_be_deleted){
        contact_graph.remove_node(id);
    }

    cerr << "Remove self edges in contact graph" << '\n';
    contact_graph.for_each_node([&](int32_t id) {
        contact_graph.remove_edge(id,id);
    });

    if (contact_graph.edge_count() == 0){
        throw runtime_error("ERROR: no inter-contig contacts detected in alignments, no usable phasing information");
    }

    cerr << "Optimizing phases..." << '\n';
    
    monte_carlo_phase_contacts(
            contact_graph,
            id_map,
            core_iterations,
            sample_size,
            n_rounds,
            n_threads,
            output_dir);

    return 0;
}