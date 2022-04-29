#include "misc.hpp"
#include <map>

using std::map;

namespace gfase{

string join(vector <string> s, char delimiter){
    string joined_string;

    for (size_t i=0; i<s.size(); i++){
        if (i < s.size()-1){
            joined_string += s[i] + delimiter;
        }
        else{
            joined_string += s[i];
        }
    }

    return joined_string;
}


void run_command(string& command){
    int exit_code = system(command.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + command);
    }
}


path align(path output_dir, path ref_path, path query_path, size_t n_threads){
    path output_path;
    string prefix = query_path.filename().replace_extension("").string();
    string suffix = ref_path.filename().replace_extension("").string();

    replace(prefix.begin(), prefix.end(), '.', '_');
    replace(suffix.begin(), suffix.end(), '.', '_');

    output_path = prefix + "_VS_" + suffix + ".sam";
    output_path = absolute(output_dir) / output_path;

    // Set up arguments
    vector<string> arguments = {"minimap2",
                                "-a",
                                "-x", "asm20",
                                "-K", "10g",                    // New parameter for large batch size, better cpu %
                                "--eqx",
                                "-t", to_string(n_threads),
                                ref_path.string(),
                                query_path.string(),
                                "-o", output_path.string()
    };

    // Convert arguments to single string
    string command = join(arguments);
    cerr << "\nRUNNING: " << command << "\n";
    run_command(command);

    return output_path;
}


path sam_to_sorted_bam(path sam_path, size_t n_threads, bool remove_sam){

    vector<string> args = {
            "samtools", "sort",
            "-@", to_string(n_threads),
            "-o", "",
            ""
    };

    auto out_path = sam_path;
    out_path.replace_extension(".sorted.bam");
    args[5] = out_path;
    args[6] = sam_path;

    auto command = join(args);

    cerr << "\nRUNNING: " << command << "\n";
    run_command(command);

    if (remove_sam){
        remove(sam_path);
    }

    return out_path;
}


void get_query_lengths_from_fasta(path fasta_path, map<string,size_t>& query_lengths){
    ifstream file(fasta_path);

    string name;
    char c;

    while (file.get(c)){
        if (c == '>'){
            name.clear();
            while (file.get(c)){
                if (std::isspace(c)){
                    query_lengths.emplace(name,0);
                    break;
                }
                else {
                    name += c;
                }
            }
        }
        else if (c != '\n'){
            query_lengths.at(name)++;
        }
    }
}


void bin_fasta_sequences(path fasta_path, const unordered_map<string,bool>& phased_contigs){
    ifstream file(fasta_path);

    string name;
    string sequence;

    // 0=pat, 1=mat, 2=both
    int bin;
    char c;

    // TODO: in progress

    while (file.get(c)){
        if (c == '>'){
            name.clear();
            while (file.get(c)){
                if (std::isspace(c)){
                    // Check if this segment is phased, set flag accordingly
                    auto result = phased_contigs.find(name);

                    if (result != phased_contigs.end()){
                        bin = int(result->second);
                    }
                    else{
                        bin = 2;
                    }
                    break;
                }
                else {
                    name += c;
                }
            }
        }
        else if (c != '\n'){
            query_lengths.at(name)++;
        }
    }
}


}