#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <iostream>

#include "improvedBnB.h"
#include "molfileParser.h"
#include "molGraph.h"

// Define benchmark parameters.
static uint16_t WARMUP_SECS = 3;
static uint16_t NUM_SAMPLES = 20;

// Parse all .mol files in the given dataset as molGraphs.
std::vector<molGraph> load_dataset_molecules(std::string& dataset) {
    std::vector<molGraph> mol_graphs;
    for (const auto& entry : std::filesystem::directory_iterator(dataset)) {
        if (entry.path().extension() == ".mol") {
            std::ifstream mol_fs(entry.path().string());
            molGraph mol_graph;
            molfileParser(mol_fs, mol_graph);
            mol_graphs.push_back(mol_graph);
        }
    }

    return mol_graphs;
}

int main(int argc, char** argv) {
    // Define the reference datasets and a container for results.
    std::vector<std::string> datasets = {
        "../../data/gdb13_1201",
        "../../data/gdb17_200",
        "../../data/checks",
        "../../data/coconut_55",
    };
    std::map<std::string, std::vector<std::chrono::nanoseconds>> results;
    
    // Run the benchmark on each reference dataset.
    for (auto& dataset : datasets) {
        // Load all molecules in the given reference dataset.
        auto mol_graphs = load_dataset_molecules(dataset);
        std::cout << "Loaded " << mol_graphs.size() << " molecules from "
            << dataset << std::endl;
        
        // Do a CPU warmup period of calculations where times aren't recorded.
        std::cout << "Warming up for " << WARMUP_SECS << " seconds..."
            << std::endl;
        auto t_start = std::chrono::steady_clock::now();
        do {
            for (auto& mol_graph : mol_graphs) {
                improvedBnB(mol_graph);
            }
        } while (std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::steady_clock::now() - t_start)
                 < std::chrono::seconds(WARMUP_SECS));
        
        // Now that things are warm, do the actual benchmark.
        for (uint16_t s = 0; s < NUM_SAMPLES; s++) {
            std::cout << "\rBenchmarking " << (s + 1) << " of " << NUM_SAMPLES
                << " samples..." << std::flush;
            auto t_start = std::chrono::steady_clock::now();
            for (auto& mol_graph : mol_graphs) {
                improvedBnB(mol_graph);
            }
            auto t_end = std::chrono::steady_clock::now();
            results[dataset].push_back(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        t_end - t_start));
        }
        std::cout << std::endl;
    }

    // TODO: Dump results into a CSV.
    
    return 0;
}
