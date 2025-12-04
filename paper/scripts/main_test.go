package main

import(
    "GoAssembly/pkg/assembly"
    "log"
    "math/rand"
    "os"
    "path/filepath"
    "testing"
)

// Always run benchmarks with the default #workers and buffer size and the
// "shortest" algorithm variant, which is the only one implemented.
var workers = 100
var bufsize = 100
var variant = "shortest"

func BenchmarkDatasets(b *testing.B) {
    // Define the list of datasets to include in this benchmark.
    // The filepaths are relative to paper/assembly_go/cmd/app/main_test.go.
    datasets := []struct {
        name string
        path string
    }{
        {"gdb13_1201", "../../../../data/gdb13_1201"},
        {"gdb17_200", "../../../../data/gdb17_200"},
        {"checks", "../../../../data/checks"},
		// {"coconut_55", "../../../../data/coconut_55"},  // SLOW! Use -count=1
    }

    // Run all sub-benchmarks.
    for _, dataset := range datasets {
        b.Run(dataset.name, func(b *testing.B) {
            // Collect all .mol filepaths in the dataset.
            entries, err := os.ReadDir(dataset.path)
            if err != nil {
                log.Fatal(err)
            }
            var mol_files []string
            for _, entry := range entries {
                if entry.IsDir() {
                    continue
                }
                if filepath.Ext(entry.Name()) == ".mol" {
                    mol_path := filepath.Join(dataset.path, entry.Name())
                    mol_files = append(mol_files, mol_path)
                }
            }

            // Make the graph data structure for each .mol file.
            var graphs []assembly.Graph
            for _, mol_file := range mol_files {
                graphs = append(graphs, assembly.MolColourGraph(mol_file))
            }

            // Shuffle the list of graphs for better statistics.
            rand.Shuffle(len(graphs), func(i, j int) {
                graphs[i], graphs[j] = graphs[j], graphs[i]
            })

            // Reset timer and start actual sub-benchmark.
            b.ResetTimer()
            for i := 0; i < b.N; i++ {
                for _, graph := range graphs {
                    var pathways []assembly.Pathway
                    pathways = assembly.Assembly(graph, workers, bufsize, variant)
                    _ = assembly.AssemblyIndex(&pathways[0], &graph)
                }
            }
        })
    }
}
