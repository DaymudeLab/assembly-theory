
#include "assemblyCalculator.h"
#include <chrono>
#include <iostream>           // for operator<<, basic_ostream, cout, ifstream
#include <string>             // for char_traits, allocator, operator+, ope...
#include <vector>             // for vector
#include "globalPrimitives.h" // for moleculeName
#include "graphio.h"          // for graphio
#include "improvedBnB.h"      // for improvedBnB
#include "molGraph.h"         // for molGraph
#include "molfileParser.h"    // for molfileParser

using namespace std;

void assemblyCalculator(string &fileName)
{
    string in = fileName + ".mol";
    ifstream molfile(in.c_str());
    if (molfile.is_open())
    {
        ofstream outputFile;
        molGraph mol_graph;
        vector<double> coords;
        molfileParser(molfile, mol_graph);
        auto t_start = chrono::steady_clock::now();
        improvedBnB(mol_graph, outputFile);
        auto t_end = chrono::steady_clock::now();
        cout << chrono::duration_cast<chrono::nanoseconds>(t_end - t_start).count();
    }
    else
    {
        ifstream graphFile(fileName.c_str());
        if (graphFile.is_open())
        {
            ofstream outputFile;
            molGraph mol_graph;
            graphio(graphFile, mol_graph);
            improvedBnB(mol_graph, outputFile);
        }
        else
            cout << "No file found\n";
    }
}
