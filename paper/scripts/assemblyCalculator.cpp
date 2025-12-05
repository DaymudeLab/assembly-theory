
#include "assemblyCalculator.h"
#include <time.h>             // for clock, clock_t
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
        moleculeName = fileName + "Pathway";
        molGraph mol_graph;
        vector<double> coords;
        molfileParser(molfile, mol_graph);
        clock_t t1 = clock();
        improvedBnB(mol_graph);
        clock_t t2 = clock();
    }
    else
    {
        ifstream graphFile(fileName.c_str());
        if (graphFile.is_open())
        {
            molGraph mol_graph;
            graphio(graphFile, mol_graph);
            moleculeName = fileName + "Pathway";
            clock_t startTime = clock();
            improvedBnB(mol_graph);
        }
        else
            cout << "No file found\n";
    }
}
