
#include "molfileParser.h"
#include <stddef.h>           // for size_t
#include <fstream>            // for operator<<, basic_ostream::operator<<
#include <iostream>           // for cout
#include <string>             // for allocator, char_traits, getline, stoi
#include <vector>             // for vector
#include <sstream>            // std::istringstream
#include "globalPrimitives.h" // for coords, removeHydrogens
#include "molGraph.h"         // for molGraph, atom

using namespace std;

string molfileParser(ifstream &molfile, molGraph &mg)
{
    string name, currLine;
    double x; // placeholder to store coordinates which are not necessary
    int totalAtoms = 0, totalBonds = 0;

    getline(molfile, name);
    for (int i = 0; i < 2; i++)
    {
        getline(molfile, currLine);
    }

    getline(molfile, currLine);
    istringstream iss(currLine);
    string molLine = iss.str(), atoms = molLine.substr(0, 3), bonds = molLine.substr(3, 3);
    totalAtoms = stoi(atoms);
    totalBonds = stoi(bonds);

    for (int i = 0; i < totalAtoms; i++)
    {
        getline(molfile, currLine);
        istringstream iss(currLine);
        for (int i = 0; i < 3; i++)
        {
            iss >> x;
            coords.push_back(x);
        }
        string s;
        iss >> s;
        mg.addAtom(s);
    }
    for (int i = 0; i < totalBonds; i++)
    {
        getline(molfile, currLine);
        istringstream iss(currLine);
        int atomA, atomB, bondOrder;
        string molLine = iss.str(), atomAs = molLine.substr(0, 3),
               atomBs = molLine.substr(3, 3), bondOrderS = molLine.substr(6, 3);
        atomA = stoi(atomAs), atomB = stoi(atomBs), bondOrder = stoi(bondOrderS);
        mg.addBond(atomA - 1, atomB - 1, bondOrder);
        string temp;
        while (iss)
        {
            string s0;
            iss >> s0;
            temp += s0 + ' ';
        }
    }
    if (removeHydrogens)
    {
        for (size_t i = 0; i < mg.mg.size(); i++)
        {
            if (mg.atype(i) == "H")
                mg.removeAtom(i);
        }
        mg.removeAndCollapse();
    }
    return name;
}
