/**
 * @file improvedBnB.h
 * @brief code relating to improved branch and bound algorithm
 */
#pragma once
#include <fstream>            // for ofstream
#include <map>                // for map
#include <vector>             // for vector
#include "globalPrimitives.h" // for standardBitset
struct assemblyState;
struct dagDuplicateSet;
struct initialDuplicateSet;
struct molGraph;

// using namespace std;

/**
 * @brief Enumerate all subgraphs during the initial phase of the pathway algorithm. See Seet et al section 4.3 Duplicate Enumeration
 *
 * @param _target The initial assembly state
 * @param stmapVector The matchings found
 * @return true if any matchings found
 * @return false if no matchings found
 */
bool initialRecursiveEnumeration(assemblyState &_target, std::vector<std::map<int, initialDuplicateSet>> &stmapVector, bool &earlyTerminate);

/**
 * @brief Enumerate all subgraphs during subsequent phases of the pathway algorithm using the DAG to speed things up.  See seet et al section 4.3 Duplicate Enumeration
 *
 * @param _target Target assembly state
 * @param stmapVector List of generated duplicate pairs
 * @return true if matches found
 * @return false otherwise
 */
int dagRecursiveEnumeration(assemblyState &_target, std::vector<std::map<int, dagDuplicateSet>> &stmapVector,
                            std::vector<std::vector<standardBitset>> &targetMasks);

/**
 * @brief This function returns the minimum possible sum of duplicate bonds that can be found if only fragments equal
 * to the size of the largest duplicatable subgraph are considered. Used to tighten the bound after fragmentation
 * @param target The target assembly state
 * @param matchMask The bitset of all graphs isomorphic to the matching
 * @param maxFragMask The bitset of all duplicatable subgraphs with the same bitset count as the matching
 * @return int the sum of duplicate bonds calculated by this cutoff function
 */
int postFragmentationCutoff(assemblyState &target, standardBitset &matchMask, standardBitset &maxFragMask);

/**
 * @brief The recursive function that enumerates duplicates and generates assembly states on all but the first pass
 * of the assembly algorithm
 *
 * @param input The input assembly state
 * @param AI The global minimum assembly index found
 * @return true if any more duplicatable substructures are found
 * @return false otherwise
 */
bool dagRecursiveAssembly(assemblyState &input, int &AI);

/**
 * @brief
 * @brief The recursive function that enumerates duplicates and generates assembly states on the first pass
 * of the assembly algorithm
 *
 * @param input The input assembly state
 * @param AI The global minimum assembly index found
 * @return true if any more duplicatable substructures are found
 * @return false otherwise
 */
bool initialRecursiveAssembly(assemblyState &input, int &AI);

/**
 * @brief Function that calls the recursive assembly function
 *
 * @param mg The target molGraph
 */
void improvedBnB(molGraph &mg);
