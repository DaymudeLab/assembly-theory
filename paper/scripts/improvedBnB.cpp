#include "improvedBnB.h"
#include <algorithm>           // for max
#include <bitset>              // for bitset, hash, operator&
#include <ctime>               // for clock, size_t
#include <fstream>             // for operator<<, char_traits, basic_ostream
#include <iostream>            // for cout
#include <map>                 // for map, operator!=, _Rb_tree_iterator
#include <unordered_map>       // for unordered_map, _Node_iterator
#include <unordered_set>       // for unordered_set
#include <utility>             // for pair
#include <vector>              // for vector
#include "assemblyState.h"     // for apWrapper, assemblyState, assemblyPath
#include "dagEnumeration.h"    // for convertDag
#include "duplicateMatching.h" // for dagDuplicateSet, initialDuplicateSet
#include "fragmentation.h"     // for fragmentAssemblyState, clearPathMap
#include "globalPrimitives.h"  // for standardBitset, bitsetHashTable, inte...
#include "graphHashes.h"       // for graphHash, canonise, graphHashMap
#include "molGraph.h"          // for molGraph, preprocessWriteback, target...
#include "pathwayGenerator.h"  // for recoverPathway2

using namespace std;

bool initialRecursiveEnumeration(assemblyState &_target, std::vector<std::map<int, initialDuplicateSet>> &stmapVector, bool &earlyTerminate)
{
    vector<std::unordered_map<standardBitset, pair<int, vector<standardBitset>>>> tempDag(2);
    int ordinal = MAX_INT;
    vector<standardBitset> &masks = _target.masks;
    bool alive = 0;
    size_t currSize = 1;
    int pr = 0;
    bool exitEnum = 0;
    vector<standardBitset> targetMask(masks.size(), 0);

    vector<initialPotentialDuplicate> *matchingList1ptr = new vector<initialPotentialDuplicate>;
    vector<initialPotentialDuplicate> &matchingList1 = *(matchingList1ptr);
    std::unordered_set<standardBitset> maskMap;
    
    // Find all one-bond duplicatable subgraphs
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (masks[i][j] != 0)
            {
                initialPotentialDuplicate m(j, masks[i], i);
                pair<int, vector<standardBitset>> p;
                p.first = -1;
                standardBitset b = 0;
                b.set(j);
                tempDag[0][b] = p;
                m.generateDAG(matchingList1, i, maskMap, tempDag);
            }
        }
    }

    bool active = 1, overweight = 0;
    vector<initialPotentialDuplicate> *currMLptr = nullptr, *prevMLptr = matchingList1ptr;
    while (active)
    {
        map<int, initialDuplicateSet> temp;
        stmapVector.push_back(temp);
        map<int, initialDuplicateSet> &stmap = stmapVector.back();
        active = 0;
        currMLptr = new vector<initialPotentialDuplicate>;
        vector<initialPotentialDuplicate> &currML = *currMLptr, &prevML = *prevMLptr;
        // Iterate through all matching subgraphs from the previous matchlist
        for (size_t i = 0; i < prevML.size(); i++)
        {
            if (clock() - startTime > runTimeMax)
            {
                interruptFlag = 1;
                return false;
            }
            if (bitsetHashTable.size() > ENUM_MAX)
            {
                exitEnum = 1;
                break;
            }
            initialPotentialDuplicate &m = prevML[i];

            int s = canonise(m.mask);
            if (s <= ordinal)
            {
                if (stmap.count(s) == 0)
                {
                    initialDuplicateSet ss(currSize + 1, masks.size());
                    ss.insert(m);
                    stmap[s] = ss;
                }
                else
                {
                    stmap[s].insert(m);
                }
            }
            else
                overweight = 1;
        }
        if (exitEnum)
        {
            earlyTerminate = 1;
            prevML.clear();
            delete prevMLptr;
            break;
        }
        if (!overweight)
            tempDag.resize(tempDag.size() + 1);
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (clock() - startTime > runTimeMax)
            {
                interruptFlag = 1;
                return false;
            }
            initialDuplicateSet &ss = it->second;
            if (ss.isValid())
            {
                active = 1;
                if (!overweight)
                    alive |= ss.dagPopulator(currML, maskMap, tempDag);
                int u = ENUM_MAX - currML.size();
                if (bitsetHashTable.size() > u)
                {
                    exitEnum = 1;
                    earlyTerminate = 1;
                    break;
                }
            }
        }
        if (overweight)
            active = 0;
        currSize++;
        prevML.clear();
        delete prevMLptr;
        prevMLptr = currMLptr;
        if (exitEnum)
            break;
    }
    currMLptr->clear();
    delete currMLptr;
    if (exitEnum)
        return false;
    convertDag(tempDag);
    return alive;
}

int dagRecursiveEnumeration(assemblyState &_target, vector<map<int, dagDuplicateSet>> &stmapVector,
                            vector<vector<standardBitset>> &targetMasks)
{
    int ordinal = MAX_INT;
    // Set the maximum index of the fragment that may be chosen
    if (bitsetHashTable.count(_target.masks.front()))
        ordinal = bitsetHashTable[_target.masks.front()].first;
    vector<standardBitset> &masks = _target.masks;
    bool alive = 0;
    size_t currSize = 1;
    int pr = 0;

    // Find all one-bond duplicatable subgraphs
    stmapVector.resize(1);
    for (size_t i = 0; i < masks.size(); i++)
    {
        for (size_t j = 0; j < univEdgeList.size(); j++)
        {
            if (masks[i][j] != 0)
            {
                standardBitset b = 0;
                b.set(j);
                potentialDuplicate m(b, i, j);
                dagGenerate(m, stmapVector[0], masks[i], currSize, ordinal, masks.size());
            }
        }
    }
    bool active = 1, overweight = 0, last = 0;
    while (active)
    {
        vector<standardBitset> targetMask(masks.size(), 0);
        active = 0;
        map<int, dagDuplicateSet> temp;
        stmapVector.push_back(temp);
        map<int, dagDuplicateSet> &stmap = stmapVector[stmapVector.size() - 2];
        // Iterate through all matching subgraphs from the previous matchlist
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            if (interruptFlag)
                return false;
            dagDuplicateSet &ss = it->second;
            if (ss.isValid())
            {
                active |= dagDuplicateGenerator(ss, stmapVector.back(), targetMask, masks, ordinal, overweight, last);
            }
        }
        if (overweight)
            last = 1;
        targetMasks.push_back(targetMask);
        currSize++;
    }
    if (stmapVector.back().size() == 0)
        stmapVector.pop_back();
    for (size_t i = 0; i < masks.size(); i++)
    {
        masks[i] &= targetMasks[0][i];
    }
    return currSize;
}

int postFragmentationCutoff(assemblyState &target, standardBitset &matchMask, standardBitset &maxFragMask)
{
    if (target.masks.size() < 2)
        return 0;
    int maxFragSize = target.masks[0].count();

    /// mainSizeList is the vector of integers corresponding to the sizes of each fragment
    /// matchSizeList is the vector of integers corresponding to the sizes of each matching fragment bitset
    /// maxFragSizeList is the vector of integers corrsponding to the sizes of each maximal size fragment bitset
    vi mainSizeList, matchSizeList, maxFragSizeList;

    mainSizeList.resize(target.masks.size());
    matchSizeList.resize(target.masks.size());
    maxFragSizeList.resize(target.masks.size());

    mainSizeList[0] = maxFragSize;
    matchSizeList[0] = maxFragSize;
    maxFragSizeList[0] = maxFragSize;
    for (size_t i = 1; i < target.masks.size(); i++)
    {
        mainSizeList[i] = target.masks[i].count();
        matchSizeList[i] = (target.masks[i] & matchMask).count();
        maxFragSizeList[i] = (target.masks[i] & maxFragMask).count();
    }
    int matchDB = target.maxDupBonds(mainSizeList, maxFragSize, matchSizeList);
    int maxFragDB = target.maxDupBonds(mainSizeList, maxFragSize, maxFragSizeList) - 1;
    return max(matchDB, maxFragDB);
}

bool dagRecursiveAssembly(assemblyState &input, int &AI)
{
    recursiveCount++;
    if (clock() - startTime > runTimeMax)
        interruptFlag = 1;
    if (interruptFlag)
        return false;
    if (input.AI() < AI)
    {
        AI = input.AI();
        minAIfound = AI;
        minAssemblyPath = input.apPtr;
    }

    vector<map<int, dagDuplicateSet>> stmapVector;
    vector<vector<standardBitset>> targetMasks;
    int maxFragSize = dagRecursiveEnumeration(input, stmapVector, targetMasks);

    if (stmapVector.size() == 0 || interruptFlag)
        return false;

    /// Find the fragment-size-specific AI lower bounds
    vi fragSizeList, fragSizeListMax, sizeList(input.masks.size());
    input.maxDupBonds(fragSizeList, maxFragSize, targetMasks);

    /// Establish the max duplicate cutoff based on the maximum fragment size
    fragSizeListMax.resize(fragSizeList.size());
    fragSizeListMax[0] = fragSizeList[0];
    for (int i = 1; i < fragSizeList.size(); i++)
    {
        if (fragSizeList[i] > fragSizeListMax[i - 1])
            fragSizeListMax[i] = fragSizeList[i];
        else
            fragSizeListMax[i] = fragSizeListMax[i - 1];
    }
    for (size_t i = 0; i < input.masks.size(); i++)
    {
        sizeList[i] = input.masks[i].count();
    }

    /// Begin iterating through the enumerated duplicatable fragments
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        map<int, dagDuplicateSet> &stmap = stmapVector[j];
        vector<standardBitset> stmapMaskList(input.masks.size(), 0);
        standardBitset maskM = 0;
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            dagDuplicateSet &ss = it->second;
            if (!ss.dead)
            {
                vector<validMatchings> matchings;
                standardBitset maskC = 0;

                for (size_t i = 0; i < ss.maskList.size(); i++)
                {
                    maskC |= ss.maskList[i];
                    stmapMaskList[i] |= ss.maskList[i];
                }
                maskM |= maskC;
                int dupBondsMaxFrag = input.maxDupBonds(sizeList, ss.size, stmapMaskList);

                int temp = fragSizeListMax[ss.size - 2] - 1;

                if (ss.size > 2)
                    temp = max(temp, fragSizeListMax[ss.size - 3]);

                /// Initial branch-and-bound before the fragmentation step
                int earlySDP = max(dupBondsMaxFrag, temp);

                int earlyAIBound = totalBonds - input.sumDupBonds - 1 - earlySDP;
                if (earlyAIBound < AI)
                {
                    if (ss.list.size() > 1)
                    {
                        ss.generateMatchings(matchings);
                    }
                    for (int i = matchings.size() - 1; i >= 0; i--)
                    {
                        assemblyState as;
                        fragmentAssemblyState(input, matchings[i], as);

                        int sumDupBonds = input.sumDupBonds + matchings[i].maxFragSize - 1;
                        as.sumDupBonds = sumDupBonds;
                        int temp = postFragmentationCutoff(as, maskC, maskM);
                        if (as.lowBoundAI(matchings[i].maxFragSize, temp) < AI)
                        {
                            apWrapper ap;
                            ap.ap = new assemblyPath;
                            ap.ap->parent = input.apPtr;
                            ap.ap->key = as.assemblyHashCalculator();

                            if (pathAssemblyMap.count(ap) == 0)
                            {
                                assemblyIx++;
                                as.ix = assemblyIx;
                                as.apPtr = ap.ap;
                                ap.ap->sumDupBonds = sumDupBonds;
                                ap.ap->match = bitsetHashTable[matchings[i].first].second;
                                ap.ap->duplicate = bitsetHashTable[matchings[i].second].second;
                                pathAssemblyMap.insert(ap);
                                dagRecursiveAssembly(as, AI);
                            }
                            else
                            {
                                auto it2 = pathAssemblyMap.find(ap);
                                delete ap.ap;
                                if (sumDupBonds > (*it2).ap->sumDupBonds)
                                {
                                    assemblyIx++;
                                    as.ix = assemblyIx;
                                    as.apPtr = (*it2).ap;
                                    (*it2).ap->sumDupBonds = sumDupBonds;
                                    (*it2).ap->match = bitsetHashTable[matchings[i].first].second;
                                    (*it2).ap->duplicate = bitsetHashTable[matchings[i].second].second;
                                    (*it2).ap->parent = input.apPtr;
                                    dagRecursiveAssembly(as, AI);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool initialRecursiveAssembly(assemblyState &input, int &AI)
{
    bool earlyTerminate = 0;
    recursiveCount++;
    if (clock() - startTime > runTimeMax)
        interruptFlag = 1;
    if (interruptFlag)
        return false;
    if (input.AI() < AI)
    {
        AI = input.AI();
        minAIfound = AI;
        minAssemblyPath = input.apPtr;
    }
    vector<map<int, initialDuplicateSet>> stmapVector;
    initialRecursiveEnumeration(input, stmapVector, earlyTerminate);
    if (earlyTerminate)
        return false;
    if (interruptFlag)
        return false;

    if (stmapVector.size() == 0)
        return false;
    /// Begin iterating through the enumerated duplicatable fragments
    for (int j = stmapVector.size() - 1; j >= 0; j--)
    {
        map<int, initialDuplicateSet> &stmap = stmapVector[j];
        for (auto it = stmap.begin(); it != stmap.end(); ++it)
        {
            initialDuplicateSet &ss = it->second;
            vector<validMatchings> matchings;
            if (ss.list.size() > 1)
            {
                ss.generateMatchings(matchings);
            }
            standardBitset maskC = 0;
            for (size_t i = 0; i < ss.maskList.size(); i++)
                maskC |= ss.maskList[i];
            for (int i = matchings.size() - 1; i >= 0; i--)
            {
                assemblyState as;
                fragmentAssemblyState(input, matchings[i], as);
                int sumDupBonds = input.sumDupBonds + matchings[i].maxFragSize - 1;
                as.sumDupBonds = sumDupBonds;
                if (as.lowBoundAI() < AI)
                {
                    apWrapper ap;
                    ap.ap = new assemblyPath;
                    ap.ap->parent = input.apPtr;
                    ap.ap->key = as.assemblyHashCalculator();

                    if (pathAssemblyMap.count(ap) == 0)
                    {
                        assemblyIx++;
                        as.ix = assemblyIx;
                        as.apPtr = ap.ap;
                        ap.ap->sumDupBonds = sumDupBonds;
                        ap.ap->match = bitsetHashTable[matchings[i].first].second;
                        ap.ap->duplicate = bitsetHashTable[matchings[i].second].second;
                        pathAssemblyMap.insert(ap);
                        dagRecursiveAssembly(as, AI);
                    }
                    else
                    {
                        auto it2 = pathAssemblyMap.find(ap);
                        delete ap.ap;
                        if (sumDupBonds > (*it2).ap->sumDupBonds)
                        {
                            assemblyIx++;
                            as.ix = assemblyIx;
                            as.apPtr = (*it2).ap;
                            (*it2).ap->sumDupBonds = sumDupBonds;
                            (*it2).ap->match = bitsetHashTable[matchings[i].first].second;
                            (*it2).ap->duplicate = bitsetHashTable[matchings[i].second].second;
                            (*it2).ap->parent = input.apPtr;
                            dagRecursiveAssembly(as, AI);
                        }
                    }
                }
            }
        }
    }
    return true;
}

void improvedBnB(molGraph &mg)
{
    startTime = clock();
    clearPathMap();
    bitsetHashTable.clear();
    graphHashMap.clear();
    vector<edgeL> removedEdges;
    totalBonds = mg.totalBonds;
    originalEdgeList = mg.writeEdgeList();
    int disjointFragments = mg.disjointFragments();
    originalMolecule = mg;
    targetMolecule = preprocessWriteback(mg, removedEdges);
    univEdgeList = targetMolecule.writeEdgeList();
    allEdges = 0;
    for (size_t i = 0; i < univEdgeList.size(); i++)
        allEdges.set(i);
    assemblyState as;
    as.masks.push_back(allEdges);
    apWrapper ap;
    ap.ap = new assemblyPath;
    ap.ap->parent = nullptr;
    ap.ap->key = as.assemblyHashCalculator();
    ap.ap->sumDupBonds = 0;
    pathAssemblyMap.insert(ap);
    as.apPtr = ap.ap;
    int AI = MAX_INT;
    initialRecursiveAssembly(as, AI);
}
