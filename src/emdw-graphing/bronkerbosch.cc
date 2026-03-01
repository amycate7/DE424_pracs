/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    September 2015
          PURPOSE: The Bron Kerbosch algorithm for finding maximal cliques
          SOURCE: Wiki page at https://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm

*******************************************************************************/

#include "bronkerbosch.hpp"

// standard headers
#include <vector>
#include <algorithm>

// emdw headers
#include "prlite_stlvecs.hpp"
#include "emdw.hpp"
#include "prlite_vector.hpp"
#include "vecset.hpp"


//#include "sortindices.hpp"

using namespace std;
using namespace emdw;

void
BronKerbosch2(
  const map<RVIdType, set<RVIdType> >& neighbours,
  const set<RVIdType>& rSet,
  set<RVIdType>& pSet,
  set<RVIdType>& xSet,
  vector< set<RVIdType> >& maxCliques) {

  if ( pSet.size() or xSet.size() ) { // not done yet

    // here we choose the first available as the pivot element. If
    // instead we chose the one with the highest number of neigbours
    // in P, we get the Tomita algorithm.
    RVIdType u;
    if ( pSet.size() ) {
      u = *( pSet.begin() );
    } // if
    else{
      u = *( xSet.begin() );
    } // else

    set<RVIdType> pNotNeighbourU;
    auto setPtr = neighbours.find(u);
    set_difference(pSet.begin(), pSet.end(),
                   setPtr->second.begin(), setPtr->second.end(),
                   inserter( pNotNeighbourU, pNotNeighbourU.end() ));
    for (auto v : pNotNeighbourU) {

      //cout << "<" << pSet.size()+xSet.size() << ">" << flush;
      set<RVIdType> newRSet(rSet);
      newRSet.insert(v);

      set<RVIdType> newPSet;
      auto setPtr = neighbours.find(v);
      set_intersection(pSet.begin(), pSet.end(),
                       setPtr->second.begin(), setPtr->second.end(),
                       inserter( newPSet, newPSet.end() ));

      set<RVIdType> newXSet;
      set_intersection(xSet.begin(), xSet.end(),
                       setPtr->second.begin(), setPtr->second.end(),
                       inserter( newXSet, newXSet.end() ));

      BronKerbosch2(neighbours, newRSet, newPSet, newXSet, maxCliques);
      //cout << "." << flush;
      pSet.erase(v);
      xSet.insert(v);

    } // for

  } // if

  else { // found a maximal clique
    //cout << "[" << maxCliques.size() << "]" << flush;
    maxCliques.push_back(rSet);
    // for (auto el : rSet) {cout << el << " ";} // for
    // cout << endl;
  } // else

} // BronKerbosch2

void
extendWithVacuousFactors(
  std::vector< rcptr<Factor> >& theFactorPtrs){

  // get all the RV adjacencies
  set<RVIdType> pSet;
  map<RVIdType, set<RVIdType> > adjacent;
  for (auto el : theFactorPtrs) {
    RVIds vars = el->getVars();
    for (size_t i = 0; i < vars.size()-1; i++) {
      pSet.insert(vars[i]);
      for (size_t j = i+1; j < vars.size(); j++) {
        adjacent[vars[i]].insert(vars[j]);
        adjacent[vars[j]].insert(vars[i]);
      } // for j
    } // for i
    pSet.insert(vars[vars.size()-1]);
  } // for el

  // for (auto el : adjacent) {
  //   cout << el.first << ": " << el.second << endl;
  // } // for

  set<RVIdType> rSet;
  set<RVIdType> xSet;
  vector< set<RVIdType> > maxCliques;
  BronKerbosch2(adjacent, rSet, pSet, xSet, maxCliques);
  //JADP16Apr_CLANG cout << maxCliques << endl;
  cout << "Found " << maxCliques.size() << " cliques";

  PRLITE_ASSERT(false,"extendWithVacuousFactors : NIY\n");
} // extendWithVacuousFactors
