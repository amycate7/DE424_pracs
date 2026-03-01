/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

#include "lbu_sweep.hpp"

using namespace std;

/**
 * @brief Calculates a message passing order that never revisits a
 * previously encountered source node.
 *
 * ALGORITHM:
 *
 * INITIALIZATION: A set of initial source nodes (sources) to propagate
 * from. All nodes are potential destinations (availDest). next
 *
 * ITERATION: For all source nodes: Remove the source node from the
 * possible destinations (can not loop back). For each source node,
 * find its neigbours that are still available as destinations. For
 * each of them record the source -> destination message. If that
 * neighbour currently is not also a source, record it in nextLayer
 * (which will be the sources next time around) and .
 */

vector< pair<int, int> > msgOrder(
  const vector< set<unsigned> >& links,
  const set<int>& initNodes) {

  vector< pair<int, int> > messages;
  set<int> nextLayer;

  set<int> sources = initNodes;
  set<int> availDest;

  // Initially all of them
  for (unsigned i = 0; i < links.size(); i++) {
    availDest.insert(i);
  } // for

  // Inside this outer loop we first fully process the current set of
  // sources and only after that, repack their neigbours as the next
  // layer of sources.
  while ( !availDest.empty() ) {

    for (auto s : sources) { // from each source node
      //cout << s << " ";
      availDest.erase(s); // can not loop back to the source

      for (auto neigh = links[s].begin(); neigh != links[s].end(); neigh++){ // all neighbouring nodes

        if (availDest.find(*neigh) != availDest.end()) { // if the neighbour is still available as destination
          //is present in set
          messages.push_back(pair<int, int>(s, *neigh));  // from to

          if (sources.find(*neigh) == sources.end()) { // only if the neighbour is not also a source
            //is not present in set
            nextLayer.insert(*neigh);
          } //if

        } // if
      } // for
    } // for
    //cout << endl;

    sources = nextLayer; // new set of sources
    nextLayer.clear(); // and rebuild the next layer from scratch
  } // while

  return messages;
} // msgOrder

/**
 * @brief Does full cluster sweep Lauritzen Spiegelhalter message
 * passing on an RIP clustergraph.
 *
 * IMPLEMENTATION: I considered using a vector for ssBeliefs (faster)
 * instead of a map. However, since the map is independent of the
 * message order, it allows us to use the query code from lbu_cg
 * directly, and also to even change the message order on subsequent
 * runs.
 */

unsigned lbu_sweep(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& ssBeliefs,
  vector< pair<int, int> >& forwardSweep,
  double dampingFactor,
  double distThresh,
  unsigned maxNoOfSweeps) {

  vector< rcptr<Factor> >& factorPtrs(cg.factorPtrs);
  map<Idx2, emdw::RVIds>& sepvecs(cg.sepvecs);

  if (forwardSweep.size() == 1) {
    forwardSweep = msgOrder( cg.linkedTo, set<int>({forwardSweep[0].first}) );
  } // if

  else if ( !forwardSweep.size() ) {
    set<int> initNode;

    // lets first search for leaf nodes
    for (unsigned srce = 0; srce < cg.linkedTo.size(); srce++) {
      if (cg.linkedTo[srce].size() == 1) {
        initNode.insert(srce);
      } // if
    } // for

    if ( !initNode.size() ) { // damn, no leaf nodes - find most informative node
      unsigned bestIdx = 0; double bestDist = cg.factorPtrs[0]->distanceFromVacuous();
      for (unsigned idx = 1; idx < cg.factorPtrs.size(); idx++) {
        double dist = cg.factorPtrs[idx]->distanceFromVacuous();
        if (dist > bestDist){
          bestIdx = idx;
          bestDist = dist;
        } // if
      } // for
      initNode.insert(bestIdx);
      //cout << "best one is " << bestIdx << " with value " << bestDist << endl;
    } // if

    forwardSweep = msgOrder(cg.linkedTo, initNode);
  } // else if

  // 'changed' records whether the particular cluster changed during
  // the most recent message passing? The .first is for changes during
  // forward sweep, .second for changes during backward sweep.  Since
  // we'll start off with a forward sweep, we initialize as if the
  // preceding backward sweep activated all clusters.
  vector< pair<bool,bool> > changed( factorPtrs.size() );
  for (auto& el : changed) {el.second = true;}

  double dist;
  double distMax;
  bool notDone = true;
  unsigned noOfAbsorbs = 0;
  unsigned noOfSweeps = 0;

  while (notDone and noOfSweeps < maxNoOfSweeps) {

    notDone = false;
    distMax = 0.0;

    // *********************
    // *** Forward sweep ***
    // *********************

    for (auto& el : changed) {el.first = false;} // clear all the forward change markers
    for (auto fSweepItr = forwardSweep.begin(); fSweepItr != forwardSweep.end(); fSweepItr++) {

      auto srce = fSweepItr->first;
      if (changed[srce].first or changed[srce].second) { // only marginalise if source has changed
        rcptr<Factor> tmp = factorPtrs[srce]->marginalize(sepvecs[*fSweepItr]);

        auto ssItr = fSweepItr->first < fSweepItr->second ?
          ssBeliefs.find({fSweepItr->first,fSweepItr->second}) :
          ssBeliefs.find({fSweepItr->second,fSweepItr->first});

        if ( ssItr != ssBeliefs.end() ) { // been here before
          dist = tmp->inplaceDampen(ssItr->second, dampingFactor);
          if (dist > distThresh) {
            auto dest = fSweepItr->second;
            factorPtrs[dest]->inplaceAbsorb( tmp->cancel(ssItr->second) );
            factorPtrs[dest]->inplaceNormalize();
            changed[dest].first = true;
            notDone = true;
            noOfAbsorbs+=2;
          } // if
          ssItr->second = tmp; // and update the message
        } // if

        else{ // first time around
          dist = tmp->distanceFromVacuous();
          if (dist > distThresh) {
            auto dest = fSweepItr->second;
            factorPtrs[dest]->inplaceAbsorb(tmp);
            factorPtrs[dest]->inplaceNormalize();
            changed[dest].first = true;
            notDone = true;
            //noOfAbsorbs+=1;
          } // if
          ssBeliefs.insert({ {min(fSweepItr->first,fSweepItr->second), max(fSweepItr->first,fSweepItr->second)}, tmp});
        } // else

        if (dist > distMax) {distMax = dist;} // if

      } // if

    } // for forward pass

      // **********************
      // *** Backward sweep ***
      // **********************

    for (auto& el : changed) {el.second = false;} // clear all the backward change markers
    for (auto bSweepItr = forwardSweep.rbegin(); bSweepItr != forwardSweep.rend(); bSweepItr++) {

      auto srce = bSweepItr->second;
      if (changed[srce].first or changed[srce].second) { // only marginalise if source has changed
        rcptr<Factor> tmp = factorPtrs[srce]->marginalize(sepvecs[*bSweepItr]);

        auto ssItr = bSweepItr->first < bSweepItr->second ?
          ssBeliefs.find({bSweepItr->first,bSweepItr->second}) :
          ssBeliefs.find({bSweepItr->second,bSweepItr->first});

        dist = tmp->inplaceDampen(ssItr->second, dampingFactor);
        if (dist > distThresh) {
          auto dest = bSweepItr->first;
          factorPtrs[dest]->inplaceAbsorb( tmp->cancel(ssItr->second) );
          factorPtrs[dest]->inplaceNormalize();
          changed[dest].second = true;
          notDone = true;
          noOfAbsorbs+=2;
        } // if

        if (dist > distMax) {distMax = dist;} // if
        ssItr->second = tmp; // and update the message

      } // if

    } // for backward pass

    cout << distMax << " ";
    noOfSweeps++;
  } // while notDone

  return noOfAbsorbs;
} // lbu_sweep
