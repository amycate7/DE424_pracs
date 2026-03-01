/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    March 2013
  PURPOSE: PGM cluster graph inference via residual Shafer-Shenoy.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"
//#include "prlite_stlvecs.hpp"  // FIXME: for some reason this must be after prlite_vector.hpp
#include "sortindices.hpp"

// emdw headers
#include "emdw.hpp"
#include "vecset.hpp"
#include "lbp_cg.hpp"

// standard headers
#include <new>  // placement new
#include <iostream>
#include <queue>

using namespace std;

//-----------------------------------------------------------------------

/**
 * needs marginalize, inplaceAbsorb and inplaceNormalize (Different
 * from the standard recipe - how you initialize can make a hell of a
 * difference!)
 */
void
setupInitialMessagesLBP_CG(
  const ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& messages,
  MessageQueue& msgQueue,
  double thresh) {

  // How far a transition is from closest leaf node.
  map< pair<unsigned,unsigned>,unsigned> hopsFromLeaf;
  queue< tuple<unsigned, unsigned, unsigned> > newOnes; // srce, dest, hops
  // initialize
  for (unsigned srce = 0; srce < cg.linkedTo.size(); srce++) {
    if (cg.linkedTo[srce].size() == 1) {
      newOnes.push({srce, *cg.linkedTo[srce].begin(), 1});
    } // if
  } // for

  // and recursively find all connected
  while ( !newOnes.empty() ) {
    auto [oldSrce, newSrce, numHops] = newOnes.front();
    newOnes.pop(); // pops oldest one from front
    if ( hopsFromLeaf.find({oldSrce, newSrce}) == hopsFromLeaf.end() ) { // been here?
      hopsFromLeaf.insert({ {oldSrce, newSrce}, numHops});
      for (auto dest : cg.linkedTo[newSrce]) { // new destinations to visit
        if (dest != oldSrce) {
          newOnes.push({newSrce,dest,numHops+1}); // newest goes to the back
        } // if
      } // for
    } // if
  } // while

  // calculate all messages
  for (auto sepItr = cg.sepvecs.begin(); sepItr != cg.sepvecs.end(); sepItr++){

    const Idx2& idx2(sepItr->first);
    unsigned srce = idx2.first;
    unsigned dest = idx2.second;
    rcptr<Factor> msgPtr;

    if ( messages.find(idx2) == messages.end() ) {

      // add the message
      msgPtr =
        cg.factorPtrs[idx2.first]->marginalize(sepItr->second, true)->normalize();
      messages[idx2] = msgPtr;

      // do not schedule messages to leaf nodes. extra priority for
      // messages from leaf nodes
      // rcptr<Factor> vacPtr = uniqptr<Factor>( msgPtr->vacuousCopy() );
      // vacPtr->inplaceNormalize();

      // double delta =
      //   ( msgPtr->distance(vacPtr) + (cg.linkedTo[srce].size() == 1) ) *
      //   (cg.linkedTo[dest].size() > 1);

      auto itr = hopsFromLeaf.find({srce,dest});
      double cudos = (itr != hopsFromLeaf.end() ? 1.0/itr->second : 0.0);
      //double cudos = (cg.linkedTo[dest].size() > 1);

      double delta = ( msgPtr->distanceFromVacuous() + cudos ) *
        (cg.linkedTo[dest].size() > 1);

      if (1 or delta > thresh) {
        msgQueue.insert(idx2, delta);
      } // if

    } // if

  } // for

} // setupInitialMessagesLBP_CG

//TODO: Remove/update, Old function; reversed direction does not work for leaves
void
scheduleChangedFactorsLBP_CG(
  const ClusterGraph& cg,
  const set<unsigned>& changedFactors,
  MessageQueue& msgQueue) {
  for (auto from : changedFactors) {
    for (auto to : cg.linkedTo[from]) {
      // NOTE: we use the reversed direction to record the *entry*
      // into the cluster - the follow-up will message will then be
      // from this cluster to somewhere else.
      Idx2 idx2(to,from);
      msgQueue.insert( idx2, 1.0 + msgQueue.priority(idx2)  );
    } // for
  } // for
} // scheduleChangedFactorsLBP_CG

//-----------------------------------------------------------------------

// Required FactorOperator's: inplaceNormalize inplaceAbsorb
// observeAndReduce marginalize absorb inplaceDampen
unsigned
loopyBP_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& messages,
  MessageQueue& msgQueue,
  double dampingFactor,
  double deltaThresh,
  unsigned maxNoOfAbsorbs,
  const std::map<emdw::RVIdType, AnyType>& observed) {

  //long ticks1 = 0; //,ticks2 = 0;
  //Timer timer1; //, timer2;

  // param controlling the iteration
  if (!maxNoOfAbsorbs) {
    maxNoOfAbsorbs = 99999999;
  } // if
  maxNoOfAbsorbs *= unsigned(sqrt(cg.factorPtrs.size()));
  unsigned noOfAbsorbs = 0;

  // set up initial messages if necessary
  if ( !msgQueue.size() ) {
    setupInitialMessagesLBP_CG(cg, messages, msgQueue, deltaThresh);
  } // if
  // cout << "Cluster graph: " << cg << endl;

  //cout << "Show initial queue\n";
  //msgQueue.show();

  ofstream errfile;
  errfile.open( "err.txt", static_cast<std::ios::openmode>(std::ios::out) );
  unsigned cnt = 0;
  while (msgQueue.size() and noOfAbsorbs < maxNoOfAbsorbs) {
    //timer1.start();

    //cout << "\nIteration " << cnt << endl;
    //cout << "\nMessage queue:\n"; msgQueue.show(); cout << endl;
    //cout << "\nExisting messages:\n";
    //writeMessagesToText(cout, messages);

    MessageQueue::MsgPriority toSend = msgQueue.pop();
    // cout << "Pop " << toSend.fromTo
    //      << " from queue (delta = " << toSend.delta
    //      << "), new size is " << msgQueue.size() << endl;
    unsigned srce = toSend.fromTo.second;

    if (cnt++ == 10000) {cnt = 0;} // if
    if (cnt == 0) {
      cout << toSend.delta << " " << flush;
      errfile << toSend.delta << " " << flush;
    } // if

    // Factor[srce] is ready to distribute messages to its destinations
    for (auto dest : cg.linkedTo[srce]) {

      // do not send back over incoming link
      if (dest == toSend.fromTo.first) {
        //cout << "not sending " << srce << " to " << dest << "\n";
        continue;
      } // if

      //cout << "Sending message from " << srce << " to " << dest << " (" << toSend.delta << ")" << endl;
      rcptr<Factor> oldMsgPtr = messages[Idx2(srce,dest)];
      rcptr<Factor> newMsgPtr;
      // this is rather expensive with highly linked factors
      //cout << "combining previous messages from: ";
      for (auto prevSrce : cg.linkedTo[srce]) {
        if (prevSrce != dest) {

          //cout << "Using previous message from " << prevSrce << " to " << srce << endl;
          // sets up the type of the new message and do first
          // absorb. NOTE: Might speedup by setting the type upfront
          if (!newMsgPtr) {
            newMsgPtr = cg.factorPtrs[srce]->absorb(messages[Idx2(prevSrce,srce)]);
            newMsgPtr->inplaceNormalize();  //added here to guard against underflow
            noOfAbsorbs++;
          } // if

          // while this one should be faster. NOTE: We are assuming
          // that the absorb will not change the type
          else {
            newMsgPtr->inplaceAbsorb(messages[Idx2(prevSrce,srce)]);
            newMsgPtr->inplaceNormalize();  //added here to guard against underflow
            noOfAbsorbs++;
          } // else

        } // if
      } // for

      if (newMsgPtr) {
        // true because sepvecs are presorted
        newMsgPtr =
          newMsgPtr->marginalize(cg.sepvecs[Idx2(srce,dest)], true);
        newMsgPtr->inplaceNormalize();
      } // if

      // if newMsgPtr is null, it is from a leaf. Normally we will ignore
      // it as static, but oracle input might have changed it.
      else {
        newMsgPtr =
          cg.factorPtrs[srce]->
          marginalize(cg.sepvecs[Idx2(srce,dest)], true);
        // no need to normalize if original factor was ok, but we cannot always make that assumption
        newMsgPtr->inplaceNormalize();
      } // else

      double delta = newMsgPtr->inplaceDampen(oldMsgPtr, dampingFactor);
      //cout << "delta: " << delta << endl;
      messages[Idx2(srce,dest)] = newMsgPtr;
      if (delta > deltaThresh) {msgQueue.insert(Idx2(srce,dest), delta);} // if
    } // for

  } // while
  errfile.close();

  //cout << "Cluster graph: " << cg << endl;
  //cout << "Final messages:\n";
  //writeMessagesToText(cout, messages);

  if ( msgQueue.size() ) { // not converged
    cout << __FILE__ << __LINE__ << ": Did not converge in the permitted " << maxNoOfAbsorbs << " iterations\n";
  } // if

  //cout << __FILE__ << __LINE__ << endl;
  //cout << ticks1 << endl;

  return noOfAbsorbs;

} //loopyBP_CG

//-----------------------------------------------------------------------

/*
 * For the moment we are going to assume that the full queryVars will
 * be contained in a single factor somewhere. At a later stage we will
 * have to generalise this to where the queryVars may be spread over
 * several factors.
 */

// Required FactorOperator's: normalize marginalize inplaceAbsorb

uniqptr<Factor>
queryLBP_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& messages,
  const emdw::RVIds& queryVars){

  std::vector<int> l2i, r2i;
  vector<unsigned> varI = cg.factorsWithRV[ queryVars[0] ];
  for (unsigned idx = 1; idx < queryVars.size(); idx++) {
    varI = sortedIntersection(varI, cg.factorsWithRV[ queryVars[idx] ],
                              l2i, r2i);
  } // for
  PRLITE_ASSERT(varI.size(),
         "Query currently can only handle cases with all the query variables in a single factor: " << queryVars);

  // varI now contains ALL the factors each of which contains ALL the
  // query vars. Arbitrarily we pick the middle one from the list.
  //unsigned factorUsed = varI[varI.size()/2];
  unsigned factorUsed = varI[0];
  uniqptr<Factor> resultPtr = uniqptr<Factor>( cg.factorPtrs[factorUsed]->copy() );

  for (auto srce : cg.linkedTo[factorUsed]) {
    resultPtr-> inplaceAbsorb(messages[Idx2(srce,factorUsed)]);
    resultPtr->inplaceNormalize();      // here to protect against underflow
    } // for

  return
    resultPtr->
    marginalize(queryVars) ->
    normalize();                        // this might be unnecessary with the above inplaceNormalize

} // queryLBP_CG

//-----------------------------------------------------------------------

double
queryLBP_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& messages,
  const emdw::RVIds& queryVars,
  const emdw::RVVals& queryVals,
  unsigned& factorUsed,
  const std::map<emdw::RVIdType, AnyType>& observed){

  emdw::RVIds unobservedVars;
  unobservedVars.reserve( queryVars.size() );
  emdw::RVVals unobservedVals;
  unobservedVals.reserve( queryVals.size() );

  for (unsigned idx = 0; idx < queryVars.size(); idx++) {

    // first check whether it is an observed RV.  we can probably do
    // this faster by first assuming that it is unobserved, and if
    // this fails, then check if it is observed.

    auto itr = observed.find( queryVars[idx] );

    if ( itr != observed.end() ) { // is observed
      factorUsed = cg.noOfFactors();
      if ( itr->second != queryVals[idx] ) {  // !!! we'll have to do this more carefully with continuous vars
        return 0.0;
      } // if
    } // if

    else { // not observed
      unobservedVars.push_back(queryVars[idx]);
      unobservedVals.push_back(queryVals[idx]);
    } // else

  } // for

  if ( !unobservedVars.size() ) {
    factorUsed = cg.noOfFactors();
    return 1.0;
  } // everything was observed and matched

    std::vector<int> l2i, r2i;
    vector<unsigned> varI = cg.factorsWithRV[ unobservedVars[0] ];
    for (unsigned idx = 1; idx < unobservedVars.size(); idx++) {
      varI = sortedIntersection(varI, cg.factorsWithRV[ unobservedVars[idx] ],
                                l2i, r2i);
    } // for
  PRLITE_ASSERT(varI.size(),
         "Query currently can only handle cases with all the query variables in a single factor: " << queryVars);

  //factorUsed = varI[varI.size()/2];
  factorUsed = varI[0];
  uniqptr<Factor> resultPtr = uniqptr<Factor>( cg.factorPtrs[factorUsed]->copy() );

  for (auto srce : cg.linkedTo[factorUsed]) {
    resultPtr->inplaceAbsorb(messages[Idx2(srce,factorUsed)]);
    resultPtr->inplaceNormalize();      // here to protect against underflow
  } // for

  return
    resultPtr->
    marginalize(unobservedVars) ->
    normalize()->                       // this might be unnecessary with the above inplaceNormalize
    potentialAt(unobservedVars, unobservedVals);

} // queryLBP_CG

//-----------------------------------------------------------------------
