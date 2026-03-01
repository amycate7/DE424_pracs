/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    March 2013
  PURPOSE: PGM cluster graph inference via residual Lauritzen-Spiegelhalter.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"
#include "sortindices.hpp"

// emdw headers
#include "emdw.hpp"
#include "vecset.hpp"
#include "lbu2_cg.hpp"

// standard headers
#include <new>  // placement new
#include <iostream>

using namespace std;

//-----------------------------------------------------------------------

/**
 * needs marginalize, inplaceAbsorb and inplaceNormalize (Different
 * from the standard recipe - how you initialize can make a hell of a
 * difference!)
 */
void
setupInitialMessagesLBU2_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& messages,
  MessageQueue& msgQueue,
  bool singleMsg,
  double thresh) {


  // calculate all messages
  for (auto sepItr = cg.sepvecs.begin(); sepItr != cg.sepvecs.end(); sepItr++){
    const Idx2& idx2(sepItr->first);
    unsigned srce = idx2.first;
    unsigned dest = idx2.second;

    auto itr = messages.find(idx2);
    if (itr  == messages.end() ) {
      // add the message
      rcptr<Factor> msgPtr =
        cg.factorPtrs[idx2.first]->marginalize(sepItr->second, true)->normalize();
      messages[idx2] = msgPtr;

      // do not schedule messages to leaf nodes. extra priority for
      // messages from leaf nodes
      // rcptr<Factor> vacPtr = uniqptr<Factor>( msgPtr->vacuousCopy() );
      // vacPtr->inplaceNormalize();
      // double delta =
      //   ( msgPtr->distance(vacPtr) + (cg.linkedTo[srce].size() == 1) ) *
      //   (cg.linkedTo[dest].size() > 1);
      double delta =
        ( msgPtr->distanceFromVacuous() + (cg.linkedTo[srce].size() == 1) ) *
        (cg.linkedTo[dest].size() > 1);
      if (1 or delta > thresh) {
        msgQueue.insert(idx2, delta);
      } // if
    } // if

  } // for

  // absorb all messages into their cliques
  unsigned nFactors = cg.factorPtrs.size();
  for (unsigned dest = 0; dest < nFactors; dest++){
    rcptr<Factor> beliefPtr = cg.factorPtrs[dest];
    for (auto srce : cg.linkedTo[dest]) {
      beliefPtr->inplaceAbsorb(messages[Idx2(srce,dest)]);
      beliefPtr->inplaceNormalize();  // temporary fix for underflow
    } // for
    //beliefPtr->inplaceNormalize();  // this is where it should be
  } // for

  // if single messages are specified, tie the two opposing messages
  if (singleMsg) {
    for (auto sepItr = cg.sepvecs.begin(); sepItr != cg.sepvecs.end(); sepItr++){
      const Idx2& idx2(sepItr->first);
      unsigned srce = idx2.first;
      unsigned dest = idx2.second;
      if (srce > dest) {
        messages[idx2] = messages[Idx2(dest,srce)];
      } // if
    } // for
  } // if

} // setupInitialMessagesLBU2_CG

//TODO: Remove/update, Old function; reversed direction does not work for leaves
void
scheduleChangedFactorsLBU2_CG(
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

// Required FactorOperator's: inplaceNormalize inplaceCancel
// inplaceAbsorb inplaceDampen marginalize observeAndReduce.

// NOTE: Results may (probably will) deviate from that op BP_CG, if an
// absorb reduces a prob to approx zero you can not recover it with a
// cancel.
unsigned
loopyBU2_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& messages,
  MessageQueue& msgQueue,
  double dampingFactor,
  double deltaThresh,
  const std::map<emdw::RVIdType, AnyType>& observed,
  bool singleMsg) {

  // params controlling the iteration
  unsigned noOfAbsorbs = 0;
  unsigned maxNoOfAbsorbs = unsigned(sqrt(cg.factorPtrs.size()))*33333333;

  // set up initial messages if necessary
  setupInitialMessagesLBU2_CG(cg, messages, msgQueue, singleMsg, deltaThresh);
  //cout << "Cluster graph: " << cg << endl;

  ofstream errfile;
  errfile.open( "err.txt", static_cast<std::ios::openmode>(std::ios::out) );
  unsigned cnt = 0;
  while (msgQueue.size() and noOfAbsorbs < maxNoOfAbsorbs) {

    //cout << "\nIteration " << cnt << endl;
    //cout << "\nMessage queue:\n"; msgQueue.show(); cout << endl;
    //cout << "\nExisting messages:\n";
    //writeMessagesToText(cout, messages);

    MessageQueue::MsgPriority toSend = msgQueue.pop();
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

      //cout << "Sending message from " << srce << " to " << dest << endl;
      Idx2 idx2(srce,dest);
      rcptr<Factor> oppMsgPtr = messages[Idx2(dest,srce)];
      rcptr<Factor> oldMsgPtr = messages[idx2];
      rcptr<Factor> newMsgPtr =
        cg.factorPtrs[srce]->
        marginalize(cg.sepvecs[idx2], true);
      newMsgPtr->inplaceCancel(oppMsgPtr);
      newMsgPtr->inplaceNormalize();
      double delta = newMsgPtr->inplaceDampen(oldMsgPtr, dampingFactor);
      if (delta > deltaThresh) {msgQueue.insert(idx2, delta);} // if
      messages[idx2] = newMsgPtr;
      cg.factorPtrs[dest]->inplaceAbsorb(newMsgPtr->cancel(oldMsgPtr));
      cg.factorPtrs[dest]->inplaceNormalize();
      noOfAbsorbs+=2;

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

} //loopyBU2_CG

//-----------------------------------------------------------------------

/*
 * For the moment we are going to assume that the full queryVars will
 * be contained in a single factor somewhere. At a later stage we will
 * have to generalise this to where the queryVars may be spread over
 * several factors.
 */

// Required FactorOperator's: normalize marginalize inplaceAbsorb

uniqptr<Factor>
queryLBU2_CG(
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
  unsigned factorUsed = varI[varI.size()/2];

  return
    cg.factorPtrs[factorUsed]->
    marginalize(queryVars) ->
    normalize();
} // queryLBU2_CG

//-----------------------------------------------------------------------
