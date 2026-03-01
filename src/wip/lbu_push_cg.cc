/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    March 2013
  PURPOSE: PGM cluster graph inference via residual Lauritzen-Spiegelhalter.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_zfstream.hpp"
#include "prlite_vector.hpp"
#include "sortindices.hpp"
//#include "timer.hpp"

// emdw headers
#include "emdw.hpp"
#include "vecset.hpp"
#include "lbu_cg.hpp"
//#include "discretetable.hpp"

// standard headers
#include <new>  // placement new
#include <iostream>

using namespace std;

//long ticksNorm=0, ticksMarg=0, ticksAbsorb=0, ticksCancel=0, ticksDampen=0;
//Timer timer;

//-----------------------------------------------------------------------

/**
 * needs marginalize, inplaceAbsorb and inplaceNormalize (Different
 * from the standard recipe - how you initialize can make a hell of a
 * difference!)
 */
void
setupInitialBeliefsLBU_CG(
  ClusterGraph& cg,
  std::map< Idx2, rcptr<Factor> >& ssBeliefs,
  MessageQueue& msgQueue,
  double thresh) {

  // calculate all messages
  std::map< Idx2, rcptr<Factor> > ssBeliefsD;
  for (auto sepItr = cg.sepvecs.begin(); sepItr != cg.sepvecs.end(); sepItr++) {
    const Idx2& idx2(sepItr->first);
    unsigned srce = idx2.first;
    unsigned dest = idx2.second;

    if (dest > srce) { // dest > srce; do both now

      auto itr = ssBeliefs.find(idx2);

      if (itr  == ssBeliefs.end() ) { // no previous message available
        // srce -> dest
        ssBeliefs[idx2] =
          cg.factorPtrs[srce]->
          marginalize(sepItr->second, true)->normalize();
        // dest -> srce
        ssBeliefsD[idx2] =
          cg.factorPtrs[dest]->
          marginalize(sepItr->second, true)->normalize();



      } // if

      else { // previous message available

      } // else

    } // if

  } // for

  //Absorb messages into their cliques
  for (auto sepItr = cg.sepvecs.begin(); sepItr != cg.sepvecs.end(); sepItr++) {
    const Idx2& idx2(sepItr->first);
    unsigned srce = idx2.first;
    unsigned dest = idx2.second;
    if (dest > srce) { // do both now

      auto itr = ssBeliefsD.find(idx2);
      if (itr  != ssBeliefsD.end() ) {
        //cout << "From: " << srce << " to: " << dest << endl;
        cg.factorPtrs[dest]->inplaceAbsorb(ssBeliefs[idx2]);
        cg.factorPtrs[dest]->inplaceNormalize();

        cg.factorPtrs[srce]->inplaceAbsorb(ssBeliefsD[idx2]);
        cg.factorPtrs[srce]->inplaceNormalize();

        ssBeliefs[idx2]->inplaceAbsorb(ssBeliefsD[idx2]);
        ssBeliefs[idx2]->inplaceNormalize();

        // There is no use in sending a vacuous message
        double delta = ssBeliefs[idx2]->distanceFromVacuous();

        // We do not exclude sending messages to leaves since it might
        // be useful for approximate inference.
        if (delta > thresh) {

          msgQueue.insert( idx2,
                           (delta + (cg.linkedTo[srce].size() == 1) ) );
          msgQueue.insert(Idx2(dest, srce),
                          (delta + (cg.linkedTo[dest].size() == 1) ) );

        } // if

        ssBeliefsD[idx2] = rcptr<Factor>(0);
      } // if
    } // if
  } // for

} // setupInitialBeliefsLBU_CG

//TODO: Remove/update, Old function; reversed direction does not work for leaves
void
scheduleChangedFactorsLBU_CG(
  const ClusterGraph& cg,
  const set<unsigned>& changedFactors,
  MessageQueue& msgQueue) {
  for (auto from : changedFactors) {
    for (auto to : cg.linkedTo[from]) {
      // NOTE: we use the reversed direction to record the *entry*
      // into the cluster - the follow-up message will then be from
      // this cluster to somewhere else.
      Idx2 idx2(to,from);
      msgQueue.insert( idx2, 1.0 + msgQueue.priority(idx2)  );
    } // for
  } // for
} // scheduleChangedFactorsLBU_CG

//-----------------------------------------------------------------------

// Required FactorOperator's: inplaceNormalize inplaceCancel
// inplaceAbsorb inplaceDampen marginalize observeAndReduce.

unsigned
loopyBU_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& ssBeliefs,
  MessageQueue& msgQueue,
  double dampingFactor,
  double deltaThresh,
  unsigned maxNoOfAbsorbs,
  const std::map<emdw::RVIdType, AnyType>& observed) {

  // param controlling the iteration
  if (!maxNoOfAbsorbs) {
    maxNoOfAbsorbs = 33333333;
  } // if
  maxNoOfAbsorbs *= unsigned(sqrt(cg.factorPtrs.size()));
  unsigned noOfAbsorbs = 0;

  // set up initial ssBeliefs if necessary
  if ( !msgQueue.size() ) {
    setupInitialBeliefsLBU_CG(cg, ssBeliefs, msgQueue, deltaThresh);
  } // if

  //cout << "Initial queue: "; msgQueue.show();
  // cout << __FILE__ << __LINE__ << ": " << "Cluster graph:\n" << cg << endl;
  // cout << __FILE__ << __LINE__ << ": " << "ssBeliefs:\n";
  // for (auto ssb : ssBeliefs) {
  //   cout << ssb.first << " " << *ssb.second << endl;
  // } // for

  ofstream errfile;
  errfile.open( "err.txt", static_cast<std::ios::openmode>(std::ios::out) );
  unsigned cnt = 0;
  while (msgQueue.size() and noOfAbsorbs < maxNoOfAbsorbs) {

    //cout << "\nIteration " << cnt << endl;
    //cout << "\nMessage queue:\n"; msgQueue.show(); cout << endl;
    //cout << "\nExisting ssBeliefs:\n";
    //writeMessagesToText(cout, ssBeliefs);

    MessageQueue::MsgPriority toSend = msgQueue.pop();
    // cout << "Pop " << toSend.fromTo
    //      << " from queue (delta = " << toSend.delta
    //      << "), new size is " << msgQueue.size() << endl;

    // Previously fromTo.first sent an individualized message to each
    // each cluster connected to it (including also
    // fromTo.second). Each of these clusters absorbed its
    // message. Now, in turn, we send messages from fromTo.second to
    // be absorbed by all clusters connected to it (excluding
    // fromTo.first).
    unsigned srce = toSend.fromTo.second;

    //if (cnt++ == 10000) {cnt = 0;} // if
    if (cnt++ == 10000) {cnt = 0;} // if
    if (cnt == 0) {

      // write out the current delta status
      cout << toSend.delta << " " << flush;
      errfile << toSend.delta << " " << flush;
      // size_t totSz = 0;
      // for (auto fPtr : cg.factorPtrs) {
      //   rcptr< DiscreteTable<unsigned> > tstDTPtr = std::dynamic_pointer_cast< DiscreteTable<unsigned> > (fPtr);
      //   if (tstDTPtr) {
      //     totSz += tstDTPtr->sparseProbs.size();
      //   } // if
      // } // for
      // cout << "Size: " << totSz << endl;
      // cout << endl
      //      << double(ticksNorm)/1.0 << " "
      //      << double(ticksMarg)/1.0 << " "
      //      << double(ticksAbsorb)/1.0 << " "
      //      << double(ticksCancel)/1.0 << " "
      //      << double(ticksDampen=0)/1.0 << endl;

      ifstream tmpfile;

      // check for request to save the current status of the cluster graph
      tmpfile.open("SAVEPGM");
      if (tmpfile) {
        prlite::gzofstream outfile("current.pgm.gz", static_cast<std::ios::openmode>(std::ios::out) );
        //ofstream outfile("current.pgm", static_cast<std::ios::openmode>(std::ios::out) );
        outfile << cg;
        outfile.close();
        int status = system("rm SAVEPGM");
        if (status) {} // just to get rid of the confounded warning
      } // if

      // check for request to terminate loopy inference
      tmpfile.open("TERMINATE");
      if (tmpfile) {
        noOfAbsorbs = maxNoOfAbsorbs;
        int status = system("rm TERMINATE");
        if (status) {} // just to get rid of the confounded warning
      } // if

    } // if

    // Factor[srce] is ready to distribute ssBeliefs to its destinations
    for (auto dest : cg.linkedTo[srce]) {

      // do not send back over incoming link
      if (dest == toSend.fromTo.first) {
        //cout << "not sending " << srce << " to " << dest << "\n";
        continue;
      } // if

      // cout << "Sending message from " << srce << " to " << dest << endl;
      Idx2 ssIdx( std::min(srce,dest), std::max(srce,dest) );
      rcptr<Factor> oldMsgPtr = ssBeliefs[ssIdx];
      //timer.start();
      rcptr<Factor> newMsgPtr =
        cg.factorPtrs[srce]->
        marginalize(cg.sepvecs[ssIdx], true);
      //ticksMarg += timer.noOfTicks();
      newMsgPtr->inplaceNormalize();
      double delta = newMsgPtr->inplaceDampen(oldMsgPtr, dampingFactor);
      if (delta > deltaThresh) {msgQueue.insert(Idx2(srce,dest), delta);} // if
      cg.factorPtrs[dest]->inplaceAbsorb(newMsgPtr->cancel(oldMsgPtr));
      cg.factorPtrs[dest]->inplaceNormalize();
      ssBeliefs[ssIdx] = newMsgPtr;
      noOfAbsorbs+=2;

    } // for

  } // while
  errfile.close();

  //cout << "Cluster graph: " << cg << endl;
  //cout << "Final ssBeliefs:\n";
  //writeMessagesToText(cout, ssBeliefs);

  if ( msgQueue.size() ) { // not converged
    cout << __FILE__ << __LINE__
         << ": Did not converge in the permitted "
         << maxNoOfAbsorbs << " iterations\n";
    //noOfAbsorbs = 0;
  } // if

  return noOfAbsorbs;

} //loopyBU_CG


unsigned
loopyBU_CG(
    ClusterGraph& cg,
    double dampingFactor,
    double deltaThresh,
    unsigned maxNoOfAbsorbs,
    const std::map<emdw::RVIdType, AnyType>& observed) {

  //provide throw-away message data structures
  map< Idx2, rcptr<Factor> > ssBeliefs;
  MessageQueue msgQueue;
  return loopyBU_CG(cg, ssBeliefs, msgQueue, dampingFactor, deltaThresh, maxNoOfAbsorbs, observed);
}



//-----------------------------------------------------------------------

/*
 * For the moment we are going to assume that the full queryVars will
 * be contained in a single factor somewhere. At a later stage we will
 * have to generalise this to where the queryVars may be spread over
 * several factors.
 */

// Required FactorOperator's: normalize marginalize inplaceAbsorb

uniqptr<Factor>
queryLBU_CG(
  ClusterGraph& cg,
  map< Idx2, rcptr<Factor> >& ssBeliefs,
  const emdw::RVIds& queryVars){

  std::vector<int> l2i, r2i;
  vector<unsigned> varI = cg.factorsWithRV[ queryVars[0] ];
  for (unsigned idx = 1; idx < queryVars.size(); idx++) {
    varI = sortedIntersection(varI, cg.factorsWithRV[ queryVars[idx] ],
                              l2i, r2i);
  } // for

  //Sanity check: If not all in one factor, are they even in any factor?
  if (not varI.size()) {
    std::stringstream errorStream;
    errorStream << "Query currently can only handle cases with all the query variables in a single factor: " << queryVars;
    vector<unsigned> orphanVars;

    for(unsigned var : queryVars) {
      if (cg.factorsWithRV[var].size() == 0) orphanVars.push_back(var);
    }

    if (orphanVars.size()) errorStream << "; also the following are not in any factor: " << orphanVars;

    PRLITE_ASSERT(varI.size(), errorStream.str());
  }

  // varI now contains ALL the factors each of which contains ALL the
  // query vars. Arbitrarily we pick the middle one from the list.
  unsigned factorUsed = varI[0];
  //unsigned factorUsed = varI[varI.size()/2];

  return
    cg.factorPtrs[factorUsed]->
    marginalize(queryVars) ->
    normalize();
} // queryLBU_CG

//-----------------------------------------------------------------------
