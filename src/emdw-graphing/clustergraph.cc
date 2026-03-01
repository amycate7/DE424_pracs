/*******************************************************************************

          AUTHOR:    JA du Preez, SF Streicher
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:      March 2013 - 2018
          PURPOSE:   ClusterGraph PGM cluster graph

*******************************************************************************/

// patrec headers
#include "prlite_vector.hpp"
#include "prlite_stlvecs.hpp"  // FIXME: for some reason this must be after prlite_vector.hpp
#include "sortindices.hpp"
#include "prlite_zfstream.hpp"

// emdw headers
#include "emdw.hpp"
#include "vecset.hpp"
#include "oddsandsods.hpp"
#include "clustergraph.hpp"
#include "mojo.hpp"
#include "discretetable.hpp"
#include "treefactor.hpp"


// standard headers
#include <new>        // placement new
#include <algorithm>  // std::set_intersection, std::sort
#include <iostream>
#include <list>
//#include <random>
#include <limits>
#include <iomanip>

using namespace std;
using emdw::RVIdType;
using emdw::RVIds;


/**************************The ClusterGraph class members***********************/

//============================ Traits and Typedefs ============================

ClusterGraph::_LTRIP_  ClusterGraph::LTRIP;
ClusterGraph::_BETHE_  ClusterGraph::BETHE;
ClusterGraph::_JTREE_  ClusterGraph::JTREE;
ClusterGraph::_LPRODS_  ClusterGraph::LPRODS;
ClusterGraph::_LBAGS_  ClusterGraph::LBAGS;
ClusterGraph::_RLJTREES_  ClusterGraph::RLJTREES;

//======================== Constructors and Destructor ========================

//

// LTRIP ctor
ClusterGraph::ClusterGraph(
  _LTRIP_,
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy,
  std::map<size_t,double> theMassOfCluster)
    : factorPtrs( theFactorPtrs.size() ) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  // use the observed variables to simplify the factors
  observeAndReduce(factorPtrs, observed, factorsWithRV);

  std::vector<std::map<size_t, emdw::RVIds> > theIntersections;
  std::map<size_t, std::set<size_t> > theReachableFrom;
  absorbSubsetFactors(factorPtrs, theIntersections, theReachableFrom, factorsWithRV);

  if ( !theMassOfCluster.size() ) {
      for (size_t m = 0; m < factorPtrs.size(); m++) {
      theMassOfCluster[m] = 1E-4 + factorPtrs[m]->distanceFromVacuous();
    } // for
  } // if

  std::map<size_t, Mojo> mojoMap(
    calcGravityMojos(factorPtrs, theIntersections, theReachableFrom, theMassOfCluster, false) );

  // for (unsigned k = 0; k < factorPtrs.size(); k++) {
  //   cout << "Factor " << k << ": " << *factorPtrs[k] << endl;
  // } // for
  //     cout << "mojoMap: " << endl;
  //     for (auto el: mojoMap) {
  //        cout << el.first << ": "; el.second.show(1);
  //     } // for

  map< Idx2, set<emdw::RVIdType> > sepsets;
  linkedTo.resize( factorPtrs.size() );

  for (auto fwRV : factorsWithRV){
    //cout << "Setting up subgraph for varId: " << fwRV.first << endl;

    if (fwRV.second.size() > 1) { // must be connected to something

      // *********************************************************************

      // extract the subset of mojo's applicable to this particular RV
      std::map<size_t, Mojo> rvMojos;
      // std::cout << "Factor set for RV= " << fwRV.first << ": " << fwRV.second << std::endl;
      for (auto to : fwRV.second) {
        auto toMojoItr = mojoMap.find(to);
        Mojo toMojo(to);
        for (auto from : fwRV.second) {
          if (from == to) {continue;}
          // cout << "Setting up mojo " << to << "<-" << from << endl;
          auto fromItr = toMojoItr->second.find(from);
          if ( fromItr == toMojoItr->second.end() ) {
            cout << __FILE__ << __LINE__ << " mighty odd, this!\n";
          } // if
          toMojo.insert(fromItr->first, fromItr->second);
        } // for
        rvMojos.insert( make_pair(to, toMojo) );
      } // for

      // cout << "rvMojos: " << endl;
      // for (auto el: rvMojos) {
      //    cout << el.first << ": "; el.second.show(1);
      // } // for

      // build an mst for these factors

      vector<Idx2> varLinks;
      primsAlgo(rvMojos, varLinks);
      // cout << __FILE__ << __LINE__ << "\n";
      // cout << "varLinks (" << varLinks.size() << "):\n";
      // for (auto el : varLinks) {
      //   cout << el.first << ":" << el.second << endl;
      // } // for

      // *********************************************************************

      for (auto aLink : varLinks){

        unsigned from = aLink.first; //fwRV.second[aLink.first];
        unsigned to = aLink.second; //fwRV.second[aLink.second];
        //cout << "Linking factor " << from << " to factor " << to << endl;
        //cout << "Linking factor " << factorPtrs[from] << " to factor " << factorPtrs[to] << endl;
        // record that there is a link - both directions
        linkedTo[from].insert(to);
        linkedTo[to].insert(from);

        // Add varId to the sepset.
        sepsets[Idx2(from, to)].insert(fwRV.first);
        sepsets[Idx2(to, from)].insert(fwRV.first);

      } // for

    } // if

  } // for


  // cout << __FILE__ << __LINE__ << endl;
  // cout << "factors after absorbing leaves:\n" << factorPtrs.size() << endl;
  // for (auto idx = 0; idx < factorPtrs.size(); idx++){
  //   cout << "factor_" << idx << " at: " << factorPtrs[idx] << endl;
  //   if (factorPtrs[idx]){
  //       cout << *factorPtrs[idx] << endl;
  //   } // if
  // } // for

  // copy the sepset into the sepvec. Note: they therefore are
  // presorted automatically.

  sepvecs.clear();
  for (auto itr = sepsets.begin(); itr != sepsets.end(); itr++){
    sepvecs[itr->first] =
      emdw::RVIds( itr->second.begin(), itr->second.end() );
  } // for

} // LTRIP ctor

// BETHE ctor
ClusterGraph::ClusterGraph(
  _BETHE_,
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy)
    : factorPtrs( theFactorPtrs.size() ) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  // use the observed variables to simplify the factors
  observeAndReduce(factorPtrs, observed, factorsWithRV);

  std::vector<std::map<size_t, emdw::RVIds> > theIntersections;
  std::map<size_t, std::set<size_t> > theReachableFrom;
  absorbSubsetFactors(factorPtrs, theIntersections, theReachableFrom, factorsWithRV);

  linkedTo.resize( factorPtrs.size() );

  for (auto& kv : factorsWithRV){
    auto N = factorPtrs.size();
    unsigned variableFactorIndex = N;
    for (auto& i : kv.second){
      if (factorPtrs[i]->noOfVars() <= 1){
        variableFactorIndex = i;
        break;
      }
    }

    if (variableFactorIndex >= N){
      // Create unit variable factor for current RV
      rcptr<Factor> variableFactor = rcptr<Factor>(
        factorPtrs[kv.second[0]]->vacuousCopy({kv.first}, true) );
      variableFactor->keepSeparate = true;
      // Add variable factor to list
      factorPtrs.push_back(variableFactor);
      linkedTo.resize( factorPtrs.size() );
      kv.second.push_back(variableFactorIndex);
      variableFactorIndex = N;
    }

    // Make links to variable factor
    for (auto& i : kv.second){
      if (i == variableFactorIndex) continue;
      linkedTo[i].insert(variableFactorIndex);
      linkedTo[variableFactorIndex].insert(i);
      sepvecs[Idx2(i, variableFactorIndex)] = {kv.first};
      sepvecs[Idx2(variableFactorIndex, i)] = {kv.first};
    }
  }

} // BETHE ctor

// Junction Tree ctor
ClusterGraph::ClusterGraph(
  _JTREE_,
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy)
    : factorPtrs( theFactorPtrs.size() ) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  if (observed.size()) { //observe if given, Junction tree takes care of subsets etc.
    observeAndReduce(factorPtrs, observed, factorsWithRV);
  }

  buildJunctionTree(factorPtrs, factorsWithRV);
} // Junction Tree ctor

//#include "lbu_cg.hpp"

// Loopy factor products ctor
ClusterGraph::ClusterGraph(
  _LPRODS_,
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy,
  size_t maxSz)
    : factorPtrs( theFactorPtrs.size() ) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  // use the observed variables to simplify the factors
  observeAndReduce(factorPtrs, observed, factorsWithRV);
  cout << "After observing " << observed.size()
       << " variables we now have " << factorPtrs.size() << " factors\n";

  std::vector<std::map<size_t, emdw::RVIds> > theIntersections;
  std::map<size_t, std::set<size_t> > theReachableFrom;
  absorbSubsetFactors(factorPtrs, theIntersections, theReachableFrom, factorsWithRV);
  cout << "After absorbing subsets we now have " << factorPtrs.size() << " factors\n";

  // Initial cluster factors, masses and variables are just those of
  // the individual factors. We are going to record these in maps
  // since the clustering process will keep on reducing them, and we
  // don't want to renumber.
  map< size_t, set<size_t> > factorsInCluster;
  map<size_t,double> massOfCluster;
  map< size_t, set<RVIdType> > varsInCluster;
  unsigned M = factorPtrs.size();
  for (size_t m = 0; m < M; m++) {
    factorsInCluster[m].insert({m});
    massOfCluster[m] = 1E-4 + factorPtrs[m]->distanceFromVacuous();
    const RVIds& theVars( factorPtrs[m]->getVars() );
    varsInCluster[m].insert( theVars.begin(), theVars.end() );
  } // for

  // record the initial gravity accelarations between the various clusters
  std::map<size_t, Mojo> mojoMap(
    calcGravityMojos(factorPtrs, theIntersections, theReachableFrom, massOfCluster, false) );

  // and cluster the factors accordingly
  clusterGravity(mojoMap, theReachableFrom, massOfCluster, factorsInCluster, varsInCluster, maxSz);
  std::cout << "Clustered factors into " << massOfCluster.size() << " groups\n";

  // ***************************************************************************************
  // ****************** Now build a factor for each cluster ********************************
  // ***************************************************************************************


  vector< rcptr<Factor> > clPtrs;
  for (auto cl : factorsInCluster) {
    vector<size_t> indices( cl.second.begin(), cl.second.end() );
    clPtrs.push_back( absorb( extract< rcptr<Factor> >(factorPtrs, indices) )->normalize() );
  } // for
  cout << endl;

  rcptr<ClusterGraph> cgPtr;
  cgPtr = uniqptr<ClusterGraph>(
    new ClusterGraph(ClusterGraph::LTRIP, clPtrs, observed, false) );
  //new ClusterGraph(ClusterGraph::LTRIP, clPtrs, observed, true, massOfCluster) );

// we should be able to find all the parameters for an inplace construction in cgPtr
  this->~ClusterGraph(); // Destroy existing

  new(this) ClusterGraph(std::move(cgPtr->factorPtrs), std::move(cgPtr->linkedTo), std::move(cgPtr->sepvecs), std::move(cgPtr->factorsWithRV), false);
} // Loopy Prods ctor

// Loopy factor products ctor
ClusterGraph::ClusterGraph(
  _LJTREES_,
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy,
  size_t maxSz)
    : factorPtrs( theFactorPtrs.size() ) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  // use the observed variables to simplify the factors
  observeAndReduce(factorPtrs, observed, factorsWithRV);
  cout << "After observing " << observed.size()
       << " variables we now have " << factorPtrs.size() << " factors\n";

  std::vector<std::map<size_t, emdw::RVIds> > theIntersections;
  std::map<size_t, std::set<size_t> > theReachableFrom;
  absorbSubsetFactors(factorPtrs, theIntersections, theReachableFrom, factorsWithRV);
  cout << "After absorbing subsets we now have " << factorPtrs.size() << " factors\n";

  // Initial cluster factors, masses and variables are just those of
  // the individual factors. We are going to record these in maps
  // since the clustering process will keep on reducing them, and we
  // don't want to renumber.
  map< size_t, set<size_t> > factorsInCluster;
  map<size_t,double> massOfCluster;
  map< size_t, set<RVIdType> > varsInCluster;
  unsigned M = factorPtrs.size();
  for (size_t m = 0; m < M; m++) {
    factorsInCluster[m].insert({m});
    massOfCluster[m] = 1E-4 + factorPtrs[m]->distanceFromVacuous();
    const RVIds& theVars( factorPtrs[m]->getVars() );
    varsInCluster[m].insert( theVars.begin(), theVars.end() );
  } // for

  // record the initial gravity accelarations between the various clusters
  std::map<size_t, Mojo> mojoMap(
    calcGravityMojos(factorPtrs, theIntersections, theReachableFrom, massOfCluster, false) );

  // and cluster the factors accordingly
  clusterGravity(mojoMap, theReachableFrom, massOfCluster, factorsInCluster, varsInCluster, maxSz);
  std::cout << "Clustered factors into " << massOfCluster.size() << " groups\n";

  // ***************************************************************************************
  // ****************** Now build a factor for each cluster ********************************
  // ***************************************************************************************


  vector< rcptr<Factor> > clPtrs;
  for (auto cl : factorsInCluster) {
    vector<size_t> indices( cl.second.begin(), cl.second.end() );
    clPtrs.push_back(
      rcptr<Factor>( new TreeFactor(extract< rcptr<Factor> >(factorPtrs, indices), false) ) );
  } // for

  rcptr<ClusterGraph> cgPtr;
  cgPtr = uniqptr<ClusterGraph>(
    new ClusterGraph(ClusterGraph::LTRIP, clPtrs, observed, false) );
  //new ClusterGraph(ClusterGraph::LTRIP, clPtrs, observed, true, massOfCluster) );

// we should be able to find all the parameters for an inplace construction in cgPtr
  this->~ClusterGraph(); // Destroy existing

  new(this) ClusterGraph(std::move(cgPtr->factorPtrs), std::move(cgPtr->linkedTo), std::move(cgPtr->sepvecs), std::move(cgPtr->factorsWithRV), false);
} // Loopy Trees ctor

// class specific ctor
ClusterGraph::ClusterGraph(
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::vector< std::set<unsigned> >& theLinkedTo,
  const std::map<Idx2, emdw::RVIds>& theSepvecs,
  const std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV,
  const bool deepCopy)
    : factorPtrs( theFactorPtrs.size() ),
      linkedTo(theLinkedTo),
      sepvecs(theSepvecs),
      factorsWithRV(theFactorsWithRV) {

  for (unsigned idx = 0; idx < theFactorPtrs.size(); idx++){
    if (deepCopy){
      factorPtrs[idx] = rcptr<Factor>( theFactorPtrs[idx]->copy() );
    } // if
    else{
      factorPtrs[idx] = theFactorPtrs[idx];
    } // else
  } // for idx

  if ( !factorsWithRV.size() ) {
    for (unsigned fIdx = 0; fIdx < factorPtrs.size(); fIdx++){
      const Factor& theFactor(*factorPtrs[fIdx]);
      unsigned nVars = theFactor.noOfVars();
      for (unsigned varNo = 0; varNo < nVars; varNo++){
        factorsWithRV[theFactor.getVar(varNo)].push_back(fIdx);
      } // for
    } // for

  } // if

} // class specific ctor


// defunct ctor
ClusterGraph::ClusterGraph(
  const std::vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy)
    : ClusterGraph(LTRIP, theFactorPtrs, observed, deepCopy) {
} // defunct ctor

// default constructor
ClusterGraph::ClusterGraph(){} // default constructor

/*
// copy constructor
ClusterGraph::ClusterGraph(const ClusterGraph& g)
    : factorPtrs(g.factorPtrs),  // shallow
      linkedTo(g.linkedTo),
      sepvecs(g.sepvecs),        // shallow
      factorsWithRV(g.factorsWithRV){
} // copy constructor
*/

ClusterGraph::~ClusterGraph(){} // destructor

//========================= Operators and Conversions =========================

/*
ClusterGraph& ClusterGraph::operator=(const ClusterGraph& g){
  if (this != &g){
    this->~ClusterGraph();              // Destroy existing
    new ( static_cast<void *>(this) ) ClusterGraph(g); // do an inplace construction
  } // if

  return *this;
} // operator=
*/

//================================= Iterators =================================

//============================= Exemplar Support ==============================

const string& ClusterGraph::isA() const {
  static const string CLASSNAME("ClusterGraph");
  return CLASSNAME;
} // isA

ClusterGraph* ClusterGraph::copy() const {
  return new ClusterGraph(*this);
} // copy

//=========================== Inplace Configuration ===========================

// class specific
unsigned
ClusterGraph::classSpecificConfigure(
  const vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  const bool deepCopy){
  this->~ClusterGraph(); // Destroy existing
  new(this) ClusterGraph(theFactorPtrs, observed, deepCopy); // and do an inplace construction
  return 1;
} // classSpecificConfigure

istream& ClusterGraph::txtRead(istream& file){ return file; } // txtRead

ostream& ClusterGraph::txtWrite(ostream& file) const {
  // first show the factors
  for (unsigned i = 0; i < factorPtrs.size(); i++){
    file << "Factor_" << i << ": " << *factorPtrs[i] << endl;
  } // for
  // Now show what is linked to/from what with what sepvec
  for (unsigned from = 0; from < linkedTo.size(); from++){
    file << from << " <==> ";
    for (auto to : linkedTo[from]){
      // bummer, the map[] is non-const, have to use find
      file << to << " ( sepvec: "
           << sepvecs.find( Idx2(from, to) )->second << ") ";
    } // for
    file << endl;
  } // for
  return file;
} // txtWrite

//===================== Required Virtual Member Functions =====================

//====================== Other Virtual Member Functions =======================

//======================= Non-virtual Member Functions ========================

// When a factor is fully observed it will be removed as a whole.
void
ClusterGraph::observeAndReduce(
  vector< rcptr<Factor> >& theFactorPtrs,
  const std::map<emdw::RVIdType, AnyType>& observed,
  map< emdw::RVIdType, vector<unsigned> >& theFactorsWithRV){

  if ( !observed.size() ){return;} // if

  // cout << __FILE__ << __LINE__ << endl;
  // for (unsigned idx = 0; idx < observedVars.size(); idx++){
  //   cout << observedVars[idx]<< ":" << unsigned(observedVals[idx]) << " " << flush;
  // } // for
  // cout << endl;

  if ( !theFactorsWithRV.size() ){ // empty, have to search for the RV's
    for (unsigned fIdx = 0; fIdx < theFactorPtrs.size(); fIdx++){
      const Factor& theFactor(*theFactorPtrs[fIdx]);
      unsigned nVars = theFactor.noOfVars();
      for (unsigned varNo = 0; varNo < nVars; varNo++){
        theFactorsWithRV[theFactor.getVar(varNo)].push_back(fIdx);
      } // for
    } // for
  } // if

  for (auto itr = observed.begin(); itr != observed.end(); itr++){

    emdw::RVIdType varId = itr->first;
    const vector<unsigned>& factorIndices(theFactorsWithRV[varId]);
    // cout << "RV " << varId << " occurs in factors " << factorIndices << endl;
    // reduce all the factors containing varId
    for (unsigned i = 0; i < factorIndices.size(); i++){
      auto fIdx = factorIndices[i];
      // cout << "Reducing factor " << fIdx << " by observing " << varId << endl;
      theFactorPtrs[fIdx] =
        theFactorPtrs[fIdx]->
        observeAndReduce(emdw::RVIds{varId},
                         emdw::RVVals{itr->second});
    } // for

    // and now we can zap the observed var from our list
    theFactorsWithRV.erase(varId);

  } // for

  // check for and remove fully observed factors, normalize the others
  auto sz = theFactorPtrs.size();
  while (sz--) {
    if ( !theFactorPtrs[sz]->getVars().size() ){ // empty factor
      theFactorPtrs[sz] = 0;
    } // if
    else {
      theFactorPtrs[sz]->inplaceNormalize();
    } // else
  } // while

  vector<unsigned> old2new;
  repackFactors(theFactorPtrs, theFactorsWithRV, old2new);
  // cout << "EO observeAndReduce\n";
  // cout << "theFactorPtrs:\n" << theFactorPtrs.size() << endl;
  // for (auto elem : theFactorPtrs){
  //   cout << elem << ": ";
  //   if (elem){ cout << elem->getVars(); }
  //   cout << endl;
  // } // for

} // observeAndReduce

void ClusterGraph::swapFactors(
  vector< rcptr<Factor> >& factors,
  size_t idx1,
  size_t idx2,
  vector<size_t>& toNew,
  vector<size_t>& toOld) {
  swap(factors[idx1], factors[idx2]);
  swap(toOld[idx1], toOld[idx2]);
  toNew[toOld[idx1]] =  idx1;
  toNew[toOld[idx2]] =  idx2;
} // newSwapFactors

void ClusterGraph::absorbSubsetFactors(
  std::vector< rcptr<Factor> >& theFactorPtrs,
  std::vector<std::map<size_t, emdw::RVIds> >& theIntersections,
  std::map<size_t, std::set<size_t> >& theReachableFrom,
  std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV) {

  size_t nFactors = theFactorPtrs.size();
  theIntersections.clear();
  theReachableFrom.clear();
  //theFactorsWithRV.clear();

  // For a large number of factors it can be very expensive to
  // directly find all theIntersections. Instead, we will make little
  // pools for which we know that all factors in that pool have at
  // least one variable in common. So lets first find out which
  // factors have which RVs.
  if ( !theFactorsWithRV.size() ){ // empty, have to search for the RV's
    for (unsigned k = 0; k < nFactors; k++) {
      for (unsigned m = 0; m < theFactorPtrs[k]->noOfVars(); m++) {
        theFactorsWithRV[theFactorPtrs[k]->getVar(m)].push_back(k);
      } // for
    } // for
  } // if

  // lets complete theReachableFrom
  std::map<size_t, std::set<size_t> > oldReachFrom;
  for (auto rv : theFactorsWithRV) {  // rv.first is the RV
    for (auto itr1 = rv.second.begin(); itr1 != rv.second.end(); itr1++) { // itr1 points to the first factor
      size_t m = *itr1;
      for (auto itr2 = std::next(itr1); itr2 != rv.second.end(); itr2++) {
        size_t n = *itr2;
        oldReachFrom[m].insert(n);
        oldReachFrom[n].insert(m); // we'll need this to find all the links to subsets
      } // for
    } // for
  } // for

  vector<size_t> toNew(nFactors); // toNew[m] where is that original factor m located now?
  vector<size_t> toOld(nFactors); // toOld[m] where was this new factor m originally?
  for (size_t idx = 0; idx < nFactors; idx++) {
    toNew[idx] = toOld[idx] = idx;
  } // for idx

  ///cout << "Factors:\n";
  ///for (size_t k = 0; k < nFactors; k++) {
  ///  cout << toOld[k] << "/" << k << ": " << theFactorPtrs[k]->getVars() << endl;
  ///} // for
  ///cout << "Reachable from:\n";
  ///for (auto mm : oldReachFrom) {
  ///  cout << mm.first << ": ";
  ///  for (auto n : oldReachFrom[mm.first]) {
  ///    cout << n << " ";
  ///  } // for
  ///  cout << endl;
  ///} // for

  size_t r = 0;  // indexes the factor currently at r
  theIntersections.resize(nFactors);

  // TRICK: We get rid of subset factors by assimilating them into
  // their supersets and then swapping the subset out with the current
  // last factor in the array, and then decrement the length of the
  // array. To avoid having to find and change previous intersection
  // indices, the calculation below is arranged such that the first
  // intersection index will never change, and the second intersection
  // index is kept at its original pre-swapping value, while tracking
  // separately where it is newly located. After all is done, these
  // second indices are updated to their final values.

  while (r < nFactors) {   // r is new index
    auto rOld = toOld[r];  // it was originally at rOld

    for (auto cOld : oldReachFrom[rOld]) { // the (old) factors that r originally could reach
      auto c = toNew[cOld]; // the new index for those factors
      if (c < r) {continue;} // we only record theIntersections with r < c
     ///cout << "\nConsidering factors " << rOld << "/" << r << " and " << cOld << "/" << c << endl;
      if (c >= nFactors) {
       ///cout << "Nope, " << cOld << "/" << c << " does not exist any more\n";
        continue;
      } // done

      vector<RVIdType> intersect;
      const RVIds& rSet(theFactorPtrs[r]->getVars());
      const RVIds& cSet(theFactorPtrs[c]->getVars());
      set_intersection(rSet.begin(),rSet.end(),cSet.begin(),cSet.end(),
                       std::inserter(intersect,intersect.begin()));
      auto iSz = intersect.size();
      if (iSz) {

        // cSet is a subset of rSet, therefore absorb it into
        // rSet. Move the last factor into position of cSet and
        // continue from c
        if (cSet.size() == iSz and !theFactorPtrs[c]->keepSeparate) {
          ///cout << cOld << "/" << c << " is a subset of " << rOld << "/" << r << endl;
          ///cout << "Absorbing factor " << cOld << "/" << c << ": " << *theFactorPtrs[c];
          ///cout << " into factor " << rOld << "/" << r; // << ": " << *theFactorPtrs[r];
          ///cout << endl;
          // absorb the subset at c into r
          theFactorPtrs[r]->inplaceAbsorb(theFactorPtrs[c]);
          theFactorPtrs[r]->inplaceNormalize();

          // one factor less after the absorb
          nFactors--;

          // swap the subset factor out to last position
          swapFactors(theFactorPtrs, c, nFactors, toNew, toOld);
          ///cout << "Swop factors " << cOld << "/" << c << " and " << nFactors
          ///     << " (LEN " << nFactors << ")" << endl;
          ///cout << "Factors:\n";
          ///for (size_t k = 0; k < nFactors; k++) {
          ///  cout << toOld[k] << "/" << k << ": " << theFactorPtrs[k]->getVars() << endl;
          ///} // for

          // continue from current c
          continue;
        } // if

        // rSet is a subset of cSet. Swap them and restart everything from this new r
        else if (rSet.size() == iSz and !theFactorPtrs[r]->keepSeparate) {
          //cout << r << " is a subset of " << c << endl;
          swapFactors(theFactorPtrs, r, c, toNew, toOld);
          ///cout << "Swop factors " << rOld << "/" << r
          ///     << " and " << cOld << "/" << c << " (LEN " << nFactors << ")" << endl;
          ///cout << "Factors:\n";
          ///for (size_t k = 0; k < nFactors; k++) {
          ///  cout << toOld[k] << "/" << k << ": " << theFactorPtrs[k]->getVars() << endl;
          ///} // for
          theIntersections[r].clear(); // start afresh
          r--; break; // restart with the new factor at r
        } // else if

        // NOTE: An important part of the trick in the above
        // arrangement is that we never swap theIntersections - that will
        // be expensive. The price to pay is that we need to
        // recalculate some theIntersections if r is an subset of
        // c. Fortunately this should be the exception.
       ///cout << "Record their intersection:" << intersect << endl;
        theIntersections[r][toOld[c]] = intersect; //iSz;
      } // if iSz

    } // for c
    r++;
  } // while r

  theFactorPtrs.resize(nFactors);
  theIntersections.resize(nFactors);

  //theReachableFrom.resize(nFactors);

  // the first index of theIntersections refers to new position, the
  // second to old position. Now translate those second indices to
  // also refer to new position, and then remove any links to factors
  // with idx >= nFactors.
  for (size_t r = 0; r < theIntersections.size(); r++) {
    map<size_t, emdw::RVIds> isects;
    for (auto toOldDest : theIntersections[r]) {
      size_t toNewDest = toNew[toOldDest.first];
      //cout << r << "->" << toOldDest.first << " becomes "  << r << "->" << toNewDest << endl;
      if (toNewDest < nFactors) {
        isects[toNewDest] = std::move(toOldDest.second);
        theReachableFrom[r].insert(toNewDest);
        theReachableFrom[toNewDest].insert(r);
      } // if
    } // for
    if ( isects.size() ) {
      theIntersections[r] = isects;
    } // if
  } // for

  ///cout << "Factors:\n";
  ///for (size_t k = 0; k < nFactors; k++) {
  ///  cout << k << ": " << theFactorPtrs[k]->getVars() << endl;
  ///} // for
  ///cout << "TheIntersections\n";
  ///for (auto k = 0; k < theIntersections.size(); k++) {
  ///  for (auto toI : theIntersections[k]) {
  ///    cout << k << "->" << toI.first << ": " << toI.second << endl;
  ///  } // for
  ///} // for

  ///cout << "RRReachable from:\n";
  ///for (auto mm : theReachableFrom) {
  ///  cout << mm.first << ": ";
  ///  for (auto n : theReachableFrom[mm.first]) {
  ///    cout << n << " ";
  ///  } // for
  ///  cout << endl;
  ///} // for
  ///exit(-1);

  // The factor numbers in theFactorsWithRV still refer to the old factor
  // positions. Lets now fix those.

  for (auto& rvNr : theFactorsWithRV) {
    size_t newSz = rvNr.second.size();
    size_t k = 0;
    while (k < newSz){
      unsigned newIdx = toNew[ rvNr.second[k] ];
      if (newIdx < nFactors) {
        rvNr.second[k] = newIdx;
        k++;
      } // if
      else {
        newSz--;
        swap(rvNr.second[k], rvNr.second[newSz]);
      } // else
    } // while
    rvNr.second.resize(newSz);
  } // for
  for (auto& rvNr : theFactorsWithRV) { // maybe inefficient?
    std::sort( rvNr.second.begin(), rvNr.second.end() );
  } // for

  // cout << "Factors:\n";
  // for (auto el : theFactorPtrs) {
  //   cout << *el << endl;
  // } // for

  // cout << "TheIntersections:\n";
  // for (auto m = 0; m < nFactors; m++) {
  //   for (auto toN : theIntersections[m]) {
  //     cout << m << "->" << toN.first << ": " << toN.second << endl;
  //   } // for
  // } // for

  // cout << "Reachable from:\n";
  // for (auto mm : theReachableFrom) {
  //   cout << mm.first << ": ";
  //   for (auto n : theReachableFrom[mm.first]) {
  //     cout << n << " ";
  //   } // for
  //   cout << endl;
  // } // for

} // absorbSubsets

// testing code for absorbSubsets
/*
  {
    typedef unsigned T;
    typedef DiscreteTable<T> DT;
    rcptr< vector<T> > bin (new vector<T>{0,1});
    vector< rcptr<Factor> > factors;
    rcptr<DT> dtPtr;

    dtPtr = uniqptr<DT> ( new DT({4,6,7},{bin,bin,bin}) ); // 0
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    dtPtr = uniqptr<DT> ( new DT({4,5},{bin,bin}) ); //1
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    dtPtr = uniqptr<DT> ( new DT({0,1,2},{bin,bin,bin}) ); //2
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    dtPtr = uniqptr<DT> ( new DT({0,1},{bin,bin}) ); //3
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    dtPtr = uniqptr<DT> ( new DT({0,1,2,3,7},{bin,bin,bin,bin,bin}) ); //4
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    dtPtr = uniqptr<DT> ( new DT({4,5,1},{bin,bin,bin}) ); //5
    dtPtr->makeDefaultZero();
    dtPtr->addNoise(1.0);
    dtPtr->Factor::inplaceNormalize();
    factors.push_back(dtPtr);

    // fliplr
    // for (size_t k = 0; k < factors.size()/2; k++) {
    //   swap(factors[k], factors[factors.size()-1-k]);
    // } // for

    rcptr<Factor> f145 = factors[5]->absorb(factors[1])->normalize();
    rcptr<Factor> f01237 = factors[4]->absorb(factors[3])->absorb(factors[2])->normalize();
    rcptr<Factor> f467 = factors[0];

    vector< map<size_t, emdw::RVIds> > theIntersections;
    std::map<size_t, std::set<size_t> > theReachableFrom;
    std::map< emdw::RVIdType, std::vector<unsigned> > theFactorsWithRV;

    absorbSubsetFactors(factors, theIntersections, theReachableFrom, theFactorsWithRV);

    if (*factors[1] != *f145) {
      cout << "Factors differ\n";
      cout << *factors[1] << endl << *f145 << endl;
    } // if
    if (*factors[2] != *f01237) {
      cout << "Factors differ\n";
      cout << *factors[1] << endl << *f01237 << endl;
    } // if
    if (*factors[0] != *f467) {
      cout << "Factors differ\n";
      cout << *factors[0] << endl << *f467 << endl;
    } // if

    cout << "theIntersections ought to be 4, 7 and 1" << endl;
    for (size_t m = 0; m < theIntersections.size(); m++) {
      //if ( theIntersections[m].size() ){
      for (auto n : theIntersections[m]) {
        cout << m << "^" << n.first << ": " << n.second << endl;
      } // for
        //} // if
    } // for
    exit(-1);
  }

 */

void ClusterGraph::buildJunctionTree(
  std::vector< rcptr<Factor> >& theFactorPtrs,
  std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV) {

  //helper function if i,j and j,i mirroring is needed
  auto get_ij_ji = [](unsigned i, unsigned j) -> vector<pair<unsigned, unsigned>>{
    vector<pair<unsigned, unsigned>> ij_ji;
    ij_ji.push_back(pair<unsigned, unsigned>(i,j));
    ij_ji.push_back(pair<unsigned, unsigned>(j,i));
    return ij_ji;
  };

  //************************************
  // Helper variables to keep track of everything (different idexes than final versions)
  vector<rcptr<Factor>> oldFactorPtrs = theFactorPtrs;
  theFactorPtrs.clear();

  vector<RVIdType> vars;                         //List of variables [v0, v1, v2, v3, ...]
  vector<vector<RVIdType >> factorVars;          //Vars in factor: f0:[v0,v2,v2], f1:[v0,v1,v2], f2:[v1,v2,v5]
  map<unsigned, set<unsigned> > junctFactors;    //Show which theFactorPtrs are in a junct. Init w/ j0:{f0}, j1:{f1}, j2:{f2}
  map<unsigned, set<RVIdType> > junctVars;       //Show vars in junction, eg: j0:{v0,v1,v2}, j1:{v1,v2,v5}, ...
  map<RVIdType, set<unsigned> > tmpVarToJunction;//Show junction contains vars, eg: v0:{j0,j3}, v1:{j0,j1,j2}, ...
  map<RVIdType, set<unsigned> > varToFactor;     //Show factor contains vars, eg: v0:{f0,f3}, v1:{f0,f1,f2}, ...
  map<unsigned, Idx2> sepsetLink;     //Sepsets are new theFactorPtrs that links two junctions: f10:(j0,j1), f11(j1,j2), ...
  auto sepsetPoint = unsigned(oldFactorPtrs.size());  //factorPtrs[sepsetPoint] onwards are sepsets.

  //***************************************
  // Section 1 - Building the junction tree
  //***************************************
  auto set_insert = [](vector<unsigned> &v, const unsigned &t){
    auto i = std::lower_bound(v.begin(), v.end(), t);
    if (i == v.end() || t < *i){
      v.insert(i, t);
    }
  };

  // Initiate the helper variables
  for (unsigned i = 0; i < oldFactorPtrs.size(); i++) {
    factorVars.push_back(oldFactorPtrs[i]->getVars());

    junctFactors[i] = set<unsigned>({i});
    junctVars[i] = set<RVIdType>();
    for (auto &var: oldFactorPtrs[i]->getVars()) {
      junctVars[i].insert(var);
    }
  }

  for (auto &fPtr: oldFactorPtrs) {
    for (auto &var: fPtr->getVars()) {
      set_insert(vars, var);
    }
  }

  for (auto &var: vars) {
    tmpVarToJunction[var] = set<unsigned>();
    varToFactor[var] = set<unsigned>();
  }

  for (unsigned fIdx = 0; fIdx < oldFactorPtrs.size(); fIdx++) {
    for (auto &var: oldFactorPtrs[fIdx]->getVars()) {
      tmpVarToJunction[var].insert(fIdx);
      varToFactor[var].insert(fIdx);
    }
  }

  //************************************
  // Calculate elimination order
  auto elemOrder = varElemOrder(oldFactorPtrs, vars);

  //************************************
  // Elimate vars, as a result creating junctions, sepsets, and links
  for (auto &varE: elemOrder) {

    auto junctionsContainVarE = tmpVarToJunction[varE]; //TODO: Is copy? we need a copy!
    auto jIdx0 = *junctionsContainVarE.begin();
    auto sepIdx = unsigned(factorVars.size());

    for (auto &jIdxI: junctionsContainVarE) {

      //pseudo: tmpVarToJunction[var assoc/w jIdxI].replace(jIdxI, sepIdx)
      //(remove eliminated factors by replacing with their sepsets)
      for (auto &var: junctVars[jIdxI]) {
        auto iterI = tmpVarToJunction[var].find(jIdxI);
        PRLITE_ASSERT(iterI != tmpVarToJunction[var].end(),
               "Oops tmpVarToJunction and juncVars have different relationship");

        tmpVarToJunction[var].erase (iterI );
        tmpVarToJunction[var].insert(sepIdx);
      }

      //Complete the link (initial idx)--[sepset jIdxI]--(jIdx0), if oldFactorPtrs[jIdxI] a prev generated sepset?
      if (sepsetLink.count(jIdxI) > 0) {
        sepsetLink[jIdxI].second = jIdx0;
      }

      if (jIdx0 != jIdxI) {
        //Add jIdxI to junction number jIdx0
        junctFactors[jIdx0].insert(jIdxI);

        //pseudo: junctVars[jIdx0] <-add-- vars assoc/w jIdxI
        junctVars[jIdx0].insert(junctVars[jIdxI].begin(), junctVars[jIdxI].end());

        //Erase
        junctFactors.erase(jIdxI);
        junctVars.erase(jIdxI);
      }
    }

    //1: tmpVarToJunction.erase(varE)                             //remove eliminated var
    //2: sepsetLink[insert @ sepIdx] = (jIdx0, None)              //add the sepset link
    //3: junctFactors[insert @ sepIdx] = {sepIdx}                 //add temp junction containing this sepset
    //4: junctVars[insert @ sepIdx]   = junctVars[jIdx0] - {varE} //sepset vars
    //5: factorVars[insert @ sepIdx] =   ^                        //sepset vars (create new factor for sepset)

    //(1)
    tmpVarToJunction.erase(varE);

    if (tmpVarToJunction.size()) {//last round, no sepset
      //(2)
      sepsetLink[sepIdx] = Idx2(jIdx0, unsigned(-1));

      //(3)
      junctFactors[sepIdx] = set<unsigned>({sepIdx});

      //(4)
      set<RVIdType> sepsetVars = junctVars[jIdx0];
      auto iterRVElim = sepsetVars.find(varE);
      sepsetVars.erase(iterRVElim);
      junctVars[sepIdx] = sepsetVars;

      //(5)
      factorVars.push_back(vector<RVIdType>(sepsetVars.size()));
      std::copy(sepsetVars.begin(), sepsetVars.end(),
                factorVars[factorVars.size() - 1].begin());
    }
  }

  //***************************************
  // Section 2 - Absorbing subset nodes
  //***************************************
  map<unsigned, map<unsigned, unsigned> > juncToJuncAndSeps; //efficient representation for graph walkthrough
  set<unsigned> loopOverNodes;

  for (auto &_: sepsetLink){
    auto sepIdx = _.first;
    auto i = _.second.first;
    auto j = _.second.second;

    juncToJuncAndSeps[i][j] = sepIdx;
    juncToJuncAndSeps[j][i] = sepIdx;

    loopOverNodes.insert(i);
    loopOverNodes.insert(j);
  }
  sepsetLink.clear();

  while(not loopOverNodes.empty()) {
    // JADP ruan se fix begin
    auto loop_over_nodes_iter = loopOverNodes.begin();
    for (uint counter = 0; counter < loopOverNodes.size(); counter++) {
      uint i = *loop_over_nodes_iter;
      loop_over_nodes_iter = loopOverNodes.erase(loop_over_nodes_iter); //remove i from merge candidate loop
      // JADP ruan se fix eindig

      for (auto &_: juncToJuncAndSeps[i]) {
        auto j = _.first;
        auto sepIdx = _.second;

        for (auto &_: get_ij_ji(i, j)) { //test i->j and j->i
          auto i = _.first;
          auto j = _.second;

          if (junctVars[j].size() == factorVars[sepIdx].size()) {

            for (auto &_: juncToJuncAndSeps[j]) { //relink all k (j neighbours) to i
              auto k = _.first;
              auto sepK = _.second;

              if (k != i) {
                juncToJuncAndSeps[i][k] = sepK; // add i <-> k
                juncToJuncAndSeps[k][i] = sepK;

                juncToJuncAndSeps[k].erase(j); //remove j <-> k
                //juncToJuncAndSeps[j].erase(k);
              }
            }
            juncToJuncAndSeps[i].erase(j); //remove i <-> j
            juncToJuncAndSeps.erase(j); //remove j

            if (loopOverNodes.find(j) != loopOverNodes.end()) { //remove node j from merge candidate loop
              loopOverNodes.erase(loopOverNodes.find(j));
            }
            loopOverNodes.insert(i); //add i back as a merge candidate

            //let junct i merge all factors and vars from from junct j
            junctFactors[i].insert(junctFactors[j].begin(), junctFactors[j].end());
            junctFactors.erase(j);

            junctVars[i].insert(junctVars[j].begin(), junctVars[j].end());
            junctVars.erase(j);

            goto breakout; //Cleanly break out of nested loops - https://stackoverflow.com/questions/1257744
          }
        }
      }
    }
    breakout:;
  }

  //*************************************************
  // SECTION 3 Make new factors from the Junctions
  //*************************************************
  theFactorPtrs.clear();
  theFactorsWithRV.clear();
  sepvecs.clear();
  linkedTo.clear();
  //map<Idx2, RVIds> sepvecs;
  //vector<set<unsigned>> linkedTo;

  map<unsigned, unsigned> idxOldToNew; //helper for junction to fa
  map<RVIdType, rcptr<Factor> > vacuousVarCache;
  map<unsigned, set<RVIdType> > junctVacuousEntries;

  auto i = unsigned(-1);
  for(auto &_: junctFactors) {
    auto j = _.first;
    auto factorsIdx = _.second;

    i++;
    idxOldToNew[j] = unsigned(i);
  }

  for(auto &_: junctFactors){
    auto jIdx = _.first;
    auto factorsIdx = _.second;

    set<RVIdType > varsFacts;
    set<RVIdType > varsSeps;
    for (auto &fIdx: factorsIdx){
      for(auto &var: factorVars[fIdx]){
        if(fIdx<sepsetPoint){ //Not a sepset
          varsFacts.insert(var);
        }else{
          varsSeps.insert(var);
        }
      }
    }

    junctVacuousEntries[jIdx] = set<unsigned>();
    set_difference(varsSeps.begin(), varsSeps.end(), varsFacts.begin(), varsFacts.end(),
                   std::inserter(junctVacuousEntries[jIdx], junctVacuousEntries[jIdx].end()));
  }

  for(auto &_: junctFactors){
    auto jIdx = _.first;
    auto factorsIdx = _.second;

    auto i = unsigned(theFactorPtrs.size());

    for(auto &fIdx: factorsIdx) {
      if (fIdx < sepsetPoint) {
        if(theFactorPtrs.size() == i) { //first entry
          theFactorPtrs.push_back(oldFactorPtrs[fIdx]);
        }else {
          if (oldFactorPtrs[fIdx]->noOfVars() > theFactorPtrs[i]->noOfVars()) {
            theFactorPtrs[i] = oldFactorPtrs[fIdx]->absorb(theFactorPtrs[i]);
          }else {
            theFactorPtrs[i] = theFactorPtrs[i]->absorb(oldFactorPtrs[fIdx]);
          }
        }
      }
    }
    for(auto &var: junctVacuousEntries[jIdx]){
      if (not vacuousVarCache.count(var)) { // build cache
        vacuousVarCache[var] = rcptr<Factor>(oldFactorPtrs[*varToFactor[var].begin()]->vacuousCopy({var}, true)
        )->marginalize({var}, true);
      }
      if(theFactorPtrs.size() == i) { //first entry
        theFactorPtrs.push_back(vacuousVarCache[var]);
      }else {
        theFactorPtrs[i] = theFactorPtrs[i]->absorb(vacuousVarCache[var]);
      }
    }
  }

  for(unsigned i = 0; i < theFactorPtrs.size(); i++){
    linkedTo.push_back(set<unsigned>());
  }

  for(auto &iItem: juncToJuncAndSeps){
    for(auto &jItem: iItem.second){
      auto iOld = iItem.first;
      auto jOld = jItem.first;
      auto idxSep = jItem.second;

      auto i = idxOldToNew[iOld];
      auto j = idxOldToNew[jOld];
      linkedTo[i].insert(j);
      sepvecs[Idx2({i,j})] = factorVars[idxSep];
    }
  }

  for(unsigned i = 0; i < theFactorPtrs.size(); i++){
    for(auto j: theFactorPtrs[i]->getVars()){
      // JADP ruan fix logic error
      if(!theFactorsWithRV.count(j)){
        theFactorsWithRV[j] = vector<unsigned>({i});
      }
      else{
        set_insert(theFactorsWithRV[j], i);
      }
    }
  }


  //*************************************************
  // SECTION 4 Remove leftover junk (like empty connections)
  //*************************************************
  vector<Idx2> emptySeps;
  for(auto &item: sepvecs){
    if(not item.second.size()){
      emptySeps.push_back(item.first);
    }
  }

  for(auto &item: emptySeps){
    sepvecs.erase(item);

    auto iter = linkedTo[item.first].find(item.second);
    if (iter != linkedTo[item.first].end()){
      linkedTo[item.first].erase(iter);
    }

    iter = linkedTo[item.second].find(item.first);
    if (iter != linkedTo[item.second].end()){
      linkedTo[item.second].erase(iter);
    }

  }

  vector<int> emptyLinks;
  for(auto i=int(linkedTo.size())-1; i>=0; i--){
    if(linkedTo[i].size()==0){
      emptyLinks.push_back(i);
    }
  }
  for(auto &i: emptyLinks){
    linkedTo.erase(linkedTo.begin()+i);
  }

}

// old2new can be uninitialized on entry. It returns by reference to
// where a specific factor has relocated.
void ClusterGraph::repackFactors(
  std::vector< rcptr<Factor> >& theFactorPtrs,
  std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV,
  std::vector<unsigned>& old2new){

  // cout << __FILE__ << __LINE__ << endl;
  // cout << "theFactorPtrs:\n" << theFactorPtrs.size() << endl;
  // for (auto elem : theFactorPtrs){
  //   cout << elem << ": ";
  //   if (elem){ cout << elem->getVars(); }
  //   cout << endl;
  // } // for

  // first pack the pointers while also creating a map
  old2new.resize( theFactorPtrs.size() );
  unsigned oldIdx = 0;
  unsigned newSz = 0;
  while ( oldIdx < theFactorPtrs.size() ){
    if (theFactorPtrs[oldIdx]){ // a valid one
      old2new[oldIdx] = newSz;
      theFactorPtrs[newSz++] = theFactorPtrs[oldIdx++];
    } // if
    else { // whoops, an empty old factor
      old2new[oldIdx++] = theFactorPtrs.size(); // this one is a gonner
    } // else
  } // for idx
  theFactorPtrs.erase( theFactorPtrs.begin() + newSz, theFactorPtrs.end() );

  // cout << __FILE__ << __LINE__ << endl;
  // cout << "old2new:\n" << old2new << endl;
  // cout << __FILE__ << __LINE__ << endl;
  // cout << "theFactorPtrs:\n" << theFactorPtrs.size() << endl;
  // for (auto idx = 0; idx < theFactorPtrs.size(); idx++){
  //   cout << "factor_" << idx << " at: " << theFactorPtrs[idx] << endl;
  //   if (theFactorPtrs[idx]){
  //     cout << *theFactorPtrs[idx] << endl;
  //   } // if
  // } // for

  // lets tackle theFactorsWithRV
  std::map< emdw::RVIdType, std::vector<unsigned> > tmpFactorsWithRV = theFactorsWithRV;
  theFactorsWithRV.clear();
  for (auto idToFactor : tmpFactorsWithRV){
    emdw::RVIdType varId = idToFactor.first;
    std::vector<unsigned> factorIdx;
    // fortunately we did not reorder the factors
    for (auto oldFactorIdx : idToFactor.second){
      unsigned newFactorIdx = old2new[oldFactorIdx];
      if ( newFactorIdx < theFactorPtrs.size() ){
        factorIdx.push_back(newFactorIdx);
      } // if
    } // for
    theFactorsWithRV[varId] = factorIdx;
    // cout << __FILE__ << __LINE__ << endl;
    // cout << varId << " was: " << idToFactor.second << endl;
    // cout << varId << " now: " << factorIdx << endl;
  } // for

  // cout << "Compacted " << old2new.size() << " factors to " << theFactorPtrs.size() << endl;

} // repackFactors

void ClusterGraph::repackGraph(
  std::vector< rcptr<Factor> >& theFactorPtrs,
  std::map< emdw::RVIdType, std::vector<unsigned> >& theFactorsWithRV,
  std::vector< std::set<unsigned> >& theLinkedTo,
  std::map<Idx2, emdw::RVIds>& theSepvecs){

  // cout << __FILE__ << __LINE__ << endl;
  // cout << "theFactorPtrs:\n" << theFactorPtrs.size() << endl;
  // for (auto elem : theFactorPtrs){
  //   cout << elem << ": ";
  //   if (elem){ cout << elem->getVars(); }
  //   cout << endl;
  // } // for

  // first pack the pointers while also creating a map
  vector<unsigned> old2new;
  repackFactors(theFactorPtrs, theFactorsWithRV, old2new);

  // cout << "theLinkedTo:\n";
  // for (unsigned idx = 0; idx < theLinkedTo.size(); idx++){
  //   cout << idx << " was linked to " << theLinkedTo[idx] << endl;
  // } // for
  std::vector< std::set<unsigned> > tmpLinkedTo = theLinkedTo;
  theLinkedTo.clear();
  theLinkedTo.resize( theFactorPtrs.size() );

  // once again luckily the order is intact
  for (unsigned oldFrom = 0; oldFrom < tmpLinkedTo.size(); oldFrom++){
    unsigned newFrom = old2new[oldFrom];
    if ( newFrom > theFactorPtrs.size() ){continue;} // if
    for (unsigned oldTo : tmpLinkedTo[oldFrom]){
      unsigned newTo = old2new[oldTo];
      if ( newTo < theFactorPtrs.size() ){
        theLinkedTo[newFrom].insert(newTo);
      } // if
    } // for
  } // for
  // cout << __FILE__ << __LINE__ << endl;
  // for (unsigned idx = 0; idx < theLinkedTo.size(); idx++){
  //   cout << idx << " now linked to " << theLinkedTo[idx] << endl;
  // } // for

  // only the sepsets remaining
  // for (auto aSepvec : theSepvecs){
  //   cout << "was between " << aSepvec.first.first
  //        << " and " << aSepvec.first.second
  //        << ": " << aSepvec.second << endl;
  // } // for
  std::map<Idx2, emdw::RVIds> tmpSepvecs = theSepvecs;
  theSepvecs.clear();
  for (auto aSepvec : tmpSepvecs){
    unsigned newFrom = old2new[aSepvec.first.first];
    unsigned newTo = old2new[aSepvec.first.second];
    if ( newFrom < theFactorPtrs.size() && newTo < theFactorPtrs.size() ){
      theSepvecs[Idx2(newFrom, newTo)] = aSepvec.second;
    } // if
  } // for
  // for (auto aSepvec : theSepvecs){
  //   cout << "now between " << aSepvec.first.first
  //        << " and " << aSepvec.first.second
  //        << ": " << aSepvec.second << endl;
  // } // for

} // repackGraph

unsigned ClusterGraph::noOfFactors() const {
  return factorPtrs.size();
} // noOfFactors

void
ClusterGraph::exportToGraphViz(
  const string& graphname,
  const std::set<emdw::RVIdType>& idsToShow) const {

  vector<unsigned> clNr( factorPtrs.size() );
  string filename = graphname + ".dot";
  ofstream gvfile( filename.c_str() );
  if (!gvfile){
    PRLITE_FAIL("Could not open file '" << filename << "'");
  } // if

  // top stuff
  gvfile << "graph " << graphname << " {" << endl;
  gvfile << "\t" << "rankdir = LR;" << endl;
  gvfile << "\t" << "orientation = landscape;" << endl;
  // gvfile << "\t" << "concentrate = true;" << endl;
  gvfile << "\t" << "size = \"7.5,10.0\";" << endl;
  gvfile << endl;

  // do all the nodes
  set<size_t> clustersToShow;

  gvfile << "\t" << "node [shape=ellipse]" << endl;
  unsigned clCnt = 0;
  for (unsigned i = 0; i < factorPtrs.size(); i++){
    bool showCluster = true;
    if ( idsToShow.size() ){
      vector<emdw::RVIdType>::iterator peItr;
      vector<emdw::RVIdType> intersection( idsToShow.size() );
      peItr = set_intersection(
        factorPtrs[i]->getVars().begin(), factorPtrs[i]->getVars().end(),
        idsToShow.begin(), idsToShow.end(), intersection.begin() );

      if ( peItr == intersection.begin() ){
        showCluster = false;
      } // if
      else {
        clustersToShow.insert(i);
      }
    } // if

    if (showCluster){
      clNr[i] = clCnt++;
      //gvfile << "\t" << "n" << clNr[i] << " [label=\"" ;   // this one excludes the clusternr
      gvfile << "\t" << "n" << clNr[i] << " [label=\"" << i << ":"; // this line prepends the clusternr also
      for (unsigned j = 0; j < factorPtrs[i]->noOfVars(); j++){
        gvfile << factorPtrs[i]->getVar(j);
        if (j < factorPtrs[i]->noOfVars() - 1){
          gvfile << ",";
        } // if
      } // for j
      gvfile << "\"]\n";
    } // if
  } // for i

  // do all the sepsets and connect them between the nodes
  gvfile << "\n\t" << "node [shape=box]" << endl;
  unsigned cnt = 0; // numbers each sepset
  for (auto sep : sepvecs){
    if (sep.first.first < sep.first.second){
      // the sepset
      bool showSepSet = true;
      auto itr1 = clustersToShow.find(sep.first.first);
      auto itr2 = clustersToShow.find(sep.first.second);
      if ( !idsToShow.size() or
           itr1 != clustersToShow.end() or
           itr2 != clustersToShow.end() ){
        gvfile << "\t" << "s" << cnt << " [label=\"";
        for (unsigned j = 0; j < sep.second.size(); j++){
          gvfile << sep.second[j];
          if (j < sep.second.size() - 1){
            gvfile << ",";
          } // if
        } // for j
        gvfile << "\"";
        if ( idsToShow.size() and ( itr1 == clustersToShow.end() or itr2 == clustersToShow.end() ) ){
          gvfile << ", style=dashed";
        } // if
        gvfile << "]\n";
      } // if
      else {
        showSepSet = false;
      } // else

      // and connect it
      if (showSepSet){

        if ( !idsToShow.size() or
             clustersToShow.find(sep.first.first) != clustersToShow.end() ){
          gvfile << "\tn" << clNr[sep.first.first] << " -- " << "s" << cnt << endl;
        } // if

        if ( !idsToShow.size() or
             clustersToShow.find(sep.first.second) != clustersToShow.end() ){
          gvfile << "\ts" << cnt << " -- " << "n" << clNr[sep.first.second] << endl;
        } // if

        cnt++;
      } // if

    } // if
  } // for

  // and terminate things
  gvfile << "}" << endl;
  gvfile.close();

  cout << "Now use 'dot -Tpdf " << filename << " > " << filename
       << ".pdf' to generate a pdf picture\n";

} // exportToGraphViz

void
ClusterGraph::plotNodeStrengths(
  const string& graphname) const {

  vector<unsigned> clNr( factorPtrs.size() );
  vector<unsigned> clStrengths( factorPtrs.size() );
  for (auto& el : clStrengths) {el = 0;} // fo

  string filename = graphname + ".dot";
  ofstream gvfile( filename.c_str() );
  if (!gvfile){
    PRLITE_FAIL("Could not open file '" << filename << "'");
  } // if

  // top stuff
  gvfile << "graph " << graphname << "_strengths {" << endl;
  gvfile << "\t" << "rankdir = LR;" << endl;
  gvfile << "\t" << "orientation = landscape;" << endl;
  // gvfile << "\t" << "concentrate = true;" << endl;
  gvfile << "\t" << "size = \"7.5,10.0\";" << endl;
  gvfile << endl;

  // do all the sepsets and add their sizes to the adjacent nodes. connect them between the nodes
  gvfile << "\n\t" << "node [shape=box]" << endl;
  for (auto sep : sepvecs){
    if (sep.first.first < sep.first.second){
      // the sepset
      clStrengths[sep.first.first] += sep.second.size();
      clStrengths[sep.first.second] += sep.second.size();
    } // if
  } // for

  // create all the factor nodes
  gvfile << "\t" << "node [shape=ellipse]" << endl;
  for (unsigned i = 0; i < factorPtrs.size(); i++){
    gvfile << "\t" << "n" << i << " [label=\"" << clStrengths[i] <<  "\"]\n";
  } // for

  // do all the sepsets and connect their adjacent nodes
  gvfile << "\n\t" << "node [shape=box]" << endl;
  for (auto sep : sepvecs){
    if (sep.first.first < sep.first.second){
      gvfile << "\tn" << sep.first.first << " -- " << "n" << sep.first.second << endl;
    } // if
  } // for


  // and terminate things
  gvfile << "}" << endl;
  gvfile.close();

  cout << "Now use 'dot -Tpdf " << filename << " > " << filename
       << ".pdf' to generate a pdf picture\n";

} // plotNodeStrengths


//================================== Friends ==================================

/// output
std::ostream& operator<<( std::ostream& file,
                          const ClusterGraph& cg ){
  return cg.txtWrite(file);
} // operator<<

/// input
// std::istream& operator>>( std::istream& file,
//                          clusterGraph& cg ){
//  return cg.txtRead(file);
//} // operator>>


//=============================== Non-Class helpers =============================

ostream& writeMessagesToText(
  ostream& file,
  const map< Idx2, rcptr<Factor> >& messages){

  for (auto msg : messages){
    file << msg.first.first << " => " << msg.first.second << endl;
    file << "Message: " << *msg.second << endl;
  } // for

  return file;
} // writeMessagesToText


emdw::RVIds varElemOrder(const std::vector<rcptr<Factor>> &theFactorPtrs,
                         const emdw::RVIds &theVars){
  auto vars = theVars;
  if(theVars.size()==0){
    for(auto &fPtr: theFactorPtrs){
      for(auto &var: fPtr->getVars()){
        sortedSetInsert(vars, var);
      }
    }
  }

  //building the MRF
  std::map<emdw::RVIdType, std::set<emdw::RVIdType> > varNeighs;
  for (auto &var: vars){ varNeighs[var] = std::set<emdw::RVIdType>(); }
  for (auto &fPtr: theFactorPtrs) {
    for (auto &var: fPtr->getVars()) {
      for (auto &varNeigh: fPtr->getVars()) {
        if(var!=varNeigh) {
          varNeighs[var].insert(varNeigh);
        }
      }
    }
  }

  //find the elimination order
  std::vector<std::pair<unsigned, unsigned>> inducedEdges_var;
  for(auto &varE: vars){
    auto nONeighs = unsigned(varNeighs[varE].size()+1);
    auto nOFullyConn = nONeighs*(nONeighs-1)/2; //n(n-1)/2 formula

    std::set<std::pair<unsigned, unsigned>> connCount({});

    for(auto &varI: varNeighs[varE]){
      for(auto &varJ: varNeighs[varI]){
        if (varNeighs[varE].count(varJ)){
          if(varI < varJ){
            connCount.insert(std::pair<unsigned, unsigned>(varI, varJ));
          }else{
            connCount.insert(std::pair<unsigned, unsigned>(varJ, varI));
          }
        }
      }
    }

    auto nOInduced = nOFullyConn-unsigned(connCount.size());
    inducedEdges_var.push_back(std::pair<unsigned, unsigned>(nOInduced, varE));
  }
  sort(inducedEdges_var.begin(), inducedEdges_var.end());

  emdw::RVIds elimOrder;
  for(auto &item: inducedEdges_var){
    elimOrder.push_back(item.second);
  }

  return elimOrder;
}


std::vector<std::pair<emdw::RVIds, std::vector<unsigned>>>
  varMarginalizationOrder(const std::vector<rcptr<Factor>> &theFactorPtrs){

  auto factorsWithRV = map<RVIdType, vector<unsigned>>();
  for(auto i=0u; i < unsigned(theFactorPtrs.size()); i++){
    for(auto &var: theFactorPtrs[i]->getVars()){
      auto iter = factorsWithRV.find(var);
      if(iter==factorsWithRV.end()){ //item not in
        factorsWithRV[var] = {i};
      }else{ //item in
        sortedSetInsert(iter->second, i);
      }
    }
  }

  auto vars = RVIds();
  for(auto &_: factorsWithRV){
    auto var = _.first;
    vars.push_back(var);
  }

  auto elimOrder = varElemOrder(theFactorPtrs, vars);
  auto marginalizationOrder = vector<pair<RVIdType, vector<unsigned>>>();

  for(auto &varE: elimOrder){
    marginalizationOrder.push_back(make_pair(varE, factorsWithRV[varE]));
  }

  for(size_t i=0; i < elimOrder.size(); i++){// go through vars
    auto junctIdx = marginalizationOrder[i].second[0];
    for(size_t _=1; _< marginalizationOrder[i].second.size(); _++){
      auto oldIdx = marginalizationOrder[i].second[_];

      for(size_t j=i+1; j < elimOrder.size(); j++){

        //Try moving oldIdx to end of vector
        auto iter = remove(marginalizationOrder[j].second.begin(),
                           marginalizationOrder[j].second.end(),
                           oldIdx);

        //Remove and replace? Cut oldIdx from vector and replace with juncIdx
        if(iter != marginalizationOrder[j].second.end()){
          //Removing end-of-vector oldIdx
          marginalizationOrder[j].second.erase(
              iter, marginalizationOrder[j].second.end());

          //Insert junction index
          sortedSetInsert(marginalizationOrder[j].second, junctIdx);
        }
      }
    }
  }

  //Have multiple var summations over the same factor as one line
  vector<pair<RVIds, vector<unsigned>>> margOrderCondense(
      {make_pair(RVIds({marginalizationOrder[0].first}), marginalizationOrder[0].second)}
  );
  auto j=0u;
  for(size_t i=1; i<marginalizationOrder.size(); i++){
    if(marginalizationOrder[i].second == margOrderCondense[j].second){
      sortedSetInsert(margOrderCondense[j].first, marginalizationOrder[i].first);
    }else{
      margOrderCondense.push_back(
          make_pair(RVIds({marginalizationOrder[i].first}), marginalizationOrder[i].second)
      );
      j++;
    }
  }

  return margOrderCondense;

}
