#include "treefactor.hpp"/*******************************************************************************


          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Multi-variate TreeFactor PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "sortindices.hpp"

//emdw headers
#include "emdw.hpp"
#include "treefactor.hpp"
#include "vecset.hpp"
#include "clustergraph.hpp"
#include "lbu_cg.hpp"
#include "hashfunctions.hpp"

//std headers
#include <functional>


using namespace std;
using namespace emdw;

/**************************The TreeFactor class members***********************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
void TreeFactor::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new TreeFactor_Marginalize() );
  } // if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new TreeFactor_InplaceNormalize() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new TreeFactor_Normalize() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new TreeFactor_InplaceAbsorb() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new TreeFactor_Absorb() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new TreeFactor_InplaceCancel() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new TreeFactor_Cancel() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new TreeFactor_ObserveAndReduce() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper = uniqptr<InplaceDamper>( new TreeFactor_InplaceDamping() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new TreeFactor_Sampler() );
  } // if

} // supplementWithDefaultOperators

// default constructor
TreeFactor::TreeFactor(
  const rcptr <ClusterGraph>& theGraphPtr,
  bool preCalibrated,
  const rcptr<Marginalizer>& theMarginalizer,
  const rcptr<InplaceNormalizer>& theInplaceNormalizer,
  const rcptr<Normalizer>& theNormalizer,
  const rcptr<InplaceAbsorber>& theInplaceAbsorber,
  const rcptr<Absorber>& theAbsorber,
  const rcptr<InplaceCanceller>& theInplaceCanceller,
  const rcptr<Canceller>& theCanceller,
  const rcptr<ObserveAndReducer>& theObserveAndReducer,
  const rcptr<InplaceDamper>& theInplaceDamper,
  const rcptr<Sampler>& theSampler)
    : graphPtr(theGraphPtr),
      marginalizer(theMarginalizer),
      inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper),
      sampler(theSampler) {

  supplementWithDefaultOperators();

  if(not theGraphPtr){
    graphPtr = make_shared<ClusterGraph>();
  }

  vars.clear();
  for(auto &item: graphPtr->factorsWithRV){
      sortedSetInsert(vars, item.first);
  }

  if(not preCalibrated){
    calibrateteWithMessagePassing();
  }

} // default constructor


// default constructor
TreeFactor::TreeFactor(
    const vector<rcptr<Factor>>& theFactorPtrs,
    bool preCalibrated,
    const rcptr<Marginalizer>& theMarginalizer,
    const rcptr<InplaceNormalizer>& theInplaceNormalizer,
    const rcptr<Normalizer>& theNormalizer,
    const rcptr<InplaceAbsorber>& theInplaceAbsorber,
    const rcptr<Absorber>& theAbsorber,
    const rcptr<InplaceCanceller>& theInplaceCanceller,
    const rcptr<Canceller>& theCanceller,
    const rcptr<ObserveAndReducer>& theObserveAndReducer,
    const rcptr<InplaceDamper>& theInplaceDamper,
    const rcptr<Sampler>& theSampler)
    : marginalizer(theMarginalizer),
      inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper),
      sampler(theSampler) {

  supplementWithDefaultOperators();

  auto _ = std::map<emdw::RVIdType, AnyType>(); //dummy
  graphPtr = make_shared<ClusterGraph>(
      ClusterGraph::JTREE,
      theFactorPtrs,
      _,
      false
  );

  vars.clear();
  for(auto &item: graphPtr->factorsWithRV){
    sortedSetInsert(vars, item.first);
  }

  if(not preCalibrated){
    calibrateteWithMessagePassing();
  }

} // default constructor


// copy constructor. We provide this to make sure we do not have gLinear deep/shallow issues
TreeFactor::TreeFactor(const TreeFactor& p)
     :vars(p.vars),
      marginalizer(p.marginalizer),
      inplaceNormalizer(p.inplaceNormalizer),
      normalizer(p.normalizer),
      inplaceAbsorber(p.inplaceAbsorber),
      absorber(p.absorber),
      inplaceCanceller(p.inplaceCanceller),
      canceller(p.canceller),
      observeAndReducer(p.observeAndReducer),
      inplaceDamper(p.inplaceDamper),
      sampler(p.sampler) {
  // NOTE: We are assuming the operators are without memory, i.e. we
  // can apply the same operator to several factor instances without
  // any harm.

    graphPtr = make_shared<ClusterGraph>(
                    p.graphPtr->factorPtrs,
                    p.graphPtr->linkedTo,
                    p.graphPtr->sepvecs,
                    p.graphPtr->factorsWithRV,
                    true);
} // cpy ctor


TreeFactor::~TreeFactor() {
  // --noOfObjects;
  // std::cout << '*' << noOfObjects << '*' << std::flush;
} // destructor

//========================= Operators and Conversions =========================

bool TreeFactor::operator==(const TreeFactor& p) const {

  auto sizeLhs = graphPtr->noOfFactors();
  auto sizeRhs = p.graphPtr->noOfFactors();

  if(sizeLhs!=sizeRhs){
      return false;
  }
  if (vars != p.vars) {return false;} // if

  //TODO: consider equal junction trees build from different factors equal !!!!

  size_t i;
  for(i = 0; i < sizeLhs; i++){

    if(graphPtr->factorPtrs[i]->getVars() != p.graphPtr->factorPtrs[i]->getVars()){
      break;
    }
    if(graphPtr->factorPtrs[i] != p.graphPtr->factorPtrs[i]){
      return false;
    }
  }
  if(i == sizeLhs) { //forloop was completed without non-equalities
    return true;
  }

  //complete the rest, but make sure their vars match up -> using sneaky sort trick

  //Lineup factors according to their variables
  vector<pair<RVIds, size_t>> compLhs;
  vector<pair<RVIds, size_t>> compRhs;

  while (i < sizeLhs) {
    compLhs.push_back(make_pair(  graphPtr->factorPtrs[i]->getVars(), i));
    compRhs.push_back(make_pair(p.graphPtr->factorPtrs[i]->getVars(), i));
    i++;
  } // while
  sort(compLhs.begin(), compLhs.end());
  sort(compRhs.begin(), compRhs.end());

  for(size_t ii = 0; ii < compLhs.size(); ii++){
    auto i = compLhs[ii].second;
    auto j = compRhs[ii].second;
    if(graphPtr->factorPtrs[i] != graphPtr->factorPtrs[j]){
      return false;
    }
  }

  return true;

} // operator==

bool TreeFactor::operator!=(const TreeFactor& p) const {
  return !operator==(p);
} // operator !=

// boiler plate code
bool TreeFactor::isEqual(const Factor* rhsPtr) const {
  const TreeFactor* dwnPtr = dynamic_cast<const TreeFactor*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual


double TreeFactor::distanceFromVacuous() const {

  //TODO: this is only a heuristic
  double dist = 0;
  for (auto &fPtr: graphPtr->factorPtrs){
    dist+=fPtr->distanceFromVacuous();
  }

  return dist;
} // distanceFromVacuous

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

unsigned TreeFactor::classSpecificConfigure(
  const rcptr<ClusterGraph>& theGraphPtr,
  const rcptr<Marginalizer>& theMarginalizer,
  const rcptr<InplaceNormalizer>& theInplaceNormalizer,
  const rcptr<Normalizer>& theNormalizer,
  const rcptr<InplaceAbsorber>& theInplaceAbsorber,
  const rcptr<Absorber>& theAbsorber,
  const rcptr<InplaceCanceller>& theInplaceCanceller,
  const rcptr<Canceller>& theCanceller,
  const rcptr<ObserveAndReducer>& theObserveAndReducer,
  const rcptr<InplaceDamper>& theInplaceDamper,
  const rcptr<Sampler>& theSampler) {

  this->~TreeFactor(); // Destroy existing

  new(this) TreeFactor(
    theGraphPtr,
    true,
    theMarginalizer,
    theInplaceNormalizer,
    theNormalizer,
    theInplaceAbsorber,
    theAbsorber,
    theInplaceCanceller,
    theCanceller,
    theObserveAndReducer,
    theInplaceDamper,
    theSampler); // and do an inplace construction

  return 1;
} // classSpecificConfigure

// !!! this still needs fixing to handle the FactorOperators
std::ostream& TreeFactor::txtWrite(std::ostream& file) const {

  cout << "TreeFactor::txtWrite : FIXME factor operators will revert to defaults\n";
  file << "Vars: " << vars << "\n";

  return file;
} // txtWrite

// !!! this still needs fixing to handle the FactorOperators
std::istream& TreeFactor::txtRead(std::istream& file) {

  cout << "TreeFactor::txtRead : FIXME factor operators will revert to defaults\n";
  std::string dummy;
  emdw::RVIds theVars;

  file >> dummy;
  PRLITE_CHECK(file >> theVars, IOError, "Error reading vars RowVector.");


  // TODO BUG losing all the factor operators
  //classSpecificConfigure(theVars, tmpval);
  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
emdw::RVVals TreeFactor::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
uniqptr<Factor> TreeFactor::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
void TreeFactor::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
uniqptr<Factor> TreeFactor::absorb(const Factor* rhsPtr,
                                   FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
void TreeFactor::inplaceAbsorb(const Factor* rhsPtr,
                               FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
uniqptr<Factor> TreeFactor::cancel(const Factor* rhsPtr,
                                   FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
void TreeFactor::inplaceCancel(const Factor* rhsPtr,
                               FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
uniqptr<Factor> TreeFactor::marginalize(const emdw::RVIds& variablesToKeep,
                                        bool presorted,
                                        const Factor* peekAheadPtr,
                                        FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep, presorted, peekAheadPtr) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep, presorted, peekAheadPtr) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
uniqptr<Factor> TreeFactor::observeAndReduce(const emdw::RVIds& variables,
                                             const emdw::RVVals& assignedVals,
                                             bool presorted,
                                             FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals, presorted) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals, presorted) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
double TreeFactor::inplaceDampen(const Factor* rhsPtr,
                                 double df,
                                 FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

TreeFactor* TreeFactor::copy(
  const emdw::RVIds& newVars,
  bool presorted) const {

  if (! newVars.size() ) { // use all variables -> efficient identical copy
    return new TreeFactor(*this);
  }
  else{
    // sort the selected variables if necessary
    const emdw::RVIds *sortedVarsPtr = &newVars;
    emdw::RVIds tmpVars; // need it here to not go out of scope too soon
    std::vector<size_t> sorted; // need it here to not go out of scope too soon
    if (!presorted) {
      sorted = sortIndices(newVars, std::less<unsigned>());
      tmpVars = extract<unsigned>(newVars, sorted);
      sortedVarsPtr = &tmpVars;
    } // if
    const emdw::RVIds &sortedVars(*sortedVarsPtr);

    // find the subset
    std::vector<size_t> l2i, r2i;
    sortedIntersection(vars, sortedVars, l2i, r2i);

    // check that there are no stragglers
    PRLITE_ASSERT(r2i.size() == newVars.size(),
           "DiscreteTable::vacuousCopy : Some selected variables are not from the original factor");

    std::vector<rcptr<Factor>> rhsFactorPtrs;
    for(auto &fPtr: graphPtr->factorPtrs){
      auto fSelVars = sortedIntersection(fPtr->getVars(), sortedVars, l2i, r2i);
      if(fSelVars.size()) {
        rhsFactorPtrs.push_back(rcptr<Factor>(fPtr->copy(fSelVars, true)));
      }
    }
    auto rhsPtr =  new TreeFactor(rhsFactorPtrs,
                              true,   //TODO: Should we assume rhs is preCalibrated? Faster, yes! But error free?
                              marginalizer,
                              inplaceNormalizer,
                              normalizer,
                              inplaceAbsorber,
                              absorber,
                              inplaceCanceller,
                              canceller,
                              observeAndReducer,
                              inplaceDamper,
                              sampler);
    return rhsPtr;
  } // if newVars

} // copy

TreeFactor*
TreeFactor::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  if ( !selVars.size() ) { // use all variables

    std::vector<rcptr<Factor>> rhsFactorPtrs;
    for(auto &fPtr: graphPtr->factorPtrs){
      rhsFactorPtrs.push_back(rcptr<Factor>(fPtr->vacuousCopy()));
    }

    auto rhsPtr =  new TreeFactor(rhsFactorPtrs,
                              true,   //TODO: Should we assume rhs is preCalibrated? Faster, yes! But error free?
                              marginalizer,
                              inplaceNormalizer,
                              normalizer,
                              inplaceAbsorber,
                              absorber,
                              inplaceCanceller,
                              canceller,
                              observeAndReducer,
                              inplaceDamper,
                              sampler);
    return rhsPtr;

  } // else
  else {
    // sort the selected variables if necessary
    const emdw::RVIds *sortedVarsPtr = &selVars;
    emdw::RVIds tmpVars; // need it here to not go out of scope too soon
    std::vector<size_t> sorted; // need it here to not go out of scope too soon
    if (!presorted) {
      sorted = sortIndices(selVars, std::less<unsigned>());
      tmpVars = extract<unsigned>(selVars, sorted);
      sortedVarsPtr = &tmpVars;
    } // if
    const emdw::RVIds &sortedVars(*sortedVarsPtr);

    // find the subset
    std::vector<size_t> l2i, r2i;
    sortedIntersection(vars, sortedVars, l2i, r2i);

    // check that there are no stragglers
    PRLITE_ASSERT(r2i.size() == selVars.size(),
           "DiscreteTable::vacuousCopy : Some selected variables are not from the original factor");

    std::vector<rcptr<Factor>> rhsFactorPtrs;
    for(auto &fPtr: graphPtr->factorPtrs){
      auto fSelVars = sortedIntersection(fPtr->getVars(), sortedVars, l2i, r2i);
      if(fSelVars.size()) {
        rhsFactorPtrs.push_back(rcptr<Factor>(fPtr->vacuousCopy(fSelVars, true)));
      }
    }

    auto rhsPtr =  new TreeFactor(rhsFactorPtrs,
                              true,   //TODO: Should we assume rhs is preCalibrated? Faster, yes! But error free?
                              marginalizer,
                              inplaceNormalizer,
                              normalizer,
                              inplaceAbsorber,
                              absorber,
                              inplaceCanceller,
                              canceller,
                              observeAndReducer,
                              inplaceDamper,
                              sampler);
    return rhsPtr;
  }

} // vacuousCopy

unsigned TreeFactor::noOfVars() const {
  return vars.size();
} // noOfVars

emdw::RVIds TreeFactor::getVars() const {
  return vars;
} // getVars

emdw::RVIdType TreeFactor::getVar(unsigned varNo) const {
  return vars[varNo];
} // getVar

//====================== Other Virtual Member Functions =======================

double TreeFactor::potentialAt(
  const emdw::RVIds& theVars,
  const emdw::RVVals& theirVals,
  bool presorted) const{

  if (theVars.size() == 1) {presorted = true;}

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  const emdw::RVVals* sortedValsPtr = &theirVals;

  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  emdw::RVVals tmpVals; // need it here to not go out of scope too soon

  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    tmpVals = extract<AnyType>(theirVals,sorted);

    sortedVarsPtr = &tmpVars;
    sortedValsPtr = &tmpVals;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);
  PRLITE_ASSERT(sortedVars == vars, "The variables do not match: " << sortedVars << " vs " << vars);


  PRLITE_ASSERT(graphPtr->noOfFactors()==1, " TODO: TreeFactor::potentialAt currently limited to single factor graphs");

  return graphPtr->factorPtrs[0]->potentialAt(*sortedVarsPtr, *sortedValsPtr, true);

} // potentialAt

double TreeFactor::distance(const Factor* rhsPtr) const {

  const TreeFactor* dwnPtr = dynamic_cast<const TreeFactor*>(rhsPtr);
  PRLITE_ASSERT(dwnPtr, " TreeFactor::distance : rhs is not a TreeFactor");
  PRLITE_ASSERT(vars == dwnPtr->vars, " TreeFactor::distance : The variables do not match");

  auto sizeLhs = graphPtr->noOfFactors();
  auto sizeRhs = dwnPtr->graphPtr->noOfFactors();

  PRLITE_ASSERT(sizeLhs==sizeRhs, " TreeFactor::distance: Can currently only handle equal graph structures.");

  double d = 0;
  size_t i;
  for(i = 0; i < sizeLhs; i++){
    auto j = i;

    auto lhsVars = graphPtr->factorPtrs[i]->getVars();
    auto rhsVars = dwnPtr->graphPtr->factorPtrs[j]->getVars();

    if(lhsVars!=rhsVars){
      break;
    }
    d += graphPtr->factorPtrs[i]->distance(dwnPtr->graphPtr->factorPtrs[j]);
  }
  if(i == sizeLhs) { //forloop was completed without miss matches
    return d;
  }

  //complete the rest, but make sure their vars match up -> using sneaky sort trick

  //Lineup factors according to their variables
  vector<pair<RVIds, size_t>> compLhs;
  vector<pair<RVIds, size_t>> compRhs;

  while (i < sizeLhs) {
    compLhs.push_back(make_pair(        graphPtr->factorPtrs[i]->getVars(), i));
    compRhs.push_back(make_pair(dwnPtr->graphPtr->factorPtrs[i]->getVars(), i));
    i++;
  } // while
  sort(compLhs.begin(), compLhs.end());
  sort(compRhs.begin(), compRhs.end());

  for(size_t ii = 0; ii < compLhs.size(); ii++){
    auto i = compLhs[ii].second;
    auto j = compRhs[ii].second;

    auto lhsVars = graphPtr->factorPtrs[i]->getVars();
    auto rhsVars = dwnPtr->graphPtr->factorPtrs[j]->getVars();
    PRLITE_ASSERT(lhsVars==rhsVars, " TreeFactor::distance: Can currently only handle equal graph structures.");

    d += graphPtr->factorPtrs[i]->distance(dwnPtr->graphPtr->factorPtrs[j]);
  }
  return d;


} // distance

//======================= Parameterizations (non-virtual) =====================


//======================= Non-virtual Member Functions ========================

void TreeFactor::calibrateteWithMessagePassing(){
  map<Idx2, rcptr<Factor> > msgs;
  MessageQueue msgQ;
  loopyBU_CG(*graphPtr, msgs, msgQ, 0.0);
}


//Initiate the cache as empty
std::map<uint64_t, std::vector<std::pair<size_t, size_t>>>
    TreeFactor::factorsLineupCache = std::map<uint64_t, std::vector<std::pair<size_t, size_t>>>();

//Generate function lineup or use cache
vector<pair<size_t, size_t>> TreeFactor::getFactorLineup(const rcptr<ClusterGraph> &lhsGPtr,
                                                         const rcptr<ClusterGraph> &rhsGPtr){
  //Generate unique hash from the lhs rhs pointers
  uint64_t hash = 0;
  for(auto &fPtr: lhsGPtr->factorPtrs){
    hash = vectorToHash(fPtr->getVars(), hash);
  }
  hash += 514404; //S1M0N as arbitrary divider
  for(auto &fPtr: rhsGPtr->factorPtrs){
    hash = vectorToHash(fPtr->getVars(), hash);
  }

  auto cacheEntry = TreeFactor::factorsLineupCache.find(hash);

  //Got cache? Return
  if(cacheEntry != TreeFactor::factorsLineupCache.end()){
    return cacheEntry->second;

  //Ain't got cache? Try generating lineup, expensive, but once and for all
  }else{

    RVIds rhsVars;
    vector<int> _; //dummy
    for(auto &fPtr: rhsGPtr->factorPtrs) {
      rhsVars = sortedUnion(rhsVars, fPtr->getVars(), _, _);
    }

    RVIds varU = sortedUnion(rhsVars, vars, _, _);
    if(varU.size() < vars.size()){ //oops, new variables introduced
      return {};
    }

    std::vector<std::pair<size_t, size_t> > lineup;
    for(size_t j = 0; j < rhsGPtr->factorPtrs.size(); j++) {
      auto jVec = rhsGPtr->factorPtrs[j]->getVars();

      if(jVec.size() == 0){ //Empty factor, how? Connect arbitrarily to i=0
        lineup.push_back(make_pair(size_t(0),j));
        continue;
      }
      auto iList = lhsGPtr->factorsWithRV[jVec[0]];
      for (size_t v = 1; v < jVec.size(); v++) {
        iList = sortedIntersection(iList, lhsGPtr->factorsWithRV[jVec[v]], _, _);
      }

      if(iList.size()==0){
        return {};
      }
      lineup.push_back(make_pair(size_t(iList[0]),j)); //yeah found a lineup
    }

    //Fully lined-up? Then good to go into cache.
    TreeFactor::factorsLineupCache[hash] = lineup;
    return lineup;
  }
}


//Initialised the marginalization sequence cache
map<uint64_t, vector<pair<RVIds, vector<unsigned>>>>
  TreeFactor::marginalizationSequenceCache = map<uint64_t, vector<pair<RVIds, vector<unsigned>>>>();

//Generate marginalization sequence or use cache
vector<pair<RVIds,vector<unsigned>>> getMarginalizationSequence(
                                          const rcptr<ClusterGraph> &lhsGPtr){
  //Generate unique hash from the lhs rhs pointers
  uint64_t hash = 0;
  for(auto &fPtr: lhsGPtr->factorPtrs){
    hash = vectorToHash(fPtr->getVars(), hash);
  }

  auto cacheEntry = TreeFactor::marginalizationSequenceCache.find(hash);

  //Got cache? Return
  if(cacheEntry != TreeFactor::marginalizationSequenceCache.end()){
    return cacheEntry->second;

  //Ain't got cache? Generate marginalization Sequence
  }else{
    auto margOrder = varMarginalizationOrder(lhsGPtr->factorPtrs);
    TreeFactor::marginalizationSequenceCache[hash] = margOrder;
    return margOrder;
  }

};


//================================== Friends ==================================

/// output
template<typename T>
std::ostream& operator<<( std::ostream& file,
                          const TreeFactor& n ){
  return n.txtWrite(file);
} // operator<<

/// input
template<typename T>
std::istream& operator>>( std::istream& file,
                          TreeFactor& n ){
  return n.txtRead(file);
} // operator>>


//=============================== Helpers ==============================

bool isTFPtr(const rcptr<Factor> &fPtr){
  if(dynamic_pointer_cast<TreeFactor>(fPtr)){
    return true;
  }
  return false;
}

rcptr<TreeFactor> toTFPtr(const rcptr<Factor> &fPtr){
  auto theTF = dynamic_pointer_cast<TreeFactor>(fPtr);
  PRLITE_ASSERT(theTF, "Illegal rcptr<Factor> downcast to rcptr<TreeFactor>.");
  return theTF;
}

vector<rcptr<Factor>> flattenTFPtrs(const vector<rcptr<Factor>> &factorPtrs){
  vector<rcptr<Factor>> flattenedPtrs;

  //Helper recursive function to extract nonTreeFactors from TreeFactors
  function<void(const rcptr<Factor>)> doFlatten= [&flattenedPtrs, &doFlatten]
                                                 (const rcptr<Factor> &fPtr){
    if(isTFPtr(fPtr)){
      for(auto &subFPtr: toTFPtr(fPtr)->graphPtr->factorPtrs){
        doFlatten(subFPtr);
      }
    }else{
      flattenedPtrs.push_back(fPtr);
    }
  };

  for(auto &fPtr: factorPtrs) {
    doFlatten(fPtr);
  }

  return flattenedPtrs;
}

//======================= Factor operators ====================================


// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample   **************************************
// ***************************************************************************

// Identify object
const std::string& TreeFactor_Sampler::isA() const {
  static const std::string CLASSNAME("TreeFactor_Sampler");
  return CLASSNAME;
} // isA

emdw::RVVals TreeFactor_Sampler::sample(const TreeFactor* lhsPtr) {
  PRLITE_ASSERT(false, " TreeFactor_Sampler::sample: not implimented yet.");
} // sample


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** TreeFactor_Normalize ************************/
/*****************************************************************************/

// Identify object
const std::string& TreeFactor_InplaceNormalize::isA() const {
  static const std::string CLASSNAME("TreeFactor_InplaceNormalize");
  return CLASSNAME;
} // isA

void TreeFactor_InplaceNormalize::inplaceProcess(TreeFactor* lhsPtr) {

  TreeFactor& lhs(*lhsPtr);
  for(auto &fPtr: lhs.graphPtr->factorPtrs){
    fPtr->inplaceNormalize();
  }

} // inplaceProcess

// Identify object
const std::string& TreeFactor_Normalize::isA() const {
  static const std::string CLASSNAME("TreeFactor_Normalize");
  return CLASSNAME;
} // isA

Factor* TreeFactor_Normalize::process(const TreeFactor* lhsPtr) {

  // make a copy of the factor
  TreeFactor* fPtr = new TreeFactor(*lhsPtr);

  // inplace normalize it
  TreeFactor_InplaceNormalize ipNorm;
  try{
    ipNorm.inplaceProcess(fPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;
} // process


/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** TreeFactor_Absorb ********************************/
/*****************************************************************************/

// Identify object
const std::string& TreeFactor_InplaceAbsorb::isA() const {
  static const std::string CLASSNAME("TreeFactor_InplaceAbsorb");
  return CLASSNAME;
} // isA

/**
 * This operator inplace absorbs another TreeFactor message, ie
 * does an inplace factor product.
 */
void TreeFactor_InplaceAbsorb::inplaceProcess(TreeFactor* lhsPtr,
                                              const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate TreeFactor*
  const TreeFactor* rhsPtr = dynamic_cast<const TreeFactor*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " TreeFactor_InplaceAbsorb::inplaceProcess : rhs not a valid TreeFactor : " << typeid(*rhsFPtr).name() );
  const TreeFactor& rhs(*rhsPtr);
  TreeFactor& lhs(*lhsPtr);

  auto lineup = lhs.getFactorLineup(lhs.graphPtr, rhs.graphPtr);

  //Lineup attempt :)
  if(lineup.size() != 0){
    for(auto &ij: lineup){
      lhs.graphPtr->factorPtrs[ij.first]->inplaceAbsorb(rhs.graphPtr->factorPtrs[ij.second]);
    }
  }
  //No lineup found :( Oh well, new graph then.
  else{
    auto lhsFactorPtrs = lhs.graphPtr->factorPtrs;
    for (auto &fPtr: rhs.graphPtr->factorPtrs) {
      lhsFactorPtrs.push_back(fPtr);
    }

    auto _ = std::map<emdw::RVIdType, AnyType>(); //dummy
    lhs.graphPtr = make_shared<ClusterGraph>(
        ClusterGraph::JTREE,
        lhsFactorPtrs,
        _,
        false
    );
    std::vector<int> __; //dummy
    lhs.vars = sortedUnion(lhs.vars, rhs.vars, __, __);
  }

  lhs.calibrateteWithMessagePassing();

} // inplaceProcess

// Identify object
const std::string& TreeFactor_Absorb::isA() const {
  static const std::string CLASSNAME("TreeFactor_Absorb");
  return CLASSNAME;
} // isA

Factor* TreeFactor_Absorb::process(const TreeFactor* lhsPtr,
                               const Factor* rhsFPtr) {

  // make a copy of the factor
  TreeFactor* fPtr = new TreeFactor(*lhsPtr);

  // inplace normalize it
  TreeFactor_InplaceAbsorb ipAbsorb;
  try{
    ipAbsorb.inplaceProcess(fPtr,rhsFPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;

} // process


/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** TreeFactor_Cancel ********************************/
/*****************************************************************************/

// Identify object
const std::string& TreeFactor_InplaceCancel::isA() const {
  static const std::string CLASSNAME("TreeFactor_InplaceCancel");
  return CLASSNAME;
} // isA

/**
 * This operator inplace cancels another TreeFactor message, ie does an
 * inplace factor division. It is assumed that we will only cancel
 * messages that were previously absorbed.
 */
void TreeFactor_InplaceCancel::inplaceProcess(TreeFactor* lhsPtr,
                                          const Factor* rhsFPtr){

  // rhsFPtr must be an appropriate TreeFactor*
  const TreeFactor* rhsPtr = dynamic_cast<const TreeFactor*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " TreeFactor_InplaceCancel::inplaceProcess : rhs not a valid TreeFactor : " << typeid(*rhsFPtr).name() );
  const TreeFactor& rhs(*rhsPtr);
  TreeFactor& lhs(*lhsPtr);

  auto lineup = lhs.getFactorLineup(lhs.graphPtr, rhs.graphPtr);

  //Lineup attempt :)
  if(lineup.size() != 0){
    for(auto &ij: lineup){
      lhs.graphPtr->factorPtrs[ij.first]->inplaceCancel(rhs.graphPtr->factorPtrs[ij.second]);
    }
  }
  //No lineup found :(
  else{
    std::vector<int> _; //dummy
    auto vecU = sortedUnion(lhs.vars, rhs.vars, _, _);
    PRLITE_ASSERT(vecU.size()==lhs.vars.size(), " TreeFactor_InplaceCancel::inplaceProcess : rhs not subset of lhs");

    //TODO: create new graph
    PRLITE_ASSERT(false, " TreeFactor_InplaceCancel::inplaceProcess : currently rhs and lhs graphs must match (TODO: fix with creating new graph)");
  }

  lhs.calibrateteWithMessagePassing();

} // inplaceProcess

// Identify object
const std::string& TreeFactor_Cancel::isA() const {
  static const std::string CLASSNAME("TreeFactor_Cancel");
  return CLASSNAME;
} // isA

Factor* TreeFactor_Cancel::process(const TreeFactor* lhsPtr,
                               const Factor* rhsFPtr) {

  // make a copy of the factor
  TreeFactor* fPtr = new TreeFactor(*lhsPtr);

  // inplace normalize it
  TreeFactor_InplaceCancel ipCancel;
  try{
    ipCancel.inplaceProcess(fPtr,rhsFPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

  // and return the result
  return fPtr;
} // process


/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** TreeFactor_Marginalize ***************************/
/*****************************************************************************/

// Identify object
const std::string& TreeFactor_Marginalize::isA() const {
  static const std::string CLASSNAME("TreeFactor_Marginalize");
  return CLASSNAME;
} // isA

Factor* TreeFactor_Marginalize::process(const TreeFactor* lhsPtr,
                                        const emdw::RVIds& toKeep,
                                        bool presorted,
                                        const Factor* peekAheadPtr) {

  const TreeFactor& lhs(*lhsPtr);
  //int origDim = int(lhs.vars.size());
  PRLITE_ASSERT(toKeep.size(), "TreeFactor_Marginalize::process : toKeep empty, cannot marginalize everything out");

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &toKeep;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( toKeep, std::less<unsigned>() );
    tmpVars = extract<unsigned>(toKeep,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);


  auto margSeq = getMarginalizationSequence(lhs.graphPtr);
  auto factorPtrs = vector<rcptr<Factor>>();
  for(auto &fPtr: lhs.graphPtr->factorPtrs){
    factorPtrs.push_back(fPtr);
  }

  for(auto &_: margSeq){
    auto vars = _.first;
    auto fIdx = _.second;

    //if(fIdx.size()) {factorPtrs[fIdx[0]] = rcptr<Factor>(factorPtrs[fIdx[0]]->copy());}

    for(size_t i=1; i<fIdx.size(); i++){
      factorPtrs[fIdx[0]] = factorPtrs[fIdx[0]]->absorb(factorPtrs[fIdx[i]]);
      factorPtrs[fIdx[i]] = rcptr<Factor>(nullptr); //zap the factor
    }

    RVIds fToZap;
    std::set_difference(vars.begin(), vars.end(),
                        sortedVars.begin(), sortedVars.end(),
                        std::inserter(fToZap, fToZap.begin()));

    if(fToZap.size()) {
      factorPtrs[fIdx[0]] = factorPtrs[fIdx[0]]->marginalizeX(fToZap, true);
    }
  }

  auto fPtrsToKeep = vector<rcptr<Factor>>();
  for(auto &fPtr: factorPtrs){
    if(fPtr){ //don't include zapped-to-null
      if(fPtr->noOfVars()) { //don't include empties
        fPtrsToKeep.push_back(fPtr);
      }
    }
  }


  //Tree with no connections (by def these factors are independent)
  vector<set<unsigned>> linkedTo;
  map<Idx2, RVIds> sepvecs;
  map<RVIdType, vector<unsigned>> factorsWithRV;
  for(auto i=0u; i<unsigned(fPtrsToKeep.size()); i++) {
    for(auto &var: fPtrsToKeep[i]->getVars()) {
      factorsWithRV[var]={i};
    }
  }
  rcptr<ClusterGraph> margGraph = make_shared<ClusterGraph>(fPtrsToKeep, linkedTo, sepvecs, factorsWithRV, false);

  auto margPtr = new TreeFactor(margGraph,
                                true, //TODO: Should we assume rhs is preCalibrated? Faster yes! But error free?
                                lhs.marginalizer,
                                lhs.inplaceNormalizer,
                                lhs.normalizer,
                                lhs.inplaceAbsorber,
                                lhs.absorber,
                                lhs.inplaceCanceller,
                                lhs.canceller,
                                lhs.observeAndReducer,
                                lhs.inplaceDamper,
                                lhs.sampler);


  return margPtr;

} // process


/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** TreeFactor_ObserveAndReduce **********************/
/*****************************************************************************/

// Identify object
const std::string& TreeFactor_ObserveAndReduce::isA() const {
  static const std::string CLASSNAME("TreeFactor_ObserveAndReduce");
  return CLASSNAME;
} // isA

Factor* TreeFactor_ObserveAndReduce::process(const TreeFactor* lhsPtr,
                                         const emdw::RVIds& theVars,
                                         const emdw::RVVals& anyVals,
                                         bool presorted) {


  // and a convenient reference for the lhs factor
  const TreeFactor& lhs(*lhsPtr);
  //int origDim = lhs.vars.size();
  if (theVars.size() == 1) {presorted = true;}

  // if observed nothing, return as is
  if (!theVars.size()) {
    return lhs.copy();
  } // if

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  const emdw::RVVals* sortedValsPtr = &anyVals;

  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  emdw::RVVals tmpVals;

  std::vector<size_t> sorted; // need it here to not go out of scope too soon

  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    tmpVals = extract<AnyType>(anyVals,sorted);

    sortedVarsPtr = &tmpVars;
    sortedValsPtr = &tmpVals;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);
  const emdw::RVVals& sortedVals(*sortedValsPtr);

  /*
  // now i need to find out where the intersection between vars and
  // sortedVars lie. Any variables outside of these need not be
  // considered.
  std::vector<int> l2i; // which indices in lhs pot are observed (i.e. y)
  std::vector<int> r2i; // with corresponding indices in rhs pot
  std::vector<int> lni; // non-observed lhs variables ARE the new variable set (i.e. x)
  std::vector<int> rni; // remaining rhs variables (i.e. observed vars not part of this factor and thus ignored)
  emdw::RVIds varI = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);
  */

  std::vector<rcptr<Factor>> rhsFactorPtrs;
  for(auto &fPtr: lhs.graphPtr->factorPtrs){
    rhsFactorPtrs.push_back(rcptr<Factor>(fPtr->observeAndReduce(sortedVars, sortedVals, true)));
  }

  auto rhsPtr = new TreeFactor(rhsFactorPtrs,
                           true,   //TODO: Should we assume rhs is preCalibrated? Faster yes! But error free?
                           lhs.marginalizer,
                           lhs.inplaceNormalizer,
                           lhs.normalizer,
                           lhs.inplaceAbsorber,
                           lhs.absorber,
                           lhs.inplaceCanceller,
                           lhs.canceller,
                           lhs.observeAndReducer,
                           lhs.inplaceDamper,
                           lhs.sampler);

  return rhsPtr;


} // process

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************

// Identify object
const std::string& TreeFactor_InplaceDamping::isA() const {
  static const std::string CLASSNAME("TreeFactor_InplaceDamping");
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and old
// parameters - this linear combination is determined by the df parameter
double
TreeFactor_InplaceDamping::inplaceProcess(
  TreeFactor* newPtr,
  const Factor* oldPtr,
  double oldWeight) { // oldWeight is the weight for rhs

  //return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);

  TreeFactor& lhs(*newPtr);
  // rhsFPtr must be an appropriate TreeFactor*
  const TreeFactor* rhsPtr = dynamic_cast<const TreeFactor*>(oldPtr);
  PRLITE_ASSERT(rhsPtr, " TreeFactor_InplaceAbsorb::inplaceProcess : rhs not a valid TreeFactor : " << typeid(*oldPtr).name() );
  const TreeFactor& rhs(*rhsPtr);


  //We will transform the old message into the treestructure of the new one
  //Expensive, but will (at the very most) happen once per sepvec.
  auto structureOldToNew = [&lhs, &rhs](){
    PRLITE_ASSERT(false, " TreeFactor_InplaceAbsorb::inplaceProcess : old and new message graphs doesn't line up")
  };


  //TODO: we strongly assume the lhs and rhs treeFactors internal factors line up
  auto lhsSize = lhs.graphPtr->factorPtrs.size();
  auto rhsSize = rhs.graphPtr->factorPtrs.size();
  if(lhsSize!=rhsSize){
    structureOldToNew();
    lhsSize = lhs.graphPtr->factorPtrs.size();
  }

  for(size_t i =0; i < lhsSize; i++){ //We can compare the actual vars, but that might end up very expensive
    if(lhs.graphPtr->factorPtrs[i]->getVars()!=rhs.graphPtr->factorPtrs[i]->getVars()){
      structureOldToNew();
      break;
    }
  }

  //TODO: May we dampen the factors seperately?
  //TODO: may we take the sum of the returns as the new return?

  double dfOld = oldWeight;//lhsSize;
  double dfNew = 0;
  for(size_t i = 0; i<lhsSize; i++){
    dfNew += lhs.graphPtr->factorPtrs[i]->inplaceDampen(rhs.graphPtr->factorPtrs[i], dfOld);
  }

  //TODO: should we calibrate? Or does a calibrateted thing stay a calibrateted thing?
  //rhs.calibrateteWithMessagePassing();

  return dfNew;

} // inplaceProcess
