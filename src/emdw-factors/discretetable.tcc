/*******************************************************************************

          AUTHOR: JA du Preez
          DATE:    March 2015
          PURPOSE: DiscreteTable PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/
// standard headers
#include <random>
#include <limits>

// patrec headers
#include "prlite_stdfun.hpp"  // realEqual

// emdw headers
#include "vecset.hpp"
#include "sortindices.hpp"
#include "combinations.hpp"
#include "oddsandsods.hpp" // almostEqual prodUnderflowProtected

// ****************************************************************************
// *********************** DiscreteTable **************************************
// ****************************************************************************

/*
  Om diskrete faktore vinnig te maak, gebruik hierdie klas 'n
  vloer-waarde waaronder enige waarskynlikhede bloot tot die
  vloerwaarde gelig word. Dit hoort die faktore yler te maak. Maar dit
  kom ongelukkig nie sonder komplikasies nie. o.a. veroorsaak dit dat
  ons nodig het om die domeins van elke toevalsveranderlike te lys.
*/


//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
template <typename AssignType>
void DiscreteTable<AssignType>::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new DiscreteTable_Marginalize<AssignType>() );
  } //if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new DiscreteTable_InplaceNormalize<AssignType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new DiscreteTable_Normalize<AssignType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new DiscreteTable_InplaceAbsorb<AssignType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new DiscreteTable_Absorb<AssignType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new DiscreteTable_InplaceCancel<AssignType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new DiscreteTable_Cancel<AssignType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new DiscreteTable_ObserveAndReduce<AssignType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper =
      uniqptr<InplaceDamper>( new DiscreteTable_InplaceDamping<AssignType>() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new DiscreteTable_Sampler<AssignType>() );
  } // if

} // supplementWithDefaultOperators

// default constructor
template <typename AssignType>
DiscreteTable<AssignType>::DiscreteTable(
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler)
    : marginalizer(theMarginalizer),
      inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper),
      sampler(theSampler){

  supplementWithDefaultOperators();
  } // default constructor


// class specific ctor
template <typename AssignType>
DiscreteTable<AssignType>::DiscreteTable(
  const emdw::RVIds& theVars,
  const std::vector< rcptr< std::vector<AssignType> > >& theDomains,
  double theDefProb,
  const std::map<std::vector<AssignType>, FProb>& theSparseProbs,
  double theMargin,
  double theRelFloor,
  bool presorted,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler)
: vars( theVars.size() ), //
  domains( theDomains.size() ),    //
  defProb(theDefProb),
  sparseProbs(),          //
  margin(theMargin),
  relFloor(theRelFloor),
  marginalizer(theMarginalizer),
  inplaceNormalizer(theInplaceNormalizer),
  normalizer(theNormalizer),
  inplaceAbsorber(theInplaceAbsorber),
  absorber(theAbsorber),
  inplaceCanceller(theInplaceCanceller),
  canceller(theCanceller),
  observeAndReducer(theObserveAndReducer),
  inplaceDamper(theInplaceDamper),
  sampler(theSampler),
  cards( theVars.size() ),
  maxSz(1){

  supplementWithDefaultOperators();
  PRLITE_ASSERT(theVars.size() == theDomains.size(),
         "DiscreteTable cs ctor : " << theVars.size()
         << " variables, but " << theDomains.size()
         << " domains");

  if (presorted or theVars.size() == 1) {
    vars = theVars;
    domains = theDomains;
    auto itr1 = cards.begin();
    auto itr2 = domains.begin();
    while ( itr1 != cards.end() ){
      *itr1 = (*itr2++)->size();
      maxSz *= *itr1++;
    } // while

    sparseProbs = theSparseProbs;
    if (sparseProbs.size() == maxSz) {defProb = 0.0;} // if
    return;
  } // if

  std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );

  // fix the RV order
  vars = extract<unsigned>(theVars, sorted);
  domains = extract< rcptr< std::vector<AssignType> > >(theDomains, sorted);

  // fix the domain dependant stuff
  auto itr1 = cards.begin();
  auto itr2 = domains.begin();
  while ( itr1 != cards.end() ){
    *itr1 = (*itr2++)->size();
    maxSz *= *itr1++;
  } // while

  // fix the RV order in the sparseProbs
  for (auto aProb : theSparseProbs) {

    // This checks that the dimension of theVars and theTable matches.
    PRLITE_ASSERT(theVars.size() == aProb.first.size(),
           "DiscreteTable cs ctor : " << theVars.size()
           << " variables, but " << aProb.first.size()
           << " assignment(s)");

    sparseProbs[extract<AssignType>(aProb.first, sorted)] = aProb.second;

  } // for
  if (sparseProbs.size() == maxSz) {defProb = 0.0;} // if

} // class specific ctor

// destructor
template <typename AssignType>
DiscreteTable<AssignType>::~DiscreteTable() {
} // destructor

//========================= Operators and Conversions =========================

// operator ==
template <typename AssignType>
bool DiscreteTable<AssignType>::operator==(const DiscreteTable<AssignType>& p) const {

  const unsigned decimalPlaces = 10;

  // Test to see whether potentials are over the same set of variables
  if (vars != p.vars) {
    FILEPOS; return false;
  } // if

  // isn't this one already implied by the vars test?
  auto itr1 = domains.begin();
  auto itr2 = p.domains.begin();
  while ( itr1 != domains.end() ) {
    if (**itr1++ != **itr2++) {
      //std::cout << **itr1 << std::endl << **itr2 << std::endl;
      FILEPOS; return false;
    } // if
  } // while

  if (!almostEqual(defProb, p.defProb, decimalPlaces) and sparseProbs.size() < maxSz and p.sparseProbs.size() < maxSz) {
    FILEPOS; return false;
  } // if

  for (auto lVal : sparseProbs) {
    auto rItr = p.sparseProbs.find(lVal.first);
    if ( rItr != p.sparseProbs.end() ) { // check against rhs sparse
      if  ( !almostEqual(lVal.second.prob, rItr->second.prob, decimalPlaces) ) {
        FILEPOS; return false;
      } // if
    } // if
    else { // check against rhs default
      if ( !almostEqual(lVal.second.prob, p.defProb, decimalPlaces) ) {
        FILEPOS; return false;
      } // if
    } // else
  } // for

  for (auto rVal : p.sparseProbs) {
    auto lItr = sparseProbs.find(rVal.first);
    if ( lItr != sparseProbs.end() ) { // check against rhs sparse
      if ( !almostEqual(rVal.second.prob, lItr->second.prob, decimalPlaces) ) {
        FILEPOS; return false;
      } // if
    } // if
    else { // check against rhs default
      if ( !almostEqual(rVal.second.prob, defProb, decimalPlaces) ) {
        FILEPOS; return false;
      } // if
    } // else
  } // for

  return true;
} // operator ==

// operator !=
template <typename AssignType>
bool DiscreteTable<AssignType>::operator!=(const DiscreteTable<AssignType>& p) const {
  return !operator==(p);
} // operator !=

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// classSpecificConfigure
template <typename AssignType>
unsigned DiscreteTable<AssignType>::classSpecificConfigure(
  const emdw::RVIds& theVars,
  const std::vector< rcptr< std::vector<AssignType> > >& theDomains,
  double theDefProb,
  const std::map<std::vector<AssignType>, FProb>& theSparseProbs,
  double theMargin,
  double theRelFloor,
  bool presorted,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper,
  rcptr<Sampler> theSampler) {

  this->~DiscreteTable<AssignType>(); // Destroy existing

  new(this) DiscreteTable<AssignType>(
    theVars,
    theDomains,
    theDefProb,
    theSparseProbs,
    theMargin,
    theRelFloor,
    presorted,
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

  return true;
} // classSpecificConfigure

// txtWrite
template <typename AssignType>
std::ostream& DiscreteTable<AssignType>::txtWrite(std::ostream& file) const {
  file << "DiscreteTable_V0\n";
  file << "Vars: " << vars;
  file << "\nDomains:\n";
  for (auto aDom : domains) {
    file << *aDom << std::endl;
  } // for
  file << "\nDefProb: " << defProb
       << " Margin: " << margin
       << " RelFloor: " << relFloor
       << "\nNoOfSparseProbs: " << sparseProbs.size() << std::endl;
  for (auto spec : sparseProbs) {
    file << spec.first << " " << spec.second.prob << std::endl;
  } // for
  file << std::endl;
  file << marginalizer->isA() << std::endl;
  file << inplaceNormalizer->isA() << std::endl;
  file << normalizer->isA() << std::endl;
  file << inplaceAbsorber->isA() << std::endl;
  file << absorber->isA() << std::endl;
  file << inplaceCanceller->isA() << std::endl;
  file << canceller->isA() << std::endl;
  file << observeAndReducer->isA() << std::endl;
  file << inplaceDamper->isA() << std::endl;
  file << sampler->isA() << std::endl;
  return file;
} // txtWrite

// txtRead
template <typename AssignType>
std::istream& DiscreteTable<AssignType>::txtRead(std::istream& file) {

  std::string dummy;

  emdw::RVIds theVars;
  double theDefProb;
  size_t nSparse;
  double theMargin;
  double theRelFloor;

  file >> dummy; // version tag
  file >> dummy >> theVars;

  std::vector< rcptr< std::vector<AssignType> > > theDomains( theVars.size() );
  file >> dummy; // to bypass the "Domains:" text (added by Rebecca 09112017)
  for (auto& aDom : theDomains) {
    aDom = rcptr< std::vector<AssignType> >(new std::vector<AssignType>() );
    file >> *aDom;
  } // for

  file >> dummy >> theDefProb
       >> dummy >> theMargin
       >> dummy >> theRelFloor
       >> dummy >> nSparse;

  std::map<std::vector<AssignType>, FProb> theSparseProbs;
  while (nSparse--) {
    std::vector<AssignType> theVals; double theProb;
    PRLITE_CHECK(file >> theVals, IOError, "Error reading parameter value.");
    PRLITE_CHECK(file >> theProb, IOError, "Error reading parameter count.");
    theSparseProbs[theVals] = theProb;
  } // while

  std::cout << "DiscreteTable::txtRead : FIXME factor operators will revert to defaults\n";
  classSpecificConfigure(theVars, theDomains, theDefProb, theSparseProbs, theMargin, theRelFloor, true);

  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
template <typename AssignType>
emdw::RVVals DiscreteTable<AssignType>::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
template <typename AssignType>
uniqptr<Factor> DiscreteTable<AssignType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
template <typename AssignType>
void DiscreteTable<AssignType>::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> DiscreteTable<AssignType>::absorb(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
template <typename AssignType>
void DiscreteTable<AssignType>::inplaceAbsorb(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> DiscreteTable<AssignType>::cancel(const Factor* rhsPtr,
                                                  FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
template <typename AssignType>
void DiscreteTable<AssignType>::inplaceCancel(const Factor* rhsPtr,
                                              FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> DiscreteTable<AssignType>::marginalize(const emdw::RVIds& variablesToKeep,
                                                       bool presorted,
                                                       const Factor* peekAheadPtr,
                                                       FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep, presorted, peekAheadPtr) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep, presorted, peekAheadPtr) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> DiscreteTable<AssignType>::observeAndReduce(const emdw::RVIds& variables,
                                                            const emdw::RVVals& assignedVals,
                                                            bool presorted,
                                                            FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals, presorted) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals, presorted) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
template <typename AssignType>
double DiscreteTable<AssignType>::inplaceDampen(const Factor* rhsPtr,
                                                double df,
                                                FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

// copy
template <typename AssignType>
DiscreteTable<AssignType>* DiscreteTable<AssignType>::copy(
  const emdw::RVIds& newVars, bool presorted) const {
  if ( newVars.size() ) {
    return new DiscreteTable<AssignType>(
      newVars,
      domains,
      defProb,
      sparseProbs,
      margin,
      relFloor,
      presorted,
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
  }
  return new DiscreteTable<AssignType>(*this);
} // copy

// vacuousCopy
template <typename AssignType>
DiscreteTable<AssignType>*
DiscreteTable<AssignType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  if ( !selVars.size() ) { // use all variables
    return new DiscreteTable<AssignType>(
      vars, domains, 1.0/maxSz,
      std::map<std::vector<AssignType>, FProb>(),
      margin, relFloor, true,
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
  } // else

  else { // use a subset of variables

    // sort the selected variables if necessary
    const emdw::RVIds* sortedVarsPtr = &selVars;
    emdw::RVIds tmpVars; // need it here to not go out of scope too soon
    std::vector<size_t> sorted; // need it here to not go out of scope too soon
    if (!presorted) {
      sorted = sortIndices( selVars, std::less<unsigned>() );
      tmpVars = extract<unsigned>(selVars,sorted);
      sortedVarsPtr = &tmpVars;
    } // if
    const emdw::RVIds& sortedVars(*sortedVarsPtr);

    // find the subset
    std::vector<size_t> l2i, r2i;
    sortedIntersection(vars, sortedVars, l2i, r2i);

    // check that there are no stragglers
    PRLITE_ASSERT(r2i.size() == selVars.size(), "DiscreteTable::vacuousCopy : Some selected variables are not from the original factor");

    // gather the domains for the subset
    std::vector< rcptr< std::vector<AssignType> > > selDoms(
      extract<rcptr< std::vector<AssignType> >>(domains, l2i) );

    // and build it
    return new DiscreteTable<AssignType>(
      sortedVars, selDoms, 1.0,
      std::map<std::vector<AssignType>, FProb>(),
      margin, relFloor, true,
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
  } // else

} // vacuousCopy

// Boiler plate code
// isEqual
template <typename AssignType>
bool DiscreteTable<AssignType>::isEqual(const Factor* rhsPtr) const {
  const DiscreteTable<AssignType>* dwnPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

// KL(p||q) = <log(p)>_p - <log(q)>_p
//          = <log(p)>_p + <log(N)>_p
//          = -H_p + log(N).
//
// I.e. KL and -H only differ from each other by a constant
// log(N). However, KL is always positive and -H is always
// negative.
//
// Interestingly, using -H instead of KL appears (???) to be slightly
// faster. Because all the loopy distances are positive (being
// symmetrical KL), using neg entropy results in first completing the
// (positive valued) effects of the chosen factor before popping an
// untouched 'vacuous' version (with its negative values) from the
// queue. Somehow this seems to be better.
template <typename AssignType>
double
DiscreteTable<AssignType>::distanceFromVacuous() const {
  return -entropy() + log2(maxSz);
} // distanceFromVacuous

// noOfVars
template <typename AssignType>
unsigned DiscreteTable<AssignType>::noOfVars() const {
  return vars.size();
} // noOfVars

// getVars
template <typename AssignType>
emdw::RVIds DiscreteTable<AssignType>::getVars() const {
  return vars;
} // getVars

// getVar
template <typename AssignType>
emdw::RVIdType DiscreteTable<AssignType>::getVar(unsigned varNo) const {
  return vars[varNo];
} // getVar

//====================== Other Virtual Member Functions =======================

// distance: one sided P||Q KL. for some reason this works better than the symmetric version
template <typename AssignType>
double DiscreteTable<AssignType>::distance(const Factor* rhsPtr) const {

  const std::map<std::vector<AssignType>, FProb>& mapP(sparseProbs);         // convenient alias
  const double defP(defProb);                              // convenient alias
  const double logDefP = defP > 0.0 ? log(defP) : -1000.0;
  const double defPLogDefP = defP > 0.0 ? defP*logDefP : 0.0;

  // Ensure rhs is also a DiscreteTable
  const DiscreteTable<AssignType>* dwnPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsPtr);
  PRLITE_ASSERT(dwnPtr, "DiscreteTable::distance : rhs is not a DiscreteTable");
  const std::map<std::vector<AssignType>, FProb>& mapQ(dwnPtr->sparseProbs); // convenient alias
  const double defQ(dwnPtr->defProb);                      // convenient alias
  const double logDefQ = defQ > 0.0 ? log(defQ) : -1000.0;
  //const double defQLogDefQ = defQ > 0.0 ? defQ*logDefQ : 0.0;

  const double defPLogDefQ = defP > 0.0 ? defP*logDefQ : 0.0;
  //const double defQLogDefP = defQ > 0.0 ? defQ*logDefP : 0.0;

  // Ensures the number of parameters of the two DiscreteTable potential functions are equal
  PRLITE_ASSERT(maxSz == dwnPtr->maxSz,
         "DiscreteTable::distance The maxSz of the two distributions differ");

  // All the variables we'll need for calculating the distance

  // to calc the number of defaults occurring in neither of the maps
  size_t nDefsP = maxSz-mapP.size();
  size_t nDefsQ = maxSz-mapQ.size();
  double dist = 0.0; // dist1 = 0.0, dist2 = 0.0;

  auto iP = mapP.begin();
  auto iQ = mapQ.begin();
  while ( iP != mapP.end() && iQ != mapQ.end() ){

    if (iP->first < iQ->first) {      // only in mapP
      nDefsQ--;

      double P(iP->second.prob);
      double logP = P > 0.0 ? log(P) : -1000.0;

      if (P > 0.0) {dist += P*(logP -logDefQ);} // if
      //if (defQ > 0.0) {dist += defQLogDefQ -defQ*logP;} // if

      iP++;
    } // if

    else if (iP->first > iQ->first) { // only in mapQ
      nDefsP--;

      double Q(iQ->second.prob);
      double logQ = Q > 0.0 ? log(Q) : -1000.0;

      //if (Q > 0.0) {dist += Q*(logQ -logDefP);} // if
      if (defP > 0.0) {dist += defPLogDefP -defP*logQ;} // if

      iQ++;
    } // else if

    else {                            //common to both maps

      double P(iP->second.prob);
      double logP = P > 0.0 ? log(P) : -1000.0;
      double Q(iQ->second.prob);
      double logQ = Q > 0.0 ? log(Q) : -1000.0;

      if (P > 0.0) {dist += P*(logP -logQ);} // if
      //if (Q > 0.0) {dist += Q*(logQ -logP);} // if

      iP++; iQ++;
    } // else

  } // while

  while ( iP != mapP.end() ){         // some extra spec's in mapP
    nDefsQ--;

    double P(iP->second.prob);
    double logP = P > 0.0 ? log(P) : -1000.0;

    if (P > 0.0) {dist += P*(logP -logDefQ);} // if
    //if (defQ > 0.0) {dist += defQLogDefQ -defQ*logP;} // if

    iP++;
  } // while

  while ( iQ != mapQ.end() ){         // some extra spec's in mapQ
    nDefsP--;

    double Q(iQ->second.prob);
    double logQ = Q > 0.0 ? log(Q) : -1000.0;

    //if (Q > 0.0) {dist += Q*(logQ -logDefP);} // if
    if (defP > 0.0) {dist += defPLogDefP -defP*logQ;} // if

    iQ++;
  } // while

  // and now all the cases simultaneously absent from both maps
  PRLITE_ASSERT(
    nDefsP == nDefsQ,
    "DiscreteTable::distance The number of cases in neither map should match, no matter from which side you look at it.");

  if (nDefsP) {
    dist += nDefsP*(defPLogDefP - defPLogDefQ); // + defQLogDefQ - defQLogDefP);
  } // if

      return dist; // /2.0;
} // distance

template <typename AssignType>
double DiscreteTable<AssignType>::potentialAt(
  const emdw::RVIds& theVars,
  const emdw::RVVals& anyVals,
  bool presorted) const{

  // A typedef for the actual internal assignment vector type
  if (theVars.size() == 1) {presorted = true;}

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);
  PRLITE_ASSERT(sortedVars == vars,
         "DiscreteTable::potentialAt : The variables do not match: "
         << sortedVars << " vs " << vars);

  std::vector<AssignType> observedVals( anyVals.size() );
  for (unsigned idx = 0; idx < anyVals.size(); idx++){
    if (presorted) {
      observedVals[idx] = static_cast<AssignType>(anyVals[idx]);
    } //if
    else {
      observedVals[idx] = static_cast<AssignType>(anyVals[ sorted[idx] ]);
    } // else
  } // for
  auto itr = sparseProbs.find(observedVals);
  return itr == sparseProbs.end() ? defProb : itr->second.prob;
} // potentialAt

template <typename AssignType>
void
DiscreteTable<AssignType>::setEntry(
  const emdw::RVIds& theVars,
  const emdw::RVVals& theirVals,
  double theProb,
  bool presorted ){

  if (theVars.size() == 1) {presorted = true;} // if

  if (presorted) {
    PRLITE_ASSERT(theVars == vars, "The variables do not match");
    std::vector<AssignType> assignedVals( theirVals.size() );
    for (unsigned idx = 0; idx < theirVals.size(); idx++) {
      assignedVals[idx] = static_cast<AssignType>(theirVals[idx]);
    } // for
    sparseProbs[assignedVals] = theProb;
    return;
  } // if

  std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );
  PRLITE_ASSERT(theVars.size() == theirVals.size(),
         theVars.size() << " variables, but " <<
         theirVals.size() << " assignment(s)");
  PRLITE_ASSERT(extract<unsigned>(theVars, sorted) == vars, "The variables do not match");
  sparseProbs[extract<AssignType>(theirVals, sorted)] = theProb;

} // setEntry

template <typename AssignType>
double DiscreteTable<AssignType>::sumZ() const {
  double sum = ( maxSz-sparseProbs.size() ) * defProb;
  for (const auto&  val : sparseProbs) {
    sum += val.second.prob;
  } // for
  return sum;
} // sumZ

//======================= Non-virtual Member Functions ========================

/// inplace absorb a Cat1EXPF. The code is mostly scavenged from the
/// corresponding DiscreteTable operator

template <typename AssignType>
void DiscreteTable<AssignType>::inplaceAbsorbCat1EXPF(const rcptr< Cat1EXPF<AssignType> >& c1Ptr) {

  const Cat1EXPF<AssignType>& rhs(*c1Ptr);
  const std::map<AssignType, double>& mapR(rhs.sparseProbs);
  const double defR = rhs.defProb;

  std::map< std::vector<AssignType>, FProb>& mapL(sparseProbs);
  double& defL = defProb;

  // where in vars can c1Ptr->varId be found
  size_t l2i = -1;
  // where in vars can all RVs except c1Ptr->varId be found
  std::vector<size_t> lni; lni.reserve(vars.size()-1);
  for (size_t idx = 0; idx < vars.size(); idx++) {
    if (vars[idx] == c1Ptr->varId) {l2i = idx;} // if
    else {lni.push_back(idx);} // else
  } // for
  PRLITE_ASSERT(l2i < vars.size(), "Variable not in cluster scope: " << c1Ptr->varId);

  // Case 3: Indien die lhs defVal > 0.0 itereer ons oor die *rhs*
  // sparse waardes. Verwerking minder as O( N*log2(N) ). Die hoop is
  // dat die vloerwaarde N en M klein genoeg kan maak om dit die
  // moeite werd te maak. Vir rhs Cat1EXPF sal ons ook soortgelyk moet doen.

  if (defProb > 0.0) { // we have a non-zero default on the lhs probs

    //std::cout << "Case 3*\n";

    // Elkeen van daardie veranderlike assignments word dan (een vir
    // een) uitgebrei na al die moontlike versoenbare lhs veranderlike
    // assignments. Daarna soek ons in die lhs vir elk van hierdie
    // assignments. Indien daar, vermenigvuldig die lhs prob met die
    // rhs prob en stel die flag in die lhs. Indien nie daar nie,
    // vermenigvuldig die rhs prob met die lhs defVal en insert 'n
    // nuwe inskrywing in die lhs sparse values, en stel die flag.

    // initialize stuff outside the iterator loop
    std::vector<size_t> newSizes(lni.size() );
    for (size_t k = 0; k < lni.size(); k++) {
      newSizes[k] = domains[lni[k]]->size();
    } // for
    std::vector<size_t> initIndices = initTuples<size_t>(newSizes);
    std::vector<size_t> indices( initIndices.size() );
    std::vector<AssignType> jointAssign(vars.size() );
    double newDef = defL*defR;

    auto iR = mapR.begin();
    while ( iR != mapR.end() ) {

      // first copy the values of the current rhs assignment
      jointAssign[ l2i ] = iR->first;

      // now extend it with the missing lhs variable assignments
      indices = initIndices;
      do {

        // copy the new variable assignments into place
        for (unsigned k = 0; k < lni.size(); k++) {
          jointAssign[ lni[k] ] = (*domains[ lni[k] ])[ indices[k] ];
        } // for

        auto result = mapL.emplace(jointAssign, defL); // leaves existing map records unchanged
        result.first->second.prob *= iR->second;
        if (std::abs(result.first->second.prob - newDef) <= margin) {
          mapL.erase(result.first);
        } // if
        else {
          result.first->second.flag = 1;
        } // else

      } while( nextCombination(newSizes, indices) );

      iR++;
    } // while

    // Daarna itereer ons deur al die lhs sparse entries. flag == 1
    // beteken dat reeds met rhs prob vermenigvuldig is. Met flag == 0
    // moet ons nog vermenigvuldig met die defR.

    for (auto& aProb : mapL) {
      if (!aProb.second.flag) { // cases not visible from mapR
        aProb.second.prob *= defR;
      } // if
      aProb.second.flag = 0;
    }  // for

    // en heel laaste moet ons die lhs default ook opdateer
    defL = newDef;

    return;
  } // if

  // Case 4: Gevalle waar die lhs defVal == 0 itereer ons oor die *lhs* sparse waardes, onttrek
  // die tersake veranderlikes in die rhs en doen dan 'n find op die
  // rhs. As ons dit vind word die lhs prob met ooreenstemmende rhs
  // prob vermenigvuldig, anders word dit met rhs defVal
  // vermenigvuldig. O( M*log2(N) ) QED

  //std::cout << "Case 4*\n";

  auto iL = mapL.begin();
  while ( iL != mapL.end() ) {
    auto iR = mapR.find( iL->first[ l2i ]);
    if ( iR != mapR.end() ) {
      iL->second.prob *= iR->second;
    } // if
    else {iL->second.prob *= defR;} // else
    iL++;
  } // while

} // inplaceAbsorbCat1EXPF

template <typename AssignType> //!!!
Cat1EXPF<AssignType>* DiscreteTable<AssignType>::marginalizeToCat1EXPF(emdw::RVIdType varId) const {

  size_t idxOfVar = -1;
  for (size_t idx = 0; idx < vars.size(); idx++) {
    if (vars[idx] == varId) {
      idxOfVar = idx;
      break;
    } // if
  } // for
  PRLITE_ASSERT(idxOfVar < vars.size(), "Variable not in cluster scope: " << varId);

  std::map<AssignType, double> theSparseProbs;

  if (defProb > 0) { // more stuff to consider

    // each value in the domain of varId must appear in this many contexts
    size_t noOfVariants = maxSz / domains[idxOfVar]->size();

    std::map<AssignType, size_t> occurrences;
    for (auto el : sparseProbs) {
      theSparseProbs[ el.first[idxOfVar] ] += el.second.prob;
      occurrences[ el.first[idxOfVar] ]++;
    } // for

    for (auto el : occurrences) {
      theSparseProbs[el.first] += defProb*(noOfVariants - el.second);
    } // if

    double theDefProb = 0.0;
    if ( occurrences.size() < domains[idxOfVar]->size() ) { // some values absent from sparseProbs
      theDefProb = defProb*noOfVariants;
    } // if

    return new Cat1EXPF<AssignType>(varId, domains[idxOfVar], theDefProb, theSparseProbs);
  } // if

  else { // lucky us, defProbs == 0 is simple
    for (auto el : sparseProbs) {
      theSparseProbs[ el.first[idxOfVar] ] += el.second.prob;
    } // for

    return new Cat1EXPF<AssignType>(varId, domains[idxOfVar], 0.0, theSparseProbs);
  } // else

} // marginalizeToCat1EXPF

template <typename AssignType>
void DiscreteTable<AssignType>::makeDefaultZero() const {

  if (defProb > 0.0) {
    //typedef std::vector<AssignType> std::vector<AssignType>;

    // iterate over all possible assignments for the newVars
    std::vector<size_t> indices = initTuples<size_t>(cards);
    do {

      std::vector<AssignType> assignVals( indices.size() );
      for (size_t k = 0; k < indices.size(); k++) {
        assignVals[k] = (*domains[k])[ indices[k] ];
      } // for

      if ( sparseProbs.find(assignVals) == sparseProbs.end() ){
        sparseProbs[assignVals] = defProb;
      } // if

    } while( nextCombination(cards, indices) );

    defProb = 0.0;
  } // if

} // makeDefaultZero

template <typename AssignType>
void DiscreteTable<AssignType>::addNoise(double upperLim){
  static std::random_device randDev{};
  static std::default_random_engine randEng{randDev()};
  static std::uniform_real_distribution<> randUnif(0.0, upperLim);
  for (auto& el : sparseProbs) {
    el.second.prob += randUnif(randEng);
  } // for
} // addNoise

template <typename AssignType>
void DiscreteTable<AssignType>::squashSmallProbs(double thresh){
  if (defProb < thresh) {defProb = 0.0;} // if
  for (auto& el : sparseProbs) {
    if (el.second.prob < thresh) {el.second.prob = 0.0;} // if
  } // for
} // squashSmallProbs

template <typename AssignType>
void DiscreteTable<AssignType>::nonZeroToUnity(){
  if (defProb > 0.0) {defProb = 1.0;} // if
  for (auto& el : sparseProbs) {
    if (el.second.prob > 0.0) {el.second.prob = 1.0;} // if
  } // for
} // nonZeroToUnity

/*
 * Calculates the entropy for the distribution, even if unnormalized.
 *
 * $h = \sum_i p_i \log(p_i)$
 * If distribution is not normalized, this changes to:
 * $h = \sum_i (p_i/k) \log(p_i/k)$ with $k = \sum_i p_i$.
 * This simplifies to:
 * $h = \log(\sum_i p_i) - \frac{\sum_i p_i \log(p_i)}{\sum_i p_i}$
 */

template <typename AssignType>
double DiscreteTable<AssignType>::entropy() const {
  double entr = 0.0;
  double sump = 0.0;

  if (defProb > 0.0) {
    size_t noOfDefs = maxSz - sparseProbs.size();
    sump = noOfDefs*defProb;
    entr = sump*log2(defProb);
  } // if

  for (auto el : sparseProbs) {
    if (el.second.prob > 0.0) {
      sump += el.second.prob;
      entr += el.second.prob*log2(el.second.prob);
    } // if
  } // for

  return log2(sump) - entr/sump;
} // entropy

template <typename AssignType>
double DiscreteTable<AssignType>::getProb(const std::vector<AssignType>& theAss) const {
  auto itr = sparseProbs.find(theAss);
  if ( itr == sparseProbs.end() ) {return defProb;} // if
  else {return itr->second.prob;} // else
} // getProb

template <typename AssignType>
double DiscreteTable<AssignType>::getDefaultProb() const {
  return defProb;
} // getDefaultProb

template <typename AssignType>
const std::map<std::vector<AssignType>, FProb>& DiscreteTable<AssignType>::getSparseProbs() const {
  return sparseProbs;
} // getSparseProbs

template <typename AssignType>
unsigned DiscreteTable<AssignType>::getCardinality(unsigned varNo) const {
  return cards[varNo];
} // getCardinality

template <typename AssignType>
void DiscreteTable<AssignType>::copyExtraSettingsFrom(const DiscreteTable<AssignType>& d) {
  margin = d.margin;
  relFloor = d.relFloor;
  marginalizer = d.marginalizer;
  inplaceNormalizer = d.inplaceNormalizer;
  normalizer = d.normalizer;
  inplaceAbsorber = d.inplaceAbsorber;
  absorber = d.absorber;
  inplaceCanceller = d.inplaceCanceller;
  canceller = d.canceller;
  observeAndReducer = d.observeAndReducer;
  inplaceDamper = d.inplaceDamper;
  sampler = d.sampler;
} // copyOperatorsFrom

//======================= Factor operators ====================================

// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample   **************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Sampler<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Sampler<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Samples a random vector from the distribution
//
// WHAMMY: For an efficient sampling process this sampler will
// propagate non-zero default into the probs table, in effect making
// it dense. This can make the probs table prohibitively large in
// tasks where the majority of probs are zero by design (such as in
// sudoku).  The memory efficient alternative will unfortunately be
// much more cpu intensive since it will have to dynamically recreate
// all the assignments corresponding to the default probs.

template <typename AssignType>
emdw::RVVals DiscreteTable_Sampler<AssignType>::sample(const DiscreteTable<AssignType>* lhsPtr) {

  if (lhsPtr->defProb > 0.0) {lhsPtr->makeDefaultZero();} // if

  //static std::random_device randDev{};
  //static std::default_random_engine randEng{randDev()};
  static std::uniform_real_distribution<> randUnif(0.0, 1.0);

  emdw::RVVals ret;
  double accum = 0.0;
  double randVal =   randUnif( emdw::randomEngine() );
  auto iL = lhsPtr->sparseProbs.begin();
  while ( iL != lhsPtr->sparseProbs.end() ){
    accum += iL->second.prob;
    if (accum >= randVal) {
      ret.assign( iL->first.begin(), iL->first.end() );
      break;
    } // if
    iL++;
  } // while

  return ret;
} // sample

// ***************************************************************************
// *************************** Family 1 **************************************
// *************************** Normalize *************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Normalizes the sum to one, applies a relFloor, make it more sparse where possible
template <typename AssignType>
void DiscreteTable_InplaceNormalize<AssignType>::inplaceProcess(DiscreteTable<AssignType>* lhsPtr) {

  double relFloor = lhsPtr->relFloor;
  double margin = lhsPtr->margin;
  double& defVal = lhsPtr->defProb;

  double theSum = 0.0;
  double theMax = defVal;
  auto itr = lhsPtr->sparseProbs.begin();
  while ( itr != lhsPtr->sparseProbs.end() ){
    if (std::abs(itr->second.prob-defVal)  <= margin) {
      itr = lhsPtr->sparseProbs.erase(itr);
    } // if
    else{
      theSum += itr->second.prob;
      if (itr->second.prob > theMax) {
        theMax = itr->second.prob;
      } // if
      itr++;
    } // else
  } // while

  double thresh = theMax*relFloor;
  if (defVal < thresh) {defVal = thresh;} // if
  size_t defSz = lhsPtr->maxSz-lhsPtr->sparseProbs.size();
  if (!defSz) {defVal = thresh;} // if : we are free to change the default
  theSum += defSz*defVal;
  if ( prlite::realEqual(theSum, 0.0) ) {
    std::string msg = "DiscreteTable_InplaceNormalize:: Empty factor with vars size "+std::to_string(lhsPtr->noOfVars())+" is fracturing the PGM";
    throw( msg.c_str() );   // just to still throw a char* exception.
  } // if
  defVal /= theSum;

  double err = 0.0;
  for (auto& aProb : lhsPtr->sparseProbs) {
    if (aProb.second.prob < thresh) {
      err += thresh - aProb.second.prob;
      aProb.second.prob = thresh;
    } // if
    aProb.second.prob /= theSum;
  } // for

  if (err > 0.0) { // we divided by a too small number
    double correction = theSum/(theSum+err);
    defVal *= correction;
    for (auto& aProb : lhsPtr->sparseProbs) { // damn, once again
      aProb.second.prob *= correction;
    } //
  } // if

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Normalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Normalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* DiscreteTable_Normalize<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr) {
  // make a copy of the factor
  DiscreteTable<AssignType>* fPtr = new DiscreteTable<AssignType>(*lhsPtr);

  // inplace normalize it
  DiscreteTable_InplaceNormalize<AssignType> ipNorm;
  ipNorm.inplaceProcess(fPtr);
  // try{
  //   ipNorm.inplaceProcess(fPtr);
  // } // try
  // catch (const char* s) {
  //   std::cout << __FILE__ << __LINE__
  //             << " call to 'inplaceProcess' failed" << std::endl;
  //   std::cout << s << std::endl;
  //   throw s;
  // } // catch

  // and return the result
  return fPtr;
} // process

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceMaxNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceMaxNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// MaxNormalizes the sum to one, applies a relFloor, make it more sparse where possible
template <typename AssignType>
void DiscreteTable_InplaceMaxNormalize<AssignType>::inplaceProcess(DiscreteTable<AssignType>* lhsPtr) {

  double relFloor = lhsPtr->relFloor;
  double margin = lhsPtr->margin;
  double& defVal = lhsPtr->defProb;

  double theMax = defVal;
  auto itr = lhsPtr->sparseProbs.begin();
  while ( itr != lhsPtr->sparseProbs.end() ){
    if (std::abs(itr->second.prob-defVal)  <= margin) {
      itr = lhsPtr->sparseProbs.erase(itr);
    } // if
    else{
      if (itr->second.prob > theMax) {
        theMax = itr->second.prob;
      } // if
      itr++;
    } // else
  } // while

  double thresh = theMax*relFloor;
  if (defVal < thresh) {defVal = thresh;} // if
  size_t defSz = lhsPtr->maxSz-lhsPtr->sparseProbs.size();
  if (!defSz) {defVal = thresh;} // if : we are free to change the default
  if ( prlite::realEqual(theMax, 0.0) ) {
    throw("DiscreteTable_InplaceNormalize:: Empty factor is fracturing the PGM");
  } // if
  defVal /= theMax;

  for (auto& aProb : lhsPtr->sparseProbs) {
    if (aProb.second.prob < thresh) {
      aProb.second.prob = thresh;
    } // if
    aProb.second.prob /= theMax;
  } // for

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_MaxNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_MaxNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* DiscreteTable_MaxNormalize<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr) {
  // make a copy of the factor
  DiscreteTable<AssignType>* fPtr = new DiscreteTable<AssignType>(*lhsPtr);

  // inplace normalize it
  DiscreteTable_InplaceMaxNormalize<AssignType> ipNorm;
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

// ***************************************************************************
// *************************** Family 2 **************************************
// *************************** Absorb ****************************************
// ***************************************************************************

/**
 * @param oldMap The map that needs to be extended.
 *
 * @param oldVars The variable ids that were used in oldMap
 *
 *
 * @param i2l Union indices for the old vars that are shared with new vars.
 *
 * @param ni2l Union indices for the unique old vars (i.e. not shared with new vars).
 *
 * @param newVars The sorted new/extra ids that oldMap needs to be extended with.
 *
 * @param ni2r Union indices for the unique new vars (i.e. not shared with old vars).
 *
 * @param newDoms The domains of newVars.
 */
template <typename AssignType>
std::map< std::vector<AssignType>, FProb>
extendSparseTable(
  const std::map< std::vector<AssignType>, FProb>& lMap,
  const emdw::RVIds& lVars,
  const emdw::RVIds& rVars,
  const std::vector< rcptr< std::vector<AssignType> > >& rDoms,
  std::vector<size_t>& l2i,
  std::vector<size_t>& lni,
  std::vector<size_t>& r2i,   // seems this one is not being used
  std::vector<size_t>& rni,
  std::vector<size_t>& i2l,
  std::vector<size_t>& ni2l,
  std::vector<size_t>& ni2r){

  //typedef std::vector<AssignType> std::vector<AssignType>;

  // std::cout << "lVars: " << lVars << std::endl;
  // std::cout << "rVars: " << rVars << std::endl;
  // std::cout << "rDoms:\n";
  // for (auto el : rDoms) {std::cout << *el << std::endl;} // for
  // std::cout << "l2i: " << l2i << std::endl;
  // std::cout << "lni: " << lni << std::endl;
  // std::cout << "r2i: " << r2i << std::endl;
  // std::cout << "rni: " << rni << std::endl;
  // std::cout << "i2l: " << i2l << std::endl;
  // std::cout << "ni2l: " << ni2l << std::endl;
  // std::cout << "ni2r: " << ni2r << std::endl;

  // array of new domain sizes
  std::vector<size_t> newSizes(rni.size() );
  for (size_t k = 0; k < rni.size(); k++) {
    newSizes[k] = rDoms[rni[k]]->size();
  } // for
  //std::cout << "newSizes: " << newSizes << std::endl;

  std::map< std::vector<AssignType>, FProb> newMap;
  std::vector<AssignType> jointAssign( i2l.size() + ni2l.size() + ni2r.size() );

  // iterate over all possible assignments for the newVars
  std::vector<size_t> indices = initTuples<size_t>(newSizes);
  do {

    // copy the new variable assignments into place
    for (unsigned k = 0; k < ni2r.size(); k++) {
      jointAssign[ ni2r[k] ] = (*rDoms[ rni[k] ])[ indices[k] ];
    } // for
    //std::cout << __LINE__ << ": " << jointAssign << std::endl;

    // iterate over all the existing oldMap records
    for (auto el : lMap) {

      // copy the old common/mutual variable assignments into place
      for (unsigned k = 0; k < i2l.size(); k++) {
        jointAssign[ i2l[k] ] = el.first[ l2i[k] ];
      } // for
      //std::cout << __LINE__ << ": " << jointAssign << std::endl;

      // copy the old unique variable assignments into place
      for (unsigned k = 0; k < ni2l.size(); k++) {
        jointAssign[ ni2l[k] ] = el.first[ lni[k] ];
      } // for
      //std::cout << __LINE__ << ": " << jointAssign << std::endl;

      // and associate the newly formed joint and its probability
      newMap[jointAssign] = el.second;

    } // for

  } while( nextCombination(newSizes, indices) );

  return newMap;
} // extendSparseTable

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceAbsorb<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceAbsorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
  Hier is twee faktore ter sprake, ons verwys na hulle as lhs en
  rhs. def is die default waarde, sparse is die detail waardes in die
  map. Groottes van die sparse gedeeltes dui ons aan as M en N. Ons
  gaan slegs gevalle waar rhs 'n subversameling van lhs is,
  hanteer. Dis al wat ons werklik effektief inplace kan doen. Ander
  gevalle moet in DiscreteTable_Absorb hanteer word.
*/

// NOTE: We should consider two further optimizations:

// 1: with relFloor > 0 cases: We can keep track of the running
// product max and remove / not insert cases into the map that, based
// on the current running max, we already know will fall below the
// threshold value.

// 2: Similarly we can also keep the cases where the product is within
// a "margin" margin of the default value, out of the sparseProbs map
// (instead of putting it in there and then rely on normalize to take
// it out again.

template <typename AssignType>
void
DiscreteTable_InplaceAbsorb<AssignType>::inplaceAbsorb(
  DiscreteTable<AssignType>& lhs,
  const DiscreteTable<AssignType>& rhs,
  const std::vector<size_t>& l2i,
  const std::vector<size_t>& lni,
  const std::vector<size_t>& i2l,
  const std::vector<size_t>& ni2l) {

  //typedef std::vector<AssignType> std::vector<AssignType>;

  const std::map< std::vector<AssignType>, FProb>& mapR(rhs.sparseProbs);
  const double defR = rhs.defProb;
  std::map< std::vector<AssignType>, FProb>& mapL(lhs.sparseProbs);
  double& defL = lhs.defProb;

  // Case 1: with a uniform rhs we can simply return the lhs.
  if ( !mapR.size() ) {
    //std::cout << "Case 1*\n";

    // one can also leave the following block out altogether, in that
    // case the result may be unnormalized
    defL *= defR;
    auto iL = mapL.begin();
    while ( iL != mapL.end() ){
      iL->second.prob *= defR;
      iL++;
    } // while
    return;
  } // if

  // Case 2: If we have exactly the same variables on both sides, the
  // double ladder trick should be fastest. Rather useful since
  // this happens a lot in the LBU algo.

  if (lhs.vars == rhs.vars){

    //std::cout << "Case 2*\n";

    // now multiply by simultaneously climbing two ladders

    // to calc the number of defaults occurring in neither of the maps
    size_t nDefsL = lhs.maxSz-mapL.size();
    size_t nDefsR = rhs.maxSz;
    double newDef = defL*defR;

    for (auto k : lni) {nDefsR *= lhs.domains[k]->size();} // for
    nDefsR -= mapR.size();

    auto iL = mapL.begin();
    auto iR = mapR.begin();
    while ( iL != mapL.end() && iR != mapR.end() ){

      if (iL->first < iR->first) {      // only in mapL
        nDefsR--;
        iL->second.prob *= defR;
        iL++;
      } // if

      else if (iL->first > iR->first) { // only in mapR
        nDefsL--;
        double prob = defL * iR->second.prob;
        if (std::abs(prob - newDef) > lhs.margin) {
          mapL.emplace_hint(iL, iR->first, prob);
        } // if
        iR++;
      } // else if

      else {                            //common to both maps
        iL->second.prob *= iR->second.prob;
        iL++; iR++;
      } // else

    } // while

    while ( iL != mapL.end() ){         // some extra spec's in mapL
      nDefsR--;
      iL->second.prob *= defR;
      iL++;
    } // while

    while ( iR != mapR.end() ){         // some extra spec's in mapR
      nDefsL--;
      double prob = defL * iR->second.prob;
      if (std::abs(prob - newDef) > lhs.margin) {
        mapL.emplace_hint(iL, iR->first, prob);
      } // if
      iR++;
    } // while

    PRLITE_ASSERT(
      nDefsL == nDefsR,
      "DiscreteTable_InplaceAbsorb : The number of cases in neither map should match, no matter from which side you look at it: " << nDefsL << "!=" << nDefsR << std::endl << lhs << std::endl << rhs);

    // there are some unused default values
    defL = newDef;

    return;
  } // if

  // Anything reaching here will have a non-uniform rhs

  // Case 3: Indien die lhs defVal > 0.0 itereer ons oor die *rhs*
  // sparse waardes. Verwerking minder as O( N*log2(N) ). Die hoop is
  // dat die vloerwaarde N en M klein genoeg kan maak om dit die
  // moeite werd te maak. Vir rhs Cat1EXPF sal ons ook soortgelyk moet doen.

  if (lhs.defProb > 0.0) { // we have a non-zero default on the lhs probs

    //std::cout << "Case 3*\n";

    // Elkeen van daardie veranderlike assignments word dan (een vir
    // een) uitgebrei na al die moontlike versoenbare lhs veranderlike
    // assignments. Daarna soek ons in die lhs vir elk van hierdie
    // assignments. Indien daar, vermenigvuldig die lhs prob met die
    // rhs prob en stel die flag in die lhs. Indien nie daar nie,
    // vermenigvuldig die rhs prob met die lhs defVal en insert 'n
    // nuwe inskrywing in die lhs sparse values, en stel die flag.

    // initialize stuff outside the iterator loop
    std::vector<size_t> newSizes(lni.size() );
    for (size_t k = 0; k < lni.size(); k++) {
      newSizes[k] = lhs.domains[lni[k]]->size();
    } // for
    std::vector<size_t> initIndices = initTuples<size_t>(newSizes);
    std::vector<size_t> indices( initIndices.size() );
    std::vector<AssignType> jointAssign(lhs.vars.size() );
    double newDef = defL*defR;

    auto iR = mapR.begin();
    while ( iR != mapR.end() ) {

      // first copy the values of the current rhs assignment
      for (size_t r = 0; r < rhs.vars.size(); r++) {
        jointAssign[ i2l[r] ] = iR->first[r];
      } // for

      // now extend it with the missing lhs variable assignments
      indices = initIndices;
      do {

        // copy the new variable assignments into place
        for (unsigned k = 0; k < ni2l.size(); k++) {
          jointAssign[ ni2l[k] ] = (*lhs.domains[ lni[k] ])[ indices[k] ];
        } // for

        auto result = mapL.emplace(jointAssign, defL); // leaves existing map records unchanged
        result.first->second.prob *= iR->second.prob;
        if (std::abs(result.first->second.prob - newDef) <= lhs.margin) {
          mapL.erase(result.first);
        } // if
        else {
          result.first->second.flag = 1;
        } // else

      } while( nextCombination(newSizes, indices) );

      iR++;
    } // while

    // Daarna itereer ons deur al die lhs sparse entries. flag == 1
    // beteken dat reeds met rhs prob vermenigvuldig is. Met flag == 0
    // moet ons nog vermenigvuldig met die defR.

    for (auto& aProb : mapL) {
      if (!aProb.second.flag) { // cases not visible from mapR
        aProb.second.prob *= defR;
      } // if
      aProb.second.flag = 0;
    }  // for

    // en heel laaste moet ons die lhs default ook opdateer
    defL = newDef;

    return;
  } // if

  // Case 4: Gevalle waar die lhs defVal == 0, hanteer ons net soos in
  // SparseTable maw ons itereer oor die *lhs* sparse waardes, onttrek
  // die tersake veranderlikes in die rhs en doen dan 'n find op die
  // rhs. As ons dit vind word die lhs prob met ooreenstemmende rhs
  // prob vermenigvuldig, anders word dit met rhs defVal
  // vermenigvuldig. O( M*log2(N) ) QED

  //std::cout << "Case 4*\n";

  auto iL = mapL.begin();
  while ( iL != mapL.end() ) {
    auto iR = mapR.find( extract<AssignType>(iL->first, l2i) );
    if ( iR != mapR.end() ) {
      iL->second.prob *= iR->second.prob;
    } // if
    else {iL->second.prob *= defR;} // else
    iL++;
  } // while
  return;

} // inplaceAbsorb

/*
  Hier is twee faktore ter sprake, ons verwys na hulle as lhs en
  rhs. def is die default waarde, sparse is die detail waardes in die
  map. Groottes van die sparse gedeeltes dui ons aan as M en N. Ons
  gaan slegs gevalle waar rhs 'n subversameling van lhs is,
  hanteer. Dis al wat ons werklik effektief inplace kan doen. Ander
  gevalle moet in DiscreteTable_Absorb hanteer word.
*/
template <typename AssignType>
void
DiscreteTable_InplaceAbsorb<AssignType>::inplaceProcess(
  DiscreteTable<AssignType>* lhsPtr,
  const Factor* rhsFPtr) {

  // rhsFPtr must be a compatible DiscreteTable*
  const DiscreteTable<AssignType>* dwnPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "DiscreteTable_InplaceAbsorb : rhs not a valid DiscreteTable : " << typeid(*dwnPtr).name() );
  // Create references to DiscreteTable instances
  const DiscreteTable<AssignType>& rhs(*dwnPtr);
  DiscreteTable<AssignType>& lhs(*lhsPtr);

  // l2i  the indices of vecL that is in the intersection
  // lni  the indices of unique vals in L (i.e. not in the intersection)
  // NOT USED: r2i   the indices of vecR that is in the intersection
  // rni  the indices of unique vals in R (i.e. not in the intersection)
  // i2l  where in the union is l2i[k] located
  // ni2l where in the union is lni[k] located
  // NOT USED: ni2r where in the union is rni[k] located

  // NOT USED: uVars

  std::vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  cupcap(lhs.vars,rhs.vars,l2i,lni,r2i,rni,i2l,ni2l,ni2r);
  PRLITE_ASSERT(!rni.size(),"DiscreteTable_InplaceAbsorb : the rhs variables must be a subset of those in the lhs");

    //std::cout << "Case X\n";

  // !!! NOTE: we could consider implementing a stripped down
  // intersection returning only l2i,lni,rni,i2l,ni2l
  inplaceAbsorb(lhs,rhs,std::move(l2i),std::move(lni),std::move(i2l),std::move(ni2l));

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Absorb<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Absorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One DiscreteTable multiplies with another
template <typename AssignType>
Factor* DiscreteTable_Absorb<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr, const Factor* rhsFPtr) {

  //typedef std::vector<AssignType> std::vector<AssignType>;

  // rhsFPtr must be a compatible DiscreteTable*
  const DiscreteTable<AssignType>* rhsPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr,
         "DiscreteTable_Absorb : rhs not a valid DiscreteTable : " << typeid(*rhsPtr).name() );

  // l2i  the indices of vecL that are in the intersection
  // lni  the indices of unique vals in L (i.e. not in the intersection)
  // r2i  the indices of vecR that are in the intersection
  // rni  the indices of unique vals in R (i.e. not in the intersection)
  // i2l  where in the union is l2i[k] located
  // ni2l where in the union is lni[k] located
  // ni2r where in the union is rni[k] located

  std::vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  std::vector<emdw::RVIdType> vecU =
    cupcap(lhsPtr->vars,rhsPtr->vars,l2i,lni,r2i,rni,i2l,ni2l,ni2r);

  if ( !rni.size() ) { // rhs vars is a subset of lhs

    //std::cout << "Case A\n";

    // make a copy of the lhs factor and an inplace absorber for it
    DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(*lhsPtr);
    DiscreteTable_InplaceAbsorb<AssignType> ipAbsorb;

    try{
      ipAbsorb.inplaceAbsorb(*lPtr, *rhsPtr, std::move(l2i), std::move(lni), std::move(i2l), std::move(ni2l));
    } // try
    catch (const char* s) {
      std::cout << __FILE__ << __LINE__
                << " call to 'inplaceAbsorb' failed" << std::endl;
      std::cout << s << std::endl;
      throw s;
    } // catch

    // and return the result
    return lPtr;
  } // if

  if ( !lni.size() ) { // lhs vars is a subset of rhs

    //std::cout << "Case B\n";

    // make a copy of the rhs factor and and inplace absorber for it
    DiscreteTable<AssignType>* rPtr = new DiscreteTable<AssignType>(*rhsPtr);
    DiscreteTable_InplaceAbsorb<AssignType> ipAbsorb;

    try{
      ipAbsorb.inplaceAbsorb(*rPtr, *lhsPtr, std::move(r2i), std::move(rni), std::move(i2l), std::move(ni2r));
    } // try
    catch (const char* s) {
      std::cout << __FILE__ << __LINE__
                << " call to 'inplaceAbsorb' failed" << std::endl;
      std::cout << s << std::endl;
      throw s;
    } // catch

    // and return the result
    return rPtr;
  } // if

  // no subsets

  // build the RV domains for the union of lhs and rhs
  std::vector< rcptr< std::vector<AssignType> > > theDomains( vecU.size() );
  for (size_t k = 0; k < l2i.size(); k++) {
    theDomains[ i2l[k] ] = lhsPtr->domains[ l2i[k] ];
  } // for
  for (size_t k = 0; k < lni.size(); k++) {
    theDomains[ ni2l[k] ] = lhsPtr->domains[ lni[k] ];
  } // for
  for (size_t k = 0; k < rni.size(); k++) {
    theDomains[ ni2r[k] ] = rhsPtr->domains[ rni[k] ];
  } // for


  std::map< std::vector<AssignType>, FProb> mapL;

  // first the case with zero default values - we only need to get the
  // product of the two sets of sparseProbs
  if ( lhsPtr->defProb == 0.0 && rhsPtr->defProb == 0.0 ) {

    //Multiply the DT of the lhs with the rhs exploiting zero-default probabilities
    for(auto &lhsIter: lhsPtr->sparseProbs){
      for(auto &rhsIter: rhsPtr->sparseProbs){
        unsigned i;
        for( i = 0; i < l2i.size(); i++) {
          //end this loop if any two elements doesn't match
          if (lhsIter.first[l2i[i]] != rhsIter.first[r2i[i]]) break;
        }
        //all lhs and rhs intersected vars match
        if(i==l2i.size()) {
          std::vector<AssignType> prodVars(vecU.size());

          for (size_t k = 0; k < l2i.size(); k++) prodVars[i2l[k] ] = lhsIter.first[l2i[k]];
          for (size_t k = 0; k < lni.size(); k++) prodVars[ni2l[k]] = lhsIter.first[lni[k]];
          for (size_t k = 0; k < rni.size(); k++) prodVars[ni2r[k]] = rhsIter.first[rni[k]];
          mapL[prodVars] = lhsIter.second.prob * rhsIter.second.prob;
        }
      }
    }

    // and rebuild a new factor for the new sparseProbs
    DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(
            std::move(vecU), std::move(theDomains),
            lhsPtr->defProb, std::move(mapL),
            lhsPtr->margin, lhsPtr->relFloor, true,
            lhsPtr->marginalizer, lhsPtr->inplaceNormalizer, lhsPtr->normalizer,
            lhsPtr->inplaceAbsorber, lhsPtr->absorber,
            lhsPtr->inplaceCanceller, lhsPtr->canceller,
            lhsPtr->observeAndReducer,lhsPtr->inplaceDamper,lhsPtr->sampler);

    return lPtr;

  } // if

  // no subsets and non-zero default - the slowest general case
  else {
    // extend lhs to include the scope of rhs
    mapL = extendSparseTable(lhsPtr->sparseProbs, lhsPtr->vars, rhsPtr->vars, rhsPtr->domains,
                              l2i, lni, r2i, rni, i2l, ni2l, ni2r);
  } // else


  // extend the lhs factor to all the variables of the rhs
  DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(
    std::move(vecU), std::move(theDomains),
    lhsPtr->defProb, std::move(mapL),
    lhsPtr->margin, lhsPtr->relFloor, true,
    lhsPtr->marginalizer, lhsPtr->inplaceNormalizer, lhsPtr->normalizer,
    lhsPtr->inplaceAbsorber, lhsPtr->absorber,
    lhsPtr->inplaceCanceller, lhsPtr->canceller,
    lhsPtr->observeAndReducer,lhsPtr->inplaceDamper,lhsPtr->sampler);

  DiscreteTable_InplaceAbsorb<AssignType> ipAbsorb;

  // l2i changes: depending on where the rhs variables fit relative to the old lhs ones
  // lni changes: depending on where the rhs variables fit relative to the old lhs ones
  // i2l unchanged: still the same positions in the intersection
  // ni2l changes: depending on where the rhs variables fit relative to the old lhs ones

  // We can either cleverly directly fix the above fields, or lazily
  // opt for ipAbsorb to figure it out. Since we suspect that this
  // part of the code will only be used rarely, we opt for the latter.

  try{
    ipAbsorb.inplaceProcess(lPtr, rhsPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

    // and return the result
  return lPtr;

} // process

// ***************************************************************************
// *************************** Family 2 **************************************
// *************************** Absorb0UF *************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceAbsorb0UF<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceAbsorb0UF<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
  Hier is twee faktore ter sprake, ons verwys na hulle as lhs en
  rhs. def is die default waarde, sparse is die detail waardes in die
  map. Groottes van die sparse gedeeltes dui ons aan as M en N. Ons
  gaan slegs gevalle waar rhs 'n subversameling van lhs is,
  hanteer. Dis al wat ons werklik effektief inplace kan doen. Ander
  gevalle moet in DiscreteTable_Absorb hanteer word.
*/

// NOTE: We should consider two further optimizations:

// 1: with relFloor > 0 cases: We can keep track of the running
// product max and remove / not insert cases into the map that, based
// on the current running max, we already know will fall below the
// threshold value.

// 2: Similarly we can also keep the cases where the product is within
// a "margin" margin of the default value, out of the sparseProbs map
// (instead of putting it in there and then rely on normalize to take
// it out again.

template <typename AssignType>
void
DiscreteTable_InplaceAbsorb0UF<AssignType>::inplaceAbsorb(
  DiscreteTable<AssignType>& lhs,
  const DiscreteTable<AssignType>& rhs,
  const std::vector<size_t>& l2i,
  const std::vector<size_t>& lni,
  const std::vector<size_t>& i2l,
  const std::vector<size_t>& ni2l) {

  //typedef std::vector<AssignType> std::vector<AssignType>;

  const std::map< std::vector<AssignType>, FProb>& mapR(rhs.sparseProbs);
  const double defR = rhs.defProb;
  std::map< std::vector<AssignType>, FProb>& mapL(lhs.sparseProbs);
  double& defL = lhs.defProb;

  // Case 1: with a uniform rhs we can simply return the lhs.
  if ( !mapR.size() ) {
    //std::cout << "Case 1*\n";

    // one can also leave the following block out altogether, in that
    // case the result may be unnormalized
    defL = prodUnderflowProtected(defL, defR);
    auto iL = mapL.begin();
    while ( iL != mapL.end() ){
      iL->second.prob = prodUnderflowProtected(iL->second.prob, defR);
      iL++;
    } // while
    return;
  } // if

  // Case 2: If we have exactly the same variables on both sides, the
  // double ladder trick should be fastest. Rather useful since
  // this happens a lot in the LBU algo.

  if (lhs.vars == rhs.vars){

    //std::cout << "Case 2*\n";

    // now multiply by simultaneously climbing two ladders

    // to calc the number of defaults occurring in neither of the maps
    size_t nDefsL = lhs.maxSz-mapL.size();
    size_t nDefsR = rhs.maxSz;
    double newDef = prodUnderflowProtected(defL, defR);

    for (auto k : lni) {nDefsR *= lhs.domains[k]->size();} // for
    nDefsR -= mapR.size();

    auto iL = mapL.begin();
    auto iR = mapR.begin();
    while ( iL != mapL.end() && iR != mapR.end() ){

      if (iL->first < iR->first) {      // only in mapL
        nDefsR--;
        iL->second.prob = prodUnderflowProtected(iL->second.prob, defR);
        iL++;
      } // if

      else if (iL->first > iR->first) { // only in mapR
        nDefsL--;
        double prob = prodUnderflowProtected(defL, iR->second.prob);
        if (std::abs(prob - newDef) > lhs.margin) {
          mapL.emplace_hint(iL, iR->first, prob);
        } // if
        iR++;
      } // else if

      else {                            //common to both maps
        iL->second.prob = prodUnderflowProtected(iL->second.prob, iR->second.prob);
        iL++; iR++;
      } // else

    } // while

    while ( iL != mapL.end() ){         // some extra spec's in mapL
      nDefsR--;
      iL->second.prob = prodUnderflowProtected(iL->second.prob, defR);
      iL++;
    } // while

    while ( iR != mapR.end() ){         // some extra spec's in mapR
      nDefsL--;
      double prob = prodUnderflowProtected(defL, iR->second.prob);
      if (std::abs(prob - newDef) > lhs.margin) {
        mapL.emplace_hint(iL, iR->first, prob);
      } // if
      iR++;
    } // while

    PRLITE_ASSERT(
      nDefsL == nDefsR,
      "DiscreteTable_InplaceAbsorb0UF : The number of cases in neither map should match, no matter from which side you look at it: " << nDefsL << "!=" << nDefsR << std::endl << lhs << std::endl << rhs);

    // there are some unused default values
    defL = newDef;

    return;
  } // if

  // Anything reaching here will have a non-uniform rhs

  // Case 3: Indien die lhs defVal > 0.0 itereer ons oor die *rhs*
  // sparse waardes. Verwerking minder as O( N*log2(N) ). Die hoop is
  // dat die vloerwaarde N en M klein genoeg kan maak om dit die
  // moeite werd te maak.

  if (lhs.defProb > 0.0) { // we have a non-zero default on the lhs probs

    //std::cout << "Case 3*\n";

    // Elkeen van daardie veranderlike assignments word dan (een vir
    // een) uitgebrei na al die moontlike versoenbare lhs veranderlike
    // assignments. Daarna soek ons in die lhs vir elk van hierdie
    // assignments. Indien daar, vermenigvuldig die lhs prob met die
    // rhs prob en stel die flag in die lhs. Indien nie daar nie,
    // vermenigvuldig die rhs prob met die lhs defVal en insert 'n
    // nuwe inskrywing in die lhs sparse values, en stel die flag.

    // initialize stuff outside the iterator loop
    std::vector<size_t> newSizes(lni.size() );
    for (size_t k = 0; k < lni.size(); k++) {
      newSizes[k] = lhs.domains[lni[k]]->size();
    } // for
    std::vector<size_t> initIndices = initTuples<size_t>(newSizes);
    std::vector<size_t> indices( initIndices.size() );
    std::vector<AssignType> jointAssign(lhs.vars.size() );
    double newDef = prodUnderflowProtected(defL, defR);

    auto iR = mapR.begin();
    while ( iR != mapR.end() ) {

      // first copy the values of the current rhs assignment
      for (size_t r = 0; r < rhs.vars.size(); r++) {
        jointAssign[ i2l[r] ] = iR->first[r];
      } // for

      // now extend it with the missing lhs variable assignments
      indices = initIndices;
      do {

        // copy the new variable assignments into place
        for (unsigned k = 0; k < ni2l.size(); k++) {
          jointAssign[ ni2l[k] ] = (*lhs.domains[ lni[k] ])[ indices[k] ];
        } // for

        auto result = mapL.emplace(jointAssign, defL); // leaves existing map records unchanged
        result.first->second.prob = prodUnderflowProtected(result.first->second.prob, iR->second.prob);
        if (std::abs(result.first->second.prob - newDef) <= lhs.margin) {
          mapL.erase(result.first);
        } // if
        else {
          result.first->second.flag = 1;
        } // else

      } while( nextCombination(newSizes, indices) );

      iR++;
    } // while

    // Daarna itereer ons deur al die lhs sparse entries. flag == 1
    // beteken dat reeds met rhs prob vermenigvuldig is. Met flag == 0
    // moet ons nog vermenigvuldig met die defR.

    for (auto& aProb : mapL) {
      if (!aProb.second.flag) { // cases not visible from mapR
        aProb.second.prob = prodUnderflowProtected(aProb.second.prob, defR);
      } // if
      aProb.second.flag = 0;
    }  // for

    // en heel laaste moet ons die lhs default ook opdateer
    defL = newDef;

    return;
  } // if

  // Case 4: Gevalle waar die lhs defVal == 0, hanteer ons net soos in
  // SparseTable maw ons itereer oor die *lhs* sparse waardes, onttrek
  // die tersake veranderlikes in die rhs en doen dan 'n find op die
  // rhs. As ons dit vind word die lhs prob met ooreenstemmende rhs
  // prob vermenigvuldig, anders word dit met rhs defVal
  // vermenigvuldig. O( M*log2(N) ) QED

  //std::cout << "Case 4*\n";

  auto iL = mapL.begin();
  while ( iL != mapL.end() ) {
    auto iR = mapR.find( extract<AssignType>(iL->first, l2i) );
    if ( iR != mapR.end() ) {
      iL->second.prob = prodUnderflowProtected(iL->second.prob, iR->second.prob);
    } // if
    else {iL->second.prob = prodUnderflowProtected(iL->second.prob, defR);} // else
    iL++;
  } // while
  return;

} // inplaceAbsorb

/*
  Hier is twee faktore ter sprake, ons verwys na hulle as lhs en
  rhs. def is die default waarde, sparse is die detail waardes in die
  map. Groottes van die sparse gedeeltes dui ons aan as M en N. Ons
  gaan slegs gevalle waar rhs 'n subversameling van lhs is,
  hanteer. Dis al wat ons werklik effektief inplace kan doen. Ander
  gevalle moet in DiscreteTable_Absorb hanteer word.
*/
template <typename AssignType>
void
DiscreteTable_InplaceAbsorb0UF<AssignType>::inplaceProcess(
  DiscreteTable<AssignType>* lhsPtr,
  const Factor* rhsFPtr) {

  // rhsFPtr must be a compatible DiscreteTable*
  const DiscreteTable<AssignType>* dwnPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "DiscreteTable_InplaceAbsorb0UF : rhs not a valid DiscreteTable : " << typeid(*dwnPtr).name() );
  // Create references to DiscreteTable instances
  const DiscreteTable<AssignType>& rhs(*dwnPtr);
  DiscreteTable<AssignType>& lhs(*lhsPtr);

  // l2i  the indices of vecL that is in the intersection
  // lni  the indices of unique vals in L (i.e. not in the intersection)
  // NOT USED: r2i   the indices of vecR that is in the intersection
  // rni  the indices of unique vals in R (i.e. not in the intersection)
  // i2l  where in the union is l2i[k] located
  // ni2l where in the union is lni[k] located
  // NOT USED: ni2r where in the union is rni[k] located

  // NOT USED: uVars

  std::vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  cupcap(lhs.vars,rhs.vars,l2i,lni,r2i,rni,i2l,ni2l,ni2r);
  PRLITE_ASSERT(!rni.size(),"DiscreteTable_InplaceAbsorb0UF : the rhs variables must be a subset of those in the lhs");

    //std::cout << "Case X\n";

  // !!! NOTE: we could consider implementing a stripped down
  // intersection returning only l2i,lni,rni,i2l,ni2l
  inplaceAbsorb(lhs,rhs,std::move(l2i),std::move(lni),std::move(i2l),std::move(ni2l));

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Absorb0UF<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Absorb0UF<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One DiscreteTable multiplies with another
template <typename AssignType>
Factor* DiscreteTable_Absorb0UF<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr, const Factor* rhsFPtr) {

  //typedef std::vector<AssignType> std::vector<AssignType>;

  // rhsFPtr must be a compatible DiscreteTable*
  const DiscreteTable<AssignType>* rhsPtr = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr,
         "DiscreteTable_Absorb0UF : rhs not a valid DiscreteTable : " << typeid(*rhsPtr).name() );

  // l2i  the indices of vecL that are in the intersection
  // lni  the indices of unique vals in L (i.e. not in the intersection)
  // r2i  the indices of vecR that are in the intersection
  // rni  the indices of unique vals in R (i.e. not in the intersection)
  // i2l  where in the union is l2i[k] located
  // ni2l where in the union is lni[k] located
  // ni2r where in the union is rni[k] located

  std::vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  std::vector<emdw::RVIdType> vecU =
    cupcap(lhsPtr->vars,rhsPtr->vars,l2i,lni,r2i,rni,i2l,ni2l,ni2r);

  if ( !rni.size() ) { // rhs vars is a subset of lhs

    //std::cout << "Case A\n";

    // make a copy of the lhs factor and an inplace absorber for it
    DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(*lhsPtr);
    DiscreteTable_InplaceAbsorb0UF<AssignType> ipAbsorb;

    try{
      ipAbsorb.inplaceAbsorb(*lPtr, *rhsPtr, std::move(l2i), std::move(lni), std::move(i2l), std::move(ni2l));
    } // try
    catch (const char* s) {
      std::cout << __FILE__ << __LINE__
                << " call to 'inplaceAbsorb' failed" << std::endl;
      std::cout << s << std::endl;
      throw s;
    } // catch

    // and return the result
    return lPtr;
  } // if

  if ( !lni.size() ) { // lhs vars is a subset of rhs

    //std::cout << "Case B\n";

    // make a copy of the rhs factor and and inplace absorber for it
    DiscreteTable<AssignType>* rPtr = new DiscreteTable<AssignType>(*rhsPtr);
    DiscreteTable_InplaceAbsorb0UF<AssignType> ipAbsorb;

    try{
      ipAbsorb.inplaceAbsorb(*rPtr, *lhsPtr, std::move(r2i), std::move(rni), std::move(i2l), std::move(ni2r));
    } // try
    catch (const char* s) {
      std::cout << __FILE__ << __LINE__
                << " call to 'inplaceAbsorb' failed" << std::endl;
      std::cout << s << std::endl;
      throw s;
    } // catch

    // and return the result
    return rPtr;
  } // if

  // no subsets

  // build the RV domains for the union of lhs and rhs
  std::vector< rcptr< std::vector<AssignType> > > theDomains( vecU.size() );
  for (size_t k = 0; k < l2i.size(); k++) {
    theDomains[ i2l[k] ] = lhsPtr->domains[ l2i[k] ];
  } // for
  for (size_t k = 0; k < lni.size(); k++) {
    theDomains[ ni2l[k] ] = lhsPtr->domains[ lni[k] ];
  } // for
  for (size_t k = 0; k < rni.size(); k++) {
    theDomains[ ni2r[k] ] = rhsPtr->domains[ rni[k] ];
  } // for


  std::map< std::vector<AssignType>, FProb> mapL;

  // first the case with zero default values - we only need to get the
  // product of the two sets of sparseProbs
  if ( lhsPtr->defProb == 0.0 && rhsPtr->defProb == 0.0 ) {

    //Multiply the DT of the lhs with the rhs exploiting zero-default probabilities
    for(auto &lhsIter: lhsPtr->sparseProbs){
      for(auto &rhsIter: rhsPtr->sparseProbs){
        unsigned i;
        for( i = 0; i < l2i.size(); i++) {
          //end this loop if any two elements doesn't match
          if (lhsIter.first[l2i[i]] != rhsIter.first[r2i[i]]) break;
        }
        //all lhs and rhs intersected vars match
        if(i==l2i.size()) {
          std::vector<AssignType> prodVars(vecU.size());

          for (size_t k = 0; k < l2i.size(); k++) prodVars[i2l[k] ] = lhsIter.first[l2i[k]];
          for (size_t k = 0; k < lni.size(); k++) prodVars[ni2l[k]] = lhsIter.first[lni[k]];
          for (size_t k = 0; k < rni.size(); k++) prodVars[ni2r[k]] = rhsIter.first[rni[k]];
          mapL[prodVars] = prodUnderflowProtected(lhsIter.second.prob, rhsIter.second.prob);
        }
      }
    }

    // and rebuild a new factor for the new sparseProbs
    DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(
            std::move(vecU), std::move(theDomains),
            lhsPtr->defProb, std::move(mapL),
            lhsPtr->margin, lhsPtr->relFloor, true,
            lhsPtr->marginalizer, lhsPtr->inplaceNormalizer, lhsPtr->normalizer,
            lhsPtr->inplaceAbsorber, lhsPtr->absorber,
            lhsPtr->inplaceCanceller, lhsPtr->canceller,
            lhsPtr->observeAndReducer,lhsPtr->inplaceDamper,lhsPtr->sampler);

    return lPtr;

  } // if

  // no subsets and non-zero default - the slowest general case
  else {
    // extend lhs to include the scope of rhs
    mapL = extendSparseTable(lhsPtr->sparseProbs, lhsPtr->vars, rhsPtr->vars, rhsPtr->domains,
                              l2i, lni, r2i, rni, i2l, ni2l, ni2r);
  } // else


  // extend the lhs factor to all the variables of the rhs
  DiscreteTable<AssignType>* lPtr = new DiscreteTable<AssignType>(
    std::move(vecU), std::move(theDomains),
    lhsPtr->defProb, std::move(mapL),
    lhsPtr->margin, lhsPtr->relFloor, true,
    lhsPtr->marginalizer, lhsPtr->inplaceNormalizer, lhsPtr->normalizer,
    lhsPtr->inplaceAbsorber, lhsPtr->absorber,
    lhsPtr->inplaceCanceller, lhsPtr->canceller,
    lhsPtr->observeAndReducer,lhsPtr->inplaceDamper,lhsPtr->sampler);

  DiscreteTable_InplaceAbsorb0UF<AssignType> ipAbsorb;

  // l2i changes: depending on where the rhs variables fit relative to the old lhs ones
  // lni changes: depending on where the rhs variables fit relative to the old lhs ones
  // i2l unchanged: still the same positions in the intersection
  // ni2l changes: depending on where the rhs variables fit relative to the old lhs ones

  // We can either cleverly directly fix the above fields, or lazily
  // opt for ipAbsorb to figure it out. Since we suspect that this
  // part of the code will only be used rarely, we opt for the latter.

  try{
    ipAbsorb.inplaceProcess(lPtr, rhsPtr);
  } // try
  catch (const char* s) {
    std::cout << __FILE__ << __LINE__
              << " call to 'inplaceProcess' failed" << std::endl;
    std::cout << s << std::endl;
    throw s;
  } // catch

    // and return the result
  return lPtr;

} // process

// ***************************************************************************
// *************************** Family 2 **************************************
// *************************** Cancel ****************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceCancel<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceCancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA


// lhs must contain rhs. returns lhs /= rhs;
template <typename AssignType>
void
DiscreteTable_InplaceCancel<AssignType>::inplaceCancel(
  DiscreteTable<AssignType>& lhs,
  const DiscreteTable<AssignType>& rhs,
  const std::vector<size_t>& l2i,
  const std::vector<size_t>& lni,
  const std::vector<size_t>& i2l,
  const std::vector<size_t>& ni2l) {

  //typedef std::vector<AssignType> std::vector<AssignType>;

  const std::map< std::vector<AssignType>, FProb>& mapR(rhs.sparseProbs);
  const double defR = rhs.defProb;
  std::map< std::vector<AssignType>, FProb>& mapL(lhs.sparseProbs);
  double& defL = lhs.defProb;

  // Case 1: with a uniform rhs we can simply return the lhs.
  if ( !mapR.size() ) {
    //std::cout << "Case 1/\n";

    // one can also leave the following block out altogether, in that
    // case the result may be unnormalized
    defL /= defR;
    auto iL = mapL.begin();
    while ( iL != mapL.end() ){
      iL->second.prob /= defR;
      iL++;
    } // while

    return;
  } // if

  // Case 2: If we have exactly the same variables on both sides, the
  // double ladder trick should be fastest. Rather useful since
  // this happens all the time in the LBU algo.
  if (lhs.vars.size() == rhs.vars.size() ){

    //std::cout << "Case 2/\n";

    // now divide by simultaneously climbing two ladders

    // to calc the number of defaults occurring in neither of the maps
    size_t nDefsL = lhs.maxSz-mapL.size();
    size_t nDefsR = rhs.maxSz;
    double newDef = defR > 0.0 ? defL/defR : defL;

    for (auto k : lni) {nDefsR *= lhs.domains[k]->size();} // for
    nDefsR -= mapR.size();

    auto iL = mapL.begin();
    auto iR = mapR.begin();
    while ( iL != mapL.end() && iR != mapR.end() ){

      if (iL->first < iR->first) {      // only in mapL
        nDefsR--;
        if(defR > 0) {
          iL->second.prob /= defR;
        } // if
        iL++;
      } // if

      else if (iL->first > iR->first) { // only in mapR
        nDefsL--;
        double prob = iR->second.prob > 0.0 ? defL / iR->second.prob : 0.0;
        if (std::abs(prob - newDef) > lhs.margin) {
          mapL.emplace_hint(iL, iR->first, prob);
        } // if
        iR++;
      } // else if

      else {                            //common to both maps
        if (iR->second.prob > 0.0) {
          iL->second.prob /= iR->second.prob;
        } // if
        iL++; iR++;
      } // else

    } // while

    while ( iL != mapL.end() ){         // some extra spec's in mapL
      nDefsR--;
      if (defR > 0.0) {
        iL->second.prob /= defR;
      } // if
      iL++;
    } // while

    while ( iR != mapR.end() ){         // some extra spec's in mapR
      nDefsL--;
      double prob = iR->second.prob > 0.0 ? defL / iR->second.prob : 0.0;
      if (std::abs(prob - newDef) > lhs.margin) {
        mapL.emplace_hint(iL, iR->first, prob);
      } // if
      iR++;
    } // while

    PRLITE_ASSERT(
      nDefsL == nDefsR,
      "DiscreteTable_InplaceCancel : The number of cases in neither map should match, no matter from which side you look at it.");

    // there are some unused default values
    defL = newDef;

    return;
  } // if

  // Anything reaching here will have a non-uniform rhs

  // Case 3: Indien die lhs defVal > 0.0 itereer ons oor die *rhs*
  // sparse waardes. Verwerking minder as O( N*log2(N) ). Die hoop is
  // dat die vloerwaarde N en M klein genoeg kan maak om dit die
  // moeite werd te maak.

  if (lhs.defProb > 0.0) { // we have a non-zero default on the lhs probs

    //std::cout << "Case 3/\n";

    // Elkeen van daardie veranderlike assignments word dan (een vir
    // een) uitgebrei na al die moontlike versoenbare veranderlike
    // assignments van die volledige stel veranderlikes (maw die
    // veranderlikes van die lhs). Daarna soek ons in die lhs vir elk
    // van hierdie assignments. Indien daar, deel die lhs prob met die
    // rhs prob en stel die flag in die lhs. Indien nie daar nie, deel
    // die rhs prob met die lhs defVal en insert 'n nuwe inskrywing in
    // die lhs sparse values, en stel die flag.

    // initialize stuff outside the iterator loop
    std::vector<size_t> newSizes(lni.size() );
    for (size_t k = 0; k < lni.size(); k++) {
      newSizes[k] = lhs.domains[lni[k]]->size();
    } // for
    std::vector<size_t> initIndices = initTuples<size_t>(newSizes);
    std::vector<size_t> indices( initIndices.size() );
    std::vector<AssignType> jointAssign(lhs.vars.size() );
    // to distinguish from freshly created entries (flag == 0)
    //for (auto& aProb : mapL) {aProb.second.flag = 1;}
    double newDef = defR > 0.0 ? defL/defR : defL;

    auto iR = mapR.begin();
    while ( iR != mapR.end() ) {

      // first copy the values of the current rhs assignment
      for (size_t r = 0; r < rhs.vars.size(); r++) {
        jointAssign[ i2l[r] ] = iR->first[r];
      } // for

      // now extend it with the missing lhs variable assignments
      indices = initIndices;
      do {

        // copy the new variable assignments into place
        for (unsigned k = 0; k < ni2l.size(); k++) {
          jointAssign[ ni2l[k] ] = (*lhs.domains[ lni[k] ])[ indices[k] ];
        } // for

        if (iR->second.prob > 0.0) {

          auto result = mapL.emplace(jointAssign, defL);  // leaves existing map records unchanged
          result.first->second.prob /= iR->second.prob;
          if (std::abs(result.first->second.prob - newDef) <= lhs.margin) {
            mapL.erase(result.first);
          } // if
          else {
            result.first->second.flag = 1;
          } // else

        } // if

        else {
          if (defL > 0.0) {
            auto& aProb = mapL[jointAssign];
            aProb.prob = 0.0;
            aProb.flag = 1;
          } // if
        } // else

      } while( nextCombination(newSizes, indices) );

      iR++;
    } // while

    // Daarna itereer ons deur al die lhs sparse entries. flag == 1
    // beteken dat reeds met rhs prob gedeel is. Met flag == 0, moet
    // ons nog deel met die defR.

    for (auto& aProb : mapL) {
      if (!aProb.second.flag) { // cases not visible from mapR
        if (defR > 0.0) {aProb.second.prob /= defR;} // if
      } // if
      aProb.second.flag = 0;
    }  // for

    // en heel laaste moet ons die lhs default ook opdateer
    defL = newDef;

    return;
  } // if

  // Case 4: Gevalle waar die lhs defVal == 0, hanteer ons net soos in
  // SparseTable maw ons itereer oor die *lhs* sparse waardes, onttrek
  // die tersake veranderlikes in die rhs en doen dan 'n find op die
  // rhs. As ons dit vind word die lhs prob deur ooreenstemmende rhs
  // prob gedeel, anders word dit deur rhs defVal
  // gedeel. O( M*log2(N) ) QED

  //std::cout << "Case 4/\n";

  auto iL = mapL.begin();
  while ( iL != mapL.end() ) {
    auto iR = mapR.find( extract<AssignType>(iL->first, l2i) );
    if ( iR != mapR.end() ) {
      if (iR->second.prob > 0.0) {
        iL->second.prob /= iR->second.prob;
      } // if
    } // if
    else if (defR > 0.0) {iL->second.prob /= defR;} // else if
    iL++;
  } // while
  return;

} // inplaceCancel

// lhs must contain rhs. returns lhs /= rhs;
template <typename AssignType>
void
DiscreteTable_InplaceCancel<AssignType>::inplaceProcess(
  DiscreteTable<AssignType>* lhsPtr,
  const Factor* rhsFPtr) {

  // rhsFPtr must be an compatible DiscreteTable*
  const DiscreteTable<AssignType>* dwnPtr =
    dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "DiscreteTable_InplaceCancel : rhs not a valid DiscreteTable : "
         << typeid(*dwnPtr).name() );
  // Create references to DiscreteTable instances
  const DiscreteTable<AssignType>& rhs(*dwnPtr);
  DiscreteTable<AssignType>& lhs(*lhsPtr);

  // l2i  the indices of vecL that is in the intersection
  // lni  the indices of unique vals in L (i.e. not in the intersection)
  // NOT USED: r2i   the indices of vecR that is in the intersection
  // rni  the indices of unique vals in R (i.e. not in the intersection)
  // i2l  where in the union is l2i[k] located
  // ni2l where in the union is lni[k] located
  // NOT USED: ni2r where in the union is rni[k] located

  // NOT USED: uVars

  std::vector<size_t> l2i,lni,r2i,rni,i2l,ni2l,ni2r;
  cupcap(lhs.vars,rhs.vars,l2i,lni,r2i,rni,i2l,ni2l,ni2r);
  PRLITE_ASSERT(!rni.size(),"DiscreteTable_InplaceCancel : the rhs variables must be a subset of those in the lhs");

    //std::cout << "Case X\n";

  // !!! NOTE: we could consider implementing a stripped down
  // intersection returning only l2i,lni,rni,i2l,ni2l
  inplaceCancel(lhs,rhs,std::move(l2i),std::move(lni),std::move(i2l),std::move(ni2l));
  return;
} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Cancel<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Cancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One DiscreteTable cancels another - this is implemented by subtracting the parameter
// of the one potential from that of the other
template <typename AssignType>
Factor* DiscreteTable_Cancel<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr, const Factor* rhsFPtr) {
  // make a copy of the factor
  DiscreteTable<AssignType>* fPtr = new DiscreteTable<AssignType>(*lhsPtr);

  // inplace cancel it
  DiscreteTable_InplaceCancel<AssignType> ipCancel;
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

// ***************************************************************************
// ******************************** Family 3 *********************************
// ******************************** Marginalize ******************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_Marginalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_Marginalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
 * 1) if defProb is 0, we can do this as for SparseTable - the cost is
 * O(N*log2(N)).
 *
 * 2) The alternative is to keep count of how many instances of the
 * toKeep variables should occur in total. After completion the
 * missing instances are added from the default values. The default
 * value itself is also multiplied by this number. The cost also is
 * roughly O(N log2(N)).
 */
template <typename AssignType>
Factor* DiscreteTable_Marginalize<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr,
                                                       const emdw::RVIds& toKeep,
                                                       bool presorted,
                                                       const Factor* peekAheadPtr) {

  //typedef std::vector<AssignType> std::vector<AssignType>;
  const DiscreteTable<AssignType>& lhs(*lhsPtr);
  if (toKeep.size() == 1) {presorted = true;}

  // fiddling about to avoid copying toKeep when there is no sorting required
  const emdw::RVIds* sortedKeepPtr = &toKeep;
  emdw::RVIds tmp; // need it here to not go out of scope too soon
  if (!presorted) {
    std::vector<size_t> sorted = sortIndices( toKeep, std::less<unsigned>() );
    tmp = extract<unsigned>(toKeep,sorted);
    sortedKeepPtr = &tmp;
  } // if
  const emdw::RVIds& sortedKeep(*sortedKeepPtr);

  // now to find out where the intersection between vars and theVars
  // lie. Any variables outside of these need not be considered.
  std::vector<size_t> l2i; // which var indices from lhs pot are also in sortedKeep
  std::vector<size_t> r2i; // the indices in sortedKeep where you will find them
  std::vector<size_t> lni; // which var indices in lhs are not in sortedKeep
  std::vector<size_t> rni; // // the indices in sortedKeep where you will find them
  emdw::RVIds theVars = sortedIntersection(lhs.vars, sortedKeep, l2i, r2i, lni, rni);

  // if theVars == vars, dont bother, nothing changes
  if ( theVars.size() == lhs.vars.size() ) {
    return lhs.copy();
  } // if

  // the number of instances we should see for each of the possible assignments
  size_t nVariants = lhs.domains[ lni[0] ]->size();
  for (size_t k = 1; k < lni.size(); k++) {
    nVariants *= lhs.domains[ lni[k] ]->size();
  } // for

  std::map<std::vector<AssignType>, FProb> theProbs;

  if (lhs.defProb > 0.0) {  // we have to take the default probs into account too

    //double newDef = lhs.defProb*nVariants;

    for (auto itr = lhs.sparseProbs.begin(); itr != lhs.sparseProbs.end(); itr++) {

      if (std::abs(itr->second.prob - lhs.defProb) > lhs.margin) { // to prevent redundant entries
        FProb& aProb = theProbs[extract<AssignType>(itr->first, l2i)];
        aProb.prob += itr->second.prob; // if never seen before this adds 0.0
        aProb.flag++;
      } // if

    } // for

    auto itr1 = theProbs.begin();
    while (itr1 != theProbs.end() ) {
      size_t nMissed = nVariants - itr1->second.flag;
      if (nMissed > 0) {itr1->second.prob += nMissed*lhs.defProb;} // if
      itr1->second.flag = 0;
      itr1++;
    } // while

  } // if

  else{ // no default probs to worry about

    for (auto itr = lhs.sparseProbs.begin(); itr != lhs.sparseProbs.end(); itr++) {
      if (itr->second.prob > 0.0) { // to prevent redundant entries
        theProbs[extract<AssignType>(itr->first, l2i)].prob += itr->second.prob;
      } // if
    } // for

  } // else

  return new DiscreteTable<AssignType>(
    theVars,
    extract< rcptr< std::vector<AssignType> > >(lhs.domains, l2i),
    lhs.defProb*nVariants,    // note that we only change the default here
    theProbs,
    lhs.margin,
    lhs.relFloor,
    true,
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

} // process

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_MaxMarginalize<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_MaxMarginalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
 * 1) if defProb is 0, we can do this as for SparseTable - the cost is
 * O(N*log2(N)).
 *
 * 2) The alternative is to keep count of how many instances of the
 * toKeep variables should occur in total. After completion the
 * missing instances are added from the default values. The default
 * value itself is also multiplied by this number. The cost also is
 * roughly O(N log2(N)), but we also have to allocate a map with the
 * counts. This will be our first implementation.
 */
template <typename AssignType>
Factor* DiscreteTable_MaxMarginalize<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr,
                                                       const emdw::RVIds& toKeep,
                                                       bool presorted,
                                                          const Factor* peekAheadPtr) {

  //typedef std::vector<AssignType> std::vector<AssignType>;
  const DiscreteTable<AssignType>& lhs(*lhsPtr);
  if (toKeep.size() == 1) {presorted = true;}

  // fiddling about to avoid copying toKeep when there is no sorting required
  const emdw::RVIds* sortedKeepPtr = &toKeep;
  emdw::RVIds tmp; // need it here to not go out of scope too soon
  if (!presorted) {
    std::vector<size_t> sorted = sortIndices( toKeep, std::less<unsigned>() );
    tmp = extract<unsigned>(toKeep,sorted);
    sortedKeepPtr = &tmp;
  } // if
  const emdw::RVIds& sortedKeep(*sortedKeepPtr);

  // now to find out where the intersection between vars and theVars
  // lie. Any variables outside of these need not be considered.
  std::vector<size_t> l2i; // which var indices from lhs pot are also in sortedKeep
  std::vector<size_t> r2i; // the indices in sortedKeep where you will find them
  std::vector<size_t> lni; // which var indices in lhs are not in sortedKeep
  std::vector<size_t> rni; // // the indices in sortedKeep where you will find them
  emdw::RVIds theVars = sortedIntersection(lhs.vars, sortedKeep, l2i, r2i, lni, rni);

  // if theVars == vars, dont bother, nothing changes
  if ( theVars.size() == lhs.vars.size() ) {
    return lhs.copy();
  } // if

  // the number of instances we should see for each of the possible assignments
  size_t nVariants = lhs.domains[ lni[0] ]->size();
  for (size_t k = 1; k < lni.size(); k++) {
    nVariants *= lhs.domains[ lni[k] ]->size();
  } // for

  std::map<std::vector<AssignType>, FProb> theProbs;

  if (lhs.defProb > 0.0) {  // we have to take the default probs into account too

    for (auto itr = lhs.sparseProbs.begin(); itr != lhs.sparseProbs.end(); itr++) {
      if (std::abs(itr->second.prob - lhs.defProb) > lhs.margin) { // to prevent redundant entries
        FProb& aProb = theProbs[extract<AssignType>(itr->first, l2i)];
        if (itr->second.prob > aProb.prob) {  // if never seen before aProb.prob will be zero
          aProb.prob = itr->second.prob;
        } // if
        aProb.flag++;
      } // if
    } // for

    auto itr1 = theProbs.begin();
    while (itr1 != theProbs.end() ) {
      size_t nMissed = nVariants - itr1->second.flag;
      if (nMissed > 0 and lhs.defProb > itr1->second.prob) {
        itr1->second.prob = lhs.defProb;
      } // if
      itr1++;
    } // while

  } // if

  else{ // no default probs to worry about

    for (auto itr = lhs.sparseProbs.begin(); itr != lhs.sparseProbs.end(); itr++) {
      if (itr->second.prob > 0.0) { // to prevent redundant entries
        std::vector<AssignType> assig( extract<AssignType>(itr->first, l2i) );
        auto itrNew = theProbs.find(assig);
        if ( itrNew != theProbs.end() ) { // have seen this before
          if (itr->second.prob > itrNew->second.prob) {itrNew->second.prob = itr->second.prob;} // if
        } // if
        else {theProbs[assig] = itr->second;} // else
      } // if
    } // for

  } // else

  return new DiscreteTable<AssignType>(
    theVars,
    extract< rcptr< std::vector<AssignType> > >(lhs.domains, l2i),
    lhs.defProb,
    theProbs,
    lhs.margin,
    lhs.relFloor,
    true,
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

} // process

// ***************************************************************************
// ******************************** Family 4 *********************************
// ******************************** ObserveAndReduce *************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_ObserveAndReduce<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_ObserveAndReduce<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* DiscreteTable_ObserveAndReduce<AssignType>::process(const DiscreteTable<AssignType>* lhsPtr,
                                                            const emdw::RVIds& theVars,
                                                            const emdw::RVVals& anyVals,
                                                            bool presorted) {

  // A typedef for the actual internal assignment vector type
  //typedef std::vector<AssignType> std::vector<AssignType>;
  const DiscreteTable<AssignType>& lhs(*lhsPtr);
  if (theVars.size() == 1) {presorted = true;}

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const emdw::RVIds& sortedVars(*sortedVarsPtr);

  // now i need to find out where the intersection between vars and
  // sortedVars lie. Any variables outside of these need not be
  // considered.
  std::vector<int> l2i; // which indices in lhs pot are shared
  std::vector<int> r2i; // with corresponding indices in rhs pot
  std::vector<int> lni; // non-observed lhs variables ARE the new variable set
  std::vector<int> rni; // remaining rhs variables (unused)
  //emdw::RVIds varI = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);
  sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);

  // now collect the underlying observed variables. if not presorted
  // we'll need some gymnastics to collect them in the proper order
  if (!presorted) {
    // now fix r2i to give the actual index in anyVals directly
    for (unsigned i = 0; i < r2i.size(); i++) {
      r2i[i] = sorted[r2i[i]];
    } // for

  } // if
  // anyVals are still of type AnyType. While extracting them in the
  // right order, we at the same time convert their type to AssignType
  std::vector<AssignType> observedVals = extract<AssignType>(anyVals, r2i);

  // the new variables and table
  emdw::RVIds newVars = extract<unsigned>(lhs.vars,lni);
  std::map<std::vector<AssignType>, FProb> theProbs;

  for (auto itr = lhs.sparseProbs.begin(); itr != lhs.sparseProbs.end(); itr++) {
    if (extract<AssignType>(itr->first, l2i) == observedVals) {
      // keep it, but on the reduced scope
      theProbs[extract<AssignType>(itr->first,lni)] = itr->second;
    } // if
  } // for

  return new DiscreteTable<AssignType>(
    newVars,
    extract< rcptr< std::vector<AssignType> > >(lhs.domains, lni),
    lhs.defProb,
    theProbs,
    lhs.margin,
    lhs.relFloor,
    true,
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

} // process

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& DiscreteTable_InplaceDamping<AssignType>::isA() const {
  static const std::string CLASSNAME("DiscreteTable_InplaceDamping<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and old
// parameters - this linear combination is determined by the df parameter
template <typename AssignType>
double
DiscreteTable_InplaceDamping<AssignType>::inplaceProcess(
  DiscreteTable<AssignType>* newPtr,
  const Factor* oldPtr,
  double oldWeight) { // oldWeight is the weight for rhs

  // A typedef for the actual internal assignment vector type
  //typedef std::vector<AssignType> std::vector<AssignType>;

  DiscreteTable<AssignType>& lhs(*newPtr);
  std::map<std::vector<AssignType>, FProb>& mapP(lhs.sparseProbs);         // convenient alias
  double& defP(lhs.defProb);                            // convenient alias

  // Ensure rhs is also a DiscreteTable
  const DiscreteTable<AssignType>* dwnPtr = dynamic_cast<const DiscreteTable<AssignType>*>(oldPtr);
  PRLITE_ASSERT(dwnPtr, "DiscreteTable::distance : rhs is not a DiscreteTable");
  const DiscreteTable<AssignType>& rhs(*dwnPtr);
  const std::map<std::vector<AssignType>, FProb>& mapQ(rhs.sparseProbs); // convenient alias
  const double defQ(rhs.defProb);                      // convenient alias

  // Ensures the number of parameters of the two DiscreteTable potential functions are equal
  PRLITE_ASSERT(lhs.maxSz == rhs.maxSz,
         "DiscreteTable::distance The maxSz of the two distributions differ");
  double dist = 0.0; // dist1 = 0.0, dist2 = 0.0;

  // to calc the number of defaults occurring in neither of the maps
  double oldWeightCompl = 1.0-oldWeight;
  double weightedDefP = oldWeightCompl*defP;
  double weightedDefQ = oldWeight*defQ;
  double newDef = defP > 0.0 ? weightedDefP + weightedDefQ : defP;

  const double logDefP = defP > 0.0 ? log(defP) : -1000.0;
  const double defPLogDefP = defP > 0.0 ? defP*logDefP : 0.0;

  const double logDefQ = defQ > 0.0 ? log(defQ) : -1000.0;
  const double defQLogDefQ = defQ > 0.0 ? defQ*logDefQ : 0.0;

  const double defPLogDefQ = defP > 0.0 ? defP*logDefQ : 0.0;
  const double defQLogDefP = defQ > 0.0 ? defQ*logDefP : 0.0;

  // All the variables we'll need for calculating the distance

  size_t nDefsP = lhs.maxSz-mapP.size();
  size_t nDefsQ = rhs.maxSz-mapQ.size();

  // we first calculate the distance of the rhs to the unsmoothed lhs,
  // and then we smooth the lhs.
  auto iP = mapP.begin();
  auto iQ = mapQ.begin();
  while ( iP != mapP.end() && iQ != mapQ.end() ){

    if (iP->first < iQ->first) {      // only in mapP
      nDefsQ--;

      double& P(iP->second.prob);
      double logP = P > 0.0 ? log(P) : -1000.0;
      if (P > 0.0) {dist += P*(logP -logDefQ);} // if
      if (defQ > 0.0) {dist += defQLogDefQ -defQ*logP;} // if

      if (P > 0.0) { // we dont smooth probs that became zero
        P = prodUnderflowProtected(P, oldWeightCompl);
        P += weightedDefQ;
      } // if

      iP++;
    } // if

    else if (iP->first > iQ->first) { // only in mapQ
      nDefsP--;

      double Q(iQ->second.prob);
      double logQ = Q > 0.0 ? log(Q) : -1000.0;

      if (Q > 0.0) {dist += Q*(logQ -logDefP);} // if
      if (defP > 0.0) {dist += defPLogDefP -defP*logQ;} // if

      if (defP > 0.0) { // we dont smooth probs that became zero
        double Q(iQ->second.prob);
        double weightedQ = prodUnderflowProtected(oldWeight, Q);
        double prob = weightedDefP + weightedQ;
        if (std::abs(prob-newDef) > lhs.margin) {
          mapP[iQ->first] = prob;
        } // if
      } // if

      iQ++;
    } // else if

    else {                            //common to both maps
      double& P(iP->second.prob);
      double logP = P > 0.0 ? log(P) : -1000.0;
      double Q(iQ->second.prob);
      double logQ = Q > 0.0 ? log(Q) : -1000.0;

      if (P > 0.0) {dist += P*(logP -logQ);} // if
      if (Q > 0.0) {dist += Q*(logQ -logP);} // if

      if (P > 0.0) { // we dont smooth probs that became zero
        double Q(iQ->second.prob);
        P = prodUnderflowProtected(P, oldWeightCompl);
        P += oldWeight*Q;
      } // if

      iP++; iQ++;
    } // else

  } // while

  while ( iP != mapP.end() ){         // some extra spec's in mapP
    nDefsQ--;

    double& P(iP->second.prob);
    double logP = P > 0.0 ? log(P) : -1000.0;

    if (P > 0.0) {dist += P*(logP -logDefQ);} // if
    if (defQ > 0.0) {dist += defQLogDefQ -defQ*logP;} // if

    if (P > 0.0) { // we dont smooth probs that became zero
      P = prodUnderflowProtected(P, oldWeightCompl);
      P += weightedDefQ;
    } // if

    iP++;
  } // while

  while ( iQ != mapQ.end() ){         // some extra spec's in mapQ
    nDefsP--;

    double Q(iQ->second.prob);
    double logQ = Q > 0.0 ? log(Q) : -1000.0;

    if (Q > 0.0) {dist += Q*(logQ -logDefP);} // if
    if (defP > 0.0) {dist += defPLogDefP -defP*logQ;} // if

    if (defP > 0.0) { // we dont smooth probs that became zero
      double Q(iQ->second.prob);
      double weightedQ = prodUnderflowProtected(oldWeight, Q);
      double prob = weightedDefP + weightedQ;
      if (std::abs(prob-newDef) > lhs.margin) {
        mapP[iQ->first] = prob;
      } // if
    } // if

    iQ++;
  } // while

  // and now all the cases simultaneously absent from both maps
  PRLITE_ASSERT(
    nDefsP == nDefsQ,
    "DiscreteTable::distance The number of cases in neither map should match, no matter from which side you look at it.");

  if (nDefsP) {

    dist += nDefsP*(defPLogDefP - defPLogDefQ + defQLogDefQ - defQLogDefP);

  } // if

  if (defP > 0.0) { // we dont smooth probs that became zero
    defP = weightedDefP + weightedDefQ;
  } // if

  return dist/2.0;

} // inplaceProcess

// DEPRECATED stuff
