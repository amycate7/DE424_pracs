/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    April 2018
          PURPOSE: DirichletSet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers

// emdw headers
#include "mathfuncs.hpp"
#include "sortindices.hpp"
#include "conditionalpolya.hpp"

/******************************************************************************/
/***************************** DirichletSet ***********************************/
/******************************************************************************/


//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
// mostly boiler plate
template <typename CondType, typename DirType>
void DirichletSet<CondType, DirType>::supplementWithDefaultOperators() {


  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new DirichletSet_Marginalize<CondType, DirType>() );
  } //if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new DirichletSet_InplaceNormalize<CondType, DirType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new DirichletSet_Normalize<CondType, DirType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new DirichletSet_InplaceAbsorb<CondType, DirType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new DirichletSet_Absorb<CondType, DirType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new DirichletSet_InplaceCancel<CondType, DirType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new DirichletSet_Cancel<CondType, DirType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new DirichletSet_ObserveAndReduce<CondType, DirType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper =
      uniqptr<InplaceDamper>( new DirichletSet_InplaceDamping<CondType, DirType>() );
  } // if

} // supplementWithDefaultOperators

// default constructor
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>::DirichletSet(
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper)
    : marginalizer(theMarginalizer),
      inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper) { //An empty dirichlet is vacuous until values are set

  supplementWithDefaultOperators();

} // default constructor


// class specific constructor
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>::DirichletSet(
  emdw::RVIdType theVarId,
  const rcptr< std::vector<CondType> >& theDomain,
  const std::map< CondType, rcptr< Dirichlet<DirType> > >& theDirPtrs,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper)
: varId(theVarId),
  domain(theDomain),
  dirPtrs(),
  marginalizer(theMarginalizer),
  inplaceNormalizer(theInplaceNormalizer),
  normalizer(theNormalizer),
  inplaceAbsorber(theInplaceAbsorber),
  absorber(theAbsorber),
  inplaceCanceller(theInplaceCanceller),
  canceller(theCanceller),
  observeAndReducer(theObserveAndReducer),
  inplaceDamper(theInplaceDamper) {

  supplementWithDefaultOperators();

  for (auto el : theDirPtrs) {
    dirPtrs[el.first] = uniqptr< Dirichlet<DirType> > ( el.second->copy() );
  } // for

} // class specific constructor

// cpy ctor
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>::DirichletSet(const DirichletSet<CondType, DirType>& p)
    : Factor(p),
      varId(p.varId),
      domain(p.domain),
      dirPtrs(),
      marginalizer(p.marginalizer),
      inplaceNormalizer(p.inplaceNormalizer),
      normalizer(p.normalizer),
      inplaceAbsorber(p.inplaceAbsorber),
      absorber(p.absorber),
      inplaceCanceller(p.inplaceCanceller),
      canceller(p.canceller),
      observeAndReducer(p.observeAndReducer),
      inplaceDamper(p.inplaceDamper) {

  for (auto el : p.dirPtrs) {
    dirPtrs[el.first] = uniqptr< Dirichlet<DirType> > ( el.second->copy() );
  } // for

} // cpy ctor

template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>::DirichletSet(DirichletSet<CondType, DirType>&& p)
    : Factor(p),
      varId(p.varId),
      dirPtrs(p.dirPtrs),  // shallow copy
      marginalizer(p.marginalizer),
      inplaceNormalizer(p.inplaceNormalizer),
      normalizer(p.normalizer),
      inplaceAbsorber(p.inplaceAbsorber),
      absorber(p.absorber),
      inplaceCanceller(p.inplaceCanceller),
      canceller(p.canceller),
      observeAndReducer(p.observeAndReducer),
      inplaceDamper(p.inplaceDamper){
} // move ctor

// destructor
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>::~DirichletSet() {
} // destructor

//========================= Operators and Conversions =========================

// operator ==
template <typename CondType, typename DirType>
bool DirichletSet<CondType, DirType>::operator==(const DirichletSet<CondType, DirType>& p) const {

  if (varId != p.varId) {return false;} // if

  // check domain
  if (*domain != *p.domain) {return false;} // if

  // Test for the same number of Dirichlets
  if ( dirPtrs.size() != p.dirPtrs.size() ) {
    return false;
  } // if

  auto itr1 = dirPtrs.begin();
  auto itr2 = p.dirPtrs.begin();
  while ( itr1 != dirPtrs.end() ) {

    // same indices
    if (itr1->first != itr2->first) {
      return false;
    } // if
    // same Dirichlets
    if (*itr1->second != *itr2->second) {
      return false;
    } // if

    itr1++; itr2++;
  } // while
  return true;
} // operator ==

// operator !=
template <typename CondType, typename DirType>
bool DirichletSet<CondType, DirType>::operator!=(const DirichletSet<CondType, DirType>& p) const {
  return !operator==(p);
} // operator !=

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// classSpecificConfigure
template <typename CondType, typename DirType>
unsigned DirichletSet<CondType, DirType>::classSpecificConfigure(
  emdw::RVIdType theVarId,
  const rcptr< std::vector<CondType> >& theDomain,
  const std::map< CondType, rcptr< Dirichlet<DirType> > >& theDirPtrs,
  rcptr<Marginalizer> theMarginalizer,
  rcptr<InplaceNormalizer> theInplaceNormalizer,
  rcptr<Normalizer> theNormalizer,
  rcptr<InplaceAbsorber> theInplaceAbsorber,
  rcptr<Absorber> theAbsorber,
  rcptr<InplaceCanceller> theInplaceCanceller,
  rcptr<Canceller> theCanceller,
  rcptr<ObserveAndReducer> theObserveAndReducer,
  rcptr<InplaceDamper> theInplaceDamper) {

  this->~DirichletSet<CondType, DirType>(); // Destroy existing

  new(this) DirichletSet<CondType, DirType>(
    theVarId,
    theDomain,
    theDirPtrs,
    theMarginalizer,
    theInplaceNormalizer,
    theNormalizer,
    theInplaceAbsorber,
    theAbsorber,
    theInplaceCanceller,
    theCanceller,
    theObserveAndReducer,
    theInplaceDamper); // and do an inplace construction

  return true;
} // classSpecificConfigure

// txtWrite
template <typename CondType, typename DirType>
std::ostream& DirichletSet<CondType, DirType>::txtWrite(std::ostream& file) const {
  file << "DirichletSet_V0\n";
  file << "VarId: " << varId << std::endl;

  file << "\nDomain: " << *domain << std::endl;
  file << "nDirs: " << dirPtrs.size() << std::endl;
  for (auto aDir : dirPtrs) {
    file << aDir.first << " " << *aDir.second << std::endl;
  } // for
  return file;
} // txtWrite

// txtRead
template <typename CondType, typename DirType>
std::istream& DirichletSet<CondType, DirType>::txtRead(std::istream& file) {
  std::string dummy;

  file >> dummy; // strip the version tag

  emdw::RVIdType theVarId;
  file >> dummy >> theVarId;

  rcptr< std::vector<CondType> > theDomain;
  theDomain = uniqptr< std::vector<CondType> >(new std::vector<CondType>());
  file >>  dummy >> *theDomain;

  size_t nDirs;
  file >> dummy >> nDirs;

  std::map< CondType, rcptr< Dirichlet<DirType> > > theDirPtrs;
  for (size_t idx = 0; idx < nDirs; idx++) {
    CondType theKey;
    file >> theKey;
    rcptr< Dirichlet<DirType> > theVal =
      uniqptr< Dirichlet<DirType> > ( new Dirichlet<DirType>() );
    theVal->txtRead(file);
    theDirPtrs[theKey] = theVal;
  } // for

  classSpecificConfigure(theVarId, theDomain, theDirPtrs);

  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 1

// Boiler plate code
// returning new factor
template <typename CondType, typename DirType>
uniqptr<Factor> DirichletSet<CondType, DirType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void DirichletSet<CondType, DirType>::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> DirichletSet<CondType, DirType>::absorb(const Factor* rhsPtr,
                                                                FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void DirichletSet<CondType, DirType>::inplaceAbsorb(const Factor* rhsPtr,
                                                            FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> DirichletSet<CondType, DirType>::cancel(const Factor* rhsPtr,
                                                                FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void DirichletSet<CondType, DirType>::inplaceCancel(const Factor* rhsPtr,
                                                            FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> DirichletSet<CondType, DirType>::marginalize(const emdw::RVIds& variablesToKeep,
                                                                     bool presorted,
                                                                     const Factor* peekAheadPtr,
                                                                     FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep, presorted, peekAheadPtr) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep, presorted, peekAheadPtr) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> DirichletSet<CondType, DirType>::observeAndReduce(const emdw::RVIds& variables,
                                                                          const emdw::RVVals& assignedVals,
                                                                          bool presorted,
                                                                          FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals, presorted) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals, presorted) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
template <typename CondType, typename DirType>
double DirichletSet<CondType, DirType>::inplaceDampen(const Factor* rhsPtr,
                                                              double df,
                                                              FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

// copy
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>* DirichletSet<CondType, DirType>::copy(
  const emdw::RVIds& newVars, bool) const {

  if ( newVars.size() ) {

    std::map< CondType, rcptr< Dirichlet<DirType> > > theDirPtrs;
    for (auto itr = dirPtrs.begin(); itr != dirPtrs.end(); itr++) {
      theDirPtrs[itr->first] = uniqptr< Dirichlet<DirType> > ( itr->second->copy() );
    } // for

    return new DirichletSet<CondType, DirType>(
      newVars[0],
      domain,
      theDirPtrs,
      marginalizer,
      inplaceNormalizer,
      normalizer,
      inplaceAbsorber,
      absorber,
      inplaceCanceller,
      canceller,
      observeAndReducer,
      inplaceDamper);
  } // if

  return new DirichletSet<CondType, DirType>(*this);
} // copy

// vacuousCopy
template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>*
DirichletSet<CondType, DirType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool) const {

  PRLITE_ASSERT(selVars.size()<2, "DirichletSet::vacuousCopy : can not select more than one variable");
  if ( selVars.size() ) {
    PRLITE_ASSERT(selVars[0] == varId, "DirichletSet::vacuousCopy : Can only select variable "
           << varId);
  }// if

  DirichletSet<CondType, DirType>* fPtr = copy();
  for (auto& aDir : fPtr->dirPtrs) {
    aDir.second = uniqptr< Dirichlet<DirType> > ( aDir.second->vacuousCopy() );
  } // for
  return fPtr;
} // vacuousCopy

// isEqual mostly boiler plate code
template <typename CondType, typename DirType>
bool DirichletSet<CondType, DirType>::isEqual(const Factor* rhsPtr) const {
  const DirichletSet<CondType, DirType>* dwnPtr = dynamic_cast<const DirichletSet<CondType, DirType>*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

template <typename CondType, typename DirType>
double
DirichletSet<CondType, DirType>::distanceFromVacuous() const {
  double dist = 0.0;
  for (const auto& dirPtr : dirPtrs) {
    dist += dirPtr.second->distanceFromVacuous();
  } // for
  return dist;
} // distanceFromVacuous

// noOfVars
template <typename CondType, typename DirType>
unsigned DirichletSet<CondType, DirType>::noOfVars() const {
  return 1u;
} // noOfVars

// getVars
template <typename CondType, typename DirType>
emdw::RVIds DirichletSet<CondType, DirType>::getVars() const {
  return emdw::RVIds{varId};
} // getVars

// getVar
template <typename CondType, typename DirType>
emdw::RVIdType DirichletSet<CondType, DirType>::getVar(unsigned varNo) const {
  PRLITE_ASSERT(varNo ==0, "Only one id in a DirichletSet");
  return varId;
} // getVar

//====================== Other Virtual Member Functions =======================

// distance
template <typename CondType, typename DirType>
double DirichletSet<CondType, DirType>::distance(const Factor* rhsPtr) const {
  const DirichletSet<CondType, DirType>* dwnPtr = dynamic_cast<const DirichletSet<CondType, DirType>*>(rhsPtr);
  // Ensure both potentials are DirichletSet potentials
  PRLITE_ASSERT(dwnPtr, "DirichletSet::distance : rhs is not a DirichletSet");
  // Ensures the number of parameters of the two DirichletSet potential functions are equal
  PRLITE_ASSERT(dirPtrs.size() == dwnPtr->dirPtrs.size(), "The number of Dirichlets differs");

  double dist = 0.0;
  auto itr1 = dirPtrs.begin(); auto itr2 = dwnPtr->dirPtrs.begin();
  while ( itr1 != dirPtrs.end() ) {
    dist += itr1->second->distance( itr2->second.get() );
    itr1++; itr2++;
  } // while

  return dist;
} // distance

//======================= Non-virtual Member Functions ========================

template <typename CondType, typename DirType>
void DirichletSet<CondType, DirType>::copyExtraSettingsFrom(const DirichletSet<CondType, DirType>& d) {
  marginalizer = d.marginalizer;
  inplaceNormalizer = d.inplaceNormalizer;
  normalizer = d.normalizer;
  inplaceAbsorber = d.inplaceAbsorber;
  absorber = d.absorber;
  inplaceCanceller = d.inplaceCanceller;
  canceller = d.canceller;
  observeAndReducer = d.observeAndReducer;
  inplaceDamper = d.inplaceDamper;
} // copyExtraSettings

template <typename CondType, typename DirType>
uniqptr< Dirichlet<DirType> > DirichletSet<CondType, DirType>::copyDirichlet(
  const CondType& dirIdx) const {
  auto itr = dirPtrs.find(dirIdx);
  PRLITE_ASSERT(itr != dirPtrs.end(), "can not return a non-existent Dirichlet");
  return uniqptr< Dirichlet<DirType> > ( itr->second->copy() );
} // copyDirichlet

template <typename CondType, typename DirType>
rcptr< Dirichlet<DirType> > DirichletSet<CondType, DirType>::aliasDirichlet(
  const CondType& dirIdx) {
  auto itr = dirPtrs.find(dirIdx);
  PRLITE_ASSERT(itr != dirPtrs.end(), "can not return a non-existent Dirichlet");
  return itr->second;
} //aliasDirichlet

template <typename CondType, typename DirType>
emdw::RVVals DirichletSet<CondType, DirType>::sampleFrom(const CondType dirIdx) const {
  auto itr = dirPtrs.find(dirIdx);
  PRLITE_ASSERT(itr != dirPtrs.end(), "Could not find a Dirichlet indexed by " << dirIdx);
  return itr->second->sample(nullptr);
} // sampleFrom

template <typename CondType, typename DirType>
size_t DirichletSet<CondType, DirType>::noOfDirichlets() const {
  return dirPtrs.size();
} // noOfDirichlets

template <typename CondType, typename DirType>
size_t DirichletSet<CondType, DirType>::noOfCategories() const {
  return dirPtrs[0]->noOfCategories();
} // noOfCategories

template <typename CondType, typename DirType>
void
DirichletSet<CondType, DirType>::addToCounts(
  const CondType& setNr,
  const DirType& catNr,
  double cnt) {
  dirPtrs[setNr]->addToCounts(catNr,cnt);
} // addToCounts

template <typename CondType, typename DirType>
DirichletSet<CondType, DirType>*
DirichletSet<CondType, DirType>::sliceAt(const DirType& theObservedWord) const {

  std::map< CondType, rcptr< Dirichlet<DirType> > > theDirPtrs;
  for (auto condOn : *domain) {
    auto itr = dirPtrs.find(condOn);
    if ( itr != dirPtrs.end() ) {
      theDirPtrs.insert (std::make_pair(condOn, itr->second->sliceAt(theObservedWord) ) );
    } // if
    else {
      PRLITE_ASSERT(false, "DirichletSet::sliceAt : unknown condition " << condOn);
    } // else
  } // for

  return new DirichletSet<CondType, DirType>(
    varId,
    domain,
    theDirPtrs,
    marginalizer,
    inplaceNormalizer,
    normalizer,
    inplaceAbsorber,
    absorber,
    inplaceCanceller,
    canceller,
    observeAndReducer,
    inplaceDamper);

} // sliceAt

/*****************************************************************************/
/***************************** Family 1 **************************************/
/***************************** Normalize *************************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_InplaceNormalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_InplaceNormalize<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename CondType, typename DirType>
void DirichletSet_InplaceNormalize<CondType, DirType>::inplaceProcess(DirichletSet<CondType, DirType>* lhsPtr) {
  for (auto& aDir : lhsPtr->dirPtrs) {
    aDir.second->Factor::inplaceNormalize();
  } // for

} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_Normalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_Normalize<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boiler plate code
template <typename CondType, typename DirType>
Factor* DirichletSet_Normalize<CondType, DirType>::process(const DirichletSet<CondType, DirType>* lhsPtr) {
  // make a copy of the factor
  DirichletSet<CondType, DirType>* fPtr = new DirichletSet<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  DirichletSet_InplaceNormalize<CondType, DirType> ipNorm;
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
/***************************** Family 2 **************************************/
/***************************** Absorb ****************************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_InplaceAbsorb<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_InplaceAbsorb<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One DirichletSet absorbs another - this is implemented by adding the parameters
// of the given DirichletSets
template <typename CondType, typename DirType>
void DirichletSet_InplaceAbsorb<CondType, DirType>::inplaceProcess(DirichletSet<CondType, DirType>* lhsPtr, const Factor* rhsFPtr) {

  //std::cout << "Multiply: " << *lhsPtr << std::endl;

  const DirichletSet<CondType, DirType>* rhsPtr1 = dynamic_cast<const DirichletSet<CondType, DirType>*>(rhsFPtr);
  if (rhsPtr1) { // ok, we've got a DirichletSet
    //std::cout << "\nWith: " << *rhsFPtr  << std::endl;
    PRLITE_ASSERT(lhsPtr->getVar(0) == rhsPtr1->getVar(0),
           "DirichletSet_InplaceAbsorb : Both DirichletSets must share the same random variable");

    auto lItr = lhsPtr->dirPtrs.begin();
    auto rItr = rhsPtr1->dirPtrs.begin();
    while ( rItr != rhsPtr1->dirPtrs.end() ) {

      lItr->second->Factor::inplaceAbsorb(rItr->second);
      lItr++; rItr++;
    } // while

    //std::cout << "\nGives: " << *lhsPtr << std::endl;

    return;
  } // if

  PRLITE_ASSERT(false,
         "DirichletSet_InplaceAbsorb : rhs is an invalid type: " << typeid(*rhsFPtr).name() );

} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_Absorb<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_Absorb<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boiler plate code
template <typename CondType, typename DirType>
Factor* DirichletSet_Absorb<CondType, DirType>::process(const DirichletSet<CondType, DirType>* lhsPtr, const Factor* rhsFPtr) {
  // make a copy of the factor
  DirichletSet<CondType, DirType>* fPtr = new DirichletSet<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  DirichletSet_InplaceAbsorb<CondType, DirType> ipAbsorb;
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
/***************************** Family 2 **************************************/
/***************************** Cancel ****************************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_InplaceCancel<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_InplaceCancel<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One DirichletSet cancels another - this is implemented by subtracting the parameter
// of the one potential from that of the other
template <typename CondType, typename DirType>
void DirichletSet_InplaceCancel<CondType, DirType>::inplaceProcess(DirichletSet<CondType, DirType>* lhsPtr, const Factor* rhsFPtr) {
  const DirichletSet<CondType, DirType>* rhsPtr1 = dynamic_cast<const DirichletSet<CondType, DirType>*>(rhsFPtr);
  if (rhsPtr1) { // ok, we've got a DirichletSet
    PRLITE_ASSERT(lhsPtr->getVar(0) == rhsPtr1->getVar(0),
           "DirichletSet_InplaceCancel : Both DirichletSets must share the same random variable");

    auto lItr = lhsPtr->dirPtrs.begin();
    auto rItr = rhsPtr1->dirPtrs.begin();
    while ( rItr != rhsPtr1->dirPtrs.end() ) {

      lItr->second->Factor::inplaceCancel(rItr->second);
      lItr++; rItr++;
    } // while

    return;
  } // if

  PRLITE_ASSERT(false,
         "DirichletSet_InplaceCancel : rhs is an invalid type: " << typeid(*rhsFPtr).name() );
} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_Cancel<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_Cancel<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boiler plate code
template <typename CondType, typename DirType>
Factor* DirichletSet_Cancel<CondType, DirType>::process(const DirichletSet<CondType, DirType>* lhsPtr, const Factor* rhsFPtr) {
  // make a copy of the factor
  DirichletSet<CondType, DirType>* fPtr = new DirichletSet<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  DirichletSet_InplaceCancel<CondType, DirType> ipCancel;
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
/********************************** Marginalize ******************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_Marginalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_Marginalize<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Given DirichletSet is already marginalized in its parametric form
// This function therefore merely returns a copy of the given potential
template <typename CondType, typename DirType>
Factor* DirichletSet_Marginalize<CondType, DirType>::process(
  const DirichletSet<CondType, DirType>* lhsPtr,
  const emdw::RVIds& toKeep,
  bool presorted,
  const Factor* peekAheadPtr) {

  PRLITE_ASSERT(toKeep[0] == lhsPtr->getVar(0),
         "DirichletSet_Marginalize bad marginalization variable");

  if (peekAheadPtr) {

    const ConditionalPolya<CondType, DirType>* dwnPtr =
        dynamic_cast<const ConditionalPolya<CondType, DirType>*> (peekAheadPtr);
    PRLITE_ASSERT(dwnPtr, "DirichletSet_Marginalize : Don't know what to do with a "
           << *peekAheadPtr);

    ValidValue<DirType> theObservedX = dwnPtr->getObservedX();
    if (theObservedX.valid) {
      return lhsPtr->sliceAt(theObservedX.value);
    } // if

  } // if

  return lhsPtr->copy();

} // process

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& DirichletSet_InplaceDamping<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("DirichletSet_InplaceDamping<"+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and old
// parameters - this linear combination is determined by the df parameter
template <typename CondType, typename DirType>
double DirichletSet_InplaceDamping<CondType, DirType>::inplaceProcess(
  DirichletSet<CondType, DirType>* lhsPtr,
  const Factor* rhsFPtr,
  double df) { //Remember df is the weight for rhs

  const DirichletSet<CondType, DirType>* rhsPtr =
    dynamic_cast<const DirichletSet<CondType, DirType>*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr,
         "DirichletSet_InplaceDamping : rhs not a valid DirichletSet : " << typeid(*rhsPtr).name() );

  // Create references to DirichletSet instances
  const DirichletSet<CondType, DirType>& rhs(*rhsPtr);
  DirichletSet<CondType, DirType>& lhs(*lhsPtr);

  PRLITE_ASSERT(lhs.getVar(0) == rhs.getVar(0),
         "DirichletSet_InplaceDamping : Both DirichletSets must share the same random variable");

  double dist = 0.0;
  auto itr1 = lhs.dirPtrs.begin(); auto itr2 = rhs.dirPtrs.begin();
  while ( itr1 != lhs.dirPtrs.end() ) {

    dist += itr1->second->Factor::inplaceDampen(itr2->second, df);
    itr1++; itr2++;
  } // while

  return dist;
} // inplaceProcess
