/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2018
          PURPOSE: Polya PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "discretetable.hpp"
#include "vecset.hpp"

// standard headers
#include <new>  // placement new
#include <math.h>
#include <typeinfo>
#include <limits>  // numeric_limits


/******************************************************************************/
/****************************** Polya ***********************************/
/******************************************************************************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
template <typename AssignType>
void Polya<AssignType>::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new Polya_Marginalize<AssignType>() );
  } //if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new Polya_InplaceNormalize<AssignType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new Polya_Normalize<AssignType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new Polya_InplaceAbsorb<AssignType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new Polya_Absorb<AssignType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new Polya_InplaceCancel<AssignType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new Polya_Cancel<AssignType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new Polya_ObserveAndReduce<AssignType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper =
      uniqptr<InplaceDamper>( new Polya_InplaceDamping<AssignType>() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new Polya_Sampler<AssignType>() );
  } // if

} // supplementWithDefaultOperators

// class specific constructor
template <typename AssignType>
Polya<AssignType>::Polya(
  emdw::RVIdType theDirVarId,
  emdw::RVIdType theCatVarId,
  rcptr< std::vector<AssignType> > theDomain,
  rcptr< Dirichlet<AssignType> > theDirPtr,
  rcptr< DiscreteTable<AssignType> > theCatPtr,
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
    : dirVarId(theDirVarId),
      catVarId(theCatVarId),
      domain(theDomain),
      dirPtr( std::move(theDirPtr) ),
      catPtr( std::move(theCatPtr) ),
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

  PRLITE_ASSERT(dirVarId < catVarId, "The Dirichlet ID must be smallest, then the Categorical ID");

  if (dirPtr) {
    PRLITE_ASSERT(
      dirPtr->getVar(0) == dirVarId,
      "Dirichlet variable should be " << dirVarId
      << ", not " << dirPtr->getVar(0) );
  } // if
  else {
    dirPtr = std::make_unique< Dirichlet<AssignType> >(dirVarId, domain);
  } // else

  if (catPtr) {
    PRLITE_ASSERT(
      catPtr->getVar(0) == catVarId,
      "DiscreteTable variable should be " << catVarId
      << ", not " << catPtr->getVar(0) );
  } // if

#ifdef CHK_IDX
  // we need something here to check that if dirPtr is set, it matches
  // the number of classes from domain
#endif // CHK_IDX

} // class specific constructor

// copy constructor is explicit since we need separate/deep dirPtr,
// catPtr (for instance with multiple links to other clusters)
template <typename AssignType>
Polya<AssignType>::Polya(const Polya<AssignType>& tcp)
    : Factor(tcp),
      dirVarId(tcp.dirVarId),
      catVarId(tcp.catVarId),
      domain(tcp.domain),
      dirPtr(0),
      catPtr(0),
      marginalizer(tcp.marginalizer),
      inplaceNormalizer(tcp.inplaceNormalizer),
      normalizer(tcp.normalizer),
      inplaceAbsorber(tcp.inplaceAbsorber),
      absorber(tcp.absorber),
      inplaceCanceller(tcp.inplaceCanceller),
      canceller(tcp.canceller),
      observeAndReducer(tcp.observeAndReducer),
      inplaceDamper(tcp.inplaceDamper),
      sampler(tcp.sampler) {
  if (tcp.dirPtr) { dirPtr = rcptr< Dirichlet<AssignType> >( tcp.dirPtr->copy() ); } // if
  if (tcp.catPtr) { catPtr = rcptr< DiscreteTable<AssignType> >( tcp.catPtr->copy() ); } // if
} // copy ctor

// move constructor
template <typename AssignType>
Polya<AssignType>::Polya(Polya<AssignType>&& tcp)
    : Factor(tcp),
      dirVarId(tcp.dirVarId),
      catVarId(tcp.catVarId),
      domain(tcp.domain),
      dirPtr(tcp.dirPtr),
      catPtr(tcp.catPtr),
      marginalizer(tcp.marginalizer),
      inplaceNormalizer(tcp.inplaceNormalizer),
      normalizer(tcp.normalizer),
      inplaceAbsorber(tcp.inplaceAbsorber),
      absorber(tcp.absorber),
      inplaceCanceller(tcp.inplaceCanceller),
      canceller(tcp.canceller),
      observeAndReducer(tcp.observeAndReducer),
      inplaceDamper(tcp.inplaceDamper),
      sampler(tcp.sampler) {
} // move ctor

// default constructor
template <typename AssignType>
Polya<AssignType>::Polya(
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
    : inplaceNormalizer(theInplaceNormalizer),
      normalizer(theNormalizer),
      inplaceAbsorber(theInplaceAbsorber),
      absorber(theAbsorber),
      inplaceCanceller(theInplaceCanceller),
      canceller(theCanceller),
      observeAndReducer(theObserveAndReducer),
      inplaceDamper(theInplaceDamper),
      sampler(theSampler) {

  supplementWithDefaultOperators();

} // default constructor

template <typename AssignType>
Polya<AssignType>::~Polya() {
} // destructor

//========================= Operators and Conversions =========================

template <typename AssignType>
bool Polya<AssignType>::operator==(const Polya<AssignType>& p) const {

  if (dirVarId != p.dirVarId) {
    return false;
  } // if

  if (catVarId != p.catVarId) {
    return false;
  } // if

    if (*domain != *p.domain) {
      return false;
    } // if

  // funky logical xor. extra nots to force them both to boolean
  if ( !dirPtr.get() != !p.dirPtr.get()  ) {
    return false;
  } // if

  if (dirPtr and *dirPtr != *p.dirPtr) {
    return false;
  } // if

  // funky logical xor. extra nots to force them both to boolean
  if ( !catPtr.get() != !p.catPtr.get() ) {
    return false;
  } // if

  if (catPtr and *catPtr != *p.catPtr) {
    return false;
  } // if

  return true;

} // operator ==

template <typename AssignType>
bool Polya<AssignType>::operator!=(const Polya<AssignType>& p) const {
  return !operator==(p);
} // operator !=

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================


// class specific
template <typename AssignType>
unsigned Polya<AssignType>::classSpecificConfigure(
  emdw::RVIdType theDirVarId,
  emdw::RVIdType theCatVarId,
  rcptr< std::vector<AssignType> > theDomain,
  rcptr< Dirichlet<AssignType> > theDirPtr,
  rcptr< DiscreteTable<AssignType> > theCatPtr,
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

  // Destroy existing
  this->~Polya<AssignType>();

  // and do an inplace construction
  new(this) Polya<AssignType>(
    theDirVarId,
    theCatVarId,
    theDomain,
    theDirPtr,
    theCatPtr,
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

// !!! this still needs fixing
template <typename AssignType>
std::ostream& Polya<AssignType>::txtWrite(std::ostream& file) const {
  file << "Polya_" << typeid(AssignType).name() << "_V0\n";

  file << "DirichletVar: " << dirVarId << " DiscreteTableVar: " << catVarId << std::endl;

  file << "\nDomain: "  << *domain << std::endl;

  file << "Dirichlet:\n";
  if (dirPtr.get() ) { file << "1 " << *dirPtr << std::endl;}
  else {file << "0\n";}

  file << "DiscreteTable:\n";
  if (catPtr.get() ) { file << "1 " << *catPtr << std::endl;}
  else {file << "0\n";}

  // if (defaultMarginalizer) {file << defaultMarginalizer->isA() << std::endl; }
  // if (defaultInplaceNormalizer) {file << defaultInplaceNormalizer->isA() << std::endl; }
  // if (defaultNormalizer) {file << defaultNormalizer->isA() << std::endl; }
  // if (defaultInplaceAbsorber) {file << defaultInplaceAbsorber->isA() << std::endl; }
  // if (defaultAbsorber) {file << defaultAbsorber->isA() << std::endl; }
  // if (defaultInplaceCanceller) {file << defaultInplaceCanceller->isA() << std::endl; }
  // if (defaultCanceller) {file << defaultCanceller->isA() << std::endl; }
  // if (defaultObserveAndReducer) {file << defaultObserveAndReducer->isA() << std::endl; }

  return file;
} // txtWrite


// !!! this still needs fixing
template <typename AssignType>
std::istream& Polya<AssignType>::txtRead(std::istream& file) {

  std::string dummy;
  file >> dummy;  // strip the version tag

  emdw::RVIdType theDirVarId;
  file >> dummy >> theDirVarId;

  emdw::RVIdType theCatVarId;
  file >> dummy >> theCatVarId;

  rcptr< std::vector<AssignType> > theDomain( new std::vector<AssignType>() );
  file >> dummy >> *theDomain; // to bypass the "Domain:" text (added by Rebecca 09112017)

  unsigned tmp;
  file >> dummy >> tmp;
  if (tmp != 0) {
    std::cerr << __FILE__ << __LINE__ << ": sorry, txtRead NIY!\n";
    exit(-1);
  } // if
  file >> dummy >> tmp;
  if (tmp != 0) {
    std::cerr << __FILE__ << __LINE__ << ": sorry, txtRead NIY!\n";
    exit(-1);
  } // if

  classSpecificConfigure(theDirVarId, theCatVarId, theDomain);
  return file;
} // txtRead

template <typename AssignType>
Polya<AssignType>* Polya<AssignType>::copy(
  const emdw::RVIds& newVars,
  bool) const {

  if ( newVars.size() ) {
    return new Polya<AssignType>(
      newVars[0],
      newVars[1],
      domain,
      dirPtr,
      catPtr,
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
  } // if

  return new Polya<AssignType>(*this);
} // copy

template <typename AssignType>
Polya<AssignType>*
Polya<AssignType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool) const {

  PRLITE_ASSERT(!selVars.size(), "Polya::vacuousCopy : variable selection NIY");
  Polya<AssignType>* fPtr = copy();

  if (dirPtr) {
    fPtr->dirPtr = uniqptr< Dirichlet<AssignType> >( dirPtr->vacuousCopy() );
  } // if

  if (catPtr) {
    fPtr->catPtr
      = uniqptr< DiscreteTable<AssignType> >( static_cast<DiscreteTable<AssignType>*> (catPtr->vacuousCopy() ) );
  } // if

  return fPtr;
} // vacuousCopy


//===================== Required Factor Operations ============================

// Boiler plate code

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
template <typename AssignType>
emdw::RVVals Polya<AssignType>::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// boiler plate code
// inplace
template <typename AssignType>
void Polya<AssignType>::inplaceNormalize(FactorOperator* procPtr) {
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

// boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Polya<AssignType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} //normalize

//------------------Family 2

// boiler plate code
// inplace
template <typename AssignType>
void Polya<AssignType>::inplaceAbsorb(const Factor* rhsPtr,
                                      FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} //inplaceAbsorb

// boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Polya<AssignType>::absorb(const Factor* rhsPtr,
                                          FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// boiler plate code
// inplace
template <typename AssignType>
void Polya<AssignType>::inplaceCancel(const Factor* rhsPtr,
                                      FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} //inplaceCancel

// boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Polya<AssignType>::cancel(const Factor* rhsPtr,
                                          FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

//------------------Family 3

// boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Polya<AssignType>::marginalize(const emdw::RVIds& variablesToKeep,
                                               bool presorted,
                                               const Factor* peekAheadPtr,
                                               FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep, presorted, peekAheadPtr) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep, presorted, peekAheadPtr) );
} //marginalize

//------------------Family 4: eg Observe

// boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Polya<AssignType>::observeAndReduce(const emdw::RVIds& variables,
                                                    const emdw::RVVals& assignedVals,
                                                    bool presorted,
                                                    FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals, presorted) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals, presorted) );
} //observeAndReduce

//------------------Family 5: eg inplaceDampen

// boiler plate code
// inplace returning a double
template <typename AssignType>
double Polya<AssignType>::inplaceDampen(const Factor* rhsPtr,
                                        double df,
                                        FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

// boiler plate code
template <typename AssignType>
bool Polya<AssignType>::isEqual(const Factor* rhsPtr) const {
    const Polya<AssignType>* dwnPtr = dynamic_cast<const Polya<AssignType>*>(rhsPtr);
    if (!dwnPtr) {return false;}
    return operator==(*dwnPtr);
} // isEqual

template <typename AssignType>
double
Polya<AssignType>::distanceFromVacuous() const {
  double dist = 0.0;
  if (dirPtr) {dist += dirPtr->distanceFromVacuous();} // if
  if (catPtr) {dist += catPtr->distanceFromVacuous();} // if
  return dist;
} // distanceFromVacuous

template <typename AssignType>
unsigned Polya<AssignType>::noOfVars() const {
  return 2;
} // noOfVars

template <typename AssignType>
emdw::RVIds Polya<AssignType>::getVars() const {

  emdw::RVIds theVars(2);
  theVars[0] = dirVarId;
  theVars[1] = catVarId;

  return theVars;
} // getVars

template <typename AssignType>
emdw::RVIdType Polya<AssignType>::getVar(unsigned varNo) const {
  if (varNo) {return catVarId;} // if
  else{return dirVarId;} // else
} // getVar

//====================== Other Virtual Member Functions =======================

template <typename AssignType>
double Polya<AssignType>::distance(const Factor* rhsPtr) const {

  const Polya<AssignType>* dwnPtr = dynamic_cast<const Polya<AssignType>*>(rhsPtr);
  PRLITE_ASSERT(dwnPtr, "Polya::distance : rhs is not a Polya");

  // will only do this from one side. could possibly miss some
  // non-zero rhs probs.

  double dist = 0.0;

  if (dirPtr) {
    if (dwnPtr->dirPtr) {dist += dirPtr->distance( dwnPtr->dirPtr.get() );}
    else {return std::numeric_limits<double>::infinity();}
  } // if

  if (catPtr) {
    if (dwnPtr->catPtr) {dist += catPtr->distance( dwnPtr->catPtr.get() );}
    else {return std::numeric_limits<double>::infinity();}
  } // if

  return dist;
} // distance

//======================= Non-virtual Member Functions ========================

template <typename AssignType>
uniqptr< Dirichlet<AssignType> > Polya<AssignType>::getDirichlet() const {
  if (dirPtr) {return uniqptr< Dirichlet<AssignType> >( dirPtr->copy() );} // if
  else {return uniqptr< Dirichlet<AssignType> >();} // else
} // getDirichlet

template <typename AssignType>
uniqptr< DiscreteTable<AssignType> > Polya<AssignType>::getDiscreteTable() const {
  if (catPtr) {return uniqptr< DiscreteTable<AssignType> >( catPtr->copy() );} // if
  else{return uniqptr< DiscreteTable<AssignType> >();} // else
} // getDiscreteTable

// template <typename AssignType>
// emdw::RVIdType Polya<AssignType>::getDirichletId() const {return dirVarId;} // getDirichletId

// template <typename AssignType>
//   emdw::RVIdType Polya<AssignType>::getCategoricalId() const {return catVarId;} // getCategoricalId

//================================== Friends ==================================

/// output
template<typename T>
std::ostream& operator<<( std::ostream& file,
                          const Polya<T>& n ){
  return n.txtWrite(file);
} // operator<<


/// input
template<typename T>
std::istream& operator>>( std::istream& file,
                          Polya<T>& n ){
  return n.txtRead(file);
} // operator>>

// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample ****************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& Polya_Sampler<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_Sampler");
  return CLASSNAME;
} // isA

// Samples a random vector from the distribution
template <typename AssignType>
emdw::RVVals Polya_Sampler<AssignType>::sample(const Polya<AssignType>* lhsPtr) {

  emdw::RVVals samples;
  PRLITE_ASSERT(lhsPtr->dirPtr, "Polya_Sampler:: Dirichlet prior uninstantiated");

  // sample the theta probabilities from the Dirichlet
  samples = lhsPtr->dirPtr->sample(nullptr);
  std::vector<double> thetas;
  for (auto el : samples) {
    thetas.push_back( double(el) );
  } // for

  // build the Categorical and sample from it.
  std::map<std::vector<AssignType>, FProb> sparseProbs;
  for (unsigned idx =0; idx < thetas.size(); idx++) {
    sparseProbs[{(*lhsPtr->domain)[idx]}] = thetas[idx];
  } // for
  rcptr<Factor> polyaCatPtr = uniqptr< DiscreteTable<AssignType> > (
    new DiscreteTable<AssignType>({1}, {lhsPtr->domain}, 0.0, sparseProbs) );
  //cout << "Categorical Polya: " << *polyaCatPtr << endl;
  samples.push_back(polyaCatPtr->sample()[0]);

  return samples;
} // sample

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Polya_InplaceNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_InplaceNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
void Polya_InplaceNormalize<AssignType>::inplaceProcess(Polya<AssignType>* lstPtr) {
  if (lstPtr->dirPtr) {lstPtr->dirPtr->inplaceNormalize(nullptr);} // if
  if (lstPtr->catPtr) {lstPtr->catPtr->inplaceNormalize(nullptr);} // if
} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Polya_Normalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_Normalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor* Polya_Normalize<AssignType>::process(const Polya<AssignType>* lhsPtr) {

  // make a copy of the factor
  Polya<AssignType>* fPtr = new Polya<AssignType>(*lhsPtr);

  // inplace normalize it
  Polya_InplaceNormalize<AssignType> ipNorm;
  ipNorm.inplaceProcess(fPtr);

  // and return the result
  return fPtr;
} // process

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Polya_InplaceAbsorb<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_InplaceAbsorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * This operator inplace absorbs another Dirichlet or DiscreteTable,
 * ie it does an inplace factor product.
 */
template <typename AssignType>
void Polya_InplaceAbsorb<AssignType>::inplaceProcess(Polya<AssignType>* lhsPtr,
                                                  const Factor* rhsFPtr){

  // check out what type of rhs Factor we've got and then handle it as needed

  // Multiplying with a Dirichlet
  const  Dirichlet<AssignType> * rhsPtr2 = dynamic_cast<const  Dirichlet<AssignType> *>(rhsFPtr);
  if (rhsPtr2) { // ok, we've got a Dirichlet
    PRLITE_ASSERT(lhsPtr->dirVarId == rhsPtr2->getVar(0),
           "Polya_InplaceAbsorb : Both Dirichlets must be defined over the same random variable ids");

    if (lhsPtr->dirPtr) {
      lhsPtr->dirPtr->inplaceAbsorb(rhsFPtr,nullptr);
    } // if
    else {
      lhsPtr->dirPtr = rcptr< Dirichlet<AssignType> >( rhsPtr2->copy() );
    } // else
    return;
  } // if

  // Multiplying with a DiscreteTable
  const DiscreteTable<AssignType>* rhsPtr3 = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  if (rhsPtr3) { // ok, we've got a DiscreteTable
    PRLITE_ASSERT(lhsPtr->catVarId == rhsPtr3->getVar(0),
           "Polya_InplaceAbsorb : The DiscreteTable must be defined over the same random variable ids as the Polya");
    if (lhsPtr->catPtr) {
      lhsPtr->catPtr->inplaceAbsorb(rhsFPtr,0);
    } // if
    else {
      lhsPtr->catPtr = rcptr< DiscreteTable<AssignType> >( static_cast<DiscreteTable<AssignType>*>( rhsPtr3->copy() ) );
    } // else
    return;
  } // if

  PRLITE_ASSERT(false, "Polya_InplaceAbsorb : rhs is an invalid type: " << typeid(*rhsFPtr).name() );

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Polya_Absorb<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_Absorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Polya_InplaceCancel<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_InplaceCancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * This operator inplace cancels either a Dirichlet or a DiscreteTable
 * message, ie it does an inplace factor division.
 */
template <typename AssignType>
void Polya_InplaceCancel<AssignType>::inplaceProcess(Polya<AssignType>* lhsPtr,
                                                  const Factor* rhsFPtr){

  // check out what type of rhs Factor we've got and then handle it as needed

  // Dividing by a Dirichlet
  const  Dirichlet<AssignType> * rhsPtr2 = dynamic_cast<const  Dirichlet<AssignType> *>(rhsFPtr);
  if (rhsPtr2) { // ok, we've got a Dirichlet
    PRLITE_ASSERT(lhsPtr->dirVarId == rhsPtr2->getVar(0),
           "Polya_InplaceCancel : Both Dirichlets must be defined over the same random variable ids");
    PRLITE_ASSERT(lhsPtr->dirPtr,
           "Polya_InplaceCancel : Can not cancel a rhs Dirichlet from a non-existing lhs");
    lhsPtr->dirPtr->inplaceCancel(rhsFPtr,nullptr);
    return;
  } // if

  // Dividing by a DiscreteTable
  const DiscreteTable<AssignType>* rhsPtr3 = dynamic_cast<const DiscreteTable<AssignType>*>(rhsFPtr);
  if (rhsPtr3) { // ok, we've got a DiscreteTable
    PRLITE_ASSERT(lhsPtr->catVarId == rhsPtr3->getVar(0),
           "Polya_InplaceCancel : Both DiscreteTables must be defined over the same random variable ids");
    PRLITE_ASSERT(lhsPtr->catPtr,
           "Polya_InplaceCancel : Can not cancel a rhs DiscreteTable from a non-existing lhs");
    lhsPtr->catPtr->inplaceCancel(rhsFPtr,nullptr);
    return;
  } // if

  PRLITE_ASSERT(false, "Polya_InplaceCancel : rhs is an invalid type: " << typeid(*rhsFPtr).name() );

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Polya_Cancel<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_Cancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

//!!!
/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Polya_Marginalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_Marginalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename AssignType>
Factor*
Polya_Marginalize<AssignType>::process(
  const Polya<AssignType>* lhsPtr,
  const emdw::RVIds& toKeep,
  bool presorted,
  const Factor* peekAheadPtr) {

  PRLITE_ASSERT(toKeep.size() > 0, "Polya_Marginalize : There are no variables to retain");

  // make sure that toKeep is sorted
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

  // ++++++++++++++++++++++++++++++
  // +++++ keep the Dirichlet +++++
  // ++++++++++++++++++++++++++++++

  if ( sortedKeep[0] == lhsPtr->dirVarId and sortedKeep.size() == 1 ) {

    // we have 4 cases to consider depending on whether or not dirPtr
    // and/or catPtr is instantiated

    if (lhsPtr->dirPtr) { // this should nearly always be the case

      // if no Categorical, return the Dirichlet
      if (!lhsPtr->catPtr) {
        return lhsPtr->dirPtr->copy();
      } // if

      // ok, so we have a Categorical and will have to update the alphas
      // ak'= ak + \frac{pk*ak}{\sum_j pj*aj}

      lhsPtr->catPtr->inplaceNormalize(nullptr);

      // to do this effectively we will use the double ladder
      // trick. To make things fast we will save the intermediate
      // results in two corresponding maps
      std::map< AssignType, double> theAj; // the aj terms
      std::map< AssignType, double> thePjAj; // the pj*aj terms
      double sumPjAj = 0.0;

      const std::map<AssignType, double>& mapL(lhsPtr->dirPtr->getSparseAlphas() );
      const std::map<std::vector<AssignType>, FProb>& mapR(lhsPtr->catPtr->getSparseProbs() );
      double defL = lhsPtr->dirPtr->getDefaultAlpha();
      double defR = lhsPtr->catPtr->getDefaultProb();

      auto iL = mapL.begin();
      auto iR = mapR.begin();
      while ( iL != mapL.end() && iR != mapR.end() ){

        if (iL->first < iR->first[0]) {      // mapL,defR
          theAj.insert( theAj.end(), std::make_pair(iL->first,iL->second) ); // hinted insert for speed
          double theirProd = iL->second*defR;
          thePjAj.insert(thePjAj.end(), std::make_pair(iL->first,theirProd) );
          sumPjAj += theirProd;
          iL++;
        } // if

        else if (iL->first > iR->first[0]) { // defL,mapR
          theAj.insert( theAj.end(), std::make_pair(iR->first[0],defL) );
          double theirProd = defL*iR->second.prob;
          thePjAj.insert(thePjAj.end(), std::make_pair(iR->first[0],theirProd) );
          sumPjAj += theirProd;
          iR++;
        } // else if

        else {                            // mapL,mapR
          theAj.insert( theAj.end(), std::make_pair(iL->first,iL->second) );
          double theirProd = iL->second*iR->second.prob;
          thePjAj.insert(thePjAj.end(), std::make_pair(iL->first,theirProd) );
          sumPjAj += theirProd;
          iL++; iR++;
        } // else

      } // while

      while ( iL != mapL.end() ){         // some extra mapL,defR
        theAj.insert( theAj.end(), std::make_pair(iL->first,iL->second) );
        double theirProd = iL->second*defR;
        thePjAj.insert(thePjAj.end(), std::make_pair(iL->first,theirProd) );
        sumPjAj += theirProd;
        iL++;
      } // while

      while ( iR != mapR.end() ){         // some extra defL,mapR
        theAj.insert( theAj.end(), std::make_pair(iR->first[0],defL) );
        double theirProd = defL*iR->second.prob;
        thePjAj.insert(thePjAj.end(), std::make_pair(iR->first[0],theirProd) );
        sumPjAj += theirProd;
        iR++;
      } // while

      // theAj.size() is the combined updated size. Subtracting that from
      // noOfCategories gives the remaining number of unseen defaults.
      double defProd = defL*defR;
      sumPjAj += defProd * ( lhsPtr->dirPtr->noOfCategories()-theAj.size() );
      if (sumPjAj > 0.0) {
        defProd /= sumPjAj;
      } // if
      else {
        defProd = 0.0;
      } // else

      auto itr1 = theAj.begin();
      auto itr2 = thePjAj.begin();
      while ( itr1 != theAj.end() ){
        if (sumPjAj > 0.0) {
          itr1->second += itr2->second/sumPjAj;
        } // if
        else {
          itr1->second = 0.0;
        } // else
        itr1++; itr2++;
      } //while

      Dirichlet<AssignType> * retPtr =
        new Dirichlet<AssignType>(lhsPtr->dirVarId, lhsPtr->domain, defL + defProd, theAj);
      retPtr->copyExtraSettingsFrom(*lhsPtr->dirPtr);
      return retPtr;

    } // if has a Dirichlet

    else { // no Dirichlet: equivalent to all alphas = 1

      if (lhsPtr->catPtr) { // has only a Categorical
        lhsPtr->catPtr->inplaceNormalize(nullptr);
        const std::map<std::vector<AssignType>, FProb>& mapR( lhsPtr->catPtr->getSparseProbs() );
        double defR = lhsPtr->catPtr->getDefaultProb();
        std::map<AssignType, double> theAj; // the aj terms

        for (auto iR = mapR.begin(); iR != mapR.end(); iR++) {
          theAj.insert(theAj.end(), std::make_pair(iR->first[0], iR->second.prob+1.0) );  // As if we received an uninformative Dirichlet
        } // for
        return new Dirichlet<AssignType>(lhsPtr->dirVarId, lhsPtr->domain, defR+1.0, theAj);
      } // if

      else { // has squat-all
        return new Dirichlet<AssignType>(lhsPtr->dirVarId, lhsPtr->domain);
      } // else

    } // else does not have a Dirichlet

  } // if keep the Dirichlet

  // ++++++++++++++++++++++++++++++++
  // +++++ keep the Categorical +++++
  // ++++++++++++++++++++++++++++++++

  if ( sortedKeep[0] == lhsPtr->catVarId and sortedKeep.size() == 1 ) {

    // we have 4 cases to consider depending on whether or not dirPtr
    // and/or catPtr is instantiated

    if (lhsPtr->catPtr) {

      lhsPtr->catPtr->inplaceNormalize(nullptr);

      // if no Dirichlet, simply return the Categorical
      if (!lhsPtr->dirPtr) {
        return lhsPtr->catPtr->copy();
      } // if

      // ok, so we have a Dirichlet and will have to update the alphas
      // pk'= \frac{pk*ak}{\sum_j pj*aj}

      // to do this effectively we will use the double ladder
      // trick. To make things fast we will save the intermediate
      // results in two corresponding maps
      std::map< std::vector<AssignType>, FProb> thePjAj; // the pj*aj terms
      double sumPjAj = 0.0;
      const std::map<AssignType, double>& mapL(lhsPtr->dirPtr->getSparseAlphas() );
      const std::map<std::vector<AssignType>, FProb>& mapR(lhsPtr->catPtr->getSparseProbs() );
      double defL = lhsPtr->dirPtr->getDefaultAlpha();
      double defR = lhsPtr->catPtr->getDefaultProb();

      auto iL = mapL.begin();
      auto iR = mapR.begin();
      while ( iL != mapL.end() && iR != mapR.end() ){

        if (iL->first < iR->first[0]) {      // mapL,defR
          double theirProd = iL->second*defR;
          thePjAj.insert(thePjAj.end(), std::make_pair(std::vector<AssignType>({iL->first}), FProb(theirProd) ) ); // hinted insert for speed
          sumPjAj += theirProd;
          iL++;
        } // if

        else if (iL->first > iR->first[0]) { // defL,mapR
          double theirProd = defL*iR->second.prob;
          thePjAj.insert(thePjAj.end(), std::make_pair(iR->first, FProb(theirProd) ) );
          sumPjAj += theirProd;
          iR++;
        } // else if

        else {                            // mapL,mapR
          double theirProd = iL->second*iR->second.prob;
          thePjAj.insert(thePjAj.end(), std::make_pair(std::vector<AssignType>({iL->first}), FProb(theirProd) ) );
          sumPjAj += theirProd;
          iL++; iR++;
        } // else

      } // while

      while ( iL != mapL.end() ){         // some extra mapL,defR
        double theirProd = iL->second*defR;
        thePjAj.insert(thePjAj.end(), std::make_pair(std::vector<AssignType>({iL->first}), FProb(theirProd) ) );
        sumPjAj += theirProd;
        iL++;
      } // while

      while ( iR != mapR.end() ){         // some extra defL,mapR
        double theirProd = defL*iR->second.prob;
        thePjAj.insert(thePjAj.end(), std::make_pair(iR->first, FProb(theirProd) ) );
        sumPjAj += theirProd;
        iR++;
      } // while

      // thePjAj.size() is the combined updated size. Subtracting that from
      // noOfCategories gives the remaining number of unseen defaults.
      double defProd = defL*defR;
      sumPjAj += defProd * ( lhsPtr->dirPtr->noOfCategories()-thePjAj.size() );
      if (sumPjAj > 0.0) {
        defProd /= sumPjAj;
      } // if
      else {
        defProd = 0.0;
      } // else

      auto itr2 = thePjAj.begin();
      while ( itr2 != thePjAj.end() ){
        if (sumPjAj > 0.0) {
          itr2->second.prob /= sumPjAj;
        } // if
        else {
          itr2->second.prob = 0.0;
        } // else
        itr2++;
      } //while
      DiscreteTable<AssignType> * retPtr =
        new DiscreteTable<AssignType>({lhsPtr->catVarId}, {lhsPtr->domain}, defProd, thePjAj, 0.0, 0.0, true);
      retPtr->copyExtraSettingsFrom(*lhsPtr->catPtr);
      return retPtr;

    } // if has a Categorical

    else { // no Categorical

      if (lhsPtr->dirPtr) { // has only a Dirichlet

        const std::map<AssignType, double>& mapL(lhsPtr->dirPtr->getSparseAlphas() );
        std::map< std::vector<AssignType>, FProb> theProbs; // the aj terms
        double sumOfAlphas = lhsPtr->dirPtr->sumAlphas();

        if (sumOfAlphas > 0.0) {
          for (auto iL = mapL.begin(); iL != mapL.end(); iL++) {
            theProbs.insert(theProbs.end(), std::make_pair( std::vector<AssignType>({iL->first}), FProb(iL->second/sumOfAlphas) ) );
          } // for

          return new DiscreteTable<AssignType>(
            {lhsPtr->catVarId}, {lhsPtr->domain}, lhsPtr->dirPtr->getDefaultAlpha()/sumOfAlphas, theProbs, 0.0, 0.0, true);
        } // if
        else {
           return new DiscreteTable<AssignType>(
            {lhsPtr->catVarId}, {lhsPtr->domain}, 0.0, theProbs, 0.0, 0.0, true);
        } // else
      } // if

      else { // has squat-all
        DiscreteTable<AssignType> * retPtr = new DiscreteTable<AssignType>({lhsPtr->catVarId}, {lhsPtr->domain});
        retPtr->inplaceNormalize(nullptr);
        return retPtr;
      } // else


    } // else no Categorical

  } // if keep the Categorical

  // ++++++++++++++++++++++++++++++++++++++++
  // +++++ WTF! This should not happen! +++++
  // ++++++++++++++++++++++++++++++++++++++++

  // Once we reached here we are in trouble. Apparently there is nothing to keep!
  PRLITE_ASSERT(false,
  "Can not marginalise " << lhsPtr->dirVarId
  << " or " << lhsPtr->catVarId
  << " to retain " << sortedKeep);

  return nullptr;
} // process

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Polya_ObserveAndReduce<AssignType>::isA() const {
  static const std::string CLASSNAME("Polya_ObserveAndReduce<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
 * PURPOSE: In principle we can observe the dirichlet variable
 * resulting in a plain DiscreteTable where the various assignments
 * will have probs as dictated by the sampled thetas. I guess that if
 * we also have an existing categorical these probs will simply
 * multiply. However, all of this seems to be somewhat
 * counter-intuitive to what this class is actually intended for. We
 * will therefore only support observing the categorical variable.
 *
 * METHOD: Rather simple. If you observe a DiscreteTable X = k, then the
 * k'th alpha count in Dirichlet should be incremented by one.
 */
template <typename AssignType>
Factor*
Polya_ObserveAndReduce<AssignType>::process(
  const Polya<AssignType>* lhsPtr,
  const emdw::RVIds& theVars,
  const emdw::RVVals& anyVals,
  bool presorted) {

  Dirichlet<AssignType> * retPtr;

  const Polya<AssignType>& lhs(*lhsPtr);
  PRLITE_ASSERT(theVars.size() == 1 and theVars[0] == lhs.catVarId, "Can only observe the categorical variable in a Polya");

  if (lhsPtr->dirPtr) {retPtr = lhsPtr->dirPtr->copy();} // if
  else {
    retPtr = new Dirichlet<AssignType> (lhsPtr->dirVarId, lhsPtr->domain );
  } // else
  retPtr->addToCounts(AssignType(anyVals[0]), 1.0);

  return retPtr;


  //**********************************************************************************


} // process


/*
  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedVarsPtr = &theVars;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( theVars, std::less<unsigned>() );
    tmpVars = extract<unsigned>(theVars,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  else{
    //!@! sorted must be set up here. also remember that the
    //!observations can also include variables that do not form part
    //!of this class. We will probably have to do some intersection
    //!thing here.
  } // else
  const emdw::RVIds& sortedVars(*sortedVarsPtr);
 */

/*
template <typename AssignType>
Factor* Polya_Absorb<AssignType>::process(const Polya<AssignType>* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  Polya<AssignType>* fPtr = new Polya<AssignType>(*lhsPtr);

  // inplace normalize it
  Polya_InplaceAbsorb<AssignType> ipAbsorb;
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
  //std::cout << __FILE__ << __LINE__ << ": " << *fPtr << std::endl;
  return fPtr;

} // process

template <typename AssignType>
Factor* Polya_Cancel<AssignType>::process(const Polya<AssignType>* lhsPtr,
                                      const Factor* rhsFPtr) {

  // make a copy of the factor
  Polya<AssignType>* fPtr = new Polya<AssignType>(*lhsPtr);

  // inplace normalize it
  Polya_InplaceCancel<AssignType> ipCancel;
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

 */
