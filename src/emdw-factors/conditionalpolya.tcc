/*******************************************************************************

          AUTHORS: JA du Preez

          DATE: October 2018

          PURPOSE: ConditionalPolya PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


//emdw headers
#include "emdw.hpp"
#include "sortindices.hpp"
#include "vecset.hpp"

// standard headers
#include <new>  // placement new
#include <math.h>
#include <random>

/**************************The ConditionalPolya class members***********************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper
// function. Boilerplate.
template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new ConditionalPolya_Marginalize<CondType, DirType>() );
  } // if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new ConditionalPolya_InplaceNormalize<CondType, DirType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new ConditionalPolya_Normalize<CondType, DirType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new ConditionalPolya_InplaceAbsorb<CondType, DirType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new ConditionalPolya_Absorb<CondType, DirType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new ConditionalPolya_InplaceCancel<CondType, DirType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new ConditionalPolya_Cancel<CondType, DirType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new ConditionalPolya_ObserveAndReduce<CondType, DirType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper = uniqptr<InplaceDamper>( new ConditionalPolya_InplaceDamping<CondType, DirType>() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new ConditionalPolya_Sampler<CondType, DirType>() );
  } // if

} // supplementWithDefaultOperators

// class specific ctor
template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>::ConditionalPolya(
  emdw::RVIdType theDVar,
  emdw::RVIdType theZVar,
  emdw::RVIdType theXVar,
  rcptr< std::vector<CondType> > theZDomain,
  rcptr< std::vector<DirType> > theXDomain,
  ValidValue<DirType> theObservedX,
  rcptr< DirichletSet<CondType, DirType> > theCDirPtr,
  rcptr< DiscreteTable<CondType> > theZCatPtr,
  rcptr< DiscreteTable<DirType> > theXCatPtr,
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
: dVar(theDVar),
  zVar(theZVar),
  xVar(theXVar),
  zDomain(theZDomain),
  xDomain(theXDomain),
  observedX(theObservedX),
  cDirPtr(0),
  zCatPtr(0),
  xCatPtr(0),
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

  PRLITE_ASSERT(dVar < zVar and zVar < xVar, "dVar=" << dVar << "; zVar=" << zVar << "; xVar=" << xVar);

  if (theCDirPtr) {
    PRLITE_ASSERT(
      theCDirPtr->getVar(0) == dVar,
      "DirichletSet variable should be " << dVar
      << ", not " << theCDirPtr->getVar(0) );
    cDirPtr = uniqptr< DirichletSet<CondType, DirType> >( theCDirPtr->copy() );
  } // if

  if (theZCatPtr) {
    PRLITE_ASSERT(theZCatPtr->noOfVars() == 1, "Only one RV per (conditioned on) z DiscreteTable");
    PRLITE_ASSERT(
      theZCatPtr->getVar(0) == zVar,
      "Conditioned on (z) variable should be " << zVar
      << ", not " << theZCatPtr->getVar(0) );
    zCatPtr = uniqptr< DiscreteTable<CondType> >( theZCatPtr->copy() );
  } // if

  if (theXCatPtr and !observedX.valid) {
    PRLITE_ASSERT(theXCatPtr->noOfVars() == 1, "Only one RV per (target) y DiscreteTable");
    PRLITE_ASSERT(theXCatPtr->getVar(0) == xVar,
           "Target  (y) variable should be " << xVar
           << ", not " << theXCatPtr->getVar(0) );
    xCatPtr = uniqptr< DiscreteTable<DirType> >( theXCatPtr->copy() );
  } // if

  updateVols();
} // class specific constructor

// default constructor
template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>::ConditionalPolya(
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
      sampler(theSampler) {

  supplementWithDefaultOperators();

} // default constructor

// copy constructor.
template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>::ConditionalPolya(const ConditionalPolya<CondType, DirType>& p)
    : dVar(p.dVar),
      zVar(p.zVar),
      xVar(p.xVar),
      zDomain(p.zDomain),
      xDomain(p.xDomain),
      observedX(p.observedX),
      cDirPtr(0),
      zCatPtr(0),
      xCatPtr(0),
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

  if (p.cDirPtr) {
    cDirPtr = uniqptr< DirichletSet<CondType, DirType> >( p.cDirPtr->copy() );
  } // if

  if (p.zCatPtr) {
    zCatPtr = uniqptr< DiscreteTable<CondType> >( p.zCatPtr->copy() );
  } // if

  if (p.xCatPtr) {
    xCatPtr = uniqptr< DiscreteTable<DirType> >( p.xCatPtr->copy() );
  } // if

  updateVols();
} // cpy ctor

// move constructor.
template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>::ConditionalPolya(ConditionalPolya<CondType, DirType>&& p)
    : dVar(p.dVar),
      zVar(p.zVar),
      xVar(p.xVar),
      zDomain(p.zDomain),
      xDomain(p.xDomain),
      observedX(p.observedX),
      cDirPtr(p.cDirPtr),
      zCatPtr(p.zCatPtr),
      xCatPtr(p.xCatPtr),
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
  updateVols();
} // move ctor

template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>::~ConditionalPolya() {
  // TIP: Typically empty, except if your class contains raw pointers
} // destructor

//========================= Operators and Conversions =========================

template <typename CondType, typename DirType>
bool ConditionalPolya<CondType, DirType>::operator==(const ConditionalPolya<CondType, DirType>& p) const {

  // variables
  if (dVar != p.dVar) {return false;} // if
  if (zVar != p.zVar) {return false;} // if
  if ( xVar != p.xVar) {return false;} // if

  // domains
  if (*zDomain != *p.zDomain) {return false;} // if
  if (*xDomain != *p.xDomain) {return false;} // if

  if (observedX != p.observedX) {return false;} // if

  // check component distributions
  // DirichletSet
  if (cDirPtr and p.cDirPtr and *cDirPtr != *p.cDirPtr) {return false;} // if
  if (cDirPtr and !p.cDirPtr and cDirPtr->distanceFromVacuous() > 1E-10) {return false;} // if
  if (!cDirPtr and p.cDirPtr and p.cDirPtr->distanceFromVacuous() > 1E-10) {return false;} // if

  // z DiscreteTable
  if (zCatPtr and p.zCatPtr and *zCatPtr != *p.zCatPtr) {return false;} // if
  if (zCatPtr and !p.zCatPtr and zCatPtr->distanceFromVacuous() > 1E-10) {return false;} // if
  if (!zCatPtr and p.zCatPtr and p.zCatPtr->distanceFromVacuous() > 1E-10) {return false;} // if

  // y DiscreteTable
  if (xCatPtr and p.xCatPtr and *xCatPtr != *p.xCatPtr) {return false;} // if
  if (xCatPtr and !p.xCatPtr and xCatPtr->distanceFromVacuous() > 1E-10) {return false;} // if
  if (!xCatPtr and p.xCatPtr and p.xCatPtr->distanceFromVacuous() > 1E-10) {return false;} // if

  return true;
} // operator==

// Boilerplate
template <typename CondType, typename DirType>
bool ConditionalPolya<CondType, DirType>::operator!=(const ConditionalPolya<CondType, DirType>& p) const {
  return !operator==(p);
} // operator !=

// Boilerplate code
template <typename CondType, typename DirType>
bool ConditionalPolya<CondType, DirType>::isEqual(const Factor* rhsPtr) const {
  const ConditionalPolya<CondType, DirType>* dwnPtr = dynamic_cast<const ConditionalPolya<CondType, DirType>*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// Boilerplate
template <typename CondType, typename DirType>
unsigned ConditionalPolya<CondType, DirType>::classSpecificConfigure(
  emdw::RVIdType theDVar,
  emdw::RVIdType theZVar,
  emdw::RVIdType theXVar,
  rcptr< std::vector<CondType> > theZDomain,
  rcptr< std::vector<DirType> > theXDomain,
  ValidValue<DirType> theObservedX,
  rcptr< DirichletSet<CondType, DirType> > theCDirPtr,
  rcptr< DiscreteTable<CondType> > theZCatPtr,
  rcptr< DiscreteTable<DirType> > theXCatPtr,
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

  this->~ConditionalPolya<CondType, DirType>(); // Destroy existing

  new(this) ConditionalPolya<CondType, DirType>(
    theDVar,
    theZVar,
    theXVar,
    theZDomain,
    theXDomain,
    theObservedX,
    theCDirPtr,
    theZCatPtr,
    theXCatPtr,
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

template <typename CondType, typename DirType>
std::ostream& ConditionalPolya<CondType, DirType>::txtWrite(std::ostream& file) const {
  file << "ConditionalPolya_V0\n";
  file << "dVarId: " << dVar << " zVarId: " << zVar << " xVarId: " << xVar << std::endl;
  file << "\nzDomain: " << *zDomain << std::endl;
  file << "\nxDomain: " << *xDomain << std::endl;

  file << "yObserved:\n";
  if (observedX.valid ) { file << "1 " << observedX.value << std::endl;}
  else {file << "0\n";}

  file << "DirichletSet:\n";
  if (cDirPtr.get() ) { file << "1 " << *cDirPtr << std::endl;}
  else {file << "0\n";}

  file << "zDiscreteTable:\n";
  if (zCatPtr.get() ) { file << "1 " << *zCatPtr << std::endl;}
  else {file << "0\n";}

  file << "yDiscreteTable:\n";
  if (xCatPtr.get() ) { file << "1 " << *xCatPtr << std::endl;}
  else {file << "0\n";}

  // TODO Remember the factor operators

  return file;
} // txtWrite

template <typename CondType, typename DirType>
std::istream& ConditionalPolya<CondType, DirType>::txtRead(std::istream& file) {
  std::string dummy;

  file >> dummy; // strip the version tag

  // TODO Remember the factor operators
  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
template <typename CondType, typename DirType>
emdw::RVVals ConditionalPolya<CondType, DirType>::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
template <typename CondType, typename DirType>
uniqptr<Factor> ConditionalPolya<CondType, DirType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> ConditionalPolya<CondType, DirType>::absorb(const Factor* rhsPtr,
                                                            FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::inplaceAbsorb(const Factor* rhsPtr,
                                                        FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> ConditionalPolya<CondType, DirType>::cancel(const Factor* rhsPtr,
                                                            FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::inplaceCancel(const Factor* rhsPtr,
                                                        FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
template <typename CondType, typename DirType>
uniqptr<Factor> ConditionalPolya<CondType, DirType>::marginalize(const emdw::RVIds& variablesToKeep,
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
uniqptr<Factor> ConditionalPolya<CondType, DirType>::observeAndReduce(const emdw::RVIds& variables,
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
double ConditionalPolya<CondType, DirType>::inplaceDampen(const Factor* rhsPtr,
                                                          double df,
                                                          FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>* ConditionalPolya<CondType, DirType>::copy(
  const emdw::RVIds& newVars,
  bool presorted) const {

  if ( newVars.size() ) { // replace the variables
    return new ConditionalPolya<CondType, DirType>(
      newVars[0],
      newVars[1],
      newVars[2],
      zDomain,
      xDomain,
      observedX,
      cDirPtr,
      zCatPtr,
      xCatPtr,
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
  } // if newVars

  // identical copy
  return new ConditionalPolya<CondType, DirType>(*this);

} // copy

// Currently we retain the status of observedX. Not sure if this is
// right.
template <typename CondType, typename DirType>
ConditionalPolya<CondType, DirType>*
ConditionalPolya<CondType, DirType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  PRLITE_ASSERT(!selVars.size(), "ConditionalPolya::vacuousCopy : variable selection not supported");
  ConditionalPolya<CondType, DirType>* fPtr = copy();

  if (fPtr->cDirPtr) {
    fPtr->cDirPtr = uniqptr< DirichletSet<CondType,DirType> >( fPtr->cDirPtr->vacuousCopy() );
  } // if
  if (fPtr->zCatPtr) {
    fPtr->zCatPtr = uniqptr< DiscreteTable<CondType> >( fPtr->zCatPtr->vacuousCopy() );
  } // if
  if (fPtr->xCatPtr) {
    fPtr->xCatPtr = uniqptr< DiscreteTable<DirType> >( fPtr->xCatPtr->vacuousCopy() );
  } // if

  return fPtr;
} // vacuousCopy

template <typename CondType, typename DirType>
unsigned ConditionalPolya<CondType, DirType>::noOfVars() const {
  return observedX.valid ? 2 : 3;
} // noOfVars

template <typename CondType, typename DirType>
emdw::RVIds ConditionalPolya<CondType, DirType>::getVars() const {
  return observedX.valid ? emdw::RVIds{dVar,zVar} : emdw::RVIds{dVar,zVar, xVar};
} // getVars

template <typename CondType, typename DirType>
emdw::RVIdType ConditionalPolya<CondType, DirType>::getVar(unsigned varNo) const {

  switch (varNo) {
      case  0: {return dVar; break;} // case 0:
      case  1: {return zVar; break;} // case 1:
      case  2: {return xVar; break;} // case 2:
      default: {PRLITE_ASSERT(false, "requested var no > 2"); break; } //  default:
  } // switch

} // getVar

//====================== Other Virtual Member Functions =======================

template <typename CondType, typename DirType>
double ConditionalPolya<CondType, DirType>::distance(const Factor* rhsPtr) const {

  const ConditionalPolya<CondType, DirType>* dwnPtr = dynamic_cast<const ConditionalPolya<CondType, DirType>*>(rhsPtr);
  PRLITE_ASSERT(dwnPtr, "ConditionalPolya<CondType, DirType>::distance : rhs is not a ConditionalPolya");
  PRLITE_ASSERT(dVar == dwnPtr->dVar, "The DirichletSet variables do not match");
  PRLITE_ASSERT(zVar == dwnPtr->zVar, "The conditioned on variables do not match");
  PRLITE_ASSERT(xVar == dwnPtr->xVar, "The target variables do not match");

  // will only do this from one side. could possibly miss some
  // non-zero rhs probs.

  double dist = 0.0;

  if (cDirPtr) {
    if (dwnPtr->cDirPtr) {dist += cDirPtr->distance( dwnPtr->cDirPtr.get() );}
    else {return std::numeric_limits<double>::infinity();}
  } // if

  if (zCatPtr) {
    if (dwnPtr->zCatPtr) {dist += zCatPtr->distance( dwnPtr->zCatPtr.get() );}
    else {return std::numeric_limits<double>::infinity();}
  } // if

  if (observedX.valid) {
    if (!dwnPtr->observedX.valid or observedX.value != dwnPtr->observedX.value) {
      return std::numeric_limits<double>::infinity();
    } // if
  } // if

  else if (xCatPtr) {
    if (!dwnPtr->observedX.valid and dwnPtr->xCatPtr) {dist += xCatPtr->distance( dwnPtr->xCatPtr.get() );}
    else {return std::numeric_limits<double>::infinity();}
  } // if

  return dist;
} // distance

//======================= Parameterizations (non-virtual) =====================


//======================= Non-virtual Member Functions ========================

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::absorbCDir(const DirichletSet<CondType, DirType>* theCDirPtr){

  PRLITE_ASSERT(dVar == theCDirPtr->getVar(0),
         "ConditionalPolya_InplaceAbsorb : Both DirichletSets must be defined over the same random variable ids");

  // QUESTION???: Is this right at all? Shouldn't we just re-initialise the DirichletSet?
  // QUESTION???: Should we be paying attention to whether y is observed or not?
  if (cDirPtr) {
    cDirPtr->inplaceAbsorb(theCDirPtr, nullptr);
  } // if
  else {
    cDirPtr = uniqptr< DirichletSet<CondType,DirType> >( theCDirPtr->copy() );
  } // else

  updateVols();
} // absorbCDir

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::cancelCDir(const DirichletSet<CondType, DirType>* theCDirPtr){
    PRLITE_ASSERT(dVar == theCDirPtr->getVar(0),
           "ConditionalPolya_InplaceCancel : Both DirichletSets must be defined over the same random variable ids");
    PRLITE_ASSERT(cDirPtr, "ConditionalPolya_InplaceCancel : Trying to cancel a non-existent DirichletSet");
    cDirPtr->inplaceCancel(theCDirPtr,nullptr);
  updateVols();
} // cancelCDir

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::absorbZCat(const DiscreteTable<CondType>* theZCatPtr){

  if (zCatPtr) {
    zCatPtr->inplaceAbsorb(theZCatPtr,nullptr);
  } // if
  else {
    zCatPtr = uniqptr< DiscreteTable<CondType> >(
      static_cast<DiscreteTable<CondType>*>( theZCatPtr->copy() ) );
  } // else

  updateVols();
} // absorbZCat

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::cancelZCat(const DiscreteTable<CondType>* theZCatPtr){
  PRLITE_ASSERT(zCatPtr, "ConditionalPolya_InplaceCancel : Trying to cancel a non-existent z DiscreteTable");
  zCatPtr->inplaceCancel(theZCatPtr,nullptr);
  updateVols();
} // cancelZCat

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::absorbXCat(const DiscreteTable<DirType>* theXCatPtr){
  if (observedX.valid) {
    xCatPtr = nullptr;
  } // if
  else if (xCatPtr) {
    xCatPtr->inplaceAbsorb(theXCatPtr,nullptr);
  } // if
  else {
    xCatPtr = uniqptr< DiscreteTable<DirType> >(
      static_cast<DiscreteTable<DirType>*>( theXCatPtr->copy() ) );
  } // else
  updateVols();
} // absorbXCat

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::cancelXCat(const DiscreteTable<DirType>* theXCatPtr){
  PRLITE_ASSERT(xCatPtr, "ConditionalPolya_InplaceCancel : Trying to cancel a non-existent y DiscreteTable");
  xCatPtr->inplaceCancel(theXCatPtr,nullptr);
  updateVols();
} // cancelXCat

template <typename CondType, typename DirType>
ValidValue<DirType> ConditionalPolya<CondType, DirType>::getObservedX() const {
  return observedX;
} // getObservedX

/// nRows, nCols
template <typename CondType, typename DirType>
std::pair<size_t,size_t> ConditionalPolya<CondType, DirType>::volDims() const {
  return std::pair<size_t,size_t>( vols.rows(), vols.cols() );
} // volDims

/// volume for z-row, y-col
template <typename CondType, typename DirType>
double ConditionalPolya<CondType, DirType>::vol(size_t rZ, size_t cY) const {
  return vols(rZ,cY);
} // vol

template <typename CondType, typename DirType>
void ConditionalPolya<CondType, DirType>::updateVols() {


  // we only bother with this if both Phi and z are multiplied in
  if (!cDirPtr or !zCatPtr) {return;} // if

  inplaceNormalize();

  double tVol = 0.0;
  if (observedX.valid) {
    vols.resize(zDomain->size(), 1);
    size_t zIdx = 0;
    for (auto zVal : *zDomain) {
      double& cell( vols(zIdx,0) );

      rcptr< Dirichlet<DirType> > dirPtr(cDirPtr ? cDirPtr->aliasDirichlet(zVal) : nullptr);
      if (dirPtr) { // has a Dirichlet
        cell = dirPtr->getAlpha(observedX.value);
        cell /= dirPtr->sumAlphas();
      } // if
      else {cell = 1.0;}

      if (zCatPtr) { // has a distribution over the conditioned-on (z) variable
        cell *= zCatPtr->getProb({zVal});
      } // if

      tVol += cell;
      zIdx++;
    } // for
  } // if

  else {
    vols.resize( zDomain->size(), xDomain->size() );
    size_t zIdx = 0;

    for (auto zVal : *zDomain) {
      rcptr< Dirichlet<DirType> > dirPtr(cDirPtr ? cDirPtr->aliasDirichlet(zVal) : nullptr);
      if (dirPtr) { } // if

      size_t yIdx = 0;
      for (auto yVal : *xDomain) {
        double& cell( vols(zIdx,yIdx) );

        if (dirPtr) { // has a Dirichlet
          cell = dirPtr->getAlpha(yVal);
          cell /= dirPtr->sumAlphas();
        } // if
        else {cell = 1.0;} // else

        if (zCatPtr) { // has a distribution over the conditioned-on (z) variable
          cell *= zCatPtr->getProb({zVal});
        } // if

        if (xCatPtr) { // has a distribution over the conditioned-on (z) variable
          cell *= xCatPtr->getProb({yVal});
        } // if

        tVol += cell;
        yIdx++;
      } // for

      zIdx++;
    } // for

  } // else

  for (int zIdx = 0; zIdx < vols.rows(); zIdx++) {
    for (int yIdx = 0; yIdx < vols.cols(); yIdx++) {
      vols(zIdx,yIdx) /= tVol;
    } // for
  } // for

} // updateVols

//================================== Friends ==================================

/// output. Boilerplate
template<typename T1, typename T2>
std::ostream& operator<<( std::ostream& file,
                          const ConditionalPolya<T1,T2>& n ){
  return n.txtWrite(file);
} // operator<<

/// input. Boilerplate
template<typename T1, typename T2>
std::istream& operator>>( std::istream& file,
                          ConditionalPolya<T1,T2>& n ){
  return n.txtRead(file);
} // operator>>


//======================= Factor operators ====================================


// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample ****************************************
// ***************************************************************************

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_Sampler<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_Sampler<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename CondType, typename DirType>
emdw::RVVals ConditionalPolya_Sampler<CondType, DirType>::sample(
  const ConditionalPolya<CondType, DirType>* lhsPtr) {

  emdw::RVVals samples;
  PRLITE_ASSERT(lhsPtr->zCatPtr, "ConditionalPolya_Sampler:: The to-condition-on (z) distribution is uninstantiated");
  PRLITE_ASSERT(lhsPtr->cDirPtr, "ConditionalPolya_Sampler:: The DirichletSet distribution is uninstantiated");
  samples.push_back( lhsPtr->zCatPtr->sample(nullptr)[0] );
  CondType zVal = CondType(samples[0]);
  emdw::RVVals anyThetas = lhsPtr->cDirPtr->sampleFrom(zVal);
  std::vector<double> thetas;
  for (auto el : anyThetas) {thetas.push_back(double(el) );} // for
  // ironically we can not push AnyType values in the array
  samples.insert( samples.end(), thetas.begin(), thetas.end() );
  if (lhsPtr->observedX.valid) {samples.push_back(lhsPtr->observedX.value);} // if
  else {
    // build the y Categorical and sample from it.
    std::map<std::vector<DirType>, FProb> sparseProbs;
    for (unsigned idx =0; idx < thetas.size(); idx++) {
      sparseProbs[{(*lhsPtr->xDomain)[idx]}] = thetas[idx];
    } // for
    rcptr<Factor> polyaCatPtr = uniqptr< DiscreteTable<DirType> > (
      new DiscreteTable<DirType>({lhsPtr->xVar}, {lhsPtr->xDomain}, 0.0, sparseProbs) );

    // guess we still need to take a pre-existing xCatPtr into
    // account? This feels just plain wrong
     // if (lhsPtr->xCatPtr) {
     //   polyaCatPtr->inplaceAbsorb(lhsPtr->xCatPtr);
     //   polyaCatPtr->inplaceNormalize();
     // } // if

    // ironically we can not push AnyType values in the array
    samples.push_back( DirType(polyaCatPtr->sample()[0]) );
  } // else

  return samples;
} // sample


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_InplaceNormalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_InplaceNormalize<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

template <typename CondType, typename DirType>
void ConditionalPolya_InplaceNormalize<CondType, DirType>::inplaceProcess(
  ConditionalPolya<CondType, DirType>* lhsPtr) {

  ConditionalPolya<CondType, DirType>& lhs(*lhsPtr);

  if (lhs.cDirPtr) {lhs.cDirPtr->inplaceNormalize(nullptr);} // if
  if (lhs.zCatPtr) {lhs.zCatPtr->inplaceNormalize(nullptr);} // if
  if (lhs.xCatPtr) {lhs.xCatPtr->inplaceNormalize(nullptr);} // if

} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_Normalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_Normalize<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Boilerplate
template <typename CondType, typename DirType>
Factor* ConditionalPolya_Normalize<CondType, DirType>::process(
  const ConditionalPolya<CondType, DirType>* lhsPtr) {

  // make a copy of the factor
  ConditionalPolya<CondType, DirType>* fPtr = new ConditionalPolya<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  ConditionalPolya_InplaceNormalize<CondType, DirType> ipNorm;
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
/********************************** Absorb ***********************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_InplaceAbsorb<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_InplaceAbsorb<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * The ConditionalPolya can get multiplied with 1) a categorical
 * distribution p(z), where z is the conditioned-on variable, 2) a
 * DirichletSet distribution describing the topic
 * distributions, or 3) a Categorical distribution p(y) where y is the
 * items/words that makes up a topic.
 */
template <typename CondType, typename DirType>
void ConditionalPolya_InplaceAbsorb<CondType, DirType>::inplaceProcess(
  ConditionalPolya<CondType, DirType>* lhsPtr,
  const Factor* rhsFPtr){

  // a convenient reference for the lhs factor
  ConditionalPolya<CondType, DirType>& lhs(*lhsPtr);

  // check out what type of rhs Factor we've got and then handle it as needed
  PRLITE_ASSERT( rhsFPtr->noOfVars() == 1, "ConditionalPolya InplaceAbsorb : Can only multiply with singel RV factors - in this case we have " << rhsFPtr->noOfVars() );

  // Multiplying with a DirichletSet
  const  DirichletSet<CondType, DirType> * rhsPtr2 =
    dynamic_cast<const  DirichletSet<CondType, DirType> *>(rhsFPtr);
  if (rhsPtr2) { // ok, we've got a DirichletSet
    lhs.absorbCDir(rhsPtr2);
    return;
  } // if

  // Multiplying with a z DiscreteTable
  const DiscreteTable<CondType>* rhsPtr3 = dynamic_cast<const DiscreteTable<CondType>*>(rhsFPtr);
  if (rhsPtr3 and lhs.zVar == rhsPtr3->getVar(0) ) { // ok, we've got a z DiscreteTable
    lhs.absorbZCat(rhsPtr3);
    return;
  } // if

  // Multiplying with a y DiscreteTable
  const DiscreteTable<DirType>* rhsPtr4 = dynamic_cast<const DiscreteTable<DirType>*>(rhsFPtr);
  if (rhsPtr4 and lhs.xVar == rhsPtr4->getVar(0) ) { // ok, we've got a y DiscreteTable
    lhs.absorbXCat(rhsPtr4);
    return;
  } // if

  PRLITE_ASSERT(
    false,
    "ConditionalPolya_InplaceAbsorb : rhs is either an invalid type: " << typeid(*rhsFPtr).name()
    << " or its id " << rhsFPtr->getVar(0) << "does not match "
    << lhs.zVar << " or " << lhs.xVar);

} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_Absorb<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_Absorb<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_InplaceCancel<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_InplaceCancel<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/**
 * The ConditionalPolya can get divided by 1) a categorical
 * distribution p(z), where z is the conditioned-on variable, 2) a
 * DirichletSet distribution describing the topic
 * distributions, or 3) a Categorical distribution p(y) where y is the
 * items/words that makes up a topic.
 */
template <typename CondType, typename DirType>
void ConditionalPolya_InplaceCancel<CondType, DirType>::inplaceProcess(
  ConditionalPolya<CondType, DirType>* lhsPtr,
  const Factor* rhsFPtr){

  // a convenient reference for the lhs factor
  ConditionalPolya<CondType, DirType>& lhs(*lhsPtr);

  // check out what type of rhs Factor we've got and then handle it as needed

  // Dividing by a DirichletSet
  const  DirichletSet<CondType, DirType> * rhsPtr2 =
    dynamic_cast<const  DirichletSet<CondType, DirType> *>(rhsFPtr);
  if (rhsPtr2) { // ok, we've got a DirichletSet
    lhs.cancelCDir(rhsPtr2);
    return;
  } // if

  // Dividing by a z DiscreteTable
  const DiscreteTable<CondType>* rhsPtr3 = dynamic_cast<const DiscreteTable<CondType>*>(rhsFPtr);
  if (rhsPtr3 and lhs.zVar == rhsPtr3->getVar(0) ) { // ok, we've got a z DiscreteTable
    lhs.cancelZCat(rhsPtr3);
    return;
  } // if

  // Dividing by a y DiscreteTable
  const DiscreteTable<DirType>* rhsPtr4 = dynamic_cast<const DiscreteTable<DirType>*>(rhsFPtr);
  if (rhsPtr4 and lhs.xVar == rhsPtr4->getVar(0) ) { // ok, we've got a y DiscreteTable
    lhs.cancelXCat(rhsPtr4);
    return;
  } // if

  PRLITE_ASSERT(
    false,
    "ConditionalPolya_InplaceCancel : rhs is either an invalid type: " << typeid(*rhsFPtr).name()
    << " or its id " << rhsFPtr->getVar(0) << "does not match either of "
    << lhs.zVar << " or " << lhs.xVar);

} // inplaceProcess

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_Cancel<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_Cancel<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_Marginalize<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_Marginalize<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
 * Special considerations: 1) y can be observed or not and 2) The
 * external y distribution might or might not exist.
 *
 * We can keep 1) The DirichletSet or 2) The target (y)
 * distribution or 3) The conditioned-on (z) distribution.  A fourth
 * possibility is to keep both of y and z - for the moment we won't
 * support that.
 */
template <typename CondType, typename DirType>
Factor* ConditionalPolya_Marginalize<CondType, DirType>::process(
  const ConditionalPolya<CondType, DirType>* lhsPtr,
  const emdw::RVIds& toKeep,
  bool presorted,
  const Factor* peekAheadPtr) {

  PRLITE_ASSERT(toKeep.size() == 1,
    "ConditionalPolya_Marginalize : We only allow keeping one of the DirichletSet, the target (y) DiscreteTable or the conditioned-on (z) DiscreteTable.");

  // a convenient reference for the lhs factor
  const ConditionalPolya<CondType, DirType>& lhs(*lhsPtr);

  if (toKeep[0] == lhs.dVar) {

    // ###############################################
    // ### Keep Phi (the DirichletSet RV)  ###
    // ###############################################

    std::map< CondType, rcptr< Dirichlet<DirType> > > dirsMap; // { {0u,refD0_DIR},{1u,refD1_DIR} };

    size_t r = 0;
    for (auto rIdx : *lhs.zDomain) {
      rcptr< Dirichlet<DirType> > theDirPtr;

      if (!lhs.cDirPtr) { // no Dirichlets: equivalent to all alphas = 1
        theDirPtr = uniqptr< Dirichlet<DirType> > ( new Dirichlet<DirType>(0u,lhs.xDomain, 1.0) );
      } // if

      else { // has got a Dirichlet and therefore also incoming message
        theDirPtr = uniqptr< Dirichlet<DirType> > ( lhs.cDirPtr->copyDirichlet(rIdx) );

        if (lhs.zCatPtr) { // properly initialised

          if (lhs.observedX.valid) {
            theDirPtr->addToCounts( lhs.observedX.value, lhs.vol(r,0) );
          } // if

          else {
            size_t c = 0;
            for (auto cIdx : *lhs.xDomain) {
              theDirPtr->addToCounts( cIdx, lhs.vol(r,c) );
              c++;
            } // for
          } // else

        } // if

      } // else

      dirsMap[rIdx] = theDirPtr;
      r++;
    } // for rIdx

    DirichletSet<CondType, DirType>* retPtr(
      new DirichletSet<CondType, DirType>(lhs.dVar,lhs.zDomain, dirsMap) );
    if (lhs.cDirPtr) {retPtr->copyExtraSettingsFrom(*lhs.cDirPtr);} // if
    return retPtr;
  } // if

  else if (toKeep[0] == lhs.zVar) {

    // ######################################
    // ### Keep z (the conditioned-on RV) ###
    // ######################################

    if (lhs.zCatPtr) { // has got a zCatPtr and therefore also incoming message

      if (lhs.cDirPtr) { // properly initialised

        std::map<std::vector<CondType>, FProb> sparseProbs;
        int r = 0;
        for (auto rIdx : *lhs.zDomain) {

          if (lhs.observedX.valid) {
            sparseProbs[{rIdx}] = lhs.vol(r,0);
          } // if

          else {
            double prob = 0;
            for (int c = 0; c < lhs.vols.cols(); c++) {
              prob += lhs.vol(r,c);
            } // for
            sparseProbs[{rIdx}] = prob;
          } // else

          r++;
        } // for r

        DiscreteTable<CondType>* retPtr = new DiscreteTable<CondType>(
          {lhs.zVar}, {lhs.zDomain}, 0.0, sparseProbs);
        retPtr->copyExtraSettingsFrom(*lhs.zCatPtr);

        return retPtr;
      } // if

      else { // no cDir -> return current zCatPtr to cancel incoming msg
        return lhs.zCatPtr->copy();
      } // else

    } // if

    else { // no zCat, therefore also no incoming message so far. Vacuous

      // !!!If we do have a DirichletSet, and the Y variable
      // is observed, we should be able to do something. However, my
      // measurements on the code below is a bit ambiguous about this.

      //*
      if (lhs.cDirPtr and lhs.observedX.valid) {

        std::map<std::vector<CondType>, FProb> sparseProbs;
        double tVol = 0.0;
        for (auto zVal : *lhs.zDomain) {
          rcptr< Dirichlet<DirType> > dirPtr( lhs.cDirPtr->aliasDirichlet(zVal) );
          double scaledProb = dirPtr->getAlpha(lhs.observedX.value)/dirPtr->sumAlphas();
          tVol += scaledProb;
          sparseProbs[{zVal}] = scaledProb;
        } // for
        for (auto& el : sparseProbs) {
          if (tVol > 0.0) {
          el.second.prob /= tVol;
          } // if
          else {
            el.second.prob = 0.0;
          } // else
        } // for

        return new DiscreteTable<CondType>({lhs.zVar},{lhs.zDomain},0.0, sparseProbs);
      } // if
      //*/

      return new DiscreteTable<CondType>({lhs.zVar},{lhs.zDomain});
    } // else

  } // else if

  else if (toKeep[0] == lhs.xVar) { //!!!

    // ##############################
    // ### Keep y (the target RV) ###
    // ##############################

    PRLITE_ASSERT(!lhs.observedX.valid, "ConditionalPolya_Marginalize : This is crazy, you can not marginalize to x if it is already observed!");

    if (lhs.cDirPtr and lhs.zCatPtr) {
      // if both z and phi are instantiated, return the marginalization from vols

      std::map<std::vector<DirType>, FProb> sparseProbs;
      int c = 0;
      for (auto cIdx : *lhs.xDomain) {

        double prob = 0;
        for (int r = 0; r < lhs.vols.rows(); r++) {
          prob += lhs.vol(r,c);
        } // for r
        sparseProbs[{cIdx}] = prob;

        c++;
      } // for c

      DiscreteTable<DirType>* retPtr = new DiscreteTable<DirType>(
        {lhs.xVar}, {lhs.xDomain}, 0.0, sparseProbs);
      retPtr->copyExtraSettingsFrom(*lhs.xCatPtr);

      return retPtr;
    } // if

    else if (lhs.xCatPtr) {
      // if we have an external xCatPtr, just return that to cancel the incoming message
      return lhs.xCatPtr->copy();
    } // else if

    else {
    // return a vacuous density
      return new DiscreteTable<DirType>({lhs.xVar},{lhs.xDomain});
    } // else

  } // else if

  PRLITE_ASSERT(
    false,
    "ConditionalPolya_Marginalize : We only allow keeping one of the DirichletSet, the target (y) DiscreteTable or the conditioned-on (z) DiscreteTable.");
  return nullptr;

} // process


/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

// Identify object
template <typename CondType, typename DirType>
const std::string& ConditionalPolya_ObserveAndReduce<CondType, DirType>::isA() const {
  static const std::string CLASSNAME("ConditionalPolya_ObserveAndReduce<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
  return CLASSNAME;
} // isA

/*
 * We will only cater for observing the value of y
 */
template <typename CondType, typename DirType>
Factor* ConditionalPolya_ObserveAndReduce<CondType, DirType>::process(const ConditionalPolya<CondType, DirType>* lhsPtr,
                                          const emdw::RVIds& theVars,
                                          const emdw::RVVals& anyVals,
                                          bool presorted) {

  // a convenient reference for the lhs factor
  const ConditionalPolya<CondType, DirType>& lhs(*lhsPtr);

  PRLITE_ASSERT(theVars.size() == 1 and theVars[0] == lhs.xVar,
         "ConditionalPolya_ObserveAndReduce : we only allow observing the target (y) value");


  ValidValue<DirType> yVal = DirType(anyVals[0]);

  return new ConditionalPolya<CondType, DirType> (
    lhs.dVar, lhs.zVar, lhs.xVar,
    lhs.zDomain, lhs.xDomain,
    yVal,
    lhs.cDirPtr, lhs.zCatPtr, nullptr,
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

/*
// Boilerplate
template <typename CondType, typename DirType>
Factor* ConditionalPolya_Absorb<CondType, DirType>::process(
  const ConditionalPolya<CondType, DirType>* lhsPtr,
  const Factor* rhsFPtr) {

  // make a copy of the factor
  ConditionalPolya<CondType, DirType>* fPtr = new ConditionalPolya<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  ConditionalPolya_InplaceAbsorb<CondType, DirType> ipAbsorb;
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


// Boilerplate
template <typename CondType, typename DirType>
Factor* ConditionalPolya_Cancel<CondType, DirType>::process(
  const ConditionalPolya<CondType, DirType>* lhsPtr,
  const Factor* rhsFPtr) {

  // make a copy of the factor
  ConditionalPolya<CondType, DirType>* fPtr = new ConditionalPolya<CondType, DirType>(*lhsPtr);

  // inplace normalize it
  ConditionalPolya_InplaceCancel<CondType, DirType> ipCancel;
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
