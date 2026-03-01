/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2018
          PURPOSE: Dirichlet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/
#include <numeric>
#include <random>

#include "mathfuncs.hpp"

/******************************************************************************/
/***************************** Dirichlet **************************************/
/******************************************************************************/


//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
template <typename AssignType>
void Dirichlet<AssignType>::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new Dirichlet_Marginalize<AssignType>() );
  } //if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new Dirichlet_InplaceNormalize<AssignType>() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new Dirichlet_Normalize<AssignType>() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new Dirichlet_InplaceAbsorb<AssignType>() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new Dirichlet_Absorb<AssignType>() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new Dirichlet_InplaceCancel<AssignType>() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new Dirichlet_Cancel<AssignType>() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new Dirichlet_ObserveAndReduce<AssignType>() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper =
      uniqptr<InplaceDamper>( new Dirichlet_InplaceWeakDamping<AssignType>() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new Dirichlet_Sampler<AssignType>() );
  } // if

} // supplementWithDefaultOperators

// default constructor
template <typename AssignType>
Dirichlet<AssignType>::Dirichlet(
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
      sampler(theSampler),
      sumOfAlphas(0.0){

  supplementWithDefaultOperators();
} // default constructor

// class specific ctor
template <typename AssignType>
Dirichlet<AssignType>::Dirichlet(
  emdw::RVIdType theVarId,
  rcptr< std::vector<AssignType> > theDomain,
  double theDefAlpha,
  const std::map<AssignType, double>& theSparseAlphas,
  double theRelMargin,
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
: varId(theVarId),
  domain(theDomain),
  defAlpha(theDefAlpha),
  sparseAlphas(theSparseAlphas),
  relMargin(theRelMargin),
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
  sumOfAlphas(0.0){

  supplementWithDefaultOperators();

  recalcSumOfAlphas(); // uses domain->size()
} // class specific ctor

// convenience ctor
template <typename AssignType>
Dirichlet<AssignType>::Dirichlet(
  rcptr< std::vector<AssignType> > theDomain,
  emdw::RVIdType theVarId,
  double nominalAlpha,
  double perturbAlpha,
  double theRelMargin,
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
: varId(theVarId),
  domain(theDomain),
  defAlpha(perturbAlpha > 0.0 ? 0.0 : nominalAlpha),
  sparseAlphas(),
  relMargin(theRelMargin),
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
  sumOfAlphas(0.0){

  supplementWithDefaultOperators();

  if (perturbAlpha > 0.0) {
    std::uniform_real_distribution<> perturb(0, perturbAlpha);
    for (auto el : *theDomain) {
      sparseAlphas[el] = nominalAlpha + perturb( emdw::randomEngine() );
    } // for
  } // if

  recalcSumOfAlphas(); // uses domain->size()
} // convenience ctor


// destructor
template <typename AssignType>
Dirichlet<AssignType>::~Dirichlet() {
} // destructor

//========================= Operators and Conversions =========================

// operator ==
template <typename AssignType>
bool Dirichlet<AssignType>::operator==(const Dirichlet<AssignType>& p) const {

  const double tol = std::numeric_limits<double>::epsilon()*10;

  // Test to see whether potentials are over the same set of variables
  if (varId != p.varId) {return false;} // if

  // check domain
  if (*domain != *p.domain) {return false;} // if

  // Check the sumOfAlphas
  if (std::abs(sumOfAlphas - p.sumOfAlphas) > tol) {return false;} // if

  // now check the corresponding values

  const std::map<AssignType, double>& mapL(sparseAlphas);         // convenient alias
  const double defL(defAlpha);                                 // convenient alias

  const std::map<AssignType, double>& mapR(p.sparseAlphas); // convenient alias
  const double defR(p.defAlpha);                         // convenient alias

  // to calc the number of defaults occurring in neither of the maps
  size_t nDefsL = domain->size()-mapL.size();
  size_t nDefsR = domain->size()-mapR.size();

  auto iL = mapL.begin();
  auto iR = mapR.begin();
  while ( iL != mapL.end() && iR != mapR.end() ){

    if (iL->first < iR->first) {      // only in mapL
      nDefsR--;
      if ( std::abs(iL->second - defR) > tol ) {return false;}
      iL++;
    } // if

    else if (iL->first > iR->first) { // only in mapR
      nDefsL--;
      if ( std::abs(iR->second - defL) > tol ) {return false;}
      iR++;
    } // else if

    else {                            //common to both maps
      if ( std::abs(iL->second - iR->second) > tol ) {return false;}
      iL++; iR++;
    } // else

  } // while

  while ( iL != mapL.end() ){         // some extra spec's in mapL
    nDefsR--;
    if ( std::abs(iL->second - defR) > tol ) {return false;}
    iL++;
  } // while

  while ( iR != mapR.end() ){         // some extra spec's in mapR
    nDefsL--;
    if ( std::abs(iR->second - defL) > tol ) {return false;}
    iR++;
  } // while

    // and now  all the cases simultaneously absent from both maps
  PRLITE_ASSERT(
    nDefsL == nDefsR,
    "The number of cases in neither map should match, no matter from which side you look at it.");

  if (nDefsL and std::abs(defL-defR) > tol){ return false;}

  return true;
} // operator ==

// operator !=
template <typename AssignType>
bool Dirichlet<AssignType>::operator!=(const Dirichlet<AssignType>& p) const {
  return !operator==(p);
} // operator !=

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

// classSpecificConfigure
template <typename AssignType>
unsigned Dirichlet<AssignType>::classSpecificConfigure(
  emdw::RVIdType theVarId,
  rcptr< std::vector<AssignType> > theDomain,
  double theDefAlpha,
  const std::map<AssignType, double>& theSparseAlphas,
  double theRelMargin,
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

  this->~Dirichlet<AssignType>(); // Destroy existing

  new(this) Dirichlet<AssignType>(
    theVarId,
    theDomain,
    theDefAlpha,
    theSparseAlphas,
    theRelMargin,
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
std::ostream& Dirichlet<AssignType>::txtWrite(std::ostream& file) const {
  file << "Dirichlet_V0\n";
  file << "VarId: " << varId << std::endl;

  file << "Domain: " << *domain << std::endl;

  file << "\nDefaultAlpha: " << defAlpha
       << " RelMarginAlpha: " << relMargin << std::endl;

  file << "\nNoOfSpecAlphas: " << sparseAlphas.size() << std::endl;
  for (auto spec : sparseAlphas) {
    file << spec.first << " " << spec.second << "   ";
  } // for
  file << std::endl;

  // we still miss the factor operators
  return file;
} // txtWrite

// txtRead
template <typename AssignType>
std::istream& Dirichlet<AssignType>::txtRead(std::istream& file) {

  std::string dummy;

  emdw::RVIdType theVarId;
  double theDefAlpha;
  size_t nSpec;
  double theRelMargin;

  file >> dummy; // version tag
  file >> dummy >> theVarId;

  rcptr< std::vector<AssignType> > theDomain;
  theDomain= uniqptr< std::vector<AssignType> >(new std::vector<AssignType>());
  file >> dummy >>  *theDomain;

  file >> dummy >> theDefAlpha;
  file >> dummy >> theRelMargin;
  file >> dummy >> nSpec;

  std::map<AssignType, double> theSparseAlphas;
  while (nSpec--) {
    AssignType theVal; double theCnt;
    PRLITE_CHECK(file >> theVal, IOError, "Error reading parameter value.");
    PRLITE_CHECK(file >> theCnt, IOError, "Error reading parameter alpha.");
    theSparseAlphas[theVal] = theCnt;
  } // while

  // we still miss the factor operators

  classSpecificConfigure(theVarId, theDomain, theDefAlpha, theSparseAlphas, theRelMargin);

  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
template <typename AssignType>
emdw::RVVals Dirichlet<AssignType>::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
template <typename AssignType>
uniqptr<Factor> Dirichlet<AssignType>::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
template <typename AssignType>
void Dirichlet<AssignType>::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Dirichlet<AssignType>::absorb(const Factor* rhsPtr,
                                              FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
template <typename AssignType>
void Dirichlet<AssignType>::inplaceAbsorb(const Factor* rhsPtr,
                                          FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Dirichlet<AssignType>::cancel(const Factor* rhsPtr,
                                              FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
template <typename AssignType>
void Dirichlet<AssignType>::inplaceCancel(const Factor* rhsPtr,
                                          FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
template <typename AssignType>
uniqptr<Factor> Dirichlet<AssignType>::marginalize(const emdw::RVIds& variablesToKeep,
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
uniqptr<Factor> Dirichlet<AssignType>::observeAndReduce(const emdw::RVIds& variables,
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
double Dirichlet<AssignType>::inplaceDampen(const Factor* rhsPtr,
                                            double df,
                                            FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

// copy
template <typename AssignType>
Dirichlet<AssignType>* Dirichlet<AssignType>::copy(
  const emdw::RVIds& newVars, bool presorted) const {
  if ( newVars.size() ) {
    return new Dirichlet<AssignType>(
               newVars[0],
               domain,
               defAlpha,
               sparseAlphas,
               relMargin,
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
  return new Dirichlet<AssignType>(*this);
} // copy

// vacuousCopy
template <typename AssignType>
Dirichlet<AssignType>*
Dirichlet<AssignType>::vacuousCopy(
  const emdw::RVIds& selVars,
  bool) const {

  PRLITE_ASSERT(selVars.size()<2, "Dirichlet::vacuousCopy : can not select more than one variable");
  if ( selVars.size() ) {
    PRLITE_ASSERT(selVars[0] == varId, "Dirichlet::vacuousCopy : Can only select variable "
           << varId);
  }// if

  Dirichlet<AssignType>* fPtr = copy();
  fPtr->defAlpha = 1.0;
  fPtr->sparseAlphas.clear();
  return fPtr;
} // vacuousCopy

// Boiler plate code
// isEqual
template <typename AssignType>
bool Dirichlet<AssignType>::isEqual(const Factor* rhsPtr) const {
  const Dirichlet<AssignType>* dwnPtr = dynamic_cast<const Dirichlet<AssignType>*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

// Symmetric KL distance from a uniform dirichlet. Several terms
// cancel due to the symmetricising.
// 2*symKL = \sum_k( (a_k-b_k)[dg(a_k)-dg(b_k)-dg(\sum_k a_k)+dg(\sum_k b_k)])

template <typename AssignType>
double Dirichlet<AssignType>::distanceFromVacuous() const {

  double dgSumAlphaP = digamma(sumAlphas());
  double dgSumAlphaQ = digamma(domain->size());
  double dg1 = digamma(1.0);

  // we first do the default alphas. remember all the betas==1
  double dist = 0.0;
  auto diffAlpha = defAlpha-1.0;
  auto nDefs = domain->size()-sparseAlphas.size();
  if (nDefs and diffAlpha != 0.0) {

    if (defAlpha > 0.0) {
      dist = nDefs*diffAlpha *
        (digamma(defAlpha)-dg1 -dgSumAlphaP + dgSumAlphaQ);
    } // if

    // if alpha_k == 0 and beta_k != 0 or vice versa, the distance is infinity
    else {return std::numeric_limits<double>::infinity();} // else
  } // if

  // and now we do the sparse alphas
  for (const auto& sCnt : sparseAlphas) {

    if (sCnt.second > 0.0){
      diffAlpha = sCnt.second-1.0;
      if (diffAlpha != 0.0) {
        dist +=  (diffAlpha)*(digamma(sCnt.second)- dg1 - dgSumAlphaP + dgSumAlphaQ);
      } // if
    } // if

    // if alpha_k == 0 and beta_k != 0 or vice versa, the distance is infinity
    else {return std::numeric_limits<double>::infinity();} // else
  } // for

  return dist/2;
} // distanceFromVacuous

// noOfVars
template <typename AssignType>
unsigned Dirichlet<AssignType>::noOfVars() const {
  return 1;
} // noOfVars

// getVars
template <typename AssignType>
emdw::RVIds Dirichlet<AssignType>::getVars() const {
  return emdw::RVIds{varId};
} // getVars

// getVar
template <typename AssignType>
emdw::RVIdType Dirichlet<AssignType>::getVar(unsigned varNo) const {
  PRLITE_ASSERT(varNo ==0, "No sense in asking for var no " << varNo << " there is only one");
  return varId;
} // getVar

//====================== Other Virtual Member Functions =======================

// Symmetric KL distance. Several terms cancel due to the symmetricising.
// 2*symKL = \sum_k( (a_k-b_k)[dg(a_k)-dg(b_k)-dg(\sum_k a_k)+dg(\sum_k b_k)])
// if alpha_k == 0 and beta_k != 0 or vice versa, the distance is infinity

template <typename AssignType>
double Dirichlet<AssignType>::distance(const Factor* rhsFPtr) const {
  // rhsFPtr must be an appropriate Dirichlet*
  const Dirichlet<AssignType>* dwnPtr = dynamic_cast<const Dirichlet<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "Dirichlet::distance : rhs not a valid Dirichlet : " << typeid(*rhsFPtr).name() );
  // Create references to Dirichlet instances
  const Dirichlet<AssignType>& rhs(*dwnPtr);
  // Ensures the number of parameters of the two Dirichlet potential functions are equal
  PRLITE_ASSERT(domain->size() == rhs.domain->size(),
         "Dirichlet_InplaceAbsorb : lhs size " << domain->size() <<
         " not equal to rhs size " << rhs.domain->size() <<
         "\nlhs vars " << getVars() << " rhs vars " << rhs.getVars() );

  // --- Variables we'll need for calculating the distance

  const double& defL(defAlpha);                             // convenient alias
  const double& defR(rhs.defAlpha);                         // convenient alias
  const std::map<AssignType, double>& mapL(sparseAlphas);     // convenient alias
  const std::map<AssignType, double>& mapR(rhs.sparseAlphas); // convenient alias
  double dgDefL = defL > 0.0 ? digamma(defL) : std::numeric_limits<double>::infinity();
  double dgDefR = defR > 0.0 ? digamma(defR) : std::numeric_limits<double>::infinity();
  double dgSumAlphaL = digamma(sumOfAlphas);
  double dgSumAlphaR = digamma(dwnPtr->sumAlphas());
  double diffDefL = dgDefL - dgSumAlphaL;
  double diffDefR = dgDefR - dgSumAlphaR;
  size_t nDefsL = domain->size()-mapL.size();
  size_t nDefsR = rhs.domain->size()-mapR.size();

  double dist = 0.0;

  // --- Lets start iterating

  auto iL = mapL.begin();
  auto iR = mapR.begin();
  while ( iL != mapL.end() && iR != mapR.end() ){

    if (iL->first < iR->first) {      // mapL,defR
      nDefsR--;
      double ak = iL->second;

      if (ak != defR){
        if (ak > 0.0 and defR > 0.0) {
          dist += (ak-defR)*(digamma(ak) - dgSumAlphaL - diffDefR);
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      iL++;
    } // if

    else if (iL->first > iR->first) { // defL,mapR
      nDefsL--;
      double bk = iR->second;

      if (bk != defL){
        if (bk > 0.0 and defL > 0.0) {
          dist += (defL-bk)*(diffDefL - (digamma(bk) -dgSumAlphaR) );
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      iR++;
    } // else if

    else {                            //mapL,mapR
      double ak = iL->second;
      double bk = iR->second;

      if (ak != bk){
        if (ak > 0.0 and bk > 0.0) {
          dist += (ak-bk)*(digamma(ak) - dgSumAlphaL - (digamma(bk) - dgSumAlphaR) );
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      iL++; iR++;
    } // else

  } // while

  while ( iL != mapL.end() ){         // some extra alphas in mapL,defR
    nDefsR--;
    double ak = iL->second;

    if (ak != defR){
      if (ak > 0.0 and defR > 0.0) {
        dist += (ak-defR)*(digamma(ak) - dgSumAlphaL - diffDefR);
      } // if
      else {return std::numeric_limits<double>::infinity();} // else
    } // if

    iL++;
  } // while

  while ( iR != mapR.end() ){         // some extra betas in mapR,defL
    nDefsL--;
    double bk = iR->second;

    if (bk != defL){
      if (bk > 0.0 and defL > 0.0) {
        dist += (defL-bk)*(diffDefL - (digamma(bk) -dgSumAlphaR) );
      } // if
      else {return std::numeric_limits<double>::infinity();} // else
    } // if

    iR++;
  } // while

  // now we may still have a number of alphas and betas in neither of the maps defL,defR
  PRLITE_ASSERT(
    nDefsL == nDefsR,
    "Dirichlet::distance : The number of cases in neither map should match, no matter from which side you look at it: " << nDefsL << "!=" << nDefsR);

  if (defL != defR){
    if (defL > 0.0 and defR > 0.0) {
      dist += nDefsL*(defL-defR)*(diffDefL - diffDefR);
    } // if
    else {return std::numeric_limits<double>::infinity();} // else
  } // if

  return dist/2.0;
} // distance

//======================= Non-virtual Member Functions ========================

template <typename AssignType>
Dirichlet<AssignType>* Dirichlet<AssignType>::sliceAt(const AssignType& theObservedWord) const {

  double newDefault = (sumAlphas() - getAlpha(theObservedWord) ) / (domain->size()-1);
  return new Dirichlet<AssignType>(
    varId,
    domain,
    newDefault,
    { {theObservedWord, getAlpha(theObservedWord)} },
    relMargin,
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

} // sliceAt

template <typename AssignType>
void Dirichlet<AssignType>::makeDefaultZero() const {

  if (defAlpha > 0.0) {

    // iterate over all possible assignments for the newVars
    for (auto val : *domain) {
      if ( sparseAlphas.find(val) == sparseAlphas.end() ){
        sparseAlphas[val] = defAlpha;
      } // if
    } // for

    defAlpha = 0.0;
  } // if

} // makeDefaultZero

template <typename AssignType>
size_t Dirichlet<AssignType>::noOfCategories() const {
  return domain->size();
} // noOfCategories

template <typename AssignType>
void Dirichlet<AssignType>::addToCounts(const AssignType& k,
                                        double cnt){
  auto itr = sparseAlphas.find(k);
  if ( itr == sparseAlphas.end() ) {sparseAlphas[k] = defAlpha + cnt;} // if
  else {itr->second += cnt;} // else

  sumOfAlphas += cnt;
} // addToCounts

template <typename AssignType>
double Dirichlet<AssignType>::getDefaultAlpha() const {
  return defAlpha;
} // getDefaultAlpha

template <typename AssignType>
const std::map<AssignType, double>& Dirichlet<AssignType>::getSparseAlphas() const {
  return sparseAlphas;
} // getSparseAlphas


template <typename AssignType>
double Dirichlet<AssignType>::getAlpha(const AssignType& k) const {
  auto itr = sparseAlphas.find(k);
  if ( itr == sparseAlphas.end() ) {return defAlpha;} // if
  else {return itr->second;} // else
} // getAlpha

template <typename AssignType>
void Dirichlet<AssignType>::copyExtraSettingsFrom(const Dirichlet<AssignType>& d) {
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

  // rescale the alphas to sum to a specified value
template <typename AssignType>
void Dirichlet<AssignType>::rescaleSumOfAlphasTo(double newSum){
  double scale = newSum/sumOfAlphas;
  defAlpha *= scale;
  for (auto& el : sparseAlphas) {
    el.second *= scale;
  } // for
  sumOfAlphas = newSum;
} // rescaleSumOfAlphasTo

  // rescale the alphas to not exceed this value
template <typename AssignType>
void Dirichlet<AssignType>::limitSumOfAlphasTo(double noMoreThan) {
  if (sumOfAlphas > noMoreThan) {
    rescaleSumOfAlphasTo(noMoreThan);
    sumOfAlphas = noMoreThan;
  } // if
} // limitSumOfAlphasTo


// sumOfAlphas
template <typename AssignType>
double Dirichlet<AssignType>::sumAlphas() const {
  return sumOfAlphas;
} // sumAlphas

template <typename AssignType>
void Dirichlet<AssignType>::recalcSumOfAlphas() const {
  sumOfAlphas =
    std::accumulate(sparseAlphas.begin(), sparseAlphas.end(), 0.0,
                    [](double sofar, const std::pair<AssignType,double>& p) { return sofar+p.second; });
  size_t nDefs = domain->size() - sparseAlphas.size();
  sumOfAlphas += nDefs*defAlpha;
  PRLITE_ASSERT(sumOfAlphas > 0.0, "Sum of Alphas must be bigger than zero: " << sumOfAlphas);
} // recalcSumOfAlphas

//======================= Factor operators ====================================

// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample ****************************************
// ***************************************************************************

// Identify object
template <typename AssignType>
const std::string& Dirichlet_Sampler<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_Sampler");
  return CLASSNAME;
} // isA

// Samples a random vector from the distribution
template <typename AssignType>
emdw::RVVals Dirichlet_Sampler<AssignType>::sample(const Dirichlet<AssignType>* lhsPtr) {

  double sum = 0.0;
  std::vector<double> randVals; randVals.reserve( lhsPtr->domain->size() );
  for ( auto ass : *(lhsPtr->domain) ) {

    // now find the appropriate alpha value
    double alpha = lhsPtr->defAlpha;
    auto itr = lhsPtr->sparseAlphas.find(ass);
    if ( itr != lhsPtr->sparseAlphas.end() ) {
      alpha = itr->second;
    } // if

    std::gamma_distribution<double> randGamma(alpha,1.0);
    double randVal = randGamma( emdw::randomEngine() );
    randVals.push_back(randVal);
    sum += randVal;
  } // for

  // and convert it into probabilities
  emdw::RVVals rvVals;
  rvVals.reserve(randVals.size() );
  for (auto el : randVals) {
    rvVals.push_back(el/sum);
  } // for

  return rvVals;
} // sample

/*****************************************************************************/
/***************************** Family 1 **************************************/
/***************************** Normalize *************************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Dirichlet_InplaceNormalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_InplaceNormalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// The Dirichlet should be implicitly normalized. However, we can use
// the opportunity to do some cleaning up to make things more sparse
// if possible.

template <typename AssignType>
void Dirichlet_InplaceNormalize<AssignType>::inplaceProcess(Dirichlet<AssignType>* lhsPtr) {

  bool specVaries = false;
  double prevVal =0.0;
  double& sumOfAlphas = lhsPtr->sumOfAlphas;
  double& defAlpha = lhsPtr->defAlpha;

  // find the max alpha
  double maxVal = defAlpha;
  for (auto itr = lhsPtr->sparseAlphas.begin(); itr != lhsPtr->sparseAlphas.end(); itr++) {
    if (itr->second > maxVal) {maxVal = itr->second;} // if
  } // for
  double theThresh = lhsPtr->relMargin*maxVal;

  // zero all alphas less than theThresh
  if (defAlpha < theThresh) {
    sumOfAlphas -= defAlpha * ( lhsPtr->domain->size()-lhsPtr->sparseAlphas.size() );
    defAlpha = 0.0;
  } // if

  auto itr = lhsPtr->sparseAlphas.begin();
  while ( itr != lhsPtr->sparseAlphas.end() ) {
    double val = itr->second;
    if (val < theThresh) {
      sumOfAlphas -= val;
      itr->second = val = 0.0;
    } // if

    // check for variation in the spec values
    if (itr != lhsPtr->sparseAlphas.begin() and std::abs(val-prevVal) > 0.0) {
      specVaries = true;
    } // if
    prevVal = val;

    // remove any values equal to the default, and increment iterator
    if (std::abs(itr->second-defAlpha) > 0.0) {
      itr++;
    } // if
    else {
      itr = lhsPtr->sparseAlphas.erase(itr); // c++11 returns the next itr
    } // else

  } // while

  if ( lhsPtr->domain->size() == lhsPtr->sparseAlphas.size() ) { // fully dense, free to change default value

    if (specVaries) { // 0.0 probably is a more sensible default
      defAlpha = 0.0;
    } // if

    else { // fully dense with no variation - bonanza
      defAlpha = prevVal;
      lhsPtr->sparseAlphas.clear();
    } // else

  } // if

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Dirichlet_Normalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_Normalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Normalizing the dirichlet distribution returns the original distribution
// since there is only a single random variable theta
template <typename AssignType>
Factor* Dirichlet_Normalize<AssignType>::process(const Dirichlet<AssignType>* lhsPtr) {
  // make a copy of the factor
  Dirichlet<AssignType>* fPtr = new Dirichlet<AssignType>(*lhsPtr);

  // inplace normalize it
  Dirichlet_InplaceNormalize<AssignType> ipNorm;
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
template <typename AssignType>
const std::string& Dirichlet_InplaceAbsorb<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_InplaceAbsorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA


// One Dirichlet absorbs another - this is implemented by adding the
// parameters of the given Dirichlets. Algebraically the product will
// no longer be normalised. However, since we do not explicitly keep
// track of the volume, we are going to assume that it is
// re-normalised automatically. Might need to revisit this.

// When one multiply two true Dirichlets alpha' = alpha1 + alpha2
// -1. However, as far as I can make out we will never multiply two
// true Dirichlets - instead we will be multiplying a true dirichlet
// (lhs) with a Dirichlet conjugate (determined by the Polya
// neighbour's marginalising). In this case alpha' = alpha1 +
// alpha2. We will implement this here.
template <typename AssignType>
void
Dirichlet_InplaceAbsorb<AssignType>::inplaceProcess(
  Dirichlet<AssignType>* lhsPtr,
  const Factor* rhsFPtr) {

  // rhsFPtr must be an appropriate Dirichlet*
  const Dirichlet<AssignType>* dwnPtr = dynamic_cast<const Dirichlet<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "Dirichlet_InplaceAbsorb : rhs not a valid Dirichlet : " << typeid(*dwnPtr).name() );
  // Create references to Dirichlet instances
  const Dirichlet<AssignType>& rhs(*dwnPtr);
  Dirichlet<AssignType>& lhs(*lhsPtr);

  std::map<AssignType, double>& mapL(lhs.sparseAlphas);         // convenient alias
  double& defL(lhs.defAlpha);                                 // convenient alias
  double& sumOfAlphas = lhs.sumOfAlphas;                      // convenient alias

  const std::map<AssignType, double>& mapR(rhs.sparseAlphas); // convenient alias
  const double& defR(rhs.defAlpha);                         // convenient alias

  // Ensures the number of parameters of the two Dirichlet potential functions are equal
  PRLITE_ASSERT(lhs.domain->size() == rhs.domain->size(),
         "Dirichlet_InplaceAbsorb : lhs size " << lhs.domain->size() <<
         " not equal to rhs size " << rhs.domain->size() <<
         "\nlhs vars " << lhs.getVars() << " rhs vars " << rhs.getVars() );

  // efficient handling of relatively sparse rhs with defAlpha == 0)
  if (defR == 0.0 and mapR.size() <= mapL.size()/10){
    for (auto r : mapR) {
      mapL[r.first] += r.second-1.00;
      sumOfAlphas += r.second-1.00;
    } // for
    return;
  } // if

  auto iL = mapL.begin();
  auto iR = mapR.begin();
  while ( iL != mapL.end() && iR != mapR.end() ){

    if (iL->first < iR->first) {      // mapL,defR
      iL->second += defR-1.00;
      sumOfAlphas += defR-1.00;
      iL++;
    } // if

    else if (iL->first > iR->first) { // defL,mapR
      mapL[iR->first] = defL + iR->second-1.00;
      sumOfAlphas += iR->second-1.00;
      iR++;
    } // else if

    else {                            // mapL,mapR
      iL->second += iR->second-1.00;
      sumOfAlphas += iR->second-1.00;
      iL++; iR++;
    } // else

  } // while

  while ( iL != mapL.end() ){         // some extra sparseL's mapL,defR
    iL->second += defR-1.00;
    sumOfAlphas += defR-1.00;
    iL++;
  } // while

  while ( iR != mapR.end() ){         // some extra sparseR's in defL,mapR
    mapL[iR->first] = defL + iR->second-1.00;
    sumOfAlphas += iR->second-1.00;
    iR++;
  } // while

  // mapL.size() is its new updated size. Subtracting that from
  // lhs.domain->size() gives the remaining number of unvisited defaults.
  sumOfAlphas += (defR-1.00)*(lhs.domain->size() - mapL.size() );
  defL += defR-1.00;

} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Dirichlet_Absorb<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_Absorb<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One Dirichlet absorbs another - this is implemented by adding the parameters
// of the given Dirichlets
template <typename AssignType>
Factor* Dirichlet_Absorb<AssignType>::process(const Dirichlet<AssignType>* lhsPtr, const Factor* rhsFPtr) {
  // make a copy of the factor
  Dirichlet<AssignType>* fPtr = new Dirichlet<AssignType>(*lhsPtr);

  // inplace absorb it
  Dirichlet_InplaceAbsorb<AssignType> ipAbsorb;
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

/*****************************************************************************/
/***************************** Family 2 **************************************/
/***************************** Cancel ****************************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Dirichlet_InplaceCancel<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_InplaceCancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA


// One Dirichlet cancels another - this is implemented by subtracting
// the parameter of the one potential from that of the
// other. Algebraically the quotient will no longer be
// normalised. However, since we do not explicitly keep track of the
// volume, we are going to assume that it is re-normalised
// automatically. Might need to revisit this.

// Check the note above on multiplying a true Dirichlet with a
// Polya-marginalised conjugate Dirichlet.
template <typename AssignType>
void
Dirichlet_InplaceCancel<AssignType>::inplaceProcess(
  Dirichlet<AssignType>* lhsPtr,
  const Factor* rhsFPtr) {

  // rhsFPtr must be an appropriate Dirichlet*
  const Dirichlet<AssignType>* dwnPtr = dynamic_cast<const Dirichlet<AssignType>*>(rhsFPtr);
  PRLITE_ASSERT(dwnPtr,
         "Dirichlet_InplaceCancel : rhs not a valid Dirichlet : " << typeid(*dwnPtr).name() );
  // Create references to Dirichlet instances
  const Dirichlet<AssignType>& rhs(*dwnPtr);
  Dirichlet<AssignType>& lhs(*lhsPtr);

  std::map<AssignType, double>& mapL(lhs.sparseAlphas);         // convenient alias
  double& defL(lhs.defAlpha);                                 // convenient alias
  double& sumOfAlphas = lhs.sumOfAlphas;                      // convenient alias

  const std::map<AssignType, double>& mapR(rhs.sparseAlphas); // convenient alias
  const double& defR(rhs.defAlpha);                         // convenient alias

  // Ensures the number of parameters of the two Dirichlet potential functions are equal
  PRLITE_ASSERT(lhs.domain->size() == rhs.domain->size(),
         "Dirichlet_InplaceCancel : The size of the Dirichlet parameters differ");

  // efficient handling of relatively sparse rhs with defAlpha == 0)
  if (defR == 0.0 and mapR.size() <= mapL.size()/10){
    for (auto r : mapR) {
      mapL[r.first] -= (r.second-1.00);
      sumOfAlphas -= (r.second-1.00);
    } // for
    return;
  } // if

  auto iL = mapL.begin();
  auto iR = mapR.begin();
  while ( iL != mapL.end() && iR != mapR.end() ){

    if (iL->first < iR->first) {      // only in mapL
      iL->second -= (defR-1.00);
      sumOfAlphas -= (defR-1.00);
      iL++;
    } // if

    else if (iL->first > iR->first) { // only in mapR
      mapL[iR->first] = defL - (iR->second-1.00);
      sumOfAlphas -= (iR->second-1.00);
      iR++;
    } // else if

    else {                            //common to both maps
      iL->second -= (iR->second-1.00);
      sumOfAlphas -= (iR->second-1.00);
      iL++; iR++;
    } // else

  } // while

  while ( iL != mapL.end() ){         // some extra spec's in mapL
    iL->second -= (defR-1.00);
    sumOfAlphas -= (defR-1.00);
    iL++;
  } // while

  while ( iR != mapR.end() ){         // some extra spec's in mapR
    mapL[iR->first] = defL - (iR->second-1.00);
    sumOfAlphas -= (iR->second-1.00);
    iR++;
  } // while

  // mapL.size() is its new updated size. Subtracting that from
  // lhs.domain->size() gives the remaining number of unvisited defaults.
  sumOfAlphas -= (defR-1.00)*(lhs.domain->size() - mapL.size() );
  defL -= (defR-1.00);
} // inplaceProcess

// Identify object
template <typename AssignType>
const std::string& Dirichlet_Cancel<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_Cancel<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// One Dirichlet cancels another - this is implemented by subtracting the parameter
// of the one potential from that of the other
template <typename AssignType>
Factor* Dirichlet_Cancel<AssignType>::process(const Dirichlet<AssignType>* lhsPtr, const Factor* rhsFPtr) {
  // make a copy of the factor
  Dirichlet<AssignType>* fPtr = new Dirichlet<AssignType>(*lhsPtr);

  // inplace normalize it
  Dirichlet_InplaceCancel<AssignType> ipCancel;
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
template <typename AssignType>
const std::string& Dirichlet_Marginalize<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_Marginalize<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Given Dirichlet is already marginalized in its parametric form
// This function therefore merely returns a copy of the given potential
template <typename AssignType>
Factor* Dirichlet_Marginalize<AssignType>::process(const Dirichlet<AssignType>* lhsPtr,
                                                   const emdw::RVIds& toKeep,
                                                   bool presorted,
                                                   const Factor* peekAheadPtr) {
  // make a copy of the factor
  PRLITE_ASSERT(toKeep == lhsPtr->getVars(),
         "Dirichlet_Marginalize : Dirichlet marginalization should not remove anything");
  return lhsPtr->copy();
} // process

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/

// Identify object
template <typename AssignType>
const std::string& Dirichlet_InplaceWeakDamping<AssignType>::isA() const {
  static const std::string CLASSNAME("Dirichlet_InplaceWeakDamping<"+std::string(typeid(AssignType).name())+std::string(">"));
  return CLASSNAME;
} // isA

// Dampens new message by taking a linear combination of the new and
// old alpha counts - this linear combination is determined by the
// df parameter
template <typename AssignType>
double
Dirichlet_InplaceWeakDamping<AssignType>::inplaceProcess(
  Dirichlet<AssignType>* lhsPtr,
  const Factor* rhsFPtr,
  double df) { // df is the weight for rhs


  //----------------------------------------------------
  // rhsFPtr must be an appropriate Dirichlet*
  const Dirichlet<AssignType>* dwnPtr = dynamic_cast<const Dirichlet<AssignType>*>(rhsFPtr);

  PRLITE_ASSERT(dwnPtr,
         "Dirichlet::distance : rhs not a valid Dirichlet : " << typeid(*rhsFPtr).name() );
  // Create references to Dirichlet instances
  const Dirichlet<AssignType>& rhs(*dwnPtr);
  Dirichlet<AssignType>& lhs(*lhsPtr);

  // Ensures the number of parameters of the two Dirichlet potential functions are equal
  PRLITE_ASSERT(lhs.domain->size() == rhs.domain->size(),
         "Dirichlet_InplaceAbsorb : lhs size " << lhs.domain->size() <<
         " not equal to rhs size " << rhs.domain->size() <<
         "\nlhs vars " << lhs.getVars() << " rhs vars " << rhs.getVars() );

  double dfCompl = 1.0 - df;

  // --- Variables we'll need for calculating the distance

  double& defL(lhs.defAlpha);                               // convenient alias
  const double& defR(rhs.defAlpha);                         // convenient alias
  std::map<AssignType, double>& mapL(lhs.sparseAlphas);       // convenient alias
  const std::map<AssignType, double>& mapR(rhs.sparseAlphas); // convenient alias
  double dgDefL = defL > 0.0 ? digamma(defL) : std::numeric_limits<double>::infinity();
  double dgDefR = defR > 0.0 ? digamma(defR) : std::numeric_limits<double>::infinity();
  double dgSumAlphaL = digamma(lhs.sumAlphas());
  double dgSumAlphaR = digamma(dwnPtr->sumAlphas());
  double diffDefL = dgDefL - dgSumAlphaL;
  double diffDefR = dgDefR - dgSumAlphaR;
  size_t nDefsL = lhs.domain->size()-mapL.size();
  size_t nDefsR = rhs.domain->size()-mapR.size();

  double dist = 0.0;

  // --- Lets start iterating

  auto iL = mapL.begin();
  auto iR = mapR.begin();
  while ( iL != mapL.end() && iR != mapR.end() ){

    if (iL->first < iR->first) {      // mapL,defR
      nDefsR--;
      double& ak = iL->second;

      if (ak != defR){
        if (ak > 0.0 and defR > 0.0) {
          dist += (ak-defR)*(digamma(ak) - dgSumAlphaL - diffDefR);
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      ak = dfCompl*ak + df*defR;
      iL++;
    } // if

    else if (iL->first > iR->first) { // defL,mapR
      nDefsL--;
      double bk = iR->second;

      if (bk != defL){
        if (bk > 0.0 and defL > 0.0) {
          dist += (defL-bk)*(diffDefL - (digamma(bk) -dgSumAlphaR) );
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      mapL[iR->first] = dfCompl*defL + df*bk;
      iR++;
    } // else if

    else {                            //mapL,mapR
      double& ak = iL->second;
      double bk = iR->second;

      if (ak != bk){
        if (ak > 0.0 and bk > 0.0) {
          dist += (ak-bk)*(digamma(ak) - dgSumAlphaL - (digamma(bk) - dgSumAlphaR) );
        } // if
        else {return std::numeric_limits<double>::infinity();} // else
      } // if

      ak = dfCompl*ak + df*bk;
      iL++; iR++;
    } // else

  } // while

  while ( iL != mapL.end() ){         // some extra alphas in mapL,defR
    nDefsR--;
    double& ak = iL->second;

    if (ak != defR){
      if (ak > 0.0 and defR > 0.0) {
        dist += (ak-defR)*(digamma(ak) - dgSumAlphaL - diffDefR);
      } // if
      else {return std::numeric_limits<double>::infinity();} // else
    } // if

    ak = dfCompl*ak + df*defR;
    iL++;
  } // while

  while ( iR != mapR.end() ){         // some extra betas in mapR,defL
    nDefsL--;
    double bk = iR->second;

    if (bk != defL){
      if (bk > 0.0 and defL > 0.0) {
        dist += (defL-bk)*(diffDefL - (digamma(bk) -dgSumAlphaR) );
      } // if
      else {return std::numeric_limits<double>::infinity();} // else
    } // if

    mapL[iR->first] = dfCompl*defL + df*bk;
    iR++;
  } // while

  // now we may still have a number of alphas and betas in neither of the maps defL,defR
  PRLITE_ASSERT(
    nDefsL == nDefsR,
    "Dirichlet::distance : The number of cases in neither map should match, no matter from which side you look at it: " << nDefsL << "!=" << nDefsR);

  if (defL != defR){
    if (defL > 0.0 and defR > 0.0) {
      dist += nDefsL*(defL-defR)*(diffDefL - diffDefR);
    } // if
    else {return std::numeric_limits<double>::infinity();} // else
  } // if

  defL = dfCompl*defL + df*defR;
  lhs.recalcSumOfAlphas();
  return dist/2.0;

} // inplaceProcess
