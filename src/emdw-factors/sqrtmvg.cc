/*******************************************************************************


          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Multi-variate SqrtMVG PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// gLinear headers
//#include "gLapackHeaders.h"
#include "gLinear/gLapackWrap.h"

// patrec headers
#include "sortindices.hpp"
#include "prlite_vector.hpp"
#include "prlite_stlvecs.hpp"
#include "prlite_stdfun.hpp"  // TestCase
#include "prlite_logging.hpp"
#include "prlite_testing.hpp" //AlmostEqual


//emdw headers
#include "emdw.hpp"
#include "sqrtmvg.hpp"
#include "vecset.hpp"
#include "matops.hpp" // subtractOuterProd
#include "lapackblas.hpp"

// standard headers
#include <new>  // placement new
#include <math.h>
#include <random>

using namespace std;
using namespace emdw;

const double HALFLOG2PI = log(2.0*M_PI)/2.0;

/**************************The SqrtMVG class members***********************/

//============================ Traits and Typedefs ============================

//======================== Constructors and Destructor ========================

// The set of default factor operators are defined in this helper function
void SqrtMVG::supplementWithDefaultOperators() {

  // provide the default version if still unspecified
  if (!marginalizer) {
    marginalizer =
      uniqptr<Marginalizer>( new SqrtMVG_Marginalize() );
  } // if

  if (!inplaceNormalizer) {
    inplaceNormalizer =
      uniqptr<InplaceNormalizer>( new SqrtMVG_InplaceNormalize() );
  } // if

  if (!normalizer) {
    normalizer =
      uniqptr<Normalizer>( new SqrtMVG_Normalize() );
  } // if

  if (!inplaceAbsorber) {
    inplaceAbsorber =
      uniqptr<InplaceAbsorber>( new SqrtMVG_InplaceAbsorb() );
  } // if

  if (!absorber) {
    absorber =
      uniqptr<Absorber>( new SqrtMVG_Absorb() );
  } // if

  if (!inplaceCanceller) {
    inplaceCanceller =
      uniqptr<InplaceCanceller>( new SqrtMVG_InplaceCancel() );
  } // if

  if (!canceller) {
    canceller =
      uniqptr<Canceller>( new SqrtMVG_Cancel() );
  } // if

  if (!observeAndReducer) {
    observeAndReducer =
      uniqptr<ObserveAndReducer>( new SqrtMVG_ObserveAndReduce() );
  } // if

  if (!inplaceDamper) {
    inplaceDamper = uniqptr<InplaceDamper>( new SqrtMVG_InplaceNoDamping() );
  } // if

  if (!sampler) {
    sampler =
      uniqptr<Sampler>( new SqrtMVG_Sampler() );
  } // if

} // supplementWithDefaultOperators

// default constructor
SqrtMVG::SqrtMVG(
  const emdw::RVIds& theVars,
  bool theAutoNormed,
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
    : vars( theVars.size() ),
      vacuous(true), //an empty gaussian is vacuous until values are set
      autoNormed(theAutoNormed),
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

  if ( theVars.size() ) {

    if (presorted or theVars.size() == 1) {
      vars = theVars;
    } // if

    else {
      std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );
      vars = extract<unsigned>(theVars, sorted);
    } // else

  } // if

} // default constructor

// copy constructor. We provide this to make sure we do not have gLinear deep/shallow issues
SqrtMVG::SqrtMVG(const SqrtMVG& p)
    : vars(p.vars),
      vacuous(p.vacuous),
      L( p.L.rows(), p.L.cols() ),
      h( p.h.size() ),
      g(p.g),
      autoNormed(p.autoNormed),
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
  L = p.L;
  h = p.h;
  // NOTE: We are assuming the operators are without memory, i.e. we
  // can apply the same operator to several factor instances without
  // any harm.
} // cpy ctor

// class specific constructor Triangular CanonicalForm
SqrtMVG::SqrtMVG(
  const emdw::RVIds& theVars,
  const prlite::ColMatrix<double>& theL,
  const prlite::ColVector<double>& theH,
  const double& theG,
  bool theAutoNormed,
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
: vars(theVars.size()),
  vacuous(false),
  L( theL.rows(),theL.cols() ),
  h( theH.size() ),
  g(theG),
  autoNormed(theAutoNormed),
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

  size_t nDims = vars.size();
  if (nDims == 1) {presorted = true;}

  supplementWithDefaultOperators();

#ifdef CHK_IDX
  // This checks that the dimension of theVars and other parameters match
  if (!vacuous) {
    PRLITE_ASSERT(nDims == (unsigned)theL.rows(),
           nDims << " variables, but " <<
           theL.rows() << " L Matrix row(s)");
    PRLITE_ASSERT(nDims == (unsigned)theL.cols(),
           nDims << " variables, but " <<
           theL.cols() << " L Matrix col(s)");
    PRLITE_ASSERT(nDims == (unsigned)theH.size(),
           nDims << " variables, but " <<
           theH.size() << " h col vector row(s)");
  } // if
#endif // CHK_IDX

  //if presorted, just save values
  if (presorted) {
    vars = theVars;
    L = theL;
    h = theH;
  } // if

  else{

    std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );
    vars = extract<unsigned>(theVars, sorted);

    //map from unsorted to sorted
    for (size_t r=0; r < nDims; r++){
      h[r] = theH[sorted[r]];
      L.row(r) = theL.row(sorted[r]); // no longer lower triangular
    } // for r
    rotateColsToLowerTriangular(L);   // therefore we fix that

  } // else

  // check posdef
  for (size_t n = 0; n < nDims; n++) {
    PRLITE_ASSERT(L(n,n) > 0, "Can not construct from a non-posdef L: " << L);
  } // for

} // class specific constructor Triangular CanonicalForm


SqrtMVG::~SqrtMVG() {
  // --noOfObjects;
  // std::cout << '*' << noOfObjects << '*' << std::flush;
} // destructor

// class specific constructor covariance form
SqrtMVG::SqrtMVG(
  const emdw::RVIds& theVars,
  const prlite::ColVector<double>& theMn,
  const prlite::RowMatrix<double>& theCov,
  bool theAutoNormed,
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
    : vars(theVars.size()),
      vacuous(false),
      L( theVars.size(), theVars.size() ),
      h(theVars.size()),
      autoNormed(theAutoNormed),
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

  int nDims = vars.size(); // deliberately int for lapack's sake
  if (nDims == 1) {presorted = true;}

  supplementWithDefaultOperators();
  prlite::ColVector<double> sortedMn;
  if (!autoNormed and !presorted) {
    sortedMn.resize(nDims);
  } // if

#ifdef CHK_IDX
  // This checks that the dimension of theVars and other parameters match
  PRLITE_ASSERT(nDims == theCov.rows(),
         nDims << " variables, but " <<
         theCov.rows() << " cov Matrix row(s)");
  PRLITE_ASSERT(nDims == theCov.cols(),
         nDims << " variables, but " <<
         theCov.cols() << " cov Matrix col(s)");
  PRLITE_ASSERT(nDims == theMn.size(),
         nDims << " variables, but " <<
         theMn.size() << " mn col vector row(s)");
#endif // CHK_IDX

  //if presorted, just save values. Will update them lower down.
  if (presorted) {
    vars = theVars;
    L = theCov;
    h = theMn;
  } // if

  else {
    std::vector<size_t> sorted = sortIndices( theVars, std::less<unsigned>() );
    vars = extract<unsigned>(theVars, sorted);

    //map from unsorted to sorted
    for (int r=0; r < theCov.rows(); r++){

      h[r] = theMn[ sorted[r] ];
      if (!autoNormed) {
        sortedMn[r] = h[r];
      } // if

      for (int c=0; c <= r; c++){  // we only need to do the lower triangle
        L(r,c) = theCov(sorted[r],sorted[c]); // this only works on symmetric matrices!
      } // for c

    } // for r

  } // else

  // --------------------Determine L so that inv(cv) = K = LL'

  // 'n Alternatief vir die prosedure hieronder sou wees om C = B*B^T
  // = \sum_n(C^{(N)})C^{(N)T} te sien en van cholUL's gebruik te
  // maak. Vermoed egter dat dit nie regtig voordele inhou nie.

  char uplo = 'L';
  int lda = nDims;
  int fail;
  // B so that BB' = cv
  dpotrf_(&uplo, &nDims, &L(0,0), &lda, &fail ); // cholesky factorize
  PRLITE_ASSERT(!fail, "Could not get the cholesky root of the covariance");

  // K = inv(cov). Let op as ons bloot L inverteer: K = inv(LL') =
  // L^{-T}L^{-1} dws upper times lower vorm. Geldig, maar
  // ongewoon. Hier vat ek die effens duurder pad deur eers K te kry,
  // en dan sy wortels.
  // K = inv(cv)
  dpotri_(&uplo, &nDims, &L(0,0), &lda, &fail ); //invert symm posdef
  PRLITE_ASSERT(!fail, "Could not invert the covariance matrix");

  //L so that K = LL'
  dpotrf_(&uplo, &nDims, &L(0,0), &lda, &fail ); // cholesky factorize
  PRLITE_ASSERT(!fail, "Could not get the cholesky root of the precision");
  // In the cholUL/DL routines we use full columns of L as vectors
  zeroUpper(L);
  // if we survived al those tests it should be posdef

  // --------------------Determine h = K * mn

  char trans = 'T';  // we want L'
  char diag = 'N';
  int incx = 1;
  // L' * mn
  dtrmv_(&uplo, &trans, &diag, &nDims, &L(0,0), &lda, &h[0], &incx); // triang * vec
  trans = 'N'; // this one uses L
  // h = LL' * mn = K * mn
  dtrmv_(&uplo, &trans, &diag, &nDims, &L(0,0), &lda, &h[0], &incx); // triang * vec

  if (autoNormed) {return;} // if

  // --------------------Determine g

  // faster than getNormedG because we have the mn available

  double sqrtDetK = L(0,0);
  for (int k = 1; k < nDims; k++) {
    sqrtDetK *= L(k,k);
  } // for
  if (presorted) {
    g = -0.5 * (theMn.transpose() * h) - nDims*HALFLOG2PI + log(sqrtDetK);
  } // if
  else{
    g = -0.5 * (sortedMn.transpose() * h) - nDims*HALFLOG2PI + log(sqrtDetK);
  } // else

} // class specific constructor CovarianceForm

//========================= Operators and Conversions =========================

bool SqrtMVG::operator==(const SqrtMVG& p) const {

  double tolerance = std::numeric_limits<double>::epsilon()*100000;

  if (vars != p.vars) {return false;} // if
  if (autoNormed != p.autoNormed){return false;} // if
  if (vacuous != p.vacuous){return false;} // if

  if (!vacuous) {

    if ( !autoNormed and  abs<double>(g - p.g) > tolerance*max( abs<double>(g), abs<double>(p.g) ) ){
      cout << "  g = " << g << endl;
      cout << "p.g = " << p.g << endl;
      return false;
    } // if

    // check h
    if ( h.size() != p.h.size() ) {return false;}
    for (int n = 0; n < h.size(); n++) {
      if ( abs<double>( h[n] - p.h[n]) > tolerance*max( abs<double>(h[n]),abs<double>(p.h[n]) ) ){
        cout << "  h[" << n << "] = " << h[n] << endl;
        cout << "p.h[" << n << "] = " << p.h[n] << endl;
        return false;
      } // if
    } // for n

    // check L
    if ( L.rows() != p.L.rows() or L.cols() != p.L.cols() ) {return false;}
    for (int r = 0; r < L.rows(); r++) {
      for (int c = 0; c <= r; c++) {
        if ( abs<double>( L(r,c) - p.L(r,c) ) > tolerance*( max( abs( L(r,c) ), abs( p.L(r,c) ) ) ) ) {
          //abs<double>( L(r,c) ) ) {
          cout << "  L(" << r << "," << c << ") = " << L(r,c) << endl;
          cout << "p.L(" << r << "," << c << ") = " << p.L(r,c) << endl;
          return false;
        } // if
      } // for c
    } // for r

  } // if ! vacuous

  return true;
} // operator==

bool SqrtMVG::operator!=(const SqrtMVG& p) const {
  return !operator==(p);
} // operator !=

// boiler plate code
bool SqrtMVG::isEqual(const Factor* rhsPtr) const {
  const SqrtMVG* dwnPtr = dynamic_cast<const SqrtMVG*>(rhsPtr);
  if (!dwnPtr) {return false;}
  return operator==(*dwnPtr);
} // isEqual

//================================= Iterators =================================

//=========================== Inplace Configuration ===========================

unsigned SqrtMVG::classSpecificConfigure(
  const emdw::RVIds& theVars,
  bool theAutoNormed,
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

  this->~SqrtMVG(); // Destroy existing

  new(this) SqrtMVG(
    theVars,
    theAutoNormed,
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

  return 1;
} // classSpecificConfigure

// class specific
unsigned SqrtMVG::classSpecificConfigure(
  const emdw::RVIds& theVars,
  const prlite::ColMatrix<double>& theL,
  const prlite::ColVector<double>& theH,
  const double& theG,
  bool theAutoNormed,
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

  this->~SqrtMVG(); // Destroy existing

  new(this) SqrtMVG(
    theVars,
    theL,
    theH,
    theG,
    theAutoNormed,
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
  return 1;
} // classSpecificConfigure

// !!! this still needs fixing to handle the FactorOperators
std::ostream& SqrtMVG::txtWrite(std::ostream& file) const {

  cout << "SqrtMVG::txtWrite : FIXME factor operators will revert to defaults\n";
  file << "Vars: " << vars << "\n";
  file << "vacuous: " << vacuous << "\n";

  if (!vacuous) {
  file << "L: " << L << "\n";
  file << "h: " << h << "\n";
  file << "g: " << g << "\n";
  file << "autoNormed: " << autoNormed << "\n";
  } // if

  return file;
} // txtWrite

// !!! this still needs fixing to handle the FactorOperators
std::istream& SqrtMVG::txtRead(std::istream& file) {

  cout << "SqrtMVG::txtRead : FIXME factor operators will revert to defaults\n";
  std::string dummy;
  emdw::RVIds theVars;
  bool theVacuous;
  prlite::ColMatrix<double> theL;
  prlite::ColVector<double> theH;
  double theG;
  bool theAutoNormed;

  file >> dummy;
  PRLITE_CHECK(file >> theVars, IOError, "Error reading vars RowVector.");

  file >> dummy;
  PRLITE_CHECK(file >> theVacuous, IOError, "Error reading vacuous bool.");

  if (theVacuous) { // do an inplace construction of a vacuous SqrtMVG
    this->~SqrtMVG(); // Destroy existing
    new(this) SqrtMVG(theVars); // and do an inplace construction
    return file;
  } // if

  file >> dummy;
  PRLITE_CHECK(file >> theL, IOError, "Error reading L matrix.");

  file >> dummy;
  PRLITE_CHECK(file >> theH, IOError, "Error reading h vector.");

  file >> dummy;
  PRLITE_CHECK(file >> theG, IOError, "Error reading g value.");

  file >> dummy;
  PRLITE_CHECK(file >> theAutoNormed, IOError, "Error reading autoNormed value.");

  // BUG losing all the factor operators
  classSpecificConfigure(theVars, theL, theH, theG, theAutoNormed);
  return file;
} // txtRead

//===================== Required Factor Operations ============================

//------------------Family 0

// Boiler plate code
// returning a sample from the distribution
emdw::RVVals SqrtMVG::sample(Sampler* procPtr) const{
  if (procPtr) return dynamicSample(procPtr, this);
  else return dynamicSample(sampler.get(), this);
} //sample

//------------------Family 1

// Boiler plate code
// returning new factor
uniqptr<Factor> SqrtMVG::normalize(FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this) );
  else return uniqptr<Factor>( dynamicApply(normalizer.get(), this) );
} // normalize

// Boiler plate code
// inplace
void SqrtMVG::inplaceNormalize(FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this);
  else dynamicInplaceApply(inplaceNormalizer.get(), this);
} //inplaceNormalize

//------------------Family 2

// Boiler plate code
// returning a new factor
uniqptr<Factor> SqrtMVG::absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(absorber.get(), this, rhsPtr) );
} //absorb

// Boiler plate code
// inplace
void SqrtMVG::inplaceAbsorb(const Factor* rhsPtr,
                            FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceAbsorber.get(), this, rhsPtr);
} // inplaceAbsorb

// Boiler plate code
// returning a new factor
uniqptr<Factor> SqrtMVG::cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, rhsPtr) );
  else return uniqptr<Factor>( dynamicApply(canceller.get(), this, rhsPtr) );
} //cancel

// Boiler plate code
// inplace
void SqrtMVG::inplaceCancel(const Factor* rhsPtr,
                            FactorOperator* procPtr){
  if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
  else dynamicInplaceApply(inplaceCanceller.get(), this, rhsPtr);
} // inplaceCancel

//------------------Family 3

// Boiler plate code
// returning a new factor
uniqptr<Factor> SqrtMVG::marginalize(const emdw::RVIds& variablesToKeep,
                                     bool presorted,
                                     const Factor* peekAheadPtr,
                                     FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variablesToKeep, presorted, peekAheadPtr) );
  else return uniqptr<Factor>( dynamicApply(marginalizer.get(), this, variablesToKeep, presorted, peekAheadPtr) );
} //marginalize

//------------------Family 4

// Boiler plate code
// returning a new factor
uniqptr<Factor> SqrtMVG::observeAndReduce(const emdw::RVIds& variables,
                                          const emdw::RVVals& assignedVals,
                                          bool presorted,
                                          FactorOperator* procPtr) const {
  if (procPtr) return uniqptr<Factor>( dynamicApply(procPtr, this, variables, assignedVals, presorted) );
  else return uniqptr<Factor>( dynamicApply(observeAndReducer.get(), this, variables, assignedVals, presorted) );
} //observeAndReduce

//------------------Family 5

// Boiler plate code
// inplace returning a double
double SqrtMVG::inplaceDampen(const Factor* rhsPtr,
                              double df,
                              FactorOperator* procPtr) {
  if (procPtr) return dynamicInplaceApply(procPtr, this, rhsPtr, df);
  else return dynamicInplaceApply(inplaceDamper.get(), this, rhsPtr, df);
} // inplaceDampen

//===================== Required Virtual Member Functions =====================

SqrtMVG* SqrtMVG::copy(
  const emdw::RVIds& newVars,
  bool presorted) const {

  if ( newVars.size() ) { // replace the variables

    if (vacuous) { // vacuous version
      return new SqrtMVG(
        newVars, presorted, autoNormed,
        marginalizer,
        inplaceNormalizer, normalizer,
        inplaceAbsorber, absorber,
        inplaceCanceller, canceller,
        observeAndReducer,
        inplaceDamper,
        sampler);
    } // if

    else { // full version
      return new SqrtMVG(
        newVars, L, h, g, autoNormed, presorted,
        marginalizer,
        inplaceNormalizer, normalizer,
        inplaceAbsorber, absorber,
        inplaceCanceller, canceller,
        observeAndReducer,
        inplaceDamper,
        sampler);
    } // else

  } // if newVars

  // identical copy
  return new SqrtMVG(*this);

} // copy

SqrtMVG*
SqrtMVG::vacuousCopy(
  const emdw::RVIds& selVars,
  bool presorted) const {

  if ( selVars.empty() ) {
    return new SqrtMVG(
      vars, presorted, autoNormed,
      marginalizer,
      inplaceNormalizer, normalizer,
      inplaceAbsorber, absorber,
      inplaceCanceller, canceller,
      observeAndReducer,
      inplaceDamper,
      sampler);
  } // if
  else {
    return new SqrtMVG(
      selVars, presorted, autoNormed,
      marginalizer,
      inplaceNormalizer, normalizer,
      inplaceAbsorber, absorber,
      inplaceCanceller, canceller,
      observeAndReducer,
      inplaceDamper,
      sampler);
  } // else

} // vacuousCopy

unsigned SqrtMVG::noOfVars() const {
  return vars.size();
} // noOfVars

emdw::RVIds SqrtMVG::getVars() const {
  return vars;
} // getVars

emdw::RVIdType SqrtMVG::getVar(unsigned varNo) const {
  return vars[varNo];
} // getVar

//====================== Other Virtual Member Functions =======================

double SqrtMVG::potentialAt(
  const emdw::RVIds& theVars,
  const emdw::RVVals& theirVals,
  bool presorted) const{

  if (vacuous) { return 0.0;} // if

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
  PRLITE_ASSERT(sortedVars == vars, "The variables do not match: " << sortedVars << " vs " << vars);

  prlite::ColVector<double> observedVals( theirVals.size() );
  for (unsigned idx = 0; idx < theirVals.size(); idx++){
    if (presorted) {
      observedVals[idx] = static_cast<double>(theirVals[idx]);
    } //if
    else {
      observedVals[idx] = static_cast<double>(theirVals[ sorted[idx] ]);
    } // else
  } // for

  //   h'*x
  int nDims = noOfVars();
  int incx = 1;
  double logf = ddot_(&nDims, &h[0], &incx, &observedVals[0], &incx);

  // -0.5*(L'*x)*(L'*x)'
  char uplo = 'L';
  char trans = 'T';  // we want L'
  char diag = 'N';
  int lda = nDims;
  // L' * delta
  dtrmv_(&uplo, &trans, &diag, &nDims, &L(0,0), &lda, &observedVals[0], &incx); // triang * vec
  // delta'*delta
  logf -=  0.5*ddot_(&nDims, &observedVals[0], &incx, &observedVals[0], &incx);

  // g
  logf += (autoNormed ? getNormedG(L,h) : g);

  return  std::exp(logf);
} // potentialAt

double SqrtMVG::distance(const Factor* rhsPtr) const {

  const SqrtMVG* dwnPtr = dynamic_cast<const SqrtMVG*>(rhsPtr);
  PRLITE_ASSERT(dwnPtr, "SqrtMVG::distance : rhs is not a SqrtMVG");
  PRLITE_ASSERT(vars == dwnPtr->vars, "The variables do not match\n" << vars << std::endl << dwnPtr->vars);

  if (vacuous and dwnPtr->vacuous) {
    return 0.0;
  } // if

  if (vacuous or dwnPtr->vacuous) {
    return 1E3;
  } // if

  return ( mahalanobis(dwnPtr->getMean()) + dwnPtr->mahalanobis(getMean()) )/sqrt( noOfVars() );
} // distance

//======================= Parameterizations (non-virtual) =====================

bool SqrtMVG::isVacuous() const {
  return vacuous;
} // isVacuous

prlite::ColMatrix<double> SqrtMVG::getL() const {

  if (vacuous) {
    prlite::ColMatrix<double> tmp( noOfVars(), noOfVars() );
    tmp.assignToAll(0.0);
    return tmp;
  } // if

  return L;
} // getL

prlite::ColVector<double> SqrtMVG::getH() const {

  if (vacuous) {
    prlite::ColVector<double> tmp( noOfVars() );
    tmp.assignToAll(0.0);
    return tmp;
  } // if

  return h;
} // getH

double SqrtMVG::getG() const {

  if (vacuous) {return 0.0;} // if
  if (autoNormed) {return getNormedG(L,h);} // if
  return g;
} // getG

/*
With K = L * L^T it follows that C = L^{-T} * L^{-1}.
Let B = L^{-T} (note upper triangular).
Then sigmapoint i = 1:dim is given by x_i = k * B^(i)
with k = sqrt(dim+0.5) and ^(i) denoting column i.
I.e. X^(1:dim) = B*D with D = kI.
Therefore L^{T} * X = D.
Solve for X. X^(0) = mn.
 */
prlite::ColMatrix<double> SqrtMVG::getSigmaPoints() const {

  PRLITE_ASSERT(!vacuous, "Can not determine sigma points from a vacuous density");

  int dim = noOfVars();
  prlite::ColMatrix<double> x(dim, 2*dim+1);
  double k = sqrt(dim+0.5);
  for (int r = 0; r < dim; r++) {
    for (int c = 1; c <= dim; c++) {
      x(r,c) = 0.0;
    } // for c
    x(r,r+1) = k;
  } // for r

  char uplo = 'L';
  char trans = 'T';
  char diag = 'N';
  int orderLHS = dim;
  int nColsRHS = orderLHS;
  int nRowsLHS = orderLHS;
  int nRowsRHS = nRowsLHS;
  int fail;

  // we compute x from col 1:dim, reserving 0 for mean
  dtrtrs_( &uplo, &trans, &diag, &orderLHS, &nColsRHS,
           &L(0,0), &nRowsLHS,
           &x(0,1), &nRowsRHS, &fail); // solves triangular system
  PRLITE_ASSERT(!fail, "Could not solve sigma points from L");

  // col(0) gets the mean
  x.col(0) = getMean();

  // the mean pm x_i crowd
  for (int r = 0; r < dim; r++) {
    for (int c = 1; c <= dim; c++) {
      x(r,c+dim) = x(r,0) - x(r,c); // do first
      x(r,c) += x(r,0);             // now overwrite
    } // for c
  } // for r

  return x;
} // getSigmaPoints

// mn = Ch = K^{-1}h = (LL^T)^{-1} = L^{-T}L^{-1}h
// Stap 1: Bepaal eers L^{-1}h = b dws Lb = h, gebruik dtrtrs om b te vind.
// Stap 2: Dan L^{-T}b = mn. dws L^{T}mn = b, gebruik dtrtrs om mn te vind.
prlite::ColVector<double> SqrtMVG::getMean() const {

  PRLITE_ASSERT(!vacuous, "Can not get the mean of a vacuous density");

  char uplo = 'L';
  char trans = 'N';
  char diag = 'N';
  int orderLHS = vars.size();
  int lda = vars.size();
  int incx = 1;

  // Step 1: first solve b in L * b = h
  prlite::ColVector<double> mn; mn = h;
  dtrsv_(&uplo, &trans, &diag, &orderLHS, &L(0,0), &lda, &mn[0], &incx);
  // should be fine with posdef L

  // Step 2: now solve mn in L^{T}mn = b
  trans = 'T';
  dtrsv_(&uplo, &trans, &diag, &orderLHS, &L(0,0), &lda, &mn[0], &incx);
  // should be fine with posdef L

  return mn;
} // getMean

// C = K^{-1} = (LL^T)^{-1} = L^{-T}L^{-1}  gebruik dtrtri vir die inverse en dtrmm vir die vermenigvuldiging. Moet redelik kophou omdat ons C as ROWDENSE gebruik.
prlite::RowMatrix<double> SqrtMVG::getCov() const {

  PRLITE_ASSERT(!vacuous, "Can not get the covariance of a vacuous density");

  // TRICK: Lapack sees ROWDENSE LOWER theCov as transposed COLDENSE UPPER
  prlite::RowMatrix<double> theCov;
  theCov = copyLowerToLowerAndZeroUpper(L);

  char uplo = 'U';
  char diag = 'N';
  int nDims = vars.size();
  int lda = nDims;
  int fail;

  dtrtri_(&uplo, &diag,  &nDims, &theCov(0,0), &lda, &fail); // invert triangular
  PRLITE_ASSERT(!fail, "Could not invert L");


  // TRICK: Lapack sees ROWDENSE LOWER theCov as transposed COLDENSE UPPER
  uplo = 'U';

  // We want to calculate ROWDENSE B^T * B. For COLDENSE Lapack this
  // is B*B^T. I.e. transpose the lhs and multiply with B^T on rhs.
  char transA = 'T'; char side = 'R';

  int nRowsB = nDims;
  int nColsB = nDims;
  double alpha = 1.0;
  int ldb = nDims;

  // CvD: START OF TEMP FIX
  prlite::RowMatrix<double> theCov_copy( noOfVars(), noOfVars() );
  theCov_copy = theCov;

  dtrmm_(&side, &uplo, &transA, &diag, &nRowsB, &nColsB, &alpha, &theCov_copy(0,0), &lda, &theCov(0,0), &ldb); // triang * mat
  // CvD: END OF TEMP FIX

  return theCov;
} // getCov

// K = LL^{T}  gebruik dtrmm vir die matriksvermenigvuldiging
prlite::RowMatrix<double> SqrtMVG::getK() const {

  prlite::RowMatrix<double> theK( noOfVars(), noOfVars() );
  if (vacuous) {
    theK.assignToAll(0.0);
    return theK;
  } // if

  theK = copyLowerToUpperAndZeroLower(L);
  // theK now is ROWDENSE UPPER L^T. COLDENSE Lapack sees it as LOWER
  // L. We want ROWDENSE L * L^T. Ask Lapack to calculate L^T * L;

  int nDims = vars.size();
  char side = 'R';
  char uplo = 'L';
  char transA = 'T';
  char diag = 'N';
  int nRowsB = nDims;
  int nColsB = nDims;
  double alpha = 1.0;
  int lda = nDims;
  int ldb = nDims;

  // CvD: START OF TEMP FIX
  prlite::RowMatrix<double> theK_copy( noOfVars(), noOfVars() );
  theK_copy = theK;

  dtrmm_(&side, &uplo, &transA, &diag, &nRowsB, &nColsB, &alpha, &theK_copy(0,0), &lda, &theK(0,0), &ldb); // triang * mat
  // CvD: END OF TEMP FIX
  
  return theK;
} // getK

bool SqrtMVG::estMnR(
  const prlite::ColMatrix<double>& samples,
  prlite::ColVector<double>& mn,
  prlite::ColMatrix<double>& R,
  bool biased) {

  size_t nSamples = samples.cols();
  size_t nDims = samples.rows();

  R.resize(nDims,nDims);
  R.assignToAll(0.0);
  mn.resize(nDims);
  mn.assignToAll(0.0);

  for (size_t n = 0; n < nSamples; n++) {
    mn += samples.col(n);
    cholUR( R, samples.col(n) );
  } // for
  mn /= nSamples;

  if (biased) {
    R /= sqrt(nSamples-1);
    if ( !cholDR( R, mn*( double(nSamples)/(nSamples-1) ) ) ){
      PRLITE_LOGERROR(PRLITE_pLog, "SqrtMVG::estMnR", "upper righthand R Matrix not positive definite");
      return false;
    } // if
  } // if

  else{
    R /= sqrt(nSamples);
    if ( !cholDR(R, mn) ) {
      PRLITE_LOGERROR(PRLITE_pLog, "SqrtMVG::estMnR", "upper righthand R Matrix not positive definite");
      return false;
    } // if
  } // if

  return true;
} // estMnR

SqrtMVG* SqrtMVG::constructFromSamples(
    const emdw::RVIds& theVars,
    const prlite::ColMatrix<double>& samples,
    bool biased,
    bool theAutoNormed,
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
    rcptr<Sampler> theSampler){

  prlite::ColMatrix<double> theR;
  prlite::ColVector<double> theH;
  prlite::ColMatrix<double> theL;
  PRLITE_ASSERT(estMnR(samples, theH, theR, biased), "constructFromSamples : could not estimate params from samples");   // theH now contains the mean

  char uplo = 'U';
  char diag = 'N';
  int nDims = theVars.size();
  int lda = nDims;
  int fail;

  dtrtri_(&uplo, &diag,  &nDims, &theR(0,0), &lda, &fail); // invert triangular
  PRLITE_ASSERT(!fail, "Could not invert R");

  theL = theR.transpose();

  uplo = 'L';
  char trans = 'T';
  int incx = 1;
  // L^T * mn
  dtrmv_(&uplo, &trans, &diag, &nDims, &theL(0,0), &lda, &theH[0], &incx); // triang * vec
  trans = 'N'; // this one uses L
  // L * L^T * mn
  dtrmv_(&uplo, &trans, &diag, &nDims, &theL(0,0), &lda, &theH[0], &incx); // triang * vec

  double theG;
  if (!theAutoNormed) {theG = SqrtMVG::getNormedG(theL, theH);} // if

  return new SqrtMVG(
    theVars, theL, theH, theG, theAutoNormed, presorted,
    theMarginalizer,
    theInplaceNormalizer, theNormalizer,
    theInplaceAbsorber, theAbsorber,
    theInplaceCanceller, theCanceller,
    theObserveAndReducer,
    theInplaceDamper,
    theSampler);
} // constructFromSamples

SqrtMVG* SqrtMVG::constructFromSigmaPoints(
    const emdw::RVIds& xVars,
    const prlite::ColMatrix<double>& xSPts,
    const emdw::RVIds& yVars,
    const prlite::ColMatrix<double>& ySPts,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed,
    bool presortedX,
    bool presortedY,
    rcptr<Marginalizer> theMarginalizer,
    rcptr<InplaceNormalizer> theInplaceNormalizer,
    rcptr<Normalizer> theNormalizer,
    rcptr<InplaceAbsorber> theInplaceAbsorber,
    rcptr<Absorber> theAbsorber,
    rcptr<InplaceCanceller> theInplaceCanceller,
    rcptr<Canceller> theCanceller,
    rcptr<ObserveAndReducer> theObserveAndReducer,
    rcptr<InplaceDamper> theInplaceDamper,
    rcptr<Sampler> theSampler){

  size_t xDim = xVars.size();
  size_t yDim = yVars.size();
  size_t nPts = xSPts.cols();


  PRLITE_ASSERT(ySPts.cols() > int(2*yDim),
         "Rank deficient Y sigma points, have : "
                << ySPts.cols() << ", should at least be: " << 2*yDim);
  PRLITE_ASSERT(int(xDim) == xSPts.rows(),
         "The dimensions of the X vars and X sigma points differs: "
         << xDim << " vs " << xSPts.rows());
  PRLITE_ASSERT(int(yDim) == ySPts.rows(),
         "The dimensions of the Y vars and Y sigma points differs: "
         << yDim << " vs " << ySPts.rows());
  PRLITE_ASSERT(ySPts.cols() == int(nPts),
         "The number of sigma points for X and Y differs: "
         << nPts << " vs " << ySPts.cols());
  PRLITE_ASSERT(noiseR.rows() == int(yDim),
         "The dimension of Y and noisediffers: "
         << yDim << " vs " << noiseR.rows());


  prlite::ColVector<double> xMn, yMn;
  prlite::ColMatrix<double> xR, yR;
  PRLITE_ASSERT(estMnR(xSPts, xMn, xR, false), "constructFromSigmaPoints : could not estimate params from X sigma points");
  PRLITE_ASSERT(estMnR(ySPts, yMn, yR, false), "constructFromSigmaPoints : could not estimate params from Y sigma points");

  // calc Cov_yx
  prlite::ColMatrix<double> yxC(yDim,xDim); yxC.assignToAll(0.0);
  for (size_t n = 0; n < nPts; n++) {
    yxC += ySPts.col(n) * xSPts.col(n).transpose();
  } // for
  yxC /= nPts;
  yxC -= yMn*xMn.transpose();

  // find R11
  prlite::ColMatrix<double> R11; R11 = noiseR;
  for (size_t n = 0; n < yDim; n++) {
    cholUR( R11, yR.col(n) );
  } // for n
  //cout << "R11: " << R11 << endl;

  // solve for R01T. We'll reuse yxC as an alias for R01T
  char uplo = 'U'; // R11 is upper
  char trans = 'N';
  char diag = 'N';
  int orderLHS = yDim;
  int nColsRHS = xDim;
  int nRowsLHS = orderLHS;
  int nRowsRHS = nRowsLHS;
  int fail;

  dtrtrs_( &uplo, &trans, &diag, &orderLHS, &nColsRHS,
           &R11(0,0), &nRowsLHS,
           &yxC(0,0), &nRowsRHS, &fail); // solves triangular system
  PRLITE_ASSERT(!fail, "Could not solve Cov_YX from R11");
  //cout << "R01T: " << yxC << endl;

  // now find  R00. We'll reuse xR as an alias for R00
  prlite::ColVector<double> rown(xDim); // need to get it contiguous in any case
  for (size_t n = 0; n < yDim; n++) {
    rown = yxC.row(n).transpose();
    PRLITE_ASSERT(cholDR(xR, rown), "constructFromSigmaPoints : could not do the downdate to find R00");
  } // for n
  //cout << "R00: " << xR << endl;

  //OK, we can now collect them in R
  int newDim = xDim + yDim;

  // we will initialize L with R^T
  prlite::ColMatrix<double> theL(newDim, newDim); theL.assignToAll(0.0);

  // copy R00T
  //cout << "R00: " << xR << endl;
  for (size_t c = 0; c < xDim; c++){
    for (size_t r = 0; r <= c; r++){
      theL(c,r) = xR(r,c);
    } // for
  } // for
  //cout << "L: " << theL << endl;

  // copy R01T
  //cout << "R01T: " << yxC << endl;
  for (size_t r = 0; r < yDim; r++){
    for (size_t c = 0; c < xDim; c++){
      theL(xDim+r,c) = yxC(r,c);
    } // for
  } // for
  //cout << "L: " << theL << endl;

  // copy R11T
  //cout << "R11: " << R11 << endl;
  for (size_t c = 0; c < yDim; c++){
    for (size_t r = 0; r <= c; r++){
      theL(xDim+c,xDim+r) = R11(r,c);
    } // for
  } // for
  //cout << "L: " << theL << endl;

  uplo = 'L';
  int lda = newDim;

  dtrtri_(&uplo, &diag,  &newDim, &theL(0,0), &lda, &fail); // invert triangular
  PRLITE_ASSERT(!fail, "Could not invert R^T");
  //cout << "L: " << theL << endl;

  // calculate the new value for h

  // first collect the new mean
  prlite::ColVector<double> theH(newDim);
  for (size_t r = 0; r < xDim; r++){
    theH[r] = xMn[r];
  } // for
  for (size_t r = 0; r < yDim; r++){
    theH[xDim+r] = yMn[r];
  } // for

  // and calculate h from that
  uplo = 'L';
  trans = 'T';
  int incx = 1;
  // L^T * mn
  dtrmv_(&uplo, &trans, &diag, &newDim, &theL(0,0), &lda, &theH[0], &incx); // triang * vec
  trans = 'N'; // this one uses L
  // L * L^T * mn
  dtrmv_(&uplo, &trans, &diag, &newDim, &theL(0,0), &lda, &theH[0], &incx); // triang * vec

  double theG;
  if (!theAutoNormed) {theG = SqrtMVG::getNormedG(theL, theH);} // if

  // Let the constructor figure out unsorted variable orders if required.
  emdw::RVIds theVars( xVars.begin(), xVars.end() );
  theVars.insert( theVars.end(), yVars.begin(), yVars.end() );
  bool presorted = presortedX and presortedY and (yVars[0] > xVars[xVars.size()-1]);

  return new SqrtMVG(
    theVars, theL, theH, theG, theAutoNormed, presorted,
    theMarginalizer,
    theInplaceNormalizer, theNormalizer,
    theInplaceAbsorber, theAbsorber,
    theInplaceCanceller, theCanceller,
    theObserveAndReducer,
    theInplaceDamper,
    theSampler);

} //constructFromSigmaPoints

// Defunct version based on constructFromSigmaPoints. Has constraint
// that yDim <= xDim.
SqrtMVG* SqrtMVG::_constructAffineGaussian(
    const SqrtMVG& xPdf,
    const prlite::ColMatrix<double>& aMat,
    const prlite::ColVector<double>& cVec,
    const emdw::RVIds& newVars,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed,
    bool presortedY,
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

  unsigned xDim = xPdf.noOfVars();
  unsigned yDim = newVars.size();
  PRLITE_ASSERT(xDim == unsigned(aMat.cols()), "The number of variables in x does not match the number of columns in aMat");
  PRLITE_ASSERT(yDim == unsigned(aMat.rows()), "The number of variables in y does not match the number of rows in aMat");
  PRLITE_ASSERT(yDim <= xDim,
                "Because we are using the constructFromSigmaPoints function, the y-dim (" << yDim << ") may not exceed the x-dim (" << xDim << ")");
  prlite::ColMatrix<double> spX; spX = xPdf.getSigmaPoints();
  prlite::ColMatrix<double> spY = aMat * spX;
  if ( cVec.size() ) { // add offset if required
    for (int c = 0; c < spY.cols(); c++) {
      spY.col(c) += cVec;
    } // for
  } // if

  return constructFromSigmaPoints(xPdf.getVars(), spX, newVars, spY, noiseR, theAutoNormed, true, presortedY,
                                  theMarginalizer,
                                  theInplaceNormalizer, theNormalizer,
                                  theInplaceAbsorber, theAbsorber,
                                  theInplaceCanceller, theCanceller,
                                  theObserveAndReducer,
                                  theInplaceDamper,
                                  theSampler);

} // _constructAffineGaussian

SqrtMVG* SqrtMVG::constructAffineGaussian(
    const SqrtMVG& xPdf,
    const prlite::ColMatrix<double>& aMat,
    const prlite::ColVector<double>& cVec,
    const emdw::RVIds& yVars,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed,
    bool presortedY,
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

  unsigned xDim = xPdf.noOfVars();
  unsigned yDim = yVars.size();
  PRLITE_ASSERT(xDim == unsigned(aMat.cols()), "The number of variables in x does not match the number of columns in aMat");
  PRLITE_ASSERT(yDim == unsigned(aMat.rows()), "The number of variables in y does not match the number of rows in aMat");

  // -------------------------------------------------------------------------------------
  // --- First we need to find upper triangular R such that R*R' = cov(x). result in xRcv
  // -------------------------------------------------------------------------------------

  // First get L' via the copyLowerToUpper. Because K = LL' ==> xCov = L^(-T)L^(-1). That means
  // inv(L^(-T)) is the COLDENSE UPPER R we require.
  char uplo = 'U';
  char diag = 'N';
  int nDims = xDim;
  int fail;

  prlite::ColMatrix<double> xRcv = copyLowerToUpperAndZeroLower(xPdf.L);

  dtrtri_(&uplo, &diag,  &nDims, &xRcv(0,0), &nDims, &fail); // inv triangular
  PRLITE_ASSERT(!fail, "Could not invert L");
  //cout << "xRcv: " << xRcv << endl;
  // it checks out NOTE: xRT.transpose is what it should be to us -
  // because lapack sees the transpose of a rowdense matrix.

  // cout << "xRcv is: " << xRcv << endl;
  // cout << "x Cov should be: " << xPdf.getCov() << endl;
  // cout << "x Cov via xRcv: " << xRcv*xRcv.transpose() << endl;
  // checks out

  // --------------------------------------------
  // --- Now we need to compute AR, save in theAR
  // --------------------------------------------

  //cout << "aMat: " << aMat << endl;
  prlite::ColMatrix<double> theAR; theAR = aMat;

  char side = 'R';   // xRcv multiplies from the right
  uplo = 'U';        // xRcv is upper triangular
  char transA = 'N'; // xRcv must not be transposed - we want AR
  diag = 'N';        // xRcv not unit triangular
  int nRowsB = theAR.rows();
  int nColsB = theAR.cols();
  double alpha = 1.0;
  int lda = xRcv.rows();
  int ldb = nRowsB;
  dtrmm_(&side, &uplo, &transA, &diag, &nRowsB, &nColsB, &alpha, &xRcv(0,0), &lda, &theAR(0,0), &ldb); // mpy triangular
  // cout << "AR should be " << aMat* xRcv << endl;
  // cout << "AR is: " << theAR << endl;
  // it checks out

  // -----------------------------
  // --- From this we can find Q11
  // -----------------------------

  // now lets find Q11
  prlite::ColMatrix<double> Q11; Q11 = noiseR;
  for (int c = 0; c < theAR.cols(); c++) { // we are using the rows, because it is the transposed
    cholUR( Q11, theAR.col(c) );
  } // for
  // cout << "Q11 is: " << Q11 << endl;
  // cout << "Q11*Q11' should be: " << noiseR*noiseR.transpose() + aMat*xRcv * (aMat*xRcv).transpose() << endl;
  // cout << "Q11*Q11' is: " << Q11*Q11.transpose() << endl;
  // checks out

  // ------------------------------------
  // --- Now we need ARR', save in theAR
  // ------------------------------------

  uplo = 'U';        // xRcv is upper triangular
  diag = 'N';        // xRcv not unit triangular
  transA = 'T';      // xRcv must now be transposed - we want ARR'
  side = 'R';        // xRcv' multiplies from the right
  alpha = 1.0;
  nRowsB = theAR.rows();
  nColsB = theAR.cols();
  lda = xRcv.rows();
  ldb = nRowsB;
  dtrmm_(&side, &uplo, &transA, &diag, &nRowsB, &nColsB, &alpha, &xRcv(0,0), &lda, &theAR(0,0), &ldb); // mpy triangular
  // cout << "ARR' should be: " << aMat*xRcv*xRcv.transpose() << endl;
  // cout << "ARR' is: " << theAR << endl;
  // checks out

  // -------------------------------------------------
  // --- Now we need to solve for Q01T, save in theAR
  // -------------------------------------------------

  // solve for Q01'. We'll reuse theAR (now containing ARR') as an alias for Q01'
  uplo = 'U';                   // Q11 is upper
  char trans = 'N';             // Q11 is not transposed
  diag = 'N';                   // Q11 is not unit triangular
  int orderLHS = yDim;
  int nColsRHS = xDim;
  int nRowsLHS = orderLHS;
  int nRowsRHS = nRowsLHS;

  dtrtrs_( &uplo, &trans, &diag, &orderLHS, &nColsRHS,
           &Q11(0,0), &nRowsLHS,
           &theAR(0,0), &nRowsRHS, &fail); // solves triangular system
  PRLITE_ASSERT(!fail, "Could not solve Cov_YX from R11");
  // cout << "Q01': " << theAR << endl;
  // cout << "Q11*Q01': " << Q11 * theAR << endl;
  // checks out

  // ---------------------------------------------------
  // --- Now we need to downdate to find Q00, save in xRcv
  // ---------------------------------------------------

  prlite::ColMatrix<double> savR; savR = xRcv; // just for checking

  // now find  Q00. We'll reuse xRcv as an alias for Q00
  prlite::ColVector<double> rown(xDim); // need to get it contiguous in any case
  for (size_t n = 0; n < yDim; n++) {
    rown = theAR.row(n).transpose();
    PRLITE_ASSERT(cholDR(xRcv, rown), "constructFromSigmaPoints : could not do the downdate to find R00");
  } // for n
  // cout << "Q00: " << xRcv << endl;
  // cout << "Q00*Q00T should be: " << savR*savR.transpose() - theAR.transpose()*theAR << endl;
  // cout << "Q00*Q00T is: " << xRcv*xRcv.transpose() << endl;
  // checks out

  // ----------------------------------------------------------------------
  // --- repack the Q's into a single transposed matrix that we can invert
  // ----------------------------------------------------------------------

  //OK, we can now collect them in R
  int newDim = xDim + yDim;

  //We now have Cov = RR' ==> K = R^(-T)R^(-1). That means we will
  //initialize L with R^T so that the inverse ends up LOWER COLDENSE
  prlite::ColMatrix<double> theL(newDim, newDim); theL.assignToAll(0.0);

  // copy Q00 - this lives in xRcv

  //cout << "Q00: " << xRcv << endl;
  for (size_t c = 0; c < xDim; c++){
    for (size_t r = 0; r <= c; r++){
      theL(c,r) = xRcv(r,c);
    } // for
  } // for
  //cout << "L: " << theL << endl;

  // copy Q01T
  for (size_t r = 0; r < yDim; r++){
    for (size_t c = 0; c < xDim; c++){
      theL(xDim+r,c) = theAR(r,c);
    } // for
  } // for
  //cout << "L: " << theL << endl;

  // copy Q11T
  //cout << "Q11: " << Q11 << endl;
  for (size_t c = 0; c < yDim; c++){
    for (size_t r = 0; r <= c; r++){
      theL(xDim+c,xDim+r) = Q11(r,c);
    } // for
  } // for
  // cout << "L: " << theL << endl;
  // cout << "Cov is: " << theL.transpose() * theL << endl;
  // checks out

  // ----------------------------------
  // --- and invert to find L, in theL
  // ----------------------------------

  uplo = 'L';
  lda = newDim;

  dtrtri_(&uplo, &diag,  &newDim, &theL(0,0), &lda, &fail); // inv triangular
  PRLITE_ASSERT(!fail, "Could not invert R^T");
  //cout << "L: " << theL << endl;
  // checks out

  // ------------------------------------
  // --- calculate the new value for h
  // ------------------------------------

  // first collect the new mean
  prlite::ColVector<double> theH(newDim);
  prlite::ColVector<double> xMn = xPdf.getMean();
  prlite::ColVector<double> yMn = aMat * xMn;
  if ( cVec.size() ) {
    yMn += cVec;
  } // if

  for (size_t r = 0; r < xDim; r++){
    theH[r] = xMn[r];
  } // for
  for (size_t r = 0; r < yDim; r++){
    theH[xDim+r] = yMn[r];
  } // for

  // and calculate h from that
  uplo = 'L';
  trans = 'T';
  int incx = 1;
  // L^T * mn
  dtrmv_(&uplo, &trans, &diag, &newDim, &theL(0,0), &lda, &theH[0], &incx); // triang * vec
  trans = 'N'; // this one uses L
  // L * L^T * mn
  dtrmv_(&uplo, &trans, &diag, &newDim, &theL(0,0), &lda, &theH[0], &incx); // triang * vec

  double theG;
  if (!theAutoNormed) {theG = SqrtMVG::getNormedG(theL, theH);} // if

  // cout << "theH: " << theH << endl;
  // cout << "theG: " << theG << endl;
  // checks out

  // Let the constructor figure out unsorted variable orders if required.
  emdw::RVIds theVars( xPdf.vars.begin(), xPdf.vars.end() );
  theVars.insert( theVars.end(), yVars.begin(), yVars.end() );
  bool presorted = presortedY and (yVars[0] > xPdf.vars[xDim-1]);

  return new SqrtMVG(
    theVars, theL, theH, theG, theAutoNormed, presorted,
    theMarginalizer,
    theInplaceNormalizer, theNormalizer,
    theInplaceAbsorber, theAbsorber,
    theInplaceCanceller, theCanceller,
    theObserveAndReducer,
    theInplaceDamper,
    theSampler);
} // constructAffineGaussian

//======================= Non-virtual Member Functions ========================

// g = -0.5 * (mn.transpose() * h) - nDims*HALFLOG2PI + log(sqrtDetK);
//
// Of in skoon kanoniese vorm:
//
// mn' * h = (L^{-1} * h)' * (L^{-1} * h)
// Laat nou (L^{-1} *h = x. Dan h = L*x. Kan x oplos met 'n dtrtrs
double SqrtMVG::getNormedG(const prlite::ColMatrix<double>& theL, const prlite::ColVector<double>& theH){

  char uplo = 'L';
  char trans = 'N';
  char diag = 'N';
  int orderLHS = theH.size();
  int lda = orderLHS;
  int incx = 1;

  // first solve x in L * x = h
  prlite::ColVector<double> x; x = theH;
  double* l00 = const_cast<double*> ( theL.getStartAddr() );
  dtrsv_(&uplo, &trans, &diag, &orderLHS, l00, &lda, &x[0], &incx);

  double sqrtDetK = theL(0,0);
  for (int k = 1; k < orderLHS; k++) {
    sqrtDetK *= theL(k,k);
  } // for

  return ( -0.5 * (x.transpose() * x) - orderLHS*HALFLOG2PI + log(sqrtDetK) );
} // getNormedG

double SqrtMVG::getMass() const {
  if (vacuous) {
    PRLITE_ASSERT(false, "Can not get mass of vacuous pdf");
    return std::numeric_limits<double>::infinity();
  } // if
  if (autoNormed) {return 1.0;} // if
  return exp( g-getNormedG(L,h) );
} // getMass

// CAREFUL, this does NOT set the total volume to newMass, it
// multiplies the existing volume with newMass.
void SqrtMVG::adjustMass(double newMass){

  PRLITE_ASSERT(!vacuous, "It does not make sense to adjust the mass of a vacuous density");
  // or should we simply do nothing thereby leaving it at infinity?

  if (autoNormed) {
    g = getNormedG(L,h);
    autoNormed = false;
  } // if
  g += log(newMass);
} // adjustMass

void SqrtMVG::makeVacuous(){
  vacuous = true;
  // unneccesary, but just to make sure
  L.resize(0,0);
  h.resize(0);
} // makeVacuous

double SqrtMVG::mahalanobis(const prlite::ColVector<double>& observedVals) const{

  if (vacuous) {return std::numeric_limits<double>::infinity();} // if

  prlite::ColVector<double> delta = observedVals-getMean();
  char uplo = 'L';
  char trans = 'T';  // we want L'
  char diag = 'N';
  int nDims = noOfVars();
  int lda = nDims;
  int incx = 1;
  // L' * delta
  dtrmv_(&uplo, &trans, &diag, &nDims, &L(0,0), &lda, &delta[0], &incx); // triang * vec
  // delta'*delta
  return ddot_(&nDims, &delta[0], &incx, &delta[0], &incx);
} // mahalanobis


/**
 * Directly calculates the symmetric KL-divergence between two MVG's as  0.5*(KL(P||Q) + KL(Q||P))
 * dist = (-2*dim + deltaU^T (K_P + K_Q) deltaU + trace(K_Q * (K_P)^-1) + trace(K_P* (K_Q)^-1)) *0.25
*/
double SqrtMVG::symmKLdivergence(const Factor *lhsPtr, const Factor *rhsPtr)
{
  using t_colVec = prlite::ColVector<double>;
  using t_colMat = prlite::ColMatrix<double>;

  auto sqrtLhsPtr = dynamic_cast< const SqrtMVG* >(lhsPtr);
  auto sqrtRhsPtr = dynamic_cast< const SqrtMVG* >(rhsPtr);

  PRLITE_ASSERT(sqrtLhsPtr->vars.size() == sqrtRhsPtr->vars.size(),"Gaussians must have the same number of dimentions");
  if (sqrtLhsPtr->isVacuous() || sqrtRhsPtr->isVacuous()) { return std::numeric_limits<double>::infinity(); }

  t_colVec dU = sqrtRhsPtr->getMean() - sqrtLhsPtr->getMean();
  t_colMat Lsum;
  Lsum = sqrtLhsPtr->getL();          //deep copy
  t_colMat rhsL = sqrtRhsPtr->getL(); //shallow copy
  uint dim = dU.size();

  /* K_P + K_Q */
  for (uint i = 0; i < dim; i++){
    cholUL(Lsum, rhsL.col(i));
  }

  /*code from mahalanobis in sqrtmvg for deltaU^T (K) deltaU*/
  char uplo = 'L';
  char trans = 'T'; // we want L'
  char diag = 'N';
  int nDims = double(dim);
  int lda = nDims;
  int incx = 1;
  // L' * delta
  dtrmv_(&uplo, &trans, &diag, &nDims, &Lsum(0, 0), &lda, &dU[0], &incx); // triang * vec
  // delta'*delta
  double dist = ddot_(&nDims, &dU[0], &incx, &dU[0], &incx);

  /*trace(K_Q* (K_P)^-1) + trace(K_P* (K_Q)^-1)*/
  auto K = sqrtLhsPtr->getK();
  auto C = sqrtRhsPtr->getCov().transpose();

  double trace = 0;
  for (uint i{0}; i < dim; i++){
    trace += K.row(i) * C.col(i);
  }

  K = sqrtRhsPtr->getK();
  C = sqrtLhsPtr->getCov().transpose();

  for (uint i{0}; i < dim; i++){
    trace += K.row(i) * C.col(i);
  }

  dist = 0.25 * (dist + trace - 2.0 * double(dim));
  return dist;
}//symmKLdivergence

void
SqrtMVG::export2DMesh(
  const std::string& filename,
  RVIdType xId, RVIdType yId,
  unsigned nPoints) const {

  PRLITE_ASSERT(!vacuous, "Can not show a 2D mesh on a vacuous pdf");

  rcptr<Factor> fPtr = marginalize({xId});
  rcptr<SqrtMVG> gPtr = std::dynamic_pointer_cast<SqrtMVG>(fPtr);
  double mnX = gPtr->getMean()[0];
  double stX = sqrt(gPtr->getCov()(0,0));
  double minX = mnX - 4.0*stX;
  double maxX = mnX + 4.0*stX;
  double delX = (maxX-minX)/(nPoints-1);
  prlite::ColVector<double> x(nPoints); x.set(minX, delX);

  fPtr = marginalize({yId});
  gPtr = std::dynamic_pointer_cast<SqrtMVG>(fPtr);
  double mnY = gPtr->getMean()[0];
  double stY = sqrt(gPtr->getCov()(0,0));
  double minY = mnY - 4.0*stY;
  double maxY = mnY + 4.0*stY;
  double delY = (maxY-minY)/(nPoints-1);
  prlite::ColVector<double> y(nPoints); y.set(minY, delY);

  if (xId != yId) fPtr = marginalize({xId,yId});
  else fPtr = marginalize({xId});
  gPtr = std::dynamic_pointer_cast<SqrtMVG>(fPtr);

  prlite::RowMatrix<double> fvals(nPoints, nPoints);
  for (int r = 0; r < int(nPoints); r++) { //y vals
    for (int c = 0; c < int(nPoints); c++) { // xvals
      if (xId != yId) fvals(r,c) = gPtr->potentialAt({xId,yId}, {x[c],y[r]});
      else  fvals(r,c) = gPtr->potentialAt({xId}, {x[c]});
    } // for c
  } // for r

  ofstream ofile( filename.c_str() );
  if (!ofile) {
    PRLITE_FAIL("Could not open file '" << filename << "'");
  } // if
  for (int c = 0; c < int(nPoints); c++) { ofile << x[c] << " "; } ofile << endl;
  for (int r = 0; r < int(nPoints); r++) { ofile << y[r] << " "; } ofile << endl;
  for (int r = 0; r < int(nPoints); r++) { //y vals
    for (int c = 0; c < int(nPoints); c++) { // xvals
      ofile << fvals(r,c) << " ";
    } // for c
    ofile << endl;
  } // for r

  ofile.close();

} // export2DMesh

//================================== Friends ==================================

/// output
std::ostream& operator<<( std::ostream& file,
                          const SqrtMVG& n ){
  return n.txtWrite(file);
} // operator<<

/// input
std::istream& operator>>( std::istream& file,
                          SqrtMVG& n ){
  return n.txtRead(file);
} // operator>>


//======================= Factor operators ====================================


// ***************************************************************************
// *************************** Family 0 **************************************
// *************************** Sample ****************************************
// ***************************************************************************

// Identify object
const std::string& SqrtMVG_Sampler::isA() const {
  static const std::string CLASSNAME("SqrtMVG_Sampler");
  return CLASSNAME;
} // isA

void SqrtMVG_Sampler::sampleCentered(
  const SqrtMVG* lhsPtr,
  prlite::ColMatrix<double>& x,
  size_t nSamples) {

  PRLITE_ASSERT(!lhsPtr->vacuous, "Can not sample from a vacuous pdf");

  static normal_distribution<> stdGauss(0.0, 1.0); // mn, stddev

  // TODO: upgrade to iterators
  size_t dim = lhsPtr->noOfVars();
  if ( !x.rows() ) {x.resize(dim, nSamples);} // if
  for (size_t c = 0; c < nSamples; c++) {
    for (size_t r = 0; r < dim; r++) {
      x(r,c) = stdGauss( emdw::randomEngine() );
    } // for r
  } // for c

  // K = L* L^T
  // C = L^{-T} * L^{-1}
  // y = L^{-T} * x
  // L^T * y = x

  char uplo = 'L';
  char trans = 'T';
  char diag = 'N';
  int orderLHS = dim;
  int nColsRHS = nSamples;
  int nRowsLHS = orderLHS;
  int nRowsRHS = nRowsLHS;
  int fail;

  dtrtrs_( &uplo, &trans, &diag, &orderLHS, &nColsRHS,
           &lhsPtr->L(0,0), &nRowsLHS,
           &x(0,0), &nRowsRHS, &fail); // solves triangular system
  PRLITE_ASSERT(!fail, "Could not solve zero mean sample from L");
} // sampleCentered

prlite::ColMatrix<double> SqrtMVG_Sampler::sample(
  const SqrtMVG* lhsPtr,
  size_t nSamples) {

  prlite::ColVector<double> theMn = lhsPtr->getMean();
  size_t dim = theMn.size();
  prlite::ColMatrix<double> x(dim,nSamples);
  sampleCentered(lhsPtr, x, nSamples);
  emdw::RVVals ret(dim);
  for (size_t idx = 0; idx < nSamples; idx++) {
    x.col(idx) += theMn;
  } // for
  return x;
} // sample

emdw::RVVals SqrtMVG_Sampler::sample(const SqrtMVG* lhsPtr) {

  prlite::ColVector<double> theMn = lhsPtr->getMean();
  size_t dim = theMn.size();
  prlite::ColMatrix<double> x(dim,1);
  sampleCentered(lhsPtr, x, 1);
  emdw::RVVals ret(dim);
  for (size_t idx = 0; idx < dim; idx++) {
    ret[idx] = x(idx,0)+theMn[idx];
  } // for

  return ret;
} // sample


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

// Identify object
const std::string& SqrtMVG_InplaceNormalize::isA() const {
  static const std::string CLASSNAME("SqrtMVG_InplaceNormalize");
  return CLASSNAME;
} // isA

void SqrtMVG_InplaceNormalize::inplaceProcess(SqrtMVG* lhsPtr) {

  SqrtMVG& lhs(*lhsPtr);

  if (lhs.vacuous) {return;}
  if (lhs.autoNormed) {return;}
  lhs.g = lhs.getNormedG(lhs.L,lhs.h);

} // inplaceProcess

// Identify object

const std::string& SqrtMVG_Normalize::isA() const {
  static const std::string CLASSNAME("SqrtMVG_Normalize");
  return CLASSNAME;
} // isA

Factor* SqrtMVG_Normalize::process(const SqrtMVG* lhsPtr) {

  // make a copy of the factor
  SqrtMVG* fPtr = new SqrtMVG(*lhsPtr);

  // inplace normalize it
  SqrtMVG_InplaceNormalize ipNorm;
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
const std::string& SqrtMVG_InplaceAbsorb::isA() const {
  static const std::string CLASSNAME("SqrtMVG_InplaceAbsorb");
  return CLASSNAME;
} // isA

void SqrtMVG_InplaceAbsorb::weightedInplaceAbsorb(SqrtMVG* lhsPtr,
                                                  const Factor* rhsFPtr,
                                                  const ValidValue<double>& weight) {

  double sqrtWeight = weight.valid ? sqrt(weight.value) : 1.0;

  // rhsFPtr must be an appropriate SqrtMVG*
  const SqrtMVG* rhsPtr = dynamic_cast<const SqrtMVG*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " SqrtMVG_InplaceAbsorb::inplaceProcess : rhs not a valid SqrtMVG : " << typeid(*rhsFPtr).name() );
  const SqrtMVG& rhs(*rhsPtr);
  SqrtMVG& lhs(*lhsPtr);

  // varU - the union of vars and rhs.vars
  // u2l  - the corresponding indices in union that came from lhs.vars
  // u2r  - the corresponding indices in union that came from rhs.vars
  std::vector<int> u2l, u2r;
  emdw::RVIds varU = sortedUnion(lhs.vars, rhs.vars, u2l, u2r);
  size_t newDim = varU.size();
  // cout << "\nvarL: " << lhs.vars
  //      << "\nvarR: " << rhs.vars
  //      << "\nu2l: " << u2l
  //      << "\nu2r: " << u2r;
  //      << "\nl2u: " << l2u
  //      << "\nr2u: " << r2u << endl;

  // CASE 0: both are vacuous, return a vacuous over the union
  // scope. It is important to do this here since we do not cater for
  // both sides vacuous later on.
  if (lhs.vacuous and rhs.vacuous){
    lhs.classSpecificConfigure(
      varU,
      lhs.autoNormed,
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
    return;
  } // if

  bool subsetR = (u2l.size() == newDim);  // rhs contributed no extra vars
  bool subsetL = (u2r.size() == newDim);  // lhs contributed no extra vars

  // CASE 1: both factors have the same scope. We choose to directly
  // work with the datastructures of the lhs instead of following the
  // more kosher inplace construction route.
  if (subsetR and subsetL) {

    if (rhs.vacuous) { return; }
    if (lhs.vacuous) { // rhs specifically NOT vacuous.

      // matching scope, lhs vacuous, rhs not vacuous

      if (weight.valid) {

        // JADP 16 March 2019
        prlite::ColMatrix<double> newL; newL = rhs.L*sqrtWeight;
        prlite::ColVector<double> newH; newH = rhs.h*weight.value;
        double newG = rhs.g*weight.value;

        lhs.classSpecificConfigure(
          rhs.vars,
          newL,
          newH,
          newG,
          lhs.autoNormed,
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
        return;

      } // if

      else { // weight not valid
        lhs.classSpecificConfigure(
          rhs.vars,
          rhs.L,
          rhs.h,
          rhs.g,
          lhs.autoNormed,
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
        return;
      } // else
    } // if lhs.vacuous

    // neither vacuous
    if (weight.valid) {
      lhs.h += rhs.h*weight.value;
      lhs.g += rhs.g*weight.value;
    } // if
    else {
      lhs.h += rhs.h;
      lhs.g += rhs.g;
    } // else

    // add K's using cholUL on the L's
    if (weight.valid) {
      for (size_t c = 0; c < newDim; c++) {
        //cout << __FILE__ << __LINE__ << " " << "col["<< c << "]: " << rhs.L.col(c);
        cholUL(lhs.L, rhs.L.col(c)*sqrtWeight);
      } // for
    } // if
    else {
      for (size_t c = 0; c < newDim; c++) {
        //cout << __FILE__ << __LINE__ << " " << "col["<< c << "]: " << rhs.L.col(c);
        cholUL( lhs.L, rhs.L.col(c) );
      } // for
    } // else

    return;
  } // if

  // CASE 2: rhs is a subset of lhs but not vice versa. This one is
  // particularly important since it happens all the time in an
  // PGM. We therefore choose to directly work with the datastructures
  // of the lhs instead of following the more kosher inplace
  // construction route.
  else if (subsetR) { // rhs is a subset of lhs but not vice versa

    // SPECIAL CASE 2: rhs is vacuous, leave lhs as it is.
    if (rhs.vacuous) {return;} // if

    // SPECIAL CASE 2: lhs is vacuous, but rhs specifically is NOT.
    if (lhs.vacuous) {

      // rhs is a subset of vacuous lhs. Result will be zero-definite!
      //  Are we ok with this? weighted version also NIY
      PRLITE_ASSERT(false, "Currently SqrtMVG does not support non-vacuous but zero-definite pdfs");

      lhs.vacuous = false;
      lhs.L.resize(newDim, newDim); lhs.L.assignToAll(0.0);
      lhs.h.resize(newDim); lhs.h.assignToAll(0.0);

      // Since lhs.K == 0 we can directly use rhs.L
      for (size_t r = 0; r < u2r.size(); r++) { //
        lhs.h[u2r[r]] = rhs.h[r];
        for (size_t c = 0; c <= r; c++) {
          lhs.L(u2r[r],u2r[c]) = rhs.L(r,c);
        } // for c
      } // for r

      return;
    } // if

    // neither vacuous with rhs subset of lhs

    // L: add rhs K to lhs by using cholUL
    //cout << "u2r: " << u2r << endl;
    //cout << "rhs.L: " << rhs.L << endl;
    prlite::ColVector<double> rhsVec(newDim);
    for (size_t c = 0; c < u2r.size(); c++) { // we can skip the null columns
      rhsVec.assignToAll(0.0);
      if (weight.valid) {
        for (size_t r = 0; r < u2r.size(); r++) {
          rhsVec[u2r[r]] = rhs.L(r,c)*sqrtWeight;
        } // for
      } // if
      else {
        for (size_t r = 0; r < u2r.size(); r++) {
          rhsVec[u2r[r]] = rhs.L(r,c);
        } // for
      } // else
      //cout << __FILE__ << __LINE__ << " " << "rhs.L_ext: " << rhsVec << endl;
      cholUL(lhs.L, rhsVec);
    } // for c

    // h and g:
    if (weight.valid) {
      for(unsigned r = 0; r < u2r.size(); r++) {
        lhs.h[u2r[r]] += rhs.h[r]*weight.value;
      } // for r
      lhs.g += rhs.g*weight.value;
    } //if
    else {
      for(unsigned r = 0; r < u2r.size(); r++) {
        lhs.h[u2r[r]] += rhs.h[r];
      } // for r
      lhs.g += rhs.g;
    } // else

    return;
  } // else if

  // with the following cases we need new lhs structures
  prlite::ColMatrix<double> newL;
  prlite::ColVector<double> newH;
  double newG = 0;


  // CASE 3: lhs is a subset of rhs
  if (subsetL) {

    // SPECIAL CASE 3: lhs is vacuous and subset of rhs. It must be updated to be the
    // variable extended version of rhs.
    if (lhs.vacuous) {

      // TODO: we'll do weighted update of  vacuous L a subset of R   later if needed
      PRLITE_ASSERT(!weight.valid, "SqrtMVG_InplaceAbsorb::weightedInplaceAbsorb : Weighted update code NIY");

      lhs.classSpecificConfigure(  // non-vacuous version
        rhs.vars, rhs.L, rhs.h, rhs.g, lhs.autoNormed, true,
        lhs.marginalizer, lhs.inplaceNormalizer, lhs.normalizer,
        lhs.inplaceAbsorber, lhs.absorber, lhs.inplaceCanceller, lhs.canceller,
        lhs.observeAndReducer, lhs.inplaceDamper, lhs.sampler);
      return;
    } // if

    // SPECIAL CASE 3: rhs is vacuous and subset lhs specifically is NOT. lhs
    // must simply be extended to include the extra variables of
    // rhs. This will be zero-definite!  weighted version also NIY
    PRLITE_ASSERT(!rhs.vacuous,
           "Currently SqrtMVG does not support non-vacuous but zero-definite pdfs");

    // neither vacuous and lhs is a subset of rhs

    if (weight.valid) { // W is OK
      // non-vacuous lhs subset of non-vacuous rhs

      // first we need to get the extended lhs L and h with interleaving zeros to match scope of rhs
      newL.resize( u2r.size(), u2r.size() ); newL.assignToAll(0.0);
      newH.resize( u2r.size() ); newH.assignToAll(0.0);
      for (size_t r = 0; r < u2l.size(); r++) {
        newH[ u2l[r] ] = lhs.h[r];
        for (size_t c = 0; c < u2l.size(); c++) {
          newL(u2l[r], u2l[c]) = lhs.L(r,c);
        } // for c
      } // for r

      // now we can update them from weighted rhs //!!!
      newH += rhs.h*weight.value;
      newG = lhs.g + rhs.g*weight.value;
      for (size_t c = 0; c < u2r.size(); c++) {
        //cout << __FILE__ << __LINE__ << " " << "col["<< c << "]: " << rhs.L.col(c);
        cholUL(newL, rhs.L.col(c)*sqrtWeight);
      } // for

    } // if

    else {

      // we work from rhs to get scope right
      newL = rhs.L;
      newH = rhs.h;
      newG = rhs.g;

      // L: add rhs K to lhs by using cholUL
      //cout << "u2l: " << u2l << endl;
      //cout << "lhs.L: " << lhs.L << endl;
      prlite::ColVector<double> lhsVec(newDim);
      for (size_t c = 0; c < u2l.size(); c++) { // we can skip the null columns
        lhsVec.assignToAll(0.0);
        for (size_t r = 0; r < u2l.size(); r++) {
          lhsVec[u2l[r]] = lhs.L(r,c);
        } // for
        //cout << __FILE__ << __LINE__ << " " << "lhs.L_ext: " << lhsVec << endl;
        cholUL(newL, lhsVec);
      } // for c

      // h:
      for(unsigned r = 0; r < u2l.size(); r++) {
        newH[u2l[r]] += lhs.h[r];
      } // for r

      // g:
      newG += lhs.g;
    } // else

  } // if

  // CASE 4: General case
  else {

    PRLITE_ASSERT(!lhs.vacuous and !rhs.vacuous,
           "Currently SqrtMVG does not support non-vacuous but zero-definite pdfs");

    newL.resize(newDim,newDim); newL.assignToAll(0.0);
    newH.resize(newDim); newH.assignToAll(0.0);
    newG = 0.0;

    // first stick an extended lhs L into newL and h into newH
    if (!lhs.vacuous) {
      for (size_t c = 0; c < u2l.size(); c++) { // we can skip the null columns
        newH[u2l[c]] = lhs.h[c];
        for (size_t r = 0; r < u2l.size(); r++) {
          newL(u2l[r],u2l[c]) = lhs.L(r,c);
        } // for r
      } // for c
      newG = lhs.g;
    } // if

    if (weight.valid) {
      // L: add rhs K to lhs by using cholUL and stick rhs h into newH
      if (!rhs.vacuous) {
        prlite::ColVector<double> rhsVec(newDim);
        for (size_t c = 0; c < u2r.size(); c++) { // we can skip the null columns
          newH[u2r[c]] += rhs.h[c]*weight.value;
          rhsVec.assignToAll(0.0);
          for (size_t r = 0; r < u2r.size(); r++) {
            rhsVec[u2r[r]] = rhs.L(r,c)*sqrtWeight;
          } // for
          cholUL(newL, rhsVec);
        } // for c
        newG += rhs.g;
      } // if
    } // if

    else {
      // L: add rhs K to lhs by using cholUL and stick rhs h into newH
      if (!rhs.vacuous) {
        prlite::ColVector<double> rhsVec(newDim);
        for (size_t c = 0; c < u2r.size(); c++) { // we can skip the null columns
          newH[u2r[c]] += rhs.h[c];
          rhsVec.assignToAll(0.0);
          for (size_t r = 0; r < u2r.size(); r++) {
            rhsVec[u2r[r]] = rhs.L(r,c);
          } // for
          cholUL(newL, rhsVec);
        } // for c
        newG += rhs.g;
      } // if
    } // else

  } // else

  lhs.classSpecificConfigure(
    varU,
    newL,
    newH,
    newG,
    lhs.autoNormed,
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


} // weightedInplaceAbsorb

// Identify object
const std::string& SqrtMVG_Absorb::isA() const {
  static const std::string CLASSNAME("SqrtMVG_Absorb");
  return CLASSNAME;
} // isA

Factor* SqrtMVG_Absorb::process(const SqrtMVG* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  SqrtMVG* fPtr = new SqrtMVG(*lhsPtr);

  // inplace normalize it
  SqrtMVG_InplaceAbsorb ipAbsorb;
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
/********************************** Cancel ***********************************/
/*****************************************************************************/

// Identify object
const std::string& SqrtMVG_InplaceCancel::isA() const {
  static const std::string CLASSNAME("SqrtMVG_InplaceCancel");
  return CLASSNAME;
} // isA

void SqrtMVG_InplaceCancel::weightedInplaceCancel(SqrtMVG* lhsPtr,
                                                  const Factor* rhsFPtr,
                                                  const ValidValue<double>& weight) {

  double sqrtWeight = weight.valid ? sqrt(weight.value) : 1.0;

  // rhsFPtr must be an appropriate SqrtMVG*
  const SqrtMVG* rhsPtr = dynamic_cast<const SqrtMVG*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " SqrtMVG_InplaceAbsorb::inplaceProcess : rhs not a valid SqrtMVG : " << typeid(*rhsFPtr).name() );
  const SqrtMVG& rhs(*rhsPtr);
  SqrtMVG& lhs(*lhsPtr);

  if (lhs.vacuous or rhs.vacuous){
    return;
  } //if

  std::vector<int> l2i; // which indices in lhs are shared
  std::vector<int> r2i; // with corresponding indices in rhs
  emdw::RVIds varI = sortedIntersection(lhs.vars, rhs.vars, l2i, r2i);  //intersection of vars
  PRLITE_ASSERT(rhs.vars.size() == varI.size(), "Factor division only supports full sharing of rhs variables");
  size_t newDim = lhs.vars.size();

  // This happens all the time in an PGM. We therefore choose to
  // directly work with the datastructures of the lhs instead of
  // following the more kosher inplace construction route.

  // L: subtract rhs K from lhs by using cholDL
  prlite::ColVector<double> rhsVec(newDim);
  for (size_t c = 0; c < varI.size(); c++) { // we can skip the null columns

    rhsVec.assignToAll(0.0);
    if (weight.valid) {
      for (size_t r = 0; r < varI.size(); r++) {
        rhsVec[l2i[r]] = rhs.L(r2i[r],c)*sqrtWeight;
      } // for
    } // if
    else {
      for (size_t r = 0; r < varI.size(); r++) {
        rhsVec[l2i[r]] = rhs.L(r2i[r],c);
      } // for
    } // else

    if ( !cholDL(lhs.L, rhsVec) ) {
      lhs.makeVacuous();
      return;
    } // if

  } // for c

    // h and g:
  if (weight.valid) {
    for(unsigned r = 0; r < varI.size(); r++) {
      lhs.h[l2i[r]] -= rhs.h[r2i[r]]*weight.value;
    } // for r
    lhs.g -= rhs.g*weight.value;
  } // if
  else {
    for(unsigned r = 0; r < varI.size(); r++) {
      lhs.h[l2i[r]] -= rhs.h[r2i[r]];
    } // for r
    lhs.g -= rhs.g;
  } // else

} // weightedInplaceCancel

// Identify object
const std::string& SqrtMVG_Cancel::isA() const {
  static const std::string CLASSNAME("SqrtMVG_Cancel");
  return CLASSNAME;
} // isA

Factor* SqrtMVG_Cancel::process(const SqrtMVG* lhsPtr,
                                       const Factor* rhsFPtr) {

  // make a copy of the factor
  SqrtMVG* fPtr = new SqrtMVG(*lhsPtr);

  // inplace normalize it
  SqrtMVG_InplaceCancel ipCancel;
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
const std::string& SqrtMVG_Marginalize::isA() const {
  static const std::string CLASSNAME("SqrtMVG_Marginalize");
  return CLASSNAME;
} // isA

Factor* SqrtMVG_Marginalize::process(const SqrtMVG* lhsPtr,
                                     const emdw::RVIds& toKeep,
                                     bool presorted,
                                     const Factor* peekAheadPtr) {

  const SqrtMVG& lhs(*lhsPtr);
  int origDim = lhs.vars.size();
  PRLITE_ASSERT(toKeep.size(), "SqrtMVG_Marginalize : Can not marginalize everything out");

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

  std::vector<int> l2i; // which indices in lhs pot are shared
  std::vector<int> r2i; // with corresponding indices in rhs pot, not used here
  std::vector<int> lni; // indices in lhs that are not in intersection, i.e. to be marginalized out
  std::vector<int> rni; // should be empty, all variables should be in lhs scope
  emdw::RVIds theVars = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);
  PRLITE_ASSERT(theVars.size(), "SqrtMVG_Marginalize : Can not marginalize everything out");
  PRLITE_ASSERT(!rni.size(), "SqrtMVG_Marginalize : Can not marginalize to variables not in scope");
  int xDim = theVars.size();
  int yDim = lhs.noOfVars()-xDim;

  // theVars are the variables we want to keep, but in sorted order
  // lni indexes the variables we want to integrate out, in sorted order

  if (xDim == origDim) { // keep everything
    return lhs.copy();
  } // if

  if (lhs.vacuous) { // return a vacuous
    return new SqrtMVG(
      theVars,
      lhs.autoNormed,
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
  } // if

  prlite::ColMatrix<double>* theLPtr;
  prlite::ColMatrix<double> tmpL; // need it here to not go out of scope too soon
  if (lni[0] < xDim) { // eish, have to do some rather expensive variable reordering

    tmpL.resize( origDim, origDim );
    theLPtr = &tmpL;

    // repack rows in new order
    size_t idx = 0;
    for (size_t n = 0; n < l2i.size(); n++) {
      tmpL.row(idx++) = lhs.L.row(l2i[n]);
    } // for
    for (size_t n = 0; n < lni.size(); n++) {
      tmpL.row(idx++) = lhs.L.row(lni[n]);
    } // for

    // and get it back to lower triangular shape again
    rotateColsToLowerTriangular(tmpL);

  } // if

  else{ // nice and fast
    theLPtr = &lhs.L;
  } // else

  prlite::ColMatrix<double>& theL(*theLPtr);

  prlite::ColVector<double> hx(xDim);
  prlite::ColVector<double> hy(yDim);
  int ySz = 0; int xSz = 0;
  for (int n = 0; n < origDim; n++)
    if (xSz < xDim and theVars[xSz] == lhs.vars[n]) {
      hx[xSz++] = lhs.h[n];
    } // if
    else{
      hy[ySz++] = lhs.h[n];
    } // else

  //NOTE: the ObserveAndReduce operator below shows how one can avoid
  //using these gLinear slices
  // Lxx = L00 (extract Lxx)
  gLinear::gIndexRange rangeX(0, xDim-1);
  gLinear::gIndexRange rangeY(xDim, lhs.noOfVars()-1);
  prlite::ColMatrix<double> Lxx; Lxx = theL.slice(rangeX, rangeX);
  prlite::ColMatrix<double> L10; L10 = theL.slice(rangeY, rangeX);
  prlite::ColMatrix<double> Lyy; Lyy = theL.slice(rangeY, rangeY); // at this stage still L11

// Kxy^T = L10*L00^T    (calculate Kxy^T directly: yDim x xDim)
  prlite::ColMatrix<double> QT; QT = L10;   // KxyT will ultimately be transformed to QT
  char side = 'R';
  char uplo = 'L';
  char transA = 'T';
  char diag = 'N';
  int nRowsB = yDim;
  int nColsB = xDim;
  double alpha = 1.0;
  int lda = xDim;
  int ldb = nRowsB;
  dtrmm_(&side, &uplo, &transA, &diag, &nRowsB, &nColsB, &alpha, &Lxx(0,0), &lda, &QT(0,0), &ldb);
  //cout << "KxyT = L10*L00^T: " << QT << endl;

  // Kyy = L11*L11^T + L10*L10^T  (calculate Lyy via cholUL)
  for (int n = 0; n < L10.cols(); n++) {
    cholUL(Lyy,  L10.col(n) );
  } // for n
  //cout << "Lyy: " << Lyy << endl;

  // Kxy*Lyy^{-T} = Q
  // Q*Lyy^T = Kxy
  // Lyy*Q^T = Kxy^T              (solve for Q^T)
  //KxyT will now get replaced as QT

  uplo = 'L';
  char trans = 'N';
  diag = 'N';
  int orderLHS = yDim;
  int nColsRHS = xDim;
  int nRowsLHS = yDim;
  int nRowsRHS = yDim;
  int fail;

  dtrtrs_( &uplo, &trans, &diag, &orderLHS, &nColsRHS,
           &Lyy(0,0), &nRowsLHS,
           &QT(0,0), &nRowsRHS, &fail); // solves triangular system
  PRLITE_ASSERT(!fail, "Could not solve Q matrix");
  // !!!! TODO Rather return vacuous???

  //cout << "QT: " << QT << endl;

  // L'*L'^T = Lxx*Lxx^T - Q*Q^T  (calculate L' via cholDL)
  for (int n = 0; n < QT.rows(); n++) {
    prlite::ColVector<double> q; q = QT.row(n).transpose();
    PRLITE_ASSERT(cholDL( Lxx, q ), "Could not do downdate to find Lxx");
    // or should we return a vacuous pdf?
  } // for n
  //cout << "L': " << Lxx << endl;

  // now we still need to get h'

  // first solve for v in Lyy v = hy
  uplo = 'L';
  trans = 'N';
  diag = 'N';
  orderLHS = yDim;
  lda = yDim;
  int incx = 1;
  dtrsv_(&uplo, &trans, &diag, &orderLHS, &Lyy(0,0), &lda, &hy[0], &incx);
  // !!! TODO need to check posdef???
  hx -= QT.transpose()*hy;
  //cout << "h': " << hx << endl;

  // lastly we need g
  double gx;
  if (!lhs.autoNormed) {
    gx = SqrtMVG::getNormedG(Lxx,hx);
  } // if

  // and thats it folks
  return new SqrtMVG(
    theVars,
    Lxx,
    hx,
    gx,
    lhs.autoNormed,
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


/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

// Identify object
const std::string& SqrtMVG_ObserveAndReduce::isA() const {
  static const std::string CLASSNAME("SqrtMVG_ObserveAndReduce");
  return CLASSNAME;
} // isA

Factor* SqrtMVG_ObserveAndReduce::process(const SqrtMVG* lhsPtr,
                                          const emdw::RVIds& theVars,
                                          const emdw::RVVals& anyVals,
                                          bool presorted) {


  // and a convenient reference for the lhs factor
  const SqrtMVG& lhs(*lhsPtr);
  int origDim = lhs.vars.size();
  if (theVars.size() == 1) {presorted = true;}

  // if observed nothing, return as is
  if (!theVars.size()) {
    return lhs.copy();
  } // if

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
  std::vector<int> l2i; // which indices in lhs pot are observed (i.e. y)
  std::vector<int> r2i; // with corresponding indices in rhs pot
  std::vector<int> lni; // non-observed lhs variables ARE the new variable set (i.e. x)
  std::vector<int> rni; // remaining rhs variables (i.e. observed vars not part of this factor and thus ignored)
  emdw::RVIds varI = sortedIntersection(lhs.vars, sortedVars, l2i, r2i, lni, rni);

  int yDim = l2i.size();

  // if observation is empty, return as is
  if (!yDim) {
    return lhs.copy();
  } // if

  int xDim = lni.size();

  // if no more rv's left, return an empty factor
  if (!xDim) {
    return new SqrtMVG(
      {},
      lhs.autoNormed,
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
  } // if

  emdw::RVIds newVars = extract<emdw::RVIdType>(lhs.vars,lni);

  if (lhs.vacuous) { // return a vacuous
    return new SqrtMVG(
      newVars,
      lhs.autoNormed,
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
  } // if

  prlite::ColMatrix<double>* theLPtr;
  prlite::ColMatrix<double> tmpL; // need it here to not go out of scope too soon

  if ( l2i[0] < xDim ) { // eish, have to do some rather expensive variable reordering
    tmpL.resize(origDim,origDim);
    theLPtr = &tmpL;

    // repack rows in new order
    size_t idx = 0;
    for (size_t n = 0; n < lni.size(); n++) {
      tmpL.row(idx++) = lhs.L.row(lni[n]);
    } // for
    for (size_t n = 0; n < l2i.size(); n++) {
      tmpL.row(idx++) = lhs.L.row(l2i[n]);
    } // for

    // and get it back to lower triangular shape again
    rotateColsToLowerTriangular(tmpL);
  } // if

  else{ // nice and fast
    theLPtr = &lhs.L;
  } // else

  prlite::ColMatrix<double>& theL(*theLPtr);

  prlite::ColVector<double> newH(xDim);
  prlite::ColVector<double> hy(yDim);
  int ySz = 0; int xSz = 0;
  for (int n = 0; n < origDim; n++) {
    if (xSz < xDim and newVars[xSz] == lhs.vars[n]) {
      newH[xSz++] = lhs.h[n];
    } // if
    else{
      hy[ySz++] = lhs.h[n];
    } // else
  } // for

  // newL = L00                    (extract newL)
  // gLinear::gIndexRange rangeX(0, xDim-1);
  // gLinear::gIndexRange rangeY(xDim, lhs.noOfVars()-1);
  // prlite::ColMatrix<double> newL; newL = theL.slice(rangeX, rangeX); // the new L
  prlite::ColMatrix<double> newL(xDim,xDim);
  for (int c = 0; c < xDim; c++){
    for (int r = 0; r < xDim; r++){
      newL(r,c) = theL(r,c);
    } // for c
  } // for r

  // prlite::ColMatrix<double> L10; L10 = theL.slice(rangeY, rangeX);
  prlite::ColMatrix<double> L10(yDim,xDim);
  for (int c = 0; c < xDim; c++){
    for (int r = 0; r < yDim; r++){
      L10(r,c) = theL(r+xDim,c);
    } // for c
  } // for r

  // prlite::ColMatrix<double> Lyy; Lyy = theL.slice(rangeY, rangeY); // at this stage still L11
  prlite::ColMatrix<double> Lyy(yDim,yDim); // at this stage still L11
  for (int c = 0; c < yDim; c++){
    for (int r = 0; r < yDim; r++){
      Lyy(r,c) = theL(r+xDim,c+xDim);
    } // for c
  } // for r

  for (int n = 0; n < L10.cols(); n++) {
    cholUL(Lyy,  L10.col(n) );
  } // for n

  // the observed values
  prlite::ColVector<double> yVals( r2i.size() );
  if (!presorted){
    for(unsigned i = 0; i < r2i.size(); i++) {
      yVals[i] = static_cast<double>(anyVals[ r2i[ sorted[i] ] ]);
    } // for
  } // if
  else{
    for(unsigned i = 0; i < r2i.size(); i++) {
      yVals[i] = static_cast<double>(anyVals[ r2i[i] ]);
    } // for
  } // else

  // now calculate the new h

  // Kxy = L00*L10^T

  // First calculate L10^T*y. Should be able to do this more efficiently
  prlite::ColVector<double> v = L10.transpose()*yVals;

  // And then do the triangular product
  char uplo = 'L';
  char trans = 'N';
  char diag = 'N';
  int incx = 1;
  dtrmv_(&uplo, &trans, &diag, &xDim, &newL(0,0), &xDim, &v[0], &incx); // triang * vec

  // and finally subtract from newH
  newH -= v;

  // now calculate the new g
  double newG = lhs.g;
  if (!lhs.autoNormed) {
    int incy = 1;
    newG += ddot_(&yDim, &hy[0], &incx, &yVals[0], &incy);

    trans = 'T'; // we need Lyy^T
    dtrmv_(&uplo, &trans, &diag, &yDim, &Lyy(0,0), &yDim, &yVals[0], &incy); // triang * vec
    newG -= 0.5*ddot_(&yDim, &yVals[0], &incy, &yVals[0], &incy);
  } // if

  // and thats it folks
  return new SqrtMVG(
    newVars,
    newL,
    newH,
    newG,
    lhs.autoNormed,
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


/*****************************************************************************/
/*************************** Family 5 ****************************************/
/*************************** Damping *****************************************/
/*****************************************************************************/

// Identify object
const std::string& SqrtMVG_InplaceNoDamping::isA() const {
  static const std::string CLASSNAME("SqrtMVG_InplaceNoDamping");
  return CLASSNAME;
} // isA

double SqrtMVG_InplaceNoDamping::inplaceProcess(
  SqrtMVG* lhsPtr,
  const Factor* rhsPtr,
  double){
  return lhsPtr->distance(rhsPtr);
} // inplaceProcess


//####################################################################################
//####################################################################################
//################# OLD STUFF THAT NEEDS UPDATING ####################################
//####################################################################################
//####################################################################################


/*****************************************************************************/
/*************************** Family 5 ****************************************/
/*************************** Damping *****************************************/
/*****************************************************************************/

/*
// Identify object
const std::string& SqrtMVG_InplaceWeakDamping::isA() const {
  static const std::string CLASSNAME("SqrtMVG_InplaceWeakDamping");
  return CLASSNAME;
} // isA

double SqrtMVG_InplaceWeakDamping::inplaceProcess(
  SqrtMVG* lhsPtr,
  const Factor* rhsPtr,
  double df){  // remember df is the weight for rhs

  const SqrtMVG* dwnPtr = dynamic_cast<const SqrtMVG*>(rhsPtr);
  ASSERT(dwnPtr, "SqrtMVG::distance : rhs is not a SqrtMVG");

  double dist = lhsPtr->distance(dwnPtr);
  //std::cout << "distance: " << dist << std::endl;

  if (df > 0.0) {
    lhsPtr->updateCovarianceForm();
    dwnPtr->updateCovarianceForm();

    // cov
    Matrix<double>  theCov = dwnPtr->cov*df + lhsPtr->cov*(1.0-df);

    // right
    ColVector<double> theMn = dwnPtr->mn*df;
    addOuterProdSym(theCov, dwnPtr->mn,theMn); // already weighted

    //left
    ColVector<double> vWeighted = lhsPtr->mn*(1.0-df);
    addOuterProdSym(theCov, lhsPtr->mn, vWeighted);
    theMn += vWeighted;

    subtractOuterProdSym(theCov,theMn);
    copyLowerToUpper(theCov);

    // std::cout << __FILE__ << __LINE__ << std::endl;
    // std::cout << "lhs::mean: " << lhsPtr->mn << std::endl
    //           << "lhs::cov:  " << lhsPtr->cov << std::endl
    //           << "rhs::mean: " << dwnPtr->mn << std::endl
    //           << "rhs::cov:  " << dwnPtr->cov << std::endl
    //           << "new::mean: " << theMn << std::endl
    //           << "new::cov:  " << theCov << std::endl;

    lhsPtr->reconfigure(theMn, theCov);
  } // if

  return dist;
} // inplaceProcess
*/
