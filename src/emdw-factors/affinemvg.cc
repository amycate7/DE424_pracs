/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    April 2022
          PURPOSE: Factor for an affine transformation of a Multivariate Gaussian

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

#include "affinemvg.hpp"
#include "sortindices.hpp"
#include "vecset.hpp"

using namespace std;
using namespace emdw;

//======================== Constructors and Destructor ========================

AffineMVG::AffineMVG(
  const RVIds& theXIds,
  const RVIds& theYIds,
  const prlite::ColMatrix<double>& theAMat,
  const prlite::ColVector<double>& theCVec,
  const prlite::ColMatrix<double>& theNoiseR,
  const emdw::RVIds& theObsvIds,
  const emdw::RVVals& theObsvVals)
    : xIds(theXIds),
      yIds(theYIds),
      aMat(theAMat),
      cVec(theCVec),
      noiseR(theNoiseR),
      obsvIds(theObsvIds),
      obsvVals(theObsvVals){

  // check sorting order

  auto iPtr = xIds.begin(); // check the x ids
  RVIdType prevId = *iPtr;
  while ( ++iPtr != xIds.end() ) {
    PRLITE_ASSERT(*iPtr > prevId, "X Ids in unsorted order: " << prevId << " vs " << *iPtr);
    prevId = *iPtr;
  } // while

  for (auto el: yIds) { // check the y ids
    PRLITE_ASSERT(el > prevId, "Y-Ids in unsorted order: " << prevId << " vs " << el);
    prevId = el;
  } // for

  // check the observed ids
  unsigned nObsvs = obsvIds.size();
  if ( nObsvs ) {
    PRLITE_ASSERT(obsvIds[0] >= yIds[0], "Can only observe Y Ids: " << obsvIds[0]);
    PRLITE_ASSERT(obsvIds[nObsvs-1] <= yIds[yIds.size()-1], "Unknown observed variable " << obsvIds[nObsvs-1]);
    for (unsigned k = 1; k < obsvIds.size(); k++) {
      PRLITE_ASSERT(obsvIds[k] > obsvIds[k-1], "Observed variable ids must be presorted" << obsvIds);
    } // for i
  } // if

  // exclude the observed y ids
  unsigned nVars = xIds.size() + yIds.size() - obsvIds.size();
  allVars.reserve(nVars);
  allVars.insert( allVars.end(), xIds.begin(), xIds.end() );

  // this code relies on both yIds and obsvVals being presorted.
  auto obsvIdPtr = obsvIds.begin();
  for (auto el : yIds) { // now we add the unobserved y ids

    if ( obsvIdPtr == obsvIds.end() ) {
      allVars.push_back(el);
    } // if
    else if (el != *obsvIdPtr) {
      allVars.push_back(el);
    } // else if
    else {
      obsvIdPtr++;
    } // else

  } // for

} // class specific ctor

//========================= Operators and Conversions =========================

bool AffineMVG::operator==(const AffineMVG& d) const {

    if (xIds != d.xIds) {return false;}
    if (yIds != d.yIds) {return false;}

    if (aMat != d.aMat) {return false;}
    if (cVec != d.cVec) {return false;}
    if (noiseR != d.noiseR) {return false;}
    if (obsvIds != d.obsvIds) {return false;}
    if (obsvVals != d.obsvVals) {return false;}

    // if xPtr exists we should check that too
    if (xPtr) {
      if (!d.xPtr) {return false;}
      if (*xPtr != *d.xPtr) {return false;}
    } // if

    return true;
  } // operator==

//=========================== Inplace Configuration ===========================

unsigned AffineMVG::classSpecificConfigure(
  const RVIds& theXIds,
  const RVIds& theYIds,
  const prlite::ColMatrix<double>& theAMat,
  const prlite::ColVector<double>& theCVec,
  const prlite::ColMatrix<double>& theNoiseR,
  const emdw::RVIds& theObsvIds,
  const emdw::RVVals& theObsvVals){

  this->~AffineMVG(); // Destroy existing
  new(this) AffineMVG(theXIds, theYIds, theAMat, theCVec, theNoiseR, theObsvIds, theObsvVals);
  return 1;
} // classSpecificConfigure

ostream& AffineMVG::txtWrite(ostream& file) const {
  file << "AffineMVG_IOv0" << endl;
  file << "xIds: " << xIds << endl;
  file << "yIds: " << yIds << endl;
  file << "aMat: " << aMat << endl;
  file << "cVec: " << cVec << endl;
  file << "noiseR: " << noiseR << endl;
  unsigned nObsv = obsvIds.size();
  file << "observed: " << nObsv << endl;
  for (unsigned idx = 0; idx < nObsv; idx++) {
    file << obsvIds[idx] << " " << double(obsvVals[idx]) << endl;
  } // for
  return file;
} // txtWrite

istream& AffineMVG::txtRead(istream& file) {

  RVIds theXIds;
  RVIds theYIds;
  prlite::ColMatrix<double> theAMat;
  prlite::ColVector<double> theCVec;
  prlite::ColMatrix<double> theNoiseR;

  string dummy;
  file >> dummy;
  PRLITE_ASSERT(dummy == "AffineMVG_IOv0", "Unknown AffineMVG format " << dummy);

  file >> dummy;
  PRLITE_CHECK(file >> theXIds, IOError, "Error reading x-vars RowVector.");
  file >> dummy;
  PRLITE_CHECK(file >> theYIds, IOError, "Error reading y-vars RowVector.");
  file >> dummy;
  PRLITE_CHECK(file >> theAMat, IOError, "Error reading the A matrix ColMatrix.");
  file >> dummy;
  PRLITE_CHECK(file >> theCVec, IOError, "Error reading the C vector ColVector.");
  file >> dummy;
  PRLITE_CHECK(file >> theNoiseR, IOError, "Error reading the noiseR matrix ColMatrix.");
  file >> dummy;
  unsigned nObsv;
  file >> dummy >> nObsv;
  RVIds theObsvIds;
  RVVals theObsvVals;

  if (nObsv) {
    theObsvIds.reserve(nObsv);
    theObsvVals.reserve(nObsv);
    RVIdType theId; double theVal;
    for (unsigned cnt = 0; cnt < nObsv; cnt++) {
      file >> theId >> theVal;
      theObsvIds.push_back(theId);
      theObsvVals.push_back(theVal);
    } // for
  } // if

  classSpecificConfigure(theXIds, theYIds, theAMat, theCVec, theNoiseR, theObsvIds, theObsvVals);
  return file;
} // txtRead

//===================== Required Virtual Member Functions =====================

AffineMVG* AffineMVG::copy(
  const RVIds& newVars,
  bool presorted) const {

  if ( newVars.size() ){ // replace the ids if requested
    RVIds theXIds,theYIds,theObsvIds;
    theXIds.reserve( xIds.size() ); theYIds.reserve( yIds.size() ); theObsvIds.reserve( obsvIds.size() );
    PRLITE_ASSERT(presorted, "At this stage we only allow AffineMVG construction with sorted X and Y RV ids");
    PRLITE_ASSERT(newVars.size() == xIds.size() + yIds.size(), "The number of new ids does not match that in the original object");

    auto obsItr = obsvIds.begin();
    for (unsigned k = 0; k < newVars.size(); k++) {

      if ( k < xIds.size() ) { // these are the new xIds;
        theXIds.push_back(newVars[k]);
      } // if

      else { // so these must be yIds
        theYIds.push_back(newVars[k]);
        if (yIds[k] == *obsItr) {
          theObsvIds.push_back(newVars[k]);
          obsItr++;
        } // if
      } // else

    } // for

    return new AffineMVG(theXIds, theYIds, aMat, cVec, noiseR, theObsvIds, obsvVals);
  } // if

  else { // same ids, just make a copy
    return new AffineMVG(*this);
  } // else

} //copy

///
unsigned AffineMVG::noOfVars() const {
  return allVars.size();
} // noOfVars

/// returns the ids of the variables the factor is defined on
RVIds AffineMVG::getVars() const {
  return allVars;
} // getVars

///
RVIdType AffineMVG::getVar(unsigned varNo) const {
  return allVars[varNo];
} // getVar

//===================== Required Factor Operations ============================

//normally these would be boilerplate that uses DDD to select the
//underlying FactorOperator. However, for the moment being we are
//going to short-cut directly implement them here.

//------------------Family 0

// returning a sample from the distribution
RVVals AffineMVG::sample(Sampler* procPtr) const{
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");
  PRLITE_ASSERT(xPtr, "It does not really make sense to sample when there is no xPtr");

  if ( obsvIds.size() ) { // a bit expensive though
    rcptr<Factor> tmpPtr = xyPtr->observeAndReduce(obsvIds, obsvVals);
    return tmpPtr->Factor::sample();
  } // if

  return xyPtr->Factor::sample();
} //sample

//------------------Family 1

// returning new factor
uniqptr<Factor> AffineMVG::normalize(FactorOperator* procPtr) const {
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");

  AffineMVG* ret = new AffineMVG(*this);
  ret->inplaceNormalize(procPtr);

  return uniqptr<AffineMVG>(ret);
} // normalize

// inplace
void AffineMVG::inplaceNormalize(FactorOperator* procPtr){
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");
  if (xPtr) {
    xPtr->inplaceNormalize();
    xyPtr->inplaceNormalize();
  } // if

} //inplaceNormalize

//------------------Family 2

void AffineMVG::inplaceAbsorb(const Factor* rhsFPtr,
                              FactorOperator* procPtr){
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");

  const  SqrtMVG* rhsPtr = dynamic_cast<const SqrtMVG*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " AffineMVG_InplaceAbsorb : rhs not a valid SqrtMVG : "
                << typeid(*rhsFPtr).name() );

  // is it an x or a y? Actually we should also be able to handle
  // anything in between, but for the sake of simplicity we don't.

  if (rhsPtr->getVars() == xIds) { // ok, we are multiplying over x
    if (xPtr) { // seems we have a previous instantiation
      xPtr->inplaceAbsorb(rhsFPtr);
    } // if
    else{
      xPtr = uniqptr<SqrtMVG>( rhsPtr->copy() );
    } // else

    xyPtr = uniqptr<SqrtMVG>( SqrtMVG::constructAffineGaussian(*xPtr, aMat, cVec, yIds, noiseR) );
  } // if

  else if (rhsPtr->getVars() == yIds) { // ok, then it should be y we are multiplying over.
    if (xyPtr) {xyPtr->inplaceAbsorb(rhsFPtr);} // if
    // do nothing if xy does not exist at this point.
  } // else if

  else{
    PRLITE_FAIL( "The distribution to multiply with does not match x or y" << rhsPtr->getVars() );
  } // else

  // we will handle obsvY when marginalising

} //inplaceAbsorb

// returning a new factor
uniqptr<Factor> AffineMVG::absorb(const Factor* rhsFPtr,
                                FactorOperator* procPtr) const {
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");

  AffineMVG* ret = new AffineMVG(*this);
  ret->inplaceAbsorb(rhsFPtr, procPtr);
  return uniqptr<AffineMVG>(ret);
} //absorb


void AffineMVG::inplaceCancel(const Factor* rhsFPtr,
                              FactorOperator* procPtr){
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");
  PRLITE_ASSERT(xPtr, "It does not really make sense to cancel when there is no xPtr");

  const  SqrtMVG* rhsPtr = dynamic_cast<const SqrtMVG*>(rhsFPtr);
  PRLITE_ASSERT(rhsPtr, " AffineMVG_InplaceCancel : rhs not a valid SqrtMVG : "
                << typeid(*rhsFPtr).name() );


  // is it an x or a y? Actually we should also be able to handle
  // anything in between, but for the sake of simplicity we don't.

  if (rhsPtr->getVars() == xIds) { // ok, we are dividing with x distribution
    xPtr->inplaceCancel(rhsFPtr);
    xyPtr = uniqptr<SqrtMVG>( SqrtMVG::constructAffineGaussian(*xPtr, aMat, cVec, yIds, noiseR) );
  } // if

  else if (rhsPtr->getVars() == yIds) { // ok, then it should be y distributiohn we are dividing with.
    xyPtr->inplaceCancel(rhsFPtr);
    // do nothing if xy does not exist at this point.
  } // else if

  else{
    PRLITE_FAIL( "The distribution to divide with does not match x or y" << rhsPtr->getVars() );
  } // else

  // we will handle obsvY when marginalising

} //inplaceCancel

// returning a new factor
uniqptr<Factor> AffineMVG::cancel(const Factor* rhsFPtr,
                                FactorOperator* procPtr) const {
  AffineMVG* ret = new AffineMVG(*this);
  ret->inplaceCancel(rhsFPtr, procPtr);
  return uniqptr<AffineMVG>(ret);
} //cancel

  //------------------Family 3 eg Marginalize

uniqptr<Factor> AffineMVG::marginalize(const RVIds& variablesToKeep,
                                          bool presorted,
                                          const Factor* peekAheadPtr,
                                          FactorOperator* procPtr) const {
  PRLITE_ASSERT(!procPtr, "Currently there is not support for multiple FactorOperators");

  // fiddling about to avoid copying when no sorting is required
  const RVIds* sortedVarsPtr = &variablesToKeep;
  RVIds tmpVars; // need it here to not go out of scope too soon
  vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( variablesToKeep, less<unsigned>() );
    tmpVars = extract<unsigned>(variablesToKeep,sorted);
    sortedVarsPtr = &tmpVars;
  } // if
  const RVIds& sortedVars(*sortedVarsPtr);

  // a bit expensive, but we have no guarantee that there are no extra variables.
  RVIds l2i;
  RVIds r2i;
  // this already takes obsvY into account.
  RVIds vecI = sortedIntersection(allVars, sortedVars, l2i, r2i);
  PRLITE_ASSERT(vecI.size(), "Attempting to marginalize over unknown variables " << variablesToKeep);

  // if no joint yet, return a vacuous distribution over the requested variables
  if (!xyPtr) {
    return make_unique<SqrtMVG>( SqrtMVG(vecI) );
  } // if

  else { // have a joint instantiated

    if ( obsvIds.size() ) {
      return xyPtr->Factor::observeAndReduce(obsvIds, obsvVals)->
        marginalize(vecI)->
        normalize();
    } // if

    else {
      return xyPtr->
        Factor::marginalize(vecI)->
        normalize();
    } // else

  } // else


} // marginalize

  //------------------Family 4: eg ObserveAndReduce

uniqptr<Factor> AffineMVG::observeAndReduce(const RVIds& variables,
                                            const RVVals& assignedVals,
                                            bool presorted,
                                            FactorOperator* procPtr) const {
  PRLITE_FAIL("observeAndReduce not implemented for AffineMVG");
  uniqptr<Factor> tmpPtr; // dummy code to stop warnings
  return tmpPtr;
} // observeAndReduce

  //------------------Family 5: eg inplaceDampen

  /**
   * NOTE the parameter is the OLD message and will be weighted with df
   */
double AffineMVG::inplaceDampen(const Factor* oldMsg,
                                double df,
                                FactorOperator* procPtr){
  PRLITE_FAIL("inplaceDampen not implemented for AffineMVG");
} // inplaceDampen

//====================== Other Virtual Member Functions =======================

double AffineMVG::potentialAt(
  const RVIds& theVars,
  const RVVals& theirVals,
  bool presorted) const{

  PRLITE_ASSERT(xyPtr, "Cannot determine the height of a non-existant pdf");
  return xyPtr->potentialAt(theVars, theirVals, presorted);
} // potentialAt

double AffineMVG::distance(const Factor* rhsPtr) const {
  PRLITE_ASSERT(xyPtr, "Cannot compare the distance with a non-existant pdf");
  return xyPtr->distance(rhsPtr);
} // distance


//================================== Friends ==================================

// /// output
// ostream& operator<<( ostream& file,
//                           const AffineMVG& n ){
//   return n.txtWrite(file);
// } // operator<<

// /// input
// istream& operator>>( istream& file,
//                           AffineMVG& n ){
//   return n.txtRead(file);
// } // operator>>
