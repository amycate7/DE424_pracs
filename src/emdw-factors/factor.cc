/*
 * Author     : JA du Preez
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
// standard headers
#include <string>  // string
#include <iostream>  // cout, endl
#include <typeinfo>

#include "factor.hpp"
#include "sortindices.hpp"

using namespace std;

//========================= Operators and Conversions =========================

  /**
   * equality.
   */
  bool Factor::operator==(const Factor& rhs) const {
    return isEqual(&rhs);
  } // operator==

  /**
   * inequality.
   */
  bool Factor::operator!=(const Factor& a) const {
    return !operator==(a);
  } // operator!=

//=========================== Inplace Configuration ===========================
//===================== Required Virtual Member Functions =====================

double
Factor::distanceFromVacuous() const {
  rcptr<Factor> vacPtr( vacuousCopy() );
  vacPtr->inplaceNormalize(); // just to be sure
  return distance(vacPtr);
} // distanceFromVacuous

emdw::RVVals
Factor::sample(
  Sampler* procPtr) const{
  PRLITE_ASSERT(false, std::string(typeid(*this).name()) + " does not (yet) support the sample member function");
} // sample

//====================== Other Virtual Member Functions =======================

double
Factor::potentialAt(
  const emdw::RVIds&,
  const emdw::RVVals&,
  bool) const {
  PRLITE_ASSERT(false, std::string(typeid(*this).name()) + " does not (yet) support the potentialAt member function");
} // potentialAt

double
Factor::distance(
  const Factor* rhsPtr) const{
  PRLITE_ASSERT(
    false,
    std::string(typeid(*this).name()) +
    " does not have the (deprecated) 'distance' member function. Check if 'distanceFromVacuous' is not what you require.  A default implementation is available from 'factor.cc'.");
} // distance

double
Factor::distance(
  const rcptr<Factor>& rhsPtr) const {
  return distance( rhsPtr.get() );
} // distance

//======================= Non-virtual Member Functions ========================

uniqptr<Factor> Factor::absorb(const std::vector< rcptr<Factor> >& facPtrs,
                               const rcptr<FactorOperator>& procPtr) const{
  uniqptr<Factor> resultPtr = uniqptr<Factor>( copy() );
  for (unsigned i = 0; i < facPtrs.size(); i++){
    resultPtr = resultPtr->absorb( facPtrs[i].get(), procPtr.get() );
  } // for
  return resultPtr;
} // apply

// these ones are just to make working with rcptrs easier

//------------------Family 0: eg Sample

emdw::RVVals Factor::sample(const rcptr<Sampler>& p1) const {
  return sample( p1.get() );
} // sample

//------------------Family 1: eg Normalize

uniqptr<Factor> Factor::normalize(const rcptr<FactorOperator>& p1) const {
  return normalize(p1.get() );
} // normalize

void Factor::inplaceNormalize(const rcptr<FactorOperator>& p1) {
  inplaceNormalize( p1.get() );
} // inplaceNormalize

//------------------Family 2: eg Product

uniqptr<Factor> Factor::absorb(const rcptr<Factor>& p2,
                              const rcptr<FactorOperator>& p1) const {
  return absorb( p2.get(), p1.get() );
} // absorb

void Factor::inplaceAbsorb(const rcptr<Factor>& p2,
                          const rcptr<FactorOperator>& p1) {
  inplaceAbsorb( p2.get(), p1.get() );
} // inplaceAbsorb

uniqptr<Factor> Factor::cancel(const rcptr<Factor>& p2,
                              const rcptr<FactorOperator>& p1) const {
  return cancel( p2.get(), p1.get() );
} // cancel

void Factor::inplaceCancel(const rcptr<Factor>& p2,
                          const rcptr<FactorOperator>& p1) {
  inplaceCancel( p2.get(), p1.get() );
} // inplaceCancel

//------------------Family 3: eg Marginalize

uniqptr<Factor> Factor::marginalize(
  const emdw::RVIds& variablesToKeep,
  bool presorted,
  const Factor* peekAheadPtr,
  const rcptr<FactorOperator>& procPtr) const {
  return marginalize( variablesToKeep, presorted, peekAheadPtr, procPtr.get() );
} //marginalize

uniqptr<Factor> Factor::marginalizeX(
  const emdw::RVIds& toZap,
  bool presorted,
  const Factor* peekAheadPtr,
  const rcptr<FactorOperator>& procPtr) const {

  if (toZap.size() < 2) {presorted = true;} // if

  // fiddling about to avoid copying when no sorting is required
  const emdw::RVIds* sortedZapPtr = &toZap;
  emdw::RVIds tmpVars; // need it here to not go out of scope too soon
  std::vector<size_t> sorted; // need it here to not go out of scope too soon
  if (!presorted) {
    sorted = sortIndices( toZap, std::less<unsigned>() );
    tmpVars = extract<unsigned>(toZap,sorted);
    sortedZapPtr = &tmpVars;
  } // if
  const emdw::RVIds& vecR(*sortedZapPtr);
  const emdw::RVIds& vecL( getVars() );
  emdw::RVIds toKeep;
  toKeep.reserve( vecL.size() );

  auto iL = vecL.begin();
  auto iR = vecR.begin();
  while ( iL != vecL.end() && iR != vecR.end() ){

    if (*iL < *iR) { // only in vecL
      toKeep.push_back(*iL++);
    } // if

    else if (*iL > *iR) {
      iR++;
    } // else if

    else {                            //common to both vecs
      iL++; iR++;
    } // else

  } // while

  if ( iL != vecL.end() ){         // some extra vars in vecL
    toKeep.insert( toKeep.end(), iL, vecL.end() );
  } // if

  return marginalize(toKeep, true, peekAheadPtr, procPtr);
} // marginalizeX

//------------------Family 4: eg Observe

uniqptr<Factor> Factor::observeAndReduce(const emdw::RVIds& variables,
                                         const emdw::RVVals& assignedVals,
                                         bool presorted,
                                         const rcptr<FactorOperator>& procPtr) const {
  return observeAndReduce( variables, assignedVals, presorted, procPtr.get() );
} //observeAndReduce

//------------------Family 5: eg inplaceDampen

double Factor::inplaceDampen(const rcptr<Factor>& oldMsg,
                             double df,
                             const rcptr<FactorOperator>& procPtr){
  return inplaceDampen( oldMsg.get(), df, procPtr.get() );
} // inplaceDampen

//================================== Friends ==================================

/// output
std::ostream& operator<<( std::ostream& file,
                          const Factor& n ){
  return n.txtWrite(file);
} // operator<<


/// input
std::istream& operator>>( std::istream& file,
                          Factor& n ){
  return n.txtRead(file);
} // operator>>

//================================== External =================================

uniqptr<Factor> absorb(const std::vector< rcptr<Factor> >& facPtrs,
                       const rcptr<FactorOperator>& procPtr) {
  uniqptr<Factor> resultPtr = uniqptr<Factor>( facPtrs[0]->copy() );
  for (unsigned i = 1; i < facPtrs.size(); i++){
    resultPtr = resultPtr->absorb( facPtrs[i].get(), procPtr.get() );
  } // for
  return resultPtr;
} // apply
