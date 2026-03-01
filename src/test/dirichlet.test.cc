/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    March 2018
          PURPOSE: Unit tests for Dirichlet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


// patrec headers
#include "prlite_testing.hpp"

// emdw headers
#include "dirichlet.hpp"

// standard headers
#include <iostream>

using namespace std;
using namespace emdw;

// v1 and v2 ito alphas
double
symmKLDist(
  vector<double> v1,
  vector<double> v2) {

  // example values for v1{1, 2, 3, 2, 5}, v2{2, 3, 3, 2, 2}  // alphas

  double sumAlpha1 = 0.0; for (auto e : v1) sumAlpha1 += e; // 13
  double sumAlpha2 = 0.0; for (auto e : v2) sumAlpha2 += e; // 12
  double dist = 0.0;

  for (unsigned k = 0; k < v1.size(); k++) {
    dist += (v1[k] - v2[k])*
      (digamma(v1[k])-digamma(sumAlpha1) - (digamma(v2[k])-digamma(sumAlpha2)));
  } // for
  return dist/2.0;
} // symmKLDist

class BaseTestDirichlet : public prlite::TestCase {
public:

  /** setup to be run before each test. */
  virtual void setup() {} // setup


  /** clean up after. */
  virtual void tearDown() {} // tearDown

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

}; // BaseTestDirichlet

class TestDirichlet : public BaseTestDirichlet {
public:

  /** a name. */
  const string& isA() {
    static const string name = "TestDirichlet";
    return name;
  } // isA

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);
    try {

#ifdef TESTING_MODE
    FILEPOS;
#endif

      // #####################################
      // ### define some generic variables ###
      // #####################################


      typedef unsigned T;
      typedef Dirichlet<T> DIR;

      rcptr< vector<T> > binDom (
        new vector<T>{0,1});
      rcptr< vector<T> > cat5Dom (
        new vector<T>{0,1,2,3,4});

      const double tol = std::numeric_limits<double>::epsilon()*10;
      double df = 0.1;

      // lets first find the ref distance using vectors of phantoms
      vector<double> v1{1, 2, 3, 2, 5};
      vector<double> v2{2, 3, 3, 2, 2};
      double refDist5 =  symmKLDist(v1, v2); // 2.3333

      vector<double> refAlpha5( v1.size() );  // 1.1 2.1 3 2 4.7
      for (unsigned k = 0; k < v1.size(); k++) {
        refAlpha5[k] = (1-df)*v1[k] + df*v2[k];
      } // for

      v1 = {3, 3};
      v2 = {2, 2};
      double refDist2 =  symmKLDist(v1, v2); // 0.05
      vector<double> refAlpha2( v1.size() );  // 2.9 2.9
      for (unsigned k = 0; k < v1.size(); k++) {
        refAlpha2[k] = (1-df)*v1[k] + df*v2[k];
      } // for

      vector<double> refAlpha, tstAlpha;
      rcptr<Factor> tstPtr;
      rcptr<DIR> tstDirPtr, dirPtr1, dirPtr2;

      // #######################################
      // ### Lets test some member functions ###
      // #######################################

      // +++++++++++++++++++++++
      // +++ sumOfAlphaCount +++
      // +++++++++++++++++++++++

      dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {0, 1.0}, {2, 3.0}, {4, 5.0} }) );

      if ( dirPtr1->sumAlphas() != 13.0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::sumAlphas test FAILED!");
        PRLITE_LOGERROR( PRLITE_pLog, "Testing", "Expected 13, but got " << dirPtr1->sumAlphas() );
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // if ( dirPtr1->prodAlphas() != 60.0) {
      //   LOGERROR(pLog, "Testing", "DirichletSet::prodAlphas FAILED!");
      //   LOGERROR( pLog, "Testing", "Expected 60.0, but got " << dirPtr1->prodAlphas() );
      //   THROW(prlite::TestFailed, "");
      // } // if

      // +++++++++++++++++++
      // +++ addToCounts +++
      // +++++++++++++++++++

      // --- fully default ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 3.0) );
      dirPtr1->addToCounts({1}, 0.5);

      refAlpha  = {3.0, 3.5};
      tstAlpha.clear();
      for (unsigned k = 0; k < dirPtr1->noOfCategories(); k++) {
        tstAlpha.push_back( dirPtr1->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // +++++++++++++++++++++++++++
      // +++ distanceFromVacuous +++
      // +++++++++++++++++++++++++++

      {
        double refDist = symmKLDist({1, 2, 3, 2, 5}, {1,1,1,1,1});  // alphas
        //cout << "refDist: " << refDist << endl;

        dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {0, 1.0}, {2, 3.0}, {4, 5.0} }) );
        double tstDist = dirPtr1->distanceFromVacuous();
        //cout << "tstDist: " << tstDist << endl;

        if (std::abs(tstDist-refDist) > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::distance FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist
                   << " but got " << tstDist);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      }

      // ++++++++++++++++
      // +++ distance +++
      // ++++++++++++++++

      dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {0, 1.0}, {2, 3.0}, {4, 5.0} }) );
      dirPtr2 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {1, 3.0}, {2, 3.0} }) );

      double tstDist = dirPtr1->distance( dirPtr2.get() );
      if (std::abs(tstDist-refDist5) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist5
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // This one is supposed to be the same as the original
      dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 1.0, { {1, 2.0}, {2, 3.0}, {3, 2.0}, {4, 5.0} }) );

      tstDist = dirPtr1->distance( dirPtr2.get() );
      if (std::abs(tstDist-refDist5) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist5
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // distance between distrib and itself should be zero
      tstDist = dirPtr1->distance( dirPtr1.get() );
      if (std::abs(tstDist) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist5
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ######################################
      // ### lets test the factor operators ###
      // ######################################


      // +++++++++++++++++
      // +++ normalize +++
      // +++++++++++++++++

      tstPtr = dirPtr1->Factor::normalize();
      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::normalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Normalize should leave the factor unchanged");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // +++++++++++++++++++
      // +++ marginalize +++
      // +++++++++++++++++++

      tstPtr = dirPtr1->Factor::marginalize({0u});
      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Marginalize should leave the factor unchanged");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++
      // +++ absorb +++
      // ++++++++++++++

      // --- both fully default ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 3.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      refAlpha  = {4.0, 4.0};
      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both fully sparse ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0, 3.0},{1, 3.0}}) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0, 2.0},{1, 2.0}}) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      refAlpha  = {4.0, 4.0};
      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully default values, second fully sparse values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 3.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,2.0},{1,2.0}}) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      refAlpha  = {4.0, 4.0};
      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully sparse values, second fully default values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,3.0},{1,3.0}}) );
      // phantoms == 1, fully default
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      refAlpha  = {4.0, 4.0};
      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both partly sparse ---

      dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {0, 1.0}, {2, 3.0}, {4, 5.0} }) );
      dirPtr2 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {1, 3.0}, {2, 3.0} }) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      refAlpha  = {2,4,5,3,6};
      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both very sparse ---

      dirPtr1 = uniqptr<DIR>( new DIR(0u, binDom, 0.0, { {1, 4.0} }) );
      dirPtr2 = uniqptr<DIR>( new DIR(0u, binDom, 0.0, { {1, 1.0} }) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if ( tstDirPtr->getAlpha(0) != -1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected -1 but got " << tstDirPtr->getAlpha(0) << "\n" << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      if ( tstDirPtr->getAlpha(1) != 4) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 4 but got " << tstDirPtr->getAlpha(1) << "\n" << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++
      // +++ cancel +++
      // ++++++++++++++

      // --- both fully default ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(dirPtr2)->normalize();

      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::cancel FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected : " << *dirPtr1 << ", but got " << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both fully sparse ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 0.0, map<unsigned,double>{{0,2.0},{1,2.0}}) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 0.0, map<unsigned,double>{{0,1.0},{1,1.0}}) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(dirPtr2)->normalize();

      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::cancel FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected : " << *dirPtr1 << ", but got " << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully default values, second fully sparse values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 0.0, map<unsigned,double>{{0,1.0},{1,1.0}}) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(dirPtr2)->normalize();

      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::cancel FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected : " << *dirPtr1 << ", but got " << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully sparse values, second fully default values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 0.0, map<unsigned,double>{{0,2.0},{1,2.0}}) );
      // phantoms == 1, fully default
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(dirPtr2)->normalize();

      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::cancel FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected : " << *dirPtr1 << ", but got " << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both partly sparse ---

      dirPtr1 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {0, 1.0}, {2, 3.0}, {4, 5.0} }) );
      dirPtr2 = uniqptr<DIR>( new DIR(0u, cat5Dom, 2.0, { {1, 3.0}, {2, 3.0} }) );

      tstPtr = dirPtr1->Factor::absorb(dirPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(dirPtr2)->normalize();

      if (*tstPtr != *dirPtr1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::cancel FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected : " << *dirPtr1 << ", but got " << *tstDirPtr);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // +++++++++++++++++++++
      // +++ InplaceDampen +++
      // +++++++++++++++++++++

      // --- both fully default ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 3.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );

      tstPtr = uniqptr<Factor>( dirPtr1->copy() );
      tstDist = tstPtr->Factor::inplaceDampen(dirPtr2, df);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if (std::abs(tstDist-refDist2) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist2
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen alphas FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha2 << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both fully sparse ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,3.0},{1,3.0}}) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,2.0},{1,2.0}}) );

      tstPtr = uniqptr<Factor>( dirPtr1->copy() );
      tstDist = tstPtr->Factor::inplaceDampen(dirPtr2, df);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if (std::abs(tstDist-refDist2) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist2
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen alphas FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha2 << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully default values, second fully sparse values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 3.0) );
      // phantoms == 1
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,2.0},{1,2.0}}) );

      tstPtr = uniqptr<Factor>( dirPtr1->copy() );
      tstDist = tstPtr->Factor::inplaceDampen(dirPtr2, df);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if (std::abs(tstDist-refDist2) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist2
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen alphas FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha2 << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- first dirichlet fully sparse values, second fully default values ---

      // phantoms == 2
      dirPtr1 = uniqptr<DIR> ( new DIR(0u, binDom, 1.0, map<unsigned,double>{{0,3.0},{1,3.0}}) );
      // phantoms == 1, fully default
      dirPtr2 = uniqptr<DIR> ( new DIR(0u, binDom, 2.0) );

      tstPtr = uniqptr<Factor>( dirPtr1->copy() );
      tstDist = tstPtr->Factor::inplaceDampen(dirPtr2, df);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if (std::abs(tstDist-refDist2) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist2
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen alphas FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha2 << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- both partly sparse ---

      dirPtr1 = uniqptr<DIR>(  // phantoms: 0,1,2,1,4
        new DIR(0u, cat5Dom, 1.0, map<unsigned,double>{{1,2.0},{2,3.0},{3,2.0},{4,5.0}}) );
      dirPtr2 = uniqptr<DIR>( // phantoms: 1,2,2,1,1
        new DIR(0u, cat5Dom, 2.0, map<unsigned,double>{{1,3.0},{2,3.0}}) );

      tstPtr = uniqptr<Factor>( dirPtr1->copy() );
      tstDist = tstPtr->Factor::inplaceDampen(dirPtr2, df);
      tstDirPtr = std::dynamic_pointer_cast< Dirichlet<T> > (tstPtr);

      if (std::abs(tstDist-refDist5) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen distance FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected distance = " << refDist2
                 << " but got " << tstDist);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstAlpha.clear();
      for (unsigned k = 0; k < tstDirPtr->noOfCategories(); k++) {
        tstAlpha.push_back( tstDirPtr->getAlpha({k}) );
      } // for
      if ( tstAlpha != refAlpha5) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet::inplaceDampen alphas FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << refAlpha2 << ", but got " << tstAlpha);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

#ifdef TESTING_MODE
    FILEPOS;
#endif

    } // try

    catch (const char* msg) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Dirichlet threw an exception!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", msg);
      PRLITE_THROW(prlite::TestFailed, "");
    } // catch


  } // runTests

}; // TestDirichlet

PRLITE_REGISTER_TESTCASE(TestDirichlet);
