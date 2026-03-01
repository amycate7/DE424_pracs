/*******************************************************************************

          AUTHOR: JA du Preez
          DATE:    March 2015
          PURPOSE: DiscreteTable PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


// patrec headers
#include "prlite_testing.hpp"

// emdw headers
#include "discretetable.hpp"
#include "cat1expf.hpp"

// standard headers
#include <iostream>

using namespace std;
using namespace emdw;

typedef unsigned AssignType;

#include "combinations.hpp"

/*
 * spoed toetse: (op lbu cnt 2000)
 *
 * ClockTicks:
 *
 * Tegniek      normalize marginalize absorb cancel dampen
 * SpTable      6         20          24     2      0
 * D0Table      4         12          31     4      0
 * D1Table    528         95        2008     6      0
 *
 * Total table sizes: (hierdie een wys die eintlike probleem)
 *
 * Def,floor == 0
-9.30411e-05 Size: 835443
29999.8 Size: 200372
999.994 Size: 25492
999.866 Size: 9823
996.875 Size: 3241
499.991 Size: 660
0.00018016 Size: 579
1.16276e-06 Size: 579
2.35875e-06 Size: 579
2.36592e-11 Size: 579

* Def,floor > 0
-0.275552 Size: 835443
439.287 Size: 635347
465.596 Size: 1530725
117.877 Size: 1891401
105.676 Size: 675841
85.3008 Size: 1978849
63.7529 Size: 1383961
130.349 Size: 1396122
98.5388 Size: 1898989
456.236 Size: 2177023

 */

double
DiscreteTable2SidedKLDist(
  vector<double> p,
  vector<double> q) {

  double dist = 0.0;
  for (unsigned k = 0; k < p.size(); k++) {
    double P = p[k];
    double logP = P > 0.0 ? log(P) : -1000.0;
    double Q = q[k];
    double logQ = Q > 0.0 ? log(Q) : -1000.0;

    if (P > 0.0) { dist += P*logP - P*logQ; } // if
    if (Q > 0.0) { dist += Q*logQ - Q*logP; } // if
  } // for

  return dist/2;
} // DiscreteTable2SidedKLDist

double
DiscreteTable1SidedKLDist(
  vector<double> p,
  vector<double> q) {

  double dist = 0.0;
  for (unsigned k = 0; k < p.size(); k++) {
    double P = p[k];
    double logP = P > 0.0 ? log(P) : -1000.0;
    double Q = q[k];
    double logQ = Q > 0.0 ? log(Q) : -1000.0;

    if (P > 0.0) { dist += P*logP - P*logQ; } // if
    //if (Q > 0.0) { dist += Q*logQ - Q*logP; } // if
  } // for

  return dist;
} // DiscreteTable1SidedKLDist

class BaseTestDiscreteTable : public prlite::TestCase {
public:

  /** setup to be run before each test. */
  virtual void setup() {} // setup


  /** clean up after. */
  virtual void tearDown() {} // tearDown

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

}; // BaseTestDiscreteTable

class TestDiscreteTable : public BaseTestDiscreteTable {
public:

  /** a name. */
  const string& isA() {
    static const string name = "TestDiscreteTable";
    return name;
  } // isA

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

    try{

#ifdef TESTING_MODE
    FILEPOS;
#endif

      // #####################################
      // ### define some generic variables ###
      // #####################################

      typedef DiscreteTable<AssignType> DT;
      const double tol = 1E-6; //std::numeric_limits<double>::epsilon()*10;

      rcptr<Factor> tstPtr, refPtr;
      rcptr<DT> tstDTPtr, dtPtr1, dtPtr2;

      // ####################################
      // ### Test cases with defVals == 0 ###
      // ####################################

      {

        rcptr< vector<AssignType> > doms0(new vector<AssignType>{0});
        rcptr< vector<AssignType> > doms1(new vector<AssignType>{1,2,3});
        rcptr< vector<AssignType> > doms2(new vector<AssignType>{1,2});
        rcptr< vector<AssignType> > doms3(new vector<AssignType>{1,2});
        rcptr< vector<AssignType> > doms5(new vector<AssignType>{1,2});
        rcptr< vector<AssignType> > doms6(new vector<AssignType>{1,2});

        std::vector< rcptr<Factor> > factorPtrs;

        // factor 0 is empty
        factorPtrs.push_back( uniqptr<DT> ( new DT({0u}, {doms0},  0.00) ) );

        // factor 1 is defined on X1 and X2
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({1,2}, {doms1,doms2},  0.0,
                   {
                     {{1,1}, 0.5},
                     {{1,2}, 0.8},
                     {{2,1}, 0.1},
                     {{3,1}, 0.3},
                     {{3,2}, 0.9}
                   }, 0.0, 0.0) ) );

        // factor 2 is defined on X2  and X3
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({2,3}, {doms2,doms3},  0.0,
                   {
                     {{1,1}, 0.0},
                     {{1,2}, 0.7},
                     {{2,1}, 0.1},
                     {{2,2}, 0.2}
                   }, 0.0, 0.0) ) );

        // factor 3 has same scope as factor 2 (defined on X2  and X3)
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({2,3}, {doms2,doms3},  0.0,
                   {
                     {{1,1}, 0.0},
                     {{1,2}, 0.7},
                     {{2,1}, 0.1},
                     {{2,2}, 0.2}
                   }, 0.0, 0.0) ) );

        // factor 4 is defined on X2 and is a subset of factorPtrs 1, 2 and 3
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({2}, {doms2},  0.0,
                   {
                     {{1}, 0.0},
                     {{2}, 0.5}
                   }, 0.0, 0.0) ) );

        // factor 5 is defined on X1, X2, X3 and is equal to the product
        // of factors 1 and 2
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({1,2,3}, {doms1,doms2,doms3},  0.0,
                   {
                     {{1,1,2}, 0.35},
                     {{1,2,1}, 0.08},
                     {{1,2,2}, 0.16},
                     {{2,1,2}, 0.07},
                     {{3,1,2}, 0.21},
                     {{3,2,1}, 0.09},
                     {{3,2,2}, 0.18}
                   }, 0.0, 0.0) ) );

        // factor 6 is defined on X5  and X6
        factorPtrs.push_back(
          uniqptr<DT>(
            new DT({5,6}, {doms5,doms6},  0.0,
                   {
                     {{1,1}, 0.5},
                     {{1,2}, 0.7},
                     {{2,1}, 0.1},
                     {{2,2}, 0.2}
                   }, 0.0, 0.0) ) );

        // *********************
        // *** inplaceDampen ***
        // *********************

        {

          // lets first find the ref distance between test vectors
          vector<double> v1{0.05, 0.375, 0.125, 0.25, 0.2, 0};
          vector<double> v2{0, 0.125, 0.25, 0.2, 0.375, 0.05};
          double refDist =  DiscreteTable2SidedKLDist(v1, v2); // 0.20599
          //cout << "refDist: " << refDist << endl;

          rcptr< vector<AssignType> > doms00(
            new vector<AssignType>{0,1,2,3,4,5});
          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms00}, 0.0,
                   {
                     {{0}, 0.05},
                     {{1}, 0.375},
                     {{2}, 0.125},
                     {{3}, 0.25},
                     {{4}, 0.2}
                   }) );
          dtPtr2 = uniqptr<DT>(
            new DT({0u}, {doms00},  0.0,
                   {
                     {{1}, 0.125},
                     {{2}, 0.25},
                     {{3}, 0.2},
                     {{4}, 0.375},
                     {{5}, 0.05},
                       }) );

          double tstDist = dtPtr1->Factor::inplaceDampen(dtPtr2, 0.1);
          //cout << "DT: " << refDist << " " << tstDist << endl << *dtPtr1 << endl;

          if (std::abs(tstDist-refDist) > tol) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable:: inplaceDampen FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected distance = " << refDist
                     << " but got " << tstDist);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

          vector<double> vRef = {0.045, 0.35, 0.1375, 0.245, 0.2175, 0.0};

          for (AssignType a = 0; a < vRef.size(); a++) {
            if (std::abs(dtPtr1->potentialAt({0u},{a}) - vRef[a]) > tol) {
              PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable:: inplaceDampen FAILED!");
              PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected prob = " << vRef[a]
                       << " but got " << dtPtr1->potentialAt({0u},{a}));
              PRLITE_THROW(prlite::TestFailed, "");
            } // if
          } // for

        }

        //========== normalize ====================

        //---------- normalize an empty factor

        // an empty table is a killer in a pgm, it fractures it into islands.

        try {
          tstPtr = factorPtrs[0]->normalize();
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize failed!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "normalizing an empty factor should fail!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // try
        catch (const char* s) {
          //if we go through here the exception worked correctly
        } // catch

        //---------- normalize a non-empty factor

        tstPtr = factorPtrs[1]->normalize();
        refPtr = uniqptr<DT>(
          new DT({1,2}, {doms1,doms2},  0.0,
                 {
                   {{1,1}, 0.5/2.6},
                   {{1,2}, 0.8/2.6},
                   {{2,1}, 0.1/2.6},
                   {{2,2}, 0.0/2.6},
                   {{3,1}, 0.3/2.6},
                   {{3,2}, 0.9/2.6}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== max normalize ====================

        //---------- normalize a non-empty factor

        rcptr<FactorOperator> maxNormPtr = uniqptr<FactorOperator>( new DiscreteTable_MaxNormalize<AssignType> );
        tstPtr = factorPtrs[1]->normalize(maxNormPtr);
        refPtr = uniqptr<DT>(
          new DT({1,2}, {doms1,doms2},  0.0,
                 {
                   {{1,1}, 0.5/0.9},
                   {{1,2}, 0.8/0.9},
                   {{2,1}, 0.1/0.9},
                   {{2,2}, 0.0/0.9},
                   {{3,1}, 0.3/0.9},
                   {{3,2}, 0.9/0.9}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable_MaxNormalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== product ======================

        //---------- two factorPtrs with identical scope

        tstPtr = factorPtrs[2]->absorb(factorPtrs[3]);
        refPtr = uniqptr<DT>(
          new DT({2,3}, {doms2,doms3},  0.0,
                 {
                   {{1,1}, 0.0},
                   {{1,2}, 0.49},
                   {{2,1}, 0.01},
                   {{2,2}, 0.04}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- two factorPtrs with identical scope no 2

        rcptr< vector<AssignType> > dm4(new vector<AssignType>{0,1,2,3});

        rcptr<Factor> f1Ptr = uniqptr<DT>(
          new DT({0}, {dm4},  0.0,
                 {
                   {{0}, 1.0},
                   {{1}, 1.0},
                   {{2}, 0.0},
                   {{3}, 0.0}
                 }, 0.0, 0.0) );

        rcptr<Factor> f2Ptr = uniqptr<DT>(
          new DT({0}, {dm4},  0.0,
                 {
                   {{0}, 1.0},
                   {{1}, 1.0},
                 }, 0.0, 0.0) );

        f1Ptr->inplaceAbsorb(f2Ptr);
        if (*f1Ptr !=  *f2Ptr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *f2Ptr << "\nGot:\n" << *f1Ptr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- two factorPtrs with identical scope (flipped)

        tstPtr = factorPtrs[3]->absorb(factorPtrs[2]);
        refPtr = uniqptr<DT>(
          new DT({2,3}, {doms2,doms3},  0.0,
                 {
                   {{1,1}, 0.0},
                   {{1,2}, 0.49},
                   {{2,1}, 0.01},
                   {{2,2}, 0.04}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- subset factorPtrs

        tstPtr = factorPtrs[3]->absorb(factorPtrs[4]);
        refPtr = uniqptr<DT>(
          new DT({2,3}, {doms2,doms3},  0.0,
                 {
                   {{1,1}, 0.0},
                   {{1,2}, 0.0},
                   {{2,1}, 0.05},
                   {{2,2}, 0.10}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- subset factorPtrs (flipped)

        tstPtr = factorPtrs[4]->absorb(factorPtrs[3]);
        refPtr = uniqptr<DT>(
          new DT({2,3}, {doms2,doms3},  0.0,
                 {
                   {{1,1}, 0.0},
                   {{1,2}, 0.0},
                   {{2,1}, 0.05},
                   {{2,2}, 0.10}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- partially overlapping factorPtrs

        tstPtr = factorPtrs[1]->absorb(factorPtrs[2]);
        if (*tstPtr !=  *factorPtrs[5]) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *factorPtrs[5] << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- partially overlapping factorPtrs (flipped)

        tstPtr = factorPtrs[2]->absorb(factorPtrs[1]);
        if (*tstPtr !=  *factorPtrs[5]) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *factorPtrs[5] << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- lhs and rhs are disjoint

        tstPtr = factorPtrs[6]->absorb(factorPtrs[2]);
        refPtr = uniqptr<DT>(
          new DT({2,3,5,6}, {doms2,doms3,doms5,doms6},  0.0,
                 {
                   {{1,1,1,1}, 0.00},
                   {{1,1,1,2}, 0.00},
                   {{1,1,2,1}, 0.00},
                   {{1,1,2,2}, 0.00},
                   {{1,2,1,1}, 0.35},
                   {{1,2,1,2}, 0.49},
                   {{1,2,2,1}, 0.07},
                   {{1,2,2,2}, 0.14},
                   {{2,1,1,1}, 0.05},
                   {{2,1,1,2}, 0.07},
                   {{2,1,2,1}, 0.01},
                   {{2,1,2,2}, 0.02},
                   {{2,2,1,1}, 0.1},
                   {{2,2,1,2}, 0.14},
                   {{2,2,2,1}, 0.02},
                   {{2,2,2,2}, 0.04}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== divide =======================

        //---------- fully overlapping non-empty factorPtrs

        tstPtr = factorPtrs[2]->absorb(factorPtrs[3]);
        tstPtr = tstPtr->cancel(factorPtrs[2]);
        if (*tstPtr !=  *factorPtrs[3]) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *factorPtrs[3] << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- partially overlapping non-empty factorPtrs

        tstPtr = factorPtrs[5]->cancel(factorPtrs[2]);
        //note, this is not factor 1, we stay on the larger scope
        refPtr = uniqptr<DT>(
          new DT({1,2,3}, {doms1,doms2,doms3},  0.0,
                 {
                   {{1,1,2}, 0.5},
                   {{1,2,1}, 0.8},
                   {{1,2,2}, 0.8},
                   {{2,1,2}, 0.1},
                   {{3,1,2}, 0.3},
                   {{3,2,1}, 0.9},
                   {{3,2,2}, 0.9}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== marginalize ==================

        //---------- marginalize to subset

        tstPtr = factorPtrs[2]->marginalize({2}, false, nullptr);
        refPtr = uniqptr<DT>(
          new DT({2}, {doms2},  0.0,
                 {
                   {{1}, 0.7},
                   {{2}, 0.3}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- marginalize to subset + extras

        tstPtr = factorPtrs[2]->marginalize({2,7}, false, nullptr);
        refPtr = uniqptr<DT>(
          new DT({2}, {doms2},  0.0,
                 {
                   {{1}, 0.7},
                   {{2}, 0.3}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== max marginalize ==============

        rcptr<FactorOperator> maxMargPtr = uniqptr<FactorOperator>( new DiscreteTable_MaxMarginalize<AssignType> );
        tstPtr = factorPtrs[5]->marginalize({2,3}, false, nullptr, maxMargPtr);
        refPtr = uniqptr<DT>(
          new DT({2,3}, {doms2,doms3},  0.0,
                 {
                   {{1,2}, 0.35},
                   {{2,1}, 0.09},
                   {{2,2}, 0.18}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable_MaxMarginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //========== observe and reduce ===========!!!

        //---------- observe empty set

        tstPtr = factorPtrs[5]->observeAndReduce({}, {}, false);

        if (*tstPtr !=  *factorPtrs[5]) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable_ObserveAndReduce FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *factorPtrs[5] << "\nGot:\n" << *tstPtr);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- observe subset

        tstPtr = factorPtrs[5]->observeAndReduce({1,3}, {AssignType(1),AssignType(2)} );
        refPtr = uniqptr<DT>(
          new DT({2}, {doms2},  0.0,
                 {
                   {{1}, 0.35},
                   {{2}, 0.16}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable_ObserveAndReduce FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- observe subset + extras

        tstPtr = factorPtrs[5]->observeAndReduce(
          {1,3,7},
          {AssignType(1),AssignType(2),AssignType(2)});
        refPtr = uniqptr<DT>(
          new DT({2}, {doms2},  0.0,
                 {
                   {{1}, 0.35},
                   {{2}, 0.16}
                 }, 0.0, 0.0) );

        if (*tstPtr !=  *refPtr) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable_ObserveAndReduce FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refPtr << "\nGot:\n" << *tstPtr);
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        //---------- observe disjoint

        tstPtr = factorPtrs[3]->observeAndReduce({7,8}, {AssignType(1),AssignType(2)});

        if (*tstPtr !=  *factorPtrs[3]) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *factorPtrs[3] << "\nGot:\n" << *tstPtr);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "LST FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
/*
 */

/*
  refPtr = uniqptr<DT>(
  new DT({}, {},  0.0,
  {
  }, 0.0, 0.0) );
*/

      } // test cases with defVal == 0

      // ###################################
      // ### Test cases with defVals > 0 ###
      // ###################################

      {

        /*
        // ****************
        // *** distance ***
        // ****************

        {

        vector<double> v1{0.0, 0.375, 0.125, 0.25, 0.25, 0};
        vector<double> v2{0, 0.125, 0.25, 0.25, 0.375, 0.0};
        double refDist =  DiscreteTable1SidedKLDist(v1, v2); // 0.20599

        rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
        dtPtr1 = uniqptr<DT>(
        new DT({0u}, {doms0},  0.25,
        { {{0}, 0.0}, {{1}, 0.375}, {{2}, 0.125}, {{5}, 0.0} }) );
        dtPtr2 = uniqptr<DT>(
        new DT({0u}, {doms0},  0.25,
        { {{0}, 0.0}, {{4}, 0.375}, {{1}, 0.125}, {{5}, 0.0} }) );
        //cout << *dtPtr1 << endl;
        //cout << *dtPtr2 << endl;

        double tstDist = dtPtr1->distance( dtPtr2.get() );
        if (std::abs(tstDist-refDist) > tol) {
        LOGERROR(pLog, "Testing", "DiscreteTable::distance FAILED!");
        LOGERROR(pLog, "Testing", "\nExpected distance = " << refDist
        << " but got " << tstDist);
        THROW(prlite::TestFailed, "");
        } // if

        }

        {
        vector<double> v1{0.0, 0.375, 0.125, 0.25, 0.25, 0};
        vector<double> v2{0, 0.125, 0.25, 0.25, 0.375, 0.0};
        double refDist =  DiscreteTable1SidedKLDist(v1, v2); // 0.20599

        rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
        dtPtr1 = uniqptr<DT>(
        new DT({0u}, {doms0},  0.25,
        { {{0}, 0.0}, {{1}, 0.375}, {{2}, 0.125}, {{5}, 0.0} }) );
        dtPtr2 = uniqptr<DT>(
        new DT({0u}, {doms0},  0.0,
        { {{4}, 0.375}, {{1}, 0.125}, {{3}, 0.25}, {{2}, 0.25} }) );
        //cout << *dtPtr1 << endl;
        //cout << *dtPtr2 << endl;

        double tstDist = dtPtr1->distance( dtPtr2.get() );
        if (std::abs(tstDist-refDist) > tol) {
        LOGERROR(pLog, "Testing", "DiscreteTable::distance FAILED!");
        LOGERROR(pLog, "Testing", "\nExpected distance = " << refDist
        << " but got " << tstDist);
        THROW(prlite::TestFailed, "");
        } // if
        }
        */

        // *****************
        // *** normalize ***
        // *****************

        { // +++ fairly plain, non-zero default +++

          rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
          refPtr = uniqptr<DT>(
            new DT({0u}, {doms0},  0.25,
                   { {{0}, 0.0}, {{1}, 0.375}, {{2}, 0.125}, {{5}, 0.0} }) );

          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms0},  2.0,
                   { {{0}, 0.0}, {{1}, 3.0}, {{2}, 1.0}, {{5}, 0.0} }) );

          tstPtr = dtPtr1->Factor::normalize();
          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ one with (inactive) flooring +++

          rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
          refPtr = uniqptr<DT>(
            new DT({0u}, {doms0},  0.05/1.1,
                   { {{1}, 0.375/1.1}, {{2}, 0.125/1.1}, {{3}, 0.25/1.1}, {{4}, 0.25/1.1} }, 1E-10, 0.1) );

          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms0},  0.05,
                   { {{1}, 0.375}, {{2}, 0.125}, {{3}, 0.25}, {{4}, 0.25} }, 1E-10, 0.1) );

          tstPtr = dtPtr1->Factor::normalize();
          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ one with (active) flooring +++

          rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
          refPtr = uniqptr<DT>(
            new DT({0u}, {doms0},  0.375*0.1/1.075,
                   { {{1}, 0.375/1.075}, {{2}, 0.125/1.075}, {{3}, 0.25/1.075}, {{4}, 0.25/1.075} }, 1E-10, 0.1) );

          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms0},  0.01,
                   { {{1}, 0.375}, {{2}, 0.125}, {{3}, 0.25}, {{4}, 0.25} }, 1E-10, 0.1) );

          tstPtr = dtPtr1->Factor::normalize();
          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ max normalize with flooring +++

          rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
          refPtr = uniqptr<DT>(
            new DT({0u}, {doms0},  0.05/0.375,
                   {
                     {{1}, 0.375/0.375},
                     {{2}, 0.125/0.375},
                     {{3}, 0.25/0.375},
                     {{4}, 0.25/0.375}
                   }, 1E-10, 0.1) );

          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms0},  0.05,
                   {
                     {{1}, 0.375},
                     {{2}, 0.125},
                     {{3}, 0.25},
                     {{4}, 0.25}
                   }, 1E-10, 0.1) );
          //cout << "First: " << *dtPtr1 << endl;

          rcptr<FactorOperator> normPtr = uniqptr<FactorOperator>(
            new DiscreteTable_MaxNormalize<AssignType> );
          tstPtr = dtPtr1->Factor::normalize(normPtr);
          //cout << "maxNormed: " << *tstPtr << endl;

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }

        //++++++++++ normalize an empty factor

        // an empty table is a killer in a pgm, it fractures it into islands.

        try {
          rcptr< vector<AssignType> > doms0(new vector<AssignType>{0,1,2,3,4,5});
          dtPtr1 = uniqptr<DT>(
            new DT({0u}, {doms0},  0.00) );
          tstPtr = dtPtr1->Factor::normalize();
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::normalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "normalizing an empty factor should fail!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // try
        catch (const char* s) {
          //if we go through here the exception worked correctly
        } // catch

        // *********************
        // *** inplaceAbsorb ***
        // *********************

        // ++++++++++++++++++++++++
        // +++ lhs contains rhs +++
        // ++++++++++++++++++++++++

        { // +++ lhs sparse probs, def val > 0, rhs no sparse probs  +++
          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2) );

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // another one
        { // --- lhs sparse probs, def val > 0, rhs no sparse probs  ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2 ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1->Factor::normalize()
                     << "Got " << *tstPtr->normalize());
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ lhs no sparse probs, rhs sparse probs, def val > 0  +++

          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2,
                                      {
                                        {{21}, 0.3},
                                        {{22}, 0.2}
                                      }) );

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);

          refPtr = uniqptr<DT>(new DT({1,2}, {doms1, doms2}, 0.0,
                                      {
                                        {{11,21}, 0.3},
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3},
                                        {{12,22}, 0.2}
                                      } ) );
          if (*tstPtr->normalize() != *refPtr->normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr->normalize()
                     << "Got " << *tstPtr->normalize());
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        { // --- lhs sparse probs, def val == 0, rhs sparse probs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{11,22,32}, 0.2},
                                        {{12,21,32}, 0.5},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,22,32}, 0.04},
                                        {{12,21,32}, 0.15},
                                        {{12,22,32}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, def val > 0, rhs sparseprobs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.02,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,21,32}, 0.03},
                                        {{12,21,31}, 0.04},
                                        {{12,21,32}, 0.03},
                                        {{12,22,32}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        { // Willie's test example
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1});

          dtPtr1 = uniqptr<DT>(new DT({5,6,7,9}, {doms2,doms2,doms2,doms2}, 1.0) );

          dtPtr2 = uniqptr<DT>(new DT({6,7,9}, {doms2,doms2,doms2}, 0.0,
                                      {
                                        {{0,0,0}, 1.0},
                                        {{1,0,1}, 1.0},
                                        {{1,1,0}, 1.0},
                                        {{1,1,1}, 1.0},
                                      } ) );

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);

          refPtr = uniqptr<DT>(new DT({5,6,7,9}, {doms2,doms2,doms2,doms2}, 0.0,
                                      {
                                        {{0,0,0,0}, 1.0},
                                        {{0,1,0,1}, 1.0},
                                        {{0,1,1,0}, 1.0},
                                        {{0,1,1,1}, 1.0},
                                        {{1,0,0,0}, 1.0},
                                        {{1,1,0,1}, 1.0},
                                        {{1,1,1,0}, 1.0},
                                        {{1,1,1,1}, 1.0},
                                      } ) );
          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        // **************
        // *** absorb ***
        // **************

        // ++++++++++++++++++++++++++++++++++++++
        // +++ lhs and rhs has identical vars +++
        // ++++++++++++++++++++++++++++++++++++++

        { // --- lhs and rhs both have sparse and def probs ---

          RVIds vars1 {1,2};
          RVIds vars2 {1,2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.7,
                                      {
                                        {{11,21}, 0.1},
                                        {{12,21}, 0.4}
                                      }) );


          tstPtr = dtPtr1->Factor::absorb(dtPtr2);

          refPtr = uniqptr<DT>(new DT({1,2}, {doms1, doms2}, 0.35,
                                      {
                                        {{11,21}, 0.05},
                                        {{11,22}, 0.14},
                                        {{12,21}, 0.12}
                                      } ) );

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs has only def prob while rhs has only sparse probs ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});

          dtPtr1 = uniqptr<DT>(new DT({111}, {doms1}, 0.5) );
          //cout << "Factor1: " << *dtPtr1 << endl;
          dtPtr2 = uniqptr<DT>(new DT({111}, {doms1}, 0.0, { {{0}, 1.0} }) );
          //cout << "Factor2: " << *dtPtr2 << endl;
          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Result: "<< *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr2->Factor::normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr2
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- rhs has only def prob while lhs has only sparse probs ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});

          dtPtr1 = uniqptr<DT>(new DT({111}, {doms1}, 0.5) );
          //cout << "Factor1: " << *dtPtr1 << endl;
          dtPtr2 = uniqptr<DT>(new DT({111}, {doms1}, 0.0, { {{0}, 1.0} }) );
          //cout << "Factor2: " << *dtPtr2 << endl;
          tstPtr = dtPtr2->Factor::absorb(dtPtr1);
          //cout << "Result: "<< *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr2->Factor::normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr2
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // ++++++++++++++++++++++++
        // +++ lhs contains rhs +++
        // ++++++++++++++++++++++++


        { // +++ lhs sparse probs, def val > 0, rhs no sparse probs  +++
          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        // another one
        { // --- lhs sparse probs, def val > 0, rhs no sparse probs  ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2 ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1->Factor::normalize()
                     << "Got " << *tstPtr->normalize());
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ lhs no sparse probs, rhs sparse probs, def val > 0  +++

          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2,
                                      {
                                        {{21}, 0.3},
                                        {{22}, 0.2}
                                      }) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);

          refPtr = uniqptr<DT>(new DT({1,2}, {doms1, doms2}, 0.0,
                                      {
                                        {{11,21}, 0.3},
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3},
                                        {{12,22}, 0.2}
                                      } ) );
          if (*tstPtr->normalize() != *refPtr->normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr->normalize()
                     << "Got " << *tstPtr->normalize());
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        { // --- lhs sparse probs, def val == 0, rhs sparse probs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{11,22,32}, 0.2},
                                        {{12,21,32}, 0.5},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,22,32}, 0.04},
                                        {{12,21,32}, 0.15},
                                        {{12,22,32}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, def val > 0, rhs sparseprobs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.02,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,21,32}, 0.03},
                                        {{12,21,31}, 0.04},
                                        {{12,21,32}, 0.03},
                                        {{12,22,32}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // ++++++++++++++++++++++++
        // +++ rhs contains lhs +++
        // ++++++++++++++++++++++++

        { // --- lhs no sparse values, rhs with sparse values ---

          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2) );

          tstPtr = dtPtr2->Factor::absorb(dtPtr1);

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize() ) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, def val > 0, rhs sparseprobs, def val > 0 ---

          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr2->Factor::absorb(dtPtr1);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.02,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,21,32}, 0.03},
                                        {{12,21,31}, 0.04},
                                        {{12,21,32}, 0.03},
                                        {{12,22,32}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }


        // +++++++++++++++++++++++
        // +++ partial overlap +++
        // +++++++++++++++++++++++

        { // --- lhs no sparse probs, rhs no sparse probs ---
          RVIds vars1 {1,2};
          RVIds vars2 {2,3};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.9) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2,doms3}, 0.2) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2)->normalize();

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1,doms2,doms3} ) )->Factor::normalize();

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, rhs no sparse probs ---
          RVIds vars1 {1,2};
          RVIds vars2 {2,3};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2,doms3}, 0.2) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2)->Factor::normalize();

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.5,
                                      {
                                        {{11,22,31}, 0.2},
                                        {{11,22,32}, 0.2},
                                        {{11,22,33}, 0.2},
                                        {{11,22,34}, 0.2},
                                        {{12,21,31}, 0.3},
                                        {{12,21,32}, 0.3},
                                        {{12,21,33}, 0.3},
                                        {{12,21,34}, 0.3}
                                      }) )->Factor::normalize();

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // another one
        { // --- lhs has sparse probs, rhs no sparse probs ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2}, {doms1, doms2}, 0.1,
                                      {
                                        {{11,21}, 0.6},
                                        {{12,22}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2 ));
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{11,21,32}, 0.6},
                                        {{11,21,33}, 0.6},
                                        {{11,21,34}, 0.6},
                                        {{12,22,31}, 0.3},
                                        {{12,22,32}, 0.3},
                                        {{12,22,33}, 0.3},
                                        {{12,22,34}, 0.3}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr->normalize()
                     << "Got " << *tstPtr->normalize() );
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // another one
        { // --- lhs sparse probs, rhs no sparse probs ---
          RVIds vars1 {1,2};
          RVIds vars2 {2,3};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2,doms3}, 0.2) );

          tstPtr = dtPtr2->Factor::absorb(dtPtr1)->Factor::normalize();

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.5,
                                      {
                                        {{11,22,31}, 0.2},
                                        {{11,22,32}, 0.2},
                                        {{11,22,33}, 0.2},
                                        {{11,22,34}, 0.2},
                                        {{12,21,31}, 0.3},
                                        {{12,21,32}, 0.3},
                                        {{12,21,33}, 0.3},
                                        {{12,21,34}, 0.3}
                                      }) )->Factor::normalize();

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, rhs sparse probs ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2}, {doms1, doms2}, 0.1,
                                      {
                                        {{11,21}, 0.6},
                                        {{12,22}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;

          refPtr = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.02,
                                      {
                                        {{11,21,31}, 0.24},
                                        {{11,21,32}, 0.18},
                                        {{11,21,33}, 0.12},
                                        {{11,21,34}, 0.12},
                                        {{12,21,31}, 0.04},
                                        {{12,21,32}, 0.03},
                                        {{12,22,31}, 0.06},
                                        {{12,22,32}, 0.06},
                                        {{12,22,33}, 0.06},
                                        {{12,22,34}, 0.06}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          if (*tstPtr->normalize() != *refPtr->normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::absorb FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        // *********************
        // *** inplaceCancel ***
        // *********************

        // ++++++++++++++++++++++++
        // +++ lhs contains rhs +++
        // ++++++++++++++++++++++++

        { // Willie's test example
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1});

          dtPtr1 = uniqptr<DT>(new DT({5,6,7,9}, {doms2,doms2,doms2,doms2}, 1.0) );

          dtPtr2 = uniqptr<DT>(new DT({6,7,9}, {doms2,doms2,doms2}, 0.0,
                                      {
                                        {{0,0,0}, 1.0},
                                        {{1,0,1}, 1.0},
                                        {{1,1,0}, 1.0},
                                        {{1,1,1}, 1.0},
                                      } ) );

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          tstPtr->Factor::inplaceAbsorb(dtPtr2);
          tstPtr->Factor::inplaceCancel(dtPtr2);

          // Careful now, the zeros in dtPtr2 will remain in the
          // result. I.e. we can not recover dtPtr1 fully.
          refPtr = uniqptr<DT>(new DT({5,6,7,9}, {doms2,doms2,doms2,doms2}, 0.0,
                                      {
                                        {{0,0,0,0}, 1.0},
                                        {{0,1,0,1}, 1.0},
                                        {{0,1,1,0}, 1.0},
                                        {{0,1,1,1}, 1.0},
                                        {{1,0,0,0}, 1.0},
                                        {{1,1,0,1}, 1.0},
                                        {{1,1,1,0}, 1.0},
                                        {{1,1,1,1}, 1.0},
                                      } ) );
          if (*tstPtr->normalize() != *refPtr->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceCancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        // **************
        // *** cancel ***
        // **************


        // ++++++++++++++++++++++++++++++++++++++
        // +++ lhs and rhs has identical vars +++
        // ++++++++++++++++++++++++++++++++++++++

        { // --- lhs and rhs both have sparse and def probs ---

          RVIds vars1 {1,2};
          RVIds vars2 {1,2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.7,
                                      {
                                        {{11,21}, 0.1},
                                        {{12,21}, 0.4}
                                      }) );


          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // ++++++++++++++++++++++++
        // +++ lhs contains rhs +++
        // ++++++++++++++++++++++++


        { // +++ lhs sparse probs, def val > 0, rhs no sparse probs  +++
          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5,
                                      {
                                        {{11,22}, 0.2},
                                        {{12,21}, 0.3}
                                      }) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // another one
        { // --- lhs sparse probs, def val > 0, rhs no sparse probs  ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2 ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // +++ lhs no sparse probs, rhs sparse probs, def val > 0  +++

          RVIds vars1 {1,2};
          RVIds vars2 {2};
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});

          dtPtr1 = uniqptr<DT>(new DT(vars1, {doms1, doms2}, 0.5) );

          dtPtr2 = uniqptr<DT>(new DT(vars2, {doms2}, 0.2,
                                      {
                                        {{21}, 0.3},
                                        {{22}, 0.2}
                                      }) );

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }


        { // --- lhs sparse probs, def val == 0, rhs sparse probs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{11,22,32}, 0.2},
                                        {{12,21,32}, 0.5},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, rhs sparse probs, both with def val ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{11,21,32}, 0.6},
                                        {{11,21,33}, 0.6},
                                        {{11,21,34}, 0.6},
                                        {{12,22,31}, 0.3},
                                        {{12,22,32}, 0.3},
                                        {{12,22,33}, 0.3},
                                        {{12,22,34}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          //cout << "Prod: " << *tstPtr << endl;
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        { // --- lhs sparse probs, def val > 0, rhs sparseprobs, def val > 0 ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{11,12});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{21,22});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{31,32,33,34});

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.1,
                                      {
                                        {{11,21,31}, 0.6},
                                        {{12,22,32}, 0.3}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          dtPtr2 = uniqptr<DT>(new DT({2,3}, {doms2,doms3}, 0.2,
                                      {
                                        {{21,31}, 0.4},
                                        {{21,32}, 0.3}
                                      } ) );
          //cout << "Second: " << *dtPtr2 << endl;

          tstPtr = dtPtr1->Factor::absorb(dtPtr2);
          tstPtr = tstPtr->cancel(dtPtr2);
          //cout << "Div: " << *tstPtr << endl;

          if (*tstPtr->normalize() != *dtPtr1->Factor::normalize()) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::cancel FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *dtPtr1
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if
        }

        // *******************
        // *** marginalize ***
        // *******************

        { // --- non-zero default value ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1,2});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{0,1});

          refPtr = uniqptr<DT>(new DT({1}, {doms1}, 42.0,
                                      {
                                        {{0}, 19.0},
                                        {{1}, 37.0}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 7.0,
                                      {
                                        {{0,0,0}, 0.0},
                                        {{0,0,1}, 1.0},
                                        {{0,1,0}, 2.0},
                                        {{0,2,0}, 4.0},
                                        {{0,2,1}, 5.0},
                                        {{1,0,0}, 6.0},
                                        {{1,2,0}, 3.0}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          tstPtr = dtPtr1->Factor::marginalize({1u});
          //cout << "Marg: " << *tstPtr << std::endl;

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }

        { // --- max marginalize, non-zero default value ---

          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1,2});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{0,1});

          rcptr<FactorOperator> margPtr = uniqptr<FactorOperator>(
            new DiscreteTable_MaxMarginalize<AssignType> );

          refPtr = uniqptr<DT>(new DT({1,2}, {doms1,doms2}, 7.0,
                                      {
                                        {{0,0}, 1.0},
                                        {{0,1}, 7.0},
                                        {{0,2}, 5.0},
                                        {{1,0}, 7.0},
                                        {{1,2}, 7.0}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 7.0,
                                      {
                                        {{0,0,0}, 0.0},
                                        {{0,0,1}, 1.0},
                                        {{0,1,0}, 2.0},
                                        {{0,2,0}, 4.0},
                                        {{0,2,1}, 5.0},
                                        {{1,0,0}, 6.0},
                                        {{1,2,0}, 3.0}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          tstPtr = dtPtr1->Factor::marginalize({1u,2u}, true, nullptr, margPtr);
          //cout << "Marg: " << *tstPtr << std::endl;

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }

        { // --- max marginalize, non-zero default value ---

          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1,2});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{0,1});

          refPtr = uniqptr<DT>(new DT({1,2}, {doms1,doms2}, 0.5,
                                      {
                                        {{0,0}, 1.0},
                                        {{0,1}, 2.0},
                                        {{0,2}, 5.0},
                                        {{1,0}, 6.0},
                                        {{1,2}, 3.0}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.5,
                                      {
                                        {{0,0,0}, 0.0},
                                        {{0,0,1}, 1.0},
                                        {{0,1,0}, 2.0},
                                        {{0,2,0}, 4.0},
                                        {{0,2,1}, 5.0},
                                        {{1,0,0}, 6.0},
                                        {{1,2,0}, 3.0}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          rcptr<FactorOperator> margPtr = uniqptr<FactorOperator>(
            new DiscreteTable_MaxMarginalize<AssignType> );
          tstPtr = dtPtr1->Factor::marginalize({1u,2u}, true, nullptr, margPtr);
          //cout << "Marg: " << *tstPtr << std::endl;

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalize FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }
        // ************************
        // *** observeAndReduce ***
        // ************************

        { // --- zero default value ---
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1});
          rcptr< vector<AssignType> > doms2 (
            new vector<AssignType>{0,1,2});
          rcptr< vector<AssignType> > doms3 (
            new vector<AssignType>{0,1});

          refPtr = uniqptr<DT>(new DT({3}, {doms3}, 0.0,
                                      {
                                        {{0}, 4.0},
                                        {{1}, 5.0}
                                      } ) );
          //cout << "Ref: " << *refPtr << endl;

          dtPtr1 = uniqptr<DT>(new DT({1,2,3}, {doms1, doms2, doms3}, 0.0,
                                      {
                                        {{0,0,0}, 0.0},
                                        {{0,0,1}, 1.0},
                                        {{0,1,0}, 2.0},
                                        {{0,2,0}, 4.0},
                                        {{0,2,1}, 5.0},
                                        {{1,0,0}, 6.0},
                                        {{1,2,0}, 3.0}
                                      } ) );
          //cout << "First: " << *dtPtr1 << endl;

          tstPtr = dtPtr1->Factor::observeAndReduce({1,2},{0u,2u});
          //cout << "Reduced: " << *tstPtr << std::endl;

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::observeAndReduce FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }


        // *********************
        // *** inplaceDampen ***
        // *********************

        // 0.9*v1+0.1*v2 = {0, 0.35,  0.1375, 0.25, 0.2625 0};
        {

          vector<double> v1{0.0, 0.375, 0.125, 0.25, 0.25, 0};
          vector<double> v2{0, 0.125, 0.25, 0.25, 0.375, 0.0};
          double refDist =  DiscreteTable2SidedKLDist(v1, v2); // 0.20599

          rcptr< vector<AssignType> > doms0(
            new vector<AssignType>{0,1});
          rcptr< vector<AssignType> > doms1 (
            new vector<AssignType>{0,1,2});

          refPtr = uniqptr<DT>(
            new DT({0,1}, {doms0,doms1},  0.0,
                   {
                     {{0,1}, 0.35},
                     {{0,2}, 0.1375},
                     {{1,0}, 0.25},
                     {{1,1}, 0.2625}
                   }) );

          dtPtr1 = uniqptr<DT>(
            new DT({0,1}, {doms0,doms1},  0.25,
                   {
                     {{0,0}, 0.0},
                     {{0,1}, 0.375},
                     {{0,2}, 0.125},
                     {{1,2}, 0.0}
                   }) );
          dtPtr2 = uniqptr<DT>(
            new DT({0,1}, {doms0,doms1},  0.0,
                   {
                     {{0,1}, 0.125},
                     {{0,2}, 0.25},
                     {{1,0}, 0.25},
                     {{1,1}, 0.375}
                   }) );

          tstPtr = uniqptr<Factor>( dtPtr1->copy() );
          double df = 0.1;
          double tstDist = tstPtr->Factor::inplaceDampen(dtPtr2, df);
          tstDTPtr = std::dynamic_pointer_cast< DiscreteTable<AssignType> > (tstPtr);

          if (std::abs(tstDist-refDist) > tol) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceDampen distance FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected distance = " << refDist
                     << " but got " << tstDist);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

          if (*tstPtr != *refPtr) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceDampen FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected " << *refPtr
                     << "Got " << *tstPtr);
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        }

      } // test cases with defVal != 0


      // ################
      // ### Cat1EXPF ###
      // ################

      {
        rcptr< vector<AssignType> > domsBin(
          new vector<AssignType>{0,1});
        rcptr< vector<AssignType> > domsTri (
          new vector<AssignType>{0,1,2});

        dtPtr1 = uniqptr<DT>(
          new DT({0,1,2}, {domsBin,domsTri, domsBin},  0.12,
                 {
                   {{0,0,0}, 0.04},
                   {{0,0,1}, 0.04},
                   {{1,0,0}, 0.04},
                   {{1,0,1}, 0.04},
                   {{0,1,1}, 0.08},
                   {{1,1,0}, 0.08},
                   {{1,1,1}, 0.08}
                 }) );
        dtPtr1->Factor::inplaceNormalize();
        //cout << *dtPtr1 << endl;

        // *****************************
        // *** marginalizeToCat1EXPF ***
        // *****************************

        rcptr< Cat1EXPF<AssignType> > cat1Ptr = uniqptr< Cat1EXPF<AssignType> >( dtPtr1->marginalizeToCat1EXPF(1) );
        //cout << *cat1Ptr << endl;


        if (std::abs(cat1Ptr->getProb(0)-0.16) > tol) { // 4/25
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalizeToCat1EXPF FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected p(x=0)  = 0.16, but got " << cat1Ptr->getProb(0));
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        if (std::abs(cat1Ptr->getProb(1)-0.36) > tol) { // 9/25
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalizeToCat1EXPF FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected p(x=1)  = 0.36, but got " << cat1Ptr->getProb(1));
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        if (std::abs(cat1Ptr->getProb(2)-0.48) > tol) { // 12/25
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::marginalizeToCat1EXPF FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected p(x=2)  = 0.48, but got " << cat1Ptr->getProb(2));
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // *****************************
        // *** inplaceAbsorbCat1EXPF ***
        // *****************************

        dtPtr1 = uniqptr<DT>(
          new DT({0,1,2}, {domsBin,domsTri, domsBin},  3,
                 {
                   {{0,0,0}, 1},
                   {{0,0,1}, 1},
                   {{1,0,0}, 1},
                   {{1,0,1}, 1},
                   {{0,1,0}, 2},
                   {{0,1,1}, 2},
                   {{1,1,0}, 2},
                   {{1,1,1}, 2}
                 }) );
        dtPtr1->Factor::inplaceNormalize();
        //cout << *dtPtr1 << endl;

        cat1Ptr = uniqptr< Cat1EXPF<AssignType> >(
          new Cat1EXPF<AssignType>(1, domsTri,  8,
                 {
                   {0, 24},
                   {1, 12}
                 }) );
        cat1Ptr->inplaceNormalize();
        //cout << *cat1Ptr << endl;
        dtPtr1->inplaceAbsorbCat1EXPF(cat1Ptr);
        dtPtr1->Factor::inplaceNormalize();
        //cout << *dtPtr1 << endl;

        if (std::abs(dtPtr1->getDefaultProb()-1.0/12) > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable::inplaceAbsorbCat1EXPF FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected default prob 0.0833, but got " << dtPtr1->getDefaultProb());
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // Cat1EXPF

      // #########################
      // ### integration tests ###
      // #########################

#ifdef TESTING_MODE
    FILEPOS;
#endif

    } // try

    catch (const char* msg) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable threw an exception!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", msg);
      PRLITE_THROW(prlite::TestFailed, "");
    } // catch

  catch (const exception& e) {
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", "DiscreteTable threw an exception!");
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", e.what() );
    throw(e);
  } // catch

  catch(...) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DiscreteTable threw an exception!");
    PRLITE_THROW(prlite::TestFailed, "");
  } // catch


  } // runTests

}; // TestDiscreteTable

PRLITE_REGISTER_TESTCASE(TestDiscreteTable);


/*
      // ****************
      // *** distance ***
      // ****************

      {

        // lets first find the ref distance between test vectors
        vector<double> v1{0.05, 0.375, 0.125, 0.25, 0.2, 0};
        vector<double> v2{0, 0.125, 0.25, 0.2, 0.375, 0.05};
        double refDist =  DiscreteTable1SidedKLDist(v1, v2); // 0.20599
        //cout << "refDist: " << refDist << endl;

        rcptr< vector<AssignType> > doms00(
          new vector<AssignType>{0,1,2,3,4,5});
        dtPtr1 = uniqptr<DT>(
          new DT({0u}, {doms00}, 0.0,
                 {
                   {{0}, 0.05},
                   {{1}, 0.375},
                   {{2}, 0.125},
                   {{3}, 0.25},
                   {{4}, 0.2}
                 }) );
        dtPtr2 = uniqptr<DT>(
          new DT({0u}, {doms00},  0.0,
                 {
                   {{1}, 0.125},
                   {{2}, 0.25},
                   {{3}, 0.2},
                   {{4}, 0.375},
                   {{5}, 0.05},
                 }) );
        //cout << *dtPtr1 << endl;
        //cout << *dtPtr2 << endl;

        double tstDist = dtPtr1->distance( dtPtr2.get() );

        if (std::abs(tstDist-refDist) > tol) {
          LOGERROR(pLog, "Testing", "DiscreteTable::distance FAILED!");
          LOGERROR(pLog, "Testing", "\nExpected distance = " << refDist
                   << " but got " << tstDist);
          THROW(prlite::TestFailed, "");
        } // if

      }

 */
