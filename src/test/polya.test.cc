/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2018
          PURPOSE: Test code for polya routines.

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

/*
Copyright (c) 1995-2006, Stellenbosch University
All rights reserved.

*/

//#ifndef NO_TESTING

// patrec headers
#include "prlite_error.hpp"
#include "oddsandsods.hpp"
#include "polya.hpp"
#include "discretetable.hpp"
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "lbu_cg.hpp"
#include "lbp_cg.hpp"
#include "prlite_testing.hpp"  // TestCase
//#include <boost/any.hpp>

// emdw headers

//#include <string>  // string

using namespace std;
using namespace emdw;

class BaseTestPolya : public prlite::TestCase {
public:

  /** setup to be run before each test. */
  virtual void setup() {} // setup


  /** clean up after. */
  virtual void tearDown() {} // tearDown

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

}; // BaseTestPolya

class TestPolya : public BaseTestPolya {
public:

  /** a name. */
  const string& isA() {
    static const string name = "TestPolya";
    return name;
  } // isA

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

    try{

#ifdef TESTING_MODE
    FILEPOS;
#endif


    typedef unsigned AssignType;
    typedef Dirichlet<AssignType> DIR;
    typedef Polya<AssignType> POL;
    double theDef = 1.0;

      rcptr< vector<AssignType> > binDom( new vector<AssignType>({0,1}) );
      rcptr< vector<AssignType> > cat3Dom( new vector<AssignType>({0,1,2}) );

      // #####################################
      // ### define some generic variables ###
      // #####################################

      // Same as dir1_DIR, but constructed differently
      rcptr<DIR> dir0_DIR = uniqptr<DIR> (new DIR(0u, binDom, 4.0, { {1, 3.0} } ) );
      dir0_DIR->keepSeparate = 1;
      //std::cout << "Dirichlet0: " << *dir0_DIR << std::endl;

      // Same as cat1_DT, but constructed differently
      map<vector<AssignType>, FProb> table0;
      // note not normalized
      table0[{1}] = 0.8;
      rcptr< DiscreteTable<AssignType> > cat0_DT = uniqptr< DiscreteTable<AssignType> > (
        new DiscreteTable<AssignType>(
          {1u}, {binDom},
          1.2, table0,
          0.0, 0.0, true) );
      //std::cout << "Categorical0: " << *cat0_DT << std::endl;

      // The first Dirichlet
      // alphas = [4,3];
      rcptr<DIR> dir1_DIR = uniqptr<DIR> (new DIR(0u, binDom, 0.0, { {0, 4.0}, {1, 3.0} } ) );
      dir1_DIR->keepSeparate = 1;
      //std::cout << "Dirichlet1: " << *dir1_DIR << std::endl;

      // The first DiscreteTable
      map<vector<AssignType>, FProb> table1;
      // note not normalized
      table1[{0u}] = 1.2;
      table1[{1}] = 2.0 - table1[{0u}].prob;
      rcptr< DiscreteTable<AssignType> > cat1_DT = uniqptr< DiscreteTable<AssignType> > (
        new DiscreteTable<AssignType>(
          {1u}, {binDom},
          0.0, table1,
          0.0, 0.0, true) );
      //std::cout << "Categorical1: " << *cat1_DT << std::endl;

      // The second Dirichlet
      map<AssignType, double> theDirPar2;
      // alphas = 2
      theDirPar2[0] = 2.0; theDirPar2[1] = 3.0;
      rcptr<DIR> dir2_DIR = uniqptr<DIR> (new DIR(0u, binDom, theDef, theDirPar2) );
      //std::cout << "Dirichlet2: " << *dir2_DIR << std::endl;

      // The second DiscreteTable
      map<vector<AssignType>, FProb> table2;
      // note not normalized
      table2[{0u}] = 0.25;
      table2[{1}] = 1.0 - table2[{0u}].prob;
      rcptr< DiscreteTable<AssignType> > cat2_DT = uniqptr< DiscreteTable<AssignType> > (
        new DiscreteTable<AssignType>(
          {1u}, {binDom},
          0.0, table2,
          0.0, 0.0, true) );
      //std::cout << "Categorical2: " << *cat2_DT << std::endl;

      // The Polya
      rcptr< Polya<AssignType> > pot1_POL, pot2_POL, tst_POL;
      rcptr<Factor> tst_FAC, tst2_FAC, tmp_FAC;
      rcptr< DiscreteTable<AssignType> > tst_DT;
      rcptr<DIR> tst_DIR;

      // ############################################
      // ### lets first test the factor operators ###
      // ############################################

      // remember that we have 4 configurations for each Polya,
      // depending on which of the the Dirichlet and DiscreteTable is
      // instantiated.


      // +++++++++++++++++
      // +++ normalize +++
      // +++++++++++++++++

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, cat1_DT) );

      tst_FAC = pot1_POL->Factor::normalize();

      // check that the categorical is normalized
      tst_POL = dynamic_pointer_cast< Polya<AssignType> > (tst_FAC);
      tst_DT = tst_POL->getDiscreteTable();
      try {
        assertAlmostEqual(1.0, tst_DT->sumZ(), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::normalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected unity sum on categorical probs, got " << tst_DT->sumZ() );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // the Dirichlet is implicitly normalized and should not change
      tst_DIR = tst_POL->getDirichlet();
      if (*tst_DIR !=  *dir1_DIR) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::normalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *dir1_DIR << endl << "Got:\n" << *tst_DIR);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if


      // +++++++++++++++++++++
      // +++ absorb/cancel +++
      // +++++++++++++++++++++

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, cat1_DT) );
      pot1_POL->inplaceNormalize(nullptr);

      // ------------------------------
      // --- absorb a DiscreteTable ---
      // ------------------------------

            tmp_FAC = uniqptr<POL>( pot1_POL->copy() );
      tmp_FAC->inplaceAbsorb(cat2_DT);

      tst_FAC = tmp_FAC->normalize();

      // check that the dirichlet alphas are correct
      tst_POL = dynamic_pointer_cast< Polya<AssignType> > (tst_FAC);
      tst_DIR = tst_POL->getDirichlet();
      if (tst_DIR->getAlpha({0}) !=  4 or tst_DIR->getAlpha({1}) != 3) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::absorb Categorical FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha counts [4,3], but got [" << tst_DIR->getAlpha({0}) << "," << tst_DIR->getAlpha({1})<< "]");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // check that the categorical probs are correct
      tst_DT = tst_POL->getDiscreteTable();
      try{
        assertAlmostEqual(1.0/3, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::absorb Categorical FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected P(0) = 0.333333 but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ------------------------------
      // --- cancel a DiscreteTable ---
      // ------------------------------

      tmp_FAC = uniqptr<Factor>( tst_FAC->copy() );
      tmp_FAC->inplaceCancel(cat2_DT);

      tst_FAC = tmp_FAC->normalize();

      if (*tst_FAC != *pot1_POL) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::cancel Categorical FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *pot1_POL << "\nbut got:\n" << *tst_FAC);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --------------------------
      // --- absorb a Dirichlet ---
      // --------------------------

            tmp_FAC = uniqptr<POL>( pot1_POL->copy() );
      tmp_FAC->inplaceAbsorb(dir2_DIR);


      tst_FAC = tmp_FAC->normalize();

      // check that the dirichlet alphas are correct
      tst_POL = dynamic_pointer_cast< Polya<AssignType> > (tst_FAC);
      tst_DIR = tst_POL->getDirichlet();
      if (tst_DIR->getAlpha({0}) !=  5 or tst_DIR->getAlpha({1}) != 5) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::absorb Dirichlet FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha [5, 5], but got [" << tst_DIR->getAlpha({0}) << ", " << tst_DIR->getAlpha({1}) << "]");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --------------------------
      // --- cancel a Dirichlet ---
      // --------------------------

            tmp_FAC = uniqptr<Factor>( tst_FAC->copy() );
      tmp_FAC->inplaceCancel(dir2_DIR);

      tst_FAC = tmp_FAC->normalize();

      if (*tst_FAC != *pot1_POL) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::cancel Dirichlet FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected:\n" << *pot1_POL << "\nbut got:\n" << *tst_FAC);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // +++++++++++++++++++
      // +++ marginalize +++
      // +++++++++++++++++++

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, cat1_DT) );

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if ( tst_DIR->getAlpha({0}) !=  (4.0+2.0/3) or tst_DIR->getAlpha({1}) != (3.0+1.0/3) ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 4.67 and 3.33 alpha counts, but got "
                 << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u})->normalize();

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(2.0/3, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 2.0/3 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // --- only the Dirichlet is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, nullptr) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(4.0/7, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 4.0/7 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (*tst_DIR != *dir1_DIR) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << *dir1_DIR << " but got " << *tst_FAC);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- Only the DiscreteTable is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, 0, cat1_DT) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 0.6 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  1.6 or tst_DIR->getAlpha({1}) != 1.4) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 1.6 and 1.4 alpha counts, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- only the Dirichlet is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, 0) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(4.0/7, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 4.0/7 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  4.0 or tst_DIR->getAlpha({1}) != 3.0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected [4.0, 3.0], but got ["
                 << tst_DIR->getAlpha({0}) << ", " << tst_DIR->getAlpha({1}) << "]");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- Neither the Dirichlet nor the the DiscreteTable is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, 0, 0) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(0.5, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::maginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 0.5 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  1.0 or tst_DIR->getAlpha({1}) != 1.0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 1.0 and 1.0 alpha counts, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- A TEST EXACTLY LIKE THE PREVIOUS SETE, BUT THE ---
      // --- DIRICHLET SAVES ONE ALPHA IN ITS DEFAULT       ---

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir0_DIR, cat0_DT) );

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if ( tst_DIR->getAlpha({0}) !=  (4.0+2.0/3) or tst_DIR->getAlpha({1}) != (3.0+1.0/3) ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 4.67 and 3.33 alpha counts, but got "
                 << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u})->normalize();

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(2.0/3, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 2.0/3 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // --- only the Dirichlet is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir0_DIR, nullptr) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(4.0/7, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 4.0/7 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (*tst_DIR != *dir0_DIR) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected " << *dir0_DIR << " but got " << *tst_FAC);
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- Only the DiscreteTable is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, 0, cat0_DT) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 0.6 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  1.6 or tst_DIR->getAlpha({1}) != 1.4) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 1.6 and 1.4 alpha counts, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- only the Dirichlet is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir0_DIR, 0) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(4.0/7, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 4.0/7 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  4.0 or tst_DIR->getAlpha({1}) != 3.0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected [4.0, 3.0], but got ["
                 << tst_DIR->getAlpha({0}) << ", " << tst_DIR->getAlpha({1}) << "]");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- Neither the Dirichlet nor the the DiscreteTable is instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, 0, 0) );

      // ... Keep the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::marginalize({1u});

      // check that the categorical probs are correct
      tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
      try{
        assertAlmostEqual(0.5, tst_DT->potentialAt({1u}, {0u}), 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::maginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = "
                 << 0.5 << " but got " << tst_DT->potentialAt({1u}, {0u}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

      // ... Keep the Dirichlet ...
      tst_FAC = pot1_POL->Factor::marginalize({0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  1.0 or tst_DIR->getAlpha({1}) != 1.0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 1.0 and 1.0 alpha counts, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if


      // ++++++++++++++++++++++++
      // +++ observeAndReduce +++
      // ++++++++++++++++++++++++

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, cat1_DT) );

      // ... We only support observing the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::observeAndReduce({1u},{0u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  5 or tst_DIR->getAlpha({1}) != 3) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::observeAndReduce FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 5 and 3 alpha counts, but got "
                 << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ... We only support observing the DiscreteTable ...
      tst_FAC = pot1_POL->Factor::observeAndReduce({1u},{1u});

      // check that the dirichlet alphas are correct
      tst_DIR = dynamic_pointer_cast<DIR> (tst_FAC);
      if (tst_DIR->getAlpha({0}) !=  4 or tst_DIR->getAlpha({1}) != 4) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::observeAndReduce FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 4 and 4 alpha counts, but got "
                 << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // +++++++++++++++++++++
      // +++ inplaceDampen +++
      // +++++++++++++++++++++

      // --- Both Dirichlet and DiscreteTable are instantiated ---
      pot1_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir1_DIR, cat1_DT) );
      pot2_POL =
        uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom, dir2_DIR, cat2_DT) );

      tst_FAC = rcptr<Factor>( pot1_POL->copy() );
      try{
        tst_FAC->inplaceDampen(pot2_POL, 0.3);

        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::inplaceDampen test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "unexpected call to non-existing inplaceDampen operation");
        PRLITE_THROW(prlite::TestFailed, "");
      } // try

      catch(const exception& e) {
        // if we go through here everything is fine
      } // catch


      // #################################
      // ### Now for integration tests ###
      // #################################


      vector< rcptr<Factor> > factors;

      // +++++++++++++++++++++++++++++++++++++
      // +++ Integration test 1:           +++
      // +++ A coin with some observations +++
      // +++++++++++++++++++++++++++++++++++++

      { // integration test 1 scope

        // the X for this Polya will not be observed. To support the
        // required multiplication this factor must also be the first
        // in our list.
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom) );
        factors.push_back(pot1_POL);
        //std::cout << factors.size()-1 << ": " << *factors[factors.size()-1] << std::endl;

        // the X's of all the next POT's are observed thereby reducing
        // them to Dirichlets.
        unsigned xId = 2u;
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, xId, binDom) );
        factors.push_back(
          pot1_POL->Factor::observeAndReduce(emdw::RVIds{xId},emdw::RVVals{0u}) );

        xId = 3u;
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, xId, binDom) );
        factors.push_back(
          pot1_POL->Factor::observeAndReduce(emdw::RVIds{xId},emdw::RVVals{0u}) );

        // !!! lyk my hierdie werk nie reg nie. die observation gaan
        // na kategorie 0 ipv kategorie 1
        xId = 4u;
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, xId, binDom) );
        factors.push_back(
          pot1_POL->Factor::observeAndReduce(emdw::RVIds{xId},emdw::RVVals{1u}) );

        // and add the Dirichlet prior to the end of the list
        factors.push_back(dir1_DIR);

        // First one is Polya, the rest are Dirichlet
                tst_FAC = uniqptr<Factor>( factors[0]->copy() );
        for (size_t k = 1; k < factors.size(); k++) {
          tst_FAC->inplaceAbsorb(factors[k]);
        } // for

        tst_POL = dynamic_pointer_cast< Polya<AssignType> > (tst_FAC);
        //std::cout << "joint: " << *tst_POL << std::endl;

        // check that the dirichlet alphas are correct
        tst_DIR = tst_POL->getDirichlet();
        if (tst_DIR->getAlpha({0}) !=  6 or tst_DIR->getAlpha({1}) != 4) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 1 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 6 and 4 alpha counts respectively, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check that the categorical distribution is absent
        tst_DT = tst_POL->getDiscreteTable();
        if (tst_DT) { // should be NULL
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 1 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable to be absent, but found it to not be so");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // this tests marginalization when there is not DiscreteTable
        // present in the Polya (marked as case 5)
        tst_FAC = tst_FAC->marginalize({1u});
        //std::cout << "marginal cat: " << *tst_FAC << std::endl;

        tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);

        // check that the categorical probs are correct
        try{
          assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
        } // try
        catch(const exception& e) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 1 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = 0.6 but got " << tst_DT->potentialAt({1u}, {0u}) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // catch

      } // integration test 1 scope

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // +++ Integration test 2:                               +++
      // +++ manual loopy belief on a coin indirectly observed +++
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if (false) { // integration test 2 scope

        factors.clear();
        // factor 0 contains the Dirichlet prior, as well as the unobserved X1
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0u, 1u, binDom) );
        factors.push_back(pot1_POL);
        factors[factors.size()-1]->inplaceAbsorb(dir1_DIR);

        //std::cout << factors.size()-1 << ": " << *factors[factors.size()-1] << std::endl;

        // Factor 1 is for X2 which will be indirectly observed as being 0
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0, 2u, binDom) );
        factors.push_back(pot1_POL);

        // Factor 2 is for X3 which will be indirectly observed as being 0
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0, 3u, binDom) );
        factors.push_back(pot1_POL);

        // Factor 3 is for X4 which will be indirectly observed as being 1
        pot1_POL = uniqptr< Polya<AssignType> > (
          new Polya<AssignType>(0, 4u, binDom) );
        factors.push_back(pot1_POL);

        // Factor 4 observes the sum X2+X3 to be 0
        map< vector<AssignType>, FProb > table;
        RVIds vars{2,3,5};
        table[{0u, 0u, 0u}] = 1.0;
        table[{0u, 1u, 1u}] = 1.0;
        table[{1u, 0u, 1u}] = 1.0;
        table[{1u, 1u, 2u}] = 1.0;
        rcptr< DiscreteTable<AssignType> > theCatPtr
          = uniqptr< DiscreteTable<AssignType> >(
            new DiscreteTable<AssignType>(
              vars, {binDom, binDom, cat3Dom},
              0.0, table,
              0.0, 0.0) );
        factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{5u},RVVals{0u}) );

        // Factor 5 observes the sum X3+X4 to be 1
        vars = {3,4,6};
        theCatPtr
          = uniqptr< DiscreteTable<AssignType> >(
            new DiscreteTable<AssignType>(
              vars, {binDom, binDom, cat3Dom},
              0.0, table,
              0.0, 0.0) );
        factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{6u},RVVals{1u}) );

                // nodes
        map<AssignType,rcptr<Factor> > nd;
        nd[0] = factors[0];
        nd[1] = factors[1];
        nd[2] = factors[2];
        nd[3] = factors[3];
        nd[4] = factors[4];
        nd[5] = factors[5];

        // sep sets
        map<vector<AssignType>, RVIds> ss;
        ss[{0,1}] = {0};
        ss[{0,2}] = {0};
        ss[{0,3}] = {0};
        ss[{1,4}] = {2};
        ss[{2,4}] = {3};
        ss[{2,5}] = {3};
        ss[{3,5}] = {4};
        for (auto elem : ss) {
          if (elem.first[0] < elem.first[1]) {
            ss[{elem.first[1],elem.first[0]}] = elem.second;
          } // if
        } // for

        // initial messages
        map<vector<AssignType>, rcptr<Factor> > msg;
        for (auto elem : ss) {
          msg[elem.first] = nd[elem.first[0]]->marginalize(elem.second);
        } // for

        for (unsigned rpt = 0; rpt < 9; rpt++) {

          // manual message passing

          // upward pass
          // std::cout << "\nUpward pass\n";
          msg[{4,1}] = nd[4]->absorb(msg[{2,4}])->marginalize(ss[{4,1}])->normalize();
          msg[{4,2}] = nd[4]->absorb(msg[{1,4}])->marginalize(ss[{4,2}])->normalize();
          msg[{5,2}] = nd[5]->absorb(msg[{3,5}])->marginalize(ss[{5,2}])->normalize();
          msg[{5,3}] = nd[5]->absorb(msg[{2,5}])->marginalize(ss[{5,3}])->normalize();
          msg[{1,0}] = nd[1]->absorb(msg[{4,1}])->marginalize(ss[{1,0}])->normalize();
          msg[{2,0}] = nd[2]->absorb(msg[{4,2}])->absorb(msg[{5,2}])->marginalize(ss[{2,0}])->normalize();
          msg[{3,0}] = nd[3]->absorb(msg[{5,3}])->marginalize(ss[{3,0}])->normalize();
          // std::cout << "\n4 to 1: " << *msg[{4,1}] << std::endl;
          // std::cout << "\n4 to 2: " << *msg[{4,2}] << std::endl;
          // std::cout << "\n5 to 2: " << *msg[{5,2}] << std::endl;
          // std::cout << "\n5 to 3: " << *msg[{5,3}] << std::endl;
          // std::cout << "\n1 to 0: " << *msg[{1,0}] << std::endl;
          // std::cout << "\n2 to 0: " << *msg[{2,0}] << std::endl;
          // std::cout << "\n3 to 0: " << *msg[{3,0}] << std::endl;

          // downward pass
          // std::cout << "\nDownward pass\n";
          msg[{0,1}] = nd[0]->absorb(msg[{2,0}])->absorb(msg[{3,0}])->marginalize(ss[{0,1}])->normalize();
          msg[{0,2}] = nd[0]->absorb(msg[{1,0}])->absorb(msg[{3,0}])->marginalize(ss[{0,2}])->normalize();
          msg[{0,3}] = nd[0]->absorb(msg[{1,0}])->absorb(msg[{2,0}])->marginalize(ss[{0,3}])->normalize();
          msg[{1,4}] = nd[1]->absorb(msg[{0,1}])->marginalize(ss[{1,4}])->normalize();
          msg[{2,4}] = nd[2]->absorb(msg[{0,2}])->absorb(msg[{5,2}])->marginalize(ss[{2,4}])->normalize();
          msg[{2,5}] = nd[2]->absorb(msg[{0,2}])->absorb(msg[{4,2}])->marginalize(ss[{2,5}])->normalize();
          msg[{3,5}] = nd[3]->absorb(msg[{0,3}])->marginalize(ss[{3,5}])->normalize();

          // std::cout << "\n0 to 1: " << *msg[{0,1}] << std::endl;
          // std::cout << "\n0 to 2: " << *msg[{0,2}] << std::endl;
          // std::cout << "\n0 to 3: " << *msg[{0,3}] << std::endl;
          // std::cout << "\n1 to 4: " << *msg[{1,4}] << std::endl;
          // std::cout << "\n2 to 4: " << *msg[{2,4}] << std::endl;
          // std::cout << "\n2 to 5: " << *msg[{2,5}] << std::endl;
          // std::cout << "\n3 to 5: " << *msg[{3,5}] << std::endl;

        } // for

        // check that the dirichlet alphas are correct
        tst_FAC =
          nd[0]->
          absorb(msg[{1,0}])->
          absorb(msg[{2,0}])->
          absorb(msg[{3,0}])->
          normalize();
        tst2_FAC = tst_FAC->marginalize({0u});
        tst_DIR = dynamic_pointer_cast<DIR>(tst2_FAC);
        if ( true and (tst_DIR->getAlpha({0}) !=  6 or tst_DIR->getAlpha({1}) != 4) ) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 2 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 6 and 4 alpha counts respectively, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check that the categorical probs are correct
        tst2_FAC = tst_FAC->marginalize({1u});
        tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst2_FAC);
        try{
          assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
        } // try
        catch(const exception& e) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 2 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = 0.6 but got " << tst_DT->potentialAt({1u}, {0u}) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // catch

      } // integration test 2 scope

      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      // +++ Integration test 3:                        +++
      // +++ loopy belief on a coin indirectly observed +++
      // ++++++++++++++++++++++++++++++++++++++++++++++++++

      if (false) { // integration test 3 scope

        try{

          factors.clear();

          factors.push_back(dir1_DIR);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0u, 1u, binDom) );
          factors.push_back( pot1_POL);
          //std::cout << factors.size()-1 << ": " << *factors[factors.size()-1] << std::endl;

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 2u, binDom) );
          factors.push_back(pot1_POL);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 3u, binDom) );
          factors.push_back(pot1_POL);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 4u, binDom) );
          factors.push_back(pot1_POL);

          map< vector<AssignType>, FProb > table;
          RVIds vars{2,3,5};
          table[{0u, 0u, 0u}] = 1.0;
          table[{0u, 1u, 1u}] = 1.0;
          table[{1u, 0u, 1u}] = 1.0;
          table[{1u, 1u, 2u}] = 1.0;
          rcptr< DiscreteTable<AssignType> > theCatPtr
            = uniqptr< DiscreteTable<AssignType> >(
              new DiscreteTable<AssignType>(
                vars, {binDom, binDom, cat3Dom},
                0.0, table,
                0.0, 0.0) );
          factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{5u},RVVals{0u}) );

          vars = {3,4,6};
          theCatPtr
            = uniqptr< DiscreteTable<AssignType> >(
              new DiscreteTable<AssignType>(
                vars, {binDom, binDom, cat3Dom},
                0.0, table,
                0.0, 0.0) );
          factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{6u},RVVals{1u}) );

          rcptr<ClusterGraph> cgPtr;
          cgPtr = uniqptr<ClusterGraph>( new ClusterGraph(factors) );
          //cgPtr->exportToGraphViz("blah");

          map< Idx2, rcptr<Factor> > msgs;
          MessageQueue msgQ;

          msgs.clear(); msgQ.clear();
          // NOTE: BU (belief update) is the Lauritzen-Spiegelhalter algo.
          // BP (belief propagation) gives the Shafer-Shenoy algo. The loopy
          // and the query MUST match on these.
          //unsigned nMsgs =
          loopyBP_CG(*cgPtr, msgs, msgQ, 0.0);
          //std::cout << "Sent " << nMsgs << " messages before convergence\n";

          // check that the dirichlet alphas are correct
          tst_FAC = queryLBP_CG( *cgPtr, msgs, {0u} );
          tst_DIR = dynamic_pointer_cast<DIR>(tst_DIR);
          if (tst_DIR->getAlpha({0}) !=  6 or tst_DIR->getAlpha({1}) != 4) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 3 FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 6 and 4 alpha counts respectively, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

          // check that the categorical probs are correct
          tst_FAC = queryLBP_CG( *cgPtr, msgs, {1u} );
          tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
          try{
            assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
          } // try
          catch(const exception& e) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 3 FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = 0.6 but got " << tst_DT->potentialAt({1u}, {0u}) );
            PRLITE_THROW(prlite::TestFailed, "");
          } // catch

        } // try

        catch(const char* msg) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 3 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Caught exception", msg);
          PRLITE_THROW(prlite::TestFailed, "");
        } //  catch

      } // integration test 3 scope

      // ++++++++++++++++++++++++++++++++++++++++++++++++++
      // +++ Integration test 4:                        +++
      // +++ loopy belief on a coin indirectly observed +++
      // ++++++++++++++++++++++++++++++++++++++++++++++++++

      { // integration test 4 scope

        try{

          factors.clear();

          factors.push_back(dir1_DIR);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0u, 1u, binDom) );
          factors.push_back( pot1_POL);
          //std::cout << factors.size()-1 << ": " << *factors[factors.size()-1] << std::endl;

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 2u, binDom) );
          factors.push_back(pot1_POL);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 3u, binDom) );
          factors.push_back(pot1_POL);

          pot1_POL = uniqptr< Polya<AssignType> > (
            new Polya<AssignType>(0, 4u, binDom) );
          factors.push_back(pot1_POL);

          map< vector<AssignType>, FProb > table;
          RVIds vars{2,3,5};
          table[{0u, 0u, 0u}] = 1.0;
          table[{0u, 1u, 1u}] = 1.0;
          table[{1u, 0u, 1u}] = 1.0;
          table[{1u, 1u, 2u}] = 1.0;
          rcptr< DiscreteTable<AssignType> > theCatPtr
            = uniqptr< DiscreteTable<AssignType> >(
              new DiscreteTable<AssignType>(
                vars, {binDom, binDom, cat3Dom},
                0.0, table,
                0.0, 0.0) );
          factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{5u},RVVals{0u}) );

          vars = {3,4,6};
          theCatPtr
            = uniqptr< DiscreteTable<AssignType> >(
              new DiscreteTable<AssignType>(
                vars, {binDom, binDom, cat3Dom},
                0.0, table,
                0.0, 0.0) );
          factors.push_back( theCatPtr->Factor::observeAndReduce(RVIds{6u},RVVals{1u}) );

          rcptr<ClusterGraph> cgPtr = uniqptr<ClusterGraph>( new ClusterGraph(factors) );
          //cgPtr->exportToGraphViz("blah");
          //std::cout << *cgPtr << std::endl;

          map< Idx2, rcptr<Factor> > msgs;
          MessageQueue msgQ;
          msgs.clear(); msgQ.clear();
          // NOTE: BU (belief update) is the Lauritzen-Spiegelhalter algo.
          // BP (belief propagation) gives the Shafer-Shenoy algo. The loopy
          // and the query MUST match on these.
          //nMsgs =
          loopyBU_CG(*cgPtr, msgs, msgQ, 0.0);

          //std::cout << "Sent " << nMsgs << " messages before convergence\n";

          // check that the dirichlet alphas are correct
          tst_FAC = queryLBU_CG( *cgPtr, msgs, {0u} );
          tst_DIR = dynamic_pointer_cast<DIR>(tst_DIR);
          if (tst_DIR->getAlpha({0}) !=  6 or tst_DIR->getAlpha({1}) != 4) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 4 FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected 6 and 4 alpha counts respectively, but got " << tst_DIR->getAlpha({0}) << " and " << tst_DIR->getAlpha({1}));
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

          // check that the categorical probs are correct
          tst_FAC = queryLBU_CG( *cgPtr, msgs, {1u} );
          tst_DT = dynamic_pointer_cast< DiscreteTable<AssignType> > (tst_FAC);
          try{
            assertAlmostEqual(0.6, tst_DT->potentialAt({1u}, {0u}), 5);
          } // try
          catch(const exception& e) {
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 4 FAILED!");
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected DiscreteTable with prob P(0) = 0.6 but got " << tst_DT->potentialAt({1u}, {0u}) );
            PRLITE_THROW(prlite::TestFailed, "");
          } // catch

        } // try

        catch(const char* msg) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya::integration test 4 FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Caught exception", msg);
          PRLITE_THROW(prlite::TestFailed, "");
        } //  catch

      } // integration test 4 scope

#ifdef TESTING_MODE
    FILEPOS;
#endif

    } // try

    catch (const char* msg) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Polya threw an exception!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", msg);
      PRLITE_THROW(prlite::TestFailed, "");
    } // catch


  } // runTests

}; // TestPolya

  PRLITE_REGISTER_TESTCASE(TestPolya);


//#endif // NO_TESTING
