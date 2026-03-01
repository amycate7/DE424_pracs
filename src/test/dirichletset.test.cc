/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    April 2018
          PURPOSE: ConditionalDirichlet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_testing.hpp"

// emdw headers
#include "dirichletset.hpp"

// standard headers
#include <iostream>

using namespace std;
using namespace emdw;

class BaseTestDirichletSet : public prlite::TestCase {
public:

  /** setup to be run before each test. */
  virtual void setup() {} // setup


  /** clean up after. */
  virtual void tearDown() {} // tearDown

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

}; // BaseTestDirichletSet

class TestDirichletSet : public BaseTestDirichletSet {
public:

  /** a name. */
  const string& isA() {
    static const string name = "TestDirichletSet";
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

      typedef string DirType;
      rcptr< vector<DirType> > wordDom (
        new vector<DirType>{"one", "zero"});

      typedef unsigned CondType;
      rcptr< vector<CondType> > topicDom (
        new vector<CondType>{0,1});

      //RVIds theDirVars;
      rcptr<Factor> tstPtr, refPtr;
      rcptr< Dirichlet<DirType> > tstDirPtr, refDirPtr;
      rcptr< DirichletSet<CondType, DirType> > potPtr1, potPtr2, tstPotPtr;

      //--------------------------------------
      //--- the first DirichletSet ---
      //--------------------------------------

      map< CondType, rcptr< Dirichlet<DirType> > > theDirPtrs1;
      // The first topic Dirichlet
      RVIdType theDirVar_a = 101u;
      map<DirType,double> theAlphas1a;
      theAlphas1a["zero"] = 3.0; theAlphas1a["one"] = 3.0;
      rcptr< Dirichlet<DirType> > theDirPtr1a = uniqptr< Dirichlet<DirType> > (
        new  Dirichlet<DirType> (theDirVar_a, wordDom, 0.0, theAlphas1a) );
      theDirPtrs1[0] = theDirPtr1a;
      //theDirVars.push_back(theDirVar_a);

      // The second topic Dirichlet
      RVIdType theDirVar_b= 100u;
      map<DirType,double> theAlphas1b;
      theAlphas1b["zero"] = 2.0; theAlphas1b["one"] = 2.0;
      rcptr< Dirichlet<DirType> > theDirPtr1b = uniqptr< Dirichlet<DirType> > (
        new  Dirichlet<DirType> (theDirVar_b, wordDom, 0.0, theAlphas1b) );
      theDirPtrs1[1] = theDirPtr1b;
      //theDirVars.push_back(theDirVar_b);

      //---------------------------------------
      //--- the second DirichletSet ---
      //---------------------------------------

      map< CondType, rcptr< Dirichlet<DirType> > > theDirPtrs2;
      // The first topic Dirichlet
      map<DirType,double> theAlphas2a;
      theAlphas2a["zero"] = 2.5; theAlphas2a["one"] = 1.5;
      rcptr< Dirichlet<DirType> > theDirPtr2a = uniqptr< Dirichlet<DirType> > (
        new  Dirichlet<DirType> (theDirVar_a, wordDom, 0.0, theAlphas2a) );
      theDirPtrs2[0] = theDirPtr2a;

      // The second topic Dirichlet
      map<DirType,double> theAlphas2b;
      theAlphas2b["zero"] = 2.2; theAlphas2b["one"] = 3.0;
      rcptr< Dirichlet<DirType> > theDirPtr2b = uniqptr< Dirichlet<DirType> > (
        new  Dirichlet<DirType> (theDirVar_b, wordDom, 0.0, theAlphas2b) );
      theDirPtrs2[1] = theDirPtr2b;

      // set<DirType> theDomain;
      // for (auto el : theAlphas1a) {theDomain.insert(el.first);} // for
      // rcptr< DirichletSet<CondType, DirType> > tmpPtr(
      //   new DirichletSet<CondType, DirType>(0u, 3, theDomain, 1.0, 0.1) );
      // cout <<  *tmpPtr << endl;
      // cout << *tmpPtr->sliceAt("zero") << endl;

      // ######################################
      // ### lets test the factor operators ###
      // ######################################

      // +++++++++++++++++
      // +++ normalize +++
      // +++++++++++++++++

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs1) );

      tstPtr = potPtr1->Factor::normalize();
      if ( !tstPtr->isEqual( potPtr1.get() ) ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::normalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Normalize should leave the factor unchanged");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // cout << *potPtr1 << endl;
      // cout << *tstPtr << endl;
      // exit(-1);

      // ++++++++++++++
      // +++ absorb +++
      // ++++++++++++++

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs1) );
      potPtr2 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs2) );

      tstPtr = potPtr1->Factor::absorb(potPtr2)->normalize();
      tstPotPtr = std::dynamic_pointer_cast< DirichletSet<CondType, DirType> > (tstPtr);
      tstDirPtr = uniqptr< Dirichlet<DirType> > (tstPotPtr->copyDirichlet({0}) );
      if (tstDirPtr->getAlpha({"one"}) !=  3.5 or tstDirPtr->getAlpha({"zero"}) != 4.5) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha counts 3.5 and 4.5 , but got " << tstDirPtr->getAlpha({"one"})
                 << " and " << tstDirPtr->getAlpha({"zero"}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstDirPtr = uniqptr< Dirichlet<DirType> > (tstPotPtr->copyDirichlet({1}) );
      if (tstDirPtr->getAlpha({"one"}) !=  4.0 or tstDirPtr->getAlpha({"zero"}) != 3.2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha counts 4.0 and 3.2 , but got " << tstDirPtr->getAlpha({"one"})
                 << " and " << tstDirPtr->getAlpha({"zero"}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // --- do a repeat test to check for shallow deep issues

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs1) );
      potPtr2 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs2) );

      tstPtr = potPtr1->Factor::absorb(potPtr2)->normalize();
      tstPotPtr = std::dynamic_pointer_cast< DirichletSet<CondType, DirType> > (tstPtr);

      tstDirPtr = uniqptr< Dirichlet<DirType> > (tstPotPtr->copyDirichlet({0}) );
      if (tstDirPtr->getAlpha({"one"}) !=  3.5 or tstDirPtr->getAlpha({"zero"}) != 4.5) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha counts 3.5 and 4.5 , but got " << tstDirPtr->getAlpha({"one"})
                 << " and " << tstDirPtr->getAlpha({"zero"}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstDirPtr = uniqptr< Dirichlet<DirType> > (tstPotPtr->copyDirichlet({1}) );
      if (tstDirPtr->getAlpha({"one"}) !=  4.0 or tstDirPtr->getAlpha({"zero"}) != 3.2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::absorb FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected alpha counts 4.0 and 3.2 , but got " << tstDirPtr->getAlpha({"one"})
                 << " and " << tstDirPtr->getAlpha({"zero"}));
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++
      // +++ cancel +++
      // ++++++++++++++

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs1) );
      potPtr2 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs2) );

      tstPtr = potPtr1->Factor::absorb(potPtr2)->normalize();
      tstPtr = tstPtr->Factor::cancel(potPtr2)->normalize();

      if ( !tstPtr->isEqual( potPtr1.get() ) ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::cancel test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Absorb followed by Cancel should leave the factor unchanged");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if


      // +++++++++++++++++++
      // +++ marginalize +++
      // +++++++++++++++++++

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom,theDirPtrs1) );

      // keep everything
      tstPtr = potPtr1->Factor::marginalize({1u})->normalize();

      if ( !tstPtr->isEqual( potPtr1.get() ) ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Keeping all variables should leave the factor unchanged");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // keeping subsets is not supported (yet)
      try{
        tstPtr = potPtr1->Factor::marginalize({2u})->normalize();
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::marginalize test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Can not marginalize to an unknown variable");
        PRLITE_THROW(prlite::TestFailed, "");
      } // try

      catch(const exception& e) {
        // if we go through here everything is fine
      } // catch

      // +++++++++++++++++++++
      // +++ inplaceDampen +++
      // +++++++++++++++++++++

      potPtr1 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs1) );
      potPtr2 =
        uniqptr< DirichletSet<CondType, DirType> > (
          new DirichletSet<CondType, DirType>(1u, topicDom, theDirPtrs2) );

      tstPtr = uniqptr<Factor>( potPtr1->copy() );
      double tstDist = tstPtr->inplaceDampen(potPtr2, 0.3);
      tstPotPtr = std::dynamic_pointer_cast< DirichletSet<CondType, DirType> > (tstPtr);

      tstDirPtr = tstPotPtr->copyDirichlet({0});
      refPtr = uniqptr<Factor>( theDirPtr1a->copy() );
      double refDist = refPtr->Factor::inplaceDampen(theDirPtr2a, 0.3);
      refDirPtr = std::dynamic_pointer_cast< Dirichlet<DirType> > (refPtr);

      if ( !tstDirPtr->isEqual( refDirPtr.get() ) ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::inplaceDampen test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Result differs from what the separate Dirichlets dictate");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      tstDirPtr = tstPotPtr->copyDirichlet({1});
      refPtr = uniqptr<Factor>( theDirPtr1b->copy() );
      refDist += refPtr->Factor::inplaceDampen(theDirPtr2b, 0.3);
      refDirPtr = std::dynamic_pointer_cast< Dirichlet<DirType> > (refPtr);

      if ( !tstDirPtr->isEqual( refDirPtr.get() ) ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::inplaceDampen test FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Result differs from what the separate Dirichlets dictate");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      try {
        assertAlmostEqual(tstDist, refDist, 5);
      } // try
      catch(const exception& e) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet::inplaceDampen FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Distance differs from what the separate Dirichlets dictate");
        PRLITE_THROW(prlite::TestFailed, "");
      } // catch

#ifdef TESTING_MODE
    FILEPOS;
#endif

    } // try

    catch (const char* msg) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "DirichletSet threw an exception!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", msg);
      PRLITE_THROW(prlite::TestFailed, "");
    } // catch

  } // runTests

}; // TestDirichletSet

PRLITE_REGISTER_TESTCASE(TestDirichletSet);
