/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    September 2018
          PURPOSE: Unit tests for ConditionalPolya

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


//#ifndef NO_TESTING

// patrec headers
#include "prlite_error.hpp"
#include "prlite_testing.hpp"  // TestCase

// emdw headers
#include "emdw.hpp"
#include "conditionalpolya.hpp"
#include "dirichletset.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <limits>

using namespace std;
using namespace emdw;

class BaseTestConditionalPolya : public prlite::TestCase
{
public:

  /** setup to be run before each test. */
  virtual void setup()
    {

    }


  /** clean up after. */
  virtual void tearDown()
    {

    }

  /** example user test. */

protected:
  // Data to use between setup and tearDown.

};

class TestConditionalPolya : public BaseTestConditionalPolya
{
public:

  /** a name. */
  const string& isA()
    {
      static const string name = "TestConditionalPolya";
      return name;
    }


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

    typedef unsigned ZTYPE;
    typedef unsigned YTYPE;
    typedef Dirichlet<YTYPE> DIR;
    typedef DiscreteTable<YTYPE> DTY;
    typedef DiscreteTable<ZTYPE> DTZ;
    typedef DirichletSet<ZTYPE,YTYPE> CDIR;
    typedef ConditionalPolya<ZTYPE,YTYPE> CPOL;
    RVIdType idCnt = 0;

    rcptr< vector<ZTYPE> > zDom (
      new vector<ZTYPE>{0,1});
    rcptr< vector<YTYPE> > yDom (
      new vector<YTYPE>{0,1});
    rcptr< vector<YTYPE> > xDom (
      new vector<YTYPE>{0,1});


    rcptr<Factor> tst_FAC;
    rcptr<CPOL> tst_CPOL;
    rcptr<CDIR> tst_CDIR;
    rcptr<DIR> tst_DIR;
    rcptr<DTZ> tst_DTZ;
    rcptr<DTY> tst_DTY;

    map<pair<ZTYPE,YTYPE>, double> expectZY;
    map<ZTYPE, double> expectZ;
    map<YTYPE, double> expectY;

    // ********************************
    // *** Build ConditionalPolya 1 ***
    // ********************************

    // ++++++++++++++++++++++++++++++
    // +++ build the DirichletSet +++
    // ++++++++++++++++++++++++++++++

    RVIdType refD0_ID = idCnt++; double b00 = 0.8; double b01 = 1.2;
    rcptr<DIR> refD0_DIR = uniqptr<DIR>(
      new DIR (refD0_ID, yDom, 0.0, map<YTYPE, double> { {0,b00},{1,b01} } ) );
    //cout << "DirichletRef0: " << *refD0_DIR << endl;

    RVIdType refD1_ID = idCnt++; double b10 = 3.1; double b11 = 0.9;
    rcptr<DIR> refD1_DIR = uniqptr<DIR>(
      new DIR (refD1_ID, yDom, 0.0, map<YTYPE, double> { {0,b10},{1,b11} } ) );
    //cout << "DirichletRef1: " << *refD1_DIR << endl;
    RVIdType refCDir_ID = idCnt++;
    map< ZTYPE, rcptr<DIR > > dirsMap { {0u,refD0_DIR},{1u,refD1_DIR} };
    rcptr<CDIR> ref_CDIR = uniqptr<CDIR>( new CDIR (refCDir_ID, yDom, dirsMap) );
    rcptr<Factor> ref_CDIR_FAC = ref_CDIR;

    //cout << *ref_CDIR << endl;

    map<vector<unsigned>, FProb> sparseProbs;

    // +++++++++++++++++++++++++++++++++++++++++++++++
    // +++ build the conditioned-on z distribution +++
    // +++++++++++++++++++++++++++++++++++++++++++++++

    RVIdType refZ_ID = idCnt++;
    sparseProbs.clear();
    sparseProbs[vector<unsigned>{0}] = 0.3;
    sparseProbs[vector<unsigned>{1}] = 0.7;
    rcptr<DTZ> refZ_DTZ =
      uniqptr<DTZ>( new DTZ({refZ_ID}, {zDom}, 0.0, sparseProbs) );
    //cout << "Categorical ZRef: " << *refZ_DTZ << endl;

    // +++++++++++++++++++++++++++++++++++++++++
    // +++ build the external Y distribution +++
    // +++++++++++++++++++++++++++++++++++++++++

    // X conditioned on Y
    RVIdType refY_ID = idCnt++;
    RVIdType refX_ID = idCnt++;
    YTYPE xVal = 1;
    ValidValue<YTYPE> theObservedX;
    sparseProbs.clear();
    sparseProbs[vector<unsigned>{0,0}] = 0.9;
    sparseProbs[vector<unsigned>{0,1}] = 0.1;
    sparseProbs[vector<unsigned>{1,0}] = 0.4;
    sparseProbs[vector<unsigned>{1,1}] = 0.6;
    rcptr<Factor> refYGivenX_FAC =
      uniqptr<DTY>( new DTY ({refY_ID,refX_ID}, {yDom,xDom}, 0.0, sparseProbs) )
      ->observeAndReduce({refX_ID},{xVal},true,nullptr)->normalize();
    rcptr<DTY> refYGivenX_DTY = dynamic_pointer_cast<DTY>(refYGivenX_FAC);
    //cout << "Categorical X|Y: " << *refYGivenX_DTY << endl;

    // +++++++++++++++++++++++++++++++++++
    // +++ build the conditional Polya +++
    // +++++++++++++++++++++++++++++++++++

    rcptr<CPOL> ref0_CPOL = uniqptr<CPOL>( new CPOL (refCDir_ID, refZ_ID, refY_ID, zDom, yDom) );
    rcptr<Factor> ref0_CPOL_FAC = ref0_CPOL;
    //cout << *ref0_CPOL << endl;

    rcptr<CPOL> ref1_CPOL = uniqptr<CPOL>(
      new CPOL (refCDir_ID, refZ_ID, refY_ID, zDom, yDom, theObservedX, ref_CDIR, refZ_DTZ, refYGivenX_DTY) );
    rcptr<Factor> ref1_CPOL_FAC = ref1_CPOL;
    //cout << *ref1_CPOL << endl;

    rcptr<CPOL> ref2_CPOL = uniqptr<CPOL>(
      new CPOL (refCDir_ID, refZ_ID, refY_ID, zDom, yDom, YTYPE(1) ) );
    rcptr<Factor> ref2_CPOL_FAC = ref2_CPOL; //!!!

    // ********************************
    // *** Build ConditionalPolya 0 ***
    // ********************************

    // #####################
    // ### CLASS MEMBERS ###
    // #####################


    // ********************
    // *** Constructors ***
    // ********************

    // *****************************
    // *** Polymorphic functions ***
    // *****************************

    // *********************
    // *** Other members ***
    // *********************


    // ########################
    // ### FACTOR OPERATORS ###
    // ########################

    // **************
    // *** sample ***
    // **************

    /*
    size_t nSamples = 100000;
    size_t cnt = 0;
    vector<RVVals> allSamples;
    while (cnt < nSamples) {
      RVVals sample = ref1_CPOL->sample(nullptr);
      YTYPE yRefVal = YTYPE(sample[3]);
      RVVals xRefSample = refYGivenX_FAC->observeAndReduce({refY_ID}, {yRefVal})->normalize()->sample();
      YTYPE xRefVal = YTYPE(xRefSample[0]);
      if (xRefVal == 1) {
        allSamples.push_back(sample);
        cnt++;
      } // if
    } // while

    double zy00,zy01,zy10,zy11;
    zy00=zy01=zy10=zy11 = 0;
    for (size_t n = 0; n < allSamples.size(); n++) {
      RVVals& sample(allSamples[n]);
      ZTYPE zVal = ZTYPE(sample[0]);
      YTYPE yVal = YTYPE(sample[sample.size()-1]);
      zy00 += (zVal == 0 and yVal == 0);
      zy01 += (zVal == 0 and yVal == 1);
      zy10 += (zVal == 1 and yVal == 0);
      zy11 += (zVal == 1 and yVal == 1);
    } // for
    cout <<  zy00 << " " << zy01 << " " << zy10<< " " << zy11 << endl;
    */

    // *****************
    // *** normalize ***
    // *****************

    // **************
    // *** absorb ***
    // **************

    rcptr<CPOL> tmp_CPOL;
    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);

    tst_FAC =
      tmp_CPOL->
      normalize();

    if (*tst_FAC !=  *ref1_CPOL_FAC) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *ref1_CPOL_FAC << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // **************
    // *** cancel ***
    // **************

    tmp_CPOL = uniqptr<CPOL>( ref1_CPOL->copy() );
    tmp_CPOL->Factor::inplaceCancel(ref_CDIR);
    tmp_CPOL->Factor::inplaceCancel(refZ_DTZ);
    tmp_CPOL->Factor::inplaceCancel(refYGivenX_DTY);

    tst_FAC =
      tmp_CPOL->
      normalize();
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *ref0_CPOL_FAC) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::cancel FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *ref0_CPOL_FAC << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // **********************
    // *** observe/reduce ***
    // **********************

    tst_FAC =
      ref1_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)});
    tst_CPOL = dynamic_pointer_cast<CPOL>(tst_FAC);
    //cout << *tst_CPOL << endl;


    auto ret = tst_CPOL->getObservedX();
    if (!ret.valid or ret.value != 1) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::observeAndReduce FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: 1; Got: " << ret.value);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ***************************************************
    // *** marginalize: Keep Phi (the DirichletSet RV) ***
    // ***************************************************

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++ y not observed and also no external y distribution +++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tst_FAC = ref0_CPOL->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 0.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *ref_CDIR) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *ref_CDIR << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 1.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);

    expectZY.clear();
    expectZY[make_pair(0,0)] = 0.92000; expectZY[make_pair(0,1)] = 1.38000;
    expectZY[make_pair(1,0)] = 3.64250; expectZY[make_pair(1,1)] = 1.05750;

    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      //cout << *tst_DIR << endl;

      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-expectZY[make_pair(zIdx,yIdx)]) > 1E-05) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = "
                   << expectZY[make_pair(zIdx,yIdx)]
                   << ", but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // +++++++++++++++++++++++++++++++
    // +++ external y distribution +++
    // +++++++++++++++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);

    tst_FAC = tmp_CPOL->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 1.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *ref_CDIR) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *ref_CDIR << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 1.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refCDir_ID});
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    //cout << *tst_FAC << endl;

    expectZY.clear();
    expectZY[make_pair(0,0)] = 0.844651; expectZY[make_pair(0,1)] = 1.60186;
    expectZY[make_pair(1,0)] = 3.30186; expectZY[make_pair(1,1)] =  1.25163;

    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if ( abs(tst_DIR->getAlpha(yIdx)-expectZY[make_pair(zIdx,yIdx)]) > 1E-05) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = " << expectZY[make_pair(zIdx,yIdx)] << ", but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for z
    } // for y

    // ++++++++++++++++++
    // +++ observed y +++
    // ++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tst_FAC = ref0_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      Factor::marginalize({refCDir_ID})->
      Factor::normalize();
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 1.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *ref_CDIR) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *ref_CDIR << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;

    // should be flat at alpha = 1.0
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);
    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-1.0) > 1E-5) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = 1.0, but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refCDir_ID});
    //cout << *tst_FAC << endl;
    tst_CDIR = dynamic_pointer_cast<CDIR>(tst_FAC);

    expectZY.clear();
    expectZY[make_pair(0,0)] = 0.8; expectZY[make_pair(0,1)] = 1.73333;
    expectZY[make_pair(1,0)] = 3.1; expectZY[make_pair(1,1)] = 1.36667;

    for (auto zIdx : *zDom) {
      tst_DIR = tst_CDIR->aliasDirichlet(zIdx);
      //cout << *tst_DIR << endl;

      for (auto yIdx : *yDom) {
        if (abs(tst_DIR->getAlpha(yIdx)-expectZY[make_pair(zIdx,yIdx)]) > 1E-05) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: alpha = "
                   << expectZY[make_pair(zIdx,yIdx)]
                   << ", but Got: " << tst_DIR->getAlpha(yIdx) );
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      } // for
    } // for

    // ************************************************************
    // *** marginalize: Keep z (the conditioned-on categorical) ***
    // ************************************************************

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++ y not observed and also no external y distribution +++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tst_FAC = ref0_CPOL->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *refZ_DTZ) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refZ_DTZ << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *refZ_DTZ) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refZ_DTZ << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // +++++++++++++++++++++++++++++++
    // +++ external y distribution +++
    // +++++++++++++++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);

    tst_FAC = tmp_CPOL->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *refZ_DTZ) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refZ_DTZ << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refZ_ID});
    tst_DTZ = dynamic_pointer_cast<DTZ>(tst_FAC);
    //cout << *tst_FAC << endl;

    expectZ.clear();
    expectZ[0] = 0.446512;
    expectZ[1] = 0.553488;
    for (auto zIdx : *zDom) {
      if ( abs(tst_DTZ->getProb({zIdx}) - expectZ[zIdx]) > 1E-05 ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: z-prob = " << expectZ[zIdx]
                 << ", but Got: " << tst_DTZ->getProb({zIdx}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // if
    } // for z

    // ++++++++++++++++++
    // +++ observed y +++
    // ++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tst_FAC = ref0_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      Factor::marginalize({refZ_ID})->
      Factor::normalize();
    //cout << *tst_FAC << endl;

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    tst_DTZ = dynamic_pointer_cast<DTZ>(tst_FAC);
    //cout << *tst_FAC << endl;

    expectZ.clear();
    expectZ[0] = 0.727273;
    expectZ[1] = 0.272727;
    for (auto zIdx : *zDom) {
      if ( abs(tst_DTZ->getProb({zIdx}) - expectZ[zIdx]) > 1E-05 ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: z-prob = " << expectZ[zIdx]
                 << ", but Got: " << tst_DTZ->getProb({zIdx}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // if
    } // for z

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *refZ_DTZ) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refZ_DTZ << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      Factor::observeAndReduce({refY_ID},{YTYPE(1)})->
      normalize()->
      Factor::marginalize({refZ_ID});
    //cout << *tst_FAC << endl;
    tst_DTZ = dynamic_pointer_cast<DTZ>(tst_FAC);

    expectZ.clear();
    expectZ[0] = 0.533333;
    expectZ[1] = 0.466667;
    for (auto zIdx : *zDom) {
      if ( abs(tst_DTZ->getProb({zIdx}) - expectZ[zIdx]) > 1E-05 ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected z-prob = " << expectZ[zIdx]
                 << ", but Got: " << tst_DTZ->getProb({zIdx}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // if
    } // for z

    // ****************************************************
    // *** marginalize: Keep y (the target categorical) ***
    // ****************************************************

    // +++++++++++++++++++++++++++++++
    // +++ external y distribution +++
    // +++++++++++++++++++++++++++++++

    // ----------------------------------
    // --- no existing z distribution ---
    // ----------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);

    tst_FAC = tmp_CPOL->
      Factor::marginalize({refY_ID});
    //cout << *tst_FAC << endl;

    // should be refY
    if (*tst_FAC != *refYGivenX_DTY) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refYGivenX_DTY << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refY_ID});
    //cout << *tst_FAC << endl;

    // should be refY
    if (*tst_FAC != *refYGivenX_DTY) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected: " << *refYGivenX_DTY << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // --------------------------------------
    // --- has an existing z distribution ---
    // --------------------------------------

    // ................................
    // ... no existing DirichletSet ...
    // ................................

        tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refY_ID});
    //cout << *tst_FAC << endl;

    if (*tst_FAC !=  *refYGivenX_DTY) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::absorb FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "\nExpected:\n" << *refYGivenX_DTY << "\nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

    // ....................................
    // ... has an existing DirichletSet ...
    // ....................................

    tmp_CPOL = uniqptr<CPOL>( ref0_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(refYGivenX_DTY);
    tmp_CPOL->Factor::inplaceAbsorb(ref_CDIR);
    tmp_CPOL->Factor::inplaceAbsorb(refZ_DTZ);

    tst_FAC = tmp_CPOL->
      normalize()->
      Factor::marginalize({refY_ID});
    tst_DTY = dynamic_pointer_cast<DTY>(tst_FAC);
    //cout << *tst_FAC << endl;

    expectY.clear();
    expectY[0] = 0.246512;
    expectY[1] = 0.753488;
    for (auto yIdx : *yDom) {
      if ( abs(tst_DTY->getProb({yIdx}) - expectY[yIdx]) > 1E-05 ){
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: y-prob = " << expectY[yIdx]
                 << ", but Got: " << tst_DTY->getProb({yIdx}) );
        PRLITE_THROW(prlite::TestFailed, "");
      } // if
    } // for y

    // ###########################
    // ### Integration testing ###
    // ###########################

    // **************************************************
    // *** Send a message from a CDIR to CPOL to CDIR ***
    // **************************************************

    rcptr<Factor> msg = ref_CDIR_FAC->marginalize( {refCDir_ID}, true, ref2_CPOL_FAC.get() );
        tmp_CPOL = uniqptr<CPOL>( ref2_CPOL->copy() );
    tmp_CPOL->Factor::inplaceAbsorb(msg);

    tst_FAC = tmp_CPOL->
      normalize()->
      marginalize({refCDir_ID})->
      cancel(msg)->
      normalize();

    // should be vacuous
    if (tst_FAC->distanceFromVacuous() > 0.0) {
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya::marginalize FAILED!");
      PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected: vacuous distribution, \nGot:\n" << *tst_FAC);
      PRLITE_THROW(prlite::TestFailed, "");
    } // if

        // #################### END OF TESTS ##################################

#ifdef TESTING_MODE
    FILEPOS;
#endif

  } // try

  catch (char msg[]) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  // catch (char const* msg) {
  //   LOGERROR(pLog, "Testing", "ConditionalPolya threw an exception!");
  //   THROW(prlite::TestFailed, msg);
  // } // catch

  catch (const string& msg) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  catch (const exception& e) {
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", "ConditionalPolya threw an exception!");
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", e.what() );
    throw(e);
  } // catch

  catch(...) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "ConditionalPolya threw an exception!");
    PRLITE_THROW(prlite::TestFailed, "");
  } // catch


  } // runTests
};

PRLITE_REGISTER_TESTCASE(TestConditionalPolya);
