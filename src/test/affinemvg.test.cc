/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    June 2022
          PURPOSE: Unit tests for AffineMVG

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


//#ifndef NO_TESTING

// patrec headers
#include "prlite_error.hpp"
#include "prlite_testing.hpp"  // TestCase
#include "prlite_genmat.hpp"  // Matrix
#include "vecset.hpp"
#include "prlite_testing.hpp"


// emdw headers
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "emdw.hpp"
#include "sqrtmvg.hpp"
#include "affinemvg.hpp"
#include "oddsandsods.hpp"
#include "matops.hpp"


// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>
#include <random>

using namespace std;
using namespace emdw;

class BaseTestAffineMVG : public prlite::TestCase
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

class TestAffineMVG : public BaseTestAffineMVG
{
public:

  /** a name. */
  const string& isA()
    {
      static const string name = "TestAffineMVG";
      return name;
    }

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

  try {

#ifdef TESTING_MODE
    FILEPOS;
#endif

    // The gist of it: All the tests will be compared to what the
    // manual result of SqrtMVG would have been.

    if (true) {

      rcptr<SqrtMVG> dwnPtr;
      double mnx = 5.0; double varx = 1.0;   // mean and variance of x
      double beta = 10.0;                    // noise precision - i.e. precision of y around a*x
      double a = 3.0;                        // y = a*x + noise
      double obsvy = 9.0;                    // observed value for y
      double mny = 3.0; double vary = 0.001;   // mean and variance of y msg
      prlite::ColMatrix<double> aMat(1,1); aMat(0,0) = a;
      prlite::ColVector<double> cVec; // keep empty
      prlite::ColMatrix<double> noiseR(1,1); noiseR(0,0) = sqrt(1.0/beta);


      // ************************************************************
      // 1) We'll start with a conditional distribution free from any
      // observations. Then we'll multipy in an x-distribution and
      // marginalise to y.
      // ************************************************************

      prlite::ColVector<double> mnX(1); mnX[0] = mnx;
      prlite::RowMatrix<double> varX(1,1); varX(0,0) = varx;
      rcptr<SqrtMVG> xMsgPtr = std::make_unique<SqrtMVG>(SqrtMVG({0},mnX,varX));

      // +++++++
      // SqrtMVG
      // +++++++

      rcptr<SqrtMVG> xyPtr( SqrtMVG::constructAffineGaussian(*xMsgPtr, aMat, cVec, {2}, noiseR) );
      rcptr<Factor> sqrtYMargPtr = xyPtr->Factor::marginalize({2});

      // +++++++++
      // AffineMVG
      // +++++++++

      rcptr<AffineMVG> affMVGPtr = make_unique<AffineMVG> (
        AffineMVG({0},{2}, aMat, cVec, noiseR) );//, ovy) );
      affMVGPtr->Factor::inplaceAbsorb(xMsgPtr);
      rcptr<Factor> affYMargPtr = affMVGPtr->marginalize({2});

      // ++++++++++++
      // Lets compare
      // ++++++++++++

      if (*sqrtYMargPtr != *affYMargPtr) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Marginalising towards y differs between SqrtMVG and AffineMVG\n" << *sqrtYMargPtr << endl << *affYMargPtr);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      //affMVGPtr->Factor::inplaceAbsorb(yMsgPtr);
      //rcptr<Factor> affFMVGPtr = affMVGPtr->marginalize({0});

      //dwnPtr = dynamic_pointer_cast<SqrtMVG>(yMargPtr);
      //cout << "y: " << dwnPtr->getMean() << dwnPtr->getCov() << endl;

      // ************************************************************
      // 2) The result will be multiplied with an y distribution and
      // then marginalised to x.
      // ************************************************************

      prlite::ColVector<double> mnY(1); mnY[0] = mny;
      prlite::RowMatrix<double> varY(1,1); varY(0,0) = vary;
      rcptr<SqrtMVG> yMsgPtr = std::make_unique<SqrtMVG>(SqrtMVG({2},mnY,varY));

      // +++++++
      // SqrtMVG
      // +++++++

      xyPtr->Factor::inplaceAbsorb(yMsgPtr);
      xyPtr->inplaceNormalize();
      rcptr<Factor> sqrtXMargPtr = xyPtr->Factor::marginalize({2})->normalize();

      // +++++++++
      // AffineMVG
      // +++++++++

      affMVGPtr->Factor::inplaceAbsorb(yMsgPtr);
      rcptr<Factor> affXMargPtr = affMVGPtr->marginalize({2})->normalize();

      // ++++++++++++
      // Lets compare
      // ++++++++++++

      if (*sqrtXMargPtr != *affXMargPtr) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Marginalising towards x differs between SqrtMVG and AffineMVG\n" << *sqrtXMargPtr << endl << *affXMargPtr);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ************************************************************
      // 3) Lastly we'll create a version with an observed y and
      // marginalise to x.
      // ************************************************************

      RVIds ovIds; ovIds.push_back(2);
      RVVals ovVals; ovVals.push_back(obsvy);

      // +++++++
      // SqrtMVG
      // +++++++

      xyPtr = uniqptr<SqrtMVG>( SqrtMVG::constructAffineGaussian(*xMsgPtr, aMat, cVec, {2}, noiseR) );
      sqrtXMargPtr = xyPtr->Factor::observeAndReduce(ovIds, ovVals)->normalize();
      //cout << *sqrtXMargPtr << endl;

      // +++++++++
      // AffineMVG
      // +++++++++

      affMVGPtr = make_unique<AffineMVG> (
        AffineMVG({0},{2}, aMat, cVec, noiseR, ovIds, ovVals) );
      affMVGPtr->Factor::inplaceAbsorb(xMsgPtr);
      //cout << "Affine MVG: " << *affMVGPtr << endl;
      affXMargPtr = affMVGPtr->marginalize({0})->normalize();
      //cout << "Marginalises to: " << *affXMargPtr << endl;

      // ++++++++++++
      // Lets compare
      // ++++++++++++

      if (*sqrtXMargPtr != *affXMargPtr) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Marginalising towards x with obsrved y differs between SqrtMVG and AffineMVG\n" << *sqrtXMargPtr << endl << *affXMargPtr);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

    } // if

    // ***************************************************************
    // 4) This one involves many operations, it really is a type of
    // integration test. Now we are going to build an AffineMVG,
    // absorb an input X-distribution, then multiply in a new
    // X-distribution while dividing out the old one, and compare the
    // result to what if we did the same multiply and divide to the
    // original joint MVG.
    // ***************************************************************

    if (true) {

      // ++++++++++++++++++++++++++++
      // The first the X-message
      // ++++++++++++++++++++++++++++

      prlite::ColVector<double> mnX1(3); mnX1[0] = mnX1[1] = mnX1[2] = 2.0;
      prlite::RowMatrix<double> cvX1(3,3);
      cvX1(0,0) = cvX1(1,1) = cvX1(2,2) = 1.0;
      cvX1(0,1) = cvX1(1,0) = 0.9;
      cvX1(0,2) = cvX1(2,0) = 0.1;
      cvX1(1,2) = cvX1(2,1) = 0.0;

      rcptr<Factor> xMsg1 = std::make_unique<SqrtMVG>(SqrtMVG({0,1,2},mnX1,cvX1));
      //cout << xMsg1->getL() << endl;

      // ++++++++++++++++++++++++++++++
      // Now we build the AffineMVG
      // ++++++++++++++++++++++++++++++

      RVIds ovIds; ovIds.push_back(4);
      RVVals ovVals; ovVals.push_back(3.0);

      prlite::ColMatrix<double> aMat(2,3);
      aMat(0,0) = aMat(0,1) = aMat(1,1) = 1.0;
      aMat(0,2) = aMat(1,0) = 0.0;
      aMat(1,2) = -1.0;
      prlite::ColVector<double> cVec(2);
      cVec[0] = 1.0; cVec[1] = 2.0;

      prlite::ColMatrix<double> noiseR(2,2);
      noiseR(0,0) = noiseR(1,1) = 0.1;
      noiseR(0,1) = noiseR(1,0) = 0.0;

      rcptr<AffineMVG> amvg1 = make_unique<AffineMVG> (
        AffineMVG({0,1,2},{3,4}, aMat, cVec, noiseR, ovIds, ovVals) );

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Multiply in the first X-message and note the joint
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++

      amvg1->Factor::inplaceAbsorb(xMsg1);
      //cout << *amvg1 << endl;
      rcptr<Factor> asFac1 = amvg1->marginalize(amvg1->getVars());
      //cout << *asFac1 << endl;
      rcptr<SqrtMVG> asMVG1 = dynamic_pointer_cast<SqrtMVG>(asFac1);
      //cout << *asMVG1 << endl;
      //cout << "msg1: " << asMVG1->getMean() << asMVG1->getCov() << endl;
      //exit(-1);

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Build a second message and then an AffineMVG based on it
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      prlite::ColVector<double> mnX2 = mnX1*2.0;
      prlite::RowMatrix<double> cvX2 = cvX1*0.5;

      rcptr<Factor> xMsg2 = std::make_unique<SqrtMVG>(SqrtMVG({0,1,2},mnX2,cvX2));

      rcptr<AffineMVG> amvg2 = make_unique<AffineMVG> (
        AffineMVG({0,1,2},{3,4}, aMat, cVec, noiseR, ovIds, ovVals) );

      amvg2->Factor::inplaceAbsorb(xMsg2);
      rcptr<Factor> asFac2 = amvg2->marginalize( amvg2->getVars() );
      //cout << *asFac << endl;
      rcptr<SqrtMVG> asMVG2 = dynamic_pointer_cast<SqrtMVG>(asFac2);
      //cout << "msg2: " << asMVG2->getMean() << asMVG2->getCov() << endl;

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // Compare this to extracting the joint SqrtMVG from the first
      // AffineMVG and then multiplying by the second X-message and
      // divided by the first
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      rcptr<Factor> asFac3 = amvg1->
        Factor::marginalize( amvg1->getVars() )->  // extract the joint
        Factor::absorb(xMsg2)->
        Factor::cancel(xMsg1);
      //cout << *asFac << endl;
      rcptr<SqrtMVG> asMVG3 = dynamic_pointer_cast<SqrtMVG>(asFac3);
      //cout << "msg3: " << asMVG3->getMean() << asMVG3->getCov() << endl;
      cout << "diff: " << asMVG3->getMean() - asMVG2->getMean()
           << asMVG3->getCov() - asMVG2->getCov() << endl;

      // Confirmed that updating an AffineGaussian by multiplying in
      // the new x distribution and THEN dividing out the old x
      // distribution, is equivalent to doing this to the original
      // (joint) SqrtMVG.

      // ++++++++++++
      // Lets compare
      // ++++++++++++

      if ( *asFac3->normalize() != *amvg2->marginalize( amvg2->getVars() )->normalize() ) {
        PRLITE_LOGERROR(
          PRLITE_pLog, "Testing",
          "Integration test differences between AffineMVG and same operations with joint SqrtMVG");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

    }

    // #################### END OF TESTS ####################

#ifdef TESTING_MODE
    FILEPOS;
#endif

  } // try

  catch (char msg[]) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  // catch (char const* msg) {
  //   LOGERROR(pLog, "Testing", "AffineMVG threw an exception!");
  //   THROW(prlite::TestFailed, msg);
  // } // catch

  catch (const string& msg) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  catch (const exception& e) {
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", "AffineMVG threw an exception!");
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", e.what() );
    throw(e);
  } // catch

  catch(...) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "AffineMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, "");
  } // catch


  } // runTests
};

PRLITE_REGISTER_TESTCASE(TestAffineMVG);
