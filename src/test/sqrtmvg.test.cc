/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Unit tests for SqrtMVG

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/


//#ifndef NO_TESTING

// patrec headers
#include "prlite_error.hpp"
#include "prlite_testing.hpp"  // TestCase
#include "prlite_genmat.hpp"  // Matrix
#include "vecset.hpp"
#include "prlite_testing.hpp"
//#include "prlite_mean_cov.hpp" //  MLEMeanAndFullCovariance_N mlemac;


// emdw headers
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "emdw.hpp"
#include "sqrtmvg.hpp"
#include "oddsandsods.hpp"
#include "lbp_cg.hpp"
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

class BaseTestSqrtMVG : public prlite::TestCase
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

class TestSqrtMVG : public BaseTestSqrtMVG
{
public:

  /** a name. */
  const string& isA()
    {
      static const string name = "TestSqrtMVG";
      return name;
    }

  //returns a matrix filled with random values chosen from a normal distribution
  prlite::ColVector<double> normVec(int sizeX, default_random_engine generator, double mean, double var){

    prlite::ColVector<double> randvec(sizeX);

    normal_distribution<double> dist(mean,sqrt(var));

    for (int i=0; i < sizeX; i++){
      randvec[i] = dist(generator);
    }
    return randvec;
  }

  /*Reads the output of a SqrtMVG and returns it's mean */
  prlite::ColVector<double> readMean(rcptr<Factor> factor){

    prlite::ColVector<double> mean;
    prlite::RowMatrix<double> cov;

    emdw::RVIds theVars;
    prlite::RowMatrix<double> K;
    prlite::ColVector<double> h;

//          cout << *factor << endl;

    stringstream file;
    file << *factor;
    string dummy;

    file >> dummy;
    PRLITE_CHECK(file >> theVars, IOError, "Error reading vars RowVector.");

    file >> dummy;
    PRLITE_CHECK(file >> K, IOError, "Error reading K.");

    file >> dummy;
    PRLITE_CHECK(file >> h, IOError, "Error reading h.");

    int fail;
    cov = inv(K,fail);

    mean = cov*h;
    return mean;
  }

  /*Compares how difference in distance vector Matrixes are*/
  double maxdif( prlite::ColVector<double>  vec1, prlite::ColVector<double>  vec2 ){

    double maxdif = 0;
    double dif;

    for(int i=0; i < vec1.size(); i++){
      dif = fabs(vec1[i] - vec2[i]) / min(vec1[i],vec2[i]);
      if (dif > maxdif) maxdif = dif;
    }

    return maxdif;
  }

  std::map<emdw::RVIdType, AnyType>
  mapObserved(
    const emdw::RVIds& observedVars,
    const emdw::RVVals& observedVals){

    std::map<emdw::RVIdType, AnyType> observed;
    for (unsigned i = 0; i < observedVars.size(); i++) {
      observed[ observedVars[i] ] = observedVals[i];
    } // for

    return observed;

  } // mapObserved

  /*
    Octave code for generating verification data:

format long e

cv1 = [ 4.0  2.0 -2.0;
        2.0  5.0 -5.0;
       -2.0 -5.0  8.0];
mn1 = [1.0 2.0 4.0]';
K1 = inv(cv1);
L1 = chol(K1)';
h1 = K1 * mn1;
g1 = -0.5*mn1'*K1*mn1-log(sqrt(det(2*pi*cv1)));

K1x = [K1 , zeros(3,1)];
K1x = [K1x; zeros(1,4)];
h1x = [h1;0];
g1x = g1;

K2 = [1.0  0.5;
      0.5  2.0];
h2 = [1.0  2.0]';
cv2 = inv(K2);
L2 = chol(K2)';
mn2 = cv2*h2;
g2 = -0.5*mn2'*K2*mn2-log(sqrt(det(2*pi*cv2)));

K2e = [K2(1,:); zeros(1,2); K2(2,:)];
K2e = [K2e(:,1), zeros(3,1), K2e(:,2)]
h2e = [h2(1), 0, h2(2)]';
g2e = g2;

K2x = [zeros(1,2); K2(1,:); zeros(1,2); K2(2,:)];
K2x = [zeros(4,1), K2x(:,1), zeros(4,1), K2x(:,2)];
h2x = [0, h2(1), 0, h2(2)]';
g2x = g2;

K3 = K1 + K2e;
h3 = h1 + h2e;
g3 = g1 + g2e;
L3 = chol(K3)';

K5 = K1x + K2x;
h5 = h1x + h2x;
g5 = g1x + g2x;
L5 = chol(K5)';

K7 = [K1, zeros(3,2);
      zeros(2,3) , K2];
h7 = [h1; h2];
g7 = g1+g2;
L7 = chol(K7)';
4  */

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

  try {

#ifdef TESTING_MODE
    FILEPOS;
#endif

    if (true) {

      double tol = 1e-9;

      // ++++++++++++++++++++ ref1{0,1,2} ++++++++++++++++++++

      int dim = 3;
      prlite::RowMatrix<double> ref1Cv(dim,dim);
      ref1Cv(0,0) = 4; ref1Cv(1,1) = 5; ref1Cv(2,2) = 8;
      ref1Cv(0,1) = ref1Cv(1,0) = 2; ref1Cv(0,2) = ref1Cv(2,0) = -2; ref1Cv(1,2) = ref1Cv(2,1) = -5;
      //cout << "ref1Cv: " << ref1Cv << endl;

      prlite::ColVector<double> ref1Mn(dim);
      ref1Mn[0] = 1; ref1Mn[1] = 2; ref1Mn[2] = 4;
      //cout << "ref1Mn: " << ref1Mn << endl;

      prlite::ColMatrix<double> ref1L(dim,dim);
      ref1L(0,0) =  0.559016994374947; ref1L(0,1) = 0.0;               ref1L(0,2) = 0.0;
      ref1L(1,0) = -0.223606797749979; ref1L(1,1) = 0.730296743340222; ref1L(1,2) = 0.0;
      ref1L(2,0) = 0.0;                ref1L(2,1) = 0.456435464587639; ref1L(2,2) = 0.353553390593274;

      prlite::ColVector<double> ref1h(3);
      ref1h[0] = 0.0625; ref1h[1] = 2.375; ref1h[2] = 2.0;

      double ref1g = -11.098666105068;

      prlite::RowMatrix<double> ref1K; ref1K = ref1L * ref1L.transpose();

      SqrtMVG ref1({0,1,2}, ref1Mn, ref1Cv, false, true);  // from covariance form

      // ++++++++++++++++++++ ref2{0,2} with scope subset of ref1{0,1,2} ++++++++++++++++++++

      prlite::ColMatrix<double> ref2L(2,2);
      ref2L(0,0) = 1.0;  ref2L(0,1) = 0.0;
      ref2L(1,0) = 0.5;  ref2L(1,1) = 1.32287565553230;

      prlite::ColVector<double> ref2h(2); ref2h[0] = 1.0; ref2h[1] = 2.0;

      double ref2g = -2.70092631529878;

      prlite::RowMatrix<double> ref2K; ref2K = ref2L * ref2L.transpose();

      SqrtMVG ref2({0,2}, ref2L, ref2h, ref2g, false, true); // from triangular canonical form

      // ++++++++++++++++++++ ref3{0,1,2} = ref1{0,1,2}*ref2{0,2}  ++++++++++++++++++++

      prlite::ColMatrix<double> ref3L(3,3);
      ref3L(0,0) =  1.14564392373896e+00;  ref3L(0,1) = 0.0;                  ref3L(0,2) = 0.0;
      ref3L(1,0) = -1.09108945117996e-01;  ref3L(1,1) = 7.55928946018455e-01; ref3L(1,2) = 0.0;
      ref3L(2,0) =  4.36435780471985e-01;  ref3L(2,1) = 5.03952630678970e-01; ref3L(2,2) = 1.37436854187255e+00;
      //cout << __FILE__ << __LINE__ << " ref3L: " << ref3L << endl;

      prlite::ColVector<double> ref3h(3);
      ref3h[0] = 1.0625; ref3h[1] = 2.375; ref3h [2] = 4.0;

      double ref3g = -1.37995924203667e+01;

      SqrtMVG ref3({0,1,2}, ref3L, ref3h, ref3g, false, true); // from triangular canonical form
      //cout << __FILE__ << __LINE__ << " ref3: " << ref3 << endl;

      // ++++++++++++++++++++ ref4{1,3} scope partially overlapping ref1{0,1,2} ++++++++++++++++++++

      // = ref2{0,2} but with changed scope to partially overlap ref1{0,1,2}
      SqrtMVG ref4({1,3}, ref2L, ref2h, ref2g, false, true); // from triangular canonical form

      // ++++++++++++++++++++ ref5{0,1,2,3} = ref1{0,1,2}*ref4{1,3}  ++++++++++++++++++++

      prlite::ColMatrix<double> ref5L(4,4);
      ref5L(0,0) =  5.59016994374947e-01; ref5L(0,1) = 0.0;                  ref5L(0,2) =  0.0;                  ref5L(0,3) = 0.0;
      ref5L(1,0) = -2.23606797749979e-01; ref5L(1,1) = 1.23827837473378e+00; ref5L(1,2) =  0.0;                  ref5L(1,3) = 0.0;
      ref5L(2,0) =  0.0;                  ref5L(2,1) = 2.69190951029083e-01; ref5L(2,2) =  5.10753918455249e-01; ref5L(2,3) = 0.0;
      ref5L(3,0) =  0.0;                  ref5L(3,1) = 4.03786426543624e-01; ref5L(3,2) = -2.12814132689687e-01; ref5L(3,3) = 1.33853153368408e+00;

      prlite::ColVector<double> ref5h(4);
      ref5h[0] = 6.25e-02;
      ref5h[1] = 3.375;
      ref5h[2] = 2.0;
      ref5h[3] = 2.0;

      double ref5g = -1.37995924203667e+01;

      SqrtMVG ref5({0,1,2,3}, ref5L, ref5h, ref5g, false, true); // from triangular canonical form

      // ++++++++++++++++++++ ref6{3,4} scope partially overlapping ref1{0,1,2} ++++++++++++++++++++

      // = ref2{0,2} but with changed scope to not overlap ref1{0,1,2}
      SqrtMVG ref6({3,4}, ref2L, ref2h, ref2g, false, true); // from triangular canonical form
      prlite::RowMatrix<double> ref6L(dim,dim); ref6L = ref6.getL();
      prlite::ColVector<double> ref6h(dim); ref6h = ref6.getH();
      double ref6g = ref6.getG();

      // ++++++++++++++++++++ ref7{0,1,2,3,4} = ref1{0,1,2} * ref6{3,4} ++++++++++++++++++++

      prlite::ColMatrix<double> ref7L(5,5);
      ref7L(0,0) =  5.59016994374947e-01; ref7L(0,1) = 0.0;                  ref7L(0,2) = 0.0;                  ref7L(0,3) = 0.0; ref7L(0,4) = 0.0;
      ref7L(1,0) = -2.23606797749979e-01; ref7L(1,1) = 7.30296743340222e-01; ref7L(1,2) = 0.0;                  ref7L(1,3) = 0.0; ref7L(1,4) = 0.0;
      ref7L(2,0) = -0.0;                  ref7L(2,1) = 4.56435464587639e-01; ref7L(2,2) = 3.53553390593274e-01; ref7L(2,3) = 0.0; ref7L(2,4) = 0.0;
      ref7L(3,0) =  0.0;                  ref7L(3,1) = 0.0;                  ref7L(3,2) = 0.0;                  ref7L(3,3) = 1.0; ref7L(3,4) = 0.0;
      ref7L(4,0) =  0.0;                  ref7L(4,1) = 0.0;                  ref7L(4,2) = 0.0;                  ref7L(4,3) = 0.5; ref7L(4,4) = 1.32287565553230e+00;

      prlite::ColVector<double> ref7h(5);
      ref7h[0] = 6.25e-02;
      ref7h[1] = 2.375;
      ref7h[2] = 2.0;
      ref7h[3] = 1.0;
      ref7h[4] = 2.0;

      double ref7g = -1.37995924203667e+01;

      SqrtMVG ref7({0,1,2,3,4}, ref7L, ref7h, ref7g, false, true); // from triangular canonical form

      // ++++++++++++++++++++ ref8{0,2} ++++++++++++++++++++

      // ref8{0,2} is the marginal of ref1{0,1,2}
      dim = 2;
      prlite::RowMatrix<double> ref8Cv(dim,dim);
      ref8Cv(0,0) = 4; ref8Cv(1,1) = 8;
      ref8Cv(0,1) = ref8Cv(1,0) = -2;
      //cout << "ref8Cv: " << ref8Cv << endl;

      prlite::ColVector<double> ref8Mn(dim);
      ref8Mn[0] = 1; ref8Mn[1] = 4;
      //cout << "ref8Mn: " << ref8Mn << endl;

      SqrtMVG ref8({0,2}, ref8Mn, ref8Cv, false, true);  // from covariance form
      prlite::RowMatrix<double> ref8L(dim,dim); ref8L = ref8.getL();
      prlite::ColVector<double> ref8h(dim); ref8h = ref8.getH();
      double ref8g = ref8.getG();

      double diff = 0.0;
      prlite::ColMatrix<double> theL;
      prlite::ColVector<double> theH;
      double theG;

      // #################### CLASS MEMBERS ####################

      // ******************** Constructors ********************

      // ++++++++++++++++++++ Construct from covariance form, unsorted ++++++++++++++++++++

      if (true){
        dim = 3;
        prlite::RowMatrix<double> tstCov(dim,dim);
        tstCov(0,0) = 8; tstCov(1,1) = 5; tstCov(2,2) = 4;
        tstCov(0,1) = tstCov(1,0) = -5; tstCov(0,2) = tstCov(2,0) = -2; tstCov(1,2) = tstCov(2,1) = 2;

        prlite::ColVector<double> tstMn(dim);
        tstMn[0] = 4; tstMn[1] = 2; tstMn[2] = 1;

        SqrtMVG ref1Alt({2,1,0}, tstMn, tstCov, false, false);
        diff = 0.0;
        theL = ref1Alt.getL();
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Constructing from unsorted covariance : L differs from sorted covariance version:\nExpected: "
                   << ref1L << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        theH = ref1Alt.getH();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref1h[k] - theH[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Constructing from unsorted covariance : h differs from sorted covariance version:\nExpected: "
                   << ref1h << "\nbut got " << theH);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        theG = ref1Alt.getG();
        diff = abs( ref1g - theG );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Constructing from unsorted covariance : g differs from sorted covariance version:\nExpected: "
                   << ref1g << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // ++++++++++++++++++++ Construct from canonical form, sorted ++++++++++++++++++++

      if (true) {
        SqrtMVG ref1Alt({0,1,2}, ref1L, ref1h, ref1g, false, true);

        if (ref1Alt != ref1) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Constructing from sorted canonical differs from sorted covariance version: \nExpected: "
                   << ref1 << "\nbut got: " << ref1Alt);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      }

      // ++++++++++++++++++++ Construct from canonical form, unsorted ++++++++++++++++++++

      {
        prlite::ColMatrix<double> tstL(dim,dim); // cholesky root of RV flipped K1
        tstL(0,0) = 0.577350269189626; tstL(0,1) =  0.0;  tstL(0,2) = 0.0;
        tstL(1,0) = 0.577350269189626; tstL(1,1) =  0.5;  tstL(1,2) = 0.0;
        tstL(2,0) = 0.00000;           tstL(2,1) = -0.25; tstL(2,2) = 0.5;

        prlite::ColVector<double> tstH(3);
        tstH[0] = 2.000000; tstH[1] = 2.375000; tstH[2] = 0.062500;

        SqrtMVG ref1Alt({2,1,0}, tstL, tstH, ref1g, false, false);

        if (ref1Alt != ref1) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Constructing from unsorted canonical differs from sorted covariance version: \nExpected: "
                   << ref1 << "\nbut got: " << ref1Alt);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
      }

      { // Blatta daaliensis - just check for not throwing exceptions
        prlite::ColMatrix<double> myL(4,4);
        myL(0,0) = 1.0; myL(0,1) = 0.0; myL(0,2) = 0.0; myL(0,3) = 0.0;
        myL(1,0) = 0.0; myL(1,1) = 1.0; myL(1,2) = 0.0; myL(1,3) = 0.0;
        myL(2,0) = 0.0; myL(2,1) = 0.0; myL(2,2) = 1.0; myL(2,3) = 0.0;
        myL(3,0) = 0.0; myL(3,1) = 0.0; myL(3,2) = 0.0; myL(3,3) = -1.0;
        prlite::ColVector<double> myH(4);
        myH[0] = myH[1] = myH[2] = myH[3] = 0.0;
        double myG = 0.0;
        SqrtMVG factor2({2,3,0,1}, myL, myH, myG);
      }

      { // another Blatta daaliensis  - just check for not throwing exceptions
        prlite::ColMatrix<double> myL(4,4);
        myL(0,0) =  1.00050; myL(0,1) =  0.00000; myL(0,2) = 0.00000; myL(0,3) = 0.00000;
        myL(1,0) =  0.00000; myL(1,1) =  1.00050; myL(1,2) = 0.00000; myL(1,3) = 0.00000;
        myL(2,0) = -0.99950; myL(2,1) =  0.00000; myL(2,2) = 0.04471; myL(2,3) = 0.00000;
        myL(3,0) =  0.00000; myL(3,1) = -0.99950; myL(3,2) = 0.00000; myL(3,3) = 0.04471;
        prlite::ColVector<double> myH(4);
        myH[0] = myH[1] = myH[2] = myH[3] = 0.0;
        double myG = 0.0;
        SqrtMVG factor2({2,3,0,1}, myL, myH, myG);
      }

      // ******************** Parameter conversions ********************

      // ++++++++++++++++++++ Check if L was calculated correctly ++++++++++++++++++++

      diff = 0.0;
      theL = ref1.getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "When constructed from covariance L deviates too far\nExpected: "
                 << ref1L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Check if h was calculated correctly ++++++++++++++++++++

      diff = 0.0;
      theH = ref1.getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref1h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "When constructed from covariance h deviates too far\nExpected: "
                 << ref1h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Check if g was calculated correctly ++++++++++++++++++++

      diff = abs(ref1.getG() - ref1g);
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "When constructed from covariance g deviates too far\nExpected: "
                 << ref1g << " but got " << ref1.getG() );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Check the sigma points ++++++++++++++++++++

      prlite::ColMatrix<double> sp; sp = ref1.getSigmaPoints();
      prlite::ColVector<double> sx(dim); sx.assignToAll(0.0);
      prlite::ColMatrix<double> sx2(dim,dim); sx2.assignToAll(0.0);

      for (int n = 0; n < sp.cols(); n++) {
        sx += sp.col(n);
        sx2 += sp.col(n)*sp.col(n).transpose();
      } // for

      prlite::ColVector<double> theMn; theMn = sx/sp.cols();
      prlite::RowMatrix<double> theCov; theCov = ( sx2- sx*(sx.transpose()/sp.cols() ) ) / sp.cols();

      diff = 0.0;
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Sigma point mean value deviates too far\nExpected: "
                 << ref1Mn << "\nbut got " << theMn);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1Cv(r,c) - theCov(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Sigma point cov value deviates too far\nExpected: "
                 << ref1Cv << "\nbut got " << theCov);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if


      /*
      prlite::ColMatrix<double> B(dim,dim); B.assignToAll(0.0);
      sx.assignToAll(0.0);
      for (int n = 0; n < sp.cols(); n++) {
        cholUR( B, sp.col(n) );
        sx += sp.col(n);
      } // for
      sx /= sqrt(dim);
      cholDR(B,sx);
      //cout << (B*B.transpose() -  << endl;
      exit(-1);
      */

      // ++++++++++++++++++++ Check if we can recover the mean ++++++++++++++++++++

      diff = 0.0;
      theMn = ref1.getMean();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Recalculated mean value deviates too far\nExpected: "
                 << ref1Mn << "\nbut got " << theMn);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Check if we can recover the covariance ++++++++++++++++++++

      diff = 0.0;
      theCov = ref1.getCov();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1Cv(r,c) - theCov(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Recalculated cov value deviates too far\nExpected: "
                 << ref1Cv << "\nbut got " << theCov);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Check if we can recover the precision ++++++++++++++++++++

      diff = 0.0;
      prlite::RowMatrix<double> theK; theK = ref1.getK();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1K(r,c) - theK(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Recalculated precision value deviates too far\nExpected: "
                 << ref1K << "\nbut got " << theK);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ******************** Virtual functions ********************

      // ++++++++++++++++++++ potentialAt ++++++++++++++++++++

      double refP = 9.16451694961552e-03;
      emdw::RVVals xVals(3);
      xVals[0] = double(1.0); xVals[1] = double(2.0); xVals[2] = double(4.0);
      //xVals[0] = ref1Mn[0]; xVals[1] = ref1Mn[1]; xVals[2] = ref1Mn[2];
      double theP = ref1.potentialAt({0,1,2}, xVals);

      if (abs<double>(theP-refP) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "potentialAt deviates too far\nExpected: "
                 << refP << "\nbut got " << theP);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      refP = 1.51324955224846e-05;
      xVals[0] = double(0.0); xVals[1] = double(0.0); xVals[2] = double(0.0);
      theP = ref1.potentialAt({0,1,2}, xVals);

      if (abs<double>(theP-refP) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "potentialAt deviates too far\nExpected: "
                 << refP << "\nbut got " << theP);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ******************** Other members ********************

      // ++++++++++++++++++++ Mahalanobis distance ++++++++++++++++++++

      double refDist = 12.8125;
      prlite::ColVector<double> x(3);
      x.assignToAll(0.0);
      double theDist = ref1.mahalanobis(x);

      if (abs<double>(theDist-refDist) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Mahalanobis distance deviates too far\nExpected: "
                 << refDist << "\nbut got " << theDist);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ symmetric Kullback-Leibler divergence +++++++++++++++
      refDist = 9.362456747;
      theDist = SqrtMVG::symmKLdivergence(&ref1,&ref3);

      if (abs<double>(theDist-refDist) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "symmetric Kullback-Leibler divergence deviates too far\nExpected: "
                 << refDist << "\nbut got " << theDist);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ++++++++++++++++++++ Lets check linear gaussian joints ++++++++++++++++++++

      if (true) {  // do a linear gaussian

        //----- build the input parameters -----

        prlite::ColMatrix<double> A(3,2); A.assignToAll(0.0);
        A(0,0) = A(1,0) = A(1,1) = 1.0; A(2,1) = -1.0;
        //cout << "A: " << A << endl;

        size_t xDim = ref1.noOfVars(); // ref1 is 3-dim
        size_t yDim = A.cols();
        size_t dim = xDim + yDim;

        //prlite::ColMatrix<double> sptsX = ref1.getSigmaPoints();
        prlite::ColMatrix<double> noiseR(2,2);
        noiseR.assignToAll(0.0); noiseR(0,0) = noiseR(1,1) = 1.0;

        //----- calculate the reference values -----

        prlite::ColMatrix<double> refCv(dim,dim);
        gLinear::gIndexRange rangeX(0, xDim-1);
        gLinear::gIndexRange rangeY(xDim, dim-1);
        refCv.slice(rangeX, rangeX) = ref1.getCov();
        refCv.slice(rangeX, rangeY) = ref1.getCov()*A;
        refCv.slice(rangeY, rangeX) = (ref1.getCov()*A).transpose();
        refCv.slice(rangeY, rangeY) = (A.transpose()*ref1.getCov()*A) + noiseR*noiseR.transpose();
        //cout << "refCv: " << refCv << endl;
        prlite::ColVector<double> refMn(dim);
        refMn.slice(rangeX) = ref1.getMean();
        refMn.slice(rangeY) = A.transpose()*ref1.getMean();
        //cout << "refMn: " << refMn << endl;

        //----- the result returned by SqrtMVG -----

        prlite::ColMatrix<double> B; B = A.transpose(); // B should be 2x3
        // NOTE: there is a weird error with A = A.transpose();

        prlite::ColVector<double> cVec; // keep empty
        SqrtMVG* tmpPtr = SqrtMVG::constructAffineGaussian(
          ref1, B, cVec, {3,4}, noiseR);

        prlite::ColVector<double> theMn; theMn = tmpPtr->getMean();
        prlite::RowMatrix<double> theCv; theCv = tmpPtr->getCov();
        //cout << "theMn: " << theMn << endl << "theCv: " << theCv << endl;
        delete tmpPtr;

        //----- and now we compare them -----

        diff = 0.0;
        for (size_t k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( refMn[k] - theMn[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Mean value deviates too far\nExpected: "
                   << refMn << "\nbut got " << theMn);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        for (size_t r = 0; r < dim; r++) {
          for (size_t c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( refCv(r,c) - theCv(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Cov value deviates too far\nExpected: "
                   << refCv << "\nbut got " << theCv);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // ++++++++++++++++++++ Lets check sigma point joints ++++++++++++++++++++

      if (true) {  // do a linear gaussian via the sigma points route

        //----- build the input parameters -----

        prlite::ColMatrix<double> A(3,2); A.assignToAll(0.0);
        A(0,0) = A(1,0) = A(1,1) = 1.0; A(2,1) = -1.0;
        //cout << "A: " << A << endl;
        prlite::ColVector<double> b(2); b[0] = 1.0; b[1] = 2.0;
        //cout << "b: " << b << endl;

        size_t xDim = ref1.noOfVars();
        size_t yDim = A.cols();
        size_t dim = xDim + yDim;

        prlite::ColMatrix<double> sptsX = ref1.getSigmaPoints();
        //cout << "sptsX: " << sptsX << endl;
        size_t nSamples = sptsX.cols();
        prlite::ColMatrix<double> sptsY; sptsY = A.transpose()*sptsX;
        for (size_t n = 0; n < nSamples; n++) {
          sptsY.col(n) += b;
        } // for
        //cout << "sptsY: " << sptsY << endl;
        prlite::ColMatrix<double> noiseR(2,2);
        noiseR.assignToAll(0.0); noiseR(0,0) = noiseR(1,1) = 1.0;

        //----- calculate the reference values -----

        prlite::ColMatrix<double> refCv(dim,dim);
        gLinear::gIndexRange rangeX(0, xDim-1);
        gLinear::gIndexRange rangeY(xDim, dim-1);
        refCv.slice(rangeX, rangeX) = ref1.getCov();
        refCv.slice(rangeX, rangeY) = ref1.getCov()*A;
        refCv.slice(rangeY, rangeX) = (ref1.getCov()*A).transpose();
        refCv.slice(rangeY, rangeY) = (A.transpose()*ref1.getCov()*A) + noiseR*noiseR.transpose();
        //cout << "refCv: " << refCv << endl;
        prlite::ColVector<double> refMn(dim);
        refMn.slice(rangeX) = ref1.getMean();
        refMn.slice(rangeY) = A.transpose()*ref1.getMean() + b;
        //cout << "refMn: " << refMn << endl;

        //----- the result returned by SqrtMVG -----

        SqrtMVG* tmpPtr = SqrtMVG::constructFromSigmaPoints(
          {0,1,2}, sptsX, {3,4}, sptsY, noiseR);
        prlite::ColVector<double> theMn; theMn = tmpPtr->getMean();
        prlite::RowMatrix<double> theCv; theCv = tmpPtr->getCov();
        //cout << "theMn: " << theMn << endl << "theCv: " << theCv << endl;
        delete tmpPtr;

        //----- and now we compare them -----

        diff = 0.0;
        for (size_t k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( refMn[k] - theMn[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Mean value deviates too far\nExpected: "
                   << refMn << "\nbut got " << theMn);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        for (size_t r = 0; r < dim; r++) {
          for (size_t c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( refCv(r,c) - theCv(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Cov value deviates too far\nExpected: "
                   << refCv << "\nbut got " << theCv);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // #################### FACTOR OPERATORS ####################

      // ******************** Lets sample from it ********************

      // -------------------- non-vacuous --------------------

      // This test might occasionally fail simply due to the
      // randomness of the sampling. Needs a fixed seed value.

      SqrtMVG_Sampler sampler;
      prlite::ColMatrix<double> samples; samples = sampler.sample(&ref1, 100000);
      rcptr<SqrtMVG> blahPtr = uniqptr<SqrtMVG>(SqrtMVG::constructFromSamples(ref1.getVars(), samples));
      theMn = blahPtr->getMean();
      theCov = blahPtr->getCov();

      diff = 0.0;
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
      } // for k

      if (diff > 5E-2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Estimated mean value deviates too far\nExpected: "
                 << ref1Mn << "\nbut got " << theMn);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1Cv(r,c) - theCov(r,c) ) );
        } // for c
      } // for r

      if (diff > 2E-1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Estimated cov value deviates too far\nExpected: "
                 << ref1Cv << "\nbut got " << theCov);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if


      for (unsigned n = 0; n < 100000; n++) {
        emdw::RVVals anySample = ref1.sample(nullptr);
        vector<double> aSample( anySample.size() );
        for (size_t idx = 0; idx < anySample.size(); idx++) {
          aSample[idx] = static_cast<double>(anySample[idx]);
        } // for
        //cout << aSample << endl;
      } // for

      // -------------------- vacuous --------------------

      // ******************** Lets check out the normalizer ********************

      // ref1.adjustMass(0.1);
      // double theMass = ref1.getMass();

      // if (abs<double>(theMass - 0.1) > tol) {
      //   LOGERROR(pLog, "Testing", "adjustMass vs getMass not consistent\nExpected: 0.1\nbut got " << theMass);
      //   LOGERROR(pLog, "Testing", "SqrtMVG Failed");
      //   THROW(prlite::TestFailed, "");
      // } // if

      // -------------------- non-vacuous --------------------

      SqrtMVG_Normalize norm;
      Factor* normedPtr = norm.process(&ref1);
      SqrtMVG* dwnPtr = dynamic_cast<SqrtMVG*>(normedPtr);
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Normalize::process did not return a SqrtMVG!");
      double theMass = dwnPtr->getMass();

      if (abs<double>(theMass - 1.0) > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Normalize returned an unnormalized mass of " << theMass);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- vacuous --------------------

      // ******************** Lets check out the absorber ********************

      SqrtMVG_Absorb absorber;

      // ++++++++++++++++++++ lhs and rhs with same scope ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      rcptr<Factor> newPtr = uniqptr<Factor> ( absorber.process(&ref1, &ref1) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Absorb::process did not return a SqrtMVG!");

      diff = 0.0;
      theK = dwnPtr->getK();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( 2.0*ref1K(r,c) - theK(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with identical scope: precision value deviates too far\nExpected: "
                 << ref1K*2.0 << "\nbut got " << theK);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( 2.0*ref1h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with identical scope: h deviates too far\nExpected: "
                 << ref1h*2.0 << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( 2.0* ref1g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with identical scope: g deviates too far\nExpected: "
                 << ref1g*2.0 << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------
      // -------------------- lhs vacuous, rhs non-vacuous --------------------

      // ++++++++++++++++++++ rhs is subset of lhs ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      newPtr = uniqptr<Factor>( absorber.process(&ref1, &ref2) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Absorb::process did not return a SqrtMVG!");

      diff = 0.0;
      theL = dwnPtr->getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref3L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with rhs subset scope: precision value deviates too far\nExpected: "
                 << ref3L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref3h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with rhs subset scope: h deviates too far\nExpected: "
                 << ref3h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( ref3g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with rhs subset scope: g deviates too far\nExpected: "
                 << ref3g << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------

      // -------------------- lhs vacuous, rhs non-vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.
      if (false) {

        SqrtMVG vac1({0,1,2});
        newPtr = uniqptr<Factor>( absorber.process(&vac1, &ref2) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );

        theL = dwnPtr->getL();
        diff = 0.0;
        diff = max<double>( diff, abs( ref2L(0,0) - theL(0,0) ) );
        diff = max<double>( diff, abs( ref2L(1,0) - theL(2,0) ) );
        diff = max<double>( diff, abs( ref2L(1,1) - theL(2,2) ) );
        diff = max<double>( diff, abs( theL(1,0) ) );
        diff = max<double>( diff, abs( theL(2,1) ) );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with rhs subset scope: precision value deviates too far\ngot " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        theH = dwnPtr->getH();
        diff = 0.0;
        diff = max<double>( diff, abs( ref2h[0] - theH[0] ) );
        diff = max<double>( diff, abs( ref2h[1] - theH[2] ) );
        diff = max<double>( diff, abs( theH[1] ) );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with rhs subset scope: h deviates too far\ngot " << theH);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // ++++++++++++++++++++ lhs is subset of rhs ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      newPtr = uniqptr<Factor>( absorber.process(&ref2, &ref1) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Absorb::process did not return a SqrtMVG!");


      diff = 0.0;
      theL = dwnPtr->getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref3L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with lhs subset scope: precision value deviates too far\nExpected: "
                 << ref3L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref3h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with lhs subset scope: h deviates too far\nExpected: "
                 << ref3h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( ref3g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with lhs subset scope: g deviates too far\nExpected: "
                 << ref3g << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.

      // -------------------- lhs vacuous, rhs non-vacuous --------------------

      // ++++++++++++++++++++ partially overlapping scopes ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      newPtr = uniqptr<Factor>( absorber.process(&ref1, &ref4) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Absorb::process did not return a SqrtMVG!");

      diff = 0.0;
      theL = dwnPtr->getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref5L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: precision value deviates too far\nExpected: "
                   << ref5L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref5h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: h deviates too far\nExpected: "
                   << ref5h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( ref5g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: g deviates too far\nExpected: "
                 << ref5g << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.

      // -------------------- lhs vacuous, rhs non-vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.

      // -------------------- lhs vacuous, rhs vacuous --------------------

      if (true) {

        SqrtMVG vac1({0,1,2});
        SqrtMVG vac2({1,2,3});
        newPtr = uniqptr<Factor>( absorber.process(&vac1, &vac2) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        RVIds expectedVars{0,1,2,3};
        PRLITE_ASSERT(dwnPtr->getVars() == expectedVars, "vars do not match expected");
        PRLITE_ASSERT(dwnPtr->isVacuous(), "The product of two vacuous densities should be vacuous");
      } // if


      // ++++++++++++++++++++ and for completeness, no overlap at all ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      newPtr = uniqptr<Factor>( absorber.process(&ref1, &ref6) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Absorb::process did not return a SqrtMVG!");

      diff = 0.0;
      theL = dwnPtr->getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref7L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: precision value deviates too far\nExpected: "
                   << ref7L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref7h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: h deviates too far\nExpected: "
                   << ref7h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( ref7g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Absorb with partially overlapping scopes: g deviates too far\nExpected: "
                 << ref7g << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.

      // -------------------- lhs vacuous, rhs non-vacuous --------------------

      // This results in a non-vacuous zero-definite
      // configuration. Currently we do not support it.

      // ******************** Lets check out the canceller ********************

      SqrtMVG_Cancel canceller;

      // ++++++++++++++++++++ rhs is subset of lhs ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      // ref3{0,1,2}/ref2{0,2} = ref1{0,1,2}

      newPtr = uniqptr<Factor>( canceller.process(&ref3, &ref2) );
      dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
      PRLITE_ASSERT(dwnPtr, "SqrtMVG_Cancel::process did not return a SqrtMVG!");

      diff = 0.0;
      theL = dwnPtr->getL();
      for (int r = 0; r < dim; r++) {
        for (int c = 0; c <= r; c++) {
          diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
        } // for c
      } // for r

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Cancel with rhs subset scope: precision value deviates too far\nExpected: "
                   << ref3L << "\nbut got " << theL);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      diff = 0.0;
      theH = dwnPtr->getH();
      for (int k = 0; k < dim; k++) {
        diff = max<double>( diff, abs( ref1h[k] - theH[k] ) );
      } // for k

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Cancel with rhs subset scope: h deviates too far\nExpected: "
                   << ref3h << "\nbut got " << theH);
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      theG = dwnPtr->getG();
      diff = abs( ref1g - theG );
      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Cancel with rhs subset scope: g deviates too far\nExpected: "
                 << ref3g << " but got " << theG );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // -------------------- lhs non-vacuous, rhs vacuous --------------------
      // -------------------- lhs vacuous, rhs non-vacuous --------------------
      // -------------------- lhs vacuous, rhs vacuous --------------------

      // ++++++++++++++++++++ this one threw a round-off exception in matops.cc ++++++++++++++++++++

      if (true) {

        prlite::ColMatrix<double> L1(4,4);
        L1(0,0) =  9.59341765932856472e-01; L1(0,1) =  0.00000000000000000e+00; L1(0,2) =  0.00000000000000000e+00; L1(0,3) = 0.00000000000000000e+00;
        L1(1,0) =  8.58410503322889817e-01; L1(1,1) =  3.55952475206043228e-01; L1(1,2) =  0.00000000000000000e+00; L1(1,3) = 0.00000000000000000e+00;
        L1(2,0) = -8.17448641259428066e-01; L1(2,1) = -1.24380987134634340e-02; L1(2,2) =  7.22536872611885084e+00; L1(2,3) = 0.00000000000000000e+00;
        L1(3,0) = -7.37365161706847072e-01; L1(3,1) = -2.86254965830731622e-01; L1(3,2) = -7.23085993512538305e+00; L1(3,3) = 1.84948545345930215e+00;

        prlite::ColMatrix<double> L2(4,4);
        L2(0,0) =  8.85098460023062250e-01; L2(0,1) =  0.00000000000000000e+00; L2(0,2) =  0.00000000000000000e+00; L2(0,3) = 0.00000000000000000e+00;
        L2(1,0) =  7.82529114830719630e-01; L2(1,1) =  3.50771824652251651e-01; L2(1,2) =  0.00000000000000000e+00; L2(1,3) = 0.00000000000000000e+00;
        L2(2,0) = -8.60834951277267058e-01; L2(2,1) = -1.58397057498013666e-01; L2(2,2) =  1.93787137041490243e+00; L2(2,3) = 0.00000000000000000e+00;
        L2(3,0) = -8.25230182052027028e-01; L2(3,1) = -1.86079888457984427e-01; L2(3,2) = -1.25530344571998276e+00; L2(3,3) = 1.53131107633148744e+00;

        prlite::ColVector<double> h(4); h.assignToAll(0.0);
        double g1 = SqrtMVG::getNormedG(L1, h);
        double g2 = SqrtMVG::getNormedG(L2, h);

        rcptr<Factor> ptr1 = uniqptr<SqrtMVG> ( new SqrtMVG({0,1,2,3}, L1, h, g1) );
        rcptr<Factor> ptr2 = uniqptr<SqrtMVG> ( new SqrtMVG({0,1,2,3}, L2, h, g2) );
        ptr1->inplaceCancel(ptr2);
        // and we do nothing, this one threw an exception

      } // if true

      // ******************** Weighted inplace absorb and cancel ********************

      SqrtMVG_InplaceAbsorb iMult;
      SqrtMVG_InplaceCancel iDiv;

      // ++++++++++++++++++++ rhs is subset of lhs ++++++++++++++++++++

      // -------------------- lhs non-vacuous, rhs non-vacuous --------------------

      SqrtMVG tst1(ref1);

      // this one we only probe in a few places
      iMult.weightedInplaceAbsorb(&tst1, &ref2, 0.5);
      diff = 0.0;
      diff = max(diff, abs(ref1K(0,0)+0.5-tst1.getK()(0,0) ) );
      diff = max(diff, abs(ref1K(0,1)-tst1.getK()(0,1) ) );
      diff = max(diff, abs(ref1K(0,2)+0.25-tst1.getK()(0,2) ) );
      diff = max(diff, abs(ref1K(2,2)+1.0-tst1.getK()(2,2) ) );
      diff = max(diff, abs(ref1h[2]+1.0-tst1.getH()[2] ) );
      diff = max(diff, abs(ref1g+ref2g/2.0-tst1.getG() ) );

      if (diff > tol) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_InplaceAbsorb::weightedInplaceAbsorb : deviates by: "
                 << diff << endl << ref1.getK() << endl << ref2.getK() << endl << tst1.getK() );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // and divide it out again
      iDiv.weightedInplaceCancel(&tst1, &ref2, 0.5);

      double tst1g = tst1.getG();
      prlite::ColVector<double> tst1h = tst1.getH();
      prlite::ColMatrix<double> tst1L = tst1.getL();

      diff = 0.0;
      diff = max( diff, abs(ref1g-tst1g) );
      for (int r = 0; r < ref1L.rows(); r++) {
        diff = max( diff, abs(ref1h[r]-tst1h[r]) );
        for (int c = 0; c <= r; c++) {
          diff = max( diff, abs( ref1L(r,c)-tst1L(r,c) ) );
        } // for c
      } // for r

      if ( diff > tol ) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_InplaceCancel::weightedInplaceCancel : deviates by: "
                 << diff << endl << ref1 << endl << tst1 );
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      // ******************** Lets check out the marginalizer ********************

      SqrtMVG_Marginalize marginalizer;

      // ++++++++++++++++++++ retrieve all the dims ++++++++++++++++++++

      // ++++++++++++++++++++ retrieve the initial dims (uncorrelated) ++++++++++++++++++++

      // -------------------- non-vacuous --------------------

      // ref7{0,1,2,3,4} = ref1{0,1,2} * ref6{3,4}

      if (true) {
        newPtr = uniqptr<Factor>( marginalizer.process(&ref7, {0,1,2}) );
        //cout << *newPtr << endl;
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Marginalize::process did not return a SqrtMVG!");
        PRLITE_ASSERT(dwnPtr->getVars() == ref1.getVars(), "vars "
               << dwnPtr->getVars() << " do not match expected " << ref1.getVars());
        dim = dwnPtr->noOfVars();

        diff = 0.0;
        theL = dwnPtr->getL();
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: precision value deviates too far\nExpected: "
                   << ref1L << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        theH = dwnPtr->getH();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref1h[k] - theH[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: h deviates too far\nExpected: "
                   << ref1h << "\nbut got " << theH);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        theG = dwnPtr->getG();
        diff = abs( ref1g - theG );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: g deviates too far\nExpected: "
                   << ref1g << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // -------------------- vacuous --------------------

      // ++++++++++++++++++++ retrieve the final dims (uncorrelated) ++++++++++++++++++++

      // -------------------- non-vacuous --------------------

      // ref7{0,1,2,3,4} = ref1{0,1,2} * ref6{3,4}

      if (true) {

        newPtr = uniqptr<Factor>( marginalizer.process(&ref7, {3,4}) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Marginalize::process did not return a SqrtMVG!");
        PRLITE_ASSERT(dwnPtr->getVars() == ref6.getVars(), "vars "
               << dwnPtr->getVars() << " do not match expected " << ref6.getVars());
        dim = dwnPtr->noOfVars();

        diff = 0.0;
        theL = dwnPtr->getL();
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref6L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: precision value deviates too far\nExpected: "
                   << ref6L << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        theH = dwnPtr->getH();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref6h[k] - theH[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: h deviates too far\nExpected: "
                   << ref6h << "\nbut got " << theH);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        theG = dwnPtr->getG();
        diff = abs( ref6g - theG );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: g deviates too far\nExpected: "
                   << ref6g << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
              } // if

      // -------------------- vacuous --------------------

      // ++++++++++++++++++++ retrieve arb dims (correlated) ++++++++++++++++++++

      // -------------------- non-vacuous --------------------

      // ref8{0,2} is marginal of ref1{0,1,2}

      if (true) {
        newPtr = uniqptr<Factor>( marginalizer.process(&ref1, {0,2}) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Marginalize::process did not return a SqrtMVG!");
        PRLITE_ASSERT(dwnPtr->getVars() == ref8.getVars(), "vars "
               << dwnPtr->getVars() << " do not match expected " << ref8.getVars());
        dim = dwnPtr->noOfVars();

        diff = 0.0;
        theL = dwnPtr->getL();
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref8L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: precision value deviates too far\nExpected: "
                   << ref8L << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        theH = dwnPtr->getH();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref8h[k] - theH[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: h deviates too far\nExpected: "
                   << ref8h << "\nbut got " << theH);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        theG = dwnPtr->getG();
        diff = abs( ref8g - theG );
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: g deviates too far\nExpected: "
                   << ref8g << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if
              } // if

      // -------------------- vacuous --------------------

      // ++++++++++++++++++++ retrieve arb dims (correlated) ++++++++++++++++++++

      // -------------------- non-vacuous --------------------

      //  marginal on {1} of ref1{0,1,2}

      if (true) {
        newPtr = uniqptr<Factor>( marginalizer.process(&ref1, {1}) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Marginalize::process did not return a SqrtMVG!");
        PRLITE_ASSERT(dwnPtr->getVars()[0] == 1, "vars "
               << dwnPtr->getVars() << " do not match expected 1  1");
        dim = dwnPtr->noOfVars();

        double refVal = ref1Cv(1,1);
        double theVal = dwnPtr->getCov()(0,0);
        diff = abs(refVal-theVal);
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: variance value deviates too far\nExpected: "
                   << refVal << "\nbut got " << theVal);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        refVal = ref1Mn[1];
        theVal = dwnPtr->getMean()[0];
        diff = abs(refVal-theVal);
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: mean deviates too far\nExpected: "
                   << refVal << "\nbut got " << theVal);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // -------------------- vacuous --------------------

      // ******************** Lets check out the observer ********************

      SqrtMVG_ObserveAndReduce observer;


      // ++++++++++++++++++++ observe independent variables ++++++++++++++++++++

      // -------------------- non-vacuous --------------------

      // ref7{0,1,2,3,4} = ref1{0,1,2} * ref6{3,4}
      // Observing X3 and X4 should return the mean and L of ref1 unchanged

      if (true) {
        newPtr = uniqptr<Factor>( observer.process(&ref7, {3,4}, {1.1, 2.2}, true) );
        //cout << *newPtr << endl;
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Observe::process did not return a SqrtMVG!");
        PRLITE_ASSERT(dwnPtr->getVars() == ref1.getVars(), "vars "
               << dwnPtr->getVars() << " do not match expected " << ref1.getVars());
        dim = dwnPtr->noOfVars();

        diff = 0.0;
        theL = dwnPtr->getL();
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : L value deviates too far\nExpected: "
                   << ref1L << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        diff = 0.0;
        theMn = dwnPtr->getMean();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : mean value deviates too far\nExpected: "
                   << ref1Mn << "\nbut got " << theMn);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // -------------------- vacuous --------------------

      // ++++++++++++++++++++ observing dependent trailing variables exactly at the mode ++++++++++++++++++++

      // ref1{0,1,2}, observe x[2] = 4.0 i.e. at the mode. Resulting
      // mean should be at same mode, and Lxx is also known

      if (true) {
        newPtr = uniqptr<Factor>( observer.process(&ref1, {2}, {ref1Mn[2]}, true) );
        //cout << *newPtr << endl;
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Observe::process did not return a SqrtMVG!");
        RVIds expectedVars{0,1};
        PRLITE_ASSERT(dwnPtr->getVars() == expectedVars, "vars "
               << dwnPtr->getVars() << " do not match expected " << expectedVars);
        dim = dwnPtr->noOfVars();

        // check L
        diff = 0.0;
        theL = dwnPtr->getL();  // should coincide with top left section of ref1L
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          gLinear::gIndexRange rangeX(0, 1);
          prlite::ColMatrix<double> Lxx; Lxx = ref1L.slice(rangeX, rangeX);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : L value deviates too far\nExpected: "
                   << Lxx << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check mean
        diff = 0.0;
        theMn = dwnPtr->getMean();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : mean value deviates too far\nExpected: "
                   << ref1Mn.subrange(0,1) << "\nbut got " << theMn);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check g
        double expectedG = ref1g + ref1h[2]*ref1Mn[2] - 0.5*ref1K(2,2)*ref1Mn[2]*ref1Mn[2];
        theG = dwnPtr->getG();
        diff = abs(expectedG - theG);
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: g deviates too far\nExpected: "
                   << expectedG << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

      } // if

      // ++++++++++++++++++++ observing dependent non-trailing variables exactly at the mode ++++++++++++++++++++

      // same as previous, but to make the observed variable
      // non-trailing, rename variables to {0 2 1} and observe
      // variable 1

      // ref1{0,1,2}, observe x[2] = 4.0 i.e. at the mode. Resulting
      // mean should be at same mode, and Lxx is also known

      if (true) {
        SqrtMVG* ref1AltPtr = ref1.copy({0,2,1});
        //cout << *ref1AltPtr << endl;
        //cout << ref1AltPtr->getMean() << endl;
        //exit(-1);

        newPtr = uniqptr<Factor>( observer.process(ref1AltPtr, {1}, {ref1Mn[2]}, true) );
        dwnPtr = dynamic_cast<SqrtMVG*>( newPtr.get() );
        PRLITE_ASSERT(dwnPtr, "SqrtMVG_Observe::process did not return a SqrtMVG!");
        RVIds expectedVars{0,2};
        PRLITE_ASSERT(dwnPtr->getVars() == expectedVars, "vars "
               << dwnPtr->getVars() << " do not match expected " << expectedVars);
        dim = dwnPtr->noOfVars();

        // check L
        diff = 0.0;
        theL = dwnPtr->getL();  // should coincide with top left section of ref1L
        for (int r = 0; r < dim; r++) {
          for (int c = 0; c <= r; c++) {
            diff = max<double>( diff, abs( ref1L(r,c) - theL(r,c) ) );
          } // for c
        } // for r

        if (diff > tol) {
          gLinear::gIndexRange rangeX(0, 1);
          prlite::ColMatrix<double> Lxx; Lxx = ref1L.slice(rangeX, rangeX);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : L value deviates too far\nExpected: "
                   << Lxx << "\nbut got " << theL);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check mean
        diff = 0.0;
        theMn = dwnPtr->getMean();
        for (int k = 0; k < dim; k++) {
          diff = max<double>( diff, abs( ref1Mn[k] - theMn[k] ) );
        } // for k

        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Observe : Observing independent components altered the covariance : mean value deviates too far\nExpected: "
                   << ref1Mn.subrange(0,1) << "\nbut got " << theMn);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        // check g
        double expectedG = ref1g + ref1h[2]*ref1Mn[2] - 0.5*ref1K(2,2)*ref1Mn[2]*ref1Mn[2];
        theG = dwnPtr->getG();
        diff = abs(expectedG - theG);
        if (diff > tol) {
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG_Marginalize with rhs subset scope: g deviates too far\nExpected: "
                   << expectedG << " but got " << theG );
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG Failed");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        delete ref1AltPtr;
      } // if

      // #################### END OF TESTS ####################

    } // if

#ifdef TESTING_MODE
    FILEPOS;
#endif

  } // try

  catch (char msg[]) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  // catch (char const* msg) {
  //   LOGERROR(pLog, "Testing", "SqrtMVG threw an exception!");
  //   THROW(prlite::TestFailed, msg);
  // } // catch

  catch (const string& msg) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, msg);
  } // catch

  catch (const exception& e) {
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", "SqrtMVG threw an exception!");
    PRLITE_LOGERROR( PRLITE_pLog, "Testing", e.what() );
    throw(e);
  } // catch

  catch(...) {
    PRLITE_LOGERROR(PRLITE_pLog, "Testing", "SqrtMVG threw an exception!");
    PRLITE_THROW(prlite::TestFailed, "");
  } // catch


  } // runTests
};

PRLITE_REGISTER_TESTCASE(TestSqrtMVG);
