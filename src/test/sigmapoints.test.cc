/*******************************************************************************

      AUTHOR: JA du Preez
      DATE: May 2014
      PURPOSE: Test code for sigmapoint routines.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

//#ifndef NO_TESTING

// patrec headers
#include "prlite_error.hpp"
#include "prlite_testing.hpp"  // TestCase

// emdw headers
#include "emdw.hpp"
#include "sigmapoints.hpp"

//#include <string>  // string

using namespace std;

class BaseTestSigmaPoints : public prlite::TestCase
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

class TestSigmaPoints : public BaseTestSigmaPoints
{
public:

  /** a name. */
  const string& isA()
  {
    static const string name = "TestSigmaPoints";
    return name;
  }

  /** run the tests. */
  virtual void runTests() {
    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);
#ifdef TESTING_MODE
    FILEPOS;
#endif


    if (true) { // conventional sigma point test

      unsigned dim = 3;
      prlite::ColVector<double> mn1(dim);
      mn1[0] = 1; mn1[1] = -3; mn1[2] = 4;
      prlite::RowMatrix<double> cv1(dim,dim);
      cv1(0,0) = 4; cv1(1,1) = 5; cv1(2,2) = 8;
      cv1(0,1) = cv1(1,0) = 2; cv1(0,2) = cv1(2,0) = -2; cv1(1,2) = cv1(2,1) = -5;
      vector< prlite::ColVector<double> > spts1;
      double w0 = 0.5/(2.0*dim+1);
      vector<double> wgts1;
      cov2sp_2Dp1(mn1, cv1, w0, spts1, wgts1);
      //for (auto elem : spts1) cout << elem << endl;

      prlite::ColVector<double> mn2; prlite::RowMatrix<double> cv2;
      sp2cov(spts1, wgts1, mn2, cv2);

      for (int i = 0; i < int(dim); i++) {

        if ( !almostEqual(mn1[i], mn2[i], 5) ){
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected mean: " << mn1[i] << "; Got: " << mn2[i]);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "EstimateSAMixturePdf_N FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        for (int j = 0; j < int(dim); j++) {
          if ( !almostEqual(cv1(i,j), cv2(i,j), 5) ){
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected covariance: " << cv1(i,j) << "; Got: " << cv2(i,j));
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "EstimateSAMixturePdf_N FAILED!");
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        } // for j
      } // for i

    } // conventional sigma point test

    if (true) { // dual sigma point test

      unsigned dim = 3;
      prlite::ColVector<double> mn1(dim);
      mn1[0] = 1; mn1[1] = -3; mn1[2] = 4;
      prlite::RowMatrix<double> cv1(dim,dim);
      cv1(0,0) = 4; cv1(1,1) = 5; cv1(2,2) = 8;
      cv1(0,1) = cv1(1,0) = 2; cv1(0,2) = cv1(2,0) = -2; cv1(1,2) = cv1(2,1) = -5;
      vector< prlite::ColVector<double> > spts1;
      vector<double> wgts1;
      cov2sp_4Dp1(mn1, cv1, spts1, wgts1);
      //for (auto elem : spts1) cout << elem << endl;

      prlite::ColVector<double> mn2; prlite::RowMatrix<double> cv2;
      sp2cov(spts1, wgts1, mn2, cv2);

      for (int i = 0; i < int(dim); i++) {

        if ( !almostEqual(mn1[i], mn2[i], 5) ){
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected mean: " << mn1[i] << "; Got: " << mn2[i]);
          PRLITE_LOGERROR(PRLITE_pLog, "Testing", "EstimateSAMixturePdf_N FAILED!");
          PRLITE_THROW(prlite::TestFailed, "");
        } // if

        for (int j = 0; j < int(dim); j++) {
          if ( !almostEqual(cv1(i,j), cv2(i,j), 5) ){
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected covariance: " << cv1(i,j) << "; Got: " << cv2(i,j));
            PRLITE_LOGERROR(PRLITE_pLog, "Testing", "EstimateSAMixturePdf_N FAILED!");
            PRLITE_THROW(prlite::TestFailed, "");
          } // if

        } // for j
      } // for i

    } // dual sigma point test


#ifdef TESTING_MODE
    FILEPOS;
#endif

  } // all tests done

};

PRLITE_REGISTER_TESTCASE(TestSigmaPoints);


//#endif // NO_TESTING
