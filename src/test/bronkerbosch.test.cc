/*******************************************************************************

      AUTHOR: JA du Preez
      DATE: May 2014
      PURPOSE: Test code for sigmapoint routines.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

//#ifndef NO_TESTING

// patrec headers
#include "prlite_logging.hpp"
#include "prlite_testing.hpp"  // TestCase

// emdw headers
#include "bronkerbosch.hpp"

//#include <string>  // string

using namespace std;
using namespace emdw;

class BaseTestBronKerbosch : public prlite::TestCase {
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

class TestBronKerbosch : public BaseTestBronKerbosch {
public:

  /** a name. */
  const string& isA() {
    static const string name = "TestBronKerbosch";
    return name;
  }

  /** run the tests. */
  virtual void runTests() {

    PRLITE_LOGDEBUG1(PRLITE_pLog, "Testing", "Testing: " << __FILE__ << ':' << __LINE__);

#ifdef TESTING_MODE
    FILEPOS;
#endif

    if (true) { // extending the three hamming parity check cliques

      map<RVIdType, set<RVIdType> > adj;

      set<RVIdType> clique0 {0,1,2,3}; // implicitly there
      set<RVIdType> clique1 {0,1,2,4}; // parity check 1
      set<RVIdType> clique2 {0,1,3,6}; // parity check 2
      set<RVIdType> clique3 {0,2,3,5}; // parity check 3

      // adjacencies extracted from cliques 1..3
      adj[0] = set<RVIdType> {1,2,4,3,5,6};
      adj[1] = set<RVIdType> {0,2,4,3,6};
      adj[2] = set<RVIdType> {0,1,4,3,5};
      adj[3] = set<RVIdType> {0,2,5,1,6};
      adj[4] = set<RVIdType> {0,1,2};
      adj[5] = set<RVIdType> {0,2,3};
      adj[6] = set<RVIdType> {0,1,3};

      set<RVIdType> pSet {0,1,2,3,4,5,6};
      set<RVIdType> rSet;
      set<RVIdType> xSet;
      vector< set<RVIdType> > maxCliques;
      BronKerbosch2(adj, rSet, pSet, xSet, maxCliques);

      if (maxCliques.size() != 4) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Expected to find 4 cliques, but got "
                 << maxCliques.size());
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "BronKerbosch FAILED!");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      if (maxCliques[0] != clique0) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Did not detect clique 0");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "BronKerbosch FAILED!");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      if (maxCliques[1] != clique1) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Did not detect clique 1");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "BronKerbosch FAILED!");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      if (maxCliques[2] != clique2) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Did not detect clique 2");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "BronKerbosch FAILED!");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

      if (maxCliques[3] != clique3) {
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "Did not detect clique 3");
        PRLITE_LOGERROR(PRLITE_pLog, "Testing", "BronKerbosch FAILED!");
        PRLITE_THROW(prlite::TestFailed, "");
      } // if

    } // extending the three hamming parity check cliques

#ifdef TESTING_MODE
    FILEPOS;
#endif

  } // all tests done

};

PRLITE_REGISTER_TESTCASE(TestBronKerbosch);


//#endif // NO_TESTING
