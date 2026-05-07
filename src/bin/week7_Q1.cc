/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"

// emdw headers
#include "emdw.hpp"
#include "discretetable.hpp"
#include "prlite_genmat.hpp"
#include "sqrtmvg.hpp"

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

typedef SqrtMVG SG; // SG now is a short-hand for SqrtMVG


int main(int, char *argv[]) {

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {

    //*********************************************************
    // Some random generator seeding. Just keep this as is
    //*********************************************************

    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout <<  seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    //*********************************************************

    // Factor 1
    prlite::ColVector<double> mn1(2);
    mn1[0] = 1 ; mn1[1] = 2;

    prlite::RowMatrix<double> cov1(2,2);
    cov1(0,1) = cov1(1,0) = 2;
    cov1(0,0) = 1;
    cov1(1,1) = 5;

    rcptr<Factor> ptr1 (uniqptr<SG>(new SG({1,2}, mn1, cov1))); // x1 identifier = 1; x2 identifier = 2

    // Factor 2
    prlite::ColVector<double> mn2(2);
    mn2[0] = 3 ; mn2[1] = 4;

    prlite::RowMatrix<double> cov2(2,2);
    cov2(0,1) = cov2(1,0) = 2;
    cov2(0,0) = 2;
    cov2(1,1) = 3;

    rcptr<Factor> ptr2 (uniqptr<SG>(new SG({2,3}, mn2, cov2))); // x2 identifier = 2; x3 identifier = 3

    rcptr<Factor> ptr3 = (*ptr1) * (*ptr2);
    rcptr<SG> dwnPtr; dwnPtr = dynamic_pointer_cast<SG>(ptr3); // Downcast to SqrtMVG to access mean and covariance

    cout << "Information matrix: " << dwnPtr->getK() << endl;
    cout << "Information vector: " << dwnPtr->getH() << endl;
    cout << "Covariance matrix: " << dwnPtr->getCov() << endl;

    // Since no entries of the covariance matrix are zero, therefore random variabes are statistically dependent. 
    // $x_1 \perp x_3 \mid x_2$.
    // MN Structure conclusion on conditional independence: Yes, because $x_2$ is the only path between $x_1$ and $x_3$ (it's a "separator"). Observing $x_2$ blocks the flow of influence.





    return 0; // tell the world that all is fine
  } // try

  catch (char msg[]) {
    cerr << msg << endl;
  } // catch

  // catch (char const* msg) {
  //   cerr << msg << endl;
  // } // catch

  catch (const string& msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const exception& e) {
    cerr << "Unhandled exception: " << e.what() << endl;
    throw e;
  } // catch

  catch(...) {
    cerr << "An unknown exception / error occurred\n";
    throw;
  } // catch

} // main
