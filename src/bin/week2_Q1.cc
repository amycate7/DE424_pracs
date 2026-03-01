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

//##################################################################
// Some example code. To compile this, go to the emdw/build
// directory and do a:
// cmake ../; make -j7 example
// To run this while in the build directory, do a:
// src/bin/example
//
// For your own stuff, make a copy of this one to start with. Then
// edit the CMakeLists.txt (also in this directory) by adding your
// new target in the same way as this example.
//##################################################################

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
    // Predefine some types and constants
    //*********************************************************

    typedef int T;                  // The type of the values that the RVs can take on
    typedef DiscreteTable<T> DT;    // DT now is a short-hand for DiscreteTable<int>
    double defprob = 0.0;           // Any unspecified probs will default to this.
    rcptr< vector<T> > triDom (     // Domain for Monty-Hall variables: three doors {0,1,2}
        new vector<T>{0,1,2});

    //*********************************************************
    // Define the RVs
    //*********************************************************

    // Define three RV ids: I (initial pick), C (car location), and M (Monty's opened door). 
    enum{I, C, M, R};

    // Define our factors using smart pointers.
    rcptr<Factor> ptrI = 
          uniqptr<DT>(
            new DT(
              {I},
              {triDom}, // Domain of I is the three doors {0,1,2}
              defprob,  // Default probs for unspecified allocations
              {
                {{0}, 1/3.0},
                {{1}, 1/3.0},
                {{2}, 1/3.0},
              } )); // The initial pick is uniform over the three doors
    
    rcptr<Factor> ptrC = 
          uniqptr<DT>(
            new DT(
              {C},
              {triDom},
              defprob,
              {
                {{0}, 1/3.0},
                {{1}, 1/3.0},
                {{2}, 1/3.0},
              } )); // The car location is uniform over the three doors

    // Define the map for M given I and C
    map<vector<T>, FProb> m_given_ic_probs;

    for (int i : {0, 1, 2}) { // Initial pick
      for (int c: {0, 1, 2}) { // Car location
        for (int m: {0, 1, 2}) { // Monty's opened door
          if (m != i && m != c) {
            if (i == c) {
              // You picked the car, so Monty can open either of the two remaining doors with equal probability
              m_given_ic_probs[{m, i, c}] = 0.5;
            } else {
              // You picked a goat, so Monty must open the only remaining door with a goat
              m_given_ic_probs[{m, i, c}] = 1.0;
            }
          }
          // By default, all other cases resolve to 0.0

        }
      }
    }

    rcptr<Factor> ptrMgIC = 
          uniqptr<DT>(
            new DT(
              {M, I, C},
              {triDom, triDom, triDom},
              defprob,
              m_given_ic_probs
            )
          ); 

    // Question 1.2 (c) - Factor p(R|M)
    map<vector<T>, FProb> r_given_m_probs;

    for (int m : {0, 1, 2}) {
      for (int r : {0, 1, 2}) {
        if (r == m) {
          r_given_m_probs[{r, m}] = 0.8; // Monty is truthful 80% of the time
        } else {
          r_given_m_probs[{r,m}] = 0.1; // Monty lies 20% of the time, and if he lies, he randomly opens one of the two remaining doors
        }
      }
    }

    rcptr<Factor> ptrRgM = 
          uniqptr<DT>(
            new DT(
              {R, M},
              {triDom, triDom},
              defprob,
              r_given_m_probs
            )
          );



    // Note the '*' before ptrXY. That is called derefencing
    // Output the defined factors
    std::cout << __FILE__ << __LINE__ << ": " << *ptrI << std::endl; // displays the factor I
    std::cout << __FILE__ << __LINE__ << ": " << *ptrC << std::endl; // displays the factor C
    std::cout << __FILE__ << __LINE__ << ": " << *ptrMgIC << std::endl; // displays the factor MgIC
    std::cout << __FILE__ << __LINE__ << ": " << *ptrRgM << std::endl; // displays the factor RgM

    // Calculate the joint probability P(M, I, C) = P(M|I,C)P(I)P(C)
    // We start with P(M|I,C) and then absorb P(I) and P(C) into it. Absorption is just multiplication.
    rcptr<Factor> ptrJoint = ptrMgIC->absorb(ptrI)->absorb(ptrC);
    std::cout << __FILE__ << __LINE__ << ": " << *ptrJoint << std::endl; // displays the joint factor

    // (b) Determine p(C|I = 0, M = 1) using the observeAndReduce and normalize method
    rcptr<Factor> ptrBelief = ptrJoint->observeAndReduce({I, M}, {0, 1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrBelief << std::endl; // displays the belief factor

    // Question 1.2 (b)
    // Setting I = 0, and the other variables unobserved. 
    // Determine p(C|I = 0). Does the marginal belief over C change when M is unobserved?
    rcptr<Factor> ptrIC = ptrJoint->marginalize({I, C}); // Marginalize out M to get the joint over I and C
    rcptr<Factor> ptrc_given_I = ptrIC->observeAndReduce({I}, {0})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrc_given_I << std::endl; // displays the belief factor

    // Distribution is the same as the prior. Thus, knowing the value of I does not change the marginal belief over C.
    // I and C are independent.
    // If the collider (M) is unobserved, the path between I and C are blocked. They are d-seperated which implied independence.

    // Question 1.2 (b)
    // Setting M = 1, and the other variables unobserved.
    rcptr<Factor> ptrC_given_M = ptrJoint->marginalize({M, C})->observeAndReduce({M}, {1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrC_given_M << std::endl; // displays the belief factor 
    
    // Setting M = 1 and I = 0, and the other variables unobserved.
    rcptr<Factor> ptrC_given_M_and_I = ptrJoint->observeAndReduce({M, I}, {1, 0})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrC_given_M_and_I << std::endl; // displays the belief factor
    
    // The belief over C changes when M and I are observed. 
    // In a collider, observing the middle node (M) opens the path between the two parent nodes (I and C), making them dependent.
    // I and C are not d-separated by M

    // Question 1.2 (c)
    // Joint distribution P(R, M, I, C) = P(R|M)P(M|I,C)P(I)P(C)
    rcptr<Factor> ptrJointFull = ptrRgM->absorb(ptrMgIC)->absorb(ptrI)->absorb(ptrC);
    std::cout << __FILE__ << __LINE__ << ": " << *ptrJointFull << std::endl; // displays the joint distr with factor R

    // Question 1.2 (e)
    // Setting I = 0 , and other variables unobserved. Determine p(C|I = 0). Does the marginal belief over C change when R and M are unobserved?
    rcptr<Factor> ptrIC2 = ptrJointFull->marginalize({I, C}); // Marginalize out M and R to get the joint over I and C
    // Test P(C|I=0) with the full joint
    rcptr<Factor> ptrC_given_I_ = ptrIC2->observeAndReduce({I}, {0})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrC_given_I_ << std::endl; // displays the belief factor
    
    // The belief over C does not change when R and M are unobserved. Without observing R or M, I and C remain independent.

    // Setting R = 1, and other variables unobserved. Determine p(C|R = 1). Does the marginal belief over C change when M and I are unobserved?
    rcptr<Factor> ptrCR = ptrJointFull->marginalize({C, R})->observeAndReduce({R}, {1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrCR << std::endl; // displays the belief factor for p(C|R=1)
    // The belief over C does change when R is observed. Observing R provides information about M, which in turn provides information about I and C.

    // Setting R = 1, and I = 0, and other variables unobserved. Determine p(C|R = 1, I = 0). Does the marginal belief over C change when M is unobserved?
    rcptr<Factor> ptrC_given_R_and_I = ptrJointFull->observeAndReduce({R, I}, {1, 0})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrC_given_R_and_I << std::endl; // displays the belief factor for p(C|R=1, I=0)
    // The belief over C does change when R and I are observed. Observing R provides information about M, which in turn provides information about I and C.
    // order in output table: mc -> probs

    // Conclusion: Observing R does change the marginal belief over C. I and C are dependent given R.
    // This matches d-seperation theory. A path between I and C is unblocked if there is a collider (M) and either the collider or uts descendants are observed.
    // R is a descendant of M, so observing R - even without directly observing M - gives us information about M, which in turn gives us information about I and C, making them dependent.
    // Influence flows from I through M and C.

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
