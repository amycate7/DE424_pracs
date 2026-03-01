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
    rcptr< vector<T> > triDom (     // Domain for Grade: three grades {1,2,3}
        new vector<T>{1,2,3});
    rcptr< vector<T> > binDom (     // Domain for D, I, S, L: binary variables
        new vector<T>{0,1});

    //*********************************************************
    // Define the RVs
    //*********************************************************

    // Define three RV ids: I (initial pick), C (car location), and M (Monty's opened door). 
    enum{D, I, G, L, S};

    // Define our factors using smart pointers.
    rcptr<Factor> ptrI = 
          uniqptr<DT>(
            new DT(
              {I},
              {binDom},
              defprob,  // Default probs for unspecified allocations
              {
                {{0}, 0.7},
                {{1}, 0.3}
              } )); // Whether student invested in time or not. 70% of students invested time, 30% did not.
    
    rcptr<Factor> ptrD = 
          uniqptr<DT>(
            new DT(
              {D},
              {binDom},
              defprob,
              {
                {{0}, 0.6},
                {{1}, 0.4}
              } ));  // Whether the exam was difficult or not. 60% of the time the exam is difficult, 40% of the time it is not.

    map<vector<T>, FProb> S_given_I_probs;

    for (int i: {0, 1}) {
      for (int s: {0, 1}) {
        if (i == 1 && s == 1) {
          S_given_I_probs[{s, i}] = 0.8; // If the student invested time, they are likely to do well in SAT
        } else if (i == 1 && s == 0) {
          S_given_I_probs[{s, i}] = 0.2;  
        } else if (i == 0 && s == 1) {
          S_given_I_probs[{s, i}] = 0.05; // If the student did not invest time, they are unlikely to do well in SAT
        } else if (i == 0 && s == 0) {
          S_given_I_probs[{s, i}] = 0.95; 
        }
      }
    }

    rcptr<Factor> ptrSgI = 
          uniqptr<DT>(
            new DT(
              {S, I},
              {binDom, binDom},
              defprob,
              S_given_I_probs
            )
          ); 

    map<vector<T>, FProb> L_given_G_probs;

    for (int G : {1, 2, 3}) {
      for (int L : {0, 1}) {
        if (G == 1 && L == 1) {
          L_given_G_probs[{L, G}] = 0.9; // If the student got a grade of 1, they are likely to be given recommendation letter
        } else if (G == 1 && L == 0) {
          L_given_G_probs[{L, G}] = 0.1;  
        } else if (G == 2 && L == 1) {
          L_given_G_probs[{L, G}] = 0.6; // If the student got a grade of 2, they are more likely to be given recommendation letter than not
        } else if (G == 2 && L == 0) {
          L_given_G_probs[{L, G}] = 0.4; 
        } else if (G == 3 && L == 1) {
          L_given_G_probs[{L, G}] = 0.01; // If the student got a grade of 3, they are unlikely to be given recommendation letter
        } else if (G == 3 && L == 0) {
          L_given_G_probs[{L, G}] = 0.99; 
        }
      }
    }

    rcptr<Factor> ptrLgG = 
          uniqptr<DT>(
            new DT(
              {L, G},
              {binDom, triDom},
              defprob,
              L_given_G_probs
            )
          );

    map<vector<T>, FProb> G_given_D_and_I_probs;

    for (int i : {0, 1}) {      // Student's invested time (0: Low, 1: High)
      for (int d : {0, 1}) {    // Difficulty of exam (0: Easy, 1: Hard)
        
        if (i == 0 && d == 0) { 
          // Scenario: Low Effort, Easy Exam
          G_given_D_and_I_probs[{1, 0, 0}] = 0.3; // Grade 1
          G_given_D_and_I_probs[{2, 0, 0}] = 0.4; // Grade 2
          G_given_D_and_I_probs[{3, 0, 0}] = 0.3; // Grade 3
        } 
        else if (i == 0 && d == 1) {
          // Scenario: Low Effort, Hard Exam
          G_given_D_and_I_probs[{1, 1, 0}] = 0.05;
          G_given_D_and_I_probs[{2, 1, 0}] = 0.25;
          G_given_D_and_I_probs[{3, 1, 0}] = 0.70;
        }
        else if (i == 1 && d == 0) {
          // Scenario: High Effort, Easy Exam
          G_given_D_and_I_probs[{1, 0, 1}] = 0.80;
          G_given_D_and_I_probs[{2, 0, 1}] = 0.15;
          G_given_D_and_I_probs[{3, 0, 1}] = 0.05;
        }
        else if (i == 1 && d == 1) {
          // Scenario: High Effort, Hard Exam
          G_given_D_and_I_probs[{1, 1, 1}] = 0.40;
          G_given_D_and_I_probs[{2, 1, 1}] = 0.40;
          G_given_D_and_I_probs[{3, 1, 1}] = 0.20;
        }
      }
    }

    // Construct the factor
    // Variables: {G, D, I}, Domains: {triDom, binDom, binDom}
    rcptr<Factor> ptrGgDI = 
          uniqptr<DT>(
            new DT(
              {G, D, I},
              {triDom, binDom, binDom},
              defprob,
              G_given_D_and_I_probs
            )
          );  


    // Note the '*' before ptrXY. That is called derefencing
    // Output the defined factors
    std::cout << __FILE__ << __LINE__ << ": " << *ptrI << std::endl; // displays the factor I
    std::cout << __FILE__ << __LINE__ << ": " << *ptrD << std::endl; // displays the factor D
    std::cout << __FILE__ << __LINE__ << ": " << *ptrSgI << std::endl; // displays the factor SgI
    std::cout << __FILE__ << __LINE__ << ": " << *ptrLgG << std::endl; // displays the factor LgG
    std::cout << __FILE__ << __LINE__ << ": " << *ptrGgDI << std::endl; // displays the factor GgDI

    // Construct the joint distribution P(G, I, D, L, S) = P(G|D,I)P(L|G)P(S|I)P(D)P(I)
    rcptr<Factor> ptrJoint = ptrGgDI->absorb(ptrLgG)->absorb(ptrSgI)->absorb(ptrD)->absorb(ptrI);
    std::cout << __FILE__ << __LINE__ << ": " << *ptrJoint << std::endl; // displays the joint factor P(G, I, D, L, S)

    // Question 2.1 (b)
    // Determine p(L|D = 1, S = 1)
    rcptr<Factor> ptrLgD_and_S = ptrJoint->marginalize({L, D, S})->observeAndReduce({D, S}, {1, 1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrLgD_and_S << std::endl; // displays the belief factor for p(L|D=1, S=1)

    // It is more likely that the student will receive a good recommendation letter (L=1), since prob =  0.551091

    // Question 2.2 (b)

    // First, marginalize out the other variables (D, L) and the common parent I to see the relationship between G and S.
    // Then, observe G and see if it changes the belief over S
    rcptr<Factor> ptrGS = ptrJoint->marginalize({G, S})->observeAndReduce({S}, {1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrGS << std::endl; // displays the belief factor for p(G|S=1)

    // Next, fix the value of the common parent I to see if it changes the belief over S.
    rcptr<Factor> ptrGSI = ptrJoint->marginalize({G, S, I})->observeAndReduce({S, I}, {1, 1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrGSI << std::endl; // displays the belief factor for p(G|S=1, I=1)

    // The above is the same as if no knowlegdge of S was given (the confounder node blocks the knowlege influence between G and S). Thus, G and S are independent given I.
    rcptr<Factor> ptrGgI = ptrJoint->marginalize({G, I})->observeAndReduce({I}, {1})->normalize();
    std::cout << __FILE__ << __LINE__ << ": " << *ptrGgI << std::endl; // displays the belief factor for p(G|I=1)



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
