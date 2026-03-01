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
    rcptr< vector<T> > sixDom (     // Domain for die
        new vector<T>{1,2,3,4,5,6});
    rcptr< vector<T> > binDom (     // Domain for binary variables: {0,1} - 0: fair; 1: loaded
        new vector<T>{0,1});

    //*********************************************************
    // Define the RVs
    //*********************************************************

    // Define RV ids for 3 time steps
    enum {z1, z2, z3, y1, y2, y3};

    // Initial state probabilities
    map<vector<T>, FProb> z1_probs;
    z1_probs[{0}] = 0.5; // 50% chance of starting fair
    z1_probs[{1}] = 0.5; // 50% chance of starting loaded

    // Emission probabilities
    map<vector<T>, FProb> emit_probs;

    // Scenario: Zt = 0 (Fair Die)
    for (int y : {1, 2, 3, 4, 5, 6}) {
        emit_probs[{y, 0}] = 1.0 / 6.0; // Uniform distribution
    }

    // Scenario: Zt = 1 (Loaded Die)
    emit_probs[{1, 1}] = 0.1;
    emit_probs[{2, 1}] = 0.1;
    emit_probs[{3, 1}] = 0.1;
    emit_probs[{4, 1}] = 0.1;
    emit_probs[{5, 1}] = 0.1;
    emit_probs[{6, 1}] = 0.5; // Highly biased toward 6

    // Transition Probabilities P(Zt | Zt-1) ---
    map<vector<T>, FProb> trans_probs;
    trans_probs[{0, 0}] = 0.95; // Stay Fair
    trans_probs[{1, 0}] = 0.05; // Switch Fair -> Loaded
    trans_probs[{1, 1}] = 0.90; // Stay Loaded
    trans_probs[{0, 1}] = 0.10; // Switch Loaded -> Fair

    // p(z1)
    rcptr<Factor> ptrZ1 = uniqptr<DT>(new DT({z1}, {binDom}, defprob, z1_probs));
    // p(z2 | z1) and p(z3 | z2)
    rcptr<Factor> ptrZ2gZ1 = uniqptr<DT>(new DT({z2, z1}, {binDom, binDom}, defprob, trans_probs));
    rcptr<Factor> ptrZ3gZ2 = uniqptr<DT>(new DT({z3, z2}, {binDom, binDom}, defprob, trans_probs));

    // p(y1 | z1), p(y2 | z2), p(y3 | z3)
    rcptr<Factor> ptrY1gZ1 = uniqptr<DT>(new DT({y1, z1}, {sixDom, binDom}, defprob, emit_probs));
    rcptr<Factor> ptrY2gZ2 = uniqptr<DT>(new DT({y2, z2}, {sixDom, binDom}, defprob, emit_probs));
    rcptr<Factor> ptrY3gZ3 = uniqptr<DT>(new DT({y3, z3}, {sixDom, binDom}, defprob, emit_probs));

    // 2. Chain the 'absorb' operations to form the full joint (HMM)
    rcptr<Factor> ptrJoint = ptrZ1->absorb(ptrZ2gZ1)
                                  ->absorb(ptrZ3gZ2)
                                  ->absorb(ptrY1gZ1)
                                  ->absorb(ptrY2gZ2)
                                  ->absorb(ptrY3gZ3); 
    // 1728 rows output   

    std::cout << __FILE__ << __LINE__ << ": " << *ptrJoint << std::endl; // displays the joint factor for the HMM

    // Question 3.2 (d)
    // 1. Set the evidence (observations y1=6, y2=5, y3=6)
    // This reduces the 1728-row table to just the 8 hidden state combinations
    rcptr<Factor> ptrPosteriorJoint = ptrJoint->observeAndReduce({y1, y2, y3}, {6, 5, 6});

    // 2. Calculate the marginal belief for each time step
    rcptr<Factor> pZ1_post = ptrPosteriorJoint->marginalize({z1})->normalize();
    rcptr<Factor> pZ2_post = ptrPosteriorJoint->marginalize({z2})->normalize();
    rcptr<Factor> pZ3_post = ptrPosteriorJoint->marginalize({z3})->normalize();

    // 3. Output the results
    cout << "P(z1 | y=6,5,6): " << *pZ1_post << endl;
    cout << "P(z2 | y=6,5,6): " << *pZ2_post << endl;
    cout << "P(z3 | y=6,5,6): " << *pZ3_post << endl;

    // Chain/cascade nodes
    // We are testing the Markov assumptions of the first-order HMM

    // First, let's see if the initial state $z_1$ tells us anything about the future state $z_3$ when we don't know what happened in between.
    rcptr<Factor> ptrZ1Z3 = ptrJoint->marginalize({z1, z3})->normalize();
    // Observe Z1 = 1 (loaded) and see if it changes our belief about Z3
    rcptr<Factor> ptrZ3_given_Z1 = ptrZ1Z3->observeAndReduce({z1}, {1})->normalize();
    cout << "P(z3 | z1=1): " << *ptrZ3_given_Z1 << endl; // We expect this to show a higher probability of z3=1 (loaded) compared to the prior belief

    // Knowing $z_1=1$ significantly increases the probability that $z_3=1$
    // $z_1$ and $z_3$ are Dependent when $z_2$ is unobserved.

    // Fix the mediator state Z2
    // 1. Marginalize to get P(z1, z2, z3)
    rcptr<Factor> ptrZ1Z2Z3 = ptrJoint->marginalize({z1, z2, z3});

    // 2. Observe z2 = 0 (Fair) and z1 = 1 (Loaded)
    rcptr<Factor> ptrZ3_given_Z1Z2 = ptrZ1Z2Z3->observeAndReduce({z1, z2}, {1, 0})->normalize();

    std::cout << "P(z3 | z1=1, z2=0): " << *ptrZ3_given_Z1Z2 << std::endl;

    // Once you know $z_2=0$, the belief about $z_3$ is determined solely by the transition $P(z_3 \mid z_2)$. The fact that $z_1$ was "Loaded" no longer matters.
    // Conclusion: $z_1$ and $z_3$ are Independent given $z_2$.

    // D-seperation theory for a chain. 
    // Path is Open: When $z_2$ is not observed. Influence "flows" from $z_1$ through $z_2$ to $z_3$. This is the basis of Markov chains: the past influences the future.
    // Path is Blocked: When $z_2$ is observed. The middle node "shields" the future from the past. This is known as the Markov Property: the future is independent of the past, given the present.



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
