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
    rcptr< vector<T> > bindom (     // Domain for binary variables
        new vector<T>{0,1});

    //*********************************************************
    // Define the RVs
    //*********************************************************

    // Define RVs 
    enum{b0, b1, b2, b3, b4, b5, b6, r0, r1, r2, r3, r4, r5, r6};

    rcptr<Factor> phi_c1 = 
          uniqptr<DT>(
            new DT(
              {b0,b1,b2,b4},
              {bindom,bindom,bindom,bindom},
              defprob,
              {
                {{0,0,0,0}, 1},
                {{1,0,0,1}, 1},
                {{0,1,0,1}, 1},
                {{0,0,1,1}, 1},
                {{1,1,0,0}, 1},
                {{0,1,1,0}, 1},
                {{1,0,1,0}, 1},
                {{1,1,1,1}, 1}
              })
          );

    rcptr<Factor> phi_c2 =
          uniqptr<DT>(
            new DT(
              {b0,b2,b3,b5},
              {bindom,bindom,bindom,bindom},
              defprob,
              {
                {{0,0,0,0}, 1},
                {{1,0,0,1}, 1},
                {{0,1,0,1}, 1},
                {{0,0,1,1}, 1},
                {{1,1,0,0}, 1},
                {{0,1,1,0}, 1},
                {{1,0,1,0}, 1},
                {{1,1,1,1}, 1}
              }
            )
          );

    rcptr<Factor> phi_c3 =
          uniqptr<DT>(
            new DT(
              {b0,b1,b3,b6},
              {bindom,bindom,bindom,bindom},
              defprob,
              {
                {{0,0,0,0}, 1},
                {{1,0,0,1}, 1},
                {{0,1,0,1}, 1},
                {{0,0,1,1}, 1},
                {{1,1,0,0}, 1},
                {{0,1,1,0}, 1},
                {{1,0,1,0}, 1},
                {{1,1,1,1}, 1}
              }
            )
          );

    rcptr<Factor> ptrR0gB0 = 
          uniqptr<DT>(
            new DT(
              {b0,r0},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );

    rcptr<Factor> ptrR1gB1 = 
          uniqptr<DT>(
            new DT(
              {b1,r1},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );

    rcptr<Factor> ptrR2gB2 = 
          uniqptr<DT>(
            new DT(
              {b2,r2},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );
     
    rcptr<Factor> ptrR3gB3 = 
          uniqptr<DT>(
            new DT(
              {b3,r3},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );    

    rcptr<Factor> ptrR4gB4 = 
          uniqptr<DT>(
            new DT(
              {b4,r4},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );   

    rcptr<Factor> ptrR5gB5 = 
          uniqptr<DT>(
            new DT(
              {b5,r5},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );

    rcptr<Factor> ptrR6gB6 = 
          uniqptr<DT>(
            new DT(
              {b6,r6},
              {bindom, bindom},
              defprob,
              {
                {{0,0}, 0.9},
                {{0,1}, 0.1},
                {{1,0}, 0.1},
                {{1,1}, 0.9}
              }
            )
          );



    // (g)
    rcptr<Factor> ptrJoint = phi_c1->absorb(phi_c2)->absorb(phi_c3)->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Joint factor: " << *ptrJoint << std::endl;

    //(h)
    rcptr<Factor> ptrBits = ptrJoint->observeAndReduce({b0,b1,b2,b3}, {1,0,1,0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Observed bits: " << *ptrBits << std::endl; 

    // Demodulator receives degraded version of the transmitted waveforms
    // (g) Create an error in r1. Transmitted: 1010 001. Received: 1110 001(parity)
    rcptr<Factor> ptrR0 = ptrR0gB0->observeAndReduce({r0}, {1})->normalize();
    rcptr<Factor> ptrR1 = ptrR1gB1->observeAndReduce({r1}, {1})->normalize();
    rcptr<Factor> ptrR2 = ptrR2gB2->observeAndReduce({r2}, {1})->normalize();
    rcptr<Factor> ptrR3 = ptrR3gB3->observeAndReduce({r3}, {0})->normalize();
    rcptr<Factor> ptrR4 = ptrR4gB4->observeAndReduce({r4}, {0})->normalize();
    rcptr<Factor> ptrR5 = ptrR5gB5->observeAndReduce({r5}, {0})->normalize();
    rcptr<Factor> ptrR6 = ptrR6gB6->observeAndReduce({r6}, {1})->normalize();

    // (g) p(b0, . . . , b6|r0, . . . , r6 = the observed values).
    rcptr<Factor> ptrReceivedJoint = ptrJoint->absorb(ptrR0)->absorb(ptrR1)->absorb(ptrR2)->absorb(ptrR3)->absorb(ptrR4)->absorb(ptrR5)->absorb(ptrR6)->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Decoded sequence: " << *ptrReceivedJoint << std::endl; 

    rcptr<Factor> ptrb0gR = ptrReceivedJoint->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Marginal belief p(b0 = 1|observed values): " << *ptrReceivedJoint << std::endl; 


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
