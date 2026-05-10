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
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "messagequeue.hpp"

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

    // Random seeding
    unsigned seedVal = emdw::randomEngine.getSeedVal();
    cout <<  seedVal << endl;
    emdw::randomEngine.setSeedVal(seedVal);

    // Types and domains
    typedef int T;                  
    typedef DiscreteTable<T> DT; // DT now is a short-hand for DiscreteTable<int>

    double defprob = 0.0; // Any unspecified probs will default to this.
    int T_max = 10; // Number of timesteps

    rcptr< vector<T> > locationDom (     // Domain location of Wumpus: 0 to 24 (5x5 grid)
        new vector<T>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24});
    rcptr< vector<T> > binDom (     // Binary RV: detection D of Wumpus: 0 (no detection), 1 (detection)
        new vector<T>{0,1});
 
    // Define the RVS
    vector<T> W_rvs(T_max); // Vector of ints to store wumpus location for each time step (instead of using enum)
    //vector<vector<RV>> D_rvs(T_max, vector<RV>(25)); // Vector to store detection RVs for each location and timestep

    int num_cells = 25; // 5x5 grid
    for (int t = 0; t < T_max; t++) {
      W_rvs[t] = t; // ID for W_t is simply t
        // for (int loc = 0; loc < num_cells; loc++) {
        //     D_rvs[t][loc] = RV(binDom); // Detection of Wumpus at location loc and time t
        // }
    }

    // Implement the factors
    // Wumpus transition factor: p(W_t | W_{t-1})
    vector<rcptr<Factor>> transitionFactors; // Vector to store our transition factor pointers for every time step

    for (int t = 1; t < T_max; t++) {
      int W_curr = W_rvs[t];
      int W_prev = W_rvs[t-1];

      map<vector<T>, FProb> transitionProbs;
      int width = 5; // Grid width
      
      // i is the index of W_prev
      for (int i = 0; i < num_cells; i++) {
        // Convert cell index to (x, y) coordinates
        int x = i % width; // x-coordinate
        int y = i / width; // y-coordinate
        double stay_prob = 0.0; // Probability of wumpus staying in the same cell

        // Define possible grid moves (up, down, left, right)
        int dx[] = {0, 0, -1, 1};
        int dy[] = {1, -1, 0, 0};

        for (int move = 0; move < 4; move++) {
          int new_x = x + dx[move]; // wumpus new x-coordinate
          int new_y = y + dy[move]; // wumpus new y-coordinate

          if (new_x >= 0 && new_x < width && new_y >=0 && new_y < width) {
            // Valid transition 
            int new_i = new_y * width + new_x; // Convert back to cell index
            transitionProbs[{i, new_i}] = 0.25; 
          } else {
            // Invalid transition
            stay_prob += 0.25; // Accumulate probability of staying in the same cell
          }
        } 
        transitionProbs[{i, i}] = stay_prob;// Probability of wumpus staying in the same cell
      }
      
      rcptr<Factor> ptrTransition = uniqptr<DT>(
        new DT(
          {W_prev, W_curr},           // Variable IDs
          {locationDom, locationDom}, // Their domains
          defprob,
          transitionProbs             // Sparse probability map
        )
      );

      transitionFactors.push_back(ptrTransition);
      cout << "Created transition factor for t=" << t << endl;
      cout << *ptrTransition << endl;

    }
     


    return 0; 
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
