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
#include <fstream> // For file reading
#include <iomanip> // For formatting file names
#include <filesystem>

using namespace std;
using namespace emdw;
namespace fs = std::filesystem;

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

    //======================================
    // Types and domains
    // =====================================
    typedef int T;                  
    typedef DiscreteTable<T> DT; // DT now is a short-hand for DiscreteTable<int>

    double defprob = 0.0; // Any unspecified probs will default to this.
    int T_max = 20; // Number of timesteps    
    int num_cells = 200; // 10x20 grid
    int R = 10; // Num rows
    int C = 20; // Num columns
    
    // Variable initialisation
    double pw = 0.7; // Initial estimates of pw
    double pc = 0.1; // Initial estimates of pc
    double epsilon = 1e-5;
    double delta = 1.0;
    int max_iters = 30;

    // Domain location of Wumpus: 0 to 200 (10x20 grid) 
    rcptr< vector<T> > locationDom(new vector<T>);
    for (int i = 0; i < num_cells; i++) {
      locationDom->push_back(i);
    }

    rcptr< vector<T> > binDom ( // Binary RV: detection D of Wumpus: 0 (no detection), 1 (detection)
      new vector<T>{0,1});

    // Create MAP operators 
    rcptr<Marginalizer> margPtr = uniqptr<Marginalizer>(new DiscreteTable_MaxMarginalize<T>);
    rcptr<InplaceNormalizer> iNormPtr = uniqptr<InplaceNormalizer>(new DiscreteTable_InplaceMaxNormalize<T>);
    rcptr<Normalizer> normPtr = uniqptr<Normalizer>(new DiscreteTable_MaxNormalize<T>);
    double margin = 0.0;
    double floor = 0.0;

    // EM requires Sum-Product, not Max-Product
    rcptr<Marginalizer> sumMarg = uniqptr<Marginalizer>(new DiscreteTable_SumMarginalize<T>);
    rcptr<Normalizer> sumNorm = uniqptr<Normalizer>(new DiscreteTable_SumNormalize<T>);
    rcptr<InplaceNormalizer> sumINorm = uniqptr<InplaceNormalizer>(new DiscreteTable_InplaceSumNormalize<T>);

    // =====================================
    // Define the RVS
    // =====================================
    vector<T> W_rvs(T_max); // Vector of ints to store wumpus location for each time step (instead of using enum)
    vector<vector<T>> D_rvs(T_max, vector<T>(200)); // Vector to store detection RVs for each location and timestep
    
    for (int t = 0; t < T_max; t++) {
      W_rvs[t] = t; // ID for W_t is simply t
      for (int loc = 0; loc < num_cells; loc++) {
        D_rvs[t][loc] = 20 + (t * 200) + loc;
      }
    }

    // ========================================================
    // Data is preloaded so that we don't reload each iteration
    // ========================================================
    map<RVIdType, AnyType> obsv;

    for (int t = 0; t < T_max; t++) {
      // First construct the file name (data_file000.txt etc...)
      stringstream ss;
      ss << "../src/dataset3/data_file" 
         << setfill('0') << setw(3) << t << ".txt"; // Ensures three digit padding (fills with zeroes)
      string filename = ss.str();
 
      ifstream dataFile(filename);

      if (!dataFile.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        continue;
      }

      // Read detection values from the file for timestep t (location tracks the cell)
      for (int loc = 0; loc < num_cells; loc++) {
        int detectionVal;
        // Attempt to read file into variable and check if the operation is successful
        if (dataFile >> detectionVal) {
          int rvID = D_rvs[t][loc]; // Obtain RV id
          obsv[rvID] = int(detectionVal); // Ensure its an int by casting
        }
      }

      dataFile.close();
      cout << "Successfully loaded evidence from " << filename << endl;
    }

    // ====================================
    // Implement the factors
    // ====================================

    // Wumpus transition factor: p(W_t | W_{t-1}): Only create once since this is not dependent on unknown parameters
    vector<rcptr<Factor>> transitionFactors; // Vector to store our transition factor pointers for every time step

    for (int t = 1; t < T_max; t++) {
      int W_curr = W_rvs[t];
      int W_prev = W_rvs[t-1];

      map<vector<T>, FProb> transitionProbs;
      
      // i is the index of W_prev
      for (int i = 0; i < num_cells; i++) {
        // Convert cell index to (x, y) coordinates
        int x = i % C; // Column index (0 to 19)
        int y = i / C; // Row index (0 to 9)
        double stay_prob = 0.0; // Probability of wumpus staying in the same cell

        // Define possible grid moves (up, down, left, right)
        int dx[] = {0, 0, -1, 1};
        int dy[] = {1, -1, 0, 0};

        for (int move = 0; move < 4; move++) {
          int new_x = x + dx[move]; // wumpus new x-coordinate
          int new_y = y + dy[move]; // wumpus new y-coordinate

          if (new_x >= 0 && new_x < C && new_y >=0 && new_y < R) {
            // Valid transition 
            int new_i = new_y * C + new_x; // Convert back to cell index
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
          transitionProbs,            // Sparse probability map
          margin, floor, false,       // Extra arguments for MAP    
          sumMarg, sumINorm, sumNorm   
        )
      );

      transitionFactors.push_back(ptrTransition);
      //cout << "Created transition factor for t=" << t << endl;
      //cout << *ptrTransition << endl;
    }

    // Build helper detection factor for MAX
    void buildDetectionFactorMax(double current_pw, double current_pc) {

    }

    // ============================
    // EM Algorithm
    // ============================
    for (int iter = 0; iter < max_iters && delta > epsilon; iter++) {
      double old_pw = pw;
      double old_pc = pc;
      vector<rcptr<Factor>> factorPtrs;

      // Push transition factors onto factorPtrs
      for (auto& ptr : transitionFactors) {
        factorPtrs.push_back(ptr);
      }

      // E-step: Construct detection factors based on current param estimates and perform inference
      // Wumpus Detection factors p(D_loc^t|W^t)
      for (int t = 0; t < T_max; t++) {
        int W_curr = W_rvs[t]; // Extract RV index of W RV at time t

        // loc is the location of the specific cell we are building a factor for
        for (int loc = 0; loc < num_cells; loc++) {
          int D_curr = D_rvs[t][loc]; // Extract index of D RV at time t and grid location loc
          map<vector<T>, FProb> detectionProbs;

          // w_pos is the latent position of the wumpus 
          for (int w_pos = 0; w_pos < num_cells; w_pos++){
            if (w_pos == loc) {
              // Detection location is equal to the wumpus position
              detectionProbs[{1, w_pos}] = current_pw; // Detection
              detectionProbs[{0, w_pos}] = 1.0 - current_pw;  // Missed detection
            } else {
              detectionProbs[{1, w_pos}] = current_pc; // Clutter 
              detectionProbs[{0, w_pos}] = 1.0 - current_pc; // No detection
            }
          } // end of wumpus location loop
          
          rcptr<Factor> ptrDetections = uniqptr<DT>(
            new DT(
              {D_curr, W_curr},
              {binDom, locationDom},
              defprob,
              detectionProbs,
              margin, floor, false,     // Ensure BP
              sumMarg, sumINorm, sumNorm
            )
          );

          factorPtrs.push_back(ptrDetections);

        } // end of grid location loop
      } // end of outer time step loop

      //Construct CG
      ClusterGraph cg(ClusterGraph::BETHE, factorPtrs, obsv);
      map<Idx2, rcptr<Factor> > msgs;
      MessageQueue msgQ;
      unsigned nMsgs = loopyBP_CG(cg, msgs, msgQ); // CG is calibrated after this pass

      // M-step: param updates
      double pw_num = 0.0;
      double pc_num = 0.0;
      for (int t = 0; t < T_max; t++) {
        rcptr<Factor> beliefT = queryLBP_CG(cg, msgs, {W_rvs[t]})->normalize(); // Query the graph
        int St = 0; // total detections at time t
        for (int loc = 0; loc < num_cells; loc++) {
          if ((int)obsv[D_rvs[t][loc]] == 1) St++;
        }

        for (int i = 0; i < num_cells; i++) {
            double gamma_ti = beliefT->potentialAt({W_rvs[t]}, {(T)i});
            int Y_ti = ((int)obsv[D_rvs[t][i]] == 1);
            if (Y_ti) pw_num += gamma_ti;
            pc_num += gamma_ti * (St - Y_ti);
        }
      }
      pw = pw_num / T_max; // Updated pw
      pc = pc_num / (T_max * (num_cells - 1)); // Updated pc
      delta = abs(pw - old_pw) + abs(pc - old_pc);
      cout << "Iteration " << iter << ": pw=" << pw << ", pc=" << pc << endl;
    }  

    // After params have converged, run one more pass with MAP inference to obtain the trajectory
    // Create a vector container to hold all the factors of the model
    vector<rcptr<Factor>> finalPtrs;

    // Rebuild transition factor with margPtr (MAX)
    vector<rcptr<Factor>> transitionFactors; // Vector to store our transition factor pointers

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
          transitionProbs,            // Sparse probability map
          margin, floor, false,       // Extra arguments for MAP    
          margPtr, iNormPtr, normPtr   
        )
      );

      transitionFactors.push_back(ptrTransition);
      //cout << "Created transition factor for t=" << t << endl;
      //cout << *ptrTransition << endl;
    }

    // Push transition factors onto factorPtrs
    for (auto& ptr : transitionFactors) {
      finalPtrs.push_back(ptr);
    }
    
    // Build Detection factors with MAX
    vector<rcptr<Factor>> detectionFactors; // Vector of pointers pointing to detection factors for each time step
    for (int t = 0; t < T_max; t++) {
      int W_curr = W_rvs[t]; // Extract RV index of W RV at time t

      // loc is the location of the specific cell we are building a factor for
      for (int loc = 0; loc < num_cells; loc++) {
        int D_curr = D_rvs[t][loc]; // Extract index of D RV at time t and grid location loc
        map<vector<T>, FProb> detectionProbs;

        // w_pos is the latent position of the wumpus 
        for (int w_pos = 0; w_pos < num_cells; w_pos++){
          if (w_pos == loc) {
            // Detection location is equal to the wumpus position
            detectionProbs[{1, w_pos}] = pw; // Detection
            detectionProbs[{0, w_pos}] = 1.0 - pw;  // Missed detection
          } else {
            detectionProbs[{1, w_pos}] = pc; // Clutter 
            detectionProbs[{0, w_pos}] = 1.0 - pc; // No detection
          }
        } // end of wumpus location loop
          
        rcptr<Factor> ptrDetections = uniqptr<DT>(
          new DT(
            {D_curr, W_curr},
            {binDom, locationDom},
            defprob,
            detectionProbs,
            margin, floor, false,     // Extra arguments for MAP 
            margPtr, iNormPtr, normPtr
          )
        );

        detectionFactors.push_back(ptrDetections);

      } // end of grid location loop
    } // end of outer time step loop

    for (auto& ptr : detectionFactors) {
      finalPtrs.push_back(ptr);
    }

    ClusterGraph cgFinal(ClusterGraph::BETHE, finalPtrs, obsv);
    map<Idx2, rcptr<Factor>> finalMsgs;
    MessageQueue finalQ;
    unsigned nMsgs = loopyBP_CG(cgFinal, finalMsgs, finalQ);

    // ================================================
    // Inferred wumpus trajectory using MAP inference
    // ================================================
        // Save MAP results to a .txt for further model performance evaluation
    ofstream outFile("wumpus_location3.txt");

    if (!outFile.is_open()) {
      cerr << "Error: Could not create wumpus_location3.txt" << endl;
    } else {
      cout << "Inferred MAP trajectory saved to file wumpus_location3.txt" << endl;
      
      for (int t = 0; t < T_max; t++) {
        rcptr<Factor> beliefT = queryLBP_CG(cgFinal, msgs, {W_rvs[t]})->normalize(); // Query the graph

        double max_likelihood = -1.0;
        int best_cell = -1;

        for (int row = 0; row < R; row++) {
          for (int col = 0; col < C; col++) {
            // Calculate the cell index and extract likelihood
            unsigned int cell_idx = row * C + col;
            double likelihood = beliefT->potentialAt({W_rvs[t]}, {(T)cell_idx});

            if (likelihood > max_likelihood) {
              max_likelihood = likelihood;
              best_cell = cell_idx;
            }

          }
        }

        // Convert cell index back to x and y coordinates
        int x = best_cell % C;
        int y = best_cell / C;

        outFile << x << " " << y << endl;

      // Convert the bestCell index back to (row, col) for clear output
      cout << "Time " << t << ": Wumpus at [" << x << ", " << y 
          << "] (Cell ID: " << best_cell << ", Likelihood: " << max_likelihood << ")" << endl;

      }
    }





    } // end of EM

    // =============================
    // Construct cluster graph (should result in junction tree, therefore exact inference)
    // =============================


    // =========================================================================================================
    // Perform MAP inference on Wumpus location for all time steps and save results to file wumpus_location1.txt
    // =========================================================================================================

    // Save MAP results to a .txt for further model performance evaluation
    ofstream outFile("wumpus_location3.txt");

    if (!outFile.is_open()) {
      cerr << "Error: Could not create wumpus_location3.txt" << endl;
    } else {
      cout << "Inferred MAP trajectory saved to file wumpus_location3.txt" << endl;
      
      for (int t = 0; t < T_max; t++) {
        rcptr<Factor> beliefT = queryLBP_CG(cg, msgs, {W_rvs[t]})->normalize(); // Query the graph

        double max_likelihood = -1.0;
        int best_cell = -1;

        for (int row = 0; row < R; row++) {
          for (int col = 0; col < C; col++) {
            // Calculate the cell index and extract likelihood
            unsigned int cell_idx = row * C + col;
            double likelihood = beliefT->potentialAt({W_rvs[t]}, {(T)cell_idx});

            if (likelihood > max_likelihood) {
              max_likelihood = likelihood;
              best_cell = cell_idx;
            }

          }
        }

        // Convert cell index back to x and y coordinates
        int x = best_cell % C;
        int y = best_cell / C;

        outFile << x << " " << y << endl;

      // Convert the bestCell index back to (row, col) for clear output
      cout << "Time " << t << ": Wumpus at [" << x << ", " << y 
          << "] (Cell ID: " << best_cell << ", Likelihood: " << max_likelihood << ")" << endl;

      }
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
