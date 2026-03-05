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

    // Prac 4 Q1 (d) -- Instead of calculating the full joint distribution, we perform a sequence of local marginalizations
    rcptr<Factor> ptr_phi4 = ptrGgDI->observeAndReduce({D}, {1});
    rcptr<Factor> ptr_phi3 = ptrSgI->observeAndReduce({S}, {1});
    rcptr<Factor> ptr_phi2 = ptrI;

    rcptr<Factor> ptr_temp1 =  ptr_phi4->absorb(ptr_phi3)->absorb(ptr_phi2);

    rcptr<Factor> ptr_phi6 = ptr_temp1->marginalize({G});
    rcptr<Factor> ptr_phi5 = ptrLgG;

    rcptr<Factor> ptr_temp2 = ptr_phi6->absorb(ptr_phi5);
    rcptr<Factor> ptr_phi7 = ptr_temp2->marginalize({L})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - VE result: " << *ptr_phi7 << std::endl; 
    
    // The output is consistent with the result obtained by calculating the full joint distribution and then marginalizing out the other variables.
    //  0  0.448909
    //  1  0.551091

    // Prac 4 Q1.2 (d) 
    rcptr<Factor> m_phi2_to_I = ptrI; // The message is the factor itself
    rcptr<Factor> m_phi3_to_I = ptrSgI->observeAndReduce({S}, {1}); 

    rcptr<Factor> m_I_to_phi4 = m_phi2_to_I->absorb(m_phi3_to_I);

    rcptr<Factor> m_phi4_to_G = ptrGgDI->observeAndReduce({D}, {1})->absorb(m_I_to_phi4)->marginalize({G});

    rcptr<Factor> m_phi5_to_L = m_phi4_to_G->absorb(ptrLgG)->marginalize({L})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - BP to determine p(L|D = 1, S = 1): " << *m_phi5_to_L << std::endl;

    // Output
    //  0  0.448909
    //  1  0.551091

    // // Prac 4 Q1.2 (e) 
    // rcptr<Factor> m_L_to_phi5 = new Factor({L});
    // m_L_to_phi5->fill(1.0);
    // rcptr<Factor> temp_phi5 = ptrLgG->absorb(m_phi4_to_G)->absorb(m_L_to_phi5);
    // rcptr<Factor> pL_Result = temp_phi5->marginalize({L})->normalize();
    // std::cout << __FILE__ << ":" << __LINE__ << " - BP to determine p(L|D = 1, S = 1) at factor: " << *pL_Result << std::endl;

    // Question 1.3

    vector< rcptr<Factor> > factorPtrs; // Create a vector of factor pointers
    factorPtrs.push_back(ptrLgG);
    factorPtrs.push_back(ptrGgDI);
    factorPtrs.push_back(ptrSgI);
    factorPtrs.push_back(ptrD);
    factorPtrs.push_back(ptrI);

    // Pass all observed RVs
    map<RVIdType, AnyType> observedRVs;
    observedRVs[D] = 1;
    observedRVs[S] = 1;

    // Create FG
    ClusterGraph cg(ClusterGraph::BETHE, factorPtrs, observedRVs);

    // Init messages and message queue
    map<Idx2, rcptr<Factor>> messages;
    MessageQueue msgQueue;

    // 1. Iterate through all pairs of nodes to find edges
    for (int i = 0; i < cg.nNodes(); ++i) {
        for (int j = 0; j < cg.nNodes(); ++j) {
            if (cg.is_edge(i, j)) {
                Idx2 ij(i, j);
                
                // 2. Initialize the message
                // If create_unit_message is missing, use create_message or 
                // manually construct using the separator variables
                messages[ij] = cg.create_message(ij); 
                messages[ij]->fill(1.0); // Use fill(1.0) or your version's equivalent

                // 3. Try push() or push_back() for the queue
                msgQueue.push(ij); 
            }
        }
    }

    unsigned nMsgs = loopyBP_CG(cg, messages, msgQueue);
    cout << "Sent " << nMsgs << " messages." << endl;

    // p(L|D = 1, S = 1)
    rcptr<Factor> belief_L = queryLBP_CG(cg, messages, {L})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - EMDW result for p(L|D = 1, S = 1): " << *belief_L << std::endl;

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
