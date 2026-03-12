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



  // Prac 3 code
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
    std::cout << __FILE__ << ":" << __LINE__ << " - Marginal belief p(b0 = 1|observed values): " << *ptrb0gR << std::endl; 

    // Prac 5 code
    // Observe and reduce the following sequence: 1110 001. The true sequence of bits is: 1010 001
    // 1.1 c.) Construct cluster potentials
    rcptr<Factor> ptrCluster0 = ptrR0->absorb(ptrR1)->absorb(ptrR2)->absorb(ptrR3); // Cluster 0 contains all the observed factors
    rcptr<Factor> ptrCluster1 = phi_c1->absorb(ptrR4);
    rcptr<Factor> ptrCluster2 = phi_c2->absorb(ptrR5);
    rcptr<Factor> ptrCluster3 = phi_c3->absorb(ptrR6);

    // 1.1 d.) Construct messages for JT
    // Inward pass
    rcptr<Factor> m1_0 = ptrCluster1->marginalize({b0,b1,b2})->normalize(); // Normalize to prevent underflow
    rcptr<Factor> m2_0 = ptrCluster2->marginalize({b0,b2,b3})->normalize();
    rcptr<Factor> m3_0 = ptrCluster3->marginalize({b0,b1,b3})->normalize();

    // Outward pass
    rcptr<Factor> m0_1 = ptrCluster0->absorb(m2_0)->absorb(m3_0)->marginalize({b0,b1,b2})->normalize();
    rcptr<Factor> m0_2 = ptrCluster0->absorb(m1_0)->absorb(m3_0)->marginalize({b0,b2,b3})->normalize();
    rcptr<Factor> m0_3 = ptrCluster0->absorb(m1_0)->absorb(m2_0)->marginalize({b0,b1,b3})->normalize();

    // 1.1 e.) Calculate the belief over b0
    rcptr<Factor> ptrC1_belief_b0 = ptrCluster1->absorb(m0_1)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from Cluster 1: " << *ptrC1_belief_b0 << std::endl;

    rcptr<Factor> ptrC2_belief_b0 = ptrCluster2->absorb(m0_2)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from Cluster 2: " << *ptrC2_belief_b0 << std::endl;

    rcptr<Factor> ptrSep03_belief_b0 = m3_0->absorb(m0_3)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from separator 0-3: " << *ptrSep03_belief_b0 << std::endl;
    // The beliefs are consistent across the clusters and separator -- calibrated

    //1.1 f.) Determine beliefs over all transmitted bits
    rcptr<Factor> beliefB1 = m1_0->absorb(m0_1)->marginalize({b1})->normalize(); // Must always normalize!
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b1: " << *beliefB1 << std::endl;

    rcptr<Factor> beliefB2 = m1_0->absorb(m0_1)->marginalize({b2})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b2: " << *beliefB2 << std::endl;

    rcptr<Factor> beliefB3 = m2_0->absorb(m0_2)->marginalize({b3})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b3: " << *beliefB3 << std::endl;

    rcptr<Factor> beliefB4 = ptrCluster1->absorb(m0_1)->marginalize({b4})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b4: " << *beliefB4 << std::endl;

    rcptr<Factor> beliefB5 = ptrCluster2->absorb(m0_2)->marginalize({b5})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b5: " << *beliefB5 << std::endl;  

    rcptr<Factor> beliefB6 = ptrCluster3->absorb(m0_3)->marginalize({b6})->normalize(); 
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b6: " << *beliefB6 << std::endl;

    // Question 1.2
    // 1.2 c.) Received bits: 1110 001. The true sequence of bits is: 1010 001
    rcptr<Factor> cl0 = phi_c1->absorb(ptrR0)->absorb(ptrR1)->absorb(ptrR2)->absorb(ptrR4);
    rcptr<Factor> cl1 = phi_c2->absorb(ptrR5);
    rcptr<Factor> cl2 = phi_c3->absorb(ptrR3)->absorb(ptrR6);

    // 1.2 d.) Initialise all CG messages to uniform distributions. The scope of each message must match the sepset of the connected clusters. 
    rcptr<Factor> m_cl0_to_cl1 = uniqptr<DT>(new DT({b0, b2}, {bindom, bindom}, 1.0));
    rcptr<Factor> m_cl1_to_cl0 = uniqptr<DT>(new DT({b0, b2}, {bindom, bindom}, 1.0));

    rcptr<Factor> m_cl0_to_cl2 = uniqptr<DT>(new DT({b0, b1}, {bindom, bindom}, 1.0));
    rcptr<Factor> m_cl2_to_cl0 = uniqptr<DT>(new DT({b0, b1}, {bindom, bindom}, 1.0));

    rcptr<Factor> m_cl1_to_cl2 = uniqptr<DT>(new DT({b3}, {bindom}, 1.0));
    rcptr<Factor> m_cl2_to_cl1 = uniqptr<DT>(new DT({b3}, {bindom}, 1.0));

    // 1.2 e.) Calculate all messages
    // 1. Update Message: Cluster 0 -> Cluster 1 (Sepset {b0, b2})
    m_cl0_to_cl1 = cl0->absorb(m_cl2_to_cl0)->marginalize({b0,b2})->normalize();

    // 2. Update Message: Cluster 1 -> Cluster 2 (Sepset {b3})
    m_cl1_to_cl2 = cl1->absorb(m_cl0_to_cl1)->marginalize({b3})->normalize();

    // 3. Update Message: Cluster 2 -> Cluster 0 (Sepset {b0, b1})
    m_cl2_to_cl0 = cl2->absorb(m_cl1_to_cl2)->marginalize({b0,b1})->normalize();

    // Reverse direction
    // 4. Update Message: Cluster 0 -> Cluster 2 (Sepset {b0, b1})
    m_cl0_to_cl2 = cl0->absorb(m_cl1_to_cl0)->marginalize({b0,b1})->normalize();

    // 5. Update Message: Cluster 2 -> Cluster 1 (Sepset {b3})
    m_cl2_to_cl1 = cl2->absorb(m_cl0_to_cl2)->marginalize({b3})->normalize();

    // 6. Update Message: Cluster 1 -> Cluster 0 (Sepset {b0, b2})
    m_cl1_to_cl0 = cl1->absorb(m_cl2_to_cl1)->marginalize({b0,b2})->normalize();

    // 1.2 f.) Iterate message-passing until convergence. 
    bool converged = false; // convergence flag
    int iteration = 0; // iteration count
    double threshold = 1e-6; // convergence threshold

    rcptr<Factor> bel01_old, bel02_old, bel12_old; // Store old sepset beliefs

    while (!converged && iteration < 100) {
      iteration++;
      // Update messages
      m_cl0_to_cl1 = cl0->absorb(m_cl2_to_cl0)->marginalize({b0,b2})->normalize();
      m_cl1_to_cl2 = cl1->absorb(m_cl0_to_cl1)->marginalize({b3})->normalize();
      m_cl2_to_cl0 = cl2->absorb(m_cl1_to_cl2)->marginalize({b0,b1})->normalize();
      m_cl0_to_cl2 = cl0->absorb(m_cl1_to_cl0)->marginalize({b0,b1})->normalize();
      m_cl2_to_cl1 = cl2->absorb(m_cl0_to_cl2)->marginalize({b3})->normalize();
      m_cl1_to_cl0 = cl1->absorb(m_cl2_to_cl1)->marginalize({b0,b2})->normalize();

      // Calculate sepset beliefs
      rcptr<Factor> bel01 = m_cl0_to_cl1->absorb(m_cl1_to_cl0)->normalize();
      rcptr<Factor> bel02 = m_cl0_to_cl2->absorb(m_cl2_to_cl0)->normalize();
      rcptr<Factor> bel12 = m_cl1_to_cl2->absorb(m_cl2_to_cl1)->normalize();

      // Check for convergence
      if (iteration > 1) { // Logic only works after one iteration
        double d01 = bel01_old->distance(bel01); // KL-divergence
        double d02 = bel02_old->distance(bel02);
        double d12 = bel12_old->distance(bel12);

        double max_change = std::max({d01, d02, d12});
        if (max_change < threshold) {
          converged = true;
        }

      }

      // Update old beliefs
      bel01_old = bel01;
      bel02_old = bel02;
      bel12_old = bel12;

    }

    // 1.2 g.) Check for calibration
    rcptr<Factor> ptrC0_bel_b0 = cl0->absorb(m_cl2_to_cl0)->absorb(m_cl1_to_cl0)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from Cluster 0: " << *ptrC0_bel_b0 << std::endl;

    rcptr<Factor> ptrC2_bel_b0 = cl2->absorb(m_cl0_to_cl2)->absorb(m_cl1_to_cl2)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from Cluster 2: " << *ptrC2_bel_b0 << std::endl;

    rcptr<Factor> ptrSep01_bel_b0 = m_cl0_to_cl1->absorb(m_cl1_to_cl0)->marginalize({b0})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b0 from separator 0-1: " << *ptrSep01_bel_b0 << std::endl;

    // The beliefs are consistent across the clusters and separator -- calibrated

    // 1.1 (h) 
    rcptr<Factor> ptrb1 = m_cl0_to_cl2->absorb(m_cl2_to_cl0)->marginalize({b1})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b1: " << *ptrb1 << std::endl;
    rcptr<Factor> ptrb2 = m_cl0_to_cl1->absorb(m_cl1_to_cl0)->marginalize({b2})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b2: " << *ptrb2 << std::endl;
    rcptr<Factor> ptrb3 = m_cl1_to_cl2->absorb(m_cl2_to_cl1)->marginalize({b3})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b3: " << *ptrb3 << std ::endl; 
    rcptr<Factor> ptrb4 = cl0->absorb(m_cl1_to_cl0)->absorb(m_cl2_to_cl0)->marginalize({b4})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b4: " << *ptrb4 << std::endl;
    rcptr<Factor> ptrb5 = cl1->absorb(m_cl0_to_cl1)->absorb(m_cl2_to_cl1)->marginalize({b5})->normalize();
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b5: " << *ptrb5 << std::endl;  
    rcptr<Factor> ptrb6 = cl2->absorb(m_cl0_to_cl2)->absorb(m_cl1_to_cl2)->marginalize({b6})->normalize(); 
    std::cout << __FILE__ << ":" << __LINE__ << " - Belief over b6: " << *ptrb6 << std::endl;

    // Question 1.3
    // a.)
    // Factors with no evidence observed
    rcptr<Factor> ptrC0 = ptrR0gB0->absorb(ptrR1gB1)->absorb(ptrR2gB2)->absorb(ptrR3gB3);
    rcptr<Factor> ptrC1 = phi_c1->absorb(ptrR4gB4);
    rcptr<Factor> ptrC2 = phi_c2->absorb(ptrR5gB5);
    rcptr<Factor> ptrC3 = phi_c3->absorb(ptrR6gB6);

    vector< rcptr<Factor> > factorPtrs;
    factorPtrs.push_back(ptrC0);
    factorPtrs.push_back(ptrC1);
    factorPtrs.push_back(ptrC2);
    factorPtrs.push_back(ptrC3);

    // 1110 001
    map<RVIdType, AnyType> obsv;
    obsv[r0] = int(1);
    obsv[r1] = int(1);
    obsv[r2] = int(1);
    obsv[r3] = int(0);
    obsv[r4] = int(0);
    obsv[r5] = int(0);
    obsv[r6] = int(1);

    ClusterGraph cg(ClusterGraph::JTREE, factorPtrs, obsv);
    cg.exportToGraphViz("errorCode_JT"); // Inspection of the PDF file shows that the JT created by emdw is identical to the one I modelled myself.

    // 1.3 c.)
    // 1. Initialize message structures
    // The compiler note suggests using 'Idx2' for the map key
    std::map<Idx2, rcptr<Factor>> msgs; 
    MessageQueue msgQ; // Remove emdw:: prefix if it's causing an error

    // 2. Perform BP message passing
    // Note: loopyBP_CG is used for both JTs and loopy graphs in this framework [cite: 273, 348]
    unsigned nMsgs = loopyBP_CG(cg, msgs, msgQ); 
    std::cout << "Sent " << nMsgs << " messages to calibrate the JT\n";

    // 3. Extract and display marginal beliefs for b0...b6
    for (uint i = b0; i <= b6; ++i) {
        // queryLBP_CG returns the marginalized distribution for variable i [cite: 343]
        rcptr<Factor> qPtr = queryLBP_CG(cg, msgs, {i})->normalize();
        std::cout << "p(b" << i << " | evidence): " << *qPtr << std::endl;
    }

    // 1.4 d/e)
    ClusterGraph cG(ClusterGraph::LTRIP, factorPtrs, obsv);
    cG.exportToGraphViz("errorCode_LTRIP");

        // 1. Initialize message structures
    // The compiler note suggests using 'Idx2' for the map key
    std::map<Idx2, rcptr<Factor>> messages; 
    MessageQueue msgQueue; // Remove emdw:: prefix if it's causing an error

    // 2. Perform BP message passing
    // Note: loopyBP_CG is used for both JTs and loopy graphs in this framework [cite: 273, 348]
    unsigned nMessages = loopyBP_CG(cG, messages, msgQueue); 
    std::cout << "Sent " << nMessages << " messages to calibrate the JT\n";

    // 3. Extract and display marginal beliefs for b0...b6
    for (uint i = b0; i <= b6; ++i) {
        // queryLBP_CG returns the marginalized distribution for variable i [cite: 343]
        rcptr<Factor> qPtr = queryLBP_CG(cG, messages, {i})->normalize();
        std::cout << "Loopy p(b" << i << " | evidence): " << *qPtr << std::endl;
    }











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
