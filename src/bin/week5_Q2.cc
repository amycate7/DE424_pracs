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
#include "textblockio.hpp"
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
// This file is adapted from example.cc and provides the framework
// for Question 2 of Week 3's practical assignment.
//
// It reads in the corrupted small image (small_img_noisy.txt)
// as a matrix and prints it to the screen.  It also creates a second
// matrix with the same size, fills it with the value 0.5, and then
// writes it out to a text file (small_img_rec.txt).
//
// It also shows you how to generate identities for the random
// variables in the model.
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
    double defProb = 0.0;           // Any unspecified probs will default to this.
    rcptr< vector<T> > pixDom (     // Lists the values that a particular RV can take on
      new vector<T>{-1,1});         // A pixel can have a value of either -1 or 1

    //*********************************************************
    // Load image from file
    //*********************************************************

    // The loadBlock function is part of the gLinear library and
    // is made available by including textblockio.hpp at the top
    // of the source file.  The contents of the text file is read
    // into the gLinear class gRowMatrix, which represents a matrix.
    // This matrix is called img_noisy here.  
    // You can then access the element at row i and column j by
    // referring to img_noisy(i,j).  You can get the number of 
    // rows by using img_noisy.rows() and the number of columns
    // by using img_noisy.cols().  
    
    gLinear::gRowMatrix<int> img_noisy;
    img_noisy = loadBlock<int>("small_img_noisy.txt");

    int rows = img_noisy.rows();
    int cols = img_noisy.cols();

    //*********************************************************
    // Print the contents of the matrix to the screen
    //*********************************************************

    std::cout << std::endl << "Contents of small_img_noisy.txt:" << std::endl;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
	std::cout << img_noisy(i,j) << " ";
      }
      std::cout << std::endl;
    }

    //*********************************************************
    // Create the RV identities
    //*********************************************************

    // Even for a 4x4 image, it would be cumbersome to specify
    // RV identities with enum{pix00, pix01, pix02, ...}.  A much
    // better approach is to create a map (or dictionary) with the
    // coordinates of the pixels in the image as they key and the
    // associated RV identity as value.  For each pixel coordinate,
    // we need two RVs: one for the latent RV x and one for the
    // observed RV y.  We therefore create two maps: xRVs and yRVs.
    // To refer to the latent RV x at row i and column j, write
    // xRVs[{i,j}], and to refer to the observed RV y, write
    // yRVs[{i,j}]. 

    emdw::RVIdType currentRVId = 0; // initial RV identity
    std::map<std::pair<int,int>, emdw::RVIdType> xRVs;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
	xRVs.insert({{i,j}, currentRVId++});
      }
    }
    std::map<std::pair<int,int>, emdw::RVIdType> yRVs;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
	yRVs.insert({{i,j}, currentRVId++});
      }
    }

    double match = 10.0;
    double mismatch = 1.0;

    // Indices: {-1,-1}, {-1,1}, {1,-1}, {1,1}
    DT potentialTable({pixDom, pixDom}, defProb, {
        {{-1, -1}, match},
        {{-1, 1},  mismatch},
        {{1, -1},  mismatch},
        {{1, 1},   match}
    });

    vector<rcptr<Factor>> factorPtrs;
    std::map<emdw::RVIdType, T> obsv; // Map for observed RV values 

    // 2. Iterate and Create Factors 
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            
            // Factor: p(xi, yi) ---
            emdw::RVIds obsIds = { xRVs[{r,c}], yRVs[{r,c}] };
            vector<rcptr<vector<T>>> obsDoms = { pixDom, pixDom };
            
            rcptr<Factor> ptrObs = uniqptr<DT>(new DT(obsIds, obsDoms, defProb, potentialTable.getMap(), 0.0, 0.0, false));
            factorPtrs.push_back(ptrObs);

            obsv[yRVs[{r,c}]] = T(img_noisy(r,c));

            // Binary/Interaction Factors: p(xi, xj)
            if (c + 1 < cols) {
                emdw::RVIds horizIds = { xRVs[{r,c}], xRVs[{r, c+1}] };
                vector<rcptr<vector<T>>> horizDoms = { pixDom, pixDom };
                
                rcptr<Factor> ptrHoriz = uniqptr<DT>(new DT(horiz_ids, horiz_doms, defProb, potentialTable.getMap(), 0.0, 0.0, false));
                factorPtrs.push_back(ptrHoriz);
            } 
            if (r + 1 < rows) {
                emdw::RVIds vertIds = { xRVs[{r,c}], xRVs[{r+1, c}] };
                vector<rcptr<vector<T>>> vertDoms = { pixDom, pixDom };

                rcptr<Factor> ptrVert = uniqptr<DT>(new DT(vert_ids, vert_doms, defProb, potentialTable.getMap(), 0.0, 0.0, false));
                factorPtrs.push_back(ptrVert);
            }
        }
    }

    // 3. Construct ClusterGraph and Perform Inference
    ClusterGraph cg(ClusterGraph::JTREE, factorPtrs, obsv);

    std::map<Idx2, rcptr<Factor>> msgs;
    MessageQueue msgQ;

    // Perform message passing to calibrate the graph 
    unsigned nMsgs = loopyBP_CG(cg, msgs, msgQ);
    std::cout << "Sent " << nMsgs << " messages to calibrate the Junction Tree." << std::endl;

    // 4. Extract Marginal Beliefs for the Reconstructed Image
    gLinear::gRowMatrix<float>img_rec(rows,cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Query the marginal for the latent variable x at (i,j) 
            rcptr<Factor> qPtr = queryLBP_CG(cg, msgs, {xRVs[{i,j}]})->normalize();
            
            // Cast and get the probability of state '1'
            rcptr<DT> dtMarginal = std::static_pointer_cast<DT>(qPtr);
            img_rec(i,j) = (float)dtMarginal->getProb({1});
        }
    }


    // Save the results for Python visualization 
    // (Followed by your existing ofstream logic to save img_rec)

    // For this practical, we explicitly construct the joint
    // distribution over all the RVs in the model by multiplying
    // all the factors together.  An easy way to do this is to
    // start with a factor that functions as a "running product"
    // and then step through the factors of the model one-by-one,
    // creating each factor and then immediately multiplying it
    // with the running product.
    // To initialise this running product, one could use one
    // of the factors of the model; however, it is easier to
    // start with an uninformative factor -- a factor that leaves
    // another factor unchanged if it is multiplied with it.
    // Such an uninformative factor can be implemented as a table
    // containing only entries of 1.  This uninformative factor
    // is created below (ptrJoint).  
    
    // rcptr<Factor> ptrJoint = uniqptr<DT>(new DT(
    //   {xRVs[{0,0}]}, // latent RV x for pixel at coordinates (0,0)
    //   {pixDom}, // the domain for each pixel value; i.e., {-1, 1}
    //   defProb,
    //   {
    //     {{-1}, 1}, // this factor is "flat" or "uninformative"
    //     {{1}, 1}, // i.e., all table entries are the same
    //   } ));

    
    //*********************************************************
    // Here you should iterate over the factors of the model,
    // create each factor, and multiply it with the running
    // product to form the joint distribution of the model.

    // You should also insert the evidence (i.e., the pixels of
    // the noisy/corrupted image) into the joint distribution.
    // It would be more efficient to insert the evidence into
    // the factor over each pixel's x and y RVs *before* you
    // multiply this factor with the running product, instead
    // of first constructing the joint distribution over all
    // the x's and y's of the model, and then inserting the
    // evidence.  
    //*********************************************************

    
    //*********************************************************
    // Create a matrix to store the reconstructed image
    //*********************************************************

    // Since the result of inference will be a probability for
    // each pixel, we create a matrix with float elements.  The
    // matrix is then filled with values of 0.5 -- you will of
    // course rather have to fill it with the probabilities that
    // result from inference.  
    
  //   gLinear::gRowMatrix<float>img_rec(rows,cols);

  //   for (int i = 0; i < rows; i++) {
  //     for (int j = 0; j < cols; j++) {
	// img_rec(i,j) = 0.5; // replace this line with the posterior
	// // probability that the pixel at coordinates (i,j) has a
	// // value of one, which you should extract from the joint
	// // distribution. 
  //     }
  //   }
    
    //*********************************************************
    // Save reconstructed image to file
    //*********************************************************

    // To write out the reconstructed image to text file, we simply
    // open the file as an output stream and pipe out the contents
    // of the matrix that contains the reconstructed image. 
    
    ofstream recFile;
    recFile.open("small_img_rec.txt");
    if (!recFile.is_open()) {
      std::cout << "Error: Could not open file to write out reconstructed image!" << std::endl;
    }
    else {
      for (int i=0; i < rows; i++) {
	for (int j=0; j < cols; j++) {
	  recFile << img_rec(i,j) << ' ';
	}
	recFile << std::endl;
      }
      recFile.close();
      std::cout << std::endl << "The reconstructed image has been written to small_img_rec.txt" << std::endl;
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