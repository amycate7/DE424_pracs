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
    // observed RV y. We therefore create two maps: xRVs and yRVs.
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

    //*********************************************************
    // Create the factors of the model and build the joint distribution
    //*********************************************************

    double match = 10.0;
    double mismatch = 1.0;

    DiscreteTable table({match, mismatch, mismatch, match});

    // 1. Initialize the joint distribution as an 'identity' factor (value 1.0)
    // We'll start it with the first latent variable to give it a scope
    DiscreteFactor joint(X[0][0], {1.0, 1.0}); 

    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            
            // --- FACTOR 1: Latent xi and Observed yi ---
            DiscreteFactor f_obs({X[r][c], Y[r][c]}, table);
            
            // Efficiency: Insert evidence (yi) before multiplying
            // This reduces f_obs from a 2x2 factor to a 1x2 factor (unary)
            int observed_val = (int)noisy_img(r, c); 
            DiscreteFactor f_unary = f_obs.observe(Y[r][c], observed_val); // Observe Y first - much faster than keeping the Y in the joint
            
            // Multiply into joint (only if it's not the very first assignment)
            if (r == 0 && c == 0) joint = f_unary;
            else joint.absorb(f_unary);

            // --- FACTOR 2: Neighbors (xi and xj) ---
            // Right neighbor
            if (c + 1 < cols) {
                joint.absorb(DiscreteFactor({X[r][c], X[r][c+1]}, table));
            }
            // Down neighbor
            if (r + 1 < rows) {
                joint.absorb(DiscreteFactor({X[r][c], X[r+1][c]}, table));
            }
        }
    }

    // 3. Normalize the final joint distribution
    joint.normalize();





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
    
    rcptr<Factor> ptrJoint = uniqptr<DT>(new DT(
      {xRVs[{0,0}]}, // latent RV x for pixel at coordinates (0,0)
      {pixDom}, // the domain for each pixel value; i.e., {-1, 1}
      defProb,
      {
        {{-1}, 1}, // this factor is "flat" or "uninformative"
        {{1}, 1}, // i.e., all table entries are the same
      } ));

    
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
    
    gLinear::gRowMatrix<float>img_rec(rows,cols);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
	img_rec(i,j) = 0.5; // replace this line with the posterior
	// probability that the pixel at coordinates (i,j) has a
	// value of one, which you should extract from the joint
	// distribution. 
      }
    }
    
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
