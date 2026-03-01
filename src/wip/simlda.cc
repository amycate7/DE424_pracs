// quick start info
// rm doc*.gz; echo "7 100 500 50 3 0.5" | ../src/wip/simlda; find doc*.gz | sort > docs.lst

/*
 * Author     :  JA du Preez (DSP Group, E&E Eng, US)
 * Created on : Febr 2015
 * Copyright  : University of Stellenbosch, all rights retained
 */

/**
 * This code generates simulated documents following the LDA model i.e
 * a document consists of a mixture of topics and a topic consists of
 * a mixture of words.
 *
 * For our words we choose the integers in the range [0:V-1].
 *
 * To make the visual identification of a topic easier, we choose the
 * words in the topic to be close to each other in numerical
 * value. (The words is firstly arranged in a circle i.e. the last one
 * is adjacent to the first one, and then a topic is centered around a
 * particular position on this circle with a Gaussian decline to both
 * sides. To this We add an additional topic, non-overlapping with the
 * others, but occurring in all documents.
 *
 * Each document is composed of two topics (although we can change
 * this of course).
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"
#include "prlite_mrandom.hpp"
#include "prlite_shuffles.hpp"
#include "prlite_stlvecs.hpp"
#include "prlite_zfstream.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>

#include <algorithm>
#include <limits>
#include <iomanip>  // setw, setprecision
#include <sstream>

using namespace std;

// ######################
// ### tuneable stuff ###
// ######################

unsigned K = 11;   // core no of topics, another common topic is added automatically.
unsigned V = 100;  // core vocab size
unsigned M = 1000; // no of docs;
unsigned N = 100; // no of words per doc
unsigned topicsPerDoc = 3; // core topics per doc
double overlap = 0.5; // weighting factor to indicate overlap of gaussians


// gaussian shape on words in topic
double gkernel(double x, double stdev) {
  double ret = exp( -x*x/(2*stdev*stdev) );
  if (ret < 1E-4) {ret = 0.0;} // if 5E-3
  return ret;
} // gkernel

// laplace shape on words in topic
double lkernel(double x, double stdev) {
  double b = stdev/sqrt(2.0);

  double ret = exp(-abs(x)/b); // /(2*b);
  if (ret < 1E-4) {ret = 0.0;} // if
  return ret;
} // lkernel

unsigned chooseOne(const vector<double> probs) {
  static prlite::MRandom rand;
  double accum = 0.0;
  double randval = rand.Flat();
  unsigned choice = 0;
  accum = probs[0];
  while (accum <= randval) {
    choice++; accum += probs[choice];
  } // while
  return choice;
} // chooseOne

int main(int argc, char *argv[]) {

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {
    // we can read in parameters as arguments or at run time

    // argument parameter option
    if (argc == 7){
      K = stoi(argv[1]);  // topics
      cout <<"nr of core topics: "<< K<<endl;
      topicsPerDoc =  stoi(argv[2]); // topics per doc
      cout <<"topics per doc: "<< topicsPerDoc<<endl;
      M = stoi(argv[3]);  // docs
      cout <<"number of docs: "<< M<<endl;
      N = stoi(argv[4]);  // words per doc
      cout <<"words per doc: "<< N<<endl;
      V = stoi(argv[5]);  // vocab size
      cout <<"vocab size: "<< V<<endl;
      overlap = atof(argv[6]);  // vocab size
      cout <<"overlap: "<< overlap<<endl;

    }
    else{
      cout << "You can supply the following arg values as parameters or input them manually: nrOfTopics topicsPerDoc nrOfDocs wordsPerDoc vocabSize overlap";
      cout << "Supply values for K, V, M, N, topicsPerDoc and overlap\n> ";
      cin >> K >> V >> M >> N;

      cout << "Simulating LDA data with K = " << K
         << "; V = " << V
         << "; M = " << M
         << "; N = " << N
         << "; topicsPerDoc = " << topicsPerDoc
         << "; overlap = " << overlap
         << endl;
    }

    double stdev = overlap*V/K; //0.25*V/K;
    unsigned nCommon = 0.05*V;   // extra vocab size

    // for doing zero-filled formatting
    unsigned widthV = unsigned(log10(V+nCommon) + 1); // no of chars in the word index
    unsigned widthM = unsigned(log10(M) + 1); // no of chars in the doc index

    // for simplicity we will take our core words to be 0, 1, ... V-1

    // We'll arrange our core words in a circle. Each core topic will
    // be placed with a focus center somewhere on the circle,
    // decreasing in word probability as you move away from this
    // center.

    // Now let's define the topic word weights (phi). Topic k
    // concentrates on the words around V/K*(k+1) // JADP !!! just recheck this

    // NOTE: There actually are K+1 topics, where the last one (index
    // K) represents a topic that is present in all documents. Same as
    // with actual text grouping function words in one topic.
    vector< vector<double> > wordsInTopic(K+1);
    for (unsigned k = 0; k < K; k++) {
      double centerAt = V/K*(k+1);
      wordsInTopic[k].resize(V+nCommon);
      double sum = 0.0;
      for (unsigned v = 0; v < V; v++) {
        double dist1 = abs(centerAt - v);
        double dist2 = V - abs(v-centerAt);
        double dist = min(dist1,dist2);
        wordsInTopic[k][v] = gkernel(dist,stdev);
        sum += wordsInTopic[k][v];
      } // for v
      for (unsigned v = 0; v < V; v++) {
        wordsInTopic[k][v] /= sum;
      } // for
      for (unsigned v = V; v < V+nCommon; v++) {
        wordsInTopic[k][v] = 0.0;
      } // for v
    } // for k
    wordsInTopic[K].resize(V+nCommon);
    for (unsigned v = 0; v < V; v++) {
      wordsInTopic[K][v] = 0.0;
    } // for v
    for (unsigned v = V; v < V+nCommon; v++) {
      wordsInTopic[K][v] = 1.0/nCommon;
    } // for

    // ok, now we define the wordsInTopic that a doc consists of (theta)
    prlite::MRandom randnum;
    vector< vector<double> > topicsInDoc(M);   // theta
    for (unsigned m = 0; m < M; m++) {
      double sum = 0.0;
      topicsInDoc[m].resize(K+1);
      for (auto& val : topicsInDoc[m]) {val = 0.0;}
      prlite::ShuffledSequence shuff(K);
      for (unsigned k = 0; k < topicsPerDoc; k++) {
        topicsInDoc[m][shuff[k]] = randnum.Flat();
        sum += topicsInDoc[m][shuff[k]];
      } // for
      topicsInDoc[m][K] = randnum.Flat();
      sum += topicsInDoc[m][K];
      for (auto& val : topicsInDoc[m]) {val /= sum;}
    } // for m

    // at last we are in the position to start generating docs
    vector< vector<string> > docs(M);
    for (unsigned m = 0; m < M; m++) { // iterate over the docs
      stringstream docNr; // left zero-padding
      docNr << setfill('0') << setw(widthM) << m;
      //cout << dicNr.str() << " ";
      string docName = "doc"+ docNr.str() + ".gz";
      prlite::gzofstream docfile(docName);
      for (unsigned n = 0; n < N; n++) { // iterate over the words in a doc
        unsigned aWrd;
        aWrd = chooseOne( wordsInTopic[ chooseOne(topicsInDoc[m]) ] );
        stringstream wordNr; // left zero-padding
        wordNr << setfill('0') << setw(widthV) << aWrd;
        //cout << wordNr.str() << " ";
        docs[m].push_back( wordNr.str() );
        docfile << wordNr.str() << " ";
      } // for n
      docfile << endl;
      docfile.close();
    } // for m
    //cout << docs << endl;

    prlite::gzofstream vocabfile("vocab.gz");
    for (unsigned v = 0; v < V+nCommon; v++) {
      stringstream wordNr; // left zero-padding
      wordNr << setfill('0') << setw(widthV) << v;
      vocabfile << wordNr.str() << " ";
    } // for v
    vocabfile.close();

    prlite::gzofstream tfile("wordsInTopic.gz");
    for (unsigned k = 0; k < K+1; k++) {
      for (unsigned v = 0; v < V+nCommon; v++) {
        tfile << wordsInTopic[k][v] << " ";
      } // for
      tfile << endl;
    } // for
    tfile.close();


    prlite::gzofstream dfile("topicsInDoc.gz");
    for (unsigned m = 0; m < M; m++) {
      dfile << m << ": ";
      for (unsigned k = 0; k < K+1; k++) {
        dfile << topicsInDoc[m][k] << " ";
      } // for
      dfile << endl;
    } // for
    dfile.close();

  } // try

  catch (const string& msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const char* msg) {
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
