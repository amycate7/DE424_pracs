// quick start info
// echo "7 1E-4 1E-4 0 1E-4 0 0 99 99" | src/bin/emdwlda

/*
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : July 2015
 * Copyright  : University of Stellenbosch, all rights retained
 */

/**
 * Does a Bayesian LDA topic model inference from training data. This
 * version saves a bit on memory by first finishing off one document
 * before it moves on to the next.
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"
#include "prlite_stlvecs.hpp"
#include "prlite_zfstream.hpp"

// emdw headers
#include "emdw.hpp"
#include "dirichlet.hpp"
#include "dirichletset.hpp"
#include "polya.hpp"
#include "conditionalpolya.hpp"

// boost headers
//#include <boost/algorithm/string.hpp>

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <sstream> // stringstream
#include <iomanip>  // setw, setprecision
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>

using namespace std;

typedef string  WordType;
typedef size_t TopicType;

#include "prlite_stringlist.hpp" // StringList


// ######################
// ### tuneable stuff ###
// ######################

double topicAlpha =   +0.0001;    // 5 //many topics per doc
double wordAlpha  =   +0.0001; // 1E-4 // only a few words per topic
double topicPerturb =  0.0;    // Positive uniform noise to add to topics per doc
double wordPerturb  = +0.0001; //1E-2;   // Positive uniform noise to add to words per topic
size_t nClusterNudges = 3;    // 5 number of words peaked per topic cluster
double nudgeWeight  =  0.0;    // 1.0; // extra weight to add to peaked words in topics

unsigned batchSz = 100;  // number of documents before updating the DirichletSet
unsigned nEpochs = 99;
unsigned useSimldaParams = 1;
size_t K = 6;

// Defines an easy set of docs to test with. Three docs with 16 single
// letter words per docs. V=5. It returns a vector of docs with each
// doc a vector of words. It also returns by reference a vector of doc
// centers where each doc center is a vector of words. Also as a side
// effect it returns the vocabulary as it was found in the set of
// docs.
vector< vector<WordType> >
getDocs(
  vector< vector<WordType> >& centers,
  rcptr< vector<WordType> >& vocabVec) {

  centers = {{"A"},{"B"},{"C"}};
  vector< vector<WordType> > docs {
    {"A","A","A","A","A","A","A","A","A","A","A","A","A","B","B","D"},
    {"B","B","B","B","B","B","B","B","B","B","B","B","B","C","C","E"},
    {"C","C","C","C","C","C","C","C","C","C","C","C","C","A","A","D"}
  };

    set<WordType> vocabSet;
    for (auto aDoc : docs) {
      for (auto aWord : aDoc) {
        vocabSet.insert(aWord);
      } // for
    } // for

    vocabVec = uniqptr< vector<WordType> > (new vector<WordType>);
    vocabVec->insert( vocabVec->begin(), vocabSet.begin(), vocabSet.end() );

    return docs;
} // getDocs

unsigned chooseOne(const vector<double> probs) {
  static std::uniform_real_distribution<> randUnif(0.0, 1.0);
  double accum = 0.0;
  double randval = randUnif( emdw::randomEngine() );
  unsigned choice = 0;
  accum = probs[0];
  while (accum <= randval) {
    choice++; accum += probs[choice];
  } // while
  return choice;
} // chooseOne

//  this routine collects the output from simlda into a set of
//  docs. You need to first list the doc*.txt files from simlda in a
//  docs.lst file. It returns a vector of docs with each doc a vector
//  of words. It also returns by reference a vector of doc centers
//  where each doc center is a vector of words. Also as a side effect
//  it returns the vocabulary as it was found in the set of docs.
vector< vector<WordType> >
loadDocs(
  const string& docsListName,
  const string& vocName,
  const string& witName,
  vector< vector<WordType> >& centers,
  rcptr< vector<WordType> >& vocabVec) {

  // The docs itself might not include all words referred to in
  // "wordsInTopics" file. Therefore explicitly load the vocabulary.

  // load the docs
  WordType aWord;
  StringList docList(docsListName);
  vector< vector<WordType> > docs( docList.size() );
  for (unsigned docNr = 0; docNr < docList.size(); docNr++) {

    string docFileName = docList[docNr];
    prlite::gzifstream dfile(docFileName);
    while (dfile) {
      dfile >> aWord;
      if (dfile) {
        docs[docNr].push_back(aWord);
      } // if
    } // while
    dfile.close();

  } // for docNr

  vocabVec = uniqptr< vector<WordType> > (new vector<WordType>);
  prlite::gzifstream vocabfile(vocName);
  while(vocabfile) {
    vocabfile >> aWord;
    if (vocabfile) {(vocabVec->push_back(aWord));} // if
  } // while
  //size_t V = vocabVec->size();


  // use the wordsInTopics file to at least mark the central word in
  // each topic. This is a bit of a cheat, typically we will need some
  // clustering procedure for this.
//  gzifstream tfile(witName);
//  vector< vector<double> > wordsInTopic;
//  size_t v = 0;
//  while (tfile) {
//    double prob;
//    tfile >> prob;
//    if (tfile) {
//
//      // new topic
//      if (!v) {
//        wordsInTopic.push_back( vector<double>() );
//      } // if
//      wordsInTopic[wordsInTopic.size()-1].push_back(prob);
//      if (++v == V) {v = 0;} // if
//    } // if
//  } // while

//  size_t K = wordsInTopic.size();
//  centers.resize(K);
//  for (size_t k = 0; k < K; k++) {
//    centers[k].resize(nClusterNudges);
//    for (size_t i = 0; i < nClusterNudges; i++) {
//      centers[k][i] = (*vocabVec)[chooseOne( wordsInTopic[k] )];
//    } // for
//  } // for

  return docs;
} // loadDocs

/// @brief this calculates the symmetric KL distance between two
/// extremely sparse Dirichlet distributions, each consisting of a
/// single non-zero phantom (NOTE: not alpha) value.
///
/// @param cnt1 phantom value from dist1, assuming all other phantoms are 0
///
/// @param cnt2 phantom value from dist2, assuming all other phantoms are 0
///
/// @param sz the number of categories in each distribution

double
symmKLDist(
  double cnt1,
  double cnt2,
  size_t sz) {

  // example values for v1{0, 1, 2, 1, 4}, v2{1, 2, 2, 1, 1}

  double sumAlpha1 = cnt1+sz;
  double sumAlpha2 = cnt2+sz;
  double sumLogGamma1 = lgamma(cnt1+1.0);
  double sumLogGamma2 = lgamma(cnt2+1.0);
  double diff = cnt1-cnt2;
  double sumDiff1 =  diff * ( digamma(cnt1+1.0) - digamma(sumAlpha1) );
  double sumDiff2 = -diff * ( digamma(cnt2+1.0) - digamma(sumAlpha2) );
  double refDist1 = lgamma(sumAlpha1) - sumLogGamma1 - lgamma(sumAlpha2) + sumLogGamma2 + sumDiff1;
  double refDist2 = lgamma(sumAlpha2) - sumLogGamma2 - lgamma(sumAlpha1) + sumLogGamma1 + sumDiff2;

  return (refDist1 + refDist2)/2.0;
} // symmKLDist

using namespace emdw;



int main(int argc, char *argv[]) {


  unsigned seedVal = emdw::randomEngine.getSeedVal();
  seedVal = 1059638474;
  cout <<  seedVal << endl;
  emdw::randomEngine.setSeedVal(seedVal);

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {


    // #################
    // ### the input ###
    // #################

    if (argc > 12){
      // you can run it like this: ../emdwlda total_topics number_of_docs words_per_doc total_vocab nudge_for_betas default_betas_word_topic default_alphas_topic_doc doc vocab output nr_of_nudges_per_topic epochs batch_size 0 (from where docs live)
      K = stoi(argv[1])+1;  // topics
      nudgeWeight = atof(argv[5]);  // word topic nudge val for betas
      wordAlpha = atof(argv[6]);  // word topic default val for betas
      topicAlpha = atof(argv[7]); // topic doc default val for alphas
      nClusterNudges = stoi(argv[11]);
      nEpochs = stoi(argv[12]);
      batchSz = stoi(argv[13]);
      useSimldaParams = stoi(argv[14]);
      if ((useSimldaParams!=0)&&(useSimldaParams!=1)){
        useSimldaParams = 1;
      }
    }
    else if (argc==9){
      // you can run it like this: ../emdwlda total_topics nudge_for_betas default_betas_word_topic default_alphas_topic_doc nr_of_nudges_per_topic epochs batch_size 0
      K = stoi(argv[1])+1;  // topics + one "distributed" topic
      nudgeWeight = atof(argv[2]);  // word topic nudge val for betas
      wordAlpha = atof(argv[3]);  // word topic default val for betas
      topicAlpha = atof(argv[4]); // topic doc default val for alphas
      nClusterNudges = stoi(argv[5]); // number of words nudged per topic
      nEpochs = stoi(argv[6]); // total number of iterations
      batchSz = stoi(argv[7]); // batch size (update CDIR every batch)
      useSimldaParams = stoi(argv[8]);
      if ((useSimldaParams!=0)&&(useSimldaParams!=1)){
        useSimldaParams = 1;
      }
      cout << nEpochs << batchSz << useSimldaParams<<endl;
    }
    else{
      cout << "This program estimates the parameters of a DirichletSet distribution\n"
           << "from the data generated by the simlda program. It must run from the same directory\n"
           << "where simlda dumped its simulated data. It also needs a 'docs.lst' file, listing\n"
           << "the names of all the simulated docs that simlda generated\n";

      cout << "How many topics K?\n> ";
      cin >> K;
      cout << "\nSupply values for topicAlpha <1E-4>, wordAlpha <1E-4>,\n"
           << "topicPerturb <0>, wordPerturb <1E-4>,\n"
           << "nNudges <0>, nudgeWeight <0>\n> ";
      cin >> topicAlpha >> wordAlpha >> topicPerturb >> wordPerturb >> nClusterNudges >> nudgeWeight;
      cout << "\nHow many message passing epochs <23>?\n> ";
      cin >> nEpochs;
      cout << "\nMini batch size (docs)? Full batch <0>; Recommended minibatch <100>\n> ";
      cin >> batchSz;
    }

    vector< vector<WordType> > centers;
    rcptr< vector<WordType> > wDOM;

    // docs[docNr][wordNr] lists for each document the words it
    // contains

    // a quick test version
    //vector< vector<WordType> > docs = getDocs(centers, wDOM);

    // this one loads a set of (simulated) documents
    vector< vector<WordType> > docs = loadDocs("docs.lst", "vocab.gz", "wordsInTopic.gz", centers, wDOM);

    // ####################
    // ### INITIALIZING ###
    // ####################

    // *****************
    // *** sizes etc ***
    // *****************

    // if (useSimldaParams){
    //   K = centers.size(); // noOfTopics
    // }
    size_t M = docs.size();
    if (batchSz < 1 or batchSz > M) {batchSz = M;} // if

    rcptr< vector<TopicType> > phiDOM ( new vector<TopicType>(K) );
    for (TopicType k = 0; k < K; k++) {
      (*phiDOM)[k] = k;
    } // for

    cout << "Inferring LDA counts with K = " << K
         << "; V = " << wDOM->size()
         << "; M = " << M
         << "; N = " << docs[0].size()
         << endl;
    cout << "Initial values for\n"
         << "topicAlpha: " << topicAlpha
         << "; wordAlpha: " << wordAlpha
         << "; topicPerturb " << topicPerturb
         << "; wordPerturb " << wordPerturb
         << endl;
    cout << "No of nudges per topic cluster: " << nClusterNudges
         << ", each with nudge weight: " << nudgeWeight
         << endl;

    RVIdType idNo = 0;

    // *****************************
    // *** set up all the RV ids ***
    // *****************************

    // Topic DirichletSet at id 0 (Phi)
    RVIdType phiID = idNo++;

    // Now the document Dirichlets (theta)
    vector<RVIdType> thetaID(M);
    map<Idx2, RVIdType> zID;
    map<Idx2, RVIdType> wID;

    for (size_t m = 0; m < M; m++) { // for each doc

      thetaID[m] = idNo++;

      for (size_t n = 0; n < docs[m].size(); n++) { // for each word in each doc
        zID[{m,n}] = idNo++;
        wID[{m,n}] = idNo++;
      } // for n

    } // for m
    cout << "We have " << idNo << " RVs here\n";

    // **************************
    // *** set up the factors ***
    // **************************

    // ++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++ the words per topic DirichletSet (c0_CDIR) +++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++

    std::map< TopicType, rcptr< Dirichlet<WordType> > > c0_DIR;
    for (TopicType k = 0; k < K; k++) {
      c0_DIR[k] = uniqptr< Dirichlet<WordType> >(
        new Dirichlet<WordType>(wDOM, k, wordAlpha, wordPerturb, 1E-7) );
      // remove the line below if you want to only use noise
      // perturbation and no extra info (such as an initial
      // clustering)
//      for (size_t i = 0; i < centers[k].size(); i++) {
//        c0_DIR[k]->addToCounts(centers[k][i], nudgeWeight);
//      } // for
    } // for

    rcptr<Factor> c0_CDIR = uniqptr<Factor>(
      new DirichletSet<TopicType, WordType>(phiID, phiDOM, c0_DIR) );
    c0_CDIR->keepSeparate = true;
    //cout << "Initial topics model: " << *c0_CDIR << endl;
    cout << "Built the DirichletSet\n";

    vector< rcptr<Factor> > c3_DIR(M); // c3_DIR[m]
    map< Idx2, rcptr<Factor> > c2_POL; // polyaPtr[{m,n}]
    map< Idx2, rcptr<Factor> > c1_CPOL;  // c1_CPOL[{m,n}]

    for (size_t m = 0; m < M; m++) { // for every doc

      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // +++ the topics per doc Dirichlets for doc m (c3_DIR[m]) +++
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      c3_DIR[m] = uniqptr< Dirichlet<TopicType> >(
        new Dirichlet<TopicType>(phiDOM, thetaID[m], topicAlpha, topicPerturb, 1E-7) );
      c3_DIR[m]->keepSeparate = true;
      //cout << *c3_DIR[m] << endl;


      for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // +++ the Polyas for the word in document m location n  (c2_POL) +++
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        c2_POL[{m,n}] = uniqptr< Polya<TopicType> >(
          new Polya<TopicType>(thetaID[m], zID[{m,n}], phiDOM) );
        //cout << "Polya: " << *c2_POL[{m,n}] << endl;

        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // +++ the ConditionalPolyas for the word in document m location n  (c1_CPOL) +++
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ValidValue<WordType> obsvW(docs[m][n]);
        c1_CPOL[{m,n}] = uniqptr< ConditionalPolya<TopicType, WordType> >(
          new ConditionalPolya<TopicType, WordType> (phiID, zID[{m,n}], wID[{m,n}], phiDOM, wDOM, obsvW) );
        //cout << m << "," << n << " ConditionalPolya: " << *c1_CPOL[{m,n}] << endl;

      } // for n

    } // for m
    cout << "Built " << c1_CPOL.size() << " ConditionalPolyas\n";
    cout << "Built " << c2_POL.size() << " Polyas\n";
    cout << "Built " << c3_DIR.size() << " Dirichlets\n";

    // *************************************************************
    // *** set up the sepsets and initial forward sweep messages ***
    // *************************************************************

    map<Idx2,RVIds> ss01;  // from c0 to c1[{m,n}]
    map<Idx2,RVIds> ss12;  // from c1[{m,n}] to c2[{m,n}]
    map<Idx2,RVIds> ss23;  // from c2[{m,n}] to c3[m]

    for (size_t m = 0; m < M; m++) { // for every doc
      for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc
        ss01[{m,n}] = RVIds{phiID};
        ss12[{m,n}] = RVIds{zID[{m,n}]};
        ss23[{m,n}] = RVIds{thetaID[m]};
      } // for n
    } // for m
    cout << "Done setting up the sepsets\n";

    // ***************************
    // *** initialize messages ***
    // ***************************

    map< Idx2,rcptr<Factor> > m01;  // from c0 to c1[{m,n}]
    map< Idx2,rcptr<Factor> > m12;  // from c1[{m,n}] to c2[{m,n}]
    map< Idx2,rcptr<Factor> > m23;  // from c2[{m,n}] to c3[m]
    rcptr<Factor> tmp_FAC;

    for (size_t m = 0; m < M; m++) { // for every doc

      // +++++++++++++++++++++
      // ++* forward sweep +++
      // +++++++++++++++++++++

      for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc

        // c0 to c1
        m01[{m,n}] = c0_CDIR->marginalize( ss01[{m,n}], true, c1_CPOL[{m,n}].get() );
        c1_CPOL[{m,n}]->inplaceAbsorb(m01[{m,n}]);
        c1_CPOL[{m,n}]->inplaceNormalize();

        // c1 to c2
        m12[{m,n}] = c1_CPOL[{m,n}]->marginalize(ss12[{m,n}], true);
        c2_POL[{m,n}]->inplaceAbsorb(m12[{m,n}]);
        c2_POL[{m,n}]->inplaceNormalize();

        // c2 to c3
        m23[{m,n}] = c2_POL[{m,n}]->marginalize(ss23[{m,n}], true);
        c3_DIR[m]->inplaceAbsorb(m23[{m,n}]);
        c3_DIR[m]->inplaceNormalize();

      } // for n

      // ++++++++++++++++++++++
      // +++ backward sweep +++
      // ++++++++++++++++++++++

      for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc

        // c3 to c2
        tmp_FAC = c3_DIR[m]->marginalize(ss23[{m,n}], true);
        c2_POL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m23[{m,n}]) );
        c2_POL[{m,n}]->inplaceNormalize();
        m23[{m,n}] = tmp_FAC;

        // c2 to c1
        tmp_FAC = c2_POL[{m,n}]->marginalize(ss12[{m,n}], true);
        c1_CPOL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m12[{m,n}]) );
        c1_CPOL[{m,n}]->inplaceNormalize();
        m12[{m,n}] = tmp_FAC;

      } // for n

    } // for m

    for (size_t m = 0; m < M; m++) { // for every doc
      for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc
        // c1 to c0
        tmp_FAC = c1_CPOL[{m,n}]->marginalize(ss01[{m,n}], true);
        c0_CDIR->inplaceAbsorb( tmp_FAC->cancel(m01[{m,n}]) );
        c0_CDIR->inplaceNormalize();
        m01[{m,n}] = tmp_FAC;

      } // for n
    } // for m
    cout << "Done setting up initial messages\n";



    // *****************************
    // *** loopy message passing ***
    // *****************************

    size_t batchCnt = 0;
    for (int iter = 0; iter < int(nEpochs); iter++) {

      for (size_t m = 0; m < M; m++) { // for every doc

        // +++++++++++++++++++++++++++++++
        // +++ forward sweep           +++
        // +++ Lauritzen Spiegelhalter +++
        // +++++++++++++++++++++++++++++++

        for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc

          // c0 to c1 (note the classical Lauritzen Spiegelhalter procedure)
          tmp_FAC = c0_CDIR->marginalize(ss01[{m,n}], true, c1_CPOL[{m,n}].get() ); // new ss_belief: simply marginalise
          c1_CPOL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m01[{m,n}]) ); // new cl_belief: cancel old and absorb new messages
          c1_CPOL[{m,n}]->inplaceNormalize(); // not sure this is necessary, check
          m01[{m,n}] = tmp_FAC; // and save new ss_belief for next time

          // c1 to c2
          tmp_FAC = c1_CPOL[{m,n}]->marginalize(ss12[{m,n}], true);
          c2_POL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m12[{m,n}]) );
          c2_POL[{m,n}]->inplaceNormalize();
          m12[{m,n}] = tmp_FAC;

          // c2 to c3
          tmp_FAC = c2_POL[{m,n}]->marginalize(ss23[{m,n}], true);
          c3_DIR[m]->inplaceAbsorb( tmp_FAC->cancel(m23[{m,n}]) );
          c3_DIR[m]->inplaceNormalize();
          m23[{m,n}] = tmp_FAC;

        } // for n

        // ++++++++++++++++++++++
        // +++ backward sweep +++
        // ++++++++++++++++++++++

        for (size_t n = 0; n < docs[m].size(); n++) { // for every word in every doc

          // c3 to c2
          tmp_FAC = c3_DIR[m]->marginalize(ss23[{m,n}], true);
          c2_POL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m23[{m,n}]) );
          c2_POL[{m,n}]->inplaceNormalize();
          m23[{m,n}] = tmp_FAC;

          // c2 to c1
          tmp_FAC = c2_POL[{m,n}]->marginalize(ss12[{m,n}], true);
          c1_CPOL[{m,n}]->inplaceAbsorb( tmp_FAC->cancel(m12[{m,n}]) );
          c1_CPOL[{m,n}]->inplaceNormalize();
          m12[{m,n}] = tmp_FAC;

        } // for n

        batchCnt++;
        if (batchCnt == batchSz) {
          batchCnt = 0;
          cout << "." << flush;
          for (int mm = int(m)+1-batchSz; mm <= int(m); mm++){
            size_t docNr = mm < 0 ? mm + M : mm;
            //cout << docNr << " " << flush;

            for (size_t n = 0; n < docs[docNr].size(); n++) { // for every word in every doc

              // c1 to c0
              tmp_FAC = c1_CPOL[{docNr,n}]->marginalize(ss01[{docNr,n}], true);
              c0_CDIR->inplaceAbsorb( tmp_FAC->cancel(m01[{docNr,n}]) );
              c0_CDIR->inplaceNormalize();
              m01[{docNr,n}] = tmp_FAC;

            } // for n

          } // for mm
          //cout << endl;

        } // if

      } // for m


      stringstream snapshot; // left zero-padding
      snapshot << "cdir" << setfill('0') << setw(2) << iter << ".gz";
      string fname = snapshot.str();
      //cout << "#" << flush;
      cout << iter << flush;
      //cout << "Iteration " << fname << " completed\n";
      prlite::gzofstream cdirfile(fname, static_cast<std::ios::openmode>(std::ios::out) );
      rcptr< DirichletSet<TopicType, WordType> > cdPtr =
        std::dynamic_pointer_cast< DirichletSet<TopicType, WordType> >(c0_CDIR);
      for (auto cond : *phiDOM) {
        rcptr< Dirichlet<WordType> > dirPtr = cdPtr->aliasDirichlet(cond);
        for (auto word : *wDOM) {
          cdirfile << dirPtr->getAlpha(word) << " ";
        } // for
        cdirfile << endl;
      } // for
      cdirfile.close();

//      cout << *c0_CDIR << endl;

      prlite::gzofstream tfile("wordsInTopicActual.gz", static_cast<std::ios::openmode>(std::ios::out) );
    rcptr< DirichletSet<TopicType, WordType> > cdPtr1 =
    std::dynamic_pointer_cast< DirichletSet<TopicType, WordType> >(c0_CDIR);
    for (auto cond : *phiDOM) {
    rcptr< Dirichlet<WordType> > dirPtr = cdPtr1->aliasDirichlet(cond);
    for (auto word : *wDOM) {
      tfile << dirPtr->getAlpha(word) << " ";
    } // for
    tfile << endl;
    } // for
    tfile.close();

    } // for iter

    cout << *c0_CDIR << endl;
  } // try

  catch (const string& msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (const char* msg) {
    cerr << msg << endl;
    throw;
  } // catch

  catch (char const * msg) {
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
