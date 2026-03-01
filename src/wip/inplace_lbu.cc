/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// #######################################################################
// # This is the same setup as in Minka's clutter problem, only much
// # more generic. In the clutter problem we assume only the mean of
// # the target to be unknown, while the clutter mean, clutter
// # variance, target variance and relative proportions are
// # known. Here we are assuming that all these quantities (except
// # for the relative proportions), are also unknown.
// # It should be relatively easy (via a Dirichlet prior and its Polya
// # dependent) to also handle an unknown relative proportion.
// #######################################################################

// patrec headers
#include "prlite_zfstream.hpp"
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"
#include "prlite_genmat.hpp"

// emdw headers
#include "emdw.hpp"
#include "dirichlet.hpp"
#include "polya.hpp"
#include "gauss1vmp.hpp"
#include "cgauss1vmp.hpp"
//#include "estcgauss1_settings.hpp"
#include "discretetable.hpp"
#include "oddsandsods.hpp"
#include "textblockio.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>  // infinity
#include <iomanip>
#include <chrono>

using namespace std;
using namespace emdw;

const unsigned rBlocks = 1;  // no of vertical subblocks
const unsigned cBlocks = 1;  // no of horizontal subblocks

// pick a block to process
const unsigned rr = min(0u,rBlocks-1);
const unsigned cc = min(0u,cBlocks-1);

const unsigned maxIters = 17;

// parameters for priors
const double p1 = 0.4; // prior guess for % black

// first the black
const double muBlackMean = 0.05;   // our guess for the mean of black x
const double muBlackStDev = 0.01; // the stdev we allow for this mean
const double gaBlackMean  = 50;    // our guess for the precision of black x
const double gaBlackCounts = 10.0; // guess based on this many phantom samples

// now the white
const double muWhiteMean = 0.75;   // our guess for the mean of white x
const double muWhiteStDev  = 0.1;  // the stdev we allow for this mean
const double gaWhiteMean  = 30;    // our guess for the precision of white x
const double gaWhiteCounts = 10.0; // guess based on this many phantom samples

using CondType = unsigned;
using DT = DiscreteTable<CondType>;
using DIR = Dirichlet<CondType>;
using POL = Polya<CondType>;
rcptr< vector<CondType> > binDom (
  new vector<CondType>{0,1});


int main(int, char *argv[]) {

  // NOTE: this activates logging and unit tests
  initLogging(argv[0]);
  prlite::TestCase::runAllTests();

  try {

    char* iName = argv[1];

    prlite::RowMatrix<double> samples;
    samples = loadBlock<double>(iName);
    unsigned nRows = samples.rows();
    unsigned nCols = samples.cols();
    cout << nRows << " " << nCols << endl;

    float rIncr = float(nRows)/rBlocks; // rows in a block
    float cIncr = float(nCols)/cBlocks; // cols in a block

    // the subblock for estimating probs
    //for (unsigned rr = 0; rr < rBlocks; rr++) {
    //for (unsigned cc = 0; cc < cBlocks; cc++) {

    unsigned cBegin = static_cast<unsigned>(cc*cIncr);
    unsigned rBegin = static_cast<unsigned>(rr*rIncr);
    float cEnd = (cc+1)*cIncr;
    float rEnd = (rr+1)*rIncr;

    // #####################
    // # Build the PGM net #
    // #####################

    // we use a simplified version where Z is fixed to P(Z=1) = w;
    // Later on we will expand this to include a Dirichlet prior.

    // *****************************
    // *** set up all the RV ids ***
    // *****************************

    RVIdType idNo = 0;

    // mu prior
    // RVIdType muID = idNo++;

    // ga prior
    // RVIdType gaID = idNo++;

    // Dirichlet prior
    RVIdType thID = idNo++;

    // Now the z ids
    prlite::RowMatrix<RVIdType> zID(nRows, nCols);
    for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
      for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {
        zID(r,c) = idNo++;
      } // for
    } // for

    // *********************************
    // *** Build the various factors ***
    // *********************************

    // +++++++++++++++++++++++
    // +++ Dirichlet prior +++
    // +++++++++++++++++++++++

    double effCnt = 200;
    cout << "sum of alphas: " << effCnt << endl;
    rcptr<Factor> c3_thDIR = uniqptr<DIR>(
      new DIR(thID, binDom, 0.0, { {0, (1-p1)*effCnt}, {1, p1*effCnt} }) );
    //cout << *c3_thDIR << endl;

    // +++++++++++++++++++++++++++++++
    // +++ Priors for the Gaussian +++
    // +++++++++++++++++++++++++++++++

    std::map< CondType, rcptr<Gauss1EXPF> > priorMuMap;
    priorMuMap[0] = std::make_unique<Gauss1EXPF>( ::muWhiteMean, 1.0/(::muWhiteStDev*::muWhiteStDev) );
    priorMuMap[1] = std::make_unique<Gauss1EXPF>( ::muBlackMean, 1.0/(::muBlackStDev*::muBlackStDev) );
    rcptr< Gauss1EXPFset<CondType> > c0_muGAU( std::make_unique< Gauss1EXPFset<CondType> >(priorMuMap) );

    std::map< CondType, rcptr<Gamma1EXPF> > priorGaMap;
    priorGaMap[0] = std::make_unique<Gamma1EXPF>(::gaWhiteMean/::gaWhiteCounts, ::gaWhiteCounts);
    priorGaMap[1] = std::make_unique<Gamma1EXPF>(::gaBlackMean/::gaBlackCounts, ::gaBlackCounts);
    rcptr< Gamma1EXPFset<CondType> > c0_gaGAM( std::make_unique< Gamma1EXPFset<CondType> >(priorGaMap));

    // +++++++++++++++++++++++++++++++++++
    // +++ Conditional Gaussians for X +++
    // +++++++++++++++++++++++++++++++++++

    prlite::RowMatrix< rcptr< CGauss1VMP<CondType> > > c1_xGAU(nRows, nCols);
    for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
      for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {
        c1_xGAU(r,c) = std::make_unique< CGauss1VMP<CondType> >  ( zID(r,c), binDom, samples(r,c) );
        if (c < 10 or c > nCols-10 or r > nRows-10) {
          c1_xGAU(r,c)->observeCategorical(0u);
        } // if
      } // for
    } // for

    // ++++++++++++++++++++
    // +++ Polyas for Z +++
    // ++++++++++++++++++++

    prlite::RowMatrix< rcptr<Factor> > c2_zPOL(nRows, nCols);
    for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
      for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {
        c2_zPOL(r,c) =
          std::make_unique<POL>( POL(thID, zID(r,c), binDom) );
      } // for
    } // for

    // ##########################
    // # Now do message passing #
    // ##########################

    // ***********************************
    // *** create the initial messages ***
    // ***********************************

    prlite::RowMatrix< rcptr< Gauss1EXPFset<CondType> > > m01_muGAU(nRows, nCols);
    prlite::RowMatrix< rcptr< Gamma1EXPFset<CondType> > > m01_gaGAM(nRows, nCols);
    prlite::RowMatrix< rcptr<Factor> > m12_zCAT(nRows, nCols);
    prlite::RowMatrix< rcptr<Factor> > m23_thDIR(nRows, nCols);

    // all basically vacuous
    for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
      for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {
        m01_muGAU(r,c) = std::unique_ptr<  Gauss1EXPFset<CondType> >( c0_muGAU->vacuousCopy() );
        m01_gaGAM(r,c) = std::unique_ptr<  Gamma1EXPFset<CondType> >( c0_gaGAM->vacuousCopy() );
        m12_zCAT(r,c) =  std::make_unique<DT>( DT({zID(r,c)}, {binDom}, 0.5) );
        m23_thDIR(r,c) = std::make_unique<DIR>( DIR(thID, binDom) );
      } // for
    } // for

    rcptr< Gauss1EXPFset<CondType> > tmp_MU;
    rcptr< Gamma1EXPFset<CondType> > tmp_GA;
    rcptr<Factor> tmp_FAC;
    rcptr<DT> tmp_DT;
    rcptr<Factor> tmp2_FAC;
    rcptr<DT> tmp2_DT;


        for (unsigned iter = 0; iter < maxIters; iter++) {

          // *********************
          // *** Forward sweep ***
          // *********************

          auto t1 = std::chrono::high_resolution_clock::now();
          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              //FILEPOS << " sample " << n << " is: " << samples[n] << std::endl;

              // ++++++++++++++++++++++++++++
              // +++ From the priors to X +++
              // +++ c0 -> c1             +++
              // ++++++++++++++++++++++++++++

              tmp_MU = unique_ptr< Gauss1EXPFset<CondType> >( c0_muGAU->copy() );  // new
              m01_muGAU(r,c)->inplaceDiv(tmp_MU);
              c1_xGAU(r,c)->inplaceCancelMeans_BP( m01_muGAU(r,c) );
              // dont need to normalize
              m01_muGAU(r,c) = tmp_MU;  // delete

              tmp_GA = unique_ptr< Gamma1EXPFset<CondType> >( c0_gaGAM->copy() );
              m01_gaGAM(r,c)->inplaceDiv(tmp_GA);
              c1_xGAU(r,c)->inplaceCancelPrecs_BP( m01_gaGAM(r,c) );
              // dont need to normalize
              m01_gaGAM(r,c) = tmp_GA;

              //FILEPOS << *c1_xGAU(r,c) << endl;

            } // for
          } // for

          auto t2 = std::chrono::high_resolution_clock::now();
          auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c0->c1: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();

          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              // +++++++++++++++++++
              // +++ From X to Z +++
              // +++ c1 -> c2    +++
              // +++++++++++++++++++

              tmp_FAC = unique_ptr<Factor>( c1_xGAU(r,c)->marginaliseToCat_BP() );
              //m12_zCAT(r,c)->inplaceCancel(tmp_FAC); !!!
              c2_zPOL(r,c)->inplaceAbsorb(tmp_FAC); // new cl_belief: cancel old and absorb new messages
              c2_zPOL(r,c)->inplaceCancel( m12_zCAT(r,c) ); // new cl_belief: cancel old and absorb new messages
              c2_zPOL(r,c)->inplaceNormalize(); // not sure this is necessary, check
              m12_zCAT(r,c) = tmp_FAC; // and save new ss_belief for next time
              //FILEPOS << *m12_zCAT(r,c) << endl;
              //FILEPOS << *c2_zPOL(r,c) << endl;

            } // for
          } // for

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c1->c2: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();
          //exit(-1);

          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              // +++++++++++++++++++++++
              // +++ From Z to theta +++
              // +++ c2 -> c3        +++
              // +++++++++++++++++++++++

              tmp_FAC = c2_zPOL(r,c)->marginalize({thID}, true);
              //m23_thDIR(r,c)->inplaceCancel(tmp_FAC); !!!
              c3_thDIR->inplaceAbsorb(tmp_FAC);
              c3_thDIR->inplaceCancel( m23_thDIR(r,c) );
              c3_thDIR->inplaceNormalize();
              m23_thDIR(r,c) = tmp_FAC;

              //FILEPOS << *c3_thDIR << endl;

            } // for
          } // for

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c2->c3: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();
          //exit(-1);
          //FILEPOS << *c3_thDIR << endl;

          // **************************************
          // *** lets not get over-enthusiastic ***
          // **************************************

          rcptr<DIR> c3_dwn = dynamic_pointer_cast<DIR>(c3_thDIR);
          double soc = c3_dwn->sumAlphas();
          double p0 = c3_dwn->getAlpha( (*binDom)[0]) / soc;
          double pthresh = 0.25;
          if (p0 < pthresh) {
            c3_thDIR = uniqptr<DIR>(
              new DIR(thID, binDom, 0.0, { {0, pthresh*soc}, {1, (1.0-pthresh)*soc} }) );
          } // if
          else if (p0 > (1.0-pthresh)) {
            c3_thDIR = uniqptr<DIR>(
              new DIR(thID, binDom, 0.0, { {0, (1.0-pthresh)*soc}, {1, pthresh*soc} }) );
          } // else if

          // **********************
          // *** Backward sweep ***
          // **********************

          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              // +++++++++++++++++++++++
              // +++ From theta to Z +++
              // +++ c3 -> c2        +++
              // +++++++++++++++++++++++

              tmp_FAC = c3_thDIR->marginalize({thID}, true);
              //m23_thDIR(r,c)->inplaceCancel(tmp_FAC); !!!
              c2_zPOL(r,c)->inplaceAbsorb(tmp_FAC);
              c2_zPOL(r,c)->inplaceCancel( m23_thDIR(r,c) );
              c2_zPOL(r,c)->inplaceNormalize();
              m23_thDIR(r,c) = tmp_FAC;

            } // for
          } // for

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c3->c2: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();

          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              // ++++++++++++++++++++++
              // +++ BP from Z to X +++
              // +++ c2 -> c1       +++
              // ++++++++++++++++++++++

              tmp_FAC = unique_ptr<Factor>( c2_zPOL(r,c)->marginalize({zID(r,c)}, true ) );
              //m12_zCAT(r,c)->inplaceCancel(tmp_FAC); !!!
              //tmp_DT = dynamic_pointer_cast<DT>(m12_zCAT(r,c));
              tmp_DT = dynamic_pointer_cast<DT>(tmp_FAC);
              c1_xGAU(r,c)->inplaceAbsorbCat_BP(tmp_DT);
              tmp_DT = dynamic_pointer_cast<DT>(m12_zCAT(r,c));
              c1_xGAU(r,c)->inplaceCancelCat_BP(tmp_DT); // new cl_belief: cancel old and absorb new messages
              m12_zCAT(r,c) = tmp_FAC; // and save new ss_belief for next time

            } // for
          } // for

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c2->c1: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();

          for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
            for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {

              // ++++++++++++++++++++++++++++
              // +++ From X to the priors +++
              // +++ c1 -> c0             +++
              // ++++++++++++++++++++++++++++

              tmp_MU = unique_ptr< Gauss1EXPFset<CondType> >( c1_xGAU(r,c)->marginaliseToMeans_BP() );
              m01_muGAU(r,c)->inplaceDiv(tmp_MU);
              //c0_muGAU->inplaceMpy(tmp_MU);
              c0_muGAU->inplaceDiv( m01_muGAU(r,c) );
              // dont need to normalize
              m01_muGAU(r,c) = tmp_MU;

              tmp_GA = unique_ptr< Gamma1EXPFset<CondType> >( c1_xGAU(r,c)->marginaliseToPrecs_BP() );
              m01_gaGAM(r,c)->inplaceDiv(tmp_GA);
              //c0_gaGAM->inplaceMpy(tmp_GA);
              c0_gaGAM->inplaceDiv( m01_gaGAM(r,c) );
              // dont need to normalize
              m01_gaGAM(r,c) = tmp_GA;

            } // for
          } // for

          t2 = std::chrono::high_resolution_clock::now();
          duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
          std::cout << "c1->c0: " << duration << std::endl;
          t1 = std::chrono::high_resolution_clock::now();

          //FILEPOS << *c0_muGAU << endl << *c0_gaGAM << endl;
        } // for iter

        for (unsigned r = rBegin; r < static_cast<unsigned>(rEnd); r++) {
          for (unsigned c = cBegin; c < static_cast<unsigned>(cEnd); c++) {
            samples(r,c) = m12_zCAT(r,c)->potentialAt({zID(r,c)},{0u});
          } // for
        } // for

        //} // for cc
    //} // for rr
    saveBlock<double>("probs.gz", samples);

  } // try

  catch (char msg[]) {
    cerr << msg << endl;
  } // catch

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
