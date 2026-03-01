/*
 * Author     :  (DSP Group, E&E Eng, US)
 * Created on :
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
#include "prlite_logging.hpp"  // initLogging
#include "prlite_testing.hpp"
#include "prlite_mrandom.hpp"
#include "shuffles.hpp"
//#include "textblockfeatures.hpp"
#include "prlite_mean_cov.hpp"
#include "gausgen.hpp"

// emdw headers
#include "emdw.hpp"
#include "sparsetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"
#include "normedgausscanonical.hpp"
#include "v2vtransform.hpp"
#include "matops.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <memory>
#include <set>
#include <map>
#include <algorithm>
#include <limits>

using namespace std;
using namespace emdw;

prlite::RowMatrix<double> corrcoef(const prlite::RowMatrix<double>& cov) {
  prlite::RowMatrix<double> cc( cov.rows(), cov.cols() );
  for (int r = 0; r < cc.rows(); r++) {
    for (int c = 0; c <= r; c++) {
      cc(r,c) = cov(r,c)/sqrt(cov(r,r)*cov(c,c));
    } // for c
  } // for r
  copyLowerToUpper(cc);
  return cc;
} //corrcoef

typedef NormedGaussCanonical GCF;
enum class Status {unknown, estimated, known};
enum class EstType {VI,NET};

template <typename T>
RVVals toRVVals(const T& x) {
  unsigned sz = x.size();
  RVVals ret(sz);
  for (unsigned i = 0; i < sz; i++) {
    ret[i] = x[i];
  } // for
  return ret;
} //toRVVals

uniqptr<GCF> buildIsotropicGaussian(
  double nomMn,
  double nomVar,
  unsigned dim,
  unsigned firstVarId = 0) {

  // the variables in this pdf
  RVIds varsMn(dim);
  setVals(varsMn, firstVarId, 1u);

  // the mean of this pdf
  prlite::ColVector<double> meanOfMn(dim); meanOfMn.assignToAll(nomMn);
  // the cov of this pdf
  prlite::RowMatrix<double> covOfMn(dim,dim);
  covOfMn.assignToAll(0.0);
  for (unsigned i = 0; i < dim; i++) {
    covOfMn(i,i) = nomVar;
  } // for i

  return uniqptr<GCF>( new GCF(varsMn, meanOfMn, covOfMn) );
} // buildIsotropicGaussian

// if the feature vector has dim D, then this prior will have dim
// D(D+1)/2. It is packed from the covariance matrix as, first the
// main diagonal, then working down the lower diagonals. The variable
// IDs are numbered in this same order. Beta is a safety factor that
// allows us to later scale the variances down without running into
// positive definiteness issues (see Gauss-Wishart distribution).
uniqptr<GCF> buildCovPriorGaussian(
  double nomVar,
  unsigned dim,
  double beta = 2,
  unsigned firstVarId = 0) {

  // the variables in this pdf
  unsigned jointDim = (dim*dim+dim)/2;
  // the first dim is the main diagonal, then the lower diagonals
  RVIds varsCov(jointDim);
  setVals(varsCov, firstVarId, 1u);

  double srhoii = min(0.25, 1.0/sqrt(dim+0.5));
  double srhoij = min(0.25, 1.0/sqrt( 0.5*dim*(dim-1.0)+0.5) );


  // first the means of this pdf
  prlite::ColVector<double> meanOfCov(jointDim);
  for (unsigned i = 0; i < dim; i++) {
    meanOfCov[i] = nomVar;
  } // for i
  for (unsigned i = dim; i < jointDim; i++) {
    meanOfCov[i] = 0.0;
  } // for i

  // now the covs of this pdf
  prlite::RowMatrix<double> covOfCov(jointDim,jointDim);
  covOfCov.assignToAll(0.0);
  for (unsigned i = 0; i < dim; i++) {
    double v1 = nomVar*srhoii/beta;
    covOfCov(i,i) = v1*v1;
  } // for i
  for (unsigned i = dim; i < jointDim; i++) {
    double v1 = nomVar*srhoij/beta;
    covOfCov(i,i) = v1*v1;
  } // for i

  return uniqptr<GCF>( new GCF(varsCov, meanOfCov, covOfCov) );

} //buildCovPriorGaussian


int main(int, char *argv[]) {

  initLogging(argv[0]);
  TestCase::runAllTests();

  try {

    MRandom randnum(7L);
    ShuffledSequence shuffle(1,7);

    // cout << *buildCovPriorGaussian(1,2,2,0);
    // exit(-1);
    // ------------ set up the basic type definitions ------------

    if (0) {  // Bayesian estimation of a Gaussian CovMat.

      // this is a no go -- the idea of a Gaussian for a CovMat just
      // has to many symmetries to it! NOT WORKING

      // --------- general parameters

      unsigned nVecs = 333;
      unsigned dim = 2;
      unsigned noOfVars = 0;
      double nomVarOfCov = 2.0;
      double qVar = 1E-4;
      double beta = 2.0;

      // --------  generate for observation data

      // the true mean of X
      prlite::ColVector<double> theKnownMn(2); theKnownMn[0] = 3.0; theKnownMn[1] = 2.0;
      // the true cov of X, we will assume this
      prlite::RowMatrix<double> theUnknownCov(2,2);
      theUnknownCov(0,0) = 1.5;
      theUnknownCov(1,1) = 2.5;
      theUnknownCov(0,1) = theUnknownCov(1,0) = 1.25;
      GaussianGenerator ggen(theKnownMn.transpose(), theUnknownCov);
      prlite::RowMatrix<double> xs(nVecs,dim);
      for (int n = 0; n < int(nVecs); n++) {
        xs[n] = ggen.generate();
        //cout << xs[n];
      } // for n

      // ---------  prior distributions

      // the measurement noise
      prlite::RowMatrix<double> qCovMat(gLinear::eye<double>(dim)*qVar);

      // first the cov

      rcptr<GCF> priorForCov = buildCovPriorGaussian(nomVarOfCov, dim, beta, noOfVars);
      //cout << "Prior Cov: " << *priorForCov << endl;
      RVIds varsCov = priorForCov->getVars();
      //cout << "varsCov: " << varsCov << endl;
      priorForCov->export2DMesh("before.txt", varsCov[0], varsCov[2], 512);
      //lets replace it with wider version so that nVec times product is original
      // ColVector<double> newH = priorForCov->getH()/nVecs;
      // Matrix<double> newK = priorForCov->getK()/nVecs;
      // priorForCov = uniqptr<GCF>( new GCF(varsCov, newK, newH) );
      noOfVars = priorForCov->noOfVars();

      // now the mean
      prlite::ColVector<double> meanOfMn(dim); meanOfMn.assignToAll(0.0);
      RVIds varsMn(dim); setVals(varsMn, noOfVars, 1u);
      rcptr<GCF> jointPrior =
        uniqptr<GCF>(
          new GCF(priorForCov.get(),
                  KnownMeanScaledUnknownCov(meanOfMn, beta),
                  varsMn, qCovMat) );
      noOfVars += dim;
      // cout << "Joint Cov: " << *jointPrior << endl;
      // exit(-1);

      // --------- observe the data and build the cluster graph

      vector< rcptr<Factor> > factors;
      RVIds varsX(dim); setVals(varsX, noOfVars, 1u);
      //cout << "varsX: " << varsX << endl;

      for (unsigned n = 0; n < nVecs; n++) {

        rcptr<GCF> joint =
          uniqptr<GCF>(
            new GCF(jointPrior.get(), UnknownMeanUnknownCov(dim), varsX, qCovMat) );
        //cout << "Joint Cov: " << *joint << endl;
        //cout << "Joint CC: " << corrcoef(joint->getCov() ) << endl;
        factors.push_back( joint->Factor::observeAndReduce( varsX, toRVVals(xs[n]) ) );
        //cout << "Obsv Cov: " << *factors[0] << endl;
        //exit(-1);
        // rcptr<GCF> qPtr = std::dynamic_pointer_cast<GCF>(factors[0]);
        // cout << *qPtr << endl;
        // qPtr->export2DMesh("tmp.txt", 0, 1, 512);

      } // for n

      rcptr<ClusterGraph> cgPtr;
      map< Idx2, rcptr<Factor> > msgs;
      MessageQueue msgQ;
      cgPtr = uniqptr<ClusterGraph>(new ClusterGraph(factors));
      cgPtr->exportToGraphViz(string("cg_cov"));
      cout << cgPtr->noOfFactors() << " factors in graph\n";

      // ---------- Do inference

      //cout << "Now doing inference\n";
      msgs.clear(); msgQ.clear();
      unsigned nMsgs = loopyBP_CG(*cgPtr, msgs, msgQ, 0.0);
      cout << "Sent " << nMsgs << " messages before convergence\n";

      // ---------- Do the query for the result

      rcptr<Factor> queryPtr = queryLBP_CG( *cgPtr, msgs, append(varsCov, varsMn) );
      //rcptr<Factor> queryPtr = queryLBP_CG( *cgPtr, msgs, varsCov);
      rcptr<GCF> queryG = std::dynamic_pointer_cast<GCF>(queryPtr);
      cout << *queryG << endl;
      queryG->export2DMesh("after.txt", varsMn[0], varsMn[1], 512);
      //queryG->export2DMesh("after.txt", varsCov[0], varsCov[2], 512);

      exit(-1);
    } // if

    if (1) {  // Bayesian estimation of a Gaussian mean.

      // --------- general parameters

      unsigned nVecs = 33;
      unsigned dim = 2;
      unsigned noOfVars = 0;
      double nomVarOfMn = 4.0;
      double qVar = 1E-4;

      // --------  generate for observation data

      // the true mean of X
      prlite::ColVector<double> theUnknownMn(2); theUnknownMn[0] = 3.0; theUnknownMn[1] = 2.0;
      // the true cov of X, we will assume this
      prlite::RowMatrix<double> theKnownCov(2,2);
      theKnownCov(0,0) = 1.5;
      theKnownCov(1,1) = 2.5;
      theKnownCov(0,1) = theKnownCov(1,0) = 1.25;
      GaussianGenerator ggen(theUnknownMn.transpose(), theKnownCov);
      prlite::RowMatrix<double> xs(nVecs,dim);
      for (int n = 0; n < int(nVecs); n++) {
        xs[n] = ggen.generate();
      } // for n

      // ---------  prior distribution for the mean

      rcptr<GCF> priorForMean = buildIsotropicGaussian(0.0, nomVarOfMn, dim, noOfVars);
      RVIds varsMn = priorForMean->getVars();
      cout << "varsMn: " << varsMn << endl;
      priorForMean->export2DMesh("before.txt", varsMn[0], varsMn[1], 512);
      //lets replace it with wider version so that nVec times product is original
      prlite::ColVector<double> newH = priorForMean->getH()/nVecs;
      prlite::RowMatrix<double> newK = priorForMean->getK()/nVecs;
      priorForMean = uniqptr<GCF>( new GCF(varsMn, newK, newH) );
      // the assumed known covariance matrix for the mean
      // Matrix<double> nomCovOfMn(dim,dim); nomCovOfMn.assignToAll(0.0);
      // for (unsigned r = 0; r < dim; r++) {nomCovOfMn(r,r) = nomVarOfMn;} // for
      noOfVars = priorForMean->noOfVars();

      // --------- observe data and build the cluster graph

      vector< rcptr<Factor> > factors;
      RVIds varsX(dim); setVals(varsX, noOfVars, 1u);
      cout << "varsX: " << varsX << endl;
      // the measurement noise
      prlite::RowMatrix<double> qCovMat(gLinear::eye<double>(dim)*qVar);


      for (unsigned n = 0; n < nVecs; n++) {

        rcptr<GCF> joint =
          uniqptr<GCF>(
            new GCF(priorForMean.get(), UnknownMeanKnownCov(theKnownCov),varsX, qCovMat) );
        factors.push_back( joint->Factor::observeAndReduce( varsX, toRVVals(xs[n]) ) );
        // rcptr<GCF> qPtr = std::dynamic_pointer_cast<GCF>(factors[0]);
        // cout << *qPtr << endl;
        // qPtr->export2DMesh("tmp.txt", 0, 1, 512);

      } // for n

      rcptr<ClusterGraph> cgPtr;
      map< Idx2, rcptr<Factor> > msgs;
      MessageQueue msgQ;
      cgPtr = uniqptr<ClusterGraph>(new ClusterGraph(factors));
      cgPtr->exportToGraphViz(string("cg_mean"));
      cout << cgPtr->noOfFactors() << " factors in graph\n";

      // ---------- Do inference

      //cout << "Now doing inference\n";
      msgs.clear(); msgQ.clear();
      unsigned nMsgs = loopyBP_CG(*cgPtr, msgs, msgQ, 0.0);
      cout << "Sent " << nMsgs << " messages before convergence\n";

      // ---------- Do the query for the result

      rcptr<Factor> queryPtr = queryLBP_CG( *cgPtr, msgs, varsMn);
      rcptr<GCF> queryG = std::dynamic_pointer_cast<GCF>(queryPtr);
      cout << *queryG << endl;
      queryG->export2DMesh("after.txt", varsMn[0], varsMn[1], 512);
      exit(-1);
    } // if

  } // try

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
