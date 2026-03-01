/*******************************************************************************

  AUTHOR: JA du Preez
  DATE: May 2014

  PURPOSE: Implements the gaussian to sigmapoint and back transformations

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// this file's header
#include "sigmapoints.hpp"

// standard headers
#include <iostream>

// gLinear headers
#include "gLinear/gLapackWrap.h"

// patrec headers
#include "prlite_genvec.hpp"
#include "prlite_genmat.hpp"
#include "prlite_outerquad_l2u.hpp"

// emdw headers
#include "matops.hpp"

using namespace std;

double sgn(double x) {
  return x < 0 ? -1.0 : 1.0;
} // sgn

prlite::RowMatrix<double> corr_coef(const prlite::RowMatrix<double>& cov) {
  prlite::RowMatrix<double> cc( cov.rows(), cov.cols() );
  for (int r = 0; r < cc.rows(); r++) {
    for (int c = 0; c <= r; c++) {
      cc(r,c) = cov(r,c)/sqrt(cov(r,r)*cov(c,c));
    } // for c
  } // for r
  prlite::copyLowerToUpper(cc);
  return cc;
} //corrcoef

void
cov2sp_2Dp1(
  const prlite::ColVector<double>& mn,
  const prlite::RowMatrix<double>& cv,
  double w0,
  vector< prlite::ColVector<double> >& spts,
  vector<double>& wgts) {

  unsigned dim = mn.size();
  unsigned nSamples = 2*dim+1;
  //double w0 = 1.0-dim; //0.5/nSamples;   // 1-dim places the sps at 1 stdev.
  double wi = (1.0-w0)/(nSamples-1);
  double scal = sqrt(dim/(1.0-w0));
  //cout << "weights: " << w0 << " " << wi << " scal: " << scal << endl;

  prlite::RowMatrix<double> L(dim,dim);
  int fail;
  cholesky_factor(cv, L, fail);

  if (fail) {
    cerr << __FILE__ << __LINE__ << " Fixing non-posdef cov!\n";
    prlite::RowMatrix<double> cv_reg(cv);
    fixCovarianceMat(cv_reg);
    cholesky_factor(cv_reg, L, fail);
  } // if
  if (fail) {cout << corr_coef(cv) << endl;}
  PRLITE_ASSERT(!fail, "Could not Cholesky factorise regularised covariance matrix");

  // the sigma points
  spts.resize(nSamples);
  wgts.resize(nSamples);
  spts[0] = mn; wgts[0] = w0;

  for (unsigned i = 1; i <= dim; i++) {
    spts[i] = mn + scal*L.col(i-1); wgts[i] = wi;
    spts[i+dim] = mn - scal*L.col(i-1); wgts[i+dim] = wi;
  } // for

} // cov2sp_2Dp1

void
cov2sp_4Dp1(
  const prlite::ColVector<double>& mn,
  const prlite::RowMatrix<double>& cv,
  vector< prlite::ColVector<double> >& spts,
  vector<double>& wgts) {

  unsigned dim = mn.size();
  unsigned nSamples = 4*dim+1;
  //double w0 = 1.0-dim; //0.5/nSamples;   // 1-dim places the sps at 1 stdev.
  double wi = 1.0; ///nSamples;
  //cout << "weights: " << w0 << " " << wi << " scal: " << scal << endl;

  prlite::RowMatrix<double> L(dim,dim);
  int fail;
  cholesky_factor(cv, L, fail);
  PRLITE_ASSERT(!fail, "Could not Cholesky factorise the covariance matrix");
  //cout << "L: " << L << endl;

  // the sigma points
  spts.resize(nSamples);
  wgts.resize(nSamples);
  spts[0] = mn; wgts[0] = wi;

  prlite::ColVector<double> tmp;
  unsigned offset;
  for (unsigned i = 1; i <= dim; i++) {
    tmp = L.col(i-1);

    spts[i]        = mn + tmp;
    wgts[i]        = wi;

    offset         = dim;
    spts[i+offset] = mn - tmp;
    wgts[i+offset] = wi;

    offset         = 2*dim;
    double scal    = sqrt(2.0*dim -0.5);
    tmp *= scal;
    spts[i+offset] = mn + tmp;
    wgts[i+offset] = wi;

    offset         = 3*dim;
    spts[i+offset] = mn - tmp;
    wgts[i+offset] = wi;
  } // for

} // cov2sp_4Dp1

void
sp2cov(
  vector< prlite::ColVector<double> >& spts,
  vector<double>& wgts,
  prlite::ColVector<double>& mn,
  prlite::RowMatrix<double>& cv) {

  // lets reestimate
  double cnt = wgts[0];
  prlite::ColVector<double> sumx = spts[0]*wgts[0];
  prlite::RowMatrix<double> sumx2 = outerProdSym<double>(spts[0], sumx);
  for (int i = 1; i < int(spts.size()); i++) {
    prlite::ColVector<double> spw(spts[i]*wgts[i]);
    sumx += spw;
    addOuterProdSym(sumx2, spts[i], spw);
    cnt += wgts[i];
  } // for i

  mn = sumx/cnt;
  cv = ( sumx2-outerProdSym<double>(mn,sumx) ) / cnt;
  prlite::copyLowerToUpper(cv);

} // sp2cov

// Bullet-proof fix for deviant covariance matrices: Limit the min
// vars and max corrs, and then up the vars to force diag dominance.
void
fixCovarianceMat(
  prlite::RowMatrix<double>& covMat,
  double minVar,
  double maxCorr,
  double diagDominancy) {
  int dim = covMat.rows();
  for (int r = 0; r < dim; r++) {
    // threshold too small variances
    covMat(r,r) = max(covMat(r,r), minVar);
    double sumOffDiag = 0.0;
    for (int c = 0; c < r; c++) {
      // limit the extreme correlations
      double maxV = sqrt( covMat(r,r)*covMat(c,c) );
      if ( abs(covMat(r,c)) > maxV*maxCorr) { // clip the correlation
        covMat(r,c) = sgn( covMat(r,c) )*maxV*maxCorr;
        covMat(c,r) = covMat(r,c);
      } // if
      sumOffDiag += abs(covMat(r,c));
    } // for c
    // and up the the variances till it's diagonally dominant
    for (int c = r+1; c < dim; c++) {
      sumOffDiag += abs(covMat(r,c));
    } // for c
    covMat(r,r) = max(covMat(r,r), sumOffDiag*diagDominancy);
  } // for r
} // fixCovarianceMat

// forces the information matrix to be diagonally dominant by reducing the off-diags
void
fixPrecisionMat(
  prlite::RowMatrix<double>& precMat,
  double maxPrec,
  double maxCorr,
  double diagDominancy) {
  //cout << precMat;
  int dim = precMat.rows();
  for (int r = 0; r < dim; r++) {
    // limit too big precisions
    precMat(r,r) = min(precMat(r,r), maxPrec);
    double sumOffDiag = 0.0;
    for (int c = 0; c < r; c++) {
      // limit the extreme correlations
      double maxV = sqrt( precMat(r,r)*precMat(c,c) );
      if ( abs(precMat(r,c)) > maxV*maxCorr) { // clip the correlation
        precMat(r,c) = sgn( precMat(r,c) )*maxV*maxCorr;
        precMat(c,r) = precMat(r,c);
      } // if
      sumOffDiag += abs(precMat(r,c));
    } // for c
    for (int c = r+1; c < dim; c++) {
      sumOffDiag += abs(precMat(r,c));
    } // for c
    if (precMat(r,r) > 0) { //soften the correlations till diag dominant
      for (int c = 0; c < dim; c++) {
        if (c == r) continue;
        precMat(r,c) *= precMat(r,r)/sumOffDiag/diagDominancy;
        precMat(c,r) = precMat(r,c);
      } // for c
    } // if
    else { // something went horribly wrong ... we try heh
      precMat(r,r) = max(precMat(r,r), sumOffDiag*diagDominancy);
    } // else
  } // for r
  //cout << precMat;
  //exit(-1);
} // fixPrecisionMat
