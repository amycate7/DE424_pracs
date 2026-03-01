#ifndef SIGMAPOINTS_HPP
#define SIGMAPOINTS_HPP

/*******************************************************************************

  AUTHOR: JA du Preez
  DATE: May 2014

  PURPOSE: Implements the gaussian to sigmapoint and back transformations

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_genvec.hpp"
#include "prlite_genmat.hpp"

// standard headers
#include <vector>

/**
 * @param w0 The weight to use for the centroid. If chosen as
 *   1.0/(2D+1) all 2D+1 sigma points will be equally weighted.
 *
 * NOTE: patrec mostly uses row vectors - will have to see how we
 * interface that.
 */

void
cov2sp_2Dp1(
  const prlite::ColVector<double>& mn,
  const prlite::RowMatrix<double>& cv,
  double w0,
  std::vector< prlite::ColVector<double> >& spts,
  std::vector<double>& wgts);

void
cov2sp_4Dp1(
  const prlite::ColVector<double>& mn,
  const prlite::RowMatrix<double>& cv,
  std::vector< prlite::ColVector<double> >& spts,
  std::vector<double>& wgts);

void
sp2cov(
  std::vector< prlite::ColVector<double> >& spts,
  std::vector<double>& wgts,
  prlite::ColVector<double>& mn,
  prlite::RowMatrix<double>& cv);

/**
 * A fairly bullet-proof fix for deviant covariance matrices: Limit
 * the minimum variances and maximum correlations, and then up the
 * variances to force diagonal dominance.
 *
 * @param minVar Variances lower than this will be upped to this.
 *
 * @param maxCorr Correlation coefficients of which the absolute value
 * exceeds this, will be reduced to this value (while of course
 * keeping its sign intact).
 *
 * @param diagDominancy No variance will be allowed to go lower than
 * this number times the absolute sum of the covariances in the same
 * row.
 *
 * NOTE: Since this does change
 * the covariance matrix it should only be used on covariance matrices
 * that need to be fixed.
 */
void
fixCovarianceMat(
  prlite::RowMatrix<double>& covMat,
  double minVar = 1E-7,
  double maxCorr = 0.75,
  double diagDominancy = 1.05);

/**
 * Only use this in Gaussian information form, in canonical form
 * changes to K will implicitly shift the mean also.
 */
void
fixPrecisionMat(
  prlite::RowMatrix<double>& precMat,
  double maxPrec=1E7,
  double maxCorr = 0.75,
  double diagDominancy=1.05);

#endif // SIGMAPOINTS_HPP
