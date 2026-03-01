/*******************************************************************************

  AUTHOR: JA du Preez
  DATE: May 2014
  PURPOSE: Implements the non-linear part for the UKF
  transforms. Transforms a ColVector<double> into another one.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

#include "v2vtransform.hpp"

// standard headers
#include <iostream>

// patrec headers
#include "prlite_genvec.hpp"
#include "prlite_genmat.hpp"
#include "prlite_outerquad_l2u.hpp"

// emdw headers
#include "sigmapoints.hpp"
#include "matops.hpp"
//#include "sortindices.hpp"
#include "vecset.hpp"

using namespace std;


// ---------------------- Exponentiation

vector< prlite::ColVector<double> > Invert::operator()(const prlite::ColVector<double>& x) const {
  vector< prlite::ColVector<double> > y(1); y[0].resize( x.size() );
  for (int i = 0; i < x.size(); i++) {
    y[0][i] = 1.0/x[i];
  } // for i
  return y;
} // operator()

vector< prlite::ColVector<double> > Exponentiate::operator()(const prlite::ColVector<double>& x) const {
  vector< prlite::ColVector<double> > y(1); y[0].resize( x.size() );
  for (int i = 0; i < x.size(); i++) {
    y[0][i] = pow(x[i],pwr);
  } // for i
  return y;
} // operator()

// ---------------------- Distances between points

vector< prlite::ColVector<double> > DistanceFromKnownPoint::operator()(const prlite::ColVector<double>& x) const {
  int dim = x.size();
  PRLITE_ASSERT(dim == x0.size(), "We need two equally sized vectors");
  vector< prlite::ColVector<double> > dist(1); dist[0].resize(1); dist[0][0] = 0.0;
  for (int i = 0; i < dim; i++) {
    double diff = x[i] - x0[i];
    dist[0][0] += diff*diff;
  } // for
  dist[0][0] = sqrt(dist[0][0]);
  return dist;
} // operator()

vector< prlite::ColVector<double> > DistancesFromKnownPoints::operator()(const prlite::ColVector<double>& x) const {
  int dim = x.size();
  unsigned nPnts = xs.size();
  vector< prlite::ColVector<double> > dist(1); dist[0].resize(nPnts); dist[0].assignToAll(0.0);
  for (unsigned n = 0; n < nPnts; n++) {
    PRLITE_ASSERT(dim == xs[n].size(), "We need two equally sized vectors");
    for (int i = 0; i < dim; i++) {
      double diff = x[i] - xs[n][i];
      dist[0][n] += diff*diff;
    } // for i
    dist[0][n] = sqrt(dist[0][n]);
    //cout << "dist is: " << dist[0][n] << endl;
  } // for n
  return dist;
} // operator()

// -------------------------------------------------------------
//                   DistanceBetweenUnknownPoints
// -------------------------------------------------------------

DistanceBetweenUnknownPoints::DistanceBetweenUnknownPoints(
    const emdw::RVIds& rvIds1,
    const emdw::RVIds& rvIds2)
    : v1Idx(0),
      v2Idx(0) {

  unsigned inDim = rvIds1.size();
  PRLITE_ASSERT(inDim == rvIds2.size(), "The two position vectors must have the same dimension");
  vector<size_t> l2u, r2u;
  emdw::RVIds vecU = sortedUnion(rvIds1, rvIds2, l2u, v1Idx, r2u, v2Idx);
  if (v2Idx.size() != inDim) { // seems we have some overlap in variables - redo v2Idx
    v2Idx.resize(inDim);
    size_t cnt = 0;
    for (size_t idx = 0; idx < vecU.size();  idx++) {
      if (rvIds2[cnt] == vecU[idx]) {
        v2Idx[cnt++] = idx;
        if (cnt == inDim) break;
      } // iff
    } // for
  } // if

  //cout << v1Idx << endl;
  //cout << v2Idx << endl;

} // ctor

vector< prlite::ColVector<double> >
DistanceBetweenUnknownPoints::operator()(
  const prlite::ColVector<double>& x) const {
  size_t dim = v1Idx.size();
  vector< prlite::ColVector<double> > dist(1); dist[0].resize(1); dist[0][0] = 0.0;
  for (size_t idx = 0; idx < dim; idx++) {
    double diff = x[ v1Idx[idx] ] - x[ v2Idx[idx] ];
    dist[0][0] += diff*diff;
  } // for
  dist[0][0] = sqrt(dist[0][0]);
  //cout << dist[0][0] << " " << flush;
  return dist;
} // operator()

// -------------------------------------------------------------

// ---------------------- Matrix vector products Ax

vector< prlite::ColVector<double> > MpyKnownAMatWithUnknownXVec::operator()(const prlite::ColVector<double>& x) const {
  vector< prlite::ColVector<double> > y(1);
  y[0] = aMat*x;
  return y;
} // operator()

vector< prlite::ColVector<double> > MpyUnknownAMatWithKnownXVec::operator()(const prlite::ColVector<double>& a) const {
  vector< prlite::ColVector<double> > y(1);
  prlite::RowMatrix<double> aMat(a.getStartAddr(), outDim, inDim);
  y[0] = aMat*x;
  return y;
} // operator()

vector< prlite::ColVector<double> > MpyUnknownAMatWithUnknownXVec::operator()(const prlite::ColVector<double>& xa) const {
  vector< prlite::ColVector<double> > y(1);
  prlite::ColVector<double> x(xa.getStartAddr(), inDim);
  prlite::RowMatrix<double> aMat(xa.getStartAddr()+inDim, outDim, inDim);
  y[0] = aMat*x;
  return y;
} // operator()

// ---------------------- Gaussian priors

// <mn> is 'n veranderbare vektor wat tipiese samples/voorbeelde van
// die VERWAGTE_WAARDE van y verskaf. (Dit word tipies bepaal deur
// klient kode in NormedGaussCanonical na aanleiding van 'n
// verspreiding/Gaussiese prior van die GEMIDDELDE van y).
//
// <cov> is die vooraf bekende kovariansie van y.
//
// Hierdie roetine verskaf tipiese samples van y deur die 2D+1
// sigmapunte op die Gaussiese verspreiding N(mn,cov) te bepaal.

vector< prlite::ColVector<double> > UnknownMeanKnownCov::operator()(const prlite::ColVector<double>&  mn) const {
  vector< prlite::ColVector<double> > ySpts;
  vector<double> wgts;
  unsigned dim = mn.size();
  cov2sp_2Dp1(mn, cov, 1.0/(2.0*dim+1), ySpts, wgts);
  return ySpts;
} // operator()

// <covAsVec> is 'n veranderbare vektor wat tipiese samples/voorbeelde
// van die KOVARIANSIE van y verskaf. (Dit word tipies bepaal deur
// klient kode in NormedGaussCanonical na aanleiding van 'n
// verspreiding/Gaussiese prior van die KOVARIANSIE van y).
//
// <mn> is die vooraf bekende gemiddelde van y.
//
// Hierdie roetine verskaf tipiese samples van y deur die 2D+1
// sigmapunte op die Gaussiese verspreiding N(mn,cov) te bepaal.

vector< prlite::ColVector<double> > KnownMeanUnknownCov::operator()(const prlite::ColVector<double>&  covAsVec) const {
  vector< prlite::ColVector<double> > ySpts;
  vector<double> wgts;
  unsigned dim = mn.size();

  prlite::RowMatrix<double> cov = packVecToLowerDiags(covAsVec, dim);
  prlite::copyLowerToUpper(cov);
  //cout << "Mean: " << mn << ", Cov: " << cov << endl;

  cov2sp_2Dp1(mn, cov, 1.0/(2.0*dim+1), ySpts, wgts);
  //cout << "Sigma points: "; for (auto elem : ySpts) cout << elem << endl;
  //exit(-1);
  return ySpts;
} // operator()

// <covYAsVec> is 'n veranderbare vektor wat tipiese samples/voorbeelde
// van die KOVARIANSIE cov(Y) verskaf. (Dit word tipies bepaal deur
// klient kode in NormedGaussCanonical na aanleiding van 'n
// verspreiding/Gaussiese prior van die KOVARIANSIE van y).
//
// <mn> is die vooraf bekende gemiddelde van die onbekende gemiddeld van y, naamlik mny.
//
// Hierdie roetine verskaf tipiese samples van y deur die 2D+1
// sigmapunte op die Gaussiese verspreiding N(mn,cov) te bepaal.

vector< prlite::ColVector<double> > KnownMeanScaledUnknownCov::operator()(const prlite::ColVector<double>&  covYAsVec) const {
  vector< prlite::ColVector<double> > ySpts;
  vector<double> wgts;
  unsigned dim = mn.size();

  prlite::RowMatrix<double> cov = packVecToLowerDiags(covYAsVec, dim);
  // scale the cov down with beta
  for (unsigned r = 0; r < dim; r++) {
    for (unsigned c = 0; c <= r; c++) {
      cov(r,c) /= beta;
    } // for c
  } // for r
  prlite::copyLowerToUpper(cov);
  //cout << "Mean: " << mn << ", Cov: " << cov << endl;

  cov2sp_2Dp1(mn, cov, 1.0/(2.0*dim+1), ySpts, wgts);
  //cout << "Sigma points: "; for (auto elem : ySpts) cout << elem << endl;
  //exit(-1);
  return ySpts;
} // operator()

//  std::vector< ColVector<double> > operator()(const ColVector<double>& covYAsVecAndMnY) cons

vector< prlite::ColVector<double> > UnknownMeanUnknownCov::operator()(const prlite::ColVector<double>&  covYAndMnYAsVec) const {
  vector< prlite::ColVector<double> > ySpts;
  vector<double> wgts;

  prlite::RowMatrix<double> cov = packVecToLowerDiags(covYAndMnYAsVec, dim);
  prlite::copyLowerToUpper(cov);
  prlite::ColVector<double> mn(covYAndMnYAsVec.getStartAddr()+dim*(dim+1)/2, dim);
  //cout << "Mean: " << mn << ", Cov: " << cov << endl;
  //exit(-1);

  cov2sp_2Dp1(mn, cov, 1.0/(2.0*dim+1), ySpts, wgts);
  //cout << "Sigma points: "; for (auto elem : ySpts) cout << elem << endl;
  //exit(-1);
  return ySpts;
} // operator()
