/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    March 2013
          PURPOSE: Various image processing routines

*******************************************************************************/

#include "inftheory.hpp"
//#include "prlite_stdfun.hpp" // I

#include "prlite_vector.hpp" // FIXME must be before prlite_stlvecs.hpp
#include "prlite_stlvecs.hpp"
#include <limits>

using namespace std;

/*****************************************************************************/

void
hist1d(
  const prlite::RowVector<float>& x,
  prlite::RowVector<float>& xc,
  prlite::RowVector<float>& p,
  unsigned nBins) {

  xc.resize(nBins);
  p.resize(nBins);
  p.assignToAll(0.0);

  float minX = x.minVal();
  float maxX = x.maxVal();
  float xStep = (maxX-minX)/nBins;
  xc.set(minX+xStep/2, xStep);

  unsigned xLen = x.size();
  for (unsigned n = 0; n < xLen; n++) {
    unsigned idx = std::min( nBins-1, unsigned( (x[n]-minX)/xStep ) );
    p[idx] += 1;
  } // for i

  p /= (xLen*xStep);
} // hist1D

void
hist2d(
  const prlite::RowVector<float>& x,
  const prlite::RowVector<float>& y,
  prlite::RowVector<float>& xc,
  prlite::RowVector<float>& yc,
  prlite::RowMatrix<float>& p,
  unsigned nXBins,
  unsigned nYBins) {

  PRLITE_ASSERT(x.size() == y.size(), "x and y must have same length");

  xc.resize(nXBins);
  yc.resize(nYBins);
  p.resize(nXBins,nYBins);
  p.assignToAll(0.0);

  float minX = x.minVal();
  float maxX = x.maxVal();
  float xStep = (maxX-minX)/nXBins;
  xc.set(minX+xStep/2, xStep);

  float minY = y.minVal();
  float maxY = y.maxVal();
  float yStep = (maxY-minY)/nYBins;
  yc.set(minY+yStep/2, yStep);

  unsigned xLen = x.size();
  for (unsigned n = 0; n < xLen; n++) {
    unsigned xIdx = std::min( nXBins-1, unsigned( (x[n]-minX)/xStep ) );
    unsigned yIdx = std::min( nYBins-1, unsigned( (y[n]-minY)/yStep ) );
    p(xIdx,yIdx) += 1;
  } // for i

  p /= (xLen*xStep*yStep);

} // hist2D

float
entropy1d(
  const prlite::RowVector<float>& p) {

  unsigned len = p.size();
  float sumP = p.sum();
  float h = 0;
  for (unsigned n = 0; n < len; n++) {
    float pn = p[n]/sumP;
    if (pn > 0) { // p log(p) == 0 with p == 0
      h -= pn*log2(pn);
    } // if
  } // for n

  return h;
} // entropy1d

float
entropy2d(
  const prlite::RowMatrix<float>& p) {

  unsigned nRows = p.rows();
  unsigned nCols = p.cols();

  float sumP = 0.0;
  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      sumP += p(r,c);
    } // for c
  } // for r

  float h = 0;
  for (unsigned r = 0; r < nRows; r++) {
    for (unsigned c = 0; c < nCols; c++) {
      float pn = p(r,c)/sumP;
      if (pn > 0) { // p log(p) == 0 with p == 0
        h -= pn*log2(pn);
      } // if
    } // for c
  } // for r

  return h;
} // entropy2d

float
KullbackLeibler(
  const prlite::RowVector<float>& p,
  const prlite::RowVector<float>& q) {

  unsigned sz = p.size();
  PRLITE_ASSERT(sz == unsigned(q.size()), "For discrete Kullback-Leibler divergence the two distributions must have the same length");
  float dist = 0.0;
  for (unsigned n = 0; n < sz; n++) {

    if ( p[n] > std::numeric_limits<float>::epsilon() and
         q[n] > std::numeric_limits<float>::epsilon() ) {
      dist += p[n]*log2(p[n]/q[n]);
    } // if

    else {
      return numeric_limits<float>::infinity();
    } // else

  }// for n

  return dist;
} // KullbackLeibler

float
mutualInformation(
  const prlite::RowVector<float>& px,
  const prlite::RowVector<float>& py,
  const prlite::RowMatrix<float>& pxy,
  float& hx,
  float& hy,
  float& hxy) {

  hx = entropy1d(px);
  hy = entropy1d(py);
  hxy = entropy2d(pxy);

  return  hx + hy - hxy;
} // mutualInformation

float
mutualInformation(
  const prlite::RowVector<float>& x,
  const prlite::RowVector<float>& y,
  float& hx,
  float& hy,
  float& hxy,
  unsigned nXBins,
  unsigned nYBins) {

  prlite::RowVector<float> xc; prlite::RowVector<float> px;
  hist1d(x, xc, px, nXBins);
  prlite::RowVector<float> yc; prlite::RowVector<float> py;
  hist1d(y, yc, py, nYBins);
  prlite::RowMatrix<float> pxy;
  // must have the same number of bins as in hist1d, if the delta
  // widths differ there are variant scaling factors in the integrals.
  hist2d(x, y, xc, yc, pxy, nXBins, nYBins);

  return mutualInformation(px,py,pxy,hx,hy,hxy);
} // mutualInformation
