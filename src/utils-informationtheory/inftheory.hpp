#ifndef INFTHEORY_HPP
#define INFTHEORY_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    October 2013
          PURPOSE: Various information theory routines

*******************************************************************************/

// standard headers
#include <iostream> // istream, ostream
#include <vector>
#include <set>
#include <map>

// patrec headers
#include "prlite_matrix.hpp"

/**
 * Determines a 1-dim histogram. xc and are the cell centerpoints. p
 * is the distribution height, normalized to integrate (not sum) to
 * unity.
 */
void
hist1d(
  const prlite::RowVector<float>& x,
  prlite::RowVector<float>& xc,
  prlite::RowVector<float>& p,
  unsigned nBins = 10);

/**
 * Determines a 2-dim histogram. xc and yc are the cell
 * centerpoints. p is the distribution height, normalized to integrate
 * (not sum) to unity.
 */
void
hist2d(
  const prlite::RowVector<float>& x,
  const prlite::RowVector<float>& y,
  prlite::RowVector<float>& xc,
  prlite::RowVector<float>& yc,
  prlite::RowMatrix<float>& p,
  unsigned nXBins = 10,
  unsigned nYBins = 10);

/**
 * How many bits of information/uncertainty in the 1-dimensional
 * discrete distribution p.
 *
 * This routine can also be used with p a histogram (i.e. continuous
 * random variable), but the result will differ from the differential
 * entropy: $-V\log V + VH(X)$ with $V$ the size/volume of a
 * bin. Note that the first term shrinks to zero as the bin size
 * shrinks to zero.
 */
float
entropy1d(
  const prlite::RowVector<float>& p);

/**
 * How many bits of information/uncertainty in the 2-dimensional
 * discrete distribution p.
 *
 * This routine can also be used with p a histogram (i.e. continuous
 * random variable), but the result will differ from the differential
 * entropy: $-V\log V + VH(X,Y)$ with $V$ the size/volume of a
 * bin. Note that the first term shrinks to zero as the bin size
 * shrinks to zero.
 */
float
entropy2d(
  const prlite::RowMatrix<float>& p);

/**
 * How many bits of information/uncertainty in the N-dimensional
 * general discrete distribution p.
 *
 * The template parameter T makes this applicable to arbitrary
 * discrete data types such as string, unsigned etc. Think of the map
 * as a generalised vector. Its first parameter (the vector<T> bit)
 * contains the particular outcome for the multiple RVs. The second
 * parameter contains the associated probability for that
 * outcome. This probability may also be unnormalized - this makes it
 * possible to use this routine directly on counts.
 */
template <typename T>
double entropy(
  const std::map<std::vector<T>, double>& p);

/**
 * Calculates a discrete Kullback-Leibler divergence / relative
 * entropy (expressed ito bits). How many extra bits will I need if I
 * encode something originating from p as if its distribution actually
 * was q. \geq 0 with equality only with p == q
 *
 * p and q are two discrete distributions of the same length.
 *
 * NOTE: this divergence is asymmetric and therefore not a true
 * distance measure. It is especially sensitive to cases where q is
 * very small compared to p. If you want to emphasize what is missing
 * (i.e. gaps not in the other) from a distribution, place it as the
 * second parameter. If you want to emphasize what is present (i.e. is
 * it supported by the other), place it as the first parameter. Also,
 * it is often made symmetric by averaging the divergences in both
 * directions.
 */
float
KullbackLeibler(
  const prlite::RowVector<float>& p,
  const prlite::RowVector<float>& q);

/**
 * Calculates how much knowing Y reduces my uncertainty
 * about X (expressed ito bits).
 *
 * @param px, py and pxy are discrete distributions describing x, y
 * and the joint (x,y). NOTE: We can also use this with histograms
 * (i.e. *volume* integrates to unity instead of values directly
 * summing to unity), as long as the same binning is applied to both
 * the scalar and joint variables.
 *
 * @param hx, hy, hxy returns by reference (as a side-effect) the
 * estimated entropies for x, y, and (x,y).
 *
 * @return the mutual information between X and Y.
 */
float
mutualInformation(
  const prlite::RowVector<float>& px,
  const prlite::RowVector<float>& py,
  const prlite::RowMatrix<float>& pxy,
  float& hx,
  float& hy,
  float& hxy);

/**
 * Calculates how much knowing Y reduces my uncertainty about X
 * (expressed ito bits). Their individual as well as joint
 * distributions will be approximated by histograms from which the MI
 * will be calculated.
 *
 * @param x,y are samples of the two continuous random variables.
 *
 * @param hx, hy, hxy returns by reference (as a side-effect) the
 * estimated entropies for x, y, and (x,y).
 *
 * @param nXBins, nYBins are the number of bins for the X and Y
 * histograms respectively (they may differ). If either X or Y are
 * discrete (with equal steps between categories) this should also
 * work -- set the number of bins to the number of discrete classes.
 *
 * @return the mutual information between X and Y.
 */
float
mutualInformation(
  const prlite::RowVector<float>& x,
  const prlite::RowVector<float>& y,
  float& hx,
  float& hy,
  float& hxy,
  unsigned nXBins=10,
  unsigned nYBins=10);

/**
 * @brief Calculates the mutual information between a set of RVs X and another set of RVs Y
 *
 * @param datInt A discrete matrix containing in each row a new
 * observation vector for the set of random variables in
 * question. I.e. each column contains the observations for a specific
 * RV. To allow easy manipulation the union of X and Y need not
 * necessarily cover all the columns, i.e. we can also work with a
 * subset of the RVs in datInt.
 *
 * @param domains the range of (discrete) values each RV can take on
 * @param varsX the indices into datInt containing the X RVs
 * @param varsY the indices into datInt containing the y RVs
 *
 * @param entY returns by reference the estimated joint entropy for
 * the y RVs. This is useful for normalising the MI to determine what
 * percentage of the uncertainty of Y is being accounted for when you
 * know X.
 *
 * @param smoothing The type of Dirichlet smoothing to apply 0: none,
 * 1: partial, 2: full. Option 1 usually is sensible, be very wary of
 * option 2 if you have many random variables or large domains (or
 * both) since the state space very soon becomes prohibitively large.
 *
 * TODO: Maybe we should extend this to calculate MI(X,Y|Z). For this
 *   we'll need a conditioned on set
 *   const std::vector<unsigned> varsZ.
 *   Should be very useful for causality explorations. With varsZ
 *   empty, it reduces to the current setup.
 * If not empty:
 * 1. find the probs for each Z configuration,
 * 2. and then condition all entropies on each of those configurations
 *    before calculating the various entropies,
 * 3. and then calculate the weighted sum of the conditional MIs.
 */
template<typename T>
double estMI(
  const std::vector< std::vector<T> >& datInt,
  const std::vector< std::vector<T> >& domains,
  const std::vector<unsigned>& varsX,
  const std::vector<unsigned>& varsY,
  double& entY,
  unsigned smoothing=1);

/**
 * @brief returns the possible value assignments that each
 * dimension/column of the data matrix can take on. This helper
 * function can calculate the domains needed by estMI (we need them to
 * do efficient Dirichlet smoothing).
 *
 * @param data is a matrix containing in its rows a sequence of (equal
 * length) data-vectors.
 *
 * @return An array of vectors where the i'th vector returns the
 * sorted sequence of values found in the i'th column of the 'data'
 * input matrix.
 */
template <typename T>
std::vector< std::vector<T> > getDomains(const std::vector< std::vector<T> >& data);

#include "inftheory.tcc"

#endif // INFTHEORY_HPP
