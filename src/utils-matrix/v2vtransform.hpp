#ifndef V2VTRANSFORM_HPP
#define V2VTRANSFORM_HPP

/*******************************************************************************

  AUTHOR: JA du Preez
  DATE: May 2014

  PURPOSE: Implements non-linear functors that transforms an
  unknown/random ColVector<double> into another one. This defines the
  non-linearity at the heart of sigma point / UKF systems. NOTE: Some
  of the transforms here are linear, but are included here for the
  convenience of a common framework.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// patrec headers
#include "prlite_genvec.hpp"
#include "prlite_genmat.hpp"

#include "emdw.hpp"

/**
 * Abstract class for functor that transforms a ColVector into
 * typically one (but could be more) ColVectors.
 */
class V2VTransform {
public:
  V2VTransform() {}
  virtual ~V2VTransform() {}
  virtual std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const = 0;

}; //V2VTransform

// ---------------------- Exponentiation

/**
 * Inverts each element of the vector (individually).
 *
 * WARNING: watch out for division by zero!
 */
class Invert : public V2VTransform{
public:
  Invert() {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const;
}; // Invert

/**
 * Exponentiates each element of the vector (individually).
 *
 * WARNING: watch out for roots of negative numbers!
 */
class Exponentiate : public V2VTransform{
public:
  Exponentiate(double thePwr) : pwr(thePwr) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const;
private:
  double pwr;
}; // Exponentiate

// ---------------------- Distances between points

/**
 * Calculates the Euclidean distance from a random/unknown vector to
 * some known reference point.
 *
 * This is useful for estimating an UKF-style joint Gaussian
 * describing the Euclidean distance between a Gaussian random vector
 * and a known point.
 *
 * NOTE 1: We can of course extend this to other norms too.
 *
 * NOTE 2: ONLY use this in UKF-style gaussians if you have informed
 * priors for those gaussians. Using it on standard zero-mean -
 * whatever-variance gaussians is a particularly bad idea since the
 * sigma point projections then operate in a very non-linear domain
 * giving senseless results.
 */
class DistanceFromKnownPoint : public V2VTransform{
public:

  /**
   * @param theX0 The known reference point.
   */
  DistanceFromKnownPoint(const prlite::ColVector<double>& theX0) : x0(theX0) {}

  /**
   * @return A vector of size 1 containing the distance between the
   * point x and the known reference point.
   *
   * @param x The point to determine the distance to. NOTE: When using
   * this function in a non-linear Gaussian, it is very important that
   * the variable ids in that Gaussian are *presorted*. Otherwise the
   * internal sorting in that Gaussian will jumble the dimensions into
   * a different order, making comparision with the original reference
   * point meaningless.
   */
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const;
private:
  prlite::ColVector<double> x0;
}; // DistanceFromKnownPoint

/**
 * Calculates the Euclidean distance from a random/unknown vector to
 * several known reference points.
 *
 * This is useful for estimating an UKF-style joint Gaussian
 * describing the Euclidean distance between a Gaussian random vector
 * and several known points.
 *
 * NOTE 1: We can of course extend this to other norms too.
 *
 * NOTE 2: ONLY use this in UKF-style gaussians if you have informed
 * priors for those gaussians. Using it on standard zero-mean -
 * whatever-variance gaussians is a particularly bad idea since the
 * sigma point projections then operate in a very non-linear domain
 * giving senseless results.
 */
class DistancesFromKnownPoints : public V2VTransform{
public:

  /**
   * @param theXs The known reference points.
   */
  DistancesFromKnownPoints(const std::vector< prlite::ColVector<double> >& theXs) : xs(theXs) {}

  /**
   * @return A vector containing the distances between the point x and
   * the known reference points.
   *
   * @param x The point to determine the distances to. NOTE: When
   * using this function in a non-linear Gaussian, it is very
   * important that the variable ids in that Gaussian are
   * *presorted*. Otherwise the internal sorting in that Gaussian will
   * jumble the dimensions into a different order, making comparision
   * with the original reference point meaningless.
   */
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const;
private:
  std::vector< prlite::ColVector<double> > xs;
}; // DistancesFromKnownPoint

/**
 * Calculates the Euclidean distance between two random/unknown
 * vectors.
 *
 * This is useful for estimating an UKF-style joint Gaussian
 * describing the Euclidean distance between two Gaussian random
 * vectors.
 *
 * NOTE 1: We can of course extend this to other norms too.
 *
 * NOTE 2: ONLY use this in UKF-style gaussians if you have informed
 * priors for those gaussians. Using it on standard zero-mean -
 * whatever-variance gaussians is a particularly bad idea since the
 * sigma point projections then operate in a very non-linear part of
 * the transformation ending in senseless results.
 */
class DistanceBetweenUnknownPoints : public V2VTransform{
public:

  /**
   * The constructor.
   *
   * @param rvIds1 The (sorted) EMDW random variable ids for the first
   * position vector. These are being supplied because when the joint
   * positional Gaussian is built from its constituent two component
   * positional Gaussians, the order of the RVs are being sorted and
   * thus may change their position in te vector. Combined with the
   * next parameter (rvIds2) these RVids are being used to keep track
   * of where in the x vector the original positional coordinates now
   * live. rvIds may even share some RVs with rvIds2. NOTE: We assume
   * that these coordinates are already sorted i.e. exactly what
   * calling the getVars() member on the Gaussian for positional
   * distribution 1 will provide.
   *
   * @param rvIds2 The (sorted) EMDW random variable ids for the
   * second position vector. These are being supplied because when the
   * joint positional Gaussian is built from its constituent two
   * component positional Gaussians, the order of the RVs are being
   * sorted and thus may change their position in te vector. Combined
   * with the previous parameter (rvIds1) these RVids are being used
   * to keep track of where in the x vector the original positional
   * coordinates now live. rvIds may even share some RVs with
   * rvIds1. NOTE: We assume that these coordinates are already sorted
   * i.e. exactly what calling the getVars() member on the Gaussian
   * for positional distribution 2 will provide.
   */
  DistanceBetweenUnknownPoints(
    const emdw::RVIds& rvIds1,
    const emdw::RVIds& rvIds2);

  /**
   * @param x The vector containing the current values of the two
   * positions to determine the distance between. The two positions
   * may even share some coordinates, i.e. this vector has a dimension
   * <= 2*dimensionality of space. Typically x is a sigma-point of the
   * joint Gaussian that describes how these two positions are jointly
   * distributed.
   *
   */
  std::vector< prlite::ColVector<double> >
  operator()(
    const prlite::ColVector<double>& x) const;

private:
  // where in x can one find the coordinates of the two vectors
  std::vector<size_t> v1Idx;
  std::vector<size_t> v2Idx;
}; // DistanceBetweenUnknownPoints

// ---------------------- Matrix vector products Ax

/**
 * y = Ax with A a known, pre-specified matrix and x
 * variable. I.e. the type of thing we will use in a linear Gaussian.
 */
class MpyKnownAMatWithUnknownXVec : public V2VTransform{
public:
  MpyKnownAMatWithUnknownXVec(const prlite::RowMatrix<double>& theAMat) : aMat(theAMat) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& x) const;
private:
  prlite::RowMatrix<double> aMat;
}; // MpyKnownAMatWithUnknownXVec

/**
 * y = Ax with x a known, pre-specified vector and A
 * variable. I.e. the type of thing we will use in a linear Gaussian.
 */
class MpyUnknownAMatWithKnownXVec : public V2VTransform{
public:
  MpyUnknownAMatWithKnownXVec(const prlite::ColVector<double>& theX, unsigned theOutDim) : x(theX), inDim(x.size()), outDim(theOutDim) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& a) const;
private:
  prlite::ColVector<double> x;
  unsigned inDim, outDim;
}; // MpyUnknownAMatWithKnownXVec

/**
 * y = Ax, but now both A and x are variables; we assume that x is the
 * first inDim variables of xa and that the further outDim by inDim
 * variables are the (row-dense) matrix A.
 *
 * The difference with MpyKnownAMatWithUnknownXVec is that A now also forms part of the
 * unknown input.
 *
 * NOTE: Use this only if you have interesting pdfs for *both* A and
 * x, otherwise this will probably not do anything useful.
 */
class MpyUnknownAMatWithUnknownXVec : public V2VTransform{
public:
  MpyUnknownAMatWithUnknownXVec(unsigned theInDim, unsigned theOutDim) : inDim(theInDim), outDim(theOutDim) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& xa) const;
private:
  unsigned inDim, outDim;
}; // MpyUnknownAMatWithUnknownXVec

// ---------------------- Gaussian priors

/**
 * Works well! The purpose of this is to provide 2D+1 predicted
 * samples/vectors of a D-dimensional random vector y given that we
 * know its covariance <cov> and the client code provides an example
 * of what its mean <mn> might be.
 */
class UnknownMeanKnownCov : public V2VTransform{
public:
  UnknownMeanKnownCov(const prlite::RowMatrix<double>& theCov) : cov(theCov) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& mn) const;
private:
  prlite::RowMatrix<double> cov;
}; // UnknownMeanKnownCov

/**
 * NOTE: DONT BOTHER with this, not useful. It turns out that the
 * predicted samples are independent of the covariances - therefore
 * one can not infer estimates for those covariances.
 *
 * The purpose of this is to provide 2D+1 predicted samples/vectors of
 * a D-dimensional random vector y given that we know its mean <mn>
 * and the client code provides an example of what its covariance,
 * packed in <covAsVec>, might be.
 *
 * <covAsVec> packs the lower triangle of the covariance matrix in a
 * vector, starting from the main diagonal and then working its way
 * down through the lower diagonals.
 */
class KnownMeanUnknownCov : public V2VTransform{
public:
  KnownMeanUnknownCov(const prlite::ColVector<double>& theMn) : mn(theMn) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& covAsVec) const;
private:
  prlite::ColVector<double> mn;
}; // UnknownMeanKnownCov


/**
 * NOTE: DONT BOTHER with this, not useful. The purpose of this is to
 * provide 2D+1 predicted samples/vectors of an unknown D-dimensional
 * MEAN_VECTOR <mnY> given that we know its own mean <mn> and the
 * client code provides an example of what cov(Y), packed in
 * <covYAsVec>, might be. Similar to how it is done in a Gauss-Wishart
 * distribution, it is assumed that the cov(mnY) is a scaled down
 * version of cov(Y) namely cov(mnY) = cov(Y)/beta (the scaling factor
 * being <beta> > 1). This of course means that the variances gets
 * smaller, and correspondingly if we model them via Gaussians, they
 * should initially be big enough to avoid problems with positive
 * definiteness.
 *
 * <covYAsVec> packs the lower triangle of the cov(Y) in a vector,
 * starting from the main diagonal and then working its way down
 * through the lower diagonals.
 */
class KnownMeanScaledUnknownCov : public V2VTransform{
public:
  KnownMeanScaledUnknownCov(const prlite::ColVector<double>& theMn, double theBeta)
      : mn(theMn), beta(theBeta) {}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& covYAsVec) const;
private:
  prlite::ColVector<double> mn;
  double beta;
}; // UnknownMeanKnownCov


// NOTE: DONT BOTHER with this, not useful.
class UnknownMeanUnknownCov : public V2VTransform{
public:
  UnknownMeanUnknownCov(unsigned theDim)
      : dim(theDim){}
  std::vector< prlite::ColVector<double> > operator()(const prlite::ColVector<double>& covYAndMnYAsVec) const;
private:
  unsigned dim;
}; // UnknownMeanUnknownCov


#endif // V2VTRANSFORM_HPP
