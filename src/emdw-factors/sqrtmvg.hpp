#ifndef SQRTMVG_HPP
#define SQRTMVG_HPP

/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Multi-variate SqrtMVG PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "mvg.hpp" // MVG
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "anytype.hpp"
#include "oddsandsods.hpp" // ValidValue

// patrec headers
#include "prlite_genmat.hpp"  // Matrix prlite::ColMatrix

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

class SqrtMVG_InplaceNormalize;
class SqrtMVG_InplaceMaxNormalize;
class SqrtMVG_InplaceAbsorb;
class SqrtMVG_Absorb;
class SqrtMVG_InplaceCancel;
class SqrtMVG_Marginalize;
class SqrtMVG_MaxMarginalize;
class SqrtMVG_ObserveAndReduce;
class SqrtMVG_InplaceWeakDamping;
class SqrtMVG_Sampler;

/******************************************************************************/
/************************ SqrtMVG **************************************/
/******************************************************************************/

/**
 *  A square root form canonical multi-variate Gaussian class.
 *
 * This code relies heavily on LaPack and BLAS. The underlying logic
 * of many routines is available in the devguide.pdf document.
 *
 * With admiration for the simple though profound insight of Wallace
 * Givens (1910 - 1993): Givens, Wallace. "Computation of plane
 * unitary rotations transforming a general matrix to triangular
 * form". J. SIAM 6(1) (1958), pp. 26–50.
 */
class SqrtMVG : public MVG {

  // we'll want to reduce these dependencies later on
  friend class SqrtMVG_InplaceNormalize;
  friend class SqrtMVG_InplaceMaxNormalize;
  friend class SqrtMVG_InplaceAbsorb;
  friend class SqrtMVG_Absorb;
  friend class SqrtMVG_InplaceCancel;
  friend class SqrtMVG_Marginalize;
  friend class SqrtMVG_MaxMarginalize;
  friend class SqrtMVG_ObserveAndReduce;
  friend class SqrtMVG_InplaceWeakDamping;
  friend class SqrtMVG_Sampler;

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:

  /**
   * @brief Default constructor.
   *
   * @param theVars Each variable in the PGM will be identified
   * with a specific integer that uniquely ids it.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  SqrtMVG(
    const emdw::RVIds& theVars = {},
    bool theAutoNormed = true,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  ///
  SqrtMVG(const SqrtMVG& p);

  /// This one should be ok
  SqrtMVG(SqrtMVG&& st) = default;

  /**
    * @brief Class Specific Constructor from square root form parameters.
    *
    * @param theL Lower triangular square root of the K matrix ie
    *   multiplied from the right with its transpose gives K. When
    *   sized 0x0 it indicates a vacuous Gaussian.
    * NOTE: If !presorted we really mean lower triangular, the upper
    *   triangle MUST be zeros.
    *
    * @param theH Canonical form potential vector. Must be same
    *   dimension as theKSqrt.
    *
    * @param theG Canonical from normalization term.
    *
    * @param theAutoNormed When this is set to true the system does
    *   not bother to calculate the g normalization constant - this
    *   implicitly assumes that g is always the value that L and h
    *   implies for a fully normalized system.
    */
  SqrtMVG(
    const emdw::RVIds& theVars,
    const prlite::ColMatrix<double>& theL,
    const prlite::ColVector<double>& theH,
    const double& theG,
    bool theAutoNormed = false,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  ///
  virtual ~SqrtMVG();

  /**
   * @brief For convenience - a covariance form constructor.
   *
   * @param theVars Each variable in the PGM will be identified
   * with a specific integer that uniquely ids it.
   *
   * @param theMn the mean vector.
   *
   * @param theCov the covariance matrix.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   * supplied if available.
   */
  SqrtMVG(
    const emdw::RVIds& theVars,
    const prlite::ColVector<double>& theMn,
    const prlite::RowMatrix<double>& theCov,
    bool theAutoNormed = false,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);


  //========================= Operators and Conversions =========================
public:

  ///
  SqrtMVG& operator=(const SqrtMVG& d) = default;
  SqrtMVG& operator=(SqrtMVG&& d) = default;

  ///
  bool operator==(const SqrtMVG& d) const;

  ///
  bool operator!=(const SqrtMVG& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /// vacuous version
  unsigned classSpecificConfigure(
    const emdw::RVIds& theVars,
    bool theAutoNormed = true,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  /// full version
  unsigned classSpecificConfigure(
    const emdw::RVIds& theVars,
    const prlite::ColMatrix<double>& theL,
    const prlite::ColVector<double>& theH,
    const double& theG,
    bool theAutoNormed = false,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  ///
  virtual std::istream& txtRead(std::istream& file);
  ///
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from (the top-level) Factor. Typically the members
  /// below would be boilerplate deflecting the request via the
  /// FactorOperators to dynamically select the particular one to use
  /// (DDD).

  //------------------Family 0: eg Sample

  virtual emdw::RVVals sample(Sampler* procPtr) const;

  //------------------Family 1: eg Normalize

  virtual void inplaceNormalize(FactorOperator* procPtr = 0);
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr = 0) const;

  //------------------Family 2: eg Absorb,Cancel

  virtual void inplaceAbsorb(const Factor* rhsPtr,
                            FactorOperator* procPtr = 0);
  /// result will inherit all the operators from the lhs factor
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  virtual void inplaceCancel(const Factor* rhsPtr,
                            FactorOperator* procPtr = 0);
  /// result will inherit all the operators from the lhs factor
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  //------------------Family 3 eg Marginalize

  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                      bool presorted = false,
                                      const Factor* peekAheadPtr = nullptr,
                                      FactorOperator* procPtr = 0) const;

  //------------------Family 4: eg ObserveAndReduce

  virtual uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                           const emdw::RVVals& assignedVals,
                                           bool presorted = false,
                                           FactorOperator* procPtr = 0) const;

  //------------------Family 5: eg inplaceDampen

  /**
   * NOTE the parameter is the OLD message and will be weighted with df
   */
  virtual double inplaceDampen(const Factor* oldMsg,
                               double df,
                               FactorOperator* procPtr = 0);

  //===================== Required Virtual Member Functions =====================
public:

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace, one-for-one AND IN
   * THIS ORDER, the original SORTED set.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual SqrtMVG* copy(
    const emdw::RVIds& newVars = {},
    bool presorted = false) const;

  /**
   * make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   *
   * NOTE: selecting a subset of variables is NIY
   */
  virtual SqrtMVG* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /// note rhs is a raw pointer
  bool isEqual(const Factor* rhsPtr) const;

  /// For the while being, use the default implementation of this member
  double distanceFromVacuous() const {return Factor::distanceFromVacuous();}

  ///
  virtual unsigned noOfVars() const;

  /// returns the ids of the variables the factor is defined on
  virtual emdw::RVIds getVars() const;

  ///
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //====================== Other Virtual Member Functions =======================
public:

  /// returns the factor potential at a given variable assignment
  virtual double potentialAt(const emdw::RVIds& theVars,
                             const emdw::RVVals& theirVals,
                             bool presorted = false) const;

  /**
   * Measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference. Currently based
   * on the Mahanalobis distance.
   *
   * We might want to refine this further, getting it right has a huge
   * influence on the order of the message passing in loopy belief
   * propagation which, in its turn, has a huge influence on
   * convergence speed.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Parameterizations (non-virtual) =====================
public:

  /**
   * @return whether this is a vacuous SqrtMVG
   */
  bool isVacuous() const;

  /**
   * @return the lower triangular root of the precision matrix. Note:
   *   the strictly upper triangle is undefined and may contain
   *   arbitrary numbers.
   */
    prlite::ColMatrix<double> getL() const;

  /// returns the current potential vector
  prlite::ColVector<double> getH() const;

  /// returns the current g value
  double getG() const;

  /// column 0 contains the mean, the other column i = 1 .. dim
  /// columns contains sqrt(dim+0.5) * (L^{-T}.col(i))
  prlite::ColMatrix<double> getSigmaPoints() const;

  ///
  prlite::ColVector<double> getMean() const;

  /// Hopefully mostly for legacy reasons. Also why we return a Matrix
  /// and not a prlite::ColMatrix
  prlite::RowMatrix<double> getCov() const;

  /// Hopefully mostly for legacy reasons. Also why we return a Matrix
  /// and not a prlite::ColMatrix
  prlite::RowMatrix<double> getK() const;

  /**
   * @brief Do an ML estimate of the mean and LHS UPPER Cholesky root
   * R of the covariance matrix from a sequence of samples.
   *
   * @param samples Each column contains a sample vector.
   *
   * @param mn The mean value estimated from samples.
   *
   * @param R The LHS UPPER Cholesky root of the estimated covariance
   * matrix i.e. cov = R*R^T with R an upper triangular matrix. This
   * unusual choice is intentional, L = inv(R^T) will be the lower
   * triangular root of the precision matrix.
   *
   * @param biased If true, normalization is with N-1 (i.e. sample
   * based). If false normalization is with N (i.e. population based).
   *
   * @return success
   */
  static bool estMnR(
    const prlite::ColMatrix<double>& samples,
    prlite::ColVector<double>& mn,
    prlite::ColMatrix<double>& R,
    bool biased = true);

  /**
   * @brief Do an ML estimate of a SqrtMVG from a sequece of samples.
   *
   * @param samples Each column contains a sample vector.
   *
   * @param biased If true, normalization is with N-1 (i.e. sample
   * based). If false normalization is with N (i.e. population based).
   */
  static  SqrtMVG* constructFromSamples(
    const emdw::RVIds& theVars,
    const prlite::ColMatrix<double>& samples,
    bool biased = true,
    bool theAutoNormed = false,
    bool presorted = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  /**
   * @brief Builds a joint Gaussian for (X,Y) where Y = g(X) +
   * noise. g(X) can be a general non-linear function but also
   * includes the so-called linear gaussian joint as a special case.
   *
   * @param xVars The RVIds of the X variables.
   * @param xSPts The sigma points of X, each in a column.
   *
   * @param yVars The RVIds of the Y variables.
   * @param ySPts The Y values that xSPts map to, each in a column.
   *
   * @param noiseR The LHS UPPER Cholesky root of the noise
   * covariance matrix. I.e. Cov_noise = noiseR * noiseR^T
   * NOTE: dpotrf_ in lapackblas.hpp returns Cov_noise = U^T*U thus
   * noiseR = U^T
   *
   * @param theAutoNormed As always, will the resultant SqrtMVG be
   * implicitly normalized or not.
   *
   * @param presortedX Are xVars sorted from lowest to highest. Save
   * yourself a lot of computation by ensuring that they are.
   *
   * @param presortedY Are yVars sorted from lowest to
   * highest. Save yourself a lot of computation by ensuring that they
   * are, and also that all of them have higher ids than the X vars.
   */
  static  SqrtMVG* constructFromSigmaPoints(
    const emdw::RVIds& xVars,
    const prlite::ColMatrix<double>& xSPts,
    const emdw::RVIds& yVars,
    const prlite::ColMatrix<double>& ySPts,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed = false,
    bool presortedX = false,
    bool presortedY = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  /**
   * @brief Linear Gaussians: we have an existing RV X with a Gaussian
   *   Factor p(X); we create an RV Y = AX + c + N where A is a linear
   *   transformation matrix, c is an offset vector and N is zero mean
   *   Gaussian noise with covariance Q; this constructor builds the
   *   joint Gaussian for p(X,Y).
   *
   * @param xPdf is an existing SqrtMVG describing RV X.
   *
   * @param aMat The matrix X is being left multiplied with. Y =
   * AX +c. NOTE: not transposed.
   *
   * @param cVec A fixed offset vector to add to Y. If length is zero,
   * it will be ignored.
   *
   * @param newVars The ids of the additional variables Y (disjoint
   *   from those already in xPdf). newVars.size() == aMat.rows()
   *   is the dimension of Y. Added to that of X it gives the total
   *   dimension of the joint.
   *
   * @param noiseR The LHS UPPER Cholesky root of the noise
   * covariance matrix. I.e. Cov_noise = noiseR * noiseR^T
   *
   * @param theAutoNormed As always, will the resultant SqrtMVG be
   * implicitly normalized or not.
   *
   * @param presortedY Set this to true if newVars are sorted from
   *   smallest to biggest and avoid sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   *   supplied if available. Take care to specify them in their
   *   proper order.
   */
  static  SqrtMVG* constructAffineGaussian(
    const SqrtMVG& xPdf,
    const prlite::ColMatrix<double>& aMat,
    const prlite::ColVector<double>& cVec,
    const emdw::RVIds& newVars,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed = false,
    bool presortedY = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  static  SqrtMVG* _constructAffineGaussian(
    const SqrtMVG& xPdf,
    const prlite::ColMatrix<double>& aMat,
    const prlite::ColVector<double>& cVec,
    const emdw::RVIds& newVars,
    const prlite::ColMatrix<double>& noiseR,
    bool theAutoNormed = false,
    bool presortedY = false,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0,
    rcptr<Sampler> theSampler = 0);

  //======================= Non-virtual Member Functions ========================
public:

  /// returns the g value that will normalize the MVG to unity volume
  static double getNormedG(const prlite::ColMatrix<double>& theL, const prlite::ColVector<double>& theH);

  /// returns the current unnormalized hyper-volume of the MVG
  double getMass() const;

  // Add this value to the existing volume of the MVG
  void adjustMass(double newMass);

  // change this pdf into a vacuous one
  void makeVacuous();

  ///
  double mahalanobis(const prlite::ColVector<double>& observedVals) const;

  /// Directly calculates the symmetric KL-divergence between two MVG's as  0.5*(KL(P||Q) + KL(Q||P))
  static double symmKLdivergence(const Factor* lhsPtr, const Factor* rhsPtr);

  /**
   * @param xId, yId: The idno's for the two random variables for
   * which the pdf must be shown.
   *
   * Writes a text file containing:
   * 1st line: the values of x where the density was evaluated
   * 2nd line: the values of y where the density was evaluated
   * Further lines: A mesh of p(x,y) with x running horizontally and y
   * running vertically (or vice versa if xId > yId)
   *
   * If xId == yId the one-dimensional result will be repeated across
   * all rows, i.e. you can also inspect one-dimensional results with
   * this function.
   */
  void
  export2DMesh(
    const std::string& filename,
    emdw::RVIdType xId, emdw::RVIdType yId,
    unsigned nPoints) const;

  //================================== Friends ==================================

public:

  /// output
  friend std::ostream& operator<<( std::ostream& file,
                                   const SqrtMVG& n );

  /// input
  friend std::istream& operator>>( std::istream& file,
                                   SqrtMVG& n );

  //=============================== Data Members ================================
private:

  emdw::RVIds vars;

  // vacuous==true:  Conceptually L=h=g=0, but physically undefined, or
  // vacuous==false: L (guaranteed posdef), h and g is properly instantiated.
  bool vacuous;

  // only defined with vacuous==false
  mutable prlite::ColMatrix<double> L; // L is the lower triangular sqrt of K. upper part undefined
  mutable prlite::ColVector<double> h; // potential vector
  mutable double g;            // normalization constant only defined with autoNormed==false
  bool autoNormed;             // no need to calculate g then

  // The set of FactorOperators in use
  rcptr<Marginalizer> marginalizer;
  rcptr<InplaceNormalizer> inplaceNormalizer;
  rcptr<Normalizer> normalizer;
  rcptr<InplaceAbsorber> inplaceAbsorber;
  rcptr<Absorber> absorber;
  rcptr<InplaceCanceller> inplaceCanceller;
  rcptr<Canceller> canceller;
  rcptr<ObserveAndReducer> observeAndReducer;
  rcptr<InplaceDamper> inplaceDamper;
  rcptr<Sampler> sampler;

  //=============================== DEPRECATED ================================

}; // SqrtMVG

// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a SqrtMVG.
 */
class SqrtMVG_Sampler: public SampleOperator< SqrtMVG > {
public:
  const std::string& isA() const;
  ///
  emdw::RVVals sample(
    const SqrtMVG* lhsPtr);
  /// Each column contains an independent sample.
  static prlite::ColMatrix<double> sample(
    const SqrtMVG* lstPtr,
    size_t nSamples);
private:
  /// A helper function doing the actual sampling, but without adding
  /// the mean. Each column contains an independent sample.
  static void sampleCentered(
     const SqrtMVG* lhsPtr,
     prlite::ColMatrix<double>& x,
     size_t nSamples=1);

}; // SqrtMVG_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/
// This operator does an inplace normalization of a SqrtMVG.

class SqrtMVG_InplaceNormalize: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  void inplaceProcess(SqrtMVG* lstPtr);
}; // SqrtMVG_InplaceNormalize

// This operator returns a normalized version of the input factor.

class SqrtMVG_Normalize: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  Factor* process(const SqrtMVG* lstPtr);
}; // SqrtMVG_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

// This operator does an inplace factor product with another
// SqrtMVG.

class SqrtMVG_InplaceAbsorb: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;

  /**
   * This operator inplace absorbs another SqrtMVG message, ie
   * does an inplace factor product.
   */
  void inplaceProcess(SqrtMVG* lhsPtr,
                      const Factor* rhsFPtr){
    weightedInplaceAbsorb(lhsPtr, rhsFPtr);
  } // inplaceProcess

  /**
   * Same as inplaceProcess, but this one puts a relative weighting
   * on rhsFPtr (defaulted to 1.0).
   */
  void weightedInplaceAbsorb(SqrtMVG* lhsPtr,
                             const Factor* rhsFPtr,
                             const ValidValue<double>& weight = ValidValue<double>() );
}; // SqrtMVG_InplaceAbsorb

// This operator returns a factor product with another
// SqrtMVG.

class SqrtMVG_Absorb: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  Factor* process(const SqrtMVG* lhsPtr,
                  const Factor* rhsFPtr);
}; // SqrtMVG_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

// This operator does an inplace factor division with another
// SqrtMVG.

class SqrtMVG_InplaceCancel: public Operator1< SqrtMVG > {
public:

  const std::string& isA() const;

/**
 * This operator inplace cancels another SqrtMVG message, ie does an
 * inplace factor division. It is assumed that we will only cancel
 * messages that were previously absorbed. If not our result might go
 * negative definite.
 */
  void inplaceProcess(SqrtMVG* lhsPtr,
                      const Factor* rhsFPtr){
    weightedInplaceCancel(lhsPtr, rhsFPtr);
  } // inplaceProcess

  /**
   * Same as inplaceProcess, but this one puts a relative weighting
   * on rhsFPtr (defaulted to 1.0).
   */
  void weightedInplaceCancel(SqrtMVG* lhsPtr,
                             const Factor* rhsFPtr,
                             const ValidValue<double>& weight = ValidValue<double>() );
}; // SqrtMVG_InplaceCancel

// This operator returns a factor division with another
// SqrtMVG.

class SqrtMVG_Cancel: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  Factor* process(const SqrtMVG* lhsPtr,
                  const Factor* rhsFPtr);
}; // SqrtMVG_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

/** This operator returns a marginalized factor, making use of the
 * covariance form representation. This version is particularly
 * efficient when combined with the Lauritzen-Spiegelhalter algorithm.
 */

class SqrtMVG_Marginalize: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  Factor* process(const SqrtMVG* lhsPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // SqrtMVG_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

/**
 * After observing some variables, the scope of the resultant factor
 * reduces to ONLY that of the unobserved variables, making further
 * computations potentially much more efficient.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 */

class SqrtMVG_ObserveAndReduce: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  Factor* process(const SqrtMVG* lhsPtr,
                const emdw::RVIds& variables,
                const emdw::RVVals& assignedVals,
                bool presorted);
}; // SqrtMVG_ObserveAndReduce

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/

/**
 * The true result of doing damping on a Gaussian, is a 2-component
 * GMM. This operator replaces that 2-component GMM with a single
 * Gaussian with corresponding first and second order moments. The
 * special case where df == 0 is also specifically handled
 * efficiently. The returned distance is determined by using the
 * standard 'distance' member on the two undamped distributions.
 *
 * new = lhs*(1-df) + rhs*df
 */

// class SqrtMVG_InplaceWeakDamping: public Operator1< SqrtMVG > {
// public:
//   const std::string& isA() const;
//   double inplaceProcess(SqrtMVG* lhsPtr,
//                         const Factor* rhsPtr,
//                         double df);
// }; // SqrtMVG_InplaceWeakDamping

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/

/**
 * DEPRECATED The true result of doing damping on a Gaussian is a
 * 2-component GMM. Unfortunately not inplace anymore. Instead, this
 * one does no damping at all. The returned distance is determined by
 * using the standard 'distance' member on the two undamped
 * distributions. This version is largely replaced by the
 * SqrtMVG_InplaceWeakDamping operator.
 */

class SqrtMVG_InplaceNoDamping: public Operator1< SqrtMVG > {
public:
  const std::string& isA() const;
  double inplaceProcess(SqrtMVG* lhsPtr,
                        const Factor* rhsPtr,
                        double df);
}; // SqrtMVG_InplaceNoDamping

#endif // SQRTMVG_HPP
