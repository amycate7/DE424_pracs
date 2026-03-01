#ifndef AFFINEMVG_HPP
#define AFFINEMVG_HPP

/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    April 2022
          PURPOSE: Factor for an affine transformation of a Multivariate Gaussian

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "mvg.hpp" // MVG
#include "sqrtmvg.hpp" // SqrtMVG
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "anytype.hpp"
//#include "oddsandsods.hpp" // ValidValue

// patrec headers
#include "prlite_genmat.hpp"  // Matrix prlite::ColMatrix

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

/******************************************************************************/
/************************ AffineMVG **************************************/
/******************************************************************************/

/**
 * A class describing the joint distribution of the affine
 * transformation y = Ax+c+noise where both x and noise are Gaussian,
 * A is a matrix and c is an offset vector.
 *
 * This also is an exact equivalent of a Gaussian with an
 * unknown mean (represented by the prior over the x variable) and a
 * known precision (represented by the noise distribution).
 *
 * The object starts out in conditional form, being dependent on a
 * message describing the distribution of x. After that
 * multiplication, it becomes a joint over both x and y.
 *
 * NOTE: Because this is a conditional distribution dependent on x, it
 * is undefined until the x-distribution has been multiplied in. This
 * makes it unsuitable to LBP type of inference, instead use LBU.
 *
 * NOTE: The current implementation implements fixed factor operations
 * and does not support variable FactorOperators (which we can later
 * implement at leasure if we so wish).
 */
class AffineMVG : public MVG{

  //============================ Traits and Typedefs ============================

  //======================== Constructors and Destructor ========================

public:

  /**
   * @brief Default constructor.
   */
  //AffineMVG();

  ///
  AffineMVG(const AffineMVG& p) = default;

  /// This one should be ok
  AffineMVG(AffineMVG&& st) = default;

  /**
   * @brief Class Specific Constructor creating a joint affine
   * gaussian with y = Ax+c+noise.
   *
   * @param theXIds: The RV ids of the X MVG. NOTE: To keep things
   * simple we require that all of these MUST be in sorted order AND
   * lower than all of the theYids.
   *
   * @param theYIds: The ids of the resultant Y MVG. To keep things
   * simple we require that all of these must be in sorted order AND
   * must be bigger than theXIds.
   *
   * @param theAMat The matrix X is being left multiplied with. NOTE:
   * not transposed.
   *
   * @param theCVec A fixed offset vector to add to Y. If length is zero,
   * it will be ignored.
   *
   * @param theNoiseR The LHS UPPER Cholesky root of the noise
   * covariance matrix. I.e. Cov_noise = noiseR * noiseR^T
   *
   * @param theObsvIds The (presorted) ids of a subset of Y for which
   * observed values are available. It does not really make sense to
   * have observed X values, you can instead just reduce the A matrix
   * correspondingly and move the fixed bits into the offset c vector.
   *
   * @param theObsvVals Observed values for the previously mentioned
   * subset of Y.
   *
   * NOTE: If necessary we can revisit the strict sorting requirement
   * at a later stage.
   */
  AffineMVG(
    const emdw::RVIds& theXIds,
    const emdw::RVIds& theYIds,
    const prlite::ColMatrix<double>& theAMat,
    const prlite::ColVector<double>& theCVec,
    const prlite::ColMatrix<double>& theNoiseR,
    const emdw::RVIds& theObsvIds = {},
    const emdw::RVVals& theObsvVals = {});

  ///
  virtual ~AffineMVG() {}

  //========================= Operators and Conversions =========================
public:

  ///
  AffineMVG& operator=(const AffineMVG& d) = default;
  AffineMVG& operator=(AffineMVG&& d) = default;

  ///
  bool operator==(const AffineMVG& d) const;
  ///
  bool operator!=(const AffineMVG& d) const {
    return !operator==(d);
  } // operator!=

  //=========================== Inplace Configuration ===========================
public:

  ///
  unsigned classSpecificConfigure(
    const emdw::RVIds& theXIds,
    const emdw::RVIds& theYIds,
    const prlite::ColMatrix<double>& theAMat,
    const prlite::ColVector<double>& theCVec,
    const prlite::ColMatrix<double>& theNoiseR,
    const emdw::RVIds& theObsvIds = {},
    const emdw::RVVals& theObsvVals = {});

  ///
  virtual std::ostream& txtWrite(std::ostream& file) const;
  virtual std::istream& txtRead(std::istream& file);

  //===================== Required Virtual Member Functions =====================
public:

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace, one-for-one AND IN
   * THIS ORDER, the original SORTED set. The sorting restrictions in
   * the class specific constructor, implies that this lot will have
   * to be pre-sorted. We might relax that at a later stage.
   *
   * @param presorted True if newVars is sorted.
   */
  virtual AffineMVG* copy(
    const emdw::RVIds& newVars = {},
    bool presorted = false) const;

  /**
   * make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   *
   * NOTE: selecting a subset of variables is NIY
   */
  virtual AffineMVG* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const {
    return copy();
  } // vacuousCopy


    /// note rhs is a raw pointer
  bool isEqual(const Factor* rhsPtr) const {
    const AffineMVG* dwnPtr = dynamic_cast<const AffineMVG*>(rhsPtr);
    if (!dwnPtr) {return false;}
    return operator==(*dwnPtr);
  } // isEqual

  /// For the while being, use the default implementation of this member
  double distanceFromVacuous() const {return Factor::distanceFromVacuous();}

    ///
  virtual unsigned noOfVars() const;

  /// returns the ids of the variables the factor is defined on
  virtual emdw::RVIds getVars() const;

  ///
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //===================== Required Factor Operations ============================

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from (the top-level) Factor. Typically the members
  /// below would be boilerplate deflecting the request via the
  /// FactorOperators to dynamically select the particular one to use
  /// (DDD). However, we have here opted to instead short circuit
  /// directly to a fixed implementation of the operation.

  //------------------Family 0: eg Sample

  virtual emdw::RVVals sample(Sampler* procPtr) const;

  //------------------Family 1: eg Normalize

  virtual void inplaceNormalize(FactorOperator* procPtr = 0);
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr = 0) const;

  //------------------Family 2: eg Absorb,Cancel

  // NOTE: Updating an AffineGaussian by multiplying in the new x
  // distribution and THEN dividing out the previous x distribution,
  // is equivalent to doing the same multiply and divide to the
  // original (joint) created by SqrtMVG::constructAffineGaussian.
  //
  // However, doing this from the other side only works
  // approximately. I.e. If you similarly multiply in from the right
  // with a y distribution, from that find the updated x-distribution
  // and then rebuild the x,y joint, the resulting (marginalised)
  // y-distribution will not quite match the one you started out with.

  /**
   * @param rhsPtr is the 'other' Factor in this operation. Must be
   * over either x or y vars.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; at the moment this is being ignored.
   */
  virtual void inplaceAbsorb(const Factor* rhsPtr,
                             FactorOperator* procPtr);

  /**
   * @param rhsPtr is the 'other' Factor in this operation. Must be
   * over either x or y vars.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; at the moment this is being ignored.
   *
   * result will inherit all the operators from the lhs factor
   */
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  /**
   * @param rhsPtr is the 'other' Factor in this operation. Must be
   * over either x or y vars.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; at the moment this is being ignored.
   */
  virtual void inplaceCancel(const Factor* rhsPtr,
                             FactorOperator* procPtr);

  /**
   * @param rhsPtr is the 'other' Factor in this operation. Must be
   * over either x or y vars.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; at the moment this is being ignored.
   *
   * result will inherit all the operators from the lhs factor
   */
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr = 0) const;

  //------------------Family 3: eg Marginalize

  /**
   * @brief Marginalises towards selected RVs taking into account
   * obsvY. If x is not instantiated, the marginal will be vacuous.
   */
  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                      bool presorted = false,
                                      const Factor* peekAheadPtr = nullptr,
                                      FactorOperator* procPtr = 0) const;

  //------------------Family 4: eg ObserveAndReduce

  /// NIY.
  virtual uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                           const emdw::RVVals& assignedVals,
                                           bool presorted = false,
                                           FactorOperator* procPtr = 0) const;

  //------------------Family 5: eg inplaceDampen

  /**
   * NIY.  NOTE the parameter is the OLD message and will be weighted
   * with df
   */
  virtual double inplaceDampen(const Factor* oldMsg,
                               double df,
                               FactorOperator* procPtr = 0);

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

  //================================== Friends ==================================

public:

  // /// output
  // friend std::ostream& operator<<( std::ostream& file,
  //                                  const AffineMVG& n );

  // /// input
  // friend std::istream& operator>>( std::istream& file,
  //                                  AffineMVG& n );

  //=============================== Data Members ================================
private:

  // definitive data members

  emdw::RVIds xIds, yIds; // includes observed RVs
  prlite::ColMatrix<double> aMat;
  prlite::ColVector<double> cVec;
  prlite::ColMatrix<double> noiseR;
  emdw::RVIds obsvIds;
  emdw::RVVals obsvVals;

  // secondary data members

  emdw::RVIds allVars; // excluding observed vars

  // This one only depends on messages coming from the x side, never
  // from the y side.
  rcptr<SqrtMVG> xPtr;

  // This joint is the actual distribution being described here. It
  // initially depends on xPtr, but may get changed by messages from
  // the y side. Marginalising (to any side) happens by integrating
  // over this distribution.
  rcptr<SqrtMVG> xyPtr;

  //=============================== DEPRECATED ================================

}; // AffineMVG

#endif // AFFINEMVG_HPP
