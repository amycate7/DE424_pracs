#ifndef TREEFACTOR_HPP
#define TREEFACTOR_HPP

/*******************************************************************************

          AUTHORS: JA du Preez
          DATE:    September 2017
          PURPOSE: Multi-variate TreeFactor PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "clustergraph.hpp"


class TreeFactor_InplaceNormalize;
class TreeFactor_InplaceMaxNormalize;
class TreeFactor_InplaceAbsorb;
class TreeFactor_Absorb;
class TreeFactor_InplaceCancel;
class TreeFactor_Marginalize;
class TreeFactor_MaxMarginalize;
class TreeFactor_ObserveAndReduce;
class TreeFactor_InplaceDamping;
class TreeFactor_Sampler;

/******************************************************************************/
/************************ TreeFactor **************************************/
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
class TreeFactor : public Factor {

  // we'll want to reduce these dependencies later on
  friend class TreeFactor_InplaceNormalize;
  friend class TreeFactor_InplaceMaxNormalize;
  friend class TreeFactor_InplaceAbsorb;
  friend class TreeFactor_Absorb;
  friend class TreeFactor_InplaceCancel;
  friend class TreeFactor_Marginalize;
  friend class TreeFactor_MaxMarginalize;
  friend class TreeFactor_ObserveAndReduce;
  friend class TreeFactor_InplaceDamping;
  friend class TreeFactor_Sampler;

  //============================ Traits and Typedefs ============================
private:
  typedef std::vector<double> DASS;

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
  TreeFactor(
    const rcptr <ClusterGraph>& theGraphPtr = 0,
    bool preCalibrated = true,
    const rcptr<Marginalizer>& theMarginalizer = 0,
    const rcptr<InplaceNormalizer>& theInplaceNormalizer = 0,
    const rcptr<Normalizer>& theNormalizer = 0,
    const rcptr<InplaceAbsorber>& theInplaceAbsorber = 0,
    const rcptr<Absorber>& theAbsorber = 0,
    const rcptr<InplaceCanceller>& theInplaceCanceller = 0,
    const rcptr<Canceller>& theCanceller = 0,
    const rcptr<ObserveAndReducer>& theObserveAndReducer = 0,
    const rcptr<InplaceDamper>& theInplaceDamper = 0,
    const rcptr<Sampler>& theSampler = 0);


  TreeFactor(
      const std::vector<rcptr<Factor>>& theFactorPtrs,
      bool preCalibrated = true,
      const rcptr<Marginalizer>& theMarginalizer = 0,
      const rcptr<InplaceNormalizer>& theInplaceNormalizer = 0,
      const rcptr<Normalizer>& theNormalizer = 0,
      const rcptr<InplaceAbsorber>& theInplaceAbsorber = 0,
      const rcptr<Absorber>& theAbsorber = 0,
      const rcptr<InplaceCanceller>& theInplaceCanceller = 0,
      const rcptr<Canceller>& theCanceller = 0,
      const rcptr<ObserveAndReducer>& theObserveAndReducer = 0,
      const rcptr<InplaceDamper>& theInplaceDamper = 0,
      const rcptr<Sampler>& theSampler = 0);

  ///
  TreeFactor(const TreeFactor& p);

  /// This one should be ok
  TreeFactor(TreeFactor&& st) = default;

  ///
  virtual ~TreeFactor();

  //========================= Operators and Conversions =========================
public:

  ///
  TreeFactor& operator=(const TreeFactor& d) = default;
  TreeFactor& operator=(TreeFactor&& d) = default;

  ///
  bool operator==(const TreeFactor& d) const;

  ///
  bool operator!=(const TreeFactor& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /// vacuous version
  unsigned classSpecificConfigure(
    const rcptr<ClusterGraph>& theGraphPtr = 0,
    const rcptr<Marginalizer>& theMarginalizer = 0,
    const rcptr<InplaceNormalizer>& theInplaceNormalizer = 0,
    const rcptr<Normalizer>& theNormalizer = 0,
    const rcptr<InplaceAbsorber>& theInplaceAbsorber = 0,
    const rcptr<Absorber>& theAbsorber = 0,
    const rcptr<InplaceCanceller>& theInplaceCanceller = 0,
    const rcptr<Canceller>& theCanceller = 0,
    const rcptr<ObserveAndReducer>& theObserveAndReducer = 0,
    const rcptr<InplaceDamper>& theInplaceDamper = 0,
    const rcptr<Sampler>& theSampler = 0);

  ///
  virtual std::istream& txtRead(std::istream& file);
  ///
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from (the top-level) Factor. They return uniqptr so
  /// that we can chain them.

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

  //------------------Family 4: eg Observe, ObserveAndReduce

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
  virtual TreeFactor* copy(
    const emdw::RVIds& newVars = {},
    bool presorted = false) const;

  /**
   * make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   *
   * NOTE: selecting a subset of variables is NIY
   */
  virtual TreeFactor* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /// note rhs is a raw pointer
  bool isEqual(const Factor* rhsPtr) const;

  /// For the while being, use the default implementation of this member
  virtual double distanceFromVacuous() const; //{return Factor::distanceFromVacuous();}

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

  //======================= Non-virtual Member Functions ========================
public:

  void calibrateteWithMessagePassing();

  std::vector<std::pair<size_t, size_t>> getFactorLineup(const rcptr<ClusterGraph> &lhsGPtr,
                                                         const rcptr<ClusterGraph> &rhsGPtr);

  std::vector<std::pair<emdw::RVIds,std::vector<unsigned>>>
    getMarginalizationSequence(const rcptr<ClusterGraph> &lhsGPtr);


  //================================== Friends ==================================

public:

  /// output
  template<typename T>
  friend std::ostream& operator<<( std::ostream& file,
                                   const TreeFactor& n );

  /// input
  template<typename T>
  friend std::istream& operator>>( std::istream& file,
                                   TreeFactor& n );

  //=============================== Data Members ================================
public:
  rcptr<ClusterGraph> graphPtr;

private:
  emdw::RVIds vars;

public:
  //Cached line-up information
  static std::map<uint64_t, //unique hash generated from lhs rhs graphs
                  std::vector<std::pair<size_t, size_t>>> factorsLineupCache;


  //Cached line-up information
  static std::map<uint64_t, //unique hash generated from lhs graph and toKeep
                  std::vector<std::pair<emdw::RVIds, std::vector<unsigned>>>
                 > marginalizationSequenceCache;

private:
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


}; // TreeFactor


//=============================== Helpers ==============================

/**
 * Test if a factor is a valid treeFactor type and return true/false
 *
 * @param fPtr. a rcptr<Factor> pointer
 */
bool isTFPtr(const rcptr<Factor> &fPtr);

/**
 * Downcast a given rcptr<Factor> to a rcptr<TreeFactor> in order
 * to access all child members.
 *
 * @param fPtr. a Factor pointer
 */
rcptr<TreeFactor> toTFPtr(const rcptr<Factor> &fPtr);

/**
 * Flatten a vector of rcptr<Factor> so that any level of TreeFactor
 * graphs are extracted into a single vector of non-TreeFactors.
 *
 * @param factorPtrs. vector of  factor pointers
 */
std::vector<rcptr<Factor>> flattenTFPtrs(
    const std::vector<rcptr<Factor>> &factorPtrs);


// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a TreeFactor.
 */
class TreeFactor_Sampler: public SampleOperator< TreeFactor > {
public:
  const std::string& isA() const;
  ///
  emdw::RVVals sample(
    const TreeFactor* lhsPtr);
}; // TreeFactor_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** TreeFactor_Normalize ************************/
/*****************************************************************************/
// This operator does an inplace normalization of a TreeFactor.

class TreeFactor_InplaceNormalize: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(TreeFactor* lstPtr);
}; // TreeFactor_InplaceNormalize

// This operator returns a normalized version of the input factor.

class TreeFactor_Normalize: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  Factor* process(const TreeFactor* lstPtr);
}; // TreeFactor_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** TreeFactor_Absorb ********************************/
/*****************************************************************************/

// This operator does an inplace factor product with another
// TreeFactor.

class TreeFactor_InplaceAbsorb: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(TreeFactor* lhsPtr,
                      const Factor* rhsFPtr);
}; // TreeFactor_InplaceAbsorb

// This operator returns a factor product with another
// TreeFactor.

class TreeFactor_Absorb: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  Factor* process(const TreeFactor* lhsPtr,
                  const Factor* rhsFPtr);
}; // TreeFactor_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** TreeFactor_Cancel ********************************/
/*****************************************************************************/

// This operator does an inplace factor division with another
// TreeFactor.

class TreeFactor_InplaceCancel: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  void inplaceProcess(TreeFactor* lhsPtr,
                      const Factor* rhsFPtr);
}; // TreeFactor_InplaceCancel

// This operator returns a factor division with another
// TreeFactor.

class TreeFactor_Cancel: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  Factor* process(const TreeFactor* lhsPtr,
                  const Factor* rhsFPtr);
}; // TreeFactor_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** TreeFactor_Marginalize ***************************/
/*****************************************************************************/

/** This operator returns a marginalized factor, making use of the
 * covariance form representation. This version is particularly
 * efficient when combined with the Lauritzen-Spiegelhalter algorithm.
 */

class TreeFactor_Marginalize: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  Factor* process(const TreeFactor* lhsPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // TreeFactor_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** TreeFactor_ObserveAndReduce **********************/
/*****************************************************************************/

/**
 * After observing some variables, the scope of the resultant factor
 * reduces to ONLY that of the unobserved variables, making further
 * computations potentially much more efficient.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 */

class TreeFactor_ObserveAndReduce: public Operator1< TreeFactor > {
public:
  const std::string& isA() const;
  Factor* process(const TreeFactor* lhsPtr,
                const emdw::RVIds& variables,
                const emdw::RVVals& assignedVals,
                bool presorted);
}; // TreeFactor_ObserveAndReduce

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************
/**
 * Dampens/smooths two corresponding TreeFactors and returns the
 * total distance between them.
 *
 * new = lhs*(1-oldWeight) + rhs*oldWeight
 */
class TreeFactor_InplaceDamping: public Operator1< TreeFactor > {
public:

  const std::string& isA() const;

  /**
   * @param newPtr is the new message
   *
   * @param oldPtr is the old message
   *
   * @param oldWeight is the weighting/damping_factor for the old message
   */
  double
  inplaceProcess(
    TreeFactor* newPtr,
    const Factor* oldPtr,
    double oldWeight);

}; // TreeFactor_InplaceDamping


#endif // TREEFACTOR_HPP
