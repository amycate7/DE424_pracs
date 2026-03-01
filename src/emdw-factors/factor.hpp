#ifndef FACTOR_HPP
#define FACTOR_HPP

/*******************************************************************************
 *                                                                             *
 * AUTHOR:  JA du Preez                                                        *
 * Created on : /02/2013                                                       *
 *                                                                             *
 * PURPOSE: The abstract base class for all Factor variants                    *
 * COPYRIGHT: University of Stellenbosch, all rights reserved.                 *
 *                                                                             *
 *******************************************************************************/

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

#include "emdw.hpp"
#include "factoroperator.hpp"

/**
 * Factors/Potentials/Probability Distributions:

 * These are defined on one or more RV's and supply plausibility
 * values for values/assignments that those RV's can take on. We can
 * think of them as generalised probability distributions. They can
 * be a function of one or several Random Variables.  Examples are
 * Discrete RV probability tables on any enumerable type, Dirichlet
 * distributions, multi-dim Gaussians, Linear Gaussians, Conditional
 * Linear Gaussians, generalised Gaussians, Gaussian Mixtures (?)
 * etc. We also often use the log-form of factors.
 *
 * Also need to think about log-linear models and features - where
 * do they fit here?
 *
 * A Factor gets operated on by a FactorOperator. There are two
 * broad versions - the first changes a factor inplace while the
 * second builds a new factor while leaving the input factor(s)
 * unchanged. These factor operations will need varying additional
 * inputs to complete their allocated task.
 */


// Top level object
class Factor {

  //======================== Constructors and Destructor ========================
public:

  virtual ~Factor() {}
  //========================= Operators and Conversions =========================
public:
  bool operator==(const Factor& rhs) const;
  bool operator!=(const Factor& a) const;

  //=========================== Inplace Configuration ===========================
public:
  ///
  virtual std::istream& txtRead(std::istream& file) = 0;
  ///
  virtual std::ostream& txtWrite(std::ostream& file) const = 0;

  //===================== Required Factor Operations ============================

  /// The following functions provide the polymorphic access to
  /// derived classes, but will probably only be called indirectly via
  /// one of the similarly named non-virtual member functions lower
  /// down.
  ///
  /// NOTE: Those that return a uniqptr do so so that we can safely
  /// chain them. With a (covariant) raw pointer return chaining will
  /// cause memory leaks.

  //------------------Family 0: eg Sample

  /**
   * @param procPtr. Can be used to specify an optional sampling
   * operator for this operation; if NULL, this routine will use the
   * operator provided by the Factor object itself.
   *
   * @return a sample drawn from the factor
   *
   * NOTE: This should be a pure virtual, but because we added it much
   * later than most Factors, we do not immediately insist on
   * it. However, it is urged that all Factors be updated to also
   * support this.
   */
  virtual emdw::RVVals sample(Sampler* procPtr) const;

  //------------------Family 1: eg Normalize

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   *
   * NOTE: Quite often a factor is implicitly not normalizable. We will
   * therefore use the semantics that this routine should normalize
   * the factor if possible, but it is not a fatal error if the result
   * does not sum/integrate to unity. StillToBeDone: include a flag
   * that signals success.
   */
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr) const = 0;

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceNormalize(FactorOperator* procPtr) = 0;

  //------------------Family 2: eg Absorb,Cancel

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceAbsorb(const Factor* rhsPtr,
                             FactorOperator* procPtr) = 0;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   *
   * Result will inherit all the operators from the lhs factor
   */
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                 FactorOperator* procPtr) const = 0;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceCancel(const Factor* rhsPtr,
                             FactorOperator* procPtr) = 0;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   *
   * Result will inherit all the operators from the lhs factor
   */
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                 FactorOperator* procPtr) const = 0;

  //------------------Family 3: eg Marginalize

  /**
   * @param variablesToKeep are the ones to remain after the rest are
   * marginalised out.
   *
   * @param presorted Set this to true if the accompanying variable
   * ids are presorted from smallest to biggest.
   *
   * @param peekAheadPtr When doing message passing it is occasionally
   * possible for the marginalize operator to boost efficiency by
   * noting some extra info at the receiving cluster. Most factors,
   * however, can simply ignore this parameter.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                      bool presorted,
                                      const Factor* peekAheadPtr,
                                      FactorOperator* procPtr) const = 0;

  //------------------Family 4: eg ObserveAndReduce

  /**
   * @param variables - The id's of the variables that were
   * observed. These may validly also include variables not part of
   * this Factor.
   *
   * @param theirVals - the values they were observed to have
   *
   * @param presorted Set this to true if the accompanying variable
   * ids are presorted from smallest to biggest.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                           const emdw::RVVals& theirVals,
                                           bool presorted,
                                           FactorOperator* procPtr) const = 0;

  //------------------Family 5: eg inplaceDampen

  /**
   * This function smooths the messages that are being passed in
   * belief-propagation algorithms.
   *
   * @param df is the damping factor. msg = df*oldMsg + (1-df)*msg (ie
   *   a first order IIR filter smoothing). The heavier the damping
   *   factor, the closer we stick to the old message.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   *
   * It is assumed that both messages are of identical type with fully
   * overlapping variables.
   *
   * returns a measure of the distance between the (pre-updated)
   * existing message and the new message.
   *
   * NOTE: Damping does not always result in the same Factor type as
   * the two messages being damped. E.g. two Gaussians turns into a
   * Mixture Gaussian. Quite often some form of approximation will be
   * necessary.
   */
  virtual double inplaceDampen(const Factor* oldMsg,
                               double df,
                               FactorOperator* procPtr) = 0;

  //===================== Required Virtual Member Functions =====================
public:

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace, one-for-one AND IN
   * THIS ORDER, the original SORTED set. This, of course, also
   * implies that your new variables must be fully compatible with the
   * ones they replace i.e. same domain etc.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual Factor* copy(const emdw::RVIds& newVars = {},
                       bool presorted = false) const = 0;

  /**
   * Make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   *
   * @param selectedVars Can specify any subset of the variables of
   * this factor to form the vacuous factor for. If empty all
   * variables will be used.
   *
   * @param presorted Set this to true if selectedVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.

   */
  virtual Factor* vacuousCopy(const emdw::RVIds& selectedVars = {},
                              bool presorted = false) const = 0;

  /// note rhs is a raw pointer
  virtual bool isEqual(const Factor* rhsPtr) const = 0;

  /**
   * Returns the distance of this potential to its uniform/flat/
   * vacuous version. (Typically we use the one-sided) KL(P||Q)
   * distance with this potential being P and the vacuous one being
   * Q.)
   *
   * NOTE: Older code can easily implement this via the 'vacuousCopy'
   * and deprecated 'distance' member functions. See 'factor.cc' for
   * an example.
   */
  virtual double distanceFromVacuous() const = 0;

  ///
  virtual unsigned noOfVars() const = 0;

  /// returns the ids of the variables the factor is defined on, in
  /// *sorted* order
  virtual emdw::RVIds getVars() const = 0;

  ///in *sorted* order
  virtual emdw::RVIdType getVar(unsigned varNo) const = 0;

  //====================== Other Virtual Member Functions =======================
public:

  /**
   * DEPRECATED. Measures the distance between two factors. zero if
   * equal, else positive and growing with degree of difference.
   */
  virtual double distance(const Factor* rhsPtr) const;

  /**
   * returns the factor potential at a given variable
   * assignment. *NOTE* the potential is supplied directly as is, no
   * normalization etc is done.
   *
   * @param theVars The RV ids in the order that the values in
   * theirVals assumes them to be.
   *
   * @param theirVals The values assigned to this random vector.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual double potentialAt(const emdw::RVIds& theVars,
                             const emdw::RVVals& theirVals,
                             bool presorted = false) const;

  /**
   * measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference.
   */
  virtual double distance(const rcptr<Factor>& rhsPtr) const;

  //======================= Non-virtual Member Functions ========================
public:

/**
 * returns lhs operator applied to all of the rhs factors
 *
 * @param procPtr. Can be used to specify an optional operator for
 * this operation; if NULL, this routine will use the operator
 * provided by the lhs Factor object itself.
 */
  uniqptr<Factor> absorb(const std::vector< rcptr<Factor> >& facPtrs,
                         const rcptr<FactorOperator>& procPtr = 0) const;

  /// In actual use factor operators will probably be activated via
  /// one or more of the following functions

  //------------------Family 0: eg Sample

  /**
   * @param procPtr. Can be used to specify an optional sample
   * operator for this operation; if NULL, this routine will use the
   * operator provided by the Factor object itself.
   *
   * @return a sample drawn from the factor
   */
  emdw::RVVals sample(const rcptr<Sampler>& procPtr = 0) const;

  //------------------Family 1: eg Normalize

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> normalize(const rcptr<FactorOperator>& procPtr = 0) const;

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  void inplaceNormalize(const rcptr<FactorOperator>& procPtr = 0);

  //------------------Family 2: eg absorb

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> absorb(const rcptr<Factor>& rhsPtr,
                         const rcptr<FactorOperator>& procPtr = 0) const;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  void inplaceAbsorb(const rcptr<Factor>& rhsPtr,
                     const rcptr<FactorOperator>& procPtr = 0);

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> cancel(const rcptr<Factor>& rhsPtr,
                         const rcptr<FactorOperator>& procPtr = 0) const;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  void inplaceCancel(const rcptr<Factor>& rhsPtr,
                     const rcptr<FactorOperator>& procPtr = 0);

  //------------------Family 3: eg Marginalize

  /**
   * @param variablesToKeep are the ones to remain after the rest are
   * marginalised out.
   *
   * @param presorted Set this to true if the accompanying variable
   * ids are presorted from smallest to biggest.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                              bool presorted = false,
                              const Factor* peekAheadPtr = nullptr,
                              const rcptr<FactorOperator>& procPtr = 0) const;

  /**
   * @param variablesToZap are the ones to marginalised out.
   *
   * @param presorted Set this to true if the accompanying variable
   * ids are presorted from smallest to biggest.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> marginalizeX(const emdw::RVIds& variablesToZap,
                               bool presorted = false,
                               const Factor* peekAheadPtr = nullptr,
                               const rcptr<FactorOperator>& procPtr = 0) const;

  //------------------Family 4: eg ObserveAndReduce

  /**
   * @param variables - the id's of the variables that were observed
   *
   * @param theirVals - the values they were observed to have
   *
   * @param presorted Set this to true if the accompanying variable
   * ids are presorted from smallest to biggest.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                   const emdw::RVVals& theirVals,
                                   bool presorted = false,
                                   const rcptr<FactorOperator>& procPtr = 0) const;

  //------------------Family 5: eg inplaceDampen

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  double inplaceDampen(const rcptr<Factor>& oldMsg,
                       double df,
                       const rcptr<FactorOperator>& procPtr = 0);


protected:

  // we use raw pointers for efficiency and covariant returns

  //------------------Family 0

  template<typename FactorType>
  inline emdw::RVVals dynamicSample(
    FactorOperator* procPtr,
    const FactorType* self) const;

  //------------------Family 1

  template<typename FactorType>
  inline Factor* dynamicApply(FactorOperator* procPtr,
                              const FactorType* self) const;

  template<typename FactorType>
  inline void dynamicInplaceApply(FactorOperator* procPtr,
                                  FactorType* self);

  //------------------Family 2

  template<typename FactorType>
  inline Factor* dynamicApply(FactorOperator* procPtr,
                              const FactorType* self,
                              const Factor* other) const;

  template<typename FactorType>
  inline void dynamicInplaceApply(FactorOperator* procPtr,
                                  FactorType* self,
                                  const Factor* other);

  //------------------Family 3

  template<typename FactorType>
  inline Factor* dynamicApply(FactorOperator* procPtr,
                              const FactorType* self,
                              const emdw::RVIds& variablesToKeep,
                              bool presorted = false,
                              const Factor* peekAheadPtr = nullptr) const;

  //------------------Family 4

  template<typename FactorType>
  inline Factor* dynamicApply(FactorOperator* procPtr,
                              const FactorType* self,
                              const emdw::RVIds& variables,
                              const emdw::RVVals& theirVals,
                              bool presorted = false) const;

  //dynamicInplaceApply(procPtr, this, variables, theirVals, presorted);
  template<typename FactorType>
  inline void dynamicInplaceApply(FactorOperator* procPtr,
                                  FactorType* self,
                                  const emdw::RVIds& variables,
                                  const emdw::RVVals& theirVals,
                                  bool presorted = false);

  //------------------Family 5

  template<typename FactorType>
  inline double dynamicInplaceApply(FactorOperator* procPtr,
                                    FactorType* self,
                                    const Factor* other,
                                    double df);

  //=============================== Data Members ================================
public:

  // A flag to indicate whether this object may be merged with others
  // during the PGM creation. Bigger values are more likely to form a
  // hub for clusters. Typically we dont consider this part of the
  // state of the factor (i.e. you wont find it in the constructors,
  // we dont save it to file etc), it is just something to manipulate
  // the cluster-graph configuration with.
  unsigned keepSeparate = 0;

  //================================== Friends ==================================
public:

  /// output
  friend std::ostream& operator<<( std::ostream& file,
                                   const Factor& n );

  /// input
  friend std::istream& operator>>( std::istream& file,
                                   Factor& n );

  friend uniqptr<Factor> absorb(
    const std::vector< rcptr<Factor> >& facPtrs,
    const rcptr<FactorOperator>& procPtr);

}; // Factor

/**
 * Apply the operator sequentially on a sequence of factors
 *
 * @param procPtr. Can be used to specify an optional operator for
 * this operation; if NULL, this routine will use the operator
 * provided by the Factor object itself.
 */
uniqptr<Factor> absorb(const std::vector< rcptr<Factor> >& facPtrs,
                       const rcptr<FactorOperator>& procPtr = 0);



#include "factor.tcc"


/**
 * FACTOR CLOSURES:
 *
 * Discrete => Discrete (always, might be in sub-families)
 *
 * Gaussian:
 *   Prod: Gaussian
 *   Marg: Gaussian
 *   Obsv: Gaussian
 *   Smth: MG   <============
 *
 * MG (Mixture Gaussian):
 *   Prod: MG
 *   Marg: MG
 *   Obsv: MG
 *
 * CG (Conditional Gaussian):
 *   Prod: CG
 *   Marg: Over Disc => CMG or MG
 *         Over Cont => CG or Disc
 *   Obsv: Disc => CG or Gaussian
 *         Cont => CG or Discrete
 *
 * CMG (Conditional Mixture Gaussian):
 *   Prod: CMG
 *   Marg: Over Disc => CMG or MG
 *         Over Cont => CMG or Discrete
 *   Obsv: Disc => CMG or MG
 *         Cont => CMG or Discrete
 *
 * LG (Linear Gaussian):
 * CLG (Conditional Linear Gaussian):
 *
 */

#endif // FACTOR_HPP
