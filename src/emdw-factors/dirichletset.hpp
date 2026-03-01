#ifndef DIRICHLETSET_HPP
#define DIRICHLETSET_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    April 2018
          PURPOSE: DirichletSet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "dirichlet.hpp"

// standard headers
#include <iosfwd>  // istream, ostream
#include <map>

// Some forward declarations
template <typename CondType, typename DirType> class DirichletSet_Marginalize;
template <typename CondType, typename DirType> class DirichletSet_InplaceNormalize;
template <typename CondType, typename DirType> class DirichletSet_InplaceAbsorb;
template <typename CondType, typename DirType> class DirichletSet_InplaceCancel;
template <typename CondType, typename DirType> class DirichletSet_InplaceDamping;
template <typename CondType, typename DirType> class ConditionalPolya;

/******************************************************************************/
/***************************** DirichletSet ***************************/
/******************************************************************************/

/**
 * This class provides a prior for the probabilities of conditional
 * Categorical (i.e. DiscreteTable) variables. Depending on the value
 * of the conditioned-on variable (of type CondType), we select a
 * particular Dirichlet prior (of type DirType). The internal
 * parameters of these Dirichlets may vary, but all other aspects
 * (such as the range of values they describe etc) are identical. Note
 * that whereas DiscreteTable requires all RVs to be of the same
 * underlying type, DirichletSet allows the conditioned-on RV
 * to have a different type (CondType) to that of the target RV
 * (DirType).
 *
 * Some examples:
 * 1. The topic models in Latent Dirichlet Allocation are each
 *    described in terms of the words occurring in it, i.e. a
 *    Dirichlet model. The DirichletSet then models the
 *    collection of these topics.
 * 2. In language modelling we are concerned with the probability of a
 *    word occurring given some preceding word (e.g. an n-gram
 *    model), i.e. once again a Dirichlet model. The collection of
 *    these n-gram models are a DirichletSet where these prior
 *    conditioned-on words indicate the particular n-gram.
 *
 * Since there essentially are no interactions between these
 * Dirichlets, their FactorOperators are straightforward.
 *
 * NOTE: We have a design choice on whether we will associate the set
 * with a single RV id (thereby hiding the ids of the underlying
 * Dirichlets), or directly use the joint ids of all the underlying
 * Dirichlets. To keep things simple (and also follow the approach
 * taken in Dirichlet) we opt for the former. This makes it possible
 * to keep the Dirichlets in a vector instead of a map, and also
 * simplifies the manipulation of variable ids considerably. On the
 * downside this unfortunately also precludes more complex scenarios
 * with multiple DirichletSets sharing subsets of their
 * variables. If the need for this arises in the future we might want
 * to reconsider/extend.
 *
 * NOTE: ConditionalPolya is the closely associated class joining the
 * the DirichletSet and Categorical variables.
 *
 * TODO: Major efficiency problem:
 * In applications where the target (y) variable is observed,
 * this class (as it currently stands) will be exceedingly
 * expensive. During the marginalisation step the whole
 * DirichletSet (all the underlying Dirichlets in all their
 * glory) will get copied (expensive) In the immediately following
 * division/cancel step the previous incoming message will get divided
 * into this - once again building a next full DirichletSet
 * (expensive). Then only after being absorbed at the ConditionalPolya
 * do we get the chance to let the observed y clean out most of this
 * DirichletSet.
 *
 * We can improve at least the second one of these inefficiencies
 * (i.e. the divide/cancel operation) by letting the
 * DirichletSet also include a field for an observed y
 * value. If two DirichletSets get multiplied, and only one
 * has an observed y, the result will still be a full Conditional
 * Dirichlet. However, if both has the same observed y, all the other
 * target parameters will resort to default. This makes of it a very
 * sparse DirichletSet with only one non-default value -
 * indexed by that observed y. If we divide a full
 * DirichletSet by a DirichletSet with observed Y, the
 * reasoning is that the message obviously is en-route to a
 * ConditionalPolya with that same observed y. One can therefore
 * short-circuit things and immediately revert to the sparse version
 * containing only that observed y.
 *
 * I am still in the dark as how to solve the first inefficiency
 * (marginalisation). Vague ideas are to somehow alias the
 * DirichletSet and only make a separate copy of it once it
 * becomes certain that it is going to be changed.
 */
template <typename CondType, typename DirType>
class DirichletSet : public Factor {

  //TODO Reduce these dependencies later?
  friend class DirichletSet_Marginalize<CondType, DirType>;
  friend class DirichletSet_InplaceNormalize<CondType, DirType>;
  friend class DirichletSet_InplaceAbsorb<CondType, DirType>;
  friend class DirichletSet_InplaceCancel<CondType, DirType>;
  friend class DirichletSet_InplaceDamping<CondType, DirType>;
  friend class ConditionalPolya<CondType, DirType>;

  //============================ Traits and Typedefs ============================
private:

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:
  /**
   * @brief Default constructor
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  DirichletSet(
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0);

  /**
   * The class specific constructor.
   *
   * @param theVarId: the (single) id to assign to this set of
   * distributions. The actual ids of the underlying Dirichlets are
   * not visible outside of this class.
   *
   * @param theDomain: the values that the Categorical RV in the
   * conditional can take on.
   *
   * @param theDirPtrs: The actual set of Dirichlet priors.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  DirichletSet(
    emdw::RVIdType theVarId,
    const rcptr< std::vector<CondType> >& theDomain,
    const std::map< CondType, rcptr< Dirichlet<DirType> > >& theDirPtrs,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0);

  /// copy and move constructors
  DirichletSet(const DirichletSet<CondType, DirType>& st);
  DirichletSet(DirichletSet<CondType, DirType>&& st);

  /// DirichletSet class destructor
  virtual ~DirichletSet();

  //========================= Operators and Conversions =========================
public:
  /// Assignment operator
  // ???JADP??? !!!As is here the individual Dirichlets will be copied shallowly
  DirichletSet& operator=(const DirichletSet<CondType, DirType>& d) = default;
  DirichletSet& operator=(DirichletSet<CondType, DirType>&& d) = default;

  /// == operator
  bool operator==(const DirichletSet<CondType, DirType>& d) const;

  /// != operator
  bool operator!=(const DirichletSet<CondType, DirType>& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /**
   * The class specific configurer.
   *
   * @param theVarId: the (single) id to assign to this set of
   * distributions. The actual ids of the underlying Dirichlets are
   * not visible outside of this class.
   *
   * @param theDomain: the values that each Categorical RV in the
   * conditional can take on.
   *
   * @param theDirPtrs: The actual set of Dirichlet priors.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  unsigned classSpecificConfigure(
    emdw::RVIdType theVarId,
    const rcptr< std::vector<CondType> >& theDomain,
    const std::map< CondType, rcptr< Dirichlet<DirType> > >& theDirPtrs,
    rcptr<Marginalizer> theMarginalizer = 0,
    rcptr<InplaceNormalizer> theInplaceNormalizer = 0,
    rcptr<Normalizer> theNormalizer = 0,
    rcptr<InplaceAbsorber> theInplaceAbsorber = 0,
    rcptr<Absorber> theAbsorber = 0,
    rcptr<InplaceCanceller> theInplaceCanceller = 0,
    rcptr<Canceller> theCanceller = 0,
    rcptr<ObserveAndReducer> theObserveAndReducer = 0,
    rcptr<InplaceDamper> theInplaceDamper = 0);

  // These functions will be used in the future to enable users to
  // write DirichletSet factors to file and to read DirichletSet
  // factors from file.
  virtual std::istream& txtRead(std::istream& file);
  virtual std::ostream& txtWrite(std::ostream& file) const;


  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from (the top-level) Factor. They return uniqptr so
  /// that we can chain them.

  //------------------Family 1: eg Normalize

  virtual void inplaceNormalize(FactorOperator* procPtr);
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr) const;

  //------------------Family 2: eg Absorb,Cancel

  virtual void inplaceAbsorb(const Factor* rhsPtr,
                            FactorOperator* procPtr);
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr) const;

  virtual void inplaceCancel(const Factor* rhsPtr,
                            FactorOperator* procPtr);
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr) const;

  //------------------Family 3 eg Marginalize

  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                      bool presorted,
                                      const Factor* peekAheadPtr,
                                      FactorOperator* procPtr) const;


  //------------------Family 4: eg Observe, ObserveAndReduce

  virtual uniqptr<Factor> observeAndReduce(const emdw::RVIds& variables,
                                           const emdw::RVVals& theirVals,
                                           bool presorted,
                                           FactorOperator* procPtr) const;

  //------------------Family 5: eg inplaceDampen

  virtual double inplaceDampen(const Factor* oldMsg,
                               double df,
                               FactorOperator* procPtr);

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
   * @param presorted Not used.
   */
  virtual DirichletSet<CondType, DirType>* copy(const emdw::RVIds& newVars = {}, bool presorted = true) const;

  /**
   * Make copy of this factor, but set the potentials to the default
   * non informative factor values. The client has to manage the
   * lifetime of this pointer - putting it in a rcptr will be a good
   * idea.
   */
  virtual DirichletSet<CondType, DirType>*
  vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /**
   * Checks if two DirichletSets are identical.
   */
  virtual bool isEqual(const Factor* rhsPtr) const ;

  /// Total KL distance to corresponding vacuous Dirichlets
  double distanceFromVacuous() const;

  /**
   * Returns 1 - the number of RV ids in a DirichletSet.
   */
  virtual unsigned noOfVars() const;

  /**
   * Returns the vector of RV ids (only 1 element).
   */
  virtual emdw::RVIds getVars() const;

  /**
   * returns the RV id of the dirichlet var. Since there is only one var,
   * the only valid assignment of varNo=0
   */
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //====================== Other Virtual Member Functions =======================
public:

  /**
   * Measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference. Currently based
   * on the Kulback Leibler distance.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Non-virtual Member Functions ========================
public:

  /**
   * Slices over all the Dirichlets in the set at a particular
   * observed value.
   *
   * @param theObservedWord The observed value at which to slice over
   * the Dirichlets in the set.
   *
   * @return A DirichletSet with all alphas zero'ed except for
   * the alpha at theObservedWord */
  DirichletSet<CondType, DirType>*
  sliceAt(const DirType& theObservedWord) const;

  /// replaces various extra settings (those following after
  /// sparseProbs in constructor) with those of d.
  void copyExtraSettingsFrom(const DirichletSet<CondType, DirType>& d);

  /// @param dirIdx The index of the particular Dirichlet to retrieve
  uniqptr< Dirichlet<DirType> > copyDirichlet(const CondType& dirIdx) const;

  /// @param dirIdx The index of the particular Dirichlet to retrieve
  rcptr< Dirichlet<DirType> > aliasDirichlet(const CondType& dirIdx);

  ///
  emdw::RVVals sampleFrom(const CondType dirIdx) const;

  ///
  size_t noOfDirichlets() const;

  /**
   * returns the number of categories per Dirichlet
   */
  size_t noOfCategories() const;

  /**
   * Adds cnt new observations to category catNr of Dirichlet setNr
   */
  void addToCounts(const CondType& setNr,
                   const DirType& catNr,
                   double cnt);


  //================================== Friends ==================================

public:

  /// output
  template<typename T>
  friend std::ostream& operator<<( std::ostream& file,
                                   const DirichletSet<CondType, DirType>& n );

  /// input
  template<typename T>
  friend std::istream& operator>>( std::istream& file,
                                   DirichletSet<CondType, DirType>& n );

  //=============================== Data Members ================================
private:

  emdw::RVIdType varId;

  // --- the values that each Categorical RV in the conditional can take on
  rcptr< std::vector<CondType> > domain;
  std::map< CondType, rcptr< Dirichlet<DirType> > > dirPtrs;

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

}; // DirichletSet


/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalization ****************************/
/*****************************************************************************/

/**
 * This operator does an inplace normalization of a DirichletSet.
 */
template <typename CondType, typename DirType>
class DirichletSet_InplaceNormalize: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(DirichletSet<CondType, DirType>* lstPtr);
}; // DirichletSet_InplaceNormalize

/**
 * This operator returns a normalized version of a DirichletSet.
 */
template <typename CondType, typename DirType>
class DirichletSet_Normalize: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const DirichletSet<CondType, DirType>* lstPtr);
}; // DirichletSet_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

/**
 * This operator inplace multiplies a DirichletSet by another DirichletSet.
 */
template <typename CondType, typename DirType>
class DirichletSet_InplaceAbsorb: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(DirichletSet<CondType, DirType>* lstPtr, const Factor*);
}; // DirichletSet_InplaceAbsorb

/**
 * This operator returns the factor product of a DirichletSet by another
 * DirichletSet.
 */
template <typename CondType, typename DirType>
class DirichletSet_Absorb: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const DirichletSet<CondType, DirType>* lstPtr, const Factor*);
}; // DirichletSet_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

/**
 * This operator does an inplace divide with the corresponding
 * Dirichlet in a second set.
 */
template <typename CondType, typename DirType>
class DirichletSet_InplaceCancel: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(DirichletSet<CondType, DirType>* lstPtr, const Factor*);
}; // DirichletSet_InplaceCancel

/**
 * This operator returns the quotient of each Dirichlet with the
 * corresponding one in a second set.
 */
template <typename CondType, typename DirType>
class DirichletSet_Cancel: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const DirichletSet<CondType, DirType>* lstPtr, const Factor*);
}; // DirichletSet_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/
/**
 * We only support keeping the full DirichletSet i.e. not reducing the
 * scope at all.  This (trivial) marginalization may be required in
 * loopy inference algorithms.
 */
template <typename CondType, typename DirType>
class DirichletSet_Marginalize: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
//  Factor* process(const DirichletSet<CondType, DirType>* lstPtr);
  Factor* process(const DirichletSet<CondType, DirType>* lstPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // DirichletSet_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/
/**
 * Currently this operator is NIY.
 *
 * TODO: Consider allowing observing the indexing (y) variable to be
 * observed, thereby doing something similar to the sliceAt operation,
 * resulting in a DirichletSet with only one non-default alpha
 * value (at the observed y index).
 */
template <typename CondType, typename DirType>
class DirichletSet_ObserveAndReduce: public Operator1< DirichletSet<CondType, DirType> > {
public:

  const std::string& isA() const {
    static const std::string CLASSNAME(
      "DirichletSet_ObserveAndReduce<"+std::string(typeid(DirType).name())+std::string(">"));
    return CLASSNAME;
  } // isA

  Factor* process(const DirichletSet<CondType, DirType>* lstPtr,
              const emdw::RVIds& variables,
              const emdw::RVVals& assignedVals,
              bool presorted) {
    PRLITE_ASSERT(
      false,
      "DirichletSet_ObserveAndReduce::process NIY since it was not expected to ever be used!");
    return (DirichletSet<CondType, DirType>*) nullptr;
  } // process

}; // DirichletSet_ObserveAndReduce

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/
/**
 * Dampens two DirichletSets. The dampingfactor 'df' is the weight of
 * the old factor, i.e. df=0 implies no damping.
 *
 * new = lhs*(1-df) + rhs*df
 */
template <typename CondType, typename DirType>
class DirichletSet_InplaceDamping: public Operator1< DirichletSet<CondType, DirType> > {
public:
  const std::string& isA() const;
  double inplaceProcess(DirichletSet<CondType, DirType>* lstPtr, const Factor*, double df);
}; // DirichletSet_InplaceDamping

#include "dirichletset.tcc"

#endif // DIRICHLETSET_HPP
