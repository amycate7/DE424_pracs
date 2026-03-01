#ifndef DIRICHLET_HPP
#define DIRICHLET_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2018
          PURPOSE: Dirichlet PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"

// standard headers
#include <iosfwd>  // istream, ostream
#include <map>

// forward declarations
template <typename AssignType> class Dirichlet_Marginalize;
template <typename AssignType> class Dirichlet_InplaceNormalize;
template <typename AssignType> class Dirichlet_InplaceAbsorb;
template <typename AssignType> class Dirichlet_InplaceCancel;
template <typename AssignType> class Dirichlet_InplaceWeakDamping;
template <typename AssignType> class Dirichlet_Sampler;
template <typename AssignType> class Cat1VMP;


/******************************************************************************/
/***************************** Dirichlet **************************************/
/******************************************************************************/
/**
 * This class provides a prior for the probabilities of a (single)
 * Categorical RV.
 * .
 *
 * A Dirichlet distributed random VECTOR describing the PROBABILITIES
 * of a K-valued SCALAR Categorial RV consists of K probability values
 * theta_0 upto theta_{K-1} that jointly sum to one. Since we are not
 * going to concern ourselves with subsets of these K probability
 * values, we will always refer to it as a whole and therefore
 * identify it with only a single RVId.
 *
 * From the above it might (erroneously) appear that the various
 * thetas are indexed with unsigned values 0..K-1. Not so - we keep
 * track of the domain of the Categorical RV and use its valid
 * assignments to index the various corresponding probabilities.
 *
 * Note: the current implementation does NOT keep track of the scaling
 * constants that keeps its volume to unity. This implies that the
 * class is implicitly normalized. This is problematic in applications
 * where one explicitly want to keep track of the normalization
 * constants.
 *
 * NOTE: Polya is the closely associated class joining the the
 * Dirichlet and Categorical variables.
 *
 * TODO: Maybe we should also provide for a relFloor similar to what
 * is done in DiscreteTable???
 */
template <typename AssignType>
class Dirichlet : public Factor {

  //TODO Reduce these dependencies later?
  friend class Dirichlet_Marginalize<AssignType>;
  friend class Dirichlet_InplaceNormalize<AssignType>;
  friend class Dirichlet_InplaceAbsorb<AssignType>;
  friend class Dirichlet_InplaceCancel<AssignType>;
  friend class Dirichlet_InplaceWeakDamping<AssignType>;
  friend class Dirichlet_Sampler<AssignType>;
  friend class Cat1VMP<AssignType>;

  //============================ Traits and Typedefs ============================
private:

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:
  /**
   * @brief Default constructor.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  Dirichlet(
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
   * @brief The class specific constructor for creating a Dirichlet
   * factor with given parameter values.
   *
   * @param theVarId This single RVID represents a vector of theta
   * random variables. Implemented in this way since we see no use in
   * referring to subsets of theta random variables and this notation
   * significantly reduces complexity.
   *
   * @param theDomain: The values the Categorical variable can take
   * on.
   *
   * @param theDefAlpha: The default alpha count (>=0.0) to use for
   * all cases not specified in theSparseAlphas. For the
   * non-informative case theDefAlpha=1 and sparseAlphas will be
   * empty. (You can think of (alpha_k - 1) as the number of prior
   * observations of Categorical variable being its k'th value.)
   *
   * @param theSparseAlphas A map containing those specific prior
   * alpha counts (>=0.0) that differ from the default value. For the
   * non-informative case this map will be empty.
   *
   * @param theRelMargin Any alpha counts smaller than
   * (theRelMargin*max(alphas) ) will be set to 0.0.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  Dirichlet(
    emdw::RVIdType theVarId,
    rcptr< std::vector<AssignType> > theDomain,
    double theDefAlpha = 1.0,
    const std::map<AssignType, double>& theSparseAlphas = std::map<AssignType, double>(),
    double theRelMargin = 0.0,
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
   * @brief A conveniency constructor to create uniform alphas with
   * some perturbation. Note, to disambiguate from other constructors,
   * the domain is the FIRST parameter here, not the second.
   *
   * @param theDomain: The values the Categorical variable can take
   * on.
   *
   * @param theVarId This single RVID represents a vector of theta
   * random variables. Implemented in this way since we see no use in
   * referring to subsets of theta random variables and this notation
   * significantly reduces complexity.
   *
   * @param nominalAlpha: A base level to use for all alphas;
   *
   * @param perturbAlpha Max value of positive uniform perturbation to
   * alphas;
   *
   * @param theRelMargin Any alpha counts smaller than
   * (theRelMargin*max(alphas) ) will be set to 0.0.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  Dirichlet(
    rcptr< std::vector<AssignType> > theDomain,
    emdw::RVIdType theVarId,
    double nominalAlpha,
    double perturbAlpha,
    double theRelMargin = 0.0,
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

  /// copy and move constructors
  Dirichlet(const Dirichlet<AssignType>& p)  = default;
  Dirichlet(Dirichlet<AssignType>&& p) = default;

  /// Dirichlet class destructor
  virtual ~Dirichlet();

  //========================= Operators and Conversions =========================
public:
  /// Assignment operator
  Dirichlet& operator=(const Dirichlet<AssignType>& d) = default;
  Dirichlet& operator=(Dirichlet<AssignType>&& d) = default;

  /// == operator
  bool operator==(const Dirichlet<AssignType>& d) const;

  /// != operator
  bool operator!=(const Dirichlet<AssignType>& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /**
   * @brief Dirichlet class specific configurer. This function
   * can be used to populate an already defined Dirichlet
   * factor.
   *
   * @param theVarId This single RVID represents a vector of theta
   * random variables. Implemented in this way since we see no use in
   * referring to subsets of theta random variables and this notation
   * significantly reduces complexity.
   *
   * @param theDomain: The values that the Categorical variable can
   * take on.
   *
   * @param theDefAlpha: The default alpha count (>=0.0) to use for
   * all cases not specified in theSparseAlphas. For the
   * non-informative case theDefAlpha=1 and sparseAlphas will be
   * empty. (You can think of (alpha_k - 1) as the number of prior
   * observations of Categorical variable being its k'th value.)
   *
   * @param theSparseAlphas A map containing those specific prior
   * alpha counts (>=0.0) that differ from the default value. For the
   * non-informative case this map will be empty.
   *
   * @param theRelMargin Any alpha counts smaller than
   * (theRelMargin*max(alphas) ) will be set to 0.0.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  unsigned classSpecificConfigure(
    emdw::RVIdType theVarId,
    rcptr< std::vector<AssignType> > theDomain,
    double theDefAlpha = 1.0,
    const std::map<AssignType, double>& theSparseAlphas = std::map<AssignType, double>(),
    double theRelMargin = 0.0,
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

  // These functions will be used in the future to enable users to
  // write Dirichlet factors to file and to read Dirichlet
  // factors from file.
  virtual std::istream& txtRead(std::istream& file);
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
   * with params it returns an identical Factor, but with the variable
   * replaced by another.
   *
   * @param newVars This variable will replace the original one.
   *
   * @param presorted This variable is necessary from the declaration
   * of the factor class, but since there is only one variable in the
   * Dirichlet class, the value of presorted will always be true.
   */
  virtual Dirichlet<AssignType>* copy(const emdw::RVIds& newVars = {}, bool presorted = true) const;

  /**
   * Make a copy of this factor, but set the potentials to the default
   * non informative factor values. The client has to manage the
   * lifetime of this pointer - putting it in a rcptr will be a good
   * idea.
   */
  virtual Dirichlet<AssignType>*
  vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /**
   * Checks if two Dirichlets are identical.
   */
  virtual bool isEqual(const Factor* rhsPtr) const ;

  /// Symmetrised KL distance to a uniform Dirichlet
  virtual double distanceFromVacuous() const;

  /**
   * Returns 1 - the number of RV ids in a Dirichlet.
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
   * on the Kulbach Leibler distance.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Non-virtual Member Functions ========================
public:

  /// Returns a very sparse Dirichlet the alpha at theObservedWord is
  /// retained and all others are replaced with a default value such
  /// that the sumOfAlphas remains unchanged.
  Dirichlet<AssignType>* sliceAt(const AssignType& theObservedWord) const;

  /** If the default alpha is not zero, fill in all the missing
   * probablities in sparseAlphas (equal to the default value) and
   * then change the default alpha to zero. Careful with this,
   * the resulting table might be huge.
   */
  void makeDefaultZero() const;

  /**
   * returns the number of categories
   */
  size_t noOfCategories() const;

  /**
   * @brief Adds cnt new observations to category k.
   *
   * @param k Is an assignment to the Categorical RV underlying this
   * Dirichlet.
   */
  void addToCounts(const AssignType& k,
                   double cnt);


  // returns the default alpha count
  double getDefaultAlpha() const;

  // returns the map with sparse alpha values
  const std::map<AssignType, double>& getSparseAlphas() const;

  /**
   * @return the alpha[k] count.
   *
   * @param k is the <AssignType> Dirichlet category, i.e. the value
   * of the underlying discrete variable.
   */
  double getAlpha(const AssignType& k) const;

  /// replaces the factor operators with those of d
  void copyExtraSettingsFrom(const Dirichlet<AssignType>& d);

  // rescale the alphas to sum to a specified value.
  void rescaleSumOfAlphasTo(double newSum);

  // rescale the alphas to not exceed this value. This is useful to
  // allow the Dirichlet to adapt in non-stationary circumstances.
  void limitSumOfAlphasTo(double noMoreThan);

  /// returns the sum of the alphas
  double sumAlphas() const;

private:
  void recalcSumOfAlphas() const;

  //================================== Friends ==================================

public:

  /// output
  template<typename T>
  friend std::ostream& operator<<( std::ostream& file,
                                   const Dirichlet<AssignType>& n );

  /// input
  template<typename T>
  friend std::istream& operator>>( std::istream& file,
                                   Dirichlet<AssignType>& n );

  //=============================== Data Members ================================
private:

  //+++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++ definitive data members ++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++

  // --- the id of the RV
  emdw::RVIdType varId;

  // --- the values that each Categorical RV can take on
  rcptr< std::vector<AssignType> > domain;

  // --- the alpha counts >= 0. All ones means no info i.e. uniform
  mutable double defAlpha;  // default value for everything not listed in sparseAlphas
  std::map<AssignType, double> sparseAlphas; // alpha counts: empty for uniform

  // --- alphas less than this times max alpha will be set to 0.0
  double relMargin = 0.0;

  // --- JADP seems we might need some normalization stuff here? For
  // --- instance, what happens when we multiply two Dirichlets?

  // --- The set of FactorOperators in use
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

  //+++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++ derived data members +++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++

  mutable double sumOfAlphas; // useful to have around to make things faster

}; // Dirichlet


// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a Dirichlet.
 */
template <typename AssignType>
class Dirichlet_Sampler: public SampleOperator< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  emdw::RVVals sample(const Dirichlet<AssignType>* lstPtr);
}; // Dirichlet_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalization ****************************/
/*****************************************************************************/

/**
 * This operator does an inplace normalization of the Dirichlets.
 */
template <typename AssignType>
class Dirichlet_InplaceNormalize: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Dirichlet<AssignType>* lstPtr);
}; // Dirichlet_InplaceNormalize

/**
 * This operator returns a normalized version of the Dirichlets.
 */
template <typename AssignType>
class Dirichlet_Normalize: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Dirichlet<AssignType>* lstPtr);
}; // Dirichlet_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

/**
 * This operator does an inplace multiply of two Dirichlet factors.
 */
template <typename AssignType>
class Dirichlet_InplaceAbsorb: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Dirichlet<AssignType>* lstPtr, const Factor*);
}; // Dirichlet_InplaceAbsorb

/**
 * This operator returns the product of two Dirichlet factors.
 */
template <typename AssignType>
class Dirichlet_Absorb: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Dirichlet<AssignType>* lstPtr, const Factor*);
}; // Dirichlet_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

/**
 * This operator inplace divides by the second Dirichlet.
 */
template <typename AssignType>
class Dirichlet_InplaceCancel: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Dirichlet<AssignType>* lstPtr, const Factor*);
}; // Dirichlet_InplaceCancel

/**
 * This operator returns the result from dividing by the second Dirichlet.
 */
template <typename AssignType>
class Dirichlet_Cancel: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Dirichlet<AssignType>* lstPtr, const Factor*);
}; // Dirichlet_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/
/**
 * We currently only support keeping the full Dirichlet i.e. not
 * reducing the scope at all.  This (trivial) marginalization may be
 * required in loopy inference algorithms.
 *
 * Note: If need be it should be easy to extend this to retain a
 * subset of the Dirichlets.
 */
template <typename AssignType>
class Dirichlet_Marginalize: public Operator1< Dirichlet<AssignType> > {
public:

  ///
  const std::string& isA() const;

  ///
  Factor* process(const Dirichlet<AssignType>* lstPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);

}; // Dirichlet_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/
/**
 * This operator is NIY, since we dont quite see it being used.
 */
template <typename AssignType>
class Dirichlet_ObserveAndReduce: public Operator1< Dirichlet<AssignType> > {
public:

  const std::string& isA() const {
    static const std::string CLASSNAME(
      "Dirichlet_ObserveAndReduce<"+std::string(typeid(AssignType).name())+std::string(">"));
    return CLASSNAME;
  } // isA

  Factor* process(const Dirichlet<AssignType>* lstPtr,
              const emdw::RVIds& variables,
              const emdw::RVVals& assignedVals,
              bool presorted) {
    PRLITE_ASSERT(
      false,
      "Dirichlet_ObserveAndReduce::process NIY since it was not expected to ever be used!");
    return (Dirichlet<AssignType>*) nullptr;
  } // process

}; // Dirichlet_ObserveAndReduce

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/
/**
 * Dampens/smooths the corresponding Dirichlets in the two sets by
 * interpolating their alpha counts.  Returns the total distance
 * between the old and non-dampened new distributions.
 *
 * new = lhs*(1-df) + rhs*df
 */
template <typename AssignType>
class Dirichlet_InplaceWeakDamping: public Operator1< Dirichlet<AssignType> > {
public:
  const std::string& isA() const;
  double inplaceProcess(Dirichlet<AssignType>* lstPtr, const Factor*, double df);
}; // Dirichlet_InplaceWeakDamping

#include "dirichlet.tcc"

#endif // DIRICHLET_HPP
