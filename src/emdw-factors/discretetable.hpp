#ifndef DISCRETETABLE_HPP
#define DISCRETETABLE_HPP

/*******************************************************************************

          AUTHOR: JA du Preez
          DATE:    March 2015
          PURPOSE: DiscreteTable PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "factor.hpp" // Factor
#include "categorical.hpp" // Categorical
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "cat1expf.hpp" // Cat1EXPF

// standard headers
#include <iosfwd>  // istream, ostream
#include <map>

template <typename AssignType> class DiscreteTable_InplaceNormalize;
template <typename AssignType> class DiscreteTable_InplaceMaxNormalize;
template <typename AssignType> class DiscreteTable_InplaceAbsorb;
template <typename AssignType> class DiscreteTable_Absorb;
template <typename AssignType> class DiscreteTable_InplaceAbsorb0UF;
template <typename AssignType> class DiscreteTable_Absorb0UF;
template <typename AssignType> class DiscreteTable_InplaceCancel;
template <typename AssignType> class DiscreteTable_Marginalize;
template <typename AssignType> class DiscreteTable_MaxMarginalize;
template <typename AssignType> class DiscreteTable_ObserveAndReduce;
template <typename AssignType> class DiscreteTable_InplaceDamping;
template <typename AssignType> class DiscreteTable_Sampler;

/**
 * A helper class, essentially there to (fairly transparently) contain
 * the sparse probabilities needed in DiscreteTable, while
 * simultaneously providing the facility to save an additional
 * flag/counts.
 */
class FProb {
public:
  FProb(double theProb=0.0) : prob(theProb) {}
  bool operator==(const FProb& p) const {return prob == p.prob;}
public:
  double prob;
  unsigned flag = 0u;
}; // FProb

// This converts a straight sparse probs map to one that uses FProb's instead
template <typename AssignType>
std::map<std::vector<AssignType>, FProb> convertToFProbMap(const std::map<AssignType, double>& sparseProbs) {
  std::map<std::vector<AssignType>, FProb> theSparseProbs;
  for (auto el : sparseProbs){
    theSparseProbs.emplace_hint(end(theSparseProbs), std::vector<AssignType>{el.first}, el.second);
  } // for
  return theSparseProbs;
} // convertToFProbMap

// ****************************************************************************
// *************************** DiscreteTable **********************************
// ****************************************************************************

/**
 * This is the standard workhorse discrete factor. It can handle a)
 * multiple random variables, b) a sparse representation c) a default
 * prob that may be larger than zero, d) a margin within which sparse
 * values will be converted to default values and e) the ability to
 * specify a relative floor value that limits how low a probability
 * may go as compared to the highest probability in the distribution.
 */
template <typename AssignType>
class DiscreteTable : public Categorical {

  friend class DiscreteTable_InplaceNormalize<AssignType>;
  friend class DiscreteTable_InplaceMaxNormalize<AssignType>;
  friend class DiscreteTable_InplaceAbsorb<AssignType>;
  friend class DiscreteTable_Absorb<AssignType>;
  friend class DiscreteTable_InplaceAbsorb0UF<AssignType>;
  friend class DiscreteTable_Absorb0UF<AssignType>;
  friend class DiscreteTable_InplaceCancel<AssignType>;
  friend class DiscreteTable_Marginalize<AssignType>;
  friend class DiscreteTable_MaxMarginalize<AssignType>;
  friend class DiscreteTable_ObserveAndReduce<AssignType>;
  friend class DiscreteTable_InplaceDamping<AssignType>;
  friend class DiscreteTable_Sampler<AssignType>;

  //============================ Traits and Typedefs ============================
private:
  //typedef std::vector<AssignType> DASS;

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:
  /**
   * @brief Default constructor.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  DiscreteTable(
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
   * The class specific constructor for creating a DiscreteTable
   * factor with given parameter values.
   *
   * @param theVars The ids of the random variables to use in this
   * factor.
   *
   * @param theDomains: The values that each one of these random
   * variables can take on (and in the same order as theVars). The
   * shared pointer is to go easy on memory, many factors will
   * probably use the same RV domains. NOTE: Would have preferred to
   * use a set instead of a vector, but unfortunately general 'absorb'
   * / 'cancel' with a non-zero lhs defProb requires random access.
   *
   * @param theDefProb: The default prob value to to use for all cases
   * not specified in theSparseProbs. For the non-informative case set
   * theDefProb > 0 (1.0 or 1.0/maxSize should do) and keep
   * sparseProbs empty. As with all the other probabilities, defProbs
   * will not be allowed to go lower than the floor value. When
   * sparseProbs specifies all possible assignments (i.e. there are no
   * constraints on defProb), defProb will automatically be set to a
   * floor value (see later), irrespective of the value specified
   * here.
   *
   * @param theSparseProbs A map containing those specific prior probs
   * that differ from the default value. For the non-informative case
   * this map will be empty.
   *
   * @param theMargin Any probs closer (or equal) than this number to
   * defProb will be grouped with the defProb and NOT with the
   * specified exceptions of the sparseProbs. Recommended value about
   * twice as big as theRelFloor.
   *
   * @param theRelFloor The lowest any probability will be allowed to
   * go relative to the highest probability in the distribution.
   * - Do not use a non-zero relFloor for factors that express some
   *   logic where certain assignments are prohibited (i.e. should
   *   have prob 0). In such a situation a non-zero relFloor will in
   *   effect legalize those prohibited assignments, thereby altering
   *   the logic of the model.
   * - However, tables where all possible assignments are allowed can
   *   benefit from lumping the small probabilities under a single
   *   lowest/floor value, thereby reducing the number of sparse
   *   probabilities to explicitly account for.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest - this avoids sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  DiscreteTable(
    const emdw::RVIds& theVars,
    const std::vector< rcptr< std::vector<AssignType> > >& theDomains,
    double theDefProb = 1.0,
    const std::map<std::vector<AssignType>, FProb>& theSparseProbs = std::map<std::vector<AssignType>, FProb>(),
    double theMargin = 0.0,
    double theRelFloor = 0.0,
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

  /// copy and move constructors
  DiscreteTable(const DiscreteTable<AssignType>& p)  = default;
  DiscreteTable(DiscreteTable<AssignType>&& p) = default;

  /// DiscreteTable class destructor
  virtual ~DiscreteTable();

  //========================= Operators and Conversions =========================
public:

  /// Assignment operator
  DiscreteTable& operator=(const DiscreteTable<AssignType>& d) = default;
  DiscreteTable& operator=(DiscreteTable<AssignType>&& d) = default;

  /// == operator
  bool operator==(const DiscreteTable<AssignType>& d) const;

  /// != operator
  bool operator!=(const DiscreteTable<AssignType>& d) const;

  // With a bit of luck we can retire the following
  // Cat1EXPF/DiscreteTable conversions somewhere in the near future.

  // ctor to convert a Cat1EXPF to a DiscreteTable.
  // NOTE: This will use the default factor operators
  DiscreteTable<AssignType>( const Cat1EXPF<AssignType>& c1EXPF)
      : DiscreteTable<AssignType>( {c1EXPF.varId},{c1EXPF.domain}, c1EXPF.defProb, ::convertToFProbMap(c1EXPF.sparseProbs) ) {}

  // conversion of a 1-dim DiscreteTable to a Cat1EXPF
  explicit operator Cat1EXPF<AssignType>() const {

    PRLITE_ASSERT(noOfVars() == 1, "Can only convert a one-dimensional DiscreteTable to a Cat1EXPF: " << noOfVars() );

    const std::map<std::vector<AssignType>, FProb>& theSparseProbs( getSparseProbs() );
    std::map<AssignType, double> sparseProbs;
    for (auto el : theSparseProbs) {
      sparseProbs.emplace_hint(end(sparseProbs), el.first[0], el.second.prob);
    } // for

    return Cat1EXPF<AssignType> (getVar(0), domains[0], getDefaultProb(), sparseProbs);
  } // convert to Cat1EXPF

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /**
   * @brief DiscreteTable class specific configurer. This function
   * can be used to populate an already defined DiscreteTable
   * factor.
   *
   * @param theVars The ids of the random variables to use in this
   * factor.
   *
   * @param theDomains: The values that each one of these random
   * variables can take on (and in the same order as theVars). The
   * shared pointer is to go easy on memory, many factors will
   * probably use the same RV domains. NOTE: Would have preferred to
   * use a set instead of a vector, but unfortunately general 'absorb'
   * / 'cancel' with a non-zero lhs defProb requires random access.
   *
   * @param theDefProb: The default prob value to to use for all cases
   * not specified in theSparseProbs. For the non-informative case set
   * theDefProb > 0 (1.0 or 1.0/maxSize should do) and keep
   * sparseProbs empty. As with all the other probabilities, defProbs
   * will not be allowed to go lower than the floor value. When
   * sparseProbs specifies all possible assignments (i.e. there are no
   * constraints on defProb), defProb will automatically be set to a
   * floor value (see later), irrespective of the value specified
   * here.
   *
   * @param theSparseProbs A map containing those specific prior probs
   * that differ from the default value. For the non-informative case
   * this map will be empty.
   *
   * @param theMargin Any probs closer (or equal) than this number to
   * defProb will be grouped with the defProb and NOT with the
   * specified exceptions of the sparseProbs. If 0.0 then this plays
   * no role, recommended value is 1E-10 (or 1.0/(maxSz*10E4) ).
   *
   * @param theRelFloor The lowest any probability will be allowed to
   * go relative to the highest probability in the distribution.
   *
   * - Use a relFloor=0 in factors expressing some logic where certain
   *   assignments are prohibited (i.e. should have prob 0). In such a
   *   situation a non-zero relFloor will in effect legalize those
   *   prohibited assignments, thereby altering the logic of the
   *   model.
   *
   * - However, tables where all possible assignments are allowed can
   *   benefit from lumping the small probabilities under a single
   *   lowest/floor value, thereby reducing the number of sparse
   *   probabilities to explicitly account for.
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest - this avoids sorting and repacking costs.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  unsigned classSpecificConfigure(
    const emdw::RVIds& theVars,
    const std::vector< rcptr< std::vector<AssignType> > >& theDomains,
    double theDefProb = 1.0,
    const std::map<std::vector<AssignType>, FProb>& theSparseProbs = std::map<std::vector<AssignType>, FProb>(),
    double theMargin = 0.0,
    double theRelFloor = 0.0,
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

  // These functions will be used in the future to enable users to
  // write DiscreteTable factors to file and to read DiscreteTable
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
  virtual DiscreteTable<AssignType>* copy(const emdw::RVIds& newVars = {}, bool presorted = false) const;

  /**
   * Make copy of this factor, but set the potentials to the default
   * non informative factor values. The client has to manage the
   * lifetime of this pointer - putting it in a rcptr will be a good
   * idea.
   */
  virtual DiscreteTable<AssignType>* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /**
   * Checks if two DiscreteTables are identical.
   */
  virtual bool isEqual(const Factor* rhsPtr) const ;

  /// one-sided KL to uniform distribution
  virtual double distanceFromVacuous() const;

  /**
   * Returns the number of RV ids in a DiscreteTable.
   */
  virtual unsigned noOfVars() const;

  /**
   * Returns the vector of RV ids.
   */
  virtual emdw::RVIds getVars() const;

  /**
   * returns the RV id of the discretetable var.
   */
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //====================== Other Virtual Member Functions =======================
public:

  /**
   * Measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference. Currently based
   * on the one-sided Kulbach Leibler distance.
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
   * Assigns a probability/potential value to a specific/single RV
   * outcome.  *NOTE* The potential value is assigned directly, it is
   * assumed to be already properly scaled.
   *
   * @param theVars The RV ids in the order that the values in
   * theirVals assumes them to be.
   *
   * @param theirVals The values assigned to this random vector.
   *
   * @param theProb The probability/potential value to assign to the
   * above random vector
   *
   * @param presorted Set this to true if theVars are sorted from
   * smallest to biggest and avoid sorting and repacking costs.
   */
  virtual void setEntry(const emdw::RVIds& theVars,
                        const emdw::RVVals& theirVals,
                        double theProb,
                        bool presorted = false);

  // returns the partition function value
  virtual double sumZ() const;

  //======================= Non-virtual Member Functions ========================
public:

  /// inplace absorb a Cat1EXPF
  void inplaceAbsorbCat1EXPF(const rcptr< Cat1EXPF<AssignType> >& c1Ptr);

  /// Marginalize to a Cat1EXPF
  Cat1EXPF<AssignType>* marginalizeToCat1EXPF(emdw::RVIdType varId) const;

  /** If the default prob is not zero, fill in all the missing
   * probablities in sparseProbs (equal to the default value) and
   * then change the default probability to zero. Careful with this,
   * the resulting table might be huge.
   */
  void makeDefaultZero() const;

  /** Sets the margin to the provided input.
   */
  void setMargin( double theMargin) {margin=theMargin;} // setMargin

  /**
   * Adds uniform noise to all elements in the sparseTable (only).
   * @param upperLim: the added noise will be uniformly spread on the interval (0:upperLim]
   */
  void addNoise(double upperLim);

  /// zero all probabilities less than the threshold
  void squashSmallProbs(double thresh);

  /// zero all probabilities less than the threshold
  void nonZeroToUnity();

  /**
   * Calculates the entropy for the distribution (even if
   * unnormalized).
   */
  double entropy() const;

  /**
   * A more direct version of potentialAt
   *
   * @param theAss The specific assignment to the random variables for
   *   which the probability must be returned. It is assumed that the
   *   corresponding random variables are sorted according to
   *   ascending RVId.
   *
   * @return The (unnormalised) probabilty at the RV assignment given
   *   by theAss.
   */
  double getProb(const std::vector<AssignType>& theAss) const;

  // returns the default probability
  double getDefaultProb() const;

  // returns the map with sparse probability values
  const std::map<std::vector<AssignType>, FProb>& getSparseProbs() const;

  /// returns the cardinality of the an RV.
  unsigned getCardinality(unsigned varNo) const;

  /// replaces various extra settings (those following after
  /// sparseProbs in constructor) with those of d.
  void copyExtraSettingsFrom(const DiscreteTable<AssignType>& d);

  //================================== Friends ==================================

public:

  /// output
  template<typename T>
  friend std::ostream& operator<<( std::ostream& file,
                                   const DiscreteTable<AssignType>& n );

  /// input
  template<typename T>
  friend std::istream& operator>>( std::istream& file,
                                   DiscreteTable<AssignType>& n );

  //=============================== Data Members ================================
private:

  //+++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++ definitive data members ++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++

  // definitive variables
  // the ids of all the RVs in use here
  emdw::RVIds vars;

public:
  // the values those RVs can take on
  std::vector< rcptr< std::vector<AssignType> > > domains;

private:
  // default probability for all assignments not in sparseProbs.
  mutable double defProb;
public:
  // all non-default probabilities. empty for uniform
  mutable std::map<std::vector<AssignType>, FProb> sparseProbs;
private:
  // probs closer than this to defProb will become defProb
  double margin;
  // the lowest a probability can go. Note that this is distinct from
  // defProb (which can be bigger than this).
  double relFloor;

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

  //+++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++ derived data members +++++++++++++++
  //+++++++++++++++++++++++++++++++++++++++++++++++++

  std::vector<size_t> cards;  // the RV cardinalities
  size_t maxSz;  // product of cardinalities of all RVs
}; // DiscreteTable


// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a DiscreteTable.
 */
template <typename AssignType>
class DiscreteTable_Sampler: public SampleOperator< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  emdw::RVVals sample(const DiscreteTable<AssignType>* lstPtr);
}; // DiscreteTable_Sampler

// ***************************************************************************
// ******************************** Family 1 *********************************
// ******************************** Normalization ****************************
// ***************************************************************************

/**
 * This operator does an inplace normalization of the DiscreteTables.
 */
template <typename AssignType>
class DiscreteTable_InplaceNormalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(DiscreteTable<AssignType>* lstPtr);
}; // DiscreteTable_InplaceNormalize

/**
 * This operator returns a normalized version of the DiscreteTables.
 */
template <typename AssignType>
class DiscreteTable_Normalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr);
}; // DiscreteTable_Normalize

/**
 * This operator does an inplace normalization of the DiscreteTables.
 */
template <typename AssignType>
class DiscreteTable_InplaceMaxNormalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(DiscreteTable<AssignType>* lstPtr);
}; // DiscreteTable_InplaceMaxNormalize

/**
 * This operator returns a normalized version of the DiscreteTables.
 */
template <typename AssignType>
class DiscreteTable_MaxNormalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr);
}; // DiscreteTable_MaxNormalize

// ***************************************************************************
// ******************************** Family 2 *********************************
// ******************************** Absorb ***********************************
// ***************************************************************************

/**
 * This operator returns the inplace product of the lhs DiscreteTable
 * with another one.
 */
template <typename AssignType>
class DiscreteTable_InplaceAbsorb: public Operator1< DiscreteTable<AssignType> > {
public:

  /// identifies this class
  const std::string& isA() const;

  /**
   * The underlying (non-virtual) multiplication engine.
   *
   * @param lhs The factor to be updated with an inplace multiply.
   *
   * @param rhs The factor we are multiplying with. Its random
   * variables must be a subset of those in lhs.
   *
   * @param l2i the indices of the RVs in lhs that are also common to
   * those in rhs.
   *
   * @param lni the indices of the RVs unique to lhs (i.e. not in the
   * intersection)
   *
   * @param i2l  where in the RV union is l2i[k] located
   *
   * @param ni2l where in the union is lni[k] located
   */
  void inplaceAbsorb(
    DiscreteTable<AssignType>& lhs,
    const DiscreteTable<AssignType>& rhs,
    const std::vector<size_t>& l2i,
    const std::vector<size_t>& lni,
    const std::vector<size_t>& i2l,
    const std::vector<size_t>& ni2l);

  /// and the version we can access polymorphically
  void inplaceProcess(DiscreteTable<AssignType>* lstPtr, const Factor*);

}; // DiscreteTable_InplaceAbsorb

/**
 * This operator returns the product of the lhs DiscreteTable with
 * another one.
 */
template <typename AssignType>
class DiscreteTable_Absorb: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr, const Factor*);
}; // DiscreteTable_Absorb

// ***************************************************************************
// ******************************** Family 2 *********************************
// ******************************** Absorb0UF ********************************
// ***************************************************************************

/**
 * This operator returns the inplace product of the lhs DiscreteTable
 * with another one.
 */
template <typename AssignType>
class DiscreteTable_InplaceAbsorb0UF: public Operator1< DiscreteTable<AssignType> > {
public:

  /// identifies this class
  const std::string& isA() const;

  /**
   * The underlying (non-virtual) multiplication engine.
   *
   * @param lhs The factor to be updated with an inplace multiply.
   *
   * @param rhs The factor we are multiplying with. Its random
   * variables must be a subset of those in lhs.
   *
   * @param l2i the indices of the RVs in lhs that are also common to
   * those in rhs.
   *
   * @param lni the indices of the RVs unique to lhs (i.e. not in the
   * intersection)
   *
   * @param i2l  where in the RV union is l2i[k] located
   *
   * @param ni2l where in the union is lni[k] located
   */
  void inplaceAbsorb(
    DiscreteTable<AssignType>& lhs,
    const DiscreteTable<AssignType>& rhs,
    const std::vector<size_t>& l2i,
    const std::vector<size_t>& lni,
    const std::vector<size_t>& i2l,
    const std::vector<size_t>& ni2l);

  /// and the version we can access polymorphically
  void inplaceProcess(DiscreteTable<AssignType>* lstPtr, const Factor*);

}; // DiscreteTable_InplaceAbsorb0UF

/**
 * This operator returns the product of the lhs DiscreteTable with
 * another one.
 */
template <typename AssignType>
class DiscreteTable_Absorb0UF: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr, const Factor*);
}; // DiscreteTable_Absorb0UF

// ***************************************************************************
// ******************************** Family 2 *********************************
// ******************************** Cancel ***********************************
// ***************************************************************************

/**
 * This operator inplace divides by the second DiscreteTable.
 */
template <typename AssignType>
class DiscreteTable_InplaceCancel: public Operator1< DiscreteTable<AssignType> > {
public:

  /// identifies this class
  const std::string& isA() const;
  /**
   * The underlying (non-virtual) division engine.
   *
   * @param lhs The factor to be updated with an inplace divide.
   *
   * @param rhs The factor we are multiplying with. Its random
   * variables must be a subset of those in lhs.
   *
   * @param l2i the indices of the RVs in lhs that are also common to
   * those in rhs.
   *
   * @param lni the indices of the RVs unique to lhs (i.e. not in the
   * intersection)
   *
   * @param i2l  where in the RV union is l2i[k] located
   *
   * @param ni2l where in the union is lni[k] located
   */
  void inplaceCancel(
    DiscreteTable<AssignType>& lhs,
    const DiscreteTable<AssignType>& rhs,
    const std::vector<size_t>& l2i,
    const std::vector<size_t>& lni,
    const std::vector<size_t>& i2l,
    const std::vector<size_t>& ni2l);

  /// and the version we can access polymorphically
  void inplaceProcess(DiscreteTable<AssignType>* lstPtr, const Factor*);

}; // DiscreteTable_InplaceCancel

/**
 * This operator returns the result from dividing by the second DiscreteTable.
 */
template <typename AssignType>
class DiscreteTable_Cancel: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr, const Factor*);
}; // DiscreteTable_Cancel

// ***************************************************************************
// ******************************** Family 3 *********************************
// ******************************** Marginalize ******************************
// ***************************************************************************
/**
 * Marginalize the distribution
 */
template <typename AssignType>
class DiscreteTable_Marginalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
//  Factor* process(const DiscreteTable<AssignType>* lstPtr);
  Factor* process(const DiscreteTable<AssignType>* lstPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // DiscreteTable_Marginalize

/**
 * Marginalize the distribution
 */
template <typename AssignType>
class DiscreteTable_MaxMarginalize: public Operator1< DiscreteTable<AssignType> > {
public:
  const std::string& isA() const;
//  Factor* process(const DiscreteTable<AssignType>* lstPtr);
  Factor* process(const DiscreteTable<AssignType>* lstPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // DiscreteTable_MaxMarginalize

// ***************************************************************************
// ******************************** Family 4 *********************************
// ******************************** ObserveAndReduce *************************
// ***************************************************************************
/**
 * After observing some variables, the scope of the resultant factor
 * reduces to ONLY that of the unobserved variables, making further
 * computations potentially much more efficient. These may validly
 * also include variables not part of this Factor.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 */
template <typename AssignType>
class DiscreteTable_ObserveAndReduce: public Operator1< DiscreteTable<AssignType> > {
public:

  const std::string& isA() const;
  Factor* process(const DiscreteTable<AssignType>* lstPtr,
                  const emdw::RVIds& variables,
                  const emdw::RVVals& assignedVals,
                  bool presorted);
}; // DiscreteTable_ObserveAndReduce

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************
/**
 * Dampens/smooths two corresponding DiscreteTables and returns the
 * total (symmetric / two-sided KL) distance between them.
 *
 * new = lhs*(1-oldWeight) + rhs*oldWeight
 */
template <typename AssignType>
class DiscreteTable_InplaceDamping: public Operator1< DiscreteTable<AssignType> > {
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
    DiscreteTable<AssignType>* newPtr,
    const Factor* oldPtr,
    double oldWeight);

}; // DiscreteTable_InplaceDamping

#include "discretetable.tcc"

#endif // DISCRETETABLE_HPP
