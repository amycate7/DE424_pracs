#ifndef POLYA_HPP
#define POLYA_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2018
          PURPOSE: Polya PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "discretetable.hpp" // DiscreteTable
#include "dirichlet.hpp" // Dirichlet
#include "factoroperator.hpp" // FactorOperator
#include "emdw.hpp"
#include "anytype.hpp"

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream
#include <map>     // map

// forward declarations
template <typename AssignType> class Polya_InplaceNormalize;
template <typename AssignType> class Polya_Marginalize;
template <typename AssignType> class Polya_InplaceAbsorb;
template <typename AssignType> class Polya_InplaceCancel;
template <typename AssignType> class Polya_ObserveAndReduce;
template <typename AssignType> class Polya_InplaceDamping;
template <typename AssignType> class Polya_Sampler;

/******************************************************************************/
/****************************** Polya ***********************************/
/******************************************************************************/

/**
 * This class models the distribution of a DiscreteTable random
 * variables given a Dirichlet prior.
 *
 *
 * From the documentation for the Dirichlet class, note that although
 * a Dirichlet actually is a random VECTOR (with a dimension/component
 * theta_0 upto theta_{K-1} per assignment value of the DiscreteTable
 * random VALUE, we will always refer to the vector as a whole and it
 * therefore is identified with only a single RVId.
 *
 * Conceptually we might also want to supply a Dirichlet prior to a
 * composition of Categorical variables (a typical DiscreteTable
 * record for instance). We have two ways of doing this:
 * 1) The first is to first uniquely map all the value combinations
 *   for those multiple Categoricals to a single Categorical variable
 *   and then model that variable in this class.
 * 2) The other option is to specify the whole collection of
 *   Categorical variables. The various joint outcomes of these
 *   Categorical RVs are then treated as if a single outcome for a
 *   composite variable.
 *
 * NOTE: It turns out that specifying multiple Categorical RVs is not
 * feasable in Conditional{Dirichlet,Polya} (elsewhere). That calls
 * our suppport for multiple RVs here in question. The same result can
 * probably be achieved while allowing only a single RV. If necessary
 * the client code can always map multiple RVs to single RVs.
 *
 * NOTE: We are ONLY going going to support Dirichlet or DiscreteTable
 * messages entering or leaving a Polya.
 */
template <typename AssignType>
class Polya : public Factor {

  // we'll want to reduce these dependencies later on
  friend class Polya_InplaceNormalize<AssignType>;
  friend class Polya_Marginalize<AssignType>;
  friend class Polya_InplaceAbsorb<AssignType>;
  friend class Polya_InplaceCancel<AssignType>;
  friend class Polya_ObserveAndReduce<AssignType>;
  friend class Polya_InplaceDamping<AssignType>;
  friend class Polya_Sampler<AssignType>;

  //============================ Traits and Typedefs ============================
private:

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

  public:

  /**
   * Default CTOR.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  Polya(
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

  /// copy ctor
  Polya(const Polya<AssignType>& st);
  /// move ctor
  Polya(Polya<AssignType>&& st);


  /**
   * Class Specific Constructor.
   *
   * @param theDirVarId: The (single) id of the Dirichlet RV. NOTE: To
   * keep things simple we require that this is a value lower than
   * theCatVarId.  Also, although a Dirichlet RV implicitly really is a
   * Random Vector, we are only allocating a single id for it as a
   * whole.
   *
   * @param theCatVarId: The id of the associated Categorical RV.
   *
   * @param theDomain: The values that the Categorical RV can take on.
   *
   * @param theDirPtr (normally not instantiated): the Dirichlet
   * distribution. This only gets instantiated as the result of
   * multiplying the conditional Polya by a Dirichlet (making it a
   * joint Polya distribution).
   *
   * @param theCatPtr (normally not instantiated): the DiscreteTable
   * distribution. This only gets instantiated as the result of
   * multiplying a Polya by a DiscreteTable resulting in a
   * ItWhoShallNotBeNamed.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
  Polya(
    emdw::RVIdType theDirVarId,  // theDirVarId < theCatVarId
    emdw::RVIdType theCatVarId,
    rcptr< std::vector<AssignType> > theDomain,
    rcptr< Dirichlet<AssignType> > theDirPtr = 0,
    rcptr< DiscreteTable<AssignType> > theCatPtr = 0,
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
  virtual ~Polya();

  //========================= Operators and Conversions =========================
public:

  ///
  // ???!!!
  Polya& operator=(const Polya<AssignType>& d) = default;
  Polya& operator=(Polya<AssignType>&& d) = default;

  ///
  bool operator==(const Polya<AssignType>& d) const;

  ///
  bool operator!=(const Polya<AssignType>& d) const;

  //================================= Iterators =================================

  //=========================== Inplace Configuration ===========================
public:

  /**
   * Class Specific Configure
   *
   * @param theDirVarId: The (single) id of the Dirichlet RV. NOTE: To
   * keep things simple we require that this is a value lower than
   * theCatVarId.  Also, although a Dirichlet RV implicitly really is a
   * Random Vector, we are only allocating a single id for it as a
   * whole.
   *
   * @param theCatVarId: The id of the associated Categorical RV.
   *
   * @param theDomain: The values that the Categorical RV can take on.
   *
   * @param theDirPtr (normally not instantiated): the Dirichlet
   * distribution. This only gets instantiated as the result of
   * multiplying the conditional Polya by a Dirichlet (making it a
   * joint Polya distribution).
   *
   * @param theCatPtr (normally not instantiated): the DiscreteTable
   * distribution. This only gets instantiated as the result of
   * multiplying a Polya by a DiscreteTable resulting in a
   * ItWhoShallNotBeNamed.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available. Take care to specify them in
   * their proper order.
   */
unsigned classSpecificConfigure(
  emdw::RVIdType theDirVarId,  // theDirVarId < theCatVarId
  emdw::RVIdType theCatVarId,
  rcptr< std::vector<AssignType> > theDomain,
  rcptr< Dirichlet<AssignType> > theDirPtr = 0,
  rcptr< DiscreteTable<AssignType> > theCatPtr = 0,
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

  /**
   * With no params it simply returns a copy of the existing Factor,
   * with params it returns an identical Factor, but with variables
   * replaced by another set.
   *
   * @param newVars These variables will replace the original
   * variables. We assume that the first is the Dirichlet variable and
   * the second is the DiscreteTable variable.
   *
   * @param presorted is irrelevant.
   */
  virtual Polya<AssignType>* copy(const emdw::RVIds& newVars = {},
                                        bool presorted = false) const;

  virtual Polya<AssignType>*
  vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from (the top-level) Factor. They return uniqptr so
  /// that we can chain them.

  //------------------Family 0: eg Sample

  virtual emdw::RVVals sample(Sampler* procPtr) const;

  //------------------Family 1 eg Normalize

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceNormalize(FactorOperator* procPtr);

  /**
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual uniqptr<Factor> normalize(FactorOperator* procPtr) const;

  //------------------Family 2 eg Absorb,Cancel

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceAbsorb(const Factor* rhsPtr,
                            FactorOperator* procPtr);

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual uniqptr<Factor> absorb(const Factor* rhsPtr,
                                FactorOperator* procPtr) const;

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual void inplaceCancel(const Factor* rhsPtr,
                            FactorOperator* procPtr);

  /**
   * @param rhsPtr is the 'other' Factor in this operation.
   *
   * @param procPtr. Can be used to specify an optional operator for
   * this operation; if NULL, this routine will use the operator
   * provided by the Factor object itself.
   */
  virtual uniqptr<Factor> cancel(const Factor* rhsPtr,
                                FactorOperator* procPtr) const;

  //------------------Family 3 eg Marginalize

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
  virtual uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep,
                                      bool presorted,
                                      const Factor* peekAheadPtr,
                                     FactorOperator* procPtr) const;

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
                                           const emdw::RVVals& assignedVals,
                                           bool presorted,
                                           FactorOperator* procPtr) const;

  //------------------Family 5: eg inplaceDampen

  /**
   * This operator is NIY, since we dont quite see it being used.
   */
  double inplaceDampen(const Factor* oldMsg,
                       double df,
                       FactorOperator* procPtr);

  //===================== Required Virtual Member Functions =====================
public:

  ///
  virtual bool isEqual(const Factor* rhsPtr) const;

  /// For the while being, use the default implementation of this member
  double distanceFromVacuous() const; // {return Factor::distanceFromVacuous();}

  ///
  virtual unsigned noOfVars() const;

  /// Note, we do not know whether the Dirichlet or the DiscreteTable
  /// comes first.
  virtual emdw::RVIds getVars() const;

  /// returns the sorted variable ids. NOTE we do not know whether the
  /// Dirichlet or the DiscreteTable variable comes first
  virtual emdw::RVIdType getVar(unsigned varNo) const;


  //====================== Other Virtual Member Functions =======================
public:

  /**
   * measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Non-virtual Member Functions ========================
public:

  uniqptr< Dirichlet<AssignType> > getDirichlet() const;
  uniqptr< DiscreteTable<AssignType> > getDiscreteTable() const;
  //emdw::RVIdType getDirichletId() const;
  //emdw::RVIdType getCategoricalId() const;

  //================================== Friends ==================================

public:

  /// output
  template<typename T>
  friend std::ostream& operator<<( std::ostream& file,
                                   const Polya<T>& n );

  /// input
  template<typename T>
  friend std::istream& operator>>( std::istream& file,
                                   Polya<T>& n );

  //=============================== Data Members ================================
private:

  // for now we are going to enforce a strict ordering on the variable ids
  emdw::RVIdType dirVarId; // dirVarId < catVarId
  emdw::RVIdType catVarId;

  // --- the values that the Categorical RV can take on
  rcptr< std::vector<AssignType> > domain;

  rcptr< Dirichlet<AssignType> > dirPtr;
  rcptr< DiscreteTable<AssignType> > catPtr;

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

}; // Polya


// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a Polya.
 */
template <typename AssignType>
class Polya_Sampler: public SampleOperator< Polya<AssignType> > {
public:
  const std::string& isA() const;
  /// first the sampled probs, then the sampled categorical
  emdw::RVVals sample(
    const Polya<AssignType>* lhsPtr);
}; // Polya_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

/**
 * This operator does an inplace normalization of a Polya.
 */
template <typename AssignType>
class Polya_InplaceNormalize: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Polya<AssignType>* lstPtr);
}; // Polya_InplaceNormalize

/**
 * This operator returns a normalized version of a Polya.
 */
template <typename AssignType>
class Polya_Normalize: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Polya<AssignType>* lstPtr);
}; // Polya_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

/**
 * This operator inplace multiplies a Polya by either
 * a) another Polya, or b) a Dirichlet, or c) a DiscreteTable.
 */
template <typename AssignType>
class Polya_InplaceAbsorb: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Polya<AssignType>* lhsPtr,
                      const Factor* rhsFPtr);
}; // Polya_InplaceAbsorb

/**
 * This operator returns the factor product of a Polya by one of
 * a) another Polya, or b) a Dirichlet, or c) a DiscreteTable.
 */
template <typename AssignType>
class Polya_Absorb: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
private:   // not available due to conditional requiring side effects
  Factor* process(const Polya<AssignType>* lhsPtr,
                  const Factor* rhsFPtr){
    PRLITE_ASSERT(
      false,
      "Polya_Absorb::process NIY since conditional requires side effects!");
    return (Polya<AssignType>*) nullptr;
  } // process

}; // Polya_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

/**
 * This operator inplace divides a Polya by either a) another
 * Polya, or b) a Dirichlet, or c) a DiscreteTable.
 */
template <typename AssignType>
class Polya_InplaceCancel: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(Polya<AssignType>* lhsPtr,
                      const Factor* rhsFPtr);
}; // Polya_InplaceCancel

/**
 * This operator returns the factor division of a Polya by one of
 * a) another Polya, or b) a Dirichlet, or c) a DiscreteTable.
 */
template <typename AssignType>
class Polya_Cancel: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
private:   // not available due to conditional requiring side effects
  Factor* process(const Polya<AssignType>* lhsPtr,
                  const Factor* rhsFPtr){
    PRLITE_ASSERT(
      false,
      "Polya_Cancel::process NIY since conditional requires side effects!");
    return (Polya<AssignType>*) nullptr;
  } // process

}; // Polya_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

/**
 * This operator marginalizes a Polya to either a) a DiscreteTable
 * or b) a Dirichlet, or c) (trivially) another Polya.
 */
template <typename AssignType>
class Polya_Marginalize: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Polya<AssignType>* lhsPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // Polya_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

/**
 * After observing some variables, the scope of the resultant factor
 * reduces to ONLY that of the unobserved variables, making further
 * computations potentially much more efficient. These may validly
 * also include variables not part of this Factor.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 *
 * Currently we only support observing the DiscreteTable variable.
 */
template <typename AssignType>
class Polya_ObserveAndReduce: public Operator1< Polya<AssignType> > {
public:
  const std::string& isA() const;
  Factor* process(const Polya<AssignType>* lhsPtr,
                const emdw::RVIds& variables,
                const emdw::RVVals& assignedVals,
                bool presorted);
}; // Polya_ObserveAndReduce

/*****************************************************************************/
/********************************** Family 5 *********************************/
/********************************** Damping **********************************/
/*****************************************************************************/

/**
 * This operator is NIY, since we dont quite see it being used.
 */
template <typename AssignType>
class Polya_InplaceDamping: public Operator1< Polya<AssignType> > {
public:

  const std::string& isA() const{
    static const std::string CLASSNAME(
      "Polya_InplaceDamping<"+std::string(typeid(AssignType).name())+std::string(">"));
    return CLASSNAME;
  } // isA

  double inplaceProcess(Polya<AssignType>* lhsPtr,
                        const Factor* rhsFPtr,
                        double df) {
    PRLITE_ASSERT(false, "Polya_InplaceDamping::inplaceProcess NIY since it was not expected to ever be used!");
    return 0.0;
  } // inplaceProcess

}; // Polya_InplaceDamping


#include "polya.tcc"

#endif // POLYA_HPP
