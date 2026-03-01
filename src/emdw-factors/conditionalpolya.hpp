#ifndef CONDITIONALPOLYA_HPP
#define CONDITIONALPOLYA_HPP

/*******************************************************************************

          AUTHORS: JA du Preez

          DATE: October 2018

          PURPOSE: ConditionalPolya PGM potential

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "emdw.hpp"
#include "factor.hpp" // Factor
#include "factoroperator.hpp" // FactorOperator
#include "anytype.hpp"
#include "dirichletset.hpp"
#include "discretetable.hpp"
#include "oddsandsods.hpp" // ValidValue
// patrec headers

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

// Some forward declarations
template <typename CondType, typename DirType> class ConditionalPolya_InplaceNormalize;
template <typename CondType, typename DirType> class ConditionalPolya_InplaceMaxNormalize;
template <typename CondType, typename DirType> class ConditionalPolya_InplaceAbsorb;
template <typename CondType, typename DirType> class ConditionalPolya_Absorb;
template <typename CondType, typename DirType> class ConditionalPolya_InplaceCancel;
template <typename CondType, typename DirType> class ConditionalPolya_Marginalize;
template <typename CondType, typename DirType> class ConditionalPolya_MaxMarginalize;
template <typename CondType, typename DirType> class ConditionalPolya_ObserveAndReduce;
template <typename CondType, typename DirType> class ConditionalPolya_InplaceDamping;
template <typename CondType, typename DirType> class ConditionalPolya_Sampler;


/************************************************************************/
/************************ ConditionalPolya ******************************/
/************************************************************************/

/**
 * A ConditionalPolya $p(x|z,phi_z)$ where $x$ may
 * possibly be latent, and and $z$ is the latent variable controlling
 * which Dirichlet to use of $x$. Where-as $x$ is a single Categorical
 * RV, we also allow $z$ to be a combination of several Categorical
 * RVs.
 *
 * NOTE: For the time being we limit the Categorical distributions to
 * DiscreteTable only.
 */
template <typename CondType, typename DirType>
class ConditionalPolya : public Factor {

  //============================ Traits and Typedefs ============================
public:

  //======================== Constructors and Destructor ========================
private:

  // private helper constructor to set up factor operators
  void supplementWithDefaultOperators();

public:

  /**
   * @brief Class Specific Constructor
   *
   * @param theDVar The (single) id for the
   * DirichletSet. Must be smaller than theZVar.
   *
   * @param theZVar The (single) id for the (conditioned on)
   * Categorical RV dictating which Dirichlet to use. Must be smaller
   * than theXVar.
   *
   * @param theXVar The (single) id for the (target) Categorical
   * RV being described by this prior.
   *
   * @param theZDomain The values that the (conditioned on) z RV can
   * take on.
   *
   * @param theXDomain The values that the (target) x RV can
   * take on.
   *
   * @param theObservedX (optional) If instantiated gives the observed
   * value of the (target) Categorical RV.
   *
   * @param theCDirPtr (normally not instantiated) The
   * DirichletSet to use as prior.
   *
   * @param theZCatPtr (normally not instantiated) The DiscreteTable
   * distribution for the z RV.
   *
   * @param theXCatPtr (normally not instantiated) The DiscreteTable
   * distribution for the y RV.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  ConditionalPolya(
    emdw::RVIdType theDVar,
    emdw::RVIdType theZVar,
    emdw::RVIdType theXVar,
    rcptr< std::vector<CondType> > theZDomain,
    rcptr< std::vector<DirType> > theXDomain,
    ValidValue<DirType> theObservedX = ValidValue<DirType>(),
    rcptr< DirichletSet<CondType, DirType> > theCDirPtr = 0,
    rcptr< DiscreteTable<CondType> > theZCatPtr = 0,
    rcptr< DiscreteTable<DirType> > theXCatPtr = 0,
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
   * @brief Default constructor.
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  ConditionalPolya(
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
  ConditionalPolya(const ConditionalPolya<CondType, DirType>& p);

  /// move ctor
  ConditionalPolya(ConditionalPolya<CondType, DirType>&& p);

  ///
  virtual ~ConditionalPolya();

  //========================= Operators and Conversions =========================
public:

  /// TIP: If you have pointers of some type in the class, you should
  /// supply these operators. Below we assume this to not be the case.
  ConditionalPolya& operator=(const ConditionalPolya<CondType, DirType>& d) = default;
  ConditionalPolya& operator=(ConditionalPolya<CondType, DirType>&& d) = default;

  ///
  bool operator==(const ConditionalPolya<CondType, DirType>& d) const;

  ///
  bool operator!=(const ConditionalPolya<CondType, DirType>& d) const;

  //================================= Iterators =================================



  //=========================== Inplace Configuration ===========================
public:

  /**
   * @brief Class Specific Configurer.
   *
   * @param theDVar The (single) id for the
   * DirichletSet. Must be smaller than theZVar.
   *
   * @param theZVar The (single) id for the (conditioned on)
   * Categorical RV dictating which Dirichlet to use. Must be smaller
   * than theXVar.
   *
   * @param theXVar The (single) id for the (target) Categorical
   * RV being described by this prior.
   *
   * @param theZDomain The values that the (conditioned on) z RV can
   * take on.
   *
   * @param theXDomain The values that the (target) y RV can
   * take on.
   *
   * @param theObservedX (optional) If instantiated gives the observed
   * value of the (target) Categorical RV.
   *
   * @param theCDirPtr (normally not instantiated) The
   * DirichletSet to use as prior.
   *
   * @param theZCatPtr (normally not instantiated) The DiscreteTable
   * distribution for the z RV.
   *
   * @param theXCatPtr (normally not instantiated) The DiscreteTable
   * distribution for the y RV
   *
   * A list of FactorOperators, if =0 a default version will be
   * automatically supplied if available.
   */
  unsigned classSpecificConfigure(
    emdw::RVIdType theDVar,
    emdw::RVIdType theZVar,
    emdw::RVIdType theXVar,
    rcptr< std::vector<CondType> > theZDomain,
    rcptr< std::vector<DirType> > theXDomain,
    ValidValue<DirType> theObservedX = ValidValue<DirType>(),
    rcptr< DirichletSet<CondType, DirType> > theCDirPtr = 0,
    rcptr< DiscreteTable<CondType> > theZCatPtr = 0,
    rcptr< DiscreteTable<DirType> > theXCatPtr = 0,
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

  /// input from stream
  virtual std::istream& txtRead(std::istream& file);
  /// output to stream
  virtual std::ostream& txtWrite(std::ostream& file) const;

  //===================== Required Factor Operations ============================
public:

  /// The preferred way to (polymorphically) access these functions is
  /// via the correspondingly named, but rcptr-based routines
  /// available from Factor (at the top-level). They return uniqptr so
  /// that we can chain them. Boilerplate.

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
   * @param newVars These variables should be sorted and will replace,
   * one-for-one AND IN THIS ORDER, the original SORTED set.
   *
   * @param presorted Not used.
   */
  virtual ConditionalPolya<CondType, DirType>* copy(
    const emdw::RVIds& newVars = {},
    bool presorted = false) const;

  /**
   * make a copy of this factor, but set the potentials to the default
   * unity factor values. The client has to manage the lifetime of
   * this pointer - putting it in a rcptr will be a good idea.
   *
   * NOTE 1: Variable selection is ignored.
   * NOTE 2: Currently we retain the status of observedX. Not sure if
   * this is right.
   */
  virtual ConditionalPolya<CondType, DirType>* vacuousCopy(
    const emdw::RVIds& selectedVars = {},
    bool presorted = false) const;

  /// note rhs is a raw pointer
  bool isEqual(const Factor* rhsPtr) const;

  /// TIP: This default implementation should do the trick, but do
  /// provide your own it if necessary.
  double distanceFromVacuous() const {return Factor::distanceFromVacuous();}

  ///
  virtual unsigned noOfVars() const;

  /// returns the ids of the variables the factor is defined on
  virtual emdw::RVIds getVars() const;

  ///
  virtual emdw::RVIdType getVar(unsigned varNo) const;

  //====================== Other Virtual Member Functions =======================
public:

  /**
   * Measures the distance between two factors. zero if equal, else
   * positive and growing with degree of difference. Typically based
   * on KL divergence.
   */
  virtual double distance(const Factor* rhsPtr) const;

  //======================= Non-virtual Member Functions ========================
public:
  void absorbCDir(const DirichletSet<CondType, DirType>* theCDirPtr);
  void cancelCDir(const DirichletSet<CondType, DirType>* theCDirPtr);
  void absorbZCat(const DiscreteTable<CondType>* theZCatPtr);
  void cancelZCat(const DiscreteTable<CondType>* theZCatPtr);
  void absorbXCat(const DiscreteTable<DirType>* theXCatPtr);
  void cancelXCat(const DiscreteTable<DirType>* theXCatPtr);

  ///
  ValidValue<DirType> getObservedX() const;

  /// nRows, nCols
  std::pair<size_t,size_t> volDims() const;

  /// volume for z-row, y-col
  double vol(size_t rZ, size_t cY) const;

private:
  void updateVols();
  //================================== Friends ==================================
public:

  // It is expedient to allow factor operators direct access to the
  // private data of its associated factor.
  friend class ConditionalPolya_InplaceNormalize<CondType, DirType>;
  friend class ConditionalPolya_InplaceMaxNormalize<CondType, DirType>;
  friend class ConditionalPolya_InplaceAbsorb<CondType, DirType>;
  friend class ConditionalPolya_Absorb<CondType, DirType>;
  friend class ConditionalPolya_InplaceCancel<CondType, DirType>;
  friend class ConditionalPolya_Marginalize<CondType, DirType>;
  friend class ConditionalPolya_MaxMarginalize<CondType, DirType>;
  friend class ConditionalPolya_ObserveAndReduce<CondType, DirType>;
  friend class ConditionalPolya_InplaceDamping<CondType, DirType>;
  friend class ConditionalPolya_Sampler<CondType, DirType>;

  /// output
  template<typename T1, typename T2>
  friend std::ostream& operator<<( std::ostream& file,
                                   const ConditionalPolya<T1,T2>& n );

  /// input
  template<typename T1, typename T2>
  friend std::istream& operator>>( std::istream& file,
                                   ConditionalPolya<T1,T2>& n );

  //=============================== Data Members ================================
private:

  // for now we are going to enforce a strict ordering on the variable ids
  emdw::RVIdType dVar;     // the DirichletSet variable
  emdw::RVIdType zVar;     // the latent variable choosing which Dirichlet to use
  emdw::RVIdType xVar;     // the variable for which the Dirichlet is a prior

  // --- the values that the z RV can take on
  rcptr< std::vector<CondType> > zDomain;

  // --- the values that the y RV can take on
  rcptr< std::vector<DirType> > xDomain;

  ValidValue<DirType> observedX;

  rcptr< DirichletSet<CondType, DirType> > cDirPtr;
  rcptr< DiscreteTable<CondType> > zCatPtr;  // only a single RV
  rcptr< DiscreteTable<DirType> > xCatPtr;   // only a single RV

  // The default set of FactorOperators in use. Boilerplate.
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

  // derived class members
  prlite::RowMatrix<double> vols;
  //=============================== DEPRECATED ================================

}; // ConditionalPolya

// ***************************************************************************
// ******************************** Family 0 *********************************
// ******************************** Sample ***********************************
// ***************************************************************************

/**
 * This operator samples from a ConditionalPolya as follows: 1) the
 * conditioned-on (z) variable (controlling the choice of Dirichlet,
 * 2) from that Dirichlet the theta probs controlling the target (y)
 * variable, and 3) using these probabilities the value of the target
 * (y) variable. Note that in this we ignore the external p(y)
 * distribution and also whether y is observed or not.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_Sampler: public SampleOperator< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  ///
  emdw::RVVals sample(
    const ConditionalPolya<CondType, DirType>* lhsPtr);
}; // ConditionalPolya_Sampler

/*****************************************************************************/
/********************************** Family 1 *********************************/
/********************************** Normalize ********************************/
/*****************************************************************************/

/// This operator normalizes the internal $p(y)$ and $p(z)$
// distributions. However, since we do not provide for explicit
// normalisation factors in the Dirichlet distributions, the total
// joint also remains unnormalised. Consult the devguide Chapter:
// ConditionalPolya Section: Normalisation for the particular volumes
// applicable.
template <typename CondType, typename DirType>
class ConditionalPolya_InplaceNormalize: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(ConditionalPolya<CondType, DirType>* lstPtr);
}; // ConditionalPolya_InplaceNormalize

/// This operator returns a normalized version of the input factor.
template <typename CondType, typename DirType>
class ConditionalPolya_Normalize: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const ConditionalPolya<CondType, DirType>* lstPtr);
}; // ConditionalPolya_Normalize

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Absorb ***********************************/
/*****************************************************************************/

/// Inplace version of the Absorb operator, see that for more detail.
template <typename CondType, typename DirType>
class ConditionalPolya_InplaceAbsorb: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(ConditionalPolya<CondType, DirType>* lhsPtr,
                      const Factor* rhsFPtr);
}; // ConditionalPolya_InplaceAbsorb

/**
 * The ConditionalPolya can get multiplied with 1) a categorical
 * distribution p(z), where z is the conditioned-on variable, 2) a
 * DirichletSet distribution describing the topic
 * distributions, or 3) a Categorical distribution p(y) where y is the
 * items/words that makes up a topic. Under the hood this operator
 * simply records these distributions so that they can be used later.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_Absorb: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
private:   // not available due to conditional requiring side effects
  Factor* process(const ConditionalPolya<CondType, DirType>* lhsPtr,
                  const Factor* rhsFPtr){
    PRLITE_ASSERT(
      false,
      "ConditionalPolya_Absorb::process NIY since conditional requires side effects!");
    return (ConditionalPolya<CondType,DirType>*) nullptr;
  } // process

}; // ConditionalPolya_Absorb

/*****************************************************************************/
/********************************** Family 2 *********************************/
/********************************** Cancel ***********************************/
/*****************************************************************************/

/// Inplace version of the Cancel operator, see that for more detail.
template <typename CondType, typename DirType>
class ConditionalPolya_InplaceCancel: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  void inplaceProcess(ConditionalPolya<CondType, DirType>* lhsPtr,
                      const Factor* rhsFPtr);
}; // ConditionalPolya_InplaceCancel

/**
 * This operator removes/cancels/divide-by a factor that was absorbed
 * at an earlier stage. See the absorb operator for more detail.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_Cancel: public Operator1<
  ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
private:   // not available due to conditional requiring side effects
  Factor* process(const ConditionalPolya<CondType,
                  DirType>* lhsPtr, const Factor* rhsFPtr){
    PRLITE_ASSERT(
      false,
      "ConditionalPolya_Cancel::process NIY since conditional requires side effects!");
    return (ConditionalPolya<CondType,DirType>*) nullptr;
  } // process

};
// ConditionalPolya_Cancel

/*****************************************************************************/
/********************************** Family 3 *********************************/
/********************************** Marginalize ******************************/
/*****************************************************************************/

/**
 * This operator returns a marginalized factor containing one of: the
 * DirichletSet, the conditioned-on (z) DiscreteTable, or the
 * target (y) DiscreteTable.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_Marginalize: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const ConditionalPolya<CondType, DirType>* lhsPtr,
                  const emdw::RVIds& variablesToKeep,
                  bool presorted = false,
                  const Factor* peekAheadPtr = nullptr);
}; // ConditionalPolya_Marginalize

/*****************************************************************************/
/********************************** Family 4 *********************************/
/********************************** ObserveAndReduce *************************/
/*****************************************************************************/

/**
 * We will only support observing the target (y) value.
 *
 * This version will ONLY yield the correct result if you apply it to
 * all of the factors containing the observed variables.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_ObserveAndReduce: public Operator1< ConditionalPolya<CondType, DirType> > {
public:
  const std::string& isA() const;
  Factor* process(const ConditionalPolya<CondType, DirType>* lhsPtr,
                const emdw::RVIds& variables,
                const emdw::RVVals& assignedVals,
                bool presorted);
}; // ConditionalPolya_ObserveAndReduce

// ***************************************************************************
// ******************************** Family 5 *********************************
// ******************************** Damping **********************************
// ***************************************************************************

/**
 * This operator is NIY, since we dont quite see it being used.
 */
template <typename CondType, typename DirType>
class ConditionalPolya_InplaceDamping: public Operator1< ConditionalPolya<CondType, DirType> > {
public:

  const std::string& isA() const{
    static const std::string CLASSNAME(
      "ConditionalPolya_InplaceDamping<"+std::string(typeid(CondType).name())+std::string(typeid(DirType).name())+std::string(">"));
    return CLASSNAME;
  } // isA

  double inplaceProcess(ConditionalPolya<CondType, DirType>* lhsPtr,
                        const Factor* rhsFPtr,
                        double df) {
    PRLITE_ASSERT(false, "ConditionalPolya_InplaceDamping::inplaceProcess NIY since it was not expected to ever be used!");
    return 0.0;
  } // inplaceProcess

}; // ConditionalPolya_InplaceDamping

#include "conditionalpolya.tcc"

#endif // CONDITIONALPOLYA_HPP
