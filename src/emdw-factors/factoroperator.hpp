#ifndef FACTOROPERATOR_HPP
#define FACTOROPERATOR_HPP

/*******************************************************************************
 *                                                                             *
 * AUTHOR:  JA du Preez                                                        *
 * Created on : /02/2013                                                       *
 *                                                                             *
 * PURPOSE: The abstract base class for all FactorOperator variants            *
 * COPYRIGHT: University of Stellenbosch, all rights reserved.                 *
 *                                                                             *
 *******************************************************************************/

// standard headers
#include <string>  // string
#include <iosfwd>  // istream, ostream

#include "emdw.hpp"

/**
 * Factor functors/operations:

 * There are a number of factor operations that all factors will
 * have to support:
 *
 * 1) Absorb/Incorporate/product of factors:
 *
 * Factor1Ptr->absorb(Factor2Ptr);
 * Factor3Ptr = Factor1Ptr->absorb(Factor2Ptr);
 *
 * Can in some cases lead to new factor types. For instance, the
 * product of two gaussians is no longer a gaussian, but a
 * generalised gaussian. Can in certain cases be done
 * inplace. Log-form: sum.
 *
 * 2) Removal/division of factors: Assuming that we will only apply
 * this to factors arising from a product, this should not introduce
 * new factors above those already associated with the original
 * product. Log-form: subtract.
 *
 * Factor1Ptr->cancel(Factor2Ptr);
 * Factor3Ptr = Factor1Ptr->cancel(Factor2Ptr);
 *
 * 3) Marginalization: Same factor type, smaller scope. Usually not
 * done inplace.
 *
 * Factor3Ptr = Factor1Ptr->marginalize(variablesToKeep, presorted);
 *
 * 4) Evidence observation: Same factor type, some probs = 0. Can be
 * done inplace. Observation itself is a std::vector<Whatever>.
 *
 * Factor1Ptr->observe(variables, values, presorted);
 * Factor3Ptr = Factor1Ptr->observe(variables, values, presorted);
 *
 *
 * 5) Normalization: Same factor type, rescaled probs. Can be done
 * inplace.
 *
 * Factor1Ptr->normalize();
 * Factor3Ptr = Factor1Ptr->normalize();
 *
 * 6) Damping: Weighted sum of messages of same type and scope. Same
 * factor type, same scope. Can be done inplace.
 *
 * Decision theory introduces a few others, we will get to them at a
 * later time.
 *
 * In some cases we can apply these operations inplace i.e. the
 * operation only changes the existing factor somewhat, but creates
 * no new factor. Mostly we will use them as a functions which
 * create new factors without changing the operand factors. These
 * new factors might be of a different kind than those it was
 * derived from; one of the challenges will be to keep the system
 * closed in terms of the types of factors that can be generated in
 * this manner.

 * The implementation details of these operations are particular to
 * the factors involved (we'll probably templatize on the
 * factors). In order to support the implementation of its own
 * functors the factors will have to provide the necessary access to
 * the functors, for instance functors operating on discrete tables
 * will have to be able to iterate over those tables. Will probably
 * be good to supply iterators for all table-based factors.
 *
 * Consideration: Do we embed the functors into the factor, or do
 * they have a own life outside the factors?
 */


/// This class collects some member functions common to all factor
/// operators. It contains no data; there is no diamond problem.
class FactorOperator{
public:
  virtual ~FactorOperator() {}
  virtual const std::string& isA() const = 0;

  /// output
  friend std::ostream& operator<<( std::ostream& file,
                                   const FactorOperator& n );


private:
  // this is just to convince the compiler that this is a base class.
  virtual void dummy() {}
}; // FactorOperator

// These classes specialize FactorOperator to distinct types in order
// to provide some type safety to constructors etc. However, this
// protection is currently only afforded by the type of pointer you
// allocate an Operator1 class to - by itself Operator 1 at this stage
// is still able to freely associate with any of these. At some point
// in the future these probably will become template parameters common
// to both Operator1 and FactorOperator, thereby truly making of them
// separate types.
class Sampler : public virtual FactorOperator{}; // Family 0
class Normalizer : public virtual FactorOperator{}; // Family 1
class InplaceNormalizer : public virtual FactorOperator{}; // Family 1
class Absorber : public virtual FactorOperator{}; // Family 2
class InplaceAbsorber : public virtual FactorOperator{}; // Family 2
class Canceller : public virtual FactorOperator{}; // Family 2
class InplaceCanceller : public virtual FactorOperator{}; // Family 2
class Marginalizer : public virtual FactorOperator{}; // Family 3
class ObserveAndReducer : public virtual FactorOperator{}; // Family 4
class InplaceDamper : public virtual FactorOperator{}; // Family 5


// forward declaration
class Factor;

// TO BE DEPRECATED sometime in the future. Mid-level abstract factor
// operator processor. This is the old / non-type safe generic factor
// operator, retained here for legacy reasons. It seems that before we
// can scrap this one, we will have to basically rid the whole system
// of the old FactorOperator versions. The problem seems to be tied up
// with the raw pointer code in the Factor classes such as
// normalize(FactorOperator* procPtr). These gets called from a rcptr
// version in Factor. Unfortunately that rcptr version can not be
// replaced with one using the detailed new classes below, since
// things become ambiguous.
template<typename FactorType>
class Operator1:
  public virtual Sampler,
  public virtual Normalizer,
  public virtual InplaceNormalizer,
  public virtual Absorber,
  public virtual InplaceAbsorber,
  public virtual Canceller,
  public virtual InplaceCanceller,
  public virtual Marginalizer,
  public virtual ObserveAndReducer,
  public virtual InplaceDamper{

public:

  //virtual const std::string& isA() const = 0;

  //------------------Family 0
  virtual emdw::RVVals sample(const FactorType* lhsPtr);

  //------------------Family 1
  virtual Factor* process(const FactorType* lhsPtr);
  virtual void inplaceProcess(FactorType* lhsPtr);

  //------------------Family 2
  virtual Factor* process(const FactorType* lhsPtr,
                          const Factor* rhsPtr);
  virtual void inplaceProcess(FactorType* lhsPtr,
                              const Factor* rhsPtr);

  //------------------Family 3
  virtual Factor* process(const FactorType* lhsPtr,
                          const emdw::RVIds& variablesToKeep,
                          bool presorted = false,
                          const Factor* peekAheadPtr = nullptr);

  //------------------Family 4
  virtual Factor* process(const FactorType* lhsPtr,
                          const emdw::RVIds& variables,
                          const emdw::RVVals& assignedVals,
                          bool presorted = false);

  //------------------Family 5
  virtual double inplaceProcess(FactorType* lhsPtr,
                                const Factor* rhsPtr,
                                double df);

}; // Operator1

// Below follows a series of operator classes that are intended to
// replace Operator1 in the long run.  However, the upgrade probably
// will be quite intrusive requiring all cases of FactorOperator and
// Operator1 to be purged from the system. On ice for the moment.

template<typename FactorType>
class InplaceDampenOperator:
  public virtual InplaceDamper{

public:

  //------------------Family 5
  virtual double inplaceProcess(FactorType* lhsPtr,
                                const Factor* rhsPtr,
                                double df) = 0;

}; // InplaceDampenOperator

template<typename FactorType>
class ObserveAndReduceOperator:
  public virtual ObserveAndReducer{

public:

  //------------------Family 4
  virtual Factor* process(const FactorType* lhsPtr,
                          const emdw::RVIds& variables,
                          const emdw::RVVals& assignedVals,
                          bool presorted = false) = 0;

}; // ObserveAndReduceOperator

template<typename FactorType>
class MarginalizeOperator:
  public virtual Marginalizer{

public:

  //------------------Family 3
  virtual Factor* process(const FactorType* lhsPtr,
                          const emdw::RVIds& variablesToKeep,
                          bool presorted = false,
                          const Factor* peekAheadPtr = nullptr) = 0;

}; // MarginalizeOperator

template<typename FactorType>
class InplaceCancelOperator:
  public virtual InplaceCanceller{

public:

  //------------------Family 2
  virtual void inplaceProcess(FactorType* lhsPtr,
                              const Factor* rhsPtr) = 0;
}; // InplaceCancelOperator

template<typename FactorType>
class CancelOperator:
  public virtual Canceller{

public:

  //------------------Family 2
  virtual Factor* process(const FactorType* lhsPtr,
                          const Factor* rhsPtr) = 0;
}; // CancelOperator

template<typename FactorType>
class InplaceAbsorbOperator:
  public virtual InplaceAbsorber{

public:

  //------------------Family 2
  virtual void inplaceProcess(FactorType* lhsPtr,
                              const Factor* rhsPtr) = 0;
}; // InplaceAbsorbOperator

template<typename FactorType>
class AbsorbOperator:
  public virtual Absorber{

public:

  //------------------Family 2
  virtual Factor* process(const FactorType* lhsPtr,
                          const Factor* rhsPtr) = 0;
}; // AbsorbOperator

template<typename FactorType>
class InplaceNormalizeOperator:
  public virtual InplaceNormalizer{

public:

  //------------------Family 1
  virtual void inplaceProcess(FactorType* lhsPtr) = 0;
}; // InplaceNormalizeOperator

template<typename FactorType>
class NormalizeOperator:
  public virtual Normalizer{

public:

  //------------------Family 1
  virtual Factor* process(const FactorType* lhsPtr) = 0;

}; // NormalizeOperator

template<typename FactorType>
class SampleOperator:
  public virtual Sampler
{

public:

  //------------------Family 0
  virtual emdw::RVVals sample(const FactorType* lhsPtr) = 0;

}; // SampleOperator


#include "factoroperator.tcc"

#endif // FACTOROPERATOR_HPP
