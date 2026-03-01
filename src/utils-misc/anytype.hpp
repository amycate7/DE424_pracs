#ifndef ANYTYPE_HPP
#define ANYTYPE_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
          (Similar to boost::any. A re-write of a class previously done by Albert Visagie)
          DATE:    March 2013
          PURPOSE: An abstract container that can save any other concrete type

*******************************************************************************/

#include "prlite_error.hpp"
//#include <boost/any.hpp>

/******************************************************************************/
/****************************** class TypeBase ********************************/
/******************************************************************************/

/**
 * This one is an abstract utility class you will never directly
 * use. It is used to keep the specific template type T out of AnyType -
 * we don't want to templatize that too.
 */
class TypeBase{
  friend class AnyType;
protected:

  //======================== Constructors and Destructor ========================

  TypeBase() {};
  TypeBase(const TypeBase&) {}
  virtual ~TypeBase() {}

  //===================== Required Virtual Member Functions =====================

  virtual TypeBase* copy() = 0;
  virtual bool isEqual(const TypeBase* rhsPtr) = 0;

}; // class TypeBase

/******************************************************************************/
/******************************** class Type<T> *******************************/
/******************************************************************************/

/**
 * And this is the derived class that is used to create all the
 * variants of AnyType.
 */
template <typename T>
class Type: public TypeBase{
  friend class AnyType;
public:
  //======================== Constructors and Destructor ========================

  Type(const T& theVal): val(theVal) {}

  //========================= Operators and Conversions =========================

  inline operator T&() {
    return val;
  }

  inline operator const T&() const {
    return val;
  }

  //===================== Required Virtual Member Functions =====================

  Type<T>* copy() {return new Type<T>(val);}

  bool isEqual(const TypeBase* rhsPtr) {
    const Type<T>* dwnPtr = dynamic_cast<const Type<T>*>(rhsPtr);
    PRLITE_ASSERT(dwnPtr, "This is odd, comparing AnyTypes of different base types");
    return val == dwnPtr->val;
  } // isEqual

  //=============================== Data Members ================================

protected:
  T val;
}; // class Type<T>

/******************************************************************************/
/******************************** class AnyType ***********************************/
/******************************************************************************/

/**
 * This is a container for arbitrary 'Type<T>' variables- you can
 * stick just about anything into an 'AnyType'. You can retrieve it again
 * via an explicit type conversion. This, of course, implies that you
 * need to know what is in there. Example:
 *
 * AnyType x(string("anytype"));
 * string retrieved = string(x);
 *
 * Saving and loading from files is a bit of an issue, generally
 * speaking you would not know what you are loading. One can address
 * this via RTTI or a whole virtual construction thing, but for the
 * moment the recommendation is to save the actual final types and
 * make sure that the loader knows what exactly to expect. We'll see
 * how this pans out in practical use.
 */
class AnyType{
public:

  //======================== Constructors and Destructor ========================
public:

  // class specific ctor allows you to stick any type into here
  template <typename T>
  AnyType(const T& theVal): basePtr(new Type<T>(theVal)) {}

  // default ctor
  AnyType() : basePtr(nullptr) {};

  // copy ctor
  AnyType(const AnyType& rhs) : basePtr( rhs.basePtr->copy() ) {}

  // move ctor
  AnyType(AnyType&& rhs) : basePtr(rhs.basePtr) {
    rhs.basePtr = nullptr;
  } // move ctor

  // dtor
  ~AnyType() {delete basePtr;}

  //========================= Operators and Conversions =========================
public:

  /**
   * Converts to actual underlying type. AnyType mismatch between the actual and
   * requested types will break badly, even seemingly innocent stuff
   * like to double from float.
   */

  template <typename T>
  inline explicit operator T() const {
// #ifdef DEBUG_MODE
    Type<T>* cstPtr = dynamic_cast<Type<T>*>(basePtr);
    PRLITE_ASSERT( cstPtr, "Request to convert AnyType to invalid type " << typeid(T).name() );
    return cstPtr->val;
// #else
//    return static_cast<Type<T>*>(basePtr)->val;
// #endif
  } // operator T()

  // assignment from final type
  template<typename T>
  AnyType& operator=(const T& theVal) {
    delete basePtr;
    basePtr = new Type<T>(theVal);
    return *this;
  } // operator=

  /**
   * Assign from another AnyType
   */
  AnyType& operator=(const AnyType& a) {
    if (a.basePtr != basePtr) {
      delete basePtr;
      basePtr = a.basePtr->copy();
    } // if
    return *this;
  } // operator=

  /**
   * Move assign from another AnyType
   */
  AnyType& operator=(AnyType&& a) {
    basePtr = a.basePtr;
    a.basePtr = nullptr;
    return *this;
  } // operator=

  /**
   * equality. Be very alert with this, some results are
   * unexpected. Example: int(1) != unsigned(1)
   */
  bool operator==(const AnyType& a) const {
    if (basePtr && a.basePtr) {
      return basePtr->isEqual(a.basePtr);
    } // if
    // check for both being zero
    return ! (basePtr || a.basePtr);
  }

  /**
   * inequality. Be very alert with this, some results are
   * unexpected. Example: int(1) != unsigned(1)
   */
  bool operator!=(const AnyType& a) const {
    return !operator==(a);
  }

  //=============================== Data Members ================================

private:
  TypeBase* basePtr;

}; // class AnyType

/******************************************************************************/
/****************************** function extract<T> ***************************/
/******************************************************************************/

/**
 * Extract, with possible type conversion, a subvector from a
 * vector. The vector must support operator[], and the returned type
 * is always std::vector<T>.
 *
 * NOTE: Handy for extracting a homogenous sub-vector from a vector of
 * AnyType, but T must EXACTLY match the actual underlying type in the
 * AnyType.  Even seemingly innocent conversions such as float to double
 * will break. Example code:

 std::vector<AnyType> any{0.0, 1.1};
 any.push_back(string("twee"));
 any.push_back(3.3);
 any.push_back(string("vier"));
 any.push_back(5.5);
 std::vector<double> doubleVec(extract<double>(any, std::vector<int>{0,1,3,5}));
 cout <<"doubleVec: " << doubleVec << endl;
 std::vector<string> stringVec(extract<string>(any, std::vector<int>{2,4}));
 cout <<"stringVec: " << stringVec << endl;

*/
/*
  template<typename T, typename VEC>
  std::vector<T> extract( const VEC& vec,
  const std::vector<int>& idxs ) {

  std::vector<T> ret( idxs.size() );
  for (unsigned i = 0; i < static_cast<unsigned>( idxs.size() ); i++) {
  ret[i] = static_cast<T>(vec[ idxs[i] ]);
  } // for i

  return ret;

  } // extract
*/

//#include "prlite_stlvecs.hpp"

/**
 * This extracts specified elements from an input vector-like
 * container fromVec (supports [] indexing) while also converting them
 * to output type TO. idxs (another vector or set of an ordinal type)
 * specifies which elements to select.
 */
template<typename TO, typename CTI, typename IDX>
std::vector<TO> extract( const CTI& fromVec,
                      const IDX& idxs ) {
  std::vector<TO> ret; ret.reserve( idxs.size() );
  for (auto i : idxs) {
    ret.push_back( static_cast<TO>(fromVec[i]) );
  } // for i

  return ret;
} // extract

/**
 * This extracts specified elements from an input vector<TI> while
 * also converting them to output type TO. idxs (another vector or
 * set of an ordinal type) specifies which elements to select.
 */
template<typename TO, typename TI, typename IDX>
std::vector<TO> defunctExtract( const std::vector<TI>& vec,
                        const IDX& idxs ) {
  std::vector<TO> ret; ret.reserve( idxs.size() );
  for (auto i : idxs) {
    ret.push_back( static_cast<TO>(vec[i]) );
  } // for i

  return ret;
} // defunctExtract

#include "anytype.tcc"

#endif // ANYTYPE_HPP
