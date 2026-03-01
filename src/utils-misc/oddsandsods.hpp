#ifndef ODDSANDSODS_HPP
#define ODDSANDSODS_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez (Copyright University of Stellenbosch, all rights reserved.)
          DATE:    26/07/96
          PURPOSE: Some utility classes / functions

*******************************************************************************/

// standard headers
#include <iostream> // istream, ostream
#include <utility>  //pair
#include <math.h>

bool almostEqual(double a, double b, unsigned decimalPlaces);

template <typename T>
bool iseven(const T& val) {
  return !(val %2);
  //return ( val == 2*floor(val/2) );
} // iseven

/**
 * A wrapper around an arbitrary value, but with a check as to whether
 * it is initialised or not.
 */
template <typename T>
class ValidValue {
public:

  /// default ctor : value not valid
  ValidValue()
      : valid(false), value() {}

  /// cs ctor: set valid value
  ValidValue(const T& theValue)
      : valid(true), value(theValue) {}

  ///
  bool operator==(const ValidValue<T>& p) const {
    if (valid != p.valid) {return false;} // if
    if (valid and value != p.value) {return false;} // if
    return true;
  } // operator==

  ///
  bool operator!=(const ValidValue<T>& p) const {
    return !operator==(p);
  } // operator!=

  /// default version should be ok for assigning from another
  ValidValue& operator=(const ValidValue<T>& d) = default;

  /// Assign a value
  ValidValue& operator=(const T& theValue) {
    valid = true;
    value = theValue;
    return *this;
  } // operator=

  void clear() {valid = false;} // clear

public:
  bool valid;
  T value;
};

/**
 * @param a value to be multiplied: c = a*b
 *
 * @param b value to be multiplied: c = a*b
 *
 */
double prodUnderflowProtected(double a, double b);

/**
 * Pulling a very thin skin over the built-in char, just to get past
 * the mess-up in the language standard with std::vector<bool>. But
 * careful, in some senses its more like a char than a bool. For
 * instance you can increment ++ this one where that won't work on a
 * bool.
 */
class Bool {
public:
  inline Bool(const char& theYesNo = 0) : yesNo(theYesNo) {}
  inline Bool (const Bool& theYesNo) : yesNo(theYesNo.yesNo) {}
  inline Bool& operator=(const char& theYesNo) {yesNo = theYesNo; return *this;}
  inline Bool& operator=(Bool theYesNo) {yesNo = theYesNo.yesNo; return *this;}
  inline operator char() const {return yesNo;}
  explicit inline operator char&() {return yesNo;}
  inline Bool& operator++() {yesNo++; return *this;}    // prefix
  inline Bool& operator++(int) {yesNo++; return *this;} // postfix
  inline Bool& operator--() {yesNo--; return *this;}    // prefix
  inline Bool& operator--(int) {yesNo--; return *this;} // postfix
  friend std::ostream& operator<<( std::ostream& file,
                                   Bool n );
  friend std::istream& operator>>( std::istream& file,
                                   Bool& n );

private:
  char yesNo;
};

#endif // ODDSANDSODS_HPP
