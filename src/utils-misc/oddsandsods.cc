/*
 * Author     : JA du Preez (DSP Group, E&E Eng, US)
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

#include "oddsandsods.hpp"

#include <limits>

double prodUnderflowProtected(double a, double b){
  if (a > 0.0 and b > 0.0) {
    // using min since denorm_min still cause fracturing.
    return std::max( (a*b), std::numeric_limits<double>::min() );
  } // if
  else {return 0.0;}
} // prodUnderflowProtected

bool almostEqual(double a, double b, unsigned decimalPlaces) {
  // check for nan
  if (a == a && b == b) {
    // not nan. check inf.
    static const double inf = std::numeric_limits<double>::infinity();
    static const double minf = -std::numeric_limits<double>::infinity();
    if ( (a == inf || a == minf) || (b == inf || b == minf) ) {
      return a == b;
    }
    else {
      // not nan or inf

      // make sure a is the smallest
      if (fabs(a)>fabs(b)) {
        return fabs(a-b) <= (fabs(b)/std::pow(10.0,static_cast<double>(decimalPlaces)));
      }
      else {
        return fabs(a-b) <= (fabs(a)/std::pow(10.0,static_cast<double>(decimalPlaces)));
      }
    }
  }
  else {
    return (!(a==a)) && (!(b==b));
  }
  // dummy for dumb compilers' warnings
  return false;
} // almostEqual

std::ostream& operator<<( std::ostream& file,
                          Bool n ){
  file << short(n.yesNo);
  return file;
} // operator<<


std::istream& operator>>( std::istream& file,
                          Bool& n ){
  file >> n.yesNo;
  return file;
} // operator>>
