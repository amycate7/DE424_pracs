/*
 * Author     : JA du Preez
 * Created on : /02/2013
 * Copyright  : University of Stellenbosch, all rights retained
 */

// patrec headers
// standard headers
#include <string>  // string
#include <iostream>  // cout, endl
//#include <typeinfo>

#include "factoroperator.hpp"

using namespace std;

//================================== Friends ==================================

/// output
std::ostream& operator<<( std::ostream& file,
                          const FactorOperator& n ){
  file << n.isA();
  return file;
} // operator<<
