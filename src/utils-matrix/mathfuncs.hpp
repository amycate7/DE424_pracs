#ifndef MATHFUNCS_HPP
#define MATHFUNCS_HPP

/*******************************************************************************

          AUTHORS: DJP Badenhorst, JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    Sept 2014
          PURPOSE: Various mathematical and probability-related functions

*******************************************************************************/

// standard headers
#include <iostream> // istream, ostream
#include <vector>
#include <set>
#include <map>

#include "prlite_matrix.hpp"


/*******************************************************************************************/
/*               Various math and prob functions                                           */
/*******************************************************************************************/

double digamma(double input);

double logSumExp(double l1, double l2);

#endif // MATHFUNCS_HPP
