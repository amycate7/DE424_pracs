#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    October 2013
          PURPOSE: Finding combinations and permutations of items

*******************************************************************************/

// standard headers
#include <iostream> // istream, ostream
#include <vector>
#include <set>
#include <map>

//*********************************************************************
//*************************** Tuples **********************************
//*********************************************************************

/**
 * @brief initTuples and nextTuple jointly generate all combinations
 * of a number of variables where each variable can vary over an
 * independent range -- initTuples generates the first of these
 * sequences after which nextTuple will generate the subsequent
 * sequences.
 *
 * @param sizes: Specifies the range over which each variable can
 * vary. With sizes[i] == Ki, the i'th variable can vary over the
 * range 0..Ki-1. The number of variables in each sequence therefore
 * is sizes.size(). NOTE: <typename T> must be an enumerable type that
 * can vary over range 0..Ki-1.
 */
template<typename T>
std::vector<T>
initTuples(
  const std::vector<size_t>& sizes);

/**
 * @brief initTuples and nextTuple jointly generate all combinations
 * of a number of variables where each variable can vary over an
 * independent range - initTuples generates the first of these
 * sequences after which nextTuple will generate the subsequent
 * sequences. The sequence is ordered such that each next vector will
 * be smaller than the previous one (comparing like indexed values
 * from the left).
 *
 * @param prevVals: On entry the previous one in this sequence of
 * variables -- this determines what the next sequence will be. On
 * exit returns by reference the next combination in this sequence.
 *
 * @return success, i.e true if the resulting prevVal is a valid next
 * combination, false if all possible combinations were exhausted.
 */
template<typename T>
bool
nextTuple(
  const std::vector<size_t>& sizes,
  std::vector<T>& prevVals);

/**
 * @return all possible N^M combinations of M unsigned variables, each
 * varying over a consecutive range. NOTE: the 'all' in the name
 * indicates that there may be repeated values in any particular
 * combination.
 *
 * @param M the number of variables
 *
 * @param N the number of values each variable can take on.
 *
 * @param offset the first value in the consecutive range of N values.
 *
 * EXAMPLE: getTuples(2,3) returns all length 2 sequences
 * that can be formed from the values {0,1,2}:
 * 2 2
 * 1 2
 * 0 2
 * 2 1
 * 1 1
 * 0 1
 * 2 0
 * 1 0
 * 0 0
 */
std::vector< std::vector<unsigned> >
getTuples(
  unsigned M,
  unsigned N,
  unsigned offset=0);

/**
 * @return all possible combinations of M variables Xi (i=0..M-1),
 * each with an individually specified domain. NOTE: the 'all' in the
 * name indicates that there may be repeated values in any particular
 * combination.
 *
 * @param assignRange[i] contains the ARBITRARY domain of the i'th
 * variable (i.e. domain(Xi) ).
 */
template<typename T>
std::vector< std::vector<T> >
getTuples(
  const std::vector< std::vector<T> >& assignRange);


//*********************************************************************
//*************************** Permutations ****************************
//*********************************************************************

/**
 * @return all N!/(N-M)! non-repeating permutations of M unsigned
 * variables, each varying over a consecutive range. NOTE: There are
 * no repeated values in any particular permutation.
 *
 * @param M the number of variables
 *
 * @param N the number of values each variable can take on.
 *
 * @param offset the first value in the consecutive range of N values.
 *
 * EXAMPLE: getPermutations(2,3) draws permutations of any 2 from
 * the set {0,1,2} (draw order makes it difference):
 * 1 2
 * 0 2
 * 2 1
 * 0 1
 * 2 0
 * 1 0
 */
std::vector< std::vector<unsigned> >
getPermutations(unsigned M,
                unsigned N,
                unsigned offset=0);

//*********************************************************************
//*************************** Combinations ****************************
//*********************************************************************

/**
 * @return all N!/(M!(N-M)!) (monotonic increasing) combinations of M
 * unsigned variables, each varying over a consecutive range.
 *
 * @param M the number of variables
 *
 * @param N the number of values each variable can take on.
 *
 * @param offset the first value in the consecutive range of N values.
 *
 * EXAMPLE: getCombinations(2,3) draws combinations of any 2
 * from the set {0,1,2} (draw order does not matter):
 * 1 2
 * 0 2
 * 0 1
 */
std::vector< std::vector<unsigned> >
getCombinations(unsigned M,
                unsigned N,
                unsigned offset=0);



/***********************************************************************/
/* these are deprecated in favour of a more standard naming convention */
/***********************************************************************/

/**
 * @brief This function creates the domains for M variables that each
 * can vary over the range [offset:N-1+offset] (a helper function used
 * by the other combination building functions).
 *
 * @param M is the number of variables for which these domains are needed
 *
 * @param N is the number of values that each variable can have
 *
 * @param offset is the first/smallest value that these variables can take on.
 */
std::vector< std::vector<unsigned> >
setAssignRange(
  unsigned M,
  unsigned N,
  unsigned offset=0);

template<typename T>
std::vector<T>
initCombinations(
  const std::vector<size_t>& sizes);

template<typename T>
bool
nextCombination(
  const std::vector<size_t>& sizes,
  std::vector<T>& prevVals);

template<typename T>
std::vector< std::vector<T> >
getAllCombinations(
  const std::vector< std::vector<T> >& assignRange);

std::vector< std::vector<unsigned> >
getAllCombinations(
  unsigned M,
  unsigned N,
  unsigned offset=0);

std::vector< std::vector<unsigned> >
getUniquePermutations(
  unsigned M,
  unsigned N,
  unsigned offset=0);

std::vector< std::vector<unsigned> >
getUniqueCombinations(
  unsigned M,
  unsigned N,
  unsigned offset=0);





#include "combinations.tcc"

#endif // COMBINATIONS_HPP
