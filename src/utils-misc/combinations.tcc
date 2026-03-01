/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    October 2014
          PURPOSE: Finding combinations and permutations of items

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

#include "prlite_stlvecs.hpp"


template<typename T>
std::vector<T>
initTuples(
  const std::vector<size_t>& sizes) {

  std::vector<T> init( sizes.size() );
  auto itr1 = sizes.begin();
  auto itr2 = init.begin();
  while ( itr1 != sizes.end() ) {
    *itr2 = *itr1 -1;
    itr1++; itr2++;
  } // while
  return init;
} // initTuples

// vector<my_class>::reverse_iterator i = my_vector.rbegin();

template<typename T>
bool
nextTuple(
  const std::vector<size_t>& sizes,
  std::vector<T>& prevVals) {

  // find the first non-zero one
  auto itr1 = sizes.rbegin();
  auto itr2 = prevVals.rbegin();
  while (itr1 != sizes.rend() and *itr2 == 0) {
    *itr2 = (*itr1)-1;
    itr1++; itr2++;
  } // while

  if ( itr1 == sizes.rend() ) {return false;}
  (*itr2)--;

  return true;
} // nextTuple

template<typename T>
std::vector< std::vector<T> >
getTuples(const std::vector< std::vector<T> >& assignRange) {

  std::vector< std::vector<T> > allAssignments;
  unsigned M = assignRange.size();

  // get the cardinalities of all the variables
  std::vector<size_t> varSizes(M);
  size_t totSz = 1;
  for (unsigned i = 0; i < M; i++) {
    varSizes[i] = assignRange[i].size();
    totSz *= varSizes[i];
  } // for i
  allAssignments.resize(totSz);

  // and pack all those combinations
  size_t idx = 0;
  std::vector<size_t> indices = initTuples<size_t>(varSizes);
  do {
    allAssignments[idx].resize(M);
    for (unsigned j = 0; j < M; j++) {
      allAssignments[idx][j] = assignRange[j][ indices[j] ];
    } // for
    idx++;
  } while( nextCombination(varSizes, indices) );

  if (idx != totSz) std::cout << "getAllCombinations : something is wrong" << std::endl;

  return allAssignments;
} // getTuples

/***********************************************************************/
/* these are deprecated in favour of a more standard naming convention */
/***********************************************************************/

// defunct
template<typename T>
std::vector<T>
initCombinations(const std::vector<size_t>& sizes) {
  std::cout << __FILE__ << __LINE__ << " ::initCombinations is deprecated. Rather use initTuples\n";
  return initTuples<T>(sizes);
} // initCombinations

template<typename T>
bool
nextCombination(const std::vector<size_t>& sizes, std::vector<T>& prevVals) {

  return nextTuple(sizes, prevVals);
} // nextCombination

template<typename T>
std::vector< std::vector<T> >
getAllCombinations(const std::vector< std::vector<T> >& assignRange) {

  std::cout << __FILE__ << __LINE__ << " ::getAllCombinations is deprecated. Rather use getTuples\n";
  return getTuples(assignRange);
} // getAllCombinations
