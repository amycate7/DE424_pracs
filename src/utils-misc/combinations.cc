/*******************************************************************************

          AUTHOR:    JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:      March 2013
          PURPOSE:   Finding combinations and permutations of items

*******************************************************************************/

#include "combinations.hpp"

using namespace std;

vector< vector<unsigned> >
getTuples(
  unsigned M,
  unsigned N,
  unsigned offset) {

  vector< vector<unsigned> > assignRange( setAssignRange(M,N,offset) );
  return getAllCombinations(assignRange);
} // getAllCombinations

vector< vector<unsigned> >
getPermutations(
  unsigned M,
  unsigned N,
  unsigned offset) {

  vector< vector<unsigned> > uniquePermutations;
  vector<size_t> varSizes(M);
  for (auto& el : varSizes) {el = N;} // for
  std::vector<unsigned> permute = initTuples<unsigned>(varSizes);

  do {

    bool allDiffer = true;
    for (unsigned i = 0; i < M; i++) {
      for (unsigned j = i+1; j < M; j++) {
        if (permute[i] == permute[j]) {
          allDiffer = false;
          break;
        } // if
      } // for j
      if (!allDiffer) {break;} // if
    } // for i

    if (allDiffer) { // we've got a legit sequence here, they all differ
      uniquePermutations.push_back(permute);
      for (auto& el : uniquePermutations[uniquePermutations.size()-1]) {
        el += offset;
      } // for
    } // if

  } while( nextCombination(varSizes, permute) );

  return uniquePermutations;
} // getPermutations

vector< vector<unsigned> >
getCombinations(
  unsigned M,
  unsigned N,
  unsigned offset) {

  vector< vector<unsigned> > uniqueCombinations;
  vector<size_t> varSizes(M);
  for (auto& el : varSizes) {el = N;} // for
  std::vector<unsigned> combine = initTuples<unsigned>(varSizes);

  do {

    bool allBigger = true;
    for (unsigned i = 0; i < M; i++) {
      for (unsigned j = i+1; j < M; j++) {
        if (combine[i] >= combine[j]) {
          allBigger = false;
          break;
        } // if
      } // for j
      if (!allBigger) {break;} // if
    } // for i

    if (allBigger) { // we've got a legit sequence here, they all differ
      uniqueCombinations.push_back(combine);
      for (auto& el : uniqueCombinations[uniqueCombinations.size()-1]) {
        el += offset;
      } // for
    } // if

  } while( nextCombination(varSizes, combine) );

  return uniqueCombinations;
} // getCombinations

/***********************************************************************/
/* these are deprecated in favour of a more standard naming convention */
/***********************************************************************/

vector< vector<unsigned> >
setAssignRange(
  unsigned M,
  unsigned N,
  unsigned offset) {

  vector< vector<unsigned> > assignRange(M);

  for (unsigned i = 0; i < M; i++) {
    assignRange[i].resize(N);
    for (unsigned j = 0; j < N; j++) {
      assignRange[i][j] = unsigned(j+offset);
    } // for j
  } // for i
  return assignRange;
} // setAssignRange


// defunct
vector< vector<unsigned> >
getAllCombinations(unsigned M,
                   unsigned N,
                   unsigned offset) {
  std::cout << __FILE__ << __LINE__ << " ::getAllCombinations is deprecated. Rather use getTuples\n";
  return getTuples(M,N,offset);
} // getAllCombinations

// defunct
vector< vector<unsigned> >
getUniquePermutations(unsigned M,
                      unsigned N,
                      unsigned offset) {

  std::cout << __FILE__ << __LINE__ << " ::getUniquePermutations is deprecated. Rather use getPermutations\n";
  return getPermutations(M,N,offset);
} // getUniquePermutations

// defunct
vector< vector<unsigned> >
getUniqueCombinations(unsigned M,
                      unsigned N,
                      unsigned offset) {

  std::cout << __FILE__ << __LINE__ << " ::getUniqueCombinations is deprecated. Rather use getCombinations\n";
  return getCombinations(M,N,offset);
} // getUniqueCombinations
