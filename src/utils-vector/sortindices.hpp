#ifndef SORTINDICES_HPP
#define SORTINDICES_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    January 2015
          PURPOSE: Sorting via the indices of a container

          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// standard headers
#include <vector>

/**
 * Finds the sequence of indices through which v is sorted.
 *
 * @param v The (const) container for the values to be sorted. v itself is not changed.
 *
 * @param comp A lambda or functor providing the relational operator
 * to base the sorting on. typically std::less for an ascending or
 * std::greater for a descending sort.
 *
 * @return a vector of size_t's giving the indices through which we
 * have a sorted version of v.
 *
 * Example:
 *
 *  vector<unsigned> blah{5,4,1,2};
 *  vector<size_t> sorted = sortIndices( blah, std::greater<unsigned>() );
 *  cout << sorted << endl;
 *
 * should give 0 1 3 2
 *
 * To sort down, instead use std::less<unsigned>()
 */
template <typename T, typename Lambda>
std::vector<size_t>
sortIndices(
  const T& v,
  Lambda&& comp);

/**
 * Finds the sequence of indices through which a subrange of a
 * container is sorted. NOTE: The container itself remains unsorted.
 *
 * @param iBegin An iterator pointing to the first element in the
 * container i.e. the one at index 0. Of course, if you set this to
 * point to the first element in your subrange, your sequence of
 * indices will start at 0.
 *
 * @param iFrom The (const) An iterator pointing to the first element
 * of the subrange to be sorted.
 *
 * @param iPastTo The (const) An iterator pointing to one past the
 * last element in the sorting subrange.
 *
 * @param comp A lambda or functor providing the relational operator
 * to base the sorting on. typically std::less for an ascending or
 * std::greater for a descending sort.
 *
 * @return a vector of size_t's giving the indices through which we
 * have a sorted version of the container range.
 *
 * Example:
 *
 *  vector<unsigned> blah{5,4,1,2};
 *  vector<size_t> sorted = sortIndices( blah, std::greater<unsigned>() );
 *  cout << sorted << endl;
 *
 * should give 0 1 3 2
 *
 * To sort down, instead use std::less<unsigned>()
 */
template <typename T, typename Lambda>
std::vector<size_t>
sortIndices(
  const T& iBegin,
  const T& iFrom,
  const T& iPastTo,
  Lambda&& comp);

#include "sortindices.tcc"

#endif // SORTINDICES_HPP
