#ifndef BRONKERBOSCH_HPP
#define BRONKERBOSCH_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          COPYRIGHT: University of Stellenbosch, all rights reserved.
          DATE:    September 2015
          PURPOSE: The Bron Kerbosch algorithm for finding maximal cliques
          SOURCE: Wiki page at https://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm

*******************************************************************************/

// emdw headers
#include "emdw.hpp"

// standard headers
#include <vector>
#include <set>
#include <map>

/**
 * @brief Find the maximal cliques in an MRF
 *
 * @param neighbours Lists all the vertices adjacent to each vertex
 * @param r These are known to be in the maximal clique (initially
 * empty)
 *
 * @param p Candidates to include in the maximal clique (initially all
 * vertices)
 *
 * @param x These are known to not be in the maximal clique (initially
 * empty)
 *
 * @param maxCliques The maximal cliques found so far, initially
 * empty, on termination contains all maximal cliques arranged from
 * biggest to smallest.
 */
void
BronKerbosch2(
  const std::map<emdw::RVIdType, std::set<emdw::RVIdType> >& neighbours,
  const std::set<emdw::RVIdType>& r,
  std::set<emdw::RVIdType>& p,
  std::set<emdw::RVIdType>& x,
  std::vector< std::set<emdw::RVIdType> >& maxCliques);

// Still to do: The BronKerbosch3 outerloop can enhance efficiency further. See Wiki page

#include "factor.hpp"

/// detects cliques additional to those implicitly implied by the
/// existing factors. Then adds additional vacuous factors to represent
/// those cliques.
void
extendWithVacuousFactors(
  std::vector< rcptr<Factor> >& theFactorPtrs);

#endif // BRONKERBOSCH_HPP
