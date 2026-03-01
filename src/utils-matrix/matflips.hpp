#ifndef MATFLIPS_HPP
#define MATFLIPS_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    March 2013
          PURPOSE: routines to flip matrices and vectors in various ways
          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

#include "prlite_matrix.hpp"

/******************************************************************************/

/**
 * This flips a vector head to tail.
 */
template <typename AssignType>
prlite::RowVector<AssignType>
flip(
  const prlite::RowVector<AssignType>& inVec);

/**
 * This physically rotates a matrix clockwise through 90 degrees.
 */
template <typename AssignType>
prlite::RowMatrix<AssignType>
rot90(
  const prlite::RowMatrix<AssignType>& inMat);

/**
 * This flips a matrix left-to-right i.e. reverses the column order.
 */
template <typename AssignType>
prlite::RowMatrix<AssignType>
fliplr(
  const prlite::RowMatrix<AssignType>& inMat);

/**
 * This flips a matrix up-to-down i.e. reverses the row order.
 */
template <typename AssignType>
prlite::RowMatrix<AssignType>
flipud(
  const prlite::RowMatrix<AssignType>& inMat);

/**
 * This does a combined lr ud flip i.e. reverses both row and column
 * orders.
 */
template <typename AssignType>
prlite::RowMatrix<AssignType>
fliplrud(
  const prlite::RowMatrix<AssignType>& inMat);

#include "matflips.tcc"

#endif // MATFLIPS_HPP
