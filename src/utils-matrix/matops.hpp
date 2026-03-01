#ifndef MATOPS_HPP
#define MATOPS_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    December 2013
          PURPOSE: some useful matrix routines
          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

#include <vector>
#include "prlite_matrix.hpp"

/******************************************************************************/

/**
 * returns by reference m -= v*v' with v a column vector
 *
 * Lower triangle only, use copyLowerToUpper to get full matrix.
 */
template<typename TOut, typename TInVec>
void
subtractOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v);

/**
 * returns by reference m += v*v' with v a column vector
 *
 * Lower triangle only, use copyLowerToUpper to get full matrix.
 */
template<typename TOut, typename TInVec>
void
addOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v);

/**
 * NOTE Symmetrical: only use this if w is a scaled version of v.
 *
 * returns v*w' with v,w column vectors
 *
 * Lower triangle only, upper is undetermined. Use copyLowerToUpper to
 * get full matrix.
 */
template<typename TOut, typename TInVec>
prlite::RowMatrix<TOut>
outerProdSym(
  const TInVec& v,
  const TInVec& w);

/**
 * NOTE Symmetrical: only use this if w is a scaled version of v and m
 * is symmetrical.
 *
 * returns by reference m += v*w' with v,w column vectors
 *
 * Lower triangle only, use copyLowerToUpper to get full matrix.
 */
template<typename TOut, typename TInVec>
void
addOuterProdSym(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v,
  const TInVec& w);

/**
 * returns by reference m += v*w' with v,w column vectors
 */
template<typename TOut, typename TInVec>
void
addOuterProd(
  prlite::RowMatrix<TOut>& m,
  const TInVec& v,
  const TInVec& w);

/**
 * Packs the lower triangle of a square matrix <mat> into a vector,
 * starting from the main diagonal and then working its way down
 * through the lower diagonals.
 */
template<typename T>
prlite::ColVector<T> packLowerDiagsToVec(const prlite::RowMatrix<T>& mat);

/**
 * Packs the content of a vector <vec> into the lower triangle of a
 * square matrix, starting from the main diagonal and then working its
 * way down through the lower diagonals. The upper triangle of the
 * matrix remains untouched.
 */
template<typename T>
prlite::RowMatrix<T> packVecToLowerDiags(const prlite::ColVector<T>& vec, unsigned dim);

///This function returns a subset of a Matrix consisting of the rows and vectors specified
template<typename T>
prlite::RowMatrix<T> subsetMat(prlite::RowMatrix<T> mat, std::vector<int> rows, std::vector<int> cols);

///This function returns a subset of a RowVector
template<typename T>
prlite::ColVector<T> subsetVec(prlite::ColVector<T> vec, std::vector<int> rows);

/**
 * @brief Does an inplace rank-one Update of a Cholesky Lower
 * triangular matrix L (root of L*L^T).
 *
 * @param L On input the original lower Cholesky matrix. Note that it
 * is column dense to support the fortran convention required by the
 * underlying BLAS routines. On output this will be L' so that L'*L'^T
 * = L*L^T + x*x^T.
 *
 * @param x the input vector to update L with.
 *
 * NOTE: if instead you require L*L^T + c*x*x^T with c a constant scalar,
 * simply do the update with x' = sqrt(c)*x.
 */
void
cholUL(
  prlite::ColMatrix<double>& L,
  const prlite::ColVector<double>& x);


/**
 * @brief Does an inplace rank-one Downdate of a Cholesky Lower
 * triangular matrix L (root of L*L^T).
 *
 * @param L On input the original lower Cholesky matrix. Note that it
 * is column dense to support the fortran convention required by the
 * underlying BLAS routines. On output this will be L' so that L'*L'^T
 * = L*L^T - x*x^T.
 *
 * @param x the input vector to downdate L with.
 *
 * @return success
 *
 * NOTE: if instead you require L*L^T - c*x*x^T with c a constant scalar,
 * simply do the downdate with x' = sqrt(c)*x.
 */
bool
cholDL(
  prlite::ColMatrix<double>& U,
  const prlite::ColVector<double>& x);

/**
 * @brief Does an inplace rank-one Update of an upper/Right triangular
 * matrix R (the root of RR^T). NOTE: We include this unusual upper
 * times lower form as a special factorisation of covariance matrices
 * which, after being inverted, will give the normal lower times upper
 * factorisation for the precision matrix.
 *
 * @param R On input the original LHS UPPER Cholesky matrix. On output
 * this will be R' so that R'*R'^T = R*R^T + x*x^T.
 *
 * @param x the input vector to update R with.
 *
 * NOTE: if instead you require R*R^T + c*x*x^T with c a constant scalar,
 * simply do the update with x' = sqrt(c)*x.
 */
void
cholUR(
  prlite::ColMatrix<double>& R,
  const prlite::ColVector<double>& x);

/**
 * @brief Does an inplace rank-one Downdate of the Cholesky
 * upper/Right triangular matrix R (the root of RR^T). NOTE: We
 * include this unusual upper times lower form as a special
 * factorisation of covariance matrices which, after being inverted,
 * will give the normal lower times upper factorisation for the
 * precision matrix.
 *
 * @param R On input the original LHS UPPER Cholesky matrix. Note that it
 * is column dense to support the fortran convention required by the
 * underlying BLAS routines. On output this will be R' so that R'*R'^T
 * = R*R^T - x*x^T.
 *
 * @param x the input vector to downdate R with.
 *
 * @return success
 *
 * NOTE: if instead you require R*R^T - c*x*x^T with c a constant scalar,
 * simply do the downdate with x' = sqrt(c)*x.
 */
bool
cholDR(
  prlite::ColMatrix<double>& R,
  const prlite::ColVector<double>& x);

/**
 * Rotate the columns of L until it is lower triangular by
 * post-multiplying (i.e. from the right) with Givens matrices. Fairly
 * expensive, (D^3)/4 ops. NOTE: Not to be confused with factorising L.
 */
void rotateColsToLowerTriangular(prlite::ColMatrix<double>& L);

void zeroUpper(prlite::ColMatrix<REAL>& mat);

prlite::RowMatrix<double> copyLowerToLowerAndZeroUpper(const prlite::ColMatrix<double>& lTri);
prlite::ColMatrix<double> copyLowerToUpperAndZeroLower(const prlite::ColMatrix<double>& lTri);

#include "matops.tcc"

#endif // MATOPS_HPP
