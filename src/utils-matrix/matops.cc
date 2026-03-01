#include "matops.hpp"
#include "lapackblas.hpp"
#include "prlite_logging.hpp"

using namespace std;

/*
  Testing data:

  Original covariance C:
  4.000000e+00 2.000000e+00 -2.000000e+00
  2.000000e+00 5.000000e+00 -5.000000e+00
  -2.000000e+00 -5.000000e+00 8.000000e+00

  Lower triangular Cholesky factor L so that C = LL^T:
  2.000000e+00 0.000000e+00 0.000000e+00
  1.000000e+00 2.000000e+00 0.000000e+00
  -1.000000e+00 -2.000000e+00 1.732051e+00

  Update vector x:
  2.000000e+00 1.000000e+00 2.000000e+00

  L' so that L'L'^T = LL^t + xx^T:
  2.828427e+00 0.000000e+00 0.000000e+00
  1.414214e+00 2.000000e+00 0.000000e+00
  7.071068e-01 -2.000000e+00 2.738613e+00

  Directly calculated C' = C + xx^T
  8.000000e+00 4.000000e+00 2.000000e+00
  4.000000e+00 6.000000e+00 -3.000000e+00
  2.000000e+00 -3.000000e+00 1.200000e+01

 TEST2:

   With original covariance C:
   1 0 0
   0 1 0
   0 0 1

   and update vector x^T
   1 1 1

   Directly calculated C' = C + xx^T:
   2 1 1
   1 2 1
   1 1 2

   R' so that R'R'^T = C':
   1.15470   0.40825   0.70711
   0.00000   1.22474   0.70711
   0.00000   0.00000   1.41421


*/

/*
  //Testing code:

  size_t dim = 3;
  gLinear::gMatrix<REAL, gLinear::COL_DENSE> cv(dim,dim);
  cv(0,0) = 4; cv(1,1) = 5; cv(2,2) = 8;
  cv(0,1) = cv(1,0) = 2; cv(0,2) = cv(2,0) = -2; cv(1,2) = cv(2,1) = -5;
  cout << "cv: " << cv << endl;

  ColVector<double> x(3);
  x[0] = 2; x[1] = 1; x[2] = 2;
  cout << "x: " << x << endl;
  cout << "C + xx^T: " << cv + x*x.transpose() << endl;

  gLinear::gMatrix<REAL, gLinear::COL_DENSE> L(dim,dim);
  int fail;
  cholesky_factor(cv, L, fail);
  cout << "L: " << L << endl;

  // Do the update
  cholUL(L,x);
  cout << "L (update): " << L << endl;
  cout << "LL^T: " << L*L.transpose() << endl;

  // and lets remove the x again using a downdate
  cholDL(L,x);
  cout << "L (downdate): " << L << endl;
  cout << "LL^T: " << L*L.transpose() << endl;
*/

/*
 * cholUL and cholDL is based on:
 *
 * P.E. Gill, G.H. Golub, W. Murray and M.A. Saunders, "Methods for
 * Modifying Matrix Factorizations", Mathematics of Computation, Vol
 * 28 nr 126, Apr 1974. Method C3 pg 521. To work with a lower
 * triangular L, I transposed the given partitioned matrices.
 *
 * Matthias Seeger, "Low Rank Updates for the Cholesky Decomposition",
 * April 5, 2008. Sections 2 and 3 have useful information on fixing
 * sign inversions etc.
 */

void
cholUL(
  prlite::ColMatrix<double>& L,
  const prlite::ColVector<double>& x) {

  // x is shallow, assign to make v deep
  // some avoidable overhead in creating the vector
  prlite::ColVector<double> v; v = x;
  size_t nDim = v.size();
  double c,s;
  int incx =1, incy=1;
  int nRots = nDim;

  for (size_t k = 0; k < nDim; k++) {

    // calculate rotation parameters plus first rotation
    drotg_(&L(k,k), &v[k], &c, &s);
    if (L(k,k) < 0.0) {
      L(k,k) = -L(k,k);
      c = -c;
      s = -s;
    } // if
    nRots--;

    // now do the other rotations
    if (nRots) {
      drot_(&nRots, &L(k+1,k), &incx, &v[k+1], &incy, &c, &s);
    } // if

  } // for k

} // cholUL

bool
cholDL(
  prlite::ColMatrix<double>& L,
  const prlite::ColVector<double>& x) {

  // prlite::ColMatrix<double> LL;
  // LL = L*L.transpose();

  // x is shallow, assign to make p deep
  // some avoidable overhead in creating the vector
  prlite::ColVector<double> p; p = x;
  int nDim = p.size();
  char uplo = 'L';
  char transp = 'N';
  char udiag = 'N';
  int incx =1, incy=1;

  dtrsv_(&uplo, &transp, &udiag, &nDim, L.getStartAddr(), &nDim,p.getStartAddr(), &incx);
  double rho2 = 1.0-(p.transpose()*p);

  if (rho2 > 0.0) { //  rho < std::numeric_limits<double>::min() ){

    double rho = sqrt(rho2);

    // some avoidable overhead in creating the vector
    std::vector<double> r(nDim); for (auto& el : r) { el = 0.0; }
    double c,s;
    int nRots = 0;

    for (int k = nDim-1; k >= 0; k--) { // k=-1 -> k = std::numeric_limits<size_t>::max()

      // calculate rotation parameters plus first rotation reducing p[k] to zero
      drotg_(&rho, &p[k], &c, &s);
      if (rho < 0.0) {
        rho = -rho;
        c = -c;
        s = -s;
      } // if
      nRots++;

      // now do the other rotations to update L and r
      drot_(&nRots, &r[k], &incx, &L(k,k), &incy, &c, &s);

    } // for k

    // L(k,k) might be negative . In that case we have to flip the sign
    // of column k in L. we can speed this up using iterators
    for (int k = 0; k < nDim; k++) {
      if (L(k,k) < 0) {
        for (int m = k; m < nDim; m++) {
          L(m,k) = -L(m,k);
        } // for
      } // if
    } // for

//  return std::abs(1.0-rho) <= 10.0*std::numeric_limits<double>::epsilon()
//         "cholDL:: final rho differs from 1.0 by " << 1.0-rho);
    return true;
  } // if

  else {
    PRLITE_LOGDEBUG2(PRLITE_pLog, "CholDL", "Result of downdate is null or negative definite; rho2: " << rho2);
    return false;
  } // else

} // cholDL

void
cholUR(
  prlite::ColMatrix<double>& R,
  const prlite::ColVector<double>& x) {

  // x is shallow, assign to make v deep
  // some avoidable overhead in creating the vector
  prlite::ColVector<double> v; v = x;
  size_t nDim = v.size();
  double c,s;
  int incx =1, incy=1;
  int nRots = nDim;

  size_t k = nDim;
  while (k--) {

    // calculate rotation parameters plus first rotation
    drotg_(&R(k,k), &v[k], &c, &s);
    if (R(k,k) < 0.0) {
      R(k,k) = -R(k,k);
      c = -c;
      s = -s;
    } // if
    nRots--;

    // now do the other rotations
    drot_(&nRots, &R(0,k), &incx, &v[0], &incy, &c, &s);

  } // for k

} // cholUR

bool
cholDR(
  prlite::ColMatrix<double>& R,
  const prlite::ColVector<double>& x) {
  // x is shallow, assign to make p deep
  // some avoidable overhead in creating the vector
  prlite::ColVector<double> p; p = x;
  int nDim = p.size();
  char uplo = 'U';
  char transp = 'N';
  char udiag = 'N';
  int incx =1, incy=1;

  dtrsv_(&uplo, &transp, &udiag, &nDim, R.getStartAddr(), &nDim, p.getStartAddr(), &incx);
  double rho2 = 1.0-(p.transpose()*p);

  if (rho2 > 0.0) { //  rho < std::numeric_limits<double>::min() ){

    double rho = sqrt(rho2);

    // some avoidable overhead in creating the vector
    std::vector<double> r(nDim); for (auto& el : r) { el = 0.0; }
    double c,s;
    int nRots = 0;

    for (int k = 0; k < nDim; k++) {

      // calculate rotation parameters plus first rotation reducing p[k] to zero
      drotg_(&rho, &p[k], &c, &s);
      if (rho < 0.0) {
        rho = -rho;
        c = -c;
        s = -s;
      } // if
      nRots++;

      // now do the other rotations to update R and r
      drot_(&nRots, &r[0], &incx, &R(0,k), &incy, &c, &s);

    } // for k

    // R(k,k) might be negative . In that case we have to flip the sign
    // of column k in R. we can speed this up using iterators
    for (int k = 0; k < nDim; k++) {
      if (R(k,k) < 0) {
        for (int m = k; m < nDim; m++) {
          R(m,k) = -R(m,k);
        } // for
      } // if
    } // for

    //  return std::abs(1.0-rho) <= 10.0*std::numeric_limits<double>::epsilon()
    //         "cholDR:: final rho differs from 1.0 by " << 1.0-rho);
    return true;
  } // if

  else {
    PRLITE_LOGDEBUG2(PRLITE_pLog, "CholDR", "Result of downdate is null or negative definite; rho2: " << rho2);
    return false;
  } // else

} // cholDR


void rotateColsToLowerTriangular(prlite::ColMatrix<double>& L) {

  size_t nDim = L.rows();
  double c,s;
  int incx =1, incy=1;
  int nRots = nDim;

  for (size_t k = 0; k < nDim; k++) {

    // if L(k,k) < 0.0, invert the sign of COLUMN k
    if (L(k,k) < 0.0) {
      for (size_t rr = k; rr < nDim; rr++) {
        L(rr,k) = -L(rr,k);
      } // for
    } // if

    nRots--;
    for (size_t c2z = k+1; c2z < nDim; c2z++) {

      if (L(k,c2z) > 0.0 or L(k,c2z) < 0.0) {
        // calculate rotation parameters plus first rotation
        drotg_(&L(k,k), &L(k,c2z), &c, &s);
        if (L(k,k) < 0.0) {
          L(k,k) = -L(k,k);
          c = -c;
          s = -s;
        } // if
        L(k,c2z) = 0.0;

        // now do the other rotations in the column
        if (nRots) {
          drot_(&nRots, &L(k+1,k), &incx, &L(k+1,c2z), &incy, &c, &s);
        } // if
      } // if

    } // for c2z

  } // for k

} //rotateColsToLowerTriangular

// zero upper
void zeroUpper(prlite::ColMatrix<REAL>& mat) {
  int nDims = mat.rows();
   for (int c = 1; c < nDims; c++) {
     for (int r = 0; r < c; r++) {
       //std::cout << c << "," << r << " <= " << r << "," << c << std::endl;
       mat(r,c) = 0.0;
     } // for
   } // for
} // zeroUpper

// copy lower to upper
void lower2upper(prlite::ColMatrix<REAL>& mat) {
  int nDims = mat.rows();
   for (int c = 1; c < nDims; c++) {
     for (int r = 0; r < c; r++) {
       //std::cout << c << "," << r << " <= " << r << "," << c << std::endl;
       mat(r,c) = mat(c,r);
     } // for
   } // for
} // lower2upper

// copy upper to lower
void upper2lower(prlite::ColMatrix<REAL>& mat) {
  int nDims = mat.rows();
   for (int c = 1; c < nDims; c++) {
     for (int r = 0; r < c; r++) {
       //std::cout << c << "," << r << " <= " << r << "," << c << std::endl;
       mat(c,r) = mat(r,c);
     } // for
   } // for
} // upper2lower

prlite::RowMatrix<double> copyLowerToLowerAndZeroUpper(const prlite::ColMatrix<double>& lTri) {
  prlite::RowMatrix<double> answ( lTri.rows(), lTri.cols() );
  // all except diagonals
  for (int r = 0; r < lTri.rows(); r++) {
    for (int c = 0; c < r ; c++) {
      answ(c,r) = 0.0;
      answ(r,c) = lTri(r,c);
    } // for c
  } // for r

  // diagonals
  for (int r = 0; r < lTri.rows(); r++) {
    answ(r,r) = lTri(r,r);
  } // for r

  return answ;
} // copyLowerToLowerAndZeroUpper

prlite::ColMatrix<double> copyLowerToUpperAndZeroLower(const prlite::ColMatrix<double>& lTri) {
  prlite::ColMatrix<double> answ( lTri.rows(), lTri.cols() );
  // all except diagonals
  for (int r = 0; r < lTri.rows(); r++) {
    for (int c = 0; c < r ; c++) {
      answ(r,c) = 0.0;
      answ(c,r) = lTri(r,c);
    } // for c
  } // for r

  // diagonals
  for (int r = 0; r < lTri.rows(); r++) {
    answ(r,r) = lTri(r,r);
  } // for r

  return answ;
} // copyLowerToUpperAndZeroLower
