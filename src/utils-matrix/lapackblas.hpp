#ifndef LAPACKBLAS_HPP
#define LAPACKBLAS_HPP

/*******************************************************************************

          AUTHOR:  JA du Preez
          DATE:    August 2017
          PURPOSE: LaPack / BLAS interfacing
          (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// http://www.netlib.org/blas/
// http://www.netlib.org/lapack/


//http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=5&t=3233
// lapack naming convention: https://www.netlib.org/lapack/lug/node24.html

//#################################################################################
//################## GIVENS ROTATIONS #############################################
//#################################################################################

// BLAS routine for calculating parameters for Givens rotation
// see www.mathkeisan.com/usersguide/man/DROTG.html
/*
NAME
       DGROTG - BLAS level one rotation subroutines

SYNOPSIS
       SUBROUTINE DROTG    ( a, b, c, s )

           DOUBLE PRECISION  a, b, c, s


DESCRIPTION
       DROTG  computes  the  elements  of  a Givens plane rotation matrix such
       that:

             _      _     _   _    _   _
             | c  s |     | a |    | r |
             |-s  c | *   | b | =  | 0 |
             -      -     -   -    -   -

       where  r = +- sqrt ( a**2 + b**2 )  and c**2 + s**2 =1.

       The Givens plane rotation can be used to introduce zero elements into
       a matrix selectively.



ARGUMENTS
       a       (input and output) DOUBLE PRECISION

               First vector component.  On input, the first component  of  the
               vector to be rotated.  On output, a is overwritten by by r, the
               first component of the vector in the rotated coordinate  system
               where:

               r = sign(sqrt(a**2 + b**2),a), if |a| > |b|

               r = sign(sqrt(a**2 + b**2),b), if |a| <= |b|

       b       (input and output) DOUBLE PRECISION
               Second vector component.

               On input, the second component of the vector to be rotated.  On
               output, b contains z, where:

                   z=s    if  |a| >  |b|
                   z=1/c  if  |a| <= |b|  and  c != 0  and r != 0
                   z=1    if  |a| <= |b|  and  c  = 0  and r != 0
                   z=0    if  r = 0

       c       (output) DOUBLE PRECISION
               Cosine of the angle of rotation:

                    c = a/r  if  r != 0
                    c = 1    if  r  = 0

       s       (output) DOUBLE PRECISION
               Sine of the angle of rotation:

                    s = b/r  if  r != 0
                    s = 0    if  r  = 0
 */
extern "C" void drotg_(double* da, double* db, double* c, double* s);


// BLAS routine for inplace 2dim rotation for x,y vector
/*
NAME
       DROT - BLAS level one, plane rotation subroutines

SYNOPSIS
       SUBROUTINE DROT    ( n, x, incx, y, incy, c, s )

           INTEGER          n, incx, incy

           DOUBLE PRECISION x, y, c, s


DESCRIPTION
       DROT   applies  a  plane  rotation matrix to a real sequence of ordered
       pairs:

            (x , y ), for all i = 1, 2, ..., n.
              i   i

ARGUMENTS
       n       INTEGER. (input)
               Number of ordered pairs (planar points in DROT) to be  rotated.
               If n <= 0, this routine returns without computation.

       x       DOUBLE PRECISION, (input and output)
               Array  of dimension (n-1) * |incx| + 1.  On input, array x con-
               tains the x-coordinate of each planar point to be rotated.   On
               output,  array x contains the x-coordinate of each rotated pla-
               nar point.

       incx    INTEGER. (input)
               Increment between elements of x.  If incx = 0, the results will
               be unpredictable.

       y       DOUBLE PRECISION, (input and output)
               array of dimension (n-1) * |incy| + 1.
               On  input,  array  y  contains  the y-coordinate of each planar
               point to be rotated.  On output, array y contains the y-coordi-
               nate of each rotated planar point.

       incy    INTEGER. (input)
               Increment between elements of y.  If incy = 0, the results will
               be unpredictable.

       c       DOUBLE PRECISION, Cosine of the angle of rotation.
               (input)

       s       DOUBLE PRECISION, Sine of the angle of rotation. (input)

NOTES
       This routine is Level 1  Basic  Linear  Algebra  Subprograms  (Level  1
       BLAS). It applies the following plane rotation to each pair of elements
       (x , y ):
                               i   i
            _      _       _      _     _      _
            | x(i) |  <--  | c  s |  .  | x(i) |
            | y(i) |       |-s  c |     | y(i) |
            -      -       -      -     -      -

            for i = 1,...,n
                                        2    2
       If coefficients c and s satisfy c  + s  = 1.0, the rotation matrix
       is orthogonal, and the transformation is called a Givens plane
       rotation. If c = 1 and s = 0, DROT returns without modifying any
       input parameters.

       To calculate the Givens coefficients c and s from a two-element
       vector to determine the angle of rotation, use SROTG(3S).

       When working backward (incx < 0 or incy < 0), each routine starts
       at the end of the vector and moves backward, as follows:

            x(1-incx * (n-1)), x(1-incx * (n-2)), ..., x(1)

            y(1-incy * (n-1)), y(1-incy * (n-2)), ..., y(1)

 */
extern "C" void drot_(int* n, double* x, int* incx, double* y, int* incy, double* c, double* s);

//#################################################################################
//################## PRODUCTS #####################################################
//#################################################################################

//*********************************************************************************
//****************** vector . vector **********************************************
//*********************************************************************************

/*
function ddot   (       integer         N,
                double precision, dimension(*)          DX,
                integer         INCX,
                double precision, dimension(*)          DY,
                integer         INCY
Purpose:

        DDOT forms the dot product of two vectors.
        uses unrolled loops for increments equal to one.

Parameters
    [in]    N
              N is INTEGER number of elements in input vector(s)

    [in]    DX
              DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )

    [in]    INCX
              INCX is INTEGER storage spacing between elements of DX

    [in]    DY
              DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )

    [in]    INCY
              INCY is INTEGER storage spacing between elements of DY

 */
extern "C" double ddot_(int* nDims, double* x, int* incx, double* y, int* incy);

//*********************************************************************************
//****************** triangular * vector ******************************************
//*********************************************************************************

/*
  subroutine dtrmv      (
  character     UPLO,
  character     TRANS,
  character     DIAG,
  integer       N,
  double precision, dimension(lda,*)    A,
  integer       LDA,
  double precision, dimension(*)        X,
  integer       INCX
  )

  DTRMV

  Purpose:

  DTRMV  performs one of the matrix-vector operations

  x := A*x,   or   x := A**T*x,

  where x is an n element vector and  A is an n by n unit, or non-unit,
  upper or lower triangular matrix.

  Parameters
  [in]  UPLO

  UPLO is CHARACTER*1
  On entry, UPLO specifies whether the matrix is an upper or
  lower triangular matrix as follows:

  UPLO = 'U' or 'u'   A is an upper triangular matrix.

  UPLO = 'L' or 'l'   A is a lower triangular matrix.

  [in]  TRANS

  TRANS is CHARACTER*1
  On entry, TRANS specifies the operation to be performed as
  follows:

  TRANS = 'N' or 'n'   x := A*x.

  TRANS = 'T' or 't'   x := A**T*x.

  TRANS = 'C' or 'c'   x := A**T*x.

  [in]  DIAG

  DIAG is CHARACTER*1
  On entry, DIAG specifies whether or not A is unit
  triangular as follows:

  DIAG = 'U' or 'u'   A is assumed to be unit triangular.

  DIAG = 'N' or 'n'   A is not assumed to be unit
  triangular.

  [in]  N

  N is INTEGER
  On entry, N specifies the order of the matrix A.
  N must be at least zero.

  [in]  A

  A is DOUBLE PRECISION array, dimension ( LDA, N )
  Before entry with  UPLO = 'U' or 'u', the leading n by n
  upper triangular part of the array A must contain the upper
  triangular matrix and the strictly lower triangular part of
  A is not referenced.
  Before entry with UPLO = 'L' or 'l', the leading n by n
  lower triangular part of the array A must contain the lower
  triangular matrix and the strictly upper triangular part of
  A is not referenced.
  Note that when  DIAG = 'U' or 'u', the diagonal elements of
  A are not referenced either, but are assumed to be unity.

  [in]  LDA

  LDA is INTEGER
  On entry, LDA specifies the first dimension of A as declared
  in the calling (sub) program. LDA must be at least
  max( 1, n ).

  [in,out]  X

  X is DOUBLE PRECISION array, dimension at least
  ( 1 + ( n - 1 )*abs( INCX ) ).
  Before entry, the incremented array X must contain the n
  element vector x. On exit, X is overwritten with the
  transformed vector x.

  [in]  INCX

  INCX is INTEGER
  On entry, INCX specifies the increment for the elements of
  X. INCX must not be zero.

*/
extern "C" void dtrmv_ (char* uplo, char* trans, char* diag, int* dim, double* triangA, int* lda, double* vecX, int* incx);

//*********************************************************************************
//****************** triangular * matrix ******************************************
//*********************************************************************************

/*
subroutine dtrmm        (       character       SIDE,
                character       UPLO,
                character       TRANSA,
                character       DIAG,
                integer         M,
                integer         N,
                double precision        ALPHA,
                double precision, dimension(lda,*)      A,
                integer         LDA,
                double precision, dimension(ldb,*)      B,
                integer         LDB
        )

Purpose:

     DTRMM  performs one of the matrix-matrix operations

        B := alpha*op( A )*B,   or   B := alpha*B*op( A ),

     where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

        op( A ) = A   or   op( A ) = A**T.

Parameters
    [in]    SIDE

              SIDE is CHARACTER*1
               On entry,  SIDE specifies whether  op( A ) multiplies B from
               the left or right as follows:

                  SIDE = 'L' or 'l'   B := alpha*op( A )*B.

                  SIDE = 'R' or 'r'   B := alpha*B*op( A ).

    [in]    UPLO

              UPLO is CHARACTER*1
               On entry, UPLO specifies whether the matrix A is an upper or
               lower triangular matrix as follows:

                  UPLO = 'U' or 'u'   A is an upper triangular matrix.

                  UPLO = 'L' or 'l'   A is a lower triangular matrix.

    [in]    TRANSA

              TRANSA is CHARACTER*1
               On entry, TRANSA specifies the form of op( A ) to be used in
               the matrix multiplication as follows:

                  TRANSA = 'N' or 'n'   op( A ) = A.

                  TRANSA = 'T' or 't'   op( A ) = A**T.

                  TRANSA = 'C' or 'c'   op( A ) = A**T.

    [in]    DIAG

              DIAG is CHARACTER*1
               On entry, DIAG specifies whether or not A is unit triangular
               as follows:

                  DIAG = 'U' or 'u'   A is assumed to be unit triangular.

                  DIAG = 'N' or 'n'   A is not assumed to be unit
                                      triangular.

    [in]    M

              M is INTEGER
               On entry, M specifies the number of rows of B. M must be at
               least zero.

    [in]    N

              N is INTEGER
               On entry, N specifies the number of columns of B.  N must be
               at least zero.

    [in]    ALPHA

              ALPHA is DOUBLE PRECISION.
               On entry,  ALPHA specifies the scalar  alpha. When  alpha is
               zero then  A is not referenced and  B need not be set before
               entry.

    [in]    A

               A is DOUBLE PRECISION array, dimension ( LDA, k ), where k is m
               when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
               Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
               upper triangular part of the array  A must contain the upper
               triangular matrix  and the strictly lower triangular part of
               A is not referenced.
               Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
               lower triangular part of the array  A must contain the lower
               triangular matrix  and the strictly upper triangular part of
               A is not referenced.
               Note that when  DIAG = 'U' or 'u',  the diagonal elements of
               A  are not referenced either,  but are assumed to be  unity.

    [in]    LDA

              LDA is INTEGER
               On entry, LDA specifies the first dimension of A as declared
               in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
               LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
               then LDA must be at least max( 1, n ).

    [in,out]    B

              B is DOUBLE PRECISION array, dimension ( LDB, N )
               Before entry,  the leading  m by n part of the array  B must
               contain the matrix  B,  and  on exit  is overwritten  by the
               transformed matrix.

    [in]    LDB

              LDB is INTEGER
               On entry, LDB specifies the first dimension of B as declared
               in  the  calling  (sub)  program.   LDB  must  be  at  least
               max( 1, m ).

 */
extern "C" void dtrmm_(char* side, char* uplo, char* transa, char* diag, int* m, int* n, double* alpha, double* A, int* lda, double* B, int* ldb);

//*********************************************************************************
//****************** matrix * matrix **********************************************
//*********************************************************************************

/*
DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
extern "C" void  dgemv_(char* trans, int* M, int* N,
                        double* alpha, double* A, int* lda,
                        double* X, int* incx,
                        double* beta, double* Y, int* incy);

//#################################################################################
//################## FACTORISATION ################################################
//#################################################################################

/*
 *  Purpose
 *  =======
 *
 *  Lapack routine for inplace calculating the Cholesky factorisation
 *  of a symmetric posdef matrix
 *
 *  DPOTRF computes the Cholesky factorization of a real symmetric
 *  positive definite matrix A.
 *
 *  The factorization has the form
 *     A = U**T * U,  if UPLO = 'U', or
 *     A = L  * L**T,  if UPLO = 'L',
 *  where U is an upper triangular matrix and L is lower triangular.
 *
 *  This is the block version of the algorithm, calling Level 3 BLAS.
 *
 *  Arguments
 *  =========
 *
 *  UPLO    (input) CHARACTER*1
 *          = 'U':  Upper triangle of A is stored;
 *          = 'L':  Lower triangle of A is stored.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
 *          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
 *          N-by-N upper triangular part of A contains the upper
 *          triangular part of the matrix A, and the strictly lower
 *          triangular part of A is not referenced.  If UPLO = 'L', the
 *          leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper
 *          triangular part of A is not referenced.
 *
 *          On exit, if INFO = 0, the factor U or L from the Cholesky
 *          factorization A = U**T*U or A = L*L**T.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *          > 0:  if INFO = i, the leading minor of order i is not
 *                positive definite, and the factorization could not be
 *                completed.
 *
 *  =====================================================================
 */
// already in gLapackHeaders.h
extern "C" void dpotrf_(char* uplo, int* orderA, double* aMat, int* nRowsA, int* info);

/*
 * Purpose
 * =======
 *
 *  DPSTRF computes the Cholesky factorization with complete
 *  pivoting of a real symmetric positive semidefinite matrix A.
 *
 *  The factorization has the form
 *     P**T * A * P = U**T * U ,  if UPLO = 'U',
 *     P**T * A * P = L  * L**T,  if UPLO = 'L',
 *  where U is an upper triangular matrix and L is lower triangular, and
 *  P is stored as vector PIV.
 *
 *  This algorithm does not attempt to check that A is positive
 *  semidefinite. This version of the algorithm calls level 3 BLAS.
 *
 * Parameters
 * ==========
 *
 *  [in]	UPLO
 *           UPLO is CHARACTER*1
 *           Specifies whether the upper or lower triangular part of the
 *           symmetric matrix A is stored.
 *           = 'U':  Upper triangular
 *           = 'L':  Lower triangular
 * [in]	N
 *           N is INTEGER
 *           The order of the matrix A.  N >= 0.
 * [in,out]	A
 *           A is DOUBLE PRECISION array, dimension (LDA,N)
 *           On entry, the symmetric matrix A.  If UPLO = 'U', the leading
 *           n by n upper triangular part of A contains the upper
 *           triangular part of the matrix A, and the strictly lower
 *           triangular part of A is not referenced.  If UPLO = 'L', the
 *           leading n by n lower triangular part of A contains the lower
 *           triangular part of the matrix A, and the strictly upper
 *           triangular part of A is not referenced.
 *
 *           On exit, if INFO = 0, the factor U or L from the Cholesky
 *           factorization as above.
 * [in]	LDA
 *           LDA is INTEGER
 *           The leading dimension of the array A.  LDA >= max(1,N).
 * [out]	PIV
 *           PIV is INTEGER array, dimension (N)
 *           PIV is such that the nonzero entries are P( PIV(K), K ) = 1.
 * [out]	RANK
 *           RANK is INTEGER
 *           The rank of A given by the number of steps the algorithm
 *           completed.
 * [in]	TOL
 *           TOL is DOUBLE PRECISION
 *           User defined tolerance. If TOL < 0, then N*U*MAX( A(K,K) )
 *           will be used. The algorithm terminates at the (K-1)st step
 *           if the pivot <= TOL.
 * [out]	WORK
 *           WORK is DOUBLE PRECISION array, dimension (2*N)
 *           Work space.
 * [out]	INFO
 *           INFO is INTEGER
 *           < 0: If INFO = -K, the K-th argument had an illegal value,
 *           = 0: algorithm completed successfully, and
 *           > 0: the matrix A is either rank deficient with computed rank
 *                as returned in RANK, or is not positive semidefinite. See
 *                Section 7 of LAPACK Working Note #161 for further
 *                information.
 */

/*
  // Example code (incomplete):
    // structures
    int dim = 3;
    gLinear::gMatrix<REAL, gLinear::COL_DENSE> A(dim,dim);
    gLinear::gMatrix<REAL, gLinear::COL_DENSE> L(dim,dim);
    prlite::ColVector<double> x(dim);

    // lets first build a rank2 choleski decomposition from rank-one updates
    L.assignToAll(0.0);
    x[0] = 1; x[1] = 3; x[2] = 4;
    cholUL(L,x);
    x[0] = 11; x[1] = 7; x[2] = 1;
    cholUL(L,x);
    A = L*L.transpose();
    // the equivalent A = L*L' matrix
    cout << "A: " << A << endl << "L: " << L << endl;
    exit(-1);
 */


//#################################################################################
//################## SOLVE LINEAR SYSTEM ##########################################
//#################################################################################


// BLAS routine for inplace solving linear system with single
// rhsides given that the lhs is triangular (in aMat).
/*
  Purpose:

  DTRSV  solves one of the systems of equations

  A*x = b,   or   A**T*x = b,

  where b and x are n element vectors and A is an n by n unit, or
  non-unit, upper or lower triangular matrix.

  No test for singularity or near-singularity is included in this
  routine. Such tests must be performed before calling this routine.

  Parameters
  [in]  UPLO

  UPLO is CHARACTER*1
  On entry, UPLO specifies whether the matrix is an upper or
  lower triangular matrix as follows:

  UPLO = 'U' or 'u'   A is an upper triangular matrix.

  UPLO = 'L' or 'l'   A is a lower triangular matrix.

  [in]  TRANS

  TRANS is CHARACTER*1
  On entry, TRANS specifies the equations to be solved as
  follows:

  TRANS = 'N' or 'n'   A*x = b.

  TRANS = 'T' or 't'   A**T*x = b.

  TRANS = 'C' or 'c'   A**T*x = b.

  [in]  DIAG

  DIAG is CHARACTER*1
  On entry, DIAG specifies whether or not A is unit
  triangular as follows:

  DIAG = 'U' or 'u'   A is assumed to be unit triangular.

  DIAG = 'N' or 'n'   A is not assumed to be unit
  triangular.

  [in]  N

  N is INTEGER
  On entry, N specifies the order of the matrix A.
  N must be at least zero.

  [in]  A

  A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  Before entry with  UPLO = 'U' or 'u', the leading n by n
  upper triangular part of the array A must contain the upper
  triangular matrix and the strictly lower triangular part of
  A is not referenced.
  Before entry with UPLO = 'L' or 'l', the leading n by n
  lower triangular part of the array A must contain the lower
  triangular matrix and the strictly upper triangular part of
  A is not referenced.
  Note that when  DIAG = 'U' or 'u', the diagonal elements of
  A are not referenced either, but are assumed to be unity.

  [in]  LDA

  LDA is INTEGER
  On entry, LDA specifies the first dimension of A as declared
  in the calling (sub) program. LDA must be at least
  max( 1, n ).

  [in,out]  X

  X is DOUBLE PRECISION array of dimension at least
  ( 1 + ( n - 1 )*abs( INCX ) ).
  Before entry, the incremented array X must contain the n
  element right-hand side vector b. On exit, X is overwritten
  with the solution vector x.

  [in]  INCX

  INCX is INTEGER
  On entry, INCX specifies the increment for the elements of
  X. INCX must not be zero.

*/
extern "C" void dtrsv_(char* uplo, char* trans, char* udiag, int* orderLHS, double* lowerTriang, int* nRowsLHS, double* rhs, int* incx);

/*
 *  Purpose
 *  =======
 *
 * Lapack routine for inplace solving linear system with multiple
 * rhsides given that the lhs is triangular (in aMat).
 *
 *  DTRTRS solves a triangular system of the form
 *
 *     A * X = B  or  A**T * X = B,
 *
 *  where A is a triangular matrix of order N, and B is an N-by-NRHS
 *  matrix.  A check is made to verify that A is nonsingular.
 *
 *  Arguments
 *  =========
 *
 *  UPLO    (input) CHARACTER*1
 *          = 'U':  A is upper triangular;
 *          = 'L':  A is lower triangular.
 *
 *  TRANS   (input) CHARACTER*1
 *          Specifies the form of the system of equations:
 *          = 'N':  A * X = B  (No transpose)
 *          = 'T':  A**T * X = B  (Transpose)
 *          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
 *
 *  DIAG    (input) CHARACTER*1
 *          = 'N':  A is non-unit triangular;
 *          = 'U':  A is unit triangular.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  NRHS    (input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the matrix B.  NRHS >= 0.
 *
 *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
 *          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
 *          upper triangular part of the array A contains the upper
 *          triangular matrix, and the strictly lower triangular part of
 *          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
 *          triangular part of the array A contains the lower triangular
 *          matrix, and the strictly upper triangular part of A is not
 *          referenced.  If DIAG = 'U', the diagonal elements of A are
 *          also not referenced and are assumed to be 1.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the right hand side matrix B.
 *          On exit, if INFO = 0, the solution matrix X.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of the array B.  LDB >= max(1,N).
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0: if INFO = -i, the i-th argument had an illegal value
 *          > 0: if INFO = i, the i-th diagonal element of A is zero,
 *               indicating that the matrix is singular and the solutions
 *               X have not been computed.
 *
 */
extern "C" void dtrtrs_(char* uplo, char* trans, char* udiag, int* orderLHS, int* nColsRHS,
                        double* lowerTriang, int* nRowsLHS, double* rhs, int* nRowsRHS, int* info);

// Lapack routine for inplace solving linear system with multiple
// rhsides given the precomputed Cholesky factorisation (in aMat).
/*
 *  Purpose
 *  =======
 *
 *  DPOTRS solves a system of linear equations A*X = B with a symmetric
 *  positive definite matrix A using the Cholesky factorization
 *  A = U**T*U or A = L*L**T computed by DPOTRF.
 *
 *  Arguments
 *  =========
 *
 *  UPLO    (input) CHARACTER*1
 *          = 'U':  Upper triangle of A is stored;
 *          = 'L':  Lower triangle of A is stored.
 *
 *  N       (input) INTEGER
 *          The order of the matrix A.  N >= 0.
 *
 *  NRHS    (input) INTEGER
 *          The number of right hand sides, i.e., the number of columns
 *          of the matrix B.  NRHS >= 0.
 *
 *  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
 *          The triangular factor U or L from the Cholesky factorization
 *          A = U**T*U or A = L*L**T, as computed by DPOTRF.
 *
 *  LDA     (input) INTEGER
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
 *          On entry, the right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *  LDB     (input) INTEGER
 *          The leading dimension of the array B.  LDB >= max(1,N).
 *
 *  INFO    (output) INTEGER
 *          = 0:  successful exit
 *          < 0:  if INFO = -i, the i-th argument had an illegal value
 *
 *  =====================================================================
 */
// already in gLapackHeaders.h
//extern "C" void dpotrs_(char* uplo, int* orderA, int* nColsB,
//double* aMat, int* nRowsA, double* bMat, int* nRowsB, int* info);

//#################################################################################
//################## INVERSES #####################################################
//#################################################################################

/*
 *  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*/
extern "C" void dtrtri_(char* uplo, char* diag,  int* n, double* aMat, int* lda, int* fail);

/*
subroutine dpotri       (       character       UPLO,
                integer         N,
                double precision, dimension( lda, * )   A,
                integer         LDA,
                integer         INFO
        )

DPOTRI

Download DPOTRI + dependencies [TGZ] [ZIP] [TXT]

Purpose:

     DPOTRI computes the inverse of a real symmetric positive definite
     matrix A using the Cholesky factorization A = U**T*U or A = L*L**T
     computed by DPOTRF.

Parameters
    [in]    UPLO

              UPLO is CHARACTER*1
              = 'U':  Upper triangle of A is stored;
              = 'L':  Lower triangle of A is stored.

    [in]    N

              N is INTEGER
              The order of the matrix A.  N >= 0.

    [in,out]    A

              A is DOUBLE PRECISION array, dimension (LDA,N)
              On entry, the triangular factor U or L from the Cholesky
              factorization A = U**T*U or A = L*L**T, as computed by
              DPOTRF.
              On exit, the upper or lower triangle of the (symmetric)
              inverse of A, overwriting the input factor U or L.

    [in]    LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

    [out]   INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, the (i,i) element of the factor U or L is
                    zero, and the inverse could not be computed.

 */
extern "C" void dpotri_(char* uplo,
                       int* dim,
                       double* A,
                       int* lda,
                       int* fail);

//#################################################################################
//#################################################################################
//#################################################################################

#endif // LAPACKBLAS_HPP
