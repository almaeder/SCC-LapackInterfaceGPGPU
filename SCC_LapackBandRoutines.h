/*
 * SCC_LapackBandRoutines.h
 *
 *  Created on: Aug 17, 2018
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright  2018 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/



#include <vector>
#include <iostream>
#include <cstring>

#include "SCC_LapackBandMatrix.h"

#ifndef SCC_LAPACK_BAND_ROUTINES_
#define SCC_LAPACK_BAND_ROUTINES_

// Headers for band matrix routines

extern "C" void dgbsvx_(char* FACT, char* TRANS, long* N, long* KL, long*  	KU, long*  	NRHS,
				double* AB, long* LDAB, double* AFB, long* LDAFB, long* IPIV, char* EQUED,
				double* R, double* C, double* B, long* LDB, double* X, long* LDX, double* RCOND, double* FERR,
				double*   BERR, double *WORK, long*	IWORK, long* INFO);



namespace SCC
{

/*
*  DGBSVX uses the LU factorization to compute the solution to a real
*  system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
*  where A is a band matrix of order N with KL subdiagonals and KU
*  superdiagonals, and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*/

class DGBSVX
{
	public:

	DGBSVX()
	{
	initialize();
	}

	DGBSVX(const DGBSVX& dgbsvx)
	{
		initialize(dgbsvx);
	}

	void initialize()
	{
		FACT  = 'E'; // E = equilibrate
		EQUED = 'B';
		RCOND = 0.0;
	    FERR  = 0.0;
	    BERR  = 0.0;

	    ABmatrix.initialize();
	    AFB.initialize();
	    IPIV.clear();
	    R.clear();
	    C.clear();
	    X.clear();
	    WORK.clear();
	    IWORK.clear();
	}
	void initialize(const DGBSVX& dgbsvx)
	{
		FACT  = dgbsvx.FACT;
		EQUED = dgbsvx.EQUED;

		RCOND = dgbsvx.RCOND;
	    FERR  = dgbsvx.FERR;
	    BERR  = dgbsvx.BERR;

	    // For caching factors

	    ABmatrix.initialize(dgbsvx.ABmatrix);
	    AFB.initialize(dgbsvx.AFB);
	    IPIV = dgbsvx.IPIV;
	    R    = dgbsvx.R;
	    C    = dgbsvx.C;
	    X    = dgbsvx.X;
	    WORK = dgbsvx.WORK;
	    IWORK= dgbsvx.IWORK;
	}

	void setEquilibration(bool val = true)
	{
	if(val) {FACT  = 'E';}
	else    {FACT  = 'N';}
	}

	void clearEquilibration()
	{
	FACT  = 'N';
	}

	 /*
          setEquilibrationType() specifies the form of equilibration that was done.
          = 'N':  No equilibration (always true if FACT = 'N').
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R).
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C).
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C).
    */


	void setEquilibrationType(char T)
	{
	EQUED = T;
	}

	void applyInverse(SCC::LapackBandMatrix& S, std::vector<double>& f)
	{
	applyInverse(S,&f[0]);
	}

	void applyInverse(SCC::LapackBandMatrix& S, double* f)
	{
	createFactors(S);
	applyInverse(f);
	return;

    //char FACT  = 'E'; // E Or N for no-equilibration

    char TRANS = 'N';

    long N     = S.N;
    long KL    = S.kl;
    long KU    = S.ku;

    long NRHS  = 1;
    double* AB = S.getDataPointer();
    long LDAB  = KL + KU + 1;
    long LDAFB = 2*KL+KU+1;


    SCC::LapackMatrix AFB(LDAFB,N);
    std::vector<long>      IPIV(N);


    std::vector<double> R(N);
    std::vector<double> C(N);

    double* Bptr = &f[0];
    long LDB     = N;

    std::vector<double> X(N); // Now an output argument

    long LDX     = N;

    //double RCOND = 1.0;
    //double FERR;
    //double BERR;

    std::vector<double> WORK(3*N);
    std::vector<long>    IWORK(N);

    long INFO = 0;

    // coeffRHS

    dgbsvx_(&FACT, &TRANS, &N, &KL, &KU, &NRHS, AB, &LDAB, AFB.getDataPointer(), &LDAFB, &IPIV[0], &EQUED,
    		&R[0], &C[0], Bptr, &LDB, &X[0], &LDX,&RCOND, &FERR,
			&BERR, &WORK[0], &IWORK[0], &INFO);

    if(INFO != 0)
    {
    	if(INFO < 0) { std::cout << -INFO << "argument to dgbsvx had an illegal value " << std::endl;}
    	if(INFO > 0) { std::cout <<  "Error in dgbsvs INFO = " << INFO << " " << " N " << N << std::endl;}
    }

    // Capture the solution

    // f = X;

    std::memcpy(&f[0],&X[0],N*sizeof(double));
	}


	void applyInverse(double* f)
	{
    //char FACT  = 'E'; // E Or N for no-equilibration

    char TRANS = 'N';

    long N     = ABmatrix.N;
    long KL    = ABmatrix.kl;
    long KU    = ABmatrix.ku;

    long NRHS  = 1;
    double* AB = ABmatrix.getDataPointer();
    long LDAB  = KL + KU + 1;
    long LDAFB = 2*KL+KU+1;

    double* Bptr = &f[0];
    long LDB     = N;
    long LDX     = N;

    long INFO = 0;

    // coeffRHS
    char FACT_TYPE = 'F';

    dgbsvx_(&FACT_TYPE, &TRANS, &N, &KL, &KU, &NRHS, AB, &LDAB, AFB.getDataPointer(), &LDAFB, &IPIV[0], &EQUED,
    		&R[0], &C[0], Bptr, &LDB, &X[0], &LDX,&RCOND, &FERR,
			&BERR, &WORK[0], &IWORK[0], &INFO);

    if(INFO != 0)
    {
    	if(INFO < 0) { std::cout << -INFO << "argument to dgbsvx had an illegal value " << std::endl;}
    	if(INFO > 0) { std::cout <<  "Error in dgbsvs INFO = " << INFO << " " << " N " << N << std::endl;}
    }

    // Capture the solution

    // f = X;

    std::memcpy(&f[0],&X[0],N*sizeof(double));
	}



	void createFactors(const SCC::LapackBandMatrix& S)
	{

    // non-default equilibration set before call to this method

    //char FACT =  'N':  The matrix A will be copied to AFB and factored.
    //char FACT  = 'E':  The matrix A will be equilibrated if necessary, then
    //                   copied to AFB and factored.

    long NRHS  = 0;  // Just factoring so no right hand sides

    char TRANS = 'N';

    long N     = S.N;
    long KL    = S.kl;
    long KU    = S.ku;

    ABmatrix.initialize(S);
    double* AB = ABmatrix.getDataPointer();
    long LDAB  = KL + KU + 1;
    long LDAFB = 2*KL+KU+1;

    AFB.initialize(LDAFB,N);
    IPIV.resize(N);
    R.resize(N);
    C.resize(N);
    X.resize(N);

    WORK.resize(3*N);
    IWORK.resize(N);

    double doubleNull = 0.0;
    double* Bptr      = &doubleNull;
    long LDB          = N;

    double xNull      = 0.0;
    double* xPtr      = &xNull;

    long LDX           = N;

    //double RCOND = 1.0;
    //double FERR;
    //double BERR;

    long INFO = 0;


    dgbsvx_(&FACT, &TRANS, &N, &KL, &KU, &NRHS, AB, &LDAB, AFB.getDataPointer(), &LDAFB, &IPIV[0], &EQUED,
    		&R[0], &C[0], Bptr, &LDB, xPtr, &LDX,&RCOND, &FERR,
			&BERR, &WORK[0], &IWORK[0], &INFO);

    if(INFO != 0)
    {
    	if(INFO < 0) { std::cout << -INFO << "argument to dgbsvx had an illegal value " << std::endl;}
    	if(INFO > 0) { std::cout <<  "Error in dgbsvs INFO = " << INFO << " " << " N " << N << std::endl;}
    }

	}


	double getReciprocalCondNumber()
	{
	return RCOND;
	}

/*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*/
	double getForwardErrEstimate()
	{
	return FERR;
	}

/*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*/
	double getBackwardErrEstimate()
	{
	return BERR;
	}


	char   FACT;
	char   EQUED;

	double RCOND;
    double FERR;
    double BERR;

    // For caching factors

    SCC::LapackBandMatrix ABmatrix;
    SCC::LapackMatrix          AFB;
    std::vector<long>         IPIV;
    std::vector<double>          R;
    std::vector<double>          C;
    std::vector<double>          X;
    std::vector<double>        WORK;
    std::vector<long>         IWORK;
};

}



/*
      The following steps are performed by this subroutine:

 1. If FACT = 'E', real scaling factors are computed to equilibrate
    the system:
       TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
       TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
       TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
    Whether or not the system will be equilibrated depends on the
    scaling of the matrix A, but if equilibration is used, A is
    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
    or diag(C)*B (if TRANS = 'T' or 'C').

 2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
    matrix A (after equilibration if FACT = 'E') as
       A = L * U,
    where L is a product of permutation and unit lower triangular
    matrices with KL subdiagonals, and U is upper triangular with
    KL+KU superdiagonals.

 3. If some U(i,i)=0, so that U is exactly singular, then the routine
    returns with INFO = i. Otherwise, the factored form of A is used
    to estimate the condition number of the matrix A.  If the
    reciprocal of the condition number is less than machine precision,
    INFO = N+1 is returned as a warning, but the routine still goes on
    to solve for X and compute error bounds as described below.

 4. The system of equations is solved for X using the factored form
    of A.

 5. Iterative refinement is applied to improve the computed solution
    matrix and calculate error bounds and backward error estimates
    for it.

 6. If equilibration was used, the matrix X is premultiplied by
    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
    that it solves the original system before equilibration.

Parameters
[in]	FACT

          FACT is CHARACTER*1
          Specifies whether or not the factored form of the matrix A is
          supplied on entry, and if not, whether the matrix A should be
          equilibrated before it is factored.
          = 'F':  On entry, AFB and IPIV contain the factored form of
                  A.  If EQUED is not 'N', the matrix A has been
                  equilibrated with scaling factors given by R and C.
                  AB, AFB, and IPIV are not modified.
          = 'N':  The matrix A will be copied to AFB and factored.
          = 'E':  The matrix A will be equilibrated if necessary, then
                  copied to AFB and factored.

[in]	TRANS

          TRANS is CHARACTER*1
          Specifies the form of the system of equations.
          = 'N':  A * X = B     (No transpose)
          = 'T':  A**T * X = B  (Transpose)
          = 'C':  A**H * X = B  (Transpose)

[in]	N

          N is INTEGER
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0.

[in]	KL

          KL is INTEGER
          The number of subdiagonals within the band of A.  KL >= 0.

[in]	KU

          KU is INTEGER
          The number of superdiagonals within the band of A.  KU >= 0.

[in]	NRHS

          NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrices B and X.  NRHS >= 0.

[in,out]	AB

          AB is DOUBLE PRECISION array, dimension (LDAB,N)
          On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
          The j-th column of A is stored in the j-th column of the
          array AB as follows:
          AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)

          If FACT = 'F' and EQUED is not 'N', then A must have been
          equilibrated by the scaling factors in R and/or C.  AB is not
          modified if FACT = 'F' or 'N', or if FACT = 'E' and
          EQUED = 'N' on exit.

          On exit, if EQUED .ne. 'N', A is scaled as follows:
          EQUED = 'R':  A := diag(R) * A
          EQUED = 'C':  A := A * diag(C)
          EQUED = 'B':  A := diag(R) * A * diag(C).

[in]	LDAB

          LDAB is INTEGER
          The leading dimension of the array AB.  LDAB >= KL+KU+1.

[in,out]	AFB

          AFB is DOUBLE PRECISION array, dimension (LDAFB,N)
          If FACT = 'F', then AFB is an input argument and on entry
          contains details of the LU factorization of the band matrix
          A, as computed by DGBTRF.  U is stored as an upper triangular
          band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
          and the multipliers used during the factorization are stored
          in rows KL+KU+2 to 2*KL+KU+1.  If EQUED .ne. 'N', then AFB is
          the factored form of the equilibrated matrix A.

          If FACT = 'N', then AFB is an output argument and on exit
          returns details of the LU factorization of A.

          If FACT = 'E', then AFB is an output argument and on exit
          returns details of the LU factorization of the equilibrated
          matrix A (see the description of AB for the form of the
          equilibrated matrix).

[in]	LDAFB

          LDAFB is INTEGER
          The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.

[in,out]	IPIV

          IPIV is INTEGER array, dimension (N)
          If FACT = 'F', then IPIV is an input argument and on entry
          contains the pivot indices from the factorization A = L*U
          as computed by DGBTRF; row i of the matrix was interchanged
          with row IPIV(i).

          If FACT = 'N', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = L*U
          of the original matrix A.

          If FACT = 'E', then IPIV is an output argument and on exit
          contains the pivot indices from the factorization A = L*U
          of the equilibrated matrix A.

[in,out]	EQUED

          EQUED is CHARACTER*1
          Specifies the form of equilibration that was done.
          = 'N':  No equilibration (always true if FACT = 'N').
          = 'R':  Row equilibration, i.e., A has been premultiplied by
                  diag(R).
          = 'C':  Column equilibration, i.e., A has been postmultiplied
                  by diag(C).
          = 'B':  Both row and column equilibration, i.e., A has been
                  replaced by diag(R) * A * diag(C).
          EQUED is an input argument if FACT = 'F'; otherwise, it is an
          output argument.

[in,out]	R

          R is DOUBLE PRECISION array, dimension (N)
          The row scale factors for A.  If EQUED = 'R' or 'B', A is
          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
          is not accessed.  R is an input argument if FACT = 'F';
          otherwise, R is an output argument.  If FACT = 'F' and
          EQUED = 'R' or 'B', each element of R must be positive.

[in,out]	C

          C is DOUBLE PRECISION array, dimension (N)
          The column scale factors for A.  If EQUED = 'C' or 'B', A is
          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
          is not accessed.  C is an input argument if FACT = 'F';
          otherwise, C is an output argument.  If FACT = 'F' and
          EQUED = 'C' or 'B', each element of C must be positive.

[in,out]	B

          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the right hand side matrix B.
          On exit,
          if EQUED = 'N', B is not modified;
          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
          diag(R)*B;
          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
          overwritten by diag(C)*B.

[in]	LDB

          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).

[out]	X

          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
          to the original system of equations.  Note that A and B are
          modified on exit if EQUED .ne. 'N', and the solution to the
          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
          and EQUED = 'R' or 'B'.

[in]	LDX

          LDX is INTEGER
          The leading dimension of the array X.  LDX >= max(1,N).

[out]	RCOND

          RCOND is DOUBLE PRECISION
          The estimate of the reciprocal condition number of the matrix
          A after equilibration (if done).  If RCOND is less than the
          machine precision (in particular, if RCOND = 0), the matrix
          is singular to working precision.  This condition is
          indicated by a return code of INFO > 0.

[out]	FERR

          FERR is DOUBLE PRECISION array, dimension (NRHS)
          The estimated forward error bound for each solution vector
          X(j) (the j-th column of the solution matrix X).
          If XTRUE is the true solution corresponding to X(j), FERR(j)
          is an estimated upper bound for the magnitude of the largest
          element in (X(j) - XTRUE) divided by the magnitude of the
          largest element in X(j).  The estimate is as reliable as
          the estimate for RCOND, and is almost always a slight
          overestimate of the true error.

[out]	BERR

          BERR is DOUBLE PRECISION array, dimension (NRHS)
          The componentwise relative backward error of each solution
          vector X(j) (i.e., the smallest relative change in
          any element of A or B that makes X(j) an exact solution).

[out]	WORK

          WORK is DOUBLE PRECISION array, dimension (3*N)
          On exit, WORK(1) contains the reciprocal pivot growth
          factor norm(A)/norm(U). The "max absolute element" norm is
          used. If WORK(1) is much less than 1, then the stability
          of the LU factorization of the (equilibrated) matrix A
          could be poor. This also means that the solution X, condition
          estimator RCOND, and forward error bound FERR could be
          unreliable. If factorization fails with 0<INFO<=N, then
          WORK(1) contains the reciprocal pivot growth factor for the
          leading INFO columns of A.

[out]	IWORK

          IWORK is INTEGER array, dimension (N)

[out]	INFO

          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, and i is
                <= N:  U(i,i) is exactly zero.  The factorization
                       has been completed, but the factor U is exactly
                       singular, so the solution and error bounds
                       could not be computed. RCOND = 0 is returned.
                = N+1: U is nonsingular, but RCOND is less than machine
                       precision, meaning that the matrix is singular
                       to working precision.  Nevertheless, the
                       solution and error bounds are computed because
                       there are a number of situations where the
                       computed solution can be more accurate than the
                       value of RCOND would suggest.

 */
#endif /* SCC_LapackBandRoutines__ */
