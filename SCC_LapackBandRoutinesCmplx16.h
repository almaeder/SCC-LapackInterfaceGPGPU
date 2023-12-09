/*
 * SCC_LapackBandRoutinesCmplx16.h
 *
 *  Created on: Dec. 3, 2023
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright  2023 Chris Anderson
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

#include "SCC_LapackHeaders.h"
#include "SCC_LapackBandMatrixCmplx16.h"

#ifndef SCC_LAPACK_BAND_ROUTINES_CMPLX_16_
#define SCC_LAPACK_BAND_ROUTINES_CMPLX_16_


namespace SCC
{

/*
 ZGBSVX uses the LU factorization to compute the solution to a complex
 system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
 where A is a band matrix of order N with KL subdiagonals and KU
 superdiagonals, and X and B are N-by-NRHS matrices.

 Error bounds on the solution and a condition estimate are also
 provided.
*/

class ZGBSVX
{
public:

	ZGBSVX()
	{
		initialize();
	}

	void initialize()
	{
    RCOND = 0.0;
	FERR.clear();
	BERR.clear();

	A.initialize();
	AF.initialize();
	}

    void applyInverse(const LapackBandMatrixCmplx16& A,std::vector <std::complex<double>>& b)
	{
    	    assert(A.sizeCheck(A.N,(long)b.size()));
	        double* bptr =  &(reinterpret_cast<double(&)[2]>(b[0])[0]);
			applyInverse(A,bptr);
	}

    void applyInverse(const LapackBandMatrixCmplx16& A,LapackMatrixCmplx16& b)
	{
    	    assert(A.sizeCheck(A.N,b.rows));
    		applyInverse(A,b.mData.dataPtr,b.cols);
	}

	void applyInverse(const LapackBandMatrixCmplx16& S, double* b, long NRHS = 1)
	{
		//assert(A.sizeCheck(A.rows,A.cols));

		char FACT  = 'E'; // Equilibrate, then factor
		char TRANS = 'N'; // No transpose

		long N     = S.N;
		long KL    = S.kl;
		long KU    = S.ku;

		//
		// Duplicate input matrix (since this zgbsvx overwrites input matrix)
		// and allocate temporaries
        //

		long LDAB  = KL + KU + 1;
		long LDAF  = 2*KL + KU + 1;

		this->A.initialize(S);
		this->AF.initialize(LDAF,N);

		double* Aptr  =  A.cmplxMdata.mData.dataPtr;
		double* AFptr =  AF.mData.dataPtr;

		std::vector <long >   IPIV(N);
		long* IPIVptr      = &IPIV[0];

		char  EQED;

		std::vector<double>   R(N);
		double* Rptr  = &R[0];

		std::vector<double>    C(N);
		double* Cptr  =  &C[0];

		std::vector<double>   B(2*N*NRHS);
		double* Bptr  =       &B[0];
	    long LDB      =          N;


		// b will be overwritten with the solution
	    // so no need to declare X separately

		double* Xptr   = b;
		long LDX       = N;

		FERR.resize(NRHS);
		BERR.resize(NRHS);

		std::vector<double>   WORK(4*N);
		double* WORKptr     = &WORK[0];

		std::vector<double>  RWORK(2*N);
		double* RWORKptr   = &RWORK[0];

		long   INFO = 0;

		// Assign right hand side to B

		for(long i = 0; i < 2*N*NRHS; i++)
		{
			Bptr[i] = b[i];
		}


		zgbsvx_(&FACT, &TRANS, &N, &KL, &KU, &NRHS, Aptr, &LDAB, AFptr, &LDAF, IPIVptr,
		        &EQED, Rptr, Cptr, Bptr,&LDB, Xptr, &LDX, &RCOND,
				&FERR[0], &BERR[0], WORKptr,RWORKptr, &INFO);


		if(INFO != 0)
        {
        std::cerr << "dgesvx  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }
	}

/*
*  RCOND is DOUBLE PRECISION
*  The estimate of the reciprocal condition number of the matrix
*  A after equilibration (if done).  If RCOND is less than the
*  machine precision (in particular, if RCOND = 0), the matrix
*  is singular to working precision.  This condition is
*  indicated by a return code of INFO > 0.
*/

	double getReciprocalConditionNumber()
	{
		return RCOND;
	}

/*
    FERR is DOUBLE PRECISION array, dimension (NRHS)
    The estimated forward error bound for each solution vector
    X(j) (the j-th column of the solution matrix X).
    If XTRUE is the true solution corresponding to X(j), FERR(j)
    is an estimated upper bound for the magnitude of the largest
    element in (X(j) - XTRUE) divided by the magnitude of the
    largest element in X(j).  The estimate is as reliable as
    the estimate for RCOND, and is almost always a slight
    overestimate of the true error.
*/

	double getSolutionErrorEst()
	{
		return FERR[0];
	}

	std::vector<double> getMultipleSolutionErrorEst()
	{
		return FERR;
	}

/*
	BERR is DOUBLE PRECISION array, dimension (NRHS)
	The componentwise relative backward error of each solution
	vector X(j) (i.e., the smallest relative change in
	any element of A or B that makes X(j) an exact solution).
*/
	double getSolutionBackwardErrorEst()
	{
		return BERR[0];
	}

	std::vector<double> getMultipleSolutionBackwardErrorEst()
	{
		return BERR;
	}



    double          RCOND;
	std::vector<double>    FERR;
	std::vector<double>    BERR;

	LapackBandMatrixCmplx16    A;
	LapackMatrixCmplx16       AF;

};


} // SCC namespace


//
// LAPACK documentation
//
/*
subroutine zgbsvx 	( 	character  	fact,
		character  	trans,
		integer  	n,
		integer  	kl,
		integer  	ku,
		integer  	nrhs,
		complex*16, dimension( ldab, * )  	ab,
		integer  	ldab,
		complex*16, dimension( ldafb, * )  	afb,
		integer  	ldafb,
		integer, dimension( * )  	ipiv,
		character  	equed,
		double precision, dimension( * )  	r,
		double precision, dimension( * )  	c,
		complex*16, dimension( ldb, * )  	b,
		integer  	ldb,
		complex*16, dimension( ldx, * )  	x,
		integer  	ldx,
		double precision  	rcond,
		double precision, dimension( * )  	ferr,
		double precision, dimension( * )  	berr,
		complex*16, dimension( * )  	work,
		double precision, dimension( * )  	rwork,
		integer  	info
	)

Purpose:

     ZGBSVX uses the LU factorization to compute the solution to a complex
     system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
     where A is a band matrix of order N with KL subdiagonals and KU
     superdiagonals, and X and B are N-by-NRHS matrices.

     Error bounds on the solution and a condition estimate are also
     provided.

Description:

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
              = 'C':  A**H * X = B  (Conjugate transpose)

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

              AB is COMPLEX*16 array, dimension (LDAB,N)
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

              AFB is COMPLEX*16 array, dimension (LDAFB,N)
              If FACT = 'F', then AFB is an input argument and on entry
              contains details of the LU factorization of the band matrix
              A, as computed by ZGBTRF.  U is stored as an upper triangular
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
              as computed by ZGBTRF; row i of the matrix was interchanged
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

              B is COMPLEX*16 array, dimension (LDB,NRHS)
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

              X is COMPLEX*16 array, dimension (LDX,NRHS)
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

              WORK is COMPLEX*16 array, dimension (2*N)

    [out]	RWORK

              RWORK is DOUBLE PRECISION array, dimension (MAX(1,N))
              On exit, RWORK(1) contains the reciprocal pivot growth
              factor norm(A)/norm(U). The "max absolute element" norm is
              used. If RWORK(1) is much less than 1, then the stability
              of the LU factorization of the (equilibrated) matrix A
              could be poor. This also means that the solution X, condition
              estimator RCOND, and forward error bound FERR could be
              unreliable. If factorization fails with 0<INFO<=N, then
              RWORK(1) contains the reciprocal pivot growth factor for the
              leading INFO columns of A.

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

Author
    Univ. of Tennessee
    Univ. of California Berkeley
    Univ. of Colorado Denver
    NAG Ltd.
*/

#endif /* SCC_LapackBandRoutines__ */
