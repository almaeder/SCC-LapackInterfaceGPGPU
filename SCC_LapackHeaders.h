/*
 * SCC_LapackHeaders.h
 *
 * LAPACK function prototypes for LAPACK routines used by the collection of classes
 * contained in LapackInterface
 *
 *
 *  Created on: Oct 25, 2017
 *      Author: anderson
 *
 *
 *  Updated : July 27, 2018 (CRA)
 *  Updated : Dec. 09, 2023 (CRA)
 */

/*
#############################################################################
#
# Copyright  2015-2018 Chris Anderson
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

#ifndef SCC_LAPACK_HEADERS_
#define SCC_LAPACK_HEADERS_

extern "C" double dlamch_(char* CMACH);

// Matrix-Vector

// double general
extern "C" void dgemv_(char* TRANS, long* M, long* N, double* alpha, double* Aptr,
				long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

// double general  banded
extern "C" void dgbmv_(char* TRANS, long* M, long* N, long* kl, long* ku, double* alpha, double* Aptr,
                       long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

// double symmetric banded
extern "C" void dsbmv_(char* UPLO, long* N, long* K, double* alpha, double* Aptr,
						long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

// complex general
extern "C" void zgemv_(char* TRANS, long* M, long* N, double* alpha, double* Aptr,
                       long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

// complex general banded
extern "C" void zgbmv_(char* TRANS, long* M, long* N, long* kl, long* ku, double* alpha, double* Aptr,
                       long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

// Matrix-Matrix multiplication

// double general
extern "C" void dgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);

// complex general

extern "C" void zgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);

// Solution of systems

// double general
extern "C" void dgesvx_(char* FACT, char* TRANS, long* N, long* NRHS, double* Aptr, long* LDA, double* AFptr, long* LDAF, long* IPIVptr,
		                char* EQED, double* Rptr, double* Cptr, double* Bptr, long* LDB, double* Xptr, long* LDX,  double* RCOND,
						double* FERR, double* BERR, double* WORKptr,
		                long* IWORKptr, long* INFO);

// double general banded
extern "C" void dgbsvx_(char* FACT, char* TRANS, long* N, long* KL, long*  	KU, long*  	NRHS,
						double* AB, long* LDAB, double* AFB, long* LDAFB, long* IPIV, char* EQUED,
						double* R, double* C, double* B, long* LDB, double* X, long* LDX, double* RCOND, double* FERR,
						double*   BERR, double *WORK, long*	IWORK, long* INFO);

// complex general banded
extern "C" void zgbsvx_(char* FACT, char* TRANS, long* N, long* KL, long* KU, long* NRHS, double* Aptr, long* LDA, double* AFptr, long* LDAF,
						long* IPIVptr,char* EQED, double* Rptr, double* Cptr, double* Bptr, long* LDB, double* Xptr, long* LDX,  double* RCOND,
						double* FERR, double* BERR, double* WORKptr, double* RWORKptr, long* INFO);

// double symmetric
extern "C"  void dsyevx_(char*JOBZ, char* RANGE, char* UPLO,long* N, double* A, long* LDA, double* VL, double* VU,
                         long*   IL, long*   IU, double*   ABSTOL, long*   M, double* W, double* Z, long* LDZ,
                         double* WORK, long * LWORK, long* IWORK, long* IFAIL, long* INFO);

// double symmetric positive definite (Choleski factorization)
extern "C" void dposv_(char* UPLO, long* N, long* NRHS, double* Aptr, long* LDA, double* Bptr, long* LDB, long* INFO);

// complex general 

extern "C" void zgesvx_(char* FACT, char* TRANS, long* N, long* NRHS, double* Aptr, long* LDA, double* AFptr, long* LDAF, long* IPIVptr,
		   char* EQED, double* Rptr, double* Cptr, double* Bptr, long* LDB, double* Xptr, long* LDX,  double* RCOND,
	       double* FERR, double* BERR, double* WORKptr,
		   double* RWORKptr, long* INFO);

// double triangular

extern "C" void dtrtrs_(char* UPLO, char* TRANS, char* DIAG, long* N, long* NRHS, double* A, long* LDA, double* B,
                        long* LDB, long* INFO);

// Least squares solution

// double general
extern "C" void dgelsy_(long* M, long* N, long* NRHS, double* APtr,
                        long* LDA, double* BPtr, long* LDB, long* JPVT, double* RCOND,
                        long* RANK, double* WORK, long* LWORK, long* INFO);

// Eigensystem computations

// double general
extern "C" void dgeevx_(char*  BALANC, char* JOBVL, char* JOBVR, char* SENSE,
                        long*  N, double* A, long* LDA, double* WR, double* WI, double* VL, long* LDVL,
                        double* VR, long*  	LDVR, long* ILO, long* IHI, double* SCALE, double* ABNRM,
                        double* RCONDE, double* RCONDV, double* WORK, long* LWORK, long* IWORK, long* INFO);
// double general
extern "C" int   dgeev_(char* JOBVL, char* JOBVR, long* N,double* A, long* LDA,double* WR, double* WI, double* VL,
                        long* LDVL, double* VR, long* LDVR, double* WORK, long* LWORK, long* INFO);
     			
             			
// double symmetric
extern "C" void dsyev_(char* JOBZ,char* UPLO, long*N, double* Aptr, long* LDA, double* Wptr,
						double* WORKptr, long* LWORK, long* INFO);

// complex Hermitian
extern "C"  void zhpevx_(char*JOBZ, char* RANGE, char* UPLO,long* N, double* AP, double* VL, double* VU,
                        long*   IL, long*   IU, double*   ABSTOL, long*   M, double* W, double* Z, long* LDZ,
             			double* WORK, double* RWORK, long* IWORK, long* IFAIL, long* INFO);

// complex general
extern "C"  void zgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, long* n,
double* A, long* 	LDA, double*  w, double*  vl, long*   ldvl, double*  vr, long* ldvr,
long* ilo, long* ihi, double* scale, double* abnrm, double* rconde, double* rcondv,
double*	work, long*	lwork, double* rwork, long* info);
             			
// double symmetric tri-diagonal
extern "C" int dsteqr_(char *compz, long *n, double *d, double *e, double *z,
		               long *ldz, double *work, long *info);

// double symmetric tri-diagonal
extern "C" int dstevx_(char *jobz, char *range, long *n, double *
                       d__, double *e, double *vl, double *vu, long *il,
                       long *iu, double *abstol, long *m, double *w,
                       double *z__, long *ldz, double *work, long *iwork,
                       long *ifail, long *info);

// double symmetric tri-diagonal
extern "C" int dstebz_(char *range, char *order, long *n, double
                       *vl, double *vu, long *il, long *iu, double *abstol,
                       double *d__, double *e, long *m, long *nsplit, double *w, long *iblock,
                       long *isplit, double *work, long *iwork, long *info);

// QR factorization

// double general (create QR factors, Q stored as elementary reflectors)
extern "C" void dgeqrf_(long* M, long* N, double* Aptr, long* LDA, double* TAU, double* WORK, long* LWORK, long* INFO);

// double general (Create Q of QR factorization from output of dgeqrf)
extern "C" void dorgqr_(long* M, long* N, long* K, double* Aptr, long* LDA, double* TAU, double* WORK, long* LWORK, long* INFO);

// double general (applies Q of QR factoriozation, Q stored as elementary reflectors)
extern "C" void dormqr_(char* SIDE, char* TRANS, long* M, long* N, long* K, double* A,
                        long* LDA, double* TAU, double* C, long* LDC, double* WORK, long* LWORK, long* INFO);

// SVD

// double general
extern "C" void dgesvd_(char* JOBU,char* JOBVT, long* M, long* N, double* APtr, long* LDA, double* SPtr, double* UPtr, long* LDU, double* VTPtr, long* LDVT,
                       double* WORKtmp, long* LWORK, long* INFO);

// LU factorization



// double general tri-diaognal (create factors)
extern "C" int dgttrf_(long* N, double* DL,double* D, double* DU, double* DU2, long* IPIV,long* INFO);

// double general tri-diagonal (create solution using factors)
extern "C" int dgttrs_(char* TRANS, long* N, long* NRHS,double* DL,double* D, double* DU, double* DU2, long* IPIV,double* B, long* LDB, long*INFO);



#endif /* SCC_LAPACKHEADERS_H_ */



