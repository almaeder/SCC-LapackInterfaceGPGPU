/*
 * SCC_LapackHeaders.h
 *
 * LAPACK function prototypes for LAPACK routines used by the collection of classes
 * contained in SCC::LapackMatrixRoutines.h
 *
 *
 *  Created on: Oct 25, 2017
 *      Author: anderson
 *
 *
 *  Last updated : July 27, 2018 (CRA)
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

extern "C" void dgemv_(char* TRANS, long* M, long* N, double* alpha, double* Aptr,
long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

extern "C" void dsyev_(char* JOBZ,char* UPLO, long*N, double* Aptr, long* LDA, double* Wptr,
double* WORKptr, long* LWORK, long* INFO);


extern "C" void dgesvx_(char* FACT, char* TRANS, long* N, long* NRHS, double* Aptr, long* LDA, double* AFptr, long* LDAF, long* IPIVptr,
		                char* EQED, double* Rptr, double* Cptr, double* Bptr, long* LDB, double* Xptr, long* LDX,  double* RCOND,
						double* FERR, double* BERR, double* WORKptr,
		                long* IWORKptr, long* INFO);

extern "C" void dgeqrf_(long* M, long* N, double* Aptr, long* LDA, double* TAU, double* WORK, long* LWORK, long* INFO);

extern "C" void dorgqr_(long* M, long* N, long* K, double* Aptr, long* LDA, double* TAU, double* WORK, long* LWORK, long* INFO);

extern "C" void dormqr_(char* SIDE, char* TRANS, long* M, long* N, long* K, double* A,
                        long* LDA, double* TAU, double* C, long* LDC, double* WORK, long* LWORK, long* INFO);

extern "C" void dtrtrs_(char* UPLO, char* TRANS, char* DIAG, long* N, long* NRHS, double* A, long* LDA, double* B,
                        long* LDB, long* INFO);


extern "C" void dgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);

extern "C" void dgesvd_(char* JOBU,char* JOBVT, long* M, long* N, double* APtr, long* LDA, double* SPtr, double* UPtr, long* LDU, double* VTPtr, long* LDVT,
                       double* WORKtmp, long* LWORK, long* INFO);

extern "C" void dgelsy_(long* M, long* N, long* NRHS, double* APtr,
                        long* LDA, double* BPtr, long* LDB, long* JPVT, double* RCOND,
                        long* RANK, double* WORK, long* LWORK, long* INFO);


extern "C" void dgeevx_(char*  BALANC, char* JOBVL, char* JOBVR, char* SENSE,
                        long*  N, double* A, long* LDA, double* WR, double* WI, double* VL, long* LDVL,
                        double* VR, long*  	LDVR, long* ILO, long* IHI, double* SCALE, double* ABNRM,
                        double* RCONDE, double* RCONDV, double* WORK, long* LWORK, long* IWORK, long* INFO);

#endif /* SCC_LAPACKHEADERS_H_ */



