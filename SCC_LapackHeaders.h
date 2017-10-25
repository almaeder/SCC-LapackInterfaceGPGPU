/*
 * SCC_LapackHeaders.h
 *
 *  Created on: Oct 25, 2017
 *      Author: anderson
 */

#ifndef SCC_LAPACKHEADERS_H_
#define SCC_LAPACKHEADERS_H_

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

extern "C" void dgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);


extern "C" void dgesvd_(char* JOBU,char* JOBVT, long* M, long* N, double* APtr, long* LDA, double* SPtr, double* UPtr, long* LDU, double* VTPtr, long* LDVT,
                       double* WORKtmp, long* LWORK, long* INFO);

extern "C" void dgelsy_(long* M, long* N, long* NRHS, double* APtr,
                        long* LDA, double* BPtr, long* LDB, long* JPVT, double* RCOND,
                        long* RANK, double* WORK, long* LWORK, long* INFO);
#endif /* SCC_LAPACKHEADERS_H_ */



