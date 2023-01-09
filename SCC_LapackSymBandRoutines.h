/*
 * SCC_LapackSymBandRoutines.h
 *
 *  Created on: Nov 29, 2022
 *      Author: anderson
 */
#include <vector>
#include <iostream>
#include <cstring>

#include "SCC_LapackSymBandMatrix.h"

#ifndef LAPACK_SYM_BAND_ROUTINES_
#define LAPACK_SYM_BAND_ROUTINES_

// Headers for symmetric band matrix routines
//
// Note: SCC_LapackSymBandMatrix stores the upper diagonal
// data of a symmetric matrix so the UPLO parameter in the
// Lapack symmetric routines should be specified as "U".
//

extern "C" void dsbmv_(char* UPLO, long* N, long* K, double* alpha, double* Aptr,
long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

namespace SCC
{
// Returns y = A*x

void applySymmetricBandMatrix(SCC::LapackSymBandMatrix& A,
std::vector<double>& x, std::vector<double>& y)
{
	long N = A.getDimension();

	if((long)y.size() != N)  {y.resize(N,0.0);}

	char UPLO    = 'U';
    long K       = A.getSuperDiagonalCount();
    double ALPHA = 1.0;
    double BETA  = 0.0;
    double*Aptr  = A.getDataPointer();
    double*Xptr  = &x[0];
    double*Yptr  = &y[0];
    long LDA     = K+1;
    long INCX    = 1;
    long INCY    = 1;

    dsbmv_(&UPLO,&N,&K,&ALPHA,Aptr,&LDA,Xptr,&INCX,&BETA,Yptr,&INCY);
}


}
#endif /* LAPACK_SYM_BAND_ROUTINES_ */
