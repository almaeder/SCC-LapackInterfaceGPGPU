/*
 * TriDiagRoutines.h
 *
 *  Created on: Jul 4, 2020
 *      Author: anderson
 */
/*
#############################################################################
#
# Copyright  2020- Chris Anderson
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
#include <stdexcept>
#include <exception>
#include <sstream>

#include "SCC_LapackMatrix.h"

#ifndef TRI_DIAG_ROUTINES_
#define TRI_DIAG_ROUTINES_

// For general tri-diagonal matrices,the diagonal component of an N x N s
// tri-diagonal matrix is specified with a double vector of size N.
//
// The off-diagonal components are specified with  double vectors of size N-1

// The diagonal component of an N x N symmetric tri-diagonal matrix is specified
// with a double vector of size N.
//
// The off-diagonal component is specified with a double vector of size N-1


namespace SCC
{
class TriDiagRoutines
{
public:

//
//##############################################################################
//                   GENERAL TRIDIAGONAL MATRIX ROUTINES
//##############################################################################
//
void realTriLUfactorization(std::vector<double>& DL, std::vector<double>& D, std::vector<double>& DU,
		                    std::vector<double>& DU2, std::vector<long>& IPIV)
{
	long N = (long)D.size();
	long INFO = 0;

	IPIV.clear();
	IPIV.resize(N,0);

	DU2.clear();
    DU2.resize(N-2,0.0);

	dgttrf_(&N, &DL[0],&D[0],&DU[0],&DU2[0], &IPIV[0],&INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nrealTriLUfactorization (dgttrf) \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }
}

std::vector<double> realTriLUsolve(std::vector<double>& DL, std::vector<double>& D, std::vector<double>& DU,
		                           std::vector<double>& DU2, std::vector<long>& IPIV,std::vector<double>& B)
{
	long N = (long)D.size();
	long INFO = 0;

    char TRANS = 'N';
    long NRHS = 1;
    long LDB  = N;

    std::vector<double> X(B);

	dgttrs_(&TRANS, &N, &NRHS, &DL[0], &D[0], &DU[0], &DU2[0], &IPIV[0], &X[0], &LDB, &INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nrealTriLUsolve (dgttrs) \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }

    return X;
}


//
//##############################################################################
//        SYMMETRIC TRIDIAGONAL MATRIX ROUTINES
//##############################################################################
//

std::vector<double> realSymTriEigenVectors(const std::vector<double> & D,
const std::vector<double> & E, SCC::LapackMatrix & Q)
{

// D is the matrix diagonal (size N)
// E is the matrix sub(or super) -diagonal (size N-1)
//
    char compz = 'I'; // eigenvalues and eigenvectors

    std::vector<double>  Dtmp = D;   // create duplicate of D and E
    std::vector<double>  Etmp = E;   //
    double* DPtr   = &Dtmp[0];
    double* EPtr   = &Etmp[0];

    long n  = (long)D.size();

    Q.initialize(n,n);
    Q.setToIdentity();

    double* QPtr    = Q.getDataPointer();

    long ldz         =   n;
    double* workPtr  = new double[2*n + 2];  // work array
    long info        = 0;

    dsteqr_(&compz, &n, DPtr, EPtr, QPtr, &ldz,workPtr, &info);
    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\nrealSymTriEigenVectors LAPACK (dsteqr) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }
//
//  clean up
//
    delete [] workPtr;
    return Dtmp;
};

std::vector<double> realSymTriEigenValues(const std::vector<double>& D,
const std::vector<double> E)
{
//
// The eigenvalues of a symmetric tridiagonal matrix are computed and returned.
//
// D is the matrix diagonal
// E is the matrix sub(or super) -diagonal.
//
    char compz = 'N'; // eigenvalues only

    std::vector<double>  Dtmp = D;   // create duplicate of D and E
    std::vector<double>  Etmp = E;
    double* DPtr   = &Dtmp[0];
    double* EPtr   = &Etmp[0];

    long n         = (long)D.size();

    double  ZPtr   =   0;
    long ldz       =   1;

    double* workPtr  = 0;   // work array not used in this call
    long info        = 0;

    dsteqr_(&compz, &n, DPtr, EPtr, &ZPtr, &ldz,workPtr, &info);
    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\nrealSymTriEigenValues LAPACK (dsteqr) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }

   return Dtmp;
};


// This routine returns the lowest nValues eigenvalues as a std::vector<double> and
// sets columns of Q to be the eigenvectors.

std::vector<double> getLowestSymTriEigSystem(long nValues, const std::vector<double>& D,
		const std::vector<double>& U, SCC::LapackMatrix& Q)
{
    long N = (long)D.size();

	if(nValues > N)
	{
		std::stringstream sout;
		sout << "\ngetLowestSymTriEigValues error \nError: Requested number of eigenvalues/eigenvectors > dimension";
		sout << "\nRequested number of eigenvalues = " << nValues;
		sout << "\nDimension of matrix             = " << N;
		throw std::runtime_error(sout.str());
	}

    std::vector<double> eVals(N);


    Q.initialize(N,nValues);

    std::vector<double> Dtmp = D;
    std::vector<double> Utmp = U;

//  Input parameters

    char jobz     =  'V';
    char range    =  'I';
    long n        =    N;
    double vLower =  0.0;
    double vUpper =  0.0;
    long   iLower =    1;       // lower computed eigenvalue index
    long   iUpper =    nValues; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &Dtmp[0];
    double* uPtr = &Utmp[0];
//
//  Output parameters
//
    long mFound = 0;  // number of eigenvalues found

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    double* vPtr   = Q.getDataPointer();    // array for the eigenvectors
    long   ldz     = N;

    double* work   = new double[5*N];   // work array
    long*  iwork   = new long[5*N];     // work array

    long*   ifail = new long[N];
    long    info  = 0;

    dstevx_(&jobz, &range, &n, dPtr,uPtr,&vLower, &vUpper, &iLower, &iUpper,
    &abstol, &mFound, ePtr, vPtr, &ldz, work, iwork, ifail, &info);

    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigSystem LAPACK (dstevx) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }

    /* extract return eigenvalues */

    std::vector<double>     eValsReturn;

    if(mFound ==  0) {return eValsReturn;}

    eValsReturn.resize(mFound);

    for(long i = 0; i < mFound; i++)
    {
    	eValsReturn[i] = eVals[i];
    }
    /* clean up */

    delete [] ifail;
    delete [] work;
    delete [] iwork;

    return eValsReturn;
}
///
/// The return value of this routine is the number M of eigenvalues found in
//  the interval (lambdaMin, lambdaMax]. The array eigVals contains the
//  eigenvalues and the first M columns of Q contain the eigenvectors.
///
long  getLowestSymTriEigSystem(double lambdaMin, double lambdaMax, const std::vector<double>& D,  const std::vector<double>& U, std::vector<double>& eigVals,
SCC::LapackMatrix& Q)
{
    long i;

    long N = (long)D.size();

    std::vector<double> eVals(N);

    Q.initialize(N,N);

    std::vector<double>  Dtmp(D);
    std::vector<double>  Utmp(U);
//
//  Input paramters
//
    char jobz     =  'V';
    char range    =  'V';
    long n        =    N;
    double vLower =  lambdaMin;
    double vUpper =  lambdaMax;
    long   iLower =    0;
    long   iUpper =    0;         // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &Dtmp[0];
    double* uPtr = &Utmp[0];
//
//  Output parameters
//
    long mFound = 0;  // number of eigenvalues found

    double* ePtr   = &eVals[0];             // array for the eigenvalues
    double* vPtr   = Q.getDataPointer();    // array for the eigenvectors
    long   ldz     = N;

    double* work   = new double[5*N];   // work array
    long*  iwork   = new long[5*N];     // work array

    long*   ifail = new long[N];
    long    info  = 0;


    dstevx_(&jobz, &range, &n, dPtr,uPtr,&vLower, &vUpper, &iLower, &iUpper,
    &abstol, &mFound, ePtr, vPtr, &ldz, work, iwork, ifail, &info);

    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigSystem LAPACK (dstevx) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }

    if(mFound == 0)
    {
    eigVals.clear();
    }
    else
    {
    eigVals.clear();
    eigVals.resize(mFound);
    for(i = 0; i < mFound; i++)
    {eigVals[i] = eVals[i];}
    }

    /* clean up */

    delete [] ifail;
    delete [] work;
    delete [] iwork;

    return mFound;
}



std::vector<double> getLowestSymTriEigValues(long nValues, std::vector<double> & D,  std::vector<double> & U)
{
    long i;

    long N = (long)D.size();

    std::vector<double>  eVals(N);

    if(nValues > N)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigValues error \nError: Requested number of eigenvalues > dimension";
    	sout << "\nRequested number of eigenvalues = " << nValues;
    	sout << "\nDimension of matrix             = " << N;
    	throw std::runtime_error(sout.str());
    }

//  Input paramters

    char range    =  'I';
    char order    =  'E';
    long n        =    N;
    double vLower =  0.0;
    double vUpper =  0.0;
    long   iLower =    1;       // lower computed eigenvalue index
    long   iUpper =    nValues; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &D[0];
    double* uPtr = &U[0];
//
//  Output parameters
//
    long mFound = 0;   // number of eigenvalues found
    long nsplit = 0;   // number of diagonal blocks

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    long*   iblock = new long[N];
    long*   isplit = new long[N];

    double* work   = new double[4*N];   // work array
    long*  iwork   = new long[3*N];     // work array

    long   info = 0;

    dstebz_(&range, &order, &n, &vLower, &vUpper, &iLower, &iUpper,
    &abstol, dPtr, uPtr, &mFound, &nsplit, ePtr, iblock, isplit, work, iwork, &info);

    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigValues LAPACK (dstebz) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }
    /* extract return eigenvalues */


    std::vector<double>  eValsReturn(nValues);

    for(i = 0; i < nValues; i++)
    {eValsReturn[i] = eVals[i];}

    /* clean up */

    delete [] iblock;
    delete [] isplit;
    delete [] work;
    delete [] iwork;

    return eValsReturn;
}


long  getLowestSymTriEigValues(double lowerLimit, double upperLimit,
std::vector<double>& D, std::vector<double>& U,std::vector<double>& eigVals)
{
    long i;

    long N = (long)D.size();

    std::vector<double>  eVals(N);
//
//  Input parameters
//
    char range    =  'V';
    char order    =  'E';
    long n        =   N;
    double vLower =  lowerLimit;
    double vUpper =  upperLimit;
    long   iLower =    1; // lower computed eigenvalue index
    long   iUpper =    1; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &D[0];
    double* uPtr = &U[0];
//
//  Output parameters
//
    long mFound = 0;   // number of eigenvalues found
    long nsplit = 1;   // number of diagonal blocks

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    long*   iblock = new long[N];
    long*   isplit = new long[N];

    double* work   = new double[4*N];   // work array
    long*  iwork   = new long[3*N];     // work array
    long   info    = 0;



    dstebz_(&range, &order, &n, &vLower, &vUpper, &iLower, &iUpper,
    &abstol, dPtr, uPtr, &mFound, &nsplit, ePtr, iblock, isplit, work, iwork,
    &info);

    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigValues LAPACK (dstebz) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }


    if(mFound == 0)
    {eigVals.clear();}
    else
    {
    eigVals.clear();
    eigVals.resize(mFound);
    for(i = 0; i < mFound; i++)
    {eigVals[i] = eVals[i];}
    }


    // clean up //

    delete [] iblock;
    delete [] isplit;
    delete [] work;
    delete [] iwork;


    return mFound;
}
};
}
#endif /* SCC_TriDiagRoutines__ */


