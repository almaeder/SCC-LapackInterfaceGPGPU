/*
 * TriDiagRoutines.h
 *
 *  Created on: Jul 4, 2020
 *      Author: anderson
 *
 *  Updated on: Dec. 11, 2023 (C.R. Anderson)
 */
//
// SCC::TriDiagRoutines
//
// A utility class for tri-diagonal matrices whose
// functionality is based upon LAPACK routines.
//
// The documentation for the each of the base LAPACK routines
// is contained at the end of this file or can be found at
//
// https://netlib.org/lapack/explore-html
//
// The member functions do not provide the complete functionality of the
// LAPACK routines upon which they are based -- only the
// functionality as needed for specific project use, functionality
// that may be updated without notice.
//
// Data mapping being used for direct invocation of
// Fortran routines
//
// C++  int    ==  Fortran LOGICAL
// C++  RC_INT   ==  Fortran INTEGER
// C++  double ==  Fortran DOUBLE PRECISION
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Current list of LAPACK routines being used by member functions
// of the TriDiagRoutines class.
//
// DGTTRF : computes an LU factorization of a real tridiagonal matrix A)
//
// DGTTRS : solves systems of equations using an LU factorization
//          (with pivoting) from DGTTRF
//
// DSTEQR : computes all eigenvalues and, optionally, eigenvectors of a
// symmetric tridiagonal matrix using the implicit QL or QR method.
//
// DSTEVX : computes selected eigenvalues and, optionally, eigenvectors
//          of a real symmetric tridiagonal matrix A.
//
// DSTEBZ : computes the eigenvalues of a symmetric tridiagonal
//          matrix T.  The user may ask for all eigenvalues, all eigenvalues
//          in the half-open interval (VL, VU], or the IL-th through IU-th
//          eigenvalues.
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
		                    std::vector<double>& DU2, std::vector<RC_INT>& IPIV)
{
	RC_INT N = (RC_INT)D.size();
	RC_INT INFO = 0;

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
		                           std::vector<double>& DU2, std::vector<RC_INT>& IPIV,std::vector<double>& B)
{
	RC_INT N = (RC_INT)D.size();
	RC_INT INFO = 0;

    char TRANS = 'N';
    RC_INT NRHS = 1;
    RC_INT LDB  = N;

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

    RC_INT n  = (RC_INT)D.size();

    Q.initialize(n,n);
    Q.setToIdentity();

    double* QPtr    = Q.getDataPointer();

    RC_INT ldz         =   n;
    double* workPtr  = new double[2*n + 2];  // work array
    RC_INT info        = 0;

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

    RC_INT n         = (RC_INT)D.size();

    double  ZPtr   =   0;
    RC_INT ldz       =   1;

    double* workPtr  = 0;   // work array not used in this call
    RC_INT info        = 0;

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

std::vector<double> getLowestSymTriEigSystem(RC_INT nValues, const std::vector<double>& D,
		const std::vector<double>& U, SCC::LapackMatrix& Q)
{
    RC_INT N = (RC_INT)D.size();

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
    RC_INT n        =    N;
    double vLower =  0.0;
    double vUpper =  0.0;
    RC_INT   iLower =    1;       // lower computed eigenvalue index
    RC_INT   iUpper =    nValues; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &Dtmp[0];
    double* uPtr = &Utmp[0];
//
//  Output parameters
//
    RC_INT mFound = 0;  // number of eigenvalues found

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    double* vPtr   = Q.getDataPointer();    // array for the eigenvectors
    RC_INT   ldz     = N;

    double* work   = new double[5*N];   // work array
    RC_INT*  iwork   = new RC_INT[5*N];     // work array

    RC_INT*   ifail = new RC_INT[N];
    RC_INT    info  = 0;

    dstevx_(&jobz, &range, &n, dPtr,uPtr,&vLower, &vUpper, &iLower, &iUpper,
    &abstol, &mFound, ePtr, vPtr, &ldz, work, iwork, ifail, &info);

    if(info != 0)
    {
    	std::stringstream sout;
    	sout << "\ngetLowestSymTriEigSystem LAPACK (dstevx) error \nError INFO = " << info << "\n";
    	throw std::runtime_error(sout.str());
    }

    /* extract return eigenvalues */

    std::vector<double> eValsReturn;

    if(mFound ==  0) {return eValsReturn;}

    eValsReturn.resize(mFound);

    for(RC_INT i = 0; i < mFound; i++)
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
RC_INT  getLowestSymTriEigSystem(double lambdaMin, double lambdaMax, const std::vector<double>& D,  const std::vector<double>& U, std::vector<double>& eigVals,
SCC::LapackMatrix& Q)
{
    RC_INT i;

    RC_INT N = (RC_INT)D.size();

    std::vector<double> eVals(N);

    Q.initialize(N,N);

    std::vector<double>  Dtmp(D);
    std::vector<double>  Utmp(U);
//
//  Input paramters
//
    char jobz     =  'V';
    char range    =  'V';
    RC_INT n        =    N;
    double vLower =  lambdaMin;
    double vUpper =  lambdaMax;
    RC_INT   iLower =    0;
    RC_INT   iUpper =    0;         // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &Dtmp[0];
    double* uPtr = &Utmp[0];
//
//  Output parameters
//
    RC_INT mFound = 0;  // number of eigenvalues found

    double* ePtr   = &eVals[0];             // array for the eigenvalues
    double* vPtr   = Q.getDataPointer();    // array for the eigenvectors
    RC_INT   ldz     = N;

    double* work   = new double[5*N];   // work array
    RC_INT*  iwork   = new RC_INT[5*N];     // work array

    RC_INT*   ifail = new RC_INT[N];
    RC_INT    info  = 0;


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



std::vector<double> getLowestSymTriEigValues(RC_INT nValues, std::vector<double> & D,  std::vector<double> & U)
{
    RC_INT i;

    RC_INT N = (RC_INT)D.size();

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
    RC_INT n        =    N;
    double vLower =  0.0;
    double vUpper =  0.0;
    RC_INT   iLower =    1;       // lower computed eigenvalue index
    RC_INT   iUpper =    nValues; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &D[0];
    double* uPtr = &U[0];
//
//  Output parameters
//
    RC_INT mFound = 0;   // number of eigenvalues found
    RC_INT nsplit = 0;   // number of diagonal blocks

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    RC_INT*   iblock = new RC_INT[N];
    RC_INT*   isplit = new RC_INT[N];

    double* work   = new double[4*N];   // work array
    RC_INT*  iwork   = new RC_INT[3*N];     // work array

    RC_INT   info = 0;

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


RC_INT  getLowestSymTriEigValues(double lowerLimit, double upperLimit,
std::vector<double>& D, std::vector<double>& U,std::vector<double>& eigVals)
{
    RC_INT i;

    RC_INT N = (RC_INT)D.size();

    std::vector<double>  eVals(N);
//
//  Input parameters
//
    char range    =  'V';
    char order    =  'E';
    RC_INT n        =   N;
    double vLower =  lowerLimit;
    double vUpper =  upperLimit;
    RC_INT   iLower =    1; // lower computed eigenvalue index
    RC_INT   iUpper =    1; // upper computed eigenvalue index

    double abstol = 1.0e-14;

    double* dPtr = &D[0];
    double* uPtr = &U[0];
//
//  Output parameters
//
    RC_INT mFound = 0;   // number of eigenvalues found
    RC_INT nsplit = 1;   // number of diagonal blocks

    double* ePtr   = &eVals[0]; // array for the eigenvalues
    RC_INT*   iblock = new RC_INT[N];
    RC_INT*   isplit = new RC_INT[N];

    double* work   = new double[4*N];   // work array
    RC_INT*  iwork   = new RC_INT[3*N];     // work array
    RC_INT   info    = 0;



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


////////////////////////////////////////////////////////////////
// DGTTRF
////////////////////////////////////////////////////////////////
/*

DGTTRF computes an LU factorization of a real tridiagonal matrix A

subroutine dgttrf	(	integer 	n,
double precision, dimension( * ) 	dl,
double precision, dimension( * ) 	d,
double precision, dimension( * ) 	du,
double precision, dimension( * ) 	du2,
integer, dimension( * ) 	ipiv,
integer 	info
)

Purpose:
 DGTTRF computes an LU factorization of a real tridiagonal matrix A
 using elimination with partial pivoting and row interchanges.

 The factorization has the form
    A = L * U
 where L is a product of permutation and unit lower bidiagonal
 matrices and U is upper triangular with nonzeros in only the main
 diagonal and first two superdiagonals.
Parameters
[in]	N
          N is INTEGER
          The order of the matrix A.
[in,out]	DL
          DL is DOUBLE PRECISION array, dimension (N-1)
          On entry, DL must contain the (n-1) sub-diagonal elements of
          A.

          On exit, DL is overwritten by the (n-1) multipliers that
          define the matrix L from the LU factorization of A.
[in,out]	D
          D is DOUBLE PRECISION array, dimension (N)
          On entry, D must contain the diagonal elements of A.

          On exit, D is overwritten by the n diagonal elements of the
          upper triangular matrix U from the LU factorization of A.
[in,out]	DU
          DU is DOUBLE PRECISION array, dimension (N-1)
          On entry, DU must contain the (n-1) super-diagonal elements
          of A.

          On exit, DU is overwritten by the (n-1) elements of the first
          super-diagonal of U.
[out]	DU2
          DU2 is DOUBLE PRECISION array, dimension (N-2)
          On exit, DU2 is overwritten by the (n-2) elements of the
          second super-diagonal of U.
[out]	IPIV
          IPIV is INTEGER array, dimension (N)
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required.
[out]	INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -k, the k-th argument had an illegal value
          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
                has been completed, but the factor U is exactly
                singular, and division by zero will occur if it is used
                to solve a system of equations.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

////////////////////////////////////////////////////////////////
// DGTTRS
////////////////////////////////////////////////////////////////
/*
DGTTRS solves systems of equations using an LU factorization
(with pivoting) from DGTTRF

dgttrs()
subroutine dgttrs	(	character 	trans,
integer 	n,
integer 	nrhs,
double precision, dimension( * ) 	dl,
double precision, dimension( * ) 	d,
double precision, dimension( * ) 	du,
double precision, dimension( * ) 	du2,
integer, dimension( * ) 	ipiv,
double precision, dimension( ldb, * ) 	b,
integer 	ldb,
integer 	info
)


Purpose:
 DGTTRS solves one of the systems of equations
    A*X = B  or  A**T*X = B,
 with a tridiagonal matrix A using the LU factorization computed
 by DGTTRF.
Parameters
[in]	TRANS
          TRANS is CHARACTER*1
          Specifies the form of the system of equations.
          = 'N':  A * X = B  (No transpose)
          = 'T':  A**T* X = B  (Transpose)
          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
[in]	N
          N is INTEGER
          The order of the matrix A.
[in]	NRHS
          NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.
[in]	DL
          DL is DOUBLE PRECISION array, dimension (N-1)
          The (n-1) multipliers that define the matrix L from the
          LU factorization of A.
[in]	D
          D is DOUBLE PRECISION array, dimension (N)
          The n diagonal elements of the upper triangular matrix U from
          the LU factorization of A.
[in]	DU
          DU is DOUBLE PRECISION array, dimension (N-1)
          The (n-1) elements of the first super-diagonal of U.
[in]	DU2
          DU2 is DOUBLE PRECISION array, dimension (N-2)
          The (n-2) elements of the second super-diagonal of U.
[in]	IPIV
          IPIV is INTEGER array, dimension (N)
          The pivot indices; for 1 <= i <= n, row i of the matrix was
          interchanged with row IPIV(i).  IPIV(i) will always be either
          i or i+1; IPIV(i) = i indicates a row interchange was not
          required.
[in,out]	B
          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the matrix of right hand side vectors B.
          On exit, B is overwritten by the solution vectors X.
[in]	LDB
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
[out]	INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
 */
////////////////////////////////////////////////////////////////
// DSTEQR
////////////////////////////////////////////////////////////////
/*
DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
symmetric tridiagonal matrix.

dsteqr()
subroutine dsteqr	(	character 	compz,
integer 	n,
double precision, dimension( * ) 	d,
double precision, dimension( * ) 	e,
double precision, dimension( ldz, * ) 	z,
integer 	ldz,
double precision, dimension( * ) 	work,
integer 	info
)


Purpose:
 DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
 symmetric tridiagonal matrix using the implicit QL or QR method.
 The eigenvectors of a full or band symmetric matrix can also be found
 if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
 tridiagonal form.
Parameters
[in]	COMPZ
          COMPZ is CHARACTER*1
          = 'N':  Compute eigenvalues only.
          = 'V':  Compute eigenvalues and eigenvectors of the original
                  symmetric matrix.  On entry, Z must contain the
                  orthogonal matrix used to reduce the original matrix
                  to tridiagonal form.
          = 'I':  Compute eigenvalues and eigenvectors of the
                  tridiagonal matrix.  Z is initialized to the identity
                  matrix.
[in]	N
          N is INTEGER
          The order of the matrix.  N >= 0.
[in,out]	D
          D is DOUBLE PRECISION array, dimension (N)
          On entry, the diagonal elements of the tridiagonal matrix.
          On exit, if INFO = 0, the eigenvalues in ascending order.
[in,out]	E
          E is DOUBLE PRECISION array, dimension (N-1)
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix.
          On exit, E has been destroyed.
[in,out]	Z
          Z is DOUBLE PRECISION array, dimension (LDZ, N)
          On entry, if  COMPZ = 'V', then Z contains the orthogonal
          matrix used in the reduction to tridiagonal form.
          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
          orthonormal eigenvectors of the original symmetric matrix,
          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
          of the symmetric tridiagonal matrix.
          If COMPZ = 'N', then Z is not referenced.
[in]	LDZ
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          eigenvectors are desired, then  LDZ >= max(1,N).
[out]	WORK
          WORK is DOUBLE PRECISION array, dimension (max(1,2*N-2))
          If COMPZ = 'N', then WORK is not referenced.
[out]	INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  the algorithm has failed to find all the eigenvalues in
                a total of 30*N iterations; if INFO = i, then i
                elements of E have not converged to zero; on exit, D
                and E contain the elements of a symmetric tridiagonal
                matrix which is orthogonally similar to the original
                matrix.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
 */
////////////////////////////////////////////////////////////////
// DSTEVX
////////////////////////////////////////////////////////////////
/*
DSTEVX computes selected eigenvalues and, optionally, eigenvectors
of a real symmetric tridiagonal matrix A.


subroutine dstevx	(	character 	jobz,
character 	range,
integer 	n,
double precision, dimension( * ) 	d,
double precision, dimension( * ) 	e,
double precision 	vl,
double precision 	vu,
integer 	il,
integer 	iu,
double precision 	abstol,
integer 	m,
double precision, dimension( * ) 	w,
double precision, dimension( ldz, * ) 	z,
integer 	ldz,
double precision, dimension( * ) 	work,
integer, dimension( * ) 	iwork,
integer, dimension( * ) 	ifail,
integer 	info
)

Purpose:
 DSTEVX computes selected eigenvalues and, optionally, eigenvectors
 of a real symmetric tridiagonal matrix A.  Eigenvalues and
 eigenvectors can be selected by specifying either a range of values
 or a range of indices for the desired eigenvalues.
Parameters
[in]	JOBZ
          JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
[in]	RANGE
          RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found.
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found.
          = 'I': the IL-th through IU-th eigenvalues will be found.
[in]	N
          N is INTEGER
          The order of the matrix.  N >= 0.
[in,out]	D
          D is DOUBLE PRECISION array, dimension (N)
          On entry, the n diagonal elements of the tridiagonal matrix
          A.
          On exit, D may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues.
[in,out]	E
          E is DOUBLE PRECISION array, dimension (max(1,N-1))
          On entry, the (n-1) subdiagonal elements of the tridiagonal
          matrix A in elements 1 to N-1 of E.
          On exit, E may be multiplied by a constant factor chosen
          to avoid over/underflow in computing the eigenvalues.
[in]	VL
          VL is DOUBLE PRECISION
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]	VU
          VU is DOUBLE PRECISION
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]	IL
          IL is INTEGER
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]	IU
          IU is INTEGER
          If RANGE='I', the index of the
          largest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]	ABSTOL
          ABSTOL is DOUBLE PRECISION
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to

                  ABSTOL + EPS *   max( |a|,|b| ) ,

          where EPS is the machine precision.  If ABSTOL is less
          than or equal to zero, then  EPS*|T|  will be used in
          its place, where |T| is the 1-norm of the tridiagonal
          matrix.

          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*DLAMCH('S').

          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3.
[out]	M
          M is INTEGER
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
[out]	W
          W is DOUBLE PRECISION array, dimension (N)
          The first M elements contain the selected eigenvalues in
          ascending order.
[out]	Z
          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M) )
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge (INFO > 0), then that
          column of Z contains the latest approximation to the
          eigenvector, and the index of the eigenvector is returned
          in IFAIL.  If JOBZ = 'N', then Z is not referenced.
          Note: the user must ensure that at least max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used.
[in]	LDZ
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).
[out]	WORK
          WORK is DOUBLE PRECISION array, dimension (5*N)
[out]	IWORK
          IWORK is INTEGER array, dimension (5*N)
[out]	IFAIL
          IFAIL is INTEGER array, dimension (N)
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge.
          If JOBZ = 'N', then IFAIL is not referenced.
[out]	INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, then i eigenvectors failed to converge.
                Their indices are stored in array IFAIL.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

////////////////////////////////////////////////////////////////
// DSTEBZ
////////////////////////////////////////////////////////////////
/*
DSTEBZ computes the eigenvalues of a symmetric tridiagonal
matrix T.  The user may ask for all eigenvalues, all eigenvalues
in the half-open interval (VL, VU], or the IL-th through IU-th
eigenvalues.

subroutine dstebz	(	character 	range,
character 	order,
integer 	n,
double precision 	vl,
double precision 	vu,
integer 	il,
integer 	iu,
double precision 	abstol,
double precision, dimension( * ) 	d,
double precision, dimension( * ) 	e,
integer 	m,
integer 	nsplit,
double precision, dimension( * ) 	w,
integer, dimension( * ) 	iblock,
integer, dimension( * ) 	isplit,
double precision, dimension( * ) 	work,
integer, dimension( * ) 	iwork,
integer 	info
)


Purpose:
 DSTEBZ computes the eigenvalues of a symmetric tridiagonal
 matrix T.  The user may ask for all eigenvalues, all eigenvalues
 in the half-open interval (VL, VU], or the IL-th through IU-th
 eigenvalues.

 To avoid overflow, the matrix must be scaled so that its
 largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
 accuracy, it should not be much smaller than that.

 See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
 Matrix", Report CS41, Computer Science Dept., Stanford
 University, July 21, 1966.
Parameters
[in]	RANGE
          RANGE is CHARACTER*1
          = 'A': ("All")   all eigenvalues will be found.
          = 'V': ("Value") all eigenvalues in the half-open interval
                           (VL, VU] will be found.
          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
                           entire matrix) will be found.
[in]	ORDER
          ORDER is CHARACTER*1
          = 'B': ("By Block") the eigenvalues will be grouped by
                              split-off block (see IBLOCK, ISPLIT) and
                              ordered from smallest to largest within
                              the block.
          = 'E': ("Entire matrix")
                              the eigenvalues for the entire matrix
                              will be ordered from smallest to
                              largest.
[in]	N
          N is INTEGER
          The order of the tridiagonal matrix T.  N >= 0.
[in]	VL
          VL is DOUBLE PRECISION

          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be returned.  VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]	VU
          VU is DOUBLE PRECISION

          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues.  Eigenvalues less than or equal
          to VL, or greater than VU, will not be returned.  VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]	IL
          IL is INTEGER

          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]	IU
          IU is INTEGER

          If RANGE='I', the index of the
          largest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]	ABSTOL
          ABSTOL is DOUBLE PRECISION
          The absolute tolerance for the eigenvalues.  An eigenvalue
          (or cluster) is considered to be located if it has been
          determined to lie in an interval whose width is ABSTOL or
          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
          will be used, where |T| means the 1-norm of T.

          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
[in]	D
          D is DOUBLE PRECISION array, dimension (N)
          The n diagonal elements of the tridiagonal matrix T.
[in]	E
          E is DOUBLE PRECISION array, dimension (N-1)
          The (n-1) off-diagonal elements of the tridiagonal matrix T.
[out]	M
          M is INTEGER
          The actual number of eigenvalues found. 0 <= M <= N.
          (See also the description of INFO=2,3.)
[out]	NSPLIT
          NSPLIT is INTEGER
          The number of diagonal blocks in the matrix T.
          1 <= NSPLIT <= N.
[out]	W
          W is DOUBLE PRECISION array, dimension (N)
          On exit, the first M elements of W will contain the
          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
          workspace.)
[out]	IBLOCK
          IBLOCK is INTEGER array, dimension (N)
          At each row/column j where E(j) is zero or small, the
          matrix T is considered to split into a block diagonal
          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
          block (from 1 to the number of blocks) the eigenvalue W(i)
          belongs.  (DSTEBZ may use the remaining N-M elements as
          workspace.)
[out]	ISPLIT
          ISPLIT is INTEGER array, dimension (N)
          The splitting points, at which T breaks up into submatrices.
          The first submatrix consists of rows/columns 1 to ISPLIT(1),
          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
          etc., and the NSPLIT-th consists of rows/columns
          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
          (Only the first NSPLIT elements will actually be used, but
          since the user cannot know a priori what value NSPLIT will
          have, N words must be reserved for ISPLIT.)
[out]	WORK
          WORK is DOUBLE PRECISION array, dimension (4*N)
[out]	IWORK
          IWORK is INTEGER array, dimension (3*N)
[out]	INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  some or all of the eigenvalues failed to converge or
                were not computed:
                =1 or 3: Bisection failed to converge for some
                        eigenvalues; these eigenvalues are flagged by a
                        negative block number.  The effect is that the
                        eigenvalues may not be as accurate as the
                        absolute and relative tolerances.  This is
                        generally caused by unexpectedly inaccurate
                        arithmetic.
                =2 or 3: RANGE='I' only: Not all of the eigenvalues
                        IL:IU were found.
                        Effect: M < IU+1-IL
                        Cause:  non-monotonic arithmetic, causing the
                                Sturm sequence to be non-monotonic.
                        Cure:   recalculate, using RANGE='A', and pick
                                out eigenvalues IL:IU.  In some cases,
                                increasing the PARAMETER "FUDGE" may
                                make things work.
                = 4:    RANGE='I', and the Gershgorin interval
                        initially used was too small.  No eigenvalues
                        were computed.
                        Probable cause: your machine has sloppy
                                        floating-point arithmetic.
                        Cure: Increase the PARAMETER "FUDGE",
                              recompile, and try again.
Internal Parameters:
  RELFAC  DOUBLE PRECISION, default = 2.0e0
          The relative tolerance.  An interval (a,b] lies within
          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
          where "ulp" is the machine precision (distance from 1 to
          the next larger floating point number.)

  FUDGE   DOUBLE PRECISION, default = 2
          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
          a value of 1 should work, but on machines with sloppy
          arithmetic, this needs to be larger.  The default for
          publicly released versions should be large enough to handle
          the worst machine around.  Note that this has no effect
          on accuracy of the solution.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.

 */
#endif /* SCC_TriDiagRoutines__ */


