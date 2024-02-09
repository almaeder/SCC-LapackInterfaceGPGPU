#include "SCC_LapackMatrix.h"
#include "SCC_LapackHeaders.h"
//
// SCC::LapackMatrixRoutines
//
// A collection of classes whose functionality is
// based upon LAPACK routines. These routine are meant
// to be used with instances of LapackMatrix. The documentation
// for the each of the base LAPACK routines is contained at the
// end of this file or can be found at
//
// https://netlib.org/lapack/explore-html
//
// These classes do not provide the complete functionality of the
// LAPACK routines upon which they are based -- only the 
// functionality as needed for specific project use, functionality
// that may be updated without notice.
//
// Author          : Chris Anderson
// Version         : Dec. 09, 2023
//
// Data mapping being used for direct invocation of
// Fortran routines
//
// C++  int    ==  Fortran LOGICAL
// C++  long   ==  Fortran INTEGER
// C++  double ==  Fortran DOUBLE PRECISION
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Current class list
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DGELSY : Created for least squares solution of A*X = B
// LAPACK base routine description:
//
// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//     minimize || A * X - B ||
// using a complete orthogonal factorization of A.  A is an M-by-N
// matrix which may be rank-deficient.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DGESVD : Created for singular value decomposition
// LAPACK base routine description:
//
// DGESVD computes the singular value decomposition (SVD) of a real
// M-by-N matrix A, optionally computing the left and/or right singular
//  vectors. The SVD is written
//
//       A = U * SIGMA * transpose(V)
//
//  where SIGMA is an M-by-N matrix which is zero except for its
//  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
//  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
//  are the singular values of A; they are real and non-negative, and
//  are returned in descending order.  The first min(m,n) columns of
//  U and V are the left and right singular vectors of A.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DSYEV : Created for eigensystem of symmetric matrix
// LAPACK base routine description:
// DSYEV computes all eigenvalues and, optionally, eigenvectors of a
// real symmetric matrix A.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DGESVX : Created for solving general  systems of equations
// LAPACK base routine description:
// DGESVX uses the LU factorization to compute the solution to a real
// system of linear equations
//      A * X = B,
// where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
// Error bounds on the solution and a condition estimate are also
// provided.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DPOSV : Created for solving symmetric positive definite systems
// of equations
// LAPACK base routine description:
// DPOSV computes the solution to a real system of linear equations
//    A * X = B,
// where A is an N-by-N symmetric positive definite matrix and X and B
// are N-by-NRHS matrices.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class NORMALEQ creates the solution of the normal equations using
// the singular value decomposition with singular value parameter cut-off
// value svdCutoff.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class QRutility created to construct solutions to
//               A x = b
// where A is a full rank M x N matrix with M >= N.
// If M > N then a least squares approximation to the
// solution is determined. QR factors are cached so
// so that solutions for multiple right hand sides
// can be computed efficiently.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class DSYEVX : Created for computing eigensystem of symmetric matricies
// LAPACK base routine description:
// DSYEVX computes selected eigenvalues and, optionally, eigenvectors
// of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
// selected by specifying either a range of values or a range of indices
// for the desired eigenvalues.
/*
#############################################################################
#
# Copyright  2015-2023 Chris Anderson
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


#include <iostream>
#include <vector>
#include <cassert>
#include <complex>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <stdexcept>

#ifndef  SCC_LAPACK_MATRIX_ROUTINES_
#define  SCC_LAPACK_MATRIX_ROUTINES_

namespace SCC
{

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DGELSY : Created for least squares solution of A*X = B
// LAPACK base routine description:
//
// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//     minimize || A * X - B ||
// using a complete orthogonal factorization of A.  A is an M-by-N
// matrix which may be rank-deficient.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class DGELSY
{
public:

    DGELSY()
    {
        initialize();
    }

    void initialize()
    {
    this->A.initialize();
    this->X.clear();
    this->WORK.clear();
    this->JPVT.clear();
    RCOND = 10.0*numLimits.epsilon();
    RANK  = 0;
    this->overwriteExtDataFlag = false;
    }


    void initialize(const DGELSY& dgelsy)
    {
    this->A    = dgelsy.A;
    this->X    = dgelsy.X;
    this->JPVT = dgelsy.JPVT;
    RCOND      = dgelsy.RCOND;
    RANK       = dgelsy.RANK;
    this->overwriteExtDataFlag = dgelsy.overwriteExtDataFlag;
    }

    long getRank()
    {
    return RANK;
    }

    std::vector<double> qrSolve(const std::vector<double>& B, const LapackMatrix& A, double rcondCutoff = -1.0)
    {
        if(rcondCutoff < 0) {RCOND = 10.0*numLimits.epsilon();}
        else                {RCOND = rcondCutoff;}

        if(not this->overwriteExtDataFlag)
        {
        // Copy A, since A is overwritten or destroyed by the routine.
            if(not equalMatrixDimensions(this->A,A))
            {
                this->A.initialize(A);
            }
            else
            {
                this->A = A;
            }
        }


        long M  = this->A.getRowDimension();
        long N  = this->A.getColDimension();

        // Capture B and increase size if N > M

        X = B;
        if(M < N) {X.resize(N,0.0);}

        JPVT.resize(N);

        long NRHS = 1;
        long LDA  = M;
        long LDB  = (long)X.size();


       // Query to obtain optimal work array size

        long LWORK = -1;
        INFO       =  0;

        double WORKtmp;

        dgelsy_(&M, &N, &NRHS, this->A.getDataPointer(), &LDA,&X[0],&LDB,
        &JPVT[0], &RCOND, &RANK, &WORKtmp, &LWORK,&INFO);


        LWORK = (long)(WORKtmp + 100);
        WORK.resize(LWORK);


        // Second call to create qr solution

        INFO       = 0;

        dgelsy_(&M, &N, &NRHS, this->A.getDataPointer(), &LDA,&X[0],&LDB,
        &JPVT[0], &RCOND, &RANK, &WORK[0], &LWORK,&INFO);

        if(INFO != 0)
        {
        std::cerr << "DGELSY  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }

        // Set X to be the right dimension
        X.resize(N);

        return X;
    }

    //
    // A convenience solve interface. It is assumed that Bptr points to contiguous
    // data of size of the number of rows of A. No bounds checking is performed.
    //

    std::vector<double> qrSolve(double* Bptr, const LapackMatrix& Ain, double rcondCutoff = -1.0)
    {
        long M     = Ain.getRowDimension();
        std::vector<double>                 B(M);
        std::memcpy(B.data(),Bptr,M*sizeof(double));
        return qrSolve(B,Ain,rcondCutoff);
    }

//
//  In this call, the matrix is passed via a pointer to the
//  matrix data assumed to be stored with the Fortran convention
//  by columns.
//
//  !!! Important: the input matrix data is overwritten during the
//  solution process. There is no bounds checking performed on
//  the input matrix.
//
//  This routine is added to avoid the need for extraneous copying
//  of input matrix data.
//
    std::vector<double> qrSolve(const std::vector<double>& B, long Arows, long Acols, double* Adata, double rcondCutoff = -1.0)
    {
        this->overwriteExtDataFlag = true;
        this->A.initialize(Arows,Acols,Adata);
        std::vector<double> X = qrSolve(B,this->A,rcondCutoff);
        this->overwriteExtDataFlag = false;
        return X;
    }

    bool equalMatrixDimensions(const LapackMatrix& A, const LapackMatrix& B)
    {
        if((A.rows != B.rows)||(A.cols != B.cols)) return false;
        return true;
    }

    LapackMatrix           A;
    std::vector<double>    X;
    std::vector<long>   JPVT;
    double             RCOND;
    long                RANK;

    long                INFO;
    std::vector<double> WORK;

    std::numeric_limits<double>   numLimits;

    bool overwriteExtDataFlag;
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DGESVD : Created for singular value decomposition
// LAPACK base routine description:
//
// DGESVD computes the singular value decomposition (SVD) of a real
// M-by-N matrix A, optionally computing the left and/or right singular
//  vectors. The SVD is written
//
//       A = U * SIGMA * transpose(V)
//
//  where SIGMA is an M-by-N matrix which is zero except for its
//  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
//  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
//  are the singular values of A; they are real and non-negative, and
//  are returned in descending order.  The first min(m,n) columns of
//  U and V are the left and right singular vectors of A.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class DGESVD
{
public:

    DGESVD()
    {
        initialize();
    }


    DGESVD(const DGESVD& dgesvd)
    {
        initialize(dgesvd);
    }

    void initialize()
    {
    this->A.initialize();
    this->U.initialize();

    this->singularValues.clear();
    this->WORK.clear();

    this->VT.initialize();
    this->svdDim = 0;
    this->overwriteExtDataFlag = false;
    }

    void initialize(const DGESVD& dgesvd)
    {
    this->A = dgesvd.A;
    this->U = dgesvd.U;

    this->singularValues       = dgesvd.singularValues;
    this->VT                   = dgesvd.VT;
    this->svdDim               = dgesvd.svdDim;
    this->overwriteExtDataFlag = dgesvd.overwriteExtDataFlag;
    }


    bool equalMatrixDimensions(const SCC::LapackMatrix& A, const SCC::LapackMatrix& B)
    {
        if((A.rows != B.rows)||(A.cols != B.cols)) return false;
        return true;
    }

    void computeSVD(const SCC::LapackMatrix& A)
    {

        if(not this->overwriteExtDataFlag)
        {
            // Copy A, since A is overwritten or destroyed by the routine.

            if(not equalMatrixDimensions(this->A,A))
            {
                this->A.initialize(A);
            }
            else
            {
                this->A = A;
            }
        }

        long M  = this->A.getRowDimension();
        long N  = this->A.getColDimension();

        // S dimension (min(M,N))

        long minMN = (M < N) ? M : N;

        if((U.rows != M) || (U.cols != M)) {U.initialize(M,M);}

        singularValues.resize(minMN,0.0);

        if((VT.rows != N) || (VT.cols != N)) {VT.initialize(N,N);}


        long  LDA = M;
        long  LDU = M;
        long LDVT = N;

        JOBU  = 'A';  //  All M columns of U are returned in array U:
        JOBVT = 'A';  //  All N rows of V**T are returned in the array VT;

       // Query to obtain optimal work array size

        long LWORK = -1;
        INFO       = 0;

        double WORKtmp;
        dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0],
                U.dataPtr, &LDU, VT.dataPtr, &LDVT, &WORKtmp, &LWORK, &INFO);

        LWORK = (long)(WORKtmp + 100);
        WORK.resize(LWORK);

        double* WORKptr =    &WORK[0];

        // Second call to create svd

        INFO       = 0;
        dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr,
                &LDU, VT.dataPtr, &LDVT, WORKptr, &LWORK, &INFO);

        if(INFO != 0)
        {
        std::cerr << "DGESVD  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }
    }

//
//  In this call, the matrix is passed via a pointer to the
//  matrix data assumed to be stored with the Fortran convention
//  by columns.
//
//  !!! Important: the input matrix data is overwritten during the
//  solution process. There is no bounds checking performed on
//  the input matrix.
//
//  This routine is added to avoid the need for extraneous copying
//  of input matrix data.
//
    void computeSVD(long Arows, long Acols, double* Adata)
    {
        this->overwriteExtDataFlag = true;
        this->A.initialize(Arows,Acols,Adata);
        computeSVD(this->A);
        this->overwriteExtDataFlag = false;
    }

    void computeThinSVD(const SCC::LapackMatrix& A)
    {
        if(not this->overwriteExtDataFlag)
        {
            // Copy A, so that input A is not overwritten or destroyed by the routine.

            if(not equalMatrixDimensions(this->A,A))
            {
                this->A.initialize(A);
            }
            else
            {
                this->A = A;
            }
        }


        long M  = this->A.getRowDimension();
        long N  = this->A.getColDimension();

        // S dimension (min(M,N))

        long minMN = (M < N) ? M : N;

        if((U.rows != M) || (U.cols != minMN)) {U.initialize(M,minMN);}

        singularValues.resize(minMN,0.0);

        if((VT.rows != minMN) || (VT.cols != N)) {VT.initialize(minMN,N);}


        long  LDA = M;
        long  LDU = M;
        long LDVT = minMN;

        JOBU  = 'S';  //  minMN columns of U are returned in array U:
        JOBVT = 'S';  //  minMN rows of V**T are returned in the array VT;

       // Query to obtain optimal work array size

        long LWORK = -1;
        INFO       = 0;

        double WORKtmp;
        dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               &WORKtmp, &LWORK, &INFO);

        LWORK = (long)(WORKtmp + 100);
        WORK.resize(LWORK);

        double* WORKptr =    &WORK[0];

        // Second call to create svd

        INFO       = 0;

        dgesvd_(&JOBU, &JOBVT, &M, &N, this->A.dataPtr, &LDA, &singularValues[0], U.dataPtr, &LDU, VT.dataPtr, &LDVT,
               WORKptr, &LWORK, &INFO);

        if(INFO != 0)
        {
        std::cerr << "THIN DGESVD  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }
    }

//
//  In this call, the matrix is passed via a pointer to the
//  matrix data assumed to be stored with the Fortran convention
//  by columns.
//
//  !!! Important: the input matrix data is overwritten during the
//  solution process. There is no bounds checking performed on
//  the input matrix.
//
//  This routine is added to avoid the need for extraneous copying
//  of input matrix data.
//
    void computeThinSVD(long Arows, long Acols, double* Adata)
    {
        this->overwriteExtDataFlag = true;
        this->A.initialize(Arows,Acols,Adata);
        computeThinSVD(this->A);
        this->overwriteExtDataFlag = false;
    }

    long getSVDdim(){return svdDim;}

    std::vector<double> applyPseudoInverse(std::vector<double>& b, double svdCutoff = -1.0)
    {

        std::vector<double> x(A.cols,0.0);

        svdDim = 0;

        if(svdCutoff < 0){svdDim = (long)singularValues.size();}
        else
        {
        for(long i = 0; i < (long)singularValues.size(); i++)
        {
            if(singularValues[i] > svdCutoff){svdDim++;}
        }
        }

        if(svdDim == 0) return x;

        //
        // Construct pseudo-inverse using components of SVD
        //

        std::vector<double> xStar;

        /*
        LapackMatrix   Ustar = U.getColSlice(0,svdDim-1);
        xStar = Ustar.applyTranspose(b);
        */

        char TRANS     = 'T';
        long Mstar     = U.rows;
        long Nstar     = svdDim;
        double ALPHA   = 1.0;
        double BETA    = 0.0;
        long LDA       = Mstar;
        long INCX      = 1;
        long INCY      = 1;

        xStar.resize(svdDim,0.0);

        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,U.dataPtr,&LDA,&b[0],&INCX,&BETA,&xStar[0],&INCY);

        for(long i = 0; i < svdDim; i++)
        {
            xStar[i] /= singularValues[i];
        }

        SCC::LapackMatrix Vstar;

        /*
        Vstar = VT.getRowSlice(0,svdDim-1);
        x = Vstar.applyTranspose(xStar);
        */

        if(svdDim == VT.rows)
        {
        TRANS     = 'T';
        Mstar     = VT.rows;
        Nstar     = VT.cols;
        ALPHA   = 1.0;
        BETA    = 0.0;
        LDA       = Mstar;
        INCX      = 1;
        INCY      = 1;
        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,VT.dataPtr,&LDA,&xStar[0],&INCX,&BETA,&x[0],&INCY);
        }
        else // Extract columns of Vstar as first svdDim rows of VT
        {
        Vstar.initialize(VT.cols,svdDim);
        for(long i = 0; i < VT.cols; i++)
        {
        for(long j = 0; j < svdDim;  j++)
        {
        Vstar(i,j) = VT(j,i);
        }}

        TRANS     = 'N';
        Mstar     = Vstar.rows;
        Nstar     = Vstar.cols;
        ALPHA   = 1.0;
        BETA    = 0.0;
        LDA       = Mstar;
        INCX      = 1;
        INCY      = 1;
        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,Vstar.dataPtr,&LDA,&xStar[0],&INCX,&BETA,&x[0],&INCX);
        }

        return x;
    }

    SCC::LapackMatrix     A;
    SCC::LapackMatrix     U;

    std::vector<double> singularValues;

    SCC::LapackMatrix     VT;

    char                 JOBU;
    char                JOBVT;
    long                 INFO;
    std::vector<double>  WORK;

    long               svdDim;
    bool overwriteExtDataFlag;
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DSYEV : eigensystem computation of real symmetric matrix.
// LAPACK base routine description:
// DSYEV computes all eigenvalues and, optionally, eigenvectors of a
// real symmetric matrix A.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class DSYEV
{
public:

    DSYEV()
    {
        JOBZ = 'N';
        UPLO = 'U';            // Using upper triangular part of A
    }

    void initialize()
    {
        U.initialize();
        eigValues.clear();
        JOBZ = 'N';
        UPLO = 'U';            // Using upper triangular part of A
    }
    void computeEigenvalues(const LapackMatrix& A, std::vector<double>& eigenValues)
    {
        assert(A.sizeCheck(A.rows,A.cols));

        U.initialize(A);
        JOBZ = 'N';
        UPLO = 'U';            // Using upper triangular part of A

        long N       = A.rows;
        double* Uptr = U.dataPtr;

        long LDA = N;

        eigenValues.resize(N);
        double*Wptr = &eigenValues[0];

        long LWORK = -1;

        double WORKtmp;

        long INFO = 0;

        // First call to get optimal workspace

        dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, &WORKtmp, &LWORK, &INFO);

        LWORK = (long)(WORKtmp + 100);
        std::vector<double>    WORK(LWORK);
        double* WORKptr =    &WORK[0];

        // Second call to create eigensystem

        dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, WORKptr, &LWORK, &INFO);

        if(INFO != 0)
        {
        std::cerr << "dsyev  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }

    }

    void computeEigensystem(const LapackMatrix& A, std::vector<double>& eigenValues, std::vector < std::vector < double> >& eigenVectors)
    {
        assert(A.sizeCheck(A.rows,A.cols));
        computeEigensystem(A,eigenValues, U);

        // Pack eigenvectors into return argument

        std::vector <double> eigVector(A.rows);

        eigenVectors.clear();
        eigenVectors.resize(A.cols,eigVector);

        for(long j = 0; j < A.cols; j++)
        {
            for(long i = 0; i < A.rows; i++)
            {
                eigenVectors[j][i] = U(i,j);
            }
        }

    }

   
    void computeEigensystem(const LapackMatrix& A, std::vector<double>& eigenValues, LapackMatrix& eigenVectors)
    {
        assert(A.sizeCheck(A.rows,A.cols));

        eigenVectors.initialize(A);
        JOBZ = 'V';
        UPLO = 'U';            // Using upper triangular part of A

        long N       = A.rows;
        double* Uptr = eigenVectors.dataPtr;

        long LDA = N;

        eigenValues.resize(N);
        double*Wptr = &eigenValues[0];

        long LWORK = -1;

        double WORKtmp;

        long INFO = 0;

        // First call to get optimal workspace

        dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, &WORKtmp, &LWORK, &INFO);

        LWORK = (long)(WORKtmp + 100);
        std::vector<double>    WORK(LWORK);
        double* WORKptr =    &WORK[0];

        // Second call to create eigensystem

        dsyev_(&JOBZ,&UPLO,&N, Uptr, &LDA, Wptr, WORKptr, &LWORK, &INFO);

        if(INFO != 0)
        {
        std::cerr << "dsyev  Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }
    }

    LapackMatrix                 U;
    std::vector<double>  eigValues;
    char            JOBZ;
    char            UPLO;

};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DGESVX : Created for solving general  systems of equations
// LAPACK base routine description:
// DGESVX uses the LU factorization to compute the solution to a real
// system of linear equations
//      A * X = B,
// where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
// Error bounds on the solution and a condition estimate are also
// provided.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class DGESVX
{
public:

    DGESVX()
    {
    initialize();
    }

    void initialize()
    {
    RCOND = 0;
    FERR.clear();
    BERR.clear();
    A.initialize();
    AF.initialize();
    }


    void applyInverse(const LapackMatrix& A,std::vector <double >& b)
    {
            applyInverse(A,&b[0]);
    }

    void applyInverse(const LapackMatrix& A,LapackMatrix& b)
    {
            applyInverse(A,b.dataPtr,b.cols);
    }

    void applyInverse(const LapackMatrix& A, double* b, long NRHS = 1)
    {
        assert(A.sizeCheck(A.rows,A.cols));

        char FACT  = 'E'; // Equilibrate, then factor
        char TRANS = 'N'; // No transpose
        long N     = A.rows;

        // Allocate temporaries

        this->A.initialize(A);
        this->AF.initialize(N,N);

        double* Aptr  =  A.dataPtr;
        double* AFptr = AF.dataPtr;

        long LDA   = N;
        long LDAF  = N;

        std::vector <long >   IPIV(N);
        long* IPIVptr = &IPIV[0];

        char  EQED;

        std::vector<double>   R(N);
        double* Rptr  = &R[0];

        std::vector<double>    C(N);
        double* Cptr  =  &C[0];

        std::vector<double>   B(N*NRHS);
        double* Bptr  =       &B[0];
        long LDB      =          N;


        // b will be overwritten with the solution
        // so no need to declare X separately

        double* Xptr = b;
        long LDX     = N;

        FERR.resize(NRHS);
        BERR.resize(NRHS);

        std::vector<double>   WORK(4*N);
        double* WORKptr = &WORK[0];

        std::vector<long>       IWORK(N);
        long* IWORKptr  = &IWORK[0];

        long   INFO = 0;


        // Assign right hand side to B

        for(long i = 0; i < N*NRHS; i++)
        {
            Bptr[i] = b[i];
        }

        dgesvx_(&FACT, &TRANS, &N, &NRHS, Aptr, &LDA, AFptr, &LDAF, IPIVptr,
                &EQED, Rptr, Cptr, Bptr,&LDB, Xptr, &LDX, &RCOND,
                &FERR[0], &BERR[0], WORKptr, IWORKptr, &INFO);


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
*  FERR is DOUBLE PRECISION array, dimension (NRHS)
*  The estimated forward error bound for each solution vector
*  X(j) (the j-th column of the solution matrix X).
*  If XTRUE is the true solution corresponding to X(j), FERR(j)
*  is an estimated upper bound for the magnitude of the largest
*  element in (X(j) - XTRUE) divided by the magnitude of the
*  largest element in X(j).  The estimate is as reliable as
*  the estimate for RCOND, and is almost always a slight
*  overestimate of the true error.
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
*  BERR is DOUBLE PRECISION array, dimension (NRHS)
*  The componentwise relative backward error of each solution
*  vector X(j) (i.e., the smallest relative change in
*   any element of A or B that makes X(j) an exact solution).
*/
    double getSolutionBackwardErrorEst()
    {
        return BERR[0];
    }

    std::vector<double> getMultipleSolutionBackwardErrorEst()
    {
        return BERR;
    }



    double              RCOND;
    std::vector<double>  FERR;
    std::vector<double>  BERR;

    LapackMatrix       A;
    LapackMatrix      AF;

};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DPOSV : Created for solving symmetric positive definite systems
// of equations
// LAPACK base routine description:
// DPOSV computes the solution to a real system of linear equations
//    A * X = B,
// where A is an N-by-N symmetric positive definite matrix and X and B
// are N-by-NRHS matrices.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class DPOSV
{
public:

    DPOSV()
    {initialize();}

    void initialize(){}

    void applyInverse(const LapackMatrix& A,std::vector <double >& b)
    {
            applyInverse(A,&b[0]);
    }

    void applyInverse(const LapackMatrix& A,LapackMatrix& b)
    {
            applyInverse(A,b.dataPtr,b.cols);
    }

    void applyInverse(const LapackMatrix& A, double* b, long NRHS = 1)
    {
        assert(A.sizeCheck(A.rows,A.cols));

        char UPLO  = 'U'; // A = U**T* U
        long N     = A.rows;

        double* Aptr  =  A.dataPtr;

        long LDA    = N;
        double* Bptr = b;
        long LDB    = N;
        long   INFO = 0;

        dposv_(&UPLO, &N, &NRHS,Aptr,&LDA,Bptr,&LDB,&INFO);

        if(INFO != 0)
        {
        std::cerr << "dposv Failed : INFO = " << INFO  << std::endl;
        exit(1);
        }
    }
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class NORMALEQ creates the solution of the normal equations with the singular value
// parameter cut-off value svdCutoff. Since the singular values of the normal
// equations are square roots of the eigenvalues of the normal equations, the
// accuracy of the solution is limited to problems where it is sufficient to
// construct a solution using singular vector components associated with
// singular values where (singular values)^2 > accuracy of the
// eigensystem routine. Typically this limits the application to problems where
// svdCutoff is approximately > 1.0e-07. The returned solution can only be
// expected to have single precision accuracy for the same reasons.
//
// In the implementation below, the solution components are only accumulated
// when the singular values are greater than the max(svdCutoff).
//
// When M < N, so the solution is underdetermined, the minimum L2
// norm solution is returned and only M approximate singular
// values are computed.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class NORMALEQ
{
public:

    NORMALEQ()
    {
    initialize();
    };

    void initialize()
    {
    bStar.clear();
    bBar.clear();
    ABnormal.initialize();

    diagonalizer.initialize();
    eigenValues.clear();
    eigenVectors.initialize();

    singularValues.clear();
    svdDim = 0;
    }

    long getSVDdim(){return svdDim;}


    // Normal equation solution of
    //
    // A*x = b
    //
    // with svdCutoff.
    //
    // Only components with singular values > svdCutoff or 1.0e-07 are accumulated.
    //
    // When A is M x N and M < N, then only M singular values are evaluated.
    //
    //
    std::vector<double> computeNormalEquationSolution(std::vector<double>& b, LapackMatrix& AB, double svdCutoff)
    {
        // Form right hand side

        bBar = AB.applyTranspose(b);

        // Form normal equations

        long M = AB.getRowDimension();
        long N = AB.getColDimension();

        std::vector<double> x(N,0.0);

        ABnormal.initialize(N,N);

        // Call dgemm to form normal equations. This is done using a direct call to dgemm so
        // that forming A^T isn't required.

        char TRANSA  = 'T';
        char TRANSB  = 'N';
        long  Mstar  = N;
        long  Nstar  = N;
        long Kstar   = M;
        double ALPHA = 1.0;
        double BETA  = 0.0;
        long LDA     = M;
        long LDB     = M;
        long LDC     = N;
        long INCX    = 1;
        long INCY    = 1;

        dgemm_(&TRANSA,&TRANSB,&Mstar,&Nstar,&Kstar,&ALPHA, AB.dataPtr,&LDA,AB.dataPtr, &LDB,&BETA,ABnormal.dataPtr,&LDC);

        // Diagonalize the normal equations

        eigenValues.resize(N,0.0);

        if(M >= N)
        {
        singularValues.resize(N,0.0);
        }
        else
        {
        singularValues.resize(M,0.0);
        }

        eigenVectors.initialize(N,N);
        diagonalizer.computeEigensystem(ABnormal, eigenValues, eigenVectors);


        // The eigenvalues from the diagonalization routine
        // are returned smallest to largest, so pack the singular values in reverse order
        // to agree with the SVD convention

        long      stopIndex = 0;
        if(M < N) stopIndex = N-M;

        long svIndex             = 0;
        long firstComponentIndex = N-1;
        long componentCount      = 0;
        for(long i= N-1; i >= stopIndex; i--)
        {
            singularValues[svIndex] = std::sqrt(std::abs(eigenValues[i]));
            if(singularValues[svIndex] > svdCutoff)
            {
              firstComponentIndex = i;
              componentCount++;
            }
            svIndex++;
        }

        /*
        std::cout << std::setprecision(15) << std::endl;
        for(long i = 0; i < eigenValues.size(); i++)
        {
            std::cout << std::sqrt(std::abs(eigenValues[i])) <<  " " << std::abs(eigenValues[i]) << std::endl;
        }
        */


        // Project solution onto selected subspace

        // A direct call to dgemv to create
        //
        // bStar = eigenVectors.applyTranspose(bBar);
        //
        // in order to avoid accumulating the values known
        // a-priori to be set to zero
        //

        bStar.resize(N,0.0);

        char TRANS = 'T';
        Mstar      =  N;
        Nstar      =  componentCount;
        ALPHA      = 1.0;
        BETA       = 0.0;
        INCX       = 1;
        INCY       = 1;
        LDA        = N;

        double* bPtr   = &bStar[firstComponentIndex];
        double* eigPtr = eigenVectors.dataPtr + firstComponentIndex*N;

        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,eigPtr,&LDA,&bBar[0],&INCX,&BETA,bPtr,&INCY);

        // Solve he diagonal system for the components being kept

        for(long i= N-1; i >= firstComponentIndex; i--)
        {
            bStar[i] /= std::abs(eigenValues[i]);
        }

        // Zero out entries in case M < N

        for(long i = firstComponentIndex-1; i >= 0; i--)
        {
            bStar[i] = 0;
        }

        //
        // Evaluate the solution. Using direct call to dgemv to avoid
        // accumulating columns with zero weights.
        //
        // x = eigenVectors*bStar;

        x.resize(M,0.0);

        TRANS      = 'N';
        Mstar      =  N;
        Nstar      =  componentCount;
        ALPHA      = 1.0;
        BETA       = 0.0;
        INCX       = 1;
        INCY       = 1;
        LDA        = N;

        bPtr   = &bStar[firstComponentIndex];
        eigPtr = eigenVectors.dataPtr + firstComponentIndex*N;

        dgemv_(&TRANS,&Mstar,&Nstar,&ALPHA,eigPtr,&LDA,bPtr,&INCX,&BETA,&x[0],&INCY);

        svdDim = componentCount;

        return x;

    }

    std::vector<double>      bStar;
    std::vector<double>       bBar;
    LapackMatrix          ABnormal;

    DSYEV                 diagonalizer;
    std::vector<double>    eigenValues;
    LapackMatrix          eigenVectors;

    std::vector<double> singularValues;
        long                    svdDim;

};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class QRutility can be used to create a solution to
//
//               A x = b
//
// where A is a full rank M x N matrix with M >= N.
// If M > N then a least squares approximation to the
// solution is determined.
//
// The solution process is decomposed into two
// steps; a create QR factors step followed by
// a create solution step. The QR factorization
// is retained internally so that, for a given
// system and multiple right hand sides, the QR
// factorization need only be computed once.
//
// The procedure used is essentially that of
// the LAPACK routine DGELS performed without scaling.
// That routine, as well as this class, creates
// a QR factorization without pivoting, so no
// rank determination is performed. For the solution
// of a system with undetermined properties, one
// might consider using SCC::DGELSY or another class
// that uses Lapack routine DGELSY to construct
// QR solutions.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
class QRutility
{
    public :

    QRutility()
    {
    initialize();
    }

    QRutility(const QRutility& Q)
    {
    initialize(Q);
    }

    void initialize()
    {
    QRfactors.initialize();
    TAU.clear();
    WORK.clear();
    bTemp.clear();
    dormqrWorkSize = -1;
    }

    void initialize(const QRutility& Q)
    {
    QRfactors.initialize(Q.QRfactors);
    TAU            = Q.TAU;
    WORK           = Q.WORK;
    bTemp          = Q.bTemp;
    dormqrWorkSize = Q.dormqrWorkSize;
    }

    std::vector<double> createQRsolution(std::vector<double>& b)
    {
    try {if(QRfactors.getRowDimension() != (long)b.size()) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "createSolution called before createQRfactors " << std::endl;
     exit(1);
    }

    // Capture right hand side

    bTemp = b;

    char SIDE    = 'L';
    char TRANS   = 'T';
    long M       = (long)b.size();
    long NRHS    = 1;
    long K       = QRfactors.getColDimension();
    long LDA     = QRfactors.getRowDimension();
    double* Aptr = QRfactors.getDataPointer();
    long LDC     = M;

    long LWORK  = -1;
    long INFO   =  0;

    double WORKDIM;

    // Obtain the optimal work size if it has not already
    // been determined

    if(dormqrWorkSize < 0)
    {
    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0], &bTemp[0],
    &LDC, &WORKDIM, &LWORK, &INFO);

    dormqrWorkSize = (long)WORKDIM+1;
    }

    LWORK = dormqrWorkSize;
    WORK.resize(LWORK);

    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0], &bTemp[0],
    &LDC, &WORK[0], &LWORK, &INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "DORMQR Failed : INFO = " << INFO  << std::endl;
     exit(1);
    }

    // Note: only using QRfactors.cols elements of bTemp.

    // Backsolve upper trignular system to obtain a solution

    char UPLO    = 'U';
    TRANS        = 'N';
    char DIAG    = 'N';
    M            = QRfactors.getColDimension();
    long N       = M;
    NRHS         = 1;
    LDA          = QRfactors.getRowDimension();
    long LDB     = M;

    dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, Aptr, &LDA,&bTemp[0],&LDB,&INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "DTRTRS detected singular matrix : INFO = " << INFO  << std::endl;
     exit(1);
    }

    return bTemp;
    }


    LapackMatrix createQRsolution(const LapackMatrix& B)
    {
    try {if(QRfactors.getRowDimension() != (long)B.rows) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "createSolution called before createQRfactors " << std::endl;
     exit(1);
    }

    LapackMatrix Btmp;

    // Capture right hand side

    Btmp.initialize(B);

    char SIDE    = 'L';
    char TRANS   = 'T';
    long M       = Btmp.rows;
    long NRHS    = Btmp.cols;
    long K       = QRfactors.getColDimension();
    long LDA     = QRfactors.getRowDimension();
    double* Aptr = QRfactors.getDataPointer();
    long LDC     = M;

    long LWORK  = -1;
    long INFO   =  0;

    double WORKDIM;

    // Obtain the optimal work size if it has not already
    // been determined

    if(dormqrWorkSize < 0)
    {
    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0], Btmp.getDataPointer(),
    &LDC, &WORKDIM, &LWORK, &INFO);

    dormqrWorkSize = (long)WORKDIM+1;
    }

    LWORK = dormqrWorkSize;
    WORK.resize(LWORK);

    dormqr_(&SIDE, &TRANS, &M, &NRHS, &K,Aptr,& LDA, &TAU[0],Btmp.getDataPointer(),
    &LDC, &WORK[0], &LWORK, &INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "DORMQR Failed : INFO = " << INFO  << std::endl;
     exit(1);
    }

    // Note: only using upper QRfactors.cols x QRfactors.cols elements of Btmp.

    // Backsolve upper trignular system to obtain a solution

    char UPLO    = 'U';
    TRANS        = 'N';
    char DIAG    = 'N';
    M            = QRfactors.getColDimension();
    long N       = M;
    NRHS         = Btmp.cols;
    LDA          = QRfactors.getRowDimension();
    long LDB     = Btmp.rows;


    dtrtrs_(&UPLO, &TRANS, &DIAG, &N, &NRHS, Aptr, &LDA,Btmp.getDataPointer(),&LDB,&INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRsolution");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "DTRTRS detected singular matrix : INFO = " << INFO  << std::endl;
     exit(1);
    }

    if(Btmp.rows == Btmp.cols)
    {
    return Btmp;
    }

    LapackMatrix Bstar(N,N);

    for(long i = 0; i < N; i++)
    {
    for(long j = 0; j < N; j++)
    {
    Bstar(i,j) = Btmp(i,j);
    }}

    return Bstar;
    }



    //
    // A convenience solve interface. It is assumed that Bptr points to contiguous
    // data of size of the number of rows of A. No bounds checking is performed.
    //

    std::vector<double> createQRsolution(double* Bptr)
    {
        long M  = QRfactors.getRowDimension();
        std::vector<double>                   B(M);
        std::memcpy(B.data(),Bptr,M*sizeof(double));
        return createQRsolution(B);
    }
//
//  This member function creates and stores internally the QR factorization
//  of the M x N input matrix A. It is assumed that M >= N and A is of
//  full rank.
//
    void createQRfactors(const SCC::LapackMatrix& A)
    {
    long M   = A.getRowDimension();
    long N   = A.getColDimension();

    try
    {
    if( M < N ) throw std::runtime_error("createQRfactors");
    }
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "Input matrix rows < cols " << '\n';
     std::cerr << "rows : " << M << " cols : " << N << std::endl;
     exit(1);
    }

    // Capture system

    QRfactors.initialize(A);

    // Create QR factors

    long LDA = M;

    TAU.resize(N);

    long LWORK  = -1;
    long INFO   =  0;

    double WORKDIM;

    // First call to obtain the optimal work size

    dgeqrf_(&M, &N, A.getDataPointer(), &LDA, &TAU[0],&WORKDIM, &LWORK, &INFO);

    LWORK = (long)WORKDIM+1;
    WORK.resize(LWORK);

    // Create QR factors

    dgeqrf_(&M, &N, QRfactors.getDataPointer(), &LDA, &TAU[0],&WORK[0], &LWORK, &INFO);

    try {if(INFO != 0) {throw std::runtime_error("createQRfactors");}}
    catch (std::runtime_error& e)
    {
     std::cerr << "Runtime exception in QRutility member function " <<  e.what() << std::endl;
     std::cerr << "DGEQRF Failed : INFO = " << INFO  << std::endl;
     exit(1);
    }

    // Reset work sizes to -1  to trigger a re-computation
    // of work storage requirements on first use
    // of createQRsolution

    dormqrWorkSize = -1;
    }

    SCC::LapackMatrix QRfactors; // QR factor component as returned by dgeqrf
    std::vector<double>          TAU; // QR factor component as returned by dgeqrf
    std::vector<double>         WORK;
    std::vector<double>        bTemp;

    long         dormqrWorkSize;
};

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class DSYEVX : Created for computing eigensystem of symmetric matricies
// LAPACK base routine description:
// DSYEVX computes selected eigenvalues and, optionally, eigenvectors
// of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
// selected by specifying either a range of values or a range of indices
// for the desired eigenvalues.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class DSYEVX
{
public :

    DSYEVX()
    {
    initialize();
    }

    void initialize()
    {
    A.initialize();
    WORK.clear();
    IWORK.clear();
    IFAIL.clear();
    }

    // Computes the eigCount algebraically smallest eigenvalues and eigenvectors.
    // The value returned is the number of eigenvalues found.

    long createAlgSmallestEigensystem(long eigCount, SCC::LapackMatrix& M, std::vector<double>& eigValues,
                                      SCC::LapackMatrix& eigVectors)
    {
    if(M.getRowDimension() != M.getColDimension())
    {
            throw std::runtime_error("\nDSYEVX : Non-square matrix input argument  \n");
    }

    long N = M.getRowDimension();

    if(eigCount > N)
    {
        std::stringstream sout;
        sout << "\nDSYEVX Error \n";
        sout << "Requested number of eigenvalues/eigenvectors exceeds system dimension. \n";
        throw std::runtime_error(sout.str());
    }
    // Copy input matrix

    A.initialize(M);

    char JOBZ   = 'V'; // Specify N for eigenvalues only
    char RANGE  = 'I'; // Specify index range of eigenvalues to be found (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    long LDA  = N;
    double VL = 0;
    double VU = 0;

    long IL = 1;        // Index of smallest eigenvalue returned
    long IU = eigCount; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long eigComputed = 0;     //  M parameter in original call = number of eigenvalues output

    eigValues.clear();         // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ   = N;
    long Mstar = (IU-IL) + 1;              // Maximal number of eigenvalues to be computed when using index specification

    eigVectors.initialize(LDZ,Mstar);      // Matrix whose columns containing the eigenvectors (Z in original call)

    long INFO = 0;

    WORK.clear();
    IWORK.clear();
    IFAIL.clear();

    // workspace query

    long LWORK = -1;
    WORK.resize(1,0.0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,&eigComputed,eigValues.data(),
    eigVectors.getDataPointer(),&LDZ,WORK.data(),&LWORK,IWORK.data(),IFAIL.data(),&INFO);

    LWORK  = (long)WORK[0];

    WORK.resize(LWORK,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,&eigComputed,eigValues.data(),
    eigVectors.getDataPointer(),&LDZ,WORK.data(),&LWORK,IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
        std::stringstream sout;
        sout << "\nDSYEVX \nError INFO = " << INFO << "\n";
        throw std::runtime_error(sout.str());
    }

    // resize the eig values array to the number of eigenvalues found

    eigValues.resize(eigComputed);
    return eigComputed;
    }


    // Computes the eigCount algebraically smallest eigenvalues and eigenvectors.
    // The value returned is the number of eigenvalues found.

    long createEigensystem(SCC::LapackMatrix& M, std::vector<double>& eigValues,
                           SCC::LapackMatrix& eigVectors)
    {
     if(M.getRowDimension() != M.getColDimension())
    {
            throw std::runtime_error("\nDSYEVX : Non-square matrix input argument  \n");
    }

    long N = M.getRowDimension();

    // Copy input matrix

    A.initialize(M);

    char JOBZ   = 'V'; // Specify N for eigenvalues only
    char RANGE  = 'A'; // Specify index range of eigenvalues to be found (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    long LDA  = N;
    double VL = 0;
    double VU = 0;

    long IL = 0; // Index of smallest eigenvalue returned
    long IU = 0; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long eigComputed = 0;     //  M parameter in original call = number of eigenvalues output

    eigValues.clear();         // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ     = N;
    eigVectors.initialize(LDZ,N);      // Matrix whose columns containing the eigenvectors (Z in original call)

    long INFO = 0;

    WORK.clear();
    IWORK.clear();
    IFAIL.clear();

    // workspace query

    long LWORK = -1;
    WORK.resize(1,0.0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,&eigComputed,eigValues.data(),
    eigVectors.getDataPointer(),&LDZ,WORK.data(),&LWORK,IWORK.data(),IFAIL.data(),&INFO);

    LWORK  = (long)WORK[0];

    WORK.resize(LWORK,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,&eigComputed,eigValues.data(),
    eigVectors.getDataPointer(),&LDZ,WORK.data(),&LWORK,IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
        std::stringstream sout;
        sout << "\nDSYEVX \nError INFO = " << INFO << "\n";
        throw std::runtime_error(sout.str());
    }

    // resize the eig values array to the number of eigenvalues found

    eigValues.resize(eigComputed);
    return eigComputed;
    }

    long createAlgSmallestEigenvalues(long eigCount, SCC::LapackMatrix& M, std::vector<double>& eigValues)
    {
        if(M.getRowDimension() != M.getColDimension())
    {
            throw std::runtime_error("\nDSYEVX : Non-square matrix input argument  \n");
    }

    long N = M.getRowDimension();

    if(eigCount > N)
    {
        std::stringstream sout;
        sout << "\nDSYEVX Error \n";
        sout << "Requested number of eigenvalues/eigenvectors exceeds system dimension. \n";
        throw std::runtime_error(sout.str());
    }
    // Copy input matrix

    A.initialize(M);

    char JOBZ   = 'N'; // Specify N for eigenvalues only
    char RANGE  = 'I'; // Specify index range of eigenvalues to be found (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    long LDA  = N;
    double VL = 0;
    double VU = 0;

    long IL = 1;        // Index of smallest eigenvalue returned
    long IU = eigCount; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long eigComputed = 0;     //  M parameter in original call = number of eigenvalues output

    eigValues.clear();         // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ     = 1;
    double ZDATA = 0.0;

    long INFO = 0;

    WORK.clear();
    IWORK.clear();
    IFAIL.clear();

    // workspace query

    long LWORK = -1;
    WORK.resize(1,0.0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,
    		&eigComputed,eigValues.data(), &ZDATA,&LDZ,WORK.data(),&LWORK,
			IWORK.data(),IFAIL.data(),&INFO);

    LWORK  = (long)WORK[0];

    WORK.resize(LWORK,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);

    dsyevx_(&JOBZ, &RANGE, &UPLO,&N,A.getDataPointer(),&LDA, &VL,&VU,&IL,&IU,&ABSTOL,&eigComputed,eigValues.data(),
    &ZDATA,&LDZ,WORK.data(),&LWORK,IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
        std::stringstream sout;
        sout << "\nDSYEVX \nError INFO = " << INFO << "\n";
        throw std::runtime_error(sout.str());
    }

    // resize the eig values array to the number of eigenvalues found

    eigValues.resize(eigComputed);
    return eigComputed;
    }


    SCC::LapackMatrix             A;
    std::vector<double>         WORK;
    std::vector<long>          IWORK;
    std::vector<long>          IFAIL;


};


} // End of SCC Namespace declaration

/////////////////////////////////////////////////////////////////////////////
// LAPACK documentation
/////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
// DGELSY
/////////////////////////////////////////////////////////////////////////////
/*
DGELSY solves overdetermined or underdetermined systems for GE matrices

subroutine dgelsy    (    integer     m,
integer     n,
integer     nrhs,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( ldb, * )     b,
integer     ldb,
integer, dimension( * )     jpvt,
double precision     rcond,
integer     rank,
double precision, dimension( * )     work,
integer     lwork,
integer     info
)

Purpose:
 DGELSY computes the minimum-norm solution to a real linear least
 squares problem:
     minimize || A * X - B ||
 using a complete orthogonal factorization of A.  A is an M-by-N
 matrix which may be rank-deficient.

 Several right hand side vectors b and solution vectors x can be
 handled in a single call; they are stored as the columns of the
 M-by-NRHS right hand side matrix B and the N-by-NRHS solution
 matrix X.

 The routine first computes a QR factorization with column pivoting:
     A * P = Q * [ R11 R12 ]
                 [  0  R22 ]
 with R11 defined as the largest leading submatrix whose estimated
 condition number is less than 1/RCOND.  The order of R11, RANK,
 is the effective rank of A.

 Then, R22 is considered to be negligible, and R12 is annihilated
 by orthogonal transformations from the right, arriving at the
 complete orthogonal factorization:
    A * P = Q * [ T11 0 ] * Z
                [  0  0 ]
 The minimum-norm solution is then
    X = P * Z**T [ inv(T11)*Q1**T*B ]
                 [        0         ]
 where Q1 consists of the first RANK columns of Q.

 This routine is basically identical to the original xGELSX except
 three differences:
   o The call to the subroutine xGEQPF has been substituted by the
     the call to the subroutine xGEQP3. This subroutine is a Blas-3
     version of the QR factorization with column pivoting.
   o Matrix B (the right hand side) is updated with Blas-3.
   o The permutation of matrix B (the right hand side) is faster and
     more simple.
Parameters
[in]    M
          M is INTEGER
          The number of rows of the matrix A.  M >= 0.
[in]    N
          N is INTEGER
          The number of columns of the matrix A.  N >= 0.
[in]    NRHS
          NRHS is INTEGER
          The number of right hand sides, i.e., the number of
          columns of matrices B and X. NRHS >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit, A has been overwritten by details of its
          complete orthogonal factorization.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
[in,out]    B
          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the M-by-NRHS right hand side matrix B.
          On exit, the N-by-NRHS solution matrix X.
          If M = 0 or N = 0, B is not referenced.
[in]    LDB
          LDB is INTEGER
          The leading dimension of the array B. LDB >= max(1,M,N).
[in,out]    JPVT
          JPVT is INTEGER array, dimension (N)
          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
          to the front of AP, otherwise column i is a free column.
          On exit, if JPVT(i) = k, then the i-th column of AP
          was the k-th column of A.
[in]    RCOND
          RCOND is DOUBLE PRECISION
          RCOND is used to determine the effective rank of A, which
          is defined as the order of the largest leading triangular
          submatrix R11 in the QR factorization with pivoting of A,
          whose estimated condition number < 1/RCOND.
[out]    RANK
          RANK is INTEGER
          The effective rank of A, i.e., the order of the submatrix
          R11.  This is the same as the order of the submatrix T11
          in the complete orthogonal factorization of A.
          If NRHS = 0, RANK = 0 on output.
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]    LWORK
          LWORK is INTEGER
          The dimension of the array WORK.
          The unblocked strategy requires that:
             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
          where MN = min( M, N ).
          The block algorithm requires that:
             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
          where NB is an upper bound on the blocksize returned
          by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR,
          and DORMRZ.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    INFO
          INFO is INTEGER
          = 0: successful exit
          < 0: If INFO = -i, the i-th argument had an illegal value.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

/////////////////////////////////////////////////////////////////////////////
// DGESVD
/////////////////////////////////////////////////////////////////////////////
/*
DGESVD computes the singular value decomposition (SVD) for GE matrices

subroutine dgesvd    (    character     jobu,
character     jobvt,
integer     m,
integer     n,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( * )     s,
double precision, dimension( ldu, * )     u,
integer     ldu,
double precision, dimension( ldvt, * )     vt,
integer     ldvt,
double precision, dimension( * )     work,
integer     lwork,
integer     info
)

Purpose:
 DGESVD computes the singular value decomposition (SVD) of a real
 M-by-N matrix A, optionally computing the left and/or right singular
 vectors. The SVD is written

      A = U * SIGMA * transpose(V)

 where SIGMA is an M-by-N matrix which is zero except for its
 min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 are the singular values of A; they are real and non-negative, and
 are returned in descending order.  The first min(m,n) columns of
 U and V are the left and right singular vectors of A.

 Note that the routine returns V**T, not V.
Parameters
[in]    JOBU
          JOBU is CHARACTER*1
          Specifies options for computing all or part of the matrix U:
          = 'A':  all M columns of U are returned in array U:
          = 'S':  the first min(m,n) columns of U (the left singular
                  vectors) are returned in the array U;
          = 'O':  the first min(m,n) columns of U (the left singular
                  vectors) are overwritten on the array A;
          = 'N':  no columns of U (no left singular vectors) are
                  computed.
[in]    JOBVT
          JOBVT is CHARACTER*1
          Specifies options for computing all or part of the matrix
          V**T:
          = 'A':  all N rows of V**T are returned in the array VT;
          = 'S':  the first min(m,n) rows of V**T (the right singular
                  vectors) are returned in the array VT;
          = 'O':  the first min(m,n) rows of V**T (the right singular
                  vectors) are overwritten on the array A;
          = 'N':  no rows of V**T (no right singular vectors) are
                  computed.

          JOBVT and JOBU cannot both be 'O'.
[in]    M
          M is INTEGER
          The number of rows of the input matrix A.  M >= 0.
[in]    N
          N is INTEGER
          The number of columns of the input matrix A.  N >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit,
          if JOBU = 'O',  A is overwritten with the first min(m,n)
                          columns of U (the left singular vectors,
                          stored columnwise);
          if JOBVT = 'O', A is overwritten with the first min(m,n)
                          rows of V**T (the right singular vectors,
                          stored rowwise);
          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
                          are destroyed.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
[out]    S
          S is DOUBLE PRECISION array, dimension (min(M,N))
          The singular values of A, sorted so that S(i) >= S(i+1).
[out]    U
          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
          if JOBU = 'S', U contains the first min(m,n) columns of U
          (the left singular vectors, stored columnwise);
          if JOBU = 'N' or 'O', U is not referenced.
[in]    LDU
          LDU is INTEGER
          The leading dimension of the array U.  LDU >= 1; if
          JOBU = 'S' or 'A', LDU >= M.
[out]    VT
          VT is DOUBLE PRECISION array, dimension (LDVT,N)
          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
          V**T;
          if JOBVT = 'S', VT contains the first min(m,n) rows of
          V**T (the right singular vectors, stored rowwise);
          if JOBVT = 'N' or 'O', VT is not referenced.
[in]    LDVT
          LDVT is INTEGER
          The leading dimension of the array VT.  LDVT >= 1; if
          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
          superdiagonal elements of an upper bidiagonal matrix B
          whose diagonal is in S (not necessarily sorted). B
          satisfies A = U * B * VT, so it has the same singular values
          as A, and singular vectors related by U and VT.
[in]    LWORK
          LWORK is INTEGER
          The dimension of the array WORK.
          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
             - PATH 1  (M much larger than N, JOBU='N')
             - PATH 1t (N much larger than M, JOBVT='N')
          LWORK >= MAX(1,3*MIN(M,N) + MAX(M,N),5*MIN(M,N)) for the other paths
          For good performance, LWORK should generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit.
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if DBDSQR did not converge, INFO specifies how many
                superdiagonals of an intermediate bidiagonal form B
                did not converge to zero. See the description of WORK
                above for details.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.

*/

/////////////////////////////////////////////////////////////////////////////
// DSYEV
/////////////////////////////////////////////////////////////////////////////
/*
DSYEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for SY matrices
dsyev()
subroutine dsyev    (    character     jobz,
character     uplo,
integer     n,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( * )     w,
double precision, dimension( * )     work,
integer     lwork,
integer     info
)

Purpose:
 DSYEV computes all eigenvalues and, optionally, eigenvectors of a
 real symmetric matrix A.
Parameters
[in]    JOBZ
          JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
[in]    UPLO
          UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.
[in]    N
          N is INTEGER
          The order of the matrix A.  N >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA, N)
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A.
          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          orthonormal eigenvectors of the matrix A.
          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
          or the upper triangle (if UPLO='U') of A, including the
          diagonal, is destroyed.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[out]    W
          W is DOUBLE PRECISION array, dimension (N)
          If INFO = 0, the eigenvalues in ascending order.
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]    LWORK
          LWORK is INTEGER
          The length of the array WORK.  LWORK >= max(1,3*N-1).
          For optimal efficiency, LWORK >= (NB+2)*N,
          where NB is the blocksize for DSYTRD returned by ILAENV.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, the algorithm failed to converge; i
                off-diagonal elements of an intermediate tridiagonal
                form did not converge to zero.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
 */

/////////////////////////////////////////////////////////////////////////////
// DGESVX
/////////////////////////////////////////////////////////////////////////////
/*
DGESVX solves linear systems of equations for GE matrices

subroutine dgesvx     (     character      fact,
        character      trans,
        integer      n,
        integer      nrhs,
        double precision, dimension( lda, * )      a,
        integer      lda,
        double precision, dimension( ldaf, * )      af,
        integer      ldaf,
        integer, dimension( * )      ipiv,
        character      equed,
        double precision, dimension( * )      r,
        double precision, dimension( * )      c,
        double precision, dimension( ldb, * )      b,
        integer      ldb,
        double precision, dimension( ldx, * )      x,
        integer      ldx,
        double precision      rcond,
        double precision, dimension( * )      ferr,
        double precision, dimension( * )      berr,
        double precision, dimension( * )      work,
        integer, dimension( * )      iwork,
        integer      info
    )

Purpose:

     DGESVX uses the LU factorization to compute the solution to a real
     system of linear equations
        A * X = B,
     where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

     Error bounds on the solution and a condition estimate are also
     provided.

Description:

     The following steps are performed:

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
           A = P * L * U,
        where P is a permutation matrix, L is a unit lower triangular
        matrix, and U is upper triangular.

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
    [in]    FACT

              FACT is CHARACTER*1
              Specifies whether or not the factored form of the matrix A is
              supplied on entry, and if not, whether the matrix A should be
              equilibrated before it is factored.
              = 'F':  On entry, AF and IPIV contain the factored form of A.
                      If EQUED is not 'N', the matrix A has been
                      equilibrated with scaling factors given by R and C.
                      A, AF, and IPIV are not modified.
              = 'N':  The matrix A will be copied to AF and factored.
              = 'E':  The matrix A will be equilibrated if necessary, then
                      copied to AF and factored.

    [in]    TRANS

              TRANS is CHARACTER*1
              Specifies the form of the system of equations:
              = 'N':  A * X = B     (No transpose)
              = 'T':  A**T * X = B  (Transpose)
              = 'C':  A**H * X = B  (Transpose)

    [in]    N

              N is INTEGER
              The number of linear equations, i.e., the order of the
              matrix A.  N >= 0.

    [in]    NRHS

              NRHS is INTEGER
              The number of right hand sides, i.e., the number of columns
              of the matrices B and X.  NRHS >= 0.

    [in,out]    A

              A is DOUBLE PRECISION array, dimension (LDA,N)
              On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
              not 'N', then A must have been equilibrated by the scaling
              factors in R and/or C.  A is not modified if FACT = 'F' or
              'N', or if FACT = 'E' and EQUED = 'N' on exit.

              On exit, if EQUED .ne. 'N', A is scaled as follows:
              EQUED = 'R':  A := diag(R) * A
              EQUED = 'C':  A := A * diag(C)
              EQUED = 'B':  A := diag(R) * A * diag(C).

    [in]    LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

    [in,out]    AF

              AF is DOUBLE PRECISION array, dimension (LDAF,N)
              If FACT = 'F', then AF is an input argument and on entry
              contains the factors L and U from the factorization
              A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then
              AF is the factored form of the equilibrated matrix A.

              If FACT = 'N', then AF is an output argument and on exit
              returns the factors L and U from the factorization A = P*L*U
              of the original matrix A.

              If FACT = 'E', then AF is an output argument and on exit
              returns the factors L and U from the factorization A = P*L*U
              of the equilibrated matrix A (see the description of A for
              the form of the equilibrated matrix).

    [in]    LDAF

              LDAF is INTEGER
              The leading dimension of the array AF.  LDAF >= max(1,N).

    [in,out]    IPIV

              IPIV is INTEGER array, dimension (N)
              If FACT = 'F', then IPIV is an input argument and on entry
              contains the pivot indices from the factorization A = P*L*U
              as computed by DGETRF; row i of the matrix was interchanged
              with row IPIV(i).

              If FACT = 'N', then IPIV is an output argument and on exit
              contains the pivot indices from the factorization A = P*L*U
              of the original matrix A.

              If FACT = 'E', then IPIV is an output argument and on exit
              contains the pivot indices from the factorization A = P*L*U
              of the equilibrated matrix A.

    [in,out]    EQUED

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

    [in,out]    R

              R is DOUBLE PRECISION array, dimension (N)
              The row scale factors for A.  If EQUED = 'R' or 'B', A is
              multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
              is not accessed.  R is an input argument if FACT = 'F';
              otherwise, R is an output argument.  If FACT = 'F' and
              EQUED = 'R' or 'B', each element of R must be positive.

    [in,out]    C

              C is DOUBLE PRECISION array, dimension (N)
              The column scale factors for A.  If EQUED = 'C' or 'B', A is
              multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
              is not accessed.  C is an input argument if FACT = 'F';
              otherwise, C is an output argument.  If FACT = 'F' and
              EQUED = 'C' or 'B', each element of C must be positive.

    [in,out]    B

              B is DOUBLE PRECISION array, dimension (LDB,NRHS)
              On entry, the N-by-NRHS right hand side matrix B.
              On exit,
              if EQUED = 'N', B is not modified;
              if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
              diag(R)*B;
              if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
              overwritten by diag(C)*B.

    [in]    LDB

              LDB is INTEGER
              The leading dimension of the array B.  LDB >= max(1,N).

    [out]    X

              X is DOUBLE PRECISION array, dimension (LDX,NRHS)
              If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
              to the original system of equations.  Note that A and B are
              modified on exit if EQUED .ne. 'N', and the solution to the
              equilibrated system is inv(diag(C))*X if TRANS = 'N' and
              EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
              and EQUED = 'R' or 'B'.

    [in]    LDX

              LDX is INTEGER
              The leading dimension of the array X.  LDX >= max(1,N).

    [out]    RCOND

              RCOND is DOUBLE PRECISION
              The estimate of the reciprocal condition number of the matrix
              A after equilibration (if done).  If RCOND is less than the
              machine precision (in particular, if RCOND = 0), the matrix
              is singular to working precision.  This condition is
              indicated by a return code of INFO > 0.

    [out]    FERR

              FERR is DOUBLE PRECISION array, dimension (NRHS)
              The estimated forward error bound for each solution vector
              X(j) (the j-th column of the solution matrix X).
              If XTRUE is the true solution corresponding to X(j), FERR(j)
              is an estimated upper bound for the magnitude of the largest
              element in (X(j) - XTRUE) divided by the magnitude of the
              largest element in X(j).  The estimate is as reliable as
              the estimate for RCOND, and is almost always a slight
              overestimate of the true error.

    [out]    BERR

              BERR is DOUBLE PRECISION array, dimension (NRHS)
              The componentwise relative backward error of each solution
              vector X(j) (i.e., the smallest relative change in
              any element of A or B that makes X(j) an exact solution).

    [out]    WORK

              WORK is DOUBLE PRECISION array, dimension (MAX(1,4*N))
              On exit, WORK(1) contains the reciprocal pivot growth
              factor norm(A)/norm(U). The "max absolute element" norm is
              used. If WORK(1) is much less than 1, then the stability
              of the LU factorization of the (equilibrated) matrix A
              could be poor. This also means that the solution X, condition
              estimator RCOND, and forward error bound FERR could be
              unreliable. If factorization fails with 0<INFO<=N, then
              WORK(1) contains the reciprocal pivot growth factor for the
              leading INFO columns of A.

    [out]    IWORK

              IWORK is INTEGER array, dimension (N)

    [out]    INFO

              INFO is INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal value
              > 0:  if INFO = i, and i is
                    <= N:  U(i,i) is exactly zero.  The factorization has
                           been completed, but the factor U is exactly
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

/////////////////////////////////////////////////////////////////////////////
// DPOSV
/////////////////////////////////////////////////////////////////////////////
/*
DPOSV computes the solution to system of linear equations A * X = B for positive definite matrices
subroutine dposv    (    character     uplo,
integer     n,
integer     nrhs,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( ldb, * )     b,
integer     ldb,
integer     info
)


Purpose:
 DPOSV computes the solution to a real system of linear equations
    A * X = B,
 where A is an N-by-N symmetric positive definite matrix and X and B
 are N-by-NRHS matrices.

 The Cholesky decomposition is used to factor A as
    A = U**T* U,  if UPLO = 'U', or
    A = L * L**T,  if UPLO = 'L',
 where U is an upper triangular matrix and L is a lower triangular
 matrix.  The factored form of A is then used to solve the system of
 equations A * X = B.
Parameters
[in]    UPLO
          UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.
[in]    N
          N is INTEGER
          The number of linear equations, i.e., the order of the
          matrix A.  N >= 0.
[in]    NRHS
          NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
          N-by-N upper triangular part of A contains the upper
          triangular part of the matrix A, and the strictly lower
          triangular part of A is not referenced.  If UPLO = 'L', the
          leading N-by-N lower triangular part of A contains the lower
          triangular part of the matrix A, and the strictly upper
          triangular part of A is not referenced.

          On exit, if INFO = 0, the factor U or L from the Cholesky
          factorization A = U**T*U or A = L*L**T.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[in,out]    B
          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the N-by-NRHS right hand side matrix B.
          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
[in]    LDB
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
          > 0:  if INFO = i, the leading principal minor of order i
                of A is not positive, so the factorization could not
                be completed, and the solution has not been computed.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/
/////////////////////////////////////////////////////////////////////////////
// LAPACK routines used by QRutility
/////////////////////////////////////////////////////////////////////////////
/*
////////////////////////////////////////////////////////////////
DGEQRF computes a QR factorization of a real M-by-N matrix A:
////////////////////////////////////////////////////////////////

subroutine dgeqrf    (    integer     m,
integer     n,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( * )     tau,
double precision, dimension( * )     work,
integer     lwork,
integer     info
)

Purpose:
 DGEQRF computes a QR factorization of a real M-by-N matrix A:

    A = Q * ( R ),
            ( 0 )
 where:
    Q is a M-by-M orthogonal matrix;
    R is an upper-triangular N-by-N matrix;
    0 is a (M-N)-by-N zero matrix, if M > N.

Parameters
[in]    M
          M is INTEGER
          The number of rows of the matrix A.  M >= 0.
[in]    N
          N is INTEGER
          The number of columns of the matrix A.  N >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA,N)
          On entry, the M-by-N matrix A.
          On exit, the elements on and above the diagonal of the array
          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          upper triangular if m >= n); the elements below the diagonal,
          with the array TAU, represent the orthogonal matrix Q as a
          product of min(m,n) elementary reflectors (see Further
          Details).
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,M).
[out]    TAU
          TAU is DOUBLE PRECISION array, dimension (min(M,N))
          The scalar factors of the elementary reflectors (see Further
          Details).
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]    LWORK
          LWORK is INTEGER
          The dimension of the array WORK.
          LWORK >= 1, if MIN(M,N) = 0, and LWORK >= N, otherwise.
          For optimum performance LWORK >= N*NB, where NB is
          the optimal blocksize.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.

////////////////////////////////////////////////////////////////
DORMQR applies the Q factor of a QR decomposition to a matrix
////////////////////////////////////////////////////////////////

subroutine dormqr    (    character     side,
character     trans,
integer     m,
integer     n,
integer     k,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( * )     tau,
double precision, dimension( ldc, * )     c,
integer     ldc,
double precision, dimension( * )     work,
integer     lwork,
integer     info
)

Purpose:
 DORMQR overwrites the general real M-by-N matrix C with

                 SIDE = 'L'     SIDE = 'R'
 TRANS = 'N':      Q * C          C * Q
 TRANS = 'T':      Q**T * C       C * Q**T

 where Q is a real orthogonal matrix defined as the product of k
 elementary reflectors

       Q = H(1) H(2) . . . H(k)

 as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
 if SIDE = 'R'.
Parameters
[in]    SIDE
          SIDE is CHARACTER*1
          = 'L': apply Q or Q**T from the Left;
          = 'R': apply Q or Q**T from the Right.
[in]    TRANS
          TRANS is CHARACTER*1
          = 'N':  No transpose, apply Q;
          = 'T':  Transpose, apply Q**T.
[in]    M
          M is INTEGER
          The number of rows of the matrix C. M >= 0.
[in]    N
          N is INTEGER
          The number of columns of the matrix C. N >= 0.
[in]    K
          K is INTEGER
          The number of elementary reflectors whose product defines
          the matrix Q.
          If SIDE = 'L', M >= K >= 0;
          if SIDE = 'R', N >= K >= 0.
[in]    A
          A is DOUBLE PRECISION array, dimension (LDA,K)
          The i-th column must contain the vector which defines the
          elementary reflector H(i), for i = 1,2,...,k, as returned by
          DGEQRF in the first k columns of its array argument A.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.
          If SIDE = 'L', LDA >= max(1,M);
          if SIDE = 'R', LDA >= max(1,N).
[in]    TAU
          TAU is DOUBLE PRECISION array, dimension (K)
          TAU(i) must contain the scalar factor of the elementary
          reflector H(i), as returned by DGEQRF.
[in,out]    C
          C is DOUBLE PRECISION array, dimension (LDC,N)
          On entry, the M-by-N matrix C.
          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
[in]    LDC
          LDC is INTEGER
          The leading dimension of the array C. LDC >= max(1,M).
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]    LWORK
          LWORK is INTEGER
          The dimension of the array WORK.
          If SIDE = 'L', LWORK >= max(1,N);
          if SIDE = 'R', LWORK >= max(1,M).
          For good performance, LWORK should generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.

////////////////////////////////////////////////////////////////
DTRTRS solves a triangular system
////////////////////////////////////////////////////////////////

dtrtrs()
subroutine dtrtrs    (    character     uplo,
character     trans,
character     diag,
integer     n,
integer     nrhs,
double precision, dimension( lda, * )     a,
integer     lda,
double precision, dimension( ldb, * )     b,
integer     ldb,
integer     info
)

Purpose:
 DTRTRS solves a triangular system of the form

    A * X = B  or  A**T * X = B,

 where A is a triangular matrix of order N, and B is an N-by-NRHS
 matrix.  A check is made to verify that A is nonsingular.
Parameters
[in]    UPLO
          UPLO is CHARACTER*1
          = 'U':  A is upper triangular;
          = 'L':  A is lower triangular.
[in]    TRANS
          TRANS is CHARACTER*1
          Specifies the form of the system of equations:
          = 'N':  A * X = B  (No transpose)
          = 'T':  A**T * X = B  (Transpose)
          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
[in]    DIAG
          DIAG is CHARACTER*1
          = 'N':  A is non-unit triangular;
          = 'U':  A is unit triangular.
[in]    N
          N is INTEGER
          The order of the matrix A.  N >= 0.
[in]    NRHS
          NRHS is INTEGER
          The number of right hand sides, i.e., the number of columns
          of the matrix B.  NRHS >= 0.
[in]    A
          A is DOUBLE PRECISION array, dimension (LDA,N)
          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
          upper triangular part of the array A contains the upper
          triangular matrix, and the strictly lower triangular part of
          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
          triangular part of the array A contains the lower triangular
          matrix, and the strictly upper triangular part of A is not
          referenced.  If DIAG = 'U', the diagonal elements of A are
          also not referenced and are assumed to be 1.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[in,out]    B
          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
          On entry, the right hand side matrix B.
          On exit, if INFO = 0, the solution matrix X.
[in]    LDB
          LDB is INTEGER
          The leading dimension of the array B.  LDB >= max(1,N).
[out]    INFO
          INFO is INTEGER
          = 0:  successful exit
          < 0: if INFO = -i, the i-th argument had an illegal value
          > 0: if INFO = i, the i-th diagonal element of A is zero,
               indicating that the matrix is singular and the solutions
               X have not been computed.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

/////////////////////////////////////////////////////////////////////////////
// DSYEVX
/////////////////////////////////////////////////////////////////////////////

/*
DSYEVX computes the eigenvalues and, optionally, the left and/or right
eigenvectors for SY matrices

subroutine dsyevx    (    character     jobz,
character     range,
character     uplo,
integer     n,
double precision, dimension( lda, * )     a,
integer     lda,
double precision     vl,
double precision     vu,
integer     il,
integer     iu,
double precision     abstol,
integer     m,
double precision, dimension( * )     w,
double precision, dimension( ldz, * )     z,
integer     ldz,
double precision, dimension( * )     work,
integer     lwork,
integer, dimension( * )     iwork,
integer, dimension( * )     ifail,
integer     info
)

Purpose:
 DSYEVX computes selected eigenvalues and, optionally, eigenvectors
 of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
 selected by specifying either a range of values or a range of indices
 for the desired eigenvalues.

Parameters
[in]    JOBZ
          JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
[in]    RANGE
          RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found.
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found.
          = 'I': the IL-th through IU-th eigenvalues will be found.
[in]    UPLO
          UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.
[in]    N
          N is INTEGER
          The order of the matrix A.  N >= 0.
[in,out]    A
          A is DOUBLE PRECISION array, dimension (LDA, N)
          On entry, the symmetric matrix A.  If UPLO = 'U', the
          leading N-by-N upper triangular part of A contains the
          upper triangular part of the matrix A.  If UPLO = 'L',
          the leading N-by-N lower triangular part of A contains
          the lower triangular part of the matrix A.
          On exit, the lower triangle (if UPLO='L') or the upper
          triangle (if UPLO='U') of A, including the diagonal, is
          destroyed.
[in]    LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[in]    VL
          VL is DOUBLE PRECISION
          If RANGE='V', the lower bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]    VU
          VU is DOUBLE PRECISION
          If RANGE='V', the upper bound of the interval to
          be searched for eigenvalues. VL < VU.
          Not referenced if RANGE = 'A' or 'I'.
[in]    IL
          IL is INTEGER
          If RANGE='I', the index of the
          smallest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]    IU
          IU is INTEGER
          If RANGE='I', the index of the
          largest eigenvalue to be returned.
          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
          Not referenced if RANGE = 'A' or 'V'.
[in]    ABSTOL
          ABSTOL is DOUBLE PRECISION
          The absolute error tolerance for the eigenvalues.
          An approximate eigenvalue is accepted as converged
          when it is determined to lie in an interval [a,b]
          of width less than or equal to

                  ABSTOL + EPS *   max( |a|,|b| ) ,

          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing A to tridiagonal form.

          Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
          If this routine returns with INFO>0, indicating that some
          eigenvectors did not converge, try setting ABSTOL to
          2*DLAMCH('S').

          See "Computing Small Singular Values of Bidiagonal Matrices
          with Guaranteed High Relative Accuracy," by Demmel and
          Kahan, LAPACK Working Note #3.
[out]    M
          M is INTEGER
          The total number of eigenvalues found.  0 <= M <= N.
          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
[out]    W
          W is DOUBLE PRECISION array, dimension (N)
          On normal exit, the first M elements contain the selected
          eigenvalues in ascending order.
[out]    Z
          Z is DOUBLE PRECISION array, dimension (LDZ, max(1,M))
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and the
          index of the eigenvector is returned in IFAIL.
          If JOBZ = 'N', then Z is not referenced.
          Note: the user must ensure that at least max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used.
[in]    LDZ
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).
[out]    WORK
          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]    LWORK
          LWORK is INTEGER
          The length of the array WORK.  LWORK >= 1, when N <= 1;
          otherwise 8*N.
          For optimal efficiency, LWORK >= (NB+3)*N,
          where NB is the max of the blocksize for DSYTRD and DORMTR
          returned by ILAENV.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]    IWORK
          IWORK is INTEGER array, dimension (5*N)
[out]    IFAIL
          IFAIL is INTEGER array, dimension (N)
          If JOBZ = 'V', then if INFO = 0, the first M elements of
          IFAIL are zero.  If INFO > 0, then IFAIL contains the
          indices of the eigenvectors that failed to converge.
          If JOBZ = 'N', then IFAIL is not referenced.
[out]    INFO
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
#endif


