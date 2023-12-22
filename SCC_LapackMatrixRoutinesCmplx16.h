
//
// SCC::LapackBandMatrixRoutinesCmplx16
//
// A collection of utility classes whose functionality is
// based upon LAPACK routines. These routines are meant
// to be used with instances of LapackMatrixCmplx16. The documentation
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
// Data mapping being used for direct invocation of
// Fortran routines
//
// C++  int                 ==  Fortran LOGICAL
// C++  long                ==  Fortran INTEGER
// C++  double              ==  Fortran DOUBLE PRECISION
// C++ std::complex<double> ==  Fortran COMPLEX*16
//
//
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Current class list
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class  ZGESVX : Created for solving A*X = B using LU factorization
//                 when A is a complex matrix.
// LAPACK base routine description:
// ZGESVX uses the LU factorization to compute the solution to a complex
// system of linear equations
//    A * X = B,
// where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
// Error bounds on the solution and a condition estimate are also
// provided.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class  ZHPEVX : Created for computing eigensystem components
//                 of a Hermitian complex matrix.
//
// LAPACK base routine description:
// ZHPEVX computes selected eigenvalues and, optionally, eigenvectors
// of a complex Hermitian matrix A in packed storage.
// Eigenvalues/vectors can be selected by specifying either a range of
// values or a range of indices for the desired eigenvalues.
//
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class  ZGEEVX : Created for computing eigensystem components
//                 of a complex matrix.
//
// LAPACK base routine description:
// ZGEEVX computes for an N-by-N complex nonsymmetric matrix A, the
// eigenvalues and, optionally, the left and/or right eigenvectors.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class   ZGEESX: Created for computing Schur decomposition of a general
//                 complex matrix.
// LAPACK base routine description:
// ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the
// eigenvalues, the Schur form T, and, optionally, the matrix of Schur
// vectors Z.
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "SCC_LapackHeaders.h"
#include "SCC_LapackMatrix.h"
#include "SCC_LapackMatrixCmplx16.h"

#include <complex>
#include <cassert>

#ifndef SCC_LAPACK_MATRIX_ROUTINES_CMPLX16
#define SCC_LAPACK_MATRIX_ROUTINES_CMPLX16

namespace SCC
{


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class  ZGESVX : Created for solving A*X = B using LU factorization
//                 when A is a complex matrix.
// LAPACK base routine description:
// ZGESVX uses the LU factorization to compute the solution to a complex
// system of linear equations
//    A * X = B,
// where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
//
// Error bounds on the solution and a condition estimate are also
// provided.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class ZGESVX
{
public:

	ZGESVX()
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

    void applyInverse(const LapackMatrixCmplx16& A,std::vector <std::complex<double>>& b)
	{
    	    assert(A.sizeCheck(A.cols,(long)b.size()));
	        double* bptr =  &(reinterpret_cast<double(&)[2]>(b[0])[0]);
			applyInverse(A,bptr);
	}

    void applyInverse(const LapackMatrixCmplx16& A,LapackMatrixCmplx16& b)
	{
    	    assert(A.sizeCheck(A.cols,b.rows));
    		applyInverse(A,b.mData.dataPtr,b.cols);
	}

	void applyInverse(const LapackMatrixCmplx16& A, double* b, long NRHS = 1)
	{
		char FACT  = 'E'; // Equilibrate, then factor
		char TRANS = 'N'; // No transpose
		long N     = A.rows;

		//
		// Duplicate input matrix (since this zgbsvx overwrites input matrix)
		// and create temporaries

		this->A.initialize(A);
		this->AF.initialize(N,N);

		double* Aptr  =  A.mData.dataPtr;
		double* AFptr = AF.mData.dataPtr;

		long LDA   = N;
		long LDAF  = N;

		std::vector <long >   IPIV(N);
		long* IPIVptr = &IPIV[0];

		char  EQED;


		std::vector<double>   R(N);
		double* Rptr  = &R[0];

		std::vector<double>    C(N);
		double* Cptr  =  &C[0];

		std::vector<double>   B(2*N*NRHS);
		double* Bptr  =      &B[0];
	    long LDB      =    N;


		// b will be overwritten with the solution
	    // so no need to declare X separately

		double* Xptr = b;
		long LDX     = N;

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


		zgesvx_(&FACT, &TRANS, &N, &NRHS, Aptr, &LDA, AFptr, &LDAF, IPIVptr,
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



    double          RCOND;
	std::vector<double>    FERR;
	std::vector<double>    BERR;

	LapackMatrixCmplx16       A;
	LapackMatrixCmplx16      AF;

};


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class  ZHPEVX : Created for computing eigensystem components
//                 of a complex matrix..
// LAPACK base routine description:
// ZHPEVX computes selected eigenvalues and, optionally, eigenvectors
// of a complex Hermitian matrix A in packed storage.
// Eigenvalues/vectors can be selected by specifying either a range of
// values or a range of indices for the desired eigenvalues.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class ZHPEVX
{
public :

	ZHPEVX(){}

	void initialize()
	{
	AP.initialize();
    WORK.initialize();
    RWORK.clear();
    IWORK.clear();
    IFAIL.clear();
	}

	// Computes the eigCount algebraically smallest eigenvalues and eigenvectors.
	// The value returned is the number of eigenvalues found.

	long createAlgSmallestEigensystem(long eigCount, SCC::LapackMatrixCmplx16& A, std::vector<double>& eigValues,
			                          SCC::LapackMatrixCmplx16& eigVectors)
	{
    if(A.getRowDimension() != A.getColDimension())
	{
			throw std::runtime_error("\nZHPEVX : Non-square matrix input argument  \n");
	}

    long N = A.getRowDimension();

    if(eigCount > N)
    {
    	std::stringstream sout;
    	sout << "\nZHPEVX Error \n";
    	sout << "Requested number of eigenvalues/eigenvectors exceeds system dimension. \n";
    	throw std::runtime_error(sout.str());
    }

    char JOBZ   = 'V'; // Specify N for eigenvalues only
    char RANGE  = 'I'; // Specify index range of eigenvalues to be find (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    AP.initialize(A.createUpperTriPacked());

    double VL = 0;
    double VU = 0;

    long IL = 1;        // Index of smallest eigenvalue returned
    long IU = eigCount; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long M = 0;                            // Number of eigenvalues output

    eigValues.clear();                     // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ   = N;
    long Mstar = (IU-IL) + 1;              // Maximal number of eigenvalues to be computed when using index specification

    eigVectors.initialize(LDZ,Mstar);      // Matrix whose columns containing the eigenvectors (Z in original call)

    long INFO = 0;

    WORK.initialize();
    RWORK.clear();
    IWORK.clear();
    IFAIL.clear();

    WORK.initialize(2*N,1);
    RWORK.resize(7*N,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);


    zhpevx_(&JOBZ, &RANGE, &UPLO,&N,AP.mData.getDataPointer(),&VL,&VU,&IL,&IU,&ABSTOL,&M,eigValues.data(),
    eigVectors.mData.getDataPointer(),&LDZ,WORK.mData.getDataPointer(),RWORK.data(),IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nZHPEVX \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }

    // resize the eig values array to the number of eigenvalues found

    eigValues.resize(M);
    return M;
	}
	
	
    // Computes the eigCount algebraically smallest eigenvalues and eigenvectors.
	// The value returned is the number of eigenvalues found.

	long createEigensystem(SCC::LapackMatrixCmplx16& A, std::vector<double>& eigValues, SCC::LapackMatrixCmplx16& eigVectors)
	{
    if(A.getRowDimension() != A.getColDimension())
	{
			throw std::runtime_error("\nZHPEVX : Non-square matrix input argument  \n");
	}

    long N = A.getRowDimension();

    char JOBZ   = 'V'; // Specify N for eigenvalues only
    char RANGE  = 'A'; // Specify index range of eigenvalues to be find (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    AP.initialize(A.createUpperTriPacked());

    double VL = 0;
    double VU = 0;

    long IL = 1; // Index of smallest eigenvalue returned
    long IU = N; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long M = 0;                            // Number of eigenvalues output

    eigValues.clear();                     // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ   = N;
    long Mstar = N;                       // Maximal number of eigenvalues to be computed when using index specification

    eigVectors.initialize(LDZ,Mstar);      // Matrix whose columns containing the eigenvectors (Z in original call)

    long INFO = 0;

    WORK.initialize();
    RWORK.clear();
    IWORK.clear();
    IFAIL.clear();

    WORK.initialize(2*N,1);
    RWORK.resize(7*N,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);


    zhpevx_(&JOBZ, &RANGE, &UPLO,&N,AP.mData.getDataPointer(),&VL,&VU,&IL,&IU,&ABSTOL,&M,eigValues.data(),
    eigVectors.mData.getDataPointer(),&LDZ,WORK.mData.getDataPointer(),RWORK.data(),IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nZHPEVX \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }

    // resize the eig values array to the number of eigenvalues found

    eigValues.resize(M);
    return M;
	}

    long createAlgSmallestEigenvalues(long eigCount, SCC::LapackMatrixCmplx16& A, std::vector<double>& eigValues)
	{
    if(A.getRowDimension() != A.getColDimension())
	{
	    throw std::runtime_error("\nZHPEVX : Non-square matrix input argument  \n");
	}

    long N = A.getRowDimension();

    if(eigCount > N)
    {
    	std::stringstream sout;
    	sout << "\nZHPEVX Error \n";
    	sout << "Requested number of eigenvalues exceeds system dimension. \n";
    	throw std::runtime_error(sout.str());
    }

    char JOBZ   = 'N'; // Specify N for eigenvalues only
    char RANGE  = 'I'; // Specify index range of eigenvalues to be find (A for all, V for interval)
    char UPLO   = 'U'; // Store complex Hermetian matrix in upper trianglar packed form

    AP.initialize(A.createUpperTriPacked());

    double VL = 0;
    double VU = 0;

    long IL = 1;        // Index of smallest eigenvalue returned
    long IU = eigCount; // Index of largest  eigenvalue returned

    char   DLAMCH_IN = 'S';
    double ABSTOL    =  2.0*(dlamch_(&DLAMCH_IN));

    long M = 0;                     // Number of eigenvalues output

    eigValues.clear();              // W parameter in original call
    eigValues.resize(N,0.0);

    long LDZ     = 1;
    double Zdata = 0.0;              // Double value to be used as reference for null eigenvector

    long INFO = 0;

    WORK.initialize();
    RWORK.clear();
    IWORK.clear();
    IFAIL.clear();

    WORK.initialize(2*N,1);
    RWORK.resize(7*N,0.0);
    IWORK.resize(5*N,0);
    IFAIL.resize(N,0);


    zhpevx_(&JOBZ, &RANGE, &UPLO,&N,AP.mData.getDataPointer(),&VL,&VU,&IL,&IU,&ABSTOL,&M,eigValues.data(),
    &Zdata,&LDZ,WORK.mData.getDataPointer(),RWORK.data(),IWORK.data(),IFAIL.data(),&INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nZHPEVX \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }

    eigValues.resize(M);
    return M;
	}


	SCC::LapackMatrixCmplx16      AP; // For storing packed matrix in packed Hermitian form

    SCC::LapackMatrixCmplx16    WORK;
    std::vector<double>        RWORK;
    std::vector<long>          IWORK;
    std::vector<long>          IFAIL;
};





//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class  ZGEEVX  : Created for computing eigensystem components
//                 of a complex matrix..
// LAPACK base routine description:
// ZGEEVX computes for an N-by-N complex general matrix A, the
// eigenvalues and, optionally, the left and/or right eigenvectors.
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

class ZGEEVX 
{
public :

	ZGEEVX (){}

	void initialize()
	{
	AS.initialize();
    WORK.clear();
    RWORK.clear();
    SCALE.clear();
    RCONDE.clear();
    RCONDV.clear();
	}

	
    // Computes the eigenvalues and left and right eigenvectors.

	void createEigensystem(SCC::LapackMatrixCmplx16& A,
	std::vector<std::complex<double>>& eigValues,
	SCC::LapackMatrixCmplx16& eigVecLeft,
	SCC::LapackMatrixCmplx16& eigVecRight)
	{
    if(A.getRowDimension() != A.getColDimension())
	{
			throw std::runtime_error("\nZHPEVX : Non-square matrix input argument  \n");
	}

	AS.initialize(A); // Copy input matrix since scaling and permuting

	char BALANC = 'B';  // Both diagonally scale and permute the copy of A.
    char JOBVL  = 'V'; 	// Compute left eigenvectors
    char JOBVR  = 'V'; 	// Compute right eigenvectors
    char SENSE  = 'B';  // Compute left and right eigenvector reciprocal condition numbers
    long N      = A.getRowDimension();
    long LDA    = N;

    eigValues.resize(N,{0.0,0.0});
    double* eigValuePtr = reinterpret_cast<double*>(const_cast< std::complex<double>* >(&eigValues[0]));

    eigVecLeft.initialize(N,N);
    eigVecRight.initialize(N,N);

    double* VLptr = eigVecLeft.mData.getDataPointer();
    long LDVL     = N;

    double* VRptr = eigVecRight.mData.getDataPointer();
    long LDVR     = N;

    long ILO      = 0;
    long IHI      = 0;

    SCALE.resize(N,0.0);

    double ABNRM;

    RCONDE.resize(N,0.0);
    RCONDV.resize(N,0.0);

    RWORK.resize(2*N,0.0);
    long INFO = 0;

    long LWORK = N*N + 3*N;
    WORK.resize(2*LWORK); // 2 X LWORK since need complex*16


    zgeevx_(&BALANC,&JOBVL, &JOBVR, &SENSE,&N,AS.mData.getDataPointer(), &LDA,
        eigValuePtr, VLptr, &LDVL, VRptr,&LDVR, &ILO, &IHI, &SCALE[0], &ABNRM,
        &RCONDE[0], &RCONDV[0], &WORK[0],&LWORK, &RWORK[0], &INFO);

    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nZGEEVX \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }

	}




	SCC::LapackMatrixCmplx16      AS; // Copy of input matrix (S for scaled)

    std::vector<double> WORK;
    std::vector<double> RWORK;
    std::vector<double> SCALE;
    std::vector<double> RCONDE;
    std::vector<double> RCONDV;
};


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Class   ZGEESX: Created for computing Schur decomposition of a general
//                 complex matrix.
// LAPACK base routine description:
// ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the
// eigenvalues, the Schur form T, and, optionally, the matrix of Schur
// vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).

// Optionally, it also orders the eigenvalues on the diagonal of the
// Schur form so that selected eigenvalues are at the top left;
// computes a reciprocal condition number for the average of the
// selected eigenvalues (RCONDE); and computes a reciprocal condition
// number for the right invariant subspace corresponding to the
// selected eigenvalues (RCONDV).  The leading columns of Z form an
// orthonormal basis for this invariant subspace.
//
// Note: Using ints for LOGICAL arguments instead of C++ bool. Since
// values are passed by pointer, this will support LOGICAL values
// used by Fortran declared LOGICAL(1) through LOGICAL(4).
//
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

//
// Required functions for ZGEESX that are
// declared external in Fortran and passed in as
// arguments to the zgeesx function call.
//

//
// Selection based on eigenvalues with positive 
// (including 0) real part. 
//
extern "C" int eigSelectRealPos(double* C)
{
    if(C[0] >= 0.0) {return true;} return false;
}

//
// Selection based on eigenvalues negative real part
//
extern "C" int eigSelectRealNeg(double* C)
{
    if(C[0] < 0.0) {return true;} return false;
}


class ZGEESX
{

public:

enum {NONE, SORT_POS, SORT_NEG};

//
// The input matrix is not altered, on output Q contains the Schur
// vectors and T the upper triangular decomposition of A with
// eigenvalues on the diagonal.
//
// Use sortType = SORT_POS to select eigenvalues with positive real
//                part to be located in upper left quadrant of
//                decomposition.
//
// Use sortType = SORT_NEG to select eigenvalues with negative real
//                part to be located in upper left quadrant of
//                decomposition.
//
// The first sortedDim columns of Q form an invariant subspace
// associated with the eigenvalues selected by sortType.
//
//
void computeSchurDecomposition(const SCC::LapackMatrixCmplx16& A,
SCC::LapackMatrixCmplx16& Q, SCC::LapackMatrixCmplx16& T,
std::vector<std::complex<double>>& eigValues,
long& sortedDim, int sortType = ZGEESX::NONE)
{
    long N        =  A.rows;

    char JOBVS    = 'V'; // Compute Schur vectors

    char SORTFLAG = 'S'; // Sort eigenvalues

    if(sortType == ZGEESX::NONE){SORTFLAG = 'N';}

    char SENSE    = 'B';
    if(sortType == ZGEESX::NONE)
    {
    SENSE     = 'N';
    RCONDE    = 0.0; // Infinite condition numbers if no sorting
    RCONDV    = 0.0; // Infinite condition numbers if no sorting
    sortedDim = 0;
    }

    long NSIZE    =  N;

    //
    // Copy in put A to T, zgeesx will overwrite T with
    // the upper triangular component of the decomposition
    //

    T.initialize(A);
    double*Aptr  = T.mData.getDataPointer();

    long LDA     = N;

    SCC::LapackMatrixCmplx16 cEigVals(N,1);

    Q.initialize(N,N);

    long LDVS = N;

    long LWORK = (N*(N+1))/2;
    std::vector<double> cWORK(2*LWORK,0.0);
    std::vector<double> RWORK(N,0.0);

    std::vector<int> BWORK(N,0); // Using int's for logical, possibly over-allocating memory
                                 // but certainly sufficient.

    long INFO = 0;

    if(sortType == ZGEESX::SORT_POS)
    {
    zgeesx_(&JOBVS,&SORTFLAG,eigSelectRealPos,&SENSE,&NSIZE,Aptr,&LDA,&sortedDim,
           cEigVals.mData.getDataPointer(),Q.mData.getDataPointer(),&LDVS,
           &RCONDE,&RCONDV,&cWORK[0],&LWORK,&RWORK[0],&BWORK[0],&INFO);
    }
    else if(sortType == ZGEESX::SORT_NEG)
    {
    zgeesx_(&JOBVS,&SORTFLAG,eigSelectRealNeg,&SENSE,&NSIZE,Aptr,&LDA,&sortedDim,
           cEigVals.mData.getDataPointer(),Q.mData.getDataPointer(),&LDVS,
           &RCONDE,&RCONDV,&cWORK[0],&LWORK,&RWORK[0],&BWORK[0],&INFO);
    }
    else
    {
    	zgeesx_(&JOBVS,&SORTFLAG,eigSelectRealPos,&SENSE,&NSIZE,Aptr,&LDA,&sortedDim,
        cEigVals.mData.getDataPointer(),Q.mData.getDataPointer(),&LDVS,
        &RCONDE,&RCONDV,&cWORK[0],&LWORK,&RWORK[0],&BWORK[0],&INFO);
    }

    eigValues.resize(N);
    for(long k = 0; k < N; k++)
    {
    eigValues[k] = cEigVals(k);
    }


    if(INFO != 0)
    {
    	std::stringstream sout;
    	sout << "\nZGEESX \nError INFO = " << INFO << "\n";
    	throw std::runtime_error(sout.str());
    }
}

// A utility for creating a matrix whose columns are those of the
// invariant subspace vectors associated with those eigenvalues
// obtained by the sorting type. The input matrix Q, and the value
// of sortedDim are those obtained from an invocation of
// computeSchurDecomposition.
//
//
void createInvariantSubpace(long sortedDim, const SCC::LapackMatrixCmplx16& Q, SCC::LapackMatrixCmplx16& V)
{
	long rows = Q.rows;
	long cols = sortedDim;

    if(sortedDim == 0)
    {
    V.initialize();
    return;
    }

	V.initialize(rows,cols);
	for(long j = 0; j < cols; j++)
	{
	for(long i = 0; i < rows; i++)
	{
	V(i,j) = Q(i,j);
	}}

}

//
//  public class variables to allow access to LAPACK return arguments
//
    double RCONDE;
    double RCONDV;

};


} // Namespace SCC

//
// LAPACK documentation
//
/////////////////////////////////////////////////////////////////////////////
// ZGESVX
/////////////////////////////////////////////////////////////////////////////
/*
        subroutine zgesvx 	( 	character  	fact,
		character  	trans,
		integer  	n,
		integer  	nrhs,
		complex*16, dimension( lda, * )  	a,
		integer  	lda,
		complex*16, dimension( ldaf, * )  	af,
		integer  	ldaf,
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

     ZGESVX uses the LU factorization to compute the solution to a complex
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
    [in]	FACT

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

    [in]	TRANS

              TRANS is CHARACTER*1
              Specifies the form of the system of equations:
              = 'N':  A * X = B     (No transpose)
              = 'T':  A**T * X = B  (Transpose)
              = 'C':  A**H * X = B  (Conjugate transpose)

    [in]	N

              N is INTEGER
              The number of linear equations, i.e., the order of the
              matrix A.  N >= 0.

    [in]	NRHS

              NRHS is INTEGER
              The number of right hand sides, i.e., the number of columns
              of the matrices B and X.  NRHS >= 0.

    [in,out]	A

              A is COMPLEX*16 array, dimension (LDA,N)
              On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
              not 'N', then A must have been equilibrated by the scaling
              factors in R and/or C.  A is not modified if FACT = 'F' or
              'N', or if FACT = 'E' and EQUED = 'N' on exit.

              On exit, if EQUED .ne. 'N', A is scaled as follows:
              EQUED = 'R':  A := diag(R) * A
              EQUED = 'C':  A := A * diag(C)
              EQUED = 'B':  A := diag(R) * A * diag(C).

    [in]	LDA

              LDA is INTEGER
              The leading dimension of the array A.  LDA >= max(1,N).

    [in,out]	AF

              AF is COMPLEX*16 array, dimension (LDAF,N)
              If FACT = 'F', then AF is an input argument and on entry
              contains the factors L and U from the factorization
              A = P*L*U as computed by ZGETRF.  If EQUED .ne. 'N', then
              AF is the factored form of the equilibrated matrix A.

              If FACT = 'N', then AF is an output argument and on exit
              returns the factors L and U from the factorization A = P*L*U
              of the original matrix A.

              If FACT = 'E', then AF is an output argument and on exit
              returns the factors L and U from the factorization A = P*L*U
              of the equilibrated matrix A (see the description of A for
              the form of the equilibrated matrix).

    [in]	LDAF

              LDAF is INTEGER
              The leading dimension of the array AF.  LDAF >= max(1,N).

    [in,out]	IPIV

              IPIV is INTEGER array, dimension (N)
              If FACT = 'F', then IPIV is an input argument and on entry
              contains the pivot indices from the factorization A = P*L*U
              as computed by ZGETRF; row i of the matrix was interchanged
              with row IPIV(i).

              If FACT = 'N', then IPIV is an output argument and on exit
              contains the pivot indices from the factorization A = P*L*U
              of the original matrix A.

              If FACT = 'E', then IPIV is an output argument and on exit
              contains the pivot indices from the factorization A = P*L*U
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
              On entry, the N-by-NRHS right hand side matrix B.
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

              RWORK is DOUBLE PRECISION array, dimension (MAX(1,2*N))
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
// ZHPEVX
/////////////////////////////////////////////////////////////////////////////

/*
ZHPEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for OTHER matrices

subroutine zhpevx	(	character 	jobz,
character 	range,
character 	uplo,
integer 	n,
complex*16, dimension( * ) 	ap,
double precision 	vl,
double precision 	vu,
integer 	il,
integer 	iu,
double precision 	abstol,
integer 	m,
double precision, dimension( * ) 	w,
complex*16, dimension( ldz, * ) 	z,
integer 	ldz,
complex*16, dimension( * ) 	work,
double precision, dimension( * ) 	rwork,
integer, dimension( * ) 	iwork,
integer, dimension( * ) 	ifail,
integer 	info
)

Purpose:
 ZHPEVX computes selected eigenvalues and, optionally, eigenvectors
 of a complex Hermitian matrix A in packed storage.
 Eigenvalues/vectors can be selected by specifying either a range of
 values or a range of indices for the desired eigenvalues.
 
Parameters
[in]	JOBZ
          JOBZ is CHARACTER*1
          = 'N':  Compute eigenvalues only;
          = 'V':  Compute eigenvalues and eigenvectors.
[in]	RANGE
          RANGE is CHARACTER*1
          = 'A': all eigenvalues will be found;
          = 'V': all eigenvalues in the half-open interval (VL,VU]
                 will be found;
          = 'I': the IL-th through IU-th eigenvalues will be found.
[in]	UPLO
          UPLO is CHARACTER*1
          = 'U':  Upper triangle of A is stored;
          = 'L':  Lower triangle of A is stored.
[in]	N
          N is INTEGER
          The order of the matrix A.  N >= 0.
[in,out]	AP
          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
          On entry, the upper or lower triangle of the Hermitian matrix
          A, packed columnwise in a linear array.  The j-th column of A
          is stored in the array AP as follows:
          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.

          On exit, AP is overwritten by values generated during the
          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
          and first superdiagonal of the tridiagonal matrix T overwrite
          the corresponding elements of A, and if UPLO = 'L', the
          diagonal and first subdiagonal of T overwrite the
          corresponding elements of A.
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

          where EPS is the machine precision.  If ABSTOL is less than
          or equal to zero, then  EPS*|T|  will be used in its place,
          where |T| is the 1-norm of the tridiagonal matrix obtained
          by reducing AP to tridiagonal form.

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
          If INFO = 0, the selected eigenvalues in ascending order.
[out]	Z
          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))
          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
          contain the orthonormal eigenvectors of the matrix A
          corresponding to the selected eigenvalues, with the i-th
          column of Z holding the eigenvector associated with W(i).
          If an eigenvector fails to converge, then that column of Z
          contains the latest approximation to the eigenvector, and
          the index of the eigenvector is returned in IFAIL.
          If JOBZ = 'N', then Z is not referenced.
          Note: the user must ensure that at least max(1,M) columns are
          supplied in the array Z; if RANGE = 'V', the exact value of M
          is not known in advance and an upper bound must be used.
[in]	LDZ
          LDZ is INTEGER
          The leading dimension of the array Z.  LDZ >= 1, and if
          JOBZ = 'V', LDZ >= max(1,N).
[out]	WORK
          WORK is COMPLEX*16 array, dimension (2*N)
[out]	RWORK
          RWORK is DOUBLE PRECISION array, dimension (7*N)
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
/*

ZGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE complex matrices

zgeevx()
subroutine zgeevx	(	character 	balanc,
character 	jobvl,
character 	jobvr,
character 	sense,
integer 	n,
complex*16, dimension( lda, * ) 	a,
integer 	lda,
complex*16, dimension( * ) 	w,
complex*16, dimension( ldvl, * ) 	vl,
integer 	ldvl,
complex*16, dimension( ldvr, * ) 	vr,
integer 	ldvr,
integer 	ilo,
integer 	ihi,
double precision, dimension( * ) 	scale,
double precision 	abnrm,
double precision, dimension( * ) 	rconde,
double precision, dimension( * ) 	rcondv,
complex*16, dimension( * ) 	work,
integer 	lwork,
double precision, dimension( * ) 	rwork,
integer 	info 
)		
 
Purpose:
 ZGEEVX computes for an N-by-N complex nonsymmetric matrix A, the
 eigenvalues and, optionally, the left and/or right eigenvectors.

 Optionally also, it computes a balancing transformation to improve
 the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
 SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
 (RCONDE), and reciprocal condition numbers for the right
 eigenvectors (RCONDV).

 The right eigenvector v(j) of A satisfies
                  A * v(j) = lambda(j) * v(j)
 where lambda(j) is its eigenvalue.
 The left eigenvector u(j) of A satisfies
               u(j)**H * A = lambda(j) * u(j)**H
 where u(j)**H denotes the conjugate transpose of u(j).

 The computed eigenvectors are normalized to have Euclidean norm
 equal to 1 and largest component real.

 Balancing a matrix means permuting the rows and columns to make it
 more nearly upper triangular, and applying a diagonal similarity
 transformation D * A * D**(-1), where D is a diagonal matrix, to
 make its rows and columns closer in norm and the condition numbers
 of its eigenvalues and eigenvectors smaller.  The computed
 reciprocal condition numbers correspond to the balanced matrix.
 Permuting rows and columns will not change the condition numbers
 (in exact arithmetic) but diagonal scaling will.  For further
 explanation of balancing, see section 4.10.2 of the LAPACK
 Users' Guide.
Parameters
[in]	BALANC	
          BALANC is CHARACTER*1
          Indicates how the input matrix should be diagonally scaled
          and/or permuted to improve the conditioning of its
          eigenvalues.
          = 'N': Do not diagonally scale or permute;
          = 'P': Perform permutations to make the matrix more nearly
                 upper triangular. Do not diagonally scale;
          = 'S': Diagonally scale the matrix, ie. replace A by
                 D*A*D**(-1), where D is a diagonal matrix chosen
                 to make the rows and columns of A more equal in
                 norm. Do not permute;
          = 'B': Both diagonally scale and permute A.

          Computed reciprocal condition numbers will be for the matrix
          after balancing and/or permuting. Permuting does not change
          condition numbers (in exact arithmetic), but balancing does.
[in]	JOBVL	
          JOBVL is CHARACTER*1
          = 'N': left eigenvectors of A are not computed;
          = 'V': left eigenvectors of A are computed.
          If SENSE = 'E' or 'B', JOBVL must = 'V'.
[in]	JOBVR	
          JOBVR is CHARACTER*1
          = 'N': right eigenvectors of A are not computed;
          = 'V': right eigenvectors of A are computed.
          If SENSE = 'E' or 'B', JOBVR must = 'V'.
[in]	SENSE	
          SENSE is CHARACTER*1
          Determines which reciprocal condition numbers are computed.
          = 'N': None are computed;
          = 'E': Computed for eigenvalues only;
          = 'V': Computed for right eigenvectors only;
          = 'B': Computed for eigenvalues and right eigenvectors.

          If SENSE = 'E' or 'B', both left and right eigenvectors
          must also be computed (JOBVL = 'V' and JOBVR = 'V').
[in]	N	
          N is INTEGER
          The order of the matrix A. N >= 0.
[in,out]	A	
          A is COMPLEX*16 array, dimension (LDA,N)
          On entry, the N-by-N matrix A.
          On exit, A has been overwritten.  If JOBVL = 'V' or
          JOBVR = 'V', A contains the Schur form of the balanced
          version of the matrix A.
[in]	LDA	
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[out]	W	
          W is COMPLEX*16 array, dimension (N)
          W contains the computed eigenvalues.
[out]	VL	
          VL is COMPLEX*16 array, dimension (LDVL,N)
          If JOBVL = 'V', the left eigenvectors u(j) are stored one
          after another in the columns of VL, in the same order
          as their eigenvalues.
          If JOBVL = 'N', VL is not referenced.
          u(j) = VL(:,j), the j-th column of VL.
[in]	LDVL	
          LDVL is INTEGER
          The leading dimension of the array VL.  LDVL >= 1; if
          JOBVL = 'V', LDVL >= N.
[out]	VR	
          VR is COMPLEX*16 array, dimension (LDVR,N)
          If JOBVR = 'V', the right eigenvectors v(j) are stored one
          after another in the columns of VR, in the same order
          as their eigenvalues.
          If JOBVR = 'N', VR is not referenced.
          v(j) = VR(:,j), the j-th column of VR.
[in]	LDVR	
          LDVR is INTEGER
          The leading dimension of the array VR.  LDVR >= 1; if
          JOBVR = 'V', LDVR >= N.
[out]	ILO	
          ILO is INTEGER
[out]	IHI	
          IHI is INTEGER
          ILO and IHI are integer values determined when A was
          balanced.  The balanced A(i,j) = 0 if I > J and
          J = 1,...,ILO-1 or I = IHI+1,...,N.
[out]	SCALE	
          SCALE is DOUBLE PRECISION array, dimension (N)
          Details of the permutations and scaling factors applied
          when balancing A.  If P(j) is the index of the row and column
          interchanged with row and column j, and D(j) is the scaling
          factor applied to row and column j, then
          SCALE(J) = P(J),    for J = 1,...,ILO-1
                   = D(J),    for J = ILO,...,IHI
                   = P(J)     for J = IHI+1,...,N.
          The order in which the interchanges are made is N to IHI+1,
          then 1 to ILO-1.
[out]	ABNRM	
          ABNRM is DOUBLE PRECISION
          The one-norm of the balanced matrix (the maximum
          of the sum of absolute values of elements of any column).
[out]	RCONDE	
          RCONDE is DOUBLE PRECISION array, dimension (N)
          RCONDE(j) is the reciprocal condition number of the j-th
          eigenvalue.
[out]	RCONDV	
          RCONDV is DOUBLE PRECISION array, dimension (N)
          RCONDV(j) is the reciprocal condition number of the j-th
          right eigenvector.
[out]	WORK	
          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]	LWORK	
          LWORK is INTEGER
          The dimension of the array WORK.  If SENSE = 'N' or 'E',
          LWORK >= max(1,2*N), and if SENSE = 'V' or 'B',
          LWORK >= N*N+2*N.
          For good performance, LWORK must generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates the optimal size of the WORK array, returns
          this value as the first entry of the WORK array, and no error
          message related to LWORK is issued by XERBLA.
[out]	RWORK	
          RWORK is DOUBLE PRECISION array, dimension (2*N)
[out]	INFO	
          INFO is INTEGER
          = 0:  successful exit
          < 0:  if INFO = -i, the i-th argument had an illegal value.
          > 0:  if INFO = i, the QR algorithm failed to compute all the
                eigenvalues, and no eigenvectors or condition numbers
                have been computed; elements 1:ILO-1 and i+1:N of W
                contain eigenvalues which have converged.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

/////////////////////////////////////////////////////////////////////////////
// ZGEESX
/////////////////////////////////////////////////////////////////////////////
/*
ZGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors general complex matrices

subroutine zgeesx	(	character 	jobvs,
character 	sort,
external 	select,
character 	sense,
integer 	n,
complex*16, dimension( lda, * ) 	a,
integer 	lda,
integer 	sdim,
complex*16, dimension( * ) 	w,
complex*16, dimension( ldvs, * ) 	vs,
integer 	ldvs,
double precision 	rconde,
double precision 	rcondv,
complex*16, dimension( * ) 	work,
integer 	lwork,
double precision, dimension( * ) 	rwork,
logical, dimension( * ) 	bwork,
integer 	info
)

 ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the
 eigenvalues, the Schur form T, and, optionally, the matrix of Schur
 vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).

 Optionally, it also orders the eigenvalues on the diagonal of the
 Schur form so that selected eigenvalues are at the top left;
 computes a reciprocal condition number for the average of the
 selected eigenvalues (RCONDE); and computes a reciprocal condition
 number for the right invariant subspace corresponding to the
 selected eigenvalues (RCONDV).  The leading columns of Z form an
 orthonormal basis for this invariant subspace.

 For further explanation of the reciprocal condition numbers RCONDE
 and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where
 these quantities are called s and sep respectively).

 A complex matrix is in Schur form if it is upper triangular.
Parameters
[in]	JOBVS
          JOBVS is CHARACTER*1
          = 'N': Schur vectors are not computed;
          = 'V': Schur vectors are computed.
[in]	SORT
          SORT is CHARACTER*1
          Specifies whether or not to order the eigenvalues on the
          diagonal of the Schur form.
          = 'N': Eigenvalues are not ordered;
          = 'S': Eigenvalues are ordered (see SELECT).
[in]	SELECT
          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument
          SELECT must be declared EXTERNAL in the calling subroutine.
          If SORT = 'S', SELECT is used to select eigenvalues to order
          to the top left of the Schur form.
          If SORT = 'N', SELECT is not referenced.
          An eigenvalue W(j) is selected if SELECT(W(j)) is true.
[in]	SENSE
          SENSE is CHARACTER*1
          Determines which reciprocal condition numbers are computed.
          = 'N': None are computed;
          = 'E': Computed for average of selected eigenvalues only;
          = 'V': Computed for selected right invariant subspace only;
          = 'B': Computed for both.
          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.
[in]	N
          N is INTEGER
          The order of the matrix A. N >= 0.
[in,out]	A
          A is COMPLEX*16 array, dimension (LDA, N)
          On entry, the N-by-N matrix A.
          On exit, A is overwritten by its Schur form T.
[in]	LDA
          LDA is INTEGER
          The leading dimension of the array A.  LDA >= max(1,N).
[out]	SDIM
          SDIM is INTEGER
          If SORT = 'N', SDIM = 0.
          If SORT = 'S', SDIM = number of eigenvalues for which
                         SELECT is true.
[out]	W
          W is COMPLEX*16 array, dimension (N)
          W contains the computed eigenvalues, in the same order
          that they appear on the diagonal of the output Schur form T.
[out]	VS
          VS is COMPLEX*16 array, dimension (LDVS,N)
          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
          vectors.
          If JOBVS = 'N', VS is not referenced.
[in]	LDVS
          LDVS is INTEGER
          The leading dimension of the array VS.  LDVS >= 1, and if
          JOBVS = 'V', LDVS >= N.
[out]	RCONDE
          RCONDE is DOUBLE PRECISION
          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
          condition number for the average of the selected eigenvalues.
          Not referenced if SENSE = 'N' or 'V'.
[out]	RCONDV
          RCONDV is DOUBLE PRECISION
          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
          condition number for the selected right invariant subspace.
          Not referenced if SENSE = 'N' or 'E'.
[out]	WORK
          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
[in]	LWORK
          LWORK is INTEGER
          The dimension of the array WORK.  LWORK >= max(1,2*N).
          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM),
          where SDIM is the number of selected eigenvalues computed by
          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also
          that an error is only returned if LWORK < max(1,2*N), but if
          SENSE = 'E' or 'V' or 'B' this may not be large enough.
          For good performance, LWORK must generally be larger.

          If LWORK = -1, then a workspace query is assumed; the routine
          only calculates upper bound on the optimal size of the
          array WORK, returns this value as the first entry of the WORK
          array, and no error message related to LWORK is issued by
          XERBLA.
[out]	RWORK
          RWORK is DOUBLE PRECISION array, dimension (N)
[out]	BWORK
          BWORK is LOGICAL array, dimension (N)
          Not referenced if SORT = 'N'.
[out]	INFO
          INFO is INTEGER
          = 0: successful exit
          < 0: if INFO = -i, the i-th argument had an illegal value.
          > 0: if INFO = i, and i is
             <= N: the QR algorithm failed to compute all the
                   eigenvalues; elements 1:ILO-1 and i+1:N of W
                   contain those eigenvalues which have converged; if
                   JOBVS = 'V', VS contains the transformation which
                   reduces A to its partially converged Schur form.
             = N+1: the eigenvalues could not be reordered because some
                   eigenvalues were too close to separate (the problem
                   is very ill-conditioned);
             = N+2: after reordering, roundoff changed values of some
                   complex eigenvalues so that leading eigenvalues in
                   the Schur form no longer satisfy SELECT=.TRUE.  This
                   could also be caused by underflow due to scaling.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
*/

#endif

