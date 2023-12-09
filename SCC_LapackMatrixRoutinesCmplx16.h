
#include "SCC_LapackHeaders.h"
#include "SCC_LapackMatrix.h"

#include "SCC_LapackMatrixCmplx16.h"

#include <complex>
#include <cassert>

// Class wrappers for Lapack routines:  ZGESVX, ZHPEVX

/*
 ZGESVX uses the LU factorization to compute the solution to a complex
 system of linear equations
    A * X = B,
 where A is an N-by-N matrix and X and B are N-by-NRHS matrices.

 Error bounds on the solution and a condition estimate are also
 provided.
*/

/*
 ZHPEVX computes selected eigenvalues and, optionally, eigenvectors
 of a complex Hermitian matrix A in packed storage.
 Eigenvalues/vectors can be selected by specifying either a range of
 values or a range of indices for the desired eigenvalues.
*/


#ifndef SCC_LAPACK_MATRIX_ROUTINES_CMPLX16
#define SCC_LAPACK_MATRIX_ROUTINES_CMPLX16


namespace SCC
{
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


} // Namespace SCC

//
// Lapack documentation
//
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
#endif

