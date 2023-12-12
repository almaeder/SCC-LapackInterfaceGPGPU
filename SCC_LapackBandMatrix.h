/*
 * SCC_LapackBandMatrix.h
 *
 *  Created on: Aug 17, 2018
 *      Author: anderson
 *
 *  Updated   : Dec. 9, 2023 (C.R. Anderson)
 */
//
// Instances of this class are N x N banded matrices with double entries
// that are specified with three parameters:
//
// kl = lower bandwidth
// ku = upper bandwidth
//  N = system size
//
//
// The indexing for the access operator()(i,j) starts at 0, consistent 
// with the indexing of the other LapackInterface matrix classes.
//
// When complied with the _DEBUG pre-processor directive, then index checking 
// is performed so that if an element is sought that is outside of the specification 
// of the band matrix structure and exception is triggered (using assert(...)). 
//
// Internally the matrix data is stored as a SCC::LapackMatrixCmplx16 
// instance that contains the matrix data stored using the Lapack band 
// storage convention, so that LAPACK band matrix routines can be
// invoked without having to transform the matrix data. 
//
// Typical use case consists of initializing an instance
// then setting values of the banded matrix using the
// standard access operator.
//
// The data pointer obtained using getDataPointer() can be
// passed to Lapack routines that assume a matrix stored
// using the Lapack band matrix storage convention.
//
// Bounds checking is only done if _DEBUG is defined
//
// Lapack routine dependencies : dgbmv_
/*
#############################################################################
#
# Copyright 2018- Chris Anderson
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

#ifdef  _DEBUG
#include <cstdio>
#else
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

#include <cmath>

#include "SCC_LapackHeaders.h"
#include "SCC_LapackMatrix.h"

#ifndef SCC_LAPACK_BAND_MATRIX_
#define SCC_LAPACK_BAND_MATRIX_

namespace SCC
{
class LapackBandMatrix
{
	public:

	LapackBandMatrix()
	{
	initialize();
	}

    LapackBandMatrix(const LapackBandMatrix& S)
	{
    	initialize(S);
	}

	LapackBandMatrix(long kl, long ku, long N)
	{
	    initialize(kl,ku,N);
	}

	void initialize()
	{
		mData.initialize();
		kl = 0;
		ku = 0;
    	N  = 0;
	}

	void initialize(const LapackBandMatrix& S)
	{
		kl = S.kl;
		ku = S.ku;
		N  = S.N;
	    mData.initialize(S.mData);

	}
    void initialize(long kl, long ku, long N)
	{
    	this->kl = kl;
    	this->ku = ku;
    	this->N  = N;
    	mData.initialize(kl + ku + 1, N);
	}


	#ifdef _DEBUG
	double& operator()(long i, long j)
	{
		assert(boundsCheck(i,j));
		return mData(ku +  (i-j),j);
	}

    const double& operator()(long i, long j) const
	{
    	assert(boundsCheck(i,j));
		return mData(ku +  (i-j),j);
	}




	#else
	inline double& operator()(long i, long j)
	{

		return mData(ku +  (i-j),j);
	}

    inline const double& operator()(long i, long j) const
	{
		return mData(ku +  (i-j),j);
	}
	#endif

	double* getDataPointer() const {return mData.dataPtr;}

//
// Algebraic operators utilize algebraic operations of underlying LapackMatric
//
    inline void operator=(const LapackBandMatrix& B)
	{
    	if(mData.isNull())
    	{
    		kl    = B.kl;
    		ku    = B.ku;
    		N     = B.N;
    		mData.initialize(B.mData);
    	}
    	else
    	{
    	assert(sizeCheck(B.kl, B.ku, B.N));
    	mData = B.mData;
    	}
	}


    inline void operator+=(const  LapackBandMatrix& B)
    {
        assert(sizeCheck(B.kl, B.ku, B.N));
    	mData += B.mData;
    }

    LapackBandMatrix operator+(const LapackBandMatrix& B)
    {
        assert(sizeCheck(B.kl, B.ku, B.N));
    	LapackBandMatrix     C(*this);
    	C.mData += B.mData;
        return C;
    }

    inline void operator-=(const  LapackBandMatrix& B)
    {
      assert(sizeCheck(B.kl, B.ku, B.N));
      mData -= B.mData;
    }

    LapackBandMatrix operator-(const LapackBandMatrix& B)
    {
        assert(sizeCheck(B.kl, B.ku, B.N));
    	LapackBandMatrix     C(*this);
    	C.mData -= B.mData;
        return C;
    }


    inline void operator*=(const double alpha)
    {
       mData *= alpha;
    }

    LapackBandMatrix operator*(const double alpha)
    {
    LapackBandMatrix R(*this);
    R *= alpha;
    return R;
    }

    friend LapackBandMatrix operator*(const double alpha, const LapackBandMatrix& B)
    {
    LapackBandMatrix R(B);
    R *= alpha;
    return R;
    }

    inline void operator/=(const double alpha)
    {
    mData /= alpha;
    }

    LapackBandMatrix operator/(const double alpha)
    {
    LapackBandMatrix R(*this);
    R /= alpha;
    return R;
    }


    bool isNull() const
    {
    return mData.isNull();
    }

    void setToValue(double val)
    {
      mData.setToValue(val);
    }

    void setToIdentity()
    {
    	setToValue(0.0);
    	for(long k = 0; k < N; k++)
    	{
    	this->operator()(k,k) = 1.0;
    	}
    }

    void setDiagonal(const std::vector<double>& diag)
    {
    	for(long k = 0; k < N; k++)
    	{
    		this->operator()(k,k) = diag[k];
    	}
    }

    double normFrobenius()
    {
    	double val    = 0.0;
    	double valSum = 0.0;

    	for(long i = 0; i < N; i++)
    	{
    	for(long j = std::max((long)0,i- kl); j <= std::min(i+ku,N-1); j++)
    	{
    		val = this->operator()(i,j);
    		valSum += val*val;
    	}}
    	return std::sqrt(valSum);
    }

    //
    // The only band matrix X matrix operation allowed
    // is a band matrix times a matrix with a single
    // column
    //
    LapackMatrix operator*(const LapackMatrix& x)
    {
    assert(sizecheckNx1(x.rows,x.cols));
	LapackMatrix y(x.rows,1);

    char TRANS     = 'N';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;


    /*
     DGBMV  performs one of the matrix-vector operations

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,

    where alpha and beta are scalars, x and y are vectors and A is an
    m by n band matrix, with kl sub-diagonals and ku super-diagonals.
   */

    dgbmv_(&TRANS,&N,&N, &kl, &ku,&ALPHA, mData.dataPtr,&mData.rows, x.dataPtr,&INCX,&BETA,y.dataPtr,&INCY);
	return y;
    }


    std::vector<double> operator*(const std::vector<double>& x)
    {
	std::vector<double> y(N,0.0);

    char TRANS     = 'N';
    double ALPHA   = 1.0;
    double BETA    = 0.0;
    long INCX      = 1;
    long INCY      = 1;

    /*
     DGBMV  performs one of the matrix-vector operations

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,

    where alpha and beta are scalars, x and y are vectors and A is an
    m by n band matrix, with kl sub-diagonals and ku super-diagonals.
   */

    dgbmv_(&TRANS,&N,&N, &kl, &ku,&ALPHA, mData.dataPtr,&mData.rows,const_cast<double*>(&x[0]),&INCX,&BETA,&y[0],&INCY);
	return y;
    }


/*!  Outputs the band matrix as a dense matrix with the (0,0) element in the upper left corner */

void printDense(std::ostream& outStream, int precision = 3)
{
	    double val;

        for(long i = 0;  i < N; i++)
        {
        for(long j = 0; j <  N; j++)
        {
          if((j > i + ku)||(j < i - kl))
          {val = 0.0;}
          else
          {val = this->operator()(i,j); }
          outStream <<   std::scientific << std::setprecision(precision) <<  std::right << std::setw(precision+7) << val << " ";
        }
        outStream << std::endl;
        }

}

/*
//
// ToDo: Decide on proper output format to support both << and >> stream operators
//       and then implements.
//
friend std::ostream& operator<<(std::ostream& outStream, const LapackBandMatrix&  V)
{
        long i; long j; double val;

        for(i = 0;  i < V.N; i++)
        {
        for(j = 0; j <  V.N; j++)
        {
          if((j > i + V.ku)||(j < i - V.kl))
          {val = 0.0;}
          else
          {val = V(i,j); }
          outStream <<   std::scientific << std::setprecision(3) <<  std::right << std::setw(10) << val << " ";
        }
        outStream << std::endl;
        }
        return outStream;
}
*/


// Fortran indexing bounds check max(1,j-ku) <= i <= min(N,j+kl)

#ifdef _DEBUG
	bool boundsCheck(long i, long j) const
	{
        long a = (j-ku > 0)   ?  j - ku : 0;
	    long b = (j+kl < N-1) ?  j + kl : N-1;
	    if((i< a) || (i > b))
	    {
	    std::cerr  <<  "Band matrix storage error " << std::endl;
	    std::cerr  <<  "kl =  " << kl  <<  std::endl;
	    std::cerr  <<  "ku =  " << ku  <<  std::endl;
	    std::cerr  <<  "N  =  " << N   <<  std::endl;
	    std::cerr  <<  "Offending indices " << "(" << i << ", " << j << ")" << std::endl;
	    return false;
	    }
	    return true;
	}


    bool sizeCheck(long dLower, long dUpper, long Msize) const
    {
    	if((dLower != kl)||(dUpper != ku)||(Msize != N))
    	{
    	std::cerr  <<  "Band matrix assignment error  "  << "\n";
	    std::cerr  <<  "kl =  " << kl   << "\n";
	    std::cerr  <<  "ku =  " << ku   << "\n";
	    std::cerr  <<  "N  =  " << N    << "\n";
	    std::cerr <<  "Offending sizes input  " << "kl : " << dLower << ", ku : " << dUpper  << ", N = " <<  Msize << "\n";
	    return false;
    	}
    	return true;
    }

    bool sizecheckNx1(long rows, long cols) const
    {
    if((rows != N) || (cols != 1))
    {
    std::cerr  <<  "LapackBandMatrix * LapackMatrix error   "  << "\n";
    std::cerr  <<  "LapackMatrix must be N x 1 matrix       "  << "\n";
    std::cerr  <<  "LapackMatrix rows : "  << rows << "\n";
    std::cerr  <<  "LapackMatrix cols : "  << cols << "\n";
    return false;
    }
    return true;
    }
#else
        bool boundsCheck(long, long) const {return true;}
        bool sizeCheck(long dLower, long dUpper, long Msize) const {return true;}
        bool sizecheckNx1(long rows, long cols)  const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    std::cerr << "LapackBandMatrix sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    std::cerr << "LapackBandMatrix sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif




    SCC::LapackMatrix mData;

	long ku;
	long kl;
	long  N;
};

} // Namespace SCC

//
// LAPACK documentation
//
////////////////////////////////////////////////////////////////
// DGBMV
////////////////////////////////////////////////////////////////
/*
subroutine dgbmv	(	character 	trans,
integer 	m,
integer 	n,
integer 	kl,
integer 	ku,
double precision 	alpha,
double precision, dimension(lda,*) 	a,
integer 	lda,
double precision, dimension(*) 	x,
integer 	incx,
double precision 	beta,
double precision, dimension(*) 	y,
integer 	incy
)

Purpose:
 DGBMV  performs one of the matrix-vector operations

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,

 where alpha and beta are scalars, x and y are vectors and A is an
 m by n band matrix, with kl sub-diagonals and ku super-diagonals.
Parameters
[in]	TRANS
          TRANS is CHARACTER*1
           On entry, TRANS specifies the operation to be performed as
           follows:

              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.

              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
[in]	M
          M is INTEGER
           On entry, M specifies the number of rows of the matrix A.
           M must be at least zero.
[in]	N
          N is INTEGER
           On entry, N specifies the number of columns of the matrix A.
           N must be at least zero.
[in]	KL
          KL is INTEGER
           On entry, KL specifies the number of sub-diagonals of the
           matrix A. KL must satisfy  0 .le. KL.
[in]	KU
          KU is INTEGER
           On entry, KU specifies the number of super-diagonals of the
           matrix A. KU must satisfy  0 .le. KU.
[in]	ALPHA
          ALPHA is DOUBLE PRECISION.
           On entry, ALPHA specifies the scalar alpha.
[in]	A
          A is DOUBLE PRECISION array, dimension ( LDA, N )
           Before entry, the leading ( kl + ku + 1 ) by n part of the
           array A must contain the matrix of coefficients, supplied
           column by column, with the leading diagonal of the matrix in
           row ( ku + 1 ) of the array, the first super-diagonal
           starting at position 2 in row ku, the first sub-diagonal
           starting at position 1 in row ( ku + 2 ), and so on.
           Elements in the array A that do not correspond to elements
           in the band matrix (such as the top left ku by ku triangle)
           are not referenced.
           The following program segment will transfer a band matrix
           from conventional full matrix storage to band storage:

                 DO 20, J = 1, N
                    K = KU + 1 - J
                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
                       A( K + I, J ) = matrix( I, J )
              10    CONTINUE
              20 CONTINUE
[in]	LDA
          LDA is INTEGER
           On entry, LDA specifies the first dimension of A as declared
           in the calling (sub) program. LDA must be at least
           ( kl + ku + 1 ).
[in]	X
          X is DOUBLE PRECISION array, dimension at least
           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
           and at least
           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
           Before entry, the incremented array X must contain the
           vector x.
[in]	INCX
          INCX is INTEGER
           On entry, INCX specifies the increment for the elements of
           X. INCX must not be zero.
[in]	BETA
          BETA is DOUBLE PRECISION.
           On entry, BETA specifies the scalar beta. When BETA is
           supplied as zero then Y need not be set on input.
[in,out]	Y
          Y is DOUBLE PRECISION array, dimension at least
           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
           and at least
           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
           Before entry, the incremented array Y must contain the
           vector y. On exit, Y is overwritten by the updated vector y.
           If either m or n is zero, then Y not referenced and the function
           performs a quick return.
[in]	INCY
          INCY is INTEGER
           On entry, INCY specifies the increment for the elements of
           Y. INCY must not be zero.
Author
Univ. of Tennessee
Univ. of California Berkeley
Univ. of Colorado Denver
NAG Ltd.
 */



#endif /* SCC_LapackBandMatrix_  */
