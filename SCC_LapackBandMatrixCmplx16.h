/*
 * SCC_LapackBandMatrixCmplx16.h
 *
 *  Created on: Dec. 3,2023
 *      Author: anderson
 */
//
// Instances of this class are N x N banded matrices with complex entries
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
/*
#############################################################################
#
# Copyright 2023- Chris Anderson
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

#include "LapackInterface/SCC_LapackMatrixCmplx16.h"

#ifndef SCC_LAPACK_BAND_MATRIX_CMPLX_16
#define SCC_LAPACK_BAND_MATRIX_CMPLX_16

namespace SCC
{
class LapackBandMatrixCmplx16
{
	public:

	LapackBandMatrixCmplx16()
	{
	initialize();
	}

    LapackBandMatrixCmplx16(const LapackBandMatrixCmplx16& S)
	{
    	initialize(S);
	}

	LapackBandMatrixCmplx16(long kl, long ku, long N)
	{
	    initialize(kl,ku,N);
	}

	void initialize()
	{
		cmplxMdata.initialize();
		kl = 0;
		ku = 0;
    	N  = 0;
	}

	void initialize(const LapackBandMatrixCmplx16& S)
	{
		kl = S.kl;
		ku = S.ku;
		N  = S.N;
	    cmplxMdata.initialize(S.cmplxMdata);

	}
    void initialize(long kl, long ku, long N)
	{
    	this->kl = kl;
    	this->ku = ku;
    	this->N  = N;
    	cmplxMdata.initialize(kl + ku + 1, N);
	}


	#ifdef _DEBUG
	 std::complex<double>& operator()(long i, long j)
	{
		assert(boundsCheck(i,j));
		return cmplxMdata(ku +  (i-j),j);
	}

    const  std::complex<double>& operator()(long i, long j) const
	{
    	assert(boundsCheck(i,j));
		return cmplxMdata(ku +  (i-j),j);
	}


	#else
	inline std::complex<double>& operator()(long i, long j)
	{
		return cmplxMdata(ku +  (i-j),j);
	}

    inline const std::complex<double>& operator()(long i, long j) const
	{
		return cmplxMdata(ku +  (i-j),j);
	}
	#endif

	double* getDataPointer() const {return cmplxMdata.mData.dataPtr;}

    void setToValue(double val)
	{
		cmplxMdata.setToValue(val);
	}

	void setToValue(const std::complex<double>& val)
	{
	cmplxMdata.setToValue(val);
	}

	bool isNull()
	{
	return cmplxMdata.isNull();
	}


    double normFrobenius()
    {
    	double valSum = 0.0;

    	for(long i = 0; i < N; i++)
    	{
    	for(long j = std::max((long)0,i- kl); j <= std::min(i+ku,N-1); j++)
    	{
    	    valSum += std::norm(this->operator()(i,j));
    	}}
    	return std::sqrt(valSum);
    }

	friend std::ostream& operator<<(std::ostream& outStream, const LapackBandMatrixCmplx16&  V)
	{
        long i; long j; std::complex<double> val;

        for(i = 0;  i < V.N; i++)
        {
        for(j = 0; j <  V.N; j++)
        {
          if((j > i + V.ku)||(j < i - V.kl))
          {val = {0.0,0.0};}
          else
          {val = V(i,j);}
          outStream <<   std::scientific << std::setprecision(3) <<  std::right << std::setw(10) << val << " ";
        }
        outStream << std::endl;
        }
        return outStream;
	}


	//  Algebraic operators

    inline void operator=(const LapackBandMatrixCmplx16& B)
	{
    	if(cmplxMdata.isNull())
    	{
    		kl    = B.kl;  ku    = B.ku; N = B.N;
    		cmplxMdata.initialize(B.cmplxMdata);
    	}
    	else
    	{
    	assert(sizeCheck(B.kl, B.ku, B.N));
        cmplxMdata = B.cmplxMdata;
    	}

	}

    inline void operator+=(const  LapackBandMatrixCmplx16& B)
    {
    	assert(sizeCheck(B.kl, B.ku, B.N));
    	cmplxMdata += B.cmplxMdata;
    }

    LapackBandMatrixCmplx16 operator+(const LapackBandMatrixCmplx16& B)
    {
    	assert(sizeCheck(B.kl, B.ku, B.N));

    	LapackBandMatrixCmplx16  C(*this);

    	C.cmplxMdata += B.cmplxMdata;

        return C;
    }

    LapackBandMatrixCmplx16 operator-(const LapackBandMatrixCmplx16& B)
    {
    	assert(sizeCheck(B.kl, B.ku, B.N));

    	LapackBandMatrixCmplx16  C(*this);

    	C.cmplxMdata -= B.cmplxMdata;
    	return C;
    }

    inline void operator-=(const  LapackBandMatrixCmplx16& D)
    {
    	assert(sizeCheck(D.kl, D.ku, D.N));
        cmplxMdata -= D.cmplxMdata;
    }

    inline void operator*=(const double alpha)
    {
    	cmplxMdata *= alpha;
    }

    inline void operator*=(const std::complex<double> alpha)
    {
    	cmplxMdata *= alpha;
    }

    LapackBandMatrixCmplx16 operator*(const double alpha)
    {
    LapackBandMatrixCmplx16 R(*this);
    R *= alpha;
    return R;
    }

    LapackBandMatrixCmplx16 operator*(const std::complex<double> alpha)
    {
    LapackBandMatrixCmplx16 R(*this);
    R *= alpha;
    return R;
    }


    friend LapackBandMatrixCmplx16 operator*(const double alpha, const LapackBandMatrixCmplx16& B)
    {
    LapackBandMatrixCmplx16 R(B);
    R *= alpha;
    return R;
    }

    friend LapackBandMatrixCmplx16 operator*(const std::complex<double> alpha, const LapackBandMatrixCmplx16& B)
    {
    LapackBandMatrixCmplx16 R(B);
    R *= alpha;
    return R;
    }


    inline void operator/=(const double alpha)
    {
    	  cmplxMdata /= alpha;
    }

    inline void operator/=(const std::complex<double> alpha)
    {
           cmplxMdata /= alpha;
    }


    LapackBandMatrixCmplx16 operator/(const double alpha)
    {
    LapackBandMatrixCmplx16 R(*this);
    R /= alpha;
    return R;
    }

    LapackBandMatrixCmplx16 operator/(const std::complex<double> alpha)
    {
    LapackBandMatrixCmplx16 R(*this);
    R /= alpha;
    return R;
    }


    bool isNull() const
    {
    	return cmplxMdata.isNull();
    }


    //
    // The only complex band matrix * matrix operation
    // allowed is a band matrix times a matrix with a
    // single column
    //
    LapackMatrixCmplx16 operator*(const LapackMatrixCmplx16& x)
    {
    assert(sizecheckNx1(x.rows,x.cols));

	LapackMatrixCmplx16 y(x.rows,1);

    char TRANS     = 'N';
    std::complex<double> ALPHA   = {1.0,0.0};
    std::complex<double> BETA    = {0.0,0.0};
    long INCX      = 1;
    long INCY      = 1;
    long LDA       = kl + ku + 1;

    zgbmv_(&TRANS,&N,&N, &kl, &ku,reinterpret_cast<double*>(const_cast< std::complex<double>* >(&ALPHA)),
    this->cmplxMdata.mData.dataPtr, &LDA, x.mData.dataPtr,&INCX,reinterpret_cast<double*>(const_cast< std::complex<double>* >(&BETA)),y.mData.dataPtr,&INCY);
	return y;
}

std::vector< std::complex<double> > operator*(const std::vector< std::complex<double> >& x)
{
	std::vector< std::complex<double> > y(N,0.0);

    char TRANS     = 'N';
    std::complex<double> ALPHA   = {1.0,0.0};
    std::complex<double> BETA    = {0.0,0.0};
    long INCX      = 1;
    long INCY      = 1;
    long LDA       = kl + ku + 1;

    zgbmv_(&TRANS,&N,&N, &kl, &ku,
    reinterpret_cast<double*>(const_cast< std::complex<double>* >(&ALPHA)),this->cmplxMdata.mData.dataPtr, &LDA,
    reinterpret_cast<double*>(const_cast< std::complex<double>* >(&x[0])),&INCX,
    reinterpret_cast<double*>(const_cast< std::complex<double>* >(&BETA)),
    reinterpret_cast<double*>(const_cast< std::complex<double>* >(&y[0])),&INCY);
	return y;
}



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
	    std::cerr  <<  "Offending sizes input  " << "kl : " << dLower << ", ku : " << dUpper  << ", N = " <<  Msize << "\n";
	    return false;
    	}
    	return true;
    }


    bool sizecheckNx1(long rows, long cols) const
    {
    if((rows != N) || (cols != 1))
    {
    std::cerr  <<  "LapackBandMatrixCmplx16 * LapackMatrixCmplx16 error   "  << "\n";
    std::cerr  <<  "LapackMatrixCmplx16 must be N x 1 matrix       "  << "\n";
    std::cerr  <<  "LapackMatrixCmplx16 rows : "  << rows << "\n";
    std::cerr  <<  "LapackMatrixCmplx16 cols : "  << cols << "\n";
    return false;
    }
    return true;
    }

#else
        bool sizeCheck(long dLower, long dUpper, long Msize) const {return true;};
        bool boundsCheck(long, long) const {return true;}
        bool sizecheckNx1(long rows, long cols)  const {return true;};
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    std::cerr << "LapackBandMatrixCmplx16 sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    std::cerr << "LapackBandMatrixCmplx16 sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif



    SCC::LapackMatrixCmplx16 cmplxMdata;

	long ku;
	long kl;
	long  N;
};

} // Namespace SCC

// LAPACK documentation
/*
ubroutine zgbmv	(	character 	trans,
integer 	m,
integer 	n,
integer 	kl,
integer 	ku,
complex*16 	alpha,
complex*16, dimension(lda,*) 	a,
integer 	lda,
complex*16, dimension(*) 	x,
integer 	incx,
complex*16 	beta,
complex*16, dimension(*) 	y,
integer 	incy
)
ZGBMV

Purpose:
 ZGBMV  performs one of the matrix-vector operations

    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or

    y := alpha*A**H*x + beta*y,

 where alpha and beta are scalars, x and y are vectors and A is an
 m by n band matrix, with kl sub-diagonals and ku super-diagonals.
Parameters
[in]	TRANS
          TRANS is CHARACTER*1
           On entry, TRANS specifies the operation to be performed as
           follows:

              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.

              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
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
          ALPHA is COMPLEX*16
           On entry, ALPHA specifies the scalar alpha.
[in]	A
          A is COMPLEX*16 array, dimension ( LDA, N )
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
          X is COMPLEX*16 array, dimension at least
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
          BETA is COMPLEX*16
           On entry, BETA specifies the scalar beta. When BETA is
           supplied as zero then Y need not be set on input.
[in,out]	Y
          Y is COMPLEX*16 array, dimension at least
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

#endif /* SCC_LapackBandMatrixCmplx16 */
