/*
 * SCC_LapackMatrixCmplx16.h
 *
 *  Modified on: June 1, 2023
 *      Author: anderson
 */
//
// A matrix class to facilitate the use of LAPACK routines for COMPLEX*16 Fortran data
// types.
//
//                       ----->  BETA version <-----
//
// The data for the matrix is stored by columns (Fortran convention) with each complex matrix
// element value stored as alternating doubles containing the real and imaginary part of that
// value.
//
// It is assumed that the data storate for std::complex<double> contains the
// real and imaginary components in adjacent memory location with the real component
// first. In addition it is assumed that the data for a std::vector<std::complex<double>
// is stored in contiguous memory locations with double values for the real and imaginary
// and compoments alternating, e.g. the storage convention used by FORTRAN for complex*16.
//
// Internally the data storage uses an SCC::LapackMatrix to faciliate the implementation
// of algebraic operations.
/*
#############################################################################
#
# Copyright 2021- Chris Anderson
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
#include "LapackInterface/SCC_LapackMatrix.h"
#include <complex>
#include <cassert>

#ifndef LAPACK_MATRIX_CMPLX_16_H_
#define LAPACK_MATRIX_CMPLX_16_H_

//
// Prototypes for the only two LAPACK BLAS routines used by this class
//
extern "C" void zgemm_(char* TRANSA,char* TRANSB,long* M, long*N ,long* K,double* ALPHA,
                       double* A,long* LDA,double* B, long* LDB,double* BETA,double* C,long* LDC);

extern "C" void zgemv_(char* TRANS, long* M, long* N, double* alpha, double* Aptr,
                       long* LDA, double* Xptr, long* INCX, double* BETA, double* Yptr, long* INCY);

namespace SCC
{
class LapackMatrixCmplx16
{
public:

    LapackMatrixCmplx16()
	{
    	this->rows = 0;
    	this->cols = 0;
	}

    LapackMatrixCmplx16(const LapackMatrixCmplx16& C)
    {
    	this->rows  = C.rows;
    	this->cols  = C.cols;
    	this->mData = C.mData;
    }

	LapackMatrixCmplx16(long M, long N)
	{
		initialize(M,N);
	}

	LapackMatrixCmplx16(const SCC::LapackMatrix& realA, const SCC::LapackMatrix& imagA)
	{
		initialize(realA,imagA);
	}

    void initialize()
	{
	    this->rows = 0;
		this->cols = 0;
		mData.initialize();
	}

	void initialize(long M, long N)
	{
	    this->rows = M;
		this->cols = N;
		mData.initialize(2*rows,cols);
		mData.setToValue(0.0);
	}

	void initialize(const LapackMatrixCmplx16& C)
    {
    	this->rows  = C.rows;
    	this->cols  = C.cols;
    	this->mData.initialize(C.mData);
    }

	void initialize(const SCC::LapackMatrix& realA, const SCC::LapackMatrix& imagA)
	{
	    this->rows = realA.getRowDimension();
		this->cols = realA.getColDimension();
		if((realA.getRowDimension() != imagA.getRowDimension())
		 ||(realA.getColDimension() != imagA.getColDimension()))
		 {
			throw std::runtime_error("\nIncorrect dimension input matrices in \nLapackMatrixCmplx (realA,imagA) constructor.\n");
		 }

		mData.initialize(2*rows,cols);

		for(long j = 0; j < cols; j++)
		{
		for(long i = 0; i < rows; i++)
		{
		mData(2*i,  j) = realA(i,j);
		mData(2*i+1,j) = imagA(i,j);
		}}
	}

	long getRowDimension() const {return rows;}
	long getColDimension() const {return cols;}

	inline void insert(long i, long j, double vReal, double vCplx)
	{
		 mData(2*i,j)      = vReal;
		 mData(2*i + 1,j)  = vCplx;
	}

	inline void extract(long i, long j, double& vReal, double& vCplx) const
	{
	     vReal = mData(2*i,j);
		 vCplx = mData(2*i + 1,j);
	}

    inline void insert(long i, long j, std::complex<double> z)
	{
		 mData(2*i,j)      = z.real();
		 mData(2*i + 1,j)  = z.imag();
	}

	inline void extract(long i, long j, std::complex<double>& z) const
	{
	     z = std::complex<double>(mData(2*i,j),mData(2*i + 1,j));
	}

    /*!
    Returns a reference to the element with index (i,j) - indexing
    starting at (0,0). Using the fact that the pointer to a complex<double> value
    is a pointer to the first of two consecutive doubles storing the
    complex value.
    */

	#ifdef _DEBUG
    std::complex<double>&  operator()(long i, long j)
    {
    assert(boundsCheck(i, 0, rows-1,1));
    assert(boundsCheck(j, 0, cols-1,2));
    return *(reinterpret_cast<std::complex<double>*>((mData.dataPtr +  (2*i) + j*(2*rows))));
    };

    const std::complex<double>&  operator()(long i, long j) const
    {
    assert(boundsCheck(i, 0, rows-1,1));
    assert(boundsCheck(j, 0, cols-1,2));
    return *(reinterpret_cast<std::complex<double>*>((mData.dataPtr +  (2*i) + j*(2*rows))));
    };
#else
    /*!
    Returns a reference to the element with index (i,j) - indexing
    starting at (0,0). Using the fact that the pointer to a complex<double> value
    is a pointer to the first of two consecutive doubles storing the
    complex value.
    */
    inline std::complex<double>&  operator()(long i, long j)
    {
    	return *(reinterpret_cast<std::complex<double>*>((mData.dataPtr +  (2*i) + j*(2*rows))));
    };

    inline const std::complex<double>&  operator()(long i, long j) const
    {
    return *(reinterpret_cast<std::complex<double>*>((mData.dataPtr +  (2*i) + j*(2*rows))));;
    };
#endif


    double normFrobenius() const
    {
	double valSum = 0.0;

	for(long j = 0; j < cols; j++)
	{
	for(long i = 0; i < rows; i++)
	{
    		valSum += std::norm(this->operator()(i,j));
    }}
    return std::sqrt(valSum);
    }

    void getColumn(long colIndex, std::vector< std::complex<double>> & Mcol)
    {
    	Mcol.resize(rows);
    	for(long i = 0; i < rows; i++)
    	{
		Mcol[i] = this->operator()(i,colIndex);
    	}
    }

    void getColumn(long colIndex, LapackMatrixCmplx16 & Mcol)
    {
    	Mcol.initialize(rows,1);
    	for(long i = 0; i < rows; i++)
    	{
		Mcol(i,0) = this->operator()(i,colIndex);
    	}
    }

	void getRealAndCmplxMatrix(LapackMatrix& realA, LapackMatrix& imagA) const
	{
		realA.initialize(rows,cols);
		imagA.initialize(rows,cols);

	    for(long j = 0; j < cols; j++)
		{
		for(long i = 0; i < rows; i++)
		{
		realA(i,j) = mData(2*i,j);
		imagA(i,j) = mData(2*i+1,j);
		}}
	}

	void getRealAndCmplxColumn(long colIndex, std::vector<double>& realCol, std::vector<double>& imagCol)
	{
		assert(boundsCheck(colIndex, 0, cols-1,2));
		realCol.resize(rows);
		imagCol.resize(rows);

	    for(long i = 0; i < rows; i++)
		{
		realCol[i] = mData(2*i,colIndex);
		imagCol[i] = mData(2*i+1,colIndex);
		}
	}

	LapackMatrixCmplx16 createUpperTriPacked()
	{
		if(rows != cols)
		{
			throw std::runtime_error("\nLapackMatrixCmplx16: No conversion of non-square matrix \nto upper traingular packed form.\n");
		}

		LapackMatrixCmplx16 AP((rows*(rows+1))/2,1);

		long     ind; long     jnd;
		double vReal; double vImag;

		for(long j = 1; j <=cols; j++)
		{
		for(long i = 1; i <= j;   i++)
		{
            ind = i-1;
            jnd = j-1;
            extract(ind,jnd,vReal,vImag);

            ind = (i + (j-1)*j/2) - 1;
            AP.insert(ind,  0,vReal,vImag);
		}}


		return AP;
	}


//  C := alpha*op( A )*op( B ) + beta*C,

LapackMatrixCmplx16 operator*(const LapackMatrixCmplx16& B) const
{
    assert(sizeCheck(this->cols,B.rows));

    LapackMatrixCmplx16 C(this->rows,B.cols);

    char TRANSA = 'N';
    char TRANSB = 'N';

    long M       = this->rows;
    long N       = B.cols;
    long K       = this->cols;

    std::complex<double> ALPHA = {1.0,0.0};
    std::complex<double> BETA  = {0.0,0.0};

    double*Aptr  = mData.getDataPointer();
    double*Bptr  = B.mData.getDataPointer();
    double*Cptr  = C.mData.getDataPointer();
    long LDA     = this->rows;
    long LDB     = B.rows;
    long LDC     = C.rows;

    zgemm_(&TRANSA,&TRANSB,&M,&N,&K,reinterpret_cast<double*>(&ALPHA), Aptr,&LDA,Bptr,&LDB,reinterpret_cast<double*>(&BETA),Cptr,&LDC);
    return C;
}


std::vector< std::complex<double> > operator*(const std::vector< std::complex<double> >& x)
{
	std::vector< std::complex<double> > y(rows,0.0);

    char TRANS     = 'N';
    std::complex<double> ALPHA = {1.0,0.0};
    std::complex<double> BETA  = {0.0,0.0};
    long INCX      = 1;
    long INCY      = 1;

    zgemv_(&TRANS,&rows,&cols,reinterpret_cast<double*>(&ALPHA),mData.getDataPointer(),&rows,
    reinterpret_cast<double*>(const_cast< std::complex<double>* >(&x[0])),&INCX,reinterpret_cast<double*>(&BETA),reinterpret_cast<double*>(&y[0]),&INCY);
	return y;
}



#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end,int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        std::cerr << "LapackMatrix index " << coordinate << " out of bounds " << std::endl;
        std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "] " << std::endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long,int) const {return true;}
#endif


#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    std::cerr << "LapackMatrixCmplx16 sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    std::cerr << "LapackMatrixCmplx16 sizes are incompatible : " << size1 << " != " << size2 << " ";
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif


	long rows;
	long cols;
	SCC::LapackMatrix mData;

};
};

#endif /* LAPACK_MATRIX_CMPLX_16_H__ */
