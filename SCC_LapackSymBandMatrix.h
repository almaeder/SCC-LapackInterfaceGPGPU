/*
 * SCC_LapackSymBandMatrix.h
 *
 *  Created on: Nov 29, 2022
 *      Author: anderson
 */

// A class for representing a dense matrix in Lapack symmetric band storage.
//
// This routine packs the upper triangular portion of the symmetric band matrix
// data so that only indices j >= i with j < i + ku are allowed where ku is
// the number of super diagonals of the matrix.
//
// Note: SCC_LapackSymBandMatrix stores the upper diagonal
// data of a symmetric matrix so the UPLO parameter in the
// Lapack symmetric routines should be specified as "U".
//
// Lapack dependencies : dsbmv_
//
/*
#############################################################################
#
# Copyright  2015-2020 Chris Anderson
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

#ifndef LAPACK_SYM_BAND_MATRIX_
#define LAPACK_SYM_BAND_MATRIX_

namespace SCC
{
class LapackSymBandMatrix
{
	public:

	LapackSymBandMatrix()
	{
	initialize();
	}

    LapackSymBandMatrix(const LapackSymBandMatrix& S)
	{
    	initialize(S);
	}

	LapackSymBandMatrix(long ku, long N)
	{
	    initialize(ku,N);
	}

	void initialize()
	{
		Sp.initialize();
		ku = 0;
    	N  = 0;
	}

	void initialize(const LapackSymBandMatrix& S)
	{
		ku = S.ku;
		N  = S.N;
	    Sp.initialize(S.Sp);

	}
    void initialize(long ku, long N)
	{
    	this->ku = ku;
    	this->N  = N;
    	Sp.initialize(ku + 1, N);
	}

	#ifdef _DEBUG
	double& operator()(long i, long j)
	{
		assert(boundsCheck(i,j));
		return Sp(ku +  (i-j),j);
	}

    const double& operator()(long i, long j) const
	{
    	assert(boundsCheck(i,j));
		return Sp(ku +  (i-j),j);
	}


	#else
	inline double& operator()(long i, long j)
	{

		return Sp(ku +  (i-j),j);
	}

    inline const double& operator()(long i, long j) const
	{
		return Sp(ku +  (i-j),j);
	}
	#endif

	double* getDataPointer() const {return Sp.dataPtr;}

    void setToValue(double val)
	{
		Sp.setToValue(val);
	}

    long getDimension()
    {
    	return N;
    }

    long getSuperDiagonalCount()
    {
    	return ku;
    }


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

// Fortran indexing bounds check

#ifdef _DEBUG
	bool boundsCheck(long i, long j) const
	{
        if((i > j)||(j > i+ku))
	    {
	    std::cerr  <<  "Symmetric and matrix storage error " << std::endl;
	    std::cerr  <<  "ku =  " << ku  <<  std::endl;
	    std::cerr  <<  "N  =  " << N   <<  std::endl;
	    std::cerr  <<  "Offending indices " << "(" << i << ", " << j << ")" << std::endl;
	    return false;
	    }
	    return true;
	}
#else
        bool boundsCheck(long, long) const {return true;}
#endif

    SCC::LapackMatrix Sp;

	long ku;
	long  N;
};

} // Namespace SCC

#endif /* LAPACK_SYM_BAND_MATRIX_ */
