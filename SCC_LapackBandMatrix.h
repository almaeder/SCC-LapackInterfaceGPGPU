/*
 * SCC_LapackBandMatrix.h
 *
 *  Created on: Aug 17, 2018
 *      Author: anderson
 */
//
// This class wraps a SCC::LapackMatrix instance that contains
// an N x N banded matrix stored using the Lapack band storage
// convention.
//
// kl = lower bandwidth
// ku = upper bandwidth
//  N = system size
//
// The indexing for the access operator(i,j) starts at 0, so that it
// is consistent with the indexing of the LapackMatrix class.
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
//
//
/*
#############################################################################
#
# Copyright 2018 Chris Anderson
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

#include <sstream>
#include <stdexcept>

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
		Sp.initialize();
		kl = 0;
		ku = 0;
    	N  = 0;
	}

	void initialize(const LapackBandMatrix& S)
	{
		kl = S.kl;
		ku = S.ku;
		N  = S.N;
	    Sp.initialize(S.Sp);

	}
    void initialize(long kl, long ku, long N)
	{
    	this->kl = kl;
    	this->ku = ku;
    	this->N  = N;
    	Sp.initialize(kl + ku + 1, N);
	}

   void operator=(const LapackBandMatrix& S)
   {
	   if(this->N == 0){initialize(S);}
	   else
	   {
	   assert(sizeCheck(S.kl, S.ku, S.N));
	   Sp.initialize(S.Sp);
	   }
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

//
// Algebraic operators utilize algebraic operations of underlying LapackMatric
//




/*!  Outputs the matrix values to the screen with the (0,0) element in the upper left corner  */

friend std::ostream& operator<<(std::ostream& outStream, const LapackMatrix&  V)
{
        long i; long j; double val;

        for(i = 0;  i < V.N; i++)
        {
        for(j = 0; j <  V.N; j++)
        {
          if((j > i + ku)||(j < i - kl))
          {val = 0.0;}
          else
          {val = V(i,j); }
          outStream <<   std::scientific << std::setprecision(3) <<  std::right << std::setw(10) << val << " ";
        }
        outStream << std::endl;
        }
        return outStream;
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
	    std::cerr <<  "Offending sizes input  " << "kl : " << dLower << ", ku : " << dUpper  << ", N = " <<  Msize << "\n";
	    return false;
    	}
    	return true;
    }
#else
        bool boundsCheck(long, long) const {return true;}
        bool sizeCheck(long dLower, long dUpper, long Msize) const {return true;}
#endif



    SCC::LapackMatrix Sp;

	long ku;
	long kl;
	long  N;
};

} // Namespace SCC


#endif /* SCC_LapackBandMatrix_  */
