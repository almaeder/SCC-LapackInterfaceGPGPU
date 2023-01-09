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
