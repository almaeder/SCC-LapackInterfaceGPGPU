/*
 * SCC_LapackHeaders.h
 *
 * LAPACK function prototypes for LAPACK routines used by the collection of classes
 * contained in LapackInterface
 *
 *
 *  Created on: Oct 25, 2017
 *      Author: anderson
 *
 *
 *  Updated : July 27, 2018 (CRA)
 *  Updated : Dec. 09, 2023 (CRA)
 */

/*
#############################################################################
#
# Copyright  2015-2018 Chris Anderson
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

#ifndef SCC_LAPACK_HEADERS_
#define SCC_LAPACK_HEADERS_

#pragma once
extern "C" {
    #ifdef USE_MKL
    #include <mkl_blas.h>
    #include <mkl_cblas.h>
    #include <mkl_lapacke.h>
    #include <mkl_lapack.h>
    #else
    #include <blas.h>
    #include <cblas.h>
    #include <lapacke.h>
    #include <lapack.h>
    #endif
}
#define RC_INT int

#endif /* SCC_LAPACKHEADERS_H_ */



