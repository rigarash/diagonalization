/*****************************************************************************
*
* ALPS/diagonalization: Full and Sparse diagonalization
*                       for quantum lattice models using parapack scheduler
*
* Copyright (C) 2009-2012 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_DIAG_PARAPACK_IETL_INTERFACE_H_
#define ALPS_DIAG_PARAPACK_IETL_INTERFACE_H_

#include <alps/config.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <ietl/traits.h>

#ifdef ALPS_HAVE_MKL
# ifdef _OPENMP
#  include <omp.h>
# endif
# include <mkl_spblas.h>
# include <vector>
#endif

namespace ietl {

template <typename T, class Gen>
inline
void
generate(boost::numeric::ublas::vector<T>& c, Gen& gen) {
    std::generate(c.begin(),c.end(),gen);
}

template <typename T>
inline
void
clear(boost::numeric::ublas::vector<T>& c) {
    c.clear();
}

template <typename T, typename S>
inline
T
dot(boost::numeric::ublas::vector<T, S> const& x,
    boost::numeric::ublas::vector<T, S> const& y) {
    return boost::numeric::ublas::inner_prod(boost::numeric::ublas::conj(x), y);
}

template <typename T>
inline
typename number_traits<T>::magnitude_type
two_norm(boost::numeric::ublas::vector<T>& x) {
    return boost::numeric::ublas::norm_2(x);
}

template <typename T>
void
copy(boost::numeric::ublas::vector<T> const& x,
     boost::numeric::ublas::vector<T>& y) {
    y.assign(x);
}

template <typename T>
inline
void
mult(boost::numeric::ublas::compressed_matrix<T, boost::numeric::ublas::row_major> const& m,
     boost::numeric::ublas::vector<T> const& x,
     boost::numeric::ublas::vector<T>& y) {
#ifdef ALPS_HAVE_MKL
    char uplo = 'N';
    int rowsize = m.size1();
    std::vector<double> value(m.value_data().size());
    std::vector<int> index1(m.index1_data().size());
    std::vector<int> index2(m.index2_data().size());
    std::vector<double> xtmp(x.size());
    std::vector<double> ytmp(m.size1());
    y.resize(ytmp.size());
    y.clear();
#pragma omp barrier
#pragma omp parallel
    {
#pragma omp for nowait
        for (std::size_t i = 0; i < value.size(); ++i) {
            value[i] = m.value_data()[i];
        }
#pragma omp for nowait
        for (std::size_t i = 0; i < index1.size(); ++i) {
            index1[i] = m.index1_data()[i];
        }
#pragma omp for nowait
        for (std::size_t i = 0; i < index2.size(); ++i) {
            index2[i] = m.index2_data()[i];
        }
#pragma omp for nowait
        for (std::size_t i = 0; i < xtmp.size(); ++i) {
            xtmp[i] = x(i);
        }
    } // omp parallel
#pragma omp barrier
    mkl_set_dynamic();
    mkl_cspblas_dcsrgemv(&uplo, &rowsize, &value[0], &index1[0], &index2[0], &xtmp[0], &ytmp[0]);
#pragma omp barrier
#pragma omp parallel for
    for (std::size_t i = 0; i < y.size(); ++i) {
        y(i) = ytmp[i];
    }
#pragma omp barrier
#else
    boost::numeric::ublas::axpy_prod(m, x, y, true);
#endif
}

template <typename M, typename T>
inline
void
mult(M const& m,
     boost::numeric::ublas::vector<T> const& x,
     boost::numeric::ublas::vector<T>& y) {
    boost::numeric::ublas::axpy_prod(m, x, y, true);
}

} // end namespace ietl

#endif // ALPS_DIAG_PARAPACK_IETL_INTERFACE_H_
