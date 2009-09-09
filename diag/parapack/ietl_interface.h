/*****************************************************************************
*
* ALPS/diag/parapack: Full and Sparse diagonalization
*                     for quantum lattice models using parapack scheduler
*
* Copyright (C) 2009-2009 by Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>
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

template <typename M, typename T>
inline
void
mult(M const& m,
     boost::numeric::ublas::vector<T> const& x,
     boost::numeric::ublas::vector<T>& y) {
#ifdef ALPS_HAVE_MKL
    typename boost::numeric::ublas::compressed_matrix<typename M::value_type, boost::numeric::ublas::row_major> mtmp = m;
    char uplo = 'N';
    int rowsize = mtmp.size1();
    const int v_size = mtmp.value_data().size();
    std::vector<double> value(v_size);
    for (std::size_t i = 0; i < v_size; ++i) {
        value[i] = mtmp.value_data()[i];
    }
    const int i1_size = mtmp.index1_data().size();
    std::vector<int> index1(i1_size);
    for (std::size_t i = 0; i < i1_size; ++i) {
        index1[i] = mtmp.index1_data()[i];
    }
    const int i2_size = mtmp.index2_data().size();
    std::vector<int> index2(i2_size);
    for (std::size_t i = 0; i < i2_size; ++i) {
        index2[i] = mtmp.index2_data()[i];
    }
    const int x_size = x.size();
    std::vector<double> xtmp(x_size);
    for (std::size_t i = 0; i < x_size; ++i) {
        xtmp[i] = x(i);
    }
    const int y_size = y.size();
    std::vector<double> ytmp(y_size);
    for (std::size_t i = 0; i < y_size; ++i) {
        ytmp[i] = y(i);
    }
    mkl_cspblas_dcsrgemv(&uplo, &rowsize, &value[0], &index1[0], &index2[0], &xtmp[0], &ytmp[0]);
    for (std::size_t i = 0; i < y_size; ++i) {
        y(i) = ytmp[i];
    }
#else
    boost::numeric::ublas::axpy_prod(m, x, y, true);
#endif
}

} // end namespace ietl

#endif // ALPS_DIAG_PARAPACK_IETL_INTERFACE_H_
