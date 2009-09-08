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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <ietl/traits.h>

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
    boost::numeric::ublas::axpy_prod(m, x, y, true);
}

} // end namespace ietl

#endif // ALPS_DIAG_PARAPACK_IETL_INTERFACE_H_
