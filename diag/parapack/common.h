/*****************************************************************************
*
* ALPS/diag/parapack: Full and Sparse diagonalization
*                     for quantum lattice models using parapack scheduler
*
* Copyright (C) 2003-2009 by Synge Todo <wistaria@comp-phys.org>,
*                            Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>
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

#ifndef ALPS_DIAG_PARAPACK_COMMON_H_
#define ALPS_DIAG_PARAPACK_COMMON_H_

#include <alps/parameter.h>
#include <alps/alea.h>
#include <alps/model.h>
#include <alps/lattice.h>

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>

#include <cmath>
#include <cstddef>

#ifdef _OPENMP
# include <omp.h>
#endif

namespace {

template <typename M, typename I, typename G>
void
add_to_matrix(
    M& matrix,
    alps::HamiltonianDescriptor<I> const& hd,
    alps::basis_states<I> const& basis_set,
    typename alps::graph_traits<G>::vertex_descriptor const& vd,
    G const& graph,
    alps::Parameters const& params)
{
    typedef typename M::value_type value_type;
    typedef alps::basis_states<I> basis_set_type;

    alps::BasisDescriptor<I> const& basis(hd.basis());

    int const t = get(alps::site_type_t(),  graph, vd);
    int const s = get(alps::site_index_t(), graph, vd);
    std::size_t const dim = basis_set.size();
    std::size_t const ds  = basis_set.basis().get_site_basis(s).num_states();

    boost::multi_array<value_type, 2> const
        site_matrix(
            alps::get_matrix(value_type(),
                             hd.site_term(t),
                             basis.site_basis(t),
                             params));

#pragma omp barrier
#pragma omp parallel for
    for (std::size_t i = 0; i < dim; ++i) {
        std::size_t const is = basis_set[i][s];
        typename basis_set_type::value_type target(basis_set[i]);
        for (std::size_t js = 0; js < ds; ++js) {
            target[s] = js;
            std::size_t const j = basis_set.index(target);
            if (j < dim) {
                matrix(i, j) += site_matrix[is][js];
            }
        }
    }
#pragma omp barrier
}

template <typename M, typename I, typename G>
void
add_to_matrix(
    M& matrix,
    alps::HamiltonianDescriptor<I> const& hd,
    alps::basis_states<I> const& basis_set,
    typename alps::graph_traits<G>::bond_descriptor const& ed,
    G const& graph,
    alps::Parameters const& params)
{
    typedef typename M::value_type value_type;
    typedef alps::basis_states<I> basis_set_type;

    alps::BasisDescriptor<I> const& basis(hd.basis());

    typename alps::graph_traits<G>::site_descriptor const& vd0 =
        alps::detail::source_wrap<G>(ed, graph);
    typename alps::graph_traits<G>::site_descriptor const& vd1 =
        alps::detail::target_wrap<G>(ed, graph);

    int const t   = get(alps::bond_type_t(), graph, ed);
    int const st0 = get(alps::site_type_t(), graph, vd0);
    int const st1 = get(alps::site_type_t(), graph, vd1);
    int const s0  = get(alps::site_index_t(), graph, vd0);
    int const s1  = get(alps::site_index_t(), graph, vd1);
    std::size_t const dim = basis_set.size();
    std::size_t const ds0 = basis_set.basis().get_site_basis(s0).num_states();
    std::size_t const ds1 = basis_set.basis().get_site_basis(s1).num_states();

    boost::multi_array<value_type, 4> const
        bond_matrix(
            alps::get_matrix(value_type(),
                             hd.bond_term(t),
                             basis.site_basis(st0),
                             basis.site_basis(st1),
                             params));

#pragma omp barrier
#pragma omp parallel for
    for (std::size_t i = 0; i < dim; ++i) {
        std::size_t const is0 = basis_set[i][s0];
        std::size_t const is1 = basis_set[i][s1];
        typename basis_set_type::value_type target(basis_set[i]);
        for (std::size_t js0 = 0; js0 < ds0; ++js0) {
            for (std::size_t js1 = 0; js1 < ds1; ++js1) {
                target[s0] = js0;
                target[s1] = js1;
                std::size_t const j = basis_set.index(target);
                if (j < dim) {
                    matrix(i, j) += bond_matrix[is0][is1][js0][js1];
                }
            }
        }
    }
#pragma omp barrier
}

template <class V>
std::pair<double, double>
static_average2(double beta, V const& evals)
{
    typedef V vector_type;
    double val  = 0.0;
    double val2 = 0.0;
    double offset = evals(0);
    BOOST_REVERSE_FOREACH(typename vector_type::value_type const& eval, evals) {
        double weight = std::exp(-beta * (eval - offset));
        val  += eval * weight;
        val2 += alps::abs2(eval) * weight;
    }
    return std::make_pair(val, val2);
}

} // end namespace

#endif // ALPS_DIAG_PARAPACK_COMMON_H_
