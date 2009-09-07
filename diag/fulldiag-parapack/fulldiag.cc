/*****************************************************************************
*
* ALPS/fulldiag-parapack: Full diagonalization for quantum lattice systems
*                         using parapack scheduler
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

#include "fulldiag.h"

#include <alps/parameter.h>
#include <alps/alea.h>

#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/matrix_traits.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include <cmath>
#include <cstddef>

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

    int t = get(alps::site_type_t(),  graph, vd);
    int s = get(alps::site_index_t(), graph, vd);
    std::size_t dim = basis_set.size();
    std::size_t ds  = basis_set.basis().get_site_basis(s).num_states();

    boost::multi_array<value_type, 2>
        site_matrix(
            alps::get_matrix(value_type(),
                             hd.site_term(t),
                             basis.site_basis(t),
                             params));

    for (std::size_t i = 0; i < dim; ++i) {
        std::size_t is = basis_set[i][s];
        for (std::size_t js = 0; js < ds; ++js) {
            typename basis_set_type::value_type target(basis_set[i]);
            target[s] = js;
            std::size_t j = basis_set.index(target);
            if (j < dim) {
                matrix(i, j) += site_matrix[is][js];
            }
        }
    }
}

template <typename M, typename I, typename G>
void
add_to_matrix(
    M& matrix,
    alps::HamiltonianDescriptor<I> const& hd,
    alps::basis_states<I> const& basis_set,
    typename alps::graph_traits<G>::bond_descriptor const& ed,
    typename alps::graph_traits<G>::site_descriptor const& vd0,
    typename alps::graph_traits<G>::site_descriptor const& vd1,
    G const& graph,
    alps::Parameters const& params)
{
    typedef typename M::value_type value_type;
    typedef alps::basis_states<I> basis_set_type;

    alps::BasisDescriptor<I> const& basis(hd.basis());

    int t   = get(alps::bond_type_t(), graph, ed);
    int st0 = get(alps::site_type_t(), graph, vd0);
    int st1 = get(alps::site_type_t(), graph, vd1);
    int s0  = get(alps::site_index_t(), graph, vd0);
    int s1  = get(alps::site_index_t(), graph, vd1);
    std::size_t dim = basis_set.size();
    std::size_t ds0 = basis_set.basis().get_site_basis(s0).num_states();
    std::size_t ds1 = basis_set.basis().get_site_basis(s1).num_states();

    boost::multi_array<value_type, 4>
        bond_matrix(
            alps::get_matrix(value_type(),
                             hd.bond_term(t),
                             basis.site_basis(st0),
                             basis.site_basis(st1),
                             params));

    for (std::size_t i = 0; i < dim; ++i) {
        std::size_t is0 = basis_set[i][s0];
        std::size_t is1 = basis_set[i][s1];
        for (std::size_t js0 = 0; js0 < ds0; ++js0) {
            for (std::size_t js1 = 0; js1 < ds1; ++js1) {
                typename basis_set_type::value_type target(basis_set[i]);
                target[s0] = js0;
                target[s1] = js1;
                std::size_t j = basis_set.index(target);
                if (j < dim) {
                    matrix(i, j) += bond_matrix[is0][is1][js0][js1];
                }
            }
        }
    }
}

template <typename T, typename R, typename A, typename V>
void
diagonalize(
    boost::numeric::ublas::matrix<T, R, A>& a,
    V& v, bool need_eigenvectors = true)
{
    using namespace boost::numeric::bindings::lapack;

    BOOST_STATIC_ASSERT((boost::is_same<typename R::orientation_category, boost::numeric::ublas::column_major_tag>::value));

    const char jobz = (need_eigenvectors ? 'V' : 'N');
    const char uplo = 'L';

    int info = syev(jobz, uplo, a, v, optimal_workspace());

    if (info != 0) {
        throw std::runtime_error("failed in syev");
    }
}

} // end namespace

namespace alps{
namespace diag{

fulldiag_worker::fulldiag_worker(alps::Parameters const& ps)
    : alps::parapack::abstract_worker(),
      alps::graph_helper<>(ps),
      alps::model_helper<>(*this, ps),
      params(ps),
      done(false)
{}

void
fulldiag_worker::init_observables(alps::Parameters const& /* params */,
                                  alps::ObservableSet& /* obs */)
{}

inline
bool
fulldiag_worker::is_thermalized() const
{ return true; }

inline
double
fulldiag_worker::progress() const
{ return done ? 1 : 0; }

void
fulldiag_worker::run(alps::ObservableSet& obs)
{
    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
    typedef boost::numeric::ublas::vector<double> diagonal_matrix_type;

    typedef boost::numeric::ublas::vector<double> vector_type;

    if (done) { return; }
    done = true;

    double beta = 1.0;
    if (params.defined("T")) {
        beta = 1.0 / alps::evaluate<double>("T", params);
    }

    // measurements
    std::map<std::string, double> m;
    m["Temperature"] = 1.0 / beta;
    m["Inverse Temperature"] = beta;
    m["Number of Sites"] = num_sites();
    m["Volume"] = volume();

    // generate basis set
    alps::basis_states<short>
        basis_set(alps::basis_states_descriptor<short>(model().basis(), graph()));
    std::size_t dim = basis_set.size();

    // generate Hamiltonian matrix
    matrix_type Hamiltonian(dim, dim);
    Hamiltonian.clear();
    BOOST_FOREACH(site_descriptor s, sites()) {
        add_to_matrix(Hamiltonian, model(), basis_set, s, graph(), params);
    }
    BOOST_FOREACH(bond_descriptor b, bonds()) {
        add_to_matrix(Hamiltonian, model(), basis_set, b, source(b, graph()), target(b, graph()), graph(), params);
    }
    diagonal_matrix_type diagonal_energy(dim);
    for (std::size_t i = 0; i < dim; ++i) {
        diagonal_energy(i) = Hamiltonian(i, i);
    }

    // partition function coefficients
    matrix_type matrix(dim, dim);
    matrix.clear();
    for (std::size_t i = 0; i < dim; ++i) {
        matrix(i, i) = 1;
    }
    double factor = 1.;
    for (std::size_t k = 0; k < 5; ++k) {
        if (k != 0) {
            matrix = prod(Hamiltonian, matrix);
            factor *= -k;
        }
        double trace = 0.;
        for (std::size_t i = 0; i < dim; ++i) {
            trace += matrix(i, i);
        }
        m["Partition Function Coefficient #" + boost::lexical_cast<std::string>(k)]
            = trace / factor;
    }

    // diagonalize Hamiltonian matrix
    vector_type evals(dim);
    std::cerr << "start diagonalization... " << std::flush;
    diagonalize(Hamiltonian, evals);
    std::cerr << "done\n";

    double E0 = evals(0);
    double Z  = 0.;
    BOOST_REVERSE_FOREACH(double eval, evals) {
        double weight = std::exp(-beta * (eval - E0));
        Z += weight;
    }
}

void
fulldiag_worker::save(alps::ODump& dump) const
{ dump << done; }

void
fulldiag_worker::load(alps::IDump& dump)
{ dump >> done; }

} // end namespace diag
} // end namespace alps

namespace {

PARAPACK_REGISTER_WORKER(alps::diag::fulldiag_worker, "Full diagonalization");

} // end namespace
