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

#include "common.h"
#include "sparsediag.h"

#include <alps/parameter.h>
#include <alps/alea.h>

#include <ietl/interface/ublas.h>
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>

#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/random.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <complex>
#include <cmath>
#include <cstddef>

namespace alps{
namespace diag{

sparsediag_worker::sparsediag_worker(alps::Parameters const& ps)
    : alps::parapack::abstract_worker(),
      alps::graph_helper<>(ps),
      alps::model_helper<>(*this, ps),
      params(ps),
      done(false)
{}

void
sparsediag_worker::init_observables(alps::Parameters const& /* params */,
                                    alps::ObservableSet& /* obs */)
{}

inline
bool
sparsediag_worker::is_thermalized() const
{ return true; }

inline
double
sparsediag_worker::progress() const
{ return done ? 1 : 0; }

void
sparsediag_worker::run(alps::ObservableSet& obs)
{
    typedef typename boost::numeric::ublas::mapped_vector_of_mapped_vector<double, boost::numeric::ublas::row_major> matrix_type;
    typedef typename boost::numeric::ublas::vector<double> vector_type;

    if (done) { return; }
    done = true;

    // measurements
    std::map<std::string, double> m;
    m["Number of Sites"] = num_sites();
    m["Volume"] = volume();

    // generate basis set
    alps::basis_states<short>
        basis_set(alps::basis_states_descriptor<short>(model().basis(), graph()));
    std::size_t dim = basis_set.size();

    // generate Hamiltonian matrix
    matrix_type Hamiltonian(dim, dim);
    // Hamiltonian.clear();
    BOOST_FOREACH(site_descriptor s, sites()) {
        add_to_matrix(Hamiltonian, model(), basis_set, s, graph(), params);
    }
    BOOST_FOREACH(bond_descriptor b, bonds()) {
        add_to_matrix(Hamiltonian, model(), basis_set, b, graph(), params);
    }

    // Lanczos diagonalization of Hamiltonian matrix
    vector_type evals;
    {
        typedef ietl::vectorspace<vector_type> vectorspace_type;
        boost::mt19937 generator;
        vectorspace_type vec(dim);
        ietl::lanczos<matrix_type, vectorspace_type> lanczos(Hamiltonian, vec);

        int max_iter = std::min(static_cast<int>(10 * dim), 1000);
        if (params.defined("MAX_ITERATIONS")) {
            max_iter = params["MAX_ITERATIONS"];
        }
        std::size_t num_eigenvalues = 1;
        if (params.defined("NUMBER_EIGENVALUES")) {
            num_eigenvalues = params["NUMBER_EIGENVALUES"];
        }
        ietl::lanczos_iteration_nlowest<double> iter(max_iter, num_eigenvalues);
        std::cerr << "Starting Lanczos... " << std::flush;
        lanczos.calculate_eigenvalues(iter, generator);
        std::cerr << "Done.\n";

        int n = std::min(num_eigenvalues, lanczos.eigenvalues().size());
        evals.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            evals[i] = lanczos.eigenvalues()[i];
        }
    }

    double E0 = evals[0];

    m["Ground State Energy"] = E0;
    m["Ground State Energy Density" ] = E0 / volume();

    // Store measurements
    typedef std::pair<std::string, double> pair_type;
    BOOST_FOREACH(pair_type const& p, m) {
        obs << alps::SimpleRealObservable(p.first);
    }
    obs.reset(true);
    BOOST_FOREACH(pair_type const& p, m) {
        obs[p.first] << p.second;
    }
}

void
sparsediag_worker::save(alps::ODump& dump) const
{ dump << done; }

void
sparsediag_worker::load(alps::IDump& dump)
{ dump >> done; }

} // end namespace diag
} // end namespace alps

namespace {

PARAPACK_REGISTER_WORKER(alps::diag::sparsediag_worker, "Sparse diagonalization");

} // end namespace
