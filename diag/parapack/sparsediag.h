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

#ifndef ALPS_DIAG_PARAPACK_SPARSEDIAG_H_
#define ALPS_DIAG_PARAPACK_SPARSEDIAG_H_

#include "matrix_worker.h"
#include "common.h"
#include <alps/parameter.h>
#include <alps/alea.h>

#include <alps/config.h>

#include "ietl_interface.h"
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

namespace alps {
namespace diag {

template <typename T, typename M1>
class sparsediag_worker
    : public matrix_worker<T, M1>
{
 public:
    typedef matrix_worker<T, M1> super_type;

    sparsediag_worker(alps::Parameters const& params)
        : super_type(params)
    {}
    void run(alps::ObservableSet& obs) {
        typedef typename boost::numeric::ublas::vector<double> vector_type;

        if (this->progress() >= 1.0) { return; }
        this->is_diagonalized_ = true;

        // measurements
        std::map<std::string, double> m;
        m["Number of Sites"] = this->num_sites();
        m["Volume"] = this->volume();

        // Lanczos diagonalization of Hamiltonian matrix
        vector_type evals;
        {
            typedef ietl::vectorspace<vector_type> vectorspace_type;
            boost::mt19937 generator;
            vectorspace_type vec(this->dimension());
#ifdef ALPS_HAVE_MKL
            typedef typename boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major> matrix2_type;
            matrix2_type Hamiltonian2 = this->matrix();
            ietl::lanczos<matrix2_type, vectorspace_type> lanczos(Hamiltonian2, vec);
#else
            ietl::lanczos<matrix_type, vectorspace_type> lanczos(matrix(), vec);
#endif
            int max_iter = std::min(static_cast<int>(10 * this->dimension()), 1000);
            if (this->params_.defined("MAX_ITERATIONS")) {
                max_iter = this->params_["MAX_ITERATIONS"];
            }
            std::size_t num_eigenvalues = 1;
            if (this->params_.defined("NUMBER_EIGENVALUES")) {
                num_eigenvalues = this->params_["NUMBER_EIGENVALUES"];
            }
            ietl::lanczos_iteration_nlowest<double> iter(max_iter, num_eigenvalues);
            boost::timer t;
            std::clog << "Starting Lanczos... " << std::flush;
            lanczos.calculate_eigenvalues(iter, generator);
            std::clog << "done. Elapsed time: " << t.elapsed() << std::endl;

            int n = std::min(num_eigenvalues, lanczos.eigenvalues().size());
            evals.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                evals[i] = lanczos.eigenvalues()[i];
            }
        }

        double E0 = evals[0];

        m["Ground State Energy"] = E0;
        m["Ground State Energy Density" ] = E0 / this->volume();

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
};


} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_SPARSEDIAG_H_
