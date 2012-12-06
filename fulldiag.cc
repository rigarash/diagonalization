/*****************************************************************************
*
* ALPS/diag/parapack: Full and Sparse diagonalization
*                     for quantum lattice models using parapack scheduler
*
* Copyright (C) 2003-2009 by Synge Todo <wistaria@comp-phys.org>,
*                            Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
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
#include "fulldiag.h"

#include <alps/parameter.h>
#include <alps/alea.h>

#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/bindings/lower.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/lapack.hpp>

#include <cmath>
#include <cstddef>

namespace {

template <typename T, typename R, typename A, typename V>
void
diagonalize(
    boost::numeric::ublas::matrix<T, R, A>& a,
    V& v, bool need_eigenvectors = true)
{
    using namespace boost::numeric::bindings::lapack;

    BOOST_STATIC_ASSERT((boost::is_same<typename R::orientation_category, boost::numeric::ublas::column_major_tag>::value));

    const char jobz = (need_eigenvectors ? 'V' : 'N');

    std::ptrdiff_t info = syev(jobz, boost::numeric::bindings::lower(a), v);

    if (info != 0) {
        throw std::runtime_error("failed in syev");
    }
}

} // end namespace

namespace alps{
namespace diag{

void
fulldiag_worker::run_subspace(alps::ObservableSet& obs)
{
    typedef boost::numeric::ublas::vector<double> diagonal_matrix_type;

    if (progress() >= 1.0) { return; }
    this->is_diagonalized_ = true;

    // partition function coefficients
    // matrix_type part(dimension(), dimension());
    // part.clear();
    // for (std::size_t i = 0; i < dimension(); ++i) {
    //     part(i, i) = 1;
    // }
    // double factor = 1.;
    // for (std::size_t k = 0; k < 5; ++k) {
    //     if (k != 0) {
    //         part = prod(matrix(), part);
    //         factor *= -k;
    //     }
    //     double trace = 0.;
    //     for (std::size_t i = 0; i < dimension(); ++i) {
    //         trace += part(i, i);
    //     }
    //     m["Partition Function Coefficient #" + boost::lexical_cast<std::string>(k)]
    //         = trace / factor;
    // }

    // diagonalize Hamiltonian matrix
    vector_type evals(dimension());
    std::cerr << "Start Diagonalization...\n" << std::flush;
    diagonalize(matrix(), evals);
    std::cerr << "Done Diagonalization.\n" << std::flush;

    BOOST_FOREACH(double d, evals) {
        eigenvalues.push_back(d);
    }
}

void
fulldiag_worker::run_measurement(alps::ObservableSet& obs) const {
    // sort eigenvalues
    std::vector<double> evv(eigenvalues);
    std::sort(evv.begin(), evv.end());
    vector_type evals(evv.size());
    for (std::size_t i = 0; i < evv.size(); ++i) {
        evals(i) = evv[i];
    }

    double beta = 1.0;
    if (params_.defined("T")) {
        beta = 1.0 / alps::evaluate<double>("T", params_);
    }

    // measurements
    std::map<std::string, double> m;
    m["Temperature"] = 1.0 / beta;
    m["Inverse Temperature"] = beta;
    m["Number of Sites"] = num_sites();
    m["Volume"] = volume();

    double Z  = 0.;
    double E0 = *std::min_element(eigenvalues.begin(), eigenvalues.end());
    BOOST_REVERSE_FOREACH(double eval, evv) {
        double weight = std::exp(-beta * (eval - E0));
        Z += weight;
    }

    double ene, ene2;
    boost::tie(ene, ene2) = static_average2(beta, evals);
    ene  = ene  / Z;
    ene2 = ene2 / Z / std::pow(std::abs(volume()), 2);
    double fene = E0 - std::log(Z) / beta;
    m["Ground State Energy"] = E0;
    m["Ground State Energy Densty" ] = E0 / volume();
    m["Energy"] = ene;
    m["Energy Density"] = ene / volume();
    m["Specific Heat"] = beta * beta * volume() * (ene2 - std::pow(std::abs(ene / volume()), 2));
    m["Free Energy"] = fene;
    m["Free Energy Density"] = fene / volume();
    m["Entropy"] = beta * (ene - fene);
    m["Entropy Density"] = beta * (ene - fene) / volume();

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

} // end namespace diag
} // end namespace alps

//namespace {

template class alps::diag::matrix_worker<double, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> >;

typedef alps::diag::matrix_evaluator<double, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> > fulldiag_evaluator;

PARAPACK_REGISTER_WORKER(alps::diag::fulldiag_worker, "Full diagonalization");
PARAPACK_REGISTER_EVALUATOR(fulldiag_evaluator, "Full diagonalization");

//} // end namespace
