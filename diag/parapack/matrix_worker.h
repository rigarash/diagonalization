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

#ifndef ALPS_DIAG_PARAPACK_MATRIX_WORKER_H_
#define ALPS_DIAG_PARAPACK_MATRIX_WORKER_H_

#include "common.h"

#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/parameter.h>
#include <alps/osiris.h>
#include <alps/alea.h>

#include <alps/parapack/serial.h>

#include <boost/timer.hpp>

namespace alps {
namespace diag {

template <typename T, typename M, typename Mtmp = M>
class matrix_worker
    : public alps::parapack::abstract_worker,
      protected alps::graph_helper<>,
      protected alps::model_helper<>
{
 public:
    typedef T value_type;
    typedef M matrix_type;
    typedef Mtmp matrix_build_type;

 public:
    typedef alps::basis_states<short> basis_states_type;
    typedef basis_states_type::value_type state_type;

    matrix_worker(alps::Parameters const& ps)
        : alps::parapack::abstract_worker(),
          alps::graph_helper<>(ps),
          alps::model_helper<>(*this, ps),
          params_(ps),
          build_matrix_(),
          built_basis_(false),
          built_matrix_(false),
          is_diagonalized_(false)
    {}
    void init_observables(alps::Parameters const& /* params */,
                          alps::ObservableSet& /* obs */) {}
    bool is_thermalized() const { return true; }
    double progress() const {
        return (built_basis_ && built_matrix_ && is_diagonalized_) ? 1.0 : 0.0;
    }
    virtual void run(alps::ObservableSet& /* obs */) = 0;
    void save(alps::ODump& odump) const
    { odump << is_diagonalized_; }
    void load(alps::IDump& idump)
    { idump >> is_diagonalized_; }

    std::size_t dimension() const {
        return basis_states().size();
    }

    basis_states_type& basis_states() {
        if (!built_basis_) {
            build_basis();
        }
        return states_;
    }
    basis_states_type const& basis_states() const {
        if (!built_basis_) {
            build_basis();
        }
        return states_;
    }
    matrix_type& matrix() {
        if (!built_matrix_) {
            build_matrix_(*this);
        }
        return matrix_;
    }
    matrix_type const& matrix() const {
        if (!built_matrix_) {
            build_matrix_(*this);
        }
        return matrix_;
    }

 protected:
    alps::Parameters params_;
    mutable alps::basis_states_descriptor<short> basis_;
    mutable basis_states_type states_;
    mutable matrix_type matrix_;

    void build_basis() const {
        boost::timer t;
        std::clog << "Building basis... " << std::flush;
        basis_  = alps::basis_states_descriptor<short>(basis(), graph());
        states_ = basis_states_type(basis_);
        std::clog << "done. Elapsed time: " << t.elapsed() << std::endl;
        built_basis_ = true;
    }

    template <typename T_, typename M_, typename Mtmp_>
    class build_matrix {
     public:
        void operator() (matrix_worker<T_, M_, Mtmp_> const& p) const {
            typedef matrix_worker<T, M, M> P;
            boost::timer t;
            std::clog << "Building matrix..." << std::flush;
            typename P::matrix_build_type mtmp = P::matrix_build_type(p.dimension(), p.dimension());
            mtmp.clear();
            BOOST_FOREACH(site_descriptor s, p.sites()) {
                add_to_matrix(mtmp, p.model(), p.basis_states(), s, p.graph(), p.params_);
            }
            BOOST_FOREACH(bond_descriptor b, p.bonds()) {
                add_to_matrix(mtmp, p.model(), p.basis_states(), b, p.graph(), p.params_);
            }
            std::clog << "done. Elapsed time: " << t.elapsed() << std::endl;
            t.restart();
            std::clog << "Copying matrix... " << std::flush;
            p.matrix_ = mtmp;
            std::clog << "done. Elapsed time: " << t.elapsed() << std::endl;
            p.built_matrix_ = true;

        }
    };
    template <typename T_, typename M_>
    class build_matrix<T_, M_, M_> {
     public:
        void operator() (matrix_worker<T_, M_, M_> const& p) const {
            typedef matrix_worker<T, M, M> P;
            boost::timer t;
            std::clog << "Building matrix..." << std::flush;
            p.matrix_ = matrix_type(p.dimension(), p.dimension());
            p.matrix_.clear();
            BOOST_FOREACH(site_descriptor s, p.sites()) {
                add_to_matrix(p.matrix_, p.model(), p.basis_states(), s, p.graph(), p.params_);
            }
            BOOST_FOREACH(bond_descriptor b, p.bonds()) {
                add_to_matrix(p.matrix_, p.model(), p.basis_states(), b, p.graph(), p.params_);
            }
            std::clog << "done. Elapsed time: " << t.elapsed() << std::endl;
            p.built_matrix_ = true;
        }
    };
    build_matrix<T, M, Mtmp> const build_matrix_;

 private:
    mutable bool built_basis_;
    mutable bool built_matrix_;

 protected:
    bool is_diagonalized_;
};

} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_MATRIX_WORKER_H_
