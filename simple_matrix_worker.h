/*****************************************************************************
*
* ALPS/diagonalization: Full and Sparse diagonalization
*                       for quantum lattice models using parapack scheduler
*
* Copyright (C) 2012-2012 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
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

#ifndef ALPS_DIAG_PARAPACK_SIMPLE_MATRIX_WORKER_H_
#define ALPS_DIAG_PARAPACK_SIMPLE_MATRIX_WORKER_H_

#include "types.h"

#include <alps/lattice/graph_helper.h>
#include <alps/model/model_helper.h>
#include <alps/parapack/worker_factory.h>

#include <boost/smart_ptr/scoped_ptr.hpp>

#include <stdexcept>

namespace alps {
namespace diag {

template <typename G = typename alps::graph_helper<>::graph_type>
class simple_matrix_worker
    : public alps::parapack::abstract_worker
{
 public:
    typedef alps::graph_helper<G> graph_helper_type;
    typedef alps::model_helper<>  model_helper_type;

    typedef alps::basis_states<short> basis_states_type;

 public:
    simple_matrix_worker(alps::Parameters const& ps)
        : alps::parapack::abstract_worker(),
          graph_(ps),
          model_(graph_, ps),
          states_ptr(),
          status_(worker_status::Undefined)
    {}
    // TODO: implement
    void init_observables(alps::Parameters const& /* ps */,
                          alps::ObservableSet& /* obs */) {
        status_ = worker_status::Ready;
    }
    // No need for thermalization for diagonalization
    bool is_thermalized() const { return true; }
    double progress() const {
        return worker_status::progress(status_);
    }
    void run(alps::ObservableSet& obs) {
        switch(status_) {
        case worker_status::Ready:
            build_basis_states();
            status_ = worker_status::Finished;
            break;
        default:
            throw std::invalid_argument("Invalid status code");
        }
    }
    void load(alps::IDump& dp) {
        int status_tmp;
        dp >> status_tmp;
        status_ = static_cast<worker_status_t>(status_tmp);
    }
    void save(alps::ODump& dp) const {
        int status_tmp = static_cast<int>(status_);
        dp << status_tmp;
    }

    worker_status_t status() const {
        return status_;
    }

    basis_states_type& basis_states() {
        if (!states_ptr) {
            build_basis_states();
        }
        return *states_ptr;
    }
    basis_states_type const& basis_states() const {
        if (!states_ptr) {
            throw std::logic_error("basis states are not built.");
        }
        return *states_ptr;
    }

    void build_basis_states() {
        if (!states_ptr) {
            states_ptr.reset(new basis_states_type(alps::basis_states_descriptor<short>(model_.basis(), graph_.graph())));
        }
    }

 private:
    graph_helper_type graph_;
    model_helper_type model_;
    boost::scoped_ptr<basis_states_type> states_ptr;

    worker_status_t status_;
};

} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_SIMPLE_MATRIX_WORKER_H_
