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

#include <alps/lattice.h>
#include <alps/parapack/worker_factory.h>

namespace alps {
namespace diag {

class simple_matrix_worker
    : public alps::parapack::abstract_worker,
      protected alps::graph_helper<>
{
 public:
    simple_matrix_worker(alps::Parameters const& ps)
        : alps::parapack::abstract_worker(),
          alps::graph_helper<>(ps)
    {}
    // TODO: implement
    void init_observables(alps::Parameters const& /* ps */,
                          alps::ObservableSet& /* obs */)
    {}
    // No need for thermalization for diagonalization
    bool is_thermalized() const { return true; }
    // TODO: implement
    double progress() const {
        return 1.0;
    }
    void run(alps::ObservableSet& obs) {
        if (progress() >= 1.0) { return; }
    }
    void load(alps::IDump& dp) {
    }
    void save(alps::ODump& dp) const {
    }
};

} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_SIMPLE_MATRIX_WORKER_H_







