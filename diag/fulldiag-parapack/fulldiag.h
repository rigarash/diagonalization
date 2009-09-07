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

#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/parapack/serial.h>

// forward declaration
class alps::Parameters;
class alps::ObservableSet;
class alps::ODump;
class alps::IDump;

namespace alps {
namespace diag {

class fulldiag_worker
    : public alps::parapack::abstract_worker,
      protected alps::graph_helper<>,
      protected alps::model_helper<>
{
 public:
    fulldiag_worker(alps::Parameters const& params);
    void init_observables(alps::Parameters const& /* params */,
                          alps::ObservableSet& /* obs */);
    bool is_thermalized() const;
    double progress() const;
    void run(alps::ObservableSet& /* obs */);
    void save(alps::ODump& /* odump */) const;
    void load(alps::IDump& /* idump */);

 private:
    bool done;
};

} // end namespace diag
} // end namespace alps
