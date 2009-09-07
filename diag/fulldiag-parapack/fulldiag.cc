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

namespace alps{
namespace diag{

fulldiag_worker::fulldiag_worker(alps::Parameters const& params)
    : alps::parapack::abstract_worker(),
      alps::graph_helper<>(params),
      alps::model_helper<>(*this, params),
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
    done = true;
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
