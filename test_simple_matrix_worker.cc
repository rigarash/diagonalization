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

#include "simple_matrix_worker.h"

#include <alps/parameter/parameters.h>
#include <alps/alea/observableset.h>

#include <boost/filesystem.hpp>

#define BOOST_TEST_MODULE test_simple_matrix_worker
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(simple_matrix_worker_h) {
    // prepare Parameters
    alps::Parameters p;
    p["LATTICE"] = "chain lattice";
    p["MODEL"]   = "fermion Hubbard";
    p["L"]       = 4;
    // prepare Observables
    alps::ObservableSet obs;
    obs.reset(true);
    // prepare dumps
    boost::filesystem::path dump(boost::filesystem::unique_path());

    // start test
    {
        alps::diag::simple_matrix_worker<> w1(p);
        w1.init_observables(p, obs);

        // precondition check
        // always thermalized
        BOOST_CHECK(w1.is_thermalized());
        // always finished
        BOOST_CHECK_GE(w1.progress(), 1.0);

        // check run
        w1.run(obs);

        // check save
        alps::OXDRFileDump odp(dump);
        w1.save(odp);
    }

    // restart test
    {
        alps::diag::simple_matrix_worker<> w2(p);
        alps::IXDRFileDump idp(dump);
        w2.load(idp);

        // postcondition check
        // always thermalized
        BOOST_CHECK(w2.is_thermalized());
        // always finished
        BOOST_CHECK_GE(w2.progress(), 1.0);
    }

    // post execution cleanup
    // remove temporary dump files
    BOOST_CHECK(boost::filesystem::remove(dump));
}
