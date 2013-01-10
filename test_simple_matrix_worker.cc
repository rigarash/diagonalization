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
#ifndef ALPS_HAS_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

struct F {
    F()
        : pp(),
          wp()
    {
        // prepare Parameters
        pp.reset(new alps::Parameters);
        (*pp)["LATTICE"] = "chain lattice";
        (*pp)["MODEL"]   = "fermion Hubbard";
        (*pp)["L"]       = 4;

        // prepare worker
        wp.reset(new alps::diag::simple_matrix_worker<matrix_type>(*pp));

        // precondition check
        BOOST_CHECK_EQUAL(wp->status(), alps::diag::worker_status::Undefined);
        BOOST_CHECK_CLOSE(wp->progress(), 0.0, 1e-4);

        // always thermalized
        BOOST_CHECK(wp->is_thermalized());
    }
    ~F()
    {
        // postcondition check
        // always thermalized
        BOOST_CHECK(wp->is_thermalized());
    }

    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
    boost::scoped_ptr<alps::Parameters> pp;
    boost::scoped_ptr<alps::diag::simple_matrix_worker<matrix_type> > wp;
};

BOOST_AUTO_TEST_CASE(simple_matrix_worker_h) {
    // prepare Observables
    alps::ObservableSet obs;
    obs.reset(true);
    // prepare dumps
    boost::filesystem::path dump(boost::filesystem::unique_path());

    // start test
    {
        F f;
        // check initialization
        f.wp->init_observables(*(f.pp), obs);
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Ready);
        BOOST_CHECK_CLOSE(f.wp->progress(), 0.1, 1e-4);

        // check run
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Ready);
        BOOST_CHECK_CLOSE(f.wp->progress(), 0.1, 1e-4);
        f.wp->run(obs);
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Finished);
        BOOST_CHECK_GE(f.wp->progress(), 1.0);

        // check save
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Finished);
        BOOST_CHECK_GE(f.wp->progress(), 1.0);
        alps::OXDRFileDump odp(dump);
        f.wp->save(odp);
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Finished);
        BOOST_CHECK_GE(f.wp->progress(), 1.0);
    }

    // restart test
    {
        F f;
        // check load
        alps::IXDRFileDump idp(dump);
        f.wp->load(idp);
        BOOST_CHECK_EQUAL(f.wp->status(), alps::diag::worker_status::Finished);
        BOOST_CHECK_GE(f.wp->progress(), 1.0);

        // postcondition check
        BOOST_CHECK_GE(f.wp->progress(), 1.0);
    }

    // post execution cleanup
    // remove temporary dump files
    BOOST_CHECK(boost::filesystem::remove(dump));
}
