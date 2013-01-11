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

#include "types.h"

#define BOOST_TEST_MODULE test_types
#ifndef ALPS_LINK_BOOST_TEST
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

using namespace alps::diag;

BOOST_AUTO_TEST_CASE(types_h) {
    worker_status_t status(worker_status::Undefined);
    BOOST_CHECK_EQUAL(status, worker_status::Undefined);
    BOOST_CHECK_CLOSE(worker_status::progress(status), 0.0, 1e-4);
    status = worker_status::Ready;
    BOOST_CHECK_EQUAL(status, worker_status::Ready);
    BOOST_CHECK_CLOSE(worker_status::progress(status), 0.1, 1e-4);
    status = worker_status::Finished;
    BOOST_CHECK_EQUAL(status, worker_status::Finished);
    BOOST_CHECK_CLOSE(worker_status::progress(status), 1.0, 1e-4);
}
