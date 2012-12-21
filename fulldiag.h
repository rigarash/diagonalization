/*****************************************************************************
*
* ALPS/diagonalization: Full and Sparse diagonalization
*                       for quantum lattice models using parapack scheduler
*
* Copyright (C) 2003-2009 by Synge Todo <wistaria@comp-phys.org>,
*               2009-2012 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
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

#ifndef ALPS_DIAG_PARAPACK_FULLDIAG_H_
#define ALPS_DIAG_PARAPACK_FULLDIAG_H_

#include "matrix_worker.h"

namespace alps {
namespace diag {

class fulldiag_worker
    : public matrix_worker<double, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> >
{
 public:
    typedef matrix_worker<double, boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> > super_type;
    typedef super_type::vector_type vector_type;

    fulldiag_worker(alps::Parameters const& params)
        : super_type(params),
          eigenvalues()
    {}

    void run_subspace(alps::ObservableSet& /* obs */);
    void run_measurement(alps::ObservableSet& /* obs */) const;

    std::vector<double> eigenvalues;
};

} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_FULLDIAG_H_
