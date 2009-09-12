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

#ifndef ALPS_DIAG_PARAPACK_SPARSEDIAG_H_
#define ALPS_DIAG_PARAPACK_SPARSEDIAG_H_

#include "matrix_worker.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace alps {
namespace diag {

class sparsediag_worker
    : public matrix_worker<double, boost::numeric::ublas::mapped_vector_of_mapped_vector<double, boost::numeric::ublas::row_major> >
{
 public:
    typedef matrix_worker<double, boost::numeric::ublas::mapped_vector_of_mapped_vector<double, boost::numeric::ublas::row_major> > super_type;

    sparsediag_worker(alps::Parameters const& params)
        : super_type(params),
          done(false)
    {}
    double progress() const;
    void run(alps::ObservableSet& /* obs */);
    void save(alps::ODump& /* odump */) const;
    void load(alps::IDump& /* idump */);

 private:
    bool done;
};

} // end namespace diag
} // end namespace alps

#endif // ALPS_DIAG_PARAPACK_SPARSEDIAG_H_
