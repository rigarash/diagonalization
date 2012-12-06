/*****************************************************************************
*
* ALPS/diag/parapack: Full and Sparse diagonalization
*                     for quantum lattice models using parapack scheduler
*
* Copyright (C) 2009-2009 by Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>
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

#include "sparsediag.h"

#include <alps/config.h>

//namespace {

typedef double value_type;
typedef boost::numeric::ublas::mapped_vector_of_mapped_vector<value_type, boost::numeric::ublas::row_major> matrix_type;
typedef boost::numeric::ublas::compressed_matrix<value_type, boost::numeric::ublas::row_major> matrix_tmp_type;

template class alps::diag::matrix_worker<value_type, matrix_type>;

#ifdef ALPS_HAVE_MKL
typedef alps::diag::sparsediag_worker<value_type, matrix_type, matrix_tmp_type> sparsediag;
template class alps::diag::sparsediag_worker<value_type, matrix_type, matrix_tmp_type>;
#else
typedef alps::diag::sparsediag_worker<value_type, matrix_type> sparsediag;
template class alps::diag::sparsediag_worker<value_type, matrix_type>;
#endif

PARAPACK_REGISTER_WORKER(sparsediag, "Sparse diagonalization");

//} // end namespace
