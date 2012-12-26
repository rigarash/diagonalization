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

#define SIMPLE_MATRIX_WORKER_VERSION "0.1"
#define SIMPLE_MATRIX_WORKER_DATE    "2012/12/25"

#define SIMPLE_MATRIX_WORKER_VERSION_STRING "ALPS/diag/simple_matrix_worker version " \
    SIMPLE_MATRIX_WORKER_VERSION " (" SIMPLE_MATRIX_WORKER_DATE ")"

#define SIMPLE_MATRIX_WORKER_COPYRIGHT SIMPLE_MATRIX_WORKER_VERSION_STRING "\n" \
    "  simple diagonalization for quantum lattice models using parapack scheduler \n" \
    "  Copyright (c) 2012-2012 Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>\n" \

PARAPACK_SET_VERSION(SIMPLE_MATRIX_WORKER_VERSION);
PARAPACK_SET_COPYRIGHT(SIMPLE_MATRIX_WORKER_COPYRIGHT);
PARAPACK_REGISTER_WORKER(alps::diag::simple_matrix_worker<>, "simple diagonalization");
PARAPACK_REGISTER_EVALUATOR(alps::parapack::simple_evaluator, "simple diagonalization");
