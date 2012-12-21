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

#define DIAG_PARAPACK_VERSION "0.2-20090908"
#define DIAG_PARAPACK_DATE    "2009/09/08"

#define DIAG_PARAPACK_VERSION_STRING "ALPS/diag/parapack version " \
    DIAG_PARAPACK_VERSION " (" DIAG_PARAPACK_DATE ")"

#define DIAG_PARAPACK_COPYRIGHT DIAG_PARAPACK_VERSION_STRING "\n" \
    "  Full and Sparse diagonalization for quantum lattice models using parapack scheduler \n" \
    "  Copyright (c) 1997-2009 by Synge Todo <wistaria@comp-phys.org>,\n" \
    "                             Ryo IGARASHI <rigarash@issp.u-tokyo.ac.jp>\n" \

#include <alps/parapack/parapack.h>

int
main(int argc, char** argv) {
    return alps::parapack::start(argc, argv);
}

PARAPACK_SET_COPYRIGHT(DIAG_PARAPACK_COPYRIGHT)
PARAPACK_SET_VERSION(DIAG_PARAPACK_VERSION_STRING)
