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

#define FULLDIAG_PARAPACK_VERSION "0.1-20090907"
#define FULLDIAG_PARAPACK_DATE    "2009/09/07"

#define FULLDIAG_PARAPACK_VERSION_STRING "ALPS/fulldiag-parapack version " \
    FULLDIAG_PARAPACK_VERSION " (" FULLDIAG_PARAPACK_DATE ")"

#define FULLDIAG_PARAPACK_COPYRIGHT FULLDIAG_PARAPACK_VERSION_STRING "\n" \
    "  Full diagonalization for quantum lattice systems using parapack scheduler \n" \
    "  Copyright (c) 1997-2009 by Synge Todo <wistaria@comp-phys.org>,\n" \
    "                             Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>\n" \

#include <alps/parapack/scheduler.h>

int
main(int argc, char** argv) {
    return alps::parapack::start(argc, argv);
}

PARAPACK_SET_COPYRIGHT(FULLDIAG_PARAPACK_COPYRIGHT)
PARAPACK_SET_VERSION(FULLDIAG_PARAPACK_VERSION_STRING)
