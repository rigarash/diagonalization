dist_noinst_HEADERS = sparsediag.h fulldiag.h common.h ietl_interface.h matrix_worker.h

SOURCES	= sparsediag.cc fulldiag.cc $(dist_noinst_HEADERS)

bin_PROGRAMS =
if BUILD_DIAG_PARAPACK
bin_PROGRAMS += diag-parapack
if HAVE_MPI
bin_PROGRAMS += diag-parapack-mpi
endif
endif

diag_parapack_SOURCES = main.cc $(SOURCES)

diag_parapack_mpi_SOURCES = main.cc $(SOURCES)
diag_parapack_mpi_LDFLAGS = @LDFLAGS_MPI@
diag_parapack_mpi_LDADD = @LIBS_MPI@

check_PROGRAMS	=
if BUILD_DIAG_PARAPACK
check_PROGRAMS	+= diag-parapack_ns
endif

dist_noinst_DATA = fulldiag_ns.ip sparsediag_ns.ip

diag_parapack_ns_SOURCES	= diag_ns.cc $(SOURCES)

.PHONY : test tests
test tests : check
