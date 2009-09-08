SOURCES	= sparsediag.cc sparsediag.h fulldiag.cc fulldiag.h common.h

bin_PROGRAMS =
if BUILD_DIAG_FULLDIAG_PARAPACK
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
if BUILD_DIAG_FULLDIAG_PARAPACK
check_PROGRAMS	+= diag-parapack_ns
endif

diag_parapack_ns_SOURCES	= diag_ns.cc $(SOURCES)

.PHONY : test tests
test tests : check
