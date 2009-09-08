SOURCES	= fulldiag.cc fulldiag.h common.h

bin_PROGRAMS =
if BUILD_DIAG_FULLDIAG_PARAPACK
bin_PROGRAMS += fulldiag-parapack
if HAVE_MPI
bin_PROGRAMS += fulldiag-parapack-mpi
endif
endif

fulldiag_parapack_SOURCES = main.cc $(SOURCES)

fulldiag_parapack_mpi_SOURCES = main.cc $(SOURCES)
fulldiag_parapack_mpi_LDFLAGS = @LDFLAGS_MPI@
fulldiag_parapack_mpi_LDADD = @LIBS_MPI@

check_PROGRAMS	=
if BUILD_DIAG_FULLDIAG_PARAPACK
check_PROGRAMS	+= fulldiag-parapack_ns
endif

fulldiag_parapack_ns_SOURCES	= fulldiag_ns.cc $(SOURCES)

.PHONY : test tests
test tests : check
