SOURCES	= main.cc fulldiag.cc fulldiag.h

bin_PROGRAMS =
if BUILD_DIAG_FULLDIAG_PARAPACK
bin_PROGRAMS += fulldiag-parapack
if HAVE_MPI
bin_PROGRAMS += fulldiag-parapack-mpi
endif
endif

fulldiag_parapack_SOURCES = $(SOURCES)

fulldiag_parapack_mpi_SOURCES = $(SOURCES)
fulldiag_parapack_mpi_LDFLAGS = @LDFLAGS_MPI@
fulldiag_parapack_mpi_LDADD = @LIBS_MPI@
