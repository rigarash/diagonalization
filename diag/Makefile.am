SUBDIRS =
if EXIST_DIAG_FULLDIAG
SUBDIRS += fulldiag
endif
if EXIST_DIAG_SPARSEDIAG
SUBDIRS += sparsediag
endif
if EXIST_DIAG_PARAPACK
SUBDIRS += parapack
endif
SUBDIRS += .
DIST_SUBDIRS = $(SUBDIRS)

noinst_HEADERS = diag.h matrix.h
