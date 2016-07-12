AM_CPPFLAGS = "-I$(srcdir)/include/"
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
ACLOCAL_AMFLAGS = -I build/m4

# The default flags are 'cru', but the 'u' causes recent versions of ar to 
# issue warnings.
ARFLAGS = cr
AR_FLAGS = cr

# Library
lib_LTLIBRARIES = \
	libsgrna.la

libsgrna_la_SOURCES = \
	src/model.cc \
	src/scoring.cc \
	src/sampling.cc \
	src/utils.cc

pkginclude_HEADERS = \
	include/model.hh \
	include/scoring.hh \
	include/sampling.hh \
	include/utils.hh

# Apps
bin_PROGRAMS = \
	mh \
	addapt

mh_SOURCES = \
	apps/mh.cc

mh_LDADD = \
	libsgrna.la

addapt_SOURCES = \
	apps/addapt.cc

addapt_LDADD = \
	libsgrna.la
	
# Tests
check_PROGRAMS = \
	run_tests

run_tests_CXXFLAGS = \
	-g

run_tests_SOURCES = \
	tests/test_model.cc \
	tests/test_scoring.cc \
	tests/test_sampling.cc \
	tests/main.cc

run_tests_LDADD = \
	libsgrna.la

TESTS = $(check_PROGRAMS)


