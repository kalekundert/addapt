AM_CPPFLAGS = -I$(srcdir)/include/
AM_CXXFLAGS = $(OPENMP_CXXFLAGS)
ACLOCAL_AMFLAGS = -I build/m4

# The default flags are 'cru', but the 'u' causes recent versions of ar to 
# issue warnings.
ARFLAGS = cr
AR_FLAGS = cr

# Library
lib_LTLIBRARIES = \
	libaddapt.la

libaddapt_la_SOURCES = \
	src/config.cc \
	src/model.cc \
	src/sampling.cc \
	src/scoring.cc \
	src/utils.cc

pkginclude_HEADERS = \
	include/config.hh \
	include/model.hh \
	include/sampling.hh \
	include/scoring.hh \
	include/utils.hh

# Apps
bin_PROGRAMS = \
	addapt

addapt_SOURCES = \
	apps/addapt.cc

addapt_LDADD = \
	libaddapt.la
	
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
	libaddapt.la

TESTS = $(check_PROGRAMS)


