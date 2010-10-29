LDADD       = @top_builddir@/pytbeaglehon/ccore/libfatbeagle.la
AM_CPPFLAGS = -I@top_srcdir@
bin_PROGRAMS = testfatbeagle testtrace
testfatbeagle_SOURCES = testfatbeagle.c
testfatbeagle_CPPFLAGS = $(AM_CPPFLAGS)

check-local: testfatbeagle
	$(abs_top_builddir)/c_example/testfatbeagle -v


testtrace_SOURCES = tracetest.c
testtrace_CPPFLAGS = $(AM_CPPFLAGS)

check-local: testtrace
	$(abs_top_builddir)/c_example/testtrace -v
