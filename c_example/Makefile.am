LDADD       = @top_builddir@/pytbeaglehon/ccore/libfatbeagle.la
AM_CPPFLAGS = -I@top_srcdir@
bin_PROGRAMS = testfatbeagle
testfatbeagle_SOURCES = testfatbeagle.c
testfatbeagle_CPPFLAGS = $(AM_CPPFLAGS)

check-local: testfatbeagle
	$(abs_top_builddir)/c_example/testfatbeagle -v
