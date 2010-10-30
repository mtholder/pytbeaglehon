LDADD       = @top_builddir@/pytbeaglehon/ccore/libfatbeagle.la
AM_CPPFLAGS = -I@top_srcdir@
noinst_PROGRAMS = testfatbeagle pytbeaglehontrace beagletrace

testfatbeagle_SOURCES = testfatbeagle.c
testfatbeagle_CPPFLAGS = $(AM_CPPFLAGS)

check-local: testfatbeagle
	$(abs_top_builddir)/c_example/testfatbeagle -v


nodist_pytbeaglehontrace_SOURCES = pytbeaglehontrace.c
nodist_beagletrace_SOURCES = beagletrace.c
BUILT_SOURCES = pytbeaglehontrace.c beagletrace.c $(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c $(abs_top_srcdir)/c_example/real_beagletrace.c

testtrace_CPPFLAGS = $(AM_CPPFLAGS)

check-local: testtrace
	$(abs_top_builddir)/c_example/testtrace -v

$(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c: Makefile
	if test -f $(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c ; then \
	echo ; \
	else \
	echo 'int main(int argc, char * argv[]) { return 0; }' > $(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c  ; \
	fi
pytbeaglehontrace.c: Makefile $(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c
	cp $(abs_top_srcdir)/c_example/real_pytbeaglehontrace.c pytbeaglehontrace.c

$(abs_top_srcdir)/c_example/real_beagletrace.c: Makefile
	if test -f $(abs_top_srcdir)/c_example/real_beagletrace.c ; then \
	echo ; \
	else \
	echo 'int main(int argc, char * argv[]) { return 0; }' > $(abs_top_srcdir)/c_example/real_beagletrace.c  ; \
	fi
beagletrace.c: Makefile $(abs_top_srcdir)/c_example/real_beagletrace.c
	cp $(abs_top_srcdir)/c_example/real_beagletrace.c beagletrace.c

