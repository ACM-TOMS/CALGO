#=======================================================================
# Makefile for files in /u/sy/beebe/jim-ball/toms/dist/laguerre/
#
# Current target list:
#	all			make everything
#	cglfd1			C test program
#	cglfd1.o
#	check			run complete validation suite
#	check-cglfd1		run C program glfd validation suite
#	check-glf		run glf validation suites
#	check-glf1
#	check-glf2
#	check-glf3
#	check-glf4
#	check-glf5
#	check-glf6
#	check-glfd		run glfd validation suites
#	check-glfd1
#	check-glfd2
#	check-glfd3
#	check-glfd4
#	check-glfd5
#	check-glfd6
#	check-qglf		run qglf validation suites
#	check-qglf1
#	check-qglf2
#	check-qglf3
#	check-qglf4
#	check-qglf5
#	check-qglf6
#	check-qglfd		run qglfd validation suites
#	check-qglfd1
#	check-qglfd2
#	check-qglfd3
#	check-qglfd4
#	check-qglfd5
#	check-qglfd6
#	check-xxx		internal target: do not use
#	clean			remove unneeded intermediate files
#	clean-test		remove unneeded test output files
#	clobber			same as distclean
#	configure		make patched configure script
#	demglq.tab		make data for Table 1 of paper 1
#	demo			run complete demonstration suite
#	demo-dp			run double-precision demonstration
#	demo-dp-1
#	demo-dp-2
#	demo-lapack		run demonstration suite for LAPACK eigensolvers
#	demo-qp			run quadruple-precision demonstration
#	demo-qp-1
#	demo-qp-2
#	demo-setup		internal target
#	distclean		remove almost everything make can rebuild
#				(paring back to the original distribution)
#	FRC			dummy target to force another target to
#				be remade
#	glqf.o
#	glqfd.o
#	glqli.o
#	install			install files on local system
#	lib			add object files to library
#	maintainer-clean	remove absolutely everything make can rebuild
#				at the author's site (do NOT use this target)
#	mostlyclean		same as clean
#	tglf1			test program for glqf
#	tglf1.o
#	tglf2
#	tglf2.o
#	tglf3
#	tglf3.o
#	tglf4
#	tglf4.o
#	tglf5
#	tglf5.o
#	tglf6
#	tglf6.o
#	tglfd.o
#	tglfd1			test program for glqfd
#	tglfd1.o
#	tglfd2
#	tglfd2.o
#	tglfd3
#	tglfd3.o
#	tglfd4
#	tglfd4.o
#	tglfd5
#	tglfd5.o
#	tglfd6
#	tglfd6.o
#	uninstall		remove files installed by "make install"
#
# [27-Mar-2000]
#=======================================================================

#=======================================================================
# Macro definitions

ABSERR		= 0

AUTOCONF	= autoconf

AR		= ar

AWK		= awk

CC		= gcc

CHMOD		= chmod

CMP		= cmp

CP		= cp

DIFF		= diff

FC		= f77

LAPACK		= -L$(prefix)/lib -llapack

LS		= ls

MKDIR		= mkdir

MV		= mv

# This script will use ndiff if it is available in the search path,
# and will otherwise fall back to using diff.
NDIFF		= `$(WHICH) ndiff`

NDIFFFLAGS	= --relerr $(RELERR) --minwidth 10 --abserr $(ABSERR)

RANLIB		= ranlib

# This is the relative error tolerance used in the check target.
# It can be redefined on the make command line if desired.
RELERR		= 1.0e-10

RM		= rm -f

SED		= sed

SHELL		= /bin/sh

TEST		= test

# We use our own private version of the `which' command, because
# system implementations of `which' have inconsistent behavior.
WHICH		= $(COMMON)/which.sh

#=======================================================================
# Compiler flags and source files

CDEFS		=

CEXTRA		= 

CFLAGS		= -DHAVE_SUNMATH_H $(COPT) $(CINC) $(CWARN) $(CDEFS) $(CEXTRA)

CDEFS		=

CINC		= -I..

CGLFD1OBJS	= cglfd1.o glqli.o glqfd.o glqrc.o $(COMMONOBJS)

CQGLFD1OBJS	= cqglfd1.o qglqli.o qglqfd.o qglqrc.o $(QCOMMONOBJS)

CHECKFILES	= unknown

CHECKGLF1       = glf101 glf102 glf103 glf104 glf105 glf106 glf107 glf108 \
		  glf109

CHECKGLF2	= glf201 glf202

CHECKGLF3	= glf301 glf302

CHECKGLF4       = glf401 glf402 glf403 glf404 glf405 glf406 glf407 glf408 \
		  glf409 glf410

CHECKGLF5       = glf501 glf502 glf503 glf504 glf505 glf506 glf507 glf508 \
		  glf509

CHECKGLF6       = glf601 glf602 glf603 glf604 glf605 glf606 glf607 glf608 \
		  glf609 glf610

CHECKGLFD1      = glfd101 glfd102 glfd103 glfd104 glfd105 glfd106 glfd107 \
		  glfd108 glfd109

CHECKGLFD2	= glfd201 glfd202

CHECKGLFD3	= glfd301 glfd302

CHECKGLFD4	= glfd401 glfd402 glfd403 glfd404 glfd405 glfd406 glfd407 \
		  glfd408 glfd409 glfd410

CHECKGLFD5	= glfd501 glfd502 glfd503 glfd504 glfd505 glfd506 glfd507 \
		  glfd508 glfd509

CHECKGLFD6	= glfd601 glfd602 glfd603 glfd604 glfd605 glfd606 glfd607 \
		  glfd608 glfd609 glfd610

CHECKCGLFD1	= cglfd101 cglfd102 cglfd103 cglfd104 cglfd105 cglfd106 \
		  cglfd107 cglfd108 cglfd109

CHECKQGLF1      = qglf101 qglf102 qglf103 qglf104 qglf105 qglf106 qglf107 \
		  qglf108 qglf109

CHECKQGLF2	= qglf201 qglf202

CHECKQGLF3	= qglf301 qglf302

CHECKQGLF4      = qglf401 qglf402 qglf403 qglf404 qglf405 qglf406 qglf407 \
		  qglf408 qglf409 qglf410

CHECKQGLF5      = qglf501 qglf502 qglf503 qglf504 qglf505 qglf506 qglf507 \
		  qglf508 qglf509

CHECKQGLF6      = qglf601 qglf602 qglf603 qglf604 qglf605 qglf606 qglf607 \
		  qglf608 qglf609 qglf610

CHECKQGLFD1     = qglfd101 qglfd102 qglfd103 qglfd104 qglfd105 qglfd106 \
		  qglfd107 qglfd108 qglfd109

CHECKQGLFD2	= qglfd201 qglfd202

CHECKQGLFD3	= qglfd301 qglfd302

CHECKQGLFD4     = qglfd401 qglfd402 qglfd403 qglfd404 qglfd405 qglfd406 \
		  qglfd407 qglfd408 qglfd409 qglfd410

CHECKQGLFD5     = qglfd501 qglfd502 qglfd503 qglfd504 qglfd505 qglfd506 \
		  qglfd507 qglfd508 qglfd509

CHECKQGLFD6     = qglfd601 qglfd602 qglfd603 qglfd604 qglfd605 qglfd606 \
		  qglfd607 qglfd608 qglfd609 qglfd610

CHECKCQGLFD1    = cqglfd101 cqglfd102 cqglfd103 cqglfd104 cqglfd105 \
		  cqglfd106 cqglfd107 cqglfd108 cqglfd109

CLIBS		= $(LIBRARY)  -lm

COMMON		= ../common

COMMONOBJS      = $(COMMON)/d2norm.o $(COMMON)/dcopy.o $(COMMON)/deps.o \
		  $(COMMON)/derbit.o $(COMMON)/dfloat.o $(COMMON)/dgamma.o \
		  $(COMMON)/dinf.o $(COMMON)/dlog2.o $(COMMON)/dnan.o \
		  $(COMMON)/dpsi.o $(COMMON)/dstore.o $(COMMON)/dstorf.o \
		  $(COMMON)/dvsum.o $(COMMON)/isdnan.o $(COMMON)/pythag.o \
		  $(COMMON)/tql1.o $(COMMON)/trapit.o

COPT		= -O

CTESTPROGS	= cglfd1

CWARN		=

DEMOPROGS       = demglq demglq-eqr demglq-erf demglq-ev demglq-evd \
		  demglq-evr demglq-evx demglq-tql1

FEXTRA		= 

FFLAGS		= -u $(FOPT) $(FINC) $(FWARN) $(FEXTRA)

FINC		=  -I$(COMMON)

FOPT		= -O

FWARN		=

GLFOBJS		= glqf.o $(GLFDOBJS)

GLFDOBJS	= glqli.o glqfd.o glqrc.o $(COMMON)/prthdr.o $(COMMONOBJS)

LDFLAGS		= 

LIBNAME		= libgjl

LIBOBJS		= glqf.o glqfd.o glqli.o glqrc.o

LIBRARY		= ../$(LIBNAME).a

LIBS		= $(LIBRARY)

QCOMMONOBJS	= $(COMMON)/q2norm.o $(COMMON)/qcopy.o $(COMMON)/qeps.o \
		  $(COMMON)/qerbit.o $(COMMON)/qfloat.o $(COMMON)/qgamma.o \
		  $(COMMON)/qinf.o $(COMMON)/qlog2.o $(COMMON)/qnan.o \
		  $(COMMON)/qpsi.o $(COMMON)/qstore.o $(COMMON)/qstorf.o \
		  $(COMMON)/qvsum.o $(COMMON)/isqnan.o $(COMMON)/qpythg.o \
		  $(COMMON)/qtql1.o $(COMMON)/trapit.o

QCTESTPROGS	= cqglfd1

QDEMOPROGS      = demqglq

QGLFOBJS	= qglqf.o $(QGLFDOBJS)

QGLFDOBJS	= qglqli.o qglqfd.o qglqrc.o $(COMMON)/qprthdr.o $(QCOMMONOBJS)

QLIBOBJS	= qglqf.o qglqfd.o qglqli.o qglqrc.o

QTESTPROGS      = tqglf1 tqglf2 tqglf3 tqglf4 tqglf5 tqglf6 tqglfd1 tqglfd2 \
		  tqglfd3 tqglfd4 tqglfd5 tqglfd6

TESTNAME	= unknown

TESTPROGS       = tglf1 tglf2 tglf3 tglf4 tglf5 tglf6 tglfd1 tglfd2 tglfd3 \
		  tglfd4 tglfd5 tglfd6

TGLF1OBJS	= tglf1.o $(GLFOBJS)

TGLF2OBJS	= tglf2.o $(GLFOBJS)

TGLF3OBJS	= tglf3.o $(GLFOBJS)

TGLF4OBJS	= tglf4.o $(GLFOBJS)

TGLF5OBJS	= tglf5.o $(GLFOBJS)

TGLF6OBJS	= tglf6.o $(GLFOBJS)

TGLFD1OBJS	= tglfd1.o $(GLFDOBJS)

TGLFD2OBJS	= tglfd2.o $(GLFDOBJS)

TGLFD3OBJS	= tglfd3.o $(GLFDOBJS)

TGLFD4OBJS	= tglfd4.o $(GLFDOBJS)

TGLFD5OBJS	= tglfd5.o $(GLFDOBJS)

TGLFD6OBJS	= tglfd6.o $(GLFDOBJS)

TQGLF1OBJS	= tqglf1.o $(QGLFOBJS)

TQGLF2OBJS	= tqglf2.o $(QGLFOBJS)

TQGLF3OBJS	= tqglf3.o $(QGLFOBJS)

TQGLF4OBJS	= tqglf4.o $(QGLFOBJS)

TQGLF5OBJS	= tqglf5.o $(QGLFOBJS)

TQGLF6OBJS	= tqglf6.o $(QGLFOBJS)

TQGLFD1OBJS	= tqglfd1.o $(QGLFDOBJS)

TQGLFD2OBJS	= tqglfd2.o $(QGLFDOBJS)

TQGLFD3OBJS	= tqglfd3.o $(QGLFDOBJS)

TQGLFD4OBJS	= tqglfd4.o $(QGLFDOBJS)

TQGLFD5OBJS	= tqglfd5.o $(QGLFDOBJS)

TQGLFD6OBJS	= tqglfd6.o $(QGLFDOBJS)

#=======================================================================
# Suffixes and extra rule(s)

.SUFFIXES:	.s .i .c .f

# Preprocessor output:
.c.i:
	$(CC) $(CFLAGS) -E $< > $@

# Assembly code:
.c.s:
	$(CC) $(CFLAGS) -S $<

# Object code (an explicit "-o $*.o" makes compilations in ../common work properly!)
.f.o:
	$(FC) $(FFLAGS) -c -o $*.o $<

# Assembly code:
.f.s:
	$(FC) $(FFLAGS) -S $<

#=======================================================================
# Targets:

all:	all-dp

all-c:	all cglfd1

all-dp:	lib-dp $(TESTPROGS)

all-qp:	lib-qp $(QTESTPROGS)

cglfd1:	$(CGLFD1OBJS)
	$(CC) $(CFLAGS) -o $@ $(CGLFD1OBJS) $(CLIBS)

cglfd1.o:	cglfd1.c ../gjl.h

cqglfd1:	$(CQGLFD1OBJS)
	$(CC) $(CFLAGS) -o $@ $(CQGLFD1OBJS) $(CLIBS)

cqglfd1.o:	cqglfd1.c ../gjl.h
	$(CC) $(CFLAGS) -DHAVE_LONG_DOUBLE -c cqglfd1.c

check:	check-dp

check-dp:	all
	@if $(TEST) -f test/glf101.in ; \
	then \
		$(MAKE) -k -s check-header check-glf check-glfd ; \
	else \
		echo ======================================================================== ; \
		echo WARNING: This seems to be the slim package distribution without the ; \
		echo extensive validation suite: running the quick demos instead. ; \
		echo ======================================================================== ; \
		$(MAKE) demo ; \
	fi

check-qp:	all
	@if $(TEST) -f test/glf101.in ; \
	then \
		$(MAKE) -k -s check-header check-qglf check-qglfd ; \
	else \
		echo ======================================================================== ; \
		echo WARNING: This seems to be the slim package distribution without the ; \
		echo extensive validation suite: running the quick demos instead. ; \
		echo ======================================================================== ; \
		$(MAKE) demo-qp ; \
	fi

check-header:
	@echo ========================================================================
	@echo There should be no output from these tests, except their names.  If
	@echo ndiff is available, there will also be reports of the maximum absolute
	@echo and relative errors found.  If small numeric differences are reported,
	@echo then YOU must decide whether they are significant.  You can do
	@echo "        " make check RELERR=xxx
	@echo to rerun the tests with a different relative error tolerance, where xxx
	@echo is either a suitable small positive floating-point number, or else a
	@echo 'value >= 1.0, meaning that multiple of the machine epsilon.'
	@echo ========================================================================

check-c:	check-cglfd1

check-cglfd1:	cglfd1
	$(MAKE) -s check-xxx TESTNAME=cglfd1 CHECKFILES="$(CHECKCGLFD1)"

check-c-qp:	check-cqglfd1

check-cqglfd1:	cqglfd1
	$(MAKE) -s check-xxx TESTNAME=cqglfd1 CHECKFILES="$(CHECKCQGLFD1)"

check-glf:	check-glf1 check-glf2 check-glf3 check-glf4 check-glf5 check-glf6

check-glf1:	tglf1
	$(MAKE) -s check-xxx TESTNAME=tglf1 CHECKFILES="$(CHECKGLF1)"

check-glf2:	tglf2
	$(MAKE) -s check-xxx TESTNAME=tglf2 CHECKFILES="$(CHECKGLF2)"

check-glf3:	tglf3
	$(MAKE) -s check-xxx TESTNAME=tglf3 CHECKFILES="$(CHECKGLF3)"

check-glf4:	tglf4
	$(MAKE) -s check-xxx TESTNAME=tglf4 CHECKFILES="$(CHECKGLF4)"

check-glf5:	tglf5
	$(MAKE) -s check-xxx TESTNAME=tglf5 CHECKFILES="$(CHECKGLF5)"

check-glf6:	tglf6
	$(MAKE) -s check-xxx TESTNAME=tglf6 CHECKFILES="$(CHECKGLF6)"

check-glfd:	check-glfd1 check-glfd2 check-glfd3 check-glfd4 check-glfd5 check-glfd6

check-glfd1:	tglfd1
	$(MAKE) -s check-xxx TESTNAME=tglfd1 CHECKFILES="$(CHECKGLFD1)"

check-glfd2:	tglfd2
	$(MAKE) -s check-xxx TESTNAME=tglfd2 CHECKFILES="$(CHECKGLFD2)"

check-glfd3:	tglfd3
	$(MAKE) -s check-xxx TESTNAME=tglfd3 CHECKFILES="$(CHECKGLFD3)"

check-glfd4:	tglfd4
	$(MAKE) -s check-xxx TESTNAME=tglfd4 CHECKFILES="$(CHECKGLFD4)"

check-glfd5:	tglfd5
	$(MAKE) -s check-xxx TESTNAME=tglfd5 CHECKFILES="$(CHECKGLFD5)"

check-glfd6:	tglfd6
	$(MAKE) -s check-xxx TESTNAME=tglfd6 CHECKFILES="$(CHECKGLFD6)"

check-qglf:	check-qglf1 check-qglf2 check-qglf3 check-qglf4 check-qglf5 check-qglf6

check-qglf1:	tqglf1
	$(MAKE) -s check-xxx TESTNAME=tqglf1 CHECKFILES="$(CHECKQGLF1)"

check-qglf2:	tqglf2
	$(MAKE) -s check-xxx TESTNAME=tqglf2 CHECKFILES="$(CHECKQGLF2)"

check-qglf3:	tqglf3
	$(MAKE) -s check-xxx TESTNAME=tqglf3 CHECKFILES="$(CHECKQGLF3)"

check-qglf4:	tqglf4
	$(MAKE) -s check-xxx TESTNAME=tqglf4 CHECKFILES="$(CHECKQGLF4)"

check-qglf5:	tqglf5
	$(MAKE) -s check-xxx TESTNAME=tqglf5 CHECKFILES="$(CHECKQGLF5)"

check-qglf6:	tqglf6
	$(MAKE) -s check-xxx TESTNAME=tqglf6 CHECKFILES="$(CHECKQGLF6)"

check-qglfd:	check-qglfd1 check-qglfd2 check-qglfd3 check-qglfd4 check-qglfd5 check-qglfd6

check-qglfd1:	tqglfd1
	$(MAKE) -s check-xxx TESTNAME=tqglfd1 CHECKFILES="$(CHECKQGLFD1)"

check-qglfd2:	tqglfd2
	$(MAKE) -s check-xxx TESTNAME=tqglfd2 CHECKFILES="$(CHECKQGLFD2)"

check-qglfd3:	tqglfd3
	$(MAKE) -s check-xxx TESTNAME=tqglfd3 CHECKFILES="$(CHECKQGLFD3)"

check-qglfd4:	tqglfd4
	$(MAKE) -s check-xxx TESTNAME=tqglfd4 CHECKFILES="$(CHECKQGLFD4)"

check-qglfd5:	tqglfd5
	$(MAKE) -s check-xxx TESTNAME=tqglfd5 CHECKFILES="$(CHECKQGLFD5)"

check-qglfd6:	tqglfd6
	$(MAKE) -s check-xxx TESTNAME=tqglfd6 CHECKFILES="$(CHECKQGLFD6)"

# This is an internal target for use ONLY by other Makefile targets.
# Invoke with TESTNAME=... CHECKFILES="..."
check-xxx:
	@echo Checking $(TESTNAME) ... ; \
	NDIFF="$(NDIFF)" ; \
	if test -n "$$NDIFF" ; \
	then \
		DIFF="$$NDIFF $(NDIFFFLAGS)" ; \
	else \
		echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ; \
		echo You should consider installing ndiff: the difference output will likely ; \
		echo be voluminous with diff.  Further information is in the ../INSTALL file. ; \
		echo !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ; \
		DIFF=$(DIFF) ; \
	fi ; \
	echo File differences are made with: $$DIFF ; \
	cd test ; \
	for f in $(CHECKFILES) ; \
	do \
		echo $$f ; \
		$(RM) output.dat ; \
		../$(TESTNAME) < $$f.in > $$f.out 2> $$f.err ; \
		$(TEST) -s $$f.err || $(RM) $$f.err ; \
		$(MV) output.dat $$f.dat ; \
		if $(CMP) -s okay/$$f.out $$f.out ; \
		then \
			$(RM) $$f.out ; \
		else \
			$$DIFF okay/$$f.out $$f.out ; \
		fi ; \
		if $(CMP) -s okay/$$f.dat $$f.dat ; \
		then \
			$(RM) $$f.dat ; \
		else \
			$$DIFF okay/$$f.dat $$f.dat ; \
		fi ; \
	done

clean:	clean-test
	-$(RM) *.dcl
	-$(RM) *.i
	-$(RM) *.lst
	-$(RM) *.o
	-$(RM) *.s
	-$(RM) *.stb
	-$(RM) *.stx
	-$(RM) *.u
	-$(RM) *~
	-$(RM) \#*
	-$(RM) core
	-$(RM) core.*
	-$(RM) demglq.tab
	-$(RM) output.dat
	-$(RM) tg_[0-9][0-9][0-9][0-9][0-9].*
	-(cd $(COMMON) ; $(MAKE) $@)

clean-exe:
	-$(RM) a.out
	-$(RM) $(CTESTPROGS)
	-$(RM) $(DEMOPROGS)
	-$(RM) $(TESTPROGS)
	-$(RM) $(QCTESTPROGS)
	-$(RM) $(QDEMOPROGS)
	-$(RM) $(QTESTPROGS)

clean-test:
	-$(RM) test/*.dat
	-$(RM) test/*.err
	-$(RM) test/*.ls?
	-$(RM) test/*.nbs
	-$(RM) test/*.out
	-$(RM) test/core
	-$(RM) test/core.*

clobber:	distclean

# We need to reverse the search order for Fortran compilers, because g77
# on IBM and HP systems adds trailing underscores to names, while the
# native Fortran compilers do not.
configure:	configure.in FRC
	$(AUTOCONF)
	$(SED) -e 's/for ac_prog in g77 f77 f2c/for ac_prog in f77 g77 f90 f95/' < $@ >$@.tmp
	$(CHMOD) 775 $@.tmp
	$(CMP) -s $@.tmp $@ || $(MV) $@.tmp $@
	-$(RM) $@.tmp

demglq:	demglq.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ demglq.o $(LDFLAGS) $(LIBS)

demglq.tab:	demglq FRC
	echo -0.9375 20 40 | \
		./demglq | \
			$(AWK) '(NF == 5) && ($$1 ~ /^[0-9]+$$/)' | \
				$(AWK) '{printf("%14d  &  %9.2e  &  %9.2e  &  %9.2e  \\\\\n", $$1, $$2, $$4, $$3)}' | \
					$(SED) -e 's/e[+]/e$$+$$/g' -e 's/e-/e$$-$$/g' -e 's/ -/  /g' > $@

demo:	demo-dp

demo-dp:	demo-info demo-dp-1 demo-dp-2

demo-dp-1:	all-dp demo-setup demglq
	@echo
	@echo "This output file should closely match results in Table 1 of the paper"
	@echo "Efficient High-Accuracy Quadrature for Two Classes of Logarithmic"
	@echo "Weight Functions"
	@echo
	echo -0.9375 20 40 | ./demglq > test/demglq.out
	@$(LS) -l test/demglq.out
	@echo
	@echo "Correct output from the development site is available here:"
	@echo
	@$(LS) -l test/okay/demglq.out
	@echo

demo-dp-2:	all-dp demo-setup demglq
	@echo
	@echo "This output file must reproduce the Gauss-Laguerre n=15 section"
	@echo "of Table 25.9 of Abramowitz and Stegen, Handbook of Mathematical"
	@echo "Functions, NBS #55, p. 923"
	@echo
	echo 0 15 0 | ./demglq > test/demglq.nbs
	@$(LS) -l test/demglq.nbs
	@echo
	@echo "Correct output from the development site is available here:"
	@echo
	@$(LS) -l test/okay/demglq.nbs
	@echo

demo-info:
	@echo "========================================================================"
	@echo "The demonstration programs demglq and demqglq read a single free-format"
	@echo "line of data on standard input with the three numeric values"
	@echo ""
	@echo "        alpha nquad pmax"
	@echo ""
	@echo "where"
	@echo ""
	@echo "        alpha = power of x in Laguerre weight function"
	@echo "        nquad = quadrature order"
	@echo "        pmax  = maximum power of x in test function"
	@echo ""
	@echo "In the demonstration targets, these values are supplied by a simple"
	@echo "echo command embedded in the Makefile.  Results are written to standard"
	@echo "output."
	@echo ========================================================================

demo-lapack:	lapack-build lapack-demo

demo-qp:	demo-info demo-qp-1 demo-qp-2

demo-qp-1:	all-qp demo-setup demqglq
	@echo
	@echo "This output file should produce quadruple-precision results analogous to"
	@echo "double-precision results in Table 1 of the paper"
	@echo "Fast Gaussian Quadrature for Two Classes of Logarithmic Weight Functions"
	@echo
	echo -0.9375 20 40 | ./demqglq > test/demqglq.out
	@$(LS) -l test/demqglq.out
	@echo
	@echo "Correct output from the development site is available here:"
	@echo
	@$(LS) -l test/okay/demqglq.out
	@echo

demo-qp-2:	all-qp demo-setup demqglq
	@echo
	@echo "This output file must reproduce the Gauss-Laguerre n=15 section"
	@echo "of Table 25.9 of Abramowitz and Stegen, Handbook of Mathematical"
	@echo "Functions, NBS #55, p. 923"
	@echo
	echo 0 15 0 | ./demqglq > test/demqglq.nbs
	@$(LS) -l test/demqglq.nbs
	@echo
	@echo "Correct output from the development site is available here:"
	@echo
	@$(LS) -l test/okay/demqglq.nbs
	@echo

demo-setup:
	@-$(TEST) -d test || $(MKDIR) test
	@-$(TEST) -d test/okay || $(MKDIR) test/okay

demqglq:	demqglq.o $(LIBRARY)
	$(FC) $(FFLAGS) -o $@ demqglq.o $(LDFLAGS) $(LIBS)

distclean:	mostlyclean clean-test
	-$(RM) config.cache config.status config.log Makefile

FRC:

install:
	@echo There is nothing to install from this directory

install-c:	install

## Executables with tql1 wrappers that call LAPACK eigenvalue solvers
lapack-build:	all-dp
	for f in eqr erf ev evd evr evx ; \
	do \
		echo demglq-$$f ; \
		$(FC) $(FFLAGS) -o demglq-$$f demglq.f ../common/tql1-$$f.f $(LDFLAGS) $(LIBS) $(LAPACK) ; \
	done

## Tests with wrapped LAPACK eigenvalue solvers
lapack-demo:	all-dp demglq
	$(CP) demglq demglq-tql1
	for f in eqr erf ev evd evr evx tql1 ; \
	do \
		echo demglq-$$f ; \
		echo -0.9375 60 120 | ./demglq-$$f > test/gjq-$$f.out 2>/dev/null ; \
		echo  0 15 30       | ./demglq-$$f > test/gjq-$$f.nbs 2>/dev/null ; \
		echo -0.5 60 120    | ./demglq-$$f > test/gjq-$$f.ls1 2>/dev/null ; \
		echo  0.0 60 120    | ./demglq-$$f > test/gjq-$$f.ls2 2>/dev/null ; \
		echo  0.5 60 120    | ./demglq-$$f > test/gjq-$$f.ls3 2>/dev/null ; \
		echo  5.0 60 120    | ./demglq-$$f > test/gjq-$$f.ls4 2>/dev/null ; \
		echo 10.0 60 120    | ./demglq-$$f > test/gjq-$$f.ls5 2>/dev/null ; \
	done
	$(LS) -l test/glq-*.*

lib-dp:	$(LIBOBJS)
	(cd ../common; $(MAKE) $@)
	$(AR) -r $(LIBRARY) $(LIBOBJS)
	$(RANLIB) $(LIBRARY) || true

lib-qp:	$(QLIBOBJS)
	(cd ../common; $(MAKE) $@)
	$(AR) -r $(LIBRARY) $(QLIBOBJS)
	$(RANLIB) $(LIBRARY) || true

maintainer-clean:	distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

mostlyclean:	clean clean-exe
	-$(RM) *.d
	-$(RM) *.tcov

tglf1:	$(TGLF1OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF1OBJS)

tglf2:	$(TGLF2OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF2OBJS)

tglf3: $(TGLF3OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF3OBJS)

tglf4: $(TGLF4OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF4OBJS)

tglf5: $(TGLF5OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF5OBJS)

tglf6: $(TGLF6OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLF6OBJS)

tglfd1:	$(TGLFD1OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD1OBJS)

tglfd2:	$(TGLFD2OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD2OBJS)

tglfd3:	$(TGLFD3OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD3OBJS)

tglfd4:	$(TGLFD4OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD4OBJS)

tglfd5:	$(TGLFD5OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD5OBJS)

tglfd6:	$(TGLFD6OBJS)
	$(FC) $(FFLAGS) -o $@ $(TGLFD6OBJS)

tqglf1:	$(TQGLF1OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF1OBJS)

tqglf2:	$(TQGLF2OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF2OBJS)

tqglf3: $(TQGLF3OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF3OBJS)

tqglf4: $(TQGLF4OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF4OBJS)

tqglf5: $(TQGLF5OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF5OBJS)

tqglf6: $(TQGLF6OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLF6OBJS)

tqglfd1:	$(TQGLFD1OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD1OBJS)

tqglfd2:	$(TQGLFD2OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD2OBJS)

tqglfd3:	$(TQGLFD3OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD3OBJS)

tqglfd4:	$(TQGLFD4OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD4OBJS)

tqglfd5:	$(TQGLFD5OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD5OBJS)

tqglfd6:	$(TQGLFD6OBJS)
	$(FC) $(FFLAGS) -o $@ $(TQGLFD6OBJS)

uninstall:

#=======================================================================
# File dependencies
# Fortran/SFTRAN3 header file dependencies

demglq.o:	$(COMMON)/maxpts.inc

demqglq.o:	$(COMMON)/maxpts.inc

glqf.o:		$(COMMON)/ecodes.inc
glqf.o:		$(COMMON)/maxpts.inc

glqfd.o:	$(COMMON)/ecodes.inc
glqfd.o:	$(COMMON)/maxpts.inc

glqli.o:	$(COMMON)/maxpts.inc

glqrc.o:	$(COMMON)/ecodes.inc
glqrc.o:	$(COMMON)/maxpts.inc

qglqf.o:	$(COMMON)/ecodes.inc
qglqf.o:	$(COMMON)/maxpts.inc

qglqfd.o:	$(COMMON)/ecodes.inc
qglqfd.o:	$(COMMON)/maxpts.inc

qglqli.o:	$(COMMON)/maxpts.inc

qglqrc.o:	$(COMMON)/ecodes.inc
qglqrc.o:	$(COMMON)/maxpts.inc

tglf1.o:	$(COMMON)/maxpts.inc
tglf1.o:	$(COMMON)/stdio.inc

tglf2.o:	$(COMMON)/maxpts.inc
tglf2.o:	$(COMMON)/stdio.inc

tglf3.o:	$(COMMON)/maxpts.inc
tglf3.o:	$(COMMON)/stdio.inc

tglf4.o:	$(COMMON)/maxpts.inc
tglf4.o:	$(COMMON)/stdio.inc

tglf5.o:	$(COMMON)/maxpts.inc
tglf5.o:	$(COMMON)/stdio.inc

tglf6.o:	$(COMMON)/maxpts.inc
tglf6.o:	$(COMMON)/stdio.inc

tglfd1.o:	$(COMMON)/maxpts.inc
tglfd1.o:	$(COMMON)/stdio.inc

tglfd2.o:	$(COMMON)/maxpts.inc
tglfd2.o:	$(COMMON)/stdio.inc

tglfd3.o:	$(COMMON)/maxpts.inc
tglfd3.o:	$(COMMON)/stdio.inc

tglfd4.o:	$(COMMON)/maxpts.inc
tglfd4.o:	$(COMMON)/stdio.inc

tglfd5.o:	$(COMMON)/maxpts.inc
tglfd5.o:	$(COMMON)/stdio.inc

tglfd6.o:	$(COMMON)/maxpts.inc
tglfd6.o:	$(COMMON)/stdio.inc

tqglf1.o:	$(COMMON)/maxpts.inc
tqglf1.o:	$(COMMON)/stdio.inc

tqglf2.o:	$(COMMON)/maxpts.inc
tqglf2.o:	$(COMMON)/stdio.inc

tqglf3.o:	$(COMMON)/maxpts.inc
tqglf3.o:	$(COMMON)/stdio.inc

tqglf4.o:	$(COMMON)/maxpts.inc
tqglf4.o:	$(COMMON)/stdio.inc

tqglf5.o:	$(COMMON)/maxpts.inc
tqglf5.o:	$(COMMON)/stdio.inc

tqglf6.o:	$(COMMON)/maxpts.inc
tqglf6.o:	$(COMMON)/stdio.inc

tqglfd1.o:	$(COMMON)/maxpts.inc
tqglfd1.o:	$(COMMON)/stdio.inc

tqglfd2.o:	$(COMMON)/maxpts.inc
tqglfd2.o:	$(COMMON)/stdio.inc

tqglfd3.o:	$(COMMON)/maxpts.inc
tqglfd3.o:	$(COMMON)/stdio.inc

tqglfd4.o:	$(COMMON)/maxpts.inc
tqglfd4.o:	$(COMMON)/stdio.inc

tqglfd5.o:	$(COMMON)/maxpts.inc
tqglfd5.o:	$(COMMON)/stdio.inc

tqglfd6.o:	$(COMMON)/maxpts.inc
tqglfd6.o:	$(COMMON)/stdio.inc

#=======================================================================
