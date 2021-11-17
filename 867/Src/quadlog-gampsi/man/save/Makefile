#=======================================================================
# Makefile for files in /u/sy/beebe/jim-ball/toms/dist
#
# This directory contains no source files itself, but this Makefile
# directs the compilation, linking, testing, et al targets in each of
# its subdirectories.
#
# Current target list:
#	all			make everything
#	check			run validation suite
#	clean			remove some intermediate files
#	clobber			same as distclean
#	configure		patch a newly-generated configure file
#	distclean		remove almost everything that make can
#				rebuild, returning the directory to the
#				state of the original distribution
#	dw			doubled-word check
#	html			same as htmlfiles
#	htmlfiles		make HTML files
#	install			install software on local system
#	install-html		install HTML files only
#	install-lib		install libraries only
#	install-man		install manual pages only
#	install-pdf		install PDF files only
#	install-ps		install PostScript files only
#	install-show		show what has been installed
#	maintainer-clean	remove absolutely everything that make
#				can rebuilds; do NOT use this: it may
#				need tools that are available only at
#				the authors' site.
#	mostlyclean		remove all intermediate files
#	pdf			same as pdffiles
#	pdffiles		make PDF files
#	ps			same as psfiles
#	psfiles			make PostScript files
#	spell			spelling check
#	uninstall		remove installed files from local system
#	uninstall-html		remove only installed HTML files
#	uninstall-lib		remove only installed libraries
#	uninstall-man		remove only installed manual pages
#	uninstall-pdf		remove only installed PDF files
#	uninstall-ps		remove only installed PostScript files
#
# [13-May-2000]
#=======================================================================

prefix		= /usr/local

htmldir		= $(prefix)/man/html/man$(manext)

mandir		= $(prefix)/man/man$(manext)

manext		= 3

pdfdir		= $(prefix)/man/pdf/man$(manext)

psdir		= $(prefix)/man/ps/man$(manext)

#=======================================================================
# Programs:

AUTOCONF	= autoconf

CMP		= cmp

CHMOD		= chmod

DEROFF		= deroff

DW		= dw

DISTILL		= ps2pdf

LS		= ls

# man2html is available at ftp://ftp.math.utah.edu/pub/sgml/man2html-x.yy.*
MAN2HTML	= man2html

MAN2PS		= ./man2ps

MV		= mv

RM		= rm -f

SED		= sed

SHELL		= /bin/sh

SORT		= sort

SPELL		= spell

#=======================================================================
# Compiler flags and source files

DWFILES		= $(MANFILES:.man=.dw)

HTMLFILES	= $(MANFILES:.man=.html)

MANFILES        = aeps.man afpmax.man ainf.man airan.man algam.man alog2.man \
		  anan.man astore.man astorf.man avsum.man dcopy.man \
		  deps.man depsln.man derbit.man dfloat.man dfpmax.man \
		  dgamma.man dinf.man diran.man dlgam.man dlog2.man dnan.man \
		  dpsi.man dpsiln.man dran.man dstore.man dstorf.man \
		  dvsum.man ffalse.man ftrue.man gamma.man gjqf.man \
		  gjqfd.man gjqrc.man glqf.man glqfd.man glqrc.man iceil.man \
		  idceil.man iqceil.man isainf.man isanan.man isdinf.man \
		  isdnan.man isqinf.man isqnan.man psi.man psiln.man \
		  pythag.man qeps.man qepsln.man qfloat.man qfpmax.man \
		  qgamma.man qgjqf.man qgjqfd.man qgjqrc.man qglqf.man \
		  qglqfd.man qglqrc.man qinf.man qiran.man qlgam.man \
		  qlog2.man qnan.man qpsi.man qpsiln.man qstore.man \
		  qstorf.man qvsum.man ran.man tql1.man

PDFFILES	= $(MANFILES:.man=.pdf)

PSFILES		= $(MANFILES:.man=.ps)

SERFILES	= $(MANFILES:.man=.ser)

#=======================================================================
# Suffixes and extra rule(s)

.SUFFIXES:	.dw .ser .i .pdf .ps .html .man

.ps.pdf:
	$(DISTILL) $<

.man.dw:
	$(DEROFF) $< | $(DW) > $@

.man.html:
	$(MAN2HTML) $<

.man.ps:
	$(MAN2PS) $< > $@

.man.ser:
	$(DEROFF) $< | $(SPELL) +$*.sok | $(SORT) -u > $@

#=======================================================================
# Targets:

all:	htmlfiles psfiles pdffiles

all-qp:	all

check:	check-dp

check-dp:
	@echo There is no validation suite in this directory

check-qp:
	@echo There is no validation suite in this directory

clean:
	-$(RM) *.dw
	-$(RM) *.ser
	-$(RM) *~
	-$(RM) \#*
	-$(RM) a.out
	-$(RM) core
	-$(RM) tg_[0-9][0-9][0-9][0-9][0-9].*

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

demo:	demo-dp

demo-dp:
	@echo There are no demonstration programs in this directory

demo-qp:
	@echo There are no demonstration programs in this directory

dw:		$(DWFILES)

distclean:	mostlyclean
	-$(RM) config.cache config.status config.log Makefile

FRC:

html htmlfiles:	$(HTMLFILES)

install:	install-html install-man install-pdf install-ps
	-$(MAKE) install-show

# This target is optional, because there is no standard place to
# install HTML versions of UNIX manual pages.  There is intentionally
# no subdirectory loop for it.  We use it only if the html directory
# already exists.
install-html:	uninstall-html
	if test -d $(htmldir) ; \
	then \
		for f in $(HTMLFILES) ; \
		do \
			$(SED) -e 's@/usr/local/@$(prefix)/@g' $$f > $(htmldir)/$$f ; \
			$(CHMOD) 664 $(htmldir)/$$f ; \
		done ; \
	fi

install-man:	uninstall-man
	for f in $(MANFILES) ; \
	do \
		g=`basename $$f .man`.$(manext) ; \
		$(SED) -e 's@/usr/local/@$(prefix)/@g' $$f > $(mandir)/$$g ; \
		$(CHMOD) 664 $(mandir)/$$g ; \
	done

# This target is optional, because there is no standard place to
# install PDF versions of UNIX manual pages.  There is intentionally
# no subdirectory loop for it.  We use it only if the html directory
# already exists.
install-pdf:	uninstall-pdf
	if test -d $(pdfdir) ; \
	then \
		for f in $(PDFFILES) ; \
		do \
			$(SED) -e 's@/usr/local/@$(prefix)/@g' $$f > $(pdfdir)/$$f ; \
			$(CHMOD) 664 $(pdfdir)/$$f ; \
		done ; \
	fi

# This target is optional, because there is no standard place to
# install PDF versions of UNIX manual pages.  There is intentionally
# no subdirectory loop for it.  We use it only if the html directory
# already exists.
install-ps:	uninstall-ps
	if test -d $(psdir) ; \
	then \
		for f in $(PSFILES) ; \
		do \
			$(SED) -e 's@/usr/local/@$(prefix)/@g' $$f > $(psdir)/$$f ; \
			$(CHMOD) 664 $(psdir)/$$f ; \
		done ; \
	fi

install-show:
	@echo 'Installed file(s)...'
	@for f in $(MANFILES) ; \
	do \
		g=`basename $$f .man`.$(manext) ; \
		$(LS) -l $(mandir)/$$g ; \
	done
	@if test -d $(htmldir) ; \
	then \
		for f in $(HTMLFILES) ; \
		do \
			$(LS) -l $(htmldir)/$$f ; \
		done ; \
	fi
	@if test -d $(pdfdir) ; \
	then \
		for f in $(PDFFILES) ; \
		do \
			$(LS) -l $(pdfdir)/$$f ; \
		done ; \
	fi
	@if test -d $(psdir) ; \
	then \
		for f in $(PSFILES) ; \
		do \
			$(LS) -l $(psdir)/$$f ; \
		done ; \
	fi

maintainer-clean:	distclean
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."
	-$(RM) $(HTMLFILES)
	-$(RM) $(PSFILES)
	-$(RM) $(PDFFILES)

mostlyclean:	clean

pdf pdffiles:	$(PDFFILES)

ps psfiles:	$(PSFILES)

spell:		$(SERFILES)

uninstall: uninstall-html uninstall-man uninstall-pdf uninstall-ps

uninstall-html:
	if test -d $(htmldir) ; \
	then \
		for f in $(HTMLFILES) ; \
		do \
			$(RM) -f $(htmldir)/$$f ; \
		done ; \
	fi

uninstall-man:
	for f in $(MANFILES) ; \
	do \
		g=`basename $$f .man`.$(manext) ; \
		$(RM) $(mandir)/$$g ; \
	done

uninstall-pdf:
	if test -d $(pdfdir) ; \
	then \
		for f in $(PDFFILES) ; \
		do \
			$(RM) -f $(pdfdir)/$$f ; \
		done ; \
	fi

uninstall-ps:
	if test -d $(psdir) ; \
	then \
		for f in $(PSFILES) ; \
		do \
			$(RM) -f $(psdir)/$$f ; \
		done ; \
	fi

#=======================================================================
