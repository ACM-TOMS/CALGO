## $Id: Makefile.am,v 1.4.2.2 2005/06/30 01:00:32 tgkolda Exp $ 
## $Source: /space/CVS-Acro/acro/packages/appspack/appspack/Makefile.am,v $ 
# @HEADER
# ************************************************************************
# 
#          APPSPACK: Asynchronous Parallel Pattern Search
#                 Copyright (2003) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#                                                                                 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA.                                                                           .
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) 
# 
# ************************************************************************
# @HEADER

## #######################################################################
## Options to automake (rarely used - don't worry about it)
## #######################################################################
AUTOMAKE_OPTIONS = foreign

## #######################################################################
## Aclocal command (rarely used - don't worry about it)
## #######################################################################
ACLOCAL	= aclocal -I config

## #######################################################################
## Additional files to be included in distribution for 'make dist'
## #######################################################################
EXTRA_DIST = \
	$(srcdir)/LICENSE \
	$(srcdir)/README \
	$(srcdir)/example-custom-evaluator/Makefile \
	$(srcdir)/example-custom-evaluator/README \
	$(srcdir)/example-custom-evaluator/evaluator.cpp \
	$(srcdir)/example-custom-evaluator/evaluator.hpp \
	$(srcdir)/example-custom-evaluator/main.cpp \
	$(srcdir)/example-custom-evaluator/custom.out \
	$(srcdir)/example-custom-executor/Makefile \
	$(srcdir)/example-custom-executor/README \
	$(srcdir)/example-custom-executor/executor.cpp \
	$(srcdir)/example-custom-executor/executor.hpp \
	$(srcdir)/example-custom-executor/main.cpp \
	$(srcdir)/example-custom-executor/master.cpp \
	$(srcdir)/example-custom-executor/msgtags.hpp \
	$(srcdir)/example-custom-executor/worker.cpp \
	$(srcdir)/example-custom-executor/custom.out

## #######################################################################
## Tools in the auxillary directory 
## #######################################################################
AUX_DIST = config/install-sh config/missing config/mkinstalldirs

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 autom4te.cache/* \
	configure config.status config.log \
	src/common/config-h.in src/common/stamp-h.in \
	$(AUX_DIST) 

## #######################################################################
## Subdirectories to be make'd recursively
## #######################################################################
SUBDIRS = src examples doc



