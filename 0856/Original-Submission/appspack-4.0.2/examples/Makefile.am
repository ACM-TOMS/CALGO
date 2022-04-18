## $Id: Makefile.am,v 1.5.2.3 2005/06/29 17:26:06 tgkolda Exp $ 
## $Source: /space/CVS-Acro/acro/packages/appspack/appspack/examples/Makefile.am,v $ 
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

if BUILD_EXECUTABLES

noinst_PROGRAMS = example1 example2 example3

endif

example1_SOURCES = $(srcdir)/example1.c 

example2_SOURCES = $(srcdir)/example2.cpp

example3_SOURCES = $(srcdir)/example3.cpp

AM_CPPFLAGS = -I$(top_srcdir)/src 

## #######################################################################
## Extra files to be distributed
## #######################################################################

EXTRA_DIST = $(srcdir)/README \
	$(srcdir)/example1.apps $(srcdir)/example2.apps $(srcdir)/example3.apps \
	$(srcdir)/example1.out $(srcdir)/example2.out $(srcdir)/example3.out

## #######################################################################
## Files to be deleted by 'make maintainer-clean'
## #######################################################################

MAINTAINERCLEANFILES = Makefile.in 

