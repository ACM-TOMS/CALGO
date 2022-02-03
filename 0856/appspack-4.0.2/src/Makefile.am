## $Id: Makefile.am,v 1.8 2004/11/23 20:40:33 tgkolda Exp $ 
## $Source: /space/CVS-Acro/acro/packages/appspack/appspack/src/Makefile.am,v $ 
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

# ------------------------------------------------------------------------
# List of all programs
# ------------------------------------------------------------------------

if BUILD_SERIAL

APPSPACK_SERIAL = appspack_serial

endif

if BUILD_MPI

APPSPACK_MPI = appspack_mpi

endif

if BUILD_EXECUTABLES

bin_PROGRAMS = $(APPSPACK_SERIAL) $(APPSPACK_MPI)

endif

# ------------------------------------------------------------------------
# Executables
# ------------------------------------------------------------------------

appspack_serial_SOURCES = $(srcdir)/APPSPACK_Main_Serial.cpp

appspack_serial_DEPENDENCIES = $(top_builddir)/src/libappspack.a 

appspack_serial_LDADD = $(appspack_serial_DEPENDENCIES)

appspack_mpi_SOURCES = $(srcdir)/APPSPACK_Main_MPI.cpp

appspack_mpi_DEPENDENCIES = $(top_builddir)/src/libappspack.a 

appspack_mpi_LDADD = $(appspack_mpi_DEPENDENCIES)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = libappspack.a

# ------------------------------------------------------------------------
# Library AR Command
# ------------------------------------------------------------------------

if USE_ALTERNATE_AR

libappspack_a_AR = $(ALTERNATE_AR)

else

libappspack_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# General CPP flags
# ------------------------------------------------------------------------

AM_CPPFLAGS = 

# ------------------------------------------------------------------------
# Library source
# ------------------------------------------------------------------------

libappspack_a_SOURCES = \
	$(srcdir)/APPSPACK_Cache_Manager.cpp \
	$(srcdir)/APPSPACK_Cache_Point.cpp \
	$(srcdir)/APPSPACK_Constraints_Bounds.cpp \
	$(srcdir)/APPSPACK_Conveyor.cpp \
	$(srcdir)/APPSPACK_Counter.cpp \
	$(srcdir)/APPSPACK_Directions.cpp \
	$(srcdir)/APPSPACK_Evaluator_SystemCall.cpp \
	$(srcdir)/APPSPACK_Executor_Serial.cpp \
	$(srcdir)/APPSPACK_Executor_MPI.cpp \
	$(srcdir)/APPSPACK_GCI.cpp \
	$(srcdir)/APPSPACK_List.cpp \
	$(srcdir)/APPSPACK_Parameter_Entry.cpp \
	$(srcdir)/APPSPACK_Parameter_List.cpp \
	$(srcdir)/APPSPACK_Point.cpp \
	$(srcdir)/APPSPACK_Print.cpp \
	$(srcdir)/APPSPACK_Solver.cpp \
	$(srcdir)/APPSPACK_Utils.cpp \
	$(srcdir)/APPSPACK_Value.cpp \
	$(srcdir)/APPSPACK_Vector.cpp

# ------------------------------------------------------------------------
# Header files to install
# ------------------------------------------------------------------------

pkginclude_HEADERS = \
	$(srcdir)/APPSPACK_Cache_Manager.hpp \
	$(srcdir)/APPSPACK_Cache_Point.hpp \
	$(srcdir)/APPSPACK_Cache_SplayTree.hpp \
	$(srcdir)/APPSPACK_Common.hpp \
	$(srcdir)/APPSPACK_Constraints_Bounds.hpp \
	$(srcdir)/APPSPACK_Constraints_Interface.hpp \
	$(srcdir)/APPSPACK_Conveyor.hpp \
	$(srcdir)/APPSPACK_Counter.hpp \
	$(srcdir)/APPSPACK_Directions.hpp \
	$(srcdir)/APPSPACK_Evaluator_Interface.hpp \
	$(srcdir)/APPSPACK_Evaluator_SystemCall.hpp \
	$(srcdir)/APPSPACK_Executor_Interface.hpp \
	$(srcdir)/APPSPACK_Executor_Serial.hpp \
	$(srcdir)/APPSPACK_Executor_MPI.hpp \
	$(srcdir)/APPSPACK_GCI.hpp \
	$(srcdir)/APPSPACK_List.hpp \
	$(srcdir)/APPSPACK_Parameter_Entry.hpp \
	$(srcdir)/APPSPACK_Parameter_List.hpp \
	$(srcdir)/APPSPACK_Point.hpp \
	$(srcdir)/APPSPACK_Print.hpp \
	$(srcdir)/APPSPACK_Solver.hpp \
	$(srcdir)/APPSPACK_Utils.hpp \
	$(srcdir)/APPSPACK_Value.hpp \
	$(srcdir)/APPSPACK_Vector.hpp 

# ------------------------------------------------------------------------
# Extra files to distribute
# ------------------------------------------------------------------------

EXTRA_DIST = 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

