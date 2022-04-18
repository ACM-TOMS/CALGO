####################################################################
#  SHHEIG main makefile                                            #
#  Top Level Makefile for generating SHHEIG Library object file,   #
#  the auxiliary library file, and linking and running the example #
#  programs.                                                       #
#  SHHEIG, Release 1.0                          ./shheig/makefile  #
#  July 23, 2013                                                   #
####################################################################
#
#  This makefile creates/updates the SHHEIG Library object file, the
#  auxiliary library, and compiles, links, and runs the example
#  programs for the SHHEIG Library. To perform all these actions, 
#  enter
#       make
#
#  To create/update the library for SHHEIG, enter 
#       make lib
#
#  To compile, link, and run the example programs, enter 
#       make example
#
#  To remove the object files for SHHEIG routines, enter
#       make cleanlib
#
#  To remove the files with the computed results (*.exa) as well as
#  the object files and executable files of the example progrmas,
#  enter
#       make cleanexample
#
#  To remove the object files for SHHEIG routines as well as the
#  files with the computed results (*.exa), enter
#       make clean
#
####################################################################

include make.inc

all: lib example

clean: cleanlib cleanexample

lib:
	( cd src; $(MAKE) )

example:
	( cd examples; $(MAKE) )

cleanlib:
	( cd src; $(MAKE) clean )

cleanexample:
	( cd examples; $(MAKE) clean )
