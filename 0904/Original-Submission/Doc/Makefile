##############################################################################
##############################################################################
##                                                                          ##
##                            SCASY Makefile                                ## 
##                                                                          ## 
## Version 1.0                                                              ##
##                                                                          ##
## Contributors: Robert Granat and Bo Kågström                              ##
##               Department of Computing Science and HPC2N                  ##
##               Umeå University, Sweden                                    ##
##                                                                          ##
## Email: {granat,bokg}@cs.umu.se                                           ##
##                                                                          ##
## Date: August 31, 2009.                                                   ##
##                                                                          ##
##                                                                          ##
##############################################################################
##############################################################################
#
include Make_include
#
all:		libscasy testscasy
#
libscasy:
		( cd ./lib/ && make all )
#
testscasy:
		( cd ./lib/ && make all )
		( cd ./test/ && make )
#
clean:
		( cd ./lib/  && make clean )
		( cd ./test/ && make clean )
#
##############################################################################
##############################################################################

