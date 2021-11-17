##############################################################################
##############################################################################
##               Makefile for test program TESTSCASY                        ##
##                                                                          ##
##############################################################################
##############################################################################
#
include ./../Make_include
#
TARGET = testscasy
#
SOURCE = TESTSCASY.f 
#
OBJECT = TESTSCASY.o
#
#
LIBS = $(RECSY) $(SLICOT) $(SCALAPACK) $(MPI_BLACS) $(BLAS) $(LAPACK) \
       $(BLAS) $(BLAS) 
#

$(TARGET):	$(OBJECT)
		$(LINKER) $(LFLAGS) -o $(TARGET) $(OBJECT) \
		-L$(SCASY)/lib -lscasy $(LIBS)
#

$(OBJECT):	$(SOURCE)
		$(FC) $(FFLAGS) $(APPFLAGS) $(INCLUDES) -c $(SOURCE)

#
clean:
		$(RM) $(TARGET) $(OBJECT)
#
##############################################################################
##############################################################################
