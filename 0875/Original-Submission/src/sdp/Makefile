
DSDPOBJ	=  dsdpstep.o printsdpa.o  sdpconevec.o \
	   sdpsss.o dsdpadddata.o dsdpadddatamat.o dsdpblock.o \
	   sdpcone.o sdpkcone.o sdpcompute.o sdpconesetup.o \
	   dsdpdsmat.o dsdpdatamat.o dsdpdualmat.o dsdpxmat.o

DCFLAGS   = -I${PP}solver

include ../../make.include

default:
	make dsdplib;

clean:
	${RM} ${DSDPOBJ}


