
DSDPOBJ  = vech.o vechu.o drowcol.o dlpack.o dufull.o \
	   sdpvec.o identity.o spds.o \
	   zeromat.o onemat.o diag.o rmmat.o \
	   cholmat.o cholmat2.o dtrsm2.o \
	   sdpmatx.o sdpnfac.o sdporder.o sdpalloc.o\
	   sdpsymb.o sdpxlist.o sdpdvct.o sdpexit.o

DCFLAGS   = -I${PP}solver -I${PP}sdp

include ../../make.include

default:
	make dsdplib;

clean:
	${RM} ${DSDPOBJ}

