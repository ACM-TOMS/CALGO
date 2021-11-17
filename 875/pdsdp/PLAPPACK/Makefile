
PLAPACK_ROOT	= /home/benson/soft/PLAPACKR30

default: all

include ../../make.include
include ${PLAPACK_ROOT}/Make.include

LINKER = mpicc

DCFLAGS =  -I${DSDPROOT}/include -I${DSDPROOT}/

OBJS	= readsdpa.o pdsdpplapack.o

LIBLAPACK = /usr/lib/liblapack.a -lblas -lg2c -lm

dsdpsdpa: ${OBJS}
	${LINKER} -o pdsdp5 ${OBJS} ${DSDPLIB} ${PLAPACKLIB} ${LIBLAPACK}
	${RM} ${OBJS}

info2:
	-@echo "======================================================"
	-@echo "Compiling PDSDP using PLAPACK parallel linear solver."
	-@echo "  Check ../../make.include to set compiler."
	-@echo "  the mpicc compiler may produce faster DSDP libraries"
	-@echo "======================================================"

clean:
	${RM} ${OBJS}

all: info info2 dsdpsdpa



