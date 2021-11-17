
default:
	make all

include ../make.include

EXAMPLE_LIBS = ${DSDPLIB} ${LAPACKBLAS}

dsdp: exinfo clean
	${MEX} -I${DSDPROOT}/include dsdp.c ${DSDPLIB}; 
	${MV} dsdp.mex* ${DSDPMATLABDIR}

dsdpsdpa: exinfo clean readsdpa.o 
	${CLINKER} readsdpa.o ${EXAMPLE_LIBS} -o dsdp5
	${MV} dsdp5 ${EXECDIR}
	${RM} readsdpa.o

maxcut: exinfo clean maxcut.o
	${CLINKER} -o maxcut maxcut.o ${EXAMPLE_LIBS}
	${MV} maxcut ${EXECDIR}
	${RM} maxcut.o

stable: exinfo clean stable.o
	${CLINKER} -o stable stable.o ${EXAMPLE_LIBS}
	${MV} stable ${EXECDIR}
	${RM} stable.o

color: exinfo clean color.o
	${CLINKER} -o color color.o ${EXAMPLE_LIBS}
	${MV} color ${EXECDIR}
	${RM} color.o

theta:	exinfo clean theta.o
	${CLINKER} -o theta theta.o ${EXAMPLE_LIBS}
	${MV} theta ${EXECDIR}
	${RM} theta.o

copysrc:
	${CP} ${DSDPROOT}/src/*/*.c allc
	${CP} ${DSDPROOT}/src/*/*.h allc
	${CP} ${DSDPROOT}/include/*.h allc
	${CP} ${DSDPROOT}/include/*.h .

# On most Linux/Unix machines
dsdpagain: copysrc
	${MEX} -DDSDPMATLAB -DDSDP_TIME dsdp.c allc/*.c
	${MV} dsdp.mex* ${DSDPMATLABDIR}
	make clean;

# On a Windows Machine
dsdpagain2: copysrc
	${MEX} -DDSDPMATLAB -DDSDP_MS_TIME -DNOUNDERBLAS dsdp.c allc/*.c
	${MV} dsdp.mex* ${DSDPMATLABDIR}
	make clean;

sdpaagain: copysrc	
	${CC} -O3 -DDSDP_TIME readsdpa.c allc/*.c ${LAPACKBLAS} -o dsdp5
	${MV} dsdp5 ${EXECDIR}
	make clean;

exinfo: chkroot
	-@echo "" 
	-@echo "Compiler version: ${CC}"
	-@echo "Compiler flags: ${CFLAGS}"
	-@echo "Using Linker: ${CLINKER}"
	-@echo "Using libraries: ${EXAMPLE_LIBS}" 
	-@echo "" 

most: exinfo
	make maxcut
	make theta
	make dsdpsdpa

all:
	make most
	make dsdp
	make clean;

clean:
	${RM} *~ *.o allc/* *.h


