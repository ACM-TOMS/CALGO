
default: all

include ../../make.include

SCALAPACKDIR=/home/benson/soft/SCALAPACK
BLACSDIR=/home/benson/soft/BLACS/LIB

#MPIHOME = /home/benson/petsc-2.3.0/externalpackages/mpich2-1.0.1/linux-gnu-O
MPIHOME = /home/petsc/soft/linux-debian_sarge-gcc-pgf90/mpich-1.2.6
#MPIHOME = /home/petsc/soft/linux-2/mpich-1.2.5

DCFLAGS = -I${SCALAPACKDIR}/PBLAS/SRC  -I${MPIHOME}/include/

SCALAPACKLIB = -L${SCALAPACKDIR}  -lscalapack 

COMMPLAT=MPI-LINUX-0
#BLACSLIB = ${BLACSDIR}/blacsCinit_${COMMPLAT}.a  ${BLACSDIR}/blacsF77init_${COMMPLAT}.a  ${BLACSDIR}/blacs_${COMMPLAT}.a
BLACSLIB = -L${BLACSDIR}/ -lblacsCinit_${COMMPLAT}  -lblacsF77init_${COMMPLAT}  -lblacs_${COMMPLAT}

#MPI = -L${MPIHOME}/lib -lfmpich -lmpich 
MPI = -L${MPIHOME}/lib -lmpich
#MPI = -L${MPIHOME}/lib -lfmpich -lmpich -lmpichf90 -lmpichf90nc -lmpichfsup
#MPI = ${MPIHOME}/lib/*.a

LAPACKBLAS = /usr/lib/liblapack.a -lblas -lf2c -lm

ALLLIBS = ${DSDPLIB} ${SCALAPACKLIB} ${BLACSLIB} ${MPI} ${LAPACKBLAS} 

FLINKER = mpicc

OBJS	= readsdpa.o pdsdpscalapack.o
dsdpsdpa: ${OBJS}
	${FLINKER} -o pdsdp5 ${OBJS} ${ALLLIBS}
	${RM} ${OBJS}

OBJS2	= maxcut.o pdsdpscalapack.o
maxcut: ${OBJS2}
	${CLINKER} -o pmaxcut ${OBJS2} ${ALLLIBS}
	${RM} ${OBJS2}

info2:
	-@echo "======================================================"
	-@echo "Compiling PDSDP using SCLAPACK parallel linear solver."
	-@echo "  Check ../../make.include to set compiler."
	-@echo "  the mpicc compiler may produce faster DSDP libraries"
	-@echo "======================================================"

clean:
	${RM} ${OBJS}

all: info info2 dsdpsdpa



