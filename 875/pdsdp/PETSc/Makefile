
BOPT=O
include ../../make.include
include ${PETSC_DIR}/bmake/common/base

ALLLIBS = ${DSDPLIB} ${LAPACKBLAS}

# Check compiler flags so PETSc and DSDP use the same compiler! g++ ?
dsdpsdpa:  readsdpa.o pdsdppetsc.o chkopts
	$(CLINKER) -o pdsdp5 readsdpa.o pdsdppetsc.o $(DSDPLIB) $(PETSC_KSP_LIB)
	rm *.o

all: dsdpsdpa


