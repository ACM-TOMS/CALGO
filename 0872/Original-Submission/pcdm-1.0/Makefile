# The makefile for PCDM
# Andrey Chernikov, January 2006

CXX                = mpiCC
CXX_FLAGS          = -DANSI_DECLARATORS -O3
PCDM_LIBRARY_FLAGS = -DPCDM_LIBRARY 

SOFTWARE     = ..
PREDICATES   = ${SOFTWARE}
METIS        = ${SOFTWARE}/metis-4.0
TRIANGLE     = ${SOFTWARE}/triangle-1.6

TRIANGLE_O   = ${TRIANGLE}/triangle.o
PREDICATES_O = ${PREDICATES}/predicates.o

INCLUDE  = -I./ 
INCLUDE += -I${TRIANGLE}
INCLUDE += -I${METIS}
INCLUDE += -I${METIS}/METISLib
INCLUDE += -I${METIS}/ParMETISLib

all : pcdm

pcdmlibrary : pcdm.cc pcdm.h pcdmio
	${CXX} ${CXX_FLAGS} ${PCDM_LIBRARY_FLAGS} ${INCLUDE} -c pcdm.cc -o pcdmlib.o

pcdm : pcdm.cc pcdm.h pcdmio
	${CXX} ${CXX_FLAGS} ${INCLUDE} -c pcdm.cc 
	${CXX} pcdm.o pcdmio.o ${PREDICATES_O} -o pcdm -L${METIS} -lmetis
        
driver : driver.cc pcdm.h pcdmlibrary pcdmio
	${CXX} ${CXX_FLAGS} ${INCLUDE} -c driver.cc
	${CXX} driver.o pcdmio.o pcdmlib.o ${PREDICATES_O} -o driver -L${METIS} -lmetis

cdt : cdt.cc pcdmio
	${CXX} ${CXX_FLAGS} ${INCLUDE} -c cdt.cc
	${CXX} cdt.o pcdmio.o ${TRIANGLE_O} -o cdt

pcdmio : pcdmio.cc pcdmio.h
	${CXX} ${CXX_FLAGS} ${INCLUDE} -c pcdmio.cc

clean :
	rm -f pcdm pcdm.o pcdmio.o pcdmlib.o driver driver.o cdt cdt.o core* *_*_*.e* *_*_*.o*
        
