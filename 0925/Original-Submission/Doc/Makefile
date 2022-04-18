
############################################################
#default environments
# you edit here along your environment
############################################################


# do not change the number of spaces before and after '='
CC       = mpicc
CXX      = mpicxx
FC       = mpif90
CFLAGS   = -O3
CXXFLAGS = -O3
FFLAGS   = -O3
AR       = ar
RANLIB   = ranlib

LAPACK_LIBS = /usr/local/lapack/lib/liblapack.a
BLAS_LIBS   = /usr/local/lapack/lib/libblas.a

SCALAPACK_LIBS  = /usr/local/scalapack/lib/libscalapack.a /usr/local/BLACS/LIB/blacs_MPI-LINUX-0.a /usr/local/BLACS/LIB/blacsCinit_MPI-LINUX-0.a /usr/local/BLACS/LIB/blacs_MPI-LINUX-0.a

FORTRAN_LIBS = -lgfortran
PTHERAD_LIBS = -lpthread

# F77_FUNC 
DEF_F77_FUNC = -DF77_FUNC\(name,NAME\)=name\ \#\#\ _ 
#DEF_F77_FUNC = -DF77_FUNC\(name,NAME\)=name\ \#\#\ __
#DEF_F77_FUNC = -DF77_FUNC\\\(name,NAME\\\)=name\\\ \\\#\\\#\\\ _
#DEF_F77_FUNC = -DF77_FUNC\\\(name,NAME\\\)=name\\\ \\\#\\\#\\\ __


############################################################
# After here, you do not need to edit
############################################################

.PHONY: all exe lib clean cleanall update_make_headers force mumps

MUMPS_DIR = ./mumps/build

EXE = sdpara
LIB = libsdpara.a
SRC = sdpa_mpist.cpp sdpa_mpicopy.cpp \
	sdpa_block.cpp sdpa_chordal.cpp sdpa_dataset.cpp \
	sdpa_dpotrf.cpp sdpa_io.cpp \
	sdpa_jordan.cpp sdpa_linear.cpp \
	sdpa_newton.cpp \
	sdpa_parts.cpp sdpa_struct.cpp \
	sdpa_tool.cpp sdpa_rpdpotrf.cpp

OBJ = $(SRC:%.cpp=%.o)
all: ${EXE} 

exe: ${EXE}
lib: ${LIB}

${EXE}: sdpa_exe.o ${LIB}
	${CXX} ${CXXFLAGS} -o $@ sdpa_exe.o ${LIB} \
		-L${MUMPS_DIR}/lib -ldmumps -lmumps_common -lpord \
		${SCALAPACK_LIBS} \
		${LAPACK_LIBS} ${BLAS_LIBS} \
		${FORTRAN_LIBS} ${PTHERAD_LIBS}

${LIB}: sdpa_solve.o sdpa_call.o ${OBJ}
	ar r $@ sdpa_solve.o sdpa_call.o ${OBJ}

.cpp.o:
	${CXX} -c ${CXXFLAGS} ${DEF_F77_FUNC} \
		-I${MUMPS_DIR}/include -o $@ $<

sdpa_solve.o: sdpa_solve.cpp
	${CXX} -c ${CXXFLAGS} ${DEF_F77_FUNC} \
		-I. -I${MUMPS_DIR}/include -o $@ $<
sdpa_call.o: sdpa_call.cpp
	${CXX} -c ${CXXFLAGS} ${DEF_F77_FUNC} \
		-I. -I${MUMPS_DIR}/include -o $@ $<
sdpa_exe.o: sdpa_exe.cpp
	${CXX} -c ${CXXFLAGS} ${DEF_F77_FUNC} \
		-I. -I${MUMPS_DIR}/include -o $@ $<

clean:
	rm -f *.o *~
cleanall: clean mumps-clean 
	rm -f ${EXE} ${LIB}

force: cleanall all

mumps:
	cd mumps; make

mumps-clean:
	cd mumps; make distclean

update_make_headers:
	g++ -MM *.cpp > make.headers
include make.headers
