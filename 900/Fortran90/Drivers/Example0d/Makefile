#########
# SHELL #
#########

SHELL = /bin/sh

##########
# MACROS #
##########

SRC = /home/torres/Desktop/ACMTOMS-RERESUBMISSION-TORRES/programs/src
BLAS = /opt/lapack/libblas.a
LAPACK = /opt/lapack/liblapack.a
BLACSF77 = /opt/blacs/LIB/blacsF77init_MPI-LINUX-0.a
BLACSMPI = /opt/blacs/LIB/blacs_MPI-LINUX-0.a
SCALAPACK = /opt/scalapack/libscalapack.a
MPICH = /opt/mpich/lib/libmpich.a

######################
# COMPILER AND FLAGS #
######################

#### MPI ( for ifort , pgf90 , g95 , gfortran ) 
FC = mpif90
FLAGS = -I/opt/mpich/include
EXT = o

##############################################################################
########################## END OF USER MODIFICATION ########################## 
##############################################################################

############
# COMMANDS #
############

CD = cd
RM = rm -fv

###########
# OBJECTS #
###########

OBJECTS_MAIN =                  example0d.o
OBJECTS_MOD =                   module_precision.o\
				module_tools.o\
				module_random.o\
				module_parameters.o\
				module_solution.o\
				module_domain.o\
				module_model.o\
				module_observations.o\
				module_parallel.o\
				module_initialize.o\
				module_ekf.o\
				module_rrsqrtkf.o\
				module_enkf.o\
				module_rrsqrtenkf.o

###########
# TARGETS #
###########

TARGET =                        example0d.exe

#########
# RULES #
#########

default:
				$(MAKE) clean
				$(MAKE) compile
				mpirun -np 1 ./example0d.exe

clean:
				$(RM) *~ *.$(EXT) *.mod *.exe *.txt *.eps *.avi
				$(CD) $(SRC) ; $(RM) *~ *.$(EXT) *.mod
				$(RM) mpif.h 

compile:
				$(MAKE) $(TARGET)

module_precision.o:		$(SRC)/module_precision.f90
				$(FC) -c $(FLAGS) $< -o $@

module_tools.o:			$(SRC)/module_tools.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_random.o: 		$(SRC)/module_random.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_parameters.o: 		module_parameters.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_solution.o: 		module_solution.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_domain.o: 		module_domain.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_model.o: 		module_model.f90 module_precision.o module_solution.o module_domain.o
				$(FC) -c $(FLAGS) $< -o $@

module_observations.o: 		module_observations.f90 module_precision.o module_random.o module_solution.o module_domain.o
				$(FC) -c $(FLAGS) $< -o $@

module_parallel.o:		$(SRC)/module_parallel.f90 module_precision.o module_model.o module_observations.o
				$(FC) -c $(FLAGS) $< -o $@

module_initialize.o: 		module_initialize.f90 module_precision.o module_random.o module_solution.o module_domain.o
				$(FC) -c $(FLAGS) $< -o $@

module_ekf.o: 			$(SRC)/module_ekf.f90 module_precision.o module_model.o module_observations.o module_initialize.o module_parallel.o
				$(FC) -c $(FLAGS) $< -o $@

module_rrsqrtkf.o: 		$(SRC)/module_rrsqrtkf.f90 module_precision.o module_tools.o module_model.o module_observations.o module_initialize.o module_parallel.o
				$(FC) -c $(FLAGS) $< -o $@

module_enkf.o: 			$(SRC)/module_enkf.f90 module_precision.o module_random.o module_model.o module_observations.o module_initialize.o module_parallel.o
				$(FC) -c $(FLAGS) $< -o $@

module_rrsqrtenkf.o: 		$(SRC)/module_rrsqrtenkf.f90 module_precision.o module_tools.o module_random.o module_model.o module_observations.o module_initialize.o module_parallel.o
				$(FC) -c $(FLAGS) $< -o $@

example0d.o: 			example0d.f90 module_precision.o module_tools.o module_random.o module_parameters.o module_solution.o module_domain.o module_model.o module_observations.o module_initialize.o module_ekf.o module_rrsqrtkf.o module_enkf.o module_rrsqrtenkf.o module_parallel.o
				$(FC) -c $(FLAGS) $< -o $@

example0d.exe:	        	$(OBJECTS_MOD) $(OBJECTS_MAIN)
				$(FC) $(FLAGS) -o $(TARGET) $(OBJECTS_MOD) $(OBJECTS_MAIN) $(SCALAPACK) $(BLACSF77) $(BLACSMPI) $(BLACSF77) $(LAPACK) $(BLAS)
