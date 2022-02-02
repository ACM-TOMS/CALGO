##################################################################################################################################
################################################ MAKEFILE TEMPLATE ###############################################################
##################################################################################################################################

#########
# SHELL #
#########

SHELL = /bin/sh

##########
# MACROS #
##########

SRC = <path to modules>
BLAS = <path to BLAS library>/libblas.a
LAPACK = <path to LAPACK library>/liblapack.a
BLACSF77 = <path to BLACS library>/blacsF77init_MPI-LINUX-0.a
BLACSMPI = <path to BLACS library>/blacs_MPI-LINUX-0.a
SCALAPACK = <path to SCALAPACK library>/libscalapack.a
MPICH = <path to MPI library>/libmpich.a

######################
# COMPILER AND FLAGS #
######################

#### MPI ( for ifort , pgf90 , g95 , gfortran ) 
FC = mpif90
FLAGS = -I<path to MPICH>/include
EXT = o

##################################################################################################################################
########################################## END OF USER MODIFICATION ##############################################################
##################################################################################################################################

############
# COMMANDS #
############

CD = cd
RM = rm -fv

###########
# OBJECTS #
###########

OBJECTS_MAIN =               	main.o

OBJECTS_MOD =                   module_precision.o\
				module_tools.o\
				module_random.o\
				module_model.o\
				module_observations.o\
				module_parallel.o\
				module_initialize.o\
				module_ekf.o\
				module_enkf.o\
				module_rrsqrtkf.o\
				module_rrsqrtenkf.o

###########
# TARGETS #
###########

TARGET =                    	 main.exe

#########
# RULES #
#########

default:
				$(MAKE) clean
				$(MAKE) compile

clean:
				$(RM) *~ *.$(EXT) *.mod *.exe *.txt *.eps *.avi
				$(CD) $(SRC) ; $(RM) *~ *.$(EXT) *.mod
				$(RM) mpif.h 

compile:
				$(MAKE) $(TARGET)

module_precision.o:		module_precision.f90
				$(FC) -c $(FLAGS) $< -o $@

module_tools.o:			module_tools.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_random.o: 		module_random.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_model.o: 		module_model.f90 module_precision.o 
				$(FC) -c $(FLAGS) $< -o $@

module_observations.o: 		module_observations.f90 module_precision.o
				$(FC) -c $(FLAGS) $< -o $@

module_parallel.o:		module_parallel.f90 module_precision.o module_model.o module_observations.o
				$(FC) -c $(FLAGS) $< -o $@

module_initialize.o: 		module_initialize.f90 module_precision.o 
				$(FC) -c $(FLAGS) $< -o $@

module_ekf.o: 			module_ekf.f90 module_precision.o module_parallel.o module_model.o module_observations.o module_initialize.o
				$(FC) -c $(FLAGS) $< -o $@

module_enkf.o: 			module_enkf.f90 module_precision.o module_random.o module_model.o module_observations.o module_parallel.o module_initialize.o
				$(FC) -c $(FLAGS) $< -o $@

module_rrsqrtkf.o: 		module_rrsqrtkf.f90 module_precision.o module_model.o module_observations.o module_parallel.o module_initialize.o
				$(FC) -c $(FLAGS) $< -o $@

module_rrsqrtenkf.o: 		module_rrsqrtenkf.f90 module_precision.o module_random.o module_model.o module_observations.o module_parallel.o module_initialize.o 
				$(FC) -c $(FLAGS) $< -o $@

main.o: 			main.f90 module_precision.o module_random.o module_tools.o module_model.o module_observations.o module_parallel.o module_initialize.o module_ekf.o module_enkf.o module_rrsqrtkf.o module_rrsqrtenkf.o
				$(FC) -c $(FLAGS) $< -o $@

main.exe:	        	$(OBJECTS_MOD) $(OBJECTS_MAIN)
				$(FC) -o $(TARGET) $(OBJECTS_MAIN) $(OBJECTS_MOD) $(SCALAPACK) $(BLACSF77) $(BLACSMPI) $(BLACSF77) $(LAPACK) $(BLAS)
