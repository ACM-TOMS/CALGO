#---------------------------------------------------#
#             Makefile for Genus                    #
#---------------------------------------------------#
# compiler
F90 = ifort -coarray
#F90 = caf

# compiler flags
F90FLAGS = -c -O3
LFLAGS = -O3

# objects alpha
INIT = .
SRCDIR = $(INIT)
OBJ=$(INIT)
OBJDIR = $(OBJ)
CALCDIR =$(INIT)

OBJECTS = $(OBJ)/mainMod.o\
$(OBJ)/genusMod.o 

#---------------------------#
#          compile          #
#---------------------------#
$(OBJDIR)/mainMod.o : $(SRCDIR)/mainMod.f90  $(SRCDIR)/Makefile
	@echo compiling mainMod.f90
	@cd $(OBJDIR); $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/mainMod.f90
#
$(OBJDIR)/genusMod.o : $(SRCDIR)/genusMod.f90 $(SRCDIR)/Makefile  
	@echo compiling genusMod.f90
	@cd $(OBJDIR); $(F90) $(F90FLAGS) -I$(SRCDIR) $(SRCDIR)/genusMod.f90

mainMod.o:  genusMod.o

genusMod.o:

#---------------------------#
#           build           #
#---------------------------#
genus_coarrays : printmsgA $(OBJECTS)
	@echo Linking...
	@$(F90) -o genus_coarrays $(OBJECTS) $(LFLAGS)
	@echo genus_coarrays built, congratulations.

#---------------------------#
#           message         #
#---------------------------#
printmsgA :
	@echo
	@echo Building genus_coarrays for $(SYSTEM)
	@echo Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)

#---------------------------#
#           clean           #
#---------------------------#
clean:
	find . \( -name '*.o' \) -exec rm {} \;
	find . \( -name '*__genmod.f90' \) -exec rm {} \;	
	find . \( -name '*.mod' \) -exec rm {} \;	
