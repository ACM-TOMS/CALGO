#---------------------------------------------------#
#             Makefile for Genus                    #
#---------------------------------------------------#
# compiler
F90 = ifort
F90 = nagfor

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
genus : printmsgA $(OBJECTS)
	@echo Linking...
	@$(F90) -o genus $(OBJECTS) $(LFLAGS)
	@echo genus  built, congratulations.

#---------------------------#
#           message         #
#---------------------------#
printmsgA :
	@echo
	@echo Building genus for $(SYSTEM)
	@echo Compiler flags : $(F90FLAGS)
	@echo Linker   flags : $(LFLAGS)

#---------------------------#
#           clean           #
#---------------------------#
clean:
	find . \( -name '*.o' \) -exec rm {} \;
	find . \( -name '*__genmod.f90' \) -exec rm {} \;	
	find . \( -name '*.mod' \) -exec rm {} \;	
