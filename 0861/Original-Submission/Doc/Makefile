%.o: %.f90
	$(F90) $(F77OPTS) -c $<
OBJS = Blanch.o Driver.o Mathieu_Zhang_Jin.o Utility.o

all : driver

driver:$(OBJS)
	$(F90) $(F77OPTS) $^ -o $@

Blanch.o : mathieu_zhang_jin.mod constants.mod 
blanch.mod : Blanch.o 
Driver.o : utility.mod mathieu_zhang_jin.mod blanch.mod constants.mod 
Mathieu_Zhang_Jin.o : constants.mod 
mathieu_zhang_jin.mod : Mathieu_Zhang_Jin.o 
Utility.o : constants.mod 
utility.mod : Utility.o 
constants.mod : constants.o 
