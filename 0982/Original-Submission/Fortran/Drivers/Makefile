.SUFFIXES:

%.o : %.mod

%.o: %.f90
	$(F95) $(F95FLAGS) -I ../Src -o $@ -c $<

vpath %.f90 ../Src
vpath %.f9h ../Src

F95= nagfor
F95FLAGS=-O

OBJ= Decay.o Dose_Factors_m.o ENDF_m.o Triangular_ODE.o Triangularize_m.o Connected_m.o

all: Decay

Dose_Factors_m.o: 

ENDF_m.o:

Connected_m.o:

Triangular_ODE.o: Compute_Solution.f9h  Compute_Z_Matrix.f9h

Decay.o: ENDF_m.o Dose_Factors_m.o Connected_m.o Triangular_ODE.o Triangularize_m.o 

Decay: $(OBJ)
	$(F95) $(F95LINKFLAGS) -o $@ $^

clean:
	rm -f Decay *.o *.mod
