F95=sunf95
F95FLAGS=

F95=gfortran
F95FLAGS=

F95=nagf95
F95FLAGS=-C=all -nan -gline -g90 -f2003 
F95FLAGS=-C=all -C=undefined -nan -gline -f2003  -g90
F95FLAGS=

%.o : %.f90
	$(F95) -c $(F95FLAGS) $< -o $@

# assume the compiler treats .f files as fixed
%.o : %.f
	$(F95) -c $(F95FLAGS) $< -o $@

# turn off the rule that attempts (on some ancient systems) to
# compile .mod files as if they contained Modula code!
%.o : %.mod

driver: cons.o air.o parabw.o testw.o
	$(F95) $(F95FLAGS) -o driver cons.o air.o parabw.o testw.o

cons.o: cons.f90 constants

air.o: air.f90 cons.o

parabw.o: parabw.f90 air.o cons.o

testw.o: testw.f90 parabw.o

clean:
	rm -rf *.mod *.o *.g90 driver
