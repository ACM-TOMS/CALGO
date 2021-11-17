#- GNU g77 compiler
#FC = g77
#FFLAGS = -pedantic -W

#- Intel ifc compiler
#FC = ifc -L/usr/local/old_lib -lc-2.3.2
#FFLAGS = -FI -72 -C

#- Nag f95 compiler
#FC = f95
#FFLAGS = -C -gline

#- Lahey LF95 compiler
#FC = lf95
#FFLAGS = --chk a,s,u,x --warn

driver1: driver1.o src.o
	$(FC) -o run_driver1 $(FFLAGS) driver1.o src.o
	run_driver1
driver1.o: driver1.f
	$(FC) -c -o driver1.o $(FFLAGS) driver1.f
src.o: src.f
	$(FC) -c -o src.o $(FFLAGS) src.f
driver2: driver2.o rnd.o
	$(FC) -o run_driver2 $(FFLAGS) driver2.o rnd.o
	run_driver2
driver2.o: driver2.f
	$(FC) -c -o driver2.o $(FFLAGS) driver2.f
rnd.o: rnd.f
	$(FC) -c -o rnd.o $(FFLAGS) rnd.f
