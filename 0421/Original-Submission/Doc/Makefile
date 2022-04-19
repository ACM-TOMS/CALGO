%.o : %.f90
	$(F77) $(F77OPTS) -c $*.f90

%.o : %.f
	$(F77) $(F77OPTS) -c $*.f

OBJS = driver.o src.o

RES: data driver
	driver < data >RES

driver: $(OBJS)
	$(F77) $(F77OPTS) -o driver $(OBJS)
