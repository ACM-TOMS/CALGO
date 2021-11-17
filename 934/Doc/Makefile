F95=nagfor
F95FLAGS=-gline -u -dcfuns -nan -C -C=undefined
%.o : %.mod

%.o: %.f95
	$(F95) $(F95FLAGS) -o $@ -c $<

%.o: %.f90
	$(F95) $(F95FLAGS) -o $@ -c $<

%: %.o
	$(F95) $(F95FLAGS) -o $@ $^ $(LIBS)

OBJS = bessel.o Blanch_Complex.o constants.o \
       feigen.o port-90.o

Driver: $(OBJS)

bessel.o: 

constants.o:

port-90.o: 

Blanch_Complex.o: constants.o feigen.o 

feigen.o: constants.o

Driver.o: bessel.o Blanch_Complex.o constants.o feigen.o port-90.o

clean:
	rm -rf *.o *.mod Driver
