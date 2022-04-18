
# Fortran compiler and linker
F77         = epcf90
F77OPTS     = -sloppy -C -d5 -g -temp=/tmp -u
F77LINK     = $(F77)
F77LINKOPTS = $(F77OPTS)

.f.o:
	$(F77) $(F77OPTS) -c $*.f

all: res1 res2 res3 res4 res5 res6 res7 res8 res9

Objs1= driver1.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

src.o:  ../Src/src.f
	$(F77) $(F77OPTS) -c $?
port.o: ../Src/port.f
	$(F77) $(F77OPTS) -c $?
blas1.o: ../Src/blas1.f
	$(F77) $(F77OPTS) -c $?
linpack.o: ../Src/linpack.f
	$(F77) $(F77OPTS) -c $?

driver1: $(Objs1)
	$(F77LINK) $(F77LINKOPTS) -o driver1 $(Objs1) $(SRCLIBS) $(Libs1)
res1: driver1  
	driver1  >temp 
	diff temp res1


Libs2= $(RANDOM)
Objs2= driver2.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver2: $(Objs2)
	$(F77LINK) $(F77LINKOPTS) -o driver2 $(Objs2) $(SRCLIBS) $(Libs2)
res2: driver2  
	driver2  >temp 
	diff temp res2


Libs3= $(RANDOM)
Objs3= driver3.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver3: $(Objs3)
	$(F77LINK) $(F77LINKOPTS) -o driver3 $(Objs3) $(SRCLIBS) $(Libs3)
res3: driver3  
	driver3  >temp 
	diff temp res3


Libs4= $(RANDOM)
Objs4= driver4.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver4: $(Objs4)
	$(F77LINK) $(F77LINKOPTS) -o driver4 $(Objs4) $(SRCLIBS) $(Libs4)
res4: driver4 data4 
	driver4 <data4 >temp 
	diff temp res4


Libs5= $(RANDOM)
Objs5= driver5.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver5: $(Objs5)
	$(F77LINK) $(F77LINKOPTS) -o driver5 $(Objs5) $(SRCLIBS) $(Libs5)
res5: driver5 data5 
	driver5 <data5 >temp 
	diff temp res5


Libs6= $(RANDOM)
Objs6= driver6.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver6: $(Objs6)
	$(F77LINK) $(F77LINKOPTS) -o driver6 $(Objs6) $(SRCLIBS) $(Libs6)
res6: driver6 data6 
	driver6 <data6 >temp 
	diff temp res6


Libs7= $(RANDOM)
Objs7= driver7.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver7: $(Objs7)
	$(F77LINK) $(F77LINKOPTS) -o driver7 $(Objs7) $(SRCLIBS) $(Libs7)
res7: driver7  
	driver7  >temp 
	diff temp res7


Libs8= $(RANDOM)
Objs8= driver8.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver8: $(Objs8)
	$(F77LINK) $(F77LINKOPTS) -o driver8 $(Objs8) $(SRCLIBS) $(Libs8)
res8: driver8  
	driver8  >temp 
	diff temp res8


Libs9= $(RANDOM)
Objs9= driver9.o dr_subs.o src.o port.o portrand.o linpack.o blas1.o

driver9: $(Objs9)
	$(F77LINK) $(F77LINKOPTS) -o driver9 $(Objs9) $(SRCLIBS) $(Libs9)
res9: driver9  
	driver9  >temp 
	diff temp res9

clean:
	rm *.o driver? temp

