%.o : %.f
	$(F77) $(F77OPTS) -c $*.f

all: res1 res2 res3 res4 res5
OBJS1 = driver1.o src.o
driver1: $(OBJS1)
	cp data1 cg.parm
	$(F77) $(F77OPTS) -o driver1 $(OBJS1)
res1: driver1 data1
	driver1 <data1  > res1

OBJS2 = driver2.o src.o
driver2: $(OBJS2)
	cp data2 cg.parm
	$(F77) $(F77OPTS) -o driver2 $(OBJS2)
res2: driver2 data2
	driver2 <data2  > res2

OBJS3 = driver3.o src.o
driver3: $(OBJS3)
	cp data3 cg.parm
	$(F77) $(F77OPTS) -o driver3 $(OBJS3)
res3: driver3 data3
	driver3 <data3 > res3

OBJS4 = driver4.o src.o
driver4: $(OBJS4)
	cp data4 cg.parm
	$(F77) $(F77OPTS) -o driver4 $(OBJS4)
res4: driver4 data4
	driver4 <data4 > res4

OBJS5 = driver5.o src.o
driver5: $(OBJS5)
	cp data5 cg.parm
	$(F77) $(F77OPTS) -o driver5 $(OBJS5)
res5: driver5 data5
	driver5 <data5 > res5

clean:
	rm -rf *.o driver1 driver2 driver3 driver4 driver5 res *.mod *.d
