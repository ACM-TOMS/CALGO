f90= frt -Am -Hasu -fi -v9
f90= f90
f90= xlf90 
f90= epcf90 -d20 -C 

INT=f90d1mach.o ivl_def.o intaritb.o testsys.o intlib.o
TEST=f90d1mach.o ivl_def.o intaritb.o intlib.o test_f90_intarith.o
 
HSL12= /rutherford/num-arcu/hsl12/source/
HSL12t= /rutherford/num-arcu/hsl12/tests/

all:	test int
 
test:	$(TEST)
	$(f90)  $(TEST)
	a.out 
	diff TEST_F90_INTARITH.OUT tinysmpl.out

int:	$(INT)
	$(f90)  $(INT)
	a.out 
	diff INTARITH.OUT sample.out

f90d1mach.o:	f90d1mach.f90
	$(f90) -c f90d1mach.f90 

ivl_def.o:	ivl_def.f90
	$(f90) -c ivl_def.f90 

intaritb.o:	intaritb.f90
	$(f90) -c intaritb.f90 

testsys.o:	testsys.f90
	$(f90) -c -o testsys.o testsys.f90 

intlib.o:	intlib.f90
	$(f90) -c intlib.f90 

test_f90_intarith.o:	test_f90_intarith.f90
	$(f90) -c test_f90_intarith.f90 

