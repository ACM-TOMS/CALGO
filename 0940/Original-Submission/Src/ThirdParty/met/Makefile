# -*- Makefile -*-
#
#  include ${HOME}/makehead/make$(EXT)

LINKOPT=-lm
LINKCOMPLEX=
FORTLINK=-lg2c

CPP=g++ -DUNDERF77
F77=g77

RANGECHECK = 
#RANGECHECK = -DRANGECHECK

INCLUDE = -Icommon
CPPI = $(CPP) $(INCLUDE) ${RANGECHECK}
#CPPI = $(CPP) $(INCLUDE) -DNUMRECIPES ${RANGECHECK}

all : xpr1Test$(EXT).out xpr1SpeedTest$(EXT).out xpr1SpeedTest2$(EXT).out xprDTest$(EXT).out xpr2Test$(EXT).out xpr2SpeedTest$(EXT).out xpr2DTest$(EXT).out # solverTest$(EXT).out leastsquareTest$(EXT).out
alltest: all cxxtest$(EXT) cxxtest2$(EXT).out cxxtest3$(EXT) cxxtest4$(EXT)

clean :
	rm xpr1Test$(EXT) xpr1SpeedTest$(EXT) xpr1SpeedTest2$(EXT) xprDTest$(EXT) xpr2Test$(EXT) xpr2SpeedTest$(EXT) xpr2DTest$(EXT)

# numrecipes$(EXT).o : numrecipes.f
# 	$(F77) -o $@ -c $?

xpr1Test$(EXT) : xpr1Test.C xpr1.h xprD.h xpr2.h vecmat.h
	$(CPPI) -o $@ xpr1Test.C $(LINKOPT)


xpr1Test$(EXT).out : xpr1Test$(EXT)
	$? > $@

xpr1SpeedTest$(EXT).o : xpr1SpeedTest.f
	$(F77) -o $@ -c $?

xpr1SpeedTest$(EXT) : xpr1SpeedTest.C xpr1.h xprD.h xpr2.h vecmat.h xpr1SpeedTest$(EXT).o
	$(CPPI) -o $@ xpr1SpeedTest.C xpr1SpeedTest$(EXT).o $(LINKOPT) ${FORTLINK}
#	$(CPPI) -DMEMCHECK -o $@ xpr1SpeedTest.C xpr1SpeedTest$(EXT).o $(LINKOPT) ${FORTLINK}

xpr1SpeedTest$(EXT).out : xpr1SpeedTest$(EXT)
	xpr1SpeedTest$(EXT) 1000 2000 4000 8000 > $@
#	xpr1SpeedTest$(EXT) 1000 2000 4000 8000 16000 > $@

xpr1SpeedTest2$(EXT).o : xpr1SpeedTest2.f
	$(F77) -o $@ -c $?

xpr1SpeedTest2$(EXT) : xpr1SpeedTest2.C xpr1.h xprD.h xpr2.h vecmat.h xpr1SpeedTest2$(EXT).o
	$(CPPI) -o $@ xpr1SpeedTest2.C xpr1SpeedTest2$(EXT).o $(LINKOPT) ${FORTLINK}
#	$(CPPI) -DMEMCHECK -o $@ xpr1SpeedTest2.C xpr1SpeedTest2$(EXT).o $(LINKOPT) ${FORTLINK}

xpr1SpeedTest2$(EXT).out : xpr1SpeedTest2$(EXT)
	xpr1SpeedTest2$(EXT) 1000 2000 4000 8000 > $@
#	xpr1SpeedTest2$(EXT) 1000 2000 4000 8000 16000 > $@

xprDTest$(EXT) : xprDTest.C xpr1.h xprD.h xpr2.h vecmat.h
	$(CPPI) -o $@ xprDTest.C $(LINKOPT)


xprDTest$(EXT).out : xprDTest$(EXT)
	$? > $@

xpr2Test$(EXT) : xpr2Test.C xpr1.h xprD.h xpr2.h vecmat.h
	$(CPPI) -o $@ xpr2Test.C $(LINKOPT)

xpr2Test$(EXT).out : xpr2Test$(EXT)
	$? > $@

#xpr2BinOpTest$(EXT) : xpr2BinOpTest.C xpr1.h xprD.h xpr2.h vecmat.h
#	$(CPPI) -o $@ xpr2BinOpTest.C $(LINKOPT)

xpr2SpeedTest$(EXT).o : xpr2SpeedTest.f
	$(F77) -o $@ -c $?

xpr2SpeedTest$(EXT) : xpr2SpeedTest.C xpr1.h xprD.h xpr2.h vecmat.h xpr2SpeedTest$(EXT).o
	$(CPPI) -o $@ xpr2SpeedTest.C xpr2SpeedTest$(EXT).o $(LINKOPT) ${FORTLINK}
#	$(CPPI) -DMEMCHECK -o $@ xpr2SpeedTest.C xpr2SpeedTest$(EXT).o $(LINKOPT) ${FORTLINK}

xpr2SpeedTest$(EXT).out : xpr2SpeedTest$(EXT)
	xpr2SpeedTest$(EXT) 100 200 400 800 > $@
#	xpr2SpeedTest$(EXT) 100 200 400 800 1600 3200 6400 > $@

xpr2DTest$(EXT) : xpr2DTest.C xpr1.h xprD.h xpr2.h vecmat.h
	$(CPPI) -o $@ xpr2DTest.C $(LINKOPT)

xpr2DTest$(EXT).out : xpr2DTest$(EXT)
	$? > $@

# solverTest$(EXT) : solverTest.C xpr1.h xprD.h xpr2.h vecmat.h numrecipes$(EXT).o solver.h
# 	$(CPPI) -o $@ solverTest.C numrecipes$(EXT).o $(LINKOPT) $(FORTLINK)

# solverTest$(EXT).out : solverTest$(EXT)
# 	$? > $@

# leastsquareTest$(EXT) : leastsquareTest.C xpr1.h xprD.h xpr2.h vecmat.h numrecipes$(EXT).o solver.h leastsquare.h
# 	$(CPPI) -o $@ leastsquareTest.C numrecipes$(EXT).o $(LINKOPT) $(FORTLINK)

# leastsquareTest$(EXT).out : leastsquareTest$(EXT)
#	$? > $@

newtest$(EXT) : newtest.C
	$(CPPI) -o $@ newtest.C

newtest$(EXT).out : newtest$(EXT)
	 newtest$(EXT) 1000 2000 3000 4000 > $@

cxxtest$(EXT) : cxxtest.C
	$(CPPI) -o $@ cxxtest.C

cxxtest2$(EXT) : cxxtest2.C
	$(CPPI) -o $@ cxxtest2.C $(LINKOPT)

cxxtest2$(EXT).out : cxxtest2$(EXT)
	cxxtest2$(EXT) > $@

cxxtest3$(EXT) : cxxtest3.C
	$(CPPI) -o $@ cxxtest3.C $(LINKOPT)

cxxtest4$(EXT) : cxxtest4.C
	$(CPPI) -o $@ cxxtest4.C $(LINKOPT)

