INCDIR = -I. -I.. -I../Src
LIBS   = -lm

all:
	@echo "To compile examples files do"
	@echo ""
	@echo "\`\`make mac'' for c++     (Mac OSX gnu)        compiler"
	@echo "\`\`make gcc'' for g++     (egcs 1.1, gcc 2.95) compiler"
	@echo "\`\`make cxx'' for cxx     (DIGITAL C++ V6.1)   compiler"
	@echo "\`\`make kcc'' for KCC     (KAI 4.3)            compiler"
	@echo "\`\`make cc''  for CC      (SGI MIPSpro V7.30)  compiler"
	@echo "\`\`make mw''  for mwccppc (Metrowerks CW 7.0)  compiler"
	@echo ""
	@echo "To run all the examples do"
	@echo ""
	@echo "\`\`make run''"
	@echo ""
	@echo "To clean up the directory do"
	@echo ""
	@echo "\`\`make clean''"
	@echo ""


GCCF="-g0 -O -ansi -pedantic -Wall -Wtraditional -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Woverloaded-virtual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -ffor-scope -fno-common -funroll-loops -finline-functions"

CCF="-g0 -O -ansi -pedantic -Wall -Wtraditional -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Woverloaded-virtual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -ffor-scope -fstrict-prototype -fno-common -funroll-loops -finline-functions"

CWF="-ansi strict -proto strict -w all -W unused -W nonotinlined, noimplicit"

mac:
	make CXX="c++" CXXFLAGS=${CCF} compile


gcc:
	make CXX="g++" CXXFLAGS=${GCCF} compile

cxx:
	make CXX="cxx" CXXFLAGS="-w2 -O4 -std strict_ansi" compile

kcc:
	make CXX="KCC" CXXFLAGS="-O3 --strict" compile

cc:
	make CXX="CC" CXXFLAGS="-O -LANG:std -OPT:Olimit=0" compile

mw:
	make CXX="mwccppc" CXXFLAGS=${CWF} LIBS="" compile


compile: p2_solver q2_solver cc_solver cv_solver

run:
	./p2_solver
	./q2_solver
	./cc_solver
	./cv_solver

diff:
	diff p2.mtv.expected p2.mtv
	diff q2.mtv.expected q2.mtv
	diff cc.mtv.expected cc.mtv
	diff cv.mtv.expected cv.mtv

p2_solver: p2_solver.cc
	$(CXX) $(INCDIR) $(CXXFLAGS) p2_solver.cc -o p2_solver

q2_solver: q2_solver.cc
	$(CXX) $(INCDIR) $(CXXFLAGS) q2_solver.cc -o q2_solver

cc_solver: cc_solver.cc eu.o
	$(CXX) $(INCDIR) $(CXXFLAGS) eu.o cc_solver.cc -o cc_solver $(LIBS)

cv_solver: cv_solver.cc eu.o
	$(CXX) $(INCDIR) $(CXXFLAGS) eu.o cv_solver.cc -o cv_solver $(LIBS)

eu.o:	eu.cc eu.hh
	$(CXX) $(INCDIR) $(CXXFLAGS) -c eu.cc -o eu.o

clean:
	rm -f p2_solver q2_solver cc_solver cv_solver *.o *.mtv *~
	rm -rf cxx* ti_* ii_*
