INCDIR = -I. -I.. -I../Src
LIBS   = -lm

GCCF="-g0 -O -ansi -pedantic -Wall -Wtraditional -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Woverloaded-virtual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -ffor-scope -fno-common -funroll-loops -finline-functions"

CCF="-g0 -O -ansi -pedantic -Wall -Wtraditional -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wconversion -Woverloaded-virtual -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -ffor-scope -fstrict-prototype -fno-common -funroll-loops -finline-functions"

CWF="-ansi strict -proto strict -w all -W unused -W nonotinlined, noimplicit"

all:
	@echo "To compile the test files do:"
	@echo ""
	@echo "\`\`make mac'' for c++     (Mac OSX gnu)        compiler"
	@echo "\`\`make gcc'  for g++     (egcs 1.1, gcc 2.95) compiler"
	@echo "\`\`make cxx'' for cxx     (DIGITAL C++ V6.1)   compiler"
	@echo "\`\`make kcc'' for KCC     (KAI 4.3)            compiler"
	@echo "\`\`make cc''  for CC      (SGI MIPSpro V7.30)  compiler"
	@echo "\`\`make mw''  for mwccppc (Metrowerks CW 7.0)  compiler"
	@echo ""
	@echo "To run the tests do:"
	@echo ""
	@echo "\`\`make run''"
	@echo ""
	@echo "To check the instantiation of the templated library"
	@echo "methods do:"
	@echo ""
	@echo "\`\`make mac.check''"
	@echo "\`\`make g++.check''"
	@echo "\`\`make cxx.check''"
	@echo "\`\`make kcc.check''"
	@echo "\`\`make cc.check''"
	@echo ""
	@echo "If all the methods are correctly instantiated no error or"
	@echo "warning messages should appear."
	@echo ""
	@echo "To clean up the directory do:"
	@echo ""
	@echo "\`\`make clean''"
	@echo ""

mac:
	make clean
	make CXX="c++" CXXFLAGS=${CCF} compile

gcc:
	make clean
	make CXX="g++" CXXFLAGS=${GCCF} compile

cxx:
	make clean
	make CXX="cxx" CXXFLAGS="-w2 -O4 -std strict_ansi" compile

kcc:
	make clean
	make CXX="KCC" CXXFLAGS="-O3 --strict" compile

cc:
	make clean
	make CXX="CC" CXXFLAGS="-O -LANG:std -OPT:Olimit=0" compile

mw:
	make clean
	make CXX="mwccppc" CXXFLAGS=${CWF} LIBS="" compile

mac.check:
	make CXX="c++" CXXFLAGS="-O" check

gcc.check:
	make CXX="g++" CXXFLAGS="-O" check

cxx.check:
	make CXX="cxx" CXXFLAGS="-w2 -O" check

kcc.check:
	make CXX="KCC" CXXFLAGS="-O --strict" check

cc.check:
	make CXX="CC" CXXFLAGS="-O -LANG:std -OPT:Olimit=0" check

mw.check:
	make CXX="mwccppc" CXXFLAGS=${CWF} LIBS="" check


compile: testall_triangle testall_quad testall_marker

run: 	
	./testall_triangle_l
	./testall_triangle 
	./testall_quad_l 
	./testall_quad 
	./testall_marker

diff: 	
	diff ./testall_triangle_l.out ./testall_triangle_l.expected
	diff ./testall_triangle.out   ./testall_triangle.expected
	diff ./testall_quad_l.out     ./testall_quad_l.expected
	diff ./testall_quad.out       ./testall_quad.expected
	diff ./testall_marker.out     ./testall_marker.expected

testall_triangle: testall_triangle.cc p2print.hh
	$(CXX) $(INCDIR) $(CXXFLAGS) -DLIST=true \
	testall_triangle.cc -o testall_triangle_l $(LIBS)
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DLIST=false \
	testall_triangle.cc -o testall_triangle $(LIBS)

testall_quad:  testall_quad.cc p2print.hh
	$(CXX) $(INCDIR) $(CXXFLAGS) -DLIST=true \
	testall_quad.cc -o testall_quad_l $(LIBS)
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DLIST=false \
	testall_quad.cc -o testall_quad $(LIBS)

testall_marker:  testall_marker.cc
	$(CXX) $(INCDIR) $(CXXFLAGS) \
	testall_marker.cc -o testall_marker $(LIBS)

check:
	make clean
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=true \
	-DP2MESH_DEBUG -DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=true \
	-DP2MESH_DEBUG -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=true \
	-DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=true \
	-c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=true \
	-DP2MESH_DEBUG -DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=true \
	-DP2MESH_DEBUG -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=true \
	-DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=true \
	 -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=false \
	-DP2MESH_DEBUG -DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=false \
	-DP2MESH_DEBUG -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=false \
	-DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=4 -DLIST=false \
	-c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=false \
	-DP2MESH_DEBUG -DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=false \
	-DP2MESH_DEBUG -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=false \
	-DP2MESH_VERBOSE -c checkall.cc
#
	$(CXX) $(INCDIR) $(CXXFLAGS) -DSIZE=3 -DLIST=false \
	-c checkall.cc
#
	make clean

clean:
	rm -f testall_triangle testall_quad testall_marker \
	testall_triangle_l testall_quad_l checkall *.out *.o cxx* *~
	rm -rf cxx* ti_* ii_*
