SRC=../../../Src

FILES=Parameters.o Manbis.o Functions.o

vpath %.cpp $(SRC)
vpath %.h $(SRC)

example1: $(FILES) main.cpp
	$(CXX) $(CXXFLAGS) $(FILES) -I $(SRC) main.cpp -o example1 $(LDFLAGS)
	./example1 -l 0 -r 10 -e 5 -q 80
	mv Results.txt ResultsRun1.txt
	./example1 -l 0 -r 5 -e 6 -q 90 -n 0 -u 1
	mv Results.txt ResultsRun2.txt
	mv SubsByLength.txt SubsByLengthRun2.txt
	mv SubsByXLeft.txt SubsByXLeftRun2.txt

Parameters.o: 

Manbis.o:

Functions.o: 

clean: 
	rm -f example1 *.o
	rm -f SubsByLengthRun2.txt
	rm -f SubsByXLeftRun2.txt
	rm -f ResultsRun*.txt
