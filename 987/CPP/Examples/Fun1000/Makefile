SRC=../../../Src

FILES=Parameters.o Manbis.o Functions.o

vpath %.cpp $(SRC)
vpath %.h $(SRC)

.PHONY: fun1000

fun1000: $(FILES) main.cpp
	$(CXX) $(CXXFLAGS) $(FILES) -I $(SRC) main.cpp -o fun1000 $(LDFLAGS)
	./fun1000 -l 0 -r 1 -e 8 -q 60
	mv Results.txt ResultsRun1.txt
	./fun1000 -l 0 -r 1 -e 8 -q 65
	mv Results.txt ResultsRun2.txt
	./fun1000 -l 0 -r 1 -e 8 -q 80
	mv Results.txt ResultsRun3.txt
	./fun1000 -l 0 -r 1 -e 8 -q 90
	mv Results.txt ResultsRun4.txt
	./fun1000 -l 0 -r 1 -e 8 -q 99
	mv Results.txt ResultsRun5.txt

Parameters.o: 

Manbis.o:

Functions.o: 

clean: 
	rm -f fun1000 *.o
	rm -f ResultsRun*.txt
