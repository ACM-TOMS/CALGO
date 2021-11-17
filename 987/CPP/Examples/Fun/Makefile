SRC=../../../Src

FILES=Parameters.o Manbis.o Functions.o

vpath %.cpp $(SRC)
vpath %.h $(SRC)

.PHONY: fun

fun: $(FILES) main.cpp
	$(CXX) $(CXXFLAGS) $(FILES) -I $(SRC) main.cpp -o fun $(LDFLAGS)
	./fun -l -200.5 -r 100 -e 8 -q 30 -s 0.5
	mv Results.txt ResultsRun1.txt
	./fun -l -200.5 -r 100 -e 8 -q 30 -s 0.4
	mv Results.txt ResultsRun2.txt
	./fun -l -200.5 -r 100 -e 8 -q 30 -s 0.3
	mv Results.txt ResultsRun3.txt
	./fun -l -200.5 -r 100 -e 8 -q 33 -s 0.5
	mv Results.txt ResultsRun4.txt
	./fun -l -200.5 -r 100 -e 8 -q 40 -s 0.5
	mv Results.txt ResultsRun5.txt
	./fun -l -200.5 -r 100 -e 8 -q 43 -s 0.5
	mv Results.txt ResultsRun6.txt
	./fun -l -200.5 -r 100 -e 8 -q 60 -s 0.5
	mv Results.txt ResultsRun7.txt
	./fun -l -200.5 -r 100 -e 8 -q 65 -s 0.5
	mv Results.txt ResultsRun8.txt
	./fun -l -200.5 -r 100 -e 8 -q 80 -s 0.5
	mv Results.txt ResultsRun9.txt
	./fun -l -200.5 -r 100 -e 8 -q 90 -s 0.5
	mv Results.txt ResultsRun10.txt
	./fun -l -200.5 -r 100 -e 8 -q 95 -s 0.5
	mv Results.txt ResultsRun11.txt
	./fun -l -200.5 -r 100 -e 8 -q 97 -s 0.5
	mv Results.txt ResultsRun12.txt
	./fun -l -200.5 -r 100 -e 8 -q 99 -s 0.5
	mv Results.txt ResultsRun13.txt

Parameters.o: 

Manbis.o:

Functions.o: 

clean: 
	rm -f fun *.o
	rm -f ResultsRun*.txt
