SRC=../../../Src

FILES=Parameters.o Manbis.o Functions.o

vpath %.cpp $(SRC)
vpath %.h $(SRC)

.PHONY: example2

example2: $(FILES) main.cpp
	$(CXX) $(CXXFLAGS) $(FILES) -I $(SRC) main.cpp -o example2 $(LDFLAGS)
	./example2 -l -5 -r 5 -e 8 -q 85 -a 16

Parameters.o: 

Manbis.o:

Functions.o: 

clean: 
	rm -f example2 *.o
	rm -f Results.txt
