CC = g++
CFLAGS = -O2 -g -Wall
SRC = common.h common.cc Interval.h Interval.cc
OBJTEST = common.o Interval.o Gabor.o input_maker.o ShiftGaborMP.o Partition.o FFT.o FFTGaborMP.o driver1.o
OBJDEMO = common.o Interval.o Gabor.o input_maker.o ShiftGaborMP.o Partition.o FFT.o FFTGaborMP.o driver2.o

LINKS = -lm -ldl

all: MPTest MPDemo


#################################################################

MPTest: $(OBJTEST)
	$(CC)  $(CFLAGS) -o MPTest $(OBJTEST) $(LINKS)

MPTest.o: $(SRC) MPTest.cc
	$(CC) $(CFLAGS) -c MPTest.cc

##################################################################

MPDemo: $(OBJDEMO)
	$(CC)  $(CFLAGS) -o MPDemo $(OBJDEMO) $(LINKS)

MPDemo.o: $(SRC) MPDemo.cc
	$(CC) $(CFLAGS) -c MPDemo.cc

###################################################################

###################################################################

common.o: common.h common.cc
	$(CC) -c $(CFLAGS) -o common.o common.cc

Interval.o: Interval.h Interval.cc common.h common.cc
	$(CC) -c $(CFLAGS) -o Interval.o Interval.cc

FFT.o: $(SRC) FFT.h FFT.cc
	$(CC) $(CFLAGS) -c FFT.cc

Gabor.o: $(SRC) Gabor.h Gabor.cc
	$(CC) $(CFLAGS) -c Gabor.cc

Partition.o: $(SRC) Partition.h Partition.cc
	$(CC) $(CFLAGS) -c Partition.cc

ShiftGaborMP.o: $(SRC) ShiftGaborMP.cc
	$(CC) $(CFLAGS) -c ShiftGaborMP.cc


FFTGaborMP.o: $(SRC) FFT.h FFT.cc FFTGaborMP.cc
	$(CC) $(CFLAGS) -c FFTGaborMP.cc

input_maker.o: $(SRC) input_maker.h input_maker.cc
	$(CC) $(CFLAGS) -c input_maker.cc


###################################################################

clean:
	rm -f results* core *.o MPTest MPDemo MPapp *~ 


