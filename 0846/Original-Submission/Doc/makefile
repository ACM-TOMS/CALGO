#compiler executable name
CC          = g++
# compile-only flag
CFLAG   = -c  -O2  -ansi
# link flag
LDFLAG  = -O2

###########################################################################
# Common macro declarations.
###########################################################################
# Main object compile rule.
#
.SUFFIXES: .cpp .o .h

.cpp.o:
	$(CC) $(CFLAG) $<

mixedvol: MixedVolDriver.o \
	PolyReader/GenericPolynomial.o PolyReader/PolynomialException.o \
	PolyReader/ExpandedPolynomial.o PolyReader/PolynomialSystemReader.o
	$(CC) $(LDFLAG) MixedVolDriver.o \
	GenericPolynomial.o PolynomialException.o \
	ExpandedPolynomial.o PolynomialSystemReader.o  -o $@

clean:
	rm -f *.o; rm -f GetOption/*.o; rm -f PolyReader/*.o; rm -f PreProcess/*.o; rm -f MixedVol/*.o;
