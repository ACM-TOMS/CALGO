# makefile for VBF (Vector Boolean Functions) Jose A. Alvarez Cubero 2006

############# User settings ######################

# set your C++ compiler and options here...
GPP=g++
LIBS=-lntl
NTLINC= -I/usr/local/include -L/usr/local/lib

############# End of user settings ###############


ex: ex.cpp VBF.h
	$(GPP) $(NTLINC) -Wall ex.cpp -o ex.exe $(LIBS)

clean:
	rm -f *.exe
