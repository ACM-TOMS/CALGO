.KEEP_STATE:

TARGET		= interpo
SRCS		= main.cc Interp.cc Net.cc
OBJS		= $(SRCS:%.cc=%.o)
# INCDIRS and LIBDIRS provide the pathnames for the netcdf include and library
#    subdirectories.  You will have to change these pathnames to coincide with
#    the area where you install this software.
INCDIRS		= /crunch/homes/minser/IMP/interp_code/code/netcdf/include 
LIBDIRS		= /usr/lib  /crunch/homes/minser/IMP/interp_code/code/netcdf/lib
LIBS		= -lm  -lnetcdf

#The C++ compiler is g++ - version 2.7.2.1
C++             = g++  -O2

%.o:	%.cc
	$(C++) $(INCDIRS:%=-I%) -c $<

$(TARGET):	$(OBJS)
#	quantify -collector=/usr/local/lib/gcc-lib/sparc-sun-sunos4.1/2.2.2/ld $(C++) $(OBJS) -o $(TARGET) $(LIBDIRS:%=-L%) $(LIBS) 
	$(C++) $(OBJS) -o $(TARGET) $(LIBDIRS:%=-L%) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS) 
