
D1       =  solver
D2       =  sdp
D3       =  lp
D4       =  vecmat
D5       =  sys
D7       =  bounds

dsdplib:
	cd ${D1}; make dsdplib;
	cd ${D2}; make dsdplib;
	cd ${D3}; make dsdplib;
	cd ${D4}; make dsdplib;
	cd ${D5}; make dsdplib;
	cd ${D7}; make dsdplib;

clean:
	cd ${D1}; make clean;
	cd ${D2}; make clean;
	cd ${D3}; make clean;
	cd ${D4}; make clean;
	cd ${D5}; make clean;
	cd ${D7}; make clean;

all: dsdplib 

default:
	make dsdplib;
