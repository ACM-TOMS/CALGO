
default: 
	make all

D1       =  src
D2       =  examples
D3	 =  matlab
D4       =  bin

include make.include

dsdplibrary: chkroot info
	cd ${D1}; make all
	${RANLIB} ${DSDPLIB}

example:
	cd ${D2}; make most

test:
	cd ${D4}; make all

dmatlab:
	cd ${D2}; make dsdp

dsdpmatlab: dsdplibrary
	make dmatlab

dsdpapi: dsdplibrary
	make example
	make test;

install:
	make dsdplibrary
	make example
	make test
	make dmatlab

all:
	make install

clean:
	cd ${D1}; make clean
	cd ${D2}; make clean
	cd ${D4}; make clean
	${RM} lib/lib* matlab/dsdp.mex*
	${RM} *~ */*~ */*/*~

htmlzip:
	zip -r DSDP5-api-html.zip dox
	${RM} -R dox

oshared: 
	-@${RM} tmp; \
	mkdir tmp; \
	cd tmp; \
	echo "building ${DSDPLIBSO}"; \
	${AR} x ${DSDPLIB} ;\
	${SH_LD} ${DSDPLIBSO} *.o -o ${DSDPLIBSO}; \
	cd ../ ; \
	${RM} tmp
