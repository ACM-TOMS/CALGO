#!/bin/bash

convert(){

	echo "-------------------------- configuring for C++ --------------------------"
	cd $2
#	CC='g++' CFLAGS='-I/home/wschrep/ArithmosRelease/Libraries/Arithmos/include -I/home/wschrep/mpfr-test/gmp-3.1.1' LIBS='-L/home/wschrep/ArithmosRelease/Libraries/Arithmos/lib -L/home/wschrep/gmp-3.1.1/lib -lArithmos -lgmp' ./configure 
#Better configuration, also works on solaris now (31/5/2006), the O2 is necessary there!!! alse -m32 might be necessary...
#	CC='g++' CFLAGS='-O2 -I/home/wschrep/ArithmosRelease/Libraries/Arithmos/include -I/home/wschrep/mpfr-test/gmp-3.1.1' LIBS='-L/home/wschrep/ArithmosRelease/Libraries/Arithmos/lib -L/home/wschrep/gmp-3.1.1/lib -lArithmos -lMpIeee -lSpecialValue -lgmp' ./configure
	CC='g++' CFLAGS='-O2 -I/home/wschrep/ArithmosRelease/Libraries/Arithmos/include -I/home/wschrep/mpfr-test/gmp-3.1.1' LIBS='' ./configure
	cd ..

	echo "-------------------------- precompile pass ------------------------------"
	./precompilearith $1 $2

	echo "-------------------- making fresh working copy of cblas -----------------"
	rm -rf "$2/cblas"
	cp -rf "$1/cblas" "$2/cblas"
	
	
	echo "------------------------- precompile cblas ------------------------------"
	./precompilecblas "$1/cblas" "$2/cblas"

	
	echo "-------------------- gsl-convert-cblas (index types) --------------------"
	cp -rf $2/cblas $2/cblas_pre
	./gsl-convert-cblas "$2/cblas_pre" "$2/cblas"
	
	echo "--------------------- converting rootdir --------------------------------"
	./precompileonedir "$1" "$2"

	echo "--------------------- converting templates  --------------------------------"
	./gsl-convert-templates $1/templates_on.h $2/templates_on.h
}


if (( $#!=2 ))
then
	echo "USAGE: $0 <gsl original> <gsl precompiled>";
	exit 1;
fi

if [ -d $2 ]
then
	echo "$2 already exists, using it..."
else
	echo "$2 does not exist, making copy..."
	cp -rf $1 $2
fi

convert $1 $2
