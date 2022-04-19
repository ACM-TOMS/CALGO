#!/bin/bash

convert(){

	echo "-------------------------- configuring for C++ --------------------------"
	cd $2
	CC='g++' CFLAGS='-I/home/wschrep/ArithmosRelease/Libraries/Arithmos/include -I/home/wschrep/mpfr-test/gmp-3.1.1' LIBS='-L/home/wschrep/ArithmosRelease/Libraries/Arithmos/lib -L/home/wschrep/gmp-3.1.1/lib -lArithmos -lgmp' ./configure 
	cd ..

	echo "-------------------------- precompile pass ------------------------------"
	./precompilearith.sh $1 $2

	echo "-------------------- making fresh working copy of cblas -----------------"
	rm -rf "$2/cblas"
	cp -rf "$1/cblas" "$2/cblas"
	
	
	echo "------------------------- precompile cblas ------------------------------"
	./precompilecblas.sh "$1/cblas" "$2/cblas"

	
	echo "-------------------- gsl-convert-cblas (index types) --------------------"
	cp -rf $2/cblas $2/cblas_pre
	./gsl-convert-cblas.sh "$2/cblas_pre" "$2/cblas"
	
	echo "--------------------- converting rootdir --------------------------------"
	./precompileonedir.sh "$1" "$2"

	echo "--------------------- converting templates  --------------------------------"
	./gsl-convert-templates.sh $1/templates_on.h $2/templates_on.h
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
