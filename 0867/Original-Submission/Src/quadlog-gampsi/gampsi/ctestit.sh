#! /bin/sh
# Compile and run one or more short C or C++ programs on any local
# architecture, with all available compilers:
#
# Usage:
#	testit.sh foo1.c [foo2.c ... foon.c]
#
# [15-Jun-2000]

libdir=`dirname $0`/../common
system=`$libdir/machid.sh`

for f in "$@"
do
	case $system in
		Apple-ppc-Linux*)
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			;;
		Apple-Power_Macintosh-Rhapsody-*)
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= c++; c++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			;;
		Compaq_DEC*)
			rm -f a.out ; echo ======= c89; c89 -ieee_with_inexact $f && ./a.out
			rm -f a.out ; echo ======= cc; cc  -ieee_with_inexact $f && ./a.out
			rm -f a.out ; echo ======= cxx; cxx -x cxx -ieee_with_inexact $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ -mieee-with-inexact $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc -mieee-with-inexact $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			;;
		HP*)
			rm -f a.out ; echo ======= c89; c89 $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= CC; CC $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			;;
		IBM*)
			rm -f a.out ; echo ======= c89; c89  $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= xlc; xlc $f && ./a.out
			rm -f a.out ; echo ======= xlC; xlC $f && ./a.out
			;;
		Intel*)
			rm -f a.out ; echo ======= CC; CC  $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			rm -f a.out ; echo ======= pgCC; pgCC $f && ./a.out
			rm -f a.out ; echo ======= pgcc; pgcc $f && ./a.out
			;;
		NeXT-*)
			rm -f a.out ; echo ======= c++; c++  $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			;;
		SGI*)
			rm -f a.out ; echo ======= CC; CC  $f && ./a.out
			rm -f a.out ; echo ======= c89; c89 $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			;;
		Sun-sparc-Linux*)
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			;;
		Sun*)
			rm -f a.out ; echo ======= CC; CC $f && ./a.out
			rm -f a.out ; echo ======= c89; c89 $f && ./a.out
			rm -f a.out ; echo ======= cc; cc $f && ./a.out
			rm -f a.out ; echo ======= g++; g++ $f && ./a.out
			rm -f a.out ; echo ======= gcc; gcc $f && ./a.out
			rm -f a.out ; echo ======= lcc; lcc $f && ./a.out
			;;
		*)
			echo Unrecognized system: $system
			exit 1
			;;
	esac
	rm -f a.out
done
