#! /bin/sh
# Compile and run one or more short Fortran program on any local
# architecture, with all available compilers:
#
# Usage:
#	testit.sh foo1.f [foo2.f ... foon.f]
#
# [12-Jun-2000]

libdir=`dirname $0`/../common
system=`$libdir/machid.sh`

for f in "$@"
do
	case $system in
		Apple-ppc-Linux*)
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			;;
		Compaq_DEC*)
			rm -f a.out ; echo ======= f77; f77 -fpe3 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 -mieee-with-inexact $f && ./a.out
			rm -f a.out ; echo ======= f90; f90 -fpe3 $f && ./a.out
			rm -f a.out ; echo ======= f95; f95 -fpe3 $f && ./a.out
			rm -f a.out ; echo ======= nagf90; nagf90 -ieee=full -fixed $f && ./a.out
			rm -f a.out ; echo ======= nagf95; nagf95 -ieee=full -fixed $f && ./a.out
			;;
		HP*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			;;
		IBM*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			rm -f a.out ; echo ======= xlf90; xlf90 -qfixed $f && ./a.out
			rm -f a.out ; echo ======= xlf95; xlf95 -qfixed $f && ./a.out
			;;
		Intel*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			rm -f a.out ; echo ======= nagf90; nagf90 -ieee=full -fixed $f && ./a.out
			rm -f a.out ; echo ======= nagf95; nagf95 -ieee=full -fixed $f && ./a.out
			rm -f a.out ; echo ======= pgf77; pgf77 $f && ./a.out
			rm -f a.out ; echo ======= pgf90; pgf90 $f && ./a.out
			;;
		SGI*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			rm -f a.out ; echo ======= f90; f90 $f && ./a.out
			rm -f a.out ; echo ======= nagf90; nagf90 -ieee=full -fixed $f && ./a.out
			rm -f a.out ; echo ======= nagf95; nagf95 -ieee=full -fixed $f && ./a.out
			;;
		Sun-sparc-Linux*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			;;
		Sun*)
			rm -f a.out ; echo ======= f77; f77 $f && ./a.out
			rm -f a.out ; echo ======= g77; g77 $f && ./a.out
			rm -f a.out ; echo ======= f90; f90 -ftrap=%none $f && ./a.out
			rm -f a.out ; echo ======= f95; f95 -ftrap=%none $f && ./a.out
			rm -f a.out ; echo ======= nagf90; nagf90 -ieee=full -fixed $f && ./a.out
			rm -f a.out ; echo ======= nagf95; nagf95 -ieee=full -fixed $f && ./a.out
			;;
		*)
			echo Unrecognized system: $system
			exit 1
			;;
	esac
	rm -f a.out
done
