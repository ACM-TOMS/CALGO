#! /bin/sh
###---------------------------------------------------------------------
### This is a script to run one or more make targets for all available
### compilers on the current system.
###
### Usage:
###	./run-tests.sh arg1 arg2 ...
###
### All of the arguments are passed to make.
###
### Normally, this script is invoked for ./run-par.sh, which does
### PARALLEL executions of this script on each test system.
###
### [06-Apr-2000]
###---------------------------------------------------------------------

### Augment the search path with a location of a private version of
### ndiff, which is otherwise not yet installed on all of the test
### systems:
PATH=$PATH:$HOME/bin

hostname=`hostname`

os=`uname -s || true`
if test -z "$os"
then
	os=unknown
fi

TARGETS="$@"

run_one()
{
	## Usage: run_one COMPILER OPTION-1 OPTION-2 ...
	COMPILER="$1"
	shift

	(cd ../common; make distclean) >/dev/null 2>&1
	make distclean >/dev/null 2>&1

	echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	./configure
	echo ========================================================================
	echo 'DATE:        ' `date`
	echo 'HOSTNAME:    ' `hostname`
	echo 'MACHINETYPE: ' `$HOME/bin/machinetype || true`
	echo 'COMPILER:    ' $COMPILER
	echo 'OPTIONS:     ' "$@"
	make -i -k RELERR=1 FC="$COMPILER" FOPT="'""$@""'" $TARGETS
	echo @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
}

case $os in
	AIX)
		run_one f77 -O3
		run_one g77 -O3
		run_one xlf90 -O3 -qfixed
		if test `uname -v``uname -r` -gt 42	# AIX 4.3 and later have xlf95
		then
			run_one xlf95 -O3 -qfixed
		fi
		;;

	OSF1)			# DEC Alpha OSF/1 4.x
		run_one f77 -O3
		run_one f90 -O3
		run_one f95 -O3
		run_one g77 -O3
		run_one nagf95 -O4 -fixed -ieee=full
		;;

	HP-UX)			# HP 9000/7xx HP-UX 10.xx
		run_one f77 +O3
		run_one g77 -O3
		run_one nagf95 -O4 -fixed -ieee=full
		;;

	IRIX|IRIX64)
		case `uname -r` in
		5*)
			run_one f77 -O2 -Wf,-I../common
			run_one g77 -O3
			run_one nagf95 -O4 -fixed -ieee=full
			;;
		6*)
			run_one f77 -O3
			run_one f90 -O3
			run_one g77 -O3
			run_one nagf95 -O4 -fixed -ieee=full
			;;
		*)
			echo Unrecognized SGI IRIX version: `uname -r`
			exit 1
		esac
		;;

	Linux)			# GNU/Linux (Alpha, Intel, PowerPC, SPARC, ...)
		run_one f77 -O3
		run_one g77 -O3
		run_one nagf95 -O4 -fixed -ieee=full
		run_one pgf77 -O3
		run_one pgf90 -O3
		run_one pghpf -O3
		;;

	SunOS)
		run_one f77 -O3
		run_one f90 -O3 -ftrap=%none
		run_one g77 -O3
		run_one nagf95 -O4 -fixed -ieee=full
		;;

	*)
		echo Unrecognized operation system: $os
		exit 1
		;;
esac
