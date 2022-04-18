#! /bin/sh
###=====================================================================
### Loop over a set of hosts defined in HOSTLIST, in parallel, copying
### the distribution tree to each remote host in a temporary build
### directory, and then launch the run-test.sh script on each to test
### all available Fortran compilers.
###
### Usage:
### 	run-par.sh [--hosts "host1 host2 ..."] target1 target2 ...
###
### [20-Apr-2000]
###=====================================================================

if test "$1" = "--hosts"
then
	shift
	hosts="$1"
	shift
else
	hosts=`cat HOSTLIST || true`
fi

if test -z "$hosts"
then
	echo There are no hosts defined in HOSTLIST
	exit 1
fi

BUILDDIR=/usr/tmp
THISDIR=`basename $PWD`
TARGETS="$@"

echo "Date: `date`"
echo
echo "Beginning parallel builds on:"
echo
for h in $hosts
do
	echo "	"$h
	LOGFILE=typescript.$h.`date +%Y.%m.%d.%H.%M.%S`
	( \
		cd ../..; tar cf - dist | \
			ssh $h "chdir $BUILDDIR ; \
				rm -rf dist ; \
				tar xf - ; \
				chdir dist/$THISDIR ; \
				./run-tests.sh $TARGETS " \
	) > $LOGFILE 2>&1 &
done
echo 
echo Waiting for completion of builds...
echo 
wait
ls -l typescript.*
echo "Date: `date`"



