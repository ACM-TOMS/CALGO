#! /bin/sh
# Move figure data and plots to LaTeX figures directory tree.
#
# Usage:
#	./movefigs.sh oldname newname subdir
#
# E.g.
#	./movefigs.sh gam101 qgamc100 dec
#
# [20-Jun-2000]

FIG=/u/sy/beebe/jim-ball/toms/latex/figures

if test -n "$1"
then
	IN=$1
else
	echo "Usage: $0 oldname newname subdir"
	exit 1
fi

if test -n "$2"
then
	OUT=$2
else
	echo "Usage: $0 oldname newname subdir"
	exit 1
fi

if test -n "$3"
then
	SYSTEM=$3
	GNUPLOT=gnuplot/$SYSTEM
	mkdir -p $FIG/$SYSTEM $FIG/$GNUPLOT $FIG/$GNUPLOT/test 2>/dev/null || true
else
	echo "Usage: $0 oldname newname subdir"
	exit 1
fi

if test -f $IN.out.ps
then
	mv test/$IN.out $FIG/$GNUPLOT/test/
	sed -e s=$IN=$OUT=g $IN.out.ps		> $FIG/$GNUPLOT/$OUT.eps
	test -f $FIG/$SYSTEM/$OUT.eps && rm $FIG/$SYSTEM/$OUT.eps
	ln $FIG/$GNUPLOT/$OUT.eps $FIG/$SYSTEM/$OUT.eps

# To save disk space, omit *.[123] files, since they can be recreated
# from the test/*.out file:
#	sed -e s=$IN=$OUT=g $IN.out.1		> $FIG/$GNUPLOT/$OUT.out.1
#	sed -e s=$IN=$OUT=g $IN.out.2		> $FIG/$GNUPLOT/$OUT.out.2
#	sed -e s=$IN=$OUT=g $IN.out.3		> $FIG/$GNUPLOT/$OUT.out.3
	sed -e s=$IN=$OUT=g $IN.out.4		> $FIG/$GNUPLOT/$OUT.out.4
	sed -e s=$IN=$OUT=g $IN.out.5		> $FIG/$GNUPLOT/$OUT.out.5
	sed -e s=$IN=$OUT=g $IN.out.gnu		> $FIG/$GNUPLOT/$OUT.out.gnu
	sed -e s=$IN=$OUT=g $IN.out.map		> $FIG/$GNUPLOT/$OUT.out.map
	sed -e s=$IN=$OUT=g $IN.out.map.out	> $FIG/$GNUPLOT/$OUT.out.map.out
	rm -f $IN.out.*

	ls -lo $FIG/$GNUPLOT/test/$IN.out $FIG/$SYSTEM/$OUT.eps $FIG/$GNUPLOT/$OUT.out.*
else
	echo cannot find $IN.out.ps
	exit 1
fi
