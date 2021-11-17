#! /bin/sh
# Run intomap on one or more test/*.out files to produce *.map files,
# and then run maple on them, in PARALLEL, filtering out the unwanted
# chatter from maple, leaving the output in files *.tmp in the current
# directory.
#
# Usage:
#	./intomap.sh test/*.out
#
# [29-Jun-2000] -- add echo command to show how to rename output files
# [04-Jun-2000]

LIBDIR=`dirname $0`

for f in "$@"
do
	g=`basename $f .out`
	echo $g
	awk -f $LIBDIR/intomap.awk < $f > $g.map
	(mapleV5.1 < $g.map | grep '^ *[#0-9]' > $g.tmp; echo Maybe mv $g.tmp test/okay/$g.out) &
done
wait
