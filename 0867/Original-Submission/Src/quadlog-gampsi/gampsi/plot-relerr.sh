#! /bin/sh
#=======================================================================
# Given a gam*.out or psi*.out filename on the command line, read the
# test/FILENAME and test/okay/FILENAME files, and split them into
# three separate files for later input into Maple for high-precision
# relative error computation, and then into gnuplot for plotting.
#
# The output files are named with the same name, but with extensions
# .1, .2, .3 and .4, for the x, fapprox(x), fexact(x), and relerr(x)
# values respectively. The Maple program for computing the relative
# error values is stored in a file with extension .map, and the
# gnuplot program for plotting the relative error is in a file with
# extension .gnu.
#
# Then run maple to generate the .4 file, then then gnuplot to produce
# both an X11 display of the results, and a corresponding PostScript
# .ps file for printing and typesetting.
#
# Usage:
#	./plot-relerr.sh file functionname \
#		[logx? [xmin:xmax [ymin:ymax [xformat [yformat [abserr [machid]]]]]]]
#
# As the script runs, it echoes the names of each output file as it is
# completed.
#
# If the third argument is omitted (or empty), the horizontal axis is
# x; otherwise, it is log10(x).  The particular value of that argument
# is irrelevant.
#
# If the fourth argument is omitted, the x range is determined
# automatically from the input data.  Otherwise it should be a
# colon-separated pair of numbers defining the x subrange to be
# selected for plotting.
#
# If the fifth argument is omitted, the y range is determined
# automatically from the input data.  Otherwise it should be a
# colon-separated pair of numbers defining the y subrange to be
# selected for plotting.
#
# If the sixth argument is omitted, the x axis is numbered with the
# default %g format, which trims trailing zeros, producing uneven
# label widths.  This argument allows specification of an alternate
# format, such as %.2f.
#
# If the seventh argument is omitted, the y axis is numbered with the
# default %g format, which trims trailing zeros, producing uneven
# label widths.  This argument allows specification of an alternate
# format, such as %.2f.
#
# If the eighth argument is omitted, then no absolute error plots
# are produced. Otherwise, it should be an arbitrary nonempty string,
# and absolute error plots will be made.
#
# If the ninth argument is provided, it is used instead of the return
# value of ../common/machid.sh; this permits processing of data from
# other systems that lack Maple and/or gnuplot.
#
# Examples:
#	./plot-relerr.sh gam01.out  dgamma
#	./plot-relerr.sh gam01.out  dgamma
#	./plot-relerr.sh psi01.out  dpsi  logx
#	./plot-relerr.sh psi01.out  dpsi  1
#	./plot-relerr.sh psi201.out dpsi
#	./plot-relerr.sh psi201.out qpsi  ''   1.3:1.5  -40:-30  %.2f  %.1f
#
# [31-May-2000]
#=======================================================================

LIBDIR=`dirname $0`

if test $# -lt 1
then
	echo Usage: $0 'file functionname [logx? [xmin:xmax [ymin:ymax [xformat [yformat [abserr]]]]]]'
	exit 1
fi

if test -f test/$1 -a -f test/okay/$1
then
	true
else
	echo 'File(s)' test/$1 and/or test/okay/$1 do not exist
	exit 1
fi


rm -f $1.[12345] $1.gnu $1.map $1.map.out $1.ps

#=======================================================================
# Quadruple-precision data cannot be filtered to prepare plots of
# relative and absolute error without a working makdat; if its build
# fails, issue a warning, and continue using sort-out.awk to convert
# (f,p) pairs to f*2**p values that are written to the $1.1 file.  The
# remaining numerical data is handled exclusively as character strings.
# makdat builds will fail on GNU/Linux systems, where the native
# compilers lack quadruple-precision support.
#=======================================================================
# rm -f makdat
if make makdat >/dev/null 2>&1
then
	# Sun f95 incorrectly outputs D exponent instead of E exponent
	FIXEXP="sed -e 's/D-/e-/g' -e 's/D[+]/e+/g'"
	FIXEXP="sed -f fixexp.sed"
	awk '{print $1, $2, $3, "\047" $5 "\047"}' < test/$1 | \
		 ./makdat | $FIXEXP | \
			 awk '{print $1}' >$1.1
	echo $1.1
	awk '{print $1, $2, $3, "\047" $5 "\047"}' < test/$1 | \
		 ./makdat | $FIXEXP | \
			 awk '{print $2}' >$1.2
	echo $1.2
	awk '{print $1, $2, $3, "\047" $5 "\047"}' < test/okay/$1 | \
		 ./makdat | $FIXEXP | \
			 awk '{print $2}' >$1.3
	echo $1.3
else
	echo 'WARNING: makdat could not be built: quadruple-precision plot data will be worthless!'
	awk -f sort-out.awk test/$1 | awk '{print $1}' >$1.1
	echo $1.1
	awk -f sort-out.awk test/$1 | awk '{print $2}' >$1.2
	echo $1.2
	awk -f sort-out.awk test/okay/$1 | awk '{print $2}' >$1.3
	echo $1.3
fi

# Check for logx? argument:
if test -n "$3"
then
	xlabel="log10(x)"
	x='log10($1)'
else
	xlabel="x"
	x='$1'
fi

# Check for xrange argument:
if test -n "$4"
then
	xrange="[$4]"
else
 	xrange='[*:*]'
fi

# Check for yrange argument:
if test -n "$5"
then
	yrange="[$5]"
else
 	yrange='[*:*]'
fi

# Check for xformat argument:
if test -n "$6"
then
	xformat="$6"
else
 	xformat="%g"
fi

# Check for yformat argument:
if test -n "$7"
then
	yformat="$7"
else
 	yformat="%g"
fi

# Check for abserr argument:
if test -n "$8"
then
	abserr="$8"
else
 	abserr=""
fi


# Check for machid argument:
if test -n "$9"
then
	title="$9"
else
	title=`../common/machid.sh`
fi
title="`echo $title | sed -e 's/-/ /g' `"
compiler=`awk '/^ *FC[ 	]*=/ {print $3}' Makefile`
title=`echo $title [$compiler]`

precision=`echo $2 | cut -c 1`
case $precision in
	d) ulp='2.0**(-52)' ;;
	q) os=`uname -s`
	   case $os in
	        AIX|IRIX64)
			ulp='2.0**(-105)'
			;;
		*)
			ulp='2.0**(-112)'
			;;
	   esac
	   ;;
	*) ulp='2.0**(-23)' ;;
esac

#-----------------------------------------------------------------------

cat <<EOF >$1.map
########################################################################
## WARNING: Do NOT edit this file.  It was created automatically
## by
##	$0 $1 $2 "$3" "$4" "$5" "$6" "$7"
## for
##	${USER}@`hostname`
## on
##	`date`
########################################################################

atof := proc(str)
		local r:
		r := sscanf(str,"%e"):
		r[1]
	end:

sizeof := proc(x)
		local e,n:
		n := 0:
		for e in x
		do
			n := n + 1
		od:
		n
	end:
Digits := 50:
tx := readdata("$1.1", string):
tfapprox := readdata("$1.2", string):
tfexact := readdata("$1.3", string):
nx := sizeof(tx):

## Sanity check on the input vectors:
if (sizeof(tfapprox) <> sizeof(tfexact)) then quit fi:
if (nx <> sizeof(tfexact)) then quit fi:

## To work around stupid limitations of Maple, we have to
## convert the temporary t*[] arrays to predeclared arrays:
x       := array(1..nx);
fapprox := array(1..nx);
fexact  := array(1..nx);

## Convert data from strings to numbers or symbols; readdata() without
## the second string argument simply drops nonnumeric data, like "NaN"
## and "Infinity", but we need to preserve them as placeholders.
for k from 1 to nx
do
	x[k] := parse(tx[k]):
	fapprox[k] := parse(tfapprox[k]):
	fexact[k] := parse(tfexact[k])
od:

## relative error file: the type checks avoid output of Infinity and NaN,
## and x=0 is eliminated to avoid problems with log(x)
for k from 1 to nx
do
	if (type(x[k],numeric) and
	    type(fapprox[k],numeric) and
	    type(fexact[k],numeric) and
	    (evalb(x[k] <> 0)))
	then
		fprintf("$1.4",
			"%.15e\t%.15e\n",
			x[k],
			(fapprox[k] - fexact[k])/fexact[k])
	fi
od:
close("$1.4"):

## absolute error file: the type checks avoid output of Infinity and NaN,
## and x=0 is eliminated to avoid problems with log(x)
for k from 1 to nx
do
	if (type(x[k],numeric) and
	    type(fapprox[k],numeric) and
	    type(fexact[k],numeric) and
	    (evalb(x[k] <> 0)))
	then
		fprintf("$1.5",
			"%.15e\t%.15e\n",
			x[k],
			(fapprox[k] - fexact[k]))
	fi
od:
close("$1.5"):
EOF
echo $1.map

#-----------------------------------------------------------------------
cat <<EOF >$1.gnu
########################################################################
## WARNING: Do NOT edit this file.  It was created automatically
## by
##	$0 $1 $2 "$3" "$4" "$5" "$6" "$7" "$8"
## for
##	${USER}@`hostname`
## on
##	`date`
########################################################################

set title "$title"
set xlabel "$xlabel"
set ylabel "log10(relative error in $2(x))"
set xrange $xrange
set yrange $yrange
set format x "$xformat"
set format y "$yformat"

plot "$1.4" using ($x):(log10(abs(\$2))) notitle, \\
     "$1.4" using ($x):(log10($ulp)) notitle with lines linewidth 5

set term postscript eps "Helvetica" 32
set output "$1.ps"
set size 1.5, 1.5
plot "$1.4" using ($x):(log10(abs(\$2))) notitle, \\
     "$1.4" using ($x):(log10($ulp)) notitle with lines linewidth 15
pause -1 "Hit return to quit"
EOF

if test -n "$abserr"
then
	cat <<EOF >>$1.gnu
set term x11
set size 1, 1
set ylabel "log10(absolute error in $2(x))"

plot "$1.5" using ($x):(log10(abs(\$2))) notitle, \\
     "$1.5" using ($x):(log10($ulp)) notitle with lines linewidth 5

set term postscript eps "Helvetica" 32
set output "$1.abs.ps"
set size 1.5, 1.5
plot "$1.5" using ($x):(log10(abs(\$2))) notitle, \\
     "$1.5" using ($x):(log10($ulp)) notitle with lines linewidth 15
pause -1 "Hit return to quit"

EOF
fi
echo $1.gnu

mapleV5.1 < $1.map > $1.map.out
echo $1.4
echo $1.5

gnuplot $1.gnu

# Optimize the output PostScript files to avoid 1000x slowdown in
# printing on Lexmark printers!
$LIBDIR/fix-gnuplot.sh $1.ps

if test -n "$abserr"
then
	$LIBDIR/fix-gnuplot.sh $1.abs.ps
fi


# The default of butt caps (0 setlinecap) in PostScript leaves ugly
# gaps in the horizontal line marking the ULP value, so patch the
# PostScript file to fix it.  The resulting line may stick out of the
# plot frame by half the linewidth, if it reaches the frame edges.
mv $1.ps $1.ps.tmp
sed -e 's/^LT1$/LT1 2 setlinecap/' < $1.ps.tmp > $1.ps
rm $1.ps.tmp
echo $1.ps

if test -n "$abserr"
then
	mv $1.abs.ps $1.abs.ps.tmp
	sed -e 's/^LT1$/LT1 2 setlinecap/' < $1.abs.ps.tmp > $1.abs.ps
	rm $1.abs.ps.tmp
	echo $1.abs.ps
fi
