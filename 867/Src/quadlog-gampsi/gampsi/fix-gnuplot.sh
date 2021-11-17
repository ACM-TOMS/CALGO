#! /bin/sh
########################################################################
# gnuplot produces excessively long paths for
#
# 	plot "qpsi105.out.4" using ($1):(log10(2.0**(-112))) notitle \
#		with lines linewidth 5
#
# which cause Lexmark printers to take up to 30 min to print.  When
# these paths are eliminated, printing takes only a few seconds.
#
# The optimization is therefore to look for sequences of commands
#	LT1 2 setlinecap
#	2088 3377 M
#	182 0 V
#	33 0 V
#	...
#	235 0 V
#	currentpoint stroke M
# with 399 V commands, and replace the V sequence by a single command.
#
# Usage:
#	./fix-gnuplot.sh eps-file(s)
#
# The input files are renamed with suffix .pre-opt, and the output
# optimized files have the original names.
#
# [04-Jul-2000]
########################################################################

fix_gnuplot()
{
    awk '
	BEGIN	{
			cpx = -2147483647
			max_cpx = cpx
			min_cpx = -cpx
			nv = 0
		}

	/^%%EndComments/ \
		{
			print 
			print "% Fixed by fix-gnuplot.sh"
			next
		}

	(NF == 3) && ($1 ~ /^[-+]?[0-9]+$/) && ($2 ~ /^[-+]?[0-9]+$/) && ($3 == "M") \
		{
			cpx = $1
			cpy = $2
			max_cpx = cpx
			min_cpx = cpx
			moveto_cpx = cpx
			moveto_cpy = cpy
		}

	(NF == 3) && ($1 ~ /^[-+]?[0-9]+$/) && ($2 == "0") && ($3 == "V") \
		{
			cpx += $1
			max_cpx = max(cpx,max_cpx)
			min_cpx = min(cpx,min_cpx)
			nv++
			Vlast = $0
			next
		}

	(nv > 0) \
		{
			# End of V sequence: replace by single line
			if (nv > 1)
			{
				if ((min_cpx != moveto_cpx) || (cpy != moveto_cpy))
					print min_cpx, cpy, "M"
				print max_cpx, cpy, "L"
			}
			else
				print Vlast
			nv = 0
			min_cpx = cpx
			max_cpx = cpx
			moveto_cpx = cpx
			moveto_cpy = cpy
		}

		{ print }

	function max(a,b)
	{
	    a += 0
	    b += 0
	    return ((a > b) ? a : b)
	}

	function min(a,b)
	{
	    a += 0
	    b += 0
	    return ((a < b) ? a : b)
	}
	' $1
}


for f in "$@"
do
	fix_gnuplot $f >$f.tmp.$$ && \
		mv -f $f $f.pre-opt && \
		mv -f $f.tmp.$$ $f
	ls -lo $f.pre-opt $f
done

