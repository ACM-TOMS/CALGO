# /u/sy/beebe/jim-ball/toms/dist/gampsi-new/gendat.awk, Tue Jul 11 05:19:57 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Take values of xmin and xmax from the command line, and output a
# candidate test file.
#
# Usage:
#	awk -f gendat.awk -v xmin=xxx -v xmax=xxx -v fun=xxx
#
# [11-Jul-2000]
# ======================================================================

BEGIN \
{
    pmax = 31
    pmin = pmax - 1
    maxint = 2^pmax - 1
    minint = 2^pmin
    if (xmin > 0)
    {
	qmin = floor(log2(abs(xmin))) - pmin
	qmax = ceil(log2(abs(xmax))) - pmax
    }
    else
    {
	qmin = ceil(log2(abs(xmin))) - pmax
	qmax = floor(log2(abs(xmax))) - pmin
    }
    sign = (xmin < 0) ? -1 : 1
    print (fun == "") ? "xxx" : fun
    print ceil(abs(10000 / (abs(qmax - qmin) + 1)))
    print sign
    if (qmin < qmax)
	print qmin, qmax
    else
	print qmax, qmin
    print minint, maxint
    if (xmin > 0)
	printf("%.5e %.5e\n", xmin, xmax)
    else
	printf("%.5e %.5e\n", -xmax, -xmin)
    print 0, 0
    print "************************************************************************"
    if (xmin > 0)
    {
	print "*   xrange = " sign*minint " * 2^(" qmin ") .. " \
	    sign*maxint " * 2^(" qmax ")"
	print "*          = " sign*minint*2^qmin " .. " \
	    sign*maxint*2^qmax " [calculated]"
    }
    else
    {
	print "*   xrange = " sign*maxint " * 2^(" qmin ") .. " \
	    sign*minint " * 2^(" qmax ")"
	print "*          = " sign*maxint*2^qmin " .. " \
	    sign*minint*2^qmax " [calculated]"
    }
    print "*          = " xmin " .. " xmax " [requested]"
    print "*   minint = 2^(" pmax - 1 ")     = " minint
    print "*   maxint = 2^(" pmax ") - 1 = " maxint
    print "************************************************************************"
    
    exit(0);    
}

function abs(x)
{
    x = x + 0
    return ((x < 0) ? -x : x)
}



function ceil(x, r)
{
    r = int(x)
    if (r < x)
	r += 1
    return (r)
}


function floor(x, r)
{
    r = int(x)
    if (r > x)
	r -= 1
    return (r)
}


function log2(x)
{
    return (log(x)/log(2.0))
}
