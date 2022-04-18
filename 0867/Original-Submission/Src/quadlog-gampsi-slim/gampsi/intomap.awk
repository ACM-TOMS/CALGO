# /u/sy/beebe/jim-ball/toms/dist/gampsi/intomap.awk, Tue Mar 28 14:36:08 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Convert a gam*.out or psi*.out file to a suitable input file for
# use in Maple V.
#
# Usage:
#	awk -f intomap.awk infile >mapfile
#       mapleV5 <mapfile >outfile
#
# [03-Aug-2000] -- update for ln(Gamma(x)) and psiln(x) support
# [14-Jul-2000] -- update for cos() support
# [01-Jun-2000] -- update for exp() and log() support
# [28-Mar-2000]
# ======================================================================

BEGIN	{
		n = 0
		fcn = "UnknownFunction"
	}

/^#/	{
		print
		if (index($0,"of cos(x)") > 0)
		    fcn = "cos"
		else if (index($0,"of gamma(x)") > 0)
		    fcn = "GAMMA"
		else if (index($0,"of Gamma(x)") > 0)
		    fcn = "GAMMA"
		else if (index($0,"of ln(Gamma(x))") > 0)
		    fcn = "lnGamma"
		else if (index($0,"of exp(x)") > 0)
		    fcn = "exp"
		else if (index($0,"of log(x)") > 0)
		    fcn = "log"
		else if (index($0,"of psi(x)") > 0)
		    fcn = "Psi"
		else if (index($0,"of psiln(x)") > 0)
		    fcn = "psiln"
		else if (index($0,"of sin(x)") > 0)
		    fcn = "sin"
		else if (index($0,"of tan(x)") > 0)
		    fcn = "tan"
	}

$1 == 2	{
		n++
		f[n] = $2
		p[n] = $3
	}

END	{
		print "interface(quiet = true):"
		if (fcn == "psiln")
		{
		    print "psiln := proc(x_)"
		    print "    local result:"
		    print "    # We need a precision increment of at least log10(xmax), where xmax"
		    print "    # is the largest argument required for the tests."
		    print "    Digits := Digits + 310:"
		    print "    result := evalf(Psi(x_) - ln(x_)):"
		    print "    Digits := Digits - 310:"
		    print "    RETURN (result):"
		    print "end:"
		}
		else if (fcn == "lnGamma")
		{
		    print "lnGamma := proc(x_)"
		    print "    local result:"
		    print "    # We need a precision increment of at least log10(xmax), where xmax"
		    print "    # is the largest argument required for the tests."
		    print "    Digits := Digits + 310:"
		    print "    result := evalf(lnGAMMA(abs(x_))):"
		    print "    Digits := Digits - 310:"
		    print "    RETURN (result):"
		    print "end:"
		}
		print "Digits := 50 :"
		print "n :=", n, ":"
		print "f := array(1.." n ") :"
		print "p := array(1.." n ") :"
		print ""
		print "f := ["
		for (i = 1; i <= n; ++i)
			print f[i] ((i < n) ? "," : "")
		print "] :"
		print ""
		print "p := ["
		for (i = 1; i <= n; ++i)
			print p[i] ((i < n) ? "," : "")
		print "] :"
		print ""
		print "for i from 1 by 1 to n do"
		print "\tprintf(\"%1d %20.0f. %5d %1d %57.49e %1d %1d\\n\","
		print "\t\t2, f[i], p[i], 0, " fcn "(f[i]*2^p[i]), 0, 0)"
		print "od :"
	}
