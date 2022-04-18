# /u/sy/beebe/jim-ball/toms/dist/jacobi/intomap.awk, Tue May  9 18:18:06 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Convert jacobi/test/*.in files to Mathematica files for creating
# high-precision *.dat files.
#
# Usage:
#	awk -f intomap.awk *.in
#
# [09-May-2000]
# ======================================================================

BEGIN {}

LastFile != FILENAME	{ new_file() }

			{ convert($1,$2,$3,$4) }

END   {
    if (OutFile != "")
	close (OutFile)
}

function convert(nquad,alpha,beta,n)
{
    print "$MinPrecision = 50;" >> OutFile
    print "alpha = ", alpha "`50;" >> OutFile
    print "beta = ", beta "`50;" >> OutFile
    print "one = 1`50;" >> OutFile
    print "n = ", n ";" >> OutFile
#    print "fp[x_] := NumberForm[x,NumberFormat -> (SequenceForm[#1,\"e\", #3]&)];" >> OutFile
    print "fp[x_] := FortranForm[x];" >> OutFile
    print "nquad = ", nquad ";" >> OutFile
    print "For \\\n[\n\tk = 0, k < n, k++," >> OutFile
    print "\tPrint \\\n\t[\n\t\tk, Tab," >> OutFile
    print "\t\talpha, Tab," >> OutFile
    print "\t\tbeta, Tab," >> OutFile
    print "\t\tnquad, Tab," >> OutFile
    if (FILENAME ~ "(gjf1|gjfd1)")
	print "\t\tfp[N[Integrate[(one + x)^k *(one - x)^alpha * (one + x)^beta * Log[one+x], {x, -one, one}]]], Tab," >> OutFile
    else
	print "\t\tfp[N[Integrate[(one + x)^k *(one - x)^alpha * (one + x)^beta, {x, -one, one}]]], Tab," >> OutFile
    print "\t\t0, Tab," >> OutFile
    print "\t\t0\n\t];" >> OutFile
    print "];" >> OutFile
}

function new_file()
{
    if (OutFile != "")
	close (OutFile)
    OutFile = FILENAME
    LastFile = OutFile
    sub("[.][a-z]+",".math",OutFile)
    DatFile = FILENAME
    sub("[.][a-z]+",".dat",DatFile)
    system("grep '^#' < okay/" DatFile " | awk '{print \"(*\", $0, \"*)\"}' > " OutFile)
}
