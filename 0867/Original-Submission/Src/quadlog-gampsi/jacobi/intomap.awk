# /u/sy/beebe/jim-ball/toms/dist/jacobi/intomap.awk, Tue May  9 18:18:06 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Convert jacobi/test/*.in files to Maple files for creating 
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
    print "alpha := ", alpha ":" >> OutFile
    print "beta := ", beta ":" >> OutFile
    print "n := ", n ":" >> OutFile
    print "nquad := ", nquad ":" >> OutFile
    print "for k from 0 to n" >> OutFile
    print "do" >> OutFile
    print "\tprintf(\"%d %4d\\t%15.8e\\t%15.8e\\t%d\\t%50.35e 0 0\\n\"," >> OutFile
    if (FILENAME ~ "(gjf1|gjfd1)")
	print "\t\t3, k, alpha, beta, nquad, evalf(int((1+x)^k * (1-x)^alpha * (1+x)^beta * ln(1+x), x = -1..+1))):" >> OutFile
    else
	print "\t\t3, k, alpha, beta, nquad, evalf(int((1+x)^k * (1-x)^alpha * (1+x)^beta, x = -1..+1))):" >> OutFile
    print "od;" >> OutFile
}

function new_file()
{
    if (OutFile != "")
	close (OutFile)
    OutFile = FILENAME
    sub("[.][a-z]+",".map",OutFile)     
    DatFile = FILENAME
    sub("[.][a-z]+",".dat",DatFile)     
    system("grep '^#' < okay/" DatFile " > " OutFile)
    print "interface(quiet = true):" >> OutFile
    print "Digits := 50:" >> OutFile
}
