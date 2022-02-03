# /u/sy/beebe/jim-ball/toms/dist/gampsi/sort-out.awk, Tue May 30 14:58:28 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Filter an gam*.out or psi*.out file, producing a file sorted by
# ascending x value, with pairs (x, f(x)) on each line.
#
# Usage:
#	awk -f sort-out.awk infile >outfile
#
# [30-May-2000]
# ======================================================================

BEGIN	{
	    sort_pipe = "/usr/local/bin/sort -k 1g" 
	    sort_pipe = "/usr/local/bin/sort +0g -1" 
	}

/^#/	{ next }

	{ x = $2 * 2.0^($3); printf("%.15e\t%s\n", x, tolower($5)) | sort_pipe }

END	{ close(sort_pipe) }
