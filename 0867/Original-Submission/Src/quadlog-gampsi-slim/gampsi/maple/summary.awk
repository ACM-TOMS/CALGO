# /u/sy/beebe/jim-ball/toms/dist/gampsi/maple/summary.awk, Fri Jul 28 06:24:21 2000
# Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
# ======================================================================
# Extract and print a summary of the Pad{\'e} results from one or more
# log files.
#
# Usage:
#	awk -f summary.awk logfile(s) >summary-file
#
# [28-Jul-2000]
# ======================================================================

# NB: This requires GNU sort's -g option for sorting floating-point numbers
BEGIN				{ sort_pipe = "/usr/local/bin/sort +7g -8 +21n -22 +0 -1" }

/^[*] *Maximum *relative *err/	{
					Line = $0
					next
				}
/^[*] *at *x *=/		{
					Line = Line " " substr($0,2)
					sub("x - x0","x-x0",Line)
					gsub(" +"," ",Line)
					LineNumber = FNR
					Filename = FILENAME
					next
				}

/^Attempting/			{
					report()
					next
				}

/^ *X *-/			{
					negative_signs++
					next
				}


END				{
					report()
					close(sort_pipe)
				}

#=======================================================================

function report( balanced,parts,positive,precision)
{
    if (Line != "")
    {
	if (match(Line,"degree \\[[0-9]+,[0-9]+\\]"))
	{
	    split(substr(Line,RSTART+8,RLENGTH-9),parts,",")
	    positive = ((negative_signs == 0) ? " POSITIVE!" : "")
	    Line = Line "\t= " (parts[1] + parts[2]) 
	    balanced = (parts[1] == parts[2]) ? "\tBALANCED!" : "\t\t "
	}
	else
	{
	    balanced = ""
	    positive = ""
	}
	split(Line,parts," ")
	if (parts[8] <= (2.0^(-112)))
	    precision = "128-bit"
	else if (parts[8] <= (2.0^(-63)))
	    precision = " 80-bit"
	else if (parts[8] <= (2.0^(-52)))
	    precision = " 64-bit"
	else if (parts[8] <= (2.0^(-23)))
	    precision = " 32-bit"
	else
	    precision = "USELESS"
	print Filename ":" sprintf("%05d",LineNumber) ":" Line "\t" precision balanced positive | sort_pipe
    }
    Line = ""
    negative_signs = 0
}
