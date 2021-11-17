###=====================================================================
### Process non-nested Fortran INCLUDE statements, replacing them by
### their contents, assumed to come from the directory ../common.
###
### Usage:
###	awk -f doincl.awk infile >outfile
###
### [14-Dec-2000]
###=====================================================================

/^ +[Ii][Nn][Cc][Ll][Uu][Dd][Ee] *'/ \
	{
		gsub(/'/,"",$2)
		print FILENAME ": " $2 > "/dev/stderr"
		incfile = ("../common/" $2)
		while ((getline x < incfile) > 0)
			print x
		close incfile
		next
	}
	{ print }
