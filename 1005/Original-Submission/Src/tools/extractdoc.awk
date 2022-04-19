# EXTRACTDOC.AWK: Extract documentation from a Fortran source file.
#
# The first comment block is extracted with the leading ! removed.
# The block begins with the first comment line and ends with the
# first non-comment line of the file.
#
# Usage: awk -f extractdoc.awk file.f90 > file.txt
#
# Compatibility: Has been tested with gawk, mawk and Apple awk.

NR==1 {
   print
   print ""
}

/^[ \t]*!/ {printon = 1}

printon    {
   if ($0 !~ /^[ \t]*!.*/) exit
   sub("^[ \t]*! ?", "")
   print
}
