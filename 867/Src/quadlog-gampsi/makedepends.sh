#! /bin/sh
# Search the *.f files for INCLUDE statements, and output a set of
# dependencies for insertion in a Makefile.in file.
#
# Usage:
#	./makedepends.sh
#
# [01-May-2000]

grep '^       *INCLUDE ' *.f | \
	awk '{printf("%-15s\t$(COMMON)/%s\n", 
		    (substr($1,1,length($1)-2) "o:"),
		    (substr($3,2,length($3)-2)))}' | \
		sort | \
			awk '{if (LAST != $1) print ""
			      LAST = $1
			      print	      
			      }' | \
				unexpand

