#! /bin/sh
###=====================================================================
### Search the PATH for the executable(s) given on the command line, and
### if found, display their full path on stdout, one per line, and exit
### with a success code (0), if at least one of them was found.
### Otherwise, display nothing, and exit with a failure code (1).
###
### Usage:
###	which.sh prog1 prog2 ... progn
###
### The reason for using this private script, rather than the system
### `which' command, is that some implementations of the latter try to
### be clever, and read user startup files, producing unpredictable
### results if a set-window-title string is produced as part of the
### output.  Experiments also show that the exit codes returned by
### system `which' commands are unreliable, and the output is too, since
### on some systems, it contains "no foobar in /a /b /c ...", sigh...
###
### [29-Jan-2000]
###=====================================================================

retcode=1

for program in "$@"
do
    IFS=: pathlist="$PATH"
    for path in $pathlist
    do
	## Some old implementations of test don't support the -x option,
	## so fall back to the -f option if -x fails.  Executable files
	## need not be readable, so -r is NOT the option to use!
	if test -x $path/$program 2>/dev/null || test -f $path/$program 2>/dev/null
	then
		echo $path/$program
		retcode=0
		## we must find only the FIRST occurrence in the path
		break
	fi
    done 2>/dev/null
done

exit $retcode
