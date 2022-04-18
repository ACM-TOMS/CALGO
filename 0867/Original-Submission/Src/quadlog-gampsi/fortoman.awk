### /u/sy/beebe/toms-test/dist/fortoman.awk, Mon Dec 18 15:03:29 2000
### Edit by Nelson H. F. Beebe <beebe@math.utah.edu>
### ====================================================================
### Filter a Fortran file with standardized documentation header (PLOT79
### style) and standardized declarations (from ftnchek), and create a
### UNIX manual page file on stdout to document the file.
###
### Usage:
###	awk -f fortoman.awk \
###		[-v AUTHORFILE=filename] \
###		[-v INCLUDEFILES="file1 file2 ..."] \
###		[-v MANSECT=manual-section ] \
###		[-v SEEALSO=filename ] \
###		Fortran-file >man-page-file
###
### The command-line variable setting of AUTHORFILE can choose an
### alternate file with author information.  The default is
### "AUTHORFILE.std".
###
### The command-line variable setting of INCLUDEFILES specifies one or
### more C header files from which suitable function prototypes can be
### extracted for inclusion in the output manual pages.  If prototypes
### are found in multiple files, the FIRST one seen is used: thus,
### order those files accordingly if it matters.
###
### The default manual section is 3 (runtime library), unless
### overridden by a command-line setting of MANSECT.
###
### The SEEALSO file contains lines with multiple words per line, for
### example,
###
### 	f77 f90 f95 hpf nagf90 nagf95 pgf77 pgf90
###	c89 CC cc cxx g++ gcc lcc pgCC pgcc
###	...
###
### The words are the basenames of related manual page files, and SEE
### ALSO entries will be generated for each of the remaining ones.
### They will be output in the order given, so the words should be
### sorted within each line.
###
### [22-Dec-2000] -- split function terminate() into separate
###                  print_xxx() functions, supply missing Declaration
###                  output, reformat the SYNOPSIS section slightly,
###                  fix .TH line output, and add support for MANSECT.
### [21-Dec-2000] -- add SEEALSO support
### [18-Dec-2000] -- original version
### ====================================================================

BEGIN	{ initialize() }

FNR == 1 {
	    Prototype = $0
	    while (($0 ~ /, *$/) && (getline > 0))
		Prototype = Prototype " " substr($0,7)
	    sub("^ *","",Prototype)
	    gsub(" +"," ", Prototype)
	 }

/^[*]     [[(].*[])] *$/ && (Summary == "") \
	{
	    Summary = substr($0,7)
	    sub("^ *[[(]","",Summary)
	    sub("[])] *$","",Summary)
	    next
	}

/^[*]     [[(][0-9]+-[A-Za-z]+-[12][0-9][0-9][0-9][])] *$/ \
	{
	    Last_Modified_Date = $0
	    sub("^[*]     [[(]","",Last_Modified_Date)
	    sub("[])] *$","",Last_Modified_Date)
	    split(Last_Modified_Date,lmd_parts,"-")
	    # Reset DMY whenever the code has its own date
	    DMY = sprintf("%02d",lmd_parts[1]) " " Month_Expansion[tolower(lmd_parts[2])] " " lmd_parts[3]
	    next
	}

/^[*][*]+$/ { StarLines++; next }

/^[*] *$/ && (StarLines == 1) \
	{
	    if (Indented_Lines)
	    {
		Description = Description ((Description == "") ? "" : "\n") ".fi\n.RE\n.PP"
		Indented_Lines = 0
	    }
	    else
	    {
		sub("^[*] *$",".PP") # convert empty comment lines to paragraph breaks
		Description = Description ((Description == "") ? "" : "\n") troff_protect($0)
	    }
	    next
	}


/^[*] /	&& (StarLines == 1) \
	{
	    sub("^[*]     ","")	       # handle nonempty comment lines
	    if ($0 ~ /^ /)
	    {
		Indented_Lines++
		if (Indented_Lines == 1)
		    Description = Description ((Description == "") ? "" : "\n") ".RS\n.nf"
		sub("^ *","")
	    }
	    else
	    {
		if (Indented_Lines)
		    Description = Description ((Description == "") ? "" : "\n") ".fi\n.RE\n.PP"
		Indented_Lines = 0
	    }

	    Description = Description ((Description == "") ? "" : "\n") troff_protect($0)
	    next
	}

/^[*] *Argument variables/, /^[*] *Local variables/ \
	{
	    if ($0 ~ /^ /)
	    {
		## Jump to END rule if we reach executable code (in case
		## there are no local variables declared) or INCLUDEd 
		## text.
		if ($0 ~ /^ *[Ii][Nn][Cc][Ll][Uu][Dd][Ee] *'/)
		    exit(0)
		else if ($0 ~ /^ *[Ii][Ff] *[(]/)
		    exit(0)
		else if ($0 ~ /^ +[Cc][Aa][Ll][Ll] *[A-Za-z]/)
		    exit(0)
		else if ($0 ~ /^ +.*=/)
		    exit(0)
		sub(/^ +/,"")
		Declarations = Declarations ((Declarations == "") ? "" : "\n") bold($0)
	    }
	}


END   { terminate() }

### ====================================================================

function bold(s)
{
    return (".B \"" s "\"")
}


function bold_inline(s)
{
    return ("\\fB" s "\\fP\\&")
}


function bold_roman(s,t)
{
    return (".BR \"" s "\" \"" t "\"")
}


function copy_file(infile, s)
{
    while ((getline s < infile) > 0)
	print s
    close(infile)
}


function do_include_files(filelist, k,parts,n)
{
    n = split(filelist,parts," ")
    for (k = 1; k <= n; ++k)
	do_one_include_file(parts[k])
}


function do_one_include_file(infile, cmd,routine_name,line,multiline_comment,n,parts,prototype)
{
    routine_name = ""
    multiline_comment = 0
    prototype = ""

    FILENAME = infile
    FNR = 0
    cmd = "expand < " infile
    n = 0
    while ((cmd | getline line) > 0)
    {
	FNR++
	if (match(line,"^ *#"))		# ignore preprocessor lines
	    continue
	else if (match(line,"^ *$"))	# ignore blank lines
	    continue
	else if (match(line,"^ */[*]"))	# ignore comment lines
	{
	    multiline_comment = !match(line,"[*]/")
	    continue
	}
	else if (multiline_comment)
	{
	    if (match(line,"[*]/"))
		multiline_comment = 0
	}
	else if (match(line,"^[ a-zA-Z0-9_*]+[(]") && !multiline_comment)
	{
	    ## warning("DEBUG: proto head = [" substr(line,RSTART,RLENGTH-1) "]")
	    n = split(substr(line,RSTART,RLENGTH-1),parts," ")
	    routine_name = trim(parts[n])
	    prototype = trim(line)
	    sub("^ *extern +","",prototype)
	    n = index(prototype,"(")
	    prototype = bold(prototype)
	}
	else if (prototype != "")
	{
	    if (n == 0)
	    {
		warning("missing prototype header line")
		n = 4
	    }
	    prototype = prototype "\n" bold(indent(n) trim(line))
	}
	if (match(line,"[)] *;"))
	{
	    if (routine_name == "")
		warning("no function name found in prototype")
	    else if (routine_name in Prototype_Table)
	    {
		if (Prototype_Table[routine_name] == prototype)
		    warning("duplicate prototype ignored for function " routine_name)
		else
		    warning("conflicting prototype ignored for function " routine_name)
	    }
	    else
	    {
		Prototype_Table[routine_name] = prototype
		Header_Table[routine_name] = FILENAME
		gsub("^.*/","",Header_Table[routine_name])
	    }
	    routine_name = ""
	    prototype = ""
	}
    }
    close(cmd)
    FILENAME = ""
    FNR = 0
}


function do_seealso_file(infile, k)
{
    while ((getline < infile) > 0)
    {
	for (k = 1; k <= NF; ++k)
	{
	    if ($k in See_Also)
		warning("already have See_Also[] entry for " $k)
	    else
		See_Also[$k] = $0
	}
    }
    close(infile)
}


function get_argument_variables(s, k,n,parts)
{
    ## Parse and store the argument variables in Argument_Table[]
    if (match(s,"[(][^)]*[)]"))  
    {
	n = split(substr(s,RSTART+1,RLENGTH-2),parts," *, *")
	for (k = 1; k <= n; ++k)
	    Argument_Table[parts[k]] = k
    }
}


function get_name(s)
{
    sub("^.*FUNCTION *", "", s)
    sub("^.*SUBROUTINE *","",s)
    sub(" *[(].*$","",s)
    return (tolower(s))
}


function gsplit(s,parts,regexp, n)
{
    ## Split s into parts[] according to regexp, and return the number of  
    ## entries in parts[].  Unlike split(), which discards the separator
    ## strings, gsplit() preserves them as members of parts[].

    split("",parts," ")		# clear parts[]
    n = 0
    while (match(s,regexp))
    {
	if (RSTART > 1)		# have NONMATCH MATCH REST
	    parts[++n] = substr(s,1,RSTART-1)
	parts[++n] = substr(s,RSTART,RLENGTH)
	s = substr(s,RSTART+RLENGTH)
    }
    if (s != "")
	parts[++n] = s
    return (n)
}


function header_file(name)
{
    return ((name in Header_Table) ? Header_Table[name] : "????")
}


function highlight_argument_variables(routine_name,s, k,n,lines)
{
    get_argument_variables(Prototype)  
    Argument_Table[routine_name] = 0 # we highlight the routine name too
    n = split(s,lines,"\n")
    s = ""
    for (k = 1; k <= n; ++k)    
	s = s ((s == "") ? "" : "\n") highlight_names(lines[k])
    return (s)
}


function highlight_embedded_names(s, k,n,parts,var_regexp)
{
    ## Handle more difficult cases of embedded names, like 
    ##
    ##    var(1..n)
    ##    var,
    ##    f(var)
    ##    f(g(h(var)))
    ##    (var1..var2)
    ##    (var1<=var2)
    ##    (var1.NE.var2)
    ##    (var1+var2)
    ##    (var1**2+var2**2)
    ##    ...
    ##
    ## These are sufficiently varied that a general solution is called
    ## for: split the string into a list of tokens that alternate
    ## between names and non-names, then supply embedded font changes
    ## for the names where needed.  The troff -man .B and .BR commands
    ## are not suited to this job, since they have a limited number of
    ## arguments, and introduce horizontal space if split into
    ## multiple such commands.

    var_regexp = "[A-Za-z][A-Za-z0-9_]*"
    n = gsplit(s,parts,var_regexp)
    s = ""
    for (k = 1; k <= n; k++)
	s = s ((parts[k] in Argument_Table) ? bold_inline(parts[k]) : parts[k])
    return (s)
}


function highlight_names(s, k,n,parts,var_regexp)
{
    n = split(s,parts," ")
    s = ""
    for (k = 1; k <= n; ++k) 
    {
	if (parts[k] in Argument_Table)
	    s = s " " bold_inline(parts[k])
	else
	    s = s " " highlight_embedded_names(parts[k])
    }
    return (trim_newlines_and_space(s))
}


function indent(n, blanks,s)
{
    s = ""
    blanks = "          "
    while (length(s) < n)
	s = s substr(blanks,1,(n <= length(blanks)) ? n : length(blanks))
    s = substr(s,1,n)
    return (s)
}


function initialize(cmd)
{
    cmd = "date"
    cmd | getline Current_Date_and_Time
    close(cmd)

    cmd = "date +'%d %B %Y'"
    cmd | getline DMY
    close (cmd)

    ## Initialize global variables to keep ``gawk --lint'' happy
    Declarations = ""
    Description = ""
    Indented_Lines = 0
    Last_Modified_Date = ""
    Manual_Section = (MANSECT == "") ? "3" : MANSECT
    Prototype = ""
    StarLines = 0
    Summary = ""

    Month_Expansion["jan"]	= "January"
    Month_Expansion["feb"]	= "February"
    Month_Expansion["mar"]	= "March"
    Month_Expansion["apr"]	= "April"
    Month_Expansion["may"]	= "May"
    Month_Expansion["jun"]	= "June"
    Month_Expansion["jul"]	= "July"
    Month_Expansion["aug"]	= "August"
    Month_Expansion["sep"]	= "September"
    Month_Expansion["oct"]	= "October"
    Month_Expansion["nov"]	= "November"
    Month_Expansion["dec"]	= "December"

    if (AUTHORFILE == "")
	AUTHORFILE = "AUTHORFILE.std"
    if (INCLUDEFILES != "")
	do_include_files(INCLUDEFILES)
    if (SEEALSO != "")
	do_seealso_file(SEEALSO)
}


function print_author()
{
    print_separator_comment()
    copy_file(AUTHORFILE)
}


function print_description(routine_name)
{
    print_separator_comment()
    print ".SH DESCRIPTION"
    print highlight_argument_variables(routine_name,Description)
}


function print_name(routine_name)
{
    print ".SH NAME"
    print routine_name " \\- " Summary
}


function print_seealso(routine_name, k,n,nout,parts)
{
    if (routine_name in See_Also)
    {
	n = split(See_Also[routine_name],parts," ")
	if (n > 1)
	{
	    print_separator_comment()
	    print ".SH \"SEE ALSO\""
	    nout = 0
	    for (k = 1; k <= n; ++k)
	    {
		if (parts[k] != routine_name)
		{
		    if (nout > 0)
			print ","
		    printf(".BR %s (%s)", parts[k], Manual_Section)
		    nout++
		}
	    }
	    if (nout > 0)
		print "."
	}
    }
}


function print_separator_comment()
{
    print ".\\\"====================================================================="
}


function print_synopsis(routine_name)
{
    print_separator_comment()
    print ".SH SYNOPSIS"
    print "Fortran (77, 90, 95, HPF):"
    print ".RS"
    print ".B f77"
    print ".I \"[ flags ] file(s) .\\|.\\|. -L/usr/local/lib -lgjl\""
    print ".RS"
    print ".nf"
    print ".B \"" Prototype "\""
    if (Declarations != "")
	print Declarations
    print ".fi"
    print ".RE"
    print ".RE"
    print "C (K&R, 89, 99), C++ (98):"
    print ".RS"
    print ".B cc"
    print ".I \"[ flags ] -I/usr/local/include file(s) .\\|.\\|. -L/usr/local/lib -lgjl\""
    print ".br"
    print "Use"
    print ".RS"
    print ".B \"#include <" header_file(routine_name) ">\""
    print ".RE"
    print "to get this prototype:"
    print ".RS"
    print ".nf"
    if ((routine_name == "") || (!(routine_name in Prototype_Table)))
    {
	warning("no C prototype for function " routine_name)
	print ".\\\" SUPPLY C PROTOTYPE HERE"
    }
    else
	print Prototype_Table[routine_name]
    print ".fi"
    print ".RE"
    print ".RE"
    print ".PP"
    print "NB: The definition of C/C++ data types"
    print ".B fortran_"
    print ".IR xxx ,"
    print "and the mapping of Fortran external names to C/C++ external names,"
    print "is handled by the C/C++ header file.  That way, the same function"
    print "or subroutine name can be used in C, C++, and Fortran code,"
    print "independent of compiler conventions for mangling of external"
    print "names in these programming languages."
    if (Last_Modified_Date != "")
    {
	print ".PP"
	print "Last code modification: " Last_Modified_Date
    }
}


function print_th(routine_name)
{
    print ".TH", toupper(routine_name), Manual_Section, ("\"" DMY "\""), "\"Version 1.00\""
    print ".\\\" WARNING: This file was produced automatically from file " FILENAME
    print ".\\\" by fortoman.awk on " Current_Date_and_Time "."
    print ".\\\" Any manual changes will be lost if this file is regenerated!"
}


function roman_bold(s,t)
{
    return (".RB \"" s "\" \"" t "\"")
}


function terminate( routine_name)
{
    routine_name = get_name(Prototype)
    print_th(routine_name)
    print_name(routine_name)
    print_synopsis(routine_name)
    print_description(routine_name)
    print_seealso(routine_name)
    print_author()
}


function trim(s)
{
    gsub(/^[ \t]+/,"",s)
    gsub(/[ \t]+$/,"",s)
    return (s)
}


function trim_newlines_and_space(s)
{
    ## Because awk and nawk (but not gawk or mawk) forbid newlines in
    ## character classes, we have to temporarily map the newline to a
    ## (almost certainly) unused character.
    gsub("\n","\001",s)
    gsub(/^[ \t\001]+/,"",s)
    gsub(/[ \t\001]+$/,"",s)
    gsub("\001","\n",s)
    return (s)
}


function troff_protect(s, quote_count)
{
    ## Protect lines that begin with dot, the troff command prefix
    ## character, provided the line is not a .PP, the only troff
    ## command which gets inserted above into s in place of an empty
    ## string.  
    ##
    ## Also, replace ellipses and quotes by the troff equivalents.

    gsub("[.][.][.]", ".\\|.\\|.", s)
    gsub("[.][.]",    ".\\|.",     s)
    
    quote_count = 0
    while (match(s,"\""))
    {
	if ((int(quote_count) % 2) == 0)
	    s = substr(s,1,RSTART-1) "``" substr(s,RSTART+RLENGTH)
	else
	    s = substr(s,1,RSTART-1) "''" substr(s,RSTART+RLENGTH)
	quote_count++
    }

    if (s == ".PP")
	return (s)
    else if (s ~ "^[.]")
	return ("\\&" s)
    else
	return (s)
}


function warning(message)
{
    print FILENAME ":" FNR ":%%" message >"/dev/stderr"
}
