
#-------------------------------------------------------------------------------#
#  NAME
#       m5, m5.awk - macro processor
#
#  SYNOPSIS
#       m5 [ -Dname ] [ -Dname=def ] [-c] [ -dp char ] [ -o file ] [
#       -sp char ] [ file ... ]
#
#       [g|n]awk -f m5.awk X [ -Dname ] [ -Dname=def ]  [-c]  [  -dp
#       char ] [ -o file ] [ -sp char ] [ file ... ]
#
#  DESCRIPTION
#       m5 is a Bourne shell script for invoking m5.awk, which actu-
#       ally   performs  the  macro  processing.   m5,  unlike  many
#       macroprocessors, does  not  directly  interpret  its  input.
#       Instead  it uses a two-pass approach in which the first pass
#       translates the input to an awk program, and the second  pass
#       executes  the  awk  program  to  produce  the  final output.
#       Details of usage are provided below.
#
#       As noted in the synopsis above, its invocation  may  require
#       specification  of  awk, gawk, or nawk, depending on the ver-
#       sion of awk  available  on  your  system.   This  choice  is
#       further  complicated  on  some systems, e.g. Sun, which have
#       both awk (original awk) and nawk (new awk).   Other  systems
#       appear to have new awk, but have named it just awk.  New awk
#       should be used, regardless of what it has been  named.   The
#       macro  processor translator will not work using original awk
#       because the former, for example, uses the built-in  function
#       match().
#
#  OPTIONS
#       The following options are supported:
#
#       -Dname        Following the cpp convention, define name as 1
#                     (one).   This  is  the  same  as if a -Dname=1
#                     appeared as an option or #name=1  appeared  as
#                     an  input  line.  Names specified using -D are
#                     awk variables  defined  just  before  main  is
#                     invoked.
#
#       -Dname=def    Define name as "def".  This is the same as  if
#                     #name="def"  appeared  as an input line. Names
#                     specified using -D are awk  variables  defined
#                     just before main is invoked.
#
#       X             Yes, that really is a capital "X".   The  ver-
#                     sion  of  nawk on Sun Solaris 2.5.1 apparently
#                     does its own argument processing before  pass-
#                     ing  the  arguments on to the awk program.  In
#                     this case, X (and all succeeding options)  are
#                     believed  by  nawk  to  be  file names and are
#                     passed on to the  macro  processor  translator
#                     (m5.awk)  for  its  own  argument processing).
#                     Without the X, Sun nawk  attempts  to  process
#                     succeeding  options  (e.g.,  -Dname)  as valid
#                     nawk  arguments  or  files,  thus  causing  an
#                     error.   This  may  not  be  a problem for all
#                     awks.
#
#       -c            Compile only.  The  output  program  is  still
#                     produced, but the final output is not.
#
#       -dp char      The directive prefix character (default is #).
#
#       -o file       The output program file (default is a.awk).
#
#       -sp char      The substitution prefix character (default  is
#                     $).
#
#  USAGE
#    Overview
#       The program that performs the  first  pass  noted  above  is
#       called  the m5 translator and is named m5.awk.  The input to
#       the translator may be either standard input or one  or  more
#       files  listed  on  the command line.  An input line with the
#       directive prefix character (# by default)  in  column  1  is
#       treated  as  a  directive  statement  in  the  MP  directive
#       language (awk).  All other input lines are processed as text
#       lines.   Simple  macros  are  created  using  awk assignment
#       statements and their values referenced using  the  substitu-
#       tion  prefix character ($ by default).  The backslash (\) is
#       the escape character; its presence forces the next character
#       to literally appear in the output.  This is most useful when
#       forcing the appearance of the  directive  prefix  character,
#       the  substitution prefix character, and the escape character
#       itself.
#
#    Macro Substitution
#       All input lines are scanned for macro  references  that  are
#       indicated  by  the  substitution prefix character.  Assuming
#       the default value of that character, macro references may be
#       of  the  form  $var, $(var), $(expr), $[str], $var[expr], or
#       $func(args).  These are replaced by  an  awk  variable,  awk
#       variable, awk expression, awk array reference to the special
#       array M[], regular awk  array  reference,  or  awk  function
#       call,  respectively.   These are, in effect, macros.  The MP
#       translator checks for proper nesting of parentheses and dou-
#       ble  quotes when translating $(expr) and $func(args) macros,
#       and checks for proper nesting of square brackets and  double
#       quotes  when translating $[expr] and $var[expr] macros.  The
#       substitution prefix character indicates a a macro  reference
#       unless  it  is (i) escaped (e.g., \$abc), (ii) followed by a
#       character other than A-Z, a-z, (, or [ (e.g., $@), or  (iii)
#       inside a macro reference (e.g., $($abc); probably an error).
#
#       An understanding of the implementation of macro substitution
#       will  help  in its proper usage. When a text line is encoun-
#       tered, it is scanned for macros, embedded in  an  awk  print
#       statement,  and  copied to the output program.  For example,
#       the input line
#
#       The quick $fox jumped over the lazy $dog.
#
#       is transformed into
#
#       print "The quick " fox " jumped over the lazy " dog "."
#
#       Obviously the use of this transformation technique relies completely
#       on the presence of the awk concatenation operator (one or more blanks).
#
#    Macros Containing Macros
#       As already noted, a macro  reference  inside  another  macro
#       reference  will not result in substitution and will probably
#       cause  an  awk   execution-time   error.    Furthermore,   a
#       substitution  prefix  character in the substituted string is
#       also generally not significant because the substitution pre-
#       fix  character  is  detected  at translation time, and macro
#       values are  assigned  at  execution  time.   However,  macro
#       references  of  the  form  $[expr]  provide  a simple nested
#       referencing capability.  For example, if $[abc] is in a text
#       line,  or  in a directive line and not on the left hand side
#       of   an   assignment   statement,   it   is   replaced    by
#       eval(M["abc"])/.   When  the output program is executed, the
#       m5 runtime routine eval()/ substitutes the value of M["abc"]
#       examining it for further macro references of the form $[str]
#       (where "str" denotes an arbitrary string).  If one is found,
#       substitution  and  scanning  proceed  recursively.  Function
#       type macro references may result in references to other mac-
#       ros,  thus  providing an additional form of nested referenc-
#       ing.
#
#    Directive Lines
#       Except for the include directive, when a directive  line  is
#       detected,  the  directive  prefix  is  removed,  the line is
#       scanned for macros, and then the line is copied to the  out-
#       put  program (as distinct from the final output).  Any valid
#       awk construct, including the function statement, is  allowed
#       in  a  directive  line.   Further information on writing awk
#       programs may be found in  Aho,  Kernighan,  and  Weinberger,
#       Dougherty and Robbins, and Robbins.
#
#    Include Directive
#       A single non-awk directive has been  provided:  the  include
#       directive.    Assuming  that  #  is  the  directive  prefix,
#       #include(filename) directs the MP translator to  immediately
#       read  from  the  indicated file, processing lines from it in
#       the normal manner.  This processing mode makes  the  include
#       directive  the  only  type  of  directive  to take effect at
#       translation time.  Nested  includes  are  allowed.   Include
#       directives  must  appear on a line by themselves.  More ela-
#       borate types of file processing may be  directly  programmed
#       using appropriate awk statements in the input file.
#
#    Main Program and Functions
#       The MP translator builds the resulting awk program in one of
#       two ways, depending on the form of the first input line.  If
#       that line begins with "function", it  is  assumed  that  the
#       user is providing one or more functions, including the func-
#       tion "main" required by m5.  If  the  first  line  does  not
#       begin  with  "function",  then  the  entire  input  file  is
#       translated  into  awk  statements  that  are  placed  inside
#       "main".   If some input lines are inside functions, and oth-
#       ers are not, awk will will detect this and complain.  The MP
#       by  design  has  little awareness of the syntax of directive
#       lines (awk statements), and as a consequence  syntax  errors
#       in directive lines are not detected until the output program
#       is executed.
#
#    Output
#       Finally, unless the -c (compile only) option is specified on
#       the  command line, the output program is executed to produce
#       the final output (directed by default to  standard  output).
#       The  version  of  awk  specified  in ARGV[0] (a built-in awk
#       variable containing the command name) is used to execute the
#       program.  If ARGV[0] is null, awk is used.
#
#  EXAMPLE
#       Understanding this example requires recognition  that  macro
#       substitution  is  a two-step process:  (i) the input text is
#       translated into an output awk  program,  and  (ii)  the  awk
#       program  is  executed  to  produce the final output with the
#       macro substitutions  actually  accomplished.   The  examples
#       below  illustrate  this  process.  # and $ are assumed to be
#       the directive  and  substitution  prefix  characters.   This
#       example  was  successfully  executed using awk on a Cray C90
#       running UNICOS 10.0.0.3, gawk on  a  Gateway  E-3200  runing
#       SuSE Linux Version 6.0, and nawk on a Sun Ultra 2 Model 2200
#       running Solaris 2.5.1.
#
#    Input Text
#       #function main() {
#
#          Example 1: Simple Substitution
#          ------------------------------
#       #  br = "brown"
#          The quick $br fox.
#
#          Example 2: Substitution inside a String
#          ---------------------------------------
#       #  r = "row"
#          The quick b$(r)n fox.
#
#          Example 3: Expression Substitution
#          ----------------------------------
#       #  a = 4
#       #  b = 3
#          The quick $(2*a + b) foxes.
#
#          Example 4: Macros References inside a Macro
#          -------------------------------------------
#       #  $[fox] = "\$[q] \$[b] \$[f]"
#       #  $[q] = "quick"
#       #  $[b] = "brown"
#       #  $[f] = "fox"
#          The $[fox].
#
#          Example 5: Array Reference Substitution
#          ---------------------------------------
#       #  x[7] = "brown"
#       #  b = 3
#          The quick $x[2*b+1] fox.
#
#          Example 6: Function Reference Substitution
#          ------------------------------------------
#          The quick $color(1,2) fox.
#
#          Example 7: Substitution of Special Characters
#          ---------------------------------------------
#       \#  The \$ quick \\ brown $# fox. $$
#       #}
#       #include(testincl.m5)
#
#    Included File testincl.m5
#       #function color(i,j) {
#          The lazy dog.
#       #  if (i == j)
#       #     return "blue"
#       #  else
#       #     return "brown"
#       #}
#
#    Output Program
#       function main() {
#          print
#          print "   Example 1: Simple Substitution"
#          print "   ------------------------------"
#          br = "brown"
#          print "   The quick " br " fox."
#          print
#          print "   Example 2: Substitution inside a String"
#          print "   ---------------------------------------"
#          r = "row"
#          print "   The quick b" r "n fox."
#          print
#          print "   Example 3: Expression Substitution"
#          print "   ----------------------------------"
#          a = 4
#          b = 3
#          print "   The quick " 2*a + b " foxes."
#          print
#          print "   Example 4: Macros References inside a Macro"
#          print "   -------------------------------------------"
#          M["fox"] = "$[q] $[b] $[f]"
#          M["q"] = "quick"
#          M["b"] = "brown"
#          M["f"] = "fox"
#          print "   The " eval(M["fox"]) "."
#          print
#          print "   Example 5: Array Reference Substitution"
#          print "   ---------------------------------------"
#          x[7] = "brown"
#          b = 3
#          print "   The quick " x[2*b+1] " fox."
#          print
#          print "   Example 6: Function Reference Substitution"
#          print "   ------------------------------------------"
#          print "   The quick " color(1,2) " fox."
#          print
#          print "   Example 7: Substitution of Special Characters"
#          print "   ---------------------------------------------"
#          print "\#  The \$ quick \\ brown $# fox. $$"
#       }
#       function color(i,j) {
#          print "   The lazy dog."
#          if (i == j)
#             return "blue"
#          else
#             return "brown"
#       }
#
#       function eval(inp   ,isplb,irb,out,name) {
#
#          splb = SP "["
#          out = ""
#
#          while( isplb = index(inp, splb) ) {
#             irb = index(inp, "]")
#             if ( irb == 0 ) {
#                out = out substr(inp,1,isplb+1)
#                inp = substr( inp, isplb+2 )
#             } else {
#                name = substr( inp, isplb+2, irb-isplb-2 )
#                sub( /^ +/, "", name )
#                sub( / +$/, "", name )
#                out = out substr(inp,1,isplb-1) eval(M[name])
#                inp = substr( inp, irb+1 )
#             }
#          }
#
#          out = out inp
#
#          return out
#       }
#       BEGIN {
#          SP = "$"
#          main()
#          exit
#       }
#
#    Final Output
#
#          Example 1: Simple Substitution
#          ------------------------------
#          The quick brown fox.
#
#          Example 2: Substitution inside a String
#          ---------------------------------------
#          The quick brown fox.
#
#          Example 3: Expression Substitution
#          ----------------------------------
#          The quick 11 foxes.
#
#          Example 4: Macros References inside a Macro
#          -------------------------------------------
#          The quick brown fox.
#
#          Example 5: Array Reference Substitution
#          ---------------------------------------
#          The quick brown fox.
#
#          Example 6: Function Reference Substitution
#          ------------------------------------------
#          The lazy dog.
#          The quick brown fox.
#
#          Example 7: Substitution of Special Characters
#          ---------------------------------------------
#       #  The $ quick \ brown $# fox. $$
#
#  FILE
#       a.awk         default output program file
#
#  SEE ALSO
#       awk(1), cpp(1), gawk(1), m4(1), nawk(1).  vi(1)
#
#  AUTHOR
#       William A. Ward, Jr., School  of  Computer  and  Information
#       Sciences, University of South Alabama, Mobile, Alabama, July
#       23, 1999.
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------#
#  GLOBAL NAMES
#
#     Awk Variable Handling
#       Awk's approach to handling variables is simple.  If a formal parameter
#       appears in the parameter list of an awk function, then it is a local
#       variable in that function; it is allocated on entry to the function
#       and deallocated on exit. It may be used even if (especially if) there
#       is no corresponding actual argument in the function invocation.  All
#       other awk variables are global.  In particular, a variable must be
#       global for its value to persist across function calls.  Furthermore,
#       there is no such thing as a named constant, only global (persistent)
#       variables that have their value set and are never changed.
#  
#     Naming Convention
#       The convention adopted here to distinguish between variable types is
#       to name global variables and constants in all upper case characters
#       and to name local variables in all lower case characters.  This is
#       consistent with awk's own convention for built-in variables.
#  
#     Global Constants (Name, Value, Description)
#       FALSE         0    Mnemonic for false
#       FS            "\n" Awk field separator is newline (awk builtin)
#       TRUE          1    Mnemonic for true
#  
#     Global Variables
#       ARGC          Number of command line arguments in ARGV[] (awk builtin)
#       ARGV[]        Array containing command line arguments (awk builtin)
#       C             Boolean flag for compile-only
#       DP            Directive prefix character
#       IMPLICIT_MAIN True if main program is implicit
#       ND            Number of -D type arguments
#       O             Name of executable output file
#       SP            Substitution prefix character
#
#------------------------------------------------------------------------------- 
#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_begin()
#
#  DESCRIPTION
#       m5_begin performs all initialization that must take place before the
#       first input line is read.  This includes (1) defining global constants,
#       (2) initializing global variables, and (3) processing command line
#       arguments in the argument array ARGV. Arguments that begin with the
#       dash (-) are treated as options while the remaining arguments are
#       assumed to be file names and are moved to the beginning of the argument
#       list array ARGV for processing by the program.  The valid arguments
#       are noted at the beginning of this file.
#
#  ARGUMENTS
#       None.
#
#  LOCAL VARIABLES
#       argc    - Number of original command line arguments.
#
#       i       - Argument counter
#
#-------------------------------------------------------------------------------

function m5_begin(	argc,i) {

   # Define global constants

   FALSE = 0					# Mnemonic for false
   TRUE = 1					# Mnemonic for true
   FS = "\n"					# awk field separator is newline

   # Set global variables

   NEED_EVAL = FALSE

   # Set default values for global variables which may be reset by arguments

   DP = "#"					# Default directive prefix
   C = FALSE					# Default for compile-only
   SP = "$"					# Default substitution prefix
   ND = 0					# Default number of -D arguments
   O = "a.awk"					# Default executable output file

   # Process the arguments in the array ARGV

   if ( ARGV[0] == "" )				# If command name is undefined
      ARGV[0] = "awk"				# ...Use "awk"
   argc = ARGC					# Save original value of ARGC
   ARGC = 1					# Reset to count file names

   for ( i=1; i<argc; i++ )			# Loop over every argument
      if ( ARGV[i] == "-D" )			# ...If arg is -D
         D[++ND] = ARGV[++i]			#    ...Save name[=def]
      else if ( substr(ARGV[i],1,2) == "-D" )	# ...Else if arg is -Dname[=def]
         D[++ND] = substr(ARGV[i],3)		#    ...Save name[=def]
      else if ( ARGV[i] == "-dp" )		# ...Else if arg is -dp <char>
         DP = ARGV[++i]				#    ...Save dir prefix char
      else if ( ARGV[i] == "-c" )		# ...Else if arg is -c
         C = TRUE				#    ...Set compile-only to true
      else if ( ARGV[i] == "-sp" )		# ...Else if arg is -sp <char>
         SP = ARGV[++i]				#    ...Save sub prefix char
      else if ( ARGV[i] == "-o" )		# ...Else if arg is -o <file>
         O = ARGV[++i]				#    ...Save name of exec file
      else if ( substr(ARGV[i],1,1) == "-" \
         && length(ARGV[i]) > 1 )		# ...Else if invalid option
            m5_prerr( "m5: Invalid argument " ARGV[i])
      else if ( i>1 || ARGV[i] != "X" )		# ...Else found a file name
         ARGV[ARGC++] = ARGV[i]			#    ...Put file name back in
						#          argument list

}						# End function m5_begin

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_dashd(d)
#
#  DESCRIPTION
#       m5_dashd processes "-D" (dashd) arguments from the command line
#       previously saved in the global array D.  This processing simply
#       requires writing an assignment statement (one per argument) in
#       the awk BEGIN section just before the call to main.
#
#  ARGUMENTS
#       d       - A string containing the macro definition, it must be of
#                 the form <name> or <name>=<value>.  If only the name is
#                 specified, the value is assumed to be 1 (one).  Other
#                 forms are treated as errors.  Macro definitions with
#                 embedded blanks are not allowed.
#  
#-------------------------------------------------------------------------------

function m5_dashd(d) {

   if ( ! match(d, /^[A-Za-z][A-Za-z0-9_]*/) )	# If invalid name
      m5_prerr("m5: Invalid argument -D " d)	# ...Display error message
   else if ( length(d) == RLENGTH )		# Else if -Dname
      print "   " substr(d,1,RLENGTH) \
         " = 1" > O				# ...Set name = 1 ala cpp
   else if ( substr(d,RLENGTH+1,1) != "=" )	# Else if next char not equal
      m5_prerr("m5: Invalid argument -D " d)	# ...Display error message
   else						# Else if -Dname=value
      print "   " substr(d,1,RLENGTH) \
         " = \"" substr(d,RLENGTH+2) "\"" > O	# ...Set name = value

}						# End function m5_dashd

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_end()
#
#  DESCRIPTION
#       m5_end performs wrap-up processing after the last input line is read.
#       First, if the last input line was still inside the function main, it
#       closes main with a "}".  Next, if the program should be executed, the
#       run-time system system and the awk BEGIN section are added to the output
#       program, and the program is executed using the awk processor specified
#       in ARGV[0] (see m5_begin).
#
#  ARGUMENTS
#       None.
#
#-------------------------------------------------------------------------------

function m5_end() {

   # If still inside main, close main with "}"

   m5_output(TRUE)

   # If the program should be executed (not compile-only), add the run-time
   # system and the awk BEGIN section to the output program, and then execute
   # the program.

   if ( ! C ) {					# If not compile-only
      if ( NEED_EVAL )
         m5_preval()				# ...Print run-time function
      print "BEGIN {" > O			# ...Print start of BEGIN sect
      print "   SP = \"" SP "\"" > O		# ...Print def of SP
      for ( i=1; i<=ND; i++)			# ...For each cmd line def
         m5_dashd(D[i])				#    ...Print cmd line def
      print "   main()" > O			# ...Print call to main
      print "   exit" > O			# ...Print exit
      print "}" > O				# ...Print end of BEGIN sect
      system( ARGV[0] " -f " O )		# ...Execute the awk program
   }						# End if not compile-only

}						# End function m5_end

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_findparen(s,rparen)
#
#  DESCRIPTION
#       m5_findparen finds the end of a substring in argument "s" containing
#       balanced parentheses and returns the position relative to the start of
#       the string s where this is the case.  The left and right parenthesis
#       characters are selectable (see below). Parentheses inside quoted strings
#       ("...") are ignored for the purpose of determining balance.  Escaped
#       quotes inside quoted strings are allowed.  Macro references which span
#       lines are not allowed.
#
#  ARGUMENTS
#       s       - A string containing the parenthesized portion of a macro
#                 invocation; if the $(ABC + DEF) is the invocation, then s
#                 contains (ABC + DEF).  This function should only be invoked
#                 when s begins with the left (open) parenthesis character as
#                 specified by argument lparen.
#
#       rparen  - The right parenthesis character.
#
#  LOCAL VARIABLES
#       ccur    - The current character extracted from s.
#
#       i       - Position of the current character in s.
#
#       escaped - True if the previous character was the escape character.
#
#       lparen  - The left parenthesis character.
#
#       nest    - The current nest depth.
#
#       quote   - True if inside a quoted string, false otherwise.
#
#-------------------------------------------------------------------------------

function m5_findparen(s,rparen		,ccur,escaped,i,lparen,nest,quoted) {

   # Initialize variables for processing the parenthesized string

   escaped = FALSE				# Turn escape flag off
   i = 1					# Init char position in string
   lparen = substr(s,1,1)			# Fetch left paren from string
   nest = 1					# Init paren nest depth
   quoted = FALSE				# Init quoted flag to off

   # Examine successive characters in sting s until the parenthesized nest
   # depth becomes zero (normal termination) or end-of-string is encountered
   # (an error condition)

   do {						# Do
      i = i + 1					# ...Move to next char position
      ccur = substr(s,i,1)			# ...Get next character
      if ( escaped )				# ...If escape flag is on
         escaped = FALSE			#    ...Skip this char
      else if ( ccur == "\\" )			# ...Else if char is escape
         escaped = TRUE				#    ...Turn escape flag on
      else if ( ccur == "\"" )			# ...Else if quote
         quoted = ! quoted			#    ...Flip quote indicator
      else if ( ccur == lparen )		# ...Else if a left paren.
         nest = nest + (! quoted)		#    ... increment nest depth
      else if ( ccur == rparen )		# ...Else if a right paren,
         nest = nest - (! quoted)		#    ... decrement nest depth
   } while( (ccur != "") && (nest > 0) )	# While not end & still nested

   # Check for error conditions

   if ( quoted )				# If quotes not closed
      m5_prerr("m5: Double quotes not closed in macro call $" s)
   else if ( nest > 0 )				# Else if still nested
      m5_prerr("m5: Unbalanced " lparen rparen "in macro call " SP s)

   return i					# Return position of right paren

}						# End function m5_findparen

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_include(inp)
#
#  DESCRIPTION
#       m5_include is called by m5_main to check for the presence of an
#       include directive.  If one is found, it is processed by extracting
#       the file name from the directive and then reading further input lines
#       from that file.  After each line is read, m5_main is recursively
#       called to process the line.
#
#  ARGUMENTS
#       inp     - The current input line.
#
#  LOCAL VARIABLES
#       file    - The name of the file to be included.
#
#       include - True if the input is an include directive, false otherwise.
#
#-------------------------------------------------------------------------------

function m5_include(inp		,file,include) {

   # Try to match input to " include ( "

   include = match( inp, /^[ \t]*include[ \t]*\([ \t]*/ )

   # If this is an include, process it by recursively calling m5_main

   if ( include ) {				# If this is an include,
      inp = substr( inp, RLENGTH+1 )		# ...Strip off "include ( "
      match( inp, /^[^) \t]+/ )			# ...Locate right paren
      if ( RSTART = 0 )				# ...If no right paren, error
         m5_prerr("m5: Badly formed include statement")
      file = substr( inp, 1, RLENGTH )		# ...Extract file name
      while( (getline < file) > 0 )		# ...While not eof on file,
         m5_main($0)				#    ...Recursively process
   }						# End if include

   return include				# Return include flag

}						# End function m5_include

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_macro(isdir)
#
#  DESCRIPTION
#       m5_macro processes a macro in the current input line.  The remaining
#       trailing portion of the input line Ris in the global variable INP.
#       The substitution prefix character has already been removed, so the
#       start of the macro reference text is the first character in INP.
#       The text to be copied to the output buffer is the returned value
#       of the function.
#
#  ARGUMENTS
#       isdir   - True (1) if this is a directive line, false (0) otherwise.
#
#  LOCAL VARIABLES
#       c       - Current input char.
#
#       iend    - Character position of the end of the macro reference.
#       
#       macro   - Contains the current macro reference with the substitution
#                 prefix character removed.
#
#-------------------------------------------------------------------------------

function m5_macro(isdir		,c,iend,macro) {

   # Get first character

   c = substr( INP, 1, 1 )			# Get first character

   # Is this $(...) ?

   if ( c == "(" ) {				# If mac ref is $(...)
      iend = m5_findparen(INP,")")		# ...Find closing )
      macro = substr( INP, 2, iend-2 )		# ...Extract macro

   # Is this $[...] ?

   } else if ( c == "[" ) {			# Else if mac ref is $[...]
      iend = m5_findparen(INP,"]")		# ...Find closing ]
      macro = substr( INP, 2, iend-2 )		# ...Extract macro
      sub( /^ +/, "", macro )			# ...Remove leading white space
      sub( / +$/, "", macro )			# ...Remove trailing white space
      if ( isdir && match( substr(INP,iend+1), \
            /[ \t]*=[^=]/ ) )			# ...If lhs of assignment stmt
         macro = "M[\"" macro "\"]"		#    ...Put M[macro]
      else {					# ...Else
         macro = "eval(M[\"" macro "\"])"	#    ...Put eval(M[macro])
         NEED_EVAL = TRUE			#    ...Will need runtime eval
      }						# ...End if

   # Is this $name, $name(...), or $name[...] ?

   } else if ( c ~ /^[A-Za-z]/ ) {		# Else if mac ref is $name...
      match( INP, /^[A-Za-z0-9_]+/ ) 		# ...Find end of name
      c = substr( INP, RLENGTH+1, 1 )		# ...Get 1st char after name
      if ( c == "(" )				# ...If function reference
         iend = RLENGTH + m5_findparen( \
            substr(INP,RLENGTH+1),")")		#    ...Find closing )
      else if ( c == "[" )			# ...If array reference
         iend = RLENGTH + m5_findparen( \
            substr(INP,RLENGTH+1),"]")		#    ...Find closing ]
      else					# ...Else
         iend = RLENGTH				#    ...Simple name
      macro = substr( INP, 1, iend )		# ...Extract macro

   # This is not a macro reference.

   } else {					# Else this is not a mac ref
      iend = 0					# ...So not to advance INP
      macro = ""				# ...Macro is null
   }						# End if

   # Adjust input string and return macro string

   INP = substr( INP, iend+1 )			# Advance input past mac ref
   return macro					# Return macro to put in output

}						# End function m5_macro

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_main(inp)
#
#  DESCRIPTION
#       To be supplied.
#
#  ARGUMENTS
#       inp     - On entry, the current input line.  inp is scanned for
#                 macros, and the macro-free prefix is appended to "out"
#                 (see below).  When inp becomes null, scanning is terminated
#                 and out is printed to the output program.
#
#  LOCAL VARIABLES
#       buf     - Macro-free text from inp accumulates in buf.  When a
#                 macro is encountered, buf is copied to out, and buf
#                 set to null.
#
#       c       - Current input char.
#
#       e       - If a text line, then the escape char, else null.
#
#       isdir   - True (1) if this is a directive line, false (0) otherwise.
#
#       out     - Text and macro references are copied to out from inp
#                 as inp is scanned; out is the final printed output.
#                 For text lines, out begins with "   print".
#
#       macro   - Contains the current macro reference with the substitution
#                 prefix character removed.
#
#       q       - If a text line, then the double quote char, else null.
#
#       s       - If a text line, then the space char, else null.
#
#-------------------------------------------------------------------------------

function m5_main(inp	,buf,c,e,isdir,macro,out,q,s) {

   # Determine input line type (directive or text) and initialization

   INP = inp					# Move inp to INP (global)
   isdir = ( substr(INP,1,1) == DP )		# Get input line type
   escaped = FALSE				# No escape prior to loop
   quoted = FALSE				# No escape prior to loop

   # If this is a directive line

   if ( isdir ) {				# If this is a directive line
      INP = substr( INP, 2 )			# ...Remove directive prefix
      if ( m5_include(INP) )			# ...If an include, process &
         return					#    ...Return
      if ( substr(INP,2,1) == " " )		# ...If 2nd char is space,
         buf = " "				#    ...Indent
      else					# ...Else
         buf = ""				#    ...No indent
      out = ""					# ...Awk command from input
      e = ""					# ...Escape char not needed
      q = ""					# ...Quote char not needed
      s = ""					# ...Space char not needed

   # Else this is a text line

   } else { 					# Else this is a text line.
      buf = ""					# ...Empty buffer
      out = "   print"				# ...Awk command is "print"
      e = "\\"					# ...Escape char is needed
      q = "\""					# ...Quote char is needed
      s = " "					# ...Space char is needed
   }						# End if

   # Loop to scan input line and process macro references.

   while ( (c = substr(INP,1,1)) != "" ) {	# While not end-of-string
      INP = substr( INP, 2 )			# ...Remove 1st char
      if ( escaped ) {				# ...If this char is escaped
         escaped = FALSE			#    ...Escape only one char
         buf = buf c				#    ...Copy current char
      } else if ( c == "\\" ) {			# ...Else if escape char
         escaped = TRUE				#    ...Turn on escape mode
         buf = buf e				#    ...Copy escape char
      } else if ( c == "\"" ) {			# ...Else if quote char
         quoted = ! quoted			#    ...Flip quote mode
         buf = buf e c				#    ...Copy escaped quote
      } else if ( c != SP ) {			# ...Else if not sub prefix
         buf = buf c				#    ...Copy current char
      } else {					# ...Else this is a macro ref
         macro = m5_macro(isdir)		#    ...Translate mac ref
         if ( macro == "" )			#    ...If not really a macro
            buf = buf c				#       ...Copy sub prefix
         else {					#    ...Else
            if ( buf != "" ) {			#       ...Buffer not empty?
               out = out s q buf q		#          ...Copy buf to out
               buf = ""				#          ...Reset buf
            }					#       ...End if	
            if ( isdir && quoted )		#       ...Mac in quote in dir?
               out = out "\" " macro " \""	#          ...Mac not quoted
            else				#       ...Else
               out = out s macro		#          ...Copy mac to out
         }					#    ...End if
      }						# ...End if
   }						# End while

   # Wrap-up processing

   if ( buf != "" )				# If buffer not empty
      out = out s q buf q			# ...Copy buffer to output
   m5_output(FALSE,out)				# Print output

}						# End function m5_main

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_output()
#
#  DESCRIPTION
#       m5_output performs all output except error messages.  It also
#       distinguishes between the situation where function main is or is not
#       user-provided.  A very simple test is used to do this: if the first word
#       on the line is "function", then it is assumed that the user has written
#       their own function main, and no special action is taken.  Otherwise, the
#       main function header and trailer are automatically inserted.
#
#  ARGUMENTS
#       eof     - True if this is end of file.
#
#       line    - The line to be printed.
#
#-------------------------------------------------------------------------------

function m5_output(eof,line) {
   NRO = NRO + 1				# Number of output records
   if ( NRO == 1 ) {				# If first line
      IMPLICIT_MAIN = \
         ! match(line, /^[ \t]*function[ \t]/)	# ...True if 1st word not func
      if ( IMPLICIT_MAIN )			# ...If main is implied
         print "function main() {" > O		#    ...Main function header
   }						# End if
   if ( eof && IMPLICIT_MAIN)			# If eof and main is implied
      print "}" > O				# ...Main function trailer
   else						# Else
      print line > O				# ...Print normal output line

}						# End function m5_output

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       m5_prerr(message)
#
#  DESCRIPTION
#       Handle fatal errors by printing the argument string "message" on
#       standard error and exiting with error code 1.
#
#  ARGUMENTS
#       message - Error message.
#
#-------------------------------------------------------------------------------

function m5_prerr(message) {

   print message | "cat >&2"			# Error message to std error.
   error = 1					# Set error condition.
   exit						# Terminate execution.

}						# End function m5_prerr

#-------------------------------------------------------------------------------
#
#  FUNCTION NAME
#       m5_preval()
#
#  DESCRIPTION
#       m5_preval adds the m5 runtime system to the output file.  Currently
#       this consists of a single function "eval" which is used to dynamically
#       evaluate macros of the form $[name] (where "$" is the substitution
#       prefix character).  eval is basically a scaled down version of m5_main.
#       Whereas m5_main replaces macro references with awk variable names at
#       translation time, eval performs true dynamic macro substitution at
#       execution time.
#
#       Why not make this a separate file?  Well, having the entire macro
#       processor in a single file has its advantages.  If the runtime system
#       were in a separate file, this program would have to be aware of the
#       file name of the separate file, and this could adversely affect the
#       portability and ease of installation of m5.
#
#  ARGUMENTS
#       None.
#
#-------------------------------------------------------------------------------

function m5_preval() {

   print "function eval(inp	,isplb,irb,out,name) {"			> O
   print ""								> O
   print "   splb = SP \"[\""						> O
   print "   out = \"\""						> O
   print ""								> O
   print "   while( isplb = index(inp, splb) ) {"			> O
   print "      irb = index(inp, \"]\")"				> O
   print "      if ( irb == 0 ) {"					> O
   print "         out = out substr(inp,1,isplb+1)"			> O
   print "         inp = substr( inp, isplb+2 )"			> O
   print "      } else {"						> O
   print "         name = substr( inp, isplb+2, irb-isplb-2 )"		> O
   print "         sub( /^ +/, \"\", name )"				> O
   print "         sub( / +$/, \"\", name )"				> O
   print "         out = out substr(inp,1,isplb-1) eval(M[name])"	> O
   print "         inp = substr( inp, irb+1 )"				> O
   print "      }"							> O
   print "   }"								> O
   print ""								> O
   print "   out = out inp"						> O
   print ""								> O
   print "   return out"						> O
   print "}"								> O

}						# End function m5_preval

#-------------------------------------------------------------------------------
#
#  FUNCTION
#       BEGIN, awk main, and END
#
#  DESCRIPTION
#       Top level procedure invocations.
#
#  ARGUMENTS
#       None.
#
#-------------------------------------------------------------------------------

BEGIN {						# Initialization...
   m5_begin()					# ...prior to reading any input
}						# End BEGIN

{						# Process each line...
   m5_main($0)					# ...to create output program
}						# End awk main loop

END {						# Wrap-up processing...
   m5_end()					# ...after all lines read
}						# End END
