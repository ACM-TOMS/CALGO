%{
// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

/* This will take an input file and split it into a configurations file and a function file */
#include "bertini.h"

// settings to switch between config & function
int configFound = 0, endOfConfig = 1, endOfFunc = 0;
// files for config & function
FILE *CONFIGFILE = NULL, *FUNCFILE = NULL;

/* %option nounput - do not warn about unput */
%}
%option nounput 

config "CONFIG"
input "INPUT"
end ("END"[ \t\r\n]+)
ends "END;"
homvargp "hom_variable_group"
vargp "variable_group"
pathvar "pathvariable"
var "variable"
param "parameter"
const "constant"
func "function"
subfunc "definedSubfunction"
rand "random"
randReal "random_real"
number ([0-9]+[\.]?([0-9]*)?([eE][-+]?[0-9]+)?)|([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?) 
I [I]
name [a-zA-Z][a-zA-Z0-9_\[\]]*
comments [%].?+[\n]
white [ \t\r\n]
eol [;]
colon [:]
other .

%%

{end}|{ends} { // end of something
             if (!endOfConfig)
             { // end config
               endOfConfig = 1;
             }
             else if (!endOfFunc)
             { // end func
               endOfFunc = 1;
               // print END;
               fprintf(FUNCFILE, "\nEND;\n\n");
             }
           }
{comments}|{white} { // ignore comments and whitespace
           }
{config}   { // start of config found - print there until reach END or INPUT
             configFound = 1;
             endOfConfig = 0;
           }
{input}    { // start of input found - print there until reach END or EOF
             endOfConfig = 1;
             endOfFunc = 0;
           }
{param}|{homvargp}|{vargp}|{var}|{pathvar}|{const}|{func}|{subfunc}|{rand}|{randReal} { // declaration
             if (endOfConfig && !endOfFunc)
             { // print to FUNCFILE
               fprintf(FUNCFILE, "%s ", splitParsetext);
             }
           }
{colon}    { // a colon
             if (!endOfConfig)
             { // add a space after a colon in configurations
               fprintf(CONFIGFILE, ": ");
             }
           }
{eol}      { // end of line
             if (!endOfConfig)
             { // print to CONFIGFILE
  	       fprintf(CONFIGFILE, ";\n");
             }
             else if (!endOfFunc)
             { // print to FUNCFILE
               fprintf(FUNCFILE, ";\n");
             }
           }
{number}|{I}|{name}|{other}   { // other symbol
             if (!endOfConfig)
             { // print to CONFIGFILE
               fprintf(CONFIGFILE, "%s", splitParsetext);
             }
             else if (!endOfFunc)
             { // print to FUNCFILE
               fprintf(FUNCFILE, "%s", splitParsetext);
             }
           }
%%

// parse a file
int splitParse(FILE *fp, char *funcName, char *configName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: splits fp into function & configurations               *
\***************************************************************/
{
  int rV = 0;

  // setup FUNCFILE & CONFIGFILE
  FUNCFILE = fopen(funcName, "w");
  CONFIGFILE = fopen(configName, "w");
  if (FUNCFILE == NULL)
  {
    printf("ERROR: The name %s is not valid!\n\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }
  if (CONFIGFILE == NULL)
  {
    printf("ERROR: The name %s is not valid!\n\n", configName);
    bexit(ERROR_INPUT_SYNTAX);
  }

  // parse the file
  splitParsein = fp;
  rV = splitParselex();

  // close FUNCFILE & CONFIGFILE
  if (!endOfFunc)
    fprintf(FUNCFILE, "\nEND;\n\n");

  fclose(FUNCFILE);
  fclose(CONFIGFILE);
  FUNCFILE = NULL;
  CONFIGFILE = NULL;

  return rV;
}

int splitParsewrap()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches EOF                         *
\***************************************************************/
{ // parser is complete
  return 1;
}

int splitParseerror(char *s)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches an error-should never occur *
\***************************************************************/
{ // error in parser
  printf("ERROR: There was an error (%s) near '%s'.\n", s, splitParsetext);
  bexit(ERROR_INPUT_SYNTAX);

  return 0;
}




