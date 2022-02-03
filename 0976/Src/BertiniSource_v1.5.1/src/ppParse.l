%{
// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler
#include "ppParse_old.h"
#include "ppParse.tab.h"

// preprocess array
preProcessArray *ppArray = NULL;

// processing line
char *ppLine = NULL;
int ppLineLen = 0;

// using a parameter homotopy
int ppUsingParamHomotopy = 0, ppAddAllNums = 0;

/* %option nounput - do not warn about unput */
%}
%option nounput

end "END"
homvargp "hom_variable_group "
vargp "variable_group "
pathvar "pathvariable "
var "variable "
param "parameter "
const "constant "
func "function "
subfunc "definedSubfunction "
number ([0-9]+[\.]?([0-9]*)?([eE][-+]?[0-9]+)?)|([0-9]*\.[0-9]+([eE][-+]?[0-9]+)?)
I [I]
Pi "Pi"
sin "sin"
cos "cos"
tan "tan"
exp "exp"
pow "pow"
name [a-zA-Z][a-zA-Z0-9_\[\]]*
eol ";"
white [ \t\r\n]
operator [-+*/^(),=]
other .

%%

{end}      { // end of file
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppEND;
           }
{homvargp} { // hom variable group
	     append_str(&ppLine, &ppLineLen, ppParsetext);
             return ppHOMVARGP;
           }
{vargp}    { // variable group
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppVARGP;
           }
{pathvar}  { // path variable
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppPATHVAR;
           }
{var}      { // variable
	     append_str(&ppLine, &ppLineLen, ppParsetext);
      	     return ppVAR;
           }
{param}    { // parameter
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppPARAM;
           }
{const}    { // const
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppCONST;
           }
{func}     { // function
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppFUNC;
           }
{subfunc}  { // defined subfunction
	     append_str(&ppLine, &ppLineLen, ppParsetext);
             return ppSUBFUNC;
           } 
{number}   { // number
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     ppParselval.name = strdup(ppParsetext);
	     return ppNUMBER;
           }
{sin}|{cos}|{tan}|{exp} { // sin, cos, tan, or exp
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     ppParselval.name = strdup(ppParsetext);
	     return ppSCE;
	   }
{pow}      { // pow 
             append_str(&ppLine, &ppLineLen, ppParsetext);
             ppParselval.name = strdup(ppParsetext);
             return ppPOW;
           }
{I}        { // I
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppI;
           }
{Pi}       { // Pi
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppPi;
           }
{name}     { // name
	     append_str(&ppLine, &ppLineLen, ppParsetext);
             ppParselval.name = strdup(ppParsetext);
	     return ppNAME;
           }
{eol}      { // end of line
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppEOL;
           }
{operator} { // -+*/^(),=
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     return ppParsetext[0];
           }
{white}    ; // ignore it
{other}    { // other symbol - error!
	     append_str(&ppLine, &ppLineLen, ppParsetext);
	     ppParseerror("invalid character");
           }

%%

void append_str(char **str1, int *lenStr1, char *str2)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: append str2 onto str1 and update its length            *
\***************************************************************/
{
  *lenStr1 += strlen(str2);
  *str1 = (char *)brealloc(*str1, (*lenStr1 + 1) * sizeof(char));
  strcat(*str1, str2);

  return;
}

int lookup_preProcessArray_type(int *loc, preProcessArrayType *Array, char *name)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - found, 0 - not found                       *
* NOTES: lookup an entry in ppType                              *
\***************************************************************/
{
  int j, rV = 0, len = strlen(name);

  *loc = -1; 

  for (j = 0; j < Array->numType && !rV; j++)
    if (len == Array->nameSize[j] && !strcmp(name, Array->name[j]))
    { // found a match!
      rV = 1;
      *loc = j;
    }

  return rV;
}

int lookup_preProcessArray(int *type, int *loc, preProcessArray *Array, char *name)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - found, 0 - not found                       *
* NOTES: lookup an entry in Array                               *
\***************************************************************/
{
  int i, rV = 0, numTypes = Array->numTypes;

  *type = *loc = -1;

  for (i = 0; i < numTypes && !rV; i++)
  {
    rV = lookup_preProcessArray_type(loc, &Array->types[i], name);
    if (rV)
    { // found a match!
      *type = i;
    }
  }

  return rV;
}

void negate_lastEntry(preProcessArrayType *Numbers)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: negate the last entry in Numbers                       *
\***************************************************************/
{
  int loc = Numbers->numType - 1;
  char *oldName = Numbers->name[loc];

  Numbers->nameSize[loc]++;
  Numbers->name[loc] = (char *)bmalloc((Numbers->nameSize[loc] + 1) * sizeof(char));
  sprintf(Numbers->name[loc], "-%s", oldName);

  free(oldName);

  return;
}


void addEntry_preProcessArray_type(preProcessArrayType *Array, char *name, int isDefined, int lineNo)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add an entry to Array                                  *
\***************************************************************/
{
  int loc = Array->numType;

  // add the entry
  Array->numType++;

  if (Array->numType > Array->allocSize)
  { // increase the allocated size
    if (Array->allocSize == 0)
      Array->allocSize = 1;
    else
      Array->allocSize *= 2;

    Array->isDefined = (int *)brealloc(Array->isDefined, Array->allocSize * sizeof(int));
    Array->lineNumber = (int *)brealloc(Array->lineNumber, Array->allocSize * sizeof(int));
    Array->nameSize = (int *)brealloc(Array->nameSize, Array->allocSize * sizeof(int));
    Array->name = (char **)brealloc(Array->name, Array->allocSize * sizeof(char *));
  }

  Array->isDefined[loc] = isDefined;
  Array->lineNumber[loc] = lineNo;
  Array->nameSize[loc] = strlen(name);
  Array->name[loc] = (char *)bmalloc((Array->nameSize[loc] + 1) * sizeof(char));
  strcpy(Array->name[loc], name);

  return;
}

void initialize_preProcessArray_type(preProcessArrayType *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes Array                                      *
\***************************************************************/
{
  Array->numType = 0;
  Array->allocSize = 0;
  Array->isDefined = NULL;
  Array->lineNumber = NULL;
  Array->nameSize = NULL;
  Array->name = NULL;

  return;
}

void initialize_preProcessArray(preProcessArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes Array                                      *
\***************************************************************/
{
  int i;

  Array->numTypes = 10;
  Array->types = (preProcessArrayType *)bmalloc(Array->numTypes * sizeof(preProcessArrayType));
  for (i = 0; i < Array->numTypes; i++)
    initialize_preProcessArray_type(&Array->types[i]);

  // setup 'I', 'Pi', '0' & '1'
  addEntry_preProcessArray_type(&Array->types[CONSTANTTYPE], "I", 1, 0);
  addEntry_preProcessArray_type(&Array->types[CONSTANTTYPE], "Pi", 1, 0);
  addEntry_preProcessArray_type(&Array->types[NUMBERTYPE], "0", 1, 0);
  addEntry_preProcessArray_type(&Array->types[NUMBERTYPE], "1", 1, 0);
 
  return;
}

void copy_preProcessArray_type(preProcessArrayType *Out, preProcessArrayType *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  int i;

  Out->numType = In->numType;

  if (In->numType > 0)
  { // setup the structures
    Out->allocSize = Out->numType;
    Out->isDefined = (int *)bmalloc(Out->numType * sizeof(int));
    Out->lineNumber = (int *)bmalloc(Out->numType * sizeof(int));
    Out->nameSize = (int *)bmalloc(Out->numType * sizeof(int));
    Out->name = (char **)bmalloc(Out->numType * sizeof(char *));

    for (i = 0; i < In->numType; i++)
    {
      Out->isDefined[i] = In->isDefined[i];
      Out->lineNumber[i] = In->lineNumber[i];
      Out->nameSize[i] = In->nameSize[i];
      Out->name[i] = (char *)bmalloc((Out->nameSize[i] + 1) * sizeof(char));
      strcpy(Out->name[i], In->name[i]);
    }
  }
  else
  { // NULL out the pointers
    Out->allocSize = 0;
    Out->isDefined = NULL;
    Out->lineNumber = NULL;
    Out->nameSize = NULL;
    Out->name = NULL;
  }

  return;
}

void copy_preProcessArray(preProcessArray *Out, preProcessArray *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  int i;

  Out->numTypes = In->numTypes;

  if (In->numTypes > 0)
  { // copyt data
    Out->types = (preProcessArrayType *)bmalloc(Out->numTypes * sizeof(preProcessArrayType));
    for (i = 0; i < In->numTypes; i++)
      copy_preProcessArray_type(&Out->types[i], &In->types[i]);
  }
  else
  { // NULL out pointer 
    Out->types = NULL;
  }

  return;
}

void clear_preProcessArray_type(preProcessArrayType *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears Array                                           *
\***************************************************************/
{
  int i;

  if (Array->numType > 0)
  {
    for (i = 0; i < Array->numType; i++)
      free(Array->name[i]);
    free(Array->isDefined);
    free(Array->lineNumber);
    free(Array->nameSize);
    free(Array->name);
  }
  
  return;
}

void clear_preProcessArray(preProcessArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears Array                                           *
\***************************************************************/
{
  int i;

  if (Array->numTypes > 0)
  {
    for (i = 0; i < Array->numTypes; i++)
      clear_preProcessArray_type(&Array->types[i]);
  
    free(Array->types);
    Array->numTypes = 0;
    Array->types = NULL;
  }

  return;
}

void verify_preProcessArray_type(preProcessArrayType *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify that everything is defined in Array             *
\***************************************************************/
{
  int i;

  for (i = 0; i < Array->numType; i++)
    if (!Array->isDefined[i])
    { // not defined
      printf("ERROR: After parsing, it appears that the symbol '%s' is not defined!\n", Array->name[i]);
      bexit(ERROR_INPUT_SYNTAX);
    }
 
  return;
}

void verify_existence_definedSubfunction(preProcessArrayType *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify that everything is defined in Array has a file  *
\***************************************************************/
{
  int i, size = 0;
  char *name = NULL;
  FILE *IN = NULL;

  for (i = 0; i < Array->numType; i++)
  {
    size = 1 + snprintf(NULL, 0, "%s.func", Array->name[i]);
    name = (char *)brealloc(name, size * sizeof(char));
    sprintf(name, "%s.func", Array->name[i]);
    IN = fopen(name, "r");
    if (IN == NULL)
    { // file does not exist
      printf("ERROR: The file '%s' for the defined subfunction '%s' does not exist!\n", name, Array->name[i]);
      bexit(ERROR_INPUT_SYNTAX);
    }
  }

  if (name != NULL)
    free(name);

  return;
}

void verify_preProcessArray(preProcessArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify that everything is defined in Array             *
\***************************************************************/
{
  int i;

  // verify the data
  for (i = 0; i < Array->numTypes; i++)
    verify_preProcessArray_type(&Array->types[i]);

  // verify that we have defined subfunction files
  verify_existence_definedSubfunction(&Array->types[DEFINEDSUBFUNCTIONTYPE]);

  // verify that there are at most 1 path variables
  if (Array->types[PATHVARIABLETYPE].numType > 1)
  { // error - too many
    printf("ERROR: After parsing, it appears that there are %d path variables defined (there must be at most 1)!\n", Array->types[PATHVARIABLETYPE].numType);
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void verify_standard_homotopy(preProcessArray *Array, variablegroupArray *vargpArray, int onlyOneGp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify settings for a standard homotopy                *
\***************************************************************/
{
  // check path variables
  if (Array->types[PATHVARIABLETYPE].numType != 0)
  { // error - too many
    printf("ERROR: After parsing, it appears that at least one path variable has been defined!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check parameters
  if (Array->types[PARAMETERTYPE].numType != 0)
  { // error - too many
    printf("ERROR: After parsing, it appears that at least one parameter has been defined!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check variable
  if (Array->types[VARIABLETYPE].numType != 0)
  { // error - too many
    printf("ERROR: After parsing, it appears that at least one 'variable' has been used!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check variable groups
  if (onlyOneGp)
  { // check that there is only 1 variable group
    if (vargpArray->numGps > 1)
    { // error - too many
      printf("ERROR: Your settings only allow for one variable group to be used (you have requested to use %d variable groups)!\n", vargpArray->numGps);
      bexit(ERROR_INPUT_SYNTAX);
    }
  }

  return;
}

void verify_parameter_homotopy(preProcessArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify settings for a standard homotopy                *
\***************************************************************/
{
  // check path variables
  if (Array->types[PATHVARIABLETYPE].numType != 0)
  { // error - too many
    printf("ERROR: After parsing, it appears that at least one path variable has been defined!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check parameters
  if (Array->types[PARAMETERTYPE].numType < 1)
  { // error - too many
    printf("ERROR: After parsing, it appears that there are %d parameters defined (there must be at least 1)!\n", Array->types[PARAMETERTYPE].numType);
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check variable
  if (Array->types[VARIABLETYPE].numType != 0)
  { // error - too many
    printf("ERROR: After parsing, it appears that at least one 'variable' has been used!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void verify_userdefined_homotopy(preProcessArray *Array, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify settings for a standard homotopy                *
\***************************************************************/
{
  // check path variables
  if (Array->types[PATHVARIABLETYPE].numType != 1)
  { // error - too many
    printf("ERROR: After parsing, it appears that there are %d path variables defined (there must be exactly 1)!\n", Array->types[PATHVARIABLETYPE].numType);
    bexit(ERROR_INPUT_SYNTAX);
  }

  // check parameters
  if (Array->types[PARAMETERTYPE].numType < 1)
  { // error - too many
    printf("ERROR: After parsing, it appears that there are %d parameters defined (there must be at least 1)!\n", Array->types[PARAMETERTYPE].numType);
    bexit(ERROR_INPUT_SYNTAX);
  }

  if (userHom == 2)
  { // check variable
    if (Array->types[VARIABLETYPE].numType != 0)
    { // error - too many
      printf("ERROR: After parsing, it appears that at least one 'variable' has been used!\n");
      printf("With this setting, you should use either 'variable_group' or 'hom_variable_group'.\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
  }
  else
  { // check variable_group
    if (Array->types[VARIABLEGROUPTYPE].numType != 0)
    { // error - too many
      printf("ERROR: After parsing, it appears that at least one 'variable_group' has been used!\n");
      printf("Standard user defined homotopies should use 'variable' rather than 'variable_group'.\n");
      bexit(ERROR_INPUT_SYNTAX);
    }

    // check hom variable_group
    if (Array->types[HOMVARIABLEGROUPTYPE].numType != 0)
    { // error - too many
      printf("ERROR: After parsing, it appears that at least one 'hom_variable_group' has been used!\n");
      printf("Standard user defined homotopies should use 'variable' rather than 'hom_variable_group'.\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
  }

  return;
}

int preProcessParse(preProcessArray *preProcArray, FILE *fp, int paramHom, int addAllNums)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: preprocess the file                                    *
\***************************************************************/
{
  int rV = 0, endFile = 0;

  // initialize
  ppLine = NULL;
  endOfFile = &endFile;

  // setup ppArray
  ppArray = preProcArray;
  initialize_preProcessArray(ppArray);

  // setup line number
  lineNum = &ppParselineno;
  *lineNum = 1;

  // setup ppUsingParamHomotopy
  ppUsingParamHomotopy = (paramHom == 1 || paramHom == 2);

  // setup ppAddAllNums
  ppAddAllNums = addAllNums;

  // parse the file
  ppParsein = fp;
  do
  { // clear out the line
    ppLineLen = 1;
    ppLine = (char *)brealloc(ppLine, ppLineLen * sizeof(char));
    ppLine[0] = '\0';

    // parse the line
    rV = ppParseparse(); 

    // increment line number
    (*lineNum)++;

  } while (!rV && !endFile);

  // verify that everything that needs defined has been defined
  verify_preProcessArray(ppArray);

  // clear memory
  ppParsein = NULL;
  ppArray = NULL;

  free(ppLine);
  ppLineLen = 0;

  return rV;
}

int ppParsewrap()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches EOF                         *
\***************************************************************/
{ // parser is complete
  return 1;
}

int ppParseerror(char *s)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when preprocessor reaches an error              *
\***************************************************************/
{ // error in parser
  if (!strcmp(s,"invalid character"))
    printf("ERROR: The '%s' is an %s in system statement %d.\n", ppParsetext, s, ppParselineno);
  else
    printf("ERROR: There was an error (%s) near '%s' in system statement %d.\n", s, ppParsetext, ppParselineno);
  errorStatement();

  return 0;
}

void errorStatement()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches an error-should never occur *
\***************************************************************/
{
  printf("Statement %d up to error: '%s'\n", ppParselineno, ppLine);
  bexit(ERROR_INPUT_SYNTAX);

  return;
}

void print_preproc_data(FILE *OUT, int numFuncs, int numVarGps, int numHomGps, int *types, int *sizes)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the required data to OUT                        *
\***************************************************************/
{
  int i, numGps = numVarGps + numHomGps;

  // print the number of functions, number of homogeneous variable groups, and number of variable groups
  fprintf(OUT, "%d %d %d\n", numFuncs, numHomGps, numVarGps);
  // print the type & size of each of the groups
  for (i = 0; i < numGps; i++)
    fprintf(OUT, "%d %d\n", types[i], sizes[i]);

  return;
}

void setup_variablegroup_numbers(preProcessArray *Array, int *numHomGps, int **sizeHomGps, int *numVarGps, int **sizeVarGps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the variable group numbers                       *
\***************************************************************/
{
  int i, j, found, lineNum, lineCount = 0, *lines = NULL, *types = NULL, currVar = 0, currHom = 0;
  int numVar = Array->types[VARIABLEGROUPTYPE].numType, numHomVar = Array->types[HOMVARIABLEGROUPTYPE].numType;

  // initialize counts
  *numHomGps = *numVarGps = 0;

  // count the number of line numbers used in var & hom_var
  for (i = 0; i < numVar; i++)
  { // see if used
    found = 0;
    lineNum = Array->types[VARIABLEGROUPTYPE].lineNumber[i];
    for (j = 0; j < lineCount && !found; j++)
      if (lines[j] == lineNum)
        found = 1;

    if (!found)
    { // add to lines & types
      lines = (int *)brealloc(lines, (lineCount + 1) * sizeof(int));
      types = (int *)brealloc(types, (lineCount + 1) * sizeof(int));
      lines[lineCount] = lineNum;
      types[lineCount] = 1; // variable group
      lineCount++;

      // increment the count of variable groups
      (*numVarGps)++;
    }
  }
  for (i = 0; i < numHomVar; i++)
  { // see if used
    found = 0;
    lineNum = Array->types[HOMVARIABLEGROUPTYPE].lineNumber[i];
    for (j = 0; j < lineCount && !found; j++)
      if (lines[j] == lineNum)
        found = 1;

    if (!found)
    { // add to lines
      lines = (int *)brealloc(lines, (lineCount + 1) * sizeof(int));
      types = (int *)brealloc(types, (lineCount + 1) * sizeof(int));
      lines[lineCount] = lineNum;
      types[lineCount] = 0; // hom variable group
      lineCount++;

      // increment the count of hom variable groups
      (*numHomGps)++;
    }
  }

  // allocate sizeHomGps & sizeVarGps
  *sizeVarGps = (int *)bmalloc(*numVarGps * sizeof(int));
  *sizeHomGps = (int *)bmalloc(*numHomGps * sizeof(int));

  // setup the group numbers
  for (i = 0; i < lineCount; i++)
  { // find the minimum still left
    lineNum = INT_MAX; 
    found = 0;
    for (j = 0; j < lineCount; j++)
      if (lines[j] < lineNum && lines[j] >= 0)
      { // store this line number
        lineNum = lines[j];
        found = j;
      }

    // determine the type and setup the line numbers and count
    if (types[found])
    { // variable group
      (*sizeVarGps)[currVar] = 0;
      for (j = 0; j < numVar; j++)
        if (Array->types[VARIABLEGROUPTYPE].lineNumber[j] == lineNum)
        {
          Array->types[VARIABLEGROUPTYPE].lineNumber[j] = i;
          (*sizeVarGps)[currVar]++;
        }

      currVar++;
    }
    else
    { // hom variable group
      (*sizeHomGps)[currHom] = 0;
      for (j = 0; j < numHomVar; j++)
        if (Array->types[HOMVARIABLEGROUPTYPE].lineNumber[j] == lineNum)
        {
          Array->types[HOMVARIABLEGROUPTYPE].lineNumber[j] = i;
          (*sizeHomGps)[currHom]++;
        }

      currHom++;
    }

    // kill off this line number
    lines[found] = -1;
  }

  // clear memory
  free(types);
  free(lines);

  return;
}

void initialize_variablegroupArray(variablegroupArray *vargpArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize vargpArray                                  *
\***************************************************************/
{
  vargpArray->numGps = vargpArray->numHomGps = vargpArray->numVarGps = vargpArray->totalNumVars = 0;
  vargpArray->types = vargpArray->sizes = vargpArray->varLoc = vargpArray->homLoc = vargpArray->groupNumber = vargpArray->nameSize = NULL;
  vargpArray->name = NULL;

  return;
}

void clear_variablegroupArray(variablegroupArray *vargpArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear vargpArray                                       *
\***************************************************************/
{
  int i;

  // clear name
  for (i = 0; i < vargpArray->totalNumVars; i++)
    free(vargpArray->name[i]);
  free(vargpArray->name);

  // clear other data
  free(vargpArray->types);
  free(vargpArray->sizes);
  free(vargpArray->varLoc);
  free(vargpArray->homLoc);
  free(vargpArray->groupNumber);
  free(vargpArray->nameSize);

  // NULL out data
  vargpArray->numGps = vargpArray->numHomGps = vargpArray->numVarGps = vargpArray->totalNumVars = 0;
  vargpArray->types = vargpArray->sizes = vargpArray->varLoc = vargpArray->homLoc = vargpArray->groupNumber = vargpArray->nameSize = NULL;
  vargpArray->name = NULL;

  return;
}

void setup_variablegroupArray(preProcessArray *Array, variablegroupArray *vargpArray, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup vargpArray using Array                           *
\***************************************************************/
{
  int i, j, found, lineNum, lineCount = 0, *lines = NULL, *types = NULL, currVar = 0, currVarGp = 0;
  int numVar = Array->types[VARIABLEGROUPTYPE].numType, numHomVar = Array->types[HOMVARIABLEGROUPTYPE].numType;

  // initialize vargpArray
  initialize_variablegroupArray(vargpArray);

  if (userHom)
  { // just using 'variable' data
    vargpArray->totalNumVars = Array->types[VARIABLETYPE].numType;
    vargpArray->numGps = vargpArray->numHomGps = vargpArray->numVarGps = 0;
    vargpArray->types = vargpArray->sizes = vargpArray->varLoc = vargpArray->homLoc = vargpArray->groupNumber = vargpArray->nameSize = NULL;
    // setup name
    vargpArray->name = (char **)bmalloc(vargpArray->totalNumVars * sizeof(char *));
    for (i = 0; i < vargpArray->totalNumVars; i++)
    {
      vargpArray->name[i] = (char *)bmalloc((Array->types[VARIABLETYPE].nameSize[i] + 1) * sizeof(char));
      strcpy(vargpArray->name[i], Array->types[VARIABLETYPE].name[i]);
    }
  }
  else
  { // using both 'variable_group' & 'hom_variable_group'

    // count the number of line numbers used in var & hom_var
    for (i = 0; i < numVar; i++)
    { // see if used
      found = 0;
      lineNum = Array->types[VARIABLEGROUPTYPE].lineNumber[i];
      for (j = 0; j < lineCount && !found; j++)
        if (lines[j] == lineNum)
          found = 1;

      if (!found)
      { // add to lines & types
        lines = (int *)brealloc(lines, (lineCount + 1) * sizeof(int));
        types = (int *)brealloc(types, (lineCount + 1) * sizeof(int));
        lines[lineCount] = lineNum;
        types[lineCount] = 1; // variable group
        lineCount++;
  
        // increment the count of variable groups
        vargpArray->numVarGps++;
      }
    }
    for (i = 0; i < numHomVar; i++)
    { // see if used
      found = 0;
      lineNum = Array->types[HOMVARIABLEGROUPTYPE].lineNumber[i];
      for (j = 0; j < lineCount && !found; j++)
        if (lines[j] == lineNum)
          found = 1;
 
      if (!found)
      { // add to lines
        lines = (int *)brealloc(lines, (lineCount + 1) * sizeof(int));
        types = (int *)brealloc(types, (lineCount + 1) * sizeof(int));
        lines[lineCount] = lineNum;
        types[lineCount] = 0; // hom variable group
        lineCount++;

        // increment the count of hom variable groups
        vargpArray->numHomGps++;
      }
    }

    // setup the total number of variable groups
    vargpArray->numGps = vargpArray->numVarGps + vargpArray->numHomGps;
 
    // allocate types & sizes
    vargpArray->types = (int *)bmalloc(vargpArray->numGps * sizeof(int)); 
    vargpArray->sizes = (int *)bmalloc(vargpArray->numGps * sizeof(int)); 
  
    // allocate varLoc & homLoc
    vargpArray->varLoc = (int *)bmalloc(numVar * sizeof(int));
    vargpArray->homLoc = (int *)bmalloc(numHomVar * sizeof(int));

    // count the total number of variables
    vargpArray->totalNumVars = vargpArray->numVarGps + numVar + numHomVar;
 
    // allocate groupNumber, nameSize, name
    vargpArray->groupNumber = (int *)bmalloc(vargpArray->totalNumVars * sizeof(int));
    vargpArray->nameSize = (int *)bmalloc(vargpArray->totalNumVars * sizeof(int));
    vargpArray->name = (char **)bmalloc(vargpArray->totalNumVars * sizeof(char *));

    // setup the new homogenizing variables
    for (i = 0; i < vargpArray->numVarGps; i++)
    { // setup nameSize[i] & name[i]
      vargpArray->nameSize[i] = snprintf(NULL, 0, "HOM_VAR_%d", i);
      vargpArray->name[i] = (char *)bmalloc((1 + vargpArray->nameSize[i]) * sizeof(char));
      sprintf(vargpArray->name[i], "HOM_VAR_%d", i);  
    }
    // move past these variables
    currVar = vargpArray->numVarGps;
 
    // setup the group numbers
    for (i = 0; i < lineCount; i++)
    { // find the minimum still left
      lineNum = INT_MAX;
      found = 0;
      for (j = 0; j < lineCount; j++)
        if (lines[j] < lineNum && lines[j] >= 0)
        { // store this line number
          lineNum = lines[j];
          found = j;
        }

      // setup type
      vargpArray->types[i] = types[found];
      // initialize size
      vargpArray->sizes[i] = 0;

      // setup the rest of the data
      if (types[found])
      { // variable group
        vargpArray->groupNumber[currVarGp] = i; // setup variable group for new homogenizing variable
        currVarGp++;

        for (j = 0; j < numVar; j++)
          if (Array->types[VARIABLEGROUPTYPE].lineNumber[j] == lineNum)
          { // setup the information about this variable
            vargpArray->varLoc[j] = currVar; // setup location
            vargpArray->groupNumber[currVar] = i; // variable group number
            // setup nameSize & name
            vargpArray->nameSize[currVar] = Array->types[VARIABLEGROUPTYPE].nameSize[j];
            vargpArray->name[currVar] = (char *)bmalloc((1 + vargpArray->nameSize[currVar]) * sizeof(char));
            strcpy(vargpArray->name[currVar], Array->types[VARIABLEGROUPTYPE].name[j]);
            // increment the current variable
            currVar++;
            // increment the size of this variable group
            vargpArray->sizes[i]++;
          }
      }
      else
      { // hom variable group
        for (j = 0; j < numHomVar; j++)
          if (Array->types[HOMVARIABLEGROUPTYPE].lineNumber[j] == lineNum)
          { // setup the information about this variable
            vargpArray->homLoc[j] = currVar; // setup location
            vargpArray->groupNumber[currVar] = i; // variable group number
            // setup nameSize & name
            vargpArray->nameSize[currVar] = Array->types[HOMVARIABLEGROUPTYPE].nameSize[j];
            vargpArray->name[currVar] = (char *)bmalloc((1 + vargpArray->nameSize[currVar]) * sizeof(char));
            strcpy(vargpArray->name[currVar], Array->types[HOMVARIABLEGROUPTYPE].name[j]);
            // increment the current variable
            currVar++;
            // increment the size of this variable group
            vargpArray->sizes[i]++;
          }
      }

      // kill off this line number
      lines[found] = -1;
    }
  
    // clear memory
    free(types);
    free(lines);
  }

  return;
} 


