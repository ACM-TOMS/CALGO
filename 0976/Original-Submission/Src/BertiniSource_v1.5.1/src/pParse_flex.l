%{
// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "ppParse.h"
#include "pParse_bison.h"

#define yylex pParse_bisonlex


variablegroupArray *pVarGpArray = NULL;
definedSubfuncStack *pActiveSFData = NULL;
parseArray *pArray = NULL;

// processing line
char *pLine = NULL;
int pLineLen = 0, pCurrNumber = 0, pAddAllNums = 0;

extern int *pendOfFile;
extern int *plineNum;

/* %option nounput - do not warn about unput */
%}
%option nounput
%option bison-bridge

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
             append_str(&pLine, &pLineLen, yytext);
	     return pEND;
           }
{homvargp} { // hom variable group
             append_str(&pLine, &pLineLen, yytext);
             return pHOMVARGP;
           }
{vargp}    { // variable group
             append_str(&pLine, &pLineLen, yytext);
	     return pVARGP;
           }
{pathvar}  { // path variable
             append_str(&pLine, &pLineLen, yytext);
	     return pPATHVAR;
           }
{var}      { // variable
             append_str(&pLine, &pLineLen, yytext);
      	     return pVAR;
           }
{param}    { // parameter
             append_str(&pLine, &pLineLen, yytext);
	     return pPARAM;
           }
{const}    { // const
             append_str(&pLine, &pLineLen, yytext);
	     return pCONST;
           }
{func}     { // function
             append_str(&pLine, &pLineLen, yytext);
	     return pFUNC;
           }
{subfunc}  { // defined subfunction
             append_str(&pLine, &pLineLen, yytext);
             return pSUBFUNC;
           } 
{number}   { // number
             append_str(&pLine, &pLineLen, yytext);
	     pParse_bisonlval.name = strdup(yytext);
	     return pNUMBER;
           }
{sin}|{cos}|{tan}|{exp} { // sin, cos, tan, or exp
             append_str(&pLine, &pLineLen, yytext);
	     pParse_bisonlval.name = strdup(yytext);
	     return pSCE;
	   }
{pow}      { // pow 
             append_str(&pLine, &pLineLen, yytext);
             pParse_bisonlval.name = strdup(yytext);
             return pPOW;
           }
{I}        { // I
             append_str(&pLine, &pLineLen, yytext);
	     return pI;
           }
{Pi}       { // Pi
             append_str(&pLine, &pLineLen, yytext);
	     return pPi;
           }
{name}     { // name
             append_str(&pLine, &pLineLen, yytext);
             pParse_bisonlval.name = strdup(yytext);
	     return pNAME;
           }
{eol}      { // end of line
             append_str(&pLine, &pLineLen, yytext);
	     return pEOL;
           }
{operator} { // -+*/^(),=
             append_str(&pLine, &pLineLen, yytext);
	     return yytext[0];
           }
{white}    ; // ignore it
{other}    { // other symbol - error!
             append_str(&pLine, &pLineLen, yytext);
	     pParse_bisonerror("invalid character");
           }

%%

void clear_definedSubfuncData(definedSubfuncData *subFunc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear subFunc                                          *
\***************************************************************/
{
  if (subFunc->numArguments > 0)
  { // clear arguments
    free(subFunc->arguments);
  }

  // NULL out memory
  subFunc->arguments = NULL;
  subFunc->subfuncLoc = subFunc->numArguments = 0;

  return;
}

void increment_subfunc_parseArray(parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add a subfuncData element to Array                     *
\***************************************************************/
{
  int currFunc = Array->numSubfuncCalls;
  // increment number of subfunction calls
  Array->numSubfuncCalls++;
  // reallocate memory
  Array->subfuncCalls = (definedSubfuncData *)brealloc(Array->subfuncCalls, Array->numSubfuncCalls * sizeof(definedSubfuncData));
  // setup new subfunction
  Array->subfuncCalls[currFunc].subfuncLoc = 0;
  Array->subfuncCalls[currFunc].numArguments = 0;
  Array->subfuncCalls[currFunc].arguments = NULL;

  return;
}

void increment_definedSubfuncStack(definedSubfuncStack *sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add a subfuncData element to sfData                    *
\***************************************************************/
{ 
  int currFunc = sfData->numActiveSubfuncCalls;
  // increment number of subfunction calls
  sfData->numActiveSubfuncCalls++;
  // reallocate memory
  sfData->subfuncCalls = (definedSubfuncData *)brealloc(sfData->subfuncCalls, sfData->numActiveSubfuncCalls * sizeof(definedSubfuncData));
  // setup new subfunction
  sfData->subfuncCalls[currFunc].subfuncLoc = 0;
  sfData->subfuncCalls[currFunc].numArguments = 0;
  sfData->subfuncCalls[currFunc].arguments = NULL;
  
  return;
}

void remove_definedSubfuncDataStack(definedSubfuncStack *sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: delete a subfuncData element from sfData               *
\***************************************************************/
{ 
  int currFunc = sfData->numActiveSubfuncCalls - 1;

  if (currFunc >= 0)
  { // clear memory
    clear_definedSubfuncData(&sfData->subfuncCalls[currFunc]);
    // decrement number of subfunction calls
    sfData->numActiveSubfuncCalls--;
    // reallocate memory
    sfData->subfuncCalls = (definedSubfuncData *)brealloc(sfData->subfuncCalls, sfData->numActiveSubfuncCalls * sizeof(definedSubfuncData));
  }
  
  return;
}

void copy_argumentData(argumentData *Out, argumentData *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  Out->argMemLoc = In->argMemLoc;
  Out->argDerivExists = In->argDerivExists;
  Out->argDerivLoc = In->argDerivLoc;

  return;
}

void copy_definedSubfuncData(definedSubfuncData *Out, definedSubfuncData *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  int i;

  Out->subfuncLoc = In->subfuncLoc;
  Out->numArguments = In->numArguments;
  Out->arguments = (argumentData *)bmalloc(Out->numArguments * sizeof(argumentData));
  for (i = 0; i < Out->numArguments; i++)
    copy_argumentData(&Out->arguments[i], &In->arguments[i]);

  return;
}

void add_arg_definedSubfuncData(definedSubfuncData *sf, int argMemLoc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add an argument to sf                                  *
\***************************************************************/
{
  int currArg = sf->numArguments;
  // increment number of arguments
  sf->numArguments++;
  // reallocate memory
  sf->arguments = (argumentData *)brealloc(sf->arguments, sf->numArguments * sizeof(argumentData));
  // setup new argument
  sf->arguments[currArg].argMemLoc = argMemLoc;
  sf->arguments[currArg].argDerivExists = sf->arguments[currArg].argDerivLoc = 0;

  return;
}

void add_op_expArrayOps(expArrayOps *exp, int memLoc, char op, int in0, int in1)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add an operation to exp                                *
\***************************************************************/
{
  int currOp = exp->numOps;
  // increment number of operations
  exp->numOps++;
  // reallocate memory
  exp->ops = (expOps *)brealloc(exp->ops, exp->numOps * sizeof(expOps));
  // setup new operation
  exp->ops[currOp].memLoc = memLoc;
  exp->ops[currOp].op = op;
  exp->ops[currOp].in[0] = in0;
  exp->ops[currOp].in[1] = in1;

  return;
}

int lookup_number_type(int *loc, preProcessArrayType *Array, char *name)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - found, 0 - not found                       *
* NOTES: lookup a number in ppType based on pCurrNumber         *
\***************************************************************/
{
  int rV = 0, len = strlen(name);

  *loc = -1;

  if (len == Array->nameSize[pCurrNumber] && !strcmp(name, Array->name[pCurrNumber]))
  { // found a match!
    rV = 1;
    *loc = pCurrNumber;
    pCurrNumber++;
  }

  return rV;
}

int lookup_memLoc(int type, int loc, memoryLocations *memLoc, variablegroupArray *varGp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the memory location                            *
\***************************************************************/
{
  int mem = 0;

  if (type == DEFINEDSUBFUNCTIONTYPE)
  { // no memory location setup for this type!
    printf("ERROR: There are no memory locations for defined subfunctions!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  else if (type == FUNCTIONTYPE)
  { // lookup memory location
    mem = memLoc->funcStart + loc;
  }
  else if (type == VARIABLETYPE)
  { // lookup memory location
    mem = memLoc->varStart + loc;
  }
  else if (type == VARIABLEGROUPTYPE)
  { // lookup memory location
    mem = memLoc->varStart + varGp->varLoc[loc];
  }
  else if (type == HOMVARIABLEGROUPTYPE)
  { // lookup memory location
    mem = memLoc->varStart + varGp->homLoc[loc];
  }
  else if (type == PARAMETERTYPE)
  { // lookup memory location
    mem = memLoc->paramStart + loc;
  }
  else if (type == PATHVARIABLETYPE)
  { // lookup memory location
    mem = memLoc->pathvarStart + loc;
  }
  else if (type == CONSTANTTYPE)
  { // lookup memory location
    mem = memLoc->constStart + loc;
  }
  else if (type == INLINESUBFUNCTIONTYPE)
  { // lookup memory location
    mem = memLoc->subfuncStart + loc;
  }
  else if (type == NUMBERTYPE)
  { // lookup memory location
    mem = memLoc->numStart + loc;
  }

  return mem;
}

void initialize_defStatement(defStatement *dState, int type, int loc, memoryLocations *memLoc, variablegroupArray *varGp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize dState                                      *
\***************************************************************/
{ // setup data about lval
  dState->lvalType = type;
  dState->lvalLoc = loc;
  dState->lvalMemLoc = lookup_memLoc(type, loc, memLoc, varGp);

  // setup data about rval
  initialize_expArrayOps(&dState->rvalExp);
  dState->rvalExp.finalMemLoc = -2;

  return;
}

void copy_expOps(expOps *Out, expOps *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  Out->memLoc = In->memLoc;
  Out->op = In->op;
  Out->in[0] = In->in[0];
  Out->in[1] = In->in[1];

  return;
}

void copy_expArrayOps(expArrayOps *Out, expArrayOps *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  int i;

  Out->numOps = In->numOps;
  Out->finalMemLoc = In->finalMemLoc;
  Out->ops = (expOps *)bmalloc(Out->numOps * sizeof(expOps));
  for (i = 0; i < Out->numOps; i++)
    copy_expOps(&Out->ops[i], &In->ops[i]);

  return;
}

void copy_defStatement(defStatement *Out, defStatement *In)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  Out->lvalType = In->lvalType;
  Out->lvalLoc = In->lvalLoc;
  Out->lvalMemLoc = In->lvalMemLoc;

  copy_expArrayOps(&Out->rvalExp, &In->rvalExp);

  return;
}

void setup_memLoc(memoryLocations *memLoc, int *firstFree, preProcessArray *preProcArray, variablegroupArray *vargpArray, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup memLoc & firstFree using preProcArray            *
\***************************************************************/
{
  int currLoc = 0;

  if (userHom)
  { // variable
    memLoc->numVars = preProcArray->types[VARIABLETYPE].numType;
  }
  else
  { // variable_group & hom_variable_group
    memLoc->numVars = vargpArray->totalNumVars;
  }

  memLoc->varStart = currLoc;
  currLoc += memLoc->numVars;

  memLoc->numPathVars = preProcArray->types[PATHVARIABLETYPE].numType;
  memLoc->pathvarStart = currLoc;
  currLoc += memLoc->numPathVars;

  memLoc->numParams = preProcArray->types[PARAMETERTYPE].numType;
  memLoc->paramStart = currLoc;
  currLoc += memLoc->numParams;
  memLoc->paramDerivStart = currLoc;
  currLoc += memLoc->numParams * memLoc->numPathVars;

  memLoc->numFuncs = preProcArray->types[FUNCTIONTYPE].numType;
  memLoc->funcStart = currLoc;
  currLoc += memLoc->numFuncs;
  memLoc->funcDerivVStart = currLoc;
  currLoc += memLoc->numFuncs * memLoc->numVars;
  memLoc->funcDerivPStart = currLoc;
  currLoc += memLoc->numFuncs * memLoc->numParams;

  memLoc->numConsts = preProcArray->types[CONSTANTTYPE].numType;
  memLoc->constStart = currLoc;
  currLoc += memLoc->numConsts;

  memLoc->numNums = preProcArray->types[NUMBERTYPE].numType;
  memLoc->numStart = currLoc;
  currLoc += memLoc->numNums;

  memLoc->numSubfuncs = preProcArray->types[INLINESUBFUNCTIONTYPE].numType;
  memLoc->subfuncStart = currLoc;
  currLoc += memLoc->numSubfuncs;
  memLoc->subfuncDerivVStart = currLoc;
  currLoc += memLoc->numSubfuncs * memLoc->numVars;
  memLoc->subfuncDerivPStart = currLoc;
  currLoc += memLoc->numSubfuncs * memLoc->numParams;

  *firstFree = currLoc;

  return;
}

void undefine_preProcessArray_type(preProcessArrayType *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: undefined this type                                    *
\***************************************************************/
{
  int i;

  for (i = 0; i < Array->numType; i++)
    Array->isDefined[i] = 0;

  return;
}

void initialize_parseArray(parseArray *parArray, preProcessArray *preProcArray, variablegroupArray *vargpArray, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize parArray using preProcArray                 *
\***************************************************************/
{
  parArray->numDefiningStatements = 0;
  parArray->definingStatements = NULL;
  parArray->numSubfuncCalls = 0;
  parArray->subfuncCalls = NULL;

  // copy preprocess info
  copy_preProcessArray(&parArray->ppArray, preProcArray);

  // mark constants, parameters, inline subfunctions & functions as undefined
  undefine_preProcessArray_type(&parArray->ppArray.types[CONSTANTTYPE]);
  undefine_preProcessArray_type(&parArray->ppArray.types[PARAMETERTYPE]);
  undefine_preProcessArray_type(&parArray->ppArray.types[INLINESUBFUNCTIONTYPE]);
  undefine_preProcessArray_type(&parArray->ppArray.types[FUNCTIONTYPE]);

  // setup memory locations and first free memory location
  setup_memLoc(&parArray->memoryLoc, &parArray->firstFreeMemLoc, &parArray->ppArray, vargpArray, userHom);

  return;
}

void initialize_expArrayOps(expArrayOps *exp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize exp                                         *
\***************************************************************/
{ // NULL out memory
  exp->numOps = exp->finalMemLoc = 0;
  exp->ops = NULL;
  
  return;
}

void clear_expArrayOps(expArrayOps *exp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear exp                                              *
\***************************************************************/
{
  if (exp->numOps > 0)
  { // release memory
    free(exp->ops);
  }

  // NULL out memory
  exp->ops = NULL;
  exp->numOps = exp->finalMemLoc = 0;

  return;
}

void clear_defStatement(defStatement *dState)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear dState                                           *
\***************************************************************/
{
  // clear rval
  clear_expArrayOps(&dState->rvalExp);
  // clear lval
  dState->lvalType = dState->lvalLoc = dState->lvalMemLoc = 0;

  return;
}

void clear_parseArray(parseArray *parArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear parArray                                         *
\***************************************************************/
{
  int i;

  // clear defining statements
  if (parArray->numDefiningStatements > 0)
  {
    for (i = 0; i < parArray->numDefiningStatements; i++)
      clear_defStatement(&parArray->definingStatements[i]); 
    free(parArray->definingStatements);
  }

  // clear subfunction calls
  if (parArray->numSubfuncCalls > 0)
  {
    for (i = 0; i < parArray->numSubfuncCalls; i++)
      clear_definedSubfuncData(&parArray->subfuncCalls[i]);
    free(parArray->subfuncCalls);
  }

  // clear preprocess array
  clear_preProcessArray(&parArray->ppArray);

  // NULL out memory
  parArray->numDefiningStatements = parArray->numSubfuncCalls = parArray->firstFreeMemLoc = 0;
  parArray->definingStatements = NULL;
  parArray->subfuncCalls = NULL;

  return;
}

void load_subFuncData(parseArray *Array, subFuncData *sfData, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the defined subfunction data using IN            *
\***************************************************************/
{
  int i, j, loc, strSize = 0, strMem = 2;
  char *str = (char *)bmalloc(strMem * sizeof(char));

  // read in data
  fscanf(IN, "TotalMem: %d;\n", &sfData->totalMem);

  // variables
  fscanf(IN, "Variables: %d;\n", &sfData->numVars);
  // move past variable names
  for (i = 0; i < sfData->numVars; i++)
    scanRestOfLine(IN);

  // constants
  fscanf(IN, "Constants: %d;\n", &sfData->numConsts);
  // setup IAddr - first constant
  sfData->IAddr = Array->memoryLoc.constStart;
  // move past constant names
  for (i = 0; i < sfData->numConsts; i++)
    scanRestOfLine(IN);

  // numbers
  fscanf(IN, "Numbers: %d;\n", &sfData->numNumbers);
  // allocate memory
  sfData->numAddr = (int *)bmalloc(sfData->numNumbers * sizeof(int));
  // read in the numbers and determine memory location
  for (i = 0; i < sfData->numNumbers; i++)
  { // read in number
    strSize = 0;
    str[strSize] = fgetc(IN);
    while(str[strSize] != '\n')
    { // increment size
      strSize++;
      // verify enough memory
      if (strSize >= strMem)
      {
        strMem *= 2;
        str = (char *)brealloc(str, strMem * sizeof(char));
      }
      // read in next char
      str[strSize] = fgetc(IN);
    }
    // setup last char to '\0'
    str[strSize] = '\0';

    // find location of this number
    if (!lookup_preProcessArray_type(&loc, &Array->ppArray.types[NUMBERTYPE], str))
    { // ERROR!!!
      printf("ERROR: Unknown symbol %s when trying to setup defined subfunctions.\n", str);
      bexit(ERROR_INPUT_SYNTAX);
    }
    else
    { // setup memory location
      sfData->numAddr[i] = Array->memoryLoc.numStart + loc;
    }
  }

  // setup function information
  sfData->derivData = (expArrayOps *)bmalloc(sfData->numVars * sizeof(expArrayOps));
  fscanf(IN, "FunctionMem: %d;\n", &sfData->funcData.finalMemLoc);
  fscanf(IN, "FunctionDerivMem:"); 
  for (i = 0; i < sfData->numVars; i++)
    fscanf(IN, "%d", &sfData->derivData[i].finalMemLoc);

  fscanf(IN, ";\n");
  fscanf(IN, "FunctionOps: %d;\n", &sfData->funcData.numOps);
  sfData->funcData.ops = (expOps *)bmalloc(sfData->funcData.numOps * sizeof(expOps));
  for (i = 0; i < sfData->funcData.numOps; i++)
  { // read in ith operation
    readInOp(&sfData->funcData.ops[i], IN);
  }
  
  // seutp derivative information
  fscanf(IN, "FunctionDerivOps:");
  for (i = 0; i < sfData->numVars; i++)
  { // read in number of operations and allocate memory
    fscanf(IN, "%d", &sfData->derivData[i].numOps);
    sfData->derivData[i].ops = (expOps *)bmalloc(sfData->derivData[i].numOps * sizeof(expOps));
  }
  fscanf(IN, ";\n");
  for (i = 0; i < sfData->numVars; i++)
    for (j = 0; j < sfData->derivData[i].numOps; j++)
    { // read in jth operation for ith derivative
      readInOp(&sfData->derivData[i].ops[j], IN);
    }   

  return;
}

void setup_subFuncData(parseArray *Array, subFuncData **sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the defined subfunction data                     *
\***************************************************************/
{
  int i, size = 0, numSub = Array->ppArray.types[DEFINEDSUBFUNCTIONTYPE].numType;
  char *name = NULL;
  FILE *IN = NULL;  

  // allocate memory
  *sfData = (subFuncData *)bmalloc(numSub * sizeof(subFuncData));

  // read in each one
  for (i = 0; i < numSub; i++)
  { // open file
    size = 1 + snprintf(NULL, 0, "%s.func", Array->ppArray.types[DEFINEDSUBFUNCTIONTYPE].name[i]);
    name = (char *)brealloc(name, size * sizeof(char));
    sprintf(name, "%s.func", Array->ppArray.types[DEFINEDSUBFUNCTIONTYPE].name[i]);
    IN = fopen(name, "r");
    if (IN == NULL)
    { // error
      printf("ERROR: The file '%s' for the defined subfunction '%s' does not exist!\n", name, Array->ppArray.types[DEFINEDSUBFUNCTIONTYPE].name[i]);
      bexit(ERROR_INPUT_SYNTAX);
    }

    // read in the data
    load_subFuncData(Array, &(*sfData)[i], IN);

    // close IN
    fclose(IN);
  }

  // free memory
  free(name);

  return;
}

void clear_subFuncData(int numType, subFuncData **sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the defined subfunction data                     *
\***************************************************************/
{
  int i, j;

  for (i = 0; i < numType; i++)
  { // clear data
    for (j = 0; j < (*sfData)[i].numVars; j++)
      clear_expArrayOps(&(*sfData)[i].derivData[j]);
    free((*sfData)[i].derivData);
    clear_expArrayOps(&(*sfData)[i].funcData);
    free((*sfData)[i].numAddr);
  }

  // release memory
  free(*sfData);
  *sfData = NULL;

  return;
}

void verify_subFuncCalls(parseArray *parArray, subFuncData *sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify the subfunction calls are okay                  *
\***************************************************************/
{
  int i, loc;

  for (i = 0; i < parArray->numSubfuncCalls; i++)
  { // setup the subfunction location
    loc = parArray->subfuncCalls[i].subfuncLoc;
    // verify number of arguments
    if (parArray->subfuncCalls[i].numArguments != sfData[loc].numVars + sfData[loc].numConsts)
    { // argument list not the same
      printf("ERROR: Subfunction call number %d, involving '%s', has an incorrect number of arguments!\n", i+1, parArray->ppArray.types[DEFINEDSUBFUNCTIONTYPE].name[loc]);
      printf("It was expecting %d argument(s), but received %d argument(s).\n", sfData[loc].numVars + sfData[loc].numConsts, parArray->subfuncCalls[i].numArguments);
      bexit(ERROR_INPUT_SYNTAX);
    }
  }

  return;
}

void initialize_definedSubfuncStack(definedSubfuncStack *sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize sfData                                      *
\***************************************************************/
{
  sfData->numActiveSubfuncCalls = 0;
  sfData->subfuncCalls = NULL;

  return;
}

void clear_definedSubfuncStack(definedSubfuncStack *sfData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear sfData                                           *
\***************************************************************/
{
  int i;

  for (i = 0; i < sfData->numActiveSubfuncCalls; i++)
    clear_definedSubfuncData(&sfData->subfuncCalls[i]);
  free(sfData->subfuncCalls);
  sfData->numActiveSubfuncCalls = 0;
  sfData->subfuncCalls = NULL;

  return;
}

int processSystem(parseArray *parArray, subFuncData **sfData, preProcessArray *preProcArray, variablegroupArray *vargpArray, int userHom, int addAllNums, FILE *fp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: process the file                                       *
\***************************************************************/
{
  int rV = 0, endFile = 0;
  definedSubfuncStack activeSFData;

  pendOfFile = &endFile;

  // setup pArray
  pArray = parArray;

  // initialize global data for lines 
  pLine = NULL;
  pLineLen = 0;

  // setup pCurrNumber -- skip by '0' & '1' so start at memory location 2
  pCurrNumber = 2;

  // setup pAddAllNums
  pAddAllNums = addAllNums;

  // setup the data inside of parArray from preProcArray
  initialize_parseArray(parArray, preProcArray, vargpArray, userHom);

  // seutp pActiveSFData
  pActiveSFData = &activeSFData;

  // initialize pActiveSFData
  initialize_definedSubfuncStack(pActiveSFData);

  // setup pVarGpArray
  pVarGpArray = vargpArray;

  // setup line number
  plineNum = &yylineno;
  *plineNum = 1;

  // process the file
  yyin = fp;
  do
  { // clear out the line
    pLineLen = 1;
    pLine = (char *)brealloc(pLine, pLineLen * sizeof(char));
    pLine[0] = '\0';

    // process the line
    rV = pParse_bisonparse(); 

    // increment line number
    (*plineNum)++;

  } while (!rV && !endFile);

  // setup defined subfunction data
  setup_subFuncData(parArray, sfData);

  // verify that the subfunctions are called with correct arguments
  verify_subFuncCalls(parArray, *sfData);

  // clear memory
  yyin = NULL;
  pArray = NULL;
  clear_definedSubfuncStack(pActiveSFData);
  pActiveSFData = NULL;

  free(pLine);
  pLineLen = 0;

  return rV;
}

int yywrap()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches EOF                         *
\***************************************************************/
{ // parser is complete
  return 1;
}

int pParse_bisonerror(char *s)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches an error-should never occur *
\***************************************************************/
{ // error in parser
  if (!strcmp(s,"invalid character"))
    printf("ERROR: The '%s' is an %s in system statement %d.\n", yytext, s, yylineno);
  else
    printf("ERROR: There was an error (%s) near '%s' in system statement %d.\n", s, yytext, yylineno);
  errorStatement_pParse();

  return 0;
}

int isUnary(char op)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - unary operation, 0 - binary operation      *
* NOTES: determine if the opeartion is unary or binary          *
\***************************************************************/
{
  if (op == '=' || op == 'N' || op == 'S' || op == 'C' || op == 'X')
    return 1;
  else
    return 0;
}

void readInOp(expOps *op, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: read in the next operation from IN                     *
\***************************************************************/
{
  fscanf(IN, "%c %d %d", &op->op, &op->memLoc, &op->in[0]);
  if (!isUnary(op->op))
  { // read in second input
    fscanf(IN, "%d", &op->in[1]);
  }
  else
  { // set to -1
    op->in[1] = -1;
  }
  fscanf(IN, "\n");

  return;
}

void load_subFuncData_numbers(preProcessArray *Array, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the numbers for defined subfunction data in IN   *
\***************************************************************/
{
  int i, tempInt, loc, strSize = 0, strMem = 2;
  char *str = (char *)bmalloc(strMem * sizeof(char));

  // read in data
  fscanf(IN, "TotalMem: %d;\n", &tempInt);

  // variables
  fscanf(IN, "Variables: %d;\n", &tempInt);
  // move past variable names
  for (i = 0; i < tempInt; i++)
    scanRestOfLine(IN);

  // constants
  fscanf(IN, "Constants: %d;\n", &tempInt);
  // move past constant names
  for (i = 0; i < tempInt; i++)
    scanRestOfLine(IN);

  // numbers
  fscanf(IN, "Numbers: %d;\n", &tempInt);
  // read in the numbers
  for (i = 0; i < tempInt; i++)
  { // read in the number
    strSize = 0;
    str[strSize] = fgetc(IN);
    while(str[strSize] != '\n')
    { // increment size
      strSize++;
      // verify enough memory
      if (strSize >= strMem)
      {
        strMem *= 2;
        str = (char *)brealloc(str, strMem * sizeof(char));
      }
      // read in next char
      str[strSize] = fgetc(IN);
    }
    // setup last char to '\0'
    str[strSize] = '\0';

    // see if this number is known
    if (!lookup_preProcessArray_type(&loc, &Array->types[NUMBERTYPE], str))
    { // add it to the list as defined
      addEntry_preProcessArray_type(&Array->types[NUMBERTYPE], str, 1, 0);
    }
  }

  // release memory
  free(str);

  return;
}

void setup_subFuncData_numbers(preProcessArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the numbers for defined subfunction data         *
\***************************************************************/
{ 
  int i, size = 0, numSub = Array->types[DEFINEDSUBFUNCTIONTYPE].numType;
  char *name = NULL;
  FILE *IN = NULL;

  // loop over each one
  for (i = 0; i < numSub; i++)
  { // open file
    size = 1 + snprintf(NULL, 0, "%s.func", Array->types[DEFINEDSUBFUNCTIONTYPE].name[i]);
    name = (char *)brealloc(name, size * sizeof(char));
    sprintf(name, "%s.func", Array->types[DEFINEDSUBFUNCTIONTYPE].name[i]);
    IN = fopen(name, "r");
    if (IN == NULL)
    { // error
      printf("ERROR: The file '%s' for the defined subfunction '%s' does not exist!\n", name, Array->types[DEFINEDSUBFUNCTIONTYPE].name[i]);
      bexit(ERROR_INPUT_SYNTAX);
    } 

    // read in the data
    load_subFuncData_numbers(Array, IN);

    // close IN
    fclose(IN);
  }

  // free memory
  free(name);

  return;
}

void errorStatement_pParse()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: called when parser reaches an error-should never occur *
\***************************************************************/
{
  printf("Statement %d up to error: '%s'\n", yylineno, pLine);
  bexit(ERROR_INPUT_SYNTAX);

  return;
}

