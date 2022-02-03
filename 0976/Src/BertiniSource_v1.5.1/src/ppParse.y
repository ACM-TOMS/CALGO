// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

%{
#include "ppParse_old.h"

int declarationType = -1, isDefined = 0, lastExpWasNumber = 0;
int definingType = -1, typeFound = -1, locFound = -1;
extern int ppUsingParamHomotopy, ppAddAllNums;
extern preProcessArray *ppArray;

%}

%union 
{
  char *name;
} // the name of the item 

%token ppEND ppHOMVARGP ppVARGP ppPATHVAR ppVAR ppPARAM ppCONST ppFUNC ppSUBFUNC ppI ppPi
%token <name> ppNUMBER ppSCE ppPOW ppNAME

// operations listed in reverse order of precendence
%left '+' '-' 
%left '*' '/'
%nonassoc ppNEG // unary operation
%right '^' // i.e. -x^4 = -(x^4)

%%

line:	/* no token */
        { // we are at the end of the file
          *endOfFile = 1;
          lastExpWasNumber = 0;
        }
|       ppEND
        { // we are at the end of the file
          *endOfFile = 1;
          lastExpWasNumber = 0;
        }
|	declaration declist 
        { // declares things
          lastExpWasNumber = 0;
        }
|	lval '=' expression 
        { // defines things
          lastExpWasNumber = 0;
        }
;

lval: // acceptable lvals must either be a declared (but undefined) function, parameter or constant, or a name that is not already used to create an inline subfunction
	ppNAME
	{ // determine if it is in ppArray
     	  if (lookup_preProcessArray(&typeFound, &locFound, ppArray, $1))
          { // found the item in the array - it must be a function that is not defined!
            if (!ppArray->types[typeFound].isDefined[locFound])
            { // declared but undefined
              if (typeFound == FUNCTIONTYPE || typeFound == PARAMETERTYPE || typeFound == CONSTANTTYPE)
  	      { // it is a function, parameter or constant - mark as defined
	        ppArray->types[typeFound].isDefined[locFound] = 1;
              }
              else
              { // not a function
                printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
                errorStatement();
              }
            }
	    else
            { // trying to redefine!
              printf("ERROR: Trying to redefine %s in system statement %d.\n", $1, *lineNum);
              errorStatement();
            }
          } 
          else
          { // add as an inline subfunction that is defined here
	    typeFound = INLINESUBFUNCTIONTYPE;
	    addEntry_preProcessArray_type(&ppArray->types[INLINESUBFUNCTIONTYPE], $1, 1, *lineNum);
          }
          // setup what we are defining
          definingType = typeFound;
          lastExpWasNumber = 0;

          // free memory
          free($1);
          $1 = NULL;
	}
|	ppNUMBER
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
          errorStatement();
          lastExpWasNumber = 1;
        }
|       ppSCE
        { // invalid lval
          printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
|	ppI
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (I) in system statement %d.\n", *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
|	ppPi
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (Pi) in system statement %d.\n", *lineNum);
    	  errorStatement();
          lastExpWasNumber = 0;
        }
;	

expression: // everything is done in the reduction phase
	expression '+' expression 
        {
          lastExpWasNumber = 0;
        }
|	expression '-' expression  
        {
          lastExpWasNumber = 0;
        }
|	expression '*' expression 
        {
          lastExpWasNumber = 0;
        }
|	expression '/' expression 
        {
          lastExpWasNumber = 0;
        }
|      	expression '^' expression 
        {
          lastExpWasNumber = 0;
        }
|       '+' expression %prec ppNEG 
        {
          lastExpWasNumber = 0;
        }
|       '-' expression %prec ppNEG 
        {
          if (ppAddAllNums && lastExpWasNumber)
          { // negate last number entered
            negate_lastEntry(&ppArray->types[NUMBERTYPE]);
          }
          lastExpWasNumber = 0;
        }
|	'(' expression ')' 
        {
          lastExpWasNumber = 0;
        }
|	elem_expression ;
|	function_call 
        {
          lastExpWasNumber = 0;
        }
;

elem_expression:
	ppNUMBER
	{ // determine if this number has been used before
          if (ppAddAllNums || !lookup_preProcessArray_type(&locFound, &ppArray->types[NUMBERTYPE], $1))
          { // add it to the list as defined
      	    addEntry_preProcessArray_type(&ppArray->types[NUMBERTYPE], $1, 1, *lineNum);
          }
	  // free memory
	  free($1);
	  $1 = NULL;
          lastExpWasNumber = 1;
	}
|	ppNAME
	{ // determine what it is
          if (!lookup_preProcessArray(&typeFound, &locFound, ppArray, $1))
          { // not found
            printf("ERROR: Attempting to use the unknown symbol '%s' in system statement %d.\n", $1, *lineNum);
    	    errorStatement();
          }

	  // verify that is has been defined
	  if (!ppArray->types[typeFound].isDefined[locFound])
          { // not defined
            printf("ERROR: Attempting to use the symbol '%s' before it has been defined in system statement %d.\n", $1, *lineNum);
    	    errorStatement();
          }

          // verify it is okay
          verify_dependent_name(definingType, typeFound, $1, *lineNum);

          // verify not defining an expression based on itself
          if (ppArray->types[typeFound].lineNumber[locFound] == *lineNum)
          { // not defined
            printf("ERROR: Attempting to define the symbol '%s' in terms of itself in system statement %d.\n", $1, *lineNum);
            errorStatement();
          }

          // free memory
          free($1);
	  $1 = NULL;
          lastExpWasNumber = 0;
	}
|	ppI 
        { // I -- ignore
         lastExpWasNumber = 0;
        }
|	ppPi 
        { // Pi -- ignore
          lastExpWasNumber = 0;
        }
;

function_call:
	ppSCE '(' expression ')' 
        { // sin, cos, tan, or exp
          lastExpWasNumber = 0;

          // free memory
          free($1);
          $1 = NULL;
        }
|       ppPOW '(' expression ',' expression ')'
        { // pow
          lastExpWasNumber = 0;

          // free memory
          free($1);
	  $1 = NULL;
        }
|       ppNAME '(' arguments ')'
        { // determine if this corresponds to a defined subfunction
          if (!lookup_preProcessArray_type(&locFound, &ppArray->types[DEFINEDSUBFUNCTIONTYPE], $1))
          { // not found
            printf("ERROR: Invalid name of a defined subfunction (%s) in system statement %d.\n", $1, *lineNum);
    	    errorStatement();
          }

          // verify it is okay
          verify_dependent_name(definingType, DEFINEDSUBFUNCTIONTYPE, $1, *lineNum);

          // free memory
          free($1);
	  $1 = NULL;
          lastExpWasNumber = 0;
        }

arguments:
	arguments ',' expression 
        {
          lastExpWasNumber = 0;
        } 
|       expression 
        {
          lastExpWasNumber = 0;
        } 
;

declaration:
	ppHOMVARGP
	{ // setup declarationType
	  declarationType = HOMVARIABLEGROUPTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
	}
|	ppVARGP
        { // setup declarationType
          declarationType = VARIABLEGROUPTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
	}
|	ppPATHVAR
        { // setup declarationType
          declarationType = PATHVARIABLETYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
|	ppVAR
        { // setup declarationType
          declarationType = VARIABLETYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
|	ppPARAM
        { // setup declarationType
          declarationType = PARAMETERTYPE;
	  isDefined = ppUsingParamHomotopy;
          lastExpWasNumber = 0;
        }
|	ppCONST
        { // setup declarationType
          declarationType = CONSTANTTYPE;
	  isDefined = 0;
          lastExpWasNumber = 0;
        }
|	ppFUNC
        { // setup declarationType
          declarationType = FUNCTIONTYPE;
	  isDefined = 0;
          lastExpWasNumber = 0;
        }
|	ppSUBFUNC
        { // setup declarationType
          declarationType = DEFINEDSUBFUNCTIONTYPE;
	  isDefined = 1;
          lastExpWasNumber = 0;
        }
;

declist: // read in the list starting at the left
	declist ',' ppNAME
	{ // determine if it is already in the list
          if (lookup_preProcessArray(&typeFound, &locFound, ppArray, $3))
          { // redefining!
            printf("ERROR: Trying to redefine %s in system statement %d.\n", $3, *lineNum);
    	    errorStatement();
          }
          else
          { // add to list
    	    addEntry_preProcessArray_type(&ppArray->types[declarationType], $3, isDefined, *lineNum);
          }

          // free memory
          free($3);
	  $3 = NULL;
          lastExpWasNumber = 0;
	}
|	ppNAME
	{ // determine if it is already in the list
          if (lookup_preProcessArray(&typeFound, &locFound, ppArray, $1))
          { // redefining!
            printf("ERROR: Trying to redefine %s in system statement %d.\n", $1, *lineNum);
    	    errorStatement();
          }
          else
          { // add to list
            addEntry_preProcessArray_type(&ppArray->types[declarationType], $1, isDefined, *lineNum);
          }

          // free memory
          free($1);
	  $1 = NULL;
          lastExpWasNumber = 0;
	}
;

%%

void verify_dependent_name(int defineType, int currType, char *name, int line)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify that 'defineType' can depend upon 'currType'    *
\***************************************************************/
{
  if (defineType == CONSTANTTYPE)
  { // must depend on other constants
    if (!(currType == CONSTANTTYPE || currType == NUMBERTYPE))
    { // not a constant
      printf("ERROR: Attempting to use the symbol '%s' when defining a constant in system statement %d.\nConstants must be defined using other constants or numbers.\n", name, line);
      errorStatement();
    }
  }
  else if (defineType == PARAMETERTYPE)
  { // must only depend on constants, pathvariables, parameters, or inline subfunction
    if (!(currType == PATHVARIABLETYPE || currType == CONSTANTTYPE || currType == NUMBERTYPE))
    { // not valid
      printf("ERROR: Attempting to use the symbol '%s' when defining a parameter in system statement %d.\nParameters must be defined using pathvariables, constants, or numbers.\n", name, line);
      errorStatement();
    }
  }
  else if (defineType == FUNCTIONTYPE)
  { // cannot depend on path variables or other functions
    if (currType == PATHVARIABLETYPE || currType == FUNCTIONTYPE)
    { // is a pathvariable
      printf("ERROR: Attempting to use the symbol '%s' when defining a function in system statement %d.\nFunctions cannot directly depend upon path variables or other functions.\n", name, line);
      errorStatement();
    }
  }

  return;
}

