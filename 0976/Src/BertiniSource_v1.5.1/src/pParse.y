// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

%{
#include "ppParse_old.h"

int pDeclarationType = -1;
int pDefiningType = -1, pTypeFound = -1, pLocFound = -1;
int pCurrDefiningStatement = -1;
extern variablegroupArray *pVarGpArray;
extern definedSubfuncStack *pActiveSFData;
extern parseArray *pArray;
extern int pAddAllNums, pCurrNumber; 

%}

%union 
{
  char *name;
  int memLoc;
} // the name of the item 

%token pEND pHOMVARGP pVARGP pPATHVAR pVAR pPARAM pCONST pFUNC pSUBFUNC pI pPi
%token <name> pNUMBER pSCE pPOW pNAME

%type <memLoc> expression elem_expression function_call

// operations listed in reverse order of precendence
%left '+' '-'
%left '*' '/'
%nonassoc pNEG // unary operation
%left '^' // i.e. -x^4 = -(x^4)

%%

line:	/* no token */
        { // we are at the end of the file
          *endOfFile = 1;
        }
|       pEND
        { // we are at the end of the file
          *endOfFile = 1;
        }
|	declaration declist ; // declares things
|	lval '=' expression // defines things
	{ // store the rval memory location
	  pArray->definingStatements[pCurrDefiningStatement].rvalExp.finalMemLoc = $3;
	}
;

lval: // acceptable lvals must either be a declared (but undefined) function, parameter or constant, or a name that is not already used to create an inline subfunction
	pNAME
	{ // determine if it is in pArray
     	  if (lookup_preProcessArray(&pTypeFound, &pLocFound, &pArray->ppArray, $1))
          { // found the item in the array - it must be a function that is not defined!
            if (!pArray->ppArray.types[pTypeFound].isDefined[pLocFound])
            { // declared but undefined
              if (pTypeFound == FUNCTIONTYPE || pTypeFound == PARAMETERTYPE || pTypeFound == CONSTANTTYPE || pTypeFound == INLINESUBFUNCTIONTYPE)
  	      { // it is a function, parameter, constant, or inline subfunction - mark as defined
	        pArray->ppArray.types[pTypeFound].isDefined[pLocFound] = 1;
		// setup the current defining statement
		pCurrDefiningStatement = pArray->numDefiningStatements;
	        // add to the number of defining statemenst
 	        pArray->numDefiningStatements++;
		// allocate memory
		pArray->definingStatements = (defStatement *)brealloc(pArray->definingStatements, pArray->numDefiningStatements * sizeof(defStatement));
		// initialize memory
		initialize_defStatement(&pArray->definingStatements[pCurrDefiningStatement], pTypeFound, pLocFound, &pArray->memoryLoc, pVarGpArray);
          	// setup what we are defining
	        pDefiningType = pTypeFound;
              }
              else
              { // not a function
                printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
    		errorStatement_pParse();
              }
            }
	    else
            { // trying to redefine!
              printf("ERROR: Trying to redefine %s in system statement %d.\n", $1, *lineNum);
    	      errorStatement_pParse();
            }
          } 
          else
          { // not found
            printf("ERROR: Unknown symbol %s in system statement %d.\n", $1, *lineNum);
            errorStatement_pParse();
          }

          // free memory
          free($1);
          $1 = NULL;
	}
|	pNUMBER
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
    	  errorStatement_pParse();
        }
|       pSCE
        { // invalid lval
          printf("ERROR: Invalid left-hand value (%s) in system statement %d.\n", $1, *lineNum);
    	  errorStatement_pParse();
        }
|	pI
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (I) in system statement %d.\n", *lineNum);
    	  errorStatement_pParse();
        }
|	pPi
	{ // invalid lval
	  printf("ERROR: Invalid left-hand value (Pi) in system statement %d.\n", *lineNum);
    	  errorStatement_pParse();
        }
;	

expression: // everything is done in the reduction phase
        expression '+' expression
        { // addition
          if ($1 == -2)
          { // just the second number
	    $$ = $3;
          }
	  else if ($3 == -2)
	  { // just the first number
	    $$ = $1;
	  }
	  else if ($1 == -1)
	  { // 1 + second number
	    if ($3 == -1)
            { // add 1 + 1 by adding an operation
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '+', pArray->memoryLoc.numStart + 1, pArray->memoryLoc.numStart + 1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
	    else
	    { // add 1 + second number by adding an operation
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '+', pArray->memoryLoc.numStart + 1, $3);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
	  else if ($3 == -1)
          { // first number + 1
            add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '+', $1, pArray->memoryLoc.numStart + 1);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
          }
	  else
	  { // add the two number
            add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '+', $1, $3);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
          }
        }
|	expression '-' expression 
        { // subtraction
          if ($1 == -2)
          { // negate the second number
            if ($3 == -2)
            { // answer is still 0
              $$ = $3;
            }
            else if ($3 == -1)
            { // negate 1 by adding an operation
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', pArray->memoryLoc.numStart + 1, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
            else
            { // negate the memory location by adding an operation
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', $3, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
          else if ($1 == -1)
          { // 1 - second number
	    if ($3 == -2)
            { // answer is just first number
              $$ = $1;
            }
	    else if ($3 == -1)
	    { // answer is 0 
	      $$ = -2;
	    }
 	    else
 	    { // perform the subtraction	
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '-', pArray->memoryLoc.numStart + 1, $3);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
    	    }
          }
          else
          { // first number is something
	    if ($3 == -2)
	    { // answer is just first number
	      $$ = $1;
	    }
	    else if ($3 == -1)
            { // perform the subtraction        
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '-', $1, pArray->memoryLoc.numStart + 1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
	    else
	    { // perform the subtraction
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '-', $1, $3);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
        }
|	expression '*' expression
        { // multiplication
          if ($1 == -2 || $3 == -2)
          { // multiplication by 0 - set equal to 0
	    $$ = -2;
          }
          else if ($1 == -1)
          { // multiplication by 1 - set equal to second number
            $$ = $3;
          }
          else if ($3 == -1)
          { // multiplication by 1 - set equal to first number
            $$ = $1;
          }
          else
          { // perform the multiplication
            add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '*', $1, $3);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
          }
        }
|	expression '/' expression 
	{ // division
          if ($3 == -2)
          { // division by 0!!!
	    printf("ERROR: Attempting to divide by 0 in system statement %d.\n", *lineNum);
            errorStatement_pParse();
          }
          else if ($3 == -1)
          { // division by 1 - set equal to numerator
            $$ = $1;
          }
          else
          { // perform the division
            add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '/', $1, $3);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
          }
	}
|      	expression '^' expression
	{ // exponetiation
	  if ($3 == -2)
	  { // exponent is 0 - set equal to 1
	    $$ = -1;
	  }
	  else if ($3 == -1)
	  { // exponent is 1 - set equal to the base
	    $$ = $1;
	  }
	  else
	  { // perform the exponentiation
	    add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '^', $1, $3);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
	  }
	}
|       '+' expression %prec pNEG 
	{ // leave expression alone
	  $$ = $2;
	};
|       '-' expression %prec pNEG
	{ // negate the expression
          if (pArray->memoryLoc.numStart <= $2 && $2 < pArray->memoryLoc.numStart + pArray->memoryLoc.numNums)
          { // expression is a number
            if (!pAddAllNums)
            { // need to negate this number
              if ($2 == -2)
              { // answer is still 0
                $$ = $2;
              }
              else if ($2 == -1)
              { // negate 1 by adding an operation
                add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', pArray->memoryLoc.numStart + 1, -1);
                $$ = pArray->firstFreeMemLoc;
                pArray->firstFreeMemLoc++;
              }
              else
              { // negate the memory location by adding an operation
                add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', $2, -1);
                $$ = pArray->firstFreeMemLoc;
                pArray->firstFreeMemLoc++;
              }
            }
            else
            { // already setup properly
              $$ = $2;
            } 
          }
          else
          { // work with the given memory location
  	    if ($2 == -2)
	    { // answer is still 0
	      $$ = $2;
	    }
	    else if ($2 == -1)
	    { // negate 1 by adding an operation
	      add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', pArray->memoryLoc.numStart + 1, -1);
	      $$ = pArray->firstFreeMemLoc;
	      pArray->firstFreeMemLoc++;	    
	    }
	    else
	    { // negate the memory location by adding an operation
	      add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'N', $2, -1);
	      $$ = pArray->firstFreeMemLoc;
	      pArray->firstFreeMemLoc++;	    
	    }
          }
  	}
|	'(' expression ')' 
	{ // leave expression alone
	  $$ = $2;
	}
|	elem_expression 
	{ // leave expression alone
	  $$ = $1;
 	}
|	function_call
	{ // leave expression alone
	  $$ = $1;
	}
;

elem_expression:
	pNUMBER
	{ // determine if this number has been used before
          if (pAddAllNums)
          { // we lookup based on the expected location
            if (!lookup_number_type(&pLocFound, &pArray->ppArray.types[NUMBERTYPE], $1))
     	    { // not found -- try again with the negative value
              char *str = (char *)bmalloc((strlen($1) + 2) * sizeof(char));
              sprintf(str, "-%s", $1);
              if (!lookup_number_type(&pLocFound, &pArray->ppArray.types[NUMBERTYPE], str))
              { // error since number and its negative are unknown!
                printf("ERROR: Unknown number %s in system statement %d.\n", $1, *lineNum);
                errorStatement_pParse();
              }
              free(str);
            }
          }
          else if (!lookup_preProcessArray_type(&pLocFound, &pArray->ppArray.types[NUMBERTYPE], $1))
 	  { // not found
            printf("ERROR: Unknown number %s in system statement %d.\n", $1, *lineNum);
            errorStatement_pParse();
          }

	  // compute memory location based on where it was found
	  $$ = lookup_memLoc(NUMBERTYPE, pLocFound, &pArray->memoryLoc, pVarGpArray);

	  // determine if either 0 or 1
	  if ($$ == pArray->memoryLoc.numStart)
	  { // == 0
	    $$ = -2;
   	  }
	  else if ($$ == pArray->memoryLoc.numStart + 1)
	  { // == 1
	    $$ = -1;
	  }

	  // free memory
	  free($1);
	  $1 = NULL;
	}
|	pNAME
	{ // determine what it is
          if (!lookup_preProcessArray(&pTypeFound, &pLocFound, &pArray->ppArray, $1))
          { // not found
            printf("ERROR: Attempting to use the unknown symbol '%s' in system statement %d.\n", $1, *lineNum);
    	    errorStatement_pParse();
          }

	  // verify that is has been defined
	  if (!pArray->ppArray.types[pTypeFound].isDefined[pLocFound])
          { // not defined
            printf("ERROR: Attempting to use the symbol '%s' before it has been defined in system statement %d.\n", $1, *lineNum);
    	    errorStatement_pParse();
          }

          // verify it is okay
          verify_dependent_name(pDefiningType, pTypeFound, $1, *lineNum);

          // verify not defining an expression based on itself
          if (pArray->ppArray.types[pTypeFound].lineNumber[pLocFound] == *lineNum)
          { // not defined
            printf("ERROR: Attempting to define the symbol '%s' in terms of itself in system statement %d.\n", $1, *lineNum);
            errorStatement_pParse();
          }

	  // lookup memory location
	  $$ = lookup_memLoc(pTypeFound, pLocFound, &pArray->memoryLoc, pVarGpArray);

          // free memory
          free($1);
	  $1 = NULL;
	}
|	pI 
	{ // return the memory location of I - the first constant 
	  $$ = pArray->memoryLoc.constStart;
	}
|	pPi 
	{ // return the memory location of Pi - the second constant 
	  $$ = pArray->memoryLoc.constStart + 1;
	}
;

function_call:
	pSCE '(' expression ')' 
	{ // determine which of sin, cos, tan, or exp it is
	  if (!strcmp($1, "sin"))
	  { // sin(expression)
	    if ($3 == -2)
	    { // sin(0) = 0
	      $$ = -2;
	    }
	    else if ($3 == -1)
	    { // sin(1)
	      add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'S', pArray->memoryLoc.numStart + 1, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
	    }
	    else
            { // sin(expression)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'S', $3, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
	  else if (!strcmp($1, "cos"))
          { // cos(expression)
            if ($3 == -2)
            { // cos(0) = 1
              $$ = -1;
            }
            else if ($3 == -1)
            { // cos(1)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'C', pArray->memoryLoc.numStart + 1, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
            else
            { // cos(expression)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'C', $3, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
          else if (!strcmp($1, "tan"))
          { // tan(expression)
            if ($3 == -2)
            { // tan(0) = 0
              $$ = -2;
            }
            else if ($3 == -1)
            { // tan(1) = sin(1)/cos(1)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'S', pArray->memoryLoc.numStart + 1, -1);
              pArray->firstFreeMemLoc++;
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'C', pArray->memoryLoc.numStart + 1, -1);
              pArray->firstFreeMemLoc++;
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '/', pArray->firstFreeMemLoc - 2, pArray->firstFreeMemLoc - 1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
            else
            { // tan(expression) = sin(expression)/cos(expression)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'S', $3, -1);
              pArray->firstFreeMemLoc++;
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'C', $3, -1);
              pArray->firstFreeMemLoc++;
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '/', pArray->firstFreeMemLoc - 2, pArray->firstFreeMemLoc - 1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
          else if (!strcmp($1, "exp"))
          { // exp(expression)
            if ($3 == -2)
            { // exp(0) = 1
              $$ = -1;
            }
            else if ($3 == -1)
            { // exp(1)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'X', pArray->memoryLoc.numStart + 1, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
            else
            { // exp(expression)
              add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, 'X', $3, -1);
              $$ = pArray->firstFreeMemLoc;
              pArray->firstFreeMemLoc++;
            }
          }
	  else
	  { // error
            printf("ERROR: Invalid name of a defined subfunction (%s) in system statement %d.\n", $1, *lineNum);
            errorStatement_pParse();
	  }

          // free memory
          free($1);
          $1 = NULL;
	}
|       pPOW '(' expression ',' expression ')'
        { // pow -> (expression) ^ (expression)
          if ($5 == -2)
          { // exponent is 0 - set equal to 1
            $$ = -1;
          }
          else if ($5 == -1)
          { // exponent is 1 - set equal to the base
            $$ = $3;
          }
          else
          { // perform the exponentiation
            add_op_expArrayOps(&pArray->definingStatements[pCurrDefiningStatement].rvalExp, pArray->firstFreeMemLoc, '^', $3, $5);
            $$ = pArray->firstFreeMemLoc;
            pArray->firstFreeMemLoc++;
          }

          // free memory
          free($1);
          $1 = NULL;
        }
|       pNAME '(' arguments ')'
        { // determine if this corresponds to a defined subfunction
          if (!lookup_preProcessArray_type(&pLocFound, &pArray->ppArray.types[DEFINEDSUBFUNCTIONTYPE], $1))
          { // not found
            printf("ERROR: Invalid name of a defined subfunction (%s) in system statement %d.\n", $1, *lineNum);
    	    errorStatement_pParse();
          }

	  // setup the actual subfunction call
	  increment_subfunc_parseArray(pArray);
	  copy_definedSubfuncData(&pArray->subfuncCalls[pArray->numSubfuncCalls - 1], &pActiveSFData->subfuncCalls[pActiveSFData->numActiveSubfuncCalls - 1]);
	  pArray->subfuncCalls[pArray->numSubfuncCalls - 1].subfuncLoc = pLocFound;

	  // this one is no longer active
	  remove_definedSubfuncDataStack(pActiveSFData);

	  // setup the return memory location - start counting at -3
	  $$ = -pArray->numSubfuncCalls - 2;

          // free memory
          free($1);
	  $1 = NULL;
        }

arguments: // read in the arguments for first to last
        arguments ',' expression
        { // setup the next argument
          add_arg_definedSubfuncData(&pActiveSFData->subfuncCalls[pActiveSFData->numActiveSubfuncCalls - 1], $3);
        }
|	expression
	{ // add an active subfunction call
	  increment_definedSubfuncStack(pActiveSFData);

          // setup the first argument
	  add_arg_definedSubfuncData(&pActiveSFData->subfuncCalls[pActiveSFData->numActiveSubfuncCalls - 1], $1);
	}
;

declaration:
	pHOMVARGP
	{ // setup pDeclarationType
	  pDeclarationType = HOMVARIABLEGROUPTYPE;
	}
|	pVARGP
        { // setup pDeclarationType
          pDeclarationType = VARIABLEGROUPTYPE;
	}
|	pPATHVAR
        { // setup pDeclarationType
          pDeclarationType = PATHVARIABLETYPE;
        }
|	pVAR
        { // setup pDeclarationType
          pDeclarationType = VARIABLETYPE;
        }
|	pPARAM
        { // setup pDeclarationType
          pDeclarationType = PARAMETERTYPE;
        }
|	pCONST
        { // setup pDeclarationType
          pDeclarationType = CONSTANTTYPE;
        }
|	pFUNC
        { // setup pDeclarationType
          pDeclarationType = FUNCTIONTYPE;
        }
|	pSUBFUNC
        { // setup pDeclarationType
          pDeclarationType = DEFINEDSUBFUNCTIONTYPE;
        }
;

declist: // read in the list starting at the left
	declist ',' pNAME
	{ // determine if it is already in the list
          if (!lookup_preProcessArray_type(&pLocFound, &pArray->ppArray.types[pDeclarationType], $3))
          { // not found!
            printf("ERROR: Unknown symbol %s in system statement %d.\n", $3, *lineNum);
    	    errorStatement_pParse();
          }

          // free memory
          free($3);
	  $3 = NULL;
	}
|	pNAME
	{ // determine if it is already in the list
          if (!lookup_preProcessArray_type(&pLocFound, &pArray->ppArray.types[pDeclarationType], $1))
          { // not found!
            printf("ERROR: Unknown symbol %s in system statement %d.\n", $1, *lineNum);
    	    errorStatement_pParse();
          }

          // free memory
          free($1);
	  $1 = NULL;
	}
;

%%


