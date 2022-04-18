// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "ppParse.h"

int isDefinedSubfunc(int *sfCall, int memoryLoc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - memoryLoc uses defined subfunction, 0-other*
* NOTES: sfCall - subfunction call number                       *
\***************************************************************/
{
  *sfCall = -(memoryLoc + 3);
  if (*sfCall >= 0)
    return 1;
  else
    return 0;
}

void setupDefinedSubfuncOps(expArrayOps *outOps, expArrayOps *inOps, int *memConv, int *firstFreeMemLoc, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: encode the operations for these operations             *
\***************************************************************/
{
  int i, numOps = inOps->numOps, numOldOps = outOps->numOps;
  int numNewOps = numOps + numOldOps;
  expOps tempOp, *tempPtrOp = NULL;

  // allocate new memory in outOps
  outOps->ops = (expOps *)brealloc(outOps->ops, numNewOps * sizeof(expOps));
  outOps->numOps = numNewOps;

  // copy over the operations
  for (i = 0; i < numOps; i++)
  { // copy the ith operation
    tempPtrOp = &inOps->ops[i];
    // setup the operation
    tempOp.op = tempPtrOp->op;
    // setup the input arguments
    tempOp.in[0] = memConv[tempPtrOp->in[0]];
    if (isUnary(tempOp.op))
      tempOp.in[1] = -1; // not used
    else
      tempOp.in[1] = memConv[tempPtrOp->in[1]];
    // setup output argument (should not be defined yet!)
    tempOp.memLoc = memConv[tempPtrOp->memLoc] = *firstFreeMemLoc;
    (*firstFreeMemLoc)++;
    // copy Op
    copyOp(&outOps->ops[i + numOldOps], &tempOp, zeroAddr, oneAddr);
    // NULL out tempPtrOp
    tempPtrOp = NULL;
  }

  return;
}

void setupDefinedSubfuncOperations(int *outNum, expArrayOps *funcOps, int SFcall, parseArray *Array, int zeroAddr, int oneAddr, subFuncData *defSFData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the proper operations for this subfunction       *
\***************************************************************/
{
  int i, currLoc, newSFnum; 
  int SFnum = Array->subfuncCalls[SFcall].subfuncLoc; // which subfunction to use
  int numSFVars = defSFData[SFnum].numVars, numSFConsts = defSFData[SFnum].numConsts, numSFNumbers = defSFData[SFnum].numNumbers;  
  int totalSFMem = defSFData[SFnum].totalMem;
  int *memConversion = (int *)bmalloc(totalSFMem * sizeof(int));

  // setup conversion
  currLoc = 0;
  for (i = 0; i < numSFVars; i++) // variable arguments
  { // determine if argument is also a defined subfunction
    if (isDefinedSubfunc(&newSFnum, Array->subfuncCalls[SFcall].arguments[i].argMemLoc))
    { // this argument is itself a defined subfunction - need to evaluate it first!
      setupDefinedSubfuncOperations(&memConversion[currLoc], funcOps, newSFnum, Array, zeroAddr, oneAddr, defSFData);    
    }
    else
    { // argument is a known memory location
      memConversion[currLoc] = Array->subfuncCalls[SFcall].arguments[i].argMemLoc; 
    }

    currLoc++;
  }
  // I
  memConversion[currLoc] = defSFData[SFnum].IAddr;
  currLoc++;
  // Pi
  memConversion[currLoc] = defSFData[SFnum].IAddr + 1;
  currLoc++;
  // constant arguments
  for (i = 0; i < numSFConsts; i++)
  { // determine if argument is also a defined subfunction
    if (isDefinedSubfunc(&newSFnum, Array->subfuncCalls[SFcall].arguments[i + numSFVars].argMemLoc))
    { // this argument is itself a defined subfunction - need to evaluate it first!
      setupDefinedSubfuncOperations(&memConversion[currLoc], funcOps, newSFnum, Array, zeroAddr, oneAddr, defSFData);
    }
    else
    { // argument is a known memory location
      memConversion[currLoc] = Array->subfuncCalls[SFcall].arguments[i + numSFVars].argMemLoc;
    }

    currLoc++;
  }
  // numbers
  for (i = 0; i < numSFNumbers; i++)
  {
    memConversion[currLoc] = defSFData[SFnum].numAddr[i];
    currLoc++;
  }
  // setup rest as unknown
  for (i = currLoc; i < totalSFMem; i++)
    memConversion[i] = -2;

  // copy function evaluations
  setupDefinedSubfuncOps(funcOps, &defSFData[SFnum].funcData, memConversion, &Array->firstFreeMemLoc, zeroAddr, oneAddr);

  // setup outNum
  *outNum = memConversion[defSFData[SFnum].funcData.finalMemLoc];

  // clear memory
  free(memConversion);

  return;
}

void addDefinedSubfuncOp(expArrayOps *funcOps, expOps *Op, parseArray *Array, int zeroAddr, int oneAddr, subFuncData *defSFData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add the operation to funcOps                           *
\***************************************************************/
{
  int isSubfunc0 = 0, isSubfunc1 = 0, sfCall0 = 0, sfCall1 = 0, unary = isUnary(Op->op), newLoc0 = 0, newLoc1 = 0;

  // determine if in[0] is a defined subfunction
  isSubfunc0 = isDefinedSubfunc(&sfCall0, Op->in[0]);
  if (!unary)
    isSubfunc1 = isDefinedSubfunc(&sfCall1, Op->in[1]);

  // error checking
  if (!isSubfunc0 && !isSubfunc1)
  { // error
    printf("ERROR: This operation does not utilize defined subfunctions!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  if (isSubfunc0)
  { // setup the operations needed to compute this value 
    setupDefinedSubfuncOperations(&newLoc0, funcOps, sfCall0, Array, zeroAddr, oneAddr, defSFData);
  }
  else
  { // leave memory location alone
    newLoc0 = Op->in[0];
  }

  if (isSubfunc1)
  { // setup the operations needed to compute this value 
    setupDefinedSubfuncOperations(&newLoc1, funcOps, sfCall1, Array, zeroAddr, oneAddr, defSFData);
  }
  else
  { // leave memory location alone
    newLoc1 = Op->in[1];
  }

  // add the input operation to funcOps
  addOp(funcOps, Op->memLoc, Op->op, newLoc0, newLoc1, zeroAddr, oneAddr);

  return;
}

void reorderFunctions(parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reorder functions: subFuncs first and then functions   *
\***************************************************************/
{
  int i, j, k, numFuncs = Array->memoryLoc.numFuncs, numSubfuncs = Array->memoryLoc.numSubfuncs, numDefStatements = Array->numDefiningStatements;
  int numOther = numDefStatements - numFuncs - numSubfuncs;
  int currFunc = 0;
  int *defOrder = (int *)bmalloc(numDefStatements * sizeof(int));
  defStatement *tempDefStatements = NULL;

  // initialize defOrder
  for (i = 0; i < numDefStatements; i++)
    defOrder[i] = i;

  // put all non-(sub)functions first
  for (i = 0; i < numOther; i++)
  { // find the next other defining statement
    k = -1;
    for (j = i; j < numDefStatements; j++)
      if (Array->definingStatements[defOrder[j]].lvalType != FUNCTIONTYPE && Array->definingStatements[defOrder[j]].lvalType != INLINESUBFUNCTIONTYPE)
      { // we have our next one
        k = j;
        break;
      }
    // copy i,..,k-1,k to k,i+1,...k-1
    currFunc = defOrder[k];
    for (j = k; j > i; j--)
      defOrder[j] = defOrder[j-1];
    defOrder[i] = currFunc;
  }

  // put all subfunctions next
  for (i = numOther; i < numSubfuncs + numOther; i++)
  { // find the next subfunction defining statement
    k = -1;
    for (j = i; j < numDefStatements; j++)
      if (Array->definingStatements[defOrder[j]].lvalType == INLINESUBFUNCTIONTYPE)
      { // we have our next one
        k = j;
        break;
      }
    // copy i,..,k-1,k to k,i+1,...k-1
    currFunc = defOrder[k];
    for (j = k; j > i; j--)
      defOrder[j] = defOrder[j-1];
    defOrder[i] = currFunc;
  }

  // order the functions by lvalLoc
  for (i = numSubfuncs + numOther; i < numDefStatements; i++)
  { // find the function of minimal degree
    k = i;
    for (j = i + 1; j < numDefStatements; j++)
    { // see if smaller
      if (Array->definingStatements[defOrder[k]].lvalLoc > Array->definingStatements[defOrder[j]].lvalLoc)
        k = j;
    }

    // copy i,..,k-1,k to k,i+1,...k-1
    currFunc = defOrder[k];
    for (j = k; j > i; j--)
      defOrder[j] = defOrder[j-1];
    defOrder[i] = currFunc;
  }

  // now that we have the order, we need to setup the definingStatements in this order
  tempDefStatements = Array->definingStatements;
  Array->definingStatements = (defStatement *)bmalloc(numDefStatements * sizeof(defStatement));
  for (i = 0; i < numDefStatements; i++)
  {
    copy_defStatement(&Array->definingStatements[i], &tempDefStatements[defOrder[i]]);
    clear_defStatement(&tempDefStatements[defOrder[i]]);
  }

  // clear memory
  free(tempDefStatements);
  free(defOrder);

  return;
}





