// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "diff.h"
#include "ppParse.h"

// This file uses forward AD to find the derivatives in a way that works with deflation
// work variable by variable, rather than function by function

int checkExp_op(func_ops *op, comp_d *mem, int oneAddr, int allowNegExp);
void performOp(func_ops *op, comp_d *mem);
void copyToOp(func_ops *fop, int memLoc, char op, int in0, int in1);

void diff_sys_fwd_arr(systemStruct *sys, char *fileName);
void diff_forward_pathVars(systemStruct *sys, int currParam, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_forward_vars_subfuncs(systemStruct *sys, int currSubFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_forward_vars_funcs(systemStruct *sys, int currFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_forward_params_subfunc(systemStruct *sys, int currSubFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_forward_params_func(systemStruct *sys, int currFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);

void diff_deflatable(int putIntoOrder, int allowNegExp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the derivative of 'finalFile.out' in a way that is*
* 'friendly' to deflation using forward AD                      *
\***************************************************************/
{ 
  systemStruct sys;
  num_t *nums = NULL;  

  // read in the system
  setupSystemStructure(&sys, "finalFile.out", putIntoOrder);

  // setup nums - double precision is okay
  setupNums(&nums, sys.numNumbers, 64, 0);

  // check exponentiation using sys
  checkExp(&sys, nums, allowNegExp);

  // differentiate and create arr.out
  diff_sys_array(&sys, "arr.out");

  // differentiate the system and create arr.out using forward AD
//  diff_sys_fwd_arr(&sys, "arr.out");

  // clear nums
  clearNums(&nums, sys.numNumbers);

  // clear sys
  clearSystemStructure(&sys);

  return;
}

void diff_sys_fwd_arr(systemStruct *sys, char *fileName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the derivative of sys & prints SLP to 'fileName' *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{ // seutp OUT
  FILE *OUT = fopen(fileName, "w");
  if (OUT == NULL)
  {
    printf("ERROR: '%s' is an invalid name!\n", fileName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  int i, j, k, end, endUpdate, endParams, endFn, endPDerivs, endJvEval, totalF = sys->numSubfuncs + sys->numFuncs;
  int *startSubfuncs = (int *)bmalloc(sys->numSubfuncs * sizeof(int)), *endSubfuncs = (int *)bmalloc(sys->numSubfuncs * sizeof(int));
  int *startFuncs = (int *)bmalloc(sys->numFuncs * sizeof(int)), *endFuncs = (int *)bmalloc(sys->numFuncs * sizeof(int));
  int *startJvsub = (int *)bmalloc(sys->numSubfuncs * sizeof(int)), *endJvsub = (int *)bmalloc(sys->numSubfuncs * sizeof(int));
  int *startJv = (int *)bmalloc(sys->numFuncs * sizeof(int)), *endJv = (int *)bmalloc(sys->numFuncs * sizeof(int));
  int totalOpCount, *memLoc = NULL, *derivAddr = NULL;
  funcStruct **derivs = NULL;

  // count the total number of things that we need the derivative of
  totalOpCount = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
    totalOpCount += sys->params[i].num_ops;
  for (i = 0; i < sys->numSubfuncs; i++)
    totalOpCount += sys->subFuncs[i].num_ops;
  for (i = 0; i < sys->numFuncs; i++)
    totalOpCount += sys->funcs[i].num_ops;

  // allocate for memLoc & derivAddr
  memLoc = (int *)bmalloc(totalOpCount * sizeof(int));
  derivAddr = (int *)bmalloc(totalOpCount * sizeof(int));

  // print the update instructions
  end = endUpdate = printOps(OUT, sys->updateOps, sys->numUpdate);

  // print the parameter instructions
  for (i = 0; i < sys->numParams; i++)
  { // print the operations for the ith parameter
    end += printOps(OUT, sys->params[i].ops, sys->params[i].num_ops);
  }
  // store the location at the end of the parameters
  endParams = end;

  // print the subfunctions & functions (in the order they were defined)
  for (i = 0; i < totalF; i++)
  { // find what is the ith subfunc/function
    for (j = 0; j < totalF; j++)
      if (j < sys->numSubfuncs && sys->subFuncOrder[j] == i)
      { // print the jth subfunction
        startSubfuncs[j] = end; // where this subfunction start
        end += printOps(OUT, sys->subFuncs[j].ops, sys->subFuncs[j].num_ops);
        endSubfuncs[j] = end; // where this subfunction ends
      }
      else if (j >= sys->numSubfuncs && sys->funcOrder[j - sys->numSubfuncs] == i)
      { // print the (j - numSubfuncs)th function
        startFuncs[j - sys->numSubfuncs] = end; // where this function start
        end += printOps(OUT, sys->funcs[j - sys->numSubfuncs].ops, sys->funcs[j - sys->numSubfuncs].num_ops);
        endFuncs[j - sys->numSubfuncs] = end; // where this function ends
      }
  }
  // store the location at the end of the functions
  endFn = end;

  ////// PARAMETERS w.r.t. PATH VARIABLES //////

  // setup memory for holding the derivatives for each parameter w.r.t. each path variable
  derivs = (funcStruct **)bmalloc(sys->numParams * sizeof(funcStruct *));
  for (i = 0; i < sys->numParams; i++)
    derivs[i] = (funcStruct *)bmalloc(sys->numPathVars * sizeof(funcStruct));

  // compute the derivatives for each parameter w.r.t. each path variable
  for (i = 0; i < sys->numPathVars; i++)
  { // initialize the data
    initialize_memLoc_derivAddr(sys, memLoc, derivAddr, totalOpCount);
    for (j = 0; j < sys->numParams; j++)
    { // find the partial derivatives for the jth parameter w.r.t. ith path variable
      diff_sys_pathVar(&derivs[j][i], sys, j, i, memLoc, derivAddr, totalOpCount);
    }
  }

  // print the derivatives for each parameter w.r.t each path variable
  for (i = 0; i < sys->numParams; i++)
  {
    for (j = 0; j < sys->numPathVars; j++)
    { // print derivs
      end += printOps(OUT, derivs[i][j].ops, derivs[i][j].num_ops);

      // clear derivs
      free(derivs[i][j].ops);
      derivs[i][j].ops = NULL;
      derivs[i][j].num_ops = 0;
    }
    // clear derivs[i]
    free(derivs[i]);
  }
  // clear derivs
  free(derivs);

  // store the location at the end of the parameter derivatives
  endPDerivs = end; 

  ////// FUNCTIONS w.r.t. VARIABLES //////

  // setup memory for holding the derivatives for each (sub)function w.r.t. each variable
  derivs = (funcStruct **)bmalloc(totalF * sizeof(funcStruct *));
  for (i = 0; i < totalF; i++)
    derivs[i] = (funcStruct *)bmalloc(sys->numVars * sizeof(funcStruct));

  // compute the derivatives for each (sub)function w.r.t. each variable
  for (i = 0; i < sys->numVars; i++)
  { // initialize the data
    initialize_memLoc_derivAddr(sys, memLoc, derivAddr, totalOpCount);
    for (j = 0; j < totalF; j++)
    { // find the partial derivatives for the jth (sub)function w.r.t. ith variable
      for (k = 0; k < totalF; k++)
        if (k < sys->numSubfuncs && sys->subFuncOrder[k] == j)
          diff_sys_var_subfunc(&derivs[k][i], sys, k, i, memLoc, derivAddr, totalOpCount);
        else if (k >= sys->numSubfuncs && sys->funcOrder[k - sys->numSubfuncs] == j)
          diff_sys_var_func(&derivs[k][i], sys, k - sys->numSubfuncs, i, memLoc, derivAddr, totalOpCount);
    }
  }

  // print the derivatives for each (sub)function w.r.t each variable
  for (i = 0; i < totalF; i++)
  { // find the (sub)function corresponding to the ith one
    for (j = 0; j < totalF; j++)
      if (j < sys->numSubfuncs && sys->subFuncOrder[j] == i)
      { // jth subfunction
        startJvsub[j] = end; // where this deriv starts
        for (k = 0; k < sys->numVars; k++)
        { // print derivs
          end += printOps(OUT, derivs[j][k].ops, derivs[j][k].num_ops);

          // clear derivs[j][k]
          free(derivs[j][k].ops);
          derivs[j][k].ops = NULL;
          derivs[j][k].num_ops = 0;
        }
        endJvsub[j] = end; // where this deriv ends

        // clear derivs[j]
        free(derivs[j]);
      }
      else if (j >= sys->numSubfuncs && sys->funcOrder[j - sys->numSubfuncs] == i)
      { // jth function
        startJv[j - sys->numSubfuncs] = end; // where this deriv starts
        for (k = 0; k < sys->numVars; k++)
        { // print derivs
          end += printOps(OUT, derivs[j][k].ops, derivs[j][k].num_ops);

          // clear derivs[j][k]
          free(derivs[j][k].ops);
          derivs[j][k].ops = NULL;
          derivs[j][k].num_ops = 0;
        }
        endJv[j - sys->numSubfuncs] = end; // where this deriv ends

        // clear derivs[j]
        free(derivs[j]);
      }
  }
  // clear derivs
  free(derivs);

  // store the location at the end of the variable derivatives
  endJvEval = end;

  ////// FUNCTIONS w.r.t. PARAMETERS //////

  // setup memory for holding the derivatives for each (sub)function w.r.t. each parameter
  derivs = (funcStruct **)bmalloc(totalF * sizeof(funcStruct *));
  for (i = 0; i < totalF; i++)
    derivs[i] = (funcStruct *)bmalloc(sys->numParams * sizeof(funcStruct));

  // compute the derivatives for each (sub)function w.r.t. each parameter
  for (i = 0; i < sys->numParams; i++)
  { // initialize the data
    initialize_memLoc_derivAddr(sys, memLoc, derivAddr, totalOpCount);
    for (j = 0; j < totalF; j++)
    { // find the partial derivatives for the jth (sub)function w.r.t. ith parameter
      for (k = 0; k < totalF; k++)
        if (k < sys->numSubfuncs && sys->subFuncOrder[k] == j)
          diff_sys_param_subfunc(&derivs[k][i], sys, k, i, memLoc, derivAddr, totalOpCount);
        else if (k >= sys->numSubfuncs && sys->funcOrder[k - sys->numSubfuncs] == j)
          diff_sys_param_func(&derivs[k][i], sys, k - sys->numSubfuncs, i, memLoc, derivAddr, totalOpCount);
    }
  }

  // print the derivatives for each (sub)function w.r.t each parameter
  for (i = 0; i < totalF; i++)
  { // find the (sub)function corresponding to the ith one
    for (j = 0; j < totalF; j++)
      if (j < sys->numSubfuncs && sys->subFuncOrder[j] == i)
      { // jth subfunction
        for (k = 0; k < sys->numParams; k++)
        { // print derivs
          end += printOps(OUT, derivs[j][k].ops, derivs[j][k].num_ops);

          // clear derivs[j][k]
          free(derivs[j][k].ops);
          derivs[j][k].ops = NULL;
          derivs[j][k].num_ops = 0;
        }

        // clear derivs[j]
        free(derivs[j]);
      }
      else if (j >= sys->numSubfuncs && sys->funcOrder[j - sys->numSubfuncs] == i)
      { // jth function
        for (k = 0; k < sys->numParams; k++)
        { // print derivs
          end += printOps(OUT, derivs[j][k].ops, derivs[j][k].num_ops);

          // clear derivs[j][k]
          free(derivs[j][k].ops);
          derivs[j][k].ops = NULL;
          derivs[j][k].num_ops = 0;
        }

        // clear derivs[j]
        free(derivs[j]);
      }
  }

  // print an 'X' at the end of the array
  fprintf(OUT, "X\n");

  // print the rest of the instructions
  fprintf(OUT, "TotalRoomNeeded %d;\n", sys->firstFreeMemLoc);
  fprintf(OUT, "NVAR %d %d;\n", sys->numVars, sys->varsAddr);
  fprintf(OUT, "NPATHVAR %d %d;\n", sys->numPathVars, sys->pathVarsAddr);
  fprintf(OUT, "NPAR %d %d %d;\n", sys->numParams, sys->paramAddr, sys->funcAddr + sys->numFuncs);
  fprintf(OUT, "NFCN %d %d %d %d;\n", sys->numFuncs, sys->funcAddr, sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars, sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars + sys->numVars * sys->numFuncs);
  fprintf(OUT, "NCON %d %d;\n", sys->numConstants, sys->constAddr);
  fprintf(OUT, "NNUM %d %d;\n", sys->numNumbers, sys->numAddr);
  fprintf(OUT, "SUBFCN %d %d %d %d;\n", sys->numSubfuncs, sys->subFuncAddr, sys->subFuncAddr + sys->numSubfuncs, sys->subFuncAddr + sys->numSubfuncs + sys->numVars * sys->numSubfuncs);
  fprintf(OUT, "NUMINST %d;\n", end);
  fprintf(OUT, "CMPLX %d;\n", sys->constAddr); // first constant is always 'I'
  fprintf(OUT, "VARGPS %d;\n", sys->numVarGps);
  for (i = 0; i < sys->numVarGps; i++)
    fprintf(OUT, " %d", sys->varGpSizes[i]);
  fprintf(OUT, ";\n");
  fprintf(OUT, "RANDINDEX %d;\n", sys->randIndex);
  fprintf(OUT, "INSTCOUNT %d %d %d %d %d;\n", endUpdate, endParams, endFn, endPDerivs, endJvEval);
  // print the start and end for the subfunctions
  for (i = 0; i < sys->numSubfuncs; i++)
    fprintf(OUT, " %d %d", startSubfuncs[i], endSubfuncs[i]);
  // print the start and end for the functions
  for (i = 0; i < sys->numFuncs; i++)
    fprintf(OUT, " %d %d", startFuncs[i], endFuncs[i]);
  // print the start and end of the derivs for subfunctions
  for (i = 0; i < sys->numSubfuncs; i++)
    fprintf(OUT, " %d %d", startJvsub[i], endJvsub[i]);
  // print the start and end of the derivs for functions
  for (i = 0; i < sys->numFuncs; i++)
    fprintf(OUT, " %d %d", startJv[i], endJv[i]);
  fprintf(OUT, ";\n");
  // print the 'matrix' of the subfunctions 'below' each function
  for (i = 0; i < sys->numFuncs; i++)
  { // print the subfunctions below the ith function
    for (j = 0; j < sys->numSubfuncs; j++)
      fprintf(OUT, " %d", sys->subFuncOrder[j] < sys->funcOrder[i] ? 1 : 0);
    if (sys->numSubfuncs > 0)
      fprintf(OUT, ";\n"); 
  }
  fprintf(OUT, "X diff complete\n\n");

  // close OUT
  fclose(OUT);

  // clear memory
  free(startSubfuncs); free(endSubfuncs);
  free(startFuncs); free(endFuncs);
  free(startJvsub); free(endJvsub);
  free(startJv); free(endJv);
  free(memLoc);  
  free(derivAddr);

  return;
}

void initialize_memLoc_derivAddr(systemStruct *sys, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup memLoc                                           *
\***************************************************************/
{
  int i, j, count;

  // initialize derivAddr (set to all zeros)
  for (i = 0; i < totalOpCount; i++)
    derivAddr[i] = -2;

  // setup memLoc for the variables, path variables, constants, numbers, parameters, subfunction & function operations
  count = 0;
  // variables
  for (i = 0; i < sys->numVars; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->varsAddr + i;
    count++;
  }
  // path variables
  for (i = 0; i < sys->numPathVars; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->pathVarsAddr + i;
    count++;
  }
  // constants
  for (i = 0; i < sys->numConstants; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->constAddr + i;
    count++;
  }
  // numbers
  for (i = 0; i < sys->numNumbers; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->numAddr + i;
    count++;
  }
  // parameters
  for (i = 0; i < sys->numParams; i++)
    for (j = 0; j < sys->params[i].num_ops; j++)
    { // setup memLoc[count]
      memLoc[count] = sys->params[i].ops[j].memLoc;
      count++;
    }
  // subfunctions
  for (i = 0; i < sys->numSubfuncs; i++)
    for (j = 0; j < sys->subFuncs[i].num_ops; j++)
    { // setup memLoc[count]
      memLoc[count] = sys->subFuncs[i].ops[j].memLoc;
      count++;
    }
  // functions
  for (i = 0; i < sys->numFuncs; i++)
    for (j = 0; j < sys->funcs[i].num_ops; j++)
    { // setup memLoc[count]
      memLoc[count] = sys->funcs[i].ops[j].memLoc;
      count++;
    }

  return;
}

void diff_sys_pathVar(funcStruct *diff, systemStruct *sys, int currParam, int currPathVar, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate param w.r.t. pathVar - store to diff     *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, addr;

  // find where to store this partial derivative - after the functions 
  addr = sys->funcAddr + sys->numFuncs + currParam * sys->numPathVars + currPathVar;

  // setup for current path variable
  for (i = 0; i < sys->numPathVars; i++)
  { // derivAddr[numVars + i] == 'one' if i == currPathVar or 'zero' if i != currPathVar
    derivAddr[sys->numVars + i] = (i == currPathVar ? -1 : -2);
  }

  // now we need to setup the derivatives for the 'currParam' parameter w.r.t. 'currPathVar' path variable
  diff_forward_pathVars(sys, currParam, addr, memLoc, derivAddr, totalOpCount, &diff->num_ops, &diff->ops);

  return;
}

void diff_forward_pathVars(systemStruct *sys, int currParam, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currParam'                              *
\***************************************************************/
{
  int i, count;

  // differentiate the 'currParam' parameter 
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i <= currParam; i++)
    count += sys->params[i].num_ops;

  // differentiate the parameter
  diff_funcStruct(&sys->params[currParam], count, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_sys_var_subfunc(funcStruct *diff, systemStruct *sys, int currSubFunc, int currVar, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate subfunc w.r.t. var - store to diff       *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, addr;

  // find where to store this partial derivative - after the subfunctions
  addr = sys->subFuncAddr + sys->numSubfuncs + currSubFunc * sys->numVars + currVar;

  // seutp for current variable
  for (i = 0; i < sys->numVars; i++)
  { // derivAddr[i] = 'one' if i == currVarr or 'zero' if i != currVar
    derivAddr[i] = (i == currVar ? -1 : -2);
  }

  // now we need to setup the derivatives for the 'currSubFunc' subfunction w.r.t. 'currVar' variable
  diff_forward_vars_subfuncs(sys, currSubFunc, addr, memLoc, derivAddr, totalOpCount, &diff->num_ops, &diff->ops);

  return;
}

void diff_forward_vars_subfuncs(systemStruct *sys, int currSubFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currSubFunc'                            *
\***************************************************************/
{
  int i, count;

  // differentiate the 'currSubFunc' subfunction w.r.t. the 'currVar' variable
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
    count += sys->params[i].num_ops;
  for (i = 0; i <= currSubFunc; i++)
    count += sys->subFuncs[i].num_ops;

  // differentiate the subfunction
  diff_funcStruct(&sys->subFuncs[currSubFunc], count, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_sys_var_func(funcStruct *diff, systemStruct *sys, int currFunc, int currVar, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate func w.r.t. var - store to diff          *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, addr;

  // find where to store this partial derivative - after the functions & param derivs
  addr = sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars + currFunc * sys->numVars + currVar;

  // seutp for current variable
  for (i = 0; i < sys->numVars; i++)
  { // derivAddr[i] = 'one' if i == currVar, 'zero' if i != currVar
    derivAddr[i] = (i == currVar ? -1 : -2);
  }

  // now we need to setup the derivatives for the 'currFunc' function w.r.t. 'currVar' variable
  diff_forward_vars_funcs(sys, currFunc, addr, memLoc, derivAddr, totalOpCount, &diff->num_ops, &diff->ops);

  return;
}

void diff_forward_vars_funcs(systemStruct *sys, int currFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currFunc' w.r.t. 'currVar'              *
\***************************************************************/
{
  int i, count;

  // find where the last operation in the function is stored
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
    count += sys->params[i].num_ops;
  for (i = 0; i < sys->numSubfuncs; i++)
    count += sys->subFuncs[i].num_ops;
  for (i = 0; i <= currFunc; i++)
    count += sys->funcs[i].num_ops;

  // differentiate the function
  diff_funcStruct(&sys->funcs[currFunc], count, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_sys_param_subfunc(funcStruct *diff, systemStruct *sys, int currSubFunc, int currParam, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate subfunc w.r.t. param - store to diff     *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, addr, count;

  // find where to store this partial derivative - after the subfunctions & subfunctions w.r.t. vars
  addr = sys->subFuncAddr + sys->numSubfuncs + sys->numVars * sys->numSubfuncs + currSubFunc * sys->numParams + currParam;

  // seutp for current parameter
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
  { // derivAddr[i] = 'one' if i == currParam, 'zero' if i != currParam
    count += sys->params[i].num_ops;
    derivAddr[count - 1] = (i == currParam ? -1 : -2);
  }

  // now we need to setup the derivatives for the 'currSubFunc' subfunction w.r.t. 'currParam' parameter
  diff_forward_params_subfunc(sys, currSubFunc, addr, memLoc, derivAddr, totalOpCount, &diff->num_ops, &diff->ops);

  return;
}

void diff_forward_params_subfunc(systemStruct *sys, int currSubFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currSubFunc' w.r.t. 'currParam'         *
\***************************************************************/
{
  int i, count;

  // differentiate the 'currSubFunc' subfunction w.r.t. the 'currParam' parameter
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
    count += sys->params[i].num_ops;
  for (i = 0; i <= currSubFunc; i++)
    count += sys->subFuncs[i].num_ops;

  // differentiate the subfunction
  diff_funcStruct(&sys->subFuncs[currSubFunc], count, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_sys_param_func(funcStruct *diff, systemStruct *sys, int currFunc, int currParam, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate funcs w.r.t. params - store to diff      *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, addr, count = 0;

  // find where to store this partial derivative - after the functions, param derivs & derivs w.r.t. vars
  addr = sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars + sys->numVars * sys->numFuncs + currFunc * sys->numParams + currParam;

  // seutp for current parameter
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
  { // derivAddr[i] = 'one' if i == currParam, 'zero' if i != currParam
    count += sys->params[i].num_ops;
    derivAddr[count - 1] = (i == currParam ? -1 : -2);
  }

  // now we need to setup the derivatives for the 'currFunc' function w.r.t. 'currParam' parameter
  diff_forward_params_func(sys, currFunc, addr, memLoc, derivAddr, totalOpCount, &diff->num_ops, &diff->ops);

  return;
}

void diff_forward_params_func(systemStruct *sys, int currFunc, int storeAddr, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currFunc' w.r.t. 'currParam'            *
\***************************************************************/
{
  int i, count;

  // differentiate the 'currFunc' function w.r.t. the 'currParam' parameter
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)
    count += sys->params[i].num_ops;
  for (i = 0; i < sys->numSubfuncs; i++)
    count += sys->subFuncs[i].num_ops;
  for (i = 0; i <= currFunc; i++)
    count += sys->funcs[i].num_ops;

  // differentiate the function
  diff_funcStruct(&sys->funcs[currFunc], count, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_op(func_ops *ops, int *memLoc, int *derivAddr, int totalOpCount, int oneAddr, int *first_free_mem_loc, int *new_ops, func_ops **diff_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: generate the derivative of ops                         *
\***************************************************************/
{
  int i, indexIn0 = 0, indexIn1 = 0, indexLoc = 0;
  char op = ops->op;

  if (op == '=')
  { // derivative does not change - update the derivAddr for the new location with the old location
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break; 
      } 
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative
    derivAddr[indexLoc] = derivAddr[indexIn0];
  }
  else if (op == 'N')
  { // derivative needs to be negated
    // find the index for in[0]
    indexIn0 = -1; 
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break; 
      } 
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }
 
    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative
    if (derivAddr[indexIn0] == -2)
    { // derivative is still zero
      derivAddr[indexLoc] = -2;
    }
    else if (derivAddr[indexIn0] == -1)
    { // derivative is -1 - setup in new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], 'N', oneAddr, -1);
    }
    else
    { // negate the old location - setup in new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], 'N', derivAddr[indexIn0], -1);
    }
  }
  else if (op == '+')
  { // derivative is sum of the derivatives
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }
    // find the index for in[1]
    indexIn1 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[1])
      { // setup indexIn1
        indexIn1 = i;
        break;
      }
    // error check
    if (indexIn1 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }
    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative
    if (derivAddr[indexIn0] == -2)
    { // derivative is just the derivative of in[1]
      derivAddr[indexLoc] = derivAddr[indexIn1];
    }
    else if (derivAddr[indexIn0] == -1)
    { // derivative is 1 + diff(in[1])
      if (derivAddr[indexIn1] == -2)
      { // derivative is just 'one'
        derivAddr[indexLoc] = -1;
      }
      else if (derivAddr[indexIn1] == -1)
      { // derivative is 1 + 1 - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '+', oneAddr, oneAddr);
      }
      else
      { // derivative is 1 + diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '+', oneAddr, derivAddr[indexIn1]);
      }
    }
    else
    { // derivative is diff(in[0]) + diff(in[1])
      if (derivAddr[indexIn1] == -2)
      { // derivative is just diff(in[0])
        derivAddr[indexLoc] = derivAddr[indexIn0];
      }
      else if (derivAddr[indexIn1] == -1)
      { // derivative is diff(in[0]) + 1 - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '+', derivAddr[indexIn0], oneAddr);
      } 
      else
      { // derivative is diff(in[0]) + diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '+', derivAddr[indexIn0], derivAddr[indexIn1]);
      }
    }
  }
  else if (op == '-')
  { // derivative is difference of the derivatives
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // find the index for in[1]
    indexIn1 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[1])
      { // setup indexIn1
        indexIn1 = i;
        break;
      }
    // error check
    if (indexIn1 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }
    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative
    if (derivAddr[indexIn1] == -2)
    { // derivative is just the derivative of in[0]
      derivAddr[indexLoc] = derivAddr[indexIn0];
    }
    else if (derivAddr[indexIn1] == -1)
    { // derivative is diff(in[0]) - 1
      if (derivAddr[indexIn0] == -2)
      { // derivative is -1 - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], 'N', oneAddr, -1);
      }
      else if (derivAddr[indexIn0] == -1)
      { // derivative is 1 - 1 = 0 = 'zero'
        derivAddr[indexLoc] = -2;
      }
      else
      { // derivative is diff(in[0]) - 1 - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '-', derivAddr[indexIn0], oneAddr);
      }
    }
    else
    { // derivative is diff(in[0]) - diff(in[1])
      if (derivAddr[indexIn0] == -2)
      { // derivative is just -diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], 'N', derivAddr[indexIn1], -1);
      }
      else if (derivAddr[indexIn0] == -1)
      { // derivative is 1 - diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '-', oneAddr, derivAddr[indexIn1]);
      }
      else
      { // derivative is diff(in[0]) - diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '-', derivAddr[indexIn0], derivAddr[indexIn1]);
      }
    }
  }
  else if (op == '/')
  { // we must have the denominator being constant and so we only need to worry about in[0]
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative
    if (derivAddr[indexIn0] == -2)
    { // derivative is still 0
      derivAddr[indexLoc] = -2;
    }
    else if (derivAddr[indexIn0] == -1)
    { // derivative is 1 / in[1] - setup to new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '/', oneAddr, ops->in[1]);
    }
    else
    { // derivative is diff(in[0]) / in[1] - setup to new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '/', derivAddr[indexIn0], ops->in[1]);
    }
  }
  else if (op == '*')
  { // use multiplication rule
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // find the index for in[1]
    indexIn1 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[1])
      { // setup indexIn1
        indexIn1 = i;
        break;
      }
    // error check
    if (indexIn1 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }
    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative = diff(in[0]) * in[1] + diff(in[1]) * in[0]
    if (derivAddr[indexIn0] == -2)
    { // derivative is in[0] * diff(in[1])
      if (derivAddr[indexIn1] == -2)
      { // derivative is 'zero'
        derivAddr[indexLoc] = -2;
      }
      else if (derivAddr[indexIn1] == -1)
      { // derivative is in[0]
        derivAddr[indexLoc] = ops->in[0];
      }
      else
      { // derivative is in[0] * diff(in[1]) - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '*', ops->in[0], derivAddr[indexIn1]);
      }
    }
    else if (derivAddr[indexIn0] == -1)
    { // derivative is in[1] + diff(in[1]) * in[0]
      if (derivAddr[indexIn1] == -2)
      { // derivative is just in[1]
        derivAddr[indexLoc] = ops->in[1];
      }
      else if (derivAddr[indexIn1] == -1)
      { // derivative is in[1] + in[0] - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '+', ops->in[0], ops->in[1]);
      }
      else
      { // derivative is in[1] + diff(in[1]) * in[0] - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 2;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[0], derivAddr[indexIn1]);
        copyToOp(&(*diff_ops)[1], derivAddr[indexLoc], '+', ops->in[1], *first_free_mem_loc);
        (*first_free_mem_loc)++;
      }
    }
    else
    { // derivative is diff(in[0]) * in[1] + diff(in[1]) * in[0]
      if (derivAddr[indexIn1] == -2)
      { // derivative is just diff(in[0]) * in[1] - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc], '*', ops->in[1], derivAddr[indexIn0]);
      }
      else if (derivAddr[indexIn1] == -1)
      { // derivative is diff(in[0]) * in[1] + in[0] - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 2;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[1], derivAddr[indexIn0]);
        copyToOp(&(*diff_ops)[1], derivAddr[indexLoc], '+', ops->in[0], *first_free_mem_loc);
        (*first_free_mem_loc)++;
      }
      else
      { // derivative is diff(in[0]) * in[1] + diff(in[1]) * in[0] - setup in new location
        derivAddr[indexLoc] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 3;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[1], derivAddr[indexIn0]);
        (*first_free_mem_loc)++;

        copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '*', ops->in[0], derivAddr[indexIn1]);
        copyToOp(&(*diff_ops)[2], derivAddr[indexLoc], '+', *first_free_mem_loc - 1, *first_free_mem_loc);
        (*first_free_mem_loc)++;
      }
    }
  }
  else if (op == '^')
  { // we must have the exponent being constant and so we only need to worry about in[0]
    // find the index for in[0]
    indexIn0 = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->in[0])
      { // setup indexIn0
        indexIn0 = i;
        break;
      }
    // error check
    if (indexIn0 == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // find the index for memLoc
    indexLoc = -1;
    for (i = 0; i < totalOpCount; i++)
      if (memLoc[i] == ops->memLoc)
      { // setup indexLoc
        indexLoc = i;
        break;
      }
    // error check
    if (indexLoc == -1)
    {
      printf("ERROR: The memory location can not be found!\n");
      bexit(ERROR_LOOP_FAILURE);
    }

    // set the derivative = in[1] * (in[0] ^ (in[1] - 1)) * diff(in[0])
    if (derivAddr[indexIn0] == -2)
    { // derivative is still 0
      derivAddr[indexLoc] = -2;
    }
    else if (derivAddr[indexIn0] == -1)
    { // derivative is in[1] * in[0] ^ (in[1] - 1) - setup in new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 3;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));

      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '-', ops->in[1], oneAddr);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '^', ops->in[0], *first_free_mem_loc - 1);
      copyToOp(&(*diff_ops)[2], derivAddr[indexLoc], '*', *first_free_mem_loc, ops->in[1]);
      (*first_free_mem_loc)++;
    }
    else
    { // derivative is in[1] * (in[0] ^ (in[1] - 1)) * diff(in[0]) - setup in new location
      derivAddr[indexLoc] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 4;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));

      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '-', ops->in[1], oneAddr);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '^', ops->in[0], *first_free_mem_loc - 1);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[2], *first_free_mem_loc, '*', *first_free_mem_loc - 1, ops->in[1]);
      copyToOp(&(*diff_ops)[3], derivAddr[indexLoc], '*', *first_free_mem_loc, derivAddr[indexIn0]);
      (*first_free_mem_loc)++;
    }
  } 
  else
  {
    printf("We have a wrong operation ('%c')!!\n", op);
    bexit(ERROR_CONFIGURATION);
  }

  return;
}

int copyInstructions(FILE *OUT, FILE *IN, char endCh)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: return the number of operations printed        *
* NOTES: read in the operation from IN and print to OUT         *
\***************************************************************/
{
  int loc, in[2], numOps = 0;
  char op;

  do
  { // read in the next operation
    op = fgetc(IN);

    if (isUnary(op))
    { // read in and print unary operation
      numOps += 3;
      fscanf(IN, "%d %d;\n", &loc, &in[0]);
      fprintf(OUT, "%d %d %d ", op, loc, in[0]);
    }
    else if (op == '+' || op == '-' || op == '*' || op == '/' || op == '^')
    { // read in the binary operation
      numOps += 4;
      fscanf(IN, "%d %d %d;\n", &loc, &in[0], &in[1]);
      fprintf(OUT, "%d %d %d %d ", op, loc, in[0], in[1]);
    }
    else if (op != 'X')
    {
      printf("ERROR: There is an invalid opertion ('%c')!\n", op);
      bexit(ERROR_CONFIGURATION);
    }
  } while (op != endCh);

  return numOps;
}

void setupSystemStructure(systemStruct *sys, char *fileName, int putIntoOrder)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: read in the system from fileName & setup the structure *
\***************************************************************/
{ // open up the file that contains the system
  FILE *FUNC = fopen(fileName, "r");
  if (FUNC == NULL)
  {
    printf("ERROR: '%s' does not exist!\n", fileName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  int i, numNames, roomNeeded, numVars, varAddr, numPathvars, pathvarAddr, numParams, paramAddr, paramDerivAddr;
  int numFuncs, funcAddr, funcDerivWRTVarsAddr, funcDerivWRTParamsAddr;
  int numConsts, constAddr, numNums, numAddr, numSubfuncs, subfuncAddr, subfuncDerivWRTVarsAddr, subfuncDerivWRTParamsAddr;
  int IAddr, zeroAddr, oneAddr, num_var_gps, rand_index;
  int tempInt1, tempInt2;

  // read in the important information
  fscanf(FUNC, "RoomBeforeTemps %d;\n", &numNames);
  fscanf(FUNC, "TotalRoomNeeded %d;\n", &roomNeeded);

  // store roomNeeded
  sys->firstFreeMemLoc = roomNeeded + 1;

  fscanf(FUNC, "NVAR %d %d;\n", &numVars, &varAddr);
  fscanf(FUNC, "NPATHVAR %d %d;\n", &numPathvars, &pathvarAddr);
  fscanf(FUNC, "NPAR %d %d %d;\n", &numParams, &paramAddr, &paramDerivAddr);
  fscanf(FUNC, "NFCN %d %d %d %d;\n", &numFuncs, &funcAddr, &funcDerivWRTVarsAddr, &funcDerivWRTParamsAddr);
  fscanf(FUNC, "NCON %d %d;\n", &numConsts, &constAddr);
  fscanf(FUNC, "NNUM %d %d;\n", &numNums, &numAddr);
  fscanf(FUNC, "SUBFCN %d %d %d %d;\n", &numSubfuncs, &subfuncAddr, &subfuncDerivWRTVarsAddr, &subfuncDerivWRTParamsAddr);

  // setup variables
  sys->numVars = numVars;
  sys->varsAddr = varAddr;
  
  // setup path variables
  sys->numPathVars = numPathvars;
  sys->pathVarsAddr = pathvarAddr;

  // allocate for parameters
  sys->numParams = numParams;
  sys->paramAddr = paramAddr;
  sys->params = (funcStruct *)bmalloc(sys->numParams * sizeof(funcStruct));
  for (i = 0; i < numParams; i++)
  {
    sys->params[i].num_ops = 1;
    sys->params[i].ops = (func_ops *)bmalloc(sys->params[i].num_ops * sizeof(func_ops));
  }

  // allocate for subfunctions
  sys->numSubfuncs = numSubfuncs;
  sys->subFuncAddr = subfuncAddr;
  sys->subFuncs = (funcStruct *)bmalloc(sys->numSubfuncs * sizeof(funcStruct));
  for (i = 0; i < numSubfuncs; i++)
  {
    sys->subFuncs[i].num_ops = 1;
    sys->subFuncs[i].ops = (func_ops *)bmalloc(sys->subFuncs[i].num_ops * sizeof(func_ops));
  }

  // allocate functions
  sys->numFuncs = numFuncs;
  sys->funcAddr = funcAddr;
  sys->funcs = (funcStruct *)bmalloc(sys->numFuncs * sizeof(funcStruct));
  for (i = 0; i < numFuncs; i++)
  {
    sys->funcs[i].num_ops = 1;
    sys->funcs[i].ops = (func_ops *)bmalloc(sys->funcs[i].num_ops * sizeof(func_ops));
  }

  fscanf(FUNC, "CMPLX %d %d %d;\n", &IAddr, &tempInt1, &tempInt2);
  fscanf(FUNC, "ZERO %d;\n", &zeroAddr);
  fscanf(FUNC, "ONE %d;\n", &oneAddr);
  fscanf(FUNC, "VARGPS %d;\n", &num_var_gps);

  // setup varaible groups in var_gp_sizes
  sys->numVarGps = num_var_gps;
  sys->varGpSizes = (int *)bmalloc(num_var_gps * sizeof(int));
  for (i = 0; i < num_var_gps; i++)
  {
    fscanf(FUNC, " %d", &tempInt1);
    sys->varGpSizes[i] = tempInt1;
  }

  fscanf(FUNC, ";\n");
  fscanf(FUNC, "RANDINDEX %d;\n", &rand_index);
  sys->randIndex = rand_index;

  // allocte update
  sys->numUpdate = 1;
  sys->updateOps = (func_ops *)bmalloc(sys->numUpdate * sizeof(func_ops));

  // setup other information
  sys->numConstants = numConsts;
  sys->constAddr = constAddr;
  sys->numNumbers = numNums;
  sys->numAddr = numAddr;

  // setup the structures
  readInOps(sys, FUNC, putIntoOrder);

  // close FUNC
  fclose(FUNC);

  return;
}

void readInOps(systemStruct *sys, FILE *FUNC, int putIntoOrder)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the operations in sys from FUNC                  *
\***************************************************************/
{
  int numParamOps = 1, numFuncOps = 1;
  func_ops *paramOps = (func_ops *)bmalloc(numParamOps * sizeof(func_ops));
  func_ops *funcOps = (func_ops *)bmalloc(numFuncOps * sizeof(func_ops));

  // start with the update steps
  fscanf(FUNC, "BEGIN UPDATE;\n");
  setupOps(&sys->updateOps, &sys->numUpdate, FUNC, 'B'); // read until B in BEGIN PARAM

  // read in parameter steps
  fscanf(FUNC, "EGIN PARAM;\n");
  setupOps(&paramOps, &numParamOps, FUNC, 'B'); // read until B in BEGIN FUNCTION

  // setup the parameters from paramOps
  setupParams(sys, paramOps, numParamOps);

  // read in function steps
  fscanf(FUNC, "EGIN FUNCTION;\n");
  setupOps(&funcOps, &numFuncOps, FUNC, 'E'); // read until E in END

  // setup the subfunctions & functions from funcOps
  setupFuncs_subFuncs(sys, funcOps, numFuncOps, putIntoOrder);

  // clear funcOps
  free(funcOps);

  return;
}

void setupParams(systemStruct *sys, func_ops *paramOps, int numParamOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decompose paramOps in the parameters                   *
\***************************************************************/
{
  int i, j, count, boolCont;

  // initialize j
  j = 0;

  // find the sizes of the parameters and allocate them
  for (i = 0; i < sys->numParams; i++)
  { // initialize
    count = 0;
    boolCont = 1;

    // determine if we have more operations to look at
    if (j + count >= numParamOps)
    {
      printf("ERROR: It appears that a user-defined parameter was not setup properly!\n");
      bexit(ERROR_INPUT_SYSTEM);
    } 

    do
    { // determine if this operation calculated the ith parameter
      if (paramOps[j + count].memLoc == sys->paramAddr + i)
      { // this is the last operation
        boolCont = 0;
      }
      // increment count
      count++;
    } while (boolCont);
 
    // store the size and allocate the ops
    sys->params[i].num_ops = count;
    sys->params[i].ops = (func_ops *)bmalloc(count * sizeof(func_ops));
 
    // store the ops
    for (count = 0; count < sys->params[i].num_ops; count++)
    {
      sys->params[i].ops[count].memLoc = paramOps[j + count].memLoc;
      sys->params[i].ops[count].op = paramOps[j + count].op;
      sys->params[i].ops[count].in[0] = paramOps[j + count].in[0];
      sys->params[i].ops[count].in[1] = paramOps[j + count].in[1];
 
      // initialize other data
      sys->params[i].ops[count].lastUsed = -1;
      sys->params[i].ops[count].numDiffInst = 0;
    }
 
    // update j
    j += count;
  }

  return;
}

void setupFuncs_subFuncs(systemStruct *sys, func_ops *funcOps, int numFuncOps, int putIntoOrder)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: decompose funcOps in the functions and subfunctions    *
\***************************************************************/
{
  int i, j, count, boolCont, isSubFunc, isFunc, subFuncNum, funcNum, total = sys->numSubfuncs + sys->numFuncs;

  // allocate orderings
  sys->subFuncOrder = (int *)bmalloc(sys->numSubfuncs * sizeof(int));
  sys->funcOrder = (int *)bmalloc(sys->numFuncs * sizeof(int));

  // initialize j
  j = 0;
  // loop over the functions and subfunctions to find them inside of funcOps
  for (i = 0; i < total; i++)
  { // find the end of the next function/subfunction
    count = isSubFunc = isFunc = 0;
    boolCont = 1;
    subFuncNum = funcNum = -1;

    do
    { // determine if this operation calculated the end of a subfunction
      if (sys->subFuncAddr <= funcOps[j + count].memLoc && funcOps[j + count].memLoc < sys->subFuncAddr + sys->numSubfuncs) 
      { // this is the last operation of a subfunction
        boolCont = 0;
        isSubFunc = 1;
        subFuncNum =  funcOps[j + count].memLoc - sys->subFuncAddr;
      }

      if (boolCont && sys->funcAddr <= funcOps[j + count].memLoc && funcOps[j + count].memLoc < sys->funcAddr + sys->numFuncs)
      { // this is the last operation of jth function
        boolCont = 0;
        isFunc = 1;
        funcNum = funcOps[j + count].memLoc - sys->funcAddr;
      }
      // increment count
      count++;
    } while (boolCont);

    // store the size and allocate the ops
    if (isSubFunc)
    { // store the size and allocate the ops for this subfunction
      sys->subFuncOrder[subFuncNum] = i;
      sys->subFuncs[subFuncNum].num_ops = count;
      sys->subFuncs[subFuncNum].ops = (func_ops *)bmalloc(count * sizeof(func_ops));

      // store the ops
      for (count = 0; count < sys->subFuncs[subFuncNum].num_ops; count++)
      {
        sys->subFuncs[subFuncNum].ops[count].memLoc = funcOps[j + count].memLoc;
        sys->subFuncs[subFuncNum].ops[count].op = funcOps[j + count].op;
        sys->subFuncs[subFuncNum].ops[count].in[0] = funcOps[j + count].in[0];
        sys->subFuncs[subFuncNum].ops[count].in[1] = funcOps[j + count].in[1];

        // initialize other data
        sys->subFuncs[subFuncNum].ops[count].lastUsed = -1;
        sys->subFuncs[subFuncNum].ops[count].numDiffInst = 0;
      }
    }
    else if (isFunc)
    { // store the size and allocate the ops for this function
      sys->funcOrder[funcNum] = i;
      sys->funcs[funcNum].num_ops = count;
      sys->funcs[funcNum].ops = (func_ops *)bmalloc(count * sizeof(func_ops));

      // store the ops
      for (count = 0; count < sys->funcs[funcNum].num_ops; count++)
      {
        sys->funcs[funcNum].ops[count].memLoc = funcOps[j + count].memLoc;
        sys->funcs[funcNum].ops[count].op = funcOps[j + count].op;
        sys->funcs[funcNum].ops[count].in[0] = funcOps[j + count].in[0];
        sys->funcs[funcNum].ops[count].in[1] = funcOps[j + count].in[1];

        // initialize other data
        sys->funcs[funcNum].ops[count].lastUsed = -1;
        sys->funcs[funcNum].ops[count].numDiffInst = 0;
      }
    }
    // update j
    j += count;
  }

  if (putIntoOrder)
  { // put the subfunctions & functions into order - subfunctions first and then functions
    for (i = 0; i < total; i++)
      if (i < sys->numSubfuncs)
      { // set to i
        sys->subFuncOrder[i] = i;
      }
      else
      { // set to i
        sys->funcOrder[i - sys->numSubfuncs] = i;
      }
  }

  return;
}

void setupOps(func_ops **ops, int *numOps, FILE *FUNC, char stopCh)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the operations in ops from FUNC                  *
\***************************************************************/
{
  int count = 0, tempInt1, tempInt2, tempInt3;
  char ch;

  do
  { // read in the next operation
    ch = fgetc(FUNC);

    if (ch != stopCh)
    { // read in the operation
      if (isUnary(ch))
      { // unary operation
        fscanf(FUNC, "%d %d;\n", &tempInt1, &tempInt2);

        (*ops)[count].memLoc = tempInt1;
        (*ops)[count].op = ch;
        (*ops)[count].in[0] = tempInt2;
        (*ops)[count].in[1] = -1; // not used

        (*ops)[count].lastUsed = -1; // intialize
      }
      else if (ch == '+' || ch == '-' || ch == '*' || ch == '/' || ch == '^')
      { // binary operation
        fscanf(FUNC, "%d %d %d;\n", &tempInt1, &tempInt2, &tempInt3);

        (*ops)[count].memLoc = tempInt1;
        (*ops)[count].op = ch;
        (*ops)[count].in[0] = tempInt2;
        (*ops)[count].in[1] = tempInt3;

        (*ops)[count].lastUsed = -1; // intialize
      }
      else
      {
        printf("ERROR: There is an invalid opertion ('%c')!\n", ch);
        bexit(ERROR_CONFIGURATION);
      }

      // increment count
      count++;

      // see if we need to increase the size of the memory
      if (count >= *numOps)
      { // double the size
        *numOps = *numOps * 2;
        *ops = (func_ops *)brealloc(*ops, *numOps * sizeof(func_ops));
      }
    }

  } while (ch != stopCh);

  // make sure that we have set the number correctly
  *numOps = count;
  *ops = (func_ops *)brealloc(*ops, *numOps * sizeof(func_ops));

  return;
}

int printOps(FILE *OUT, func_ops *funcOps, int numOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of instructions printed                 *
* NOTES: prints the operations to OUT                           *
\***************************************************************/
{
  int i, count = 0;

  for (i = 0; i < numOps; i++)
    if (isUnary(funcOps[i].op))
    {
      count += 3;
      fprintf(OUT, "%d %d %d ", funcOps[i].op, funcOps[i].memLoc, funcOps[i].in[0]);
    }
    else if (funcOps[i].op == '+' || funcOps[i].op == '-' || funcOps[i].op == '*' || funcOps[i].op == '/' || funcOps[i].op == '^')
    {
      count += 4;
      fprintf(OUT, "%d %d %d %d ", funcOps[i].op, funcOps[i].memLoc, funcOps[i].in[0], funcOps[i].in[1]);
    }
    else
    {
      printf("ERROR: '%c' is an invalid operation!\n", funcOps[i].op);
      bexit(ERROR_CONFIGURATION);
    }

  return count;
}

void clearFuncStruct(funcStruct *funcs, int numFuncs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear funcs                                            *
\***************************************************************/
{
  int i;

  for (i = numFuncs - 1; i >= 0; i--)
  { // clear ops
    free(funcs[i].ops);
  }
  free(funcs);

  return;
}

void clearVarOps(var_ops *vars, int numVars)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear funcs                                            *
\***************************************************************/
{
  int i;

  for (i = numVars - 1; i >= 0; i--)
  { // clear ops
    free(vars[i].derivAddr);
    free(vars[i].varUsed);
  }

  return;
}

void clearSystemStructure(systemStruct *sys)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear sys                                              *
\***************************************************************/
{
  // clear varGpSizes
  free(sys->varGpSizes);

  // clear updateOps
  free(sys->updateOps);

  // clear params
  clearFuncStruct(sys->params, sys->numParams);

  // clear subFuncs
  clearFuncStruct(sys->subFuncs, sys->numSubfuncs);

  // clear funcs
  clearFuncStruct(sys->funcs, sys->numFuncs);

  // clear orders
  free(sys->subFuncOrder);
  free(sys->funcOrder);

  return;
}

void checkExp(systemStruct *sys, num_t *nums, int allowNegExp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: check that the exponents are what they should be!      *
\***************************************************************/
{
  int i, j, retVal, totalF = sys->numSubfuncs + sys->numFuncs, memSize = sys->firstFreeMemLoc, oneAddr = sys->numAddr + 1; // '1' is the second number
  comp_d *mem = (comp_d *)bmalloc(memSize * sizeof(comp_d));

  // set mem to zero
  for (i = 0; i < memSize; i++)
    set_zero_d(mem[i]);

  // setup the numbers
  for (i = 0; i < sys->numNumbers; i++)
    mem[i + sys->numAddr]->r = mpq_get_d(nums[i].rat);

  // setup I
  set_double_d(mem[sys->constAddr], 0, 1); // first constant is always I
  // setup Pi
  set_double_d(mem[sys->constAddr + 1], M_PI, 0); // second constant is always Pi

  // setup variables (random)
  for (i = 0; i < sys->numVars; i++)
  {
    get_comp_rand_d(mem[i + sys->varsAddr]);
  }

  // setup path variables (random)
  for (i = 0; i < sys->numPathVars; i++)
  {
    get_comp_rand_d(mem[i + sys->pathVarsAddr]);
  }

  // start going through the evaluation instructions looking for bad exponents!

  // update instructions
  retVal = checkExp_ops(sys->updateOps, sys->numUpdate, mem, oneAddr, allowNegExp);

  if (retVal == 1)
  { // print error message
    printf("ERROR: When seting up the functions, there was a bad exponent!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  else if (retVal == 2)
  { // print error message
    printf("ERROR: When setting up the functions, an exponent was found that was not an integer!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  else if (retVal == 3)
  { // print error message
    printf("ERROR: When setting up the functions, an exponent was found that was negative!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // parameter instructions
  for (i = 0; i < sys->numParams; i++)
  { // check ith parameter
    retVal = checkExp_ops(sys->params[i].ops, sys->params[i].num_ops, mem, oneAddr, allowNegExp);

    if (retVal == 1)
    { // print error message
      printf("ERROR: When seting up the functions, there was a bad exponent!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
    else if (retVal == 2)
    { // print error message
      printf("ERROR: When setting up the functions, an exponent was found that was not an integer!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
    else if (retVal == 3)
    { // print error message
      printf("ERROR: When setting up the functions, an exponent was found that was negative!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
  }

  // subfunctions & functions - evaluate in order they were defined
  for (i = 0; i < totalF; i++)
  { // find what is the ith subfunc/function
    for (j = 0; j < totalF; j++)
      if (j < sys->numSubfuncs && sys->subFuncOrder[j] == i)
      { // check the jth subfunction
        retVal = checkExp_ops(sys->subFuncs[j].ops, sys->subFuncs[j].num_ops, mem, oneAddr, allowNegExp);

        if (retVal == 1)
        { // print error message
          printf("ERROR: When seting up the functions, there was a bad exponent!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
        else if (retVal == 2)
        { // print error message
          printf("ERROR: When setting up the functions, an exponent was found that was not an integer!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
        else if (retVal == 3)
        { // print error message
          printf("ERROR: When setting up the functions, an exponent was found that was negative!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
      }
      else if (j >= sys->numSubfuncs && sys->funcOrder[j - sys->numSubfuncs] == i)
      { // check the (j - numSubfuncs)th function
        retVal = checkExp_ops(sys->funcs[j - sys->numSubfuncs].ops, sys->funcs[j - sys->numSubfuncs].num_ops, mem, oneAddr, allowNegExp);

        if (retVal == 1)
        { // print error message
          printf("ERROR: When seting up the functions, there was a bad exponent!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
        else if (retVal == 2)
        { // print error message
          printf("ERROR: When setting up the functions, an exponent was found that was not an integer!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
        else if (retVal == 3)
        { // print error message
          printf("ERROR: When setting up the functions, an exponent was found that was negative!\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
      }
  }

  // clear mem
  free(mem);

  return;
}

int checkExp_ops(func_ops *ops, int num_ops, comp_d *mem, int oneAddr, int allowNegExp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: check that the exponents are what they should be - ops *
\***************************************************************/
{
  int i, retVal = 0;

  // loop through the operations
  for (i = 0; i < num_ops && !retVal; i++)
  { // see if the operation is '^'
    if (ops[i].op == '^')
    { // see if this exponentiation operation is okay
      retVal = checkExp_op(&ops[i], mem, oneAddr, allowNegExp);
    }

    if (!retVal)
    { // perform the operation
      performOp(&ops[i], mem);
    }
  }

  return retVal;
}

int checkExp_op(func_ops *op, comp_d *mem, int oneAddr, int allowNegExp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: update the exponentiation operation if it can be fixed!*
\***************************************************************/
{ // make sure that the operation is really an exponentiation
  if (op->op != '^')
    return 0;

  // check this exp. operation
  int retVal = 0;

  // make sure that exp->i is '0'
  verifyRealExponent(mem[op->in[1]]);

  // verify exponent is integer if base is not real
  verifyBaseExponent(mem[op->in[0]], mem[op->in[1]]);

  // if okay, see if exp is 0 or 1 - correct it if it is!
  if (!retVal && mem[op->in[1]]->r == 0)
  { // set the location to 1 since the exponent is zero
    op->op = '='; // leave memLoc alone
    op->in[0] = oneAddr;
    op->in[1] = -1; // unary operation
  }
  else if (!retVal && mem[op->in[1]]->r == 1)
  { // set the location to the base since the exponent is one
    op->op = '='; // leave memLoc & in[0] alone
    op->in[1] = -1; // unary operation
  }

  return retVal;
}

void performOp(func_ops *op, comp_d *mem)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform the operation                                  *
\***************************************************************/
{
  if (op->op == '=')
  {
    set_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else if (op->op == '*')
  {
    mul_d(mem[op->memLoc], mem[op->in[0]], mem[op->in[1]]);
  }
  else if (op->op == '+')
  {
    add_d(mem[op->memLoc], mem[op->in[0]], mem[op->in[1]]);
  }
  else if (op->op == '-')
  {
    sub_d(mem[op->memLoc], mem[op->in[0]], mem[op->in[1]]);
  }
  else if (op->op == '^')
  { // verify real exponent
    verifyRealExponent(mem[op->in[1]]);

    // verify exponent is integer if base is not real
    verifyBaseExponent(mem[op->in[0]], mem[op->in[1]]);

    // perform the operation
    exp_d(mem[op->memLoc], mem[op->in[0]], mem[op->in[1]]->r);
  }
  else if (op->op == 'N')
  {
    neg_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else if (op->op == '/')
  {
    div_d(mem[op->memLoc], mem[op->in[0]], mem[op->in[1]]);
  }
  else if (op->op == 'C')
  {
    cos_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else if (op->op == 'S')
  {
    sin_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else if (op->op == 'X')
  {
    exponential_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else
  { // immediately exit if something is wrong with the program
    printf("\nERROR: Not a valid operation (%d)!\n", op->op);
    bexit(ERROR_CONFIGURATION);
  }

  return;
}

void copyToOp(func_ops *fop, int memLoc, char op, int in0, int in1)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy data to fop                                       *
\***************************************************************/
{
  fop->memLoc = memLoc;
  fop->op = op;
  fop->in[0] = in0;
  fop->in[1] = in1;
  
  return;
}

void diff_funcStruct(funcStruct *func, int endFuncCount, int storeAddr, int zeroAddr, int oneAddr, int *firstFreeMemLoc, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate func                                     *
\***************************************************************/
{
  int i, j, curr_loc, *op_count = (int *)bmalloc(func->num_ops * sizeof(int));
  func_ops **temp_ops = (func_ops **)bmalloc(func->num_ops * sizeof(func_ops *));

  // initialize derivCount
  *derivCount = 0;

  // differentiate each operation for this function w.r.t. this variable
  for (i = 0; i < func->num_ops; i++)
  { // setup op_count[i] & temp_ops[i]
    op_count[i] = 0;
    temp_ops[i] = NULL;

    diff_op(&func->ops[i], memLoc, derivAddr, totalOpCount, oneAddr, firstFreeMemLoc, &op_count[i], &temp_ops[i]);

    // update the total derivative count
    *derivCount += op_count[i];
  }

  // increment the derivative count by 1 - storing to the correct location
  *derivCount += 1;

  // setup deriv_ops
  curr_loc = 0;
  *deriv_ops = (func_ops *)bmalloc(*derivCount * sizeof(func_ops));
  for (i = 0; i < func->num_ops; i++)
  {
    for (j = 0; j < op_count[i]; j++)
    { // copy to deriv_ops
      copyToOp(&(*deriv_ops)[curr_loc], temp_ops[i][j].memLoc, temp_ops[i][j].op, temp_ops[i][j].in[0], temp_ops[i][j].in[1]);

      // update curr_loc
      curr_loc++;
    }
    // clear temp_ops[i]
    free(temp_ops[i]);
  }

  // store the partial derivative to the correct location
  if (derivAddr[endFuncCount - 1] == -2)
  { // derivative is zero
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', zeroAddr, -1);
  }
  else if (derivAddr[endFuncCount - 1] == -1)
  { // derivative is one
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', oneAddr, -1);
  }
  else
  { // derivative is at derivAddr
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', derivAddr[endFuncCount - 1], -1);
  }

  // clear temp_ops & op_count
  free(temp_ops);
  free(op_count);

  return;
}

int diff_funcStruct_test(funcStruct *func, int endFuncCount, int storeAddr, int zeroAddr, int oneAddr, int firstFreeTemp, int *memLoc, int *derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -2, -1 or something else                       *
* NOTES: determine if deriv is 0, 1, or something else          *
\***************************************************************/
{
  int i, *op_count = (int *)bmalloc(func->num_ops * sizeof(int)), retVal = 0, firstFreeMemLoc = firstFreeTemp;
  func_ops **temp_ops = (func_ops **)bmalloc(func->num_ops * sizeof(func_ops *));

  // differentiate each operation for this function w.r.t. this variable
  for (i = 0; i < func->num_ops; i++)
  { // setup op_count[i] & temp_ops[i]
    op_count[i] = 0;
    temp_ops[i] = NULL;

    diff_op(&func->ops[i], memLoc, derivAddr, totalOpCount, oneAddr, &firstFreeMemLoc, &op_count[i], &temp_ops[i]);
  }

  // determine what the derivative is
  retVal = derivAddr[endFuncCount - 1];
  if (retVal != -2 && retVal != -1)
    retVal = storeAddr;

  // clear temp_ops & op_count
  for (i = 0; i < func->num_ops; i++)
    free(temp_ops[i]);
  free(temp_ops);
  free(op_count);

  return retVal;
}

void setup_derivAddr_subfunc_vals(systemStruct *sys, int *derivAddr, int currVar, int **subFuncVals)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup derivAddr for the current 'variable'             *
\***************************************************************/
{
  int i, count = 0;

  // move past variables, path variables, constants, numbers, and parameters
  count = sys->numVars + sys->numPathVars + sys->numConstants + sys->numNumbers;
  for (i = 0; i < sys->numParams; i++)  
    count += sys->params[i].num_ops;

  // setup subfunctions
  for (i = 0; i < sys->numSubfuncs; i++)
  { // move past ops
    count += sys->subFuncs[i].num_ops;
    // setup derivAddr
    derivAddr[count-1] = subFuncVals[i][currVar];
  }

  return;
}

// NEW FUNCTIONS

void diff_file(char *fileName, int putIntoOrder, int allowNegExp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: turns file into "arr.out"                              *
\***************************************************************/
{
  systemStruct sys;
  num_t *nums = NULL;

  // read in the system
  setupSystemStructure(&sys, fileName, putIntoOrder);

  // setup nums - double precision is okay
  setupNums(&nums, sys.numNumbers, 64, 0);

  // check exponentiation using sys
  checkExp(&sys, nums, allowNegExp);

  // differentiate and create arr.out
  diff_sys_array(&sys, "arr.out");

  // clear nums
  clearNums(&nums, sys.numNumbers);

  // clear sys
  clearSystemStructure(&sys);

  return;
}

void setup_memoryLocations_sys(memoryLocations *memoryLoc, systemStruct *sys)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup memoryLoc from sys                               *
\***************************************************************/
{
  memoryLoc->numVars = sys->numVars;
  memoryLoc->varStart = sys->varsAddr;

  memoryLoc->numPathVars = sys->numPathVars;
  memoryLoc->pathvarStart = sys->pathVarsAddr;

  memoryLoc->numParams = sys->numParams;
  memoryLoc->paramStart = sys->paramAddr;
  memoryLoc->paramDerivStart = sys->funcAddr + sys->numFuncs;

  memoryLoc->numFuncs = sys->numFuncs;
  memoryLoc->funcStart = sys->funcAddr;
  memoryLoc->funcDerivVStart = sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars;
  memoryLoc->funcDerivPStart = sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars + sys->numVars * sys->numFuncs;

  memoryLoc->numConsts = sys->numConstants;
  memoryLoc->constStart = sys->constAddr;

  memoryLoc->numNums = sys->numNumbers;
  memoryLoc->numStart = sys->numAddr;

  memoryLoc->numSubfuncs = sys->numSubfuncs;
  memoryLoc->subfuncStart = sys->subFuncAddr;
  memoryLoc->subfuncDerivVStart = sys->subFuncAddr + sys->numSubfuncs;
  memoryLoc->subfuncDerivPStart = sys->subFuncAddr + sys->numSubfuncs + sys->numVars * sys->numSubfuncs;

  return;
}

void copy_func_op_expOps(expOps *op, func_ops *funcOp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy funcOp to op                                      *
\***************************************************************/
{
  op->memLoc = funcOp->memLoc;
  op->op = funcOp->op;
  op->in[0] = funcOp->in[0];
  op->in[1] = funcOp->in[1];

  return;
}

void copy_func_ops_expArrayOps(expArrayOps *eOps, func_ops *funcOps, int numOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy funcOps to eOps                                   *
\***************************************************************/
{
  int i;

  initialize_expArrayOps(eOps);
  eOps->numOps = numOps;
  eOps->ops = (expOps *)bmalloc(numOps * sizeof(expOps));
  for (i = 0; i < numOps; i++)
    copy_func_op_expOps(&eOps->ops[i], &funcOps[i]);

  return;
}

void setupParamOps_func_ops(expArrayOps *paramOps, expArrayOps *paramDiffOps, func_ops *ops, int numOps, memoryLocations *memoryLoc, int *firstFreeMemLoc, int *memLocDerivs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup paramOps & paramDiffOps from ops                 *
\***************************************************************/
{
  int i, numOldOps = paramOps->numOps, zeroAddr = memoryLoc->numStart, oneAddr = memoryLoc->numStart + 1;
  int *outDeriv = NULL, deriv0, deriv1;
  expOps *tempOp = NULL;

  // setup memory for new operations in paramOps
  paramOps->numOps += numOps;
  paramOps->ops = (expOps *)brealloc(paramOps->ops, paramOps->numOps * sizeof(expOps));

  // loop over the operations - add on the operations and compute the deriv
  for (i = 0; i < numOps; i++)
  { // copy over the operation
    tempOp = &paramOps->ops[numOldOps + i];
    copy_func_op_expOps(tempOp, &ops[i]);
    // compute the derivative of this operation
    outDeriv = &memLocDerivs[tempOp->memLoc];
    deriv0 = memLocDerivs[tempOp->in[0]];
    if (isUnary(tempOp->op))
    { // setup deriv1
      deriv1 = memLocDerivs[tempOp->in[1]];
    }
    else
    { // setup deriv1
      deriv1 = -2;
    }
    derivMemLoc(paramDiffOps, tempOp, outDeriv, deriv0, deriv1, zeroAddr, oneAddr, firstFreeMemLoc);
  }

  // null out memory
  outDeriv = NULL;
  tempOp = NULL;

  return;
}

void setupParamOps_funcStruct(expArrayOps *paramOps, expArrayOps *paramDiffOps, funcStruct *paramStruct, int numParams, memoryLocations *memoryLoc, int *firstFreeMemLoc, int *memLocDerivs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup paramOps & paramDiffOps from paramStruct         *
\***************************************************************/
{
  int i;

  initialize_expArrayOps(paramOps);
  initialize_expArrayOps(paramDiffOps);

  for (i = 0; i < numParams; i++)
  { // initialize memLocDerivs
    initialize_memLocDerivs(memLocDerivs, memoryLoc);
    // set deriv of the path variable to 1
    memLocDerivs[memoryLoc->pathvarStart] = -1;
    // copy the ops for this param to paramOps & setup paramDiffOps
    setupParamOps_func_ops(paramOps, paramDiffOps, paramStruct[i].ops, paramStruct[i].num_ops, memoryLoc, firstFreeMemLoc, memLocDerivs);

  }

  return;
}

void setupFuncOps_funcStruct(expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, int *funcOrder, int *subFuncOrder, int **startSubfuncs, int **endSubfuncs, int **startFuncs, int **endFuncs, int **startJvSub, int **endJvSub, int **startJv, int **endJv, funcStruct *funcs, int numFuncs, funcStruct *subFuncs, int numSubfuncs, memoryLocations *memoryLoc, int *firstFreeMemLoc, int *memLocDerivs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup funcOps, JvOps & JpOps from (sub)funcStruct      *
\***************************************************************/
{
  int i, j, currMemFunc = 0, currFuncOp = 0, currMemJv = 0, currJvOp = 0, currMemJp = 0, currJpOp = 0;
  int numVars = memoryLoc->numVars, numParams = memoryLoc->numParams;
  int isJv, totalJvOps = 0, totalJpOps = 0;
  int zeroAddr = memoryLoc->numStart, oneAddr = memoryLoc->numStart + 1;
  expArrayOps *func_temp = (expArrayOps *)bmalloc(numFuncs * sizeof(expArrayOps));
  expArrayOps *subFunc_temp = (expArrayOps *)bmalloc(numSubfuncs * sizeof(expArrayOps));
  expArrayOps **Jv_func_temp = (expArrayOps **)bmalloc(numFuncs * sizeof(expArrayOps *));
  expArrayOps **Jv_subFunc_temp = (expArrayOps **)bmalloc(numSubfuncs * sizeof(expArrayOps *));
  expArrayOps **Jp_func_temp = (expArrayOps **)bmalloc(numFuncs * sizeof(expArrayOps *));
  expArrayOps **Jp_subFunc_temp = (expArrayOps **)bmalloc(numSubfuncs * sizeof(expArrayOps *));

  // verify ordering
  for (i = 0; i < numSubfuncs; i++)
    if (subFuncOrder[i] != i)
    {
      printf("ERROR: Invalid ordering of the subfunctions!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
  for (i = 0; i < numFuncs; i++)
    if (funcOrder[i] != i + numSubfuncs)
    {
      printf("ERROR: Invalid ordering of the functions!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }

  // initialize
  initialize_expArrayOps(funcOps);
  initialize_expArrayOps(JvOps);
  initialize_expArrayOps(JpOps);

  for (i = 0; i < numFuncs; i++)
  {
    initialize_expArrayOps(&func_temp[i]);
    Jv_func_temp[i] = (expArrayOps *)bmalloc(numVars * sizeof(expArrayOps));
    for (j = 0; j < numVars; j++)
      initialize_expArrayOps(&Jv_func_temp[i][j]);
    Jp_func_temp[i] = (expArrayOps *)bmalloc(numParams * sizeof(expArrayOps));
    for (j = 0; j < numParams; j++)
      initialize_expArrayOps(&Jp_func_temp[i][j]);
  }
  for (i = 0; i < numSubfuncs; i++)
  {
    initialize_expArrayOps(&subFunc_temp[i]);
    Jv_subFunc_temp[i] = (expArrayOps *)bmalloc(numVars * sizeof(expArrayOps));
    for (j = 0; j < numVars; j++)
      initialize_expArrayOps(&Jv_subFunc_temp[i][j]);
    Jp_subFunc_temp[i] = (expArrayOps *)bmalloc(numParams * sizeof(expArrayOps));
    for (j = 0; j < numParams; j++)
      initialize_expArrayOps(&Jp_subFunc_temp[i][j]);
  }

  // allocate other data
  *startSubfuncs = (int *)bmalloc(numSubfuncs * sizeof(int));
  *endSubfuncs = (int *)bmalloc(numSubfuncs * sizeof(int));
  *startFuncs = (int *)bmalloc(numFuncs * sizeof(int));
  *endFuncs = (int *)bmalloc(numFuncs * sizeof(int));
  *startJvSub = (int *)bmalloc(numSubfuncs * sizeof(int));
  *endJvSub = (int *)bmalloc(numSubfuncs * sizeof(int));
  *startJv = (int *)bmalloc(numFuncs * sizeof(int));
  *endJv = (int *)bmalloc(numFuncs * sizeof(int));

  // setup funcOps - first subfunctions & then functions
  for (i = 0; i < numSubfuncs; i++)
  { // copy this subfunction to funcOps
    (*startSubfuncs)[i] = currMemFunc;
    funcOps->numOps += subFuncs[i].num_ops;
    funcOps->ops = (expOps *)brealloc(funcOps->ops, funcOps->numOps * sizeof(expOps));
    for (j = 0; j < subFuncs[i].num_ops; j++)
      copy_func_op_expOps(&funcOps->ops[currFuncOp + j], &subFuncs[i].ops[j]);

    // udpate memory count
    (*endSubfuncs)[i] = currMemFunc += countMem(&funcOps->ops[currFuncOp], subFuncs[i].num_ops);
    currFuncOp = funcOps->numOps;
  }
  for (i = 0; i < numFuncs; i++)
  { // copy this function to funcOps
    (*startFuncs)[i] = currMemFunc;
    funcOps->numOps += funcs[i].num_ops;
    funcOps->ops = (expOps *)brealloc(funcOps->ops, funcOps->numOps * sizeof(expOps));
    for (j = 0; j < funcs[i].num_ops; j++)
      copy_func_op_expOps(&funcOps->ops[currFuncOp + j], &funcs[i].ops[j]);

    // udpate memory count
    (*endFuncs)[i] = currMemFunc += countMem(&funcOps->ops[currFuncOp], funcs[i].num_ops);
    currFuncOp = funcOps->numOps;
  }

  // point to the proper operations using func_temp & subFunc_temp
  currFuncOp = 0;
  for (i = 0; i < numSubfuncs; i++)
  { // subfunction
    subFunc_temp[i].ops = &funcOps->ops[currFuncOp];
    // increment
    currFuncOp += subFunc_temp[i].numOps = subFuncs[i].num_ops;
  }
  for (i = 0; i < numFuncs; i++)
  { // function
    func_temp[i].ops = &funcOps->ops[currFuncOp];
    // increment
    currFuncOp += func_temp[i].numOps = funcs[i].num_ops;
  }

  // loop over the variables to compute the derivatives
  isJv = 1;
  for (j = 0; j < numVars; j++)
  { // initialize derivs
    initialize_memLocDerivs(memLocDerivs, memoryLoc);
    // jth variable
    memLocDerivs[memoryLoc->varStart + j] = -1;

    // loop over the subfunctions
    for (i = 0; i < numSubfuncs; i++)
    { // compute the derivative of ith subfunction w.r.t. jth variable
      setupDiffOps(&Jv_subFunc_temp[i][j], &subFunc_temp[i], memoryLoc, firstFreeMemLoc, memLocDerivs, INLINESUBFUNCTIONTYPE, i, isJv, j, 1);

      // increment the number of operations
      totalJvOps += Jv_subFunc_temp[i][j].numOps;
    }

    // loop over the functions
    for (i = 0; i < numFuncs; i++)
    { // compute the derivative of ith function w.r.t. jth variable
      setupDiffOps(&Jv_func_temp[i][j], &func_temp[i], memoryLoc, firstFreeMemLoc, memLocDerivs, FUNCTIONTYPE, i, isJv, j, 1);

      // increment the number of operations
      totalJvOps += Jv_func_temp[i][j].numOps;
    }
  }

  // setup Jv
  JvOps->numOps = totalJvOps;
  JvOps->ops = (expOps *)bmalloc(totalJvOps * sizeof(expOps));
  for (i = 0; i < numSubfuncs; i++)
  { // copy over the operations for this subfunction
    (*startJvSub)[i] = currMemJv;

    for (j = 0; j < numVars; j++)
    { // copy Jv_subFunc_temp[i][j] to JvOps - update currJvOp & currMemJv
      currMemJv += copyOps(&JvOps->ops[currJvOp], Jv_subFunc_temp[i][j].ops, Jv_subFunc_temp[i][j].numOps, zeroAddr, oneAddr);

      currJvOp += Jv_subFunc_temp[i][j].numOps;
      // clear Jv_subFunc_temp[i][j]
      clear_expArrayOps(&Jv_subFunc_temp[i][j]);
    }
    // free Jv_subFunc_temp[i]
    free(Jv_subFunc_temp[i]);

    // store the end of this function
    (*endJvSub)[i] = currMemJv;
  }
  for (i = 0; i < numFuncs; i++)
  { // copy over the operations for this function
    (*startJv)[i] = currMemJv;

    for (j = 0; j < numVars; j++)
    { // copy Jv_func_temp[i][j] to JvOps - update currJvOp & currMemJv
      currMemJv += copyOps(&JvOps->ops[currJvOp], Jv_func_temp[i][j].ops, Jv_func_temp[i][j].numOps, zeroAddr, oneAddr);
      currJvOp += Jv_func_temp[i][j].numOps;
      // clear Jv_func_temp[i][j]
      clear_expArrayOps(&Jv_func_temp[i][j]);
    }
    // free Jv_func_temp[i]
    free(Jv_func_temp[i]);

    // store the end of this function
    (*endJv)[i] = currMemJv;
  }

  // loop over the parameters to compute the derivatives
  isJv = 0;
  for (j = 0; j < numParams; j++)
  { // initialize derivs
    initialize_memLocDerivs(memLocDerivs, memoryLoc);
    // jth parameter
    memLocDerivs[memoryLoc->paramStart + j] = -1;

    // loop over the subfunctions
    for (i = 0; i < numSubfuncs; i++)
    { // compute the derivative of ith subfunction w.r.t. jth parameter
      setupDiffOps(&Jp_subFunc_temp[i][j], &subFunc_temp[i], memoryLoc, firstFreeMemLoc, memLocDerivs, INLINESUBFUNCTIONTYPE, i, isJv, j, 1);

      // increment the number of operations
      totalJpOps += Jp_subFunc_temp[i][j].numOps;
    }

    // loop over the functions
    for (i = 0; i < numFuncs; i++)
    { // compute the derivative of ith function w.r.t. jth parameter
      setupDiffOps(&Jp_func_temp[i][j], &func_temp[i], memoryLoc, firstFreeMemLoc, memLocDerivs, FUNCTIONTYPE, i, isJv, j, 1);

      // increment the number of operations
      totalJpOps += Jp_func_temp[i][j].numOps;
    }
  }

  // setup Jp
  JpOps->numOps = totalJpOps;
  JpOps->ops = (expOps *)bmalloc(totalJpOps * sizeof(expOps));
  for (i = 0; i < numSubfuncs; i++)
  { // copy over the operations for this subfunction
    for (j = 0; j < numParams; j++)
    { // copy Jp_subFunc_temp[i][j] to JpOps - update currJpOp & currMemJp
      currMemJp += copyOps(&JpOps->ops[currJpOp], Jp_subFunc_temp[i][j].ops, Jp_subFunc_temp[i][j].numOps, zeroAddr, oneAddr);
      currJpOp += Jp_subFunc_temp[i][j].numOps;
      // clear Jp_subFunc_temp[i][j]
      clear_expArrayOps(&Jp_subFunc_temp[i][j]);
    }
    // free Jp_subFunc_temp[i]
    free(Jp_subFunc_temp[i]);
  }
  for (i = 0; i < numFuncs; i++)
  { // copy over the operations for this function
    for (j = 0; j < numParams; j++)
    { // copy Jp_func_temp[i][j] to JpOps - update currJpOp & currMemJp
      currMemJp += copyOps(&JpOps->ops[currJpOp], Jp_func_temp[i][j].ops, Jp_func_temp[i][j].numOps, zeroAddr, oneAddr);
      currJpOp += Jp_func_temp[i][j].numOps;
      // clear Jp_func_temp[i][j]
      clear_expArrayOps(&Jp_func_temp[i][j]);
    }
    // free Jp_func_temp[i]
    free(Jp_func_temp[i]);
  }

  // clear func_temp & subFunc_temp
  for (i = 0; i < numFuncs; i++)
    func_temp[i].ops = NULL;
  free(func_temp);
  for (i = 0; i < numSubfuncs; i++)
    subFunc_temp[i].ops = NULL;
  free(subFunc_temp);

  // clear J*_func_temp & J*_subFunc_temp
  free(Jv_func_temp);
  free(Jv_subFunc_temp);
  free(Jp_func_temp);
  free(Jp_subFunc_temp);

  return;
}

void diff_sys_array(systemStruct *sys, char *fileName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: turns sys into "arr.out"                               *
\***************************************************************/
{
  int firstFreeMemLoc = sys->firstFreeMemLoc, *memLocDerivs = (int *)bmalloc(sys->firstFreeMemLoc * sizeof(int));
  int *startSubfuncs = NULL, *endSubfuncs = NULL, *startFuncs = NULL, *endFuncs = NULL, *startJvSub = NULL, *endJvSub = NULL, *startJv = NULL, *endJv = NULL;
  memoryLocations memoryLoc;
  expArrayOps constOps, paramOps, paramDiffOps, funcOps, JvOps, JpOps;

  FILE *OUT = fopen("arr.out", "w");

  // setup memoryLoc from sys
  setup_memoryLocations_sys(&memoryLoc, sys);

  // setup constOps
  copy_func_ops_expArrayOps(&constOps, sys->updateOps, sys->numUpdate);

  // setup paramOps &  paramDiffOps
  setupParamOps_funcStruct(&paramOps, &paramDiffOps, sys->params, sys->numParams, &memoryLoc, &firstFreeMemLoc, memLocDerivs);

  // setup funcOps, JvOps, JpOps & all other data associated with it
  setupFuncOps_funcStruct(&funcOps, &JvOps, &JpOps, sys->funcOrder, sys->subFuncOrder, &startSubfuncs, &endSubfuncs, &startFuncs, &endFuncs, &startJvSub, &endJvSub, &startJv, &endJv, sys->funcs, sys->numFuncs, sys->subFuncs, sys->numSubfuncs, &memoryLoc, &firstFreeMemLoc, memLocDerivs);

  // setup arr.out
  setupArrOut(OUT, &constOps, &paramOps, &paramDiffOps, &funcOps, &JvOps, &JpOps, sys->funcOrder, sys->subFuncOrder, startSubfuncs, endSubfuncs, startFuncs, endFuncs, startJvSub, endJvSub, startJv, endJv, &memoryLoc, firstFreeMemLoc, sys->numVarGps, sys->varGpSizes, sys->randIndex);

  fclose(OUT);

  // clear memory
  free(memLocDerivs);
  free(startSubfuncs); free(endSubfuncs);
  free(startFuncs); free(endFuncs);
  free(startJvSub); free(endJvSub);
  free(startJv); free(endJv);
  clear_expArrayOps(&constOps);
  clear_expArrayOps(&paramOps);
  clear_expArrayOps(&paramDiffOps);
  clear_expArrayOps(&funcOps);
  clear_expArrayOps(&JvOps);
  clear_expArrayOps(&JpOps);

  return;
}





