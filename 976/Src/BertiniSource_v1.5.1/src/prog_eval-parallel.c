// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bertini.h"

// this file contains the functions needed to evaluate the straight line program

// global variables
_comp_d  **mem_d  = NULL;
_comp_mp **mem_mp = NULL;
int *size_d  = NULL;  // size of mem_d
int *size_mp = NULL;  // size of mem_mp
int *mem_needs_init_d  = NULL; // determine if mem_d has been initialized
int *mem_needs_init_mp = NULL; // determine if mem_mp has been initialized

// initialize the varaibles used in evalProg - must be called in a serial (w.r.t OpenMP) portion of the code - i.e. not in an OpenMP parallel region
void initEvalProg(int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes the sturctures based on memSize & MPType   *
\***************************************************************/
{
  int i, max = max_threads();

  // allocate memory based on the max number of threads & the MPType
  if (MPType == 0 || MPType == 2)
  { // double precision needed
    mem_needs_init_d = (int *)bmalloc(max * sizeof(int));
    size_d = (int *)bmalloc(max * sizeof(int));
    mem_d  = (_comp_d **)bmalloc(max * sizeof(_comp_d *));

    for (i = 0; i < max; i++)
    {
      mem_d[i] = NULL;
      size_d[i] = 0; // the mem_d[i] structure has no elements in it
      mem_needs_init_d[i] = 1; // saying that mem_d[i] has not been initialized by evalProg_d
    }
  }
  // notice that this is not an if-else structure based on AMP needing both
  if (MPType == 1 || MPType == 2)
  { // multi precision needed
    mem_needs_init_mp = (int *)bmalloc(max * sizeof(int));
    size_mp = (int *)bmalloc(max * sizeof(int));
    mem_mp  = (_comp_mp **)bmalloc(max * sizeof(_comp_mp *));

    for (i = 0; i < max; i++)
    {
      mem_mp[i] = NULL;
      size_mp[i] = 0; // the mem_mp[i] structure has no elements in it
      mem_needs_init_mp[i] = 1; // saying that mem_mp[i] has not been initialized by evalProg_mp
    }
  }

  return;
}

void freeEvalProg(int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: releases the sturctures based on MPType                *
\***************************************************************/
{
  int i, j, max = max_threads();

  if (MPType == 0 || MPType == 2)
  { // clear double precision, if needed

    if (mem_d != NULL)
    {
      for (i = max - 1; i >= 0; i--)
      {
        free(mem_d[i]);
      }
      free(mem_d);
      mem_d = NULL;
    }
    if (mem_needs_init_d != NULL)
    {
      free(mem_needs_init_d);
      mem_needs_init_d = NULL;
    }
    if (size_d != NULL)
    {
      free(size_d);
      size_d = NULL; 
    }
  }
  // notice that this is not an if-else structure based on AMP needing both
  if (MPType == 1 || MPType == 2)
  { // clear multi precision, if needed

    if (mem_mp != NULL)
    {
      for (i = max - 1; i >= 0; i--)
      {
        if (mem_needs_init_mp != NULL && !mem_needs_init_mp[i])
        { // mem_mp[i] has been initialized so it needs cleared
          for (j = size_mp[i] - 1; j >= 0; j--)
          {
            clear_mp(&mem_mp[i][j]);
          }
        }
        free(mem_mp[i]);
      }
      free(mem_mp);
      mem_mp = NULL;
    }
    if (mem_needs_init_mp != NULL)
    {
      free(mem_needs_init_mp);
      mem_needs_init_mp = NULL;
    }
    if (size_mp != NULL)
    {
      free(size_mp);
      size_mp = NULL;
    }
  }

  return;
}

int change_prec_prog(void const *ED, int new_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: changes the precision on Prog                          *
\***************************************************************/
{
  prog_t *Prog = (prog_t *)ED;

  Prog->precision = new_prec;

  Prog = NULL;

  return 0;
}

void evalInsts_d(_comp_d *mem_d, int *prog, int startOp, int endOp, int oid)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do the instructions from startOp to endOp              *
\***************************************************************/
{
  int i = startOp, j, k, l;
  double tempD, tempD2;

  while (i < endOp)
  {
    j = prog[i];
    if (j == '*')
    {
      j = prog[i+1]; k = prog[i+2]; l = prog[i+3];
      mul_d2(&mem_d[j], &mem_d[k], &mem_d[l], tempD);
      i += 4;
    }
    else if (j == '+')
    {
      add_d(&mem_d[prog[i+1]], &mem_d[prog[i+2]], &mem_d[prog[i+3]]);
      i += 4;
    }
    else if (j == '=')
    {
      set_d(&mem_d[prog[i+1]], &mem_d[prog[i+2]]);
      i += 3;
    }
    else if (j == '-')
    {
      sub_d(&mem_d[prog[i+1]], &mem_d[prog[i+2]], &mem_d[prog[i+3]]);
      i += 4;
    }
    else if (j == 'N')
    {
      neg_d(&mem_d[prog[i+1]], &mem_d[prog[i+2]]);
      i += 3;
    }
    else if (j == '^')
    {
      j = prog[i+1]; k = prog[i+2]; l = prog[i+3];
      exp_d2(&mem_d[j], &mem_d[k], mem_d[l].r, tempD);
      i += 4;
    }
    else if (j == 'C')
    {
      cos_d2(&mem_d[prog[i+1]], &mem_d[prog[i+2]], tempD, tempD2);
      i += 3;
    }
    else if (j == 'S')
    {
      sin_d2(&mem_d[prog[i+1]], &mem_d[prog[i+2]], tempD, tempD2);
      i += 3;
    }
    else if (j == 'X')
    {
      exponential_d2(&mem_d[prog[i+1]], &mem_d[prog[i+2]], tempD);
      i += 3;
    }
    else if (j == '/')
    {
      div_d2(&mem_d[prog[i+1]], &mem_d[prog[i+2]], &mem_d[prog[i+3]], tempD, tempD2);
      i += 4;
    }
    else
    { // immediately exit if something is wrong with the program
      printf("\nERROR: Not a valid operation (%d) in evalProg_d!\n", j);
      bexit(ERROR_CONFIGURATION);
    }
  }

  return;
}

int evalProg_d_void(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate evalProg but with void at end                 *
\***************************************************************/
{
  prog_t *Prog = (prog_t *)ED;

  evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, Prog);

  Prog = NULL;

  return 0;
}

int evalProg_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, oid = thread_num();

  // check to see if the memory has been initialized before
  if (mem_needs_init_d[oid])
  { // allocate memory
    mem_d[oid] = (_comp_d *)bmalloc(Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // initialize the memory
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
    }

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }
  else if (size_d[oid] != Prog->memSize) // check to see if the memory is the same
  { // reallocate memory
    mem_d[oid] = (_comp_d *)brealloc(mem_d[oid], Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // initialize the memory
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
    }

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }

  // setup the numbers
  begin = Prog->numAddr;
  end = begin + Prog->numNums;
  for (i = begin; i < end; i++)
  {
    set_double_d(&mem_d[oid][i], mpq_get_d(Prog->nums[i - begin].rat), 0);
  }
  // setup I
  set_double_d(&mem_d[oid][Prog->IAddr], 0, 1);
  // setup Pi
  set_double_d(&mem_d[oid][Prog->IAddr + 1], M_PI, 0);

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_d(&mem_d[oid][i], &vars->coord[i - begin]);
  }

  // copy in the value of the path variable
  set_d(&mem_d[oid][Prog->inpPathVars], pathVars);

  // do the evaluations
  evalInsts_d(mem_d[oid], Prog->prog, 0, Prog->size, oid);

  // Gather the output
  change_size_vec_d(funcVals, Prog->numFuncs);
  change_size_point_d(parVals, Prog->numPars);
  change_size_vec_d(parDer, Prog->numPars);
  change_size_mat_d(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_d(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows     = Jp->rows = Prog->numFuncs;
  parVals->size  = parDer->size = Jp->cols = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_d(&funcVals->coord[i], &mem_d[oid][Prog->evalFuncs + i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_d(&Jv->entry[i][j], &mem_d[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_d(&Jp->entry[i][j], &mem_d[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_d(&parVals->coord[i], &mem_d[oid][Prog->evalPars+i]);
    set_d(&parDer->coord[i], &mem_d[oid][Prog->evalDPars+i]);
  }

  return 0;
}

int evalProg_eff_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog, int startFuncNum, int endFuncNum, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, skipUpdate = 1, oid = thread_num();
  int *subFuncsToEval = (int *)bmalloc(Prog->numSubfuncs * sizeof(int));

  // check to see if the memory has been initialized before
  if (mem_needs_init_d[oid])
  { // need to do update steps
    skipUpdate = 0;

    // allocate memory
    mem_d[oid] = (_comp_d *)bmalloc(Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // set all location to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mem_d[oid][i].r = mpq_get_d(Prog->nums[i - begin].rat);
    }

    // setup I
    set_double_d(&mem_d[oid][Prog->IAddr], 0, 1);
    // setup Pi
    set_double_d(&mem_d[oid][Prog->IAddr + 1], M_PI, 0);

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }
  else if (size_d[oid] != Prog->memSize) // check to see if the memory is the same
  { // need to do update steps
    skipUpdate = 0;

    // reallocate memory
    mem_d[oid] = (_comp_d *)brealloc(mem_d[oid], Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // set all location to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mem_d[oid][i].r = mpq_get_d(Prog->nums[i - begin].rat);
    }

    // setup I
    set_double_d(&mem_d[oid][Prog->IAddr], 0, 1);
    // setup Pi
    set_double_d(&mem_d[oid][Prog->IAddr + 1], M_PI, 0);

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_d(&mem_d[oid][i], &vars->coord[i - begin]);
  }
 
  // copy in the value of the path variables
  begin = Prog->inpPathVars;
  end = begin + Prog->numPathVars;
  for (i = begin; i < end; i++)
  {
    set_d(&mem_d[oid][i], pathVars);
  }

  // setup the subfunctions that need to be evaluated
  for (i = 0; i < Prog->numSubfuncs; i++)
  {
    subFuncsToEval[i] = 0;
    for (j = startFuncNum; j < endFuncNum; j++)
      if (subFuncsBelow[j][i])
      { // subfunction i needs to be evaluated for function j
        subFuncsToEval[i] = 1;
        j = endFuncNum; // end loop
      }
  }

  // see where we need to start at (do update or not)
  if (skipUpdate)
    i = Prog->numInstAtEndUpdate; 
  else
    i = 0;

  // evaluate up to Params
  evalInsts_d(mem_d[oid], Prog->prog, i, Prog->numInstAtEndParams, oid);

  // evaluate the subfunctions that are needed
  for (j = 0; j < Prog->numSubfuncs; j++)
    if (subFuncsToEval[j])
    { // evaluate the jth subfunction
      evalInsts_d(mem_d[oid], Prog->prog, startSub[j], endSub[j], oid);
    }

  // evaluate the functions that are needed
  for (j = startFuncNum; j < endFuncNum; j++)
  { // evaluate the jth function
    evalInsts_d(mem_d[oid], Prog->prog, startFunc[j], endFunc[j], oid);
  }

  // evaluate up to Param Derivs
  evalInsts_d(mem_d[oid], Prog->prog, Prog->numInstAtEndFnEval, Prog->numInstAtEndPDeriv, oid);

  // evaluate the subfunction derivatives that are needed
  for (j = 0; j < Prog->numSubfuncs; j++)
    if (subFuncsToEval[j])
    { // evaluate the jth subfunction derivatives
      evalInsts_d(mem_d[oid], Prog->prog, startJvsub[j], endJvsub[j], oid);
    }
  
  // evaluate the function derivatives that are needed
  for (j = startFuncNum; j < endFuncNum; j++)
  { // evaluate the jth function derivatives
    evalInsts_d(mem_d[oid], Prog->prog, startJv[j], endJv[j], oid);
  }

  // evaluate the param derivatives
  evalInsts_d(mem_d[oid], Prog->prog, Prog->numInstAtEndJvEval, Prog->size, oid);

  // Gather the output.
  change_size_vec_d(funcVals, Prog->numFuncs);
  change_size_point_d(parVals, Prog->numPars);
  change_size_vec_d(parDer, Prog->numPars);
  change_size_mat_d(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_d(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows     = Jp->rows = Prog->numFuncs;
  parVals->size  = parDer->size = Jp->cols = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_d(&funcVals->coord[i], &mem_d[oid][Prog->evalFuncs + i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_d(&Jv->entry[i][j], &mem_d[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_d(&Jp->entry[i][j], &mem_d[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_d(&parVals->coord[i], &mem_d[oid][Prog->evalPars+i]);
    set_d(&parDer->coord[i], &mem_d[oid][Prog->evalDPars+i]);
  }

  free(subFuncsToEval);

  return 0;
}

int evalProg_d_std(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, skipUpdate = 1, oid = thread_num();

  // check to see if the memory has been initialized before
  if (mem_needs_init_d[oid])
  { // need to do update steps
    skipUpdate = 0;

    // allocate memory
    mem_d[oid] = (_comp_d *)bmalloc(Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // initialize the memory
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mem_d[oid][i].r = mpq_get_d(Prog->nums[i - begin].rat);
    }
    // setup I
    set_double_d(&mem_d[oid][Prog->IAddr], 0, 1);
    // setup Pi
    set_double_d(&mem_d[oid][Prog->IAddr + 1], M_PI, 0);

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }
  else if (size_d[oid] != Prog->memSize) // check to see if the memory is the same
  { // need to do update steps
    skipUpdate = 0;

    // reallocate memory
    mem_d[oid] = (_comp_d *)brealloc(mem_d[oid], Prog->memSize * sizeof(_comp_d));
    size_d[oid] = Prog->memSize; // size of memory

    // initialize the memory
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // set to 0
      set_zero_d(&mem_d[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mem_d[oid][i].r = mpq_get_d(Prog->nums[i - begin].rat);
    }
    // setup I
    set_double_d(&mem_d[oid][Prog->IAddr], 0, 1);
    // setup Pi
    set_double_d(&mem_d[oid][Prog->IAddr + 1], M_PI, 0);

    mem_needs_init_d[oid] = 0; // memory has been initialized
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_d(&mem_d[oid][i], &vars->coord[i - begin]);
  }

  // copy in the value of the path variables
  set_d(&mem_d[oid][i], pathVars);

  // Run the program.
  if (skipUpdate)
  { // we are able to skip the update steps since they have been done before
    i = Prog->numInstAtEndUpdate;
  }
  else
  { // need to do the update steps since we have made changes
    i = 0;
  }

  // do the evaluations
  evalInsts_d(mem_d[oid], Prog->prog, i, Prog->size, oid);

  // Gather the output
  change_size_vec_d(funcVals, Prog->numFuncs);
  change_size_point_d(parVals, Prog->numPars);
  change_size_vec_d(parDer, Prog->numPars);
  change_size_mat_d(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_d(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows     = Jp->rows = Prog->numFuncs;
  parVals->size  = parDer->size = Jp->cols = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_d(&funcVals->coord[i], &mem_d[oid][Prog->evalFuncs + i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_d(&Jv->entry[i][j], &mem_d[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_d(&Jp->entry[i][j], &mem_d[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_d(&parVals->coord[i], &mem_d[oid][Prog->evalPars+i]);
    set_d(&parDer->coord[i], &mem_d[oid][Prog->evalDPars+i]);
  }

  return 0;
}

void evalInsts_mp(_comp_mp *mem_mp, int *prog, int startOp, int endOp, int oid)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: do the instructions from startOp to endOp              *
\***************************************************************/
{
  int i = startOp, j;
  while (i < endOp)
  {
    j = prog[i];
    if (j == '*')
    {
      mul_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], &mem_mp[prog[i+3]], oid);
      i += 4;
    }
    else if (j == '+')
    {
      add_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], &mem_mp[prog[i+3]]);
      i += 4;
    }
    else if (j == '=')
    {
      set_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]]);
      i += 3;
    }
    else if (j == '-')
    {
      sub_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], &mem_mp[prog[i+3]]);
      i += 4;
    } 
    else if (j == 'N')
    {
      neg_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]]);
      i += 3;
    }
    else if (j == '^')
    {
      exp_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], mem_mp[prog[i+3]].r, oid);
      i += 4;
    }
    else if (j == 'C')
    {
      cos_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], oid);
      i += 3;
    }
    else if (j == 'S')
    {
      sin_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], oid);
      i += 3;
    }
    else if (j == 'X')
    {
      exponential_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], oid);
      i += 3;
    }
    else if (j == '/')
    {
      div_omp_mp(&mem_mp[prog[i+1]], &mem_mp[prog[i+2]], &mem_mp[prog[i+3]], oid);
      i += 4;
    }
    else
    { // immediately exit if something is wrong with the program
      printf("\nERROR: Not a valid operation (%d) in evalProg_mp!\n", prog[i]);
      bexit(ERROR_CONFIGURATION);
    }
  }

  return;
}

int evalProg_mp_void(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate evalProg but with void at end                 *
\***************************************************************/
{
  prog_t *Prog = (prog_t *)ED;

  evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, Prog);

  Prog = NULL;

  return 0;
}

int evalProg_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, oid = thread_num(), curr_prec = Prog->precision;

  // check to see if the memory has been initialized before
  if (mem_needs_init_mp[oid])
  { // allocate memory
    mem_mp[oid] = (_comp_mp *)bmalloc(Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory & set to 0
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
    }

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (size_mp[oid] != Prog->memSize)  // check to see if the memory is the same
  { // clear old memory
    for (i = 0; i < size_mp[oid]; i++)
      clear_mp(&mem_mp[oid][i]);

    // allocate more memory
    mem_mp[oid] = (_comp_mp *)brealloc(mem_mp[oid], Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory & set to 0
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
    }

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (mpfr_get_prec(mem_mp[oid][0].r) != curr_prec) // precision is not correct
  { // set precision & set to zero
    for (i = 0; i < Prog->memSize; i++)
    { // set the precision
      setprec_mp(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
    }
  }

  // setup the numbers
  begin = Prog->numAddr;
  end = begin + Prog->numNums;
  for (i = begin; i < end; i++)
  { // setup real parts
    mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);

    // imaginary parts are always zero
    mpf_set_ui(mem_mp[oid][i].i, 0);
  }
  
  // setup I
  mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
  mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
  // setup Pi
  mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
  mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_mp(&mem_mp[oid][i], &vars->coord[i - begin]);
  }
  // copy in the value of the path variable
  set_mp(&mem_mp[oid][Prog->inpPathVars], pathVars);

  // do the evaluations
  evalInsts_mp(mem_mp[oid], Prog->prog, 0, Prog->size, oid);

  // Gather the output.
  change_size_vec_mp(funcVals, Prog->numFuncs);
  change_size_point_mp(parVals, Prog->numPars);
  change_size_vec_mp(parDer, Prog->numPars);
  change_size_mat_mp(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_mp(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows = Jp->rows     = Prog->numFuncs;
  parVals->size  = Jp->cols = parDer->size = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_mp(&funcVals->coord[i], &mem_mp[oid][Prog->evalFuncs+i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_mp(&Jv->entry[i][j], &mem_mp[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_mp(&Jp->entry[i][j], &mem_mp[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_mp(&parVals->coord[i], &mem_mp[oid][Prog->evalPars+i]);
    set_mp(&parDer->coord[i], &mem_mp[oid][Prog->evalDPars+i]);
  }

  return 0;
}

int evalProg_eff_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog, int startFuncNum, int endFuncNum, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, skipUpdate = 1, oid = thread_num(), curr_prec = Prog->precision;
  int *subFuncsToEval = (int *)bmalloc(Prog->numSubfuncs * sizeof(int));

  // check to see if the memory has been initialized before
  if (mem_needs_init_mp[oid])
  { // need to do update steps
    skipUpdate = 0;

    // allocate memory
    mem_mp[oid] = (_comp_mp *)bmalloc(Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (size_mp[oid] != Prog->memSize)  // check to see if the memory is the same
  { // need to do update steps
    skipUpdate = 0;

    // clear old memory
    for (i = 0; i < size_mp[oid]; i++)
      clear_mp(&mem_mp[oid][i]);

    // allocate more memory
    mem_mp[oid] = (_comp_mp *)brealloc(mem_mp[oid], Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (mpfr_get_prec(mem_mp[oid][0].r) != curr_prec) // precision is not correct
  { // need to do update steps
    skipUpdate = 0;

    // set the precision, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      setprec_mp(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_mp(&mem_mp[oid][i], &vars->coord[i - begin]);
  }

  // copy in the value of the path variables
  begin = Prog->inpPathVars;
  end = begin + Prog->numPathVars;
  for (i = begin; i < end; i++)
  {
    set_mp(&mem_mp[oid][i], pathVars);
  }

  // setup the subfunctions that need to be evaluated
  for (i = 0; i < Prog->numSubfuncs; i++)
  {
    subFuncsToEval[i] = 0;
    for (j = startFuncNum; j < endFuncNum; j++)
      if (subFuncsBelow[j][i])
      { // subfunction i needs to be evaluated for function j
        subFuncsToEval[i] = 1;
        j = endFuncNum; // end loop
      }
  }

  // see where we need to start at (do update or not)
  if (skipUpdate)
    i = Prog->numInstAtEndUpdate;
  else
    i = 0;

  // evaluate up to Params
  evalInsts_mp(mem_mp[oid], Prog->prog, i, Prog->numInstAtEndParams, oid);

  // evaluate the subfunctions that are needed
  for (j = 0; j < Prog->numSubfuncs; j++)
    if (subFuncsToEval[j])
    { // evaluate the jth subfunction
      evalInsts_mp(mem_mp[oid], Prog->prog, startSub[j], endSub[j], oid);
    }

  // evaluate the functions that are needed
  for (j = startFuncNum; j < endFuncNum; j++)
  { // evaluate the jth function
    evalInsts_mp(mem_mp[oid], Prog->prog, startFunc[j], endFunc[j], oid);
  }

  // evaluate up to Param Derivs
  evalInsts_mp(mem_mp[oid], Prog->prog, Prog->numInstAtEndFnEval, Prog->numInstAtEndPDeriv, oid);

  // evaluate the subfunction derivatives that are needed
  for (j = 0; j < Prog->numSubfuncs; j++)
    if (subFuncsToEval[j])
    { // evaluate the jth subfunction derivatives
      evalInsts_mp(mem_mp[oid], Prog->prog, startJvsub[j], endJvsub[j], oid);
    }

  // evaluate the function derivatives that are needed
  for (j = startFuncNum; j < endFuncNum; j++)
  { // evaluate the jth function derivatives
    evalInsts_mp(mem_mp[oid], Prog->prog, startJv[j], endJv[j], oid);
  }

  // evaluate the param derivatives
  evalInsts_mp(mem_mp[oid], Prog->prog, Prog->numInstAtEndJvEval, Prog->size, oid);

  // Gather the output.
  change_size_vec_mp(funcVals, Prog->numFuncs);
  change_size_point_mp(parVals, Prog->numPars);
  change_size_vec_mp(parDer, Prog->numPars);
  change_size_mat_mp(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_mp(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows     = Jp->rows = Prog->numFuncs;
  parVals->size  = parDer->size = Jp->cols = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_mp(&funcVals->coord[i], &mem_mp[oid][Prog->evalFuncs+i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_mp(&Jv->entry[i][j], &mem_mp[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_mp(&Jp->entry[i][j], &mem_mp[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_mp(&parVals->coord[i], &mem_mp[oid][Prog->evalPars+i]);
    set_mp(&parDer->coord[i], &mem_mp[oid][Prog->evalDPars+i]);
  }

  free(subFuncsToEval);

  return 0;
}

int evalProg_mp_std(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, begin, end, skipUpdate = 1, oid = thread_num(), curr_prec = Prog->precision;

  // check to see if the memory has been initialized before
  if (mem_needs_init_mp[oid])
  { // need to do update steps
    skipUpdate = 0;

    // allocate memory
    mem_mp[oid] = (_comp_mp *)bmalloc(Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0 
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (size_mp[oid] != Prog->memSize)  // check to see if the memory is the same
  { // need to do update steps
    skipUpdate = 0;

    // clear old memory
    for (i = 0; i < size_mp[oid]; i++)
      clear_mp(&mem_mp[oid][i]);

    // allocate more memory
    mem_mp[oid] = (_comp_mp *)brealloc(mem_mp[oid], Prog->memSize * sizeof(_comp_mp));
    size_mp[oid] = Prog->memSize; // size of memory

    // initialize the memory, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      init_mp2(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);

    mem_needs_init_mp[oid] = 0; // memory has been initialized
  }
  else if (mpfr_get_prec(mem_mp[oid][0].r) != curr_prec) // precision is not correct
  { // need to do update steps
    skipUpdate = 0;

    // set the precision, set to 0 & setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSize; i++)
    { // initialize
      setprec_mp(&mem_mp[oid][i], curr_prec);
      // set to 0
      set_zero_mp(&mem_mp[oid][i]);
      // see if we need to setup a number
      if (i >= begin && i < end)
        mpf_set_q(mem_mp[oid][i].r, Prog->nums[i - begin].rat);
    }

    // setup I
    mpf_set_ui(mem_mp[oid][Prog->IAddr].r, 0);
    mpf_set_ui(mem_mp[oid][Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(mem_mp[oid][Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(mem_mp[oid][Prog->IAddr + 1].i, 0);
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_mp(&mem_mp[oid][i], &vars->coord[i - begin]);
  }
  // copy in the value of the path variable
  set_mp(&mem_mp[oid][Prog->inpPathVars], pathVars);

  // Run the program.
  if (skipUpdate)
  { // we are able to skip the update steps since they have been done before
    i = Prog->numInstAtEndUpdate;
  }
  else
  { // need to do the update steps since we have made changes
    i = 0;
  }

  // do the evaluations
  evalInsts_mp(mem_mp[oid], Prog->prog, i, Prog->size, oid);

  // Gather the output.
  change_size_vec_mp(funcVals, Prog->numFuncs);
  change_size_point_mp(parVals, Prog->numPars);
  change_size_vec_mp(parDer, Prog->numPars);
  change_size_mat_mp(Jv, Prog->numFuncs, Prog->numVars);
  change_size_mat_mp(Jp, Prog->numFuncs, Prog->numPars);
  funcVals->size = Jv->rows = Jp->rows     = Prog->numFuncs;
  parVals->size  = Jp->cols = parDer->size = Prog->numPars;
  Jv->cols       = Prog->numVars;

  for (i = 0; i < Prog->numFuncs; i++)
  {
    set_mp(&funcVals->coord[i], &mem_mp[oid][Prog->evalFuncs+i]);
    for (j = 0; j < Prog->numVars; j++)
    {
      set_mp(&Jv->entry[i][j], &mem_mp[oid][Prog->evalJVars + i*(Prog->numVars) + j]);
    }
    for (j = 0; j < Prog->numPars; j++)
    {
      set_mp(&Jp->entry[i][j], &mem_mp[oid][Prog->evalJPars + i*(Prog->numPars) + j]);
    }
  }

  for (i = 0; i < Prog->numPars; i++)
  {
    set_mp(&parVals->coord[i], &mem_mp[oid][Prog->evalPars+i]);
    set_mp(&parDer->coord[i], &mem_mp[oid][Prog->evalDPars+i]);
  }

  return 0;
}


