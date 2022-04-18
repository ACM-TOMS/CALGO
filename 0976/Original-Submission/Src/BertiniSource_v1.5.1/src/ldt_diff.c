// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "diff.h"

void diff_forward_vars_subfuncs_old(systemStruct *sys, int currSubFunc, int currVar, int storeAddr, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_forward_vars_funcs_old(systemStruct *sys, int currFunc, int currVar, int storeAddr, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);

void diff_op_old(func_ops *ops, int diff_var, int *memLoc, int **derivAddr, int totalOpCount, int oneAddr, int *first_free_mem_loc, int *new_ops, func_ops **diff_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: generate the derivative of ops w.r.t. diff_var         *
\***************************************************************/
{
  int i, indexIn0, indexIn1, indexLoc;
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
    derivAddr[indexLoc][diff_var] = derivAddr[indexIn0][diff_var];
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still zero
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is -1 - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', oneAddr, -1);
    }
    else
    { // negate the old location - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', derivAddr[indexIn0][diff_var], -1);
    }
  }
  else if (op == 'S')
  { // derivative is cos(in[0])*deriv0
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still zero
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is cos(in[0]) - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'C', ops->in[0], -1);
    }
    else
    { // cos(in[0])*deriv0 - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 2;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, 'C', ops->in[0], -1);
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '*', *first_free_mem_loc, derivAddr[indexIn0][diff_var]);
      (*first_free_mem_loc)++;
    }
  }
  else if (op == 'C')
  { // derivative is -sin(in[0])*deriv0
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still zero
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is -sin(in[0]) - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 2;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, 'S', ops->in[0], -1);
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', *first_free_mem_loc, -1);
      (*first_free_mem_loc)++;
    }
    else
    { // -sin(in[0])*deriv0 - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 3;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, 'S', ops->in[0], -1);
      (*first_free_mem_loc)++;
      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', *first_free_mem_loc - 1, derivAddr[indexIn0][diff_var]);
      (*first_free_mem_loc)++;
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', *first_free_mem_loc - 1, -1);
    }
  }
  else if (op == 'X')
  { // derivative is exp(in[0])*deriv0 = memLoc*deriv0
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still zero
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is exp(in[0]) = memLoc
      derivAddr[indexLoc][diff_var] = ops->memLoc;
    }
    else
    { // exp(in[0])*deriv0 = memLoc*deriv0 - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '*', ops->memLoc, derivAddr[indexIn0][diff_var]);
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is just the derivative of in[1]
      derivAddr[indexLoc][diff_var] = derivAddr[indexIn1][diff_var];
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is 1 + diff(in[1])
      if (derivAddr[indexIn1][diff_var] == -2)
      { // derivative is just 'one'
        derivAddr[indexLoc][diff_var] = -1;
      }
      else if (derivAddr[indexIn1][diff_var] == -1)
      { // derivative is 1 + 1 - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '+', oneAddr, oneAddr);
      }
      else
      { // derivative is 1 + diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '+', oneAddr, derivAddr[indexIn1][diff_var]);
      }
    }
    else
    { // derivative is diff(in[0]) + diff(in[1])
      if (derivAddr[indexIn1][diff_var] == -2)
      { // derivative is just diff(in[0])
        derivAddr[indexLoc][diff_var] = derivAddr[indexIn0][diff_var];
      }
      else if (derivAddr[indexIn1][diff_var] == -1)
      { // derivative is diff(in[0]) + 1 - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '+', derivAddr[indexIn0][diff_var], oneAddr);
      }
      else
      { // derivative is diff(in[0]) + diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '+', derivAddr[indexIn0][diff_var], derivAddr[indexIn1][diff_var]);
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
    if (derivAddr[indexIn1][diff_var] == -2)
    { // derivative is just the derivative of in[0]
      derivAddr[indexLoc][diff_var] = derivAddr[indexIn0][diff_var];
    }
    else if (derivAddr[indexIn1][diff_var] == -1)
    { // derivative is diff(in[0]) - 1
      if (derivAddr[indexIn0][diff_var] == -2)
      { // derivative is -1 - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', oneAddr, -1);
      }
      else if (derivAddr[indexIn0][diff_var] == -1)
      { // derivative is 1 - 1 = 0 = 'zero'
        derivAddr[indexLoc][diff_var] = -2;
      }
      else
      { // derivative is diff(in[0]) - 1 - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '-', derivAddr[indexIn0][diff_var], oneAddr);
      }
    }
    else
    { // derivative is diff(in[0]) - diff(in[1])
      if (derivAddr[indexIn0][diff_var] == -2)
      { // derivative is just -diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], 'N', derivAddr[indexIn1][diff_var], -1);
      }
      else if (derivAddr[indexIn0][diff_var] == -1)
      { // derivative is 1 - diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '-', oneAddr, derivAddr[indexIn1][diff_var]);
      }
      else
      { // derivative is diff(in[0]) - diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '-', derivAddr[indexIn0][diff_var], derivAddr[indexIn1][diff_var]);
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still 0
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is 1 / in[1] - setup to new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '/', oneAddr, ops->in[1]);
    }
    else
    { // derivative is diff(in[0]) / in[1] - setup to new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 1;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
      copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '/', derivAddr[indexIn0][diff_var], ops->in[1]);
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is in[0] * diff(in[1])
      if (derivAddr[indexIn1][diff_var] == -2)
      { // derivative is 'zero'
        derivAddr[indexLoc][diff_var] = -2;
      }
      else if (derivAddr[indexIn1][diff_var] == -1)
      { // derivative is in[0]
        derivAddr[indexLoc][diff_var] = ops->in[0];
      }
      else
      { // derivative is in[0] * diff(in[1]) - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '*', ops->in[0], derivAddr[indexIn1][diff_var]);
      }
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is in[1] + diff(in[1]) * in[0]
      if (derivAddr[indexIn1][diff_var] == -2)
      { // derivative is just in[1]
        derivAddr[indexLoc][diff_var] = ops->in[1];
      }
      else if (derivAddr[indexIn1][diff_var] == -1)
      { // derivative is in[1] + in[0] - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '+', ops->in[0], ops->in[1]);
      }
      else
      { // derivative is in[1] + diff(in[1]) * in[0] - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 2;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[0], derivAddr[indexIn1][diff_var]);
        copyToOp(&(*diff_ops)[1], derivAddr[indexLoc][diff_var], '+', ops->in[1], *first_free_mem_loc);
        (*first_free_mem_loc)++;
      }
    }
    else
    { // derivative is diff(in[0]) * in[1] + diff(in[1]) * in[0]
      if (derivAddr[indexIn1][diff_var] == -2)
      { // derivative is just diff(in[0]) * in[1] - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 1;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], derivAddr[indexLoc][diff_var], '*', ops->in[1], derivAddr[indexIn0][diff_var]);
      }
      else if (derivAddr[indexIn1][diff_var] == -1)
      { // derivative is diff(in[0]) * in[1] + in[0] - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 2;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[1], derivAddr[indexIn0][diff_var]);
        copyToOp(&(*diff_ops)[1], derivAddr[indexLoc][diff_var], '+', ops->in[0], *first_free_mem_loc);
        (*first_free_mem_loc)++;
      }
      else
      { // derivative is diff(in[0]) * in[1] + diff(in[1]) * in[0] - setup in new location
        derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
        (*first_free_mem_loc)++;

        *new_ops = 3;
        *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));
        copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '*', ops->in[1], derivAddr[indexIn0][diff_var]);
        (*first_free_mem_loc)++;

        copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '*', ops->in[0], derivAddr[indexIn1][diff_var]);
        copyToOp(&(*diff_ops)[2], derivAddr[indexLoc][diff_var], '+', *first_free_mem_loc - 1, *first_free_mem_loc);
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
    if (derivAddr[indexIn0][diff_var] == -2)
    { // derivative is still 0
      derivAddr[indexLoc][diff_var] = -2;
    }
    else if (derivAddr[indexIn0][diff_var] == -1)
    { // derivative is in[1] * in[0] ^ (in[1] - 1) - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 3;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));

      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '-', ops->in[1], oneAddr);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '^', ops->in[0], *first_free_mem_loc - 1);
      copyToOp(&(*diff_ops)[2], derivAddr[indexLoc][diff_var], '*', *first_free_mem_loc, ops->in[1]);
      (*first_free_mem_loc)++;
    }
    else
    { // derivative is in[1] * (in[0] ^ (in[1] - 1)) * diff(in[0]) - setup in new location
      derivAddr[indexLoc][diff_var] = *first_free_mem_loc;
      (*first_free_mem_loc)++;

      *new_ops = 4;
      *diff_ops = (func_ops *)bmalloc(*new_ops * sizeof(func_ops));

      copyToOp(&(*diff_ops)[0], *first_free_mem_loc, '-', ops->in[1], oneAddr);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[1], *first_free_mem_loc, '^', ops->in[0], *first_free_mem_loc - 1);
      (*first_free_mem_loc)++;

      copyToOp(&(*diff_ops)[2], *first_free_mem_loc, '*', *first_free_mem_loc - 1, ops->in[1]);
      copyToOp(&(*diff_ops)[3], derivAddr[indexLoc][diff_var], '*', *first_free_mem_loc, derivAddr[indexIn0][diff_var]);
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

void diff_funcStruct_old(funcStruct *func, int endFuncCount, int currVar, int storeAddr, int zeroAddr, int oneAddr, int *firstFreeMemLoc, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate func w.r.t. 'currVar'                    *
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

    diff_op_old(&func->ops[i], currVar, memLoc, derivAddr, totalOpCount, oneAddr, firstFreeMemLoc, &op_count[i], &temp_ops[i]);

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
  if (derivAddr[endFuncCount - 1][currVar] == -2)
  { // derivative is zero
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', zeroAddr, -1);
  }
  else if (derivAddr[endFuncCount - 1][currVar] == -1)
  { // derivative is one
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', oneAddr, -1);
  }
  else
  { // derivative is at derivAddr
    copyToOp(&(*deriv_ops)[curr_loc], storeAddr, '=', derivAddr[endFuncCount - 1][currVar], -1);
  }

  // clear temp_ops & op_count
  free(temp_ops);
  free(op_count);

  return;
}

void diff_sys_vars_subfuncs_old(funcStruct *diff, systemStruct *sys, int currSubFunc, int *memLoc, int **derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate subfuncs w.r.t. vars - store to diff     *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, j, addr, count = 0;

  // setup memLoc for the variables, path variables, constants, numbers, parameters, subfunction & function operations
  // also setup derivAddr for variables (path variables, constants and numbers have 'zero' deriv already)
  count = 0;
  // variables
  for (i = 0; i < sys->numVars; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->varsAddr + i;
    for (j = 0; j < sys->numVars; j++)
    { // derivAddr[count][j + numPathVars] = 'one' if i == j, 'zero' if i != j
      derivAddr[count][j + sys->numPathVars] = (i == j ? -1 : -2);
    }
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

  // find the partial derivatives w.r.t. each variable
  for (i = 0; i < sys->numVars; i++)
  { // find where to store this partial derivative - after the subfunctions
    addr = sys->subFuncAddr + sys->numSubfuncs + currSubFunc * sys->numVars + i;

    // now we need to setup the derivatives for the 'currSubFunc' subfunction w.r.t. ith variable
    diff_forward_vars_subfuncs_old(sys, currSubFunc, i + sys->numPathVars, addr, memLoc, derivAddr, totalOpCount, &diff[i].num_ops, &diff[i].ops);
  }

  return;
}

void diff_forward_vars_subfuncs_old(systemStruct *sys, int currSubFunc, int currVar, int storeAddr, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate 'currSubFunc' w.r.t. 'currVar'           *
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
  diff_funcStruct_old(&sys->subFuncs[currSubFunc], count, currVar, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}

void diff_sys_vars_funcs_old(funcStruct *diff, systemStruct *sys, int currFunc, int *memLoc, int **derivAddr, int totalOpCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate funcs w.r.t. vars - store to diff        *
*  -2 means zero, -1 means one, >= 0 means that memory location *
\***************************************************************/
{
  int i, j, addr, count = 0;

  // setup memLoc for the variables, path variables, constants, numbers, parameters, subfunction & function operations
  // also setup derivAddr for variables (path variables, constants and numbers have 'zero' deriv already)
  count = 0;
  // variables
  for (i = 0; i < sys->numVars; i++)
  { // memLoc[count] = addr + i
    memLoc[count] = sys->varsAddr + i;
    for (j = 0; j < sys->numVars; j++)
    { // derivAddr[count][j + sys->numPathVars] = 'one' if i == j, 'zero' if i != j
      derivAddr[count][j + sys->numPathVars] = (i == j ? -1 : -2);
    }
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

  // find the partial derivatives w.r.t. each path variable
  for (i = 0; i < sys->numVars; i++)
  { // find where to store this partial derivative - after the functions & param derivs
    addr = sys->funcAddr + sys->numFuncs + sys->numParams * sys->numPathVars + currFunc * sys->numVars + i;

    // now we need to setup the derivatives for the 'currFunc' function w.r.t. ith variable
    diff_forward_vars_funcs_old(sys, currFunc, i + sys->numPathVars, addr, memLoc, derivAddr, totalOpCount, &diff[i].num_ops, &diff[i].ops);
  }

  return;
}

void diff_forward_vars_funcs_old(systemStruct *sys, int currFunc, int currVar, int storeAddr, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops)
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
  diff_funcStruct_old(&sys->funcs[currFunc], count, currVar, storeAddr, sys->numAddr, sys->numAddr + 1, &sys->firstFreeMemLoc, memLoc, derivAddr, totalOpCount, derivCount, deriv_ops);

  return;
}



