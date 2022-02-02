// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "ppParse.h"
#include "parallel.h"

void performOps(expArrayOps *ops, comp_d *mem);
void determineDependencies_func(int **sfDepend, int **fDepend, parseArray *Array, expArrayOps *func, expArrayOps *subFunc);
void reallocate_degArray(int ***degArray, int *degMem, int newMem, int numGps);
void degreeBound(double *degBound, int numFuncs, int funcStart, int numVarGps, int **degArray);
void coeffBound(double *cBound, expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, parseArray *Array);
void verify_exponents(expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, parseArray *Array);

void verifyRealExponent(comp_d exp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify real number                                     *
\***************************************************************/
{ // check this exp. operation
  double tol = 1e-15;

  // make sure that exp->i is '0'
  if (fabs(exp->i) > tol)
  { // print error message
    printf("ERROR: An exponent (%0.15e + %0.15e*I) was found that was nonreal!\n", exp->r, exp->i);
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void verifyBaseExponent(comp_d base, comp_d exp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify base is real if exp is not an integer           *
\***************************************************************/
{ // check this exp. operation
  int intExp = (int) exp->r;
  double tol = 1e-15;

  // make sure that exp->i is '0'
  if (fabs(exp->i) > tol)
  { // print error message
    printf("ERROR: An exponent (%0.15e + %0.15e*I) was found that was nonreal!\n", exp->r, exp->i);
    bexit(ERROR_INPUT_SYNTAX);
  }

  // determine if exp->r is an integer
  if (fabs(exp->r - intExp) > tol)
  { // exponent is not an integer - base must be real and nonnegative
    if (fabs(base->i) > tol || base->r < 0)
    { // print error message
      printf("ERROR: Invalid base with exponent %0.15e!\n", exp->r);
      bexit(ERROR_INPUT_SYNTAX);
    } 
  }

  return;
}

void verifyIntegerExponent(comp_d exp, char *funcName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify nonnegative integer                             *
\***************************************************************/
{ // check this exp. operation
  int retVal = 0, realExp = (int) exp->r;
  double tol = 1e-15;

  // make sure that exp->i is '0'
  if (fabs(exp->i) > tol)
  { // bad exponent!
    retVal = 1;
  }

  // make sure that exp->r is an integer
  if (!retVal && fabs(exp->r - realExp) > tol)
  { // not an integer
    retVal = 2;
  }

  // if needed, make sure that exp->r is non-negative
  if (!retVal && exp->r < 0)
  { // negative integer
    retVal = 3;
  }

  if (retVal == 1)
  { // print error message
    printf("ERROR: When seting up %s, an exponent was found that was nonreal!\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }
  else if (retVal == 2)
  { // print error message
    printf("ERROR: When setting up %s, an exponent was found that was not an integer!\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }
  else if (retVal == 3)
  { // print error message
    printf("ERROR: When setting up %s, an exponent was found that was negative!\n", funcName);
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void initialize_memLocDerivs(int *memLocDerivs, memoryLocations *memoryLoc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initailize memLocDerivs                                *
\***************************************************************/
{
  int i, currLoc = 0;
  int numVars = memoryLoc->numVars, numPathVars = memoryLoc->numPathVars, numParams = memoryLoc->numParams, numFuncs = memoryLoc->numFuncs;
  int numConsts = memoryLoc->numConsts, numNums = memoryLoc->numNums, numSubfuncs = memoryLoc->numSubfuncs;

  // variables
  currLoc = memoryLoc->varStart;
  for (i = 0; i < numVars; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // path variables
  currLoc = memoryLoc->pathvarStart;
  for (i = 0; i < numPathVars; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // parameters
  currLoc = memoryLoc->paramStart;
  for (i = 0; i < numParams; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // functions
  currLoc = memoryLoc->funcStart;
  for (i = 0; i < numFuncs; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // constants
  currLoc = memoryLoc->constStart;
  for (i = 0; i < numConsts; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // numbers
  currLoc = memoryLoc->numStart;
  for (i = 0; i < numNums; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  // subfunctions
  currLoc = memoryLoc->subfuncStart;
  for (i = 0; i < numSubfuncs; i++)
  {
    memLocDerivs[currLoc] = -2;
    currLoc++;
  }

  return;
}

int opUsesSubFunc(expOps *Op)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - not use defined subfunction, 1 - otherwise *
* NOTES: determine if Op uses defined subfunction call          *
\***************************************************************/
{
  int isSubfunc0 = 0, isSubfunc1 = 0, sfCall0 = 0, sfCall1 = 0, unary = isUnary(Op->op);

  // determine if in[0] is a defined subfunction
  isSubfunc0 = isDefinedSubfunc(&sfCall0, Op->in[0]);
  if (!unary)
    isSubfunc1 = isDefinedSubfunc(&sfCall1, Op->in[1]);

  return (isSubfunc0 || isSubfunc1);
}

void setupOp(expOps *Op, int memLoc, char op, int in0, int in1, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the operation in Op                              *
\***************************************************************/
{
  Op->memLoc = memLoc;
  Op->op = op;

  // verify in[0]
  if (in0 == -2)
    Op->in[0] = zeroAddr;
  else if (in0 == -1)
    Op->in[0] = oneAddr;
  else
    Op->in[0] = in0;

  if (!isUnary(op))
  {
    if (in1 == -2)
      Op->in[1] = zeroAddr;
    else if (in1 == -1)
      Op->in[1] = oneAddr;
    else
      Op->in[1] = in1;
  }
  else
    Op->in[1] = in1;

  return;
}

void copyOp(expOps *Out, expOps *In, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  Out->memLoc = In->memLoc;
  Out->op = In->op;

  // veryify in[0]
  if (In->in[0] == -2)
    Out->in[0] = zeroAddr;
  else if (In->in[0] == -1)
    Out->in[0] = oneAddr;
  else
    Out->in[0] = In->in[0];

  if (!isUnary(In->op))
  {
    if (In->in[1] == -2)
      Out->in[1] = zeroAddr;
    else if (In->in[1] == -1)
      Out->in[1] = oneAddr;
    else
      Out->in[1] = In->in[1];
  }
  else
    Out->in[1] = In->in[1];

  return;
}

int copyOps(expOps *Out, expOps *In, int numOps, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: memory locations used                          *
* NOTES: copy In to Out                                         *
\***************************************************************/
{
  int i, count = 0;
  
  for (i = 0; i < numOps; i++)
  { // copy op
    copyOp(&Out[i], &In[i], zeroAddr, oneAddr);
    // increment count
    count += 3 + !isUnary(Out[i].op);
  }

  return count;
}

void addOp(expArrayOps *Array, int memLoc, char op, int in0, int in1, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add operation to Array                                 *
\***************************************************************/
{
  int num = Array->numOps;

  Array->numOps++;
  Array->ops = (expOps *)brealloc(Array->ops, Array->numOps * sizeof(expOps));
  setupOp(&Array->ops[num], memLoc, op, in0, in1, zeroAddr, oneAddr);

  return;
}

int countMem(expOps *Ops, int numOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of memory locations used                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, count = 0;

  for (i = 0; i < numOps; i++)
    count += 3 + !isUnary(Ops[i].op);

  return count;
}

void derivMemLoc(expArrayOps *derivArray, expOps *operation, int *outDeriv, int deriv0, int deriv1, int zeroAddr, int oneAddr, int *first_free_mem_loc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the derivative of operation using memLocDerivs *
\***************************************************************/
{
  char op = operation->op;

  if (opUsesSubFunc(operation))
  { // error - this is only for standard operations
    printf("ERROR: This operation utilizes defined subfunctions!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  if (op == '=')
  { // deriv does not change
    *outDeriv = deriv0;
  }
  else if (op == 'N')
  { // deriv is just negative of deriv0
    if (deriv0 == -2)
    { // deriv is still 0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is -1
      addOp(derivArray, *first_free_mem_loc, 'N', oneAddr, -1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv is negative of deriv0
      addOp(derivArray, *first_free_mem_loc, 'N', deriv0, -1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == 'S')
  { // deriv is cos(in[0])*deriv0
    if (deriv0 == -2)
    { // deriv is 0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is cos(in[0])
      addOp(derivArray, *first_free_mem_loc, 'C', operation->in[0], -1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv is cos(in[0])*deriv0
      addOp(derivArray, *first_free_mem_loc, 'C', operation->in[0], -1, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '*', *first_free_mem_loc - 1, deriv0, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == 'C')
  { // deriv is -sin(in[0])*deriv0
    if (deriv0 == -2)
    { // deriv is 0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is -sin(in[0])
      addOp(derivArray, *first_free_mem_loc, 'S', operation->in[0], -1, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, 'N', *first_free_mem_loc - 1, -1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv is -sin(in[0])*deriv0
      addOp(derivArray, *first_free_mem_loc, 'S', operation->in[0], -1, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '*', *first_free_mem_loc - 1, deriv0, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, 'N', *first_free_mem_loc - 1, -1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == 'X')
  { // deriv is memLoc * deriv0, where memLoc = exp(in[0])
    if (deriv0 == -2)
    { // deriv is 0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is memLoc
      *outDeriv = operation->memLoc;
    }
    else
    { // deriv is memLoc * deriv0
      addOp(derivArray, *first_free_mem_loc, '*', operation->memLoc, deriv0, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == '+')
  { // deriv is deriv0 + deriv1
    if (deriv0 == -2)
    { // deriv is just deriv1
      *outDeriv = deriv1;
    }
    else if (deriv1 == -2)
    { // deriv is just deriv0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is 1 + deriv1
      if (deriv1 == -1)
      { // deriv is 1 + 1
        addOp(derivArray, *first_free_mem_loc, '+', oneAddr, oneAddr, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is 1 + deriv1
        addOp(derivArray, *first_free_mem_loc, '+', oneAddr, deriv1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else if (deriv1 == -1)
    { // deriv is deriv0 + 1
      addOp(derivArray, *first_free_mem_loc, '+', deriv0, oneAddr, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv0 + deriv1
      addOp(derivArray, *first_free_mem_loc, '+', deriv0, deriv1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == '-')
  { // deriv is deriv0 - deriv1
    if (deriv1 == -2)
    { // deriv is just deriv0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -2)
    { // deriv is just -deriv1
      if (deriv1 == -1)
      { // deriv is -1
        addOp(derivArray, *first_free_mem_loc, 'N', oneAddr, -1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is -deriv1
        addOp(derivArray, *first_free_mem_loc, 'N', deriv1, -1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else if (deriv0 == -1)
    { // deriv is 1 - deriv1
      if (deriv1 == -1)
      { // deriv is 1 - 1 = 0
        *outDeriv = -2;
      }
      else
      { // deriv is 1 - deriv1
        addOp(derivArray, *first_free_mem_loc, '-', oneAddr, deriv1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else if (deriv1 == -1)
    { // deriv is deriv0 - 1
      addOp(derivArray, *first_free_mem_loc, '-', deriv0, oneAddr, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv0 - deriv1
      addOp(derivArray, *first_free_mem_loc, '-', deriv0, deriv1, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else if (op == '*')
  { // deriv is deriv0 * in[1] + deriv1 * in[0]
    if (deriv0 == -2)
    { // deriv is deriv1 * in[0]
      if (deriv1 == -2)
      { // deriv is 0
        *outDeriv = -2;
      }
      else if (deriv1 == -1)
      { // deriv is in[0]
        *outDeriv = operation->in[0];
      }
      else
      { // deriv is deriv1 * in[0]
        addOp(derivArray, *first_free_mem_loc, '*', deriv1, operation->in[0], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else if (deriv0 == -1)
    { // deriv is in[1] + deriv1 * in[0]
      if (deriv1 == -2)
      { // deriv is in[1]
        *outDeriv = operation->in[1];
      }
      else if (deriv1 == -1)
      { // deriv is in[0] + in[1]
        addOp(derivArray, *first_free_mem_loc, '+', operation->in[0], operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is in[1] + deriv1 * in[0]
        addOp(derivArray, *first_free_mem_loc, '*', deriv1, operation->in[0], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '+', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
       (*first_free_mem_loc)++;
      }
    }
    else
    { // deriv is deriv0 * in[1] + deriv1 * in[0]
      if (deriv1 == -2)
      { // deriv is deriv0 * in[1]
        addOp(derivArray, *first_free_mem_loc, '*', deriv0, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else if (deriv1 == -1)
      { // deriv is deriv0 * in[1] + in[0]
        addOp(derivArray, *first_free_mem_loc, '*', deriv0, operation->in[1], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '+', *first_free_mem_loc - 1, operation->in[0], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is deriv0 * in[1] + deriv1 * in[0]
        addOp(derivArray, *first_free_mem_loc, '*', deriv0, operation->in[1], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '*', deriv1, operation->in[0], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '+', *first_free_mem_loc - 2, *first_free_mem_loc - 1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
  }
  else if (op == '/')
  { // deriv is (deriv0 - memLoc * deriv1) / in[1], where memLoc = in[0] / in[1]
    if (deriv1 == -2)
    { // deriv is just deriv0 / in[1]
      if (deriv0 == -2)
      { // deriv is still 0
        *outDeriv = deriv0;
      }
      else if (deriv0 == -1)
      { // deriv is just 1 / in[1]
        addOp(derivArray, *first_free_mem_loc, '/', oneAddr, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is deriv0 / in[1]
        addOp(derivArray, *first_free_mem_loc, '/', deriv0, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else if (deriv1 == -1)
    { // deriv is (deriv0 - memLoc) / in[1]
      if (deriv0 == -2)
      { // deriv is -memLoc / in[1]
        addOp(derivArray, *first_free_mem_loc, '/', operation->memLoc, operation->in[1], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, 'N', *first_free_mem_loc - 1, -1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else if (deriv0 == -1)
      { // deriv is (1 - memLoc) / in[1]
        addOp(derivArray, *first_free_mem_loc, '-', oneAddr, operation->memLoc, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '/', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is (deriv0 - memLoc) / in[1]
        addOp(derivArray, *first_free_mem_loc, '-', deriv0, operation->memLoc, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '/', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
    else 
    { // deriv is (deriv0 - memLoc * deriv1) / in[1]
      if (deriv0 == -2)
      { // deriv is -memLoc * deriv1 / in[1]
        addOp(derivArray, *first_free_mem_loc, '*', operation->memLoc, deriv1, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '/', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, 'N', *first_free_mem_loc - 1, -1, zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else if (deriv0 == -1)
      { // deriv is (1 - memLoc * deriv1) / in[1]
        addOp(derivArray, *first_free_mem_loc, '*', operation->memLoc, deriv1, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '-', oneAddr, *first_free_mem_loc - 1, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '/', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
      else
      { // deriv is (deriv0 - memLoc * deriv1) / in[1]
        addOp(derivArray, *first_free_mem_loc, '*', operation->memLoc, deriv1, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '-', deriv0, *first_free_mem_loc - 1, zeroAddr, oneAddr);
        (*first_free_mem_loc)++;
        addOp(derivArray, *first_free_mem_loc, '/', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
        *outDeriv = *first_free_mem_loc;
        (*first_free_mem_loc)++;
      }
    }
  } 
  else if (op == '^')
  { // verify that the deriv of exponent is -2 == 0!!
    if (deriv1 != -2)
    { // error!
      printf("ERROR: It appears that there is a nonconstant exponent!!\n");
      bexit(ERROR_INPUT_SYNTAX);
    }
    // deriv is in[1] * in[0] ^ (in[1] - 1) * deriv0
    if (deriv0 == -2)
    { // deriv is still 0
      *outDeriv = deriv0;
    }
    else if (deriv0 == -1)
    { // deriv is in[1] * in[0] ^ (in[1] - 1)
      addOp(derivArray, *first_free_mem_loc, '-', operation->in[1], oneAddr, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '^', operation->in[0], *first_free_mem_loc - 1, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '*', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
    else
    { // deriv is in[1] * in[0] ^ (in[1] - 1) * deriv0
      addOp(derivArray, *first_free_mem_loc, '-', operation->in[1], oneAddr, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '^', operation->in[0], *first_free_mem_loc - 1, zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '*', *first_free_mem_loc - 1, operation->in[1], zeroAddr, oneAddr);
      (*first_free_mem_loc)++;
      addOp(derivArray, *first_free_mem_loc, '*', *first_free_mem_loc - 1, deriv0, zeroAddr, oneAddr);
      *outDeriv = *first_free_mem_loc;
      (*first_free_mem_loc)++;
    }
  }
  else
  {
    printf("ERROR: Invalid operation ('%c')!!\n", op);
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void addOps(expArrayOps *Ops, defStatement *defSt, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: add the operations onto Ops                            *
\***************************************************************/
{
  int i, numOps = defSt->rvalExp.numOps;
  int numNewOps = numOps + 1, numOldOps = Ops->numOps; // +1 for '=' operation
  int rVal = defSt->rvalExp.finalMemLoc;

  // setup memory for the new operations
  Ops->numOps += numNewOps;
  Ops->ops = (expOps *)brealloc(Ops->ops, Ops->numOps * sizeof(expOps));

  // add on the operations
  for (i = 0; i < numOps; i++)
  { // copy over the operation
    copyOp(&Ops->ops[numOldOps + i], &defSt->rvalExp.ops[i], zeroAddr, oneAddr);
  }

  // setup '=' operation
  setupOp(&Ops->ops[numOldOps + numOps], defSt->lvalMemLoc, '=', rVal, -1, zeroAddr, oneAddr); 

  return;
}

void setupConstOps(expArrayOps *constOps, parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the contant defining operations                  *
\***************************************************************/
{
  int i;
  int zeroAddr = Array->memoryLoc.numStart, oneAddr = Array->memoryLoc.numStart + 1;

  // initialize
  initialize_expArrayOps(constOps);

  // loop over the defining statements finding the ones that define constants
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == CONSTANTTYPE)
    { // this defines a constant - add it to the list
      addOps(constOps, &Array->definingStatements[i], zeroAddr, oneAddr);
    }

  return;
}

void setupParamDiffOps(expArrayOps *paramOps, expArrayOps *paramDerivOps, parseArray *Array, int defStateNum, int *memLocDerivs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the parameter and its deriv w.r.t. pathvar       *
\***************************************************************/
{
  int i, numOps = Array->definingStatements[defStateNum].rvalExp.numOps;
  int numNewOps = numOps + 1, numOldOps = paramOps->numOps; // +1 for '=' operation
  int zeroAddr = Array->memoryLoc.numStart, oneAddr = Array->memoryLoc.numStart + 1;
  int paramNum = Array->definingStatements[defStateNum].lvalLoc;
  int paramDiffLoc = Array->memoryLoc.paramDerivStart + paramNum, paramAddr = Array->memoryLoc.paramStart + paramNum;
  int rVal = Array->definingStatements[defStateNum].rvalExp.finalMemLoc;
  int *outDeriv = NULL, deriv0, deriv1;
  expOps *tempOp = NULL;

  // setup memory for new operations in paramOps
  paramOps->numOps += numNewOps;
  paramOps->ops = (expOps *)brealloc(paramOps->ops, paramOps->numOps * sizeof(expOps));

  // loop over the operations - add on the operations and compute the deriv
  for (i = 0; i < numOps; i++)
  { // copy over the operation
    tempOp = &paramOps->ops[numOldOps + i];
    copyOp(tempOp, &Array->definingStatements[defStateNum].rvalExp.ops[i], zeroAddr, oneAddr);
    // compute the derivative of this operation
    outDeriv = &memLocDerivs[tempOp->memLoc];
    deriv0 = memLocDerivs[tempOp->in[0]];
    if (!isUnary(tempOp->op))
    { // setup deriv1
      deriv1 = memLocDerivs[tempOp->in[1]];
    }
    else
    { // setup deriv1
      deriv1 = -2;
    }
    derivMemLoc(paramDerivOps, tempOp, outDeriv, deriv0, deriv1, zeroAddr, oneAddr, &Array->firstFreeMemLoc);
  }

  // setup '=' operation
  tempOp = &paramOps->ops[numOldOps + numOps];
  setupOp(tempOp, paramAddr, '=', rVal, -1, zeroAddr, oneAddr);

  // differentiate the '=' operation
  outDeriv = &memLocDerivs[tempOp->memLoc];
  deriv0 = memLocDerivs[tempOp->in[0]];
  deriv1 = -2;
  derivMemLoc(paramDerivOps, tempOp, outDeriv, deriv0, deriv1, zeroAddr, oneAddr, &Array->firstFreeMemLoc);

  // setup the '=' operation
  addOp(paramDerivOps, paramDiffLoc, '=', memLocDerivs[paramAddr], -1, zeroAddr, oneAddr);

  return;
}

void setupParamOps(expArrayOps *paramOps, expArrayOps *paramDerivOps, parseArray *Array, int *memLocDerivs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the parameters and their derivs w.r.t. pathvars  *
\***************************************************************/
{
  int i;

  // initialize
  initialize_expArrayOps(paramOps);
  initialize_expArrayOps(paramDerivOps);

  // loop over the defining statements finding the ones that define parameters
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == PARAMETERTYPE)
    { // this defines a parameter - initialize memory, add it to the list and differentiate it
      initialize_memLocDerivs(memLocDerivs, &Array->memoryLoc);
      // set deriv of the path variable to 1
      memLocDerivs[Array->memoryLoc.pathvarStart] = -1; 

      setupParamDiffOps(paramOps, paramDerivOps, Array, i, memLocDerivs);
    }

  return;
}

void diffOp(expArrayOps *diffOps, expOps *Op, int *firstFreeMemLoc, int *memLocDerivs, int zeroAddr, int oneAddr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate the operation                            *
\***************************************************************/
{ // determine if operation uses a defined subfunction or no
  if (opUsesSubFunc(Op))
  { // setup the operations for this defined subfunction
    printf("ERROR: This function does not differentiate defined subfunctions!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  else
  { // differentiate normally
    int *outDeriv = &memLocDerivs[Op->memLoc], deriv0 = memLocDerivs[Op->in[0]], deriv1 = isUnary(Op->op) ? -1 : memLocDerivs[Op->in[1]];
    derivMemLoc(diffOps, Op, outDeriv, deriv0, deriv1, zeroAddr, oneAddr, firstFreeMemLoc);
    outDeriv = NULL;
  }

  return;
}

void copyFuncOp(expArrayOps *funcOps, expOps *Op, parseArray *Array, int zeroAddr, int oneAddr, subFuncData *defSFData)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the function operation                           *
\***************************************************************/
{ // determine if this operation uses a defined subfunction or not
  if (opUsesSubFunc(Op))
  { // setup the operations using defined subfunction
    addDefinedSubfuncOp(funcOps, Op, Array, zeroAddr, oneAddr, defSFData);
  }
  else
  { // add the operation to funcOps
    addOp(funcOps, Op->memLoc, Op->op, Op->in[0], Op->in[1], zeroAddr, oneAddr);
  }

  return;
}

void homogenizeOp(expArrayOps *outOps, expOps *inOp, parseArray *Array, variablegroupArray *vargpArray, int **degArray, char *funcName, comp_d *memArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: homogenize inOp and add to outOps                      *
\***************************************************************/
{
  int i, j, zeroAddr = Array->memoryLoc.numStart, oneAddr = Array->memoryLoc.numStart + 1;
  int numGps = vargpArray->numGps, newDeg, currVarGp = 0;
  int deg0, deg1, newLoc0, newLoc1;
  char op = inOp->op;

  if (op == '=' || op == 'N')
  { // leave the degrees alone
    for (i = 0; i < numGps; i++)
      degArray[inOp->memLoc][i] = degArray[inOp->in[0]][i];

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, inOp->in[0], inOp->in[1], zeroAddr, oneAddr);
  }
  else if (op == 'S' || op == 'C' || op == 'X')
  { // verify degrees are all zero
    for (i = 0; i < numGps; i++)
    { // input must be constant w.r.t. variables
      if (degArray[inOp->in[0]][i] != 0)
      {
        printf("ERROR: It appears that %s uses %s with an input argument that depends upon a variable.\n", funcName, op == 'S' ? "sin" : (op == 'C' ? "cos" : "exp"));
        bexit(ERROR_INPUT_SYNTAX);
      }

      // set degree to 0
      degArray[inOp->memLoc][i] = 0;
    }

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, inOp->in[0], inOp->in[1], zeroAddr, oneAddr);
  }
  else if (op == '+' || op == '-')
  { // make the degrees are equal
    newLoc0 = inOp->in[0];
    newLoc1 = inOp->in[1];
    for (i = 0; i < numGps; i++)
    { // see what type of variable group this is
      deg0 = degArray[inOp->in[0]][i];
      deg1 = degArray[inOp->in[1]][i];

      if (vargpArray->types[i])
      { // variable_group - make sure they are equal
        if (deg0 < deg1)
        { // increase the degree on the first input
          for (j = deg0; j < deg1; j++)
          { // multiply by correct homogenizing variable
            addOp(outOps, Array->firstFreeMemLoc, '*', Array->memoryLoc.varStart + currVarGp, newLoc0, zeroAddr, oneAddr);
            newLoc0 = Array->firstFreeMemLoc;
            Array->firstFreeMemLoc++;
          }

          // store the new degree
          newDeg = deg1;
        }
        else if (deg0 > deg1)
        { // increase the degree on the second input 
          for (j = deg1; j < deg0; j++)
          { // multiply by correct homogenizing variable
            addOp(outOps, Array->firstFreeMemLoc, '*', Array->memoryLoc.varStart + currVarGp, newLoc1, zeroAddr, oneAddr);
            newLoc1 = Array->firstFreeMemLoc;
            Array->firstFreeMemLoc++;
          }
  
          // store the new degree
          newDeg = deg0;
        }
        else
        { // both are the same
          newDeg = deg0;
        }

        // increment currVarGp
        currVarGp++;
      }
      else
      { // hom_variable_group - must be equal!
        if (deg0 != deg1)
        {
          printf("ERROR: It appears that %s is not properly homogenized in the variables:\n", funcName);
          for (j = 0; j < vargpArray->totalNumVars; j++)
            if (vargpArray->groupNumber[j] == i)
              printf("%s ", vargpArray->name[j]);
          printf("\n");
          bexit(ERROR_INPUT_SYNTAX);
        }
        // both are the same
        newDeg = deg0;
      }

      // store the new degree
      degArray[inOp->memLoc][i] = newDeg;
    }

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, newLoc0, newLoc1, zeroAddr, oneAddr);
  }
  else if (op == '*')
  { // add on to the degrees
    for (i = 0; i < numGps; i++)
      degArray[inOp->memLoc][i] = degArray[inOp->in[0]][i] + degArray[inOp->in[1]][i];

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, inOp->in[0], inOp->in[1], zeroAddr, oneAddr);
  }
  else if (op == '/')
  { // degree of the numerator - denom must not depend on variables
    for (i = 0; i < numGps; i++)
    { // denom must be constant w.r.t. variables
      if (degArray[inOp->in[1]][i] != 0)
      {
        printf("ERROR: It appears that %s has a nonconstant denominator!\n", funcName);
        bexit(ERROR_INPUT_SYNTAX);
      }

      // set degree to that of numerator
      degArray[inOp->memLoc][i] = degArray[inOp->in[0]][i];
    }

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, inOp->in[0], inOp->in[1], zeroAddr, oneAddr);
  }
  else if (op == '^')
  { // compute the total degree of the base and verify exponent is constant
    j = 0; 
    for (i = 0; i < numGps; i++)
    {
      j += degArray[inOp->in[0]][i];

      if (degArray[inOp->in[1]][i] != 0)
      {
        printf("ERROR: It appears that %s has a nonconstant exponent!\n", funcName);
        bexit(ERROR_INPUT_SYNTAX);
      }
    }

    if (j != 0)
    { // base is not constant - exponent must be nonnegative integer
      verifyIntegerExponent(memArray[inOp->in[1]], funcName);

      // mutlipy the degree by the exponent - a nonnegative integer
      newDeg = (int) memArray[inOp->in[1]]->r;
      for (i = 0; i < numGps; i++)
      { // set degree of the result to base * exponent
        degArray[inOp->memLoc][i] = degArray[inOp->in[0]][i] * newDeg;
      }
    }
    else
    { // result is constant
      for (i = 0; i < numGps; i++)
        degArray[inOp->memLoc][i] = 0;
    }

    // copy the operation to outOps
    addOp(outOps, inOp->memLoc, inOp->op, inOp->in[0], inOp->in[1], zeroAddr, oneAddr);
  }
  else
  {
    printf("ERROR: Invalid operation ('%c')!!\n", op);
    bexit(ERROR_INPUT_SYNTAX);
  }

  return;
}

void homogenizeOps(expArrayOps *outOps, expArrayOps *inOps, parseArray *Array, variablegroupArray *vargpArray, int **degArray, char *funcName, comp_d *memArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: homogenize inOps and add to outOps                     *
\***************************************************************/
{
  int i, numOps = inOps->numOps;

  for (i = 0; i < numOps; i++)
  { // homogenize the ith operation of inOps
    homogenizeOp(outOps, &inOps->ops[i], Array, vargpArray, degArray, funcName, memArray);
  }

  return;
}

void copyFuncOps(expArrayOps *funcOps, parseArray *Array, int defStateNum, subFuncData *defSFData, variablegroupArray *vargpArray, int userHom, int ***degArray, int *degSize, comp_d **memArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the function operations                          *
\***************************************************************/
{
  int i, numOps = Array->definingStatements[defStateNum].rvalExp.numOps;
  int zeroAddr = Array->memoryLoc.numStart, oneAddr = Array->memoryLoc.numStart + 1;
  int storeAddr = 0, storeType = Array->definingStatements[defStateNum].lvalType, storeNum = Array->definingStatements[defStateNum].lvalLoc;
  char *funcName = Array->ppArray.types[storeType].name[storeNum];
  expOps tempOp;

  // setup storeAddr
  if (storeType == FUNCTIONTYPE)
  { // function
    storeAddr = Array->memoryLoc.funcStart + storeNum;
  }
  else if (storeType == INLINESUBFUNCTIONTYPE)
  { // subfunction
    storeAddr = Array->memoryLoc.subfuncStart + storeNum;
  }
  else
  { // error
    printf("ERROR: The defining statement does not define a function or subfunction!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // determine if we need to homogenize after the operations are setup
  if (userHom)
  { // loop over the operations
    for (i = 0; i < numOps; i++)
    { // setup the operation
      copyFuncOp(funcOps, &Array->definingStatements[defStateNum].rvalExp.ops[i], Array, zeroAddr, oneAddr, defSFData);
    }

    // setup '=' operation in tempOp
    setupOp(&tempOp, storeAddr, '=', Array->definingStatements[defStateNum].rvalExp.finalMemLoc, -1, zeroAddr, oneAddr);

    // setup the '=' operation in funcOps
    copyFuncOp(funcOps, &tempOp, Array, zeroAddr, oneAddr, defSFData);
  }
  else // need to homogenize
  { // setup operations in a temporary structure and then homogenize by putting into funcOps
    expArrayOps tempOps;
    initialize_expArrayOps(&tempOps);

    // loop over the operations
    for (i = 0; i < numOps; i++)
    { // setup the operation
      copyFuncOp(&tempOps, &Array->definingStatements[defStateNum].rvalExp.ops[i], Array, zeroAddr, oneAddr, defSFData);
    }

    // setup '=' operation in tempOp
    setupOp(&tempOp, storeAddr, '=', Array->definingStatements[defStateNum].rvalExp.finalMemLoc, -1, zeroAddr, oneAddr);

    // setup the '=' operation in tempOps
    copyFuncOp(&tempOps, &tempOp, Array, zeroAddr, oneAddr, defSFData);

    // reallocate memArray & degArray - could have added memory locations for defined subfunctions
    *memArray = (comp_d *)brealloc(*memArray, Array->firstFreeMemLoc * sizeof(comp_d));    
    reallocate_degArray(degArray, degSize, Array->firstFreeMemLoc, vargpArray->numGps);

    // perform the operations
    performOps(&tempOps, *memArray);

    // homogenize the operations and setup in funcOps
    homogenizeOps(funcOps, &tempOps, Array, vargpArray, *degArray, funcName, *memArray);

    // clear tempOps
    clear_expArrayOps(&tempOps);
  }

  funcName = NULL;

  return;
}

void setupDiffOps(expArrayOps *diffOps, expArrayOps *funcOps, memoryLocations *memoryLoc, int *firstFreeMemLoc, int *memLocDerivs, int storeType, int storeNum, int isJv, int varOrParamNumber, int depend)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the derivs w.r.t. memLocDerivs                   *
\***************************************************************/
{
  int i, numOps = funcOps->numOps;
  int zeroAddr = memoryLoc->numStart, oneAddr = memoryLoc->numStart + 1;
  int numVars = memoryLoc->numVars, numParams = memoryLoc->numParams;
  int storeAddr = 0, storeLoc = 0;

  // setup storeLoc
  if (storeType == FUNCTIONTYPE)
  { // function
    storeAddr = memoryLoc->funcStart + storeNum;
    if (isJv)
      storeLoc = memoryLoc->funcDerivVStart + storeNum * numVars + varOrParamNumber;
    else
      storeLoc = memoryLoc->funcDerivPStart + storeNum * numParams + varOrParamNumber;
  }
  else if (storeType == INLINESUBFUNCTIONTYPE)
  { // subfunction
    storeAddr = memoryLoc->subfuncStart + storeNum;
    if (isJv)
      storeLoc = memoryLoc->subfuncDerivVStart + storeNum * numVars + varOrParamNumber;
    else
      storeLoc = memoryLoc->subfuncDerivPStart + storeNum * numParams + varOrParamNumber;
  }
  else
  { // error
    printf("ERROR: The defining statement does not define a function or subfunction!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  // determine if this has dependence
  if (depend)
  { // loop over the operations to compute deriv
    for (i = 0; i < numOps; i++)
    { // differentiate the operation
      diffOp(diffOps, &funcOps->ops[i], firstFreeMemLoc, memLocDerivs, zeroAddr, oneAddr);
    } 
  }
  else
  { // set deriv to 0
    memLocDerivs[storeAddr] = -2;
  }

  // store the deriv to diffOps
  addOp(diffOps, storeLoc, '=', memLocDerivs[storeAddr], -1, zeroAddr, oneAddr);

  return;
}

void setupFuncOps(expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, int **funcOrder, int **subfuncOrder, int **startSubfuncs, int **endSubfuncs, int **startFuncs, int **endFuncs, int **startJvSub, int **endJvSub, int **startJv, int **endJv, parseArray *Array, subFuncData *defSFData, variablegroupArray *vargpArray, int userHom, int ***degArray, int *degSize, comp_d **memArray, int orderFunc, int useParallelDiff, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the functions and derivs w.r.t. vars & params    *
\***************************************************************/
{
  int i, j, isJv, numVars = Array->memoryLoc.numVars, numParams = Array->memoryLoc.numParams;
  int numVarsParams = numVars + numParams;
  int numFuncs = Array->memoryLoc.numFuncs, numSubfuncs = Array->memoryLoc.numSubfuncs;
  int count = 0, currFunc = 0, currSubFunc = 0, currFuncOp = 0, currJvOp = 0, currJpOp = 0, currMemFunc = 0, currMemJv = 0;
  int totalJvOps = 0, totalJpOps = 0;
  int totalMem = 0, *memLocDerivs = NULL;
  int zeroAddr = Array->memoryLoc.numStart, oneAddr = Array->memoryLoc.numStart + 1;
  int **sfDepend = (int **)bmalloc(numSubfuncs * sizeof(int *)), **fDepend = (int **)bmalloc(numFuncs * sizeof(int *));
  int *numSubFuncOps = (int *)bmalloc(numSubfuncs * sizeof(int)), *numFuncOps = (int *)bmalloc(numFuncs * sizeof(int));
  expArrayOps *func_temp = (expArrayOps *)bmalloc(numFuncs * sizeof(expArrayOps));
  expArrayOps *subFunc_temp = (expArrayOps *)bmalloc(numSubfuncs * sizeof(expArrayOps));
  expArrayOps **Jv_func_temp = (expArrayOps **)bmalloc(numFuncs * sizeof(expArrayOps *));
  expArrayOps **Jv_subFunc_temp = (expArrayOps **)bmalloc(numSubfuncs * sizeof(expArrayOps *));
  expArrayOps **Jp_func_temp = (expArrayOps **)bmalloc(numFuncs * sizeof(expArrayOps *));
  expArrayOps **Jp_subFunc_temp = (expArrayOps **)bmalloc(numSubfuncs * sizeof(expArrayOps *));

  // initialize temp ops, sfDepend & fDepend
  for (i = 0; i < numFuncs; i++)
  { // allocate and initialize
    initialize_expArrayOps(&func_temp[i]);
    Jv_func_temp[i] = (expArrayOps *)bmalloc(numVars * sizeof(expArrayOps));
    for (j = 0; j < numVars; j++)
      initialize_expArrayOps(&Jv_func_temp[i][j]);
    Jp_func_temp[i] = (expArrayOps *)bmalloc(numParams * sizeof(expArrayOps));
    for (j = 0; j < numParams; j++)
      initialize_expArrayOps(&Jp_func_temp[i][j]);
    fDepend[i] = (int *)bmalloc(numVarsParams * sizeof(int));
    for (j = 0; j < numVarsParams; j++)
      fDepend[i][j] = 0;
  }
  for (i = 0; i < numSubfuncs; i++)
  { // allocate and initialize
    initialize_expArrayOps(&subFunc_temp[i]);
    Jv_subFunc_temp[i] = (expArrayOps *)bmalloc(numVars * sizeof(expArrayOps));
    for (j = 0; j < numVars; j++)
      initialize_expArrayOps(&Jv_subFunc_temp[i][j]);
    Jp_subFunc_temp[i] = (expArrayOps *)bmalloc(numParams * sizeof(expArrayOps));
    for (j = 0; j < numParams; j++)
      initialize_expArrayOps(&Jp_subFunc_temp[i][j]);
    sfDepend[i] = (int *)bmalloc(numVarsParams * sizeof(int));
    for (j = 0; j < numVarsParams; j++)
      sfDepend[i][j] = 0;
  }

  // initialize funcOps, JvOps & JpOps
  initialize_expArrayOps(funcOps);
  initialize_expArrayOps(JvOps);
  initialize_expArrayOps(JpOps);

  // allocate other data
  *funcOrder = (int *)bmalloc(numFuncs * sizeof(int));
  *subfuncOrder = (int *)bmalloc(numSubfuncs * sizeof(int));
  *startSubfuncs = (int *)bmalloc(numSubfuncs * sizeof(int));
  *endSubfuncs = (int *)bmalloc(numSubfuncs * sizeof(int));
  *startFuncs = (int *)bmalloc(numFuncs * sizeof(int));
  *endFuncs = (int *)bmalloc(numFuncs * sizeof(int));
  *startJvSub = (int *)bmalloc(numSubfuncs * sizeof(int));
  *endJvSub = (int *)bmalloc(numSubfuncs * sizeof(int));
  *startJv = (int *)bmalloc(numFuncs * sizeof(int));
  *endJv = (int *)bmalloc(numFuncs * sizeof(int));

  // reorder defining statements of functions, if needed
  if (orderFunc)
    reorderFunctions(Array);

  // setup funcOps
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
    { // copy this function to funcOps
      currFunc = Array->definingStatements[i].lvalLoc;
      (*funcOrder)[currFunc] = count;
      (*startFuncs)[currFunc] = currMemFunc;

      // copy to funcOps
      copyFuncOps(funcOps, Array, i, defSFData, vargpArray, userHom, degArray, degSize, memArray);

      // setup the number of operations for this function
      numFuncOps[currFunc] = funcOps->numOps - currFuncOp;

      // udpate memory count
      (*endFuncs)[currFunc] = currMemFunc += countMem(&funcOps->ops[currFuncOp], numFuncOps[currFunc]);
      currFuncOp = funcOps->numOps;

      // increment the total count of functions & subfunctions
      count++;
    }
    else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
    { // this defines a subfunction - setup the data needed
      currSubFunc = Array->definingStatements[i].lvalLoc;
      (*subfuncOrder)[currSubFunc] = count;
      (*startSubfuncs)[currSubFunc] = currMemFunc;

      // copy to funcOps
      copyFuncOps(funcOps, Array, i, defSFData, vargpArray, userHom, degArray, degSize, memArray);

      // setup the number of operations for this subfunction
      numSubFuncOps[currSubFunc] = funcOps->numOps - currFuncOp;

      // udpate memory count
      (*endSubfuncs)[currSubFunc] = currMemFunc += countMem(&funcOps->ops[currFuncOp], numSubFuncOps[currSubFunc]);
      currFuncOp = funcOps->numOps;

      // increment the total count of functions & subfunctions
      count++;
    }

  // point to the proper operations using func_temp & subFunc_temp
  currFuncOp = 0;
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
    { // function
      currFunc = Array->definingStatements[i].lvalLoc;
      // setup func_temp
      func_temp[currFunc].ops = &funcOps->ops[currFuncOp];
      // increment
      currFuncOp += func_temp[currFunc].numOps = numFuncOps[currFunc];
    }
    else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
    { // subfunction
      currSubFunc = Array->definingStatements[i].lvalLoc;
      // setup subFunc_temp
      subFunc_temp[currSubFunc].ops = &funcOps->ops[currFuncOp];
      // increment
      currFuncOp += subFunc_temp[currSubFunc].numOps = numSubFuncOps[currSubFunc];
    }

  // determine the dependencies now that the subfunction & function operations are completely setup
  determineDependencies_func(sfDepend, fDepend, Array, func_temp, subFunc_temp);

  // setup totalMem
  totalMem = Array->firstFreeMemLoc;

#ifdef _HAVE_MPI
  if (num_processes > 1 && useParallelDiff)
  { // send the totalMem
    MPI_Bcast(&totalMem, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // send defining statements type & loc
    int *defType = (int *)bmalloc(Array->numDefiningStatements * sizeof(int)), *defLoc = (int *)bmalloc(Array->numDefiningStatements * sizeof(int));
    for (i = 0; i < Array->numDefiningStatements; i++)
    {
      defType[i] = Array->definingStatements[i].lvalType;
      defLoc[i] = Array->definingStatements[i].lvalLoc;
    }

    // send the defining statements type & loc
    bcast_definingStatement_type_loc(&defType, &defLoc, &Array->numDefiningStatements, my_id, num_processes, headnode);
    free(defType);
    free(defLoc);

    // send memoryLocations
    bcast_memoryLocations(&Array->memoryLoc, my_id, num_processes, headnode);

    // send func & subFunc
    for (i = 0; i < numFuncs; i++)
      bcast_expArrayOps(&func_temp[i], my_id, num_processes, headnode);
    for (i = 0; i < numSubfuncs; i++)
      bcast_expArrayOps(&subFunc_temp[i], my_id, num_processes, headnode);

    // send depenedencies
    bcast_dependencies(fDepend, sfDepend, numFuncs, numSubfuncs, numVarsParams, my_id, num_processes, headnode);
  }
#endif

  isJv = 1;
  if (num_processes == 1 || !useParallelDiff)
  { // serial processing 

    // setup memLocDerivs now that the functions are setup
    memLocDerivs = (int *)bmalloc(totalMem * sizeof(int));

    // loop over the variables to compute the derivatives
    for (j = 0; j < numVars; j++)
    { // initialize derivs
      initialize_memLocDerivs(memLocDerivs, &Array->memoryLoc);
      // jth variable
      memLocDerivs[Array->memoryLoc.varStart + j] = -1;

      // loop over the defining statements finding the ones that define functions or inline subfunctions
      for (i = 0; i < Array->numDefiningStatements; i++)
        if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
        { // defines a function - computes its derivative w.r.t. jth variable
          currFunc = Array->definingStatements[i].lvalLoc;

          // compute deriv of func_temp[currFunc] w.r.t. jth variable
          setupDiffOps(&Jv_func_temp[currFunc][j], &func_temp[currFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, FUNCTIONTYPE, currFunc, isJv, j, fDepend[currFunc][j]);

          // increment the number of operations
          totalJvOps += Jv_func_temp[currFunc][j].numOps;
        }
        else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
        { // defines a subfunction - compute its derivatives w.r.t. jth variable
          currSubFunc = Array->definingStatements[i].lvalLoc;

          // compute deriv of subFunc_temp[currSubFunc] w.r.t. jth variable
          setupDiffOps(&Jv_subFunc_temp[currSubFunc][j], &subFunc_temp[currSubFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, INLINESUBFUNCTIONTYPE, currSubFunc, isJv, j, sfDepend[currSubFunc][j]);

          // increment the number of operations
          totalJvOps += Jv_subFunc_temp[currSubFunc][j].numOps;
        }
    }
  }
  else 
  { // parallel processing
#ifdef _HAVE_MPI
    parallel_diff_vars_head(Jv_func_temp, Jv_subFunc_temp, func_temp, subFunc_temp, Array, sfDepend, fDepend, isJv, totalMem, my_id, num_processes, headnode);

    // compute the total number of operations
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < numVars; j++)
        totalJvOps += Jv_func_temp[i][j].numOps;
    for (i = 0; i < numSubfuncs; i++)
      for (j = 0; j < numVars; j++)
        totalJvOps += Jv_subFunc_temp[i][j].numOps;
#endif
  }

  // setup Jv
  JvOps->numOps = totalJvOps;
  JvOps->ops = (expOps *)bmalloc(totalJvOps * sizeof(expOps));
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
    { // copy over the operations for this function
      currFunc = Array->definingStatements[i].lvalLoc;
      (*startJv)[currFunc] = currMemJv;

      for (j = 0; j < numVars; j++)
      { // copy Jv_func_temp[currFunc][j] to JvOps - update currJvOp & currMemJv
        currMemJv += copyOps(&JvOps->ops[currJvOp], Jv_func_temp[currFunc][j].ops, Jv_func_temp[currFunc][j].numOps, zeroAddr, oneAddr);
        currJvOp += Jv_func_temp[currFunc][j].numOps;
        // clear Jv_func_temp[currFunc][j]
        clear_expArrayOps(&Jv_func_temp[currFunc][j]);
      }
      // free Jv_func_temp[currFunc]
      free(Jv_func_temp[currFunc]);

      // store the end of this function
      (*endJv)[currFunc] = currMemJv;
    }
    else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
    { // copy over the operations for this subfunction
      currSubFunc = Array->definingStatements[i].lvalLoc;
      (*startJvSub)[currSubFunc] = currMemJv;

      for (j = 0; j < numVars; j++)
      { // copy Jv_subFunc_temp[currSubFunc][j] to JvOps - update currJvOp & currMemJv
        currMemJv += copyOps(&JvOps->ops[currJvOp], Jv_subFunc_temp[currSubFunc][j].ops, Jv_subFunc_temp[currSubFunc][j].numOps, zeroAddr, oneAddr);
        currJvOp += Jv_subFunc_temp[currSubFunc][j].numOps;
        // clear Jv_subFunc_temp[currSubFunc][j]
        clear_expArrayOps(&Jv_subFunc_temp[currSubFunc][j]);
      }
      // free Jv_subFunc_temp[currSubFunc]
      free(Jv_subFunc_temp[currSubFunc]);

      // store the end of this function
      (*endJvSub)[currSubFunc] = currMemJv; 
    }

  // clear Jv_func_temp & Jv_subFunc_temp
  free(Jv_func_temp);
  free(Jv_subFunc_temp);

  isJv = 0;
  if (num_processes == 1 || !useParallelDiff)
  { // serial processing 

    // loop over the parameters to compute the derivatives
    for (j = 0; j < numParams; j++)
    { // initialize derivs
      initialize_memLocDerivs(memLocDerivs, &Array->memoryLoc);
      // jth parameter
      memLocDerivs[Array->memoryLoc.paramStart + j] = -1;

      // loop over the defining statements finding the ones that define functions or inline subfunctions
      for (i = 0; i < Array->numDefiningStatements; i++)
        if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
        { // defines a function - computes its derivative w.r.t. jth parameter
          currFunc = Array->definingStatements[i].lvalLoc;

          // compute deriv of func_temp[currFunc] w.r.t. jth parameter
          setupDiffOps(&Jp_func_temp[currFunc][j], &func_temp[currFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, FUNCTIONTYPE, currFunc, isJv, j, fDepend[currFunc][numVars + j]);

          // increment the number of operations
          totalJpOps += Jp_func_temp[currFunc][j].numOps;
        }
        else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
        { // defines a subfunction - compute its derivatives w.r.t. jth parameter
          currSubFunc = Array->definingStatements[i].lvalLoc;

          // compute deriv of subFunc_temp[currSubFunc] w.r.t. jth parameter
          setupDiffOps(&Jp_subFunc_temp[currSubFunc][j], &subFunc_temp[currSubFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, INLINESUBFUNCTIONTYPE, currSubFunc, isJv, j, sfDepend[currSubFunc][numVars + j]);

          // increment the number of operations
          totalJpOps += Jp_subFunc_temp[currSubFunc][j].numOps;
        }
    }
  }
  else
  { // parallel processing
#ifdef _HAVE_MPI
    parallel_diff_vars_head(Jp_func_temp, Jp_subFunc_temp, func_temp, subFunc_temp, Array, sfDepend, fDepend, isJv, totalMem, my_id, num_processes, headnode);

    // compute the total number of operations
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < numParams; j++)
        totalJpOps += Jp_func_temp[i][j].numOps;
    for (i = 0; i < numSubfuncs; i++)
      for (j = 0; j < numParams; j++)
        totalJpOps += Jp_subFunc_temp[i][j].numOps;
#endif
  }

  // setup Jp
  JpOps->numOps = totalJpOps;
  JpOps->ops = (expOps *)bmalloc(totalJpOps * sizeof(expOps));
  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
    { // copy over the operations for this function
      currFunc = Array->definingStatements[i].lvalLoc;

      for (j = 0; j < numParams; j++)
      { // copy Jp_func_temp[currFunc][j] to JpOps - update currJpOp
        copyOps(&JpOps->ops[currJpOp], Jp_func_temp[currFunc][j].ops, Jp_func_temp[currFunc][j].numOps, zeroAddr, oneAddr);
        currJpOp += Jp_func_temp[currFunc][j].numOps;
        // clear Jp_func_temp[currFunc][j]
        clear_expArrayOps(&Jp_func_temp[currFunc][j]);
      }
      // free Jp_func_temp[currFunc]
      free(Jp_func_temp[currFunc]);
    }
    else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
    { // copy over the operations for this subfunction
      currSubFunc = Array->definingStatements[i].lvalLoc;

      for (j = 0; j < numParams; j++)
      { // copy Jp_func_temp[currFunc][j] to JpOps - update currJpOp
        copyOps(&JpOps->ops[currJpOp], Jp_subFunc_temp[currSubFunc][j].ops, Jp_subFunc_temp[currSubFunc][j].numOps, zeroAddr, oneAddr);
        currJpOp += Jp_subFunc_temp[currSubFunc][j].numOps;
        // clear Jp_subFunc_temp[currSubFunc][j]
        clear_expArrayOps(&Jp_subFunc_temp[currSubFunc][j]);
      }
      // free Jp_subFunc_temp[currSubFunc]
      free(Jp_subFunc_temp[currSubFunc]);
    }

  // clear Jp_func_temp & Jp_subFunc_temp
  free(Jp_func_temp);
  free(Jp_subFunc_temp);

  // clear func_temp & subFunc_temp, fDepend & sfDepend
  for (i = 0; i < numFuncs; i++)
  {
    func_temp[i].ops = NULL;
    free(fDepend[i]);
  }
  free(func_temp);
  free(fDepend);
  for (i = 0; i < numSubfuncs; i++)
  {
    subFunc_temp[i].ops = NULL;
    free(sfDepend[i]);
  }
  free(subFunc_temp);
  free(sfDepend);

  // clear other memory
  free(numSubFuncOps);
  free(numFuncOps);
  free(memLocDerivs);

  return;
}

void computeNumDenom(char **numer, char **denom, char *s)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the numerator and denominator for s          *
\***************************************************************/
{
  char *exp = NULL, ch;
  int indexIN = 0;  // Keeps track of where we are in the string (including dec point).
  int indexOUT = 0; // Keeps track of where we are in the array to be printed (excluding dec point).
  int i, ptSeen = 0, powerOfTen = 0, exp_int;  // the exponent as an integer.
  int size_exp;

  // original allocation to remove decimal point and the scientific notation - so originally setup to be at most the size of the string + 1
  *numer = (char *)bmalloc((strlen(s) + 1) * sizeof(char));

  ch = s[0];
  while ((ch != '\0') && (ch != 'e') && (ch != 'E'))
  {
    if (ch == '.')
    {
      ptSeen = 1;
      indexIN++;
    }
    else
    {
      (*numer)[indexOUT++] = ch;
      if (ptSeen == 1)
        powerOfTen++;
      indexIN++;
    }
    ch = s[indexIN];
  }
  (*numer)[indexOUT] = '\0';

  if ((ch == 'e') || (ch == 'E'))  //In other words, if we are using scientific notation.
  { // move past 'e' or 'E' and find the next char
    ch = s[++indexIN];
    // store the exponent
    size_exp = 2;
    exp = (char *)bmalloc(size_exp * sizeof(char));
    i = 0;
    while (ch != '\0')
    { // store the char
      exp[i++] = ch;
      // find the next char
      ch = s[++indexIN];
      // see if we need to increase the size of exp
      if (i + 1 >= size_exp)
      { // increase the size
        exp = (char *)brealloc(exp, 2 * size_exp * sizeof(char));
        size_exp *= 2;
      }
    }
    exp[i] = '\0';
    exp_int = atoi(exp);
    if (exp_int < 0)  //I.e., if the exponent of 10 if negative, we just increase powerOfTen (the denom) by the correct amount.
      powerOfTen -= exp_int;
    else  //If the exponent of 10 is positive:
    {
      if (powerOfTen > exp_int)  //If our denominator has more zeros than the exponent of 10, we just take some off.
        powerOfTen -= exp_int;
      else //Otherwise, the denominator should be 1 (so that powerOfTen = 0) and the numerator needs to be made larger!
      {
        exp_int -= powerOfTen;
        powerOfTen = 0;
        // increase to the correct size of noDecPt
        *numer = (char *)brealloc(*numer, (strlen(*numer) + exp_int + 1) * sizeof(char));
        for (i=0; i<exp_int; i++)  //Here we add on the appropriate number of zeros.
          (*numer)[indexOUT++] = '0';
        (*numer)[indexOUT] = '\0';
      }
    }
  }

  // allocate denom
  *denom = (char *)bmalloc((powerOfTen + 2) * sizeof(char));
  (*denom)[0] = '1';
  for (i=1; i<=powerOfTen; i++)
    (*denom)[i] = '0';
  (*denom)[powerOfTen+1] = '\0';

  if (exp != NULL)
    free(exp);

  return;
}

int printOps2(FILE *OUT, expOps *funcOps, int numOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of instructions printed                 *
* NOTES: prints the operations to OUT                           *
\***************************************************************/
{
  int i, count = 0;
  char op;

  for (i = 0; i < numOps; i++)
  {
    op = funcOps[i].op;
    if (isUnary(op))
    {
      count += 3;
      fprintf(OUT, "%d %d %d ", op, funcOps[i].memLoc, funcOps[i].in[0]);
    }
    else if (op == '+' || op == '-' || op == '*' || op == '/' || op == '^')
    {
      count += 4;
      fprintf(OUT, "%d %d %d %d ", op, funcOps[i].memLoc, funcOps[i].in[0], funcOps[i].in[1]);
    }
    else
    {
      printf("ERROR: '%c' is an invalid operation!\n", op);
      bexit(ERROR_CONFIGURATION);
    }
  }

  return count;
}

void setupArrOut(FILE *OUT, expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, int *funcOrder, int *subfuncOrder, int *startSubfuncs, int *endSubfuncs, int *startFuncs, int *endFuncs, int *startJvSub, int *endJvSub, int *startJv, int *endJv, memoryLocations *memoryLoc, int totalRoomNeeded, int numVarGps, int *varGpSizes, int randIndex)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup arr.out                                          *
\***************************************************************/
{
  int i, j, end, endUpdate, endParams, endFn, endPDerivs, endJvEval, endJpEval;
  int numSubfuncs = memoryLoc->numSubfuncs;
  int numFuncs = memoryLoc->numFuncs;

  // print the constant expressions
  end = endUpdate = printOps2(OUT, constOps->ops, constOps->numOps);

  // print the parameter expressions
  endParams = end += printOps2(OUT, paramOps->ops, paramOps->numOps);

  // print the functions
  endFn = end += printOps2(OUT, funcOps->ops, funcOps->numOps);

  // print the derivatives of the parameters w.r.t. path variables
  endPDerivs = end += printOps2(OUT, paramDerivOps->ops, paramDerivOps->numOps);

  // print the derivatives of the functions w.r.t. variables
  endJvEval = end += printOps2(OUT, JvOps->ops, JvOps->numOps);

  // print the derivatives of the functions w.r.t. parameters
  endJpEval = end += printOps2(OUT, JpOps->ops, JpOps->numOps);

  // print an 'X' at the end of the array
  fprintf(OUT, "X\n");

  fprintf(OUT, "TotalRoomNeeded %d;\n", totalRoomNeeded);
  fprintf(OUT, "NVAR %d %d;\n", memoryLoc->numVars, memoryLoc->varStart);
  fprintf(OUT, "NPATHVAR %d %d;\n", memoryLoc->numPathVars, memoryLoc->pathvarStart);
  fprintf(OUT, "NPAR %d %d %d;\n", memoryLoc->numParams, memoryLoc->paramStart, memoryLoc->paramDerivStart);
  fprintf(OUT, "NFCN %d %d %d %d;\n", memoryLoc->numFuncs, memoryLoc->funcStart, memoryLoc->funcDerivVStart, memoryLoc->funcDerivPStart);
  fprintf(OUT, "NCON %d %d;\n", memoryLoc->numConsts, memoryLoc->constStart);
  fprintf(OUT, "NNUM %d %d;\n", memoryLoc->numNums, memoryLoc->numStart);
  fprintf(OUT, "SUBFCN %d %d %d %d;\n", memoryLoc->numSubfuncs, memoryLoc->subfuncStart, memoryLoc->subfuncDerivVStart, memoryLoc->subfuncDerivPStart);
  fprintf(OUT, "NUMINST %d;\n", end);
  fprintf(OUT, "CMPLX %d;\n", memoryLoc->constStart); // first constant is always 'I'
  fprintf(OUT, "VARGPS %d;\n", numVarGps);
  for (i = 0; i < numVarGps; i++)
    fprintf(OUT, " %d", varGpSizes[i]);
  fprintf(OUT, ";\n");
  fprintf(OUT, "RANDINDEX %d;\n", randIndex);
  fprintf(OUT, "INSTCOUNT %d %d %d %d %d;\n", endUpdate, endParams, endFn, endPDerivs, endJvEval);
  // print the start and end for the subfunctions
  for (i = 0; i < numSubfuncs; i++)
    fprintf(OUT, " %d %d", endParams + startSubfuncs[i], endParams + endSubfuncs[i]);
  // print the start and end for the functions
  for (i = 0; i < numFuncs; i++)
    fprintf(OUT, " %d %d", endParams + startFuncs[i], endParams + endFuncs[i]);
  // print the start and end of the derivs for subfunctions
  for (i = 0; i < numSubfuncs; i++)
    fprintf(OUT, " %d %d", endPDerivs + startJvSub[i], endPDerivs + endJvSub[i]);
  // print the start and end of the derivs for functions
  for (i = 0; i < numFuncs; i++)
    fprintf(OUT, " %d %d", endPDerivs + startJv[i], endPDerivs + endJv[i]);
  fprintf(OUT, ";\n");
  // print the 'matrix' of the subfunctions 'below' each function
  for (i = 0; i < numFuncs; i++)
  { // print the subfunctions below the ith function
    for (j = 0; j < numSubfuncs; j++)
      fprintf(OUT, " %d", subfuncOrder[j] < funcOrder[i] ? 1 : 0);
    if (numSubfuncs > 0)
      fprintf(OUT, ";\n");
  }
  fprintf(OUT, "X diff complete\n\n");

  return;
}

void setupSLPNumOut(int numNums, char **nums, char *numName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup num.out                                          *
\***************************************************************/
{
  int i;
  char *numer = NULL, *denom = NULL;
  FILE *OUT = fopen(numName, "w");
  if (OUT == NULL)
  { // error
    printf("ERROR: '%s' is an invalid name!\n", numName);
    bexit(ERROR_CONFIGURATION);
  }

  for (i = 0; i < numNums; i++)
  {
    computeNumDenom(&numer, &denom, nums[i]);
    fprintf(OUT, "%s/%s ;\n", numer, denom);
    free(numer);
    free(denom);
    numer = denom = NULL;
  }
  fclose(OUT);

  return;
}

void setupSLPOutputFiles(variablegroupArray *vargpArray, expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, int *funcOrder, int *subfuncOrder, int *startSubfuncs, int *endSubfuncs, int *startFuncs, int *endFuncs, int *startJvSub, int *endJvSub, int *startJv, int *endJv, parseArray *Array, int numVarGps, int *varGpSizes, int randIndex, int **degArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup names.out, num.out, arr.out & deg.out            *
\***************************************************************/
{
  int i, j, numVars = Array->memoryLoc.numVars, numNums = Array->memoryLoc.numNums, numFuncs = Array->memoryLoc.numFuncs, funcStart = Array->memoryLoc.funcStart;
  FILE *OUT = NULL;

  // print the variable names to names.out
  OUT = fopen("names.out", "w");
  for (i = 0; i < numVars; i++)
    fprintf(OUT, "%s\n", vargpArray->name[i]);
  fclose(OUT);

  // print the numbers to num.out
  setupSLPNumOut(numNums, Array->ppArray.types[NUMBERTYPE].name, "num.out");

  // print the SLP to arr.out
  OUT = fopen("arr.out", "w");
  setupArrOut(OUT, constOps, paramOps, paramDerivOps, funcOps, JvOps, JpOps, funcOrder, subfuncOrder, startSubfuncs, endSubfuncs, startFuncs, endFuncs, startJvSub, endJvSub, startJv, endJv, &Array->memoryLoc, Array->firstFreeMemLoc, numVarGps, varGpSizes, randIndex);
  fclose(OUT);

  // print the degrees to deg.out
  OUT = fopen("deg.out", "w");
  for (i = 0; i < numFuncs; i++)
  { // print the degrees for the ith function
    for (j = 0; j < numVarGps; j++)
      fprintf(OUT, "%d\n", degArray[funcStart + i][j]);
    fprintf(OUT, "\n");
  }
  fclose(OUT);

  OUT = NULL;

  return;
}

void initialize_degArray(int ***degArray, int *degSize, int totalMem, variablegroupArray *vargpArray, memoryLocations *memLoc, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, numGps = vargpArray->numGps;

  // initialize degSize
  *degSize = 0;

  if (!userHom)
  { // allocate based on totalMem & numGps
    *degSize = totalMem;
    *degArray = (int **)bmalloc(totalMem * sizeof(int *));
    for (i = 0; i < totalMem; i++)
    {
      (*degArray)[i] = (int *)bmalloc(numGps * sizeof(int));
      for (j = 0; j < numGps; j++)
        (*degArray)[i][j] = 0;
    }

    // setup based for the variables
    for (i = 0; i < vargpArray->totalNumVars; i++)
    { // put a 1 in the memory location & variable group corresponding to the ith variable
      (*degArray)[memLoc->varStart + i][vargpArray->groupNumber[i]] = 1;    
    }
  }

  return;
}

void reallocate_degArray(int ***degArray, int *degSize, int newSize, int numGps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;

  *degArray = (int **)brealloc(*degArray, newSize * sizeof(int *));
  for (i = *degSize; i < newSize; i++)
  {
    (*degArray)[i] = (int *)bmalloc(numGps * sizeof(int));
    for (j = 0; j < numGps; j++)
      (*degArray)[i][j] = 0;
  }  

  *degSize = newSize;

  return;
}

void clear_degArray(int ***degArray, int degSize, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  if (!userHom)
  { // allocate based on totalMem & numGps
    for (i = 0; i < degSize; i++)
      free((*degArray)[i]);
    free(*degArray);
  }

  *degArray = NULL;

  return;
}

void init_memArray(comp_d **memArray, int totalMem, parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  *memArray = (comp_d *)bmalloc(totalMem * sizeof(comp_d));
  for (i = 0; i < totalMem; i++)
    set_zero_d((*memArray)[i]);

  // setup variables as random
  for (i = 0; i < Array->memoryLoc.numVars; i++)
    get_comp_rand_d((*memArray)[Array->memoryLoc.varStart + i]);

  // setup pathvars as random
  for (i = 0; i < Array->memoryLoc.numPathVars; i++)
    get_comp_rand_d((*memArray)[Array->memoryLoc.pathvarStart + i]);

  // setup parameters as random
  for (i = 0; i < Array->memoryLoc.numParams; i++)
    get_comp_rand_d((*memArray)[Array->memoryLoc.paramStart + i]);

  // setup numbers
  for (i = 0; i < Array->memoryLoc.numNums; i++)
    set_double_d((*memArray)[Array->memoryLoc.numStart + i], atof(Array->ppArray.types[NUMBERTYPE].name[i]), 0);

  // setup I
  set_double_d((*memArray)[Array->memoryLoc.constStart], 0, 1);

  // setup Pi
  set_double_d((*memArray)[Array->memoryLoc.constStart + 1], M_PI, 0);

  return;
}

void initialize_memArray(comp_d **memArray, int totalMem, parseArray *Array, int userHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  if (!userHom)
  {
    init_memArray(memArray, totalMem, Array);
  }
  else
  {
    *memArray = NULL;
  }

  return;
}

void clear_memArray(comp_d **memArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  if (*memArray != NULL)
  {
    free(*memArray);
    *memArray = NULL;
  }

  return;
}

void performExpOp(expOps *op, comp_d *mem)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
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
  else if (op->op == 'S')
  {
    sin_d(mem[op->memLoc], mem[op->in[0]]);
  }
  else if (op->op == 'C')
  {
    cos_d(mem[op->memLoc], mem[op->in[0]]);
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

void performOps(expArrayOps *ops, comp_d *mem)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, numOps = ops->numOps;

  for (i = 0; i < numOps; i++)
  { // perform the operation
    performExpOp(&ops->ops[i], mem);
  }
  
  return;
}

void compute_SLP(double *degBound, double *coBound, int computeBounds, parseArray *Array, subFuncData *defSFData, variablegroupArray *vargpArray, int userHom, int orderFunc, int useParallelDiff, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int totalMem = Array->firstFreeMemLoc;
  int *memLocDerivs = (int *)bmalloc(totalMem * sizeof(int));
  int *funcOrder = NULL, *subfuncOrder = NULL;
  int *startSubfuncs = NULL, *endSubfuncs = NULL, *startFuncs = NULL, *endFuncs = NULL, *startJvSub = NULL, *endJvSub = NULL, *startJv = NULL, *endJv = NULL;
  int degSize = 0, **degArray = NULL;
  comp_d *memArray = NULL;
  expArrayOps constOps, paramOps, paramDerivOps, funcOps, JvOps, JpOps;

  // setup degArray & memArray
  initialize_degArray(&degArray, &degSize, totalMem, vargpArray, &Array->memoryLoc, userHom);
  initialize_memArray(&memArray, totalMem, Array, userHom);

  // setup the constant statements
  setupConstOps(&constOps, Array);

  // setup parameter statements and derivs w.r.t. pathvariables
  setupParamOps(&paramOps, &paramDerivOps, Array, memLocDerivs);

  if (!userHom)
  { // perform constOps
    performOps(&constOps, memArray);
    // perform paramOps
    performOps(&paramOps, memArray);
  }

  // setup subfunction & function statements and derivs w.r.t. variables & parameters
  setupFuncOps(&funcOps, &JvOps, &JpOps, &funcOrder, &subfuncOrder, &startSubfuncs, &endSubfuncs, &startFuncs, &endFuncs, &startJvSub, &endJvSub, &startJv, &endJv, Array, defSFData, vargpArray, userHom, &degArray, &degSize, &memArray, orderFunc, useParallelDiff, my_id, num_processes, headnode);

  // setup degree & coeff bound, if needed
  *degBound = *coBound = 0;
  if (computeBounds)
  { // compute degree & coeff bound
    degreeBound(degBound, Array->memoryLoc.numFuncs, Array->memoryLoc.funcStart, vargpArray->numGps, degArray);
    coeffBound(coBound, &constOps, &paramOps, &paramDerivOps, &funcOps, &JvOps, &JpOps, Array);
  }
  else
  { // verify all exponents (this is also done by coeffBound)
    verify_exponents(&constOps, &paramOps, &paramDerivOps, &funcOps, &JvOps, &JpOps, Array);
  }

  // setup output files now that everything is known
  setupSLPOutputFiles(vargpArray, &constOps, &paramOps, &paramDerivOps, &funcOps, &JvOps, &JpOps, funcOrder, subfuncOrder, startSubfuncs, endSubfuncs, startFuncs, endFuncs, startJvSub, endJvSub, startJv, endJv, Array, vargpArray->numGps, vargpArray->sizes, 0, degArray);

  // clear Ops memory
  clear_expArrayOps(&constOps);
  clear_expArrayOps(&paramOps);
  clear_expArrayOps(&paramDerivOps);
  clear_expArrayOps(&funcOps);
  clear_expArrayOps(&JvOps);
  clear_expArrayOps(&JpOps);

  // clear degArray & memArray
  clear_degArray(&degArray, degSize, userHom);
  clear_memArray(&memArray);

  // clear other memory
  free(funcOrder); free(subfuncOrder);
  free(startSubfuncs); free(endSubfuncs);
  free(startFuncs); free(endFuncs);
  free(startJvSub); free(endJvSub);
  free(startJv); free(endJv);

  free(memLocDerivs);

  return;
}

void determineDependencies_mem(int *depend, int **sfDepend, parseArray *Array, int memLoc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the dependencies of memLoc                   *
\***************************************************************/
{
  int i, temp;
  int varStart = Array->memoryLoc.varStart, numVars = Array->memoryLoc.numVars;
  int varEnd = varStart + numVars;
  int paramStart = Array->memoryLoc.paramStart, numParams = Array->memoryLoc.numParams;
  int paramEnd = paramStart + numParams, numVarsParams = numVars + numParams;
  int subfuncStart = Array->memoryLoc.subfuncStart, numSubfuncs = Array->memoryLoc.numSubfuncs;
  int subfuncEnd = subfuncStart + numSubfuncs;

  if (varStart <= memLoc && memLoc < varEnd)
  { // this is a variable
    depend[memLoc - varStart] = 1;
  }
  else if (paramStart <= memLoc && memLoc < paramEnd)
  { // this is a parameter
    depend[memLoc - paramStart + numVars] = 1;
  }
  else if (subfuncStart <= memLoc && memLoc < subfuncEnd)
  { // this is a subfunction - copy over dependencies
    temp = memLoc - subfuncStart; // which subfunction
    for (i = 0; i < numVarsParams; i++)
      depend[i] = depend[i] || sfDepend[temp][i]; // or
  }

  return;
}

void determineDependencies_op(int *depend, int **sfDepend, parseArray *Array, expOps *Op)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the dependencies of Op                       *
\***************************************************************/
{
  // determine dependencies for the first input value
  determineDependencies_mem(depend, sfDepend, Array, Op->in[0]);

  if (!isUnary(Op->op))
  { // determine dependencies for the second input value
    determineDependencies_mem(depend, sfDepend, Array, Op->in[1]);
  }

  return;
}

void determineDependencies_ops(int *depend, int **sfDepend, parseArray *Array, expArrayOps *Ops)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the dependencies of Ops                      *
\***************************************************************/
{
  int i, numOps = Ops->numOps;

  // loop over the operations
  for (i = 0; i < numOps; i++)
  { // determine the dependencies for this operation
    determineDependencies_op(depend, sfDepend, Array, &Ops->ops[i]);
  }

  return;
}

void determineDependencies_func(int **sfDepend, int **fDepend, parseArray *Array, expArrayOps *func, expArrayOps *subFunc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the dependencies                             *
\***************************************************************/
{
  int i, currFunc, currSubFunc;

  for (i = 0; i < Array->numDefiningStatements; i++)
    if (Array->definingStatements[i].lvalType == FUNCTIONTYPE)
    { // function
      currFunc = Array->definingStatements[i].lvalLoc;
      // setup the dependencies for this function
      determineDependencies_ops(fDepend[currFunc], sfDepend, Array, &func[currFunc]);
    }
    else if (Array->definingStatements[i].lvalType == INLINESUBFUNCTIONTYPE)
    { // subfunction
      currSubFunc = Array->definingStatements[i].lvalLoc;
      // setup the dependencies for this subfunction
      determineDependencies_ops(sfDepend[currSubFunc], sfDepend, Array, &subFunc[currSubFunc]);
    }

  return;
}

void degreeBound(double *degBound, int numFuncs, int funcStart, int numVarGps, int **degArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the degree & coeff bound                     *
\***************************************************************/
{
  int i, j, currDeg;

  // initialize degBound
  *degBound = 1;

  // compute degree bound
  for (i = 0; i < numFuncs; i++)
  { // compute the total degree of the ith function
    currDeg = 0;
    for (j = 0; j < numVarGps; j++)
      currDeg += degArray[funcStart + i][j];
    if (currDeg > *degBound)
    *degBound = currDeg;
  }

  return;
}

void coeffBound(double *cBound, expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximate coeff bound                                *
\***************************************************************/
{
  int i, j, tempInt, numFuncs = Array->memoryLoc.numFuncs, numVars = Array->memoryLoc.numVars, numParams = Array->memoryLoc.numParams;
  int funcStart = Array->memoryLoc.funcStart, JvStart = Array->memoryLoc.funcDerivVStart, JpStart = Array->memoryLoc.funcDerivPStart;
  double tempD;
  comp_d *memArray = NULL;

  // initialize cBound
  *cBound = 1;

  // compute coeff bound
  init_memArray(&memArray, Array->firstFreeMemLoc, Array);
  performOps(constOps, memArray);
  performOps(paramOps, memArray);
  performOps(funcOps, memArray);
  performOps(paramDerivOps, memArray);
  performOps(JvOps, memArray);
  performOps(JpOps, memArray);

  // find maximum of functions & derivatives of functions
  for (i = 0; i < numFuncs; i++)
  {
    tempD = d_abs_d(memArray[funcStart + i]);
    if (tempD > *cBound)
      *cBound = tempD;

    for (j = 0; j < numVars; j++)
    {
      tempInt = JvStart + i*numVars + j;
      tempD = d_abs_d(memArray[tempInt]);
      if (tempD > *cBound)
        *cBound = tempD;
    }

    for (j = 0; j < numParams; j++)
    {
      tempInt = JpStart + i*numParams + j;
      tempD = d_abs_d(memArray[tempInt]);
      if (tempD > *cBound)
        *cBound = tempD;
    }
  }

  // clear memArray
  clear_memArray(&memArray);

  return;
}

void verify_exponents(expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: verify exponents -> real numbers                       *
\***************************************************************/
{
  comp_d *memArray = NULL;

  // initialize memArray
  init_memArray(&memArray, Array->firstFreeMemLoc, Array);
  performOps(constOps, memArray);
  performOps(paramOps, memArray);
  performOps(funcOps, memArray);
  performOps(paramDerivOps, memArray);
  performOps(JvOps, memArray);
  performOps(JpOps, memArray);

  // clear memArray
  clear_memArray(&memArray);

  return;
}


