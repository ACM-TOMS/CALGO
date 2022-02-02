// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "ppParse.h"
#include "parallel.h"

// Differentiate a system in parallel - use a variable based distribution

#define FIRST_FREE_START -1000000000

#ifdef _HAVE_MPI

void postProcessOps_array(expArrayOps *func, parseArray *Array, int firstOffset);
int recvDerivs_expArrayOps(expArrayOps **J_func, expArrayOps **J_subFunc, int numFuncs, int numSubfuncs);
void sendDerivs_expArrayOps(expArrayOps *J_func, expArrayOps *J_subFunc, int numFuncs, int numSubfuncs, int varNum, int proc);
void send_recv_expArrayOps(expArrayOps *Array, int proc, int isSending);
void create_mpi_expOps(MPI_Datatype *mpi_expOps);
void parallel_diff_vars_worker(parseArray *Array, expArrayOps *func, expArrayOps *subFunc, int **fDepend, int **sfDepend, int *memLocDerivs, int *defType, int *defLoc, int numDefStatement, int my_id, int num_processes, int headnode);
void readFuncOps_depend(char *funcOps, char *funcDepend, expArrayOps *func, expArrayOps *subFunc, int **fDepend, int **sfDepend, int numFuncs, int numSubfuncs, int numVarsParams);

void parallel_diff_vars_head(expArrayOps **J_func, expArrayOps **J_subFunc, expArrayOps *func, expArrayOps *subFunc, parseArray *Array, int **sfDepend, int **fDepend, int isJv, int totalMem, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate Array in parallel                        *
\***************************************************************/
{
  int i, j, k, firstOffset, recvProc, count = 0, numDefStatement = Array->numDefiningStatements;
  int numVars = Array->memoryLoc.numVars, numParams = Array->memoryLoc.numParams, numFuncs = Array->memoryLoc.numFuncs, numSubfuncs = Array->memoryLoc.numSubfuncs;
  int numCount, currFunc, currSubFunc;
  int *lastSize = (int *)bmalloc(num_processes * sizeof(int));
  
  for (i = 0; i < num_processes; i++)
    lastSize[i] = 0;

  // setup numCount
  numCount = isJv ? numVars : numParams;

  // send isJv
  MPI_Bcast(&isJv, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // send out the initial variables
  for (i = 0; i < num_processes; i++)
    if (i != my_id && count < numCount)
    { // send out the variable
      MPI_Send(&count, 1, MPI_INT, i, TAG_NUM_FUNC, MPI_COMM_WORLD);

      // update the counts
      lastSize[i] = 1;
      count += lastSize[i];
    }
    else if (i != my_id)
    { // tell the worker not to work
      k = -1;
      MPI_Send(&k, 1, MPI_INT, i, TAG_NUM_FUNC, MPI_COMM_WORLD);

      // update the counts
      lastSize[i] = 0;
      count += lastSize[i];
    }

  // loop until all sent out
  while (count < numCount)
  { // recv one deriv back
    recvProc = recvDerivs_expArrayOps(J_func, J_subFunc, numFuncs, numSubfuncs);

    // send out the next variable
    MPI_Send(&count, 1, MPI_INT, recvProc, TAG_NUM_FUNC, MPI_COMM_WORLD);

    // update the counts
    lastSize[recvProc] = 1;
    count += lastSize[recvProc];
  }

  // count the number still out
  count = 0;
  for (i = 0; i < num_processes; i++)
    if (i != headnode && lastSize[i] > 0)
      count++;

  // loop to recv all the packets back
  while (count > 0)
  { // recv deriv back
    recvProc = recvDerivs_expArrayOps(J_func, J_subFunc, numFuncs, numSubfuncs);

    // tell the worker that this level is complete
    k = -1;
    MPI_Send(&k, 1, MPI_INT, recvProc, TAG_NUM_FUNC, MPI_COMM_WORLD);

    // decrement count
    count--;
  }

  // post-process the operations 
  for (i = 0; i < numCount; i++)
  { // setup the offset
    firstOffset = Array->firstFreeMemLoc;

    // set the numbers correctly
    for (j = 0; j < numDefStatement; j++)
      if (Array->definingStatements[j].lvalType == FUNCTIONTYPE)
      { // defines a function
        currFunc = Array->definingStatements[j].lvalLoc;

        postProcessOps_array(&J_func[currFunc][i], Array, firstOffset); 
      }
      else if (Array->definingStatements[j].lvalType == INLINESUBFUNCTIONTYPE)
      { // defines an inline subfunction
        currSubFunc = Array->definingStatements[j].lvalLoc;

        postProcessOps_array(&J_subFunc[currSubFunc][i], Array, firstOffset); 
      }
  }

  // clear memory
  free(lastSize);  

  return;
}

void postProcessOps_array(expArrayOps *func, parseArray *Array, int firstOffset)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: post process the operations                            *
\***************************************************************/
{
  int i, num_ops = func->numOps;

  for (i = 0; i < num_ops; i++)
  {
    if (func->ops[i].memLoc < 0)
    { // update memLoc
      func->ops[i].memLoc = firstOffset + func->ops[i].memLoc - FIRST_FREE_START;
      // see if we need to update firstFreeMemLoc
      if (func->ops[i].memLoc >= Array->firstFreeMemLoc)
        Array->firstFreeMemLoc = func->ops[i].memLoc + 1;
    }

    if (func->ops[i].in[0] < 0)
    { // update the first input 
      func->ops[i].in[0] = firstOffset + func->ops[i].in[0] - FIRST_FREE_START;
      // see if we need to update firstFreeMemLoc (should never happen here!)
      if (func->ops[i].in[0] >= Array->firstFreeMemLoc)
        Array->firstFreeMemLoc = func->ops[i].in[0] + 1;
    }

    if (!isUnary(func->ops[i].op))
    { // binary operation
      if (func->ops[i].in[1] < 0)
      { // update the second input
        func->ops[i].in[1] = firstOffset + func->ops[i].in[1] - FIRST_FREE_START;
        // see if we need to update firstFreeMemLoc (should never happen here!)
        if (func->ops[i].in[1] >= Array->firstFreeMemLoc)
          Array->firstFreeMemLoc = func->ops[i].in[1] + 1;
      }
    }
  }

  return;
}

// worker functions

void parallel_diff_worker(int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate in parallel                              *
\***************************************************************/
{
  int i, needToDiff, totalMem, numFuncs, numSubfuncs, numVarsParams, **fDepend = NULL, **sfDepend = NULL;
  int numDefStatement, *defType = NULL, *defLoc = NULL, *memLocDerivs = NULL;
  parseArray Array;
  expArrayOps *func = NULL, *subFunc = NULL;

  // recv needToDiff
  MPI_Bcast(&needToDiff, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  if (needToDiff)
  { // recv the totalMem
    MPI_Bcast(&totalMem, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // setup memLocDerivs
    memLocDerivs = (int *)bmalloc(totalMem * sizeof(int));

    // recv the defining statements type & loc
    bcast_definingStatement_type_loc(&defType, &defLoc, &numDefStatement, my_id, num_processes, headnode);

    // send memoryLocations
    bcast_memoryLocations(&Array.memoryLoc, my_id, num_processes, headnode);
    numFuncs = Array.memoryLoc.numFuncs;
    numSubfuncs = Array.memoryLoc.numSubfuncs;
    numVarsParams = Array.memoryLoc.numVars + Array.memoryLoc.numParams;

    // allocate func, subFunc, fDepend & sfDepend
    func = (expArrayOps *)bmalloc(numFuncs * sizeof(expArrayOps));
    subFunc = (expArrayOps *)bmalloc(numSubfuncs * sizeof(expArrayOps));
    fDepend = (int **)bmalloc(numFuncs * sizeof(int *));
    sfDepend = (int **)bmalloc(numSubfuncs * sizeof(int *));
    for (i = 0; i < numFuncs; i++)
    {
      initialize_expArrayOps(&func[i]);
      fDepend[i] = (int *)bmalloc(numVarsParams * sizeof(int));
    }
    for (i = 0; i < numSubfuncs; i++)
    {
      initialize_expArrayOps(&subFunc[i]);
      sfDepend[i] = (int *)bmalloc(numVarsParams * sizeof(int));
    }

    /// recv func & subFunc
    for (i = 0; i < numFuncs; i++)
      bcast_expArrayOps(&func[i], my_id, num_processes, headnode);
    for (i = 0; i < numSubfuncs; i++)
      bcast_expArrayOps(&subFunc[i], my_id, num_processes, headnode);
    
    // recv dependencies
    bcast_dependencies(fDepend, sfDepend, numFuncs, numSubfuncs, numVarsParams, my_id, num_processes, headnode);

    // compute derivs w.r.t. vars
    parallel_diff_vars_worker(&Array, func, subFunc, fDepend, sfDepend, memLocDerivs, defType, defLoc, numDefStatement, my_id, num_processes, headnode);
    // compute derivs w.r.t. params
    parallel_diff_vars_worker(&Array, func, subFunc, fDepend, sfDepend, memLocDerivs, defType, defLoc, numDefStatement, my_id, num_processes, headnode);

    // clear memory
    for (i = 0; i < numFuncs; i++)
    {
      clear_expArrayOps(&func[i]);
      free(fDepend[i]);
    }
    for (i = 0; i < numSubfuncs; i++)
    {
      clear_expArrayOps(&subFunc[i]);
      free(sfDepend[i]);
    }
    free(func);
    free(subFunc);
    free(fDepend);
    free(sfDepend);
    free(memLocDerivs);
    free(defType);
    free(defLoc);
  }

  // no need to clear Array

  return;
}

void parallel_diff_vars_worker(parseArray *Array, expArrayOps *func, expArrayOps *subFunc, int **fDepend, int **sfDepend, int *memLocDerivs, int *defType, int *defLoc, int numDefStatement, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: differentiate in parallel                              *
\***************************************************************/
{
  int i, varNum, isJv, currFunc, currSubFunc, numVars, numParams, numVarsParams, numFuncs, numSubfuncs;
  MPI_Status status;
  expArrayOps *J_func = NULL, *J_subFunc = NULL;

  // recv isJv
  MPI_Bcast(&isJv, 1, MPI_INT, headnode, MPI_COMM_WORLD);

  // allocate Jv_func, Jv_subFunc
  numVars = Array->memoryLoc.numVars;
  numParams = Array->memoryLoc.numParams;
  numVarsParams = numVars + numParams;
  numFuncs = Array->memoryLoc.numFuncs;
  numSubfuncs = Array->memoryLoc.numSubfuncs;
  J_func = (expArrayOps *)bmalloc(numFuncs * sizeof(expArrayOps));
  J_subFunc = (expArrayOps *)bmalloc(numSubfuncs * sizeof(expArrayOps));

  // recv a variable
  MPI_Recv(&varNum, 1, MPI_INT, headnode, TAG_NUM_FUNC, MPI_COMM_WORLD, &status);

  // loop over the variables
  while (varNum >= 0)
  { // setup derivs
    for (i = 0; i < numFuncs; i++)
      initialize_expArrayOps(&J_func[i]);
    for (i = 0; i < numSubfuncs; i++)
      initialize_expArrayOps(&J_subFunc[i]);

    // initialize the data
    initialize_memLocDerivs(memLocDerivs, &Array->memoryLoc);
    Array->firstFreeMemLoc = FIRST_FREE_START;

    // setup variable/parameter
    if (isJv)
      memLocDerivs[Array->memoryLoc.varStart + varNum] = -1;
    else
      memLocDerivs[Array->memoryLoc.paramStart + varNum] = -1;

    // loop over the defining statements computing the derivatives of the functions & subfunctions
    for (i = 0; i < numDefStatement; i++)
      if (defType[i] == FUNCTIONTYPE)
      { // defines a function
        currFunc = defLoc[i];

        // compute its derivative w.r.t. varNum
        setupDiffOps(&J_func[currFunc], &func[currFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, FUNCTIONTYPE, currFunc, isJv, varNum, fDepend[currFunc][varNum]);
      }
      else if (defType[i] == INLINESUBFUNCTIONTYPE)
      { // defines an inline subfunction
        currSubFunc = defLoc[i];

        // compute its derivative w.r.t. varNum
        setupDiffOps(&J_subFunc[currSubFunc], &subFunc[currSubFunc], &Array->memoryLoc, &Array->firstFreeMemLoc, memLocDerivs, INLINESUBFUNCTIONTYPE, currSubFunc, isJv, varNum, sfDepend[currSubFunc][varNum]);
      }

    // send derivs back
    sendDerivs_expArrayOps(J_func, J_subFunc, numFuncs, numSubfuncs, varNum, headnode);

    // clear derivs
    for (i = 0; i < numFuncs; i++)
      clear_expArrayOps(&J_func[i]);
    for (i = 0; i < numSubfuncs; i++)
      clear_expArrayOps(&J_subFunc[i]);

    // recv the next variable
    MPI_Recv(&varNum, 1, MPI_INT, headnode, TAG_NUM_FUNC, MPI_COMM_WORLD, &status);
  }

  // clear memory 
  free(J_func);
  free(J_subFunc);

  return;
}

// other functions

int recvDerivs_expArrayOps(expArrayOps **J_func, expArrayOps **J_subFunc, int numFuncs, int numSubfuncs)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: recv J_func & J_subFunc                                *
\***************************************************************/
{
  int i, recvProc, varNum;
  MPI_Status status;

  // recv the varNum
  MPI_Recv(&varNum, 1, MPI_INT, MPI_ANY_SOURCE, TAG_VAR_NUM, MPI_COMM_WORLD, &status);
  // save the processor number
  recvProc = status.MPI_SOURCE;

  // recv J_func & J_subFunc
  for (i = 0; i < numFuncs; i++)
    send_recv_expArrayOps(&J_func[i][varNum], recvProc, 0);
  for (i = 0; i < numSubfuncs; i++)
    send_recv_expArrayOps(&J_subFunc[i][varNum], recvProc, 0);

  return recvProc;
}

void sendDerivs_expArrayOps(expArrayOps *J_func, expArrayOps *J_subFunc, int numFuncs, int numSubfuncs, int varNum, int proc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: send J_func & J_subFunc                                *
\***************************************************************/
{
  int i;

  // send varNum
  MPI_Send(&varNum, 1, MPI_INT, proc, TAG_VAR_NUM, MPI_COMM_WORLD);

  // send J_func & J_subFunc
  for (i = 0; i < numFuncs; i++)
    send_recv_expArrayOps(&J_func[i], proc, 1);
  for (i = 0; i < numSubfuncs; i++)
    send_recv_expArrayOps(&J_subFunc[i], proc, 1);

  return;
}

void bcast_definingStatement_type_loc(int **type, int **loc, int *numDefStatement, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcast type & loc for defining statements           *
\***************************************************************/
{
  if (my_id == headnode)
  { // send the number of defining statements
    MPI_Bcast(numDefStatement, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // send type & loc
    MPI_Bcast(*type, *numDefStatement, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(*loc, *numDefStatement, MPI_INT, headnode, MPI_COMM_WORLD);
  }
  else
  { // recv the number of defining statements
    MPI_Bcast(numDefStatement, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // allocate memory
    *type = (int *)bmalloc(*numDefStatement * sizeof(int));
    *loc = (int *)bmalloc(*numDefStatement * sizeof(int));
    // recv type & loc
    MPI_Bcast(*type, *numDefStatement, MPI_INT, headnode, MPI_COMM_WORLD);
    MPI_Bcast(*loc, *numDefStatement, MPI_INT, headnode, MPI_COMM_WORLD);
  }

  return;
}

void bcast_memoryLocations(memoryLocations *memLoc, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcast memoryLocations                              *
\***************************************************************/
{
  int *memInt = (int *)bmalloc(19 * sizeof(int));

  if (my_id == headnode)
  { // copy data to memInt
    memInt[0] = memLoc->numVars;
    memInt[1] = memLoc->varStart;

    memInt[2] = memLoc->numPathVars;
    memInt[3] = memLoc->pathvarStart;

    memInt[4] = memLoc->numParams;
    memInt[5] = memLoc->paramStart;
    memInt[6] = memLoc->paramDerivStart;

    memInt[7] = memLoc->numFuncs;
    memInt[8] = memLoc->funcStart;
    memInt[9] = memLoc->funcDerivVStart;
    memInt[10] = memLoc->funcDerivPStart;

    memInt[11] = memLoc->numConsts;
    memInt[12] = memLoc->constStart;

    memInt[13] = memLoc->numNums;
    memInt[14] = memLoc->numStart;

    memInt[15] = memLoc->numSubfuncs;
    memInt[16] = memLoc->subfuncStart;
    memInt[17] = memLoc->subfuncDerivVStart;
    memInt[18] = memLoc->subfuncDerivPStart;

    // send it out
    MPI_Bcast(memInt, 19, MPI_INT, headnode, MPI_COMM_WORLD);
  }
  else
  { // recv data
    MPI_Bcast(memInt, 19, MPI_INT, headnode, MPI_COMM_WORLD);

    // setup data
    memLoc->numVars = memInt[0];
    memLoc->varStart = memInt[1];

    memLoc->numPathVars = memInt[2];
    memLoc->pathvarStart = memInt[3];

    memLoc->numParams = memInt[4];
    memLoc->paramStart = memInt[5];
    memLoc->paramDerivStart = memInt[6];

    memLoc->numFuncs = memInt[7];
    memLoc->funcStart = memInt[8];
    memLoc->funcDerivVStart = memInt[9];
    memLoc->funcDerivPStart = memInt[10];

    memLoc->numConsts = memInt[11];
    memLoc->constStart = memInt[12];

    memLoc->numNums = memInt[13];
    memLoc->numStart = memInt[14];

    memLoc->numSubfuncs = memInt[15];
    memLoc->subfuncStart = memInt[16];
    memLoc->subfuncDerivVStart = memInt[17];
    memLoc->subfuncDerivPStart = memInt[18];
  }

  // free memInt
  free(memInt);

  return;
}

void send_recv_expArrayOps(expArrayOps *Array, int proc, int isSending)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: send expArrayOps                                       *
\***************************************************************/
{
  MPI_Datatype mpi_expOps;
  create_mpi_expOps(&mpi_expOps);

  if (isSending)
  { // send the number of operations
    MPI_Send(&Array->numOps, 1, MPI_INT, proc, TAG_RETURN_NUM_OPS, MPI_COMM_WORLD);
    // send the operations
    MPI_Send(Array->ops, Array->numOps, mpi_expOps, proc, TAG_RETURN_OPS, MPI_COMM_WORLD);
    // send the finalMemLoc
    MPI_Send(&Array->finalMemLoc, 1, MPI_INT, proc, TAG_RETURN_NUM_OPS, MPI_COMM_WORLD);
  }
  else
  { // recv the number of operations
    MPI_Status status;
    MPI_Recv(&Array->numOps, 1, MPI_INT, proc, TAG_RETURN_NUM_OPS, MPI_COMM_WORLD, &status);
    // allocate the operations
    Array->ops = (expOps *)bmalloc(Array->numOps * sizeof(expOps));
    // recv the operations
    MPI_Recv(Array->ops, Array->numOps, mpi_expOps, proc, TAG_RETURN_OPS, MPI_COMM_WORLD, &status);
    // recv the finalMemLoc
    MPI_Recv(&Array->finalMemLoc, 1, MPI_INT, proc, TAG_RETURN_NUM_OPS, MPI_COMM_WORLD, &status);
  }

  MPI_Type_free(&mpi_expOps);

  return;
}

void bcast_expArrayOps(expArrayOps *Array, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: broadcast expArrayOps                                  *
\***************************************************************/
{
  MPI_Datatype mpi_expOps;
  create_mpi_expOps(&mpi_expOps);

  if (my_id == headnode)
  { // send the number of operations
    MPI_Bcast(&Array->numOps, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // send the operations
    MPI_Bcast(Array->ops, Array->numOps, mpi_expOps, headnode, MPI_COMM_WORLD);
    // send the finalMemLoc
    MPI_Bcast(&Array->finalMemLoc, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  }
  else
  { // recv the number of operations
    MPI_Bcast(&Array->numOps, 1, MPI_INT, headnode, MPI_COMM_WORLD);
    // allocate the operations
    Array->ops = (expOps *)bmalloc(Array->numOps * sizeof(expOps));
    // recv the operations
    MPI_Bcast(Array->ops, Array->numOps, mpi_expOps, headnode, MPI_COMM_WORLD);
    // recv the finalMemLoc
    MPI_Bcast(&Array->finalMemLoc, 1, MPI_INT, headnode, MPI_COMM_WORLD);
  }

  MPI_Type_free(&mpi_expOps);

  return;
}

void create_mpi_expOps(MPI_Datatype *mpi_expOps)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: create the expOps data structure                       *
\***************************************************************/
{
  // used to calculate displacements for datatypes
  MPI_Aint startaddress, address;

  expOps tempStruct;

  // arrays for length, displacement and datatypes in mpi_endgame_data_t
  int systemStruct_length[3] = {1, 1, 2};
  MPI_Datatype systemStruct_datatypes[3] = {MPI_INT, MPI_CHAR, MPI_INT};
  MPI_Aint systemStruct_displacements[3];

  // calculate displacements
  systemStruct_displacements[0] = 0;
  MPI_Address(&tempStruct.memLoc, &startaddress);
  MPI_Address(&tempStruct.op, &address);
  systemStruct_displacements[1] = address - startaddress;
  MPI_Address(&tempStruct.in, &address);
  systemStruct_displacements[2] = address - startaddress;

  // build the mpi datatype mpi_expOps
  MPI_Type_struct(3, systemStruct_length, systemStruct_displacements, systemStruct_datatypes, mpi_expOps);
  MPI_Type_commit(mpi_expOps);

  return;
}

void bcast_dependencies(int **fDepend, int **sfDepend, int numFuncs, int numSubfuncs, int numVarsParams, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: send dependencies                                      *
\***************************************************************/
{
  int i, j, count = 0, *depend = (int *)bmalloc((numFuncs + numSubfuncs) * numVarsParams * sizeof(int));

  if (my_id == headnode)
  { // setup depend and send it
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < numVarsParams; j++)
        depend[count++] = fDepend[i][j];
    for (i = 0; i < numSubfuncs; i++)
      for (j = 0; j < numVarsParams; j++)
        depend[count++] = sfDepend[i][j];

    MPI_Bcast(depend, (numFuncs + numSubfuncs) * numVarsParams, MPI_INT, headnode, MPI_COMM_WORLD);
  }
  else
  { // recv depend and setup fDepend & sfDepend
    MPI_Bcast(depend, (numFuncs + numSubfuncs) * numVarsParams, MPI_INT, headnode, MPI_COMM_WORLD);
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < numVarsParams; j++)
        fDepend[i][j] = depend[count++];
    for (i = 0; i < numSubfuncs; i++)
      for (j = 0; j < numVarsParams; j++)
        sfDepend[i][j] = depend[count++];
  }

  // free memory
  free(depend);

  return;
}

void printFuncOp(FILE *FUNC, expOps *op)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print op data                                          *
\***************************************************************/
{
  fprintf(FUNC, "%d %c %d %d\n", op->memLoc, op->op, op->in[0], op->in[1]);

  return;
}

void scanFuncOp(FILE *FUNC, expOps *op)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: scan op data                                           *
\***************************************************************/
{
  fscanf(FUNC, "%d %c %d %d\n", &op->memLoc, &op->op, &op->in[0], &op->in[1]);

  return;
}

void printFuncOps_depend(char *funcOps, char *funcDepend, expArrayOps *func, expArrayOps *subFunc, int **fDepend, int **sfDepend, parseArray *Array)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: create funcOps & funcDepend files                      *
\***************************************************************/
{
  int i, j, numOps;
  int numVars = Array->memoryLoc.numVars, numParams = Array->memoryLoc.numParams, numFuncs = Array->memoryLoc.numFuncs, numSubfuncs = Array->memoryLoc.numSubfuncs;
  int numVarsParams = numVars + numParams;
  FILE *FUNC = NULL;
  
  // create func ops
  FUNC = fopen(funcOps, "w");
  fprintf(FUNC, "%d %d\n", numFuncs, numSubfuncs);
  for (i = 0; i < numFuncs; i++)
    fprintf(FUNC, "%d\n", func[i].numOps);
  for (i = 0; i < numSubfuncs; i++)
    fprintf(FUNC, "%d\n", subFunc[i].numOps);
  for (i = 0; i < numFuncs; i++)
  {
    numOps = func[i].numOps;
    for (j = 0; j < numOps; j++)
      printFuncOp(FUNC, &func[i].ops[j]);
  }
  for (i = 0; i < numSubfuncs; i++)
  {
    numOps = subFunc[i].numOps;
    for (j = 0; j < numOps; j++)
      printFuncOp(FUNC, &subFunc[i].ops[j]);
  }
  fclose(FUNC);

  // create func depend
  FUNC = fopen(funcDepend, "w");
  fprintf(FUNC, "%d %d %d\n", numFuncs, numSubfuncs, numVarsParams);
  for (i = 0; i < numFuncs; i++)
  {
    for (j = 0; j < numVarsParams; j++)
      fprintf(FUNC, "%d ", fDepend[i][j]);
    fprintf(FUNC, "\n");
  }
  for (i = 0; i < numSubfuncs; i++)
  {
    for (j = 0; j < numVarsParams; j++)
      fprintf(FUNC, "%d ", sfDepend[i][j]);
    fprintf(FUNC, "\n");
  }
  fclose(FUNC);

  return;
}

void readFuncOps_depend(char *funcOps, char *funcDepend, expArrayOps *func, expArrayOps *subFunc, int **fDepend, int **sfDepend, int numFuncs, int numSubfuncs, int numVarsParams)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: read in funcOps & funcDepend files                     *
\***************************************************************/
{
  int i, j, numOps, tempFuncs, tempSubfuncs, tempVarsParams;
  FILE *FUNC = NULL;

  // read in func ops
  FUNC = fopen(funcOps, "r");
  fscanf(FUNC, "%d %d\n", &tempFuncs, &tempSubfuncs);
  // verify data
  if (numFuncs != tempFuncs || numSubfuncs != tempSubfuncs)
  {
    printf("ERROR: The number of functions is incorrect!!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  for (i = 0; i < numFuncs; i++)
    fscanf(FUNC, "%d\n", &func[i].numOps);
  for (i = 0; i < numSubfuncs; i++)
    fscanf(FUNC, "%d\n", &subFunc[i].numOps);
  for (i = 0; i < numFuncs; i++)
  {
    numOps = func[i].numOps;
    func[i].ops = (expOps *)bmalloc(numOps * sizeof(expOps));
    for (j = 0; j < numOps; j++)
      scanFuncOp(FUNC, &func[i].ops[j]);
  }
  for (i = 0; i < numSubfuncs; i++)
  {
    numOps = subFunc[i].numOps;
    subFunc[i].ops = (expOps *)bmalloc(numOps * sizeof(expOps));
    for (j = 0; j < numOps; j++)
      scanFuncOp(FUNC, &subFunc[i].ops[j]);
  }
  fclose(FUNC);

  // read in func depend
  FUNC = fopen(funcDepend, "r");
  fscanf(FUNC, "%d %d %d\n", &tempFuncs, &tempSubfuncs, &tempVarsParams);
  // verify data
  if (numFuncs != tempFuncs || numSubfuncs != tempSubfuncs)
  {
    printf("ERROR: The number of functions is incorrect!!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  // verify data
  if (tempVarsParams != numVarsParams)
  {
    printf("ERROR: The number of variables is incorrect!!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }
  for (i = 0; i < numFuncs; i++)
  {
    for (j = 0; j < numVarsParams; j++)
      fscanf(FUNC, "%d ", &fDepend[i][j]);
    fscanf(FUNC, "\n");
  }
  for (i = 0; i < numSubfuncs; i++)
  {
    for (j = 0; j < numVarsParams; j++)
      fscanf(FUNC, "%d ", &sfDepend[i][j]);
    fscanf(FUNC, "\n");
  }
  fclose(FUNC);

  return;
}





#endif

