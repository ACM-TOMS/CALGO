// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

typedef struct
{
  int  memLoc;        // memory location
  char op;            // operation
  int  in[2];         // memory locations for the two input arguments
  int  lastUsed;      // array index of where it was last used
  int  numDiffInst;   // the number of instructions needed to create derivative
  char chDiffInst[3]; // the operations needed to create derivative
  int  diffInst[3][3];// the memory locations for the operations to create derivative
} func_ops;

typedef struct
{
  int memLoc;       // memory location
  int *derivAddr;   // only used with subfunctions - where derivative is stored w.r.t. vars[i]
  int useCount;     // number of times that it is used
  int varUsedSize;  // memory size of varUsed
  int *varUsed;     // array of array indices of where it is used in the function
} var_ops;

typedef struct
{
  int num_ops;    // number of operations needed to compute the function
  func_ops *ops;
} funcStruct;

typedef struct
{
  int numVars;          // number of variables
  int varsAddr;         // address of the first variable
  
  int numPathVars;      // number of path variables
  int pathVarsAddr;     // address of the first path variable
    
  int numUpdate;        // number of update operations
  func_ops *updateOps;  // the update operations

  int numParams;        // number of parameters
  int paramAddr;        // address of the first parameter
  funcStruct *params;   // the parameter operations 
  
  int numSubfuncs;      // number of subfunctions 
  int subFuncAddr;      // address of the first subfunction
  funcStruct *subFuncs; // the subfunctions
  
  int numFuncs;         // number of functions
  int funcAddr;         // address of the first function
  funcStruct *funcs;    // the functions

  int *subFuncOrder;    // orderings of the subfunctions
  int *funcOrder;       // orderings of the functions
  
  // other important information
  int firstFreeMemLoc;  // address of the first free memory location
  int numConstants;     // number of constants
  int constAddr;        // address of the first constant
  int numNumbers;       // number of numbers
  int numAddr;          // address of the first number
  
  // other marginally useful information
  int numVarGps;        // number of variable groups
  int *varGpSizes;      // sizes of the variable groups
  int randIndex;        // rand index
} systemStruct;

// diff_deflatable.c

void checkExp(systemStruct *sys, num_t *nums, int allowNegExp);
void diff_sys(systemStruct *sys, char *paramName, char *jvName, char *jpName); // differentiate sys and print derivative instructions to the files
void diff_forward(systemStruct *sys, int *memLoc, int **derivAddr, int totalOpCount, FILE *PDERIVS, FILE *JV, FILE *JP); // differentiates sys using memLoc & derivAddr

void diff_op(func_ops *ops, int *memLoc, int *derivAddr, int totalOpCount, int oneAddr, int *first_free_mem_loc, int *new_ops, func_ops **diff_ops); // finds the derivative of the operation w.r.t. diff_var
int copyInstructions(FILE *OUT, FILE *IN, char endCh); // copy IN to OUT until reach endCh
void setupSystemStructure(systemStruct *sys, char *fileName, int putIntoOrder); // reads in the syteme from 'fileName' and sets it up in sys
void readInOps(systemStruct *sys, FILE *FUNC, int putIntoOrder); // reads in the operations
void setupParams(systemStruct *sys, func_ops *paramOps, int numParamOps); // sets up the parameters
void setupFuncs_subFuncs(systemStruct *sys, func_ops *funcOps, int numFuncOps, int putIntoOrder); // sets up the functions and subfunctions
void setupOps(func_ops **ops, int *numOps, FILE *FUNC, char stopCh); // setup the operations
int  printOps(FILE *OUT, func_ops *funcOps, int numOps); // print the operations
void clearFuncStruct(funcStruct *funcs, int numFuncs); // clear funcs
void clearVarOps(var_ops *vars, int numVars); // clear vars
void clearSystemStructure(systemStruct *sys); // clear sys

void diff_sys_pathVar(funcStruct *diff, systemStruct *sys, int currParam, int currPathVar, int *memLoc, int *derivAddr, int totalOpCount);
void diff_sys_var_subfunc(funcStruct *diff, systemStruct *sys, int currSubFunc, int currVar, int *memLoc, int *derivAddr, int totalOpCount);
void diff_sys_var_func(funcStruct *diff, systemStruct *sys, int currFunc, int currVar, int *memLoc, int *derivAddr, int totalOpCount);
void diff_sys_param_subfunc(funcStruct *diff, systemStruct *sys, int currSubFunc, int currParam, int *memLoc, int *derivAddr, int totalOpCount);
void diff_sys_param_func(funcStruct *diff, systemStruct *sys, int currFunc, int currParam, int *memLoc, int *derivAddr, int totalOpCount);

void initialize_memLoc_derivAddr(systemStruct *sys, int *memLoc, int *derivAddr, int totalOpCount);
void diff_file(char *fileName, int putIntoOrder, int allowNegExp);
void diff_sys_array(systemStruct *sys, char *fileName);

int checkExp_ops(func_ops *ops, int num_ops, comp_d *mem, int oneAddr, int allowNegOps);
void copyToOp(func_ops *fop, int memLoc, char op, int in0, int in1);
int checkExp_op(func_ops *op, comp_d *mem, int oneAddr, int allowNegExp);
void performOp(func_ops *op, comp_d *mem);

void diff_funcStruct(funcStruct *func, int endFuncCount, int storeAddr, int zeroAddr, int oneAddr, int *firstFreeMemLoc, int *memLoc, int *derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops); // differentiate a funcStruct
void setup_derivAddr_subfunc_vals(systemStruct *sys, int *derivAddr, int currVar, int **subFuncVals);

int diff_funcStruct_test(funcStruct *func, int endFuncCount, int storeAddr, int zeroAddr, int oneAddr, int firstFreeTemp, int *memLoc, int *derivAddr, int totalOpCount);

// regen_diff.c
void diff_sys_rev_arr(systemStruct *sys, char *fileName); // differentiate sys and print SLP array 'fileName' - use reverseAD
void findPartialDerivs_subfuncs(funcStruct *diff, systemStruct *sys, int currSubFunc, int startAddr);
void findPartialDerivs_funcs(funcStruct *diff, systemStruct *sys, int currFunc, int startAddr);

// diff_Jon.c
void traceFuncDown(func_ops *function, int size_function, var_ops *vars, int num_vars, int zeroAddr, int oneAddr);
int dFunc_dVars(func_ops *function, int size_function, var_ops *vars, int num_vars, int vars_index, FILE *OUT, int zeroAddr, int oneAddr, int numVars, int *first_free_mem_loc, int derivLoc);

