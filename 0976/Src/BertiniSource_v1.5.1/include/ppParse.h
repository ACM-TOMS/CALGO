// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// end of line
#define pEOL 0
#define ppEOL 0

#define DEFINEDSUBFUNCTIONTYPE 0
#define FUNCTIONTYPE 1
#define VARIABLETYPE 2
#define VARIABLEGROUPTYPE 3
#define HOMVARIABLEGROUPTYPE 4
#define PARAMETERTYPE 5
#define PATHVARIABLETYPE 6
#define CONSTANTTYPE 7
#define INLINESUBFUNCTIONTYPE 8
#define NUMBERTYPE 9

#define TAG_NUM_FUNC       1 // tag for sending the function number
#define TAG_RETURN_NUM_OPS 2 // tag for sending the number of operations
#define TAG_RETURN_OPS     3 // tag for sending the operations        
#define TAG_VAR_NUM        4 // tag for sending variable number     

typedef struct
{
  int numType;     // number of elements of this type
  int allocSize;   // number of allocated elements
  int *isDefined;  // if the element has been defined already
  int *lineNumber; // line number where it is defined
  int *nameSize;   // size of each name
  char **name;     // name of the elements
} preProcessArrayType;

typedef struct
{
  int numTypes;               // number of types
  preProcessArrayType *types; // the data about each type
} preProcessArray;

typedef struct
{
  int memLoc; // memory location of answer
  char op;   // operation to perform
  int in[2];  // input memory locations
} expOps;

typedef struct
{
  int numOps;  // number of operations
  expOps *ops; // list of operations
  int finalMemLoc; // final memory location of the expression
} expArrayOps;

typedef struct
{
  int lvalType; // type of object being defined
  int lvalLoc;  // location of object being defined
  int lvalMemLoc; // memory location of the object being defined
  expArrayOps rvalExp; // expression on the right-hand side
} defStatement;

typedef struct
{
  int numVars;
  int varStart;

  int numPathVars;
  int pathvarStart;

  int numParams;
  int paramStart;
  int paramDerivStart;

  int numFuncs;
  int funcStart;
  int funcDerivVStart;
  int funcDerivPStart;

  int numConsts;
  int constStart;

  int numNums;
  int numStart;

  int numSubfuncs;
  int subfuncStart;
  int subfuncDerivVStart;
  int subfuncDerivPStart;
} memoryLocations; // information regarding memory locations

typedef struct
{
  int argMemLoc;      // memory location of argument
  int argDerivExists; // determine if the derivative w.r.t. this argument has been already computed
  int argDerivLoc;    // location of the derivative w.r.t. this argument
} argumentData;

typedef struct
{
  int subfuncLoc;          // location of defined subfunction
  int numArguments;        // number of arguments
  argumentData *arguments; // list of argument information
} definedSubfuncData;

typedef struct
{
  int numDefiningStatements; // number of defining statements
  defStatement *definingStatements; // defining statements
  int numSubfuncCalls; // number of defined subfunction calls used
  definedSubfuncData *subfuncCalls; // defined subfunction call data
  preProcessArray ppArray; // preprocess array
  memoryLocations memoryLoc; // memory locations
  int firstFreeMemLoc; // first free memory location
} parseArray;

typedef struct
{
  int numActiveSubfuncCalls; // number of active subfunction calls
  definedSubfuncData *subfuncCalls; // defined subfunction call data
} definedSubfuncStack; 

typedef struct
{
  int totalMem;   // total number of memory locations used (variables + constants + temporary locations)
  int numVars;    // number of variables
  int numConsts;  // number of constants
  int numNumbers; // number of numbers
  int IAddr;      // address of 'I'
  int *numAddr;   // address of the numbers used
  expArrayOps funcData; // function data
  expArrayOps *derivData; // derivative data
} subFuncData;

typedef struct
{
  int numGps;    // total number of groups
  int numHomGps; // number of hom_var groups
  int numVarGps; // number of var groups
  int *types;    // type of variable group: 0 - hom_var, 1 - var
  int *sizes;    // size of variable groups

  int *varLoc;   // change 'loc' of VARIABLEGROUPTYPE to variable number
  int *homLoc;   // change 'loc' of HOMVARIABLEGROUPTYPE to variable number

  int totalNumVars; // total number of variables
  int *groupNumber; // group number for each variable
  int *nameSize;    // size of each name
  char **name;      // name of the elements
} variablegroupArray;

// function declarations in ppParse.l
void append_str(char **str1, int *lenStr1, char *str2);
void print_preproc_data(FILE *OUT, int numFuncs, int numVarGps, int numHomGps, int *types, int *sizes);
int lookup_preProcessArray_type(int *loc, preProcessArrayType *Array, char *name);
int lookup_preProcessArray(int *type, int *loc, preProcessArray *ppArray, char *name);
void addEntry_preProcessArray_type(preProcessArrayType *Array, char *name, int isDefined, int lineNo);
void initialize_preProcessArray_type(preProcessArrayType *Array);
void initialize_preProcessArray(preProcessArray *ppArray);
void clear_preProcessArray_type(preProcessArrayType *Array);
void clear_preProcessArray(preProcessArray *ppArray);
void copy_preProcessArray_type(preProcessArrayType *Out, preProcessArrayType *In);
void copy_preProcessArray(preProcessArray *Out, preProcessArray *In);
void verify_standard_homotopy(preProcessArray *Array, variablegroupArray *vargpArray, int onlyOneGp);
void verify_parameter_homotopy(preProcessArray *Array);
void verify_userdefined_homotopy(preProcessArray *Array, int userHom);
void verify_preProcessArray_type(preProcessArrayType *Array);
void verify_preProcessArray(preProcessArray *Array);
int preProcessParse(preProcessArray *preProcArray, FILE *fp, int paramHom, int addAllNums);
void setup_variablegroup_numbers(preProcessArray *Array, int *numHomGps, int **sizeHomGps, int *numVarGps, int **sizeVarGps);
void setup_variablegroupArray(preProcessArray *Array, variablegroupArray *vargpArray, int userHom);
void clear_variablegroupArray(variablegroupArray *vargpArray);
void errorStatement();

// function declarations in ppParse.y
void verify_dependent_name(int defineType, int currType, char *name, int line);

// function declarations in pParse.l
void increment_subfunc_parseArray(parseArray *Array);
void increment_definedSubfuncStack(definedSubfuncStack *sfData);
void remove_definedSubfuncDataStack(definedSubfuncStack *sfData);
void copy_definedSubfuncData(definedSubfuncData *Out, definedSubfuncData *In);
void add_arg_definedSubfuncData(definedSubfuncData *sf, int argMemLoc);
void add_op_expArrayOps(expArrayOps *exp, int memLoc, char op, int in0, int in1);
void negate_lastEntry(preProcessArrayType *Numbers);
int lookup_number_type(int *loc, preProcessArrayType *Array, char *name);
int lookup_memLoc(int type, int loc, memoryLocations *memLoc, variablegroupArray *varGp);
void initialize_defStatement(defStatement *dState, int type, int loc, memoryLocations *memLoc, variablegroupArray *varGp);
void copy_defStatement(defStatement *Out, defStatement *In);
void copy_expArrayOps(expArrayOps *Out, expArrayOps *In);
void initialize_expArrayOps(expArrayOps *exp);
void clear_expArrayOps(expArrayOps *exp);
void clear_parseArray(parseArray *parArray);
void clear_subFuncData(int numType, subFuncData **sfData);
void clear_defStatement(defStatement *dState);
int processSystem(parseArray *parArray, subFuncData **sfData, preProcessArray *preProcArray, variablegroupArray *vargpArray, int userHom, int paramHom, FILE *fp);
extern int pParseparse();
extern int pParselex();
extern int pParsewrap();
extern int pParseerror(char *s);
void readInOp(expOps *op, FILE *IN);
void setup_subFuncData_numbers(preProcessArray *Array);
void errorStatement_pParse();

// function declarations in setupSLP.c
void verifyRealExponent(comp_d exp);
void verifyBaseExponent(comp_d base, comp_d exp);
int opUsesSubFunc(expOps *Op);
void setupOp(expOps *Op, int memLoc, char op, int in0, int in1, int zeroAddr, int oneAddr);
void copyOp(expOps *Out, expOps *In, int zeroAddr, int oneAddr);
int copyOps(expOps *Out, expOps *In, int numOps, int zeroAddr, int oneAddr);
void addOp(expArrayOps *Array, int memLoc, char op, int in0, int in1, int zeroAddr, int oneAddr);
void compute_SLP(double *degBound, double *coBound, int computeBounds, parseArray *Array, subFuncData *defSFData, variablegroupArray *vargpArray, int userHom, int orderFunc, int useParallelDiff, int my_id, int num_processes, int headnode);
void setupSLPNumOut(int numNums, char **nums, char *numName);
void definedSubfuncDiffOp(expArrayOps *diffOps, expOps *Op, parseArray *Array, int *memLocDerivs, int *depend, int zeroAddr, int oneAddr, subFuncData *defSFData);
void addDefinedSubfuncOp(expArrayOps *funcOps, expOps *Op, parseArray *Array, int zeroAddr, int oneAddr, subFuncData *defSFData);
void derivMemLoc(expArrayOps *derivArray, expOps *operation, int *outDeriv, int deriv0, int deriv1, int zeroAddr, int oneAddr, int *first_free_mem_loc);
int isDefinedSubfunc(int *sfCall, int memoryLoc);
void reorderFunctions(parseArray *Array);
int countMem(expOps *Ops, int numOps);
void initialize_memLocDerivs(int *memLocDerivs, memoryLocations *memoryLoc);
void setupDiffOps(expArrayOps *diffOps, expArrayOps *funcOps, memoryLocations *memoryLoc, int *firstFreeMemLoc, int *memLocDerivs, int storeType, int storeNum, int isJv, int varOrParamNumber, int depend);
void setupArrOut(FILE *OUT, expArrayOps *constOps, expArrayOps *paramOps, expArrayOps *paramDerivOps, expArrayOps *funcOps, expArrayOps *JvOps, expArrayOps *JpOps, int *funcOrder, int *subfuncOrder, int *startSubfuncs, int *endSubfuncs, int *startFuncs, int *endFuncs, int *startJvSub, int *endJvSub, int *startJv, int *endJv, memoryLocations *memoryLoc, int totalRoomNeeded, int numVarGps, int *varGpSizes, int randIndex);

// function declarations in parallel_diff.c
void parallel_diff_vars_head(expArrayOps **J_func, expArrayOps **J_subFunc, expArrayOps *func, expArrayOps *subFunc, parseArray *Array, int **sfDepend, int **fDepend, int isJv, int totalMem, int my_id, int num_processes, int headnode);
void parallel_diff_worker(int my_id, int num_processes, int headnode);
void bcast_dependencies(int **fDepend, int **sfDepend, int numFuncs, int numSubfuncs, int numVarsParams, int my_id, int num_processes, int headnode);
void bcast_expArrayOps(expArrayOps *Array, int my_id, int num_processes, int headnode);
void bcast_definingStatement_type_loc(int **type, int **loc, int *numDefStatement, int my_id, int num_processes, int headnode);
void bcast_memoryLocations(memoryLocations *memLoc, int my_id, int num_processes, int headnode);

// function declarations in paramHomotopy.c
void setupParameterHomotopy(FILE *OUT, FILE* IN, int paramHom, int readInOld, int maxPrec, preProcessArray *ppArray, variablegroupArray *vargpArray);

// extern definitions required for bison/flex to play nicely with the autotools system.
extern int pParse_bisonerror(char *s);
extern int pParse_bisonlex();

extern int ppParse_bisonerror(char *s);
extern int ppParse_bisonlex();

