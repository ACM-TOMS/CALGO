// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _LOCALDIM_H
#define _LOCALDIM_H

#include "bertini.h"
#include "diff.h"

// SLP for evaluating all partial derivs up to a given order
typedef struct {
  int *prog;          // program instructions
  int  size;          // size of prog
  int  memSizeNeeded; // amount of memory needed for workspace

  _comp_d  *mem_d;  // double precision workspace array 
  _comp_mp *mem_mp; // higher precision workspace array
  int  memSize_d;   // size of allocated array in double precision
  int  memSize_mp;  // size of allocated array in higher precision

  num_t *nums;    // array of real numbers
  int  precision; // current precision to utlize

  int numLinears; // number of extra linear equations (patches or slices)
  int linearType; // -1 - not initialized, 0 - double, 1 - MP, 2 - AMP
  mat_d linearCoeff_d; // coefficients for the linear equations
  mat_mp linearCoeff_mp;
  mpq_t ***linearCoeff_rat;
  vec_d linearConst_d; // value of the constant (patches == 1, slices == 0)
  vec_mp linearConst_mp;
  mpq_t **linearConst_rat;

  int  order; // top order of the partial derivatives that are calculated (order >= 1)

  int  numInstAtEndUpdate; // instruction number at end of update

  int  numVars;     // number of variables
  int  numNums;     // number of real numbers
  int  numConsts;   // number of constants 
  int  numFuncs;    // number of functions
  int  numSubfuncs; // number of subfunctions
  int *numDerivs;   // number of derivs of each order
  int *numSubDerivs;// number of subderivs of each order

  int  inpVars;  // where input variables are stored
  int  IAddr;    // where 'I' is stored
  int  numAddr;  // whre first of nums is stored
  int  constAddr;// where first constant is stored

  int  evalFuncs;  // where the first function is stored
  int  evalSubs;   // where the first subfunction is stored 
  int *evalJFuncs; // where the first derivative of the first function of a given order is stored
  int *evalJSubs;  // where the first derivative of the first subfunction of a given order is stored

  systemStruct sys; // used to hold the function operations
  int totalOpCount; // total number of operations
  int *memLoc;      // memory locations
  int **derivAddr;  // derivatives
  funcStruct **subFuncDerivs; // hold the derivatives of the subfunctions [order][numSubDerivs]
  funcStruct **funcDerivs; // hold the derivatives of the functions [order][numDerivs]

} prog_deriv_t;

// function declarations for ldt_eval.c - evaluate derivatives
int evalDeriv_d(point_d funcVals, point_d derivVals, point_d linVals, point_d linDerivVals, point_d vars, prog_deriv_t *Prog);
int evalDeriv_mp(point_mp funcVals, point_mp derivVals, point_mp linVals, point_mp linDerivVals, point_mp vars, prog_deriv_t *Prog);
void setup_deriv_from_SLP(prog_deriv_t *deriv, prog_t *SLP);
void clear_deriv(prog_deriv_t *deriv);
void setupNext_derivs(prog_deriv_t *deriv);
void add_linears_to_deriv(prog_deriv_t *deriv, mat_d linears_d, vec_d const_d, mat_mp linears_mp, vec_mp const_mp, mpq_t ***linears_rat, mpq_t **const_rat, int linearsMP);
void add_slices_patch_to_deriv(prog_deriv_t *deriv, mat_d slices_d, mat_mp slices_mp, mpq_t ***slices_rat, vec_d patch_d, vec_mp patch_mp, mpq_t **patch_rat, int MP);
void add_patches_to_deriv(prog_deriv_t *deriv, mat_d patch_d, mat_mp patch_mp, mpq_t ***patch_rat, int patchMP);
int change_prec_prog_deriv(void const *ED, int prec);

// function declarations for ldt_mult_mat.c - setup the multiplicity matrix
void setup_multiplicity_matrix_d(mat_d MM, int order, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int numFuncs, int numLinears);
void combine_func_deriv_vals_d(vec_d vals, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int max_diff_order, int *diff_size);
void setup_MM_entry_d(mat_d MM, int row_num, int col_num, vec_d vals, int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size);
void setup_multiplicity_matrix_mp(mat_mp MM, int order, vec_mp funcVals, vec_mp derivVals, vec_mp linVals, vec_mp linDerivVals, int numVars, int numFuncs, int numLinears);
void combine_func_deriv_vals_mp(vec_mp vals, vec_mp funcVals, vec_mp derivVals, vec_mp linVals, vec_mp linDerivVals, int numVars, int max_diff_order, int *diff_size);
void setup_MM_entry_mp(mat_mp MM, int row_num, int col_num, vec_mp vals, int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size);

void start_array(int *array, int length, int weight);
void start_array_rev(int *array, int length, int weight);
int advance_array(int *array, int length);
int advance_array_rev(int *array, int length);
int array_compare(int *array1, int *array2, int N);
int find_partial_index(int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size);

// function declarations for ldt_main.c - main LDT function
int is_isolated(int *mult, prog_deriv_t *deriv, point_d pt1_d, point_mp pt1_mp, int pt1_prec, point_d pt2_d, point_mp pt2_mp, int pt2_prec, int mult_bound, tracker_config_t *T, int printHilbert);
int is_isolated_d(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_d pt1_d, point_d pt2_d, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert);
int is_isolated_mp(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_mp pt1_mp, point_mp pt2_mp, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert);
int is_isolated_amp(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_d pt1_d, point_mp pt1_mp, int pt1_prec, point_d pt2_d, point_mp pt2_mp, int pt2_prec, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert);

// function declarations for ldt_rank.c - rank determination functions
int rank_MM_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d MM1, mat_d MM2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs, double max_CN, double max_SV_ratio, double SV_tol);
int rank_MM_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp MM1, mat_mp MM2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_mp **houseHolderVecs, double max_CN, double max_SV_ratio, double SV_tol);
int rank_MM_amp(double *CN, double *smallest_nonzero, double *largest_zero, mat_d MM1_d, mat_mp MM1_mp, int prec1, mat_d MM2_d, mat_mp MM2_mp, int prec2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs_d, vec_mp **houseHolderVecs_mp, double max_CN, double max_SV_ratio);

// function declarations for ldt_dual_basis.c - compute dual basis
void create_dual_basis(vec_d **dual_d, vec_mp **dual_mp, int *dual_prec, int multiplicity, mat_d MM_d, mat_mp MM_mp, int MM_prec, int MPType, mat_d topMat_d, mat_mp topMat_mp, vec_d **randVec_d, vec_mp **randVec_mp, int setupTopMat);
void find_hilbert_func_d(int **hilFn, int *hilOrder, int *reg, int mult, vec_d *dual1_d, vec_d *dual2_d, int *colSize, int colOrder, tracker_config_t *T);
void find_hilbert_func_mp(int **hilFn, int *hilOrder, int *reg, int mult, vec_mp *dual1_mp, vec_mp *dual2_mp, int *colSize, int colOrder, tracker_config_t *T);

// ldt_diff.c
void diff_funcStruct_old(funcStruct *func, int endFuncCount, int currVar, int storeAddr, int zeroAddr, int oneAddr, int *firstFreeMemLoc, int *memLoc, int **derivAddr, int totalOpCount, int *derivCount, func_ops **deriv_ops);
void diff_sys_vars_subfuncs_old(funcStruct *diff, systemStruct *sys, int currSubFunc, int *memLoc, int **derivAddr, int totalOpCount);
void diff_sys_vars_funcs_old(funcStruct *diff, systemStruct *sys, int currFunc, int *memLoc, int **derivAddr, int totalOpCount);

// testing
void first_order_dual_space_d(int *num_first, vec_d **first_order_duals, int ***first_dual_support, int ***first_monomial_support, prog_deriv_t *deriv, point_d pt1_d, point_d pt2_d, tracker_config_t *T);
double compute_num_monomials(int num_vars, int deg);
void next_order_dual_space_d(int *num_next, vec_d **next_order_duals, int ***next_dual_support, int ***next_monomial_support, int **next_num_monomials, int next_order, int num_prev, vec_d *prev_order_duals, int **prev_dual_support, int **prev_monomial_support, prog_deriv_t *deriv, point_d pt1_d, point_d pt2_d, tracker_config_t *T);
void setup_multiplicity_matrix_mon_d(mat_d MM, int order, int **monomial_support, int *diff_size, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int numFuncs, int numLinears);

#endif

