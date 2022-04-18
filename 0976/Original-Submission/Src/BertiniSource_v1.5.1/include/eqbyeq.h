// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _EQBYEQ_H
#define _EQBYEQ_H

#include "bertini.h"

// Endpoint classification codes - only need to know where to move this on to next level or not
#define UNCLASSIFIED 0
#define MOVE_TO_NEXT 10        // this endpoint will be moved on
#define DO_NOT_MOVE_TO_NEXT 15 // this endpoint will NOT be moved on
// < 0 -> bad path (code = retVal returned from path)

/*****  Data structures *****/
typedef struct
{
  int startFunction;        // The starting function for witness data.
  int depth;                // The number of functions whose witness data will be found starting with 'startFunction'
  int num_paths;            // The number of paths stored.
  int num_sing;             // The number of singular endpoints.
  int num_nonsing;          // The number of nonsingular endpoints.
  int num_inf;              // The number of infinite endpoints.
  int num_higher_dim;       // The number of endpoints that lie on higher dimensional components.
  int num_bad;              // The number of other bad paths.

  // start system information
  int num_linears;          // The total number of linears used in the start system
  comp_d **startSystemCoeff;// The coefficients for the start system - [linear][variable]

  // information for intrinsic slicing - witness data always found using intrinsic slicing
  mat_d B;                  // The matrix for the null space
  vec_d p;                  // One vector that satisfies the original linears

  point_d *startPts;        // The startpoints of each path.
  point_d *endPts_in;       // The endpoints - intrinsic.
  point_d *endPts;          // The endpoints - extrinsic

  comp_d  *finalTs;         // The final T values (non-0 in case of failure).
  double  *condition_nums;  // The final observed condition number of each path, for classifying endpoint types.
  int *endPt_retVals;       // The retVal from endgame for each path.
  int *endPt_types;         // The type of each endpoint - must be == MOVE_TO_NEXT to move to the next level
  int *higherDim;           // Determine if point lies on higher dim'l component
} eqWitnessData_d;          // Called witnessData[#] inside EqD.

typedef struct
{
  int startFunction;         // The starting function for witness data.
  int depth;                 // The number of functions whose witness data will be found starting with 'startFunction'
  int num_paths;             // The number of paths stored.
  int num_sing;              // The number of singular endpoints.
  int num_nonsing;           // The number of nonsingular endpoints.
  int num_inf;               // The number of infinite endpoints.
  int num_higher_dim;        // The number of endpoints that lie on higher dimensional components.
  int num_bad;               // The number of other bad paths.

  // start system information
  int num_linears;           // The total number of linears used in the start system
  comp_mp **startSystemCoeff;// The coefficients for the start system - [linear][variable]
  mpq_t ***startSystemCoeff_rat;

  // information for intrinsic slicing - witness data always found using intrinsic slicing
  mat_mp B;                  // The matrix for the null space
  mpq_t ***B_rat;           
  vec_mp p;                  // One vector that satisfies the original linears
  mpq_t **p_rat;

  point_mp *startPts;        // The startpoints of each path.
  point_mp *endPts_in;       // The endpoints - intrinisc
  point_mp *endPts;          // The endpoints - extrinsic

  comp_mp  *finalTs;         // The final T values (non-0 in case of failure).
  double  *condition_nums;   // The final observed condition number of each path, for classifying endpoint types.
  int *endPt_retVals;        // The retVal from endgame for each path.
  int *endPt_types;          // The type of each endpoint - must be == MOVE_TO_NEXT to move to the next level
  int *higherDim;            // Determine if point lies on higher dim'l component
} eqWitnessData_mp;          // Called witnessData[#] inside EqD.

typedef struct
{
  int depth_x;             // the number of functions, starting at 0, that the x-variables satisfy
  int depth_y;             // the number of functions, starting at depth_x, that the y-variables satisfy
  int num_paths;           // The number of paths stored.
  int num_sing;            // The number of singular endpoints.
  int num_nonsing;         // The number of nonsingular endpoints.
  int num_inf;             // The number of infinite endpoints.
  int num_higher_dim;      // The number of endpoints that lie on higher dimensional components.
  int num_bad;             // The number of other bad paths.

  // information for intrinsic slicing
  int useIntrinsicSlice;   // determine if using intrinsic slicing
  mat_d B;                 // The matrix for the null space for the 'fixed linears'
  vec_d p;                 // One vector that satisfies the original 'fixed linears'
  vec_d p1;                // One vector that satisfies the linears at t = 1
  vec_d p0;                // One vector that satisfies the linears at t = 0
  mat_d B1;                // The matrix for the null space for the linears at t = 1
  mat_d B0;                // The matrix for the null space for the linears at t = 0

  point_d *startPts;       // The startpoints of each path.
  point_d *endPts_in;      // The endpoints - intrinsic.
  point_d *endPts;         // The endpoints - extrinsic

  comp_d  *finalTs;        // The final T values (non-0 in case of failure).
  double  *condition_nums; // The final observed condition number of each path, for classifying endpoint types.
  int *endPt_retVals;      // The retVal from endgame for each path.
  int *endPt_types;        // The type of each endpoint - must be == MOVE_TO_NEXT to move to the next level
  int *higherDim;          // Determine if point lies on higher dim'l component
} eqStageData_d;

typedef struct
{
  int depth_x;             // the number of functions, starting at 0, that the x-variables satisfy
  int depth_y;             // the number of functions, starting at depth_x, that the y-variables satisfy
  int num_paths;           // The number of paths stored.
  int num_sing;            // The number of singular endpoints.
  int num_nonsing;         // The number of nonsingular endpoints.
  int num_inf;             // The number of infinite endpoints.
  int num_higher_dim;      // The number of endpoints that lie on higher dimensional components.
  int num_bad;             // The number of other bad paths.

  // information for intrinsic slicing
  int useIntrinsicSlice;   // determine if using intrinsic slicing
  mat_mp B;                // The matrix for the null space for the 'fixed linears'
  mpq_t ***B_rat;
  vec_mp p;                // One vector that satisfies the original 'fixed linears'
  mpq_t **p_rat;
  vec_mp p1;               // One vector that satisfies the linears at t = 1
  mpq_t **p1_rat;
  vec_mp p0;               // One vector that satisfies the linears at t = 0
  mpq_t **p0_rat;
  mat_mp B1;               // The matrix for the null space for the linears at t = 1
  mpq_t ***B1_rat;
  mat_mp B0;               // The matrix for the null space for the linears at t = 0
  mpq_t ***B0_rat;

  point_mp *startPts;      // The startpoints of each path.
  point_mp *endPts_in;     // The endpoints - intrinsic
  point_mp *endPts;        // The endpoints - extrinisc

  comp_mp  *finalTs;       // The final T values (non-0 in case of failure).
  double  *condition_nums; // The final observed condition number of each path, for classifying endpoint types.
  int *endPt_retVals;      // The retVal from endgame for each path.
  int *endPt_types;        // The type of each endpoint - must be == MOVE_TO_NEXT to move to the next level
  int *higherDim;          // Determine if point lies on higher dim'l component
} eqStageData_mp;

typedef struct
{
  eqWitnessData_d  *witnessData_d;  // Stores all of the witness points for the subsystem (1 per subsystem).
  eqWitnessData_mp *witnessData_mp;

  eqStageData_d  *stageData_d;      // Stores all of the data for each stage
  eqStageData_mp *stageData_mp;

  comp_d  gamma_d;                  // the random gamma for the homotopy
  mpq_t   gamma_rat[2];
  comp_mp gamma_mp;
  comp_d  **coeff_d;                // the coefficients for the random linears [func_num][var_num]
  mpq_t  ***coeff_rat;
  comp_mp **coeff_mp;

  int noChanges;                    // whether we can use the SLP to get a square system to use instCounts
  int numSubFuncs;                  // number of subfunctions
  int *startSub;                    // where the subfunctions start
  int *endSub;                      // where the subfunctions end
  int *startFunc;                   // where the functions start
  int *endFunc;                     // where the functions end
  int *startJvsub;                  // where the derivs of subfunctions start
  int *endJvsub;                    // where the derivs of subfunctions end
  int *startJv;                     // where the derivs of functions start
  int *endJv;                       // where the derivs of functions end
  int **subFuncsBelow;              // matrix of which subfunctions are below each function  

  int curr_precision;               // the current precision for the multiprecision numbers - used only in AMP
  int increase_witness_prec;        // 0 - do not increase the precision for the witness data since it is no longer needed
                                    // 1 - do increase the precision for the witness data

  int **degrees;                    // mhom degrees of the function - [func_num][var_gp_num]
  int num_funcs;                    // The number of functions - we can & do assume that we are dealing with a squared-up system
  int num_subsystems;               // The number of subsystems that are used
  int num_var_gps;                  // The number of variable groups
  int num_vars;                     // The total number of variables
  int curr_stage_num;               // The array index of the current stage - used in evaluation of the function
  int curr_path_num;                // The current path number being tracked
} eqData_t;  // Usually called EqD.

// function declarations
int change_eqbyeq_eval_prec(void const *ED, int prec);

int eqbyeq_witness_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
int eqbyeq_stage_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// eqbyeq_eval.c
int witness_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int witness_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int standard_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int standard_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int stage_sort_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int stage_sort_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int eqbyeq_square_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED, int startFunc, int endFunc);
int eqbyeq_square_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED, int startFunc, int endFunc);

// clear functions
void eqbyeq_clear_d(eqData_t *EqD, int MPType);
void eqbyeq_clear_mp(eqData_t *EqD);

// eqbyeq_output.c
void eqbyeqOutputChart_d(eqData_t *EqD, FILE *fp, int infRemoved);
void eqbyeqOutputChart_mp(eqData_t *EqD, FILE *fp, int infRemoved);


#endif
