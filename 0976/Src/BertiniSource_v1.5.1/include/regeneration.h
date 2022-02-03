// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _REGENERATION_H
#define _REGENERATION_H

#include "bertini.h"

// Endpoint classification codes - only need to know where to move this on to next level or not
#define UNCLASSIFIED 0
#define MOVE_TO_NEXT 10        // this endpoint will be moved on to next level
#define DO_NOT_MOVE_TO_NEXT 15 // this endpoint will NOT be moved on to next level
// < 0 -> bad path (code = retVal returned from path)

typedef struct 
{
  int level;                // The regeneration starting level of this data.
  int depth;                // the number of functions that will be solved for at this level starting with 'level'

  int num_paths;            // The number of paths needed to be tracked.
  int num_sing;             // The number of singular endpoints.
  int num_nonsing;          // The number of nonsingular endpoints.
  int num_inf;              // The number of infinite endpoints.
  int num_higher_dim;       // The number of endpoints that lie on higher dimensional components.
  int num_bad;              // The number of other bad paths.

  // information for intrinsic slicing - either use B & p (intrinsic) or linears described by coeff (extrinsic)
  int useIntrinsicSlice;    // determine if this level uses intrinsic slicing
  mat_d B_d;                // The matrix for the null space
  mat_mp B_mp;
  mpq_t ***B_rat;
  vec_d p_d;                // One vector that satisfies the original linears
  vec_mp p_mp;
  mpq_t **p_rat;

  // all data related to the start points and end points will be handled locally with mass data storage using files rather than memory
} regenLevel_t;

typedef struct 
{
  int num_levels;        // The number of regeneration levels used
  int num_funcs;         // number of functions
  int num_variables;     // number of variables
  int num_var_gps;       // number of variable groups

  preproc_data PPD;      // preprocessing data

  // since we need a square system, lets be general as to how we describe this square system
  void const *square_d;
  void const *square_mp;
  int (*square_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *);
  int (*square_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
  int (*change_square_prec)(void const *, int);

  int noChanges;       // whether we can use the SLP to get a square system to use instCounts
  int numSubFuncs;     // number of subfunctions
  int *startSub;       // where the subfunctions start
  int *endSub;         // where the subfunctions end
  int *startFunc;      // where the functions start
  int *endFunc;        // where the functions end
  int *startJvsub;     // where the derivs of subfunctions start
  int *endJvsub;       // where the derivs of subfunctions end
  int *startJv;        // where the derivs of functions start
  int *endJv;          // where the derivs of functions end
  int **subFuncsBelow; // matrix of which subfunctions are below each function

  // coefficients for the patches - all are of size num_var_gps x num_variables
  mat_d patchCoeff_d;
  mat_mp patchCoeff_mp;
  mpq_t ***patchCoeff_rat;
  
  // other things used for regeneration
  int **degrees;          // 'mhom' degrees of the function - [func_num][var_gp_num]
  comp_d  gamma_d;        // the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t   gamma_rat[2];
  mat_d mainCoeff_d;      // the coefficients for the main support regeneration linears [func_num][var_num] 
  mat_mp mainCoeff_mp;
  mpq_t  ***mainCoeff_rat;
  comp_d  ****coeff_d;    // the coefficients for the random regeneration linears [func_num][var_gp_num][deg_num][var_num]
  comp_mp ****coeff_mp;
  mpq_t  *****coeff_rat;

  // vector for new patch with main homogeneous variable
  vec_d main_homVar_d;     // numVars
  vec_mp main_homVar_mp;
  mpq_t **main_homVar_rat;

  int curr_precision;     // the current precision for the multiprecision numbers - used only in AMP
  int curr_level_num;     // The array index of the current level for the regeneration - used in function evaluation
  int *curr_linear;       // current linear - used in function evaluation
  int *curr_linear_degree;// current degree of the linear - used in function evaluation

  regenLevel_t *level;    // stores the data for each level of the regeneration
} regen_t;

// regen_track.c
void printRegenFooter(regen_t *regen, int level_num, endgame_data_t *endPt, int rankDef, point_d orig_d, point_mp orig_mp, point_d dehom_d, point_mp dehom_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, tracker_config_t *T, trackingStats *trackCount, int updateTrackCount, prog_t *origProg);
void regen_track_seq(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, mat_d finalPatch_d, mat_mp finalPatch_mp, mpq_t ***finalPatch_rat, prog_t *origProg);
void regenTrackLevel(int pathMod, regen_t *regen, tracker_config_t *T, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, int level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg);
void regenTrackPath(regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *MIDOUT, FILE *FAIL, FILE *NONSOLN, point_data_d *startPt_d, point_data_mp *startPt_mp, int path_num, int level_num, int *curr_linear, int *curr_linear_degree, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *), prog_t *origProg);
void regenSortLevel(int max_threads, int pathMod, regen_t *regen, tracker_config_t T[], FILE **OUT, FILE *RAWOUT, int level_num, FILE *ENDPTS, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void regenSortEndpoint(int *rankDef, int *finite, int *higherDim, double condNum, regen_t *regen, int level_num, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, int *curr_linear, int *curr_linear_degree, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int determineRegenFinite(int remove_inf, double maxNorm, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, regen_t *regen, int level_num);
int determineRegenHigherDim(double tol, double ratio, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, regen_t *regen, int level_num);
int regenPrepareNextLevel(int pathMod, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *START, char *nextStartFile, int curr_level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int regenMoveNextLevel(int max_threads, int pathMod, regen_t *regen, tracker_config_t T[], FILE **OUT, FILE **RAWOUT, FILE **MIDOUT, FILE **FAIL, char *nextStartFile, int curr_level_num, point_d *movePts_d, point_mp *movePts_mp, int *Pts_prec, trackingStats trackCount[], int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int regenLinearTrackPath(int path_num, point_data_d *startPt_d, point_data_mp *startPt_mp, int *curr_linear, int *curr_linear_degree, regen_t *regen, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *NEXTSTARTPTS, int next_level_num, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void regenRemoveBadPaths(FILE *NEXTSTARTPTS, FILE *INPUT, int *num_new_paths, int num_funcs, int num_vars, int MPType, int numOrigPts, int numMovePts, int *startPathNum, int num_bad, int *badPaths);
void regenOutputChart(regen_t *regen, FILE *fp, int infRemoved);
void regenSetupFinalSoln(regen_t *regen, tracker_config_t *T, FILE *NONSINGIN, FILE *RAWIN, FILE *FINALSOLN, mat_d patch_d, mat_mp patch_mp, mpq_t ***patch_rat);
int computeRegenRetval(regen_t *regen, int level_num, endgame_data_t *endPt, point_d orig_d, point_mp orig_mp, point_d orig_last_d, point_mp orig_last_mp, point_d dehom_d, point_mp dehom_mp, tracker_config_t *T, int junkCheck, prog_t *origProg);
void regenSortEndpoint_basic(int checkRankDef, int *rankDef, int *finite, int *higherDim, double condNum, regen_t *regen, int level_num, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int *curr_linear, int *curr_linear_degree, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void printRegenTrackFooter(int retVal, regen_t *regen, int level_num, endgame_data_t *endPt, int rankDef, int finite, int higherDim, point_d orig_d, point_mp orig_mp, point_d dehom_d, point_mp dehom_mp, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *FAIL, FILE *NONSOLN, tracker_config_t *T, trackingStats *trackCount, int updateTrackCount);
int printRegenLinearTrackFooter(int retVal, regen_t *regen, int level_num, int path_num, endgame_data_t *endPt, int rankDef, int finite, int higherDim, point_d orig_d, point_mp orig_mp, point_d dehom_d, point_mp dehom_mp, FILE *OUT, FILE *RAWOUT, FILE *NEXTSTARTPTS, FILE *FAIL, tracker_config_t *T);

// dehomogenize the regeneration point
int regen_dehom_stack(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
int regen_dehom_hom(point_d hom_d, point_mp hom_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
int regen_dehom2(point_d nonhom_d, point_mp nonhom_mp, point_d hom_d, point_mp hom_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
int regen_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// regen_setup.c
int change_regen_prec(void const *RED, int prec);
void setupRegenLevels(regen_t *regen, tracker_config_t *T, char *depthName, int intrinsicCutoff);
void setupFirstLevel(regen_t *regen, tracker_config_t *T, char *startName);
void setupRegenIntrinsicSlice(regen_t *regen, int MPType, int level_num, int max_prec);
int createFirstRegenStartPoints(regen_t *regen, tracker_config_t *T, char *startName);
int createRegenStartPts(regen_t *regen, int MPType, int count, int **decomp, FILE *START);
void setupRegenRestart(regen_t *regen, tracker_config_t *T, char *startName, int curr_level);
void setupRegenLevelStructures(regen_t *regen, int MPType, int level_num, int max_prec);
void clearRegenLevelStructures(int max_threads, regen_t *regen, int MPType, int level_num);
void copyRegenLevelStructures(regen_t *regenIn, int max_threads, regen_t *regen_copy, int MPType, int level_num);
void copyRegenRandom(regen_t *regen, regen_t *regenIn, int MPType, void const *square_d, void const *square_mp);
void printRegenSummaryData(regen_t *regen, int finished_level, FILE *FP);
void printRegenRelevantData(regen_t *regen, int MPType, int curr_level, FILE *FP);
void printPointLinearDegree(FILE *FP, point_d Pt_d, point_mp Pt_mp, int prec, int *linear, int *degree, int size);
void createDecomp(int ***goodLoc, int *good_size, regen_t *regen, int startFunc, int endFunc);

// regen_eval.c
int regen_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int regen_square_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int regen_moving_linear_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);

int regen_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int regen_square_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int regen_moving_linear_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

// regen_parallel.h
void regenClassifyLevel(regen_t *regen, int curr_level, double final_tol, int MPType, FILE *RAWOUT, FILE *RAWIN);
int head_regenMoveOrigPts(int minPacketSize, int maxPacketSize, int numOrigPts, int pathsTracked, int pathMod, regen_t *regen, tracker_config_t *T, int **new_linear, int **new_linear_degree, int new_count, mat_d B_transpose_d, mat_mp B_transpose_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *START, FILE *NEXTSTARTPTS, int level_num, trackingStats *trackCount, int (*change_prec)(void const *, int), int my_id, int headnode, int num_processes);

int regen_worker_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

#endif
