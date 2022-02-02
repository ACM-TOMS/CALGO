// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _REGEN_POS_DIM_H
#define _REGEN_POS_DIM_H

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"

// Endpoint classification codes
#define UNCLASSIFIED 0
#define SOLUTION_AND_NONSING 10
#define SOLUTION_AND_SING 15
#define NONSOLUTION_AND_NONSING 20
#define NONSOLUTION_AND_SING 25
// < 0 -> bad path (code = retVal returned from path)

typedef struct
{
  int codim;                // current codimension

  int num_paths;            // The number of paths needed to be tracked.
  int num_superset;         // number of points in the witness superset
  int num_nonsing;          // number of points in witness superset that are non-singular
  int num_sing;             // number of points in witness superset that are singular
  int num_nonsolns;         // number of nonsolutions
  int num_inf;              // number of infinite endpoints
  int num_bad;              // number of other bad endpoints

  // information for intrinsic slicing - either use B & p (intrinsic) or linears described by coeff (extrinsic)
  int useIntrinsicSlice;    // determine if this level uses intrinsic slicing
  mat_d B_d;                // The matrix for the null space
  mat_mp B_mp;
  mpq_t ***B_rat;
  vec_d p_d;                // One vector that satisfies the original linears
  vec_mp p_mp;
  mpq_t **p_rat;

  // all data related to the start points and end points will be handled locally with mass data storage using files rather than memory
} regenCodim_t;

typedef struct
{
  prog_t *Prog;          // SLP used to evaluate f
  preproc_data PPD;      // preprocessing data
 
  int *orig_degrees;     // degrees of f_i
  int *new_degrees;      // degrees of f_P[i]
  int *P;                // permutations of f_i's used to minimize paths

  int system_rank;       // rank of the original function f

  int orig_variables;    // number of variables used in 'Prog'
  int new_variables;     // number of variables after slicing
  mat_d C_d;             // matrix to convert new variables to orig variables
  mat_mp C_mp;
  mpq_t ***C_rat;

  int ***W;              // exponents of homogenizing variable to maintain the 1-hom of randomized system
  vec_d H_d;             // random vector - used to maintain the 1-hom of randomized system - (H*x + homVarConst)^W will maintain 1-hom
  vec_mp H_mp;
  mpq_t **H_rat;

  comp_d homVarConst_d;  // the constant term in the homogeneous variable expression
  comp_mp homVarConst_mp;
  mpq_t homVarConst_rat[2];

  comp_d  gamma_d;       // the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t   gamma_rat[2];

  comp_d  ***coeff_d;    // the coefficients for the random regeneration linears [codim][deg_num][var_num]
  comp_mp ***coeff_mp;
  mpq_t  ****coeff_rat;

  vec_d patchCoeff_d;    // coefficients for the patch
  vec_mp patchCoeff_mp;
  mpq_t **patchCoeff_rat;

  int sameA;             // determine if A has the same rows as A grows - avoid 2 randomizations
  mat_d *A_d;            // matrices that randomize the system - A[i].*W - used to randomize to the correct size for codim i
  mat_mp *A_mp;          // just assumed to be upper triangular 
  mpq_t ****A_rat; 

  int num_funcs;         // number of functions
  int num_codim;         // number of codimensions
  int curr_precision;    // the current precision for the multiprecision numbers - used only in AMP
  int curr_codim;        // current codimension - for function evaluation
  int moving_degree;     // degree where the linear is moving to - for function evaluation

  regenCodim_t *codim;   // stores the data for each codim
} regen_pos_dim_t;

// regen_pos_dim_track.c
int regen_pos_dim_main(witness_t *witnessSuperset, int startLevel, int maxCodim, int specificCodim, tracker_config_t *T, int pathMod, double midpoint_tol, double intrinsicCutoffMultiplier, int my_id, int num_processes, int headnode);
void printRPDFooter(regen_pos_dim_t *RPD, int codim_index, endgame_data_t *endPt, int corank, FILE *OUT, FILE *RAWOUT, FILE *FAIL, tracker_config_t *T, trackingStats *trackCount, double smallest_nonzero_SV, double largest_zero_SV);
void regen_pos_dim_SortEndpoint(int *corank, int *finite, int *soln, double condNum, regen_pos_dim_t *RPD, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT, point_d endPt_d, point_mp endPt_mp, comp_d endTime_d, comp_mp endTime_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int regen_pos_dim_setupNextStart(regen_pos_dim_t *RPD, int codim_index, tracker_config_t *T, FILE *START, FILE *NEXTSTARTPTS, FILE *PREPAREPTS, int **startPathNum);
void regen_pos_dim_RemoveBadPaths(FILE *NEXTSTARTPTS, FILE *INPUT, int *num_new_paths, int num_vars, int MPType, int numOrigPts, int numMovePts, int *startPathNum, int num_bad, int *badPaths);
int regen_pos_dim_PrepareNextCodim(int pathMod, regen_pos_dim_t *RPD, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, FILE *START, char *nextStartFile, int codim_index, int num_codim, trackingStats *trackCount, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void regen_pos_dim_seq_track(int startCodimIndex, int maxCodim, trackingStats *trackCount, int pathMod, double midpoint_tol, tracker_config_t *T, regen_pos_dim_t *RPD, char *startName, char *witName);
void regen_pos_dim_OutputChart(regen_pos_dim_t *RPD, FILE *fp, int maxCodim);
// dehomogenize the regen pos dim point
int regen_pos_dim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// regen_pos_dim_setup.c
void regen_pos_dim_setup(int startLevel, int *maxCodim, int specificCodim, tracker_config_t *T, regen_pos_dim_t *RPD, char *preprocFile, char *degreeFile, char *startFileName, double intrinsicCutoffMultiplier);
void setupRegenCodimData(regen_pos_dim_t *RPD, int codim_index, int codim, int MPType, int max_prec, int intrinsicCutoff);
int change_regen_pos_dim_prec(void const *ED, int prec);
void printRPDSummaryData(regen_pos_dim_t *RPD, int codim_index, FILE *FP);
void printRPDRelevantData(regen_pos_dim_t *RPD, int MPType, int codim_index, FILE *FP);

// regen_pos_dim_eval.c
int regen_pos_dim_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int regen_pos_dim_moving_linear_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int regen_pos_dim_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int regen_pos_dim_moving_linear_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

// regen_pos_dim_copy.c
int regen_pos_dim_copyWitness_clear(witness_t *witnessSuperset, regen_pos_dim_t *RPD, char *witName, int MPType, int max_prec, int specificCodim);
void regen_pos_dim_clear(regen_pos_dim_t *RPD, int MPType);

// regen_pos_dim_extend.c
void regenExtendMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode);
void worker_regenExtend(int my_id, int num_processes, int headnode, int dataType);

#endif

