// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _DIMBYDIM_H
#define _DIMBYDIM_H

#include "bertini.h"
#include "cascade.h"

// Endpoint classification codes
#define UNCLASSIFIED 0
#define NON_SINGULAR 10
#define SINGULAR 15
// < 0 -> bad path (code = retVal returned from path)

// Data structures
typedef struct
{
  int codim;                     // current codimension

  int num_paths;                 // number of paths to be tracked for this codim
  int num_superset;              // number of points in the witness superset
  int num_nonsing;               // number of points in witness superset that are non-singular
  int num_sing;                  // number of points in witness superset that are singular
  int num_nonsolns;              // number of nonsolutions
  int num_inf;                   // number of infinite endpoints
  int num_bad;                   // number of other bad endpoints

  mat_d A_d;                     // random matrix - randomize the system to the proper codim
  mat_mp A_mp;
  mpq_t ***A_rat;
  int **W;                       // exponents of homogenizing variable to maintain the 1-hom of randomized system

  vec_d H_d;                     // random vector - used to maintain the 1-hom of randomized system - (H*x + homVarConst)^W will maintain 1-hom
  vec_mp H_mp;
  mpq_t **H_rat;

  comp_d homVarConst_d;          // the constant term in the homogeneous variable expression
  comp_mp homVarConst_mp;
  mpq_t homVarConst_rat[2];

  int useIntrinsicSlice;         // determine if this codim uses intrinsic slicing
  // the setup of B & p depend upon if this codim uses intrinsic or extrinsic slices
  // either way, they will be used
  mat_d B_d;                     // random matrix - 'coefficients' (ext) or 'null space' (int)
  mat_mp B_mp;
  mpq_t ***B_rat;

  vec_d p_d;                     // random vector - 'patch coefficients' (ext) or 'solution vector' (int)
  vec_mp p_mp;
  mpq_t **p_rat;

  point_d *startPts_d;           // start point for each path - double precision
  point_mp *startPts_mp;         // start point for each path - multi precision

  endpoint_data_d *endPts_d;     // end point data for each path - double precision
  endpoint_data_mp *endPts_mp;   // end point data for each path - multi precision
  endpoint_data_amp *endPts_amp; // end point data for each path - AMP

  int *endPt_types;              // end point types
} codimData_t;

typedef struct
{
  prog_t *Prog;         // slp to evaluate the original function f
  preproc_data PPD;     // preprocessing data

  int *orig_degrees;    // degrees of f_i
  int *new_degrees;     // degrees of f_P[i]
  int *P;               // permutations of f_i's used to minimize paths

  int system_rank;      // rank of the original function f

  int orig_variables;   // number of variables used in 'Prog'
  int new_variables;    // number of variables after slicing
  mat_d C_d;            // matrix to convert new variables to orig variables
  mat_mp C_mp;  
  mpq_t ***C_rat;
  
  comp_d gamma_d;       // the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t gamma_rat[2];

  int curr_precision;   // current precision

  int num_funcs;        // number of functions
  int num_codim;        // number of codimensions
  int curr_codim_index; // current codimension

  codimData_t *codim;   // stores the codimension data
} codim_t;

// dimbydim_track.c
int printDimbyDimFooter_d(codim_t *CD, int codim_index, int path_num, point_data_d *endPoint, point_d orig_vars, point_d dehomP, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount);
int printDimbyDimFooter_mp(codim_t *CD, int codim_index, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp dehomP, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

int determineDimbyDimFinite_d(double max_norm, point_d point, codim_t *CD, int codim_index);
int determineDimbyDimFinite_mp(double max_norm, point_mp point, int prec, codim_t *CD, int codim_index);

int determineDimbyDimSoln_d(double tol, double ratio, point_d point, point_d last_point, comp_d time, codim_t *CD, int codim_index);
int determineDimbyDimSoln_mp(double tol, double ratio, point_mp point, point_mp last_point, comp_mp time, int prec, codim_t *CD, int codim_index);

void dimbydim_classifyCodim(codim_t *CD, int codim_index, double final_tol, int MPType);

void dimbydimSortEndpoint(int *rankDef, int *finite, int *soln, double condNum, codim_t *CD, int codim_index, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int (*change_prec)(void const *, int));

// dehomogenize the dimbydim point
int dimbydim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// dimbydim_parallel.c
void dimbydim_par_track(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, double midpoint_tol, tracker_config_t *T, codim_t *CD, int my_id, int num_processes, int headnode);
int store_dimbydim_endPoint(endgame_data_t *endPt, int corank, double smallest, double largest, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, void const *CD_d, void const *CD_mp);

// dimbydim_setup.c
void dimbydim_setup(FILE **OUT, char *outName, FILE **RAWOUT, char *rawName, FILE **MIDOUT, char *midName, FILE **FAIL, char *failName, tracker_config_t *T, codim_t *CD, char *preprocFile, char *degreeFile, double intrinsicCutoffMultiplier, int maxCodim, int specificCodim);

void dimbydim_clear(codim_t *CD, int MPType);
void clearCodimData(codimData_t *CD, int MPType);

void setup_omp_codim_t(codim_t *CD, codim_t *CD_in, int MPType);

void setup_dimbydim_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, tracker_config_t **T_copy, tracker_config_t *T, codim_t **CD_copy, codim_t *CD);

void clear_omp_codim_t(codim_t *CD, int MPType);

void clear_dimbydim_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, tracker_config_t **T_copy, codim_t **CD_copy);

void dimbydim_clear_start_points(codim_t *CD, int codim_index, int MPType);

int change_dimbydim_prec(void const *ED, int prec);

// dimbydim_eval.c
int standard_dimbydim_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int standard_dimbydim_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

// dimbydim_output.c
void dimbydimOutputChart(codim_t *CD, FILE *fp, int MPType);

void dimbydimFindOrigVarsDehom_d(point_d orig_vars_d, point_d dehom_d, point_d P_d, codim_t *CD, int codim_index);
void dimbydimFindOrigVarsDehom_mp(point_mp orig_vars_mp, point_mp dehom_mp, point_mp P_mp, codim_t *CD, int codim_index);

// parallel_send_functions.c
void bcast_codim_t(codim_t *CD, int MPType, int my_id, int headnode);
void bcast_codimData_t(codimData_t *CD, int curr_prec, int MPType, int my_id, int headnode);

// copy_functions.c
void cp_codim_int(void *Out, void *In, int MPType, char **codimStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType);
void cp_codimData_int(void *Out, void *In, int curr_prec, int MPType, char **codimStr, int freeStr, comp_d **coeff_d, int **W_send, int inType);

#endif
