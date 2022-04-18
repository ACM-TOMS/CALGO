// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _PARALLEL_H
#define _PARALLEL_H

#include "bertini.h"
#include "cascade.h"
#include "regeneration.h"
#include "regen_pos_dim.h"
#include "eqbyeq.h"

#ifdef _HAVE_MPI
#include <mpi.h>
#endif

#define STOPCODE       -199 // tells the worker process that it can stop
#define ZERO_DIM_D     -198 // tells the worker process that it needs to do standard zero dimensional tracking in either double or AMP
#define USER_HOM_D     -197 // tells the worker process that it needs to do user defined homotopy zero dimensional tracking in either double or AMP
#define ZERO_DIM_MP    -196 // tells the worker process that it needs to do standard zero dimensional tracking in fixed multi precision
#define USER_HOM_MP    -195 // tells the worker process that it needs to do user defined homotopy zero dimensional tracking in fixed multi precision
#define ZERO_DIM_REGEN -194 // tells the worker process that it needs to do zero dimensionial tracking using regeneration 
#define EQBYEQ_D       -193 // tells the worker process that it needs to do equation-by-equation in either double or adaptive precision
#define EQBYEQ_MP      -192 // tells the worker process that it needs to do equation-by-equation in fixed multi precision
#define IRREDDECOMP    -191 // tells the worker process that it needs to do numerical irreducible decomposition
#define DIMBYDIM       -190 // tells the worker process that it needs to do dimension-by-dimension tracking
#define CASCADE        -189 // tells the worker process that it needs to do cascade tracking
#define REGEN_POS_DIM  -188 // tells the worker process that it needs to do pos-dim regeneration
#define PARAM_HOM_D    -187 // tells the worker process that it needs to do parameter homtopy
#define PARAM_HOM_MP   -186 // tells the worker process that it needs to do parameter homtopy
#define REGEN_EXTEND   -185 // tells the worker process that it needs to do regeneration extension 

#define TAG_NUM_PTS        1 // tag for sending the number of points
#define TAG_POINT_DATA_D   2 // tag for sending point_data_d
#define TAG_ENDGAME_DATA_T 3 // tag for sending endgame_data_t_int
#define TAG_STR            4 // tag for sending string
#define TAG_ENDGAME_STR    TAG_STR // tag for sending the string for endgame_data_t_int
#define TAG_COEFF          5 // tag for sending coefficients
#define TAG_ENDGAME_COEFF  TAG_COEFF // tag for sending the coefficients for endgame_data_t_int
#define TAG_RANK_TYPE      6 // tag for sending the rank type
#define TAG_POINT_D_INT    7 // tag for sending point_d_int
#define TAG_POINT_MP_INT   8 // tag for sending point_mp_int
#define TAG_TRACKBACK_T    9 // tag for sending trackBack_samples_t_int
#define TAG_DOUBLE        10 // tag for sending doubles
#define TAG_INT           11 // tag for sending ints
#define TAG_OTHER_DATA    12 // tag for sending other data

typedef struct { // this can be expanded when more information about the task needs to be sent to the worker process
  int dataType;  // tells the worker process which parallel task it needs to run
} worker_info;

typedef struct
{
  int prec;        // current precision
  int totalLength; // string length with '\0' at the end
} mpf_int;

typedef struct
{
  int prec;         // current precision
  int totalLength;  // string length for r & i with a '\0' in between and at end
} comp_mp_int;

typedef struct
{
  int length[2];    // string length for r & i with a '\0' in between and at end
} comp_rat_int;

typedef struct
{
  int size;         // size of the point
} point_d_int; 

typedef struct 
{
  int prec;         // current precision
  int size;         // size of the point
  int totalLength;  // string length to reconstruct each comp_mp - separated by '\0'
} point_mp_int;

typedef struct
{
  int size;         // size of the vec
  int totalLength;  // string length to reconstruct each comp_mp - separated by '\0'
} point_rat_int;

typedef struct
{
  point_d_int point_int;  // integer information to reconstruct 'point' - time is put at bottom of comp_d array for sending point
  comp_d time;            // time
  int cycle_num;          // cycle number
} point_data_d_int;

typedef struct 
{
  point_mp_int point_int; // integer information to reconstruct 'point'
  comp_mp_int  time_int;  // integer information to reconstruct 'time'
  int cycle_num;          // cycle number
  int totalLength;        // total number of characters for both point & time
} point_data_mp_int;

typedef struct
{
  int rows;        // number of rows
  int cols;        // number of cols
} mat_d_int;

typedef struct
{
  int prec;        // current precision
  int rows;        // number of rows
  int cols;        // number of cols
  int totalLength; // string length to reconstruct each comp_mp - separated by '\0'
} mat_mp_int;

typedef struct
{
  int rows;        // number of rows
  int cols;        // number of cols
  int totalLength; // string length to reconstruct each comp_mp - separated by '\0'
} mat_rat_int;

typedef struct
{   
  int prec;   
  point_data_d_int PD_d_int;
  point_data_mp_int PD_mp_int;

  int last_approx_prec;       // precision of the last approximation
  point_d_int last_approx_d_int;      // last approximation to the end point
  point_mp_int last_approx_mp_int;    // last approximation to the end point

  int retVal;   
  int pathNum;   
  int codim;
  double first_increase;
  double condition_number;
  double  function_residual_d;
  mpf_int function_residual_mp_int;
  double  latest_newton_residual_d;   
  mpf_int latest_newton_residual_mp_int;   
  double  t_val_at_latest_sample_point_d;   
  mpf_int t_val_at_latest_sample_point_mp_int;
  double  error_at_latest_sample_point_d;
  mpf_int error_at_latest_sample_point_mp_int; 
  int totalLength; // string length to reconstruct all MP datatypes
} endgame_data_t_int; 

typedef struct
{
  endgame_data_t_int EG_int;

  int numSamples;
  int samplePts_prec;
  int midPt_prec;
  int pointSize;

  int num_double;
  int num_comp_d;
  int totalLength;
} trackBack_samples_t_int;

// only the relevant data in tracker_config_t
typedef struct {
  int numVars;
  int numPathVars;
  int numParams;
  int numFuncs;

  double maxStepSize;
  double minStepSizeBeforeEndGame;
  double minStepSizeDuringEndGame;

  double minTrackT; 
  double basicNewtonTol; 
  double endgameNewtonTol;

  int cSecInc;   
  int maxNewtonIts;
  int MPType; 
  int Precision; 
  int outputLevel;
  int screenOut; 
  double targetT;
  double endgameBoundary; 

  double goingToInfinity;
  int maxNumSteps;
  int endgameNumber;

  double power_series_sample_factor;  
  int cycle_num_max; 
  int num_PSEG_sample_points; 

  double final_tolerance; 
  double real_threshold;  
  int endgameOnly; 

  double AMP_bound_on_abs_vals_of_coeffs; 
  double AMP_bound_on_degree; 
  double AMP_eps; 
  double AMP_Phi; 
  double AMP_Psi;
  int AMP_safety_digits_1;
  int AMP_safety_digits_2;
  int AMP_max_prec;

  double sing_val_zero_tol;  
  double cond_num_threshold;

  double step_fail_factor;
  double step_success_factor;

  int max_num_pts_for_trace;
  int max_num_mon_linears;
  int max_num_bad_loops_in_mon;

  double final_tol_multiplier;
  double final_tol_times_mult;

  int sharpenOnly;
  int sharpenDigits; 

  int regen_remove_inf;
  int regen_higher_dim_check;
  double sliceBasicNewtonTol;
  double sliceEndgameNewtonTol;
  double sliceFinalTol;

  int minCycleTrackBack;       
  int junkRemovalTest;         
  int maxDepthLDT;
  int odePredictor;            

  int securityLevel;
  double securityMaxNorm;

  double cutoffCycleTime;
  double cutoffRatioTime; 

  double finiteThreshold;

  double funcResTol;
  double ratioTol;

  int maxStepsBeforeNewton;

} tracker_config_t_relevant;

typedef struct 
{
  int size; 
  int memSize;
  int precision;

  int num_var_gps; 
  int index_of_first_number_for_proj_trans; 

  int numInstAtEndUpdate; 
  int numInstAtEndParams;
  int numInstAtEndPDeriv;
  int numInstAtEndFnEval; 
  int numInstAtEndJvEval; 

  int numVars;
  int numPathVars;
  int numNums;    
  int numConsts;

  int numPars; 
  int numFuncs;
  int numSubfuncs;

  int inpVars;
  int inpPathVars;
  int IAddr;     
  int numAddr;  
  int constAddr;

  int evalPars;
  int evalDPars;
  int evalFuncs;
  int evalJVars;
  int evalJPars;
  int evalSubs; 
  int evalJSubsV;
  int evalJSubsP;

  int totalLength; // total length for the strings to create the nums
} prog_t_int;

typedef struct
{
  int num_funcs;
  int num_hom_var_gp;
  int num_var_gp;
} preproc_data_int;

typedef struct
{
  int num_patches;
  int patchCoeff_rows;
  int patchCoeff_cols;
} patch_eval_data_d_int;

typedef struct
{
  int prec;
  int num_patches;
  int patchCoeff_rows;
  int patchCoeff_cols;
  int totalLength;
} patch_eval_data_mp_int; 

typedef struct
{
  int startSystemType; 
  int size_r;
  int max_degree; 
  int coeff_cols;
  int num_coeff;
  comp_d gamma;
} start_system_eval_data_d_int; 

typedef struct
{
  int startSystemType; 
  int size_r; 
  int max_degree; 
  int coeff_cols;
  int coeff_strLength;
  int prec;
  int totalLength;
} start_system_eval_data_mp_int;

typedef struct
{ 
  prog_t_int Prog_int;
  int size_f;
  int B_rows;
  int B_cols;
  int B_perp_rows;
  int B_perp_cols;
  int noChanges;
  int max_of_W;
  int A_rows;
  int A_cols;
  int size_r;
  int num_comp_d;    
} square_system_eval_data_d_int; 

typedef struct
{ 
  prog_t_int Prog_int;
  int size_f;
  int B_rows;
  int B_cols;
  int B_strLength;
  int B_perp_rows;
  int B_perp_cols;
  int B_perp_strLength; 
  int noChanges;
  int max_of_W;
  int A_rows; 
  int A_cols;
  int A_strLength;
  int size_r;
  int prec;
  int totalLength;
} square_system_eval_data_mp_int; 

typedef struct
{
  square_system_eval_data_mp_int squareSystem_int;
  patch_eval_data_mp_int patch_int;
  start_system_eval_data_mp_int startSystem_int;
  preproc_data_int preProcData_int;
} basic_eval_data_mp_int;

typedef struct
{
  square_system_eval_data_d_int squareSystem_int;
  patch_eval_data_d_int patch;
  start_system_eval_data_d_int startSystem_int;
  preproc_data_int preProcData_int;
  int MPType; // if == 2, need to send BED_mp
} basic_eval_data_d_int;

typedef struct
{
  comp_d gamma_d; 

  int curr_precision;    // the current precision for the multiprecision numbers - used only in AMP

  int num_funcs;         // The number of functions - we can & do assume that we are dealing with a squared-up system
  int num_subsystems;    // The number of subsystems used
  int num_var_gps;       // The number of variable groups
  int num_vars;          // The total number of variables

  int num_coeff;         // the total number of coeff that are to be sent - used only in double precision
  int totalLength;       // the total length of the string needed to send the rest of the data

  int noChanges;         // whether we can use the SLP to get a square system to use instCounts
  int numSubFuncs;       // number of subfunctions
  int numInts;           // number of integers needed to setup information about location of functions/subfunctions
} eqData_t_int;

typedef struct
{
  int curr_precision;// the current precision for the multiprecision numbers - used only in AMP

  int startFunction; // the starting function for witness data
  int depth;         // the number of functions whose witness data will be found starting with 'startFunction'
  int num_paths;     // The number of paths stored.

  int num_linears;   // The total number of linears used in the start system

  int B_rows;        // The number of rows in B
  int B_cols;        // The number of cols in B
  int p_size;        // The number of entries in p

  int num_comp_d;    // the total number of comp_d that are to be sent - used only in double precision
  int totalLength;   // the total length of the string needed to send the rest of the data
} witnessData_t_int;

typedef struct
{
  int curr_precision;    // the current precision for the multiprecision numbers - used only in AMP

  int depth_x;           // the number of functions, starting at 0, that the x-variables satisfy
  int depth_y;           // the number of functions, starting at depth_x, that the y-variables satisfy
  int num_paths;         // The number of paths stored.

  int useIntrinsicSlice; // determine if this level uses intrinsic slicing
  int B_rows;            // The number of rows in B
  int B_cols;            // The number of cols in B
  int p_size;            // The number of entries in p
  int Bt_rows;           // The number of rows in B1 & B0
  int Bt_cols;           // The number of cols in B1 & B0
  int pt_size;            // The number of entries in p1 & p0

  int num_comp_d;        // the total number of comp_d that are to be sent - used only in double precision
  int totalLength;       // the total length of the string needed to send the rest of the data
} stageData_t_int;

typedef struct
{
  int curr_precision;       // the current precision for the multiprecision numbers

  prog_t_int Prog_int;      // the slp used to evaluate the original system f
  preproc_data_int PPD_int; // preprocessing data

  int num_funcs;            // number of functions in f
  int system_rank;          // rank of f
  int num_codim;            // number of codimensions 

  int orig_variables;       // number of original variables used to evaluate f
  int new_variables;        // number of variables after slicing

  int num_comp_d;           // number of comp_d that are to be sent - used only in double preicison
  int totalLength;          // total length of the string needed to send the rest of the data
} codim_t_int;
  
typedef struct
{
  int curr_precision;       // the current precision for the multiprecision numbers

  int codim;                // the codimension this data represents
  int A_W_rows;             // number of rows in A & W
  int A_W_cols;             // number of cols in A & W
 
  int H_size;               // size of H

  int useIntrinsicSlice;    // determine if this codim uses intrinsic slicing

  int B_rows;               // number of rows in B
  int B_cols;               // number of cols in B
  int p_size;               // size of p

  int num_comp_d;           // number of comp_d that are to be sent - used only in double preicison
  int totalLength;          // total length of the string needed to send the rest of the data
} codimData_t_int;

typedef struct
{
  int curr_precision;       // the current precision for the multiprecision numbers

  prog_t_int Prog_int;      // the slp used to evaluate the original system f
  preproc_data_int PPD_int; // preprocessing data

  int num_funcs;            // number of functions in f
  int system_rank;          // rank of f
  int num_codim;            // number of codimensions

  int orig_variables;       // number of original variables used to evaluate f
  int new_variables;        // number of variables after slicing

  int A_W_rows;             // number of rows in A & W
  int A_W_cols;             // number of cols in A & W

  int H_size;               // size of H

  int C_rows;               // number of rows in C
  int C_cols;               // number of cols in C

  int R_rows;               // number of rows in R
  int R_cols;               // number of cols in R

  int T_size;               // size of T

  int B_rows;               // number of rows in B
  int B_cols;               // number of cols in B
  int p_size;               // size of p

  int num_int;              // number of int that are to be sent
  int num_comp_d;           // number of comp_d that are to be sent - used only in double preicison
  int totalLength;          // total length of the string needed to send the rest of the data
} cascade_t_int;

typedef struct
{
  int codim;                // the codimension this data represents
} cascadeCodim_t_int;

typedef struct
{
  prog_t_int Prog_int; // slp to evaluate f
  
  int curr_codim; // codimension that this represents
  int orig_variables; // original number of variables (extrinsically)
  int curr_precision; // current precision

  int A_rows;  // random matrix - randomize the system to the proper codim
  int A_cols;

  int B_rows;  // random matrix - coefficients for linear slices
  int B_cols;

  int p_size;  // random vector - patch coefficients

  int startSliceVec_size; // vector of where to start at
  int startSliceVec_init; // whether startSliceVec has been initialized

  int targetSliceVec_size; // vector of where to finish at
  int targetSliceVec_init; // whether targetSliceVec has been initialized

  int K_rows; // kernel of the linear space defined by B & p
  int K_cols;

  int num_comp_d;  // number of comp_d that are to be sent - used only in double precision
  int totalLength; // total length of the string needed to send the rest of the data
} membership_slice_moving_t_int;

typedef struct
{
  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value

  int pt_size;                // number of coordinates
  int last_approx_size;       // number of coordinates
  int num_comp_d;             // number of comp_d that are to be sent
} endpoint_data_d_int;

typedef struct
{
  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value

  int pt_size;                // number of coordinates
  int last_approx_size;       // number of coordinates
  int totalLength;            // total length of the string needed to send the rest of the data
} endpoint_data_mp_int;

typedef struct
{
  int curr_prec;              // current precision
  int last_approx_prec;       // precision of the last approximation

  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value

  int pt_size;                // number of coordinates
  int last_approx_size;       // number of coordinates
  int num_comp_d;             // number of comp_d that are to be sent
  int totalLength;            // total length of the string needed to send the rest of the data
} endpoint_data_amp_int;

typedef struct
{
  int codim;                         // current codimension
  int num_set;                       // number of points in the witness set
  int num_nonsing;                   // number of points in witness set that are non-singular
  int num_sing;                      // number of points in witness set that are singular

  int A_rows;                        // random matrix - randomize the system to the proper codim
  int A_cols;

  int H_size;                        // random vector - used to maintain the 1-hom of randomized system - (H*x + homVarConst)^W will maintain 1-hom

  int B_rows;                        // random matrix - coefficients for linear slices
  int B_cols;

  int p_size;                        // random vector - patch coefficients

  int curr_prec;
  int num_comp_d;
  int totalLength;

} witnessCodim_t_int;

typedef struct
{
  prog_t_int Prog_int;      // the slp used to evaluate the original system f
  preproc_data_int PPD_int; // preprocessing data

  int num_funcs;         // number of functions
  int num_codim;         // number of codimensions
  int curr_precision;    // current precision

  int system_rank;       // rank of the original function f

  int orig_variables;    // number of variables used in 'Prog'
  int new_variables;     // number of variables after slicing

  int num_comp_d;
  int totalLength;
} witness_t_int;

// main tracking functions to do MPI tracking - parallel_main.c
void head_zero_dim_endgame(int useTrackBack, int minCycleTrackBack, trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode);

void head_zero_dim_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int my_id, int num_processes, int headnode);

void head_zero_dim_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int my_id, int num_processes, int headnode);

void head_eqbyeq_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode);
void head_eqbyeq_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_mp *ED, int my_id, int num_processes, int headnode);

void worker_process_main(int my_id, int num_processes, int headnode);

void worker_zero_dim_endgame(int my_id, int num_processes, int headnode, int dataType);

void worker_eqbyeq_d(int my_id, int num_processes, int headnode, int dataType);
void worker_eqbyeq_mp(int my_id, int num_processes, int headnode, int dataType);

void delete_parallel_files(char *fileName, int headnode, int num_processes);
void combine_midpath_data(char *midFile, char *parFile, int delete_files, int headnode, int num_processes);
int parallel_midpoint_checking(char *midFile, char *parFile, int delete_files, int num_paths, int numVars, double midpoint_tol, int my_id, int headnode, int num_processes);

// positive dimensional parallel worker functions
void worker_numericalIrredDecomp(int my_id, int num_processes, int headnode, int dataType);
void worker_dimbydim(int my_id, int num_processes, int headnode, int dataType);
void worker_cascade(int my_id, int num_processes, int headnode, int dataType);

void packetSize_maker(int *packetSizes, int num_procs, int headnode, int totalToDistribute, int currProc, int minSize, int maxSize, int *max);

void cp_patch_d_int(void *PED_out, void *PED_in, comp_d **coeff, int inType);
void cp_patch_mp_int(void *PED_out, void *PED_in, char **patchStr, int freeStr, int inType);

void head_trackPaths(char *strJob, int num_paths, int minPacketSize, int maxPacketSize, trackingStats *trackCount, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*setup_startPoint)(endgame_data_t *, int, int, void const *, void const *), int (*store_endPoint)(endgame_data_t *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, void const *, void const *));

void worker_trackPaths(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*worker_track)(endgame_data_t *, endgame_data_t *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *)));

void worker_sortPaths(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*setup_worker_track)(int, int, void const *, void const *), int (*classifyEndpoint)(endgame_data_t *, tracker_config_t *, FILE *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int)));

int worker_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int worker_track_path_rank(int rankType, int *rankDef, int *corank, double *smallest_nonzero_sv, double *largest_zero_sv, endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int worker_useSharpener(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp);
int worker_setup(int pathNum, int MPType, void const *ED_d, void const *ED_mp);

#ifdef _HAVE_MPI // the following include MPI_Datatype and so only need to be declared when MPI is being used

// functions to create the MPI datatypes - parallel_datatypes.c
void create_worker_info(MPI_Datatype *mpi_worker_info);
void create_comp_d(MPI_Datatype *mpi_comp_d);
void create_comp_rat_int(MPI_Datatype *mpi_comp_rat_int);
void create_point_d_int(MPI_Datatype *mpi_point_d_int);
void create_comp_point_d_int(MPI_Datatype *mpi_comp_d, MPI_Datatype *mpi_point_d_int);
void create_point_data_d_int(MPI_Datatype *mpi_point_data_d_int);
void create_mat_d_int(MPI_Datatype *mpi_mat_d_int);
void create_mpf_int(MPI_Datatype *mpi_mpf_int);
void create_comp_mp_int(MPI_Datatype *mpi_comp_mp_int);
void create_point_mp_int(MPI_Datatype *mpi_point_mp_int);
void create_comp_point_mp_int(MPI_Datatype *mpi_comp_mp_int, MPI_Datatype *mpi_point_mp_int);
void create_point_data_mp_int(MPI_Datatype *mpi_point_data_mp_int);
void create_mat_mp_int(MPI_Datatype *mpi_mat_mp_int);
void create_point_rat_int(MPI_Datatype *mpi_point_rat_int);
void create_mat_rat_int(MPI_Datatype *mpi_mat_rat_int);
void create_endgame_data_t_int(MPI_Datatype *mpi_endgame_data_t);
void create_trackBack_samples_t_int(MPI_Datatype *mpi_trackBack_samples_t);
void create_tracker_config_t(MPI_Datatype *mpi_tracker_config_t_relevant);
void create_prog_t_int(MPI_Datatype *mpi_prog_t_int);
void create_preproc_data_int(MPI_Datatype *mpi_preproc_data_int);
void create_patch_eval_data_d_int(MPI_Datatype *mpi_patch_eval_data_d_int);
void create_patch_eval_data_mp_int(MPI_Datatype *mpi_patch_eval_data_mp_int);
void create_start_system_eval_data_d_int(MPI_Datatype *mpi_ssed_d_int);
void create_start_system_eval_data_mp_int(MPI_Datatype *mpi_ssed_mp_int);
void create_square_system_eval_data_d_int(MPI_Datatype *mpi_ssed_d_int);
void create_square_system_eval_data_mp_int(MPI_Datatype *mpi_ssed_mp_int);
void create_basic_eval_data_d_int(MPI_Datatype *mpi_bed_d_int);
void create_basic_eval_data_mp_int(MPI_Datatype *mpi_bed_mp_int);
void create_eqData_t_int(MPI_Datatype *mpi_eqd_t_int);
void create_witnessData_t_int(MPI_Datatype *mpi_wd_t_int);
void create_stageData_t_int(MPI_Datatype *mpi_sd_t_int);
void create_codim_t_int(MPI_Datatype *mpi_codim_t_int);
void create_codimData_t_int(MPI_Datatype *mpi_codimData_t_int);
void create_cascade_t_int(MPI_Datatype *mpi_cascade_t_int);
void create_cascadeCodim_t_int(MPI_Datatype *mpi_cascadeCodim_t_int);

void create_regen_int(MPI_Datatype *mpi_regen_int);
void create_regenLevel_int(MPI_Datatype *mpi_regenLevel_int);

void create_membership_slice_moving_t_int(MPI_Datatype *mpi_slice_moving_int);
void create_endpoint_data_d_int(MPI_Datatype *mpi_endpoint_d_int);
void create_endpoint_data_mp_int(MPI_Datatype *mpi_endpoint_mp_int);
void create_endpoint_data_amp_int(MPI_Datatype *mpi_endpoint_amp_int);
void create_witnessCodim_t_int(MPI_Datatype *mpi_witnessCodim_int);
void create_witness_t_int(MPI_Datatype *mpi_witness_int);

#endif

//////// REGENERATION SECTION /////////////

typedef struct
{
  int num_levels;
  int num_funcs;
  int num_variables;
  int num_var_gps;

  preproc_data_int PPD_int;

  int patch_rows;
  int patch_cols;

  int curr_precision;

  int num_comp_d;  // total number of comp_d that are to be sent - double precision
  int totalLength; // total length of the string needed to setup the data

  int noChanges;   // whether we can use the SLP to get a square system to use instCounts
  int numSubFuncs; // number of subfunctions
  int numInts;     // number of integers
} regen_t_int;

typedef struct
{
  int level;
  int depth;

  int useIntrinsicSlice;
  int B_rows;
  int B_cols;
  int p_size;

  int num_comp_d;
  int totalLength;
} regenLevel_t_int;

void head_trackPaths2(char *strJob, int randomize_paths, int num_paths, int minPacketSize, int maxPacketSize, trackingStats *trackCount, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int), FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, int my_id, int headnode, int num_processes, int (*create_send_packet)(int, int, FILE *, endgame_data_t *, int *, int, int, void const *, void const *, char *, int, int), int (*recv_store_packet)(endgame_data_t **, int *, trackingStats *, tracker_config_t *, FILE *, FILE *, FILE *, FILE *, FILE *, int *, void const *, void const *, int (*change_prec)(void const *, int)));
void worker_trackPaths2(trackingStats *trackCount, tracker_config_t *T, void const *ED_d, void const *ED_mp, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *), int (*recv_track_send_packet)(int, endgame_data_t **, int *, endgame_data_t **, int *, FILE *, FILE *, tracker_config_t *, void const *, void const *, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *)));

// regen_parallel.h
void head_regen_track_zero_dim(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int my_id, int num_processes, int headnode);
void head_regen_track(int startLevel, regen_t *regen, tracker_config_t *T, int pathMod, double midpoint_tol, char *startName, FILE *OUT, char *midFile, FILE *FAIL, FILE *FINALSOLN, trackingStats *trackCount, mat_d finalPatch_d, mat_mp finalPatch_mp, mpq_t ***finalPatch_rat, int my_id, int num_processes, int headnode);
void head_regenTrackLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int curr_level, int pathMod, tracker_config_t *T, regen_t *regen, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *RAWSORT, FILE *FAIL, FILE *NONSOLN, int my_id, int headnode, int num_processes);
void head_regenSortLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int curr_level, int pathMod, tracker_config_t *T, regen_t *regen, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
int head_regenPrepareNextLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int curr_level, int pathMod, tracker_config_t *T, regen_t *regen, FILE *START, char *nextStartFile, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void head_regenMoveNextLevel(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int next_level_num, int pathMod, tracker_config_t *T, regen_t *regen, FILE *NEXTSTARTPTS, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int newPaths, FILE *START, int *retVals, int my_id, int headnode, int num_processes);
void worker_regen_zero_dim(int my_id, int num_processes, int headnode, int dataType);
void worker_regen_track(regen_t *regen, tracker_config_t *T, trackingStats *trackCount, int my_id, int num_processes, int headnode);
void worker_regenTrackLevel(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_t *regen, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void worker_regenSortLevel(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_t *regen, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void worker_regenPrepareLevel(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_t *regen, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void cp_regen_int(void *regen_out, void *regen_in, int MPType, int **Degrees, int **PPD_type, int **PPD_size, char **regenStr, comp_d **coeff_d, int freeStr, int **instCount, int inType);
void cp_regenLevel_int(void *rL_out, void *rL_in, int MPType, int curr_prec, char **rlStr, comp_d **coeff_d, int freeStr, int inType);
void bcast_regen_t(regen_t *regen, int MPType, int my_id, int headnode);
void bcast_regenLevel_t(regenLevel_t *rL, int MPType, int curr_prec, int my_id, int headnode);
int regen_create_send_packet_track(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int regen_create_send_packet_move(int startNum, int size, FILE *START, endgame_data_t *sendPts, int *pathNum, int MPType, int pathMod, void const *ED_d, void const *ED_mp, char *jobName, int totalPaths, int sendProc);
int regen_recv_store_packet_track(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int regen_recv_store_packet_move(endgame_data_t **recvPts, int *numRecvPts, trackingStats *trackCount, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *OTHER, FILE *OTHER2, int *rV, void const *ED_d, void const *ED_mp, int (*change_prec)(void const *, int));
int regen_useSharpener(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp);
int regen_useSharpener_move(int retVal, int sharpenDigits, void const *ED_d, void const *ED_mp);
int regen_recv_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));
int regen_recv_linear_track_send_packet(int headnode, endgame_data_t **startPts, int *numStartPts, endgame_data_t **endPts, int *numEndPts, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));
int regen_worker_linear_track_path(endgame_data_t *startPt, endgame_data_t *endPt, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *CD_d, void const *CD_mp, int (*eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*useSharpener)(int, int, void const *, void const *));

//////////////////////////////////////////////

// functions to send the MPI datatypes - parallel_send_functions.c
void bcast_mat_d(mat_d A, int my_id, int headnode);
void bcast_mat_mp(mat_mp A, int my_id, int headnode);
void bcast_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int my_id, int headnode);
void bcast_vec_d(vec_d b, int my_id, int headnode);
void bcast_vec_mp(vec_mp b, int my_id, int headnode);
void bcast_vec_rat(mpq_t ***b, int size, int my_id, int headnode);
void bcast_comp_d(comp_d c, int my_id, int headnode);
void bcast_comp_num_d(comp_d *c, int num, int my_id, int headnode);
void bcast_comp_mp(comp_mp c, int my_id, int headnode);
void bcast_comp_num_mp(comp_mp *c, int num, int my_id, int headnode);
void bcast_comp_rat(mpq_t c[2], int my_id, int headnode);
void bcast_comp_num_rat(mpq_t c[][2], int num, int my_id, int headnode);
void bcast_worker_info(worker_info *WI, int my_id, int headnode);
void bcast_tracker_config_t(tracker_config_t *T, int my_id, int headnode);
void bcast_basic_eval_data_d(basic_eval_data_d *BED, int MPType, int my_id, int headnode);
void bcast_basic_eval_data_mp(basic_eval_data_mp *BED, int setupProg, int my_id, int headnode);
void bcast_basic_eval_data_amp(basic_eval_data_d *BED, int my_id, int headnode);
void bcast_eqData_t(eqData_t *EqD, int MPType, int my_id, int headnode);
void bcast_witnessData(eqData_t *EqD, int stage, int MPType, int my_id, int headnode);
void bcast_stageData(eqData_t *EqD, int stage, int MPType, int my_id, int headnode);
void bcast_patch_data_d_amp(basic_eval_data_d *BED, int MPType, int my_id, int headnode);
void bcast_prog_t(prog_t *Prog, int MPType, int my_id, int headnode);
void send_recv_point_data_d(point_data_d **PD, int *numPts, int targetNum, int isSending);
int  send_recv_endgame_data_t(endgame_data_t **EG, int *numPts, int MPType, int targetNum, int isSending);
int  send_recv_corank_data(int *corank, double *sm, double *lg, int numPts, int targetNum, int isSending);
int  send_recv_last_approximations(int *prec, point_d *pt_d, point_mp *pt_mp, int *corank, double *sm, double *lg, int numPts, int MPType, int targetNum, int isSending);
int  send_recv_trackBack_samples_t(trackBack_samples_t **TB, int *numPts, int MPType, int targetNum, int isSending);

void cp_witnessCodim_t_int(void *Out, void *In, comp_d **witComp, char **witStr, int **expW, int **witTypes, int **mult, endpoint_data_d_int **wit_d, endpoint_data_mp_int **wit_mp, endpoint_data_amp_int **wit_amp, int curr_prec, int freeData, int MPType, int inType);
void cp_witness_t_int(void *Out, void *In, int **progInst, int **gpSizes, char **progStr, int **PPDtype, int **PPDsize, int **origDegs, int **newDegs, int **P, comp_d **witComp, char **witStr, int freeData, int MPType, int inType);

/////////////////////////// REGEN POS DIM SECTION ///////////////////////////////

typedef struct
{
  prog_t_int Prog_int;      // the slp used to evaluate the original system f
  preproc_data_int PPD_int; // preprocessing data

  int system_rank;       // rank of the original function f
  int orig_variables;    // number of variables used in 'Prog'
  int new_variables;     // number of variables after slicing

  int C_rows;
  int C_cols;
  int H_size;
  int patchCoeff_size;

  int num_funcs;         // number of functions
  int num_codim;         // number of codimensions
  int curr_precision;    // current precision

  int sameA;             // if A is the same for the codims

  int num_int;           
  int num_comp_d;
  int totalLength;
} regen_pos_dim_t_int;

typedef struct
{
  int codim;                         // current codimension

  int useIntrinsicSlice;
  int B_rows;
  int B_cols;
  int p_size;

  int num_comp_d;
  int totalLength;
} regenCodim_t_int;

#ifdef _HAVE_MPI

void create_regen_pos_dim_t_int(MPI_Datatype *mpi_regen_pos_dim_t_int);
void create_regenCodim_t_int(MPI_Datatype *mpi_regenCodim_t_int);

#endif

void cp_regen_pos_dim_int(void *Out, void *In, int MPType, char **rpdStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType);
void cp_regenCodim_int(void *Out, void *In, int MPType, char **rpdStr, comp_d **coeff_d, int curr_prec, int freeStr, int inType);

void bcast_regen_pos_dim_t(regen_pos_dim_t *RPD, int MPType, int my_id, int headnode);
void bcast_regenCodim_t(regenCodim_t *RCD, int curr_prec, int MPType, int my_id, int headnode);

// regen_pos_dim_parallel.c
void regen_pos_dim_par_track(int startCodimIndex, int maxCodim, trackingStats *trackCount, int pathMod, double midpoint_tol, tracker_config_t *T, regen_pos_dim_t *RPD, char *startName, char *witName, int my_id, int num_processes, int headnode);
void worker_regen_pos_dim(int my_id, int num_processes, int headnode, int dataType);
int head_regen_pos_dim_PrepareNextCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *nextStartFile, int my_id, int headnode, int num_processes);
void head_regen_pos_dim_TrackCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void head_regen_pos_dim_SortCodim(int minPacketSize, int maxPacketSize, trackingStats *trackCount, int codim_index, int num_codim, int pathMod, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *START, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *WITSUPER, int my_id, int headnode, int num_processes);
void worker_regen_pos_dim_PrepareNextCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void worker_regen_pos_dim_TrackCodim(trackingStats *trackCount, int curr_level, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);
void worker_regen_pos_dim_SortCodim(trackingStats *trackCount, int codim_index, tracker_config_t *T, regen_pos_dim_t *RPD, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *FAIL, int my_id, int headnode, int num_processes);


	
#endif



