// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

/***********************************************\ 
* bertini.h is the main library for Bertini.    * 
\***********************************************/

#ifndef _BERTINI_H
#define _BERTINI_H

#define BERTINI_VERSION_STRING "1.5.1"
#define BERTINI_DATE_STRING "August 29, 2016"
#define BERTINI_QUIT_MESSAGE "Bertini will now exit due to this error.\n"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>
#include <mpf2mpfr.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _HAVE_MPI
#define bclock_t double
#define bclock(_t) *(_t) = MPI_Wtime()
#define totalTime(_t, _s, _f) *(_t) = _f - _s
#else
#define bclock_t clock_t
#define bclock(_t) *(_t) = clock()
#define totalTime(_t, _s, _f) *(_t) = (_f - _s)/((double) CLOCKS_PER_SEC)
#endif

#define MAX(_a,_b) ((_a)>(_b) ? (_a) : (_b))  
#define MIN(_a,_b) ((_a)<(_b) ? (_a) : (_b))

#define retVal_max_prec_reached -100

// these 2 mean the same thing!!!
#define retVal_reached_minTrackT -50
#define retVal_EG_failed_to_converge -50

#define retVal_cycle_num_too_high -200
#define retVal_PSEG_failed -15
#define retVal_going_to_infinity -2
#define retVal_security_max -4
#define retVal_step_size_too_small -3
#define retVal_too_many_steps -10
#define retVal_refining_failed -20
#define retVal_higher_prec_needed 100
#define retVal_NAN -99
#define retVal_Bertini_Junk -98
#define retVal_Failed_to_converge -97 // this is used in Newton iterations

#define retVal_sharpening_singular_endpoint -22 // this is used when the sharpening sees that the endpoint is singular and cannot sharpen it
#define retVal_sharpening_failed -21 // this is used when the sharpening of an endpoint does not reach the desired tolerance

#define retVal_higher_dim -23 // this is used in regeneration when an endpoint lies on a higher dimensional component

// Exit Error Codes
#define ERROR_OTHER 1
#define ERROR_WRITE_PRIVILEGE 2   // write privilege problems - e.g. not able to create and write to files
#define ERROR_FILE_NOT_EXIST 3    // file exist problems - e.g. expected files do not exist
#define ERROR_INVALID_SIZE 4      // size problems - e.g. expected size not the same as the current size
#define ERROR_MEMORY_ALLOCATION 5 // memory problems - e.g. unable to allocate memory
#define ERROR_CONFIGURATION 6     // configuration problems - e.g. function calls called with wrong input values
#define ERROR_INPUT_SYSTEM 7      // input errors - e.g. degenerate system
#define ERROR_INPUT_SYNTAX 8      // syntax errors for input system
#define ERROR_LOOP_FAILURE 9      // loop failed to exit properly - e.g. fail-safe for 'infinite loops'

typedef struct
{
  int num_funcs;
  int num_hom_var_gp;
  int num_var_gp;
  int *type; // 0 - hom_var_gp, 1 - var_gp
  int *size; // size of the group of the user listed variables (total size = size + type)
} preproc_data;

/*** low-level data types. ***/
typedef struct 
{
  double r, i;
} _comp_d;

typedef struct 
{
  mpf_t r, i;
} _comp_mp;

typedef struct 
{
  _comp_d *coord;
  int alloc_size; // allocated size
  int size;       // size of the point
} _point_d;

typedef struct 
{
  _comp_mp *coord;
  int alloc_size; // allocated size
  int curr_prec;  // current precision
  int size;       // size of the point
} _point_mp;

typedef struct 
{
  _comp_d **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _mat_d;

typedef struct 
{
  _comp_mp **entry;
  int alloc_rows; // allocated number of rows
  int alloc_cols; // allocated number of cols
  int curr_prec;  // current precision
  int rows;       // number of rows for the matrix
  int cols;       // number of cols for the matrix
} _mat_mp;


/* 
   These are the most basic types to be used. 
   `_d'  ==> double precision.
   `_mp' ==> multi-precision.
   `_t'  ==> both/either.
*/
typedef _comp_d   comp_d[1];   /* complex number */
typedef _point_d  point_d[1];  /* complex point */
typedef  point_d  vec_d;       /* complex vector (same structure as a point - useful to have two names) */
typedef _mat_d    mat_d[1];    /* complex matrix */

typedef _comp_mp  comp_mp[1];
typedef _point_mp point_mp[1];
typedef  point_mp vec_mp; 
typedef _mat_mp   mat_mp[1];

// these contain all of the data structures that are returned when doing function evaluation
typedef struct
{
  point_d funcVals;
  point_d parVals;
  vec_d parDer;
  mat_d Jv;
  mat_d Jp;
} eval_struct_d;

typedef struct
{
  point_mp funcVals;
  point_mp parVals;
  vec_mp parDer;
  mat_mp Jv;
  mat_mp Jp;
} eval_struct_mp;

/* A real number - rat holds a rational representation, real holds the real approximation (to currPrec, the current precision).
This data structure is not used much (only in function evaluation. */
typedef struct  
{		
  mpq_t rat;
  mpf_t real;
  int currPrec;
} num_t;

/*** The straight-line program structure.  This is the way that polynomials are stored internally. ***/
typedef struct {
  int *prog;     /*  The program instructions. (a big integer array)  */
  int  size;     /*  size of the instruction program.   */   
  int  memSize;  /* Amount of memory it needs in workspace (for temp and final results).*/   
  num_t *nums;   /* The array of real numbers. */
  int precision; /* The precision at which evaluation should occur */
  
  /* INFO NEEDED FOR M-HOM: */
  int num_var_gps;  /* The total number of variable groups (i.e., m from m-hom).*/
  int *var_gp_sizes;  /* The size of each of the groups. */
  int index_of_first_number_for_proj_trans;  /* The address of the first number used in the projective transformation polynomials. */
                                                                                    
  /* STOP LOCATIONS: */
  int  numInstAtEndUpdate; /* instruction number at end of update. i.e. i = 0; while (i < numInstAtEndUpdate) .. */
  int  numInstAtEndParams; /* instruction number at end of params. i.e. i = numInstAtEndUpdate; while (i < numInstAtEndParams) .. */
  int  numInstAtEndFnEval; /* instruction number at end of function eval. i.e. i = numInstAtEndParams; while (i < numInstAtEndFnEval) .. */
  int  numInstAtEndPDeriv; /* instruction number at end of param diff. i.e. i = numInstAtEndFnEval; while (i < numInstAtEndPDeriv) .. */
  int  numInstAtEndJvEval; /* instruction number at end of Jv eval. i.e. i = numInstAtEndPDeriv; while (i < numInstAtEndJvEval) .. */
                           /* for Jp eval: i = numInstAtEndJvEval; while (i < size) .. */ 

  /* INPUT AMOUNTS: */
  int  numVars;  /*  Number of variables in the function being computed.   */
  int  numPathVars;  /*  Number of path variables.  Ought to be 1 usually.   */
  int  numNums;  /*  Number of real numbers used in evaluation.  */
  int  numConsts;  /*  Number of constants.  */   
                                                                                   
  /* OUTPUT AMOUNTS: */
  int  numPars;  /*  Number of parameters   */
  int  numFuncs; /*  Number of coordinate functions in the homotopy.   */
  int  numSubfuncs;  /*  Number of subfunctions.  */
                                                                                      
  /* INPUT LOCATIONS: */
  int  inpVars;  /*  Where the input variable values are stored.   */
  int  inpPathVars;  /*  Where the values of the path variables are stored.   */
  int  IAddr;  /*  Where the constant I is stored.  */
  int  numAddr;  /*  Where the first num_t is stored.  */
  int  constAddr;  /*  Where the first constant is stored.  */
                                                                                      
  /* OUTPUT LOCATIONS: */
  int  evalPars;  /*  Where U(t), for given t, is stored.   */
  int  evalDPars;  /*  Where the derivatives of the parameters are stored.   */
  int  evalFuncs;  /*  Where H(x,t) is stored.   */
  int  evalJVars;  /*  Where the Jacobian w.r.t. vars is stored.   */
  int  evalJPars;  /*  Where the Jacobian w.r.t. pars is stored.   */
  int  evalSubs;  /*  Where the subfunctions are stored */
  int  evalJSubsV;  /*  Where the derivatives of the subfunctions w.r.t. vars are stored.  */
  int  evalJSubsP;  /*  Where the derivatives of the subfunctions w.r.t. pars are stored.  */
} prog_t;

/***  The point_data_t structures hold the data for a given path.  It is basic to the path tracker.  
The "param" field was not being used (as of 6/21/06), so it has been commented out. ***/
typedef struct {
  point_d point; 
//  point_d param;
  comp_d  time;
  int     cycle_num;
} point_data_d;

typedef struct {
  point_mp point;
//  point_mp param;
  comp_mp  time;
  int      cycle_num;
} point_data_mp;         

/*** The midpoint_data_t structures are used only in the post-processor to compare the path values at t=tEndgame. ***/

typedef struct {
  point_d point;
  double norm;
  int path_num;
} midpoint_data_d;

typedef struct {
  point_mp point;
  mpf_t norm;
  int path_num;
} midpoint_data_mp;

typedef struct {
  int prec;
  point_d  Pt_d;
  point_mp Pt_mp;
} midpoint_data_t;

// the post_process_t structure is used in post-processing //
typedef struct 
{
  int path_num;     // path number of the solution
  int sol_num;      // solution number
  comp_d  *sol_d;   // solution
  comp_mp *sol_mp;
  int sol_prec;     // precision of the solution
  int size_sol;     // the number of entries in sol
  double function_resid_d;  // the function residual
  mpf_t  function_resid_mp; 
  double cond_est;  // the estimate of the condition number
  double newton_resid_d;    // the newton residual
  mpf_t  newton_resid_mp;
  double final_t;   // the final value of time
  double accuracy_estimate; // accuracy estimate between extrapolations
  double first_increase;    // time value of the first increase in precision
  int cycle_num;    // cycle number used in extrapolations
  int success;      // success flag 
  int multiplicity; // multiplicity
  int isReal;       // real flag:  0 - not real, 1 - real
  int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
  int isSing;       // singular flag: 0 - non-sigular, 1 - singular
} post_process_t; 

/***  The tracker_config_t structure is key.  It holds the basic bits of data necessary to track, except for the SLP. ***/
/* Not all of this data will be used on a given run, e.g. latest_cond_num_exp will not be used if using fixed precision. */
typedef struct {
  int numVars;
  int numPathVars;  /* This is always 1 -> could be removed!!!??? */
  int numParams;
  int numFuncs;
  double maxStepSize;  /* The maximum step size (in t) allowed by the user. */
  double minStepSizeBeforeEndGame;  /* Threshold for declaring failure due to adaptive step size becoming too small, before endgame. */
  double minStepSizeDuringEndGame;  /* Threshold for declaring failure due to adaptive step size becoming too small, during endgame. */
  double minStepSize;  /* Will generally be set to one of the previous two (eliminates many unattractive if/else's). */ 
  double currentStepSize; /* Useful when exiting and re-entering tracker (always starting tracker w/ maxStepSize may waste time if small step is needed). */
  int first_step_of_path;  /* Indicates whether the first step on the current path has been taken. */
  double minTrackT;    /* Smallest time interval over which tracking may occur - may be obsolete!!!???!!! */
  double basicNewtonTol;  /* The Newton tolerance before the endgame. */
  double endgameNewtonTol;  /* The Newton tolerance during the endgame. */
  double currentNewtonTol;  /* The Newton tolerance currently being used - may be one of the two above or other (if using adaptive tolerances). */
  int cSecInc;      /* # of consecutive successful steps before trying to increase deltaT. */
  int maxNewtonIts; /* Max number of Newton iterations in one step.  Set this low to avoid path-crossing. */
  int MPType;  /* 0 for machine precision, 1 for multiprecision, 2 for adaptiveMP. */
  int Precision;  /* Current precision (not used if MPType = 0). */
  int outputLevel;  /* 0 for just start and end points. */
                    /* 1 = 0 + all times. */
                    /* 2 = 1 + all residuals of Newton. */
                    /* 3 = 2 + all function evals, etc.  (EVERYTHING!) */
  int    screenOut;  /* 0 for no output to screen, 1 for all output to screen (as well as file). */
  double targetT;  /* Final desired value of path variable (not always 0, could be endgameBoundary (below) or some complex number for tracking around 0. */
  double endgameBoundary;  /* Path variable value at which switch to endgame occurs. */
  int endgameSwitch;  /* Set to 0 if tracking before endgameBoundary, 1 if in endgame. Useful in tracker. */
  double goingToInfinity;  /* If Newton residual exceeds this, the path is declared a failure (going to infinity). */
  int maxNumSteps;  /* Maximum number of steps to be taken during the endgame. */
  int endgameNumber;  /* Indicates which endgame to use. See config for options. */
  int latest_cond_num_exp;  /* Latest estimate of the exponent of the condition number (for AMP). */
  int steps_since_last_CN;  /* Number of Euler steps since last approximation of the condition number (for AMP). */
  double power_series_sample_factor;  /* The factor used to choose values of t at which to sample for the power series endgame. */
  int cycle_num_max;  /* The largest cycle number to be considered in the power series endgame. */
  int num_PSEG_sample_points; // The number of sample points used in the interpolation for PSEG
  double latest_newton_residual_d; // Keeps track of the most recent newton residual
  mpf_t  latest_newton_residual_mp;
  double t_val_at_latest_sample_point;  /* Keeps track of the most recent value of t in the power series. */
  double error_at_latest_sample_point;  /* Keeps track of the error between consecutive extrapolations to the origin in the PS EG. */
  double final_tolerance;  /* Residual at which we are satisfied with the approximation by the power series endgame. */
  double real_threshold;  /* Threshold at which point we call imaginary part of complex number 0 (i.e., we call it a real number if im < real_threshold). */
  int endgameOnly;  /* Set internally (not via config) to 1 if only the endgame should be used, 0 if basic tracking should be used up to endgame breakpoing. */
  double AMP_bound_on_abs_vals_of_coeffs;  /* User-defined bound on the sum of the abs vals of the coeffs for any polynomial in the system (for adaptive precision). */
  double AMP_bound_on_degree;  /* User-set bound on degrees of polynomials in the system - tricky to compute for factored polys, subfuncs, etc. (for adaptive precision). */
  double AMP_eps;  /* Bound on \epsilon (an error bound) from adaptive precision paper - see paper for details. */
  double AMP_Phi;  /* Bound on \Phi (an error bound) from adaptive precision paper - see paper for details. */
  double AMP_Psi;  /* Bound on \Psi (an error bound) from adaptive precision paper - see paper for details. */

  int AMP_safety_digits_1;
  int AMP_safety_digits_2;
  int AMP_max_prec;

  double sing_val_zero_tol;  /* User-defined (eventually) setting for threshold between numerically zero and nonzero singular values. */ 
  double cond_num_threshold;

  double step_fail_factor;
  double step_success_factor;

  int max_num_pts_for_trace;
  int max_num_mon_linears;
  int max_num_bad_loops_in_mon;

  double final_tol_multiplier;
  double final_tol_times_mult;

  // sharpen information
  int sharpenOnly;
  int sharpenDigits; // number of digits to sharpen to - non-zero means that we utilize sharpening

  // regen information
  int regen_remove_inf;        // determine whether to keep tracking or remove the infinite endpoints
  int regen_higher_dim_check;  // determine whether to use the higher dimensional check in regeneration
  double sliceBasicNewtonTol;  // newton tolerance before the endgame for moving slices in regeneration
  double sliceEndgameNewtonTol;// newton tolerance during the endgame for moving slices in regeneration
  double sliceFinalTol;        // final tolerance for moving slices in regeneration

  // other information
  int minCycleTrackBack;       // minimum cycle for using track back endgame
  int junkRemovalTest;         // use either LDT or membership test
  int maxDepthLDT;             // maximum depth to use LDT before switching to membership test
  int odePredictor;            // which ODE predictor to use

  // security level
  int securityLevel;           // security level for path tracking
  double securityMaxNorm;      // max norm for the security level

  // endgame cutoff values
  double cutoffCycleTime;      // time cutoff for cycle agreement for 'pre-Cauchy endgame'
  double cutoffRatioTime;      // time cutoff for ratio agreement

  // finite threshold
  double finiteThreshold;      // maximum norm to be considered finite point

  // functtion residual tolerance
  double funcResTol;           // tolerance for function evaluation

  // ratio tolerance
  double ratioTol;             // tolerance for ratios when deciding zero

  // AMP3 tracking 
  int maxStepsBeforeNewton;   // maximum number of steps before using a Newton iteration when using AMP3 tracking

} tracker_config_t;

// structure to hold the counts for the tracking statistics
typedef struct {
  int numPoints;     // number of points that are tracked
  int successes;     // number of successes
  int failures;      // number of failures
  int junkCount;     // number of Bertini junk
  int nanCount;      // number that are not a number
  int infCount;      // number that exceeded the max norm
  int securityCount; // number that security truncated
  int sizeCount;     // number that had step size drop below min
  int PSEGCount;     // number that failed the Power Series EndGame
  int precCount;     // number that exceeded the max precision
  int cycleCount;    // number that had a cycle number that was too high
  int stepCount;     // number that took too many steps
  int refineCount;   // number that failed to refine
  int otherCount;    // all the other failures
} trackingStats;

// structure to hold PSEG sample points
typedef struct {
  int num_samples;
  int mem_count;
  point_data_d *samples;
  vec_d *dX;
  _point_d *Z_rev;
} PSEG_samples_struct_d;

typedef struct {
  int curr_prec;
  int num_samples;
  int mem_count;
  point_data_mp *samples;
  vec_mp *dX;
  _point_mp *Z_rev;
} PSEG_samples_struct_mp;

typedef struct {
  int *max_prec; // maximum precision that we know the sample point to
  int num_samples; // number of samples actually found
  int mem_count; // number of memory locations
  point_data_d  *samples_d; // samples in double precision
  point_data_mp *samples_mp; // samples in multi precision
  vec_d  *dX_d;
  vec_mp *dX_mp;
  _point_d  *Z_rev_d;
  _point_mp *Z_rev_mp;
} PSEG_samples_struct_amp;

typedef struct
{
  int prec;
  point_data_d PD_d;
  point_data_mp PD_mp;

  int last_approx_prec;       // precision of the last approximation
  point_d last_approx_d;      // last approximation to the end point
  point_mp last_approx_mp;    // last approximation to the end point

  int retVal;
  int pathNum;
  int codim;
  double first_increase;
  double condition_number;
  double function_residual_d;
  mpf_t  function_residual_mp;
  double latest_newton_residual_d;
  mpf_t  latest_newton_residual_mp;
  double t_val_at_latest_sample_point_d;
  mpf_t  t_val_at_latest_sample_point_mp;
  double error_at_latest_sample_point_d;
  mpf_t  error_at_latest_sample_point_mp;
} endgame_data_t;

typedef struct
{
  endgame_data_t endPt;

  int numSamples;
  int samplePts_prec;
  double *normSamples_d;
  mpf_t *normSamples_mp;
  point_d *samplePts_d;
  point_mp *samplePts_mp;
  int *samplePts_notUsed;
} endgame_samples_t; 

typedef struct
{
  endgame_data_t endPt;

  int numSamples;
  int samplePts_prec;
  double *normSamples_d;
  mpf_t *normSamples_mp;
  point_d *samplePts_d;
  point_mp *samplePts_mp;
  int *samplePts_notUsed;

  int midPt_prec;
  point_d *midPt_d;
  point_mp *midPt_mp;

} trackBack_samples_t;

/*** eval_data structures and eval_t prototypes ***/
/* Each form of evaluation requires special sets of data, so each form of evaluation gets its own structure. */

/* This is the basic, eval_data structure.  It is assumed that regular SLP evaluation is all that is needed, so we only include one SLP. */
/* The cascade structures and prototypes may be found in cascade.h. */
typedef struct {
prog_t *Prog;
}  basic_eval_data;

/*** Most prototypes. ***/
/* Some prototypes for minor functions appear in the file in which they are needed. */

/* It's VERY IMPORTANT to call the following before you use any of the _mp functions or macros!!! */
void  initMP(int prec);  /* Initializes all global instances of MP data types and sets precision to prec. */ 
void  clearMP();

double d_abs_mp(comp_mp a);  /* Absolute value of a MP complex number - returned as a double. */
void   pow_rdouble_d(comp_d res, comp_d a, double e); /* Store a^e in res (double power of a complex number. */
void   pow_rdouble_mp(comp_mp res, comp_mp a, double e);  /* same as above */
void   pow_rmpf_mp(comp_mp res, comp_mp base, mpf_t e); // a^e where e is multiprecision

double d_vec_abs_d(vec_d v); /* One-norm of the vector v!!! */
double d_vec_abs_mp(vec_mp v); /* same as above. */
void   print_d(FILE *fp, int digits, comp_d Z); /* Print a complex number to file pointed at by fp */
void   printMat_d(FILE *fp, int digits, mat_d M); /* Print matrix to file */
void   printVec_d(FILE *fp, int digits, vec_d v); /* Print vector to file */
void   printPoint_d(FILE *fp, int digits, point_d v); /* Print point to file */

void print_Matlab_d(FILE *fp, int digits, comp_d Z);
void printMat_Matlab_d(FILE *fp, int digits, mat_d M);
void printVec_Matlab_d(FILE *fp, int digits, vec_d v);
void printPoint_Matlab_d(FILE *fp, int digits, point_d v);

void   print_rat(FILE *fp, mpq_t *Z);

void   print_mp(FILE *fp, int digits, comp_mp Z);
void   printMat_mp(FILE *fp, int digits, mat_mp M);
void   printVec_mp(FILE *fp, int digits, vec_mp v);
void   printPoint_mp(FILE *fp, int digits, point_mp v);

void print_Matlab_mp(FILE *fp, int digits, comp_mp Z);
void printMat_Matlab_mp(FILE *fp, int digits, mat_mp M);
void printVec_Matlab_mp(FILE *fp, int digits, vec_mp v);
void printPoint_Matlab_mp(FILE *fp, int digits, point_mp v);

void mat_perp_d(mat_d Cp, mat_d C);
void mat_perp_mp(mat_mp Cp, mat_mp C, int prec);
void mat_perp_rat(mat_d Cp_d, mat_mp Cp_mp, mpq_t ***Cp_rat, mpq_t ***C_rat, int rows, int cols, int curr_prec, int max_prec, int init_rat);

void print_comp_out_d(FILE *FP, comp_d A);
void print_vec_out_d(FILE *FP, vec_d A);
void print_mat_out_d(FILE *FP, mat_d A);
void setup_comp_in_d(comp_d A, FILE *FP);
void setup_vec_in_d(vec_d A, FILE *FP);
void setup_mat_in_d(mat_d A, FILE *FP);
void print_comp_out_mp(FILE *FP, comp_mp A);
void print_vec_out_mp(FILE *FP, vec_mp A);
void print_mat_out_mp(FILE *FP, mat_mp A);
void setup_comp_in_mp(comp_mp A, FILE *FP);
void setup_vec_in_mp(vec_mp A, FILE *FP);
void setup_mat_in_mp(mat_mp A, FILE *FP);
void print_comp_out_rat(FILE *FP, mpq_t *A);
void print_vec_out_rat(FILE *FP, mpq_t **A, int size);
void print_mat_out_rat(FILE *FP, mpq_t ***A, int rows, int cols);
void setup_comp_in_rat(mpq_t *A, FILE *FP);
void setup_vec_in_rat(mpq_t ***A, FILE *FP, int *size);
void setup_mat_in_rat(mpq_t ****A, FILE *FP, int *rows, int *cols);

/* Prototypes from zero_dim_main.c: */
int  setupProg(prog_t *P, int precision, int MPType); /* Reads straight-line program instructions from arr.out and numbers from num.out into P. */
int setupProg_count(prog_t *P, int precision, int MPType, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow);
void setupNums(num_t **nums, int numNums, int precision, int MPType);
void clearProg(prog_t *P, int MPType, int clearEvalProg);  // clears P
void clearNums(num_t **nums, int numNums); // clear nums
void setupStart_d(tracker_config_t *T, point_data_d *PD, FILE *StartPts);  /* Grabs a start point (for a path) from file "start" */
void setupStart_mp(tracker_config_t *T, point_data_mp *PD, FILE *StartPts);  /* same as above */
void zero_dim_main(int MPType, double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);
int  zero_dim_main_d(int MPType, double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);
int  zero_dim_main_mp(double parse_time, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);
void remove_output_files(int trackType, int sharpenOnly, int removeRawData);
void nonsolutions_check_compare_d(int *isZero, point_d f1, point_d f2, int startFunc, int endFunc, double zero_threshold, double max_ratio);
int nonsolutions_check_d(int size_f, int size_r, point_d x, point_d last_x, comp_d time, double zero_threshold, double max_ratio, prog_t *Prog);
void nonsolutions_check_compare_mp(int *isZero, point_mp f1, point_mp f2, int startFunc, int endFunc, double zero_threshold, double max_ratio);
int nonsolutions_check_mp(int size_f, int size_r, point_mp x, point_mp last_x, comp_mp time, double zero_threshold, double max_ratio, prog_t *Prog);
int bertini_junk_check_d(int size_f, int size_r, point_d x, comp_d time, double zero_threshold, prog_t *Prog);
int bertini_junk_check_mp(int size_f, int size_r, point_mp x, comp_mp time, double zero_threshold, prog_t *Prog);
void remove_temp_files();
void getDehomPoint_comp_d(point_d dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_d *sol, int num_vars, preproc_data *PPD, double accuracyEstimate);
void getDehomPoint_comp_mp(point_mp dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_mp *sol, int num_vars, preproc_data *PPD, double accuracyEstimate);
void initialize_mpi(int argc, char *args[], int *num_of_processes, int *my_id);
void finalize_mpi(void);

/* Prototypes from parse_input.c: */
void parse_input(char *inputName, int *trackType, int *MPType, int *genType, int *userHom, unsigned int *randomSeed, int *sharpenOnly, int *needToDiff, int *remove_temp, int useParallelDiff, int my_id, int num_processes, int headnode);
int setupConfig(tracker_config_t *T, double *midpointTol, int *userHom, int *useRegen, int *regenStartLevel, int *maxCodim, int *specificCodim, int *printMod, double *intrinsicCutoffMultiplier, int *reducedOnly, int *constructWitnessSet, int *supersetOnly, int *paramHom, int MPType); // Opens the file "config" and stores data in tracker_config_t.
void printConfigValues(FILE *OUT, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom);
int setupRandomValues(FILE *OUT, FILE *IN, int createNewValues, int maxPrec);

/* Prototypes from matrixSolve.c: */
int matrixSolve_d(vec_d x, mat_d A, vec_d b); // solves A*x = b
int matrixSolve_cond_num_norms_d(vec_d x, mat_d A, vec_d b, double *cond_num, double *norm_A, double *norm_A_inv); // solves A*x = b & finds the
                                                                             // condition number of A along with the norm of A and A^-1

int matrixSolve_LU_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange); // solves A*x = b for x using LU decompostiion on a square matrix A
int matrixSolve2_LU_d(vec_d x1, vec_d x2, mat_d A, vec_d b1, vec_d b2, double tol, double largeChange); // solves A*x1 = b1 & A*x2 = b2 for x1 & x2 simultaneously

int matrixSolve_LU_cond_num_norms_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); // solves A*x = b for x using LU decomp 
int matrixSolve_svd_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); // finds the least squares solution
// to A*x = b using the svd factorization of A & finds the condition number of A along with the norm of A and A^-1
int matrixSolve_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); // finds the least squares solution 
// to A*x = b using the QR factorization of A & finds the condition number of A along with the norm of A and A^-1

int LU_matrixSolve_d(vec_d x, mat_d LU, int **rwnm, int *sign, mat_d A, vec_d b, double tol, double largeChange); 
int LU_matrixSolve_mp(vec_mp x, mat_mp LU, int **rwnm, int *sign, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange);
int matrixSolve_from_LU_d(mat_d X, mat_d LU, int *rownum, mat_d B); 
int matrixSolve_from_LU_mp(mat_mp X, mat_mp LU, int *rownum, mat_mp B);

int matrixSolve_mp(vec_mp x, mat_mp A, vec_mp b); // solves A*x = b for x
int matrixSolve_cond_num_norms_mp(vec_mp x, mat_mp A, vec_mp b, double *cond_num, double *norm_A, double *norm_A_inv); // solves A*x = b & finds the condition number of A along with the norm of A and A^-1

int matrixSolve_LU_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange); // solves A*x = b for x using LU decompostiion on a square matrix A
int matrixSolve2_LU_mp(vec_mp x1, vec_mp x2, mat_mp A, vec_mp b1, vec_mp b2, double tol, double largeChange); // solves A*x1 = b1 & A*x2 = b2 for x1 & x2 simultaneously
int matrixSolve_LU_cond_num_norms_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); // solves A*x = b for x using LU decomp

// finds the least squares solution to A*x = b using the svd factorization of A & finds the condition number of A along with the norm of A and A^-1
int matrixSolve_svd_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); 
int matrixSolve_svd_Least_Squares_mp2(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange, double *cond_num, double *norm_A, double *norm_A_inv);

// finds the least squares solution to A*x = b using the QR factorization of A & finds the condition number of A along with the norm of A and A^-1
int matrixSolve_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv); 
int matrixSolve_Least_Squares_mp2(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange, double *cond_num, double *norm_A, double *norm_A_inv);


int matrixSolve_Hessenberg_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange);
void gengr_d(comp_d f, comp_d g, comp_d top, mat_d A, int r, int c, double tol_sign);
int matrixSolve_Hessenberg_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange);
void gengr_mp(comp_mp f, comp_mp g, comp_mp top, mat_mp A, int r, int c, mpf_t tol_sign);

/* Prototypes from svd_from_QR.c */
int svd_corank_d(mat_d A, double rank_tol, double largeChange);
int svd_corank_analyze_d(mat_d E, double rank_tol, double largeChange);
int svd_jacobi_d_prec(mat_d U, mat_d E, mat_d V, mat_d InputMat, double rankTol);
int svd_jacobi_d(mat_d U, mat_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange);
int svd_R_jacobi_d(mat_d U, mat_d E, mat_d V, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange);
int jacobi_d(mat_d U, mat_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange);
int checkGood_d(double currSV, double prevSV, double tooSmall, double largeChange);
// the following 2 functions are for finding only the singular values
int svd_jacobi_E_d(vec_d E, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange);
int jacobi_E_d(vec_d E, mat_d U, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange);
// the following 3 functions are for finding only singular values and right singular vectors
int svd_jacobi_EV_d(vec_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange);
int svd_R_jacobi_UE_d(mat_d U, vec_d E, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange);
int jacobi_UE_d(mat_d U, vec_d E, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange);

int svd_corank_mp(mat_mp A, mpf_t rank_tol, mpf_t largeChange);
int svd_corank_analyze_mp(mat_mp E, mpf_t rank_tol, mpf_t largeChange);
int svd_jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange);
int svd_jacobi_mp_prec(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, double rankTol, int curr_prec);
int svd_R_jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange);
int jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange);
int checkGood_mp(mpf_t currSV, mpf_t prevSV, mpf_t tooSmall, mpf_t largeChange);
// the following 2 functions are for finding only the singular values
int svd_jacobi_E_mp(vec_mp E, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange);
int jacobi_E_mp(vec_mp E, mat_mp U, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange);
// the following 3 functions are for finding only singular values and right singular vectors
int svd_jacobi_EV_mp(vec_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange);
int svd_R_jacobi_UE_mp(mat_mp U, vec_mp E, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange);
int jacobi_UE_mp(mat_mp U, vec_mp E, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange);

int corank_rrv_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0, mat_d mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio, double FV_tol);
int corank_rrv_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, mat_mp mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio, double FV_tol);

/* Prototypes from QR.c */
int QR_d_prec(mat_d Q, mat_d R, mat_d P, mat_d A, int preSortRows);
int QR_d(mat_d Q, mat_d R, mat_d P, mat_d A, double tol_pivot, double tol_sign, double largeChange, int preSortRows);
void convertToP_d(mat_d P, int *perm, int n);
void convertToPerm_d(int *perm, mat_d P, int n);
void genhh_d(vec_d u, mat_d A, int sr, int c, double norm, int *colnum, double tol);
void findZ_mat_d(vec_d z, mat_d C, vec_d u, int sr, int sc, int *colnum);
void findZ_vec_d(comp_d z, vec_d C, vec_d u, int sr);
void apphh_Z_mat_d(mat_d C, vec_d u, vec_d z, int sr, int er, int sc, int ec, int u_sc, int u_ec, int z_sc, int z_ec, int *colnum);
void apphh_Z_vec_d(vec_d C, vec_d u, comp_d z, int sr, int er, int u_sc, int u_ec);
double sign_d(comp_d s, comp_d in, double tol);
double cond_est_d(vec_d y, vec_d x, double norm_x, mat_d C, int sr, int er, int col, comp_d gamma, int *colnum, double tol);
// do QR factorization without finding Q
int QR_R_d(mat_d R, mat_d P, mat_d A, double tol_pivot, double tol_sign, double largeChange, int preSortRows);
int matrixSolve_from_QR_d(vec_d x, mat_d Q_trans, mat_d R, int *perm, vec_d b);

int QR_mp_prec(mat_mp Q, mat_mp R, mat_mp P, mat_mp A, int preSortRows, int curr_prec);
int QR_mp(mat_mp Q, mat_mp R, mat_mp P, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange, int preSortRows);
void convertToP_mp(mat_mp P, int *perm, int n);
void convertToPerm_mp(int *perm, mat_mp P, int n);
void genhh_mp(vec_mp u, mat_mp A, int sr, int c, mpf_t norm, int *colnum, mpf_t tol);
void findZ_mat_mp(vec_mp z, mat_mp C, vec_mp u, int sr, int sc, int *colnum);
void findZ_vec_mp(comp_mp z, vec_mp C, vec_mp u, int sr);
void apphh_Z_mat_mp(mat_mp C, vec_mp u, vec_mp z, int sr, int er, int sc, int ec, int u_sc, int u_ec, int z_sc, int z_ec, int *colnum);
void apphh_Z_vec_mp(vec_mp C, vec_mp u, comp_mp z, int sr, int er, int u_sc, int u_ec);
double sign_mp(comp_mp s, comp_mp in, double tol);
double sign_mp2(comp_mp s, comp_mp in, mpf_t tol);
void sign_mp3(mpf_t norm, comp_mp s, comp_mp in, mpf_t tol);
void cond_est_mp(vec_mp y, vec_mp x, mpf_t norm_x_y, mat_mp C, int sr, int er, int col, comp_mp gamma, int *colnum, mpf_t tol);
// do QR factorization without finding Q
int QR_R_mp_prec(mat_mp R, mat_mp P, mat_mp A, int preSortRows, int curr_prec);
int QR_R_mp(mat_mp R, mat_mp P, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange, int preSortRows);
int matrixSolve_from_QR_mp(vec_mp x, mat_mp Q_trans, mat_mp R, int *perm, vec_mp b);

/* Prototypes from QLP.c: */
int QLP_L_d(mat_d L, mat_d A, double tol_pivot, double tol_sign, double largeChange);
int QLP_L_mp_prec(mat_mp L, mat_mp A, int curr_prec);
int QLP_L_mp(mat_mp L, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange);

/* Prototypes from misc.c: */
void init_trackingStats(trackingStats *S);
void add_trackingStats(trackingStats *Total, trackingStats *Input, int num_input);
void reproduceInputFile(char *outputName, char *inName, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom);
double maxDiffMat_d(mat_d A, mat_d B);
double infNormVec_d(vec_d X); /* Infinity norm of a vector. */
double frobNormMat_d(mat_d M);
double infNormMatRow_d(mat_d A, int i); /* Infinity norm of the row or a matrix. */
double infNormMat_d(mat_d M); /* Infinity norm of a matrix. */    
double conditionNumber_d(mat_d M); /* ESTIMATE of the condition number of a matrix. */
double conditionNumber_and_norms_d(mat_d M, double *norm_M, double *norm_M_inv);  /* ESTIMATE of the condition number of a matrix AND the norms of M and M inverse. */
void mat_d_to_mp_rat(mat_mp MP, mpq_t ***RAT, mat_d D, int digits, int prec, int need_to_init_mp, int need_to_init_rat); // _d matrix to _mp & rational
void comp_d_to_mp_rat(comp_mp MP, mpq_t *RAT, comp_d D, int digits, int prec, int need_to_init_mp, int need_to_init_rat); // comp_d to rational 'num' & mp MP
void d_to_mp2(comp_mp m, comp_d d);
double maxDiffMat_mp(mat_mp A, mat_mp B);
double infNormVec_mp(vec_mp X); /* See above _d version */
void infNormVec_mp2(mpf_t norm, vec_mp X);
double infNormMatRow_mp(mat_mp A, int i); /* See above _d version */
double infNormMat_mp(mat_mp M); /* See above _d version */
void infNormMat_mp2(mpf_t norm, mat_mp M);
double frobNormMat_mp(mat_mp M);
double conditionNumber_mp(mat_mp M); /* See above _d version */
double conditionNumber_and_norms_mp(mat_mp M, double *norm_M, double *norm_M_inv);
void   mul_mat_vec_d(vec_d Res, mat_d M, vec_d X);  /* Multiply a matrix, M, by a vector, X. */
void   mul_mat_vec_mp(vec_mp Res, mat_mp M, vec_mp X); /* same */
void vec_mat_mul_d(vec_d Res, vec_d b, mat_d A);
void vec_mat_mul_mp(vec_mp Res, vec_mp b, mat_mp A);
void vec_mat_mul_rat(mpq_t **Res, mpq_t **b, mpq_t ***A, int b_size, int A_rows, int A_cols, int init_res);
void   add_vec_d(vec_d Res, vec_d x, vec_d y); /* Add two vectors (componentwise). */
void   add_vec_mp(vec_mp Res, vec_mp x, vec_mp y); /* same */
void   hermiteInterpCW_d(point_d Res, _comp_d *T, _point_d *Y, _point_d *dHdT, comp_d T0, int n); /* Interpolate the n function values Y and derivs dHdT */
												  /* at t-values T to get function value Res at T0. */
void   hermiteInterpCW_mp(point_mp Res, _comp_mp *T, _point_mp *Y, _point_mp *dHdT, comp_mp T0, int n); /* same */
int    track_d(point_data_d *Final, point_data_d *Start, comp_d FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); /* Track from Start.time to FinalT, if possible, */
												       /* using T (containing Program), with result in */
												       /* Final.  Report output to file OUT. */
int    track_mp(point_data_mp *Final, point_data_mp *Start, comp_mp FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));  /* same */
void   twoNormVec_d(vec_d V, double *retVal); /* The 2-norm of a vector (note that d_vec_abs_d is the 1-norm). */
void   twoNormVec_mp(vec_mp V, comp_mp retVal); /* same */
void twoNormVec_mp2(vec_mp V, mpf_t retVal);
void   outerProduct_d(mat_d Res, vec_d V, vec_d W); /* Outer product of two vectors. */
void   outerProduct_mp(mat_mp Res, vec_mp V, vec_mp W); /* same */
void   transpose_d(mat_d Res, mat_d M); /* Stores the transpose of M in Res. */
void   transpose_mp(mat_mp Res, mat_mp M); /* same */
void   find_opposite_phase_d(comp_d opposite_phase, comp_d number); /* Finds a complex number, opposite_phase, so that product with "number" is real. */ 
void   find_opposite_phase_mp(comp_mp opposite_phase, comp_mp number); /* same */
void make_matrix_ID_d(mat_d A, int rows, int cols); /* Sets A to identity matrix. */
void make_matrix_ID_mp(mat_mp A, int rows, int cols); 
void make_matrix_ID_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int rows, int cols, int prec, int max_prec, int init_mp, int init_rat);
void make_square_matrix_orth_d(mat_d A, int size); // sets A to a size x size conjugate-orthonormal matrix
void make_square_matrix_orth_mp(mat_mp A, int size, int prec); // sets A to a size x size conjugate-orthonormal matrix
void make_vec_random_d(vec_d p, int size);  // unit vector in inf-norm
void make_vec_random_d2(vec_d p, int size); // unit vector in 2-norm
void make_vec_random_mp(vec_mp p, int size);
void make_vec_random_rat(vec_d p_d, vec_mp p_mp, mpq_t **p_rat, int size, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat);
void make_matrix_random_d(mat_d A, int rows, int cols); /* Sets A to a random matrix. */
void make_matrix_random_mp(mat_mp A, int rows, int cols, int prec); /* Sets A to a random matrix. */
void make_matrix_random_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int rows, int cols, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat);

void make_matrix_random_real_d(mat_d A, int rows, int cols); /* Sets A to a random real matrix. */
void make_matrix_random_real_mp(mat_mp A, int rows, int cols, int prec); /* Sets A to a random real matrix. */
void make_square_matrix_orth_real_d(mat_d A, int size); // sets A to a size x size orthonormal real matrix
void make_square_matrix_orth_real_mp(mat_mp A, int size, int prec); // sets A to a size x size orthonormal real matrix

void   make_elem_d(mat_d A, int num1, int num2);  /* Sets A to identity matrix, but with rows num1 and num2 and cols num1 and num2 swapped. */
void   make_elem_mp(mat_mp A, int num1, int num2); /* same */
void   mat_mul_d(mat_d Res, mat_d A, mat_d B); /* Multiply two matrices. */
void   mat_mul_mp(mat_mp Res, mat_mp A, mat_mp B); /* same */
void mat_mul_rat(mpq_t ***Res, mpq_t ***A, mpq_t ***B, int A_rows, int A_cols, int B_rows, int B_cols, int max_prec);
void   Gram_Schmidt_d(mat_d A);  /* Perform Gram-Schmidt orthogonalization on A. */
int    get_rand_int(int *x);  /* Used to grab a random integer. */
void create_random_number_str(char *str, int num_digits);

double get_comp_rand_d(comp_d x);
double get_comp_rand_real_d(comp_d x);

double get_comp_rand_mp(comp_mp x);
void   get_comp_rand_mp2(mpf_t mod, comp_mp x);
double get_comp_rand_real_mp(comp_mp x);
void   get_comp_rand_real_mp2(mpf_t mod, comp_mp x);

double get_comp_rand_rat(comp_d x_d, comp_mp x_mp, mpq_t *x_rat, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat);
double get_comp_rand_real_rat(comp_d x_d, comp_mp x_mp, mpq_t *x_rat, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat);

void   convert_point_data_d_to_mp(point_data_mp *dataMP, point_data_d *data); /* Changes path data from d to mp. */
void   convert_point_data_mp_to_d(point_data_d *data, point_data_mp *dataMP); /* Vice versa. */
void   point_data_cp_d(point_data_d *dest, point_data_d *src);  /* Copies all path data of src into dest. */
void   point_data_cp_mp(point_data_mp *dest, point_data_mp *src); /* same */

char* mpf_to_str(mpf_t MPF, int base); // MPF to string
size_t outStr_to_frac_size(size_t *numer_zeros, size_t *denom_zeros, char *strMan, long int exp); // find sizes to create string
void   outStr_to_frac(char *strFrac, char *strMan, size_t numer_zeros, long int exp, size_t denom_zeros); // setup the fraction string
void   mpf_t_to_rat(mpq_t RAT, mpf_t MP); // convert a multiprecision floating-point value to a rational number
void change_prec_mp(comp_mp x, int new_prec);
void change_prec_mp2(comp_mp x, int new_prec);

void   mypause(); /* Manual debugging function (esp. for memory testing) - causes Bertini to wait for (meaningless) integer input from user. */
void *bmalloc(size_t size); // Bertini's malloc with error checking
void *bcalloc(size_t num, size_t size); // Bertini's calloc with error checking
void *brealloc(void *ptr, size_t size); // Bertini's realloc with error checking
void bexit(int errorCode); // Bertini's exit function

double amp_criterion_A(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi);
double amp_criterion_B(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi, double tol, double residual, int maxIts, int it_number);
double amp_criterion_B2(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi, double tol, double proportion, int maxIts);
double amp_criterion_C(int safety_digits, double norm_J_inv, double Psi, double tol, double size);

void d_get_str(char *strOut, long int *exp, int digits, double D);  // converts D to a string for mantissa and exponent like mpf_get_str

void findFunctionResidual_conditionNumber_d(double *func_residual, double *cond_num, point_data_d *PD, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
void findFunctionResidual_d(double *func_residual, point_data_d *PD, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
void findFunctionResidual_conditionNumber_mp(mpf_t func_residual, double *cond_num, point_data_mp *PD, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
void findFunctionResidual_mp(mpf_t func_residual, point_data_mp *PD, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void init_all_mat_mp(int num, ...);
void init_all_mat_mp2(int prec, int num, ...);
void init_all_vec_mp(int num, ...);
void init_all_vec_mp2(int prec, int num, ...);
void init_all_point_mp(int num, ...);
void clear_all_mat_mp(int num, ...);
void clear_all_vec_mp(int num, ...);
void clear_all_point_mp(int num, ...);

int checkWritePrivilege();
int prec_to_digits(int prec);
int digits_to_prec(int digits);

void extrinsicToIntrinsic_d(point_d out, point_d in, mat_d B_transpose, vec_d p);
void intrinsicToExtrinsic_d(point_d out, point_d in, mat_d B, vec_d p);
void extrinsicToIntrinsic_mp(point_mp out, point_mp in, mat_mp B_transpose, vec_mp p);
void intrinsicToExtrinsic_mp(point_mp out, point_mp in, mat_mp B, vec_mp p);

int isSamePoint(point_d endPt1_d, point_mp endPt1_mp, int prec1, point_d endPt2_d, point_mp endPt2_mp, int prec2, double tol);
void findDiff_point(mpf_t norm_diff, point_d endPt1_d, point_mp endPt1_mp, int prec1, point_d endPt2_d, point_mp endPt2_mp, int prec2);

void intrinsicToExtrinsicSlices_d(mat_d B_out, vec_d p_out, mat_d B_in, vec_d p_in);
void intrinsicToExtrinsicSlices_mp(mat_mp B_out, vec_mp p_out, mat_mp B_in, vec_mp p_in);
void intrinsicToExtrinsicSlices_rat(mpq_t ****B_out, int *B_out_rows, int *B_out_cols, mpq_t ***p_out, int *p_out_size, mpq_t ***B_in, mpq_t **p_in, int B_in_rows, int B_in_cols, int p_in_size, int curr_prec, int max_prec);

void intrinsicToExtrinsicMat_d(mat_d B_out, mat_d B_in);
void intrinsicToExtrinsicMat_mp(mat_mp B_out, mat_mp B_in);
void intrinsicToExtrinsicMat_rat(mpq_t ****B_out, int *B_out_rows, int *B_out_cols, mpq_t ***B_in, int B_in_rows, int B_in_cols, int curr_prec, int max_prec);

void printRawMat(FILE *OUT, mat_d M_d, mat_mp M_mp, mpq_t ***M_rat, int MPType);
void printRawVec(FILE *OUT, vec_d V_d, vec_mp V_mp, mpq_t **V_rat, int MPType);
void printRawComp(FILE *OUT, comp_d C_d, comp_mp C_mp, mpq_t *C_rat, int MPType);

void setupRawComp(FILE *IN, comp_d C_d, comp_mp C_mp, mpq_t *C_rat, int C_Type, int inputType, int mp_need_init, int rat_need_init);
void setupRawVec(FILE *IN, vec_d V_d, vec_mp V_mp, mpq_t ***V_rat, int V_Type, int inputType);
void setupRawMat(FILE *IN, mat_d M_d, mat_mp M_mp, mpq_t ****M_rat, int M_Type, int inputType);

int sort_order_d(const void *vp, const void *vq);
int sort_order_mp(const void *vp, const void *vq);

double factorial(int n);
double factorial_array(int *array, int length);
void factorial2(mpf_t factorial, int n);
void factorial_array2(mpf_t rV, int *array, int length);
double combination(int d, int n);
int residual_to_digits_d(double residual);
int residual_to_digits_mp(mpf_t residual, int curr_prec);
int scanRestOfLine(FILE *IN);
int printRestOfLine(FILE *OUT, FILE *IN);
void setupInput(char *outputName, char *inputName);
void printVersion(FILE *OUT);
void normalize_vec_d(vec_d x1, vec_d x);
void normalize_vec_mp(vec_mp x1, vec_mp x);
double point_cp_d_norm(point_d outPt, point_d inPt);
void point_cp_mp_norm(mpf_t norm, point_mp outPt, point_mp inPt, int currPrec);
void point_cp_mp_norm2(mpf_t norm, point_mp inOutPt, int newPrec);
void point_d_to_mp_norm(mpf_t norm, point_mp outPt, point_d inPt, int currPrec);

/* Prototypes from diff.c: */
void diff(); /* Reads the function evaluation instructions from finalFile.out (from parser) and creates instructions for derivatives. */
               /* Writes to the file arr.out. */
int determine_dependencies(int numInst);  /* Goes through the numInst instructions for the eval of a func and determines which vars, subfuncs appear. */

/* Prototypes from diff_deflatable.c */
void diff_deflatable(int putIntoOrder, int allowNegExp); // create an SLP that is 'friendly' with deflation

/* Prototypes from endgame.c: */
int  endgame_d(int endgameNumber, int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *)); /* Sends data to appropriate endgame. */
int endgame_rank_d(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int  endgame_mp(int endgameNumber, int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int endgame_rank_mp(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int  endgame_amp(int endgameNumber, int pathNum, int *prec, double *first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int endgame_rank_amp(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

// Prototypes from CauchyEG.c - the Cauchy Integral endgame
double find_closed_loop_tol(double min_tol, double max_tol, point_data_d *x_d, point_data_mp *x_mp, int x_prec, tracker_config_t *T, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int CauchyEG_d(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d **endSamples, int *cycle, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_main_d(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int CauchyEG_rank_d(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void find_Cauchy_preEG_approx_d(point_d preEG_approx, comp_d finalTime, point_data_d *samples_d, int num_samples, int cycle_num, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int find_Cauchy_samples_d(double closed_loop_max_error, point_data_d **Samples, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_d *Start, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
void find_Cauchy_approx_d(point_d approx, point_data_d Samples[], int cycle_num, int num_samples);
int Cauchy_rank_main_d(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, point_data_d *Final, point_d last_approx, comp_d finalTime, point_data_d *samplePts, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int circle_track_d(point_data_d Final[], int Final_digits[], point_data_d *Start, int M, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int check_closed_loop_d(double closed_loop_max_error, point_data_d *Pt0, int *digitsCorrect0, point_data_d *Pt1, int *digitsCorrect1, comp_d time, eval_struct_d *e_d, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int CauchyEG_mp(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp **endSamples, int *cycle, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_main_mp(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int CauchyEG_rank_mp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void find_Cauchy_preEG_approx_mp(point_mp preEG_approx, comp_mp finalTime, point_data_mp *samples_mp, int num_samples, int cycle_num, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int find_Cauchy_samples_mp(double closed_loop_max_error, point_data_mp **Samples, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_mp *Start, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
void find_Cauchy_approx_mp(point_mp approx, point_data_mp Samples[], int cycle_num, int num_samples);
int Cauchy_rank_main_mp(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, point_data_mp *Final, point_mp last_approx, comp_mp finalTime, point_data_mp *samplePt, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int circle_track_mp(point_data_mp Final[], int Final_digits[], point_data_mp *Start, int M, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int check_closed_loop_mp(double closed_loop_max_error, point_data_mp *Pt0, int *digitsCorrect0, point_data_mp *Pt1, int *digitsCorrect1, comp_mp time, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int CauchyEG_amp(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endSamples_d, point_data_mp **endSamples_mp, int *cycle, int *samples_per_loop, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_main_amp(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, double *time_first_increase, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int CauchyEG_rank_amp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void find_Cauchy_preEG_approx_amp(point_d preEG_approx_d, point_mp preEG_approx_mp, int *preEG_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *samples_d, point_data_mp *samples_mp, int *samples_prec, int num_samples, int cycle_num, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int find_Cauchy_samples_amp(double closed_loop_max_error, point_data_d **Samples_d, point_data_mp **Samples_mp, int *Samples_prec, int *cycle_num, int num_samples, int *curr_cycle_size, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
void find_Cauchy_approx_amp(point_d approx_d, point_mp approx_mp, int *approx_prec, point_data_d Samples_d[], point_data_mp Samples_mp[], int prec_in, int cycle_num, int num_samples);
int Cauchy_rank_main_amp(double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int Cauchy_retVal, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *samplePt_d, point_data_mp *samplePt_mp, int prec_in, int cycle_num, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int circle_track_amp(int M, point_data_d Final_d[], point_data_mp Final_mp[], int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int check_closed_loop_amp(double closed_loop_max_error, point_data_d *Pt0_d, point_data_mp *Pt0_mp, int *prec0, int *digitsCorrect0, point_data_d *Pt1_d, point_data_mp *Pt1_mp, int *prec1, int *digitsCorrect1, comp_d time, eval_struct_d *e_d, eval_struct_mp *e_mp, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int Cauchy_corank_d(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_d approx0, point_d approx1, comp_d finalTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int Cauchy_corank_mp(double *CN, double *smallest_nonzero_SV, double *largest_zero_SV, point_mp approx0, point_mp approx1, comp_mp finalTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e_mp, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int Cauchy_corank_amp(double *CN, double *smallest_nonzero_SV, double *larget_zero_SV, point_d approx0_d, point_mp approx0_mp, int prec0, point_d approx1_d, point_mp approx1_mp, int prec1, comp_d finalTime_d, comp_mp finalTime_mp, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

// CauchyEG_new.c
int CauchyEG_main_d2(point_data_d *Final, point_d last_approx, point_data_d **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_d finalTime, point_data_d *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_main_mp2(point_data_mp *Final, point_mp last_approx, point_data_mp **endLoopSamples, int *cycle_num, int *samples_per_loop, comp_mp finalTime, point_data_mp *Start, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int CauchyEG_main_amp2(int *prec_out, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d **endLoopSamples_d, point_data_mp **endLoopSamples_mp, int *cycle_num, int *samples_per_loop, double *time_first_increase, comp_d finalTime_d, comp_mp finalTime_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

// trackBackEG.c
void find_minimum_separation(int *numPoints, double *minDist_d, double **startPoint_norm_d, point_d **startPts_d, mpf_t minDist_mp, mpf_t **startPoint_norm_mp, point_mp **startPts_mp, tracker_config_t *T, FILE *START);
int check_point_trackBack(int *indexI, int *indexJ, double tol, trackBack_samples_t **EGsamples, int numSamples, double norm_d, mpf_t norm_mp, point_d testPt_d, point_mp testPt_mp, int testPt_prec);
void zero_dim_trackBack_path_d(int pathNum, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void zero_dim_trackBack_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void zero_dim_trackBack_path_rank_d(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void zero_dim_trackBack_path_mp(int pathNum, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void zero_dim_trackBack_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
void zero_dim_trackBack_path_rank_mp(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

// Prototypes from NoEG.c - no endgame - straight tracking
int NoEG_d(point_data_d *Final, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int NoEG_mp(point_data_mp *Final, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int NoEG_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int NoEG_reactive_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int NoEG_amp2(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

// Prototypes from OneStopEG.c - tracking to t = endgameBoundary and then switch tolerances and track to t = targetT
int OneStopEG_d(point_data_d *Final, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int OneStopEG_mp(point_data_mp *Final, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int OneStopEG_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int OneStopEG_reactive_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// Prototypes from PSEG.c - the power series endgame
double compare_approximations_d(point_d approx0, point_d approx1);
double compare_approximations_mp(point_mp approx0, point_mp approx1);

int find_cycle_num_d(point_data_d *Pts, int start, double sample_factor);
int find_cycle_num_mp(point_data_mp *Pts, int start, double sample_factor);

void init_PSEG_samples_struct_d(PSEG_samples_struct_d *S, int num_points);
void clear_PSEG_samples_struct_d(PSEG_samples_struct_d *S);

void init_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S, int num_points, int prec);
void clear_PSEG_samples_struct_mp(PSEG_samples_struct_mp *S);

int PSEG_struct_d(point_data_d *Final, point_d last_approx, PSEG_samples_struct_d *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_struct_mp(point_data_mp *Final, point_mp last_approx, PSEG_samples_struct_mp *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int PSEG_samples_d(point_data_d *Final, int *goodSamples, int numSamples, point_data_d *samplesT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int PSEG_samples_mp(point_data_mp *Final, int *goodSamples, int numSamples, point_data_mp *samplesT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int PS_prediction_d(point_d predictedVal, comp_d predictedT, int *cycle_num, point_data_d *samples, int n, int maxCycleNum, int screenOut, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
void cubicInterp_d(point_d Res, _comp_d *T, _point_d *Y, _point_d *dHdT, comp_d T0);
int PSEG_d(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PS_prediction_mp(point_mp predictedVal, comp_mp predictedT, int *cycle_num, point_data_mp *samples, int n, int maxCycleNum, int screenOut, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
void cubicInterp_mp(point_mp Res, _comp_mp *T, _point_mp *Y, _point_mp *dHdT, comp_mp T0);
int PSEG_mp(int pathNum, point_data_mp *Final, point_mp last_approx_mp, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_rank_mp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_mp2(int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *lastSample, int *cycle_num, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int PSEG_rank_d(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_d2(int pathNum, point_data_d *Final, point_d last_approx, point_data_d *lastSample, int *cycle_num, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

int PSEG_reactive_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// Prototypes from PSEG_amp.c
void init_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S, int num_points, int prec);
void clear_PSEG_samples_struct_amp(PSEG_samples_struct_amp *S);

int PSEG_struct_amp(int *final_prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, PSEG_samples_struct_amp *PSEG_samples, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_amp(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_amp2(int pathNum, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *lastSample_d, point_data_mp *lastSample_mp, int *cycle_num, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));
int PSEG_rank_amp(int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_d last_approx_d, point_mp last_approx_mp, int *last_approx_prec, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

// Prototypes from ZeroOrderEG.c - the zero order endgame
int ZeroOrderEG_d(point_data_d *Final, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int ZeroOrderEG_samples_d(point_data_d *Final, int *goodSamples, int numSamples, point_data_d *samplesT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int ZeroOrderEG_mp(point_data_mp *Final, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int ZeroOrderEG_samples_mp(point_data_mp *Final, int *goodSamples, int numSamples, point_data_mp *samplesT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int ZeroOrderEG_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int ZeroOrderEG_reactive_amp(int *prec, double *time_first_increase, point_data_d *Final_d, point_data_mp *Final_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

/* Prototypes from adaptiveMP.c: */
int    AMPtrack_d(point_data_d *Final, point_data_d *Start, comp_d FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));  /* Basic tracker for AMP.  Makes use of the */
int    AMPtrack_mp(point_data_mp *Final, point_data_mp *Start, comp_mp FinalT, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *)); /* same */

/* Prototypes from svd.c: */
void   convert_bidiag_to_real_d(mat_d A, mat_d L, mat_d R);  /* Uses elem. matrices to convert bidiagonal complex A into real. */
void   convert_bidiag_to_real_mp(mat_mp A, mat_mp L, mat_mp R); /* same */
void   svd_d(mat_d A, mat_d L, mat_d R);  /* Stores singular vectors in L and R and overwrites A with diagonal matrix of singular values. */
void   svd_mp(mat_mp A, mat_mp L, mat_mp R);  /* same */
double condition_number_d(mat_d A);  /* Returns condition number of A without damaging the contents of A. */
double condition_number_mp(mat_mp A);  /* same */
int    find_num_rank_d(mat_d A, double tolerance);  /* Returns numerical rank and does not change A. */
int    find_num_rank_mp(mat_mp A, double tolerance);  /* same */
int    find_num_rank_and_sing_vecs_d(mat_d A, mat_d L, mat_d R, double tolerance);  /* Stores singular vectors in L and R, does not change A, and returns numerical rank. */
int    find_num_rank_and_sing_vecs_mp(mat_mp A, mat_mp L, mat_mp R, double tolerance);  /* same */
int    bidiag_svd_d(mat_d A, mat_d L, mat_d R);  /* Finds the SVD of a bidiagonal matrix. Uses algorithm from G. W. Stewart. */
int    bidiag_svd_mp(mat_mp A, mat_mp L, mat_mp R);  /* same */
void   prune_zeros_d(mat_d A, int *zero_rows, int *zero_cols); /* Shrinks A by cutting out rows and columns of all 0's. */
void   prune_zeros_mp(mat_mp A, int *zero_rows, int *zero_cols);  /* same */
void   reinsert_zeros_d(mat_d A, mat_d L, mat_d R, int *zero_rows, int *zero_cols); /* Plugs zero rows and columns back into matrix after SVD. */
void   reinsert_zeros_mp(mat_mp A, mat_mp L, mat_mp R, int *zero_rows, int *zero_cols);  /* same */

void   svd_via_newton(mat_mp M);  /* For the SVD/Newton project - not working yet (as of 3/1/05)!!!???!!! */

/* Prototypes from track.c: */
int	step_d(point_data_d *P, comp_d newTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); /* Called by track_d, tries to advance P to newTime using T, euler, and newton. */
int 	euler_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); /* Simple Euler's method. */
int 	constant_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); 
int 	newton_d(vec_d newParam, vec_d newPoint, tracker_config_t *T, FILE *OUT, comp_d t1, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); /* Simple Newton's method. */
int newton_iteration_d(double *residual, int return_norm, double *norm_J, double *norm_J_inv, vec_d P, comp_d time, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int newton_residual_d(vec_d dX, double *residual, int return_norm, double *norm_J, double *norm_J_inv, vec_d P, comp_d time, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int     refine_d(point_data_d *P, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *)); /* Uses Newton's method to refine to full precision available. */
int 	step_mp(point_data_mp *P, comp_mp newTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));  /* Similar to above _d version. */
int 	euler_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *)); /* Similar to above _d version. */
int 	constant_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *)); 
int 	newton_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp t1, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *)); /* Similar to above _d version. */
int newton_iteration_mp(mpf_t residual, int return_norm, double *norm_J, double *norm_J_inv, vec_mp P, comp_mp time, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int newton_residual_mp(vec_mp dX, mpf_t residual, int return_norm, double *norm_J, double *norm_J_inv, vec_mp P, comp_mp time, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int     refine_mp(point_data_mp *P, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *)); /* Similar to above _d version. */

void step_size_increase_d(double predictorError, int *consecSuccess, tracker_config_t *T);
void step_size_increase_mp(mpf_t predictorError, mpf_t currentStepSize, int *consecSuccess, tracker_config_t *T);

int step_error_d(double *predictorError, int consecSuccess, point_data_d *P, comp_d newTime, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int step_error_mp(mpf_t predictorError, int consecSuccess, point_data_mp *P, comp_mp newTime, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// ode.c
int modified_euler_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int heun_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int runge_kutta_d(vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int modified_euler_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int heun_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int runge_kutta_mp(vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// predictor methods with error estimates
int heun_error_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkn34_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk34_methods_d(double *err, double a[5], double a_minus_b[5], double c[5], double d[5][4], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkf45_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkck45_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk45_methods_d(double *err, double a[6], double a_minus_b[6], double c[6], double d[6][5], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkdp56_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk56_methods_d(double *err, double a[8], double a_minus_b[8], double c[8], double d[8][7], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkv67_d(double *err, vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk67_methods_d(double *err, double a[10], double a_minus_b[10], double c[10], double d[10][9], vec_d newParam, vec_d newPoint, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int heun_error_mp(mpf_t error, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rkn34_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rk34_methods_mp(mpf_t err, mpf_t a[5], mpf_t a_minus_b[5], mpf_t c[5], mpf_t d[5][4], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rkf45_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rkck45_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rk45_methods_mp(mpf_t err, mpf_t a[6], mpf_t a_minus_b[6], mpf_t c[6], mpf_t d[6][5], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rkdp56_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rk56_methods_mp(mpf_t err, mpf_t a[8], mpf_t a_minus_b[8], mpf_t c[8], mpf_t d[8][7], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rkv67_mp(mpf_t err, vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int rk67_methods_mp(mpf_t err, mpf_t a[10], mpf_t a_minus_b[10], mpf_t c[10], mpf_t d[10][9], vec_mp newParam, vec_mp newPoint, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// svd_analyze.c
int jacobian_analyze(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, int MPType, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int determineRankDef(double *CN, double rank_tol_d, mpf_t rank_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp,  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int prec_needed(int MPType, mpf_t min_sv, mpf_t func_resid, point_data_mp *in_mp, int prec_in, eval_struct_mp *e, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

/* Prototypes from prog_eval.c: */
/* These functions evaluate the polynomials and derivatives, storing function values in funcVals, parameter values in parVals (since the 
parameters depend upon the value of the path variable), derivatives of the parameters in parDer, jacobian w.r.t vars in Jv, jacobian w.r.t 
parameters in Jp, on input vars (the point), pathVars (the t-value), and Prog (the straight-line program). */
int eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void initEvalProg(int MPType);
void freeEvalProg(int MPType);
int evalProg_d_void(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int evalProg_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog);
int evalProg_eff_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog, int startFuncNum, int endFuncNum, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow);
int evalProg_d_std(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, prog_t *Prog);
void evalInsts_d(_comp_d *mem_d, int *prog, int startOp, int endOp, int oid);
int evalProg_mp_void(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int evalProg_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog);
int evalProg_eff_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog, int startFuncNum, int endFuncNum, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow);
int evalProg_mp_std(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, prog_t *Prog);
void evalInsts_mp(_comp_mp *mem_mp, int *prog, int startOp, int endOp, int oid);
int change_prec_prog(void const *ED, int new_prec);

/* Prototypes from prog_manip.c: */
int     write_out_prog(prog_t *Prog);  /* Writes straight-line program in "input" file format - uses a lot of memory, many size problems! */

/* Prototypes from min_svd.c: */
int min_svd_d(double *min_sv, mat_d A, double error_tol);
void min_svd_mp_prec(mpf_t min_sv, mat_mp Jv, int prec);
int min_svd_mp(mpf_t min_sv, mat_mp A, mpf_t error_tol);
void approx_min_svd_d(double *min_sv, mat_d A);
void approx_min_svd_mp(mpf_t min_sv, mat_mp A);
void approx_min_svd_mp_prec(mpf_t min_sv, mat_mp A, int prec);

/* Prototypes from corank.c: */
int rankDef_d(double *CN, double minSV0, double minSV1, double mat_norm, double max_CN, double max_SV_ratio, double SV_tol);
int corank_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0, mat_d mat1, double max_CN, double max_SV_ratio, double SV_tol);
int rankDef_mp(double *CN, mpf_t minSV0, mpf_t minSV1, mpf_t mat_norm, int curr_prec, double max_CN, double max_SV_ratio, double SV_tol);
int corank_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, mat_mp mat1, int curr_prec, double max_CN, double max_SV_ratio, double SV_tol);
int rankDef_amp(double *CN, mpf_t minSV0, int prec0, mpf_t minSV1, int prec1, mpf_t mat_norm, double max_SV_ratio);
int corank_amp(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0_d, mat_mp mat0_mp, int prec0, mat_d mat1_d, mat_mp mat1_mp, int prec1, double max_SV_ratio);

/* Prototypes from deflation.c: */
int deflation(int *deflations_needed, prog_t *deflatedProg, point_data_d *output_PD_d, point_data_mp *output_PD_mp, int *output_PD_prec, prog_t *origProg, int orig_corank, double smallest_nonzero_SV, double largest_zero_SV, point_data_d *input_PD_d, point_data_mp *input_PD_mp, int input_PD_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, mpq_t ***K, int K_rows, int K_cols, tracker_config_t *T, FILE *OUT, int deflation_iteration_bound);
void deflator(prog_t *new_prog, prog_t *old_prog, mpq_t ***K, int K_rows, int K_cols, mpq_t ***B, int B_rows, int B_cols); // single-stage deflator
void add_vec_patch_SLP(prog_t *new_prog, prog_t *old_prog, mpq_t **patch, int patch_size, mpq_t *patch_rhs);
void randomize_SLP(prog_t *new_prog, prog_t *old_prog, mpq_t ***A, int A_rows, int A_cols);

/* Prototypes from function_main.c */
int function_eval_main(int printJacobian, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);

/* Prototypes from newton_main.c */
int newton_eval_main(int computeCN, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);

// Basic OpenMP functions that are needed throughout
#ifdef _OPENMP
  #define max_threads() omp_get_max_threads()
  #define thread_num()  omp_get_thread_num()
#else
  #define max_threads() 1
  #define thread_num()  0
#endif

/* Basic arithmetic functions. */
/*******************************************************************************\ 
* A Preceeding 'd_' means that it should return a double regardless		*
* of the underlying precision. A trailing '_d' means that the arguments		*
* are expected to have an underlying type of double. A trailing '_mp' means	*
* that the underlying data types are multi-precision.				*
\*******************************************************************************/

// absolute value
#define d_abs_d(_a)  (sqrt((_a)->r*(_a)->r + (_a)->i*(_a)->i))
#define mpf_abs_mp(_a, _x) { int _oid = thread_num(); \
                             mpf_mul(_tempMPF1[_oid], (_x)->r, (_x)->r); mpf_mul(_a, (_x)->i, (_x)->i); \
                             mpf_add(_a, _a, _tempMPF1[_oid]); mpf_sqrt(_a, _a); }
#define abs_mp(_a, _x) { mpf_mul((_a)->r, (_x)->r, (_x)->r); mpf_mul((_a)->i, (_x)->i, (_x)->i); \
                         mpf_add((_a)->r, (_a)->r, (_a)->i); mpf_sqrt((_a)->r, (_a)->r); mpf_set_ui((_a)->i, 0); }

// one norm
#define d_oneNorm_d(_a) (fabs((_a)->r) + fabs((_a)->i))
#define d_oneNorm_mp(_a) (fabs(mpf_get_d((_a)->r)) + fabs(mpf_get_d((_a)->i)))
#define mp_oneNorm_mp(_n,_a) { int _oid = thread_num();  mpf_abs(_n, (_a)->r); \
  mpf_abs(_tempMPF1[_oid], (_a)->i); mpf_add((_n), (_n), _tempMPF1[_oid]); }

// addition
#define add_d(_r,_a,_b)  { (_r)->r = (_a)->r + (_b)->r; (_r)->i = (_a)->i + (_b)->i; }
#define add_mp(_r,_a,_b) { mpf_add((_r)->r, (_a)->r, (_b)->r); mpf_add((_r)->i, (_a)->i, (_b)->i); }

// subtraction
#define sub_d(_r,_a,_b)  { (_r)->r = (_a)->r - (_b)->r; (_r)->i = (_a)->i - (_b)->i; }
#define sub_mp(_r,_a,_b) { mpf_sub((_r)->r, (_a)->r, (_b)->r); mpf_sub((_r)->i, (_a)->i, (_b)->i); }

// mutliply
#define mul_d(_r,_a,_b) { double _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; \
                         (_r)->i = (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r = _s;}
#define mul_d2(_r,_a,_b,_s) { _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i = (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r = _s;}

#define mul_mp(_r, _a, _b) { int _oid = thread_num(); \
           mpf_mul(_tempMPF1[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_b)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
           mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->i); mpf_mul((_r)->i, (_a)->i, (_b)->r); mpf_add((_r)->i, (_r)->i, _tempMPF2[_oid]); \
           mpf_set((_r)->r, _tempMPF1[_oid]);}
#define mul_omp_mp(_r, _a, _b, _oid) { \
           mpf_mul(_tempMPF1[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_b)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
           mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->i); mpf_mul((_r)->i, (_a)->i, (_b)->r); mpf_add((_r)->i, (_r)->i, _tempMPF2[_oid]); \
           mpf_set((_r)->r, _tempMPF1[_oid]); }
#define mul_rat(_r, _a, _b) { mpq_t _t1, _t2; mpq_init(_t1); mpq_init(_t2); \
  mpq_mul(_t1, (_a)[0], (_b)[0]); mpq_mul(_t2, (_a)[1], (_b)[1]); mpq_sub(_t1, _t1, _t2); \
  mpq_mul(_t2, (_a)[0], (_b)[1]); mpq_mul((_r)[1], (_a)[1], (_b)[0]); mpq_add((_r)[1], (_r)[1], _t2); \
  mpq_set((_r)[0], _t1); mpq_clear(_t1); mpq_clear(_t2); }

#define mul_rdouble_d(_r,_a,_d)  { (_r)->r = (_a)->r * (_d); (_r)->i = (_a)->i * (_d); }
#define mul_rdouble_mp(_r,_a,_d)  {  int _oid = thread_num(); mpf_set_d(_tempMPF1[_oid], _d); \
                                     mpf_mul((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_mul((_r)->i, (_a)->i, _tempMPF1[_oid]); }
#define mul_rmpf_mp(_r,_a,_d)  { mpf_mul((_r)->r, (_a)->r, _d); mpf_mul((_r)->i, (_a)->i, _d); }

#define mul_ui_mp(_r,_a,_d)  { mpf_mul_ui((_r)->r, (_a)->r, _d); mpf_mul_ui((_r)->i, (_a)->i, _d); }
#define div_ui_mp(_r,_a,_d)  { mpf_div_ui((_r)->r, (_a)->r, _d); mpf_div_ui((_r)->i, (_a)->i, _d); }


// divide
#define div_d(_r,_a,_b) { double _s = 1 / ((_b)->r*(_b)->r + (_b)->i*(_b)->i);\
                          double _rt = _s * ((_a)->r*(_b)->r + (_a)->i*(_b)->i);\
                             (_r)->i = _s * ((_a)->i*(_b)->r - (_a)->r*(_b)->i); (_r)->r = _rt; }
#define div_d2(_r,_a,_b, _s, _rt) { _s = 1 / ((_b)->r*(_b)->r + (_b)->i*(_b)->i); _rt = _s * ((_a)->r*(_b)->r + (_a)->i*(_b)->i);\
                             (_r)->i = _s * ((_a)->i*(_b)->r - (_a)->r*(_b)->i); (_r)->r = _rt; }

#define div_mp(_r, _a, _b) { int _oid = thread_num(); \
  mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
    mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
    mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
  mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
  mpf_set((_r)->r, _tempMPF2[_oid]); }
#define div_omp_mp(_r, _a, _b, _oid) { \
  mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
    mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
    mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
  mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
  mpf_set((_r)->r, _tempMPF2[_oid]); }
#define div_mp2(_r, _a, _b) { int _oid = thread_num(), _prec = mpf_get_prec((_a)->r), _prec2 = mpf_get_prec((_b)->r), _curr_prec[3]; \
 if (_prec > _prec2) _prec = _prec2;\
   _curr_prec[0] = mpf_get_prec(_tempMPF1[_oid]); _curr_prec[1] = mpf_get_prec(_tempMPF2[_oid]); _curr_prec[2] = mpf_get_prec(_tempMPF3[_oid]);\
 if (_prec > _curr_prec[0]) mpf_set_prec(_tempMPF1[_oid], _prec);\
 if (_prec > _curr_prec[1]) mpf_set_prec(_tempMPF2[_oid], _prec);\
 if (_prec > _curr_prec[2]) mpf_set_prec(_tempMPF3[_oid], _prec);\
 mpf_mul(_tempMPF1[_oid], (_b)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_b)->i); mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
   mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
 mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->i); \
   mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); \
 mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_mul((_r)->i, (_a)->r, (_b)->i); mpf_sub((_r)->i, _tempMPF3[_oid], (_r)->i); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); \
 mpf_set((_r)->r, _tempMPF2[_oid]); \
 if (_prec > _curr_prec[0]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[0]);\
 if (_prec > _curr_prec[1]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[1]);\
 if (_prec > _curr_prec[2]) mpf_set_prec(_tempMPF1[_oid], _curr_prec[2]); }

// reciprocate
#define recip_d(_r, _a) { double _s = 1 / ((_a)->r*(_a)->r + (_a)->i*(_a)->i); (_r)->r = (_a)->r * _s; (_r)->i = - (_a)->i * _s; }
#define recip_d2(_r, _a, _s) { _s = 1 / ((_a)->r*(_a)->r + (_a)->i*(_a)->i); (_r)->r = (_a)->r * _s; (_r)->i = - (_a)->i * _s; }
#define recip_mp(_r, _a) { int _oid = thread_num(); \
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i); \
  mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_ui_div(_tempMPF1[_oid], 1, _tempMPF1[_oid]); \
  mpf_mul((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_neg(_tempMPF1[_oid], _tempMPF1[_oid]); mpf_mul((_r)->i, (_a)->i, _tempMPF1[_oid]); }
#define recip_rat(_r, _a) { mpq_t _t1, _t2; mpq_init(_t1); mpq_init(_t2); \
  mpq_mul(_t1, (_a)[0], (_a)[0]); mpq_mul(_t1, (_a)[1], (_a)[1]); mpq_add(_t1, _t1, _t1); mpq_inv(_t1, _t1); \
  mpq_mul((_r)[0], (_r)[0], _t1); mpq_neg(_t1, _t1); mpq_mul((_r)[1], (_a)[1], _t1); mpq_clear(_t1); mpq_clear(_t2); }

// exponentiate
#define sqr_d(_r, _a) { double _s = (_a)->r*(_a)->r - (_a)->i*(_a)->i; (_r)->i = 2*(_a)->r*(_a)->i; (_r)->r = _s; }
#define sqr_d2(_r, _a, _s) { _s = (_a)->r*(_a)->r - (_a)->i*(_a)->i; (_r)->i = 2*(_a)->r*(_a)->i; (_r)->r = _s; }

#define sqr_mp(_r, _a) { int _oid = thread_num();\
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_a)->i); mpf_add((_r)->i, _tempMPF2[_oid], _tempMPF2[_oid]); mpf_set((_r)->r, _tempMPF1[_oid]);}
#define sqr_omp_mp(_r, _a, _oid) { \
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); \
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_a)->i); mpf_add((_r)->i, _tempMPF2[_oid], _tempMPF2[_oid]); mpf_set((_r)->r, _tempMPF1[_oid]);}

#define cube_d(_r, _a) { double _s = (_a)->r * ((_a)->r*(_a)->r - 3*(_a)->i*(_a)->i); (_r)->i = (_a)->i * (3*(_a)->r*(_a)->r - (_a)->i*(_a)->i); \
                         (_r)->r = _s; }
#define cube_d2(_r, _a, _s) { _s = (_a)->r * ((_a)->r*(_a)->r - 3*(_a)->i*(_a)->i); (_r)->i = (_a)->i * (3*(_a)->r*(_a)->r - (_a)->i*(_a)->i); \
                         (_r)->r = _s; }

#define cube_mp(_r, _a) { int _oid = thread_num(); \
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i);\
  mpf_mul_ui(_tempMPF3[_oid], _tempMPF1[_oid], 3); mpf_sub(_tempMPF3[_oid], _tempMPF3[_oid], _tempMPF2[_oid]);\
  mpf_mul_ui(_tempMPF2[_oid], _tempMPF2[_oid], 3); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]);\
  mpf_mul((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_mul((_r)->i, (_a)->i, _tempMPF3[_oid]); }
#define cube_omp_mp(_r, _a, _oid) { \
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_a)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_a)->i);\
  mpf_mul_ui(_tempMPF3[_oid], _tempMPF1[_oid], 3); mpf_sub(_tempMPF3[_oid], _tempMPF3[_oid], _tempMPF2[_oid]);\
  mpf_mul_ui(_tempMPF2[_oid], _tempMPF2[_oid], 3); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]);\
  mpf_mul((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_mul((_r)->i, (_a)->i, _tempMPF3[_oid]); }

// used for rational functions (negative integer exponents)
#define exp_d_neg(_r, _a, _d, _s) { int _dd = -_d; \
           if (_dd == 2)  { sqr_d2(_r, _a, _s); }\
      else if (_dd == 3)  { cube_d2(_r, _a, _s); } \
      else if (_dd == 4)  { sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s); }\
      else if (_dd == 5)  { comp_d _t; set_d(_t, _a); sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s);}\
      else if (_dd == 6)  { cube_d2(_r, _a, _s); sqr_d2(_r, _r, _s); }\
      else if (_dd == 7)  { comp_d _t; set_d(_t, _a); cube_d2(_r, _a, _s); sqr_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s); }\
      else if (_dd == 8)  { sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s);  sqr_d2(_r, _r, _s); }\
      else if (_dd == 9)  { cube_d2(_r, _a, _s); cube_d2(_r, _r, _s); }\
      else if (_dd == 10) { comp_d _t; set_d(_t, _a); cube_d2(_r, _a, _s); cube_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s); }\
      else if (_dd == 1)  { set_d(_r, _a); }\
      else if (_dd == 0)  { set_one_d(_r); }\
      else if (_dd > 10)  { int _i,_size; comp_d _t; set_d(_t,_a); set_one_d(_r); \
       mpz_t _pow_int; mpz_init_set_ui(_pow_int, _dd); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
       for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { mul_d2(_r,_r,_t,_s); } sqr_d2(_t,_t,_s); } \
       mpz_clear(_pow_int); free(_base2); } \
      else { printf("\nERROR: Invalid power (%d) in exp_d_neg.\n",(int) _d); bexit(0); } \
      recip_d2(_r,_r,_s); }

#define exp_mp_neg(_r, _a, _d, _oid) { int _dd = -_d; \
           if (_dd == 2)  { sqr_omp_mp(_r, _a, _oid); }\
      else if (_dd == 3)  { cube_omp_mp(_r, _a, _oid); } \
      else if (_dd == 4)  { sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_dd == 5)  { mpf_t _t; mpf_init2(_t, mpf_get_prec(_tempMPF1[_oid])); mpf_set(_t, (_a)->r); mpf_set(_tempMPF3[_oid], (_a)->i); \
                           sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); \
                           mpf_mul(_tempMPF2[_oid], (_r)->r, _tempMPF3[_oid]); mpf_mul((_r)->r, (_r)->r, _t); mpf_mul(_tempMPF1[_oid], (_r)->i, _tempMPF3[_oid]); \
                           mpf_sub((_r)->r, (_r)->r, _tempMPF1[_oid]); mpf_mul((_r)->i, (_r)->i, _t); mpf_add((_r)->i, (_r)->i, _tempMPF2[_oid]); mpf_clear(_t); } \
      else if (_dd == 6)  { cube_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_dd == 7)  { comp_mp _t; int _prec = mpf_get_prec(_tempMPF1[_oid]); init_mp2(_t, _prec); set_mp(_t, _a); \
                           cube_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); mul_omp_mp(_r, _r, _t, _oid); clear_mp(_t); }\
      else if (_dd == 8)  { sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_dd == 9)  { cube_omp_mp(_r, _a, _oid); cube_omp_mp(_r, _r, _oid); }\
      else if (_dd == 10) { comp_mp _t; int _prec = mpf_get_prec(_tempMPF1[_oid]); init_mp2(_t, _prec); set_mp(_t, _a); \
                           cube_omp_mp(_r, _a, _oid); cube_omp_mp(_r, _r, _oid); mul_omp_mp(_r, _r, _t, _oid); clear_mp(_t); }\
      else if (_dd == 1)  { set_mp(_r, _a); }\
      else if (_dd == 0)  { set_one_mp(_r); }\
      else if (_dd > 10)  { int _i,_size,_j = mpf_get_prec(_tempMPF1[_oid]); comp_mp _t; init_mp2(_t,_j); set_mp(_t,_a); set_one_mp(_r); \
        mpz_t _pow_int; mpz_init_set_ui(_pow_int,_dd); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
        for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { mul_omp_mp(_r,_r,_t,_oid); } sqr_omp_mp(_t,_t,_oid); } \
        mpz_clear(_pow_int); free(_base2); clear_mp(_t); } \
      else { printf("\nERROR: Invalid power (%d) in exp_mp_neg.\n",(int) _d); bexit(0); } \
      recip_mp(_r, _r); }

// general exponentiation functions
#define exp_d(_r, _a, _d) { double _temp; exp_d2(_r,_a,_d,_temp); }

#define exp_d2(_r, _a, _d, _s) { /* determine if exponent is integer or floating point */  \
      if (fabs(_d - ((int) _d)) > 0) { /* assume base is real */ (_r)->r = pow((_a)->r,_d); (_r)->i = 0; } \
      else if (_d == 2)  { sqr_d2(_r, _a, _s); }\
      else if (_d == 3)  { cube_d2(_r, _a, _s); } \
      else if (_d == 4)  { sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s); }\
      else if (_d == 5)  { comp_d _t; set_d(_t, _a); sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s);}\
      else if (_d == 6)  { cube_d2(_r, _a, _s); sqr_d2(_r, _r, _s); }\
      else if (_d == 7)  { comp_d _t; set_d(_t, _a); cube_d2(_r, _a, _s); sqr_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s); }\
      else if (_d == 8)  { sqr_d2(_r, _a, _s); sqr_d2(_r, _r, _s);  sqr_d2(_r, _r, _s); }\
      else if (_d == 9)  { cube_d2(_r, _a, _s); cube_d2(_r, _r, _s); }\
      else if (_d == 10) { comp_d _t; set_d(_t, _a); cube_d2(_r, _a, _s); cube_d2(_r, _r, _s); mul_d2(_r, _r, _t, _s); }\
      else if (_d == 1)  { set_d(_r, _a); }\
      else if (_d == 0)  { set_one_d(_r); }\
      else if (_d > 10)  { int _i,_size; comp_d _t; set_d(_t,_a); set_one_d(_r); \
       mpz_t _pow_int; mpz_init_set_ui(_pow_int, _d); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
       for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { mul_d2(_r,_r,_t,_s); } sqr_d2(_t,_t,_s); } \
       mpz_clear(_pow_int); free(_base2); } \
      else if (_d < 0) { /* use negative exponentiation */ exp_d_neg(_r,_a,_d,_s); } \
      else { printf("\nERROR: Invalid power (%d) in exp_d2.\n", (int) _d); bexit(0); } }

#define exp_mp(_r, _a, _d) { int _oid = thread_num(); exp_omp_mp(_r, _a, _d, _oid); }
#define exp_mp_int(_r, _a, _d) { int _oid = thread_num(); exp_omp_mp_int(_r, _a, _d, _oid); }

#define exp_omp_mp_int(_r, _a, _d, _oid) { /* _d is an integer */ \
           if (_d == 2)  { sqr_omp_mp(_r, _a, _oid); }\
      else if (_d == 3)  { cube_omp_mp(_r, _a, _oid); } \
      else if (_d == 4)  { sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_d == 5)  { mpf_t _t; mpf_init2(_t, mpf_get_prec(_tempMPF1[_oid])); mpf_set(_t, (_a)->r); mpf_set(_tempMPF3[_oid], (_a)->i); \
                           sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); \
                           mpf_mul(_tempMPF2[_oid], (_r)->r, _tempMPF3[_oid]); mpf_mul((_r)->r, (_r)->r, _t); mpf_mul(_tempMPF1[_oid], (_r)->i, _tempMPF3[_oid]); \
                           mpf_sub((_r)->r, (_r)->r, _tempMPF1[_oid]); mpf_mul((_r)->i, (_r)->i, _t); mpf_add((_r)->i, (_r)->i, _tempMPF2[_oid]); mpf_clear(_t); } \
      else if (_d == 6)  { cube_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_d == 7)  { comp_mp _t; int _prec = mpf_get_prec(_tempMPF1[_oid]); init_mp2(_t, _prec); set_mp(_t, _a); \
                           cube_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); mul_omp_mp(_r, _r, _t, _oid); clear_mp(_t); }\
      else if (_d == 8)  { sqr_omp_mp(_r, _a, _oid); sqr_omp_mp(_r, _r, _oid); sqr_omp_mp(_r, _r, _oid); }\
      else if (_d == 9)  { cube_omp_mp(_r, _a, _oid); cube_omp_mp(_r, _r, _oid); }\
      else if (_d == 10) { comp_mp _t; int _prec = mpf_get_prec(_tempMPF1[_oid]); init_mp2(_t, _prec); set_mp(_t, _a); \
                           cube_omp_mp(_r, _a, _oid); cube_omp_mp(_r, _r, _oid); mul_omp_mp(_r, _r, _t, _oid); clear_mp(_t); }\
      else if (_d == 1)  { set_mp(_r, _a); }\
      else if (_d == 0)  { mpf_set_ui((_r)->r, 1); mpf_set_ui((_r)->i, 0); }\
      else if (_d > 10)  { int _i,_size,_j = mpf_get_prec(_tempMPF1[_oid]); comp_mp _t; init_mp2(_t,_j); set_mp(_t,_a); set_one_mp(_r); \
        mpz_t _pow_int; mpz_init_set_ui(_pow_int,_d); char *_base2 = mpz_get_str(NULL, 2, _pow_int); _size = strlen(_base2); \
        for (_i = _size - 1; _i >= 0; _i--) { if (_base2[_i] == '1') { mul_omp_mp(_r,_r,_t,_oid); } sqr_omp_mp(_t,_t,_oid); } \
        mpz_clear(_pow_int); free(_base2); clear_mp(_t); } \
      else if (_d < 0) { /* use negative exponentiation */ exp_mp_neg(_r,_a,_d,_oid); } \
      else { printf("\nERROR: Invalid power (%d) in exp_mp.\n", _d); bexit(0); } }

#define exp_omp_mp(_r, _a, _d, _oid) { /* determine if exponent is integer or floating point */  \
      if (mpfr_integer_p(_d)) { /* _d is an integer */ int _D = (int) mpf_get_d(_d); exp_omp_mp_int(_r, _a, _D, _oid); } \
      else { /* _d is a floating point - assume base is real */ mpfr_pow((_r)->r, (_a)->r, _d, __gmp_default_rounding_mode); mpf_set_ui((_r)->i, 0); } }

// e^x
#define exponential_d(_r, _a) { double _s = exp((_a)->r); (_r)->r = _s * cos((_a)->i); (_r)->i = _s * sin((_a)->i); }
#define exponential_d2(_r,_a,_s) { _s = exp((_a)->r); (_r)->r = _s * cos((_a)->i); (_r)->i = _s * sin((_a)->i); }

#define exponential_mp(_r, _a) { int _oid = thread_num(); mpfr_exp(_tempMPF1[_oid], (_a)->r, __gmp_default_rounding_mode); \
	mpfr_cos((_r)->r, (_a)->i, __gmp_default_rounding_mode); mpf_mul((_r)->r, (_r)->r, _tempMPF1[_oid]); \
	mpfr_sin((_r)->i, (_a)->i, __gmp_default_rounding_mode); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); }
#define exponential_omp_mp(_r, _a, _oid) { mpfr_exp(_tempMPF1[_oid], (_a)->r, __gmp_default_rounding_mode); \
        mpfr_cos((_r)->r, (_a)->i, __gmp_default_rounding_mode); mpf_mul((_r)->r, (_r)->r, _tempMPF1[_oid]); \
        mpfr_sin((_r)->i, (_a)->i, __gmp_default_rounding_mode); mpf_mul((_r)->i, (_r)->i, _tempMPF1[_oid]); }

// cos(x)
#define cos_d(_r, _a) { double _s = exp((_a)->i); double _sinv = 1/_s; (_r)->r = 0.5 * cos((_a)->r) * (_s + _sinv); (_r)->i = 0.5 * sin((_a)->r) * (_sinv - _s); }
#define cos_d2(_r, _a, _s, _sinv) {   _s = exp((_a)->i); _sinv = 1/_s; (_r)->r = 0.5 * cos((_a)->r) * (_s + _sinv); (_r)->i = 0.5 * sin((_a)->r) * (_sinv - _s); }

#define cos_mp(_r, _a) { int _oid = thread_num(); mpfr_exp(_tempMPF1[_oid], (_a)->i, __gmp_default_rounding_mode); mpf_ui_div(_tempMPF2[_oid], 1, _tempMPF1[_oid]); \
	mpfr_sin((_r)->i, (_a)->r, __gmp_default_rounding_mode); mpfr_cos((_r)->r, (_a)->r, __gmp_default_rounding_mode); \
	mpf_div_ui((_r)->r, (_r)->r, 2); mpf_div_ui((_r)->i, (_r)->i, 2); \
	mpf_add(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->r, (_r)->r, _tempMPF3[_oid]); \
	mpf_sub(_tempMPF3[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); mpf_mul((_r)->i, (_r)->i, _tempMPF3[_oid]); }
#define cos_omp_mp(_r, _a, _oid) { mpfr_exp(_tempMPF1[_oid], (_a)->i, __gmp_default_rounding_mode); mpf_ui_div(_tempMPF2[_oid], 1, _tempMPF1[_oid]); \
        mpfr_sin((_r)->i, (_a)->r, __gmp_default_rounding_mode); mpfr_cos((_r)->r, (_a)->r, __gmp_default_rounding_mode); \
        mpf_div_ui((_r)->r, (_r)->r, 2); mpf_div_ui((_r)->i, (_r)->i, 2); \
        mpf_add(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->r, (_r)->r, _tempMPF3[_oid]); \
        mpf_sub(_tempMPF3[_oid], _tempMPF2[_oid], _tempMPF1[_oid]); mpf_mul((_r)->i, (_r)->i, _tempMPF3[_oid]); }

// sin(x)
#define sin_d(_r, _a) { double _s = exp((_a)->i); double _sinv = 1/_s; (_r)->r = 0.5 * sin((_a)->r) * (_s + _sinv); (_r)->i = 0.5 * cos((_a)->r) * (_s - _sinv); }
#define sin_d2(_r, _a, _s, _sinv) {   _s = exp((_a)->i); _sinv = 1/_s; (_r)->r = 0.5 * sin((_a)->r) * (_s + _sinv); (_r)->i = 0.5 * cos((_a)->r) * (_s - _sinv); }

#define sin_mp(_r, _a) { int _oid = thread_num(); mpfr_exp(_tempMPF1[_oid], (_a)->i, __gmp_default_rounding_mode); mpf_ui_div(_tempMPF2[_oid], 1, _tempMPF1[_oid]); \
        mpfr_cos((_r)->i, (_a)->r, __gmp_default_rounding_mode); mpfr_sin((_r)->r, (_a)->r, __gmp_default_rounding_mode); \
        mpf_div_ui((_r)->r, (_r)->r, 2); mpf_div_ui((_r)->i, (_r)->i, 2); \
        mpf_add(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->r, (_r)->r, _tempMPF3[_oid]); \
        mpf_sub(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->i, (_r)->i, _tempMPF3[_oid]); }
#define sin_omp_mp(_r, _a, _oid) { mpfr_exp(_tempMPF1[_oid], (_a)->i, __gmp_default_rounding_mode); mpf_ui_div(_tempMPF2[_oid], 1, _tempMPF1[_oid]); \
        mpfr_cos((_r)->i, (_a)->r, __gmp_default_rounding_mode); mpfr_sin((_r)->r, (_a)->r, __gmp_default_rounding_mode); \
        mpf_div_ui((_r)->r, (_r)->r, 2); mpf_div_ui((_r)->i, (_r)->i, 2); \
        mpf_add(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->r, (_r)->r, _tempMPF3[_oid]); \
        mpf_sub(_tempMPF3[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul((_r)->i, (_r)->i, _tempMPF3[_oid]); }

// _r += _a * _b
#define sum_mul_d(_r, _a, _b) { double _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i += (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r += _s; }
#define sum_mul_d2(_r, _a, _b, _s) { _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i += (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r += _s; }

#define sum_mul_mp(_r, _a, _b) { int _oid = thread_num();\
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_b)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]);\
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->i); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]);\
  mpf_add((_r)->r, (_r)->r, _tempMPF1[_oid]); mpf_add((_r)->i, (_r)->i, _tempMPF2[_oid]); }

// _r -= _a * _b
#define sub_mul_d(_r, _a, _b) { double _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i -= (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r -= _s; }
#define sub_mul_d2(_r, _a, _b, _s) { _s = (_a)->r * (_b)->r - (_a)->i * (_b)->i; (_r)->i -= (_a)->r * (_b)->i + (_a)->i * (_b)->r; (_r)->r -= _s; }

#define sub_mul_mp(_r, _a, _b) { int _oid = thread_num();\
  mpf_mul(_tempMPF1[_oid], (_a)->r, (_b)->r); mpf_mul(_tempMPF2[_oid], (_a)->i, (_b)->i); mpf_sub(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]);\
  mpf_mul(_tempMPF2[_oid], (_a)->r, (_b)->i); mpf_mul(_tempMPF3[_oid], (_a)->i, (_b)->r); mpf_add(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]);\
  mpf_sub((_r)->r, (_r)->r, _tempMPF1[_oid]); mpf_sub((_r)->i, (_r)->i, _tempMPF2[_oid]); }

// _r += conj(_a) * _b
#define sum_conj_mul_d(_r, _a, _b) { double _s = (_a)->r * (_b)->r + (_a)->i * (_b)->i; (_r)->i += (_a)->r * (_b)->i - (_a)->i * (_b)->r; (_r)->r += _s; }
#define sum_conj_mul_d2(_r, _a, _b, _s) { _s = (_a)->r * (_b)->r + (_a)->i * (_b)->i; (_r)->i += (_a)->r * (_b)->i - (_a)->i * (_b)->r; (_r)->r += _s; }

// _r = _a - 2 * _b * _c
#define householder_mult_left_d(_r, _a, _b, _c) { double _s = (_b)->r * (_c)->r - (_b)->i * (_c)->i; \
    (_r)->i = (_a)->i - 2 * ((_b)->i * (_c)->r + (_b)->r * (_c)->i); (_r)->r = (_a)->r - 2 * _s; }
#define householder_mult_left_d2(_r, _a, _b, _c, _s) { _s = (_b)->r * (_c)->r - (_b)->i * (_c)->i; \
    (_r)->i = (_a)->i - 2 * ((_b)->i * (_c)->r + (_b)->r * (_c)->i); (_r)->r = (_a)->r - 2 * _s; }

// _r = _a - 2 * _b * conj(_c)
#define householder_mult_right_d(_r, _a, _b, _c) { double _s = (_b)->r * (_c)->r + (_b)->i * (_c)->i; \
    (_r)->i = (_a)->i - 2 * ((_b)->i * (_c)->r - (_b)->r * (_c)->i); (_r)->r = (_a)->r - 2 * _s; }
#define householder_mult_right_d2(_r, _a, _b, _c, _s) { _s = (_b)->r * (_c)->r + (_b)->i * (_c)->i; \
    (_r)->i = (_a)->i - 2 * ((_b)->i * (_c)->r - (_b)->r * (_c)->i); (_r)->r = (_a)->r - 2 * _s; }

#define householder_mult_right_mp(_r, _a, _b, _c) { int _oid = thread_num(); \
    mpf_mul(_tempMPF1[_oid], (_b)->r, (_c)->r); mpf_mul(_tempMPF2[_oid], (_b)->i, (_c)->i); \
    mpf_add(_tempMPF1[_oid], _tempMPF1[_oid], _tempMPF2[_oid]); mpf_mul_ui(_tempMPF1[_oid], _tempMPF1[_oid], 2); \
    mpf_mul(_tempMPF2[_oid], (_b)->i, (_c)->r); mpf_mul(_tempMPF3[_oid], (_b)->r, (_c)->i); \
    mpf_sub(_tempMPF2[_oid], _tempMPF2[_oid], _tempMPF3[_oid]); mpf_mul_ui(_tempMPF2[_oid], _tempMPF2[_oid], 2); \
    mpf_sub((_r)->r, (_a)->r, _tempMPF1[_oid]); mpf_sub((_r)->i, (_a)->i, _tempMPF2[_oid]); }

// initialize
#define init_d(_a)  { (_a)->r = (_a)->i = 0; }
#define init_mp(_a) { mpfr_init((_a)->r); mpfr_init((_a)->i); }
#define init_mp2(_a, _prec) {mpfr_init2((_a)->r, _prec), mpfr_init2((_a)->i, _prec); }
#define init_rat(_a) { mpq_init((_a)[0]); mpq_init((_a)[1]); }

// set precision
#define setprec_mp(_a,_prec) {mpfr_set_prec((_a)->r,_prec); mpfr_set_prec((_a)->i,_prec); }
#define setprec_mp_rat(_a, _prec, _rat) { mpfr_set_prec((_a)->r, _prec); mpf_set_q((_a)->r, (_rat)[0]); mpfr_set_prec((_a)->i, _prec); mpf_set_q((_a)->i, (_rat)[1]); }
#define change_prec_mp_rat setprec_mp_rat

// copy value
#define set_d(_r,_a)   { (_r)->r = (_a)->r; (_r)->i = (_a)->i; }
#define set_mp(_r,_a)  { mpf_set((_r)->r, (_a)->r); mpf_set((_r)->i, (_a)->i); }
#define set_rat(_r,_a) { mpq_set((_r)[0], (_a)[0]); mpq_set((_r)[1], (_a)[1]); }

// negate value
#define neg_d(_r,_a)  { (_r)->r = -(_a)->r; (_r)->i = -(_a)->i; }
#define neg_mp(_r,_a) { mpf_neg((_r)->r, (_a)->r); mpf_neg((_r)->i, (_a)->i); }
#define neg_rat(_r,_a) { mpq_neg((_r)[0], (_a)[0]); mpq_neg((_r)[1], (_a)[1]); }

// conjugate value
#define conjugate_d(_r, _a) { (_r)->r = (_a)->r, (_r)->i = -(_a)->i; }
#define neg_conjugate_d(_r, _a) { (_r)->r = -(_a)->r, (_r)->i = (_a)->i; }
#define conjugate_mp(_r, _a) { mpf_set((_r)->r, (_a)->r); mpf_neg((_r)->i, (_a)->i); }
#define neg_conjugate_mp(_r, _a) { mpf_neg((_r)->r, (_a)->r); mpf_set((_r)->i, (_a)->i); }

// set real & imag values
#define set_double_d(_r, _x, _y) { (_r)->r = _x; (_r)->i = _y; }
#define set_double_mp(_r, _x, _y) { mpf_set_d((_r)->r, _x); mpf_set_d((_r)->i, _y); }
#define set_double_mp2(_r, _x, _y) { mpf_set_d((_r)->r, _x); mpf_set_d((_r)->i, _y); }

// set zero
#define set_zero_d(_r) { (_r)->r = (_r)->i = 0; }
#define set_zero_mp(_r) { mpf_set_ui((_r)->r, 0); mpf_set_ui((_r)->i, 0); }
#define set_zero_rat(_r) { mpq_set_ui((_r)[0], 0, 1); mpq_set_ui((_r)[1], 0, 1); }

// set one
#define set_one_d(_r) { (_r)->r = 1; (_r)->i = 0; }
#define set_one_mp(_r) { mpf_set_ui((_r)->r, 1); mpf_set_ui((_r)->i, 0); }
#define set_one_rat(_r) { mpq_set_ui((_r)[0], 1, 1); mpq_set_ui((_r)[1], 0, 1); }

// set neg one
#define set_neg_one_d(_r) { (_r)->r = -1; (_r)->i = 0; }
#define set_neg_one_mp(_r) { mpf_set_si((_r)->r, -1); mpf_set_ui((_r)->i, 0); }
#define set_neg_one_rat(_r) { mpq_set_si((_r)[0], -1, 1); mpq_set_ui((_r)[1], 0, 1); }

// conversion functions (d to mp & mp to d & rat to d & rat to mp)
#define d_to_mp(_m, _d)  { mpf_set_d((_m)->r, (_d)->r); mpf_set_d((_m)->i, (_d)->i); }
#define mp_to_d(_d, _m)  { (_d)->r = mpf_get_d((_m)->r); (_d)->i = mpf_get_d((_m)->i); }
#define mp_to_rat(_r, _m)  { mpf_t_to_rat((_r)[0], (_m)->r); mpf_t_to_rat((_r)[1], (_m)->i); }
#define d_to_rat(_r, _d)  { mpq_set_d((_r)[0], (_d)->r); mpq_set_d((_r)[1], (_d)->i); }
#define rat_to_d(_d, _r) { (_d)->r = mpq_get_d((_r)[0]); (_d)->i = mpq_get_d((_r)[1]); }
#define rat_to_mp(_m, _r){ mpf_set_q((_m)->r, (_r)[0]);  mpf_set_q((_m)->i, (_r)[1]); }

#define norm_sqr_d(_d) (_d)->r*(_d)->r + (_d)->i*(_d)->i
#define norm_sqr_mp(_n,_m) { int _oid = thread_num(); mpf_mul(_tempMPF1[_oid], (_m)->r, (_m)->r); mpf_mul(_n, (_m)->i, (_m)->i); mpf_add(_n, _n, _tempMPF1[_oid]); }

// clear
#define clear_d(_a)  { (_a)->r = (_a)->i = 0; }
#define clear_mp(_a) { mpf_clear((_a)->r); mpf_clear((_a)->i); }
#define clear_rat(_a) { mpq_clear((_a)[0]); mpq_clear((_a)[1]); }
#define clear_d_mp_rat(_d, _mp, _rat, _MPType) { \
if (_MPType == 0 || _MPType == 2) { clear_d(_d); } \
if (_MPType == 1 || _MPType == 2) { clear_mp(_mp); } \
if (_MPType == 2) { clear_rat(_rat); } }

////////////////////////////
// Basic vector operations - initialize, clear, etc..

// initialize
#define init_point_d(_a, _s) { (_a)->coord = (_comp_d *)bmalloc((_s) * sizeof(_comp_d)); (_a)->alloc_size = _s; (_a)->size = 0; }
#define init_point_mp2(_a, _s, _prec) { int _i; (_a)->coord = (_comp_mp *)bmalloc((_s) * sizeof(_comp_mp)); for (_i = 0; _i < _s; _i++) init_mp2(&(_a)->coord[_i], _prec); \
  (_a)->curr_prec = _prec; (_a)->alloc_size = _s; (_a)->size = 0; }
#define init_point_mp(_a, _s) { int _p = mpf_get_default_prec(); init_point_mp2(_a, _s, _p); }
#define init_point_rat(_a, _size) { int _i; (_a) = (mpq_t **)bmalloc(_size * sizeof(mpq_t *)); \
  for (_i = 0; _i < _size; _i++) { (_a)[_i] = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); init_rat((_a)[_i]); } }
#define init_vec_d   init_point_d
#define init_vec_mp  init_point_mp
#define init_vec_mp2 init_point_mp2
#define init_vec_rat init_point_rat

// increase size
#define increase_size_point_d(_a, _new_size) { if ((_a)->alloc_size < _new_size) { (_a)->coord = (_comp_d *)brealloc((_a)->coord, (_new_size) * sizeof(_comp_d)); (_a)->alloc_size = _new_size; }}
#define increase_size_point_mp(_a, _new_size) { if ((_a)->alloc_size < _new_size) { int _i; (_a)->coord = (_comp_mp *)brealloc((_a)->coord, (_new_size) * sizeof(_comp_mp)); \
  for (_i = (_a)->alloc_size; _i < _new_size; _i++) init_mp2(&(_a)->coord[_i], (_a)->curr_prec); (_a)->alloc_size = _new_size; }}
#define increase_size_vec_d increase_size_point_d
#define increase_size_vec_mp increase_size_point_mp

// change size
#define change_size_point_d(_a, _new_size) { if ((_a)->alloc_size != _new_size) { \
  (_a)->coord = (_comp_d *)brealloc((_a)->coord, (_new_size) * sizeof(_comp_d)); (_a)->size = (_a)->alloc_size = _new_size; }}
#define change_size_point_mp(_a, _new_size) { if ((_a)->alloc_size != _new_size) { \
  int _i; for (_i = (_a)->alloc_size - 1; _i >= 0; _i--) clear_mp(&(_a)->coord[_i]); \
  (_a)->coord = (_comp_mp *)brealloc((_a)->coord, (_new_size) * sizeof(_comp_mp)); \
  for (_i = _new_size - 1; _i >= 0; _i--) init_mp2(&(_a)->coord[_i], (_a)->curr_prec); (_a)->size = (_a)->alloc_size = _new_size; }}
#define change_size_vec_d change_size_point_d
#define change_size_vec_mp change_size_point_mp

// set precision
#define setprec_point_mp(_a, _prec) { if ((_a)->curr_prec != _prec) { int _i; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_size - 1; _i >= 0; _i--) setprec_mp(&(_a)->coord[_i], _prec); } }
#define setprec_vec_mp setprec_point_mp

// change precision
#define change_prec_point_mp(_a, _prec) { if ((_a)->curr_prec != _prec) { int _i; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_size - 1; _i >= 0; _i--) change_prec_mp(&(_a)->coord[_i], _prec); }}
#define change_prec_point_mp_rat(_a, _prec, _rat) { if ((_a)->curr_prec != _prec) { int _i; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_size - 1; _i >= 0; _i--) if (_i < (_a)->alloc_size) { change_prec_mp_rat(&(_a)->coord[_i], _prec, (_rat)[_i]); } \
    else { setprec_mp(&(_a)->coord[_i], _prec); }}}
#define change_prec_vec_mp change_prec_point_mp
#define change_prec_vec_mp_rat change_prec_point_mp_rat

// copy functions
#define point_cp_d(_d, _s) { int _i, _size = (_s)->size; increase_size_point_d(_d, _size); \
  (_d)->size = _size; for (_i = 0; _i < _size; _i++) set_d(&(_d)->coord[_i], &(_s)->coord[_i]); }
#define point_cp_mp(_d,_s) { int _i, _size = (_s)->size; increase_size_point_mp(_d, _size); \
  (_d)->size = _size; for (_i = 0; _i < _size; _i++) set_mp(&(_d)->coord[_i], &(_s)->coord[_i]); }
#define point_cp_rat(_d,_s,_r) { int _i; for (_i = 0; _i < _r; _i++) { mpq_set((_d)[_i][0], (_s)[_i][0]); mpq_set((_d)[_i][1], (_s)[_i][1]); } }

#define vec_cp_d point_cp_d
#define vec_cp_mp point_cp_mp
#define vec_cp_rat point_cp_rat

// conversion functions (d to mp & mp to d)
#define point_d_to_mp(_m, _d) { int _i, _size = (_d)->size; increase_size_point_mp(_m, _size); \
  (_m)->size = _size; for (_i = 0; _i < _size; _i++) d_to_mp(&(_m)->coord[_i], &(_d)->coord[_i]); }
#define point_d_to_rat(_r, _d) { int _i, _size = (_d)->size; for (_i = 0; _i < _size; _i++) d_to_rat((_r)[_i], &(_d)->coord[_i]); }
#define point_mp_to_d(_d, _m) { int _i, _size = (_m)->size; increase_size_point_d(_d, _size); \
  (_d)->size = _size; for (_i = 0; _i < _size; _i++) mp_to_d(&(_d)->coord[_i], &(_m)->coord[_i]); }
#define point_mp_to_rat(_r, _m) { int _i, _size = (_m)->size; for (_i = 0; _i < _size; _i++) mp_to_rat((_r)[_i], &(_m)->coord[_i]); }
#define point_rat_to_d(_d, _r, _size) { int _i; increase_size_point_d(_d, _size); (_d)->size = _size; \
  for (_i = 0; _i < _size; _i++) rat_to_d(&(_d)->coord[_i], (_r)[_i]); }
#define point_rat_to_mp(_m, _r, _size) { int _i; increase_size_point_mp(_m, _size); (_m)->size = _size; \
  for (_i = 0; _i < _size; _i++) rat_to_mp(&(_m)->coord[_i], (_r)[_i]); }

#define vec_d_to_mp point_d_to_mp
#define vec_d_to_rat point_d_to_rat
#define vec_mp_to_d point_mp_to_d
#define vec_mp_to_rat point_mp_to_rat
#define vec_rat_to_d point_rat_to_d
#define vec_rat_to_mp point_rat_to_mp

// scale vectors
#define vec_mulrdouble_d(_r, _x, _t) { int _i, _size = (_x)->size; increase_size_vec_d(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) mul_rdouble_d(&(_r)->coord[_i], &(_x)->coord[_i], _t); }
#define vec_mulrdouble_mp(_r, _x, _t) { int _i, _size = (_x)->size; increase_size_vec_mp(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) mul_rdouble_mp(&(_r)->coord[_i], &(_x)->coord[_i], _t); }
#define vec_mulrmpf_mp(_r, _x, _t) { int _i, _size = (_x)->size; increase_size_vec_mp(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) mul_rmpf_mp(&(_r)->coord[_i], &(_x)->coord[_i], _t); }
#define vec_mulcomp_d(_r, _x, _t) { int _i, _size = (_x)->size; increase_size_vec_d(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) mul_d(&(_r)->coord[_i], &(_x)->coord[_i], _t); }
#define vec_mulcomp_mp(_r, _x, _t) { int _i, _size = (_x)->size; increase_size_vec_mp(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) mul_mp(&(_r)->coord[_i], &(_x)->coord[_i], _t); }

// add vectors
#define vec_add_d(_r, _u, _v) { int _i, _size = (_u)->size; increase_size_vec_d(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) add_d(&(_r)->coord[_i], &(_u)->coord[_i], &(_v)->coord[_i]); }
#define vec_add_mp(_r, _u, _v) { int _i, _size = (_v)->size; increase_size_vec_mp(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) add_mp(&(_r)->coord[_i], &(_u)->coord[_i], &(_v)->coord[_i]); }

// subtract vector
#define vec_sub_d(_r, _u, _v) { int _i, _size = (_u)->size; increase_size_vec_d(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) sub_d(&(_r)->coord[_i], &(_u)->coord[_i], &(_v)->coord[_i]); }
#define vec_sub_mp(_r, _u, _v) { int _i, _size = (_u)->size; increase_size_vec_mp(_r, _size); \
  (_r)->size = _size; for (_i = 0; _i < _size; _i++) sub_mp(&(_r)->coord[_i], &(_u)->coord[_i], &(_v)->coord[_i]); }

// clear
#define clear_point_d(_a)  { free((_a)->coord); (_a)->coord = NULL; (_a)->alloc_size = (_a)->size = 0; }
#define clear_point_mp(_a) { int _i; for (_i = (_a)->alloc_size - 1; _i >= 0; _i--) clear_mp(&(_a)->coord[_i]); free((_a)->coord); (_a)->coord = NULL; (_a)->alloc_size = (_a)->size = 0; }
#define clear_point_rat(_a, _size) { int _i; for (_i = _size - 1; _i >= 0; _i--) { clear_rat((_a)[_i]); free((_a)[_i]); } free(_a); }
#define clear_point(_d, _mp, _rat, _MPType) { int _size = _MPType == 2 ? (_d)->size : 0; \
  if (_MPType == 0 || _MPType == 2) { clear_point_d(_d); } \
  if (_MPType == 1 || _MPType == 2) { clear_point_mp(_mp); } \
  if (_MPType == 2) { clear_point_rat(_rat, _size); } }

#define clear_vec_d   clear_point_d
#define clear_vec_mp  clear_point_mp
#define clear_vec_rat clear_point_rat
#define clear_vec clear_point

// initialize point_data_*
#define init_point_data_d(_pd, _s) { init_point_d((_pd)->point, _s); init_d((_pd)->time); (_pd)->cycle_num = 0; }
#define init_point_data_mp(_pd, _s) { init_point_mp((_pd)->point, _s); init_mp((_pd)->time); (_pd)->cycle_num = 0; }
#define init_point_data_mp2(_pd, _s, _prec) { init_point_mp2((_pd)->point, _s, _prec); init_mp2((_pd)->time, _prec); (_pd)->cycle_num = 0; }

// change precision
#define change_prec_point_data_mp(_pd, _prec) { change_prec_point_mp((_pd)->point, _prec); change_prec_mp((_pd)->time, _prec); }
#define change_prec_point_data_mp2(_pd, _prec) { change_prec_point_mp((_pd)->point, _prec); change_prec_mp2((_pd)->time, _prec); }

// set precision
#define setprec_point_data_mp(_pd, _prec) { setprec_point_mp((_pd)->point, _prec); setprec_mp((_pd)->time, _prec); }

// clear point_data_*
#define clear_point_data_d(_pd) { clear_point_d((_pd)->point); clear_d((_pd)->time); (_pd)->cycle_num = 0; }
#define clear_point_data_mp(_pd) { clear_point_mp((_pd)->point); clear_mp((_pd)->time); (_pd)->cycle_num = 0; }

// initialize endgame_data_t
#define init_endgame_data(_eg, _prec) { int _p = MAX(64, _prec); (_eg)->prec = (_eg)->last_approx_prec = _prec; \
  init_point_data_d(&(_eg)->PD_d, 0); init_point_data_mp2(&(_eg)->PD_mp, 0, _p); \
  init_point_d((_eg)->last_approx_d, 0); init_point_mp((_eg)->last_approx_mp, 0); \
  (_eg)->retVal = (_eg)->pathNum = (_eg)->codim = (_eg)->first_increase = (_eg)->condition_number = 0; \
  (_eg)->function_residual_d = (_eg)->latest_newton_residual_d = (_eg)->t_val_at_latest_sample_point_d = (_eg)->error_at_latest_sample_point_d = 0; \
  mpf_init2((_eg)->function_residual_mp, _p); mpf_init2((_eg)->latest_newton_residual_mp, _p); \
  mpf_init2((_eg)->t_val_at_latest_sample_point_mp, _p); mpf_init2((_eg)->error_at_latest_sample_point_mp, _p); }

// clear endgame_data_t
#define clear_endgame_data(_eg) { clear_point_data_d(&(_eg)->PD_d); clear_point_data_mp(&(_eg)->PD_mp); \
  clear_point_d((_eg)->last_approx_d); clear_point_mp((_eg)->last_approx_mp); \
  (_eg)->retVal = (_eg)->pathNum = (_eg)->codim = (_eg)->first_increase = (_eg)->condition_number = 0; \
  (_eg)->function_residual_d = (_eg)->latest_newton_residual_d = (_eg)->t_val_at_latest_sample_point_d = (_eg)->error_at_latest_sample_point_d = 0; \
  mpf_clear((_eg)->function_residual_mp); mpf_clear((_eg)->latest_newton_residual_mp); \
  mpf_clear((_eg)->t_val_at_latest_sample_point_mp); mpf_clear((_eg)->error_at_latest_sample_point_mp); }

// initialize trackBack_samples_t
#define init_trackBack_sample(_eg, _prec) { init_endgame_data(&(_eg)->endPt, _prec); \
  (_eg)->numSamples = (_eg)->samplePts_prec = (_eg)->midPt_prec = 0; \
  (_eg)->normSamples_d = NULL; (_eg)->normSamples_mp = NULL; (_eg)->midPt_d = (_eg)->samplePts_d = NULL; (_eg)->midPt_mp = (_eg)->samplePts_mp = NULL; \
  (_eg)->samplePts_notUsed = NULL; }

// clear trackBack_samples_t
#define clear_trackBack_sample(_tb) { int _j; clear_endgame_data(&(_tb)->endPt); \
  if ((_tb)->samplePts_prec < 64) { for (_j = (_tb)->numSamples - 1; _j >= 0; _j--) { clear_point_d((_tb)->samplePts_d[_j]); } \
    free((_tb)->samplePts_d); free((_tb)->normSamples_d); } \
  else { for (_j = (_tb)->numSamples - 1; _j >= 0; _j--) { mpf_clear((_tb)->normSamples_mp[_j]); clear_point_mp((_tb)->samplePts_mp[_j]); } \
    free((_tb)->samplePts_mp); free((_tb)->normSamples_mp); } \
  free((_tb)->samplePts_notUsed); \
  if ((_tb)->midPt_prec < 64) { for (_j = (_tb)->numSamples - 1; _j >= 0; _j--) { clear_point_d((_tb)->midPt_d[_j]); } free((_tb)->midPt_d); } \
  else { for (_j = (_tb)->numSamples - 1; _j >= 0; _j--) { clear_point_mp((_tb)->midPt_mp[_j]); }  free((_tb)->midPt_mp); } \
  (_tb)->numSamples = (_tb)->samplePts_prec = (_tb)->midPt_prec = 0; \
  (_tb)->normSamples_d = NULL; (_tb)->normSamples_mp = NULL; (_tb)->midPt_d = (_tb)->samplePts_d = NULL; (_tb)->midPt_mp = (_tb)->samplePts_mp = NULL; \
  (_tb)->samplePts_notUsed = NULL; }

/////////////////////////////////////
// Basic matrix operations - initialize, clear, etc..

// initialize
#define init_mat_d(_a, _r, _c) { int _i; (_a)->entry = (_comp_d **)bmalloc((_r) * sizeof(_comp_d *)); \
  for (_i = 0; _i < _r; _i++) (_a)->entry[_i] = (_comp_d *)bmalloc((_c) * sizeof(_comp_d)); (_a)->alloc_rows = _r; (_a)->alloc_cols = _c; (_a)->rows = (_a)->cols = 0; }
#define init_mat_mp2(_a, _r, _c, _prec) { int _i, _j; (_a)->entry = (_comp_mp **)bmalloc((_r) * sizeof(_comp_mp *)); \
  for (_i = 0; _i < _r; _i++) { (_a)->entry[_i] = (_comp_mp *)bmalloc((_c) * sizeof(_comp_mp)); for (_j = 0; _j < _c; _j++) init_mp2(&(_a)->entry[_i][_j], _prec); } \
  (_a)->curr_prec = _prec; (_a)->alloc_rows = _r; (_a)->alloc_cols = _c; (_a)->rows = (_a)->cols = 0; }
#define init_mat_mp(_a, _r, _c) { int _p = mpf_get_default_prec(); init_mat_mp2(_a, _r, _c, _p); }
#define init_mat_rat(_a, _rows, _cols) { int _i, _j; (_a) = (mpq_t ***)bmalloc((_rows) * sizeof(mpq_t **)); \
  for (_i = 0; _i < _rows; _i++) { (_a)[_i] = (mpq_t **)bmalloc((_cols) * sizeof(mpq_t *)); \
  for (_j = 0; _j < _cols; _j++) { (_a)[_i][_j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); init_rat((_a)[_i][_j]); } } }

// increase rows
#define increase_rows_mat_d(_a, _new_rows) { if ((_a)->alloc_rows < _new_rows) { int _i; \
  (_a)->entry = (_comp_d **)brealloc((_a)->entry, (_new_rows) * sizeof(_comp_d *)); \
  for (_i = (_a)->alloc_rows; _i < _new_rows; _i++) (_a)->entry[_i] = (_comp_d *)bmalloc((_a)->alloc_cols * sizeof(_comp_d)); (_a)->alloc_rows = _new_rows; }}
#define increase_rows_mat_mp(_a, _new_rows) { if ((_a)->alloc_rows < _new_rows) { int _i, _j; \
  (_a)->entry = (_comp_mp **)brealloc((_a)->entry, (_new_rows) * sizeof(_comp_mp *)); \
  for (_i = (_a)->alloc_rows; _i < _new_rows; _i++) { (_a)->entry[_i] = (_comp_mp *)bmalloc((_a)->alloc_cols * sizeof(_comp_mp)); \
    for (_j = 0; _j < (_a)->alloc_cols; _j++) init_mp2(&(_a)->entry[_i][_j], (_a)->curr_prec); } (_a)->alloc_rows = _new_rows; }}

// increase cols
#define increase_cols_mat_d(_a, _new_cols) { if ((_a)->alloc_cols < _new_cols) { int _i; for (_i = 0; _i < (_a)->alloc_rows; _i++) \
  (_a)->entry[_i] = (_comp_d *)brealloc((_a)->entry[_i], (_new_cols) * sizeof(_comp_d)); (_a)->alloc_cols = _new_cols; }}
#define increase_cols_mat_mp(_a, _new_cols) { if ((_a)->alloc_cols < _new_cols) { int _i,_j; for (_i = 0; _i < (_a)->alloc_rows; _i++) { \
  (_a)->entry[_i] = (_comp_mp *)brealloc((_a)->entry[_i], (_new_cols) * sizeof(_comp_mp)); for (_j = (_a)->alloc_cols; _j < _new_cols; _j++) \
  init_mp2(&(_a)->entry[_i][_j], (_a)->curr_prec); } (_a)->alloc_cols = _new_cols; }}

// increase size
#define increase_size_mat_d(_a, _new_rows, _new_cols) { increase_rows_mat_d(_a, _new_rows); increase_cols_mat_d(_a, _new_cols); }
#define increase_size_mat_mp(_a, _new_rows, _new_cols) { increase_rows_mat_mp(_a, _new_rows); increase_cols_mat_mp(_a, _new_cols); }

// decrease rows
#define decrease_rows_mat_d(_a, _new_rows) { if ((_a)->alloc_rows > _new_rows) { int _i; /* clear extra rows */ \
  for (_i = _new_rows; _i < (_a)->alloc_rows; _i++) if ((_a)->entry[_i] != NULL) free((_a)->entry[_i]); \
  (_a)->entry = (_comp_d **)brealloc((_a)->entry, (_new_rows) * sizeof(_comp_d *)); (_a)->alloc_rows = _new_rows; }}
#define decrease_rows_mat_mp(_a, _new_rows) { if ((_a)->alloc_rows > _new_rows) { int _i,_j; /* clear extra rows */ \
  for (_i = _new_rows; _i < (_a)->alloc_rows; _i++) if ((_a)->entry[_i] != NULL) { \
    for (_j = 0; _j < (_a)->alloc_cols; _j++) clear_mp(&(_a)->entry[_i][_j]); free((_a)->entry[_i]); } \
  (_a)->entry = (_comp_mp **)brealloc((_a)->entry, (_new_rows) * sizeof(_comp_mp *)); (_a)->alloc_rows = _new_rows; }}

// decrease cols
#define decrease_cols_mat_d(_a, _new_cols) { if ((_a)->alloc_cols > _new_cols) { int _i; /* clear extra cols */ \
  for (_i = 0; _i < (_a)->alloc_rows; _i++) (_a)->entry[_i] = (_comp_d *)brealloc((_a)->entry[_i], (_new_cols) * sizeof(_comp_d)); \
  (_a)->alloc_cols = _new_cols; }}
#define decrease_cols_mat_mp(_a, _new_cols) { if ((_a)->alloc_cols > _new_cols) { int _i,_j; /* clear extra cols */ \
  for (_i = 0; _i < (_a)->alloc_rows; _i++) { for (_j = _new_cols; _j < (_a)->alloc_cols; _j++) clear_mp(&(_a)->entry[_i][_j]); \
    (_a)->entry[_i] = (_comp_mp *)brealloc((_a)->entry[_i], (_new_cols) * sizeof(_comp_mp)); } (_a)->alloc_cols = _new_cols; }}

// change size
#define change_size_mat_d(_a, _new_rows, _new_cols) { \
  if ((_a)->alloc_rows < _new_rows) { /* increase the rows first and then change the size of the columns */ \
    increase_rows_mat_d(_a, _new_rows); \
    if ((_a)->alloc_cols < _new_cols) { increase_cols_mat_d(_a, _new_cols); } \
    else if ((_a)->alloc_cols > _new_cols) { decrease_cols_mat_d(_a, _new_cols); }} \
  else { /* change the size of the columns first and then change the rows */ \
    if ((_a)->alloc_cols < _new_cols) { increase_cols_mat_d(_a, _new_cols); } \
    else if ((_a)->alloc_cols > _new_cols) { decrease_cols_mat_d(_a, _new_cols); }\
    if ((_a)->alloc_rows > _new_rows) { decrease_rows_mat_d(_a, _new_rows); }}}
  
#define change_size_mat_mp(_a, _new_rows, _new_cols) { \
  if ((_a)->alloc_rows < _new_rows) { /* increase the rows first and then change the size of the columns */ \
    increase_rows_mat_mp(_a, _new_rows); \
    if ((_a)->alloc_cols < _new_cols) { increase_cols_mat_mp(_a, _new_cols); } \
    else if ((_a)->alloc_cols > _new_cols) { decrease_cols_mat_mp(_a, _new_cols); }} \
  else { /* change the size of the columns first and then change the rows */ \
    if ((_a)->alloc_cols < _new_cols) { increase_cols_mat_mp(_a, _new_cols); } \
    else if ((_a)->alloc_cols > _new_cols) { decrease_cols_mat_mp(_a, _new_cols); }\
    if ((_a)->alloc_rows > _new_rows) { decrease_rows_mat_mp(_a, _new_rows); }}}

// set precision
#define setprec_mat_mp(_a, _prec) { if ((_a)->curr_prec != _prec) { int _i, _j; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) for (_j = (_a)->alloc_cols - 1; _j >= 0; _j--) setprec_mp(&(_a)->entry[_i][_j], _prec); } }

// change precision
#define change_prec_mat_mp(_a, _prec) { if ((_a)->curr_prec != _prec) { int _i, _j; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) for (_j = (_a)->alloc_cols - 1; _j >= 0; _j--) change_prec_mp(&(_a)->entry[_i][_j], _prec); }}
#define change_prec_mat_mp_rat(_a, _prec, _rat) { if ((_a)->curr_prec != _prec) { int _i,_j; (_a)->curr_prec = _prec; \
  for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) for (_j = (_a)->alloc_cols - 1; _j >= 0; _j--) \
    if (_i < (_a)->rows && _j < (_a)->cols) { change_prec_mp_rat(&(_a)->entry[_i][_j], _prec, (_rat)[_i][_j]); } \
    else { setprec_mp(&(_a)->entry[_i][_j], _prec); }}} 

// copy functions
#define mat_cp_d(_d,_s)  { int _i, _j, _rows = (_s)->rows, _cols = (_s)->cols; increase_size_mat_d(_d, _rows, _cols); \
  (_d)->rows = _rows; (_d)->cols = _cols; for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) set_d(&(_d)->entry[_i][_j], &(_s)->entry[_i][_j]); }
#define mat_cp_mp(_d,_s) { int _i, _j, _rows = (_s)->rows, _cols = (_s)->cols; increase_size_mat_mp(_d, _rows, _cols); \
  (_d)->rows = _rows; (_d)->cols = _cols; for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) set_mp(&(_d)->entry[_i][_j], &(_s)->entry[_i][_j]); }
#define mat_cp_rat(_A, _B, _r, _c) { int _i, _j; for (_i = _r - 1; _i >= 0; _i--) for (_j = _c - 1; _j >= 0; _j--) \
    { mpq_set((_A)[_i][_j][0], (_B)[_i][_j][0]); mpq_set((_A)[_i][_j][1], (_B)[_i][_j][1]); } }
#define mat_cp_mp_rat(_A, _A_rat, _B, _B_rat) { int _i, _j; _A->rows = _B->rows; _A->cols = _B->cols; \
  for (_i = (_B)->rows - 1; _i >= 0; _i--) for (_j = (_B)->cols - 1; _j >= 0; _j--) \
    { set_mp(&(_A)->entry[_i][_j], &(_B)->entry[_i][_j]); \
      mpq_set((_A_rat)[_i][_j][0], (_B_rat)[_i][_j][0]); mpq_set((_A_rat)[_i][_j][1], (_B_rat)[_i][_j][1]); } }

// conversion functions (d to mp & mp to d)
#define mat_d_to_mp(_m, _d) { int _i, _j, _r = (_d)->rows, _c = (_d)->cols; increase_size_mat_mp(_m, _r, _c); \
  (_m)->rows = _r; (_m)->cols = _c; for (_i = 0; _i < _r; _i++) for (_j = 0; _j < _c; _j++) d_to_mp(&(_m)->entry[_i][_j], &(_d)->entry[_i][_j]); }
#define mat_d_to_rat(_r, _d) { int _i, _j, _rows = (_d)->rows, _cols = (_d)->cols; \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) d_to_rat((_r)[_i][_j], &(_d)->entry[_i][_j]); }
#define mat_mp_to_d(_d, _m) { int _i, _j, _r = (_m)->rows, _c = (_m)->cols; increase_size_mat_d(_d, _r, _c); \
  (_d)->rows = _r; (_d)->cols = _c; for (_i = 0; _i < _r; _i++) for (_j = 0; _j < _c; _j++) mp_to_d(&(_d)->entry[_i][_j], &(_m)->entry[_i][_j]); }
#define mat_mp_to_rat(_r, _m) { int _i, _j, _rows = (_m)->rows, _cols = (_m)->cols; \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) mp_to_rat((_r)[_i][_j], &(_m)->entry[_i][_j]); }
#define mat_rat_to_d(_d, _r, _rows, _cols) { int _i, _j; increase_size_mat_d(_d, _rows, _cols); (_d)->rows = _rows; (_d)->cols = _cols; \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) rat_to_d(&(_d)->entry[_i][_j], (_r)[_i][_j]); }
#define mat_rat_to_mp(_m, _r, _rows, _cols) { int _i, _j; increase_size_mat_mp(_m, _rows, _cols); (_m)->rows = _rows; (_m)->cols = _cols; \
  for (_i = 0; _i < _rows; _i++) for (_j = 0; _j < _cols; _j++) rat_to_mp(&(_m)->entry[_i][_j], (_r)[_i][_j]); }

// negate
#define mat_neg_d(_m) { int _i, _j; for (_i = 0; _i < (_m)->rows; _i++) for (_j = 0; _j < (_m)->cols; _j++) neg_d(&(_m)->entry[_i][_j]); }
#define mat_neg_mp(_m) { int _i, _j; for (_i = 0; _i < (_m)->rows; _i++) for (_j = 0; _j < (_m)->cols; _j++) neg_mp(&(_m)->entry[_i][_j]); }

// clear
#define clear_mat_d(_a)  { int _i; for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) { free((_a)->entry[_i]); } free((_a)->entry); (_a)->entry = NULL; \
  (_a)->alloc_rows = (_a)->alloc_cols = (_a)->rows = (_a)->cols = 0; }
#define clear_mat_mp(_a) { int _i, _j; for (_i = (_a)->alloc_rows - 1; _i >= 0; _i--) { for (_j = (_a)->alloc_cols - 1; _j >= 0; _j--) clear_mp(&(_a)->entry[_i][_j]); \
  free((_a)->entry[_i]); } free((_a)->entry); (_a)->entry = NULL; (_a)->alloc_rows = (_a)->alloc_cols = (_a)->rows = (_a)->cols = 0; }
#define clear_mat_rat(_a, _rows, _cols) { int _i, _j; \
    for (_i = _rows - 1; _i >= 0; _i--) { for (_j = _cols - 1; _j >= 0; _j--) { \
    clear_rat((_a)[_i][_j]); free((_a)[_i][_j]); } free((_a)[_i]); } free(_a); }
#define clear_mat(_d, _mp, _rat, _MPType) { \
  int _rows = _MPType == 2 ? (_d)->rows : 0, _cols = _MPType == 2 ? (_d)->cols : 0; \
  if (_MPType == 0 || _MPType == 2) { clear_mat_d(_d); } \
  if (_MPType == 1 || _MPType == 2) { clear_mat_mp(_mp); } \
  if (_MPType == 2) { clear_mat_rat(_rat, _rows, _cols); } }

// eval struct
#define init_eval_struct_d(_e, _num_funcs, _num_vars, _num_params) { init_mat_d((_e).Jv, _num_funcs, _num_vars); init_mat_d((_e).Jp, _num_funcs, _num_params); \
  init_vec_d((_e).funcVals, _num_funcs); init_vec_d((_e).parVals, _num_params); init_vec_d((_e).parDer, _num_params); }
#define init_eval_struct_mp(_e, _num_funcs, _num_vars, _num_params) { init_mat_mp((_e).Jv, _num_funcs, _num_vars); init_mat_mp((_e).Jp, _num_funcs, _num_params); \
  init_vec_mp((_e).funcVals, _num_funcs); init_vec_mp((_e).parVals, _num_params); init_vec_mp((_e).parDer, _num_params); }
#define init_eval_struct_mp2(_e, _num_funcs, _num_vars, _num_params, _pr) { init_mat_mp2((_e).Jv, _num_funcs, _num_vars, _pr); init_mat_mp2((_e).Jp, _num_funcs, _num_params, _pr); \
  init_vec_mp2((_e).funcVals, _num_funcs, _pr); init_vec_mp2((_e).parVals, _num_params, _pr); init_vec_mp2((_e).parDer, _num_params, _pr); }
#define setprec_eval_struct_mp(_e, _pr) { setprec_mat_mp((_e).Jv, _pr); setprec_mat_mp((_e).Jp, _pr); \
  setprec_vec_mp((_e).funcVals, _pr); setprec_vec_mp((_e).parVals, _pr); setprec_vec_mp((_e).parDer, _pr); }
#define clear_eval_struct_d(_e) { clear_mat_d((_e).Jv); clear_mat_d((_e).Jp); clear_vec_d((_e).funcVals); clear_vec_d((_e).parVals); clear_vec_d((_e).parDer); }
#define clear_eval_struct_mp(_e) { clear_mat_mp((_e).Jv); clear_mat_mp((_e).Jp); clear_vec_mp((_e).funcVals); clear_vec_mp((_e).parVals); clear_vec_mp((_e).parDer); }

// variables used in basic multiprecsion arithmetic defined above which are actually allocated & initialized in misc.c
extern mpf_t *_tempMPF1, *_tempMPF2, *_tempMPF3;

// adaptiveMP2.c
int AMPtrack(point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, comp_d fT_d, comp_mp fT_mp, int time_prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

int refine_d_basic(point_data_d *P, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int refine_mp_basic(point_data_mp *P, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int refine_digits_amp(int outputLevel, int digitsCorrect, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, int digitsCorrect_in, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, comp_d time, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));
int refine_digits_d(int outputLevel, int digitsCorrect, double *latest_newton_residual_d, int digitsCorrect_in, point_data_d *out, point_data_d *in, comp_d time, FILE *OUT, eval_struct_d *e_d, void const *ED_d, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int refine_digits_mp(int outputLevel, int digitsCorrect, mpf_t latest_newton_residual_mp, int digitsCorrect_in, point_data_mp *out, point_data_mp *in, int prec_in, comp_mp time, FILE *OUT, eval_struct_mp *e_mp, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void AMP2_update(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, int currPrec, double P0, int digits_C, double newtonTol, int maxNewtonIts, double finalTol, int checkFinalTol, double eta_minSS, double eta_maxSS, comp_d currTime_d, comp_mp currTime_mp);
int newton_d_amp(double *norm_J, double *norm_J_inv, int num_digits, point_d P, comp_d t, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int newton_mp_amp(double *norm_J, double *norm_J_inv, int num_digits, point_mp P, comp_mp t, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int minDigits_from_eta(double eta, comp_d currTime_d, comp_mp currTime_mp, int prec_in);
double AMP2_cost(int prec, double eta);
int AMP2_convergence_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T);
int AMP2_convergence_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T);

double minTime_d(int digits, comp_d t_d, comp_mp t_mp, int prec_in);
void minTime_mp(mpf_t minPrecTime, int digits, comp_d t_d, comp_mp t_mp, int prec_in);

// adaptiveMP_error.c
int AMPtrack_error(point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, comp_d fT_d, comp_mp fT_mp, int time_prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

// other 
void setup_omp_file(FILE ***Fptr, FILE *Orig, char *name, int max);
void combine_omp_file(FILE *F, FILE ***Fptr, int max);
void delete_omp_file(int max_threads, char *fName);
void setup_basic_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, tracker_config_t **T_copy, tracker_config_t *T);
void clear_basic_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *midName, char *failName, tracker_config_t **T_copy);

int find_Cauchy_double_samples_amp(point_data_d *Samples_d, point_data_mp *Samples_mp, int *samples_prec, point_data_d *samplePts_d, point_data_mp *samplePts_mp, int prec_in, int cycle_num, int samples_per_loop, double *time_first_increase, tracker_config_t *T, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

// QR_block.c
int QR_R_block_d(mat_d R, mat_d P, mat_d A, int rowTop, int colTop, double tol_pivot, double tol_sign, double largeChange);
int QLP_block_pair_d(mat_d L0, mat_d L1, mat_d A0, mat_d A1, int rowTop, int colTop, double tol_pivot, double tol_sign, double largeChange);
int QLP_block_pair_mp_prec(mat_mp L0, mat_mp L1, mat_mp A0, mat_mp A1, int rowTop, int colTop, int curr_prec);
int QLP_block_pair_mp(mat_mp L0, mat_mp L1, mat_mp A0, mat_mp A1, int rowTop, int colTop, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange);
int QLP_block_pair_amp_mp_d_prec(mat_mp L0, mat_d L1, mat_mp A0, int prec0, mat_d A1, int rowTop, int colTop);
int QLP_block_pair_amp_mp_d(mat_mp L0, mat_d L1, mat_mp A0, int prec0, mat_d A1, int rowTop, int colTop, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange);

void setup_frobenius_d(double **FV, mat_d m0, mat_d m1);

int isUnary(char op);

// splitParse.l
int splitParse(FILE *fp, char *funcName, char *configName); // split an input file into the function and configurations

#endif

