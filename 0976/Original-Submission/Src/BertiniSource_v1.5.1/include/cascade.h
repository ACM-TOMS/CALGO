// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _CASCADE_H
#define _CASCADE_H

#include "bertini.h"
#include "regeneration.h"
#include "eqbyeq.h"

typedef struct 
{
  int num_patches;
  mat_d patchCoeff; // of size num_patches x num_variables
} patch_eval_data_d;  // the structure to evaluate the patch equation in terms of the y variables, where (By = x)
                      // this will allow us to use the y input variables to easy calculate the patches and read off the Jacobian w.r.t. y

typedef struct
{
  int num_patches;
  mat_mp patchCoeff; // of size num_patches x num_variables
  mpq_t ***patchCoeff_rat; // num_patches x num_variables x {real, imaginary}
  int curr_prec;
} patch_eval_data_mp;  // the structure to evaluate the patch equation in terms of the y variables, where (By = x)
                       // this will allow us to use the y input variables to easy calculate the patches and read off the Jacobian w.r.t. y

typedef struct 
{
  int startSystemType; // 0 if total degree, 1 if mhom
  int size_r;     // the number of entries in *degrees
  int *degrees;   // the "total" degree of each function - if mhom, we store the sum of the mhom degrees - this will allow
                  // the multiplication of the linears described by coeff to be easy
  int max_degree; // the maxium "total" degree - this is NEEDED - we allocate memory based on this number!
  int coeff_cols; // the number of columns in coeff
  comp_d **coeff; // [linear number][variable number], where the linear number goes from 0 to sum(degrees[i], i = 0,..,num_funcs)
                  // this will be setup properly so that we can simply do multiplication by the vector of variables
                  // and have the linears be in the correct hom var group.
  comp_d gamma;
} start_system_eval_data_d; // this has the structure to store the degree of the functions in F to create the total degree start system

typedef struct
{
  int startSystemType; // 0 if total degree, 1 if mhom
  int size_r;     // the number of entries in *degrees
  int *degrees;   // the "total" degree of each function - if mhom, we store the sum of the mhom degrees - this will allow
                  // the multiplication of the linears described by coeff to be easy
  int max_degree; // the maxium "total" degree - this is NEEDED - we allocate memory based on this number!
  int coeff_cols; // the number of columns in coeff
  comp_mp **coeff; // [linear number][variable number], where the linear number goes from 0 to sum(degrees[i], i = 0,..,num_funcs)
                  // this will be setup properly so that we can simply do multiplication by the vector of variables
                  // and have the linears be in the correct hom var group.
  mpq_t ***coeff_rat;
  comp_mp gamma;  // the random gamma for the homotopy
  mpq_t *gamma_rat;
  int curr_prec;
} start_system_eval_data_mp; // this has the structure to store the degree of the functions in F to create the total degree start system

typedef struct
{ // main idea: we want to evaluate F = [I A.*W][P]f(By), where f is described in the straight-line program
  prog_t *Prog;       // s.l.p. to evaluate the original function f
  int size_f;         // the number of functions in f
  int *orig_degrees;  // the degrees of the functions in f
  int *new_degrees;   // the degrees of F - the square system
  mat_d B;            // the matrix that converts the input y variables to x variables that Prog needs (x = By)
  mat_d B_perp;       // orthogonal complement of B
  int noChanges;      // == 0 if the square system is something other then Prog, otherwise the square system is just Prog
  int *P;             // the vector that describes the permuation of the f_i's
  int **W;            // the matrix that contains the exponents of the homogenizing variable to maintain the 1-homogenity of the square system F
  int max_of_W;       // the maximum entry of W - this is not implemented now, but could be later for efficiency
  mat_d A;            // the matrix of random numbers used in creating the random linear combinations in F
  int size_r;         // the size of F
} square_system_eval_data_d;  // this is the structure to evaluate the square system F

typedef struct
{ // main idea: we want to evaluate F = [I A.*W][P]f(By), where f is described in the straight-line program
  prog_t *Prog;       // s.l.p. to evaluate the original function f 
  int size_f;         // the number of functions in f
  int *orig_degrees;  // the degrees of the functions in f
  int *new_degrees;   // the degrees of F - the square system
  mat_mp B;           // the matrix that converts the input y variables to x variables that Prog needs (x = By)
  mat_mp B_perp;      // orthogonal complement of B
  mpq_t ***B_rat;
  mpq_t ***B_perp_rat;
  int noChanges;      // == 0 if the square system is something other then Prog, otherwise the square system is just Prog
  int *P;             // the vector that describes the permuation of the f_i's
  int **W;            // the matrix that contains the exponents of the homogenizing variable to maintain the 1-homogenity of the square system F
  int max_of_W;       // the maximum entry of W - this is not implemented now, but could be later for efficiency
  mat_mp A;           // the matrix of random numbers used in creating the random linear combinations in F
  mpq_t ***A_rat;
  int size_r;         // the size of F
  int curr_prec;
} square_system_eval_data_mp;  // this is the structure to evaluate the square system F

typedef struct
{
  square_system_eval_data_mp squareSystem;
  patch_eval_data_mp patch;
  start_system_eval_data_mp startSystem;
  preproc_data preProcData;
  eqData_t *EqD;   // for equation-by-equation
} basic_eval_data_mp;

typedef struct 
{
  square_system_eval_data_d squareSystem;
  patch_eval_data_d patch;
  start_system_eval_data_d startSystem;
  preproc_data preProcData;
  basic_eval_data_mp *BED_mp; // used only for AMP
  eqData_t *EqD;   // for equation-by-equation
} basic_eval_data_d;

typedef struct
{
  int path_num;
  double norm;
} sortStruct_d;

typedef struct
{
  int path_num;
  mpf_t norm;
} sortStruct_mp;

// declarations from eval_functions.c
int eval_d(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int eval_mp(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int patch_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int patch_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int start_system_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int start_system_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int square_system_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int square_system_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int basic_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int basic_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int userHom_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int userHom_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int paramHom_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int paramHom_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

// declarations from setup_functions.c
int change_basic_eval_prec(void const *, int);
int change_square_prec(void const *Sq, int prec);

void setupTD_startPoints_d(char pointsIN[], char pointsOUT[], int size, int *degs, patch_eval_data_d *PED);
void setupTD_startPoints_mp(char pointsIN[], char pointsOUT[], int size, int *degs, patch_eval_data_mp *PED);

void setupPreProcData(char preprocFile[], preproc_data *PPD);

void setupSquareSystem_d(prog_t *Prog, int finalSize, preproc_data *PPD, char degreeFile[], square_system_eval_data_d *SSED, int adjustDegrees);
void setupSquareSystem_d_to_mp(square_system_eval_data_d *SSED_d, square_system_eval_data_mp *SSED_mp, int digits, int prec);
void setupSquareSystem_mp(prog_t *Prog, int finalSize, preproc_data *PPD, char degreeFile[], square_system_eval_data_mp *SSED, int digits, int prec, int adjustDegrees);

void setupCoeffInSS_d(char degreeFile[], int *P, preproc_data *PPD, start_system_eval_data_d *SSED);
void setupCoeffInSS_mp(char degreeFile[], int *P, preproc_data *PPD, start_system_eval_data_mp *SSED, int digits, int prec);

void setupStartSystem_d(int SSType, int size, int *deg, int *P, preproc_data *PPD, char degreeFile[], start_system_eval_data_d *SSED);
void setupStartSystem_d_to_mp(start_system_eval_data_d *SSED_d, start_system_eval_data_mp *SSED_mp, int digits, int prec);
void setupStartSystem_mp(int SSType, int size, int *deg, int *P, preproc_data *PPD, char degreeFile[], start_system_eval_data_mp *SSED, int digits, int prec);

void setupPatch_d(int PatchType, patch_eval_data_d *PED, void const *ptr1, void const *ptr2);
void setupPatch_d_to_mp(patch_eval_data_d *PED_d, patch_eval_data_mp *PED_mp, int digits, int prec, int patchType, int numVars, void const *ptr1, void const *ptr2);
void setupPatch_mp(int PatchType, patch_eval_data_mp *PED, int digits, int prec, void const *ptr1, void const *ptr2);
void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED);

void setupBasicEval_d(char preprocFile[], char degreeFile[], prog_t *dummyProg, int squareSize, int patchType, int ssType, int trackType, void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4, basic_eval_data_d *BED, int adjustDegrees);
void setupBasicEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg, int squareSize, int patchType, int ssType, int prec, void const *ptr1, void const *ptr2, basic_eval_data_mp *BED, int adjustDegrees);
void changeBasicEvalPrec_mp(int new_prec, basic_eval_data_mp *BED);

void setupBEDUsingUserHom_d(prog_t *dummyProg, int MPType, basic_eval_data_d *BED);
void setupBEDUsingUserHom_mp(prog_t *dummyProg, int MPType, basic_eval_data_mp *BED);

void setupBEDUsingParamHom_d(prog_t *dummyProg, char *preprocFile, int MPType, basic_eval_data_d *BED);
void setupBEDUsingParamHom_mp(prog_t *dummyProg, char *preprocFile, int MPType, basic_eval_data_mp *BED);

// declarations from startpoint_maker.c
int checkLoc(int *loc, int *P, int **mhomDeg, preproc_data *PPD);

void TDstartMaker_d(int *degs, int num_funcs);
void TDstartMaker_mp(int *degs, int num_funcs);

void MHstartMaker_d(preproc_data *PPD, int *P, comp_d **coeff, mat_d patchCoeff);
void MHstartMaker_mp(preproc_data *PPD, int *P, comp_mp **coeff, mat_mp patchCoeff);

/* Prototypes from rank_finder.c: */
int rank_finder_d(preproc_data *PPD, prog_t *orig_function, tracker_config_t *T, int num_variables);
int rank_finder_mp(preproc_data *PPD, prog_t *orig_function, tracker_config_t *T, int num_variables);

/* Prototypes from post_processing.c: */
void sort_points(int num_crossings, int *convergence_failures, int *sharpening_failures, int *sharpening_singular, char *inputName, int num_sols, int num_vars, double midpoint_tol, double final_tol, tracker_config_t *T, preproc_data *PPD, int regenToggle, int eqbyeqToggle);
void midpoint_checker(int num_paths, int num_vars, double tol, int *num_crossings);
void printFailureSummary(trackingStats *tC, int convergence_failures, int sharpening_failures, int sharpening_singular);
void printSharpeningFailureSummary(int total, int sharpening_failures, int sharpening_singular, int finiteSuccess, int infiniteSuccess, int nonfiniteSuccess);
int setupPostProcess(int *orig_prec, FILE *IN, post_process_t *endPoint, int size, int MPType);
void zeroDimPostProcess(FILE *OUT, post_process_t *endpoints, int num_sols, int num_vars, double final_tol, tracker_config_t *T, preproc_data *PPD, int num_crossings, int convergence_failures, char *inputName, int regenToggle, int eqbyeqToggle);
int checkForReal_d(point_d Pt, double realTol);
int checkForReal_mp(point_mp Pt, double realTol);
void findFiniteSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxNorm);
void findRealSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double realTol);
void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle);
void findMultSol(post_process_t *endPoints, int num_sols, int num_vars, preproc_data *PPD, double finalTol);

// zero_dim_main.c
void getDehomPoint_d(point_d dehomPoint, point_d inPoint, int num_vars, preproc_data *PPD);
void getDehomPoint_mp(point_mp dehomPoint, point_mp inPoint, int num_vars, preproc_data *PPD);

// zero_dim_setup.c
int userHom_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int userHom_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int paramHom_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, int moveStartPts, char *startName);
int paramHom_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, int moveStartPts, char *startName);

int zero_dim_basic_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, int findStartPts, char *pointsIN, char *pointsOUT);

int zero_dim_basic_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, int findStartPts, char *pointsIN, char *pointsOUT);

void printZeroDimRelevantData(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);
void updateZeroDimRelevantData(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP);

// zero_dim_track.c
void printSuccess_d(FILE *RAWOUT, point_d orig_vars, int cycle_num, int path_num, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, int success_flag);
void printSuccess_mp(FILE *RAWOUT, point_mp orig_vars, int cycle_num, int path_num, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, int success_flag);

void zero_dim_track_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void zero_dim_track_path_rank_d(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void zero_dim_track_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, basic_eval_data_mp *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void zero_dim_track_path_rank_mp(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, endgame_data_t *EG_out, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void setup_zero_dim_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_d **BED_copy, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp);

void setup_zero_dim_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, tracker_config_t *T, basic_eval_data_mp **BED_copy, basic_eval_data_mp *ED);

void clear_zero_dim_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, basic_eval_data_d **BED_copy);

void clear_zero_dim_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, basic_eval_data_mp **BED_copy);

void zero_dim_track_path_d(int pathNum, endgame_data_t *EG_out, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void zero_dim_track_path_mp(int pathNum, endgame_data_t *EG_out, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

void printPathHeader_d(FILE *OUT, point_data_d *PD, tracker_config_t *T, int pathNum, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

void printPathHeader_mp(FILE *OUT, point_data_mp *PD, tracker_config_t *T, int pathNum, void const *ED, int (*eval_f)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void printPathFooter_d(trackingStats *trackCount, endgame_data_t *EG, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, void const *ED);

void printPathFooter_mp(trackingStats *trackCount, endgame_data_t *EG, tracker_config_t *T, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, void const *ED);

void printPathFooterOut_d(FILE *OUT, FILE *RAWOUT, int success, int pathNum, point_data_d *PD, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, point_d dehomP, tracker_config_t *T, prog_t *Prog, int print_dehom, int eval_orig_f);

void printBasicFooter_d(FILE *OUT, point_data_d *PD, tracker_config_t *T, double func_residual);

void printPathFooterOut_mp(FILE *OUT, FILE *RAWOUT, int success, int pathNum, point_data_mp *PD, double cond_num, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, double first_increase, point_mp dehomP, tracker_config_t *T, prog_t *Prog, int print_dehom, int eval_orig_f);

void printBasicFooter_mp(FILE *OUT, point_data_mp *PD, tracker_config_t *T, mpf_t func_residual);

void printFailureMsg_d(FILE *FAIL, point_data_d *PD, point_d dehomP, int pathNum, int retVal, int isNumber, int isJunk, trackingStats *trackCount, tracker_config_t *T);
void printFailureMsg_mp(FILE *FAIL, point_data_mp *PD, point_mp dehomP, int pathNum, int retVal, int isNumber, int isJunk, trackingStats *trackCount, tracker_config_t *T);

void printResultOfPath(FILE *OUT, int retVal, tracker_config_t *T);

// dehomogenize the zero-dim point
int zero_dim_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
// set to size 0
int zero_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// copy_functions.c

void cp_tracker_config_t(tracker_config_t *T, const tracker_config_t *T_input);
void cp_preproc_data(preproc_data *PPD, preproc_data *PPD_input);
void cp_prog_t(prog_t *Prog, prog_t *Prog_input);

void cp_basic_eval_data_d(basic_eval_data_d *BED, basic_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType);
void cp_square_system_d(square_system_eval_data_d *SSED, square_system_eval_data_d *SSED_input);
void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d *PED_input);
void cp_start_system_d(start_system_eval_data_d *SSED, start_system_eval_data_d *SSED_input);

void cp_basic_eval_data_mp(basic_eval_data_mp *BED, basic_eval_data_mp *BED_input);
void cp_square_system_mp(square_system_eval_data_mp *SSED, square_system_eval_data_mp *SSED_input, int copy_SLP, prog_t *Prog_ptr);
void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp *PED_input);
void cp_start_system_mp(start_system_eval_data_mp *SSED, start_system_eval_data_mp *SSED_input);

void cp_comp_mp_int(void *Out, void *In, char **str, int freeStr, int inType);
void cp_comp_rat_int(void *Out, void *In, char **str, int freeStr, int inType);

void cp_point_d_int(void *Out, void *In, comp_d **coeff, int freeCoeff, int initPoint, int inType);
void cp_point_data_d_int(void *PD_out, void *PD_in, comp_d **coeff, int freeCoeff, int inType);
void cp_point_mp_int(void *Out, void *In, char **str, int freeStr, int initPoint, int inType);

void cp_mat_d_int(void *Out, void *In, comp_d **coeff, int inType);
void cp_mat_mp_int(void *Out, void *In, char **str, int freeStr, int inType);
void cp_mat_rat_int(void *Out, void *In, char **str, int rows, int cols, int freeStr, int inType);

void cp_vec_rat_char(void *Out, void *In, int *length, int size, int freeStr, int inType);
void cp_mat_rat_char(void *Out, void *In, int *length, int rows, int cols, int freeStr, int inType);

void cp_prog_t_int(void *Prog_out, void *Prog_in, int **instructions, int **gp_sizes, char **progStr, int freeStr, int inType);

void cp_endgame_data_t(endgame_data_t *EG_out, endgame_data_t *EG_in, int initEG_out);
void cp_trackBack_samples_t(trackBack_samples_t *TB_out, trackBack_samples_t *TB_in, int initTB_out);

#ifdef _HAVE_MPI // parallel copy functions

void cp_tracker_config_t_relevant(void *T_out, void *T_in, int inType);
void cp_basic_eval_data_d_int(void *BED_out, void *BED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int freeStr, int **orig_deg, int **new_deg, int **P, int **W, int **startDeg, comp_d **st_coeff, comp_d **sq_coeff, comp_d **patch_coeff, int **ppd_type, int **ppd_size, int MPType, int inType);
void cp_basic_eval_data_mp_int(void *BED_out, void *BED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int setupProg, int **orig_deg, int **new_deg, int **P, int **W, char **sqStr, char **patchStr, int **startDeg, char **startStr, int **ppd_type, int **ppd_size, int freeStr, int inType);
void cp_preproc_data_int(void *PPD_out, void *PPD_in, int **type, int **size, int inType);
void cp_endgame_data_t_int(void *EG_out, void *EG_in, char **egStr, comp_d **coeff, int freeStr, int freeCoeff, int inType);

void cp_eqData_int(void *EqD_out, void *EqD_in, int MPType, char **eqdStr, int freeStr, comp_d **coeff_d, int **degrees, int **instCount, int inType);
void cp_witnessData_int(void *Out, void *In, int stage, int MPType, char **wdStr, int freeStr, comp_d **coeff_d, int inType);
void cp_stageData_int(void *Out, void *In, int stage, int MPType, char **sdStr, int freeStr, comp_d **coeff_d, int inType);
void cp_trackBack_samples_t_int(void *TB_out, void *TB_in, char **egStr, comp_d **egCoeff, double **tbDouble, comp_d **tbCoeff, char **tbStr, int **tbInt, int freeStr, int freeCoeff, int inType);

#endif

void pos_dim_main(int trackType, int genType, int MPType, unsigned int currentSeed, char *startName, int my_id, int num_processes, int headnode);

// clear_functions.c prototypes
void tracker_config_clear(tracker_config_t *T);

void preproc_data_clear(preproc_data *PPD);

void basic_eval_clear_d(basic_eval_data_d *ED, int clearRegen, int MPType);
void basic_eval_clear_mp(basic_eval_data_mp *ED, int clearRegen, int clrProg);

void square_system_eval_data_clear_d(square_system_eval_data_d *SSED, int MPType);
void square_system_eval_data_clear_mp(square_system_eval_data_mp *SSED, int clrProg);

void start_system_eval_data_clear_d(start_system_eval_data_d *SSED);
void start_system_eval_data_clear_mp(start_system_eval_data_mp *SSED);

void patch_eval_data_clear_d(patch_eval_data_d *PED);
void patch_eval_data_clear_mp(patch_eval_data_mp *PED);

// Prototypes for sharpen.c
void sharpen_process_main(int MPType, int trackType, unsigned int currentSeed, int my_id, int num_processes, int headnode);
void create_raw_data_from_endPoints(post_process_t *endPoints, int num_endPoints, int num_variables, FILE *rawOUT);

int sharpen_zero_main(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void sharpen_endpoint_endgame(endgame_data_t *endPoint, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

// misc.c
void printPatchCoeff(FILE *OUT, int MPType, basic_eval_data_d *BED_d, basic_eval_data_mp *BED_mp);
void readInPatch(FILE *IN, int MPType, int old_MPType, patch_eval_data_d *PED_d, patch_eval_data_mp *PED_mp);
void move_to_patch_mat_d(point_d outPt, point_d inPt, mat_d patch_d, preproc_data *PPD);
void move_to_patch_d(point_d outPt, point_d inPt, basic_eval_data_d *BED);
void move_to_patch_mat_mp(point_mp outPt, point_mp inPt, mat_mp patch_mp, preproc_data *PPD);
void move_to_patch_mp(point_mp outPt, point_mp inPt, basic_eval_data_mp *BED);
void move_to_patch_vec_d(point_d outPt, point_d inPt, vec_d patch);
void move_to_patch_vec_mp(point_mp outPt, point_mp inPt, vec_mp patch);

// eqbyeq_setup.c
int eqbyeq_setup_d(FILE **OUT, char *outName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, char *pointsIN, char *pointsOUT, char *depthFile, double intrinsicCutoffMultiplier);
int eqbyeq_setup_mp(FILE **OUT, char *outName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, char *pointsIN, char *pointsOUT, char *depthFile, double intrinsicCutoffMultiplier);

void setup_omp_eqbyeq_d(basic_eval_data_d *BED, eqData_t *EqD_in, int MPType);
void setup_omp_eqbyeq_mp(basic_eval_data_mp *BED, eqData_t *EqD_in);

void setupEqbyEqFirstStage_d(eqData_t *EqD, int MPType);
void setupEqbyEqFirstStage_mp(eqData_t *EqD);

void clearEqbyEqFirstWitnessData_d(eqData_t *EqD, int MPType);
void clearEqbyEqFirstWitnessData_mp(eqData_t *EqD);

void setupEqbyEqNextStage_d(basic_eval_data_d *ED, int stage, int MPType, int max_prec);
void setupEqbyEqNextStage_mp(basic_eval_data_mp *ED, int stage);

void clearEqbyEqStageData_d(eqData_t *EqD, int stage, int MPType);
void clearEqbyEqStageData_mp(eqData_t *EqD, int stage);

void clearEqbyEqWitnessData_d(eqData_t *EqD, int subsystem, int MPType);
void clearEqbyEqWitnessData_mp(eqData_t *EqD, int subsystem);

void setup_omp_eqbyeq_first_stage_d(int max_threads, basic_eval_data_d ED_copy[], eqData_t *EqD, int MPType);
void setup_omp_eqbyeq_first_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], eqData_t *EqD);

void setup_omp_eqbyeq_next_stage_d(int max_threads, basic_eval_data_d ED_copy[], eqData_t *EqD, int stage, int MPType);
void setup_omp_eqbyeq_next_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], eqData_t *EqD, int stage);

void clear_omp_eqbyeq_stage_d(int max_threads, basic_eval_data_d ED_copy[], int stage, int MPType);
void clear_omp_eqbyeq_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], int stage);

void clear_omp_eqbyeq_witness_d(int max_threads, basic_eval_data_d ED_copy[], int subsystem, int MPType);
void clear_omp_eqbyeq_witness_mp(int max_threads, basic_eval_data_mp ED_copy[], int subsystem);

void clear_omp_eqData_d(basic_eval_data_d *BED, int MPType);
void clear_omp_eqData_mp(basic_eval_data_mp *BED);

// eqbyeq_track.c
void eqbyeq_track_d(FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, trackingStats *trackCount);
void eqbyeq_track_mp(FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, tracker_config_t *T, double midpoint_tol, double target_tol, basic_eval_data_mp *ED, trackingStats *trackCount);

int printWitnessFooter_d(basic_eval_data_d *ED, int subsystem_num, int path_num, point_data_d *endPoint, point_d orig_vars, point_d orig_last_d, point_d dehom, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

int printWitnessFooter_mp(basic_eval_data_mp *ED, int subsystem_num, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp orig_last_mp, point_mp dehomP, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

int printStageFooter_d(basic_eval_data_d *ED, int stage_num, int path_num, point_data_d *endPoint, point_d orig_vars, point_d orig_last_d, point_d dehom_d, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

int printStageFooter_mp(basic_eval_data_mp *ED, int stage_num, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp orig_last_mp, point_mp dehom_mp, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, FILE *NONSOLN, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

void setup_regen_from_zero_dim_seq(int max, regen_t *regen, char *startName, int startLevel, double intrinsicCutoffMultiplier, char *depthName, char *degreeName, tracker_config_t *T, basic_eval_data_d *BED_d, basic_eval_data_mp *BED_mp, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow);
void clearRegenRandom_zero_dim(int max, regen_t *regen, int MPType);

void eqbyeqWitnessSortEndpoint_d(int *rankDef, int *finite, int higherDim, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int subsystem, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void eqbyeqWitnessSortEndpoint_mp(int *rankDef, int *finite, int higherDim, basic_eval_data_mp *ED, int subsystem, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_mp *endPt_mp, int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

void eqbyeqStageSortEndpoint_d(int *rankDef, int *finite, int higherDim, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int stage_num, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void eqbyeqStageSortEndpoint_mp(int *rankDef, int *finite, int higherDim, basic_eval_data_mp *ED, int stage_num, int path_num, double *condNum, tracker_config_t *T, FILE *OUT, point_data_mp *endPt_mp, int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int determineEqbyEqHigherDim(double tol, double ratio, point_data_d *Pt_d, point_data_mp *Pt_mp, int Pt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, basic_eval_data_d *ED, basic_eval_data_mp *ED_mp, int num, int isStage);

// New cascade

// Endpoint classification codes:
#define UNCLASSIFIED 0
#define SOLUTION_AND_NONSING 10
#define SOLUTION_AND_SING 15
#define NONSOLUTION_AND_NONSING 20
#define NONSOLUTION_AND_SING 25
// < 0 -> bad path (code = retVal returned from path)

typedef struct
{
  point_d endPt;              // end point for the path
  point_d last_approx;        // last approximation to the end point
  comp_d finalT;              // final T value
  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value
} endpoint_data_d;

typedef struct
{
  point_mp endPt;             // end point for the path
  point_mp last_approx;       // last approximation to the end point
  comp_mp finalT;             // final T value
  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value
} endpoint_data_mp;

typedef struct
{
  int curr_prec;              // current precision
  point_d endPt_d;            // end point for the path - double precision
  comp_d finalT_d;            // final T value - double precision
  point_mp endPt_mp;          // end point for the path - mulit precision
  comp_mp finalT_mp;          // final T value - multi precision

  int last_approx_prec;       // precision of the last approximation
  point_d last_approx_d;      // last approximation to the end point
  point_mp last_approx_mp;    // last approximation to the end point

  double cond_num;            // condition number
  int corank;                 // corank at the endpoint
  double smallest_nonzero_SV; // smallest singular value that is non-zero
  double largest_zero_SV;     // largest singular value that is zero
  int retVal;                 // return value
} endpoint_data_amp;

typedef struct
{
  int codim;                     // current codimension
  int num_paths;                 // number of paths to be tracked for this codim
  int num_superset;              // number of points in the witness superset
  int num_nonsing;               // number of points in witness superset that are non-singular
  int num_sing;                  // number of points in witness superset that are singular
  int num_nonsolns;              // number of nonsolutions
  int num_inf;                   // number of infinite endpoints
  int num_other;                 // number of other endpoints that can be removed
  int num_bad;                   // number of bad endpoints

  point_d *startPts_d;           // start point for each path - double precision
  point_mp *startPts_mp;         // start point for each path - multi precision

  endpoint_data_d *endPts_d;     // end point data for each path - double precision
  endpoint_data_mp *endPts_mp;   // end point data for each path - multi precision
  endpoint_data_amp *endPts_amp; // end point data for each path - AMP

  int *endPt_types;              // end point types
} cascadeCodim_t;

typedef struct
{
  prog_t *Prog;           // slp to evaluate the original function f
  preproc_data PPD;       // preprocessing data

  mat_d A_d;              // random matrix - randomize the original system to 'square system'
  mat_mp A_mp;
  mpq_t ***A_rat;
  int **W;                // exponents of homogenizing variable to maintain the 1-hom of randomized system

  vec_d H_d;              // random vector - used to maintain the 1-hom of randomized system - (H*x)^W will maintain 1-hom
  vec_mp H_mp;
  mpq_t **H_rat;

  comp_d homVarConst_d;   // the constant term in the homogeneous variable expression
  comp_mp homVarConst_mp;
  mpq_t homVarConst_rat[2];

  int *orig_degrees;      // degrees of f_i
  int *new_degrees;       // degrees of f_P[i]
  int *P;                 // permutations of f_i's used to minimize paths

  int system_rank;        // rank of the original function f

  int orig_variables;     // number of variables used in 'Prog'
  int new_variables;      // number of variables after slicing
  mat_d C_d;              // matrix to convert new variables to orig variables
  mat_mp C_mp;
  mpq_t ***C_rat;

  comp_d gamma_d;         // the random gamma for the initial homotopy
  comp_mp gamma_mp;
  mpq_t gamma_rat[2];

  int curr_precision;     // current precision

  int num_funcs;          // number of functions
  int num_codim;          // number of codimensions
  int curr_codim_index;   // current codimension

  mat_d R_d;              // random matrix - randomize the linears into the 'square system'
  mat_mp R_mp;
  mpq_t ***R_rat;

  vec_d T_d;              // vector to t-values that control the cascade homotopy by turning off linear slices
  vec_mp T_mp;

  mat_d B_d;              // random matrix - 'coefficients' for linear slices
  mat_mp B_mp;
  mpq_t ***B_rat;

  vec_d p_d;              // random vector - 'patch coefficients'
  vec_mp p_mp;
  mpq_t **p_rat;

  int *W_prime;           // exponents of homogenizing variable to maintin the 1-hom of cascade system

  cascadeCodim_t *codim;  // store the codimension data
} cascade_t;

// cascade_setup.c
void cascade_setup(FILE **OUT, char *outName, FILE **RAWOUT, char *rawName, FILE **MIDOUT, char *midName, FILE **FAIL, char *failName, tracker_config_t *T, cascade_t *CD, char *preprocFile, char *degreeFile, int maxCodim, int specificCodim);

void cascade_clear(cascade_t *CD, int MPType);
void clearCascadeCodim(cascade_t *CD, int codim_index, int MPType);

void allocateCascadeCodim(cascade_t *CD, int codim_index, int codim, int num_paths, int MPType);

void cascade_clear_start_points(cascade_t *CD, int codim_index, int MPType);

void setup_omp_cascade_t(cascade_t *CD, cascade_t *CD_in, int MPType);
void clear_omp_cascade_t(cascade_t *CD, int MPType);

int change_cascade_prec(void const *ED, int prec);

// cascade_track.c
void setup_cascade_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, char *outName, FILE ***RAWOUT_copy, FILE *RAWOUT, char *rawName, FILE ***MIDOUT_copy, FILE *MIDOUT, char *midName, FILE ***FAIL_copy, FILE *FAIL, char *failName, tracker_config_t **T_copy, tracker_config_t *T, cascade_t **CD_copy, cascade_t *CD);
void clear_cascade_omp(int max_threads, trackingStats **trackCount_copy, trackingStats *trackCount, char *outName, char *rawName, char *midName, char *failName, tracker_config_t **T_copy, cascade_t **CD_copy);

void cascadePrepareNextCodim(cascade_t *CD, int curr_codim_index, int MPType);

void cascadeFindOrigVarsDehom_d(point_d orig_vars_d, point_d dehom_d, point_d P_d, cascade_t *CD);
void cascadeFindOrigVarsDehom_mp(point_mp orig_vars_mp, point_mp dehom_mp, point_mp P_mp, cascade_t *CD);

int printCascadeFooter_d(cascade_t *CD, int codim_index, int path_num, point_data_d *endPoint, point_d orig_vars, point_d dehomP, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

int printCascadeFooter_mp(cascade_t *CD, int codim_index, int path_num, point_data_mp *endPoint, point_mp orig_vars, point_mp dehomP, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, FILE *RAWOUT, FILE *FAIL, int retVal_in, tracker_config_t *T, trackingStats *trackCount);

void findDiff_endpoint_data(mpf_t norm_diff, endpoint_data_amp endPt1, endpoint_data_amp endPt2);

void cascade_classifyCodim(cascade_t *CD, int codim_index, double final_tol, int MPType);

int determineCascadeFinite_d(double max_norm, point_d point, cascade_t *CD);
int determineCascadeFinite_mp(double max_norm, point_mp point, int prec, cascade_t *CD);

int determineCascadeSoln_d(double tol, double ratio, point_d point, point_d last_approx, comp_d time, cascade_t *CD, int codim_index);
int determineCascadeSoln_mp(double tol, double ratio, point_mp point, point_mp last_approx, comp_mp time, int prec, cascade_t *CD, int codim_index);

void cascadeSortEndpoint(int *rankDef, int *finite, int *soln, double condNum, cascade_t *CD, int codim_index, int path_num, tracker_config_t *T, FILE *OUT, point_data_d *endPt_d, point_data_mp *endPt_mp, int endPt_prec, point_d last_approx_d, point_mp last_approx_mp, int last_approx_prec, int (*change_prec)(void const *, int));

// dehomogenize the cascade point
int cascade_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);

// cascade_eval.c
int initial_codim_cascade_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int standard_codim_cascade_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);

int initial_codim_cascade_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int standard_codim_cascade_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

// cascade_output.c
void cascadeOutputChart(cascade_t *CD, FILE *fp, int MPType);

// copy_functions.c
void cp_cascade_int(void *Out, void *In, int MPType, char **cascadeStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType);
void cp_cascadeCodim_int(void *Out, void *In, int inType);

void cp_endpoint_data_d_int(void *Out, void *In, comp_d **coeff, int freeCoeff, int initPoint, int inType);
void cp_endpoint_data_mp_int(void *Out, void *In, char **epStr, int freeStr, int initPoint, int inType);
void cp_endpoint_data_amp_int(void *Out, void *In, comp_d **coeff, char **epStr, int freeData, int initPoint, int inType);

// parallel_send_functions.c
void bcast_cascade_t(cascade_t *CD, int MPType, int my_id, int headnode);
void bcast_cascadeCodim_t(cascadeCodim_t *CD, int my_id, int headnode);

// cascade_parallel.c
void cascade_par_track(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *FAIL, char *midFile, int pathMod, double midpoint_tol, tracker_config_t *T, cascade_t *CD, int my_id, int num_processes, int headnode);

// regen_pos_dim_track.c
void print_endpoint_data_d(FILE *fp, endpoint_data_d *endPt);
void setup_endpoint_data_d(endpoint_data_d *endPt, FILE *fp);
void print_endpoint_data_mp(FILE *fp, endpoint_data_mp *endPt, int prec);
void setup_endpoint_data_mp(endpoint_data_mp *endPt, FILE *fp);
void print_endpoint_data_amp(FILE *fp, endpoint_data_amp *endPt);
void setup_endpoint_data_amp(endpoint_data_amp *endPt, FILE *fp);

// initialize endpoint_data_*
#define init_endpoint_data_d(_r) { init_point_d((_r)->endPt, 0); init_point_d((_r)->last_approx, 0); init_d((_r)->finalT); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }
#define init_endpoint_data_mp(_r) { init_point_mp((_r)->endPt, 0); init_point_mp((_r)->last_approx, 0); init_mp((_r)->finalT); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }
#define init_endpoint_data_mp2(_r, _prec) { init_point_mp2((_r)->endPt, 0, _prec); init_point_mp2((_r)->last_approx, 0, _prec); init_mp2((_r)->finalT, _prec); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }
#define init_endpoint_data_amp(_r, _prec, _last_prec) { int _p = MAX(_prec, 64), _l_p = MAX(_last_prec, 64); \
  (_r)->curr_prec = _prec; (_r)->last_approx_prec = _last_prec; \
  init_point_d((_r)->endPt_d, 0); init_d((_r)->finalT_d); init_point_mp2((_r)->endPt_mp, 0, _p); init_mp2((_r)->finalT_mp, _p); \
  init_point_d((_r)->last_approx_d, 0); init_point_mp2((_r)->last_approx_mp, 0, _l_p); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }

// copy endpoint_data_*
#define endpoint_data_cp_d(_r, _a) { point_cp_d((_r)->endPt, (_a)->endPt); point_cp_d((_r)->last_approx, (_a)->last_approx); \
  set_d((_r)->finalT, (_a)->finalT); (_r)->cond_num = (_a)->cond_num; (_r)->corank = (_a)->corank; (_r)->smallest_nonzero_SV = (_a)->smallest_nonzero_SV; \
  (_r)->largest_zero_SV = (_a)->largest_zero_SV; (_r)->retVal = (_a)->retVal; }
#define endpoint_data_cp_mp(_r, _a) { point_cp_mp((_r)->endPt, (_a)->endPt); point_cp_mp((_r)->last_approx, (_a)->last_approx); \
  set_mp((_r)->finalT, (_a)->finalT); (_r)->cond_num = (_a)->cond_num; (_r)->corank = (_a)->corank; (_r)->smallest_nonzero_SV = (_a)->smallest_nonzero_SV; \
  (_r)->largest_zero_SV = (_a)->largest_zero_SV; (_r)->retVal = (_a)->retVal; }
#define endpoint_data_cp_amp(_r, _a) { int _p = (_r)->curr_prec = (_a)->curr_prec; \
  if (_p < 64) { point_cp_d((_r)->endPt_d, (_a)->endPt_d); set_d((_r)->finalT_d, (_a)->finalT_d); } \
  else { setprec_point_mp((_r)->endPt_mp, _p); setprec_mp((_r)->finalT_mp, _p); point_cp_mp((_r)->endPt_mp, (_a)->endPt_mp); set_mp((_r)->finalT_mp, (_a)->finalT_mp); }\
  _p = (_r)->last_approx_prec = (_a)->last_approx_prec; \
  if (_p < 64) { point_cp_d((_r)->last_approx_d, (_a)->last_approx_d); } \
  else { setprec_point_mp((_r)->last_approx_mp, _p); point_cp_mp((_r)->last_approx_mp, (_a)->last_approx_mp); } \
  (_r)->cond_num = (_a)->cond_num; (_r)->corank = (_a)->corank; (_r)->smallest_nonzero_SV = (_a)->smallest_nonzero_SV; \
  (_r)->largest_zero_SV = (_a)->largest_zero_SV; (_r)->retVal = (_a)->retVal; }

// clear endpoint_data_* 
#define clear_endpoint_data_d(_r) { clear_point_d((_r)->endPt); clear_point_d((_r)->last_approx); clear_d((_r)->finalT); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; } 
#define clear_endpoint_data_mp(_r) { clear_point_mp((_r)->endPt); clear_point_mp((_r)->last_approx); clear_mp((_r)->finalT); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }
#define clear_endpoint_data_amp(_r) { clear_point_d((_r)->endPt_d); clear_d((_r)->finalT_d); clear_point_mp((_r)->endPt_mp); clear_mp((_r)->finalT_mp); \
  clear_point_d((_r)->last_approx_d); clear_point_mp((_r)->last_approx_mp); \
  (_r)->cond_num = (_r)->corank = (_r)->smallest_nonzero_SV = (_r)->largest_zero_SV = (_r)->retVal = 0; }

#endif

