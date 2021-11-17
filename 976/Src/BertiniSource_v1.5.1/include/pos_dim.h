// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#ifndef _POS_DIM_H
#define _POS_DIM_H

#include "bertini.h"
#include "cascade.h"
#include "dimbydim.h"

#define NON_SINGULAR 10
#define SINGULAR 15

// data structure for witness (super) set
typedef struct
{
  int codim;                         // current codimension
  int num_set;                       // number of points in the witness set
  int num_nonsing;                   // number of points in witness set that are non-singular
  int num_sing;                      // number of points in witness set that are singular

  mat_d A_d;                         // random matrix - randomize the system to the proper codim
  mat_mp A_mp;
  mpq_t ***A_rat;
  int A_rows; 
  int A_cols;

  int **W;                           // exponents of homogenizing variable to maintain the 1-hom of randomized system

  vec_d H_d;                         // random vector - used to maintain the 1-hom of randomized system - (H*x + homVarConst)^W will maintain 1-hom
  vec_mp H_mp;
  mpq_t **H_rat;

  comp_d homVarConst_d;              // the constant term in the homogeneous variable expression
  comp_mp homVarConst_mp;
  mpq_t homVarConst_rat[2];

  mat_d B_d;                         // random matrix - coefficients for linear slices
  mat_mp B_mp;
  mpq_t ***B_rat;

  vec_d p_d;                         // random vector - patch coefficients
  vec_mp p_mp;
  mpq_t **p_rat;

  endpoint_data_d *witnessPts_d;     // witness point data - double precision
  endpoint_data_mp *witnessPts_mp;   // witness point data - multi precision
  endpoint_data_amp *witnessPts_amp; // witness point data - AMP

  int *witnessPt_types;              // type of witness points (either singular or non-singular)

  int num_components;                // number of irreducible components for this pure-dimensional witness set
  int *component_nums;               // describe which points go with which components
  int *multiplicities;               // describe the multiplicities for the points - > 0 means the multiplicity, otherwise the negative of duplicate path numb  
  int *deflations_needed;            // number of deflations needed for this point to become non-singular

} witnessCodim_t;

typedef struct
{
  prog_t *Prog;          // slp to evaluate the original function f
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

  comp_d gamma_d;        // the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t gamma_rat[2];

  int curr_precision;    // current precision

  int num_funcs;         // number of functions
  int num_codim;         // number of codimensions
  int curr_codim_index;  // current codimension

  // used for junk removal
  int targetSliceInit;   // determine if the target slices have been initialized or not
  mat_d targetSliceMat_d;// random matrix - 'coefficients' (ext) or 'null space' (int)
  mat_mp targetSliceMat_mp;
  mpq_t ***targetSliceMat_rat;

  vec_d targetSliceVec_d;// random vector - 'patch coefficients' (ext) or 'solution vector' (int)
  vec_mp targetSliceVec_mp;
  mpq_t **targetSliceVec_rat;

  witnessCodim_t *codim; // holds the witness data
} witness_t;

typedef struct
{
  prog_t *Prog; // slp to evaluate f

  int curr_codim; // codimension that this represents
  int orig_variables; // original number of variables (extrinsically)
  int curr_precision; // current precision

  comp_d gamma_d;// the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t gamma_rat[2];

  mat_d A_d;    // random matrix - randomize the system to the proper codim
  mat_mp A_mp;
  mpq_t ***A_rat;

  mat_d B_d;    // random matrix - coefficients for linear slices
  mat_mp B_mp;
  mpq_t ***B_rat;

  vec_d p_d;    // random vector - patch coefficients
  vec_mp p_mp;
  mpq_t **p_rat;

  vec_d startSliceVec_d;// vector of where to start at
  vec_mp startSliceVec_mp;
  mpq_t **startSliceVec_rat;
  int startSliceVec_init; // whether startSliceVec has been initialized

  vec_d targetSliceVec_d;// vector of where to finish at
  vec_mp targetSliceVec_mp;
  mpq_t **targetSliceVec_rat;
  int targetSliceVec_init; // whether targetSliceVec has been initialized

  mpq_t ***K_rat;// kernel of the linear space defined by B & p
  int K_rows;
  int K_cols;

} membership_slice_moving_t;

typedef struct
{
  prog_t *Prog;     // slp to evaluate f
  preproc_data PPD; // preprocessing data

  int curr_codim;     // codimension that this represents
  int orig_variables; // original number of variables (extrinsically)
  int num_funcs;      // number of functions
  int curr_precision; // current precision
  int *P;             // permutations of f_i's used to minimize paths

  comp_d gamma_d; // the random gamma for the homotopy
  comp_mp gamma_mp;
  mpq_t gamma_rat[2];

  mat_d A_d;  // randomize the system to the proper codim
  mat_mp A_mp;
  mpq_t ***A_rat;
  int A_rows;
  int A_cols;

  int **W;    // exponents of homogenizing variable to maintain the 1-hom of randomized system

  vec_d H_d;  // used to maintain the 1-hom of randomized system - (H*x + homVarConst)^W will maintain 1-hom
  vec_mp H_mp;
  mpq_t **H_rat;

  comp_d homVarConst_d; // the constant term in the homogeneous variable expression
  comp_mp homVarConst_mp;
  mpq_t homVarConst_rat[2];

  mat_d B_start_d;    // coefficients for the starting linear slices
  mat_mp B_start_mp;
  mpq_t ***B_start_rat;

  mat_d B_end_d;    // coefficients for the ending linear slices
  mat_mp B_end_mp;
  mpq_t ***B_end_rat;

  vec_d patch_start_d;    // patch coefficients
  vec_mp patch_start_mp;
  mpq_t **patch_start_rat;

  vec_d patch_end_d;    // patch coefficients
  vec_mp patch_end_mp;
  mpq_t **patch_end_rat;

} general_slice_moving_t;

// sharpen.c
int pos_dim_sharpening_menu(tracker_config_t *T);

// pos_dim.c
void sampleComponent(unsigned int currentSeed, int MPType, int useSharpeningMenu, int my_id, int num_processes, int headnode);
void sampleComponentMenu(witness_t *W, tracker_config_t *T, int pathMod);
void printSamplePoints(point_d *samples_d, point_mp *samples_mp, int *samples_prec, int num_samples, FILE *fp, int useMatlab);
int generateSamplePoints(point_d *samples_d, point_mp *samples_mp, int *samples_prec, witness_t *W, tracker_config_t *T, int codim_index, int component_number, int num_samples, int pathMod);
int sampleTrack(point_d endPt_d, point_mp endPt_mp, int *endPt_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);
void setupWitnessDataFromFile(char *witnessName, char *newWitnessName, char *preprocFile, char *degreeFile, witness_t *W, tracker_config_t *T, int needToSetupProg);
void witness_clear(witness_t *W, int MPType);
void witness_clear_codim(witnessCodim_t *WC, int MPType);
int change_witness_prec(void const *ED, int prec);
void dimbydim_copyWitness_clear(witness_t *witnessSuperset, codim_t *CD, int MPType, int max_prec);
void dimbydim_copyWitness_clear_codim(witnessCodim_t *witCodim, codimData_t *CD, int MPType, int curr_prec, int max_prec);
int cascade_copyWitness_clear(witness_t *witnessSuperset, cascade_t *CD, int MPType, int max_prec, int specificCodim);
void cascade_copyWitness_clear_codim(witnessCodim_t *witCodim, cascade_t *CD, int codim_index, int MPType, int curr_prec, int max_prec);
void printCodimWitnessStructures(FILE *OUT, witnessCodim_t *WC, int MPType);
void printCodimWitnessSet(FILE *OUT, witnessCodim_t *WC, int MPType);
void setupDegrees_orig_new_perm(int **orig_degrees, int **new_degrees, int **perm, int num_funcs, int num_var_gps, char *degreeFile);

// numericalIrredDecomp.c
void numericalIrredDecomp(unsigned int currentSeed, int MPType, int genType, int my_id, int num_processes, int headnode);
void remove_junk_points(witness_t *W, int codim_index, int MPType, int *isJunk);
void witnessSetOutputChart(witness_t *W, FILE *fp, int MPType, int reducedOnly);
void multiplicity_witness(witness_t *W, int codim_index, int MPType, double tol);
void sort_endpoint_data(int *mult, endpoint_data_d *Pts_d, endpoint_data_mp *Pts_mp, endpoint_data_amp *Pts_amp, int numPoints, int MPType, double tol);
void witnessSupersetOutput(witness_t *W, int MPType);
void setupWitnessTotallyExtrinisic(witness_t *W, int MPType, int max_prec);
void displayDeflationSummary(int **fullRankProgInfo, witness_t *W);

void numIrredDecompChart(witness_t *W, FILE *fp, int MPType, int reducedOnly);
void numIrredDecompOutput(witness_t *W, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom);

void witnessFindDehom_d(point_d dehom_d, point_d P_d, witness_t *W, int codim_index);
void witnessFindDehom_mp(point_mp dehom_mp, point_mp P_mp, witness_t *W, int codim_index, int P_prec);
void deflate_for_junkRemoval(prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int component_number, tracker_config_t *T, FILE *OUT);
void clear_sliceMover_fullRankProgs(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, int MPType);
void clear_sliceMover_fullRankProgs_codim(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, int num_set, int MPType);
void clear_fullRankProg_endPt(prog_t **fullRankProg, int fullRankProgInfo, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, int MPType);
void sharpen_deflated(witness_t *W, int codim_index, prog_t **fullRankProgs, int *fullRankProgInfo, membership_slice_moving_t *sliceMover, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT);
void numIrredDecompWitnessData(char *witnessName, witness_t *W, int MPType);

// membership.c
void membershipMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode);
int pureDimMembershipTest(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, double midpoint_tol, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode);
int printMembershipFooter_d(point_data_d *endPoint, double cond_num, double func_residual, double newton_error, double t_val_sample, double error_sample, FILE *OUT, int retVal_in, tracker_config_t *T);
int printMembershipFooter_mp(point_data_mp *endPoint, double cond_num, double first_increase, mpf_t func_residual, mpf_t newton_error, mpf_t t_val_sample, mpf_t error_sample, FILE *OUT, int retVal_in, tracker_config_t *T);

void junkRemoval(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);
void junkRemoval_par(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);
void junkRemoval_seq(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);

int junkRemoval_mem(witness_t *W, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode);
int junkRemoval_ldt(witness_t *W, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int my_id, int num_processes, int headnode);

// membership_junk.c
int junkRemoval_membershipTest(witness_t *W, int pure_codim_index, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t **fullRankProgs, int *fullRankProgInfo, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode);
void basic_setup_slice_moving(membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int MPType, int max_prec);
int slice_mover_change_prec(void const *ED, int prec);
void setup_slice_moving_slice(membership_slice_moving_t *sliceMover, point_d endPt_d, point_mp endPt_mp, int endPt_prec, int MPType, int max_prec);
int membership_slice_moving_track(int *isMember, point_d testPt_d, point_mp testPt_mp, int testPt_prec, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, double tol, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);
void setup_random_slice_moving(membership_slice_moving_t *sliceMover, int MPType, int max_prec);
void clear_slice_mover(membership_slice_moving_t *SM, int MPType);


// membership_eval.c
// the membership slice evaluators do things intrinsically from orig_variables to new_variables
// while the true slice evaluators do things extrinisically in orig_variables
int standard_witness_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int membership_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int decomposition_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int standard_witness_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int membership_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int decomposition_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);

int slice_mover_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int slice_mover_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
int slice_moving_track(endgame_data_t *endPt, membership_slice_moving_t *sliceMover, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, int checkSharpen, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);
void final_setup_slice_moving(membership_slice_moving_t *sliceMover, prog_t *fullRankProg, int MPType, int max_prec, int setupGamma);
void initialize_slice_moving_sliceVec(membership_slice_moving_t *sliceMover, int size, int MPType);

// cascade_track.c
int cascade_main(witness_t *witnessSuperset, int maxCodim, int specificCodim, tracker_config_t *T, int pathMod, double midpoint_tol, double intrinsicCutoffMultiplier, int my_id, int num_processes, int headnode);

// dimbydim_track.c
void dimbydim_main(witness_t *witnessSuperset, int maxCodim, int specificCodim, tracker_config_t *T, int pathMod, double midpoint_tol, double intrinsicCutoffMultiplier, int my_id, int num_processes, int headnode);

// pureDecomp.c
void pureDecomp(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode);
void pureDecomp_par(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode);
void pureDecomp_seq(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgsInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, int my_id, int num_processes, int headnode);

void calculateTrace(comp_d trace_d, comp_mp trace_mp, int *trace_prec, membership_slice_moving_t *sliceMover, prog_t *fullRankProg, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, int codim, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, vec_d proj_d, vec_mp proj_mp, mpq_t **proj_rat, vec_d v_d, vec_mp v_mp, mpq_t **v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], int setupGamma);

void calculateGroupTrace_d(comp_d gp_trace, int num_points, int *gp_nums, int curr_gp, comp_d *traces);
void calculateGroupTrace_mp(comp_mp gp_trace, int num_points, int *gp_nums, int curr_gp, comp_mp *traces);
void checkCompleteComponents(witness_t *W, int codim_index, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec);

void traceDecomposition(witness_t *W, int codim_index, double tol, int *num_to_classify, int **isClassified, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec);
int updateTempGroupsFromLoops(int num_paths, int num_to_classify, int *isClassified, int *loop_results, int *num_temp_gps, int **temp_gp_nums, comp_d **temp_gp_trace_d, comp_mp **temp_gp_trace_mp, int trace_prec);

int basicMonodromyLoop(endgame_data_t *monodromyPt, membership_slice_moving_t *sliceMover, prog_t *fullRankProgs, point_d Pt_d, point_mp Pt_mp, int Pt_prec, witness_t *W, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat);
int monodromyLoop(membership_slice_moving_t *sliceMover, prog_t *fullRankProgs, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, witness_t *W, int codim_index, int pathNum, tracker_config_t *T, FILE *OUT, FILE *MIDOUT, double tol, int *isClassified, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t *gamma_out_rat, comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t *gamma_in_rat, vec_d vec_out_d, vec_mp vec_out_mp, mpq_t **vec_out_rat, vec_d vec_in_d, vec_mp vec_in_mp, mpq_t **vec_in_rat, double *norms);

// parallel related items
void bcast_membership_slice_moving_t(membership_slice_moving_t *slice, int sendProg, int MPType, int my_id, int headnode);
void bcast_witnessCodim_t(witnessCodim_t *witCodim, int MPType, int curr_prec, int my_id, int headnode);
void bcast_witness_t(witness_t *W, int MPType, int my_id, int headnode);
void bcast_witness_structures(prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, int MPType, witness_t *W, int codim_index, int my_id, int headnode);
void bcast_trace_structures(vec_d proj_d, vec_mp proj_mp, mpq_t ***proj_rat, vec_d v_d, vec_mp v_mp, mpq_t ***v_rat, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], comp_d gamma_d, comp_mp gamma_mp, mpq_t gamma_rat[2], int MPType, int curr_prec, int my_id, int headnode);
void bcast_monodromy_structures(int *continueMonodromy, vec_d v_out_d, vec_mp v_out_mp, mpq_t ***v_out_rat, vec_d v_in_d, vec_mp v_in_mp, mpq_t ***v_in_rat, comp_d gamma_out_d, comp_mp gamma_out_mp, mpq_t gamma_out_rat[2], comp_d gamma_in_d, comp_mp gamma_in_mp, mpq_t gamma_in_rat[2], int MPType, int curr_prec, int my_id, int headnode);

void cp_slice_moving_int(void *Out, char **sliceStr, char **progStr, int **prog_inst, int **prog_gp_sizes, comp_d **coeff_d, void *In, int sendProg, int MPType, int freeStr, int inType);
void worker_witness_superset_decomposition(int my_id, int num_processes, int headnode);
void worker_junkRemoval(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int specificCodim, int topDimension, int my_id, int num_processes, int headnode);

// witnessGeneration.c
void witnessGeneration(unsigned int currentSeed, int MPType, char *startName, int my_id, int num_processes, int headnode);

// isosingular.c
int isosingularDimTest(membership_slice_moving_t *sliceMover, point_d *endPts_d, point_mp *endPts_mp, int *endPt_prec, int nullity, int numPoints, point_d *Pts_d, point_mp *Pts_mp, int prec, tracker_config_t *T, witness_t *W);

// witnessProjection.c
void witnessProjectionMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode);

// printWitness.c
void printWitnessMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode);
void printLinearSystem(witness_t *W, int codim, int component_num, char *pointsName, char *linearName, int MPType, int max_prec);

// generalSliceMoving.c
int general_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
int general_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED);
void initialize_setup_general_slice(general_slice_moving_t *M, witness_t *W_start, int start_codim_index, witness_t *W_end, int end_codim_index, int MPType);
int change_general_slice_prec(void const *ED, int prec);
void clear_general_slice(general_slice_moving_t *M, int MPType);
int general_slice_moving_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);
int general_slice_moving_track(endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, endgame_data_t *endGame, general_slice_moving_t *M, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, int checkSharpen, tracker_config_t *T, FILE *OUT, FILE *MIDOUT);
int projection_sortEndpoint(point_data_d *PD_d, point_data_mp *PD_mp, int Pt_prec, int retVal_in, preproc_data *PPD, trackingStats *trackCount, FILE *FAIL, int pathNum, tracker_config_t *T);

// diag_main.c
void diagIntersectionMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode);
void runDiagIntersection(witness_t *witnessA, int codim_indexA, int componentA, witness_t *witnessB, int codim_indexB, int componentB, tracker_config_t *T, int pathMod);

#endif

