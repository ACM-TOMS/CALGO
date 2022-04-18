// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"
#include "regeneration.h"
#include "regen_pos_dim.h"
#include "dimbydim.h"
#include "pos_dim.h"

///////////////////// COPY FUNCTIONS /////////////////////////////////

void cp_tracker_config_t(tracker_config_t *T, const tracker_config_t *T_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of T to T_input                          *
\***************************************************************/
{
  T->numVars = T_input->numVars;
  T->numPathVars = T_input->numPathVars;
  T->numParams = T_input->numParams;
  T->numFuncs = T_input->numFuncs;

  T->maxStepSize = T_input->maxStepSize;
  T->minStepSizeBeforeEndGame = T_input->minStepSizeBeforeEndGame;
  T->minStepSizeDuringEndGame = T_input->minStepSizeDuringEndGame;
  T->minStepSize = T_input->minStepSize;
  T->currentStepSize = T_input->currentStepSize;
  T->first_step_of_path = T_input->first_step_of_path;
  T->minTrackT = T_input->minTrackT;

  T->basicNewtonTol = T_input->basicNewtonTol;
  T->endgameNewtonTol = T_input->endgameNewtonTol;
  T->final_tolerance = T_input->final_tolerance;

  T->cSecInc = T_input->cSecInc;
  T->maxNewtonIts = T_input->maxNewtonIts;
  T->MPType = T_input->MPType;
  T->Precision = T_input->Precision;
  T->outputLevel = T_input->outputLevel;
  T->screenOut = T_input->screenOut;

  T->targetT = T_input->targetT;
  T->endgameBoundary = T_input->endgameBoundary;
  T->endgameSwitch = T_input->endgameSwitch;

  T->goingToInfinity = T_input->goingToInfinity;
  T->maxNumSteps = T_input->maxNumSteps;
  T->endgameNumber = T_input->endgameNumber;

  T->latest_cond_num_exp = T_input->latest_cond_num_exp;
  T->steps_since_last_CN = T_input->steps_since_last_CN;

  T->power_series_sample_factor = T_input->power_series_sample_factor;
  T->cycle_num_max = T_input->cycle_num_max;
  T->num_PSEG_sample_points = T_input->num_PSEG_sample_points;

  T->latest_newton_residual_d = T_input->latest_newton_residual_d;
  if (T->MPType == 1)
  { // initialize to the current precision
    mpf_init2(T->latest_newton_residual_mp, T_input->Precision);
  }
  else if (T->MPType == 2)
  { // initialize to the maximum precision
    mpf_init2(T->latest_newton_residual_mp, T_input->AMP_max_prec);
  }

  T->t_val_at_latest_sample_point = T_input->t_val_at_latest_sample_point;
  T->error_at_latest_sample_point = T_input->error_at_latest_sample_point;
  T->final_tolerance = T_input->final_tolerance;

  T->real_threshold = T_input->real_threshold;
  T->endgameOnly = T_input->endgameOnly;

  T->AMP_bound_on_abs_vals_of_coeffs = T_input->AMP_bound_on_abs_vals_of_coeffs;
  T->AMP_bound_on_degree = T_input->AMP_bound_on_degree;
  T->AMP_eps = T_input->AMP_eps;
  T->AMP_Phi = T_input->AMP_Phi;
  T->AMP_Psi = T_input->AMP_Psi;

  T->AMP_safety_digits_1 = T_input->AMP_safety_digits_1;
  T->AMP_safety_digits_2 = T_input->AMP_safety_digits_2;
  T->AMP_max_prec = T_input->AMP_max_prec;

  T->sing_val_zero_tol = T_input->sing_val_zero_tol;
  T->cond_num_threshold = T_input->cond_num_threshold;

  T->step_fail_factor = T_input->step_fail_factor;
  T->step_success_factor = T_input->step_success_factor;
 
  T->max_num_pts_for_trace = T_input->max_num_pts_for_trace;
  T->max_num_mon_linears = T_input->max_num_mon_linears;
  T->max_num_bad_loops_in_mon = T_input->max_num_bad_loops_in_mon;

  T->final_tol_multiplier = T_input->final_tol_multiplier;
  T->final_tol_times_mult = T_input->final_tol_times_mult;
 
  T->sharpenDigits = T_input->sharpenDigits;
  T->sharpenOnly = T_input->sharpenOnly;

  T->regen_remove_inf = T_input->regen_remove_inf;
  T->regen_higher_dim_check = T_input->regen_higher_dim_check;
  T->sliceBasicNewtonTol = T_input->sliceBasicNewtonTol;
  T->sliceEndgameNewtonTol = T_input->sliceEndgameNewtonTol;
  T->sliceFinalTol = T_input->sliceFinalTol;

  T->minCycleTrackBack = T_input->minCycleTrackBack;
  T->junkRemovalTest = T_input->junkRemovalTest;
  T->maxDepthLDT = T_input->maxDepthLDT;
  T->odePredictor = T_input->odePredictor;

  T->securityLevel = T_input->securityLevel;
  T->securityMaxNorm = T_input->securityMaxNorm;

  T->cutoffCycleTime = T_input->cutoffCycleTime;
  T->cutoffRatioTime = T_input->cutoffRatioTime;
  T->finiteThreshold = T_input->finiteThreshold;
  T->funcResTol = T_input->funcResTol;
  T->ratioTol = T_input->ratioTol;
  T->maxStepsBeforeNewton = T_input->maxStepsBeforeNewton;

  return;
}

void cp_basic_eval_data_d(basic_eval_data_d *BED, basic_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of BED_(t)_input to BED                  *
\***************************************************************/
{
  cp_preproc_data(&BED->preProcData, &BED_d_input->preProcData);
  cp_square_system_d(&BED->squareSystem, &BED_d_input->squareSystem);
  cp_patch_d(&BED->patch, &BED_d_input->patch);
  cp_start_system_d(&BED->startSystem, &BED_d_input->startSystem);

  if (MPType == 2)
  { // need to also setup MP versions since using AMP
    BED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));

    cp_preproc_data(&BED->BED_mp->preProcData, &BED_mp_input->preProcData);
    // simply point to the SLP that was setup in BED
    cp_square_system_mp(&BED->BED_mp->squareSystem, &BED_mp_input->squareSystem, 0, BED->squareSystem.Prog);
    cp_patch_mp(&BED->BED_mp->patch, &BED_mp_input->patch);
    cp_start_system_mp(&BED->BED_mp->startSystem, &BED_mp_input->startSystem);
  }
  else
    BED->BED_mp = NULL;

  return;
}

void cp_preproc_data(preproc_data *PPD, preproc_data *PPD_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PPD_input to PPD                      *
\***************************************************************/
{
  int i, total_gp;

  PPD->num_funcs = PPD_input->num_funcs;
  PPD->num_hom_var_gp = PPD_input->num_hom_var_gp;
  PPD->num_var_gp = PPD_input->num_var_gp;

  total_gp = PPD->num_hom_var_gp + PPD->num_var_gp;

  PPD->size = (int *)bmalloc(total_gp * sizeof(int));
  PPD->type = (int *)bmalloc(total_gp * sizeof(int));

  for (i = 0; i < total_gp; i++)
  {
    PPD->size[i] = PPD_input->size[i];
    PPD->type[i] = PPD_input->type[i];
  }

  return;
}

void cp_square_system_d(square_system_eval_data_d *SSED, square_system_eval_data_d *SSED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_input to SSED                    *
\***************************************************************/
{
  int i, j;

  // copy SLP
  SSED->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(SSED->Prog, SSED_input->Prog);

  SSED->size_f = SSED_input->size_f;
  SSED->size_r = SSED_input->size_r;
  SSED->max_of_W = SSED_input->max_of_W;
  SSED->noChanges = SSED_input->noChanges;

  SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
  SSED->P = (int *)bmalloc(SSED->size_f * sizeof(int));
  for (i = 0; i < SSED->size_f; i++)
  {
    SSED->orig_degrees[i] = SSED_input->orig_degrees[i];
    SSED->P[i] = SSED_input->P[i];
  }

  SSED->new_degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
  SSED->W = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
  {
    SSED->new_degrees[i] = SSED_input->new_degrees[i];
    SSED->W[i] = (int *)bmalloc((SSED->size_f - SSED->size_r) * sizeof(int));
    for (j = 0; j < (SSED->size_f - SSED->size_r); j++)
      SSED->W[i][j] = SSED_input->W[i][j];
  }

  // copy the matrices
  init_mat_d(SSED->B, SSED_input->B->rows, SSED_input->B->cols);
  mat_cp_d(SSED->B, SSED_input->B);
  init_mat_d(SSED->B_perp, SSED_input->B_perp->rows, SSED_input->B_perp->cols);
  mat_cp_d(SSED->B_perp, SSED_input->B_perp);
  init_mat_d(SSED->A, SSED_input->A->rows, SSED_input->A->cols);
  mat_cp_d(SSED->A, SSED_input->A);

  return;
}

void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d *PED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PED_input to PED                      *
\***************************************************************/
{
  PED->num_patches = PED_input->num_patches;
  init_mat_d(PED->patchCoeff, PED_input->patchCoeff->rows, PED_input->patchCoeff->cols);
  mat_cp_d(PED->patchCoeff, PED_input->patchCoeff);

  return;
}

void cp_start_system_d(start_system_eval_data_d *SSED, start_system_eval_data_d *SSED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_input to SSED                    *
\***************************************************************/
{
  int i, j, total_deg = 0;

  SSED->startSystemType = SSED_input->startSystemType;
  SSED->size_r = SSED_input->size_r;
  SSED->max_degree = SSED_input->max_degree;
  SSED->coeff_cols = SSED_input->coeff_cols;

  SSED->degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
  for (i = 0; i < SSED->size_r; i++)
    total_deg += SSED->degrees[i] = SSED_input->degrees[i];

  set_d(SSED->gamma, SSED_input->gamma);

  if (SSED->startSystemType == 1)
  { // if we are using an mhom start system, we need to copy over the coeff
    SSED->coeff = (comp_d **)bmalloc(total_deg * sizeof(comp_d *));
    for (i = 0; i < total_deg; i++)
    {
      SSED->coeff[i] = (comp_d *)bmalloc(SSED->coeff_cols * sizeof(comp_d));
      for (j = 0; j < SSED->coeff_cols; j++)
      {
        set_d(SSED->coeff[i][j], SSED_input->coeff[i][j]);
      }
    }
  }

  return;
}

void cp_prog_t(prog_t *P, prog_t *P_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of P_input to P                          *
\***************************************************************/
{
  int i;

  P->size = P_input->size;
  P->prog = (int *)bmalloc(P->size * sizeof(int));
  for (i = 0; i < P->size; i++)
  {
    P->prog[i] = P_input->prog[i];
  }
  P->memSize = P_input->memSize; 
  P->numNums = P_input->numNums;
  P->nums = (num_t *)bmalloc(P->numNums * sizeof(num_t));
  for (i = 0; i < P->numNums; i++)
  {
    mpq_init(P->nums[i].rat);
    mpq_set(P->nums[i].rat, P_input->nums[i].rat);
    P->nums[i].currPrec = P_input->nums[i].currPrec;
    mpf_init2(P->nums[i].real, P->nums[i].currPrec);
    mpf_set_q(P->nums[i].real, P->nums[i].rat);
  }
  P->precision = P_input->precision;
  P->num_var_gps = P_input->num_var_gps;
  P->var_gp_sizes = (int *)bmalloc(P->num_var_gps * sizeof(int));
  for (i = 0; i < P->num_var_gps; i++)
  {
    P->var_gp_sizes[i] = P_input->var_gp_sizes[i];
  }
  P->index_of_first_number_for_proj_trans = P_input->index_of_first_number_for_proj_trans;

  P->numInstAtEndUpdate = P_input->numInstAtEndUpdate;
  P->numInstAtEndParams = P_input->numInstAtEndParams;
  P->numInstAtEndFnEval = P_input->numInstAtEndFnEval;
  P->numInstAtEndPDeriv = P_input->numInstAtEndPDeriv;
  P->numInstAtEndJvEval = P_input->numInstAtEndJvEval;

  P->numVars = P_input->numVars;
  P->numPathVars = P_input->numPathVars;
  P->numConsts = P_input->numConsts;
  P->numPars = P_input->numPars;
  P->numFuncs = P_input->numFuncs;
  P->numSubfuncs = P_input->numSubfuncs;

  P->inpVars = P_input->inpVars;
  P->inpPathVars = P_input->inpPathVars;
  P->IAddr = P_input->IAddr;
  P->numAddr = P_input->numAddr;
  P->constAddr = P_input->constAddr;

  P->evalPars = P_input->evalPars;
  P->evalDPars = P_input->evalDPars;
  P->evalFuncs = P_input->evalFuncs;
  P->evalJVars = P_input->evalJVars;
  P->evalJPars = P_input->evalJPars;
  P->evalSubs = P_input->evalSubs;
  P->evalJSubsV = P_input->evalJSubsV;
  P->evalJSubsP = P_input->evalJSubsP;

  return;
}

void cp_basic_eval_data_mp(basic_eval_data_mp *BED, basic_eval_data_mp *BED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of BED_input to BED                      *
\***************************************************************/
{
  cp_preproc_data(&BED->preProcData, &BED_input->preProcData);
  // create a new SLP in BED
  cp_square_system_mp(&BED->squareSystem, &BED_input->squareSystem, 1, NULL);
  cp_patch_mp(&BED->patch, &BED_input->patch);
  cp_start_system_mp(&BED->startSystem, &BED_input->startSystem);

  return;
}

void cp_square_system_mp(square_system_eval_data_mp *SSED, square_system_eval_data_mp *SSED_input, int copy_SLP, prog_t *Prog_ptr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_input to SSED                    *
\***************************************************************/
{
  int i, j;

  if (copy_SLP)
  { // copy SLP
    SSED->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    cp_prog_t(SSED->Prog, SSED_input->Prog);
  }
  else
  { // point to Prog_ptr - used in AMP so only 1 version of the SLP exists
    SSED->Prog = Prog_ptr;
  }

  SSED->size_f = SSED_input->size_f;
  SSED->size_r = SSED_input->size_r;
  SSED->max_of_W = SSED_input->max_of_W;
  SSED->noChanges = SSED_input->noChanges;

  SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
  SSED->P = (int *)bmalloc(SSED->size_f * sizeof(int));
  for (i = 0; i < SSED->size_f; i++)
  {
    SSED->orig_degrees[i] = SSED_input->orig_degrees[i];
    SSED->P[i] = SSED_input->P[i];
  }

  SSED->new_degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
  SSED->W = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
  {
    SSED->new_degrees[i] = SSED_input->new_degrees[i];
    SSED->W[i] = (int *)bmalloc((SSED->size_f - SSED->size_r) * sizeof(int));
    for (j = 0; j < (SSED->size_f - SSED->size_r); j++)
      SSED->W[i][j] = SSED_input->W[i][j];
  }

  // set the current precision
  SSED->curr_prec = SSED_input->curr_prec;

  // initialize all matrices to this preicision
  init_mat_mp2(SSED->B, SSED_input->B->rows, SSED_input->B->cols, SSED->curr_prec);
  init_mat_rat(SSED->B_rat, SSED_input->B->rows, SSED_input->B->cols);
  init_mat_mp2(SSED->B_perp, SSED_input->B_perp->rows, SSED_input->B_perp->cols, SSED->curr_prec);
  init_mat_rat(SSED->B_perp_rat, SSED_input->B_perp->rows, SSED_input->B_perp->cols);
  init_mat_mp2(SSED->A, SSED_input->A->rows, SSED_input->A->cols, SSED->curr_prec);
  init_mat_rat(SSED->A_rat, SSED_input->A->rows, SSED_input->A->cols);

  // setup matrix B
  mat_cp_mp_rat(SSED->B, SSED->B_rat, SSED_input->B, SSED_input->B_rat);  

  // setup matrix B_perp
  mat_cp_mp_rat(SSED->B_perp, SSED->B_perp_rat, SSED_input->B_perp, SSED_input->B_perp_rat);

  // setup matrix A
  mat_cp_mp_rat(SSED->A, SSED->A_rat, SSED_input->A, SSED_input->A_rat);

  return;
}

void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp *PED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PED_input to PED                      *
\***************************************************************/
{
  PED->num_patches = PED_input->num_patches;

  // set the current precision
  PED->curr_prec = PED_input->curr_prec;

  // initialize patchCoeff to this preicision
  init_mat_mp2(PED->patchCoeff, PED_input->patchCoeff->rows, PED_input->patchCoeff->cols, PED->curr_prec);
  init_mat_rat(PED->patchCoeff_rat, PED_input->patchCoeff->rows, PED_input->patchCoeff->cols);

  // setup patchCoeff
  mat_cp_mp_rat(PED->patchCoeff, PED->patchCoeff_rat, PED_input->patchCoeff, PED_input->patchCoeff_rat);

  return;
}

void cp_start_system_mp(start_system_eval_data_mp *SSED, start_system_eval_data_mp *SSED_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_input to SSED                    *
\***************************************************************/
{
  int i, j, total_deg = 0;

  SSED->startSystemType = SSED_input->startSystemType;
  SSED->size_r = SSED_input->size_r;
  SSED->max_degree = SSED_input->max_degree;
  SSED->coeff_cols = SSED_input->coeff_cols;

  SSED->degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
  for (i = 0; i < SSED->size_r; i++)
    total_deg += SSED->degrees[i] = SSED_input->degrees[i];

  // set the current precision
  SSED->curr_prec = SSED_input->curr_prec;

  // setup (mp) gamma 
  init_mp2(SSED->gamma, SSED->curr_prec);
  set_mp(SSED->gamma, SSED_input->gamma);

  // setup (rat) gamma
  SSED->gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
  init_rat(SSED->gamma_rat);
  set_rat(SSED->gamma_rat, SSED_input->gamma_rat);

  if (SSED->startSystemType == 1)
  { // if we are using an mhom start system, we need to copy over the coeff

    SSED->coeff = (comp_mp **)bmalloc(total_deg * sizeof(comp_mp *));
    SSED->coeff_rat = (mpq_t ***)bmalloc(total_deg * sizeof(mpq_t **));
    for (i = 0; i < total_deg; i++)
    {
      SSED->coeff[i] = (comp_mp *)bmalloc(SSED->coeff_cols * sizeof(comp_mp));
      SSED->coeff_rat[i] = (mpq_t **)bmalloc(SSED->coeff_cols * sizeof(mpq_t *));
      for (j = 0; j < SSED->coeff_cols; j++)
      { // setup mp 
        init_mp2(SSED->coeff[i][j], SSED->curr_prec);
        set_mp(SSED->coeff[i][j], SSED_input->coeff[i][j]);       

        // setup rat
        SSED->coeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
        mpq_init(SSED->coeff_rat[i][j][0]);
        mpq_init(SSED->coeff_rat[i][j][1]);
        mpq_set(SSED->coeff_rat[i][j][0], SSED_input->coeff_rat[i][j][0]);
        mpq_set(SSED->coeff_rat[i][j][1], SSED_input->coeff_rat[i][j][1]);
      }
    }
  }

  return;
}

//////////////////// COPY FUNCTIONS FOR MPI THINGS ////////////////////

#ifdef _HAVE_MPI

void cp_tracker_config_t_relevant(void *T_out, void *T_in, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - T_out is _t_relevant, T_in is _t       *
*                otherwise - T_out is _t, T_in is _t_relevant   *
* RETURN VALUES:                                                *
* NOTES: stores a copy of T_in to T_out                         *
\***************************************************************/
{
  if (inType == 0)
  { // _t to _t_relevant so that it can be sent
    tracker_config_t *T_right = (tracker_config_t *)T_in;
    tracker_config_t_relevant *T_left = (tracker_config_t_relevant *)T_out;

    T_left->numVars = T_right->numVars;
    T_left->numPathVars = T_right->numPathVars;
    T_left->numParams = T_right->numParams;
    T_left->numFuncs = T_right->numFuncs;

    T_left->maxStepSize = T_right->maxStepSize;
    T_left->minStepSizeBeforeEndGame = T_right->minStepSizeBeforeEndGame;
    T_left->minStepSizeDuringEndGame = T_right->minStepSizeDuringEndGame;

    T_left->minTrackT = T_right->minTrackT;
    T_left->basicNewtonTol = T_right->basicNewtonTol;
    T_left->endgameNewtonTol = T_right->endgameNewtonTol;

    T_left->cSecInc = T_right->cSecInc;
    T_left->maxNewtonIts = T_right->maxNewtonIts;
    T_left->MPType = T_right->MPType;
    T_left->Precision = T_right->Precision;
    T_left->outputLevel = T_right->outputLevel;
    T_left->screenOut = T_right->screenOut;
    T_left->targetT = T_right->targetT;
    T_left->endgameBoundary = T_right->endgameBoundary;

    T_left->goingToInfinity = T_right->goingToInfinity;
    T_left->maxNumSteps = T_right->maxNumSteps;
    T_left->endgameNumber = T_right->endgameNumber;

    T_left->power_series_sample_factor = T_right->power_series_sample_factor;
    T_left->cycle_num_max = T_right->cycle_num_max;
    T_left->num_PSEG_sample_points = T_right->num_PSEG_sample_points;

    T_left->final_tolerance = T_right->final_tolerance;
    T_left->real_threshold = T_right->real_threshold;
    T_left->endgameOnly = T_right->endgameOnly;

    T_left->AMP_bound_on_abs_vals_of_coeffs = T_right->AMP_bound_on_abs_vals_of_coeffs;
    T_left->AMP_bound_on_degree = T_right->AMP_bound_on_degree;
    T_left->AMP_eps = T_right->AMP_eps;
    T_left->AMP_Phi = T_right->AMP_Phi;
    T_left->AMP_Psi = T_right->AMP_Psi;
    T_left->AMP_safety_digits_1 = T_right->AMP_safety_digits_1;
    T_left->AMP_safety_digits_2 = T_right->AMP_safety_digits_2;
    T_left->AMP_max_prec = T_right->AMP_max_prec;

    T_left->sing_val_zero_tol = T_right->sing_val_zero_tol;
    T_left->cond_num_threshold = T_right->cond_num_threshold;

    T_left->step_fail_factor = T_right->step_fail_factor;
    T_left->step_success_factor = T_right->step_success_factor;

    T_left->max_num_pts_for_trace = T_right->max_num_pts_for_trace;
    T_left->max_num_mon_linears = T_right->max_num_mon_linears;
    T_left->max_num_bad_loops_in_mon = T_right->max_num_bad_loops_in_mon;

    T_left->final_tol_multiplier = T_right->final_tol_multiplier;
    T_left->final_tol_times_mult = T_right->final_tol_times_mult;

    T_left->sharpenOnly = T_right->sharpenOnly;
    T_left->sharpenDigits = T_right->sharpenDigits;

    T_left->regen_remove_inf = T_right->regen_remove_inf;
    T_left->regen_higher_dim_check = T_right->regen_higher_dim_check;
    T_left->sliceBasicNewtonTol = T_right->sliceBasicNewtonTol;
    T_left->sliceEndgameNewtonTol = T_right->sliceEndgameNewtonTol;
    T_left->sliceFinalTol = T_right->sliceFinalTol;

    T_left->minCycleTrackBack = T_right->minCycleTrackBack;
    T_left->junkRemovalTest = T_right->junkRemovalTest;
    T_left->maxDepthLDT = T_right->maxDepthLDT;
    T_left->odePredictor = T_right->odePredictor;

    T_left->securityLevel = T_right->securityLevel;
    T_left->securityMaxNorm = T_right->securityMaxNorm;

    T_left->cutoffCycleTime = T_right->cutoffCycleTime;
    T_left->cutoffRatioTime = T_right->cutoffRatioTime;
    T_left->finiteThreshold = T_right->finiteThreshold;
    T_left->funcResTol = T_right->funcResTol;
    T_left->ratioTol = T_right->ratioTol;
    T_left->maxStepsBeforeNewton = T_right->maxStepsBeforeNewton;
  }
  else
  { // _t_relevant to _t so that it can be used by the workers
    tracker_config_t *T_left = (tracker_config_t *)T_out;
    tracker_config_t_relevant *T_right = (tracker_config_t_relevant *)T_in;

    T_left->numVars = T_right->numVars;
    T_left->numPathVars = T_right->numPathVars;
    T_left->numParams = T_right->numParams;
    T_left->numFuncs = T_right->numFuncs;

    T_left->maxStepSize = T_right->maxStepSize;
    T_left->minStepSizeBeforeEndGame = T_right->minStepSizeBeforeEndGame;
    T_left->minStepSizeDuringEndGame = T_right->minStepSizeDuringEndGame;

    T_left->minTrackT = T_right->minTrackT;
    T_left->basicNewtonTol = T_right->basicNewtonTol;
    T_left->endgameNewtonTol = T_right->endgameNewtonTol;

    T_left->cSecInc = T_right->cSecInc;
    T_left->maxNewtonIts = T_right->maxNewtonIts;
    T_left->MPType = T_right->MPType;
    T_left->Precision = T_right->Precision;
    T_left->outputLevel = T_right->outputLevel;
    T_left->screenOut = T_right->screenOut;
    T_left->targetT = T_right->targetT;
    T_left->endgameBoundary = T_right->endgameBoundary;

    T_left->goingToInfinity = T_right->goingToInfinity;
    T_left->maxNumSteps = T_right->maxNumSteps;
    T_left->endgameNumber = T_right->endgameNumber;

    T_left->power_series_sample_factor = T_right->power_series_sample_factor;
    T_left->cycle_num_max = T_right->cycle_num_max;
    T_left->num_PSEG_sample_points = T_right->num_PSEG_sample_points;

    T_left->final_tolerance = T_right->final_tolerance;
    T_left->real_threshold = T_right->real_threshold;
    T_left->endgameOnly = T_right->endgameOnly;

    T_left->AMP_bound_on_abs_vals_of_coeffs = T_right->AMP_bound_on_abs_vals_of_coeffs;
    T_left->AMP_bound_on_degree = T_right->AMP_bound_on_degree;
    T_left->AMP_eps = T_right->AMP_eps;
    T_left->AMP_Phi = T_right->AMP_Phi;
    T_left->AMP_Psi = T_right->AMP_Psi;
    T_left->AMP_safety_digits_1 = T_right->AMP_safety_digits_1;
    T_left->AMP_safety_digits_2 = T_right->AMP_safety_digits_2;
    T_left->AMP_max_prec = T_right->AMP_max_prec;

    T_left->sing_val_zero_tol = T_right->sing_val_zero_tol;
    T_left->cond_num_threshold = T_right->cond_num_threshold;

    T_left->step_fail_factor = T_right->step_fail_factor;
    T_left->step_success_factor = T_right->step_success_factor;

    T_left->max_num_pts_for_trace = T_right->max_num_pts_for_trace;
    T_left->max_num_mon_linears = T_right->max_num_mon_linears;
    T_left->max_num_bad_loops_in_mon = T_right->max_num_bad_loops_in_mon;

    T_left->final_tol_multiplier = T_right->final_tol_multiplier;
    T_left->final_tol_times_mult = T_right->final_tol_times_mult;

    T_left->sharpenOnly = T_right->sharpenOnly;
    T_left->sharpenDigits = T_right->sharpenDigits;

    T_left->regen_remove_inf = T_right->regen_remove_inf;
    T_left->regen_higher_dim_check = T_right->regen_higher_dim_check;
    T_left->sliceBasicNewtonTol = T_right->sliceBasicNewtonTol;
    T_left->sliceEndgameNewtonTol = T_right->sliceEndgameNewtonTol;
    T_left->sliceFinalTol = T_right->sliceFinalTol;

    T_left->minCycleTrackBack = T_right->minCycleTrackBack;
    T_left->junkRemovalTest = T_right->junkRemovalTest;
    T_left->maxDepthLDT = T_right->maxDepthLDT;
    T_left->odePredictor = T_right->odePredictor;

    T_left->securityLevel = T_right->securityLevel;
    T_left->securityMaxNorm = T_right->securityMaxNorm;

    T_left->cutoffCycleTime = T_right->cutoffCycleTime;
    T_left->cutoffRatioTime = T_right->cutoffRatioTime;
    T_left->finiteThreshold = T_right->finiteThreshold;
    T_left->funcResTol = T_right->funcResTol;
    T_left->ratioTol = T_right->ratioTol;
    T_left->maxStepsBeforeNewton = T_right->maxStepsBeforeNewton;

    // initialize the rest of T_left
    T_left->latest_newton_residual_d = 0;
    T_left->t_val_at_latest_sample_point = 0;
    T_left->error_at_latest_sample_point = 0;
    if (T_left->MPType == 1)
    { // fixed MP
      mpf_init2(T_left->latest_newton_residual_mp, T_left->Precision);
      mpf_set_ui(T_left->latest_newton_residual_mp, 0);
    }
    else if (T_left->MPType == 2)
    { // AMP
      mpf_init2(T_left->latest_newton_residual_mp, T_left->AMP_max_prec);
      mpf_set_ui(T_left->latest_newton_residual_mp, 0);
    }
  }

  return;
}

void cp_prog_str(void *Out, void *In, int *length, int numNums, int prec, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is char, In is num_t               *
*                otherwise - Out is num_t, In is char           *
* RETURN VALUES:                                                *
* NOTES: stores num_t to string or string to num_t              *
\***************************************************************/
{
  int i, strLength, size, base = 10;
  char *temp = NULL;

  if (inType == 0)
  { // num_t to char so that it can be sent
    num_t **Right = (num_t **)In;
    char **Left = (char **)Out;

    // initialize
    *length = 0;
    *Left = NULL;

    // setup Left
    for (i = 0; i < numNums; i++)
    { // convert mpq to string
      temp = mpq_get_str(NULL, base, (*Right)[i].rat); 
      strLength = strlen(temp) + 1; // +1 for '\0'
      // reallocate enough space
      *Left = (char *)brealloc(*Left, (*length + strLength) * sizeof(char));
      // copy on to end of *Left
      strcpy(&(*Left)[*length], temp);
      // update the length
      *length = *length + strLength; 
      // free temp
      free(temp);
      temp = NULL;
    }   
  }
  else
  { // char to num_t so that it can be used by the workers
    num_t **Left = (num_t **)Out;
    char **Right = (char **)In;

    // initialize
    *Left = (num_t *)bmalloc(numNums * sizeof(num_t));
    strLength = 0;

    // setup each num
    for (i = 0; i < numNums; i++)
    { // setup rat
      mpq_init((*Left)[i].rat);
      mpq_set_str((*Left)[i].rat, &(*Right)[strLength], base);
      mpq_canonicalize((*Left)[i].rat);
      // setup real
      (*Left)[i].currPrec = prec;
      mpf_init2((*Left)[i].real, prec);
      mpf_set_q((*Left)[i].real, (*Left)[i].rat);
      // update strLength
      size = strlen(&(*Right)[strLength]);
      strLength += size + 1;
    }
    // clear Right
    if (freeStr)
      free(*Right);
  }

  return;
}

void cp_prog_t_int(void *Prog_out, void *Prog_in, int **instructions, int **gp_sizes, char **progStr, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Prog_out is _t_int, Prog_in is _t      *
*                otherwise - Prog_out is _t, Prog_in is _t_int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of Prog_in to Prog_out                   *
*  either sets up instructions, gp_sizes and progStr or use them*
\***************************************************************/
{
  int i;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    prog_t *Prog_right = (prog_t *)Prog_in;
    prog_t_int *Prog_left = (prog_t_int *)Prog_out;

    Prog_left->size = Prog_right->size;
    Prog_left->memSize = Prog_right->memSize;
    Prog_left->precision = Prog_right->precision;

    Prog_left->num_var_gps = Prog_right->num_var_gps;
    Prog_left->index_of_first_number_for_proj_trans = Prog_right->index_of_first_number_for_proj_trans;

    Prog_left->numInstAtEndUpdate = Prog_right->numInstAtEndUpdate;
    Prog_left->numInstAtEndParams = Prog_right->numInstAtEndParams;
    Prog_left->numInstAtEndPDeriv = Prog_right->numInstAtEndPDeriv;
    Prog_left->numInstAtEndFnEval = Prog_right->numInstAtEndFnEval;
    Prog_left->numInstAtEndJvEval = Prog_right->numInstAtEndJvEval;

    Prog_left->numVars = Prog_right->numVars;
    Prog_left->numPathVars = Prog_right->numPathVars;
    Prog_left->numNums = Prog_right->numNums;
    Prog_left->numConsts = Prog_right->numConsts;

    Prog_left->numPars = Prog_right->numPars;
    Prog_left->numFuncs = Prog_right->numFuncs;
    Prog_left->numSubfuncs = Prog_right->numSubfuncs;

    Prog_left->inpVars = Prog_right->inpVars;
    Prog_left->inpPathVars = Prog_right->inpPathVars;
    Prog_left->IAddr = Prog_right->IAddr;
    Prog_left->numAddr = Prog_right->numAddr;
    Prog_left->constAddr = Prog_right->constAddr;

    Prog_left->evalPars = Prog_right->evalPars;
    Prog_left->evalDPars = Prog_right->evalDPars;
    Prog_left->evalFuncs = Prog_right->evalFuncs;
    Prog_left->evalJVars = Prog_right->evalJVars;
    Prog_left->evalJPars = Prog_right->evalJPars;
    Prog_left->evalSubs = Prog_right->evalSubs;
    Prog_left->evalJSubsV = Prog_right->evalJSubsV;
    Prog_left->evalJSubsP = Prog_right->evalJSubsP;

    // setup instructions
    *instructions = (int *)bmalloc(Prog_right->size * sizeof(int));
    for (i = Prog_right->size - 1; i >= 0; i--)
      (*instructions)[i] = Prog_right->prog[i];

    // setup gp_sizes
    *gp_sizes = (int *)bmalloc(Prog_right->num_var_gps * sizeof(int));
    for (i = Prog_right->num_var_gps - 1; i >= 0; i--)
      (*gp_sizes)[i] = Prog_right->var_gp_sizes[i];

    // setup progStr
    cp_prog_str(progStr, &Prog_right->nums, &Prog_left->totalLength, Prog_right->numNums, Prog_right->precision, freeStr, inType);
  }
  else
  { // _t_int to _t so that it can be used by the workers
    prog_t *Prog_left = (prog_t *)Prog_out;
    prog_t_int *Prog_right = (prog_t_int *)Prog_in;

    Prog_left->size = Prog_right->size;
    Prog_left->memSize =  Prog_right->memSize;
    Prog_left->precision = Prog_right->precision;

    Prog_left->num_var_gps = Prog_right->num_var_gps;
    Prog_left->index_of_first_number_for_proj_trans = Prog_right->index_of_first_number_for_proj_trans;

    Prog_left->numInstAtEndUpdate = Prog_right->numInstAtEndUpdate;
    Prog_left->numInstAtEndParams = Prog_right->numInstAtEndParams;
    Prog_left->numInstAtEndPDeriv = Prog_right->numInstAtEndPDeriv;
    Prog_left->numInstAtEndFnEval = Prog_right->numInstAtEndFnEval;
    Prog_left->numInstAtEndJvEval = Prog_right->numInstAtEndJvEval;

    Prog_left->numVars = Prog_right->numVars;
    Prog_left->numPathVars = Prog_right->numPathVars;
    Prog_left->numNums = Prog_right->numNums;
    Prog_left->numConsts = Prog_right->numConsts;

    Prog_left->numPars = Prog_right->numPars;
    Prog_left->numFuncs = Prog_right->numFuncs;
    Prog_left->numSubfuncs = Prog_right->numSubfuncs;

    Prog_left->inpVars = Prog_right->inpVars;
    Prog_left->inpPathVars = Prog_right->inpPathVars;
    Prog_left->IAddr = Prog_right->IAddr;
    Prog_left->numAddr = Prog_right->numAddr;
    Prog_left->constAddr = Prog_right->constAddr;

    Prog_left->evalPars = Prog_right->evalPars;
    Prog_left->evalDPars = Prog_right->evalDPars;
    Prog_left->evalFuncs = Prog_right->evalFuncs;
    Prog_left->evalJVars = Prog_right->evalJVars;
    Prog_left->evalJPars = Prog_right->evalJPars;
    Prog_left->evalSubs = Prog_right->evalSubs;
    Prog_left->evalJSubsV = Prog_right->evalJSubsV;
    Prog_left->evalJSubsP = Prog_right->evalJSubsP;

    // setup prog
    Prog_left->prog = (int *)bmalloc(Prog_right->size * sizeof(int));
    for (i = Prog_right->size - 1; i >= 0; i--)
      Prog_left->prog[i] = (*instructions)[i];
    // clear instructions
    free(*instructions);

    // setup var_gp_sizes
    Prog_left->var_gp_sizes = (int *)bmalloc(Prog_right->num_var_gps * sizeof(int));
    for (i = Prog_right->num_var_gps - 1; i >= 0; i--)
      Prog_left->var_gp_sizes[i] = (*gp_sizes)[i];
    // clear gp_sizes
    free(*gp_sizes);

    // setup progStr
    cp_prog_str(&Prog_left->nums, progStr, &Prog_right->totalLength, Prog_right->numNums, Prog_right->precision, freeStr, inType);
  }

  return;
}

void cp_square_system_eval_data_d_int(void *SSED_out, void *SSED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int freeStr, int **orig_deg, int **new_deg, int **P, int **W, comp_d **coeff, int MPType, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - SSED_out is _d_int, SSED_in is _d      *
*                otherwise - SSED_out is _d, SSED_in is _d_int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_in to SSED_out                   *
\***************************************************************/
{
  int i, j, num_cols;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    square_system_eval_data_d *SSED_right = (square_system_eval_data_d *)SSED_in;
    square_system_eval_data_d_int *SSED_left = (square_system_eval_data_d_int *)SSED_out;

    // copy Prog
    cp_prog_t_int(&SSED_left->Prog_int, SSED_right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);

    // setup orig_deg
    *orig_deg = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      (*orig_deg)[i] = SSED_right->orig_degrees[i];

    // setup new_deg
    *new_deg = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      (*new_deg)[i] = SSED_right->new_degrees[i];

    // setup P
    *P = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      (*P)[i] = SSED_right->P[i];

    // setup W
    num_cols = SSED_right->size_f - SSED_right->size_r;
    *W = (int *)bmalloc(SSED_right->size_r * num_cols * sizeof(int));
    for (i = 0; i < SSED_right->size_r; i++)
      for (j = 0; j < num_cols; j++)
        (*W)[num_cols * i + j] = SSED_right->W[i][j];

    // setup B, B_perp & A
    SSED_left->B_rows = SSED_right->B->rows;
    SSED_left->B_cols = SSED_right->B->cols;
    SSED_left->B_perp_rows = SSED_right->B_perp->rows;
    SSED_left->B_perp_cols = SSED_right->B_perp->cols;
    SSED_left->A_rows = SSED_right->A->rows;
    SSED_left->A_cols = SSED_right->A->cols;

    // setup coeff
    SSED_left->num_comp_d = SSED_left->B_rows * SSED_left->B_cols + SSED_left->B_perp_rows * SSED_left->B_perp_cols + SSED_left->A_rows * SSED_left->A_cols;
    *coeff = (comp_d *)bmalloc(SSED_left->num_comp_d * sizeof(comp_d));
    num_cols = 0;
    for (i = 0; i < SSED_left->B_rows; i++)
      for (j = 0; j < SSED_left->B_cols; j++)
      {
        set_d((*coeff)[num_cols], &SSED_right->B->entry[i][j]);
        num_cols++;
      }
    for (i = 0; i < SSED_left->B_perp_rows; i++)
      for (j = 0; j < SSED_left->B_perp_cols; j++)
      {
        set_d((*coeff)[num_cols], &SSED_right->B_perp->entry[i][j]);
        num_cols++;
      }
    for (i = 0; i < SSED_left->A_rows; i++)
      for (j = 0; j < SSED_left->A_cols; j++)
      {
        set_d((*coeff)[num_cols], &SSED_right->A->entry[i][j]);
        num_cols++;
      }

    // setup the rest
    SSED_left->size_f = SSED_right->size_f;
    SSED_left->noChanges = SSED_right->noChanges;
    SSED_left->max_of_W = SSED_right->max_of_W;
    SSED_left->size_r = SSED_right->size_r;
  }
  else
  { // _d_int to _d so that it can be used by the workers
    square_system_eval_data_d *SSED_left = (square_system_eval_data_d *)SSED_out;
    square_system_eval_data_d_int *SSED_right = (square_system_eval_data_d_int *)SSED_in;

    // allocate memory for Prog
    SSED_left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t)); 

    // copy Prog
    cp_prog_t_int(SSED_left->Prog, &SSED_right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    // initialize evalProg
    initEvalProg(MPType);

    // setup orig_degrees
    SSED_left->orig_degrees = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      SSED_left->orig_degrees[i] = (*orig_deg)[i];
    // clear orig_deg
    free(*orig_deg);

    // setup new_degrees
    SSED_left->new_degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      SSED_left->new_degrees[i] = (*new_deg)[i];
    // clear new_deg
    free(*new_deg);

    // setup P
    SSED_left->P = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      SSED_left->P[i] = (*P)[i];
    // clear P
    free(*P);

    // setup W
    num_cols = SSED_right->size_f - SSED_right->size_r;
    SSED_left->W = (int **)bmalloc(SSED_right->size_r * sizeof(int *));
    for (i = 0; i < SSED_right->size_r; i++)
    {
      SSED_left->W[i] = (int *)bmalloc(num_cols * sizeof(int));
      for (j = 0; j < num_cols; j++)
        SSED_left->W[i][j] = (*W)[num_cols * i + j];
    }
    // clear W
    free(*W);

    num_cols = 0;
    // setup B
    init_mat_d(SSED_left->B, SSED_right->B_rows, SSED_right->B_cols);
    SSED_left->B->rows = SSED_right->B_rows;
    SSED_left->B->cols = SSED_right->B_cols;
    for (i = 0; i < SSED_right->B_rows; i++)
      for (j = 0; j < SSED_right->B_cols; j++)
      {
        set_d(&SSED_left->B->entry[i][j], (*coeff)[num_cols]);
        num_cols++;
      }

    // setup B_perp
    init_mat_d(SSED_left->B_perp, SSED_right->B_perp_rows, SSED_right->B_perp_cols);
    SSED_left->B_perp->rows = SSED_right->B_perp_rows;
    SSED_left->B_perp->cols = SSED_right->B_perp_cols;
    for (i = 0; i < SSED_right->B_perp_rows; i++)
      for (j = 0; j < SSED_right->B_perp_cols; j++)
      {
        set_d(&SSED_left->B_perp->entry[i][j], (*coeff)[num_cols]);
        num_cols++;
      }

    // setup A
    init_mat_d(SSED_left->A, SSED_right->A_rows, SSED_right->A_cols);
    SSED_left->A->rows = SSED_right->A_rows;
    SSED_left->A->cols = SSED_right->A_cols;
    for (i = 0; i < SSED_right->A_rows; i++)
      for (j = 0; j < SSED_right->A_cols; j++)
      {
        set_d(&SSED_left->A->entry[i][j], (*coeff)[num_cols]);
        num_cols++;
      }

    // clear coeff
    free(*coeff);

    // setup the rest
    SSED_left->size_f = SSED_right->size_f;
    SSED_left->noChanges = SSED_right->noChanges;
    SSED_left->max_of_W = SSED_right->max_of_W;
    SSED_left->size_r = SSED_right->size_r;
  }


  return;
}

void cp_start_system_eval_data_d_int(void *SSED_out, void *SSED_in, int **degrees, comp_d **coeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - SSED_out is _d_int, SSED_in is _d      *
*                otherwise - SSED_out is _d, SSED_in is _d_int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_in to SSED_out                   *
\***************************************************************/
{
  int i, j, num_cols, totalDeg = 0;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    start_system_eval_data_d *SSED_right = (start_system_eval_data_d *)SSED_in;
    start_system_eval_data_d_int *SSED_left = (start_system_eval_data_d_int *)SSED_out;

    // setup degrees
    *degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      totalDeg += (*degrees)[i] = SSED_right->degrees[i];

    // setup coeff - if needed
    if (SSED_right->startSystemType == 1)
    {
      num_cols = SSED_right->coeff_cols;
      SSED_left->num_coeff = totalDeg * num_cols;
      *coeff = (comp_d *)bmalloc(SSED_left->num_coeff * sizeof(comp_d));
      for (i = 0; i < totalDeg; i++)
        for (j = 0; j < num_cols; j++)
        {
          set_d((*coeff)[num_cols * i + j] , SSED_right->coeff[i][j]);
        }
    }
    else
    {
      *coeff = NULL;
      SSED_left->num_coeff = 0;
    }
    
    // setup the rest
    SSED_left->startSystemType = SSED_right->startSystemType;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->max_degree = SSED_right->max_degree;
    SSED_left->coeff_cols = SSED_right->coeff_cols;
    set_d(SSED_left->gamma, SSED_right->gamma);
  }
  else
  { // _d_int to _d so that it can be used by the workers
    start_system_eval_data_d *SSED_left = (start_system_eval_data_d *)SSED_out;
    start_system_eval_data_d_int *SSED_right = (start_system_eval_data_d_int *)SSED_in;

    // setup degrees
    SSED_left->degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      totalDeg += SSED_left->degrees[i] = (*degrees)[i];
    // clear degrees
    free(*degrees);

    // setup coeff - if needed
    if (SSED_right->startSystemType == 1)
    {
       num_cols = SSED_right->coeff_cols;
       SSED_left->coeff = (comp_d **)bmalloc(totalDeg * sizeof(comp_d *));
       for (i = 0; i < totalDeg; i++)
       {
         SSED_left->coeff[i] = (comp_d *)bmalloc(SSED_right->coeff_cols * sizeof(comp_d));
         for (j = 0; j < num_cols; j ++)
         {
           set_d(SSED_left->coeff[i][j], (*coeff)[num_cols * i + j]);
         }
       }
       // clear coeff
       free(*coeff);
    }

    // setup the rest
    SSED_left->startSystemType = SSED_right->startSystemType;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->max_degree = SSED_right->max_degree;
    SSED_left->coeff_cols = SSED_right->coeff_cols;
    set_d(SSED_left->gamma, SSED_right->gamma);
  }

  return;
}

void cp_patch_d_int(void *PED_out, void *PED_in, comp_d **coeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - PED_out is _d_int, PED_in is _d        *
*                otherwise - PED_out is _d, PED_in is _d_int    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PED_in to PED_out                     *
\***************************************************************/
{
  int i, j, num;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    patch_eval_data_d *PED_right = (patch_eval_data_d *)PED_in;
    patch_eval_data_d_int *PED_left = (patch_eval_data_d_int *)PED_out;

    // setup sizes
    PED_left->num_patches = PED_right->num_patches;
    PED_left->patchCoeff_rows = PED_right->patchCoeff->rows;
    PED_left->patchCoeff_cols = PED_right->patchCoeff->cols;

    // setup coeff
    *coeff = (comp_d *)bmalloc(PED_left->patchCoeff_rows * PED_left->patchCoeff_cols * sizeof(comp_d));
    num = 0;
    for (i = 0; i < PED_left->patchCoeff_rows; i++)
      for (j = 0; j < PED_left->patchCoeff_cols; j++)
      {
        set_d((*coeff)[num], &PED_right->patchCoeff->entry[i][j]);
        num++;
      }
  }
  else
  { // _d_int to _d so that it can be used by the workers
    patch_eval_data_d *PED_left = (patch_eval_data_d *)PED_out;
    patch_eval_data_d_int *PED_right = (patch_eval_data_d_int *)PED_in;

    // initialize patchCoeff
    init_mat_d(PED_left->patchCoeff, PED_right->patchCoeff_rows, PED_right->patchCoeff_cols);

    // setup patchCoeff
    PED_left->patchCoeff->rows = PED_right->patchCoeff_rows;
    PED_left->patchCoeff->cols = PED_right->patchCoeff_cols;
    num = 0;
    for (i = 0; i < PED_right->patchCoeff_rows; i++)
      for (j = 0; j < PED_right->patchCoeff_cols; j++)
      {
        set_d(&PED_left->patchCoeff->entry[i][j], (*coeff)[num]);
        num++;
      }

    // setup the rest
    PED_left->num_patches = PED_right->num_patches;

    // clear coeff
    free(*coeff);
  }

  return;
}

void cp_preproc_data_int(void *PPD_out, void *PPD_in, int **type, int **size, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - PPD_out is _int, PPD_in is _           *
*                otherwise - PPD_out is _, PPD_in is _int       
* RETURN VALUES:                                                *
* NOTES: stores a copy of PPD_in to PPD_out                     *
\***************************************************************/
{
  int i, totalSize = 0;

  if (inType == 0)
  { // _t to _int so that it can be sent
    preproc_data *PPD_right = (preproc_data *)PPD_in;
    preproc_data_int *PPD_left = (preproc_data_int *)PPD_out;

    totalSize = PPD_right->num_hom_var_gp + PPD_right->num_var_gp;

    // setup type & size
    *type = (int *)bmalloc(totalSize * sizeof(int));
    *size = (int *)bmalloc(totalSize * sizeof(int));
    for (i = totalSize - 1; i >= 0; i--)
    {
      (*type)[i] = PPD_right->type[i];
      (*size)[i] = PPD_right->size[i];
    }

    // setup the rest
    PPD_left->num_funcs = PPD_right->num_funcs;
    PPD_left->num_hom_var_gp = PPD_right->num_hom_var_gp;
    PPD_left->num_var_gp = PPD_right->num_var_gp;
  }
  else
  { // _int to _ so that it can be used by the workers
    preproc_data *PPD_left = (preproc_data *)PPD_out;
    preproc_data_int *PPD_right = (preproc_data_int *)PPD_in;

    totalSize = PPD_right->num_hom_var_gp + PPD_right->num_var_gp;

    // setup type & size
    PPD_left->type = (int *)bmalloc(totalSize * sizeof(int));
    PPD_left->size = (int *)bmalloc(totalSize * sizeof(int));
    for (i = totalSize - 1; i >= 0; i--)
    {
      PPD_left->type[i] = (*type)[i];
      PPD_left->size[i] = (*size)[i];
    }
    // clear type & size
    free(*type);
    free(*size);

    // setup the rest
    PPD_left->num_funcs = PPD_right->num_funcs;
    PPD_left->num_hom_var_gp = PPD_right->num_hom_var_gp;
    PPD_left->num_var_gp = PPD_right->num_var_gp;
  }
 
  return;
}

void cp_basic_eval_data_d_int(void *BED_out, void *BED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int freeStr, int **orig_deg, int **new_deg, int **P, int **W, int **startDeg, comp_d **st_coeff, comp_d **sq_coeff, comp_d **patch_coeff, int **ppd_type, int **ppd_size, int MPType, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - BED_out is _d_int, BED_in is _d        *
*                otherwise - BED_out is _d, BED_in is _d_int    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of BED_in to BED_out                     *
\***************************************************************/
{
  if (inType == 0)
  { // _d to _d_int so that it can be sent
    basic_eval_data_d *BED_right = (basic_eval_data_d *)BED_in;
    basic_eval_data_d_int *BED_left = (basic_eval_data_d_int *)BED_out;

    // copy square_system
    cp_square_system_eval_data_d_int(&BED_left->squareSystem_int, &BED_right->squareSystem, prog_inst, prog_gp_sizes, progStr, freeStr, orig_deg, new_deg, P, W, sq_coeff, MPType, inType);
    // copy patch
    cp_patch_d_int(&BED_left->patch, &BED_right->patch, patch_coeff, inType);
    // copy start_system
    cp_start_system_eval_data_d_int(&BED_left->startSystem_int, &BED_right->startSystem, startDeg, st_coeff, inType);
    // copy preproc_data
    cp_preproc_data_int(&BED_left->preProcData_int, &BED_right->preProcData, ppd_type, ppd_size, inType);
    // setup the MPType
    BED_left->MPType = MPType;
  }
  else
  { // _d_int to _d so that it can be used by the workers
    basic_eval_data_d *BED_left = (basic_eval_data_d *)BED_out;
    basic_eval_data_d_int *BED_right = (basic_eval_data_d_int *)BED_in;

    // copy square system
    cp_square_system_eval_data_d_int(&BED_left->squareSystem, &BED_right->squareSystem_int, prog_inst, prog_gp_sizes, progStr, freeStr, orig_deg, new_deg, P, W, sq_coeff, MPType, inType);

    // copy patch
    cp_patch_d_int(&BED_left->patch, &BED_right->patch, patch_coeff, inType);
    // copy start system
    cp_start_system_eval_data_d_int(&BED_left->startSystem, &BED_right->startSystem_int, startDeg, st_coeff, inType);
    // copy preproc_data
    cp_preproc_data_int(&BED_left->preProcData, &BED_right->preProcData_int, ppd_type, ppd_size, inType);

    if (BED_right->MPType == 2)
    { // need BED_mp allocated
      BED_left->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
    }
  }

  return;
}

////////////////////// MP FUNCTIONS ////////////////////////////////

void cp_mat_mp_int(void *Out, void *In, char **str, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _int, In is _mp                 *
*                otherwise - Out is _mp, In is _int             *
* RETURN VALUES:                                                *
* NOTES: copies mat_mp to string and _int                       *
\***************************************************************/
{
  int i, j, strLoc, rows, cols, sizeR, sizeI, base = 10;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // mat_mp to string and _int so that it can be sent
    _mat_mp *Right = (_mat_mp *)In;
    mat_mp_int *Left = (mat_mp_int *)Out;

    // setup Left
    rows = Left->rows = (*Right).rows;
    cols = Left->cols = (*Right).cols;
    Left->prec = mpf_get_prec((*Right).entry[0][0].r);

    // initialize
    *str = NULL;
    strLoc = 0;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      { // real part
        strR = mpf_to_str((*Right).entry[i][j].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str((*Right).entry[i][j].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update str
        *str = (char *)brealloc(*str, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*str)[strLoc], strR);
        strcpy(&(*str)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
    // set totalLength
    Left->totalLength = strLoc;
  }
  else
  { // string and _int to mat_mp
    _mat_mp *Left = (_mat_mp *)Out;
    mat_mp_int *Right = (mat_mp_int *)In;

    rows = (*Left).rows = Right->rows;
    cols = (*Left).cols = Right->cols;

    strLoc = 0;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      { // setup real part
        mpf_set_str((*Left).entry[i][j].r, &(*str)[strLoc], base);
        strLoc += 1 + strlen(&(*str)[strLoc]);
        // setup imag part
        mpf_set_str((*Left).entry[i][j].i, &(*str)[strLoc], base);
        strLoc += 1 + strlen(&(*str)[strLoc]);
      }
    if (freeStr)
      free(*str);
  }
 
  return;
}

void cp_mat_rat_int(void *Out, void *In, char **str, int rows, int cols, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _int, In is _rat                *
*                otherwise - Out is _rat, In is _int            *
* RETURN VALUES:                                                *
* NOTES: copies mpq_t *** to string and _int                    *
\***************************************************************/
{
  int i, j, strLoc, sizeR, sizeI, base = 10;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // mat_rat to string and _int so that it can be sent
    mpq_t ***Right = (mpq_t ***)In;
    mat_rat_int *Left = (mat_rat_int *)Out;

    // setup Left
    Left->rows = rows;
    Left->cols = cols;

    // initialize
    *str = NULL;
    strLoc = 0;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      { // real part
        strR = mpq_get_str(NULL, base, Right[i][j][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right[i][j][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update str
        *str = (char *)brealloc(*str, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*str)[strLoc], strR);
        strcpy(&(*str)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
    // set totalLength
    Left->totalLength = strLoc;
  }
  else
  { // string & _int to mat_rat
    mpq_t ***Left = (mpq_t ***)Out;

    strLoc = 0;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      { // setup real part
        mpq_set_str(Left[i][j][0], &(*str)[strLoc], base);
        mpq_canonicalize(Left[i][j][0]);
        strLoc += 1 + strlen(&(*str)[strLoc]);
        // setup imag part
        mpq_set_str(Left[i][j][1], &(*str)[strLoc], base);
        mpq_canonicalize(Left[i][j][1]);
        strLoc += 1 + strlen(&(*str)[strLoc]);
      }
    if (freeStr)
      free(*str);
  }

  return;
}

void cp_vec_rat_char(void *Out, void *In, int *length, int size, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is char, In is vec_rat             *
*                otherwise - Out is vec_rat, In is char         *
* RETURN VALUES:                                                *
* NOTES: stores vec_rat to string or string to vec_rat          *
\***************************************************************/
{
  int i, k, strLength, strSize, base = 10;
  char *temp = NULL;

  if (inType == 0)
  { // vec_rat to char so that it can be sent
    mpq_t ***Right = (mpq_t ***)In;
    char **Left = (char **)Out;

    // initialize
    *length = 0;
    *Left = NULL;

    // setup Left
    for (i = 0; i < size; i++)
      for (k = 0; k < 2; k++) // for real & imag
      { // convert mpq to string
        temp = mpq_get_str(NULL, base, (*Right)[i][k]);
        strLength = strlen(temp) + 1; // +1 for '\0'
        // reallocate enough space
        *Left = (char *)brealloc(*Left, (*length + strLength) * sizeof(char));
        // copy on to end of *Left
        strcpy(&(*Left)[*length], temp);
        // update the length
        *length = *length + strLength;
        // free temp
        free(temp);
        temp = NULL;
      }
  }
  else
  { // char to mpq_t so that it can be used by the workers
    mpq_t ***Left = (mpq_t ***)Out;
    char **Right = (char **)In;

    // initialize
    strLength = 0;

    // setup Left
    *Left = (mpq_t **)bmalloc(size * sizeof(mpq_t *));
    for (i = 0; i < size; i++)
    {
      (*Left)[i] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
      for (k = 0; k < 2; k++) // for real & imag
      { // setup rat
        mpq_init((*Left)[i][k]);
        mpq_set_str((*Left)[i][k], &(*Right)[strLength], base);
        mpq_canonicalize((*Left)[i][k]);
        // update strLength
        strSize = strlen(&(*Right)[strLength]);
        strLength += strSize + 1;
      }
    }
    // free Right
    if (freeStr)
      free(*Right);
  }

  return;
}

void cp_mat_rat_char(void *Out, void *In, int *length, int rows, int cols, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is char, In is mat_rat             *
*                otherwise - Out is mat_rat, In is char         *
* RETURN VALUES:                                                *
* NOTES: stores mat_rat to string or string to mat_rat          *
\***************************************************************/
{
  int i, j, k, strLength, size, base = 10;
  char *temp = NULL;

  if (inType == 0)
  { // mat_rat to char so that it can be sent
    mpq_t ****Right = (mpq_t ****)In;
    char **Left = (char **)Out;

    // initialize
    *length = 0;
    *Left = NULL;

    // setup Left
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
        for (k = 0; k < 2; k++) // for real & imag
        { // convert mpq to string
          temp = mpq_get_str(NULL, base, (*Right)[i][j][k]);
          strLength = strlen(temp) + 1; // +1 for '\0'
          // reallocate enough space
          *Left = (char *)brealloc(*Left, (*length + strLength) * sizeof(char));
          // copy on to end of *Left
          strcpy(&(*Left)[*length], temp);
          // update the length
          *length = *length + strLength;
          // free temp
          free(temp);
          temp = NULL;
        }
  }
  else
  { // char to mpq_t so that it can be used by the workers
    mpq_t ****Left = (mpq_t ****)Out;
    char **Right = (char **)In;

    // initialize
    strLength = 0;

    // setup Left
    init_mat_rat(*Left, rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
        for (k = 0; k < 2; k++) // for real & imag
        { // setup rat
          mpq_set_str((*Left)[i][j][k], &(*Right)[strLength], base);
          mpq_canonicalize((*Left)[i][j][k]);
          // update strLength
          size = strlen(&(*Right)[strLength]);
          strLength += size + 1;
        }

    // free Right
    if (freeStr)
      free(*Right);
  }

  return;
}

void cp_square_system_eval_data_mp_int(void *SSED_out, void *SSED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int setupProg, int **orig_deg, int **new_deg, int **P, int **W, char **sqStr, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - SSED_out is _d_int, SSED_in is _d      *
*                otherwise - SSED_out is _d, SSED_in is _d_int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_in to SSED_out                   *
* setupProg says whether Prog needs setup or not
\***************************************************************/
{
  int i, j, num_cols;
  char *strA = NULL, *strB = NULL, *strB_perp = NULL;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    square_system_eval_data_mp *SSED_right = (square_system_eval_data_mp *)SSED_in;
    square_system_eval_data_mp_int *SSED_left = (square_system_eval_data_mp_int *)SSED_out;

    if (setupProg)
    { // copy Prog
      cp_prog_t_int(&SSED_left->Prog_int, SSED_right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    }
    else
    { // set to NULL
      *prog_inst = NULL;
      *prog_gp_sizes = NULL;
      *progStr = NULL;
      SSED_left->Prog_int.totalLength = 0;
    }

    // setup orig_deg
    *orig_deg = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      (*orig_deg)[i] = SSED_right->orig_degrees[i];

    // setup new_deg
    *new_deg = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      (*new_deg)[i] = SSED_right->new_degrees[i];

    // setup P
    *P = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      (*P)[i] = SSED_right->P[i];

    // setup W
    num_cols = SSED_right->size_f - SSED_right->size_r;
    *W = (int *)bmalloc(SSED_right->size_r * num_cols * sizeof(int));
    for (i = 0; i < SSED_right->size_r; i++)
      for (j = 0; j < num_cols; j++)
        (*W)[num_cols * i + j] = SSED_right->W[i][j];

    // setup strB, strB_perp, strA
    cp_mat_rat_char(&strB, &SSED_right->B_rat, &SSED_left->B_strLength, SSED_right->B->rows, SSED_right->B->cols, freeStr, inType);
    cp_mat_rat_char(&strB_perp, &SSED_right->B_perp_rat, &SSED_left->B_perp_strLength, SSED_right->B_perp->rows, SSED_right->B_perp->cols, freeStr, inType);
    cp_mat_rat_char(&strA, &SSED_right->A_rat, &SSED_left->A_strLength, SSED_right->A->rows, SSED_right->A->cols, freeStr, inType);

    // allocate and setup sqStr
    SSED_left->totalLength = SSED_left->B_strLength + SSED_left->B_perp_strLength + SSED_left->A_strLength;
    *sqStr = (char *)bmalloc(SSED_left->totalLength * sizeof(char));
    // cannot use strcpy since we have multiple NULL characters throughout the 'strings'
    for (i = SSED_left->totalLength - 1; i >= 0; i--)
    {
      if (i < SSED_left->B_strLength)
        (*sqStr)[i] = strB[i];
      else if (i < SSED_left->B_strLength + SSED_left->B_perp_strLength)
        (*sqStr)[i] = strB_perp[i - SSED_left->B_strLength];
      else
        (*sqStr)[i] = strA[i - SSED_left->B_strLength - SSED_left->B_perp_strLength];
    }

    // free strB, strB_perp & strA
    free(strB);
    free(strB_perp);
    free(strA);

    // setup the rest
    SSED_left->size_f = SSED_right->size_f;
    SSED_left->B_rows = SSED_right->B->rows;
    SSED_left->B_cols = SSED_right->B->cols;
    SSED_left->B_perp_rows = SSED_right->B_perp->rows;
    SSED_left->B_perp_cols = SSED_right->B_perp->cols;
    SSED_left->noChanges = SSED_right->noChanges;
    SSED_left->max_of_W = SSED_right->max_of_W;
    SSED_left->A_rows = SSED_right->A->rows;
    SSED_left->A_cols = SSED_right->A->cols;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->prec = SSED_right->curr_prec;
    SSED_left->totalLength = SSED_left->B_strLength + SSED_left->B_perp_strLength + SSED_left->A_strLength;
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    square_system_eval_data_mp *SSED_left = (square_system_eval_data_mp *)SSED_out;
    square_system_eval_data_mp_int *SSED_right = (square_system_eval_data_mp_int *)SSED_in;

    if (setupProg)
    { // allocate memory for Prog
      SSED_left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
 
      // copy Prog
      cp_prog_t_int(SSED_left->Prog, &SSED_right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
      // initialize evalProg
      initEvalProg(1); // MPType is 1
    }

    // setup orig_degrees
    SSED_left->orig_degrees = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      SSED_left->orig_degrees[i] = (*orig_deg)[i];
    // clear orig_deg
    free(*orig_deg);
 
    // setup new_degrees
    SSED_left->new_degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      SSED_left->new_degrees[i] = (*new_deg)[i];
    // clear new_deg
    free(*new_deg);
 
    // setup P
    SSED_left->P = (int *)bmalloc(SSED_right->size_f * sizeof(int));
    for (i = SSED_right->size_f - 1; i >= 0; i--)
      SSED_left->P[i] = (*P)[i];
    // clear P
    free(*P);
 
    // setup W
    num_cols = SSED_right->size_f - SSED_right->size_r;
    SSED_left->W = (int **)bmalloc(SSED_right->size_r * sizeof(int *));
    for (i = 0; i < SSED_right->size_r; i++)
    {
      SSED_left->W[i] = (int *)bmalloc(num_cols * sizeof(int));
      for (j = 0; j < num_cols; j++)
        SSED_left->W[i][j] = (*W)[num_cols * i + j];
    }
    // clear W
    free(*W);

    // setup B
    strB = &(*sqStr)[0];
    cp_mat_rat_char(&SSED_left->B_rat, &strB, &SSED_right->B_strLength, SSED_right->B_rows, SSED_right->B_cols, 0, inType);
    init_mat_mp2(SSED_left->B, SSED_right->B_rows, SSED_right->B_cols, SSED_right->prec);
    SSED_left->B->rows = SSED_right->B_rows;
    SSED_left->B->cols = SSED_right->B_cols;
    for (i = SSED_right->B_rows - 1; i >= 0; i--)
      for (j = SSED_right->B_cols - 1; j >= 0; j--)
      {
        mpf_set_q(SSED_left->B->entry[i][j].r, SSED_left->B_rat[i][j][0]);
        mpf_set_q(SSED_left->B->entry[i][j].i, SSED_left->B_rat[i][j][1]);
      }
    strB = NULL;

    // setup B_perp
    strB_perp = &(*sqStr)[SSED_right->B_strLength];
    cp_mat_rat_char(&SSED_left->B_perp_rat, &strB_perp, &SSED_right->B_perp_strLength, SSED_right->B_perp_rows, SSED_right->B_perp_cols, 0, inType);
    init_mat_mp2(SSED_left->B_perp, SSED_right->B_perp_rows, SSED_right->B_perp_cols, SSED_right->prec);
    SSED_left->B_perp->rows = SSED_right->B_perp_rows;
    SSED_left->B_perp->cols = SSED_right->B_perp_cols;
    for (i = SSED_right->B_perp_rows - 1; i >= 0; i--)
      for (j = SSED_right->B_perp_cols - 1; j >= 0; j--)
      {
        mpf_set_q(SSED_left->B_perp->entry[i][j].r, SSED_left->B_perp_rat[i][j][0]);
        mpf_set_q(SSED_left->B_perp->entry[i][j].i, SSED_left->B_perp_rat[i][j][1]);
      }
    strB_perp = NULL;

    // setup A
    strA = &(*sqStr)[SSED_right->B_strLength + SSED_right->B_perp_strLength];
    cp_mat_rat_char(&SSED_left->A_rat, &strA, &SSED_right->A_strLength, SSED_right->A_rows, SSED_right->A_cols, 0, inType);
    init_mat_mp2(SSED_left->A, SSED_right->A_rows, SSED_right->A_cols, SSED_right->prec);
    SSED_left->A->rows = SSED_right->A_rows;
    SSED_left->A->cols = SSED_right->A_cols;
    for (i = SSED_right->A_rows - 1; i >= 0; i--)
      for (j = SSED_right->A_cols - 1; j >= 0; j--)
      {
        mpf_set_q(SSED_left->A->entry[i][j].r, SSED_left->A_rat[i][j][0]);
        mpf_set_q(SSED_left->A->entry[i][j].i, SSED_left->A_rat[i][j][1]);
      }
    strA = NULL;

    // clear sqStr
    if (freeStr)
      free(*sqStr);

    // setup the rest
    SSED_left->size_f = SSED_right->size_f;
    SSED_left->noChanges = SSED_right->noChanges;
    SSED_left->max_of_W = SSED_right->max_of_W;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->curr_prec = SSED_right->prec;
  }
  
  return;
}

void cp_patch_mp_int(void *PED_out, void *PED_in, char **patchStr, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - PED_out is _mp_int, PED_in is _mp      *
*                otherwise - PED_out is _mp, PED_in is _mp_int  *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PED_in to PED_out                     *
\***************************************************************/
{
  int i, j;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    patch_eval_data_mp *PED_right = (patch_eval_data_mp *)PED_in;
    patch_eval_data_mp_int *PED_left = (patch_eval_data_mp_int *)PED_out;

    // setup sqStr
    cp_mat_rat_char(patchStr, &PED_right->patchCoeff_rat, &PED_left->totalLength, PED_right->patchCoeff->rows, PED_right->patchCoeff->cols, freeStr, inType);
    
    // setup the rest
    PED_left->num_patches = PED_right->num_patches;
    PED_left->patchCoeff_rows = PED_right->patchCoeff->rows;
    PED_left->patchCoeff_cols = PED_right->patchCoeff->cols;
    PED_left->prec = PED_right->curr_prec;
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    patch_eval_data_mp *PED_left = (patch_eval_data_mp *)PED_out;
    patch_eval_data_mp_int *PED_right = (patch_eval_data_mp_int *)PED_in; 

    // setup patchCoeff
    cp_mat_rat_char(&PED_left->patchCoeff_rat, patchStr, &PED_right->totalLength, PED_right->patchCoeff_rows, PED_right->patchCoeff_cols, freeStr, inType);

    init_mat_mp2(PED_left->patchCoeff, PED_right->patchCoeff_rows, PED_right->patchCoeff_cols, PED_right->prec);
    PED_left->patchCoeff->rows = PED_right->patchCoeff_rows;
    PED_left->patchCoeff->cols = PED_right->patchCoeff_cols;
    for (i = PED_right->patchCoeff_rows - 1; i >= 0; i--)
      for (j = PED_right->patchCoeff_cols - 1; j >= 0; j--)
      {
        mpf_set_q(PED_left->patchCoeff->entry[i][j].r, PED_left->patchCoeff_rat[i][j][0]);
        mpf_set_q(PED_left->patchCoeff->entry[i][j].i, PED_left->patchCoeff_rat[i][j][1]);
      }

    // setup the rest
    PED_left->num_patches = PED_right->num_patches;
    PED_left->curr_prec = PED_right->prec;
  }

  return;
}

void cp_start_system_eval_data_mp_int(void *SSED_out, void *SSED_in, int **degrees, char **startStr, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - SSED_out is _mp_int, SSED_in is _mp    *
*                otherwise - SSED_out is _mp, SSED_in is _mp_int*
* RETURN VALUES:                                                *
* NOTES: stores a copy of SSED_in to SSED_out                   *
\***************************************************************/
{
  int i, j, totalDeg = 0, strLength = 0, base = 10;
  char *temp = NULL;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    start_system_eval_data_mp *SSED_right = (start_system_eval_data_mp *)SSED_in;
    start_system_eval_data_mp_int *SSED_left = (start_system_eval_data_mp_int *)SSED_out;

    // setup degrees
    *degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      totalDeg += (*degrees)[i] = SSED_right->degrees[i];

    // setup the coeff in startStr - if needed
    if (SSED_right->startSystemType == 1)
    {
      cp_mat_rat_char(startStr, &SSED_right->coeff_rat, &SSED_left->coeff_strLength, totalDeg, SSED_right->coeff_cols, freeStr, inType);
    }
    else
    {
      *startStr = NULL;
      SSED_left->coeff_strLength = 0;
    }

    // initialize the total length
    SSED_left->totalLength = SSED_left->coeff_strLength;

    // setup gamma
    for (i = 0; i < 2; i++) // for real & imag
    {
      temp = mpq_get_str(NULL, base, SSED_right->gamma_rat[i]);
      strLength = strlen(temp) + 1; // +1 for '\-'
      // reallocate enough space
      *startStr = (char *)brealloc(*startStr, (SSED_left->totalLength + strLength) * sizeof(char));
      // copy on to end of *Left
      strcpy(&(*startStr)[SSED_left->totalLength], temp);
      // update the length
      SSED_left->totalLength += strLength;
      // free temp
      free(temp);
      temp = NULL;
    }

    // setup the rest
    SSED_left->startSystemType = SSED_right->startSystemType;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->max_degree = SSED_right->max_degree;
    SSED_left->coeff_cols = SSED_right->coeff_cols;
    SSED_left->prec = SSED_right->curr_prec;
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    start_system_eval_data_mp *SSED_left = (start_system_eval_data_mp *)SSED_out;
    start_system_eval_data_mp_int *SSED_right = (start_system_eval_data_mp_int *)SSED_in;

    // setup degrees
    SSED_left->degrees = (int *)bmalloc(SSED_right->size_r * sizeof(int));
    for (i = SSED_right->size_r - 1; i >= 0; i--)
      totalDeg += SSED_left->degrees[i] = (*degrees)[i];
    // clear degrees
    free(*degrees);

    // setup the coeff in startStr - if needed
    if (SSED_right->startSystemType == 1)
    { // setup coeff_rat
      cp_mat_rat_char(&SSED_left->coeff_rat, startStr, &SSED_right->coeff_strLength, totalDeg, SSED_right->coeff_cols, 0, inType);
      // setup coeff
      SSED_left->coeff = (comp_mp **)bmalloc(totalDeg * sizeof(comp_mp *));
      for (i = 0; i < totalDeg; i++)
      {
        SSED_left->coeff[i] = (comp_mp *)bmalloc(SSED_right->coeff_cols * sizeof(comp_mp));
        for (j = 0; j < SSED_right->coeff_cols; j++)
        {
          init_mp2(SSED_left->coeff[i][j], SSED_right->prec);
          mpf_set_q(SSED_left->coeff[i][j]->r, SSED_left->coeff_rat[i][j][0]);
          mpf_set_q(SSED_left->coeff[i][j]->i, SSED_left->coeff_rat[i][j][1]);
        }
      }
    }

    // initialize where to start gamma
    strLength = SSED_right->coeff_strLength;

    // setup gamma_rat
    SSED_left->gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
    init_rat(SSED_left->gamma_rat);
    for (i = 0; i < 2; i++) // for real & imag
    { // setup rat
      mpq_set_str(SSED_left->gamma_rat[i], &(*startStr)[strLength], base);
      mpq_canonicalize(SSED_left->gamma_rat[i]);
      // update strLength
      totalDeg = strlen(&(*startStr)[strLength]);
      strLength += totalDeg + 1;
    }
    // setup gamma
    init_mp2(SSED_left->gamma, SSED_right->prec);
    mpf_set_q(SSED_left->gamma->r, SSED_left->gamma_rat[0]);
    mpf_set_q(SSED_left->gamma->i, SSED_left->gamma_rat[1]);
    // free startStr
    if (freeStr)
      free(*startStr);

    // setup the rest
    SSED_left->startSystemType = SSED_right->startSystemType;
    SSED_left->size_r = SSED_right->size_r;
    SSED_left->max_degree = SSED_right->max_degree;
    SSED_left->coeff_cols = SSED_right->coeff_cols;
    SSED_left->curr_prec = SSED_right->prec;
  }

  return;
}

void cp_basic_eval_data_mp_int(void *BED_out, void *BED_in, int **prog_inst, int **prog_gp_sizes, char **progStr, int setupProg, int **orig_deg, int **new_deg, int **P, int **W, char **sqStr, char **patchStr, int **startDeg, char **startStr, int **ppd_type, int **ppd_size, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - BED_out is _d_int, BED_in is _d        *
*                otherwise - BED_out is _d, BED_in is _d_int    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of BED_in to BED_out                     *
\***************************************************************/
{
  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    basic_eval_data_mp *BED_right = (basic_eval_data_mp *)BED_in;
    basic_eval_data_mp_int *BED_left = (basic_eval_data_mp_int *)BED_out;

    // copy square_system
    cp_square_system_eval_data_mp_int(&BED_left->squareSystem_int, &BED_right->squareSystem, prog_inst, prog_gp_sizes, progStr, setupProg, orig_deg, new_deg, P, W, sqStr, freeStr, inType);
    // copy patch
    cp_patch_mp_int(&BED_left->patch_int, &BED_right->patch, patchStr, freeStr, inType);
    // copy start_system
    cp_start_system_eval_data_mp_int(&BED_left->startSystem_int, &BED_right->startSystem, startDeg, startStr, freeStr, inType);
    // copy preproc_data
    cp_preproc_data_int(&BED_left->preProcData_int, &BED_right->preProcData, ppd_type, ppd_size, inType);
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    basic_eval_data_mp *BED_left = (basic_eval_data_mp *)BED_out;
    basic_eval_data_mp_int *BED_right = (basic_eval_data_mp_int *)BED_in;

    // copy square system
    cp_square_system_eval_data_mp_int(&BED_left->squareSystem, &BED_right->squareSystem_int, prog_inst, prog_gp_sizes, progStr, setupProg, orig_deg, new_deg, P, W, sqStr, freeStr, inType);
    // copy start system
    cp_start_system_eval_data_mp_int(&BED_left->startSystem, &BED_right->startSystem_int, startDeg, startStr, freeStr, inType);
    // copy preproc_data
    cp_preproc_data_int(&BED_left->preProcData, &BED_right->preProcData_int, ppd_type, ppd_size, inType);
    // copy patch
    cp_patch_mp_int(&BED_left->patch, &BED_right->patch_int, patchStr, freeStr, inType);
  }

  return;
}

void cp_mpf_int(void *Out, void *In, char **str, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int,In is _t                 *
*                otherwise - Out is _t, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    mpf_t *Right = (mpf_t *)In;
    mpf_int *Left = (mpf_int *)Out;

    // setup str & totalLength
    *str = mpf_to_str(*Right, base);
    Left->totalLength = strlen(*str) + 1; // +1 for '\0'

    // setup prec
    Left->prec = mpf_get_prec(*Right);
  }
  else
  { // _t_int to _t so that it can be used by the workers
    mpf_t *Left = (mpf_t *)Out;
    mpf_int *Right = (mpf_int *)In;

    // initialize
    mpf_init2(*Left, Right->prec);

    // setup Left
    mpf_set_str(*Left, *str, base);

    // free str
    if (freeStr)
      free(*str);
  }

  return;
}

void cp_comp_mp_int(void *Out, void *In, char **str, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _mp_int,In is _mp               *
*                otherwise - Out is _mp, In is _mp_int          *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, size, base = 10;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    comp_mp *Right = (comp_mp *)In;
    comp_mp_int *Left = (comp_mp_int *)Out;

    int sizeR, sizeI;
    char *strR = NULL, *strI = NULL;

    // real part
    strR = mpf_to_str((*Right)->r, base);
    sizeR = strlen(strR) + 1; // +1 for '\0'
    // imag part
    strI = mpf_to_str((*Right)->i, base);
    sizeI = strlen(strI) + 1; // +1 for '\0' 

    // find totalLength and setup str
    Left->totalLength = sizeR + sizeI;
    *str = (char *)bmalloc(Left->totalLength * sizeof(char));
    // copy to str - cannot use strcpy
    for (i = 0; i < Left->totalLength; i++)
    {
      if (i < sizeR)
        (*str)[i] = strR[i];
      else
        (*str)[i] = strI[i - sizeR];
    }

    // setup prec
    Left->prec = mpf_get_prec((*Right)->r);

    // free strR, strI
    free(strR);
    free(strI);
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    comp_mp *Left = (comp_mp *)Out;
    comp_mp_int *Right = (comp_mp_int *)In; 

    // initialize
    init_mp2((*Left), Right->prec);

    // setup real part
    mpf_set_str((*Left)->r, *str, base);

    // find where the imag part starts
    size = strlen(*str) + 1;

    // setup imag part
    mpf_set_str((*Left)->i, &(*str)[size], base);

    // free str
    if (freeStr)
      free(*str);
  }
  
  return;
}

void cp_comp_rat_int(void *Out, void *In, char **str, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _rat_int,In is _rat             *
*                otherwise - Out is _rat, In is _rat_int        *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, base = 10;

  if (inType == 0)
  { // _rat to _rat_int so that it can be sent
    mpq_t *Right = (mpq_t *)In;
    comp_rat_int *Left = (comp_rat_int *)Out;

    int sizeR, sizeI;
    char *strR = NULL, *strI = NULL;

    // real part
    strR = mpq_get_str(NULL, base, Right[0]);
    sizeR = Left->length[0] = strlen(strR) + 1; // +1 for '\0'
    // imag part
    strI = mpq_get_str(NULL, base, Right[1]);
    sizeI = Left->length[1] = strlen(strI) + 1; // +1 for '\0'

    // setup str
    *str = (char *)bmalloc((sizeR + sizeI) * sizeof(char));
    // copy to str - cannot use strcpy
    for (i = 0; i < sizeR; i++)
      (*str)[i] = strR[i];
    for (i = 0; i < sizeI; i++)
      (*str)[i + sizeR] = strI[i];

    // free strR, strI
    free(strR);
    free(strI);
  }
  else
  { // _rat_int to _rat so that it can be used by the workers
    mpq_t *Left = (mpq_t *)Out;
    comp_rat_int *Right = (comp_rat_int *)In;

    // initialize
    init_rat(Left);

    // setup real part
    mpq_set_str(Left[0], &(*str)[0], base);
    mpq_canonicalize(Left[0]);
    // setup imag part
    mpq_set_str(Left[1], &(*str)[Right->length[0]], base);
    mpq_canonicalize(Left[1]);

    // free str
    if (freeStr)
      free(*str);
  }

  return;
}

void cp_point_d_int(void *Out, void *In, comp_d **coeff, int freeCoeff, int initPoint, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _d_int,In is _d                 *
*                otherwise - Out is _d, In is _d_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, size;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    _point_d *Right = (_point_d *)In;
    point_d_int *Left = (point_d_int *)Out;

    // initialize
    size = Left->size = (*Right).size;

    *coeff = (comp_d *)bmalloc(size * sizeof(comp_d));
    // setup coeff
    for (i = 0; i < size; i++)
    {
      set_d((*coeff)[i], &(*Right).coord[i]);
    }
  }
  else
  { // _d_int to _d so that it can be sent
    point_d_int *Right = (point_d_int *)In;
    _point_d *Left = (_point_d *)Out;

    // setup size
    size = Right->size;

    // initialize
    if (initPoint)
    {
      init_point_d(&(*Left), size);
    }
    else
    {
      change_size_point_d(&(*Left), size);
    }
    (*Left).size = size;

    // setup each coord
    for (i = 0; i < size; i++)
    {
      set_d(&(*Left).coord[i], (*coeff)[i]);
    }

    // clear coeff
    if (freeCoeff)
      free(*coeff);
  }

  return;
}

void cp_point_mp_int(void *Out, void *In, char **str, int freeStr, int initPoint, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _mp_int,In is _mp               *
*                otherwise - Out is _mp, In is _mp_int          *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, size, strLoc, base = 10;
  
  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    _point_mp *Right = (_point_mp *)In;
    point_mp_int *Left = (point_mp_int *)Out;

    int sizeR, sizeI;
    char *strR = NULL, *strI = NULL;

    // initialize
    size = Left->size = (*Right).size;
    strLoc = 0;
    *str = NULL;

    // setup str
    for (i = 0; i < size; i++)
    { // real part
      strR = mpf_to_str((*Right).coord[i].r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str((*Right).coord[i].i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0' 

      // update str
      *str = (char *)brealloc(*str, (strLoc + sizeR + sizeI) * sizeof(char));
      for (j = sizeR + sizeI - 1; j >= 0; j--)
      {
        if (j < sizeR)
          (*str)[j + strLoc] = strR[j];
        else
          (*str)[j + strLoc] = strI[j - sizeR];
      }
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);
    }
    // setup totalLength
    Left->totalLength = strLoc;

    // setup prec
    Left->prec = (*Right).curr_prec;
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    _point_mp *Left = (_point_mp *)Out;
    point_mp_int *Right = (point_mp_int *)In;

    // initialize
    if (initPoint) 
    {
      init_point_mp2(&(*Left), Right->size, Right->prec);
    }
    else
    {
      change_prec_point_mp(&(*Left), Right->prec);
      change_size_point_mp(&(*Left), Right->size);
    }
    strLoc = 0;

    // setup size
    (*Left).size = Right->size;

    // setup each coord
    for (i = 0; i < Right->size; i++)
    { // setup real part
      mpf_set_str((*Left).coord[i].r, &(*str)[strLoc], base);

      // find where the imag part starts
      size = strlen(&(*str)[strLoc]) + 1;
      strLoc += size;

      // setup imag part
      mpf_set_str((*Left).coord[i].i, &(*str)[strLoc], base);

      // find where the next real part starts
      size = strlen(&(*str)[strLoc]) + 1;
      strLoc += size;
    }

    // free str
    if (freeStr)
      free(*str);
  }

  return;
}

void cp_point_data_d_int(void *PD_out, void *PD_in, comp_d **coeff, int freeCoeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - PD_out is _d_int, PD_in is _d          *
*                otherwise - PD_out is _d, PD_in is _d_int      *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PD_in to PD_out                       *
\***************************************************************/
{
  if (inType == 0)
  { // _d to _d_int so that it can be sent
    point_data_d *PD_right = (point_data_d *)PD_in;
    point_data_d_int *PD_left = (point_data_d_int *)PD_out;

    // setup point
    cp_point_d_int(&PD_left->point_int, &PD_right->point, coeff, freeCoeff, 0, inType);

    // setup time
    set_d(PD_left->time, PD_right->time);

    // setup cycle_num
    PD_left->cycle_num = PD_right->cycle_num;
  }
  else
  { // _d_int to _d so that it can be used by the workers
    point_data_d *PD_left = (point_data_d *)PD_out;
    point_data_d_int *PD_right = (point_data_d_int *)PD_in;

    // setup point
    cp_point_d_int(&PD_left->point, &PD_right->point_int, coeff, freeCoeff, 0, inType);
    // setup time
    set_d(PD_left->time, PD_right->time);
    // setup cycle_num
    PD_left->cycle_num = PD_right->cycle_num;

    // free coeff
    if (freeCoeff)
      free(*coeff);
  }

  return;
}

void cp_point_data_mp_int(void *PD_out, void *PD_in, char **str, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - PD_out is _mp_int, PD_in is _mp        *
*                otherwise - PD_out is _mp, PD_in is _mp_int    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of PD_in to PD_out                       *
\***************************************************************/
{
  int i;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    point_data_mp *PD_right = (point_data_mp *)PD_in;
    point_data_mp_int *PD_left = (point_data_mp_int *)PD_out;

    char *strPoint = NULL, *strTime = NULL;

    // setup point
    cp_point_mp_int(&PD_left->point_int, PD_right->point, &strPoint, 0, 0, inType);

    // setup time
    cp_comp_mp_int(&PD_left->time_int, &PD_right->time, &strTime, 0, inType);

    // setup str & totalLength
    PD_left->totalLength = PD_left->point_int.totalLength + PD_left->time_int.totalLength;
    *str = (char *)bmalloc(PD_left->totalLength * sizeof(char));
    for (i = PD_left->totalLength - 1; i >= 0; i--)
    {
      if (i < PD_left->point_int.totalLength)
        (*str)[i] = strPoint[i];
      else
        (*str)[i] = strTime[i - PD_left->point_int.totalLength];
    }
    
    // setup cycle_num
    PD_left->cycle_num = PD_right->cycle_num;

    // clear memory
    free(strPoint);
    free(strTime);
  }
  else
  { // _mp_int to _mp so that it can be used by the workers
    point_data_mp *PD_left = (point_data_mp *)PD_out;
    point_data_mp_int *PD_right = (point_data_mp_int *)PD_in;

    // find where time starts
    char *strTime = &(*str)[PD_right->point_int.totalLength];

    // setup point 
    cp_point_mp_int(PD_left->point, &PD_right->point_int, str, 0, 0, inType);
 
    // setup time 
    cp_comp_mp_int(&PD_left->time, &PD_right->time_int, &strTime, 0, inType);

    // setup cycle_num
    PD_left->cycle_num = PD_right->cycle_num;

    // free str
    if (freeStr)
      free(*str);
  }

  return;
}

void cp_mat_d_int(void *Out, void *In, comp_d **coeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _d_int,In is _d                 *
*                otherwise - Out is _d, In is _d_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, r, c, count;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    _mat_d *Right = (_mat_d *)In;
    mat_d_int *Left = (mat_d_int *)Out;

    // initialize
    r = Left->rows = (*Right).rows;
    c = Left->cols = (*Right).cols;

    *coeff = (comp_d *)bmalloc(r * c * sizeof(comp_d));
    // setup coeff
    count = 0;
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      {
        set_d((*coeff)[count], &(*Right).entry[i][j]);
        count++;
      }
  }
  else
  { // _d_int to _d so that it can be sent
    mat_d_int *Right = (mat_d_int *)In;
    _mat_d *Left = (_mat_d *)Out;

    // setup r & c
    r = Right->rows;
    c = Right->cols;

    // initialize
    init_mat_d(&(*Left), r, c);
    (*Left).rows = r;
    (*Left).cols = c;

    // setup each entry
    count = 0;
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      {
        set_d(&(*Left).entry[i][j], (*coeff)[count]);
        count++;
      }

    // clear coeff
    free(*coeff);
  }

  return;
}

void cp_endgame_data_t_int(void *EG_out, void *EG_in, char **egStr, comp_d **coeff, int freeStr, int freeCoeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - EG_out is _t_int, EG_in is _t          *
*                otherwise - EG_out is _t, EG_in is _t_int      *
* RETURN VALUES:                                                *
* NOTES: stores a copy of EG_in to EG_out                       *
\***************************************************************/
{
  int i;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    endgame_data_t *EG_right = (endgame_data_t *)EG_in;
    endgame_data_t_int *EG_left = (endgame_data_t_int *)EG_out;

    if (EG_right->prec == 52)
    { // copy over double precision values
      comp_d *PD = NULL, *Last = NULL;
      cp_point_data_d_int(&EG_left->PD_d_int, &EG_right->PD_d, &PD, freeCoeff, inType);

      EG_left->function_residual_d = EG_right->function_residual_d;
      EG_left->latest_newton_residual_d = EG_right->latest_newton_residual_d;
      EG_left->t_val_at_latest_sample_point_d = EG_right->t_val_at_latest_sample_point_d;
      EG_left->error_at_latest_sample_point_d = EG_right->error_at_latest_sample_point_d;

      if (EG_right->last_approx_prec == 52)
      { // setup Last
        cp_point_d_int(&EG_left->last_approx_d_int, &EG_right->last_approx_d, &Last, freeCoeff, 0, inType);

        // combine PD & Last
        *coeff = (comp_d *)bmalloc((EG_right->PD_d.point->size + EG_right->last_approx_d->size) * sizeof(comp_d));
        for (i = 0; i < EG_right->PD_d.point->size + EG_right->last_approx_d->size; i++)
        {
          if (i < EG_right->PD_d.point->size)
          {
            set_d((*coeff)[i], PD[i]);
          }
          else
          {
            set_d((*coeff)[i], Last[i - EG_right->PD_d.point->size]);
          }
        }

        free(PD);
        free(Last);

        *egStr = NULL;
        EG_left->totalLength = 0;
      }
      else
      { // setup str & totalLength for last_approx_mp
        cp_point_mp_int(&EG_left->last_approx_mp_int, &EG_right->last_approx_mp, egStr, 0, 0, inType);
        EG_left->totalLength = EG_left->last_approx_mp_int.totalLength;

        // setup coeff
        *coeff = (comp_d *)bmalloc(EG_right->PD_d.point->size * sizeof(comp_d));
        for (i = 0; i < EG_right->PD_d.point->size; i++)
        {
          set_d((*coeff)[i], PD[i]);
        }

        free(PD);
      }
    }
    else
    { // need to setup for MP
      char *strPD = NULL, *strFunc = NULL, *strNewton = NULL, *strTval = NULL, *strError = NULL, *strLast = NULL;
      cp_point_data_mp_int(&EG_left->PD_mp_int, &EG_right->PD_mp, &strPD, 0, inType);
      cp_mpf_int(&EG_left->function_residual_mp_int, &EG_right->function_residual_mp, &strFunc, 0, inType);
      cp_mpf_int(&EG_left->latest_newton_residual_mp_int, &EG_right->latest_newton_residual_mp, &strNewton, 0, inType);
      cp_mpf_int(&EG_left->t_val_at_latest_sample_point_mp_int, &EG_right->t_val_at_latest_sample_point_mp, &strTval, 0, inType);
      cp_mpf_int(&EG_left->error_at_latest_sample_point_mp_int, &EG_right->error_at_latest_sample_point_mp, &strError, 0, inType);

      if (EG_right->last_approx_prec == 52)
      { // setup coeff
        cp_point_d_int(&EG_left->last_approx_d_int, &EG_right->last_approx_d, coeff, freeCoeff, 0, inType);
        EG_left->last_approx_mp_int.totalLength = 0;
      }
      else
      { // setup strLast
        cp_point_mp_int(&EG_left->last_approx_mp_int, &EG_right->last_approx_mp, &strLast, 0, 0, inType);
        *coeff = NULL;
      }

      // setup str & totalLength
      EG_left->totalLength = EG_left->PD_mp_int.totalLength + EG_left->last_approx_mp_int.totalLength + EG_left->function_residual_mp_int.totalLength + EG_left->latest_newton_residual_mp_int.totalLength + EG_left->t_val_at_latest_sample_point_mp_int.totalLength + EG_left->error_at_latest_sample_point_mp_int.totalLength;
      *egStr = (char *)bmalloc(EG_left->totalLength * sizeof(char));
      for (i = EG_left->totalLength - 1; i >= 0; i--)
      {
        if (i < EG_left->PD_mp_int.totalLength)
          (*egStr)[i] = strPD[i];
        else if (i < EG_left->PD_mp_int.totalLength + EG_left->last_approx_mp_int.totalLength)
          (*egStr)[i] = strLast[i - EG_left->PD_mp_int.totalLength];
        else if (i < EG_left->PD_mp_int.totalLength + EG_left->last_approx_mp_int.totalLength + EG_left->function_residual_mp_int.totalLength)
          (*egStr)[i] = strFunc[i - EG_left->PD_mp_int.totalLength - EG_left->last_approx_mp_int.totalLength];
        else if (i < EG_left->PD_mp_int.totalLength + EG_left->last_approx_mp_int.totalLength + EG_left->function_residual_mp_int.totalLength + EG_left->latest_newton_residual_mp_int.totalLength)
          (*egStr)[i] = strNewton[i - EG_left->PD_mp_int.totalLength - EG_left->last_approx_mp_int.totalLength - EG_left->function_residual_mp_int.totalLength];
        else if (i < EG_left->PD_mp_int.totalLength + EG_left->last_approx_mp_int.totalLength + EG_left->function_residual_mp_int.totalLength + EG_left->latest_newton_residual_mp_int.totalLength + EG_left->t_val_at_latest_sample_point_mp_int.totalLength) 
          (*egStr)[i] = strTval[i - EG_left->PD_mp_int.totalLength - EG_left->last_approx_mp_int.totalLength - EG_left->function_residual_mp_int.totalLength - EG_left->latest_newton_residual_mp_int.totalLength];
        else
          (*egStr)[i] = strError[i - EG_left->PD_mp_int.totalLength - EG_left->last_approx_mp_int.totalLength - EG_left->function_residual_mp_int.totalLength - EG_left->latest_newton_residual_mp_int.totalLength - EG_left->t_val_at_latest_sample_point_mp_int.totalLength];
      }

      // free strPD, strFunc, strNewton, strTval, strError
      free(strPD);
      free(strLast);
      free(strFunc);
      free(strNewton);
      free(strTval);
      free(strError);
    }

    // setup the rest
    EG_left->prec = EG_right->prec;
    EG_left->last_approx_prec = EG_right->last_approx_prec;
    EG_left->retVal = EG_right->retVal;
    EG_left->pathNum = EG_right->pathNum;
    EG_left->codim = EG_right->codim;
    EG_left->first_increase = EG_right->first_increase;
    EG_left->condition_number = EG_right->condition_number;
  }
  else
  { // _t_int to _t so that it can be used by the workers
    endgame_data_t *EG_left = (endgame_data_t *)EG_out;
    endgame_data_t_int *EG_right = (endgame_data_t_int *)EG_in;

    if (EG_right->prec == 52)
    { // copy over double precision values
      comp_d *Last = NULL;

      // setup PD
      cp_point_data_d_int(&EG_left->PD_d, &EG_right->PD_d_int, coeff, 0, inType);

      EG_left->function_residual_d = EG_right->function_residual_d;
      EG_left->latest_newton_residual_d = EG_right->latest_newton_residual_d;
      EG_left->t_val_at_latest_sample_point_d = EG_right->t_val_at_latest_sample_point_d;
      EG_left->error_at_latest_sample_point_d = EG_right->error_at_latest_sample_point_d;

      if (EG_right->last_approx_prec == 52)
      { // setup Last
        Last = &(*coeff)[EG_right->PD_d_int.point_int.size];

        // setup last_approx 
        cp_point_d_int(&EG_left->last_approx_d, &EG_right->last_approx_d_int, &Last, 0, 0, inType);
      }
      else
      { // setup form str -- this will clear egStr if needed
        cp_point_mp_int(&EG_left->last_approx_mp, &EG_right->last_approx_mp_int, egStr, freeStr, 0, inType);
      }

      // clear coeff
      if (freeCoeff)
        free(*coeff);
    }
    else
    { // setup pointers to the values
      char *strLast = &(*egStr)[EG_right->PD_mp_int.totalLength];
      char *strFunc = &(*egStr)[EG_right->PD_mp_int.totalLength + EG_right->last_approx_mp_int.totalLength];
      char *strNewton = &(*egStr)[EG_right->PD_mp_int.totalLength + EG_right->last_approx_mp_int.totalLength + EG_right->function_residual_mp_int.totalLength];
      char *strTval = &(*egStr)[EG_right->PD_mp_int.totalLength + EG_right->last_approx_mp_int.totalLength + EG_right->function_residual_mp_int.totalLength + EG_right->latest_newton_residual_mp_int.totalLength];
      char *strError = &(*egStr)[EG_right->PD_mp_int.totalLength + EG_right->last_approx_mp_int.totalLength + EG_right->function_residual_mp_int.totalLength + EG_right->latest_newton_residual_mp_int.totalLength + EG_right->error_at_latest_sample_point_mp_int.totalLength];

      // setup each of the MP values
      cp_point_data_mp_int(&EG_left->PD_mp, &EG_right->PD_mp_int, egStr, 0, inType);
      cp_mpf_int(&EG_left->function_residual_mp, &EG_right->function_residual_mp_int, &strFunc, 0, inType);
      cp_mpf_int(&EG_left->latest_newton_residual_mp, &EG_right->latest_newton_residual_mp_int, &strNewton, 0, inType);
      cp_mpf_int(&EG_left->t_val_at_latest_sample_point_mp, &EG_right->t_val_at_latest_sample_point_mp_int, &strTval, 0, inType);
      cp_mpf_int(&EG_left->error_at_latest_sample_point_mp, &EG_right->error_at_latest_sample_point_mp_int, &strError, 0, inType);

      if (EG_right->last_approx_prec == 52)
      { // setup from coeff -- this will clear coeff if needed
        cp_point_d_int(&EG_left->last_approx_d, &EG_right->last_approx_d_int, coeff, freeCoeff, 0, inType);
      }
      else
      { // setup from strLast
        cp_point_mp_int(&EG_left->last_approx_mp, &EG_right->last_approx_mp_int, &strLast, 0, 0, inType);
      }
  
      // set pointer to NULL
      strLast = strFunc = strNewton = strTval = strError = NULL;

      // free egStr
      if (freeStr)
        free(*egStr);
    }

    // setup the rest
    EG_left->prec = EG_right->prec;
    EG_left->last_approx_prec = EG_right->last_approx_prec;
    EG_left->retVal = EG_right->retVal;
    EG_left->pathNum = EG_right->pathNum;
    EG_left->codim = EG_right->codim;
    EG_left->first_increase = EG_right->first_increase;
    EG_left->condition_number = EG_right->condition_number;
  }

  return;
}

void cp_eqData_int(void *EqD_out, void *EqD_in, int MPType, char **eqdStr, int freeStr, comp_d **coeff_d, int **degrees, int **instCount, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - EqD_out is _t_int, EqD_in is _t        *
*                otherwise - EqD_out is _t, EqD_in is _t_int    *
* RETURN VALUES:                                                *
* NOTES: stores a copy of EqD_in to EqD_out                     *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    eqData_t *EqD_right = (eqData_t *)EqD_in;
    eqData_t_int *EqD_left = (eqData_t_int *)EqD_out;

    if (MPType == 0)
    { // count the total number of coeff_d
      EqD_left->num_coeff = EqD_right->num_funcs * EqD_right->num_vars;

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(EqD_left->num_coeff * sizeof(comp_d));
      count = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
        for (j = 0; j < EqD_right->num_vars; j++)
        {
          set_d((*coeff_d)[count], EqD_right->coeff_d[i][j]);
          count++;
        }

      // setup gamma
      set_d(EqD_left->gamma_d, EqD_right->gamma_d);

      // setup totalLength
      EqD_left->totalLength = 0;
    }
    else if (MPType == 1)
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      EqD_left->num_coeff = 0; // no coeff setup in coeff_d

      // setup eqdStr
      *eqdStr = NULL;
      strLoc = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
        for (j = 0; j < EqD_right->num_vars; j++)
        { // real part
          strR = mpf_to_str(EqD_right->coeff_mp[i][j]->r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(EqD_right->coeff_mp[i][j]->i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update eqdStr
          *eqdStr = (char *)brealloc(*eqdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*eqdStr)[strLoc], strR);
          strcpy(&(*eqdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup gamma
      // real part
      strR = mpf_to_str(EqD_right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(EqD_right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update eqdStr
      *eqdStr = (char *)brealloc(*eqdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*eqdStr)[strLoc], strR);
      strcpy(&(*eqdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup totalLength
      EqD_left->totalLength = strLoc;
    }
    else // MPType == 2
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      EqD_left->num_coeff = 0; // no coeff setup in coeff_d

      // setup eqdStr
      *eqdStr = NULL;
      strLoc = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
        for (j = 0; j < EqD_right->num_vars; j++)
        { // real part
          strR = mpq_get_str(NULL, base, EqD_right->coeff_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, EqD_right->coeff_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update eqdStr
          *eqdStr = (char *)brealloc(*eqdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*eqdStr)[strLoc], strR);
          strcpy(&(*eqdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup gamma
      // real part
      strR = mpq_get_str(NULL, base, EqD_right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, EqD_right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update eqdStr
      *eqdStr = (char *)brealloc(*eqdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*eqdStr)[strLoc], strR);
      strcpy(&(*eqdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup totalLength
      EqD_left->totalLength = strLoc;
    }

    // setup degrees
    *degrees = (int *)bmalloc(EqD_right->num_funcs * EqD_right->num_var_gps * sizeof(int));
    for (i = 0; i < EqD_right->num_funcs; i++)
      for (j = 0; j < EqD_right->num_var_gps; j++)
        (*degrees)[i * EqD_right->num_var_gps + j] = EqD_right->degrees[i][j];

    // setup instCount, if needed
    EqD_left->noChanges = EqD_right->noChanges;
    EqD_left->numSubFuncs = EqD_right->numSubFuncs;
    if (EqD_left->noChanges)
    { // count the number of things to send
      EqD_left->numInts = 4 * (EqD_right->numSubFuncs + EqD_right->num_funcs) + EqD_right->numSubFuncs * EqD_right->num_funcs;
      *instCount = (int *)bmalloc(EqD_left->numInts * sizeof(int));
      // setup instCount
      for (i = 0; i < EqD_right->numSubFuncs; i++)
      { // copy the things related to subfunctions
        (*instCount)[i] = EqD_right->startSub[i];
        (*instCount)[i + EqD_right->numSubFuncs] = EqD_right->endSub[i];
        (*instCount)[i + 2 * EqD_right->numSubFuncs] = EqD_right->startJvsub[i];
        (*instCount)[i + 3 * EqD_right->numSubFuncs] = EqD_right->endJvsub[i];
      }
      for (i = 4 * EqD_right->numSubFuncs; i < 4 * EqD_right->numSubFuncs + EqD_right->num_funcs; i++)
      { // copy the things releted to functions
        (*instCount)[i] = EqD_right->startFunc[i - 4 * EqD_right->numSubFuncs];
        (*instCount)[i + EqD_right->num_funcs] = EqD_right->endFunc[i - 4 * EqD_right->numSubFuncs];
        (*instCount)[i + 2 * EqD_right->num_funcs] = EqD_right->startJv[i - 4 * EqD_right->numSubFuncs];
        (*instCount)[i + 3 * EqD_right->num_funcs] = EqD_right->endJv[i - 4 * EqD_right->numSubFuncs];
      }
      for (i = 0; i < EqD_right->num_funcs; i++)
        for (j = 0; j < EqD_right->numSubFuncs; j++)
          (*instCount)[4 * (EqD_right->num_funcs + EqD_right->numSubFuncs) + i * EqD_right->numSubFuncs + j] = EqD_right->subFuncsBelow[i][j];
    }
    else
    { // set to NULL
      EqD_left->numInts = 0;
      *instCount = NULL;
    }

    // copy over the rest of the values
    EqD_left->curr_precision = EqD_right->curr_precision;
    EqD_left->num_funcs = EqD_right->num_funcs;
    EqD_left->num_subsystems = EqD_right->num_subsystems;
    EqD_left->num_var_gps = EqD_right->num_var_gps;
    EqD_left->num_vars = EqD_right->num_vars;
  }
  else
  { // _t_int to _t so that it can be used by the workers
    eqData_t *EqD_left = (eqData_t *)EqD_out;
    eqData_t_int *EqD_right = (eqData_t_int *)EqD_in;

    // setup the degrees
    EqD_left->degrees = (int **)bmalloc(EqD_right->num_funcs * sizeof(int *));
    for (i = 0; i < EqD_right->num_funcs; i++)
    {
      EqD_left->degrees[i] = (int *)bmalloc(EqD_right->num_var_gps * sizeof(int));
      for (j = 0; j < EqD_right->num_var_gps; j++)
        EqD_left->degrees[i][j] = (*degrees)[i * EqD_right->num_var_gps + j];
    }
    free(*degrees);

    // setup instCount
    EqD_left->noChanges = EqD_right->noChanges;
    EqD_left->numSubFuncs = EqD_right->numSubFuncs;
    if (EqD_left->noChanges)
    { // allocate
      EqD_left->startSub = (int *)bmalloc(EqD_right->numSubFuncs * sizeof(int));
      EqD_left->endSub = (int *)bmalloc(EqD_right->numSubFuncs * sizeof(int));
      EqD_left->startFunc = (int *)bmalloc(EqD_right->num_funcs * sizeof(int));
      EqD_left->endFunc = (int *)bmalloc(EqD_right->num_funcs * sizeof(int));
      EqD_left->startJvsub = (int *)bmalloc(EqD_right->numSubFuncs * sizeof(int));
      EqD_left->endJvsub = (int *)bmalloc(EqD_right->numSubFuncs * sizeof(int));
      EqD_left->startJv = (int *)bmalloc(EqD_right->num_funcs * sizeof(int));
      EqD_left->endJv = (int *)bmalloc(EqD_right->num_funcs * sizeof(int));

      // setup
      for (i = 0; i < EqD_right->numSubFuncs; i++)
      { // copy the things related to subfunctions
        EqD_left->startSub[i] = (*instCount)[i];
        EqD_left->endSub[i] = (*instCount)[i + EqD_right->numSubFuncs];
        EqD_left->startJvsub[i] = (*instCount)[i + 2 * EqD_right->numSubFuncs];
        EqD_left->endJvsub[i] = (*instCount)[i + 3 * EqD_right->numSubFuncs];
      }
      for (i = 0; i < EqD_right->num_funcs; i++)
      { // copy the things releted to functions
        EqD_left->startFunc[i] = (*instCount)[i + 4 * EqD_right->numSubFuncs];
        EqD_left->endFunc[i] = (*instCount)[i + EqD_right->num_funcs + 4 * EqD_right->numSubFuncs];
        EqD_left->startJv[i] = (*instCount)[i + 2 * EqD_right->num_funcs + 4 * EqD_right->numSubFuncs];
        EqD_left->endJv[i] = (*instCount)[i + 3 * EqD_right->num_funcs + 4 * EqD_right->numSubFuncs];
      }
      if (EqD_right->numSubFuncs > 0)
      { // setup subFuncsBelow
        EqD_left->subFuncsBelow = (int **)bmalloc(EqD_right->num_funcs * sizeof(int *));
        for (i = 0; i < EqD_right->num_funcs; i++)
          EqD_left->subFuncsBelow[i] = (int *)bmalloc(EqD_right->numSubFuncs * sizeof(int));

        for (i = 0; i < EqD_right->num_funcs; i++)
          for (j = 0; j < EqD_right->numSubFuncs; j++)
            EqD_left->subFuncsBelow[i][j] = (*instCount)[4 * (EqD_right->num_funcs + EqD_right->numSubFuncs) + i * EqD_right->numSubFuncs + j];
      }
      else
        EqD_left->subFuncsBelow = NULL;

      // clear instCount
      free(*instCount);
    }
    else
    { // set to NULL
      EqD_left->startSub = EqD_left->endSub = EqD_left->startFunc = EqD_left->endFunc = EqD_left->startJvsub = EqD_left->endJvsub = EqD_left->startJv = EqD_left->endJv = NULL;
      EqD_left->subFuncsBelow = NULL;
    }

    if (MPType == 0)
    { // setup coeff_d
      EqD_left->coeff_d = (comp_d **)bmalloc(EqD_right->num_funcs * sizeof(comp_d *));
      count = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
      {
        EqD_left->coeff_d[i] = (comp_d *)bmalloc(EqD_right->num_vars * sizeof(comp_d));
        for (j = 0; j < EqD_right->num_vars; j++)
        {
          set_d(EqD_left->coeff_d[i][j], (*coeff_d)[count]);
          count++;
        }
      }
      free(*coeff_d);

      // setup gamma
      set_d(EqD_left->gamma_d, EqD_right->gamma_d);
    }
    else if (MPType == 1)
    { // setup coeff_mp
      EqD_left->coeff_mp = (comp_mp **)bmalloc(EqD_right->num_funcs * sizeof(comp_mp *));
      strLoc = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
      {
        EqD_left->coeff_mp[i] = (comp_mp *)bmalloc(EqD_right->num_vars * sizeof(comp_mp));
        for (j = 0; j < EqD_right->num_vars; j++)
        {
          init_mp2(EqD_left->coeff_mp[i][j], EqD_right->curr_precision);
          // setup real part
          mpf_set_str(EqD_left->coeff_mp[i][j]->r, &(*eqdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*eqdStr)[strLoc]);
          // setup imag part
          mpf_set_str(EqD_left->coeff_mp[i][j]->i, &(*eqdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*eqdStr)[strLoc]);
        }
      }
      // setup gamma
      init_mp2(EqD_left->gamma_mp, EqD_right->curr_precision);
      // real part
      mpf_set_str(EqD_left->gamma_mp->r, &(*eqdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*eqdStr)[strLoc]);
      // imag part
      mpf_set_str(EqD_left->gamma_mp->i, &(*eqdStr)[strLoc], base);

      // free str
      if (freeStr)
        free(*eqdStr);
    }
    else // MPType == 2
    { // setup coeff_rat, then coeff_d & coeff_mp
      EqD_left->coeff_d = (comp_d **)bmalloc(EqD_right->num_funcs * sizeof(comp_d *));
      EqD_left->coeff_mp = (comp_mp **)bmalloc(EqD_right->num_funcs * sizeof(comp_mp *));
      EqD_left->coeff_rat = (mpq_t ***)bmalloc(EqD_right->num_funcs * sizeof(mpq_t **));
      strLoc = 0;
      for (i = 0; i < EqD_right->num_funcs; i++)
      {
        EqD_left->coeff_d[i] = (comp_d *)bmalloc(EqD_right->num_vars * sizeof(comp_d));
        EqD_left->coeff_mp[i] = (comp_mp *)bmalloc(EqD_right->num_vars * sizeof(comp_mp));
        EqD_left->coeff_rat[i] = (mpq_t **)bmalloc(EqD_right->num_vars * sizeof(mpq_t *));
        for (j = 0; j < EqD_right->num_vars; j++)
        { // setup coeff_rat
          EqD_left->coeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          // setup real part
          mpq_init(EqD_left->coeff_rat[i][j][0]);
          mpq_set_str(EqD_left->coeff_rat[i][j][0], &(*eqdStr)[strLoc], base);
          mpq_canonicalize(EqD_left->coeff_rat[i][j][0]);
          strLoc += 1 + strlen(&(*eqdStr)[strLoc]);
          // setup imag part
          mpq_init(EqD_left->coeff_rat[i][j][1]);
          mpq_set_str(EqD_left->coeff_rat[i][j][1], &(*eqdStr)[strLoc], base);
          mpq_canonicalize(EqD_left->coeff_rat[i][j][1]);
          strLoc += 1 + strlen(&(*eqdStr)[strLoc]);

          // setup coeff_mp & coeff_d
          init_mp2(EqD_left->coeff_mp[i][j], EqD_right->curr_precision);
          mpf_set_q(EqD_left->coeff_mp[i][j]->r, EqD_left->coeff_rat[i][j][0]);
          EqD_left->coeff_d[i][j]->r = mpq_get_d(EqD_left->coeff_rat[i][j][0]);

          mpf_set_q(EqD_left->coeff_mp[i][j]->i, EqD_left->coeff_rat[i][j][1]);
          EqD_left->coeff_d[i][j]->i = mpq_get_d(EqD_left->coeff_rat[i][j][1]);
        }
      }
      // setup gamma_rat
      mpq_init(EqD_left->gamma_rat[0]);
      mpq_set_str(EqD_left->gamma_rat[0], &(*eqdStr)[strLoc], base);
      mpq_canonicalize(EqD_left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*eqdStr)[strLoc]);
      mpq_init(EqD_left->gamma_rat[1]);
      mpq_set_str(EqD_left->gamma_rat[1], &(*eqdStr)[strLoc], base);
      mpq_canonicalize(EqD_left->gamma_rat[1]);

      // setup gamma_mp & gamma_d
      init_mp2(EqD_left->gamma_mp, EqD_right->curr_precision);
      mpf_set_q(EqD_left->gamma_mp->r, EqD_left->gamma_rat[0]);
      EqD_left->gamma_d->r = mpq_get_d(EqD_left->gamma_rat[0]);

      mpf_set_q(EqD_left->gamma_mp->i, EqD_left->gamma_rat[1]);
      EqD_left->gamma_d->i = mpq_get_d(EqD_left->gamma_rat[1]);

      // free str
      if (freeStr)
        free(*eqdStr);
    }

    // copy over the rest of the values
    EqD_left->curr_precision = EqD_right->curr_precision;
    EqD_left->num_funcs = EqD_right->num_funcs;
    EqD_left->num_subsystems = EqD_right->num_subsystems;
    EqD_left->num_var_gps = EqD_right->num_var_gps;
    EqD_left->num_vars = EqD_right->num_vars;
  }

  return;
}

void cp_witnessData_int(void *Out, void *In, int stage, int MPType, char **wdStr, int freeStr, comp_d **coeff_d, int inType)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is EqD               *
*                otherwise - Out is EqD, In is _t_int           *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    eqData_t *Right = (eqData_t *)In;
    witnessData_t_int *Left = (witnessData_t_int *)Out;

    if (MPType == 0)
    { // count the total number of comp_d to be sent
      Left->num_comp_d = Right->witnessData_d[stage].num_linears * Right->witnessData_d[stage].depth + Right->witnessData_d[stage].B->rows * Right->witnessData_d[stage].B->cols + Right->witnessData_d[stage].p->size;

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
      count = 0;
      for (i = 0; i < Right->witnessData_d[stage].num_linears; i++)
        for (j = 0; j < Right->witnessData_d[stage].depth; j++)
        {
          set_d((*coeff_d)[count], Right->witnessData_d[stage].startSystemCoeff[i][j]);
          count++;
        }
      for (i = 0; i < Right->witnessData_d[stage].B->rows; i++)
        for (j = 0; j < Right->witnessData_d[stage].B->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->witnessData_d[stage].B->entry[i][j]);
          count++;
        }
      for (i = 0; i < Right->witnessData_d[stage].p->size; i++)
      {
        set_d((*coeff_d)[count], &Right->witnessData_d[stage].p->coord[i]);
        count++;
      }

      // setup totalLength
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup wdStr
      *wdStr = NULL;
      strLoc = 0;
      for (i = 0; i < Right->witnessData_mp[stage].num_linears; i++)
        for (j = 0; j < Right->witnessData_mp[stage].depth; j++)
        { // real part
          strR = mpf_to_str(Right->witnessData_mp[stage].startSystemCoeff[i][j]->r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->witnessData_mp[stage].startSystemCoeff[i][j]->i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update wStr
          *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*wdStr)[strLoc], strR);
          strcpy(&(*wdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      for (i = 0; i < Right->witnessData_mp[stage].B->rows; i++)
        for (j = 0; j < Right->witnessData_mp[stage].B->cols; j++)
        { // real part
          strR = mpf_to_str(Right->witnessData_mp[stage].B->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->witnessData_mp[stage].B->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update wStr
          *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*wdStr)[strLoc], strR);
          strcpy(&(*wdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      for (i = 0; i < Right->witnessData_mp[stage].p->size; i++)
      { // real part
        strR = mpf_to_str(Right->witnessData_mp[stage].p->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->witnessData_mp[stage].p->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update wStr
        *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*wdStr)[strLoc], strR);
        strcpy(&(*wdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else // MPType == 2
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup wdStr
      *wdStr = NULL;
      strLoc = 0;
      for (i = 0; i < Right->witnessData_mp[stage].num_linears; i++)
        for (j = 0; j < Right->witnessData_mp[stage].depth; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->witnessData_mp[stage].startSystemCoeff_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->witnessData_mp[stage].startSystemCoeff_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update wStr
          *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*wdStr)[strLoc], strR);
          strcpy(&(*wdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      for (i = 0; i < Right->witnessData_mp[stage].B->rows; i++)
        for (j = 0; j < Right->witnessData_mp[stage].B->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->witnessData_mp[stage].B_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->witnessData_mp[stage].B_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update wStr
          *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*wdStr)[strLoc], strR);
          strcpy(&(*wdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      for (i = 0; i < Right->witnessData_mp[stage].p->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->witnessData_mp[stage].p_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->witnessData_mp[stage].p_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update wStr
        *wdStr = (char *)brealloc(*wdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*wdStr)[strLoc], strR);
        strcpy(&(*wdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      Left->totalLength = strLoc;
    }

    // copy over the rest of the values
    Left->curr_precision = Right->curr_precision;
    if (MPType == 0 || MPType == 2)
    {
      Left->startFunction = Right->witnessData_d[stage].startFunction;
      Left->depth = Right->witnessData_d[stage].depth;
      Left->num_paths = Right->witnessData_d[stage].num_paths;
      Left->num_linears = Right->witnessData_d[stage].num_linears;
      Left->B_rows = Right->witnessData_d[stage].B->rows;
      Left->B_cols = Right->witnessData_d[stage].B->cols;
      Left->p_size = Right->witnessData_d[stage].p->size;
    }
    else
    {
      Left->startFunction = Right->witnessData_mp[stage].startFunction;
      Left->depth = Right->witnessData_mp[stage].depth;
      Left->num_paths = Right->witnessData_mp[stage].num_paths;
      Left->num_linears = Right->witnessData_mp[stage].num_linears;
      Left->B_rows = Right->witnessData_mp[stage].B->rows;
      Left->B_cols = Right->witnessData_mp[stage].B->cols;
      Left->p_size = Right->witnessData_mp[stage].p->size;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    eqData_t *Left = (eqData_t *)Out;
    witnessData_t_int *Right = (witnessData_t_int *)In;

    if (MPType == 0)
    { // setup all of the comp_d
      count = 0;
      Left->witnessData_d[stage].startSystemCoeff = (comp_d **)bmalloc(Right->num_linears * sizeof(comp_d *));
      for (i = 0; i < Right->num_linears; i++)
      {
        Left->witnessData_d[stage].startSystemCoeff[i] = (comp_d *)bmalloc(Right->depth * sizeof(comp_d));
        for (j = 0; j < Right->depth; j++)
        {
          set_d(Left->witnessData_d[stage].startSystemCoeff[i][j], (*coeff_d)[count]);
          count++;
        }
      }
      init_mat_d(Left->witnessData_d[stage].B, Right->B_rows, Right->B_cols);
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)        
        {
          set_d(&Left->witnessData_d[stage].B->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      init_vec_d(Left->witnessData_d[stage].p, Right->p_size);
      for (i = 0; i < Right->p_size; i++)
      {
        set_d(&Left->witnessData_d[stage].p->coord[i], (*coeff_d)[count]);
        count++;
      }
      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup all of the comp_mp
      Left->witnessData_mp[stage].startSystemCoeff = (comp_mp **)bmalloc(Right->num_linears * sizeof(comp_mp *));
      strLoc = 0;
      for (i = 0; i < Right->num_linears; i++)
      {
        Left->witnessData_mp[stage].startSystemCoeff[i] = (comp_mp *)bmalloc(Right->depth * sizeof(comp_mp));
        for (j = 0; j < Right->depth; j++)
        {
          init_mp2(Left->witnessData_mp[stage].startSystemCoeff[i][j], Right->curr_precision);
          // setup real part
          mpf_set_str(Left->witnessData_mp[stage].startSystemCoeff[i][j]->r, &(*wdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->witnessData_mp[stage].startSystemCoeff[i][j]->i, &(*wdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
        }
      }
      init_mat_mp2(Left->witnessData_mp[stage].B, Right->B_rows, Right->B_cols, Right->curr_precision);
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)
        { // setup real part
          mpf_set_str(Left->witnessData_mp[stage].B->entry[i][j].r, &(*wdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->witnessData_mp[stage].B->entry[i][j].i, &(*wdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
        }
      init_vec_mp2(Left->witnessData_mp[stage].p, Right->p_size, Right->curr_precision);
      for (i = 0; i < Right->p_size; i++)
      { // setup real part
        mpf_set_str(Left->witnessData_mp[stage].p->coord[i].r, &(*wdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*wdStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->witnessData_mp[stage].p->coord[i].i, &(*wdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*wdStr)[strLoc]);
      }

      // free str
      if (freeStr)
        free(*wdStr);
    }
    else // MPType == 2
    { // setup for AMP tracking
      Left->witnessData_d[stage].startSystemCoeff = (comp_d **)bmalloc(Right->num_linears * sizeof(comp_d *));
      Left->witnessData_mp[stage].startSystemCoeff = (comp_mp **)bmalloc(Right->num_linears * sizeof(comp_mp *));
      Left->witnessData_mp[stage].startSystemCoeff_rat = (mpq_t ***)bmalloc(Right->num_linears * sizeof(mpq_t **));
      strLoc = 0;
      for (i = 0; i < Right->num_linears; i++)
      {
        Left->witnessData_d[stage].startSystemCoeff[i] = (comp_d *)bmalloc(Right->depth * sizeof(comp_d));
        Left->witnessData_mp[stage].startSystemCoeff[i] = (comp_mp *)bmalloc(Right->depth * sizeof(comp_mp));
        Left->witnessData_mp[stage].startSystemCoeff_rat[i] = (mpq_t **)bmalloc(Right->depth * sizeof(mpq_t *));
        for (j = 0; j < Right->depth; j++)
        {
          Left->witnessData_mp[stage].startSystemCoeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          // setup real part
          mpq_init(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][0]);
          mpq_set_str(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][0], &(*wdStr)[strLoc], base);
          mpq_canonicalize(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][0]);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
          // setup imag part
          mpq_init(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][1]);
          mpq_set_str(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][1], &(*wdStr)[strLoc], base);
          mpq_canonicalize(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][1]);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);

          // setup _mp & _d
          init_mp2(Left->witnessData_mp[stage].startSystemCoeff[i][j], Right->curr_precision);
          mpf_set_q(Left->witnessData_mp[stage].startSystemCoeff[i][j]->r, Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][0]);
          Left->witnessData_d[stage].startSystemCoeff[i][j]->r = mpq_get_d(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][0]);

          mpf_set_q(Left->witnessData_mp[stage].startSystemCoeff[i][j]->i, Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][1]);
          Left->witnessData_d[stage].startSystemCoeff[i][j]->i = mpq_get_d(Left->witnessData_mp[stage].startSystemCoeff_rat[i][j][1]);
        }
      }
      init_mat_d(Left->witnessData_d[stage].B, Right->B_rows, Right->B_cols);
      init_mat_mp2(Left->witnessData_mp[stage].B, Right->B_rows, Right->B_cols, Right->curr_precision);
      init_mat_rat(Left->witnessData_mp[stage].B_rat, Right->B_rows, Right->B_cols);
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)
        { // setup real part
          mpq_init(Left->witnessData_mp[stage].B_rat[i][j][0]);
          mpq_set_str(Left->witnessData_mp[stage].B_rat[i][j][0], &(*wdStr)[strLoc], base);
          mpq_canonicalize(Left->witnessData_mp[stage].B_rat[i][j][0]);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);
          // setup imag part
          mpq_init(Left->witnessData_mp[stage].B_rat[i][j][1]);
          mpq_set_str(Left->witnessData_mp[stage].B_rat[i][j][1], &(*wdStr)[strLoc], base);
          mpq_canonicalize(Left->witnessData_mp[stage].B_rat[i][j][1]);
          strLoc += 1 + strlen(&(*wdStr)[strLoc]);

          // setup _mp & _d
          mpf_set_q(Left->witnessData_mp[stage].B->entry[i][j].r, Left->witnessData_mp[stage].B_rat[i][j][0]);
          Left->witnessData_d[stage].B->entry[i][j].r = mpq_get_d(Left->witnessData_mp[stage].B_rat[i][j][0]);

          mpf_set_q(Left->witnessData_mp[stage].B->entry[i][j].i, Left->witnessData_mp[stage].B_rat[i][j][1]);
          Left->witnessData_d[stage].B->entry[i][j].i = mpq_get_d(Left->witnessData_mp[stage].B_rat[i][j][1]);
        }
      init_vec_d(Left->witnessData_d[stage].p, Right->p_size);
      init_vec_mp2(Left->witnessData_mp[stage].p, Right->p_size, Right->curr_precision);
      init_vec_rat(Left->witnessData_mp[stage].p_rat, Right->p_size);
      for (i = 0; i < Right->p_size; i++)
      { // setup real part
        mpq_init(Left->witnessData_mp[stage].p_rat[i][0]);
        mpq_set_str(Left->witnessData_mp[stage].p_rat[i][0], &(*wdStr)[strLoc], base);
        mpq_canonicalize(Left->witnessData_mp[stage].p_rat[i][0]);
        strLoc += 1 + strlen(&(*wdStr)[strLoc]);
        // setup imag part
        mpq_init(Left->witnessData_mp[stage].p_rat[i][1]);
        mpq_set_str(Left->witnessData_mp[stage].p_rat[i][1], &(*wdStr)[strLoc], base);
        mpq_canonicalize(Left->witnessData_mp[stage].p_rat[i][1]);
        strLoc += 1 + strlen(&(*wdStr)[strLoc]);

        // setup _mp & _d
        mpf_set_q(Left->witnessData_mp[stage].p->coord[i].r, Left->witnessData_mp[stage].p_rat[i][0]);
        Left->witnessData_d[stage].p->coord[i].r = mpq_get_d(Left->witnessData_mp[stage].p_rat[i][0]);

        mpf_set_q(Left->witnessData_mp[stage].p->coord[i].i, Left->witnessData_mp[stage].p_rat[i][1]);
        Left->witnessData_d[stage].p->coord[i].i = mpq_get_d(Left->witnessData_mp[stage].p_rat[i][1]);
      }

      // free str
      if (freeStr)
        free(*wdStr);
    }

    // copy over the rest of the values
    if (MPType == 0 || MPType == 2)
    { // copy to _d
      Left->witnessData_d[stage].startFunction = Right->startFunction;
      Left->witnessData_d[stage].depth = Right->depth;
      Left->witnessData_d[stage].num_paths = Right->num_paths;
      Left->witnessData_d[stage].num_linears = Right->num_linears;
      Left->witnessData_d[stage].B->rows = Right->B_rows;
      Left->witnessData_d[stage].B->cols = Right->B_cols;
      Left->witnessData_d[stage].p->size = Right->p_size;

      // set all other pointers to NULL
      Left->witnessData_d[stage].startPts = NULL;
      Left->witnessData_d[stage].endPts_in = NULL;
      Left->witnessData_d[stage].endPts = NULL;
      Left->witnessData_d[stage].finalTs = NULL;
      Left->witnessData_d[stage].condition_nums = NULL;
      Left->witnessData_d[stage].endPt_retVals = NULL;
      Left->witnessData_d[stage].endPt_types = NULL;
      Left->witnessData_d[stage].higherDim = NULL;
    }

    if (MPType == 1 || MPType == 2)
    { // copy to _mp
      Left->witnessData_mp[stage].startFunction = Right->startFunction;
      Left->witnessData_mp[stage].depth = Right->depth;
      Left->witnessData_mp[stage].num_paths = Right->num_paths;
      Left->witnessData_mp[stage].num_linears = Right->num_linears;
      Left->witnessData_mp[stage].B->rows = Right->B_rows;
      Left->witnessData_mp[stage].B->cols = Right->B_cols;
      Left->witnessData_mp[stage].p->size = Right->p_size;

      // set all other pointers to NULL
      Left->witnessData_mp[stage].startPts = NULL;
      Left->witnessData_mp[stage].endPts_in = NULL;
      Left->witnessData_mp[stage].endPts = NULL;
      Left->witnessData_mp[stage].finalTs = NULL;
      Left->witnessData_mp[stage].condition_nums = NULL;
      Left->witnessData_mp[stage].endPt_retVals = NULL;
      Left->witnessData_mp[stage].endPt_types = NULL;
      Left->witnessData_mp[stage].higherDim = NULL;
    }
  }
  return;
}

void cp_stageData_int(void *Out, void *In, int stage, int MPType, char **sdStr, int freeStr, comp_d **coeff_d, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is EqD               *
*                otherwise - Out is EqD, In is _t_int           *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    eqData_t *Right = (eqData_t *)In;
    stageData_t_int *Left = (stageData_t_int *)Out;

    // setup gamma & slice
    if (MPType == 0)
    {  // find the number of comp_d needed
      if (Right->stageData_d[stage].useIntrinsicSlice)
        Left->num_comp_d = 1 + Right->stageData_d[stage].B->rows * Right->stageData_d[stage].B->cols + Right->stageData_d[stage].p->size + Right->stageData_d[stage].B1->rows * Right->stageData_d[stage].B1->cols + Right->stageData_d[stage].p1->size + Right->stageData_d[stage].B0->rows * Right->stageData_d[stage].B0->cols + Right->stageData_d[stage].p0->size;
      else
        Left->num_comp_d = 1; // only gamma

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));

      count = 0;
      // setup gamma in coeff_d
      set_d((*coeff_d)[count], Right->gamma_d);
      count++;

      // setup slice, if needed
      if (Right->stageData_d[stage].useIntrinsicSlice)
      { // setup B in coeff_d
        for (i = 0; i < Right->stageData_d[stage].B->rows; i++)
          for (j = 0; j < Right->stageData_d[stage].B->cols; j++)
          {
            set_d((*coeff_d)[count], &Right->stageData_d[stage].B->entry[i][j]);
            count++;
          }
        // setup p in coeff_d
        for (i = 0; i < Right->stageData_d[stage].p->size; i++)
        {
          set_d((*coeff_d)[count], &Right->stageData_d[stage].p->coord[i]);
          count++;
        }

        // setup B1 in coeff_d
        for (i = 0; i < Right->stageData_d[stage].B1->rows; i++)
          for (j = 0; j < Right->stageData_d[stage].B1->cols; j++)
          {
            set_d((*coeff_d)[count], &Right->stageData_d[stage].B1->entry[i][j]);
            count++;
          }
        // setup p1 in coeff_d
        for (i = 0; i < Right->stageData_d[stage].p1->size; i++)
        {
          set_d((*coeff_d)[count], &Right->stageData_d[stage].p1->coord[i]);
          count++;
        }

        // setup B0 in coeff_d
        for (i = 0; i < Right->stageData_d[stage].B0->rows; i++)
          for (j = 0; j < Right->stageData_d[stage].B0->cols; j++)
          {
            set_d((*coeff_d)[count], &Right->stageData_d[stage].B0->entry[i][j]);
            count++;
          }
        // setup p0 in coeff_d
        for (i = 0; i < Right->stageData_d[stage].p0->size; i++)
        {
          set_d((*coeff_d)[count], &Right->stageData_d[stage].p0->coord[i]);
          count++;
        }
      }
      // setup totalLength
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup sdStr
      *sdStr = NULL;
      strLoc = 0;
      // setup gamma in sdStr
      // real part of gamma
      strR = mpf_to_str(Right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part of gamma
      strI = mpf_to_str(Right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update sdStr
      *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*sdStr)[strLoc], strR);
      strcpy(&(*sdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      if (Right->stageData_mp[stage].useIntrinsicSlice)
      { // setup B in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B->cols; j++)
          { // real part
            strR = mpf_to_str(Right->stageData_mp[stage].B->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->stageData_mp[stage].B->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'
  
            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;
  
            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p->size; i++)
        { // real part
          strR = mpf_to_str(Right->stageData_mp[stage].p->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->stageData_mp[stage].p->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'
  
          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

        // setup B1 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B1->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B1->cols; j++)
          { // real part
            strR = mpf_to_str(Right->stageData_mp[stage].B1->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->stageData_mp[stage].B1->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p1 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p1->size; i++)
        { // real part
          strR = mpf_to_str(Right->stageData_mp[stage].p1->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->stageData_mp[stage].p1->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

        // setup B0 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B0->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B0->cols; j++)
          { // real part
            strR = mpf_to_str(Right->stageData_mp[stage].B0->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->stageData_mp[stage].B0->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p0 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p0->size; i++)
        { // real part
          strR = mpf_to_str(Right->stageData_mp[stage].p0->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->stageData_mp[stage].p0->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else if (MPType == 2)
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup sdStr
      *sdStr = NULL;
      strLoc = 0;
      // setup gamma in sdStr
      // real part of gamma
      strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part of gamma
      strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update sdStr
      *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*sdStr)[strLoc], strR);
      strcpy(&(*sdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      if (Right->stageData_d[stage].useIntrinsicSlice)
      { // setup B in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B->cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].B_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].B_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p->size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].p_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].p_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;
  
          // free strR, strI
          free(strR);
          free(strI);
        }

        // setup B1 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B1->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B1->cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].B1_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].B1_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p1 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p1->size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].p1_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].p1_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

        // setup B0 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].B0->rows; i++)
          for (j = 0; j < Right->stageData_mp[stage].B0->cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].B0_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].B0_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update sdStr
            *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*sdStr)[strLoc], strR);
            strcpy(&(*sdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // setup p0 in sdStr
        for (i = 0; i < Right->stageData_mp[stage].p0->size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->stageData_mp[stage].p0_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->stageData_mp[stage].p0_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sdStr
          *sdStr = (char *)brealloc(*sdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sdStr)[strLoc], strR);
          strcpy(&(*sdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }

    // copy over the rest of the values
    Left->curr_precision = Right->curr_precision;
    if (MPType == 0 || MPType == 2)
    {
      Left->depth_x = Right->stageData_d[stage].depth_x;
      Left->depth_y = Right->stageData_d[stage].depth_y;
      Left->num_paths = Right->stageData_d[stage].num_paths;
      Left->useIntrinsicSlice = Right->stageData_d[stage].useIntrinsicSlice;
      Left->B_rows = Right->stageData_d[stage].B->rows;
      Left->B_cols = Right->stageData_d[stage].B->cols;
      Left->p_size = Right->stageData_d[stage].p->size;
      Left->Bt_rows = Right->stageData_d[stage].B1->rows;
      Left->Bt_cols = Right->stageData_d[stage].B1->cols;
      Left->pt_size = Right->stageData_d[stage].p1->size;
    }
    else
    {
      Left->depth_x = Right->stageData_mp[stage].depth_x;
      Left->depth_y = Right->stageData_mp[stage].depth_y;
      Left->num_paths = Right->stageData_mp[stage].num_paths;
      Left->useIntrinsicSlice = Right->stageData_mp[stage].useIntrinsicSlice;
      Left->B_rows = Right->stageData_mp[stage].B->rows;
      Left->B_cols = Right->stageData_mp[stage].B->cols;
      Left->p_size = Right->stageData_mp[stage].p->size;
      Left->Bt_rows = Right->stageData_mp[stage].B1->rows;
      Left->Bt_cols = Right->stageData_mp[stage].B1->cols;
      Left->pt_size = Right->stageData_mp[stage].p1->size;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    eqData_t *Left = (eqData_t *)Out;
    stageData_t_int *Right = (stageData_t_int *)In;

    if (MPType == 0)
    { // setup all of the comp_d
      count = 0;
      // setup gamma_d
      set_d(Left->gamma_d, (*coeff_d)[count]);
      count++;

      if (Right->useIntrinsicSlice)
      { // setup B
        init_mat_d(Left->stageData_d[stage].B, Right->B_rows, Right->B_cols);
        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          {
            set_d(&Left->stageData_d[stage].B->entry[i][j], (*coeff_d)[count]);
            count++;
          }
        // setup p
        init_vec_d(Left->stageData_d[stage].p, Right->p_size);
        for (i = 0; i < Right->p_size; i++)
        {
          set_d(&Left->stageData_d[stage].p->coord[i], (*coeff_d)[count]);
          count++;
        }

        // setup B1
        init_mat_d(Left->stageData_d[stage].B1, Right->Bt_rows, Right->Bt_cols);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          {
            set_d(&Left->stageData_d[stage].B1->entry[i][j], (*coeff_d)[count]);
            count++;
          }
        // setup p1
        init_vec_d(Left->stageData_d[stage].p1, Right->pt_size);
        for (i = 0; i < Right->pt_size; i++)
        {
          set_d(&Left->stageData_d[stage].p1->coord[i], (*coeff_d)[count]);
          count++;
        }

        // setup B0
        init_mat_d(Left->stageData_d[stage].B0, Right->Bt_rows, Right->Bt_cols);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          {
            set_d(&Left->stageData_d[stage].B0->entry[i][j], (*coeff_d)[count]);
            count++;
          }
        // setup p1
        init_vec_d(Left->stageData_d[stage].p0, Right->pt_size);
        for (i = 0; i < Right->pt_size; i++)
        {
          set_d(&Left->stageData_d[stage].p0->coord[i], (*coeff_d)[count]);
          count++;
        }
      }  
      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup all of the comp_mp
      strLoc = 0;
      // setup gamma_mp
      // setup real part
      mpf_set_str(Left->gamma_mp->r, &(*sdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*sdStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->gamma_mp->i, &(*sdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*sdStr)[strLoc]);

      if (Right->useIntrinsicSlice)
      { // setup B
        init_mat_mp2(Left->stageData_mp[stage].B, Right->B_rows, Right->B_cols, Right->curr_precision);
        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          { // setup real part
            mpf_set_str(Left->stageData_mp[stage].B->entry[i][j].r, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->stageData_mp[stage].B->entry[i][j].i, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          }
        // setup p
        init_vec_mp2(Left->stageData_mp[stage].p, Right->p_size, Right->curr_precision);
        for (i = 0; i < Right->p_size; i++)
        { // setup real part
          mpf_set_str(Left->stageData_mp[stage].p->coord[i].r, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->stageData_mp[stage].p->coord[i].i, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
        }

        // setup B1
        init_mat_mp2(Left->stageData_mp[stage].B1, Right->Bt_rows, Right->Bt_cols, Right->curr_precision);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          { // setup real part
            mpf_set_str(Left->stageData_mp[stage].B1->entry[i][j].r, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->stageData_mp[stage].B1->entry[i][j].i, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          }
        // setup p1
        init_vec_mp2(Left->stageData_mp[stage].p1, Right->pt_size, Right->curr_precision);
        for (i = 0; i < Right->pt_size; i++)
        { // setup real part
          mpf_set_str(Left->stageData_mp[stage].p1->coord[i].r, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->stageData_mp[stage].p1->coord[i].i, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
        }

        // setup B0
        init_mat_mp2(Left->stageData_mp[stage].B0, Right->Bt_rows, Right->Bt_cols, Right->curr_precision);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          { // setup real part
            mpf_set_str(Left->stageData_mp[stage].B0->entry[i][j].r, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->stageData_mp[stage].B0->entry[i][j].i, &(*sdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          }
        // setup p0
        init_vec_mp2(Left->stageData_mp[stage].p0, Right->pt_size, Right->curr_precision);
        for (i = 0; i < Right->pt_size; i++)
        { // setup real part
          mpf_set_str(Left->stageData_mp[stage].p0->coord[i].r, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->stageData_mp[stage].p0->coord[i].i, &(*sdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
        }
      }

      // free str
      if (freeStr)
        free(*sdStr);
    }
    else // MPType == 2
    { // setup for AMP tracking
      strLoc = 0;
      // setup gamma
      // setup real part
      mpq_set_str(Left->gamma_rat[0], &(*sdStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*sdStr)[strLoc]);
      // setup imag part
      mpq_set_str(Left->gamma_rat[1], &(*sdStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[1]);
      strLoc += 1 + strlen(&(*sdStr)[strLoc]);
      // setup _mp & _d
      mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
      Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);
      mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
      Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);

      if (Right->useIntrinsicSlice)
      { // setup B
        init_mat_d(Left->stageData_d[stage].B, Right->B_rows, Right->B_cols);
        init_mat_mp2(Left->stageData_mp[stage].B, Right->B_rows, Right->B_cols, Right->curr_precision);
        init_mat_rat(Left->stageData_mp[stage].B_rat, Right->B_rows, Right->B_cols);
        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          { // setup real part
            mpq_set_str(Left->stageData_mp[stage].B_rat[i][j][0], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B_rat[i][j][0]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->stageData_mp[stage].B_rat[i][j][1], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B_rat[i][j][1]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
  
            // setup _mp & _d
            mpf_set_q(Left->stageData_mp[stage].B->entry[i][j].r, Left->stageData_mp[stage].B_rat[i][j][0]);
            Left->stageData_d[stage].B->entry[i][j].r = mpq_get_d(Left->stageData_mp[stage].B_rat[i][j][0]);
  
            mpf_set_q(Left->stageData_mp[stage].B->entry[i][j].i, Left->stageData_mp[stage].B_rat[i][j][1]);
            Left->stageData_d[stage].B->entry[i][j].i = mpq_get_d(Left->stageData_mp[stage].B_rat[i][j][1]);
          }
        // setup p
        init_vec_d(Left->stageData_d[stage].p, Right->p_size);
        init_vec_mp2(Left->stageData_mp[stage].p, Right->p_size, Right->curr_precision);
        init_vec_rat(Left->stageData_mp[stage].p_rat, Right->p_size);
        for (i = 0; i < Right->p_size; i++)
        { // setup real part
          mpq_set_str(Left->stageData_mp[stage].p_rat[i][0], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p_rat[i][0]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->stageData_mp[stage].p_rat[i][1], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p_rat[i][1]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
  
          // setup _mp & _d
          mpf_set_q(Left->stageData_mp[stage].p->coord[i].r, Left->stageData_mp[stage].p_rat[i][0]);
          Left->stageData_d[stage].p->coord[i].r = mpq_get_d(Left->stageData_mp[stage].p_rat[i][0]);
  
          mpf_set_q(Left->stageData_mp[stage].p->coord[i].i, Left->stageData_mp[stage].p_rat[i][1]);
          Left->stageData_d[stage].p->coord[i].i = mpq_get_d(Left->stageData_mp[stage].p_rat[i][1]);
        }

        // setup B1
        init_mat_d(Left->stageData_d[stage].B1, Right->Bt_rows, Right->Bt_cols);
        init_mat_mp2(Left->stageData_mp[stage].B1, Right->Bt_rows, Right->Bt_cols, Right->curr_precision);
        init_mat_rat(Left->stageData_mp[stage].B1_rat, Right->Bt_rows, Right->Bt_cols);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          { // setup real part
            mpq_set_str(Left->stageData_mp[stage].B1_rat[i][j][0], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B1_rat[i][j][0]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->stageData_mp[stage].B1_rat[i][j][1], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B1_rat[i][j][1]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);

            // setup _mp & _d
            mpf_set_q(Left->stageData_mp[stage].B1->entry[i][j].r, Left->stageData_mp[stage].B1_rat[i][j][0]);
            Left->stageData_d[stage].B1->entry[i][j].r = mpq_get_d(Left->stageData_mp[stage].B1_rat[i][j][0]);

            mpf_set_q(Left->stageData_mp[stage].B1->entry[i][j].i, Left->stageData_mp[stage].B1_rat[i][j][1]);
            Left->stageData_d[stage].B1->entry[i][j].i = mpq_get_d(Left->stageData_mp[stage].B1_rat[i][j][1]);
          }
        // setup p1
        init_vec_d(Left->stageData_d[stage].p1, Right->pt_size);
        init_vec_mp2(Left->stageData_mp[stage].p1, Right->pt_size, Right->curr_precision);
        init_vec_rat(Left->stageData_mp[stage].p1_rat, Right->pt_size);
        for (i = 0; i < Right->pt_size; i++)
        { // setup real part
          mpq_set_str(Left->stageData_mp[stage].p1_rat[i][0], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p1_rat[i][0]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->stageData_mp[stage].p1_rat[i][1], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p1_rat[i][1]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);

          // setup _mp & _d
          mpf_set_q(Left->stageData_mp[stage].p1->coord[i].r, Left->stageData_mp[stage].p1_rat[i][0]);
          Left->stageData_d[stage].p1->coord[i].r = mpq_get_d(Left->stageData_mp[stage].p1_rat[i][0]);

          mpf_set_q(Left->stageData_mp[stage].p1->coord[i].i, Left->stageData_mp[stage].p1_rat[i][1]);
          Left->stageData_d[stage].p1->coord[i].i = mpq_get_d(Left->stageData_mp[stage].p1_rat[i][1]);
        }

        // setup B0
        init_mat_d(Left->stageData_d[stage].B0, Right->Bt_rows, Right->Bt_cols);
        init_mat_mp2(Left->stageData_mp[stage].B0, Right->Bt_rows, Right->Bt_cols, Right->curr_precision);
        init_mat_rat(Left->stageData_mp[stage].B0_rat, Right->Bt_rows, Right->Bt_cols);
        for (i = 0; i < Right->Bt_rows; i++)
          for (j = 0; j < Right->Bt_cols; j++)
          { // setup real part
            mpq_set_str(Left->stageData_mp[stage].B0_rat[i][j][0], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B0_rat[i][j][0]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->stageData_mp[stage].B0_rat[i][j][1], &(*sdStr)[strLoc], base);
            mpq_canonicalize(Left->stageData_mp[stage].B0_rat[i][j][1]);
            strLoc += 1 + strlen(&(*sdStr)[strLoc]);

            // setup _mp & _d
            mpf_set_q(Left->stageData_mp[stage].B0->entry[i][j].r, Left->stageData_mp[stage].B0_rat[i][j][0]);
            Left->stageData_d[stage].B0->entry[i][j].r = mpq_get_d(Left->stageData_mp[stage].B0_rat[i][j][0]);

            mpf_set_q(Left->stageData_mp[stage].B0->entry[i][j].i, Left->stageData_mp[stage].B0_rat[i][j][1]);
            Left->stageData_d[stage].B0->entry[i][j].i = mpq_get_d(Left->stageData_mp[stage].B0_rat[i][j][1]);
          }
        // setup p0
        init_vec_d(Left->stageData_d[stage].p0, Right->pt_size);
        init_vec_mp2(Left->stageData_mp[stage].p0, Right->pt_size, Right->curr_precision);
        init_vec_rat(Left->stageData_mp[stage].p0_rat, Right->pt_size);
        for (i = 0; i < Right->pt_size; i++)
        { // setup real part
          mpq_set_str(Left->stageData_mp[stage].p0_rat[i][0], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p0_rat[i][0]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->stageData_mp[stage].p0_rat[i][1], &(*sdStr)[strLoc], base);
          mpq_canonicalize(Left->stageData_mp[stage].p0_rat[i][1]);
          strLoc += 1 + strlen(&(*sdStr)[strLoc]);

          // setup _mp & _d
          mpf_set_q(Left->stageData_mp[stage].p0->coord[i].r, Left->stageData_mp[stage].p0_rat[i][0]);
          Left->stageData_d[stage].p0->coord[i].r = mpq_get_d(Left->stageData_mp[stage].p0_rat[i][0]);

          mpf_set_q(Left->stageData_mp[stage].p0->coord[i].i, Left->stageData_mp[stage].p0_rat[i][1]);
          Left->stageData_d[stage].p0->coord[i].i = mpq_get_d(Left->stageData_mp[stage].p0_rat[i][1]);
        }
      }

      // free str
      if (freeStr)
        free(*sdStr);
    }

    // copy over the rest of the values
    if (MPType == 0 || MPType == 2)
    { // copy to _d
      Left->stageData_d[stage].depth_x = Right->depth_x;
      Left->stageData_d[stage].depth_y = Right->depth_y;
      Left->stageData_d[stage].num_paths = Right->num_paths;
      Left->stageData_d[stage].useIntrinsicSlice = Right->useIntrinsicSlice;
      Left->stageData_d[stage].B->rows = Right->B_rows;
      Left->stageData_d[stage].B->cols = Right->B_cols;
      Left->stageData_d[stage].p->size = Right->p_size;
      Left->stageData_d[stage].B1->rows = Left->stageData_d[stage].B0->rows = Right->Bt_rows;
      Left->stageData_d[stage].B1->cols = Left->stageData_d[stage].B0->cols = Right->Bt_cols;
      Left->stageData_d[stage].p1->size = Left->stageData_d[stage].p0->size = Right->pt_size;

      // set all other pointers to NULL
      Left->stageData_d[stage].startPts = NULL;
      Left->stageData_d[stage].endPts_in = NULL;
      Left->stageData_d[stage].endPts = NULL;
      Left->stageData_d[stage].finalTs = NULL;
      Left->stageData_d[stage].condition_nums = NULL;
      Left->stageData_d[stage].endPt_retVals = NULL;
      Left->stageData_d[stage].endPt_types = NULL;
      Left->stageData_d[stage].higherDim = NULL;
    }

    if (MPType == 1 || MPType == 2)
    { // copy to _mp
      Left->stageData_mp[stage].depth_x = Right->depth_x;
      Left->stageData_mp[stage].depth_y = Right->depth_y;
      Left->stageData_mp[stage].num_paths = Right->num_paths;
      Left->stageData_mp[stage].useIntrinsicSlice = Right->useIntrinsicSlice;
      Left->stageData_mp[stage].B->rows = Right->B_rows;
      Left->stageData_mp[stage].B->cols = Right->B_cols;
      Left->stageData_mp[stage].p->size = Right->p_size;
      Left->stageData_mp[stage].B1->rows = Left->stageData_mp[stage].B0->rows = Right->Bt_rows;
      Left->stageData_mp[stage].B1->cols = Left->stageData_mp[stage].B0->cols = Right->Bt_cols;
      Left->stageData_mp[stage].p1->size = Left->stageData_mp[stage].p0->size = Right->pt_size;

      // set all other pointers to NULL
      Left->stageData_mp[stage].startPts = NULL;
      Left->stageData_mp[stage].endPts_in = NULL;
      Left->stageData_mp[stage].endPts = NULL;
      Left->stageData_mp[stage].finalTs = NULL;
      Left->stageData_mp[stage].condition_nums = NULL;
      Left->stageData_mp[stage].endPt_retVals = NULL;
      Left->stageData_mp[stage].endPt_types = NULL;
      Left->stageData_mp[stage].higherDim = NULL;
    }
  }
  return;
}

void cp_codim_int(void *Out, void *In, int MPType, char **codimStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    codim_t *Right = (codim_t *)In;
    codim_t_int *Left = (codim_t_int *)Out;

    // setup Prog
    cp_prog_t_int(&Left->Prog_int, Right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD_int, &Right->PPD, ppd_type, ppd_size, inType);

    // setup degrees - orig_degrees, new_degrees & P
    *degrees = (int *)bmalloc(3 * Right->num_funcs * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      (*degrees)[i] = Right->orig_degrees[i];
      (*degrees)[Right->num_funcs + i] = Right->new_degrees[i];
      (*degrees)[2 * Right->num_funcs + i] = Right->P[i];
    }

    // setup other values
    Left->curr_precision = Right->curr_precision;
    Left->num_funcs = Right->num_funcs;
    Left->system_rank = Right->system_rank;
    Left->num_codim = Right->num_codim;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup the rest
    if (MPType == 0)
    { // setup coeff_d
      if (Right->orig_variables != Right->new_variables)
      { // need to send gamma_d & C_d
        Left->num_comp_d = 1 + Right->C_d->rows * Right->C_d->cols;

        // setup coeff_d
        *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
        count = 0;
        // setup gamma first
        set_d((*coeff_d)[count], Right->gamma_d);
        count++;
        // setup C_d
        for (i = 0; i < Right->C_d->rows; i++)
          for (j = 0; j < Right->C_d->cols; j++)
          {
            set_d((*coeff_d)[count], &Right->C_d->entry[i][j]);
            count++;
          }
      }
      else
      { // need to only send gamma_d
        Left->num_comp_d = 1;

        // setup coeff_d
        *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
        count = 0;
        set_d((*coeff_d)[count], Right->gamma_d);
        count++;
      }

      // set totalLength = 0
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup codimStr
      *codimStr = NULL;
      strLoc = 0;

      // setup gamma
      // real part
      strR = mpf_to_str(Right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update codimStr
      *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*codimStr)[strLoc], strR);
      strcpy(&(*codimStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup C_mp, if needed
      if (Right->orig_variables != Right->new_variables)
      { // setup C_mp
        for (i = 0; i < Right->C_mp->rows; i++)
          for (j = 0; j < Right->C_mp->cols; j++)
          { // real part
            strR = mpf_to_str(Right->C_mp->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->C_mp->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update codimStr
            *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*codimStr)[strLoc], strR);
            strcpy(&(*codimStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else // MPType == 2
    {
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      // setup codimStr
      *codimStr = NULL;
      strLoc = 0;

      // setup gamma
      // real part
      strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update codimStr
      *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*codimStr)[strLoc], strR);
      strcpy(&(*codimStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup C_rat, if needed
      if (Right->orig_variables != Right->new_variables)
      { // setup C_rat
        for (i = 0; i < Right->C_mp->rows; i++)
          for (j = 0; j < Right->C_mp->cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->C_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->C_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update codimStr
            *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*codimStr)[strLoc], strR);
            strcpy(&(*codimStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    codim_t *Left = (codim_t *)Out;
    codim_t_int *Right = (codim_t_int *)In;

    // setup Prog
    Left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    cp_prog_t_int(Left->Prog, &Right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    // initialize evalProg
    initEvalProg(MPType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD, &Right->PPD_int, ppd_type, ppd_size, inType);

    // setup other values
    Left->curr_precision = Right->curr_precision;
    Left->num_funcs = Right->num_funcs;
    Left->system_rank = Right->system_rank;
    Left->num_codim = Right->num_codim;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup orig_degrees, new_degrees & P from degrees
    Left->orig_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->new_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->P = (int *)bmalloc(Right->num_funcs * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      Left->orig_degrees[i] = (*degrees)[i];
      Left->new_degrees[i] = (*degrees)[Right->num_funcs + i];
      Left->P[i] = (*degrees)[2 * Right->num_funcs + i];
    }

    free(*degrees);

    // setup the rest
    if (MPType == 0)
    { // setup gamma_d
      count = 0;
      set_d(Left->gamma_d, (*coeff_d)[count]);
      count++;

      // see if C_d needs setup
      if (Left->orig_variables != Left->new_variables)
      { // setup C_d
        init_mat_d(Left->C_d, Left->orig_variables, Left->new_variables);
        Left->C_d->rows = Left->orig_variables;
        Left->C_d->cols = Left->new_variables;
        for (i = 0; i < Left->C_d->rows; i++)
          for (j = 0; j < Left->C_d->cols; j++)
          {
            set_d(&Left->C_d->entry[i][j], (*coeff_d)[count]);
            count++;
          }
      }
      
      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup gamma_mp
      strLoc = 0;
      init_mp2(Left->gamma_mp, Left->curr_precision);
      // setup real part
      mpf_set_str(Left->gamma_mp->r, &(*codimStr)[strLoc], base);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->gamma_mp->i, &(*codimStr)[strLoc], base);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);

      // see if C_mp needs setup
      if (Left->orig_variables != Left->new_variables)
      { // setup C_mp
        init_mat_mp2(Left->C_mp, Left->orig_variables, Left->new_variables, Left->curr_precision);
        Left->C_mp->rows = Left->orig_variables;
        Left->C_mp->cols = Left->new_variables;
        for (i = 0; i < Left->C_mp->rows; i++)
          for (j = 0; j < Left->C_mp->cols; j++)
          { // setup real part
            mpf_set_str(Left->C_mp->entry[i][j].r, &(*codimStr)[strLoc], base);
            strLoc += 1 + strlen(&(*codimStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->C_mp->entry[i][j].i, &(*codimStr)[strLoc], base);
            strLoc += 1 + strlen(&(*codimStr)[strLoc]);
          }
      }

      // free str
      if (freeStr)
        free(*codimStr);
    }
    else // MPType == 2
    { // setup gamma_rat
      strLoc = 0;
      mpq_init(Left->gamma_rat[0]);
      mpq_set_str(Left->gamma_rat[0], &(*codimStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      mpq_init(Left->gamma_rat[1]);
      mpq_set_str(Left->gamma_rat[1], &(*codimStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[1]);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);

      // setup gamma_mp & gamma_d
      init_mp2(Left->gamma_mp, Left->curr_precision);
      mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
      Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);

      mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
      Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);

      // see if C needs setup
      if (Left->orig_variables != Left->new_variables)
      { // setup C_d, C_mp & C_rat
        init_mat_d(Left->C_d, Left->orig_variables, Left->new_variables);
        init_mat_mp2(Left->C_mp, Left->orig_variables, Left->new_variables, Left->curr_precision);
        init_mat_rat(Left->C_rat, Left->orig_variables, Left->new_variables);
        Left->C_d->rows = Left->C_mp->rows = Left->orig_variables;
        Left->C_d->cols = Left->C_mp->cols = Left->new_variables;
        for (i = 0; i < Left->C_mp->rows; i++)
          for (j = 0; j < Left->C_mp->cols; j++)
          { // setup real part
            mpq_set_str(Left->C_rat[i][j][0], &(*codimStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][0]);
            strLoc += 1 + strlen(&(*codimStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->C_rat[i][j][1], &(*codimStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][1]);
            strLoc += 1 + strlen(&(*codimStr)[strLoc]);

            // setup C_d & C_mp
            mpf_set_q(Left->C_mp->entry[i][j].r, Left->C_rat[i][j][0]);
            Left->C_d->entry[i][j].r = mpq_get_d(Left->C_rat[i][j][0]);

            mpf_set_q(Left->C_mp->entry[i][j].i, Left->C_rat[i][j][1]);
            Left->C_d->entry[i][j].i = mpq_get_d(Left->C_rat[i][j][1]);
          }
      }
      // free str
      if (freeStr)
        free(*codimStr);
    }
  }

  return;
}

void cp_codimData_int(void *Out, void *In, int curr_prec, int MPType, char **codimStr, int freeStr, comp_d **coeff_d, int **W_send, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    codimData_t *Right = (codimData_t *)In;
    codimData_t_int *Left = (codimData_t_int *)Out;

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    { // use _d to get the sizes
      Left->A_W_rows = Right->A_d->rows;
      Left->A_W_cols = Right->A_d->cols;
      Left->H_size = Right->H_d->size;
      Left->B_rows = Right->B_d->rows;
      Left->B_cols = Right->B_d->cols;
      Left->p_size = Right->p_d->size;
    }
    else
    { // use _mp to get the sizes
      Left->A_W_rows = Right->A_mp->rows;
      Left->A_W_cols = Right->A_mp->cols;
      Left->H_size = Right->H_mp->size;
      Left->B_rows = Right->B_mp->rows;
      Left->B_cols = Right->B_mp->cols;
      Left->p_size = Right->p_mp->size;
    }

    // setup other values
    Left->curr_precision = curr_prec;
    Left->codim = Right->codim;
    Left->useIntrinsicSlice = Right->useIntrinsicSlice;

    // setup W_send
    count = Left->A_W_rows * Left->A_W_cols;
    *W_send = (int *)bmalloc(count * sizeof(int));
    count = 0;
    for (i = 0; i < Left->A_W_rows; i++)
      for (j = 0; j < Left->A_W_cols; j++)
      {
        (*W_send)[count] = Right->W[i][j];
        count++;
      }

    // setup codimStr & coeff_d
    if (MPType == 0)
    { // setup num_comp_d
      Left->num_comp_d = 1 + Left->A_W_rows * Left->A_W_cols + Left->H_size + Left->B_rows * Left->B_cols + Left->p_size;
      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
      count = 0;
      // setup homVarConst_d
      set_d((*coeff_d)[count], Right->homVarConst_d);
      count++;
      // setup A_d
      for (i = 0; i < Right->A_d->rows; i++)
        for (j = 0; j < Right->A_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->A_d->entry[i][j]);
          count++;
        }
      // setup H_d
      for (i = 0; i < Right->H_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->H_d->coord[i]);
        count++;
      }
      // setup B_d
      for (i = 0; i < Right->B_d->rows; i++)
        for (j = 0; j < Right->B_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->B_d->entry[i][j]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Right->p_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->p_d->coord[i]);
        count++;
      }

      // set totalLength = 0
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    { // setup codimStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      *codimStr = NULL;
      strLoc = 0;

      // setup homVarConst_mp
      // real part
      strR = mpf_to_str(Right->homVarConst_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->homVarConst_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update codimStr
      *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*codimStr)[strLoc], strR);
      strcpy(&(*codimStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_mp
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->A_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->A_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update codimStr
          *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*codimStr)[strLoc], strR);
          strcpy(&(*codimStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup H_mp
      for (i = 0; i < Right->H_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->H_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->H_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update codimStr
        *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*codimStr)[strLoc], strR);
        strcpy(&(*codimStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup B_mp
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->B_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->B_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update codimStr
          *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*codimStr)[strLoc], strR);
          strcpy(&(*codimStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_mp
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->p_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->p_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update codimStr
        *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*codimStr)[strLoc], strR);
        strcpy(&(*codimStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else // MPType == 2
    { // setup codimStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      *codimStr = NULL;
      strLoc = 0;

      // setup homVarConst
      // real part
      strR = mpq_get_str(NULL, base, Right->homVarConst_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->homVarConst_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update codimStr
      *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*codimStr)[strLoc], strR);
      strcpy(&(*codimStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_rat
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->A_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->A_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update codimStr
          *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*codimStr)[strLoc], strR);
          strcpy(&(*codimStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup H_rat
      for (i = 0; i < Right->H_mp->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->H_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->H_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update codimStr
        *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*codimStr)[strLoc], strR);
        strcpy(&(*codimStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup B_rat
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->B_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->B_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update codimStr
          *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*codimStr)[strLoc], strR);
          strcpy(&(*codimStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_rat
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->p_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->p_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update codimStr
        *codimStr = (char *)brealloc(*codimStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*codimStr)[strLoc], strR);
        strcpy(&(*codimStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    codimData_t *Left = (codimData_t *)Out;
    codimData_t_int *Right = (codimData_t_int *)In;

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    { // setup _d sizes
      init_mat_d(Left->A_d, Right->A_W_rows, Right->A_W_cols);
      Left->A_d->rows = Right->A_W_rows;
      Left->A_d->cols = Right->A_W_cols;
      init_vec_d(Left->H_d, Right->H_size);
      Left->H_d->size = Right->H_size;
      init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
      Left->B_d->rows = Right->B_rows;
      Left->B_d->cols = Right->B_cols;
      init_vec_d(Left->p_d, Right->p_size);
      Left->p_d->size = Right->p_size;
    }

    if (MPType == 1 || MPType == 2)
    { // setup _mp sizes
      init_mat_mp2(Left->A_mp, Right->A_W_rows, Right->A_W_cols, Right->curr_precision);
      init_vec_mp2(Left->H_mp, Right->H_size, Right->curr_precision);
      init_mat_mp2(Left->B_mp, Right->B_rows, Right->B_cols, Right->curr_precision);
      init_vec_mp2(Left->p_mp, Right->p_size, Right->curr_precision);

      Left->A_mp->rows = Right->A_W_rows;
      Left->A_mp->cols = Right->A_W_cols;
      Left->H_mp->size = Right->H_size;
      Left->B_mp->rows = Right->B_rows;
      Left->B_mp->cols = Right->B_cols;
      Left->p_mp->size = Right->p_size;
    }

    // setup other values
    Left->codim = Right->codim;
    Left->useIntrinsicSlice = Right->useIntrinsicSlice;
    Left->num_paths = Left->num_superset = Left->num_nonsing = Left->num_sing = Left->num_nonsolns = Left->num_inf = Left->num_bad = 0;

    // NULL out other pointers
    Left->startPts_d = NULL;
    Left->startPts_mp = NULL;
    Left->endPts_d = NULL;
    Left->endPts_mp = NULL;
    Left->endPts_amp = NULL;
    Left->endPt_types = NULL;

    // setup W
    count = 0;
    Left->W = (int **)bmalloc(Right->A_W_rows * sizeof(int *));
    for (i = 0; i < Right->A_W_rows; i++)
    {
      Left->W[i] = (int *)bmalloc(Right->A_W_cols * sizeof(int));
      for (j = 0; j < Right->A_W_cols; j++)
      {
        Left->W[i][j] = (*W_send)[count];
        count++;
      }
    }

    // setup the other structures
    if (MPType == 0)
    { // initialize count
      count = 0;
      // setup homVarConst_d
      set_d(Left->homVarConst_d, (*coeff_d)[count]);
      count++;
      // setup A_d
      for (i = 0; i < Left->A_d->rows; i++)
        for (j = 0; j < Left->A_d->cols; j++)
        {
          set_d(&Left->A_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup H_d
      for (i = 0; i < Left->H_d->size; i++)
      {
        set_d(&Left->H_d->coord[i], (*coeff_d)[count]);
        count++;
      }
      // setup B_d
      for (i = 0; i < Left->B_d->rows; i++)
        for (j = 0; j < Left->B_d->cols; j++)
        {
          set_d(&Left->B_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Left->p_d->size; i++)
      {
        set_d(&Left->p_d->coord[i], (*coeff_d)[count]);
        count++;
      }

      free(*coeff_d);
    }
    else if (MPType == 1)
    { // intialize strLoc
      strLoc = 0;
      // setup homVarConst_mp
      init_mp2(Left->homVarConst_mp, Right->curr_precision);
      // setup real part
      mpf_set_str(Left->homVarConst_mp->r, &(*codimStr)[strLoc], base);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->homVarConst_mp->i, &(*codimStr)[strLoc], base);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);

      // setup A_mp
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->A_mp->entry[i][j].r, &(*codimStr)[strLoc], base);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->A_mp->entry[i][j].i, &(*codimStr)[strLoc], base);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        }
      // setup H_mp
      for (i = 0; i < Left->H_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->H_mp->coord[i].r, &(*codimStr)[strLoc], base);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->H_mp->coord[i].i, &(*codimStr)[strLoc], base);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      }
      // setup B_mp
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->B_mp->entry[i][j].r, &(*codimStr)[strLoc], base);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->B_mp->entry[i][j].i, &(*codimStr)[strLoc], base);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        }
      // setup p_mp
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->p_mp->coord[i].r, &(*codimStr)[strLoc], base);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->p_mp->coord[i].i, &(*codimStr)[strLoc], base);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      }

      if (freeStr)
        free(*codimStr);
    }
    else // MPType == 2
    { // intialize strLoc
      strLoc = 0;
      // setup homVarConst
      mpq_init(Left->homVarConst_rat[0]);
      mpq_set_str(Left->homVarConst_rat[0], &(*codimStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[0]);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);
      mpq_init(Left->homVarConst_rat[1]);
      mpq_set_str(Left->homVarConst_rat[1], &(*codimStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[1]);
      strLoc += 1 + strlen(&(*codimStr)[strLoc]);

      // setup homVarConst_mp & homVarConst_d
      init_mp2(Left->homVarConst_mp, Right->curr_precision);
      mpf_set_q(Left->homVarConst_mp->r, Left->homVarConst_rat[0]);
      Left->homVarConst_d->r = mpq_get_d(Left->homVarConst_rat[0]);

      mpf_set_q(Left->homVarConst_mp->i, Left->homVarConst_rat[1]);
      Left->homVarConst_d->i = mpq_get_d(Left->homVarConst_rat[1]);

      // setup A_d, A_mp & A_rat
      init_mat_rat(Left->A_rat, Left->A_mp->rows, Left->A_mp->cols);
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->A_rat[i][j][0], &(*codimStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][0]);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->A_rat[i][j][1], &(*codimStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][1]);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);

          // setup A_d & A_mp
          mpf_set_q(Left->A_mp->entry[i][j].r, Left->A_rat[i][j][0]);
          Left->A_d->entry[i][j].r = mpq_get_d(Left->A_rat[i][j][0]);

          mpf_set_q(Left->A_mp->entry[i][j].i, Left->A_rat[i][j][1]);
          Left->A_d->entry[i][j].i = mpq_get_d(Left->A_rat[i][j][1]);
        }

      // setup H_d, H_mp & H_rat
      init_vec_rat(Left->H_rat, Left->H_mp->size);
      for (i = 0; i < Left->H_mp->size; i++)
      { // setup real part
        mpq_set_str(Left->H_rat[i][0], &(*codimStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][0]);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->H_rat[i][1], &(*codimStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][1]);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);

        // setup H_d & H_mp
        mpf_set_q(Left->H_mp->coord[i].r, Left->H_rat[i][0]);
        Left->H_d->coord[i].r = mpq_get_d(Left->H_rat[i][0]);

        mpf_set_q(Left->H_mp->coord[i].i, Left->H_rat[i][1]);
        Left->H_d->coord[i].i = mpq_get_d(Left->H_rat[i][1]);
      }

      // setup B_d, B_mp & B_rat
      init_mat_rat(Left->B_rat, Left->B_mp->rows, Left->B_mp->cols);
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->B_rat[i][j][0], &(*codimStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][0]);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->B_rat[i][j][1], &(*codimStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][1]);
          strLoc += 1 + strlen(&(*codimStr)[strLoc]);

          // setup B_d & B_mp
          mpf_set_q(Left->B_mp->entry[i][j].r, Left->B_rat[i][j][0]);
          Left->B_d->entry[i][j].r = mpq_get_d(Left->B_rat[i][j][0]);

          mpf_set_q(Left->B_mp->entry[i][j].i, Left->B_rat[i][j][1]);
          Left->B_d->entry[i][j].i = mpq_get_d(Left->B_rat[i][j][1]);
        }

      // setup p_d, p_mp & p_rat
      init_vec_rat(Left->p_rat, Left->p_mp->size);
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpq_set_str(Left->p_rat[i][0], &(*codimStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][0]);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->p_rat[i][1], &(*codimStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][1]);
        strLoc += 1 + strlen(&(*codimStr)[strLoc]);

        // setup p_d & p_mp
        mpf_set_q(Left->p_mp->coord[i].r, Left->p_rat[i][0]);
        Left->p_d->coord[i].r = mpq_get_d(Left->p_rat[i][0]);

        mpf_set_q(Left->p_mp->coord[i].i, Left->p_rat[i][1]);
        Left->p_d->coord[i].i = mpq_get_d(Left->p_rat[i][1]);
      }
      // free str
      if (freeStr)
        free(*codimStr);
    }
  }

  return;
}

void cp_cascade_int(void *Out, void *In, int MPType, char **cascadeStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    cascade_t *Right = (cascade_t *)In;
    cascade_t_int *Left = (cascade_t_int *)Out;

    // setup Prog
    cp_prog_t_int(&Left->Prog_int, Right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD_int, &Right->PPD, ppd_type, ppd_size, inType);

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    { // use _d to get the sizes
      Left->A_W_rows = Right->A_d->rows;
      Left->A_W_cols = Right->A_d->cols;
      Left->H_size = Right->H_d->size;
      Left->C_rows = Right->C_d->rows;
      Left->C_cols = Right->C_d->cols;
      Left->R_rows = Right->R_d->rows;
      Left->R_cols = Right->R_d->cols;
      Left->T_size = Right->T_d->size;
      Left->B_rows = Right->B_d->rows;
      Left->B_cols = Right->B_d->cols;
      Left->p_size = Right->p_d->size;
    }
    else
    { // use _mp to get the sizes
      Left->A_W_rows = Right->A_mp->rows;
      Left->A_W_cols = Right->A_mp->cols;
      Left->H_size = Right->H_mp->size;
      Left->C_rows = Right->C_mp->rows;
      Left->C_cols = Right->C_mp->cols;
      Left->R_rows = Right->R_mp->rows;
      Left->R_cols = Right->R_mp->cols;
      Left->T_size = Right->T_mp->size;
      Left->B_rows = Right->B_mp->rows;
      Left->B_cols = Right->B_mp->cols;
      Left->p_size = Right->p_mp->size;
    }

    // setup other values
    Left->curr_precision = Right->curr_precision;
    Left->num_funcs = Right->num_funcs;
    Left->system_rank = Right->system_rank;
    Left->num_codim = Right->num_codim;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup degrees - orig_degrees, new_degrees, P, W, & W_prime
    Left->num_int = 3 * Right->num_funcs + Left->A_W_rows * Left->A_W_cols + Right->system_rank;
    *degrees = (int *)bmalloc(Left->num_int * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    { // copy over orig_degrees, new_degrees & P
      (*degrees)[i] = Right->orig_degrees[i];
      (*degrees)[Right->num_funcs + i] = Right->new_degrees[i];
      (*degrees)[2 * Right->num_funcs + i] = Right->P[i];
    }
    count = 3 * Right->num_funcs;
    for (i = 0; i < Left->A_W_rows; i++)
      for (j = 0; j < Left->A_W_cols; j++)
      { // copy over W
        (*degrees)[count] = Right->W[i][j];
        count++;
      }
    for (i = 0; i < Right->system_rank; i++)
    { // copy over W_prime
      (*degrees)[count] = Right->W_prime[i];
      count++;
    }

    // setup cascadeStr & coeff_d
    if (MPType == 0)
    { // setup num_comp_d
      Left->num_comp_d = 2 + Left->A_W_rows * Left->A_W_cols + Left->H_size + Left->C_rows * Left->C_cols + Left->R_rows * Left->R_cols + Left->T_size + Left->B_rows * Left->B_cols + Left->p_size;

      // set totalLength = 0
      Left->totalLength = 0;

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
      count = 0;
      // setup gamma_d
      set_d((*coeff_d)[count], Right->gamma_d);
      count++;
      // setup homVarConst_d
      set_d((*coeff_d)[count], Right->homVarConst_d);
      count++;
      // setup A_d
      for (i = 0; i < Right->A_d->rows; i++)
        for (j = 0; j < Right->A_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->A_d->entry[i][j]);
          count++;
        }
      // setup H_d
      for (i = 0; i < Right->H_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->H_d->coord[i]);
        count++;
      }
      // setup C_d
      for (i = 0; i < Right->C_d->rows; i++)
        for (j = 0; j < Right->C_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->C_d->entry[i][j]);
          count++;
        }
      // setup R_d
      for (i = 0; i < Right->R_d->rows; i++)
        for (j = 0; j < Right->R_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->R_d->entry[i][j]);
          count++;
        }
      // setup T_d
      for (i = 0; i < Right->T_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->T_d->coord[i]);
        count++;
      }
      // setup B_d
      for (i = 0; i < Right->B_d->rows; i++)
        for (j = 0; j < Right->B_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->B_d->entry[i][j]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Right->p_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->p_d->coord[i]);
        count++;
      }

      // set totalLength = 0
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    { // setup cascadeStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      *cascadeStr = NULL;
      strLoc = 0;

      // setup gamma_mp
      // real part
      strR = mpf_to_str(Right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update cascadeStr
      *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*cascadeStr)[strLoc], strR);
      strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup homVarConst_mp
      // real part
      strR = mpf_to_str(Right->homVarConst_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->homVarConst_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update cascadeStr
      *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*cascadeStr)[strLoc], strR);
      strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_mp
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->A_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->A_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup H_mp
      for (i = 0; i < Right->H_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->H_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->H_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup C_mp
      for (i = 0; i < Right->C_mp->rows; i++)
        for (j = 0; j < Right->C_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->C_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->C_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup R_mp
      for (i = 0; i < Right->R_mp->rows; i++)
        for (j = 0; j < Right->R_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->R_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->R_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup T_mp
      for (i = 0; i < Right->T_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->T_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->T_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup B_mp
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->B_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->B_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_mp
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->p_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->p_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else // MPType == 2
    { // setup cascadeStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      Left->num_comp_d = 0; // no coeff setup in coeff_d

      *cascadeStr = NULL;
      strLoc = 0;

      // setup gamma
      // real part
      strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update cascadeStr
      *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*cascadeStr)[strLoc], strR);
      strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup homVarConst
      // real part
      strR = mpq_get_str(NULL, base, Right->homVarConst_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->homVarConst_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update cascadeStr
      *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*cascadeStr)[strLoc], strR);
      strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_rat
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->A_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->A_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup H_rat
      for (i = 0; i < Right->H_mp->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->H_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->H_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup C_rat
      for (i = 0; i < Right->C_mp->rows; i++)
        for (j = 0; j < Right->C_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->C_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->C_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup R_rat
      for (i = 0; i < Right->R_mp->rows; i++)
        for (j = 0; j < Right->R_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->R_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->R_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup T_mp 
      for (i = 0; i < Right->T_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->T_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->T_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup B_rat
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->B_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->B_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update cascadeStr
          *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*cascadeStr)[strLoc], strR);
          strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_rat
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->p_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->p_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update cascadeStr
        *cascadeStr = (char *)brealloc(*cascadeStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*cascadeStr)[strLoc], strR);
        strcpy(&(*cascadeStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    cascade_t *Left = (cascade_t *)Out;
    cascade_t_int *Right = (cascade_t_int *)In;

    // setup Prog
    Left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    cp_prog_t_int(Left->Prog, &Right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    // initialize evalProg
    initEvalProg(MPType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD, &Right->PPD_int, ppd_type, ppd_size, inType);

    // setup other values
    Left->curr_precision = Right->curr_precision;
    Left->num_funcs = Right->num_funcs;
    Left->system_rank = Right->system_rank;
    Left->num_codim = Right->num_codim;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    { // setup _d sizes
      init_mat_d(Left->A_d, Right->A_W_rows, Right->A_W_cols);
      init_vec_d(Left->H_d, Right->H_size);
      init_mat_d(Left->C_d, Right->C_rows, Right->C_cols);
      init_mat_d(Left->R_d, Right->R_rows, Right->R_cols);
      init_vec_d(Left->T_d, Right->T_size);
      init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
      init_vec_d(Left->p_d, Right->p_size);

      Left->A_d->rows = Right->A_W_rows;
      Left->A_d->cols = Right->A_W_cols;
      Left->H_d->size = Right->H_size;
      Left->C_d->rows = Right->C_rows;
      Left->C_d->cols = Right->C_cols;
      Left->R_d->rows = Right->R_rows;
      Left->R_d->cols = Right->R_cols;
      Left->T_d->size = Right->T_size;
      Left->B_d->rows = Right->B_rows;
      Left->B_d->cols = Right->B_cols;
      Left->p_d->size = Right->p_size;
    }

    if (MPType == 1 || MPType == 2)
    { // setup _mp sizes
      init_mat_mp2(Left->A_mp, Right->A_W_rows, Right->A_W_cols, Right->curr_precision);
      init_vec_mp2(Left->H_mp, Right->H_size, Right->curr_precision);
      init_mat_mp2(Left->C_mp, Right->C_rows, Right->C_cols, Right->curr_precision);
      init_mat_mp2(Left->R_mp, Right->R_rows, Right->R_cols, Right->curr_precision);
      init_vec_mp2(Left->T_mp, Right->T_size, Right->curr_precision);
      init_mat_mp2(Left->B_mp, Right->B_rows, Right->B_cols, Right->curr_precision);
      init_vec_mp2(Left->p_mp, Right->p_size, Right->curr_precision);

      Left->A_mp->rows = Right->A_W_rows;
      Left->A_mp->cols = Right->A_W_cols;
      Left->H_mp->size = Right->H_size;
      Left->C_mp->rows = Right->C_rows;
      Left->C_mp->cols = Right->C_cols;
      Left->R_mp->rows = Right->R_rows;
      Left->R_mp->cols = Right->R_cols;
      Left->T_mp->size = Right->T_size;
      Left->B_mp->rows = Right->B_rows;
      Left->B_mp->cols = Right->B_cols;
      Left->p_mp->size = Right->p_size;
    }

    // setup orig_degrees, new_degrees, P, W, & W_prime
    Left->orig_degrees = (int *)bmalloc(Left->num_funcs * sizeof(int));
    Left->new_degrees = (int *)bmalloc(Left->num_funcs * sizeof(int));
    Left->P = (int *)bmalloc(Left->num_funcs * sizeof(int));
    Left->W_prime = (int *)bmalloc(Left->system_rank * sizeof(int));
    Left->W = (int **)bmalloc(Right->A_W_rows * sizeof(int *));
    for (i = 0; i < Right->A_W_rows; i++)
      Left->W[i] = (int *)bmalloc(Right->A_W_cols * sizeof(int));
    for (i = 0; i < Left->num_funcs; i++)
    { // setup orig_degrees, new_degrees & P
      Left->orig_degrees[i] = (*degrees)[i];
      Left->new_degrees[i] = (*degrees)[Left->num_funcs + i];
      Left->P[i] = (*degrees)[2 * Left->num_funcs + i];
    }
    count = 3 * Left->num_funcs;
    for (i = 0; i < Right->A_W_rows; i++)
      for (j = 0; j < Right->A_W_cols; j++)
      { // setup W
        Left->W[i][j] = (*degrees)[count];
        count++;
      }
    for (i = 0; i < Left->system_rank; i++)
    { // setup W_prime
      Left->W_prime[i] = (*degrees)[count];
      count++;
    }

    // free degrees
    free(*degrees);

    // setup the other structures
    if (MPType == 0)
    { // initialize count
      count = 0;
      // setup gamma_d
      set_d(Left->gamma_d, (*coeff_d)[count]);
      count++;
      // setup homVarConst_d
      set_d(Left->homVarConst_d, (*coeff_d)[count]);
      count++;
      // setup A_d
      for (i = 0; i < Left->A_d->rows; i++)
        for (j = 0; j < Left->A_d->cols; j++)
        {
          set_d(&Left->A_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup H_d
      for (i = 0; i < Left->H_d->size; i++)
      {
        set_d(&Left->H_d->coord[i], (*coeff_d)[count]);
        count++;
      }
      // setup C_d
      for (i = 0; i < Left->C_d->rows; i++)
        for (j = 0; j < Left->C_d->cols; j++)
        {
          set_d(&Left->C_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup R_d
      for (i = 0; i < Left->R_d->rows; i++)
        for (j = 0; j < Left->R_d->cols; j++)
        {
          set_d(&Left->R_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup T_d
      for (i = 0; i < Left->T_d->size; i++)
      {
        set_d(&Left->T_d->coord[i], (*coeff_d)[count]);
        count++;
      }
      // setup B_d
      for (i = 0; i < Left->B_d->rows; i++)
        for (j = 0; j < Left->B_d->cols; j++)
        {
          set_d(&Left->B_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Left->p_d->size; i++)
      {
        set_d(&Left->p_d->coord[i], (*coeff_d)[count]);
        count++;
      }

      free(*coeff_d);
    }
    else if (MPType == 1)
    { // intialize strLoc
      strLoc = 0;
      // setup gamma_mp
      init_mp2(Left->gamma_mp, Right->curr_precision);
      // setup real part
      mpf_set_str(Left->gamma_mp->r, &(*cascadeStr)[strLoc], base);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->gamma_mp->i, &(*cascadeStr)[strLoc], base);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

      // setup homVarConst_mp
      init_mp2(Left->homVarConst_mp, Right->curr_precision);
      // setup real part
      mpf_set_str(Left->homVarConst_mp->r, &(*cascadeStr)[strLoc], base);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->homVarConst_mp->i, &(*cascadeStr)[strLoc], base);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

      // setup A_mp
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->A_mp->entry[i][j].r, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->A_mp->entry[i][j].i, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        }
      // setup H_mp
      for (i = 0; i < Left->H_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->H_mp->coord[i].r, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->H_mp->coord[i].i, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      }
      // setup C_mp
      for (i = 0; i < Left->C_mp->rows; i++)
        for (j = 0; j < Left->C_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->C_mp->entry[i][j].r, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->C_mp->entry[i][j].i, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        }
      // setup R_mp
      for (i = 0; i < Left->R_mp->rows; i++)
        for (j = 0; j < Left->R_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->R_mp->entry[i][j].r, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->R_mp->entry[i][j].i, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        }
      // setup T_mp
      for (i = 0; i < Left->T_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->T_mp->coord[i].r, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->T_mp->coord[i].i, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      }
      // setup B_mp
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->B_mp->entry[i][j].r, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->B_mp->entry[i][j].i, &(*cascadeStr)[strLoc], base);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        }
      // setup p_mp
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->p_mp->coord[i].r, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->p_mp->coord[i].i, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      }

      // free string
      if (freeStr)
        free(*cascadeStr);
    }
    else // MPType == 2
    { // intialize strLoc
      strLoc = 0;
      // setup gamma
      mpq_init(Left->gamma_rat[0]);
      mpq_set_str(Left->gamma_rat[0], &(*cascadeStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      mpq_init(Left->gamma_rat[1]);
      mpq_set_str(Left->gamma_rat[1], &(*cascadeStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[1]);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

      // setup gamma_mp & gamma_d
      init_mp2(Left->gamma_mp, Right->curr_precision);
      mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
      Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);

      mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
      Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);

      // setup homVarConst
      mpq_init(Left->homVarConst_rat[0]);
      mpq_set_str(Left->homVarConst_rat[0], &(*cascadeStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[0]);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
      mpq_init(Left->homVarConst_rat[1]);
      mpq_set_str(Left->homVarConst_rat[1], &(*cascadeStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[1]);
      strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

      // setup homVarConst_mp & homVarConst_d
      init_mp2(Left->homVarConst_mp, Right->curr_precision);
      mpf_set_q(Left->homVarConst_mp->r, Left->homVarConst_rat[0]);
      Left->homVarConst_d->r = mpq_get_d(Left->homVarConst_rat[0]);

      mpf_set_q(Left->homVarConst_mp->i, Left->homVarConst_rat[1]);
      Left->homVarConst_d->i = mpq_get_d(Left->homVarConst_rat[1]);

      // setup A_d, A_mp & A_rat
      init_mat_rat(Left->A_rat, Left->A_mp->rows, Left->A_mp->cols);
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->A_rat[i][j][0], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][0]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->A_rat[i][j][1], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][1]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

          // setup A_d & A_mp
          mpf_set_q(Left->A_mp->entry[i][j].r, Left->A_rat[i][j][0]);
          Left->A_d->entry[i][j].r = mpq_get_d(Left->A_rat[i][j][0]);

          mpf_set_q(Left->A_mp->entry[i][j].i, Left->A_rat[i][j][1]);
          Left->A_d->entry[i][j].i = mpq_get_d(Left->A_rat[i][j][1]);
        }
      // setup H_d, H_mp & H_rat
      init_vec_rat(Left->H_rat, Left->H_mp->size);
      for (i = 0; i < Left->H_mp->size; i++)
      { // setup real part
        mpq_set_str(Left->H_rat[i][0], &(*cascadeStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][0]);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->H_rat[i][1], &(*cascadeStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][1]);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

        // setup H_d & H_mp
        mpf_set_q(Left->H_mp->coord[i].r, Left->H_rat[i][0]);
        Left->H_d->coord[i].r = mpq_get_d(Left->H_rat[i][0]);

        mpf_set_q(Left->H_mp->coord[i].i, Left->H_rat[i][1]);
        Left->H_d->coord[i].i = mpq_get_d(Left->H_rat[i][1]);
      }
      // setup C_d, C_mp & C_rat
      init_mat_rat(Left->C_rat, Left->C_mp->rows, Left->C_mp->cols);
      for (i = 0; i < Left->C_mp->rows; i++)
        for (j = 0; j < Left->C_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->C_rat[i][j][0], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->C_rat[i][j][0]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->C_rat[i][j][1], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->C_rat[i][j][1]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

          // setup C_d & C_mp
          mpf_set_q(Left->C_mp->entry[i][j].r, Left->C_rat[i][j][0]);
          Left->C_d->entry[i][j].r = mpq_get_d(Left->C_rat[i][j][0]);

          mpf_set_q(Left->C_mp->entry[i][j].i, Left->C_rat[i][j][1]);
          Left->C_d->entry[i][j].i = mpq_get_d(Left->C_rat[i][j][1]);
        }
      // setup R_d, R_mp & R_rat
      init_mat_rat(Left->R_rat, Left->R_mp->rows, Left->R_mp->cols);
      for (i = 0; i < Left->R_mp->rows; i++)
        for (j = 0; j < Left->R_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->R_rat[i][j][0], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->R_rat[i][j][0]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->R_rat[i][j][1], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->R_rat[i][j][1]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

          // setup R_d & R_mp
          mpf_set_q(Left->R_mp->entry[i][j].r, Left->R_rat[i][j][0]);
          Left->R_d->entry[i][j].r = mpq_get_d(Left->R_rat[i][j][0]);

          mpf_set_q(Left->R_mp->entry[i][j].i, Left->R_rat[i][j][1]);
          Left->R_d->entry[i][j].i = mpq_get_d(Left->R_rat[i][j][1]);
        }
      // setup T_d & T_mp
      for (i = 0; i < Left->T_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->T_mp->coord[i].r, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->T_mp->coord[i].i, &(*cascadeStr)[strLoc], base);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

        // setup T_d
        mp_to_d(&Left->T_d->coord[i], &Left->T_mp->coord[i]);
      }
      // setup B_d, B_mp & B_rat
      init_mat_rat(Left->B_rat, Left->B_mp->rows, Left->B_mp->cols);
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->B_rat[i][j][0], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][0]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->B_rat[i][j][1], &(*cascadeStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][1]);
          strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

          // setup B_d & B_mp
          mpf_set_q(Left->B_mp->entry[i][j].r, Left->B_rat[i][j][0]);
          Left->B_d->entry[i][j].r = mpq_get_d(Left->B_rat[i][j][0]);

          mpf_set_q(Left->B_mp->entry[i][j].i, Left->B_rat[i][j][1]);
          Left->B_d->entry[i][j].i = mpq_get_d(Left->B_rat[i][j][1]);
        }
      // setup p_d, p_mp & p_rat
      init_vec_rat(Left->p_rat, Left->p_mp->size);
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpq_set_str(Left->p_rat[i][0], &(*cascadeStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][0]);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->p_rat[i][1], &(*cascadeStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][1]);
        strLoc += 1 + strlen(&(*cascadeStr)[strLoc]);

        // setup p_d & p_mp
        mpf_set_q(Left->p_mp->coord[i].r, Left->p_rat[i][0]);
        Left->p_d->coord[i].r = mpq_get_d(Left->p_rat[i][0]);

        mpf_set_q(Left->p_mp->coord[i].i, Left->p_rat[i][1]);
        Left->p_d->coord[i].i = mpq_get_d(Left->p_rat[i][1]);
      }

      // free string
      if (freeStr)
        free(*cascadeStr);
    }
  }

  return;
}

void cp_cascadeCodim_int(void *Out, void *In, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  if (inType == 0)
  { // _t to _t_int so that it can be sent
    cascadeCodim_t *Right = (cascadeCodim_t *)In;
    cascadeCodim_t_int *Left = (cascadeCodim_t_int *)Out;

    Left->codim = Right->codim;
  }
  else
  { // _t_int to _t so that it can be used by the workers
    cascadeCodim_t *Left = (cascadeCodim_t *)Out;
    cascadeCodim_t_int *Right = (cascadeCodim_t_int *)In;

    // setup values
    Left->codim = Right->codim;
    Left->num_paths = Left->num_superset = Left->num_nonsing = Left->num_sing = Left->num_nonsolns = Left->num_inf = Left->num_bad = 0;

    // NULL out other pointers
    Left->startPts_d = NULL;
    Left->startPts_mp = NULL;
    Left->endPts_d = NULL;
    Left->endPts_mp = NULL;
    Left->endPts_amp = NULL;
    Left->endPt_types = NULL;
  }

  return;
}

void cp_slice_moving_int(void *Out, char **sliceStr, char **progStr, int **prog_inst, int **prog_gp_sizes, comp_d **coeff_d, void *In, int sendProg, int MPType, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is slice             *
*                otherwise - Out is slice, In is _t_int         *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, sizeR = 0, sizeI = 0, strLoc = 0, base = 10;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    membership_slice_moving_t *Right = (membership_slice_moving_t *)In;
    membership_slice_moving_t_int *Left = (membership_slice_moving_t_int *)Out;

    // setup Prog
    if (sendProg)
    { // setup the structures
      cp_prog_t_int(&Left->Prog_int, Right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    }
    else
    { // NULL out the structures
      *prog_inst = NULL;
      *prog_gp_sizes = NULL;
      *progStr = NULL;
      Left->Prog_int.totalLength = 0;
    }

    // copy other data
    Left->curr_codim = Right->curr_codim;
    Left->orig_variables = Right->orig_variables;
    Left->curr_precision = Right->curr_precision;

    Left->startSliceVec_init = Right->startSliceVec_init;
    Left->targetSliceVec_init = Right->targetSliceVec_init;
    Left->K_rows = Right->K_rows;
    Left->K_cols = Right->K_cols;

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    {
      Left->A_rows = Right->A_d->rows;
      Left->A_cols = Right->A_d->cols;
      Left->B_rows = Right->B_d->rows;
      Left->B_cols = Right->B_d->cols;
      Left->p_size = Right->p_d->size;
      Left->startSliceVec_size = Left->startSliceVec_init ? Right->startSliceVec_d->size : 0;
      Left->targetSliceVec_size = Left->targetSliceVec_init ? Right->targetSliceVec_d->size : 0;
    }
    else
    {
      Left->A_rows = Right->A_mp->rows;
      Left->A_cols = Right->A_mp->cols;
      Left->B_rows = Right->B_mp->rows;
      Left->B_cols = Right->B_mp->cols;
      Left->p_size = Right->p_mp->size;
      Left->startSliceVec_size = Left->startSliceVec_init ? Right->startSliceVec_mp->size : 0;
      Left->targetSliceVec_size = Left->targetSliceVec_init ? Right->targetSliceVec_mp->size : 0;
    }

    // setup the first part of sliceStr using K_rat
    for (i = 0; i < Right->K_rows; i++)
      for (j = 0; j < Right->K_cols; j++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->K_rat[i][j][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->K_rat[i][j][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update sliceStr
        *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*sliceStr)[strLoc], strR);
        strcpy(&(*sliceStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

    // setup the rest of sliceStr & coeff_d
    if (MPType == 0)
    { // setup num_comp_d
      Left->num_comp_d = 2 + Left->A_rows * Left->A_cols + Left->B_rows * Left->B_cols + Left->p_size + Left->startSliceVec_size + Left->targetSliceVec_size;

      // setup totalLength
      Left->totalLength = strLoc;

      // setup coeff_d
      *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
      count = 0;

      // setup gamma_d
      set_d((*coeff_d)[count], Right->gamma_d);
      count++;
      // setup A_d
      for (i = 0; i < Right->A_d->rows; i++)
        for (j = 0; j < Right->A_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->A_d->entry[i][j]);
          count++;
        }
      // setup B_d
      for (i = 0; i < Right->B_d->rows; i++)
        for (j = 0; j < Right->B_d->cols; j++)
        {
          set_d((*coeff_d)[count], &Right->B_d->entry[i][j]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Right->p_d->size; i++)
      {
        set_d((*coeff_d)[count], &Right->p_d->coord[i]);
        count++;
      }
      // setup startSliceVec_d
      if (Right->startSliceVec_init)
      {
        for (i = 0; i < Right->startSliceVec_d->size; i++)
        {
          set_d((*coeff_d)[count], &Right->startSliceVec_d->coord[i]);
          count++;
        }
      }
      // setup targetSliceVec_d
      if (Right->targetSliceVec_init)
      {
        for (i = 0; i < Right->targetSliceVec_d->size; i++)
        {
          set_d((*coeff_d)[count], &Right->targetSliceVec_d->coord[i]);
          count++;
        }
      }
    }
    else if (MPType == 1)
    { // setup num_comp_d
      Left->num_comp_d = 0;

      // setup gamma_mp
      // real part
      strR = mpf_to_str(Right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update sliceStr
      *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*sliceStr)[strLoc], strR);
      strcpy(&(*sliceStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_mp
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->A_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->A_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup B_mp
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpf_to_str(Right->B_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->B_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_mp
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpf_to_str(Right->p_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->p_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update sliceStr
        *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*sliceStr)[strLoc], strR);
        strcpy(&(*sliceStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup startSliceVec_mp
      if (Right->startSliceVec_init)
      {
        for (i = 0; i < Right->startSliceVec_mp->size; i++)
        { // real part
          strR = mpf_to_str(Right->startSliceVec_mp->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->startSliceVec_mp->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }
      // setup targetSliceVec_mp
      if (Right->targetSliceVec_init)
      {
        for (i = 0; i < Right->targetSliceVec_mp->size; i++)
        { // real part
          strR = mpf_to_str(Right->targetSliceVec_mp->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->targetSliceVec_mp->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
    else // MPType == 2
    { // setup num_comp_d
      Left->num_comp_d = 0;

      // setup gamma
      // real part
      strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update sliceStr
      *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*sliceStr)[strLoc], strR);
      strcpy(&(*sliceStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // setup A_rat
      for (i = 0; i < Right->A_mp->rows; i++)
        for (j = 0; j < Right->A_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->A_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->A_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup B_rat
      for (i = 0; i < Right->B_mp->rows; i++)
        for (j = 0; j < Right->B_mp->cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->B_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->B_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // setup p_rat
      for (i = 0; i < Right->p_mp->size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->p_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->p_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update sliceStr
        *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*sliceStr)[strLoc], strR);
        strcpy(&(*sliceStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // setup startSliceVec_rat
      if (Right->startSliceVec_init)
      {
        for (i = 0; i < Right->startSliceVec_mp->size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->startSliceVec_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->startSliceVec_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }
      // setup startSliceVec_rat
      if (Right->targetSliceVec_init)
      {
        for (i = 0; i < Right->targetSliceVec_mp->size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->targetSliceVec_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->targetSliceVec_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update sliceStr
          *sliceStr = (char *)brealloc(*sliceStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*sliceStr)[strLoc], strR);
          strcpy(&(*sliceStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      }

      // setup totalLength
      Left->totalLength = strLoc;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    membership_slice_moving_t *Left = (membership_slice_moving_t *)Out;
    membership_slice_moving_t_int *Right = (membership_slice_moving_t_int *)In;

    // setup Prog
    if (sendProg)
    { // setup Prog
      Left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
      cp_prog_t_int(Left->Prog, &Right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    }
    else
    { // NULL it out
      Left->Prog = NULL;
    }
    // intialize evalProg
    initEvalProg(MPType);

    // setup other values
    Left->curr_codim = Right->curr_codim;
    Left->orig_variables = Right->orig_variables;
    Left->curr_precision = Right->curr_precision;

    Left->startSliceVec_init = Right->startSliceVec_init;
    Left->targetSliceVec_init = Right->targetSliceVec_init;
    Left->K_rows = Right->K_rows;
    Left->K_cols = Right->K_cols;

    // setup K
    init_mat_rat(Left->K_rat, Left->K_rows, Left->K_cols);
    for (i = 0; i < Left->K_rows; i++)
      for (j = 0; j < Left->K_cols; j++)
      { // setup real part
        mpq_set_str(Left->K_rat[i][j][0], &(*sliceStr)[strLoc], base);
        mpq_canonicalize(Left->K_rat[i][j][0]);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->K_rat[i][j][1], &(*sliceStr)[strLoc], base);
        mpq_canonicalize(Left->K_rat[i][j][1]);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
      }

    // setup the sizes
    if (MPType == 0 || MPType == 2)
    { // setup _d sizes
      init_mat_d(Left->A_d, Right->A_rows, Right->A_cols);
      init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
      init_vec_d(Left->p_d, Right->p_size);
      if (Left->startSliceVec_init)
        init_vec_d(Left->startSliceVec_d, Right->startSliceVec_size);
      if (Left->targetSliceVec_init)
        init_vec_d(Left->targetSliceVec_d, Right->targetSliceVec_size);

      Left->A_d->rows = Right->A_rows;
      Left->A_d->cols = Right->A_cols;
      Left->B_d->rows = Right->B_rows;
      Left->B_d->cols = Right->B_cols;
      Left->p_d->size = Right->p_size;
      Left->startSliceVec_d->size = Right->startSliceVec_init ? Right->startSliceVec_size : 0;
      Left->targetSliceVec_d->size = Right->targetSliceVec_init ? Right->targetSliceVec_size : 0;
    }

    if (MPType == 1 || MPType == 2)
    { // setup _mp sizes
      init_mat_mp2(Left->A_mp, Right->A_rows, Right->A_cols, Right->curr_precision);
      init_mat_mp2(Left->B_mp, Right->B_rows, Right->B_cols, Right->curr_precision);
      init_vec_mp2(Left->p_mp, Right->p_size, Right->curr_precision);
      if (Left->startSliceVec_init)
        init_vec_mp2(Left->startSliceVec_mp, Right->startSliceVec_size, Right->curr_precision);
      if (Left->targetSliceVec_init)
        init_vec_mp2(Left->targetSliceVec_mp, Right->targetSliceVec_size, Right->curr_precision);

      Left->A_mp->rows = Right->A_rows;
      Left->A_mp->cols = Right->A_cols;
      Left->B_mp->rows = Right->B_rows;
      Left->B_mp->cols = Right->B_cols;
      Left->p_mp->size = Right->p_size;
      Left->startSliceVec_mp->size = Right->startSliceVec_init ? Right->startSliceVec_size : 0;
      Left->targetSliceVec_mp->size = Right->targetSliceVec_init ? Right->targetSliceVec_size : 0;
    }

    if (MPType == 2)
    { // setup _rat 
      init_mat_rat(Left->A_rat, Right->A_rows, Right->A_cols);
      init_mat_rat(Left->B_rat, Right->B_rows, Right->B_cols);
      init_vec_rat(Left->p_rat, Right->p_size);
      if (Left->startSliceVec_init)
        init_vec_rat(Left->startSliceVec_rat, Right->startSliceVec_size);
      if (Left->targetSliceVec_init)
        init_vec_rat(Left->targetSliceVec_rat, Right->targetSliceVec_size);
    }

    // setup the other structures
    if (MPType == 0)
    { // setup gamma_d
      set_d(Left->gamma_d, (*coeff_d)[count]);
      count++;
      // setup A_d
      for (i = 0; i < Left->A_d->rows; i++)
        for (j = 0; j < Left->A_d->cols; j++)
        {
          set_d(&Left->A_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup B_d
      for (i = 0; i < Left->B_d->rows; i++)
        for (j = 0; j < Left->B_d->cols; j++)
        {
          set_d(&Left->B_d->entry[i][j], (*coeff_d)[count]);
          count++;
        }
      // setup p_d
      for (i = 0; i < Left->p_d->size; i++)
      {
        set_d(&Left->p_d->coord[i], (*coeff_d)[count]);
        count++;
      }
      // setup startSliceVec_d
      if (Left->startSliceVec_init)
      {
        for (i = 0; i < Left->startSliceVec_d->size; i++)
        {
          set_d(&Left->startSliceVec_d->coord[i], (*coeff_d)[count]);
          count++;
        }
      }
      // setup targetSliceVec_d
      if (Left->targetSliceVec_init)
      {
        for (i = 0; i < Left->targetSliceVec_d->size; i++)
        {
          set_d(&Left->targetSliceVec_d->coord[i], (*coeff_d)[count]);
          count++;
        }
      }

      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup gamma_mp
      init_mp2(Left->gamma_mp, Right->curr_precision);
      // setup real part
      mpf_set_str(Left->gamma_mp->r, &(*sliceStr)[strLoc], base);
      strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->gamma_mp->i, &(*sliceStr)[strLoc], base);
      strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

      // setup A_mp
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->A_mp->entry[i][j].r, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->A_mp->entry[i][j].i, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        }
      // setup B_mp
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpf_set_str(Left->B_mp->entry[i][j].r, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->B_mp->entry[i][j].i, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        }
      // setup p_mp
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpf_set_str(Left->p_mp->coord[i].r, &(*sliceStr)[strLoc], base);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->p_mp->coord[i].i, &(*sliceStr)[strLoc], base);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
      }
      // setup startSliceVec_mp
      if (Left->startSliceVec_init)
      {
        for (i = 0; i < Left->startSliceVec_mp->size; i++)
        { // setup real part
          mpf_set_str(Left->startSliceVec_mp->coord[i].r, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->startSliceVec_mp->coord[i].i, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        }
      }
      // setup targetSliceVec_mp
      if (Left->targetSliceVec_init)
      {
        for (i = 0; i < Left->targetSliceVec_mp->size; i++)
        { // setup real part
          mpf_set_str(Left->targetSliceVec_mp->coord[i].r, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->targetSliceVec_mp->coord[i].i, &(*sliceStr)[strLoc], base);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        }
      }

      // free string
      if (freeStr)
        free(*sliceStr);
    }
    else // MPType == 2
    { // setup gamma
      mpq_init(Left->gamma_rat[0]);
      mpq_set_str(Left->gamma_rat[0], &(*sliceStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
      mpq_init(Left->gamma_rat[1]);
      mpq_set_str(Left->gamma_rat[1], &(*sliceStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[1]);
      strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

      // setup gamma_mp & gamma_d
      init_mp2(Left->gamma_mp, Right->curr_precision);
      mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
      Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);

      mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
      Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);

      // setup A_d, A_mp & A_rat
      for (i = 0; i < Left->A_mp->rows; i++)
        for (j = 0; j < Left->A_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->A_rat[i][j][0], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][0]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->A_rat[i][j][1], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][1]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

          // setup A_d & A_mp
          mpf_set_q(Left->A_mp->entry[i][j].r, Left->A_rat[i][j][0]);
          Left->A_d->entry[i][j].r = mpq_get_d(Left->A_rat[i][j][0]);

          mpf_set_q(Left->A_mp->entry[i][j].i, Left->A_rat[i][j][1]);
          Left->A_d->entry[i][j].i = mpq_get_d(Left->A_rat[i][j][1]);
        }
      // setup B_d, B_mp & B_rat
      for (i = 0; i < Left->B_mp->rows; i++)
        for (j = 0; j < Left->B_mp->cols; j++)
        { // setup real part
          mpq_set_str(Left->B_rat[i][j][0], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][0]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->B_rat[i][j][1], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][1]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

          // setup B_d & B_mp
          mpf_set_q(Left->B_mp->entry[i][j].r, Left->B_rat[i][j][0]);
          Left->B_d->entry[i][j].r = mpq_get_d(Left->B_rat[i][j][0]);

          mpf_set_q(Left->B_mp->entry[i][j].i, Left->B_rat[i][j][1]);
          Left->B_d->entry[i][j].i = mpq_get_d(Left->B_rat[i][j][1]);
        }
      // setup p_d, p_mp & p_rat
      for (i = 0; i < Left->p_mp->size; i++)
      { // setup real part
        mpq_set_str(Left->p_rat[i][0], &(*sliceStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][0]);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->p_rat[i][1], &(*sliceStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][1]);
        strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

        // setup p_d & p_mp
        mpf_set_q(Left->p_mp->coord[i].r, Left->p_rat[i][0]);
        Left->p_d->coord[i].r = mpq_get_d(Left->p_rat[i][0]);

        mpf_set_q(Left->p_mp->coord[i].i, Left->p_rat[i][1]);
        Left->p_d->coord[i].i = mpq_get_d(Left->p_rat[i][1]);
      }

      // setup startSliceVec_d,_mp,_rat
      if (Left->startSliceVec_init)
      {
        for (i = 0; i < Left->startSliceVec_mp->size; i++)
        { // setup real part
          mpq_set_str(Left->startSliceVec_rat[i][0], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->startSliceVec_rat[i][0]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->startSliceVec_rat[i][1], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->startSliceVec_rat[i][1]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

          // setup startSliceVec_d & startSliceVec_mp
          mpf_set_q(Left->startSliceVec_mp->coord[i].r, Left->startSliceVec_rat[i][0]);
          Left->startSliceVec_d->coord[i].r = mpq_get_d(Left->startSliceVec_rat[i][0]);

          mpf_set_q(Left->startSliceVec_mp->coord[i].i, Left->startSliceVec_rat[i][1]);
          Left->startSliceVec_d->coord[i].i = mpq_get_d(Left->startSliceVec_rat[i][1]);
        }
      }
      // setup targetSliceVec_d,_mp,_rat
      if (Left->targetSliceVec_init)
      {
        for (i = 0; i < Left->targetSliceVec_mp->size; i++)
        { // setup real part
          mpq_set_str(Left->targetSliceVec_rat[i][0], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->targetSliceVec_rat[i][0]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->targetSliceVec_rat[i][1], &(*sliceStr)[strLoc], base);
          mpq_canonicalize(Left->targetSliceVec_rat[i][1]);
          strLoc += 1 + strlen(&(*sliceStr)[strLoc]);

          // setup targetSliceVec_d & targetSliceVec_mp
          mpf_set_q(Left->targetSliceVec_mp->coord[i].r, Left->targetSliceVec_rat[i][0]);
          Left->targetSliceVec_d->coord[i].r = mpq_get_d(Left->targetSliceVec_rat[i][0]);

          mpf_set_q(Left->targetSliceVec_mp->coord[i].i, Left->targetSliceVec_rat[i][1]);
          Left->targetSliceVec_d->coord[i].i = mpq_get_d(Left->targetSliceVec_rat[i][1]);
        }
      }

      // free string
      if (freeStr)
        free(*sliceStr);
    }
  }

  return;
}

void cp_endpoint_data_d_int(void *Out, void *In, comp_d **coeff, int freeCoeff, int initPoint, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _d_int,In is _d                 *
*                otherwise - Out is _d, In is _d_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, count = 0, size = 0;

  if (inType == 0)
  { // _d to _d_int so that it can be sent
    endpoint_data_d *Right = (endpoint_data_d *)In;
    endpoint_data_d_int *Left = (endpoint_data_d_int *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // compute number of comp_d needed
    Left->num_comp_d = 1 + Right->endPt->size + Right->last_approx->size;

    // setup coeff
    *coeff = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
    set_d((*coeff)[count], Right->finalT);
    count++;

    size = Left->pt_size = Right->endPt->size;
    for (i = 0; i < size; i++)
    {
      set_d((*coeff)[count], &Right->endPt->coord[i]);
      count++;
    }

    size = Left->last_approx_size = Right->last_approx->size;
    for (i = 0; i < size; i++)
    {
      set_d((*coeff)[count], &Right->last_approx->coord[i]);
      count++;
    }
  }
  else
  { // _d_int to _d so that it can be used
    endpoint_data_d_int *Right = (endpoint_data_d_int *)In;
    endpoint_data_d *Left = (endpoint_data_d *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // initialize
    if (initPoint)
    {
      init_endpoint_data_d(Left);
    }

    change_size_point_d(Left->endPt, Right->pt_size);
    Left->endPt->size = Right->pt_size;
    change_size_point_d(Left->last_approx, Right->last_approx_size);
    Left->last_approx->size = Right->last_approx_size;

    // setup values
    set_d(Left->finalT, (*coeff)[count]);
    count++;
    size = Right->pt_size;
    for (i = 0; i < size; i++)
    {
      set_d(&Left->endPt->coord[i], (*coeff)[count]);
      count++;
    }
    size = Right->last_approx_size;
    for (i = 0; i < size; i++)
    {
      set_d(&Left->last_approx->coord[i], (*coeff)[count]);
      count++;
    }

    // clear coeff
    if (freeCoeff)
      free(*coeff);
  }

  return;
}

void cp_endpoint_data_mp_int(void *Out, void *In, char **epStr, int freeStr, int initPoint, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _mp_int,In is _mp               *
*                otherwise - Out is _mp, In is _mp_int          *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, size = 0, strLoc = 0, base = 10, sizeR = 0, sizeI = 0;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    endpoint_data_mp *Right = (endpoint_data_mp *)In;
    endpoint_data_mp_int *Left = (endpoint_data_mp_int *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // setup epStr
    *epStr = NULL;

    // copy finalT
    // real part
    strR = mpf_to_str(Right->finalT->r, base);
    sizeR = strlen(strR) + 1; // +1 for '\0'
    // imag part
    strI = mpf_to_str(Right->finalT->i, base);
    sizeI = strlen(strI) + 1; // +1 for '\0'
    // update epStr
    *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
    strcpy(&(*epStr)[strLoc], strR);
    strcpy(&(*epStr)[strLoc + sizeR], strI);
    // update strLoc
    strLoc += sizeR + sizeI;
    // free strR, strI
    free(strR);
    free(strI);

    // copy endPt
    size = Left->pt_size = Right->endPt->size;
    for (i = 0; i < size; i++)
    { // real part
      strR = mpf_to_str(Right->endPt->coord[i].r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->endPt->coord[i].i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update epStr
      *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*epStr)[strLoc], strR);
      strcpy(&(*epStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);
    }

    size = Left->last_approx_size = Right->last_approx->size;
    // copy last_approx
    for (i = 0; i < size; i++)
    { // real part
      strR = mpf_to_str(Right->last_approx->coord[i].r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->last_approx->coord[i].i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update epStr
      *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*epStr)[strLoc], strR);
      strcpy(&(*epStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);
    }

    // seutp totalLength
    Left->totalLength = strLoc;
  }
  else
  { // _mp_int to _mp so that it can be used
    endpoint_data_mp_int *Right = (endpoint_data_mp_int *)In;
    endpoint_data_mp *Left = (endpoint_data_mp *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // initialize
    if (initPoint)
    {
      init_endpoint_data_mp(Left);
    }

    change_size_point_mp(Left->endPt, Right->pt_size);
    Left->endPt->size = Right->pt_size;
    change_size_point_mp(Left->last_approx, Right->last_approx_size);
    Left->last_approx->size = Right->last_approx_size;

    // setup finalT
    // setup real part
    mpf_set_str(Left->finalT->r, &(*epStr)[strLoc], base);
    strLoc += 1 + strlen(&(*epStr)[strLoc]);
    // setup imag part
    mpf_set_str(Left->finalT->i, &(*epStr)[strLoc], base);
    strLoc += 1 + strlen(&(*epStr)[strLoc]);

    // setup endPt
    size = Right->pt_size;
    for (i = 0; i < size; i++)
    { // setup real part
      mpf_set_str(Left->endPt->coord[i].r, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->endPt->coord[i].i, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);
    }
    // setup last_approx
    size = Right->last_approx_size;
    for (i = 0; i < size; i++)
    { // setup real part
      mpf_set_str(Left->last_approx->coord[i].r, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->last_approx->coord[i].i, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);
    }

    // clear str
    if (freeStr)
      free(*epStr);
  }

  return;
}

void cp_endpoint_data_amp_int(void *Out, void *In, comp_d **coeff, char **epStr, int freeData, int initPoint, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _amp_int,In is _amp             *
*                otherwise - Out is _amp, In is _amp_int        *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, count = 0, size = 0, strLoc = 0, base = 10, sizeR = 0, sizeI = 0;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // _mp to _mp_int so that it can be sent
    endpoint_data_amp *Right = (endpoint_data_amp *)In;
    endpoint_data_amp_int *Left = (endpoint_data_amp_int *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // store the precision
    Left->curr_prec = Right->curr_prec;
    Left->last_approx_prec = Right->last_approx_prec;
 
    // setup coeff & epStr
    *coeff = NULL;
    *epStr = NULL;

    // copy finalT & endPt
    if (Right->curr_prec < 64)
    { // setup using _d
      *coeff = (comp_d *)brealloc(*coeff, (count + Right->endPt_d->size + 1) * sizeof(comp_d));

      set_d((*coeff)[count], Right->finalT_d);
      count++;

      size = Left->pt_size = Right->endPt_d->size;
      for (i = 0; i < size; i++)
      {
        set_d((*coeff)[count], &Right->endPt_d->coord[i]);
        count++;
      }
    }
    else
    { // setup using _mp
      // real part
      strR = mpf_to_str(Right->finalT_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->finalT_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'
      // update epStr
      *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*epStr)[strLoc], strR);
      strcpy(&(*epStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;
      // free strR, strI
      free(strR);
      free(strI);

      // copy endPt
      size = Left->pt_size = Right->endPt_mp->size;
      for (i = 0; i < size; i++)
      { // real part
        strR = mpf_to_str(Right->endPt_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->endPt_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update epStr
        *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*epStr)[strLoc], strR);
        strcpy(&(*epStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;
  
        // free strR, strI
        free(strR);
        free(strI);
      }
    }
    
    // copy finalT & endPt
    // copy last_approx
    if (Right->last_approx_prec < 64)
    { // setup using _d
      *coeff = (comp_d *)brealloc(*coeff, (count + Right->last_approx_d->size) * sizeof(comp_d));

      size = Left->last_approx_size = Right->last_approx_d->size;
      for (i = 0; i < size; i++)
      {
        set_d((*coeff)[count], &Right->last_approx_d->coord[i]);
        count++;
      }
    }
    else
    { // setup using _mp
      size = Left->last_approx_size = Right->last_approx_mp->size;
      for (i = 0; i < size; i++)
      { // real part
        strR = mpf_to_str(Right->last_approx_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->last_approx_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update epStr
        *epStr = (char *)brealloc(*epStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*epStr)[strLoc], strR);
        strcpy(&(*epStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
    }

    // seutp totalLength & num_comp_d
    Left->totalLength = strLoc;
    Left->num_comp_d = count;
  }
  else
  { // _amp_int to _amp so that it can be used
    endpoint_data_amp_int *Right = (endpoint_data_amp_int *)In;
    endpoint_data_amp *Left = (endpoint_data_amp *)Out;

    // setup other vaules
    Left->cond_num = Right->cond_num;
    Left->corank = Right->corank;
    Left->smallest_nonzero_SV = Right->smallest_nonzero_SV;
    Left->largest_zero_SV = Right->largest_zero_SV;
    Left->retVal = Right->retVal;

    // store the precision
    Left->curr_prec = Right->curr_prec;
    Left->last_approx_prec = Right->last_approx_prec;

    // initialize
    if (initPoint)
    {
      init_endpoint_data_amp(Left, Left->curr_prec, Left->last_approx_prec);
    }

    if (Left->curr_prec < 64)
    {
      change_size_point_d(Left->endPt_d, Right->pt_size);
      Left->endPt_d->size = Right->pt_size;
    }
    else
    {
      change_prec_mp(Left->finalT_mp, Left->curr_prec);
      change_prec_point_mp(Left->endPt_mp, Left->curr_prec);
      change_size_point_mp(Left->endPt_mp, Right->pt_size);
      Left->endPt_mp->size = Right->pt_size;
    }

    if (Left->last_approx_prec < 64)
    {
      change_size_point_d(Left->last_approx_d, Right->last_approx_size);
      Left->last_approx_d->size = Right->last_approx_size;
    }
    else
    {
      change_prec_point_mp(Left->last_approx_mp, Left->last_approx_prec);
      change_size_point_mp(Left->last_approx_mp, Right->last_approx_size);
      Left->last_approx_mp->size = Right->last_approx_size;
    }

    // setup finalT & endPt
    if (Right->curr_prec < 64)
    { // setup using _d
      set_d(Left->finalT_d, (*coeff)[count]);
      count++;

      size = Right->pt_size;
      for (i = 0; i < size; i++)
      {
        set_d(&Left->endPt_d->coord[i], (*coeff)[count]);
        count++;
      }
    }
    else
    { // setup using _mp
      // setup real part
      mpf_set_str(Left->finalT_mp->r, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->finalT_mp->i, &(*epStr)[strLoc], base);
      strLoc += 1 + strlen(&(*epStr)[strLoc]);

      // setup endPt
      size = Right->pt_size;
      for (i = 0; i < size; i++)
      { // setup real part
        mpf_set_str(Left->endPt_mp->coord[i].r, &(*epStr)[strLoc], base);
        strLoc += 1 + strlen(&(*epStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->endPt_mp->coord[i].i, &(*epStr)[strLoc], base);
        strLoc += 1 + strlen(&(*epStr)[strLoc]);
      }
    }

    // copy last_approx
    if (Right->last_approx_prec < 64)
    { // setup using _d
      size = Right->last_approx_size;
      for (i = 0; i < size; i++)
      {
        set_d(&Left->last_approx_d->coord[i], (*coeff)[count]);
        count++;
      }
    }
    else
    { // setup using _mp
      size = Right->last_approx_size;
      for (i = 0; i < size; i++)
      { // setup real part
        mpf_set_str(Left->last_approx_mp->coord[i].r, &(*epStr)[strLoc], base);
        strLoc += 1 + strlen(&(*epStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->last_approx_mp->coord[i].i, &(*epStr)[strLoc], base);
        strLoc += 1 + strlen(&(*epStr)[strLoc]);
      }
    }

    // clear data
    if (freeData)
    {
      free(*coeff);
      free(*epStr);
    }
  }

  return;
}

void cp_witnessCodim_t_int(void *Out, void *In, comp_d **witComp, char **witStr, int **expW, int **witTypes, int **mult, endpoint_data_d_int **wit_d, endpoint_data_mp_int **wit_mp, endpoint_data_amp_int **wit_amp, int curr_prec, int freeData, int MPType, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int,In is _t                 *
*                otherwise - Out is _t, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10, sizeR = 0, sizeI = 0;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    witnessCodim_t *Right = (witnessCodim_t *)In;
    witnessCodim_t_int *Left = (witnessCodim_t_int *)Out;

    // setup curr_prec
    Left->curr_prec = curr_prec;

    // setup A_rows & A_cols
    Left->A_rows = Right->A_rows;
    Left->A_cols = Right->A_cols;

    // setup expW 
    *expW = (int *)bmalloc(Left->A_rows * Left->A_cols * sizeof(int));
    for (i = 0; i < Left->A_rows; i++)
      for (j = 0; j < Left->A_cols; j++)
        (*expW)[i * Left->A_cols + j] = Right->W[i][j];

    // setup witTypes & mult
    *witTypes = (int *)bmalloc(Right->num_set * sizeof(int));
    *mult = (int *)bmalloc(Right->num_set * sizeof(int));
    for (i = 0; i < Right->num_set; i++)
    {
      (*witTypes)[i] = Right->witnessPt_types[i];
      (*mult)[i] = Right->multiplicities[i];
    }

    // setup other vaules
    Left->codim = Right->codim;
    Left->num_set = Right->num_set;
    Left->num_nonsing = Right->num_nonsing;
    Left->num_sing = Right->num_sing;

    if (MPType == 0)
    { // setup _d
      int num_comp_codim_d = Right->A_rows * Right->A_cols + Right->H_d->size + 1 + Right->B_d->rows * Right->B_d->cols + Right->p_d->size;
      *witComp = (comp_d *)bmalloc(num_comp_codim_d * sizeof(comp_d));

      // setup the sizes
      Left->H_size = Right->H_d->size;
      Left->B_rows = Right->B_d->rows;  
      Left->B_cols = Right->B_d->cols;
      Left->p_size = Right->p_d->size;

      // copy A
      for (i = 0; i < Left->A_rows; i++)
        for (j = 0; j < Left->A_cols; j++)
        {
          set_d((*witComp)[count], &Right->A_d->entry[i][j]);
          count++;
        }
      // copy H
      for (i = 0; i < Left->H_size; i++)
      {
        set_d((*witComp)[count], &Right->H_d->coord[i]);
        count++;
      }
      // copy homVarConst
      set_d((*witComp)[count], Right->homVarConst_d);
      count++;
      // copy B
      for (i = 0; i < Left->B_rows; i++)
        for (j = 0; j < Left->B_cols; j++)
        {
          set_d((*witComp)[count], &Right->B_d->entry[i][j]);
          count++;
        }
      // copy p
      for (i = 0; i < Left->p_size; i++)
      {
        set_d((*witComp)[count], &Right->p_d->coord[i]);
        count++;
      }

      // setup the witness points
      comp_d *wit_comp_d = NULL;
      *wit_d = (endpoint_data_d_int *)bmalloc(Left->num_set * sizeof(endpoint_data_d_int));
      *wit_mp = NULL;
      *wit_amp = NULL;

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        cp_endpoint_data_d_int(&(*wit_d)[i], &Right->witnessPts_d[i], &wit_comp_d, freeData, 0, inType);

        *witComp = (comp_d *)brealloc(*witComp, (count + (*wit_d)[i].num_comp_d) * sizeof(comp_d));
        for (j = 0; j < (*wit_d)[i].num_comp_d; j++)
        {
          set_d((*witComp)[count], wit_comp_d[j]);
          count++;
        }
        free(wit_comp_d);
      }

      // store the sizes
      Left->num_comp_d = count;
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    { // setup _mp

      // setup the sizes
      Left->H_size = Right->H_mp->size;
      Left->B_rows = Right->B_mp->rows;
      Left->B_cols = Right->B_mp->cols;
      Left->p_size = Right->p_mp->size;

      // copy A
      for (i = 0; i < Left->A_rows; i++)
        for (j = 0; j < Left->A_cols; j++)
        { // real part
          strR = mpf_to_str(Right->A_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->A_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update witStr
          *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*witStr)[strLoc], strR);
          strcpy(&(*witStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy H
      for (i = 0; i < Left->H_size; i++)
      { // real part
        strR = mpf_to_str(Right->H_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->H_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // copy homVarConst
      // real part
      strR = mpf_to_str(Right->homVarConst_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->homVarConst_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update witStr
      *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*witStr)[strLoc], strR);
      strcpy(&(*witStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // copy B
      for (i = 0; i < Left->B_rows; i++)
        for (j = 0; j < Left->B_cols; j++)
        { // real part
          strR = mpf_to_str(Right->B_mp->entry[i][j].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->B_mp->entry[i][j].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update witStr
          *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*witStr)[strLoc], strR);
          strcpy(&(*witStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy p
      for (i = 0; i < Left->p_size; i++)
      { // real part
        strR = mpf_to_str(Right->p_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->p_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup the witness points
      char *str = NULL;
      *wit_mp = (endpoint_data_mp_int *)bmalloc(Left->num_set * sizeof(endpoint_data_mp_int));
      *wit_d = NULL;
      *wit_amp = NULL;

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        cp_endpoint_data_mp_int(&(*wit_mp)[i], &Right->witnessPts_mp[i], &str, freeData, 0, inType);
        *witStr = (char *)brealloc(*witStr, (strLoc + (*wit_mp)[i].totalLength) * sizeof(char));
        for (j = 0; j < (*wit_mp)[i].totalLength; j++)
        {
          (*witStr)[strLoc] = str[j];
          strLoc++;
        }
        free(str);
      }

      // store the sizes
      *witComp = NULL;
      Left->num_comp_d = 0;
      Left->totalLength = strLoc;
    }
    else
    { // setup _amp

      // setup the sizes
      Left->H_size = Right->H_d->size;
      Left->B_rows = Right->B_d->rows;
      Left->B_cols = Right->B_d->cols;
      Left->p_size = Right->p_d->size;

      // copy A
      for (i = 0; i < Left->A_rows; i++)
        for (j = 0; j < Left->A_cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->A_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->A_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update witStr
          *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*witStr)[strLoc], strR);
          strcpy(&(*witStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }
      // copy H
      for (i = 0; i < Left->H_size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->H_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->H_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }
      // copy homVarConst
      // real part
      strR = mpq_get_str(NULL, base, Right->homVarConst_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->homVarConst_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update witStr
      *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*witStr)[strLoc], strR);
      strcpy(&(*witStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // copy B
      for (i = 0; i < Left->B_rows; i++)
        for (j = 0; j < Left->B_cols; j++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->B_rat[i][j][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->B_rat[i][j][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update witStr
          *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*witStr)[strLoc], strR);
          strcpy(&(*witStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

      // copy p
      for (i = 0; i < Left->p_size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->p_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->p_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup the witness points
      comp_d *wit_comp_d = NULL;
      char *str = NULL;
      *wit_amp = (endpoint_data_amp_int *)bmalloc(Left->num_set * sizeof(endpoint_data_amp_int));
      *wit_d = NULL;
      *wit_mp = NULL;

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        cp_endpoint_data_amp_int(&(*wit_amp)[i], &Right->witnessPts_amp[i], &wit_comp_d, &str, freeData, 0, inType);

        *witComp = (comp_d *)brealloc(*witComp, (count + (*wit_amp)[i].num_comp_d) * sizeof(comp_d));
        for (j = 0; j < (*wit_amp)[i].num_comp_d; j++)
        {
          set_d((*witComp)[count], wit_comp_d[j]);
          count++;
        }
        free(wit_comp_d);

        *witStr = (char *)brealloc(*witStr, (strLoc + (*wit_amp)[i].totalLength) * sizeof(char));
        for (j = 0; j < (*wit_amp)[i].totalLength; j++)
        {
          (*witStr)[strLoc] = str[j];
          strLoc++;
        }
        free(str);
      }

      // store the sizes
      Left->num_comp_d = count;
      Left->totalLength = strLoc;
    }
  }
  else
  { // _t_int to _t so that it can be used
    witnessCodim_t_int *Right = (witnessCodim_t_int *)In;
    witnessCodim_t *Left = (witnessCodim_t *)Out;

    // setup A_rows & A_cols
    Left->A_rows = Right->A_rows;
    Left->A_cols = Right->A_cols;

    // setup W
    Left->W = (int **)bmalloc(Left->A_rows * sizeof(int *));
    for (i = 0; i < Left->A_rows; i++)
    {
      Left->W[i] = (int *)bmalloc(Left->A_cols * sizeof(int));
      for (j = 0; j < Left->A_cols; j++)
        Left->W[i][j] = (*expW)[i * Left->A_cols + j];
    }

    // setup witnessPt_types & multiplicities
    Left->witnessPt_types = (int *)bmalloc(Right->num_set * sizeof(int));
    Left->multiplicities = (int *)bmalloc(Right->num_set * sizeof(int));
    for (i = 0; i < Right->num_set; i++)
    {
      Left->witnessPt_types[i] = (*witTypes)[i];
      Left->multiplicities[i] = (*mult)[i];
    }

    // setup other vaules
    Left->codim = Right->codim;
    Left->num_set = Right->num_set;
    Left->num_nonsing = Right->num_nonsing;
    Left->num_sing = Right->num_sing;

    Left->num_components = 0;
    Left->component_nums = Left->deflations_needed = NULL;

    if (MPType == 0)
    { // setup A
      init_mat_d(Left->A_d, Right->A_rows, Right->A_cols);
      Left->A_d->rows = Right->A_rows;
      Left->A_d->cols = Right->A_cols;
      for (i = 0; i < Right->A_rows; i++)
        for (j = 0; j < Right->A_cols; j++)
        {
          set_d(&Left->A_d->entry[i][j], (*witComp)[count]);
          count++;
        }
      // setup H
      init_vec_d(Left->H_d, Right->H_size);
      Left->H_d->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      {
        set_d(&Left->H_d->coord[i], (*witComp)[count]);
        count++;
      }
      // setup homVarConst
      set_d(Left->homVarConst_d, (*witComp)[count]);
      count++;
      // setup B
      init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
      Left->B_d->rows = Right->B_rows;
      Left->B_d->cols = Right->B_cols;
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)
        {
          set_d(&Left->B_d->entry[i][j], (*witComp)[count]);
          count++;
        }
      // setup p
      init_vec_d(Left->p_d, Right->p_size);
      Left->p_d->size = Right->p_size;
      for (i = 0; i < Right->p_size; i++)
      {
        set_d(&Left->p_d->coord[i], (*witComp)[count]);
        count++;
      }

      // setup the witness points
      comp_d *wit_comp_d = NULL;
      Left->witnessPts_d = (endpoint_data_d *)bmalloc(Left->num_set * sizeof(endpoint_data_d));
      Left->witnessPts_mp = NULL;
      Left->witnessPts_amp = NULL;

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        wit_comp_d = &(*witComp)[count];
        cp_endpoint_data_d_int(&Left->witnessPts_d[i], &(*wit_d)[i], &wit_comp_d, 0, 1, inType);
        count += (*wit_d)[i].num_comp_d;
      }
      wit_comp_d = NULL;
    }
    else if (MPType == 1)
    { // setup A
      init_mat_mp(Left->A_mp, Right->A_rows, Right->A_cols);
      Left->A_mp->rows = Right->A_rows;
      Left->A_mp->cols = Right->A_cols;
      for (i = 0; i < Right->A_rows; i++)
        for (j = 0; j < Right->A_cols; j++)
        { // setup real part
          mpf_set_str(Left->A_mp->entry[i][j].r, &(*witStr)[strLoc], base);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->A_mp->entry[i][j].i, &(*witStr)[strLoc], base);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
        }
      // setup H
      init_vec_mp(Left->H_mp, Right->H_size);
      Left->H_mp->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      { // setup real part
        mpf_set_str(Left->H_mp->coord[i].r, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->H_mp->coord[i].i, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
      }
      // setup homVarConst
      init_mp(Left->homVarConst_mp);
      // setup real part
      mpf_set_str(Left->homVarConst_mp->r, &(*witStr)[strLoc], base);
      strLoc += 1 + strlen(&(*witStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->homVarConst_mp->i, &(*witStr)[strLoc], base);
      strLoc += 1 + strlen(&(*witStr)[strLoc]);

      // setup B
      init_mat_mp(Left->B_mp, Right->B_rows, Right->B_cols);
      Left->B_mp->rows = Right->B_rows;
      Left->B_mp->cols = Right->B_cols;
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)
        { // setup real part
          mpf_set_str(Left->B_mp->entry[i][j].r, &(*witStr)[strLoc], base);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->B_mp->entry[i][j].i, &(*witStr)[strLoc], base);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
        }
      // setup p
      init_vec_mp(Left->p_mp, Right->p_size);
      Left->p_mp->size = Right->p_size;
      for (i = 0; i < Right->p_size; i++)
      { // setup real part
        mpf_set_str(Left->p_mp->coord[i].r, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->p_mp->coord[i].i, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
      }

      // setup the witness points
      char *str = NULL;
      Left->witnessPts_d = NULL;
      Left->witnessPts_mp = (endpoint_data_mp *)bmalloc(Left->num_set * sizeof(endpoint_data_mp));
      Left->witnessPts_amp = NULL;

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        str = &(*witStr)[strLoc];
        cp_endpoint_data_mp_int(&Left->witnessPts_mp[i], &(*wit_mp)[i], &str, 0, 1, inType);
        strLoc += (*wit_mp)[i].totalLength;
      }
      str = NULL;
    }
    else
    { // setup A
      init_mat_d(Left->A_d, Right->A_rows, Right->A_cols);
      init_mat_mp2(Left->A_mp, Right->A_rows, Right->A_cols, Right->curr_prec);
      init_mat_rat(Left->A_rat, Right->A_rows, Right->A_cols);
      Left->A_d->rows = Left->A_mp->rows = Right->A_rows;
      Left->A_d->cols = Left->A_mp->cols = Right->A_cols;
      for (i = 0; i < Right->A_rows; i++)
        for (j = 0; j < Right->A_cols; j++)
        { // setup real part
          mpq_set_str(Left->A_rat[i][j][0], &(*witStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][0]);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->A_rat[i][j][1], &(*witStr)[strLoc], base);
          mpq_canonicalize(Left->A_rat[i][j][1]);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);

          // setup A_d & A_mp
          mpf_set_q(Left->A_mp->entry[i][j].r, Left->A_rat[i][j][0]);
          Left->A_d->entry[i][j].r = mpq_get_d(Left->A_rat[i][j][0]);

          mpf_set_q(Left->A_mp->entry[i][j].i, Left->A_rat[i][j][1]);
          Left->A_d->entry[i][j].i = mpq_get_d(Left->A_rat[i][j][1]);
        }
      // setup H
      init_vec_d(Left->H_d, Right->H_size);
      init_vec_mp2(Left->H_mp, Right->H_size, Right->curr_prec);
      init_vec_rat(Left->H_rat, Right->H_size);
      Left->H_d->size = Left->H_mp->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      { // setup real part
        mpq_set_str(Left->H_rat[i][0], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][0]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->H_rat[i][1], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][1]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);

        // setup H_d & H_mp
        mpf_set_q(Left->H_mp->coord[i].r, Left->H_rat[i][0]);
        Left->H_d->coord[i].r = mpq_get_d(Left->H_rat[i][0]);

        mpf_set_q(Left->H_mp->coord[i].i, Left->H_rat[i][1]);
        Left->H_d->coord[i].i = mpq_get_d(Left->H_rat[i][1]);
      }
      // setup homVarConst
      init_mp2(Left->homVarConst_mp, Right->curr_prec);
      init_rat(Left->homVarConst_rat);
      // setup real part
      mpq_set_str(Left->homVarConst_rat[0], &(*witStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[0]);
      strLoc += 1 + strlen(&(*witStr)[strLoc]);
      // setup imag part
      mpq_set_str(Left->homVarConst_rat[1], &(*witStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[1]);
      strLoc += 1 + strlen(&(*witStr)[strLoc]);
      // setup homVarConst_d & _mp
      mpf_set_q(Left->homVarConst_mp->r, Left->homVarConst_rat[0]);
      Left->homVarConst_d->r = mpq_get_d(Left->homVarConst_rat[0]);
      mpf_set_q(Left->homVarConst_mp->i, Left->homVarConst_rat[1]);
      Left->homVarConst_d->i = mpq_get_d(Left->homVarConst_rat[1]);
      // setup B
      init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
      init_mat_mp2(Left->B_mp, Right->B_rows, Right->B_cols, Right->curr_prec);
      init_mat_rat(Left->B_rat, Right->B_rows, Right->B_cols);
      Left->B_d->rows = Left->B_mp->rows = Right->B_rows;
      Left->B_d->cols = Left->B_mp->cols = Right->B_cols;
      for (i = 0; i < Right->B_rows; i++)
        for (j = 0; j < Right->B_cols; j++)
        { // setup real part
          mpq_set_str(Left->B_rat[i][j][0], &(*witStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][0]);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->B_rat[i][j][1], &(*witStr)[strLoc], base);
          mpq_canonicalize(Left->B_rat[i][j][1]);
          strLoc += 1 + strlen(&(*witStr)[strLoc]);

          // setup B_d & B_mp
          mpf_set_q(Left->B_mp->entry[i][j].r, Left->B_rat[i][j][0]);
          Left->B_d->entry[i][j].r = mpq_get_d(Left->B_rat[i][j][0]);

          mpf_set_q(Left->B_mp->entry[i][j].i, Left->B_rat[i][j][1]);
          Left->B_d->entry[i][j].i = mpq_get_d(Left->B_rat[i][j][1]);
        }
      // setup p
      init_vec_d(Left->p_d, Right->p_size);
      init_vec_mp2(Left->p_mp, Right->p_size, Right->curr_prec);
      init_vec_rat(Left->p_rat, Right->p_size);
      Left->p_d->size = Left->p_mp->size = Right->p_size;
      for (i = 0; i < Right->p_size; i++)
      { // setup real part
        mpq_set_str(Left->p_rat[i][0], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][0]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->p_rat[i][1], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->p_rat[i][1]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);

        // setup p_d & p_mp
        mpf_set_q(Left->p_mp->coord[i].r, Left->p_rat[i][0]);
        Left->p_d->coord[i].r = mpq_get_d(Left->p_rat[i][0]);

        mpf_set_q(Left->p_mp->coord[i].i, Left->p_rat[i][1]);
        Left->p_d->coord[i].i = mpq_get_d(Left->p_rat[i][1]);
      }

      // setup the witness points
      char *str = NULL;
      comp_d *wit_comp_d = NULL;
      Left->witnessPts_d = NULL;
      Left->witnessPts_mp = NULL;
      Left->witnessPts_amp = (endpoint_data_amp *)bmalloc(Left->num_set * sizeof(endpoint_data_amp));

      for (i = 0; i < Left->num_set; i++)
      { // setup the ith witness point
        wit_comp_d = &(*witComp)[count];
        str = &(*witStr)[strLoc];
        cp_endpoint_data_amp_int(&Left->witnessPts_amp[i], &(*wit_amp)[i], &wit_comp_d, &str, 0, 1, inType);
        count += (*wit_amp)[i].num_comp_d;
        strLoc += (*wit_amp)[i].totalLength;
      }
      str = NULL;
      wit_comp_d = NULL;
    }

    // clear data
    if (freeData)
    {
      free(*witComp);
      free(*witStr);
      free(*expW);
      free(*witTypes);
      free(*mult);
      if (MPType == 0)
        free(*wit_d);
      else if (MPType == 1)
        free(*wit_mp);
      else
        free(*wit_amp);
    }
  }

  return;
}

void cp_witness_t_int(void *Out, void *In, int **progInst, int **gpSizes, char **progStr, int **PPDtype, int **PPDsize, int **origDegs, int **newDegs, int **P, comp_d **witComp, char **witStr, int freeData, int MPType, int inType) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int,In is _t                 *
*                otherwise - Out is _t, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10, sizeR = 0, sizeI = 0;
  char *strR = NULL, *strI = NULL;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    witness_t *Right = (witness_t *)In;
    witness_t_int *Left = (witness_t_int *)Out;

    // copy prog
    cp_prog_t_int(&Left->Prog_int, Right->Prog, progInst, gpSizes, progStr, freeData, inType);
    // copy PPD
    cp_preproc_data_int(&Left->PPD_int, &Right->PPD, PPDtype, PPDsize, inType);

    // setup other data
    Left->num_funcs = Right->num_funcs;
    Left->num_codim = Right->num_codim;
    Left->curr_precision = Right->curr_precision;
    Left->system_rank = Right->system_rank;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup origDegs, newDegs & P
    *origDegs = (int *)bmalloc(Right->num_funcs * sizeof(int));
    *newDegs = (int *)bmalloc(Right->num_funcs * sizeof(int));
    *P = (int *)bmalloc(Right->num_funcs * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      (*origDegs)[i] = Right->orig_degrees[i];
      (*newDegs)[i] = Right->new_degrees[i];
      (*P)[i] = Right->P[i];
    }

    // setup the rest
    if (Left->new_variables != Left->orig_variables)
    { // need to setup C & gamma
      if (MPType == 0)
      {
        Left->num_comp_d = 1 + Left->new_variables * Left->orig_variables;
        *witComp = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));

        // copy C
        for (i = 0; i < Left->orig_variables; i++)
          for (j = 0; j < Left->new_variables; j++)
          {
            set_d((*witComp)[count], &Right->C_d->entry[i][j]);
            count++;
          }
        // copy gamma
        set_d((*witComp)[count], Right->gamma_d);
        count++;

        *witStr = NULL;
        Left->totalLength = 0;
      }
      else if (MPType == 1)
      { // copy C
        for (i = 0; i < Left->orig_variables; i++)
          for (j = 0; j < Left->new_variables; j++)
          { // real part
            strR = mpf_to_str(Right->C_mp->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->C_mp->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update witStr
            *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*witStr)[strLoc], strR);
            strcpy(&(*witStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;
 
            // free strR, strI
            free(strR);
            free(strI);
          }
        // copy gamma
        // real part
        strR = mpf_to_str(Right->gamma_mp->r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->gamma_mp->i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);

        *witComp = NULL; 
        Left->num_comp_d = 0;
        Left->totalLength = strLoc;
      }
      else
      { // copy C
        for (i = 0; i < Left->orig_variables; i++)
          for (j = 0; j < Left->new_variables; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->C_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->C_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update witStr
            *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*witStr)[strLoc], strR);
            strcpy(&(*witStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        // copy gamma
        // real part
        strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);

        *witComp = NULL;
        Left->num_comp_d = 0;
        Left->totalLength = strLoc;
      }
    }
    else
    { // need to gamma
      if (MPType == 0)
      {
        Left->num_comp_d = 1;
        *witComp = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));

        // copy gamma
        set_d((*witComp)[count], Right->gamma_d);
        count++;

        *witStr = NULL;
        Left->totalLength = 0;
      }
      else if (MPType == 1)
      { // copy gamma
        // real part
        strR = mpf_to_str(Right->gamma_mp->r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->gamma_mp->i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);

        *witComp = NULL;
        Left->num_comp_d = 0;
        Left->totalLength = strLoc;
      }
      else
      { // copy gamma
        // real part
        strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update witStr
        *witStr = (char *)brealloc(*witStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*witStr)[strLoc], strR);
        strcpy(&(*witStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);

        *witComp = NULL;
        Left->num_comp_d = 0;
        Left->totalLength = strLoc;
      }
    }
  }
  else
  { // _t_int to _t so that it can be used
    witness_t_int *Right = (witness_t_int *)In;
    witness_t *Left = (witness_t *)Out;

    // copy prog
    Left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    cp_prog_t_int(Left->Prog, &Right->Prog_int, progInst, gpSizes, progStr, freeData, inType);     

    // copy PPD
    cp_preproc_data_int(&Left->PPD, &Right->PPD_int, PPDtype, PPDsize, inType);

    // setup other data
    Left->num_funcs = Right->num_funcs;
    Left->num_codim = Right->num_codim;
    Left->curr_precision = Right->curr_precision;
    Left->system_rank = Right->system_rank;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;

    // setup origDegs, newDegs & P
    Left->orig_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->new_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->P = (int *)bmalloc(Right->num_funcs * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      Left->orig_degrees[i] = (*origDegs)[i];
      Left->new_degrees[i] = (*newDegs)[i];
      Left->P[i] = (*P)[i];
    }

    // setup the rest
    if (Left->new_variables != Left->orig_variables)
    { // need to setup C & gamma
      if (MPType == 0)
      { // setup C
        init_mat_d(Left->C_d, Right->orig_variables, Right->new_variables);
        Left->C_d->rows = Right->orig_variables;
        Left->C_d->cols = Right->new_variables;
        for (i = 0; i < Right->orig_variables; i++)
          for (j = 0; j < Right->new_variables; j++)
          {
            set_d(&Left->C_d->entry[i][j], (*witComp)[count]);
            count++;
          }
        // setup gamma
        set_d(Left->gamma_d, (*witComp)[count]);
        count++;
      }
      else if (MPType == 1)
      { // setup C
        init_mat_mp(Left->C_mp, Right->orig_variables, Right->new_variables);
        Left->C_mp->rows = Right->orig_variables;
        Left->C_mp->cols = Right->new_variables;
        for (i = 0; i < Right->orig_variables; i++)
          for (j = 0; j < Right->new_variables; j++)
          { // setup real part
            mpf_set_str(Left->C_mp->entry[i][j].r, &(*witStr)[strLoc], base);
            strLoc += 1 + strlen(&(*witStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->C_mp->entry[i][j].i, &(*witStr)[strLoc], base);
            strLoc += 1 + strlen(&(*witStr)[strLoc]);
          }
        // setup gamma
        init_mp(Left->gamma_mp);
        // setup real part
        mpf_set_str(Left->gamma_mp->r, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->gamma_mp->i, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
      }
      else
      { // setup C
        init_mat_d(Left->C_d, Right->orig_variables, Right->new_variables);
        init_mat_mp2(Left->C_mp, Right->orig_variables, Right->new_variables, Right->curr_precision);
        init_mat_rat(Left->C_rat, Right->orig_variables, Right->new_variables);
        Left->C_d->rows = Left->C_mp->rows = Right->orig_variables;
        Left->C_d->cols = Left->C_mp->cols = Right->new_variables;
        for (i = 0; i < Right->orig_variables; i++)
          for (j = 0; j < Right->new_variables; j++)
          { // setup real part
            mpq_set_str(Left->C_rat[i][j][0], &(*witStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][0]);
            strLoc += 1 + strlen(&(*witStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->C_rat[i][j][1], &(*witStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][1]);
            strLoc += 1 + strlen(&(*witStr)[strLoc]);

            // setup C_d & C_mp
            mpf_set_q(Left->C_mp->entry[i][j].r, Left->C_rat[i][j][0]);
            Left->C_d->entry[i][j].r = mpq_get_d(Left->C_rat[i][j][0]);

            mpf_set_q(Left->C_mp->entry[i][j].i, Left->C_rat[i][j][1]);
            Left->C_d->entry[i][j].i = mpq_get_d(Left->C_rat[i][j][1]);
          }
        // setup gamma
        init_mp2(Left->gamma_mp, Right->curr_precision);
        init_rat(Left->gamma_rat);
        // setup real part
        mpq_set_str(Left->gamma_rat[0], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->gamma_rat[0]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->gamma_rat[1], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->gamma_rat[1]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup gamma_d & _mp
        mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
        Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);
        mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
        Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);
      }
    }
    else
    { // need to gamma
      if (MPType == 0)
      { // setup gamma
        set_d(Left->gamma_d, (*witComp)[count]);
        count++;
      }
      else if (MPType == 1)
      { // setup gamma
        init_mp(Left->gamma_mp);
        // setup real part
        mpf_set_str(Left->gamma_mp->r, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->gamma_mp->i, &(*witStr)[strLoc], base);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
      }
      else
      { // setup gamma
        init_mp2(Left->gamma_mp, Right->curr_precision);
        init_rat(Left->gamma_rat);
        // setup real part
        mpq_set_str(Left->gamma_rat[0], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->gamma_rat[0]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->gamma_rat[1], &(*witStr)[strLoc], base);
        mpq_canonicalize(Left->gamma_rat[1]);
        strLoc += 1 + strlen(&(*witStr)[strLoc]);
        // setup gamma_d & _mp
        mpf_set_q(Left->gamma_mp->r, Left->gamma_rat[0]);
        Left->gamma_d->r = mpq_get_d(Left->gamma_rat[0]);
        mpf_set_q(Left->gamma_mp->i, Left->gamma_rat[1]);
        Left->gamma_d->i = mpq_get_d(Left->gamma_rat[1]);
      }
    }

    Left->targetSliceInit = Left->curr_codim_index = 0;

    // clear data
    if (freeData)
    {
      free(*origDegs);
      free(*newDegs);
      free(*P);
      free(*witStr);
    }
  }

  return;
}

void cp_trackBack_samples_t_int(void *TB_out, void *TB_in, char **egStr, comp_d **egCoeff, double **tbDouble, comp_d **tbCoeff, char **tbStr, int **tbInt, int freeStr, int freeCoeff, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - TB_out is _t_int, TB_in is _t          *
*                otherwise - TB_out is _t, TB_in is _t_int      *
* RETURN VALUES:                                                *
* NOTES: stores a copy of TB_in to TB_out                       *
\***************************************************************/
{
  int i, j, k, count = 0, strLoc = 0, base = 10;
  char *tempStr = NULL;  

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    trackBack_samples_t *TB_right = (trackBack_samples_t *)TB_in;
    trackBack_samples_t_int *TB_left = (trackBack_samples_t_int *)TB_out;

    // setup endgame data information
    cp_endgame_data_t_int(&TB_left->EG_int, &TB_right->endPt, egStr, egCoeff, 0, 0, inType);

    // setup the rest of the information
    TB_left->numSamples = TB_right->numSamples;
    TB_left->samplePts_prec = TB_right->samplePts_prec;
    TB_left->midPt_prec = TB_right->midPt_prec;

    // setup sample points
    if (TB_left->numSamples > 0)
    {
      if (TB_right->samplePts_prec < 64)
      { // NULL out tbStr
        *tbStr = NULL;
        TB_left->totalLength = 0;

        // setup tbDouble
        TB_left->num_double = TB_right->numSamples;
        *tbDouble = (double *)bmalloc(TB_right->numSamples * sizeof(double));
        for (i = 0; i < TB_right->numSamples; i++)
          (*tbDouble)[i] = TB_right->normSamples_d[i];

        // setup sample points
        TB_left->pointSize = TB_right->samplePts_d[0]->size;
       TB_left->num_comp_d = TB_left->pointSize * TB_right->numSamples;
         *tbCoeff = (comp_d *)bmalloc(TB_left->num_comp_d * sizeof(comp_d));
        for (i = 0; i < TB_right->numSamples; i++)
          for (j = 0; j < TB_left->pointSize; j++)
          {
            set_d((*tbCoeff)[count], &TB_right->samplePts_d[i]->coord[j]);
            count++;
          }
      }
      else
      { // NULL out tbDouble
        TB_left->num_double = 0;
        *tbDouble = NULL;

        // NULL out tbCoeff
        TB_left->num_comp_d = 0;
        *tbCoeff = NULL;

        // setup normSamples_mp
        *tbStr = NULL;
        for (i = 0; i < TB_right->numSamples; i++)
        {
          tempStr = mpf_to_str(TB_right->normSamples_mp[i], base);
          count = strLoc + strlen(tempStr) + 1;
          // reallocate tbStr
          *tbStr = (char *)brealloc(*tbStr, count * sizeof(char));
          // update tbStr
          for (j = strLoc; j < count; j++)
            (*tbStr)[j] = tempStr[j - strLoc];
          free(tempStr);
          // update strLoc
          strLoc = count;
        }

        // setup sample points
        TB_left->pointSize = TB_right->samplePts_mp[0]->size;
        for (i = 0; i < TB_right->numSamples; i++)
          for (j = 0; j < TB_left->pointSize; j++)
          { // copy real part
            tempStr = mpf_to_str(TB_right->samplePts_mp[i]->coord[j].r, base);
            count = strLoc + strlen(tempStr) + 1;
            // reallocate tbStr
            *tbStr = (char *)brealloc(*tbStr, count * sizeof(char));
            // update tbStr
            for (k = strLoc; k < count; k++)
              (*tbStr)[k] = tempStr[k - strLoc];
            free(tempStr);
            // update strLoc
            strLoc = count;

            // copy imag part
            tempStr = mpf_to_str(TB_right->samplePts_mp[i]->coord[j].i, base);
            count = strLoc + strlen(tempStr) + 1;
            // reallocate tbStr
            *tbStr = (char *)brealloc(*tbStr, count * sizeof(char));
            // update tbStr
            for (k = strLoc; k < count; k++)
              (*tbStr)[k] = tempStr[k - strLoc];
            free(tempStr);
            // update strLoc
            strLoc = count;
          }

        // store the size
        TB_left->totalLength = strLoc;
      }

      // setup mid points
      if (TB_right->midPt_prec < 64)
      { // add on to tbCoeff
        *tbCoeff = (comp_d *)brealloc(*tbCoeff, (TB_left->num_comp_d + TB_left->pointSize * TB_right->numSamples) * sizeof(comp_d));
        count = TB_left->num_comp_d;
        TB_left->num_comp_d += TB_left->pointSize * TB_right->numSamples;
        for (i = 0; i < TB_right->numSamples; i++)
          for (j = 0; j < TB_left->pointSize; j++)
          {
            set_d((*tbCoeff)[count], &TB_right->midPt_d[i]->coord[j]);
            count++;
          }
      }
      else
      { // add on to tbStr
        for (i = 0; i < TB_right->numSamples; i++)
          for (j = 0; j < TB_left->pointSize; j++)
          { // copy real part
            tempStr = mpf_to_str(TB_right->midPt_mp[i]->coord[j].r, base);
            count = strLoc + strlen(tempStr) + 1;
            // reallocate tbStr
            *tbStr = (char *)brealloc(*tbStr, count * sizeof(char));
            // update tbStr
            for (k = strLoc; k < count; k++)
              (*tbStr)[k] = tempStr[k - strLoc];
            free(tempStr);
            // update strLoc
            strLoc = count;

            // copy imag part
            tempStr = mpf_to_str(TB_right->midPt_mp[i]->coord[j].i, base);
            count = strLoc + strlen(tempStr) + 1;
            // reallocate tbStr
            *tbStr = (char *)brealloc(*tbStr, count * sizeof(char));
            // update tbStr
            for (k = strLoc; k < count; k++)
              (*tbStr)[k] = tempStr[k - strLoc];
            free(tempStr);
            // update strLoc
            strLoc = count;
          }

        // store the size
        TB_left->totalLength = strLoc;
      }

      // setup tbInt
      *tbInt = (int *)bmalloc(TB_right->numSamples * sizeof(int));
      for (i = 0; i < TB_right->numSamples; i++)
        (*tbInt)[i] = TB_right->samplePts_notUsed[i];
    }
    else
    { // NULL out everything and set to 0
      *tbDouble = NULL;
      *tbCoeff = NULL;
      *tbStr = NULL;
      *tbInt = NULL;

      TB_left->pointSize = TB_left->num_double = TB_left->num_comp_d = TB_left->totalLength = 0;
    }
  }
  else
  { // _t_int to _t so that it can be used by the workers
    trackBack_samples_t *TB_left = (trackBack_samples_t *)TB_out;
    trackBack_samples_t_int *TB_right = (trackBack_samples_t_int *)TB_in;

    // setup endgame data information
    cp_endgame_data_t_int(&TB_left->endPt, &TB_right->EG_int, egStr, egCoeff, freeStr, freeCoeff, inType);

    // setup the rest of the information
    TB_left->numSamples = TB_right->numSamples;
    TB_left->samplePts_prec = TB_right->samplePts_prec;
    TB_left->midPt_prec = TB_right->midPt_prec;

    // setup samplePts_notUsed
    TB_left->samplePts_notUsed = (int *)bmalloc(TB_right->numSamples * sizeof(int));
    for (i = 0; i < TB_right->numSamples; i++)
      TB_left->samplePts_notUsed[i] = (*tbInt)[i];

    // setup sample points
    if (TB_right->samplePts_prec < 64)
    { // setup normSamples_d
      TB_left->normSamples_d = (double *)bmalloc(TB_right->numSamples * sizeof(double));
      for (i = 0; i < TB_right->numSamples; i++)
        TB_left->normSamples_d[i] = (*tbDouble)[i];

      // setup samplePts_d
      TB_left->samplePts_d = (point_d *)bmalloc(TB_right->numSamples * sizeof(point_d));
      for (i = 0; i < TB_right->numSamples; i++)
      { // initialize
        init_point_d(TB_left->samplePts_d[i], TB_right->pointSize);
        TB_left->samplePts_d[i]->size = TB_right->pointSize;
        // setup
        for (j = 0; j < TB_right->pointSize; j++)
        {
          set_d(&TB_left->samplePts_d[i]->coord[j], (*tbCoeff)[count]);
          count++;
        }
      }
    }
    else
    { // setup normSamples_mp
      TB_left->normSamples_mp = (mpf_t *)bmalloc(TB_right->numSamples * sizeof(mpf_t));
      for (i = 0; i < TB_right->numSamples; i++)
      { // initialize
        mpf_init2(TB_left->normSamples_mp[i], TB_right->samplePts_prec);
        mpf_set_str(TB_left->normSamples_mp[i], &(*tbStr)[strLoc], base);
        strLoc += strlen(&(*tbStr)[strLoc]) + 1;
      }

      // setup samplePts_mp
      TB_left->samplePts_mp = (point_mp *)bmalloc(TB_right->numSamples * sizeof(point_mp));
      for (i = 0; i < TB_right->numSamples; i++)
      { // initialize
        init_point_mp2(TB_left->samplePts_mp[i], TB_right->pointSize, TB_right->samplePts_prec);
        TB_left->samplePts_mp[i]->size = TB_right->pointSize;
        // setup
        for (j = 0; j < TB_right->pointSize; j++)
        { // setup real
          mpf_set_str(TB_left->samplePts_mp[i]->coord[j].r, &(*tbStr)[strLoc], base);
          strLoc += strlen(&(*tbStr)[strLoc]) + 1;
          // setup imag
          mpf_set_str(TB_left->samplePts_mp[i]->coord[j].i, &(*tbStr)[strLoc], base);
          strLoc += strlen(&(*tbStr)[strLoc]) + 1;
        }
      }
    }

    // setup mid points
    if (TB_right->midPt_prec < 64)
    { // setup midPt_d
      TB_left->midPt_d = (point_d *)bmalloc(TB_right->numSamples * sizeof(point_d));
      for (i = 0; i < TB_right->numSamples; i++)
      { // initialize
        init_point_d(TB_left->midPt_d[i], TB_right->pointSize);
        TB_left->midPt_d[i]->size = TB_right->pointSize;
        // setup
        for (j = 0; j < TB_right->pointSize; j++)
        {
          set_d(&TB_left->midPt_d[i]->coord[j], (*tbCoeff)[count]);
          count++;
        }
      }
    }
    else
    { // setup midPt_mp
      TB_left->midPt_mp = (point_mp *)bmalloc(TB_right->numSamples * sizeof(point_mp));
      for (i = 0; i < TB_right->numSamples; i++)
      { // initialize
        init_point_mp2(TB_left->midPt_mp[i], TB_right->pointSize, TB_right->midPt_prec);
        TB_left->midPt_mp[i]->size = TB_right->pointSize;
        // setup
        for (j = 0; j < TB_right->pointSize; j++)
        { // setup real
          mpf_set_str(TB_left->midPt_mp[i]->coord[j].r, &(*tbStr)[strLoc], base);
          strLoc += strlen(&(*tbStr)[strLoc]) + 1;
          // setup imag
          mpf_set_str(TB_left->midPt_mp[i]->coord[j].i, &(*tbStr)[strLoc], base);
          strLoc += strlen(&(*tbStr)[strLoc]) + 1;
        }
      }
    }

    // clear
    if (freeStr)
      free(*tbStr);
    if (freeCoeff)
    {
      free(*tbDouble);
      free(*tbCoeff);
      free(*tbInt);
    }
  }

  return;
}

void cp_endgame_data_t(endgame_data_t *EG_out, endgame_data_t *EG_in, int initEG_out)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  if (initEG_out)
  { // initialize to the correct precision
    init_endgame_data(EG_out, EG_in->prec);
  }
  else
  { // setup to the correct precision
    clear_endgame_data(EG_out);
    init_endgame_data(EG_out, EG_in->prec);
  }

  // copy values
  EG_out->prec = EG_in->prec;
  EG_out->retVal = EG_in->retVal;
  EG_out->pathNum = EG_in->pathNum;
  EG_out->codim = EG_in->codim;
  EG_out->first_increase = EG_in->first_increase;
  EG_out->condition_number = EG_in->condition_number;
  if (EG_out->prec < 64)
  { // copy _d
    point_data_cp_d(&EG_out->PD_d, &EG_in->PD_d);
    EG_out->function_residual_d = EG_in->function_residual_d;
    EG_out->latest_newton_residual_d = EG_in->latest_newton_residual_d;
    EG_out->t_val_at_latest_sample_point_d = EG_in->t_val_at_latest_sample_point_d;
    EG_out->error_at_latest_sample_point_d = EG_in->error_at_latest_sample_point_d;
  }
  else
  { // copy _mp
    point_data_cp_mp(&EG_out->PD_mp, &EG_in->PD_mp);
    mpf_set(EG_out->function_residual_mp, EG_in->function_residual_mp);
    mpf_set(EG_out->latest_newton_residual_mp, EG_in->latest_newton_residual_mp);
    mpf_set(EG_out->t_val_at_latest_sample_point_mp, EG_in->t_val_at_latest_sample_point_mp);
    mpf_set(EG_out->error_at_latest_sample_point_mp, EG_in->error_at_latest_sample_point_mp);
  }

  EG_out->last_approx_prec = EG_in->last_approx_prec;
  if (EG_out->last_approx_prec < 64)
  {
    point_cp_d(EG_out->last_approx_d, EG_in->last_approx_d);
  }
  else
  {
    setprec_point_mp(EG_out->last_approx_mp, EG_out->last_approx_prec);
    point_cp_mp(EG_out->last_approx_mp, EG_in->last_approx_mp);
  }

  return;
}

void cp_trackBack_samples_t(trackBack_samples_t *TB_out, trackBack_samples_t *TB_in, int initTB_out)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  if (initTB_out)
  { // initialize to the correct precision
    init_trackBack_sample(TB_out, TB_in->endPt.prec);
  }
  else
  { // setup to the correct precision
    clear_trackBack_sample(TB_out);
    init_trackBack_sample(TB_out, TB_in->endPt.prec);
  }

  // copy endPt
  cp_endgame_data_t(&TB_out->endPt, &TB_in->endPt, 0);

  // copy other data
  TB_out->numSamples = TB_in->numSamples;
  TB_out->samplePts_prec = TB_in->samplePts_prec;
  TB_out->midPt_prec = TB_in->midPt_prec;

  // setup samplePts
  if (TB_out->samplePts_prec < 64)
  { // setup _d
    TB_out->normSamples_d = (double *)bmalloc(TB_out->numSamples * sizeof(double));
    TB_out->samplePts_d = (point_d *)bmalloc(TB_out->numSamples * sizeof(point_d)); 
    TB_out->samplePts_notUsed = (int *)bmalloc(TB_out->numSamples * sizeof(int));
    for (i = 0; i < TB_out->numSamples; i++)
    { // copy values
      TB_out->samplePts_notUsed[i] = TB_in->samplePts_notUsed[i];
      TB_out->normSamples_d[i] = TB_in->normSamples_d[i];
      init_point_d(TB_out->samplePts_d[i], 0);
      point_cp_d(TB_out->samplePts_d[i], TB_in->samplePts_d[i]);
    }
  }
  else
  { // setup _mp
    TB_out->normSamples_mp = (mpf_t *)bmalloc(TB_out->numSamples * sizeof(mpf_t));
    TB_out->samplePts_mp = (point_mp *)bmalloc(TB_out->numSamples * sizeof(point_mp));
    TB_out->samplePts_notUsed = (int *)bmalloc(TB_out->numSamples * sizeof(int));
    for (i = 0; i < TB_out->numSamples; i++)
    { // copy values
      TB_out->samplePts_notUsed[i] = TB_in->samplePts_notUsed[i];
      mpf_init2(TB_out->normSamples_mp[i], TB_out->samplePts_prec);
      init_point_mp2(TB_out->samplePts_mp[i], 0, TB_out->samplePts_prec);
      mpf_set(TB_out->normSamples_mp[i], TB_in->normSamples_mp[i]);
      point_cp_mp(TB_out->samplePts_mp[i], TB_in->samplePts_mp[i]);
    }
  }

  // setup midPt
  if (TB_out->midPt_prec < 64)
  { // setup _d
    TB_out->midPt_d = (point_d *)bmalloc(TB_out->numSamples * sizeof(point_d));
    for (i = 0; i < TB_out->numSamples; i++)
    { // copy values
      init_point_d(TB_out->midPt_d[i], 0);
      point_cp_d(TB_out->midPt_d[i], TB_in->midPt_d[i]);
    }
  }
  else
  { // setup _mp
    TB_out->midPt_mp = (point_mp *)bmalloc(TB_out->numSamples * sizeof(point_mp));
    for (i = 0; i < TB_out->numSamples; i++)
    { // copy values
      init_point_mp2(TB_out->midPt_mp[i], 0, TB_out->midPt_prec);
      point_cp_mp(TB_out->midPt_mp[i], TB_in->midPt_mp[i]);
    }
  }

  return;
}

void cp_regen_pos_dim_int(void *Out, void *In, int MPType, char **rpdStr, char **progStr, int freeStr, comp_d **coeff_d, int **degrees, int **ppd_type, int **ppd_size, int **prog_inst, int **prog_gp_sizes, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, k, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    regen_pos_dim_t *Right = (regen_pos_dim_t *)In;
    regen_pos_dim_t_int *Left = (regen_pos_dim_t_int *)Out;

    // setup Prog
    cp_prog_t_int(&Left->Prog_int, Right->Prog, prog_inst, prog_gp_sizes, progStr, freeStr, inType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD_int, &Right->PPD, ppd_type, ppd_size, inType);

    // setup other values
    Left->system_rank = Right->system_rank;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;
    Left->num_funcs = Right->num_funcs;
    Left->num_codim = Right->num_codim;
    Left->curr_precision = Right->curr_precision;
    Left->sameA = Right->sameA;

    // count the number of integers being sent - orig_degrees, new_degrees, P, W
    Left->num_int = 3 * Right->num_funcs;
    for (i = 0; i < Right->num_codim; i++)
      for (j = 0; j <= i; j++)
        Left->num_int += Right->num_funcs - j - 1;

    // setup degrees
    count = 0;
    *degrees = (int *)bmalloc(Left->num_int * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      (*degrees)[count] = Right->orig_degrees[i];
      count++;
    }
    for (i = 0; i < Right->num_funcs; i++)
    {
      (*degrees)[count] = Right->new_degrees[i];
      count++;
    }
    for (i = 0; i < Right->num_funcs; i++)
    {
      (*degrees)[count] = Right->P[i];
      count++;
    }
    for (i = 0; i < Right->num_codim; i++)
      for (j = 0; j <= i; j++)
        for (k = 0; k < Right->num_funcs - j - 1; k++)
        {
          (*degrees)[count] = Right->W[i][j][k];
          count++;
        }

    // setup coeff_d & rpdStr for C, H, homVarConst, gamma, coeff, patchCoeff & A
    count = strLoc = 0;
    *coeff_d = NULL;
    *rpdStr = NULL;
    if (MPType == 0)
    { // setup coeff_d
      if (Right->orig_variables != Right->new_variables)
      { // need to send C_d
        Left->C_rows = Right->C_d->rows;
        Left->C_cols = Right->C_d->cols;
        Left->num_comp_d = Left->C_rows * Left->C_cols;

        // setup coeff_d with C
        *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
        for (i = 0; i < Left->C_rows; i++)
          for (j = 0; j < Left->C_cols; j++)
          {
            set_d((*coeff_d)[count], &Right->C_d->entry[i][j]);
            count++;
          }
      }
      else
      { // set to 0
        Left->num_comp_d = Left->C_rows = Left->C_cols = 0;
      }

      // setup H
      Left->num_comp_d += Left->H_size = Right->H_d->size;
      *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
      for (i = 0; i < Left->H_size; i++)
      {
        set_d((*coeff_d)[count], &Right->H_d->coord[i]);
        count++;
      }

      // setup homVarConst
      Left->num_comp_d += 1;
      *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
      set_d((*coeff_d)[count], Right->homVarConst_d);
      count++;

      // setup gamma
      Left->num_comp_d += 1;
      *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
      set_d((*coeff_d)[count], Right->gamma_d);
      count++;

      // setup patchCoeff
      Left->num_comp_d += Left->patchCoeff_size = Right->patchCoeff_d->size;
      *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
      for (i = 0; i < Left->patchCoeff_size; i++)
      {
        set_d((*coeff_d)[count], &Right->patchCoeff_d->coord[i]);
        count++;
      }

      // setup A
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        Left->num_comp_d += (i + 1) * Right->num_funcs;
        *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          {
            set_d((*coeff_d)[count], &Right->A_d[i]->entry[j][k]);
            count++;
          }
      }

      // setup coeff
      for (i = 0; i < Right->num_codim; i++)
      {
        Left->num_comp_d += Right->new_degrees[i] * Right->new_variables;
        *coeff_d = (comp_d *)brealloc(*coeff_d, Left->num_comp_d * sizeof(comp_d));
        for (j = 0; j < Right->new_degrees[i]; j++)
          for (k = 0; k < Right->new_variables; k++)
          {
            set_d((*coeff_d)[count], Right->coeff_d[i][j][k]);
            count++;
          }
      }

      // set totalLength = 0
      Left->totalLength = 0;
    }
    else if (MPType == 1)
    { // setup rpdStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      // setup C, if needed
      if (Right->orig_variables != Right->new_variables)
      { // need to send C_mp
        Left->C_rows = Right->C_mp->rows;
        Left->C_cols = Right->C_mp->cols;

        for (i = 0; i < Left->C_rows; i++)
          for (j = 0; j < Left->C_cols; j++)
          { // real part
            strR = mpf_to_str(Right->C_mp->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->C_mp->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }
      else
      { // set to 0
        Left->C_rows = Left->C_cols = 0;
      }

      // setup H
      Left->H_size = Right->H_mp->size;
      for (i = 0; i < Left->H_size; i++)
      { // real part
        strR = mpf_to_str(Right->H_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->H_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update rpdStr
        *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*rpdStr)[strLoc], strR);
        strcpy(&(*rpdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup homVarConst
      // real part
      strR = mpf_to_str(Right->homVarConst_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->homVarConst_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update rpdStr
      *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*rpdStr)[strLoc], strR);
      strcpy(&(*rpdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // setup gamma
      // real part
      strR = mpf_to_str(Right->gamma_mp->r, base);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpf_to_str(Right->gamma_mp->i, base);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update rpdStr
      *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*rpdStr)[strLoc], strR);
      strcpy(&(*rpdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // setup patchCoeff
      Left->patchCoeff_size = Right->patchCoeff_mp->size;
      for (i = 0; i < Left->patchCoeff_size; i++)
      { // real part
        strR = mpf_to_str(Right->patchCoeff_mp->coord[i].r, base);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpf_to_str(Right->patchCoeff_mp->coord[i].i, base);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update rpdStr
        *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*rpdStr)[strLoc], strR);
        strcpy(&(*rpdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup A
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          { // real part
            strR = mpf_to_str(Right->A_mp[i]->entry[j][k].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->A_mp[i]->entry[j][k].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }

      // setup coeff
      for (i = 0; i < Right->num_codim; i++)
        for (j = 0; j < Right->new_degrees[i]; j++)
          for (k = 0; k < Right->new_variables; k++)
          { // real part
            strR = mpf_to_str(Right->coeff_mp[i][j][k]->r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->coeff_mp[i][j][k]->i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }

      // store the size
      Left->totalLength = strLoc;

      // set num_comp_d = 0
      Left->num_comp_d = 0;
    } 
    else // MPType == 2
    { // setup rpdStr
      int sizeR, sizeI;
      char *strR = NULL, *strI = NULL;

      // setup C, if needed
      if (Right->orig_variables != Right->new_variables)
      { // need to send C_mp
        Left->C_rows = Right->C_mp->rows;
        Left->C_cols = Right->C_mp->cols;

        for (i = 0; i < Left->C_rows; i++)
          for (j = 0; j < Left->C_cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->C_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->C_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }
      else
      { // set to 0
        Left->C_rows = Left->C_cols = 0;
      }

      // setup H
      Left->H_size = Right->H_mp->size;
      for (i = 0; i < Left->H_size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->H_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->H_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update rpdStr
        *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*rpdStr)[strLoc], strR);
        strcpy(&(*rpdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup homVarConst
      // real part
      strR = mpq_get_str(NULL, base, Right->homVarConst_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->homVarConst_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update rpdStr
      *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*rpdStr)[strLoc], strR);
      strcpy(&(*rpdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // setup gamma
      // real part
      strR = mpq_get_str(NULL, base, Right->gamma_rat[0]);
      sizeR = strlen(strR) + 1; // +1 for '\0'
      // imag part
      strI = mpq_get_str(NULL, base, Right->gamma_rat[1]);
      sizeI = strlen(strI) + 1; // +1 for '\0'

      // update rpdStr
      *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
      strcpy(&(*rpdStr)[strLoc], strR);
      strcpy(&(*rpdStr)[strLoc + sizeR], strI);
      // update strLoc
      strLoc += sizeR + sizeI;

      // free strR, strI
      free(strR);
      free(strI);

      // setup patchCoeff
      Left->patchCoeff_size = Right->patchCoeff_mp->size;
      for (i = 0; i < Left->patchCoeff_size; i++)
      { // real part
        strR = mpq_get_str(NULL, base, Right->patchCoeff_rat[i][0]);
        sizeR = strlen(strR) + 1; // +1 for '\0'
        // imag part
        strI = mpq_get_str(NULL, base, Right->patchCoeff_rat[i][1]);
        sizeI = strlen(strI) + 1; // +1 for '\0'

        // update rpdStr
        *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
        strcpy(&(*rpdStr)[strLoc], strR);
        strcpy(&(*rpdStr)[strLoc + sizeR], strI);
        // update strLoc
        strLoc += sizeR + sizeI;

        // free strR, strI
        free(strR);
        free(strI);
      }

      // setup A
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->A_rat[i][j][k][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->A_rat[i][j][k][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
      }

      // setup coeff
      for (i = 0; i < Right->num_codim; i++)
        for (j = 0; j < Right->new_degrees[i]; j++)
          for (k = 0; k < Right->new_variables; k++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->coeff_rat[i][j][k][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->coeff_rat[i][j][k][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }

      // store the size
      Left->totalLength = strLoc;

      // set num_comp_d = 0
      Left->num_comp_d = 0;
    }
  }
  else
  { // _t_int to _t so that it can be used
    regen_pos_dim_t_int *Right = (regen_pos_dim_t_int *)In;
    regen_pos_dim_t *Left = (regen_pos_dim_t *)Out;

    // setup Prog
    Left->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
    cp_prog_t_int(Left->Prog, &Right->Prog_int, prog_inst, prog_gp_sizes, progStr, freeStr, inType);
    // initialize evalProg
    initEvalProg(MPType);

    // setup PPD
    cp_preproc_data_int(&Left->PPD, &Right->PPD_int, ppd_type, ppd_size, inType);

    // setup other values
    Left->system_rank = Right->system_rank;
    Left->orig_variables = Right->orig_variables;
    Left->new_variables = Right->new_variables;
    Left->num_funcs = Right->num_funcs;
    Left->num_codim = Right->num_codim;
    Left->curr_precision = Right->curr_precision;
    Left->sameA = Right->sameA;

    // setup orig_degrees, new_degrees & P
    Left->orig_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->new_degrees = (int *)bmalloc(Right->num_funcs * sizeof(int));
    Left->P = (int *)bmalloc(Right->num_funcs * sizeof(int));
    for (i = 0; i < Right->num_funcs; i++)
    {
      Left->orig_degrees[i] = (*degrees)[i];
      Left->new_degrees[i] = (*degrees)[Right->num_funcs + i];
      Left->P[i] = (*degrees)[2 * Right->num_funcs + i];
    }

    // setup W 
    count = 3 * Right->num_funcs;
    Left->W = (int ***)bmalloc(Right->num_codim * sizeof(int **));
    for (i = 0; i < Right->num_codim; i++)
    {
      Left->W[i] = (int **)bmalloc((i + 1) * sizeof(int *));
      for (j = 0; j <= i; j++)
      {
        Left->W[i][j] = (int *)bmalloc((Right->num_funcs - j - 1) * sizeof(int));
        for (k = 0; k < Right->num_funcs - j - 1; k++)
        {
          Left->W[i][j][k] = (*degrees)[count];
          count++;
        }
      }
    }

    // clear degrees
    free(*degrees);

    // setup the rest
    count = strLoc = 0;
    if (MPType == 0)
    { // setup _d
      if (Right->orig_variables != Right->new_variables)
      { // setup C_d
        init_mat_d(Left->C_d, Right->C_rows, Right->C_cols);
        Left->C_d->rows = Right->C_rows;
        Left->C_d->cols = Right->C_cols;

        for (i = 0; i < Right->C_rows; i++)
          for (j = 0; j < Right->C_cols; j++)
          {
            set_d(&Left->C_d->entry[i][j], (*coeff_d)[count]);
            count++;
          }
      }

      // setup H
      init_vec_d(Left->H_d, Right->H_size);
      Left->H_d->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      {
        set_d(&Left->H_d->coord[i], (*coeff_d)[count]);
        count++;
      }

      // setup homVarConst
      set_d(Left->homVarConst_d, (*coeff_d)[count]);
      count++;

      // setup gamma
      set_d(Left->gamma_d, (*coeff_d)[count]);
      count++;

      // setup patchCoeff
      init_vec_d(Left->patchCoeff_d, Right->patchCoeff_size);
      Left->patchCoeff_d->size = Right->patchCoeff_size;
      for (i = 0; i < Right->patchCoeff_size; i++)
      {
        set_d(&Left->patchCoeff_d->coord[i], (*coeff_d)[count]);
        count++;
      }

      // setup A
      Left->A_d = (mat_d *)bmalloc(Right->num_codim * sizeof(mat_d));
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        init_mat_d(Left->A_d[i], i + 1, Right->num_funcs);
        Left->A_d[i]->rows = i + 1;
        Left->A_d[i]->cols = Right->num_funcs;
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          {
            set_d(&Left->A_d[i]->entry[j][k], (*coeff_d)[count]);
            count++;
          }
      }

      // setup coeff_d
      Left->coeff_d = (comp_d ***)bmalloc(Right->num_codim * sizeof(comp_d **));
      for (i = 0; i < Right->num_codim; i++)
      {
        Left->coeff_d[i] = (comp_d **)bmalloc(Left->new_degrees[i] * sizeof(comp_d *));
        for (j = 0; j < Left->new_degrees[i]; j++)
        {
          Left->coeff_d[i][j] = (comp_d *)bmalloc(Left->new_variables * sizeof(comp_d));
          for (k = 0; k < Left->new_variables; k++)
          {
            set_d(Left->coeff_d[i][j][k], (*coeff_d)[count]);
            count++;
          }
        }
      }

      // clear coeff_d
      free(*coeff_d);
    }
    else if (MPType == 1)
    { // setup _mp
      if (Right->orig_variables != Right->new_variables)
      { // setup C_mp
        init_mat_mp(Left->C_mp, Right->C_rows, Right->C_cols);
        Left->C_mp->rows = Right->C_rows;
        Left->C_mp->cols = Right->C_cols;

        for (i = 0; i < Right->C_rows; i++)
          for (j = 0; j < Right->C_cols; j++)
          { // setup real part
            mpf_set_str(Left->C_mp->entry[i][j].r, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->C_mp->entry[i][j].i, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          }
      }

      // setup H
      init_vec_mp(Left->H_mp, Right->H_size);
      Left->H_mp->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      { // setup real part
        mpf_set_str(Left->H_mp->coord[i].r, &(*rpdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->H_mp->coord[i].i, &(*rpdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      }

      // setup homVarConst
      init_mp(Left->homVarConst_mp);
      // setup real part
      mpf_set_str(Left->homVarConst_mp->r, &(*rpdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->homVarConst_mp->i, &(*rpdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

      // setup gamma
      init_mp(Left->gamma_mp);
      // setup real part
      mpf_set_str(Left->gamma_mp->r, &(*rpdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      // setup imag part
      mpf_set_str(Left->gamma_mp->i, &(*rpdStr)[strLoc], base);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

      // setup patchCoeff
      init_vec_mp(Left->patchCoeff_mp, Right->patchCoeff_size);
      Left->patchCoeff_mp->size = Right->patchCoeff_size;
      for (i = 0; i < Right->patchCoeff_size; i++)
      { // setup real part
        mpf_set_str(Left->patchCoeff_mp->coord[i].r, &(*rpdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
        // setup imag part
        mpf_set_str(Left->patchCoeff_mp->coord[i].i, &(*rpdStr)[strLoc], base);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      }

      // setup A
      Left->A_mp = (mat_mp *)bmalloc(Right->num_codim * sizeof(mat_mp));
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        init_mat_mp(Left->A_mp[i], i + 1, Right->num_funcs);
        Left->A_mp[i]->rows = i + 1;
        Left->A_mp[i]->cols = Right->num_funcs;
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          { // setup real part
            mpf_set_str(Left->A_mp[i]->entry[j][k].r, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->A_mp[i]->entry[j][k].i, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          }
      }

      // setup coeff_mp
      Left->coeff_mp = (comp_mp ***)bmalloc(Right->num_codim * sizeof(comp_mp **));
      for (i = 0; i < Right->num_codim; i++)
      {
        Left->coeff_mp[i] = (comp_mp **)bmalloc(Left->new_degrees[i] * sizeof(comp_mp *));
        for (j = 0; j < Left->new_degrees[i]; j++)
        {
          Left->coeff_mp[i][j] = (comp_mp *)bmalloc(Left->new_variables * sizeof(comp_mp));
          for (k = 0; k < Left->new_variables; k++)
          {
            init_mp(Left->coeff_mp[i][j][k]);
            // setup real part
            mpf_set_str(Left->coeff_mp[i][j][k]->r, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->coeff_mp[i][j][k]->i, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          }
        }
      }
    }
    else // MPType == 2  
    { // setup _d, _mp, _rat
      if (Right->orig_variables != Right->new_variables)
      { // setup C_mp
        init_mat_d(Left->C_d, Right->C_rows, Right->C_cols);
        init_mat_mp2(Left->C_mp, Right->C_rows, Right->C_cols, Left->curr_precision);
        init_mat_rat(Left->C_rat, Right->C_rows, Right->C_cols);
        Left->C_d->rows = Left->C_mp->rows = Right->C_rows;
        Left->C_d->cols = Left->C_mp->cols = Right->C_cols;
        for (i = 0; i < Right->C_rows; i++)
          for (j = 0; j < Right->C_cols; j++)
          { // setup real part
            mpq_set_str(Left->C_rat[i][j][0], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][0]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->C_rat[i][j][1], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->C_rat[i][j][1]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

            // setup C_d & C_mp
            rat_to_d(&Left->C_d->entry[i][j], Left->C_rat[i][j]);
            rat_to_mp(&Left->C_mp->entry[i][j], Left->C_rat[i][j]);
          }
      }

      // setup H
      init_vec_d(Left->H_d, Right->H_size);
      init_vec_mp2(Left->H_mp, Right->H_size, Left->curr_precision);
      init_vec_rat(Left->H_rat, Right->H_size);
      Left->H_d->size = Left->H_mp->size = Right->H_size;
      for (i = 0; i < Right->H_size; i++)
      { // setup real part
        mpq_set_str(Left->H_rat[i][0], &(*rpdStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][0]);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->H_rat[i][1], &(*rpdStr)[strLoc], base);
        mpq_canonicalize(Left->H_rat[i][1]);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

        // setup H_d & H_mp
        rat_to_d(&Left->H_d->coord[i], Left->H_rat[i]);
        rat_to_mp(&Left->H_mp->coord[i], Left->H_rat[i]);
      }

      // setup homVarConst
      init_d(Left->homVarConst_d);
      init_mp2(Left->homVarConst_mp, Left->curr_precision);
      init_rat(Left->homVarConst_rat);
      // setup real part
      mpq_set_str(Left->homVarConst_rat[0], &(*rpdStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[0]);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      // setup imag part
      mpq_set_str(Left->homVarConst_rat[1], &(*rpdStr)[strLoc], base);
      mpq_canonicalize(Left->homVarConst_rat[1]);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

      // setup homVarConst_d & homVarConst_mp
      rat_to_d(Left->homVarConst_d, Left->homVarConst_rat);
      rat_to_mp(Left->homVarConst_mp, Left->homVarConst_rat);

      // setup gamma
      init_d(Left->gamma_d);
      init_mp2(Left->gamma_mp, Left->curr_precision);
      init_rat(Left->gamma_rat);
      // setup real part
      mpq_set_str(Left->gamma_rat[0], &(*rpdStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[0]);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
      // setup imag part
      mpq_set_str(Left->gamma_rat[1], &(*rpdStr)[strLoc], base);
      mpq_canonicalize(Left->gamma_rat[1]);
      strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

      // setup gamma_d & gamma_mp
      rat_to_d(Left->gamma_d, Left->gamma_rat);
      rat_to_mp(Left->gamma_mp, Left->gamma_rat);

      // setup patchCoeff
      init_vec_d(Left->patchCoeff_d, Right->patchCoeff_size);
      init_vec_mp2(Left->patchCoeff_mp, Right->patchCoeff_size, Left->curr_precision);
      init_vec_rat(Left->patchCoeff_rat, Right->patchCoeff_size);
      Left->patchCoeff_d->size = Left->patchCoeff_mp->size = Right->patchCoeff_size;
      for (i = 0; i < Right->patchCoeff_size; i++)
      { // setup real part
        mpq_set_str(Left->patchCoeff_rat[i][0], &(*rpdStr)[strLoc], base);
        mpq_canonicalize(Left->patchCoeff_rat[i][0]);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
        // setup imag part
        mpq_set_str(Left->patchCoeff_rat[i][1], &(*rpdStr)[strLoc], base);
        mpq_canonicalize(Left->patchCoeff_rat[i][1]);
        strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

        // setup patchCoeff_d & patchCoeff_mp
        rat_to_d(&Left->patchCoeff_d->coord[i], Left->patchCoeff_rat[i]);
        rat_to_mp(&Left->patchCoeff_mp->coord[i], Left->patchCoeff_rat[i]);
      }

      // setup A
      Left->A_d = (mat_d *)bmalloc(Right->num_codim * sizeof(mat_d));
      Left->A_mp = (mat_mp *)bmalloc(Right->num_codim * sizeof(mat_mp));
      Left->A_rat = (mpq_t ****)bmalloc(Right->num_codim * sizeof(mpq_t ***));
      for (i = 0; i < Right->num_codim; i++)
      { // setup A[i]
        init_mat_d(Left->A_d[i], i + 1, Right->num_funcs);
        init_mat_mp2(Left->A_mp[i], i + 1, Right->num_funcs, Left->curr_precision);
        init_mat_rat(Left->A_rat[i], i + 1, Right->num_funcs);
        Left->A_d[i]->rows = Left->A_mp[i]->rows = i + 1;
        Left->A_d[i]->cols = Left->A_mp[i]->cols = Right->num_funcs;
        for (j = 0; j <= i; j++)
          for (k = 0; k < Right->num_funcs; k++)
          { // setup real part
            mpq_set_str(Left->A_rat[i][j][k][0], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->A_rat[i][j][k][0]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->A_rat[i][j][k][1], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->A_rat[i][j][k][1]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

            // setup A_d & A_mp
            rat_to_d(&Left->A_d[i]->entry[j][k], Left->A_rat[i][j][k]);
            rat_to_mp(&Left->A_mp[i]->entry[j][k], Left->A_rat[i][j][k]);
          }
      }

      // setup coeff_mp
      Left->coeff_d = (comp_d ***)bmalloc(Right->num_codim * sizeof(comp_mp **));
      Left->coeff_mp = (comp_mp ***)bmalloc(Right->num_codim * sizeof(comp_mp **));
      Left->coeff_rat = (mpq_t ****)bmalloc(Right->num_codim * sizeof(mpq_t ***));
      for (i = 0; i < Right->num_codim; i++)
      {
        Left->coeff_d[i] = (comp_d **)bmalloc(Left->new_degrees[i] * sizeof(comp_d *));
        Left->coeff_mp[i] = (comp_mp **)bmalloc(Left->new_degrees[i] * sizeof(comp_mp *));
        Left->coeff_rat[i] = (mpq_t ***)bmalloc(Left->new_degrees[i] * sizeof(mpq_t **));
        for (j = 0; j < Left->new_degrees[i]; j++)
        {
          Left->coeff_d[i][j] = (comp_d *)bmalloc(Left->new_variables * sizeof(comp_d));
          Left->coeff_mp[i][j] = (comp_mp *)bmalloc(Left->new_variables * sizeof(comp_mp));
          Left->coeff_rat[i][j] = (mpq_t **)bmalloc(Left->new_variables * sizeof(mpq_t *));
          for (k = 0; k < Left->new_variables; k++)
          {
            Left->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            init_d(Left->coeff_d[i][j][k]);
            init_mp2(Left->coeff_mp[i][j][k], Left->curr_precision);
            init_rat(Left->coeff_rat[i][j][k]);
            // setup real part
            mpq_set_str(Left->coeff_rat[i][j][k][0], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->coeff_rat[i][j][k][0]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->coeff_rat[i][j][k][1], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->coeff_rat[i][j][k][1]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

            // setup coeff_d & coeff_mp
            rat_to_d(Left->coeff_d[i][j][k], Left->coeff_rat[i][j][k]);
            rat_to_mp(Left->coeff_mp[i][j][k], Left->coeff_rat[i][j][k]);
          }
        }
      }
    }

    // free str
    if (freeStr)
      free(*rpdStr);
  }

  return;
}

void cp_regenCodim_int(void *Out, void *In, int MPType, char **rpdStr, comp_d **coeff_d, int curr_prec, int freeStr, int inType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: inType: 0 - Out is _t_int, In is CD                *
*                otherwise - Out is CD, In is _t_int            *
* RETURN VALUES:                                                *
* NOTES: stores a copy of In to Out                             *
\***************************************************************/
{
  int i, j, count = 0, strLoc = 0, base = 10;

  if (inType == 0)
  { // _t to _t_int so that it can be sent
    regenCodim_t *Right = (regenCodim_t *)In;
    regenCodim_t_int *Left = (regenCodim_t_int *)Out;

    // setup values
    Left->codim = Right->codim;
    Left->useIntrinsicSlice = Right->useIntrinsicSlice;

    // setup other values
    *rpdStr = NULL;
    *coeff_d = NULL;
    if (Right->useIntrinsicSlice)
    { // copy the values
      if (MPType == 0)
      { // copy to coeff_d
        Left->B_rows = Right->B_d->rows;
        Left->B_cols = Right->B_d->cols;
        Left->p_size = Right->p_d->size;

        Left->num_comp_d = Left->B_rows * Left->B_cols * Left->p_size;
        *coeff_d = (comp_d *)bmalloc(Left->num_comp_d * sizeof(comp_d));
        for (i = 0; i < Left->B_rows; i++)
          for (j = 0; j < Left->B_cols; j++)
          {
            set_d((*coeff_d)[count], &Right->B_d->entry[i][j]);
            count++;
          }
        for (i = 0; i < Left->p_size; i++)
        {
          set_d((*coeff_d)[count], &Right->p_d->coord[i]);
          count++;
        }

        // set length to 0
        Left->totalLength = 0;
      }
      else if (MPType == 1)
      { // setup rpdStr from _mp
        int sizeR, sizeI;
        char *strR = NULL, *strI = NULL;
 
        Left->B_rows = Right->B_mp->rows;
        Left->B_cols = Right->B_mp->cols;
        Left->p_size = Right->p_mp->size;

        for (i = 0; i < Left->B_rows; i++)
          for (j = 0; j < Left->B_cols; j++)
          { // real part
            strR = mpf_to_str(Right->B_mp->entry[i][j].r, base);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpf_to_str(Right->B_mp->entry[i][j].i, base);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        for (i = 0; i < Left->p_size; i++)
        { // real part
          strR = mpf_to_str(Right->p_mp->coord[i].r, base);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpf_to_str(Right->p_mp->coord[i].i, base);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update rpdStr
          *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*rpdStr)[strLoc], strR);
          strcpy(&(*rpdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

        // store the size
        Left->totalLength = strLoc;
        Left->num_comp_d = 0;
      }
      else
      { // setup rpdStr from _rat
        int sizeR, sizeI;
        char *strR = NULL, *strI = NULL;

        Left->B_rows = Right->B_d->rows;
        Left->B_cols = Right->B_d->cols;
        Left->p_size = Right->p_d->size;

        for (i = 0; i < Left->B_rows; i++)
          for (j = 0; j < Left->B_cols; j++)
          { // real part
            strR = mpq_get_str(NULL, base, Right->B_rat[i][j][0]);
            sizeR = strlen(strR) + 1; // +1 for '\0'
            // imag part
            strI = mpq_get_str(NULL, base, Right->B_rat[i][j][1]);
            sizeI = strlen(strI) + 1; // +1 for '\0'

            // update rpdStr
            *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
            strcpy(&(*rpdStr)[strLoc], strR);
            strcpy(&(*rpdStr)[strLoc + sizeR], strI);
            // update strLoc
            strLoc += sizeR + sizeI;

            // free strR, strI
            free(strR);
            free(strI);
          }
        for (i = 0; i < Left->p_size; i++)
        { // real part
          strR = mpq_get_str(NULL, base, Right->p_rat[i][0]);
          sizeR = strlen(strR) + 1; // +1 for '\0'
          // imag part
          strI = mpq_get_str(NULL, base, Right->p_rat[i][1]);
          sizeI = strlen(strI) + 1; // +1 for '\0'

          // update rpdStr
          *rpdStr = (char *)brealloc(*rpdStr, (strLoc + sizeR + sizeI) * sizeof(char));
          strcpy(&(*rpdStr)[strLoc], strR);
          strcpy(&(*rpdStr)[strLoc + sizeR], strI);
          // update strLoc
          strLoc += sizeR + sizeI;

          // free strR, strI
          free(strR);
          free(strI);
        }

        // store the size
        Left->totalLength = strLoc;
        Left->num_comp_d = 0;
      }
    }
    else
    { // set to 0
      Left->B_rows = Left->B_cols = Left->p_size = Left->num_comp_d = Left->totalLength = 0;
    }
  }
  else
  { // _t_int to _t so that it can be used
    regenCodim_t_int *Right = (regenCodim_t_int *)In;
    regenCodim_t *Left = (regenCodim_t *)Out;

    // setup values
    Left->codim = Right->codim;
    Left->useIntrinsicSlice = Right->useIntrinsicSlice;

    // setup other values
    if (Right->useIntrinsicSlice)
    { // copy the values
      if (MPType == 0)
      { // setup _d
        init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
        init_vec_d(Left->p_d, Right->p_size);
        Left->B_d->rows = Right->B_rows;
        Left->B_d->cols = Right->B_cols;
        Left->p_d->size = Right->p_size;

        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          {
            set_d(&Left->B_d->entry[i][j], (*coeff_d)[count]);
            count++;
          }
        for (i = 0; i < Right->p_size; i++)
        {
          set_d(&Left->p_d->coord[i], (*coeff_d)[count]);
          count++;
        }

        // clear coeff_d
        free(*coeff_d);
      }
      else if (MPType == 1)
      { // setup _mp
        init_mat_mp(Left->B_mp, Right->B_rows, Right->B_cols);
        init_vec_mp(Left->p_mp, Right->p_size);
        Left->B_mp->rows = Right->B_rows;
        Left->B_mp->cols = Right->B_cols;
        Left->p_mp->size = Right->p_size;

        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          { // setup real part
            mpf_set_str(Left->B_mp->entry[i][j].r, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpf_set_str(Left->B_mp->entry[i][j].i, &(*rpdStr)[strLoc], base);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          }
        for (i = 0; i < Right->p_size; i++)
        { // setup real part
          mpf_set_str(Left->p_mp->coord[i].r, &(*rpdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          // setup imag part
          mpf_set_str(Left->p_mp->coord[i].i, &(*rpdStr)[strLoc], base);
          strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
        }
      }
      else // MPType == 2
      { // setup _d, _mp, _rat
        init_mat_d(Left->B_d, Right->B_rows, Right->B_cols);
        init_mat_mp2(Left->B_mp, Right->B_rows, Right->B_cols, curr_prec);
        init_mat_rat(Left->B_rat, Right->B_rows, Right->B_cols);
        init_vec_d(Left->p_d, Right->p_size);
        init_vec_mp2(Left->p_mp, Right->p_size, curr_prec);
        init_vec_rat(Left->p_rat, Right->p_size);
        Left->B_d->rows = Left->B_mp->rows = Right->B_rows;
        Left->B_d->cols = Left->B_mp->cols = Right->B_cols;
        Left->p_d->size = Left->p_mp->size = Right->p_size;

        for (i = 0; i < Right->B_rows; i++)
          for (j = 0; j < Right->B_cols; j++)
          { // setup real part
            mpq_set_str(Left->B_rat[i][j][0], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->B_rat[i][j][0]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
            // setup imag part
            mpq_set_str(Left->B_rat[i][j][1], &(*rpdStr)[strLoc], base);
            mpq_canonicalize(Left->B_rat[i][j][1]);
            strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

            // setup B_d & B_mp
            rat_to_d(&Left->B_d->entry[i][j], Left->B_rat[i][j]);
            rat_to_mp(&Left->B_mp->entry[i][j], Left->B_rat[i][j]);
          }
        for (i = 0; i < Right->p_size; i++)
        { // setup real part
          mpq_set_str(Left->p_rat[i][0], &(*rpdStr)[strLoc], base);
          mpq_canonicalize(Left->p_rat[i][0]);
          strLoc += 1 + strlen(&(*rpdStr)[strLoc]);
          // setup imag part
          mpq_set_str(Left->p_rat[i][1], &(*rpdStr)[strLoc], base);
          mpq_canonicalize(Left->p_rat[i][1]);
          strLoc += 1 + strlen(&(*rpdStr)[strLoc]);

          // setup p_d & p_mp
          rat_to_d(&Left->p_d->coord[i], Left->p_rat[i]);
          rat_to_mp(&Left->p_mp->coord[i], Left->p_rat[i]);
        }
      }

      // clear rpdStr
      if (freeStr)
        free(*rpdStr);
    }
    else
    { // NULL out
      Left->B_rat = NULL;
      Left->p_rat = NULL;
    }
  }

  return;
}

#endif

