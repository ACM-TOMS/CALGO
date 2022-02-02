// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "regeneration.h"

void tracker_config_clear(tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears the allocated memory of T                       *
\***************************************************************/
{
  if (T->MPType == 1 || T->MPType == 2)
  { // clear latest_newton_residual_mp
    mpf_clear(T->latest_newton_residual_mp);
  }

  return;
}

void preproc_data_clear(preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear PPD                                              *
\***************************************************************/
{
  free(PPD->type);
  free(PPD->size);

  return;
}

void patch_eval_data_clear_d(patch_eval_data_d *PED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear PED                                              *
\***************************************************************/
{
  clear_mat_d(PED->patchCoeff);

  return;
}

void patch_eval_data_clear_mp(patch_eval_data_mp *PED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear PED                                              *
\***************************************************************/
{
  clear_mat_rat(PED->patchCoeff_rat, PED->patchCoeff->rows, PED->patchCoeff->cols);
  clear_mat_mp(PED->patchCoeff);

  return;
}

void start_system_eval_data_clear_d(start_system_eval_data_d *SSED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear SSED                                             *
\***************************************************************/
{
  int i, j, total_deg = 0;
  
  if (SSED->startSystemType == 1)
  { // find the total degree   
    for (i = 0; i < SSED->size_r; i++)
      total_deg += SSED->degrees[i];

    // clear the coeff
    for (i = total_deg - 1; i >= 0; i--)
    {
      for (j = SSED->coeff_cols - 1; j >= 0; j--)
      {
        clear_d(SSED->coeff[i][j]);
      }
      free(SSED->coeff[i]);
    }
    free(SSED->coeff);   
  }
  free(SSED->degrees);
  clear_d(SSED->gamma);

  return;
}

void start_system_eval_data_clear_mp(start_system_eval_data_mp *SSED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear SSED                                             *
\***************************************************************/
{
  int i, j, total_deg = 0;

  if (SSED->startSystemType == 1)
  { // find the total degree
    for (i = 0; i < SSED->size_r; i++)
      total_deg += SSED->degrees[i];

    // clear the coeff
    for (i = total_deg - 1; i >= 0; i--)
    {
      for (j = SSED->coeff_cols - 1; j >= 0; j--)
      {
        mpq_clear(SSED->coeff_rat[i][j][0]);
        mpq_clear(SSED->coeff_rat[i][j][1]);

        clear_mp(SSED->coeff[i][j]);

        free(SSED->coeff_rat[i][j]);
      }
      free(SSED->coeff[i]);
      free(SSED->coeff_rat[i]);
    }
    free(SSED->coeff);
    free(SSED->coeff_rat);
  }
  free(SSED->degrees);
  clear_mp(SSED->gamma);
  clear_rat(SSED->gamma_rat);
  free(SSED->gamma_rat);

  return;
}

void square_system_eval_data_clear_d(square_system_eval_data_d *SSED, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear SSED                                             *
\***************************************************************/
{
  int i;

  clearProg(SSED->Prog, MPType, 1);

  for (i = SSED->size_r - 1; i >= 0; i--)
    free(SSED->W[i]);
  free(SSED->W);
  free(SSED->P);
  free(SSED->orig_degrees);
  free(SSED->new_degrees);
  clear_mat_d(SSED->B);
  clear_mat_d(SSED->B_perp);
  clear_mat_d(SSED->A);  

  return;
}

void square_system_eval_data_clear_mp(square_system_eval_data_mp *SSED, int clrProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear SSED                                             *
\***************************************************************/
{
  int i;

  if (clrProg)
    clearProg(SSED->Prog, 1, 1); // MPType is 1 if it needs to be cleared

  for (i = SSED->size_r - 1; i >= 0; i--)
    free(SSED->W[i]);
  free(SSED->W);
  free(SSED->P);
  free(SSED->orig_degrees);
  free(SSED->new_degrees);

  clear_mat_rat(SSED->B_rat, SSED->B->rows, SSED->B->cols);
  clear_mat_mp(SSED->B);

  clear_mat_rat(SSED->B_perp_rat, SSED->B_perp->rows, SSED->B_perp->cols);
  clear_mat_mp(SSED->B_perp);

  clear_mat_rat(SSED->A_rat, SSED->A->rows, SSED->A->cols);
  clear_mat_mp(SSED->A);

  return;
}

void basic_eval_clear_mp(basic_eval_data_mp *ED, int clearRegen, int clrProg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear ED                                               *
\***************************************************************/
{
  square_system_eval_data_clear_mp(&ED->squareSystem, clrProg);
  start_system_eval_data_clear_mp(&ED->startSystem);
  patch_eval_data_clear_mp(&ED->patch);
  preproc_data_clear(&ED->preProcData);

  if (clearRegen == -59)
  {
    eqbyeq_clear_mp(ED->EqD);
  }

  return;
}
  
void basic_eval_clear_d(basic_eval_data_d *ED, int clearRegen, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear ED                                               *
\***************************************************************/
{
  square_system_eval_data_clear_d(&ED->squareSystem, MPType);
  start_system_eval_data_clear_d(&ED->startSystem);
  patch_eval_data_clear_d(&ED->patch);
  preproc_data_clear(&ED->preProcData);

  if (MPType == 2)
  { // ED->BED_mp->RD is just a pointer to ED->RD, so it will be cleared below - and thus 0 for not regen and 0 for not clear Prog since this was done above
    basic_eval_clear_mp(ED->BED_mp, 0, 0);
  }

  if (clearRegen == -59)
  {
    eqbyeq_clear_d(ED->EqD, MPType);
  }

  return;
}

///////////// EQ-BY-EQ CLEARING FUNCTIONS ////////////////////////

void eqbyeq_clear_mp(eqData_t *EqD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear EqD                                              *
\***************************************************************/
{
  int i, j;

  // clear only the top stage - all other cleared during eq-by-eq
  if (EqD->num_subsystems > 1)
    clearEqbyEqStageData_mp(EqD, EqD->num_subsystems - 1);
  free(EqD->stageData_mp);
  free(EqD->witnessData_mp);

  // clear gamma
  clear_mp(EqD->gamma_mp);

  // clear coeff_mp & degrees
  for (i = EqD->num_funcs - 1; i >= 0; i--)
  {
    for (j = EqD->num_vars - 1; j >= 0; j--)
    {
      clear_mp(EqD->coeff_mp[i][j]);
    }
    free(EqD->coeff_mp[i]);
    free(EqD->degrees[i]);
  }
  free(EqD->coeff_mp);
  free(EqD->degrees);

  free(EqD->startSub);
  free(EqD->endSub);
  free(EqD->startFunc);
  free(EqD->endFunc);
  free(EqD->startJvsub);
  free(EqD->endJvsub);
  free(EqD->startJv);
  free(EqD->endJv);

  if (EqD->numSubFuncs > 0)
  { // clear subFuncsBelow
    for (j = EqD->num_funcs - 1; j >= 0; j--)
      free(EqD->subFuncsBelow[j]);
    free(EqD->subFuncsBelow);
  }

  return;
}

void eqbyeq_clear_d(eqData_t *EqD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear EqD                                              *
\***************************************************************/
{
  int i, j;

  // clear only the top stage - all other cleared during eq-by-eq
  if (EqD->num_subsystems > 1)
    clearEqbyEqStageData_d(EqD, EqD->num_subsystems - 1, MPType);

  free(EqD->stageData_d);
  free(EqD->witnessData_d);

  // clear gamma
  clear_d(EqD->gamma_d);

  if (MPType == 2)
  { // clear the stages
    free(EqD->stageData_mp);
    free(EqD->witnessData_mp);

    // clear gamma
    clear_rat(EqD->gamma_rat);
    clear_mp(EqD->gamma_mp);

    // clear coeff_mp & coeff_rat
    for (i = EqD->num_funcs - 1; i >= 0; i--)
    {
      for (j = EqD->num_vars - 1; j >= 0; j--)
      {
        mpq_clear(EqD->coeff_rat[i][j][0]);
        mpq_clear(EqD->coeff_rat[i][j][1]);
        free(EqD->coeff_rat[i][j]);
        clear_mp(EqD->coeff_mp[i][j]);
      }
      free(EqD->coeff_rat[i]);
      free(EqD->coeff_mp[i]);
    }
    free(EqD->coeff_rat);
    free(EqD->coeff_mp);
  }

  // clear coeff_d & degrees
  for (i = EqD->num_funcs - 1; i >= 0; i--)
  {
    for (j = EqD->num_vars - 1; j >= 0; j--)
    {
      clear_d(EqD->coeff_d[i][j]);
    }
    free(EqD->coeff_d[i]);
    free(EqD->degrees[i]);
  }
  free(EqD->coeff_d);
  free(EqD->degrees);

  free(EqD->startSub);
  free(EqD->endSub);
  free(EqD->startFunc);
  free(EqD->endFunc);
  free(EqD->startJvsub);
  free(EqD->endJvsub);
  free(EqD->startJv);
  free(EqD->endJv);

  if (EqD->numSubFuncs > 0)
  { // clear subFuncsBelow
    for (j = EqD->num_funcs - 1; j >= 0; j--)
      free(EqD->subFuncsBelow[j]);
    free(EqD->subFuncsBelow);
  }

  return;
}

