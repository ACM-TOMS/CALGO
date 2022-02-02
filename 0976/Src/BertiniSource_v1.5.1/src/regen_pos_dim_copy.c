// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"

// COPY RPD TO witnessSuperset & clear it out

void regen_pos_dim_clearCodim(regenCodim_t *RC, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear RC                                               *
\***************************************************************/
{ // clear B & p
  if (RC->useIntrinsicSlice)
  { // see what is setup
    clear_mat(RC->B_d, RC->B_mp, RC->B_rat, MPType);
    clear_vec(RC->p_d, RC->p_mp, RC->p_rat, MPType);
  }

  // clear out other data
  RC->num_paths = RC->num_superset = RC->num_nonsing = RC->num_sing = RC->num_nonsolns = RC->num_inf = RC->num_bad = RC->useIntrinsicSlice = 0;  

  return;
}

void regen_pos_dim_copyWitness_clear_codim(witnessCodim_t *witCodim, regen_pos_dim_t *RPD, int codim_index, FILE *WITPTS, int MPType, int curr_prec, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy RPD[codim_index] to witCodim & clear it           *
\***************************************************************/
{
  int i, j, codim = RPD->codim[codim_index].codim;

  // codim, num_set, num_nonsing, num_sing
  witCodim->codim = RPD->codim[codim_index].codim;
  witCodim->num_set = RPD->codim[codim_index].num_superset;
  witCodim->num_nonsing = RPD->codim[codim_index].num_nonsing;
  witCodim->num_sing = RPD->codim[codim_index].num_sing;

  // setup W
  witCodim->W = (int **)bmalloc(codim * sizeof(int *));
  for (i = 0; i < codim; i++)
  {
    witCodim->W[i] = (int *)bmalloc((RPD->num_funcs - codim) * sizeof(int));
    for (j = 0; j < RPD->num_funcs - codim; j++)
      witCodim->W[i][j] = RPD->new_degrees[i] - RPD->new_degrees[j + codim];
  }

  // allocate witnessPt_types
  witCodim->witnessPt_types = (int *)bmalloc(witCodim->num_set * sizeof(int));

  // setup A
  witCodim->A_rows = codim;
  witCodim->A_cols = RPD->num_funcs - codim;
  if (MPType == 0)
  { // setup A_d & A_rat (random)
    init_mat_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);
    make_matrix_random_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);

    // setup A_rat
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpq_set_d(witCodim->A_rat[i][j][0], witCodim->A_d->entry[i][j].r);
        mpq_set_d(witCodim->A_rat[i][j][1], witCodim->A_d->entry[i][j].i);
      } 
  }
  else if (MPType == 1)
  { // seutp A_mp & A_rat (random)
    init_mat_mp2(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);
    make_matrix_random_mp(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);

    // setup A_rat
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    for (i = 0; i < witCodim->A_rows; i++)
      for (j = 0; j < witCodim->A_cols; j++)
      {
        mpf_t_to_rat(witCodim->A_rat[i][j][0], witCodim->A_mp->entry[i][j].r);
        mpf_t_to_rat(witCodim->A_rat[i][j][1], witCodim->A_mp->entry[i][j].i);
      }
  }
  else
  { // setup A_d, A_mp & A_rat
    // allocate for A_d, A_mp, A_rat
    init_mat_d(witCodim->A_d, witCodim->A_rows, witCodim->A_cols);
    init_mat_mp2(witCodim->A_mp, witCodim->A_rows, witCodim->A_cols, curr_prec);
    init_mat_rat(witCodim->A_rat, witCodim->A_rows, witCodim->A_cols);
    // setup A_d, A_mp, A_rat
    make_matrix_random_rat(witCodim->A_d, witCodim->A_mp, witCodim->A_rat, witCodim->A_rows, witCodim->A_cols, curr_prec, max_prec, 0, 0);
  }

  // copy over H & homVarConst
  if (MPType == 0 || MPType == 2)
  { // setup H_d & homVarConst_d
    init_vec_d(witCodim->H_d, RPD->H_d->size);
    vec_cp_d(witCodim->H_d, RPD->H_d);
    set_d(witCodim->homVarConst_d, RPD->homVarConst_d);
  }
  if (MPType == 1 || MPType == 2)
  { // setup H_mp & homVarConst_mp
    init_vec_mp2(witCodim->H_mp, RPD->H_mp->size, curr_prec);
    init_mp2(witCodim->homVarConst_mp, curr_prec);
    vec_cp_mp(witCodim->H_mp, RPD->H_mp);
    set_mp(witCodim->homVarConst_mp, RPD->homVarConst_mp);
  }
  if (MPType == 2)
  { // setup H_rat & homVarConst_rat

    // allocate, initialize and set H_rat
    init_vec_rat(witCodim->H_rat, witCodim->H_mp->size);
    for (i = 0; i < witCodim->H_mp->size; i++)
    {
      set_rat(witCodim->H_rat[i], RPD->H_rat[i]);
    }

    // initialize and set homVarConst_rat
    init_rat(witCodim->homVarConst_rat);
    set_rat(witCodim->homVarConst_rat, RPD->homVarConst_rat);
  }

  // setup B - coefficients for linear slices & p - patch coefficients
  if (MPType == 0 || MPType == 2)
  { // setup B_d & p_d
    init_mat_d(witCodim->B_d, RPD->new_variables - codim - 1, RPD->new_variables);
    init_vec_d(witCodim->p_d, RPD->new_variables);
    witCodim->B_d->rows = RPD->new_variables - codim - 1;
    witCodim->B_d->cols = witCodim->p_d->size = RPD->new_variables;

    for (j = 0; j < RPD->new_variables; j++)
    { // setup p[i]
      set_d(&witCodim->p_d->coord[j], &RPD->patchCoeff_d->coord[j]);
      for (i = 0; i < witCodim->B_d->rows; i++)
      { // setup B[i][j]
        set_d(&witCodim->B_d->entry[i][j], RPD->coeff_d[i + codim][0][j]); 
      }
    }
  }
  if (MPType == 1 || MPType == 2)
  { // setup B_mp & p_mp
    init_mat_mp(witCodim->B_mp, RPD->new_variables - codim - 1, RPD->new_variables);
    init_vec_mp(witCodim->p_mp, RPD->new_variables);
    witCodim->B_mp->rows = RPD->new_variables - codim - 1;
    witCodim->B_mp->cols = witCodim->p_mp->size = RPD->new_variables;
 
    for (j = 0; j < RPD->new_variables; j++)
    { // setup p[i]
      set_mp(&witCodim->p_mp->coord[j], &RPD->patchCoeff_mp->coord[j]);
      for (i = 0; i < witCodim->B_mp->rows; i++)
      { // setup B[i][j]
        set_mp(&witCodim->B_mp->entry[i][j], RPD->coeff_mp[i + codim][0][j]);
      } 
    } 
  }
  if (MPType == 2)
  { // setup B_rat & p_rat
    init_mat_rat(witCodim->B_rat, RPD->new_variables - codim - 1, RPD->new_variables);
    init_vec_rat(witCodim->p_rat, RPD->new_variables);
  
    for (j = 0; j < RPD->new_variables; j++)
    { // setup p[i]
      set_rat(witCodim->p_rat[j], RPD->patchCoeff_rat[j]);
      for (i = 0; i < witCodim->B_mp->rows; i++) 
      { // setup B[i][j] 
        set_rat(witCodim->B_rat[i][j], RPD->coeff_rat[i + codim][0][j]);
      }  
    }  
  }  
  
  // setup witnessPts & witnessPt_types
  if (MPType == 0)
  { // allocate witnessPts_d
    witCodim->witnessPts_d = (endpoint_data_d *)bmalloc(witCodim->num_set * sizeof(endpoint_data_d));
    witCodim->witnessPts_mp = NULL;
    witCodim->witnessPts_amp = NULL;

    for (i = 0; i < witCodim->num_set; i++)
    { // setup the ith endpoint from WITPTS
      init_endpoint_data_d(&witCodim->witnessPts_d[i]);
      setup_endpoint_data_d(&witCodim->witnessPts_d[i], WITPTS);
      // setup type
      if (witCodim->witnessPts_d[i].corank > 0)
        witCodim->witnessPt_types[i] = SOLUTION_AND_SING;
      else
        witCodim->witnessPt_types[i] = SOLUTION_AND_NONSING;
    }    
  }  
  else if (MPType == 1)
  { // allocate witnessPts_mp
    witCodim->witnessPts_mp = (endpoint_data_mp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_mp));
    witCodim->witnessPts_d = NULL;
    witCodim->witnessPts_amp = NULL;

    for (i = 0; i < witCodim->num_set; i++)
    { // setup the ith endpoint from WITPTS
      init_endpoint_data_mp(&witCodim->witnessPts_mp[i]);
      setup_endpoint_data_mp(&witCodim->witnessPts_mp[i], WITPTS);
      // setup type
      if (witCodim->witnessPts_mp[i].corank > 0)
        witCodim->witnessPt_types[i] = SOLUTION_AND_SING;
      else
        witCodim->witnessPt_types[i] = SOLUTION_AND_NONSING;
    } 
  }
  else
  { // allocate witnessPts_amp
    witCodim->witnessPts_amp = (endpoint_data_amp *)bmalloc(witCodim->num_set * sizeof(endpoint_data_amp));
    witCodim->witnessPts_d = NULL;
    witCodim->witnessPts_mp = NULL;

    for (i = 0; i < witCodim->num_set; i++)
    { // setup the ith endpoint from WITPTS
      init_endpoint_data_amp(&witCodim->witnessPts_amp[i], 64, 64);
      setup_endpoint_data_amp(&witCodim->witnessPts_amp[i], WITPTS);
      // setup type
      if (witCodim->witnessPts_amp[i].corank > 0)
        witCodim->witnessPt_types[i] = SOLUTION_AND_SING;
      else
        witCodim->witnessPt_types[i] = SOLUTION_AND_NONSING;
    } 
  }

  // intialize other values in witCodim
  witCodim->num_components = 0;
  witCodim->multiplicities = witCodim->component_nums = witCodim->deflations_needed = NULL;

  // clear codim
  regen_pos_dim_clearCodim(&RPD->codim[codim_index], MPType);

  return;
}

int regen_pos_dim_copyWitness_clear(witness_t *witnessSuperset, regen_pos_dim_t *RPD, char *witName, int MPType, int max_prec, int specificCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if specificCodim is top dimensional            *
* NOTES: copy RPD to witnessSuperset and clear RPD              *
\***************************************************************/
{
  int i, j, k, count, retVal = 0;
  FILE *WITPTS = NULL;
  char *witnessPts = NULL;
  size_t size;

  // target slices have not been initialized
  witnessSuperset->targetSliceInit = 0;

  // Prog
  witnessSuperset->Prog = RPD->Prog;

  // PPD
  cp_preproc_data(&witnessSuperset->PPD, &RPD->PPD);

  // orig_degrees, new_degrees, P
  witnessSuperset->orig_degrees = RPD->orig_degrees;
  witnessSuperset->new_degrees = RPD->new_degrees;
  witnessSuperset->P = RPD->P;

  // system_rank, orig_variables, new_variables, curr_precision, num_funcs
  witnessSuperset->system_rank = RPD->system_rank;
  witnessSuperset->orig_variables = RPD->orig_variables;
  witnessSuperset->new_variables = RPD->new_variables;
  witnessSuperset->curr_precision = RPD->curr_precision;
  witnessSuperset->num_funcs = RPD->num_funcs;

  // C (if needed)
  if (RPD->orig_variables != RPD->new_variables)
  {
    if (MPType == 0 || MPType == 2)
    { // copy and clear C_d
      init_mat_d(witnessSuperset->C_d, RPD->C_d->rows, RPD->C_d->cols);
      mat_cp_d(witnessSuperset->C_d, RPD->C_d);
      clear_mat_d(RPD->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // copy and clear C_mp
      init_mat_mp2(witnessSuperset->C_mp, RPD->C_mp->rows, RPD->C_mp->cols, witnessSuperset->curr_precision);
      mat_cp_mp(witnessSuperset->C_mp, RPD->C_mp);
      clear_mat_mp(RPD->C_mp);
    }

    if (MPType == 2)
    { // copy and clear C_rat
      witnessSuperset->C_rat = RPD->C_rat;
      RPD->C_rat = NULL;
    }
  }

  // gamma
  if (MPType == 0 || MPType == 2)
  { // copy and clear gamma_d
    set_d(witnessSuperset->gamma_d, RPD->gamma_d);
    clear_d(RPD->gamma_d);
  }
  if (MPType == 1 || MPType == 2)
  { // copy and clear gamma_mp
    init_mp2(witnessSuperset->gamma_mp, witnessSuperset->curr_precision);
    set_mp(witnessSuperset->gamma_mp, RPD->gamma_mp);
    clear_mp(RPD->gamma_mp);
  }
  if (MPType == 2)
  { // copy and clear gamma_rat
    init_rat(witnessSuperset->gamma_rat);
    set_rat(witnessSuperset->gamma_rat, RPD->gamma_rat);
    clear_rat(RPD->gamma_rat);
  }

  // count the number of codim that have witness points
  count = 0;
  retVal = 1;
  for (i = 0; i < RPD->num_codim; i++)
    if (RPD->codim[i].num_superset > 0)
    {
      if (specificCodim == 0 || RPD->codim[i].codim == specificCodim)
        count++;
      else if (specificCodim > 0)
        retVal = 0;
    }


  // set num_codim to count
  witnessSuperset->num_codim = count;

  // allocate memory for the 'num_codim' witness sets
  witnessSuperset->codim = (witnessCodim_t *)bmalloc(witnessSuperset->num_codim * sizeof(witnessCodim_t));

  // codim
  count = 0;
  for (i = 0; i < RPD->num_codim; i++)
  { // setup the string for the name of the file containing the witness superset points
    size = 1 + snprintf(NULL, 0, "%s_%d", witName, RPD->codim[i].codim);
    witnessPts = (char *)brealloc(witnessPts, size * sizeof(char));
    sprintf(witnessPts, "%s_%d", witName, RPD->codim[i].codim);

    if ((specificCodim == 0 || RPD->codim[i].codim == specificCodim) && RPD->codim[i].num_superset > 0)
    { // open up the file containing the witness superset points for the codim
      WITPTS = fopen(witnessPts, "r");
      // make sure WITPTS exits
      if (WITPTS == NULL)
      {
        printf("ERROR: The file to contain the witness superset points, '%s', does not exist!\n", witnessPts);
        bexit(ERROR_FILE_NOT_EXIST);
      }

      // copy and clear the ith codim
      regen_pos_dim_copyWitness_clear_codim(&witnessSuperset->codim[count], RPD, i, WITPTS, MPType, witnessSuperset->curr_precision, max_prec);

      // close WITPTS
      fclose(WITPTS);
 
      // increment count
      count++;
    }
    else
    { // just clear this codim 
      regen_pos_dim_clearCodim(&RPD->codim[i], MPType);
    }
  }

  // clear W & A
  for (i = RPD->num_codim - 1; i >= 0; i--)
  { // clear W
    for (j = 0; j <= i; j++)
      free(RPD->W[i][j]);
    free(RPD->W[i]);

    // clear A
    if (MPType == 0)
    { // clear A_d
      clear_mat_d(RPD->A_d[i]);
    }
    else if (MPType == 1)
    { // clear A_mp
      clear_mat_mp(RPD->A_mp[i]);
    }
    else
    { // clear A_d, A_mp & A_rat
      clear_mat(RPD->A_d[i], RPD->A_mp[i], RPD->A_rat[i], MPType);
    }
  }
  free(RPD->W);
  if (MPType == 0)
    free(RPD->A_d);
  else if (MPType == 1)
    free(RPD->A_mp);
  else
    free(RPD->A_rat);
  RPD->A_d = NULL;
  RPD->A_mp = NULL;
  RPD->A_rat = NULL;

  // clear H, homVarConst, patchCoeff
  clear_vec(RPD->H_d, RPD->H_mp, RPD->H_rat, MPType);
  clear_d_mp_rat(RPD->homVarConst_d, RPD->homVarConst_mp, RPD->homVarConst_rat, MPType);
  clear_vec(RPD->patchCoeff_d, RPD->patchCoeff_mp, RPD->patchCoeff_rat, MPType);

  // clear coeff
  if (MPType == 0)
  { // clear coeff_d
    for (i = RPD->num_codim - 1; i >= 0; i--)
    { 
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
        free(RPD->coeff_d[i][j]);
      free(RPD->coeff_d[i]);
    }
    free(RPD->coeff_d);
  }
  else if (MPType == 1)
  { // clear coeff_mp
    for (i = RPD->num_codim - 1; i >= 0; i--)
    {
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
      {
        for (k = RPD->new_variables - 1; k >= 0; k--)
        {
          clear_mp(RPD->coeff_mp[i][j][k]);
        }
        free(RPD->coeff_mp[i][j]);
      }
      free(RPD->coeff_mp[i]);
    }
    free(RPD->coeff_mp);
  }
  else
  { // clear coeff_d, coeff_mp & coeff_rat
    for (i = RPD->num_codim - 1; i >= 0; i--)
    {
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
      {
        for (k = RPD->new_variables - 1; k >= 0; k--)
        {
          clear_d_mp_rat(RPD->coeff_d[i][j][k], RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k], MPType);
        }
        free(RPD->coeff_d[i][j]);
        free(RPD->coeff_mp[i][j]);
        free(RPD->coeff_rat[i][j]);
      }
      free(RPD->coeff_d[i]);
      free(RPD->coeff_mp[i]);
      free(RPD->coeff_rat[i]);
    }
    free(RPD->coeff_d);
    free(RPD->coeff_mp);
    free(RPD->coeff_rat);
  }
  RPD->coeff_d = NULL;
  RPD->coeff_mp = NULL;
  RPD->coeff_rat = NULL;

  // clear other data
  RPD->Prog = NULL;
  preproc_data_clear(&RPD->PPD);
  RPD->orig_degrees = RPD->new_degrees = RPD->P = NULL;
  free(RPD->codim);

  // clear memory
  free(witnessPts);

  return retVal;
}

void regen_pos_dim_clear(regen_pos_dim_t *RPD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear RPD                                              *
\***************************************************************/
{
  int i, j, k;

  // clear codim
  for (i = RPD->num_codim - 1; i >= 0; i--)
    regen_pos_dim_clearCodim(&RPD->codim[i], MPType);
  free(RPD->codim);

  // clear Prog
  clearProg(RPD->Prog, MPType, 0);
  free(RPD->Prog);

  // clear PPD
  preproc_data_clear(&RPD->PPD);

  // clear C (if needed)
  if (RPD->orig_variables != RPD->new_variables)
  {
    clear_mat(RPD->C_d, RPD->C_mp, RPD->C_rat, MPType);
  }

  // clear gamma
  clear_d_mp_rat(RPD->gamma_d, RPD->gamma_mp, RPD->gamma_rat, MPType);

  // clear H, homVarConst, patchCoeff
  clear_vec(RPD->H_d, RPD->H_mp, RPD->H_rat, MPType);
  clear_d_mp_rat(RPD->homVarConst_d, RPD->homVarConst_mp, RPD->homVarConst_rat, MPType);
  clear_vec(RPD->patchCoeff_d, RPD->patchCoeff_mp, RPD->patchCoeff_rat, MPType);

  // clear W & A
  for (i = RPD->num_codim - 1; i >= 0; i--)
  { // clear W
    for (j = 0; j <= i; j++)
      free(RPD->W[i][j]);
    free(RPD->W[i]);

    if (MPType == 0)
    { // clear A_d
      clear_mat_d(RPD->A_d[i]);
    }
    else if (MPType == 1)
    { // clear A_mp
      clear_mat_mp(RPD->A_mp[i]);
    }
    else
    { // clear A_d, A_mp & A_rat
      clear_mat(RPD->A_d[i], RPD->A_mp[i], RPD->A_rat[i], MPType);
    }
  }
  free(RPD->W);
  if (MPType == 0)
    free(RPD->A_d);
  else if (MPType == 1)
    free(RPD->A_mp);
  else
    free(RPD->A_rat);
  RPD->A_d = NULL;
  RPD->A_mp = NULL;
  RPD->A_rat = NULL;

  // clear coeff
  if (MPType == 0)
  { // clear coeff_d
    for (i = RPD->num_codim - 1; i >= 0; i--)
    {
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
        free(RPD->coeff_d[i][j]);
      free(RPD->coeff_d[i]);
    }
    free(RPD->coeff_d);
  }
  else if (MPType == 1)
  { // clear coeff_mp
    for (i = RPD->num_codim - 1; i >= 0; i--)
    {
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
      {
        for (k = RPD->new_variables - 1; k >= 0; k--)
        {
          clear_mp(RPD->coeff_mp[i][j][k]);
        }
        free(RPD->coeff_mp[i][j]);
      }
      free(RPD->coeff_mp[i]);
    }
    free(RPD->coeff_mp);
  }
  else
  { // clear coeff_d, coeff_mp & coeff_rat
    for (i = RPD->num_codim - 1; i >= 0; i--)
    {
      for (j = RPD->new_degrees[i] - 1; j >= 0; j--)
      {
        for (k = RPD->new_variables - 1; k >= 0; k--)
        {
          clear_d_mp_rat(RPD->coeff_d[i][j][k], RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k], MPType);
        }
        free(RPD->coeff_d[i][j]);
        free(RPD->coeff_mp[i][j]);
        free(RPD->coeff_rat[i][j]);
      }
      free(RPD->coeff_d[i]);
      free(RPD->coeff_mp[i]);
      free(RPD->coeff_rat[i]);
    }
    free(RPD->coeff_d);
    free(RPD->coeff_mp);
    free(RPD->coeff_rat);
  }
  RPD->coeff_d = NULL;
  RPD->coeff_mp = NULL;
  RPD->coeff_rat = NULL;

  // clear orig_degrees, new_degrees, P
  free(RPD->new_degrees);
  free(RPD->orig_degrees);
  free(RPD->P);

  return;
}












