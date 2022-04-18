// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"

void createFirstRegenCodim_startPoints(regen_pos_dim_t *RPD, int MPType, char *startFileName);
void writeRegenCodimStartPts(regen_pos_dim_t *RPD, int MPType, int count, FILE *START);
void setupRegenCodimIntrinsicSlice(regen_pos_dim_t *RPD, int MPType, int codim_index, int max_prec);
void setupRPDRestart(regen_pos_dim_t *RPD, int MPType, int max_prec, int codim_index, FILE *FP);

void regen_pos_dim_setup(int startCodim, int *maxCodim, int specificCodim, tracker_config_t *T, regen_pos_dim_t *RPD, char *preprocFile, char *degreeFile, char *startName, double intrinsicCutoffMultiplier)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for positive dimensional tracking - regeneration *
\***************************************************************/
{
  int i, j, k, intrinsicCutoff;
  size_t size;
  char *startFileName = NULL;

  // setup the startFileName
  k = startCodim <= 0 ? 1 : startCodim;
  size = 1 + snprintf(NULL, 0, "%s_%d", startName, k);
  startFileName = (char *)brealloc(startFileName, size * sizeof(char));
  sprintf(startFileName, "%s_%d", startName, k);

  // store the precision
  RPD->curr_precision = T->Precision;

  // setup the SLP
  RPD->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  RPD->orig_variables = setupProg(RPD->Prog, T->Precision, T->MPType);

  // error checking
  if (RPD->Prog->numPathVars > 0)
  { // path variable present
    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (RPD->Prog->numPars > 0)
  { // parameter present
    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup pre proc data
  setupPreProcData(preprocFile, &RPD->PPD);

  // verify that we are using only 1 homogenous variable group
  if (RPD->PPD.num_hom_var_gp + RPD->PPD.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: Positive dimensional regeneration is implemented for systems with only one variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find the rank
  if (T->MPType == 0 || T->MPType == 2)
    RPD->system_rank = rank_finder_d(&RPD->PPD, RPD->Prog, T, RPD->orig_variables);
  else
    RPD->system_rank = rank_finder_mp(&RPD->PPD, RPD->Prog, T, RPD->orig_variables);

  // setup the number of variables that will be used for tracking
  T->numVars = RPD->new_variables = RPD->system_rank + RPD->PPD.num_var_gp + RPD->PPD.num_hom_var_gp; // add on the number of patches attached - should be 1!

  // setup the number of functions and codimension
  RPD->num_funcs = RPD->PPD.num_funcs;
  RPD->num_codim = RPD->system_rank;
  *maxCodim = *maxCodim > 0 ? MIN(*maxCodim, RPD->system_rank) : RPD->system_rank; 

  // error checking on specific codimension
  if (specificCodim > 0 && specificCodim > RPD->system_rank)
  {
    printf("NOTE: Based on the computed system rank (%d), this system has no components of codimension %d!\n", RPD->system_rank, specificCodim);
    bexit(ERROR_INPUT_SYSTEM);
  }

  // determine where to switch from intrinsic to extrinsic
  intrinsicCutoff = floor(intrinsicCutoffMultiplier * RPD->new_variables);

  // setup orig_degrees, new_degrees & P
  setupDegrees_orig_new_perm(&RPD->orig_degrees, &RPD->new_degrees, &RPD->P, RPD->num_funcs, RPD->PPD.num_var_gp + RPD->PPD.num_hom_var_gp, degreeFile);

  // setup W
  RPD->W = (int ***)bmalloc(RPD->num_codim * sizeof(int **));
  for (i = 0; i < RPD->num_codim; i++)
  { // setup for codimension i + 1
    RPD->W[i] = (int **)bmalloc((i + 1) * sizeof(int *));
    for (j = 0; j <= i; j++)
    {
      RPD->W[i][j] = (int *)bmalloc((RPD->num_funcs - j - 1) * sizeof(int));
      for (k = 0; k < RPD->num_funcs - j - 1; k++)
      {
        RPD->W[i][j][k] = RPD->new_degrees[j] - RPD->new_degrees[k + j + 1];
      }
    }
  }

  // setup gamma
  if (T->MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(RPD->gamma_d);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    init_mp(RPD->gamma_mp);
    get_comp_rand_mp(RPD->gamma_mp);
  }
  else
  { // setup gamma_d, gamma_mp & gamma_rat
    get_comp_rand_rat(RPD->gamma_d, RPD->gamma_mp, RPD->gamma_rat, RPD->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // allocate codim
  RPD->codim = (regenCodim_t *)bmalloc(RPD->num_codim * sizeof(regenCodim_t));

  // initialize
  RPD->sameA = 0;

  if (startCodim <= 0)
  { // we are starting from scratch

    // setup C, if needed
    if (RPD->orig_variables != RPD->new_variables)
    { // we will need to convert between the old and new variables (orig_vars = C * new_vars, we take C = [[I];[RAND]])
      if (T->MPType == 0)
      { // only setup C_d
        mat_d tempMat_d;
        init_mat_d(tempMat_d, RPD->orig_variables - RPD->new_variables, RPD->new_variables);
        make_matrix_random_d(tempMat_d, RPD->orig_variables - RPD->new_variables, RPD->new_variables);

        // setup C_d 
        init_mat_d(RPD->C_d, RPD->orig_variables, RPD->new_variables);
        RPD->C_d->rows = RPD->orig_variables;
        RPD->C_d->cols = RPD->new_variables;
        // set the top to be Identity and bottom to be tempMat
        for (i = 0; i < RPD->orig_variables; i++)
          if (i < RPD->new_variables)
          { // Identity 
            for (j = 0; j < RPD->new_variables; j++)
              if (i == j)
              {
                set_one_d(&RPD->C_d->entry[i][j]);
              }
              else
              {
                set_zero_d(&RPD->C_d->entry[i][j]);
              }
          }
          else
          { // tempMat
            for (j = 0; j < RPD->new_variables; j++)
            {
              set_d(&RPD->C_d->entry[i][j], &tempMat_d->entry[i - RPD->new_variables][j]);
            }
          }

        clear_mat_d(tempMat_d);
      }
      else if (T->MPType == 1)
      { // only setup C_mp
        mat_mp tempMat_mp;
        init_mat_mp(tempMat_mp, RPD->orig_variables - RPD->new_variables, RPD->new_variables);
        make_matrix_random_mp(tempMat_mp, RPD->orig_variables - RPD->new_variables, RPD->new_variables, T->Precision);

        // setup C_mp
        init_mat_mp(RPD->C_mp, RPD->orig_variables, RPD->new_variables);
        RPD->C_mp->rows = RPD->orig_variables;
        RPD->C_mp->cols = RPD->new_variables;
        // set the top to be Identity and bottom to be tempMat
        for (i = 0; i < RPD->orig_variables; i++)
          if (i < RPD->new_variables)
          { // Identity
            for (j = 0; j < RPD->new_variables; j++)
              if (i == j)
              {
                set_one_mp(&RPD->C_mp->entry[i][j]);
              }
              else
              {
                set_zero_mp(&RPD->C_mp->entry[i][j]);
              }
          }
          else
          { // tempMat
            for (j = 0; j < RPD->new_variables; j++)
            {
              set_mp(&RPD->C_mp->entry[i][j], &tempMat_mp->entry[i - RPD->new_variables][j]);
            }
          }

        clear_mat_mp(tempMat_mp);
      }
      else
      { // setup C_d, C_mp & C_rat
        mat_d tempMat_d;
        mat_mp tempMat_mp;
        mpq_t ***tempMat_rat = NULL;
        init_mat_d(tempMat_d, RPD->orig_variables - RPD->new_variables, RPD->new_variables);
        init_mat_mp2(tempMat_mp, RPD->orig_variables - RPD->new_variables, RPD->new_variables, RPD->curr_precision);
        init_mat_rat(tempMat_rat, RPD->orig_variables - RPD->new_variables, RPD->new_variables);
        make_matrix_random_rat(tempMat_d, tempMat_mp, tempMat_rat, RPD->orig_variables - RPD->new_variables, RPD->new_variables, RPD->curr_precision, T->AMP_max_prec, 0, 0);

        // allocate for C_d, C_mp & C_rat
        init_mat_d(RPD->C_d, RPD->orig_variables, RPD->new_variables);
        init_mat_mp2(RPD->C_mp, RPD->orig_variables, RPD->new_variables, RPD->curr_precision);
        init_mat_rat(RPD->C_rat, RPD->orig_variables, RPD->new_variables);
        RPD->C_d->rows = RPD->C_mp->rows = RPD->orig_variables;
        RPD->C_d->cols = RPD->C_mp->cols = RPD->new_variables;
        // set the top to be Identity and bottom to be tempMat
        for (i = 0; i < RPD->orig_variables; i++)
          if (i < RPD->new_variables)
          { // Identity
            for (j = 0; j < RPD->new_variables; j++)
              if (i == j)
              {
                set_one_d(&RPD->C_d->entry[i][j]);
                set_one_mp(&RPD->C_mp->entry[i][j]);
                set_one_rat(RPD->C_rat[i][j]);
              }
              else
              {
                set_zero_d(&RPD->C_d->entry[i][j]);
                set_zero_mp(&RPD->C_mp->entry[i][j]);
                set_zero_rat(RPD->C_rat[i][j]);
              }
          }
          else
          { // tempMat
            for (j = 0; j < RPD->new_variables; j++)
            {
              set_d(&RPD->C_d->entry[i][j], &tempMat_d->entry[i - RPD->new_variables][j]);
              set_mp(&RPD->C_mp->entry[i][j], &tempMat_mp->entry[i - RPD->new_variables][j]);
              set_rat(RPD->C_rat[i][j], tempMat_rat[i - RPD->new_variables][j]);
            }
          }

        clear_mat(tempMat_d, tempMat_mp, tempMat_rat, T->MPType);
      }
    }
    else
    { // set to 0
      RPD->C_d->rows = RPD->C_d->cols = RPD->C_mp->rows = RPD->C_mp->cols = 0;
    }

    // setup H, homVarConst, patchCoeff, A
    if (T->MPType == 0)
    { // setup H_d & homVarConst_d
      init_vec_d(RPD->H_d, RPD->new_variables);
      RPD->H_d->size = RPD->new_variables;
      if (RPD->PPD.num_var_gp > 0)
      { // using a variable group that was originally not homogenized
        if (RPD->orig_variables != RPD->new_variables)
        { // H_d is first row of C
          for (i = 0; i < RPD->new_variables; i++)
          {
            set_d(&RPD->H_d->coord[i], &RPD->C_d->entry[0][i]);
          }
        }
        else
        { // H_d = [1,0..0]
          set_one_d(&RPD->H_d->coord[0]);
          for (i = 1; i < RPD->new_variables; i++)
          {
            set_zero_d(&RPD->H_d->coord[i]);
          }
        }

        // setup homVarConst_d to be 0
        set_zero_d(RPD->homVarConst_d);
      }
      else
      { // using a homogeneous variable group

        // setup H_d to be random
        make_vec_random_d(RPD->H_d, RPD->new_variables);

        // setup homVarConst_d to be random
        get_comp_rand_d(RPD->homVarConst_d);
      }

      // setup patchCoeff_d
      init_vec_d(RPD->patchCoeff_d, RPD->new_variables);
      make_vec_random_d(RPD->patchCoeff_d, RPD->new_variables);

      // setup A_d
      RPD->sameA = 1;
      RPD->A_d = (mat_d *)bmalloc(RPD->num_codim * sizeof(mat_d));
      // make them in reverse order
      for (i = RPD->num_codim - 1; i >= 0; i--)
      { // setup A_d[i]
        init_mat_d(RPD->A_d[i], i + 1, RPD->num_funcs); 
        RPD->A_d[i]->rows = i + 1;
        RPD->A_d[i]->cols = RPD->num_funcs;

        if (i == RPD->num_codim - 1)
        { // we generate the main matrix
          make_matrix_random_d(RPD->A_d[i], i + 1, RPD->num_funcs);
          // make upper triangular with 1's on diagonal
          for (j = 0; j <= i; j++)
            for (k = 0; k <= j; k++)
              if (j == k)
              { // set to 1
                set_one_d(&RPD->A_d[i]->entry[j][k]);
              }
              else
              { // set to 0
                set_zero_d(&RPD->A_d[i]->entry[j][k]);
              }
        }
        else
        { // copy the top of A_d[i+1]
          for (j = 0; j <= i; j++)
            for (k = 0; k < RPD->num_funcs; k++)
            {
              set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
            }
        }
      }
    }
    else if (T->MPType == 1)
    { // setup H_mp & homVarConst_mp
      init_vec_mp(RPD->H_mp, RPD->new_variables);
      RPD->H_mp->size = RPD->new_variables;
      init_mp(RPD->homVarConst_mp);
      if (RPD->PPD.num_var_gp > 0)
      { // using a variable group that was originally not homogenized
        if (RPD->orig_variables != RPD->new_variables)
        { // H_mp is first row of C
          for (i = 0; i < RPD->new_variables; i++)
          {
            set_mp(&RPD->H_mp->coord[i], &RPD->C_mp->entry[0][i]);
          }
        }
        else
        { // H_d = [1,0..0]
          set_one_mp(&RPD->H_mp->coord[0]);
          for (i = 1; i < RPD->new_variables; i++)
          {
            set_zero_mp(&RPD->H_mp->coord[i]);
          }
        }

        // setup homVarConst_mp to be 0
        set_zero_mp(RPD->homVarConst_mp);
      }
      else
      { // using a homogeneous variable group

        // setup H_mp to be random
        make_vec_random_mp(RPD->H_mp, RPD->new_variables);

        // setup homVarConst_mp to be random
        get_comp_rand_mp(RPD->homVarConst_mp);
      }

      // setup patchCoeff_mp
      init_vec_mp(RPD->patchCoeff_mp, RPD->new_variables);
      make_vec_random_mp(RPD->patchCoeff_mp, RPD->new_variables);

      // setup A_mp
      RPD->sameA = 1;
      RPD->A_mp = (mat_mp *)bmalloc(RPD->num_codim * sizeof(mat_mp));
      // make them in reverse order
      for (i = RPD->num_codim - 1; i >= 0; i--)
      { // setup A_mp[i]
        init_mat_mp(RPD->A_mp[i], i + 1, RPD->num_funcs);
        RPD->A_mp[i]->rows = i + 1;
        RPD->A_mp[i]->cols = RPD->num_funcs;

        if (i == RPD->num_codim - 1)
        { // we generate the main matrix
          make_matrix_random_mp(RPD->A_mp[i], i + 1, RPD->num_funcs, T->Precision);
          // make upper triangular with 1's on diagonal
          for (j = 0; j <= i; j++)
            for (k = 0; k <= j; k++)
              if (j == k)
              { // set to 1
                set_one_mp(&RPD->A_mp[i]->entry[j][k]);
              }
              else
              { // set to 0
                set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
              }
        }
        else
        { // copy the top of A_mp[i+1]
          for (j = 0; j <= i; j++)
            for (k = 0; k < RPD->num_funcs; k++)
            {
              set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
            }
        }
      }
    }
    else
    { // setup H_d, H_mp, H_rat & homVarConst_d, homVarConst_mp, homVarConst_rat
      init_vec_d(RPD->H_d, RPD->new_variables);
      init_vec_mp2(RPD->H_mp, RPD->new_variables, RPD->curr_precision);
      init_vec_rat(RPD->H_rat, RPD->new_variables);
      RPD->H_d->size = RPD->H_mp->size = RPD->new_variables;
      init_mp2(RPD->homVarConst_mp, RPD->curr_precision);
      init_rat(RPD->homVarConst_rat);

      if (RPD->PPD.num_var_gp > 0)
      { // using a variable group that was originally not homogenized
        if (RPD->orig_variables != RPD->new_variables)
        { // H is first row of C
          for (i = 0; i < RPD->new_variables; i++)
          {
            set_d(&RPD->H_d->coord[i], &RPD->C_d->entry[0][i]);
            set_mp(&RPD->H_mp->coord[i], &RPD->C_mp->entry[0][i]);
            set_rat(RPD->H_rat[i], RPD->C_rat[0][i]);
          }
        }
        else
        { // H_d = [1,0..0]
          set_one_d(&RPD->H_d->coord[0]);
          set_one_mp(&RPD->H_mp->coord[0]);
          set_one_rat(RPD->H_rat[0]);
          for (i = 1; i < RPD->new_variables; i++)
          {
            set_zero_d(&RPD->H_d->coord[i]);
            set_zero_mp(&RPD->H_mp->coord[i]);
            set_zero_rat(RPD->H_rat[i]);
          }
        }

        // setup homVarConst to be 0
        set_zero_d(RPD->homVarConst_d);
        set_zero_mp(RPD->homVarConst_mp);
        set_zero_rat(RPD->homVarConst_rat);
      }
      else
      { // using a homogeneous variable group

        // setup H to be random
        make_vec_random_rat(RPD->H_d, RPD->H_mp, RPD->H_rat, RPD->new_variables, RPD->curr_precision, T->AMP_max_prec, 0, 0);

        // setup homVarConst to be random
        get_comp_rand_rat(RPD->homVarConst_d, RPD->homVarConst_mp, RPD->homVarConst_rat, RPD->curr_precision, T->AMP_max_prec, 0, 0);
      }

      // setup patchCoeff
      init_vec_d(RPD->patchCoeff_d, RPD->new_variables);
      init_vec_mp2(RPD->patchCoeff_mp, RPD->new_variables, RPD->curr_precision);
      init_vec_rat(RPD->patchCoeff_rat, RPD->new_variables);
      make_vec_random_rat(RPD->patchCoeff_d, RPD->patchCoeff_mp, RPD->patchCoeff_rat, RPD->new_variables, RPD->curr_precision, T->AMP_max_prec, 0, 0);

      // setup A
      RPD->sameA = 1;
      RPD->A_d = (mat_d *)bmalloc(RPD->num_codim * sizeof(mat_d));
      RPD->A_mp = (mat_mp *)bmalloc(RPD->num_codim * sizeof(mat_mp));
      RPD->A_rat = (mpq_t ****)bmalloc(RPD->num_codim * sizeof(mpq_t ***));
      // make them in reverse order
      for (i = RPD->num_codim - 1; i >= 0; i--)
      { // setup A[i]
        init_mat_d(RPD->A_d[i], i + 1, RPD->num_funcs);
        init_mat_mp2(RPD->A_mp[i], i + 1, RPD->num_funcs, RPD->curr_precision);
        init_mat_rat(RPD->A_rat[i], i + 1, RPD->num_funcs);
        RPD->A_d[i]->rows = RPD->A_mp[i]->rows = i + 1;
        RPD->A_d[i]->cols = RPD->A_mp[i]->cols = RPD->num_funcs;

        if (i == RPD->num_codim - 1)
        { // we generate the main matrix
          make_matrix_random_rat(RPD->A_d[i], RPD->A_mp[i], RPD->A_rat[i], i + 1, RPD->num_funcs, RPD->curr_precision, T->AMP_max_prec, 0, 0);
          // make upper triangular with 1's on diagonal
          for (j = 0; j <= i; j++)
            for (k = 0; k <= j; k++)
              if (j == k)
              { // set to 1
                set_one_d(&RPD->A_d[i]->entry[j][k]);
                set_one_mp(&RPD->A_mp[i]->entry[j][k]);
                set_one_rat(RPD->A_rat[i][j][k]);
              }
              else
              { // set to 0
                set_zero_d(&RPD->A_d[i]->entry[j][k]);
                set_zero_mp(&RPD->A_mp[i]->entry[j][k]);
                set_zero_rat(RPD->A_rat[i][j][k]);
              }
        }
        else
        { // copy the top of A_d[i+1]
          for (j = 0; j <= i; j++)
            for (k = 0; k < RPD->num_funcs; k++)
            {
              set_d(&RPD->A_d[i]->entry[j][k], &RPD->A_d[i+1]->entry[j][k]);
              set_mp(&RPD->A_mp[i]->entry[j][k], &RPD->A_mp[i+1]->entry[j][k]);
              set_rat(RPD->A_rat[i][j][k], RPD->A_rat[i+1][j][k]);
            }
        }
      }
    }

    // setup coeff
    if (T->MPType == 0)
    { // setup coeff_d
      RPD->coeff_mp = NULL;
      RPD->coeff_rat = NULL;

      // allocate coeff_d make each one random
      RPD->coeff_d = (comp_d ***)bmalloc(RPD->num_codim * sizeof(comp_d **));
      for (i = 0; i < RPD->num_codim; i++)
      { // allocate for degree
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          for (k = 0; k < RPD->new_variables; k++)
          { // make random
            get_comp_rand_d(RPD->coeff_d[i][j][k]);
          }
        }
      }
    }
    else if (T->MPType == 1)
    { // setup coeff_mp
      RPD->coeff_d = NULL;
      RPD->coeff_rat = NULL;

      // allocate coeff_mp make each one random
      RPD->coeff_mp = (comp_mp ***)bmalloc(RPD->num_codim * sizeof(comp_mp **));
      for (i = 0; i < RPD->num_codim; i++)
      { // allocate for degree
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          for (k = 0; k < RPD->new_variables; k++)
          { // make random
            init_mp(RPD->coeff_mp[i][j][k]);
            get_comp_rand_mp(RPD->coeff_mp[i][j][k]);
          }
        }
      }
    }
    else
    { // setup coeff
      RPD->coeff_d = (comp_d ***)bmalloc(RPD->num_codim * sizeof(comp_d **));
      RPD->coeff_mp = (comp_mp ***)bmalloc(RPD->num_codim * sizeof(comp_mp **));
      RPD->coeff_rat = (mpq_t ****)bmalloc(RPD->num_codim * sizeof(mpq_t ***));
      for (i = 0; i < RPD->num_codim; i++)
      { // allocate for degree
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        RPD->coeff_rat[i] = (mpq_t ***)bmalloc(RPD->new_degrees[i] * sizeof(mpq_t **));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          RPD->coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
          for (k = 0; k < RPD->new_variables; k++)
          { // make random
            RPD->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
            init_rat(RPD->coeff_rat[i][j][k]);
            get_comp_rand_rat(RPD->coeff_d[i][j][k], RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k], RPD->curr_precision, T->AMP_max_prec, 0, 0);
          }
        }
      }
    }

    // setup the codimensions
    for (i = 0; i < RPD->num_codim; i++)
    { // setup the ith codimData
      setupRegenCodimData(RPD, i, i + 1, T->MPType, T->AMP_max_prec, intrinsicCutoff);
    }

    // setup the first level
    createFirstRegenCodim_startPoints(RPD, T->MPType, startFileName);
  }
  else if (startCodim <= RPD->num_codim + 1)
  { // we are starting with data from a previous run
    FILE *FP = fopen(startFileName, "r");
    // make sure FP exits
    if (FP == NULL)
    {
      printf("ERROR: The file to contain the start points, '%s', does not exist!\n", startFileName);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // setup from the file
    setupRPDRestart(RPD, T->MPType, T->AMP_max_prec, startCodim - 1, FP);

    // close FP
    fclose(FP);
  }
  else
  { // codim is too large
    printf("ERROR: The starting codimension (%d) is larger than the number of codimensions (%d)!\n", startCodim, RPD->num_codim);
    bexit(ERROR_CONFIGURATION);
  }

  // print message about codimensions
  if (specificCodim > 0)
  { // print a message
    printf("NOTE: You have requested to compute only codimension %d.\n", specificCodim);
  }
  else if (*maxCodim < RPD->num_codim)
  { // print a message 
    printf("NOTE: You have requested a maximum codimension of %d.\n", *maxCodim);
  }

  free(startFileName);

  return;
}

void createFirstRegenCodim_startPoints(regen_pos_dim_t *RPD, int MPType, char *startFileName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the start points for the first codim             *
\***************************************************************/
{
  FILE *START = fopen(startFileName, "w");
  int count;

  // the number of start points is the degree of the first function
  RPD->codim[0].num_paths = count = RPD->new_degrees[0];

  // print the number of start points to START
  fprintf(START, "%d\n\n", count);

  // create the start points
  writeRegenCodimStartPts(RPD, MPType, count, START);

  // print the relevant data to START so that it can be used to rerun this exact problem
  printRPDRelevantData(RPD, MPType, 0, START);

  // close START
  fclose(START);

  return;
}

void writeRegenCodimStartPts(regen_pos_dim_t *RPD, int MPType, int count, FILE *START)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: solves for the start points and prints them to START   *
\***************************************************************/
{
  int i, j, k, num_vars = RPD->new_variables, num_codim = RPD->num_codim;

  if (MPType == 0 || MPType == 2)
  { // find & print using double precision
    mat_d A, B_transpose;
    vec_d b, p, tempPt;

    init_mat_d(A, num_vars, num_vars); init_mat_d(B_transpose, 0, 0);
    init_vec_d(b, num_vars); init_vec_d(p, 0); init_vec_d(tempPt, 0);

    A->rows = A->cols = b->size = num_vars; // == num_codim + 1

    // see if we need to setup for intrinsic conversion
    if (RPD->codim[0].useIntrinsicSlice)
    { // setup B_transpose & p
      transpose_d(B_transpose, RPD->codim[0].B_d);
      vec_cp_d(p, RPD->codim[0].p_d);
    }

    // setup the patch
    if (RPD->PPD.num_var_gp)
    { // rhs is 1
      set_one_d(&b->coord[num_codim]);
      // coefficient as expected
      for (j = 0; j < num_vars; j++)
      {
        set_d(&A->entry[num_codim][j], &RPD->patchCoeff_d->coord[j]);
      }
    }
    else
    { // setup rhs
      set_d(&b->coord[num_codim], RPD->homVarConst_d);
      // coefficients are patchCoeff - H_d
      for (j = 0; j < num_vars; j++)
      {
        sub_d(&A->entry[num_codim][j], &RPD->patchCoeff_d->coord[j], &RPD->H_d->coord[j]);
      }
    }

    // copy over the linear coefficients for the bottom codim
    for (j = 1; j < num_codim; j++)
      for (k = 0; k < num_vars; k++)
        set_d(&A->entry[j][k], RPD->coeff_d[j][0][k]);

    // loop over the number of start points
    for (i = 0; i < count; i++)
    { // copy over the linears (i == degree of first linear)
      for (j = 0; j < num_codim; j++)
      { // b[j] == 0
        set_zero_d(&b->coord[j]);

        if (j == 0)
        { // setup the coeff for this point
          for (k = 0; k < num_vars; k++)
            set_d(&A->entry[j][k], RPD->coeff_d[j][i][k]);
        }
      }

      // solve for the start point
      if (matrixSolve_d(tempPt, A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      if (RPD->codim[0].useIntrinsicSlice)
      { // covert tempPt to intrinsic coordinates 
        extrinsicToIntrinsic_d(tempPt, tempPt, B_transpose, p);
      }

      // print the point to START
      for (j = 0; j < tempPt->size; j++)
        fprintf(START, "%.15e %.15e;\n", tempPt->coord[j].r, tempPt->coord[j].i);
      fprintf(START, "\n");
    }

    clear_mat_d(A); clear_mat_d(B_transpose);
    clear_vec_d(b); clear_vec_d(p); clear_vec_d(tempPt);
  }
  else
  { // find & print using fixed multi precision
    mat_mp A, B_transpose;
    vec_mp b, p, tempPt;

    init_mat_mp(A, num_vars, num_vars); init_mat_mp(B_transpose, 0, 0);
    init_vec_mp(b, num_vars); init_vec_mp(p, 0); init_vec_mp(tempPt, 0);

    A->rows = A->cols = b->size = num_vars; // == num_codim + 1

    // see if we need to setup for intrinsic conversion
    if (RPD->codim[0].useIntrinsicSlice)
    { // setup B_transpose & p
      transpose_mp(B_transpose, RPD->codim[0].B_mp);
      vec_cp_mp(p, RPD->codim[0].p_mp);
    }

    // setup the patch
    if (RPD->PPD.num_var_gp)
    { // rhs is 1
      set_one_mp(&b->coord[num_codim]);
      // coefficient as expected
      for (j = 0; j < num_vars; j++)
      {
        set_mp(&A->entry[num_codim][j], &RPD->patchCoeff_mp->coord[j]);
      }
    }
    else
    { // setup rhs
      set_mp(&b->coord[num_codim], RPD->homVarConst_mp);
      // coefficients are patchCoeff - H_d
      for (j = 0; j < num_vars; j++)
      {
        sub_mp(&A->entry[num_codim][j], &RPD->patchCoeff_mp->coord[j], &RPD->H_mp->coord[j]);
      }
    }

    // copy over the linear coefficients for the bottom codim
    for (j = 1; j < num_codim; j++)
      for (k = 0; k < num_vars; k++)
        set_mp(&A->entry[j][k], RPD->coeff_mp[j][0][k]);

    // loop over the number of start points
    for (i = 0; i < count; i++)
    { // copy over the linears (i == degree of first linear)
      for (j = 0; j < num_codim; j++)
      { // b[j] == 0
        set_zero_mp(&b->coord[j]);

        if (j == 0)
        { // setup the coeff for this point
          for (k = 0; k < num_vars; k++)
            set_mp(&A->entry[j][k], RPD->coeff_mp[j][i][k]);
        }
      }

      // solve for the start point
      if (matrixSolve_mp(tempPt, A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      if (RPD->codim[0].useIntrinsicSlice)
      { // covert tempPt to intrinsic coordinates
        extrinsicToIntrinsic_mp(tempPt, tempPt, B_transpose, p);
      }

      // print the point to START
      for (j = 0; j < tempPt->size; j++)
      {
        print_mp(START, 0, &tempPt->coord[j]);
        fprintf(START, ";\n");
      }
      fprintf(START, "\n");
    }

    clear_mat_mp(A); clear_mat_mp(B_transpose);
    clear_vec_mp(b); clear_vec_mp(p); clear_vec_mp(tempPt);
  }

  return;
}

void setupRegenCodimData(regen_pos_dim_t *RPD, int codim_index, int codim, int MPType, int max_prec, int intrinsicCutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the codim                                        *
\***************************************************************/
{
  RPD->codim[codim_index].codim = codim;
 
  // set counts to 0
  RPD->codim[codim_index].num_paths = RPD->codim[codim_index].num_superset = RPD->codim[codim_index].num_nonsing = RPD->codim[codim_index].num_sing =
    RPD->codim[codim_index].num_nonsolns = RPD->codim[codim_index].num_inf = RPD->codim[codim_index].num_bad = 0;

  // see if we are using an intrinsic slice
  if (codim > intrinsicCutoff)
  { // use extrinsic slice
    RPD->codim[codim_index].useIntrinsicSlice = 0;
  }
  else
  { // use intrinsic slice
    RPD->codim[codim_index].useIntrinsicSlice = 1;
  }

  if (RPD->codim[codim_index].useIntrinsicSlice)
  { // setup the intrinsic slice
    setupRegenCodimIntrinsicSlice(RPD, MPType, codim_index, max_prec);
  }
  else
  { // NULL out the pointers
    RPD->codim[codim_index].B_rat = NULL;
    RPD->codim[codim_index].p_rat = NULL;
  }

  return;
}

void setupRegenCodimIntrinsicSlice(regen_pos_dim_t *RPD, int MPType, int codim_index, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup B & p for using an intrinsic slice               *
\***************************************************************/
{
  int i, j, start, finish, count;

  // setup start & finish
  start = RPD->codim[codim_index].codim; // where the linears start
  finish = RPD->num_codim; // where the linears finish

  if (MPType == 0)
  { // setup B_d & p_d
    init_mat_d(RPD->codim[codim_index].B_d, 0, 0);
    init_vec_d(RPD->codim[codim_index].p_d, 0);

    vec_d b;
    mat_d tempMat, A, Q, R, P;
    double tol_pivot = 1e-15, tol_sign = 1e-20, largeChange = 1e14;

    init_vec_d(b, 0);
    init_mat_d(tempMat, 0, 0); init_mat_d(A, 0, 0); 
    init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

    // setup A
    count = finish - start + 1; // add the patch
    change_size_mat_d(A, count, RPD->new_variables);
    A->rows = count;
    A->cols = RPD->new_variables;

    // setup tempMat & b
    change_size_mat_d(tempMat, count, count);
    change_size_vec_d(b, count);
    tempMat->rows = tempMat->cols = b->size = count;

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < RPD->new_variables; j++)
      { // function i degree 0 variable j
        set_d(&A->entry[count][j], RPD->coeff_d[i][0][j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_d(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
      // set b[count] to 0
      set_zero_d(&b->coord[count]);
      count++;
    }
    // setup the patch
    if (RPD->PPD.num_var_gp)
    { // rhs is 1
      set_one_d(&b->coord[count]);
      // coefficient as expected
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_d(&A->entry[count][j], &RPD->patchCoeff_d->coord[j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_d(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }
    else
    { // setup rhs
      set_d(&b->coord[count], RPD->homVarConst_d);
      // coefficients are patchCoeff - H_d
      for (j = 0; j < RPD->new_variables; j++)
      {
        sub_d(&A->entry[count][j], &RPD->patchCoeff_d->coord[j], &RPD->H_d->coord[j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_d(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }

    // setup p - putting 0 in the extra positions
    // by doing it this way, we have a standard way of constructing p from coeff
    if (matrixSolve_d(RPD->codim[codim_index].p_d, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // correctly setup p
    count = tempMat->cols + start;
    increase_size_vec_d(RPD->codim[codim_index].p_d, count);
    RPD->codim[codim_index].p_d->size = count;
    for (j = tempMat->cols; j < count; j++)
    {
      set_zero_d(&RPD->codim[codim_index].p_d->coord[j]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_d(A, A);
    QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    change_size_mat_d(RPD->codim[codim_index].B_d, RPD->new_variables, start);
    RPD->codim[codim_index].B_d->rows = RPD->new_variables;
    RPD->codim[codim_index].B_d->cols = start;
    count = finish - start + 1;

    for (i = 0; i < RPD->new_variables; i++)
      for (j = 0; j < start; j++)
      {
        set_d(&RPD->codim[codim_index].B_d->entry[i][j], &Q->entry[i][j+count]);
      }

    clear_vec_d(b);
    clear_mat_d(tempMat); clear_mat_d(A); 
    clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
  }
  else if (MPType == 1)
  { // setup B_mp & p_mp
    init_mat_mp(RPD->codim[codim_index].B_mp, 0, 0);
    init_vec_mp(RPD->codim[codim_index].p_mp, 0);

    int num_digits = prec_to_digits(RPD->curr_precision);
    size_t size;
    char *str = NULL;

    vec_mp b;
    mat_mp tempMat, A, Q, R, P;
    mpf_t tol_pivot, tol_sign, largeChange;

    // initialize MP
    init_vec_mp(b, 0);
    init_mat_mp(tempMat, 0, 0); init_mat_mp(A, 0, 0);
    init_mat_mp(Q, 0, 0); init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0);
    mpf_init(tol_pivot); mpf_init(tol_sign); mpf_init(largeChange);

    // setup tol_pivot
    size = 1 + snprintf(NULL, 0, "1e-%d", num_digits-1);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", num_digits-1);
    mpf_set_str(tol_pivot, str, 10);
    // setup tol_sign
    size = 1 + snprintf(NULL, 0, "1e-%d", 2 * num_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", 2 * num_digits);
    mpf_set_str(tol_sign, str, 10);
    // setup largeChange
    size = 1 + snprintf(NULL, 0, "1e%d", num_digits-2);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", num_digits-2);
    mpf_set_str(largeChange, str, 10);

    // setup A
    count = finish - start + 1; // add the patch
    change_size_mat_mp(A, count, RPD->new_variables);
    A->rows = count; 
    A->cols = RPD->new_variables;

    // setup tempMat & b
    change_size_mat_mp(tempMat, count, count);
    change_size_vec_mp(b, count);
    tempMat->rows = tempMat->cols = b->size = count;

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < RPD->new_variables; j++)
      { // function i degree 0 variable j
        set_mp(&A->entry[count][j], RPD->coeff_mp[i][0][j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // setup the patch
    if (RPD->PPD.num_var_gp)
    { // rhs is 1
      set_one_mp(&b->coord[count]);
      // coefficient as expected
      for (j = 0; j < RPD->new_variables; j++)
      {
        set_mp(&A->entry[count][j], &RPD->patchCoeff_mp->coord[j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }
    else
    { // setup rhs
      set_mp(&b->coord[count], RPD->homVarConst_mp);
      // coefficients are patchCoeff - H_d
      for (j = 0; j < RPD->new_variables; j++)
      {
        sub_mp(&A->entry[count][j], &RPD->patchCoeff_mp->coord[j], &RPD->H_mp->coord[j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }

    // setup p - putting 0 in the extra positions
    // by doing it this way, we have a standard way of constructing p from coeff
    if (matrixSolve_mp(RPD->codim[codim_index].p_mp, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // correctly setup p
    count = tempMat->cols + start;
    increase_size_vec_mp(RPD->codim[codim_index].p_mp, count);
    RPD->codim[codim_index].p_mp->size = count;
    for (j = tempMat->cols; j < count; j++)
    {
      set_zero_mp(&RPD->codim[codim_index].p_mp->coord[j]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    change_size_mat_mp(RPD->codim[codim_index].B_mp, RPD->new_variables, start);
    RPD->codim[codim_index].B_mp->rows = RPD->new_variables;
    RPD->codim[codim_index].B_mp->cols = start;
    count = finish - start + 1;

    for (i = 0; i < RPD->new_variables; i++)
      for (j = 0; j < start; j++)
      {
        set_mp(&RPD->codim[codim_index].B_mp->entry[i][j], &Q->entry[i][j+count]);
      }

    clear_vec_mp(b);
    clear_mat_mp(tempMat); clear_mat_mp(A);
    clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);
    mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
    free(str);
  }
  else
  { // setup B & p
    init_mat_d(RPD->codim[codim_index].B_d, 0, 0);
    init_vec_d(RPD->codim[codim_index].p_d, 0);
    init_mat_mp(RPD->codim[codim_index].B_mp, 0, 0);
    init_vec_mp(RPD->codim[codim_index].p_mp, 0);

    int num_digits = prec_to_digits(max_prec);
    size_t size;
    char *str = NULL;

    comp_mp tempComp;
    vec_mp tempVec, b;
    mat_mp tempMat, A, Q, R, P;
    mpf_t tol_pivot, tol_sign, largeChange;

    // initialize MP
    init_mp2(tempComp, max_prec);
    init_vec_mp2(tempVec, 0, max_prec); init_vec_mp2(b, 0, max_prec);
    init_mat_mp2(tempMat, 0, 0, max_prec); init_mat_mp2(A, 0, 0, max_prec);
    init_mat_mp2(Q, 0, 0, max_prec); init_mat_mp2(R, 0, 0, max_prec); init_mat_mp2(P, 0, 0, max_prec);
    mpf_init2(tol_pivot, max_prec); mpf_init2(tol_sign, max_prec); mpf_init2(largeChange, max_prec);

    // setup tol_pivot
    size = 1 + snprintf(NULL, 0, "1e-%d", num_digits-1);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", num_digits-1);
    mpf_set_str(tol_pivot, str, 10);
    // setup tol_sign
    size = 1 + snprintf(NULL, 0, "1e-%d", 2 * num_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", 2 * num_digits);
    mpf_set_str(tol_sign, str, 10);
    // setup largeChange
    size = 1 + snprintf(NULL, 0, "1e%d", num_digits-2);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", num_digits-2);
    mpf_set_str(largeChange, str, 10);

    // setup A
    count = finish - start + 1; // add the patch
    change_size_mat_mp(A, count, RPD->new_variables);
    A->rows = count;
    A->cols = RPD->new_variables;

    // setup tempMat & b
    change_size_mat_mp(tempMat, count, count);
    change_size_vec_mp(b, count);
    tempMat->rows = tempMat->cols = b->size = count;

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < RPD->new_variables; j++)
      { // function i degree 0 variable j
        rat_to_mp(&A->entry[count][j], RPD->coeff_rat[i][0][j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }

    // set the global prec to the maximum precision
    initMP(max_prec);

    // setup the patch
    if (RPD->PPD.num_var_gp)
    { // rhs is 1
      set_one_mp(&b->coord[count]);
      // coefficient as expected
      for (j = 0; j < RPD->new_variables; j++)
      {
        rat_to_mp(&A->entry[count][j], RPD->patchCoeff_rat[j]);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }
    else
    { // setup rhs
      rat_to_mp(&b->coord[count], RPD->homVarConst_rat);
      // coefficients are patchCoeff - H_d
      for (j = 0; j < RPD->new_variables; j++)
      {
        rat_to_mp(&A->entry[count][j], RPD->patchCoeff_rat[j]);
        rat_to_mp(tempComp, RPD->H_rat[j]);
        sub_mp(&A->entry[count][j], &A->entry[count][j], tempComp);
        if (j < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
        }
      }
    }

    // solve for tempVec
    if (matrixSolve_mp(tempVec, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }
    // correctly setup tempVec
    count = tempMat->cols + start;
    increase_size_vec_mp(tempVec, count);
    tempVec->size = count;
    for (j = tempMat->cols; j < count; j++)
    {
      set_zero_mp(&tempVec->coord[j]);
    }

    // copy tempVec to p
    change_size_vec_d(RPD->codim[codim_index].p_d, tempVec->size);
    change_size_vec_mp(RPD->codim[codim_index].p_mp, tempVec->size);
    init_vec_rat(RPD->codim[codim_index].p_rat, tempVec->size);
    RPD->codim[codim_index].p_d->size = RPD->codim[codim_index].p_mp->size = tempVec->size;
    for (j = 0; j < tempVec->size; j++)
    { // copy tempVec to p_d, p_mp & p_rat
      mpf_t_to_rat(RPD->codim[codim_index].p_rat[j][0], tempVec->coord[j].r);
      mpf_t_to_rat(RPD->codim[codim_index].p_rat[j][1], tempVec->coord[j].i);
      rat_to_mp(&RPD->codim[codim_index].p_mp->coord[j], RPD->codim[codim_index].p_rat[j]);
      rat_to_d(&RPD->codim[codim_index].p_d->coord[j], RPD->codim[codim_index].p_rat[j]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);

    // allocate space and then setup B
    change_size_mat_d(RPD->codim[codim_index].B_d, RPD->new_variables, start);
    change_size_mat_mp(RPD->codim[codim_index].B_mp, RPD->new_variables, start);
    init_mat_rat(RPD->codim[codim_index].B_rat, RPD->new_variables, start);
    RPD->codim[codim_index].B_d->rows = RPD->codim[codim_index].B_mp->rows = RPD->new_variables;
    RPD->codim[codim_index].B_d->cols = RPD->codim[codim_index].B_mp->cols = start;
    count = finish - start + 1;

    for (i = 0; i < RPD->codim[codim_index].B_d->rows; i++)
      for (j = 0; j < RPD->codim[codim_index].B_d->cols; j++)
      { // copy Q to B_d, B_mp & B_rat
        mpf_t_to_rat(RPD->codim[codim_index].B_rat[i][j][0], Q->entry[i][j+count].r);
        mpf_t_to_rat(RPD->codim[codim_index].B_rat[i][j][1], Q->entry[i][j+count].i);
        rat_to_mp(&RPD->codim[codim_index].B_mp->entry[i][j], RPD->codim[codim_index].B_rat[i][j]);
        rat_to_d(&RPD->codim[codim_index].B_d->entry[i][j], RPD->codim[codim_index].B_rat[i][j]);
      }

    // set back to the current precision
    initMP(RPD->curr_precision);

    // clear MP
    mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
    clear_mp(tempComp);
    clear_vec_mp(tempVec); clear_vec_mp(b);
    clear_mat_mp(A); clear_mat_mp(tempMat);
    clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

    free(str);
  }

  return;
}

////// CHANGE PRECISION ///////

int change_regen_pos_dim_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for RPD                               *
\***************************************************************/
{
  int i, j, k;

  // cast ED as RPD
  regen_pos_dim_t *RPD = (regen_pos_dim_t *)ED;

  // set the SLP to the correct precision
  RPD->Prog->precision = prec;

  if (prec != RPD->curr_precision)
  { // need to change the precision
    RPD->curr_precision = prec;

    // change the precision for C, if needed
    if (RPD->new_variables != RPD->orig_variables)
    {
      for (i = 0; i < RPD->C_mp->rows; i++)
        for (j = 0; j < RPD->C_mp->cols; j++)
        {
          setprec_mp(&RPD->C_mp->entry[i][j], prec);
          rat_to_mp(&RPD->C_mp->entry[i][j], RPD->C_rat[i][j]);
        }
    }

    // change the precision for H
    for (i = 0; i < RPD->H_mp->size; i++)
    {
      setprec_mp(&RPD->H_mp->coord[i], prec);
      rat_to_mp(&RPD->H_mp->coord[i], RPD->H_rat[i]);
    }

    // change the precision for homVarConst
    setprec_mp(RPD->homVarConst_mp, prec);
    rat_to_mp(RPD->homVarConst_mp, RPD->homVarConst_rat);

    // change the precision for gamma
    setprec_mp(RPD->gamma_mp, prec);
    rat_to_mp(RPD->gamma_mp, RPD->gamma_rat);

    // change the precision for coeff
    for (i = 0; i < RPD->num_codim; i++)
      for (j = 0; j < RPD->new_degrees[i]; j++)
        for (k = 0; k < RPD->new_variables; k++)
        {
          setprec_mp(RPD->coeff_mp[i][j][k], prec);
          rat_to_mp(RPD->coeff_mp[i][j][k], RPD->coeff_rat[i][j][k]);
        }

    // change the precision for patchCoeff
    for (i = 0; i < RPD->patchCoeff_mp->size; i++)
    {
      setprec_mp(&RPD->patchCoeff_mp->coord[i], prec);
      rat_to_mp(&RPD->patchCoeff_mp->coord[i], RPD->patchCoeff_rat[i]);
    }

    for (k = 0; k < RPD->num_codim; k++)
    { // change the precision for A
      for (i = 0; i < RPD->A_mp[k]->rows; i++)
        for (j = 0; j < RPD->A_mp[k]->cols; j++)
        {
          setprec_mp(&RPD->A_mp[k]->entry[i][j], prec);
          rat_to_mp(&RPD->A_mp[k]->entry[i][j], RPD->A_rat[k][i][j]);
        }

      if (RPD->codim[k].useIntrinsicSlice)
      { // change the precision for B
        for (i = 0; i < RPD->codim[k].B_mp->rows; i++)
          for (j = 0; j < RPD->codim[k].B_mp->cols; j++)
          {
            setprec_mp(&RPD->codim[k].B_mp->entry[i][j], prec);
            rat_to_mp(&RPD->codim[k].B_mp->entry[i][j], RPD->codim[k].B_rat[i][j]);
          }

        // change the precision for p
        for (i = 0; i < RPD->codim[k].p_mp->size; i++)
        {
          setprec_mp(&RPD->codim[k].p_mp->coord[i], prec);
          rat_to_mp(&RPD->codim[k].p_mp->coord[i], RPD->codim[k].p_rat[i]);
        }
      }
    }
  }

  return 0;
}

void printRPDSummaryData(regen_pos_dim_t *RPD, int codim_index, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the summary data to FP                          *
\***************************************************************/
{
  int i;

  // print the num_codim, num_funcs, orig_variables, new_variables, system_rank
  fprintf(FP, "%d %d %d %d %d\n", RPD->num_codim, RPD->num_funcs, RPD->orig_variables, RPD->new_variables, RPD->system_rank);

  // print the current codimension
  fprintf(FP, "%d\n", RPD->codim[codim_index].codim);

  // print the info about each codim
  for (i = 0; i < RPD->num_codim; i++)
    if (i <= codim_index)
    { // print info about how the codim was
      fprintf(FP, "%d %d %d %d %d %d %d %d\n", RPD->codim[i].codim, RPD->codim[i].num_paths, RPD->codim[i].num_superset, RPD->codim[i].num_nonsing, RPD->codim[i].num_sing, RPD->codim[i].num_nonsolns, RPD->codim[i].num_inf, RPD->codim[i].num_bad);
    }
    else
    { // print info about future codims
      fprintf(FP, "%d %d\n", RPD->codim[i].codim, RPD->codim[i].useIntrinsicSlice);
    }
  fprintf(FP, "\n");

  return;
}

///// PRINT RELEVANT DATA TO FILE ///////

void printRPDRelevantData(regen_pos_dim_t *RPD, int MPType, int codim_index, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the relevant data to FP so that we can begin    *
* exactly right here in case of failure or change tolerances    *
\***************************************************************/
{
  int i, j, k, maxCodim = RPD->system_rank;

  // print an X signifying that this is extra info
  fprintf(FP, "X\n");

  // print the MPType, num_codim, num_funcs, orig_variables, new_variables, system_rank
  fprintf(FP, "%d %d %d %d %d %d\n", MPType, RPD->num_codim, RPD->num_funcs, RPD->orig_variables, RPD->new_variables, RPD->system_rank);

  // print the current codim_index
  fprintf(FP, "%d\n", codim_index);

  // print the info about each codim
  for (i = 0; i < RPD->num_codim; i++)
    if (i < codim_index)
    { // print info about how the codim was
      fprintf(FP, "%d %d %d %d %d %d %d %d\n", RPD->codim[i].codim, RPD->codim[i].num_paths, RPD->codim[i].num_superset, RPD->codim[i].num_nonsing, RPD->codim[i].num_sing, RPD->codim[i].num_nonsolns, RPD->codim[i].num_inf, RPD->codim[i].num_bad);
    }
    else if (i == codim_index)
    { // print info about the current codim
      fprintf(FP, "%d %d %d\n", RPD->codim[i].codim, RPD->codim[i].num_paths, RPD->codim[i].useIntrinsicSlice);
    }
    else
    { // print info about future codims
      fprintf(FP, "%d %d\n", RPD->codim[i].codim, RPD->codim[i].useIntrinsicSlice);
    }
  fprintf(FP, "\n");

  // print C, H, homVarConst, patchCoeff
  if (MPType == 0)
  { // print _d
    print_mat_out_d(FP, RPD->C_d);
    print_vec_out_d(FP, RPD->H_d);
    print_comp_out_d(FP, RPD->homVarConst_d);
    print_vec_out_d(FP, RPD->patchCoeff_d);
  }
  else if (MPType == 1)
  { // print _mp
    print_mat_out_mp(FP, RPD->C_mp);
    print_vec_out_mp(FP, RPD->H_mp);
    print_comp_out_mp(FP, RPD->homVarConst_mp);
    print_vec_out_mp(FP, RPD->patchCoeff_mp);
  }
  else
  { // print _rat
    print_mat_out_rat(FP, RPD->C_rat, RPD->C_d->rows, RPD->C_d->cols);
    print_vec_out_rat(FP, RPD->H_rat, RPD->H_d->size);
    print_comp_out_rat(FP, RPD->homVarConst_rat);
    print_vec_out_rat(FP, RPD->patchCoeff_rat, RPD->patchCoeff_d->size);
  }
  
  // print coeff
  if (MPType == 0)
  { // print coeff_d
    for (i = 0; i < maxCodim; i++)
      for (j = 0; j < RPD->new_degrees[i]; j++)
        for (k = 0; k < RPD->new_variables; k++)
          print_comp_out_d(FP, RPD->coeff_d[i][j][k]);
  }
  else if (MPType == 1)
  { // print coeff_mp
    for (i = 0; i < maxCodim; i++)
      for (j = 0; j < RPD->new_degrees[i]; j++)
        for (k = 0; k < RPD->new_variables; k++)
          print_comp_out_mp(FP, RPD->coeff_mp[i][j][k]);
  }
  else
  { // print coeff_rat
    for (i = 0; i < maxCodim; i++)
      for (j = 0; j < RPD->new_degrees[i]; j++)
        for (k = 0; k < RPD->new_variables; k++)
          print_comp_out_rat(FP, RPD->coeff_rat[i][j][k]);
  }

  // print A
  fprintf(FP, "%d\n", RPD->sameA);
  if (MPType == 0)
  { // print A_d
    for (i = 0; i < maxCodim; i++)
      print_mat_out_d(FP, RPD->A_d[i]);
  }
  else if (MPType == 1)
  { // print A_mp
    for (i = 0; i < maxCodim; i++)
      print_mat_out_mp(FP, RPD->A_mp[i]);
  }
  else
  { // print A_rat
    for (i = 0; i < maxCodim; i++)
      print_mat_out_rat(FP, RPD->A_rat[i], RPD->A_d[i]->rows, RPD->A_d[i]->cols);
  }
  fprintf(FP, "\n");

  return;
}

///// READ IN RELEVANT DATA FROM FILE ///////

void setupRPDRestart(regen_pos_dim_t *RPD, int MPType, int max_prec, int codim_index, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup RPD from the data in FP                          *
\***************************************************************/
{
  int i, j, k, file_codim_index, file_MPType, num_codim, num_funcs, orig_variables, new_variables, system_rank;
  char ch;
  int H_size = 0, patchCoeff_size = 0, C_rows = 0, C_cols = 0, sameA = 0, *A_rows = NULL, *A_cols = NULL;
  comp_d homVarConst_d, ***coeff_d = NULL;
  comp_mp homVarConst_mp, ***coeff_mp = NULL;
  mpq_t homVarConst_rat[2], ****coeff_rat = NULL;
  vec_d H_d, patchCoeff_d;
  vec_mp H_mp, patchCoeff_mp;
  mpq_t **H_rat = NULL, **patchCoeff_rat = NULL;
  mat_d C_d, *A_d = NULL;
  mat_mp C_mp, *A_mp = NULL;
  mpq_t ***C_rat = NULL, ****A_rat = NULL;

  // move past all the other data and find the 'X'
  do
  {
    ch = fgetc(FP);
  }
  while (ch != 'X');

  // read in MPType, num_codim, num_funcs, orig_variables, new_variables, system_rank
  fscanf(FP, "%d%d%d%d%d%d\n", &file_MPType, &num_codim, &num_funcs, &orig_variables, &new_variables, &system_rank);
  // read in the current codimension index
  fscanf(FP, "%d\n", &file_codim_index);

  // make sure we have agreement on num_funcs, orig_variables, new_variables & system_rank
  if (num_codim != RPD->num_codim)
  {
    printf("ERROR: The number of codimensions (%d vs %d) is not correct!\n", RPD->num_codim, num_codim);
    bexit(ERROR_CONFIGURATION);
  }
  if (num_funcs != RPD->num_funcs)
  {
    printf("ERROR: The number of functions (%d vs %d) is not correct!\n", RPD->num_funcs, num_funcs);
    bexit(ERROR_CONFIGURATION);
  }
  if (orig_variables != RPD->orig_variables)
  {
    printf("ERROR: The number of original variables (%d vs %d) is not correct!\n", RPD->orig_variables, orig_variables);
    bexit(ERROR_CONFIGURATION);
  }
  if (new_variables != RPD->new_variables)
  {
    printf("ERROR: The number of tracking variables (%d vs %d) is not correct!\n", RPD->new_variables, new_variables);
    bexit(ERROR_CONFIGURATION);
  }
  if (system_rank != RPD->system_rank)
  {
    printf("ERROR: The rank of the system (%d vs %d) is not correct!\n", RPD->system_rank, system_rank);
    bexit(ERROR_CONFIGURATION);
  }
  if (file_codim_index != codim_index)
  {
    printf("ERROR: The current codimension index (%d vs %d) is not correct!\n", file_codim_index, codim_index);
    bexit(ERROR_CONFIGURATION);
  }

  // so we can assume that we probably have the correct file

  // NOTE: codim has been allocated, but not yet setup

  // read in info about each codim
  for (i = 0; i < num_codim; i++)
    if (i < codim_index)
    { // read in info about this codim
      fscanf(FP, "%d%d%d%d%d%d%d%d\n", &RPD->codim[i].codim, &RPD->codim[i].num_paths, &RPD->codim[i].num_superset, &RPD->codim[i].num_nonsing, &RPD->codim[i].num_sing, &RPD->codim[i].num_nonsolns, &RPD->codim[i].num_inf, &RPD->codim[i].num_bad);
    }
    else if (i == codim_index)
    { // read in info about the current codim
      fscanf(FP, "%d%d%d\n", &RPD->codim[i].codim, &RPD->codim[i].num_paths, &RPD->codim[i].useIntrinsicSlice);
    }
    else
    { // read in info about future codims
      fscanf(FP, "%d%d\n", &RPD->codim[i].codim, &RPD->codim[i].useIntrinsicSlice);
    }

  // read in C, H, homVarConst, patchCoeff
  if (file_MPType == 0)
  { // read in _d
    init_mat_d(C_d, 0, 0);
    setup_mat_in_d(C_d, FP);
    init_vec_d(H_d, 0);
    setup_vec_in_d(H_d, FP);
    init_d(homVarConst_d);
    setup_comp_in_d(homVarConst_d, FP);
    init_vec_d(patchCoeff_d, 0);
    setup_vec_in_d(patchCoeff_d, FP);    
  }
  else if (file_MPType == 1)
  { // read in _mp
    init_mat_mp(C_mp, 0, 0);
    setup_mat_in_mp(C_mp, FP);
    init_vec_mp(H_mp, 0);
    setup_vec_in_mp(H_mp, FP);
    init_mp(homVarConst_mp);
    setup_comp_in_mp(homVarConst_mp, FP);
    init_vec_mp(patchCoeff_mp, 0);
    setup_vec_in_mp(patchCoeff_mp, FP);
  }
  else
  { // read in _rat
    setup_mat_in_rat(&C_rat, FP, &C_rows, &C_cols);
    setup_vec_in_rat(&H_rat, FP, &H_size);
    setup_comp_in_rat(homVarConst_rat, FP);
    setup_vec_in_rat(&patchCoeff_rat, FP, &patchCoeff_size);
  }

  // read in coeff
  if (file_MPType == 0)
  { // read in coeff_d
    coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
    for (i = 0; i < system_rank; i++)
    {
      coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      {
        coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
        for (k = 0; k < RPD->new_variables; k++)
        {
          init_d(coeff_d[i][j][k]);
          setup_comp_in_d(coeff_d[i][j][k], FP);
        }
      }
    }
  }
  else if (file_MPType == 1)
  { // read in coeff_mp
    coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
    for (i = 0; i < system_rank; i++)
    {
      coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      {
        coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
        for (k = 0; k < RPD->new_variables; k++)
        {
          init_mp(coeff_mp[i][j][k]);
          setup_comp_in_mp(coeff_mp[i][j][k], FP);
        }
      }
    }
  }
  else
  { // read in coeff_rat
    coeff_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
    for (i = 0; i < system_rank; i++)
    {
      coeff_rat[i] = (mpq_t ***)bmalloc(RPD->new_degrees[i] * sizeof(mpq_t **));
      for (j = 0; j < RPD->new_degrees[i]; j++)
      {
        coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
        for (k = 0; k < RPD->new_variables; k++)
        {
          coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          setup_comp_in_rat(coeff_rat[i][j][k], FP);
        }
      }
    }
  }

  // read in A
  fscanf(FP, "%d\n", &sameA);
  if (file_MPType == 0)
  { // read in A_d
    A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
    for (i = 0; i < system_rank; i++)
    {
      init_mat_d(A_d[i], 0, 0);
      setup_mat_in_d(A_d[i], FP);
    }
  }
  else if (file_MPType == 1)
  { // read in A_mp
    A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
    for (i = 0; i < system_rank; i++)
    {
      init_mat_mp(A_mp[i], 0, 0);
      setup_mat_in_mp(A_mp[i], FP);
    }
  }
  else
  { // read in A_rat
    A_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
    A_rows = (int *)bmalloc(system_rank * sizeof(int));
    A_cols = (int *)bmalloc(system_rank * sizeof(int));
    for (i = 0; i < system_rank; i++)
      setup_mat_in_rat(&A_rat[i], FP, &A_rows[i], &A_cols[i]);
  }

  // setup RPD and clear structures, namely C, H, homVarConst, patchCoeff, coeff & A
  if (MPType == 0)
  { // setup _d 
    if (file_MPType == 0)
    { // copy _d to _d
      init_mat_d(RPD->C_d, 0, 0);
      mat_cp_d(RPD->C_d, C_d);
      clear_mat_d(C_d);

      init_vec_d(RPD->H_d, 0);
      vec_cp_d(RPD->H_d, H_d);
      clear_vec_d(H_d);

      init_d(RPD->homVarConst_d);
      set_d(RPD->homVarConst_d, homVarConst_d);
      clear_d(homVarConst_d);

      init_vec_d(RPD->patchCoeff_d, 0);
      vec_cp_d(RPD->patchCoeff_d, patchCoeff_d);
      clear_vec_d(patchCoeff_d);

      RPD->coeff_d = coeff_d;
      coeff_d = NULL;
      RPD->coeff_mp = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_d = A_d;
      A_d = NULL;
      RPD->A_mp = NULL;
      RPD->A_rat = NULL;
    }
    else if (file_MPType == 1)
    { // convert _mp to _d
      init_mat_d(RPD->C_d, 0, 0);
      mat_mp_to_d(RPD->C_d, C_mp);
      clear_mat_mp(C_mp);

      init_vec_d(RPD->H_d, 0);
      vec_mp_to_d(RPD->H_d, H_mp);
      clear_vec_mp(H_mp);

      init_d(RPD->homVarConst_d);
      mp_to_d(RPD->homVarConst_d, homVarConst_mp);
      clear_mp(homVarConst_mp);

      init_vec_d(RPD->patchCoeff_d, 0);
      vec_mp_to_d(RPD->patchCoeff_d, patchCoeff_mp);
      clear_vec_mp(patchCoeff_mp);

      RPD->coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
      for (i = 0; i < system_rank; i++)
      {
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        {
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          for (k = 0; k < RPD->new_variables; k++) 
          {
            init_d(RPD->coeff_d[i][j][k]);
            mp_to_d(RPD->coeff_d[i][j][k], coeff_mp[i][j][k]);
            clear_mp(coeff_mp[i][j][k]);
          }
          free(coeff_mp[i][j]);
        }
        free(coeff_mp[i]);
      }
      free(coeff_mp);
      RPD->coeff_mp = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_d(RPD->A_d[i], 0, 0);
        mat_mp_to_d(RPD->A_d[i], A_mp[i]);
        clear_mat_mp(A_mp[i]);
      }
      free(A_mp);
      RPD->A_mp = NULL;
      RPD->A_rat = NULL;
    }
    else 
    { // convert _rat to _d
      init_mat_d(RPD->C_d, 0, 0);
      mat_rat_to_d(RPD->C_d, C_rat, C_rows, C_cols);
      clear_mat_rat(C_rat, C_rows, C_cols);

      init_vec_d(RPD->H_d, 0);
      vec_rat_to_d(RPD->H_d, H_rat, H_size);
      clear_vec_rat(H_rat, H_size);

      init_d(RPD->homVarConst_d);
      rat_to_d(RPD->homVarConst_d, homVarConst_rat);
      clear_rat(homVarConst_rat);

      init_vec_d(RPD->patchCoeff_d, 0);
      vec_rat_to_d(RPD->patchCoeff_d, patchCoeff_rat, patchCoeff_size);
      clear_vec_rat(patchCoeff_rat, patchCoeff_size);

      RPD->coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
      for (i = 0; i < system_rank; i++)
      {
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        {
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          for (k = 0; k < RPD->new_variables; k++)
          {
            init_d(RPD->coeff_d[i][j][k]);
            rat_to_d(RPD->coeff_d[i][j][k], coeff_rat[i][j][k]);
            clear_rat(coeff_rat[i][j][k]);
            free(coeff_rat[i][j][k]);
          }
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
      RPD->coeff_mp = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_d(RPD->A_d[i], 0, 0);
        mat_rat_to_d(RPD->A_d[i], A_rat[i], A_rows[i], A_cols[i]);
        clear_mat_rat(A_rat[i], A_rows[i], A_cols[i]);
      }
      free(A_rat);
      free(A_rows);
      free(A_cols);
      RPD->A_mp = NULL;
      RPD->A_rat = NULL;
    }
  }
  else if (MPType == 1)
  { // setup _mp
    if (file_MPType == 0)
    { // copy _d to _mp
      init_mat_mp(RPD->C_mp, 0, 0);
      mat_d_to_mp(RPD->C_mp, C_d);
      clear_mat_d(C_d);

      init_vec_mp(RPD->H_mp, 0);
      vec_d_to_mp(RPD->H_mp, H_d);
      clear_vec_d(H_d);

      init_mp(RPD->homVarConst_mp);
      d_to_mp(RPD->homVarConst_mp, homVarConst_d);
      clear_d(homVarConst_d);

      init_vec_mp(RPD->patchCoeff_mp, 0);
      vec_d_to_mp(RPD->patchCoeff_mp, patchCoeff_d);
      clear_vec_d(patchCoeff_d);

      RPD->coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
      for (i = 0; i < system_rank; i++)
      {
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        {
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          for (k = 0; k < RPD->new_variables; k++)
          {
            init_mp(RPD->coeff_mp[i][j][k]);
            d_to_mp(RPD->coeff_mp[i][j][k], coeff_d[i][j][k]);
            clear_d(coeff_d[i][j][k]);
          }
          free(coeff_d[i][j]);
        }
        free(coeff_d[i]);
      }
      free(coeff_d);
      RPD->coeff_d = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_mp(RPD->A_mp[i], 0, 0);
        mat_d_to_mp(RPD->A_mp[i], A_d[i]);
        clear_mat_d(A_d[i]);
      }
      free(A_d);
      RPD->A_d = NULL;
      RPD->A_rat = NULL;
    }
    else if (file_MPType == 1)
    { // copy _mp to _mp
      init_mat_mp(RPD->C_mp, 0, 0);
      mat_cp_mp(RPD->C_mp, C_mp);
      clear_mat_mp(C_mp);

      init_vec_mp(RPD->H_mp, 0);
      vec_cp_mp(RPD->H_mp, H_mp);
      clear_vec_mp(H_mp);

      init_mp(RPD->homVarConst_mp);
      set_mp(RPD->homVarConst_mp, homVarConst_mp);
      clear_mp(homVarConst_mp);

      init_vec_mp(RPD->patchCoeff_mp, 0);
      vec_cp_mp(RPD->patchCoeff_mp, patchCoeff_mp);
      clear_vec_mp(patchCoeff_mp);

      RPD->coeff_mp = coeff_mp;
      coeff_mp = NULL;
      RPD->coeff_d = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_mp = A_mp;
      A_mp = NULL;
      RPD->A_d = NULL;
      RPD->A_rat = NULL;
    }
    else
    { // convert _rat to _mp
      init_mat_mp(RPD->C_mp, 0, 0);
      mat_rat_to_mp(RPD->C_mp, C_rat, C_rows, C_cols);
      clear_mat_rat(C_rat, C_rows, C_cols);

      init_vec_mp(RPD->H_mp, 0);
      vec_rat_to_mp(RPD->H_mp, H_rat, H_size);
      clear_vec_rat(H_rat, H_size);

      init_mp(RPD->homVarConst_mp);
      rat_to_mp(RPD->homVarConst_mp, homVarConst_rat);
      clear_rat(homVarConst_rat);

      init_vec_mp(RPD->patchCoeff_mp, 0);
      vec_rat_to_mp(RPD->patchCoeff_mp, patchCoeff_rat, patchCoeff_size);
      clear_vec_rat(patchCoeff_rat, patchCoeff_size);

      RPD->coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
      for (i = 0; i < system_rank; i++)
      {
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        {
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          for (k = 0; k < RPD->new_variables; k++)
          {
            init_mp(RPD->coeff_mp[i][j][k]);
            rat_to_mp(RPD->coeff_mp[i][j][k], coeff_rat[i][j][k]);
            clear_rat(coeff_rat[i][j][k]);
            free(coeff_rat[i][j][k]);
          }
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
      RPD->coeff_d = NULL;
      RPD->coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_mp(RPD->A_mp[i], 0, 0);
        mat_rat_to_mp(RPD->A_mp[i], A_rat[i], A_rows[i], A_cols[i]);
        clear_mat_rat(A_rat[i], A_rows[i], A_cols[i]);
      }
      free(A_rat);
      free(A_rows);
      free(A_cols);
      RPD->A_d = NULL;
      RPD->A_rat = NULL;
    }
  }
  else
  { // setup _d, _mp & _rat
    if (file_MPType == 0)
    { // copy _d to _d, _mp & _rat
      init_mat_d(RPD->C_d, C_d->rows, C_d->cols);
      init_mat_mp2(RPD->C_mp, C_d->rows, C_d->cols, RPD->curr_precision);
      init_mat_rat(RPD->C_rat, C_d->rows, C_d->cols);
      mat_cp_d(RPD->C_d, C_d);
      mat_d_to_mp(RPD->C_mp, C_d);
      mat_d_to_rat(RPD->C_rat, C_d);
      clear_mat_d(C_d);

      init_vec_d(RPD->H_d, H_d->size);
      init_vec_mp2(RPD->H_mp, H_d->size, RPD->curr_precision);
      init_vec_rat(RPD->H_rat, H_d->size);
      vec_cp_d(RPD->H_d, H_d);
      vec_d_to_mp(RPD->H_mp, H_d);
      vec_d_to_rat(RPD->H_rat, H_d);
      clear_vec_d(H_d);

      init_d(RPD->homVarConst_d);
      init_mp2(RPD->homVarConst_mp, RPD->curr_precision);
      init_rat(RPD->homVarConst_rat);
      set_d(RPD->homVarConst_d, homVarConst_d);
      d_to_mp(RPD->homVarConst_mp, homVarConst_d);
      d_to_rat(RPD->homVarConst_rat, homVarConst_d);
      clear_d(homVarConst_d);

      init_vec_d(RPD->patchCoeff_d, patchCoeff_d->size);
      init_vec_mp2(RPD->patchCoeff_mp, patchCoeff_d->size, RPD->curr_precision);
      init_vec_rat(RPD->patchCoeff_rat, patchCoeff_d->size);
      vec_cp_d(RPD->patchCoeff_d, patchCoeff_d);
      vec_d_to_mp(RPD->patchCoeff_mp, patchCoeff_d);
      vec_d_to_rat(RPD->patchCoeff_rat, patchCoeff_d);
      clear_vec_d(patchCoeff_d);

      RPD->coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
      RPD->coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
      RPD->coeff_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
      for (i = 0; i < system_rank; i++)
      { // allocate for degree
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        RPD->coeff_rat[i] = (mpq_t ***)bmalloc(RPD->new_degrees[i] * sizeof(mpq_t **));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          RPD->coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
          for (k = 0; k < RPD->new_variables; k++)
          {
            RPD->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            init_d(RPD->coeff_d[i][j][k]);
            init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
            init_rat(RPD->coeff_rat[i][j][k]);
            set_d(RPD->coeff_d[i][j][k], coeff_d[i][j][k]);
            d_to_mp(RPD->coeff_mp[i][j][k], coeff_d[i][j][k]);
            d_to_rat(RPD->coeff_rat[i][j][k], coeff_d[i][j][k]);
            clear_d(coeff_d[i][j][k]);
          }
          free(coeff_d[i][j]);
        }
        free(coeff_d[i]);
      }
      free(coeff_d);
    
      RPD->sameA = sameA;
      RPD->A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
      RPD->A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
      RPD->A_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_d(RPD->A_d[i], A_d[i]->rows, A_d[i]->cols);
        init_mat_mp2(RPD->A_mp[i], A_d[i]->rows, A_d[i]->cols, RPD->curr_precision);
        init_mat_rat(RPD->A_rat[i], A_d[i]->rows, A_d[i]->cols);
        mat_cp_d(RPD->A_d[i], A_d[i]);
        mat_d_to_mp(RPD->A_mp[i], A_d[i]);
        mat_d_to_rat(RPD->A_rat[i], A_d[i]);
        clear_mat_d(A_d[i]);
      }
      free(A_d);
    }
    else if (file_MPType == 1)
    { // copy _mp to _d, _mp & _rat
      init_mat_d(RPD->C_d, C_mp->rows, C_mp->cols);
      init_mat_mp2(RPD->C_mp, C_mp->rows, C_mp->cols, RPD->curr_precision);
      init_mat_rat(RPD->C_rat, C_mp->rows, C_mp->cols);
      mat_mp_to_d(RPD->C_d, C_mp);
      mat_cp_mp(RPD->C_mp, C_mp);
      mat_mp_to_rat(RPD->C_rat, C_mp);
      clear_mat_mp(C_mp);

      init_vec_d(RPD->H_d, H_mp->size);
      init_vec_mp2(RPD->H_mp, H_mp->size, RPD->curr_precision);
      init_vec_rat(RPD->H_rat, H_mp->size);
      vec_mp_to_d(RPD->H_d, H_mp);
      vec_cp_mp(RPD->H_mp, H_mp);
      vec_mp_to_rat(RPD->H_rat, H_mp);
      clear_vec_mp(H_mp);

      init_d(RPD->homVarConst_d);
      init_mp2(RPD->homVarConst_mp, RPD->curr_precision);
      init_rat(RPD->homVarConst_rat);
      mp_to_d(RPD->homVarConst_d, homVarConst_mp);
      set_mp(RPD->homVarConst_mp, homVarConst_mp);
      mp_to_rat(RPD->homVarConst_rat, homVarConst_mp);
      clear_mp(homVarConst_mp);

      init_vec_d(RPD->patchCoeff_d, patchCoeff_mp->size);
      init_vec_mp2(RPD->patchCoeff_mp, patchCoeff_mp->size, RPD->curr_precision);
      init_vec_rat(RPD->patchCoeff_rat, patchCoeff_mp->size);
      vec_mp_to_d(RPD->patchCoeff_d, patchCoeff_mp);
      vec_cp_mp(RPD->patchCoeff_mp, patchCoeff_mp);
      vec_mp_to_rat(RPD->patchCoeff_rat, patchCoeff_mp);
      clear_vec_mp(patchCoeff_mp);

      RPD->coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
      RPD->coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
      RPD->coeff_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
      for (i = 0; i < system_rank; i++)
      { // allocate for degree
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        RPD->coeff_rat[i] = (mpq_t ***)bmalloc(RPD->new_degrees[i] * sizeof(mpq_t **));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          RPD->coeff_rat[i][j] = (mpq_t **)bmalloc(RPD->new_variables * sizeof(mpq_t *));
          for (k = 0; k < RPD->new_variables; k++)
          {
            RPD->coeff_rat[i][j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            init_d(RPD->coeff_d[i][j][k]);
            init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
            init_rat(RPD->coeff_rat[i][j][k]);
            mp_to_d(RPD->coeff_d[i][j][k], coeff_mp[i][j][k]);
            set_mp(RPD->coeff_mp[i][j][k], coeff_mp[i][j][k]);
            mp_to_rat(RPD->coeff_rat[i][j][k], coeff_mp[i][j][k]);
            clear_mp(coeff_mp[i][j][k]);
          }
          free(coeff_mp[i][j]);
        }
        free(coeff_mp[i]);
      }
      free(coeff_mp);

      RPD->sameA = sameA;
      RPD->A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
      RPD->A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
      RPD->A_rat = (mpq_t ****)bmalloc(system_rank * sizeof(mpq_t ***));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_d(RPD->A_d[i], A_mp[i]->rows, A_mp[i]->cols);
        init_mat_mp2(RPD->A_mp[i], A_mp[i]->rows, A_mp[i]->cols, RPD->curr_precision);
        init_mat_rat(RPD->A_rat[i], A_mp[i]->rows, A_mp[i]->cols);
        mat_mp_to_d(RPD->A_d[i], A_mp[i]);
        mat_cp_mp(RPD->A_mp[i], A_mp[i]);
        mat_mp_to_rat(RPD->A_rat[i], A_mp[i]);
        clear_mat_mp(A_mp[i]);
      }
      free(A_mp);
    }
    else 
    { // copy _rat to _d, _mp & _rat
      init_mat_d(RPD->C_d, C_rows, C_cols);
      init_mat_mp2(RPD->C_mp, C_rows, C_cols, RPD->curr_precision);
      mat_rat_to_d(RPD->C_d, C_rat, C_rows, C_cols);
      mat_rat_to_mp(RPD->C_mp, C_rat, C_rows, C_cols);
      RPD->C_rat = C_rat;
      C_rat = NULL;

      init_vec_d(RPD->H_d, H_size);
      init_vec_mp2(RPD->H_mp, H_size, RPD->curr_precision);
      vec_rat_to_d(RPD->H_d, H_rat, H_size);
      vec_rat_to_mp(RPD->H_mp, H_rat, H_size);
      RPD->H_rat = H_rat;
      H_rat = NULL;

      init_d(RPD->homVarConst_d);
      init_mp2(RPD->homVarConst_mp, RPD->curr_precision);
      init_rat(RPD->homVarConst_rat);
      rat_to_d(RPD->homVarConst_d, homVarConst_rat);
      rat_to_mp(RPD->homVarConst_mp, homVarConst_rat);
      set_rat(RPD->homVarConst_rat, homVarConst_rat);
      clear_rat(homVarConst_rat);

      init_vec_d(RPD->patchCoeff_d, patchCoeff_size);
      init_vec_mp2(RPD->patchCoeff_mp, patchCoeff_size, RPD->curr_precision);
      vec_rat_to_d(RPD->patchCoeff_d, patchCoeff_rat, patchCoeff_size);
      vec_rat_to_mp(RPD->patchCoeff_mp, patchCoeff_rat, patchCoeff_size);
      RPD->patchCoeff_rat = patchCoeff_rat;
      patchCoeff_rat = NULL;

      RPD->coeff_d = (comp_d ***)bmalloc(system_rank * sizeof(comp_d **));
      RPD->coeff_mp = (comp_mp ***)bmalloc(system_rank * sizeof(comp_mp **));
      for (i = 0; i < system_rank; i++)
      { // allocate for degree
        RPD->coeff_d[i] = (comp_d **)bmalloc(RPD->new_degrees[i] * sizeof(comp_d *));
        RPD->coeff_mp[i] = (comp_mp **)bmalloc(RPD->new_degrees[i] * sizeof(comp_mp *));
        for (j = 0; j < RPD->new_degrees[i]; j++)
        { // allocate for variables
          RPD->coeff_d[i][j] = (comp_d *)bmalloc(RPD->new_variables * sizeof(comp_d));
          RPD->coeff_mp[i][j] = (comp_mp *)bmalloc(RPD->new_variables * sizeof(comp_mp));
          for (k = 0; k < RPD->new_variables; k++)
          {
            init_d(RPD->coeff_d[i][j][k]);
            init_mp2(RPD->coeff_mp[i][j][k], RPD->curr_precision);
            rat_to_d(RPD->coeff_d[i][j][k], coeff_rat[i][j][k]);
            rat_to_mp(RPD->coeff_mp[i][j][k], coeff_rat[i][j][k]);
          }
        }
      }
      RPD->coeff_rat = coeff_rat;
      coeff_rat = NULL;

      RPD->sameA = sameA;
      RPD->A_d = (mat_d *)bmalloc(system_rank * sizeof(mat_d));
      RPD->A_mp = (mat_mp *)bmalloc(system_rank * sizeof(mat_mp));
      for (i = 0; i < system_rank; i++)
      {
        init_mat_d(RPD->A_d[i], A_rows[i], A_cols[i]);
        init_mat_mp2(RPD->A_mp[i], A_rows[i], A_cols[i], RPD->curr_precision);
        mat_rat_to_d(RPD->A_d[i], A_rat[i], A_rows[i], A_cols[i]);
        mat_rat_to_mp(RPD->A_mp[i], A_rat[i], A_rows[i], A_cols[i]);
      }
      RPD->A_rat = A_rat;
      A_rat = NULL;
      free(A_rows);
      free(A_cols);
    }
  }

  // now that we have the other structures, finish setting up the codim
  for (i = 0; i < num_codim; i++)
    if (i < codim_index)
    { // setup other structures - so they clear properly
      RPD->codim[i].useIntrinsicSlice = 0;
      RPD->codim[i].B_rat = NULL;
      RPD->codim[i].p_rat = NULL;
    }
    else if (i == codim_index)
    { // setup other structures
      RPD->codim[i].num_superset = RPD->codim[i].num_nonsing = RPD->codim[i].num_sing = RPD->codim[i].num_nonsolns = RPD->codim[i].num_inf = RPD->codim[i].num_bad = 0;
      if (RPD->codim[i].useIntrinsicSlice)
      { // setup the intrinsic slice
        setupRegenCodimIntrinsicSlice(RPD, MPType, i, max_prec);
      }
      else
      { // NULL out the pointers
        RPD->codim[i].B_rat = NULL;
        RPD->codim[i].p_rat = NULL;
      }
    }
    else
    { // setup other structures
      RPD->codim[i].num_superset = RPD->codim[i].num_nonsing = RPD->codim[i].num_sing = RPD->codim[i].num_nonsolns = RPD->codim[i].num_inf = RPD->codim[i].num_bad = 0;
      if (RPD->codim[i].useIntrinsicSlice)
      { // setup the intrinsic slice
        setupRegenCodimIntrinsicSlice(RPD, MPType, i, max_prec);
      }
      else
      { // NULL out the pointers
        RPD->codim[i].B_rat = NULL;
        RPD->codim[i].p_rat = NULL;
      }
    }

  // all memory should be cleared 
 
  return;
}



























