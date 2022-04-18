// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"

void setupCascadeFirstCodim(cascade_t *CD, int MPType);

void cascade_setup(FILE **OUT, char *outName, FILE **RAWOUT, char *rawName, FILE **MIDOUT, char *midName, FILE **FAIL, char *failName, tracker_config_t *T, cascade_t *CD, char *preprocFile, char *degreeFile, int maxCodim, int specificCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for positive dimensional tracking - cascade      *
\***************************************************************/
{
  int i, j, actualMaxCodim;

  // store the precision
  CD->curr_precision = T->Precision;

  // setup OUT, RAWOUT, MIDOUT & FAIL
  *OUT = fopen(outName, "w");
  *RAWOUT = fopen(rawName, "w");
  *MIDOUT = fopen(midName, "w");
  *FAIL = fopen(failName, "w");

  // setup the slp
  CD->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  CD->orig_variables = setupProg(CD->Prog, T->Precision, T->MPType);

  // error checking
  if (CD->Prog->numPathVars > 0)
  { // path variable present
    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (CD->Prog->numPars > 0)
  { // parameter present
    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup CD.PPD
  setupPreProcData(preprocFile, &CD->PPD);

  // verify that we are using only 1 homogenous variable group
  if (CD->PPD.num_hom_var_gp + CD->PPD.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: The cascade method is implemented for systems with only one variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find the rank
  if (T->MPType == 0 || T->MPType == 2)
    CD->system_rank = rank_finder_d(&CD->PPD, CD->Prog, T, CD->orig_variables);
  else
    CD->system_rank = rank_finder_mp(&CD->PPD, CD->Prog, T, CD->orig_variables);

  // setup the number of variables that will be used for tracking
  T->numVars = CD->new_variables = CD->system_rank + CD->PPD.num_var_gp + CD->PPD.num_hom_var_gp; // add on the number of patches attached

  // setup the number of functions and the number of codimensions
  CD->num_funcs = CD->PPD.num_funcs;
  CD->num_codim = CD->system_rank;
  actualMaxCodim = maxCodim > 0 ? MIN(maxCodim, CD->system_rank) : CD->system_rank;

  // error checking on specific codimension
  if (specificCodim > 0 && specificCodim > CD->system_rank)
  {
    printf("NOTE: Based on the computed system rank (%d), this system has no components of codimension %d!\n", CD->system_rank, specificCodim);
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup orig_degrees, new_degrees & P
  setupDegrees_orig_new_perm(&CD->orig_degrees, &CD->new_degrees, &CD->P, CD->num_funcs, CD->PPD.num_var_gp + CD->PPD.num_hom_var_gp, degreeFile);

  // setup gamma
  if (T->MPType == 0)
  { // only setup gamma_d
    get_comp_rand_d(CD->gamma_d);
  }
  else if (T->MPType == 1)
  { // only setup gamma_mp
    init_mp(CD->gamma_mp);
    get_comp_rand_mp(CD->gamma_mp);
  }
  else
  { // setup gamma_rat, gamma_mp & gamma_d
    get_comp_rand_rat(CD->gamma_d, CD->gamma_mp, CD->gamma_rat, CD->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // setup C, if needed
  if (CD->orig_variables != CD->new_variables)
  { // we will need to convert between the old and new variables (orig_vars = C * new_vars)
    if (T->MPType == 0)
    { // only setup C_d
      init_mat_d(CD->C_d, CD->orig_variables, CD->new_variables);
      make_matrix_random_d(CD->C_d, CD->orig_variables, CD->new_variables);
    }
    else if (T->MPType == 1)
    { // only setup C_mp
      init_mat_mp(CD->C_mp, CD->orig_variables, CD->new_variables);
      make_matrix_random_mp(CD->C_mp, CD->orig_variables, CD->new_variables, CD->curr_precision);
    }
    else
    { // allocate for C_rat
      init_mat_d(CD->C_d, CD->orig_variables, CD->new_variables); 
      init_mat_mp2(CD->C_mp, CD->orig_variables, CD->new_variables, CD->curr_precision); 
      init_mat_rat(CD->C_rat, CD->orig_variables, CD->new_variables); 

      // setup C_rat, C_mp & C_d
      make_matrix_random_rat(CD->C_d, CD->C_mp, CD->C_rat, CD->orig_variables, CD->new_variables, CD->curr_precision, T->AMP_max_prec, 0, 0);
    }
  }
  else
  { // set sizes to 0
    CD->C_d->rows = CD->C_d->cols = CD->C_mp->rows = CD->C_mp->cols = 0;
  }

  // setup A, if needed
  if (CD->num_funcs != CD->system_rank)
  { // we need to square the system
    if (T->MPType == 0)
    { // only setup A_d
      init_mat_d(CD->A_d, CD->system_rank, CD->num_funcs - CD->system_rank);
      make_matrix_random_d(CD->A_d, CD->system_rank, CD->num_funcs - CD->system_rank);
    }
    else if (T->MPType == 1)
    { // only setup A_mp
      init_mat_mp(CD->A_mp, CD->system_rank, CD->num_funcs - CD->system_rank);
      make_matrix_random_mp(CD->A_mp, CD->system_rank, CD->num_funcs - CD->system_rank, CD->curr_precision);
    }
    else
    { // allocate for A_rat
      init_mat_d(CD->A_d, CD->system_rank, CD->num_funcs - CD->system_rank);
      init_mat_mp2(CD->A_mp, CD->system_rank, CD->num_funcs - CD->system_rank, CD->curr_precision);
      init_mat_rat(CD->A_rat, CD->system_rank, CD->num_funcs - CD->system_rank);

      // setup A_rat, A_mp & A_d
      make_matrix_random_rat(CD->A_d, CD->A_mp, CD->A_rat, CD->system_rank, CD->num_funcs - CD->system_rank, CD->curr_precision, T->AMP_max_prec, 0, 0);
    }
  }
  else
  { // set sizes to 0
    CD->A_d->rows = CD->A_d->cols = CD->A_mp->rows = CD->A_mp->cols = 0;
  }

  // setup W, if needed
  if (CD->num_funcs != CD->system_rank)
  { // we need to square the system
    CD->W = (int **)bmalloc(CD->system_rank * sizeof(int *));
    for (i = 0; i < CD->system_rank; i++)
    {
      CD->W[i] = (int *)bmalloc((CD->num_funcs - CD->system_rank) * sizeof(int));
      for (j = 0; j < CD->num_funcs - CD->system_rank; j++)
        CD->W[i][j] = CD->new_degrees[i] - CD->new_degrees[CD->system_rank + j];
    }
  }

  // setup R - make square so that standard_eval works
  if (T->MPType == 0)
  { // only setup R_d
    init_mat_d(CD->R_d, CD->system_rank, CD->system_rank);
    make_matrix_random_d(CD->R_d, CD->system_rank, CD->system_rank);
  }
  else if (T->MPType == 1)
  { // only setup R_mp
    init_mat_mp(CD->R_mp, CD->system_rank, CD->system_rank);
    make_matrix_random_mp(CD->R_mp, CD->system_rank, CD->system_rank, CD->curr_precision);
  }
  else
  { // allocate for R_rat
    init_mat_d(CD->R_d, CD->system_rank, CD->system_rank);
    init_mat_mp2(CD->R_mp, CD->system_rank, CD->system_rank, CD->curr_precision);
    init_mat_rat(CD->R_rat, CD->system_rank, CD->system_rank);

    // setup R_rat, R_mp & R_d
    make_matrix_random_rat(CD->R_d, CD->R_mp, CD->R_rat, CD->system_rank, CD->system_rank, CD->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup T - initialize to all ones
  if (T->MPType == 0)
  { // only setup T_d
    init_vec_d(CD->T_d, CD->system_rank - 1);
    CD->T_d->size = CD->system_rank - 1;
    for (i = 0; i < CD->T_d->size; i++)
    {
      set_one_d(&CD->T_d->coord[i]);
    }
  }
  else if (T->MPType == 1)
  { // only setup T_mp
    init_vec_mp(CD->T_mp, CD->system_rank - 1);
    CD->T_mp->size = CD->system_rank - 1;
    for (i = 0; i < CD->T_mp->size; i++)
    {
      set_one_mp(&CD->T_mp->coord[i]);
    }
  }
  else
  { // setup both T_d & T_mp
    init_vec_d(CD->T_d, CD->system_rank - 1);
    init_vec_mp2(CD->T_mp, CD->system_rank - 1, CD->curr_precision);
    CD->T_d->size = CD->T_mp->size = CD->system_rank - 1;
    for (i = 0; i < CD->T_d->size; i++)
    {
      set_one_d(&CD->T_d->coord[i]);
      set_one_mp(&CD->T_mp->coord[i]);
    }
  }

  // setup B
  if (T->MPType == 0)
  { // only setup B_d
    init_mat_d(CD->B_d, CD->system_rank - 1, CD->new_variables);
    make_matrix_random_d(CD->B_d, CD->system_rank - 1, CD->new_variables);
  }
  else if (T->MPType == 1)
  { // only setup B_mp
    init_mat_mp(CD->B_mp, CD->system_rank - 1, CD->new_variables);
    make_matrix_random_mp(CD->B_mp, CD->system_rank - 1, CD->new_variables, CD->curr_precision);
  }
  else
  { // allocate B_rat
    init_mat_d(CD->B_d, CD->system_rank - 1, CD->new_variables);
    init_mat_mp2(CD->B_mp, CD->system_rank - 1, CD->new_variables, CD->curr_precision);
    init_mat_rat(CD->B_rat, CD->system_rank - 1, CD->new_variables);

    // setup B_rat, B_mp & B_d
    make_matrix_random_rat(CD->B_d, CD->B_mp, CD->B_rat, CD->system_rank - 1, CD->new_variables, CD->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup p
  if (T->MPType == 0)
  { // only setup p_d
    init_vec_d(CD->p_d, CD->new_variables);
    make_vec_random_d(CD->p_d, CD->new_variables);
  }
  else if (T->MPType == 1)
  { // only setup p_mp
    init_vec_mp(CD->p_mp, CD->new_variables);
    make_vec_random_mp(CD->p_mp, CD->new_variables);
  }
  else
  { // allocate p_rat
    init_vec_d(CD->p_d, CD->new_variables);
    init_vec_mp2(CD->p_mp, CD->new_variables, CD->curr_precision);
    init_vec_rat(CD->p_rat, CD->new_variables);

    // setup p_rat, p_mp & p_d
    make_vec_random_rat(CD->p_d, CD->p_mp, CD->p_rat, CD->new_variables, CD->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup W_prime
  CD->W_prime = (int *)bmalloc(CD->system_rank * sizeof(int));
  for (i = 0; i < CD->system_rank; i++)
    CD->W_prime[i] = CD->new_degrees[i] - 1; // the - 1 for the linear having degree 1

  // setup H & homVarConst
  if (CD->PPD.num_var_gp)
  { // using a variable group that was originally not homogenized

    // set H to be [1,0..0] and homVarConst to be 0 if we are in the original variables, otherwise H = first row of C so that x_0 = H * new_vars
    if (T->MPType == 0)
    { // homVarConst_d
      set_zero_d(CD->homVarConst_d);

      // H_d 
      init_vec_d(CD->H_d, CD->new_variables);
      CD->H_d->size = CD->new_variables;
      if (CD->orig_variables != CD->new_variables)
      { // first row of C
        for (i = 0; i < CD->H_d->size; i++)
        {
          set_d(&CD->H_d->coord[i], &CD->C_d->entry[0][i]);
        }
      }
      else
      { // [1,0,..,0]
        set_one_d(&CD->H_d->coord[0]);
        for (i = 1; i < CD->H_d->size; i++)
        {
          set_zero_d(&CD->H_d->coord[i]);
        }
      }
    }
    else if (T->MPType == 1)
    { // homVarConst_mp
      init_mp(CD->homVarConst_mp);
      set_zero_mp(CD->homVarConst_mp);

      // H_mp
      init_vec_mp(CD->H_mp, CD->new_variables);
      CD->H_mp->size = CD->new_variables;
      if (CD->orig_variables != CD->new_variables)
      { // first row of C
        for (i = 0; i < CD->H_mp->size; i++)
        {
          set_mp(&CD->H_mp->coord[i], &CD->C_mp->entry[0][i]);
        }
      }
      else
      { // [1,0,..,0]
        set_one_mp(&CD->H_mp->coord[0]);
        for (i = 1; i < CD->H_mp->size; i++)
        {
          set_zero_mp(&CD->H_mp->coord[i]);
        }
      }
    }
    else
    { // setup H_d, H_mp, H_rat, homVarConst_d, homVarConst_mp & homVarcConst_rat

      // setup homVarConst
      init_mp2(CD->homVarConst_mp, CD->curr_precision);
      mpq_init(CD->homVarConst_rat[0]); mpq_init(CD->homVarConst_rat[1]);

      set_zero_d(CD->homVarConst_d);
      set_zero_mp(CD->homVarConst_mp);
      set_zero_rat(CD->homVarConst_rat);

      // setup H
      init_vec_d(CD->H_d, CD->new_variables);
      init_vec_mp2(CD->H_mp, CD->new_variables, CD->curr_precision);
      init_vec_rat(CD->H_rat, CD->new_variables);
      CD->H_d->size = CD->H_mp->size = CD->new_variables;
      if (CD->orig_variables != CD->new_variables)
      { // first row of C
        for (i = 0; i < CD->new_variables; i++)
        {
          set_rat(CD->H_rat[i], CD->C_rat[0][i]);
          set_d(&CD->H_d->coord[i], &CD->C_d->entry[0][i]);
          set_mp(&CD->H_mp->coord[i], &CD->C_mp->entry[0][i]);
        }
      }
      else
      { // [1,0..,0]
        for (i = 0; i < CD->new_variables; i++)
        {
          if (i == 0)
          { // set to 1
            set_one_d(&CD->H_d->coord[i]);
            set_one_mp(&CD->H_mp->coord[i]);
            set_one_rat(CD->H_rat[i]);
          }
          else
          { // set to 0
            set_zero_d(&CD->H_d->coord[i]);
            set_zero_mp(&CD->H_mp->coord[i]);
            set_zero_rat(CD->H_rat[i]);
          }
        }
      }
    }
  }
  else
  { // using a homogeneous variable group

    // set H and homVarConst to be random
    if (T->MPType == 0)
    { // only setup H_d and homVarConst_d
      init_vec_d(CD->H_d, CD->new_variables);
      make_vec_random_d(CD->H_d, CD->new_variables);

      get_comp_rand_d(CD->homVarConst_d);
    }
    else if (T->MPType == 1)
    { // only setup H_mp and homVarConst_mp
      init_vec_mp(CD->H_mp, CD->new_variables);
      make_vec_random_mp(CD->H_mp, CD->new_variables);

      init_mp(CD->homVarConst_mp);
      get_comp_rand_mp(CD->homVarConst_mp);
    }
    else
    { // setup H_d, H_mp, H_rat, homVarConst_d, homVarConst_mp & homVarcConst_rat
      init_vec_d(CD->H_d, CD->new_variables);
      init_vec_mp2(CD->H_mp, CD->new_variables, CD->curr_precision);
      init_vec_rat(CD->H_rat, CD->new_variables);

      // setup H_rat, H_mp & H_d
      make_vec_random_rat(CD->H_d, CD->H_mp, CD->H_rat, CD->new_variables, CD->curr_precision, T->AMP_max_prec, 0, 0);

      // setup homVarConst
      get_comp_rand_rat(CD->homVarConst_d, CD->homVarConst_mp, CD->homVarConst_rat, CD->curr_precision, T->AMP_max_prec, 1, 1);
    }
  }

  // allocate for the codimensions
  CD->codim = (cascadeCodim_t *)bmalloc(actualMaxCodim * sizeof(cascadeCodim_t));

  // setup the first codimension
  setupCascadeFirstCodim(CD, T->MPType);

  // print message about codimensions
  if (specificCodim > 0)
  { // print a message
    CD->num_codim = actualMaxCodim; 
    printf("NOTE: You have requested to compute only codimension %d.\n", specificCodim);
  }
  else if (actualMaxCodim < CD->num_codim)
  { // change the number of codim and print a message 
    CD->num_codim = actualMaxCodim;
    printf("NOTE: You have requested a maximum codimension of %d.\n", actualMaxCodim);
  }

  return;
}

void setupCascadeFirstCodim(cascade_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the intial cascade codimension                   *
\***************************************************************/
{
  int i, j, num_paths;
  FILE *tempFile = NULL;

  // setup the total degree start points to 'nonhom_start'
  if (MPType == 0 || MPType == 2)
  { // setup in double precision
    TDstartMaker_d(CD->new_degrees, CD->system_rank);
  }
  else
  { // setup the start points in multi precision
    TDstartMaker_mp(CD->new_degrees, CD->system_rank);
  }

  // open up nonhom_start
  tempFile = fopen("nonhom_start", "r");
  if (tempFile == NULL)
  {
    printf("ERROR: 'nonhom_start' does not exist!!!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in the number of paths
  fscanf(tempFile, "%d\n\n", &num_paths);

  // allocate the appropriate memory
  allocateCascadeCodim(CD, 0, 1, num_paths, MPType);

  if (MPType == 0 || MPType == 2)
  { // setup the start points in double precision
    comp_d tempComp, tempComp2;
    vec_d tempVec;

    init_vec_d(tempVec, CD->new_variables);
    tempVec->size = CD->new_variables; 

    for (i = 0; i < num_paths; i++)
    { // read in the ith point 
      set_one_d(&tempVec->coord[0]);
      for (j = 1; j < tempVec->size; j++)
      { // read in the jth coordinate
        fscanf(tempFile, "%lf %lf;\n", &tempVec->coord[j].r, &tempVec->coord[j].i);
      }
      fscanf(tempFile, "\n");

      // find the patch normalizer
      set_zero_d(tempComp);
      if (CD->PPD.num_var_gp)
      { // p * vars
        for (j = 0; j < tempVec->size; j++)
        {
          sum_mul_d(tempComp, &CD->p_d->coord[j], &tempVec->coord[j]);
        }
        // tempComp = 1 / tempComp
        recip_d(tempComp, tempComp);
      }
      else
      { // (p - H) * vars
        for (j = 0; j < tempVec->size; j++)
        {
          sub_d(tempComp2, &CD->p_d->coord[j], &CD->H_d->coord[j]);
          sum_mul_d(tempComp, tempComp2, &tempVec->coord[j]);
        }
        // tempComp = homVarConst / tempComp
        div_d(tempComp, CD->homVarConst_d, tempComp);
      }

      // store the start point after it is adjusted to the patch
      change_size_vec_d(CD->codim[0].startPts_d[i], CD->new_variables);
      CD->codim[0].startPts_d[i]->size = CD->new_variables;
      for (j = 0; j < tempVec->size; j++) 
      {
        mul_d(&CD->codim[0].startPts_d[i]->coord[j], &tempVec->coord[j], tempComp);
      }
    }

    clear_vec_d(tempVec);
  }
  else
  { // setup the start points in multi precision
    comp_mp tempComp, tempComp2;
    vec_mp tempVec;
 
    init_mp(tempComp); init_mp(tempComp2);
    init_vec_mp(tempVec, CD->new_variables);
    tempVec->size = CD->new_variables;

    for (i = 0; i < num_paths; i++)
    { // read in the ith point
      set_one_mp(&tempVec->coord[0]);
      for (j = 1; j < tempVec->size; j++)
      { // read in the jth coordinate
        mpf_inp_str(tempVec->coord[j].r, tempFile, 10);
        mpf_inp_str(tempVec->coord[j].i, tempFile, 10);
        fscanf(tempFile, ";\n");
      }
      fscanf(tempFile, "\n");

      // find the patch normalizer
      set_zero_mp(tempComp);
      if (CD->PPD.num_var_gp)
      { // p * vars
        for (j = 0; j < tempVec->size; j++)
        {
          sum_mul_mp(tempComp, &CD->p_mp->coord[j], &tempVec->coord[j]);
        }
        // tempComp = 1 / tempComp
        recip_mp(tempComp, tempComp);
      }
      else
      { // (p - H) * vars
        for (j = 0; j < tempVec->size; j++)
        {
          sub_mp(tempComp2, &CD->p_mp->coord[j], &CD->H_mp->coord[j]);
          sum_mul_mp(tempComp, tempComp2, &tempVec->coord[j]);
        }
        // tempComp = homVarConst / tempComp
        div_mp(tempComp, CD->homVarConst_mp, tempComp);
      }

      // store the start point after it is adjusted to the patch
      change_size_vec_mp(CD->codim[0].startPts_mp[i], CD->new_variables);
      CD->codim[0].startPts_mp[i]->size = CD->new_variables;
      for (j = 0; j < tempVec->size; j++)
      {
        mul_mp(&CD->codim[0].startPts_mp[i]->coord[j], &tempVec->coord[j], tempComp);
      }
    }

    clear_mp(tempComp); clear_mp(tempComp2);
    clear_vec_mp(tempVec);
  }

  return;
}

void allocateCascadeCodim(cascade_t *CD, int codim_index, int codim, int num_paths, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: allocate the structures inside of 'codim_index'        *
\***************************************************************/
{
  int i;

  // setup codim & num_paths
  CD->codim[codim_index].codim = codim;
  CD->codim[codim_index].num_paths = num_paths;

  // set other counts to 0
  CD->codim[codim_index].num_superset = CD->codim[codim_index].num_nonsing = CD->codim[codim_index].num_sing = CD->codim[codim_index].num_nonsolns 
    = CD->codim[codim_index].num_inf = CD->codim[codim_index].num_other = CD->codim[codim_index].num_bad = 0;

  // allocate endPt_types
  CD->codim[codim_index].endPt_types = (int *)bmalloc(num_paths * sizeof(int));

  // allocate startPts and endPts
  if (MPType == 0)
  { // allocate startPts_d
    CD->codim[codim_index].startPts_d = (point_d *)bmalloc(num_paths * sizeof(point_d));
    // allocate endPts_d
    CD->codim[codim_index].endPts_d = (endpoint_data_d *)bmalloc(num_paths * sizeof(endpoint_data_d));

    // NULL out other structures
    CD->codim[codim_index].startPts_mp = NULL;
    CD->codim[codim_index].endPts_mp = NULL;
    CD->codim[codim_index].endPts_amp = NULL;

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_d(CD->codim[codim_index].startPts_d[i], 0);
      init_endpoint_data_d(&CD->codim[codim_index].endPts_d[i]);
    }
  }
  else if (MPType == 1)
  { // allocate startPts_mp
    CD->codim[codim_index].startPts_mp = (point_mp *)bmalloc(num_paths * sizeof(point_mp));
    // allocate endPts_mp
    CD->codim[codim_index].endPts_mp = (endpoint_data_mp *)bmalloc(num_paths * sizeof(endpoint_data_mp));

    // NULL out other structures
    CD->codim[codim_index].startPts_d = NULL;
    CD->codim[codim_index].endPts_d = NULL;
    CD->codim[codim_index].endPts_amp = NULL;

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_mp(CD->codim[codim_index].startPts_mp[i], 0);
      init_endpoint_data_mp(&CD->codim[codim_index].endPts_mp[i]);
    }
  }
  else
  { // allocate startPts_d
    CD->codim[codim_index].startPts_d = (point_d *)bmalloc(num_paths * sizeof(point_d));
    // allocate endPts_mp
    CD->codim[codim_index].endPts_amp = (endpoint_data_amp *)bmalloc(num_paths * sizeof(endpoint_data_amp));

    // NULL out other structures
    CD->codim[codim_index].startPts_mp = NULL;
    CD->codim[codim_index].endPts_d = NULL;
    CD->codim[codim_index].endPts_mp = NULL;

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_d(CD->codim[codim_index].startPts_d[i], 0);
      init_endpoint_data_amp(&CD->codim[codim_index].endPts_amp[i], 64, 64); 
    }
  }

  return;
}   

void cascade_clear(cascade_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear CD - positive dimensional tracking - cascade     *
\***************************************************************/
{
  int i;

  // clear the slp
  clearProg(CD->Prog, MPType, 0);
  free(CD->Prog);

  // clear PPD
  preproc_data_clear(&CD->PPD);

  // clear orig_degrees
  free(CD->orig_degrees);

  // clear P
  free(CD->P);

  // clear new_degrees
  free(CD->new_degrees);

  // clear gamma
  clear_d_mp_rat(CD->gamma_d, CD->gamma_mp, CD->gamma_rat, MPType);

  // clear C, if needed
  if (CD->orig_variables != CD->new_variables)
  { // clear C
    clear_mat(CD->C_d, CD->C_mp, CD->C_rat, MPType);
  }

  // clear A, if needed
  if (CD->num_funcs != CD->system_rank)
  { // clear A
    clear_mat(CD->A_d, CD->A_mp, CD->A_rat, MPType);
  }

  // clear W, if needed
  if (CD->num_funcs != CD->system_rank)
  { // we need to square the system
    for (i = CD->system_rank - 1; i >= 0; i--)
      free(CD->W[i]);
    free(CD->W);
  }

  // clear R
  clear_mat(CD->R_d, CD->R_mp, CD->R_rat, MPType);

  // clear T
  if (MPType == 0 || MPType == 2)
  { // clear T_d
    clear_vec_d(CD->T_d);
  }
  if (MPType == 1 || MPType == 2)
  { // clear T_mp
    clear_vec_mp(CD->T_mp);
  }

  // clear B
  clear_mat(CD->B_d, CD->B_mp, CD->B_rat, MPType);

  // clear p
  clear_vec(CD->p_d, CD->p_mp, CD->p_rat, MPType);

  // clear W_prime
  free(CD->W_prime);

  // clear H
  clear_vec(CD->H_d, CD->H_mp, CD->H_rat, MPType);

  // clear homVarConst
  clear_d_mp_rat(CD->homVarConst_d, CD->homVarConst_mp, CD->homVarConst_rat, MPType);

  // clear the codimensions
  for (i = CD->num_codim - 1; i >= 0; i--)
  { // clear codim i
    clearCascadeCodim(CD, i, MPType);
  }

  // free the codim
  free(CD->codim);

  return;
}

void clearCascadeCodim(cascade_t *CD, int codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears codimensin 'codim_index'                        *
\***************************************************************/
{
  int i, num_paths = CD->codim[codim_index].num_paths;

  // free endPt_types
  free(CD->codim[codim_index].endPt_types);

  if (MPType == 0)
  { // clear the memory
    for (i = num_paths - 1; i >= 0; i--) 
    {
      if (CD->codim[codim_index].startPts_d != NULL)  
      { // clear startPts_d[i]
        clear_point_d(CD->codim[codim_index].startPts_d[i]);
      }

      if (CD->codim[codim_index].endPts_d != NULL)
      { // clear endPts_d[i]
        clear_endpoint_data_d(&CD->codim[codim_index].endPts_d[i]);
      } 
    }
     // free startPts_d, if needed
    if (CD->codim[codim_index].startPts_d != NULL)
      free(CD->codim[codim_index].startPts_d);
    // free endPts_d, if needed
    if (CD->codim[codim_index].endPts_d != NULL)
      free(CD->codim[codim_index].endPts_d);
  }
  else if (MPType == 1)
  { // clear the memory
    for (i = num_paths - 1; i >= 0; i--)
    {
      if (CD->codim[codim_index].startPts_mp != NULL)
      { // clear startPts_mp[i]
        clear_point_mp(CD->codim[codim_index].startPts_mp[i]);
      }

      if (CD->codim[codim_index].endPts_mp != NULL)
      { // clear finalT
        clear_endpoint_data_mp(&CD->codim[codim_index].endPts_mp[i]);
      }
    }
    // free startPts_mp, if needed
    if (CD->codim[codim_index].startPts_mp != NULL)
      free(CD->codim[codim_index].startPts_mp);
    // free endPts_mp, if needed
    if (CD->codim[codim_index].endPts_mp != NULL)
      free(CD->codim[codim_index].endPts_mp);
  }
  else
  { // clear the memory
    for (i = num_paths - 1; i >= 0; i--)
    {
      if (CD->codim[codim_index].startPts_d != NULL)
      { // clear startPts_d[i]
        clear_point_d(CD->codim[codim_index].startPts_d[i]);
      }
     
      if (CD->codim[codim_index].endPts_amp != NULL)
      { // clear endPts_amp[i]
        clear_endpoint_data_amp(&CD->codim[codim_index].endPts_amp[i]);
      }
    }
    // free startPts_d, if needed
    if (CD->codim[codim_index].startPts_d != NULL)
      free(CD->codim[codim_index].startPts_d);
    // free endPts_amp, if needed
    if (CD->codim[codim_index].endPts_amp != NULL)
      free(CD->codim[codim_index].endPts_amp);
  }  

  return;
}

void cascade_clear_start_points(cascade_t *CD, int codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears the start points from 'codim_index'             *
\***************************************************************/
{
  int i, num_paths = CD->codim[codim_index].num_paths;

  if (MPType == 0 || MPType == 2)
  { // clear the start points
    for (i = num_paths - 1; i >= 0; i--)
      clear_point_d(CD->codim[codim_index].startPts_d[i]);

    // free the start points
    free(CD->codim[codim_index].startPts_d);
    // NULL out _d
    CD->codim[codim_index].startPts_d = NULL;
  }
  else if (MPType == 1)
  { // clear the start points
    for (i = num_paths - 1; i >= 0; i--)
      clear_point_mp(CD->codim[codim_index].startPts_mp[i]);

    // free the start points
    free(CD->codim[codim_index].startPts_mp);
    // NULL out _mp
    CD->codim[codim_index].startPts_mp = NULL;
  }

  return;
}

////// OPENMP COPY FUNCTIONS ///////

void setup_omp_cascade_t(cascade_t *CD, cascade_t *CD_in, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copies CD_in to CD for use in OpenMP tracking          *
\***************************************************************/
{
  // setup the current precision
  CD->curr_precision = CD_in->curr_precision;

  // setup a new prog
  CD->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(CD->Prog, CD_in->Prog);

  // setup PPD
  CD->PPD.num_funcs = CD_in->PPD.num_funcs;
  CD->PPD.num_hom_var_gp = CD_in->PPD.num_hom_var_gp;
  CD->PPD.num_var_gp = CD_in->PPD.num_var_gp;
  CD->PPD.type = CD_in->PPD.type;
  CD->PPD.size = CD_in->PPD.size;

  // point to other structures
  CD->orig_degrees = CD_in->orig_degrees;
  CD->new_degrees = CD_in->new_degrees;
  CD->P = CD_in->P;

  // copy other values
  CD->system_rank = CD_in->system_rank;
  CD->orig_variables = CD_in->orig_variables;
  CD->new_variables = CD_in->new_variables;

  // setup gamma
  if (MPType == 0 || MPType == 2)
  { // copy gamma_d
    set_d(CD->gamma_d, CD_in->gamma_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy gamma_mp
    init_mp2(CD->gamma_mp, CD->curr_precision);
    set_mp(CD->gamma_mp, CD_in->gamma_mp);
  }

  if (MPType == 2)
  { // setup gamma_rat
    init_rat(CD->gamma_rat);
    set_rat(CD->gamma_rat, CD_in->gamma_rat);
  }

  // setup C, if needed
  if (CD->orig_variables != CD->new_variables)
  {
    if (MPType == 0 || MPType == 2)
    { // copy C_d
      init_mat_d(CD->C_d, CD_in->C_d->rows, CD_in->C_d->cols);
      mat_cp_d(CD->C_d, CD_in->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // copy C_mp
      init_mat_mp2(CD->C_mp, CD_in->C_mp->rows, CD_in->C_mp->cols, CD->curr_precision);
      mat_cp_mp(CD->C_mp, CD_in->C_mp);
    }

    if (MPType == 2)
    { // point to C_rat
      CD->C_rat = CD_in->C_rat;
    }
  }

  // copy other values
  CD->num_funcs = CD_in->num_funcs;
  CD->num_codim = CD_in->num_codim;
  CD->curr_codim_index = CD_in->curr_codim_index;

  // setup A, if needed
  if (CD->num_funcs != CD->system_rank)
  { // we need to square the system
    if (MPType == 0 || MPType == 2)
    { // copy A_d
      init_mat_d(CD->A_d, CD_in->A_d->rows, CD_in->A_d->cols);
      mat_cp_d(CD->A_d, CD_in->A_d);
    }

    if (MPType == 1 || MPType == 2)
    { // copy A_mp
      init_mat_mp2(CD->A_mp, CD_in->A_mp->rows, CD_in->A_mp->cols, CD->curr_precision);
      mat_cp_mp(CD->A_mp, CD_in->A_mp);
    }

    if (MPType == 2)
    { // point to A_rat
      CD->A_rat = CD_in->A_rat;
    }
  }

  // point to W
  CD->W = CD_in->W;

  // setup R
  if (MPType == 0 || MPType == 2)
  { // copy R_d
    init_mat_d(CD->R_d, CD_in->R_d->rows, CD_in->R_d->cols);
    mat_cp_d(CD->R_d, CD_in->R_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy R_mp
    init_mat_mp2(CD->R_mp, CD_in->R_mp->rows, CD_in->R_mp->cols, CD->curr_precision);
    mat_cp_mp(CD->R_mp, CD_in->R_mp);
  }

  if (MPType == 2)
  { // point to R_rat
    CD->R_rat = CD_in->R_rat;
  }

  // setup T
  if (MPType == 0 || MPType == 2)
  { // copy T_d
    init_vec_d(CD->T_d, CD_in->T_d->size);
    vec_cp_d(CD->T_d, CD_in->T_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy T_mp
    init_vec_mp2(CD->T_mp, CD_in->T_mp->size, CD->curr_precision);
    vec_cp_mp(CD->T_mp, CD_in->T_mp);
  }

  // setup B
  if (MPType == 0 || MPType == 2)
  { // copy B_d
    init_mat_d(CD->B_d, CD_in->B_d->rows, CD_in->B_d->cols);
    mat_cp_d(CD->B_d, CD_in->B_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy B_mp
    init_mat_mp2(CD->B_mp, CD_in->B_mp->rows, CD_in->B_mp->cols, CD->curr_precision);
    mat_cp_mp(CD->B_mp, CD_in->B_mp);
  }

  if (MPType == 2)
  { // point to B_rat
    CD->B_rat = CD_in->B_rat;
  }

  // setup p
  if (MPType == 0 || MPType == 2)
  { // copy p_d
    init_vec_d(CD->p_d, CD_in->p_d->size);
    vec_cp_d(CD->p_d, CD_in->p_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy p_mp
    init_vec_mp2(CD->p_mp, CD_in->p_mp->size, CD->curr_precision);
    vec_cp_mp(CD->p_mp, CD_in->p_mp);
  }

  if (MPType == 2)
  { // point to p_rat
    CD->p_rat = CD_in->p_rat;
  }

  // setup W_prime
  CD->W_prime = CD_in->W_prime;

  // setup H
  if (MPType == 0 || MPType == 2)
  { // copy H_d
    init_vec_d(CD->H_d, CD_in->H_d->size);
    vec_cp_d(CD->H_d, CD_in->H_d);
  }
  
  if (MPType == 1 || MPType == 2)
  { // copy H_mp
    init_vec_mp2(CD->H_mp, CD_in->H_mp->size, CD->curr_precision);
    vec_cp_mp(CD->H_mp, CD_in->H_mp);
  }
  
  if (MPType == 2)
  { // point to H_rat
    CD->H_rat = CD_in->H_rat;
  }

  // setup homVarConst
  if (MPType == 0 || MPType == 2)
  { // copy homVarConst_d
    set_d(CD->homVarConst_d, CD_in->homVarConst_d);
  }

  if (MPType == 1 || MPType == 2)
  { // copy homVarConst_mp
    init_mp2(CD->homVarConst_mp, CD->curr_precision);
    set_mp(CD->homVarConst_mp, CD_in->homVarConst_mp);
  }

  if (MPType == 2)
  { // copy homVarConst_rat
    init_rat(CD->homVarConst_rat);
    set_rat(CD->homVarConst_rat, CD_in->homVarConst_rat);
  }
   
  // since there is nothing of interst in the codimensions, just point to them
  CD->codim = CD_in->codim;

  return;
}

void clear_omp_cascade_t(cascade_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears CD where is was used in OpenMP tracking         *
\***************************************************************/
{
  // clear Prog
  clearProg(CD->Prog, MPType, 0);
  free(CD->Prog);

  // clear PPD
  CD->PPD.type = NULL;
  CD->PPD.size = NULL;

  // clear pointers to other structures
  CD->orig_degrees = NULL;
  CD->new_degrees = NULL;
  CD->P = NULL;

  // clear gamma
  clear_d_mp_rat(CD->gamma_d, CD->gamma_mp, CD->gamma_rat, MPType);

  // clear C, if used
  if (CD->orig_variables != CD->new_variables)
  {
    if (MPType == 0 || MPType == 2)
    { // clear C_d
      clear_mat_d(CD->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // clear C_mp
      clear_mat_mp(CD->C_mp);
    }

    if (MPType == 2)
    { // clear C_rat
      CD->C_rat = NULL;
    }
  }

  // clear A, if needed
  if (CD->num_funcs != CD->system_rank)
  { 
    if (MPType == 0 || MPType == 2)
    { // clear A_d
      clear_mat_d(CD->C_d);
    }

    if (MPType == 1 || MPType == 2)
    { // clear A_mp
      clear_mat_mp(CD->A_mp);
    }

    if (MPType == 2)
    { // clear A_rat
      CD->A_rat = NULL;
    }
  }

  // clear to W
  CD->W = NULL;

  // clear R
  if (MPType == 0 || MPType == 2)
  { // clear R_d
    clear_mat_d(CD->R_d);
  }

  if (MPType == 1 || MPType == 2)
  { // clear R_mp
    clear_mat_mp(CD->R_mp);
  }

  if (MPType == 2)
  { // clear R_rat
    CD->R_rat = NULL;
  }

  // clear T
  if (MPType == 0 || MPType == 2)
  { // clear T_d
    clear_vec_d(CD->T_d);
  }

  if (MPType == 1 || MPType == 2)
  { // clear T_mp
    clear_vec_mp(CD->T_mp);
  }

  // clear B
  if (MPType == 0 || MPType == 2)
  { // clear B_d
    clear_mat_d(CD->B_d);
  }

  if (MPType == 1 || MPType == 2)
  { // clear B_mp
    clear_mat_mp(CD->B_mp);
  }

  if (MPType == 2)
  { // clear B_rat
    CD->B_rat = NULL;
  }

  // clear p
  if (MPType == 0 || MPType == 2)
  { // clear p_d
    clear_vec_d(CD->p_d);
  }

  if (MPType == 1 || MPType == 2)
  { // clear p_mp
    clear_vec_mp(CD->p_mp);
  }

  if (MPType == 2)
  { // clear p_rat
    CD->p_rat = NULL;
  }

  // clear W_prime
  CD->W_prime = NULL;

  // clear H
  if (MPType == 0 || MPType == 2)
  { // clear H_d
    clear_vec_d(CD->H_d);
  }

  if (MPType == 1 || MPType == 2)
  { // clear H_mp
    clear_vec_mp(CD->H_mp);
  }

  if (MPType == 2)
  { // clear H_rat
    CD->H_rat = NULL;
  }

  // clear homVarConst
  clear_d_mp_rat(CD->homVarConst_d, CD->homVarConst_mp, CD->homVarConst_rat, MPType);

  // clear codim
  CD->codim = NULL;

  return;
}

///// CHANGE PRECISION /////

int change_cascade_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for cascade                           *
\***************************************************************/
{
  int i, j;

  // cast ED as CD
  cascade_t *CD = (cascade_t *)ED;

  // set the SLP to the correct precision
  CD->Prog->precision = prec;

  if (prec != CD->curr_precision)
  { // need to change the precision
    CD->curr_precision = prec;

    // change the precision for gamma
    setprec_mp(CD->gamma_mp, prec);
    mpf_set_q(CD->gamma_mp->r, CD->gamma_rat[0]);
    mpf_set_q(CD->gamma_mp->i, CD->gamma_rat[1]);

    // change the precision for C, if needed
    if (CD->new_variables != CD->orig_variables)
    {
      for (i = 0; i < CD->C_mp->rows; i++)
        for (j = 0; j < CD->C_mp->cols; j++)
        {
          setprec_mp(&CD->C_mp->entry[i][j], prec);
          mpf_set_q(CD->C_mp->entry[i][j].r, CD->C_rat[i][j][0]);
          mpf_set_q(CD->C_mp->entry[i][j].i, CD->C_rat[i][j][1]);
        }
    }

    // change the precision for A, if needed
    if (CD->num_funcs != CD->system_rank)
    {
      for (i = 0; i < CD->A_mp->rows; i++)
        for (j = 0; j < CD->A_mp->cols; j++)
        {
          setprec_mp(&CD->A_mp->entry[i][j], prec);
          mpf_set_q(CD->A_mp->entry[i][j].r, CD->A_rat[i][j][0]);
          mpf_set_q(CD->A_mp->entry[i][j].i, CD->A_rat[i][j][1]);
        }
    }

    // change the precision for R
    for (i = 0; i < CD->R_mp->rows; i++)
      for (j = 0; j < CD->R_mp->cols; j++)
      {
        setprec_mp(&CD->R_mp->entry[i][j], prec);
        mpf_set_q(CD->R_mp->entry[i][j].r, CD->R_rat[i][j][0]);
        mpf_set_q(CD->R_mp->entry[i][j].i, CD->R_rat[i][j][1]);
      }

    // change the precision for T
    for (i = 0; i < CD->T_mp->size; i++)
    {
      change_prec_mp(&CD->T_mp->coord[i], prec);
    }

    // change the precision for B
    for (i = 0; i < CD->B_mp->rows; i++)
      for (j = 0; j < CD->B_mp->cols; j++)
      {
        setprec_mp(&CD->B_mp->entry[i][j], prec);
        mpf_set_q(CD->B_mp->entry[i][j].r, CD->B_rat[i][j][0]);
        mpf_set_q(CD->B_mp->entry[i][j].i, CD->B_rat[i][j][1]);
      }

    // change the precision for p
    for (i = 0; i < CD->p_mp->size; i++)
    {
      setprec_mp(&CD->p_mp->coord[i], prec);
      mpf_set_q(CD->p_mp->coord[i].r, CD->p_rat[i][0]);
      mpf_set_q(CD->p_mp->coord[i].i, CD->p_rat[i][1]);
    }

    // change the precision for H
    for (i = 0; i < CD->H_mp->size; i++)
    {
      setprec_mp(&CD->H_mp->coord[i], prec);
      mpf_set_q(CD->H_mp->coord[i].r, CD->H_rat[i][0]);
      mpf_set_q(CD->H_mp->coord[i].i, CD->H_rat[i][1]);
    }

    // change the precision for homVarConst
    setprec_mp(CD->homVarConst_mp, prec);
    mpf_set_q(CD->homVarConst_mp->r, CD->homVarConst_rat[0]);
    mpf_set_q(CD->homVarConst_mp->i, CD->homVarConst_rat[1]);
  }

  return 0;
}

