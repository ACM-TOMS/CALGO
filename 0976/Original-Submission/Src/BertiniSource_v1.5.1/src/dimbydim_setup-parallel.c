// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"
#include "pos_dim.h"
#include "dimbydim.h"

// provides the functions to setup for using dimension-by-dimension witness superset generation
void setupCodimData(codim_t *CD, int codim_index, int codim, int MPType, int max_prec, int intrinsicCutoff);
void setupCodimStartPts_d(codim_t *CD, int codim_index);
void setupCodimStartPts_mp(codim_t *CD, int codim_index);

void dimbydim_setup(FILE **OUT, char *outName, FILE **RAWOUT, char *rawName, FILE **MIDOUT, char *midName, FILE **FAIL, char *failName, tracker_config_t *T, codim_t *CD, char *preprocFile, char *degreeFile, double intrinsicCutoffMultiplier, int maxCodim, int specificCodim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for positive dimensional tracking - dim-by-dim   *
\***************************************************************/
{
  int i, intrinsicCutoff, actualMaxCodim;

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
    printf("ERROR: The dimension-by-dimension method is implemented for systems with only one variable group.\n");
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

  // determine where to switch from intrinsic to extrinsic
  intrinsicCutoff = floor(intrinsicCutoffMultiplier * CD->new_variables);

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

  if (specificCodim > 0)
  { // only setup the correct codimension
    CD->codim = (codimData_t *)bmalloc(1 * sizeof(codimData_t));
    setupCodimData(CD, 0, specificCodim, T->MPType, T->AMP_max_prec, intrinsicCutoff);
    CD->num_codim = 1;
  }
  else
  { // allocate codim
    CD->codim = (codimData_t *)bmalloc(actualMaxCodim * sizeof(codimData_t));

    // setup the codimensions
    for (i = 0; i < actualMaxCodim; i++)
    { // setup the ith codimData
      setupCodimData(CD, i, i + 1, T->MPType, T->AMP_max_prec, intrinsicCutoff);
    }

    CD->num_codim = actualMaxCodim;
  }

  // print message about codimensions
  if (specificCodim > 0)
  { // print a message
    printf("NOTE: You have requested to compute only codimension %d.\n", specificCodim);
  }
  else if (actualMaxCodim < CD->num_codim)
  { // change the number of codim and print a message 
    printf("NOTE: You have requested a maximum codimension of %d.\n", actualMaxCodim);
  }

  return;
}

void setupCodimData(codim_t *CD, int codim_index, int codim, int MPType, int max_prec, int intrinsicCutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for the codim                                    *
\***************************************************************/
{
  int i, j, num_paths;

  // setup codim
  CD->codim[codim_index].codim = codim;
  
  // initialize the counts
  CD->codim[codim_index].num_superset = CD->codim[codim_index].num_nonsing = CD->codim[codim_index].num_sing = CD->codim[codim_index].num_nonsolns
      = CD->codim[codim_index].num_inf = CD->codim[codim_index].num_bad = 0;

  // find the number of paths for this codim
  num_paths = 1;
  for (i = 0; i < codim; i++)
  {
    num_paths *= CD->new_degrees[i];
  }

  // setup the number of paths
  CD->codim[codim_index].num_paths = num_paths;

  // see if using intrinsic slice
  if (codim > intrinsicCutoff)
  { // use extrinsic slice
    CD->codim[codim_index].useIntrinsicSlice = 0;
  }
  else
  { // use intrinsic slice
    CD->codim[codim_index].useIntrinsicSlice = 1;
  }

  // setup endPt_types
  CD->codim[codim_index].endPt_types = (int *)bmalloc(num_paths * sizeof(int));

  // setup W
  CD->codim[codim_index].W = (int **)bmalloc(codim * sizeof(int *));
  for (i = 0; i < codim; i++)
  { 
    CD->codim[codim_index].W[i] = (int *)bmalloc((CD->num_funcs - codim) * sizeof(int));
    for (j = 0; j < CD->num_funcs - codim; j++)
      CD->codim[codim_index].W[i][j] = CD->new_degrees[i] - CD->new_degrees[j+codim];
  }

  if (MPType == 0)
  { // setup the _d structures

    // setup H_d & homVarConst_d
    init_vec_d(CD->codim[codim_index].H_d, CD->new_variables);
    CD->codim[codim_index].H_d->size = CD->new_variables;
    if (CD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (CD->orig_variables != CD->new_variables)
      { // H_d is first row of C
        for (i = 0; i < CD->new_variables; i++)
        {
          set_d(&CD->codim[codim_index].H_d->coord[i], &CD->C_d->entry[0][i]);
        }
      }
      else
      { // H_d = [1,0..0]
        set_one_d(&CD->codim[codim_index].H_d->coord[0]);
        for (i = 1; i < CD->new_variables; i++)
        {
          set_zero_d(&CD->codim[codim_index].H_d->coord[i]);
        }
      }

      // setup homVarConst_d to be 0
      set_zero_d(CD->codim[codim_index].homVarConst_d);
    }
    else
    { // using a homogeneous variable group

      // setup H_d to be random
      make_vec_random_d(CD->codim[codim_index].H_d, CD->new_variables);

      // setup homVarConst_d to be random
      get_comp_rand_d(CD->codim[codim_index].homVarConst_d);
    }

    // setup A_d
    init_mat_d(CD->codim[codim_index].A_d, codim, CD->num_funcs - codim);
    make_matrix_random_d(CD->codim[codim_index].A_d, codim, CD->num_funcs - codim);

    // setup B_d
    if (CD->codim[codim_index].useIntrinsicSlice)
    { // setup for intrinsic slicing
      init_mat_d(CD->codim[codim_index].B_d, CD->new_variables, codim);
      make_matrix_random_d(CD->codim[codim_index].B_d, CD->new_variables, codim);
    }
    else
    { // setup for extrinsic slicing
      init_mat_d(CD->codim[codim_index].B_d, CD->new_variables - codim - 1, CD->new_variables);
      make_matrix_random_d(CD->codim[codim_index].B_d, CD->new_variables - codim - 1, CD->new_variables);
    }

    // adjust B if using intrinsic slice to the form [[I][B]]
    if (CD->codim[codim_index].useIntrinsicSlice)
    {
      for (i = 0; i < codim; i++)
        for (j = 0; j < codim; j++)
          if (i == j)
          {
            set_one_d(&CD->codim[codim_index].B_d->entry[i][j]);
          }
          else
          {
            set_zero_d(&CD->codim[codim_index].B_d->entry[i][j]);
          }
    }

    // setup p_d
    init_vec_d(CD->codim[codim_index].p_d, CD->new_variables);
    make_vec_random_d(CD->codim[codim_index].p_d, CD->new_variables);

    // allocate startPts_d
    CD->codim[codim_index].startPts_d = (point_d *)bmalloc(num_paths * sizeof(point_d));

    // allocate endPts_d
    CD->codim[codim_index].endPts_d = (endpoint_data_d *)bmalloc(num_paths * sizeof(endpoint_data_d));

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_d(CD->codim[codim_index].startPts_d[i], 0);
      init_endpoint_data_d(&CD->codim[codim_index].endPts_d[i]);
    }

    // setup other pointers to NULL
    CD->codim[codim_index].startPts_mp = NULL;
    CD->codim[codim_index].endPts_mp = NULL;
    CD->codim[codim_index].endPts_amp = NULL;
  }
  else if (MPType == 1)
  { // setup the _mp structures

    // initialize H_mp & homVarConst_mp
    init_vec_mp(CD->codim[codim_index].H_mp, CD->new_variables);
    CD->codim[codim_index].H_mp->size = CD->new_variables;
    init_mp(CD->codim[codim_index].homVarConst_mp);

    // setup H_mp & homVarConst_mp
    if (CD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (CD->orig_variables != CD->new_variables)
      { // H_mp is first row of C
        for (i = 0; i < CD->new_variables; i++)
        {
          set_mp(&CD->codim[codim_index].H_mp->coord[i], &CD->C_mp->entry[0][i]);
        }
      }
      else
      { // H_mp = [1,0..0]
        set_one_mp(&CD->codim[codim_index].H_mp->coord[0]);
        for (i = 1; i < CD->new_variables; i++)
        {
          set_zero_mp(&CD->codim[codim_index].H_mp->coord[i]);
        }
      }

      // setup homVarConst_mp to be 0
      set_zero_mp(CD->codim[codim_index].homVarConst_mp);
    }
    else
    { // using a homogeneous variable group

      // setup H_mp to be random
      make_vec_random_mp(CD->codim[codim_index].H_mp, CD->new_variables);

      // setup homVarConst_mp to be random
      get_comp_rand_mp(CD->codim[codim_index].homVarConst_mp);
    }

    // setup A_mp
    init_mat_mp(CD->codim[codim_index].A_mp, codim, CD->num_funcs - codim);
    make_matrix_random_mp(CD->codim[codim_index].A_mp, codim, CD->num_funcs - codim, CD->curr_precision);

    // setup B_mp
    if (CD->codim[codim_index].useIntrinsicSlice)
    { // setup for intrinsic slicing
      init_mat_mp(CD->codim[codim_index].B_mp, CD->new_variables, codim);
      make_matrix_random_mp(CD->codim[codim_index].B_mp, CD->new_variables, codim, CD->curr_precision);
    }
    else
    { // setup for extrinsic slicing
      init_mat_mp(CD->codim[codim_index].B_mp, CD->new_variables - codim - 1, CD->new_variables);
      make_matrix_random_mp(CD->codim[codim_index].B_mp, CD->new_variables - codim - 1, CD->new_variables, CD->curr_precision);
    }

    // adjust B if using intrinsic slice to the form [[I][B]]
    if (CD->codim[codim_index].useIntrinsicSlice)
    {
      for (i = 0; i < codim; i++)
        for (j = 0; j < codim; j++)
          if (i == j)
          {
            set_one_mp(&CD->codim[codim_index].B_mp->entry[i][j]);
          }
          else
          {
            set_zero_mp(&CD->codim[codim_index].B_mp->entry[i][j]);
          }
    }

    // seutp p_mp
    init_vec_mp(CD->codim[codim_index].p_mp, CD->new_variables);
    make_vec_random_mp(CD->codim[codim_index].p_mp, CD->new_variables);

    // allocate startPts_mp
    CD->codim[codim_index].startPts_mp = (point_mp *)bmalloc(num_paths * sizeof(point_mp));

    // allocate endPts_mp
    CD->codim[codim_index].endPts_mp = (endpoint_data_mp *)bmalloc(num_paths * sizeof(endpoint_data_mp));

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_mp(CD->codim[codim_index].startPts_mp[i], 0);
      init_endpoint_data_mp(&CD->codim[codim_index].endPts_mp[i]); 
    }

    // setup other pointers to NULL
    CD->codim[codim_index].startPts_d = NULL;
    CD->codim[codim_index].endPts_d = NULL;
    CD->codim[codim_index].endPts_amp = NULL;
  }
  else // setup for AMP
  { // setup _d, _mp & _rat structures

    // initialize H & homVarConst
    init_vec_d(CD->codim[codim_index].H_d, CD->new_variables);
    init_vec_mp2(CD->codim[codim_index].H_mp, CD->new_variables, CD->curr_precision);
    init_vec_rat(CD->codim[codim_index].H_rat, CD->new_variables);
    CD->codim[codim_index].H_d->size = CD->codim[codim_index].H_mp->size = CD->new_variables;

    init_mp2(CD->codim[codim_index].homVarConst_mp, CD->curr_precision);
    init_rat(CD->codim[codim_index].homVarConst_rat);

    // setup H & homVarConst
    if (CD->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized
      if (CD->orig_variables != CD->new_variables)
      { // H is first row of C
        for (i = 0; i < CD->new_variables; i++)
        {
          mpq_set(CD->codim[codim_index].H_rat[i][0], CD->C_rat[0][i][0]);
          mpq_set(CD->codim[codim_index].H_rat[i][1], CD->C_rat[0][i][1]);

          set_d(&CD->codim[codim_index].H_d->coord[i], &CD->C_d->entry[0][i]);
          set_mp(&CD->codim[codim_index].H_mp->coord[i], &CD->C_mp->entry[0][i]);
        }
      }
      else
      { // H = [1,0..0]
        set_one_d(&CD->codim[codim_index].H_d->coord[0]);
        set_one_mp(&CD->codim[codim_index].H_mp->coord[0]);
        set_one_rat(CD->codim[codim_index].H_rat[0]);
        for (i = 1; i < CD->new_variables; i++)
        {
          set_zero_d(&CD->codim[codim_index].H_d->coord[i]);
          set_zero_mp(&CD->codim[codim_index].H_mp->coord[i]);
          set_zero_rat(CD->codim[codim_index].H_rat[i]);
        }
      }

      // setup homVarConst to be 0
      set_zero_d(CD->codim[codim_index].homVarConst_d);
      set_zero_mp(CD->codim[codim_index].homVarConst_mp);
      set_zero_rat(CD->codim[codim_index].homVarConst_rat);
    }
    else
    { // using a homogeneous variable group

      // setup H to be random
      make_vec_random_rat(CD->codim[codim_index].H_d, CD->codim[codim_index].H_mp, CD->codim[codim_index].H_rat, CD->new_variables, CD->curr_precision, max_prec, 0, 0);

      // setup homVarConst to be random
      get_comp_rand_rat(CD->codim[codim_index].homVarConst_d, CD->codim[codim_index].homVarConst_mp, CD->codim[codim_index].homVarConst_rat, CD->curr_precision, max_prec, 0, 0);
    }

    // setup A_d, A_mp, A_rat
    init_mat_d(CD->codim[codim_index].A_d, codim, CD->num_funcs - codim);
    init_mat_mp2(CD->codim[codim_index].A_mp, codim, CD->num_funcs - codim, CD->curr_precision);
    init_mat_rat(CD->codim[codim_index].A_rat, codim, CD->num_funcs - codim);
    make_matrix_random_rat(CD->codim[codim_index].A_d, CD->codim[codim_index].A_mp, CD->codim[codim_index].A_rat, codim, CD->num_funcs - codim, CD->curr_precision, max_prec, 0, 0);

    // setup B
    if (CD->codim[codim_index].useIntrinsicSlice)
    { // setup for intrinsic slicing
      init_mat_d(CD->codim[codim_index].B_d, CD->new_variables, codim);
      init_mat_mp2(CD->codim[codim_index].B_mp, CD->new_variables, codim, CD->curr_precision);
      init_mat_rat(CD->codim[codim_index].B_rat, CD->new_variables, codim);
      // setup B_d, B_mp & B_rat
      make_matrix_random_rat(CD->codim[codim_index].B_d, CD->codim[codim_index].B_mp, CD->codim[codim_index].B_rat, CD->new_variables, codim, CD->curr_precision, max_prec, 0, 0);
    }
    else
    { // setup for extrinsic slicing
      init_mat_d(CD->codim[codim_index].B_d, CD->new_variables - codim - 1, CD->new_variables);
      init_mat_mp2(CD->codim[codim_index].B_mp, CD->new_variables - codim - 1, CD->new_variables, CD->curr_precision);
      init_mat_rat(CD->codim[codim_index].B_rat, CD->new_variables - codim - 1, CD->new_variables);
      // setup B_d, B_mp & B_rat
      make_matrix_random_rat(CD->codim[codim_index].B_d, CD->codim[codim_index].B_mp, CD->codim[codim_index].B_rat, CD->new_variables - codim - 1, CD->new_variables, CD->curr_precision, max_prec, 0, 0);
    }

    // adjust B if using intrinsic slice to the form [[I][B]]
    if (CD->codim[codim_index].useIntrinsicSlice)
    {
      for (i = 0; i < codim; i++)
        for (j = 0; j < codim; j++)
        {
          if (i == j) 
          { 
            set_one_rat(CD->codim[codim_index].B_rat[i][j]); 
          }
          else 
          { 
            set_zero_rat(CD->codim[codim_index].B_rat[i][j]); 
          }
          mpf_set_q(CD->codim[codim_index].B_mp->entry[i][j].r, CD->codim[codim_index].B_rat[i][j][0]);
          mpf_set_q(CD->codim[codim_index].B_mp->entry[i][j].i, CD->codim[codim_index].B_rat[i][j][1]);
          CD->codim[codim_index].B_d->entry[i][j].r = mpq_get_d(CD->codim[codim_index].B_rat[i][j][0]);
          CD->codim[codim_index].B_d->entry[i][j].i = mpq_get_d(CD->codim[codim_index].B_rat[i][j][1]);
        }
    }

    // allocate for p_rat
    init_vec_d(CD->codim[codim_index].p_d, CD->new_variables);
    init_vec_mp2(CD->codim[codim_index].p_mp, CD->new_variables, CD->curr_precision);
    init_vec_rat(CD->codim[codim_index].p_rat, CD->new_variables);
    // setup p_d, p_mp & p_rat
    make_vec_random_rat(CD->codim[codim_index].p_d, CD->codim[codim_index].p_mp, CD->codim[codim_index].p_rat, CD->new_variables, CD->curr_precision, max_prec, 0, 0);

    // allocate startPts_d
    CD->codim[codim_index].startPts_d = (point_d *)bmalloc(num_paths * sizeof(point_d));

    // allocate endPts_amp
    CD->codim[codim_index].endPts_amp = (endpoint_data_amp *)bmalloc(num_paths * sizeof(endpoint_data_amp));    

    // initialize memory
    for (i = 0; i < num_paths; i++)
    {
      init_point_d(CD->codim[codim_index].startPts_d[i], 0);
      init_endpoint_data_amp(&CD->codim[codim_index].endPts_amp[i], 64, 64);
    }

    // setup other pointers to NULL
    CD->codim[codim_index].startPts_mp = NULL;
    CD->codim[codim_index].endPts_d = NULL;
    CD->codim[codim_index].endPts_mp = NULL;
  }

  // setup start points for this codimension
  if (MPType == 0 || MPType == 2)
    setupCodimStartPts_d(CD, codim_index);
  else
    setupCodimStartPts_mp(CD, codim_index);

  return;
}

void setupCodimStartPts_d(codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for the start points for 'codim_index'           *
\***************************************************************/
{
  int i, j, num_paths = CD->codim[codim_index].num_paths, codim = CD->codim[codim_index].codim;
  int *curr_degs = (int *)bmalloc(codim * sizeof(int));
  mat_d A;
  vec_d b;

  init_mat_d(A, 0, 0);
  init_vec_d(b, 0);

  // the start system has the form [x_1^d_1 - x_0^d_0, .., x_codim^d_codim - x_0^d_codim] 
  // with extrinsic slices: B*x = 0 and patch p*x = 1 
  // with intrinsic slices: x = p + B*y, where y is tracking variables

  // initialize curr_degs
  for (i = 0; i < codim; i++)
    curr_degs[i] = 0;

  if (CD->codim[codim_index].useIntrinsicSlice)
  { // setup the start points for intrinsic tracking

    // setup A to originally be [B**,-1] ** - of the correct size (the extra column is for the homogenenous multiplier since we have variables y_1, .. y_codim, lambda)
    change_size_mat_d(A, codim + 1, codim + 1);
    A->rows = A->cols = codim + 1;
    for (i = 0; i < A->rows; i++)
    { // copy over B to the first 'codim' columns
      for (j = 0; j < codim; j++)
      {
        set_d(&A->entry[i][j], &CD->codim[codim_index].B_d->entry[i][j]);
      }
      // set last column to [-1,..-1]
      set_one_d(&A->entry[i][codim]);
      neg_d(&A->entry[i][codim], &A->entry[i][codim]);
    }

    // setup b to be [-p_0, .. -p_codim]
    change_size_vec_d(b, codim + 1);
    b->size = codim + 1;
    for (i = 0; i < b->size; i++)
    {
      neg_d(&b->coord[i], &CD->codim[codim_index].p_d->coord[i]);
    }

    // loop through to find the start points
    for (i = 0; i < num_paths; i++)
    { // setup A for this start point
      for (j = 0; j < codim; j++)
      {
        set_double_d(&A->entry[j+1][codim], -cos(2 * M_PI * curr_degs[j] / CD->new_degrees[j]), -sin(2 * M_PI * curr_degs[j] / CD->new_degrees[j]));
      }

      // solve for the start point
      if (matrixSolve_d(CD->codim[codim_index].startPts_d[i], A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      // remove the bottom coordinate (the homogeneous multiplier)
      CD->codim[codim_index].startPts_d[i]->size -= 1;

      // update curr_degs
      for (j = 0; j < codim; j++)
        if (curr_degs[j] + 1 == CD->new_degrees[j])
        { // set back to 0
          curr_degs[j] = 0;
        }
        else
        { // increment and then exit loop
          curr_degs[j]++;
          j = codim;
        }
    }
  }
  else
  { // setup the start point for extrinsic tracking

    // setup A to originally be [p][0, -I, 0][B]
    change_size_mat_d(A, CD->new_variables, CD->new_variables);
    A->rows = A->cols = CD->new_variables;
    for (i = 0; i < A->rows; i++)
      if (i < 1)
      { // copy p to the first row
        for (j = 0; j < A->cols; j++)
        {
          set_d(&A->entry[i][j], &CD->codim[codim_index].p_d->coord[j]);
        }
      }
      else if (i <= codim)
      { // set entry to be \delta_ij
        for (j = 0; j < A->cols; j++)
        {
          set_zero_d(&A->entry[i][j]);
        }
        A->entry[i][i].r = -1;
      }
      else
      { // copy over B
        for (j = 0; j < A->cols; j++)
        {
          set_d(&A->entry[i][j], &CD->codim[codim_index].B_d->entry[i - codim - 1][j]);
        }
      }       

    // setup b to be [1,0..0]
    change_size_vec_d(b, CD->new_variables);
    b->size = CD->new_variables;
    for (i = 0; i < b->size; i++)
      if (i > 0)
      {
        set_zero_d(&b->coord[i]);
      }
      else
      {
        set_one_d(&b->coord[i]);
      }

    // adjust A & b if using hom variables
    if (CD->PPD.num_hom_var_gp)
    { // first row becomes p - h
      for (j = 0; j < A->cols; j++)
      {
        sub_d(&A->entry[0][j], &A->entry[0][j], &CD->codim[codim_index].H_d->coord[j]);
      }
      // first entry of b is homVarConst
      set_d(&b->coord[0], CD->codim[codim_index].homVarConst_d);
    }
      
    // loop through to find the start points    
    for (i = 0; i < num_paths; i++)
    { // setup A for this start point
      for (j = 0; j < codim; j++)
      {
        set_double_d(&A->entry[j+1][0], cos(2 * M_PI * curr_degs[j] / CD->new_degrees[j]), sin(2 * M_PI * curr_degs[j] / CD->new_degrees[j]));
      }

      // solve for the start point
      if (matrixSolve_d(CD->codim[codim_index].startPts_d[i], A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      // update curr_degs
      for (j = 0; j < codim; j++)
        if (curr_degs[j] + 1 == CD->new_degrees[j])
        { // set back to 0
          curr_degs[j] = 0;
        }
        else
        { // increment and then exit loop
          curr_degs[j]++;
          j = codim;
        }
    }
  }

  clear_mat_d(A);
  clear_vec_d(b);

  return;
}

void setupCodimStartPts_mp(codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for the start points for 'codim_index'           *
\***************************************************************/
{
  int i, j, num_paths = CD->codim[codim_index].num_paths, codim = CD->codim[codim_index].codim;
  int *curr_degs = (int *)bmalloc(codim * sizeof(int));
  mpf_t tempMPF, two_pi;
  mat_mp A;
  vec_mp b;

  mpf_init(tempMPF); mpf_init(two_pi);
  init_mat_mp(A, 0, 0);
  init_vec_mp(b, 0);

  // find 2*PI
  mpfr_const_pi(two_pi, __gmp_default_rounding_mode);
  mpf_add(two_pi, two_pi, two_pi);
  mpfr_free_cache(); // free the cache needed to create PI

  // the start system has the form [x_1^d_1 - x_0^d_0, .., x_codim^d_codim - x_0^d_codim]
  // with extrinsic slices: B*x = 0 and patch p*x = 1
  // with intrinsic slices: x = p + B*y, where y is tracking variables

  // initialize curr_degs
  for (i = 0; i < codim; i++)
    curr_degs[i] = 0;

  if (CD->codim[codim_index].useIntrinsicSlice)
  { // setup the start points for intrinsic tracking

    // setup A to originally be [B**,-1] ** - of the correct size (the extra column is for the homogenenous multiplier since we have variables y_1, .. y_codim, lambda)
    change_size_mat_mp(A, codim + 1, codim + 1);
    A->rows = A->cols = codim + 1;
    for (i = 0; i < A->rows; i++)
    { // copy over B to the first 'codim' columns
      for (j = 0; j < codim; j++)
      {
        set_mp(&A->entry[i][j], &CD->codim[codim_index].B_mp->entry[i][j]);
      }
      // set last column to [-1,..-1]
      mpf_set_si(A->entry[i][codim].r, -1);
      mpf_set_ui(A->entry[i][codim].i, 0);
    }

    // setup b to be [-p_0, .. -p_codim]
    change_size_vec_mp(b, codim + 1);
    b->size = codim + 1;
    for (i = 0; i < b->size; i++)
    {
      neg_mp(&b->coord[i], &CD->codim[codim_index].p_mp->coord[i]);
    }

    // loop through to find the start points
    for (i = 0; i < num_paths; i++)
    { // setup A for this start point
      for (j = 0; j < codim; j++)
      { // find the angle = 2*PI*curr_degs/new_degrees
        mpf_mul_ui(tempMPF, two_pi, curr_degs[j]);
        mpf_div_ui(tempMPF, tempMPF, CD->new_degrees[j]);

        // find -sin & -cos of angle
        mpfr_sin_cos(A->entry[j+1][codim].i, A->entry[j+1][codim].r, tempMPF, __gmp_default_rounding_mode);
        neg_mp(&A->entry[j+1][codim], &A->entry[j+1][codim]);
      }

      // solve for the start point
      if (matrixSolve_mp(CD->codim[codim_index].startPts_mp[i], A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      // remove the bottom coordinate (the homogeneous multiplier)
      CD->codim[codim_index].startPts_mp[i]->size -= 1;

      // update curr_degs
      for (j = 0; j < codim; j++)
        if (curr_degs[j] + 1 == CD->new_degrees[j])
        { // set back to 0
          curr_degs[j] = 0;
        }
        else
        { // increment and then exit loop
          curr_degs[j]++;
          j = codim;
        }
    }
  }
  else
  { // setup the start point for extrinsic tracking

    // setup A to originally be [p][0, -I, 0][B]
    change_size_mat_mp(A, CD->new_variables, CD->new_variables);
    A->rows = A->cols = CD->new_variables;
    for (i = 0; i < A->rows; i++)
      if (i < 1)
      { // copy p to the first row
        for (j = 0; j < A->cols; j++)
        {
          set_mp(&A->entry[i][j], &CD->codim[codim_index].p_mp->coord[j]);
        }
      }
      else if (i <= codim)
      { // set entry to be -\delta_ij
        for (j = 0; j < A->cols; j++)
        {
          set_zero_mp(&A->entry[i][j]);
        }
        mpf_set_si(A->entry[i][i].r, -1);
      }
      else
      { // copy over B
        for (j = 0; j < A->cols; j++)
        {
          set_mp(&A->entry[i][j], &CD->codim[codim_index].B_mp->entry[i - codim - 1][j]);
        }
      }

    // setup b to be [1,0..0]
    change_size_vec_mp(b, CD->new_variables);
    b->size = CD->new_variables;
    for (i = 0; i < b->size; i++)
      if (i > 0)
      {
        set_zero_mp(&b->coord[i]);
      }
      else
      {
        set_one_mp(&b->coord[i]);
      }

    // loop through to find the start points
    for (i = 0; i < num_paths; i++)
    { // setup A for this start point
      for (j = 0; j < codim; j++)
      { // find the angle = 2*PI*curr_degs/new_degrees
        mpf_mul_ui(tempMPF, two_pi, curr_degs[j]);
        mpf_div_ui(tempMPF, tempMPF, CD->new_degrees[j]);

        // find sin & cos of angle
        mpfr_sin_cos(A->entry[j+1][0].i, A->entry[j+1][0].r, tempMPF, __gmp_default_rounding_mode);
      }

      // solve for the start point
      if (matrixSolve_mp(CD->codim[codim_index].startPts_mp[i], A, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      // update curr_degs
      for (j = 0; j < codim; j++)
        if (curr_degs[j] + 1 == CD->new_degrees[j])
        { // set back to 0
          curr_degs[j] = 0;
        }
        else
        { // increment and then exit loop
          curr_degs[j]++;
          j = codim;
        }
    }
  }

  // clear MP
  mpf_clear(tempMPF); mpf_clear(two_pi);
  clear_mat_mp(A);
  clear_vec_mp(b);

  return;
}

void dimbydim_clear(codim_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear CD - positive dimensional tracking - dim-by-dim  *
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

  // clear the codimensions
  for (i = 0; i < CD->num_codim; i++)
  { // clear the ith codimData
    clearCodimData(&CD->codim[i], MPType);
  }

  // clear codim
  free(CD->codim);

  return;
}

void dimbydim_clear_start_points(codim_t *CD, int codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears the start points from 'codim_index'             *
\***************************************************************/
{
  int i, num_paths = CD->codim[codim_index].num_paths;

  if (MPType == 0 || MPType == 2)
  { // free the start points
    for (i = num_paths - 1; i >= 0; i--)
      clear_point_d(CD->codim[codim_index].startPts_d[i]);
    free(CD->codim[codim_index].startPts_d);
  }
  else if (MPType == 1)
  { // clear the start points
    for (i = num_paths - 1; i >= 0; i--)
      clear_point_mp(CD->codim[codim_index].startPts_mp[i]);
    free(CD->codim[codim_index].startPts_mp);
  }

  // NULL out _d & _mp
  CD->codim[codim_index].startPts_d = NULL;
  CD->codim[codim_index].startPts_mp = NULL;

  return;
}

void clearCodimData(codimData_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear for the codim                                    *
\***************************************************************/
{
  int i, num_paths = CD->num_paths;

  // clear endPt_types
  free(CD->endPt_types);

  // clear H
  clear_vec(CD->H_d, CD->H_mp, CD->H_rat, MPType);

  // clear homVarConst
  clear_d_mp_rat(CD->homVarConst_d, CD->homVarConst_mp, CD->homVarConst_rat, MPType);

  // clear A
  clear_mat(CD->A_d, CD->A_mp, CD->A_rat, MPType);

  // clear B
  clear_mat(CD->B_d, CD->B_mp, CD->B_rat, MPType);

  // clear p
  clear_vec(CD->p_d, CD->p_mp, CD->p_rat, MPType);

  if (MPType == 0)
  { // clear _d structures
    for (i = num_paths - 1; i >= 0; i--)
    { // clear startPts_d, if needed 
      if (CD->startPts_d != NULL)
        clear_point_d(CD->startPts_d[i]);
      // clear endPts_d, if needed 
      if (CD->endPts_d != NULL)
        clear_endpoint_data_d(&CD->endPts_d[i]);
    }
    if (CD->startPts_d != NULL)
      free(CD->startPts_d);
    if (CD->endPts_d != NULL)
      free(CD->endPts_d);
  }
  else if (MPType == 1)
  { // clear _mp structures
    for (i = num_paths - 1; i >= 0; i--)
    { // clear startPts_mp, if needed
      if (CD->startPts_mp != NULL)
        clear_point_mp(CD->startPts_mp[i]);
      // clear endPts_mp, if needed
      if (CD->endPts_mp != NULL)
        clear_endpoint_data_mp(&CD->endPts_mp[i]);
    }
    if (CD->startPts_mp != NULL)
      free(CD->startPts_mp);
    if (CD->endPts_mp != NULL)
      free(CD->endPts_mp);
  }
  else
  { // clear _d, _mp & _rat
    for (i = num_paths - 1; i >= 0; i--)
    { // clear startPts_d, if needed
      if (CD->startPts_d != NULL)
        clear_point_d(CD->startPts_d[i]);
      // clear endPts_amp, if needed
      if (CD->endPts_amp != NULL)
        clear_endpoint_data_amp(&CD->endPts_amp[i]);
    }
    if (CD->startPts_d != NULL)
      free(CD->startPts_d);
    if (CD->endPts_amp != NULL)
      free(CD->endPts_amp);
  }

  // NULL out the pointers
  CD->startPts_d = NULL;
  CD->startPts_mp = NULL;
  CD->endPts_d = NULL;
  CD->endPts_mp = NULL;
  CD->endPts_amp = NULL;

  return;
}

//// OPENMP COPY FUNCTIONS /////
void setup_omp_codimData(codim_t *CD, codim_t *CD_in, int codim_index, int MPType);
void clear_omp_codimData(codim_t *CD, int codim_index, int MPType);

void setup_omp_codim_t(codim_t *CD, codim_t *CD_in, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copies CD_in to CD for use in OpenMP tracking          *
\***************************************************************/
{
  int i;

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
    mpq_init(CD->gamma_rat[0]); mpq_init(CD->gamma_rat[1]);
    mpq_set(CD->gamma_rat[0], CD_in->gamma_rat[0]); mpq_set(CD->gamma_rat[1], CD_in->gamma_rat[1]);
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

  // allocate codim_t
  CD->codim = (codimData_t *)bmalloc(CD->num_codim * sizeof(codimData_t));

  // setup the codimensions
  for (i = 0; i < CD->num_codim; i++)
    setup_omp_codimData(CD, CD_in, i, MPType);

  return;
}

void setup_omp_codimData(codim_t *CD, codim_t *CD_in, int codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  CD->codim[codim_index].codim = CD_in->codim[codim_index].codim;
  CD->codim[codim_index].num_paths = CD_in->codim[codim_index].num_paths;
  CD->codim[codim_index].num_superset = CD_in->codim[codim_index].num_superset;
  CD->codim[codim_index].num_nonsing = CD_in->codim[codim_index].num_nonsing;
  CD->codim[codim_index].num_sing = CD_in->codim[codim_index].num_sing;
  CD->codim[codim_index].num_nonsolns = CD_in->codim[codim_index].num_nonsolns;
  CD->codim[codim_index].num_inf = CD_in->codim[codim_index].num_inf;
  CD->codim[codim_index].num_bad = CD_in->codim[codim_index].num_bad;
  CD->codim[codim_index].useIntrinsicSlice = CD_in->codim[codim_index].useIntrinsicSlice;

  // point to W
  CD->codim[codim_index].W = CD_in->codim[codim_index].W;

  // point to endPt_types
  CD->codim[codim_index].endPt_types = CD_in->codim[codim_index].endPt_types;

  // point to startPts
  CD->codim[codim_index].startPts_d = CD_in->codim[codim_index].startPts_d;
  CD->codim[codim_index].startPts_mp = CD_in->codim[codim_index].startPts_mp;

  // point to endPts
  CD->codim[codim_index].endPts_d = CD_in->codim[codim_index].endPts_d;
  CD->codim[codim_index].endPts_mp = CD_in->codim[codim_index].endPts_mp;
  CD->codim[codim_index].endPts_amp = CD_in->codim[codim_index].endPts_amp;

  if (MPType == 0)
  { // setup _d structures

    // copy H_d
    init_vec_d(CD->codim[codim_index].H_d, CD_in->codim[codim_index].H_d->size);
    vec_cp_d(CD->codim[codim_index].H_d, CD_in->codim[codim_index].H_d);

    // copy homVarConst_d
    set_d(CD->codim[codim_index].homVarConst_d, CD_in->codim[codim_index].homVarConst_d);

    // copy A_d
    init_mat_d(CD->codim[codim_index].A_d, CD_in->codim[codim_index].A_d->rows, CD_in->codim[codim_index].A_d->cols);
    mat_cp_d(CD->codim[codim_index].A_d, CD_in->codim[codim_index].A_d);

    // copy B_d
    init_mat_d(CD->codim[codim_index].B_d, CD_in->codim[codim_index].B_d->rows, CD_in->codim[codim_index].B_d->cols);
    mat_cp_d(CD->codim[codim_index].B_d, CD_in->codim[codim_index].B_d);

    // copy p_d
    init_vec_d(CD->codim[codim_index].p_d, CD_in->codim[codim_index].p_d->size);
    vec_cp_d(CD->codim[codim_index].p_d, CD_in->codim[codim_index].p_d);
  }
  else if (MPType == 1)
  { // setup _mp structures

    // copy H_mp
    init_vec_mp2(CD->codim[codim_index].H_mp, CD_in->codim[codim_index].H_mp->size, CD->curr_precision);
    vec_cp_mp(CD->codim[codim_index].H_mp, CD_in->codim[codim_index].H_mp);

    // copy homVarConst_mp
    init_mp2(CD->codim[codim_index].homVarConst_mp, CD->curr_precision);
    set_mp(CD->codim[codim_index].homVarConst_mp, CD_in->codim[codim_index].homVarConst_mp);

    // copy A_mp
    init_mat_mp2(CD->codim[codim_index].A_mp, CD_in->codim[codim_index].A_mp->rows, CD_in->codim[codim_index].A_mp->cols, CD->curr_precision);
    mat_cp_mp(CD->codim[codim_index].A_mp, CD_in->codim[codim_index].A_mp);

    // copy B_mp
    init_mat_mp2(CD->codim[codim_index].B_mp, CD_in->codim[codim_index].B_mp->rows, CD_in->codim[codim_index].B_mp->cols, CD->curr_precision);
    mat_cp_mp(CD->codim[codim_index].B_mp, CD_in->codim[codim_index].B_mp);

    // copy p_mp
    init_vec_mp2(CD->codim[codim_index].p_mp, CD_in->codim[codim_index].p_mp->size, CD->curr_precision);
    vec_cp_mp(CD->codim[codim_index].p_mp, CD_in->codim[codim_index].p_mp);
  }
  else if (MPType == 2)
  { // setup _d, _mp & _rat

    // copy H
    init_vec_d(CD->codim[codim_index].H_d, CD_in->codim[codim_index].H_d->size);
    init_vec_mp2(CD->codim[codim_index].H_mp, CD_in->codim[codim_index].H_mp->size, CD->curr_precision);
    vec_cp_d(CD->codim[codim_index].H_d, CD_in->codim[codim_index].H_d);
    vec_cp_mp(CD->codim[codim_index].H_mp, CD_in->codim[codim_index].H_mp);
    CD->codim[codim_index].H_rat = CD_in->codim[codim_index].H_rat;

    // copy homVarConst
    init_mp2(CD->codim[codim_index].homVarConst_mp, CD->curr_precision);
    init_rat(CD->codim[codim_index].homVarConst_rat);
    set_d(CD->codim[codim_index].homVarConst_d, CD_in->codim[codim_index].homVarConst_d);
    set_mp(CD->codim[codim_index].homVarConst_mp, CD_in->codim[codim_index].homVarConst_mp);
    set_rat(CD->codim[codim_index].homVarConst_rat, CD_in->codim[codim_index].homVarConst_rat);

    // copy A
    init_mat_d(CD->codim[codim_index].A_d, CD_in->codim[codim_index].A_d->rows, CD_in->codim[codim_index].A_d->cols);
    init_mat_mp2(CD->codim[codim_index].A_mp, CD_in->codim[codim_index].A_mp->rows, CD_in->codim[codim_index].A_mp->cols, CD->curr_precision);
    mat_cp_d(CD->codim[codim_index].A_d, CD_in->codim[codim_index].A_d);
    mat_cp_mp(CD->codim[codim_index].A_mp, CD_in->codim[codim_index].A_mp);
    CD->codim[codim_index].A_rat = CD_in->codim[codim_index].A_rat;

    // copy B
    init_mat_d(CD->codim[codim_index].B_d, CD_in->codim[codim_index].B_d->rows, CD_in->codim[codim_index].B_d->cols);
    init_mat_mp2(CD->codim[codim_index].B_mp, CD_in->codim[codim_index].B_mp->rows, CD_in->codim[codim_index].B_mp->cols, CD->curr_precision);
    mat_cp_d(CD->codim[codim_index].B_d, CD_in->codim[codim_index].B_d);
    mat_cp_mp(CD->codim[codim_index].B_mp, CD_in->codim[codim_index].B_mp);
    CD->codim[codim_index].B_rat = CD_in->codim[codim_index].B_rat;

    // copy p
    init_vec_d(CD->codim[codim_index].p_d, CD_in->codim[codim_index].p_d->size);
    init_vec_mp2(CD->codim[codim_index].p_mp, CD_in->codim[codim_index].p_mp->size, CD->curr_precision);
    vec_cp_d(CD->codim[codim_index].p_d, CD_in->codim[codim_index].p_d);
    vec_cp_mp(CD->codim[codim_index].p_mp, CD_in->codim[codim_index].p_mp);
    CD->codim[codim_index].p_rat = CD_in->codim[codim_index].p_rat;
  }

  return;
}

void clear_omp_codim_t(codim_t *CD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears CD where is was used in OpenMP tracking         *
\***************************************************************/
{
  int i;

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

  // clear the codimensions
  for (i = 0; i < CD->num_codim; i++)
    clear_omp_codimData(CD, i, MPType);

  // free codim
  free(CD->codim);

  return;
}

void clear_omp_codimData(codim_t *CD, int codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  // clear W
  CD->codim[codim_index].W = NULL;

  // clear endPt_types
  CD->codim[codim_index].endPt_types = NULL;

  // clear startPts
  CD->codim[codim_index].startPts_d = NULL;
  CD->codim[codim_index].startPts_mp = NULL;

  // clear endPts
  CD->codim[codim_index].endPts_d = NULL;
  CD->codim[codim_index].endPts_mp = NULL;
  CD->codim[codim_index].endPts_amp = NULL;

  if (MPType == 0)
  { // clear _d structures

    // clear H_d
    clear_vec_d(CD->codim[codim_index].H_d);

    // clear homVarConst_d
    clear_d(CD->codim[codim_index].homVarConst_d);

    // clear A_d
    clear_mat_d(CD->codim[codim_index].A_d);

    // clear B_d
    clear_mat_d(CD->codim[codim_index].B_d);

    // clear p_d
    clear_vec_d(CD->codim[codim_index].p_d);
  }
  else if (MPType == 1)
  { // clear _mp structures

    // clear H_mp
    clear_vec_mp(CD->codim[codim_index].H_mp);

    // clear homVarConst_mp
    clear_mp(CD->codim[codim_index].homVarConst_mp);

    // clear A_mp
    clear_mat_mp(CD->codim[codim_index].A_mp);

    // clear B_mp
    clear_mat_mp(CD->codim[codim_index].B_mp);

    // clear p_mp
    clear_vec_mp(CD->codim[codim_index].p_mp);
  }
  else if (MPType == 2)
  { // clear _d, _mp & _rat structures

    // clear H
    clear_vec_d(CD->codim[codim_index].H_d);
    clear_vec_mp(CD->codim[codim_index].H_mp);
    CD->codim[codim_index].H_rat = NULL;

    // clear homVarConst
    clear_d(CD->codim[codim_index].homVarConst_d);
    clear_mp(CD->codim[codim_index].homVarConst_mp);
    clear_rat(CD->codim[codim_index].homVarConst_rat);

    // clear A
    clear_mat_d(CD->codim[codim_index].A_d);
    clear_mat_mp(CD->codim[codim_index].A_mp);
    CD->codim[codim_index].A_rat = NULL;

    // clear B
    clear_mat_d(CD->codim[codim_index].B_d);
    clear_mat_mp(CD->codim[codim_index].B_mp);
    CD->codim[codim_index].B_rat = NULL;

    // clear p
    clear_vec_d(CD->codim[codim_index].p_d);
    clear_vec_mp(CD->codim[codim_index].p_mp);
    CD->codim[codim_index].p_rat = NULL;
  }

  return;
}

///// CHANGE PRECISION /////

int change_dimbydim_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for dim-by-dim                        *
\***************************************************************/
{
  int i, j, k;

  // cast ED as CD
  codim_t *CD = (codim_t *)ED;

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

    // change the precision for all of A_mp, B_mp & p_mp in codim
    for (k = 0; k < CD->num_codim; k++)
    { // change the precision for H
      for (i = 0; i < CD->codim[k].H_mp->size; i++)
      {
        setprec_mp(&CD->codim[k].H_mp->coord[i], prec);
        mpf_set_q(CD->codim[k].H_mp->coord[i].r, CD->codim[k].H_rat[i][0]);
        mpf_set_q(CD->codim[k].H_mp->coord[i].i, CD->codim[k].H_rat[i][1]);
      }

      // change the precision for homVarConst
      setprec_mp(CD->codim[k].homVarConst_mp, prec);
      mpf_set_q(CD->codim[k].homVarConst_mp->r, CD->codim[k].homVarConst_rat[0]);
      mpf_set_q(CD->codim[k].homVarConst_mp->i, CD->codim[k].homVarConst_rat[1]);

      // change the precision for A
      for (i = 0; i < CD->codim[k].A_mp->rows; i++)
        for (j = 0; j < CD->codim[k].A_mp->cols; j++)
        {
          setprec_mp(&CD->codim[k].A_mp->entry[i][j], prec);
          mpf_set_q(CD->codim[k].A_mp->entry[i][j].r, CD->codim[k].A_rat[i][j][0]);
          mpf_set_q(CD->codim[k].A_mp->entry[i][j].i, CD->codim[k].A_rat[i][j][1]);
        }

      // change the precision for B
      for (i = 0; i < CD->codim[k].B_mp->rows; i++)
        for (j = 0; j < CD->codim[k].B_mp->cols; j++)
        {
          setprec_mp(&CD->codim[k].B_mp->entry[i][j], prec);
          mpf_set_q(CD->codim[k].B_mp->entry[i][j].r, CD->codim[k].B_rat[i][j][0]);
          mpf_set_q(CD->codim[k].B_mp->entry[i][j].i, CD->codim[k].B_rat[i][j][1]);
        }
 
      // change the precision for p
      for (i = 0; i < CD->codim[k].p_mp->size; i++)
      {
        setprec_mp(&CD->codim[k].p_mp->coord[i], prec);
        mpf_set_q(CD->codim[k].p_mp->coord[i].r, CD->codim[k].p_rat[i][0]);
        mpf_set_q(CD->codim[k].p_mp->coord[i].i, CD->codim[k].p_rat[i][1]);
      }
    }
  }

  return 0;
}

