// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"
#include "eqbyeq.h"

// provides the functions to setup for using the equation-by-equation method of Sommese, Verschelde and Wampler

void setupEqbyEqRandom_d(eqData_t *EqD, tracker_config_t *T, char *degreeFile, preproc_data *PPD, square_system_eval_data_d *SqD, int *new_degrees, int curr_prec);
void setupEqbyEqRandom_mp(eqData_t *EqD, tracker_config_t *T, char *degreeFile, preproc_data *PPD, int *new_degrees);
void setupEqbyEqStages_d(eqData_t *EqD, tracker_config_t *T, char *depthFile, int intrinsicCutoff);
void setupEqbyEqStages_mp(eqData_t *EqD, tracker_config_t *T, char *depthFile, int intrinsicCutoff);
void setupEqbyEqWitnessData_d(eqData_t *EqD, patch_eval_data_d *PD_d);
void setupEqbyEqWitnessData_mp(eqData_t *EqD, patch_eval_data_mp *PD_mp);
void setupEqbyEqWitnessData_amp(eqData_t *EqD, patch_eval_data_mp *PD_mp, int max_prec);

void changeWitnessDataPrec(eqWitnessData_mp *EqD, int prec);
void changeStageDataPrec(eqStageData_mp *EqD, int prec);

int eqbyeq_setup_d(FILE **OUT, char *outName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, char *pointsIN, char *pointsOUT, char *depthFile, double intrinsicCutoffMultiplier)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of original variables                   *
* NOTES: setup for zero dimensional tracking using eq-by-eq     *
\***************************************************************/
{
  int num_orig_vars, intrinsicCutoff;
  FILE *midOUT = NULL;

  // allocate space
  ED->EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));

  // setup the function
  num_orig_vars = zero_dim_basic_setup_d(OUT, outName, &midOUT, "midpath_data", T, ED, dummyProg, &ED->EqD->startSub, &ED->EqD->endSub, &ED->EqD->startFunc, &ED->EqD->endFunc, &ED->EqD->startJvsub, &ED->EqD->endJvsub, &ED->EqD->startJv, &ED->EqD->endJv, &ED->EqD->subFuncsBelow, eval_d, eval_mp, preprocFile, degreeFile, 0, pointsIN, pointsOUT);

  // setup the point if using AMP
  if (T->MPType == 2)
    ED->BED_mp->EqD = ED->EqD;

  // midOUT does not need setup
  fclose(midOUT);

  // verify that we are using only 1 homogenous variable group
  if (ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: The equation-by-equation method of Sommese, Verschelde and Wampler is only implemented in Bertini for systems with\nonly one variable group.");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup whether there are changes & the number of subfunctions
  ED->EqD->noChanges = ED->squareSystem.noChanges;
  ED->EqD->numSubFuncs = ED->squareSystem.Prog->numSubfuncs;

  // setup the random numbers needed for eq-by-eq
  setupEqbyEqRandom_d(ED->EqD, T, degreeFile, &ED->preProcData, &ED->squareSystem, ED->squareSystem.new_degrees, 64);

  intrinsicCutoff = floor(intrinsicCutoffMultiplier * num_orig_vars);

  // allocate space for the witness data and the stage data
  setupEqbyEqStages_d(ED->EqD, T, depthFile, intrinsicCutoff);

  // setup the witness data
  setupEqbyEqWitnessData_d(ED->EqD, &ED->patch);
  if (T->MPType == 2)
  { // setup for AMP as well
    setupEqbyEqWitnessData_amp(ED->EqD, &ED->BED_mp->patch, T->AMP_max_prec);
  }

  return num_orig_vars;
}

void setupEqbyEqRandom_d(eqData_t *EqD, tracker_config_t *T, char *degreeFile, preproc_data *PPD, square_system_eval_data_d *SqD, int *new_degrees, int curr_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the random numbers for eq-by-eq                  *
\***************************************************************/
// ASSUME 1-hom!
{
  int i, j, k, beg, digits = 16, num_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  FILE *degIN = fopen(degreeFile, "r"); // open the file to read in the degrees
  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  EqD->increase_witness_prec = 1; // need to increase precision on witness data since this information is still needed
  EqD->num_funcs = T->numVars - num_var_gps; // calculate the true number of functions (we can ignore the projective transformations)
  EqD->num_var_gps = num_var_gps; // set the number of variable groups
  EqD->num_vars = T->numVars;     // set the number of variables

  // setup the degrees and allocate space for the random numbers
  EqD->degrees = (int **)bmalloc(EqD->num_funcs * sizeof(int *));
  EqD->coeff_d = (comp_d **)bmalloc(EqD->num_funcs * sizeof(comp_d *));

  // read in the degrees and allocate the space for the random numbers - the patches will be handled by the function evaluator and not in eq-by-eq
  for (i = 0; i < EqD->num_funcs; i++)
  { // setup degrees
    EqD->degrees[i] = (int *)bmalloc(num_var_gps * sizeof(int));

    if (num_var_gps == 1)
    { // the functions could have been permuted so we use new_degrees
      EqD->degrees[i][0] = new_degrees[i];
    }
    else
    { // need to read in all of the m-hom degrees from degIN
      for (j = 0; j < num_var_gps; j++)
        fscanf(degIN, "%d\n", &EqD->degrees[i][j]);
      fscanf(degIN, "\n"); // extra "new line" character
    }

    // allocate based on the total number of variables & initialize to zero
    EqD->coeff_d[i] = (comp_d *)bmalloc(T->numVars * sizeof(comp_d));
    for (j = 0; j < T->numVars; j++)
      set_zero_d(EqD->coeff_d[i][j]);
  }
  // close degreeFile
  fclose(degIN);

  // now that we have the space, initialize it and generate the random numbers

  // generate the random gamma
  get_comp_rand_d(EqD->gamma_d);

  // now do random numbers for the appropriate coeff
  for (i = 0; i < EqD->num_funcs; i++)
  { // need to make the linears truly linears in only the variable groups - zero for the other coeff
    if (PPD->type[0])
    { // this is a regular variable group - thus it has a homogenous coordinate that was created
      beg = 0;

      get_comp_rand_d(EqD->coeff_d[i][beg]);  // get random for hom coord for this variable group
    }
    beg = PPD->num_var_gp;

    for (k = 0; k < PPD->size[0]; k++)
      get_comp_rand_d(EqD->coeff_d[i][k + beg]); // get random for variables that are in this variable group
  }

  if (T->MPType == 2)
  { // convert all to mp and rationals
    T->Precision = curr_prec;
    initMP(curr_prec);
    EqD->curr_precision = curr_prec;

    // gamma
    get_comp_rand_rat(EqD->gamma_d, EqD->gamma_mp, EqD->gamma_rat, curr_prec, T->AMP_max_prec, 1, 1);

    // coeff
    EqD->coeff_mp = (comp_mp **)bmalloc(EqD->num_funcs * sizeof(comp_mp *));
    EqD->coeff_rat = (mpq_t ***)bmalloc(EqD->num_funcs * sizeof(mpq_t **));
    for (i = 0; i < EqD->num_funcs; i++)
    {
      EqD->coeff_mp[i] = (comp_mp *)bmalloc(T->numVars * sizeof(comp_mp));
      EqD->coeff_rat[i] = (mpq_t **)bmalloc(T->numVars * sizeof(mpq_t *));

      for (k = 0; k < T->numVars; k++)
      { // allocate for real & imaginary rational numbers
        EqD->coeff_rat[i][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
        // check to see if coeff_d is zero or not
        if (EqD->coeff_d[i][k]->r == 0 && EqD->coeff_d[i][k]->i == 0)
        { // just copy over 0 to _mp & _rat
          comp_d_to_mp_rat(EqD->coeff_mp[i][k], EqD->coeff_rat[i][k], EqD->coeff_d[i][k], digits, curr_prec, 1, 1);
        }
        else
        { // setup coeff_d, coeff_mp & coeff_rat to 'extended precision'
          get_comp_rand_rat(EqD->coeff_d[i][k], EqD->coeff_mp[i][k], EqD->coeff_rat[i][k], curr_prec, T->AMP_max_prec, 1, 1);
        }
      }
    }
  }

  return;
}

void setupEqbyEqStages_d(eqData_t *EqD, tracker_config_t *T, char *depthFile, int intrinsicCutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: allocate memory for the witness data and the stages    *
\***************************************************************/
{ // NOTE: If using AMP, all start points and end points are stored as doubles - for now
  int i, j, *depths = NULL;
  FILE *DEPTH = fopen(depthFile, "r");

  if (DEPTH == NULL)
  { // setup with each depth as 1
    EqD->num_subsystems = EqD->num_funcs;

    // allocate space
    EqD->witnessData_d = (eqWitnessData_d *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_d));
    EqD->stageData_d = (eqStageData_d *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_d));

    if (T->MPType == 2)
    { // allocate MP as well
      EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_mp));
      EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_mp));
    }

    // setup the basic info for the witness data
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      EqD->witnessData_d[i].startFunction = i;
      EqD->witnessData_d[i].depth = 1;

      if (T->MPType == 2)
      {
        EqD->witnessData_mp[i].startFunction = EqD->witnessData_d[i].startFunction;
        EqD->witnessData_mp[i].depth = EqD->witnessData_d[i].depth;
      }
    }
  }
  else
  { // read in the depths
    EqD->num_subsystems = 1;
    fscanf(DEPTH, "%d\n", &EqD->num_subsystems); // read in the number of subsystems that will be needed

    // error checking
    if (EqD->num_subsystems <= 0)
    {
      printf("ERROR: The number of subsystems must be > 0!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // setup depths
    depths = (int *)bmalloc(EqD->num_subsystems * sizeof(int));
    j = 0;
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      depths[i] = 1;
      fscanf(DEPTH, "%d\n", &depths[i]);
      if (depths[i] <= 0)
      {
        printf("ERROR: The depths must be > 0!\n");
        bexit(ERROR_CONFIGURATION);
      }
      j += depths[i]; // sum up the total depths
    }

    if (j != EqD->num_funcs) // the total depth is not correct
    {
      printf("ERROR: The number of functions (%d) is not equal to the total depth (%d)!\n", EqD->num_funcs, j);
      bexit(ERROR_CONFIGURATION);
    }

    // allocate space
    EqD->witnessData_d = (eqWitnessData_d *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_d));
    EqD->stageData_d = (eqStageData_d *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_d));

    if (T->MPType == 2)
    { // allocate MP as well
      EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_mp));
      EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_mp));
    }

    // setup the basic info for the witness data
    j = 0;
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      EqD->witnessData_d[i].startFunction = j;
      EqD->witnessData_d[i].depth = depths[i];

      if (T->MPType == 2)
      {
        EqD->witnessData_mp[i].startFunction = EqD->witnessData_d[i].startFunction;
        EqD->witnessData_mp[i].depth = EqD->witnessData_d[i].depth;
      }
      j += depths[i];
    }

    free(depths);
    fclose(DEPTH);
  }

  // setup the basic info for the stages
  for (i = 0; i < EqD->num_subsystems; i++)
  { // setup depth_x
    if (i == 0)
      EqD->stageData_d[i].depth_x = 0;
    else
      EqD->stageData_d[i].depth_x = EqD->stageData_d[i-1].depth_x + EqD->stageData_d[i-1].depth_y;
    // setup depth_y
    EqD->stageData_d[i].depth_y = EqD->witnessData_d[i].depth;

    // determine if we are using intrinsic slice
    if (i > 0 && EqD->stageData_d[i].depth_x + EqD->stageData_d[i].depth_y > intrinsicCutoff)
      EqD->stageData_d[i].useIntrinsicSlice = 0;
    else
      EqD->stageData_d[i].useIntrinsicSlice = 1;

    if (T->MPType == 2)
    {
      EqD->stageData_mp[i].depth_x = EqD->stageData_d[i].depth_x;
      EqD->stageData_mp[i].depth_y = EqD->stageData_d[i].depth_y;
      EqD->stageData_mp[i].useIntrinsicSlice = EqD->stageData_d[i].useIntrinsicSlice;
    }
  }

  return;
}

void setupEqbyEqWitnessData_d(eqData_t *EqD, patch_eval_data_d *PD_d)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the internal structures of witness data          *
\***************************************************************/
// ASSUME 1-hom!!!
{
  int i, j, k, r, start, *curr_deg = NULL;
  vec_d b;
  mat_d tempMat, A, Q, R, P;
  double tol_pivot = 1e-15, tol_sign = 1e-20, largeChange = 1e14;

  init_vec_d(b, 0);
  init_mat_d(tempMat, 0, 0); init_mat_d(A, 0, 0); init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

  for (i = 0; i < EqD->num_subsystems; i++)
  { // setup the ith witness data

    // initialize other values
    EqD->witnessData_d[i].num_sing = EqD->witnessData_d[i].num_nonsing = EqD->witnessData_d[i].num_inf = EqD->witnessData_d[i].num_bad = 0;

    // setup the number of paths
    EqD->witnessData_d[i].num_linears = 0; // the sum of the degrees
    EqD->witnessData_d[i].num_paths = 1; // the multiplication of the degrees
    start = EqD->witnessData_d[i].startFunction;
    for (j = 0; j < EqD->witnessData_d[i].depth; j++)
    {
      EqD->witnessData_d[i].num_linears += EqD->degrees[start + j][0];
      EqD->witnessData_d[i].num_paths *= EqD->degrees[start + j][0];
    }

    // setup the start system information - we will be using intrinsic slices to make it only 'depth' number of variables to track
    EqD->witnessData_d[i].startSystemCoeff = (comp_d **)bmalloc(EqD->witnessData_d[i].num_linears * sizeof(comp_d *));
    for (j = 0; j < EqD->witnessData_d[i].num_linears; j++)
    { // setup the jth linear - assuming here we are in the 1-hom case
      EqD->witnessData_d[i].startSystemCoeff[j] = (comp_d *)bmalloc(EqD->witnessData_d[i].depth * sizeof(comp_d));
      for (k = 0; k < EqD->witnessData_d[i].depth; k++)
      { // get a random number
        get_comp_rand_d(EqD->witnessData_d[i].startSystemCoeff[j][k]);
      }
    }

    // setup A
    k = EqD->num_funcs - EqD->witnessData_d[i].depth + PD_d->num_patches;
    change_size_mat_d(A, k, EqD->num_vars);
    A->rows = k;
    A->cols = EqD->num_vars; // == EqD->num_funcs + PD_d->num_patches
    // setup tempMat & b
    change_size_mat_d(tempMat, EqD->num_vars, EqD->num_vars);
    change_size_vec_d(b, EqD->num_vars);
    tempMat->rows = tempMat->cols = b->size = EqD->num_vars;

    // setup to find the intrinsic slice
    r = 0;
    for (j = 0; j < EqD->num_funcs; j++)
      if (j < EqD->witnessData_d[i].startFunction || j >= EqD->witnessData_d[i].startFunction + EqD->witnessData_d[i].depth)
      { // need to copy over the coefficients
        for (k = 0; k < EqD->num_vars; k++)
        {
          set_d(&A->entry[r][k], EqD->coeff_d[j][k]); // assume 1-hom
          set_d(&tempMat->entry[r][k], &A->entry[r][k]);
        }
        // set b[r] to 0
        set_zero_d(&b->coord[r]);
        r++;
      } 
    // copy over the patch coefficients
    for (j = 0; j < PD_d->num_patches; j++)
    { // set b[r] to 1
      set_one_d(&b->coord[r]);
      // copy patch
      for (k = 0; k < EqD->num_vars; k++)
      {
        set_d(&A->entry[r][k], &PD_d->patchCoeff->entry[j][k]);
        set_d(&tempMat->entry[r][k], &A->entry[r][k]);
      }
      r++;
    }
    // setup others as random
    for (r = r; r < EqD->num_vars; r++)
    { // set b[r] to random
      get_comp_rand_d(&b->coord[r]);
      for (k = 0; k < EqD->num_vars; k++)
        get_comp_rand_d(&tempMat->entry[r][k]);
    }

    // compute p
    init_vec_d(EqD->witnessData_d[i].p, 0);
    if (matrixSolve_d(EqD->witnessData_d[i].p, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_d(A, A); 
    QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    init_mat_d(EqD->witnessData_d[i].B, EqD->num_vars, EqD->witnessData_d[i].depth);
    EqD->witnessData_d[i].B->rows = EqD->num_vars;
    EqD->witnessData_d[i].B->cols = EqD->witnessData_d[i].depth;
    start = A->cols;

    for (j = 0; j < EqD->witnessData_d[i].B->rows; j++)
      for (k = 0; k < EqD->witnessData_d[i].B->cols; k++)
      {
        set_d(&EqD->witnessData_d[i].B->entry[j][k], &Q->entry[j][k+start]);
      }

    // setup other structures
    EqD->witnessData_d[i].startPts = (point_d *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(point_d));
    EqD->witnessData_d[i].endPts_in = (point_d *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(point_d));
    EqD->witnessData_d[i].endPts = (point_d *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(point_d));
    EqD->witnessData_d[i].finalTs = (comp_d *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(comp_d));
    EqD->witnessData_d[i].condition_nums = (double *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(double));
    EqD->witnessData_d[i].endPt_retVals = (int *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(int));
    EqD->witnessData_d[i].endPt_types = (int *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(int));
    EqD->witnessData_d[i].higherDim = (int *)bmalloc(EqD->witnessData_d[i].num_paths * sizeof(int));
    for (j = 0; j < EqD->witnessData_d[i].num_paths; j++)
    { // initialize the memory
      init_point_d(EqD->witnessData_d[i].startPts[j], EqD->witnessData_d[i].depth);
      init_point_d(EqD->witnessData_d[i].endPts_in[j], EqD->witnessData_d[i].depth);
      EqD->witnessData_d[i].startPts[j]->size = EqD->witnessData_d[i].depth;
      EqD->witnessData_d[i].endPts_in[j]->size = EqD->witnessData_d[i].depth;

      init_point_d(EqD->witnessData_d[i].endPts[j], EqD->num_vars);
      EqD->witnessData_d[i].endPts[j]->size = EqD->num_vars;

      set_zero_d(EqD->witnessData_d[i].finalTs[j]);
      EqD->witnessData_d[i].condition_nums[j] = EqD->witnessData_d[i].endPt_retVals[j] = EqD->witnessData_d[i].endPt_types[j] = EqD->witnessData_d[i].higherDim[j] = 0;
    }

    // setup the start points
    curr_deg = (int *)brealloc(curr_deg, EqD->witnessData_d[i].depth * sizeof(int));
    for (j = 0; j < EqD->witnessData_d[i].depth; j++)
      curr_deg[j] = 0;

    change_size_mat_d(tempMat, EqD->witnessData_d[i].depth, EqD->witnessData_d[i].depth);
    change_size_vec_d(b, EqD->witnessData_d[i].depth);
    tempMat->rows = tempMat->cols = b->size = EqD->witnessData_d[i].depth;
    // setup b
    for (k = 0; k < b->size; k++)
    {
      set_one_d(&b->coord[k]);
    }
    for (j = 0; j < EqD->witnessData_d[i].num_paths; j++)
    { // setup tempMat
      start = 0;
      for (k = 0; k < tempMat->rows; k++)
      { // copy over the entries
        for (r = 0; r < tempMat->cols; r++)
        {
          set_d(&tempMat->entry[k][r], EqD->witnessData_d[i].startSystemCoeff[start + curr_deg[k]][r]);
        }
        // update start
        start += EqD->degrees[EqD->witnessData_d[i].startFunction + k][0];
      }
      // solve for the start point
      r = matrixSolve_d(EqD->witnessData_d[i].startPts[j], tempMat, b);
      if (r)
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      } 

      // update curr_deg
      for (k = EqD->witnessData_d[i].depth - 1; k >= 0; k--)
      { // check to see if we are at the top
        if (curr_deg[k] == EqD->degrees[EqD->witnessData_d[i].startFunction + k][0] - 1)
        { // set to 0
          curr_deg[k] = 0;
        }
        else
        { // increment this location and exit loop
          curr_deg[k]++;
          k = -10;
        }
      }
    }  
  }    

  // clear memory
  clear_vec_d(b);
  clear_mat_d(tempMat); clear_mat_d(A); clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);

  free(curr_deg);

  return;
}

void setupEqbyEqFirstStage_d(eqData_t *EqD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first stage of eq-by-eq                      *
\***************************************************************/
{
  // do not worry about increasing the precision for the witness data structures
  EqD->increase_witness_prec = 0;

  // copy over the path information
  EqD->stageData_d[0].num_paths = EqD->witnessData_d[0].num_paths;
  EqD->stageData_d[0].num_sing = EqD->witnessData_d[0].num_sing;
  EqD->stageData_d[0].num_nonsing = EqD->witnessData_d[0].num_nonsing;
  EqD->stageData_d[0].num_inf = EqD->witnessData_d[0].num_inf;
  EqD->stageData_d[0].num_bad = EqD->witnessData_d[0].num_bad;

  // copy B and p
  EqD->stageData_d[0].useIntrinsicSlice = 1;
  init_mat_d(EqD->stageData_d[0].B, EqD->witnessData_d[0].B->rows, EqD->witnessData_d[0].B->cols);
  init_vec_d(EqD->stageData_d[0].p, EqD->witnessData_d[0].p->size);
  mat_cp_d(EqD->stageData_d[0].B, EqD->witnessData_d[0].B);
  vec_cp_d(EqD->stageData_d[0].p, EqD->witnessData_d[0].p);

  // point to the structures inside of witnessData[0]
  EqD->stageData_d[0].startPts = EqD->witnessData_d[0].startPts;
  EqD->stageData_d[0].endPts_in = EqD->witnessData_d[0].endPts_in;
  EqD->stageData_d[0].endPts = EqD->witnessData_d[0].endPts;
  EqD->stageData_d[0].finalTs = EqD->witnessData_d[0].finalTs;
  EqD->stageData_d[0].condition_nums = EqD->witnessData_d[0].condition_nums;
  EqD->stageData_d[0].endPt_retVals = EqD->witnessData_d[0].endPt_retVals;
  EqD->stageData_d[0].endPt_types = EqD->witnessData_d[0].endPt_types;
  EqD->stageData_d[0].higherDim = EqD->witnessData_d[0].higherDim;

  // clear out B1, p1, B0, p0
  init_mat_d(EqD->stageData_d[0].B1, 0, 0);
  init_vec_d(EqD->stageData_d[0].p1, 0);
  init_mat_d(EqD->stageData_d[0].B0, 0, 0);
  init_vec_d(EqD->stageData_d[0].p0, 0);

  if (MPType == 2)
  { // copy over B and p 
    init_mat_mp2(EqD->stageData_mp[0].B, EqD->witnessData_mp[0].B->rows, EqD->witnessData_mp[0].B->cols, EqD->curr_precision);
    mat_cp_mp(EqD->stageData_mp[0].B, EqD->witnessData_mp[0].B);
    init_vec_mp2(EqD->stageData_mp[0].p, EqD->witnessData_mp[0].p->size, EqD->curr_precision);
    vec_cp_mp(EqD->stageData_mp[0].p, EqD->witnessData_mp[0].p);

    // point to B_rat & p_rat
    EqD->stageData_mp[0].useIntrinsicSlice = 1;
    EqD->stageData_mp[0].B_rat = EqD->witnessData_mp[0].B_rat;
    EqD->stageData_mp[0].p_rat = EqD->witnessData_mp[0].p_rat;

    // NULL out B1_rat, p1_rat, B0_rat, p0_rat
    EqD->stageData_mp[0].p1_rat = EqD->stageData_mp[0].p0_rat = NULL;
    EqD->stageData_mp[0].B1_rat = EqD->stageData_mp[0].B0_rat = NULL;

    init_mat_mp(EqD->stageData_mp[0].B1, 0, 0);
    init_vec_mp(EqD->stageData_mp[0].p1, 0);
    init_mat_mp(EqD->stageData_mp[0].B0, 0, 0);
    init_vec_mp(EqD->stageData_mp[0].p0, 0);

    // set all other pointers to NULL
    EqD->stageData_mp[0].startPts = NULL;
    EqD->stageData_mp[0].endPts_in = NULL;
    EqD->stageData_mp[0].endPts = NULL;
    EqD->stageData_mp[0].finalTs = NULL;
    EqD->stageData_mp[0].condition_nums = NULL;
    EqD->stageData_mp[0].endPt_retVals = NULL;
    EqD->stageData_mp[0].endPt_types = NULL;
    EqD->stageData_mp[0].higherDim = NULL;
  }

  return;
}

void clearEqbyEqFirstWitnessData_d(eqData_t *EqD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the first witness data information               *
\***************************************************************/
{
  // set pointers to NULL since stage[0] uses this information
  EqD->witnessData_d[0].startPts = NULL;
  EqD->witnessData_d[0].endPts_in = NULL;
  EqD->witnessData_d[0].endPts = NULL;
  EqD->witnessData_d[0].finalTs = NULL;
  EqD->witnessData_d[0].condition_nums = NULL;
  EqD->witnessData_d[0].endPt_retVals = NULL;
  EqD->witnessData_d[0].endPt_types = NULL;
  EqD->witnessData_d[0].higherDim = NULL;

  // clear B and p
  clear_mat_d(EqD->witnessData_d[0].B);
  clear_vec_d(EqD->witnessData_d[0].p);

  if (MPType == 2)
  { // clear B and p and set all other pointers to NULL
    clear_mat_mp(EqD->witnessData_mp[0].B);
    clear_vec_mp(EqD->witnessData_mp[0].p);

    // stage[0] uses this information
    EqD->witnessData_mp[0].B_rat = NULL;
    EqD->witnessData_mp[0].p_rat = NULL;

    // not used
    EqD->witnessData_mp[0].startPts = NULL;
    EqD->witnessData_mp[0].endPts_in = NULL;
    EqD->witnessData_mp[0].endPts = NULL;
    EqD->witnessData_mp[0].finalTs = NULL;
    EqD->witnessData_mp[0].condition_nums = NULL;
    EqD->witnessData_mp[0].endPt_retVals = NULL;
    EqD->witnessData_mp[0].endPt_types = NULL;
    EqD->witnessData_mp[0].higherDim = NULL;
  }

  return;
}

void setupInt_d(basic_eval_data_d *ED, int stage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the things needed to do intrinisic tracking      *
* using only double precision                                   *
\***************************************************************/
{ 
  mat_d tempMat, A, Q, R, P;
  vec_d b, p_temp;
  double tol_pivot = 1e-15, tol_sign = 1e-20, largeChange = 1e14;
  int j, k, depth_x, depth_y, depth_sum, count;

  init_mat_d(tempMat, 0, 0); init_mat_d(A, 0, 0);
  init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);
  init_vec_d(b, 0); init_vec_d(p_temp, 0);

  // setup depth_x & depth_y
  depth_x = ED->EqD->stageData_d[stage].depth_x;
  depth_y = ED->EqD->stageData_d[stage].depth_y;
  depth_sum = depth_x + depth_y;

  // find B & p for the 'fixed linears'
  k = ED->EqD->num_funcs - depth_sum + ED->patch.num_patches;
  change_size_mat_d(A, k, ED->EqD->num_vars);
  A->rows = k;
  A->cols = ED->EqD->num_vars; // == ED->EqD->num_funcs + ED->patch.num_patches
  // setup tempMat & b
  change_size_mat_d(tempMat, ED->EqD->num_vars, ED->EqD->num_vars);
  change_size_vec_d(b, ED->EqD->num_vars);  
  tempMat->rows = tempMat->cols = b->size = ED->EqD->num_vars;

  // setup to find the intrinsic slice
  count = 0;
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients of the linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], ED->EqD->coeff_d[j][k]); // assume 1-hom
      set_d(&tempMat->entry[count][k], &A->entry[count][k]);
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_d(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
      set_d(&tempMat->entry[count][k], &A->entry[count][k]);
    }
    count++;
  }
  // setup rest as random
  for (count = count; count < ED->EqD->num_vars; count++)
  { // setup b[count] as random
    get_comp_rand_d(&b->coord[count]);
    for (k = 0; k < ED->EqD->num_vars; k++)
      get_comp_rand_d(&tempMat->entry[count][k]);
  }
 
  // setup p
  init_vec_d(ED->EqD->stageData_d[stage].p, tempMat->cols);
  if (matrixSolve_d(ED->EqD->stageData_d[stage].p, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }

  // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
  transpose_d(A, A);
  QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
  init_mat_d(ED->EqD->stageData_d[stage].B, ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B->rows = ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B->cols = depth_sum;
  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_d[stage].B->rows; j++)
    for (k = 0; k < ED->EqD->stageData_d[stage].B->cols; k++)
    {
      set_d(&ED->EqD->stageData_d[stage].B->entry[j][k], &Q->entry[j][k+count]);
    }

  // setup B1,p1,B0,p0
  k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
  change_size_mat_d(A, k, 2 * ED->EqD->num_vars);
  A->rows = k;
  A->cols = 2 * ED->EqD->num_vars; // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
  // setup tempMat & b
  change_size_mat_d(tempMat, 2 * ED->EqD->num_vars, 2 * ED->EqD->num_vars);
  change_size_vec_d(b, 2 * ED->EqD->num_vars);
  tempMat->rows = tempMat->cols = b->size = 2 * ED->EqD->num_vars; 

  // setup to find the intrinsic slice
  count = 0;
  // setup linears in the x variables
  for (j = depth_x; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients of the linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], ED->EqD->coeff_d[j][k]); // assume 1-hom
      set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for x variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_d(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
      set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
    }
    count++;
  }
  // setup linears in the y variables
  for (j = 0; j < depth_x; j++)
  { // copy over the coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    { 
      set_zero_d(&A->entry[count][k]);
      set_d(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_d[j][k]); // assume 1-hom
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_d(&A->entry[count][k]);
      set_d(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_d[j][k]); // assume 1-hom
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for y variables 
  for (j = 0; j < ED->patch.num_patches; j++) 
  { // set b[count] to 1 
    set_one_d(&b->coord[count]); 
    // copy patch 
    for (k = 0; k < ED->EqD->num_vars; k++) 
    { 
      set_zero_d(&A->entry[count][k]);
      set_d(&A->entry[count][k + ED->EqD->num_vars], &ED->patch.patchCoeff->entry[j][k]);
    } 
    count++; 
  }

  // setup tempMat to find p
  for (count = 0; count < 2 * ED->EqD->num_vars; count++)
    if (count < A->rows)
    { // copy A
      for (k = 0; k < A->cols; k++)
        set_d(&tempMat->entry[count][k], &A->entry[count][k]);
    }
    else
    { // setup rest as random
      get_comp_rand_d(&b->coord[count]);
      for (k = 0; k < A->cols; k++)
        get_comp_rand_d(&tempMat->entry[count][k]);
    }

  // setup p1
  init_vec_d(ED->EqD->stageData_d[stage].p1, 2 * ED->EqD->num_vars);
  if (matrixSolve_d(ED->EqD->stageData_d[stage].p1, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }

  // do QR decomposition on A^T and B1 is extra columns of Q
  transpose_d(A, A);
  QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
  init_mat_d(ED->EqD->stageData_d[stage].B1, 2 * ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B1->rows = 2 * ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B1->cols = depth_sum;

  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_d[stage].B1->rows; j++)
    for (k = 0; k < ED->EqD->stageData_d[stage].B1->cols; k++)
    {
      set_d(&ED->EqD->stageData_d[stage].B1->entry[j][k], &Q->entry[j][k+count]);
    }

  // now we do the same procedure to find p0 & B0
  k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
  change_size_mat_d(A, k, 2 * ED->EqD->num_vars);
  A->rows = k;
  A->cols = 2 * ED->EqD->num_vars; // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
  // setup tempMat & b
  change_size_mat_d(tempMat, k, k);
  change_size_vec_d(b, k);
  tempMat->rows = tempMat->cols = b->size = k;

  count = 0;
  // setup x_j - y_j for j = depth_x to depth_sum
  for (j = depth_x; j < depth_sum; j++)
  {
    for (k = 0; k < ED->EqD->num_vars; k++)
      if (j + 1 == k)
      {
        set_one_d(&A->entry[count][k]);
        neg_d(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
      }
      else
      {
        set_zero_d(&A->entry[count][k]);
        set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
      }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // setup linears in the x variables
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients of the linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], ED->EqD->coeff_d[j][k]); // assume 1-hom
      set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for x variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_d(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_d(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
      set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
    }
    count++;
  }
  // setup x_j - y_j for j = 0 to depth_x
  for (j = 0; j < depth_x; j++)
  {
    for (k = 0; k < ED->EqD->num_vars; k++)
      if (j + 1 == k)
      {
        set_one_d(&A->entry[count][k]);
        neg_d(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
      }
      else
      {
        set_zero_d(&A->entry[count][k]);
        set_zero_d(&A->entry[count][k + ED->EqD->num_vars]);
      }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // setup linears in the y variables
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_d(&A->entry[count][k]);
      set_d(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_d[j][k]); // assume 1-hom
    }
    // set b[count] to 0
    set_zero_d(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for y variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_d(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_d(&A->entry[count][k]);
      set_d(&A->entry[count][k + ED->EqD->num_vars], &ED->patch.patchCoeff->entry[j][k]);
    }
    count++;
  }

  // setup tempMat to find p using a randomization process
  make_matrix_random_d(R, A->cols, A->rows);
  mat_mul_d(tempMat, A, R);

  // setup p_temp
  if (matrixSolve_d(p_temp, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }
  // undo the randomization to find p0
  init_vec_d(ED->EqD->stageData_d[stage].p0, 0);
  mul_mat_vec_d(ED->EqD->stageData_d[stage].p0, R, p_temp);

  // do QR decomposition on A^T and B0 is extra columns of Q
  transpose_d(A, A);
  QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
  init_mat_d(ED->EqD->stageData_d[stage].B0, 2 * ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B0->rows = 2 * ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B0->cols = depth_sum;

  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_d[stage].B0->rows; j++)
    for (k = 0; k < ED->EqD->stageData_d[stage].B0->cols; k++)
    {
      set_d(&ED->EqD->stageData_d[stage].B0->entry[j][k], &Q->entry[j][k+count]);
    }

  // clear memory
  clear_mat_d(tempMat); clear_mat_d(A);
  clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
  clear_vec_d(b); clear_vec_d(p_temp);

  return;
}

void setupInt_amp(basic_eval_data_d *ED, int stage, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the things needed to do intrinisic tracking      *
* using adaptive multi precision                                *
\***************************************************************/
{
  int num_digits = prec_to_digits(max_prec);
  size_t size;
  char *str = NULL;
  vec_mp tempVec, b;
  mat_mp tempMat, A, Q, R, P;
  mpf_t tol_pivot, tol_sign, largeChange;

  int j, k, depth_x, depth_y, depth_sum, count;

  // setup depth_x & depth_y
  depth_x = ED->EqD->stageData_d[stage].depth_x;
  depth_y = ED->EqD->stageData_d[stage].depth_y;
  depth_sum = depth_x + depth_y;

  // initialize MP
  mpf_init2(tol_pivot, max_prec); mpf_init2(tol_sign, max_prec); mpf_init2(largeChange, max_prec);
  init_vec_mp2(tempVec, 0, max_prec); init_vec_mp2(b, 0, max_prec);
  init_mat_mp2(tempMat, 0, 0, max_prec); init_mat_mp2(A, 0, 0, max_prec); init_mat_mp2(Q, 0, 0, max_prec);
  init_mat_mp2(R, 0, 0, max_prec); init_mat_mp2(P, 0, 0, max_prec);

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

  // find B & p for the 'fixed linears'
  k = ED->EqD->num_funcs - depth_sum + ED->patch.num_patches;
  change_size_mat_mp(A, k, ED->EqD->num_vars);
  A->rows = k;
  A->cols = ED->EqD->num_vars; // == ED->EqD->num_funcs + num_patches
  // setup tempMat & b
  change_size_mat_mp(tempMat, k, k);
  change_size_vec_mp(b, k);
  tempMat->rows = tempMat->cols = b->size = k;

  // setup to find the intrinsic slice
  count = 0;
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->EqD->coeff_rat[j][k][1]);
      if (k < tempMat->cols)
      { // copy to tempMat
        set_mp(&tempMat->entry[count][k], &A->entry[count][k]);
      }
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients
  for (j = 0; j < ED->BED_mp->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_mp(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->BED_mp->patch.patchCoeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->BED_mp->patch.patchCoeff_rat[j][k][1]);
      if (k < tempMat->cols)
      { // copy to tempMat
        set_mp(&tempMat->entry[count][k], &A->entry[count][k]);
      }
    }
    count++;
  }

  // set the global prec to the maximum precision
  initMP(max_prec);

  // setup p - putting 0 in the extra positions
  // by doing it this way, we have a standard way of constructing p from coeff
  if (matrixSolve_mp(tempVec, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }
  // add 0's
  count = tempMat->cols + depth_sum;
  increase_size_vec_mp(tempVec, count);
  tempVec->size = count;
  for (j = tempMat->cols; j < count; j++)
  {
    set_zero_mp(&tempVec->coord[j]);
  }

  // copy tempVec to p
  init_vec_d(ED->EqD->stageData_d[stage].p, tempVec->size);
  init_vec_mp2(ED->EqD->stageData_mp[stage].p, tempVec->size, ED->EqD->curr_precision);
  init_vec_rat(ED->EqD->stageData_mp[stage].p_rat, tempVec->size);
  ED->EqD->stageData_d[stage].p->size = ED->EqD->stageData_mp[stage].p->size = tempVec->size;
  for (j = 0; j < tempVec->size; j++)
  { // copy tempVec to p & p_rat & copy back to p in double
    mpf_t_to_rat(ED->EqD->stageData_mp[stage].p_rat[j][0], tempVec->coord[j].r);
    mpf_t_to_rat(ED->EqD->stageData_mp[stage].p_rat[j][1], tempVec->coord[j].i);
    mpf_set_q(ED->EqD->stageData_mp[stage].p->coord[j].r, ED->EqD->stageData_mp[stage].p_rat[j][0]);
    mpf_set_q(ED->EqD->stageData_mp[stage].p->coord[j].i, ED->EqD->stageData_mp[stage].p_rat[j][1]);
    ED->EqD->stageData_d[stage].p->coord[j].r = mpq_get_d(ED->EqD->stageData_mp[stage].p_rat[j][0]);
    ED->EqD->stageData_d[stage].p->coord[j].i = mpq_get_d(ED->EqD->stageData_mp[stage].p_rat[j][1]);
  }

  // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
  transpose_mp(A, A);
  QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
  init_mat_d(ED->EqD->stageData_d[stage].B, ED->EqD->num_vars, depth_sum);
  init_mat_mp(ED->EqD->stageData_mp[stage].B, ED->EqD->num_vars, depth_sum);
  init_mat_rat(ED->EqD->stageData_mp[stage].B_rat, ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B->rows = ED->EqD->stageData_mp[stage].B->rows = ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B->cols = ED->EqD->stageData_mp[stage].B->cols = depth_sum;
  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_mp[stage].B->rows; j++)
    for (k = 0; k < ED->EqD->stageData_mp[stage].B->cols; k++)
    { // copy to B & B_rat and then copy back to B in double
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B_rat[j][k][0], Q->entry[j][k+count].r);
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B_rat[j][k][1], Q->entry[j][k+count].i);
      mpf_set_q(ED->EqD->stageData_mp[stage].B->entry[j][k].r, ED->EqD->stageData_mp[stage].B_rat[j][k][0]);
      mpf_set_q(ED->EqD->stageData_mp[stage].B->entry[j][k].i, ED->EqD->stageData_mp[stage].B_rat[j][k][1]);
      ED->EqD->stageData_d[stage].B->entry[j][k].r = mpq_get_d(ED->EqD->stageData_mp[stage].B_rat[j][k][0]);
      ED->EqD->stageData_d[stage].B->entry[j][k].i = mpq_get_d(ED->EqD->stageData_mp[stage].B_rat[j][k][1]);
    }

  // set back to the current precision
  initMP(ED->EqD->curr_precision);

  // setup B1 & p1
  k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
  change_size_mat_mp(A, k, 2 * ED->EqD->num_vars);
  A->rows = k;
  A->cols = 2 * ED->EqD->num_vars;  // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
  // setup tempMat & b
  change_size_mat_mp(tempMat, k, k);
  change_size_vec_mp(b, k);
  tempMat->rows = tempMat->cols = b->size = k;

  count = 0;
  // setup linears in the x variables
  for (j = depth_x; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients of the linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->EqD->coeff_rat[j][k][1]);
      set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for x variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_mp(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->BED_mp->patch.patchCoeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->BED_mp->patch.patchCoeff_rat[j][k][1]);
      set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
    }
    count++;
  }
  // setup linears in the y variables
  for (j = 0; j < depth_x; j++)
  { // copy over the coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_mp(&A->entry[count][k]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].i, ED->EqD->coeff_rat[j][k][1]);
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_mp(&A->entry[count][k]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].i, ED->EqD->coeff_rat[j][k][1]);
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for y variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_mp(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_mp(&A->entry[count][k]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].r, ED->BED_mp->patch.patchCoeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].i, ED->BED_mp->patch.patchCoeff_rat[j][k][1]);
    }
    count++;
  }

  // setup tempMat to find p1
  count = 0;
  // copy over the first n+p-depth_x cols related to the x variables
  for (j = 0; j < ED->EqD->num_funcs + ED->patch.num_patches - depth_x; j++)
  {
    for (k = 0; k < tempMat->rows; k++)
    {
      set_mp(&tempMat->entry[k][count], &A->entry[k][j]);
    }
    count++;
  }
  // copy over the first n+p-depth_y cols related to the y variables
  for (j = ED->EqD->num_vars; j < ED->EqD->num_vars + ED->EqD->num_funcs + ED->patch.num_patches - depth_y; j++)
  {
    for (k = 0; k < tempMat->rows; k++)
    {
      set_mp(&tempMat->entry[k][count], &A->entry[k][j]);
    }
    count++;
  }

  // set the global prec to the maximum precision
  initMP(max_prec);

  // setup tempVec
  if (matrixSolve_mp(tempVec, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }

  // setup p1
  init_vec_d(ED->EqD->stageData_d[stage].p1, 2 * ED->EqD->num_vars);
  init_vec_mp2(ED->EqD->stageData_mp[stage].p1, 2 * ED->EqD->num_vars, ED->EqD->curr_precision);
  init_vec_rat(ED->EqD->stageData_mp[stage].p1_rat, 2 * ED->EqD->num_vars);
  ED->EqD->stageData_d[stage].p1->size = ED->EqD->stageData_mp[stage].p1->size = 2 * ED->EqD->num_vars;
  count = 0;
  for (j = 0; j < 2 * ED->EqD->num_vars; j++)
    if (j < ED->EqD->num_funcs + ED->patch.num_patches - depth_x)
    {
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].p1_rat[j][0], tempVec->coord[count].r);
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].p1_rat[j][1], tempVec->coord[count].i);
      mpf_set_q(ED->EqD->stageData_mp[stage].p1->coord[j].r, ED->EqD->stageData_mp[stage].p1_rat[j][0]);
      mpf_set_q(ED->EqD->stageData_mp[stage].p1->coord[j].i, ED->EqD->stageData_mp[stage].p1_rat[j][1]);
      ED->EqD->stageData_d[stage].p1->coord[j].r = mpq_get_d(ED->EqD->stageData_mp[stage].p1_rat[j][0]);
      ED->EqD->stageData_d[stage].p1->coord[j].i = mpq_get_d(ED->EqD->stageData_mp[stage].p1_rat[j][1]);
      count++;
    }
    else if (j < ED->EqD->num_vars)
    {
      set_zero_rat(ED->EqD->stageData_mp[stage].p1_rat[j]);
      set_zero_mp(&ED->EqD->stageData_mp[stage].p1->coord[j]);
      set_zero_d(&ED->EqD->stageData_d[stage].p1->coord[j]);
    }
    else if (j < ED->EqD->num_vars + ED->EqD->num_funcs + ED->patch.num_patches - depth_y)
    {
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].p1_rat[j][0], tempVec->coord[count].r);
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].p1_rat[j][1], tempVec->coord[count].i);
      mpf_set_q(ED->EqD->stageData_mp[stage].p1->coord[j].r, ED->EqD->stageData_mp[stage].p1_rat[j][0]);
      mpf_set_q(ED->EqD->stageData_mp[stage].p1->coord[j].i, ED->EqD->stageData_mp[stage].p1_rat[j][1]);
      ED->EqD->stageData_d[stage].p1->coord[j].r = mpq_get_d(ED->EqD->stageData_mp[stage].p1_rat[j][0]);
      ED->EqD->stageData_d[stage].p1->coord[j].i = mpq_get_d(ED->EqD->stageData_mp[stage].p1_rat[j][1]);
      count++;
    }
    else
    {
      set_zero_rat(ED->EqD->stageData_mp[stage].p1_rat[j]);
      set_zero_mp(&ED->EqD->stageData_mp[stage].p1->coord[j]);
      set_zero_d(&ED->EqD->stageData_d[stage].p1->coord[j]);
    }

  // now, to find B1, do a QR decomposition on A^T and B is the extra columns of Q  
  transpose_mp(A, A);
  QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);

  // allocate space for B1_rat and then setup B1 & B1_rat
  init_mat_d(ED->EqD->stageData_d[stage].B1, 2 * ED->EqD->num_vars, depth_sum);
  init_mat_mp2(ED->EqD->stageData_mp[stage].B1, 2 * ED->EqD->num_vars, depth_sum, ED->EqD->curr_precision);
  init_mat_rat(ED->EqD->stageData_mp[stage].B1_rat, 2 * ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B1->rows = ED->EqD->stageData_mp[stage].B1->rows = 2 * ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B1->cols = ED->EqD->stageData_mp[stage].B1->cols = depth_sum;

  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_mp[stage].B1->rows; j++)
    for (k = 0; k < ED->EqD->stageData_mp[stage].B1->cols; k++)
    { // copy to B & B_rat and then copy back to B in double
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B1_rat[j][k][0], Q->entry[j][k+count].r);
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B1_rat[j][k][1], Q->entry[j][k+count].i);
      mpf_set_q(ED->EqD->stageData_mp[stage].B1->entry[j][k].r, ED->EqD->stageData_mp[stage].B1_rat[j][k][0]);
      mpf_set_q(ED->EqD->stageData_mp[stage].B1->entry[j][k].i, ED->EqD->stageData_mp[stage].B1_rat[j][k][1]);
      ED->EqD->stageData_d[stage].B1->entry[j][k].r = mpq_get_d(ED->EqD->stageData_mp[stage].B1_rat[j][k][0]);
      ED->EqD->stageData_d[stage].B1->entry[j][k].i = mpq_get_d(ED->EqD->stageData_mp[stage].B1_rat[j][k][1]);
    }

  // set back to the current precision
  initMP(ED->EqD->curr_precision);

  // now we do the same procedure to find p0 & B0
  k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
  change_size_mat_mp(A, k, 2 * ED->EqD->num_vars);
  A->rows = k;
  A->cols = 2 * ED->EqD->num_vars; // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
  // setup tempMat & b
  change_size_mat_mp(tempMat, k, k);
  change_size_vec_mp(b, k);
  tempMat->rows = tempMat->cols = b->size = k;

  count = 0;
  // setup x_j - y_j for j = depth_x to depth_sum
  for (j = depth_x; j < depth_sum; j++)
  {
    for (k = 0; k < ED->EqD->num_vars; k++)
      if (j == k)
      {
        set_one_mp(&A->entry[count][k]);
        neg_mp(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
      }
      else
      {
        set_zero_mp(&A->entry[count][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // setup linears in the x variables
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over the coefficients of the linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->EqD->coeff_rat[j][k][1]);
      set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for x variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_mp(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      mpf_set_q(A->entry[count][k].r, ED->BED_mp->patch.patchCoeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k].i, ED->BED_mp->patch.patchCoeff_rat[j][k][1]);
      set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
    }
    count++;
  }
  // setup x_j - y_j for j = 0 to depth_x
  for (j = 0; j < depth_x; j++)
  {
    for (k = 0; k < ED->EqD->num_vars; k++)
      if (j == k)
      {
        set_one_mp(&A->entry[count][k]);
        neg_mp(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
      }
      else
      {
        set_zero_mp(&A->entry[count][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  for (j = depth_sum; j < ED->EqD->num_funcs; j++)
  { // copy over coefficients of linear j
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_mp(&A->entry[count][k]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].r, ED->EqD->coeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].i, ED->EqD->coeff_rat[j][k][1]);
    }
    // set b[count] to 0
    set_zero_mp(&b->coord[count]);
    count++;
  }
  // copy over the patch coefficients for y variables
  for (j = 0; j < ED->patch.num_patches; j++)
  { // set b[count] to 1
    set_one_mp(&b->coord[count]);
    // copy patch
    for (k = 0; k < ED->EqD->num_vars; k++)
    {
      set_zero_mp(&A->entry[count][k]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].r, ED->BED_mp->patch.patchCoeff_rat[j][k][0]);
      mpf_set_q(A->entry[count][k + ED->EqD->num_vars].i, ED->BED_mp->patch.patchCoeff_rat[j][k][1]);
    }
    count++;
  }

  // set the global prec to the maximum precision
  initMP(max_prec);

  // setup tempMat to find p using a randomization process
  make_matrix_random_mp(R, A->cols, A->rows, max_prec);
  mat_mul_mp(tempMat, A, R);

  // find tempVec
  if (matrixSolve_mp(tempVec, tempMat, b))
  { // this should never happen!
    printf("ERROR: Problem solving a random matrix\n");
    bexit(ERROR_OTHER);
  }

  // undo the randomization
  mul_mat_vec_mp(tempVec, R, tempVec);

  // seutp p0
  init_vec_d(ED->EqD->stageData_d[stage].p0, 2 * ED->EqD->num_vars);
  init_vec_mp2(ED->EqD->stageData_mp[stage].p0, 2 * ED->EqD->num_vars, ED->EqD->curr_precision);
  init_vec_rat(ED->EqD->stageData_mp[stage].p0_rat, 2 * ED->EqD->num_vars);
  ED->EqD->stageData_d[stage].p0->size = ED->EqD->stageData_mp[stage].p0->size = 2 * ED->EqD->num_vars;
  for (j = 0; j < tempVec->size; j++)
  { // setup p0[j]
    mpf_t_to_rat(ED->EqD->stageData_mp[stage].p0_rat[j][0], tempVec->coord[j].r);
    mpf_t_to_rat(ED->EqD->stageData_mp[stage].p0_rat[j][1], tempVec->coord[j].i);
    mpf_set_q(ED->EqD->stageData_mp[stage].p0->coord[j].r, ED->EqD->stageData_mp[stage].p0_rat[j][0]);
    mpf_set_q(ED->EqD->stageData_mp[stage].p0->coord[j].i, ED->EqD->stageData_mp[stage].p0_rat[j][1]);
    ED->EqD->stageData_d[stage].p0->coord[j].r = mpq_get_d(ED->EqD->stageData_mp[stage].p0_rat[j][0]);
    ED->EqD->stageData_d[stage].p0->coord[j].i = mpq_get_d(ED->EqD->stageData_mp[stage].p0_rat[j][1]);
  }

  // now, to find B0, do a QR decomposition on A^T and B is the extra columns of Q
  transpose_mp(A, A);
  QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);

  // allocate space for B0
  init_mat_d(ED->EqD->stageData_d[stage].B0, 2 * ED->EqD->num_vars, depth_sum);
  init_mat_mp2(ED->EqD->stageData_mp[stage].B0, 2 * ED->EqD->num_vars, depth_sum, ED->EqD->curr_precision);
  init_mat_rat(ED->EqD->stageData_mp[stage].B0_rat, 2 * ED->EqD->num_vars, depth_sum);
  ED->EqD->stageData_d[stage].B0->rows = ED->EqD->stageData_mp[stage].B0->rows = 2 * ED->EqD->num_vars;
  ED->EqD->stageData_d[stage].B0->cols = ED->EqD->stageData_mp[stage].B0->cols = depth_sum;

  count = A->cols;
  for (j = 0; j < ED->EqD->stageData_mp[stage].B0->rows; j++)
    for (k = 0; k < ED->EqD->stageData_mp[stage].B0->cols; k++)
    { // copy to B & B_rat and then copy back to B in double
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B0_rat[j][k][0], Q->entry[j][k+count].r);
      mpf_t_to_rat(ED->EqD->stageData_mp[stage].B0_rat[j][k][1], Q->entry[j][k+count].i);
      mpf_set_q(ED->EqD->stageData_mp[stage].B0->entry[j][k].r, ED->EqD->stageData_mp[stage].B0_rat[j][k][0]);
      mpf_set_q(ED->EqD->stageData_mp[stage].B0->entry[j][k].i, ED->EqD->stageData_mp[stage].B0_rat[j][k][1]);
      ED->EqD->stageData_d[stage].B0->entry[j][k].r = mpq_get_d(ED->EqD->stageData_mp[stage].B0_rat[j][k][0]);
      ED->EqD->stageData_d[stage].B0->entry[j][k].i = mpq_get_d(ED->EqD->stageData_mp[stage].B0_rat[j][k][1]);
    }

  // set back to the current precision
  initMP(ED->EqD->curr_precision);

  // clear MP
  mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_mp(tempVec); clear_vec_mp(b);
  clear_mat_mp(tempMat); clear_mat_mp(A);
  clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

  free(str);

  return;
}

void setupEqbyEqNextStage_d(basic_eval_data_d *ED, int stage, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the 'stage'th stage of eq-by-eq                  *
* Setup a new gamma, patch, move the points to that patch and   *
* then setup B, p and the start points                          *
\***************************************************************/
{
  int i, j, k, depth_x, depth_y, depth_sum, count, indexI, indexJ, num_vars = ED->EqD->num_vars, num_paths;
  vec_d new_vars;
  mat_d B_transpose;

  init_vec_d(new_vars, 0);
  init_mat_d(B_transpose, 0, 0);

  // error checking
  if (stage <= 0 || stage >= ED->EqD->num_subsystems)
  {
    printf("ERROR: The stage to setup is incorrect!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup gamma
  if (MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(ED->EqD->gamma_d);
  }
  else
  { // setup gamma_d, gamma_mp & gamma_rat
    get_comp_rand_rat(ED->EqD->gamma_d, ED->EqD->gamma_mp, ED->EqD->gamma_rat, ED->EqD->curr_precision, max_prec, 0, 0);
  }

  // setup a new patch
  indexI = ED->patch.patchCoeff->rows;
  indexJ = ED->patch.patchCoeff->cols;

  // setup the patch coefficients - only change the ones that are not 0
  for (i = 0; i < indexI; i++)
    for (j = 0; j < indexJ; j++)
      if (ED->patch.patchCoeff->entry[i][j].r != 0 || ED->patch.patchCoeff->entry[i][j].i != 0)
      {
        if (MPType == 0)
        { // setup _d patch coefficients
          get_comp_rand_d(&ED->patch.patchCoeff->entry[i][j]);
        }
        else
        { // setup _d, _mp & _rat patch coefficients - only change the ones that are not 0
          get_comp_rand_rat(&ED->patch.patchCoeff->entry[i][j], &ED->BED_mp->patch.patchCoeff->entry[i][j], ED->BED_mp->patch.patchCoeff_rat[i][j], ED->BED_mp->patch.curr_prec, max_prec, 0, 0);
        }
      }

  // move the previous stage points to this new patch
  indexI = ED->EqD->stageData_d[stage-1].num_paths;
  for (i = 0; i < indexI; i++)
  {
    move_to_patch_d(ED->EqD->stageData_d[stage-1].endPts[i], ED->EqD->stageData_d[stage-1].endPts[i], ED);
  }
  // move the witness data points to this new patch
  indexI = ED->EqD->witnessData_d[stage].num_paths;
  for (i = 0; i < indexI; i++)
  {
    move_to_patch_d(ED->EqD->witnessData_d[stage].endPts[i], ED->EqD->witnessData_d[stage].endPts[i], ED);
  }

  // setup depth_x & depth_y
  depth_x = ED->EqD->stageData_d[stage].depth_x;
  depth_y = ED->EqD->stageData_d[stage].depth_y;
  depth_sum = depth_x + depth_y; 

  // calculate the number of variables for this stage
  if (ED->EqD->stageData_d[stage].useIntrinsicSlice)
    num_vars = depth_sum;
  else
    num_vars = 2 * ED->EqD->num_vars;

  // calculate number of paths for this stage
  num_paths = ED->EqD->stageData_d[stage].num_paths = ED->EqD->stageData_d[stage - 1].num_nonsing * ED->EqD->witnessData_d[stage].num_nonsing;
  
  // initialize the other counts
  ED->EqD->stageData_d[stage].num_sing = ED->EqD->stageData_d[stage].num_nonsing = ED->EqD->stageData_d[stage].num_inf = ED->EqD->stageData_d[stage].num_bad = 0;

  // allocate memory
  ED->EqD->stageData_d[stage].startPts = (point_d *)bmalloc(num_paths * sizeof(point_d));
  if (ED->EqD->stageData_d[stage].useIntrinsicSlice)
    ED->EqD->stageData_d[stage].endPts_in = (point_d *)bmalloc(num_paths * sizeof(point_d));
  else
    ED->EqD->stageData_d[stage].endPts_in = NULL;
  ED->EqD->stageData_d[stage].endPts = (point_d *)bmalloc(num_paths * sizeof(point_d));
  ED->EqD->stageData_d[stage].finalTs = (comp_d *)bmalloc(num_paths * sizeof(comp_d));
  ED->EqD->stageData_d[stage].condition_nums = (double *)bmalloc(num_paths * sizeof(double));
  ED->EqD->stageData_d[stage].endPt_retVals = (int *)bmalloc(num_paths * sizeof(int));
  ED->EqD->stageData_d[stage].endPt_types = (int *)bmalloc(num_paths * sizeof(int));
  ED->EqD->stageData_d[stage].higherDim = (int *)bmalloc(num_paths * sizeof(int));
  for (j = 0; j < num_paths; j++)
  { // initialize the memory
    init_point_d(ED->EqD->stageData_d[stage].startPts[j], num_vars);
    ED->EqD->stageData_d[stage].startPts[j]->size = num_vars;
    if (ED->EqD->stageData_d[stage].useIntrinsicSlice)
    {
      init_point_d(ED->EqD->stageData_d[stage].endPts_in[j], num_vars);
      ED->EqD->stageData_d[stage].endPts_in[j]->size = num_vars;
    }
    init_point_d(ED->EqD->stageData_d[stage].endPts[j], ED->EqD->num_vars);
    ED->EqD->stageData_d[stage].endPts[j]->size = ED->EqD->num_vars;

    set_zero_d(ED->EqD->stageData_d[stage].finalTs[j]);
    ED->EqD->stageData_d[stage].condition_nums[j] = ED->EqD->stageData_d[stage].endPt_retVals[j] = ED->EqD->stageData_d[stage].endPt_types[j] = ED->EqD->stageData_d[stage].higherDim[j] = 0;
  }

  if (MPType == 2)
  { // setup the first set of values
    ED->EqD->stageData_mp[stage].num_paths = ED->EqD->stageData_d[stage].num_paths;
    ED->EqD->stageData_mp[stage].num_sing = ED->EqD->stageData_mp[stage].num_nonsing = ED->EqD->stageData_mp[stage].num_inf = ED->EqD->stageData_mp[stage].num_bad = 0;

    // NULL out all other structures that will not be used
    ED->EqD->stageData_mp[stage].startPts = NULL;
    ED->EqD->stageData_mp[stage].endPts_in = NULL;
    ED->EqD->stageData_mp[stage].endPts = NULL;
    ED->EqD->stageData_mp[stage].finalTs = NULL;
    ED->EqD->stageData_mp[stage].condition_nums = NULL;
    ED->EqD->stageData_mp[stage].endPt_retVals = NULL;
    ED->EqD->stageData_mp[stage].endPt_types = NULL;
    ED->EqD->stageData_mp[stage].higherDim = NULL;
  }

  // setup to find the intrinsic slice, if needed
  if (ED->EqD->stageData_d[stage].useIntrinsicSlice)
  { // setup for intrinisic tracking
    if (MPType == 2)
    { // setup for intrinsic slicing using AMP
      setupInt_amp(ED, stage, max_prec);
    }
    else
    { // setup for intrinsic slicing in only double precision
      setupInt_d(ED, stage);
    }

    // now that we have the slice, setup the start points on this slice
    transpose_d(B_transpose, ED->EqD->stageData_d[stage].B1);

    // setup the start points
    count = indexI = indexJ = 0;
    num_vars = ED->EqD->num_vars;
    for (i = 0; i < ED->EqD->stageData_d[stage - 1].num_nonsing; i++)
    { // find the next good path from the previous stage
      while (ED->EqD->stageData_d[stage - 1].endPt_types[indexI] != MOVE_TO_NEXT)
        indexI++;

      // setup the top of new_vars
      point_cp_d(new_vars, ED->EqD->stageData_d[stage-1].endPts[indexI]);
      increase_size_vec_d(new_vars, 2 * num_vars);
      new_vars->size = 2 * num_vars;
    
      indexJ = 0;
      for (j = 0; j < ED->EqD->witnessData_d[stage].num_nonsing; j++)
      { // find the next good path from the next subsystem
        while (ED->EqD->witnessData_d[stage].endPt_types[indexJ] != MOVE_TO_NEXT)
          indexJ++;

        // setup the bottom of new_vars
        for (k = num_vars; k < 2 * num_vars; k++)
        {
          set_d(&new_vars->coord[k], &ED->EqD->witnessData_d[stage].endPts[indexJ]->coord[k - num_vars]);
        }
        new_vars->size = 2 * num_vars;

        // turn new_vars into the intrinsic coordinates
        extrinsicToIntrinsic_d(ED->EqD->stageData_d[stage].startPts[count], new_vars, B_transpose, ED->EqD->stageData_d[stage].p1);

        // increment indexJ
        indexJ++;
        // increment count
        count++;
      }
      // increment indexI
      indexI++;
    }
  }
  else
  { // setup for extrinsic tracking
    count = indexI = indexJ = 0;
    for (i = 0; i < ED->EqD->stageData_d[stage - 1].num_nonsing; i++)
    { // find the next good path from the previous stage
      while (ED->EqD->stageData_d[stage - 1].endPt_types[indexI] != MOVE_TO_NEXT)
        indexI++;

      indexJ = 0;
      for (j = 0; j < ED->EqD->witnessData_d[stage].num_nonsing; j++)
      { // find the next good path from the next subsystem
        while (ED->EqD->witnessData_d[stage].endPt_types[indexJ] != MOVE_TO_NEXT)
          indexJ++;

        // setup the start point
        for (k = 0; k < ED->EqD->num_vars; k++)
        {
          set_d(&ED->EqD->stageData_d[stage].startPts[count]->coord[k], &ED->EqD->stageData_d[stage-1].endPts[indexI]->coord[k]);
          set_d(&ED->EqD->stageData_d[stage].startPts[count]->coord[k+ED->EqD->num_vars], &ED->EqD->witnessData_d[stage].endPts[indexJ]->coord[k]);
        }

        // increment indexJ
        indexJ++;
        // increment count
        count++;
      }
      // increment indexI
      indexI++;
    }
  }

  clear_vec_d(new_vars);
  clear_mat_d(B_transpose);

  return;
}

void clearEqbyEqStageData_d(eqData_t *EqD, int stage, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'stage' stage data information               *
\***************************************************************/
{
  int i, num_paths = EqD->stageData_d[stage].num_paths;

  // clear memory
  for (i = num_paths - 1; i >= 0; i--)
  {
    if (EqD->stageData_d[stage].startPts != NULL)
      clear_point_d(EqD->stageData_d[stage].startPts[i]);
    if (EqD->stageData_d[stage].useIntrinsicSlice && EqD->stageData_d[stage].endPts_in != NULL)
      clear_point_d(EqD->stageData_d[stage].endPts_in[i]);
    if (EqD->stageData_d[stage].endPts != NULL)
      clear_point_d(EqD->stageData_d[stage].endPts[i]);
  }
  free(EqD->stageData_d[stage].startPts);
  if (EqD->stageData_d[stage].useIntrinsicSlice)
    free(EqD->stageData_d[stage].endPts_in);
  free(EqD->stageData_d[stage].endPts);
  free(EqD->stageData_d[stage].finalTs);
  free(EqD->stageData_d[stage].condition_nums);
  free(EqD->stageData_d[stage].endPt_retVals);
  free(EqD->stageData_d[stage].endPt_types);
  free(EqD->stageData_d[stage].higherDim);

  if (MPType == 2)
  { 
    if (EqD->stageData_d[stage].useIntrinsicSlice)
    { // clear B
      clear_mat(EqD->stageData_d[stage].B, EqD->stageData_mp[stage].B, EqD->stageData_mp[stage].B_rat, MPType);
      // clear p
      clear_vec(EqD->stageData_d[stage].p, EqD->stageData_mp[stage].p, EqD->stageData_mp[stage].p_rat, MPType);

      // clear B1, B0, p1, p0, if needed
      if (EqD->stageData_mp[stage].B1_rat != NULL)
      { // clear B1
        clear_mat(EqD->stageData_d[stage].B1, EqD->stageData_mp[stage].B1, EqD->stageData_mp[stage].B1_rat, MPType);
      }
      if (EqD->stageData_mp[stage].B0_rat != NULL)
      { // clear B0
        clear_mat(EqD->stageData_d[stage].B0, EqD->stageData_mp[stage].B0, EqD->stageData_mp[stage].B0_rat, MPType);
      }
      if (EqD->stageData_mp[stage].p1_rat != NULL)
      { // clear p1
        clear_vec(EqD->stageData_d[stage].p1, EqD->stageData_mp[stage].p1, EqD->stageData_mp[stage].p1_rat, MPType);
      }
      if (EqD->stageData_mp[stage].p0_rat != NULL)
      { // clear p0
        clear_vec(EqD->stageData_d[stage].p0, EqD->stageData_mp[stage].p0, EqD->stageData_mp[stage].p0_rat, MPType);
      }
    }

    // set other pointers to NULL
    EqD->stageData_mp[stage].startPts = NULL;
    EqD->stageData_mp[stage].endPts_in = NULL;
    EqD->stageData_mp[stage].endPts = NULL;
    EqD->stageData_mp[stage].finalTs = NULL;
    EqD->stageData_mp[stage].condition_nums = NULL;
    EqD->stageData_mp[stage].endPt_retVals = NULL;
    EqD->stageData_mp[stage].endPt_types = NULL;
    EqD->stageData_mp[stage].higherDim = NULL;
  }
  else if (EqD->stageData_d[stage].useIntrinsicSlice)
  { // clear B
    clear_mat_d(EqD->stageData_d[stage].B);
    // clear p
    clear_vec_d(EqD->stageData_d[stage].p);

    // clear B1, B0, p1, p0, if needed
    clear_mat_d(EqD->stageData_d[stage].B1);
    clear_mat_d(EqD->stageData_d[stage].B0);
    clear_vec_d(EqD->stageData_d[stage].p1);
    clear_vec_d(EqD->stageData_d[stage].p0);
  }


  return;
}

void clearEqbyEqWitnessData_d(eqData_t *EqD, int subsystem, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'subsystem' witness data information         *
\***************************************************************/
{
  int i, j, num_paths = EqD->witnessData_d[subsystem].num_paths;

  // clear startSystemCoeff
  for (i = EqD->witnessData_d[subsystem].num_linears - 1; i >= 0; i--)
    free(EqD->witnessData_d[subsystem].startSystemCoeff[i]);
  free(EqD->witnessData_d[subsystem].startSystemCoeff);  

  if (EqD->witnessData_d[subsystem].startPts != NULL)
  {
    for (i = num_paths - 1; i >= 0; i--)
    {
      clear_point_d(EqD->witnessData_d[subsystem].startPts[i]);
      clear_point_d(EqD->witnessData_d[subsystem].endPts_in[i]);
      clear_point_d(EqD->witnessData_d[subsystem].endPts[i]);
    }
  }
  // clear the other memory 
  free(EqD->witnessData_d[subsystem].startPts);
  free(EqD->witnessData_d[subsystem].endPts_in);
  free(EqD->witnessData_d[subsystem].endPts);
  free(EqD->witnessData_d[subsystem].finalTs);
  free(EqD->witnessData_d[subsystem].condition_nums);
  free(EqD->witnessData_d[subsystem].endPt_retVals);
  free(EqD->witnessData_d[subsystem].endPt_types);
  free(EqD->witnessData_d[subsystem].higherDim);

  if (MPType == 2)
  { // clear MP structures
    for (i = EqD->witnessData_mp[subsystem].num_linears - 1; i >= 0; i--)
    {
      for (j = EqD->witnessData_mp[subsystem].depth - 1; j >= 0; j--)
      {  
        mpq_clear(EqD->witnessData_mp[subsystem].startSystemCoeff_rat[i][j][0]);
        mpq_clear(EqD->witnessData_mp[subsystem].startSystemCoeff_rat[i][j][1]);
        free(EqD->witnessData_mp[subsystem].startSystemCoeff_rat[i][j]);

        clear_mp(EqD->witnessData_mp[subsystem].startSystemCoeff[i][j]);
      }
      free(EqD->witnessData_mp[subsystem].startSystemCoeff_rat[i]);
      free(EqD->witnessData_mp[subsystem].startSystemCoeff[i]);
    }
    free(EqD->witnessData_mp[subsystem].startSystemCoeff_rat);
    free(EqD->witnessData_mp[subsystem].startSystemCoeff);

    // clear B
    clear_mat(EqD->witnessData_d[subsystem].B, EqD->witnessData_mp[subsystem].B, EqD->witnessData_mp[subsystem].B_rat, MPType);

    // clear p
    clear_vec(EqD->witnessData_d[subsystem].p, EqD->witnessData_mp[subsystem].p, EqD->witnessData_mp[subsystem].p_rat, MPType);

    // set other pointers to NULL
    EqD->witnessData_mp[subsystem].startPts = NULL;
    EqD->witnessData_mp[subsystem].endPts_in = NULL;
    EqD->witnessData_mp[subsystem].endPts = NULL;
    EqD->witnessData_mp[subsystem].finalTs = NULL;
    EqD->witnessData_mp[subsystem].condition_nums = NULL;
    EqD->witnessData_mp[subsystem].endPt_retVals = NULL;
    EqD->witnessData_mp[subsystem].endPt_types = NULL;
    EqD->witnessData_mp[subsystem].higherDim = NULL;
  }
  else
  { // clear B
    clear_mat_d(EqD->witnessData_d[subsystem].B);

    // clear p
    clear_vec_d(EqD->witnessData_d[subsystem].p);
  }
 
  return;
}

void setupEqbyEqWitnessData_amp(eqData_t *EqD, patch_eval_data_mp *PD_mp, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the internal structures of witness data          *
\***************************************************************/
// ASSUME 1-hom!!!
{
  int i, j, k, r, start, num_digits = prec_to_digits(max_prec), *curr_deg = NULL;
  size_t size;
  char *str = NULL;
  vec_d b_d;
  mat_d tempMat_d;
  vec_mp tempVec, b;
  mat_mp tempMat, A, Q, R, P;
  mpf_t tol_pivot, tol_sign, largeChange;

  // initialize MP
  mpf_init2(tol_pivot, max_prec); mpf_init2(tol_sign, max_prec); mpf_init2(largeChange, max_prec);

  // initialize
  init_vec_d(b_d, 0);
  init_mat_d(tempMat_d, 0, 0);
  init_vec_mp2(tempVec, 0, max_prec); init_vec_mp2(b, 0, max_prec);
  init_mat_mp2(tempMat, 0, 0, max_prec); init_mat_mp2(A, 0, 0, max_prec);
  init_mat_mp2(Q, 0, 0, max_prec); init_mat_mp2(R, 0, 0, max_prec); init_mat_mp2(P, 0, 0, max_prec);

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

  for (i = 0; i < EqD->num_subsystems; i++)
  { // setup the ith witness data
    EqD->witnessData_mp[i].num_paths = EqD->witnessData_d[i].num_paths;

    // setup the linears
    EqD->witnessData_mp[i].num_linears = EqD->witnessData_d[i].num_linears;
    EqD->witnessData_mp[i].startSystemCoeff = (comp_mp **)bmalloc(EqD->witnessData_mp[i].num_linears * sizeof(comp_mp *));
    EqD->witnessData_mp[i].startSystemCoeff_rat = (mpq_t ***)bmalloc(EqD->witnessData_mp[i].num_linears * sizeof(mpq_t **));

    for (j = 0; j < EqD->witnessData_mp[i].num_linears; j++)
    { // setup the jth linear
      EqD->witnessData_mp[i].startSystemCoeff[j] = (comp_mp *)bmalloc(EqD->witnessData_mp[i].depth * sizeof(comp_mp));
      EqD->witnessData_mp[i].startSystemCoeff_rat[j] = (mpq_t **)bmalloc(EqD->witnessData_mp[i].depth * sizeof(mpq_t *));
      for (k = 0; k < EqD->witnessData_mp[i].depth; k++)
      { // allocate for real & imag
        EqD->witnessData_mp[i].startSystemCoeff_rat[j][k] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
        // get a random number
        get_comp_rand_rat(EqD->witnessData_d[i].startSystemCoeff[j][k], EqD->witnessData_mp[i].startSystemCoeff[j][k], EqD->witnessData_mp[i].startSystemCoeff_rat[j][k], EqD->curr_precision, max_prec, 1, 1);
      }
    }

    // setup A
    k = EqD->num_funcs - EqD->witnessData_mp[i].depth + PD_mp->num_patches;
    change_size_mat_mp(A, k, EqD->num_vars);
    A->rows = k;
    A->cols = EqD->num_vars;  // == EqD->num_funcs + PD_mp->num_patches
    // setup tempMat & b
    change_size_mat_mp(tempMat, k, k);
    change_size_vec_mp(b, k);
    tempMat->rows = tempMat->cols = b->size = k;

    // setup to find the intrinsic slice
    r = 0;
    for (j = 0; j < EqD->num_funcs; j++)
      if (j < EqD->witnessData_mp[i].startFunction || j >= EqD->witnessData_mp[i].startFunction + EqD->witnessData_mp[i].depth)
      { // need to copy over the coefficients
        for (k = 0; k < EqD->num_vars; k++)
        {
          mpf_set_q(A->entry[r][k].r, EqD->coeff_rat[j][k][0]);
          mpf_set_q(A->entry[r][k].i, EqD->coeff_rat[j][k][1]);
          if (k < tempMat->cols)
          { // copy to tempMat
            set_mp(&tempMat->entry[r][k], &A->entry[r][k]);
          }
        }
        // set b[r] to 0
        set_zero_mp(&b->coord[r]);
        r++;
      }
    // copy over the patch coefficients
    for (j = 0; j < PD_mp->num_patches; j++)
    { // set b[r] to 1
      mpf_set_ui(b->coord[r].r, 1);
      mpf_set_ui(b->coord[r].i, 0);
      // copy patch
      for (k = 0; k < EqD->num_vars; k++)
      {
        mpf_set_q(A->entry[r][k].r, PD_mp->patchCoeff_rat[j][k][0]);
        mpf_set_q(A->entry[r][k].i, PD_mp->patchCoeff_rat[j][k][1]);
        if (k < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[r][k], &A->entry[r][k]);
        }
      }
      r++;
    }

    // set the global prec to the maximum precision
    initMP(max_prec);

    // setup p - putting 0 in the extra positions
    if (matrixSolve_mp(tempVec, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // correctly setup tempVec
    start = tempMat->cols + EqD->witnessData_mp[i].depth;
    increase_size_vec_mp(tempVec, start);
    tempVec->size = start;
    for (j = tempMat->cols; j < start; j++)
    {
      set_zero_mp(&tempVec->coord[j]);
    }

    // setup p
    init_vec_mp2(EqD->witnessData_mp[i].p, start, EqD->curr_precision);
    init_vec_rat(EqD->witnessData_mp[i].p_rat, start);
    EqD->witnessData_mp[i].p->size = start;
    for (j = 0; j < EqD->witnessData_mp[i].p->size; j++)
    { // copy tempVec to p & p_rat & copy back to p in double
      mpf_t_to_rat(EqD->witnessData_mp[i].p_rat[j][0], tempVec->coord[j].r);
      mpf_t_to_rat(EqD->witnessData_mp[i].p_rat[j][1], tempVec->coord[j].i);
      mpf_set_q(EqD->witnessData_mp[i].p->coord[j].r, EqD->witnessData_mp[i].p_rat[j][0]);
      mpf_set_q(EqD->witnessData_mp[i].p->coord[j].i, EqD->witnessData_mp[i].p_rat[j][1]);
      EqD->witnessData_d[i].p->coord[j].r = mpq_get_d(EqD->witnessData_mp[i].p_rat[j][0]);
      EqD->witnessData_d[i].p->coord[j].i = mpq_get_d(EqD->witnessData_mp[i].p_rat[j][1]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);

    // setup B
    init_mat_mp2(EqD->witnessData_mp[i].B, EqD->num_vars, EqD->witnessData_mp[i].depth, EqD->curr_precision);
    init_mat_rat(EqD->witnessData_mp[i].B_rat, EqD->num_vars, EqD->witnessData_mp[i].depth);
    EqD->witnessData_mp[i].B->rows = EqD->num_vars;
    EqD->witnessData_mp[i].B->cols = EqD->witnessData_mp[i].depth;
    start = A->cols;

    for (j = 0; j < EqD->witnessData_mp[i].B->rows; j++)
      for (k = 0; k < EqD->witnessData_mp[i].B->cols; k++)
      { // copy to B & B_rat and then copy back to B in double
        mpf_t_to_rat(EqD->witnessData_mp[i].B_rat[j][k][0], Q->entry[j][k+start].r);
        mpf_t_to_rat(EqD->witnessData_mp[i].B_rat[j][k][1], Q->entry[j][k+start].i);
        mpf_set_q(EqD->witnessData_mp[i].B->entry[j][k].r, EqD->witnessData_mp[i].B_rat[j][k][0]);
        mpf_set_q(EqD->witnessData_mp[i].B->entry[j][k].i, EqD->witnessData_mp[i].B_rat[j][k][1]);
        EqD->witnessData_d[i].B->entry[j][k].r = mpq_get_d(EqD->witnessData_mp[i].B_rat[j][k][0]);
        EqD->witnessData_d[i].B->entry[j][k].i = mpq_get_d(EqD->witnessData_mp[i].B_rat[j][k][1]);
      }

    // set back to the current precision
    initMP(EqD->curr_precision);

    // NULL out all other structures
    EqD->witnessData_mp[i].startPts = NULL;
    EqD->witnessData_mp[i].endPts_in = NULL;
    EqD->witnessData_mp[i].endPts = NULL;
    EqD->witnessData_mp[i].finalTs = NULL;
    EqD->witnessData_mp[i].condition_nums = NULL;
    EqD->witnessData_mp[i].endPt_retVals = NULL;
    EqD->witnessData_mp[i].endPt_types = NULL;
    EqD->witnessData_mp[i].higherDim = NULL;

    // setup the start points in double precision again since the coefficients have changed
    curr_deg = (int *)brealloc(curr_deg, EqD->witnessData_d[i].depth * sizeof(int));
    for (j = 0; j < EqD->witnessData_d[i].depth; j++)
      curr_deg[j] = 0;

    change_size_mat_d(tempMat_d, EqD->witnessData_d[i].depth, EqD->witnessData_d[i].depth);
    change_size_vec_d(b_d, EqD->witnessData_d[i].depth);
    tempMat_d->rows = tempMat_d->cols = b_d->size = EqD->witnessData_d[i].depth;
    // setup b
    for (k = 0; k < b_d->size; k++)
    {
      set_one_d(&b_d->coord[k]);

    }
    for (j = 0; j < EqD->witnessData_d[i].num_paths; j++)
    { // setup tempMat
      start = 0;
      for (k = 0; k < tempMat_d->rows; k++)
      { // copy over the entries
        for (r = 0; r < tempMat_d->cols; r++)
        {
          set_d(&tempMat_d->entry[k][r], EqD->witnessData_d[i].startSystemCoeff[start + curr_deg[k]][r]);
        }
        // update start
        start += EqD->degrees[EqD->witnessData_d[i].startFunction + k][0];
      }
      // solve for the start point
      if (matrixSolve_d(EqD->witnessData_d[i].startPts[j], tempMat_d, b_d))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }
      // update curr_deg
      for (k = EqD->witnessData_d[i].depth - 1; k >= 0; k--)
      { // check to see if we are at the top
        if (curr_deg[k] == EqD->degrees[EqD->witnessData_d[i].startFunction + k][0] - 1)
        { // set to 0
          curr_deg[k] = 0;
        }
        else
        { // increment this location and exit loop
          curr_deg[k]++;
          k = -10;
        }
      }
    }
  }

  // release memory
  free(str);
  free(curr_deg);

  // clear MP
  mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_d(b_d);
  clear_mat_d(tempMat_d);
  clear_vec_mp(tempVec); clear_vec_mp(b);
  clear_mat_mp(tempMat); clear_mat_mp(A);
  clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

  return;
}

/////// MULTI PRECISION ////////////

int eqbyeq_setup_mp(FILE **OUT, char *outName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, char *pointsIN, char *pointsOUT, char *depthFile, double intrinsicCutoffMultiplier)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of original variables                   *
* NOTES: setup for zero dimensional tracking using eq-by-eq     *
\***************************************************************/
{
  int num_orig_vars, intrinsicCutoff;
  FILE *midOUT = NULL;

  // allocate space
  ED->EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));

  // setup the function
  num_orig_vars = zero_dim_basic_setup_mp(OUT, outName, &midOUT, "midpath_data", T, ED, dummyProg, &ED->EqD->startSub, &ED->EqD->endSub, &ED->EqD->startFunc, &ED->EqD->endFunc, &ED->EqD->startJvsub, &ED->EqD->endJvsub, &ED->EqD->startJv, &ED->EqD->endJv, &ED->EqD->subFuncsBelow, eval_mp, preprocFile, degreeFile, 0, pointsIN, pointsOUT);

  // midOUT does not need setup
  fclose(midOUT);

  // verify that we are using only 1 homogenous variable group
  if (ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: The equation-by-equation method of Sommese, Verschelde and Wampler is only implemented in Bertini for systems with\nonly one variable group.");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup whether there are changes & the number of subfunctions
  ED->EqD->noChanges = ED->squareSystem.noChanges;
  ED->EqD->numSubFuncs = ED->squareSystem.Prog->numSubfuncs;

  // setup the random numbers needed for eq-by-eq
  setupEqbyEqRandom_mp(ED->EqD, T, degreeFile, &ED->preProcData, ED->squareSystem.new_degrees);

  intrinsicCutoff = floor(intrinsicCutoffMultiplier * num_orig_vars);

  // allocate space for the witness data and the stage data
  setupEqbyEqStages_mp(ED->EqD, T, depthFile, intrinsicCutoff);

  // setup the witness data
  setupEqbyEqWitnessData_mp(ED->EqD, &ED->patch);

  return num_orig_vars;
}

void setupEqbyEqRandom_mp(eqData_t *EqD, tracker_config_t *T, char *degreeFile, preproc_data *PPD, int *new_degrees)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the random numbers for eq-by-eq                  *
\***************************************************************/
// ASSUME 1-hom!
{
  int i, j, k, beg, num_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;

  FILE *degIN = fopen(degreeFile, "r"); // open the file to read in the degrees
  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  EqD->curr_precision = T->Precision; // the fixed precision
  EqD->num_funcs = T->numVars - num_var_gps; // calculate the true number of functions (we can ignore the projective transformations)
  EqD->num_var_gps = num_var_gps; // set the number of variable groups
  EqD->num_vars = T->numVars;     // set the number of variables

  // setup the degrees and allocate space for the random numbers
  EqD->degrees = (int **)bmalloc(EqD->num_funcs * sizeof(int *));
  EqD->coeff_mp = (comp_mp **)bmalloc(EqD->num_funcs * sizeof(comp_mp *));

  // read in the degrees and allocate the space for the random numbers - the patches will be handled by the function evaluator and not in eq-by-eq
  for (i = 0; i < EqD->num_funcs; i++)
  { // setup degrees
    EqD->degrees[i] = (int *)bmalloc(num_var_gps * sizeof(int));

    if (num_var_gps == 1)
    { // the functions could have been permuted so we use new_degrees
      EqD->degrees[i][0] = new_degrees[i];
    }
    else
    { // need to read in all of the m-hom degrees from degIN
      for (j = 0; j < num_var_gps; j++)
        fscanf(degIN, "%d\n", &EqD->degrees[i][j]);
      fscanf(degIN, "\n"); // extra "new line" character
    }

    // allocate based on the total number of variables
    EqD->coeff_mp[i] = (comp_mp *)bmalloc(T->numVars * sizeof(comp_mp));

  }
  // close degreeFile
  fclose(degIN);

  // now that we have the space, initialize it and generate the random numbers

  // generate the random gamma
  init_mp(EqD->gamma_mp);
  get_comp_rand_mp(EqD->gamma_mp);

  // initialize all coeff to zero
  for (i = 0; i < EqD->num_funcs; i++)
    for (j = 0; j < T->numVars; j++)
    {
      init_mp(EqD->coeff_mp[i][j]);
      set_zero_mp(EqD->coeff_mp[i][j]);
    }

  // now do random numbers for the appropriate coeff
  for (i = 0; i < EqD->num_funcs; i++)
  { // need to make the linears truly linears in only the variable groups - zero for the other coeff
    if (PPD->type[0])
    { // this is a regular variable group - thus it has a homogenous coordinate that was created
      beg = 0;

      get_comp_rand_mp(EqD->coeff_mp[i][beg]);  // get random for hom coord for this variable group
    }
    beg = PPD->num_var_gp;

    for (k = 0; k < PPD->size[0]; k++)
      get_comp_rand_mp(EqD->coeff_mp[i][k + beg]); // get random for variables that are in this variable group
  }

  return;
}

void setupEqbyEqStages_mp(eqData_t *EqD, tracker_config_t *T, char *depthFile, int intrinsicCutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: allocate memory for the witness data and the stages    *
\***************************************************************/
{
  int i, j, *depths = NULL;
  FILE *DEPTH = fopen(depthFile, "r");

  if (DEPTH == NULL)
  { // setup with each depth as 1
    EqD->num_subsystems = EqD->num_funcs;

    // allocate space
    EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_mp));
    EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_mp));

    // setup the basic info for the witness data
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      EqD->witnessData_mp[i].startFunction = i;
      EqD->witnessData_mp[i].depth = 1;
    }
  }
  else
  { // read in the depths
    EqD->num_subsystems = 1;
    fscanf(DEPTH, "%d\n", &EqD->num_subsystems); // read in the number of subsystems that will be needed

    // error checking
    if (EqD->num_subsystems <= 0)
    {
      printf("ERROR: The number of subsystems must be > 0!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // setup depths
    depths = (int *)bmalloc(EqD->num_subsystems * sizeof(int));
    j = 0;
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      depths[i] = 1;
      fscanf(DEPTH, "%d\n", &depths[i]);
      if (depths[i] <= 0)
      {
        printf("ERROR: The depths must be > 0!\n");
        bexit(ERROR_CONFIGURATION);
      }
      j += depths[i]; // sum up the total depths
    }

    if (j != EqD->num_funcs) // the total depth is not correct
    {
      printf("ERROR: The number of functions (%d) is not equal to the total depth (%d)!\n", EqD->num_funcs, j);
      bexit(ERROR_CONFIGURATION);
    }

    // allocate space
    EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqWitnessData_mp));
    EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD->num_subsystems * sizeof(eqStageData_mp));

    // setup the basic info for the witness data
    j = 0;
    for (i = 0; i < EqD->num_subsystems; i++)
    {
      EqD->witnessData_mp[i].startFunction = j;
      EqD->witnessData_mp[i].depth = depths[i];

      j += depths[i];
    }

    free(depths);
    fclose(DEPTH);
  }

  // setup the basic info for the stages
  for (i = 0; i < EqD->num_subsystems; i++)
  { // setup depth_x
    if (i == 0)
      EqD->stageData_mp[i].depth_x = 0;
    else
      EqD->stageData_mp[i].depth_x = EqD->stageData_mp[i-1].depth_x + EqD->stageData_mp[i-1].depth_y;
    // setup depth_y
    EqD->stageData_mp[i].depth_y = EqD->witnessData_mp[i].depth;

    // determine if we are using intrinsic slice
    if (i > 0 && EqD->stageData_mp[i].depth_x + EqD->stageData_mp[i].depth_y > intrinsicCutoff)
      EqD->stageData_mp[i].useIntrinsicSlice = 0;
    else
      EqD->stageData_mp[i].useIntrinsicSlice = 1;
  }

  return;
}

void setupEqbyEqWitnessData_mp(eqData_t *EqD, patch_eval_data_mp *PD_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the internal structures of witness data          *
\***************************************************************/
// ASSUME 1-hom!!!
{
  int i, j, k, r, start, num_digits = prec_to_digits(EqD->curr_precision), *curr_deg = NULL;
  size_t size;
  char *str = NULL;
  vec_mp b;
  mat_mp tempMat, A, Q, R, P;
  mpf_t tol_pivot, tol_sign, largeChange;

  // initialize MP
  mpf_init(tol_pivot); mpf_init(tol_sign); mpf_init(largeChange);
  init_vec_mp(b, 0);
  init_mat_mp(tempMat, 0, 0); init_mat_mp(A, 0, 0); init_mat_mp(Q, 0, 0);
  init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0);

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

  for (i = 0; i < EqD->num_subsystems; i++)
  { // setup the ith witness data

    // initialize other values
    EqD->witnessData_mp[i].num_sing = EqD->witnessData_mp[i].num_nonsing = EqD->witnessData_mp[i].num_inf = EqD->witnessData_mp[i].num_bad = 0;

    // setup the number of paths
    EqD->witnessData_mp[i].num_linears = 0; // the sum of the degrees
    EqD->witnessData_mp[i].num_paths = 1; // the multiplication of the degrees
    start = EqD->witnessData_mp[i].startFunction;
    for (j = 0; j < EqD->witnessData_mp[i].depth; j++)
    {
      EqD->witnessData_mp[i].num_linears += EqD->degrees[start + j][0];
      EqD->witnessData_mp[i].num_paths *= EqD->degrees[start + j][0];
    }

    // setup the start system information - we will be using intrinsic slices to make it only 'depth' number of variables to track
    EqD->witnessData_mp[i].startSystemCoeff = (comp_mp **)bmalloc(EqD->witnessData_mp[i].num_linears * sizeof(comp_mp *));
    for (j = 0; j < EqD->witnessData_mp[i].num_linears; j++)
    { // setup the jth linear - assuming here we are in the 1-hom case
      EqD->witnessData_mp[i].startSystemCoeff[j] = (comp_mp *)bmalloc(EqD->witnessData_mp[i].depth * sizeof(comp_mp));
      for (k = 0; k < EqD->witnessData_mp[i].depth; k++)
      { // get a random number
        init_mp(EqD->witnessData_mp[i].startSystemCoeff[j][k]);
        get_comp_rand_mp(EqD->witnessData_mp[i].startSystemCoeff[j][k]);
      }
    }

    // setup A
    k = EqD->num_funcs - EqD->witnessData_mp[i].depth + PD_mp->num_patches;
    change_size_mat_mp(A, k, EqD->num_vars);
    A->rows = k;
    A->cols = EqD->num_vars; // == EqD->num_funcs + PD_mp->num_patches
    // setup tempMat & b
    change_size_mat_mp(tempMat, k, k);
    change_size_vec_mp(b, k);
    tempMat->rows = tempMat->cols = b->size = k;

    // setup to find the intrinsic slice
    r = 0;
    for (j = 0; j < EqD->num_funcs; j++)
      if (j < EqD->witnessData_mp[i].startFunction || j >= EqD->witnessData_mp[i].startFunction + EqD->witnessData_mp[i].depth)
      { // need to copy over the coefficients
        for (k = 0; k < EqD->num_vars; k++)
        {
          set_mp(&A->entry[r][k], EqD->coeff_mp[j][k]);
          if (k < tempMat->cols)
          { // copy to tempMat
            set_mp(&tempMat->entry[r][k], &A->entry[r][k]);
          }
        }
        // set b[r] to 0
        set_zero_mp(&b->coord[r]);
        r++;
      }
    // copy over the patch coefficients
    for (j = 0; j < PD_mp->num_patches; j++)
    { // set b[r] to 1
      set_one_mp(&b->coord[r]);
      // copy patch
      for (k = 0; k < EqD->num_vars; k++)
      {
        set_mp(&A->entry[r][k], &PD_mp->patchCoeff->entry[j][k]);
        if (k < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[r][k], &A->entry[r][k]);
        }
      }
      r++;
    }

    // setup p - putting 0 in the extra positions
    // by doing it this way, we have a standard way of constructing p from coeff
    init_vec_mp(EqD->witnessData_mp[i].p, 0);
    if (matrixSolve_mp(EqD->witnessData_mp[i].p, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // correctly setup p
    start = tempMat->cols + EqD->witnessData_mp[i].depth;
    increase_size_vec_mp(EqD->witnessData_mp[i].p, start);
    EqD->witnessData_mp[i].p->size = start;
    for (j = tempMat->cols; j < start; j++)
    {
      set_zero_mp(&EqD->witnessData_mp[i].p->coord[j]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    init_mat_mp(EqD->witnessData_mp[i].B, EqD->num_vars, EqD->witnessData_mp[i].depth);
    EqD->witnessData_mp[i].B->rows = EqD->num_vars;
    EqD->witnessData_mp[i].B->cols = EqD->witnessData_mp[i].depth;
    start = A->cols;

    for (j = 0; j < EqD->witnessData_mp[i].B->rows; j++)
      for (k = 0; k < EqD->witnessData_mp[i].B->cols; k++)
      { // copy to B
        set_mp(&EqD->witnessData_mp[i].B->entry[j][k], &Q->entry[j][k+start]);
      }

    // setup other structures
    EqD->witnessData_mp[i].startPts = (point_mp *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(point_mp));
    EqD->witnessData_mp[i].endPts_in = (point_mp *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(point_mp));
    EqD->witnessData_mp[i].endPts = (point_mp *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(point_mp));
    EqD->witnessData_mp[i].finalTs = (comp_mp *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(comp_mp));
    EqD->witnessData_mp[i].condition_nums = (double *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(double));
    EqD->witnessData_mp[i].endPt_retVals = (int *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(int));
    EqD->witnessData_mp[i].endPt_types = (int *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(int));
    EqD->witnessData_mp[i].higherDim = (int *)bmalloc(EqD->witnessData_mp[i].num_paths * sizeof(int));
    for (j = 0; j < EqD->witnessData_mp[i].num_paths; j++)
    { // initialize the memory
      init_mp(EqD->witnessData_mp[i].finalTs[j]);
      init_point_mp(EqD->witnessData_mp[i].startPts[j], EqD->witnessData_mp[i].depth);
      init_point_mp(EqD->witnessData_mp[i].endPts_in[j], EqD->witnessData_mp[i].depth);
      init_point_mp(EqD->witnessData_mp[i].endPts[j], EqD->num_vars);

      EqD->witnessData_mp[i].startPts[j]->size = EqD->witnessData_mp[i].depth;
      EqD->witnessData_mp[i].endPts_in[j]->size = EqD->witnessData_mp[i].depth;
      EqD->witnessData_mp[i].endPts[j]->size = EqD->num_vars;

      set_zero_mp(EqD->witnessData_mp[i].finalTs[j]);
      EqD->witnessData_mp[i].condition_nums[j] = EqD->witnessData_mp[i].endPt_retVals[j] = EqD->witnessData_mp[i].endPt_types[j] = EqD->witnessData_mp[i].higherDim[j] = 0;
    }

    // setup the start points
    curr_deg = (int *)brealloc(curr_deg, EqD->witnessData_mp[i].depth * sizeof(int));
    for (j = 0; j < EqD->witnessData_mp[i].depth; j++)
      curr_deg[j] = 0;

    change_size_mat_mp(tempMat, EqD->witnessData_mp[i].depth, EqD->witnessData_mp[i].depth);
    change_size_vec_mp(b, EqD->witnessData_mp[i].depth);
    tempMat->rows = tempMat->cols = b->size = EqD->witnessData_mp[i].depth;
    // setup b
    for (k = 0; k < b->size; k++)
    {
      set_one_mp(&b->coord[k]);
    }
    for (j = 0; j < EqD->witnessData_mp[i].num_paths; j++)
    { // setup tempMat
      start = 0;
      for (k = 0; k < tempMat->rows; k++)
      { // copy over the entries
        for (r = 0; r < tempMat->cols; r++)
        {
          set_mp(&tempMat->entry[k][r], EqD->witnessData_mp[i].startSystemCoeff[start + curr_deg[k]][r]);
        }
        // update start
        start += EqD->degrees[EqD->witnessData_mp[i].startFunction + k][0];
      }
      // solve for the start point
      if (matrixSolve_mp(EqD->witnessData_mp[i].startPts[j], tempMat, b))
      { // this should never happen!
        printf("ERROR: Problem solving a random matrix\n");
        bexit(ERROR_OTHER);
      }

      // update curr_deg
      for (k = EqD->witnessData_mp[i].depth - 1; k >= 0; k--)
      { // check to see if we are at the top
        if (curr_deg[k] == EqD->degrees[EqD->witnessData_mp[i].startFunction + k][0] - 1)
        { // set to 0
          curr_deg[k] = 0;
        }
        else
        { // increment this location and exit loop
          curr_deg[k]++;
          k = -10;
        }
      }
    }
  }

  // release memory
  free(str);
  free(curr_deg);

  // clear MP
  mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_mp(b);
  clear_mat_mp(tempMat); clear_mat_mp(A); clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

  return;
}

void setupEqbyEqFirstStage_mp(eqData_t *EqD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first stage of eq-by-eq                      *
\***************************************************************/
{
  // copy over the path information
  EqD->stageData_mp[0].num_paths = EqD->witnessData_mp[0].num_paths;
  EqD->stageData_mp[0].num_sing = EqD->witnessData_mp[0].num_sing;
  EqD->stageData_mp[0].num_nonsing = EqD->witnessData_mp[0].num_nonsing;
  EqD->stageData_mp[0].num_inf = EqD->witnessData_mp[0].num_inf;
  EqD->stageData_mp[0].num_bad = EqD->witnessData_mp[0].num_bad;

  // copy B and p
  EqD->stageData_mp[0].useIntrinsicSlice = 1;
  init_mat_mp(EqD->stageData_mp[0].B, EqD->witnessData_mp[0].B->rows, EqD->witnessData_mp[0].B->cols);
  mat_cp_mp(EqD->stageData_mp[0].B, EqD->witnessData_mp[0].B);
  init_vec_mp(EqD->stageData_mp[0].p, EqD->witnessData_mp[0].p->size);
  vec_cp_mp(EqD->stageData_mp[0].p, EqD->witnessData_mp[0].p);

  // point to the structures inside of witnessData[0]
  EqD->stageData_mp[0].startPts = EqD->witnessData_mp[0].startPts;
  EqD->stageData_mp[0].endPts_in = EqD->witnessData_mp[0].endPts_in;
  EqD->stageData_mp[0].endPts = EqD->witnessData_mp[0].endPts;
  EqD->stageData_mp[0].finalTs = EqD->witnessData_mp[0].finalTs;
  EqD->stageData_mp[0].condition_nums = EqD->witnessData_mp[0].condition_nums;
  EqD->stageData_mp[0].endPt_retVals = EqD->witnessData_mp[0].endPt_retVals;
  EqD->stageData_mp[0].endPt_types = EqD->witnessData_mp[0].endPt_types;
  EqD->stageData_mp[0].higherDim = EqD->witnessData_mp[0].higherDim;

  return;
}

void clearEqbyEqFirstWitnessData_mp(eqData_t *EqD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the first witness data information               *
\***************************************************************/
{
  // set pointers to NULL since stage[0] uses this information
  EqD->witnessData_mp[0].startPts = NULL;
  EqD->witnessData_mp[0].endPts_in = NULL;
  EqD->witnessData_mp[0].endPts = NULL;
  EqD->witnessData_mp[0].finalTs = NULL;
  EqD->witnessData_mp[0].condition_nums = NULL;
  EqD->witnessData_mp[0].endPt_retVals = NULL;
  EqD->witnessData_mp[0].endPt_types = NULL;
  EqD->witnessData_mp[0].higherDim = NULL;

  // clear B and p
  clear_mat_mp(EqD->witnessData_mp[0].B);
  clear_vec_mp(EqD->witnessData_mp[0].p);

  return;
}

void setupEqbyEqNextStage_mp(basic_eval_data_mp *ED, int stage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the 'stage'th stage of eq-by-eq                  *
* Setup a new gamma, patch, move the points to that patch and   *
* then setup B, p and the start points                          *
\***************************************************************/
{
  int i, j, k, depth_x, depth_y, depth_sum, count, indexI, indexJ, num_vars, num_paths, num_digits = prec_to_digits(ED->EqD->curr_precision);
  vec_mp new_vars, tempVec, b;
  mat_mp B_transpose, tempMat, A, Q, R, P;
  mpf_t tol_pivot, tol_sign, largeChange;
  size_t size;
  char *str = NULL;

  // initialize MP
  mpf_init(tol_pivot); mpf_init(tol_sign); mpf_init(largeChange);
  init_vec_mp(new_vars, 0); init_vec_mp(tempVec, 0); init_vec_mp(b, 0);
  init_mat_mp(B_transpose, 0, 0); init_mat_mp(tempMat, 0, 0); init_mat_mp(A, 0, 0);
  init_mat_mp(Q, 0, 0); init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0);

  // error checking
  if (stage <= 0 || stage >= ED->EqD->num_subsystems)
  {
    printf("ERROR: The stage to setup is incorrect!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // setup gamma_mp
  get_comp_rand_mp(ED->EqD->gamma_mp);

  // setup a new patch
  indexI = ED->patch.patchCoeff->rows;
  indexJ = ED->patch.patchCoeff->cols;

  // setup the patch coefficients - only change the ones that are not 0
  for (i = 0; i < indexI; i++)
    for (j = 0; j < indexJ; j++)
      if (!mpfr_zero_p(ED->patch.patchCoeff->entry[i][j].r) && !mpfr_zero_p(ED->patch.patchCoeff->entry[i][j].i))
      {
        get_comp_rand_mp(&ED->patch.patchCoeff->entry[i][j]);
      }

  // move the previous stage points to this new patch
  indexI = ED->EqD->stageData_mp[stage-1].num_paths;
  for (i = 0; i < indexI; i++)
  {
    move_to_patch_mp(ED->EqD->stageData_mp[stage-1].endPts[i], ED->EqD->stageData_mp[stage-1].endPts[i], ED);
  }
  // move the witness data points to this new patch
  indexI = ED->EqD->witnessData_mp[stage].num_paths;
  for (i = 0; i < indexI; i++)
  {
    move_to_patch_mp(ED->EqD->witnessData_mp[stage].endPts[i], ED->EqD->witnessData_mp[stage].endPts[i], ED);
  }

  // setup depth_x & depth_y
  depth_x = ED->EqD->stageData_mp[stage].depth_x;
  depth_y = ED->EqD->stageData_mp[stage].depth_y;
  depth_sum = depth_x + depth_y;

  // calculate the number of variables for this stage
  if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
    num_vars = depth_sum;
  else
    num_vars = 2 * ED->EqD->num_vars;

  // calculate number of paths for this stage
  num_paths = ED->EqD->stageData_mp[stage].num_paths = ED->EqD->stageData_mp[stage - 1].num_nonsing * ED->EqD->witnessData_mp[stage].num_nonsing;

  // initialize the other counts
  ED->EqD->stageData_mp[stage].num_sing = ED->EqD->stageData_mp[stage].num_nonsing = ED->EqD->stageData_mp[stage].num_inf = ED->EqD->stageData_mp[stage].num_bad = 0;

  // allocate memory
  ED->EqD->stageData_mp[stage].startPts = (point_mp *)bmalloc(num_paths * sizeof(point_mp));
  if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
    ED->EqD->stageData_mp[stage].endPts_in = (point_mp *)bmalloc(num_paths * sizeof(point_mp));
  else
    ED->EqD->stageData_mp[stage].endPts_in = NULL;
  ED->EqD->stageData_mp[stage].endPts = (point_mp *)bmalloc(num_paths * sizeof(point_mp));
  ED->EqD->stageData_mp[stage].finalTs = (comp_mp *)bmalloc(num_paths * sizeof(comp_mp));
  ED->EqD->stageData_mp[stage].condition_nums = (double *)bmalloc(num_paths * sizeof(double));
  ED->EqD->stageData_mp[stage].endPt_retVals = (int *)bmalloc(num_paths * sizeof(int));
  ED->EqD->stageData_mp[stage].endPt_types = (int *)bmalloc(num_paths * sizeof(int));
  ED->EqD->stageData_mp[stage].higherDim = (int *)bmalloc(num_paths * sizeof(int));
  for (j = 0; j < num_paths; j++)
  { // initialize the memory
    init_point_mp(ED->EqD->stageData_mp[stage].startPts[j], num_vars);
    ED->EqD->stageData_mp[stage].startPts[j]->size = num_vars;
    if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
    {
      init_point_mp(ED->EqD->stageData_mp[stage].endPts_in[j], num_vars);
      ED->EqD->stageData_mp[stage].endPts_in[j]->size = num_vars; 
    }
    init_point_mp(ED->EqD->stageData_mp[stage].endPts[j], ED->EqD->num_vars);
    ED->EqD->stageData_mp[stage].endPts[j]->size = ED->EqD->num_vars;

    init_mp(ED->EqD->stageData_mp[stage].finalTs[j]);
    set_zero_mp(ED->EqD->stageData_mp[stage].finalTs[j]);
    ED->EqD->stageData_mp[stage].condition_nums[j] = ED->EqD->stageData_mp[stage].endPt_retVals[j] = ED->EqD->stageData_mp[stage].endPt_types[j] = ED->EqD->stageData_mp[stage].higherDim[j] = 0;
  }

  // setup to find the intrinsic slice, if needed
  if (ED->EqD->stageData_mp[stage].useIntrinsicSlice)
  { // setup tol_pivot
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

    // find B & p for the 'fixed linears'
    k = ED->EqD->num_funcs - depth_sum + ED->patch.num_patches;
    change_size_mat_mp(A, k, ED->EqD->num_vars);
    A->rows = k;
    A->cols = ED->EqD->num_vars; // == ED->EqD->num_funcs + ED->patch.num_patches
    // setup tempMat & b
    change_size_mat_mp(tempMat, k, k);
    change_size_vec_mp(b, k);  
    tempMat->rows = tempMat->cols = b->size = k;

    // setup to find the intrinsic slice
    count = 0;
    for (j = depth_sum; j < ED->EqD->num_funcs; j++)
    { // copy over the coefficients
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], ED->EqD->coeff_mp[j][k]);
        if (k < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][k], &A->entry[count][k]);
        }
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients
    for (j = 0; j < ED->patch.num_patches; j++)
    { // set b[count] to 1
      set_one_mp(&b->coord[count]);
      // copy patch
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
        if (k < tempMat->cols)
        { // copy to tempMat
          set_mp(&tempMat->entry[count][k], &A->entry[count][k]);
        }
      }
      count++;
    }

    // setup p - putting 0 in the extra positions
    // by doing it this way, we have a standard way of constructing p from coeff
    init_vec_mp(ED->EqD->stageData_mp[stage].p, 0);
    if (matrixSolve_mp(ED->EqD->stageData_mp[stage].p, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }

    // correctly setup p
    count = tempMat->cols + depth_sum;
    increase_size_vec_mp(ED->EqD->stageData_mp[stage].p, count);
    ED->EqD->stageData_mp[stage].p->size = count;
    for (j = tempMat->cols; j < count; j++)
    {
      set_zero_mp(&ED->EqD->stageData_mp[stage].p->coord[j]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    init_mat_mp(ED->EqD->stageData_mp[stage].B, ED->EqD->num_vars, depth_sum);
    ED->EqD->stageData_mp[stage].B->rows = ED->EqD->num_vars;
    ED->EqD->stageData_mp[stage].B->cols = depth_sum;
    count = A->cols;
    for (j = 0; j < ED->EqD->stageData_mp[stage].B->rows; j++)
      for (k = 0; k < ED->EqD->stageData_mp[stage].B->cols; k++)
      { // copy to B
        set_mp(&ED->EqD->stageData_mp[stage].B->entry[j][k], &Q->entry[j][k+count]);
      }

    // setup B1 & p1
    k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
    change_size_mat_mp(A, k, 2 * ED->EqD->num_vars);
    A->rows = k;
    A->cols = 2 * ED->EqD->num_vars; // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
    // setup tempMat & b
    change_size_mat_mp(tempMat, k, k);
    change_size_vec_mp(b, k);
    tempMat->rows = tempMat->cols = b->size = k;

    count = 0;
    // setup linears in the x variables
    for (j = depth_x; j < ED->EqD->num_funcs; j++)
    { // copy over the coefficients of the linear j
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], ED->EqD->coeff_mp[j][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients for x variables
    for (j = 0; j < ED->patch.num_patches; j++)
    { // set b[count] to 1
      set_one_mp(&b->coord[count]);
      // copy patch
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
      count++;
    }
    // setup linears in the y variables
    for (j = 0; j < depth_x; j++)
    { // copy over the coefficients of linear j
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_zero_mp(&A->entry[count][k]);
        set_mp(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_mp[j][k]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    for (j = depth_sum; j < ED->EqD->num_funcs; j++)
    { // copy over coefficients of linear j
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_zero_mp(&A->entry[count][k]);
        set_mp(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_mp[j][k]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients for y variables
    for (j = 0; j < ED->patch.num_patches; j++)
    { // set b[count] to 1
      set_one_mp(&b->coord[count]);
      // copy patch
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_zero_mp(&A->entry[count][k]);
        set_mp(&A->entry[count][k + ED->EqD->num_vars], &ED->patch.patchCoeff->entry[j][k]);
      }
      count++;
    }

    // setup tempMat to find p using a randomization process
    make_matrix_random_mp(R, A->cols, A->rows, ED->EqD->curr_precision);
    mat_mul_mp(tempMat, A, R);

    // setup p1
    init_vec_mp(ED->EqD->stageData_mp[stage].p1, 0);
    if (matrixSolve_mp(ED->EqD->stageData_mp[stage].p1, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }
    // undo the randomization to find p1
    mul_mat_vec_mp(ED->EqD->stageData_mp[stage].p1, R, ED->EqD->stageData_mp[stage].p1);

    // now, to find B1, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    init_mat_mp(ED->EqD->stageData_mp[stage].B1, 2 * ED->EqD->num_vars, depth_sum);
    ED->EqD->stageData_mp[stage].B1->rows = 2 * ED->EqD->num_vars;
    ED->EqD->stageData_mp[stage].B1->cols = depth_sum;

    count = A->cols;
    for (j = 0; j < ED->EqD->stageData_mp[stage].B1->rows; j++)
      for (k = 0; k < ED->EqD->stageData_mp[stage].B1->cols; k++)
      { // copy to B
        set_mp(&ED->EqD->stageData_mp[stage].B1->entry[j][k], &Q->entry[j][k+count]);
      }

    // setup B0 & p0
    k = 2 * (ED->EqD->num_funcs + ED->patch.num_patches) - depth_sum;
    change_size_mat_mp(A, k, 2 * ED->EqD->num_vars);
    A->rows = k;
    A->cols = 2 * ED->EqD->num_vars; // == 2 * (ED->EqD->num_funcs + ED->patch.num_patches)
    // setup tempMat & b
    change_size_mat_mp(tempMat, k, k);
    change_size_vec_mp(b, k);
    tempMat->rows = tempMat->cols = b->size = k;

    count = 0;
    // setup x_j - y_j for j = depth_x to depth_sum
    for (j = depth_x; j < depth_sum; j++)
    {
      for (k = 0; k < ED->EqD->num_vars; k++)
        if (j == k)
        {
          set_one_mp(&A->entry[count][k]);
          neg_mp(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
        }
        else
        {
          set_zero_mp(&A->entry[count][k]);
          set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
        }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // setup linears in the x variables
    for (j = depth_sum; j < ED->EqD->num_funcs; j++)
    { // copy over the coefficients of the linear j
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], ED->EqD->coeff_mp[j][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients for x variables
    for (j = 0; j < ED->patch.num_patches; j++)
    { // set b[count] to 1
      set_one_mp(&b->coord[count]);
      // copy patch
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_mp(&A->entry[count][k], &ED->patch.patchCoeff->entry[j][k]);
        set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
      }
      count++;
    }
    // setup x_j - y_j for j = 0 to depth_x
    for (j = 0; j < depth_x; j++)
    {
      for (k = 0; k < ED->EqD->num_vars; k++)
        if (j == k)
        {
          set_one_mp(&A->entry[count][k]);
          neg_mp(&A->entry[count][k + ED->EqD->num_vars], &A->entry[count][k]);
        }
        else
        {
          set_zero_mp(&A->entry[count][k]);
          set_zero_mp(&A->entry[count][k + ED->EqD->num_vars]);
        }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    for (j = depth_sum; j < ED->EqD->num_funcs; j++)
    { // copy over coefficients of linear j
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_zero_mp(&A->entry[count][k]);
        set_mp(&A->entry[count][k + ED->EqD->num_vars], ED->EqD->coeff_mp[j][k]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients for y variables
    for (j = 0; j < ED->patch.num_patches; j++)
    { // set b[count] to 1
      set_one_mp(&b->coord[count]);
      // copy patch
      for (k = 0; k < ED->EqD->num_vars; k++)
      {
        set_zero_mp(&A->entry[count][k]);
        set_mp(&A->entry[count][k + ED->EqD->num_vars], &ED->patch.patchCoeff->entry[j][k]);
     }
      count++;
    }

    // setup tempMat to find p using a randomization process
    make_matrix_random_mp(R, A->cols, A->rows, ED->EqD->curr_precision);
    mat_mul_mp(tempMat, A, R);

    // setup p0
    init_vec_mp(ED->EqD->stageData_mp[stage].p0, 0);
    if (matrixSolve_mp(ED->EqD->stageData_mp[stage].p0, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem solving a random matrix\n");
      bexit(ERROR_OTHER);
    }
    // undo the randomization to find p0
    mul_mat_vec_mp(ED->EqD->stageData_mp[stage].p0, R, ED->EqD->stageData_mp[stage].p0);

    // now, to find B0, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    init_mat_mp(ED->EqD->stageData_mp[stage].B0, 2 * ED->EqD->num_vars, depth_sum);
    ED->EqD->stageData_mp[stage].B0->rows = 2 * ED->EqD->num_vars;
    ED->EqD->stageData_mp[stage].B0->cols = depth_sum;

    count = A->cols;
    for (j = 0; j < ED->EqD->stageData_mp[stage].B0->rows; j++)
      for (k = 0; k < ED->EqD->stageData_mp[stage].B0->cols; k++)
      { // copy to B0
        set_mp(&ED->EqD->stageData_mp[stage].B0->entry[j][k], &Q->entry[j][k+count]);
      }

    // now that we have the slice, setup the start points on this slice
    transpose_mp(B_transpose, ED->EqD->stageData_mp[stage].B1);
    // setup the start points
    count = indexI = indexJ = 0;
    num_vars = ED->EqD->num_vars;
    for (i = 0; i < ED->EqD->stageData_mp[stage - 1].num_nonsing; i++)
    { // find the next good path from the previous stage
      while (ED->EqD->stageData_mp[stage - 1].endPt_types[indexI] != MOVE_TO_NEXT)
        indexI++;

      // setup the top of new_vars
      point_cp_mp(new_vars, ED->EqD->stageData_mp[stage-1].endPts[indexI]);
      increase_size_vec_mp(new_vars, 2 * num_vars);
      new_vars->size = 2 * num_vars;

      indexJ = 0;
      for (j = 0; j < ED->EqD->witnessData_mp[stage].num_nonsing; j++)
      { // find the next good path from the next subsystem
        while (ED->EqD->witnessData_mp[stage].endPt_types[indexJ] != MOVE_TO_NEXT)
          indexJ++;

        // setup the bottom of new_vars
        for (k = num_vars; k < 2 * num_vars; k++)
        {
          set_mp(&new_vars->coord[k], &ED->EqD->witnessData_mp[stage].endPts[indexJ]->coord[k - num_vars]);
        }
        new_vars->size = 2 * num_vars;

        // turn new_vars into the intrinsic coordinates
        extrinsicToIntrinsic_mp(ED->EqD->stageData_mp[stage].startPts[count], new_vars, B_transpose, ED->EqD->stageData_mp[stage].p1);

        // increment indexJ
        indexJ++;
        // increment count
        count++;
      }
      // increment indexI
      indexI++;
    }
  }
  else
  { // setup for extrinsic tracking
    count = indexI = indexJ = 0;
    for (i = 0; i < ED->EqD->stageData_mp[stage-1].num_nonsing; i++)
    { // find the next good path from the previous stage
      while (ED->EqD->stageData_mp[stage-1].endPt_types[indexI] != MOVE_TO_NEXT)
        indexI++;

      indexJ = 0;
      for (j = 0; j < ED->EqD->witnessData_mp[stage].num_nonsing; j++)
      { // find the next good path from the next subsystem
        while (ED->EqD->witnessData_mp[stage].endPt_types[indexJ] != MOVE_TO_NEXT)
          indexJ++;

        // setup the start point
        for (k = 0; k < ED->EqD->num_vars; k++)
        {
          set_mp(&ED->EqD->stageData_mp[stage].startPts[count]->coord[k], &ED->EqD->stageData_mp[stage-1].endPts[indexI]->coord[k]);
          set_mp(&ED->EqD->stageData_mp[stage].startPts[count]->coord[k+ED->EqD->num_vars], &ED->EqD->witnessData_mp[stage].endPts[indexJ]->coord[k]);
        }

        // increment indexJ
        indexJ++;
        // increment count
        count++;
      }
      // increment indexI
      indexI++;
    }
  }

  free(str);

  // clear MP
  mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_mp(new_vars); clear_vec_mp(tempVec); clear_vec_mp(b); 
  clear_mat_mp(B_transpose); clear_mat_mp(tempMat); clear_mat_mp(A);
  clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

  return;
}

void clearEqbyEqStageData_mp(eqData_t *EqD, int stage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'stage' stage data information               *
\***************************************************************/
{ 
  int i, num_paths = EqD->stageData_mp[stage].num_paths;
  
  for (i = num_paths - 1; i >= 0; i--)
  {
    if (EqD->stageData_mp[stage].finalTs != NULL)
      clear_mp(EqD->stageData_mp[stage].finalTs[i]);
    if (EqD->stageData_mp[stage].startPts != NULL)
      clear_point_mp(EqD->stageData_mp[stage].startPts[i]);
    if (EqD->stageData_mp[stage].useIntrinsicSlice && EqD->stageData_mp[stage].endPts_in != NULL)
      clear_point_mp(EqD->stageData_mp[stage].endPts_in[i]);
    if (EqD->stageData_mp[stage].endPts != NULL)
      clear_point_mp(EqD->stageData_mp[stage].endPts[i]);
  }
  // clear memory
  free(EqD->stageData_mp[stage].startPts);
  if (EqD->stageData_mp[stage].useIntrinsicSlice);
    free(EqD->stageData_mp[stage].endPts_in);
  free(EqD->stageData_mp[stage].endPts);
  free(EqD->stageData_mp[stage].finalTs);
  free(EqD->stageData_mp[stage].condition_nums);
  free(EqD->stageData_mp[stage].endPt_retVals);
  free(EqD->stageData_mp[stage].endPt_types);
  free(EqD->stageData_mp[stage].higherDim);
  
  if (EqD->stageData_mp[stage].useIntrinsicSlice)
  { // clear B
    clear_mat_mp(EqD->stageData_mp[stage].B);
    // clear p
    clear_vec_mp(EqD->stageData_mp[stage].p);

    // clear B1,B0,p1,p0 if stage != 0
    if (stage != 0)
    { // clear B1
      clear_mat_mp(EqD->stageData_mp[stage].B1);
      // clear p1
      clear_vec_mp(EqD->stageData_mp[stage].p1);
      // clear B0
      clear_mat_mp(EqD->stageData_mp[stage].B0);
      // clear p0
      clear_vec_mp(EqD->stageData_mp[stage].p0);
    }
  }

  return;
}

void clearEqbyEqWitnessData_mp(eqData_t *EqD, int subsystem)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'subsystem' witness data information         *
\***************************************************************/
{
  int i, j, depth = EqD->witnessData_mp[subsystem].depth, num_paths = EqD->witnessData_mp[subsystem].num_paths;

  // clear startSystemCoeff
  for (i = EqD->witnessData_mp[subsystem].num_linears - 1; i >= 0; i--)
  {
    for (j = depth - 1; j >= 0; j--)
    {
      clear_mp(EqD->witnessData_mp[subsystem].startSystemCoeff[i][j]);
    }
    free(EqD->witnessData_mp[subsystem].startSystemCoeff[i]);
  }
  free(EqD->witnessData_mp[subsystem].startSystemCoeff);

  // clear the other memory, if needed
  if (EqD->witnessData_mp[subsystem].startPts != NULL)
  {
    for (i = num_paths - 1; i >= 0; i--)
    {
      clear_mp(EqD->witnessData_mp[subsystem].finalTs[i]);
      clear_point_mp(EqD->witnessData_mp[subsystem].startPts[i]);
      clear_point_mp(EqD->witnessData_mp[subsystem].endPts_in[i]);
      clear_point_mp(EqD->witnessData_mp[subsystem].endPts[i]);
    }
  }
  free(EqD->witnessData_mp[subsystem].startPts);
  free(EqD->witnessData_mp[subsystem].endPts_in);
  free(EqD->witnessData_mp[subsystem].endPts);
  free(EqD->witnessData_mp[subsystem].finalTs);
  free(EqD->witnessData_mp[subsystem].condition_nums);
  free(EqD->witnessData_mp[subsystem].endPt_retVals);
  free(EqD->witnessData_mp[subsystem].endPt_types);
  free(EqD->witnessData_mp[subsystem].higherDim);

  // clear B
  clear_mat_mp(EqD->witnessData_mp[subsystem].B);
  // clear p
  clear_vec_mp(EqD->witnessData_mp[subsystem].p);

  return;
}

/////// CHANGE PRECISION ///////

int change_eqbyeq_eval_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for equation-by-equation              *
\***************************************************************/
{
  int i, k, num_funcs, num_vars;

  // cast ED as BED
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;

  // change the precision on the main parts of BED
  changeBasicEvalPrec_mp(prec, BED);

  // change the precision on the eq-by-eq parts of BED, if needed
  if (BED->EqD->curr_precision != prec)
  { // need to change the random numbers to the correct precision
    BED->EqD->curr_precision = prec;

    // change precision for gamma_mp
    setprec_mp(BED->EqD->gamma_mp, prec);
    mpf_set_q(BED->EqD->gamma_mp->r, BED->EqD->gamma_rat[0]);
    mpf_set_q(BED->EqD->gamma_mp->i, BED->EqD->gamma_rat[1]);

    // change precision for coeff_mp
    num_funcs = BED->EqD->num_funcs;
    num_vars = BED->EqD->num_vars;
    for (i = 0; i < num_funcs; i++)
      for (k = 0; k < num_vars; k++)
      {
        setprec_mp(BED->EqD->coeff_mp[i][k], prec);
        mpf_set_q(BED->EqD->coeff_mp[i][k]->r, BED->EqD->coeff_rat[i][k][0]);
        mpf_set_q(BED->EqD->coeff_mp[i][k]->i, BED->EqD->coeff_rat[i][k][1]);
      }

    if (BED->EqD->increase_witness_prec)
    { // change precision on the witness data
      for (i = 0; i < BED->EqD->num_subsystems; i++)
      {
        changeWitnessDataPrec(&BED->EqD->witnessData_mp[i], prec);
      }
    }
    else
    { // change precision on the current stage
      if (0 <= BED->EqD->curr_stage_num && BED->EqD->curr_stage_num < BED->EqD->num_subsystems)
      {
        changeStageDataPrec(&BED->EqD->stageData_mp[BED->EqD->curr_stage_num], prec);
      }
    }
  }

  return 0;
}

void changeWitnessDataPrec(eqWitnessData_mp *EqD, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for the witness data                  *
\***************************************************************/
{
  int i, j, depth = EqD->depth;
  
  // start system coefficients
  for (i = 0; i < EqD->num_linears; i++)
    for (j = 0; j < depth; j++)
    {
      setprec_mp(EqD->startSystemCoeff[i][j], prec);
      mpf_set_q(EqD->startSystemCoeff[i][j]->r, EqD->startSystemCoeff_rat[i][j][0]);
      mpf_set_q(EqD->startSystemCoeff[i][j]->i, EqD->startSystemCoeff_rat[i][j][1]);
    }

  // null space matrix B
  for (i = 0; i < EqD->B->rows; i++)
    for (j = 0; j < EqD->B->cols; j++)
    {
      setprec_mp(&EqD->B->entry[i][j], prec);
      mpf_set_q(EqD->B->entry[i][j].r, EqD->B_rat[i][j][0]);
      mpf_set_q(EqD->B->entry[i][j].i, EqD->B_rat[i][j][1]);
    }

  // vector p
  for (i = 0; i < EqD->p->size; i++)
  {
    setprec_mp(&EqD->p->coord[i], prec);
    mpf_set_q(EqD->p->coord[i].r, EqD->p_rat[i][0]);
    mpf_set_q(EqD->p->coord[i].i, EqD->p_rat[i][1]);
  }
 
  return;
}

void changeStageDataPrec(eqStageData_mp *EqD, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for the stage data                    *
\***************************************************************/
{
  int i, j;

  if (EqD->useIntrinsicSlice)
  { // null space matrix B
    for (i = 0; i < EqD->B->rows; i++)
      for (j = 0; j < EqD->B->cols; j++)
      {
        setprec_mp(&EqD->B->entry[i][j], prec);
        mpf_set_q(EqD->B->entry[i][j].r, EqD->B_rat[i][j][0]);
        mpf_set_q(EqD->B->entry[i][j].i, EqD->B_rat[i][j][1]);
      }

    // vector p
    for (i = 0; i < EqD->p->size; i++)
    {
      setprec_mp(&EqD->p->coord[i], prec);
      mpf_set_q(EqD->p->coord[i].r, EqD->p_rat[i][0]);
      mpf_set_q(EqD->p->coord[i].i, EqD->p_rat[i][1]);
    }

    // null space matrices B1 & B0
    for (i = 0; i < EqD->B1->rows; i++)
      for (j = 0; j < EqD->B1->cols; j++)
      {
        setprec_mp(&EqD->B1->entry[i][j], prec);
        mpf_set_q(EqD->B1->entry[i][j].r, EqD->B1_rat[i][j][0]);
        mpf_set_q(EqD->B1->entry[i][j].i, EqD->B1_rat[i][j][1]);

        setprec_mp(&EqD->B0->entry[i][j], prec);
        mpf_set_q(EqD->B0->entry[i][j].r, EqD->B0_rat[i][j][0]);
        mpf_set_q(EqD->B0->entry[i][j].i, EqD->B0_rat[i][j][1]);
      }

    // vector p1 & p0
    for (i = 0; i < EqD->p1->size; i++)
    {
      setprec_mp(&EqD->p1->coord[i], prec);
      mpf_set_q(EqD->p1->coord[i].r, EqD->p1_rat[i][0]);
      mpf_set_q(EqD->p1->coord[i].i, EqD->p1_rat[i][1]);

      setprec_mp(&EqD->p0->coord[i], prec);
      mpf_set_q(EqD->p0->coord[i].r, EqD->p0_rat[i][0]);
      mpf_set_q(EqD->p0->coord[i].i, EqD->p0_rat[i][1]);
    }
  }

  return;
}

//// OPENMP SETUP FUNCTIONS /////

void setup_omp_eqbyeq_d(basic_eval_data_d *BED, eqData_t *EqD_in, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copies EqD_in to BED for use in OpenMP tracking        *
*  Since function evaluators can change precision of gamma_mp & *
*  coeff_mp each copy needs to have its own version of these    *
\***************************************************************/
{
  int i, j, k;

  // allocate space
  BED->EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));
  BED->EqD->witnessData_d = (eqWitnessData_d *)bmalloc(EqD_in->num_subsystems * sizeof(eqWitnessData_d));
  BED->EqD->stageData_d = (eqStageData_d *)bmalloc(EqD_in->num_subsystems * sizeof(eqStageData_d));
  if (MPType == 2)
  { // setup for AMP
    BED->EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD_in->num_subsystems * sizeof(eqWitnessData_mp));
    BED->EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD_in->num_subsystems * sizeof(eqStageData_mp));
    BED->BED_mp->EqD = BED->EqD;
  }

  // setup the items in EqD
  BED->EqD->curr_precision = EqD_in->curr_precision;
  BED->EqD->increase_witness_prec = EqD_in->increase_witness_prec;
  BED->EqD->num_funcs = EqD_in->num_funcs;
  BED->EqD->num_subsystems = EqD_in->num_subsystems;
  BED->EqD->num_var_gps = EqD_in->num_var_gps;
  BED->EqD->num_vars = EqD_in->num_vars;
  BED->EqD->curr_stage_num = BED->EqD->curr_path_num = 0;

  // simply point to degrees
  BED->EqD->degrees = EqD_in->degrees;

  // copy over gamma_d
  set_d(BED->EqD->gamma_d, EqD_in->gamma_d);

  // setup InstCount
  BED->EqD->noChanges = EqD_in->noChanges;
  BED->EqD->numSubFuncs = EqD_in->numSubFuncs;
  BED->EqD->startSub = EqD_in->startSub;
  BED->EqD->endSub = EqD_in->endSub;
  BED->EqD->startFunc = EqD_in->startFunc;
  BED->EqD->endFunc = EqD_in->endFunc;
  BED->EqD->startJvsub = EqD_in->startJvsub;
  BED->EqD->endJvsub = EqD_in->endJvsub;
  BED->EqD->startJv = EqD_in->startJv;
  BED->EqD->endJv = EqD_in->endJv;
  BED->EqD->subFuncsBelow = EqD_in->subFuncsBelow;

  // point to coeff_d
  BED->EqD->coeff_d = EqD_in->coeff_d;
  if (MPType == 2)
  { // initialize gamma_mp & gamma_rat
    init_mp2(BED->EqD->gamma_mp, BED->EqD->curr_precision);
    mpq_init(BED->EqD->gamma_rat[0]);
    mpq_init(BED->EqD->gamma_rat[1]);
    // copy over gamma_mp & gamma_rat
    set_mp(BED->EqD->gamma_mp, EqD_in->gamma_mp);
    mpq_set(BED->EqD->gamma_rat[0], EqD_in->gamma_rat[0]);
    mpq_set(BED->EqD->gamma_rat[1], EqD_in->gamma_rat[1]);
    // point to coeff_rat
    BED->EqD->coeff_rat = EqD_in->coeff_rat;
    // setup coeff_mp
    BED->EqD->coeff_mp = (comp_mp **)bmalloc(EqD_in->num_funcs * sizeof(comp_mp *));
    for (i = 0; i < EqD_in->num_funcs; i++)
    {
      BED->EqD->coeff_mp[i] = (comp_mp *)bmalloc(EqD_in->num_vars * sizeof(comp_mp));
      for (k = 0; k < EqD_in->num_vars; k++)
      {
        init_mp2(BED->EqD->coeff_mp[i][k], BED->EqD->curr_precision);
        set_mp(BED->EqD->coeff_mp[i][k], EqD_in->coeff_mp[i][k]);
      }
    }
  }
  // setup the witness data
  for (i = 0; i < EqD_in->num_subsystems; i++)
  {
    BED->EqD->witnessData_d[i].startFunction = EqD_in->witnessData_d[i].startFunction;
    BED->EqD->witnessData_d[i].depth = EqD_in->witnessData_d[i].depth;

    BED->EqD->witnessData_d[i].num_paths = EqD_in->witnessData_d[i].num_paths;
    BED->EqD->witnessData_d[i].num_sing = BED->EqD->witnessData_d[i].num_nonsing = BED->EqD->witnessData_d[i].num_inf = BED->EqD->witnessData_d[i].num_bad = 0;

    BED->EqD->witnessData_d[i].num_linears = EqD_in->witnessData_d[i].num_linears;
    BED->EqD->witnessData_d[i].startSystemCoeff = EqD_in->witnessData_d[i].startSystemCoeff;
    init_mat_d(BED->EqD->witnessData_d[i].B, EqD_in->witnessData_d[i].B->rows, EqD_in->witnessData_d[i].B->cols);
    mat_cp_d(BED->EqD->witnessData_d[i].B, EqD_in->witnessData_d[i].B);
    init_vec_d(BED->EqD->witnessData_d[i].p, EqD_in->witnessData_d[i].p->size);
    vec_cp_d(BED->EqD->witnessData_d[i].p, EqD_in->witnessData_d[i].p);

    BED->EqD->witnessData_d[i].startPts = EqD_in->witnessData_d[i].startPts;
    BED->EqD->witnessData_d[i].endPts_in = EqD_in->witnessData_d[i].endPts_in;
    BED->EqD->witnessData_d[i].endPts = EqD_in->witnessData_d[i].endPts;

    BED->EqD->witnessData_d[i].finalTs = EqD_in->witnessData_d[i].finalTs;
    BED->EqD->witnessData_d[i].condition_nums = EqD_in->witnessData_d[i].condition_nums;

    BED->EqD->witnessData_d[i].endPt_retVals = EqD_in->witnessData_d[i].endPt_retVals;
    BED->EqD->witnessData_d[i].endPt_types = EqD_in->witnessData_d[i].endPt_types;
    BED->EqD->witnessData_d[i].higherDim = EqD_in->witnessData_d[i].higherDim;

    if (MPType == 2)
    { // setup witnessData_mp[i]
      BED->EqD->witnessData_mp[i].startFunction = EqD_in->witnessData_mp[i].startFunction;
      BED->EqD->witnessData_mp[i].depth = EqD_in->witnessData_mp[i].depth;

      BED->EqD->witnessData_mp[i].num_paths = EqD_in->witnessData_mp[i].num_paths;
      BED->EqD->witnessData_mp[i].num_sing = BED->EqD->witnessData_mp[i].num_nonsing = BED->EqD->witnessData_mp[i].num_inf = BED->EqD->witnessData_mp[i].num_bad = 0;
      BED->EqD->witnessData_mp[i].startPts = EqD_in->witnessData_mp[i].startPts;
      BED->EqD->witnessData_mp[i].endPts_in = EqD_in->witnessData_mp[i].endPts_in;
      BED->EqD->witnessData_mp[i].endPts = EqD_in->witnessData_mp[i].endPts;

      BED->EqD->witnessData_mp[i].finalTs = EqD_in->witnessData_mp[i].finalTs;
      BED->EqD->witnessData_mp[i].condition_nums = EqD_in->witnessData_mp[i].condition_nums;

      BED->EqD->witnessData_mp[i].endPt_retVals = EqD_in->witnessData_mp[i].endPt_retVals;
      BED->EqD->witnessData_mp[i].endPt_types = EqD_in->witnessData_mp[i].endPt_types;
      BED->EqD->witnessData_mp[i].higherDim = EqD_in->witnessData_mp[i].higherDim;

      // setup the coefficients, B and p
      BED->EqD->witnessData_mp[i].num_linears = EqD_in->witnessData_mp[i].num_linears;
      // point to coeff_rat
      BED->EqD->witnessData_mp[i].startSystemCoeff_rat = EqD_in->witnessData_mp[i].startSystemCoeff_rat;
      // setup coeff
      BED->EqD->witnessData_mp[i].startSystemCoeff = (comp_mp **)bmalloc(EqD_in->witnessData_mp[i].num_linears * sizeof(comp_mp *));
      for (j = 0; j < EqD_in->witnessData_mp[i].num_linears; j++)
      {
        BED->EqD->witnessData_mp[i].startSystemCoeff[j] = (comp_mp *)bmalloc(EqD_in->witnessData_mp[i].depth * sizeof(comp_mp));
        for (k = 0; k < EqD_in->witnessData_mp[i].depth; k++)
        {
          init_mp2(BED->EqD->witnessData_mp[i].startSystemCoeff[j][k], BED->EqD->curr_precision);
          set_mp(BED->EqD->witnessData_mp[i].startSystemCoeff[j][k], EqD_in->witnessData_mp[i].startSystemCoeff[j][k]);
        }
      }
      // point to B_rat
      BED->EqD->witnessData_mp[i].B_rat = EqD_in->witnessData_mp[i].B_rat;
      // setup B
      init_mat_mp2(BED->EqD->witnessData_mp[i].B, EqD_in->witnessData_mp[i].B->rows, EqD_in->witnessData_mp[i].B->cols, BED->EqD->curr_precision);
      mat_cp_mp(BED->EqD->witnessData_mp[i].B, EqD_in->witnessData_mp[i].B);
      // point to p_rat
      BED->EqD->witnessData_mp[i].p_rat = EqD_in->witnessData_mp[i].p_rat;
      // setup p
      init_vec_mp2(BED->EqD->witnessData_mp[i].p, EqD_in->witnessData_mp[i].p->size, BED->EqD->curr_precision);
      vec_cp_mp(BED->EqD->witnessData_mp[i].p, EqD_in->witnessData_mp[i].p);
    }
  }
  // setup the stage data
  for (i = 0; i < EqD_in->num_subsystems; i++)
  {
    BED->EqD->stageData_d[i].depth_x = EqD_in->stageData_d[i].depth_x;
    BED->EqD->stageData_d[i].depth_y = EqD_in->stageData_d[i].depth_y;
    BED->EqD->stageData_d[i].useIntrinsicSlice = EqD_in->stageData_d[i].useIntrinsicSlice;

    if (MPType == 2)
    {
      BED->EqD->stageData_mp[i].depth_x = EqD_in->stageData_mp[i].depth_x;
      BED->EqD->stageData_mp[i].depth_y = EqD_in->stageData_mp[i].depth_y;
      BED->EqD->stageData_mp[i].useIntrinsicSlice = EqD_in->stageData_mp[i].useIntrinsicSlice;
    }
  }

  return;
}

void setup_omp_eqbyeq_mp(basic_eval_data_mp *BED, eqData_t *EqD_in)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copies EqD_in to BED for use in OpenMP tracking        *
* Since fixed precision, can simply point to coeff_mp           *
\***************************************************************/
{
  int i;
  
  // allocate space
  BED->EqD = (eqData_t *)bmalloc(1 * sizeof(eqData_t));
  BED->EqD->witnessData_mp = (eqWitnessData_mp *)bmalloc(EqD_in->num_subsystems * sizeof(eqWitnessData_mp));
  BED->EqD->stageData_mp = (eqStageData_mp *)bmalloc(EqD_in->num_subsystems * sizeof(eqStageData_mp));

  // setup the items in EqD
  BED->EqD->curr_precision = EqD_in->curr_precision;
  BED->EqD->increase_witness_prec = EqD_in->increase_witness_prec;
  BED->EqD->num_funcs = EqD_in->num_funcs;
  BED->EqD->num_subsystems = EqD_in->num_subsystems;
  BED->EqD->num_var_gps = EqD_in->num_var_gps;
  BED->EqD->num_vars = EqD_in->num_vars;
  BED->EqD->curr_stage_num = BED->EqD->curr_path_num = 0;

  // simply point to degrees
  BED->EqD->degrees = EqD_in->degrees;

  // setup gamma_mp
  init_mp(BED->EqD->gamma_mp);
  set_mp(BED->EqD->gamma_mp, EqD_in->gamma_mp);

  // setup InstCount
  BED->EqD->noChanges = EqD_in->noChanges;
  BED->EqD->numSubFuncs = EqD_in->numSubFuncs;
  BED->EqD->startSub = EqD_in->startSub;
  BED->EqD->endSub = EqD_in->endSub;
  BED->EqD->startFunc = EqD_in->startFunc;
  BED->EqD->endFunc = EqD_in->endFunc;
  BED->EqD->startJvsub = EqD_in->startJvsub;
  BED->EqD->endJvsub = EqD_in->endJvsub;
  BED->EqD->startJv = EqD_in->startJv;
  BED->EqD->endJv = EqD_in->endJv;
  BED->EqD->subFuncsBelow = EqD_in->subFuncsBelow;

  // setup coeff_mp
  BED->EqD->coeff_mp = EqD_in->coeff_mp;

  // setup the witness data
  for (i = 0; i < EqD_in->num_subsystems; i++)
  { // setup witnessData_mp[i]
    BED->EqD->witnessData_mp[i].startFunction = EqD_in->witnessData_mp[i].startFunction;
    BED->EqD->witnessData_mp[i].depth = EqD_in->witnessData_mp[i].depth;

    BED->EqD->witnessData_mp[i].num_paths = EqD_in->witnessData_mp[i].num_paths;
    BED->EqD->witnessData_mp[i].num_sing = BED->EqD->witnessData_mp[i].num_nonsing = BED->EqD->witnessData_mp[i].num_inf = BED->EqD->witnessData_mp[i].num_bad = 0;

    BED->EqD->witnessData_mp[i].num_linears = EqD_in->witnessData_mp[i].num_linears;
    BED->EqD->witnessData_mp[i].startSystemCoeff = EqD_in->witnessData_mp[i].startSystemCoeff;

    // setup B
    init_mat_mp(BED->EqD->witnessData_mp[i].B, EqD_in->witnessData_mp[i].B->rows, EqD_in->witnessData_mp[i].B->cols);
    mat_cp_mp(BED->EqD->witnessData_mp[i].B, EqD_in->witnessData_mp[i].B);
    // setup p
    init_vec_mp(BED->EqD->witnessData_mp[i].p, EqD_in->witnessData_mp[i].p->size);
    vec_cp_mp(BED->EqD->witnessData_mp[i].p, EqD_in->witnessData_mp[i].p);

    BED->EqD->witnessData_mp[i].startPts = EqD_in->witnessData_mp[i].startPts;
    BED->EqD->witnessData_mp[i].endPts_in = EqD_in->witnessData_mp[i].endPts_in;
    BED->EqD->witnessData_mp[i].endPts = EqD_in->witnessData_mp[i].endPts;

    BED->EqD->witnessData_mp[i].finalTs = EqD_in->witnessData_mp[i].finalTs;
    BED->EqD->witnessData_mp[i].condition_nums = EqD_in->witnessData_mp[i].condition_nums;

    BED->EqD->witnessData_mp[i].endPt_retVals = EqD_in->witnessData_mp[i].endPt_retVals;
    BED->EqD->witnessData_mp[i].endPt_types = EqD_in->witnessData_mp[i].endPt_types;
    BED->EqD->witnessData_mp[i].higherDim = EqD_in->witnessData_mp[i].higherDim;
  }
  // setup the stage data
  for (i = 0; i < EqD_in->num_subsystems; i++)
  {
    BED->EqD->stageData_mp[i].depth_x = EqD_in->stageData_mp[i].depth_x;
    BED->EqD->stageData_mp[i].depth_y = EqD_in->stageData_mp[i].depth_y;
    BED->EqD->stageData_mp[i].useIntrinsicSlice = EqD_in->stageData_mp[i].useIntrinsicSlice;
  }

  return;
}

void setup_omp_eqbyeq_first_stage_d(int max_threads, basic_eval_data_d ED_copy[], eqData_t *EqD, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first stage of eq-by-eq for use with OpenMP  *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // do not worry about increasing the precision for the witness data structures
      ED_copy[i].EqD->increase_witness_prec = 0;

      // copy over the path information
      ED_copy[i].EqD->stageData_d[0].num_paths = EqD->stageData_d[0].num_paths;
      ED_copy[i].EqD->stageData_d[0].num_sing = EqD->stageData_d[0].num_sing;
      ED_copy[i].EqD->stageData_d[0].num_nonsing = EqD->stageData_d[0].num_nonsing;
      ED_copy[i].EqD->stageData_d[0].num_inf = EqD->stageData_d[0].num_inf;
      ED_copy[i].EqD->stageData_d[0].num_bad = EqD->stageData_d[0].num_bad;

      // copy B and p
      ED_copy[i].EqD->stageData_d[0].useIntrinsicSlice = EqD->stageData_d[0].useIntrinsicSlice;
      init_mat_d(ED_copy[i].EqD->stageData_d[0].B, EqD->stageData_d[0].B->rows, EqD->stageData_d[0].B->cols);
      mat_cp_d(ED_copy[i].EqD->stageData_d[0].B, EqD->stageData_d[0].B);
      init_vec_d(ED_copy[i].EqD->stageData_d[0].p, EqD->stageData_d[0].p->size);
      vec_cp_d(ED_copy[i].EqD->stageData_d[0].p, EqD->stageData_d[0].p);

      // point to the structures inside of ED
      ED_copy[i].EqD->stageData_d[0].startPts = EqD->stageData_d[0].startPts;
      ED_copy[i].EqD->stageData_d[0].endPts_in = EqD->stageData_d[0].endPts_in;
      ED_copy[i].EqD->stageData_d[0].endPts = EqD->stageData_d[0].endPts;
      ED_copy[i].EqD->stageData_d[0].finalTs = EqD->stageData_d[0].finalTs;
      ED_copy[i].EqD->stageData_d[0].condition_nums = EqD->stageData_d[0].condition_nums;
      ED_copy[i].EqD->stageData_d[0].endPt_retVals = EqD->stageData_d[0].endPt_retVals;
      ED_copy[i].EqD->stageData_d[0].endPt_types = EqD->stageData_d[0].endPt_types;
      ED_copy[i].EqD->stageData_d[0].higherDim = EqD->stageData_d[0].higherDim;

      if (MPType == 2)
      { 
        ED_copy[i].EqD->stageData_mp[0].useIntrinsicSlice = EqD->stageData_mp[0].useIntrinsicSlice;
        if (EqD->stageData_d[0].useIntrinsicSlice)
        { // copy B and p
          init_mat_mp2(ED_copy[i].EqD->stageData_mp[0].B, EqD->stageData_mp[0].B->rows, EqD->stageData_mp[0].B->cols, ED_copy[i].EqD->curr_precision);
          mat_cp_mp(ED_copy[i].EqD->stageData_mp[0].B, EqD->stageData_mp[0].B);
          init_vec_mp2(ED_copy[i].EqD->stageData_mp[0].p, EqD->stageData_mp[0].p->size, ED_copy[i].EqD->curr_precision);
          vec_cp_mp(ED_copy[i].EqD->stageData_mp[0].p, EqD->stageData_mp[0].p);

          // point to B_rat & p_rat
          ED_copy[i].EqD->stageData_mp[0].B_rat = EqD->stageData_mp[0].B_rat;
          ED_copy[i].EqD->stageData_mp[0].p_rat = EqD->stageData_mp[0].p_rat;

          // NULL out B1_rat, p1_rat, B0_rat, p0_rat
          ED_copy[i].EqD->stageData_mp[0].p1_rat = ED_copy[i].EqD->stageData_mp[0].p0_rat = NULL;
          ED_copy[i].EqD->stageData_mp[0].B1_rat = ED_copy[i].EqD->stageData_mp[0].B0_rat = NULL;
        }
        else
        { // set to B_rat & p_rat to NULL
          ED_copy[i].EqD->stageData_mp[0].B_rat = NULL;
          ED_copy[i].EqD->stageData_mp[0].p_rat = NULL;
        }

        // set all other pointers to NULL
        ED_copy[i].EqD->stageData_mp[0].startPts = NULL;
        ED_copy[i].EqD->stageData_mp[0].endPts_in = NULL;
        ED_copy[i].EqD->stageData_mp[0].endPts = NULL;
        ED_copy[i].EqD->stageData_mp[0].finalTs = NULL;
        ED_copy[i].EqD->stageData_mp[0].condition_nums = NULL;
        ED_copy[i].EqD->stageData_mp[0].endPt_retVals = NULL;
        ED_copy[i].EqD->stageData_mp[0].endPt_types = NULL;
        ED_copy[i].EqD->stageData_mp[0].higherDim = NULL;
      }
    }
  }

  return;
}

void setup_omp_eqbyeq_first_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], eqData_t *EqD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first stage of eq-by-eq for use with OpenMP  *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // do not worry about increasing the precision for the witness data structures
      ED_copy[i].EqD->increase_witness_prec = 0;

      // copy over the path information
      ED_copy[i].EqD->stageData_mp[0].num_paths = EqD->stageData_mp[0].num_paths;
      ED_copy[i].EqD->stageData_mp[0].num_sing = EqD->stageData_mp[0].num_sing;
      ED_copy[i].EqD->stageData_mp[0].num_nonsing = EqD->stageData_mp[0].num_nonsing;
      ED_copy[i].EqD->stageData_mp[0].num_inf = EqD->stageData_mp[0].num_inf;
      ED_copy[i].EqD->stageData_mp[0].num_bad = EqD->stageData_mp[0].num_bad;

      ED_copy[i].EqD->stageData_mp[0].useIntrinsicSlice = EqD->stageData_mp[0].useIntrinsicSlice;
      if (EqD->stageData_mp[0].useIntrinsicSlice)
      { // copy B and p
        init_mat_mp(ED_copy[i].EqD->stageData_mp[0].B, EqD->stageData_mp[0].B->rows, EqD->stageData_mp[0].B->cols);
        mat_cp_mp(ED_copy[i].EqD->stageData_mp[0].B, EqD->stageData_mp[0].B);
        init_vec_mp(ED_copy[i].EqD->stageData_mp[0].p, EqD->stageData_mp[0].p->size);
        vec_cp_mp(ED_copy[i].EqD->stageData_mp[0].p, EqD->stageData_mp[0].p);
      }

      // point to the structures inside of ED
      ED_copy[i].EqD->stageData_mp[0].startPts = EqD->stageData_mp[0].startPts;
      ED_copy[i].EqD->stageData_mp[0].endPts_in = EqD->stageData_mp[0].endPts_in;
      ED_copy[i].EqD->stageData_mp[0].endPts = EqD->stageData_mp[0].endPts;
      ED_copy[i].EqD->stageData_mp[0].finalTs = EqD->stageData_mp[0].finalTs;
      ED_copy[i].EqD->stageData_mp[0].condition_nums = EqD->stageData_mp[0].condition_nums;
      ED_copy[i].EqD->stageData_mp[0].endPt_retVals = EqD->stageData_mp[0].endPt_retVals;
      ED_copy[i].EqD->stageData_mp[0].endPt_types = EqD->stageData_mp[0].endPt_types;
      ED_copy[i].EqD->stageData_mp[0].higherDim = EqD->stageData_mp[0].higherDim;
    }
  }

  return;
}

void setup_omp_eqbyeq_next_stage_d(int max_threads, basic_eval_data_d ED_copy[], eqData_t *EqD, int stage, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the next stage of eq-by-eq for use with OpenMP   *
* Copy over gamma since that was changed                        *
* Even though the patch was changed, the copies do not need to  *
* know this since all the data is contained in B & p            *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // copy over the path information
      ED_copy[i].EqD->stageData_d[stage].num_paths = EqD->stageData_d[stage].num_paths;
      ED_copy[i].EqD->stageData_d[stage].num_sing = EqD->stageData_d[stage].num_sing;
      ED_copy[i].EqD->stageData_d[stage].num_nonsing = EqD->stageData_d[stage].num_nonsing;
      ED_copy[i].EqD->stageData_d[stage].num_inf = EqD->stageData_d[stage].num_inf;
      ED_copy[i].EqD->stageData_d[stage].num_bad = EqD->stageData_d[stage].num_bad;

      // copy gamma_d
      set_d(ED_copy[i].EqD->gamma_d, EqD->gamma_d);

      // copy B and p
      ED_copy[i].EqD->stageData_d[stage].useIntrinsicSlice = EqD->stageData_d[stage].useIntrinsicSlice;
      if (ED_copy[i].EqD->stageData_d[stage].useIntrinsicSlice)
      {
        init_mat_d(ED_copy[i].EqD->stageData_d[stage].B, EqD->stageData_d[stage].B->rows, EqD->stageData_d[stage].B->cols);
        mat_cp_d(ED_copy[i].EqD->stageData_d[stage].B, EqD->stageData_d[stage].B);
        init_vec_d(ED_copy[i].EqD->stageData_d[stage].p, EqD->stageData_d[stage].p->size);
        vec_cp_d(ED_copy[i].EqD->stageData_d[stage].p, EqD->stageData_d[stage].p);

        init_mat_d(ED_copy[i].EqD->stageData_d[stage].B1, EqD->stageData_d[stage].B1->rows, EqD->stageData_d[stage].B1->cols);
        mat_cp_d(ED_copy[i].EqD->stageData_d[stage].B1, EqD->stageData_d[stage].B1);
        init_vec_d(ED_copy[i].EqD->stageData_d[stage].p1, EqD->stageData_d[stage].p1->size);
        vec_cp_d(ED_copy[i].EqD->stageData_d[stage].p1, EqD->stageData_d[stage].p1);
        init_mat_d(ED_copy[i].EqD->stageData_d[stage].B0, EqD->stageData_d[stage].B0->rows, EqD->stageData_d[stage].B0->cols);
        mat_cp_d(ED_copy[i].EqD->stageData_d[stage].B0, EqD->stageData_d[stage].B0);
        init_vec_d(ED_copy[i].EqD->stageData_d[stage].p0, EqD->stageData_d[stage].p0->size);
        vec_cp_d(ED_copy[i].EqD->stageData_d[stage].p0, EqD->stageData_d[stage].p0);
      }

      // point to the structures inside of ED
      ED_copy[i].EqD->stageData_d[stage].startPts = EqD->stageData_d[stage].startPts;
      ED_copy[i].EqD->stageData_d[stage].endPts_in = EqD->stageData_d[stage].endPts_in;
      ED_copy[i].EqD->stageData_d[stage].endPts = EqD->stageData_d[stage].endPts;
      ED_copy[i].EqD->stageData_d[stage].finalTs = EqD->stageData_d[stage].finalTs;
      ED_copy[i].EqD->stageData_d[stage].condition_nums = EqD->stageData_d[stage].condition_nums;

      ED_copy[i].EqD->stageData_d[stage].endPt_retVals = EqD->stageData_d[stage].endPt_retVals;
      ED_copy[i].EqD->stageData_d[stage].endPt_types = EqD->stageData_d[stage].endPt_types;
      ED_copy[i].EqD->stageData_d[stage].higherDim = EqD->stageData_d[stage].higherDim;

      if (MPType == 2)
      { 
        ED_copy[i].EqD->stageData_mp[stage].useIntrinsicSlice = EqD->stageData_mp[stage].useIntrinsicSlice;
        if (ED_copy[i].EqD->stageData_mp[stage].useIntrinsicSlice)
        { // copy B and p
          init_mat_mp2(ED_copy[i].EqD->stageData_mp[stage].B, EqD->stageData_mp[stage].B->rows, EqD->stageData_mp[stage].B->cols, ED_copy[i].EqD->curr_precision);
          mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B, EqD->stageData_mp[stage].B);
          init_vec_mp2(ED_copy[i].EqD->stageData_mp[stage].p, EqD->stageData_mp[stage].p->size, ED_copy[i].EqD->curr_precision);
          vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p, EqD->stageData_mp[stage].p);

          // point to B_rat & p_rat
          ED_copy[i].EqD->stageData_mp[stage].B_rat = EqD->stageData_mp[stage].B_rat;
          ED_copy[i].EqD->stageData_mp[stage].p_rat = EqD->stageData_mp[stage].p_rat;

          // copy B1, B0, p1 & p0
          init_mat_mp2(ED_copy[i].EqD->stageData_mp[stage].B1, EqD->stageData_mp[stage].B1->rows, EqD->stageData_mp[stage].B1->cols, ED_copy[i].EqD->curr_precision);
          mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B1, EqD->stageData_mp[stage].B1);
          init_vec_mp2(ED_copy[i].EqD->stageData_mp[stage].p1, EqD->stageData_mp[stage].p1->size, ED_copy[i].EqD->curr_precision);
          vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p1, EqD->stageData_mp[stage].p1);

          init_mat_mp2(ED_copy[i].EqD->stageData_mp[stage].B0, EqD->stageData_mp[stage].B0->rows, EqD->stageData_mp[stage].B0->cols, ED_copy[i].EqD->curr_precision);
          mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B0, EqD->stageData_mp[stage].B0);
          init_vec_mp2(ED_copy[i].EqD->stageData_mp[stage].p0, EqD->stageData_mp[stage].p0->size, ED_copy[i].EqD->curr_precision);
          vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p0, EqD->stageData_mp[stage].p0);

          // point to B1_rat, B0_rat, p1_rat & p0_rat
          ED_copy[i].EqD->stageData_mp[stage].B1_rat = EqD->stageData_mp[stage].B1_rat;
          ED_copy[i].EqD->stageData_mp[stage].p1_rat = EqD->stageData_mp[stage].p1_rat;
          ED_copy[i].EqD->stageData_mp[stage].B0_rat = EqD->stageData_mp[stage].B0_rat;
          ED_copy[i].EqD->stageData_mp[stage].p0_rat = EqD->stageData_mp[stage].p0_rat;
        }
        else
        { // set to B_rat & p_rat to NULL
          ED_copy[i].EqD->stageData_mp[stage].B_rat = NULL;
          ED_copy[i].EqD->stageData_mp[stage].p_rat = NULL;
        }

        // copy gamma_mp & gamma_rat
        set_mp(ED_copy[i].EqD->gamma_mp, EqD->gamma_mp);
        mpq_set(ED_copy[i].EqD->gamma_rat[0], EqD->gamma_rat[0]);
        mpq_set(ED_copy[i].EqD->gamma_rat[1], EqD->gamma_rat[1]);
     
        // set all other pointers to NULL
        ED_copy[i].EqD->stageData_mp[stage].startPts = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPts_in = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPts = NULL;
        ED_copy[i].EqD->stageData_mp[stage].finalTs = NULL;
        ED_copy[i].EqD->stageData_mp[stage].condition_nums = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPt_retVals = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPt_types = NULL;
        ED_copy[i].EqD->stageData_mp[stage].higherDim = NULL;
      }
    }
  }

  return;
}

void setup_omp_eqbyeq_next_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], eqData_t *EqD, int stage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the next stage of eq-by-eq for use with OpenMP   *
* Copy over gamma since that was changed                        *
* Even though the patch was changed, the copies do not need to  *
* know this since all the data is contained in B & p            *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // copy over the path information
      ED_copy[i].EqD->stageData_mp[stage].num_paths = EqD->stageData_mp[stage].num_paths;
      ED_copy[i].EqD->stageData_mp[stage].num_sing = EqD->stageData_mp[stage].num_sing;
      ED_copy[i].EqD->stageData_mp[stage].num_nonsing = EqD->stageData_mp[stage].num_nonsing;
      ED_copy[i].EqD->stageData_mp[stage].num_inf = EqD->stageData_mp[stage].num_inf;
      ED_copy[i].EqD->stageData_mp[stage].num_bad = EqD->stageData_mp[stage].num_bad;

      // copy gamma_mp & gamma_rat
      set_mp(ED_copy[i].EqD->gamma_mp, EqD->gamma_mp);

      ED_copy[i].EqD->stageData_mp[stage].useIntrinsicSlice = EqD->stageData_mp[stage].useIntrinsicSlice;
      if (EqD->stageData_mp[stage].useIntrinsicSlice)
      { // copy B and p
        init_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B, EqD->stageData_mp[stage].B->rows, EqD->stageData_mp[stage].B->cols);
        mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B, EqD->stageData_mp[stage].B);
        init_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p, EqD->stageData_mp[stage].p->size);
        vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p, EqD->stageData_mp[stage].p);

        // copy B1, B0, p1 and p0
        init_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B1, EqD->stageData_mp[stage].B1->rows, EqD->stageData_mp[stage].B1->cols);
        mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B1, EqD->stageData_mp[stage].B1);
        init_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p1, EqD->stageData_mp[stage].p1->size);
        vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p1, EqD->stageData_mp[stage].p1);

        init_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B0, EqD->stageData_mp[stage].B0->rows, EqD->stageData_mp[stage].B0->cols);
        mat_cp_mp(ED_copy[i].EqD->stageData_mp[stage].B0, EqD->stageData_mp[stage].B0);
        init_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p0, EqD->stageData_mp[stage].p0->size);
        vec_cp_mp(ED_copy[i].EqD->stageData_mp[stage].p0, EqD->stageData_mp[stage].p0);
      }

      // point to the structures inside of ED
      ED_copy[i].EqD->stageData_mp[stage].startPts = EqD->stageData_mp[stage].startPts;
      ED_copy[i].EqD->stageData_mp[stage].endPts_in = EqD->stageData_mp[stage].endPts_in;
      ED_copy[i].EqD->stageData_mp[stage].endPts = EqD->stageData_mp[stage].endPts;
      ED_copy[i].EqD->stageData_mp[stage].finalTs = EqD->stageData_mp[stage].finalTs;
      ED_copy[i].EqD->stageData_mp[stage].condition_nums = EqD->stageData_mp[stage].condition_nums;

      ED_copy[i].EqD->stageData_mp[stage].endPt_retVals = EqD->stageData_mp[stage].endPt_retVals;
      ED_copy[i].EqD->stageData_mp[stage].endPt_types = EqD->stageData_mp[stage].endPt_types;
      ED_copy[i].EqD->stageData_mp[stage].higherDim = EqD->stageData_mp[stage].higherDim;
    }
  }

  return;
}

void clear_omp_eqbyeq_stage_d(int max_threads, basic_eval_data_d ED_copy[], int stage, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'stage' stage data information               *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // set the pointers to NULL
      ED_copy[i].EqD->stageData_d[stage].startPts = NULL;
      ED_copy[i].EqD->stageData_d[stage].endPts_in = NULL;
      ED_copy[i].EqD->stageData_d[stage].endPts = NULL;
      ED_copy[i].EqD->stageData_d[stage].finalTs = NULL;
      ED_copy[i].EqD->stageData_d[stage].condition_nums = NULL;
      ED_copy[i].EqD->stageData_d[stage].endPt_retVals = NULL;
      ED_copy[i].EqD->stageData_d[stage].endPt_types = NULL;
      ED_copy[i].EqD->stageData_d[stage].higherDim = NULL;

      if (MPType == 2)
      { 
        if (ED_copy[i].EqD->stageData_mp[stage].useIntrinsicSlice)
        { // clear B and p
          clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B);
          clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p);

          // clear B1, B0, p1 & p0, if needed
          if (ED_copy[i].EqD->stageData_mp[stage].B1_rat != NULL)
          { // clear B1
            clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B1);
          }
          if (ED_copy[i].EqD->stageData_mp[stage].B0_rat != NULL)
          { // clear B0
            clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B0);
          }
          if (ED_copy[i].EqD->stageData_mp[stage].p1_rat != NULL)
          { // clear p1
            clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p1);
          }
          if (ED_copy[i].EqD->stageData_mp[stage].p0_rat != NULL)
          { // clear p0
            clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p0);
          }
        }
        // set all other pointers to NULL
        ED_copy[i].EqD->stageData_mp[stage].B_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].p_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].B1_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].p1_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].B0_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].p0_rat = NULL;
        ED_copy[i].EqD->stageData_mp[stage].startPts = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPts_in = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPts = NULL;
        ED_copy[i].EqD->stageData_mp[stage].finalTs = NULL;
        ED_copy[i].EqD->stageData_mp[stage].condition_nums = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPt_retVals = NULL;
        ED_copy[i].EqD->stageData_mp[stage].endPt_types = NULL;
        ED_copy[i].EqD->stageData_mp[stage].higherDim = NULL;
      }
    }
  }

  return;
}

void clear_omp_eqbyeq_stage_mp(int max_threads, basic_eval_data_mp ED_copy[], int stage)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'stage' stage data information               *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // set the pointers to NULL
      ED_copy[i].EqD->stageData_mp[stage].startPts = NULL;
      ED_copy[i].EqD->stageData_mp[stage].endPts_in = NULL;
      ED_copy[i].EqD->stageData_mp[stage].endPts = NULL;
      ED_copy[i].EqD->stageData_mp[stage].finalTs = NULL;
      ED_copy[i].EqD->stageData_mp[stage].condition_nums = NULL;
      ED_copy[i].EqD->stageData_mp[stage].endPt_retVals = NULL;
      ED_copy[i].EqD->stageData_mp[stage].endPt_types = NULL;
      ED_copy[i].EqD->stageData_mp[stage].higherDim = NULL;

      if (ED_copy[i].EqD->stageData_mp[stage].useIntrinsicSlice)
      { // clear B and p
        clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B);
        clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p);

        if (stage != 0)
        { // clear B1, B0, p1 & p0
          clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B1);
          clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p1);
          clear_mat_mp(ED_copy[i].EqD->stageData_mp[stage].B0);
          clear_vec_mp(ED_copy[i].EqD->stageData_mp[stage].p0);
        }
      }
    }
  }

  return;
}

void clear_omp_eqbyeq_witness_d(int max_threads, basic_eval_data_d ED_copy[], int subsystem, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'subsystem' witness data information         *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // set the pointers to NULL
      ED_copy[i].EqD->witnessData_d[subsystem].startPts = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].endPts_in = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].endPts = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].finalTs = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].condition_nums = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].endPt_retVals = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].endPt_types = NULL;
      ED_copy[i].EqD->witnessData_d[subsystem].higherDim = NULL;

      if (MPType == 2)
      { // clear B and p and set all other pointer to NULL
        clear_mat_mp(ED_copy[i].EqD->witnessData_mp[subsystem].B);
        clear_vec_mp(ED_copy[i].EqD->witnessData_mp[subsystem].p);

        ED_copy[i].EqD->witnessData_mp[subsystem].B_rat = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].p_rat = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].startPts = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].endPts_in = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].endPts = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].finalTs = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].condition_nums = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].endPt_retVals = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].endPt_types = NULL;
        ED_copy[i].EqD->witnessData_mp[subsystem].higherDim = NULL;
      }
    }
  }
  
  return;
}

void clear_omp_eqbyeq_witness_mp(int max_threads, basic_eval_data_mp ED_copy[], int subsystem)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the 'subsystem' witness data information         *
\***************************************************************/
{
  int i;

  if (max_threads > 1)
  {
    for (i = 0; i < max_threads; i++)
    { // set the pointers to NULL
      ED_copy[i].EqD->witnessData_mp[subsystem].startPts = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].endPts_in = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].endPts = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].finalTs = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].condition_nums = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].endPt_retVals = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].endPt_types = NULL;
      ED_copy[i].EqD->witnessData_mp[subsystem].higherDim = NULL;

      // clear B and p
      clear_mat_mp(ED_copy[i].EqD->witnessData_mp[subsystem].B);
      clear_vec_mp(ED_copy[i].EqD->witnessData_mp[subsystem].p);
    }
  }

  return;
}

void clear_omp_eqData_d(basic_eval_data_d *BED, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear BED where BED was used for OpenMP tracking       *
*  Since each copy has its own version of gamma_mp & coeff_mp,  *
*  these need to be cleared - otherwise set pointers to NULL    *
\***************************************************************/
{
  int i, j;

  // clear MP data first, if needed
  if (MPType == 2)
  { // clear gamma_rat, gamma_mp & set coeff_rat to NULL
    clear_mp(BED->EqD->gamma_mp);
    clear_rat(BED->EqD->gamma_rat);
    BED->EqD->coeff_rat = NULL;

    // clear coeff_mp
    for (i = BED->EqD->num_funcs - 1; i >= 0; i--)
    {
      for (j = BED->EqD->degrees[i][0] - 1; j >= 0; j--)
      {
        clear_mp(BED->EqD->coeff_mp[i][j]);
      }
      free(BED->EqD->coeff_mp[i]);
    }
    free(BED->EqD->coeff_mp);

    // clear the MP witness data and stage data
    free(BED->EqD->witnessData_mp);
    free(BED->EqD->stageData_mp);
  }

  // set degrees and coeff_d to NULL
  BED->EqD->degrees = NULL;
  BED->EqD->coeff_d = NULL;

  // setup InstCount to NULL
  BED->EqD->startSub = BED->EqD->endSub = BED->EqD->startFunc = BED->EqD->endFunc = BED->EqD->startJvsub = BED->EqD->endJvsub = BED->EqD->startJv = BED->EqD->endJv = NULL;
  BED->EqD->subFuncsBelow = NULL;

  // clear the _d witness data and stage data
  free(BED->EqD->witnessData_d);
  free(BED->EqD->stageData_d);

  // free EqD
  free(BED->EqD);
  
  return;
}

void clear_omp_eqData_mp(basic_eval_data_mp *BED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear BED where BED was used for OpenMP tracking       *
*  Since each copy has its own version of gamma_mp, clear it    *
\***************************************************************/
{
  // clear gamma_mp
  clear_mp(BED->EqD->gamma_mp);

  BED->EqD->coeff_mp = NULL;

  // clear the witness data and stage data
  free(BED->EqD->witnessData_mp);
  free(BED->EqD->stageData_mp);

  // set degrees to NULL
  BED->EqD->degrees = NULL;

  // setup InstCount to NULL
  BED->EqD->startSub = BED->EqD->endSub = BED->EqD->startFunc = BED->EqD->endFunc = BED->EqD->startJvsub = BED->EqD->endJvsub = BED->EqD->startJv = BED->EqD->endJv = NULL;
  BED->EqD->subFuncsBelow = NULL;

  // free EqD
  free(BED->EqD);

  return;
}

