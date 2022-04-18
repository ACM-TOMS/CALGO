// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

void setupTD_startPoints_d(char pointsIN[], char pointsOUT[], int size, int *degs, patch_eval_data_d *PED)
{ // reads points from pointsIN and moves them to points on the patch and writes them to pointsOUT
  // ASSUME ONLY 1 PATCH!!!!!!!

  int i, j, num_points;
  comp_d multiplier;
  point_d tempPoint;
  vec_d patchVals;
  FILE *IN, *OUT;

  // initialize
  init_point_d(tempPoint, size + 1);
  tempPoint->size = size + 1; // extra coord at the beginning since patch depends on r + 1 variables
  init_vec_d(patchVals, 0);

  // create the start points
  TDstartMaker_d(degs, size);

  // open the files
  IN = fopen(pointsIN, "r");
  OUT = fopen(pointsOUT, "w");
  if (IN == NULL)
  {
    printf("ERROR: '%s' does not exist!!!\n\n", pointsIN);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  
  // read in the number of points
  fscanf(IN, "%d\n\n", &num_points);
  // write the number of points
  fprintf(OUT, "%d\n\n", num_points);
  
  for (i = 0; i < num_points; i++)
  { // read in the ith point
    set_double_d(&tempPoint->coord[0], 1.0, 0.0);
    for (j = 1; j < tempPoint->size; j++)
      fscanf(IN, "%lf %lf;\n", &tempPoint->coord[j].r, &tempPoint->coord[j].i);
    fscanf(IN, "\n");

    // calculate the value of the patch at this point
    mul_mat_vec_d(patchVals, PED->patchCoeff, tempPoint);
    
    // adjust tempPoint based on patchVals[0] and write to OUT
    recip_d(multiplier, &patchVals->coord[0]);
    for (j = 0; j < tempPoint->size; j++)
    {
      mul_d(&tempPoint->coord[j], &tempPoint->coord[j], multiplier);
      fprintf(OUT, "%.15e %.15e;\n", tempPoint->coord[j].r, tempPoint->coord[j].i);
    }
    fprintf(OUT, "\n");
  }

  // close the files
  fclose(IN);
  fclose(OUT);

  // clear
  clear_point_d(tempPoint);
  clear_vec_d(patchVals);

  return;
}
    
void setupPreProcData(char preprocFile[], preproc_data *PPD)
{ // setup PPD with the information stored in preprocFile

  int i, total_gp;
  FILE *ppIN = fopen(preprocFile, "r");
  if (ppIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", preprocFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in the first line of data
  fscanf(ppIN, "%d %d %d\n", &PPD->num_funcs, &PPD->num_hom_var_gp, &PPD->num_var_gp);

  total_gp = PPD->num_hom_var_gp + PPD->num_var_gp;

  PPD->size = (int *)bmalloc(total_gp * sizeof(int));
  PPD->type = (int *)bmalloc(total_gp * sizeof(int));

  // read in the type & sizes of the variable groups
  for (i = 0; i < total_gp; i++)
    fscanf(ppIN, "%d %d\n", &PPD->type[i], &PPD->size[i]);

  fclose(ppIN);

  return;
}

void setupSquareSystem_d(prog_t *Prog, int finalSize, preproc_data *PPD, char degreeFile[], square_system_eval_data_d *SSED, int adjustDegrees)
{ // setup the square system - finalSize is the final size of the 'square system' - without the patch

  int i, j, total_vars, max, max_loc = 0;
  int *tempInts = NULL;

  FILE *degIN = fopen(degreeFile, "r");
  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // calcualate total number of variables
  total_vars = 0;
  for (i = 0; i < (PPD->num_var_gp + PPD->num_hom_var_gp); i++)
    total_vars += PPD->size[i] + PPD->type[i];

  SSED->Prog = Prog;
  SSED->size_f = PPD->num_funcs;
  SSED->size_r = finalSize;
  SSED->noChanges = 0;

  SSED->P = (int *)bmalloc(SSED->size_f * sizeof(int));
  SSED->W = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
    SSED->W[i] = (int *)bmalloc((SSED->size_f - SSED->size_r) * sizeof(int));

  if (PPD->num_hom_var_gp + PPD->num_var_gp == 1)
  { // 1-hom case!

    // setup the original degrees
    SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
    for (i = 0; i < SSED->size_f; i++)
      fscanf(degIN, "%d\n\n", &SSED->orig_degrees[i]); 

    // find B & B_perp
    if (total_vars == SSED->size_r + 1)
    { // initialize B & B_perp
      init_mat_d(SSED->B, total_vars, total_vars);
      SSED->B->rows = SSED->B->cols = total_vars; // number of x variables
      init_mat_d(SSED->B_perp, total_vars, 0);
      SSED->B_perp->rows = total_vars;
      SSED->B_perp->cols = 0;

      make_matrix_ID_d(SSED->B, total_vars, total_vars);  // if we are not changing the number of variables, then we can stick with the same ones!!
    }
    else
    { // initialize B to a square matrix
      init_mat_d(SSED->B, total_vars, total_vars);
      make_matrix_random_d(SSED->B, total_vars, total_vars); // number of x variables - have it be square to make the matrix at first

      // now find B & B_perp
      init_mat_d(SSED->B_perp, total_vars, total_vars - (SSED->size_r + 1));
      SSED->B->cols = SSED->size_r + 1; // number of y variables
      SSED->B_perp->rows = total_vars;
      SSED->B_perp->cols = total_vars - (SSED->size_r + 1);

      for (i = 0; i < total_vars; i++)
        for (j = 0; j < SSED->B_perp->cols; j++)
        {
          conjugate_d(&SSED->B_perp->entry[i][j], &SSED->B->entry[i][j + SSED->B->cols]); // conjugate needed since G.S. returns a matrix
                                                                                          // whose rows are 'conjugate dot product' perpendicular
        }
    }

    if (!adjustDegrees && SSED->size_f == SSED->size_r)
    { // the calling function does not want the degrees adjusted and so since it is square, there is no need to adjust list based on degree
      SSED->noChanges = 1;
      for (i = 0; i < SSED->size_f; i++)
        SSED->P[i] = i;
    }
    else
    { // to create P, we need to find the location of the ith largest degree and store as P[i]
      tempInts = (int *)bmalloc(SSED->size_f * sizeof(int));
      for (i = 0; i < SSED->size_f; i++)
        tempInts[i] = 1; // will be = 0 when the function has been used

      for (i = 0; i < SSED->size_f; i++)
      {
        max = -1;
        for (j = 0; j < SSED->size_f; j++)
        {
          if ((max < SSED->orig_degrees[j]) && (tempInts[j] == 1))
          {
            max = SSED->orig_degrees[j];
            max_loc = j;
          }
        }
        SSED->P[i] = max_loc;
        tempInts[max_loc] = 0;
      }
    }

    // setup new_degrees
    SSED->new_degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
    for (i = 0; i < SSED->size_r; i++)
      SSED->new_degrees[i] = SSED->orig_degrees[SSED->P[i]];
    
    // setup A
    init_mat_d(SSED->A, SSED->size_r, SSED->size_f - SSED->size_r);
    make_matrix_random_d(SSED->A, SSED->size_r, SSED->size_f - SSED->size_r);

    // to setup W, we need to find the difference in the exponents
    for (i = 0; i < SSED->size_r; i++)
      for (j = 0; j < SSED->size_f - SSED->size_r; j++)
        SSED->W[i][j] = SSED->orig_degrees[SSED->P[i]] - SSED->orig_degrees[SSED->P[j+SSED->size_r]];

    // to find max_of_W, we simply move through W and find the maximum
    SSED->max_of_W = 0;
    for (i = 0; i < SSED->size_r; i++)
      for (j = 0; j < SSED->size_f - SSED->size_r; j++)
        if (SSED->max_of_W < SSED->W[i][j])
          SSED->max_of_W = SSED->W[i][j];
  }
  else if (PPD->num_hom_var_gp + PPD->num_var_gp > 1)
  { // mhom case
    if (finalSize != SSED->size_f)
    {
      printf("ERROR: To use a multihomogenous homotopy, the system must be square (finalSize (%d) != size_f (%d)).\n", finalSize, SSED->size_f);
      bexit(ERROR_INVALID_SIZE);
    }

    // with mhom, there are no changes are done
    SSED->noChanges = 1;

    // change of variable matrix is identity
    init_mat_d(SSED->B, total_vars, total_vars);
    make_matrix_ID_d(SSED->B, total_vars, total_vars);

    // permutation is the identity
    for (i = 0; i < SSED->size_f; i++)
      SSED->P[i] = i;

    // W & A & B_perp have no columns!!
    init_mat_d(SSED->A, SSED->size_f, 0);
    init_mat_d(SSED->B_perp, total_vars, 0);
    SSED->B_perp->rows = total_vars;
    SSED->A->rows = SSED->size_f;
    SSED->A->cols = SSED->B_perp->cols = 0;

    // setup orig_degrees to be the total of the mhom degrees
    SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
    SSED->new_degrees  = (int *)bmalloc(SSED->size_f * sizeof(int));
    for (i = 0; i < SSED->size_f; i++)
    {
      SSED->orig_degrees[i] = 0;
      for (j = 0; j < PPD->num_hom_var_gp + PPD->num_var_gp; j++)
      {
        fscanf(degIN, "%d\n", &max);
        SSED->orig_degrees[i] += max;
      }
      fscanf(degIN, "\n");
      SSED->new_degrees[i] = SSED->orig_degrees[i];
    }

    SSED->max_of_W = 0;
  }
  else
  {
    printf("ERROR: In setupSquareSystem, there are %d variables groups!\n", PPD->num_hom_var_gp + PPD->num_var_gp);
    bexit(ERROR_INVALID_SIZE);
  }

  free(tempInts);

  return;
}

void setupCoeffInSS_d(char degreeFile[], int *P, preproc_data *PPD, start_system_eval_data_d *SSED)
{ // setup the mhom start system - i.e. setup the coeffs in the start system

// SINCE THIS IS ONLY USED IN MHOM, the total number of variables described in PPD is the total number of variables used throughout

  int i, j, k, total_vars, total_deg, curr_loc, var_gp_count;
  int **mhomDeg, *linearGroup, **varGroupArray;
  FILE *degIN = fopen(degreeFile, "r");

  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // calculate the total number of variables that are needed
  total_vars = 0;
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
    total_vars += PPD->size[i] + PPD->type[i];

  SSED->coeff_cols = total_vars;

  mhomDeg = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
    mhomDeg[i] = (int *)bmalloc((PPD->num_hom_var_gp + PPD->num_var_gp) * sizeof(int));

  // read in the mhom degrees - and find the total number of linears that will be needed
  total_deg = 0;
  for (i = 0; i < SSED->size_r; i++)
  {
    for (j = 0; j < (PPD->num_hom_var_gp + PPD->num_var_gp); j++)
    {
      fscanf(degIN, "%d\n", &mhomDeg[i][j]);
      total_deg += mhomDeg[i][j];
    }
    fscanf(degIN, "\n"); // extra 
  }

  // convert to a list of the variable groups that are needed for the linears - use the Permuation P to get the list correct!
  linearGroup = (int *)bmalloc(total_deg * sizeof(int));
  curr_loc = 0;
  for (i = 0; i < SSED->size_r; i++)
    for (j = 0; j < (PPD->num_hom_var_gp + PPD->num_var_gp); j++)
      for (k = 0; k < mhomDeg[P[i]][j]; k++)
      {
        linearGroup[curr_loc] = j;
        curr_loc++;
      }

  // setup an array to tell which variables are in which group
  varGroupArray = (int **)bmalloc((PPD->num_hom_var_gp + PPD->num_var_gp) * sizeof(int *));
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
    varGroupArray[i] = (int *)bmalloc(total_vars * sizeof(int));

  var_gp_count = 0;
  curr_loc = PPD->num_var_gp;
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
  {
    for (j = 0; j < total_vars; j++)
      // initialize to zero
      varGroupArray[i][j] = 0;

    // put 1 in the appropriate places
    if (PPD->type[i] == 1)
    { // needed for the extra homogenizing variable that was added 
      varGroupArray[i][var_gp_count] = 1;
      var_gp_count++;
    }
    for (j = 0; j < PPD->size[i]; j++)
    {
      varGroupArray[i][curr_loc] = 1;
      curr_loc++;
    }
  }

  SSED->coeff = (comp_d **)bmalloc(total_deg * sizeof(comp_d *));
  for (i = 0; i < total_deg; i++)
    SSED->coeff[i] = (comp_d *)bmalloc(total_vars * sizeof(comp_d));

  // setup the coeff
  for (i = 0; i < total_deg; i++)
    // at linear i, the variable group that this corresponds to is linearGroup[i]
    for (j = 0; j < total_vars; j++)
      if (varGroupArray[linearGroup[i]][j] == 1)
      {
        get_comp_rand_d(SSED->coeff[i][j]);
      }
      else
      {
        set_zero_d(SSED->coeff[i][j]);
      }

  // close files & release memory
  fclose(degIN);

  for (i = (PPD->num_hom_var_gp + PPD->num_var_gp - 1); i >= 0; i--)
    free(varGroupArray[i]);
  free(varGroupArray);

  for (i = SSED->size_r - 1; i >= 0; i--)
    free(mhomDeg[i]);
  free(mhomDeg);

  free(linearGroup);

  return;
}

void setupStartSystem_d(int SSType, int size, int *deg, int *P, preproc_data *PPD, char degreeFile[], start_system_eval_data_d *SSED)
{ // if SSType == 0, we are setuping up for total degree, else, for mhom

  int i;

  // setup the startSystem
  SSED->startSystemType = SSType;
  SSED->size_r = size;
  SSED->degrees = (int *)bmalloc(size * sizeof(int));

  get_comp_rand_d(SSED->gamma);

  SSED->max_degree = -1;
  for (i = 0; i < size; i++)
  {
    SSED->degrees[i] = deg[i];
    if (deg[i] > SSED->max_degree)
      SSED->max_degree = deg[i];
  }

  if (SSType == 1) // need to setup the coeff
    setupCoeffInSS_d(degreeFile, P, PPD, SSED);

  return;
}

void setupPatch_d(int PatchType, patch_eval_data_d *PED, void const *ptr1, void const *ptr2)
{ // if PatchType == 1, we are setting up the patches by Q*B, 
  // else if PatchType == 2, we are setting up for 1-hom
  // else, we generate randomly based on variable groups

  int i, j, curr_loc, total_vars, var_gp_count;
  comp_d tempComp;

  if (PatchType == 1)
  {
    vec_d *vecQ = (vec_d *)ptr1;
    mat_d *matB = (mat_d *)ptr2;

    init_mat_d(PED->patchCoeff, 1, (*matB)->cols);
    PED->num_patches = PED->patchCoeff->rows = 1;
    PED->patchCoeff->cols = (*matB)->cols;

    // = Q*B
    for (j = 0; j < PED->patchCoeff->cols; j++)
    {
      set_zero_d(&PED->patchCoeff->entry[0][j]);
      for (i = 0; i < (*vecQ)->size; i++)
      {
        mul_d(tempComp, &(*vecQ)->coord[i], &(*matB)->entry[i][j]);
        add_d(&PED->patchCoeff->entry[0][j], &PED->patchCoeff->entry[0][j], tempComp);
      }
    }
  }
  else if (PatchType == 2)
  {
    int *numVars = (int *)ptr1;

    // one patch
    init_mat_d(PED->patchCoeff, 1, *numVars);
    PED->num_patches = PED->patchCoeff->rows = 1;
    PED->patchCoeff->cols = *numVars;

    for (i = 0; i < PED->patchCoeff->cols; i++)
      get_comp_rand_d(&PED->patchCoeff->entry[0][i]);
  }
  else
  {
    preproc_data *PPD = (preproc_data *)ptr1;

    // calculate the total number of variables
    total_vars = 0;
    for (i = 0; i < PPD->num_hom_var_gp + PPD->num_var_gp; i++)
      total_vars += PPD->type[i] + PPD->size[i];

    // one patch for each variable group
    PED->num_patches = PPD->num_hom_var_gp + PPD->num_var_gp;
    init_mat_d(PED->patchCoeff, PED->num_patches, total_vars);
    PED->patchCoeff->rows = PED->num_patches;
    PED->patchCoeff->cols = total_vars; 

    // generate a patch for each of these variable groups
    var_gp_count = 0; // the number of variable_groups that we have been through
    curr_loc = PPD->num_var_gp; // current location of the variables
    for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
    { // zero out the row
      for (j = 0; j < total_vars; j++)
      { 
        set_zero_d(&PED->patchCoeff->entry[i][j]);
      }
      // put random numbers in the appropriate places
      if (PPD->type[i] == 1)
      { // needed for the extra homogenizing variable that was added
        get_comp_rand_d(&PED->patchCoeff->entry[i][var_gp_count]);
        var_gp_count++;
      }
      for (j = 0; j < PPD->size[i]; j++)
      {
        get_comp_rand_d(&PED->patchCoeff->entry[i][curr_loc]);
        curr_loc++;
      }
    }
  }

  return;
}

void setupBasicEval_d(char preprocFile[], char degreeFile[], prog_t *dummyProg, int squareSize, int patchType, int ssType, int MPType, void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4, basic_eval_data_d *BED, int adjustDegrees)
{
  int i;

  setupPreProcData(preprocFile, &BED->preProcData);
  setupSquareSystem_d(dummyProg, squareSize, &BED->preProcData, degreeFile, &BED->squareSystem, adjustDegrees); // NOTE: squareSystem must be setup before the patch!!!
  setupPatch_d(patchType, &BED->patch, ptr1, ptr2);

  // find the degrees of the 'square system'
  int *deg = (int *)bmalloc(BED->squareSystem.size_r * sizeof(int));
  for (i = 0; i < BED->squareSystem.size_r; i++)
    deg[i] = BED->squareSystem.orig_degrees[BED->squareSystem.P[i]];

  setupStartSystem_d(ssType, BED->squareSystem.size_r, deg, BED->squareSystem.P, &BED->preProcData, degreeFile, &BED->startSystem);

  if (MPType == 2)
  { // using AMP - initialize using 16 digits & 64-bit precison 
    int digits = 16, prec = 64;

    // setup preProcData
    setupPreProcData(preprocFile, &BED->BED_mp->preProcData);
    // setup the square system
    setupSquareSystem_d_to_mp(&BED->squareSystem, &BED->BED_mp->squareSystem, digits, prec);
    // setup the patch
    setupPatch_d_to_mp(&BED->patch, &BED->BED_mp->patch, digits, prec, patchType, dummyProg->numVars, ptr3, ptr4); 
    // setup the start system
    setupStartSystem_d_to_mp(&BED->startSystem, &BED->BED_mp->startSystem, digits, prec);
  }

  free(deg);

  return;
}

void setupBEDUsingUserHom_d(prog_t *dummyProg, int MPType, basic_eval_data_d *BED)
{ // set everything to nothing expect for the s.l.p.!!

  // patch things
  init_mat_d(BED->patch.patchCoeff, 0, 0);
  BED->patch.num_patches = BED->patch.patchCoeff->rows = BED->patch.patchCoeff->cols = 0;

  // start system things
  BED->startSystem.startSystemType = -1;
  BED->startSystem.degrees = NULL;
  BED->startSystem.size_r = BED->startSystem.max_degree = BED->startSystem.coeff_cols = 0;
  get_comp_rand_d(BED->startSystem.gamma);
  
  // square system things
  BED->squareSystem.Prog = dummyProg;
  BED->squareSystem.orig_degrees = BED->squareSystem.new_degrees = BED->squareSystem.P = NULL;
  BED->squareSystem.W = NULL;
  init_mat_d(BED->squareSystem.B, 0, 0);
  BED->squareSystem.size_f = BED->squareSystem.B->rows = BED->squareSystem.B->cols = 0;
  init_mat_d(BED->squareSystem.B_perp, 0, 0);
  BED->squareSystem.B_perp->rows = BED->squareSystem.B_perp->cols = 0;
  init_mat_d(BED->squareSystem.A, 0, 0);
  BED->squareSystem.max_of_W = BED->squareSystem.A->rows = BED->squareSystem.A->cols = 0;
  BED->squareSystem.size_r = 0;

  // preproc data things
  BED->preProcData.num_funcs = BED->preProcData.num_hom_var_gp = BED->preProcData.num_var_gp = 0;
  BED->preProcData.type = BED->preProcData.size = NULL;

  if (MPType == 2)
  { // using AMP - initialize using 16 digits & 64-bit precision
    int digits = 16, prec = 64;

    // setup preProcData
    BED->BED_mp->preProcData.num_funcs = BED->BED_mp->preProcData.num_hom_var_gp = BED->BED_mp->preProcData.num_var_gp = 0;
    BED->BED_mp->preProcData.type = BED->BED_mp->preProcData.size = NULL;
    // setup the square system
    setupSquareSystem_d_to_mp(&BED->squareSystem, &BED->BED_mp->squareSystem, digits, prec);
    // setup the patch
    setupPatch_d_to_mp(&BED->patch, &BED->BED_mp->patch, digits, prec, 0, 0, NULL, NULL);
   // setup the start system
    setupStartSystem_d_to_mp(&BED->startSystem, &BED->BED_mp->startSystem, digits, prec);
  }

  return;
}

void setupBEDUsingParamHom_d(prog_t *dummyProg, char *preprocFile, int MPType, basic_eval_data_d *BED)
{ // set SLP & patch only

  // setup preProcData
  setupPreProcData(preprocFile, &BED->preProcData);

  // setup patch
  setupPatch_d(0, &BED->patch, &BED->preProcData, NULL);

  // start system things
  BED->startSystem.startSystemType = -1;
  BED->startSystem.degrees = NULL;
  BED->startSystem.size_r = BED->startSystem.max_degree = BED->startSystem.coeff_cols = 0;
  get_comp_rand_d(BED->startSystem.gamma);

  // square system things
  BED->squareSystem.Prog = dummyProg;
  BED->squareSystem.orig_degrees = BED->squareSystem.new_degrees = BED->squareSystem.P = NULL;
  BED->squareSystem.W = NULL;
  init_mat_d(BED->squareSystem.B, 0, 0);
  BED->squareSystem.size_f = BED->squareSystem.B->rows = BED->squareSystem.B->cols = 0;
  init_mat_d(BED->squareSystem.B_perp, 0, 0);
  BED->squareSystem.B_perp->rows = BED->squareSystem.B_perp->cols = 0;
  init_mat_d(BED->squareSystem.A, 0, 0);
  BED->squareSystem.max_of_W = BED->squareSystem.A->rows = BED->squareSystem.A->cols = 0;
  BED->squareSystem.size_r = 0;

  if (MPType == 2)
  { // using AMP - initialize using 16 digits & 64-bit precision
    int digits = 16, prec = 64;

    // setup preProcData
    setupPreProcData(preprocFile, &BED->BED_mp->preProcData);
    // setup the square system
    setupSquareSystem_d_to_mp(&BED->squareSystem, &BED->BED_mp->squareSystem, digits, prec);
    // setup the patch
    setupPatch_d_to_mp(&BED->patch, &BED->BED_mp->patch, digits, prec, 0, 0, NULL, NULL);
   // setup the start system
    setupStartSystem_d_to_mp(&BED->startSystem, &BED->BED_mp->startSystem, digits, prec);
  }

  return;
}

////// AMP (CON)VERSIONS /////////

void setupPatch_d_to_mp(patch_eval_data_d *PED_d, patch_eval_data_mp *PED_mp, int digits, int prec, int patchType, int numVars, void const *ptr1, void const *ptr2)
{ // convert _d patch to _mp patch
  // if PatchType == 1, we are setting up the patches by Q*B,
  // else we create new extended precision random numbers when the _d entry is not 0

  int i, j, max_prec = 1024, rows = PED_d->patchCoeff->rows, cols = PED_d->patchCoeff->cols;

  // initialize & set to current precision
  init_mat_mp2(PED_mp->patchCoeff, rows, cols, prec);
  PED_mp->curr_prec = prec;
  PED_mp->num_patches = PED_d->num_patches;
  PED_mp->patchCoeff->rows = rows;
  PED_mp->patchCoeff->cols = cols;

  // allocate for the coeff
  init_mat_rat(PED_mp->patchCoeff_rat, rows, cols);

  if (patchType == 1)
  {
    mpq_t ***Q_rat = (mpq_t ***)ptr1;
    mpq_t ****B_rat = (mpq_t ****)ptr2;

    vec_mat_mul_rat(PED_mp->patchCoeff_rat[0], *Q_rat, *B_rat, numVars, numVars, cols, 1);

    // setup the patch coefficients
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpf_set_q(PED_mp->patchCoeff->entry[i][j].r, PED_mp->patchCoeff_rat[i][j][0]);
        mpf_set_q(PED_mp->patchCoeff->entry[i][j].i, PED_mp->patchCoeff_rat[i][j][1]);
        PED_d->patchCoeff->entry[i][j].r = mpq_get_d(PED_mp->patchCoeff_rat[i][j][0]);
        PED_d->patchCoeff->entry[i][j].i = mpq_get_d(PED_mp->patchCoeff_rat[i][j][1]);
      }
  }
  else
  { // setup the patch coefficients
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
        if (PED_d->patchCoeff->entry[i][j].r == 0 && PED_d->patchCoeff->entry[i][j].i == 0)
        { // just copy over 0 to _mp & _rat
          comp_d_to_mp_rat(&PED_mp->patchCoeff->entry[i][j], PED_mp->patchCoeff_rat[i][j], &PED_d->patchCoeff->entry[i][j], digits, prec, 0, 1);
        }
        else
        { // setup to 'extended precision'
          get_comp_rand_rat(&PED_d->patchCoeff->entry[i][j], &PED_mp->patchCoeff->entry[i][j], PED_mp->patchCoeff_rat[i][j], prec, max_prec, 0, 1);
        }
  }

  return;
}

void setupSquareSystem_d_to_mp(square_system_eval_data_d *SSED_d, square_system_eval_data_mp *SSED_mp, int digits, int prec)
{ // convert _d square system to _mp square system
  int i, j;

  // copy over the 'regular' information
  SSED_mp->Prog = SSED_d->Prog;
  SSED_mp->size_f = SSED_d->size_f;
  SSED_mp->size_r = SSED_d->size_r;
  SSED_mp->max_of_W = SSED_d->max_of_W;
  SSED_mp->noChanges = SSED_d->noChanges;

  SSED_mp->orig_degrees = (int *)bmalloc(SSED_mp->size_f * sizeof(int));
  SSED_mp->P = (int *)bmalloc(SSED_mp->size_f * sizeof(int));
  for (i = 0; i < SSED_mp->size_f; i++)
  {
    SSED_mp->orig_degrees[i] = SSED_d->orig_degrees[i];
    SSED_mp->P[i] = SSED_d->P[i];
  }

  SSED_mp->new_degrees = (int *)bmalloc(SSED_mp->size_r * sizeof(int));
  SSED_mp->W = (int **)bmalloc(SSED_mp->size_r * sizeof(int *));
  for (i = 0; i < SSED_mp->size_r; i++)
  {
    SSED_mp->new_degrees[i] = SSED_d->new_degrees[i];
    SSED_mp->W[i] = (int *)bmalloc((SSED_mp->size_f - SSED_mp->size_r) * sizeof(int));
    for (j = 0; j < (SSED_mp->size_f - SSED_mp->size_r); j++)
      SSED_mp->W[i][j] = SSED_d->W[i][j];
  }
  
  SSED_mp->curr_prec = prec;

  // setup the matrix B
  init_mat_mp2(SSED_mp->B, SSED_d->B->rows, SSED_d->B->cols, prec);
  init_mat_rat(SSED_mp->B_rat, SSED_d->B->rows, SSED_d->B->cols);
  mat_d_to_mp_rat(SSED_mp->B, SSED_mp->B_rat, SSED_d->B, digits, prec, 0, 0);

  // setup the matrix B_perp
  init_mat_mp2(SSED_mp->B_perp, SSED_d->B_perp->rows, SSED_d->B_perp->cols, prec);
  init_mat_rat(SSED_mp->B_perp_rat, SSED_d->B_perp->rows, SSED_d->B_perp->cols);
  mat_d_to_mp_rat(SSED_mp->B_perp, SSED_mp->B_perp_rat, SSED_d->B_perp, digits, prec, 0, 0);

  // setup the matrix A
  init_mat_mp2(SSED_mp->A, SSED_d->A->rows, SSED_d->A->cols, prec);
  init_mat_rat(SSED_mp->A_rat, SSED_d->A->rows, SSED_d->A->cols);
  mat_d_to_mp_rat(SSED_mp->A, SSED_mp->A_rat, SSED_d->A, digits, prec, 0, 0);

  return;
}

void setupStartSystem_d_to_mp(start_system_eval_data_d *SSED_d, start_system_eval_data_mp *SSED_mp, int digits, int prec)
{ // convert _d patch to _mp start system
  int i, j, total_deg = 0, max_prec = 1024;

  SSED_mp->startSystemType = SSED_d->startSystemType;
  SSED_mp->size_r = SSED_d->size_r;
  SSED_mp->max_degree = SSED_d->max_degree;
  SSED_mp->coeff_cols = SSED_d->coeff_cols;

  SSED_mp->degrees = (int *)bmalloc(SSED_mp->size_r * sizeof(int));
  for (i = 0; i < SSED_mp->size_r; i++)
    total_deg += SSED_mp->degrees[i] = SSED_d->degrees[i];

  SSED_mp->curr_prec = prec;

  SSED_mp->gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); // for real & imaginary
  get_comp_rand_rat(SSED_d->gamma, SSED_mp->gamma, SSED_mp->gamma_rat, prec, max_prec, 1, 1);

  if (SSED_mp->startSystemType == 1)
  { // if we are using an mhom start system, we need to copy over the coeff
    SSED_mp->coeff = (comp_mp **)bmalloc(total_deg * sizeof(comp_mp *));
    SSED_mp->coeff_rat = (mpq_t ***)bmalloc(total_deg * sizeof(mpq_t **));
    for (i = 0; i < total_deg; i++)
    {
      SSED_mp->coeff[i] = (comp_mp *)bmalloc(SSED_mp->coeff_cols * sizeof(comp_mp));
      SSED_mp->coeff_rat[i] = (mpq_t **)bmalloc(SSED_mp->coeff_cols * sizeof(mpq_t *));
      for (j = 0; j < SSED_mp->coeff_cols; j++)
      {
        SSED_mp->coeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));  // for real & imaginary
        if (SSED_d->coeff[i][j]->r == 0 && SSED_d->coeff[i][j]->i == 0)
        { // just copy over 0 to _mp & _rat
          comp_d_to_mp_rat(SSED_mp->coeff[i][j], SSED_mp->coeff_rat[i][j], SSED_d->coeff[i][j], digits, prec, 1, 1);
        }
        else
        { // setup to 'extended precision'
          get_comp_rand_rat(SSED_d->coeff[i][j], SSED_mp->coeff[i][j], SSED_mp->coeff_rat[i][j], prec, max_prec, 1, 1);
        }
      }
    }
  }

  return;
}

////// CHANGE PRECISION ////////

void changePatchPrec_mp(int new_prec, patch_eval_data_mp *PED)
{ // convert _d patch to _mp patch

  int i, j;

  if (new_prec != PED->curr_prec)
  { // we need to change the precision
    int rows = PED->patchCoeff->rows, cols = PED->patchCoeff->cols;

    PED->curr_prec = new_prec;

    // change precision for patch coeff
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        setprec_mp(&PED->patchCoeff->entry[i][j], PED->curr_prec);
        mpf_set_q(PED->patchCoeff->entry[i][j].r, PED->patchCoeff_rat[i][j][0]);
        mpf_set_q(PED->patchCoeff->entry[i][j].i, PED->patchCoeff_rat[i][j][1]);
      }
  }

  return;
}

void changeStartPrec_mp(int new_prec, start_system_eval_data_mp *SSED)
{
  int i, j, total_deg = 0;

  if (new_prec != SSED->curr_prec)
  { // we need to change the precision
 
    SSED->curr_prec = new_prec;

    // change precision for gamma
    setprec_mp(SSED->gamma, SSED->curr_prec);
    mpf_set_q(SSED->gamma->r, SSED->gamma_rat[0]);
    mpf_set_q(SSED->gamma->i, SSED->gamma_rat[1]);

    if (SSED->startSystemType == 1)
    { // change precision for coeff
      for (i = 0; i < SSED->size_r; i++)
        total_deg += SSED->degrees[i];
      for (i = 0; i < total_deg; i++)
        for (j = 0; j < SSED->coeff_cols; j++)
        {
          setprec_mp(SSED->coeff[i][j], SSED->curr_prec);
          mpf_set_q(SSED->coeff[i][j]->r, SSED->coeff_rat[i][j][0]);
          mpf_set_q(SSED->coeff[i][j]->i, SSED->coeff_rat[i][j][1]);
        }
    }
  }

  
  return;
}

void changeSquarePrec_mp(int new_prec, square_system_eval_data_mp *SSED)
{
  int i, j;

  // set the SLP to the correct precision
  SSED->Prog->precision = new_prec;

  if (new_prec != SSED->curr_prec)
  { // we need to change the precision

    SSED->curr_prec = new_prec;

    // change the precision for A
    for (i = 0; i < SSED->A->rows; i++)
      for (j = 0; j < SSED->A->cols; j++)
      {
        setprec_mp(&SSED->A->entry[i][j], SSED->curr_prec);
        mpf_set_q(SSED->A->entry[i][j].r, SSED->A_rat[i][j][0]);
        mpf_set_q(SSED->A->entry[i][j].i, SSED->A_rat[i][j][1]);
      }
  
    // change the precision for B
    for (i = 0; i < SSED->B->rows; i++)
      for (j = 0; j < SSED->B->cols; j++)
      {
        setprec_mp(&SSED->B->entry[i][j], SSED->curr_prec);
        mpf_set_q(SSED->B->entry[i][j].r, SSED->B_rat[i][j][0]);
        mpf_set_q(SSED->B->entry[i][j].i, SSED->B_rat[i][j][1]);
      }

    // change the precision for B_perp
    for (i = 0; i < SSED->B_perp->rows; i++)
      for (j = 0; j < SSED->B_perp->cols; j++)
      {
        setprec_mp(&SSED->B_perp->entry[i][j], SSED->curr_prec);
        mpf_set_q(SSED->B_perp->entry[i][j].r, SSED->B_perp_rat[i][j][0]);
        mpf_set_q(SSED->B_perp->entry[i][j].i, SSED->B_perp_rat[i][j][1]);
      }
  }

  return;
}

int change_square_prec(void const *Sq, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  changeSquarePrec_mp(prec, (square_system_eval_data_mp *)Sq);
  return 0;
}

int change_basic_eval_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for standard zero dimensional solving *
\***************************************************************/
{
  changeBasicEvalPrec_mp(prec, (basic_eval_data_mp *)ED);
  return 0;
}

void changeBasicEvalPrec_mp(int new_prec, basic_eval_data_mp *BED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for the main part of BED              *
\***************************************************************/
{
  // change the precision for the patch
  changePatchPrec_mp(new_prec, &BED->patch);
  // change the precision for the start system
  changeStartPrec_mp(new_prec, &BED->startSystem);
  // change the precision for the square system
  changeSquarePrec_mp(new_prec, &BED->squareSystem);

  return;
}

////// MP VERSIONS ///////////

void setupPatch_mp(int PatchType, patch_eval_data_mp *PED, int digits, int prec, void const *ptr1, void const *ptr2)
{ // if PatchType == 1, we are setuping up the patches by Q*B - only for cascade which is not yet implemented in _mp!
  // else if PatchType == 2, we are setuping up for 1-hom
  // else, we generate randomly based on variable groups

  int i, j, total_vars, curr_loc, var_gp_count;

  PED->curr_prec = prec;

  if (PatchType == 1)
  { // patch = Q*B
    printf("ERROR: trying to make patch = Q*B - a cascade thing - using multiprecion!!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (PatchType == 2)
  {
    int *numVars = (int *)ptr1;

    // one patch
    init_mat_mp2(PED->patchCoeff, 1, *numVars, prec);
    PED->num_patches = PED->patchCoeff->rows = 1;
    PED->patchCoeff->cols = *numVars;

    for (i = 0; i < PED->patchCoeff->cols; i++)
      get_comp_rand_mp(&PED->patchCoeff->entry[0][i]);

    init_mat_rat(PED->patchCoeff_rat, PED->patchCoeff->rows, PED->patchCoeff->cols);
    for (i = 0; i < PED->patchCoeff->rows; i++)
      for (j = 0; j < PED->patchCoeff->cols; j++)
        mp_to_rat(PED->patchCoeff_rat[i][j], &PED->patchCoeff->entry[i][j]);
  }
  else
  { // generate randomly based on variable groups - should only be used if using "true" MP tracking
    preproc_data *PPD = (preproc_data *)ptr1;

    // one patch for each variable group
    PED->num_patches = PPD->num_hom_var_gp + PPD->num_var_gp;

    // calculate the total number of variables that are needed
    total_vars = 0;
    for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
      total_vars += PPD->size[i] + PPD->type[i];

    init_mat_mp2(PED->patchCoeff, PED->num_patches, total_vars, prec);
    init_mat_rat(PED->patchCoeff_rat, PED->num_patches, total_vars);
    PED->patchCoeff->rows = PED->num_patches;
    PED->patchCoeff->cols = total_vars;

    // generate a patch for each of these variable groups
    var_gp_count = 0; // the number of variable_groups that we have been through
    curr_loc = PPD->num_var_gp; // current location of the variables
    for (i = 0; i < PED->num_patches; i++)
    { // zero out the row
      for (j = 0; j < total_vars; j++)
      {
        set_zero_mp(&PED->patchCoeff->entry[i][j]);
      }
      // put random numbers in the appropriate places
      if (PPD->type[i] == 1)
      { // needed for the extra homogenizing variable that was added
        get_comp_rand_mp(&PED->patchCoeff->entry[i][var_gp_count]);
        var_gp_count++;
      }
      for (j = 0; j < PPD->size[i]; j++)
      {
        get_comp_rand_mp(&PED->patchCoeff->entry[i][curr_loc]);
        curr_loc++;
      }
    }

    for (i = 0; i < PED->num_patches; i++)
      for (j = 0; j < total_vars; j++)
        mp_to_rat(PED->patchCoeff_rat[i][j], &PED->patchCoeff->entry[i][j]);
  }

  return;
}

void setupCoeffInSS_mp(char degreeFile[], int *P, preproc_data *PPD, start_system_eval_data_mp *SSED, int digits, int prec)
{ // setup the mhom start system - i.e. setup the coeffs in the start system

// SINCE THIS IS ONLY USED IN MHOM, the total number of variables described in PPD is the total number of variables used throughout

  int i, j, k, total_vars, total_deg, curr_loc, var_gp_count;
  int **mhomDeg, *linearGroup, **varGroupArray;
  comp_d tempComp;
  FILE *degIN = fopen(degreeFile, "r");

  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // calculate the total number of variables that are needed
  total_vars = 0;
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
    total_vars += PPD->size[i] + PPD->type[i];

  SSED->coeff_cols = total_vars;

  mhomDeg = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
    mhomDeg[i] = (int *)bmalloc((PPD->num_hom_var_gp + PPD->num_var_gp) * sizeof(int));

  // read in the mhom degrees - and find the total number of linears that will be needed
  total_deg = 0;
  for (i = 0; i < SSED->size_r; i++)
  {
    for (j = 0; j < (PPD->num_hom_var_gp + PPD->num_var_gp); j++)
    {
      fscanf(degIN, "%d\n", &mhomDeg[i][j]);
      total_deg += mhomDeg[i][j];
    }
    fscanf(degIN, "\n"); // extra
  }

  // convert to a list of the variable groups that are needed for the linears - use the Permuation P to get the list correct!
  linearGroup = (int *)bmalloc(total_deg * sizeof(int));
  curr_loc = 0;
  for (i = 0; i < SSED->size_r; i++)
    for (j = 0; j < (PPD->num_hom_var_gp + PPD->num_var_gp); j++)
      for (k = 0; k < mhomDeg[P[i]][j]; k++)
      {
        linearGroup[curr_loc] = j;
        curr_loc++;
      }

  // setup an array to tell which variables are in which group
  varGroupArray = (int **)bmalloc((PPD->num_hom_var_gp + PPD->num_var_gp) * sizeof(int *));
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
    varGroupArray[i] = (int *)bmalloc(total_vars * sizeof(int));

  var_gp_count = 0;
  curr_loc = PPD->num_var_gp;
  for (i = 0; i < (PPD->num_hom_var_gp + PPD->num_var_gp); i++)
  {
    for (j = 0; j < total_vars; j++)
      // initialize to zero
      varGroupArray[i][j] = 0;

    // put 1 in the appropriate places
    if (PPD->type[i] == 1)
    { // needed for the extra homogenizing variable that was added
      varGroupArray[i][var_gp_count] = 1;
      var_gp_count++;
    }
    for (j = 0; j < PPD->size[i]; j++)
    {
      varGroupArray[i][curr_loc] = 1;
      curr_loc++;
    }
  }

  // setup coeff
  SSED->coeff = (comp_mp **)bmalloc(total_deg * sizeof(comp_mp *));
  SSED->coeff_rat = (mpq_t ***)bmalloc(total_deg * sizeof(mpq_t **));
  for (i = 0; i < total_deg; i++)
  {
    SSED->coeff[i] = (comp_mp *)bmalloc(total_vars * sizeof(comp_mp));
    SSED->coeff_rat[i] = (mpq_t **)bmalloc(total_vars * sizeof(mpq_t *));

    // at linear i, the variable group that this corresponds to is linearGroup[i]
    for (j = 0; j < total_vars; j++)
    {
      SSED->coeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); // real & imaginary parts
      if (varGroupArray[linearGroup[i]][j] == 1)
      { // setup coeff & coeff_rat to random number
        get_comp_rand_rat(tempComp, SSED->coeff[i][j], SSED->coeff_rat[i][j], prec, prec, 1, 1);
      }
      else
      { // setup coeff & coeff_rat to 0
        set_zero_d(tempComp);
        comp_d_to_mp_rat(SSED->coeff[i][j], SSED->coeff_rat[i][j], tempComp, digits, prec, 1, 1);
      }
    }
  }

  // close files & release memory
  fclose(degIN);

  for (i = (PPD->num_hom_var_gp + PPD->num_var_gp - 1); i >= 0; i--)
    free(varGroupArray[i]);
  free(varGroupArray);

  for (i = SSED->size_r - 1; i >= 0; i--)
    free(mhomDeg[i]);
  free(mhomDeg);

  free(linearGroup);

  return;
}

void setupStartSystem_mp(int SSType, int size, int *deg, int *P, preproc_data *PPD, char degreeFile[], start_system_eval_data_mp *SSED, int digits, int prec)
{ // if SSType == 0, we are setuping up for total degree, else, for mhom

  int i;
  comp_d tempGamma;

  SSED->curr_prec = prec;

  // setup the startSystem
  SSED->startSystemType = SSType;
  SSED->size_r = size;
  SSED->degrees = (int *)bmalloc(size * sizeof(int));

  SSED->gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); // for real & imaginary
  get_comp_rand_rat(tempGamma, SSED->gamma, SSED->gamma_rat, prec, prec, 1, 1);

  SSED->max_degree = -1;
  for (i = 0; i < size; i++)
  {
    SSED->degrees[i] = deg[i];
    if (deg[i] > SSED->max_degree)
      SSED->max_degree = deg[i];
  }

  if (SSType == 1) // need to setup the coeff
    setupCoeffInSS_mp(degreeFile, P, PPD, SSED, digits, prec);

  return;
}

void setupSquareSystem_mp(prog_t *Prog, int finalSize, preproc_data *PPD, char degreeFile[], square_system_eval_data_mp *SSED, int digits, int prec, int adjustDegrees)
{ // setup the square system - finalSize is the final size of the 'square system' - without the patch

  int i, j, total_vars, max, max_loc = 0;
  int *tempInts = NULL;
  mat_d tempMat, tempMat2;

  init_mat_d(tempMat, 0, 0);
  init_mat_d(tempMat2, 0, 0);

  FILE *degIN = fopen(degreeFile, "r");
  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeFile);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // calcualate total number of variables
  total_vars = 0;
  for (i = 0; i < (PPD->num_var_gp + PPD->num_hom_var_gp); i++)
    total_vars += PPD->size[i] + PPD->type[i];

  SSED->Prog = Prog;
  SSED->size_f = PPD->num_funcs;
  SSED->size_r = finalSize;
  SSED->curr_prec = prec;
  SSED->noChanges = 0; // initialize to 0

  SSED->P = (int *)bmalloc(SSED->size_f * sizeof(int));
  SSED->W = (int **)bmalloc(SSED->size_r * sizeof(int *));
  for (i = 0; i < SSED->size_r; i++)
    SSED->W[i] = (int *)bmalloc((SSED->size_f - SSED->size_r) * sizeof(int));

  if (PPD->num_hom_var_gp + PPD->num_var_gp == 1)
  { // 1-hom case!

    // setup the original degrees
    SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
    for (i = 0; i < SSED->size_f; i++)
      fscanf(degIN, "%d\n\n", &SSED->orig_degrees[i]);

    // find B & B_perp
    if (total_vars == SSED->size_r + 1)
    { // setup B
      init_mat_mp(SSED->B, total_vars, total_vars);
      change_size_mat_d(tempMat, total_vars, total_vars);
      SSED->B->rows = SSED->B->cols = total_vars; // number of x variables

      make_matrix_ID_d(tempMat, total_vars, total_vars);  // if we are not changing the number of variables, then we can stick with the same ones!!

      // setup the matrix B
      init_mat_rat(SSED->B_rat, total_vars, total_vars);
      mat_d_to_mp_rat(SSED->B, SSED->B_rat, tempMat, digits, prec, 0, 0); // store tempMat in B & B_rat

      // setup the matrix B_perp
      change_size_mat_d(tempMat, total_vars, 0);
      init_mat_mp(SSED->B_perp, total_vars, 0);
      init_mat_rat(SSED->B_perp_rat, total_vars, 0);
      SSED->B_perp->rows = tempMat->rows = total_vars;
      SSED->B_perp->cols = tempMat->cols = 0;

      mat_d_to_mp_rat(SSED->B_perp, SSED->B_perp_rat, tempMat, digits, prec, 0, 0); // this won't do anything!
    }
    else
    {
      change_size_mat_d(tempMat, total_vars, total_vars);
      make_matrix_random_d(tempMat, total_vars, total_vars);

      // put to correct sizes
      change_size_mat_d(tempMat2, total_vars, total_vars - (SSED->size_r + 1));
      init_mat_mp(SSED->B, total_vars, SSED->size_r + 1);
      init_mat_mp(SSED->B_perp, total_vars, total_vars - (SSED->size_r + 1));
      SSED->B->rows = tempMat->rows = SSED->B_perp->rows = tempMat2->rows = total_vars;
      SSED->B->cols = tempMat->cols = SSED->size_r + 1;
      SSED->B_perp->cols = tempMat2->cols = total_vars - (SSED->size_r + 1);

      for (i = 0; i < tempMat2->rows; i++)
        for (j = 0; j < tempMat2->cols; j++)
        {
          set_d(&tempMat2->entry[i][j], &tempMat->entry[i][j+tempMat->cols]); // store the perpendicular part to tempMat2
        }

      // now we can just copy tempMat to B & tempMat2 to B_perp
      init_mat_rat(SSED->B_rat, SSED->B->rows, SSED->B->cols);
      mat_d_to_mp_rat(SSED->B, SSED->B_rat, tempMat, digits, prec, 0, 0);

      init_mat_rat(SSED->B_perp_rat, SSED->B_perp->rows, SSED->B_perp->cols);
      mat_d_to_mp_rat(SSED->B_perp, SSED->B_perp_rat, tempMat2, digits, prec, 0, 0);
    }

    if (!adjustDegrees && SSED->size_f == SSED->size_r)
    { // the calling function does not want the degrees adjusted and so since it is square, there is no need to adjust list based on degree
      SSED->noChanges = 1;
      for (i = 0; i < SSED->size_f; i++)
        SSED->P[i] = i;
    }
    else
    { // to create P, we need to find the location of the ith largest degree and store as P[i]
      tempInts = (int *)bmalloc(SSED->size_f * sizeof(int));
      for (i = 0; i < SSED->size_f; i++)
        tempInts[i] = 1; // will be = 0 when the function has been used

      for (i = 0; i < SSED->size_f; i++)
      {
        max = -1;
        for (j = 0; j < SSED->size_f; j++)
        {
          if ((max < SSED->orig_degrees[j]) && (tempInts[j] == 1))
          {
            max = SSED->orig_degrees[j];
            max_loc = j;
          }
        }
        SSED->P[i] = max_loc;
        tempInts[max_loc] = 0;
      }
    }

    // setup new_degrees
    SSED->new_degrees = (int *)bmalloc(SSED->size_r * sizeof(int));
    for (i = 0; i < SSED->size_r; i++)
      SSED->new_degrees[i] = SSED->orig_degrees[SSED->P[i]];

    // create A
    change_size_mat_d(tempMat, SSED->size_r, SSED->size_f - SSED->size_r);
    make_matrix_random_d(tempMat, SSED->size_r, SSED->size_f - SSED->size_r);

    init_mat_mp2(SSED->A, tempMat->rows, tempMat->cols, prec);
    init_mat_rat(SSED->A_rat, tempMat->rows, tempMat->cols);
    mat_d_to_mp_rat(SSED->A, SSED->A_rat, tempMat, digits, prec, 0, 0);

    // to setup W, we need to find the difference in the exponents
    for (i = 0; i < SSED->size_r; i++)
      for (j = 0; j < SSED->size_f - SSED->size_r; j++)
        SSED->W[i][j] = SSED->orig_degrees[SSED->P[i]] - SSED->orig_degrees[SSED->P[j+SSED->size_r]];

    // to find max_of_W, we simply move through W and find the maximum
    SSED->max_of_W = 0;
    for (i = 0; i < SSED->size_r; i++)
      for (j = 0; j < SSED->size_f - SSED->size_r; j++)
        if (SSED->max_of_W < SSED->W[i][j])
          SSED->max_of_W = SSED->W[i][j];
  }
  else if (PPD->num_hom_var_gp + PPD->num_var_gp > 1)
  { // mhom case
    if (finalSize != SSED->size_f)
    {
      printf("setupSquareSystem: finalSize (%d) != size_f (%d) in mhom case!\n", finalSize, SSED->size_f);
      bexit(ERROR_INVALID_SIZE);
    }
    // with mhom, there are no changes are done
    SSED->noChanges = 1;

    // change of variable matrix is identity
    change_size_mat_d(tempMat, total_vars, total_vars);
    make_matrix_ID_d(tempMat, total_vars, total_vars); // if we are not changing the number of variables, then we can stick with the same ones!!

    // setup the matrix B
    init_mat_mp2(SSED->B, total_vars, total_vars, prec);
    init_mat_rat(SSED->B_rat, total_vars, total_vars);
    mat_d_to_mp_rat(SSED->B, SSED->B_rat, tempMat, digits, prec, 0, 0);

    // B_perp has no columns!!
    change_size_mat_d(tempMat, total_vars, 0);
    tempMat->rows = total_vars;
    tempMat->cols = 0;

    init_mat_mp2(SSED->B_perp, total_vars, 0, prec);
    init_mat_rat(SSED->B_perp_rat, total_vars, 0);
    mat_d_to_mp_rat(SSED->B_perp, SSED->B_perp_rat, tempMat, digits, prec, 0, 0); // this won't do anything!

    // permutation is the identity
    for (i = 0; i < SSED->size_f; i++)
      SSED->P[i] = i;

    // A has no columns!!
    change_size_mat_d(tempMat, SSED->size_r, 0);
    tempMat->rows = SSED->size_r;
    tempMat->cols = 0;

    init_mat_mp2(SSED->A, SSED->size_r, 0, prec);
    init_mat_rat(SSED->A_rat, SSED->size_r, 0); 
    mat_d_to_mp_rat(SSED->A, SSED->A_rat, tempMat, digits, prec, 0, 0);

    // setup orig_degrees to be the total of the mhom degrees
    SSED->orig_degrees = (int *)bmalloc(SSED->size_f * sizeof(int));
    SSED->new_degrees  = (int *)bmalloc(SSED->size_f * sizeof(int));
    for (i = 0; i < SSED->size_f; i++)
    {
      SSED->orig_degrees[i] = 0;
      for (j = 0; j < PPD->num_hom_var_gp + PPD->num_var_gp; j++)
      {
        fscanf(degIN, "%d\n", &max);
        SSED->orig_degrees[i] += max;
      }
      fscanf(degIN, "\n");
      SSED->new_degrees[i] = SSED->orig_degrees[i];
    }

    // W has no columns
    SSED->max_of_W = 0;
  }
  else
  {
    printf("ERROR: In setupSquareSystem, there are %d variables groups!\n", PPD->num_hom_var_gp + PPD->num_var_gp);
    bexit(ERROR_INVALID_SIZE);
  }

  free(tempInts);
  clear_mat_d(tempMat);
  clear_mat_d(tempMat2);

  return;
}

void setupBasicEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg, int squareSize, int patchType, int ssType, int prec, void const *ptr1, void const *ptr2, basic_eval_data_mp *BED, int adjustDegrees)
{
  int i, digits = prec_to_digits(mpf_get_default_prec());

  setupPreProcData(preprocFile, &BED->preProcData);
  setupSquareSystem_mp(dummyProg, squareSize, &BED->preProcData, degreeFile, &BED->squareSystem, digits, prec, adjustDegrees); 
  // NOTE: squareSystem must be setup before the patch!!!
  setupPatch_mp(patchType, &BED->patch, digits, prec, ptr1, ptr2);

  // find the degrees of the 'square system'
  int *deg = (int *)bmalloc(BED->squareSystem.size_r * sizeof(int));
  for (i = 0; i < BED->squareSystem.size_r; i++)
    deg[i] = BED->squareSystem.orig_degrees[BED->squareSystem.P[i]];

  setupStartSystem_mp(ssType, BED->squareSystem.size_r, deg, BED->squareSystem.P, &BED->preProcData, degreeFile, &BED->startSystem, digits, prec);

  free(deg);

  return;
}

void setupBEDUsingUserHom_mp(prog_t *dummyProg, int MPType, basic_eval_data_mp *BED)
{ // set everything to nothing expect for the s.l.p.!!
  comp_d tempGamma;
  int prec = dummyProg->precision; // precision is the one used by Prog

  // patch things
  init_mat_mp(BED->patch.patchCoeff, 0, 0);
  BED->patch.num_patches = BED->patch.patchCoeff->rows = BED->patch.patchCoeff->cols = 0;
  BED->patch.patchCoeff_rat = NULL;
  BED->patch.curr_prec = prec;

  // start system things
  BED->startSystem.startSystemType = -1;
  BED->startSystem.degrees = NULL;
  BED->startSystem.size_r = BED->startSystem.max_degree = BED->startSystem.coeff_cols = 0;
  BED->startSystem.gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); // for real & imaginary
  get_comp_rand_rat(tempGamma, BED->startSystem.gamma, BED->startSystem.gamma_rat, prec, prec, 1, 1);
  BED->startSystem.curr_prec = prec;

  // square system things
  BED->squareSystem.Prog = dummyProg;
  BED->squareSystem.orig_degrees = BED->squareSystem.new_degrees = BED->squareSystem.P = NULL;
  BED->squareSystem.W = NULL;
  init_mat_mp(BED->squareSystem.B, 0, 0);
  BED->squareSystem.size_f = BED->squareSystem.B->rows = BED->squareSystem.B->cols = 0;
  init_mat_mp(BED->squareSystem.B_perp, 0, 0);
  BED->squareSystem.B_perp->rows = BED->squareSystem.B_perp->cols = 0;
  init_mat_mp(BED->squareSystem.A, 0, 0);
  BED->squareSystem.max_of_W = BED->squareSystem.A->rows = BED->squareSystem.A->cols = 0;
  BED->squareSystem.B_rat = BED->squareSystem.B_perp_rat = BED->squareSystem.A_rat = NULL;
  BED->squareSystem.size_r = 0;
  BED->squareSystem.curr_prec = prec;

  // preproc data things
  BED->preProcData.num_funcs = BED->preProcData.num_hom_var_gp = BED->preProcData.num_var_gp = 0;
  BED->preProcData.type = BED->preProcData.size = NULL;

  return;
}

void setupBEDUsingParamHom_mp(prog_t *dummyProg, char *preprocFile, int MPType, basic_eval_data_mp *BED)
{ // setup SLP & patch
  comp_d tempGamma;
  int prec = dummyProg->precision; // precision is the one used by Prog
  int digits = prec_to_digits(prec); // digits

  // setup preProcData
  setupPreProcData(preprocFile, &BED->preProcData);

  // setup patch
  setupPatch_mp(0, &BED->patch, digits, prec, &BED->preProcData, NULL);

  // start system things
  BED->startSystem.startSystemType = -1;
  BED->startSystem.degrees = NULL;
  BED->startSystem.size_r = BED->startSystem.max_degree = BED->startSystem.coeff_cols = 0;
  BED->startSystem.gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t)); // for real & imaginary
  get_comp_rand_rat(tempGamma, BED->startSystem.gamma, BED->startSystem.gamma_rat, prec, prec, 1, 1);
  BED->startSystem.curr_prec = prec;

  // square system things
  BED->squareSystem.Prog = dummyProg;
  BED->squareSystem.orig_degrees = BED->squareSystem.new_degrees = BED->squareSystem.P = NULL;
  BED->squareSystem.W = NULL;
  init_mat_mp(BED->squareSystem.B, 0, 0);
  BED->squareSystem.size_f = BED->squareSystem.B->rows = BED->squareSystem.B->cols = 0;
  init_mat_mp(BED->squareSystem.B_perp, 0, 0);
  BED->squareSystem.B_perp->rows = BED->squareSystem.B_perp->cols = 0;
  init_mat_mp(BED->squareSystem.A, 0, 0);
  BED->squareSystem.max_of_W = BED->squareSystem.A->rows = BED->squareSystem.A->cols = 0;
  BED->squareSystem.B_rat = BED->squareSystem.B_perp_rat = BED->squareSystem.A_rat = NULL;
  BED->squareSystem.size_r = 0;
  BED->squareSystem.curr_prec = prec;

  return;
}

void setupTD_startPoints_mp(char pointsIN[], char pointsOUT[], int size, int *degs, patch_eval_data_mp *PED)
{ // reads points from pointsIN and moves them to points on the patch and writes them to pointsOUT
  // ASSUME ONLY 1 PATCH!!!!!!!

  int i, j, num_points;
  FILE *IN, *OUT;
  comp_mp tempComp;
  vec_mp patchVals, tempVec;

  // initialize MP
  init_mp(tempComp);
  init_vec_mp(patchVals, 0);
  init_vec_mp(tempVec, size + 1);
  tempVec->size = size + 1; // extra coord at the beginning since patch depends on r + 1 variables

  // create the start points
  TDstartMaker_mp(degs, size);

  // open the files
  IN = fopen(pointsIN, "r");
  OUT = fopen(pointsOUT, "w");
  if (IN == NULL)
  {
    printf("ERROR: '%s' does not exist!!!\n\n", pointsIN);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in the number of points
  fscanf(IN, "%d\n\n", &num_points);
  // write the number of points
  fprintf(OUT, "%d\n\n", num_points);

  for (i = 0; i < num_points; i++)
  { // read in the ith point
    set_zero_mp(&tempVec->coord[0]);
    mpf_set_ui(tempVec->coord[0].r, 1); // hom coord is 1
    for (j = 1; j < tempVec->size; j++)
    { // read in the next point
      mpf_inp_str(tempVec->coord[j].r, IN, 10);
      mpf_inp_str(tempVec->coord[j].i, IN, 10);
      fscanf(IN, ";\n");
    }
    fscanf(IN, "\n");

    // calculate the value of the patch at this point
    mul_mat_vec_mp(patchVals, PED->patchCoeff, tempVec);

    // adjust tempPoint based on patchVals[0] & write to OUT
    recip_mp(tempComp, &patchVals->coord[0]);
    for (j = 0; j < tempVec->size; j++)
    {
      mul_mp(&tempVec->coord[j], &tempVec->coord[j], tempComp);
      print_mp(OUT, 0, &tempVec->coord[j]);
      fprintf(OUT, ";\n");
    }
    fprintf(OUT, "\n");
  }

  // close the files
  fclose(IN);
  fclose(OUT);

  // clear MP
  clear_mp(tempComp);
  clear_vec_mp(patchVals);
  clear_vec_mp(tempVec);

  return;
}

