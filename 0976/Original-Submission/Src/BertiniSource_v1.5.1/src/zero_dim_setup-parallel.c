// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

int setupArr(prog_t *P, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow);

int setupProg(prog_t *P, int precision, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: numvars - number of variables                  *
* NOTES:                                                        *
\***************************************************************/
{ // just call setupProg_count
  int i, numVars, *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;

  // setup prog
  numVars = setupProg_count(P, precision, MPType, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow);

  // clear memory
  if (P->numSubfuncs > 0)
  { // clear subFuncsBelow
    for (i = P->numFuncs - 1; i >= 0; i--)
      free(subFuncsBelow[i]);
    free(subFuncsBelow);
  }
  free(startSub);
  free(endSub);
  free(startFunc);
  free(endFunc);
  free(startJvsub);
  free(endJvsub);
  free(startJv);
  free(endJv);

  return numVars;
}

int setupProg_count(prog_t *P, int precision, int MPType, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: numvars - number of variables                  *
* NOTES:                                                        *
\***************************************************************/
{ 
  int numvars = 0;

  // initialize MP to the desired precision & setup the basic things to evaluate the program
  initMP(precision);
  initEvalProg(MPType);

  // setup precision in P
  P->precision = precision;

  // setup the part in P dealing with arr.out
  numvars = setupArr(P, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

  // setup the nums in P
  setupNums(&P->nums, P->numNums, precision, MPType);

  return numvars;
}  

int setupArr(prog_t *P, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: numvars - number of variables                  *
* NOTES:                                                        *
\***************************************************************/
{
  FILE *arrIN = NULL;
  int i, j, num_var_gps, workspaceSize, numvars, varAddr, numpathvars, pathvarAddr, numpars;
  int parAddr, parderAddr, numfuncs, funcAddr, jacVAddr, jacPAddr;
  int numconsts, constAddr, numNums, numAddr, numsubfuncs, subfuncAddr;
  int subfuncDerivWRTVarsStart, subfuncDerivWRTParamsStart, rand_index, numInst, IAddr;
  char ch;

  // open arr.out
  arrIN = fopen("arr.out", "r");
  if (arrIN == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n", "arr.out");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // move to the 'X' that is after the instructions
  do
  {
    ch = fgetc(arrIN);
  } while (ch != 'X');
  fscanf(arrIN, "\n");

  // read in the postamble
  fscanf(arrIN, "TotalRoomNeeded %d;\n", &workspaceSize);
  fscanf(arrIN, "NVAR %d %d;\n", &numvars, &varAddr);
  fscanf(arrIN, "NPATHVAR %d %d;\n", &numpathvars, &pathvarAddr);
  fscanf(arrIN, "NPAR %d %d %d;\n", &numpars, &parAddr, &parderAddr);
  fscanf(arrIN, "NFCN %d %d %d %d;\n", &numfuncs, &funcAddr, &jacVAddr, &jacPAddr);
  fscanf(arrIN, "NCON %d %d;\n", &numconsts, &constAddr);
  fscanf(arrIN, "NNUM %d %d;\n", &numNums, &numAddr);
  fscanf(arrIN, "SUBFCN %d %d %d %d;\n", &numsubfuncs, &subfuncAddr, &subfuncDerivWRTVarsStart, &subfuncDerivWRTParamsStart);
  fscanf(arrIN, "NUMINST %d;\n", &numInst);
  fscanf(arrIN, "CMPLX %d;\n", &IAddr);
  fscanf(arrIN, "VARGPS %d;\n", &num_var_gps);
  // setup number of variable groups in P
  P->num_var_gps = num_var_gps;
  P->var_gp_sizes = (int *)bcalloc(num_var_gps, sizeof(int));
  for (i = 0; i < num_var_gps; i++)
    fscanf(arrIN, " %d", &P->var_gp_sizes[i]);
  fscanf(arrIN, ";\n");
  fscanf(arrIN, "RANDINDEX %d;\n", &rand_index);

  // read in INSTCOUNT - endUpdate, endParams, endFnEval, endPDeriv, endJvEval
  fscanf(arrIN, "INSTCOUNT %d %d %d %d %d;\n", &P->numInstAtEndUpdate, &P->numInstAtEndParams, &P->numInstAtEndFnEval, &P->numInstAtEndPDeriv, &P->numInstAtEndJvEval);
  // next line contains start/end of subfunctions/functions and their derivs
  *startSub = (int *)bmalloc(numsubfuncs * sizeof(int));
  *endSub = (int *)bmalloc(numsubfuncs * sizeof(int));
  *startFunc = (int *)bmalloc(numfuncs * sizeof(int));
  *endFunc = (int *)bmalloc(numfuncs * sizeof(int));
  *startJvsub = (int *)bmalloc(numsubfuncs * sizeof(int));
  *endJvsub = (int *)bmalloc(numsubfuncs * sizeof(int));
  *startJv = (int *)bmalloc(numfuncs * sizeof(int));
  *endJv = (int *)bmalloc(numfuncs * sizeof(int));
  for (i = 0; i < numsubfuncs; i++)
    fscanf(arrIN, "%d %d", &(*startSub)[i], &(*endSub)[i]);
  for (i = 0; i < numfuncs; i++)
    fscanf(arrIN, "%d %d", &(*startFunc)[i], &(*endFunc)[i]);
  for (i = 0; i < numsubfuncs; i++)
    fscanf(arrIN, "%d %d", &(*startJvsub)[i], &(*endJvsub)[i]);
  for (i = 0; i < numfuncs; i++)
    fscanf(arrIN, "%d %d", &(*startJv)[i], &(*endJv)[i]);
  scanRestOfLine(arrIN);

  // next is the 'matrix' of which subfunctions are below each function
  if (numsubfuncs > 0)
  { // read in the matrix
    *subFuncsBelow = (int **)bmalloc(numfuncs * sizeof(int *));
    for (i = 0; i < numfuncs; i++)
      (*subFuncsBelow)[i] = (int *)bmalloc(numsubfuncs * sizeof(int));

    for (i = 0; i < numfuncs; i++)
    { // read in the values
      for (j = 0; j < numsubfuncs; j++)
        fscanf(arrIN, "%d", &(*subFuncsBelow)[i][j]);
      scanRestOfLine(arrIN);
    }
  }
  else
  { // set to NULL
    *subFuncsBelow = NULL;
  }

  // Now we read in the instructions from arr.out:
  rewind(arrIN);
  P->prog = (int *)bcalloc(numInst, sizeof(int));
  for (i = 0; i < numInst; i++)
    fscanf(arrIN, "%d", &P->prog[i]);

  //Finally, we set up the Prog in T:
    //Basic sizes.
  P->size = numInst;
  P->memSize = workspaceSize;
    //# of various types of names.
  P->numVars = numvars;
  P->numPathVars = numpathvars;  //This, of course, should always be 1.
  P->numPars = numpars;
  P->numFuncs = numfuncs;
  P->numNums = numNums;
  P->numConsts = numconsts;
  P->numSubfuncs = numsubfuncs;
    //Locations of results.
  P->inpVars = varAddr;
  P->inpPathVars = pathvarAddr;
  P->evalPars = parAddr;
  P->numAddr = numAddr;
  P->constAddr = constAddr;
  P->evalDPars = parderAddr;
  P->evalFuncs = funcAddr;
  P->evalJVars = jacVAddr;
  P->evalJPars = jacPAddr;
  P->evalSubs = subfuncAddr;
  P->evalJSubsV = subfuncDerivWRTVarsStart;
  P->evalJSubsP = subfuncDerivWRTParamsStart;
  P->IAddr = IAddr;
  P->index_of_first_number_for_proj_trans = rand_index;

  // close arrIN
  fclose(arrIN);

  return numvars;
}

void setupNums(num_t **nums, int numNums, int precision, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i;
  
  // open num.out
  FILE *numIN = fopen("num.out", "r");
  if (numIN == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n", "num.out");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read the numbers of num.out into nums:
  *nums = (num_t *)bmalloc(numNums * sizeof(num_t));
  for (i = 0; i < numNums; i++)
  { // initialize
    mpq_init((*nums)[i].rat);
    // read in the value
    mpq_inp_str((*nums)[i].rat, numIN, 10);
    scanRestOfLine(numIN);
    // remove common factors
    mpq_canonicalize((*nums)[i].rat);
    // make the first approximation here:
    (*nums)[i].currPrec = precision;
    mpf_init2((*nums)[i].real, (*nums)[i].currPrec);
    mpf_set_q((*nums)[i].real, (*nums)[i].rat);
  }

  // close the file
  fclose(numIN);

  return;
}

int userHom_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for user defined homotopy tracking               *
\***************************************************************/
{ // the user defined a homotopy and we will simply track it
  int numOrigVars;

  *eval_d = &userHom_eval_d;
  *eval_mp = &userHom_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  if (T->MPType == 2) // using AMP - need to allocate space to store BED_mp
    ED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
  else
    ED->BED_mp = NULL;

  T->numVars = numOrigVars = setupProg(dummyProg, T->Precision, T->MPType);  // setup a straight-line program, using the file(s) created by the parser.

  setupBEDUsingUserHom_d(dummyProg, T->MPType, ED); // setup the main structure - ED

  return numOrigVars;
}

int paramHom_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, int moveStartPts, char *startName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for a parameter homotopy tracking                *
\***************************************************************/
{ // a parameter homotopy setup in the SLP with added patches
  int i, j, numOrigVars = 0, numVars = 0, numGps = 0, numPts = 0;

  *eval_d = &paramHom_eval_d;
  *eval_mp = &paramHom_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  if (T->MPType == 2) // using AMP - need to allocate space to store BED_mp
    ED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
  else
    ED->BED_mp = NULL;

  T->numVars = numOrigVars = setupProg(dummyProg, T->Precision, T->MPType);  // setup a straight-line program, using the file(s) created by the parser.

  setupBEDUsingParamHom_d(dummyProg, preprocFile, T->MPType, ED); // setup the main structure - ED

  if (moveStartPts)
  { // move the start points to the patch
    point_d inPt, outPt;
    FILE *INPTS = NULL, *OUTPTS = NULL;

    // check that name is different
    if (strcmp(startName, "start_param_hom") == 0)
    {
      printf("\nERROR: Due to your settings, the start file can not be named 'start_param_hom'!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // setup files
    INPTS = fopen(startName, "r");
    OUTPTS = fopen("start_param_hom", "w");

    if (INPTS == NULL)
    {
      printf("\n\nERROR: '%s' does not exist!!!\n\n", startName);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // count the number of variables for the points in start
    numGps = ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp;
    for (i = 0; i < numGps; i++)
      numVars += ED->preProcData.size[i];

    init_point_d(inPt, numOrigVars);
    init_point_d(outPt, numOrigVars);
    inPt->size = outPt->size = numOrigVars;

    // read in the number of points and print
    fscanf(INPTS, "%d", &numPts);
    fprintf(OUTPTS, "%d\n\n", numPts);
  
    // more the points to the patch
    for (i = 0; i < numPts; i++)
    { // setup point
      for (j = 0; j < numOrigVars - numVars; j++)
      { // hom coordinate is 1
        set_one_d(&inPt->coord[j]);
      }
      for (j = numOrigVars - numVars; j < numOrigVars; j++)
      {
        fscanf(INPTS, "%lf%lf", &inPt->coord[j].r, &inPt->coord[j].i);
        scanRestOfLine(INPTS);
      }
      // move to patch
      move_to_patch_mat_d(outPt, inPt, ED->patch.patchCoeff, &ED->preProcData);
      for (j = 0; j < numOrigVars; j++)
      {
        print_d(OUTPTS, 0, &outPt->coord[j]);
        fprintf(OUTPTS, "\n");
     }
      fprintf(OUTPTS, "\n");
    }
  
    fclose(INPTS);
    fclose(OUTPTS);

    clear_point_d(inPt);
    clear_point_d(outPt);
  }

  return numOrigVars;
}

int zero_dim_basic_setup_d(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_d *ED, prog_t *dummyProg, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow, int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, int findStartPts, char *pointsIN, char *pointsOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of original variables                   *
* NOTES: setup for zero dimensional tracking                    *
\***************************************************************/
{ // need to create the homotopy
  int rank, patchType, ssType, numOrigVars, adjustDegrees, numGps;

  *eval_d = &basic_eval_d;
  *eval_mp = &basic_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  if (T->MPType == 2) // using AMP - need to allocate space to store BED_mp
    ED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
  else
    ED->BED_mp = NULL;

  // setup a straight-line program, using the file(s) created by the parser
  T->numVars = numOrigVars = setupProg_count(dummyProg, T->Precision, T->MPType, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

  // setup preProcData
  setupPreProcData(preprocFile, &ED->preProcData);
  numGps = ED->preProcData.num_var_gp + ED->preProcData.num_hom_var_gp;

  // find the rank
  rank = rank_finder_d(&ED->preProcData, dummyProg, T, T->numVars);

  // check to make sure that it is possible to have a zero dimensional component
  if (T->numVars > rank + numGps)
  {
    printf("The system has no zero dimensional solutions based on its rank!\n");
    printf("The rank of the system including the patches is %d while the total number of variables is %d.\n\n", rank + numGps, T->numVars);
    bexit(ERROR_INPUT_SYSTEM);
  }
  // adjust the number of variables based on the rank
  T->numVars = rank + numGps;

  // now that we know the rank, we can setup the rest of ED
  if (numGps == 1)
  { // 1-hom
    patchType = 2; // 1-hom patch
    ssType = 0;    // with 1-hom, we use total degree start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setupBasicEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &T->numVars, NULL, NULL, NULL, ED, adjustDegrees);
  }
  else
  { // m-hom, m > 1
    patchType = 0; // random patch based on m-hom variable structure
    ssType = 1;    // with m-hom, we use the mhom structure for start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setupBasicEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &ED->preProcData, NULL, NULL, NULL, ED, adjustDegrees);
  }

  // create the start points, if needed
  if (findStartPts)
  {
    if (ssType == 0)
    { // setup total degree start points
      setupTD_startPoints_d(pointsIN, pointsOUT, ED->startSystem.size_r, ED->startSystem.degrees, &ED->patch);
    }
    else
    { // setup m-hom start points
      MHstartMaker_d(&ED->preProcData, ED->squareSystem.P, ED->startSystem.coeff, ED->patch.patchCoeff);
    }
  }

  return numOrigVars;
}

/////////////////////// MP VERSIONS //////////////////////////////

int userHom_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for user defined homotopy tracking               *
\***************************************************************/
{ // the user defined a homotopy and we will simply track it
  int numOrigVars;

  *eval = &userHom_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  T->numVars = numOrigVars = setupProg(dummyProg, T->Precision, T->MPType);  // setup a straight-line program, using the file(s) created by the parser.

  setupBEDUsingUserHom_mp(dummyProg, T->MPType, ED); // setup the main structure - ED

  return numOrigVars;
}

int paramHom_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, int moveStartPts, char *startName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup for a parameter homotopy tracking                *
\***************************************************************/
{ // a parameter homotopy setup in the SLP with added patches
  int i, j, numOrigVars = 0, numVars = 0, numGps = 0, numPts = 0;

  *eval = &paramHom_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  T->numVars = numOrigVars = setupProg(dummyProg, T->Precision, T->MPType);  // setup a straight-line program, using the file(s) created by the parser.

  setupBEDUsingParamHom_mp(dummyProg, preprocFile, T->MPType, ED); // setup the main structure - ED

  if (moveStartPts)
  { // move the start points to the patch
    point_mp inPt, outPt;
    FILE *INPTS = NULL, *OUTPTS = NULL;

    // check that name is different
    if (strcmp(startName, "start_param_hom") == 0)
    {
      printf("\nERROR: Due to your settings, the start file can not be named 'start_param_hom'!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // setup files
    INPTS = fopen(startName, "r");
    OUTPTS = fopen("start_param_hom", "w");

    if (INPTS == NULL)
    {
      printf("\nERROR: '%s' does not exist!!!\n\n", startName);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // count the number of variables for the points in start
    numGps = ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp;
    for (i = 0; i < numGps; i++)
      numVars += ED->preProcData.size[i];

    init_point_mp(inPt, numOrigVars);
    init_point_mp(outPt, numOrigVars);
    inPt->size = outPt->size = numOrigVars;
 
    // read in the number of points and print
    fscanf(INPTS, "%d", &numPts);
    fprintf(OUTPTS, "%d\n\n", numPts);

    // more the points to the patch
    for (i = 0; i < numPts; i++)
    { // setup point
      for (j = 0; j < numOrigVars - numVars; j++)
      { // hom coordinate is 1
        set_one_mp(&inPt->coord[j]);
      }
      for (j = numOrigVars - numVars; j < numOrigVars; j++)
      {
        mpf_inp_str(inPt->coord[j].r, INPTS, 10);
        mpf_inp_str(inPt->coord[j].i, INPTS, 10);
        scanRestOfLine(INPTS);
      }
      // move to patch
      move_to_patch_mat_mp(outPt, inPt, ED->patch.patchCoeff, &ED->preProcData);
      for (j = 0; j < numOrigVars; j++)
      {
        print_mp(OUTPTS, 0, &outPt->coord[j]);
        fprintf(OUTPTS, "\n");
      }
      fprintf(OUTPTS, "\n");
    }

    fclose(INPTS);
    fclose(OUTPTS);

    clear_point_mp(inPt);
    clear_point_mp(outPt);
  }

  return numOrigVars;
}

int zero_dim_basic_setup_mp(FILE **OUT, char *outName, FILE **midOUT, char *midName, tracker_config_t *T, basic_eval_data_mp *ED, prog_t *dummyProg, int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow, int (**eval)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), char *preprocFile, char *degreeFile, int findStartPts, char *pointsIN, char *pointsOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of original variables                   *
* NOTES: setup for zero dimensional tracking                    *
\***************************************************************/
{ // need to create the homotopy
  int rank, patchType, ssType, numOrigVars, adjustDegrees;

  *eval = &basic_eval_mp;

  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");

  // setup a straight-line program, using the file(s) created by the parser
  T->numVars = numOrigVars = setupProg_count(dummyProg, T->Precision, T->MPType, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);

  // setup preProcData
  setupPreProcData("preproc_data", &ED->preProcData);

  // find the rank
  rank = rank_finder_mp(&ED->preProcData, dummyProg, T, T->numVars);

  // check to make sure that it is possible to have a zero dimensional component
  if (T->numVars > rank + ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp)
  {
    printf("The system has no zero dimensional solutions based on its rank!\n");
    printf("The rank of the system including the patches is %d while the total number of variables is %d.\n\n", rank + ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp, T->numVars);
    bexit(ERROR_INPUT_SYSTEM);
  }
  // adjust the number of variables based on the rank
  T->numVars = rank + ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp;

  // now that we know the rank, we can setup the rest of ED
  if (ED->preProcData.num_hom_var_gp + ED->preProcData.num_var_gp == 1)
  { // 1-hom
    patchType = 2; // 1-hom patch
    ssType = 0;    // with 1-hom, we use total degree start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setupBasicEval_mp(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->Precision, &T->numVars, NULL, ED, adjustDegrees);
  }
  else
  { // m-hom, m > 1
    patchType = 0; // random patch based on m-hom variable structure
    ssType = 1;    // with m-hom, we use the mhom structure for start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setupBasicEval_mp(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->Precision, &ED->preProcData, NULL, ED, adjustDegrees);
  }

  // create the start points, if needed
  if (findStartPts)
  {
    if (ssType == 0)
    { // setup total degree start points
      setupTD_startPoints_mp(pointsIN, pointsOUT, ED->startSystem.size_r, ED->startSystem.degrees, &ED->patch);
    }
    else
    { // setup m-hom start points
      MHstartMaker_mp(&ED->preProcData, ED->squareSystem.P, ED->startSystem.coeff, ED->patch.patchCoeff);
    }
  }

  return numOrigVars;
}

void printStartSystem(FILE *FP, int MPType, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the relevant data to FP so that we can recover  *
* the system if needed for sharpening                           *
\***************************************************************/
{
  int i, j, totalDeg;
  
  if (MPType == 0)
  { // print the start system type
    fprintf(FP, "%d\n", ED_d->startSystem.startSystemType);

    // print gamma
    print_d(FP, 0, ED_d->startSystem.gamma);
    fprintf(FP, "\n");

    if (ED_d->startSystem.startSystemType == 1)
    { // print the coefficients
      totalDeg = 0;
      for (i = 0; i < ED_d->startSystem.size_r; i++)
        totalDeg += ED_d->startSystem.degrees[i];

      fprintf(FP, "%d %d\n", totalDeg, ED_d->startSystem.coeff_cols);
      for (i = 0; i < totalDeg; i++)
        for (j = 0; j < ED_d->startSystem.coeff_cols; j++)
        {
          print_d(FP, 0, ED_d->startSystem.coeff[i][j]);
          fprintf(FP, "\n");
        }
    }
  }
  else
  { // print the start system type
    fprintf(FP, "%d\n", ED_mp->startSystem.startSystemType);

    // print gamma
    mpq_out_str(FP, 10, ED_mp->startSystem.gamma_rat[0]);
    fprintf(FP, " ");
    mpq_out_str(FP, 10, ED_mp->startSystem.gamma_rat[1]);
    fprintf(FP, "\n");

    if (ED_mp->startSystem.startSystemType == 1)
    { // print the coefficients
      totalDeg = 0;
      for (i = 0; i < ED_mp->startSystem.size_r; i++)
        totalDeg += ED_mp->startSystem.degrees[i];

      fprintf(FP, "%d %d\n", totalDeg, ED_mp->startSystem.coeff_cols);
      for (i = 0; i < totalDeg; i++)
        for (j = 0; j < ED_mp->startSystem.coeff_cols; j++)
        {
          mpq_out_str(FP, 10, ED_mp->startSystem.coeff_rat[i][j][0]);
          fprintf(FP, " ");
          mpq_out_str(FP, 10, ED_mp->startSystem.coeff_rat[i][j][1]);
          fprintf(FP, "\n");
        }
    }
  }

  return;
}

void printSquareSystem(FILE *FP, int MPType, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the relevant data to FP so that we can recover  *
* the system if needed for sharpening                           *
\***************************************************************/
{
  int i, j, rows, cols;

  if (MPType == 0)
  { // print B
    rows = ED_d->squareSystem.B->rows;
    cols = ED_d->squareSystem.B->cols;
    fprintf(FP, "%d %d\n", rows, cols); 
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_d(FP, 0, &ED_d->squareSystem.B->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");

    // print B_perp
    rows = ED_d->squareSystem.B_perp->rows;
    cols = ED_d->squareSystem.B_perp->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_d(FP, 0, &ED_d->squareSystem.B_perp->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
    
    // print A
    rows = ED_d->squareSystem.A->rows;
    cols = ED_d->squareSystem.A->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_d(FP, 0, &ED_d->squareSystem.A->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }
  else
  { // print B
    rows = ED_mp->squareSystem.B->rows;
    cols = ED_mp->squareSystem.B->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(FP, 10, ED_mp->squareSystem.B_rat[i][j][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, ED_mp->squareSystem.B_rat[i][j][1]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");

    // print B_perp
    rows = ED_mp->squareSystem.B_perp->rows;
    cols = ED_mp->squareSystem.B_perp->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(FP, 10, ED_mp->squareSystem.B_perp_rat[i][j][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, ED_mp->squareSystem.B_perp_rat[i][j][1]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");

    // print A
    rows = ED_mp->squareSystem.A->rows;
    cols = ED_mp->squareSystem.A->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(FP, 10, ED_mp->squareSystem.A_rat[i][j][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, ED_mp->squareSystem.A_rat[i][j][1]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }

  return;
}

void printZeroDimRelevantData(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the relevant data to FP so that we can recover  *
* the system if needed for sharpening                           *
\***************************************************************/
{
  // print the MPType and if an eq-by-eq method (diagona/regen) was used
  fprintf(FP, "%d %d\n", MPType, eqbyeqMethod);

  // print the patch
  printPatchCoeff(FP, MPType, ED_d, ED_mp);

  // print the start system
  printStartSystem(FP, MPType, ED_d, ED_mp);

  // print the square system
  printSquareSystem(FP, MPType, ED_d, ED_mp);

  return;
}

void readInStartSystem(FILE *IN, int MPType, int old_MPType, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads the start system from IN and stores to BED       *
\***************************************************************/
{
  int i, j, totalDeg, inputType;

  // read in the start system type
  fscanf(IN, "%d", &inputType);

  // verify we have agreement
  if ((MPType == 1 && (inputType != ED_mp->startSystem.startSystemType)) || ((MPType == 0 || MPType == 2) && (inputType != ED_d->startSystem.startSystemType)))
  {
    printf("ERROR: The type of start system is incorrect!\n\n");
    bexit(ERROR_CONFIGURATION);
  }

  // read in gamma
  if (old_MPType == 0)
  { // read in _d
    comp_d gamma;
    fscanf(IN, "%lf%lf", &gamma->r, &gamma->i); 

    // setup in ED_d
    if (MPType == 0 || MPType == 2)
    { 
      set_d(ED_d->startSystem.gamma, gamma);
    }

    // setup in ED_mp
    if (MPType == 1 || MPType == 2)
    { // set to _mp & _rat
      comp_d_to_mp_rat(ED_mp->startSystem.gamma, ED_mp->startSystem.gamma_rat, gamma, 15, 1024, 0, 0);
    }
  }
  else 
  { // read in _rat
    mpq_t gamma[2];
    init_rat(gamma);
    mpq_inp_str(gamma[0], IN, 10);
    mpq_canonicalize(gamma[0]);
    mpq_inp_str(gamma[1], IN, 10);
    mpq_canonicalize(gamma[1]);

    // setup in ED_d
    if (MPType == 0 || MPType == 2)
    {
      rat_to_d(ED_d->startSystem.gamma, gamma);
    }

    // setup in ED_mp
    if (MPType == 1 || MPType == 2)
    { // set to _mp & _rat
      mpq_set(ED_mp->startSystem.gamma_rat[0], gamma[0]);
      mpq_set(ED_mp->startSystem.gamma_rat[1], gamma[1]);
      rat_to_mp(ED_mp->startSystem.gamma, ED_mp->startSystem.gamma_rat);
    }

    clear_rat(gamma);
  }

  if (inputType == 1)
  { // setup coeff
    fscanf(IN, "%d%d", &totalDeg, &inputType);

    // verify input
    j = 0;
    if (MPType == 1)
    {
      for (i = 0; i < ED_mp->startSystem.size_r; i++)
        j += ED_mp->startSystem.degrees[i];
    }
    else
    {
      for (i = 0; i < ED_d->startSystem.size_r; i++)
        j += ED_d->startSystem.degrees[i];
    }

    if (j != totalDeg)
    {
      printf("ERROR: The total degree is incorrect!\n\n");
      bexit(ERROR_CONFIGURATION);
    }
    else if ((MPType == 1 && ED_mp->startSystem.coeff_cols != inputType) || ((MPType == 0 || MPType == 2) && ED_d->startSystem.coeff_cols != inputType))
    {
      printf("ERROR: The number of variables is incorrect!\n\n");
      bexit(ERROR_CONFIGURATION);
    }

    if (old_MPType == 0)
    { // read in _d
      comp_d **coeff_d = (comp_d **)bmalloc(totalDeg * sizeof(comp_d *));
      for (i = 0; i < totalDeg; i++)
      {
        coeff_d[i] = (comp_d *)bmalloc(inputType * sizeof(comp_d));
        for (j = 0; j < inputType; j++)
          fscanf(IN, "%lf%lf", &coeff_d[i][j]->r, &coeff_d[i][j]->i);
      }

      // setup in ED_d
      if (MPType == 0 || MPType == 2)
      {
        for (i = 0; i < totalDeg; i++)
          for (j = 0; j < inputType; j++)
          {
            set_d(ED_d->startSystem.coeff[i][j], coeff_d[i][j]);
          }
      }

      // setup in ED_mp
      if (MPType == 1 || MPType == 2)
      { // set to _mp & _rat
        for (i = 0; i < totalDeg; i++)
          for (j = 0; j < inputType; j++)
          {
            comp_d_to_mp_rat(ED_mp->startSystem.coeff[i][j], ED_mp->startSystem.coeff_rat[i][j], coeff_d[i][j], 15, 1024, 0, 0);
          }
      }

      // clear
      for (i = 0; i < totalDeg; i++)
        free(coeff_d[i]);
      free(coeff_d);
    }
    else
    { // read in _rat
      mpq_t ***coeff_rat = (mpq_t ***)bmalloc(totalDeg * sizeof(mpq_t **));
      for (i = 0; i < totalDeg; i++)
      {
        coeff_rat[i] = (mpq_t **)bmalloc(inputType * sizeof(mpq_t *));
        for (j = 0; j < inputType; j++)
        {
          coeff_rat[i][j] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
          init_rat(coeff_rat[i][j]);
          mpq_inp_str(coeff_rat[i][j][0], IN, 10);
          mpq_canonicalize(coeff_rat[i][j][0]);
          mpq_inp_str(coeff_rat[i][j][1], IN, 10);
          mpq_canonicalize(coeff_rat[i][j][1]);
        }
      }

      // setup in ED_d
      if (MPType == 0 || MPType == 2)
      {
        for (i = 0; i < totalDeg; i++)
          for (j = 0; j < inputType; j++)
          {
            rat_to_d(ED_d->startSystem.coeff[i][j], coeff_rat[i][j]);
          }
      }

      // setup in ED_mp
      if (MPType == 1 || MPType == 2)
      { // set to _mp & _rat
        for (i = 0; i < totalDeg; i++)
          for (j = 0; j < inputType; j++)
          {
            mpq_set(ED_mp->startSystem.coeff_rat[i][j][0], coeff_rat[i][j][0]);
            mpq_set(ED_mp->startSystem.coeff_rat[i][j][1], coeff_rat[i][j][1]);
            rat_to_mp(ED_mp->startSystem.coeff[i][j], ED_mp->startSystem.coeff_rat[i][j]);
          }
      }

      // clear
      for (i = 0; i < totalDeg; i++)
      {
        for (j = 0; j < inputType; j++)
        {
          clear_rat(coeff_rat[i][j]);
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
    }
  }

  return;
}

void readInMat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int outputMPType, int inputMPType, int expectedRows, int expectedCols, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads a matrix and saves it accordingly                *
\***************************************************************/
{
  int i, j, rows, cols;

  // read in size
  fscanf(IN, "%d%d", &rows, &cols);

  if (rows != expectedRows)
  {
    printf("ERROR: The matrix has the wrong number of rows!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (cols != expectedCols)
  {
    printf("ERROR: The matrix has the wrong number of columns!\n\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (inputMPType == 0)
  { // read in _d
    mat_d in;
    init_mat_d(in, rows, cols);
    in->rows = rows;
    in->cols = cols;

    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        fscanf(IN, "%lf%lf", &in->entry[i][j].r, &in->entry[i][j].i);
      }

    // setup
    if (outputMPType == 0 || outputMPType == 2)
    { // copy to _d
      mat_cp_d(A_d, in);
    }
    
    // setup
    if (outputMPType == 1 || outputMPType == 2)
    { // copy to _mp & _rat
      mat_d_to_mp(A_mp, in);
      mat_d_to_rat(A_rat, in);
    }

    // clear
    clear_mat_d(in);
  }
  else
  { // read in _rat
    mpq_t ***in;
    init_mat_rat(in, rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_inp_str(in[i][j][0], IN, 10);
        mpq_canonicalize(in[i][j][0]);
        mpq_inp_str(in[i][j][1], IN, 10);
        mpq_canonicalize(in[i][j][1]);
      }

    // setup
    if (outputMPType == 0 || outputMPType == 2)
    { // copy to _d
      mat_rat_to_d(A_d, in, rows, cols);
    }

    // setup
    if (outputMPType == 1 || outputMPType == 2)
    { // copy to _mp & _rat
      mat_rat_to_mp(A_mp, in, rows, cols);
      mat_cp_rat(A_rat, in, rows, cols);
    }

    // clear
    clear_mat_rat(in, rows, cols);
  }

  return;
}

void readInSquareSystem(FILE *IN, int MPType, int old_MPType, basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads the square system from IN and stores to BED      *
\***************************************************************/
{
  // setup B
  if (MPType == 0)
    readInMat(ED_d->squareSystem.B, NULL, NULL, MPType, old_MPType, ED_d->squareSystem.B->rows, ED_d->squareSystem.B->cols, IN);
  else if (MPType == 1)
    readInMat(NULL, ED_mp->squareSystem.B, ED_mp->squareSystem.B_rat, MPType, old_MPType, ED_mp->squareSystem.B->rows, ED_mp->squareSystem.B->cols, IN);
  else
    readInMat(ED_d->squareSystem.B, ED_mp->squareSystem.B, ED_mp->squareSystem.B_rat, MPType, old_MPType, ED_d->squareSystem.B->rows, ED_d->squareSystem.B->cols, IN);

  // setup B_perp
  if (MPType == 0)
    readInMat(ED_d->squareSystem.B_perp, NULL, NULL, MPType, old_MPType, ED_d->squareSystem.B_perp->rows, ED_d->squareSystem.B_perp->cols, IN);
  else if (MPType == 1)
    readInMat(NULL, ED_mp->squareSystem.B_perp, ED_mp->squareSystem.B_perp_rat, MPType, old_MPType, ED_mp->squareSystem.B_perp->rows, ED_mp->squareSystem.B_perp->cols, IN);
  else
    readInMat(ED_d->squareSystem.B_perp, ED_mp->squareSystem.B_perp, ED_mp->squareSystem.B_perp_rat, MPType, old_MPType, ED_d->squareSystem.B_perp->rows, ED_d->squareSystem.B_perp->cols, IN);

  // setup A
  if (MPType == 0)
    readInMat(ED_d->squareSystem.A, NULL, NULL, MPType, old_MPType, ED_d->squareSystem.A->rows, ED_d->squareSystem.A->cols, IN);
  else if (MPType == 1)
    readInMat(NULL, ED_mp->squareSystem.A, ED_mp->squareSystem.A_rat, MPType, old_MPType, ED_mp->squareSystem.A->rows, ED_mp->squareSystem.A->cols, IN);
  else
    readInMat(ED_d->squareSystem.A, ED_mp->squareSystem.A, ED_mp->squareSystem.A_rat, MPType, old_MPType, ED_d->squareSystem.A->rows, ED_d->squareSystem.A->cols, IN);

  return;
}

void updateZeroDimRelevantData(basic_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: update the relevant data to ED from FP  so that we can *
* recover the system if needed for sharpening                   *
\***************************************************************/
{
  int inputMPType, inputEqbyeq;

  fscanf(FP, "%d%d", &inputMPType, &inputEqbyeq);

  // verify same type of method
  if (inputEqbyeq && !eqbyeqMethod)
  {
    printf("ERROR: The system was originally solved using an equation-by-equation approach!\n\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!inputEqbyeq && eqbyeqMethod)
  {
    printf("ERROR: You are trying to use an equation-by-equation approach on a system that was originally solved using a standard homotopy!\n\n");
    bexit(ERROR_CONFIGURATION);
  }

  // read in the patch
  if (MPType == 0)
    readInPatch(FP, MPType, inputMPType, &ED_d->patch, NULL);
  else if (MPType == 1)
    readInPatch(FP, MPType, inputMPType, NULL, &ED_mp->patch);
  else
    readInPatch(FP, MPType, inputMPType, &ED_d->patch, &ED_mp->patch);

  // read in the start system
  readInStartSystem(FP, MPType, inputMPType, ED_d, ED_mp);

  // read in the square system
  readInSquareSystem(FP, MPType, inputMPType, ED_d, ED_mp);

  return;
}





