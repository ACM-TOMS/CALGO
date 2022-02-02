// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

// Given a general point on a generically reduced irreducible component 
// of dimension k and compute a witness set

int witnessSetGeneration(witness_t *W, tracker_config_t *T, int pathMod, char *witnessName, int genWitnessSet);
void initializeWitnessSet(int numPoints, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W);
int generateWitnessSet(int numPoints, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W);
int monodromyAndTraces(int *num_points, int *amp_trace_prec, int **traces_prec, comp_d **traces_d, comp_mp **traces_mp, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], vec_d trace_proj_d, vec_mp trace_proj_mp, mpq_t **trace_proj_rat, vec_d trace_slice_d, vec_mp trace_slice_mp, mpq_t **trace_slice_rat, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W, FILE *OUT, FILE *MIDOUT);

void witnessGeneration(unsigned int currentSeed, int MPType, char *startName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: main control function for generating witness set from  *
*   a given point                                               *
\***************************************************************/
{
  int rV, userHom = 0, useRegen = 0, regenStartLevel = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, maxCodim = 0, specificCodim = 0, paramHom = 0;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  tracker_config_t T;
  witness_t W;

  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType);

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

#ifdef _OPENMP
  #pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }

  // quit if we are using MPI
  if (num_processes > 1)
  { // exit since parallel witness generation is not implemented
#ifdef _HAVE_MPI
    printf("ERROR: Parallel witness set generation is not implemented. Please use sequential version!\n");
    bexit(ERROR_OTHER);
#endif
  }

  // initialize witnessSet
  W.curr_precision = T.Precision;
 
  // setup the slp
  W.Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  W.orig_variables = W.new_variables = setupProg(W.Prog, T.Precision, T.MPType);

  // setup W.PPD
  setupPreProcData("preproc_data", &W.PPD);

  // verify that we are using only 1 homogenous variable group
  if (W.PPD.num_hom_var_gp + W.PPD.num_var_gp > 1)
  { // exit immediately
    printf("ERROR: Positive dimensional setup is implemented for systems with only one variable group.\n");
    printf("  Please change the input file so that the variables are listed as a single variable group.\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find the rank
  if (T.MPType == 0 || T.MPType == 2)
    W.system_rank = rank_finder_d(&W.PPD, W.Prog, &T, W.orig_variables);
  else
    W.system_rank = rank_finder_mp(&W.PPD, W.Prog, &T, W.orig_variables);

  // setup the number of functions
  W.num_funcs = W.PPD.num_funcs;

  // setup orig_degrees, new_degrees & P
  setupDegrees_orig_new_perm(&W.orig_degrees, &W.new_degrees, &W.P, W.num_funcs, W.PPD.num_var_gp + W.PPD.num_hom_var_gp, "deg.out");

  // initialize codim
  W.num_codim = 0;
  W.codim = NULL;

  // setup gamma
  if (T.MPType == 0)
  { // only setup gamma_d
    get_comp_rand_d(W.gamma_d);
  }
  else if (T.MPType == 1)
  { // only setup gamma_mp
    init_mp(W.gamma_mp);
    get_comp_rand_mp(W.gamma_mp);
  }
  else
  { // setup gamma_rat, gamma_mp & gamma_d
    get_comp_rand_rat(W.gamma_d, W.gamma_mp, W.gamma_rat, W.curr_precision, T.AMP_max_prec, 1, 1);
  }

  // target slice info is not setup
  W.targetSliceInit = 0;

  // setup the number of variables
  T.numVars = W.Prog->numVars;

  // setup the rest of T, if needed
  if (T.MPType == 1)
  { // initialize latest_newton_residual_mp
    mpf_init(T.latest_newton_residual_mp);
  }
  else if (T.MPType == 2)
  { // initialize latest_newton_residual_mp
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);

    // setup eps, Phi & Psi
    T.AMP_eps = (double) T.numVars * T.numVars;
    T.AMP_Phi = T.AMP_bound_on_degree * (T.AMP_bound_on_degree - 1) * T.AMP_bound_on_abs_vals_of_coeffs;
    T.AMP_Psi = T.AMP_bound_on_degree * T.AMP_bound_on_abs_vals_of_coeffs;
  }

  // now that everything is initialized, do the actual witness set generation
  rV = witnessSetGeneration(&W, &T, pathMod, startName, constructWitnessSet);

  if (constructWitnessSet)
  { // check for success
    if (!rV)
    { // display decomposition chart
      numIrredDecompChart(&W, stdout, T.MPType, reducedOnly);

      // create output files
      numIrredDecompOutput(&W, &T, 6, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 6
    }
    else if (rV < 0)
    { // print error message
      printf("\nNOTE: Bertini was unable to compute a complete witness set.\n");
      printf("      Try increasing MaxNumMonoLoops.\n\n");
    }
  }

  // clear witnessSet
  witness_clear(&W, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void setupPerturbationPoint(point_d testPoint_d, point_mp testPoint_mp, point_d orig_d, point_mp orig_mp, int prec, double final_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perturb the given point                                *
\***************************************************************/
{
  int i;
  double tolDigits = -log10(final_tol), precDigits = prec_to_digits(prec) - 1;

  if (prec < 64)
  { // use _d
    double tempD, perturbationSize;
    comp_d tempComp;

    // setup the perturbationSize
    tempD = MAX(precDigits - tolDigits, 0.5);
    perturbationSize = tolDigits + (rand() / (RAND_MAX + 1.0)) * tempD; // random number 'between' tolDigits & precDigits
    perturbationSize = pow(10, -perturbationSize);

    // setup testPoint_d
    init_point_d(testPoint_d, orig_d->size);
    testPoint_d->size = orig_d->size;
    for (i = 0; i < orig_d->size; i++)
    {
      get_comp_rand_d(tempComp);
      mul_rdouble_d(tempComp, tempComp, perturbationSize);
      add_d(&testPoint_d->coord[i], &orig_d->coord[i], tempComp);
    }
  }
  else
  { // use _mp
    comp_mp tempComp;
    init_mp2(tempComp, prec);

    // setup testPoint_mp
    init_point_mp2(testPoint_mp, orig_mp->size, prec);
    testPoint_mp->size = orig_mp->size;
    for (i = 0; i < orig_mp->size; i++)
    {
      get_comp_rand_mp(tempComp);
      mul_rdouble_mp(tempComp, tempComp, final_tol);
      add_mp(&testPoint_mp->coord[i], &orig_mp->coord[i], tempComp);
    }
   
    clear_mp(tempComp);   
  }

  return;
}

void setupPoints(int *numPoints, point_d **p_d, point_mp **p_mp, int prec, int numVars, int num_var_gp, char *fileName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: read in a point from a file                            *
\***************************************************************/
{
  int i, j, rV, base = 10;
  FILE *IN = fopen(fileName, "r");

  // check that IN exists
  if (IN == NULL)
  { // file does not exist
    printf("\n\nERROR: '%s' does not exist!!!\n\n", fileName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // read in the number of points
  *numPoints = 0;
  fscanf(IN, "%d", numPoints);
  if (*numPoints < 1)
  { // error
    printf("\n\nERROR: The number of points in '%s' must be positive!!\n\n", fileName);
    bexit(ERROR_CONFIGURATION);
  }

  if (prec < 64)
  { // setup p_d
    *p_d = (point_d *)bmalloc(*numPoints * sizeof(point_d));
    for (j = 0; j < *numPoints; j++)
    { // intialize
      init_point_d((*p_d)[j], numVars);
      (*p_d)[j]->size = numVars;

      if (num_var_gp)
      { // Bertini homogenized the variable group
        set_one_d(&(*p_d)[j]->coord[0]); // set hom coord == 1
        for (i = 1; i < numVars; i++)
        { // read in coordinates
          rV = fscanf(IN, "%lf%lf", &(*p_d)[j]->coord[i].r, &(*p_d)[j]->coord[i].i);
          if (rV < 0)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", fileName);
            bexit(ERROR_INVALID_SIZE);
          }
          // scan in rest of line
          rV = scanRestOfLine(IN);
          if (rV && i + 1 < numVars)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", fileName);
            bexit(ERROR_INVALID_SIZE);
          }
        }
      }
      else
      { // already homogenized variable group
        for (i = 0; i < numVars; i++)
        { // read in coordinates
          rV = fscanf(IN, "%lf%lf", &(*p_d)[j]->coord[i].r, &(*p_d)[j]->coord[i].i);
          if (rV < 0)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", fileName);
            bexit(ERROR_INVALID_SIZE);
          }
          // scan in rest of line
          rV = scanRestOfLine(IN);
          if (rV && i + 1 < numVars)
          { // we are at EOF
            printf("\n\nERROR: There are not enough coordinates in %s!!!\n\n", fileName);
            bexit(ERROR_INVALID_SIZE);
          }
        }
      }
    }
  }
  else
  { // setup p_mp
    *p_mp = (point_mp *)bmalloc(*numPoints * sizeof(point_mp));
    for (j = 0; j < *numPoints; j++)
    { // intialize
      init_point_mp2((*p_mp)[j], numVars, prec);
      (*p_mp)[j]->size = numVars;

      if (num_var_gp)
      { // Bertini homogenized the variable group
        set_one_mp(&(*p_mp)[j]->coord[0]); // set hom coord == 1
        for (i = 1; i < numVars; i++)
        {
          mpf_inp_str((*p_mp)[j]->coord[i].r, IN, base);
          mpf_inp_str((*p_mp)[j]->coord[i].i, IN, base);
          // scan rest of line
          scanRestOfLine(IN);
        }
      }
      else
      { // already homogenized variable group
        for (i = 0; i < numVars; i++)
        {
          mpf_inp_str((*p_mp)[j]->coord[i].r, IN, base);
          mpf_inp_str((*p_mp)[j]->coord[i].i, IN, base);
          // scan rest of line
          scanRestOfLine(IN);
        }
      }
    }
  }

  // close file
  fclose(IN);

  // test for distinct points


  return;
}

int computeJacobianNull(point_d p_d, point_d test_d, point_mp p_mp, point_mp test_mp, int prec, tracker_config_t *T, prog_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: dimension of null space of Jacobian            *
* NOTES: compute nullity of the Jacobian using 2 approximations *
\***************************************************************/
{ 
  int nullity = 0;
  double CN, s, l, max_CN, max_SV_ratio, SV_tol;

  if (prec < 64)
  { // use _d
    comp_d time;
    mat_d J;
    eval_struct_d e;

    // initialize
    max_CN = 1e13;
    max_SV_ratio = T->ratioTol; 
    SV_tol = MAX(T->sing_val_zero_tol, 1e-15);
    set_zero_d(time);
    init_mat_d(J, 0, 0);
    init_eval_struct_d(e, 0, 0, 0);

    // evaluate at p_d
    evalProg_d(e.funcVals, e.parVals, e.parDer, J, e.Jp, p_d, time, Prog);

    // verify that the point sufficiently satisfies the polynomial system
    CN = infNormVec_d(e.funcVals);
    if (CN > T->final_tol_times_mult)
    {
      printf("ERROR: The point does not sufficiently satisfy the original system (residual: %e, tolerance: %e).\n", CN, T->final_tol_times_mult);
      bexit(ERROR_INPUT_SYSTEM);
    }

    // evaluate at test_d
    evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, test_d, time, Prog);

    // compute nullity
    nullity = corank_rrv_d(&CN, &s, &l, J, e.Jv, 0, 0, max_CN, max_SV_ratio, SV_tol);

    // clear memory
    clear_mat_d(J);
    clear_eval_struct_d(e);
  }
  else 
  { // use _mp
    comp_mp time;
    mat_mp J;
    eval_struct_mp e;

    // set the precision correctly
    initMP(prec);
    Prog->precision = prec;

    // initialize
    nullity = prec_to_digits(T->Precision) - 4;
    max_CN = MIN(1e150, pow(10, nullity));
    max_SV_ratio = T->ratioTol; 
    SV_tol = MAX(T->sing_val_zero_tol, pow(10, -nullity - 2));
    init_mp(time); 
    set_zero_mp(time);
    init_mat_mp(J, 0, 0);
    init_eval_struct_mp(e, 0, 0, 0);

    // evaluate at p_mp
    evalProg_mp(e.funcVals, e.parVals, e.parDer, J, e.Jp, p_mp, time, Prog);

    // verify that the point sufficiently satisfies the polynomial system
    CN = infNormVec_mp(e.funcVals);
    if (CN > T->final_tol_times_mult)
    {
      printf("ERROR: The point does not sufficiently satisfy the original system (residual: %e, tolerance: %e).\n", CN, T->final_tol_times_mult);
      bexit(ERROR_INPUT_SYSTEM);
    }
  
    // evaluate at test_mp
    evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, test_mp, time, Prog);

    // compute corank
    nullity = corank_rrv_mp(&CN, &s, &l, J, e.Jv, 0, 0, max_CN, max_SV_ratio, SV_tol);

    // clear memory
    clear_mp(time);
    clear_mat_mp(J);
    clear_eval_struct_mp(e);
  }

  return nullity;
}

int witnessSetGeneration(witness_t *W, tracker_config_t *T, int pathMod, char *witnessName, int genWitnessSet)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - success, -1 - error, 1 - successful test   *
* NOTES: generate a witness set from a given point              *
\***************************************************************/
{
  int i, rV, nullity = 0, numVars = T->numVars, prec = T->Precision, *newPoints_prec = NULL, numPoints = 0, dimension = 0;
  point_d *witnessPoints_d = NULL, testPoint_d, *newPoints_d = NULL;
  point_mp *witnessPoints_mp = NULL, testPoint_mp, *newPoints_mp = NULL;
  membership_slice_moving_t sliceMover;
  FILE *F = NULL;

  // setup the precision to use
  if (T->MPType == 0)
  { // use D
    prec = 52;
  }
  else if (T->MPType == 1)
  { // use MP
    prec = T->Precision;
  }
  else
  { // use enough precision based on final tolerance
    prec = -floor(log10(T->final_tolerance) - 0.5);
    prec = digits_to_prec(prec);
  }

  // setup witnessPoints - verify numPoints > 0
  setupPoints(&numPoints, &witnessPoints_d, &witnessPoints_mp, prec, numVars, W->PPD.num_var_gp, witnessName);

  for (i = 0; i < numPoints; i++)
  { // setup perturbation of witnessPoint
    if (prec < 64)
      setupPerturbationPoint(testPoint_d, testPoint_mp, witnessPoints_d[i], NULL, prec, T->final_tolerance);
    else
      setupPerturbationPoint(testPoint_d, testPoint_mp, NULL, witnessPoints_mp[i], prec, T->final_tolerance);

    // compute the dimension of the null space of the Jacobian at this point
    if (prec < 64)
      rV = computeJacobianNull(witnessPoints_d[i], testPoint_d, NULL, testPoint_mp, prec, T, W->Prog);
    else
      rV = computeJacobianNull(NULL, testPoint_d, witnessPoints_mp[i], testPoint_mp, prec, T, W->Prog);

    if (i == 0)
      nullity = rV;
    else if (nullity != rV)
    {
      printf("ERROR: All points must have the same null space dimension!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
  }

  // setup dimension
  dimension = nullity - W->PPD.num_var_gp  - W->PPD.num_hom_var_gp;

  // error checking on the nullity - verify that nullity >= 0 
  if (dimension < 0)
  {
    printf("ERROR: The Jacobian must have a nonnegative dimensional null space!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // allocate memory for newPoints
  newPoints_prec = (int *)bmalloc(numPoints * sizeof(int));
  newPoints_d = (point_d *)bmalloc(numPoints * sizeof(point_d));
  newPoints_mp = (point_mp *)bmalloc(numPoints * sizeof(point_mp));
  for (i = 0; i < numPoints; i++)
    newPoints_prec[i] = 0;

  // test that there does exist a generically reduced component (isosingular test)
  printf("Testing for a component of dimension %d.\n\n", dimension); 
  rV = isosingularDimTest(&sliceMover, newPoints_d, newPoints_mp, newPoints_prec, nullity, numPoints, witnessPoints_d, witnessPoints_mp, prec, T, W);

  // create an "isosingular_summary" file
  F = fopen("isosingular_summary", "w");
  fprintf(F, "%d %d\n\n", dimension, !rV);
  fclose(F);

  if (rV)
  { // print error message
    printf("Bertini was unable to verify that a witness point is a smooth\npoint on a %d dimensional generically reduced component.\n\n", dimension);
    rV = 1;
  }
  else if (genWitnessSet)
  { // setup the basic witness set
    initializeWitnessSet(numPoints, newPoints_d, newPoints_mp, newPoints_prec, &sliceMover, T, W);

    // use monodromy and trace test to complete the witness set
    rV = -generateWitnessSet(numPoints, newPoints_d, newPoints_mp, newPoints_prec, &sliceMover, T, W);
  }
  else
  { // print success message
    printf("Bertini verified that a witness point is a smooth\npoint on a %d dimensional generically reduced component.\n\n", dimension); 
    rV = 1;
  }

  // clear witnessPoint & testPoint
  if (prec < 64)
  { // clear _d
    for (i = 0; i < numPoints; i++) 
      clear_point_d(witnessPoints_d[i]);
    clear_point_d(testPoint_d);
  }
  else
  { // clear _mp
    for (i = 0; i < numPoints; i++)
      clear_point_mp(witnessPoints_mp[i]);
    clear_point_mp(testPoint_mp);
  }
  free(witnessPoints_d);
  free(witnessPoints_mp);

  // clear newPoints
  for (i = 0; i < numPoints; i++)
    if (0 < newPoints_prec[i] && newPoints_prec[i] < 64)
    { // clear _d
      clear_point_d(newPoints_d[i]);
    }
    else if (newPoints_prec[i] >= 64)
    { // clear _mp
      clear_point_mp(newPoints_mp[i]);
    }
  free(newPoints_d);
  free(newPoints_mp);
  free(newPoints_prec);

  // clear sliceMover
  clear_slice_mover(&sliceMover, T->MPType); 

  return rV;
}

void initializeWitnessSet(int numPoints, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: basic setup for the witness set                        *
\***************************************************************/
{
  int i, j, codim = sliceMover->curr_codim;
  
  // allocate a codim in W
  W->num_codim = 1;
  W->curr_codim_index = 0;
  W->codim = (witnessCodim_t *)bmalloc(1 * sizeof(witnessCodim_t)); 

  // setup the codim
  W->codim[0].codim = codim;
  W->codim[0].num_set = 1;
  W->codim[0].num_nonsing = 1;
  W->codim[0].num_sing = 0;

  // setup A using sliceMover
  if (T->MPType == 0)
  { // setup A_d
    W->codim[0].A_rows = sliceMover->A_d->rows;
    W->codim[0].A_cols = sliceMover->A_d->cols;
    init_mat_d(W->codim[0].A_d, W->codim[0].A_rows, W->codim[0].A_cols);
    mat_cp_d(W->codim[0].A_d, sliceMover->A_d);

    // setup A_rat
    init_mat_rat(W->codim[0].A_rat, W->codim[0].A_rows, W->codim[0].A_cols);
    for (i = 0; i < W->codim[0].A_rows; i++)
      for (j = 0; j < W->codim[0].A_cols; j++)
      {
        mpq_set_d(W->codim[0].A_rat[i][j][0], W->codim[0].A_d->entry[i][j].r);
        mpq_set_d(W->codim[0].A_rat[i][j][1], W->codim[0].A_d->entry[i][j].i);
      } 
  }
  else if (T->MPType == 1)
  { // setup A_mp
    W->codim[0].A_rows = sliceMover->A_mp->rows;
    W->codim[0].A_cols = sliceMover->A_mp->cols;
    init_mat_mp(W->codim[0].A_mp, W->codim[0].A_rows, W->codim[0].A_cols);
    mat_cp_mp(W->codim[0].A_mp, sliceMover->A_mp);

    // setup A_rat
    init_mat_rat(W->codim[0].A_rat, W->codim[0].A_rows, W->codim[0].A_cols);
    for (i = 0; i < W->codim[0].A_rows; i++)
      for (j = 0; j < W->codim[0].A_cols; j++)
      {
        mpf_t_to_rat(W->codim[0].A_rat[i][j][0], W->codim[0].A_mp->entry[i][j].r);
        mpf_t_to_rat(W->codim[0].A_rat[i][j][1], W->codim[0].A_mp->entry[i][j].i);
      }
  }
  else
  { // setup A_d, A_mp, A_rat
    W->codim[0].A_rows = sliceMover->A_d->rows;
    W->codim[0].A_cols = sliceMover->A_d->cols;
    init_mat_d(W->codim[0].A_d, W->codim[0].A_rows, W->codim[0].A_cols);
    init_mat_mp2(W->codim[0].A_mp, W->codim[0].A_rows, W->codim[0].A_cols, W->curr_precision);
    init_mat_rat(W->codim[0].A_rat, W->codim[0].A_rows, W->codim[0].A_cols);

    W->codim[0].A_rows = sliceMover->A_d->rows;
    W->codim[0].A_cols = sliceMover->A_d->cols;
    mat_cp_d(W->codim[0].A_d, sliceMover->A_d);
    mat_cp_mp(W->codim[0].A_mp, sliceMover->A_mp);
    for (i = 0; i < W->codim[0].A_rows; i++)
      for (j = 0; j < W->codim[0].A_cols; j++)
      {
        mpf_t_to_rat(W->codim[0].A_rat[i][j][0], W->codim[0].A_mp->entry[i][j].r);
        mpf_t_to_rat(W->codim[0].A_rat[i][j][1], W->codim[0].A_mp->entry[i][j].i);
      }
  }

  // setup W
  W->codim[0].W = (int **)bmalloc(codim * sizeof(int *));
  for (i = 0; i < codim; i++)
  {
    W->codim[0].W[i] = (int *)bmalloc((sliceMover->Prog->numFuncs - codim) * sizeof(int));
    for (j = 0; j < sliceMover->Prog->numFuncs - codim; j++)
      W->codim[0].W[i][j] = W->new_degrees[i] - W->new_degrees[j + codim];
  }

  // setup H & homVarConst
  if (T->MPType == 0)
  { // setup H_d & homVarConst_d
    init_vec_d(W->codim[0].H_d, W->orig_variables);
    W->codim[0].H_d->size = W->orig_variables;
    if (W->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized: H_d = [1,0..0]
      set_one_d(&W->codim[0].H_d->coord[0]);
      for (i = 1; i < W->orig_variables; i++)
      {
        set_zero_d(&W->codim[0].H_d->coord[i]);
      }
      // setup homVarConst_d to be 0
      set_zero_d(W->codim[0].homVarConst_d);
    }
    else
    { // using a homogeneous variable group

      // setup H_d to be random
      make_vec_random_d(W->codim[0].H_d, W->orig_variables);

      // setup homVarConst_d to be random
      get_comp_rand_d(W->codim[0].homVarConst_d);
    }
  }
  else if (T->MPType == 1)
  { // setup H_mp & homVarConst_mp
    init_vec_mp(W->codim[0].H_mp, W->orig_variables);
    W->codim[0].H_mp->size = W->orig_variables;
    init_mp(W->codim[0].homVarConst_mp);
    if (W->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized: H_mp = [1,0..0]
      set_one_mp(&W->codim[0].H_mp->coord[0]);
      for (i = 1; i < W->orig_variables; i++)
      {
        set_zero_mp(&W->codim[0].H_mp->coord[i]);
      }
      // setup homVarConst_mp to be 0
      set_zero_mp(W->codim[0].homVarConst_mp);
    }
    else
    { // using a homogeneous variable group

      // setup H_mp to be random
      make_vec_random_mp(W->codim[0].H_mp, W->orig_variables);

      // setup homVarConst_mp to be random
      get_comp_rand_mp(W->codim[0].homVarConst_mp);
    }
  }
  else
  { // setup _d,_mp,_rat
    init_vec_d(W->codim[0].H_d, W->orig_variables);
    init_vec_mp(W->codim[0].H_mp, W->orig_variables);
    init_vec_rat(W->codim[0].H_rat, W->orig_variables);
    init_mp(W->codim[0].homVarConst_mp);
    init_rat(W->codim[0].homVarConst_rat);
    W->codim[0].H_d->size = W->orig_variables;
    W->codim[0].H_mp->size = W->orig_variables;
    if (W->PPD.num_var_gp > 0)
    { // using a variable group that was originally not homogenized: H = [1,0..0]
      set_one_d(&W->codim[0].H_d->coord[0]);
      set_one_mp(&W->codim[0].H_mp->coord[0]);
      set_one_rat(W->codim[0].H_rat[0]);
      for (i = 1; i < W->orig_variables; i++)
      {
        set_zero_d(&W->codim[0].H_d->coord[i]);
        set_zero_mp(&W->codim[0].H_mp->coord[i]);
        set_zero_rat(W->codim[0].H_rat[i]);
      }
      // setup homVarConst to be 0
      set_zero_d(W->codim[0].homVarConst_d);
      set_zero_mp(W->codim[0].homVarConst_mp);
      set_zero_rat(W->codim[0].homVarConst_rat);
    }
    else
    { // using a homogeneous variable group

      // setup H to be random
      make_vec_random_rat(W->codim[0].H_d, W->codim[0].H_mp, W->codim[0].H_rat, W->orig_variables, W->curr_precision, T->AMP_max_prec, 0, 0);

      // setup homVarConst to be random
      get_comp_rand_rat(W->codim[0].homVarConst_d, W->codim[0].homVarConst_mp, W->codim[0].homVarConst_rat, W->curr_precision, T->AMP_max_prec, 0, 0);
    }
  }

  // setup B & p
  if (T->MPType == 0)
  { // setup B_d & p_d
    init_mat_d(W->codim[0].B_d, sliceMover->B_d->rows, sliceMover->B_d->cols);
    mat_cp_d(W->codim[0].B_d, sliceMover->B_d);

    init_vec_d(W->codim[0].p_d, sliceMover->p_d->size);
    vec_cp_d(W->codim[0].p_d, sliceMover->p_d);
  }
  else if (T->MPType == 1)
  { // setup B_mp & p_mp
    init_mat_mp(W->codim[0].B_mp, sliceMover->B_mp->rows, sliceMover->B_mp->cols);
    mat_cp_mp(W->codim[0].B_mp, sliceMover->B_mp);

    init_vec_mp(W->codim[0].p_mp, sliceMover->p_mp->size);
    vec_cp_mp(W->codim[0].p_mp, sliceMover->p_mp);
  }
  else
  { // setup B & p
    init_mat_d(W->codim[0].B_d, sliceMover->B_d->rows, sliceMover->B_d->cols);
    init_mat_mp2(W->codim[0].B_mp, sliceMover->B_d->rows, sliceMover->B_d->cols, sliceMover->curr_precision);
    init_mat_rat(W->codim[0].B_rat, sliceMover->B_d->rows, sliceMover->B_d->cols);

    mat_cp_d(W->codim[0].B_d, sliceMover->B_d);
    mat_cp_mp(W->codim[0].B_mp, sliceMover->B_mp);
    mat_cp_rat(W->codim[0].B_rat, sliceMover->B_rat, sliceMover->B_d->rows, sliceMover->B_mp->cols);

    init_vec_d(W->codim[0].p_d, sliceMover->p_d->size);
    init_vec_mp2(W->codim[0].p_mp, sliceMover->p_d->size, sliceMover->curr_precision);
    init_vec_rat(W->codim[0].p_rat, sliceMover->p_d->size);    

    vec_cp_d(W->codim[0].p_d, sliceMover->p_d);
    vec_cp_mp(W->codim[0].p_mp, sliceMover->p_mp);
    vec_cp_rat(W->codim[0].p_rat, sliceMover->p_rat, sliceMover->p_d->size);
  }

  // setup 1 component
  W->codim[0].num_components = 1;
  W->codim[0].num_set = numPoints;
  W->codim[0].num_nonsing = numPoints;
  W->codim[0].num_sing = 0;

  // setup numPoints points
  W->codim[0].witnessPt_types = (int *)bmalloc(numPoints * sizeof(int));
  W->codim[0].component_nums = (int *)bmalloc(numPoints * sizeof(int));
  W->codim[0].multiplicities = (int *)bmalloc(numPoints * sizeof(int));
  W->codim[0].deflations_needed = (int *)bmalloc(numPoints * sizeof(int));
  for (i = 0; i < numPoints; i++)
  {
    W->codim[0].witnessPt_types[i] = NON_SINGULAR;
    W->codim[0].component_nums[i] = 0;
    W->codim[0].multiplicities[i] = 1;
    W->codim[0].deflations_needed[i] = 0;
  }

  // setup points
  if (T->MPType == 0)
  { // setup witnessPts_d
    W->codim[0].witnessPts_d = (endpoint_data_d *)bmalloc(numPoints * sizeof(endpoint_data_d));
    for (i = 0; i < numPoints; i++)
    {
      init_endpoint_data_d(&W->codim[0].witnessPts_d[i]);
      vec_cp_d(W->codim[0].witnessPts_d[i].endPt, currPoints_d[i]);
      vec_cp_d(W->codim[0].witnessPts_d[i].last_approx, currPoints_d[i]);
      set_zero_d(W->codim[0].witnessPts_d[i].finalT);
      W->codim[0].witnessPts_d[i].cond_num = 1;
      W->codim[0].witnessPts_d[i].corank = 0;
      W->codim[0].witnessPts_d[i].smallest_nonzero_SV = 1;
      W->codim[0].witnessPts_d[i].largest_zero_SV = 0;
      W->codim[0].witnessPts_d[i].retVal = 0;
    }
  }
  else if (T->MPType == 1)
  { // setup witnessPts_mp
    W->codim[0].witnessPts_mp = (endpoint_data_mp *)bmalloc(numPoints * sizeof(endpoint_data_mp));
    for (i = 0; i < numPoints; i++)
    {
      init_endpoint_data_mp(&W->codim[0].witnessPts_mp[i]);
      vec_cp_mp(W->codim[0].witnessPts_mp[i].endPt, currPoints_mp[i]);
      vec_cp_mp(W->codim[0].witnessPts_mp[i].last_approx, currPoints_mp[i]);
      set_zero_mp(W->codim[0].witnessPts_mp[i].finalT);
      W->codim[0].witnessPts_mp[i].cond_num = 1;
      W->codim[0].witnessPts_mp[i].corank = 0;
      W->codim[0].witnessPts_mp[i].smallest_nonzero_SV = 1;
      W->codim[0].witnessPts_mp[i].largest_zero_SV = 0;
      W->codim[0].witnessPts_mp[i].retVal = 0;
    }
  }
  else 
  { // setup witnessPts_amp
    W->codim[0].witnessPts_amp = (endpoint_data_amp *)bmalloc(numPoints * sizeof(endpoint_data_amp));
    for (i = 0; i < numPoints; i++)
    {
      init_endpoint_data_amp(&W->codim[0].witnessPts_amp[i], currPoints_prec[i], currPoints_prec[i]);
      W->codim[0].witnessPts_amp[i].curr_prec = W->codim[0].witnessPts_amp[i].last_approx_prec = currPoints_prec[i];
      if (currPoints_prec[i] < 64)
      {
        vec_cp_d(W->codim[0].witnessPts_amp[i].endPt_d, currPoints_d[i]);
        set_zero_d(W->codim[0].witnessPts_amp[i].finalT_d);
        vec_cp_d(W->codim[0].witnessPts_amp[i].last_approx_d, currPoints_d[i]);
      }
      else
      {
        vec_cp_mp(W->codim[0].witnessPts_amp[i].endPt_mp, currPoints_mp[i]);
        set_zero_mp(W->codim[0].witnessPts_amp[i].finalT_mp);
        vec_cp_mp(W->codim[0].witnessPts_amp[i].last_approx_mp, currPoints_mp[i]);
      }
      W->codim[0].witnessPts_amp[i].cond_num = 1;
      W->codim[0].witnessPts_amp[i].corank = 0;
      W->codim[0].witnessPts_amp[i].smallest_nonzero_SV = 1;
      W->codim[0].witnessPts_amp[i].largest_zero_SV = 0;
      W->codim[0].witnessPts_amp[i].retVal = 0;
    }
  }

  return;
}

int generateWitnessSet(int numPoints, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - success, 1 - failure                       *
* NOTES: try to generate the witness set                        *
\***************************************************************/
{
  int i, amp_trace_prec = 52, trace_num_computed = 0, isComplete = 0;
  int *traces_prec = T->MPType == 2 ? (int *)bmalloc(numPoints * sizeof(int)) : NULL;
  comp_d s_d[2];
  comp_mp s_mp[2];
  mpq_t s_rat[2][2];
  vec_d trace_proj_d, trace_slice_d;
  vec_mp trace_proj_mp, trace_slice_mp;
  mpq_t **trace_proj_rat = NULL, **trace_slice_rat = NULL;
  comp_d sum_traces_d, *traces_d = (T->MPType == 0 || T->MPType == 2) ? (comp_d *)bmalloc(numPoints * sizeof(comp_d)) : NULL;
  comp_mp sum_traces_mp, *traces_mp = (T->MPType == 1 || T->MPType == 2) ? (comp_mp *)bmalloc(numPoints * sizeof(comp_mp)) : NULL;
  endpoint_data_d endPt_d;
  endpoint_data_mp endPt_mp;
  endpoint_data_amp endPt_amp;
  FILE *OUT = fopen("output_witness", "w"), *MIDOUT = fopen("midout", "w");

  // setup projection vector, slice vector, s, and gamma
  if (T->MPType == 0)
  { // setup trace_proj_d
    init_vec_d(trace_proj_d, sliceMover->orig_variables);
    make_vec_random_d(trace_proj_d, sliceMover->orig_variables);
    // setup trace_slice_d
    init_vec_d(trace_slice_d, sliceMover->B_d->rows);
    make_vec_random_d(trace_slice_d, sliceMover->B_d->rows);
    // setup s_d
    get_comp_rand_d(s_d[0]);
    get_comp_rand_d(s_d[1]);
  }
  else if (T->MPType == 1)
  { // setup trace_proj_mp
    init_vec_mp(trace_proj_mp, sliceMover->orig_variables);
    make_vec_random_mp(trace_proj_mp, sliceMover->orig_variables);
    // setup trace_slice_mp
    init_vec_mp(trace_slice_mp, sliceMover->B_mp->rows);
    make_vec_random_mp(trace_slice_mp, sliceMover->B_mp->rows);
    // setup s_mp
    init_mp(s_mp[0]);
    init_mp(s_mp[1]);
    get_comp_rand_mp(s_mp[0]);
    get_comp_rand_mp(s_mp[1]);
  }
  else
  { // setup trace_proj
    init_vec_d(trace_proj_d, sliceMover->orig_variables);
    init_vec_mp(trace_proj_mp, sliceMover->orig_variables);
    init_vec_rat(trace_proj_rat, sliceMover->orig_variables);
    make_vec_random_rat(trace_proj_d, trace_proj_mp, trace_proj_rat, sliceMover->orig_variables, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);
    // setup trace_slice
    init_vec_d(trace_slice_d, sliceMover->B_d->rows);
    init_vec_mp(trace_slice_mp, sliceMover->B_d->rows);
    init_vec_rat(trace_slice_rat, sliceMover->B_d->rows);
    make_vec_random_rat(trace_slice_d, trace_slice_mp, trace_slice_rat, sliceMover->B_d->rows, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);
    // setup s_d
    get_comp_rand_rat(s_d[0], s_mp[0], s_rat[0], sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
    get_comp_rand_rat(s_d[1], s_mp[1], s_rat[1], sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // initialize memory
  init_mp(sum_traces_mp);
  if (T->MPType == 0)
  {
    init_endpoint_data_d(&endPt_d);
  }
  else if (T->MPType == 1)
  {
    init_endpoint_data_mp(&endPt_mp);
  }
  else
  {
    init_endpoint_data_amp(&endPt_amp, 64, 64);
  }

  // compute the traces
  printf("Calculating trace%s for %d point%s\n", numPoints == 1 ? "" : "s", numPoints, numPoints == 1 ? "" : "s");
  for (i = 0; i < numPoints; i++)
  {
    if (T->MPType == 0)
    { // compute traces_d[0]
      point_cp_d(endPt_d.endPt, currPoints_d[i]);
      set_zero_d(traces_d[i]);
      // calculate using double precision
      calculateTrace(traces_d[i], NULL, NULL, sliceMover, sliceMover->Prog, &endPt_d, NULL, NULL, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);
    }
    else if (T->MPType == 1)
    { // compute traces_mp[0]
      point_cp_mp(endPt_mp.endPt, currPoints_mp[i]);
      init_mp(traces_mp[i]);
      set_zero_mp(traces_mp[i]);
      // calculate using multi precision
      calculateTrace(NULL, traces_mp[i], NULL, sliceMover, sliceMover->Prog, NULL, &endPt_mp, NULL, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);
    }
    else
    { // compute traces_d/_mp
      endPt_amp.curr_prec = currPoints_prec[i];
      if (currPoints_prec[i] < 64)
      {
        point_cp_d(endPt_amp.endPt_d, currPoints_d[i]);
      }
      else
      {
        setprec_point_mp(endPt_amp.endPt_mp, endPt_amp.curr_prec);
        point_cp_mp(endPt_amp.endPt_mp, currPoints_mp[i]);
      }

      traces_prec[i] = 52;
      set_zero_d(traces_d[i]);
      init_mp(traces_mp[i]);
      // calculate using AMP
      calculateTrace(traces_d[i], traces_mp[i], &traces_prec[i], sliceMover, sliceMover->Prog, NULL, NULL, &endPt_amp, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);

      // update amp_trace_prec, if needed
      if (amp_trace_prec < traces_prec[i])
        amp_trace_prec = traces_prec[i];
    }
    // increment the number of traces computed
    trace_num_computed++;
  }

  // setup trace_sum and determine if zero
  if (T->MPType == 0)
  {
    set_zero_d(sum_traces_d);
    for (i = 0; i < numPoints; i++)
      add_d(sum_traces_d, sum_traces_d, traces_d[i]);

    printf("Trace residual for %d point%s: %e\n\n", numPoints, numPoints == 1 ? "" : "s", d_abs_d(sum_traces_d));
    if (d_abs_d(sum_traces_d) < T->final_tol_times_mult)
      isComplete = 1;
  }
  else if (T->MPType == 1)
  {
    set_zero_mp(sum_traces_mp);
    for (i = 0; i < numPoints; i++)
      add_mp(sum_traces_mp, sum_traces_mp, traces_mp[i]);

    printf("Trace residual for %d point%s: %e\n\n", numPoints, numPoints == 1 ? "" : "s", d_abs_mp(sum_traces_mp));
    if (d_abs_mp(sum_traces_mp) < T->final_tol_times_mult)
      isComplete = 1;
  }
  else
  {
    if (amp_trace_prec < 64)
    {
      set_zero_d(sum_traces_d);
      for (i = 0; i < numPoints; i++)
        add_d(sum_traces_d, sum_traces_d, traces_d[i]);

      printf("Trace residual for %d point%s: %e\n\n", numPoints, numPoints == 1 ? "" : "s", d_abs_d(sum_traces_d));
      if (d_abs_d(sum_traces_d) < T->final_tol_times_mult)
        isComplete = 1;
    }
    else
    {
      setprec_mp(sum_traces_mp, amp_trace_prec);
      set_zero_mp(sum_traces_mp);
      for (i = 0; i < numPoints; i++)
      {
        if (traces_prec[i] < 64)
          d_to_mp(traces_mp[i], traces_d[i]);

        add_mp(sum_traces_mp, sum_traces_mp, traces_mp[i]);
      }

      printf("Trace residual for %d point%s: %e\n\n", numPoints, numPoints == 1 ? "" : "s", d_abs_mp(sum_traces_mp));
      if (d_abs_mp(sum_traces_mp) < T->final_tol_times_mult)
        isComplete = 1;
    }
  }

  if (!isComplete)
  { // run monodromy with the trace test
    isComplete = !monodromyAndTraces(&trace_num_computed, &amp_trace_prec, &traces_prec, &traces_d, &traces_mp, s_d, s_mp, s_rat, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, currPoints_d, currPoints_mp, currPoints_prec, sliceMover, T, W, OUT, MIDOUT);
  }

  // close files and delete midout
  fclose(OUT);
  fclose(MIDOUT);
  remove("midout"); 

  // clear memory
  clear_vec(trace_proj_d, trace_proj_mp, trace_proj_rat, T->MPType);
  clear_vec(trace_slice_d, trace_slice_mp, trace_slice_rat, T->MPType);
  clear_d_mp_rat(s_d[0], s_mp[0], s_rat[0], T->MPType);
  clear_d_mp_rat(s_d[1], s_mp[1], s_rat[1], T->MPType);
  clear_mp(sum_traces_mp);

  if (T->MPType == 0)
  { // clear traces_d & endPt_d
    free(traces_d);
    clear_endpoint_data_d(&endPt_d);
  }
  else if (T->MPType == 1)
  { // clear traces_mp & endPt_mp
    for (i = 0; i < trace_num_computed; i++)
    {
      clear_mp(traces_mp[i]);
    }
    free(traces_mp);
    clear_endpoint_data_mp(&endPt_mp);
  }
  else
  { // clear
    for (i = 0; i < trace_num_computed; i++)
    {
      clear_mp(traces_mp[i]);
    }
    free(traces_d);
    free(traces_mp);
    free(traces_prec);
    clear_endpoint_data_amp(&endPt_amp);
  }

  return !isComplete;
}

int monodromyAndTraces(int *num_points, int *amp_trace_prec, int **traces_prec, comp_d **traces_d, comp_mp **traces_mp, comp_d s_d[2], comp_mp s_mp[2], mpq_t s_rat[2][2], vec_d trace_proj_d, vec_mp trace_proj_mp, mpq_t **trace_proj_rat, vec_d trace_slice_d, vec_mp trace_slice_mp, mpq_t **trace_slice_rat, point_d *currPoints_d, point_mp *currPoints_mp, int *currPoints_prec, membership_slice_moving_t *sliceMover, tracker_config_t *T, witness_t *W, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - success, 1 - failure                       *
* NOTES: try to generate the witness set                        *
\***************************************************************/
{
  int i, j, rV, matchingPathNum, size, its = 0, isComplete = 0, maxIts = T->max_num_mon_linears;
  double tol = T->final_tol_times_mult;
  endgame_data_t monodromyPt;
  mpf_t normDiff;
  FILE *WITPTS = fopen("witness_points", "w");
  fpos_t witPositionTop, witPositionBottom;

  init_endgame_data(&monodromyPt, T->Precision);
  mpf_init(normDiff);

  if (T->MPType == 0)
  { // initialize
    double currNorm, *norms = (double *)bmalloc(*num_points * sizeof(double));
    comp_d gamma_out_d, gamma_in_d, sum_trace;
    vec_d v_out_d, v_in_d;
    point_d *Pts = (point_d *)bmalloc(*num_points * sizeof(point_d)), dehom;
    endpoint_data_d endPt;

    // initialize
    size = sliceMover->B_d->rows;
    init_vec_d(v_out_d, size);
    init_vec_d(v_in_d, size);
    init_point_d(dehom, 0);
    init_endpoint_data_d(&endPt);

    // set top position of file
    fgetpos(WITPTS, &witPositionTop);
    
    // print number of points
    fprintf(WITPTS, "%d                                    \n\n", *num_points);

    // setup norms & Pts
    for (i = 0; i < *num_points; i++)
    {
      init_point_d(Pts[i], 0);
      point_cp_d(Pts[i], currPoints_d[i]);
      norms[i] = infNormVec_d(Pts[i]);

      // print dehomogenized points
      witnessFindDehom_d(dehom, Pts[i], W, 0);
      for (j = 0; j < dehom->size; j++)
        print_comp_out_d(WITPTS, &dehom->coord[j]);
      fprintf(WITPTS, "\n"); 
    }

    // loop 
    while (its < maxIts && !isComplete)
    { // display message
      printf("Performing monodromy loops starting with %d point%s\n", *num_points, *num_points == 1 ? "" : "s");

      // generate random data
      get_comp_rand_d(gamma_out_d);
      get_comp_rand_d(gamma_in_d);
      v_out_d->size = v_in_d->size = size;
      make_vec_random_d(v_out_d, size);
      for (i = 0; i < size; i++)
        set_zero_d(&v_in_d->coord[i]);
      
      // do a loop for each point
      for (i = 0; i < *num_points; i++)
      { // finish the setup for sliceMover
        initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
        final_setup_slice_moving(sliceMover, sliceMover->Prog, T->MPType, T->AMP_max_prec, 0);

        // perform the monodromy loop
        rV = basicMonodromyLoop(&monodromyPt, sliceMover, sliceMover->Prog, Pts[i], NULL, 52, W, i, T, OUT, MIDOUT, gamma_out_d, NULL, NULL, gamma_in_d, NULL, NULL, v_out_d, NULL, NULL, v_in_d, NULL, NULL);

        // determine if this is a new point
        if (!rV)
        { // compute its norm
          currNorm = infNormVec_d(monodromyPt.PD_d.point);

          // compare with other points to determine if new
          matchingPathNum = -1;
          for (j = 0; j < *num_points; j++)
            if (fabs(currNorm - norms[j]) < tol)
            { // compute difference
              findDiff_point(normDiff, monodromyPt.PD_d.point, NULL, 52, Pts[j], NULL, 52);
              if (mpf_cmp_d(normDiff, tol) < 0)
              { // found a match
                matchingPathNum = j;
                break;
              }
            }
  
          // see if new
          if (matchingPathNum < 0)
          { // update structures with new point
            Pts = (point_d *)brealloc(Pts, (*num_points + 1) * sizeof(point_d));
            *traces_d = (comp_d *)brealloc(*traces_d, (*num_points + 1) * sizeof(comp_d));
            norms = (double *)brealloc(norms, (*num_points + 1) * sizeof(double));

            init_point_d(Pts[*num_points], 0);
            point_cp_d(Pts[*num_points], monodromyPt.PD_d.point); 
            norms[*num_points] = currNorm;
            point_cp_d(endPt.endPt, monodromyPt.PD_d.point);    

            // print to dehom point to file
            witnessFindDehom_d(dehom, Pts[*num_points], W, 0);
            for (j = 0; j < dehom->size; j++)
              print_comp_out_d(WITPTS, &dehom->coord[j]);
            fprintf(WITPTS, "\n");

            // set the bottom position of file
            fgetpos(WITPTS, &witPositionBottom);

            // move to the top and print the number of points
            fsetpos(WITPTS, &witPositionTop);
            fprintf(WITPTS, "%d", *num_points + 1);

            // move back to the bottom
            fsetpos(WITPTS, &witPositionBottom);            

            // compute the new trace and increment the number of points
            calculateTrace((*traces_d)[*num_points], NULL, NULL, sliceMover, sliceMover->Prog, &endPt, NULL, NULL, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);
            (*num_points)++;

            // add up the traces
            set_zero_d(sum_trace);
            for (j = 0; j < *num_points; j++)
              add_d(sum_trace, sum_trace, (*traces_d)[j]);

            printf("Trace residual for %d point%s: %e\n", *num_points, *num_points == 1 ? "" : "s", d_abs_d(sum_trace));
            // check for completeness
            if (d_abs_d(sum_trace) < tol)
            { // complete!
              isComplete = 1;
              break;
            }
          }
        }
      }

      // increment the number of monodromy loops    
      its++;  
    }

    // save points to W
    if (*num_points > W->codim[0].num_set)
    { // add new points 
      W->codim[0].witnessPts_d = (endpoint_data_d *)brealloc(W->codim[0].witnessPts_d, *num_points * sizeof(endpoint_data_d));
      for (i = W->codim[0].num_set; i < *num_points; i++)
      {
        init_endpoint_data_d(&W->codim[0].witnessPts_d[i]);
        vec_cp_d(W->codim[0].witnessPts_d[i].endPt, Pts[i]);
        vec_cp_d(W->codim[0].witnessPts_d[i].last_approx, Pts[i]);
        set_zero_d(W->codim[0].witnessPts_d[i].finalT);
        W->codim[0].witnessPts_d[i].cond_num = 1;
        W->codim[0].witnessPts_d[i].corank = 0;
        W->codim[0].witnessPts_d[i].smallest_nonzero_SV = 1;
        W->codim[0].witnessPts_d[i].largest_zero_SV = 0;
        W->codim[0].witnessPts_d[i].retVal = 0;
      }

      // setup other data
      W->codim[0].num_set = *num_points;
      W->codim[0].num_nonsing = *num_points;
      W->codim[0].num_sing = 0;
      W->codim[0].num_components = 1;
      W->codim[0].witnessPt_types = (int *)brealloc(W->codim[0].witnessPt_types, *num_points * sizeof(int));
      W->codim[0].component_nums = (int *)brealloc(W->codim[0].component_nums, *num_points * sizeof(int));
      W->codim[0].multiplicities = (int *)brealloc(W->codim[0].multiplicities, *num_points * sizeof(int));
      W->codim[0].deflations_needed = (int *)brealloc(W->codim[0].deflations_needed, *num_points * sizeof(int));
      for (i = 0; i < *num_points; i++)
      {
        W->codim[0].witnessPt_types[i] = NON_SINGULAR;
        W->codim[0].component_nums[i] = W->codim[0].deflations_needed[i] = 0;
        W->codim[0].multiplicities[i] = 1;
      }
    }

    // clear
    clear_endpoint_data_d(&endPt);
    clear_vec_d(v_out_d);
    clear_vec_d(v_in_d);    
    for (i = 0; i < *num_points; i++)
      clear_point_d(Pts[i]);
    clear_point_d(dehom);
    free(Pts);
    free(norms);
  }
  else if (T->MPType == 1)
  { // initialize
    mpf_t currNorm, *norms = (mpf_t *)bmalloc(*num_points * sizeof(mpf_t));
    comp_mp gamma_out_mp, gamma_in_mp, sum_trace;
    vec_mp v_out_mp, v_in_mp;
    point_mp *Pts = (point_mp *)bmalloc(*num_points * sizeof(point_mp)), dehom;
    endpoint_data_mp endPt;

    // initialize
    mpf_init(currNorm);
    init_mp(gamma_out_mp);
    init_mp(gamma_in_mp);
    init_mp(sum_trace);
    size = sliceMover->B_mp->rows;
    init_vec_mp(v_out_mp, size);
    init_vec_mp(v_in_mp, size);
    init_point_mp(dehom, 0);
    init_endpoint_data_mp(&endPt);

    // set top position of file
    fgetpos(WITPTS, &witPositionTop);

    // print number of points
    fprintf(WITPTS, "%d                                    \n\n", *num_points);

    // setup norms & Pts
    for (i = 0; i < *num_points; i++)
    {
      init_point_mp(Pts[i], 0);
      point_cp_mp(Pts[i], currPoints_mp[i]);
      mpf_init(norms[i]);
      infNormVec_mp2(norms[i], Pts[i]);

      // print dehomogenized points
      witnessFindDehom_mp(dehom, Pts[i], W, 0, T->Precision);
      for (j = 0; j < dehom->size; j++)
        print_comp_out_mp(WITPTS, &dehom->coord[j]);
      fprintf(WITPTS, "\n");
    }
    
    // loop 
    while (its < maxIts && !isComplete)
    { // display message
      printf("Performing monodromy loops starting with %d point%s\n", *num_points, *num_points == 1 ? "" : "s");

      // generate random data
      get_comp_rand_mp(gamma_out_mp);
      get_comp_rand_mp(gamma_in_mp);
      v_out_mp->size = v_in_mp->size = size;
      make_vec_random_mp(v_out_mp, size);
      for (i = 0; i < size; i++)
        set_zero_mp(&v_in_mp->coord[i]);
      
      // do a loop for each point
      for (i = 0; i < *num_points; i++)
      { // finish the setup for sliceMover
        initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
        final_setup_slice_moving(sliceMover, sliceMover->Prog, T->MPType, T->AMP_max_prec, 0);

        // perform the monodromy loop
        rV = basicMonodromyLoop(&monodromyPt, sliceMover, sliceMover->Prog, NULL, Pts[i], T->Precision, W, i, T, OUT, MIDOUT, NULL, gamma_out_mp, NULL, NULL, gamma_in_mp, NULL, NULL, v_out_mp, NULL, NULL, v_in_mp, NULL);

        // determine if this is a new point
        if (!rV)
        { // compute its norm
          infNormVec_mp2(currNorm, monodromyPt.PD_mp.point);

          // compare with other points to determine if new
          matchingPathNum = -1;
          for (j = 0; j < *num_points; j++)
          {
            mpf_sub(normDiff, currNorm, norms[j]);
            mpf_abs(normDiff, normDiff);
            if (mpf_cmp_d(normDiff, tol) < 0)
            { // compute difference
              findDiff_point(normDiff, NULL, monodromyPt.PD_mp.point, T->Precision, NULL, Pts[j], T->Precision);
              if (mpf_cmp_d(normDiff, tol) < 0)
              { // found a match
                matchingPathNum = j;
                break;
              }
            }
          }
  
          // see if new
          if (matchingPathNum < 0)
          { // update structures with new point
            Pts = (point_mp *)brealloc(Pts, (*num_points + 1) * sizeof(point_mp));
            *traces_mp = (comp_mp *)brealloc(*traces_mp, (*num_points + 1) * sizeof(comp_mp));
            norms = (mpf_t *)brealloc(norms, (*num_points + 1) * sizeof(mpf_t));

            mpf_init(norms[*num_points]);
            init_mp((*traces_mp)[*num_points]);
            init_point_mp(Pts[*num_points], 0);
            point_cp_mp(Pts[*num_points], monodromyPt.PD_mp.point); 
            mpf_set(norms[*num_points], currNorm);
            point_cp_mp(endPt.endPt, monodromyPt.PD_mp.point);    

            // print to dehom point to file
            witnessFindDehom_mp(dehom, Pts[*num_points], W, 0, T->Precision);
            for (j = 0; j < dehom->size; j++)
              print_comp_out_mp(WITPTS, &dehom->coord[j]);
            fprintf(WITPTS, "\n");

            // set the bottom position of file
            fgetpos(WITPTS, &witPositionBottom);

            // move to the top and print the number of points
            fsetpos(WITPTS, &witPositionTop);
            fprintf(WITPTS, "%d", *num_points + 1);

            // move back to the bottom
            fsetpos(WITPTS, &witPositionBottom);

            // compute the new trace and increment the number of points
            calculateTrace(NULL, (*traces_mp)[*num_points], NULL, sliceMover, sliceMover->Prog, NULL, &endPt, NULL, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);
            (*num_points)++;

            // add up the traces
            set_zero_mp(sum_trace);
            for (j = 0; j < *num_points; j++)
              add_mp(sum_trace, sum_trace, (*traces_mp)[j]);
            
            printf("Trace residual for %d point%s: %e\n", *num_points, *num_points == 1 ? "" : "s", d_abs_mp(sum_trace));
            // check for completeness
            if (d_abs_mp(sum_trace) < tol)
            { // complete!
              isComplete = 1;
              break;
            }
          }
        }
      }

      // increment the number of monodromy loops    
      its++;  
    }

    // save points to W
    if (*num_points > W->codim[0].num_set)
    { // add new points 
      W->codim[0].witnessPts_mp = (endpoint_data_mp *)brealloc(W->codim[0].witnessPts_mp, *num_points * sizeof(endpoint_data_mp));
      for (i = 0; i < *num_points; i++)
      {
        init_endpoint_data_mp(&W->codim[0].witnessPts_mp[i]);
        vec_cp_mp(W->codim[0].witnessPts_mp[i].endPt, Pts[i]);
        vec_cp_mp(W->codim[0].witnessPts_mp[i].last_approx, Pts[i]);
        set_zero_mp(W->codim[0].witnessPts_mp[i].finalT);
        W->codim[0].witnessPts_mp[i].cond_num = 1;
        W->codim[0].witnessPts_mp[i].corank = 0;
        W->codim[0].witnessPts_mp[i].smallest_nonzero_SV = 1;
        W->codim[0].witnessPts_mp[i].largest_zero_SV = 0;
        W->codim[0].witnessPts_mp[i].retVal = 0;
      }

      // setup other data
      W->codim[0].num_set = *num_points;
      W->codim[0].num_nonsing = *num_points;
      W->codim[0].num_sing = 0;
      W->codim[0].num_components = 1;
      W->codim[0].witnessPt_types = (int *)brealloc(W->codim[0].witnessPt_types, *num_points * sizeof(int));
      W->codim[0].component_nums = (int *)brealloc(W->codim[0].component_nums, *num_points * sizeof(int));
      W->codim[0].multiplicities = (int *)brealloc(W->codim[0].multiplicities, *num_points * sizeof(int));
      W->codim[0].deflations_needed = (int *)brealloc(W->codim[0].deflations_needed, *num_points * sizeof(int));
      for (i = 0; i < *num_points; i++)
      {
        W->codim[0].witnessPt_types[i] = NON_SINGULAR;
        W->codim[0].component_nums[i] = W->codim[0].deflations_needed[i] = 0;
        W->codim[0].multiplicities[i] = 1;
      }
    }

    // clear
    clear_endpoint_data_mp(&endPt);
    mpf_clear(currNorm);
    clear_mp(gamma_out_mp);
    clear_mp(gamma_in_mp);
    clear_mp(sum_trace);
    clear_vec_mp(v_out_mp);
    clear_vec_mp(v_in_mp);    
    for (i = 0; i < *num_points; i++)
    {
      mpf_clear(norms[i]);
      clear_point_mp(Pts[i]);
    }
    clear_point_mp(dehom);
    free(Pts);
    free(norms);
  }
  else 
  { // initialize
    int *prec = (int *)bmalloc(*num_points * sizeof(int));
    double currNorm_d = 0, *norms_d = (double *)bmalloc(*num_points * sizeof(double));
    mpf_t currNorm_mp, *norms_mp = (mpf_t *)bmalloc(*num_points * sizeof(mpf_t));
    comp_d gamma_out_d, gamma_in_d, sum_trace_d;
    comp_mp gamma_out_mp, gamma_in_mp, sum_trace_mp;
    mpq_t gamma_out_rat[2], gamma_in_rat[2];
    vec_d v_out_d, v_in_d;
    vec_mp v_out_mp, v_in_mp;
    mpq_t **v_out_rat, **v_in_rat;
    point_d *Pts_d = (point_d *)bmalloc(*num_points * sizeof(point_d)), dehom_d;
    point_mp *Pts_mp = (point_mp *)bmalloc(*num_points * sizeof(point_mp)), dehom_mp;
    endpoint_data_amp endPt;

    // initialize
    mpf_init2(currNorm_mp, T->Precision);
    init_mp2(gamma_out_mp, T->Precision);
    init_mp2(gamma_in_mp, T->Precision);
    init_mp2(sum_trace_mp, T->Precision);
    init_rat(gamma_out_rat);
    init_rat(gamma_in_rat);
    size = sliceMover->B_d->rows;
    init_vec_d(v_out_d, size);
    init_vec_d(v_in_d, size);
    init_vec_mp2(v_out_mp, size, T->Precision);
    init_vec_mp2(v_in_mp, size, T->Precision);
    init_vec_rat(v_out_rat, size);
    init_vec_rat(v_in_rat, size);
    init_endpoint_data_amp(&endPt, T->Precision, T->Precision);
    init_point_d(dehom_d, 0);
    init_point_mp(dehom_mp, 0);

    // set top position of file
    fgetpos(WITPTS, &witPositionTop);

    // print number of points
    fprintf(WITPTS, "%d                                    \n\n", *num_points);

    // setup norms & Pts
    for (i = 0; i < *num_points; i++)
    {
      prec[i] = currPoints_prec[i];
      if (prec[i] < 64)
      { // setup _d
        init_point_d(Pts_d[i], 0);
        point_cp_d(Pts_d[i], currPoints_d[i]);
        norms_d[i] = infNormVec_d(Pts_d[i]);

        // print dehomogenized points
        witnessFindDehom_d(dehom_d, Pts_d[i], W, 0);
        for (j = 0; j < dehom_d->size; j++)
          print_comp_out_d(WITPTS, &dehom_d->coord[j]);
        fprintf(WITPTS, "\n");
      }
      else
      { // setup _mp
        mpf_init2(norms_mp[i], prec[i]);
        init_point_mp2(Pts_mp[i], 0, prec[i]);
        point_cp_mp(Pts_mp[i], currPoints_mp[i]);
        infNormVec_mp2(norms_mp[i], Pts_mp[i]);

        // print dehomogenized points
        setprec_point_mp(dehom_mp, prec[i]);
        witnessFindDehom_mp(dehom_mp, Pts_mp[i], W, 0, prec[i]);
        for (j = 0; j < dehom_mp->size; j++)
          print_comp_out_mp(WITPTS, &dehom_mp->coord[j]);
        fprintf(WITPTS, "\n");
      }
    }
    
    // loop 
    while (its < maxIts && !isComplete)
    { // display message
      printf("Performing monodromy loops starting with %d point%s\n", *num_points, *num_points == 1 ? "" : "s");

      // generate random data
      get_comp_rand_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, T->Precision, T->AMP_max_prec, 0, 0);
      get_comp_rand_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, T->Precision, T->AMP_max_prec, 0, 0);
      v_out_d->size = v_in_d->size = v_out_mp->size = v_in_mp->size = size;
      make_vec_random_rat(v_out_d, v_out_mp, v_out_rat, size, T->Precision, T->AMP_max_prec, 0, 0);
      for (i = 0; i < size; i++)
      {
        set_zero_d(&v_in_d->coord[i]);
        set_zero_mp(&v_in_mp->coord[i]);
        set_zero_rat(v_in_rat[i]);
      }
      
      // do a loop for each point
      for (i = 0; i < *num_points; i++)
      { // finish the setup for sliceMover
        initialize_slice_moving_sliceVec(sliceMover, size, T->MPType);
        final_setup_slice_moving(sliceMover, sliceMover->Prog, T->MPType, T->AMP_max_prec, 0);

        // perform the monodromy loop
        rV = basicMonodromyLoop(&monodromyPt, sliceMover, sliceMover->Prog, Pts_d[i], Pts_mp[i], prec[i], W, i, T, OUT, MIDOUT, gamma_out_d, gamma_out_mp, gamma_out_rat, gamma_in_d, gamma_in_mp, gamma_in_rat, v_out_d, v_out_mp, v_out_rat, v_in_d, v_in_mp, v_in_rat);

        // determine if this is a new point
        if (!rV)
        { // determine precision
          if (monodromyPt.prec < 64) 
          { // use _d
            currNorm_d = infNormVec_d(monodromyPt.PD_d.point);
          }
          else
          { // use _mp
            mpfr_set_prec(currNorm_mp, monodromyPt.prec);
            infNormVec_mp2(currNorm_mp, monodromyPt.PD_mp.point);
          }

          // compare with other points to determine if new
          matchingPathNum = -1;
          for (j = 0; j < *num_points; j++)
            if (prec[j] < 64 && monodromyPt.prec < 64)
            { // all in _d
              if (fabs(currNorm_d - norms_d[j]) < tol)
              { // compute difference
                findDiff_point(normDiff, monodromyPt.PD_d.point, monodromyPt.PD_mp.point, monodromyPt.prec, Pts_d[j], Pts_mp[j], prec[j]);
                if (mpf_cmp_d(normDiff, tol) < 0)
                { // found a match
                  matchingPathNum = j;
                  break;
                }
              }
            }
            else if (prec[j] < 64) // monodromyPt.prec >= 64
            { // increase to _mp
              mpfr_set_prec(normDiff, monodromyPt.prec);
              mpf_set_d(normDiff, norms_d[j]);
              mpf_sub(normDiff, normDiff, currNorm_mp);
              mpf_abs(normDiff, normDiff);

              if (mpf_cmp_d(normDiff, tol) < 0)
              { // compute difference
                findDiff_point(normDiff, monodromyPt.PD_d.point, monodromyPt.PD_mp.point, monodromyPt.prec, Pts_d[j], Pts_mp[j], prec[j]);
                if (mpf_cmp_d(normDiff, tol) < 0)
                { // found a match
                  matchingPathNum = j;
                  break;
                }
              }
            }
            else if (monodromyPt.prec < 64) // prec[j] >= 64
            { // increase to _mp
              mpfr_set_prec(normDiff, prec[j]);
              mpf_set_d(normDiff, currNorm_d);
              mpf_sub(normDiff, normDiff, norms_mp[j]);
              mpf_abs(normDiff, normDiff);

              if (mpf_cmp_d(normDiff, tol) < 0)
              { // compute difference
                findDiff_point(normDiff, monodromyPt.PD_d.point, monodromyPt.PD_mp.point, monodromyPt.prec, Pts_d[j], Pts_mp[j], prec[j]);
                if (mpf_cmp_d(normDiff, tol) < 0)
                { // found a match
                  matchingPathNum = j;
                  break;
                }
              }
            }
            else
            { // both in _mp
              mpfr_set_prec(normDiff, MAX(prec[j], monodromyPt.prec));
              mpf_sub(normDiff, norms_mp[j], currNorm_mp);
              mpf_abs(normDiff, normDiff);

              if (mpf_cmp_d(normDiff, tol) < 0)
              { // compute difference
                findDiff_point(normDiff, monodromyPt.PD_d.point, monodromyPt.PD_mp.point, monodromyPt.prec, Pts_d[j], Pts_mp[j], prec[j]);
                if (mpf_cmp_d(normDiff, tol) < 0)
                { // found a match
                  matchingPathNum = j;
                  break;
                }
              }
            }

          // see if new
          if (matchingPathNum < 0)
          { // update structures with new point
            Pts_d = (point_d *)brealloc(Pts_d, (*num_points + 1) * sizeof(point_d));
            Pts_mp = (point_mp *)brealloc(Pts_mp, (*num_points + 1) * sizeof(point_mp));
            *traces_prec = (int *)brealloc(*traces_prec, (*num_points + 1) * sizeof(int));
            *traces_d = (comp_d *)brealloc(*traces_d, (*num_points + 1) * sizeof(comp_d));
            *traces_mp = (comp_mp *)brealloc(*traces_mp, (*num_points + 1) * sizeof(comp_mp));
            prec = (int *)brealloc(prec, (*num_points + 1) * sizeof(int));
            norms_d = (double *)brealloc(norms_d, (*num_points + 1) * sizeof(double));
            norms_mp = (mpf_t *)brealloc(norms_mp, (*num_points + 1) * sizeof(mpf_t));

            init_mp((*traces_mp)[*num_points]);

            prec[*num_points] = endPt.curr_prec = monodromyPt.prec;
            if (prec[*num_points] < 64)
            { // save to _d
              init_point_d(Pts_d[*num_points], 0);
              norms_d[*num_points] = currNorm_d;
              point_cp_d(Pts_d[*num_points], monodromyPt.PD_d.point);
              point_cp_d(endPt.endPt_d, monodromyPt.PD_d.point);

              // print to dehom point to file
              witnessFindDehom_d(dehom_d, Pts_d[*num_points], W, 0);
              for (j = 0; j < dehom_d->size; j++)
                print_comp_out_d(WITPTS, &dehom_d->coord[j]);
              fprintf(WITPTS, "\n");
            }
            else
            { // save to _mp
              mpf_init2(norms_mp[*num_points], prec[*num_points]);
              init_point_mp2(Pts_mp[*num_points], 0, prec[*num_points]);
              setprec_point_mp(endPt.endPt_mp, prec[*num_points]);
              point_cp_mp(Pts_mp[*num_points], monodromyPt.PD_mp.point); 
              mpf_set(norms_mp[*num_points], currNorm_mp);
              point_cp_mp(endPt.endPt_mp, monodromyPt.PD_mp.point);    

              // print to dehom point to file
              witnessFindDehom_mp(dehom_mp, Pts_mp[*num_points], W, 0, prec[*num_points]);
              for (j = 0; j < dehom_mp->size; j++)
                print_comp_out_mp(WITPTS, &dehom_mp->coord[j]);
              fprintf(WITPTS, "\n");
            }

            // set the bottom position of file
            fgetpos(WITPTS, &witPositionBottom);

            // move to the top and print the number of points
            fsetpos(WITPTS, &witPositionTop);
            fprintf(WITPTS, "%d", *num_points + 1);

            // move back to the bottom
            fsetpos(WITPTS, &witPositionBottom);

            // compute the new trace and increment the number of points
            calculateTrace((*traces_d)[*num_points], (*traces_mp)[*num_points], &(*traces_prec)[*num_points], sliceMover, sliceMover->Prog, NULL, NULL, &endPt, W->codim[0].codim, 0, T, OUT, MIDOUT, trace_proj_d, trace_proj_mp, trace_proj_rat, trace_slice_d, trace_slice_mp, trace_slice_rat, s_d, s_mp, s_rat, 1);

            (*num_points)++;

            *amp_trace_prec = MAX(*amp_trace_prec, (*traces_prec)[*num_points - 1]);

            // compute sum
            if (*amp_trace_prec < 64)
            { // all in _d
              set_zero_d(sum_trace_d);
              for (j = 0; j < *num_points; j++)
                add_d(sum_trace_d, sum_trace_d, (*traces_d)[j]);

              printf("Trace residual for %d point%s: %e\n", *num_points, *num_points == 1 ? "" : "s", d_abs_d(sum_trace_d));
              // check for completeness
              if (d_abs_d(sum_trace_d) < tol)
              { // complete!
                isComplete = 1;
                break;
              }
            }
            else
            { // increase all precision and add them up
              setprec_mp(sum_trace_mp, *amp_trace_prec);
              set_zero_mp(sum_trace_mp);
              for (j = 0; j < *num_points; j++)
              { // increase prec
                if ((*traces_prec)[j] < 64)
                { // increase to _mp
                  setprec_mp((*traces_mp)[j], *amp_trace_prec);
                  d_to_mp((*traces_mp)[j], (*traces_d)[j]);
                  (*traces_prec)[j] = *amp_trace_prec;
                }
                else
                { // set correct precision
                  change_prec_mp((*traces_mp)[j], *amp_trace_prec);
                }
                add_mp(sum_trace_mp, sum_trace_mp, (*traces_mp)[j]);
              }

              printf("Trace residual for %d point%s: %e\n", *num_points, *num_points == 1 ? "" : "s", d_abs_mp(sum_trace_mp));
              // check for completeness
              if (d_abs_mp(sum_trace_mp) < tol)
              { // complete!
                isComplete = 1;
                break;
              }
            }
          }
        }
      }

      // increment the number of monodromy loops    
      its++;  
    }

    // save points to W
    if (*num_points > W->codim[0].num_set)
    { // add new points 
      W->codim[0].witnessPts_amp = (endpoint_data_amp *)brealloc(W->codim[0].witnessPts_amp, *num_points * sizeof(endpoint_data_amp));
      for (i = W->codim[0].num_set; i < *num_points; i++)
      {
        init_endpoint_data_amp(&W->codim[0].witnessPts_amp[i], prec[i], prec[i]);
        if (prec[i] < 64)
        { // set to _d
          vec_cp_d(W->codim[0].witnessPts_amp[i].endPt_d, Pts_d[i]);
          vec_cp_d(W->codim[0].witnessPts_amp[i].last_approx_d, Pts_d[i]);
          set_zero_d(W->codim[0].witnessPts_amp[i].finalT_d);
        }
        else
        { // set to _mp
          vec_cp_mp(W->codim[0].witnessPts_amp[i].endPt_mp, Pts_mp[i]);
          vec_cp_mp(W->codim[0].witnessPts_amp[i].last_approx_mp, Pts_mp[i]);
          set_zero_mp(W->codim[0].witnessPts_amp[i].finalT_mp);
        }

        W->codim[0].witnessPts_amp[i].cond_num = 1;
        W->codim[0].witnessPts_amp[i].corank = 0;
        W->codim[0].witnessPts_amp[i].smallest_nonzero_SV = 1;
        W->codim[0].witnessPts_amp[i].largest_zero_SV = 0;
        W->codim[0].witnessPts_amp[i].retVal = 0;
      }

      // setup other data
      W->codim[0].num_set = *num_points;
      W->codim[0].num_nonsing = *num_points;
      W->codim[0].num_sing = 0;
      W->codim[0].num_components = 1;
      W->codim[0].witnessPt_types = (int *)brealloc(W->codim[0].witnessPt_types, *num_points * sizeof(int));
      W->codim[0].component_nums = (int *)brealloc(W->codim[0].component_nums, *num_points * sizeof(int));
      W->codim[0].multiplicities = (int *)brealloc(W->codim[0].multiplicities, *num_points * sizeof(int));
      W->codim[0].deflations_needed = (int *)brealloc(W->codim[0].deflations_needed, *num_points * sizeof(int));
      for (i = 0; i < *num_points; i++)
      {
        W->codim[0].witnessPt_types[i] = NON_SINGULAR;
        W->codim[0].component_nums[i] = W->codim[0].deflations_needed[i] = 0;
        W->codim[0].multiplicities[i] = 1;
      }
    }

    // clear
    clear_endpoint_data_amp(&endPt);
    mpf_clear(currNorm_mp);
    clear_d_mp_rat(gamma_out_d, gamma_out_mp, gamma_out_rat, T->MPType);
    clear_d_mp_rat(gamma_in_d, gamma_in_mp, gamma_in_rat, T->MPType);
    clear_mp(sum_trace_mp);
    clear_vec(v_out_d, v_out_mp, v_out_rat, T->MPType);
    clear_vec(v_in_d, v_in_mp, v_in_rat, T->MPType);    
    for (i = 0; i < *num_points; i++)
      if (prec[i] < 64)
      { 
        clear_point_d(Pts_d[i]);
      }
      else
      {
        mpf_clear(norms_mp[i]);
        clear_point_mp(Pts_mp[i]);
      }
    free(Pts_d);
    free(Pts_mp);
    free(norms_d);
    free(norms_mp);
    free(prec);
  }  

  clear_endgame_data(&monodromyPt);
  mpf_clear(normDiff);
  fclose(WITPTS);

  return !isComplete;
}

