// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

// Given a general witness set for a generically reduced irreducible component
// and a projection, compute a witness set for the projection.

// Given a witness set for a projection, perform a membership test.

void witnessProjectionMenu(witness_t *W, tracker_config_t *T, int pathMod);
void witnessProjection(witness_t *W, tracker_config_t *T, int pathMod, int codim_index, int component_number, char *outName, char *projName);
int computeFiberDim(witness_t *W, tracker_config_t *T, int codim_index, int component_number, int pathNum, int *projection);
void computeProjectionWitnessSet(witness_t *W_proj, witness_t *W, tracker_config_t *T, int codim_index, int component_number, int degree, int *pathNums, int *projection, int proj_dim, int fiber_dim, int *proj_deg, int *fiber_deg, int pathMod, trackingStats *trackCount);

void witnessProjectionMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: main control function for generating a witness set for *
*  the projection of an irreducible & gen. reduced component    *
\***************************************************************/
{
  int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  tracker_config_t T;
  witness_t witnessSet;

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
  { // exit since parallel witness projection is not implemented
#ifdef _HAVE_MPI
    printf("ERROR: Parallel witness set projection is not implemented. Please use sequential version!\n");
    bexit(ERROR_OTHER);
#endif
  }

  // setup witnessSet
  setupWitnessDataFromFile("witness_data", "witness_data_old", "preproc_data", "deg.out", &witnessSet, &T, 1);

  // error checking
  if (witnessSet.Prog->numPathVars > 0)
  { // path variable present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (witnessSet.Prog->numPars > 0)
  { // parameter present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (witnessSet.PPD.num_hom_var_gp)
  { // projective coordinates not allowe
    rename("witness_data_old", "witness_data");

    printf("Witness set projection is only implemented for affine varieties (variable_group)!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup the number of variables
  T.numVars = witnessSet.Prog->numVars;

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

  // create output files now that everything is setup
  numIrredDecompOutput(&witnessSet, &T, 5, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 5

  // ask the user which component to project and do the actual projection
  witnessProjectionMenu(&witnessSet, &T, pathMod);

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void witnessProjectionMenu(witness_t *W, tracker_config_t *T, int pathMod)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: displays a menu of components to project and calls the *
*  function to do the projection                                *
\***************************************************************/
{
  int i, j, codim_index, min_deg, max_deg, count, rV, dim_number, component_number, selection_made, gen_reduced;
  int *degrees = NULL, *dim = (int *)bmalloc(W->num_codim * sizeof(int)), *codim_good = (int *)bmalloc(W->num_codim * sizeof(int));
  char ch, *tempStr = NULL, *outputFile = NULL;

  // find the dimension and good dimensions, and make sure one exists
  rV = 0;
  for (codim_index = 0; codim_index < W->num_codim; codim_index++)
  { // see if there are classified components for this codim
    codim_good[codim_index] = (W->codim[codim_index].num_components > 0 ? 1 : 0);

    if (codim_good[codim_index])
      rV = 1;

    // determine the dimension of this codim
    dim[codim_index] = W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
  }

  if (!rV)
  { // there are no classified components!!
    printf("\nThere are no classified components to sample!\n\n");
    free(dim);
    free(codim_good);
    return;
  }

  // so we have atleast one classified component
  do 
  { // initialize
    selection_made = 0;

    // print title
    printf("\n\n*************** Components to Project ****************\n\n");

    // display a catalog of the available components in each codim
    for (codim_index = 0; codim_index < W->num_codim; codim_index++)
    { // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      // find the minimum and maximum degree
      max_deg = 0;
      min_deg = W->codim[codim_index].num_set + 1;
      for (i = 0; i < W->codim[codim_index].num_components; i++)
      { // find maximum degree
        if (max_deg < degrees[i])
          max_deg = degrees[i];

        // find minimum degree
        if (min_deg > degrees[i])
          min_deg = degrees[i];
      }

      if (W->codim[codim_index].num_components > 0)
      {
        printf("Dimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
        if (W->codim[codim_index].num_components == 1)
          printf("\n");
        else
          printf("s\n");
        printf("-----------------------------------------------------\n");

        // display the summary
        for (i = min_deg; i <= max_deg; i++)
        { // count the number that have degree == i
          count = 0;
          for (j = 0; j < W->codim[codim_index].num_components; j++)
            if (degrees[j] == i)
              count++;

          if (count > 0)
          { // display the number
            printf("   degree %d: %d component", i, count);
            if (count == 1)
              printf("\n");
            else
              printf("s\n");
          }
        }
        printf("\n");
      }
    }

    printf("\nPlease select a dimension to project (-1 to quit): ");
    rV = scanf("%d", &dim_number);

    if (rV < 0)
    { // at EOF - so we need to quit
      dim_number = -1;
      selection_made = 1;
    }
    else
    { // we are not at EOF - flush the buffer
      do
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");
        selection_made = 0;
      }
      else if (dim_number == -1)
      { // dim_number is valid
        selection_made = 1;
      }
      else
      { // verify that the dimension is one of the good ones
        selection_made = 0;
        for (j = 0; j < W->num_codim; j++)
          if (codim_good[j] && dim[j] == dim_number)
          { // selection is valid
            selection_made = 1;
            break;
          }
          
        if (!selection_made) 
        { // dim_number is not valid
          printf("\nThe dimension %d is not valid!\n", dim_number);
        }
      }
    }
  } while (!selection_made);

  // so, either dim_number == -1 OR dim_number is one of the ones that have a classified component
  if (dim_number != -1)
  { // find the codim_index for the selected dimension
    codim_index = 0;
    for (i = 0; i < W->num_codim; i++)
    { // see if the dimension agrees
      if (dim[i] == dim_number)
      { // store the index and exit loop
        codim_index = i;
        break;
      }
    }

    do
    { // initialize
      selection_made = 0;

      // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      printf("\nDimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
      if (W->codim[codim_index].num_components == 1)
        printf("\n");
      else
        printf("s\n");
      printf("-----------------------------------------------------\n");

      // display the summary of components
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        printf("   component %d has degree %d\n", i, degrees[i]);
      printf("\n");

      printf("\nPlease select a component to project (-1 to quit): ");
      rV = scanf("%d", &component_number);

      if (rV < 0)
      { // at EOF - so we need to quit
        component_number = -1;
        selection_made = 1;
      }
      else
      { // we are not at EOF - flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (rV == 0)
        { // invalid input
          printf("\nThe input was not read in correctly!\n");
          selection_made = 0;
        }
        else if (component_number != -1 && (component_number < 0 || component_number >= W->codim[codim_index].num_components))
        { // component_number is not valid
          printf("\nThe component %d is not valid!\n", component_number);
          selection_made = 0;
        }
        else
        { // component_number is valid
          selection_made = 1;
        }
      }
    } while (!selection_made);

    // so, either component_number == -1 OR 0 <= component_number < num_components
    if (0 <= component_number && component_number < W->codim[codim_index].num_components)
    { // verify that every point on the component has multiplicity 1 & required no deflations (should be equivalent!)
      gen_reduced = 1;
      for (i = 0; i < W->codim[codim_index].num_set && gen_reduced; i++)
        if (W->codim[codim_index].component_nums[i] == component_number)
        { // this point lies on component - verify multiplicity and deflations
          if (W->codim[codim_index].multiplicities[i] != 1)
          { // error
            printf("\nComponent %d of dimension %d is not generically reduced!\n", component_number, W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp);
            gen_reduced = 0;
          }
          else if (W->codim[codim_index].deflations_needed[i] != 0)
          { // error
            printf("\nComponent %d of dimension %d is not generically reduced!\n", component_number, W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp);
            gen_reduced = 0;
          }
        }

      if (gen_reduced)
      { // perform the projection
        witnessProjection(W, T, pathMod, codim_index, component_number, "witness_data_projection", "projection");
      }
    }
  }

  free(tempStr);
  free(outputFile);

  return;
}

void witnessProjection(witness_t *W, tracker_config_t *T, int pathMod, int codim_index, int component_number, char *outName, char *projName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform the projection - Assume the codim exists!      *
\***************************************************************/
{
  int i, rV, degree = 0, num_set = W->codim[codim_index].num_set, dim = W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp;
  int num_fiber_vars = 0, num_affine_vars = W->orig_variables - W->PPD.num_var_gp;
  int fiber_dim = 0, projection_dim = 0, fiber_degree = 0, projection_degree = 0;
  int *pathNums = NULL, *projVars = (int *)bmalloc(num_affine_vars * sizeof(int));
  witness_t W_proj;
  FILE *PROJ = NULL;
  trackingStats trackCount;

  init_trackingStats(&trackCount); // initialize trackCount to all 0

  // verify the component is affine
  if (W->PPD.num_hom_var_gp > 0)
  { // error
    rename("witness_data_old", "witness_data");

    printf("Witness set projection is only implemented for affine varieties (variable_group)!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // compute the degree of the component and the points which correspond to it
  for (i = 0; i < num_set; i++)
    if (W->codim[codim_index].component_nums[i] == component_number)
    { // increment degree and store the path number
      degree++;
      pathNums = (int *)brealloc(pathNums, degree * sizeof(int));
      pathNums[degree - 1] = i;

      // verify multiplicity and number of deflations
      if (W->codim[codim_index].multiplicities[i] != 1 || W->codim[codim_index].deflations_needed[i] != 0)
      { // error
        rename("witness_data_old", "witness_data");

        printf("\nERROR: Component %d of dimension %d is not generically reduced!\n", component_number, dim);
        bexit(ERROR_CONFIGURATION);
      }
    }

  // verify that degree > 0
  if (degree == 0)
  { // error
    rename("witness_data_old", "witness_data");

    printf("\nERROR: Component %d of dimension %d does not exist!\n", component_number, dim);
    bexit(ERROR_CONFIGURATION);
  }

  // read in the projection 
  PROJ = fopen(projName, "r");
  if (PROJ == NULL)
  { // file does not exist
    rename("witness_data_old", "witness_data");

    printf("\n\nERROR: '%s' does not exist!!!\n\n", projName);
    bexit(ERROR_FILE_NOT_EXIST);
  }
  for (i = 0; i < num_affine_vars; i++)
  {
    rV = fscanf(PROJ, "%d", &projVars[i]);
    if (rV < 0 || projVars[i] < 0 || projVars[i] > 1)
    { // we are at EOF
      rename("witness_data_old", "witness_data");

      printf("\n\nERROR: The projection is improperly defined in '%s'!\n\n", projName);
      bexit(ERROR_INVALID_SIZE);
    }
    if (!projVars[i])
      num_fiber_vars++;
  }
  // close PROJ
  fclose(PROJ);

  // verify that this is a nontrivial projection
  if (num_fiber_vars == 0 || num_fiber_vars == num_affine_vars)
  { // error
    rename("witness_data_old", "witness_data");

    printf("\n\nERROR: '%s' describes a trivial projection!\n", projName);
    bexit(ERROR_CONFIGURATION);
  }

  // use the first witness point to compute the projection dimension and fiber dimension
  fiber_dim = computeFiberDim(W, T, codim_index, component_number, pathNums[0], projVars);
  projection_dim = dim - fiber_dim;
  printf("\n\n");

  // compute a witness set for the projection
  computeProjectionWitnessSet(&W_proj, W, T, codim_index, component_number, degree, pathNums, projVars, projection_dim, fiber_dim, &projection_degree, &fiber_degree, pathMod, &trackCount);

  // print a message about the dimension
  printf("\nDimensions\n  Projection: %d\n       Fiber: %d\n", projection_dim, fiber_dim);

  // compute the projection degree and fiber degree by looking at the points in W_proj
  printf("\nDegrees\n  Projection: %d\n       Fiber: %d\n", projection_degree, fiber_degree);

  // print the witness_data file for the projection
  numIrredDecompWitnessData(outName, &W_proj, T->MPType);
  printf("\nPrinted a witness set to '%s'.\n\n", outName);

  // print summary
  printFailureSummary(&trackCount, 0, 0, 0);

  // clear memory
  witness_clear(&W_proj, T->MPType);
  free(pathNums);
  free(projVars);

  return;
}

int computeFiberDim(witness_t *W, tracker_config_t *T, int codim_index, int component_number, int pathNum, int *projection)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute the fiber dimension for a projection           *
\***************************************************************/
{
  int i, j, k, fiberDim = 0, num_fiber_vars = 0, num_affine_vars = W->orig_variables - W->PPD.num_var_gp;
  double CN, s, l, max_CN, max_SV_ratio, SV_tol;

  // count the number of variables in the fiber
  for (i = 0; i < num_affine_vars; i++)
    if (!projection[i])
      num_fiber_vars++;

  // verify that this is a nontrivial projection
  if (num_fiber_vars == 0 || num_fiber_vars == num_affine_vars)
  { // error
    rename("witness_data_old", "witness_data");

    printf("\n\nERROR: The projection is trivial!\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (T->MPType == 0)
  { // use _d
    comp_d time;
    mat_d Jv_orig, Jv_other, Jv_temp;
    eval_struct_d e;

    // initialize
    max_CN = 1e13;
    max_SV_ratio = T->ratioTol; 
    SV_tol = MAX(T->sing_val_zero_tol, 1e-15);
    set_zero_d(time);
    init_mat_d(Jv_orig, 0, 0);
    init_mat_d(Jv_other, 0, 0);
    init_mat_d(Jv_temp, 0, 0);
    init_eval_struct_d(e, 0, 0, 0);

    // evalute at the point
    evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, W->codim[codim_index].witnessPts_d[pathNum].endPt, time, W->Prog);

    // evaluate at the other point
    evalProg_d(e.funcVals, e.parVals, e.parDer, Jv_temp, e.Jp, W->codim[codim_index].witnessPts_d[pathNum].last_approx, time, W->Prog);

    // copy over the columns corresponding to the fiber variables
    change_size_mat_d(Jv_orig, e.Jv->rows, num_fiber_vars);
    change_size_mat_d(Jv_other, e.Jv->rows, num_fiber_vars);
    Jv_orig->rows = Jv_other->rows = e.Jv->rows;
    Jv_orig->cols = Jv_other->cols = num_fiber_vars;

    k = 0;
    for (j = 0; j < num_affine_vars; j++)
      if (!projection[j])
      { // copy over to column k
        for (i = 0; i < e.Jv->rows; i++)
        {
          set_d(&Jv_orig->entry[i][k], &e.Jv->entry[i][j + W->PPD.num_var_gp]);
          set_d(&Jv_other->entry[i][k], &Jv_temp->entry[i][j + W->PPD.num_var_gp]);
        }
        // increment k
        k++;
      }

    // compute the fiber dimension == nullity
    fiberDim = corank_rrv_d(&CN, &s, &l, Jv_orig, Jv_other, 0, 0, max_CN, max_SV_ratio, SV_tol); 

    // clear
    clear_mat_d(Jv_orig);
    clear_mat_d(Jv_other);
    clear_mat_d(Jv_temp);
    clear_eval_struct_d(e);
  }
  else if (T->MPType == 1)
  { // use _mp
    comp_mp time;
    mat_mp Jv_orig, Jv_other, Jv_temp;
    eval_struct_mp e;

    // initialize
    i = prec_to_digits(T->Precision) - 4;
    max_CN = MIN(1e150, pow(10, i));
    max_SV_ratio = T->ratioTol; 
    SV_tol = MAX(T->sing_val_zero_tol, pow(10, -i - 2));
    init_mp(time); 
    set_zero_mp(time);
    init_mat_mp(Jv_orig, 0, 0);
    init_mat_mp(Jv_other, 0, 0);
    init_mat_mp(Jv_temp, 0, 0);
    init_eval_struct_mp(e, 0, 0, 0);

    // evalute at the point
    evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, W->codim[codim_index].witnessPts_mp[pathNum].endPt, time, W->Prog);

    // evaluate at the other point
    evalProg_mp(e.funcVals, e.parVals, e.parDer, Jv_temp, e.Jp, W->codim[codim_index].witnessPts_mp[pathNum].last_approx, time, W->Prog);

    // copy over the columns corresponding to the fiber variables
    change_size_mat_mp(Jv_orig, e.Jv->rows, num_fiber_vars);
    change_size_mat_mp(Jv_other, e.Jv->rows, num_fiber_vars);
    Jv_orig->rows = Jv_other->rows = e.Jv->rows;
    Jv_orig->cols = Jv_other->cols = num_fiber_vars;

    k = 0;
    for (j = 0; j < num_affine_vars; j++)
      if (!projection[j])
      { // copy over to column k
        for (i = 0; i < e.Jv->rows; i++)
        {
          set_mp(&Jv_orig->entry[i][k], &e.Jv->entry[i][j + W->PPD.num_var_gp]);
          set_mp(&Jv_other->entry[i][k], &Jv_temp->entry[i][j + W->PPD.num_var_gp]);
        }
        // increment k
        k++;
      }

    // compute the fiber dimension == nullity
    fiberDim = corank_rrv_mp(&CN, &s, &l, Jv_orig, Jv_other, 0, 0, max_CN, max_SV_ratio, SV_tol); 

    // clear
    clear_mp(time);
    clear_mat_mp(Jv_orig);
    clear_mat_mp(Jv_other);
    clear_mat_mp(Jv_temp);
    clear_eval_struct_mp(e);
  }
  else
  { // use _amp
    if (W->codim[codim_index].witnessPts_amp[pathNum].curr_prec < 64 && W->codim[codim_index].witnessPts_amp[pathNum].last_approx_prec < 64)
    { // use _d
      comp_d time;
      mat_d Jv_orig, Jv_other, Jv_temp;
      eval_struct_d e;

      // initialize
      max_CN = 1e13;
      max_SV_ratio = T->ratioTol; 
      SV_tol = MAX(T->sing_val_zero_tol, 1e-15);
      set_zero_d(time);
      init_mat_d(Jv_orig, 0, 0);
      init_mat_d(Jv_other, 0, 0);
      init_mat_d(Jv_temp, 0, 0);
      init_eval_struct_d(e, 0, 0, 0);

      // evalute at the point
      evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, W->codim[codim_index].witnessPts_amp[pathNum].endPt_d, time, W->Prog);

      // evaluate at the other point
      evalProg_d(e.funcVals, e.parVals, e.parDer, Jv_temp, e.Jp, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_d, time, W->Prog);

      // copy over the columns corresponding to the fiber variables
      change_size_mat_d(Jv_orig, e.Jv->rows, num_fiber_vars);
      change_size_mat_d(Jv_other, e.Jv->rows, num_fiber_vars);
      Jv_orig->rows = Jv_other->rows = e.Jv->rows;
      Jv_orig->cols = Jv_other->cols = num_fiber_vars;

      k = 0;
      for (j = 0; j < num_affine_vars; j++)
        if (!projection[j])
        { // copy over to column k
          for (i = 0; i < e.Jv->rows; i++)
          {
            set_d(&Jv_orig->entry[i][k], &e.Jv->entry[i][j + W->PPD.num_var_gp]);
            set_d(&Jv_other->entry[i][k], &Jv_temp->entry[i][j + W->PPD.num_var_gp]);
          }
          // increment k
          k++;
        }

      // compute the fiber dimension == nullity
      fiberDim = corank_rrv_d(&CN, &s, &l, Jv_orig, Jv_other, 0, 0, max_CN, max_SV_ratio, SV_tol); 

      // clear
      clear_mat_d(Jv_orig);
      clear_mat_d(Jv_other);
      clear_mat_d(Jv_temp);
      clear_eval_struct_d(e);
    }
    else 
    { // use _mp
      int prec = MAX(W->codim[codim_index].witnessPts_amp[pathNum].curr_prec, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_prec);
      comp_mp time;
      point_mp pt;
      mat_mp Jv_orig, Jv_other, Jv_temp;
      eval_struct_mp e;

      // set precision correctly
      initMP(prec);
      W->Prog->precision = prec;

      // initialize
      i = prec_to_digits(prec) - 4;
      max_CN = MIN(1e150, pow(10, i));
      max_SV_ratio = T->ratioTol; 
      SV_tol = MAX(T->sing_val_zero_tol, pow(10, -i - 2));
      init_mp(time); 
      set_zero_mp(time);
      init_point_mp(pt, 0);
      init_mat_mp(Jv_orig, 0, 0);
      init_mat_mp(Jv_other, 0, 0);
      init_mat_mp(Jv_temp, 0, 0);
      init_eval_struct_mp(e, 0, 0, 0);

      // setup pt
      if (W->codim[codim_index].witnessPts_amp[pathNum].curr_prec < 64)
      { // copy _d to _mp
        point_d_to_mp(pt, W->codim[codim_index].witnessPts_amp[pathNum].endPt_d);
      }
      else
      { // copy _mp to _mp
        point_cp_mp(pt, W->codim[codim_index].witnessPts_amp[pathNum].endPt_mp);
      }

      // evalute at the point
      evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, pt, time, W->Prog);

      // setup pt
      if (W->codim[codim_index].witnessPts_amp[pathNum].last_approx_prec < 64)
      { // copy _d to _mp
        point_d_to_mp(pt, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_d);
      }
      else
      { // copy _mp to _mp
        point_cp_mp(pt, W->codim[codim_index].witnessPts_amp[pathNum].last_approx_mp);
      }

      // evaluate at the other point
      evalProg_mp(e.funcVals, e.parVals, e.parDer, Jv_temp, e.Jp, pt, time, W->Prog);

      // copy over the columns corresponding to the fiber variables
      change_size_mat_mp(Jv_orig, e.Jv->rows, num_fiber_vars);
      change_size_mat_mp(Jv_other, e.Jv->rows, num_fiber_vars);
      Jv_orig->rows = Jv_other->rows = e.Jv->rows;
      Jv_orig->cols = Jv_other->cols = num_fiber_vars;

      k = 0;
      for (j = 0; j < num_affine_vars; j++)
        if (!projection[j])
        { // copy over to column k
          for (i = 0; i < e.Jv->rows; i++)
          {
            set_mp(&Jv_orig->entry[i][k], &e.Jv->entry[i][j + W->PPD.num_var_gp]);
            set_mp(&Jv_other->entry[i][k], &Jv_temp->entry[i][j + W->PPD.num_var_gp]);
          }
          // increment k
          k++;
        }

      // compute the fiber dimension == nullity
      fiberDim = corank_rrv_mp(&CN, &s, &l, Jv_orig, Jv_other, 0, 0, max_CN, max_SV_ratio, SV_tol); 

      // clear
      clear_mp(time);
      clear_point_mp(pt);
      clear_mat_mp(Jv_orig);
      clear_mat_mp(Jv_other);
      clear_mat_mp(Jv_temp);
      clear_eval_struct_mp(e);
    }
  }

  return fiberDim;
} 

void setupProjectionWitnessSet(witness_t *W_proj, witness_t *W, tracker_config_t *T, int codim_index, int *projection, int proj_dim, int fiber_dim)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the basic parts of W_proj                        *
\***************************************************************/
{
  int i, j, num_affine_vars = W->orig_variables - W->PPD.num_var_gp;

  // setup the basic parts of W_proj
  // copy Prog
  W_proj->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(W_proj->Prog, W->Prog);
  // copy PPD
  cp_preproc_data(&W_proj->PPD, &W->PPD);
  // copy other data
  W_proj->num_funcs = W->num_funcs;
  W_proj->curr_precision = W->curr_precision;
  W_proj->system_rank = W->system_rank;
  W_proj->orig_variables = W->orig_variables;
  W_proj->new_variables = W->new_variables;
  W_proj->orig_degrees = (int *)bmalloc(W_proj->num_funcs * sizeof(int));
  W_proj->new_degrees = (int *)bmalloc(W_proj->num_funcs * sizeof(int));
  W_proj->P = (int *)bmalloc(W_proj->num_funcs * sizeof(int));
  for (i = 0; i < W_proj->num_funcs; i++)
  {
    W_proj->orig_degrees[i] = W->orig_degrees[i];
    W_proj->new_degrees[i] = W->new_degrees[i];
    W_proj->P[i] = W->P[i];
  }

  // setup C
  if (W_proj->new_variables != W_proj->orig_variables)
  { // setup C
    if (T->MPType == 0)
    { // copy C_d
      init_mat_d(W_proj->C_d, W_proj->orig_variables, W_proj->new_variables);
      mat_cp_d(W_proj->C_d, W->C_d);
    }
    else if (T->MPType == 1)
    { // copy C_mp
      init_mat_mp(W_proj->C_mp, W_proj->orig_variables, W_proj->new_variables);
      mat_cp_mp(W_proj->C_mp, W->C_mp);
    }
    else
    { // copy C_d,_mp,_rat
      init_mat_d(W_proj->C_d, W_proj->orig_variables, W_proj->new_variables);
      init_mat_mp2(W_proj->C_mp, W_proj->orig_variables, W_proj->new_variables, W_proj->curr_precision);
      init_mat_rat(W_proj->C_rat, W_proj->orig_variables, W_proj->new_variables);
      mat_cp_d(W_proj->C_d, W->C_d);
      mat_cp_mp(W_proj->C_mp, W->C_mp);
      mat_cp_rat(W_proj->C_rat, W->C_rat, W_proj->orig_variables, W_proj->new_variables);
    }
  }

  // setup gamma
  if (T->MPType == 0)
  { // copy gamma_d
    set_d(W_proj->gamma_d, W->gamma_d);
  }
  else if (T->MPType == 1)
  { // copy gamma_mp
    init_mp(W_proj->gamma_mp);
    set_mp(W_proj->gamma_mp, W->gamma_mp);
  }
  else
  { // copy gamma_d,_mp,_rat
    init_mp2(W_proj->gamma_mp, W_proj->curr_precision);
    init_rat(W_proj->gamma_rat);
    set_d(W_proj->gamma_d, W->gamma_d);
    set_mp(W_proj->gamma_mp, W->gamma_mp);
    set_rat(W_proj->gamma_rat, W->gamma_rat);
  }
  W_proj->num_codim = 1;
  W_proj->curr_codim_index = 0;
  W_proj->targetSliceInit = 0;

  // setup codim
  W_proj->codim = (witnessCodim_t *)bmalloc(1 * sizeof(witnessCodim_t));
  W_proj->codim[0].codim = W->codim[codim_index].codim;
  // setup A
  W_proj->codim[0].A_rows = W->codim[codim_index].A_rows;
  W_proj->codim[0].A_cols = W->codim[codim_index].A_cols;
  if (T->MPType == 0)
  { // setup A_d
    init_mat_d(W_proj->codim[0].A_d, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols);
    mat_cp_d(W_proj->codim[0].A_d, W->codim[codim_index].A_d);
  }
  else if (T->MPType == 1)
  { // setup A_mp
    init_mat_mp(W_proj->codim[0].A_mp, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols);
    mat_cp_mp(W_proj->codim[0].A_mp, W->codim[codim_index].A_mp);
  }
  else
  { // setup A_d,_mp,_rat
    init_mat_d(W_proj->codim[0].A_d, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols);
    init_mat_mp2(W_proj->codim[0].A_mp, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols, W_proj->curr_precision);
    init_mat_rat(W_proj->codim[0].A_rat, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols);
    mat_cp_d(W_proj->codim[0].A_d, W->codim[codim_index].A_d);
    mat_cp_mp(W_proj->codim[0].A_mp, W->codim[codim_index].A_mp);
    mat_cp_rat(W_proj->codim[0].A_rat, W->codim[codim_index].A_rat, W_proj->codim[0].A_rows, W_proj->codim[0].A_cols);
  }
  // setup W
  W_proj->codim[0].W = (int **)bmalloc(W_proj->codim[0].A_rows * sizeof(int *));
  for (i = 0; i < W_proj->codim[0].A_rows; i++)
  {
    W_proj->codim[0].W[i] = (int *)bmalloc(W_proj->codim[0].A_cols * sizeof(int));
    for (j = 0; j < W_proj->codim[0].A_cols; j++)
      W_proj->codim[0].W[i][j] = W->codim[codim_index].W[i][j];
  }
  // setup H
  if (T->MPType == 0)
  { // setup H_d
    init_vec_d(W_proj->codim[0].H_d, W->codim[codim_index].H_d->size);
    vec_cp_d(W_proj->codim[0].H_d, W->codim[codim_index].H_d);
  }
  else if (T->MPType == 1)
  { // setup H_mp
    init_vec_mp(W_proj->codim[0].H_mp, W->codim[codim_index].H_mp->size);
    vec_cp_mp(W_proj->codim[0].H_mp, W->codim[codim_index].H_mp);
  }
  else
  { // setup H_rat
    init_vec_d(W_proj->codim[0].H_d, W->codim[codim_index].H_d->size);
    init_vec_mp2(W_proj->codim[0].H_mp, W->codim[codim_index].H_mp->size, W_proj->curr_precision);
    init_vec_rat(W_proj->codim[0].H_rat, W->codim[codim_index].H_mp->size);
    vec_cp_d(W_proj->codim[0].H_d, W->codim[codim_index].H_d);
    vec_cp_mp(W_proj->codim[0].H_mp, W->codim[codim_index].H_mp);
    vec_cp_rat(W_proj->codim[0].H_rat, W->codim[codim_index].H_rat, W->codim[codim_index].H_mp->size);
  } 
  // setup homVarConst
  if (T->MPType == 0)
  { // copy homVarConst_d
    set_d(W_proj->codim[0].homVarConst_d, W->codim[codim_index].homVarConst_d);
  }
  else if (T->MPType == 1)
  { // copy homVarConst_mp
    init_mp(W_proj->codim[0].homVarConst_mp);
    set_mp(W_proj->codim[0].homVarConst_mp, W->codim[codim_index].homVarConst_mp);
  }
  else
  { // copy homVarConst_d,_mp,_rat
    init_mp2(W_proj->codim[0].homVarConst_mp, W_proj->curr_precision);
    init_rat(W_proj->codim[0].homVarConst_rat);
    set_d(W_proj->codim[0].homVarConst_d, W->codim[codim_index].homVarConst_d);
    set_mp(W_proj->codim[0].homVarConst_mp, W->codim[codim_index].homVarConst_mp);
    set_rat(W_proj->codim[0].homVarConst_rat, W->codim[codim_index].homVarConst_rat);
  }
  // setup p
  if (T->MPType == 0)
  { // setup p_d
    init_vec_d(W_proj->codim[0].p_d, W->codim[codim_index].p_d->size);
    vec_cp_d(W_proj->codim[0].p_d, W->codim[codim_index].p_d);
  }
  else if (T->MPType == 1)
  { // setup p_mp
    init_vec_mp(W_proj->codim[0].p_mp, W->codim[codim_index].p_mp->size);
    vec_cp_mp(W_proj->codim[0].p_mp, W->codim[codim_index].p_mp);
  }
  else
  { // setup p_rat
    init_vec_d(W_proj->codim[0].p_d, W->codim[codim_index].p_d->size);
    init_vec_mp2(W_proj->codim[0].p_mp, W->codim[codim_index].p_mp->size, W_proj->curr_precision);
    init_vec_rat(W_proj->codim[0].p_rat, W->codim[codim_index].p_mp->size);
    vec_cp_d(W_proj->codim[0].p_d, W->codim[codim_index].p_d);
    vec_cp_mp(W_proj->codim[0].p_mp, W->codim[codim_index].p_mp);
    vec_cp_rat(W_proj->codim[0].p_rat, W->codim[codim_index].p_rat, W->codim[codim_index].p_mp->size);
  }

  // setup B - based on projection & dimensions
  if (T->MPType == 0)
  { // setup B_d
    init_mat_d(W_proj->codim[0].B_d, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    mat_cp_d(W_proj->codim[0].B_d, W->codim[codim_index].B_d);
    for (j = 0; j < num_affine_vars; j++)
      if (projection[j])
      { // clear out fiber linears
        for (i = proj_dim; i < proj_dim + fiber_dim; i++)
        {
          set_zero_d(&W_proj->codim[0].B_d->entry[i][j + W_proj->PPD.num_var_gp]);
        }
      }
      else
      { // clear out proj linears
        for (i = 0; i < proj_dim; i++)
        {
          set_zero_d(&W_proj->codim[0].B_d->entry[i][j + W_proj->PPD.num_var_gp]);
        }
      }
  }
  else if (T->MPType == 1)
  { // setup B_mp
    init_mat_mp(W_proj->codim[0].B_mp, W->codim[codim_index].B_mp->rows, W->codim[codim_index].B_mp->cols);
    mat_cp_mp(W_proj->codim[0].B_mp, W->codim[codim_index].B_mp);
    for (j = 0; j < num_affine_vars; j++)
      if (projection[j])
      { // clear out fiber linears
        for (i = proj_dim; i < proj_dim + fiber_dim; i++)
        {
          set_zero_mp(&W_proj->codim[0].B_mp->entry[i][j + W_proj->PPD.num_var_gp]);
        }
      }
      else
      { // clear out proj linears
        for (i = 0; i < proj_dim; i++)
        {
          set_zero_mp(&W_proj->codim[0].B_mp->entry[i][j + W_proj->PPD.num_var_gp]);
        }
      }
  }
  else
  { // setup B_d,_mp,_rat
    init_mat_d(W_proj->codim[0].B_d, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    init_mat_mp2(W_proj->codim[0].B_mp, W->codim[codim_index].B_mp->rows, W->codim[codim_index].B_mp->cols, W_proj->curr_precision);
    init_mat_rat(W_proj->codim[0].B_rat, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    mat_cp_d(W_proj->codim[0].B_d, W->codim[codim_index].B_d);
    mat_cp_mp(W_proj->codim[0].B_mp, W->codim[codim_index].B_mp);
    mat_cp_rat(W_proj->codim[0].B_rat, W->codim[codim_index].B_rat, W->codim[codim_index].B_d->rows, W->codim[codim_index].B_d->cols);
    for (j = 0; j < num_affine_vars; j++)
      if (projection[j])
      { // clear out fiber linears
        for (i = proj_dim; i < proj_dim + fiber_dim; i++)
        {
          set_zero_d(&W_proj->codim[0].B_d->entry[i][j + W_proj->PPD.num_var_gp]);
          set_zero_mp(&W_proj->codim[0].B_mp->entry[i][j + W_proj->PPD.num_var_gp]);
          set_zero_rat(W_proj->codim[0].B_rat[i][j + W_proj->PPD.num_var_gp]);
        }
      }
      else
      { // clear out proj linears
        for (i = 0; i < proj_dim; i++)
        {
          set_zero_d(&W_proj->codim[0].B_d->entry[i][j + W_proj->PPD.num_var_gp]);
          set_zero_mp(&W_proj->codim[0].B_mp->entry[i][j + W_proj->PPD.num_var_gp]);
          set_zero_rat(W_proj->codim[0].B_rat[i][j + W_proj->PPD.num_var_gp]);
        }
      }
  }

  // initialize other things
  W_proj->codim[0].num_set = W_proj->codim[0].num_nonsing = W_proj->codim[0].num_sing = 0;
  W_proj->codim[0].witnessPts_d = NULL;
  W_proj->codim[0].witnessPts_mp = NULL;
  W_proj->codim[0].witnessPts_amp = NULL;
  W_proj->codim[0].witnessPt_types = W_proj->codim[0].component_nums = W_proj->codim[0].multiplicities = W_proj->codim[0].deflations_needed = NULL;
  W_proj->codim[0].num_components = 0;

  return;
}

int projection_sortEndpoint(point_data_d *PD_d, point_data_mp *PD_mp, int Pt_prec, int retVal_in, preproc_data *PPD, trackingStats *trackCount, FILE *FAIL, int pathNum, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - good, 1 - bad                              *
* NOTES: determine the correct retVal & update trackCount       *
\***************************************************************/
{
  int rV = 0;
  double norm = 0;

  if (Pt_prec < 64)
  { // use _d
    point_d dehom_d;
    init_point_d(dehom_d, 0);

    // compute dehom_d
    getDehomPoint_d(dehom_d, PD_d->point, PD_d->point->size, PPD);
    norm = infNormVec_d(dehom_d);

    // check to see if finite
    if (T->finiteThreshold < norm)
      retVal_in = retVal_going_to_infinity;

    if (retVal_in != 0)
    { // set rV and update the number of failures
      rV = 1;
      trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_d(FAIL, PD_d, dehom_d, pathNum, retVal_in, 1, 0, trackCount, T);
    }
    else
    { // set rV and update the number of successes
      rV = 0;
      trackCount->successes++;
    }

    clear_point_d(dehom_d);
  }
  else
  { // use _mp
    comp_d tempDouble;
    point_mp dehom_mp;
    init_point_mp2(dehom_mp, 0, Pt_prec);
  
    mp_to_d(tempDouble, PD_mp->time);

    // compute dehom_mp
    getDehomPoint_mp(dehom_mp, PD_mp->point, PD_mp->point->size, PPD);
    norm = infNormVec_mp(dehom_mp);

    // check to see if finite
    if (T->finiteThreshold < norm)
      retVal_in = retVal_going_to_infinity;

    if (retVal_in != 0)
    { // set rV and update the number of failures
      rV = 1;
      trackCount->failures++;

      // print the path number, error message, time and point to FAIL
      printFailureMsg_mp(FAIL, PD_mp, dehom_mp, pathNum, retVal_in, 1, 0, trackCount, T);
    }
    else
    { // set rV and update the number of successes
      rV = 0;
      trackCount->successes++;
    }
    
    clear_point_mp(dehom_mp);
  }

  return rV;
}

void completeProjectionWitnessSet(witness_t *W, tracker_config_t *T, int *projection, int *proj_deg, int *fiber_deg)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: completes a witness set for the projection             *
\***************************************************************/
{
  int i, j, count, numVars = W->orig_variables, numAffineVars = W->orig_variables - W->PPD.num_var_gp, numSet = 0, projSize = 0;
  witness_t W_proj;

  // compute multiplicities of the points and remove extra
  multiplicity_witness(W, 0, T->MPType, T->final_tol_times_mult);
  numSet = W->codim[0].num_set;

  // compute projection size
  for (i = 0; i < numAffineVars; i++)
    projSize += projection[i];

  // setup W_proj
  W_proj.codim = (witnessCodim_t *)bmalloc(1 * sizeof(witnessCodim_t));
  W_proj.codim[0].multiplicities = (int *)bmalloc(numSet * sizeof(int));
  W_proj.codim[0].num_set = W->codim[0].num_set;
  if (T->MPType == 0)
  { // allocate and setup witnessPts_d
    W_proj.codim[0].witnessPts_d = (endpoint_data_d *)bmalloc(numSet * sizeof(endpoint_data_d));
    for (i = 0; i < numSet; i++)
    {
      init_endpoint_data_d(&W_proj.codim[0].witnessPts_d[i]);
      endpoint_data_cp_d(&W_proj.codim[0].witnessPts_d[i], &W->codim[0].witnessPts_d[i]);
      // convert to original coordinates
      getDehomPoint_d(W_proj.codim[0].witnessPts_d[i].endPt, W->codim[0].witnessPts_d[i].endPt, numVars, &W->PPD);
      // perform the projection
      count = 0;
      for (j = 0; j < numVars; j++)
        if (projection[j])
        {
          set_d(&W_proj.codim[0].witnessPts_d[i].endPt->coord[count], &W_proj.codim[0].witnessPts_d[i].endPt->coord[j]);
          count++;
        }
      W_proj.codim[0].witnessPts_d[i].endPt->size = projSize;
    }
  }
  else if (T->MPType == 1)
  { // allocate and setup witnessPts_mp
    W_proj.codim[0].witnessPts_mp = (endpoint_data_mp *)bmalloc(numSet * sizeof(endpoint_data_mp));
    for (i = 0; i < numSet; i++)
    {
      init_endpoint_data_mp(&W_proj.codim[0].witnessPts_mp[i]);
      endpoint_data_cp_mp(&W_proj.codim[0].witnessPts_mp[i], &W->codim[0].witnessPts_mp[i]);
      // convert to original coordinates
      getDehomPoint_mp(W_proj.codim[0].witnessPts_mp[i].endPt, W->codim[0].witnessPts_mp[i].endPt, numVars, &W->PPD);
      // perform the projection
      count = 0;
      for (j = 0; j < numVars; j++)
        if (projection[j])
        {
          set_mp(&W_proj.codim[0].witnessPts_mp[i].endPt->coord[count], &W_proj.codim[0].witnessPts_mp[i].endPt->coord[j]);
          count++;
        }
      W_proj.codim[0].witnessPts_mp[i].endPt->size = projSize;
    }
  }
  else
  { // allocate and setup witnessPts_amp 
    W_proj.codim[0].witnessPts_amp = (endpoint_data_amp *)bmalloc(numSet * sizeof(endpoint_data_amp));
    for (i = 0; i < numSet; i++)
    {
      init_endpoint_data_amp(&W_proj.codim[0].witnessPts_amp[i], W->codim[0].witnessPts_amp[i].curr_prec, W->codim[0].witnessPts_amp[i].last_approx_prec);
      endpoint_data_cp_amp(&W_proj.codim[0].witnessPts_amp[i], &W->codim[0].witnessPts_amp[i]);
      // convert to original coordinates
      if (W_proj.codim[0].witnessPts_amp[i].curr_prec < 64)
        getDehomPoint_d(W_proj.codim[0].witnessPts_amp[i].endPt_d, W->codim[0].witnessPts_amp[i].endPt_d, numVars, &W->PPD);
      else
        getDehomPoint_mp(W_proj.codim[0].witnessPts_amp[i].endPt_mp, W->codim[0].witnessPts_amp[i].endPt_mp, numVars, &W->PPD);

      // perform the projection
      count = 0;
      if (W_proj.codim[0].witnessPts_amp[i].curr_prec < 64)
      { // endPt_d
        for (j = 0; j < numVars; j++)
          if (projection[j])
          {
            set_d(&W_proj.codim[0].witnessPts_amp[i].endPt_d->coord[count], &W_proj.codim[0].witnessPts_amp[i].endPt_d->coord[j]);
            count++;
          }
        W_proj.codim[0].witnessPts_amp[i].endPt_d->size = projSize;
      }
      else
      { // endPt_mp
        for (j = 0; j < numVars; j++)
          if (projection[j])
          {
            set_mp(&W_proj.codim[0].witnessPts_amp[i].endPt_mp->coord[count], &W_proj.codim[0].witnessPts_amp[i].endPt_mp->coord[j]);
            count++;
          }
        W_proj.codim[0].witnessPts_amp[i].endPt_mp->size = projSize;
      }
    }
  } 

  // compute the number of points over each point
  sort_endpoint_data(W_proj.codim[0].multiplicities, W_proj.codim[0].witnessPts_d, W_proj.codim[0].witnessPts_mp, W_proj.codim[0].witnessPts_amp, numSet, T->MPType, T->final_tol_times_mult);

  // count the number with positive multiplicity and verify that all have the same multiplicity
  *proj_deg = 0;
  *fiber_deg = 0;
  for (i = 0; i < numSet; i++)
    if (W_proj.codim[0].multiplicities[i] > 0)
    { // increment proj_deg
      (*proj_deg)++;
      // setup fiber_deg or check the multiplicities
      if (*fiber_deg == 0)
      { // set fiber_deg
        *fiber_deg = W_proj.codim[0].multiplicities[i];
      }
      else if (*fiber_deg != W_proj.codim[0].multiplicities[i])
      { // error!
        printf("ERROR: The number of points in the fibers are not constant!\n");
        bexit(ERROR_CONFIGURATION);
      }
    }

  // setup other structures
  W->codim[0].component_nums = (int *)bmalloc(W->codim[0].num_set * sizeof(int));
  W->codim[0].deflations_needed = (int *)bmalloc(W->codim[0].num_set * sizeof(int));

  // deflations needed to 0 & component number to 0
  for (i = 0; i < W->codim[0].num_set; i++)
    W->codim[0].component_nums[i] = W->codim[0].deflations_needed[i] = 0;

  // clear W_proj
  if (T->MPType == 0)
  {
    for (i = 0; i < W_proj.codim[0].num_set; i++)
      clear_endpoint_data_d(&W_proj.codim[0].witnessPts_d[i]);
    free(W_proj.codim[0].witnessPts_d);
  }
  else if (T->MPType == 1)
  {
    for (i = 0; i < W_proj.codim[0].num_set; i++)
      clear_endpoint_data_mp(&W_proj.codim[0].witnessPts_mp[i]);
    free(W_proj.codim[0].witnessPts_mp);
  }
  else 
  {
    for (i = 0; i < W_proj.codim[0].num_set; i++)
      clear_endpoint_data_amp(&W_proj.codim[0].witnessPts_amp[i]);
    free(W_proj.codim[0].witnessPts_amp);
  }
  free(W_proj.codim[0].multiplicities);
  free(W_proj.codim);

  return;
}

void computeProjectionWitnessSet(witness_t *W_proj, witness_t *W, tracker_config_t *T, int codim_index, int component_number, int degree, int *pathNums, int *projection, int proj_dim, int fiber_dim, int *proj_deg, int *fiber_deg, int pathMod, trackingStats *trackCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: compute a witness set for the projection               *
\***************************************************************/
{
  int i, count;
  general_slice_moving_t M;
  endgame_data_t endGame;
  FILE *OUT = fopen("output_projection", "w"), *MIDOUT = fopen("midout_projection", "w"), *FAIL = fopen("failed_paths", "w");

  // initialize
  init_endgame_data(&endGame, T->Precision);
  trackCount->numPoints = degree;

  // setup the basic parst of W_proj
  setupProjectionWitnessSet(W_proj, W, T, codim_index, projection, proj_dim, fiber_dim);

  // setup M using W & W_proj
  initialize_setup_general_slice(&M, W, codim_index, W_proj, 0, T->MPType);

  // move the points to these slices
  if (T->MPType == 0)
  { // allocate endPts_d
    endpoint_data_d *endPts_d = (endpoint_data_d *)bmalloc(degree * sizeof(endpoint_data_d));
    // loop over to compute endPts
    for (i = 0; i < degree; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, degree);

      // initialize endPts_d
      init_endpoint_data_d(&endPts_d[i]);

      // track the path
      general_slice_moving_track(&endPts_d[i], NULL, NULL, &endGame, &M, W->codim[codim_index].witnessPts_d[pathNums[i]].endPt, NULL, 52, i, 1, T, OUT, MIDOUT);

      // sort the endpoint
      endPts_d[i].retVal = projection_sortEndpoint(&endGame.PD_d, &endGame.PD_mp, endGame.prec, endPts_d[i].retVal, &W_proj->PPD, trackCount, FAIL, i, T);
    } 

    // setup structures in W_proj
    W_proj->codim[0].num_set = trackCount->successes;
    W_proj->codim[0].witnessPts_d = (endpoint_data_d *)bmalloc(trackCount->successes * sizeof(endpoint_data_d));
    W_proj->codim[0].witnessPt_types = (int *)bmalloc(trackCount->successes * sizeof(int));

    count = 0;
    for (i = 0; i < degree; i++)
    { // determine if we need to copy over
      if (!endPts_d[i].retVal)
      { // copy
        init_endpoint_data_d(&W_proj->codim[0].witnessPts_d[count]);
        endpoint_data_cp_d(&W_proj->codim[0].witnessPts_d[count], &endPts_d[i]);
        // determine is singular or nonsingular
        if (endPts_d[i].corank > 0)
        { // singular
          W_proj->codim[0].num_sing++;
          W_proj->codim[0].witnessPt_types[count] = SINGULAR;
        }
        else
        { // nonsingular
          W_proj->codim[0].num_nonsing++;
          W_proj->codim[0].witnessPt_types[count] = NON_SINGULAR;
        }

        // increment count
        count++;
      }
      // clear endPts_d
      clear_endpoint_data_d(&endPts_d[i]);
    }

    // free endPts_d
    free(endPts_d);
  }
  else if (T->MPType == 1)
  { // allocate endPts_mp
    endpoint_data_mp *endPts_mp = (endpoint_data_mp *)bmalloc(degree * sizeof(endpoint_data_mp));
    // loop over to compute endPts
    for (i = 0; i < degree; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, degree);

      // initialize endPts_mp
      init_endpoint_data_mp(&endPts_mp[i]);

      // track the path
      general_slice_moving_track(NULL, &endPts_mp[i], NULL, &endGame, &M, NULL, W->codim[codim_index].witnessPts_mp[pathNums[i]].endPt, T->Precision, i, 1, T, OUT, MIDOUT);

      // sort the endpoint
      endPts_mp[i].retVal = projection_sortEndpoint(&endGame.PD_d, &endGame.PD_mp, endGame.prec, endPts_mp[i].retVal, &W_proj->PPD, trackCount, FAIL, i, T);
    }

    // setup structures in W_proj
    W_proj->codim[0].num_set = trackCount->successes;
    W_proj->codim[0].witnessPts_mp = (endpoint_data_mp *)bmalloc(trackCount->successes * sizeof(endpoint_data_mp));
    W_proj->codim[0].witnessPt_types = (int *)bmalloc(trackCount->successes * sizeof(int));

    count = 0;
    for (i = 0; i < degree; i++)
    { // determine if we need to copy over
      if (!endPts_mp[i].retVal)
      { // copy
        init_endpoint_data_mp(&W_proj->codim[0].witnessPts_mp[count]);
        endpoint_data_cp_mp(&W_proj->codim[0].witnessPts_mp[count], &endPts_mp[i]);
        // determine is singular or nonsingular
        if (endPts_mp[i].corank > 0)
        { // singular
          W_proj->codim[0].num_sing++;
          W_proj->codim[0].witnessPt_types[count] = SINGULAR;
        }
        else
        { // nonsingular
          W_proj->codim[0].num_nonsing++;
          W_proj->codim[0].witnessPt_types[count] = NON_SINGULAR;
        }

        // increment count
        count++;
      }
      // clear endPts_mp
      clear_endpoint_data_mp(&endPts_mp[i]);
    }

    // free endPts_mp
    free(endPts_mp);
  }
  else
  { // allocate endPts_amp
    endpoint_data_amp *endPts_amp = (endpoint_data_amp *)bmalloc(degree * sizeof(endpoint_data_amp));
    // loop over to compute endPts
    for (i = 0; i < degree; i++)
    { // print the path number if needed
      if (pathMod > 0 && !(i % pathMod))
        printf("Tracking path %d of %d\n", i, degree);

      // initialize endPts_amp
      init_endpoint_data_amp(&endPts_amp[i], T->Precision, T->Precision);

      // track the path
      general_slice_moving_track(NULL, NULL, &endPts_amp[i], &endGame, &M, W->codim[codim_index].witnessPts_amp[pathNums[i]].endPt_d, W->codim[codim_index].witnessPts_amp[pathNums[i]].endPt_mp, W->codim[codim_index].witnessPts_amp[pathNums[i]].curr_prec, i, 1, T, OUT, MIDOUT);

      // sort the endpoint
      endPts_amp[i].retVal = projection_sortEndpoint(&endGame.PD_d, &endGame.PD_mp, endGame.prec, endPts_amp[i].retVal, &W_proj->PPD, trackCount, FAIL, i, T);
    }

    // setup structures in W_proj
    W_proj->codim[0].num_set = trackCount->successes;
    W_proj->codim[0].witnessPts_amp = (endpoint_data_amp *)bmalloc(trackCount->successes * sizeof(endpoint_data_amp));
    W_proj->codim[0].witnessPt_types = (int *)bmalloc(trackCount->successes * sizeof(int));

    count = 0;
    for (i = 0; i < degree; i++)
    { // determine if we need to copy over
      if (!endPts_amp[i].retVal)
      { // copy
        init_endpoint_data_amp(&W_proj->codim[0].witnessPts_amp[count], endPts_amp[i].curr_prec, endPts_amp[i].last_approx_prec);
        endpoint_data_cp_amp(&W_proj->codim[0].witnessPts_amp[count], &endPts_amp[i]);
        // determine is singular or nonsingular
        if (endPts_amp[i].corank > 0)
        { // singular
          W_proj->codim[0].num_sing++;
          W_proj->codim[0].witnessPt_types[count] = SINGULAR;
        }
        else
        { // nonsingular
          W_proj->codim[0].num_nonsing++;
          W_proj->codim[0].witnessPt_types[count] = NON_SINGULAR;
        }

        // increment count
        count++;
      }
      // clear endPts_mp
      clear_endpoint_data_amp(&endPts_amp[i]);
    }

    // clear endPts_amp
    free(endPts_amp);
  }

  // complete the witness set computation
  completeProjectionWitnessSet(W_proj, T, projection, proj_deg, fiber_deg);

  // close files
  fclose(OUT);
  fclose(MIDOUT);
  remove("midout_projection");

  // clear endGame
  clear_endgame_data(&endGame);

  // clear M
  clear_general_slice(&M, T->MPType);

  return;
}


