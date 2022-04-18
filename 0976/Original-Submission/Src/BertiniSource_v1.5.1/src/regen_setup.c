// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regeneration.h"
#include "cascade.h"

#define ZERO_HOM_COORD_D 1e-12
#define ZERO_HOM_COORD_MP 1e-12

void setupRegenRandom_zero_dim(regen_t *regen, tracker_config_t *T, char *degreeName, preproc_data *PPD, square_system_eval_data_d *SSED_d, patch_eval_data_d *patch_d, square_system_eval_data_mp *SSED_mp, patch_eval_data_mp *patch_mp, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow);

void setup_regen_from_zero_dim_seq(int max, regen_t *regen, char *startName, int startLevel, double intrinsicCutoffMultiplier, char *depthName, char *degreeName, tracker_config_t *T, basic_eval_data_d *BED_d, basic_eval_data_mp *BED_mp, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup regen for doing basic zero dimensional tracking  *
* using standard sequential operation                           *
\***************************************************************/
{
  int i, intrinsicCutoff;
  square_system_eval_data_d *square_copy_d = (square_system_eval_data_d *)bmalloc((max - 1) * sizeof(square_system_eval_data_d));
  square_system_eval_data_mp *square_copy_mp = (square_system_eval_data_mp *)bmalloc((max - 1) * sizeof(square_system_eval_data_mp));
  char *str = NULL;
  size_t size;

  // setup the random numbers in regen
  if (T->MPType == 0)
  { // setup using only BED_d
    setupRegenRandom_zero_dim(&regen[0], T, degreeName, &BED_d->preProcData, &BED_d->squareSystem, &BED_d->patch, NULL, NULL, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
  }
  else if (T->MPType == 1)
  { // setup using only BED_mp
    setupRegenRandom_zero_dim(&regen[0], T, degreeName, &BED_mp->preProcData, NULL, NULL, &BED_mp->squareSystem, &BED_mp->patch, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
  }
  else
  { // setup using BED_d & BED_mp
    setupRegenRandom_zero_dim(&regen[0], T, degreeName, &BED_d->preProcData, &BED_d->squareSystem, &BED_d->patch, &BED_mp->squareSystem, &BED_mp->patch, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
  }

  // setup the name of the start file
  size = 1 + snprintf(NULL, 0, "%s_%d", startName, startLevel);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "%s_%d", startName, startLevel);

  // determine if we are starting from scratch or restarting from a previous run
  if (startLevel == 0)
  { // starting from scratch

    // setup the intrinsic cutoff
    intrinsicCutoff = floor(intrinsicCutoffMultiplier * regen[0].num_variables);

    // allocate space for the levels needed for the regeneration
    setupRegenLevels(&regen[0], T, depthName, intrinsicCutoff);

    // setup the first level
    setupFirstLevel(&regen[0], T, str);
  }
  else
  { // we are restarting a previous run

    // read in the data from 'str' and allocate space for the levels and setup the starting level
    setupRegenRestart(&regen[0], T, str, startLevel);
  }

  // setup the other copies
  if (max > 1)
  { // setup square_d &/or square_mp
    for (i = 1; i < max; i++)
    {
      if (T->MPType == 0 || T->MPType == 2)
      { // setup square_d
        cp_square_system_d(&square_copy_d[i-1], &BED_d->squareSystem);
      }

      if (T->MPType == 1)
      { // setup square_mp using fixed precision
        cp_square_system_mp(&square_copy_mp[i-1], &BED_mp->squareSystem, 1, NULL);
      }
      else if (T->MPType == 2)
      { // setup square_mp using AMP
        cp_square_system_mp(&square_copy_mp[i-1], &BED_mp->squareSystem, 0, square_copy_d[i-1].Prog);
      }

      // setup the ith copy
      copyRegenRandom(&regen[i], &regen[0], T->MPType, &square_copy_d[i-1], &square_copy_mp[i-1]);
    }

    // setup the starting level for the copies
    copyRegenLevelStructures(&regen[0], max - 1, &regen[1], T->MPType, startLevel);
  }

  // NULL out the pointers since they are now being pointed to inside of regen
  square_copy_d = NULL;
  square_copy_mp = NULL;

  free(str);

  return;
}

void regen_setup_coeff(int MPType, int AMP_max_prec, regen_t *regen, int num_var_gps, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the coefficients for regen                       *
\***************************************************************/
{
  int i, j, k, l, beg, numFuncs = regen->num_funcs, numVars = regen->num_variables;
  int **varUsed = (int **)bmalloc(numFuncs * sizeof(int *));
  for (i = 0; i < numFuncs; i++)
    varUsed[i] = (int *)bmalloc(numVars * sizeof(int));

  if (MPType == 0)
  { // setup mainCoeff_d & coeff_d
    regen->mainCoeff_rat = NULL;
    regen->coeff_mp = NULL;
    regen->coeff_rat = NULL;

    // compute the Jacobian matrix at a random point
    comp_d time_zero;
    point_d rand_point;
    mat_d Jv;
    eval_struct_d e;

    init_point_d(rand_point, 0);
    init_mat_d(Jv, 0, 0);
    init_eval_struct_d(e, 0, 0, 0);

    // setup time_zero & rand_point
    set_zero_d(time_zero);
    make_vec_random_d(rand_point, numVars - 1); // remove extra variable

    // evaluate system
    if (regen->noChanges)
    { // evaluate all of the functions
      square_system_eval_data_d *SSED = (square_system_eval_data_d *)regen->square_d;
      evalProg_d(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, SSED->Prog);
      SSED = NULL;
    }
    else
    { // evaluate the square system using square_d
      eval_d(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, regen->square_d, regen->square_eval_d);
    }

    // allocate and setup mainCoeff_d
    init_mat_d(regen->mainCoeff_d, numFuncs, numVars);
    regen->mainCoeff_d->rows = numFuncs;
    regen->mainCoeff_d->cols = numVars;
    for (i = 0; i < numFuncs; i++)
    { // setup the coefficient for this function
      for (j = 0; j < numVars; j++)
        if (j == numVars - 1)
        { // setup for new hom variables
          get_comp_rand_d(&regen->mainCoeff_d->entry[i][j]);
        }
        else if (Jv->entry[i][j].r != 0 || Jv->entry[i][j].i != 0)
        { // variable is used
          varUsed[i][j] = 1;
          get_comp_rand_d(&regen->mainCoeff_d->entry[i][j]);
        }
        else
        { // variable is not used
          varUsed[i][j] = 0;
          set_zero_d(&regen->mainCoeff_d->entry[i][j]);
        }
    }

    // allocate coeff_d and initialize to zero
    regen->coeff_d = (comp_d ****)bmalloc(regen->num_funcs * sizeof(comp_d ***));
    for (i = 0; i < regen->num_funcs; i++)
    {
      regen->coeff_d[i] = (comp_d ***)bmalloc(num_var_gps * sizeof(comp_d **));
      for (j = 0; j < num_var_gps; j++)
      {
        regen->coeff_d[i][j] = (comp_d **)bmalloc(regen->degrees[i][j] * sizeof(comp_d *));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          regen->coeff_d[i][j][k] = (comp_d *)bmalloc(numVars * sizeof(comp_d));
          for (l = 0; l < numVars; l++)
            set_zero_d(regen->coeff_d[i][j][k][l]);
        }
      }
    }

    // now setup the appropriate coeff
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
        { // need to make the linears truly linears in only the variable groups - zero for the other coeff
          // this linear is for variable group j
          if (PPD->type[j])
          { // this is a regular variable group - thus it has a homogenous coordinate that was created
            beg = 0;
            for (l = 0; l < j; l++)
              beg += PPD->type[l];

            // check to see if this variable is used
            if (varUsed[i][beg])
            { // get random for hom coord for this variable group that is used
              get_comp_rand_d(regen->coeff_d[i][j][k][beg]);  
            }

            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              { // get random for variables that are in this variable group that are used
                get_comp_rand_d(regen->coeff_d[i][j][k][l+beg]);
              }
          }
          else
          { // this is a homogenous variable group - thus no homogenous coordinate was created
            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              { // get random for variables that are in this variable group that are used
                get_comp_rand_d(regen->coeff_d[i][j][k][l+beg]); 
              }
          }

          // setup for new hom variables
          get_comp_rand_d(regen->coeff_d[i][j][k][numVars - 1]);
        }

    // we enforce that if using 1-hom, main slice and first slice are the same
    if (num_var_gps == 1)
    { // setup first slice as main slice
      for (i = 0; i < numFuncs; i++)
        for (j = 0; j < numVars; j++)
        {
          set_d(regen->coeff_d[i][0][0][j], &regen->mainCoeff_d->entry[i][j]);
        }
    }

    // clear memory
    clear_point_d(rand_point);
    clear_mat_d(Jv);
    clear_eval_struct_d(e);
  }
  else if (MPType == 1)
  { // setup mainCoeff_mp & coeff_mp
    regen->mainCoeff_rat = NULL;
    regen->coeff_d = NULL;
    regen->coeff_rat = NULL;

    // compute the Jacobian matrix at a random point
    comp_mp time_zero;
    point_mp rand_point;
    mat_mp Jv;
    eval_struct_mp e;

    init_mp(time_zero);
    init_point_mp(rand_point, 0);
    init_mat_mp(Jv, 0, 0);
    init_eval_struct_mp(e, 0, 0, 0);

    // setup time_zero & rand_point
    set_zero_mp(time_zero);
    make_vec_random_mp(rand_point, numVars - 1); // remove extra variable

    if (regen->noChanges)
    { // evaluate all of the functions
      square_system_eval_data_mp *SSED = (square_system_eval_data_mp *)regen->square_mp;
      evalProg_mp(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, SSED->Prog);
      SSED = NULL;
    }
    else
    { // evaluate the square system using square_mp
      eval_mp(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, regen->square_mp, regen->square_eval_mp);
    }

    // allocate and setup mainCoeff_mp
    init_mat_mp(regen->mainCoeff_mp, numFuncs, numVars);
    regen->mainCoeff_mp->rows = numFuncs;
    regen->mainCoeff_mp->cols = numVars;
    for (i = 0; i < numFuncs; i++)
    { // setup the coefficient for this function
      for (j = 0; j < numVars; j++)
        if (j == numVars - 1)
        { // setup for new hom coordinate
          get_comp_rand_mp(&regen->mainCoeff_mp->entry[i][j]);
        }
        else if (!mpfr_zero_p(Jv->entry[i][j].r) || !mpfr_zero_p(Jv->entry[i][j].i))
        { // variable is used
          varUsed[i][j] = 1;
          get_comp_rand_mp(&regen->mainCoeff_mp->entry[i][j]);
        }
        else
        { // variable is not used
          varUsed[i][j] = 0;
          set_zero_mp(&regen->mainCoeff_mp->entry[i][j]);
        }
    }

    // allocate coeff_mp and initialize to zero
    regen->coeff_mp = (comp_mp ****)bmalloc(regen->num_funcs * sizeof(comp_mp ***));
    for (i = 0; i < regen->num_funcs; i++)
    {
      regen->coeff_mp[i] = (comp_mp ***)bmalloc(num_var_gps * sizeof(comp_mp **));
      for (j = 0; j < num_var_gps; j++)
      {
        regen->coeff_mp[i][j] = (comp_mp **)bmalloc(regen->degrees[i][j] * sizeof(comp_mp *));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          regen->coeff_mp[i][j][k] = (comp_mp *)bmalloc(numVars * sizeof(comp_mp));
          for (l = 0; l < numVars; l++)
          {
            init_mp(regen->coeff_mp[i][j][k][l]);
            set_zero_mp(regen->coeff_mp[i][j][k][l]);
          }
        }
      }
    }

    // now setup the appropriate coeff
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
        { // need to make the linears truly linears in only the variable groups - zero for the other coeff
          // this linear is for variable group j
          if (PPD->type[j])
          { // this is a regular variable group - thus it has a homogenous coordinate that was created
            beg = 0;
            for (l = 0; l < j; l++)
              beg += PPD->type[l];

            if (varUsed[i][beg])
            { // get random for hom coord for this variable group that is used
              get_comp_rand_mp(regen->coeff_mp[i][j][k][beg]); 
            }

            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              { // get random for variables that are in this variable group that are used
                get_comp_rand_mp(regen->coeff_mp[i][j][k][l+beg]);
              }
          }
          else
          { // this is a homogenous variable group - thus no homogenous coordinate was created
            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              { // get random for variables that are in this variable group that are used
                get_comp_rand_mp(regen->coeff_mp[i][j][k][l+beg]);
              }
          }

          // setup for new hom variables
          get_comp_rand_mp(regen->coeff_mp[i][j][k][numVars - 1]);
        }

    // we enforce that if using 1-hom, main slice and first slice are the same
    if (num_var_gps == 1)
    { // setup first slice as main slice
      for (i = 0; i < numFuncs; i++)
        for (j = 0; j < numVars; j++)
        {
          set_mp(regen->coeff_mp[i][0][0][j], &regen->mainCoeff_mp->entry[i][j]);
        }
    }

    // clear memory
    clear_mp(time_zero);
    clear_point_mp(rand_point);
    clear_mat_mp(Jv);
    clear_eval_struct_mp(e);
  }
  else
  { // setup mainCoeff_d, mainCoeff_mp, mainCoeff_rat, coeff_d, coeff_mp, coeff_rat

    // compute the Jacobian matrix at a random point
    comp_d time_zero;
    point_d rand_point;
    mat_d Jv;
    eval_struct_d e;

    init_point_d(rand_point, 0);
    init_mat_d(Jv, 0, 0);
    init_eval_struct_d(e, 0, 0, 0);

    // setup time_zero & rand_point
    set_zero_d(time_zero);
    make_vec_random_d(rand_point, numVars - 1); // remove extra variable

    if (regen->noChanges)
    { // evaluate all of the functions
      square_system_eval_data_d *SSED = (square_system_eval_data_d *)regen->square_d;
      evalProg_d(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, SSED->Prog);
      SSED = NULL;
    }
    else
    { // evaluate the square system using square_d
      eval_d(e.funcVals, e.parVals, e.parDer, Jv, e.Jp, rand_point, time_zero, regen->square_d, regen->square_eval_d);
    }

    // allocate and setup mainCoeff_d, mainCoeff_mp & mainCoeff_rat
    init_mat_d(regen->mainCoeff_d, numFuncs, numVars);
    init_mat_mp(regen->mainCoeff_mp, numFuncs, numVars);
    init_mat_rat(regen->mainCoeff_rat, numFuncs, numVars);
    regen->mainCoeff_d->rows = regen->mainCoeff_mp->rows = numFuncs;
    regen->mainCoeff_d->cols = regen->mainCoeff_mp->cols = numVars;
    for (i = 0; i < numFuncs; i++)
    { // setup the coefficient for this function
      for (j = 0; j < numVars; j++)
        if (j == numVars - 1)
        { // setup for new hom coordinate
          get_comp_rand_rat(&regen->mainCoeff_d->entry[i][numVars - 1], &regen->mainCoeff_mp->entry[i][numVars - 1], regen->mainCoeff_rat[i][numVars - 1], regen->curr_precision, AMP_max_prec, 0, 0);
        }
        else if (Jv->entry[i][j].r != 0 || Jv->entry[i][j].i != 0)
        { // variable is used
          varUsed[i][j] = 1;
          get_comp_rand_rat(&regen->mainCoeff_d->entry[i][j], &regen->mainCoeff_mp->entry[i][j], regen->mainCoeff_rat[i][j], regen->curr_precision, AMP_max_prec, 0, 0);
        }
        else
        { // variable is not used
          varUsed[i][j] = 0;
          set_zero_d(&regen->mainCoeff_d->entry[i][j]);
          set_zero_mp(&regen->mainCoeff_mp->entry[i][j]);
          set_zero_rat(regen->mainCoeff_rat[i][j]);
        }
    }

    // allocate coeff_d, coeff_mp & coeff_rat and initialize to zero
    regen->coeff_d = (comp_d ****)bmalloc(regen->num_funcs * sizeof(comp_d ***));
    regen->coeff_mp = (comp_mp ****)bmalloc(regen->num_funcs * sizeof(comp_mp ***));
    regen->coeff_rat = (mpq_t *****)bmalloc(regen->num_funcs * sizeof(mpq_t ****));
    for (i = 0; i < regen->num_funcs; i++)
    {
      regen->coeff_d[i] = (comp_d ***)bmalloc(num_var_gps * sizeof(comp_d **));
      regen->coeff_mp[i] = (comp_mp ***)bmalloc(num_var_gps * sizeof(comp_mp **));
      regen->coeff_rat[i] = (mpq_t ****)bmalloc(num_var_gps * sizeof(mpq_t ***));
      for (j = 0; j < num_var_gps; j++)
      {
        regen->coeff_d[i][j] = (comp_d **)bmalloc(regen->degrees[i][j] * sizeof(comp_d *));
        regen->coeff_mp[i][j] = (comp_mp **)bmalloc(regen->degrees[i][j] * sizeof(comp_mp *));
        regen->coeff_rat[i][j] = (mpq_t ***)bmalloc(regen->degrees[i][j] * sizeof(mpq_t **));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          regen->coeff_d[i][j][k] = (comp_d *)bmalloc(numVars * sizeof(comp_d));
          regen->coeff_mp[i][j][k] = (comp_mp *)bmalloc(numVars * sizeof(comp_mp));
          regen->coeff_rat[i][j][k] = (mpq_t **)bmalloc(numVars * sizeof(mpq_t *));
          for (l = 0; l < numVars; l++)
          {
            regen->coeff_rat[i][j][k][l] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            mpq_init(regen->coeff_rat[i][j][k][l][0]);
            mpq_init(regen->coeff_rat[i][j][k][l][1]);
            set_zero_rat(regen->coeff_rat[i][j][k][l]);

            init_mp(regen->coeff_mp[i][j][k][l]);
            set_zero_mp(regen->coeff_mp[i][j][k][l]);

            set_zero_d(regen->coeff_d[i][j][k][l]);
          }
        }
      }
    }

    // now setup the appropriate coeff
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
        { // need to make the linears truly linears in only the variable groups - zero for the other coeff
          // this linear is for variable group j
          if (PPD->type[j])
          { // this is a regular variable group - thus it has a homogenous coordinate that was created
            beg = 0;
            for (l = 0; l < j; l++)
              beg += PPD->type[l];

            // check to see if this variable is used
            if (varUsed[i][beg])
            { // get random for hom coord for this variable group
              get_comp_rand_rat(regen->coeff_d[i][j][k][beg], regen->coeff_mp[i][j][k][beg], regen->coeff_rat[i][j][k][beg], regen->curr_precision, AMP_max_prec, 0, 0);
            }

            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            // get random for variables that are in this variable group that are used
            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              {
                get_comp_rand_rat(regen->coeff_d[i][j][k][l + beg], regen->coeff_mp[i][j][k][l + beg], regen->coeff_rat[i][j][k][l + beg], regen->curr_precision, AMP_max_prec, 0, 0);
              }
          }
          else
          { // this is a homogenous variable group - thus no homogenous coordinate was created
            beg = PPD->num_var_gp;
            for (l = 0; l < j; l++)
              beg += PPD->size[l]; // add sizes to find the beginning of the variables for this variable group

            // get random for variables that are in this variable group that are used
            for (l = 0; l < PPD->size[j]; l++)
              if (varUsed[i][l+beg])
              {
                get_comp_rand_rat(regen->coeff_d[i][j][k][l + beg], regen->coeff_mp[i][j][k][l + beg], regen->coeff_rat[i][j][k][l + beg], regen->curr_precision, AMP_max_prec, 0, 0);
              }
          }

          // setup for new hom coordinate
          get_comp_rand_rat(regen->coeff_d[i][j][k][numVars - 1], regen->coeff_mp[i][j][k][numVars - 1], regen->coeff_rat[i][j][k][numVars - 1], regen->curr_precision, AMP_max_prec, 0, 0);
        }

    // we enforce that if using 1-hom, main slice and first slice are the same
    if (num_var_gps == 1)
    { // setup first slice as main slice
      for (i = 0; i < numFuncs; i++)
        for (j = 0; j < numVars; j++)
        {
          set_d(regen->coeff_d[i][0][0][j], &regen->mainCoeff_d->entry[i][j]);
          set_mp(regen->coeff_mp[i][0][0][j], &regen->mainCoeff_mp->entry[i][j]);
          set_rat(regen->coeff_rat[i][0][0][j], regen->mainCoeff_rat[i][j]);
        }
    }

    // clear memory
    clear_point_d(rand_point);
    clear_mat_d(Jv);
    clear_eval_struct_d(e);
  }

  // clear varUsed
  for (i = 0; i < numFuncs; i++)
    free(varUsed[i]);
  free(varUsed);

  return;
}

void setupRegenRandom_zero_dim(regen_t *regen, tracker_config_t *T, char *degreeName, preproc_data *PPD, square_system_eval_data_d *SSED_d, patch_eval_data_d *patch_d, square_system_eval_data_mp *SSED_mp, patch_eval_data_mp *patch_mp, int *startSub, int *endSub, int *startFunc, int *endFunc, int *startJvsub, int *endJvsub, int *startJv, int *endJv, int **subFuncsBelow)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the random numbers inside of regen for doing     *
* basic zero dimensional tracking                               *
\***************************************************************/
{
  int i, j, num_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  FILE *degIN = fopen(degreeName, "r"); // open the file to read in the degrees
  if (degIN == NULL)
  {
    printf("ERROR: '%s' does not exist!!\n", degreeName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // cp PPD into regen
  cp_preproc_data(&regen->PPD, PPD);

  // setup the basics in regen
  regen->curr_precision = T->Precision;
  regen->num_variables = T->numVars + 1; // adding a new hom variable
  regen->num_var_gps = num_var_gps;
  regen->num_funcs = T->numVars - num_var_gps;

  // verify the number of functions is correct
  if (((T->MPType == 0 || T->MPType == 2) && regen->num_funcs != SSED_d->size_r) || (T->MPType == 1 && regen->num_funcs != SSED_mp->size_r))
  { 
    printf("ERROR: The number of functions does not match when setting up regeneration!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // setup square system data
  regen->square_d = SSED_d;
  regen->square_mp = SSED_mp;
  regen->square_eval_d = square_system_eval_d;
  regen->square_eval_mp = square_system_eval_mp;
  regen->change_square_prec = change_square_prec;

  // setup noChanges
  if (T->MPType == 0 || T->MPType == 2)
    regen->noChanges = SSED_d->noChanges;
  else
    regen->noChanges = SSED_mp->noChanges;

  // setup instruction counts
  if (regen->noChanges)
  { // setup counts
    if (T->MPType == 0 || T->MPType == 2)
      regen->numSubFuncs = SSED_d->Prog->numSubfuncs;
    else
      regen->numSubFuncs = SSED_mp->Prog->numSubfuncs;

    regen->startSub = (int *)bmalloc(regen->numSubFuncs * sizeof(int));
    regen->endSub = (int *)bmalloc(regen->numSubFuncs * sizeof(int));
    regen->startFunc = (int *)bmalloc(regen->num_funcs * sizeof(int));
    regen->endFunc = (int *)bmalloc(regen->num_funcs * sizeof(int));
    regen->startJvsub = (int *)bmalloc(regen->numSubFuncs * sizeof(int));
    regen->endJvsub = (int *)bmalloc(regen->numSubFuncs * sizeof(int));
    regen->startJv = (int *)bmalloc(regen->num_funcs * sizeof(int));
    regen->endJv = (int *)bmalloc(regen->num_funcs * sizeof(int));
    for (i = 0; i < regen->numSubFuncs; i++)
    {
      regen->startSub[i] = startSub[i];
      regen->endSub[i] = endSub[i];
      regen->startJvsub[i] = startJvsub[i];
      regen->endJvsub[i] = endJvsub[i];
    }
    for (i = 0; i < regen->num_funcs; i++)
    {
      regen->startFunc[i] = startFunc[i];
      regen->endFunc[i] = endFunc[i];
      regen->startJv[i] = startJv[i];
      regen->endJv[i] = endJv[i];
    }

    if (regen->numSubFuncs > 0)
    { // setup subFuncsBelow
      regen->subFuncsBelow = (int **)bmalloc(regen->num_funcs * sizeof(int *));
      for (i = 0; i < regen->num_funcs; i++)
        regen->subFuncsBelow[i] = (int *)bmalloc(regen->numSubFuncs * sizeof(int));

      for (i = 0; i < regen->num_funcs; i++)
        for (j = 0; j < regen->numSubFuncs; j++)
          regen->subFuncsBelow[i][j] = subFuncsBelow[i][j];
    }
    else
      regen->subFuncsBelow = NULL;
  }
  else
  { // set to NULL
    regen->startSub = regen->endSub = regen->startFunc = regen->endFunc = regen->startJvsub = regen->endJvsub = regen->startJv = regen->endJv = NULL;
    regen->subFuncsBelow = NULL;
    regen->numSubFuncs = 0;
  }

  // setup patch
  if (T->MPType == 0)
  { // setup patchCoeff_d
    init_mat_d(regen->patchCoeff_d, num_var_gps, regen->num_variables);
    regen->patchCoeff_d->rows = num_var_gps; 
    regen->patchCoeff_d->cols = regen->num_variables;
    for (i = 0; i < num_var_gps; i++)
    { // [patch, -1]
      for (j = 0; j < patch_d->patchCoeff->cols; j++)
        set_d(&regen->patchCoeff_d->entry[i][j], &patch_d->patchCoeff->entry[i][j]);
      set_neg_one_d(&regen->patchCoeff_d->entry[i][j]);
    }
  }
  else if (T->MPType == 1)
  { // setup patchCoeff_mp
    init_mat_mp(regen->patchCoeff_mp, num_var_gps, regen->num_variables);
    regen->patchCoeff_mp->rows = num_var_gps;
    regen->patchCoeff_mp->cols = regen->num_variables;
    for (i = 0; i < num_var_gps; i++)
    { // [patch, -1]
      for (j = 0; j < patch_mp->patchCoeff->cols; j++)
        set_mp(&regen->patchCoeff_mp->entry[i][j], &patch_mp->patchCoeff->entry[i][j]);
      set_neg_one_mp(&regen->patchCoeff_mp->entry[i][j]);
    }
  }
  else
  { // setup patchCoeff_d, patchCoeff_mp & patchCoeff_rat
    init_mat_d(regen->patchCoeff_d, num_var_gps, regen->num_variables);
    init_mat_mp(regen->patchCoeff_mp, num_var_gps, regen->num_variables);
    init_mat_rat(regen->patchCoeff_rat, num_var_gps, regen->num_variables);

    regen->patchCoeff_d->rows = regen->patchCoeff_mp->rows = num_var_gps;
    regen->patchCoeff_d->cols = regen->patchCoeff_mp->cols = regen->num_variables;

    for (i = 0; i < num_var_gps; i++)
    { // [patch, -1]
      for (j = 0; j < patch_d->patchCoeff->cols; j++)
      {
        set_rat(regen->patchCoeff_rat[i][j], patch_mp->patchCoeff_rat[i][j]);
        rat_to_mp(&regen->patchCoeff_mp->entry[i][j], regen->patchCoeff_rat[i][j]);
        rat_to_d(&regen->patchCoeff_d->entry[i][j], regen->patchCoeff_rat[i][j]);
      }
      set_neg_one_rat(regen->patchCoeff_rat[i][j]);
      set_neg_one_mp(&regen->patchCoeff_mp->entry[i][j]);
      set_neg_one_d(&regen->patchCoeff_d->entry[i][j]);
    }
  }

  // setup main homogeneous variables
  if (T->MPType == 0)
  { // setup _d
    init_vec_d(regen->main_homVar_d, regen->num_variables);
    make_vec_random_d(regen->main_homVar_d, regen->num_variables);
  }
  else if (T->MPType == 1)
  { // setup _mp
    init_vec_mp(regen->main_homVar_mp, regen->num_variables);
    make_vec_random_mp(regen->main_homVar_mp, regen->num_variables);
  }
  else
  { // setup _d, _mp, _rat
    init_vec_d(regen->main_homVar_d, regen->num_variables);
    init_vec_mp(regen->main_homVar_mp, regen->num_variables);
    init_vec_rat(regen->main_homVar_rat, regen->num_variables);
    make_vec_random_rat(regen->main_homVar_d, regen->main_homVar_mp, regen->main_homVar_rat, regen->num_variables, regen->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // setup the degrees
  regen->degrees = (int **)bmalloc(regen->num_funcs * sizeof(int *));
  for (i = 0; i < regen->num_funcs; i++)
  { // allocate degrees[i]
    regen->degrees[i] = (int *)bmalloc(num_var_gps * sizeof(int));

    if (num_var_gps == 1)
    { // the functions could have been permuted, so we use SSED to get the degree
      if (regen->noChanges)
      { // simply use the degrees from SSED
        if (T->MPType == 0 || T->MPType == 2)
          regen->degrees[i][0] = SSED_d->new_degrees[i];
        else
          regen->degrees[i][0] = SSED_mp->new_degrees[i];
      }
      else
      { // sort from bottom to top on the functions to use
        if (T->MPType == 0 || T->MPType == 2)
          regen->degrees[i][0] = SSED_d->new_degrees[regen->num_funcs - 1 - i];
        else
          regen->degrees[i][0] = SSED_mp->new_degrees[regen->num_funcs - 1 - i];
      } 
    }
    else
    { // read in the m-hom degrees from degIN
      for (j = 0; j < num_var_gps; j++)
        fscanf(degIN, "%d\n", &regen->degrees[i][j]);
      fscanf(degIN, "\n"); // extra "new line" character
    }
  }
  // close degIN
  fclose(degIN);

  // setup gamma
  if (T->MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(regen->gamma_d);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    init_mp(regen->gamma_mp);
    get_comp_rand_mp(regen->gamma_mp);
  }
  else
  { // setup gamma_d, gamma_mp, gamma_rat
    get_comp_rand_rat(regen->gamma_d, regen->gamma_mp, regen->gamma_rat, regen->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // setup coeff
  regen_setup_coeff(T->MPType, T->AMP_max_prec, regen, num_var_gps, PPD);

  return;
}

void setupRegenLevels(regen_t *regen, tracker_config_t *T, char *depthName, int intrinsicCutoff)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the levels for the regeneration                  *
\***************************************************************/
{
  int i, j, *depths = NULL;
  FILE *DEPTH = fopen(depthName, "r");

  if (DEPTH == NULL)
  { // setup with each depth as 1
    regen->num_levels = regen->num_funcs;

    // allocate space
    regen->level = (regenLevel_t *)bmalloc(regen->num_levels * sizeof(regenLevel_t));

    for (i = 0; i < regen->num_levels; i++)
    { // setup level & depth
      regen->level[i].level = i;
      regen->level[i].depth = 1;
      // determine if we are using intrinsic slicing for this level
      if (regen->level[i].level + regen->level[i].depth > intrinsicCutoff)
      { // use extrinsic
        regen->level[i].useIntrinsicSlice = 0;
      }
      else
      { // use intrinsic
        regen->level[i].useIntrinsicSlice = 1;
      }
    }
  }
  else
  { // read in the depths
    regen->num_levels = 1;
    fscanf(DEPTH, "%d\n", &regen->num_levels); // read in the number of levels that will be needed

    // error checking
    if (regen->num_levels <= 0)
    {
      printf("ERROR: The number of levels (%d) must be > 0!\n", regen->num_levels);
      bexit(ERROR_CONFIGURATION);
    }

    // setup depths
    depths = (int *)bmalloc(regen->num_levels * sizeof(int));
    j = 0;
    for (i = 0; i < regen->num_levels; i++)
    {
      depths[i] = 1;
      fscanf(DEPTH, "%d\n", &depths[i]);
      if (depths[i] <= 0)
      {
        printf("ERROR: Each depth must be > 0 (%d)!\n", depths[i]);
        bexit(ERROR_CONFIGURATION);
      }
      j += depths[i]; // sum up the total depths
    }

    // verify the total depth is correct
    if (j != regen->num_funcs) 
    {
      printf("ERROR: The number of functions (%d) is not equal to the total depth (%d)!\n", regen->num_funcs, j);
      bexit(ERROR_CONFIGURATION);
    }

    // allocate space
    regen->level = (regenLevel_t *)bmalloc(regen->num_levels * sizeof(regenLevel_t));

    j = 0;
    for (i = 0; i < regen->num_levels; i++)
    { // setup level & depth
      regen->level[i].level = j;
      regen->level[i].depth = depths[i];
      // determine if we are using intrinsic slicing for this level
      if (regen->level[i].level + regen->level[i].depth > intrinsicCutoff)
      { // use extrinsic
        regen->level[i].useIntrinsicSlice = 0;
      }
      else
      { // use intrinsic
        regen->level[i].useIntrinsicSlice = 1;
      }

      j += depths[i];
    }

    free(depths);
    fclose(DEPTH);
  }

  return;
}

void setupFirstLevel(regen_t *regen, tracker_config_t *T, char *startName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first level for the regeneration             *
\***************************************************************/
{
  // setup the structures for the first level
  setupRegenLevelStructures(regen, T->MPType, 0, T->AMP_max_prec);
 
  // create the start points for the first level
  regen->level[0].num_paths = createFirstRegenStartPoints(regen, T, startName);

  return;
}

void setupRegenLevelStructures(regen_t *regen, int MPType, int level_num, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the structures for level_num of regen            *
\***************************************************************/
{
  if (regen->level[level_num].useIntrinsicSlice)
  { // need to setup B & p
    setupRegenIntrinsicSlice(regen, MPType, level_num, max_prec);
  }
  else
  { // do not need to setup B & p
    regen->level[level_num].B_d->rows = regen->level[level_num].B_d->cols = regen->level[level_num].B_mp->rows = regen->level[level_num].B_mp->cols = 0;
    regen->level[level_num].B_rat = NULL;
    regen->level[level_num].p_d->size = regen->level[level_num].p_mp->size = 0;
    regen->level[level_num].p_rat = NULL;
  }

  // initialize other values
  regen->level[level_num].num_sing = regen->level[level_num].num_nonsing = regen->level[level_num].num_inf 
    = regen->level[level_num].num_higher_dim = regen->level[level_num].num_bad = 0;

  return;
}

void setupRegenIntrinsicSlice(regen_t *regen, int MPType, int level_num, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup B & p for using an intrinsic slice               *
\***************************************************************/
{
  int i, j, start, finish, count;
  int numVars = regen->num_variables, num_var_gps = regen->num_var_gps;

  // setup start & finish
  start = regen->level[level_num].level + regen->level[level_num].depth;
  finish = regen->num_funcs;

  // all slices are of the form coeff*x - hom_var

  if (MPType == 0)
  { // setup B_d & p_d
    init_mat_d(regen->level[level_num].B_d, 0, 0);
    init_vec_d(regen->level[level_num].p_d, 0);

    vec_d b;
    mat_d tempMat, A, Q, R, P;
    double tol_pivot = 1e-15, tol_sign = 1e-20, largeChange = 1e14;

    init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

    // setup A
    count = finish - start + num_var_gps + 1; 
    init_mat_d(A, count, numVars);
    A->rows = count; 
    A->cols = numVars;
    // setup tempMat & b
    init_mat_d(tempMat, numVars, numVars);
    init_vec_d(b, numVars);
    tempMat->rows = tempMat->cols = b->size = numVars;

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < numVars; j++)
      { // function i, variable j of mainCoeff_d
        set_d(&A->entry[count][j], &regen->mainCoeff_d->entry[i][j]);
        // copy to tempMat
        set_d(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      // set b[count] to 0
      set_zero_d(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients
    for (i = 0; i < num_var_gps; i++)
    { // setup b[count] to 0
      set_zero_d(&b->coord[count]);
      for (j = 0; j < numVars; j++)
      { // copy patch coeff
        set_d(&A->entry[count][j], &regen->patchCoeff_d->entry[i][j]);
        // copy to tempMat
        set_d(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      count++;
    }
    // copy over main_homVar
    set_one_d(&b->coord[count]);
    for (j = 0; j < numVars; j++)
    { // copy coeff
      set_d(&A->entry[count][j], &regen->main_homVar_d->coord[j]);
      // copy to tempMat
      set_d(&tempMat->entry[count][j], &A->entry[count][j]);
    }
    count++;
    // put random values below
    for (count = count; count < numVars; count++)
    { // set b[count] to random value
      get_comp_rand_d(&b->coord[count]);
      for (j = 0; j < numVars; j++)
        get_comp_rand_d(&tempMat->entry[count][j]);
    }

    // compute p
    if (matrixSolve_d(regen->level[level_num].p_d, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem setting up the intrinsic slice.\n");
      bexit(ERROR_OTHER);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_d(A, A);
    QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    change_size_mat_d(regen->level[level_num].B_d, numVars, start);
    regen->level[level_num].B_d->rows = numVars;
    regen->level[level_num].B_d->cols = start;
    count = finish - start + num_var_gps + 1;

    for (i = 0; i < numVars; i++)
      for (j = 0; j < start; j++)
      {
        set_d(&regen->level[level_num].B_d->entry[i][j], &Q->entry[i][j+count]);
      }

    clear_vec_d(b);
    clear_mat_d(A); clear_mat_d(tempMat);
    clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
  }
  else if (MPType == 1)
  { // setup B_mp & p_mp
    init_mat_mp(regen->level[level_num].B_mp, 0, 0);
    init_vec_mp(regen->level[level_num].p_mp, 0);

    int num_digits = prec_to_digits(regen->curr_precision);
    size_t size;
    char *str = NULL;

    vec_mp b;
    mat_mp tempMat, A, Q, R, P;
    mpf_t tol_pivot, tol_sign, largeChange;

    // initialize MP
    mpf_init(tol_pivot); mpf_init(tol_sign); mpf_init(largeChange);
    init_mat_mp(Q, 0, 0); init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0);

    // setup A
    count = finish - start + num_var_gps + 1;
    init_mat_mp(A, count, numVars);
    A->rows = count;
    A->cols = numVars;
    // setup tempMat & b
    init_mat_mp(tempMat, numVars, numVars);
    init_vec_mp(b, numVars);
    tempMat->rows = tempMat->cols = b->size = numVars;

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

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < numVars; j++)
      { // function i, variable j of mainCoeff
        set_mp(&A->entry[count][j], &regen->mainCoeff_mp->entry[i][j]); 
        // copy to tempMat
        set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients
    for (i = 0; i < num_var_gps; i++)
    { // setup b[count] to 0
      set_zero_mp(&b->coord[count]);
      for (j = 0; j < numVars; j++)
      { // copy patch coeff
        set_mp(&A->entry[count][j], &regen->patchCoeff_mp->entry[i][j]);
        // copy to tempMat
        set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      count++;
    }
    // copy over main_homVar
    set_one_mp(&b->coord[count]);
    for (j = 0; j < numVars; j++)
    { // copy coeff
      set_mp(&A->entry[count][j], &regen->main_homVar_mp->coord[j]);
      // copy to tempMat
      set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
    }
    count++;
    // put random values below
    for (count = count; count < numVars; count++)
    { // set b[count] to random value
      get_comp_rand_mp(&b->coord[count]);
      for (j = 0; j < numVars; j++)
        get_comp_rand_mp(&tempMat->entry[count][j]);
    }

    // compute p 
    if (matrixSolve_mp(regen->level[level_num].p_mp, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem setting up the intrinsic slice.\n");
      bexit(ERROR_OTHER);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);
    change_size_mat_mp(regen->level[level_num].B_mp, numVars, start);
    regen->level[level_num].B_mp->rows = numVars; 
    regen->level[level_num].B_mp->cols = start;
    count = finish - start + num_var_gps + 1;

    for (i = 0; i < regen->level[level_num].B_mp->rows; i++)
      for (j = 0; j < regen->level[level_num].B_mp->cols; j++)
      {
        set_mp(&regen->level[level_num].B_mp->entry[i][j], &Q->entry[i][j+count]);
      }

    // clear MP
    mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
    clear_vec_mp(b);
    clear_mat_mp(A); clear_mat_mp(tempMat); 
    clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

    free(str);
  }
  else
  { // setup B_d, B_mp, B_rat, p_d, p_mp & p_rat
    init_mat_d(regen->level[level_num].B_d, 0, 0);
    init_vec_d(regen->level[level_num].p_d, 0);
    init_mat_mp(regen->level[level_num].B_mp, 0, 0);
    init_vec_mp(regen->level[level_num].p_mp, 0);

    int num_digits = prec_to_digits(max_prec);
    size_t size;
    char *str = NULL;

    vec_mp tempVec, b;
    mat_mp tempMat, A, Q, R, P;
    mpf_t tol_pivot, tol_sign, largeChange;

    // initialize MP
    mpf_init2(tol_pivot, max_prec); mpf_init2(tol_sign, max_prec); mpf_init2(largeChange, max_prec);
    init_vec_mp2(tempVec, 0, max_prec);
    init_mat_mp2(Q, 0, 0, max_prec); init_mat_mp2(R, 0, 0, max_prec); init_mat_mp2(P, 0, 0, max_prec);

    // setup A
    count = finish - start + num_var_gps + 1;
    init_mat_mp2(A, count, numVars, max_prec);
    A->rows = count; 
    A->cols = numVars;
    // setup tempMat & b
    init_mat_mp2(tempMat, numVars, numVars, max_prec);
    init_vec_mp2(b, numVars, max_prec);
    tempMat->rows = tempMat->cols = b->size = numVars;
  
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

    // copy over the coefficients
    count = 0;
    for (i = start; i < finish; i++)
    { // copy over the coefficients for linear i
      for (j = 0; j < numVars; j++)
      { // function i, variable j of mainCoeff
        rat_to_mp(&A->entry[count][j], regen->mainCoeff_rat[i][j]);
        // copy to tempMat
        set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      count++;
    }
    // copy over the patch coefficients
    for (i = 0; i < num_var_gps; i++)
    { // set b[count] to 0
      set_zero_mp(&b->coord[count]);
      for (j = 0; j < numVars; j++)
      { // copy patch coeff
        rat_to_mp(&A->entry[count][j], regen->patchCoeff_rat[i][j]);
        // copy to tempMat
        set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
      }
      count++;
    }
    // copy over main_homVar
    set_one_mp(&b->coord[count]);
    for (j = 0; j < numVars; j++)
    { // copy coeff
      rat_to_mp(&A->entry[count][j], regen->main_homVar_rat[j]);
      // copy to tempMat
      set_mp(&tempMat->entry[count][j], &A->entry[count][j]);
    }
    count++;
    // put random values below
    comp_d tempComp;
    mpq_t tempMPQ[2];
    mpq_init(tempMPQ[0]); mpq_init(tempMPQ[1]);
    for (count = count; count < numVars; count++)
    { // set b[count] to random value
      get_comp_rand_rat(tempComp, &b->coord[count], tempMPQ, max_prec, max_prec, 0, 0);
      for (j = 0; j < numVars; j++)
        get_comp_rand_rat(tempComp, &tempMat->entry[count][j], tempMPQ, max_prec, max_prec, 0, 0);
    }
    mpq_clear(tempMPQ[0]); mpq_clear(tempMPQ[1]);

    // set the global prec to the maximum precision
    initMP(max_prec);

    // solve for tempVec
    if (matrixSolve_mp(tempVec, tempMat, b))
    { // this should never happen!
      printf("ERROR: Problem setting up the intrinsic slice.\n");
      bexit(ERROR_OTHER);
    }

    // copy tempVec to p
    change_size_vec_d(regen->level[level_num].p_d, tempVec->size);
    change_size_vec_mp(regen->level[level_num].p_mp, tempVec->size);
    init_vec_rat(regen->level[level_num].p_rat, tempVec->size);
    regen->level[level_num].p_d->size = regen->level[level_num].p_mp->size = tempVec->size;
    for (j = 0; j < tempVec->size; j++)
    { // copy tempVec to p_d, p_mp & p_rat
      mpf_t_to_rat(regen->level[level_num].p_rat[j][0], tempVec->coord[j].r);
      mpf_t_to_rat(regen->level[level_num].p_rat[j][1], tempVec->coord[j].i);
      mpf_set_q(regen->level[level_num].p_mp->coord[j].r, regen->level[level_num].p_rat[j][0]);
      mpf_set_q(regen->level[level_num].p_mp->coord[j].i, regen->level[level_num].p_rat[j][1]);
      regen->level[level_num].p_d->coord[j].r = mpq_get_d(regen->level[level_num].p_rat[j][0]);
      regen->level[level_num].p_d->coord[j].i = mpq_get_d(regen->level[level_num].p_rat[j][1]);
    }

    // now, to find B, do a QR decomposition on A^T and B is the extra columns of Q
    transpose_mp(A, A);
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 0);

    // allocate space and then setup B
    change_size_mat_d(regen->level[level_num].B_d, numVars, start);
    change_size_mat_mp(regen->level[level_num].B_mp, numVars, start);
    init_mat_rat(regen->level[level_num].B_rat, numVars, start);
    regen->level[level_num].B_d->rows = regen->level[level_num].B_mp->rows = numVars;
    regen->level[level_num].B_d->cols = regen->level[level_num].B_mp->cols = start;
    count = finish - start + num_var_gps + 1;

    for (i = 0; i < regen->level[level_num].B_d->rows; i++)
      for (j = 0; j < regen->level[level_num].B_d->cols; j++)
      { // copy Q to B_d, B_mp & B_rat 
        mpf_t_to_rat(regen->level[level_num].B_rat[i][j][0], Q->entry[i][j+count].r);
        mpf_t_to_rat(regen->level[level_num].B_rat[i][j][1], Q->entry[i][j+count].i);
        mpf_set_q(regen->level[level_num].B_mp->entry[i][j].r, regen->level[level_num].B_rat[i][j][0]);
        mpf_set_q(regen->level[level_num].B_mp->entry[i][j].i, regen->level[level_num].B_rat[i][j][1]);
        regen->level[level_num].B_d->entry[i][j].r = mpq_get_d(regen->level[level_num].B_rat[i][j][0]);
        regen->level[level_num].B_d->entry[i][j].i = mpq_get_d(regen->level[level_num].B_rat[i][j][1]);
      }

    // set back to the current precision
    initMP(regen->curr_precision);

    // clear MP
    mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
    clear_vec_mp(tempVec); clear_vec_mp(b);
    clear_mat_mp(A); clear_mat_mp(tempMat); 
    clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

    free(str);
  }

  return;
}

int isCompatible(int *loc, int *P, int size, int **mhomDeg, int num_var_gps, int *var_gp_sizes, int *tempCount)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - compatible, 0 - otherwise                  *
* NOTES: determine if loc is compatible with the function       *
\***************************************************************/
{
  int i, j;

  // initialize tempCount to 1 - patch
  for (i = 0; i < num_var_gps; i++)
    tempCount[i] = 1;

  // count the number of entries for each variable group
  for (i = 0; i < size; i++)
  {
    j = loc[i];

    if (mhomDeg[P[i]][j] == 0)
    { // does not depend on this variable group!
      return 0;
    }

    // increment this variable group
    tempCount[j]++;

    // make sure that the number of entries in this spot <= var_gp_sizes
    if (tempCount[j] > var_gp_sizes[j])
      return 0;
  }

  return 1;
}

int createFirstRegenStartPoints(regen_t *regen, tracker_config_t *T, char *startName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of start points for first level         *
* NOTES: setup the start points for the first level of regen    *
\***************************************************************/
{
  FILE *START = fopen(startName, "w");
  int j, goodCount, depth = regen->level[0].depth, good_size = 0;
  int **goodLoc = NULL;

  // create the possible decompostions
  createDecomp(&goodLoc, &good_size, regen, 0, depth);

  // now, we have all of the good partitions in goodLoc - setup the start points

  // leave room for the number of start points in START
  fprintf(START, "                                                 \n\n");

  // setup the start points based on goodLoc
  goodCount = createRegenStartPts(regen, T->MPType, good_size, goodLoc, START);

  // setup the number of start points
  regen->level[0].num_paths = goodCount;

  // print the relevant data to START that can be used to rerun this exact same problem
  printRegenRelevantData(regen, T->MPType, 0, START);

  // rewind START
  rewind(START);

  // print the number of start points to START
  fprintf(START, "%d", goodCount);

  // close START
  fclose(START);

  // release memory
  for (j = good_size - 1; j >= 0; j--)
    free(goodLoc[j]);
  free(goodLoc);

  return goodCount;
}

void createDecompExtensions(int ****decompExt, int **numExt, regen_t *regen, int **decomp, int numDecomp, int decompSize, int newSize)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: takes a decomposition and computes its extensions      *
\***************************************************************/
{
  int i, j, cont, depth = newSize - decompSize;
  int num_var_gps = regen->num_var_gps;

  // initialize
  *decompExt = (int ***)bmalloc(numDecomp * sizeof(int **));
  *numExt = (int *)bmalloc(numDecomp * sizeof(int));
  for (i = 0; i < numDecomp; i++)
  {
    (*decompExt)[i] = NULL;
    (*numExt)[i] = 0;
  }

  if (num_var_gps == 1)
  { // each one only creates 1 new decomposition!
    for (i = 0; i < numDecomp; i++)
    { // setup numExt
      (*numExt)[i] = 1;
      // setup decompExt
      (*decompExt)[i] = (int **)bmalloc(1 * sizeof(int *));
      (*decompExt)[i][0] = (int *)bmalloc(newSize * sizeof(int));
      for (j = 0; j < newSize; j++)
        (*decompExt)[i][0][j] = 0; // 1st variable group is only variable group!
    }
  }
  else
  { // setup the basic data
    int *func_gp_count = (int *)bmalloc(depth * sizeof(int)); // number of variable groups used by each function
    int *curr_func_gp_count = (int *)bmalloc(depth * sizeof(int)); // current index for the function
    int **func_var_gps = (int **)bmalloc(depth * sizeof(int *)); // the variable groups used by the functions
    int *loc = (int *)bmalloc(newSize * sizeof(int)), *perm = (int *)bmalloc(newSize * sizeof(int));
    int *tempCount = (int *)bmalloc(num_var_gps * sizeof(int)), *var_gp_sizes = (int *)bmalloc(num_var_gps * sizeof(int));

    // setup sizes of variable groups
    for (i = 0; i < num_var_gps; i++)
      var_gp_sizes[i] = regen->PPD.size[i] + regen->PPD.type[i]; 

    // setup perm to identity
    for (i = 0; i < newSize; i++)
      perm[i] = i;

    // look through degrees so that we can see which variable groups were used in which functions
    for (i = 0; i < depth; i++)
    { // find the number of groups used for this function
      func_gp_count[i] = 0;
      for (j = 0; j < num_var_gps; j++)
        if (regen->degrees[i + decompSize][j] > 0)
          func_gp_count[i]++;

      // allocate func_var_gps
      func_var_gps[i] = (int *)bmalloc(func_gp_count[i] * sizeof(int));
      curr_func_gp_count[i] = 0;
      for (j = 0; j < num_var_gps; j++)
        if (regen->degrees[i + decompSize][j] > 0)
        {
          func_var_gps[i][curr_func_gp_count[i]] = j;
          curr_func_gp_count[i]++;
        }
    }

    // loop over the decompositions
    for (i = 0; i < numDecomp; i++)
    { // setup the top of loc
      for (j = 0; j < decompSize; j++)
        loc[j] = decomp[i][j];

      // initialize for loops
      for (j = 0; j < depth; j++)
      { // setup for first variable group used
        curr_func_gp_count[j] = 0;
        loc[j + decompSize] = func_var_gps[j][curr_func_gp_count[j]];
      }

      // loop through the partitions
      cont = 1;
      while (cont)
      {
        if (isCompatible(loc, perm, newSize, regen->degrees, num_var_gps, var_gp_sizes, tempCount))
        { // loc is good - add to decompExt
          (*decompExt)[i] = (int **)brealloc((*decompExt)[i], ((*numExt)[i] + 1) * sizeof(int *));
          (*decompExt)[i][(*numExt)[i]] = (int *)bmalloc(newSize * sizeof(int));
          for (j = 0; j < newSize; j++)
            (*decompExt)[i][(*numExt)[i]][j] = loc[j];

          // increment the number
          (*numExt)[i]++;
        }
        // update loc to 'next' partition
        for (j = depth - 1; j >= 0; j--)
        { // check to see if we are at the top
          if (curr_func_gp_count[j] == func_gp_count[j] - 1)
          { // set to 0
            curr_func_gp_count[j] = 0;
            loc[j + decompSize] = func_var_gps[j][curr_func_gp_count[j]];
          }
          else
          { // increment this location and exit loop
            curr_func_gp_count[j]++;
            loc[j + decompSize] = func_var_gps[j][curr_func_gp_count[j]];
            j = -10;
          }
        }
        // check to see if all are now back ot zero
        if (j == -1)
          cont = 0; // exit loop
      }
    }

    // clear memory
    for (i = depth - 1; i >= 0; i--)
      free(func_var_gps[i]);
    free(func_var_gps);
    free(func_gp_count);
    free(curr_func_gp_count);
    free(loc); free(perm);
    free(var_gp_sizes); free(tempCount); 
  } 

  return;
}

void createDecomp(int ***goodLoc, int *good_size, regen_t *regen, int startFunc, int endFunc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the decomp for [startFunc,endFunc)             *
\***************************************************************/
{
  int i, j, cont = 1, num_var_gps = regen->num_var_gps;
  int depth = endFunc - startFunc;
  int *var_gp_sizes = (int *)bmalloc(num_var_gps * sizeof(int));
  int *perm = (int *)bmalloc(depth * sizeof(int));
  int *tempCount = (int *)bmalloc(num_var_gps * sizeof(int));
  int *loc = (int *)bmalloc(depth * sizeof(int));

  // seutp var_gp_sizes
  for (i = 0; i < num_var_gps; i++)
   var_gp_sizes[i] = regen->PPD.size[i] + regen->PPD.type[i];

  for (i = 0; i < depth; i++)
  {
    perm[i] = i; // identity permutation
    loc[i] = 0;  // initialize to all 0s
  }

  // initialize goodLoc & good_size
  *goodLoc = NULL;
  *good_size = 0;

  if (num_var_gps == 1)
  { // only 1 way to partition things
    *good_size = 1;
    *goodLoc = (int **)bmalloc(1 * sizeof(int *));
    (*goodLoc)[0] = (int *)bmalloc(depth * sizeof(int));
    for (i = 0; i < depth; i++)
      (*goodLoc)[0][i] = 0;
  }
  else
  { // multihomogeneous
    int *func_gp_count = (int *)bmalloc(depth * sizeof(int));
    int *curr_func_gp_count = (int *)bmalloc(depth * sizeof(int));
    int **func_var_gps = (int **)bmalloc(depth * sizeof(int *));
    
    // see which variable groups are used in which function
    for (i = 0; i < depth; i++)
    { // find the number of groups used for this funcion
      func_gp_count[i] = 0;
      for (j = 0; j < num_var_gps; j++)
        if (regen->degrees[startFunc + i][j] > 0)
          func_gp_count[i]++;

      // allocate func_var_gps
      func_var_gps[i] = (int *)bmalloc(func_gp_count[i] * sizeof(int));
      curr_func_gp_count[i] = 0;
      for (j = 0; j < num_var_gps; j++)
        if (regen->degrees[startFunc + i][j] > 0)
        {
          func_var_gps[i][curr_func_gp_count[i]] = j;
          curr_func_gp_count[i]++;
        }

      // setup for first variable group used
      curr_func_gp_count[i] = 0;
      loc[i] = func_var_gps[i][curr_func_gp_count[i]];
    }

    // loop through the partitions
    while (cont)
    {
      if (isCompatible(loc, perm, depth, &regen->degrees[startFunc], num_var_gps, var_gp_sizes, tempCount))
      { // loc is good - store to goodLoc
        *goodLoc = (int **)brealloc(*goodLoc, (*good_size + 1) * sizeof(int *));
        (*goodLoc)[*good_size] = (int *)bmalloc(depth * sizeof(int));
        for (i = 0; i < depth; i++)
          (*goodLoc)[*good_size][i] = loc[i];

        // increment the number of good locations
        (*good_size)++;
      }
      // update loc to "next" partition
      for (j = depth - 1; j >= 0; j--)
      { // check to see if we are at the top
        if (curr_func_gp_count[j] == func_gp_count[j] - 1)
        { // set to 0
          curr_func_gp_count[j] = 0;
          loc[j] = func_var_gps[j][curr_func_gp_count[j]];
        }
        else
        { // increment this location and exit loop
          curr_func_gp_count[j]++;
          loc[j] = func_var_gps[j][curr_func_gp_count[j]];
          j = -10;
        }
      }
      // check to see if all are now back to zero
      if (j == -1)
        cont = 0; // exit loop
    }

    // free the memory
    for (i = depth - 1; i >= 0; i--)
      free(func_var_gps[i]);
    free(func_var_gps);
    free(func_gp_count);
    free(curr_func_gp_count);
  }

  // clear memory
  free(var_gp_sizes);
  free(perm);
  free(tempCount);
  free(loc);

  return;
}

void verifyDecompExtensions(int ****decomp, int **num_decomp, regen_t *regen, int MPType, int num, int size, mat_d A_d, mat_mp A_mp, vec_d b_d, vec_mp b_mp, vec_d tempVec_d, vec_mp tempVec_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: update decomp with the ones that can actually lead to  *
* solutions                                                     *
\***************************************************************/
{
  int i, j, k, l, var_gp, num_vars = regen->num_variables, num_funcs = regen->num_funcs;
  int ***good_decomp = (int ***)bmalloc(num * sizeof(int **));
  int *good_count = (int *)bmalloc(num * sizeof(int));

  for (i = 0; i < num; i++)
  {
    good_decomp[i] = NULL;
    good_count[i] = 0;
  }

  if (MPType == 0 || MPType == 2)
  { // verify in double precision

    // verify size of A & b
    change_size_mat_d(A_d, num_vars, num_vars);
    change_size_vec_d(b_d, num_vars);
    A_d->rows = A_d->cols = b_d->size = num_vars;

    // the bottom of A & b is fixed - main slices and patches
    for (i = size; i < num_funcs; i++)
    { // seutp main slices
      set_d(&b_d->coord[i], &regen->mainCoeff_d->entry[i][num_vars]); // rhs
      for (j = 0; j < num_vars; j++)
        set_d(&A_d->entry[i][j], &regen->mainCoeff_d->entry[i][j]);
    }
    for (i = 0; i < regen->num_var_gps; i++)
    { // setup patches
      k = i + num_funcs;
      set_one_d(&b_d->coord[k]);
      for (j = 0; j < num_vars; j++)
        set_d(&A_d->entry[k][j], &regen->patchCoeff_d->entry[i][j]);    
    }

    // loop over the decompositions
    for (i = 0; i < num; i++)
    { // initialize
      good_count[i] = 0;
 
      for (j = 0; j < (*num_decomp)[i]; j++)
      { // setup A & b for the decomp with degree 0
        for (k = 0; k < size; k++)
        { // setup coefficient slices
          var_gp = (*decomp)[i][j][k];
          set_d(&b_d->coord[k], regen->coeff_d[k][var_gp][0][num_vars]);
          for (l = 0; l < num_vars; l++)
            set_d(&A_d->entry[k][l], regen->coeff_d[k][var_gp][0][l]);
        }

        if (!matrixSolve_d(tempVec_d, A_d, b_d))
        { // this decomposition is good
          l = good_count[i];
          good_decomp[i] = (int **)brealloc(good_decomp[i], (l + 1) * sizeof(int *));
          good_decomp[i][l] = (int *)bmalloc(size * sizeof(int));
          for (k = 0; k < size; k++)
            good_decomp[i][l][k] = (*decomp)[i][j][k];

          // increment the number
          good_count[i]++;
        }
        // clear decomp[i][j]
        free((*decomp)[i][j]);
      }
      // clear decomp[i]
      free((*decomp)[i]);
    }
    // clear decomp & num_decomp
    free((*decomp));
    free((*num_decomp));
  }
  else
  { // verify in multi precision

    // verify size of A & b
    change_size_mat_mp(A_mp, num_vars, num_vars);
    change_size_vec_mp(b_mp, num_vars);
    A_mp->rows = A_mp->cols = b_mp->size = num_vars;

    // the bottom of A & b is fixed - main slices and patches
    for (i = size; i < num_funcs; i++)
    { // seutp main slices
      set_mp(&b_mp->coord[i], &regen->mainCoeff_mp->entry[i][num_vars]); // rhs
      for (j = 0; j < num_vars; j++)
        set_mp(&A_mp->entry[i][j], &regen->mainCoeff_mp->entry[i][j]);
    }
    for (i = 0; i < regen->num_var_gps; i++)
    { // setup patches
      k = i + num_funcs;
      set_one_mp(&b_mp->coord[k]);
      for (j = 0; j < num_vars; j++)
        set_mp(&A_mp->entry[k][j], &regen->patchCoeff_mp->entry[i][j]);
    }

    // loop over the decompositions
    for (i = 0; i < num; i++)
    { // initialize
      good_count[i] = 0;

      for (j = 0; j < (*num_decomp)[i]; j++)
      { // setup A & b for the decomp with degree 0
        for (k = 0; k < size; k++)
        { // setup coefficient slices
          var_gp = (*decomp)[i][j][k];
          set_mp(&b_mp->coord[k], regen->coeff_mp[k][var_gp][0][num_vars]);
          for (l = 0; l < num_vars; l++)
            set_mp(&A_mp->entry[k][l], regen->coeff_mp[k][var_gp][0][l]);
        }

        if (!matrixSolve_mp(tempVec_mp, A_mp, b_mp))
        { // this decomposition is good
          l = good_count[i];
          good_decomp[i] = (int **)brealloc(good_decomp[i], (l + 1) * sizeof(int *));
          good_decomp[i][l] = (int *)bmalloc(size * sizeof(int));
          for (k = 0; k < size; k++)
            good_decomp[i][l][k] = (*decomp)[i][j][k];

          // increment the number
          good_count[i]++;
        }
        // clear decomp[i][j]
        free((*decomp)[i][j]);
      }
      // clear decomp[i]
      free((*decomp)[i]);
    }
    // clear decomp & num_decomp
    free((*decomp));
    free((*num_decomp));
  }

  // setup decomp & num_decmop as good_decomp & good_count
  *decomp = good_decomp;
  *num_decomp = good_count;
  good_decomp = NULL;
  good_count = NULL;

  return;
}

int createRegenStartPtsDecomp_d(regen_t *regen, int *decomp, int depth, mat_d A, vec_d b, vec_d tempPt, mat_d B_transpose, FILE *START)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of start points actually printed        *
* NOTES: solves for the start points and prints them to START   *
\***************************************************************/
{
  int i, j, k, var_gp, count = 0, num_vars = regen->num_variables;
  int degMult = 1, *degrees = (int *)bmalloc(depth * sizeof(int));

  // setup A & b for the decomp with degree = 0
  for (j = 0; j < depth; j++)
  { // setup coefficient slices
    var_gp = decomp[j];
    degrees[j] = 0;
    degMult *= regen->degrees[j][var_gp];
    set_zero_d(&b->coord[j]);
    for (k = 0; k < num_vars; k++)
      set_d(&A->entry[j][k], regen->coeff_d[j][var_gp][degrees[j]][k]);
  }

  // solve for the start point
  if (!matrixSolve_d(tempPt, A, b))
  { // we have a possible start point and decomp - verify new hom coord is non-zero
    if (d_abs_d(&tempPt->coord[num_vars-1]) >= ZERO_HOM_COORD_D)
    { // we have a good start point and a good decomp
      count++;

      if (regen->level[0].useIntrinsicSlice)
      { // convert tempPt to intrinsic coordinates
        extrinsicToIntrinsic_d(tempPt, tempPt, B_transpose, regen->level[0].p_d);
      }

      // print to START
      printPointLinearDegree(START, tempPt, NULL, 52, decomp, degrees, depth);

      // loop over the degrees that are possible
      for (i = 1; i < degMult; i++)
      { // update to next degree
        for (j = depth - 1; j >= 0; j--)
          if (degrees[j] == regen->degrees[j][decomp[j]] - 1)
          { // set to 0
            degrees[j] = 0;
          }
          else
          { // increment and exit
            degrees[j]++;
            j = -10;
          }

        // seutp A & b
        for (j = 0; j < depth; j++)
        { // setup coefficient slices
          var_gp = decomp[j];
          set_zero_d(&b->coord[j]);
          for (k = 0; k < num_vars; k++)
            set_d(&A->entry[j][k], regen->coeff_d[j][var_gp][degrees[j]][k]);
        }

        if (!matrixSolve_d(tempPt, A, b))
        { // we have a good start point and a good decomp
          count++;

          if (regen->level[0].useIntrinsicSlice)
          { // convert tempPt to intrinsic coordinates
            extrinsicToIntrinsic_d(tempPt, tempPt, B_transpose, regen->level[0].p_d);
          }

          // print to START
          printPointLinearDegree(START, tempPt, NULL, 52, decomp, degrees, depth);
        }
        else
        { // error - should never happen!
          printf("ERROR: Problem setting up a start point.\n");
          bexit(ERROR_OTHER);
        }
      }
    }
  }

  free(degrees);

  return count;
}

int createRegenStartPtsDecomp_mp(regen_t *regen, int *decomp, int depth, mat_mp A, vec_mp b, vec_mp tempPt, mat_mp B_transpose, FILE *START)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of start points actually printed        *
* NOTES: solves for the start points and prints them to START   *
\***************************************************************/
{
  int i, j, k, var_gp, count = 0, num_vars = regen->num_variables;
  int degMult = 1, *degrees = (int *)bmalloc(depth * sizeof(int));

  // setup A & b for the decomp with degree = 0
  for (j = 0; j < depth; j++)
  { // setup coefficient slices
    var_gp = decomp[j];
    degrees[j] = 0;
    degMult *= regen->degrees[j][var_gp];
    set_zero_mp(&b->coord[j]);
    for (k = 0; k < num_vars; k++)
      set_mp(&A->entry[j][k], regen->coeff_mp[j][var_gp][degrees[j]][k]);
  }

  // solve for the start point
  if (!matrixSolve_mp(tempPt, A, b))
  { // we have a possible start point and decomp - verify new hom coord is non-zero
    if (d_abs_mp(&tempPt->coord[num_vars-1]) >= ZERO_HOM_COORD_MP)
    { // we have a good start point and a good decomp
      count++;

      if (regen->level[0].useIntrinsicSlice)
      { // convert tempPt to intrinsic coordinates
        extrinsicToIntrinsic_mp(tempPt, tempPt, B_transpose, regen->level[0].p_mp);
      }

      // print to START
      printPointLinearDegree(START, NULL, tempPt, regen->curr_precision, decomp, degrees, depth);

      // loop over the degrees that are possible
      for (i = 1; i < degMult; i++)
      { // update to next degree
        for (j = depth - 1; j >= 0; j--)
          if (degrees[j] == regen->degrees[j][decomp[j]] - 1)
          { // set to 0
            degrees[j] = 0;
          }
          else
          { // increment and exit
            degrees[j]++;
            j = -10;
          }

        // seutp A & b
        for (j = 0; j < depth; j++)
        { // setup coefficient slices
          var_gp = decomp[j];
          set_zero_mp(&b->coord[j]);
          for (k = 0; k < num_vars; k++)
            set_mp(&A->entry[j][k], regen->coeff_mp[j][var_gp][degrees[j]][k]);
        }

        if (!matrixSolve_mp(tempPt, A, b))
        { // we have a good start point and a good decomp
          count++;

          if (regen->level[0].useIntrinsicSlice)
          { // convert tempPt to intrinsic coordinates
            extrinsicToIntrinsic_mp(tempPt, tempPt, B_transpose, regen->level[0].p_mp);
          }

          // print to START
          printPointLinearDegree(START, NULL, tempPt, regen->curr_precision, decomp, degrees, depth);
        }
        else
        { // error - should never happen!
          printf("ERROR: Problem setting up a start point.\n");
          bexit(ERROR_OTHER);
        }
      }
    }
  }

  free(degrees);

  return count;
}

int createRegenStartPts(regen_t *regen, int MPType, int count, int **decomp, FILE *START)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: number of start points actually printed        *
* NOTES: solves for the start points and prints them to START   *
\***************************************************************/
{
  int i, j, k, depth = regen->level[0].depth;
  int num_funcs = regen->num_funcs, num_vars = regen->num_variables;
  int actualCount = 0;

  if (MPType == 0 || MPType == 2)
  { // find & print the start points in double precision
    mat_d A, B_transpose;
    vec_d b, tempPt;

    init_mat_d(A, num_vars, num_vars);
    init_vec_d(b, num_vars);
    init_vec_d(tempPt, num_vars);
    A->rows = A->cols = b->size = tempPt->size = num_vars;

    if (regen->level[0].useIntrinsicSlice)
    { // setup B_transpose
      init_mat_d(B_transpose, regen->level[0].B_d->cols, regen->level[0].B_d->rows);
      transpose_d(B_transpose, regen->level[0].B_d);
    }

    // the bottom of A & b is fixed - main slices and patches
    for (j = depth; j < num_funcs; j++)
    { // setup main slices
      set_zero_d(&b->coord[j]); // rhs
      for (k = 0; k < num_vars; k++)
        set_d(&A->entry[j][k], &regen->mainCoeff_d->entry[j][k]);
    }
    for (j = 0; j < regen->num_var_gps; j++)
    { // copy patches
      i = j + num_funcs;
      set_zero_d(&b->coord[i]); // rhs
      for (k = 0; k < num_vars; k++)
        set_d(&A->entry[i][k], &regen->patchCoeff_d->entry[j][k]);
    }
    // new hom patch
    i = num_funcs + regen->num_var_gps;
    set_one_d(&b->coord[i]);
    for (k = 0; k < num_vars; k++)
      set_d(&A->entry[i][k], &regen->main_homVar_d->coord[k]);

    // loop over the decompositions
    for (i = 0; i < count; i++)
    { // compute the start points related to this decomposition
      actualCount += createRegenStartPtsDecomp_d(regen, decomp[i], depth, A, b, tempPt, B_transpose, START);
    }

    clear_mat_d(A); 
    clear_vec_d(b); clear_vec_d(tempPt);
    if (regen->level[0].useIntrinsicSlice)
      clear_mat_d(B_transpose);
  }
  else
  { // find & print the start points in fixed multi precision
    mat_mp A, B_transpose;
    vec_mp b, tempPt;

    init_mat_mp(A, num_vars, num_vars); 
    init_vec_mp(b, num_vars); 
    init_vec_mp(tempPt, num_vars);
    A->rows = A->cols = b->size = tempPt->size = num_vars;

    if (regen->level[0].useIntrinsicSlice)
    { // setup B_transpose
      init_mat_mp(B_transpose, regen->level[0].B_mp->cols, regen->level[0].B_mp->rows);
      transpose_mp(B_transpose, regen->level[0].B_mp);
    }

    // the bottom of A & b is fixed - main slices and patches
    for (j = depth; j < num_funcs; j++)
    { // setup main slices
      set_zero_mp(&b->coord[j]); // rhs
      for (k = 0; k < num_vars; k++)
        set_mp(&A->entry[j][k], &regen->mainCoeff_mp->entry[j][k]);
    }
    for (j = 0; j < regen->num_var_gps; j++)
    { // copy patches
      i = j + num_funcs;
      set_zero_mp(&b->coord[i]); // rhs
      for (k = 0; k < num_vars; k++)
        set_mp(&A->entry[i][k], &regen->patchCoeff_mp->entry[j][k]);
    }
    // new hom patch
    i = num_funcs + regen->num_var_gps;
    set_one_mp(&b->coord[i]);
    for (j = 0; j < num_vars; j++)
      set_mp(&A->entry[i][j], &regen->main_homVar_mp->coord[j]);

    // loop over the decompositions
    for (i = 0; i < count; i++)
    { // compute the start points related to this decomposition
      actualCount += createRegenStartPtsDecomp_mp(regen, decomp[i], depth, A, b, tempPt, B_transpose, START);
    }

    clear_mat_mp(A);
    clear_vec_mp(b); clear_vec_mp(tempPt);
    if (regen->level[0].useIntrinsicSlice)
      clear_mat_mp(B_transpose);
  }

  return actualCount;
}

void printPointLinearDegree(FILE *FP, point_d Pt_d, point_mp Pt_mp, int prec, int *linear, int *degree, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the point, linear & degree to FP                *
\***************************************************************/
{
  int i;

  // print Pt
  if (prec < 64)
  { // print Pt_d
    for (i = 0; i < Pt_d->size; i++)
      fprintf(FP, "%.15e %.15e;\n", Pt_d->coord[i].r, Pt_d->coord[i].i);
  }
  else
  { // print Pt_mp
    for (i = 0; i < Pt_mp->size; i++)
    {
      print_mp(FP, 0, &Pt_mp->coord[i]);
      fprintf(FP, ";\n");
    }
  }

  // print linear & degree
  for (i = 0; i < size; i++)
    fprintf(FP, "%d %d\n", linear[i], degree[i]);
  fprintf(FP, "\n");

  return;
}

void printRegenSummaryData(regen_t *regen, int finished_level, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the summary data to FP                          *
\***************************************************************/
{
  int i;

  // print the num_levels, num_funcs, num_variables, num_var_gps
  fprintf(FP, "%d %d %d %d\n", regen->num_levels, regen->num_funcs, regen->num_variables, regen->num_var_gps);

  // print the level just finished
  fprintf(FP, "%d\n", finished_level);

  // print the info about each level
  for (i = 0; i < regen->num_levels; i++)
    if (i <= finished_level)
    { // print info about how the level was
      fprintf(FP, "%d %d %d %d %d %d %d %d\n", regen->level[i].level, regen->level[i].depth, regen->level[i].num_paths, regen->level[i].num_sing, regen->level[i].num_nonsing, regen->level[i].num_inf, regen->level[i].num_higher_dim, regen->level[i].num_bad);
    }
    else
    { // print info about future levels
      fprintf(FP, "%d %d %d\n", regen->level[i].level, regen->level[i].depth, regen->level[i].useIntrinsicSlice);
    }
  fprintf(FP, "\n");

  return;
}

void printRegenRelevantData(regen_t *regen, int MPType, int curr_level, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the relevant data to FP so that we can begin    *
* exactly right here in case of failure or change tolerances    *
\***************************************************************/
{
  int i, j, k, l, rows, cols;

  // print an X signifying that this is extra info
  fprintf(FP, "X\n");

  // print the MPType, num_levels, num_funcs, num_variables, num_var_gps
  fprintf(FP, "%d %d %d %d %d\n", MPType, regen->num_levels, regen->num_funcs, regen->num_variables, regen->num_var_gps);

  // print the current level
  fprintf(FP, "%d\n", curr_level);

  // print the info about each level
  for (i = 0; i < regen->num_levels; i++)
    if (i < curr_level)
    { // print info about how the level was
      fprintf(FP, "%d %d %d %d %d %d %d %d\n", regen->level[i].level, regen->level[i].depth, regen->level[i].num_paths, regen->level[i].num_sing, regen->level[i].num_nonsing, regen->level[i].num_inf, regen->level[i].num_higher_dim, regen->level[i].num_bad);
    }
    else if (i == curr_level)
    { // print info about the current level
      fprintf(FP, "%d %d %d %d\n", regen->level[i].level, regen->level[i].depth, regen->level[i].num_paths, regen->level[i].useIntrinsicSlice);
    }
    else
    { // print info about future levels
      fprintf(FP, "%d %d %d\n", regen->level[i].level, regen->level[i].depth, regen->level[i].useIntrinsicSlice);
    }
  fprintf(FP, "\n");

  // print the patch
  if (MPType == 0)
  { // print _d
    rows = regen->patchCoeff_d->rows;
    cols = regen->patchCoeff_d->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_d(FP, 0, &regen->patchCoeff_d->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }
  else if (MPType == 1)
  { // print _mp
    rows = regen->patchCoeff_mp->rows;
    cols = regen->patchCoeff_mp->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_mp(FP, 0, &regen->patchCoeff_mp->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }
  else
  { // print _rat
    rows = regen->patchCoeff_d->rows;
    cols = regen->patchCoeff_d->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(FP, 10, regen->patchCoeff_rat[i][j][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, regen->patchCoeff_rat[i][j][1]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }

  // print the main coefficients
  if (MPType == 0)
  { // print _d
    rows = regen->mainCoeff_d->rows;
    cols = regen->mainCoeff_d->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_d(FP, 0, &regen->mainCoeff_d->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }
  else if (MPType == 1)
  { // print _mp
    rows = regen->mainCoeff_mp->rows;
    cols = regen->mainCoeff_mp->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        print_mp(FP, 0, &regen->mainCoeff_mp->entry[i][j]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }
  else
  { // print _rat
    rows = regen->mainCoeff_d->rows;
    cols = regen->mainCoeff_d->cols;
    fprintf(FP, "%d %d\n", rows, cols);
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(FP, 10, regen->mainCoeff_rat[i][j][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, regen->mainCoeff_rat[i][j][1]);
        fprintf(FP, "\n");
      }
    fprintf(FP, "\n");
  }

  // print the coefficients
  if (MPType == 0)
  { // print _d
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < regen->num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
          for (l = 0; l < regen->num_variables; l++)
          {
            print_d(FP, 0, regen->coeff_d[i][j][k][l]);
            fprintf(FP, "\n");
          }
    fprintf(FP, "\n");
  }
  else if (MPType == 1)
  { // print _mp
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < regen->num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
          for (l = 0; l < regen->num_variables; l++)
          {
            print_mp(FP, 0, regen->coeff_mp[i][j][k][l]);
            fprintf(FP, "\n");
          }
    fprintf(FP, "\n");
  }
  else
  { // print _rat
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < regen->num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
          for (l = 0; l < regen->num_variables; l++)
          {
            mpq_out_str(FP, 10, regen->coeff_rat[i][j][k][l][0]);
            fprintf(FP, " ");
            mpq_out_str(FP, 10, regen->coeff_rat[i][j][k][l][1]);
            fprintf(FP, "\n");
          }
    fprintf(FP, "\n");
  }

  // print main_homVar
  if (MPType == 0)
  { // print _d
    rows = regen->main_homVar_d->size;
    fprintf(FP, "%d\n", rows);
    for (i = 0; i < rows; i++)
    {
      print_d(FP, 0, &regen->main_homVar_d->coord[i]);
      fprintf(FP, "\n");
    }
    fprintf(FP, "\n");
  }
  else if (MPType == 1)
  { // print _mp
    rows = regen->main_homVar_mp->size;
    fprintf(FP, "%d\n", rows);
    for (i = 0; i < rows; i++)
    {
      print_mp(FP, 0, &regen->main_homVar_mp->coord[i]);
      fprintf(FP, "\n");
    }
    fprintf(FP, "\n");
  }
  else
  { // print _rat
    rows = regen->main_homVar_d->size;
    fprintf(FP, "%d\n", rows);
    for (i = 0; i < rows; i++)
    {
      mpq_out_str(FP, 10, regen->main_homVar_rat[i][0]);
      fprintf(FP, " ");
      mpq_out_str(FP, 10, regen->main_homVar_rat[i][1]);
      fprintf(FP, "\n");
    }
    fprintf(FP, "\n");
  }

  // print the intrinsic slice, if needed
  if (regen->level[curr_level].useIntrinsicSlice)
  { // print B
    if (MPType == 0)
    { // print _d
      rows = regen->level[curr_level].B_d->rows;
      cols = regen->level[curr_level].B_d->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          print_d(FP, 0, &regen->level[curr_level].B_d->entry[i][j]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
    }
    else if (MPType == 1)
    { // print _mp
      rows = regen->level[curr_level].B_mp->rows;
      cols = regen->level[curr_level].B_mp->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          print_mp(FP, 0, &regen->level[curr_level].B_mp->entry[i][j]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
    } 
    else 
    { // print _rat
      rows = regen->level[curr_level].B_d->rows;
      cols = regen->level[curr_level].B_d->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_out_str(FP, 10, regen->level[curr_level].B_rat[i][j][0]);
          fprintf(FP, " ");
          mpq_out_str(FP, 10, regen->level[curr_level].B_rat[i][j][1]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
    } 

    // print p
    if (MPType == 0)
    { // print _d
      rows = regen->level[curr_level].p_d->size;
      fprintf(FP, "%d\n", rows);
      for (i = 0; i < rows; i++)
      {
        print_d(FP, 0, &regen->level[curr_level].p_d->coord[i]);
        fprintf(FP, "\n");
      }
      fprintf(FP, "\n");
    }
    else if (MPType == 1)
    { // print _mp
      rows = regen->level[curr_level].p_mp->size;
      fprintf(FP, "%d\n", rows);
      for (i = 0; i < rows; i++)
      {
        print_mp(FP, 0, &regen->level[curr_level].p_mp->coord[i]);
        fprintf(FP, "\n");
      }
      fprintf(FP, "\n");
    }
    else
    { // print _rat
      rows = regen->level[curr_level].p_d->size;
      fprintf(FP, "%d\n", rows);
      for (i = 0; i < rows; i++)
      {
        mpq_out_str(FP, 10, regen->level[curr_level].p_rat[i][0]);
        fprintf(FP, " ");
        mpq_out_str(FP, 10, regen->level[curr_level].p_rat[i][1]);
        fprintf(FP, "\n");
      }
      fprintf(FP, "\n");
    }
  }
  fprintf(FP, "\n");

  // print the squaring matrix, if needed
  if (!regen->noChanges)
  { // print the matrix A
    if (MPType == 0)
    { // print _d
      square_system_eval_data_d *SSED = (square_system_eval_data_d *)regen->square_d;
      rows = SSED->A->rows;
      cols = SSED->A->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          print_d(FP, 0, &SSED->A->entry[i][j]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
      SSED = NULL;
    }
    else if (MPType == 1)
    { // print _mp
      square_system_eval_data_mp *SSED = (square_system_eval_data_mp *)regen->square_mp;
      rows = SSED->A->rows;
      cols = SSED->A->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          print_mp(FP, 0, &SSED->A->entry[i][j]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
    }
    else
    { // print _rat
      square_system_eval_data_mp *SSED = (square_system_eval_data_mp *)regen->square_mp;
      rows = SSED->A->rows;
      cols = SSED->A->cols;
      fprintf(FP, "%d %d\n", rows, cols);
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_out_str(FP, 10, SSED->A_rat[i][j][0]);
          fprintf(FP, " ");
          mpq_out_str(FP, 10, SSED->A_rat[i][j][1]);
          fprintf(FP, "\n");
        }
      fprintf(FP, "\n");
    }
  }

  return;
}

int change_regen_prec(void const *RED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for regeneration                      *
\***************************************************************/
{
  int i, j, k, l;

  regen_t *regen = (regen_t *)RED;

  // change the precision on the square system
  regen->change_square_prec(regen->square_mp, prec);

  // change the precision on the other structures, if needed
  if (regen->curr_precision != prec)
  { // change precision on patchCoeff
    change_prec_mat_mp_rat(regen->patchCoeff_mp, prec, regen->patchCoeff_rat);

    // change precision on mainCoeff
    change_prec_mat_mp_rat(regen->mainCoeff_mp, prec, regen->mainCoeff_rat);

    // change precision on main_homVar
    change_prec_vec_mp_rat(regen->main_homVar_mp, prec, regen->main_homVar_rat);

    // change the precision on gamma
    change_prec_mp_rat(regen->gamma_mp, prec, regen->gamma_rat);

    // change the precision on coeff
    for (i = 0; i < regen->num_funcs; i++)
      for (j = 0; j < regen->num_var_gps; j++)
        for (k = 0; k < regen->degrees[i][j]; k++)
          for (l = 0; l < regen->num_variables; l++)
          {
            change_prec_mp_rat(regen->coeff_mp[i][j][k][l], prec, regen->coeff_rat[i][j][k][l]);
          }

    // change the precision for each the current level
    k = regen->curr_level_num;
    if (0 <= k && k < regen->num_levels)
    { // see if the B & p are used
      if (regen->level[k].useIntrinsicSlice)
      { // change the precision on B
        if (regen->level[k].B_rat != NULL)
          change_prec_mat_mp_rat(regen->level[k].B_mp, prec, regen->level[k].B_rat);
 
        // change the precision on p
        if (regen->level[k].p_rat != NULL)
          change_prec_vec_mp_rat(regen->level[k].p_mp, prec, regen->level[k].p_rat);
      }
    }
  }
  regen->curr_precision = prec;

  return 0;
}

void copyRegenLevelStructures(regen_t *regenIn, int max_threads, regen_t *regen_copy, int MPType, int level_num)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy the structures for level_num of regen             *
\***************************************************************/
{
  int i, j, k, r, c;

  for (i = 0; i < max_threads; i++)
  { // initilialize values
    regen_copy[i].level[level_num].level = regenIn->level[level_num].level;
    regen_copy[i].level[level_num].depth = regenIn->level[level_num].depth;
    regen_copy[i].level[level_num].num_paths = regenIn->level[level_num].num_paths;
    regen_copy[i].level[level_num].num_sing = regenIn->level[level_num].num_sing;
    regen_copy[i].level[level_num].num_nonsing = regenIn->level[level_num].num_nonsing;
    regen_copy[i].level[level_num].num_inf = regenIn->level[level_num].num_inf;
    regen_copy[i].level[level_num].num_higher_dim = regenIn->level[level_num].num_higher_dim;
    regen_copy[i].level[level_num].num_bad = regenIn->level[level_num].num_bad;

    regen_copy[i].level[level_num].useIntrinsicSlice = regenIn->level[level_num].useIntrinsicSlice;

    // setup B & p, if needed
    if (regen_copy[i].level[level_num].useIntrinsicSlice)
    { // setup needed
      if (MPType == 0)
      { // setup B_d & p_d
        r = regenIn->level[level_num].B_d->rows;
        c = regenIn->level[level_num].B_d->cols;
        init_mat_d(regen_copy[i].level[level_num].B_d, r, c);
        init_vec_d(regen_copy[i].level[level_num].p_d, regenIn->level[level_num].p_d->size);

        mat_cp_d(regen_copy[i].level[level_num].B_d, regenIn->level[level_num].B_d);
        vec_cp_d(regen_copy[i].level[level_num].p_d, regenIn->level[level_num].p_d);
      }
      else if (MPType == 1)
      { // setup B_mp & p_mp
        r = regenIn->level[level_num].B_mp->rows;
        c = regenIn->level[level_num].B_mp->cols;
        init_mat_mp(regen_copy[i].level[level_num].B_mp, r, c);
        init_vec_mp(regen_copy[i].level[level_num].p_mp, regenIn->level[level_num].p_mp->size);

        mat_cp_mp(regen_copy[i].level[level_num].B_mp, regenIn->level[level_num].B_mp);
        vec_cp_mp(regen_copy[i].level[level_num].p_mp, regenIn->level[level_num].p_mp);
      }
      else
      { // setup B
        r = regenIn->level[level_num].B_d->rows;
        c = regenIn->level[level_num].B_d->cols; 
        init_mat_d(regen_copy[i].level[level_num].B_d, r, c);
        init_mat_mp2(regen_copy[i].level[level_num].B_mp, r, c, regen_copy[i].curr_precision);
        // point to B_rat
        regen_copy[i].level[level_num].B_rat = regenIn->level[level_num].B_rat;
        // setup B_d & B_mp
        regen_copy[i].level[level_num].B_d->rows = regen_copy[i].level[level_num].B_mp->rows = r;
        regen_copy[i].level[level_num].B_d->cols = regen_copy[i].level[level_num].B_mp->cols = c;
        for (j = 0; j < r; j++)
          for (k = 0; k < c; k++)
          {
            mpf_set_q(regen_copy[i].level[level_num].B_mp->entry[j][k].r, regen_copy[i].level[level_num].B_rat[j][k][0]);
            mpf_set_q(regen_copy[i].level[level_num].B_mp->entry[j][k].i, regen_copy[i].level[level_num].B_rat[j][k][1]);
            regen_copy[i].level[level_num].B_d->entry[j][k].r = mpq_get_d(regen_copy[i].level[level_num].B_rat[j][k][0]);
            regen_copy[i].level[level_num].B_d->entry[j][k].i = mpq_get_d(regen_copy[i].level[level_num].B_rat[j][k][1]);
          }

        // setup p
        r = regenIn->level[level_num].p_d->size;
        init_vec_d(regen_copy[i].level[level_num].p_d, r);
        init_vec_mp2(regen_copy[i].level[level_num].p_mp, r, regen_copy[i].curr_precision);
        // point to p_rat
        regen_copy[i].level[level_num].p_rat = regenIn->level[level_num].p_rat;
        // setup p_d & p_mp
        regen_copy[i].level[level_num].p_d->size = regen_copy[i].level[level_num].p_mp->size = r;

        for (j = 0; j < r; j++)
        {
          mpf_set_q(regen_copy[i].level[level_num].p_mp->coord[j].r, regen_copy[i].level[level_num].p_rat[j][0]);
          mpf_set_q(regen_copy[i].level[level_num].p_mp->coord[j].i, regen_copy[i].level[level_num].p_rat[j][1]);
          regen_copy[i].level[level_num].p_d->coord[j].r = mpq_get_d(regen_copy[i].level[level_num].p_rat[j][0]);
          regen_copy[i].level[level_num].p_d->coord[j].i = mpq_get_d(regen_copy[i].level[level_num].p_rat[j][1]);
        }
      }
    }
  }

  return;
}

void clearRegenLevelStructures(int max_threads, regen_t *regen, int MPType, int level_num)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the structures for level_num of regen            *
\***************************************************************/
{
  int i;

  for (i = 0; i < max_threads; i++)
  {
    if (regen[i].level[level_num].useIntrinsicSlice)
    { // need to clear B & p
      if (MPType == 0)
      { // clear B_d & p_d
        clear_mat_d(regen[i].level[level_num].B_d);
        clear_vec_d(regen[i].level[level_num].p_d);
      }
      else if (MPType == 1)
      { // clear B_mp & p_mp
        clear_mat_mp(regen[i].level[level_num].B_mp);
        clear_vec_mp(regen[i].level[level_num].p_mp);
      }
      else
      { // clear B & p
        if (i == 0)
        { // the rationals only need cleared for the first one since all others point to these
          clear_mat_rat(regen[i].level[level_num].B_rat, regen[i].level[level_num].B_d->rows, regen[i].level[level_num].B_d->cols);
          clear_vec_rat(regen[i].level[level_num].p_rat, regen[i].level[level_num].p_d->size);
        }
        clear_mat_d(regen[i].level[level_num].B_d);
        clear_vec_d(regen[i].level[level_num].p_d);
        clear_mat_mp(regen[i].level[level_num].B_mp);
        clear_vec_mp(regen[i].level[level_num].p_mp);
      }
    }

    regen[i].level[level_num].B_d->rows = regen[i].level[level_num].B_d->cols = regen[i].level[level_num].B_mp->rows = regen[i].level[level_num].B_mp->cols = 0;
    regen[i].level[level_num].B_rat = NULL;
    regen[i].level[level_num].p_d->size = regen[i].level[level_num].p_mp->size = 0;
    regen[i].level[level_num].p_rat = NULL;
  }

  return;
}

void copyRegenRandom(regen_t *regen, regen_t *regenIn, int MPType, void const *square_d, void const *square_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: copy the random numbers to the other copies of regen   *
\***************************************************************/
{
  int j, k, l, m;

  // initialize the values
  regen->num_levels = regenIn->num_levels;
  regen->num_funcs = regenIn->num_funcs;
  regen->num_variables = regenIn->num_variables;
  regen->num_var_gps = regenIn->num_var_gps;

  // setup PPD
  regen->PPD.num_funcs = regenIn->PPD.num_funcs;
  regen->PPD.num_hom_var_gp = regenIn->PPD.num_hom_var_gp;
  regen->PPD.num_var_gp = regenIn->PPD.num_var_gp;
  regen->PPD.type = regenIn->PPD.type;
  regen->PPD.size = regenIn->PPD.size;

  // setup 'square'
  regen->square_d = square_d;
  regen->square_mp = square_mp;
  regen->square_eval_d = regenIn->square_eval_d;
  regen->square_eval_mp = regenIn->square_eval_mp;
  regen->change_square_prec = regenIn->change_square_prec;

  // setup noChanges & inst counts
  regen->noChanges = regenIn->noChanges;
  regen->numSubFuncs = regenIn->numSubFuncs;
  regen->startSub = regenIn->startSub;
  regen->endSub = regenIn->endSub;
  regen->startFunc = regenIn->startFunc;
  regen->endFunc = regenIn->endFunc;
  regen->startJvsub = regenIn->startJvsub;
  regen->endJvsub = regenIn->endJvsub;
  regen->startJv = regenIn->startJv;
  regen->endJv = regenIn->endJv;
  regen->subFuncsBelow = regenIn->subFuncsBelow;

  // setup patchCoeff
  if (MPType == 0 || MPType == 2)
  { // setup _d
    init_mat_d(regen->patchCoeff_d, regenIn->patchCoeff_d->rows, regenIn->patchCoeff_d->cols);
    mat_cp_d(regen->patchCoeff_d, regenIn->patchCoeff_d);
  }
  if (MPType == 1 || MPType == 2)
  { // setup _mp
    init_mat_mp2(regen->patchCoeff_mp, regenIn->patchCoeff_mp->rows, regenIn->patchCoeff_mp->cols, regenIn->curr_precision);
    mat_cp_mp(regen->patchCoeff_mp, regenIn->patchCoeff_mp);
  }
  if (MPType == 2)
  { // setup _rat
    regen->patchCoeff_rat = regenIn->patchCoeff_rat;
  }

  // setup degrees
  regen->degrees = regenIn->degrees;

  // setup gamma
  if (MPType == 0 || MPType == 2)
  { // copy gamma_d
    set_d(regen->gamma_d, regenIn->gamma_d);
  }
  if (MPType == 1 || MPType == 2)
  { // copy gamma_mp
    init_mp2(regen->gamma_mp, regenIn->curr_precision);
    set_mp(regen->gamma_mp, regenIn->gamma_mp);
  }
  if (MPType == 2)
  { // setup gamma_rat
    init_rat(regen->gamma_rat);
    set_rat(regen->gamma_rat, regenIn->gamma_rat);
  }

  // setup coeff
  if (MPType == 0)
  { // point to coeff_d since the precision is not changing (fixed in double)
    regen->coeff_d = regenIn->coeff_d;
  }
  else if (MPType == 1)
  { // point to coeff_mp since the precision is not changing (fixed mulit precision)
    regen->coeff_mp = regenIn->coeff_mp;
  }
  else
  { // point to coeff_d since the precision is not changing (fixed in double)
    regen->coeff_d = regenIn->coeff_d;
    // point to coeff_rat
    regen->coeff_rat = regenIn->coeff_rat;
    // setup coeff_mp since we can change the precision
    regen->coeff_mp = (comp_mp ****)bmalloc(regenIn->num_funcs * sizeof(comp_mp ***));
    for (j = 0; j < regenIn->num_funcs; j++)
    {
      regen->coeff_mp[j] = (comp_mp ***)bmalloc(regenIn->num_var_gps * sizeof(comp_mp **));
      for (k = 0; k < regenIn->num_var_gps; k++)
      {
        regen->coeff_mp[j][k] = (comp_mp **)bmalloc(regenIn->degrees[j][k] * sizeof(comp_mp *));
        for (l = 0; l < regenIn->degrees[j][k]; l++)
        {
          regen->coeff_mp[j][k][l] = (comp_mp *)bmalloc((regenIn->num_variables + 1) * sizeof(comp_mp));
          for (m = 0; m <= regenIn->num_variables; m++)
          {
            init_mp2(regen->coeff_mp[j][k][l][m], regenIn->curr_precision);
            set_mp(regen->coeff_mp[j][k][l][m], regenIn->coeff_mp[j][k][l][m]);
          }
        }
      }
    }
  }

  // setup other values
  regen->curr_precision = regenIn->curr_precision;
  regen->curr_level_num = regenIn->curr_level_num;
  regen->curr_linear = regen->curr_linear_degree = NULL;

  // allocate for the levels
  regen->level = (regenLevel_t *)bmalloc(regen->num_levels * sizeof(regenLevel_t));

  return;
}

void clearRegenRandom_zero_dim(int max, regen_t *regen, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear the random number structures in regen            *
\***************************************************************/
{
  int i, j, k, l, m;

  for (i = max - 1; i >= 0; i--)
  { // clear regen[i]

    // clear PPD
    if (i == 0)
    { // release memory
      free(regen[i].PPD.type);
      free(regen[i].PPD.size);
    }
    // NULL out pointers
    regen[i].PPD.type = NULL;
    regen[i].PPD.size = NULL;

    // clear 'square' for the copies - do not clear the original, just NULL it out
    if (i > 0)
    {
      if (MPType == 0 || MPType == 2)
      { // clear square_d
        square_system_eval_data_clear_d((square_system_eval_data_d *)regen[i].square_d, MPType);
      }
     if (MPType == 1 || MPType == 2)
     { // clear square_mp
       square_system_eval_data_clear_mp((square_system_eval_data_mp *)regen[i].square_mp, MPType == 1);
      }
    }
    regen[i].square_d = NULL;
    regen[i].square_mp = NULL;
    regen[i].square_eval_d = NULL;
    regen[i].square_eval_mp = NULL;
    regen[i].change_square_prec = NULL;

    // clear instCount
    if (i == 0)
    {
      free(regen[i].startSub);
      free(regen[i].endSub);
      free(regen[i].startFunc);
      free(regen[i].endFunc);
      free(regen[i].startJvsub);
      free(regen[i].endJvsub);
      free(regen[i].startJv);
      free(regen[i].endJv);

      if (regen[i].numSubFuncs > 0)
      { // clear subFuncsBelow
        for (j = regen[i].num_funcs - 1; j >= 0; j--)
          free(regen[i].subFuncsBelow[j]);
        free(regen[i].subFuncsBelow);
      }
    }
    regen[i].startSub = regen[i].endSub = regen[i].startFunc = regen[i].endFunc = regen[i].startJvsub = regen[i].endJvsub = regen[i].startJv = regen[i].endJv = NULL;
    regen[i].subFuncsBelow = NULL;

    // clear patchCoeff
    if (i == 0 && MPType == 2)
    { // clear _rat first
      clear_mat_rat(regen[i].patchCoeff_rat, regen[i].patchCoeff_d->rows, regen[i].patchCoeff_d->cols);
    }
    regen[i].patchCoeff_rat = NULL;
    if (MPType == 0 || MPType == 2)
    { // clear _d
      clear_mat_d(regen[i].patchCoeff_d);
    }
    if (MPType == 1 || MPType == 2)
    { // clear _mp
      clear_mat_mp(regen[i].patchCoeff_mp);
    }

    // clear mainCoeff, main_homVar
    if (i == 0 && MPType == 2)
    { // clear _rat
      clear_mat_rat(regen[i].mainCoeff_rat, regen[i].mainCoeff_d->rows, regen[i].mainCoeff_d->cols);
      clear_vec_rat(regen[i].main_homVar_rat, regen[i].main_homVar_d->size);
    }
    regen[i].mainCoeff_rat = NULL;
    regen[i].main_homVar_rat = NULL;
    if (MPType == 0 || MPType == 2)
    { // clear _d
      clear_mat_d(regen[i].mainCoeff_d);
      clear_vec_d(regen[i].main_homVar_d);
    }
    if (MPType == 1 || MPType == 2)
    { // clear _mp
      clear_mat_mp(regen[i].mainCoeff_mp);
      clear_vec_mp(regen[i].main_homVar_mp);
    }

    // clear gamma
    if (MPType == 0 || MPType == 2)
    { // clear gamma_d
      clear_d(regen[i].gamma_d);
    }
    if (MPType == 1 || MPType == 2)
    { // clear gamma_mp
      clear_mp(regen[i].gamma_mp);
    }
    if (MPType == 2)
    { // clear gamma_rat
      mpq_clear(regen[i].gamma_rat[0]); 
      mpq_clear(regen[i].gamma_rat[1]);
    }

    // clear coeff
    if ((i > 0 && MPType == 2) || (i == 0 && (MPType == 1 || MPType == 2)))
    { // clear coeff_mp
      for (j = regen[i].num_funcs - 1; j >= 0; j--)
      {
        for (k = regen[i].num_var_gps - 1; k >= 0; k--)
        {
          for (l = regen[i].degrees[j][k] - 1; l >= 0; l--)
          {
            for (m = regen[i].num_variables - 1; m >= 0; m--)
            {
              clear_mp(regen[i].coeff_mp[j][k][l][m]);
            }
            free(regen[i].coeff_mp[j][k][l]);
          }
          free(regen[i].coeff_mp[j][k]);
        }
        free(regen[i].coeff_mp[j]);
      }
      free(regen[i].coeff_mp);
    }
    if (i == 0 && MPType == 2)
    { // clear coeff_rat
      for (j = regen[i].num_funcs - 1; j >= 0; j--)
      {
        for (k = regen[i].num_var_gps - 1; k >= 0; k--)
        {
          for (l = regen[i].degrees[j][k] - 1; l >= 0; l--)
          {
            for (m = regen[i].num_variables - 1; m >= 0; m--)
            {
              clear_rat(regen[i].coeff_rat[j][k][l][m]);
              free(regen[i].coeff_rat[j][k][l][m]);
            }
            free(regen[i].coeff_rat[j][k][l]);
          }
          free(regen[i].coeff_rat[j][k]);
        }
        free(regen[i].coeff_rat[j]);
      }
      free(regen[i].coeff_rat);
    }
    if (i == 0 && (MPType == 0 || MPType == 2))
    { // clear coeff_d
      for (j = regen[i].num_funcs - 1; j >= 0; j--)
      {
        for (k = regen[i].num_var_gps - 1; k >= 0; k--)
        {
          for (l = regen[i].degrees[j][k] - 1; l >= 0; l--)
          {
            free(regen[i].coeff_d[j][k][l]);
          }
          free(regen[i].coeff_d[j][k]);
        }
        free(regen[i].coeff_d[j]);
      }
      free(regen[i].coeff_d);
    }
    regen[i].coeff_d = NULL;
    regen[i].coeff_mp = NULL;
    regen[i].coeff_rat = NULL;

    // clear degrees
    if (i == 0)
    { // release memory
      free(regen[i].degrees);
    }
    // NULL out pointer
    regen[i].degrees = NULL;

    // NULL out curr_linear & curr_linear_degree
    regen[i].curr_linear = regen[i].curr_linear_degree = NULL;

    // clear level
    free(regen[i].level);
    regen[i].level = NULL;
  }

  return;
}

void setupRegenRestart(regen_t *regen, tracker_config_t *T, char *startName, int curr_level)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the regeneration using the data in 'startName'   *
\***************************************************************/
{
  int i, j, k, l, rows, cols, inputMPType = -1, num_levels = -1, num_funcs = -1, num_variables = -1, num_var_gps = -1, inputCurrLevel = -1;
  char ch;
  comp_d ****coeff_d = NULL;
  comp_mp ****coeff_mp = NULL;
  mpq_t *****coeff_rat = NULL;
  vec_d p_d;
  vec_mp p_mp;
  mpq_t **p_rat = NULL;
  mat_d patch_d;
  mat_mp patch_mp;
  mpq_t ***patch_rat = NULL;
  FILE *INPUT = fopen(startName, "r");
  if (INPUT == NULL)
  {
    printf("ERROR: The file used to restart the regeneration, '%s', does not exist!\n", startName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // initialize patch_d & patch_mp, p_d & p_mp
  init_mat_d(patch_d, 0, 0);
  init_mat_mp(patch_mp, 0, 0);
  init_vec_d(p_d, 0);
  init_vec_mp(p_mp, 0);

  // move past all the other data and find the 'X'
  do
  {
    ch = fgetc(INPUT);
  } while (ch != 'X');


  // read in the MPType, num_levels, num_funcs, num_variables, num_var_gps
  fscanf(INPUT, "%d %d %d %d %d\n", &inputMPType, &num_levels, &num_funcs, &num_variables, &num_var_gps);
  // read in the current level
  fscanf(INPUT, "%d\n", &inputCurrLevel);

  // error checking
  if (num_funcs != regen->num_funcs)
  {
    printf("ERROR: The number of functions (%d vs %d) is not correct!\n", regen->num_funcs, num_funcs);
    bexit(ERROR_CONFIGURATION);
  }
  if (num_levels <= 0 || num_levels > num_funcs)
  {
    printf("ERROR: The number of levels (%d) must be > 0 and <= %d!\n", num_levels, num_funcs);
    bexit(ERROR_CONFIGURATION);
  }
  if (num_variables != regen->num_variables)
  {
    printf("ERROR: The number of variables (%d vs %d) is not correct!\n", regen->num_variables, num_variables);
    bexit(ERROR_CONFIGURATION);
  }
  if (num_var_gps != regen->num_var_gps)
  {
    printf("ERROR: The number of variable groups (%d vs %d) is not correct!\n", regen->num_var_gps, num_var_gps);
    bexit(ERROR_CONFIGURATION);
  }
  if (inputCurrLevel != curr_level)
  {
    printf("ERROR: The current level number (%d vs %d) is not correct!\n", inputCurrLevel, curr_level);
    bexit(ERROR_CONFIGURATION);
  }
  if (curr_level < 0 || curr_level >= num_levels)
  {
    printf("ERROR: The current level (%d) must be >= 0 and < %d!\n", curr_level, num_levels);
    bexit(ERROR_CONFIGURATION);
  }

  // so we can assume that we probably have the correct file

  // setup the number of levels (other data is already known)
  regen->num_levels = num_levels; 

  // allocate level
  regen->level = (regenLevel_t *)bmalloc(num_levels * sizeof(regenLevel_t));

  // read in the info about each level
  for (i = 0; i < num_levels; i++)
    if (i < curr_level)
    { // read info about how this level was
      fscanf(INPUT, "%d %d %d %d %d %d %d %d\n", &regen->level[i].level, &regen->level[i].depth, &regen->level[i].num_paths, &regen->level[i].num_sing, &regen->level[i].num_nonsing, &regen->level[i].num_inf, &regen->level[i].num_higher_dim, &regen->level[i].num_bad);
    }
    else if (i == curr_level)
    { // read info about the current level
      fscanf(INPUT, "%d %d %d %d\n", &regen->level[i].level, &regen->level[i].depth, &regen->level[i].num_paths, &regen->level[i].useIntrinsicSlice);
    }
    else
    { // read info about future levels
      fscanf(INPUT, "%d %d %d\n", &regen->level[i].level, &regen->level[i].depth, &regen->level[i].useIntrinsicSlice);
    }
  fscanf(INPUT, "\n");

  // read in the patch
  if (inputMPType == 0)
  { // use _d
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_d(patch_d, rows, cols);
    patch_d->rows = rows;
    patch_d->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        fscanf(INPUT, "%lf %lf\n", &patch_d->entry[i][j].r, &patch_d->entry[i][j].i);
      }
    fscanf(INPUT, "\n");
  }
  else if (inputMPType == 1)
  { // use _mp
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_mp(patch_mp, rows, cols);
    patch_mp->rows = rows;
    patch_mp->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpf_inp_str(patch_mp->entry[i][j].r, INPUT, 10);
        mpf_inp_str(patch_mp->entry[i][j].i, INPUT, 10);
        fscanf(INPUT, "\n");
      }
    fscanf(INPUT, "\n");
  }    
  else
  { // use _rat
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_d(patch_d, rows, cols);
    change_size_mat_mp(patch_mp, rows, cols);
    init_mat_rat(patch_rat, rows, cols);
    patch_d->rows = patch_mp->rows = rows;
    patch_d->cols = patch_mp->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_inp_str(patch_rat[i][j][0], INPUT, 10);
        mpq_canonicalize(patch_rat[i][j][0]);
        mpq_inp_str(patch_rat[i][j][1], INPUT, 10);
        mpq_canonicalize(patch_rat[i][j][1]);
        fscanf(INPUT, "\n");
      }
    fscanf(INPUT, "\n");
  }

  // setup patch
  if (T->MPType == 0)
  { // update patchCoeff_d
    if (inputMPType == 0)
    { // copy _d
      mat_cp_d(regen->patchCoeff_d, patch_d);
    }
    else if (inputMPType == 1)
    { // copy _mp to _d
      mat_mp_to_d(regen->patchCoeff_d, patch_mp);
    }
    else
    { // copy _rat to _d
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          regen->patchCoeff_d->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
          regen->patchCoeff_d->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  }
  else if (T->MPType == 1)
  { // update patchCoeff_mp
    if (inputMPType == 0)
    { // copy _d to _mp
      mat_d_to_mp(regen->patchCoeff_mp, patch_d);
    }
    else if (inputMPType == 1)
    { // copy _mp
      mat_cp_mp(regen->patchCoeff_mp, patch_mp);
    }
    else
    { // copy _rat to _mp
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_set_q(regen->patchCoeff_mp->entry[i][j].r, patch_rat[i][j][0]);
          mpf_set_q(regen->patchCoeff_mp->entry[i][j].i, patch_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  }
  else
  { // update _d, _mp & _rat
    if (inputMPType == 0)
    { // copy _d to _rat
      comp_mp tempComp;
      init_mp(tempComp);
   
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          comp_d_to_mp_rat(tempComp, regen->patchCoeff_rat[i][j], &patch_d->entry[i][j], 16, regen->curr_precision, 0, 0);

          mpf_set_q(regen->patchCoeff_mp->entry[i][j].r, regen->patchCoeff_rat[i][j][0]);
          mpf_set_q(regen->patchCoeff_mp->entry[i][j].i, regen->patchCoeff_rat[i][j][1]);
          regen->patchCoeff_d->entry[i][j].r = mpq_get_d(regen->patchCoeff_rat[i][j][0]);
          regen->patchCoeff_d->entry[i][j].i = mpq_get_d(regen->patchCoeff_rat[i][j][1]);
        }
      clear_mp(tempComp);
    }
    else if (inputMPType == 1)
    { // copy _mp to _rat
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_t_to_rat(regen->patchCoeff_rat[i][j][0], patch_mp->entry[i][j].r);
          mpf_t_to_rat(regen->patchCoeff_rat[i][j][1], patch_mp->entry[i][j].i);

          mpf_set_q(regen->patchCoeff_mp->entry[i][j].r, regen->patchCoeff_rat[i][j][0]);
          mpf_set_q(regen->patchCoeff_mp->entry[i][j].i, regen->patchCoeff_rat[i][j][1]);
          regen->patchCoeff_d->entry[i][j].r = mpq_get_d(regen->patchCoeff_rat[i][j][0]);
          regen->patchCoeff_d->entry[i][j].i = mpq_get_d(regen->patchCoeff_rat[i][j][1]);
        }
    }
    else
    { // copy _rat to _rat
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_set(regen->patchCoeff_rat[i][j][0], patch_rat[i][j][0]);
          mpq_set(regen->patchCoeff_rat[i][j][1], patch_rat[i][j][1]);

          mpf_set_q(regen->patchCoeff_mp->entry[i][j].r, regen->patchCoeff_rat[i][j][0]);
          mpf_set_q(regen->patchCoeff_mp->entry[i][j].i, regen->patchCoeff_rat[i][j][1]);
          regen->patchCoeff_d->entry[i][j].r = mpq_get_d(regen->patchCoeff_rat[i][j][0]);
          regen->patchCoeff_d->entry[i][j].i = mpq_get_d(regen->patchCoeff_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  }

  // read in the main coefficients
  if (inputMPType == 0)
  { // use _d
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_d(patch_d, rows, cols);
    patch_d->rows = rows;
    patch_d->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        fscanf(INPUT, "%lf %lf\n", &patch_d->entry[i][j].r, &patch_d->entry[i][j].i);
      }
    fscanf(INPUT, "\n");
  }
  else if (inputMPType == 1)
  { // use _mp
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_mp(patch_mp, rows, cols);
    patch_mp->rows = rows;
    patch_mp->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpf_inp_str(patch_mp->entry[i][j].r, INPUT, 10);
        mpf_inp_str(patch_mp->entry[i][j].i, INPUT, 10);
        fscanf(INPUT, "\n");
      }
    fscanf(INPUT, "\n");
  }
  else
  { // use _rat
    fscanf(INPUT, "%d %d\n", &rows, &cols);
    change_size_mat_d(patch_d, rows, cols);
    change_size_mat_mp(patch_mp, rows, cols);
    init_mat_rat(patch_rat, rows, cols);
    patch_d->rows = patch_mp->rows = rows;
    patch_d->cols = patch_mp->cols = cols;
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_inp_str(patch_rat[i][j][0], INPUT, 10);
        mpq_canonicalize(patch_rat[i][j][0]);
        mpq_inp_str(patch_rat[i][j][1], INPUT, 10);
        mpq_canonicalize(patch_rat[i][j][1]);
        fscanf(INPUT, "\n");
      }
    fscanf(INPUT, "\n");
  }

  // setup main coefficients
  if (T->MPType == 0)
  { // update mainCoeff_d
    if (inputMPType == 0)
    { // copy _d
      mat_cp_d(regen->mainCoeff_d, patch_d);
    }
    else if (inputMPType == 1)
    { // copy _mp to _d
      mat_mp_to_d(regen->mainCoeff_d, patch_mp);
    }
    else
    { // copy _rat to _d
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          regen->mainCoeff_d->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
          regen->mainCoeff_d->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  } 
  else if (T->MPType == 1)
  { // update mainCoeff_mp
    if (inputMPType == 0)
    { // copy _d to _mp
      mat_d_to_mp(regen->mainCoeff_mp, patch_d);
    }
    else if (inputMPType == 1)
    { // copy _mp
      mat_cp_mp(regen->mainCoeff_mp, patch_mp);
    }
    else
    { // copy _rat to _mp
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_set_q(regen->mainCoeff_mp->entry[i][j].r, patch_rat[i][j][0]);
          mpf_set_q(regen->mainCoeff_mp->entry[i][j].i, patch_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  }
  else
  { // update _d, _mp & _rat
    if (inputMPType == 0)
    { // copy _d to _rat
      comp_mp tempComp;
      init_mp(tempComp);

      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          comp_d_to_mp_rat(tempComp, regen->mainCoeff_rat[i][j], &patch_d->entry[i][j], 16, regen->curr_precision, 0, 0);

          mpf_set_q(regen->mainCoeff_mp->entry[i][j].r, regen->mainCoeff_rat[i][j][0]);
          mpf_set_q(regen->mainCoeff_mp->entry[i][j].i, regen->mainCoeff_rat[i][j][1]);
          regen->mainCoeff_d->entry[i][j].r = mpq_get_d(regen->mainCoeff_rat[i][j][0]);
          regen->mainCoeff_d->entry[i][j].i = mpq_get_d(regen->mainCoeff_rat[i][j][1]);
        }
      clear_mp(tempComp);
    }
    else if (inputMPType == 1)
    { // copy _mp to _rat
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_t_to_rat(regen->mainCoeff_rat[i][j][0], patch_mp->entry[i][j].r);
          mpf_t_to_rat(regen->mainCoeff_rat[i][j][1], patch_mp->entry[i][j].i);

          mpf_set_q(regen->mainCoeff_mp->entry[i][j].r, regen->mainCoeff_rat[i][j][0]);
          mpf_set_q(regen->mainCoeff_mp->entry[i][j].i, regen->mainCoeff_rat[i][j][1]);
          regen->mainCoeff_d->entry[i][j].r = mpq_get_d(regen->mainCoeff_rat[i][j][0]);
          regen->mainCoeff_d->entry[i][j].i = mpq_get_d(regen->mainCoeff_rat[i][j][1]);
        }
    }
    else
    { // copy _rat to _rat
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_set(regen->mainCoeff_rat[i][j][0], patch_rat[i][j][0]);
          mpq_set(regen->mainCoeff_rat[i][j][1], patch_rat[i][j][1]);

          mpf_set_q(regen->mainCoeff_mp->entry[i][j].r, regen->mainCoeff_rat[i][j][0]);
          mpf_set_q(regen->mainCoeff_mp->entry[i][j].i, regen->mainCoeff_rat[i][j][1]);
          regen->mainCoeff_d->entry[i][j].r = mpq_get_d(regen->mainCoeff_rat[i][j][0]);
          regen->mainCoeff_d->entry[i][j].i = mpq_get_d(regen->mainCoeff_rat[i][j][1]);
        }
      // clear patch_rat
      clear_mat_rat(patch_rat, rows, cols);
    }
  }

  // read in the coeff
  if (inputMPType == 0)
  { // use _d
    coeff_d = (comp_d ****)bmalloc(regen->num_funcs * sizeof(comp_d ***));
    for (i = 0; i < regen->num_funcs; i++)
    {
      coeff_d[i] = (comp_d ***)bmalloc(regen->num_var_gps * sizeof(comp_d **));
      for (j = 0; j < regen->num_var_gps; j++)
      {
        coeff_d[i][j] = (comp_d **)bmalloc(regen->degrees[i][j] * sizeof(comp_d *));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          coeff_d[i][j][k] = (comp_d *)bmalloc(regen->num_variables * sizeof(comp_d));
          for (l = 0; l < regen->num_variables; l++)
          { 
            fscanf(INPUT, "%lf %lf\n", &coeff_d[i][j][k][l]->r, &coeff_d[i][j][k][l]->i);
          }
        }
      }
    }
  }
  else if (inputMPType == 1)
  { // use _mp
    coeff_mp = (comp_mp ****)bmalloc(regen->num_funcs * sizeof(comp_mp ***));
    for (i = 0; i < regen->num_funcs; i++)
    {
      coeff_mp[i] = (comp_mp ***)bmalloc(regen->num_var_gps * sizeof(comp_mp **));
      for (j = 0; j < regen->num_var_gps; j++)
      {
        coeff_mp[i][j] = (comp_mp **)bmalloc(regen->degrees[i][j] * sizeof(comp_mp *));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          coeff_mp[i][j][k] = (comp_mp *)bmalloc(regen->num_variables * sizeof(comp_mp));
          for (l = 0; l < regen->num_variables; l++)
          {
            init_mp(coeff_mp[i][j][k][l]);
            mpf_inp_str(coeff_mp[i][j][k][l]->r, INPUT, 10);
            mpf_inp_str(coeff_mp[i][j][k][l]->i, INPUT, 10);
            fscanf(INPUT, "\n");
          }
        }
      }
    }
  }
  else
  { // use _rat
    coeff_rat = (mpq_t *****)bmalloc(regen->num_funcs * sizeof(mpq_t ****));
    for (i = 0; i < regen->num_funcs; i++)
    {
      coeff_rat[i] = (mpq_t ****)bmalloc(regen->num_var_gps * sizeof(mpq_t ***));
      for (j = 0; j < regen->num_var_gps; j++)
      {
        coeff_rat[i][j] = (mpq_t ***)bmalloc(regen->degrees[i][j] * sizeof(mpq_t **));
        for (k = 0; k < regen->degrees[i][j]; k++)
        {
          coeff_rat[i][j][k] = (mpq_t **)bmalloc(regen->num_variables * sizeof(mpq_t *));
          for (l = 0; l < regen->num_variables; l++)
          {
            coeff_rat[i][j][k][l] = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
            init_rat(coeff_rat[i][j][k][l]);

            mpq_inp_str(coeff_rat[i][j][k][l][0], INPUT, 10);
            mpq_canonicalize(coeff_rat[i][j][k][l][0]);

            mpq_inp_str(coeff_rat[i][j][k][l][1], INPUT, 10);
            mpq_canonicalize(coeff_rat[i][j][k][l][1]);

            fscanf(INPUT, "\n");
          }
        }
      }
    }
  }

  // setup coeff
  if (T->MPType == 0)
  { // update coeff_d
    if (inputMPType == 0)
    { // copy _d
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              set_d(regen->coeff_d[i][j][k][l], coeff_d[i][j][k][l]);
            }
            free(coeff_d[i][j][k]);
          }
          free(coeff_d[i][j]);
        }
        free(coeff_d[i]);
      }
      free(coeff_d);
    }
    else if (inputMPType == 1)
    { // copy _mp to _d
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              mp_to_d(regen->coeff_d[i][j][k][l], coeff_mp[i][j][k][l]);
              clear_mp(coeff_mp[i][j][k][l]);
            }
            free(coeff_mp[i][j][k]);
          }
          free(coeff_mp[i][j]);
        }
        free(coeff_mp[i]);
      }
      free(coeff_mp);
    }
    else
    { // copy _rat to _d
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              rat_to_d(regen->coeff_d[i][j][k][l], coeff_rat[i][j][k][l]);
              clear_rat(coeff_rat[i][j][k][l]); 
              free(coeff_rat[i][j][k][l]);
            }
            free(coeff_rat[i][j][k]);
          }
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
    }
  }
  else if (T->MPType == 1)
  { // update coeff_mp
    if (inputMPType == 0)
    { // copy _d to _mp
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              d_to_mp(regen->coeff_mp[i][j][k][l], coeff_d[i][j][k][l]);
            }
            free(coeff_d[i][j][k]);
          }
          free(coeff_d[i][j]);
        }
        free(coeff_d[i]);
      }
      free(coeff_d);
    }
    else if (inputMPType == 1)
    { // copy _mp
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              set_mp(regen->coeff_mp[i][j][k][l], coeff_mp[i][j][k][l]);
              clear_mp(coeff_mp[i][j][k][l]);
            }
            free(coeff_mp[i][j][k]);
          }
          free(coeff_mp[i][j]);
        }
        free(coeff_mp[i]);
      }
      free(coeff_mp);
    }
    else
    { // copy _rat to _mp
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              rat_to_mp(regen->coeff_mp[i][j][k][l], coeff_rat[i][j][k][l]);
              clear_rat(coeff_rat[i][j][k][l]);
              free(coeff_rat[i][j][k][l]);
            }
            free(coeff_rat[i][j][k]);
          }
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
    }
  }
  else
  { // update _d, _mp & _rat
    if (inputMPType == 0)
    { // copy _d to _rat
      comp_mp tempComp;
      init_mp(tempComp);

      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              comp_d_to_mp_rat(tempComp, regen->coeff_rat[i][j][k][l], coeff_d[i][j][k][l], 16, regen->curr_precision, 0, 0);
              rat_to_mp(regen->coeff_mp[i][j][k][l], regen->coeff_rat[i][j][k][l]);
              rat_to_d(regen->coeff_d[i][j][k][l], regen->coeff_rat[i][j][k][l]);
            }
            free(coeff_d[i][j][k]);
          }
          free(coeff_d[i][j]);
        }
        free(coeff_d[i]);
      }
      free(coeff_d);

      clear_mp(tempComp);
    }
    else if (inputMPType == 1)
    { // copy _mp to _rat
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              mpf_t_to_rat(regen->coeff_rat[i][j][k][l][0], coeff_mp[i][j][k][l]->r);
              mpf_t_to_rat(regen->coeff_rat[i][j][k][l][1], coeff_mp[i][j][k][l]->i);

              rat_to_mp(regen->coeff_mp[i][j][k][l], regen->coeff_rat[i][j][k][l]);
              rat_to_d(regen->coeff_d[i][j][k][l], regen->coeff_rat[i][j][k][l]);

              clear_mp(coeff_mp[i][j][k][l]);
            }
            free(coeff_mp[i][j][k]);
          }
          free(coeff_mp[i][j]);
        }
        free(coeff_mp[i]);
      }
      free(coeff_mp);
    }
    else
    { // copy _rat to _rat
      for (i = regen->num_funcs - 1; i >= 0; i--)
      {
        for (j = regen->num_var_gps - 1; j >= 0; j--)
        {
          for (k = regen->degrees[i][j] - 1; k >= 0; k--)
          {
            for (l = regen->num_variables - 1; l >= 0; l--)
            {
              set_rat(regen->coeff_rat[i][j][k][l], coeff_rat[i][j][k][l]);
              rat_to_mp(regen->coeff_mp[i][j][k][l], regen->coeff_rat[i][j][k][l]);
              rat_to_d(regen->coeff_d[i][j][k][l], regen->coeff_rat[i][j][k][l]);
              clear_rat(coeff_rat[i][j][k][l]);
              free(coeff_rat[i][j][k][l]);
            }
            free(coeff_rat[i][j][k]);
          }
          free(coeff_rat[i][j]);
        }
        free(coeff_rat[i]);
      }
      free(coeff_rat);
    }
  }

  // read in p
  if (inputMPType == 0)
  { // use _d
    fscanf(INPUT, "%d\n", &rows);
    change_size_vec_d(p_d, rows);
    p_d->size = rows;
    for (i = 0; i < rows; i++)
    {
      fscanf(INPUT, "%lf %lf\n", &p_d->coord[i].r, &p_d->coord[i].i);
    }
    fscanf(INPUT, "\n");
  }
  else if (inputMPType == 1)
  { // use _mp
    fscanf(INPUT, "%d\n", &rows);
    change_size_vec_mp(p_mp, rows);
    p_mp->size = rows;
    for (i = 0; i < rows; i++)
    {
      mpf_inp_str(p_mp->coord[i].r, INPUT, 10);
      mpf_inp_str(p_mp->coord[i].i, INPUT, 10);
      fscanf(INPUT, "\n");
    }
    fscanf(INPUT, "\n");
  }
  else
  { // use _rat
    fscanf(INPUT, "%d\n", &rows);
    change_size_vec_d(p_d, rows);
    change_size_vec_mp(p_mp, rows);
    init_vec_rat(p_rat, rows);
    p_d->size = p_mp->size = rows;
    for (i = 0; i < rows; i++)
    {
      mpq_inp_str(p_rat[i][0], INPUT, 10);
      mpq_canonicalize(p_rat[i][0]);
      mpq_inp_str(p_rat[i][1], INPUT, 10);
      mpq_canonicalize(p_rat[i][1]);
      fscanf(INPUT, "\n");
    }
    fscanf(INPUT, "\n");
  }

  // setup p
  if (T->MPType == 0)
  { // setup p_d
    init_vec_d(regen->main_homVar_d, rows);
    regen->main_homVar_d->size = rows;
    if (inputMPType == 0)
    { // copy _d
      vec_cp_d(regen->main_homVar_d, p_d);
    }
    else if (inputMPType == 1)
    { // copy _mp to _d
      vec_mp_to_d(regen->main_homVar_d, p_mp);
    }
    else
    { // copy _rat to _d
      for (i = 0; i < rows; i++)
      {
        regen->main_homVar_d->coord[i].r = mpq_get_d(p_rat[i][0]);
        regen->main_homVar_d->coord[i].i = mpq_get_d(p_rat[i][1]);
      }
      // clear p_rat
      clear_vec_rat(p_rat, rows);
    }
  }
  else if (T->MPType == 1)
  { // setup p_mp
    init_vec_mp(regen->main_homVar_mp, rows);
    regen->main_homVar_mp->size = rows;
    if (inputMPType == 0)
    { // copy _d to _mp
      vec_d_to_mp(regen->main_homVar_mp, p_d);
    }
    else if (inputMPType == 1)
    { // copy _mp
      vec_cp_mp(regen->main_homVar_mp, p_mp);
    }
    else
    { // copy _rat to _mp
      for (i = 0; i < rows; i++)
      {
        mpf_set_q(regen->main_homVar_mp->coord[i].r, p_rat[i][0]);
        mpf_set_q(regen->main_homVar_mp->coord[i].i, p_rat[i][1]);
      }
      // clear p_rat
      clear_vec_rat(p_rat, rows);
    }
  }
  else
  { // setup p_d, p_mp & p_rat
    init_vec_d(regen->main_homVar_d, rows);
    init_vec_mp(regen->main_homVar_mp, rows);
    init_vec_rat(regen->main_homVar_rat, rows);
    regen->main_homVar_d->size = regen->main_homVar_mp->size = rows;
    if (inputMPType == 0)
    { // copy _d to _rat
      comp_mp tempComp;
      init_mp(tempComp);

      for (i = 0; i < rows; i++)
      {
        comp_d_to_mp_rat(tempComp, regen->main_homVar_rat[i], &p_d->coord[i], 16, regen->curr_precision, 0, 0);
        rat_to_d(&regen->main_homVar_d->coord[i], regen->main_homVar_rat[i]);
        rat_to_mp(&regen->main_homVar_mp->coord[i], regen->main_homVar_rat[i]);
      }
      clear_mp(tempComp);
    }
    else if (inputMPType == 1)
    { // copy _mp to _rat
      for (i = 0; i < rows; i++)
      {
        mp_to_rat(regen->main_homVar_rat[i], &p_mp->coord[i]);
        rat_to_d(&regen->main_homVar_d->coord[i], regen->main_homVar_rat[i]);
        rat_to_mp(&regen->main_homVar_mp->coord[i], regen->main_homVar_rat[i]);
      }
    }
    else
    { // copy _rat to _rat
      for (i = 0; i < rows; i++)
      {
        set_rat(regen->main_homVar_rat[i], p_rat[i]);
        rat_to_d(&regen->main_homVar_d->coord[i], regen->main_homVar_rat[i]);
        rat_to_mp(&regen->main_homVar_mp->coord[i], regen->main_homVar_rat[i]);
      }
      // clear p_rat
      clear_vec_rat(p_rat, rows);
    }
  }

  // setup the intrinsic slice, if needed
  if (regen->level[curr_level].useIntrinsicSlice)
  { // read in B
    if (inputMPType == 0)
    { // use _d
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_d(patch_d, rows, cols);
      patch_d->rows = rows;
      patch_d->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
          fscanf(INPUT, "%lf %lf\n", &patch_d->entry[i][j].r, &patch_d->entry[i][j].i);
      fscanf(INPUT, "\n");
    } 
    else if (inputMPType == 1)
    { // use _mp
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_mp(patch_mp, rows, cols);
      patch_mp->rows = rows;
      patch_mp->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_inp_str(patch_mp->entry[i][j].r, INPUT, 10);
          mpf_inp_str(patch_mp->entry[i][j].i, INPUT, 10);
          fscanf(INPUT, "\n");
        }
      fscanf(INPUT, "\n");
    } 
    else
    { // use _rat
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_d(patch_d, rows, cols);
      change_size_mat_mp(patch_mp, rows, cols);
      init_mat_rat(patch_rat, rows, cols);
      patch_d->rows = patch_mp->rows = rows;
      patch_d->cols = patch_mp->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_inp_str(patch_rat[i][j][0], INPUT, 10);
          mpq_canonicalize(patch_rat[i][j][0]);
          mpq_inp_str(patch_rat[i][j][1], INPUT, 10);
          mpq_canonicalize(patch_rat[i][j][1]);
          fscanf(INPUT, "\n");
        }
      fscanf(INPUT, "\n");
    } 

    // setup B
    if (T->MPType == 0)
    { // setup B_d
      init_mat_d(regen->level[curr_level].B_d, rows, cols);
      regen->level[curr_level].B_d->rows = rows;
      regen->level[curr_level].B_d->cols = cols;
      if (inputMPType == 0)
      { // copy _d
        mat_cp_d(regen->level[curr_level].B_d, patch_d);
      }
      else if (inputMPType == 1)
      { // copy _mp to _d
        mat_mp_to_d(regen->level[curr_level].B_d, patch_mp);
      }
      else
      { // copy _rat to _d
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            regen->level[curr_level].B_d->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
            regen->level[curr_level].B_d->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }
    else if (T->MPType == 1)
    { // setup B_mp
      init_mat_mp(regen->level[curr_level].B_mp, rows, cols);
      regen->level[curr_level].B_mp->rows = rows;
      regen->level[curr_level].B_mp->cols = cols;
      if (inputMPType == 0)
      { // copy _d to _mp
        mat_d_to_mp(regen->level[curr_level].B_mp, patch_d);
      }
      else if (inputMPType == 1)
      { // copy _mp
        mat_cp_mp(regen->level[curr_level].B_mp, patch_mp);
      }
      else
      { // copy _rat to _mp
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].r, patch_rat[i][j][0]);
            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].i, patch_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }
    else
    { // setup B_d, B_mp & B_rat
      init_mat_d(regen->level[curr_level].B_d, rows, cols);
      init_mat_mp(regen->level[curr_level].B_mp, rows, cols);
      init_mat_rat(regen->level[curr_level].B_rat, rows, cols);
      regen->level[curr_level].B_d->rows = regen->level[curr_level].B_mp->rows = rows;
      regen->level[curr_level].B_d->cols = regen->level[curr_level].B_mp->cols = cols;
      if (inputMPType == 0)
      { // copy _d to _rat
        comp_mp tempComp;
        init_mp(tempComp);

        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            comp_d_to_mp_rat(tempComp, regen->level[curr_level].B_rat[i][j], &patch_d->entry[i][j], 16, regen->curr_precision, 0, 0);

            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].r, regen->level[curr_level].B_rat[i][j][0]);
            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].i, regen->level[curr_level].B_rat[i][j][1]);
            regen->level[curr_level].B_d->entry[i][j].r = mpq_get_d(regen->level[curr_level].B_rat[i][j][0]);
            regen->level[curr_level].B_d->entry[i][j].i = mpq_get_d(regen->level[curr_level].B_rat[i][j][1]);
          }
        clear_mp(tempComp);
      }
      else if (inputMPType == 1)
      { // copy _mp to _rat
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpf_t_to_rat(regen->level[curr_level].B_rat[i][j][0], patch_mp->entry[i][j].r);
            mpf_t_to_rat(regen->level[curr_level].B_rat[i][j][1], patch_mp->entry[i][j].i);

            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].r, regen->level[curr_level].B_rat[i][j][0]);
            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].i, regen->level[curr_level].B_rat[i][j][1]);
            regen->level[curr_level].B_d->entry[i][j].r = mpq_get_d(regen->level[curr_level].B_rat[i][j][0]);
            regen->level[curr_level].B_d->entry[i][j].i = mpq_get_d(regen->level[curr_level].B_rat[i][j][1]);
          }
      }
      else
      { // copy _rat to _rat
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpq_set(regen->level[curr_level].B_rat[i][j][0], patch_rat[i][j][0]);
            mpq_set(regen->level[curr_level].B_rat[i][j][1], patch_rat[i][j][1]);

            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].r, regen->level[curr_level].B_rat[i][j][0]);
            mpf_set_q(regen->level[curr_level].B_mp->entry[i][j].i, regen->level[curr_level].B_rat[i][j][1]);
            regen->level[curr_level].B_d->entry[i][j].r = mpq_get_d(regen->level[curr_level].B_rat[i][j][0]);
            regen->level[curr_level].B_d->entry[i][j].i = mpq_get_d(regen->level[curr_level].B_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }

    // read in p
    if (inputMPType == 0)
    { // use _d
      fscanf(INPUT, "%d\n", &rows);
      change_size_vec_d(p_d, rows);
      p_d->size = rows;
      for (i = 0; i < rows; i++)
      {
        fscanf(INPUT, "%lf %lf\n", &p_d->coord[i].r, &p_d->coord[i].i);
      }
      fscanf(INPUT, "\n");
    }
    else if (inputMPType == 1)
    { // use _mp
      fscanf(INPUT, "%d\n", &rows);
      change_size_vec_mp(p_mp, rows);
      p_mp->size = rows;
      for (i = 0; i < rows; i++)
      {
        mpf_inp_str(p_mp->coord[i].r, INPUT, 10);
        mpf_inp_str(p_mp->coord[i].i, INPUT, 10);
        fscanf(INPUT, "\n");
      }
      fscanf(INPUT, "\n");
    }
    else
    { // use _rat
      fscanf(INPUT, "%d\n", &rows);
      change_size_vec_d(p_d, rows);
      change_size_vec_mp(p_mp, rows);
      init_vec_rat(p_rat, rows);
      p_d->size = p_mp->size = rows;
      for (i = 0; i < rows; i++)
      {
        mpq_inp_str(p_rat[i][0], INPUT, 10);
        mpq_canonicalize(p_rat[i][0]);
        mpq_inp_str(p_rat[i][1], INPUT, 10);
        mpq_canonicalize(p_rat[i][1]);
        fscanf(INPUT, "\n");
      }
      fscanf(INPUT, "\n");
    }

    // setup p
    if (T->MPType == 0)
    { // setup p_d
      init_vec_d(regen->level[curr_level].p_d, rows);
      regen->level[curr_level].p_d->size = rows;
      if (inputMPType == 0)
      { // copy _d
        vec_cp_d(regen->level[curr_level].p_d, p_d);
      }
      else if (inputMPType == 1)
      { // copy _mp to _d
        vec_mp_to_d(regen->level[curr_level].p_d, p_mp);
      }
      else
      { // copy _rat to _d
        for (i = 0; i < rows; i++)
        {
          regen->level[curr_level].p_d->coord[i].r = mpq_get_d(p_rat[i][0]);
          regen->level[curr_level].p_d->coord[i].i = mpq_get_d(p_rat[i][1]);
        }
        // clear p_rat
        clear_vec_rat(p_rat, rows);
      }
    }
    else if (T->MPType == 1)
    { // setup p_mp
      init_vec_mp(regen->level[curr_level].p_mp, rows);
      regen->level[curr_level].p_mp->size = rows;
      if (inputMPType == 0)
      { // copy _d to _mp
        vec_d_to_mp(regen->level[curr_level].p_mp, p_d);
      }
      else if (inputMPType == 1)
      { // copy _mp
        vec_cp_mp(regen->level[curr_level].p_mp, p_mp);
      }
      else
      { // copy _rat to _mp
        for (i = 0; i < rows; i++)
        {
          mpf_set_q(regen->level[curr_level].p_mp->coord[i].r, p_rat[i][0]);
          mpf_set_q(regen->level[curr_level].p_mp->coord[i].i, p_rat[i][1]);
        }
        // clear p_rat
        clear_vec_rat(p_rat, rows);
      }
    }
    else
    { // setup p_d, p_mp & p_rat
      init_vec_d(regen->level[curr_level].p_d, rows);
      init_vec_mp(regen->level[curr_level].p_mp, rows);
      init_vec_rat(regen->level[curr_level].p_rat, rows);
      regen->level[curr_level].p_d->size = regen->level[curr_level].p_mp->size = rows;
      if (inputMPType == 0)
      { // copy _d to _rat
        comp_mp tempComp;
        init_mp(tempComp);

        for (i = 0; i < rows; i++)
        {
          comp_d_to_mp_rat(tempComp, regen->level[curr_level].p_rat[i], &p_d->coord[i], 16, regen->curr_precision, 0, 0);

          mpf_set_q(regen->level[curr_level].p_mp->coord[i].r, regen->level[curr_level].p_rat[i][0]);
          mpf_set_q(regen->level[curr_level].p_mp->coord[i].i, regen->level[curr_level].p_rat[i][1]);
          regen->level[curr_level].p_d->coord[i].r = mpq_get_d(regen->level[curr_level].p_rat[i][0]);
          regen->level[curr_level].p_d->coord[i].i = mpq_get_d(regen->level[curr_level].p_rat[i][1]);
        }
        clear_mp(tempComp);
      }
      else if (inputMPType == 1)
      { // copy _mp to _rat
        for (i = 0; i < rows; i++)
        {
          mpf_t_to_rat(regen->level[curr_level].p_rat[i][0], p_mp->coord[i].r);
          mpf_t_to_rat(regen->level[curr_level].p_rat[i][1], p_mp->coord[i].i);

          mpf_set_q(regen->level[curr_level].p_mp->coord[i].r, regen->level[curr_level].p_rat[i][0]);
          mpf_set_q(regen->level[curr_level].p_mp->coord[i].i, regen->level[curr_level].p_rat[i][1]);
          regen->level[curr_level].p_d->coord[i].r = mpq_get_d(regen->level[curr_level].p_rat[i][0]);
          regen->level[curr_level].p_d->coord[i].i = mpq_get_d(regen->level[curr_level].p_rat[i][1]);
        }
      }
      else
      { // copy _rat to _rat
        for (i = 0; i < rows; i++)
        {
          mpq_set(regen->level[curr_level].p_rat[i][0], p_rat[i][0]);
          mpq_set(regen->level[curr_level].p_rat[i][1], p_rat[i][1]);

          mpf_set_q(regen->level[curr_level].p_mp->coord[i].r, regen->level[curr_level].p_rat[i][0]);
          mpf_set_q(regen->level[curr_level].p_mp->coord[i].i, regen->level[curr_level].p_rat[i][1]);
          regen->level[curr_level].p_d->coord[i].r = mpq_get_d(regen->level[curr_level].p_rat[i][0]);
          regen->level[curr_level].p_d->coord[i].i = mpq_get_d(regen->level[curr_level].p_rat[i][1]);
        }
        // clear p_rat
        clear_vec_rat(p_rat, rows);
      }
    }
  }
  else
  { // do not need to setup B & p
    regen->level[curr_level].B_d->rows = regen->level[curr_level].B_d->cols = regen->level[curr_level].B_mp->rows = regen->level[curr_level].B_mp->cols = 0;
    regen->level[curr_level].B_rat = NULL;
    regen->level[curr_level].p_d->size = regen->level[curr_level].p_mp->size = 0;
    regen->level[curr_level].p_rat = NULL;
  }

  if (!regen->noChanges)
  { // setup the squaring matrix A
    if (inputMPType == 0)
    { // use _d
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_d(patch_d, rows, cols);
      patch_d->rows = rows;
      patch_d->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          fscanf(INPUT, "%lf %lf\n", &patch_d->entry[i][j].r, &patch_d->entry[i][j].i);
        }
      fscanf(INPUT, "\n");
    }
    else if (inputMPType == 1)
    { // use _mp
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_mp(patch_mp, rows, cols);
      patch_mp->rows = rows;
      patch_mp->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_inp_str(patch_mp->entry[i][j].r, INPUT, 10);
          mpf_inp_str(patch_mp->entry[i][j].i, INPUT, 10);
          fscanf(INPUT, "\n");
        }
      fscanf(INPUT, "\n");
    }
    else
    { // use _rat
      fscanf(INPUT, "%d %d\n", &rows, &cols);
      change_size_mat_d(patch_d, rows, cols);
      change_size_mat_mp(patch_mp, rows, cols);
      init_mat_rat(patch_rat, rows, cols);
      patch_d->rows = patch_mp->rows = rows;
      patch_d->cols = patch_mp->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpq_inp_str(patch_rat[i][j][0], INPUT, 10);
          mpq_canonicalize(patch_rat[i][j][0]);
          mpq_inp_str(patch_rat[i][j][1], INPUT, 10);
          mpq_canonicalize(patch_rat[i][j][1]);
          fscanf(INPUT, "\n");
        }
      fscanf(INPUT, "\n");
    }

    // setup A
    if (T->MPType == 0)
    { // update square_d
      square_system_eval_data_d *SSED = (square_system_eval_data_d *)regen->square_d;
      if (inputMPType == 0)
      { // copy _d
        mat_cp_d(SSED->A, patch_d);
      }
      else if (inputMPType == 1)
      { // copy _mp to _d
        mat_mp_to_d(SSED->A, patch_mp);
      }
      else
      { // copy _rat to _d
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            SSED->A->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
            SSED->A->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }
    else if (T->MPType == 1)
    { // update square_mp
      square_system_eval_data_mp *SSED = (square_system_eval_data_mp *)regen->square_mp;
      if (inputMPType == 0)
      { // copy _d to _mp
        mat_d_to_mp(SSED->A, patch_d);
      }
      else if (inputMPType == 1)
      { // copy _mp
        mat_cp_mp(SSED->A, patch_mp);
      }
      else
      { // copy _rat to _mp
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpf_set_q(SSED->A->entry[i][j].r, patch_rat[i][j][0]);
            mpf_set_q(SSED->A->entry[i][j].i, patch_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }
    else
    { // update _d, _mp & _rat
      square_system_eval_data_d *SSED_d = (square_system_eval_data_d *)regen->square_d;
      square_system_eval_data_mp *SSED_mp = (square_system_eval_data_mp *)regen->square_mp;
      if (inputMPType == 0)
      { // copy _d to _rat
        comp_mp tempComp;
        init_mp(tempComp);

        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            comp_d_to_mp_rat(tempComp, SSED_mp->A_rat[i][j], &patch_d->entry[i][j], 16, regen->curr_precision, 0, 0);

            mpf_set_q(SSED_mp->A->entry[i][j].r, SSED_mp->A_rat[i][j][0]);
            mpf_set_q(SSED_mp->A->entry[i][j].i, SSED_mp->A_rat[i][j][1]);
            SSED_d->A->entry[i][j].r = mpq_get_d(SSED_mp->A_rat[i][j][0]);
            SSED_d->A->entry[i][j].i = mpq_get_d(SSED_mp->A_rat[i][j][1]);
          }
        clear_mp(tempComp);
      }
      else if (inputMPType == 1)
      { // copy _mp to _rat
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpf_t_to_rat(SSED_mp->A_rat[i][j][0], patch_mp->entry[i][j].r);
            mpf_t_to_rat(SSED_mp->A_rat[i][j][1], patch_mp->entry[i][j].i);

            mpf_set_q(SSED_mp->A->entry[i][j].r, SSED_mp->A_rat[i][j][0]);
            mpf_set_q(SSED_mp->A->entry[i][j].i, SSED_mp->A_rat[i][j][1]);
            SSED_d->A->entry[i][j].r = mpq_get_d(SSED_mp->A_rat[i][j][0]);
            SSED_d->A->entry[i][j].i = mpq_get_d(SSED_mp->A_rat[i][j][1]);
          }
      }
      else
      { // copy _rat to _rat
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpq_set(SSED_mp->A_rat[i][j][0], patch_rat[i][j][0]);
            mpq_set(SSED_mp->A_rat[i][j][1], patch_rat[i][j][1]);

            mpf_set_q(SSED_mp->A->entry[i][j].r, SSED_mp->A_rat[i][j][0]);
            mpf_set_q(SSED_mp->A->entry[i][j].i, SSED_mp->A_rat[i][j][1]);
            SSED_d->A->entry[i][j].r = mpq_get_d(SSED_mp->A_rat[i][j][0]);
            SSED_d->A->entry[i][j].i = mpq_get_d(SSED_mp->A_rat[i][j][1]);
          }
        // clear patch_rat
        clear_mat_rat(patch_rat, rows, cols);
      }
    }
  }

  // initialize other values
  regen->level[curr_level].num_sing = regen->level[curr_level].num_nonsing = regen->level[curr_level].num_inf
    = regen->level[curr_level].num_higher_dim = regen->level[curr_level].num_bad = 0;

  // close INPUT
  fclose(INPUT);

  // clear patch & p
  clear_mat_d(patch_d);
  clear_mat_mp(patch_mp);
  clear_vec_d(p_d);
  clear_vec_mp(p_mp);

  return;
}

