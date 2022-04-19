// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include <stdio.h>
#include <math.h>
#include "bertini.h"
#include "cascade.h"

int matrixIsZero_d(mat_d A);
int matrixIsZero_mp(mat_mp A);

int rank_finder_d(preproc_data *PPD, prog_t *orig_function, tracker_config_t *T, int num_variables)  
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: rank of the system                             *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, numGps = PPD->num_var_gp + PPD->num_hom_var_gp, currVar = PPD->num_var_gp, currVarGp = 0, final_rank = -1;
  eval_struct_d e;
  mat_d newMat;
  point_d test_point;
  comp_d test_time;

  init_eval_struct_d(e, orig_function->numFuncs, orig_function->numVars, orig_function->numPars);
  init_mat_d(newMat, 0, 0);

  // random point
  init_point_d(test_point, num_variables);
  make_vec_random_d(test_point, num_variables);
  set_one_d(test_time);

  // Evaluate the original function to get the Jacobian at the point  
  evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, test_point, test_time, orig_function);

  if (infNormVec_d(e.funcVals) < T->final_tolerance)
  {
    printf("ERROR: The system is numerically 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  //Form a new matrix [[Jv];[patches]]
  change_size_mat_d(newMat, e.Jv->rows + numGps, e.Jv->cols);
  newMat->rows = e.Jv->rows + numGps;
  newMat->cols = e.Jv->cols;
  // top is Jv
  for (i = 0; i < e.Jv->rows; i++)
    for (j = 0; j < e.Jv->cols; j++)
      set_d(&newMat->entry[i][j], &e.Jv->entry[i][j]);
  // zero out bottom of Jv
  for (i = e.Jv->rows; i < newMat->rows; i++)
    for (j = 0; j < e.Jv->cols; j++)
      set_zero_d(&newMat->entry[i][j]);
  // setup patches
  for (i = 0; i < numGps; i++)
  { // put random in correct location
    if (PPD->type[i])
    { // extra homogenizing variable
      get_comp_rand_d(&newMat->entry[e.Jv->rows + i][currVarGp]);
      currVarGp++;
    }
    for (j = 0; j < PPD->size[i]; j++)
    { // variables in group
      get_comp_rand_d(&newMat->entry[e.Jv->rows + i][currVar]);
      currVar++;
    }
  }

  //Now we get the rank of the new matrix - do a simple error check to avoid stopping svd!
  if (matrixIsZero_d(newMat))
  {
    printf("ERROR: The Jacobian with respect to the variables is 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  else
  { // perform the svd on the new matrix to find the corank
    j = svd_corank_d(newMat, T->sing_val_zero_tol, 1e12);

    // Now we just determine the rank
    final_rank = MIN(newMat->rows, newMat->cols) - j - numGps;
  }

  if (final_rank == 0)
  {
    printf("ERROR: The rank of the system is 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);  
  }

  clear_eval_struct_d(e);
  clear_mat_d(newMat);
  clear_point_d(test_point);

  return final_rank;
}

int rank_finder_mp(preproc_data *PPD, prog_t *orig_function, tracker_config_t *T, int num_variables)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: rank of the system                             *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, numGps = PPD->num_var_gp + PPD->num_hom_var_gp, currVar = PPD->num_var_gp, currVarGp = 0, final_rank = -1;
  mpf_t tol_sing_val, largeChange;
  mat_mp newMat;
  point_mp test_point;
  eval_struct_mp e;
  comp_mp test_time;

  // initialize
  mpf_init(tol_sing_val); mpf_init(largeChange);
  init_eval_struct_mp(e, orig_function->numFuncs, orig_function->numVars, orig_function->numPars);
  init_mp(test_time);
  init_mat_mp(newMat, 0, 0);

  // random point
  init_point_mp(test_point, num_variables);
  make_vec_random_mp(test_point, num_variables);
  set_one_mp(test_time);

  // Evaluate the original function to get the Jacobian at the point
  evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, test_point, test_time, orig_function);

  if (infNormVec_mp(e.funcVals) < T->final_tolerance)
  {
    printf("ERROR: The system is numerically 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  //Form a new matrix [[Jv];[patches]]
  change_size_mat_mp(newMat, e.Jv->rows + numGps, e.Jv->cols);
  newMat->rows = e.Jv->rows + numGps;
  newMat->cols = e.Jv->cols;
  // top is Jv
  for (i = 0; i < e.Jv->rows; i++)
    for (j = 0; j < e.Jv->cols; j++)
      set_mp(&newMat->entry[i][j], &e.Jv->entry[i][j]);
  // zero out bottom of Jv
  for (i = e.Jv->rows; i < newMat->rows; i++)
    for (j = 0; j < e.Jv->cols; j++)
      set_zero_mp(&newMat->entry[i][j]);
  // setup patches
  for (i = 0; i < numGps; i++)
  { // put random in correct location
    if (PPD->type[i])
    { // extra homogenizing variable
      get_comp_rand_mp(&newMat->entry[e.Jv->rows + i][currVarGp]);
      currVarGp++;
    }
    for (j = 0; j < PPD->size[i]; j++)
    { // variables in group
      get_comp_rand_mp(&newMat->entry[e.Jv->rows + i][currVar]);
      currVar++;
    }
  }

  //Now we get the rank of the new matrix:
  if (matrixIsZero_mp(newMat))
  {
    printf("ERROR: The Jacobian with respect to the variables is 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  else
  { // perform the svd on the new matrix

    // setup tol_sing_val & largeChange
    mpf_set_d(tol_sing_val, T->sing_val_zero_tol);
    mpf_set_d(largeChange, 1e17);

    j = svd_corank_mp(newMat, tol_sing_val, largeChange);

    // Now we just determine the rank
    final_rank = MIN(newMat->rows, newMat->cols) - j - numGps;
  }

  if (final_rank == 0)
  {
    printf("ERROR: The rank of the system is 0!  Please input a non-degenerate system.\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // clear
  mpf_clear(tol_sing_val); mpf_clear(largeChange);
  clear_point_mp(test_point);
  clear_eval_struct_mp(e);
  clear_mp(test_time);
  clear_mat_mp(newMat);

  return final_rank;
}

int matrixIsZero_d(mat_d A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: checks to see if A is all zeros                        *
\***************************************************************/
{
  int i, j;

  for (i = 0; i < A->rows; i++)
    for (j = 0; j < A->cols; j++)
      if (d_abs_d(&A->entry[i][j]) > 0)
        return 0;

  return 1;
}

int matrixIsZero_mp(mat_mp A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: checks to see if A is all zeros                        *
\***************************************************************/
{ 
  int i, j;
 
  for (i = 0; i < A->rows; i++)
    for (j = 0; j < A->cols; j++)
      if (d_abs_mp(&A->entry[i][j]) > 0)
        return 0;
 
  return 1;
}

