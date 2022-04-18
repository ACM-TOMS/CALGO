// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

int rankDef_newton_mp(mpf_t newton_resid, int *newton_digits, mpf_t func_resid, int *func_digits, point_data_mp *in_mp, int curr_prec, tracker_config_t *T, eval_struct_mp *e, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// does an svd analysis of the Jacobian to determine its corank

int jacobian_analyze_d(double current_tol, point_data_d *in, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int relative_svd_change_d(double percent_change, mat_d E1, mat_d E2);

int jacobian_analyze_mp(mpf_t current_tol, point_data_mp *in, int prec_in, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int relative_svd_change_mp(mpf_t percent_change, mat_mp E1, mat_mp E2);

int jacobian_analyze_amp(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

// does newton iterations to determine if rank deficient

int determineRankDef_d(double *CN_ret, double current_tol, point_data_d *in_d, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int determineRankDef_mp(double *CN_ret, mpf_t current_tol, point_data_mp *in_mp, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int determineRankDef_amp(double *CN_ret, double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int));

void findSingVals_d(mat_d singVals, double *funcResid, double *condNum, point_data_d *in_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
void findSingVals_mp(mat_mp singVals, mpf_t funcResid, double *condNum, point_data_mp *in_mp, int curr_prec, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int jacobian_analyze(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, int MPType, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank:  == 0 means non-singular               *
* NOTES: verifies that the Jacobian is non-singular             *
\***************************************************************/
{
  int retVal;
  double temp_tol_d = 0;
  mpf_t temp_tol_mp;

  if (MPType == 0 || MPType == 1)
  { // find temp_tol_d & temp_tol_mp
    if (tol_prec > 52)
    { // tolerance input as higher precision
      mpf_init2(temp_tol_mp, tol_prec);
      mpf_set(temp_tol_mp, current_tol_mp);
      temp_tol_d = mpf_get_d(temp_tol_mp);
    }
    else
    { // tolerance input as double precision
      tol_prec = MAX(prec_in, 64); // atleast need 64-bit preicision
      mpf_init2(temp_tol_mp, tol_prec);
      mpf_set_d(temp_tol_mp, current_tol_d);
      temp_tol_d = current_tol_d;
      tol_prec = 52; // set back 
    }
  }

  if (MPType == 0)
  { // we only have double precision evaluators
    if (prec_in == 52)
    { // input as double precision
      retVal = jacobian_analyze_d(temp_tol_d, in_d, ED_d, eval_func_d);
    }
    else
    { // we need to convert in_mp to double precision since we can only evaluate in double precision
      point_data_d tempPD_d;
      init_point_data_d(&tempPD_d, in_mp->point->size);
      convert_point_data_mp_to_d(&tempPD_d, in_mp);

      retVal = jacobian_analyze_d(temp_tol_d, &tempPD_d, ED_d, eval_func_d);

      clear_point_data_d(&tempPD_d);
    }
  }
  else if (MPType == 1)
  { // we only have multiprecision evaluators

    if (prec_in > 52)
    { // input as multiprecision
      retVal = jacobian_analyze_mp(temp_tol_mp, in_mp, prec_in, ED_mp, eval_func_mp);
    }
    else
    { // we need to convert in_d to fixed multiprecision since we can only evaluate in multiprecision
      // find the current precision
      prec_in = (int) mpf_get_default_prec();
      point_data_mp tempPD_mp;
      init_point_data_mp2(&tempPD_mp, in_d->point->size, prec_in);
      convert_point_data_d_to_mp(&tempPD_mp, in_d);

      retVal = jacobian_analyze_mp(temp_tol_mp, &tempPD_mp, prec_in, ED_mp, eval_func_mp);

      clear_point_data_mp(&tempPD_mp);
    }
  }
  else 
  { // use AMP - both evaluators are known
    retVal = jacobian_analyze_amp(current_tol_d, current_tol_mp, tol_prec, in_d, in_mp, prec_in, ED_d, ED_mp, eval_func_d, eval_func_mp);
  }

  if (MPType == 0 || MPType == 1)
    mpf_clear(temp_tol_mp);

  return retVal;
}

int jacobian_analyze_d(double current_tol, point_data_d *in, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank:  == 0 means non-singular               *
* NOTES: analyzes the Jacobian                                  *
\***************************************************************/
{
  int retVal = 0, cont = 1, its = 10, svd_its = 100, init_corank = 0, jac_corank = 0, analyze_corank = 0;
  double residual, largeChange, tol_conv, tol_sign, percent_change;
  mat_d U, E1, E2, V;
  eval_struct_d e;
  point_d newPoint;
 
  init_mat_d(U, 0, 0); init_mat_d(E1, 0, 0); init_mat_d(E2, 0, 0); init_mat_d(V, 0, 0);
  init_eval_struct_d(e, 0, 0, 0);

  // copy to newPoint
  init_point_d(newPoint, in->point->size);
  point_cp_d(newPoint, in->point);

  // setup largeChange
  if (current_tol >= 1e-15) // current_tol cannot be better than the precision can handle
    largeChange = 0.1 / current_tol;
  else
  { // no tolerance information so set to a drop in 14 orders of magnitude for double precision
    current_tol = 1e-14;
    largeChange = 1e14; 
  }

  // setup tol_conv & tol_sign
  tol_conv = 1e-15; // tolerance for convergence in jacobi iterations
  tol_sign = 1e-20; // tolerance for finding the sign of a number

  // setup percent_change
  percent_change = 0.20; // a 20% change in singular value - implies that it is close to 0

  // find Jv
  eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, in->point, in->time, ED, eval_func);

  // find svd of Jv
  init_corank = svd_jacobi_d(U, E1, V, e.Jv, svd_its, current_tol, tol_conv, tol_conv, tol_sign, largeChange);

  if (init_corank < 0)
  { // error in finding svd - quit
    retVal = init_corank;
    cont = 0;
  }
  else
  { // continue on with analysis
    cont = 1;
  }

  while (cont)
  { // perform a newton iteration on newPoint
    retVal = newton_iteration_d(&residual, 0, NULL, NULL, newPoint, in->time, &e, ED, eval_func);

    if (retVal)
    { // problem with matrixSolve in newton iteration
      cont = 0;
    }
    else
    { // evaluate the new Jacobian
      eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, newPoint, in->time, ED, eval_func);

      // find svd of Jv
      jac_corank = svd_jacobi_d(U, E2, V, e.Jv, svd_its, current_tol, tol_conv, tol_conv, tol_sign, largeChange);

      // compare E1 & E2
      analyze_corank = relative_svd_change_d(percent_change, E1, E2);

      if (init_corank == jac_corank && jac_corank == analyze_corank)
      { // everythings agrees on corank
        retVal = init_corank;
        cont = 0;
      }
      else if (init_corank != jac_corank)
      { // update init_corank and try again
        init_corank = jac_corank;
      }
    }

    // fail safe
    if (cont > 0)
      cont++;
    if (cont > its)
    {
      if (jac_corank > 0 || analyze_corank > 0)
        retVal = MAX(jac_corank, analyze_corank); // return maximum corank when one of them is positive
      else
        retVal = MIN(jac_corank, analyze_corank); // retun the minimum if both are non-positive (-1 means failure)

      cont = 0;
    }
  }

  clear_mat_d(U); clear_mat_d(E1); clear_mat_d(E2); clear_mat_d(V);
  clear_eval_struct_d(e);
  clear_point_d(newPoint);

  return retVal;
}

int relative_svd_change_d(double percent_change, mat_d E1, mat_d E2)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank based on percent_change                 *
* NOTES: determines if there are changes in the diagonal entries* 
*           E1 & E2 based on percent_change                     *
\***************************************************************/
{
  int cont, corank, i, rows = E1->rows, cols = E1->cols;
  double absVal_diff, absVal;
  vec_d diffVec;

  // make sure E1 & E2 are the same size
  if (E1->rows != E2->rows || E1->cols != E2->cols)
  {
    printf("ERROR: The matrix sizes in relative_svd_change_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (percent_change <= 0)
  {
    printf("ERROR: The tolerances must be positive in relative_svd_change_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // put rows to be the minimum of rows & cols so that rows is in the number of diagonal entries to compare
  rows = MIN(rows, cols);

  // initialize cont & corank
  cont = 1;
  corank = rows;

  // initialize diffVec
  init_vec_d(diffVec, rows);

  // compare the difference & relative difference in the diagonal entries
  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < rows; i++)
  { // find difference
    sub_d(&diffVec->coord[i], &E1->entry[i][i], &E2->entry[i][i]);

    // find the magnitude of the difference
    absVal_diff = d_abs_d(&diffVec->coord[i]);

    // find the relative difference
    absVal = d_abs_d(&E1->entry[i][i]);
    if (absVal == 0)
    { // try E2[i][i] since E1[i][i] is 0
      absVal = d_abs_d(&E2->entry[i][i]);

      if (absVal != 0)
      { // E2[i][i] is not 0 - divide by it
        absVal = absVal_diff / absVal;
      }
      else
      { // both are 0 so this is okay!
        absVal = 0;
      }
    }
    else
    { // E1[i][i] is not 0 - divide by it
      absVal = absVal_diff / absVal;
    }

    // compare relative difference
    if (absVal > percent_change)
    { // too big of a relative change
      cont = 0;
      i = rows;
    }
    else
    { // this one is okay
      corank--;
    }
  }

  clear_vec_d(diffVec);

  return corank;
}

////// MULTIPRECISION //////

int jacobian_analyze_mp(mpf_t current_tol, point_data_mp *in, int prec_in, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank:  == 0 means non-singular               *
* NOTES: analyzes the Jacobian                                  *
\***************************************************************/
{
  int retVal = 0, cont = 1, its = 10, svd_its = 100, init_corank = 0, jac_corank = 0, analyze_corank = 0;
  int digits = (int) floor(prec_in * log10(2.0) - 0.5);
  char *str = NULL;
  size_t size;

  mpf_t residual, largeChange, tol_conv, tol_sign, percent_change;
  eval_struct_mp e;
  mat_mp U, E1, E2, V;
  point_mp newPoint;

  // initialize to this precision
  mpf_init2(residual, prec_in); mpf_init2(largeChange, prec_in); mpf_init2(tol_conv, prec_in); mpf_init2(tol_sign, prec_in); mpf_init2(percent_change, prec_in);
  init_eval_struct_mp2(e, 0, 0, 0, prec_in);
  init_mat_mp2(U, 0, 0, prec_in); init_mat_mp2(V, 0, 0, prec_in); init_mat_mp2(E1, 0, 0, prec_in); init_mat_mp2(E2, 0, 0, prec_in);
  init_point_mp2(newPoint, in->point->size, prec_in);

  // copy to newPoint
  point_cp_mp(newPoint, in->point);

  // find what the precision can handle
  size = 1 + snprintf(NULL, 0, "1e%d", -(digits - 2));
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", -(digits - 2));
  mpf_set_str(largeChange, str, 10);

  // setup largeChange
  if (mpf_cmp(current_tol, largeChange) >= 0) // current_tol cannot be better than the precision can handle
  {
    mpf_mul_ui(largeChange, current_tol, 10);
    mpf_ui_div(largeChange, 1, largeChange);
  }
  else
  { // no tolerance information so set to a drop of (digits - 2) orders of magnitude for this precision
    size = 1 + snprintf(NULL, 0, "1e%d", -(digits - 2));
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", -(digits - 3));
    mpf_set_str(current_tol, str, 10);

    size = 1 + snprintf(NULL, 0, "1e%d", digits - 2);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", digits - 3);
    mpf_set_str(largeChange, str, 10);
  }

  // set tol_conv & tol_sign
  size = 1 + snprintf(NULL, 0, "1e%d", -digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", -digits);
  mpf_set_str(tol_conv, str, 10); // tolerance for convergence in jacobi iterations
  size = 1 + snprintf(NULL, 0, "1e%d", -2 * (digits - 2)); 
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", -2 * (digits - 2));
  mpf_set_str(tol_sign, str, 10); // tolerance for finding the sign of a number

  // setup percent_change
  mpf_set_ui(percent_change, 1);
  mpf_div_ui(percent_change, percent_change, 5); // a 20% change in singular value - implies that it is close to 0

  // find Jv
  eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, in->point, in->time, ED, eval_func);

  // find svd of Jv
  init_corank = svd_jacobi_mp(U, E1, V, e.Jv, svd_its, current_tol, tol_conv, tol_conv, tol_sign, largeChange);

  if (init_corank < 0)
  { // svd error - quit
    retVal = init_corank;
    cont = 0;
  }
  else
  { // svd was okay so we need to continue
    cont = 1;
  }

  while (cont)
  {
    // perform a newton iteration on newPoint
    retVal = newton_iteration_mp(residual, 0, NULL, NULL, newPoint, in->time, &e, ED, eval_func);

    if (retVal)
    { // problem with matrixSolve in newton iteration - use higher precision???
      cont = 0;
    }
    else
    { 
      // evaluate the new Jacobian
      eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, newPoint, in->time, ED, eval_func);

      // find svd of Jv
      jac_corank = svd_jacobi_mp(U, E2, V, e.Jv, svd_its, current_tol, tol_conv, tol_conv, tol_sign, largeChange);

      // compare E1 & E2
      analyze_corank = relative_svd_change_mp(percent_change, E1, E2);

      if (init_corank == jac_corank && jac_corank == analyze_corank)
      { // everythings agrees on corank
        retVal = init_corank;
        cont = 0;
      }
      else if (init_corank != jac_corank)
      { // update init_corank and try again
        init_corank = jac_corank;
      }
    }

    // fail safe
    if (cont > 0)
      cont++;
    if (cont > its)
    {
      if (jac_corank > 0 || analyze_corank > 0)
        retVal = MAX(jac_corank, analyze_corank); // return maximum corank when one of them is positive
      else
        retVal = MIN(jac_corank, analyze_corank); // retun the minimum if both are non-positive (-1 means failure)

      cont = 0;
    }
  }

  // clear MP
  mpf_clear(residual); mpf_clear(largeChange); mpf_clear(tol_conv); mpf_clear(tol_sign); mpf_clear(percent_change);
  clear_eval_struct_mp(e);
  clear_mat_mp(U); clear_mat_mp(V); clear_mat_mp(E1); clear_mat_mp(E2);
  clear_point_mp(newPoint);

  // free str
  free(str);

  return retVal;
}

int relative_svd_change_mp(mpf_t percent_change, mat_mp E1, mat_mp E2)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank based on percent_change                 *
* NOTES: determines if there are changes in the diagonal entries*
*           E1 & E2 based on percent_change                     *
\***************************************************************/
{
  int cont, corank, i, rows = E1->rows, cols = E1->cols;
  mpf_t absVal_diff, absVal;
  vec_mp diffVec;

  // make sure E1 & E2 are the same size
  if (E1->rows != E2->rows || E1->cols != E2->cols)
  {
    printf("ERROR: The matrix sizes in relative_svd_change_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (mpf_cmp_ui(percent_change, 0) <= 0)
  {
    printf("ERROR: The tolerances must be positive in relative_svd_change_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // find the minimum precision
  i = mpf_get_prec(E1->entry[0][0].r);
  cont = mpf_get_prec(E2->entry[0][0].r);
  i = MIN(i, cont);

  // initialize MP to this precision
  mpf_init2(absVal_diff, i);
  mpf_init2(absVal, i);

  // put rows to be the minimum of rows & cols so that rows is in the number of diagonal entries to compare
  rows = MIN(rows, cols);

  // initialize cont & corank
  cont = 1;
  corank = rows;

  // setup diffVec
  init_vec_mp2(diffVec, rows, i);

  // compare the difference & relative difference in the diagonal entries
  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < rows; i++)
  { // find difference
    sub_mp(&diffVec->coord[i], &E1->entry[i][i], &E2->entry[i][i]);

    // find the magnitude of the difference
    mpf_abs_mp(absVal_diff, &diffVec->coord[i]);

    // find the relative difference
    mpf_abs_mp(absVal, &E1->entry[i][i]);
    if (mpfr_zero_p(absVal))
    { // try E2[i][i] since E1[i][i] is 0
      mpf_abs_mp(absVal, &E2->entry[i][i]);

      if (!mpfr_zero_p(absVal))
      { // E2[i][i] is not 0 - divide by it
        mpf_div(absVal, absVal_diff, absVal);
      }
      else
      { // both are 0 so this is okay!
        mpf_set_ui(absVal, 0);
      }
    }
    else
    { // E1[i][i] is not 0 - divide by it
      mpf_div(absVal, absVal_diff, absVal);
    }

    // compare relative difference
    if (mpf_cmp(absVal,percent_change) > 0)
    { // too big of a relative change
      cont = 0;
      i = rows;
    }
    else
    { // this one is okay
      corank--;
    }
  }

  // clear MP 
  mpf_clear(absVal_diff);
  mpf_clear(absVal);
  clear_vec_mp(diffVec);

  return corank;
}

////// ADAPTIVE MULTIPRECISION ///////

int jacobian_analyze_amp(double current_tol_d, mpf_t current_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank:  == 0 means non-singular               *
* NOTES: analyzes the Jacobian using adaptive precision         *
*  double precision is not checked since this does not cost much*
*  more in 64-bit precision and for ease of programming         *
\***************************************************************/
{
  int retVal = -1, curr_prec, tol_digits, its, max_its = 20, size;
  char *str = NULL;
  double temp_tol_d;
  comp_mp temp_tol_mp; // comp_mp used to have the ability to change precisions easily without losing the information
  
  // since we are changing precision, we need access to do it
  basic_eval_data_d *BED = (basic_eval_data_d *)ED_d;

  curr_prec = MAX(tol_prec, 64); // atleast need 64 bits
  init_mp2(temp_tol_mp, curr_prec);
  set_zero_mp(temp_tol_mp);  

  if (tol_prec == 52)
  { // use current_tol_d to get an idea as to how many correct digits there should be currently
    if (current_tol_d <= 0)
    { // the tolerance is 0 so we set the number of correct digits based on the current precision
      if (prec_in == 52)
        tol_digits = 14;
      else
        tol_digits = (int) floor(prec_in * log10(2) - 1.5);

      size = 1 + snprintf(NULL, 0, "1e%d", -tol_digits);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "1e%d", -tol_digits);
      mpf_set_str(temp_tol_mp->r, str, 10);

      temp_tol_d = mpf_get_d(temp_tol_mp->r);
    }
    else
    { // use current_tol_d as it is
      temp_tol_d = current_tol_d;
      mpf_set_d(temp_tol_mp->r, temp_tol_d);
      tol_digits = (int) floor(-log10(current_tol_d));
    }
  }
  else
  { // user current_tol_mp to get an idea as to how many correct digits there should be currently
    if (mpf_cmp_ui(current_tol_mp, 0) <= 0)
    { // the tolerance is 0 so we set the number of correct digits based on the current precision
      if (prec_in == 52)
        tol_digits = 14;
      else
        tol_digits = (int) floor(prec_in * log10(2) - 1.5);

      size = 1 + snprintf(NULL, 0, "1e%d", -tol_digits);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "1e%d", -tol_digits);
      mpf_set_str(temp_tol_mp->r, str, 10);

      temp_tol_d = mpf_get_d(temp_tol_mp->r);
    }
    else
    { // use current_tol_mp as it is
      mpf_set(temp_tol_mp->r, current_tol_mp);
      temp_tol_d = mpf_get_d(temp_tol_mp->r);

      mpfr_log10(temp_tol_mp->i, temp_tol_mp->r, __gmp_default_rounding_mode);
      mpf_neg(temp_tol_mp->i, temp_tol_mp->i);
      mpfr_floor(temp_tol_mp->i, temp_tol_mp->i);
      tol_digits = (int) mpf_get_si(temp_tol_mp->i);
    }
  }
  // convert digits to bits of precision
  if (tol_digits <= 14)
  { // double precision suffices
    tol_digits = 52;
  }
  else
  { // calculate the minimum precision required based on tol_digits
    tol_digits = (int) ceil((tol_digits + 0.5)/ (32 * log10(2)));
    if (tol_digits < 2)
    { // set to use 64-bit
      tol_digits = 64;
    }
    else
      tol_digits *= 32;
  }
  // so, tol_digits is the minimum precision required
  curr_prec = MAX(tol_digits, prec_in);

  // see if we can use double precision
  if (curr_prec == 52)
  { // check double precision analysis
    retVal = jacobian_analyze_d(temp_tol_d, in_d, ED_d, eval_func_d);
  }
  else
  { // setup to go into the MP loop
    retVal = -1;
  }  

  // check for success
  if (retVal < 0)
  { // we need to do the MP loop
    point_data_mp tempPD;

    curr_prec = MAX(curr_prec, 64); // need atleast 64-bits

    // set global things to curr_prec
    initMP(curr_prec);
    changeBasicEvalPrec_mp(curr_prec, BED->BED_mp);

    // initialize tempPD to this precision
    init_point_data_mp2(&tempPD, 0, curr_prec);
  
    // change temp_tol_mp to this precision
    change_prec_mp2(temp_tol_mp, curr_prec);

    // copy to tempPD
    if (prec_in == 52)
    {
      convert_point_data_d_to_mp(&tempPD, in_d);
    }
    else
    { 
      point_data_cp_mp(&tempPD, in_mp);
    }

    // loop over higher multiprecision until one gives a good answer or we took too many iterations
    its = 0;
    while (retVal < 0 && its < max_its)
    { // increment the number of iterations
      its++;

      // try this precision
      retVal = jacobian_analyze_mp(temp_tol_mp->r, &tempPD, curr_prec, ED_mp, eval_func_mp);

      if (retVal < 0)
      { // find new precision
        curr_prec += 32;

        if (its < max_its)
        { // increase precision so we can try again
          initMP(curr_prec);
          changeBasicEvalPrec_mp(curr_prec, BED->BED_mp);
          change_prec_point_data_mp(&tempPD, curr_prec);
          change_prec_mp2(temp_tol_mp, curr_prec);
        }
      }
    }
  
    if (prec_in == 52)
    { // set things back to 64-bits
      initMP(64);
      changeBasicEvalPrec_mp(64, BED->BED_mp);
    }
    else
    { // set things back to prec_in
      initMP(prec_in);
      changeBasicEvalPrec_mp(prec_in, BED->BED_mp);
    }
  
    // clear MP
    clear_point_data_mp(&tempPD);
  }

  // release memory
  free(str);
  clear_mp(temp_tol_mp);

  return retVal;
}

///////////////// USED TO DETERMINE RANK DEFICIENCY //////////////////

void findSingVals_d(mat_d singVals, double *funcResid, double *condNum, point_data_d *in_d, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the singular values of the jacobian, function    *
* residual and condition number                                 *
\***************************************************************/
{
  int retVal, size, its = 100, num_digits = -14; // -14 = -(int) floor(52 * log10(2.0) - 1.5)
  double rank_tol = pow(10, num_digits + 1), tol_prec = pow(10, num_digits), tol_sign = pow(10, 2 * num_digits), largeChange = pow(10, -(num_digits + 1));
  // rank_tol & largeChange are not really needed since we are only interested in the singular values
  eval_struct_d e;
  mat_d U, V;

  init_eval_struct_d(e, 0, 0, 0);
  init_mat_d(U, 0, 0);
  init_mat_d(V, 0, 0);

  // evaluate the function
  eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, in_d->point, in_d->time, ED, eval_func);
  // find the function residual
  *funcResid = infNormVec_d(e.funcVals);

  // perform the svd on Jv
  retVal = svd_jacobi_d(U, singVals, V, e.Jv, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);

  if (retVal < 0)
  { // SVD was not successful - adjust tolerances to see if we can get convergence - this should rarely, if ever, happen
    // if we don't have convergence, it still returns a good approximation of the SVD
    tol_prec *= 100;
    tol_sign /= 100;
    
    svd_jacobi_d(U, singVals, V, e.Jv, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);
  }

  // find the last singular value
  if (singVals->rows > singVals->cols)
    size = singVals->cols - 1;
  else
    size = singVals->rows - 1;

  // find the condition number
  if (singVals->entry[size][size].r == 0)
  { // report condition number as -1
    *condNum = -1;
  }
  else
  { // divide to get condition number
    *condNum = singVals->entry[0][0].r / singVals->entry[size][size].r;
  }

  clear_eval_struct_d(e);
  clear_mat_d(U);
  clear_mat_d(V);

  return;
}

void findSingVals_mp(mat_mp singVals, mpf_t funcResid, double *condNum, point_data_mp *in_mp, int curr_prec, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the singular values of the jacobian, function    *
* residual and condition number                                 *
\***************************************************************/
{
  int retVal, its = 100, num_digits = -(int) floor(curr_prec * log10(2.0) - 1.5); 
  size_t size;
  char *str = NULL;
  mpf_t rank_tol, tol_prec, tol_sign, largeChange;
  eval_struct_mp e;
  mat_mp U, V;

  // initialize to curr_prec
  mpf_init2(rank_tol, curr_prec); mpf_init2(tol_prec, curr_prec); mpf_init2(tol_sign, curr_prec); mpf_init2(largeChange, curr_prec);
  init_eval_struct_mp2(e, 0, 0, 0,curr_prec);
  init_mat_mp2(U, 0, 0, curr_prec); init_mat_mp2(V, 0, 0, curr_prec);

  // setup tolerances
  size = 1 + snprintf(NULL, 0, "1e%d", num_digits + 2);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", num_digits);
  mpf_set_str(rank_tol, str, 10); // not really used

  size = 1 + snprintf(NULL, 0, "1e%d", num_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", num_digits);
  mpf_set_str(tol_prec, str, 10);

  size = 1 + snprintf(NULL, 0, "1e%d", 2 * num_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", 2 * num_digits);
  mpf_set_str(tol_sign, str, 10);

  size = 1 + snprintf(NULL, 0, "1e%d", -(num_digits + 2));
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", -(num_digits + 2));
  mpf_set_str(largeChange, str, 10); // not really used

  // evaluate the function
  eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, in_mp->point, in_mp->time, ED, eval_func);

  // find the function residual
  infNormVec_mp2(funcResid, e.funcVals);

  // perform the svd on Jv
  retVal = svd_jacobi_mp(U, singVals, V, e.Jv, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);

  if (retVal < 0)
  { // SVD was not successful - adjust tolerances to see if we can get convergence - this should rarely, if ever, happen
    // if we don't have convergence, it still returns a good approximation of the SVD
    mpf_mul_ui(tol_prec, tol_prec, 100);
    mpf_div_ui(tol_sign, tol_sign, 100);

    svd_jacobi_mp(U, singVals, V, e.Jv, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);
  }

  // find the last singular value
  if (singVals->rows > singVals->cols)
    size = singVals->cols - 1;
  else
    size = singVals->rows - 1;

  // find the condition number
  if (mpfr_zero_p(singVals->entry[size][size].r))
  { // report condition number as -1
    *condNum = -1;
  }
  else
  { // divide to get condition number
    mpf_div(largeChange, singVals->entry[0][0].r, singVals->entry[size][size].r);
    *condNum = mpf_get_d(largeChange);
  }

  // clear memory
  mpf_clear(rank_tol); mpf_clear(tol_prec); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_eval_struct_mp(e);
  clear_mat_mp(U); clear_mat_mp(V);

  free(str);

  return;
}

int determineRankDef(double *CN, double rank_tol_d, mpf_t rank_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp,  int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: == 0 means non-singular, otherwise singular    *
*   Also returns condition number                               *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal;

  if (T->MPType == 0)
  { // setup to use the double precision function
    double temp_tol_d;
    point_data_d PD_d;
    init_point_data_d(&PD_d, 0);

    if (tol_prec > 52)
    { // tolerance is in rank_tol_mp
      temp_tol_d = mpf_get_d(rank_tol_mp);
    }
    else
    { // tolerance is in rank_tol_d
      temp_tol_d = rank_tol_d;
    }

    if (prec_in > 52)
    { // point in in_mp
      convert_point_data_mp_to_d(&PD_d, in_mp);
    }
    else
    { // point in in_d
      point_data_cp_d(&PD_d, in_d);
    }
    retVal = determineRankDef_d(CN, temp_tol_d, &PD_d, T, OUT, ED_d, eval_func_d);

    clear_point_data_d(&PD_d);
  }
  else if (T->MPType == 1)
  { // setup to use the fixed multiprecision function
    mpf_t temp_tol_mp;
    point_data_mp PD_mp;

    mpf_init2(temp_tol_mp, T->Precision);
    init_point_data_mp2(&PD_mp, 0, T->Precision);

    if (tol_prec > 52)
    { // tolerance is in rank_tol_mp
      mpf_set(temp_tol_mp, rank_tol_mp);
    }
    else
    { // tolerance is in rank_tol_d
      mpf_set_d(temp_tol_mp, rank_tol_d);
    } 

    if (prec_in > 52)
    { // point in in_mp
      point_data_cp_mp(&PD_mp, in_mp);
    }
    else
    { // point in in_d
      convert_point_data_d_to_mp(&PD_mp, in_d);
    }
    retVal = determineRankDef_mp(CN, temp_tol_mp, &PD_mp, T, OUT, ED_mp, eval_func_mp);

    clear_point_data_mp(&PD_mp);
    mpf_clear(temp_tol_mp);
  }
  else // MPType == 2
  { // call AMP function as is
    retVal = determineRankDef_amp(CN, rank_tol_d, rank_tol_mp, tol_prec, in_d, in_mp, prec_in, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  }

  return retVal;
}

int determineRankDef_d(double *CN_ret, double rank_tol, point_data_d *in_d, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: == 0 means non-singular, otherwise singular    *
*   Also returns condition number                               *
* NOTES: analyzes the Jacobian using double precision           *
\***************************************************************/
{
  int retVal = 0, copy_back = 1, its, maxIts = MAX(1 + T->maxNewtonIts, 2), func_digits_for_correction = 10, newton_digits_for_correction = 10;
  int min_digit_change = maxIts; // each iteration should provide more correct digits
  int max_sv_digit_change = 3; // the minimum singular might change a little on each iteration - but since we are increasing the precision at the beginning to find the
                               // minimum singular value to full precision, it should not change all that much
  double max_CN = 1e13; // max CN allowed based on precision
  double resid_tol = MAX(T->endgameNewtonTol, 1e-13); // the residual can not be smaller than the input precision
  point_data_d tempPD;
  eval_struct_d e;

  double initial_min_sv, final_min_sv, initial_log_sv, final_log_sv, initial_func_resid;

  mat_d *Jv = (mat_d *)bmalloc(maxIts * sizeof(mat_d));
  double *newton_resid = (double *)bmalloc(maxIts * sizeof(double)), *func_resid = (double *)bmalloc(maxIts * sizeof(double));
  int *newton_digits = (int *)bmalloc(maxIts * sizeof(int)), *func_digits = (int *)bmalloc(maxIts * sizeof(int));

  init_eval_struct_d(e, 0, 0, 0);
  for (its = 0; its < maxIts; its++)
    init_mat_d(Jv[its], 0, 0);

  // copy in_d to tempPD
  init_point_data_d(&tempPD, in_d->point->size);
  point_data_cp_d(&tempPD, in_d);

  // find the initial CN & function residual
  // evaluate the function
  eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, tempPD.point, tempPD.time, ED, eval_func);
  // find the minimum singular value
  approx_min_svd_d(&initial_min_sv, e.Jv);
  initial_func_resid = infNormVec_d(e.funcVals);

  // find log10(1/min_sv)
  initial_log_sv = -log10(initial_min_sv);

  // do the main loop of newton iterations
  retVal = 0;
  for (its = 0; its < maxIts && retVal == 0; its++)
  { // initialize values
    newton_resid[its] = func_resid[its] = 0;
    // do a newton iteration
    retVal = newton_iteration_d(&newton_resid[its], 0, NULL, NULL, tempPD.point, tempPD.time, &e, ED, eval_func);
    // evaluate the function to find the new Jacobian and funcVals
    eval_d(e.funcVals, e.parVals, e.parDer, Jv[its], e.Jp, tempPD.point, tempPD.time, ED, eval_func);
    // find the function residual
    func_resid[its] = infNormVec_d(e.funcVals);    
    fprintf(OUT, "newton: %e func: %e\n", newton_resid[its], func_resid[its]);

    if (retVal == 0)
    { // look at function residuals
      if (func_resid[its] > resid_tol)
        retVal = 1;
      else if (func_resid[its] == 0)
      { // function residual is zero - set digits to the number of digits of precision
        func_digits[its] = 15;
      }
      else
      { // calculate the number of digits in the function residual
        func_digits[its] = (int) ceil(-log10(func_resid[its]));
      }
    }

    if (retVal == 0)
    { // look at newton residuals
      if (newton_resid[its] > resid_tol)
        retVal = 1;
      else if (newton_resid[its] == 0)
      { // newton residual is zero - set digits to the number of digits of precision
        newton_digits[its] = 15;
      }
      else
      { // calculate the number of digits correct based on newton residual
        newton_digits[its] = (int) ceil(-log10(newton_resid[its]));
      }
    }

    // check to see if we have refined to full precision
    if (its > 1 && retVal == 0 && func_digits[its] >= func_digits_for_correction && newton_digits[its] >= newton_digits_for_correction)
    { // increment its and then break loop since this has been refined to full precision
      its++;
      break;
    }
  }
  // store the number of iterations completed
  maxIts = its;

  // check for newton iteration error
  if (retVal)
    retVal = 1;

  // make sure the function digits correct is atleast min_digit_change more than the first digits correct OR we have more than func_digits_for_correction
//  if (retVal == 0 && func_digits[0] + min_digit_change > func_digits[maxIts - 1] && func_digits[maxIts - 1] < func_digits_for_correction)
//    retVal = 1;

  // make sure the final digits correct is atleast min_digit_change more than the first digits correct OR we have more than newton_digits_for_correction
  if (retVal == 0 && newton_digits[0] + min_digit_change > newton_digits[maxIts - 1] && newton_digits[maxIts - 1] < newton_digits_for_correction)
    retVal = 1;

  // calculate CN and determine if minimum singular value has not changed too much
  if (retVal == 0)
  { // find minimum singular value for the last iteration
    approx_min_svd_d(&final_min_sv, Jv[maxIts - 1]);

    // setup CN_ret
    *CN_ret = infNormMat_d(Jv[maxIts - 1]) / final_min_sv;
    if (*CN_ret < 1)
      *CN_ret = 1;

    // find the log(1 / minimum singluar value)
    final_log_sv = -log10(final_min_sv);
    fprintf(OUT, "min sv: %e\n", final_min_sv);

    // determine if the minimum singular value has changed too much
    if (fabs(initial_log_sv - final_log_sv) > max_sv_digit_change || *CN_ret > max_CN)
      retVal = 1;
  }
  else
  { // calculate CN_ret based on initial minimum singular value
    *CN_ret = infNormMat_d(Jv[maxIts - 1]) / initial_min_sv;
    if (*CN_ret < 1)
      *CN_ret = 1;

    fprintf(OUT, "min sv: %e\n", initial_min_sv);
  }

  // check to see if we should copy tempPD back to in_mp
  copy_back = 1;
  if (retVal == 0 || (maxIts > 1 && initial_func_resid > func_resid[maxIts - 1]))
  {
    for (its = 0; its < maxIts && copy_back; its++)
      if (newton_resid[its] > resid_tol)
        copy_back = 0;

    if (copy_back)
    { // copy back to in_d
      point_data_cp_d(in_d, &tempPD);
    }
  }

  // clear the memory
  clear_eval_struct_d(e);
  for (its = 0; its < maxIts; its++)
    clear_mat_d(Jv[its]);
  clear_point_data_d(&tempPD);
  free(Jv);
  free(newton_resid);
  free(func_resid);
  free(newton_digits);
  free(func_digits);

  return retVal;
}

int determineRankDef_mp(double *CN_ret, mpf_t rank_tol, point_data_mp *in_mp, tracker_config_t *T, FILE *OUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: == 0 means non-singular, otherwise singular    *
*   Also returns condition number                               *
* NOTES: analyzes the Jacobian using fixed multiprecision       *
\***************************************************************/
{
  int retVal = 0, copy_back = 1, its, maxIts = MAX(1 + T->maxNewtonIts, 2), func_digits_for_correction = (int) floor(T->Precision * log10(2.0) - 4.5), newton_digits_for_correction = (int) floor(T->Precision * log10(2.0) - 4.5);
  int min_digit_change = maxIts; // each iteration should provide more correct digits
  int max_sv_digit_change = 3; // the minimum singular might change a little on each iteration - but since we are increasing the precision at the beginning to find the
                               // minimum singular value to full precision, it should not change all that much
  double max_CN = func_digits_for_correction < 300 ? pow(10, func_digits_for_correction) : 1e300; // max CN allowed based on precision
  double resid_tol = newton_digits_for_correction < 300 ? pow(10, -newton_digits_for_correction + 2) : 1e-300;
  resid_tol = MAX(resid_tol, T->endgameNewtonTol); // the residual can not be smaller than the input precision
  mpf_t tempMPF;
  point_data_mp tempPD;
  eval_struct_mp e;

  double min_sv_change;
  mpf_t initial_func_resid, initial_min_sv, final_min_sv;

  mpf_t *newton_resid = (mpf_t *)bmalloc(maxIts * sizeof(mpf_t)), *func_resid = (mpf_t *)bmalloc(maxIts * sizeof(mpf_t));
  int *newton_digits = (int *)bmalloc(maxIts * sizeof(int)), *func_digits = (int *)bmalloc(maxIts * sizeof(int));

  mpf_init(tempMPF);
  init_eval_struct_mp(e, 0, 0, 0);

  // copy in_mp to tempPD
  init_point_data_mp(&tempPD, in_mp->point->size);
  point_data_cp_mp(&tempPD, in_mp);

  // find the initial CN & function residual
  mpf_init(initial_func_resid); mpf_init(initial_min_sv); mpf_init(final_min_sv);
  // evaluate the function
  eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, tempPD.point, tempPD.time, ED, eval_func);
  // find the minimum singular value
  approx_min_svd_mp_prec(initial_min_sv, e.Jv, T->Precision);
  infNormVec_mp2(initial_func_resid, e.funcVals);

  // do the main loop of newton iterations
  retVal = 0;
  for (its = 0; its < maxIts && retVal == 0; its++)
  { // initialize
    mpf_init(func_resid[its]); mpf_set_ui(func_resid[its], 0);
    mpf_init(newton_resid[its]); mpf_set_ui(newton_resid[its], 0);

    // perform the iteration steps
    retVal = rankDef_newton_mp(newton_resid[its], &newton_digits[its], func_resid[its], &func_digits[its], &tempPD, T->Precision, T, &e, ED, eval_func);
    fprintf(OUT, "newton: %e func: %e\n", mpf_get_d(newton_resid[its]), mpf_get_d(func_resid[its]));

    // check to see if newton_resid or func_resid is too big
    if (retVal == 0 && (mpf_get_d(newton_resid[its]) > resid_tol || mpf_get_d(func_resid[its]) > resid_tol))
      retVal = 1;

    // check to see if we have refined to full precision
    if (its > 0 && retVal == 0 && func_digits[its] >= func_digits_for_correction && newton_digits[its] >= newton_digits_for_correction)
    { // increment its and then break loop since this has been refined to full precision
      its++;
      break;
    }
  }
  // store the number of iterations completed
  maxIts = its;

  // check for newton iteration error -> singular
  if (retVal)  
    retVal = 1;

  // make sure the function digits correct is atleast min_digit_change more than the first digits correct OR we have more than num_digits_for_correction
//  if (retVal == 0 && func_digits[0] + min_digit_change > func_digits[maxIts - 1] && func_digits[maxIts - 1] < func_digits_for_correction)
//    retVal = 1;

  // make sure the final digits correct is atleast min_digit_change more than the first digits correct OR we have more than num_digits_for_correction
  if (retVal == 0 && newton_digits[0] + min_digit_change > newton_digits[maxIts - 1] && newton_digits[maxIts - 1] < newton_digits_for_correction)
    retVal = 1;

  // calculate CN and determine if minimum singular value has not changed too much
  if (retVal == 0)
  { // setup CN_ret
    approx_min_svd_mp_prec(final_min_sv, e.Jv, T->Precision);

    fprintf(OUT, "min sv: %e\n", mpf_get_d(final_min_sv));
    *CN_ret = infNormMat_mp(e.Jv) / mpf_get_d(final_min_sv);

    if (*CN_ret < 1)
      *CN_ret = 1;

    // determine the differences in the the minimum singular values
    mpfr_log10(initial_min_sv, initial_min_sv, __gmp_default_rounding_mode);
    mpfr_log10(final_min_sv, final_min_sv, __gmp_default_rounding_mode);
    mpf_neg(initial_min_sv, initial_min_sv);
    mpf_neg(final_min_sv, final_min_sv);
    mpf_sub(final_min_sv, final_min_sv, initial_min_sv);
    mpf_abs(final_min_sv, final_min_sv);
    min_sv_change = mpf_get_d(final_min_sv);

    // determine if the minimum singular value has changed too much
    if (min_sv_change >= max_sv_digit_change || *CN_ret > max_CN)
      retVal = 1;
  }
  else
  { // calculate CN_ret based on initial minimum singular value
    *CN_ret = infNormMat_mp(e.Jv) / mpf_get_d(initial_min_sv);
    if (*CN_ret < 1)
      *CN_ret = 1;
    fprintf(OUT, "min sv: %e\n", mpf_get_d(initial_min_sv));
  }

  // check to see if we should copy tempPD back to in_mp
  copy_back = 1;
  if (retVal == 0 || (maxIts > 1 && mpf_cmp(initial_func_resid, func_resid[maxIts - 1]) > 0))
  {
    for (its = 0; its < maxIts && copy_back; its++)
      if (mpf_get_d(newton_resid[its]) > resid_tol)
        copy_back = 0;

    if (copy_back)
    { // copy back to in_mp
      point_data_cp_mp(in_mp, &tempPD);
    }
  }

  // clear the memory
  mpf_clear(initial_func_resid); mpf_clear(initial_min_sv); mpf_clear(final_min_sv);
  for (its = maxIts - 1; its >= 0; its--)
  {
    mpf_clear(func_resid[its]);
    mpf_clear(newton_resid[its]);
  }
  free(func_resid);
  free(newton_resid);
  free(newton_digits);
  free(func_digits);
 
  mpf_clear(tempMPF); 
  clear_point_data_mp(&tempPD);
  clear_eval_struct_mp(e);

  return retVal;
}

int determineRankDef_amp(double *CN_ret, double rank_tol_d, mpf_t rank_tol_mp, int tol_prec, point_data_d *in_d, point_data_mp *in_mp, int prec_in, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: == 0 means non-singular, otherwise singular    *
*   Also returns condition number                               *
* NOTES: analyzes the Jacobian using adaptive precision         *
\***************************************************************/
{
  int retVal = 0, copy_back = 1, its, maxIts = MAX(1 + T->maxNewtonIts, 2), curr_prec, tempInt;
  int min_digit_change = 3 * maxIts; // each iteration should provide more correct digits
  int max_sv_digit_change = 3; // the minimum singular might change a little on each iteration - but since we are increasing the precision at the beginning to find the
                               // minimum singular value to full precision, it should not change all that much
  double resid_tol = 200 * T->final_tol_times_mult;
  mpf_t tempMPF;
  point_data_mp tempPD_mp;
  eval_struct_mp e;

  double min_sv_change;
  mpf_t initial_min_sv, final_min_sv, initial_func_resid;

  mpf_t *newton_resid = (mpf_t *)bmalloc(maxIts * sizeof(mpf_t)), *func_resid = (mpf_t *)bmalloc(maxIts * sizeof(mpf_t));
  int *newton_digits = (int *)bmalloc(maxIts * sizeof(int)), *func_digits = (int *)bmalloc(maxIts * sizeof(int));

  curr_prec = MAX(tol_prec, 64); // atleast need 64 bits
  curr_prec = MAX(curr_prec, prec_in);
  init_eval_struct_mp2(e, 0, 0, 0, curr_prec);
  init_point_data_mp2(&tempPD_mp, 0, curr_prec);

  // setup over to tempPD_mp
  if (prec_in == 52)
  { // copy in_d to tempPD_mp
    convert_point_data_d_to_mp(&tempPD_mp, in_d);
  }
  else
  { // copy in_mp to tempPD_mp
    point_data_cp_mp(&tempPD_mp, in_mp);
  }

  // find the precision that is needed based on the smallest singular value
  mpf_init2(initial_min_sv, curr_prec); mpf_init2(initial_func_resid, curr_prec); mpf_init2(final_min_sv, curr_prec);
  curr_prec = prec_needed(T->MPType, initial_min_sv, initial_func_resid, &tempPD_mp, curr_prec, &e, ED_mp, eval_func_mp, change_prec);

  // set everything to this precision
  T->Precision = curr_prec;
  initMP(curr_prec);
  change_prec(ED_mp, curr_prec);
  setprec_eval_struct_mp(e, curr_prec);
  setprec_point_data_mp(&tempPD_mp, curr_prec);
  mpf_init2(tempMPF, curr_prec);

  // setup tempPD_mp with the new precision
  if (prec_in == 52)
  { // copy in_d to tempPD_mp
    convert_point_data_d_to_mp(&tempPD_mp, in_d);
  }
  else
  { // copy in_mp to tempPD_mp
    point_data_cp_mp(&tempPD_mp, in_mp);
  }

  // setup resid_tol
  if (prec_in == 52 && resid_tol < 1e-12)
  { // the residual can not be smaller than the input precision
    resid_tol = 1e-12;
  }
  else
  { // make sure the residual is not smaller than the input precision
    double tempD;
    tempInt = prec_to_digits(prec_in);
    tempD = tempInt < 300 ? pow(10, -tempInt + 6) : 1e-300;
    if (resid_tol < tempD)
    { // the residual can not be smaller than the input precision
      resid_tol = tempD;
    }
  }

  // do the main loop of newton iterations
  retVal = 0;
  for (its = 0; its < maxIts && retVal == 0; its++)
  { // increase precision
    T->Precision = curr_prec = 32 + curr_prec;

    initMP(curr_prec);
    change_prec(ED_mp, curr_prec);
    setprec_eval_struct_mp(e, curr_prec);
    change_prec_point_data_mp(&tempPD_mp, curr_prec);

    // initialize to this precision
    mpf_init2(func_resid[its], curr_prec);
    mpf_init2(newton_resid[its], curr_prec);

    // perform the iteration step
    retVal = rankDef_newton_mp(newton_resid[its], &newton_digits[its], func_resid[its], &func_digits[its], &tempPD_mp, curr_prec, T, &e, ED_mp, eval_func_mp);
    fprintf(OUT, "prec: %d newton: %e func: %e\n", curr_prec, mpf_get_d(newton_resid[its]), mpf_get_d(func_resid[its]));

    // check to see if newton_resid or func_resid is too big
    if (retVal == 0 && (mpf_get_d(newton_resid[its]) > resid_tol || mpf_get_d(func_resid[its]) > resid_tol))
      retVal = 1; 
  }
  // store the number of iterations completed
  maxIts = its;

  // check for newton iteration error -> singular
  if (retVal)
    retVal = 1;

  // make sure the function digits correct is atleast min_digit_change more than the first digits correct
  if (retVal == 0 && func_digits[0] + min_digit_change > func_digits[maxIts - 1] && !mpfr_zero_p(func_resid[maxIts - 1]))
    retVal = 1;

  // make sure the final digits correct is atleast min_digit_change more than the first digits correct
  if (retVal == 0 && newton_digits[0] + min_digit_change > newton_digits[maxIts - 1] && !mpfr_zero_p(newton_resid[maxIts - 1]))
    retVal = 1;

  // calculate CN and determine if minimum singular value has not changed too much
  if (retVal == 0)
  { // setup CN_ret
    mpf_set_prec(final_min_sv, curr_prec);
    approx_min_svd_mp_prec(final_min_sv, e.Jv, curr_prec);

    fprintf(OUT, "min sv: %e\n", mpf_get_d(final_min_sv));
    *CN_ret = infNormMat_mp(e.Jv) / mpf_get_d(final_min_sv);
    if (*CN_ret < 1)
      *CN_ret = 1;

    // determine the differences in the the minimum singular values
    mpfr_log10(initial_min_sv, initial_min_sv, __gmp_default_rounding_mode);
    mpfr_log10(final_min_sv, final_min_sv, __gmp_default_rounding_mode);
    mpf_neg(initial_min_sv, initial_min_sv);
    mpf_neg(final_min_sv, final_min_sv);
    mpf_sub(final_min_sv, final_min_sv, initial_min_sv);
    mpf_abs(final_min_sv, final_min_sv); 
    min_sv_change = mpf_get_d(final_min_sv);

    // determine if the minimum singular value has changed too much
    if (min_sv_change >= max_sv_digit_change)
      retVal = 1;
  }
  else
  { // calculate CN_ret based on initial minimum singular value
    *CN_ret = infNormMat_mp(e.Jv) / mpf_get_d(initial_min_sv);
    if (*CN_ret < 1)
      *CN_ret = 1;
    fprintf(OUT, "min sv: %e\n", mpf_get_d(initial_min_sv));
  }

  // check to see if we should copy tempPD_mp back to in_d or in_mp
  copy_back = 1;
  if (retVal == 0 || (maxIts > 1 && mpf_cmp(initial_func_resid, func_resid[maxIts - 1]) > 0))
  {
    for (its = 0; its < maxIts && copy_back; its++)
      if (mpf_get_d(newton_resid[its]) > resid_tol)
        copy_back = 0;

    if (copy_back)
    { // copy back to the correct location
      if (prec_in == 52)
      {
        convert_point_data_mp_to_d(in_d, &tempPD_mp);
      }
      else
      {
        point_data_cp_mp(in_mp, &tempPD_mp);
      }
    }
  }

  // clear the memory
  mpf_clear(initial_min_sv);
  mpf_clear(final_min_sv);
  mpf_clear(initial_func_resid);
  for (its = maxIts - 1; its >= 0; its--)
  {
    mpf_clear(func_resid[its]);
    mpf_clear(newton_resid[its]);
  }
  free(func_resid);
  free(newton_resid);
  free(newton_digits);
  free(func_digits);
  
  mpf_clear(tempMPF);
  clear_point_data_mp(&tempPD_mp);
  clear_eval_struct_mp(e);

  return retVal;
}

int prec_needed(int MPType, mpf_t min_sv, mpf_t func_resid, point_data_mp *in_mp, int prec_in, eval_struct_mp *e, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the precision needed to find see the minimum   *
*  singular value as a non-zero number                          *
* NOTES:                                                        *
\***************************************************************/
{ 
  int its = 0, maxIts = 10, cont = 1, size, curr_prec = prec_in, curr_digits = prec_to_digits(prec_in);
  point_data_mp tempPD;
  char *str = NULL;
  mpf_t min_sv_allowed;

  init_point_data_mp(&tempPD, 0);
  mpf_init(min_sv_allowed);

  do
  { 
    if (MPType == 2)
    { // make sure everything is set to the correct precision
      initMP(curr_prec);
      change_prec(ED_mp, curr_prec);
      setprec_eval_struct_mp(*e, curr_prec);
      setprec_point_mp(tempPD.point, curr_prec);
      setprec_mp(tempPD.time, curr_prec);
      mpf_set_prec(min_sv_allowed, curr_prec);

      mpf_set_prec(min_sv, curr_prec);
      mpf_set_prec(func_resid, curr_prec);
    }

    // setup min_sv_allowed
    size = 1 + snprintf(NULL, 0, "1e-%d", curr_digits - 3);
    str = (char *)bmalloc(size * sizeof(char));
    sprintf(str, "1e-%d", curr_digits - 3);
    mpf_set_str(min_sv_allowed, str, 10);

    // copy in_mp to tempPD - tempPD changes precision while in_mp does not!
    point_cp_mp(tempPD.point, in_mp->point);
    set_mp(tempPD.time, in_mp->time);

    // evaluate the function
    eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, tempPD.point, tempPD.time, ED_mp, eval_func_mp);

    // find the minimum singular value
    approx_min_svd_mp_prec(min_sv, e->Jv, curr_prec);

    // determine if min_sv is too small compared to the current precision
    if (MPType == 2 && mpf_cmp(min_sv, min_sv_allowed) < 0)
    { // min_sv is still too small - increase precision and try again
      curr_prec += 32;
      curr_digits = prec_to_digits(curr_prec);
      cont = 1;
    }
    else
    { // find the function residual and quit
      infNormVec_mp2(func_resid, e->funcVals);
      cont = 0;
    }

    // check to see if we have done too many iterations
    its++;
    if (cont && its > maxIts)
    { // find the function residual and quit since we have done too many iterations
      infNormVec_mp2(func_resid, e->funcVals);
      cont = 0;
    }

  } while (cont);

  // clear MP
  clear_point_data_mp(&tempPD);
  mpf_clear(min_sv_allowed);
 
  free(str);

  return curr_prec;
}

int rankDef_newton_mp(mpf_t newton_resid, int *newton_digits, mpf_t func_resid, int *func_digits, point_data_mp *in_mp, int curr_prec, tracker_config_t *T, eval_struct_mp *e, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: == 0 means iteration was successful            *
* NOTES: does a newton iteration and finds newton/function resid* 
\***************************************************************/
{ 
  int retVal;
  mpf_t tempMPF;

  mpf_init2(tempMPF, curr_prec);

  // do newton iteration
  retVal = newton_iteration_mp(newton_resid, 0, NULL, NULL, in_mp->point, in_mp->time, e, ED_mp, eval_func_mp);

  // check for successful newton iteration
  if (retVal == 0)
  { // evaluate the function to find the new Jacobian and funcVals
    eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, in_mp->point, in_mp->time, ED_mp, eval_func_mp);

    // find function residual
    infNormVec_mp2(func_resid, e->funcVals);

    // find newton_digits
    if (mpfr_zero_p(newton_resid))
    { // minimum singular value is zero - set digits to the number of digits of precision
      *newton_digits = prec_to_digits(curr_prec);
    }
    else
    { // calculate the number of digits in the function residual
      mpfr_log10(tempMPF, newton_resid, __gmp_default_rounding_mode);
      mpf_neg(tempMPF, tempMPF);
      mpfr_floor(tempMPF, tempMPF);
      *newton_digits = (int) mpf_get_si(tempMPF);

      // make sure that this is not more than the current precision should be able to handle
      *newton_digits = MIN(*newton_digits, prec_to_digits(curr_prec));
    }

    // find func_digits
    if (mpfr_zero_p(func_resid))
    { // minimum singular value is zero - set digits to the number of digits of precision
      *func_digits = prec_to_digits(curr_prec);
    }
    else
    { // calculate the number of digits in the function residual
      mpfr_log10(tempMPF, func_resid, __gmp_default_rounding_mode);
      mpf_neg(tempMPF, tempMPF);
      mpfr_floor(tempMPF, tempMPF);
      *func_digits = (int) mpf_get_si(tempMPF);

      // make sure that this is not more than the current precision should be able to handle
      *func_digits = MIN(*func_digits, prec_to_digits(curr_prec));
    }
  }

  mpf_clear(tempMPF);

  return retVal;
}

