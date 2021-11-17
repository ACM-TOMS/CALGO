// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

int QLP_L_d(mat_d L, mat_d A, double tol_pivot, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: Takes a matrix and produces a QLP^T decomposition      *
* where L is lower triangular matrix and Q & P would be unitary *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES:                                                        *
\***************************************************************/
{ 
  int rank, m = A->rows, n = A->cols;
  mat_d tempMat, tempR;

  // error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QLP_L_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QLP_L_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  if (m == 0 || n == 0)
  { // return
    change_size_mat_d(L, m, n);
    mat_cp_d(L, A);
    return 0;
  }

  // so we have m >= n >= 1
  init_mat_d(tempR, n, n);
  init_mat_d(tempMat, m, n);

  // FIRST: do a QR decomposition for A
  QR_R_d(tempR, tempMat, A, tol_pivot, tol_sign, largeChange, 0); 

  // now that we have an upper triangular matrix, transpose it and do a pivoted QR decomposition
  transpose_d(L, tempR);

  // remove the extra columns, if any
  L->cols = L->rows;

  // SECOND: do a pivoted QR decomposition of L 
  rank = QR_R_d(tempR, tempMat, L, tol_pivot, tol_sign, largeChange, 0);

  // setup L
  transpose_d(L, tempR);

  // free the memory
  clear_mat_d(tempR); clear_mat_d(tempMat);

  return (n - rank);
}

/////// multi precision ////////

int QLP_L_mp_prec(mat_mp L, mat_mp A, int curr_prec)
/***************************************************************\
* USAGE: Takes a matrix and produces a QLP^T decomposition      *
* where L is lower triangular matrix and Q & P would be unitary *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, num_digits = prec_to_digits(curr_prec);
  size_t size;
  char *str = NULL;
  mpf_t tol_pivot, tol_sign, largeChange;

  mpf_init2(tol_pivot, curr_prec); 
  mpf_init2(tol_sign, curr_prec);
  mpf_init2(largeChange, curr_prec);

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

  retVal = QLP_L_mp(L, A, tol_pivot, tol_sign, largeChange);

  mpf_clear(tol_pivot);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);
  free(str);

  return retVal;
}

int QLP_L_mp(mat_mp L, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: Takes a matrix and produces a QLP^T decomposition      *
* where L is lower triangular matrix and Q & P would be unitary *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES:                                                        *
\***************************************************************/
{
  int rank, m = A->rows, n = A->cols;
  mat_mp tempMat, tempR;

  // error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QLP_L_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (mpf_cmp_ui(tol_sign, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QLP_L_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  if (m == 0 || n == 0)
  { // return
    change_size_mat_mp(L, m, n);
    mat_cp_mp(L, A);
    return 0;
  }

  // so we have m >= n >= 1
  init_mat_mp(tempR, n, n);
  init_mat_mp(tempMat, m, n);

  // FIRST: do a QR decomposition for A
  QR_R_mp(tempR, tempMat, A, tol_pivot, tol_sign, largeChange, 0);

  // now that we have an upper triangular matrix, transpose it and do a pivoted QR decomposition
  transpose_mp(L, tempR);

  // remove the extra columns, if any
  L->cols = L->rows;

  // SECOND: do a pivoted QR decomposition of L
  rank = QR_R_mp(tempR, tempMat, L, tol_pivot, tol_sign, largeChange, 0);

  // setup L
  transpose_mp(L, tempR);

  // free the memory
  clear_mat_mp(tempR); clear_mat_mp(tempMat);

  return (n - rank);
}

