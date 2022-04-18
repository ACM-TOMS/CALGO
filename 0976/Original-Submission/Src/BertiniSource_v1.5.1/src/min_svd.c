// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// Find the minimum singular value of the given matrix in an efficient manner

///// DOUBLE PRECISION /////

void approx_min_svd_d(double *min_sv, mat_d A)
/***************************************************************\
* USAGE: approximates the minimum singular value of A using QLP *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int n;
  double tol_sign = 1e-20, tol_pivot = 1e-14, largeChange = 1e12;
  mat_d L;

  init_mat_d(L, 0, 0);

  // compute a QLP decomposition
  QLP_L_d(L, A, tol_pivot, tol_sign, largeChange);

  // min_sv is approximated by last diagonal entry
  n = MIN(L->rows, L->cols);
  if (n > 0)
  {
    *min_sv = d_abs_d(&L->entry[n-1][n-1]);
  }
  else
  {
    *min_sv = 0;
  }

  clear_mat_d(L);

  return;
}

int min_svd_d(double *min_sv, mat_d A, double error_tol)
/***************************************************************\
* USAGE: finds the minimum singular value of A using Li/Zeng    *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if convergence to error_tol, -1 otherwise    *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, its = 0, maxIts = 50, retVal = -1, m = A->rows, n = A->cols;
  double norm, t, currSV = 1e300, prevSV = 1e300, tol_pivot = MIN(1e-14, error_tol / 100), tol_sign = 1e-20, largeChange = 1e13;
  vec_d x, z, b;
  mat_d P, Q, R;

  // make sure m > 0 && n > 0
  if (m == 0 || n == 0)
  { // a dimension is 0
    *min_sv = 0;
    retVal = 0;
    return retVal;
  }

  init_vec_d(x, 0); init_vec_d(z, 0); init_vec_d(b, 0);
  init_mat_d(P, 0, 0); init_mat_d(Q, 0, 0); init_mat_d(R, 0, 0);

  // need m >= n
  if (m < n)
  { // work with the transpose
    mat_d tempMat;
    init_mat_d(tempMat, 0, 0);

    transpose_d(tempMat, A);
    m = tempMat->rows;
    n = tempMat->cols;
    
    // find QR-decomposition of A': tempMat * P = Q * R
    QR_d(Q, R, P, tempMat, tol_pivot, tol_sign, largeChange, 1);
  
    clear_mat_d(tempMat);
  }
  else
  { // find QR-decomposition of A: A * P = Q * R
    QR_d(Q, R, P, A, tol_pivot, tol_sign, largeChange, 1);
  }
  // so we have R being an m x n upper triangular matrix with m >= n >= 1
  // since m >= n & upper triangular, the last (m-n) rows are 0 - so we can treat R as n x n upper triangular
  R->rows = R->cols = n;

  // find the norm of R[0][0]
  t = d_abs_d(&R->entry[0][0]);

  // normalize R - an n x n upper triangular matrix
  if (t > 1)
  {
    norm = 1 / t;
    for (i = 0; i < n; i++)
      for (j = i; j < n; j++)
      {
        mul_rdouble_d(&R->entry[i][j], &R->entry[i][j], norm);
      }
  }

  // generate a random unit vector x
  increase_size_vec_d(x, n);
  x->size = n;
  norm = 0;
  for (i = 0; i < n; i++)
  {
    get_comp_rand_d(&x->coord[i]);
    norm += x->coord[i].r * x->coord[i].r + x->coord[i].i * x->coord[i].i;
  }
  norm = 1 / sqrt(norm);
  for (i = 0; i < n; i++)
  {
    mul_rdouble_d(&x->coord[i], &x->coord[i], norm);
  }

  // setup the size of Q & b
  increase_size_mat_d(Q, n + 1, n);
  increase_size_vec_d(b, n + 1);
  Q->rows = b->size = n + 1;
  Q->cols = n;

  // setup the bottom of Q - always == R
  for (i = 1; i <= n; i++)
    for (j = 0; j < n; j++)
    {
      set_d(&Q->entry[i][j], &R->entry[i-1][j]);
    }

  // main loop - tau is not used since we normalized R
  while (its < maxIts)
  { // setup Q & b
    for (i = 0; i <= n; i++)
      if (i == 0)
      { // setup top row of Q & find x' * x
        norm = 0;
        for (j = 0; j < n; j++)
        { // notice the conjugation
          Q->entry[i][j].r = 2 * x->coord[j].r;
          Q->entry[i][j].i = -2 * x->coord[j].i;
          norm += x->coord[j].r * x->coord[j].r + x->coord[j].i * x->coord[j].i;
        }
        
        // setup top entry in b
        set_double_d(&b->coord[0], norm - 1, 0);
      }
      else
      { // find ith entry of b = (i-1)st entry of R * x
        set_zero_d(&b->coord[i]);
        for (j = 0; j < n; j++)
        {
          sum_mul_d(&b->coord[i], &R->entry[i-1][j], &x->coord[j]);
        }
      }

    // find the least squares solution z to Q*z = b
    i = matrixSolve_Hessenberg_Least_Squares_d(z, Q, b, tol_pivot, largeChange);
    if (i < 0)
    { // error in solving
      break;
    }

    // update x = x - z
    for (i = 0; i < n; i++)
    {
      sub_d(&x->coord[i], &x->coord[i], &z->coord[i]);
    }        
    
    // copy currSV to prevSV
    prevSV = currSV;

    // find currSV = norm of (R * x) / norm of x
    currSV = 0;
    norm = 0;
    for (i = 0; i < n; i++)
    {
      set_zero_d(&z->coord[i]);
      for (j = 0; j < n; j++)
      {
        sum_mul_d(&z->coord[i], &R->entry[i][j], &x->coord[j]);
      }
      currSV += z->coord[i].r * z->coord[i].r + z->coord[i].i * z->coord[i].i;
      norm += x->coord[i].r * x->coord[i].r + x->coord[i].i * x->coord[i].i;
    }

    // make sure we can divide
    if (norm > 0)
    {
      currSV /= norm;
      currSV = sqrt(currSV);

      // check for convergence
      if (fabs(currSV - prevSV) < error_tol)
      {
        retVal = 0;
        break;
      }
    } 
    else
    { // we have to stop since norm(x) = 0
      currSV = prevSV;

      retVal = 0;
      break;
    }

    // increment the number of iterations
    its++;
  }

  // set min_sv
  if (t > 1)
  {
    *min_sv = currSV * t;
  }
  else
  {
    *min_sv = currSV;
  }

  clear_vec_d(x); clear_vec_d(z); clear_vec_d(b);
  clear_mat_d(P); clear_mat_d(Q); clear_mat_d(R);

  // return
  return retVal;
}

//// GIVENS ROTATIONS ///////

void gengr_d(comp_d f, comp_d g, comp_d top, mat_d A, int r, int c, double tol_sign)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: [f,g] is the Givens rotation for [A(r-1,c),A(r,c)], r>0*
\***************************************************************/
{
  double norm_sqr[2], norm[2], tempD[2];
  comp_d tempComp;  

  // error checking
  if (r <= 0)
  {
    printf("ERROR: The element to clear cannot be in the first row to generate a Givens rotation!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (r >= A->rows || c >= A->cols || c < 0)
  {
    printf("ERROR: Cannot generate a Givens rotations because the location does not exist in the matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have 1 <= r < A->rows && 0 <= c < A->cols
  
  // find the square of the norms and the norms
  norm_sqr[0] = A->entry[r-1][c].r * A->entry[r-1][c].r + A->entry[r-1][c].i * A->entry[r-1][c].i;
  norm_sqr[1] = A->entry[r][c].r * A->entry[r][c].r + A->entry[r][c].i * A->entry[r][c].i;

  norm[0] = sqrt(norm_sqr[0]);
  norm[1] = sqrt(norm_sqr[1]);

  if (norm[0] < tol_sign && norm[1] < tol_sign)
  { // both are really small
    if (norm[1] < norm[0])
    { // f = 1, g = 0, top = (r-1,c)
      set_double_d(f, 1, 0);
      set_zero_d(g);
      set_d(top, &A->entry[r-1][c]);
    }
    else
    { // f = 0, g = 1, top = (r, c)
      set_zero_d(f);
      set_double_d(g, 1, 0);
      set_d(top, &A->entry[r][c]);
    }
  }
  else if (norm[1] < norm[0])
  { // (r-1,c) is larger than (r,c)
    if (norm[1] < tol_sign)
    { // (r,c) is very small - f = conj(sign((r-1,c))), g = 0, top = norm[0]
      sign_d(f, &A->entry[r-1][c], tol_sign);
      conjugate_d(f, f);
      set_zero_d(g);
      set_double_d(top, norm[0], 0);
    }
    else
    { // f = norm[0] / sqrt(norm_sqr[0] + norm_sqr[1]), g = sign((r-1,c)) * conj((r,c))/ sqrt(norm_sqr[0] + norm_sqr[1]), top = sign((r-1,c)) * sqrt(norm_sqr[0] + norm_sqr[1])
      tempD[0] = norm_sqr[0] + norm_sqr[1];
      tempD[1] = 1 / sqrt(norm_sqr[0] * tempD[0]);
      set_double_d(f, norm_sqr[0] * tempD[1], 0);

      tempD[0] *= tempD[1];
      mul_rdouble_d(tempComp, &A->entry[r-1][c], tempD[0]);

      mul_rdouble_d(g, &A->entry[r-1][c], tempD[1]);

      // setup top after (r-1,c) value is not needed anymore
      set_d(top, tempComp);

      conjugate_d(tempComp, &A->entry[r][c]); 
      mul_d(g, g, tempComp);
    }
  }
  else
  { // (r,c) is larger than (r-1,c)
    if (norm[0] < tol_sign)
    { // (r-1,c) is very small - f = 0, g = conj(sign((r,c))), top = norm[1]
      set_zero_d(f);
      sign_d(g, &A->entry[r][c], tol_sign);
      conjugate_d(g, g);
      set_double_d(top, norm[1], 0);
    }
    else
    { // f = sign((r,c)) * conj((r-1,c)) / sqrt(norm_sqr[0] + norm_sqr[1]), g = norm[1] / sqrt(norm_sqr[0] + norm_sqr[1]), top = sign((r,c)) * sqrt(norm_sqr[0] + norm_sqr[1])
      tempD[0] = norm_sqr[0] + norm_sqr[1];
      tempD[1] = 1 / sqrt(norm_sqr[1] * tempD[0]);
      set_double_d(g, norm_sqr[1] * tempD[1], 0);

      conjugate_d(tempComp, &A->entry[r-1][c]);

      // setup top after (r-1,c) value is not needed anymore
      tempD[0] *= tempD[1];
      mul_rdouble_d(top, &A->entry[r][c], tempD[0]);

      mul_rdouble_d(f, &A->entry[r][c], tempD[1]);
      mul_d(f, f, tempComp);
    }
  }

  return;
}

/////// MULTI PRECISION ///////

void approx_min_svd_mp(mpf_t min_sv, mat_mp A)
/***************************************************************\
* USAGE: approximates the minimum singular value of A using QLP *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int prec = mpf_get_prec(A->entry[0][0].r);

  approx_min_svd_mp_prec(min_sv, A, prec);

  return;
}

void approx_min_svd_mp_prec(mpf_t min_sv, mat_mp A, int prec)
/***************************************************************\
* USAGE: approximates the minimum singular value of A using QLP *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int n;
  mat_mp L;

  init_mat_mp2(L, 0, 0, prec);

  // compute a QLP decomposition
  QLP_L_mp_prec(L, A, prec);

  // min_sv is approximated by last diagonal entry
  n = MIN(L->rows, L->cols);
  if (n > 0)
  {
    mpf_abs_mp(min_sv, &L->entry[n-1][n-1]);
  }
  else
  {
    mpf_set_ui(min_sv, 0);
  }

  clear_mat_mp(L);

  return;
}

void min_svd_mp_prec(mpf_t min_sv, mat_mp Jv, int prec)
/***************************************************************\
* USAGE: finds the minimum singular value of A using Li/Zeng    *
* where the tolerance is set to the epsilon of the precision    *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if convergence to error_tol, -1 otherwise    *
* NOTES:                                                        *
\***************************************************************/
{ 
  int curr_digits = prec_to_digits(prec) - 2;
  size_t size;
  char *str = NULL;
  mpf_t error_tol;

  // setup error_tol
  mpf_init2(error_tol, prec);
  size = 1 + snprintf(NULL, 0, "1e-%d", curr_digits);
  str = (char *)bmalloc(size * sizeof(char));
  sprintf(str, "1e-%d", curr_digits);
  mpf_set_str(error_tol, str, 10);

  // find the minimum singular value
  min_svd_mp(min_sv, Jv, error_tol);

  // clear error_tol & str
  mpf_clear(error_tol);
  free(str);

  return;
}

int min_svd_mp(mpf_t min_sv, mat_mp A, mpf_t error_tol)
/***************************************************************\
* USAGE: finds the minimum singular value of A using Li/Zeng    *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if convergence to error_tol, -1 otherwise    *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, its = 0, maxIts = 50, retVal = -1, m = A->rows, n = A->cols, num_digits = - (int) floor(mpf_get_prec(A->entry[0][0].r) * log10(2.0) - 2.5);
  char *str = NULL;
  mpf_t norm, t, tempMPF, currSV, prevSV, tol_pivot, tol_sign, largeChange;
  vec_mp x, z, b;
  mat_mp P, Q, R;

  // error checking
  if (mpf_cmp_ui(error_tol, 0) < 0)
  {
    printf("ERROR: The error tolerance for min_svd_mp needs to be non-negative!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // make sure m > 0 && n > 0
  if (m == 0 || n == 0)
  { // a dimension is 0
    mpf_set_ui(min_sv, 0);
    retVal = 0;
    return retVal;
  }

  // initialize MP
  mpf_init(norm); mpf_init(t); mpf_init(tempMPF); mpf_init(currSV); 
  mpf_init(prevSV); mpf_init(tol_pivot); mpf_init(tol_sign); mpf_init(largeChange);
  init_vec_mp(x, 0); init_vec_mp(z, 0); init_vec_mp(b, 0);
  init_mat_mp(P, 0, 0); init_mat_mp(Q, 0, 0); init_mat_mp(R, 0, 0);

  // initialize values
  mpf_set_d(currSV, 1e300); mpf_set_d(prevSV, 1e300);
  // set tol_pivot
  i = 1 + snprintf(NULL, 0, "1e%d", num_digits);
  str = (char *)brealloc(str, i * sizeof(char));
  sprintf(str, "1e%d", num_digits);
  mpf_set_str(tol_pivot, str, 10);
  // set tol_sign
  i = 1 + snprintf(NULL, 0, "1e%d", 2 * num_digits + 2);
  str = (char *)brealloc(str, i * sizeof(char));
  sprintf(str, "1e%d", 2 * num_digits + 2);
  mpf_set_str(tol_sign, str, 10);
  // set largeChange
  i = 1 + snprintf(NULL, 0, "1e%d", -num_digits - 2);
  str = (char *)brealloc(str, i * sizeof(char));
  sprintf(str, "1e%d", -num_digits - 2);
  mpf_set_str(largeChange, str, 10);

  // need m >= n for the R in A = QR
  if (m < n)
  { // work with the transpose
    mat_mp tempMat;
    init_mat_mp(tempMat, 0, 0);

    transpose_mp(tempMat, A);
    m = tempMat->rows;
    n = tempMat->cols;

    // find QR-decomposition of A' - tempMat * P = Q * R
    QR_mp(Q, R, P, tempMat, tol_pivot, tol_sign, largeChange, 1);

    clear_mat_mp(tempMat);
  }
  else
  { // find QR-decomposition of A - A * P = Q * R
    QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, 1);
  }
  // so we have R being an m x n upper triangular matrix with m >= n >= 1
  // since m >= n & upper triangular, the last (m-n) rows are 0 - so we can treat R as n x n upper triangular
  R->rows = R->cols = n;

  // find the norm of the R[0][0]
  mpf_abs_mp(t, &R->entry[0][0]);

  // normalize R - an n x n upper triangular matrix
  if (mpf_cmp_ui(t, 1) > 0)
  {
    mpf_ui_div(norm, 1, t);
    for (i = 0; i < n; i++)
      for (j = i; j < n; j++)
      {
        mul_rmpf_mp(&R->entry[i][j], &R->entry[i][j], norm);
      }
  }

  // generate a random unit vector x
  increase_size_vec_mp(x, n);
  x->size = n;
  mpf_set_ui(norm, 0);
  for (i = 0; i < n; i++)
  {
    get_comp_rand_mp(&x->coord[i]);
    mpf_mul(tempMPF, x->coord[i].r, x->coord[i].r);
    mpf_add(norm, norm, tempMPF);
    mpf_mul(tempMPF, x->coord[i].i, x->coord[i].i);
    mpf_add(norm, norm, tempMPF);
  }
  mpf_sqrt(norm, norm);
  mpf_ui_div(norm, 1, norm);
  for (i = 0; i < n; i++)
  {
    mul_rmpf_mp(&x->coord[i], &x->coord[i], norm);
  }

  // setup the size of Q & b
  increase_size_mat_mp(Q, n + 1, n);
  increase_size_vec_mp(b, n + 1);
  Q->rows = b->size = n + 1;
  Q->cols = n;

  // setup the bottom of Q - always == R
  for (i = 1; i <= n; i++)
    for (j = 0; j < n; j++)
    {
      set_mp(&Q->entry[i][j], &R->entry[i-1][j]);
    }

  // main loop - tau is not used since we normalized R
  while (its < maxIts)
  { // setup Q & b
    for (i = 0; i <= n; i++)
      if (i == 0)
      { // setup top row of Q & find x' * x
        mpf_set_ui(norm, 0);
        for (j = 0; j < n; j++)
        { // notice the conjugation
          mpf_mul_ui(Q->entry[i][j].r, x->coord[j].r, 2);
          mpf_mul_ui(Q->entry[i][j].i, x->coord[j].i, 2);
          mpf_neg(Q->entry[i][j].i, Q->entry[i][j].i);
          mpf_mul(tempMPF, x->coord[j].r, x->coord[j].r);
          mpf_add(norm, norm, tempMPF);
          mpf_mul(tempMPF, x->coord[j].i, x->coord[j].i);
          mpf_add(norm, norm, tempMPF);
        }

        // setup top entry in b
        mpf_sub_ui(b->coord[0].r, norm, 1);
        mpf_set_ui(b->coord[0].i, 0);
      }
      else
      { // find ith entry of b = (i-1)st entry of R * x
        set_zero_mp(&b->coord[i]);
        for (j = 0; j < n; j++)
        {
          sum_mul_mp(&b->coord[i], &R->entry[i-1][j], &x->coord[j]);
        }
      }

    // find the least squares solution z to Q*z = b
    i = matrixSolve_Hessenberg_Least_Squares_mp(z, Q, b, tol_pivot, largeChange);
    if (i < 0)
    { // error in solving
      break;
    }

    // update x = x - z
    for (i = 0; i < n; i++)
    {
      sub_mp(&x->coord[i], &x->coord[i], &z->coord[i]);
    }

    // copy currSV to prevSV
    mpf_set(prevSV, currSV);

    // find currSV = norm of R * x / norm of x
    mpf_set_ui(currSV, 0);
    mpf_set_ui(norm, 0);
    for (i = 0; i < n; i++)
    {
      set_zero_mp(&z->coord[i]);
      for (j = 0; j < n; j++)
      {
        sum_mul_mp(&z->coord[i], &R->entry[i][j], &x->coord[j]);
      }
      mpf_mul(tempMPF, z->coord[i].r, z->coord[i].r);
      mpf_add(currSV, currSV, tempMPF);
      mpf_mul(tempMPF, z->coord[i].i, z->coord[i].i);
      mpf_add(currSV, currSV, tempMPF);

      mpf_mul(tempMPF, x->coord[i].r, x->coord[i].r);
      mpf_add(norm, norm, tempMPF);
      mpf_mul(tempMPF, x->coord[i].i, x->coord[i].i);
      mpf_add(norm, norm, tempMPF);
    }

    // make sure we can divide
    if (mpf_cmp_ui(norm, 0) > 0)
    {
      mpf_div(currSV, currSV, norm);
      mpf_sqrt(currSV, currSV);

      // check for convergence
      mpf_sub(tempMPF, currSV, prevSV);
      mpf_abs(tempMPF, tempMPF);
      if (mpf_cmp(tempMPF, error_tol) < 0)
      {
        retVal = 0;
        break;
      }
    }
    else
    { // we have to stop since norm(x) = 0
      mpf_set(currSV, prevSV);

      retVal = 0;
      break;
    }

    // increment the number of iterations
    its++;
  }

  // set min_sv
  if (mpf_cmp_ui(t, 1) > 0)
  {
    mpf_mul(min_sv, currSV, t);
  }
  else
  {
    mpf_set(min_sv, currSV);
  }

  // clear memory & MP
  free(str);
  mpf_clear(norm); mpf_clear(t); mpf_clear(tempMPF); mpf_clear(currSV);
  mpf_clear(prevSV); mpf_clear(tol_pivot); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_mp(x); clear_vec_mp(z); clear_vec_mp(b);
  clear_mat_mp(P); clear_mat_mp(Q); clear_mat_mp(R);

  // return
  return retVal;
}

//// GIVENS ROTATIONS ///////

void gengr_mp(comp_mp f, comp_mp g, comp_mp top, mat_mp A, int r, int c, mpf_t tol_sign)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: [f,g] is the Givens rotation for [A(r-1,c),A(r,c)], r>0*
\***************************************************************/
{
  mpf_t norm_sqr[2], norm[2], tempMPF[2];
  comp_mp tempComp;

  // error checking
  if (r <= 0)
  {
    printf("ERROR: The element to clear cannot be in the first row to generate a Givens rotation!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (r >= A->rows || c >= A->cols || c < 0)
  {
    printf("ERROR: Cannot generate a Givens rotations because the location does not exist in the matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have 1 <= r < A->rows && 0 <= c < A->cols

  // initialize MP
  mpf_init(norm_sqr[0]); mpf_init(norm_sqr[1]); 
  mpf_init(norm[0]); mpf_init(norm[1]);
  mpf_init(tempMPF[0]); mpf_init(tempMPF[1]);
  init_mp(tempComp);

  // find the square of the norms and the norms
  mpf_mul(tempMPF[0], A->entry[r-1][c].r, A->entry[r-1][c].r);
  mpf_mul(tempMPF[1], A->entry[r-1][c].i, A->entry[r-1][c].i);
  mpf_add(norm_sqr[0], tempMPF[0], tempMPF[1]);

  mpf_mul(tempMPF[0], A->entry[r][c].r, A->entry[r][c].r);
  mpf_mul(tempMPF[1], A->entry[r][c].i, A->entry[r][c].i);
  mpf_add(norm_sqr[1], tempMPF[0], tempMPF[1]);

  mpf_sqrt(norm[0], norm_sqr[0]);
  mpf_sqrt(norm[1], norm_sqr[1]);

  if (mpf_cmp(norm[0], tol_sign) < 0 && mpf_cmp(norm[1], tol_sign) < 0)
  { // both are really small
    if (mpf_cmp(norm[1], norm[0]) < 0)
    { // f = 1, g = 0, top = (r-1,c)
      mpf_set_ui(f->r, 1); 
      mpf_set_ui(f->i, 0);
      set_zero_mp(g);
      set_mp(top, &A->entry[r-1][c]);
    }
    else
    { // f = 0, g = 1, top = (r, c)
      set_zero_mp(f);
      mpf_set_ui(g->r, 1);
      mpf_set_ui(g->i, 0);
      set_mp(top, &A->entry[r][c]);
    }
  }
  else if (mpf_cmp(norm[1], norm[0]) < 0)
  { // (r-1,c) is larger than (r,c)
    if (mpf_cmp(norm[1], tol_sign) < 0)
    { // (r,c) is very small - f = conj(sign((r-1,c))), g = 0, top = norm[0]
      sign_mp2(f, &A->entry[r-1][c], tol_sign);
      conjugate_mp(f, f);
      set_zero_mp(g);
      mpf_set(top->r, norm[0]); 
      mpf_set_ui(top->i, 0);
    }
    else
    { // f = norm[0] / sqrt(norm_sqr[0] + norm_sqr[1]), g = sign((r-1,c)) * conj((r,c))/ sqrt(norm_sqr[0] + norm_sqr[1]), top = sign((r-1,c)) * sqrt(norm_sqr[0] + norm_sqr[1])
      mpf_add(tempMPF[0], norm_sqr[0], norm_sqr[1]);
      mpf_mul(tempMPF[1], norm_sqr[0], tempMPF[0]);
      mpf_sqrt(tempMPF[1], tempMPF[1]);
      mpf_ui_div(tempMPF[1], 1, tempMPF[1]);

      mpf_mul(f->r, norm_sqr[0], tempMPF[1]);
      mpf_set_ui(f->i, 0);

      mul_rmpf_mp(g, &A->entry[r-1][c], tempMPF[1]);

      // setup top after (r-1,c) value is not needed anymore
      mpf_mul(tempMPF[0], tempMPF[0], tempMPF[1]);
      mul_rmpf_mp(top, &A->entry[r-1][c], tempMPF[0]);

      conjugate_mp(tempComp, &A->entry[r][c]);
      mul_mp(g, g, tempComp);
    }
  }
  else
  { // (r,c) is larger than (r-1,c)
    if (mpf_cmp(norm[0], tol_sign) < 0)
    { // (r-1,c) is very small - f = 0, g = conj(sign((r,c))), top = norm[1]
      set_zero_mp(f);
      sign_mp2(g, &A->entry[r][c], tol_sign);
      conjugate_mp(g, g);
      mpf_set(top->r, norm[1]);
      mpf_set_ui(top->i, 0);
    }
    else
    { // f = sign((r,c)) * conj((r-1,c)) / sqrt(norm_sqr[0] + norm_sqr[1]), g = norm[1] / sqrt(norm_sqr[0] + norm_sqr[1]), top = sign((r,c)) * sqrt(norm_sqr[0] + norm_sqr[1])
      mpf_add(tempMPF[0], norm_sqr[0], norm_sqr[1]);
      mpf_mul(tempMPF[1], norm_sqr[1], tempMPF[0]);
      mpf_sqrt(tempMPF[1], tempMPF[1]);
      mpf_ui_div(tempMPF[1], 1, tempMPF[1]);

      mpf_mul(g->r, norm_sqr[1], tempMPF[1]);
      mpf_set_ui(g->i, 0);

      conjugate_mp(tempComp, &A->entry[r-1][c]);

      // setup top after (r-1,c) value is not needed anymore
      mpf_mul(tempMPF[0], tempMPF[0], tempMPF[1]);
      mul_rmpf_mp(top, &A->entry[r][c], tempMPF[0]);

      mul_rmpf_mp(f, &A->entry[r][c], tempMPF[1]);
      mul_mp(f, f, tempComp);
    }
  }

  mpf_clear(norm_sqr[0]); mpf_clear(norm_sqr[1]);
  mpf_clear(norm[0]); mpf_clear(norm[1]);
  mpf_clear(tempMPF[0]); mpf_clear(tempMPF[1]);
  clear_mp(tempComp);

  return;
}
