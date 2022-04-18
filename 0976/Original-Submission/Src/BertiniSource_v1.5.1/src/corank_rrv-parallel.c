// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// uses two different approximations to approximate the corank of the (block triangular) matrix

/////////// double precision ////////////

void setup_frobenius_d(double **FV, mat_d m0, mat_d m1)
/***************************************************************\
* USAGE: setup the Frobenius norms for the two matrices         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume m0 & m1 are n x n matrices                      *
\***************************************************************/
{
  int i, j, n = m0->rows;

  // error checking
  if (m0->cols != n || m1->rows != n || m1->cols != n)
  {
    printf("ERROR: To find the Frobenius norms, the matrices need to be square!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (n > 0)
  { // setup the last one
    FV[0][n-1] = norm_sqr_d(&m0->entry[n-1][n-1]);
    FV[1][n-1] = norm_sqr_d(&m1->entry[n-1][n-1]);

    // setup the rest of them
    for (i = n - 2; i >= 0; i--)
    { // initialize
      FV[0][i] = FV[0][i+1] + norm_sqr_d(&m0->entry[i][i]);
      FV[1][i] = FV[1][i+1] + norm_sqr_d(&m1->entry[i][i]);

      for (j = i + 1; j < n; j++)
      { // add on
        FV[0][i] += norm_sqr_d(&m0->entry[j][i]) + norm_sqr_d(&m0->entry[i][j]);
        FV[1][i] += norm_sqr_d(&m1->entry[j][i]) + norm_sqr_d(&m1->entry[i][j]);
      }

      // take the previous sqrt
      FV[0][i+1] = sqrt(FV[0][i+1]);
      FV[1][i+1] = sqrt(FV[1][i+1]);
    }
    // take the final sqrt
    FV[0][0] = sqrt(FV[0][0]);
    FV[1][0] = sqrt(FV[1][0]);
  }

  return;
}

int corank_rrv_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0, mat_d mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio, double FV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of null space                             *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero FV*
\***************************************************************/
{ // assume mat0 is the more accurate one
  int i, j, num_FV, bad_loop, corank = 0;
  double ratio, normalization_factor, tol_conv = 1e-15, tol_sign = 1e-20, largeChange = 1e13, **FV = (double **)bmalloc(2 * sizeof(double *));
  mat_d L0, L1;

  // make sure that mat0 & mat1 are the same size
  if (mat0->rows != mat1->rows || mat0->cols != mat1->cols)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_mat_d(L0, 0, 0);
  init_mat_d(L1, 0, 0);

  // compute the size of the null space
  if (mat0->rows >= mat0->cols)
  { // size of null space == corank

    // do the rank revealing decomposition normally
    QLP_block_pair_d(L0, L1, mat0, mat1, rowTop, colTop, tol_conv, tol_sign, largeChange);
  }
  else // make the matrices tall (more rows than columns)
  { // size of null space == corank + (cols - rows) == corank of [[mat][0]]
    increase_size_mat_d(L0, mat0->cols, mat0->cols);
    increase_size_mat_d(L1, mat1->cols, mat1->cols);
    for (j = 0; j < mat0->cols; j++)
    { // copy the top rows
      for (i = 0; i < mat0->rows; i++)
      {
        set_d(&L0->entry[i][j], &mat0->entry[i][j]);
        set_d(&L1->entry[i][j], &mat1->entry[i][j]);
      }
      // set bottom row to 0
      for (i = mat0->rows; i < mat0->cols; i++)
      {
        set_zero_d(&L0->entry[i][j]);
        set_zero_d(&L1->entry[i][j]);
      }
    }
    L0->rows = L1->rows = L0->cols = L1->cols = mat0->cols;

    // do the rank revealing decomposition
    QLP_block_pair_d(L0, L1, L0, L1, rowTop, colTop, tol_conv, tol_sign, largeChange);
  }

  // find the number of values
  num_FV = L0->rows;

  // allocate for the Frobenius norms of the matrices
  FV[0] = (double *)bmalloc(num_FV * sizeof(double));
  FV[1] = (double *)bmalloc(num_FV * sizeof(double));

  // setup FV
  setup_frobenius_d(FV, L0, L1);

  // estimate the condition number
  if (FV[0][num_FV - 1] == 0 || FV[1][num_FV - 1] == 0)
  { // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    *CN = FV[0][0] / FV[0][num_FV - 1];
  }

  // look to normalize the entries so that SV[0][0] & SV[1][0] <= 1
  if (FV[0][0] > 1 || FV[1][0] > 1)
  { // find normalization_factor = max_entry
    normalization_factor = MAX(FV[0][0], FV[1][0]);

    // normalize SV[0] & SV[1] by multiplying by 1 / normalization_factor
    ratio = 1 / normalization_factor;
    for (i = 0; i < num_FV; i++)
    {
      FV[0][i] *= ratio;
      FV[1][i] *= ratio;
    }
  }
  else
  { // normalization_factor = 1
    normalization_factor = 1;
  }

  // initialize corank & bad_loop
  corank = num_FV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_FV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (FV[0][i] == 0 || FV[1][i] == 0)
      bad_loop = 1;
    else
    { // find the ratio of the norms (min / max)
      if (FV[0][i] < FV[1][i])
        ratio = FV[0][i] / FV[1][i];
      else
        ratio = FV[1][i] / FV[0][i];

      // determine if this is okay or not
      if (ratio <= max_FV_ratio || FV[0][i] < FV_tol || FV[1][i] < FV_tol)
        bad_loop = 1;
    }      

    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    *smallest_nonzero = MAX(FV[0][num_FV - 1], FV[1][num_FV - 1]) * normalization_factor;
    *largest_zero = 0;
  }
  else if (corank == num_FV)
  { // all are zero
    *smallest_nonzero = 0;
    *largest_zero = MAX(FV[0][0], FV[1][0]) * normalization_factor;
  }
  else // 0 < corank < num_FV
  { // num_FV - corank - 1 is smallest non-zero & num_FV - corank is largest zero
    *smallest_nonzero = MAX(FV[0][num_FV - corank - 1], FV[1][num_FV - corank - 1]) * normalization_factor;
    *largest_zero = MAX(FV[0][num_FV - corank], FV[1][num_FV - corank]) * normalization_factor;
  }

  // release memory
  free(FV[0]); free(FV[1]); free(FV);
  clear_mat_d(L0); clear_mat_d(L1); 

  return corank;
}

/////////// multi precision ////////////

void setup_frobenius_mp(mpf_t **FV, mat_mp m0, mat_mp m1)
/***************************************************************\
* USAGE: setup the Frobenius norms for the two matrices         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume m0 & m1 are n x n matrices                      *
\***************************************************************/
{
  int i, j, n = m0->rows;
  mpf_t tempMPF;

  mpf_init(tempMPF);

  // error checking
  if (m0->cols != n || m1->rows != n || m1->cols != n)
  {
    printf("ERROR: To find the Frobenius norms, the matrices need to be square!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (n > 0)
  { // setup the last one
    mpf_init(FV[0][n-1]); mpf_init(FV[1][n-1]);

    norm_sqr_mp(FV[0][n-1], &m0->entry[n-1][n-1]);
    norm_sqr_mp(FV[1][n-1], &m1->entry[n-1][n-1]);

    // setup the rest of them
    for (i = n - 2; i >= 0; i--)
    { // initialize
      mpf_init(FV[0][i]); mpf_init(FV[1][i]);

      mpf_set(FV[0][i], FV[0][i+1]);
      mpf_set(FV[1][i], FV[1][i+1]);

      norm_sqr_mp(tempMPF, &m0->entry[i][i]);
      mpf_add(FV[0][i], FV[0][i], tempMPF);

      norm_sqr_mp(tempMPF, &m1->entry[i][i]);
      mpf_add(FV[1][i], FV[1][i], tempMPF);

      for (j = i + 1; j < n; j++)
      { // add on
        norm_sqr_mp(tempMPF, &m0->entry[j][i]);
        mpf_add(FV[0][i], FV[0][i], tempMPF);
        norm_sqr_mp(tempMPF, &m0->entry[i][j]);
        mpf_add(FV[0][i], FV[0][i], tempMPF);

        norm_sqr_mp(tempMPF, &m1->entry[j][i]);
        mpf_add(FV[1][i], FV[1][i], tempMPF);
        norm_sqr_mp(tempMPF, &m1->entry[i][j]);
        mpf_add(FV[1][i], FV[1][i], tempMPF);
      }

      // take the previous sqrt
      mpf_sqrt(FV[0][i+1], FV[0][i+1]);
      mpf_sqrt(FV[1][i+1], FV[1][i+1]);
    }
    // take the final sqrt
    mpf_sqrt(FV[0][0], FV[0][0]);
    mpf_sqrt(FV[1][0], FV[1][0]);
  }

  mpf_clear(tempMPF);

  return;
}

int corank_rrv_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, mat_mp mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio, double FV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of null space                             *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero FV*
\***************************************************************/
{ // assume mat0 is the more accurate one
  int i, j, num_FV, bad_loop, corank = 0, curr_prec = mpf_get_default_prec();
  mpf_t ratio, normalization_factor, **FV = (mpf_t **)bmalloc(2 * sizeof(mpf_t *));
  mat_mp L0, L1;

  // make sure that mat0 & mat1 are the same size
  if (mat0->rows != mat1->rows || mat0->cols != mat1->cols)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  mpf_init(ratio); mpf_init(normalization_factor);
  init_mat_mp(L0, 0, 0);
  init_mat_mp(L1, 0, 0);

  // compute the size of the null space
  if (mat0->rows >= mat0->cols)
  { // size of null space == corank

    // do the rank revealing decomposition normally
    QLP_block_pair_mp_prec(L0, L1, mat0, mat1, rowTop, colTop, curr_prec);
  }
  else // make the matrices tall (more rows than columns)
  { // size of null space == corank + (cols - rows) == corank of [[mat][0]]
    increase_size_mat_mp(L0, mat0->cols, mat0->cols);
    increase_size_mat_mp(L1, mat1->cols, mat1->cols);
    for (j = 0; j < mat0->cols; j++)
    { // copy the top rows
      for (i = 0; i < mat0->rows; i++)
      {
        set_mp(&L0->entry[i][j], &mat0->entry[i][j]);
        set_mp(&L1->entry[i][j], &mat1->entry[i][j]);
      }
      // set bottom row to 0
      for (i = mat0->rows; i < mat0->cols; i++)
      {
        set_zero_mp(&L0->entry[i][j]);
        set_zero_mp(&L1->entry[i][j]);
      }
    }
    L0->rows = L1->rows = L0->cols = L1->cols = mat0->cols;

    // do the rank revealing decomposition
    QLP_block_pair_mp_prec(L0, L1, L0, L1, rowTop, colTop, curr_prec);
  }

  // find the number of values
  num_FV = L0->rows;

  // allocate for the Frobenius norms of the matrices
  FV[0] = (mpf_t *)bmalloc(num_FV * sizeof(mpf_t));
  FV[1] = (mpf_t *)bmalloc(num_FV * sizeof(mpf_t));

  // setup FV
  setup_frobenius_mp(FV, L0, L1);

  // estimate the condition number
  if (mpfr_zero_p(FV[0][num_FV - 1]))
  { // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, FV[0][0], FV[0][num_FV - 1]);
    *CN = mpf_get_d(ratio);
  }

  // look to normalize the entries so that FV[0][0] & FV[1][0] <= 1
  if (mpf_cmp_ui(FV[0][0], 1) > 0 || mpf_cmp_ui(FV[1][0], 1) > 0)
  { // find normalization_factor = max_entry
    mpfr_max(normalization_factor, FV[0][0], FV[1][0], __gmp_default_rounding_mode);

    // normalize FV[0] & FV[1] by multiplying by 1 / normalization_factor
    mpf_ui_div(ratio, 1, normalization_factor);
    for (i = 0; i < num_FV; i++)
    {
      mpf_mul(FV[0][i], FV[0][i], ratio);
      mpf_mul(FV[1][i], FV[1][i], ratio);
    }
  }
  else
  { // normalization_factor = 1
    mpf_set_ui(normalization_factor, 1);
  }

  // initialize corank & bad_loop
  corank = num_FV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_FV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (mpfr_zero_p(FV[0][i]) || mpfr_zero_p(FV[1][i]))
      bad_loop = 1;
    else
    { // find the ratio of the singular values (min / max)
      if (mpf_cmp(FV[0][i], FV[1][i]) < 0)
        mpf_div(ratio, FV[0][i], FV[1][i]);
      else 
        mpf_div(ratio, FV[1][i], FV[0][i]);

      // determine if this is okay or not
      if (mpf_cmp_d(ratio, max_FV_ratio) <= 0 || mpf_cmp_d(FV[0][i], FV_tol) < 0 || mpf_cmp_d(FV[1][i], FV_tol) < 0)
        bad_loop = 1;
    }

    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    mpfr_max(ratio, FV[0][num_FV - 1], FV[1][num_FV - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    *largest_zero = 0;
  }
  else if (corank == num_FV)
  { // all are zero
    *smallest_nonzero = 0;

    mpfr_max(ratio, FV[0][0], FV[1][0], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }
  else // 0 < corank < num_FV
  { // num_FV - corank - 1 is smallest non-zero & num_FV - corank is largest zero
    mpfr_max(ratio, FV[0][num_FV - corank - 1], FV[1][num_FV - corank - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    mpfr_max(ratio, FV[0][num_FV - corank], FV[1][num_FV - corank], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }

  // release memory
  for (i = num_FV - 1; i >= 0 ; i--)
  {
    mpf_clear(FV[0][i]);
    mpf_clear(FV[1][i]);
  }
  free(FV[0]); free(FV[1]); free(FV);
  mpf_clear(ratio); mpf_clear(normalization_factor);
  clear_mat_mp(L0); clear_mat_mp(L1); 

  return corank;
}

/////////// adaptive multi precision ////////////

void setup_frobenius_amp_mp_d(mpf_t **FV, mat_mp m0, mat_d m1)
/***************************************************************\
* USAGE: setup the Frobenius norms for the two matrices         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume m0 & m1 are n x n matrices                      *
\***************************************************************/
{
  int i, j, n = m0->rows;
  double tempD;
  mpf_t tempMPF;

  mpf_init(tempMPF);

  // error checking
  if (m0->cols != n || m1->rows != n || m1->cols != n)
  {
    printf("ERROR: To find the Frobenius norms, the matrices need to be square!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (n > 0)
  { // setup the last one
    mpf_init(FV[0][n-1]); mpf_init(FV[1][n-1]);

    norm_sqr_mp(FV[0][n-1], &m0->entry[n-1][n-1]);
    tempD = norm_sqr_d(&m1->entry[n-1][n-1]);
    mpf_set_d(FV[1][n-1], tempD);

    // setup the rest of them
    for (i = n - 2; i >= 0; i--)
    { // initialize
      mpf_init(FV[0][i]); mpf_init(FV[1][i]);

      mpf_set(FV[0][i], FV[0][i+1]);
      mpf_set(FV[1][i], FV[1][i+1]);

      norm_sqr_mp(tempMPF, &m0->entry[i][i]);
      mpf_add(FV[0][i], FV[0][i], tempMPF);

      tempD = norm_sqr_d(&m1->entry[i][i]);
      mpf_set_d(tempMPF, tempD);
      mpf_add(FV[1][i], FV[1][i], tempMPF);

      for (j = i + 1; j < n; j++)
      { // add on
        norm_sqr_mp(tempMPF, &m0->entry[j][i]);
        mpf_add(FV[0][i], FV[0][i], tempMPF);
        norm_sqr_mp(tempMPF, &m0->entry[i][j]);
        mpf_add(FV[0][i], FV[0][i], tempMPF);

        tempD = norm_sqr_d(&m1->entry[j][i]);
        mpf_set_d(tempMPF, tempD);
        mpf_add(FV[1][i], FV[1][i], tempMPF);
        tempD = norm_sqr_d(&m1->entry[i][j]);
        mpf_set_d(tempMPF, tempD);
        mpf_add(FV[1][i], FV[1][i], tempMPF);
      }

      // take the previous sqrt
      mpf_sqrt(FV[0][i+1], FV[0][i+1]);
      mpf_sqrt(FV[1][i+1], FV[1][i+1]);
    }
    // take the final sqrt
    mpf_sqrt(FV[0][0], FV[0][0]);
    mpf_sqrt(FV[1][0], FV[1][0]);
  }

  mpf_clear(tempMPF);

  return;
}

int corank_rrv_amp_mp_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, int prec0, mat_d mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of null space                             *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero FV*
\***************************************************************/
{ // assume mat0 is the more accurate one
  int i, j, num_FV, bad_loop, corank = 0, curr_prec = prec0;
  mpf_t ratio, normalization_factor, **FV = (mpf_t **)bmalloc(2 * sizeof(mpf_t *));
  mat_mp L0;
  mat_d L1;

  // make sure that mat0 & mat1 are the same size
  if (mat0->rows != mat1->rows || mat0->cols != mat1->cols)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  mpf_init(ratio); mpf_init(normalization_factor);
  init_mat_mp2(L0, 0, 0, curr_prec);
  init_mat_d(L1, 0, 0);

  // compute the size of the null space
  if (mat0->rows >= mat0->cols)
  { // size of null space == corank

    // do the rank revealing decomposition normally
    QLP_block_pair_amp_mp_d_prec(L0, L1, mat0, prec0, mat1, rowTop, colTop);
  }
  else // make the matrices tall (more rows than columns)
  { // size of null space == corank + (cols - rows) == corank of [[mat][0]]
    increase_size_mat_mp(L0, mat0->cols, mat0->cols);
    increase_size_mat_d(L1, mat1->cols, mat1->cols);
    for (j = 0; j < mat0->cols; j++)
    { // copy the top rows
      for (i = 0; i < mat0->rows; i++)
      {
        set_mp(&L0->entry[i][j], &mat0->entry[i][j]);
        set_d(&L1->entry[i][j], &mat1->entry[i][j]);
      }
      // set bottom row to 0
      for (i = mat0->rows; i < mat0->cols; i++)
      {
        set_zero_mp(&L0->entry[i][j]);
        set_zero_d(&L1->entry[i][j]);
      }
    }
    L0->rows = L1->rows = L0->cols = L1->cols = mat0->cols;

    // do the rank revealing decomposition
    QLP_block_pair_amp_mp_d_prec(L0, L1, L0, prec0, L1, rowTop, colTop);
  }

  // find the number of values
  num_FV = L0->rows;

  // allocate for the Frobenius norms of the matrices
  FV[0] = (mpf_t *)bmalloc(num_FV * sizeof(mpf_t));
  FV[1] = (mpf_t *)bmalloc(num_FV * sizeof(mpf_t));

  // setup FV
  setup_frobenius_amp_mp_d(FV, L0, L1);

  // estimate the condition number
  if (mpfr_zero_p(FV[0][num_FV - 1]))
  { // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, FV[0][0], FV[0][num_FV - 1]);
    *CN = mpf_get_d(ratio);
  }

  // look to normalize the entries so that FV[0][0] & FV[1][0] <= 1
  if (mpf_cmp_ui(FV[0][0], 1) > 0 || mpf_cmp_ui(FV[1][0], 1) > 0)
  { // find normalization_factor = max_entry
    mpfr_max(normalization_factor, FV[0][0], FV[1][0], __gmp_default_rounding_mode);

    // normalize FV[0] & FV[1] by multiplying by 1 / normalization_factor
    mpf_ui_div(ratio, 1, normalization_factor);
    for (i = 0; i < num_FV; i++)
    {
      mpf_mul(FV[0][i], FV[0][i], ratio);
      mpf_mul(FV[1][i], FV[1][i], ratio);
    }
  }
  else
  { // normalization_factor = 1
    mpf_set_ui(normalization_factor, 1);
  }

  // initialize corank & bad_loop
  corank = num_FV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_FV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (mpfr_zero_p(FV[0][i]) || mpfr_zero_p(FV[1][i]))
      bad_loop = 1;
    else
    { // find the ratio of the singular values (min / max)
      if (mpf_cmp(FV[0][i], FV[1][i]) < 0)
        mpf_div(ratio, FV[0][i], FV[1][i]);
      else 
        mpf_div(ratio, FV[1][i], FV[0][i]);

      // determine if this is okay or not
      if (mpf_cmp_d(ratio, max_FV_ratio) <= 0)
        bad_loop = 1;
    }

    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    mpfr_max(ratio, FV[0][num_FV - 1], FV[1][num_FV - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    *largest_zero = 0;
  }
  else if (corank == num_FV)
  { // all are zero
    *smallest_nonzero = 0;

    mpfr_max(ratio, FV[0][0], FV[1][0], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }
  else // 0 < corank < num_FV
  { // num_FV - corank - 1 is smallest non-zero & num_FV - corank is largest zero
    mpfr_max(ratio, FV[0][num_FV - corank - 1], FV[1][num_FV - corank - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    mpfr_max(ratio, FV[0][num_FV - corank], FV[1][num_FV - corank], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }

  // release memory
  for (i = num_FV - 1; i >= 0 ; i--)
  {
    mpf_clear(FV[0][i]);
    mpf_clear(FV[1][i]);
  }
  free(FV[0]); free(FV[1]); free(FV);
  mpf_clear(ratio); mpf_clear(normalization_factor);
  clear_mat_mp(L0); clear_mat_d(L1); 

  return corank;
}

