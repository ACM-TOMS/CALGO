// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// uses two different approximations to find either if the matrix is rank deficient or the corank of the matrix

/////////// double precision ////////////

int rankDef_d(double *CN, double minSV0, double minSV1, double mat_norm, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if the matrix   *
* is rank deficient or not, and finds the condition number      *
\***************************************************************/
{
  int rankDef = 0;
  double ratio;

  // check to see if either are exactly 0
  if (minSV0 == 0 || minSV1 == 0)
  { // we know we have rank deficiency
    rankDef = 1;

    // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    *CN = mat_norm / minSV1;

    // find the ratio of the minimum singular values (min / max)
    if (minSV0 < minSV1)
      ratio = minSV0 / minSV1;
    else
      ratio = minSV1 / minSV0;

    // look to normalize (used to compare with SV_tol)
    if (mat_norm > 1)
    { // normalize minSV0 & minSV1
      mat_norm = 1 / mat_norm;

      minSV0 *= mat_norm;
      minSV1 *= mat_norm;
    }

    // determine if we have rank deficiency or not
    if (ratio <= max_SV_ratio || minSV0 < SV_tol || minSV1 < SV_tol || *CN > max_CN)
      rankDef = 1;
    else
      rankDef = 0;
  }

  return rankDef;
}

int corank_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0, mat_d mat1, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero SV*
\***************************************************************/
{
  int i, num_SV, bad_loop, corank = 0, svd_its = 100;
  double ratio, normalization_factor, tol_conv = 1e-15, tol_sign = 1e-20, largeChange = 1e13, **SV = (double **)bmalloc(2 * sizeof(double *));
  vec_d E;

  // make sure that mat0 & mat1 are the same size
  if (mat0->rows != mat1->rows || mat0->cols != mat1->cols)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_vec_d(E, 0);

  // find the singular values for the first matrix
  svd_jacobi_E_d(E, mat0, svd_its, SV_tol, tol_conv, tol_conv, tol_sign, largeChange);

  // find the number of singular values
  num_SV = E->size;

  // setup SV[0] - storing the singular values of mat0
  SV[0] = (double *)bmalloc(num_SV * sizeof(double));
  for (i = 0; i < num_SV; i++)
    SV[0][i] = E->coord[i].r;

  // find the singular values for the second matrix
  svd_jacobi_E_d(E, mat1, svd_its, SV_tol, tol_conv, tol_conv, tol_sign, largeChange);

  // setup SV[1] - storing the singular values of mat1
  SV[1] = (double *)bmalloc(num_SV * sizeof(double));
  for (i = 0; i < num_SV; i++)
    SV[1][i] = E->coord[i].r;

  // find the condition number
  if (SV[0][num_SV - 1] == 0 || SV[1][num_SV - 1] == 0)
  { // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    *CN = SV[1][0] / SV[1][num_SV - 1];
  }

  // look to normalize the entries so that SV[0][0] & SV[1][0] <= 1
  if (SV[0][0] > 1 || SV[1][0] > 1)
  { // find normalization_factor = max_entry
    normalization_factor = MAX(SV[0][0], SV[1][0]);

    // normalize SV[0] & SV[1] by multiplying by 1 / normalization_factor
    ratio = 1 / normalization_factor;
    for (i = 0; i < num_SV; i++)
    {
      SV[0][i] *= ratio;
      SV[1][i] *= ratio;
    }
  }
  else
  { // normalization_factor = 1
    normalization_factor = 1;
  }

  // initialize corank & bad_loop
  corank = num_SV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_SV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (SV[0][i] == 0 || SV[1][i] == 0)
      bad_loop = 1;
    else
    { // find the ratio of the singular values (min / max)
      if (SV[0][i] < SV[1][i])
        ratio = SV[0][i] / SV[1][i];
      else
        ratio = SV[1][i] / SV[0][i];

      // determine if this is okay or not
     if (ratio <= max_SV_ratio || (i == 0 && (!checkGood_d(SV[0][i], 1, SV_tol, largeChange) || !checkGood_d(SV[1][i], 1, SV_tol, largeChange))) || (i > 0 && (!checkGood_d(SV[0][i], SV[0][i-1], SV_tol, largeChange) || !checkGood_d(SV[1][i], SV[1][i-1], SV_tol, largeChange))))
        bad_loop = 1;
    }      

    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    *smallest_nonzero = MAX(SV[0][num_SV - 1], SV[1][num_SV - 1]) * normalization_factor;
    *largest_zero = 0;
  }
  else if (corank == num_SV)
  { // all are zero
    *smallest_nonzero = 0;
    *largest_zero = MAX(SV[0][0], SV[1][0]) * normalization_factor;
  }
  else // 0 < corank < num_SV
  { // num_SV - corank - 1 is smallest non-zero & num_SV - corank is largest zero
    *smallest_nonzero = MAX(SV[0][num_SV - corank - 1], SV[1][num_SV - corank - 1]) * normalization_factor;
    *largest_zero = MAX(SV[0][num_SV - corank], SV[1][num_SV - corank]) * normalization_factor;
  }

  // release memory
  free(SV[0]); free(SV[1]); free(SV);
  clear_vec_d(E);
 
  return corank;
}

/////////// fixed multi precision ////////////

int rankDef_mp(double *CN, mpf_t minSV0, mpf_t minSV1, mpf_t mat_norm, int curr_prec, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if the matrix   *
* is rank deficient or not, and finds the condition number      *
\***************************************************************/
{
  int rankDef = 0;
  mpf_t ratio;

  // initialize
  mpf_init(ratio);

  // check to see if either are exactly 0
  if (mpfr_zero_p(minSV0) || mpfr_zero_p(minSV1))
  { // we know we have rank deficiency
    rankDef = 1;
  
    // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, mat_norm, minSV1);
    *CN = mpf_get_d(ratio);

    // find the ratio of the minimum singular values (min / max)
    if (mpf_cmp(minSV0, minSV1) < 0)
      mpf_div(ratio, minSV0, minSV1);
    else
      mpf_div(ratio, minSV1, minSV0);

    // look to normalize (used to compare with SV_tol)
    if (mpf_cmp_ui(mat_norm, 1) > 0)
    { // normalize minSV0 & minSV1
      mpf_ui_div(mat_norm, 1, mat_norm);
   
      mpf_mul(minSV0, minSV0, mat_norm);
      mpf_mul(minSV1, minSV1, mat_norm);
    }

    // determine if we have rank deficiency or not
    if (mpf_cmp_d(ratio, max_SV_ratio) <= 0 || mpf_cmp_d(minSV0, SV_tol) < 0 || mpf_cmp_d(minSV1, SV_tol) < 0 || *CN > max_CN)
      rankDef = 1;
    else
      rankDef = 0;
  }

  // clear
  mpf_clear(ratio);

  return rankDef;
}

int corank_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, mat_mp mat1, int curr_prec, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero SV*
\***************************************************************/
{
  int i, num_SV, bad_loop, corank = 0, svd_its = 100, prec_digits = prec_to_digits(curr_prec) - 3;
  size_t size;
  char *str = NULL;
  mpf_t one, ratio, normalization_factor, SV_tol_mp, tol_conv, tol_sign, largeChange, **SV = (mpf_t **)bmalloc(2 * sizeof(mpf_t *));
  vec_mp E;

  // make sure that mat0 & mat1 are the same size
  if (mat0->rows != mat1->rows || mat0->cols != mat1->cols)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize MP
  mpf_init(one); mpf_init(ratio); mpf_init(normalization_factor); mpf_init(SV_tol_mp);
  mpf_init(tol_conv); mpf_init(tol_sign); mpf_init(largeChange);
  mpf_set_ui(one, 1);
  init_vec_mp(E, 0); 

  // setup the tolerances 
  mpf_set_d(SV_tol_mp, SV_tol);

  size = 1 + snprintf(NULL, 0, "1e-%d", prec_digits + 1);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", prec_digits + 1);
  mpf_set_str(tol_conv, str, 10);

  size = 1 + snprintf(NULL, 0, "1e-%d", 2 * prec_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", 2 * prec_digits);
  mpf_set_str(tol_sign, str, 10);

  size = 1 + snprintf(NULL, 0, "1e%d", prec_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", prec_digits);
  mpf_set_str(largeChange, str, 10); 

  // find the singular values for the first matrix
  svd_jacobi_E_mp(E, mat0, svd_its, SV_tol_mp, tol_conv, tol_conv, tol_sign, largeChange);

  // find the number of singular values
  num_SV = E->size;

  // setup SV[0] - storing the singular values of mat0
  SV[0] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
  for (i = 0; i < num_SV; i++)
  {
    mpf_init(SV[0][i]);
    mpf_set(SV[0][i], E->coord[i].r);
  }

  // find the singular values for the second matrix  
  svd_jacobi_E_mp(E, mat1, svd_its, SV_tol_mp, tol_conv, tol_conv, tol_sign, largeChange);

  // setup SV[1] - storing the singular values of mat1
  SV[1] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
  for (i = 0; i < num_SV; i++)
  {
    mpf_init(SV[1][i]);
    mpf_set(SV[1][i], E->coord[i].r);
  }

  // find the condition number
  if (mpfr_zero_p(SV[0][num_SV - 1]) || mpfr_zero_p(SV[1][num_SV - 1]))
  { // setup CN
    *CN = max_CN * 10;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, SV[1][0], SV[1][num_SV - 1]);
    *CN = mpf_get_d(ratio);
  }

  // look to normalize the entries so that SV[0][0] & SV[1][0] <= 1
  if (mpf_cmp_ui(SV[0][0], 1) > 0 || mpf_cmp_ui(SV[1][0], 1) > 0)
  { // find normalization_factor = max_entry
    mpfr_max(normalization_factor, SV[0][0], SV[1][0], __gmp_default_rounding_mode);

    // normalize SV[0] & SV[1] by multiplying by 1 / normalization_factor
    mpf_div(ratio, one, normalization_factor);
    for (i = 0; i < num_SV; i++)
    {
      mpf_mul(SV[0][i], SV[0][i], ratio);
      mpf_mul(SV[1][i], SV[1][i], ratio);
    }
  }
  else
  { // normalization_factor = 1
    mpf_set(normalization_factor, one);
  }

  // initialize corank & bad_loop
  corank = num_SV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_SV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (mpfr_zero_p(SV[0][i]) || mpfr_zero_p(SV[1][i]))
      bad_loop = 1;
    else
    { // find the ratio of the singular values (min / max)
      if (mpf_cmp(SV[0][i], SV[1][i]) < 0)
        mpf_div(ratio, SV[0][i], SV[1][i]);
      else 
        mpf_div(ratio, SV[1][i], SV[0][i]);

      // determine if this is okay or not
      if (mpf_cmp_d(ratio, max_SV_ratio) <= 0 || (i == 0 && (!checkGood_mp(SV[0][i], one, SV_tol_mp, largeChange) || !checkGood_mp(SV[1][i], one, SV_tol_mp, largeChange))) || (i > 0 && (!checkGood_mp(SV[0][i], SV[0][i-1], SV_tol_mp, largeChange) || !checkGood_mp(SV[1][i], SV[1][i-1], SV_tol_mp, largeChange))))
        bad_loop = 1;
    }

    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    mpfr_max(ratio, SV[0][num_SV - 1], SV[1][num_SV - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    *largest_zero = 0;
  }
  else if (corank == num_SV)
  { // all are zero
    *smallest_nonzero = 0;

    mpfr_max(ratio, SV[0][0], SV[1][0], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }
  else // 0 < corank < num_SV
  { // num_SV - corank - 1 is smallest non-zero & num_SV - corank is largest zero
    mpfr_max(ratio, SV[0][num_SV - corank - 1], SV[1][num_SV - corank - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    mpfr_max(ratio, SV[0][num_SV - corank], SV[1][num_SV - corank], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }

  // clear MP
  for (i = num_SV - 1; i >= 0; i--)
  {
    mpf_clear(SV[0][i]);
    mpf_clear(SV[1][i]);
  }
  free(SV[0]); free(SV[1]);
  free(SV);

  mpf_clear(one); mpf_clear(ratio); mpf_clear(normalization_factor); mpf_clear(SV_tol_mp);
  mpf_clear(tol_conv); mpf_clear(tol_sign); mpf_clear(largeChange);
  clear_vec_mp(E);

  free(str);

  return corank;
}

/////////// adaptive precision ////////////

int rankDef_amp(double *CN, mpf_t minSV0, int prec0, mpf_t minSV1, int prec1, mpf_t mat_norm, double max_SV_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - rank deficient, 0 - not rank deficient     *
* NOTES: uses the 2 approximations to determine if the matrix   *
* is rank deficient or not, and finds the condition number      *
\***************************************************************/
{
  int rankDef = 0, max_prec = MAX(prec0, prec1);
  mpf_t ratio;

  mpf_init2(ratio, max_prec);

  // check to see if either are exactly 0
  if (mpfr_zero_p(minSV0) || mpfr_zero_p(minSV1))
  { // we know we hav rank deficiency
    rankDef = 1;

    // setup CN
    *CN = 1e99;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, mat_norm, minSV1);
    *CN = mpf_get_d(ratio);

    // find the ratio of the minimum singular values (min / max)
    if (mpf_cmp(minSV0, minSV1) < 0)
      mpf_div(ratio, minSV0, minSV1);
    else
      mpf_div(ratio, minSV1, minSV0);

    // determine if the minimum singular value has changed too much
    if (mpf_cmp_d(ratio, max_SV_ratio) <= 0)
      rankDef = 1;
    else
      rankDef = 0;
  }

  mpf_clear(ratio);

  return rankDef;
}

int corank_amp(double *CN, double *smallest_nonzero, double *largest_zero, mat_d mat0_d, mat_mp mat0_mp, int prec0, mat_d mat1_d, mat_mp mat1_mp, int prec1, double max_SV_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: corank = number of 'zero' singular values      *
* NOTES: uses the 2 approximations to determine the corank for  *
* the matrix, also find CN & smallest non-zero & largest zero SV*
\***************************************************************/
{
  int i, num_SV, bad_loop, mat0r, mat0c, mat1r, mat1c, corank = 0, svd_its = 100, prec_digits = 0, max_prec = MAX(prec0, prec1);
  size_t size;
  char *str = NULL;
  double tol_conv_d = 1e-15, tol_sign_d = 1e-20, largeChange_d = 1e13;
  mpf_t ratio, normalization_factor, tol_conv, tol_sign, largeChange, **SV = (mpf_t **)bmalloc(2 * sizeof(mpf_t *));
  mpf_t one, good_tol, good_largeChange;
  vec_d E_d;
  vec_mp E_mp;

  // make sure mat_prec >= 64
  max_prec = MAX(max_prec, 64);

  // setup mat0 sizes
  if (prec0 < 64)
  {
    mat0r = mat0_d->rows;
    mat0c = mat0_d->cols;
  } 
  else
  {
    mat0r = mat0_mp->rows;
    mat0c = mat0_mp->cols;
  } 

  // setup mat1 sizes
  if (prec1 < 64)
  {
    mat1r = mat1_d->rows;
    mat1c = mat1_d->cols;
  }
  else
  {
    mat1r = mat1_mp->rows;
    mat1c = mat1_mp->cols;
  }

  // make sure that mat0 & mat1 are the same size
  if (mat0r != mat1r || mat0c != mat1c)
  {
    printf("ERROR: To find corank, the matrix sizes need to be the same!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  mpf_init2(ratio, max_prec); mpf_init2(normalization_factor, max_prec);
  mpf_init2(tol_conv, max_prec); mpf_init2(tol_sign, max_prec); mpf_init2(largeChange, max_prec);
  mpf_init2(one, max_prec); mpf_init2(good_tol, max_prec); mpf_init2(good_largeChange, max_prec);
  mpf_set_ui(one, 1);
  init_vec_d(E_d, 0); init_vec_mp2(E_mp, 0, max_prec);

  // find the singular values for the first matrix
  if (prec0 < 64)
  { // compute using mat0_d
    mpf_set_d(good_tol, tol_conv_d);
    mpf_set_d(good_largeChange, largeChange_d);
    svd_jacobi_E_d(E_d, mat0_d, svd_its, tol_conv_d, tol_conv_d, tol_conv_d, tol_sign_d, largeChange_d);

    // find the number of singular values
    num_SV = E_d->size;

    // setup SV[0] - storing the singular values of mat0
    SV[0] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
    for (i = 0; i < num_SV; i++)
    {
      mpf_init2(SV[0][i], max_prec);
      mpf_set_d(SV[0][i], E_d->coord[i].r);
    }
  }
  else
  { // compute using mat0_mp
    initMP(prec0);
    prec_digits = prec_to_digits(prec0) - 3;

    // set the tolerances
    size = 1 + snprintf(NULL, 0, "1e-%d", prec_digits + 1);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", prec_digits + 1);
    mpf_set_str(tol_conv, str, 10);

    size = 1 + snprintf(NULL, 0, "1e-%d", 2 * prec_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", 2 * prec_digits);
    mpf_set_str(tol_sign, str, 10);

    size = 1 + snprintf(NULL, 0, "1e%d", prec_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", prec_digits);
    mpf_set_str(largeChange, str, 10);

    mpf_set(good_tol, tol_conv);
    mpf_set(good_largeChange, largeChange);

    svd_jacobi_E_mp(E_mp, mat0_mp, svd_its, tol_conv, tol_conv, tol_conv, tol_sign, largeChange);

    // find the number of singular values
    num_SV = E_mp->size;

    // setup SV[0] - storing the singular values of mat0
    SV[0] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
    for (i = 0; i < num_SV; i++)
    {
      mpf_init2(SV[0][i], max_prec);
      mpf_set(SV[0][i], E_mp->coord[i].r);
    }
  }

  // find the singular values for the second matrix 
  if (prec1 < 64)
  { // compute using mat1_d
    svd_jacobi_E_d(E_d, mat1_d, svd_its, tol_conv_d, tol_conv_d, tol_conv_d, tol_sign_d, largeChange_d);

    // setup SV[1] - storing the singular values of mat0
    SV[1] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
    for (i = 0; i < num_SV; i++)
    {
      mpf_init2(SV[1][i], max_prec);
      mpf_set_d(SV[1][i], E_d->coord[i].r);
    }
  }
  else
  { // compute using mat1_mp
    initMP(prec1);
    prec_digits = prec_to_digits(prec1) - 3;

    // set the tolerances
    size = 1 + snprintf(NULL, 0, "1e-%d", prec_digits + 1);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", prec_digits + 1);
    mpf_set_str(tol_conv, str, 10);

    size = 1 + snprintf(NULL, 0, "1e-%d", 2 * prec_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e-%d", 2 * prec_digits);
    mpf_set_str(tol_sign, str, 10);

    size = 1 + snprintf(NULL, 0, "1e%d", prec_digits);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "1e%d", prec_digits);
    mpf_set_str(largeChange, str, 10);

    svd_jacobi_E_mp(E_mp, mat1_mp, svd_its, tol_conv, tol_conv, tol_conv, tol_sign, largeChange);

    // setup SV[1] - storing the singular values of mat1
    SV[1] = (mpf_t *)bmalloc(num_SV * sizeof(mpf_t));
    for (i = 0; i < num_SV; i++)
    {
      mpf_init2(SV[1][i], max_prec);
      mpf_set(SV[1][i], E_mp->coord[i].r);
    }
  }

  // finish setup if prec1 is bigger
  if (prec0 < prec1)
  { // prec1 >= 64 so use mpf
    mpf_set(good_tol, tol_conv);
    mpf_set(good_largeChange, largeChange);
  }

  // find the condition number
  if (mpfr_zero_p(SV[0][num_SV - 1]) || mpfr_zero_p(SV[1][num_SV - 1]))
  { // setup CN
    *CN = 1e199;
  }
  else
  { // calculate CN as expected
    mpf_div(ratio, SV[1][0], SV[1][num_SV - 1]);
    *CN = mpf_get_d(ratio);
  }

  // look to normalize the entries so that SV[0][0] & SV[1][0] <= 1
  if (mpf_cmp(SV[0][0], one) > 0 || mpf_cmp(SV[1][0], one) > 0)
  { // find normalization_factor = max_entry
    mpfr_max(normalization_factor, SV[0][0], SV[1][0], __gmp_default_rounding_mode);

    // normalize SV[0] & SV[1] by multiplying by 1 / normalization_factor
    mpf_div(ratio, one, normalization_factor);
    for (i = 0; i < num_SV; i++)
    {
      mpf_mul(SV[0][i], SV[0][i], ratio);
      mpf_mul(SV[1][i], SV[1][i], ratio);
    }
  }
  else
  { // normalization_factor = 1
    mpf_set(normalization_factor, one);
  }

  // initialize corank & bad_loop
  corank = num_SV;
  bad_loop = 0;

  // once a bad one has been found - all the rest below it are bad as well
  for (i = 0; i < num_SV && !bad_loop; i++)
  { // make sure that they are non-zero
    if (mpfr_zero_p(SV[0][i]) || mpfr_zero_p(SV[1][i]))
      bad_loop = 1;
    else
    { // find the ratio of the singular values (min / max)
      if (mpf_cmp(SV[0][i], SV[1][i]) > 0)
        mpf_div(ratio, SV[1][i], SV[0][i]);
      else
        mpf_div(ratio, SV[0][i], SV[1][i]);

      // determine if this is okay or not
      if (mpf_cmp_d(ratio, max_SV_ratio) <= 0 || (i == 0 && (!checkGood_mp(SV[0][i], one, good_tol, good_largeChange) || !checkGood_mp(SV[1][i], one, good_tol, good_largeChange))) || (i > 0 && (!checkGood_mp(SV[0][i], SV[0][i-1], good_tol, good_largeChange) || !checkGood_mp(SV[1][i], SV[1][i-1], good_tol, good_largeChange))))
        bad_loop = 1;
    }
    // check to see if we had a good loop
    if (!bad_loop)
      corank--;
  }

  // setup smallest_nonzero & largest_zero
  if (corank == 0)
  { // all are non-zero
    mpfr_max(ratio, SV[0][num_SV - 1], SV[1][num_SV - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    *largest_zero = 0;
  }
  else if (corank == num_SV)
  { // all are zero
    *smallest_nonzero = 0;

    mpfr_max(ratio, SV[0][0], SV[1][0], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }
  else // 0 < corank < num_SV
  { // num_SV - corank - 1 is smallest non-zero & num_SV - corank is largest zero
    mpfr_max(ratio, SV[0][num_SV - corank - 1], SV[1][num_SV - corank - 1], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *smallest_nonzero = mpf_get_d(ratio);

    mpfr_max(ratio, SV[0][num_SV - corank], SV[1][num_SV - corank], __gmp_default_rounding_mode);
    mpf_mul(ratio, ratio, normalization_factor);
    *largest_zero = mpf_get_d(ratio);
  }

  // clear
  for (i = num_SV - 1; i >= 0; i--)
  {
    mpf_clear(SV[0][i]);
    mpf_clear(SV[1][i]);
  }
  free(SV[0]); free(SV[1]);
  free(SV);

  mpf_clear(ratio); mpf_clear(normalization_factor); 
  mpf_clear(tol_conv); mpf_clear(tol_sign); mpf_clear(largeChange);
  mpf_clear(one); mpf_clear(good_tol); mpf_clear(good_largeChange);
  clear_vec_d(E_d); clear_vec_mp(E_mp);

  free(str);

  return corank;
}


