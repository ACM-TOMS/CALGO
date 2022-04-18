// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

// this file contains all of the different types of matrixSolve functions that bertini can use
// matrixSolve_d,  matrixSolve2_d,  matrixSolve_cond_num_norms_d,  LU_matrixSolve_d,  matrixSolve_from_LU_d
// matrixSolve_mp, matrixSolve2_mp, matrixSolve_cond_num_norms_mp, LU_matrixSolve_mp, matrixSolve_from_LU_mp

#define MS_NOSOLUTION -1
#define _VERYSMALL 1e-14 // double precision can be taken as accurate to 15 digits and with round off errors, we can assume that we should be accurate within 14 digits
#define _VERYLARGECHANGE 1e11  // if the pivot entries change by a factor of this, we call the matrix singular - looking for a sharp decline in singular values

/****** DOUBLE PRECISION ********/

int matrixSolve_d(vec_d x, mat_d A, vec_d b)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, m = A->rows, n = A->cols;
  double c, nA, nA_inv;

  // check that we have atleast as many rows as columns and that the size of b is correct
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // if square, try to use LU solving first
  if (m == n)
  {
    retVal = matrixSolve_LU_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE);
    if (retVal) // if LU solving failed, try QR solving to see if really rank deficient
    {
      retVal = matrixSolve_Least_Squares_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE, &c, &nA, &nA_inv); // retVal == corank of the matrix

      if (retVal != 0) // if we do not have a full rank matrix, we return no solution
        retVal = MS_NOSOLUTION;
    }
  }
  else // m > n so we need to use least squares solution
  {
    retVal = matrixSolve_Least_Squares_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE, &c, &nA, &nA_inv); // retVal == corank of the matrix      
  }

  return retVal;
}   

int matrixSolve_cond_num_norms_d(vec_d x, mat_d A, vec_d b, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES: does matrixSolve & condition number together!!         *
\***************************************************************/
{
  int retVal, m = A->rows, n = A->cols;

  // check that we have atleast as many rows as columns and that the size of b is correct
  if (m < n)
  {
    printf("ERROR: The matrix has more columns (%d) than rows (%d) in matrixSolve_cond_num_norms_d!!\n", n, m);
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_cond_num_norms_d (%d x %d and %d) do not match!!\n", A->rows, A->cols, b->size);
    bexit(ERROR_INVALID_SIZE);
  }

  // if square, try to use LU solving first
  if (m == n)
  {
    retVal = matrixSolve_LU_cond_num_norms_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE, cond_num, norm_A, norm_A_inv);
    if (retVal) // if LU solving failed, try QR solving to see if really rank deficient
    {
      retVal = matrixSolve_Least_Squares_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE, cond_num, norm_A, norm_A_inv); // retVal == corank of the matrix
      if (retVal != 0)
        retVal = MS_NOSOLUTION;
    }
  }
  else // m > n so we need to use least squares solution
  {
    retVal = matrixSolve_Least_Squares_d(x, A, b, _VERYSMALL, _VERYLARGECHANGE, cond_num, norm_A, norm_A_inv); // retVal == corank of the matrix
  }

  return retVal;
}

int matrixSolve_LU_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/*
 (JDH - 06/19/06) We shall assume that A is n by n and b is n by 1.
 We solve for x where A*x = b.  This is done using an implementation of
 Gauss Elimination with Scaled Partial Pivoting as described in
 Numerical Methods, Software, and Analysis, 2nd Ed, John R. Rice,
 algorithms 6.2.1 and 6.2.4.
*/
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL;
  double max, tempD, tempD2, prevNorm = 1e300, *scale = NULL;

  vec_d tempVec;
  mat_d tempMat;
  comp_d multiplier, tempComp;

  // error checking
  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in matrixSolve_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup rownum, scale
  rownum = (int *)bmalloc(n * sizeof(int));
  scale = (double *)bmalloc(n * sizeof(double));

  // setup x
  change_size_vec_d(x, n);
  x->size = n;

  // copy A to tempMat & b to tempVec
  init_mat_d(tempMat, n, n);
  init_vec_d(tempVec, n);
  tempMat->rows = tempMat->cols = tempVec->size = n;

  // find scale factors for each row and assign row numbers
  prevNorm = 0;
  for (i = 0; i < n; i++)
  { // setup tempVec
    set_d(&tempVec->coord[i], &b->coord[i]);
    // assign row numbers
    rownum[i] = i;
    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    set_d(&tempMat->entry[i][0], &A->entry[i][0]);
    max = d_oneNorm_d(&tempMat->entry[i][0]);
    for (j = 1; j < n; j++)
    {
      set_d(&tempMat->entry[i][j], &A->entry[i][j]);
      tempD = d_oneNorm_d(&tempMat->entry[i][j]);
      if (max < tempD)
        max = tempD;
    }
    if (max < tol)
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }

    scale[i] = 1 / max; // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find sacled maximum element in column k
    max = d_oneNorm_d(&tempMat->entry[rownum[k]][k]) * scale[rownum[k]];
    for (i = k + 1; i < n; i++)
    {
      tempD = d_oneNorm_d(&tempMat->entry[rownum[i]][k]) * scale[rownum[i]];
      if (max < tempD)
      {
        pivot = i;
        max = tempD;
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    if (max < tol || (k > 0 && prevNorm > max * largeChange))
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }
    else
      prevNorm = max;

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    // update elements in "bottom" block of tempMat
    pivot = rownum[k];
    recip_d2(tempComp, &tempMat->entry[pivot][k], tempD);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];
      // compute muliplier
      mul_d2(multiplier, &tempMat->entry[l][k], tempComp, tempD);
      // update tempVec (RHS)
      sub_mul_d2(&tempVec->coord[l], multiplier, &tempVec->coord[pivot], tempD);
      // update elements of row i
      for (j = k + 1; j < n; j++)
        sub_mul_d2(&tempMat->entry[l][j], multiplier, &tempMat->entry[pivot][j], tempD);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  max = d_oneNorm_d(&tempMat->entry[pivot][l]) * scale[pivot];
  if (max < tol || prevNorm > max * largeChange)
  { // clear
    clear_mat_d(tempMat);
    clear_vec_d(tempVec);
    free(rownum);
    free(scale);

    return MS_NOSOLUTION;
  }

  // calculate x->coord[n-1]
  div_d2(&x->coord[l], &tempVec->coord[pivot], &tempMat->entry[pivot][l], tempD, tempD2);
  for (i = n - 2; i >= 0; i--)
  { // calculate x->coord[i]
    k = rownum[i];
    set_d(tempComp, &tempVec->coord[k]);
    for (j = i + 1; j < n; j++)
      sub_mul_d2(tempComp, &tempMat->entry[k][j], &x->coord[j], tempD);
    div_d2(&x->coord[i], tempComp, &tempMat->entry[k][i], tempD, tempD2);
  }

  // check answer
  for (i = 0; i < n; i++)
    if (isnan(x->coord[i].r) || isnan(x->coord[i].i))
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }

  // clear
  clear_mat_d(tempMat);
  clear_vec_d(tempVec);
  free(rownum);
  free(scale);

  return 0;
}

int matrixSolve2_LU_d(vec_d x1, vec_d x2, mat_d A, vec_d b1, vec_d b2, double tol, double largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
// Solve A * x = b, for both x = x1, b = b1 and x = x2, b = b2 where A - n x n, b - n x 1, x - n x 1
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL;
  double max, tempD, tempD2, prevNorm = 1e300, *scale = NULL;

  vec_d tempVec1, tempVec2;
  mat_d  tempMat;
  comp_d multiplier, tempComp1, tempComp2, tempComp3;

  // error checking
  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in matrixSolve2_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b2->size || A->rows != b2->size)
  {
    printf("ERROR: The sizes in matrixSolve2_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup rownum, scale
  rownum = (int *)bmalloc(n * sizeof(int));
  scale = (double *)bmalloc(n * sizeof(double));

  // setup x1 & x2
  change_size_vec_d(x1, n);
  change_size_vec_d(x2, n);
  x1->size = x2->size = n;

  // copy A to tempMat & b1 to tempVec1 & b2 to tempVec2
  init_mat_d(tempMat, n, n);
  init_vec_d(tempVec1, n);
  init_vec_d(tempVec2, n);
  tempMat->rows = tempMat->cols = tempVec1->size = tempVec2->size = n;

  // find scale factors for each row and assign row numbers
  for (i = 0; i < n; i++)
  { // setup tempVec1 & tempVec2
    set_d(&tempVec1->coord[i], &b1->coord[i]);
    set_d(&tempVec2->coord[i], &b2->coord[i]);

    // assign row numbers
    rownum[i] = i;

    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    set_d(&tempMat->entry[i][0], &A->entry[i][0]);
    max = d_oneNorm_d(&tempMat->entry[i][0]);
    for (j = 1; j < n; j++)
    {
      set_d(&tempMat->entry[i][j], &A->entry[i][j]);
      tempD = d_oneNorm_d(&tempMat->entry[i][j]);
      if (max < tempD)
        max = tempD;
    }
    if (max < tol)
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec1);
      clear_vec_d(tempVec2);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }

    scale[i] = 1 / max; // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;
    max = d_oneNorm_d(&tempMat->entry[rownum[k]][k]) * scale[rownum[k]];
    // find scaled maximum element in column k
    for (i = k + 1; i < n; i++)
    {
      tempD = d_oneNorm_d(&tempMat->entry[rownum[i]][k]) * scale[rownum[i]];
      if (max < tempD)
      {
        pivot = i;
        max = tempD;
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    if (max < tol || (k > 0 && prevNorm > max * largeChange))
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec1);
      clear_vec_d(tempVec2);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }
    else
      prevNorm = max;

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    pivot = rownum[k];
    recip_d2(tempComp3, &tempMat->entry[pivot][k], tempD);
    // update elements in "bottom" block of tempMat
    for (i = k + 1; i < n; i++)
    { // find the row number 
      l = rownum[i];
      // compute the multiplier
      mul_d2(multiplier, &tempMat->entry[l][k], tempComp3, tempD);
      // update tempVec1 & tempVec2 (the RHS)
      sub_mul_d2(&tempVec1->coord[l], multiplier, &tempVec1->coord[pivot], tempD);
      sub_mul_d2(&tempVec2->coord[l], multiplier, &tempVec2->coord[pivot], tempD);
      // update elements of row l
      for (j = k + 1; j < n; j++)
        sub_mul_d2(&tempMat->entry[l][j], multiplier, &tempMat->entry[pivot][j], tempD);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  max = d_oneNorm_d(&tempMat->entry[pivot][l]) * scale[pivot];
  if (max < tol || prevNorm > max * largeChange)
  { // clear
    clear_mat_d(tempMat);
    clear_vec_d(tempVec1);
    clear_vec_d(tempVec2);
    free(rownum);
    free(scale);

    return MS_NOSOLUTION;
  }

  // calculate x1->coord[n-1] & x2->coord[n-1]
  recip_d2(tempComp3, &tempMat->entry[pivot][l], tempD);
  mul_d2(&x1->coord[l], &tempVec1->coord[pivot], tempComp3, tempD);
  mul_d2(&x2->coord[l], &tempVec2->coord[pivot], tempComp3, tempD);
  // compute the rest of the solution vectors
  for (i = n - 2; i >= 0; i--)
  { // calculate x1->coord[i] & x2->coord[i]
    k = rownum[i];
    set_d(tempComp1, &tempVec1->coord[k]);
    set_d(tempComp2, &tempVec2->coord[k]);
    for (j = i + 1; j < n; j++)
    {
      sub_mul_d2(tempComp1, &tempMat->entry[k][j], &x1->coord[j], tempD);
      sub_mul_d2(tempComp2, &tempMat->entry[k][j], &x2->coord[j], tempD);
    }
    div_d2(&x1->coord[i], tempComp1, &tempMat->entry[k][i], tempD, tempD2);
    div_d2(&x2->coord[i], tempComp2, &tempMat->entry[k][i], tempD, tempD2);
  }

  // check answer
  for (i = 0; i < n; i++)
    if (isnan(x1->coord[i].r) || isnan(x1->coord[i].i) || isnan(x2->coord[i].r) || isnan(x2->coord[i].i))
    { // clear
      clear_mat_d(tempMat);
      clear_vec_d(tempVec1);
      clear_vec_d(tempVec2);
      free(rownum);
      free(scale);

      return MS_NOSOLUTION;
    }

  // clear
  clear_mat_d(tempMat);
  clear_vec_d(tempVec1);
  clear_vec_d(tempVec2);
  free(rownum);
  free(scale);

  return 0;
}

int matrixSolve_LU_cond_num_norms_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does matrixSolve & condition number together!!         *
\***************************************************************/
{
  int i, numCols = A->cols, retVal = 0;;
  double maxNorm, tempNorm;

  if ((A->rows == 1) && (A->cols == 1))
  {
    *norm_A = d_abs_d(&A->entry[0][0]);
    *norm_A_inv = 1.0 / (*norm_A);
    *cond_num = 1.0;
    retVal = matrixSolve_d(x, A, b);
  }
  else
  {
    vec_d randVec;
    init_vec_d(randVec, numCols);
    randVec->size = numCols;

    // obtain norm of A
    *norm_A = infNormMat_d(A);

    // Estimate ||A^{-1}|| by :
    //   ||A^{-1}|| ~ 2||A^{-1}v|| , for a random unit vector v.

    maxNorm = 0.0;
    for (i = 0; i < numCols; i++)
    {
      // obtain complex random number and its modulus
      tempNorm = get_comp_rand_d(&randVec->coord[i]);

      // determine if it is the largest modulus so far
      if (maxNorm < tempNorm)
        maxNorm = tempNorm;
    }

    // update maxNorm to be 1 / max modulus = 1 / infNormVec_d(randVec)
    maxNorm = 1.0 / maxNorm;

    // solve for x and 'randVec' where A*x = b and A * 'randVec' = randVec
    retVal = matrixSolve2_LU_d(x, randVec, A, b, randVec, tol, largeChange);

    if (retVal)
    { // matrixSolve2_d failed
      *norm_A_inv = -1.0;
      *cond_num = -1.0;
    }
    else
    { // correct result by multiplying by maxNorm - makes the input randVec a unit vector
      *norm_A_inv = 2.0 * infNormVec_d(randVec) * maxNorm;
      *cond_num = (*norm_A) * (*norm_A_inv);
    }

    clear_vec_d(randVec);
  }

  return retVal;
}

int matrixSolve_svd_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using SVD!!                                 *
\***************************************************************/
{
  int i, j, retVal, rank, its = 50, m = A->rows, n = A->cols;
  double tempD, tol_prec = 1e-15, tol_sign = 1e-18;
  comp_d tempComp;
  mat_d U, E, V;
  vec_d tempVec;

  // do error checking
  if (tol < 0)
  {
    printf("ERROR: The tolerance for matrixSolve_svd_Least_Squares_d needs to be non-negative!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (largeChange < 0)
  {
    printf("ERROR: The change tolerance for matrixSolve_svd_Least_Squares_d needs to be non-negative!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_svd_Least_Squares_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_svd_Least_Squares_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup the size of x
  change_size_vec_d(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    *cond_num = *norm_A = *norm_A_inv = 1;
    return 0;
  }
  // so we can assume that m >= n >= 1

  init_vec_d(tempVec, 0);
  init_mat_d(U, 0, 0);
  init_mat_d(E, 0, 0);
  init_mat_d(V, 0, 0);

  // find the SVD of A using at most its iterations with tolerances of
  // tol, tol, tol_sign, and largeChange
  retVal = svd_jacobi_d(U, E, V, A, its, tol, tol_prec, tol_prec, tol_sign, largeChange);

  if (retVal >= 0)
  { // the svd is good with corank == retVal, so we find x = V * (E^-1)^H * U^H * b

    // find the rank
    rank = n - retVal;

    // set the size of tempVec
    change_size_vec_d(tempVec, rank);

    // find (E^-1)^H * U^H * b - take into account the rank to save operations
    for (j = 0; j < rank; j++)
    {
      set_zero_d(&tempVec->coord[j]);
      for (i = 0; i < m; i++)
      {
        conjugate_d(tempComp, &U->entry[i][j]);
        sum_mul_d(&tempVec->coord[j], tempComp, &b->coord[i]);
      }
      tempD = 1 / E->entry[j][j].r ; // E[j][j] is real and large enough to be inverted
      mul_rdouble_d(&tempVec->coord[j], &tempVec->coord[j], tempD);
    }

    // find V * tempVec - take into account the rank to save operations
    for (i = 0; i < n; i++)
    {
      set_zero_d(&x->coord[i]);
      for (j = 0; j < rank; j++)
      {
        sum_mul_d(&x->coord[i], &V->entry[i][j], &tempVec->coord[j]);
      }
    }

    // find the norm of A - the largest singular value
    *norm_A = E->entry[0][0].r;
    // find the norm of A^-1 - the inverse of the smallest singular value
    tempD = E->entry[n-1][n-1].r;
    if (tempD > 0)
    {
      *norm_A_inv = 1 / tempD; 
      *cond_num = (*norm_A) * (*norm_A_inv);
    }
    else
    { // smallest singular value is 0
      *cond_num = *norm_A_inv = -1;
    }
  }
  else // give up
  {
    *norm_A = infNormMat_d(A);
    *cond_num = *norm_A_inv = -1;
    retVal = MS_NOSOLUTION;
  }

  // clear
  clear_vec_d(tempVec);
  clear_mat_d(U);
  clear_mat_d(E);
  clear_mat_d(V);

  return retVal;
}

int matrixSolve_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using QR!!                                  *
\***************************************************************/
{
  int i, j, k, pivot, m = A->rows, n = A->cols, *perm = NULL;
  double max, norm_y, prevNorm, *norm = NULL;
  comp_d tempComp;
  vec_d tempVec, u, z, y;
  mat_d tempMat;

  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_Least_Squares_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_Least_Squares_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup the size of x
  change_size_vec_d(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    *cond_num = *norm_A = *norm_A_inv = 1;
    return 0;
  }
  // so we can assume that m >= n >= 1

  init_vec_d(u, 0);
  init_vec_d(z, 0);
  init_vec_d(y, 0);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));

  // copy A to tempMat & b to tempVec, find norms for each col and find the one that is the maximum
  init_vec_d(tempVec, m);
  init_mat_d(tempMat, m, n);
  tempMat->rows = tempVec->size = m;
  tempMat->cols = n;
  max = pivot = 0;
  for (j = 0; j < m; j++)
  {
    // copy b to tempVec
    set_d(&tempVec->coord[j], &b->coord[j]);

    if (j < n)
    { // initialize perm
      perm[j] = j;

      // find the twoNorm for each column
      norm[j] = norm_y = 0;
      for (i = 0; i < m; i++)
      { // copy A to tempMat
        set_d(&tempMat->entry[i][j], &A->entry[i][j]);
        // find the size to this entry
        norm[j] += A->entry[i][j].r * A->entry[i][j].r + A->entry[i][j].i * A->entry[i][j].i;
      }
      norm[j] = sqrt(norm[j]);
      // check to see if it is the maximum twoNorm thus far
      if (max < norm[j])
      {
        pivot = j;
        max = norm[j];
      }
    }
  }
  prevNorm = max;
  // take the norm_A as the sqrt(n) * maximum of the twoNorms, which a good enough estimate of the norm of A
  *norm_A = sqrt(n) * max;

  // initialize the I.C.E.
  y->size = 0;
  norm_y = 0;

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  { 
    if (norm[perm[pivot]] > 0)
    { // setup tempComp for I.C.E.
      sign_d(tempComp, &tempMat->entry[k][perm[pivot]], tol);
      mul_rdouble_d(tempComp, tempComp, -norm[perm[pivot]]);
      // find the I.C.E.
      norm_y = cond_est_d(y, y, norm_y, tempMat, 0, k, pivot, tempComp, perm, tol);
    }
    else
    { // next pivot is 0 and so there is no use trying to find the next singular value since it is 0!! - the next if statement will stop the loop
      norm_y = -1; // meaning that the matrix A is singulare & so norm_A_inv = infinity!
    }

    // check to make sure pivot element is not too small
    if (norm[perm[pivot]] <= tol || (k > 0 && prevNorm > norm[perm[pivot]] * largeChange))
      break; 
    else
      prevNorm = norm[perm[pivot]];

    // swap the columns, if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_d(u, tempMat, k, k, norm[perm[k]], perm, tol);
    // find the z associated to u & tempMat
    findZ_mat_d(z, tempMat, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of tempMat using every entry in u & every entry in z
    apphh_Z_mat_d(tempMat, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);
    // find the z associated to u & tempVec
    findZ_vec_d(tempComp, tempVec, u, k);
    // apply tempComp to rows k to m tempVec using every entry in u
    apphh_Z_vec_d(tempVec, u, tempComp, k, m, 0, u->size);

    // update the norms and find the one that is the maximum
    pivot = k+1;
    max = 0;
    for (j = k+1; j < n; j++)
    {
      norm[perm[j]] = 0;
      for (i = k+1; i < m; i++)
        norm[perm[j]] += tempMat->entry[i][perm[j]].r * tempMat->entry[i][perm[j]].r + tempMat->entry[i][perm[j]].i * tempMat->entry[i][perm[j]].i;
      norm[perm[j]] = sqrt(norm[perm[j]]);
      if (max < norm[perm[j]])
      {
        max = norm[perm[j]];
        pivot = j;
      }
    }
  }
  // so the rank of A is k

  // now do back substitutions to find the least squares solution for x
  for (i = n-1; i >= 0; i--)
  {
    if (i < k)
    { // calculate x->coord[i]
      set_zero_d(tempComp);
      for (j = i + 1; j < k; j++)
      {
        sum_mul_d(tempComp, &tempMat->entry[i][perm[j]], &x->coord[perm[j]]);
      }
      sub_d(tempComp, &tempVec->coord[i], tempComp);
      div_d(&x->coord[perm[i]], tempComp, &tempMat->entry[i][perm[i]]);
    }
    else
    {
      set_zero_d(&x->coord[perm[i]]);
    }
  }

  // norm_A_inv = norm_y
  *norm_A_inv = norm_y;
  // calculate cond_num
  if (norm_y < 0)
  {
    *cond_num = -1;
  }
  else
  {
    *cond_num = (*norm_A_inv) * (*norm_A);
  }

  // check answer
  for (i = 0; i < n; i++)
    if (isnan(x->coord[i].r) || isnan(x->coord[i].i))
    { // clear
      clear_vec_d(tempVec); clear_vec_d(u); clear_vec_d(z); clear_vec_d(y);
      clear_mat_d(tempMat);
      free(norm);
      free(perm);

      return MS_NOSOLUTION;
    }

  clear_vec_d(tempVec); clear_vec_d(u); clear_vec_d(z); clear_vec_d(y);
  clear_mat_d(tempMat);
  free(norm);
  free(perm);

  return (n - k);
}

int matrixSolve_Hessenberg_Least_Squares_d(vec_d x, mat_d A, vec_d b, double tol, double largeChange)
/***************************************************************\
* USAGE: Least squares solving using QR (Givens rotations) on an*
* upper Hessenberg matrix                                       *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is least squares solution to Ax = b*
* NOTES: does least squares solution to Ax = b                  *
\***************************************************************/
{
  int i, j, k, m = A->rows, n = A->cols;
  double tempD, prevD;
  comp_d c, s, c1, s1, tempComp;
  vec_d tempVec;
  mat_d tempMat;

  // error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_Hessenberg_Least_Squares_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_Hessenberg_Least_Squares_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  // setup the size of x
  change_size_vec_d(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    return 0;
  }
  // so we can assume that m >= n >= 1

  // copy A to tempMat & b to tempVec
  init_vec_d(tempVec, m);
  init_mat_d(tempMat, m, n);
  tempMat->rows = tempVec->size = m;
  tempMat->cols = n;
  // copy first column of A & b
  for (i = 0; i < m; i++)
  { // copy b to tempVec
    set_d(&tempVec->coord[i], &b->coord[i]);
  
    // copy A to tempMat
    for (j = 0; j < n; j++)
    {
      set_d(&tempMat->entry[i][j], &A->entry[i][j]);
    }
  }

  // main algorithm - do Givens rotation for entry (k+1,k)
  for (k = 0; k < n; k++)
  { // generate the Givens rotation and apply to the kth column
    gengr_d(c, s, &tempMat->entry[k][k], tempMat, k+1, k, tol);
    set_zero_d(&tempMat->entry[k+1][k]);

    // find conj(c) & conj(s)
    conjugate_d(c1, c);
    conjugate_d(s1, s);

    // apply the Givens rotation to the kth & (k+1)st rows and (k+1:m) columns of tempMat
    for (j = k + 1; j < n; j++)
    { // apply to jth column
      mul_d(tempComp, s1, &tempMat->entry[k][j]); // conj(s) * a
      mul_d(&tempMat->entry[k][j], c, &tempMat->entry[k][j]); // c * a
      sum_mul_d(&tempMat->entry[k][j], s, &tempMat->entry[k+1][j]); // c*a + s*b
      mul_d(&tempMat->entry[k+1][j], c1, &tempMat->entry[k+1][j]); // conj(c) * b
      sub_d(&tempMat->entry[k+1][j], &tempMat->entry[k+1][j], tempComp); // conj(c)*b - conj(s)*a
    }

    // apply the Givens rotation to the kth & (k+1)st entries of tempVec
    mul_d(tempComp, s1, &tempVec->coord[k]); // conj(s) * a
    mul_d(&tempVec->coord[k], c, &tempVec->coord[k]); // c * a
    sum_mul_d(&tempVec->coord[k], s, &tempVec->coord[k+1]); // c*a + s*b
    mul_d(&tempVec->coord[k+1], c1, &tempVec->coord[k+1]); // conj(c) * b
    sub_d(&tempVec->coord[k+1], &tempVec->coord[k+1], tempComp); // conj(c)*b - conj(s)*a
  }

  // find the rank of A by looking at the diagonal entries
  prevD = 1; // initialize to something
  for (k = 0; k < n; k++)
  { // check to make sure pivot element is not too small
    tempD = d_abs_d(&tempMat->entry[k][k]);
    if (checkGood_d(tempD, prevD, tol, largeChange))
    { // this one is good
      prevD = tempD;
    }
    else
    { // this one is bad
      break;
    }
  }
  // so the rank of A is k

  // now do back substitutions to find the least squares solution for x
  for (i = n - 1; i >= 0; i--)
  {
    if (i < k)
    { // calculate x->coord[i]
      set_zero_d(tempComp);
      for (j = i + 1; j < k; j++)
      {
        sum_mul_d(tempComp, &tempMat->entry[i][j], &x->coord[j]);
      }
      sub_d(tempComp, &tempVec->coord[i], tempComp);
      div_d(&x->coord[i], tempComp, &tempMat->entry[i][i]);
    }
    else
    {
      set_zero_d(&x->coord[i]);
    }
  }

  // check answer
  for (i = 0; i < n; i++)
    if (isnan(x->coord[i].r) || isnan(x->coord[i].i))
    { // clear
      clear_vec_d(tempVec);
      clear_mat_d(tempMat);

      return MS_NOSOLUTION;
    }

  clear_vec_d(tempVec);
  clear_mat_d(tempMat);

  return (n - k);
}

/********** MULTI PRECISION ******************/

int matrixSolve_mp(vec_mp x, mat_mp A, vec_mp b)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, m = A->rows, n = A->cols, num_digits = - (int) floor(mpf_get_prec(A->entry[0][0].r) * log10(2.0) - 2.5); // Computes the number of digits being used(-1.5 for round off errors)
  double c, nA,nA_inv, VerySmall = num_digits < -307 ? 1e-307 : pow(10, num_digits), VeryLargeChange = num_digits < -310 ? 1e307 : pow(10, -num_digits - 3);

  // check that we have atleast as many rows as columns and that the size of b is correct
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // if square, try to use LU solving first
  if (m == n)
  {
    retVal = matrixSolve_LU_mp(x, A, b, VerySmall, VeryLargeChange);

    if (retVal) // if LU solving failed, try QR solving to see if really rank deficient
    {
      retVal = matrixSolve_Least_Squares_mp(x, A, b, VerySmall, VeryLargeChange, &c, &nA, &nA_inv); // retVal == corank of the matrix
      if (retVal != 0)
        retVal = MS_NOSOLUTION;
    }
  }
  else // m > n so we need to use least squares solution
  {
    retVal = matrixSolve_Least_Squares_mp(x, A, b, VerySmall, VeryLargeChange, &c, &nA, &nA_inv); // retVal == corank of the matrix
  }

  return retVal;
}

int matrixSolve_cond_num_norms_mp(vec_mp x, mat_mp A, vec_mp b, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES: does matrixSolve & condition number together!!         *
\***************************************************************/
{
  int retVal, m = A->rows, n = A->cols, num_digits = - (int) floor(mpf_get_prec(A->entry[0][0].r) * log10(2.0) - 1.5); // Computes the number of digits being used(-1.5 for round off errors)
  double VerySmall = num_digits < -307 ? 1e-307 : pow(10, num_digits), VeryLargeChange = num_digits < -310 ? 1e307 : pow(10, -num_digits - 3);

  // check that we have atleast as many rows as columns and that the size of b is correct
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_cond_num_norms_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_cond_num_norms_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // if square, try to use LU solving first
  if (m == n)
  {
    retVal = matrixSolve_LU_cond_num_norms_mp(x, A, b, VerySmall, VeryLargeChange, cond_num, norm_A, norm_A_inv);
    if (retVal) // if LU solving failed, try QR solving to see if really rank deficient
    {
      retVal = matrixSolve_svd_Least_Squares_mp(x, A, b, VerySmall, VeryLargeChange, cond_num, norm_A, norm_A_inv); // retVal == corank of the matrix
      if (retVal != 0)
        retVal = MS_NOSOLUTION;
    }
  }
  else // m > n so we need to use least squares solution
  {
    retVal = matrixSolve_Least_Squares_mp(x, A, b, VerySmall, VeryLargeChange, cond_num, norm_A, norm_A_inv); // retVal == corank of the matrix
  }

  return retVal;
}

int matrixSolve_LU_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/*
 (JDH - 06/19/06) We shall assume that A is n by n and b is n by 1.
 We solve for x where A*x = b.  This is done using an implementation of
 Gauss Elimination with Scaled Partial Pivoting as described in
 Numerical Methods, Software, and Analysis, 2nd Ed, John R. Rice,
 algorithms 6.2.1 and 6.2.4.
*/
{
  int i, j, k, l, pivot, *rownum = NULL, n = A->rows;
  double max, tempDouble, *scale = NULL, prevNorm = 1e300;
  comp_mp tempComp1, tempComp2;
  mat_mp tempA;
  vec_mp tempb;

  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in matrixSolve_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup rownum, scale
  rownum = (int *)bmalloc(n * sizeof(int));
  scale = (double *)bmalloc(n * sizeof(double));

  // setup x
  change_size_vec_mp(x, n);
  x->size = n;

  // initialize
  init_mp(tempComp1); init_mp(tempComp2);
  // initialize and copy A to tempA & b to tempb
  init_mat_mp(tempA, n, n);
  init_vec_mp(tempb, n);
  tempA->rows = tempA->cols = tempb->size = n;

  // find scale factors for each row and assign row numbers
  for (i = 0; i < n; i++)
  { // copy tempb
    set_mp(&tempb->coord[i], &b->coord[i]);
    // assign row numbers
    rownum[i] = i;
    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    max = 0;
    for (j = 0; j < n; j++)
    {
      set_mp(&tempA->entry[i][j], &A->entry[i][j]);
      tempDouble = d_oneNorm_mp(&tempA->entry[i][j]);
      if (max < tempDouble)
        max = tempDouble;
    }

    if (max < tol)
    { // clear mp
      clear_mp(tempComp1); clear_mp(tempComp2);
      clear_mat_mp(tempA); clear_vec_mp(tempb);
      free(rownum); free(scale);

      return MS_NOSOLUTION;
    }

    scale[i] = 1 / max; // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find scaled maximum element in column k
    max = d_oneNorm_mp(&tempA->entry[rownum[k]][k]) * scale[rownum[k]];
    for (i = k + 1; i < n; i++)
    {
      tempDouble = d_oneNorm_mp(&tempA->entry[rownum[i]][k]) * scale[rownum[i]];
      if (max < tempDouble)
      {
        pivot = i;
        max = tempDouble;
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    if (max < tol || (k > 0 && prevNorm > max * largeChange))
    { // clear mp
      clear_mp(tempComp1); clear_mp(tempComp2);
      clear_mat_mp(tempA); clear_vec_mp(tempb);
      free(rownum); free(scale);

      return MS_NOSOLUTION;
    }
    else
      prevNorm = max;

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    // update elements in "bottom" block of A
    pivot = rownum[k];
    recip_mp(tempComp2, &tempA->entry[pivot][k]);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];
      // compute muliplier
      mul_mp(tempComp1, &tempA->entry[l][k], tempComp2);
      // update tempb (RHS)
      sub_mul_mp(&tempb->coord[l], tempComp1, &tempb->coord[pivot]);
      // update elements of row i
      for (j = k + 1; j < n; j++)
        sub_mul_mp(&tempA->entry[l][j], tempComp1, &tempA->entry[pivot][j]);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  max = d_oneNorm_mp(&tempA->entry[pivot][l]) * scale[pivot];
  if (max < tol || prevNorm > max * largeChange)
  { // clear mp
    clear_mp(tempComp1); clear_mp(tempComp2);
    clear_mat_mp(tempA); clear_vec_mp(tempb);
    free(rownum); free(scale);

    return MS_NOSOLUTION;
  }

  // calculate x->coord[n-1]
  div_mp(&x->coord[l], &tempb->coord[pivot], &tempA->entry[pivot][l]);
  for (i = n - 2; i >= 0; i--)
  { // calculate x->coord[i]
    k = rownum[i];
    set_mp(tempComp1, &tempb->coord[k]);
    for (j = i + 1; j < n; j++)
      sub_mul_mp(tempComp1, &tempA->entry[k][j], &x->coord[j]);
    div_mp(&x->coord[i], tempComp1, &tempA->entry[k][i]);
  }

  // clear mp
  clear_mp(tempComp1); clear_mp(tempComp2);
  clear_mat_mp(tempA); clear_vec_mp(tempb);
  free(rownum); free(scale);

  return 0;
}

int matrixSolve2_LU_mp(vec_mp x1, vec_mp x2, mat_mp A, vec_mp b1, vec_mp b2, double tol, double largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, l, pivot, *rownum = NULL, n = A->rows;
  double max, tempDouble, prevNorm = 1e300, *scale = NULL;

  comp_mp tempComp1, tempComp2, tempComp3;
  mat_mp tempA;
  vec_mp tempb1, tempb2;

  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in matrixSolve2_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b1->size || A->rows != b2->size)
  {
    printf("ERROR: The sizes in matrixSolve2_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup rownum, scale
  rownum = (int *)bmalloc(n * sizeof(int));
  scale = (double *)bmalloc(n * sizeof(double));

  // setup x1 & x2
  change_size_vec_mp(x1, n);
  change_size_vec_mp(x2, n);
  x1->size = x2->size = n;

  // initialize
  init_mp(tempComp1); init_mp(tempComp2); init_mp(tempComp3);

  // initialize and copy A to tempA & b1 to tempb1 & b2 to tempb2
  init_mat_mp(tempA, n, n);
  init_vec_mp(tempb1, n);
  init_vec_mp(tempb2, n);
  tempA->rows = tempA->cols = tempb1->size = tempb2->size = n;

  // find scale factors for each row and assign row numbers
  for (i = 0; i < n; i++)
  { // setup tempb1 & tempb2
    set_mp(&tempb1->coord[i], &b1->coord[i]);
    set_mp(&tempb2->coord[i], &b2->coord[i]);
    // assign row numbers
    rownum[i] = i;
    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    max = 0;
    for (j = 0; j < n; j++)
    {
      set_mp(&tempA->entry[i][j], &A->entry[i][j]);
      tempDouble = d_oneNorm_mp(&tempA->entry[i][j]);
      if (max < tempDouble)
        max = tempDouble;
    }
    if (max < tol)
    { // clear mp
      clear_mp(tempComp1); clear_mp(tempComp2); clear_mp(tempComp3);
      clear_mat_mp(tempA); clear_vec_mp(tempb1); clear_vec_mp(tempb2);
      free(rownum); free(scale);

      return MS_NOSOLUTION;
    }

    scale[i] = 1 / max; // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find sacled maximum element in column k
    max = d_oneNorm_mp(&tempA->entry[rownum[k]][k]) * scale[rownum[k]];
    for (i = k + 1; i < n; i++)
    {
      tempDouble = d_oneNorm_mp(&tempA->entry[rownum[i]][k]) * scale[rownum[i]];
      if (max < tempDouble)
      {
        pivot = i;
        max = tempDouble;
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    if (max < tol || (k > 0 && prevNorm > max * largeChange))
    { // clear mp
      clear_mp(tempComp1); clear_mp(tempComp2); clear_mp(tempComp3);
      clear_mat_mp(tempA); clear_vec_mp(tempb1); clear_vec_mp(tempb2);
      free(rownum); free(scale);

      return MS_NOSOLUTION;
    }
    else
      prevNorm = max;

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
    }

    // update elements in "bottom" block of A
    pivot = rownum[k];
    recip_mp(tempComp3, &tempA->entry[pivot][k]);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];
      // compute muliplier
      mul_mp(tempComp1, &tempA->entry[l][k], tempComp3);
      // update the RHS
      sub_mul_mp(&tempb1->coord[l], tempComp1, &tempb1->coord[pivot]);
      sub_mul_mp(&tempb2->coord[l], tempComp1, &tempb2->coord[pivot]);
      // update elements of row i
      for (j = k + 1; j < n; j++)
        sub_mul_mp(&tempA->entry[l][j], tempComp1, &tempA->entry[pivot][j]);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  max = d_oneNorm_mp(&tempA->entry[pivot][l]) * scale[pivot];
  if (max < tol || prevNorm > max * largeChange)
  { // clear mp
    clear_mp(tempComp1); clear_mp(tempComp2); clear_mp(tempComp3);
    clear_mat_mp(tempA); clear_vec_mp(tempb1); clear_vec_mp(tempb2);
    free(rownum); free(scale);

    return MS_NOSOLUTION;
  }

  // calculate x1->coord[n-1], x2->coord[n-1]
  recip_mp(tempComp3, &tempA->entry[pivot][l]);
  mul_mp(&x1->coord[l], &tempb1->coord[pivot], tempComp3);
  mul_mp(&x2->coord[l], &tempb2->coord[pivot], tempComp3);
  // compute the rest of the solution vectors
  for (i = n - 2; i >= 0; i--)
  { // calculate x1->coord[i], x2->coord[i]
    k = rownum[i];
    set_mp(tempComp1, &tempb1->coord[k]);
    set_mp(tempComp2, &tempb2->coord[k]);
    for (j = i + 1; j < n; j++)
    {
      sub_mul_mp(tempComp1, &tempA->entry[k][j], &x1->coord[j]);
      sub_mul_mp(tempComp2, &tempA->entry[k][j], &x2->coord[j]);
    }
    div_mp(&x1->coord[i], tempComp1, &tempA->entry[k][i]);
    div_mp(&x2->coord[i], tempComp2, &tempA->entry[k][i]);
  }

  // clear mp
  clear_mp(tempComp1); clear_mp(tempComp2); clear_mp(tempComp3);
  clear_mat_mp(tempA); clear_vec_mp(tempb1); clear_vec_mp(tempb2);
  free(rownum); free(scale);

  return 0;
}

int matrixSolve_LU_cond_num_norms_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does matrix solve & condition number together          *
\***************************************************************/
{
  int  i, numCols = A->cols, retVal = 0;
  double maxNorm, tempNorm;

  if ((A->rows == 1) && (A->cols == 1))
  {
    *norm_A = d_abs_mp(&A->entry[0][0]);
    *norm_A_inv = 1 / (*norm_A);
    *cond_num = 1;
    retVal = matrixSolve_mp(x, A, b);
  }
  else
  {
    vec_mp tempVec;
    init_vec_mp(tempVec, numCols);
    tempVec->size = numCols;

    // obtain norm of A
    *norm_A = infNormMat_mp(A);

    // Estimate ||M^{-1}|| by :
    //   ||M^{-1}|| ~ 2||M^{-1}v|| , for a random unit vector v.

    maxNorm = 0;
    for (i = 0; i < numCols; i++)
    {
      // obtain complex random number and its modulus
      tempNorm = get_comp_rand_mp(&tempVec->coord[i]);

      // determine if it is the largest modulus so far
      if (maxNorm < tempNorm)
        maxNorm = tempNorm;
    }
    // update maxNorm to be 1 / max modulus = 1 / infNormVec_mp(tempVec)
    maxNorm = 1.0 / maxNorm;

    // solve for x and 'tempVec' where A*x = b and A * 'tempVec3' = tempVec3
    retVal = matrixSolve2_LU_mp(x, tempVec, A, b, tempVec, tol, largeChange);

    if (retVal)
    { // matrixSolve2_mp failed
      *norm_A_inv = -1.0;
      *cond_num = -1.0;
    }
    else
    { // correct result by multiplying by maxNorm - makes the input tempVec a unit vector
      *norm_A_inv = 2.0 * infNormVec_mp(tempVec) * maxNorm;
      *cond_num = (*norm_A) * (*norm_A_inv);
    }
    
    clear_vec_mp(tempVec);
  }

  return retVal;
}

int matrixSolve_svd_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using SVD!!                                 *
\***************************************************************/
// INPUT AS DOUBLES - convert to MP and run mS
{
  int retVal;
  mpf_t tol_mp, largeChange_mp;

  mpf_init_set_d(tol_mp, tol);
  mpf_init_set_d(largeChange_mp, largeChange);

  retVal = matrixSolve_svd_Least_Squares_mp2(x, A, b, tol_mp, largeChange_mp, cond_num, norm_A, norm_A_inv);

  mpf_clear(tol_mp);
  mpf_clear(largeChange_mp);

  return retVal;
}

int matrixSolve_svd_Least_Squares_mp2(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using SVD!!                                 *
\***************************************************************/
{
  int i, j, retVal, rank, its = 50, m = A->rows, n = A->cols, prec = mpf_get_prec(A->entry[0][0].r);
  int num_digits = - (int) floor(prec * log10(2.0) - 0.5); // Computes the number of digits being used

  char *str = NULL;
  size_t size;

  mpf_t tol_prec, tol_sign, tempMPF;
  comp_mp tempComp;
  vec_mp tempVec;
  mat_mp U, E, V;

  // do error checking
  if (mpf_cmp_ui(tol, 0) < 0)
  {
    printf("ERROR: The tolerance for matrixSolve_svd_Least_Squares_mp needs to be non-negative!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(largeChange, 0) < 0)
  {
    printf("ERROR: The change tolerance for matrixSolve_svd_Least_Squares_mp needs to be non-negative!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_svd_Least_Squares_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_svd_Least_Squares_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup the size of x
  change_size_vec_mp(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    *cond_num = *norm_A = *norm_A_inv = 1;
    return 0;
  }
  // so we can assume that m >= n >= 1

  // initialize MP
  mpf_init2(tol_prec, prec); mpf_init2(tol_sign, prec); mpf_init2(tempMPF, prec); 
  init_mp2(tempComp, prec);
  init_vec_mp2(tempVec, 0, prec);
  init_mat_mp2(U, 0, 0, prec); init_mat_mp2(E, 0, 0, prec); init_mat_mp2(V, 0, 0, prec);

  // set tol_prec
  size = 1 + snprintf(NULL, 0, "1e%d", num_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", num_digits);
  mpf_set_str(tol_prec, str, 10);
  // set tol_sign
  size = 1 + snprintf(NULL, 0, "1e%d", 2 * num_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", 2 * num_digits);
  mpf_set_str(tol_sign, str, 10);

  // find the SVD of A using at most its iterations with tolerances of
  // tol, tol, tol_sign, and largeChange
  retVal = svd_jacobi_mp(U, E, V, A, its, tol, tol_prec, tol_prec, tol_sign, largeChange);

  if (retVal >= 0)
  { // the svd is good with corank == retVal, so we find x = V * (E^-1)^H * U^H * b

    // find the rank 
    rank = n - retVal;

    // setup the size of tempVec
    change_size_vec_mp(tempVec, rank);

    // find (E^-1)^H * U^H * b - take into account the rank to save operations
    for (j = 0; j < rank; j++)
    {
      set_zero_mp(&tempVec->coord[j]);
      for (i = 0; i < m; i++)
      {
        conjugate_mp(tempComp, &U->entry[i][j]);
        sum_mul_mp(&tempVec->coord[j], tempComp, &b->coord[i]);
      }
      mpf_ui_div(tempMPF, 1, E->entry[j][j].r);  // E[j][j] is real and large enough to be inverted
      mpf_mul(tempVec->coord[j].r, tempVec->coord[j].r, tempMPF);
      mpf_mul(tempVec->coord[j].i, tempVec->coord[j].i, tempMPF);
    }

    // find V * tempVec - take into account the rank to save operations
    for (i = 0; i < n; i++)
    {
      set_zero_mp(&x->coord[i]);
      for (j = 0; j < rank; j++)
      {
        sum_mul_mp(&x->coord[i], &V->entry[i][j], &tempVec->coord[j]);
      }
    }

    // find the norm of A - the largest singular value
    *norm_A = mpf_get_d(E->entry[0][0].r);
    // find the norm of A^-1 - the inverse of the smallest singular value
    if (mpf_cmp_ui(E->entry[n-1][n-1].r, 0) > 0)
    {
      mpf_ui_div(tempMPF, 1, E->entry[n-1][n-1].r);
      *norm_A_inv = mpf_get_d(tempMPF);
      *cond_num = (*norm_A) * (*norm_A_inv);
    }
    else
    { // smallest singular value is 0
      *cond_num = *norm_A_inv = -1;
    }
  }
  else // give up
  {
    *norm_A = infNormMat_mp(A);
    *cond_num = *norm_A_inv = -1;
    retVal = MS_NOSOLUTION;
  }

  // free str
  free(str);

  // clear MP
  mpf_clear(tol_prec); mpf_clear(tol_sign); mpf_clear(tempMPF); 
  clear_mp(tempComp);
  clear_vec_mp(tempVec);
  clear_mat_mp(U); clear_mat_mp(E); clear_mat_mp(V);

  return retVal;
}

int matrixSolve_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, double tol, double largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using QR!!                                  *
\***************************************************************/
// INPUT AS DOUBLES - convert to MP and run mS
{
  int retVal;
  mpf_t tol_mp, largeChange_mp;

  mpf_init_set_d(tol_mp, tol);
  mpf_init_set_d(largeChange_mp, largeChange);

  retVal = matrixSolve_Least_Squares_mp2(x, A, b, tol_mp, largeChange_mp, cond_num, norm_A, norm_A_inv);

  mpf_clear(tol_mp);
  mpf_clear(largeChange_mp);

  return retVal;
}

int matrixSolve_Least_Squares_mp2(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange, double *cond_num, double *norm_A, double *norm_A_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is solution to Ax=b                *
* NOTES: does least squares solution to Ax = b & calculates the *
* condition number, using QR!!                                  *
\***************************************************************/
{
  int i, j, k, pivot, *perm = NULL, m = A->rows, n = A->cols, prec = mpf_get_prec(A->entry[0][0].r);
  comp_mp tempComp;
  mpf_t tempMPF, max, prevNorm, norm_y, *norm = NULL;
  vec_mp tempVec, u, z, y;
  mat_mp tempMat;
  
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_Least_Squares_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_Least_Squares_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup the size of x
  change_size_vec_mp(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    *cond_num = *norm_A = *norm_A_inv = 1;
    return 0;
  }
  // so we can assume that m >= n >= 1

  // initialize MP
  init_mp2(tempComp, prec);
  mpf_init2(tempMPF, prec); mpf_init2(max, prec); mpf_init2(prevNorm, prec); mpf_init2(norm_y, prec);
  init_vec_mp(u, 0);
  init_vec_mp(z, 0);
  init_vec_mp(y, 0);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (i = 0; i < n; i++)
    mpf_init2(norm[i], prec);

  // copy A to tempMat & b to tempVec, find norms for each col and find the one that is the maximum
  init_vec_mp(tempVec, m);
  init_mat_mp(tempMat, m, n);
  tempMat->rows = tempVec->size = m;
  tempMat->cols = n;
  pivot = 0;
  mpf_set_ui(max, 0);
  for (j = 0; j < m; j++)
  { // copy b to tempVec
    set_mp(&tempVec->coord[j], &b->coord[j]);

    if (j < n)
    { // initialize perm
      perm[j] = j;

      // find the twoNorm for each column
      mpf_set_ui(norm[j], 0);
      for (i = 0; i < m; i++)
      { // copy A to tempMat
        set_mp(&tempMat->entry[i][j], &A->entry[i][j]);
        // find the size to this entry
        mpf_mul(tempComp->r, A->entry[i][j].r, A->entry[i][j].r);
        mpf_mul(tempComp->i, A->entry[i][j].i, A->entry[i][j].i);
        mpf_add(tempComp->r, tempComp->r, tempComp->i);
        mpf_add(norm[j], norm[j], tempComp->r);
      }
      mpf_sqrt(norm[j], norm[j]);
      // check to see if it is the maximum twoNorm thus far
      if (mpf_cmp(max, norm[j]) < 0)
      {
        pivot = j;
        mpf_set(max, norm[j]);
      }
    }
  }
  mpf_set(prevNorm, max);
  // take the norm_A as the sqrt(n) * maximum of the twoNorms, which a good enough estimate of the norm of A
  *norm_A = sqrt(n) * mpf_get_d(max);

  // initialize the I.C.E.
  y->size = 0;
  mpf_set_ui(norm_y, 0);

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  { 
    if (mpf_sgn(norm[perm[pivot]]) > 0)
    { // setup tempComp for I.C.E.
      sign_mp2(tempComp, &tempMat->entry[k][perm[pivot]], tol);
      mpf_neg(tempMPF, norm[perm[pivot]]);
      mpf_mul(tempComp->r, tempComp->r, tempMPF);
      mpf_mul(tempComp->i, tempComp->i, tempMPF);
      // find the I.C.E.
      cond_est_mp(y, y, norm_y, tempMat, 0, k, pivot, tempComp, perm, tol);
    }
    else
    { // next pivot is 0 and so there is no use trying to find the next singular value since it is 0!! - the next if statement will stop the loop
      mpf_set_si(norm_y, -1);
    }

    // check to make sure pivot element is not too small
    mpf_mul(tempMPF, norm[perm[pivot]], largeChange);
    if (mpf_cmp(norm[perm[pivot]], tol) < 0 || (k > 0 && mpf_cmp(prevNorm, tempMPF) > 0))
      break;
    else
      mpf_set(prevNorm, norm[perm[pivot]]);

    // swap the columns, if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_mp(u, tempMat, k, k, norm[perm[k]], perm, tol);
    // find the z associated to u & tempMat
    findZ_mat_mp(z, tempMat, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of tempMat using every entry in u & every entry in z
    apphh_Z_mat_mp(tempMat, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);
    // find the z associated to u & tempVec
    findZ_vec_mp(tempComp, tempVec, u, k);
    // apply tempComp to rows k to m of tempVec using every entry in u
    apphh_Z_vec_mp(tempVec, u, tempComp, k, m, 0, u->size);

    // update the norms and find the one that is the maximum
    pivot = k+1;
    mpf_set_ui(max, 0);
    for (j = k+1; j < n; j++)
    {
      mpf_set_ui(norm[perm[j]], 0);
      for (i = k+1; i < m; i++)
      {
        mpf_mul(tempComp->r, tempMat->entry[i][perm[j]].r, tempMat->entry[i][perm[j]].r);
        mpf_mul(tempComp->i, tempMat->entry[i][perm[j]].i, tempMat->entry[i][perm[j]].i);
        mpf_add(tempComp->r, tempComp->r, tempComp->i);
        mpf_add(norm[perm[j]], norm[perm[j]], tempComp->r);
      }
      mpf_sqrt(norm[perm[j]], norm[perm[j]]);
      if (mpf_cmp(max, norm[perm[j]]) < 0)
      {
        mpf_set(max, norm[perm[j]]);
        pivot = j;
      }
    }
  }
  // so the rank of A is k

  // now do back substitutions to find the least squares solution for x
  for (i = n-1; i >= 0; i--)
  {
    if (i < k)
    { // calculate x->coord[i]
      set_zero_mp(tempComp);
      for (j = i + 1; j < k; j++)
      {
        sum_mul_mp(tempComp, &tempMat->entry[i][perm[j]], &x->coord[perm[j]]);
      }
      sub_mp(tempComp, &tempVec->coord[i], tempComp);
      div_mp(&x->coord[perm[i]], tempComp, &tempMat->entry[i][perm[i]]);
    }
    else
    {
      set_zero_mp(&x->coord[perm[i]]);
    }
  }

  // norm_A_inv = norm_y
  *norm_A_inv = mpf_get_d(norm_y);
  // calculate cond_num
  if (*norm_A_inv < 0)
  {
    *cond_num = -1;
  }
  else
  {
    *cond_num = (*norm_A_inv) * (*norm_A);
  }

  // clear MP
  clear_mp(tempComp);
  mpf_clear(tempMPF); mpf_clear(max); mpf_clear(prevNorm); mpf_clear(norm_y);
  clear_vec_mp(tempVec); clear_vec_mp(u); clear_vec_mp(z); clear_vec_mp(y);
  clear_mat_mp(tempMat);
  for (i = 0; i < n; i++)
    mpf_clear(norm[i]);
  free(norm);

  return (n - k);
}

int matrixSolve_Hessenberg_Least_Squares_mp(vec_mp x, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange)
/***************************************************************\
* USAGE: Least squares solving using QR (Givens rotations) on an*
* upper Hessenberg matrix                                       *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> x is least squares solution to Ax = b*
* NOTES: does least squares solution to Ax = b                  *
\***************************************************************/
{
  int i, j, k, m = A->rows, n = A->cols;
  mpf_t tempMPF, prevMPF;
  comp_mp c, s, c1, s1, tempComp;
  vec_mp tempVec;
  mat_mp tempMat;

  // error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in matrixSolve_Hessenberg_Least_Squares_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != b->size)
  {
    printf("ERROR: The sizes in matrixSolve_Hessenberg_Least_Squares_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  // setup the size of x
  change_size_vec_mp(x, n);
  x->size = n;

  if (m == 0 || n == 0) // if there are no rows/cols, then just return
  {
    return 0;
  }
  // so we can assume that m >= n >= 1

  // initialize 
  mpf_init(tempMPF); mpf_init(prevMPF);
  init_mp(c); init_mp(s); init_mp(c1); init_mp(s1); init_mp(tempComp);

  // copy A to tempMat & b to tempVec
  init_vec_mp(tempVec, m);
  init_mat_mp(tempMat, m, n);
  tempMat->rows = tempVec->size = m;
  tempMat->cols = n;
  for (i = 0; i < m; i++)
  { // copy b to tempVec
    set_mp(&tempVec->coord[i], &b->coord[i]);

    // copy A to tempMat
    for (j = 0; j < n; j++)
    { // copy A to tempMat
      set_mp(&tempMat->entry[i][j], &A->entry[i][j]);
    }
  }

  // main algorithm - do Givens rotation for entry (k+1,k)
  for (k = 0; k < n; k++)
  { // generate the Givens rotation and apply to the kth column
    gengr_mp(c, s, &tempMat->entry[k][k], tempMat, k+1, k, tol);
    set_zero_mp(&tempMat->entry[k+1][k]);

    // find conj(c) & conj(s)
    conjugate_mp(c1, c);
    conjugate_mp(s1, s);

    // apply the Givens rotation to the kth & (k+1)st rows and (k+1:m) columns of tempMat
    for (j = k+1; j < n; j++)
    { // apply to jth column
      mul_mp(tempComp, s1, &tempMat->entry[k][j]); // conj(s) * a
      mul_mp(&tempMat->entry[k][j], c, &tempMat->entry[k][j]); // c * a
      sum_mul_mp(&tempMat->entry[k][j], s, &tempMat->entry[k+1][j]); // c*a + s*b
      mul_mp(&tempMat->entry[k+1][j], c1, &tempMat->entry[k+1][j]); // conj(c) * b
      sub_mp(&tempMat->entry[k+1][j], &tempMat->entry[k+1][j], tempComp); // conj(c)*b - conj(s)*a
    }

    // apply the Givens rotation to the kth & (k+1)st entries of tempVec
    mul_mp(tempComp, s1, &tempVec->coord[k]); // conj(s) * a
    mul_mp(&tempVec->coord[k], c, &tempVec->coord[k]); // c * a
    sum_mul_mp(&tempVec->coord[k], s, &tempVec->coord[k+1]); // c*a + s*b
    mul_mp(&tempVec->coord[k+1], c1, &tempVec->coord[k+1]); // conj(c) * b
    sub_mp(&tempVec->coord[k+1], &tempVec->coord[k+1], tempComp); // conj(c)*b - conj(s)*a
  }

  // find the rank of A by looking at the diagonal entries
  mpf_set_ui(prevMPF, 1); // initialize to something
  for (k = 0; k < n; k++)
  { // check to make sure pivot element is not too small
    mpf_abs_mp(tempMPF, &tempMat->entry[k][k]);
    if (checkGood_mp(tempMPF, prevMPF, tol, largeChange))
    { // this one is good
      mpf_set(prevMPF, tempMPF);
    }
    else
    { // this one is bad
      break;
    }
  }
  // so the rank of A is k

  // now do back substitutions to find the least squares solution for x
  for (i = n-1; i >= 0; i--)
  {
    if (i < k)
    { // calculate x->coord[i]
      set_zero_mp(tempComp);
      for (j = i + 1; j < k; j++)
      {
        sum_mul_mp(tempComp, &tempMat->entry[i][j], &x->coord[j]);
      }
      sub_mp(tempComp, &tempVec->coord[i], tempComp);
      div_mp(&x->coord[i], tempComp, &tempMat->entry[i][i]);
    }
    else
    {
      set_zero_mp(&x->coord[i]);
    }
  }

  // clear memory
  mpf_clear(tempMPF); mpf_clear(prevMPF);
  clear_mp(c); clear_mp(s); clear_mp(c1); clear_mp(s1); clear_mp(tempComp);
  clear_vec_mp(tempVec); clear_mat_mp(tempMat);

  return (n - k);
}

////// EXTRA FUNCTIONS //////////////

int LU_matrixSolve_d(vec_d x, mat_d LU, int **rwnm, int *sign, mat_d A, vec_d b, double tol, double largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: LU - matrix containing the LU-decompositon of A*
* Also, rownum which "defines" the row-swap matrix assoc w/LU   *
*  and sign defines the sign of this permutation                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL; // rownum is simply a pointer to *rwnm (to avoid dereferencing)
  double max, tempD, tempD2, prevNorm = 1e300, *scale = NULL;
  vec_d tempVec;
  comp_d tempComp;

  // error checking
  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in LU_matrixSolve_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b->size)
  {
    printf("ERROR: The sizes in LU_matrixSolve_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // setup rwnm (rownum), scale, & sign
  rownum = *rwnm = (int *)brealloc(*rwnm, n * sizeof(int));
  scale = (double *)bmalloc(n * sizeof(double));
  *sign = 1;

  // setup x
  change_size_vec_d(x, n);
  x->size = n;

  // copy A to LU & b to tempVec
  change_size_mat_d(LU, n, n);
  init_vec_d(tempVec, n);
  LU->rows = LU->cols = tempVec->size = n;

  // find scale factors for each row and assign row numbers
  prevNorm = 0;
  for (i = 0; i < n; i++)
  { // setup tempVec
    set_d(&tempVec->coord[i], &b->coord[i]);
    // assign row numbers
    rownum[i] = i;
    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    set_d(&LU->entry[i][0], &A->entry[i][0]);
    max = d_oneNorm_d(&LU->entry[i][0]);
    for (j = 1; j < n; j++)
    {
      set_d(&LU->entry[i][j], &A->entry[i][j]);
      tempD = d_oneNorm_d(&LU->entry[i][j]);
      if (max < tempD)
        max = tempD;
    }
    if (max < tol)
    { // clear
      clear_vec_d(tempVec);
      free(scale);
      rownum = NULL;

      return MS_NOSOLUTION;
    }

    scale[i] = 1 / max; // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find sacled maximum element in column k
    max = d_oneNorm_d(&LU->entry[rownum[k]][k]) * scale[rownum[k]];
    for (i = k + 1; i < n; i++)
    {
      tempD = d_oneNorm_d(&LU->entry[rownum[i]][k]) * scale[rownum[i]];
      if (max < tempD)
      {
        pivot = i;
        max = tempD;
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    if (max < tol || (k > 0 && prevNorm > max * largeChange))
    { // clear
      clear_vec_d(tempVec);
      free(scale);
      rownum = NULL;

      return MS_NOSOLUTION;
    }

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    { 
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
      // update the sign of the permuation
      *sign = -*sign;
    }

    // update elements in "bottom" block of LU
    pivot = rownum[k];
    recip_d2(tempComp, &LU->entry[pivot][k], tempD);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];

      // compute muliplier
      mul_d2(&LU->entry[l][k], &LU->entry[l][k], tempComp, tempD);
      // update tempVec (RHS)
      sub_mul_d2(&tempVec->coord[l], &LU->entry[l][k], &tempVec->coord[pivot], tempD);

      // update elements of row i
      for (j = k + 1; j < n; j++)
        sub_mul_d2(&LU->entry[l][j], &LU->entry[l][k], &LU->entry[pivot][j], tempD);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  max = d_oneNorm_d(&LU->entry[pivot][l]) * scale[pivot];
  if (max < tol || prevNorm > max * largeChange)
  { // clear
    clear_vec_d(tempVec);
    free(scale);
    rownum = NULL;

    return MS_NOSOLUTION;
  }

  // calculate x->coord[n-1]
  div_d2(&x->coord[l], &tempVec->coord[pivot], &LU->entry[pivot][l], tempD, tempD2);

  for (i = n - 2; i >= 0; i--)
  { // calculate x->coord[i]
    k = rownum[i];
    set_d(tempComp, &tempVec->coord[k]);
    for (j = i + 1; j < n; j++)
      sub_mul_d2(tempComp, &LU->entry[k][j], &x->coord[j], tempD);
    div_d2(&x->coord[i], tempComp, &LU->entry[k][i], tempD, tempD2);
  }

  // check answer
  for (i = 0; i < n; i++)
    if (isnan(x->coord[i].r) || isnan(x->coord[i].i))
    { // clear
      clear_vec_d(tempVec);
      free(scale);
      rownum = NULL;

      return MS_NOSOLUTION;
    }

  // clear
  clear_vec_d(tempVec);
  free(scale);
  rownum = NULL;

  return 0;
}

int matrixSolve_from_LU_d(mat_d X, mat_d LU, int *rownum, mat_d B)
/***************************************************************\
* USAGE: used after finding the LU-decomp by LU_matrixSovle_d   *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: solves (L*U)X = R*B, where R is stored in rownum       *
* and L and U are stored in LU                                  *
\***************************************************************/
{
  int i, j, k, l, m, pivot, n = LU->rows, p = B->cols;
  mat_d Y;
  comp_d tempComp;
  double tempD, tempD2;
  
  if (LU->rows != LU->cols)
  { 
    printf("ERROR: The matrix is not square in matrixSolve_from_LU_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (LU->rows != B->rows)
  { 
    printf("ERROR: The sizes in matrixSolve_from_LU_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  
  // setup X & Y - same size as B
  change_size_mat_d(X, n, p);
  init_mat_d(Y, n, p);
  X->rows = Y->rows = n;
  X->cols = Y->cols = p;
  
  // first compute LY = RB
  for (k = 0; k < p; k++)
  { // calculate the kth column of Y
    set_d(&Y->entry[0][k], &B->entry[rownum[0]][k]);

    for (i = 1; i < n; i++)
    { // calculate Y[i][k]
      l = rownum[i];
      neg_d(&Y->entry[i][k], &B->entry[l][k]);
      for (j = 0; j < i; j++)
        sum_mul_d2(&Y->entry[i][k], &LU->entry[l][j], &Y->entry[j][k], tempD);
      neg_d(&Y->entry[i][k], &Y->entry[i][k]);
    }
  }
  
  // then compute UX = Y
  l = n - 1;
  pivot = rownum[l]; 
  for (k = 0; k < p; k++)
  { // calculate the kth column of X
    div_d2(&X->entry[l][k], &Y->entry[l][k], &LU->entry[pivot][l], tempD, tempD2);
    for (i = n - 2; i >= 0; i--)
    {
      set_d(tempComp, &Y->entry[i][k]);
      m = rownum[i];
      for (j = i + 1; j < n; j++)
        sub_mul_d2(tempComp, &LU->entry[m][j], &X->entry[j][k], tempD);
      div_d2(&X->entry[i][k], tempComp, &LU->entry[m][i], tempD, tempD2);
    }
  }
  
  clear_mat_d(Y);
  
  return 0;
}

int LU_matrixSolve_mp(vec_mp x, mat_mp LU, int **rwnm, int *sign, mat_mp A, vec_mp b, mpf_t tol, mpf_t largeChange)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: LU - matrix containing the LU-decompositon of A*
* Also, rownum which "defines" the row-swap matrix assoc w/LU   *
*  and sign defines the sign of this permutation                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, l, pivot, n = A->rows, *rownum = NULL; // rownum is simply a pointer to *rwnm (to avoid dereferencing) 
  mpf_t max, tempMPF, prevNorm, *scale = NULL;
  vec_mp tempVec;
  comp_mp tempComp;

  // error checking
  if (A->rows != A->cols)
  {
    printf("ERROR: The matrix is not square in LU_matrixSolve_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (A->rows != b->size)
  {
    printf("ERROR: The sizes in LU_matrixSolve_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize memory
  mpf_init(max); mpf_init(tempMPF); mpf_init(prevNorm);
  init_mp(tempComp);

  // setup rwnm (rownum), scale, and sign
  rownum = *rwnm = (int *)brealloc(*rwnm, n * sizeof(int));
  scale = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (i = 0; i < n; i++)
  {
    mpf_init(scale[i]);
    mpf_set_ui(scale[i], 0);
  }
  *sign = 1;

  // setup x
  change_size_vec_mp(x, n);
  x->size = n;

  // copy A to LU & b to tempVec
  change_size_mat_mp(LU, n, n);
  init_vec_mp(tempVec, n);
  LU->rows = LU->cols = tempVec->size = n;

  // find scale factors for each row and assign row numbers
  mpf_set_ui(prevNorm, 0);
  for (i = 0; i < n; i++)
  { // setup tempVec
    set_mp(&tempVec->coord[i], &b->coord[i]);
    // assign row numbers
    rownum[i] = i;
    // using the oneNorm on each of the entries to find the infinity norm for the row - avoid sqrt
    set_mp(&LU->entry[i][0], &A->entry[i][0]);
    mp_oneNorm_mp(max, &LU->entry[i][0]);
    for (j = 1; j < n; j++)
    {
      set_mp(&LU->entry[i][j], &A->entry[i][j]);
      mp_oneNorm_mp(tempMPF, &LU->entry[i][j]);
      if (mpf_cmp(max, tempMPF) < 0)
        mpf_set(max, tempMPF);
    }
    if (mpf_cmp(max, tol) < 0)
    { // clear
      clear_vec_mp(tempVec);
      for (i = 0; i < n; i++)
        mpf_clear(scale[i]);
      free(scale);
      rownum = NULL;
      mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
      clear_mp(tempComp);

      return MS_NOSOLUTION;
    }

    mpf_ui_div(scale[i], 1, max); // scale = 1 / scale - turn later / into *
  }

  // do the pivoting for column k
  for (k = 0; k < n - 1; k++)
  {
    pivot = k;

    // find sacled maximum element in column k
    mp_oneNorm_mp(max, &LU->entry[rownum[k]][k]);
    mpf_mul(max, max, scale[rownum[k]]);
    for (i = k + 1; i < n; i++)
    {
      mp_oneNorm_mp(tempMPF, &LU->entry[rownum[i]][k]);
      mpf_mul(tempMPF, tempMPF, scale[rownum[i]]);
      if (mpf_cmp(max, tempMPF) < 0)
      {
        pivot = i;
        mpf_set(max, tempMPF);
      }
    }

    // check to make sure the scaled pivot element is not too small or too much of a change that reveals rank deficiency
    mpf_mul(tempMPF, max, largeChange);
    if (mpf_cmp(max, tol) < 0 || (k > 0 && mpf_cmp(prevNorm, tempMPF) > 0))
    { // clear
      clear_vec_mp(tempVec);
      for (i = 0; i < n; i++)
        mpf_clear(scale[i]);
      free(scale);
      rownum = NULL;
      mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
      clear_mp(tempComp);

      return MS_NOSOLUTION;
    }

    // swap the rows, if needed
    if (rownum[pivot] != rownum[k])
    {
      i = rownum[pivot];
      rownum[pivot] = rownum[k];
      rownum[k] = i;
      // update the sign of the permuation
      *sign = -*sign;
    }

    // update elements in "bottom" block of LU
    pivot = rownum[k];
    recip_mp(tempComp, &LU->entry[pivot][k]);
    for (i = k + 1; i < n; i++)
    { // find the row number
      l = rownum[i];

      // compute muliplier
      mul_mp(&LU->entry[l][k], &LU->entry[l][k], tempComp);
      // update tempVec (RHS)
      sub_mul_mp(&tempVec->coord[l], &LU->entry[l][k], &tempVec->coord[pivot]);

      // update elements of row i
      for (j = k + 1; j < n; j++)
        sub_mul_mp(&LU->entry[l][j], &LU->entry[l][k], &LU->entry[pivot][j]);
    }
  }

  // check to make sure the last scaled pivot element is not too small or too much of a change that reveals rank deficiency
  l = n - 1;
  pivot = rownum[l];
  mp_oneNorm_mp(max, &LU->entry[pivot][l]);
  mpf_mul(max, max, scale[pivot]);
  mpf_mul(tempMPF, max, largeChange);
  if (mpf_cmp(max, tol) < 0 || mpf_cmp(prevNorm, tempMPF) > 0)
  { // clear
    clear_vec_mp(tempVec);
    for (i = 0; i < n; i++)
      mpf_clear(scale[i]);
    free(scale);
    rownum = NULL;
    mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
    clear_mp(tempComp);

    return MS_NOSOLUTION;
  }

  // calculate x->coord[n-1]
  div_mp(&x->coord[l], &tempVec->coord[pivot], &LU->entry[pivot][l]);

  for (i = n - 2; i >= 0; i--)
  { // calculate x->coord[i]
    k = rownum[i];
    set_mp(tempComp, &tempVec->coord[k]);
    for (j = i + 1; j < n; j++)
      sub_mul_mp(tempComp, &LU->entry[k][j], &x->coord[j]);
    div_mp(&x->coord[i], tempComp, &LU->entry[k][i]);
  }

  clear_vec_mp(tempVec);
  for (i = 0; i < n; i++)
    mpf_clear(scale[i]);
  free(scale);
  rownum = NULL;
  mpf_clear(max); mpf_clear(tempMPF); mpf_clear(prevNorm);
  clear_mp(tempComp);

  return 0;
}

int matrixSolve_from_LU_mp(mat_mp X, mat_mp LU, int *rownum, mat_mp B)
/***************************************************************\
* USAGE: used after finding the LU-decomp by LU_matrixSovle_mp  *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: solves (L*U)X = R*B, where R is stored in rownum       *
* and L and U are stored in LU                                  *
\***************************************************************/
{
  int i, j, k, l, m, pivot, n = LU->rows, p = B->cols;
  mat_mp Y;
  comp_mp tempComp;

  if (LU->rows != LU->cols)
  {
    printf("ERROR: The matrix is not square in matrixSolve_from_LU_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (LU->rows != B->rows)
  {
    printf("ERROR: The sizes in matrixSolve_from_LU_mp do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_mp(tempComp);

  // setup X & Y - same size as B
  change_size_mat_mp(X, n, p);
  init_mat_mp(Y, n, p);
  X->rows = Y->rows = n;
  X->cols = Y->cols = p;

  // first compute LY = RB
  for (k = 0; k < p; k++)
  { // calculate the kth column of Y
    set_mp(&Y->entry[0][k], &B->entry[rownum[0]][k]);

    for (i = 1; i < n; i++)
    { // calculate Y[i][k]
      l = rownum[i];
      neg_mp(&Y->entry[i][k], &B->entry[l][k]);
      for (j = 0; j < i; j++)
        sum_mul_mp(&Y->entry[i][k], &LU->entry[l][j], &Y->entry[j][k]);
      neg_mp(&Y->entry[i][k], &Y->entry[i][k]);
    }
  }

  // then compute UX = Y
  l = n - 1;
  pivot = rownum[l];
  for (k = 0; k < p; k++)
  { // calculate the kth column of X
    div_mp(&X->entry[l][k], &Y->entry[l][k], &LU->entry[pivot][l]);
    for (i = n - 2; i >= 0; i--)
    {
      set_mp(tempComp, &Y->entry[i][k]);
      m = rownum[i];
      for (j = i + 1; j < n; j++)
        sub_mul_mp(tempComp, &LU->entry[m][j], &X->entry[j][k]);
      div_mp(&X->entry[i][k], tempComp, &LU->entry[m][i]);
    }
  }

  clear_mat_mp(Y);
  clear_mp(tempComp);

  return 0;
}

