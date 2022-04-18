// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

#define MS_NOSOLUTION -1

// This file contains the functions that the QR decomposition needs

///////// DOUBLE PRECISION /////////////

int QR_d_prec(mat_d Q, mat_d R, mat_d P, mat_d A, int preSortRows)
/***************************************************************\
* USAGE: finds A * P = Q * R                                    *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order and           *
* incorporate this into Q                                       *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
{
  int retVal;

  retVal = QR_d(Q, R, P, A, 1e-14, 1e-20, 1e13, preSortRows);

  return retVal;
}

int QR_d(mat_d Q, mat_d R, mat_d P, mat_d A, double tol_pivot, double tol_sign, double largeChange, int preSortRows)
/***************************************************************\
* USAGE: finds A * P = Q * R                                    *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order and           *
* incorporate this into Q                                       *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
// if presorting the rows, we want max{|A[1][j]|^2 : j} >= max{|A[2][j]|^2 : j} >= .. >= max{|A[m][j]|^2 : j}
{
  int i, j, k, rank, pivot, m = A->rows, n = A->cols;
  double max, prevNorm, tempD;
  comp_d tempComp, tempComp2;

  int *perm = NULL, *rowPerm = NULL;
  double *norm = NULL, *rowMax = NULL;
  vec_d u, z;

  // setup Q, R, P
  change_size_mat_d(Q, m, m);
  change_size_mat_d(R, m, n);
  change_size_mat_d(P, n, n);
  Q->rows = Q->cols = R->rows = m;
  R->cols = P->rows = P->cols = n;

  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (R == A)
  {
    printf("ERROR: The matrix R cannot be the matrix A in QR_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m == 0 || n == 0) 
  { // if there are no rows/cols, then just return
    return 0;
  }
  // so we can assume that m >= n >= 1

  init_vec_d(u, m);
  init_vec_d(z, m);

  // find A * P = Q * R, where P is n x n, Q is m x m & R is m x n

  // initialize perm, rowPerm, norm & rowMax
  perm = (int *)bmalloc(n * sizeof(int));
  rowPerm = (int *)bmalloc(m * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));
  rowMax = (double *)bmalloc(m * sizeof(double));
  for (i = 0; i < m; i++)
  {
    rowPerm[i] = i;
    rowMax[i] = 0;
    if (i < n)
    {
      perm[i] = i;
      norm[i] = 0;
    }
  }

  // find norms for each col, infinite norm of each row and find the row that is the maximum
  max = pivot = 0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    { // find the size of this [i][j] entry
      tempD = norm_sqr_d(&A->entry[i][j]);
      norm[j] += tempD;

      if (tempD > rowMax[i])
        rowMax[i] = tempD;      
    }
    norm[j] = sqrt(norm[j]);
    // check to see if it is the maximum twoNorm thus far
    if (max < norm[j])
    {
      pivot = j;
      max = norm[j];
    }
  }
  prevNorm = max; // initialize prevNorm to max so that we have an idea of the size of numbers that we are dealing with

  if (preSortRows)
  { // adjust rowPerm based on sorting rowMax
    for (i = 0; i < m; i++)
      for (j = i+1; j < m; j++)
        if (rowMax[i] < rowMax[j])
        { // swap i & j
          k = rowPerm[i];
          rowPerm[i] = rowPerm[j];
          rowPerm[j] = k;
          
          tempD = rowMax[i];
          rowMax[i] = rowMax[j];
          rowMax[j] = tempD;
        }
  }
  // initialize Q to the permutation matrix associated with rowPerm - if we are not doing row swaps, Q will be Id
  convertToP_d(Q, rowPerm, m);
  
  // setup R so that R = Q^T * A 
  for (i = 0; i < m; i++)
  {
    k = rowPerm[i];
    for (j = 0; j < n; j++)
    {
      set_d(&R->entry[i][j], &A->entry[k][j]);
    }
  }

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  {
    // check to make sure pivot element is not too small
    if (norm[perm[pivot]] < tol_pivot || prevNorm > norm[perm[pivot]] * largeChange)
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
    genhh_d(u, R, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & R
    findZ_mat_d(z, R, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in u & every entry in z
    apphh_Z_mat_d(R, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

    // update Q = Q * U_k
    // to improve performance, we find z = Q*u, where we take into account the size of u
    increase_size_vec_d(z, m);
    for (i = 0; i < m; i++)
    {
      set_zero_d(&z->coord[i]);
      for (j = k; j < m; j++)
      {
        sum_mul_d(&z->coord[i], &Q->entry[i][j], &u->coord[j-k]);
      }
    }

    // then update Q = Q - z*u^H = Q + z * (-u^H)
    for (i = 0; i < m; i++)
      for (j = k; j < m; j++)
      { // tempComp = -conj(u[j-k])
        tempComp->r = -u->coord[j-k].r;
        tempComp->i = u->coord[j-k].i;
        // Q[i][j] += z[i] * tempComp
        sum_mul_d(&Q->entry[i][j], &z->coord[i], tempComp);
      }

    // update the norms and find the one that is the maximum
    pivot = k+1;
    max = 0;
    for (j = k+1; j < n; j++)
    {
      rank = perm[j];
      norm[rank] = 0;
      for (i = k+1; i < m; i++)
        norm[rank] += norm_sqr_d(&R->entry[i][rank]);
      norm[rank] = sqrt(norm[rank]);
      if (max < norm[rank])
      {
        max = norm[rank];
        pivot = j;
      }
    }
  }
  // so the rank of A is k and thus corank is n - k
  rank = k;

  // now we need to convert perm to P
  convertToP_d(P, perm, n);

  // convert R to upper triangular matrix using the swaps in perm 
  for (j = 0; j < n; j++)
  { // find which column should be at j
    for (k = 0; k < n; k++)
      if (j == perm[k])
      {
        pivot = k;
        k = n;
       }

    // swap the columns of R if needed - perm[pivot] == j moves to perm[j]
    if (pivot != j)
    {
      for (i = 0; i < m; i++)
      {
        set_d(tempComp, &R->entry[i][j]);
        set_d(&R->entry[i][j], &R->entry[i][perm[j]]);
        set_d(&R->entry[i][perm[j]], tempComp);
      }
      perm[pivot] = perm[j];
      perm[j] = j;
    }
    
    // make the diagonal entries of R real and positive by multiply Q by sign and R by conj(sign)
    if (j < rank)
    {
      tempD = sign_d(tempComp, &R->entry[j][j], tol_sign);
      conjugate_d(tempComp2, tempComp);
      set_double_d(&R->entry[j][j], tempD, 0);
      for (i = 0; i < m; i++)
      {
        mul_d(&Q->entry[i][j], &Q->entry[i][j], tempComp);
        if (i > j && i < n)
        {
          mul_d(&R->entry[j][i], &R->entry[j][i], tempComp2);
        }
      }
    }
  }

  // free the memory
  free(perm);
  free(rowPerm);
  free(norm);
  free(rowMax);

  clear_vec_d(u);
  clear_vec_d(z);

  return (n - rank);
}

void convertToP_d(mat_d P, int *perm, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: convert perm to the permutation matrix P               *
\***************************************************************/
{
  int i, j;

  change_size_mat_d(P, n, n);
  P->rows = P->cols = n;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_double_d(&P->entry[i][j], (i == perm[j]), 0);
    }

  return;
}

void convertToPerm_d(int *perm, mat_d P, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: convert the permutation matrix P to perm               *
\***************************************************************/
{
  int i, j;

  for (i = 0; i < n; i++)
  {
    // find j so that P[i][j] == 1
    for (j = 0; j < n; j++)
      if (P->entry[i][j].r == 1)
      {
        perm[i] = j;
        j = n;
      }
  }

  return;
}

double sign_d(comp_d s, comp_d in, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: return |in|                                    *
* NOTES: s = in / |in| if |in| >= tol, else s = 1               *
\***************************************************************/
{
  double norm, retVal = d_abs_d(in);

  if (retVal < tol)
  {
    set_one_d(s);
  }
  else
  {
    norm = 1 / retVal;
    mul_rdouble_d(s, in, norm);
  }

  return retVal; 
}

void genhh_d(vec_d u, mat_d A, int sr, int c, double norm, int *colnum, double tol_sign)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: u is the householder vector associated to A(sr:rows,c) *
\***************************************************************/
{
  int i, col, size;
  double d;
  comp_d tempComp;

  size = A->rows - sr;
  change_size_vec_d(u, size);
  u->size = size;

  if (size > 0)
  {
    col = colnum[c];
    sign_d(tempComp, &A->entry[sr][col], tol_sign);
    mul_rdouble_d(tempComp, tempComp, norm);

    // find u & update A
    add_d(&u->coord[0], &A->entry[sr][col], tempComp);
    neg_d(&A->entry[sr][col], tempComp);
    d = norm_sqr_d(&u->coord[0]);
    for (i = 1; i < size; i++)
    {
      sr++;
      set_d(&u->coord[i], &A->entry[sr][col]);
      set_zero_d(&A->entry[sr][col]);
      d += norm_sqr_d(&u->coord[i]);
    }

    // now we need to normalize u
    if (d == 0)
    {
      printf("Unable to generate next Householder vector since this would require dividing by 0!\n");
      bexit(ERROR_CONFIGURATION);
    }

    d = sqrt(2.0) / sqrt(d);
    for (i = 0; i < size; i++)
    {
      mul_rdouble_d(&u->coord[i], &u->coord[i], d);
    }
  }

  return;
}

void findZ_mat_d(vec_d z, mat_d C, vec_d u, int sr, int sc, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: z = C(sr:m, sc:n)H * u                                 *
\***************************************************************/
{
  int i, j, col, row, size, rows = u->size;
  comp_d tempComp;
 
  if (u->size != C->rows - sr)
  {
    printf("The sizes of matrices do not match in findZ_mat_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  size = C->cols - sc;
  change_size_vec_d(z, size);
  z->size = size;

  for (j = 0; j < size; j++)
  {
    row = sr;
    col = colnum[sc+j];
    set_zero_d(&z->coord[j]);
    for (i = 0; i < rows; i++, row++)
    {
      conjugate_d(tempComp, &C->entry[row][col]); // doing conjugate transpose
      sum_mul_d(&z->coord[j], tempComp, &u->coord[i]);
    }
  }

  return;
}

void findZ_vec_d(comp_d z, vec_d C, vec_d u, int sr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: z = C(sr:m)H * u                                       *
\***************************************************************/
{
  int j, rows = u->size;
  comp_d tempComp;

  if (u->size != C->size - sr)
  {
    printf("The sizes of matrices do not match in findZ_vec_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  set_zero_d(z);
  for (j = 0; j < rows; j++)
  {
    conjugate_d(tempComp, &C->coord[sr]);
    sum_mul_d(z, tempComp, &u->coord[j]);
    sr++;
  }
  
  return;
}

void apphh_Z_mat_d(mat_d C, vec_d u, vec_d z, int sr, int er, int sc, int ec, int u_sc, int u_ec, int z_sc, int z_ec, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: C = C + u*-zH                                          *
\***************************************************************/
{
  int i, j, ucoord, col, row, rows = er - sr, cols = ec - sc;

  if ((er - sr != u_ec - u_sc) || (z_ec - z_sc != ec - sc))
  {
    printf("The sizes of matrices do not match in apphh_Z_mat_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  for (j = 0; j < cols; j++, z_sc++)
  {
    row = sr;
    col = colnum[sc+j];
    ucoord = u_sc;
    z->coord[z_sc].r = -z->coord[z_sc].r; // conj(neg(z[z_sc+j]))
    for (i = 0; i < rows; i++, row++, ucoord++)
    {
      sum_mul_d(&C->entry[row][col], &u->coord[ucoord], &z->coord[z_sc]);
    }
  }

  return;
}

void apphh_Z_vec_d(vec_d C, vec_d u, comp_d z, int sr, int er, int u_sc, int u_ec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: C = C + u*-zH                                          *
\***************************************************************/
{
  int i, ccoord = sr, ucoord = u_sc, rows = er - sr;

  if (er - sr != u_ec - u_sc)
  {
    printf("The sizes of vecotrs do not match in apphh_Z_vec_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  z->r = -z->r; // conj(neg(z))
  for (i = 0; i < rows; i++, ccoord++, ucoord++)
  {
    sum_mul_d(&C->coord[ccoord], &u->coord[ucoord], z);
  }
  
  return;
}

double cond_est_d(vec_d y, vec_d x, double norm_x, mat_d C, int sr, int er, int col, comp_d gamma, int *colnum, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the norm of y ~= 1 / min sing value of C       *
* NOTES: does Bischof's incremental condtion number to estimate *
* the smallest singular value of C                              *
\***************************************************************/
{
  int i, size = x->size;
  double norm_y, norm_x_sqr, norm_alpha, b, size_gamma_sqr, tempDouble;
  comp_d alpha, tempComp, s, c;

  if (er - sr != size)
  {
    printf("The sizes of vectors do not match in cond_est_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  increase_size_vec_d(y, size + 1);
  y->size = size + 1;

  // calculate alpha = vH*x
  set_zero_d(alpha);
  for (i = 0; i < size; i++)
  {
    conjugate_d(tempComp, &C->entry[sr+i][colnum[col]]);
    sum_mul_d(alpha, tempComp, &x->coord[i]);
  }
  // find |alpha|
  norm_alpha = d_abs_d(alpha);

  // find ||x||^2
  norm_x_sqr = norm_x * norm_x;

  // find |gamma|^2
  size_gamma_sqr = norm_sqr_d(gamma);

  if (size_gamma_sqr == 0)
  { // gamma is too small, so we have a failure
    set_zero_d(&y->coord[size]);
 
    return -1;
  }

  if (norm_alpha < tol)
  { // assume alpha == 0
    if (norm_x_sqr * size_gamma_sqr > 1)
    { // the max is at s=1,c=0 with the norm_y being norm_x
      for (i = 0; i < size; i++)
      {
        set_d(&y->coord[i], &x->coord[i]);
      }
      set_zero_d(&y->coord[size]);
      norm_y = norm_x;
    }
    else // norm_x_sqr * size_gamma_sqr <= 1
    { // the max is at s=0,c=1 with the norm_y = 1 / |gamma|
      for (i = 0; i < size; i++)
      {
        set_zero_d(&y->coord[i]);
      }
      recip_d(&y->coord[size], gamma); // y[size] = 1 / gamma
      norm_y = 1 / sqrt(size_gamma_sqr);
    }
  }
  else // alpha != 0 - this should be the majority of cases!!!
  {
    // calculate b = (|gamma|^2 * norm_x_sqr + |alpha|^2 - 1) / 2
    b = 0.5 * (size_gamma_sqr * norm_x_sqr + norm_alpha * norm_alpha - 1);

    // mu = 1/alpha * b + sign(conj(alpha))*sqrt((b/norm(alpha))^2 + 1) 
    tempDouble = sqrt((b*b) / (norm_alpha*norm_alpha) + 1);
    conjugate_d(c, alpha);
    sign_d(c, c, tol);
    mul_rdouble_d(c, c, tempDouble);
    recip_d(s, alpha);
    mul_rdouble_d(s, s, b);
    add_d(c, c, s);

    // s = mu / sqrt(|mu|^2 + 1)
    tempDouble = 1 / sqrt(norm_sqr_d(c) + 1);
    mul_rdouble_d(s, c, tempDouble);

    // lambda = alpha * mu + 1 & norm_y = sqrt(lambda) / |gamma|
    mul_d(c, alpha, c); // should be real!!
    norm_y = sqrt((c->r + 1) / size_gamma_sqr);

    // c = -1 / sqrt(|mu|^2 + 1)
    set_double_d(c, -tempDouble, 0.0);

    // setup y
    for (i = 0; i < size; i++)
    {
      mul_d(&y->coord[i], s, &x->coord[i]);
    }
    mul_d(&y->coord[size], s, alpha);
    sub_d(&y->coord[size], c, &y->coord[size]);
    div_d(&y->coord[size], &y->coord[size], gamma);
  }

  return norm_y;
}

///////// DO QR WITHOUT COMPUTING Q ////////////

int QR_R_d(mat_d R, mat_d P, mat_d A, double tol_pivot, double tol_sign, double largeChange, int preSortRows)
/***************************************************************\
* USAGE: finds A * P = Q * R - without finding Q explicitly     *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order               *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
// if presorting the rows, we want max{|A[1][j]|^2 : j} >= max{|A[2][j]|^2 : j} >= .. >= max{|A[m][j]|^2 : j}
{
  int i, j, k, rank, pivot, m = A->rows, n = A->cols;
  double max, prevNorm, tempD;
  comp_d tempComp;

  int *perm = NULL, *rowPerm = NULL;
  double *norm = NULL, *rowMax = NULL;
  vec_d u, z;

  // setup R, P
  change_size_mat_d(R, m, n);
  change_size_mat_d(P, n, n);
  R->rows = m;
  R->cols = P->rows = P->cols = n;

  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_R_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (R == A)
  {
    printf("ERROR: The matrix R cannot be the matrix A in QR_R_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m == 0 || n == 0)
  { // if there are no rows/cols, then just return
    return 0;
  }
  // so we can assume that m >= n >= 1

  init_vec_d(u, m);
  init_vec_d(z, m);

  // find A * P = Q * R, where P is n x n, Q would be m x m & R is m x n

  // initialize perm, rowPerm, norm & rowMax
  perm = (int *)bmalloc(n * sizeof(int));
  rowPerm = (int *)bmalloc(m * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));
  rowMax = (double *)bmalloc(m * sizeof(double));
  for (i = 0; i < m; i++)
  {
    rowPerm[i] = i;
    rowMax[i] = 0;
    if (i < n)
    {
      perm[i] = i;
      norm[i] = 0;
    }
  }

  // find norms for each col, infinite norm of each row and find the row that is the maximum
  max = pivot = 0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    { // find the size of this [i][j] entry
      tempD = norm_sqr_d(&A->entry[i][j]);
      norm[j] += tempD;

      if (tempD > rowMax[i])
        rowMax[i] = tempD;
    }
    norm[j] = sqrt(norm[j]);
    // check to see if it is the maximum twoNorm thus far
    if (max < norm[j])
    {
      pivot = j;
      max = norm[j];
    }
  }
  prevNorm = max; // initialize prevNorm to max so that we have an idea of the size of numbers that we are dealing with

  if (preSortRows)
  { // adjust rowPerm based on sorting rowMax
    for (i = 0; i < m; i++)
      for (j = i+1; j < m; j++)
        if (rowMax[i] < rowMax[j])
        { // swap i & j
          k = rowPerm[i];
          rowPerm[i] = rowPerm[j];
          rowPerm[j] = k;

          tempD = rowMax[i];
          rowMax[i] = rowMax[j];
          rowMax[j] = tempD;
        }
  }

  // setup R
  for (i = 0; i < m; i++)
  {
    k = rowPerm[i];
    for (j = 0; j < n; j++)
    {
      set_d(&R->entry[i][j], &A->entry[k][j]);
    }
  }

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  {
    // check to make sure pivot element is not too small
    if (norm[perm[pivot]] < tol_pivot || prevNorm > norm[perm[pivot]] * largeChange)
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
    genhh_d(u, R, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & R
    findZ_mat_d(z, R, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in u & every entry in z
    apphh_Z_mat_d(R, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

    // update the norms and find the one that is the maximum
    pivot = k+1;
    max = 0;
    for (j = k+1; j < n; j++)
    {
      norm[perm[j]] = 0;
      for (i = k+1; i < m; i++)
        norm[perm[j]] += norm_sqr_d(&R->entry[i][perm[j]]);
      norm[perm[j]] = sqrt(norm[perm[j]]);
      if (max < norm[perm[j]])
      {
        max = norm[perm[j]];
        pivot = j;
      }
    }
  }
  // so the rank of A is k and thus corank is n - k
  rank = k;

  // now we need to convert perm to P
  convertToP_d(P, perm, n);

  // convert R to upper triangular matrix using the swaps in perm
  for (j = 0; j < n; j++)
  { // find which column should be at j
    for (k = 0; k < n; k++)
      if (j == perm[k])
      {
        pivot = k;
        k = n;
      }

    // swap the columns of R if needed - perm[pivot] == j moves to perm[j]
    if (pivot != j)
    {
      for (i = 0; i < m; i++)
      {
        set_d(tempComp, &R->entry[i][j]);
        set_d(&R->entry[i][j], &R->entry[i][perm[j]]);
        set_d(&R->entry[i][perm[j]], tempComp);
      }
      perm[pivot] = perm[j];
      perm[j] = j;
    }

    // make the diagonal entries of R real and positive by multiply R by conj(sign)
    if (j < rank)
    {
      tempD = sign_d(tempComp, &R->entry[j][j], tol_sign);
      conjugate_d(tempComp, tempComp);
      set_double_d(&R->entry[j][j], tempD, 0);
      for (i = j+1; i < n; i++)
      {
        mul_d(&R->entry[j][i], &R->entry[j][i], tempComp);
      }
    }
  }

  // free the memory
  free(perm);
  free(rowPerm);
  free(norm);
  free(rowMax);

  clear_vec_d(u);
  clear_vec_d(z);

  return (n - rank);
}

///////// MULTI PRECISION /////////////

int QR_mp_prec(mat_mp Q, mat_mp R, mat_mp P, mat_mp A, int preSortRows, int curr_prec)
/***************************************************************\
* USAGE: finds A * P = Q * R                                    *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order and           *
* incorporate this into Q                                       *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
{
  int retVal, num_digits = prec_to_digits(curr_prec);
  size_t size;
  char *str = NULL;
  mpf_t tol_pivot, tol_sign, largeChange;

  mpf_init(tol_pivot); 
  mpf_init(tol_sign);
  mpf_init(largeChange);

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

  retVal = QR_mp(Q, R, P, A, tol_pivot, tol_sign, largeChange, preSortRows);

  mpf_clear(tol_pivot);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);
  free(str);

  return retVal;
}

int QR_mp(mat_mp Q, mat_mp R, mat_mp P, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange, int preSortRows)
/***************************************************************\
* USAGE: finds A * P = Q * R                                    *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order and           *
* incorporate this into Q                                       *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
// if presorting the rows, we want max{|A[1][j]|^2 : j} >= max{|A[2][j]|^2 : j} >= .. >= max{|A[m][j]|^2 : j}
{
  int i, j, k, rank, pivot, m = A->rows, n = A->cols;
  int *perm = NULL, *rowPerm = NULL;
  comp_mp tempComp, tempComp2;
  mpf_t tempMPF1, tempMPF2, max, prevNorm;
  mpf_t *norm = NULL, *rowMax = NULL;
  vec_mp u, z;

  // setup Q, R, P
  change_size_mat_mp(Q, m, m);
  change_size_mat_mp(R, m, n);
  change_size_mat_mp(P, n, n);
  Q->rows = Q->cols = R->rows = m;
  R->cols = P->rows = P->cols = n;

  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (R == A)
  {
    printf("ERROR: The matrix R cannot be the matrix A in QR_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m == 0 || n == 0) 
  { // if there are no rows/cols, then just return
    return 0;
  }
  // so we can assume that m >= n >= 1

  // initialize MP
  init_mp(tempComp); init_mp(tempComp2);
  mpf_init(tempMPF1); mpf_init(tempMPF2); mpf_init(max); mpf_init(prevNorm);
  init_vec_mp(u, m); init_vec_mp(z, m);

  // find A * P = Q * R, where P is n x n, Q is m x m & R is m x n

  // initialize perm, rowPerm, norm & rowMax
  perm = (int *)bmalloc(n * sizeof(int));
  rowPerm = (int *)bmalloc(m * sizeof(int));
  rowMax = (mpf_t *)bmalloc(m * sizeof(mpf_t));
  norm = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (i = 0; i < m; i++)
  {
    rowPerm[i] = i;
    mpf_init_set_ui(rowMax[i], 0);
    if (i < n)
    {
      perm[i] = i;
      mpf_init_set_ui(norm[i], 0);
    }
  }

  // find norms for each col, infinite norm of each row and find the row that is the maximum
  pivot = 0;
  mpf_set_ui(max, 0);
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    { // find the size of this [i][j] entry
      mpf_mul(tempMPF1, A->entry[i][j].r, A->entry[i][j].r);
      mpf_mul(tempMPF2, A->entry[i][j].i, A->entry[i][j].i);
      mpf_add(tempMPF1, tempMPF1, tempMPF2);
      mpf_add(norm[j], norm[j], tempMPF1);

      if (mpf_cmp(tempMPF1, rowMax[i]) > 0)
        mpf_set(rowMax[i], tempMPF1);
    }
    mpf_sqrt(norm[j], norm[j]);
    // check to see if it is the maximum twoNorm thus far
    if (mpf_cmp(max, norm[j]) < 0)
    {
      pivot = j;
      mpf_set(max, norm[j]);
    }
  }
  mpf_set(prevNorm, max); // initialize prevNorm to max so that we have an idea of the size of numbers that we are dealing with

  if (preSortRows)
  { // adjust rowPerm based on sorting rowMax
    for (i = 0; i < m; i++)
      for (j = i+1; j < m; j++)
        if (rowMax[i] < rowMax[j])
        { // swap i & j
          k = rowPerm[i];
          rowPerm[i] = rowPerm[j];
          rowPerm[j] = k;

          mpf_set(tempMPF1, rowMax[i]);
          mpf_set(rowMax[i], rowMax[j]);
          mpf_set(rowMax[j], tempMPF1);
        }
  }
  // initialize Q to the permutation matrix associated with rowPerm - if we are not doing row swaps, Q will be Id
  convertToP_mp(Q, rowPerm, m);

  // setup R so that R = Q^T * A
  for (i = 0; i < m; i++)
  {
    k = rowPerm[i];
    for (j = 0; j < n; j++)
    {
      set_mp(&R->entry[i][j], &A->entry[k][j]);
    }
  }

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  {
    // check to make sure pivot element is not too small
    if (checkGood_mp(norm[perm[pivot]], prevNorm, tol_pivot, largeChange))
      mpf_set(prevNorm, norm[perm[pivot]]);
    else
      break;

    // swap the columns, if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_mp(u, R, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & R
    findZ_mat_mp(z, R, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in u & every entry in z
    apphh_Z_mat_mp(R, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

    // update Q = Q * U_k
    // to improve performance, we find z = Q*u, where we take into account the size of u
    increase_size_vec_mp(z, m);
    for (i = 0; i < m; i++)
    {
      set_zero_mp(&z->coord[i]);
      for (j = k; j < m; j++)
      {
        sum_mul_mp(&z->coord[i], &Q->entry[i][j], &u->coord[j-k]);
      }
    }

    // then update Q = Q - z*u^H = Q + z * (-u^H)
    for (i = 0; i < m; i++)
      for (j = k; j < m; j++)
      { // tempComp = -conj(u[j-k])
        mpf_neg(tempComp->r, u->coord[j-k].r);
        mpf_set(tempComp->i, u->coord[j-k].i);
        // Q[i][j] += z[i] * tempComp
        sum_mul_mp(&Q->entry[i][j], &z->coord[i], tempComp);
      }

    // update the norms and find the one that is the maximum
    pivot = k+1;
    mpf_set_ui(max, 0);
    for (j = k+1; j < n; j++)
    {
      mpf_set_ui(norm[perm[j]], 0);
      for (i = k+1; i < m; i++)
      {
        mpf_mul(tempMPF1, R->entry[i][perm[j]].r, R->entry[i][perm[j]].r);
        mpf_mul(tempMPF2, R->entry[i][perm[j]].i, R->entry[i][perm[j]].i);
        mpf_add(tempMPF1, tempMPF1, tempMPF2);
        mpf_add(norm[perm[j]], norm[perm[j]], tempMPF1);
      }
      mpf_sqrt(norm[perm[j]], norm[perm[j]]);
      // check to see if it is the maximum twoNorm thus far
      if (mpf_cmp(max, norm[perm[j]]) < 0)
      {
        mpf_set(max, norm[perm[j]]);
        pivot = j;
      }
    }
  }
  // so the rank of A is k and thus corank is n - k
  rank = k;

  // now we need to convert perm to P
  convertToP_mp(P, perm, n);

  // convert R to the upper triangular matrix using the swaps in perm
  for (j = 0; j < n; j++)
  { // find which column should be at j
    for (k = 0; k < n; k++)
    { // find which column should be at j
      if (j == perm[k])
      {
        pivot = k;
        k = n;
      }
    }

    // swap the columns of R if needed - perm[pivot] == j moves to perm[j]
    if (pivot != j)
    {
      for (i = 0; i < n; i++)
      {
        set_mp(tempComp, &R->entry[i][j]);
        set_mp(&R->entry[i][j], &R->entry[i][perm[j]]);
        set_mp(&R->entry[i][perm[j]], tempComp);
      }
      perm[pivot] = perm[j];
      perm[j] = j;
    }

    // make the diagonal entries of R real and positive by multiply Q by sign and R by conj(sign)
    if (j < rank)
    {
      sign_mp2(tempComp, &R->entry[j][j], tol_sign);
      // normalze R[j][j]
      mpf_mul(tempMPF1, R->entry[j][j].r, R->entry[j][j].r);
      mpf_mul(tempMPF2, R->entry[j][j].i, R->entry[j][j].i);
      mpf_add(tempMPF1, tempMPF1, tempMPF2);
      mpf_sqrt(R->entry[j][j].r, tempMPF1);
      mpf_set_ui(R->entry[j][j].i, 0);

      // find conj(tempComp)
      conjugate_mp(tempComp2, tempComp);

      // update Q & R
      for (i = 0; i < m; i++)
      {
        mul_mp(&Q->entry[i][j], &Q->entry[i][j], tempComp);
        if (i > j && i < n)
        {
          mul_mp(&R->entry[j][i], &R->entry[j][i], tempComp2);
        }
      }
    }
  }

  // clear MP
  clear_mp(tempComp); clear_mp(tempComp2);
  mpf_clear(tempMPF1); mpf_clear(tempMPF2); mpf_clear(max); mpf_clear(prevNorm);
  clear_vec_mp(u); clear_vec_mp(z);

  // clear norm & rowMax
  for (i = m - 1; i >= 0; i--)
  {
    mpf_clear(rowMax[i]);
    if (i < n)
    {
      mpf_clear(norm[i]);
    }
  }
  free(norm);
  free(rowMax);

  free(perm);
  free(rowPerm);

  return (n - rank);
}

void convertToP_mp(mat_mp P, int *perm, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: convert perm to the permutation matrix P               *
\***************************************************************/
{
  int i, j;
  
  change_size_mat_mp(P, n, n);
  P->rows = P->cols = n;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_double_mp(&P->entry[i][j], (i == perm[j]), 0);
    }

  return;
}

void convertToPerm_mp(int *perm, mat_mp P, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: convert the permutation matrix P to perm               *
\***************************************************************/
{
  int i, j;

  for (i = 0; i < n; i++)
  {
    // find j so that P[i][j] == 1
    for (j = 0; j < n; j++)
      if (mpf_cmp_ui(P->entry[i][j].r, 1) == 0)
      {
        perm[i] = j;
        j = n;
      }
  }

  return;
}

double sign_mp(comp_mp s, comp_mp in, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: return |in|                                    *
* NOTES: s = in / |in| if |in| >= tol, else s = 1               *
\***************************************************************/
{
  double retVal;
  mpf_t t1, t2;

  mpf_init(t1); mpf_init(t2);

  mpf_mul(t1, in->r, in->r);
  mpf_mul(t2, in->i, in->i);
  mpf_add(t1, t1, t2);
  mpf_sqrt(t1, t1);

  retVal = mpf_get_d(t1);

  if (retVal < tol)
  {
    set_one_mp(s);
  }
  else
  {
    mpf_ui_div(t1, 1, t1);
    mpf_mul(s->r, in->r, t1);
    mpf_mul(s->i, in->i, t1);
  }

  mpf_clear(t1); mpf_clear(t2);

  return retVal;
}

double sign_mp2(comp_mp s, comp_mp in, mpf_t tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: return |in|                                    *
* NOTES: s = in / |in| if |in| >= tol, else s = 1               *
\***************************************************************/
{
  double retVal;
  mpf_t t1, t2;

  mpf_init(t1); mpf_init(t2);
 
  mpf_mul(t1, in->r, in->r);
  mpf_mul(t2, in->i, in->i);
  mpf_add(t1, t1, t2);
  mpf_sqrt(t1, t1);

  retVal = mpf_get_d(t1);

  if (mpf_cmp(t1, tol) < 0)
  {
    set_one_mp(s);
  }
  else
  {
    mpf_ui_div(t1, 1, t1);
    mpf_mul(s->r, in->r, t1);
    mpf_mul(s->i, in->i, t1);
  }

  mpf_clear(t1); mpf_clear(t2);

  return retVal;
}

void sign_mp3(mpf_t norm, comp_mp s, comp_mp in, mpf_t tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: return |in|                                    *
* NOTES: s = in / |in| if |in| >= tol, else s = 1               *
\***************************************************************/
{
  mpf_t t;

  mpf_init(t);

  mpf_mul(t, in->r, in->r);
  mpf_mul(t, in->i, in->i);
  mpf_add(t, t, t);
  mpf_sqrt(norm, t);

  if (mpf_cmp(norm, tol) < 0)
  {
    set_one_mp(s);
  }
  else
  {
    mpf_ui_div(t, 1, norm);
    mpf_mul(s->r, in->r, t);
    mpf_mul(s->i, in->i, t);
  }

  mpf_clear(t);

  return;
}

void genhh_mp(vec_mp u, mat_mp A, int sr, int c, mpf_t norm, int *colnum, mpf_t tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: u is the householder vector associated to A(sr:rows,c) *
\***************************************************************/
{
  int i, size;

  size = A->rows - sr;
  change_size_vec_mp(u, size);
  u->size = size;

  if (size > 0)
  {
    mpf_t tempMPF;
    comp_mp tempComp;

    mpf_init(tempMPF);
    init_mp(tempComp);

    sign_mp2(tempComp, &A->entry[sr][colnum[c]], tol);
    mpf_mul(tempComp->r, tempComp->r, norm);
    mpf_mul(tempComp->i, tempComp->i, norm);

    // find u & update A
    add_mp(&u->coord[0], &A->entry[sr][colnum[c]], tempComp);
    neg_mp(&A->entry[sr][colnum[c]], tempComp);

    mpf_mul(tempComp->r, u->coord[0].r, u->coord[0].r);
    mpf_mul(tempComp->i, u->coord[0].i, u->coord[0].i);
    mpf_add(tempMPF, tempComp->r, tempComp->i);
    for (i = 1; i < size; i++)
    {
      set_mp(&u->coord[i], &A->entry[sr+i][colnum[c]]);
      set_zero_mp(&A->entry[sr+i][colnum[c]]);

      mpf_mul(tempComp->r, u->coord[i].r, u->coord[i].r);
      mpf_add(tempMPF, tempMPF, tempComp->r);
      mpf_mul(tempComp->r, u->coord[i].i, u->coord[i].i);
      mpf_add(tempMPF, tempMPF, tempComp->r);
    }

    // now we need to normalize u
    if (mpfr_zero_p(tempMPF))
    { // maybe we should think about returning u as the zero vector ??
      printf("Unable to generate next Householder vector since this would require dividing by 0!\n");
      bexit(ERROR_CONFIGURATION);
    }

    // d = sqrt(2) / norm(u)
    mpf_ui_div(tempMPF, 2, tempMPF);
    mpf_sqrt(tempMPF, tempMPF);
    for (i = 0; i < size; i++)
    {
      mpf_mul(u->coord[i].r, u->coord[i].r, tempMPF);
      mpf_mul(u->coord[i].i, u->coord[i].i, tempMPF);
    }

    mpf_clear(tempMPF);
    clear_mp(tempComp);
  }

  return;
}

void findZ_mat_mp(vec_mp z, mat_mp C, vec_mp u, int sr, int sc, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: z = C(sr:m, sc:n)H * u                                 *
\***************************************************************/
{
  int i, j, size, rows = u->size;
  comp_mp tempComp;

  if (u->size != C->rows - sr)
  {
    printf("The sizes of matrices do not match in findZ_mat_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_mp(tempComp);

  size = C->cols - sc;
  change_size_vec_mp(z, size);
  z->size = size;

  for (j = 0; j < size; j++)
  {
    set_zero_mp(&z->coord[j]);
    for (i = 0; i < rows; i++)
    {
      conjugate_mp(tempComp, &C->entry[sr+i][colnum[sc+j]]); // doing conjugate transpose
      sum_mul_mp(&z->coord[j], tempComp, &u->coord[i]);
    }
  }

  clear_mp(tempComp);

  return;
}

void findZ_vec_mp(comp_mp z, vec_mp C, vec_mp u, int sr)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: z = C(sr:m)H * u                                       *
\***************************************************************/
{
  int j, rows = u->size;
  comp_mp tempComp;

  if (u->size != C->size - sr)
  {
    printf("The sizes of matrices do not match in findZ_vec_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_mp(tempComp);

  set_zero_mp(z);
  for (j = 0; j < rows; j++)
  {
    conjugate_mp(tempComp, &C->coord[sr+j]);
    sum_mul_mp(z, tempComp, &u->coord[j]);
  }

  clear_mp(tempComp);

  return;
}

void apphh_Z_mat_mp(mat_mp C, vec_mp u, vec_mp z, int sr, int er, int sc, int ec, int u_sc, int u_ec, int z_sc, int z_ec, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: C = C + u*-zH                                          *
\***************************************************************/
{
  int i, j, rows = er - sr, cols = ec - sc;

  if ((er - sr != u_ec - u_sc) || (z_ec - z_sc != ec - sc))
  {
    printf("The sizes of matrices do not match in apphh_Z_mat_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  for (j = 0; j < cols; j++)
  {
    mpf_neg(z->coord[z_sc+j].r, z->coord[z_sc+j].r); // conj(neg(z[z_sc+j]))
    for (i = 0; i < rows; i++)
    {
      sum_mul_mp(&C->entry[sr+i][colnum[sc+j]], &u->coord[u_sc+i], &z->coord[z_sc+j]);
    }
  }

  return;
}

void apphh_Z_vec_mp(vec_mp C, vec_mp u, comp_mp z, int sr, int er, int u_sc, int u_ec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: C = C + u*-zH                                          *
\***************************************************************/
{
  int i, rows = er - sr;

  if (er - sr != u_ec - u_sc)
  {
    printf("The sizes of vecotrs do not match in apphh_Z_vec_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  mpf_neg(z->r, z->r); // conj(neg(z))
  for (i = 0; i < rows; i++)
  {
    sum_mul_mp(&C->coord[sr+i], &u->coord[u_sc+i], z);
  }

  return;
}

void cond_est_mp(vec_mp y, vec_mp x, mpf_t norm_x_y, mat_mp C, int sr, int er, int col, comp_mp gamma, int *colnum, mpf_t tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS: norm_x_y is the norm of x upon entering and the    *
* norm of y upon exiting                                        *
* RETURN VALUES: the norm of y ~= 1 / min sing value of C       *
* NOTES: does Bischof's incremental condtion number to estimate *
* the smallest singular value of C                              *
\***************************************************************/
{
  int i, size = x->size, prec = mpf_get_prec(C->entry[0][0].r);
  mpf_t norm_x_sqr, norm_alpha_sqr, b, size_gamma_sqr, tempD; // since these are used in calculations, they need to be multiprecision
  comp_mp alpha, tempComp, s, c;

  init_mp2(alpha, prec);
  init_mp2(tempComp, prec);
  init_mp2(s, prec);
  init_mp2(c, prec);

  mpf_init2(norm_x_sqr, prec);
  mpf_init2(norm_alpha_sqr, prec);
  mpf_init2(b, prec);
  mpf_init2(size_gamma_sqr, prec);
  mpf_init2(tempD, prec);

  if (er - sr != size)
  {
    printf("The sizes of vectors do not match in cond_est_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  increase_size_vec_mp(y, size + 1);
  y->size = size + 1;

  // calculate alpha = vH*x
  set_zero_mp(alpha);
  for (i = 0; i < size; i++)
  {
    conjugate_mp(tempComp, &C->entry[sr+i][colnum[col]]);
    sum_mul_mp(alpha, tempComp, &x->coord[i]);
  }
  // find |alpha|^2
  mpf_mul(tempComp->r, alpha->r, alpha->r);
  mpf_mul(tempComp->i, alpha->i, alpha->i);
  mpf_add(norm_alpha_sqr, tempComp->r, tempComp->i);

  // find |x|^2
  mpf_mul(norm_x_sqr, norm_x_y, norm_x_y);

  // find |gamma|^2
  mpf_mul(tempComp->r, gamma->r, gamma->r);
  mpf_mul(tempComp->i, gamma->i, gamma->i);
  mpf_add(size_gamma_sqr, tempComp->r, tempComp->i);

  if (mpfr_zero_p(size_gamma_sqr))
  { // gamma is too small, so we have a failure
    set_zero_mp(&y->coord[size]);
    mpf_set_si(norm_x_y, -1);    

    // clear MP
    clear_mp(alpha);
    clear_mp(tempComp);
    clear_mp(s);
    clear_mp(c);

    mpf_clear(norm_x_sqr);
    mpf_clear(norm_alpha_sqr);
    mpf_clear(b);
    mpf_clear(size_gamma_sqr);
    mpf_clear(tempD);

    return;
  }

  // make sure that alpha is large enough
  mpf_sqrt(tempD, norm_alpha_sqr);
  if (mpf_cmp(tempD, tol) < 0)
  {
    mpf_mul(tempD, norm_x_sqr, size_gamma_sqr);
    if (mpf_cmp_ui(tempD, 1) > 0)
    { // the max is at s=1,c=0 with the norm_y being norm_x
      for (i = 0; i < size; i++)
      {
        set_mp(&y->coord[i], &x->coord[i]);
      }
      set_zero_mp(&y->coord[size]);
      // the norm of y is the same as the norm of x, so we leave norm_x_y alone
    }
    else // norm_x_sqr * size_gamma_sqr <= 1
    { // the max is at s=0,c=1 with the norm_y = 1 / |gamma|
      for (i = 0; i < size; i++)
      {
        set_zero_mp(&y->coord[i]);
      }
      recip_mp(&y->coord[size], gamma); // y[size] = 1 / gamma
      mpf_sqrt(tempD, size_gamma_sqr);
      mpf_ui_div(norm_x_y, 1, tempD);
    }
  }
  else // alpha != 0 - this should be the majority of cases!!!
  {
    // calculate b = (|gamma|^2 * norm_x_sqr + |alpha|^2 - 1) / 2
    mpf_mul(b, size_gamma_sqr, norm_x_sqr);
    mpf_add(b, b, norm_alpha_sqr);
    mpf_sub_ui(b, b, 1);
    mpf_div_ui(b, b, 2);

    // mu = 1/alpha * b + sign(conj(alpha))*sqrt((b/norm(alpha))^2 + 1)
    mpf_mul(tempD, b, b);
    mpf_div(tempD, tempD, norm_alpha_sqr);
    mpf_add_ui(tempD, tempD, 1);
    mpf_sqrt(tempD, tempD);

    conjugate_mp(c, alpha);
    sign_mp2(c, c, tol);
    mpf_mul(c->r, c->r, tempD);
    mpf_mul(c->i, c->i, tempD);

    recip_mp(s, alpha);
    mpf_mul(s->r, s->r, b);
    mpf_mul(s->i, s->i, b);

    add_mp(c, c, s);

    // s = mu / sqrt(|mu|^2 + 1)
    mpf_mul(tempComp->r, c->r, c->r);
    mpf_mul(tempComp->i, c->i, c->i);
    mpf_add(tempD, tempComp->r, tempComp->i);
    mpf_add_ui(tempD, tempD, 1);
    mpf_sqrt(tempD, tempD);
    mpf_ui_div(tempD, 1, tempD); // tempD = 1 / sqrt(|mu|^2 + 1)
    mpf_mul(s->r, c->r, tempD);
    mpf_mul(s->i, c->i, tempD);

    // lambda = alpha * mu + 1 & norm_y = sqrt(lambda) / |gamma|
    mul_mp(c, alpha, c); // this should be real
    mpf_add_ui(tempComp->r, c->r, 1);
    mpf_div(norm_x_y, tempComp->r, size_gamma_sqr);
    mpf_sqrt(norm_x_y, norm_x_y);

    // c = -1 / sqrt(|mu|^2 + 1)
    mpf_neg(c->r, tempD);
    mpf_set_ui(c->i, 0);

    // setup y
    for (i = 0; i < size; i++)
    {
      mul_mp(&y->coord[i], s, &x->coord[i]);
    }
    // (c - s*alpha) / gamma
    mul_mp(&y->coord[size], s, alpha);
    sub_mp(&y->coord[size], c, &y->coord[size]);
    div_mp(&y->coord[size], &y->coord[size], gamma);
  }

  // clear MP
  clear_mp(alpha);
  clear_mp(tempComp);
  clear_mp(s);
  clear_mp(c);

  mpf_clear(norm_x_sqr);
  mpf_clear(norm_alpha_sqr);
  mpf_clear(b);
  mpf_clear(size_gamma_sqr);
  mpf_clear(tempD);

  return;
}

///////// DO QR WITHOUT COMPUTING Q ////////////

int QR_R_mp_prec(mat_mp R, mat_mp P, mat_mp A, int preSortRows, int curr_prec)
/***************************************************************\
* USAGE: finds A * P = Q * R - without finding Q explicitly     *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order               *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
{
  int retVal, num_digits = prec_to_digits(curr_prec);
  size_t size;
  char *str = NULL;
  mpf_t tol_pivot, tol_sign, largeChange;

  mpf_init(tol_pivot);
  mpf_init(tol_sign);
  mpf_init(largeChange);

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

  retVal = QR_R_mp(R, P, A, tol_pivot, tol_sign, largeChange, preSortRows);

  mpf_clear(tol_pivot);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);
  free(str);

  return retVal;
}

int QR_R_mp(mat_mp R, mat_mp P, mat_mp A, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange, int preSortRows)
/***************************************************************\
* USAGE: finds A * P = Q * R - without finding Q explicitly     *
* ARGUMENTS: preSortRows == 0 - do not pre sort the rows,       *
* otherwise, presort the rows in decreasing order               *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* that is, retVal == 0 <=> A is full rank                       *
* NOTES: Assume that R != A                                     *
\***************************************************************/
// if presorting the rows, we want max{|A[1][j]|^2 : j} >= max{|A[2][j]|^2 : j} >= .. >= max{|A[m][j]|^2 : j}
{
  int i, j, k, rank, pivot, m = A->rows, n = A->cols;
  int *perm = NULL, *rowPerm = NULL;
  comp_mp tempComp, tempComp2;
  mpf_t tempMPF1, tempMPF2, max, prevNorm;
  mpf_t *norm = NULL, *rowMax = NULL;
  vec_mp u, z;

  // setup R, P
  change_size_mat_mp(R, m, n);
  change_size_mat_mp(P, n, n);
  R->rows = m;
  R->cols = P->rows = P->cols = n;

  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (R == A)
  {
    printf("ERROR: The matrix R cannot be the matrix A in QR_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m == 0 || n == 0)
  { // if there are no rows/cols, then just return
    return 0;
  }
  // so we can assume that m >= n >= 1

  // initialize MP
  init_mp(tempComp); init_mp(tempComp2);
  mpf_init(tempMPF1); mpf_init(tempMPF2); mpf_init(max); mpf_init(prevNorm);
  init_vec_mp(u, m); init_vec_mp(z, m);

  // find A * P = Q * R, where P is n x n, Q would be m x m & R is m x n

  // initialize perm, rowPerm, norm & rowMax
  perm = (int *)bmalloc(n * sizeof(int));
  rowPerm = (int *)bmalloc(m * sizeof(int));
  rowMax = (mpf_t *)bmalloc(m * sizeof(mpf_t));
  norm = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (i = 0; i < m; i++)
  {
    rowPerm[i] = i;
    mpf_init_set_ui(rowMax[i], 0);
    if (i < n)
    {
      perm[i] = i;
      mpf_init_set_ui(norm[i], 0);
    }
  }

  // find norms for each col, infinite norm of each row and find the row that is the maximum
  pivot = 0;
  mpf_set_ui(max, 0);
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    { // find the size of this [i][j] entry
      mpf_mul(tempMPF1, A->entry[i][j].r, A->entry[i][j].r);
      mpf_mul(tempMPF2, A->entry[i][j].i, A->entry[i][j].i);
      mpf_add(tempMPF1, tempMPF1, tempMPF2);
      mpf_add(norm[j], norm[j], tempMPF1);

      if (mpf_cmp(tempMPF1, rowMax[i]) > 0)
        mpf_set(rowMax[i], tempMPF1);
    }
    mpf_sqrt(norm[j], norm[j]);
    // check to see if it is the maximum twoNorm thus far
    if (mpf_cmp(max, norm[j]) < 0)
    {
      pivot = j;
      mpf_set(max, norm[j]);
    }
  }
  mpf_set(prevNorm, max); // initialize prevNorm to max so that we have an idea of the size of numbers that we are dealing with

  if (preSortRows)
  { // adjust rowPerm based on sorting rowMax
    for (i = 0; i < m; i++)
      for (j = i+1; j < m; j++)
        if (rowMax[i] < rowMax[j])
        { // swap i & j
          k = rowPerm[i];
          rowPerm[i] = rowPerm[j];
          rowPerm[j] = k;

          mpf_set(tempMPF1, rowMax[i]);
          mpf_set(rowMax[i], rowMax[j]);
          mpf_set(rowMax[j], tempMPF1);
        }
  }

  // setup R
  for (i = 0; i < m; i++)
  {
    k = rowPerm[i];
    for (j = 0; j < n; j++)
    {
      set_mp(&R->entry[i][j], &A->entry[k][j]);
    }
  }

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  {
    // check to make sure pivot element is not too small
    if (checkGood_mp(norm[perm[pivot]], prevNorm, tol_pivot, largeChange))
      mpf_set(prevNorm, norm[perm[pivot]]);
    else
      break;

    // swap the columns, if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_mp(u, R, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & R
    findZ_mat_mp(z, R, u, k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in u & every entry in z
    apphh_Z_mat_mp(R, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

    // update the norms and find the one that is the maximum
    pivot = k+1;
    mpf_set_ui(max, 0);
    for (j = k+1; j < n; j++)
    {
      mpf_set_ui(norm[perm[j]], 0);
      for (i = k+1; i < m; i++)
      {
        mpf_mul(tempMPF1, R->entry[i][perm[j]].r, R->entry[i][perm[j]].r);
        mpf_mul(tempMPF2, R->entry[i][perm[j]].i, R->entry[i][perm[j]].i);
        mpf_add(tempMPF1, tempMPF1, tempMPF2);
        mpf_add(norm[perm[j]], norm[perm[j]], tempMPF1);
      }
      mpf_sqrt(norm[perm[j]], norm[perm[j]]);
      // check to see if it is the maximum twoNorm thus far
      if (mpf_cmp(max, norm[perm[j]]) < 0)
      {
        mpf_set(max, norm[perm[j]]);
        pivot = j;
      }
    }
  }
  // so the rank of A is k and thus corank is n - k
  rank = k;

  // now we need to convert perm to P
  convertToP_mp(P, perm, n);

  // convert R to the upper triangular matrix using the swaps in perm
  for (j = 0; j < n; j++)
  { // find which column should be at j
    for (k = 0; k < n; k++)
    { // find which column should be at j
      if (j == perm[k])
      {
        pivot = k;
        k = n;
      }
    }

    // swap the columns of R if needed - perm[pivot] == j moves to perm[j]
    if (pivot != j)
    {
      for (i = 0; i < n; i++)
      {
        set_mp(tempComp, &R->entry[i][j]);
        set_mp(&R->entry[i][j], &R->entry[i][perm[j]]);
        set_mp(&R->entry[i][perm[j]], tempComp);
      }
      perm[pivot] = perm[j];
      perm[j] = j;
    }

    // make the diagonal entries of R real and positive by multiply R by conj(sign)
    if (j < rank)
    {
      sign_mp2(tempComp, &R->entry[j][j], tol_sign);
      // normalze R[j][j]
      mpf_mul(tempMPF1, R->entry[j][j].r, R->entry[j][j].r);
      mpf_mul(tempMPF2, R->entry[j][j].i, R->entry[j][j].i);
      mpf_add(tempMPF1, tempMPF1, tempMPF2);
      mpf_sqrt(R->entry[j][j].r, tempMPF1);
      mpf_set_ui(R->entry[j][j].i, 0);

      // find conj(tempComp)
      conjugate_mp(tempComp2, tempComp);

      // update R
      for (i = 0; i < m; i++)
        if (i > j && i < n)
        {
          mul_mp(&R->entry[j][i], &R->entry[j][i], tempComp2);
        }
    }
  }

  // clear MP
  clear_mp(tempComp); clear_mp(tempComp2);
  mpf_clear(tempMPF1); mpf_clear(tempMPF2); mpf_clear(max); mpf_clear(prevNorm);
  clear_vec_mp(u); clear_vec_mp(z);

  // clear norm & rowMax
  for (i = m - 1; i >= 0; i--)
  {
    mpf_clear(rowMax[i]);
    if (i < n)
    {
      mpf_clear(norm[i]);
    }
  }
  free(norm);
  free(rowMax);
  free(perm);
  free(rowPerm);

  return (n - rank);
}


/////// matrix solving using QR decompostion //////////

int matrixSolve_from_QR_d(vec_d x, mat_d Q_trans, mat_d R, int *perm, vec_d b)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, retVal = 0, m = R->rows, n = R->cols;
  comp_d sum;
  vec_d Qb, y;

  // check that we have  m == n & the size of b is correct
  if (m != n)
  {
    printf("ERROR: The matrix is not square in matrixSolve_from_QR_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != Q_trans->rows)
  {
    printf("ERROR: The sizes in matrixSolve_from_QR_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // m == n
  
  // make sure m == n > 0
  if (n == 0)
  {
    x->size = 0;
    return 0;
  }

  // initialize Qb & y
  init_vec_d(Qb, m);
  init_vec_d(y, n);
  Qb->size = m;
  y->size = n;

  // compute Qb = Q_trans * b
  mul_mat_vec_d(Qb, Q_trans, b);

  // compute y such that R*y = Qb
  for (i = n - 1; i >= 0; i--)
  { // compute y[i]
    set_zero_d(sum);
    for (j = n - 1; j > i; j--)
    {
      sum_mul_d(sum, &R->entry[i][j], &y->coord[j]);
    }
    sub_d(sum, &Qb->coord[i], sum);
    div_d(&y->coord[i], sum, &R->entry[i][i]);
  }

  // setup x
  change_size_vec_d(x, n);
  x->size = n;
  for (i = 0; i < n; i++)
  { // compute x[i] = y[perm[i]]
    set_d(&x->coord[i], &y->coord[perm[i]]); 

    // check for errors
    if (isnan(x->coord[i].r) || isnan(x->coord[i].i))
    { // error
      retVal = MS_NOSOLUTION;
    }
  }

  // clear memory
  clear_vec_d(y);

  return retVal;
}

int matrixSolve_from_QR_mp(vec_mp x, mat_mp Q_trans, mat_mp R, int *perm, vec_mp b)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if acceptable answer, non-zero otherwise     *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, retVal = 0, m = R->rows, n = R->cols;
  comp_mp sum;
  vec_mp Qb, y;

  // check that we have  m == n & the size of b is correct
  if (m != n)
  {
    printf("ERROR: The matrix is not square in matrixSolve_from_QR_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m != Q_trans->rows)
  {
    printf("ERROR: The sizes in matrixSolve_from_QR_d do not match!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // m == n

  // make sure m == n > 0
  if (n == 0)
  {
    x->size = 0;
    return 0;
  }

  init_mp(sum);

  // initialize Qb & y
  init_vec_mp(Qb, m);
  init_vec_mp(y, n);
  Qb->size = m;
  y->size = n;

  // compute Qb = Q_trans * b
  mul_mat_vec_mp(Qb, Q_trans, b);

  // compute y such that R*y = Qb
  for (i = n - 1; i >= 0; i--)
  { // compute y[i]
    set_zero_mp(sum);
    for (j = n - 1; j > i; j--)
    {
      sum_mul_mp(sum, &R->entry[i][j], &y->coord[j]);
    }
    sub_mp(sum, &Qb->coord[i], sum);
    div_mp(&y->coord[i], sum, &R->entry[i][i]);
  }

  // setup x
  change_size_vec_mp(x, n);
  x->size = n;
  for (i = 0; i < n; i++)
  { // compute x[i] = y[perm[i]]
    set_mp(&x->coord[i], &y->coord[perm[i]]);

    // check for errors
    if (!(mpfr_number_p(x->coord[i].r) && mpfr_number_p(x->coord[i].i)))
    { // error
      retVal = MS_NOSOLUTION;
    }
  }

  // clear memory
  clear_vec_d(y);

  return retVal;
}


