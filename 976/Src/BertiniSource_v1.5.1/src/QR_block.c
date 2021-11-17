// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

int QR_R_block_d(mat_d R, mat_d P, mat_d A, int rowTop, int colTop, double tol_pivot, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: Takes a block triangular matrix A and turns it into a  *
* upper triangular matrix R using Householder transformations   *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES:                                                        *
\***************************************************************/
{ // assume the top block has full column rank (i.e. colTop <= rowTop and has colTop nonzero pivots)
  int i, j, k, l, rank, pivot, m = A->rows, n = A->cols, *perm = NULL;
  double prevNorm, max, tempD, *norm = NULL;
  comp_d tempComp;
  vec_d u, z;
  mat_d tempR;

  // error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_R_block_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QR_R_block_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= rowTop && rowTop <= A->rows))
  {
    printf("ERROR: The top block is not setup properly in QR_R_block_d!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= colTop && colTop <= A->rows))
  {
    printf("ERROR: The top block is not setup properly in QR_R_block_d!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (colTop > rowTop)
  {
    printf("ERROR: The top block is not setup properly in QR_R_block_d!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  if (m == 0 || n == 0)
  { // copy A to R and return
    mat_cp_d(R, A);
    return 0;
  }

  // so we have m >= n >= 1
  init_vec_d(u, m);
  init_vec_d(z, m);
  init_mat_d(tempR, m, n);

  // copy A to tempR
  mat_cp_d(tempR, A);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));
  max = pivot = 0;
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    norm[k] = 0;
    if (k < colTop)
    { // setup norm and find the first pivot
      for (j = 0; j < rowTop; j++)
        norm[k] += norm_sqr_d(&tempR->entry[j][k]);
      norm[k] = sqrt(norm[k]);
      // check to see if it is the maximum twoNorm thus far
      if (max < norm[k])
      {
        max = norm[k];
        pivot = k;
      }
    }
  }

  // do QR decomposition for top block
  tempR->rows = rowTop;
  for (k = 0; k < colTop; k++)
  { // swap the columns if needed
    if (pivot != k)
    { 
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    if (norm[perm[k]] > tol_sign)
    { // generate the Householder quantities
      genhh_d(u, tempR, k, k, norm[perm[k]], perm, tol_sign);
      // find the z associated to u & tempR
      findZ_mat_d(z, tempR, u, k, k+1, perm);
      // apply z to rows k to rowTop & cols k+1 to n of tempR using every entry in u & every entry in z
      apphh_Z_mat_d(tempR, u, z, k, rowTop, k+1, n, 0, u->size, 0, z->size, perm);
    }

    // update the norms and find the one that is the maximum, if needed
    pivot = k + 1;
    max = 0;
    for (j = k + 1; j < colTop; j++)
    {
      l = perm[j];
      norm[l] = 0;
      for (i = k + 1; i < rowTop; i++)
        norm[l] += norm_sqr_d(&tempR->entry[i][l]);
      norm[l] = sqrt(norm[l]);
      if (max < norm[l])
      {
        max = norm[l];
        pivot = j;
      }
    }
  }
  // setup tempR back to full size
  tempR->rows = m;

  // now we need to a normal QR decomposition on the bottom
  // find the norms for each col and find the col that is the maximum remaining
  max = pivot = 0;
  for (k = colTop; k < n; k++)
  { // find the norm for the bottom of the kth column
    norm[k] = 0;
    for (j = rowTop; j < m; j++)
      norm[k] += norm_sqr_d(&tempR->entry[j][k]);
    norm[k] = sqrt(norm[k]);
    // check to see if it is the maximum twoNorm thus far
    if (max < norm[k])
    {
      max = norm[k];
      pivot = k;
    }
  }
  prevNorm = max;

  // do QR decomposition for bottom block
  for (k = colTop; k < n; k++)
  { // check to make sure pivot element is not too small
    if (norm[perm[pivot]] < tol_pivot || prevNorm > norm[perm[pivot]] * largeChange)
      break;
    else
      prevNorm = norm[perm[pivot]];

    // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_d(u, tempR, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & tempR
    findZ_mat_d(z, tempR, u, k, k+1, perm);
    // apply z to rows k to rowTop & cols k+1 to n of tempR using every entry in u & every entry in z
    apphh_Z_mat_d(tempR, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

    // update the norms and find the one that is the maximum
    pivot = k + 1;
    max = 0;
    for (j = k + 1; j < n; j++)
    {
      l = perm[j];
      norm[l] = 0;
      for (i = k + 1; i < m; i++)
        norm[l] += norm_sqr_d(&tempR->entry[i][l]);
      norm[l] = sqrt(norm[l]);
      if (max < norm[l])
      {
        max = norm[l];
        pivot = j;
      }
    }
  }
  // so the rank of A is k and thus corank is n - k
  rank = k;

  // now we need to convert perm to P
  convertToP_d(P, perm, n);

  // convert tempR to R - upper traingular matrix using the swaps in perm
  change_size_mat_d(R, m, n);
  for (j = 0; j < n; j++)
  { // setup column j
    for (i = 0; i < m; i++)
    {
      set_d(&R->entry[i][j], &tempR->entry[i][perm[j]]);
    }

    // make the diagonal entries of R real and positive
    if (j < rank)
    {
      tempD = sign_d(tempComp, &R->entry[j][j], tol_sign);
      conjugate_d(tempComp, tempComp);
      set_double_d(&R->entry[j][j], tempD, 0);
      for (i = j + 1; i < n; i++)
      {
        mul_d(&tempR->entry[j][perm[i]], &tempR->entry[j][perm[i]], tempComp);
      }
    }
  }

  // clear memory
  free(perm);
  free(norm);
  clear_vec_d(u); clear_vec_d(z);
  clear_mat_d(tempR);

  return (n - rank);
}

int QLP_block_pair_d(mat_d L0, mat_d L1, mat_d A0, mat_d A1, int rowTop, int colTop, double tol_pivot, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: Takes a pair of block triangular matrices and produces *
* a QLP^T decomposition where L is lower triangular matrix and  *
* Q & P would be unitary                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES: Work with A0 and update A1 (i.e. A0 is more accurate)  *
\***************************************************************/
{ 
  int i, j, k, l, rank, pivot, m = A0->rows, n = A0->cols, *perm = NULL;
  double prevNorm, max, *norm = NULL;
  vec_d u, z;
  mat_d tempR0, tempR1;

  // error checking
  if (A0->rows != A1->rows || A0->cols != A1->cols)
  {
    printf("ERROR: The matrices need to be the same size in QLP_block_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QLP_block_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QLP_block_pair_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= rowTop && rowTop <= m))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= colTop && colTop <= n))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  // copy A0 to L0 & A1 to L1
  change_size_mat_d(L0, m, n); 
  change_size_mat_d(L1, m, n);
  L0->rows = L1->rows = m;
  L0->cols = L1->cols = n;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      set_d(&L0->entry[i][j], &A0->entry[i][j]);
      set_d(&L1->entry[i][j], &A1->entry[i][j]);
    }

  if (m == 0 || n == 0)
  { // return
    return 0;
  }

  // so we have m >= n >= 1
  init_vec_d(u, m);
  init_vec_d(z, m);
  init_mat_d(tempR0, n, n);
  init_mat_d(tempR1, n, n);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
  }

  // FIRST: do a QR decomposition for L0 & L1 (without pivoting)

  // do QR decomposition for top block
  L0->rows = L1->rows = rowTop;
  for (k = 0; k < colTop; k++)
  { // find the norm for this column
    norm[k] = 0;
    for (j = k; j < rowTop; j++)
      norm[k] += norm_sqr_d(&L0->entry[j][k]);
    norm[k] = sqrt(norm[k]);

    if (norm[k] > tol_sign)
    { // generate the Householder quantities
      genhh_d(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_d(z, L0, u, k, k+1, perm);
      // apply z to rows k to rowTop & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_d(L0, u, z, k, rowTop, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & L1
      findZ_mat_d(z, L1, u, k, k, perm);
      // apply z to rows k to rowTop & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_d(L1, u, z, k, rowTop, k, n, 0, u->size, 0, z->size, perm);
    }
  }
  // setup L0 & L1 back to full size
  L0->rows = L1->rows = m;

  // now we need to a normal QR decomposition on the bottom
  for (k = colTop; k < n; k++)
  { // find the norm for this column
    norm[k] = 0; 
    for (j = k; j < m; j++) 
      norm[k] += norm_sqr_d(&L0->entry[j][k]);
    norm[k] = sqrt(norm[k]);

    if (norm[k] > tol_sign)
    { // generate the Householder quantities
      genhh_d(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_d(z, L0, u, k, k+1, perm);
      // apply z to rows k to m & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_d(L0, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & L1
      findZ_mat_d(z, L1, u, k, k, perm);
      // apply z to rows k to m & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_d(L1, u, z, k, m, k, n, 0, u->size, 0, z->size, perm);
    }
  }

  // now that we have an upper triangular matrix, we remove the bottom 0's, transpose it and do a pivoted QR decomposition
  change_size_mat_d(tempR0, n, n);
  change_size_mat_d(tempR1, n, n);
  tempR0->rows = tempR0->cols = tempR1->rows = tempR1->cols = L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      conjugate_d(&tempR0->entry[j][i], &L0->entry[i][j]);
      conjugate_d(&tempR1->entry[j][i], &L1->entry[i][j]);
    }

  // SECOND: do a pivoted QR decomposition of tempR

  // setup the norms and find the first pivot
  max = pivot = 0;
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    norm[k] = 0;
    for (j = k; j < n; j++) // take advantage of lower triangular
      norm[k] += norm_sqr_d(&tempR0->entry[j][k]);
    // check to see if it is the maximum thus far
    if (max < norm[k])
    {
      max = norm[k];
      pivot = k;
    }
  }
  // take the sqrt for the pivot column
  prevNorm = norm[pivot] = sqrt(norm[pivot]);

  // do QR decomposition for tempR0
  for (k = 0; k < n; k++)
  { // make sure that the pivot is large enough
    if (norm[perm[pivot]] < tol_pivot || prevNorm > norm[perm[pivot]] * largeChange)
      break;
    else
      prevNorm = norm[perm[pivot]];

    // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_d(u, tempR0, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & tempR0
    findZ_mat_d(z, tempR0, u, k, k+1, perm);
    // apply z to rows k to n & cols k+1 to n of tempR0 using every entry in u & every entry in z
    apphh_Z_mat_d(tempR0, u, z, k, n, k+1, n, 0, u->size, 0, z->size, perm);

    // find the z associated to u & tempR1
    findZ_mat_d(z, tempR1, u, k, k, perm);
    // apply z to rows k to n & cols k to n of tempR1 using every entry in u & every entry in z
    apphh_Z_mat_d(tempR1, u, z, k, n, k, n, 0, u->size, 0, z->size, perm);

    // update the norms and find the one that is the maximum, if needed
    if (k + 1 < n)
    {
      pivot = k + 1;
      max = 0;
      for (j = k + 1; j < n; j++)
      {
        l = perm[j];
        norm[l] = 0;
        for (i = k + 1; i < n; i++)
          norm[l] += norm_sqr_d(&tempR0->entry[i][l]);
        if (max < norm[l])
        {
          max = norm[l];
          pivot = j;
        }
      }
      // take the sqrt for the pivot column
      norm[perm[pivot]] = sqrt(norm[perm[pivot]]);
    }
  }
  // so rank = k
  rank = k;

  // convert tempR0 to L0 & tempR1 to L1 - lower traingular matrix using the swaps in perm
  change_size_mat_d(L0, n, n);
  change_size_mat_d(L1, n, n);
  L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (j = 0; j < n; j++)
  { // setup column j
    for (i = 0; i < n; i++)
    {
      k = perm[j];
      conjugate_d(&L0->entry[j][i], &tempR0->entry[i][k]);
      conjugate_d(&L1->entry[j][i], &tempR1->entry[i][k]);
    }
  }

  // free the memory
  free(perm);
  free(norm);
  clear_vec_d(u); clear_vec_d(z);
  clear_mat_d(tempR0); clear_mat_d(tempR1);

  return (n - rank);
}

/////// multi precision ////////

int QLP_block_pair_mp_prec(mat_mp L0, mat_mp L1, mat_mp A0, mat_mp A1, int rowTop, int colTop, int curr_prec)
/***************************************************************\
* USAGE: runs QLP with the tolerances based on curr_prec        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES: Work with A0 and update A1 (i.e. A0 is more accurate)  *
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

  retVal = QLP_block_pair_mp(L0, L1, A0, A1, rowTop, colTop, tol_pivot, tol_sign, largeChange);

  mpf_clear(tol_pivot);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);
  free(str);

  return retVal;
}

int QLP_block_pair_mp(mat_mp L0, mat_mp L1, mat_mp A0, mat_mp A1, int rowTop, int colTop, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: Takes a pair of block triangular matrices and produces *
* a QLP^T decomposition where L is lower triangular matrix and  *
* Q & P would be unitary                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES: Work with A0 and update A1 (i.e. A0 is more accurate)  *
\***************************************************************/
{ 
  int i, j, k, l, rank, pivot, m = A0->rows, n = A0->cols, *perm = NULL;
  mpf_t tempMPF, prevNorm, max, *norm = NULL;
  vec_mp u, z;
  mat_mp tempR0, tempR1;

  // error checking
  if (A0->rows != A1->rows || A0->cols != A1->cols)
  {
    printf("ERROR: The matrices need to be the same size in QLP_block_pair_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QLP_block_pair_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QLP_block_pair_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= rowTop && rowTop <= m))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= colTop && colTop <= n))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  // copy A0 to L0 & A1 to L1
  change_size_mat_mp(L0, m, n); 
  change_size_mat_mp(L1, m, n);
  L0->rows = L1->rows = m;
  L0->cols = L1->cols = n;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      set_mp(&L0->entry[i][j], &A0->entry[i][j]);
      set_mp(&L1->entry[i][j], &A1->entry[i][j]);
    }

  if (m == 0 || n == 0)
  { // return
    return 0;
  }

  // so we have m >= n >= 1
  mpf_init(prevNorm); mpf_init(max); mpf_init(tempMPF);
  init_vec_mp(u, m);
  init_vec_mp(z, m);
  init_mat_mp(tempR0, n, n);
  init_mat_mp(tempR1, n, n);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    mpf_init(norm[k]);
  }

  // FIRST: do a QR decomposition for L0 & L1 (without pivoting)

  // do QR decomposition for top block
  L0->rows = L1->rows = rowTop;
  for (k = 0; k < colTop; k++)
  { // find the norm for this column
    mpf_set_ui(norm[k], 0);
    for (j = k; j < rowTop; j++)
    {
      norm_sqr_mp(tempMPF, &L0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    mpf_sqrt(norm[k], norm[k]);

    if (mpf_cmp(norm[k], tol_sign) > 0)
    { // generate the Householder quantities
      genhh_mp(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_mp(z, L0, u, k, k+1, perm);
      // apply z to rows k to rowTop & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_mp(L0, u, z, k, rowTop, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & L1
      findZ_mat_mp(z, L1, u, k, k, perm);
      // apply z to rows k to rowTop & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_mp(L1, u, z, k, rowTop, k, n, 0, u->size, 0, z->size, perm);
    }
  }
  // setup L0 & L1 back to full size
  L0->rows = L1->rows = m;

  // now we need to a normal QR decomposition on the bottom
  for (k = colTop; k < n; k++)
  { // find the norm for this column
    mpf_set_ui(norm[k], 0);
    for (j = k; j < m; j++) 
    {
      norm_sqr_mp(tempMPF, &L0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    mpf_sqrt(norm[k], norm[k]);

    if (mpf_cmp(norm[k], tol_sign) > 0)
    { // generate the Householder quantities
      genhh_mp(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_mp(z, L0, u, k, k+1, perm);
      // apply z to rows k to m & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_mp(L0, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & L1
      findZ_mat_mp(z, L1, u, k, k, perm);
      // apply z to rows k to m & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_mp(L1, u, z, k, m, k, n, 0, u->size, 0, z->size, perm);
    }
  }

  // now that we have an upper triangular matrix, we remove the bottom 0's, transpose it and do a pivoted QR decomposition
  change_size_mat_mp(tempR0, n, n);
  change_size_mat_mp(tempR1, n, n);
  tempR0->rows = tempR0->cols = tempR1->rows = tempR1->cols = L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      conjugate_mp(&tempR0->entry[j][i], &L0->entry[i][j]);
      conjugate_mp(&tempR1->entry[j][i], &L1->entry[i][j]);
    }

  // SECOND: do a pivoted QR decomposition of tempR

  // setup the norms and find the first pivot
  pivot = 0;
  mpf_set_ui(max, 0);
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    mpf_set_ui(norm[k], 0);
    for (j = k; j < n; j++) // take advantage of lower triangular
    {
      norm_sqr_mp(tempMPF, &tempR0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    // check to see if it is the maximum thus far
    if (mpf_cmp(max, norm[k]) < 0)
    {
      mpf_set(max, norm[k]);
      pivot = k;
    }
  }
  // take the sqrt for the pivot column
  mpf_sqrt(norm[pivot], norm[pivot]);
  mpf_set(prevNorm, norm[pivot]);

  // do QR decomposition for tempR0
  for (k = 0; k < n; k++)
  { // make sure that the pivot is large enough
    if (checkGood_mp(norm[perm[pivot]], prevNorm, tol_pivot, largeChange))
      mpf_set(prevNorm, norm[perm[pivot]]);
    else
      break;

    // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_mp(u, tempR0, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & tempR0
    findZ_mat_mp(z, tempR0, u, k, k+1, perm);
    // apply z to rows k to n & cols k+1 to n of tempR0 using every entry in u & every entry in z
    apphh_Z_mat_mp(tempR0, u, z, k, n, k+1, n, 0, u->size, 0, z->size, perm);

    // find the z associated to u & tempR1
    findZ_mat_mp(z, tempR1, u, k, k, perm);
    // apply z to rows k to n & cols k to n of tempR1 using every entry in u & every entry in z
    apphh_Z_mat_mp(tempR1, u, z, k, n, k, n, 0, u->size, 0, z->size, perm);

    // update the norms and find the one that is the maximum, if needed
    if (k + 1 < n)
    {
      pivot = k + 1;
      mpf_set_ui(max, 0);
      for (j = k + 1; j < n; j++)
      {
        l = perm[j];
        mpf_set_ui(norm[l], 0);
        for (i = k + 1; i < n; i++)
        {
          norm_sqr_mp(tempMPF, &tempR0->entry[i][l]);
          mpf_add(norm[l], norm[l], tempMPF);
        }
        if (mpf_cmp(max, norm[l]) < 0)
        {
          mpf_set(max, norm[l]);
          pivot = j;
        }
      }
      // take the sqrt for the pivot column
      mpf_sqrt(norm[perm[pivot]], norm[perm[pivot]]);
    }
  }
  // so rank = k
  rank = k;

  // convert tempR0 to L0 & tempR1 to L1 - lower traingular matrix using the swaps in perm
  change_size_mat_mp(L0, n, n);
  change_size_mat_mp(L1, n, n);
  L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (j = 0; j < n; j++)
  { // setup column j
    for (i = 0; i < n; i++)
    {
      k = perm[j];
      conjugate_mp(&L0->entry[j][i], &tempR0->entry[i][k]);
      conjugate_mp(&L1->entry[j][i], &tempR1->entry[i][k]);
    }
  }

  // free the memory
  mpf_clear(tempMPF); mpf_clear(prevNorm); mpf_clear(max);
  for (j = n - 1; j >= 0; j--)
    mpf_clear(norm[j]);
  free(perm); free(norm);
  clear_vec_mp(u); clear_vec_mp(z);
  clear_mat_mp(tempR0); clear_mat_mp(tempR1);

  return (n - rank);
}

int QLP_block_pair_amp_mp_d_prec(mat_mp L0, mat_d L1, mat_mp A0, int prec0, mat_d A1, int rowTop, int colTop)
/***************************************************************\
* USAGE: runs QLP with the tolerances based on curr_prec        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES: Work with A0 and update A1 (i.e. A0 is more accurate)  *
\***************************************************************/
{
  int retVal, num_digits = prec_to_digits(prec0);
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

  retVal = QLP_block_pair_amp_mp_d(L0, L1, A0, prec0, A1, rowTop, colTop, tol_pivot, tol_sign, largeChange);

  mpf_clear(tol_pivot);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);
  free(str);

  return retVal;
}

int QLP_block_pair_amp_mp_d(mat_mp L0, mat_d L1, mat_mp A0, int prec0, mat_d A1, int rowTop, int colTop, mpf_t tol_pivot, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: Takes a pair of block triangular matrices and produces *
* a QLP^T decomposition where L is lower triangular matrix and  *
* Q & P would be unitary                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank of A - i.e. rank deficiency *
* NOTES: Work with A0 and update A1 (i.e. A0 is more accurate)  *
\***************************************************************/
{ 
  int i, j, k, l, rank, pivot, m = A0->rows, n = A0->cols, *perm = NULL;
  mpf_t tempMPF, prevNorm, max, *norm = NULL;
  vec_d u_d, z_d;
  vec_mp u, z;
  mat_d tempR1;
  mat_mp tempR0;

  // error checking
  if (A0->rows != A1->rows || A0->cols != A1->cols)
  {
    printf("ERROR: The matrices need to be the same size in QLP_block_pair_amp_mp_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QLP_block_pair_amp_mp_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QLP_block_pair_amp_mp_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= rowTop && rowTop <= m))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_amp_mp_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (!(0 <= colTop && colTop <= n))
  {
    printf("ERROR: The top block is not setup properly in QLP_block_pair_amp_mp_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  // copy A0 to L0 & A1 to L1
  change_size_mat_mp(L0, m, n); 
  change_size_mat_d(L1, m, n);
  L0->rows = L1->rows = m;
  L0->cols = L1->cols = n;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      set_mp(&L0->entry[i][j], &A0->entry[i][j]);
      set_d(&L1->entry[i][j], &A1->entry[i][j]);
    }

  if (m == 0 || n == 0)
  { // return
    return 0;
  }

  // so we have m >= n >= 1
  mpf_init2(prevNorm, prec0); mpf_init2(max, prec0); mpf_init2(tempMPF, prec0);
  init_vec_d(u_d, m); init_vec_d(z_d, m);
  init_vec_mp2(u, m, prec0); init_vec_mp2(z, m, prec0);
  init_mat_mp2(tempR0, n, n, prec0);
  init_mat_d(tempR1, n, n);

  // setup perm & norm
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    mpf_init2(norm[k], prec0);
  }

  // FIRST: do a QR decomposition for L0 & L1 (without pivoting)

  // do QR decomposition for top block
  L0->rows = L1->rows = rowTop;
  for (k = 0; k < colTop; k++)
  { // find the norm for this column
    mpf_set_ui(norm[k], 0);
    for (j = k; j < rowTop; j++)
    {
      norm_sqr_mp(tempMPF, &L0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    mpf_sqrt(norm[k], norm[k]);

    if (mpf_cmp(norm[k], tol_sign) > 0)
    { // generate the Householder quantities
      genhh_mp(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_mp(z, L0, u, k, k+1, perm);
      // apply z to rows k to rowTop & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_mp(L0, u, z, k, rowTop, k+1, n, 0, u->size, 0, z->size, perm);

      // convert u to u_d
      vec_mp_to_d(u_d, u);
      // find the z associated to u & L1
      findZ_mat_d(z_d, L1, u_d, k, k, perm);
      // apply z to rows k to rowTop & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_d(L1, u_d, z_d, k, rowTop, k, n, 0, u_d->size, 0, z_d->size, perm);
    }
  }
  // setup L0 & L1 back to full size
  L0->rows = L1->rows = m;

  // now we need to a normal QR decomposition on the bottom
  for (k = colTop; k < n; k++)
  { // find the norm for this column
    mpf_set_ui(norm[k], 0);
    for (j = k; j < m; j++) 
    {
      norm_sqr_mp(tempMPF, &L0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    mpf_sqrt(norm[k], norm[k]);

    if (mpf_cmp(norm[k], tol_sign) > 0)
    { // generate the Householder quantities
      genhh_mp(u, L0, k, k, norm[k], perm, tol_sign);
      // find the z associated to u & L0
      findZ_mat_mp(z, L0, u, k, k+1, perm);
      // apply z to rows k to m & cols k+1 to n of L0 using every entry in u & every entry in z
      apphh_Z_mat_mp(L0, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

      // convert u to u_d
      vec_mp_to_d(u_d, u);
      // find the z associated to u & L1
      findZ_mat_d(z_d, L1, u_d, k, k, perm);
      // apply z to rows k to m & cols k to n of L1 using every entry in u & every entry in z
      apphh_Z_mat_d(L1, u_d, z_d, k, m, k, n, 0, u_d->size, 0, z_d->size, perm);
    }
  }

  // now that we have an upper triangular matrix, we remove the bottom 0's, transpose it and do a pivoted QR decomposition
  change_size_mat_mp(tempR0, n, n);
  change_size_mat_d(tempR1, n, n);
  tempR0->rows = tempR0->cols = tempR1->rows = tempR1->cols = L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      conjugate_mp(&tempR0->entry[j][i], &L0->entry[i][j]);
      conjugate_d(&tempR1->entry[j][i], &L1->entry[i][j]);
    }

  // SECOND: do a pivoted QR decomposition of tempR

  // setup the norms and find the first pivot
  pivot = 0;
  mpf_set_ui(max, 0);
  for (k = 0; k < n; k++)
  {
    perm[k] = k;
    mpf_set_ui(norm[k], 0);
    for (j = k; j < n; j++) // take advantage of lower triangular
    {
      norm_sqr_mp(tempMPF, &tempR0->entry[j][k]);
      mpf_add(norm[k], norm[k], tempMPF);
    }
    // check to see if it is the maximum thus far
    if (mpf_cmp(max, norm[k]) < 0)
    {
      mpf_set(max, norm[k]);
      pivot = k;
    }
  }
  // take the sqrt for the pivot column
  mpf_sqrt(norm[pivot], norm[pivot]);
  mpf_set(prevNorm, norm[pivot]);

  // do QR decomposition for tempR0
  for (k = 0; k < n; k++)
  { // make sure that the pivot is large enough
    if (checkGood_mp(norm[perm[pivot]], prevNorm, tol_pivot, largeChange))
      mpf_set(prevNorm, norm[perm[pivot]]);
    else
      break;

    // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // generate the Householder quantities
    genhh_mp(u, tempR0, k, k, norm[perm[k]], perm, tol_sign);
    // find the z associated to u & tempR0
    findZ_mat_mp(z, tempR0, u, k, k+1, perm);
    // apply z to rows k to n & cols k+1 to n of tempR0 using every entry in u & every entry in z
    apphh_Z_mat_mp(tempR0, u, z, k, n, k+1, n, 0, u->size, 0, z->size, perm);

    // convert u to u_d
    vec_mp_to_d(u_d, u);
    // find the z associated to u & tempR1
    findZ_mat_d(z_d, tempR1, u_d, k, k, perm);
    // apply z to rows k to n & cols k to n of tempR1 using every entry in u & every entry in z
    apphh_Z_mat_d(tempR1, u_d, z_d, k, n, k, n, 0, u_d->size, 0, z_d->size, perm);

    // update the norms and find the one that is the maximum, if needed
    if (k + 1 < n)
    {
      pivot = k + 1;
      mpf_set_ui(max, 0);
      for (j = k + 1; j < n; j++)
      {
        l = perm[j];
        mpf_set_ui(norm[l], 0);
        for (i = k + 1; i < n; i++)
        {
          norm_sqr_mp(tempMPF, &tempR0->entry[i][l]);
          mpf_add(norm[l], norm[l], tempMPF);
        }
        if (mpf_cmp(max, norm[l]) < 0)
        {
          mpf_set(max, norm[l]);
          pivot = j;
        }
      }
      // take the sqrt for the pivot column
      mpf_sqrt(norm[perm[pivot]], norm[perm[pivot]]);
    }
  }
  // so rank = k
  rank = k;

  // convert tempR0 to L0 & tempR1 to L1 - lower traingular matrix using the swaps in perm
  change_size_mat_mp(L0, n, n);
  change_size_mat_d(L1, n, n);
  L0->rows = L0->cols = L1->rows = L1->cols = n;
  for (j = 0; j < n; j++)
  { // setup column j
    for (i = 0; i < n; i++)
    {
      k = perm[j];
      conjugate_mp(&L0->entry[j][i], &tempR0->entry[i][k]);
      conjugate_d(&L1->entry[j][i], &tempR1->entry[i][k]);
    }
  }

  // free the memory
  mpf_clear(tempMPF); mpf_clear(prevNorm); mpf_clear(max);
  for (j = n - 1; j >= 0; j--)
    mpf_clear(norm[j]);
  free(perm); free(norm);
  clear_vec_d(u_d); clear_vec_d(z_d);
  clear_vec_mp(u); clear_vec_mp(z);
  clear_mat_mp(tempR0); 
  clear_mat_d(tempR1);

  return (n - rank);
}


// QR decomposition for a pair of matrices

void findZ_mat_left_d(vec_d z, mat_d C, vec_d u, int sr, int sc, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: z = - C(sr:m, sc:n) * u                                *
\***************************************************************/
{
  int i, j, col, row, size, cols = u->size;

  if (u->size != C->cols - sc)
  {
    printf("The sizes of matrices do not match in findZ_mat_left_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  size = C->rows - sr;
  change_size_vec_d(z, size);
  z->size = size;

  for (j = 0; j < size; j++)
  {
    row = j;
    set_zero_d(&z->coord[j]);
    for (i = 0; i < cols; i++)
    {
      col = colnum[sc+i];
      sum_mul_d(&z->coord[j], &C->entry[row][col], &u->coord[i]);
    }
    neg_d(&z->coord[j], &z->coord[j]);
  }

  return;
}

void apphh_Z_mat_left_d(mat_d C, vec_d u, vec_d z, int sr, int er, int sc, int ec, int u_sc, int u_ec, int z_sc, int z_ec, int *colnum)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: C = C + z*uH                                           *
\***************************************************************/
{
  int i, j, zcoord, col, row, rows = er - sr, cols = ec - sc;

  if ((er - sr != z_ec - z_sc) || (u_ec - u_sc != ec - sc))
  {
    printf("The sizes of matrices do not match in apphh_Z_mat_left_d!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  for (j = 0; j < cols; j++, sc++, u_sc++)
  { // update jth column of C
    row = sr;
    col = colnum[sc];
    zcoord = z_sc;
    for (i = 0; i < rows; i++, row++, zcoord++)
    { // update C[i][j] -= conj(u[j])*z[i]
      sum_conj_mul_d(&C->entry[row][col], &u->coord[u_sc], &z->coord[zcoord]);
    }
  }

  return;
}

void QR_pair_d(mat_d Q, mat_d R0, mat_d R1, mat_d P, mat_d A0, mat_d A1, double tol_pivot, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: Takes a block triangular matrix A and turns it into a  *
* upper triangular matrix R using Householder transformations   *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computes Q, R & P                                      *
\***************************************************************/
{
  int i, j, k, l, pivot, m = A0->rows, n = A0->cols, *perm = NULL, *id_perm = NULL;
  double prevNorm, max, tempD, *norm = NULL;
  comp_d tempComp, tempComp_conj, tempComp2;
  vec_d u, z;
  mat_d tempR0, tempR1;

  // error checking
  if (A0->rows != A1->rows || A0->cols != A1->cols)
  {
    printf("ERROR: The matrices need to be the same size in QR_R_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_R_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QR_R_pair_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  if (m == 0 || n == 0)
  { // copy A to R and return
    mat_cp_d(R0, A0);
    mat_cp_d(R1, A1);
    return;
  }

  // we have m >= n >= 1

  // initialize
  init_vec_d(u, m); init_vec_d(z, m);
  init_mat_d(tempR0, m, n); init_mat_d(tempR1, m, n);
  perm = (int *)bmalloc(n * sizeof(int));
  id_perm = (int *)bmalloc(m * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));

  // setup id_perm
  for (i = 0; i < m; i++)
    id_perm[i] = i;

  // setup tempR0 & tempR1
  max = pivot = 0;
  tempR0->rows = tempR1->rows = m;
  tempR0->cols = tempR1->cols = n;
  for (j = 0; j < n; j++)
  { // setup perm, id_perm & norm
    perm[j] = j;
    norm[j] = 0;

    for (i = 0; i < m; i++)
    { // setup tempR0 & tempR1
      set_d(&tempR0->entry[i][j], &A0->entry[i][j]);
      norm[j] += norm_sqr_d(&tempR0->entry[i][j]);
      set_d(&tempR1->entry[i][j], &A1->entry[i][j]);
    }
    norm[j] = sqrt(norm[j]);
    if (norm[j] > max)
    {
      max = norm[j];
      pivot = j;
    }
  }
  prevNorm = max;

  // set Q to be identity matrix
  make_matrix_ID_d(Q, m, m);

  // do a QR decomposition with pivoting
  for (k = 0; k < n; k++)
  { // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // make sure the pivot is large enough
    if (norm[perm[k]] > tol_pivot && prevNorm < norm[perm[k]] * largeChange)
    { // update prevNorm
      prevNorm = norm[perm[k]];

      // generate the Householder quantities
      genhh_d(u, tempR0, k, k, norm[perm[k]], perm, tol_sign);
      // find the z associated to u & tempR0
      findZ_mat_d(z, tempR0, u, k, k+1, perm);
      // apply z to rows k to m & cols k+1 to n of tempR0
      apphh_Z_mat_d(tempR0, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & tempR1
      findZ_mat_d(z, tempR1, u, k, k, perm);
      // apply z to rows k to m & cols k to n of tempR1
      apphh_Z_mat_d(tempR1, u, z, k, m, k, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & Q
      findZ_mat_left_d(z, Q, u, 0, k, id_perm);
      // apply z to rows 0 to m & cols k to m of Q
      apphh_Z_mat_left_d(Q, u, z, 0, m, k, m, 0, u->size, 0, z->size, id_perm);
    }

    // update the norms and find the one that is the maximum, if needed
    pivot = k+1;
    max = 0;
    for (j = k+1; j < n; j++)
    { // determine the column
      l = perm[j];
      norm[l] = 0;
      for (i = k+1; i < m; i++)
        norm[l] += norm_sqr_d(&tempR0->entry[i][l]);
      norm[l] = sqrt(norm[l]);
      if (max < norm[l])
      {
         max = norm[l];
         pivot = j;
      }
    }
  }

  // convert perm to P
  convertToP_d(P, perm, n);

  // setup R0 & R1
  change_size_mat_d(R0, m, n);
  change_size_mat_d(R1, m, n);
  R0->rows = R1->rows = m;
  R0->cols = R1->cols = n;
  for (j = 0; j < n; j++)
  { // setup col j
    for (i = 0; i < m; i++)
    {
      set_d(&R0->entry[i][j], &tempR0->entry[i][perm[j]]);
      set_d(&R1->entry[i][j], &tempR1->entry[i][perm[j]]);
    }

    // make the diagonal entries of R0 real & positive
    tempD = sign_d(tempComp, &R0->entry[j][j], tol_sign);
    set_d(tempComp_conj, tempComp);
    conjugate_d(tempComp, tempComp);
    set_double_d(&R0->entry[j][j], tempD, 0);

    tempD = sign_d(tempComp2, &R1->entry[j][j], tol_sign);
    conjugate_d(tempComp2, tempComp2);
    set_double_d(&R1->entry[j][j], tempD, 0);

    for (i = 0; i < m; i++)
    {
      mul_d(&Q->entry[i][j], &Q->entry[i][j], tempComp_conj);
      if (i > j && i < n)
      {
        mul_d(&tempR0->entry[j][perm[i]], &tempR0->entry[j][perm[i]], tempComp);
        mul_d(&tempR1->entry[j][perm[i]], &tempR1->entry[j][perm[i]], tempComp2);
      }
    }
  }

  // clear memory
  clear_vec_d(u); clear_vec_d(u);
  clear_mat_d(tempR0); clear_mat_d(tempR1);
  free(perm); free(id_perm); free(norm);

  return;
}

void QR_R_pair_d(mat_d R0, mat_d R1, mat_d P, mat_d A0, mat_d A1, double tol_pivot, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: Takes a block triangular matrix A and turns it into a  *
* upper triangular matrix R using Householder transformations   *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computes R & P (does not compute Q)                    *
\***************************************************************/
{ 
  int i, j, k, l, pivot, m = A0->rows, n = A0->cols, *perm = NULL;
  double prevNorm, max, tempD, *norm = NULL;
  comp_d tempComp, tempComp2;
  vec_d u, z;
  mat_d tempR0, tempR1;

  // error checking
  if (A0->rows != A1->rows || A0->cols != A1->cols)
  {
    printf("ERROR: The matrices need to be the same size in QR_R_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  } 
  else if (m < n)
  {
    printf("ERROR: The matrix has more columns than rows in QR_R_pair_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in QR_R_pair_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m >= n

  if (m == 0 || n == 0)
  { // copy A to R and return
    mat_cp_d(R0, A0);
    mat_cp_d(R1, A1);
    return;
  }

  // we have m >= n >= 1

  // initialize
  init_vec_d(u, m); init_vec_d(z, m);
  init_mat_d(tempR0, m, n); init_mat_d(tempR1, m, n);
  perm = (int *)bmalloc(n * sizeof(int));
  norm = (double *)bmalloc(n * sizeof(double));
  
  // setup tempR0 & tempR1
  max = pivot = 0;
  tempR0->rows = tempR1->rows = m;
  tempR0->cols = tempR1->cols = n;
  for (j = 0; j < n; j++)
  { // setup perm & norm
    perm[j] = j;
    norm[j] = 0;

    for (i = 0; i < m; i++)
    { // setup tempR0 & tempR1
      set_d(&tempR0->entry[i][j], &A0->entry[i][j]);
      norm[j] += norm_sqr_d(&tempR0->entry[i][j]);
      set_d(&tempR1->entry[i][j], &A1->entry[i][j]);
    }
    norm[j] = sqrt(norm[j]);
    if (norm[j] > max)
    {
      max = norm[j];
      pivot = j;
    }
  }
  prevNorm = max;

  // do a QR decomposition with pivoting
  for (k = 0; k < n; k++)
  { // swap the columns if needed
    if (pivot != k)
    {
      j = perm[pivot];
      perm[pivot] = perm[k];
      perm[k] = j;
    }

    // make sure the pivot is large enough
    if (norm[perm[k]] > tol_pivot && prevNorm < norm[perm[k]] * largeChange)
    { // update prevNorm
      prevNorm = norm[perm[k]];

      // generate the Householder quantities
      genhh_d(u, tempR0, k, k, norm[perm[k]], perm, tol_sign);
      // find the z associated to u & tempR0
      findZ_mat_d(z, tempR0, u, k, k+1, perm);
      // apply z to rows k to m & cols k+1 to n of tempR0
      apphh_Z_mat_d(tempR0, u, z, k, m, k+1, n, 0, u->size, 0, z->size, perm);

      // find the z associated to u & tempR1
      findZ_mat_d(z, tempR1, u, k, k, perm);
      // apply z to rows k to m & cols k to n of tempR1
      apphh_Z_mat_d(tempR1, u, z, k, m, k, n, 0, u->size, 0, z->size, perm);
    }

    // update the norms and find the one that is the maximum, if needed
    pivot = k+1;
    max = 0;
    for (j = k+1; j < n; j++)
    { // determine the column
      l = perm[j];
      norm[l] = 0;
      for (i = k+1; i < m; i++)
        norm[l] += norm_sqr_d(&tempR0->entry[i][l]);
      norm[l] = sqrt(norm[l]);
      if (max < norm[l])
      {
         max = norm[l];
         pivot = j;
      }
    }
  }

  // convert perm to P
  convertToP_d(P, perm, n);

  // setup R0 & R1
  change_size_mat_d(R0, m, n);
  change_size_mat_d(R1, m, n);
  R0->rows = R1->rows = m;
  R0->cols = R1->cols = n;
  for (j = 0; j < n; j++)
  { // setup col j
    for (i = 0; i < m; i++)
    {
      set_d(&R0->entry[i][j], &tempR0->entry[i][perm[j]]);
      set_d(&R1->entry[i][j], &tempR1->entry[i][perm[j]]);
    }

    // make the diagonal entries of R0 real & positive
    tempD = sign_d(tempComp, &R0->entry[j][j], tol_sign);
    conjugate_d(tempComp, tempComp);
    set_double_d(&R0->entry[j][j], tempD, 0);

    tempD = sign_d(tempComp2, &R1->entry[j][j], tol_sign);
    conjugate_d(tempComp2, tempComp2);
    set_double_d(&R1->entry[j][j], tempD, 0);

    for (i = j+1; i < n; i++)
    {
      mul_d(&tempR0->entry[j][perm[i]], &tempR0->entry[j][perm[i]], tempComp);
      mul_d(&tempR1->entry[j][perm[i]], &tempR1->entry[j][perm[i]], tempComp2);
    }
  }

  // clear memory
  clear_vec_d(u); clear_vec_d(u);
  clear_mat_d(tempR0); clear_mat_d(tempR1);
  free(perm); free(norm);

  return;
}
  

