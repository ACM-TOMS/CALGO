// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "localdim.h"

// take in 2 similar matrices and determine the rank (number of linealy independent columns) efficiently
void QR_houseHolder_d(vec_d *houseHolderVecs, mat_d A, double tol_sign);
void apply_HH_d(mat_d A, int num, vec_d *HHVecs);
void apply_HH_pair_d(mat_d A1, mat_d A2, int num, vec_d *HHVecs);

int rank_MM_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_d MM1, mat_d MM2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of the new null space                     *
* NOTES: Assume MM1 is more accurate                            *
\***************************************************************/
{
  int i, j, k, l, r, c, nullSize = 0, m = MM1->rows, n = MM1->cols;
  double tol_sign = 1e-20;
  mat_d tempR1, tempR2, randMat;

  // error checking
  if (MM1->rows != MM2->rows || MM1->cols != MM2->cols)
  {
    printf("ERROR: The matrices need to be the same size in rank_MM_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m == 0 || n == 0)
  {
    printf("ERROR: The matrices need to have both columns and rows in rank_MM_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (order < 1)
  {
    printf("ERROR: The order of the multiplicity matrix must be atleast 1 in rank_MM_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize memory
  init_mat_d(tempR1, 0, 0);
  init_mat_d(tempR2, 0, 0);
  init_mat_d(randMat, 0, 0);

  // determine if we have a multiplicity matrix of order == 1, == 2 or > 2
  if (order == 1)
  { // we compute the size of the null space normally
    nullSize = corank_rrv_d(CN, smallest_nonzero, largest_zero, MM1, MM2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }
  else if (order == 2)
  { // we randomize the first block of columns to the correct rank and then do a QR decomposition on that block
    change_size_mat_d(randMat, colSize[0], ranks[0]);
    make_matrix_random_d(randMat, colSize[0], ranks[0]);

    // setup tempR1 to be the randomized columns of the first block of MM1
    r = rowSize[0];
    c = ranks[0];
    change_size_mat_d(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = colSize[0];
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_d(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_d(&tempR1->entry[i][j], &MM1->entry[i][k], &randMat->entry[k][j]);
        }
      } 

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs = ranks[0];
    *houseHolderVecs = (vec_d *)bmalloc(*numHHVecs * sizeof(vec_d));
    QR_houseHolder_d(*houseHolderVecs, tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[0];
    change_size_mat_d(tempR1, r, c);
    change_size_mat_d(tempR2, r, c);
    tempR1->rows = tempR2->rows = r;
    tempR1->cols = tempR2->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[0] + j;
      for (i = 0; i < r; i++)
      {
        set_d(&tempR1->entry[i][j], &MM1->entry[i][l]);
        set_d(&tempR2->entry[i][j], &MM2->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_d(tempR1, tempR2, *numHHVecs, *houseHolderVecs); 

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_d(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_d(&tempR2->entry[i][j], &tempR2->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2->rows = r;

    // find the size of the new null space 
    nullSize = corank_rrv_d(CN, smallest_nonzero, largest_zero, tempR1, tempR2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }
  else // order >= 3
  { // setup tempR1 to be last known block of MM1

    r = rowSize[order - 2];
    c = colSize[order - 2] - colSize[order - 3];
    change_size_mat_d(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;
    for (j = 0; j < c; j++)
    {
      l = j + colSize[order - 3];
      for (i = 0; i < r; i++)
      {
        set_d(&tempR1->entry[i][j], &MM1->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempR1
    apply_HH_d(tempR1, *numHHVecs, *houseHolderVecs);

    // remove the top part of tempR1
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    change_size_mat_d(tempR2, r, c);
    tempR2->rows = r;
    tempR2->cols = c;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_d(&tempR2->entry[i][j], &tempR1->entry[l][j]);
      }
    }

    // setup a random matrix to randomize to the correct rank
    change_size_mat_d(randMat, c, ranks[order - 2]);
    make_matrix_random_d(randMat, c, ranks[order - 2]);

    // setup tempR1 to be the randomized columns of the last known block of MM1
    r = tempR2->rows;
    c = ranks[order - 2];
    change_size_mat_d(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = randMat->rows;
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_d(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_d(&tempR1->entry[i][j], &tempR2->entry[i][k], &randMat->entry[k][j]);
        }
      }

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs += c;
    *houseHolderVecs = (vec_d *)brealloc(*houseHolderVecs, *numHHVecs * sizeof(vec_d));
    QR_houseHolder_d(&(*houseHolderVecs)[*numHHVecs - c], tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[order - 2];
    change_size_mat_d(tempR1, r, c);
    change_size_mat_d(tempR2, r, c);
    tempR1->rows = tempR2->rows = r;
    tempR1->cols = tempR2->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[order - 2] + j;
      for (i = 0; i < r; i++)
      {
        set_d(&tempR1->entry[i][j], &MM1->entry[i][l]);
        set_d(&tempR2->entry[i][j], &MM2->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_d(tempR1, tempR2, *numHHVecs, *houseHolderVecs);

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_d(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_d(&tempR2->entry[i][j], &tempR2->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2->rows = r;

    // find the size of the new null space
    nullSize = corank_rrv_d(CN, smallest_nonzero, largest_zero, tempR1, tempR2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }

  // clear memory
  clear_mat_d(tempR1);
  clear_mat_d(tempR2);
  clear_mat_d(randMat);

  return nullSize;
}

// compute a QR decomposition on A without computing Q and without pivoting, but does store the Householder vectors
void QR_houseHolder_d(vec_d *houseHolderVecs, mat_d A, double tol_sign)
/***************************************************************\
* USAGE: finds A = Q * A - without finding Q explicitly         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume A has full column rank                          *
\***************************************************************/
{
  int i, k, m = A->rows, n = A->cols, *perm = (int *)bmalloc(n * sizeof(int));
  double norm;
  vec_d z;

  // initialize z
  init_vec_d(z, m);

  // setup perm to be the identity
  for (i = 0; i < n; i++)
    perm[i] = i;

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  { // compute the norm of the kth column
    norm = 0;
    for (i = k; i < m; i++)
      norm += norm_sqr_d(&A->entry[i][k]);
    norm = sqrt(norm);

    // generate the Householder vector
    init_vec_d(houseHolderVecs[k], 0);
    genhh_d(houseHolderVecs[k], A, k, k, norm, perm, tol_sign);
    // find the z associated to the Householder vector & A
    findZ_mat_d(z, A, houseHolderVecs[k], k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in the Householder vector & every entry in z
    apphh_Z_mat_d(A, houseHolderVecs[k], z, k, m, k+1, n, 0, houseHolderVecs[k]->size, 0, z->size, perm);
  }

  // clear memory
  free(perm);
  clear_vec_d(z);

  return;
}

// apply Householder vectors to a matrices
void apply_HH_d(mat_d A, int num, vec_d *HHVecs)
/***************************************************************\
* USAGE: apply the Householder vectors to A                     *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int k, m = A->rows, n = A->cols, *perm = (int *)bmalloc(n * sizeof(int));
  vec_d z;

  // initialize z
  init_vec_d(z, m);

  // setup perm to be the identity
  for (k = 0; k < n; k++)
    perm[k] = k;

  // perform the Householder updates
  for (k = 0; k < num; k++)
  { // set the sizes that are needed for this Householder vector
    A->rows = HHVecs[k]->size + k;

    // find the z associated to HHVecs[k] & A
    findZ_mat_d(z, A, HHVecs[k], k, 0, perm);
    // apply z to rows k to rows & cols 0 to cols of A using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_d(A, HHVecs[k], z, k, A->rows, 0, A->cols, 0, HHVecs[k]->size, 0, z->size, perm);
  }

  // set the sizes back
  A->rows = m;
  A->cols = n;

  // clear memory
  free(perm);
  clear_vec_d(z);

  return;
}

// apply Householder vectors to a pair of matrices
void apply_HH_pair_d(mat_d A1, mat_d A2, int num, vec_d *HHVecs)
/***************************************************************\
* USAGE: apply the Householder vectors to A1 & A2               *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int k, m = A1->rows, n = A1->cols, *perm = (int *)bmalloc(n * sizeof(int));
  vec_d z;

  // initialize z
  init_vec_d(z, m);

  // setup perm to be the identity
  for (k = 0; k < n; k++)
    perm[k] = k;

  // perform the Householder updates
  for (k = 0; k < num; k++)
  { // set the sizes that are needed for this Householder vector
    A1->rows = A2->rows = HHVecs[k]->size + k;

    // find the z associated to HHVecs[k] & A1
    findZ_mat_d(z, A1, HHVecs[k], k, 0, perm);
    // apply z to rows k to rows & cols 0 to cols of A1 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_d(A1, HHVecs[k], z, k, A1->rows, 0, A1->cols, 0, HHVecs[k]->size, 0, z->size, perm);

    // find the z associated to HHVecs[k] & A2
    findZ_mat_d(z, A2, HHVecs[k], k, 0, perm);
    // apply z to rows k to m & cols 0 to n of A2 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_d(A2, HHVecs[k], z, k, A2->rows, 0, A2->cols, 0, HHVecs[k]->size, 0, z->size, perm);
  }

  // set the sizes back
  A1->rows = A2->rows = m;
  A1->cols = A2->cols = n;

  // clear memory
  free(perm);
  clear_vec_d(z);

  return;
}

/////// MULTI PRECISION ///////////

// take in 2 similar matrices and determine the rank (number of linealy independent columns) efficiently
void QR_houseHolder_mp(vec_mp *houseHolderVecs, mat_mp A, mpf_t tol_sign);
void apply_HH_mp(mat_mp A, int num, vec_mp *HHVecs);
void apply_HH_pair_mp(mat_mp A1, mat_mp A2, int num, vec_mp *HHVecs);

int rank_MM_mp(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp MM1, mat_mp MM2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_mp **houseHolderVecs, double max_CN, double max_SV_ratio, double SV_tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of the new null space                     *
* NOTES: Assume MM1 is more accurate                            *
\***************************************************************/
{
  int i, j, k, l, r, c, prec_digits, nullSize = 0, m = MM1->rows, n = MM1->cols, curr_prec = mpf_get_default_prec();
  size_t size;
  char *str = NULL;
  mpf_t tol_sign;
  mat_mp tempR1, tempR2, randMat;

  // setup prec_digits
  prec_digits = prec_to_digits(curr_prec) - 3;

  // setup tol_sign
  mpf_init(tol_sign);
  size = 1 + snprintf(NULL, 0, "1e-%d", 2 * prec_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", 2 * prec_digits);
  mpf_set_str(tol_sign, str, 10);  
  
  // error checking
  if (MM1->rows != MM2->rows || MM1->cols != MM2->cols)
  {
    printf("ERROR: The matrices need to be the same size in rank_MM_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m == 0 || n == 0)
  {
    printf("ERROR: The matrices need to have both columns and rows in rank_MM_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (order < 1)
  {
    printf("ERROR: The order of the multiplicity matrix must be atleast 1 in rank_MM_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize memory
  init_mat_mp(tempR1, 0, 0);
  init_mat_mp(tempR2, 0, 0);
  init_mat_mp(randMat, 0, 0);

  // determine if we have a multiplicity matrix of order == 1, == 2 or > 2
  if (order == 1)
  { // we compute the size of the null space normally
    nullSize = corank_rrv_mp(CN, smallest_nonzero, largest_zero, MM1, MM2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }
  else if (order == 2)
  { // we randomize the first block of columns to the correct rank and then do a QR decomposition on that block
    change_size_mat_mp(randMat, colSize[0], ranks[0]);
    make_matrix_random_mp(randMat, colSize[0], ranks[0], mpf_get_default_prec());

    // setup tempR1 to be the randomized columns of the first block of MM1
    r = rowSize[0];
    c = ranks[0];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = colSize[0];
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_mp(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_mp(&tempR1->entry[i][j], &MM1->entry[i][k], &randMat->entry[k][j]);
        }
      } 

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs = ranks[0];
    *houseHolderVecs = (vec_mp *)bmalloc(*numHHVecs * sizeof(vec_mp));
    QR_houseHolder_mp(*houseHolderVecs, tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[0];
    change_size_mat_mp(tempR1, r, c);
    change_size_mat_mp(tempR2, r, c);
    tempR1->rows = tempR2->rows = r;
    tempR1->cols = tempR2->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[0] + j;
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1->entry[i][l]);
        set_mp(&tempR2->entry[i][j], &MM2->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_mp(tempR1, tempR2, *numHHVecs, *houseHolderVecs); 

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_mp(&tempR2->entry[i][j], &tempR2->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2->rows = r;

    // find the size of the new null space 
    nullSize = corank_rrv_mp(CN, smallest_nonzero, largest_zero, tempR1, tempR2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }
  else // order >= 3
  { // setup tempR1 to be last known block of MM1

    r = rowSize[order - 2];
    c = colSize[order - 2] - colSize[order - 3];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;
    for (j = 0; j < c; j++)
    {
      l = j + colSize[order - 3];
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempR1
    apply_HH_mp(tempR1, *numHHVecs, *houseHolderVecs);

    // remove the top part of tempR1
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    change_size_mat_mp(tempR2, r, c);
    tempR2->rows = r;
    tempR2->cols = c;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR2->entry[i][j], &tempR1->entry[l][j]);
      }
    }

    // setup a random matrix to randomize to the correct rank
    change_size_mat_mp(randMat, c, ranks[order - 2]);
    make_matrix_random_mp(randMat, c, ranks[order - 2], mpf_get_default_prec());

    // setup tempR1 to be the randomized columns of the last known block of MM1
    r = tempR2->rows;
    c = ranks[order - 2];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = randMat->rows;
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_mp(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_mp(&tempR1->entry[i][j], &tempR2->entry[i][k], &randMat->entry[k][j]);
        }
      }

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs += c;
    *houseHolderVecs = (vec_mp *)brealloc(*houseHolderVecs, *numHHVecs * sizeof(vec_mp));
    QR_houseHolder_mp(&(*houseHolderVecs)[*numHHVecs - c], tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[order - 2];
    change_size_mat_mp(tempR1, r, c);
    change_size_mat_mp(tempR2, r, c);
    tempR1->rows = tempR2->rows = r;
    tempR1->cols = tempR2->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[order - 2] + j;
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1->entry[i][l]);
        set_mp(&tempR2->entry[i][j], &MM2->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_mp(tempR1, tempR2, *numHHVecs, *houseHolderVecs);

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_mp(&tempR2->entry[i][j], &tempR2->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2->rows = r;

    // find the size of the new null space
    nullSize = corank_rrv_mp(CN, smallest_nonzero, largest_zero, tempR1, tempR2, 0, 0, max_CN, max_SV_ratio, SV_tol);
  }

  // clear memory
  free(str);
  mpf_clear(tol_sign);
  clear_mat_mp(tempR1);
  clear_mat_mp(tempR2);
  clear_mat_mp(randMat);

  return nullSize;
}

// compute a QR decomposition on A without computing Q and without pivoting, but does store the Householder vectors
void QR_houseHolder_mp(vec_mp *houseHolderVecs, mat_mp A, mpf_t tol_sign)
/***************************************************************\
* USAGE: finds A = Q * A - without finding Q explicitly         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume A has full column rank                          *
\***************************************************************/
{
  int i, k, m = A->rows, n = A->cols, *perm = (int *)bmalloc(n * sizeof(int));
  mpf_t tempMPF, norm;
  vec_mp z;

  // initialize
  mpf_init(tempMPF); mpf_init(norm);
  init_vec_mp(z, m);

  // setup perm to be the identity
  for (i = 0; i < n; i++)
    perm[i] = i;

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  { // compute the norm of the kth column
    mpf_set_ui(norm, 0);
    for (i = k; i < m; i++)
    {
      norm_sqr_mp(tempMPF, &A->entry[i][k]);
      mpf_add(norm, norm, tempMPF);
    }
    mpf_sqrt(norm, norm);

    // generate the Householder vector
    init_vec_mp(houseHolderVecs[k], 0);
    genhh_mp(houseHolderVecs[k], A, k, k, norm, perm, tol_sign);
    // find the z associated to the Householder vector & A
    findZ_mat_mp(z, A, houseHolderVecs[k], k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in the Householder vector & every entry in z
    apphh_Z_mat_mp(A, houseHolderVecs[k], z, k, m, k+1, n, 0, houseHolderVecs[k]->size, 0, z->size, perm);
  }

  // clear memory
  free(perm);
  mpf_clear(tempMPF); mpf_clear(norm);
  clear_vec_mp(z);

  return;
}

// apply Householder vectors to a matrices
void apply_HH_mp(mat_mp A, int num, vec_mp *HHVecs)
/***************************************************************\
* USAGE: apply the Householder vectors to A                     *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int k, m = A->rows, n = A->cols, *perm = (int *)bmalloc(n * sizeof(int));
  vec_mp z;

  // initialize z
  init_vec_mp(z, m);

  // setup perm to be the identity
  for (k = 0; k < n; k++)
    perm[k] = k;

  // perform the Householder updates
  for (k = 0; k < num; k++)
  { // set the sizes that are needed for this Householder vector
    A->rows = HHVecs[k]->size + k;

    // find the z associated to HHVecs[k] & A
    findZ_mat_mp(z, A, HHVecs[k], k, 0, perm);
    // apply z to rows k to rows & cols 0 to cols of A using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_mp(A, HHVecs[k], z, k, A->rows, 0, A->cols, 0, HHVecs[k]->size, 0, z->size, perm);
  }

  // set the sizes back
  A->rows = m;
  A->cols = n;

  // clear memory
  free(perm);
  clear_vec_mp(z);

  return;
}

// apply Householder vectors to a pair of matrices
void apply_HH_pair_mp(mat_mp A1, mat_mp A2, int num, vec_mp *HHVecs)
/***************************************************************\
* USAGE: apply the Householder vectors to A1 & A2               *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int k, m = A1->rows, n = A1->cols, *perm = (int *)bmalloc(n * sizeof(int));
  vec_mp z;

  // initialize z
  init_vec_mp(z, m);

  // setup perm to be the identity
  for (k = 0; k < n; k++)
    perm[k] = k;

  // perform the Householder updates
  for (k = 0; k < num; k++)
  { // set the sizes that are needed for this Householder vector
    A1->rows = A2->rows = HHVecs[k]->size + k;

    // find the z associated to HHVecs[k] & A1
    findZ_mat_mp(z, A1, HHVecs[k], k, 0, perm);
    // apply z to rows k to rows & cols 0 to cols of A1 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_mp(A1, HHVecs[k], z, k, A1->rows, 0, A1->cols, 0, HHVecs[k]->size, 0, z->size, perm);

    // find the z associated to HHVecs[k] & A2
    findZ_mat_mp(z, A2, HHVecs[k], k, 0, perm);
    // apply z to rows k to m & cols 0 to n of A2 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_mp(A2, HHVecs[k], z, k, A2->rows, 0, A2->cols, 0, HHVecs[k]->size, 0, z->size, perm);
  }

  // set the sizes back
  A1->rows = A2->rows = m;
  A1->cols = A2->cols = n;

  // clear memory
  free(perm);
  clear_vec_mp(z);

  return;
}

/////// ADAPTIVE MULTI PRECISION ///////////

// take in 2 similar matrices and determine the rank (number of linealy independent columns) efficiently
int rank_MM_amp_mp_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp MM1_mp, int prec1, mat_d MM2_d, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs_d, vec_mp **houseHolderVecs_mp, double max_CN, double max_SV_ratio);
int corank_rrv_amp_mp_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp mat0, int prec0, mat_d mat1, int rowTop, int colTop, double max_CN, double max_FV_ratio);
void QR_houseHolder_amp(vec_d *houseHolderVecs_d, vec_mp *houseHolderVecs_mp, mat_mp A, mpf_t tol_sign);
void apply_HH_pair_amp(mat_mp A1, mat_d A2, int num, vec_d *HHVecs_d, vec_mp *HHVecs_mp);

int rank_MM_amp(double *CN, double *smallest_nonzero, double *largest_zero, mat_d MM1_d, mat_mp MM1_mp, int prec1, mat_d MM2_d, mat_mp MM2_mp, int prec2, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs_d, vec_mp **houseHolderVecs_mp, double max_CN, double max_SV_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of the new null space                     *
* NOTES: Assume MM1 is more accurate                            *
\***************************************************************/
{
  int tempInt, max_prec, nullSize = 0;
  double SV_tol; 

  // max sure that prec1 is greater than prec2 (more accurate) - otherwise, swap them
  if (prec1 < prec2)
  { // swap MM1 & MM2 so that prec1 >= prec2
    if (prec1 < 64)
    { // copy MM1_d to MM2_d and MM2_mp to MM1_mp
      mat_cp_d(MM2_d, MM1_d);
      setprec_mat_mp(MM1_mp, prec2);
      mat_cp_mp(MM1_mp, MM2_mp);
    }
    else
    { // both are using _mp - swap MM1_mp & MM2_mp
      mat_mp randMat;
      init_mat_mp2(randMat, 0, 0, prec1);

      mat_cp_mp(randMat, MM1_mp);
      setprec_mat_mp(MM1_mp, prec2);
      mat_cp_mp(MM1_mp, MM2_mp);
      setprec_mat_mp(MM2_mp, prec1);
      mat_cp_mp(MM2_mp, randMat);

      clear_mat_mp(randMat);
    }

    // swap prec1 & prec2
    tempInt = prec1;
    prec1 = prec2;
    prec2 = tempInt;
  }

  // so we have prec1 >= prec2
  max_prec = prec1;

  // see if both are in double precision
  if (max_prec < 64)
  { // use double precision decomposition
    SV_tol = 1e-14;
    nullSize = rank_MM_d(CN, smallest_nonzero, largest_zero, MM1_d, MM2_d, order, rowSize, colSize, ranks, numHHVecs, houseHolderVecs_d, max_CN, max_SV_ratio, SV_tol);
  }
  else 
  { // we have atleast MM1 using _mp - set the default precision
    initMP(max_prec);

    if (prec2 < 64)
    { // MM1 using _mp & MM2 using _d
      nullSize = rank_MM_amp_mp_d(CN, smallest_nonzero, largest_zero, MM1_mp, prec1, MM2_d, order, rowSize, colSize, ranks, numHHVecs, houseHolderVecs_d, houseHolderVecs_mp, max_CN, max_SV_ratio);
    }
    else
    { // both in MP
      tempInt = prec_to_digits(max_prec) - 4;
      SV_tol = pow(10, -tempInt - 2);
      nullSize = rank_MM_mp(CN, smallest_nonzero, largest_zero, MM1_mp, MM2_mp, order, rowSize, colSize, ranks, numHHVecs, houseHolderVecs_mp, max_CN, max_SV_ratio, SV_tol);
    }
  }

  return nullSize;
}

int rank_MM_amp_mp_d(double *CN, double *smallest_nonzero, double *largest_zero, mat_mp MM1_mp, int prec1, mat_d MM2_d, int order, int *rowSize, int *colSize, int *ranks, int *numHHVecs, vec_d **houseHolderVecs_d, vec_mp **houseHolderVecs_mp, double max_CN, double max_SV_ratio)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: size of the new null space                     *
* NOTES: Assume MM1 is more accurate                            *
\***************************************************************/
{
  int i, j, k, l, r, c, prec_digits, nullSize = 0, m = MM1_mp->rows, n = MM1_mp->cols, curr_prec = prec1;
  size_t size;
  char *str = NULL;
  mpf_t tol_sign;
  mat_d tempR2_d;
  mat_mp tempR1, tempR2, randMat;

  // setup prec_digits
  prec_digits = prec_to_digits(curr_prec) - 3;

  // setup tol_sign
  mpf_init(tol_sign);
  size = 1 + snprintf(NULL, 0, "1e-%d", 2 * prec_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", 2 * prec_digits);
  mpf_set_str(tol_sign, str, 10);  
  
  // error checking
  if (MM1_mp->rows != MM2_d->rows || MM1_mp->cols != MM2_d->cols)
  {
    printf("ERROR: The matrices need to be the same size in rank_MM_amp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (m == 0 || n == 0)
  {
    printf("ERROR: The matrices need to have both columns and rows in rank_MM_amp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (order < 1)
  {
    printf("ERROR: The order of the multiplicity matrix must be atleast 1 in rank_MM_amp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize memory
  init_mat_d(tempR2_d, 0, 0);
  init_mat_mp2(tempR1, 0, 0, curr_prec);
  init_mat_mp2(tempR2, 0, 0, curr_prec);
  init_mat_mp2(randMat, 0, 0, curr_prec);

  // determine if we have a multiplicity matrix of order == 1, == 2 or > 2
  if (order == 1)
  { // we compute the size of the null space normally
    nullSize = corank_rrv_amp_mp_d(CN, smallest_nonzero, largest_zero, MM1_mp, prec1, MM2_d, 0, 0, max_CN, max_SV_ratio);
  }
  else if (order == 2)
  { // we randomize the first block of columns to the correct rank and then do a QR decomposition on that block
    change_size_mat_mp(randMat, colSize[0], ranks[0]);
    make_matrix_random_mp(randMat, colSize[0], ranks[0], curr_prec);

    // setup tempR1 to be the randomized columns of the first block of MM1
    r = rowSize[0];
    c = ranks[0];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = colSize[0];
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_mp(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_mp(&tempR1->entry[i][j], &MM1_mp->entry[i][k], &randMat->entry[k][j]);
        }
      } 

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs = ranks[0];
    *houseHolderVecs_d = (vec_d *)bmalloc(*numHHVecs * sizeof(vec_d));
    *houseHolderVecs_mp = (vec_mp *)bmalloc(*numHHVecs * sizeof(vec_mp));
    QR_houseHolder_amp(*houseHolderVecs_d, *houseHolderVecs_mp, tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[0];
    change_size_mat_mp(tempR1, r, c);
    change_size_mat_d(tempR2_d, r, c);
    tempR1->rows = tempR2_d->rows = r;
    tempR1->cols = tempR2_d->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[0] + j;
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1_mp->entry[i][l]);
        set_d(&tempR2_d->entry[i][j], &MM2_d->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_amp(tempR1, tempR2_d, *numHHVecs, *houseHolderVecs_d, *houseHolderVecs_mp); 

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_d(&tempR2_d->entry[i][j], &tempR2_d->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2_d->rows = r;

    // find the size of the new null space 
    nullSize = corank_rrv_amp_mp_d(CN, smallest_nonzero, largest_zero, tempR1, prec1, tempR2_d, 0, 0, max_CN, max_SV_ratio);
  }
  else // order >= 3
  { // setup tempR1 to be last known block of MM1

    r = rowSize[order - 2];
    c = colSize[order - 2] - colSize[order - 3];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;
    for (j = 0; j < c; j++)
    {
      l = j + colSize[order - 3];
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1_mp->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempR1
    apply_HH_mp(tempR1, *numHHVecs, *houseHolderVecs_mp);

    // remove the top part of tempR1
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    change_size_mat_mp(tempR2, r, c);
    tempR2->rows = r;
    tempR2->cols = c;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR2->entry[i][j], &tempR1->entry[l][j]);
      }
    }

    // setup a random matrix to randomize to the correct rank
    change_size_mat_mp(randMat, c, ranks[order - 2]);
    make_matrix_random_mp(randMat, c, ranks[order - 2], curr_prec);

    // setup tempR1 to be the randomized columns of the last known block of MM1
    r = tempR2->rows;
    c = ranks[order - 2];
    change_size_mat_mp(tempR1, r, c);
    tempR1->rows = r;
    tempR1->cols = c;

    l = randMat->rows;
    for (i = 0; i < r; i++)
      for (j = 0; j < c; j++)
      { // compute tempR1[i][j]
        set_zero_mp(&tempR1->entry[i][j]);
        for (k = 0; k < l; k++)
        { // [i][j] += [i][k] * [k][j]
          sum_mul_mp(&tempR1->entry[i][j], &tempR2->entry[i][k], &randMat->entry[k][j]);
        }
      }

    // we need to do a QR decomposition for tempR1 to find the associated Householder vectors
    *numHHVecs += c;
    *houseHolderVecs_d = (vec_d *)brealloc(*houseHolderVecs_d, *numHHVecs * sizeof(vec_d));
    *houseHolderVecs_mp = (vec_mp *)brealloc(*houseHolderVecs_mp, *numHHVecs * sizeof(vec_mp));
    QR_houseHolder_amp(&(*houseHolderVecs_d)[*numHHVecs - c], &(*houseHolderVecs_mp)[*numHHVecs - c], tempR1, tol_sign);

    // setup tempRi to be the last part of MMi
    r = m;
    c = n - colSize[order - 2];
    change_size_mat_mp(tempR1, r, c);
    change_size_mat_d(tempR2_d, r, c);
    tempR1->rows = tempR2_d->rows = r;
    tempR1->cols = tempR2_d->cols = c;

    for (j = 0; j < c; j++)
    { // copy over the columns
      l = colSize[order - 2] + j;
      for (i = 0; i < r; i++)
      {
        set_mp(&tempR1->entry[i][j], &MM1_mp->entry[i][l]);
        set_d(&tempR2_d->entry[i][j], &MM2_d->entry[i][l]);
      }
    }

    // apply the Householder vectors to tempRi
    apply_HH_pair_amp(tempR1, tempR2_d, *numHHVecs, *houseHolderVecs_d, *houseHolderVecs_mp);

    // remove the top part of tempRi
    r = tempR1->rows - *numHHVecs;
    c = tempR1->cols;
    for (i = 0; i < r; i++)
    { // move up row i + numHHVecs to i
      l = i + *numHHVecs;
      for (j = 0; j < c; j++)
      {
        set_mp(&tempR1->entry[i][j], &tempR1->entry[l][j]);
        set_d(&tempR2_d->entry[i][j], &tempR2_d->entry[l][j]);
      }
    }
    // reset the size
    tempR1->rows = tempR2_d->rows = r;

    // find the size of the new null space
    nullSize = corank_rrv_amp_mp_d(CN, smallest_nonzero, largest_zero, tempR1, prec1, tempR2_d, 0, 0, max_CN, max_SV_ratio);
  }

  // clear memory
  free(str);
  mpf_clear(tol_sign);
  clear_mat_d(tempR2_d);
  clear_mat_mp(tempR1);
  clear_mat_mp(tempR2);
  clear_mat_mp(randMat);

  return nullSize;
}

// compute a QR decomposition on A without computing Q and without pivoting, but does store the Householder vectors
void QR_houseHolder_amp(vec_d *houseHolderVecs_d, vec_mp *houseHolderVecs, mat_mp A, mpf_t tol_sign)
/***************************************************************\
* USAGE: finds A = Q * A - without finding Q explicitly         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume A has full column rank                          *
* setup HH_mp and then convert to HH_d                          *
\***************************************************************/
{
  int i, k, m = A->rows, n = A->cols, *perm = (int *)bmalloc(n * sizeof(int));
  mpf_t tempMPF, norm;
  vec_mp z;

  // initialize
  mpf_init(tempMPF); mpf_init(norm);
  init_vec_mp(z, m);

  // setup perm to be the identity
  for (i = 0; i < n; i++)
    perm[i] = i;

  // main algorithm - do pivoting for column k
  for (k = 0; k < n; k++)
  { // compute the norm of the kth column
    mpf_set_ui(norm, 0);
    for (i = k; i < m; i++)
    {
      norm_sqr_mp(tempMPF, &A->entry[i][k]);
      mpf_add(norm, norm, tempMPF);
    }
    mpf_sqrt(norm, norm);

    // generate the Householder vector
    init_vec_mp(houseHolderVecs[k], 0);
    genhh_mp(houseHolderVecs[k], A, k, k, norm, perm, tol_sign);
    // find the z associated to the Householder vector & A
    findZ_mat_mp(z, A, houseHolderVecs[k], k, k+1, perm);
    // apply z to rows k to m & cols k+1 to n of R using every entry in the Householder vector & every entry in z
    apphh_Z_mat_mp(A, houseHolderVecs[k], z, k, m, k+1, n, 0, houseHolderVecs[k]->size, 0, z->size, perm);

    // convert to _d
    init_vec_d(houseHolderVecs_d[k], houseHolderVecs[k]->size);
    vec_mp_to_d(houseHolderVecs_d[k], houseHolderVecs[k]);
  }

  // clear memory
  free(perm);
  mpf_clear(tempMPF); mpf_clear(norm);
  clear_vec_mp(z);

  return;
}

// apply Householder vectors to a pair of matrices
void apply_HH_pair_amp(mat_mp A1, mat_d A2, int num, vec_d *HHVecs_d, vec_mp *HHVecs_mp)
/***************************************************************\
* USAGE: apply the Householder vectors to A1 & A2               *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int k, m = A1->rows, n = A1->cols, *perm = (int *)bmalloc(n * sizeof(int));
  vec_d z_d;
  vec_mp z_mp;

  // initialize z
  init_vec_d(z_d, m);
  init_vec_mp(z_mp, m);

  // setup perm to be the identity
  for (k = 0; k < n; k++)
    perm[k] = k;

  // perform the Householder updates
  for (k = 0; k < num; k++)
  { // set the sizes that are needed for this Householder vector
    A1->rows = A2->rows = HHVecs_mp[k]->size + k;

    // find the z associated to HHVecs_mp[k] & A1
    findZ_mat_mp(z_mp, A1, HHVecs_mp[k], k, 0, perm);
    // apply z to rows k to rows & cols 0 to cols of A1 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_mp(A1, HHVecs_mp[k], z_mp, k, A1->rows, 0, A1->cols, 0, HHVecs_mp[k]->size, 0, z_mp->size, perm);

    // find the z associated to HHVecs_d[k] & A2
    findZ_mat_d(z_d, A2, HHVecs_d[k], k, 0, perm);
    // apply z to rows k to m & cols 0 to n of A2 using every entry in HHVecs[k] & every entry in z
    apphh_Z_mat_d(A2, HHVecs_d[k], z_d, k, A2->rows, 0, A2->cols, 0, HHVecs_d[k]->size, 0, z_d->size, perm);
  }

  // set the sizes back
  A1->rows = A2->rows = m;
  A1->cols = A2->cols = n;

  // clear memory
  free(perm);
  clear_vec_d(z_d);
  clear_vec_mp(z_mp);

  return;
}





