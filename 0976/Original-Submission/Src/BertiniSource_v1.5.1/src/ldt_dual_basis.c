// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "localdim.h"

void generate_dual_basis_d(vec_d *dual_d, int multiplicity, mat_d fullRankMat, vec_d *randVec);
void generate_dual_basis_mp(vec_mp *dual_mp, int multiplicity, mat_mp fullRankMat, vec_mp *randVec, int curr_prec);

void create_dual_basis(vec_d **dual_d, vec_mp **dual_mp, int *dual_prec, int multiplicity, mat_d MM_d, mat_mp MM_mp, int MM_prec, int MPType, mat_d topMat_d, mat_mp topMat_mp, vec_d **randVec_d, vec_mp **randVec_mp, int setupTopMat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the dual basis since the multiplicity is known *
\***************************************************************/
{
  int i, j, rows, cols;

  // verify that multiplicity >= 1
  if (multiplicity < 1)
  {
    printf("ERROR: The multiplicity must be positive when creating a dual basis!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (MPType == 0 || (MPType == 2 && MM_prec < 64))
  { // compute the dual basis in double precision
    mat_d fullRankMat;
    init_mat_d(fullRankMat, 0, 0);

    // allocate dual_d & setup dual_prec
    *dual_prec = 52;
    *dual_d = (vec_d *)bmalloc(multiplicity * sizeof(vec_d));
    for (i = 0; i < multiplicity; i++)
      init_vec_d((*dual_d)[i], 0);

    // setup rows & cols
    rows = MM_d->rows;
    cols = MM_d->cols;

    if (setupTopMat)
    { // create a random matrix to kill off the rank deficiency
      make_matrix_random_d(topMat_d, multiplicity, cols);

      // generate random vectors
      *randVec_d = (vec_d *)bmalloc(multiplicity * sizeof(vec_d));
      for (i = 0; i < multiplicity; i++)
      {
        init_vec_d((*randVec_d)[i], multiplicity);
        make_vec_random_d((*randVec_d)[i], multiplicity);
      }
    }

    // setup fullRankMat == [[topMat];[MM]];
    change_size_mat_d(fullRankMat, rows + multiplicity, cols);
    fullRankMat->rows = rows + multiplicity;
    fullRankMat->cols = cols;
    for (i = rows + multiplicity - 1; i >= 0; i--)
      if (i < multiplicity) 
      { // copy topMat to top of fullRankMat
        for (j = 0; j < cols; j++)
          set_d(&fullRankMat->entry[i][j], &topMat_d->entry[i][j]);
      }
      else
      { // copy MM_d to bottom of fullRankMat
        for (j = 0; j < cols; j++)
          set_d(&fullRankMat->entry[i][j], &MM_d->entry[i - multiplicity][j]);
      }

    // generate the dual basis vectors
    generate_dual_basis_d(*dual_d, multiplicity, fullRankMat, *randVec_d);

    // clear memory
    clear_mat_d(fullRankMat);
  }
  else
  { // compute the dual basis using multi precision
    
    // set the overall precision
    initMP(MM_prec);

    mat_mp fullRankMat;
    init_mat_mp(fullRankMat, 0, 0);

    // allocate dual_mp & setup dual_prec
    *dual_prec = MM_prec;
    *dual_mp = (vec_mp *)bmalloc(multiplicity * sizeof(vec_mp));
    for (i = 0; i < multiplicity; i++)
      init_vec_mp((*dual_mp)[i], 0);

    // setup rows & cols
    rows = MM_mp->rows;
    cols = MM_mp->cols;

    if (setupTopMat)
    { // create a random matrix to kill off the rank deficiency
      make_matrix_random_mp(topMat_mp, multiplicity, cols, MM_prec);

      // generate random vectors
      *randVec_mp = (vec_mp *)bmalloc(multiplicity * sizeof(vec_mp));
      for (i = 0; i < multiplicity; i++)
      {
        init_vec_mp((*randVec_mp)[i], multiplicity);
        make_vec_random_mp((*randVec_mp)[i], multiplicity);
      }
    }

    // setup fullRankMat == [[topMat];[MM]];
    change_size_mat_mp(fullRankMat, rows + multiplicity, cols);
    fullRankMat->rows = rows + multiplicity;
    fullRankMat->cols = cols;
    for (i = rows + multiplicity - 1; i >= 0; i--)
      if (i < multiplicity)
      { // copy topMat to top of fullRankMat
        for (j = 0; j < cols; j++)
          set_mp(&fullRankMat->entry[i][j], &topMat_mp->entry[i][j]);
      }
      else
      { // copy MM_d to bottom of fullRankMat
        for (j = 0; j < cols; j++)
          set_mp(&fullRankMat->entry[i][j], &MM_mp->entry[i - multiplicity][j]);
      }

    // generate the dual basis vectors
    generate_dual_basis_mp(*dual_mp, multiplicity, fullRankMat, *randVec_mp, MM_prec);

    // clear memory
    clear_mat_mp(fullRankMat);
  }

  return;
}

void generate_dual_basis_d(vec_d *dual_d, int multiplicity, mat_d fullRankMat, vec_d *randVec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the dual basis since the multiplicity is known *
\***************************************************************/
{
  int i, j, rows, cols, *perm = NULL;
  vec_d b;
  mat_d Q, Q_trans, R, P;

  init_vec_d(b, 0);
  init_mat_d(Q, 0, 0); init_mat_d(Q_trans, 0, 0);
  init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

  // find a QR decomposition of fullRankMat
  QR_d_prec(Q, R, P, fullRankMat, 0);

  // make R a full rank square matrix by removing extra zeros from the bottom
  R->rows = R->cols = fullRankMat->cols;
 
  // setup Q to the correct size
  Q->rows = fullRankMat->rows;
  Q->cols = fullRankMat->cols;

  // transpose Q
  transpose_d(Q_trans, Q);

  // convert P to perm
  perm = (int *)bmalloc(P->rows * sizeof(int));
  convertToPerm_d(perm, P, P->rows);

  // setup rows & cols
  rows = fullRankMat->rows;
  cols = fullRankMat->cols;

  // solve for the dual basis
  for (i = 0; i < multiplicity; i++)
  { // setup b = [[rand_vec];[0]];
    change_size_vec_d(b, rows);
    b->size = rows;
    for (j = 0; j < rows; j++) 
      if (j < multiplicity)
      { // setup top to be randVec
        set_d(&b->coord[j], &randVec[i]->coord[j]);
      }
      else
      { // setup bottom to be 0
        set_zero_d(&b->coord[j]);
      }

    // setup the ith dual basis vector
    matrixSolve_from_QR_d(dual_d[i], Q_trans, R, perm, b);

    // normalize the ith dual basis vector
    normalize_vec_d(dual_d[i], dual_d[i]);
  }

  // clear memory
  free(perm);
  clear_vec_d(b);
  clear_mat_d(Q); clear_mat_d(Q_trans);
  clear_mat_d(R); clear_mat_d(P);

  return;
}

void generate_dual_basis_mp(vec_mp *dual_mp, int multiplicity, mat_mp fullRankMat, vec_mp *randVec, int curr_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the dual basis since the multiplicity is known *
\***************************************************************/
{
  int i, j, rows, cols, *perm = NULL;
  vec_mp b;
  mat_mp Q, Q_trans, R, P;

  init_vec_mp(b, 0);
  init_mat_mp(Q, 0, 0); init_mat_mp(Q_trans, 0, 0);
  init_mat_mp(R, 0, 0); init_mat_mp(P, 0, 0);

  // find a QR decomposition of fullRankMat
  QR_mp_prec(Q, R, P, fullRankMat, 0, curr_prec);

  // make R a full rank square matrix by removing extra zeros from the bottom
  R->rows = R->cols = fullRankMat->cols;

  // setup Q to the correct size
  Q->rows = fullRankMat->rows;
  Q->cols = fullRankMat->cols;

  // transpose Q
  transpose_mp(Q_trans, Q);

  // convert P to perm
  perm = (int *)bmalloc(P->rows * sizeof(int));
  convertToPerm_mp(perm, P, P->rows);

  // setup rows & cols
  rows = fullRankMat->rows;
  cols = fullRankMat->cols;

  // solve for the dual basis
  for (i = 0; i < multiplicity; i++)
  { // setup b = [[rand_vec];[0]];
    change_size_vec_mp(b, rows);
    b->size = rows;
    for (j = 0; j < rows; j++)
      if (j < multiplicity)
      { // setup top to be randVec
        set_mp(&b->coord[j], &randVec[i]->coord[j]);
      }
      else
      { // setup bottom to be 0
        set_zero_mp(&b->coord[j]);
      }

    // setup the ith dual basis vector
    matrixSolve_from_QR_mp(dual_mp[i], Q_trans, R, perm, b);

    // normalize the ith dual basis vector
    normalize_vec_mp(dual_mp[i], dual_mp[i]);
  }

  // clear memory
  free(perm);
  clear_vec_mp(b);
  clear_mat_mp(Q); clear_mat_mp(Q_trans);
  clear_mat_mp(R); clear_mat_mp(P);

  return;
}

// compute the hilbert function from a dual basis

void find_hilbert_func_d(int **hilFn, int *hilOrder, int *reg, int mult, vec_d *dual1_d, vec_d *dual2_d, int *colSize, int colOrder, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the hilbert function of the ideal from dual basis *
\***************************************************************/
{
  int i, j, rows, cols, corank, curr_order, numHH = 0;
  int *ranks = NULL, *rowSize = NULL;
  double max_CN, max_SV_ratio, SV_tol, CN, smallest_nonzero, largest_zero;
  vec_d *HH = NULL;
  mat_d A1, A2, MM1, MM2;

  // verify that multiplicity >= 1
  if (mult < 1)
  {
    printf("ERROR: The multiplicity must be positive when creating a hilbert function!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize A1, A2, MM1, MM2
  init_mat_d(A1, 0, 0);
  init_mat_d(A2, 0, 0);
  init_mat_d(MM1, 0, 0);
  init_mat_d(MM2, 0, 0);

  // setup the tolerances
  max_CN = 1e13;
  max_SV_ratio = T->ratioTol; 
  SV_tol = MAX(T->sing_val_zero_tol, 1e-15);

  // setup hilbertFn
  *hilOrder = 1;
  *hilFn = (int *)bmalloc(*hilOrder * sizeof(int));
  (*hilFn)[*hilOrder - 1] = 1;

  // initialize reg
  *reg = -1;

  // setup rows & cols
  rows = mult;
  cols = dual1_d[0]->size;

  // copy the dual basis vectors to A1 & A2 - a perturbation of A1
  change_size_mat_d(A1, rows, cols);
  change_size_mat_d(A2, rows, cols);
  A1->rows = A2->rows = rows;
  A1->cols = A2->cols = cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      set_d(&A1->entry[i][j], &dual1_d[i]->coord[j]);
      set_d(&A2->entry[i][j], &dual2_d[i]->coord[j]);
    }

  // find the rank of A[i]
  corank = curr_order = 0;
  do
  { // increase the order
    curr_order++;

    // increase the size of the hilbert function
    *hilOrder += 1;
    *hilFn = (int *)brealloc(*hilFn, *hilOrder * sizeof(int));

    // setup MM1 & MM2
    A1->cols = A2->cols = colSize[curr_order - 1];
    mat_cp_d(MM1, A1);
    mat_cp_d(MM2, A2);

    // setup rowSize
    rowSize = (int *)brealloc(rowSize, curr_order * sizeof(int));
    rowSize[curr_order - 1] = rows;

    // find the corank
    corank = rank_MM_d(&CN, &smallest_nonzero, &largest_zero, MM1, MM2, curr_order, rowSize, colSize, ranks, &numHH, &HH, max_CN, max_SV_ratio, SV_tol);

    // store the ranks
    ranks = (int *)brealloc(ranks, curr_order * sizeof(int));
    ranks[curr_order - 1] = curr_order == 1 ? colSize[curr_order - 1] - corank : colSize[curr_order - 1] - colSize[curr_order - 2] - corank;

    // find hilFn - new rank
    (*hilFn)[*hilOrder - 1] = *hilOrder == 2 ? ranks[curr_order - 1] - (*hilFn)[*hilOrder - 2] : ranks[curr_order - 1];

    // see if we have reached the regularity
    if ((*hilFn)[*hilOrder - 1] == 0 && *reg < 0)
      *reg = *hilOrder - 1;

  } while ((*hilFn)[*hilOrder - 1] > 0 && curr_order < colOrder);
 
  // clear memory
  clear_mat_d(A1);
  clear_mat_d(A2);
  clear_mat_d(MM1);
  clear_mat_d(MM2);

  for (curr_order = numHH - 1; curr_order >= 0; curr_order--)
    clear_vec_d(HH[curr_order]);
  free(HH);
  free(rowSize);  free(ranks);

  return;
}

void find_hilbert_func_mp(int **hilFn, int *hilOrder, int *reg, int mult, vec_mp *dual1_mp, vec_mp *dual2_mp, int *colSize, int colOrder, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the hilbert function of the ideal from dual basis *
\***************************************************************/
{
  int i, j, rows, cols, corank, curr_order, numHH = 0;
  int *ranks = NULL, *rowSize = NULL;
  double max_CN, max_SV_ratio, SV_tol, CN, smallest_nonzero, largest_zero;
  vec_mp *HH = NULL;
  mat_mp A1, A2, MM1, MM2;

  // verify that multiplicity >= 1
  if (mult < 1)
  {
    printf("ERROR: The multiplicity must be positive when creating a hilbert function!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize A1, A2, MM1, MM2
  init_mat_mp(A1, 0, 0);
  init_mat_mp(A2, 0, 0);
  init_mat_mp(MM1, 0, 0);
  init_mat_mp(MM2, 0, 0);

  // setup the tolerances
  corank = prec_to_digits(T->Precision) - 4;
  max_CN = MIN(1e150, pow(10, corank));
  max_SV_ratio = T->ratioTol; 
  SV_tol = MAX(T->sing_val_zero_tol, pow(10, -corank - 2));

  // setup hilbertFn
  *hilOrder = 1;
  *hilFn = (int *)bmalloc(*hilOrder * sizeof(int));
  (*hilFn)[*hilOrder - 1] = 1;

  // initialize reg
  *reg = -1;

  // setup rows & cols
  rows = mult;
  cols = dual1_mp[0]->size;

  // copy the dual basis vectors to A1 & A2 - a perturbation of A1
  change_size_mat_mp(A1, rows, cols);
  change_size_mat_mp(A2, rows, cols);
  A1->rows = A2->rows = rows;
  A1->cols = A2->cols = cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      set_mp(&A1->entry[i][j], &dual1_mp[i]->coord[j]);
      set_mp(&A2->entry[i][j], &dual2_mp[i]->coord[j]);
    }

  // find the rank of A[i]
  corank = curr_order = 0;
  do
  { // increase the order
    curr_order++;

    // increase the size of the hilbert function
    *hilOrder += 1;
    *hilFn = (int *)brealloc(*hilFn, *hilOrder * sizeof(int));

    // setup MM1 & MM2
    A1->cols = A2->cols = colSize[curr_order - 1];
    mat_cp_mp(MM1, A1);
    mat_cp_mp(MM2, A2);

    // setup rowSize
    rowSize = (int *)brealloc(rowSize, curr_order * sizeof(int));
    rowSize[curr_order - 1] = rows;

    // find the corank
    corank = rank_MM_mp(&CN, &smallest_nonzero, &largest_zero, MM1, MM2, curr_order, rowSize, colSize, ranks, &numHH, &HH, max_CN, max_SV_ratio, SV_tol);

    // store the ranks
    ranks = (int *)brealloc(ranks, curr_order * sizeof(int));
    ranks[curr_order - 1] = curr_order == 1 ? colSize[curr_order - 1] - corank : colSize[curr_order - 1] - colSize[curr_order - 2] - corank;

    // find hilFn - new rank
    (*hilFn)[*hilOrder - 1] = *hilOrder == 2 ? ranks[curr_order - 1] - (*hilFn)[*hilOrder - 2] : ranks[curr_order - 1];

    // see if we have reached the regularity
    if ((*hilFn)[*hilOrder - 1] == 0 && *reg < 0)
      *reg = *hilOrder - 1;

  } while ((*hilFn)[*hilOrder - 1] > 0 && curr_order < colOrder);

  // clear memory
  clear_mat_mp(A1);
  clear_mat_mp(A2);
  clear_mat_mp(MM1);
  clear_mat_mp(MM2);

  for (curr_order = numHH - 1; curr_order >= 0; curr_order--)
    clear_vec_mp(HH[curr_order]);
  free(HH);
  free(rowSize);  free(ranks);

  return;
}

