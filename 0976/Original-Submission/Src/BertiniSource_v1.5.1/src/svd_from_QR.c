// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"

void expand_unitary_d(mat_d U, int startCol, int endCol);
void expand_unitary_mp(mat_mp U, int startCol, int endCol);
 
///////DOUBLE PRECISION ////////

int svd_corank_d(mat_d A, double rank_tol, double largeChange)
/***************************************************************\
* USAGE: finds the corank of A using its SVD                    *
* ARGUMENTS:                                                    *
* RETURN VALUES: the corank of A if successful, else -1 if not  *
* NOTES: 				                        * 
\***************************************************************/
{ 
  int retVal, its = 50, num_digits = -14; // -14 = -(int) floor(52 * log10(2.0) - 1.5), where 52 is the number of bits for double precision
  // QR does not need to have too tight of a tolerance
  // while tol_sign needs to be relatively small for the best results
  double tol_prec = pow(10, num_digits), tol_sign = pow(10, 2 * num_digits);
  vec_d E;
  init_vec_d(E, 0);

  if (rank_tol < tol_prec)
  { // ideally, we would like tol_prec a couple orders of magnitude smaller than rank_tol 
    printf("WARNING: The tolerance to determine the rank is smaller than the precision tolerance for SVD!\nBy default, the rank tolerance will be adjusted.\n");
    rank_tol = tol_prec;
  }

  // try to find the singular values with the tolerances
  retVal = svd_jacobi_E_d(E, A, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);

  // check to see if successful
  if (retVal < 0)
  { // SVD was not successful
    // adjust tolerances to see if we can get convergence - this should very rarely happen!
    its = 100;
    // loosen convergence criterion tol_Jacobi by making larger and loosen tol_sign by making smaller
    tol_prec *= 10;
    tol_sign /= 100;

    retVal = svd_jacobi_E_d(E, A, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);
  }

  clear_vec_d(E);

  return retVal; 
}

int svd_corank_analyze_d(mat_d E, double rank_tol, double largeChange)
/***************************************************************\
* USAGE: analyzes E-diagonal matrix from SVD-to find its corank *
* ARGUMENTS:                                                    *
* RETURN VALUES: the corank of E                                *
* NOTES: Assume diagonal entries of E are decreasing,real & >= 0*
\***************************************************************/
{ 
  int i, rank, m = E->rows, n = E->cols;
  double prevSV, currentSV;

  // do error checking  
  if (rank_tol < 0 || largeChange < 0) 
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_rank_analyze_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // set n = min(m, n) so that E has n diagonal entries to look at
  if (m < n)
    n = m;

  if (n == 0) // no singular values has rank 0
    return 0;
  // so we can assume that n > 0

  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  rank = 0;
  currentSV = E->entry[0][0].r; // initialize to something
  for (i = 0; i < n; i++)
  {
    prevSV = currentSV; // the previous singular value
    currentSV = E->entry[i][i].r; // the current singular value

    if (checkGood_d(currentSV, prevSV, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // this one is good!
      rank++; // this one is good!
    }
    else
    { // this one is no good!
      i = n; // all the rest of them are not good as well
    }
  }

  return (n - rank); 
}

int svd_jacobi_d_prec(mat_d U, mat_d E, mat_d V, mat_d InputMat, double rankTol)
/***************************************************************\
* USAGE: finds the SVD of InputMat using Jacobi iterations      *
* ARGUMENTS: InputMat = U*E*V^H                                 *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  return svd_jacobi_d(U, E, V, InputMat, 100, rankTol, 1e-15, 1e-14, 1e-20, 1e13);
}

int svd_jacobi_d(mat_d U, mat_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: finds the SVD of InputMat using Jacobi iterations      *
* ARGUMENTS: InputMat = U*E*V^H                                 *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, k, retVal, transposed, m = InputMat->rows, n = InputMat->cols;
  int *perm = NULL;
  double norm, norm_inv;
  comp_d tempComp;
  mat_d H, Q, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_mat_d(U, m, n);
    change_size_mat_d(E, n, n);
    change_size_mat_d(V, n, n);
    U->rows = m;
    U->cols = E->rows = E->cols = V->rows = V->cols = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv_Jacobi < 0 || tol_QR < 0 || tol_sign < 0 || largeChange < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // for best results, we should normalize InputMat based on its infinity norm if it is > 1
  norm = infNormMat_d(InputMat);

  // check to make sure that InputMat is 'skinny'
  if (m >= n)
  { // we have more rows than columns so we just proceed
    init_mat_d(H, m, n); init_mat_d(Q, m, m); init_mat_d(R, m, n); init_mat_d(P, n, n);
    H->rows = m;
    H->cols = n;

    if (norm > 1)
    { // normalize InputMat and store as H
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          mul_rdouble_d(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // cp InputMat to H
      mat_cp_d(H, InputMat);
    }

    // did not transpose the matrix
    transposed = 0;
  }
  else
  { // we have more columns than rows so we must work with the transpose
    init_mat_d(H, n, m); init_mat_d(Q, n, n); init_mat_d(R, n, m); init_mat_d(P, m, m);
    H->rows = n;
    H->cols = m;

    if (norm > 1)
    { // normalize InputMat and store transposed as H
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          conjugate_d(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rdouble_d(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // H = InputMat^T
      transpose_d(H, InputMat);
    }

    // transposed the matrix
    transposed = 1;
    m = H->rows;
    n = H->cols;
  }
  // so we have m >= n >= 1

  // do the QR decomposition - 1 at the end so that QR presorts the rows ( H * P = Q * R, where P is n x n, Q is m x m & R is m x n)
  QR_d(Q, R, P, H, tol_QR, tol_sign, largeChange, 1);

  // decompose R using jacobi_d ( R = U * E * V^H, U is m x m, E is m x n and V is n x n)
  // NOTE: jacobi will always return its approximation of the svd even if it does not converge within the desired number of iterations
  retVal = svd_R_jacobi_d(U, E, V, R, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

  if (norm > 1)
  { // undo the normalization by multiplying the entries of E by norm
    for (i = 0; i < n; i++)
    {
      E->entry[i][i].r *= norm;
    }
  }

  // update U based on Q - U = Q * U
  mat_mul_d(U, Q, U);

  // update V based on P -> V = P * V -> do this based no row swaps for efficiency
  perm = (int *)bmalloc(n * sizeof(int));
  convertToPerm_d(perm, P, n);
  for (i = 0; i < n; i++)
    if (perm[i] != i)
    { // swap i & perm[i]
      k = perm[i];
      for (j = 0; j < n; j++)
      { // update V
        set_d(tempComp, &V->entry[i][j]);
        set_d(&V->entry[i][j], &V->entry[k][j]);
        set_d(&V->entry[k][j], tempComp);
      }
      // update perm
      for (j = 0; j < n; j++)
        if (perm[j] == i)
        {
          perm[j] = perm[i];
          perm[i] = i;
          j = n;
        }
    }
  // so we have H = U * E * V^H, thus we need to check if H = InputMat or H = InputMat^H

  if (transposed)
  { // (U * E * V^H) ^ H = (V * E^H * U^H) = V * E * U^H, so we just need to swap V & U
    mat_cp_d(H, V);
    mat_cp_d(V, U);
    mat_cp_d(U, H);
    transpose_d(E, E);
  }

  // clear
  free(perm);
  clear_mat_d(H); clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
 
  return retVal;
}

int svd_R_jacobi_d(mat_d U, mat_d E, mat_d V, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: finds the SVD of R - upper triangular matrix -  using  *
* Jacobi iterations                                             *
* ARGUMENTS: R = U*E*V^H, U - m x m, E - m x n, V - n x n       *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, retVal, m = R->rows, n = R->cols;
  mat_d tempR;

  // do error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more rows than columns in svd_R_jacobi_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  if (n == 0) // if R has no cols - immediately exit
  {
    change_size_mat_d(U, m, m);
    change_size_mat_d(E, m, n);
    change_size_mat_d(V, n, n);
    U->rows = U->cols = E->rows = m;
    E->cols = V->rows = V->cols = n;
    make_matrix_ID_d(U, m, m);
    return 0;
  }

  // so m >= n >= 1
 
  // copy the main n x n part of R to tempR
  init_mat_d(tempR, n, n);
  tempR->rows = tempR->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_d(&tempR->entry[i][j], &R->entry[i][j]);
    }

  // do jacobi iterations on tempR
  retVal = jacobi_d(U, E, V, tempR, its, rank_tol, tol_conv, tol_sign, largeChange);

  // adjust U & E based on m  - U = [[U, 0][0, I]] and E = [[E][0]]
  increase_size_mat_d(U, m, m);
  increase_size_mat_d(E, m, E->cols);
  U->rows = U->cols = E->rows = m;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i < n)
      {
        if (j >= n)
        { // U[i][j] = 0
          set_zero_d(&U->entry[i][j]);
        }
      }
      else // i >= n
      {
        if (j < n)
        { // U[i][j] & E[i][j] = 0
          set_zero_d(&U->entry[i][j]);
          set_zero_d(&E->entry[i][j]);
        }
        else // j >= n
        { // U[i][j] is (i == j)
          set_double_d(&U->entry[i][j], i == j, 0);
        }
      }

  clear_mat_d(tempR);

  return retVal;
}

int jacobi_d(mat_d U, mat_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on InputMat to find its SVD     *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: InputMat has to be square - should be R of a QR decomp *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = InputMat->rows, n = InputMat->cols;
  double da, db, sizeC_sqr, sizeC, stoppingError, prevError1, prevError2, tempD;
  comp_d c, t, tempComp;

  int *perm = NULL;
  double *sv = NULL;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix is not square in jacobi_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv < 0 || tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_mat_d(U, n, n);
    change_size_mat_d(E, n, n);
    change_size_mat_d(V, n, n);
    U->rows = U->cols = E->rows = E->cols = V->rows = V->cols = n;
    return 0;
  }

  // so we can assume that m == n >= 1

  // set U to InputMat, E to 0 and V to Id
  change_size_mat_d(U, n, n);
  change_size_mat_d(E, n, n);
  change_size_mat_d(V, n, n);
  U->rows = U->cols = E->rows = E->cols = V->rows = V->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    { // U
      set_d(&U->entry[i][j], &InputMat->entry[i][j]);
      // E
      set_zero_d(&E->entry[i][j]);
      // V
      if (i < n)
      {
        set_double_d(&V->entry[i][j], i == j, 0);
      }
    }

  // main loop
  count = 0;
  stoppingError = prevError1 = prevError2 = 0; // initialize to something
  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    prevError2 = prevError1;
    prevError1 = stoppingError;
    stoppingError = 0; // reset stoppingError to 0
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for (i = 0; i < n-1; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        da = db = 0;
        set_zero_d(c);
        for (k = 0; k < n; k++)
        {
          da += U->entry[k][i].r * U->entry[k][i].r + U->entry[k][i].i * U->entry[k][i].i;
          db += U->entry[k][j].r * U->entry[k][j].r + U->entry[k][j].i * U->entry[k][j].i;
          conjugate_d(t, &U->entry[k][i]);
          sum_mul_d(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        sizeC_sqr = c->r * c->r + c->i * c->i;
        sizeC = sqrt(sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (da > tol_sign && db > tol_sign && sizeC > tol_sign)
        { // we can safely find the stopping condition estimate
          tempD = sizeC / (sqrt(da) * sqrt(db));
          if (tempD > stoppingError)
            stoppingError = tempD;
        }

        if (sizeC > tol_sign)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          db = db / (2 * sizeC_sqr)  - da / (2 * sizeC_sqr);
          mul_rdouble_d(tempComp, c, db);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_d(t, tempComp, tol_sign);
          db = tempComp->r * tempComp->r + tempComp->i * tempComp->i;
          da = sqrt(db);
          db = 1 / (da + sqrt(1 + db));
          mul_rdouble_d(t, t, db);
        }
        else
        {
          set_zero_d(t);
        }

        // find db = 1 / sqrt(1 + |t|^2)
        db = 1 / sqrt(1 + t->r * t->r + t->i * t->i);
        // find t = db * t
        mul_rdouble_d(t, t, db);

        // tempComp = -conj(t)
        tempComp->r = -t->r;
        tempComp->i = t->i;

        // update columns i & j of U & V
        for (k = 0; k < n; k++)
        { // update U
          set_d(c, &U->entry[k][i]);

          mul_rdouble_d(&U->entry[k][i], c, db);
          sum_mul_d(&U->entry[k][i], tempComp, &U->entry[k][j]);
 
          mul_rdouble_d(&U->entry[k][j], &U->entry[k][j], db);
          sum_mul_d(&U->entry[k][j], c, t);

          // update V
          set_d(c, &V->entry[k][i]);

          mul_rdouble_d(&V->entry[k][i], c, db);
          sum_mul_d(&V->entry[k][i], tempComp, &V->entry[k][j]);

          mul_rdouble_d(&V->entry[k][j], &V->entry[k][j], db);
          sum_mul_d(&V->entry[k][j], c, t);
        }
      }

    // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
    if (stoppingError <= tol_conv || (count > 2 && fabs(stoppingError - prevError1) <= tol_conv && fabs(prevError1 - prevError2) <= tol_conv))
      convergence = 1;
    else
      convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  // initialize perm & sv and see if they are in order
  perm = (int *)bmalloc(n * sizeof(int));
  sv = (double *)bmalloc(n * sizeof(double));
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = sv[j]
    sv[j] = 0;
    for (i = 0; i < n; i++)
      sv[j] += U->entry[i][j].r * U->entry[i][j].r + U->entry[i][j].i * U->entry[i][j].i;
    E->entry[j][j].r = sv[j] = sqrt(sv[j]);
    E->entry[j][j].i = 0;

    perm[j] = j;
    if (j > 0 && sv[j-1] < sv[j])
      order = 0;
  }

  if (!order)
  { // we need to find the order and update U & V

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (sv[j] > sv[i])
        { // swap j & i
          tempD = sv[j];
          sv[j] = sv[i];
          sv[i] = tempD;

          k = perm[j];
          perm[j] = perm[i];
          perm[i] = k;
        }

    // update U = U * P & V = V * P
    for (j = 0; j < n; j++)
    { // check to see if U & V need columns swapped
      if (perm[j] != j)
      {
        k = perm[j];
        // update E
        set_d(tempComp, &E->entry[j][j]);
        set_d(&E->entry[j][j], &E->entry[k][k]);
        set_d(&E->entry[k][k], tempComp);
  
        for (i = 0; i < n; i++)
        { // update U
          set_d(tempComp, &U->entry[i][j]);
          set_d(&U->entry[i][j], &U->entry[i][k]);
          set_d(&U->entry[i][k], tempComp);
          // update V
          set_d(tempComp, &V->entry[i][j]);
          set_d(&V->entry[i][j], &V->entry[i][k]);
          set_d(&V->entry[i][k], tempComp);
        }
        // update perm
        for (i = 0; i < n; i++)
          if (perm[i] == j)
          {
            perm[i] = perm[j];
            perm[j] = j;
            i = n;
          }
      }
    }
  }

  // update U - columns need normalized
  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  da = E->entry[0][0].r; // initialize to something
  for (j = 0; j < n; j++)
  {
    // the 2-norms of the jth column of U is E->entry[j][j].r
    db = da; // the previous singular value
    da = E->entry[j][j].r; // the current singular value

    if (checkGood_d(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column so we need to normalize this column of U
      count++; // this one is good!

      tempD = 1 / da;
      for (i = 0; i < n; i++)
      { // update U
        mul_rdouble_d(&U->entry[i][j], &U->entry[i][j], tempD);
      }
    }
    else
    { // column is the zero column - expand to unitary matrix
      expand_unitary_d(U, j, n);
      // exit loop
      j = n;
    }
  }

  free(perm);
  free(sv);

  // return the correct value based on convergence - either corank or -1
  if (convergence)
    return (n - count);
  else
    return -1;
}

int checkGood_d(double currSV, double prevSV, double tooSmall, double largeChange)
/***************************************************************\
* USAGE: checks to see if the current singular value should be  *
*   considered as acceptable                                    *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns 1 if good, otherwise 0                 *
* NOTES:                                                        *
\***************************************************************/
{
  if (currSV > tooSmall && prevSV < currSV * largeChange)
    return 1;
  else
    return 0;
}

//////// FIND SINGULAR VALUES & V ONLY ///////////

int svd_jacobi_EV_d(vec_d E, mat_d V, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: finds the E & V of SVD of InputMat using Jacobi its    *
* ARGUMENTS: InputMat = U*E*V^H - only return diag(E) & V       *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, k, retVal, m = InputMat->rows, n = InputMat->cols;
  int *perm = NULL;
  double norm, norm_inv;
  mat_d H, Q, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_vec_d(E, n);
    change_size_mat_d(V, n, n);
    E->size = V->rows = V->cols = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv_Jacobi < 0 || tol_QR < 0 || tol_sign < 0 || largeChange < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // for best results, we should normalize InputMat based on its infinity norm if it is > 1
  norm = infNormMat_d(InputMat);

  // determine how to proceed
  if (m >= n)
  { // we have more rows than columns so we just proceed
    init_mat_d(H, m, n); init_mat_d(R, m, n); init_mat_d(P, n, n);
    H->rows = m;
    H->cols = n;

    if (norm > 1)
    { // normalize InputMat & store to H
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          mul_rdouble_d(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // cp InputMat to H
      mat_cp_d(H, InputMat);
    }

    // compute H * P = Q * R, without computing Q - 1 at the end so that QR presorts the rows (P is n x n, Q would be m x m & R is m x n)
    QR_R_d(R, P, H, tol_QR, tol_sign, largeChange, 1);

    // compute R^T = H = U * E * V^T without computing V
    transpose_d(H, R);
    retVal = svd_R_jacobi_UE_d(R, E, H, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

    // so InputMat = (Q*V) * E * (P*U)^T, up to scaling factor
    // compute V = P * 'U' -> do this based on row swaps for efficiency
    // rescale E, if needed
    change_size_mat_d(V, n, n);
    V->rows = V->cols = n;
    perm = (int *)bmalloc(n * sizeof(int));
    convertToPerm_d(perm, P, n);
    for (i = 0; i < n; i++)
    { // setup row i
      k = perm[i];
      for (j = 0; j < n; j++)
      { // setup V[i][j]
        set_d(&V->entry[i][j], &R->entry[k][j]);
      }

      // rescale E, if needed
      if (norm > 1)
      { // undo the normalization by multiplying the entries of E by norm
        E->coord[i].r *= norm;
      }
    }

    // clear
    free(perm);
    clear_mat_d(H); clear_mat_d(R); clear_mat_d(P);
  }
  else // m < n
  { // we have more columns than rows so we must work with the transpose
    init_mat_d(H, n, m); init_mat_d(R, n, m); init_mat_d(P, m, m); init_mat_d(Q, n, n);
    H->rows = n;
    H->cols = m;

    if (norm > 1)
    { // normalize InputMat & store transpose to H
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          conjugate_d(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rdouble_d(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // H = InputMat^T
      transpose_d(H, InputMat);
    }

    // compute H * P = Q * R - 1 at the end so that QR presorts the rows (P is m x m, Q is n x n & R is n x m)
    QR_d(Q, R, P, H, tol_QR, tol_sign, largeChange, 1);

    // compute R = U * E * V^T without computing V
    retVal = svd_R_jacobi_UE_d(H, E, R, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

    // then InputMat = (P*V) * E * (Q*U)^T, up to scaling factor
    // compute V = Q * 'U'
    // rescale E, if needed
    mat_mul_d(V, Q, H);
    if (norm > 1)
    { // rescale E
      for (i = 0; i < m; i++)
      { // undo the normalization by multiplying the entries of E by norm
        E->coord[i].r *= norm;
      }
    }

    // clear
    clear_mat_d(H); clear_mat_d(Q); clear_mat_d(R); clear_mat_d(P);
  }

  return retVal;
}

int svd_R_jacobi_UE_d(mat_d U, vec_d E, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: finds the SVD of R - upper triangular matrix -  using  *
* Jacobi iterations                                             *
* ARGUMENTS: R = U*E*V^H, U - m x m, E - m x n, V would be n x n*
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, retVal, m = R->rows, n = R->cols;
  mat_d tempR;

  // do error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more rows than columns in svd_R_jacobi_UE_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  if (n == 0) // if R has no cols - immediately exit
  {
    change_size_mat_d(U, m, m);
    change_size_vec_d(E, n);
    U->rows = U->cols = m;
    E->size = n;
    make_matrix_ID_d(U, m, m);
    return 0;
  }

  // so m >= n >= 1
 
  // copy the main n x n part of R to tempR
  init_mat_d(tempR, n, n);
  tempR->rows = tempR->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_d(&tempR->entry[i][j], &R->entry[i][j]);
    }

  // do jacobi iterations on tempR
  retVal = jacobi_UE_d(U, E, tempR, its, rank_tol, tol_conv, tol_sign, largeChange);

  // adjust U & E based on m  - U = [[U, 0][0, I]] 
  increase_size_mat_d(U, m, m);
  U->rows = U->cols = m;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i < n)
      {
        if (j >= n)
        { // U[i][j] = 0
          set_zero_d(&U->entry[i][j]);
        }
      }
      else // i >= n
      {
        if (j < n)
        { // U[i][j] & E[i][j] = 0
          set_zero_d(&U->entry[i][j]);
        }
        else // j >= n
        { // U[i][j] is (i == j)
          set_double_d(&U->entry[i][j], i == j, 0);
        }
      }

  clear_mat_d(tempR);

  return retVal;
}

int jacobi_UE_d(mat_d U, vec_d E, mat_d R, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on R to find its singular values*
* and left singular values                                      *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: R has to be square                                     *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = R->rows, n = R->cols;
  int *perm = NULL;
  double da, db, sizeC_sqr, sizeC, stoppingError, prevError1, prevError2, tempD;
  double *sv = NULL;
  comp_d c, t, tempComp;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix is not square in jacobi_UE_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_UE_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv < 0 || tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_UE_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_mat_d(U, n, n);
    change_size_vec_d(E, n);
    U->rows = U->cols = E->size = n;
    return 0;
  }
  // so we can assume that m == n >= 1

  // copy R to U
  mat_cp_d(U, R);

  // main loop
  count = 0;
  stoppingError = prevError1 = prevError2 = 0; // initialize to something
  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    prevError2 = prevError1;
    prevError1 = stoppingError;
    stoppingError = 0; // reset stoppingError to 0
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        da = db = 0;
        set_zero_d(c);
        for (k = 0; k < n; k++)
        {
          da += U->entry[k][i].r * U->entry[k][i].r + U->entry[k][i].i * U->entry[k][i].i;
          db += U->entry[k][j].r * U->entry[k][j].r + U->entry[k][j].i * U->entry[k][j].i;
          conjugate_d(t, &U->entry[k][i]);
          sum_mul_d(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        sizeC_sqr = c->r * c->r + c->i * c->i;
        sizeC = sqrt(sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (da > tol_sign && db > tol_sign && sizeC > tol_sign)
        { // we can safely find the stopping condition estimate
          tempD = sizeC / (sqrt(da) * sqrt(db));
          if (tempD > stoppingError)
            stoppingError = tempD;
        }

        if (sizeC > tol_sign)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          db = (db - da) / 2;
          mul_rdouble_d(tempComp, c, db);
          db = 1 / sizeC_sqr;
          mul_rdouble_d(tempComp, tempComp, db);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_d(t, tempComp, tol_sign);
          db = tempComp->r * tempComp->r + tempComp->i * tempComp->i;
          da = sqrt(db);
          db = 1 / (da + sqrt(1 + db));
          mul_rdouble_d(t, t, db);
        }
        else
        {
          set_zero_d(t);
        }

        if (t->r != 0 || t->i != 0)
        { // we need to do the rotation
          // find db = 1 / sqrt(1 + |t|^2)
          db = 1 / sqrt(1 + t->r * t->r + t->i * t->i);
          // find t = db * t
          mul_rdouble_d(t, t, db);

          // tempComp = -conj(t)
          tempComp->r = -t->r;
          tempComp->i = t->i;

          // update columns i & j of U
          for (k = 0; k < n; k++)
          { // update U[k][i] & U[k][j]
            set_d(c, &U->entry[k][i]);

            mul_rdouble_d(&U->entry[k][i], c, db);
            sum_mul_d(&U->entry[k][i], tempComp, &U->entry[k][j]);

            mul_rdouble_d(&U->entry[k][j], &U->entry[k][j], db);
            sum_mul_d(&U->entry[k][j], c, t);
          }
        }
      }

    // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
    if (stoppingError <= tol_conv || (count > 2 && fabs(stoppingError - prevError1) <= tol_conv && fabs(prevError1 - prevError2) <= tol_conv))
      convergence = 1;
    else
      convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  // initialize perm & sv and see if they are in order
  perm = (int *)bmalloc(n * sizeof(int));
  sv = (double *)bmalloc(n * sizeof(double));
  change_size_vec_d(E, n);
  E->size = n;
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = sv[j]
    sv[j] = 0;
    for (i = 0; i < n; i++)
      sv[j] += U->entry[i][j].r * U->entry[i][j].r + U->entry[i][j].i * U->entry[i][j].i;
    E->coord[j].r = sv[j] = sqrt(sv[j]);
    E->coord[j].i = 0;

    perm[j] = j;
    if (j > 0 && sv[j-1] < sv[j])
      order = 0;
  }

  if (!order)
  { // we need to find the order and update U

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (sv[j] > sv[i])
        { // swap j & i
          tempD = sv[j];
          sv[j] = sv[i];
          sv[i] = tempD;

          k = perm[j];
          perm[j] = perm[i];
          perm[i] = k;
        }

    // update U = U * P
    for (j = 0; j < n; j++)
    { // check to see if U need columns swapped
      if (perm[j] != j)
      {
        k = perm[j];
        // update E
        set_d(tempComp, &E->coord[j]);
        set_d(&E->coord[j], &E->coord[k]);
        set_d(&E->coord[k], tempComp);
  
        for (i = 0; i < n; i++)
        { // update U
          set_d(tempComp, &U->entry[i][j]);
          set_d(&U->entry[i][j], &U->entry[i][k]);
          set_d(&U->entry[i][k], tempComp);
        }
        // update perm
        for (i = 0; i < n; i++)
          if (perm[i] == j)
          {
            perm[i] = perm[j];
            perm[j] = j;
            i = n;
          }
      }
    }
  }

  // update U - columns need normalized
  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  da = E->coord[0].r; // initialize to something
  for (j = 0; j < n; j++)
  { // the 2-norms of the jth column of U is E->coord[j].r
    db = da; // the previous singular value
    da = E->coord[j].r; // the current singular value

    if (checkGood_d(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column so we need to normalize this column of U
      count++; // this one is good!

      tempD = 1 / da;
      for (i = 0; i < n; i++)
      { // update U
        mul_rdouble_d(&U->entry[i][j], &U->entry[i][j], tempD);
      }
    }
    else
    { // column is the zero column - expand to unitary matrix
      expand_unitary_d(U, j, n);
      // exit loop
      j = n;
    }
  }

  free(perm);
  free(sv);

  // return the correct value based on convergence - either corank or -1
  if (convergence)
    return (n - count);
  else
    return -1;
}

//////// FIND SINGULAR VALUES ONLY //////////

int svd_jacobi_E_d(vec_d E, mat_d InputMat, int its, double rank_tol, double tol_conv_Jacobi, double tol_QR, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: finds the singular values of InputMat using Jacobi its *
* ARGUMENTS: InputMat = U*E*V^H - only return diag(E)           *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the singular values done          *
\***************************************************************/
{
  int i, j, retVal, m = InputMat->rows, n = InputMat->cols;
  double norm, norm_inv;
  mat_d H, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_vec_d(E, n);
    E->size = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv_Jacobi < 0 || tol_QR < 0 || tol_sign < 0 || largeChange < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  init_mat_d(H, 0, 0); init_mat_d(R, 0, 0); init_mat_d(P, 0, 0);

  // for best results, we should normalize InputMat based on its infinity norm if it is > 1
  norm = infNormMat_d(InputMat);

  // normalize and make 'skinny'
  if (m < n)
  { // setup H = (normalized) InputMat^T
    change_size_mat_d(H, n, m);
    H->rows = n;
    H->cols = m;
    if (norm > 1)
    { // need to normalize and transpose
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        { // normalize and transpose
          conjugate_d(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rdouble_d(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // just transpose
      transpose_d(H, InputMat);
    }
    // setup m & n
    m = H->rows;
    n = H->cols;
  }
  else
  { // setup H = (normalized) InputMat
    change_size_mat_d(H, m, n);
    H->rows = m;
    H->cols = n;
    if (norm > 1)
    { // need to normalize
      norm_inv = 1 / norm;
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        { // normalize
          mul_rdouble_d(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // just copy
      mat_cp_d(H, InputMat);
    }
  }
  // so we have m >= n >= 1

  // do the QR decomposition - 1 at the end so that QR presorts the rows ( H * P = Q * R, where P is n x n, Q would be m x m & R is m x n)
  QR_R_d(R, P, H, tol_QR, tol_sign, largeChange, 1);

  // since m >= n, we want to chop off the zeros on the bottom of R and then run the SVD on R^T
  R->rows = R->cols = n;
  transpose_d(H, R);

  // find the singular values of H using jacobi iterations
  // NOTE: jacobi will always return its approximation of the svd even if it does not converge within the desired number of iterations
  retVal = jacobi_E_d(E, H, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

  if (norm > 1)
  { // undo the normalization by multiplying the entries of E by norm
    for (i = 0; i < n; i++)
      E->coord[i].r *= norm;
  }

  // clear
  clear_mat_d(H); clear_mat_d(R); clear_mat_d(P);

  return retVal;
}

int jacobi_E_d(vec_d E, mat_d U, int its, double rank_tol, double tol_conv, double tol_sign, double largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on U to find its singular values*
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: U has to be square                                     *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = U->rows, n = U->cols;
  double da, db, sizeC_sqr, sizeC, stoppingError, prevError1, prevError2, tempD;
  comp_d c, t, tempComp;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix is not square in jacobi_E_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_E_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (rank_tol < 0 || tol_conv < 0 || tol_sign < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_E_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_vec_d(E, n);
    E->size = n;
    return 0;
  }
  // so we can assume that m == n >= 1

  // main loop
  count = 0;
  stoppingError = prevError1 = prevError2 = 0; // initialize to something
  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    prevError2 = prevError1;
    prevError1 = stoppingError;
    stoppingError = 0; // reset stoppingError to 0
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        da = db = 0;
        set_zero_d(c);
        for (k = 0; k < n; k++)
        {
          da += U->entry[k][i].r * U->entry[k][i].r + U->entry[k][i].i * U->entry[k][i].i;
          db += U->entry[k][j].r * U->entry[k][j].r + U->entry[k][j].i * U->entry[k][j].i;
          conjugate_d(t, &U->entry[k][i]);
          sum_mul_d(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        sizeC_sqr = c->r * c->r + c->i * c->i;
        sizeC = sqrt(sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (da > tol_sign && db > tol_sign && sizeC > tol_sign)
        { // we can safely find the stopping condition estimate
          tempD = sizeC / (sqrt(da) * sqrt(db));
          if (tempD > stoppingError)
            stoppingError = tempD;
        }

        if (sizeC > tol_sign)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          db = (db - da) / 2;
          mul_rdouble_d(tempComp, c, db);
          db = 1 / sizeC_sqr;
          mul_rdouble_d(tempComp, tempComp, db);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_d(t, tempComp, tol_sign);
          db = tempComp->r * tempComp->r + tempComp->i * tempComp->i;
          da = sqrt(db);
          db = 1 / (da + sqrt(1 + db));
          mul_rdouble_d(t, t, db);
        }
        else
        {
          set_zero_d(t);
        }

        if (t->r != 0 || t->i != 0)
        { // we need to do the rotation
          // find db = 1 / sqrt(1 + |t|^2)
          db = 1 / sqrt(1 + t->r * t->r + t->i * t->i);
          // find t = db * t
          mul_rdouble_d(t, t, db);

          // tempComp = -conj(t)
          tempComp->r = -t->r;
          tempComp->i = t->i;

          // update columns i & j of U
          for (k = 0; k < n; k++)
          { // update U[k][i] & U[k][j]
            set_d(c, &U->entry[k][i]);

            mul_rdouble_d(&U->entry[k][i], c, db);
            sum_mul_d(&U->entry[k][i], tempComp, &U->entry[k][j]);

            mul_rdouble_d(&U->entry[k][j], &U->entry[k][j], db);
            sum_mul_d(&U->entry[k][j], c, t);
          }
        }
      }

    // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
    if (stoppingError <= tol_conv || (count > 2 && fabs(stoppingError - prevError1) <= tol_conv && fabs(prevError1 - prevError2) <= tol_conv))
      convergence = 1;
    else
      convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  change_size_vec_d(E, n);
  E->size = n;
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = E[j]
    set_zero_d(&E->coord[j]);
    for (i = 0; i < n; i++)
      E->coord[j].r += U->entry[i][j].r * U->entry[i][j].r + U->entry[i][j].i * U->entry[i][j].i;
    E->coord[j].r = sqrt(E->coord[j].r);

    if (j > 0 && E->coord[j-1].r < E->coord[j].r)
      order = 0;
  }

  if (!order)
  { // we need to find the order

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (E->coord[j].r > E->coord[i].r)
        { // swap j & i
          tempD = E->coord[j].r;
          E->coord[j].r = E->coord[i].r;
          E->coord[i].r = tempD;
        }
  }

  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  da = E->coord[0].r; // initialize to something
  for (j = 0; j < n; j++)
  { // the 2-norms of the jth column of U is E->entry[j][j].r
    db = da; // the previous singular value
    da = E->coord[j].r; // the current singular value

    if (checkGood_d(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column so we need to normalize this column of U
      count++; // this one is good!
    }
    else
    { // exit loop
      j = n;
    }
  }

  // return the correct value based on convergence - either corank or -1
  if (convergence)
    return (n - count);
  else
    return -1;
}

//////// EXPAND TO A UNITARY MATRIX ///////////

void expand_unitary_d(mat_d U, int startCol, int endCol)
/***************************************************************\
* USAGE: puts in random orthonormal columns                     *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, m = U->rows, n = U->cols;
  double norm = 0;
  comp_d proj, temp;

  // error checking - make sure 0 <= startCol < endCol <= m
  if (!(0 <= startCol && startCol < endCol && endCol <= m))
  {
    printf("ERROR: The columns are not correct in expand_unitary_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m < n) // need to make sure that we have a square or tall matrix (otherwise cannot make lin. indep.)
  {
    printf("ERROR: The matrix needs to have atleast as many rows as columns in expand_unitary_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
 
  for (j = startCol; j < endCol; j++)
  { // create a random column
    for (i = 0; i < m; i++)
      get_comp_rand_d(&U->entry[i][j]);

    // loop over the other columns to create an orthogonal vector
    for (k = 0; k < j; k++)
    { // set proj to 0
      set_zero_d(proj);
      // find U[:][k]^T * U[:][j]
      for (i = 0; i < m; i++)
      { // proj += conj(U[i][k]) * U[i][j]
        conjugate_d(temp, &U->entry[i][k]);
        sum_mul_d(proj, temp, &U->entry[i][j]);
      }

      // negate proj
      neg_d(proj, proj);

      // update U[:][j] += proj * U[:][k] & find its norm
      for (i = 0; i < m; i++)
      {
        sum_mul_d(&U->entry[i][j], proj, &U->entry[i][k]);
      }
    }

    // find the norm
    norm = 0;
    for (i = 0; i < m; i++)
      norm += U->entry[i][j].r * U->entry[i][j].r + U->entry[i][j].i * U->entry[i][j].i;
    norm = sqrt(norm);

    // invert norm
    norm = 1 / norm;

    // normalize U[:][j]
    for (i = 0; i < m; i++)
    {
      mul_rdouble_d(&U->entry[i][j], &U->entry[i][j], norm);
    }
  }

  return;
}

//////// MULTI PRECISION ///////////

int svd_corank_mp(mat_mp A, mpf_t rank_tol, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the corank of A using its SVD                    *
* ARGUMENTS:                                                    *
* RETURN VALUES: the corank of A if successful, else -1 if not  *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal, its = 75, num_digits = -prec_to_digits(mpf_get_prec(A->entry[0][0].r));  // Computes the number of digits being used
  // QR does not need to have too tight of a tolerance
  // while tol_sign needs to be relatively small for the best results
  char *str = NULL;
  size_t size;
  mpf_t tol_prec, tol_sign;
  vec_mp E;

  // initialize MP
  mpf_init(tol_prec); mpf_init(tol_sign);
  init_vec_mp(E, 0);

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
 
  // check tolerances
  if (mpf_cmp(rank_tol, tol_prec) < 0)
  { // ideally, we would like tol_prec a couple orders of magnitude smaller than rank_tol
    printf("WARNING: The tolerance to determine the rank is smaller than the precision tolerance for SVD!\nBy default, the rank tolerance will be adjusted.\n");
    mpf_set(rank_tol, tol_prec);
  }

  // try to find the singular values with the tolerances
  retVal = svd_jacobi_E_mp(E, A, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);

  // check to see if successful
  if (retVal < 0)
  { // SVD was not successful
    // adjust tolerances to see if we can get convergence - this should very rarely happen!
    its = 100;
    // loosen convergence criterion tol_Jacobi by making larger and loosen tol_sign by making smaller
    mpf_mul_ui(tol_prec, tol_prec, 10);
    mpf_div_ui(tol_sign, tol_sign, 100);

    retVal = svd_jacobi_E_mp(E, A, its, rank_tol, tol_prec, tol_prec, tol_sign, largeChange);
  }

  // clear MP
  mpf_clear(tol_prec); mpf_clear(tol_sign);
  clear_vec_mp(E);

  // free str
  free(str);

  return retVal;
}

int svd_corank_analyze_mp(mat_mp E, mpf_t rank_tol, mpf_t largeChange)
/***************************************************************\
* USAGE: analyzes E-diagonal matrix from SVD-to find its corank *
* ARGUMENTS:                                                    *
* RETURN VALUES: the corank of E                                *
* NOTES: Assume diagonal entries of E are decreasing,real & >= 0*
\***************************************************************/
{
  int i, rank, m = E->rows, n = E->cols;

  // do error checking
  if (rank_tol < 0 || largeChange < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_rank_analyze_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // set n = min(m, n) so that E has n diagonal entries to look at
  if (m < n)
    n = m;

  if (n == 0) // no singular values has rank 0
    return 0;
  // so we can assume that n > 0

  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  rank = 0;
  // check the first one
  if (checkGood_mp(E->entry[0][0].r, E->entry[0][0].r, rank_tol, largeChange))
  { // the first one is good
    rank++;

    // check the rest of them
    for (i = 1; i < n; i++)
    {
      if (checkGood_mp(E->entry[i][i].r, E->entry[i-1][i-1].r, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
      { // this one is good!
        rank++; // this one is good!
      }
      else
      { // this one is no good!
        i = n; // all the rest of them are not good as well
      }
    }
  }

  return (n - rank);
}

int svd_jacobi_mp_prec(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, double rankTol, int curr_prec)
/***************************************************************\
* USAGE: finds the SVD of InputMat using Jacobi iterations      *
* ARGUMENTS: InputMat = U*E*V^H                                 *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int retVal, its = 100, num_digits = prec_to_digits(curr_prec);
  size_t size;
  char *str = NULL;
  mpf_t rankTol_mp, tol_conv_mp, tol_QR_mp, tol_sign_mp, largeChange_mp; 

  mpf_init(rankTol_mp);
  mpf_init(tol_conv_mp);
  mpf_init(tol_QR_mp);
  mpf_init(tol_sign_mp);
  mpf_init(largeChange_mp);

  // setup rankTol
  mpf_set_d(rankTol_mp, rankTol);
  // setup tol_conv
  size = 1 + snprintf(NULL, 0, "1e-%d", num_digits-3);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", num_digits-3);
  mpf_set_str(tol_conv_mp, str, 10);
  // setup tol_pivot
  size = 1 + snprintf(NULL, 0, "1e-%d", num_digits-1);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", num_digits-1);
  mpf_set_str(tol_QR_mp, str, 10);
  // setup tol_sign
  size = 1 + snprintf(NULL, 0, "1e-%d", 2 * num_digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", 2 * num_digits);
  mpf_set_str(tol_sign_mp, str, 10);
  // setup largeChange
  size = 1 + snprintf(NULL, 0, "1e%d", num_digits-2);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", num_digits-2);
  mpf_set_str(largeChange_mp, str, 10);

  retVal = svd_jacobi_mp(U, E, V, InputMat, its, rankTol_mp, tol_conv_mp, tol_QR_mp, tol_sign_mp, largeChange_mp);

  mpf_clear(rankTol_mp);
  mpf_clear(tol_conv_mp);
  mpf_clear(tol_QR_mp);
  mpf_clear(tol_sign_mp);
  mpf_clear(largeChange_mp);
  free(str);

  return retVal;
}

int svd_jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the SVD of InputMat using Jacobi iterations      *
* ARGUMENTS: InputMat = U*E*V^H                                 *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, k, retVal, transposed, m = InputMat->rows, n = InputMat->cols;
  int *perm = NULL;
  mpf_t norm, norm_inv, dr, di;
  comp_mp tempComp;
  mat_mp H, Q, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_mat_mp(U, m, n);
    change_size_mat_mp(E, n, n);
    change_size_mat_mp(V, n, n);
    U->rows = m;
    U->cols = E->rows = E->cols = V->rows = V->cols = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(tol_conv_Jacobi, 0) < 0 || mpf_cmp_ui(tol_QR, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0 || mpf_cmp_ui(largeChange, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize
  mpf_init(norm); mpf_init(norm_inv); mpf_init(dr); mpf_init(di);
  init_mp(tempComp);

  // for best results, we should normalize H based on its infinity norm if it is > 1
  infNormMat_mp2(norm, InputMat);

  // check to make sure that InputMat is 'skinny'
  if (m >= n)
  { // we have more rows than columns so we just proceed
    init_mat_mp(H, m, n); init_mat_mp(Q, m, m); init_mat_mp(R, m, n); init_mat_mp(P, n, n);
    H->rows = m;
    H->cols = n;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // normalize InputMat and store to H
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          mul_rmpf_mp(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // cp InputMat to H
      mat_cp_mp(H, InputMat);
    }

    // did not transpose the matrix
    transposed = 0;
  }
  else
  { // we have more columns than rows so we must work with the transpose
    init_mat_mp(H, n, m); init_mat_mp(Q, n, n); init_mat_mp(R, n, m); init_mat_mp(P, m, m);
    H->rows = n;
    H->cols = m;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // normalize InputMat and store transpose to H
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          conjugate_mp(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rmpf_mp(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // H = InputMat^T
      transpose_mp(H, InputMat);
    }

    // transposed the matrix
    transposed = 1;
    m = H->rows;
    n = H->cols;
  }
  // so we have m >= n >= 1

  // do the QR decomposition - 1 at the end so that QR presorts the rows ( H * P = Q * R )
  QR_mp(Q, R, P, H, tol_QR, tol_sign, largeChange, 1);

  // decompose R using jacobi_mp ( R = U * E * V^H, U is m x m, E is m x n and V is n x n)
  // NOTE: jacobi will always return its approximation of the svd even if it does not converge within the desired number of iterations
  retVal = svd_R_jacobi_mp(U, E, V, R, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

  if (mpf_cmp_ui(norm, 1) > 0)
  { // undo the normalization by multiplying the entries of E by norm
    for (i = 0; i < n; i++)
    {
      mpf_mul(E->entry[i][i].r, E->entry[i][i].r, norm);
    }
  }

  // update U based on Q - U = Q * U
  mat_mul_mp(U, Q, U);

  // update V based on P - V = P * V - do this based no row swaps for efficiency
  perm = (int *)bmalloc(n * sizeof(int));
  convertToPerm_mp(perm, P, n);
  for (i = 0; i < n; i++)
    if (perm[i] != i)
    { // swap i & perm[i]
      k = perm[i];
      for (j = 0; j < n; j++)
      { // update V
        set_mp(tempComp, &V->entry[i][j]);
        set_mp(&V->entry[i][j], &V->entry[k][j]);
        set_mp(&V->entry[k][j], tempComp);
      }
      // update perm
      for (j = 0; j < n; j++)
        if (perm[j] == i)
        {
          perm[j] = perm[i];
          perm[i] = i;
          j = n;
        }
    }
  // so we have H = U * E * V^H, thus we need to check if H = InputMat or H = InputMat^H

  if (transposed)
  { // (U * E * V^H) ^ H = (V * E^H * U^H) = V * E * U^H, so we just need to swap V & U
    mat_cp_mp(H, V);
    mat_cp_mp(V, U);
    mat_cp_mp(U, H);
    transpose_mp(E, E);
  }

  // clear
  free(perm);
  mpf_clear(norm); mpf_clear(norm_inv); mpf_clear(dr); mpf_clear(di);
  clear_mp(tempComp);
  clear_mat_mp(H); clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);

  return retVal;
}

int svd_R_jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the SVD of R - upper triangular matrix -  using  *
* Jacobi iterations                                             *
* ARGUMENTS: R = U*E*V^H, U - m x m, E - m x n, V - n x n       *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, retVal, m = R->rows, n = R->cols;
  mat_mp tempR;

  // do error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more rows than columns in svd_R_jacobi_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  if (n == 0) // if R has no cols - immediately exit
  {
    change_size_mat_mp(U, m, m);
    change_size_mat_mp(E, m, n);
    change_size_mat_mp(V, n, n);
    U->rows = U->cols = E->rows = m;
    E->cols = V->rows = V->cols = n;
    make_matrix_ID_mp(U, m, m);
    return 0;
  }

  // so m >= n >= 1

  // copy the main n x n part of R to tempR
  init_mat_mp(tempR, n, n);
  tempR->rows = tempR->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_mp(&tempR->entry[i][j], &R->entry[i][j]);
    }

  // do jacobi iterations on tempR
  retVal = jacobi_mp(U, E, V, tempR, its, rank_tol, tol_conv, tol_sign, largeChange);

  // adjust U & E based on m  - U = [[U, 0][0, I]] and E = [[E][0]]
  increase_size_mat_mp(U, m, m);
  increase_size_mat_mp(E, m, E->cols);
  U->rows = U->cols = E->rows = m;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i < n)
      {
        if (j >= n)
        { // U[i][j] = 0
          set_zero_mp(&U->entry[i][j]);
        }
      }
      else // i >= n
      {
        if (j < n)
        { // U[i][j] & E[i][j] = 0
          set_zero_mp(&U->entry[i][j]);
          set_zero_mp(&E->entry[i][j]);
        }
        else // j >= n
        { // U[i][j] is (i == j)
          mpf_set_ui(U->entry[i][j].r, i == j);
          mpf_set_ui(U->entry[i][j].i, 0);
        }
      }

  clear_mat_mp(tempR);

  return retVal;
}

int jacobi_mp(mat_mp U, mat_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on InputMat to find its SVD     *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: InputMat has to be square - should be R of a QR decomp *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = InputMat->rows, n = InputMat->cols;
  mpf_t stoppingError, prevError1, prevError2;
  mpf_t da, db, dr, di, sizeC_sqr, sizeC;
  comp_mp c, t, tempComp;

  int *perm = NULL;
  mpf_t *sv = NULL;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix has needs to be square in jacobi_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(rank_tol, 0) < 0 || mpf_cmp_ui(tol_conv, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_mat_mp(U, n, n);
    change_size_mat_mp(E, n, n);
    change_size_mat_mp(V, n, n);
    U->rows = U->cols = E->rows = E->cols = V->rows = V->cols = n;
    return 0;
  }
  // so we can assume that m == n > 0

  // initialize MP
  mpf_init(stoppingError); mpf_init(prevError1); mpf_init(prevError2);
  mpf_init(da); mpf_init(db); mpf_init(dr); mpf_init(di); mpf_init(sizeC_sqr); mpf_init(sizeC);
  init_mp(c); init_mp(t); init_mp(tempComp);

  // initialiaze U to InputMat, E to 0 and V to Id
  change_size_mat_mp(U, n, n);
  change_size_mat_mp(E, n, n);
  change_size_mat_mp(V, n, n);
  U->rows = U->cols = E->rows = E->cols = V->rows = V->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    { // U
      set_mp(&U->entry[i][j], &InputMat->entry[i][j]);
      // E
      set_zero_mp(&E->entry[i][j]);
      // V
      set_double_mp(&V->entry[i][j], i == j, 0);
    }

  // main loop
  count = 0;
  // initialize errors to something
  mpf_set_ui(stoppingError, 0);
  mpf_set_ui(prevError1, 0);
  mpf_set_ui(prevError2, 0);

  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    mpf_set(prevError2, prevError1);
    mpf_set(prevError1, stoppingError);
    // reset stoppingError to 0
    mpf_set_ui(stoppingError, 0);
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for(i = 0; i < n-1; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        mpf_set_ui(da, 0);
        mpf_set_ui(db, 0);
        set_zero_mp(c);
        for (k = 0; k < n; k++)
        {
          mpf_mul(dr, U->entry[k][i].r, U->entry[k][i].r);
          mpf_mul(di, U->entry[k][i].i, U->entry[k][i].i);
          mpf_add(dr, dr, di);
          mpf_add(da, da, dr);
          mpf_mul(dr, U->entry[k][j].r, U->entry[k][j].r);
          mpf_mul(di, U->entry[k][j].i, U->entry[k][j].i);
          mpf_add(dr, dr, di);
          mpf_add(db, db, dr);
          conjugate_mp(t, &U->entry[k][i]);
          sum_mul_mp(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        mpf_mul(dr, c->r, c->r);
        mpf_mul(di, c->i, c->i);
        mpf_add(sizeC_sqr, dr, di);
        mpf_sqrt(sizeC, sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (mpf_cmp(da, tol_sign) > 0 && mpf_cmp(db, tol_sign) > 0 && mpf_cmp(sizeC, tol_sign) > 0)
        { // we can safely find the stopping condition estimate - sizeC / sqrt(da * db)
          mpf_mul(dr, da, db);
          mpf_sqrt(dr, dr);
          mpf_div(dr, sizeC, dr);
          if (mpf_cmp(dr, stoppingError) > 0)
            mpf_set(stoppingError, dr);
        }

        if (mpf_cmp(sizeC, tol_sign) > 0)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          mpf_mul_ui(dr, sizeC_sqr, 2);
          mpf_sub(di, db, da);
          mpf_div(db, di, dr);
          mpf_mul(tempComp->r, c->r, db);
          mpf_mul(tempComp->i, c->i, db);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_mp2(t, tempComp, tol_sign);
          mpf_mul(dr, tempComp->r, tempComp->r);
          mpf_mul(di, tempComp->i, tempComp->i);
          mpf_add(db, dr, di);
          mpf_sqrt(da, db);
          mpf_add_ui(db, db, 1);
          mpf_sqrt(db, db);
          mpf_add(db, da, db);
          mpf_ui_div(db, 1, db);

          mpf_mul(t->r, t->r, db);
          mpf_mul(t->i, t->i, db);
        }
        else
        {
          set_zero_mp(t);
        }

        // find db = 1 / sqrt(1 + |t|^2)
        mpf_mul(dr, t->r, t->r);
        mpf_mul(di, t->i, t->i);
        mpf_add(db, dr, di);
        mpf_add_ui(db, db, 1);
        mpf_sqrt(db, db);
        mpf_ui_div(db, 1, db);

        // find t = db * t
        mpf_mul(t->r, t->r, db);
        mpf_mul(t->i, t->i, db);

        // tempComp = -conj(t)
        mpf_neg(tempComp->r, t->r);
        mpf_set(tempComp->i, t->i);

        // update columns i & j of U & V
        for (k = 0; k < n; k++)
        { // update U
          set_mp(c, &U->entry[k][i]);

          mpf_mul(U->entry[k][i].r, c->r, db);
          mpf_mul(U->entry[k][i].i, c->i, db);
          sum_mul_mp(&U->entry[k][i], tempComp, &U->entry[k][j]);

          mpf_mul(U->entry[k][j].r, U->entry[k][j].r, db);
          mpf_mul(U->entry[k][j].i, U->entry[k][j].i, db);
          sum_mul_mp(&U->entry[k][j], c, t);

          // update V
          set_mp(c, &V->entry[k][i]);

          mpf_mul(V->entry[k][i].r, c->r, db);
          mpf_mul(V->entry[k][i].i, c->i, db);
          sum_mul_mp(&V->entry[k][i], tempComp, &V->entry[k][j]);

          mpf_mul(V->entry[k][j].r, V->entry[k][j].r, db);
          mpf_mul(V->entry[k][j].i, V->entry[k][j].i, db);
          sum_mul_mp(&V->entry[k][j], c, t);
        }
      }

      // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
      mpf_sub(da, stoppingError, prevError1);
      mpf_abs(da, da);
      mpf_sub(db, prevError1, prevError2);
      mpf_abs(db, db);
      if (mpf_cmp(stoppingError, tol_conv) <= 0 || (count > 2 && mpf_cmp(da, tol_conv) <= 0 && mpf_cmp(db, tol_conv) <= 0))
        convergence = 1;
      else
        convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  // initialize perm & sv and see if they are in order
  perm = (int *)bmalloc(n * sizeof(int));
  sv = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = sv[j]
    mpf_set_ui(da, 0);
    for (i = 0; i < n; i++)
    {
      mpf_mul(dr, U->entry[i][j].r, U->entry[i][j].r);
      mpf_mul(di, U->entry[i][j].i, U->entry[i][j].i);
      mpf_add(dr, dr, di);
      mpf_add(da, da, dr);
    }
    mpf_sqrt(E->entry[j][j].r, da);
    mpf_set_ui(E->entry[j][j].i, 0);

    mpf_init_set(sv[j], E->entry[j][j].r);

    perm[j] = j;
    if (j > 0 && mpf_cmp(sv[j-1], sv[j]) < 0)
      order = 0;
  }

  if (!order)
  { // we need to find the order and update U, E & V

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (mpf_cmp(sv[j], sv[i]) > 0)
        { // swap j & i
          mpf_set(da, sv[j]);
          mpf_set(sv[j], sv[i]);
          mpf_set(sv[i], da);

          k = perm[j];
          perm[j] = perm[i];
          perm[i] = k;
        }

    // update U = U * P & V = V * P
    for (j = 0; j < n; j++)
    { // check to see if U & V need columns swapped
      if (perm[j] != j)
      {
        k = perm[j];
        // update E
        set_mp(tempComp, &E->entry[j][j]);
        set_mp(&E->entry[j][j], &E->entry[k][k]);
        set_mp(&E->entry[k][k], tempComp);

        for (i = 0; i < n; i++)
        { // update U
          set_mp(tempComp, &U->entry[i][j]);
          set_mp(&U->entry[i][j], &U->entry[i][k]);
          set_mp(&U->entry[i][k], tempComp);
          // update V
          set_mp(tempComp, &V->entry[i][j]);
          set_mp(&V->entry[i][j], &V->entry[i][k]);
          set_mp(&V->entry[i][k], tempComp);
        }
        // update perm
        for (i = 0; i < n; i++)
          if (perm[i] == j)
          {
            perm[i] = perm[j];
            perm[j] = j;
            i = n;
          }
      }
    }
  }

  // update U - columns need normalized
  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  mpf_set(da, E->entry[0][0].r); // initialize to something
  for (j = 0; j < n; j++)
  { 
    // the 2-norms of the jth column of U is E->entry[j][j].r
    mpf_set(db, da); // the previous singular value
    mpf_set(da, E->entry[j][j].r); // the current singular value

    if (checkGood_mp(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column so we need to normalize this column of U
      count++; // this one is good!

      mpf_ui_div(dr, 1, da);
      for (i = 0; i < n; i++)
      { // update U
        mpf_mul(U->entry[i][j].r, U->entry[i][j].r, dr);
        mpf_mul(U->entry[i][j].i, U->entry[i][j].i, dr);
      }
    }
    else
    { // column is the zero column - expand to unitary matrix
      expand_unitary_mp(U, j, n);
      // exit loop
      j = n;
    }
  }

  // clear MP
  for (j = n - 1; j >= 0; j--)
    mpf_clear(sv[j]);
  mpf_clear(stoppingError); mpf_clear(prevError1); mpf_clear(prevError2);
  mpf_clear(da); mpf_clear(db); mpf_clear(dr); mpf_clear(di); mpf_clear(sizeC_sqr); mpf_clear(sizeC);
  clear_mp(c); clear_mp(t); clear_mp(tempComp);

  // free memory
  free(sv);
  free(perm);

  // return the correct value based on convergence - either corank or -1 
  if (convergence)
    return (n - count);
  else
    return -1;
}

int checkGood_mp(mpf_t currSV, mpf_t prevSV, mpf_t tooSmall, mpf_t largeChange)
/***************************************************************\
* USAGE: checks to see if the current singular value should be  *
*   considered as acceptable                                    *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns 1 if good, otherwise 0                 *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal;
  mpf_t tempMPF;
  mpf_init(tempMPF);

  mpf_mul(tempMPF, currSV, largeChange);

  if (mpf_cmp(currSV, tooSmall) > 0 && mpf_cmp(prevSV, tempMPF) < 0)
    retVal = 1;
  else
    retVal = 0;

  mpf_clear(tempMPF);

  return retVal;
}

//////// FIND SINGULAR VALUES & V ONLY ///////////

int svd_jacobi_EV_mp(vec_mp E, mat_mp V, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the E & V of SVD of InputMat using Jacobi its    *
* ARGUMENTS: InputMat = U*E*V^H - only return diag(E) & V       *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, k, retVal, m = InputMat->rows, n = InputMat->cols;
  int *perm = NULL;
  mpf_t norm, norm_inv;
  mat_mp H, Q, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_vec_mp(E, n);
    change_size_mat_mp(V, n, n);
    E->size = V->rows = V->cols = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(tol_conv_Jacobi, 0) < 0 || mpf_cmp_ui(tol_QR, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0 || mpf_cmp_ui(largeChange, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize
  mpf_init(norm); mpf_init(norm_inv);

  // for best results, we should normalize H based on its infinity norm if it is > 1
  infNormMat_mp2(norm, InputMat);

  // determine how to proceed
  if (m >= n)
  { // we have more rows than columns so we just proceed
    init_mat_mp(H, m, n); init_mat_mp(R, m, n); init_mat_mp(P, n, n);
    H->rows = m;
    H->cols = n;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // normalize InputMat & store to H
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          mul_rmpf_mp(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // cp InputMat to H
      mat_cp_mp(H, InputMat);
    }

    // compute H * P = Q * R, without computing Q - 1 at the end so that QR presorts the rows (P is n x n, Q would be m x m & R is m x n)
    QR_R_mp(R, P, H, tol_QR, tol_sign, largeChange, 1);

    // compute R^T = H = U * E * V^T without computing V
    transpose_mp(H, R);
    retVal = svd_R_jacobi_UE_mp(R, E, H, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

    // so InputMat = (Q*V) * E * (P*U)^T, up to scaling factor
    // compute V = P * 'U' -> do this based on row swaps for efficiency
    // rescale E, if needed
    change_size_mat_mp(V, n, n);
    V->rows = V->cols = n;
    perm = (int *)bmalloc(n * sizeof(int));
    convertToPerm_mp(perm, P, n);
    for (i = 0; i < n; i++)
    { // setup row i
      k = perm[i];
      for (j = 0; j < n; j++)
      { // setup V[i][j]
        set_mp(&V->entry[i][j], &R->entry[k][j]);
      }

      // rescale E, if needed
      if (mpf_cmp_ui(norm, 1) > 0)
      { // undo the normalization by multiplying the entries of E by norm
        mpf_mul(E->coord[i].r, E->coord[i].r, norm);
      }
    }

    // clear
    free(perm);
    clear_mat_mp(H); clear_mat_mp(R); clear_mat_mp(P);
  }
  else // m < n
  { // we have more columns than rows so we must work with the transpose
    init_mat_mp(H, n, m); init_mat_mp(R, n, m); init_mat_mp(P, m, m); init_mat_mp(Q, n, n);
    H->rows = n;
    H->cols = m;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // normalize InputMat & store transpose to H
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        {
          conjugate_mp(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rmpf_mp(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // H = InputMat^T
      transpose_mp(H, InputMat);
    }

    // compute H * P = Q * R - 1 at the end so that QR presorts the rows (P is m x m, Q is n x n & R is n x m)
    QR_mp(Q, R, P, H, tol_QR, tol_sign, largeChange, 1);

    // compute R = U * E * V^T without computing V
    retVal = svd_R_jacobi_UE_mp(H, E, R, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

    // then InputMat = (P*V) * E * (Q*U)^T, up to scaling factor
    // compute V = Q * 'U'
    // rescale E, if needed
    mat_mul_mp(V, Q, H);
    if (mpf_cmp_ui(norm, 1) > 0)
    { // rescale E
      for (i = 0; i < m; i++)
      { // undo the normalization by multiplying the entries of E by norm
        mpf_mul(E->coord[i].r, E->coord[i].r, norm);
      }
    }

    // clear
    clear_mat_mp(H); clear_mat_mp(Q); clear_mat_mp(R); clear_mat_mp(P);
  }

  // clear
  mpf_clear(norm); mpf_clear(norm_inv);

  return retVal;
}

int svd_R_jacobi_UE_mp(mat_mp U, vec_mp E, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the SVD of R - upper triangular matrix -  using  *
* Jacobi iterations                                             *
* ARGUMENTS: R = U*E*V^H, U - m x m, E - m x n, V would be n x n*
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the SVD that it has done          *
\***************************************************************/
{
  int i, j, retVal, m = R->rows, n = R->cols;
  mat_mp tempR;

  // do error checking
  if (m < n)
  {
    printf("ERROR: The matrix has more rows than columns in svd_R_jacobi_UE_d!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  // so we have m >= n

  if (n == 0) // if R has no cols - immediately exit
  {
    change_size_mat_mp(U, m, m);
    change_size_vec_mp(E, n);
    U->rows = U->cols = m;
    E->size = n;
    make_matrix_ID_mp(U, m, m);
    return 0;
  }

  // so m >= n >= 1

  // copy the main n x n part of R to tempR
  init_mat_mp(tempR, n, n);
  tempR->rows = tempR->cols = n;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      set_mp(&tempR->entry[i][j], &R->entry[i][j]);
    }

  // do jacobi iterations on tempR
  retVal = jacobi_UE_mp(U, E, tempR, its, rank_tol, tol_conv, tol_sign, largeChange);

  // adjust U & E based on m  - U = [[U, 0][0, I]]
  increase_size_mat_mp(U, m, m);
  U->rows = U->cols = m;
  for (i = 0; i < m; i++)
    for (j = 0; j < m; j++)
      if (i < n)
      {
        if (j >= n)
        { // U[i][j] = 0
          set_zero_mp(&U->entry[i][j]);
        }
      }
      else // i >= n
      {
        if (j < n)
        { // U[i][j] & E[i][j] = 0
          set_zero_mp(&U->entry[i][j]);
        }
        else // j >= n
        { // U[i][j] is (i == j)
          set_double_mp(&U->entry[i][j], i == j, 0);
        }
      }

  clear_mat_mp(tempR);

  return retVal;
}

int jacobi_UE_mp(mat_mp U, vec_mp E, mat_mp R, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on R to find its singular values*
* and left singular values                                      *
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: R has to be square                                     *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = R->rows, n = R->cols;
  mpf_t stoppingError, prevError1, prevError2;
  mpf_t da, db, dr, di, sizeC_sqr, sizeC;
  comp_mp c, t, tempComp;

  int *perm = NULL;
  mpf_t *sv = NULL;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix has needs to be square in jacobi_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(rank_tol, 0) < 0 || mpf_cmp_ui(tol_conv, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_mat_mp(U, n, n);
    change_size_vec_mp(E, n);
    U->rows = U->cols = E->size = n;
    return 0;
  }
  // so we can assume that m == n > 0

  // initialize MP
  mpf_init(stoppingError); mpf_init(prevError1); mpf_init(prevError2);
  mpf_init(da); mpf_init(db); mpf_init(dr); mpf_init(di); mpf_init(sizeC_sqr); mpf_init(sizeC);
  init_mp(c); init_mp(t); init_mp(tempComp);

  // copy R to U
  mat_cp_mp(U, R);

  // main loop
  count = 0;
  // initialize errors to something
  mpf_set_ui(stoppingError, 0);
  mpf_set_ui(prevError1, 0);
  mpf_set_ui(prevError2, 0);

  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    mpf_set(prevError2, prevError1);
    mpf_set(prevError1, stoppingError);
    // reset stoppingError to 0
    mpf_set_ui(stoppingError, 0);
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for(i = 0; i < n-1; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        mpf_set_ui(da, 0);
        mpf_set_ui(db, 0);
        set_zero_mp(c);
        for (k = 0; k < n; k++)
        {
          mpf_mul(dr, U->entry[k][i].r, U->entry[k][i].r);
          mpf_mul(di, U->entry[k][i].i, U->entry[k][i].i);
          mpf_add(dr, dr, di);
          mpf_add(da, da, dr);
          mpf_mul(dr, U->entry[k][j].r, U->entry[k][j].r);
          mpf_mul(di, U->entry[k][j].i, U->entry[k][j].i);
          mpf_add(dr, dr, di);
          mpf_add(db, db, dr);
          conjugate_mp(t, &U->entry[k][i]);
          sum_mul_mp(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        mpf_mul(dr, c->r, c->r);
        mpf_mul(di, c->i, c->i);
        mpf_add(sizeC_sqr, dr, di);
        mpf_sqrt(sizeC, sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (mpf_cmp(da, tol_sign) > 0 && mpf_cmp(db, tol_sign) > 0 && mpf_cmp(sizeC, tol_sign) > 0)
        { // we can safely find the stopping condition estimate - sizeC / sqrt(da * db)
          mpf_mul(dr, da, db);
          mpf_sqrt(dr, dr);
          mpf_div(dr, sizeC, dr);
          if (mpf_cmp(dr, stoppingError) > 0)
            mpf_set(stoppingError, dr);
        }

        if (mpf_cmp(sizeC, tol_sign) > 0)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          mpf_mul_ui(dr, sizeC_sqr, 2);
          mpf_sub(di, db, da);
          mpf_div(db, di, dr);
          mpf_mul(tempComp->r, c->r, db);
          mpf_mul(tempComp->i, c->i, db);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_mp2(t, tempComp, tol_sign);
          mpf_mul(dr, tempComp->r, tempComp->r);
          mpf_mul(di, tempComp->i, tempComp->i);
          mpf_add(db, dr, di);
          mpf_sqrt(da, db);
          mpf_add_ui(db, db, 1);
          mpf_sqrt(db, db);
          mpf_add(db, da, db);
          mpf_ui_div(db, 1, db);

          mpf_mul(t->r, t->r, db);
          mpf_mul(t->i, t->i, db);
        }
        else
        {
          set_zero_mp(t);
        }

        if (!(mpfr_zero_p(t->r) && mpfr_zero_p(t->i)))
        { // we need to do the rotation
          // find db = 1 / sqrt(1 + |t|^2)
          mpf_mul(dr, t->r, t->r);
          mpf_mul(di, t->i, t->i);
          mpf_add(db, dr, di);
          mpf_add_ui(db, db, 1);
          mpf_sqrt(db, db);
          mpf_ui_div(db, 1, db);
 
          // find t = db * t
          mpf_mul(t->r, t->r, db);
          mpf_mul(t->i, t->i, db);

          // tempComp = -conj(t)
          mpf_neg(tempComp->r, t->r);
          mpf_set(tempComp->i, t->i);

          // update columns i & j of U
          for (k = 0; k < n; k++)
          { // update U
            set_mp(c, &U->entry[k][i]);

            mpf_mul(U->entry[k][i].r, c->r, db);
            mpf_mul(U->entry[k][i].i, c->i, db);
            sum_mul_mp(&U->entry[k][i], tempComp, &U->entry[k][j]);

            mpf_mul(U->entry[k][j].r, U->entry[k][j].r, db);
            mpf_mul(U->entry[k][j].i, U->entry[k][j].i, db);
            sum_mul_mp(&U->entry[k][j], c, t);
          }
        }
      }

      // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
      mpf_sub(da, stoppingError, prevError1);
      mpf_abs(da, da);
      mpf_sub(db, prevError1, prevError2);
      mpf_abs(db, db);
      if (mpf_cmp(stoppingError, tol_conv) <= 0 || (count > 2 && mpf_cmp(da, tol_conv) <= 0 && mpf_cmp(db, tol_conv) <= 0))
        convergence = 1;
      else
        convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  // initialize perm & sv and see if they are in order
  perm = (int *)bmalloc(n * sizeof(int));
  sv = (mpf_t *)bmalloc(n * sizeof(mpf_t));
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = sv[j]
    mpf_set_ui(da, 0);
    for (i = 0; i < n; i++)
    {
      mpf_mul(dr, U->entry[i][j].r, U->entry[i][j].r);
      mpf_mul(di, U->entry[i][j].i, U->entry[i][j].i);
      mpf_add(dr, dr, di);
      mpf_add(da, da, dr);
    }
    mpf_sqrt(E->coord[j].r, da);
    mpf_set_ui(E->coord[j].i, 0);

    mpf_init_set(sv[j], E->coord[j].r);

    perm[j] = j;
    if (j > 0 && mpf_cmp(sv[j-1], sv[j]) < 0)
      order = 0;
  }

  if (!order)
  { // we need to find the order and update U, E & V

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (mpf_cmp(sv[j], sv[i]) > 0)
        { // swap j & i
          mpf_set(da, sv[j]);
          mpf_set(sv[j], sv[i]);
          mpf_set(sv[i], da);

          k = perm[j];
          perm[j] = perm[i];
          perm[i] = k;
        }

    // update U = U * P
    for (j = 0; j < n; j++)
    { // check to see if U & V need columns swapped
      if (perm[j] != j)
      {
        k = perm[j];
        // update E
        set_mp(tempComp, &E->coord[j]);
        set_mp(&E->coord[j], &E->coord[k]);
        set_mp(&E->coord[k], tempComp);

        for (i = 0; i < n; i++)
        { // update U
          set_mp(tempComp, &U->entry[i][j]);
          set_mp(&U->entry[i][j], &U->entry[i][k]);
          set_mp(&U->entry[i][k], tempComp);
        }
        // update perm
        for (i = 0; i < n; i++)
          if (perm[i] == j)
          {
            perm[i] = perm[j];
            perm[j] = j;
            i = n;
          }
      }
    }
  }

  // update U - columns need normalized
  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  mpf_set(da, E->coord[0].r); // initialize to something
  for (j = 0; j < n; j++)
  {
    // the 2-norms of the jth column of U is E->entry[j][j].r
    mpf_set(db, da); // the previous singular value
    mpf_set(da, E->coord[j].r); // the current singular value

    if (checkGood_mp(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column so we need to normalize this column of U
      count++; // this one is good!

      mpf_ui_div(dr, 1, da);
      for (i = 0; i < n; i++)
      { // update U
        mpf_mul(U->entry[i][j].r, U->entry[i][j].r, dr);
        mpf_mul(U->entry[i][j].i, U->entry[i][j].i, dr);
      }
    }
    else
    { // column is the zero column - expand to unitary matrix
      expand_unitary_mp(U, j, n);
      // exit loop
      j = n;
    }
  }

  // clear MP
  for (j = n - 1; j >= 0; j--)
    mpf_clear(sv[j]);
  mpf_clear(stoppingError); mpf_clear(prevError1); mpf_clear(prevError2);
  mpf_clear(da); mpf_clear(db); mpf_clear(dr); mpf_clear(di); mpf_clear(sizeC_sqr); mpf_clear(sizeC);
  clear_mp(c); clear_mp(t); clear_mp(tempComp);

  // free memory
  free(sv);
  free(perm);

  // return the correct value based on convergence - either corank or -1
  if (convergence)
    return (n - count);
  else
    return -1;
}

//////// FIND SINGULAR VALUES ONLY //////////

int svd_jacobi_E_mp(vec_mp E, mat_mp InputMat, int its, mpf_t rank_tol, mpf_t tol_conv_Jacobi, mpf_t tol_QR, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: finds the singular values of InputMat using Jacobi its *
* ARGUMENTS: InputMat = U*E*V^H                                 *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: Even if jacobi does not converge with its, we still    *
* return the approximation of the singular values done          *
\***************************************************************/
{
  int i, j, retVal, transposed, m = InputMat->rows, n = InputMat->cols;
  mpf_t norm, norm_inv;
  mat_mp H, R, P;

  // do error checking
  if (m == 0 || n == 0) // there are no rows/cols so we can exit immediately
  {
    change_size_vec_mp(E, n);
    E->size = n;
    return 0;
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(tol_conv_Jacobi, 0) < 0 || mpf_cmp_ui(tol_QR, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0 || mpf_cmp_ui(largeChange, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in svd_jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize
  mpf_init(norm); mpf_init(norm_inv); 

  // for best results, we should normalize InputMat based on its infinity norm if it is > 1
  infNormMat_mp2(norm, InputMat);

  // normalize and make 'skinny'
  if (m < n)
  { // setup H = (normalized) InputMat^T
    init_mat_mp(H, n, m); init_mat_mp(R, n, m); init_mat_mp(P, m, m);
    H->rows = n;
    H->cols = m;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // need to normalize and transpose
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        { // normalize and transpose
          conjugate_mp(&H->entry[j][i], &InputMat->entry[i][j]);
          mul_rmpf_mp(&H->entry[j][i], &H->entry[j][i], norm_inv);
        }
    }
    else
    { // just transpose
      transpose_mp(H, InputMat);
    }
    // setup m & n
    m = H->rows;
    n = H->cols;
    transposed = 1;
  }
  else
  { // setup H = (normalized) InputMat
    init_mat_mp(H, m, n); init_mat_mp(R, m, n); init_mat_mp(P, n, n);
    H->rows = m;
    H->cols = n;

    if (mpf_cmp_ui(norm, 1) > 0)
    { // need to normalize
      mpf_ui_div(norm_inv, 1, norm);
      for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
        { // normalize
          mul_rmpf_mp(&H->entry[i][j], &InputMat->entry[i][j], norm_inv);
        }
    }
    else
    { // just copy
      mat_cp_mp(H, InputMat);
    }
    transposed = 0;
  }
  // so we have m >= n >= 1

  // do the QR decomposition - 1 at the end so that QR presorts the rows ( H * P = Q * R, wher P is n x n, Q would be m x m & R is m x n)
  QR_R_mp(R, P, H, tol_QR, tol_sign, largeChange, 1);

  // since m >= n, we want to chop off the zeros on the bottom of R and then run the SVD on R^T
  R->rows = R->cols = n;
  transpose_mp(H, R);

  // find the singular values of H using jacobi iterations
  // NOTE: jacobi will always return its approximation of the svd even if it does not converge within the desired number of iterations
  retVal = jacobi_E_mp(E, H, its, rank_tol, tol_conv_Jacobi, tol_sign, largeChange);

  if (mpf_cmp_ui(norm, 1) > 0)
  { // undo the normalization by multiplying the entries of E by norm
    for (i = 0; i < n; i++)
      mpf_mul(E->coord[i].r, E->coord[i].r, norm);
  }

  // clear
  mpf_clear(norm); mpf_clear(norm_inv);
  clear_mat_mp(H); clear_mat_mp(R); clear_mat_mp(P);

  return retVal;
}

int jacobi_E_mp(vec_mp E, mat_mp U, int its, mpf_t rank_tol, mpf_t tol_conv, mpf_t tol_sign, mpf_t largeChange)
/***************************************************************\
* USAGE: does Jacobi iterations on U to find its singular values*
* ARGUMENTS:                                                    *
* RETURN VALUES: returns the corank if successful (>= 0)        *
*                returns -1 if not successful                   *
* NOTES: U has to be square                                     *
* Also, even if it does not converge with its, we still return  *
* the approximation of the SVD that it has done                 *
\***************************************************************/
{
  int i, j, k, convergence, count, order, m = U->rows, n = U->cols;
  mpf_t stoppingError, prevError1, prevError2;
  mpf_t da, db, dr, di, sizeC_sqr, sizeC;
  comp_mp c, t, tempComp;

  // do error checking
  if (m != n)
  {
    printf("ERROR: The matrix has needs to be square in jacobi_mp!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (its < 1)
  {
    printf("ERROR: The number of iterations needs to be positive in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (mpf_cmp_ui(rank_tol, 0) < 0 || mpf_cmp_ui(tol_conv, 0) < 0 || mpf_cmp_ui(tol_sign, 0) < 0)
  {
    printf("ERROR: The tolerance needs to be non-negative in jacobi_mp!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  // so we have m == n

  if (n == 0) // if InputMat has no rows/cols - immediately exit
  {
    change_size_vec_mp(E, n);
    E->size = n;
    return 0;
  }
  // so we can assume that m == n > 0

  // initialize MP
  mpf_init(stoppingError); mpf_init(prevError1); mpf_init(prevError2);
  mpf_init(da); mpf_init(db); mpf_init(dr); mpf_init(di); mpf_init(sizeC_sqr); mpf_init(sizeC);
  init_mp(c); init_mp(t); init_mp(tempComp);

  // main loop
  count = 0;
  // initialize errors to something
  mpf_set_ui(stoppingError, 0);
  mpf_set_ui(prevError1, 0);
  mpf_set_ui(prevError2, 0);

  do
  {
    convergence = 1; // initialize convergence to 1
    // update previous error
    mpf_set(prevError2, prevError1);
    mpf_set(prevError1, stoppingError);
    // reset stoppingError to 0
    mpf_set_ui(stoppingError, 0);
    count++; // increment the number of iterations that we have done

    // loop over all pairs (i,j) with 0 <= i < j < n
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
      { // find a, b, & c
        mpf_set_ui(da, 0);
        mpf_set_ui(db, 0);
        set_zero_mp(c);
        for (k = 0; k < n; k++)
        {
          mpf_mul(dr, U->entry[k][i].r, U->entry[k][i].r);
          mpf_mul(di, U->entry[k][i].i, U->entry[k][i].i);
          mpf_add(dr, dr, di);
          mpf_add(da, da, dr);

          mpf_mul(dr, U->entry[k][j].r, U->entry[k][j].r);
          mpf_mul(di, U->entry[k][j].i, U->entry[k][j].i);
          mpf_add(dr, dr, di);
          mpf_add(db, db, dr);

          conjugate_mp(t, &U->entry[k][i]);
          sum_mul_mp(c, t, &U->entry[k][j]);
        }

        // find sizeC_sqr & sizeC
        mpf_mul(dr, c->r, c->r);
        mpf_mul(di, c->i, c->i);
        mpf_add(sizeC_sqr, dr, di);
        mpf_sqrt(sizeC, sizeC_sqr);

        // update stoppingError
        // NOTE: da, db & sizeC are positive real numbers
        if (mpf_cmp(da, tol_sign) > 0 && mpf_cmp(db, tol_sign) > 0 && mpf_cmp(sizeC, tol_sign) > 0)
        { // we can safely find the stopping condition estimate - sizeC / sqrt(da * db)
          mpf_sqrt(dr, da);
          mpf_sqrt(di, db);
          mpf_mul(dr, dr, di);
          mpf_div(dr, sizeC, dr);
          if (mpf_cmp(dr, stoppingError) > 0)
            mpf_set(stoppingError, dr);
        }

        if (mpf_cmp(sizeC, tol_sign) > 0)
        { // since we have to invert conj(c), we need to make sure that it is large enough

          // find tempComp = (db - da) / (2 * conj(c)) = c * (db - da) / (2 * |c|^2)
          mpf_sub(dr, db, da);
          mul_rmpf_mp(tempComp, c, dr);

          mpf_mul_ui(dr, sizeC_sqr, 2);
          mpf_ui_div(dr, 1, dr);
          mul_rmpf_mp(tempComp, tempComp, dr);

          // find t = sign(tempComp) / (|tempComp| + sqrt(1 + |tempComp|^2))
          sign_mp2(t, tempComp, tol_sign);
          mpf_mul(dr, tempComp->r, tempComp->r);
          mpf_mul(di, tempComp->i, tempComp->i);
          mpf_add(db, dr, di);
          mpf_sqrt(da, db);

          mpf_add_ui(db, db, 1);
          mpf_sqrt(db, db);
          mpf_add(db, da, db);
          mpf_ui_div(db, 1, db);

          mul_rmpf_mp(t, t, db);
        }
        else
        {
          set_zero_mp(t);
        }

        if (!mpfr_zero_p(t->r) || !mpfr_zero_p(t->i))
        { // we need to do the rotation
          // find db = 1 / sqrt(1 + |t|^2)
          mpf_mul(dr, t->r, t->r);
          mpf_mul(di, t->i, t->i);
          mpf_add(db, dr, di);
          mpf_add_ui(db, db, 1);
          mpf_sqrt(db, db);
          mpf_ui_div(db, 1, db);

          // find t = db * t
          mul_rmpf_mp(t, t, db);

          // tempComp = -conj(t)
          mpf_neg(tempComp->r, t->r);
          mpf_set(tempComp->i, t->i);

          // update columns i & j of U
          for (k = 0; k < n; k++)
          { // update U
            set_mp(c, &U->entry[k][i]);

            mul_rmpf_mp(&U->entry[k][i], c, db);
            sum_mul_mp(&U->entry[k][i], tempComp, &U->entry[k][j]);

            mul_rmpf_mp(&U->entry[k][j], &U->entry[k][j], db);
            sum_mul_mp(&U->entry[k][j], c, t);
          }
        }
      }

      // check for the stopping criterion - quit if the error is small or the error has not changed too much in the previous 2 loops
      mpf_sub(da, stoppingError, prevError1);
      mpf_abs(da, da);
      mpf_sub(db, prevError1, prevError2);
      mpf_abs(db, db);
      if (mpf_cmp(stoppingError, tol_conv) <= 0 || (count > 2 && mpf_cmp(da, tol_conv) <= 0 && mpf_cmp(db, tol_conv) <= 0))
        convergence = 1;
      else
        convergence = 0;

  } while (convergence == 0 && count < its);

  // it could happen that the norms of the columns of U are not in decreasing order - i.e. the singular values - so we need to make sure of this
  change_size_vec_mp(E, n);
  order = 1;
  for (j = 0; j < n; j++)
  { // find the 2-norms of the jth column of U = E[j]
    mpf_set_ui(da, 0);
    for (i = 0; i < n; i++)
    {
      mpf_mul(dr, U->entry[i][j].r, U->entry[i][j].r);
      mpf_mul(di, U->entry[i][j].i, U->entry[i][j].i);
      mpf_add(dr, dr, di);
      mpf_add(da, da, dr);
    }
    mpf_sqrt(E->coord[j].r, da);
    mpf_set_ui(E->coord[j].i, 0);

    if (j > 0 && mpf_cmp(E->coord[j-1].r, E->coord[j].r) < 0)
      order = 0;
  }

  if (!order)
  { // we need to find the order

    // bubble sort to find the correct ordering
    for (i = 0; i < n; i++)
      for (j = i+1; j < n; j++)
        if (mpf_cmp(E->coord[j].r, E->coord[i].r) > 0)
        { // swap j & i
          mpf_set(da, E->coord[j].r);
          mpf_set(E->coord[j].r, E->coord[i].r);
          mpf_set(E->coord[i].r, da);
        }
  }

  // since we have an ordering on E, once we have a bad singular value, all the ones after it are bad as well!
  count = 0;
  mpf_set(da, E->coord[0].r); // initialize to something
  for (j = 0; j < n; j++)
  { // the 2-norms of the jth column of U is E->coord[j].r
    mpf_set(db, da); // the previous singular value
    mpf_set(da, E->coord[j].r); // the current singular value

    if (checkGood_mp(da, db, rank_tol, largeChange)) // check to make sure that the sing value is not too small and not too much difference than the last one
    { // column is not the zero column 
      count++; // this one is good!
    }
    else
    { // exit loop
      j = n;
    }
  }

  mpf_clear(stoppingError); mpf_clear(prevError1); mpf_clear(prevError2);
  mpf_clear(da); mpf_clear(db); mpf_clear(dr); mpf_clear(di); mpf_clear(sizeC_sqr); mpf_clear(sizeC);
  clear_mp(c); clear_mp(t); clear_mp(tempComp);

  // return the correct value based on convergence - either corank or -1
  if (convergence)
    return (n - count);
  else
    return -1;
}
 
//////// EXPAND TO A UNITARY MATRIX ///////////

void expand_unitary_mp(mat_mp U, int startCol, int endCol)
/***************************************************************\
* USAGE: puts in random orthonormal columns                     *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, m = U->rows, n = U->cols;
  mpf_t norm;
  comp_mp proj, temp;

  // error checking - make sure 0 <= startCol < endCol <= m
  if (!(0 <= startCol && startCol < endCol && endCol <= m))
  {
    printf("ERROR: The columns are not correct in expand_unitary_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (m < n) // need to make sure that we have a square or tall matrix (otherwise cannot make lin. indep.)
  {
    printf("ERROR: The matrix needs to have atleast as many rows as columns in expand_unitary_d!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // initialize MP
  mpf_init(norm);
  init_mp(proj); init_mp(temp);

  for (j = startCol; j < endCol; j++)
  { // create a random column
    for (i = 0; i < m; i++)
      get_comp_rand_mp(&U->entry[i][j]);

    // loop over the other columns to create an orthogonal vector
    for (k = 0; k < j; k++)
    { // set proj to 0
      set_zero_mp(proj);
      // find U[:][k]^T * U[:][j]
      for (i = 0; i < m; i++)
      { // proj += conj(U[i][k]) * U[i][j]
        conjugate_mp(temp, &U->entry[i][k]);
        sum_mul_mp(proj, temp, &U->entry[i][j]);
      }

      // negate proj
      neg_mp(proj, proj);

      // update U[:][j] += proj * U[:][k] & find its norm
      for (i = 0; i < m; i++)
      {
        sum_mul_mp(&U->entry[i][j], proj, &U->entry[i][k]);
      }
    }

    // find the norm
    mpf_set_ui(norm, 0);
    for (i = 0; i < m; i++)
    {
      mpf_mul(temp->r, U->entry[i][j].r, U->entry[i][j].r);
      mpf_mul(temp->i, U->entry[i][j].i, U->entry[i][j].i);
      mpf_add(norm, norm, temp->r);
      mpf_add(norm, norm, temp->i);
    }
    mpf_sqrt(norm, norm);

    // invert norm
    mpf_ui_div(norm, 1, norm);

    // normalize U[:][j]
    for (i = 0; i < m; i++)
    {
      mul_rmpf_mp(&U->entry[i][j], &U->entry[i][j], norm);
    }
  }

  return;
}

