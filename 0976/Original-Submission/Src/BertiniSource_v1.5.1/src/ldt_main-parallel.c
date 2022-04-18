// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "localdim.h"

// determine if the point is isolated and if it is, find its multiplicity

int is_isolated(int *mult, prog_deriv_t *deriv, point_d pt1_d, point_mp pt1_mp, int pt1_prec, point_d pt2_d, point_mp pt2_mp, int pt2_prec, int mult_bound, tracker_config_t *T, int printHilbert)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if it is isolated, and if is, its multiplicity *
* NOTES: assume that deriv is setup for atleast 1st order derivs*
\***************************************************************/
{
  int retVal, order, *hilbertFn = NULL, *rowSize = NULL, *colSize = NULL;

  // error checking
  if (deriv->order < 1)
  {
    printf("ERROR: The derivatives are not setup properly when doing the local dimension test!!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (mult_bound < 1)
  {
    printf("ERROR: The multiplicity bound must be a positive integer for the local dimension test!!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  if (T->MPType == 0)
  { // performing double precision local dim test
    retVal = is_isolated_d(mult, &hilbertFn, &order, deriv, pt1_d, pt2_d, mult_bound, T, &rowSize, &colSize, printHilbert);
  }
  else if (T->MPType == 1)
  { // performing fixed multi precision local dim test
    retVal = is_isolated_mp(mult, &hilbertFn, &order, deriv, pt1_mp, pt2_mp, mult_bound, T, &rowSize, &colSize, printHilbert);
  }
  else
  { // performing AMP local dim test
    retVal = is_isolated_amp(mult, &hilbertFn, &order, deriv, pt1_d, pt1_mp, pt1_prec, pt2_d, pt2_mp, pt2_prec, mult_bound, T, &rowSize, &colSize, printHilbert);
  }

  // clear memory
  free(hilbertFn);

  return retVal;
}

int is_isolated_d(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_d pt1_d, point_d pt2_d, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if it is isolated, and if is, its multiplicity *
* NOTES: assume that deriv is setup for atleast 1st order derivs*
\***************************************************************/
{
  int retVal, curr_order, corank, numHH = 0;
  int *ranks = NULL;
  double max_CN, max_SV_ratio, SV_tol, CN, smallest_nonzero, largest_zero;
  vec_d fn_d, fnDeriv_d, lin_d, linDeriv_d, *HH = NULL;
  mat_d MM1_d, MM2_d;

  // initialize
  init_vec_d(fn_d, 0); init_vec_d(fnDeriv_d, 0); init_vec_d(lin_d, 0); init_vec_d(linDeriv_d, 0);
  init_mat_d(MM1_d, 0, 0); init_mat_d(MM2_d, 0, 0);

  // setup the tolerances
  max_CN = 1e13;
  max_SV_ratio = T->ratioTol; 
  SV_tol = MAX(T->sing_val_zero_tol, 1e-15);

  // setup hilbertFn
  *hilbertOrder = 1;
  *hilbertFn = (int *)bmalloc(*hilbertOrder * sizeof(int));
  (*hilbertFn)[*hilbertOrder - 1] = 1;

  // initialize retVal & mult
  retVal = 0;
  *mult = 1;

  if (printHilbert)
    printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  // loop until we have stability of the size of the null space or we have exceeded the multiplicity bound
  corank = curr_order = 0;

  do
  { // increase the current order
    curr_order++;

    // increase the size of hilbertFn
    *hilbertOrder += 1;
    *hilbertFn = (int *)brealloc(*hilbertFn, *hilbertOrder * sizeof(int));

    // see if deriv is setup to enough this order
    while (curr_order > deriv->order)
    { // setup the next order
      setupNext_derivs(deriv);
    }

    // evaluate enough derivatives to setup the multiplicity matrix for this order for pt1
    evalDeriv_d(fn_d, fnDeriv_d, lin_d, linDeriv_d, pt1_d, deriv);
    // setup the multiplicity matrix of this order for pt1
    setup_multiplicity_matrix_d(MM1_d, curr_order, fn_d, fnDeriv_d, lin_d, linDeriv_d, deriv->numVars, deriv->numFuncs, deriv->numLinears);

    // evaluate enough derivatives to setup the multiplicity matrix for this order for pt2
    evalDeriv_d(fn_d, fnDeriv_d, lin_d, linDeriv_d, pt2_d, deriv);
    // setup the multiplicity matrix of this order for pt2
    setup_multiplicity_matrix_d(MM2_d, curr_order, fn_d, fnDeriv_d, lin_d, linDeriv_d, deriv->numVars, deriv->numFuncs, deriv->numLinears);

    // store the size of the old matrix
    *rowSize = (int *)brealloc(*rowSize, curr_order * sizeof(int));
    *colSize = (int *)brealloc(*colSize, curr_order * sizeof(int));
    (*rowSize)[curr_order - 1] = MM1_d->rows;
    (*colSize)[curr_order - 1] = MM1_d->cols;

    // find the corank
    corank = rank_MM_d(&CN, &smallest_nonzero, &largest_zero, MM1_d, MM2_d, curr_order, *rowSize, *colSize, ranks, &numHH, &HH, max_CN, max_SV_ratio, SV_tol);

    // store the ranks
    ranks = (int *)brealloc(ranks, curr_order * sizeof(int));
    ranks[curr_order - 1] = curr_order == 1 ? (*colSize)[curr_order - 1] - corank : (*colSize)[curr_order - 1] - (*colSize)[curr_order - 2] - corank;

    // find hilbertFn = number of new null space vectors
    (*hilbertFn)[*hilbertOrder - 1] = *hilbertOrder == 2 ? corank - (*hilbertFn)[*hilbertOrder - 2] : corank;
    // update mult
    *mult += (*hilbertFn)[*hilbertOrder - 1];

    if (printHilbert)
      printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  } while ((*hilbertFn)[*hilbertOrder - 1] > 0 && *mult <= mult_bound && curr_order < T->maxDepthLDT);

  // determine if we know an answer
  if (*mult > mult_bound)
    retVal = 0; // not isolated
  else if (curr_order < T->maxDepthLDT)
    retVal = 1; // is isolated
  else
    retVal = -1; // don't know!

  // clear memory
  for (curr_order = numHH - 1; curr_order >= 0; curr_order--)
    clear_vec_d(HH[curr_order]);
  free(HH);
  free(ranks);
  clear_vec_d(fn_d); clear_vec_d(fnDeriv_d); clear_vec_d(lin_d); clear_vec_d(linDeriv_d);
  clear_mat_d(MM1_d); clear_mat_d(MM2_d);

  return retVal;
}

int is_isolated_mp(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_mp pt1_mp, point_mp pt2_mp, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if it is isolated, and if is, its multiplicity *
* NOTES: assume that deriv is setup for atleast 1st order derivs*
\***************************************************************/
{
  int retVal, curr_order, corank, numHH = 0;
  int *ranks = NULL;
  double max_CN, max_SV_ratio, SV_tol, CN, smallest_nonzero, largest_zero;
  vec_mp fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, *HH = NULL;
  mat_mp MM1_mp, MM2_mp;

  // initialize
  init_vec_mp(fn_mp, 0); init_vec_mp(fnDeriv_mp, 0); init_vec_mp(lin_mp, 0); init_vec_mp(linDeriv_mp, 0);
  init_mat_mp(MM1_mp, 0, 0); init_mat_mp(MM2_mp, 0, 0);

  // setup the tolerances
  corank = prec_to_digits(T->Precision) - 4;
  max_CN = MIN(1e150, pow(10, corank));
  max_SV_ratio = T->ratioTol; 
  SV_tol = MAX(T->sing_val_zero_tol, pow(10, -corank - 2));

  // setup hilbertFn
  *hilbertOrder = 1;
  *hilbertFn = (int *)bmalloc(*hilbertOrder * sizeof(int));
  (*hilbertFn)[*hilbertOrder - 1] = 1;

  // initialize retVal & mult
  retVal = 0;
  *mult = 1;

  if (printHilbert)
    printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  // loop until we have stability of the size of the null space or we have exceeded the multiplicity bound
  corank = curr_order = 0;
  do
  { // increase the current order
    curr_order++;

    // increase the size of hilbertFn
    *hilbertOrder += 1;
    *hilbertFn = (int *)brealloc(*hilbertFn, *hilbertOrder * sizeof(int));

    // see if deriv is setup to enough this order
    while (curr_order > deriv->order)
    { // setup the next order
      setupNext_derivs(deriv);
    }

    // evaluate enough derivatives to setup the multiplicity matrix for this order for pt1
    evalDeriv_mp(fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, pt1_mp, deriv);
    // setup the multiplicity matrix of this order for pt1
    setup_multiplicity_matrix_mp(MM1_mp, curr_order, fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, deriv->numVars, deriv->numFuncs, deriv->numLinears);

    // evaluate enough derivatives to setup the multiplicity matrix for this order for pt2
    evalDeriv_mp(fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, pt2_mp, deriv);
    // setup the multiplicity matrix of this order for pt2
    setup_multiplicity_matrix_mp(MM2_mp, curr_order, fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, deriv->numVars, deriv->numFuncs, deriv->numLinears);

    // store the size of the old matrix
    *rowSize = (int *)brealloc(*rowSize, curr_order * sizeof(int));
    *colSize = (int *)brealloc(*colSize, curr_order * sizeof(int));
    (*rowSize)[curr_order - 1] = MM1_mp->rows;
    (*colSize)[curr_order - 1] = MM1_mp->cols;

    // find the corank
    corank = rank_MM_mp(&CN, &smallest_nonzero, &largest_zero, MM1_mp, MM2_mp, curr_order, *rowSize, *colSize, ranks, &numHH, &HH, max_CN, max_SV_ratio, SV_tol);

    // store the ranks
    ranks = (int *)brealloc(ranks, curr_order * sizeof(int));
    ranks[curr_order - 1] = curr_order == 1 ? (*colSize)[curr_order - 1] - corank : (*colSize)[curr_order - 1] - (*colSize)[curr_order - 2] - corank;

    // find hilbertFn = number of new null space vectors
    (*hilbertFn)[*hilbertOrder - 1] = *hilbertOrder == 2 ? corank - (*hilbertFn)[*hilbertOrder - 2] : corank;

    // update mult
    *mult += (*hilbertFn)[*hilbertOrder - 1];

    if (printHilbert)
      printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  } while ((*hilbertFn)[*hilbertOrder - 1] > 0 && *mult <= mult_bound && curr_order < T->maxDepthLDT);

  // determine if we know an answer
  if (*mult > mult_bound)
    retVal = 0; // not isolated
  else if (curr_order < T->maxDepthLDT)
    retVal = 1; // is isolated
  else
    retVal = -1; // don't know!

  // clear memory
  for (curr_order = numHH - 1; curr_order >= 0; curr_order--)
    clear_vec_mp(HH[curr_order]);
  free(HH);
  free(ranks);
  clear_vec_mp(fn_mp); clear_vec_mp(fnDeriv_mp); clear_vec_mp(lin_mp); clear_vec_mp(linDeriv_mp);
  clear_mat_mp(MM1_mp); clear_mat_mp(MM2_mp);

  return retVal;
}

int is_isolated_amp(int *mult, int **hilbertFn, int *hilbertOrder, prog_deriv_t *deriv, point_d pt1_d, point_mp pt1_mp, int pt1_prec, point_d pt2_d, point_mp pt2_mp, int pt2_prec, int mult_bound, tracker_config_t *T, int **rowSize, int **colSize, int printHilbert)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: if it is isolated, and if is, its multiplicity *
* NOTES: assume that deriv is setup for atleast 1st order derivs*
\***************************************************************/
{
  int retVal, curr_order, corank, numHH = 0;
  int *ranks = NULL;
  double CN, smallest_nonzero, largest_zero, max_CN = 1e300, max_SV_ratio = T->ratioTol; 
  vec_d fn_d, fnDeriv_d, lin_d, linDeriv_d, *HH_d = NULL;
  vec_mp fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, *HH_mp = NULL;
  mat_d MM1_d, MM2_d;
  mat_mp MM1_mp, MM2_mp;

  // initialize
  init_vec_mp(fn_mp, 0); init_vec_mp(fnDeriv_mp, 0); init_vec_mp(lin_mp, 0); init_vec_mp(linDeriv_mp, 0);
  init_mat_mp(MM1_mp, 0, 0); init_mat_mp(MM2_mp, 0, 0);

  // initialize
  init_vec_d(fn_d, 0); init_vec_d(fnDeriv_d, 0); init_vec_d(lin_d, 0); init_vec_d(linDeriv_d, 0);
  init_vec_mp(fn_mp, 0); init_vec_mp(fnDeriv_mp, 0); init_vec_mp(lin_mp, 0); init_vec_mp(linDeriv_mp, 0);
  init_mat_d(MM1_d, 0, 0); init_mat_d(MM2_d, 0, 0);
  init_mat_mp(MM1_mp, 0, 0); init_mat_mp(MM2_mp, 0, 0);

  // setup hilbertFn
  *hilbertOrder = 1;
  *hilbertFn = (int *)bmalloc(*hilbertOrder * sizeof(int));
  (*hilbertFn)[*hilbertOrder - 1] = 1;

  // initialize retVal & mult
  retVal = 0;
  *mult = 1;

  if (printHilbert)
    printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  // loop until we have stability of the size of the null space or we have exceeded the multiplicity bound
  curr_order = 0;
  do
  { // increase the current order
    curr_order++;

    // increase the size of hilbertFn
    *hilbertOrder += 1;
    *hilbertFn = (int *)brealloc(*hilbertFn, *hilbertOrder * sizeof(int));

    // see if deriv is setup to enough this order
    while (curr_order > deriv->order)
    { // setup the next order
      setupNext_derivs(deriv);
    }

    if (pt1_prec < 64)
    { // setup MM1_d
      // evaluate enough derivatives to setup the multiplicity matrix for this order for pt1
      evalDeriv_d(fn_d, fnDeriv_d, lin_d, linDeriv_d, pt1_d, deriv);
      // setup the multiplicity matrix of this order for pt1
      setup_multiplicity_matrix_d(MM1_d, curr_order, fn_d, fnDeriv_d, lin_d, linDeriv_d, deriv->numVars, deriv->numFuncs, deriv->numLinears);
    }
    else
    { // setup MM1_mp using the correct precision
      initMP(pt1_prec);
      change_prec_vec_mp(fn_mp, pt1_prec); change_prec_vec_mp(fnDeriv_mp, pt1_prec);
      change_prec_vec_mp(lin_mp, pt1_prec); change_prec_vec_mp(linDeriv_mp, pt1_prec);
      change_prec_mat_mp(MM1_mp, pt1_prec);
      deriv->precision = pt1_prec;

      // evaluate enough derivatives to setup the multiplicity matrix for this order for pt1
      evalDeriv_mp(fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, pt1_mp, deriv);
      // setup the multiplicity matrix of this order for pt1
      setup_multiplicity_matrix_mp(MM1_mp, curr_order, fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, deriv->numVars, deriv->numFuncs, deriv->numLinears);
    }

    if (pt2_prec < 64)
    { // setup MM2_d
      // evaluate enough derivatives to setup the multiplicity matrix for this order for pt2
      evalDeriv_d(fn_d, fnDeriv_d, lin_d, linDeriv_d, pt2_d, deriv);
      // setup the multiplicity matrix of this order for pt2
      setup_multiplicity_matrix_d(MM2_d, curr_order, fn_d, fnDeriv_d, lin_d, linDeriv_d, deriv->numVars, deriv->numFuncs, deriv->numLinears);
    }
    else
    { // setup MM2_mp using the correct precision
      initMP(pt2_prec);
      change_prec_vec_mp(fn_mp, pt2_prec); change_prec_vec_mp(fnDeriv_mp, pt2_prec);
      change_prec_vec_mp(lin_mp, pt2_prec); change_prec_vec_mp(linDeriv_mp, pt2_prec);
      change_prec_mat_mp(MM2_mp, pt2_prec);
      deriv->precision = pt2_prec;

      // evaluate enough derivatives to setup the multiplicity matrix for this order for pt2
      evalDeriv_mp(fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, pt2_mp, deriv);
      // setup the multiplicity matrix of this order for pt2
      setup_multiplicity_matrix_mp(MM2_mp, curr_order, fn_mp, fnDeriv_mp, lin_mp, linDeriv_mp, deriv->numVars, deriv->numFuncs, deriv->numLinears);
    }

    // store the size of the old matrix
    *rowSize = (int *)brealloc(*rowSize, curr_order * sizeof(int));
    *colSize = (int *)brealloc(*colSize, curr_order * sizeof(int));
    (*rowSize)[curr_order - 1] = pt1_prec < 64 ? MM1_d->rows : MM1_mp->rows;
    (*colSize)[curr_order - 1] = pt1_prec < 64 ? MM1_d->cols : MM1_mp->cols;

    // compute the corank
    corank = rank_MM_amp(&CN, &smallest_nonzero, &largest_zero, MM1_d, MM1_mp, pt1_prec, MM2_d, MM2_mp, pt2_prec, curr_order, *rowSize, *colSize, ranks, &numHH, &HH_d, &HH_mp, max_CN, max_SV_ratio);

    // store the ranks
    ranks = (int *)brealloc(ranks, curr_order * sizeof(int));
    ranks[curr_order - 1] = curr_order == 1 ? (*colSize)[curr_order - 1] - corank : (*colSize)[curr_order - 1] - (*colSize)[curr_order - 2] - corank;

    // find hilbertFn = number of new null space vectors
    (*hilbertFn)[*hilbertOrder - 1] = *hilbertOrder == 2 ? corank - (*hilbertFn)[*hilbertOrder - 2] : corank;

    // update mult
    *mult += (*hilbertFn)[*hilbertOrder - 1];

    if (printHilbert)
      printf("h(%d) = %d,  mult(%d) = %d\n", *hilbertOrder - 1, (*hilbertFn)[*hilbertOrder - 1], *hilbertOrder - 1, *mult);

  } while ((*hilbertFn)[*hilbertOrder - 1] > 0 && *mult <= mult_bound && curr_order < T->maxDepthLDT);

  // determine if we know an answer
  if (*mult > mult_bound)
    retVal = 0; // not isolated
  else if ((*hilbertFn)[*hilbertOrder -1] <= 0)
    retVal = 1; // is isolated since we have stabilized
  else
    retVal = -1; // don't know!

  // clear memory
  for (curr_order = numHH - 1; curr_order >= 0; curr_order--)
  { 
    if (HH_d != NULL)
      clear_vec_d(HH_d[curr_order]);
    if (HH_mp != NULL)
      clear_vec_mp(HH_mp[curr_order]);
  }
  if (HH_d != NULL)
    free(HH_d);
  if (HH_mp != NULL)
    free(HH_mp);
  free(ranks);
  clear_vec_d(fn_d); clear_vec_d(fnDeriv_d); clear_vec_d(lin_d); clear_vec_d(linDeriv_d);
  clear_vec_mp(fn_mp); clear_vec_mp(fnDeriv_mp); clear_vec_mp(lin_mp); clear_vec_mp(linDeriv_mp);
  clear_mat_d(MM1_d); clear_mat_d(MM2_d);
  clear_mat_mp(MM1_mp); clear_mat_mp(MM2_mp);

  return retVal;
}

