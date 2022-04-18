// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "localdim.h"

// setup the multiplicy matrix for a given order

///////////////////// DOUBLE PRECISION //////////////////////////////

void setup_multiplicity_matrix_d(mat_d MM, int order, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int numFuncs, int numLinears)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the mult matrix of 'order' from func & derivVals *
\***************************************************************/
{ // for the 'order' multiplicity matrix, the monomials are of maximal degree 'order - 1' while the derivatives are of 'order'
  int i, num_mons, num_diffs, diff_num, mon_num, mon_degree, func_num, diff_order, row_num, col_num, max_mon_degree = order - 1, max_diff_order = order;
  int totalFuncs = numFuncs + numLinears, num_rows = 0, num_cols = 0;
  int *mon_array = (int *)bmalloc(numVars * sizeof(int)), *diff_array = (int *)bmalloc(numVars * sizeof(int));
  int *mon_size = (int *)bmalloc((max_mon_degree + 1) * sizeof(int)), *diff_size = (int *)bmalloc((max_diff_order + 1) * sizeof(int));
  vec_d vals;
  init_vec_d(vals, 0);

  // error checking
  if (order <= 0)
  {
    printf("ERROR: The order of the multiplicity matrix must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (funcVals->size != numFuncs)
  {
    printf("ERROR: The number of the functions is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (linVals->size != numLinears)
  {
    printf("ERROR: The number of the linears is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  //First, we set up the size of the multiplicity matrix 
  for (i = 0; i <= max_mon_degree; i++)
  { // setup mon_size and update num_rows
    mon_size[i] = combination(numVars + i - 1, i);
    num_rows += totalFuncs * mon_size[i];
  }
  for (i = 0; i <= max_diff_order; i++)
  { // setup diff_size and update num_cols
    diff_size[i] = combination(numVars + i - 1, i);
    num_cols += diff_size[i];
  }

  // setup MM to the correct size
  change_size_mat_d(MM, num_rows, num_cols);
  MM->rows = num_rows;
  MM->cols = num_cols;

  // setup vals from funcVals, linVals, derivVals & linDerivVals
  combine_func_deriv_vals_d(vals, funcVals, derivVals, linVals, linDerivVals, numVars, max_diff_order, diff_size);

  // main loop
  row_num = 0;
  for (mon_degree = 0; mon_degree <= max_mon_degree; mon_degree++)  //for each monomial degree
  { //determine the number of monomials of this degree
    num_mons = mon_size[mon_degree];

    for (func_num = 0; func_num < totalFuncs; func_num++) //for each function
    { //setup mon_array as the first monomial of this degree:
      start_array(mon_array, numVars, mon_degree);

      for (mon_num = 0; mon_num < num_mons; mon_num++)  //for each monomial of the given degree (the last row index)
      { // reset col_num back to 0 and loop over the columns
        col_num = 0;
        for (diff_order = 0; diff_order <= max_diff_order; diff_order++)  //for each order of the partials
        {
          num_diffs = diff_size[diff_order]; //determine the number of partials of this order

          //setup diff_array as the first diff of this order
          start_array(diff_array, numVars, diff_order);

          for (diff_num = 0; diff_num < num_diffs; diff_num++)  //for each diff of this order
          { // setup the [row_num][col_num] entry of MM
            setup_MM_entry_d(MM, row_num, col_num, vals, mon_array, diff_array, func_num, totalFuncs, numVars, diff_size);

            col_num++;  //Advance to the next col.
            advance_array(diff_array, numVars);  //move diff_array to the next partial.
          }
        }
        row_num++;  //Advance to the next row.
        advance_array(mon_array, numVars); //move mon_array to the next monomial.
      }
    }
  }

  // clear memory
  clear_vec_d(vals);
  free(mon_array); free(diff_array);
  free(mon_size);  free(diff_size);

  return;
}

void combine_func_deriv_vals_d(vec_d vals, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int max_diff_order, int *diff_size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: combines funcVals & derivVals into vals with scaling   *
\***************************************************************/
{
  int i, j, derivLoc, linDerivLoc, order, currLoc, size, numFuncs = funcVals->size, numLinears = linVals->size;
  int totalFuncs = numFuncs + numLinears;

  // setup size
  size = 0;
  for (i = 0; i <= max_diff_order; i++)
    size += diff_size[i];
  size *= totalFuncs;

  // set vals to the correct size
  change_size_vec_d(vals, size);
  vals->size = size;

  // setup top of vals - [funcVals, linVals]
  for (currLoc = 0; currLoc < totalFuncs; currLoc++)
    if (currLoc < numFuncs)
    {
      set_d(&vals->coord[currLoc], &funcVals->coord[currLoc]);
    }
    else
    {
      set_d(&vals->coord[currLoc], &linVals->coord[currLoc - numFuncs]);
    }

  // setup the derivatives
  derivLoc = 0;
  for (order = 1; order <= max_diff_order; order++)
  { // setup the top of this order with function derivatives
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < diff_size[order]; j++)
      {
        set_d(&vals->coord[currLoc], &derivVals->coord[derivLoc]);
        currLoc++;
        derivLoc++;
      }

    // setup the bottom of this order with linear derivatives
    if (order == 1)
    { // copy linear derivatives
      linDerivLoc = 0;
      for (i = 0; i < numLinears; i++)
        for (j = 0; j < diff_size[order]; j++)
        {
          set_d(&vals->coord[currLoc], &linDerivVals->coord[linDerivLoc]);
          currLoc++;
          linDerivLoc++;
        }
    }
    else
    { // all are zero
      for (i = 0; i < numLinears; i++)
        for (j = 0; j < diff_size[order];j++)
        {
          set_zero_d(&vals->coord[currLoc]);
          currLoc++;
        }
    }
  }

  if (max_diff_order > 1)
  { // we need to scale them accordingly (1 / alpha! * d_alpha(f))
    int *array = (int *)bmalloc(numVars * sizeof(int));
    double factor;

    currLoc = totalFuncs + totalFuncs * numVars; // functions/linears & first order derivatives do not need scaled
    for (order = 2; order <= max_diff_order; order++)
      for (i = 0; i < totalFuncs; i++)
      { // setup array to be [order,0,..,0]
        start_array(array, numVars, order);

        // loop over the partial derivatives of this order
        for (j = 0; j < diff_size[order]; j++)
        { // find array!
          factor = factorial_array(array, numVars);

          if (factor != 1)
          { // this value needs scaled by 1 / array!
            factor = 1 / factor;
            mul_rdouble_d(&vals->coord[currLoc], &vals->coord[currLoc], factor);
          }

          // advance the array
          advance_array(array, numVars);
          // increment currLoc
          currLoc++;
        }
      }

    // clear memory
    free(array);
  }

  return;
}

void setup_MM_entry_d(mat_d MM, int row_num, int col_num, vec_d vals, int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup MM[row_num][col_num]                             *
\***************************************************************/
{
  //See Jon's sheet on the multivariate Leibnitz rule for details.
  //Essentially, depending on the two arrays (mon_array and diff_array), it is very easy to evaluate the derivative.
  //If any one entry of diff_array is greater than the corresp. entry of mon_array, then the entry is 0.
  //Otherwise, it is just the mon_array-diff_array partial derivative of the function at hand.

  if (array_compare(mon_array, diff_array, numVars) == -1)  //if mon_array is not smaller, the entry is just 0.
  {
    set_zero_d(&MM->entry[row_num][col_num]); //Could be dropped if MM is initialized to all 0's.  (efficiency!!!)
  }
  else //otherwise (if mon_array is smaller), it is just a particular partial.
  { // find the partial index and copy this to the correct location
    int coord = find_partial_index(mon_array, diff_array, func_num, numFuncs, numVars, diff_size);
    set_d(&MM->entry[row_num][col_num], &vals->coord[coord]);
  }

  return;
}

// use monomial support to setup a multiplicity matrix with only the functions and their derivs (no terms of the form x^a*f)

void setup_multiplicity_matrix_mon_d(mat_d MM, int order, int **monomial_support, int *diff_size, vec_d funcVals, vec_d derivVals, vec_d linVals, vec_d linDerivVals, int numVars, int numFuncs, int numLinears)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the mult matrix of 'order' from func & derivVals *
\***************************************************************/
{
  int i, j, curr_order, startCol = 0, currCol = 0, derivLoc = 0, linDerivLoc = 0, totalFuncs = numFuncs + numLinears, num_rows = 0, num_cols = 0;
  int *non_support = (int *)bmalloc((order + 1) * sizeof(int)), *array = (int *)bmalloc(numVars * sizeof(int));
  double multiplier;

  // error checking
  if (order <= 0)
  {
    printf("ERROR: The order of the multiplicity matrix must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (funcVals->size != numFuncs)
  {
    printf("ERROR: The number of the functions is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (linVals->size != numLinears)
  {
    printf("ERROR: The number of the linears is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // First, we set up the size of the multiplicity matrix 
  num_rows = totalFuncs;
  num_cols = 0;
  for (i = 0; i <= order; i++)
  { // count the number not in the support
    non_support[i] = diff_size[i];
    for (j = 0; j < diff_size[i]; j++)
      non_support[i] -= monomial_support[i][j];
    num_cols += diff_size[i] - non_support[i];
  }

  // setup MM to the correct size
  change_size_mat_d(MM, num_rows, num_cols);
  MM->rows = num_rows;
  MM->cols = num_cols;

  // setup the evaluations
  startCol = currCol = 0;
  if (monomial_support[0][0])
  { // set the values
    for (i = 0; i < totalFuncs; i++)
      if (i < numFuncs)
      { // set function evaluation
        set_d(&MM->entry[i][currCol], &funcVals->coord[i]);
      }
      else
      { // set linear evaluation
        set_d(&MM->entry[i][currCol], &linVals->coord[i - numFuncs]);
      }
  }

  // setup the deriv entries of MM: row - function i, column - derivative j
  startCol = currCol = diff_size[0] - non_support[0];
  derivLoc = linDerivLoc = 0;
  for (curr_order = 1; curr_order <= order; curr_order++)
  { // setup the function derivative entries
    for (i = 0; i < numFuncs; i++)
    { // reset the counter
      currCol = startCol;
      for (j = 0; j < diff_size[curr_order]; j++)
      { // set value if monomial is used
        if (monomial_support[curr_order][j]) 
        { // set the value 
          set_d(&MM->entry[i][currCol], &derivVals->coord[derivLoc]);
          // increment counter
          currCol++;
        }
        derivLoc++;  
      }             
    }

    // setup linears
    if (curr_order == 1)
    { // linear derivs
      linDerivLoc = 0;
      for (i = 0; i < numLinears; i++)
      { // reset the counter
        currCol = startCol;
        for (j = 0; j < diff_size[curr_order]; j++)
        { // set the value if monomial is used
          if (monomial_support[curr_order][j]) 
          { // set the value
            set_d(&MM->entry[i + numFuncs][currCol], &linDerivVals->coord[linDerivLoc]);
            // increment counters
            currCol++;
          }
          linDerivLoc++;
        }
      }
    }
    else
    { // all zeros
      for (i = 0; i < numLinears; i++)
      { // reset the counter
        currCol = startCol;
        for (j = diff_size[curr_order] - non_support[curr_order]; j > 0; j--)
        { // set the value
          set_zero_d(&MM->entry[i + numFuncs][currCol]);
          // increment counter
          currCol++;
        }
      }
    }

    // update the startCol
    startCol += diff_size[curr_order] - non_support[curr_order];
  }

  // normalize the columns for the functions

  startCol = currCol = diff_size[0] - non_support[0] + diff_size[1] - non_support[1]; // start at order 2
  for (curr_order = 2; curr_order <= order; curr_order++)
  { // determine which monomials where used
    start_array(array, numVars, curr_order); // monomials of this degree
    currCol = startCol;
    for (j = 0; j < diff_size[curr_order]; j++)
    { // loop over the monomials for this degree
      if (monomial_support[curr_order][j])
      { // compute multiplier and normalize this column
        multiplier = factorial_array(array, numVars);
        if (multiplier != 1)
        { // these values need scaled by 1 / array!
          multiplier = 1 / multiplier;
          for (i = 0; i < numFuncs; i++)
            mul_rdouble_d(&MM->entry[i][currCol], &MM->entry[i][currCol], multiplier);
        }

        // increment column
        currCol++;
      }
      // advance the array
      advance_array(array, numVars);
    }

    // update the startCol
    startCol += diff_size[curr_order] - non_support[curr_order];
  }
  
  free(non_support); 
  free(array);

  return;
}

///////////////////// MULTI PRECISION //////////////////////////////

void setup_multiplicity_matrix_mp(mat_mp MM, int order, vec_mp funcVals, vec_mp derivVals, vec_mp linVals, vec_mp linDerivVals, int numVars, int numFuncs, int numLinears)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the mult matrix of 'order' from func & derivVals *
\***************************************************************/
{ // for the 'order' multiplicity matrix, the monomials are of maximal degree 'order - 1' while the derivatives are of 'order'
  int i, num_mons, num_diffs, diff_num, mon_num, mon_degree, func_num, diff_order, row_num, col_num, max_mon_degree = order - 1, max_diff_order = order;
  int totalFuncs = numFuncs + numLinears, num_rows = 0, num_cols = 0;
  int *mon_array = (int *)bmalloc(numVars * sizeof(int)), *diff_array = (int *)bmalloc(numVars * sizeof(int));
  int *mon_size = (int *)bmalloc((max_mon_degree + 1) * sizeof(int)), *diff_size = (int *)bmalloc((max_diff_order + 1) * sizeof(int));
  vec_mp vals;
  init_vec_mp(vals, 0);

  // error checking
  if (order <= 0)
  {
    printf("ERROR: The order of the multiplicity matrix must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (funcVals->size != numFuncs)
  {
    printf("ERROR: The number of the functions is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (linVals->size != numLinears)
  {
    printf("ERROR: The number of the linears is not correct when setting up the multiplicity matrix!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  //First, we set up the size of the multiplicity matrix
  for (i = 0; i <= max_mon_degree; i++)
  { // setup mon_size and update num_rows
    mon_size[i] = combination(numVars + i - 1, i);
    num_rows += totalFuncs * mon_size[i];
  }
  for (i = 0; i <= max_diff_order; i++)
  { // setup diff_size and update num_cols
    diff_size[i] = combination(numVars + i - 1, i);
    num_cols += diff_size[i];
  }

  // setup MM to the correct size
  change_size_mat_mp(MM, num_rows, num_cols);
  MM->rows = num_rows;
  MM->cols = num_cols;

  // setup vals from funcVals & derivVals
  combine_func_deriv_vals_mp(vals, funcVals, derivVals, linVals, linDerivVals, numVars, max_diff_order, diff_size);

  // main loop
  row_num = 0;
  for (mon_degree = 0; mon_degree <= max_mon_degree; mon_degree++)  //for each monomial degree
  { //determine the number of monomials of this degree
    num_mons = mon_size[mon_degree];

    for (func_num = 0; func_num < totalFuncs; func_num++)  //for each function
    { //setup mon_array as the first monomial of this degree:
      start_array(mon_array, numVars, mon_degree);

      for (mon_num = 0; mon_num < num_mons; mon_num++)  //for each monomial of the given degree (the last row index)
      { // reset col_num back to 0 and loop over the columns
        col_num = 0;
        for (diff_order = 0; diff_order <= max_diff_order; diff_order++)  //for each order of the partials
        { //determine the number of partials of this order
          num_diffs = diff_size[diff_order];

          //setup diff_array as the first diff of this order
          start_array(diff_array, numVars, diff_order);

          for (diff_num = 0; diff_num < num_diffs; diff_num++)  //for each diff of this order
          { // setup the MM[row_num][col_num]
            setup_MM_entry_mp(MM, row_num, col_num, vals, mon_array, diff_array, func_num, totalFuncs, numVars, diff_size);

            col_num++;  //Advance to the next col.
            advance_array(diff_array, numVars);  //move diff_array to the next partial.
          }
        }
        row_num++;  //Advance to the next row.
        advance_array(mon_array, numVars); //move mon_array to the next monomial.
      }
    }
  }

  // clear memory
  clear_vec_mp(vals);
  free(mon_array); free(diff_array);
  free(mon_size);  free(diff_size);

  return;
}

void combine_func_deriv_vals_mp(vec_mp vals, vec_mp funcVals, vec_mp derivVals, vec_mp linVals, vec_mp linDerivVals, int numVars, int max_diff_order, int *diff_size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: combines funcVals & derivVals into vals with scaling   *
\***************************************************************/
{
  int i, j, derivLoc, linDerivLoc, order, currLoc, size, numFuncs = funcVals->size, numLinears = linVals->size;
  int totalFuncs = numFuncs + numLinears;

  // setup size
  size = 0;
  for (i = 0; i <= max_diff_order; i++)
    size += diff_size[i];
  size *= totalFuncs;

  // make sure vals is of the correct size
  change_size_vec_mp(vals, size);
  vals->size = size;

  // setup top of vals - [funcVals, linVals]
  for (currLoc = 0; currLoc < totalFuncs; currLoc++)
    if (currLoc < numFuncs)
    {
      set_mp(&vals->coord[currLoc], &funcVals->coord[currLoc]);
    }
    else
    {
      set_mp(&vals->coord[currLoc], &linVals->coord[currLoc - numFuncs]);
    }

  // setup the derivatives
  derivLoc = 0;
  for (order = 1; order <= max_diff_order; order++)
  { // setup the top of this order with function derivatives
    for (i = 0; i < numFuncs; i++)
      for (j = 0; j < diff_size[order]; j++)
      {
        set_mp(&vals->coord[currLoc], &derivVals->coord[derivLoc]);
        currLoc++;
        derivLoc++;
      }

    // setup the bottom of this order with linear derivatives
    if (order == 1)
    { // copy linear derivatives
      linDerivLoc = 0;
      for (i = 0; i < numLinears; i++)
        for (j = 0; j < diff_size[order]; j++)
        {
          set_mp(&vals->coord[currLoc], &linDerivVals->coord[linDerivLoc]);
          currLoc++;
          linDerivLoc++;
        }
    }
    else
    { // all are zero
      for (i = 0; i < numLinears; i++)
        for (j = 0; j < diff_size[order];j++)
        {
          set_zero_mp(&vals->coord[currLoc]);
          currLoc++;
        }
    }
  }

  if (max_diff_order > 1)
  { // we need to scale them accordingly (1 / alpha! * d_alpha(f))
    int *array = (int *)bmalloc(numVars * sizeof(int));
    mpf_t factor;
    mpf_init(factor);

    currLoc = totalFuncs + totalFuncs * numVars; // functions/linears & first order derivatives do not need scaled
    for (order = 2; order <= max_diff_order; order++)
      for (i = 0; i < totalFuncs; i++)
      { // setup array to be [order,0,..,0]
        start_array(array, numVars, order);

        // loop over the partial derivatives of this order
        for (j = 0; j < diff_size[order]; j++)
        { // find array!
          factorial_array2(factor, array, numVars);

          if (mpf_cmp_ui(factor, 1) > 0)
          { // this value needs scaled by 1 / array!
            mpf_ui_div(factor, 1, factor);
            mul_rmpf_mp(&vals->coord[currLoc], &vals->coord[currLoc], factor);
          }

          // advance the array
          advance_array(array, numVars);
          // increment currLoc
          currLoc++;
        }
      }

    // clear memory
    mpf_clear(factor);
    free(array);
  }

  return;
}

void setup_MM_entry_mp(mat_mp MM, int row_num, int col_num, vec_mp vals, int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup MM[row_num][col_num]                             *
\***************************************************************/
{
  //See Jon's sheet on the multivariate Leibnitz rule for details.
  //Essentially, depending on the two arrays (mon_array and diff_array), it is very easy to evaluate the derivative.
  //If any one entry of diff_array is greater than the corresp. entry of mon_array, then the entry is 0.
  //Otherwise, it is just the mon_array-diff_array partial derivative of the function at hand.

  if (array_compare(mon_array, diff_array, numVars) == -1)  //if mon_array is not smaller, the entry is just 0.
  {
    set_zero_mp(&MM->entry[row_num][col_num]); //Could be dropped if MM is initialized to all 0's.  (efficiency!!!)
  }
  else //otherwise (if mon_array is smaller), it is just a particular partial.
  { // find the partial index and copy this to the correct location
    int coord = find_partial_index(mon_array, diff_array, func_num, numFuncs, numVars, diff_size);
    set_mp(&MM->entry[row_num][col_num], &vals->coord[coord]);
  }

  return;
}

///////////////////// GENERAL FUNCTIONS //////////////////////////////

void start_array(int *array, int length, int weight)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup array to be the starting one in the ordering     *
\***************************************************************/
{
  int i;

  // error checking
  if (length <= 0)
  {
    printf("ERROR: The length of the array must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (weight < 0)
  {
    printf("ERROR: The weight for the array must be non-negative!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  array[0] = weight;
  for (i = 1; i < length; i++)
    array[i] = 0;

  return;
}

void start_array_rev(int *array, int length, int weight)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup array to be the starting one in reverse ordering *
\***************************************************************/
{
  int i;

  // error checking
  if (length <= 0)
  {
    printf("ERROR: The length of the array must be positive!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (weight < 0)
  {
    printf("ERROR: The weight for the array must be non-negative!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  for (i = 0; i < length; i++)
    array[i] = 0;
  array[length - 1] = weight;

  return;
}

int advance_array(int *array, int length)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: advances the array to the next one in the ordering     *
\***************************************************************/
{
  //The idea here is that we have a fixed way of cycling through all arrays of a given length with a given total weight (sum of entries).
  //This function just takes the array and tries to move to the next array.
  //Returns 0 upon success, 1 if it is the last array for this weight.

  //We start with [weight, 0, ..., 0] (setup by the calling function).
  //To move to the next, we start at the next to last entry and go towards the first until we find the first nonzero entry (index "target").
  //Once we find one (** if we do **), we move both one from this entry to the entry "target" and also all of the final entry to entry "target."
  //If we do not find a nonzero entry before the final entry, then we are done (return 1).

  //Find the target:
  int target = length - 2;  //length-1 is actually the last entry.
  while (target >= 0)
  { // check to see if target is non-zero
    if (array[target] != 0)
      break;
    target--;
  }

  if (target < 0)  //Didn't find a nonzero entry before the last one - must have been the last array of this length and weight.
    return 1;

  //At this point, we have found the "target" entry and can act on it:
  //move one from this entry to the next:
  array[target]--;
  array[target+1]++;

  //move all of final entry (index length-1) to target+1 - if they are not the same location
  if (target + 1 != length - 1)
  {
    array[target+1] += array[length-1];
    array[length-1] = 0;
  }

  return 0;
}

int advance_array_rev(int *array, int length)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: advances the array to the next one in reverse ordering *
\***************************************************************/
{
  //This function just takes the array and tries to move to the next array.
  //Returns 0 upon success, 1 if it is the last array for this weight.

  //We start with [0,...,0,weight] (setup by the calling function).
  //To move to the next, we start at the second entry and go towards the end until we find the first nonzero entry (index "target").
  //Once we find one (** if we do **), we move both one from this entry to the entry "target" and also all of the final entry to entry "target."
  //If we do not find a nonzero entry before the final entry, then we are done (return 1).

  //Find the target:
  int target = 1;
  while (target < length)
  { // check to see if target is non-zero
    if (array[target] != 0)
      break;
    target++;
  }

  if (target >= length)  //Didn't find a nonzero entry before the last one - must have been the last array of this length and weight.
    return 1;

  //At this point, we have found the "target" entry and can act on it:
  //move one from this entry to the previous:
  array[target]--;
  array[target-1]++;

  //move all of first entry to target-1 - if they are not the same location
  if (target - 1 != 0)
  {
    array[target-1] += array[0];
    array[0] = 0;
  }

  return 0;
}

int array_compare(int *array1, int *array2, int N)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 !(<=), 0 ==, 1 <= & (!=)                    *
* NOTES: determine if array1 <= array2                          *
\***************************************************************/
{
  //We just check to see if array1 is less than or equal to array2, i.e., if array1[i] <= array2[i] for all i.
  //return -1 if array1 !<= array2, 0 if array1 = array2, 1 if array1 < array2.
  int i, ret_val = 0;

  for (i = 0; i < N; i++)
    if (array1[i] > array2[i])
    { // array1 is not less than or equal to array2 - immediately exit
      ret_val = -1;
      break;
    }
    else if (array1[i] < array2[i])
    { // array1 <= array2 and array1 != array2
      ret_val = 1;
    }

  return ret_val;
}

int find_partial_index(int *mon_array, int *diff_array, int func_num, int numFuncs, int numVars, int *diff_size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determine the index associated w/ diff_array-mon_array *
\***************************************************************/
{ 
  int i, retVal = 0, norm = 0, *diff = (int *)bmalloc(numVars * sizeof(int)), *test = (int *)bmalloc(numVars * sizeof(int));

  // find (diff_array - mon_array) and its norm
  for (i = 0; i < numVars; i++)
  { // setup diff = diff_array - mon_array
    diff[i] = diff_array[i] - mon_array[i];
    norm += diff[i];
  }

  // now we move past everything of smaller norm
  for (i = 0; i < norm; i++)
    retVal += numFuncs * diff_size[i];

  // move past the functions before this one
  retVal += func_num * diff_size[norm];

  // setup test to be the first array of this norm
  start_array(test, numVars, norm);

  // increment test until it is equal to diff
  while (array_compare(test, diff, numVars) != 0)
  { // increment test
    advance_array(test, numVars);
    // increment retVal
    retVal++;
  }

  // clear memory
  free(diff); free(test);

  return retVal;
}


