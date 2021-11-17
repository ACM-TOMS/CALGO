// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include <stdio.h>
#include <math.h>
#include "bertini.h"
#include "cascade.h"

// these are functions that are needed to create the mhom start points - local to here
int generateFromPartition_d(preproc_data *PPD, comp_d **coeff, mat_d patchCoeff, int *good_loc, int *P, int **mhomDeg, FILE *OUT);
int generateFromPartition_mp(preproc_data *PPD, comp_mp **coeff, mat_mp patchCoeff, int *good_loc, int *P, int **mhomDeg, FILE *OUT);

void TDstartMaker_d(int *degs, int num_funcs)
{ // total degree start point maker - output to 'nonhom_start'
  FILE *OUT;
  int  i, j, index, Done, total_deg = 1;
  int *curr_degs;

  OUT = fopen("nonhom_start", "w");

  for (i = 0; i < num_funcs; i++)
    total_deg = total_deg * degs[i];

  // Print out number of startpoints.
  fprintf(OUT, "%d\n\n", total_deg);
  
  // We'll use curr_degs to move through all possible degrees:
  curr_degs = (int *)bcalloc(num_funcs, sizeof(int));

  for (i = 0; i < num_funcs; i++)
    curr_degs[i] = 0;

  for (i = 0; i < total_deg; i++)
  { // Make start point for the current degrees
    for (j = 0; j < num_funcs; j++)
      fprintf(OUT, "%.15e %.15e;\n", cos(2*M_PI*curr_degs[j]/degs[j]), sin(2*M_PI*curr_degs[j]/degs[j]));
    fprintf(OUT, "\n");
          
    // Update curr_degs
    index = num_funcs - 1;
    Done = 0;
    while (!Done)
    {
      if (curr_degs[index] < degs[index] - 1)
      {
        curr_degs[index]++;
        for (j = index + 1; j < num_funcs; j++)
          curr_degs[j] = 0;
        Done = 1;
      }
      else
      {
        index--;
        if (index < 0)
          Done = 1;
      }
    }
  }

  fclose(OUT);
  free(curr_degs);

  return;
}

void setupMHstructs(preproc_data *PPD, int **func_gp_count, int ***mhomDeg, FILE *degIN)
{ // read in mhomDeg & func_gp_count from degIN
  int i, j, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp, *fn = *func_gp_count, **deg = *mhomDeg;

  for (i = 0; i < PPD->num_funcs; i++)
  { // initialize fn[i] to 0
    fn[i] = 0;
    // setup deg[i]
    deg[i] = (int *)bmalloc(total_var_gps * sizeof(int));
    for (j = 0; j < total_var_gps; j++)
    {
      fscanf(degIN, "%d\n", &deg[i][j]);
      if (deg[i][j] > 0)
        fn[i]++;
    }
    fscanf(degIN, "\n");
  }

  fn = NULL;
  deg = NULL;

  return;
}

int twoHomStartMaker_d(int *func_gp_count, int **mhomDeg, preproc_data *PPD, int *P, comp_d **coeff, mat_d patchCoeff, FILE *OUT)
{
  int i, j, count = 0, numOnly1 = 0, numOnly2 = 0, numBoth = 0, num_funcs = PPD->num_funcs, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;

  // make sure that we only have 2 variable groups!
  if (total_var_gps != 2)
  {
    printf("ERROR: The number of variable groups is incorrect!\n");
    bexit(ERROR_CONFIGURATION);
  }

  int funcCount[2] = {PPD->size[0] + PPD->type[0] - 1, PPD->size[1] + PPD->type[1] - 1};
  int *only1 = NULL, *only2 = NULL, *both = NULL, *both_curr_loc = NULL;
  int *curr_loc = (int *)bmalloc(num_funcs * sizeof(int));

  // find the number that only involve the 1st or 2nd variable group and the number that involve both and save the function number
  for (i = 0; i < num_funcs; i++)
  {
    if (mhomDeg[i][1] == 0)
    { // this only involves the 1st var group
      only1 = (int *)brealloc(only1, (numOnly1 + 1) * sizeof(int));
      only1[numOnly1] = i;

      numOnly1++;
    }
    else if (mhomDeg[i][0] > 0)
    { // this involves both var groups
      both = (int *)brealloc(both, (numBoth + 1) * sizeof(int));
      both[numBoth] = i;

      numBoth++;
    }
    else
    { // this involves only 2nd var group
      only2 = (int *)brealloc(only2, (numOnly2 + 1) * sizeof(int));
      only2[numOnly2] = i;

      numOnly2++;
    }
  }

  // make sure we have the correct numbers
  if (numOnly1 + numBoth < funcCount[0] || numOnly2 + numBoth < funcCount[1] || funcCount[0] < numOnly1 || funcCount[1] < numOnly2)
  {
    printf("ERROR: Based on the degrees, the system has no isolated solutions!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // update funcCount
  funcCount[0] -= numOnly1;
  funcCount[1] -= numOnly2;

  // initialize both_curr_loc
  both_curr_loc = (int *)bmalloc(numBoth * sizeof(int));
  for (i = 0; i < numBoth; i++)
    if (i < funcCount[0])
      both_curr_loc[i] = 0;
    else
      both_curr_loc[i] = 1;

  // initialize curr_loc
  for (i = 0; i < num_funcs; i++)
    curr_loc[i] = 0;
  for (i = 0; i < numOnly2; i++)
    curr_loc[only2[i]] = 1;
  for (i = 0; i < numBoth; i++)
    curr_loc[both[i]] = both_curr_loc[i];

  // loop over the ones that will generate start points
  while (1)
  { // we have a good partition - so we need to move through all of the degrees relating to this partition
    count += generateFromPartition_d(PPD, coeff, patchCoeff, curr_loc, P, mhomDeg, OUT);

    // find first 0 on the left 
    for (i = 0; i < numBoth; i++)
      if (both_curr_loc[i] == 0)
        break;

    // find the next 1 
    for (j = i + 1; j < numBoth; j++)
      if (both_curr_loc[j] == 1)
        break;

    // see if we have something to move
    if (j >= numBoth)
      break;

    // change j & j - 1
    both_curr_loc[j - 1] = 1;
    both_curr_loc[j] = 0;

    // update curr_loc to the next partition
    curr_loc[both[j - 1]] = both_curr_loc[j - 1];
    curr_loc[both[j]] = both_curr_loc[j];
  }

  // free memory
  free(curr_loc);
  free(both_curr_loc);
  free(only1);
  free(only2);
  free(both);

  return count;
}






int choose_col_in_row(int **mhomDeg, int *var_gp_ctr, int row, int column, int m)
{
    int col = column + 1;  //We assume the current column is done and we need to increment by at least one.
    int done = 0;
    
    if (col-1 > -1)  //If we are coming off of a good partition, we need to remember to increment var_gp_ctr as we move away from that column.
      var_gp_ctr[col-1] = var_gp_ctr[col-1] + 1;
    
    while (!done)
    {
        if (col == m) //got to the right end of the degree matrix!
            done = 1;
        else
        {
            if ((mhomDeg[row][col] == 0) || (var_gp_ctr[col] == 0))  //bad choice, either way!
                col = col + 1;
            else  //means mhomDeg[row][col] > 0, col <= m, and var_gp_ctr[col] > 0 --> good choice of column!
            {
                done = 1;
                var_gp_ctr[col] = var_gp_ctr[col]-1;
            }
        }
    }
    return col;
}


int multiHomStartMaker2_d(int *func_gp_count, int **mhomDeg, preproc_data *PPD, int *P, comp_d **coeff, mat_d patchCoeff, FILE *OUT)
{ //Finds good multihom startpoints intelligently, without going through all possible combinations (like multiHomStartMaker_d()).
  
  int i, row=0, m=PPD->num_hom_var_gp+PPD->num_var_gp;  //m=# var gps; we cycle through row and col in degree table mhomDeg.
  int count=0;
  int old_current_part_row = -1;
  int bad_choice = 0;
    
  int *current_part = (int *)bmalloc(PPD->num_funcs * sizeof(int));  //holds our current choice of partition
//  int *curr_func_gp_count = (int *)bmalloc(PPD->num_funcs * sizeof(int));  //stores # var gps with nonzero degree for each func
//  int **func_var_gps = (int **)bmalloc(PPD->num_funcs * sizeof(int *));  //stores which variable groups have nonzero degree (col) for each function (row)
  int *var_gp_ctr = (int *)bmalloc(m * sizeof(int));  //tracks # of each var gp type chosen so far
  int *K = (int *)bmalloc(m * sizeof(int));  //size of each variable group
    
  for (i = 0; i < m; i++)  // cycle through all var gps
  {
      K[i] = PPD->size[i] + PPD->type[i] - 1;  //# vars in gp i
      var_gp_ctr[i] = K[i]; //start with var_gp_ctr = K and decrement as we go
  }
    
    //printf("K %d %d\nvar_gp_ctr %d %d\n", K[0], K[1], var_gp_ctr[0], var_gp_ctr[1]);
    
  // setup func_var_gps, curr_func_gp_count and initialize curr_loc
  for (i = 0; i < PPD->num_funcs; i++)  //for each function
  {
/*
 func_var_gps[i] = (int *)bmalloc(func_gp_count[i] * sizeof(int));  //preps row i of func_var_gps matrix to take nonzero degrees (and only the nonzero degrees!)
    curr_func_gp_count[i] = 0;    //stores # var gps in which function i has nonzero degree
    for (j = 0; j < m; j++)  //for each var gp
      if (mhomDeg[i][j] > 0)  //if degree of func i in var gp j > 0, we record that:
      {
        func_var_gps[i][curr_func_gp_count[i]] = j;  // we add j to the list of var gps for which func i has nonzero degree
        curr_func_gp_count[i]++;  // we remember that function i has another var gp with nonzero degree
      }
*/
      current_part[i] = -1;
  }

    
  while (row > -1)  // Algorithm will move up and down rows, kicking out to row=-1 at end
  {
    old_current_part_row = current_part[row];  //Hang on to previous choice of column for this row, in case we are done with this row.
    current_part[row] = choose_col_in_row(mhomDeg,var_gp_ctr,row,current_part[row],m);  //Pick next column (var gp) for the current row (func)
      
    if (current_part[row] == m) // means we have exhausted all good columns for the current row, so we go back up a row
    {
//      var_gp_ctr[old_current_part_row] = var_gp_ctr[old_current_part_row]+1;  //done with this row, so increment counter from previously chosen column for this row
      row = row-1;  //go back up a row
      bad_choice = 1;
    }
    else  //found a good choice of column for this row!
    {
      row = row+1;  //move on to next row!
      if (row < PPD->num_funcs)
        current_part[row] = -1; //since we are starting a new row, we start with the left-most entry (choose_col_in_row() first increments col)
    }
     
    if ((row == PPD->num_funcs) && (!bad_choice))  //We have reached the final row with a good partition!
    {
      count += generateFromPartition_d(PPD, coeff, patchCoeff, current_part, P, mhomDeg, OUT);  //This line generates all startpoints for this partition
     // var_gp_ctr[current_part[row-1]] = var_gp_ctr[current_part[row-1]] + 1;  //done with this row, so increment counter from previously chosen column for this row
      row = row - 1; //put the counter back on the last row to try to move to the next column
    }
    bad_choice=0;
  }
    
  // free memory
  free(current_part);
  free(var_gp_ctr);
  free (K);
//  free(curr_func_gp_count);
//  for (i = PPD->num_funcs - 1; i >= 0; i--)
//    free(func_var_gps[i]);
//  free(func_var_gps);

  return count;
}







int multiHomStartMaker_d(int *func_gp_count, int **mhomDeg, preproc_data *PPD, int *P, comp_d **coeff, mat_d patchCoeff, FILE *OUT)
{
  int i, j, count = 0, cont = 1, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  int *curr_loc = (int *)bmalloc(PPD->num_funcs * sizeof(int)), *curr_func_gp_count = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int **func_var_gps = (int **)bmalloc(PPD->num_funcs * sizeof(int *));

  // setup func_var_gps, curr_func_gp_count and initialize curr_loc
  for (i = 0; i < PPD->num_funcs; i++)
  { // setup func_var_gps
    func_var_gps[i] = (int *)bmalloc(func_gp_count[i] * sizeof(int));
    curr_func_gp_count[i] = 0;
    for (j = 0; j < total_var_gps; j++)
      if (mhomDeg[i][j] > 0)
      {
        func_var_gps[i][curr_func_gp_count[i]] = j;
        curr_func_gp_count[i]++;
      }
    curr_func_gp_count[i] = 0;
    // initialize curr_loc[i]
    curr_loc[i] = func_var_gps[i][curr_func_gp_count[i]];
  }

    
    
  while (cont)
  { // check to see if the partition is good
    if (!checkLoc(curr_loc, P, mhomDeg, PPD))
    { // we have a good partition - so we need to move through all of the degrees relating to this partition
      count += generateFromPartition_d(PPD, coeff, patchCoeff, curr_loc, P, mhomDeg, OUT);
    }

    // update to next partition
    for (j = PPD->num_funcs - 1; j >= 0; j--)
    { // check to see if we are at the top

      if (curr_func_gp_count[j] == func_gp_count[j] - 1)
      { // set to 0
        curr_func_gp_count[j] = 0;
        curr_loc[j] = func_var_gps[j][curr_func_gp_count[j]];
      }
      else
      { // increment this position and exit this for loop
        curr_func_gp_count[j]++;
        curr_loc[j] = func_var_gps[j][curr_func_gp_count[j]];
        j = -10;
      }
    }
    // check to see if we need to continue
    if (j == -1)
    { // all were set back to zero so we are done
      cont = 0;
    }
  }

  // free memory
  free(curr_loc);
  free(curr_func_gp_count);
  for (i = PPD->num_funcs - 1; i >= 0; i--)
    free(func_var_gps[i]);
  free(func_var_gps);

  return count;
}

void MHstartMaker_d(preproc_data *PPD, int *P, comp_d **coeff, mat_d patchCoeff)
{
  int i, count = 0;
  int *func_gp_count = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int **mhomDeg = (int **)bmalloc(PPD->num_funcs * sizeof(int *));

  char ch;
  FILE *OUT = fopen("start_temp", "w");
  FILE *degIN = fopen("deg.out", "r");
  if (degIN == NULL)
  {
    printf("ERROR: 'deg.out' does not exist!!!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  
  // setup mhomDeg & func_gp_count
  setupMHstructs(PPD, &func_gp_count, &mhomDeg, degIN);

  // close degIN
  fclose(degIN);

  // generate the start points
  count = multiHomStartMaker2_d(func_gp_count, mhomDeg, PPD, P, coeff, patchCoeff, OUT);

  // close OUT
  fclose(OUT);

  // now we need to make 'start'
  degIN = fopen("start_temp", "r");
  OUT = fopen("start", "w");
  fprintf(OUT, "%d\n\n", count);
  ch = fgetc(degIN);
  while (ch != EOF)
  {
    fputc(ch, OUT);
    ch = fgetc(degIN);
  }
  fclose(OUT);
  fclose(degIN);
  remove("start_temp");

  // free memory
  free(func_gp_count);
  for (i = PPD->num_funcs - 1; i >= 0; i--)
    free(mhomDeg[i]);
  free(mhomDeg);

  return;
}

int generateFromPartition_d(preproc_data *PPD, comp_d **coeff, mat_d patchCoeff, int *good_loc, int *P, int **mhomDeg, FILE *OUT)
{ // returns number of start points created
  int i, j, k, rowNum, num_vars, cont = 1, total_deg = 0, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  int *curr_degs = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int *deg_offset = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  mat_d A;
  vec_d b, x;

  // find num_vars
  num_vars = 0;
  for (i = 0; i < total_var_gps; i++)
    num_vars += PPD->size[i] + PPD->type[i];

  // initialize A, b & x
  init_mat_d(A, PPD->num_funcs + total_var_gps, num_vars);
  init_vec_d(b, PPD->num_funcs + total_var_gps);
  init_vec_d(x, num_vars);
  A->rows = b->size = PPD->num_funcs + total_var_gps; // should be square!
  A->cols = x->size = num_vars;

  for (i = 0; i < PPD->num_funcs; i++)
  {
    curr_degs[i] = deg_offset[i] = 0;
    if (i > 0)
    {
      deg_offset[i] = deg_offset[i-1];
      for (j = 0; j < total_var_gps; j++)
        deg_offset[i] += mhomDeg[P[i-1]][j];
    }
  }

  // setup A & b
  for (i = 0; i < PPD->num_funcs; i++)
  { // top of b is 0
    set_zero_d(&b->coord[i]);
  }
  for (i = 0; i < patchCoeff->rows; i++)
  { // bottom of b is 1
    set_double_d(&b->coord[i + PPD->num_funcs], 1.0, 0.0);
    for (j = 0; j < num_vars; j++)
    { // put patches at bottom of A
      set_d(&A->entry[i + PPD->num_funcs][j], &patchCoeff->entry[i][j]);
    }
  }

  // loop through all of the possible degree configurations
  while (cont)
  { // increment the number of start points made
    total_deg++;

    // Make start point for the current degrees
    for (j = 0; j < PPD->num_funcs; j++)
    { // find the row number of the linear described by good_loc[j]
      rowNum = deg_offset[j] + curr_degs[j]; 
      for (k = 0; k < good_loc[j]; k++)
        rowNum += mhomDeg[P[j]][k];

      // copy this to A
      for (k = 0; k < num_vars; k++)
      {
        set_d(&A->entry[j][k], coeff[rowNum][k]);
      }
    }
    // perform matrix solve and print output
    matrixSolve_d(x, A, b);
    // print to OUT
    for (j = 0; j < x->size; j++)
      fprintf(OUT, "%.15e %.15e;\n", x->coord[j].r, x->coord[j].i);
    fprintf(OUT, "\n");

    // update curr_degs
    for (j = PPD->num_funcs - 1; j >= 0; j--)
    { // set to 0 and continue loop if at the top
      if (curr_degs[j] == mhomDeg[P[j]][good_loc[j]] - 1)
        curr_degs[j] = 0;
      else
      { // increment this poisition and exit this for loop
        curr_degs[j] ++;
        j = -10;
      }
    }
    // check to see if we need to continue
    if (j == -1)
    { // all were set back to zero so we are done
      cont = 0;
    }
  }
  // clear memory
  free(curr_degs);
  free(deg_offset);
  clear_mat_d(A);
  clear_vec_d(b); clear_vec_d(x);

  // return the number of points printed to OUT
  return total_deg;
}

int checkLoc(int *loc, int *P, int **mhomDeg, preproc_data *PPD)
{
  int i, num_var_gps = PPD->num_var_gp + PPD->num_hom_var_gp;
  int *count = (int *)bmalloc(num_var_gps * sizeof(int));
  for (i = 0; i < num_var_gps; i++)
    count[i] = 0;
 
  // count the number of entries in each spot
  for (i = 0; i < PPD->num_funcs; i++)
  {
    if (mhomDeg[P[i]][loc[i]] == 0) // meaning that the 'loc[i]' variable group in the ith function has no degree!
    {
      free(count);
      return -1;
    }
 
    count[loc[i]]++;
  }

  // check to make sure that the number of entries in each spot == var_gp_sizes - 1
  for (i = 0; i < num_var_gps; i++)
    if (count[i] != (PPD->size[i] + PPD->type[i] - 1))
    {
      free(count);
      return -1;
    }
 
  free(count);
  return 0;
}

/////////////////// MP VERSIONS ///////////////////////////

void TDstartMaker_mp(int *degs, int num_funcs)
{ // total degree start point maker - output to 'nonhom_start'
  FILE *OUT;
  int  i, j, index, Done, total_deg = 1;
  int *curr_degs;
  comp_mp tempComp;
  mpf_t two_pi, tempMPF;

  // initialize MP
  init_mp(tempComp);
  mpf_init(two_pi); mpf_init(tempMPF);
  mpfr_const_pi(two_pi, __gmp_default_rounding_mode);
  mpf_mul_ui(two_pi, two_pi, 2);
  mpfr_free_cache(); // free the cache needed to create PI

  OUT = fopen("nonhom_start", "w");

  for (i = 0; i < num_funcs; i++)
    total_deg = total_deg * degs[i];

  // Print out number of startpoints.
  fprintf(OUT, "%d\n\n", total_deg);

  // We'll use curr_degs to move through all possible degrees:
  curr_degs = (int *)bcalloc(num_funcs, sizeof(int));

  for (i = 0; i < num_funcs; i++)
    curr_degs[i] = 0;

  for (i = 0; i < total_deg; i++)
  {
    // Make start point for the current degrees
    for (j = 0; j < num_funcs; j++)
    { // find the angle - 2*PI*curr[j]/deg[j]
      mpf_mul_ui(tempMPF, two_pi, curr_degs[j]);
      mpf_div_ui(tempMPF, tempMPF, degs[j]);
      // find sin & cos of angle
      mpfr_sin_cos(tempComp->i, tempComp->r, tempMPF, __gmp_default_rounding_mode);

      // print to OUT
      print_mp(OUT, 0, tempComp);
      fprintf(OUT, ";\n");
    }
    fprintf(OUT, "\n");

    // Update curr_degs
    index = num_funcs - 1;
    Done = 0;
    while (!Done)
    {
      if (curr_degs[index] < degs[index] - 1)
      {
        curr_degs[index]++;
        for (j = index + 1; j < num_funcs; j++)
          curr_degs[j] = 0;
        Done = 1;
      }
      else
      {
        index--;
        if (index < 0)
          Done = 1;
      }
    }
  }

  fclose(OUT);
  free(curr_degs);

  // clear MP
  clear_mp(tempComp);
  mpf_clear(two_pi); mpf_clear(tempMPF);

  return;
}

int twoHomStartMaker_mp(int *func_gp_count, int **mhomDeg, preproc_data *PPD, int *P, comp_mp **coeff, mat_mp patchCoeff, FILE *OUT)
{
  int i, j, count = 0, numOnly1 = 0, numOnly2 = 0, numBoth = 0, num_funcs = PPD->num_funcs, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;

  // make sure that we only have 2 variable groups!
  if (total_var_gps != 2)
  {
    printf("ERROR: The number of variable groups is incorrect!\n");
    bexit(ERROR_CONFIGURATION);
  }

  int funcCount[2] = {PPD->size[0] + PPD->type[0] - 1, PPD->size[1] + PPD->type[1] - 1};
  int *only1 = NULL, *only2 = NULL, *both = NULL, *both_curr_loc = NULL;
  int *curr_loc = (int *)bmalloc(num_funcs * sizeof(int));

  // find the number that only involve the 1st or 2nd variable group and the number that involve both and save the function number
  for (i = 0; i < num_funcs; i++)
  {
    if (mhomDeg[i][1] == 0)
    { // this only involves the 1st var group
      only1 = (int *)brealloc(only1, (numOnly1 + 1) * sizeof(int));
      only1[numOnly1] = i;

      numOnly1++;
    }
    else if (mhomDeg[i][0] > 0)
    { // this involves both var groups
      both = (int *)brealloc(both, (numBoth + 1) * sizeof(int));
      both[numBoth] = i;

      numBoth++;
    }
    else
    { // this involves only 2nd var group
      only2 = (int *)brealloc(only2, (numOnly2 + 1) * sizeof(int));
      only2[numOnly2] = i;

      numOnly2++;
    }
  }

  // make sure we have the correct numbers
  if (numOnly1 + numBoth < funcCount[0] || numOnly2 + numBoth < funcCount[1] || funcCount[0] < numOnly1 || funcCount[1] < numOnly2)
  {
    printf("ERROR: Based on the degrees, the system has no isolated solutions!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // update funcCount
  funcCount[0] -= numOnly1;
  funcCount[1] -= numOnly2;

  // initialize both_curr_loc
  both_curr_loc = (int *)bmalloc(numBoth * sizeof(int));
  for (i = 0; i < numBoth; i++)
    if (i < funcCount[0])
      both_curr_loc[i] = 0;
    else
      both_curr_loc[i] = 1;

  // initialize curr_loc
  for (i = 0; i < num_funcs; i++)
    curr_loc[i] = 0;
  for (i = 0; i < numOnly2; i++)
    curr_loc[only2[i]] = 1;
  for (i = 0; i < numBoth; i++)
    curr_loc[both[i]] = both_curr_loc[i];

  // loop over the ones that will generate start points
  while (1)
  { // we have a good partition - so we need to move through all of the degrees relating to this partition
    count += generateFromPartition_mp(PPD, coeff, patchCoeff, curr_loc, P, mhomDeg, OUT);

    // find first 0 on the left
    for (i = 0; i < numBoth; i++)
      if (both_curr_loc[i] == 0)
        break;

    // find the next 1
    for (j = i + 1; j < numBoth; j++)
      if (both_curr_loc[j] == 1)
        break;

    // see if we have something to move
    if (j >= numBoth)
      break;

    // change j & j - 1
    both_curr_loc[j - 1] = 1;
    both_curr_loc[j] = 0;

    // update curr_loc to the next partition
    curr_loc[both[j - 1]] = both_curr_loc[j - 1];
    curr_loc[both[j]] = both_curr_loc[j];
  }

  // free memory
  free(curr_loc);
  free(both_curr_loc);
  free(only1);
  free(only2);
  free(both);

  return count;
}

int multiHomStartMaker_mp(int *func_gp_count, int **mhomDeg, preproc_data *PPD, int *P, comp_mp **coeff, mat_mp patchCoeff, FILE *OUT)
{
  int i, j, count = 0, cont = 1, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  int *curr_loc = (int *)bmalloc(PPD->num_funcs * sizeof(int)), *curr_func_gp_count = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int **func_var_gps = (int **)bmalloc(PPD->num_funcs * sizeof(int *));

  // setup func_var_gps, curr_func_gp_count and initialize curr_loc
  for (i = 0; i < PPD->num_funcs; i++)
  { // setup func_var_gps
    func_var_gps[i] = (int *)bmalloc(func_gp_count[i] * sizeof(int));
    curr_func_gp_count[i] = 0;
    for (j = 0; j < total_var_gps; j++)
      if (mhomDeg[i][j] > 0)
      {
        func_var_gps[i][curr_func_gp_count[i]] = j;
        curr_func_gp_count[i]++;
      }
    curr_func_gp_count[i] = 0;

    // initialize curr_loc[i]
    curr_loc[i] = func_var_gps[i][curr_func_gp_count[i]];
  }

  while (cont)
  { // check to see if the partition is good
    if (!checkLoc(curr_loc, P, mhomDeg, PPD))
    { // we have a good partition - so we need to move through all of the degrees relating to this partition
      count += generateFromPartition_mp(PPD, coeff, patchCoeff, curr_loc, P, mhomDeg, OUT);
    }

    // update to next partition
    for (j = PPD->num_funcs - 1; j >= 0; j--)
    { // check to see if we are at the top

      if (curr_func_gp_count[j] == func_gp_count[j] - 1)
      { // set to 0
        curr_func_gp_count[j] = 0;
        curr_loc[j] = func_var_gps[j][curr_func_gp_count[j]];
      }
      else
      { // increment this position and exit this for loop
        curr_func_gp_count[j]++;
        curr_loc[j] = func_var_gps[j][curr_func_gp_count[j]];
        j = -10;
      }
    }
    // check to see if we need to continue
    if (j == -1)
    { // all were set back to zero so we are done
      cont = 0;
    }
  }

  // free memory
  free(curr_loc);
  free(curr_func_gp_count);
  for (i = PPD->num_funcs - 1; i >= 0; i--)
    free(func_var_gps[i]);
  free(func_var_gps);

  return count;
}

void MHstartMaker_mp(preproc_data *PPD, int *P, comp_mp **coeff, mat_mp patchCoeff)
{
  int i, count = 0;
  int *func_gp_count = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int **mhomDeg = (int **)bmalloc(PPD->num_funcs * sizeof(int *));

  char ch;
  FILE *OUT = fopen("start_temp", "w");
  FILE *degIN = fopen("deg.out", "r");
  if (degIN == NULL)
  {
    printf("ERROR: 'deg.out' does not exist!!!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // setup mhomDeg & func_gp_count
  setupMHstructs(PPD, &func_gp_count, &mhomDeg, degIN);

  // close degIN
  fclose(degIN);

  // generate the start points
  count = multiHomStartMaker_mp(func_gp_count, mhomDeg, PPD, P, coeff, patchCoeff, OUT);

  // close OUT
  fclose(OUT);

  // now we need to make 'start'
  degIN = fopen("start_temp", "r");
  OUT = fopen("start", "w");
  fprintf(OUT, "%d\n\n", count);
  ch = fgetc(degIN);
  while (ch != EOF)
  {
    fputc(ch, OUT);
    ch = fgetc(degIN);
  }
  fclose(OUT);
  fclose(degIN);
  remove("start_temp");

  // free memory
  free(func_gp_count);
  for (i = PPD->num_funcs - 1; i >= 0; i--)
    free(mhomDeg[i]);
  free(mhomDeg);

  return;
}

int generateFromPartition_mp(preproc_data *PPD, comp_mp **coeff, mat_mp patchCoeff, int *good_loc, int *P, int **mhomDeg, FILE *OUT)
{ // returns number of start points created
  int i, j, k, rowNum, num_vars, cont = 1, total_deg = 0, total_var_gps = PPD->num_hom_var_gp + PPD->num_var_gp;
  int *curr_degs = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  int *deg_offset = (int *)bmalloc(PPD->num_funcs * sizeof(int));
  mat_mp A;
  vec_mp b, x;

  // find num_vars
  num_vars = 0;
  for (i = 0; i < total_var_gps; i++)
    num_vars += PPD->size[i] + PPD->type[i];

  // initialize A, b, & x
  init_mat_mp(A, PPD->num_funcs + total_var_gps, num_vars);
  init_vec_mp(b, PPD->num_funcs + total_var_gps);
  init_vec_mp(x, num_vars);
  A->rows = b->size = PPD->num_funcs + total_var_gps; // should be square!
  A->cols = x->size = num_vars;

  for (i = 0; i < PPD->num_funcs; i++)
  {
    curr_degs[i] = deg_offset[i] = 0;
    if (i > 0)
    {
      deg_offset[i] = deg_offset[i-1];
      for (j = 0; j < total_var_gps; j++)
        deg_offset[i] += mhomDeg[P[i-1]][j];
    }
  }

  // setup A & b
  for (i = 0; i < PPD->num_funcs; i++)
  { // top of b is 0
    set_zero_mp(&b->coord[i]);
  }
  for (i = 0; i < patchCoeff->rows; i++)
  { // bottom of b is 1
    mpf_set_ui(b->coord[i+PPD->num_funcs].r, 1);
    mpf_set_ui(b->coord[i+PPD->num_funcs].i, 0);
    for (j = 0; j < num_vars; j++)
    { // put patches at bottom of A
      set_mp(&A->entry[i + PPD->num_funcs][j], &patchCoeff->entry[i][j]);
    }
  }

  // loop through all of the possible degree configurations
  while (cont)
  { // increment the number of start points made
    total_deg++;

    // Make start point for the current degrees
    for (j = 0; j < PPD->num_funcs; j++)
    { // find the row number of the linear described by good_loc[j]
      rowNum = deg_offset[j] + curr_degs[j];
      for (k = 0; k < good_loc[j]; k++)
        rowNum += mhomDeg[P[j]][k];

      // copy this to A
      for (k = 0; k < num_vars; k++)
      {
        set_mp(&A->entry[j][k], coeff[rowNum][k]);
      }
    }

    // perform matrix solve and print output
    matrixSolve_mp(x, A, b);

    // print to OUT
    for (j = 0; j < x->size; j++)
    {
      print_mp(OUT, 0, &x->coord[j]);
      fprintf(OUT, ";\n");
    }
    fprintf(OUT, "\n");

    // update curr_degs
    for (j = PPD->num_funcs - 1; j >= 0; j--)
    { // set to 0 and continue loop if at the top
      if (curr_degs[j] == mhomDeg[P[j]][good_loc[j]] - 1)
        curr_degs[j] = 0;
      else
      { // increment this position and exit this for loop
        curr_degs[j]++;
        j = -10;
      }
    }
    // check to see if we need to continue
    if (j == -1)
    { // all were set back to 0 so we are done
      cont = 0;
    }
  }
  free(curr_degs);
  free(deg_offset);

  clear_mat_mp(A);
  clear_vec_mp(b); 
  clear_vec_mp(x);

  // return the number of points printed to OUT
  return total_deg;
}

