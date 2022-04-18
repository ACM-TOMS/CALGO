// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bertini.h"
#include "cascade.h"
#include "regeneration.h"
#include "parallel.h"

//  USAGE:  This file will contain functions to create the output of the solutions found in a given Bertini run.

void multiplicityPointSummary(int numStartPts, post_process_t *endPoints, int finiteToggle, int regenToggle, int eqbyeqToggle);
void singularPointSummary(int numStartPts, post_process_t *endPoints, int finiteToggle, int regenToggle, int eqbyeqToggle);
void zeroDimOutputChart(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, char **name_table, int num_sols, int num_vars, preproc_data *PPD, double maxNorm, double realTol, double finalTol, double maxCondNum, int regenToggle, int eqbyeqToggle);
void findRealSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double realTol);
int checkForReal_d(point_d Pt, double realTol);
int checkForReal_mp(point_mp Pt, double realTol);
void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle);
void findMultSol(post_process_t *endPoints, int num_sols, int num_vars, preproc_data *PPD, double finalTol);
void findFiniteSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxNorm);
void createRawSoln(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD);

void printFileHeader(FILE *OUT, int num_vars, preproc_data *PPD, char **name_table)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header of the file                          *
\***************************************************************/
{
  int i;

  fprintf(OUT, "Number of variables: %d\n", num_vars - PPD->num_var_gp);
  fprintf(OUT, "Variables: ");
  for (i = PPD->num_var_gp; i < num_vars; i++)
    fprintf(OUT, " %s", name_table[i]);
  fprintf(OUT, "\nRank: %d\n", num_vars - PPD->num_hom_var_gp - PPD->num_var_gp);

  return;
}

void printMainDataPointHeader_d(FILE *OUT, int sol_num, int path_num, double cond_num, double function_resid, double newton_resid, double final_T, double function_error, int prec, double first_increase, int cycle_num, int success, int origErrorIsInf, double origErrorEst)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header of a path for main_data              *
\***************************************************************/
{
  fprintf(OUT, "-------------------------\nSolution %d (path number %d", sol_num, path_num);
  if (success == retVal_sharpening_failed || success == retVal_sharpening_singular_endpoint)
    fprintf(OUT, "&");
  else if (success != 1)
    fprintf(OUT, "#");
  fprintf(OUT, ")\n");
  fprintf(OUT, "Estimated condition number: %.15e\nFunction residual: %.15e\nLatest Newton residual: %.15e\n", cond_num, function_resid, newton_resid);
  fprintf(OUT, "T value at final sample point: %.15e\nMaximum precision utilized: %d\n", final_T, prec);
  fprintf(OUT, "T value of first precision increase: %.15e\n", first_increase);
  fprintf(OUT, "Accuracy estimate, internal coordinates (difference of last two endpoint estimates):  %.15e\n", function_error);
  fprintf(OUT, "Accuracy estimate, user's coordinates (after dehomogenization, if applicable):  ");
  if (origErrorIsInf)
    fprintf(OUT, "infinity\n");
  else
    fprintf(OUT, "%.15e\n", origErrorEst);
  fprintf(OUT, "Cycle number: %d\n", cycle_num);

  return;
}

void printMainDataPointHeader_mp(FILE *OUT, int sol_num, int path_num, double cond_num, mpf_t function_resid, mpf_t newton_resid, double final_T, double function_error, int prec, double first_increase, int cycle_num, int success, int origErrorIsInf, double origErrorEst)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header of a path for main_data              *
\***************************************************************/
{
  fprintf(OUT, "-------------------------\nSolution %d (path number %d", sol_num, path_num);
  if (success == retVal_sharpening_failed || success == retVal_sharpening_singular_endpoint)
    fprintf(OUT, "&");
  else if (success != 1)
    fprintf(OUT, "#");
  fprintf(OUT, ")\n");
  fprintf(OUT, "Estimated condition number: %.15e\nFunction residual: ", cond_num);
  mpf_out_str(OUT, 10, 15, function_resid);
  fprintf(OUT, "\nLatest Newton residual: ");
  mpf_out_str(OUT, 10, 15, newton_resid); 
  fprintf(OUT, "\nT value at final sample point: %.15e\nMaximum precision utilized: %d\n", final_T, prec);
  fprintf(OUT, "T value of first precision increase: %.15e\n", first_increase);
  fprintf(OUT, "Accuracy estimate, internal coordinates (difference of last two endpoint estimates):  %.15e\n", function_error);
  fprintf(OUT, "Accuracy estimate, user's coordinates (after dehomogenization, if applicable):  ");
  if (origErrorIsInf)
    fprintf(OUT, "infinity\n");
  else
    fprintf(OUT, "%.15e\n", origErrorEst);
  fprintf(OUT, "Cycle number: %d\n", cycle_num);

  return;
}

void printMainDataPointFooter(FILE *OUT, int curr_sol, int num_sols, post_process_t *endPoints)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the header of a path for main_data              *
\***************************************************************/
{
  int j;

  fprintf(OUT, "Paths with the same endpoint, to the prescribed tolerance:  ");

  //  Now we compare it to the other points:
  if (endPoints[curr_sol].multiplicity > 1)
  {
    for (j = 0; j < num_sols; j++)
      if (j != curr_sol && endPoints[j].sol_num == endPoints[curr_sol].sol_num)
      {
        fprintf(OUT, "%d", endPoints[j].path_num);
        if (endPoints[j].success == retVal_sharpening_failed || endPoints[j].success == retVal_sharpening_singular_endpoint)
          fprintf(OUT, "&");
        else if (endPoints[j].success != 1)
          fprintf(OUT, "#");

        fprintf(OUT, "  ");
      }
  }

  fprintf(OUT, "\nMultiplicity: %d\n", endPoints[curr_sol].multiplicity);

  return;
}

double local_error_d(int *isInf, double u, double v, double ruv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximate the local error when converting back to    *
*   original coordinates                                        *
\***************************************************************/
{
  double error = 0;

  // verify all are nonnegative
  if (u < 0 || v < 0 || ruv < 0)
  {
    printf("ERROR: The accuracy estimates must be nonnegative!\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (ruv == 0)
  { // local error is zero up to machine precision
    ruv = 1e-16;
  }

  // compare v and ruv
  if (v < ruv)
  { // infinite error
    *isInf = 1;
  } 
  else 
  { // not infinite error
    *isInf = 0; 
    error = (ruv*ruv*(u/v + 1) + ruv*(u+v)) / (v*v - ruv*ruv);
  }

  return error;
}

double local_error_mp(int *isInf, mpf_t u, mpf_t v, mpf_t ruv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: approximate the local error when converting back to    *
*   original coordinates                                        *
\***************************************************************/
{
  double error = 0;
  mpf_t tempRUV, tempMPF1, tempMPF2;

  mpf_init(tempRUV);
  mpf_init(tempMPF1);
  mpf_init(tempMPF2);

  // verify all are nonnegative
  if (mpfr_sgn(u) < 0 || mpfr_sgn(v) < 0 || mpfr_sgn(ruv) < 0)
  {
    printf("ERROR: The accuracy estimates must be nonnegative!\n");
    bexit(ERROR_CONFIGURATION);
  }

  if (mpfr_zero_p(ruv))
  { // local error is zero up to machine precision
    int num_digits = prec_to_digits(mpfr_get_prec(u));
    mpfr_ui_pow_ui(tempRUV, 10, num_digits, __gmp_default_rounding_mode);
    mpf_ui_div(tempRUV, 1, tempRUV);
  }
  else
    mpf_set(tempRUV, ruv);

  // compare v and ruv
  if (mpf_cmp(v, tempRUV) < 0)
  { // infinite error
    *isInf = 1;
  }
  else
  { // not infinite error
    *isInf = 0;

    // ruv*ruv*(u/v + 1)
    mpf_div(tempMPF1, u, v);
    mpf_set_ui(tempMPF2, 1);
    mpf_add(tempMPF1, tempMPF1, tempMPF2);
    mpf_mul(tempMPF1, tempRUV, tempMPF1);
    mpf_mul(tempMPF1, tempRUV, tempMPF1);

    // ruv*(u+v)
    mpf_add(tempMPF2, u, v);
    mpf_mul(tempMPF2, tempRUV, tempMPF2);

    // ruv*ruv*(u/v + 1) + ruv*(u+v)
    mpf_add(tempMPF1, tempMPF1, tempMPF2);

    // v*v  - ruv*ruv
    mpf_mul(tempMPF2, v, v);
    mpf_mul(tempRUV, tempRUV, tempRUV);
    mpf_sub(tempMPF2, tempMPF2, tempRUV);

    // error = (ruv*ruv*(u/v + 1) + ruv*(u+v)) / (v*v - ruv*ruv);
    mpf_div(tempMPF1, tempMPF1, tempMPF2);
    error = mpf_get_d(tempMPF1);
  }

  mpf_clear(tempRUV);
  mpf_clear(tempMPF1);
  mpf_clear(tempMPF2);

  return error;
}

void getDehomPoint_comp_d(point_d dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_d *sol, int num_vars, preproc_data *PPD, double accuracyEstimate)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the dehom point                             *
\***************************************************************/
{
  int i, j, largeCoord, currVarGp = 0, var_gp = PPD->num_var_gp, hom_var_gp = PPD->num_hom_var_gp;
  int num_gps = var_gp + hom_var_gp, size = num_vars - var_gp;
  double infNorm, tempD, localErrorEst, numNorm, denomNorm;
  comp_d recip, largest;
 
  // setup the final size 
  change_size_point_d(dehomPoint, size); 
  dehomPoint->size = size;

  // initialize 
  *origErrorEst = *origErrorIsInf = 0;

  if (num_gps == 0)
  { // no variable groups - simply copy the point 
    for (i = 0; i < num_vars; i++)
    {
      set_d(&dehomPoint->coord[i], sol[i]);
    }
    *origErrorEst = accuracyEstimate;
  }
  else
  { // loop over the variable groups
    for (i = 0; i < num_gps; i++) 
    { // determine if this group is affine or homogeneous
      if (PPD->type[i] == 0) 
      { // homogeneous variable group - normalize by largest coordinate
        infNorm = largeCoord = -1;
        set_one_d(largest);
        for (j = 0; j < PPD->size[i]; j++)
        { // compute the modulus of the coordinate
          tempD = d_abs_d(sol[var_gp + j]);
          if (tempD > infNorm)
          {
            infNorm = tempD;
            set_d(largest, sol[var_gp + j]);
            largeCoord = j;
          }
        }

        // setup the norm of the denominator
        denomNorm = infNorm;

        // verify largest is nonzero - should never be zero!
        if (largest->r == 0 && largest->i == 0) 
        { // generate a random perturbation
          get_comp_rand_d(recip);
          mul_rdouble_d(recip, recip, 1e-16);
          recip_d(recip, recip);
        }
        else
        { // reciprocate the largest
          recip_d(recip, largest);
        }

        // normalize the coordinates
        for (j = 0; j < PPD->size[i]; j++)
        {
          if (j == largeCoord)
          { // approximate error, if needed
            if (*origErrorIsInf == 0)
            {
              numNorm = denomNorm;
              localErrorEst = local_error_d(origErrorIsInf, numNorm, denomNorm, accuracyEstimate);

              // update error
              if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
                *origErrorEst = localErrorEst;
            }

            // set to 1
            set_one_d(&dehomPoint->coord[var_gp - PPD->num_var_gp]);
          }
          else
          { // approximate error, if needed
            if (*origErrorIsInf == 0)
            {
              numNorm = d_abs_d(sol[var_gp]); 
              localErrorEst = local_error_d(origErrorIsInf, numNorm, denomNorm, accuracyEstimate);

              // update error
              if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
                *origErrorEst = localErrorEst;
            }
            
            // multiplity by recip
            mul_d(&dehomPoint->coord[var_gp - PPD->num_var_gp], sol[var_gp], recip);
          }

          // update
          var_gp++;
        }
      }
      else
      { // variable group - dehomogenize
        if (sol[currVarGp]->r == 0 && sol[currVarGp]->i == 0)
        { // generate a random perturbation
          get_comp_rand_d(recip);
          mul_rdouble_d(recip, recip, 1e-16);
          recip_d(recip, recip);
        }
        else
        { // reciprocate
          recip_d(recip, sol[currVarGp]);
        }

        // setup the norm of the denominator
        denomNorm = d_abs_d(sol[currVarGp]);

        // dehomogenize the coordinates
        for (j = 0; j < PPD->size[i]; j++) 
        { // approximate error, if needed
          if (*origErrorIsInf == 0)
          {
            numNorm = d_abs_d(sol[var_gp]);
            localErrorEst = local_error_d(origErrorIsInf, numNorm, denomNorm, accuracyEstimate);

            if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
              *origErrorEst = localErrorEst;
          }

          // dehomogenize
          mul_d(&dehomPoint->coord[var_gp - PPD->num_var_gp], sol[var_gp], recip);

          // update
          var_gp++;
        }

        // increment since we have seen a homogenizing coordinate
        currVarGp++;
      }
    }
  }

  return;
}

void countSingRealPoints(int *singReal, int *nonSingReal, int *singNonReal, int *nonSingNonReal, int numStartPts, post_process_t *endPoints, int finiteToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: counts the number of singular & real endpoints         *
\***************************************************************/
{
  int i;

  // initialize all to 0
  *singReal = *nonSingReal = *singNonReal = *nonSingNonReal = 0;

  for (i = 0; i < numStartPts; i++)
    if (endPoints[i].multiplicity > 0 && endPoints[i].isFinite == finiteToggle)
    {
      if (endPoints[i].isReal && endPoints[i].isSing)
        *singReal = *singReal + 1;
      else if (endPoints[i].isReal && !endPoints[i].isSing)
        *nonSingReal = *nonSingReal + 1;
      else if (!endPoints[i].isReal && endPoints[i].isSing)
        *singNonReal = *singNonReal + 1;
      else
        *nonSingNonReal = *nonSingNonReal + 1;
    }

  return;
}

void singularPointSummary(int numStartPts, post_process_t *endPoints, int finiteToggle, int regenToggle, int eqbyeqToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether anything was printed or not            *
* NOTES: prints the singular point summary                      *
\***************************************************************/
{
  int printInfo = 0, singReal = 0, singNonReal = 0, nonSingReal = 0, nonSingNonReal = 0;

  // get a count
  countSingRealPoints(&singReal, &nonSingReal, &singNonReal, &nonSingNonReal, numStartPts, endPoints, finiteToggle);

  if (regenToggle || eqbyeqToggle)
  { // only care about non singular ones
    if (nonSingReal + nonSingNonReal > 0)
      printInfo = 1;
  }
  else
  { // care about all of them
    if (singReal + nonSingReal + singNonReal + nonSingNonReal > 0)
      printInfo = 1;
  }

  // print the information
  if (printInfo)
  {
    if (finiteToggle == 1)
    {
      if (regenToggle == 1 || eqbyeqToggle == 1)
        printf("Non-singular Finite Solution Summary\n");
      else
        printf("\nFinite Solution Summary\n");
    }
    else if (finiteToggle == 0)
      printf("\nInfinite Solution Summary\n");
    else if (finiteToggle == -1)
      printf("\nSolution Summary\n");

    if (!regenToggle && !eqbyeqToggle)
      printf("\nNOTE: nonsingular vs singular is based on condition number and identical endpoints\n");
    else
      printf("\nNOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n");

    printf("\n\t\t| Number of real solns\t|  Number of non-real solns\t|  Total\n");
    printf("------------------------------------------------------------------------------------------\n");
    printf("Non-singular\t|\t%d\t\t|\t\t%d\t\t|   %d\n", nonSingReal, nonSingNonReal, nonSingReal + nonSingNonReal);
    if (!(regenToggle || eqbyeqToggle))
      printf("Singular\t|\t%d\t\t|\t\t%d\t\t|   %d\n", singReal, singNonReal, singReal + singNonReal);
    printf("------------------------------------------------------------------------------------------\n");
    if (regenToggle || eqbyeqToggle)
      printf("Total\t\t|\t%d\t\t|\t\t%d\t\t|   %d\n\n", nonSingReal, nonSingNonReal, nonSingReal + nonSingNonReal);
    else
      printf("Total\t\t|\t%d\t\t|\t\t%d\t\t|   %d\n\n", (nonSingReal + singReal), (nonSingNonReal + singNonReal), (nonSingReal + singReal + nonSingNonReal + singNonReal));
  }

  return;    
}

void multiplicityPointSummary(int numStartPts, post_process_t *endPoints, int finiteToggle, int regenToggle, int eqbyeqToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the multiplicity point summary                  *
\***************************************************************/
{
  int i, j, realCount, nonRealCount, maxMult = 0;

  if (!(regenToggle || eqbyeqToggle))
  { // find the maximum multiplicity of all of them
    for (i = 0; i < numStartPts; i++)
      if (maxMult < endPoints[i].multiplicity && endPoints[i].isFinite == finiteToggle)
        maxMult = endPoints[i].multiplicity;

    if (maxMult > 0)
    {
      if (finiteToggle == 1)
        printf("\nFinite Multiplicity Summary\n");
      else if (finiteToggle == 0)
        printf("\nInfinite Multiplicity Summary\n");
      else if (finiteToggle == -1)
        printf("\nMultiplicity Summary\n");

      printf("\n  Multiplicity\t|  Number of real solns\t|  Number of non-real solns\n"); 
      printf("------------------------------------------------------------------------------------------\n");

      for (i = 1; i <= maxMult; i++)
      {
        realCount = nonRealCount = 0;
 
        for (j = 0; j < numStartPts; j++)
          if (endPoints[j].multiplicity == i && endPoints[j].isFinite == finiteToggle)
          {  
            if (endPoints[j].isReal)
              realCount++;
            else
              nonRealCount++;
          }

        if ((realCount > 0) || (nonRealCount > 0))
          printf("\t%d\t|\t%d\t\t|\t%d\n", i, realCount, nonRealCount); 
      }
      printf("\n");
    }
  }

  if (finiteToggle == 0 || finiteToggle == -1)
  {
    printf("\n------------------------------------------------------------------------------------------\n");
    printf("The following files may be of interest to you:\n\n");
    printf("main_data:             A human-readable version of the solutions - main output file.\n");
    printf("raw_solutions:         A list of the solutions with the corresponding path numbers.\n");
  }

  if (finiteToggle == 0)
  {
    printf("raw_data:              Similar to the previous, but with the points in Bertini's homogeneous\n");
    printf("                         coordinates along with more information about the solutions.\n");
    printf("real_finite_solutions: A list of all real finite solutions.\n");
    printf("finite_solutions:      A list of all finite solutions.\n");
    printf("nonsingular_solutions: A list of all nonsingular solutions.\n");
  
    if (!regenToggle && !eqbyeqToggle)
      printf("singular_solutions:    A list of all singular solutions.\n");

    printf("------------------------------------------------------------------------------------------\n");
  }
  else if (finiteToggle == -1)
  {
    printf("raw_data:              Similar to the previous, but with more information about the solutions.\n");
    printf("real_solutions:        A list of all real solutions.\n");
    printf("nonsingular_solutions: A list of all nonsingular solutions.\n");

    if (!regenToggle)
      printf("singular_solutions:    A list of all singular solutions.\n");

    printf("------------------------------------------------------------------------------------------\n");
  }

  return;
}

int norm_order_d(const void *vp, const void *vq)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: provides the order function for qsort                  *
\***************************************************************/
{
  midpoint_data_d *a = (midpoint_data_d *) vp;
  midpoint_data_d *b = (midpoint_data_d *) vq;
  double x;
 
  x = a->norm - b->norm;

  // return the sign of x 
  if (x > 0)
    return 1;
  else if (x < 0)
    return -1;
  else
    return 0;
}

int norm_order_mp(const void *vp, const void *vq)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: provides the order function for qsort                  *
\***************************************************************/
{
  int sign;
  midpoint_data_mp *a = (midpoint_data_mp *) vp;
  midpoint_data_mp *b = (midpoint_data_mp *) vq;
  mpf_t x;
  mpf_init(x);

  mpf_sub(x, a->norm, b->norm);

  sign = mpfr_sgn(x);
 
  mpf_clear(x);

  // return the sign of x
  return sign;
}


void midpoint_checker(int num_paths, int num_vars, double tol, int *num_crossings)
/***************************************************************\
* RETURN VALUES:                                                *
* NOTES: checks for path crossing                               *
\***************************************************************/
{
  FILE *midIN;
  int i, j, index;

  if (num_paths < 2)
    return;

  if (tol > 1e-14)
  { // check using double precision
    point_d test_vec;
    midpoint_data_d *midpoint_data = NULL;

    init_vec_d(test_vec, 0);

    // read in the data and set up the structure that we'll be using to sort.
    midIN = fopen("midpath_data", "r");
    if (midIN == NULL)
    {
      printf("ERROR: 'midpath_data' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }

    midpoint_data = (midpoint_data_d *)bcalloc(num_paths, sizeof(midpoint_data_d));
    for (i = 0; i < num_paths; i++)
    {
      init_point_d(midpoint_data[i].point, num_vars);
      midpoint_data[i].point->size = num_vars;
      fscanf(midIN, "%d\n", &midpoint_data[i].path_num);
      for (j = 0; j < num_vars; j++)
        fscanf(midIN, "%lf%lf\n", &midpoint_data[i].point->coord[j].r, &midpoint_data[i].point->coord[j].i);
      midpoint_data[i].norm = infNormVec_d(midpoint_data[i].point);
    }

    fclose(midIN);

    // sort
    qsort(midpoint_data, num_paths, sizeof(midpoint_data_d), norm_order_d);
  
    //Finally, we cycle through all of the midpoints and compare to the ones following it whose norms are within the specified tolerance.
    //For each midpoint whose norm is within the specified tolerance, we check to see how close the points are (using the infinity norm).
    //If there's no problem, we move on.  Otherwise, we print a warning (no automatic reruns for now!).
    for (i = num_paths - 1; i >= 0; i--)  //Start at num_paths-1 since there is nothing to compare the last point to!
    {
      index = i + 1;
      while ((index < num_paths) && (midpoint_data[index].norm - midpoint_data[i].norm < tol)) //Terminates if norm of the indexth midpt greater than norm of ith + tol.
      {
        vec_sub_d(test_vec, midpoint_data[index].point, midpoint_data[i].point);
        if (infNormVec_d(test_vec) < tol)
        {
          printf("!!!WARNING!!!  Paths %d and %d may have crossed!!!\n", midpoint_data[i].path_num, midpoint_data[index].path_num);
          *num_crossings = *num_crossings + 1;  //At the end, this will tell us how many crossings happened (NOTE:  WE DON'T ASSUME THIS STARTS AT 0!!!).
        }

        index++;
      }
    }

    // clear memory
    for (i = 0; i < num_paths; i++)
      clear_point_d(midpoint_data[i].point);
    free(midpoint_data);
    clear_vec_d(test_vec);
  }
  else
  { // check using multi precision
    int curr_prec = (int) mpf_get_default_prec(), new_prec = digits_to_prec(2-log10(tol));
    mpf_t norm_diff, tol_mp;
    point_mp test_vec;
    midpoint_data_mp *midpoint_data = NULL;

    // set new precision
    initMP(new_prec);

    mpf_init(norm_diff);
    mpf_init(tol_mp);
    mpf_set_d(tol_mp, tol);
    init_vec_mp(test_vec, 0);

    // read in the data and set up the structure that we'll be using to sort.
    midIN = fopen("midpath_data", "r");
    if (midIN == NULL)
    {
      printf("ERROR: 'midpath_data' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }

    midpoint_data = (midpoint_data_mp *)bcalloc(num_paths, sizeof(midpoint_data_mp));
    for (i = 0; i < num_paths; i++)
    {
      mpf_init(midpoint_data[i].norm);
      init_point_mp(midpoint_data[i].point, num_vars);
      midpoint_data[i].point->size = num_vars;
      fscanf(midIN, "%d\n", &midpoint_data[i].path_num);
      for (j = 0; j < num_vars; j++)
      {
        mpf_inp_str(midpoint_data[i].point->coord[j].r, midIN, 10);
        mpf_inp_str(midpoint_data[i].point->coord[j].i, midIN, 10);
      }

      infNormVec_mp2(midpoint_data[i].norm, midpoint_data[i].point);
    }

    fclose(midIN);

    // sort
    qsort(midpoint_data, num_paths, sizeof(midpoint_data_mp), norm_order_mp);

    //Finally, we cycle through all of the midpoints and compare to the ones following it whose norms are within the specified tolerance.
    //For each midpoint whose norm is within the specified tolerance, we check to see how close the points are (using the infinity norm).
    //If there's no problem, we move on.  Otherwise, we print a warning (no automatic reruns for now!).
    for (i = num_paths - 1; i >= 0; i--)  //Start at num_paths-1 since there is nothing to compare the last point to!
    {
      index = i + 1;
      while (index < num_paths)
      { // compute norm_diff
        mpf_sub(norm_diff, midpoint_data[index].norm, midpoint_data[i].norm);
        if (mpf_cmp(norm_diff, tol_mp) < 0)
        {
          vec_sub_mp(test_vec, midpoint_data[index].point, midpoint_data[i].point);
          infNormVec_mp2(norm_diff, test_vec);
          if (mpf_cmp(norm_diff, tol_mp) < 0)
          {
            printf("!!!WARNING!!!  Paths %d and %d may have crossed!!!\n", midpoint_data[i].path_num, midpoint_data[index].path_num);
            *num_crossings = *num_crossings + 1;  //At the end, this will tell us how many crossings happened (NOTE:  WE DON'T ASSUME THIS STARTS AT 0!!!).
          }

          index++;
        }
        else
        {
          index = num_paths; // break out of loop
        }
      }
    }

    // clear memory
    for (i = 0; i < num_paths; i++)
      clear_point_mp(midpoint_data[i].point);
    free(midpoint_data);
    clear_vec_mp(test_vec);
    mpf_clear(norm_diff);
    mpf_clear(tol_mp);

    // return to old precision
    initMP(curr_prec);
  }

  return;
}

void printFailureSummary(trackingStats *tC, int convergence_failures, int sharpening_failures, int sharpening_singular)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the failure summary                             *
\***************************************************************/
{
  int trueFailures = tC->failures - tC->securityCount;

  printf("\nPaths Tracked: %d\n", tC->numPoints);

  if (tC->securityCount > 0)
    printf(" Truncated infinite paths: %d - try adjusting SecurityMaxNorm or set SecurityLevel to 1 in the input file\n", tC->securityCount);

  if (trueFailures > 0)
  {
    printf(" Number of failures: %d\n", trueFailures);
    if (tC->junkCount > 0)
      printf("     Nonsolutions: %d - see 'nonsolutions'\n", tC->junkCount);
    if (tC->nanCount > 0)
      printf("     Not a valid number: %d - verify the configurations in the input file\n", tC->nanCount);
    if (tC->infCount > 0)
      printf("     Norm exceeded PathTruncationThreshold: %d - try adjusting PathTruncationThreshold in the input file\n", tC->infCount);
    if (tC->sizeCount > 0)
      printf("     Step size too small: %d - try adjusting MinStepSizeBeforeEG or MinStepSizeDuringEG in the input file\n", tC->sizeCount);
    if (tC->PSEGCount > 0)
      printf("     Power series endgame failure: %d\n", tC->PSEGCount);
    if (tC->precCount > 0)
      printf("     Maximum precision reached: %d - try adjusting TrackToldBeforeEG or TrackTolDuringEG in the input file\n", tC->precCount);
    if (tC->cycleCount > 0)
      printf("     Cycle Number too high: %d - try adjusting TrackTolDuringEG in the input file\n", tC->cycleCount);
    if (tC->stepCount > 0)
      printf("     Too many steps: %d - try adjusting MaxNumberSteps in the input file\n", tC->stepCount);
    if (tC->refineCount > 0)
      printf("     Refining failed: %d - try adjusting TrackTolBeforeEG or TrackTolDuringEG in the input file\n", tC->refineCount);
    if (tC->otherCount > 0)
      printf("     Other failures: %d\n", tC->otherCount);
  }

  if (tC->failures > 0)
    printf("   Please see 'failed_paths' for more information about these paths.\n\n");

  if (convergence_failures > 0)
  {
    printf(" Number that failed to converge: %d - try adjusting FinalTol or NbhdRadius in the input file\n", convergence_failures);
    printf("   Please see the paths marked with '#' in 'main_data' for more information about their accuracy.\n\n");
  }

  if (sharpening_singular + sharpening_failures > 0)
  {
    printf(" Number that failed to sharpen: %d\n", sharpening_singular + sharpening_failures);
    if (sharpening_singular > 0)
      printf("     Singular endpoints: %d\n", sharpening_singular);
    if (sharpening_failures > 0)
    {
      printf("     Sharpening convergence failures: %d - consider using adaptive precision\n", sharpening_failures);
    }
    printf("   Please see the paths marked with '&' in 'main_data' for more information about their accuracy.\n\n");
  }
}

void printSharpeningFailureSummary(int total, int sharpening_failures, int sharpening_singular, int finiteSuccess, int infiniteSuccess, int nonfiniteSuccess)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the failure summary                             *
\***************************************************************/
{
  printf("\nEndpoints Considered for Sharpening: %d\n", total);

  if (finiteSuccess > 0)
    printf(" Finite Sharpening Successes: %d\n", finiteSuccess);

  if (infiniteSuccess > 0)
    printf(" Infinite Sharpening Successes: %d\n", infiniteSuccess);

  if (nonfiniteSuccess > 0)
    printf(" Sharpening Successes: %d\n", nonfiniteSuccess);

  if (sharpening_singular + sharpening_failures > 0)
  {
    printf(" Sharpening Failures: %d\n", sharpening_singular + sharpening_failures);
    if (sharpening_singular > 0)
      printf("     Singular endpoints: %d\n", sharpening_singular);
    if (sharpening_failures > 0)
    {
      printf("     Sharpening convergence failures: %d - try adjusting the number of sharpening digits\n", sharpening_failures);
      printf("       (option 5 on the main menu) or consider using adaptive precision by adding 'MPTYPE: 2;' to the input file.\n");
    }
    printf("   Please see the paths marked with '&' in 'main_data' for more information about their accuracy.\n\n");
  }
}

void getDehomPoint_comp_mp(point_mp dehomPoint, int *origErrorIsInf, double *origErrorEst, comp_mp *sol, int num_vars, preproc_data *PPD, double accuracyEstimate)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: calculates the dehom point                             *
\***************************************************************/
{
  int i, j, largeCoord, currVarGp = 0, var_gp = PPD->num_var_gp, hom_var_gp = PPD->num_hom_var_gp, num_digits = prec_to_digits(dehomPoint->curr_prec);
  int num_gps = var_gp + hom_var_gp, size = num_vars - var_gp;
  double localErrorEst;
  mpf_t epsilon, infNorm, numNorm, denomNorm, accuracy;
  comp_mp recip, largest;

  mpf_init(epsilon);
  mpf_init(infNorm);
  mpf_init(numNorm);
  mpf_init(denomNorm);
  mpf_init(accuracy);
  init_mp(recip);
  init_mp(largest);

  // setup the final size
  change_size_point_mp(dehomPoint, size);
  dehomPoint->size = size;

  // initialize
  *origErrorEst = *origErrorIsInf = 0;
  mpf_set_d(accuracy, accuracyEstimate); 

  if (num_gps == 0)
  { // no variable groups - simply copy the point
    for (i = 0; i < num_vars; i++)
    {
      set_mp(&dehomPoint->coord[i], sol[i]);
    }
    *origErrorEst = accuracyEstimate;
  }
  else
  { // loop over the variable groups
    for (i = 0; i < num_gps; i++)
    { // determine if this group is affine or homogeneous
      if (PPD->type[i] == 0)
      { // homogeneous variable group - normalize by largest coordinate
        largeCoord = -1;
        mpf_set_si(infNorm, -1);
        set_one_mp(largest);
        for (j = 0; j < PPD->size[i]; j++)
        { // compute the modulus of the coordinate
          norm_sqr_mp(epsilon, sol[var_gp + j]);
          if (mpf_cmp(epsilon, infNorm) > 0)
          {
            mpf_set(infNorm, epsilon);
            set_mp(largest, sol[var_gp + j]);
            largeCoord = j;
          }
        }

        // setup the norm of the denominator
        mpf_sqrt(denomNorm, infNorm);

        // verify largest is nonzero - should never be zero!
        if (mpfr_zero_p(largest->r) && mpfr_zero_p(largest->i))
        { // generate a random perturbation
          get_comp_rand_mp(recip);
          mpfr_ui_pow_ui(epsilon, 10, num_digits, __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(recip, recip, epsilon);
          recip_mp(recip, recip);
        }
        else
        { // reciprocate the largest
          recip_mp(recip, largest);
        }

        // normalize the coordinates
        for (j = 0; j < PPD->size[i]; j++)
        {
          if (j == largeCoord)
          { // approximate error, if needed
            if (*origErrorIsInf == 0)
            {
              mpf_set(numNorm, denomNorm);
              localErrorEst = local_error_mp(origErrorIsInf, numNorm, denomNorm, accuracy);

              // update error
              if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
                *origErrorEst = localErrorEst;
            }

            // set to 1
            set_one_mp(&dehomPoint->coord[var_gp - PPD->num_var_gp]);
          }
          else
          { // approximate error, if needed
            if (*origErrorIsInf == 0)
            {
              mpf_abs_mp(numNorm, sol[var_gp]);
              localErrorEst = local_error_mp(origErrorIsInf, numNorm, denomNorm, accuracy);

              // update error
              if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
                *origErrorEst = localErrorEst;
            }

            // multiplity by recip
            mul_mp(&dehomPoint->coord[var_gp - PPD->num_var_gp], sol[var_gp], recip);
          }
          var_gp++;
        }
      }
      else
      { // variable group - dehomogenize
        if (mpfr_zero_p(sol[currVarGp]->r) && mpfr_zero_p(sol[currVarGp]->i))
        { // generate a random perturbation
          get_comp_rand_mp(recip);
          mpfr_ui_pow_ui(epsilon, 10, num_digits, __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(recip, recip, epsilon);
          recip_mp(recip, recip);
        }
        else
        { // reciprocate
          recip_mp(recip, sol[currVarGp]);
        }

        // setup the norm of the denominator
        mpf_abs_mp(denomNorm, sol[currVarGp]); 

        // dehomogenize the coordinates
        for (j = 0; j < PPD->size[i]; j++)
        { // approximate error, if needed
          if (*origErrorIsInf == 0)
          {
            mpf_abs_mp(numNorm, sol[var_gp]);
            localErrorEst = local_error_mp(origErrorIsInf, numNorm, denomNorm, accuracy);

            if (*origErrorIsInf == 0 && localErrorEst > *origErrorEst)
              *origErrorEst = localErrorEst;
          }

          // dehomogenize
          mul_mp(&dehomPoint->coord[var_gp - PPD->num_var_gp], sol[var_gp], recip);

          // update
          var_gp++;
        }

        // increment since we have seen a homogenizing coordinate
        currVarGp++;
      }
    }
  }

  mpf_clear(epsilon);
  mpf_clear(infNorm);
  mpf_clear(numNorm);
  mpf_clear(denomNorm);
  mpf_clear(accuracy);
  clear_mp(recip);
  clear_mp(largest);

  return;
}

void printRefinedPoint_d(FILE *OUT, point_d Pt) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the point in a 'refined' way                    *
\***************************************************************/
{
  int i, size = Pt->size;

  for (i = 0; i < size; i++)
  {
    print_d(OUT, 16, &Pt->coord[i]);
    fprintf(OUT, "\n");
  }

  return;
}

void printRefinedPoint_mp(FILE *OUT, point_mp Pt) 
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the point in a 'refined' way                    *
\***************************************************************/
{
  int i, size = Pt->size;
  
  for (i = 0; i < size; i++)
  { 
    print_mp(OUT, 0, &Pt->coord[i]);
    fprintf(OUT, "\n");
  }

  return;
}

void sort_points(int num_crossings, int *convergence_failures, int *sharpening_failures, int *sharpening_singular, char *inputName, int num_sols, int num_vars, double midpoint_tol, double final_tol, tracker_config_t *T, preproc_data *PPD, int regenToggle, int eqbyeqToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: final_tol is the final tolerance - all relaxation of   *
* this tolerance needs to be done before calling sort_points    *
\***************************************************************/
{
  // initialize counters
  *convergence_failures = *sharpening_failures = *sharpening_singular = 0;

  // check to make sure that we are supposed to create "main_data"
  if (T->outputLevel <= -1)
  { // delete all files except raw_data
    remove_output_files(0, 0, 0);  // track type is 0 - zero dimensional tracking
    remove(inputName);

    // print message to screen about path crossing since this has been done
    if (num_crossings > 0)
    {
      printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame\n\n", num_crossings);
      printf("   Try adjusting TRACKTOLBEFOREEG in the input file\n\n");
    }
  }
  else
  { // create "main_data" 
    int i, j;
    FILE *IN = NULL, *OUT = NULL;
    post_process_t *endPoints = (post_process_t *)bmalloc(num_sols * sizeof(post_process_t));

    // open 'raw_data'
    IN = fopen("raw_data", "r");
    if (IN == NULL)
    {
      printf("ERROR: The file 'raw_data' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }
    // move past the number of variables and the dimension
    fscanf(IN, "%d\n%d\n", &i, &j); 

    // open 'main_data'
    OUT = fopen("main_data", "w");

    // read in the data
    for (i = 0; i < num_sols; i++)
    { // read in the next solution and relevant data
      setupPostProcess(&j, IN, &endPoints[i], num_vars, T->MPType);

      // check success 
      if (endPoints[i].success == retVal_sharpening_singular_endpoint)
        *sharpening_singular = *sharpening_singular + 1;
      else if (endPoints[i].success == retVal_sharpening_failed)
        *sharpening_failures = *sharpening_failures + 1;
      else if (endPoints[i].success != 1)
        *convergence_failures = *convergence_failures + 1;
    }
    // close 'raw_data'
    fclose(IN);

    // do the post-processing
    zeroDimPostProcess(OUT, endPoints, num_sols, num_vars, final_tol, T, PPD, num_crossings, *convergence_failures, inputName, regenToggle, eqbyeqToggle);

    // close 'main_data' 
    fclose(OUT);

    // clear memory
    for (i = num_sols - 1; i >= 0; i--)
    {
      if (endPoints[i].sol_prec >= 64)
      { // clear _mp
        mpf_clear(endPoints[i].function_resid_mp);
        mpf_clear(endPoints[i].newton_resid_mp);
        for (j = 0; j < num_vars; j++)
        {
          clear_mp(endPoints[i].sol_mp[j]);
        }
      }
      free(endPoints[i].sol_mp);
      free(endPoints[i].sol_d);
    }
    free(endPoints);
  }

  return;
}

void zeroDimPostProcess(FILE *OUT, post_process_t *endPoints, int num_sols, int num_vars, double final_tol, tracker_config_t *T, preproc_data *PPD, int num_crossings, int convergence_failures, char *inputName, int regenToggle, int eqbyeqToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the actual post processing for a zero dim run     *
\***************************************************************/
{
  int i, num_good_sols = 0, size = 0;
  FILE *NAMES = NULL;
  char ch, **name_table = (char **)bmalloc(num_vars * sizeof(char *));
  int *origErrorIsInf = (int *)bmalloc(num_sols * sizeof(int)); 
  double *origErrorEst = (double *)bmalloc(num_sols * sizeof(double));
  point_d *dehomPoints_d = T->MPType == 1 ? NULL : (point_d *)bmalloc(num_sols * sizeof(point_d));
  point_mp *dehomPoints_mp = T->MPType == 0 ? NULL : (point_mp *)bmalloc(num_sols * sizeof(point_mp));

  // Reading in the names of the variables.
  NAMES = fopen("names.out", "r");
  if (NAMES == NULL)
  {
    printf("ERROR: 'names.out' does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  for (i = 0; i < num_vars; i++)
  { // initial allocation
    size = 1;
    name_table[i] = (char *)bmalloc(size * sizeof(char));
    // read in name
    while ((name_table[i][size - 1] = fgetc(NAMES)) != '\n')
    {
      size++;
      name_table[i] = (char *)brealloc(name_table[i], size * sizeof(char));
    }
    name_table[i][size - 1] = '\0';
  }
  fclose(NAMES);

  // setup the dehomogenized points and error estimates
  for (i = 0; i < num_sols; i++)
  { // setup dehom
    if (endPoints[i].sol_prec < 64)
    { // setup _d
      init_point_d(dehomPoints_d[i], 0);
      getDehomPoint_comp_d(dehomPoints_d[i], &origErrorIsInf[i], &origErrorEst[i], endPoints[i].sol_d, num_vars, PPD, endPoints[i].accuracy_estimate);
    }
    else
    { // setup _mp
      initMP(endPoints[i].sol_prec);
      init_point_mp(dehomPoints_mp[i], 0);
      getDehomPoint_comp_mp(dehomPoints_mp[i], &origErrorIsInf[i], &origErrorEst[i], endPoints[i].sol_mp, num_vars, PPD, endPoints[i].accuracy_estimate);
    }
  }

  // print top of main_data
  printFileHeader(OUT, num_vars, PPD, name_table);

  // create various output from this data
  zeroDimOutputChart(endPoints, dehomPoints_d, dehomPoints_mp, name_table, num_sols, num_vars, PPD, T->finiteThreshold, T->real_threshold, final_tol, T->cond_num_threshold, regenToggle, eqbyeqToggle);
  createRawSoln(endPoints, dehomPoints_d, dehomPoints_mp, num_sols, num_vars, PPD);

  // Now we sort the points - and write the body of main_data
  if (regenToggle)
    fprintf(OUT, "\nNOTE: Since regeneration is being used, only non-singular solutions are printed.\n\n");
  else if (eqbyeqToggle)
    fprintf(OUT, "\nNOTE: Since equation-by-equation is being used, only non-singular solutions are printed.\n\n");

  for (i = 0; i < num_sols; i++)
    if (endPoints[i].multiplicity > 0 && (!(regenToggle || eqbyeqToggle) || ((regenToggle || eqbyeqToggle) && endPoints[i].isSing == 0)))  // print only for regen/eqbyeq if non-singular
    {
      if (endPoints[i].sol_prec < 64)
      { // print header for the solution
        printMainDataPointHeader_d(OUT, endPoints[i].sol_num, endPoints[i].path_num, endPoints[i].cond_est, endPoints[i].function_resid_d, endPoints[i].newton_resid_d, endPoints[i].final_t, endPoints[i].accuracy_estimate, endPoints[i].sol_prec, endPoints[i].first_increase, endPoints[i].cycle_num, endPoints[i].success, origErrorIsInf[i], origErrorEst[i]);
        // print the point
        printRefinedPoint_d(OUT, dehomPoints_d[i]);
      }
      else
      { // print header for the solution
        printMainDataPointHeader_mp(OUT, endPoints[i].sol_num, endPoints[i].path_num, endPoints[i].cond_est, endPoints[i].function_resid_mp, endPoints[i].newton_resid_mp, endPoints[i].final_t, endPoints[i].accuracy_estimate, endPoints[i].sol_prec, endPoints[i].first_increase, endPoints[i].cycle_num, endPoints[i].success, origErrorIsInf[i], origErrorEst[i]);
        // print the point
        printRefinedPoint_mp(OUT, dehomPoints_mp[i]); 
      }

      // print footer for the solution
      printMainDataPointFooter(OUT, i, num_sols, endPoints);

      num_good_sols++;
    }

  fprintf(OUT, "-------------------------\n");

  if (num_good_sols == 1)
    fprintf(OUT, "At tol=%.12e, there appears to be a unique solution.\n", final_tol);
  else
    fprintf(OUT, "At tol=%.12e, there appear to be %d solutions.\n", final_tol, num_good_sols);

  if (num_crossings > 0)
  {
    printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame\n\n", num_crossings);
    fprintf(OUT, "\nIt appears that %d path crossing(s) occurred prior to t=tEndgame\n\n", num_crossings);
    printf("   Try adjusting TRACKTOLBEFOREEG in the input file\n\n");
  }

  if (convergence_failures == 1)
  {
    fprintf(OUT, "\nThere is 1 path that was tracked near the target time which did not converge. Try adjusting FinalTol or NbhdRadius in the input file.\n");
    fprintf(OUT, "This path is marked with '#' after its path number.\n");
  }
  else if (convergence_failures > 1)
  {
    fprintf(OUT, "\nThere are %d paths that were tracked near the target time which did not converge. Try adjusting FinalTol or NbhdRadius in the input file.\n", convergence_failures);
    fprintf(OUT, "These paths are marked with '#' after the path numbers.\n");
  }

  // print the input onto the bottom of refined solutions
  NAMES = fopen(inputName, "r");
  if (NAMES == NULL)
  {
    printf("ERROR: '%s' does not exist!\n", inputName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  while((ch = fgetc(NAMES)) != EOF)
    fputc(ch, OUT);
  fclose(NAMES);
  remove(inputName);

  // clear memory
  for (i = num_vars - 1; i >= 0; i--)
    free(name_table[i]);
  free(name_table);
  for (i = 0; i < num_sols; i++)
  {
    if (endPoints[i].sol_prec < 64)
    {
      clear_point_d(dehomPoints_d[i]);
    }
    else
    {
      clear_point_mp(dehomPoints_mp[i]);
    }
  }
  if (T->MPType != 1)
    free(dehomPoints_d);
  if (T->MPType != 0)
    free(dehomPoints_mp);
  free(origErrorEst);
  free(origErrorIsInf);

  return;
}

void zeroDimOutputChart(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, char **name_table, int num_sols, int num_vars, preproc_data *PPD, double maxNorm, double realTol, double finalTol, double maxCondNum, int regenToggle, int eqbyeqToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the output charts                               *
\***************************************************************/
{
  int i;

  // find solution numbers & multiplicity
  findMultSol(endPoints, num_sols, num_vars, PPD, finalTol);

  if (PPD->num_var_gp > 0)
  { // find the finite solutions
    findFiniteSol(endPoints, dehomPoints_d, dehomPoints_mp, num_sols, num_vars, PPD, maxNorm);
  }
  else
  { // set finite as -1
    for (i = 0; i < num_sols; i++)
      endPoints[i].isFinite = -1;
  }

  // determine which ones are real
  findRealSol(endPoints, dehomPoints_d, dehomPoints_mp, num_sols, num_vars, PPD, realTol);

  // determine which ones are singular
  findSingSol(endPoints, dehomPoints_d, dehomPoints_mp, num_sols, num_vars, PPD, maxCondNum, finalTol, regenToggle);

  if (PPD->num_var_gp > 0)
  { // print the finite summary
    singularPointSummary(num_sols, endPoints, 1, regenToggle, eqbyeqToggle);
    multiplicityPointSummary(num_sols, endPoints, 1, regenToggle, eqbyeqToggle);
    // print the infinite summary
    singularPointSummary(num_sols, endPoints, 0, regenToggle, eqbyeqToggle);
    multiplicityPointSummary(num_sols, endPoints, 0, regenToggle, eqbyeqToggle);
  }
  else
  { // homogenized or user defined homotopy - we don't know what infinity is so we just put everything together!
    singularPointSummary(num_sols, endPoints, -1, regenToggle, eqbyeqToggle);
    multiplicityPointSummary(num_sols, endPoints, -1, regenToggle, eqbyeqToggle);
  }

  return;
}

void findMultSol(post_process_t *endPoints, int num_sols, int num_vars, preproc_data *PPD, double finalTol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the multiplicities & solution numbers            *
\***************************************************************/
{
  int i, j, k, curr_sol = 0, max_prec = 52, indexI = 0, indexJ = 0, cont = 1;
  double tempTot_d;

  // initialize multiplicity & solution and find the maximum precision used
  for (i = 0; i < num_sols; i++)
  { // initialize multiplicity & solution number
    endPoints[i].multiplicity = 1;
    endPoints[i].sol_num = -1;
    
    if (max_prec < endPoints[i].sol_prec)
      max_prec = endPoints[i].sol_prec;
  }

  if (max_prec < 64)
  { // sort using double precision
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(num_sols * sizeof(sortStruct_d));
    comp_d tempComp;

    for (i = 0; i < num_sols; i++)
    { // initialize sortPts
      sortPts[i].path_num = i;
      // find norm using double prec
      sortPts[i].norm = 0;
      for (j = 0; j < num_vars; j++)
      {
        sortPts[i].norm += norm_sqr_d(endPoints[i].sol_d[j]);
      }
      sortPts[i].norm = sqrt(sortPts[i].norm);
    }
 
    // sort
    qsort(sortPts, num_sols, sizeof(sortStruct_d), sort_order_d);
  
    // loop through to classify the paths
    for (i = 0; i < num_sols; i++)
    {
      indexI = sortPts[i].path_num;
      if (endPoints[indexI].sol_num == -1)
      { // set the solution number correctly
        endPoints[indexI].sol_num = curr_sol;
        curr_sol++;
 
        // loop through to find the multiplicity
        cont = 1;
        j = i;
        do 
        { // increment the counter (start at i + 1)
          j++;
  
          if (j < num_sols)
          { // check norm diff
            indexJ = sortPts[j].path_num;
            if (sortPts[j].norm - sortPts[i].norm > finalTol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we can continue & this one has not yet been classified
          if (cont && endPoints[indexJ].sol_num == -1)
          {
            tempTot_d = 0;
            for (k = 0; k < num_vars; k++)
            { // find the difference
              sub_d(tempComp, endPoints[indexJ].sol_d[k], endPoints[indexI].sol_d[k]);
              tempTot_d += tempComp->r * tempComp->r + tempComp->i * tempComp->i;
            }
            tempTot_d = sqrt(tempTot_d);
            if (tempTot_d < finalTol)
            { // lump solution j in with solution
              endPoints[indexJ].multiplicity = -1;
              endPoints[indexI].multiplicity++;
              endPoints[indexJ].sol_num = endPoints[indexI].sol_num;
            }
          }
        } while (cont);
      }
    }
    // clear memory
    free(sortPts);
  }
  else
  { // sort using MP
    sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(num_sols * sizeof(sortStruct_mp));
    comp_mp tempComp1, tempComp2;
    mpf_t tempMPF, tempTot_mp;

    init_mp2(tempComp1, max_prec);
    init_mp2(tempComp2, max_prec);
    mpf_init2(tempMPF, max_prec);
    mpf_init2(tempTot_mp, max_prec);

    for (i = 0; i < num_sols; i++)
    { // initialize sortPts
      sortPts[i].path_num = i;
      mpf_init2(sortPts[i].norm, max_prec);

      if (endPoints[i].sol_prec < 64)
      { // find norm using double prec     
        tempTot_d = 0;
        for (j = 0; j < num_vars; j++)
        {
          tempTot_d += norm_sqr_d(endPoints[i].sol_d[j]);
        }
        tempTot_d = sqrt(tempTot_d);
        mpf_set_d(sortPts[i].norm, tempTot_d);
      }
      else
      { // find norm using MP
        mpf_set_ui(sortPts[i].norm, 0);
        for (j = 0; j < num_vars; j++)
        {
          mpf_mul(tempMPF, endPoints[i].sol_mp[j]->r, endPoints[i].sol_mp[j]->r);
          mpf_add(sortPts[i].norm, sortPts[i].norm, tempMPF);
          mpf_mul(tempMPF, endPoints[i].sol_mp[j]->i, endPoints[i].sol_mp[j]->i);
          mpf_add(sortPts[i].norm, sortPts[i].norm, tempMPF);
        }
        mpf_sqrt(sortPts[i].norm, sortPts[i].norm);
      }
    }

    // sort
    qsort(sortPts, num_sols, sizeof(sortStruct_mp), sort_order_mp);

    // loop through to classify the paths
    for (i = 0; i < num_sols; i++)
    {
      indexI = sortPts[i].path_num;
      if (endPoints[indexI].sol_num == -1)
      { // set the solution number correctly
        endPoints[indexI].sol_num = curr_sol;
        curr_sol++;

        // loop through to find the multiplicity
        cont = 1;
        j = i;
        do 
        { // increment the counter (start at i + 1)
          j++; 

          if (j < num_sols)
          { // check norm diff
            indexJ = sortPts[j].path_num;
            mpf_sub(tempMPF, sortPts[j].norm, sortPts[i].norm);
            if (mpf_get_d(tempMPF) > finalTol)
              cont = 0;
          }
          else
            cont = 0;
 
          // check to see if we can continue & this one has not yet been classified
          if (cont && endPoints[indexJ].sol_num == -1)
          {
            mpf_set_ui(tempTot_mp, 0);
            for (k = 0; k < num_vars; k++)
            { // find the difference
              if (endPoints[indexI].sol_prec < 64)
              {
                d_to_mp(tempComp1, endPoints[indexI].sol_d[k]);
              }
              else
              {
                set_mp(tempComp1, endPoints[indexI].sol_mp[k]);
              }

              if (endPoints[indexJ].sol_prec < 64)
              {
                d_to_mp(tempComp2, endPoints[indexJ].sol_d[k]);
              }
              else
              {
                set_mp(tempComp2, endPoints[indexJ].sol_mp[k]);
              }

              sub_mp(tempComp1, tempComp1, tempComp2);
              mpf_mul(tempMPF, tempComp1->r, tempComp1->r);
              mpf_add(tempTot_mp, tempTot_mp, tempMPF);
              mpf_mul(tempMPF, tempComp1->i, tempComp1->i);
              mpf_add(tempTot_mp, tempTot_mp, tempMPF);
            }
            mpf_sqrt(tempTot_mp, tempTot_mp);

            if (mpf_get_d(tempTot_mp) < finalTol)
            { // lump solution j in with solution
              endPoints[indexJ].multiplicity = -1;
              endPoints[indexI].multiplicity++;
              endPoints[indexJ].sol_num = endPoints[indexI].sol_num;
            }
          }
        } while (cont);
      }
    }

    // clear memory
    for (i = num_sols - 1; i >= 0; i--)
    {
      mpf_clear(sortPts[i].norm);
    }
    free(sortPts);
    clear_mp(tempComp1);
    clear_mp(tempComp2);
    mpf_clear(tempMPF);
    mpf_clear(tempTot_mp);
  }

  return;
}

void findRealSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double realTol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the real solutions                               *
\***************************************************************/
{
  int i, j, k, l, real_count = 0;
  FILE *REAL = NULL;

  if (PPD->num_var_gp > 0)
    REAL = fopen("real_finite_solutions", "w");
  else
    REAL = fopen("real_solutions", "w");

  // leave room for the number of solutions
  fprintf(REAL, "                                                 \n");

  // initialize
  for (i = 0; i < num_sols; i++)
    endPoints[i].isReal = 0;

  for (i = 0; i < num_sols; i++)
   if (endPoints[i].multiplicity > 0)
    { // check for real solution
      if (endPoints[i].sol_prec < 64)
      {
        endPoints[i].isReal = checkForReal_d(dehomPoints_d[i], realTol);
      }
      else
      {
        endPoints[i].isReal = checkForReal_mp(dehomPoints_mp[i], realTol);
      }

      if (endPoints[i].isReal && (PPD->num_var_gp == 0 || endPoints[i].isFinite))
      { // print info to REAL
  
        // increment the number of singular solutions printed
        real_count++;

        // move to the next line
        fprintf(REAL, "\n");

        if (endPoints[i].sol_prec < 64)
        { // print information about point to REAL in double precision
          for (j = 0; j < dehomPoints_d[i]->size; j++)
          {
            print_d(REAL, 16, &dehomPoints_d[i]->coord[j]);
            fprintf(REAL, "\n");
          }
        } 
        else // prec >= 64
        { // print information about point to REAL in multi precision
          for (j = 0; j < dehomPoints_mp[i]->size; j++)
          {
            print_mp(REAL, 0, &dehomPoints_mp[i]->coord[j]);
            fprintf(REAL, "\n");
          }
        }

        l = 0;
        // print the other points with the same solution number
        for (k = 1; k < endPoints[i].multiplicity; k++)
        { // increment the number of real solutions printed
          real_count++;
          
          // move to the next line
          fprintf(REAL, "\n");

          // find the point with the same solution number
          while (l == i || endPoints[l].sol_num != endPoints[i].sol_num)
            l++;

          if (endPoints[l].sol_prec < 64)
          { // print information about point to REAL in double precision
            for (j = 0; j < dehomPoints_d[l]->size; j++)
            {
              print_d(REAL, 16, &dehomPoints_d[l]->coord[j]);
              fprintf(REAL, "\n");
            }
          }
          else // prec >= 64
          { // print information about point to REAL in multi precision
            for (j = 0; j < dehomPoints_mp[l]->size; j++)
            {
              print_mp(REAL, 0, &dehomPoints_mp[l]->coord[j]);
              fprintf(REAL, "\n");
            }
          }
          l++;
        }
      }
    }

  // print an extra line at the bottom
  fprintf(REAL, "\n");

  // rewind to the beginning and print the number of solutions printed
  rewind(REAL);
  fprintf(REAL, "%d", real_count);

  // close the files
  fclose(REAL);

  return;
}

int checkForReal_d(point_d Pt, double realTol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines is Pt is real or not                        *
\***************************************************************/
{
  int j, size = Pt->size, isReal = 1;

  for (j = 0; j < size && isReal; j++)
    if (fabs(Pt->coord[j].i) > realTol)
      isReal = 0;

  return isReal;
}

int checkForReal_mp(point_mp Pt, double realTol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines is Pt is real or not                        *
\***************************************************************/
{
  int j, isReal = 1, size = Pt->size;

  for (j = 0; j < size && isReal; j++)
    if (fabs(mpf_get_d(Pt->coord[j].i)) > realTol)
      isReal = 0;

  return isReal;
}

void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the nonsingular/singular sols                   *
\***************************************************************/
{
  int i, j, k, l, sing_count = 0, nonsing_count = 0;
  FILE *NONSING = fopen("nonsingular_solutions", "w"), *SINGU = fopen("singular_solutions", "w");

  // leave room for the number of solutions
  fprintf(SINGU, "                                                 \n");
  fprintf(NONSING, "                                                 \n");

  for (i = 0; i < num_sols; i++)
    if (endPoints[i].multiplicity > 0)
    { // regeneration endpoints are always non-singular
      if (!regenToggle && ((endPoints[i].cond_est > maxCondNum) || (endPoints[i].cond_est < 0.0) || (endPoints[i].multiplicity > 1)))
        endPoints[i].isSing = 1;
      else
        endPoints[i].isSing = 0;

      if (endPoints[i].isSing)
      { // print info to SINGU

        // increment the number of singular solutions printed
        sing_count++;

        // move to the next line
        fprintf(SINGU, "\n");

        // print the point
        if (endPoints[i].sol_prec < 64)
        { // print information about point to SINGU in double precision in the user defined variables
          for (j = 0; j < dehomPoints_d[i]->size; j++)
          {
            print_d(SINGU, 16, &dehomPoints_d[i]->coord[j]);
            fprintf(SINGU, "\n");
          }
        }
        else // prec >= 64
        { // print information about point to SINGU in multi precision
          for (j = 0; j < dehomPoints_mp[i]->size; j++)
          {
            print_mp(SINGU, 0, &dehomPoints_mp[i]->coord[j]);
            fprintf(SINGU, "\n");
          }
        }

        l = 0;
        // print the other points with the same solution number
        for (k = 1; k < endPoints[i].multiplicity; k++)
        { // increment the number of singular solutions printed
          sing_count++;

          // move to the next line
          fprintf(SINGU, "\n");

          // find the point with the same solution number
          while (l == i || endPoints[l].sol_num != endPoints[i].sol_num)
            l++;

          if (endPoints[l].sol_prec < 64)
          { // print information about point to SINGU in double precision
            for (j = 0; j < dehomPoints_d[l]->size; j++)
            {
              print_d(SINGU, 16, &dehomPoints_d[l]->coord[j]);
              fprintf(SINGU, "\n");
            }
          }
          else // prec >= 64
          { // print information about point to SINGU in multi precision
            for (j = 0; j < dehomPoints_mp[l]->size; j++)
            {
              print_mp(SINGU, 0, &dehomPoints_mp[l]->coord[j]);
              fprintf(SINGU, "\n");
            }
          }
          l++;
        }
      }
      else
      { // print info to NONSING

        // increment the number of nonsingular solutions printed
        nonsing_count++;

        // move to the next line
        fprintf(NONSING, "\n");

        if (endPoints[i].sol_prec < 64)
        { // print information about point to NONSING in double precision in the user-defined variables
          for (j = 0; j < dehomPoints_d[i]->size; j++)
          {
            print_d(NONSING, 16, &dehomPoints_d[i]->coord[j]);
            fprintf(NONSING, "\n");
          }
        }
        else // prec >= 64
        { // print information about point to NONSING in multi precision
          for (j = 0; j < dehomPoints_mp[i]->size; j++)
          {
            print_mp(NONSING, 0, &dehomPoints_mp[i]->coord[j]);
            fprintf(NONSING, "\n");
          }
        }

        // there are no other points with the same solution number since this point is NONSINGULAR!
      }
    }

  // print an extra line at the bottom
  fprintf(SINGU, "\n");
  fprintf(NONSING, "\n");

  // rewind to the beginning and print the number of solutions printed
  rewind(SINGU);
  rewind(NONSING);
  fprintf(SINGU, "%d", sing_count);
  fprintf(NONSING, "%d", nonsing_count);

  // close the files
  fclose(SINGU);
  fclose(NONSING);

  return;
}

void findFiniteSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxNorm)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the finite sols                                 *
\***************************************************************/
{
  int i, j, k, l, finite_count = 0;
  FILE *FINITE = fopen("finite_solutions", "w");

  // leave room for the number of solutions
  fprintf(FINITE, "                                                 \n");

  for (i = 0; i < num_sols; i++)
    if (endPoints[i].multiplicity > 0)
    { // find the dehom point
      if (endPoints[i].sol_prec < 64)
      { // use double precision
        if (infNormVec_d(dehomPoints_d[i]) < maxNorm)
        { // print information about point to FINITE
          endPoints[i].isFinite = 1;

          // increment the number of finite solutions printed
          finite_count++;

          // move to the next line
          fprintf(FINITE, "\n");

          // print the point
          for (j = 0; j < dehomPoints_d[i]->size; j++)
          {
            print_d(FINITE, 16, &dehomPoints_d[i]->coord[j]);
            fprintf(FINITE, "\n");
          }
        }
        else
        {
          endPoints[i].isFinite = 0;
        }
      }
      else
      { // use higher precision
        if (infNormVec_mp(dehomPoints_mp[i]) < maxNorm )
        { // print information about point to FINITE
          endPoints[i].isFinite = 1;

          // increment the number of finite solutions printed
          finite_count++;

          // move to the next line
          fprintf(FINITE, "\n");

          // print the point
          for (j = 0; j < dehomPoints_mp[i]->size; j++)
          {
            print_mp(FINITE, 0, &dehomPoints_mp[i]->coord[j]);
            fprintf(FINITE, "\n");
          }
        }
        else
        {
          endPoints[i].isFinite = 0;
        }
      }

      if (endPoints[i].isFinite)
      { // print the other points with the same solution number

        l = 0;
        for (k = 1; k < endPoints[i].multiplicity; k++)
        { // increment the number of singular solutions printed
          finite_count++;

          // move to the next line
          fprintf(FINITE, "\n");

          // find the point with the same solution number
          while (l == i || endPoints[l].sol_num != endPoints[i].sol_num)
            l++;

          // setup finite
          endPoints[l].isFinite = 1;

          if (endPoints[l].sol_prec < 64)
          { // print information about point to FINITE in double precision
            for (j = 0; j < dehomPoints_d[l]->size; j++)
            {
              print_d(FINITE, 16, &dehomPoints_d[l]->coord[j]);
              fprintf(FINITE, "\n");
            }
          }
          else // Sols_prec[l] >= 64
          { // print information about point to FINITE in multi precision
            for (j = 0; j < dehomPoints_mp[l]->size; j++)
            {
              print_mp(FINITE, 0, &dehomPoints_mp[l]->coord[j]);
              fprintf(FINITE, "\n");
            }
          }
          l++;
        }
      }
    }

  // print an extra line at the bottom
  fprintf(FINITE, "\n");
  
  // rewind to the beginning and print the number of solutions printed
  rewind(FINITE);
  fprintf(FINITE, "%d", finite_count);

  // close the files
  fclose(FINITE);

  return;
}

void createRawSoln(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates raw_solutions                                  *
\***************************************************************/
{
  int i, j;
  FILE *OUT = fopen("raw_solutions", "w");

  // print the number of points at the top
  fprintf(OUT, "%d\n\n", num_sols);

  for (i = 0; i < num_sols; i++)
  { // path number
    fprintf(OUT, "%d\n", endPoints[i].path_num);

    // print the dehom point 
    if (endPoints[i].sol_prec < 64)
    { // use double precision
      for (j = 0; j < dehomPoints_d[i]->size; j++)
      {
        print_d(OUT, 16, &dehomPoints_d[i]->coord[j]);
        fprintf(OUT, "\n");
      }
      fprintf(OUT, "\n");
    }
    else
    { // use higher precision
      for (j = 0; j < dehomPoints_mp[i]->size; j++)
      {
        print_mp(OUT, 0, &dehomPoints_mp[i]->coord[j]);
        fprintf(OUT, "\n");
      }
      fprintf(OUT, "\n");
    }
  }

  // close the file
  fclose(OUT);

  return;
}

int setupPostProcess(int *orig_prec, FILE *IN, post_process_t *endPoint, int size, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 if successful & -1 otherwise                 *
* NOTES: reads in the data from IN to setup endPoint            *
\***************************************************************/
{
  int i, retVal = 0;

  // read in the path number
  endPoint->path_num = -1;

  fscanf(IN, "%d\n", &endPoint->path_num);

  if (endPoint->path_num == -1)
  { // we are the bottom so we return an error
    retVal = -1;
  }
  else
  { // continue reading in everything else
    fscanf(IN, "%d\n", orig_prec);

    if (MPType == 0 || (*orig_prec < 64 && MPType == 2))
    { // answer in double precision
      endPoint->sol_prec = 52;
      endPoint->size_sol = size;
      endPoint->sol_d  = (comp_d *)bmalloc(size * sizeof(comp_d));
      endPoint->sol_mp = NULL;

      for (i = 0; i < size; i++)
        fscanf(IN, "%lf %lf\n", &endPoint->sol_d[i]->r, &endPoint->sol_d[i]->i);

      fscanf(IN, "%lf\n", &endPoint->function_resid_d);
      fscanf(IN, "%lf\n", &endPoint->cond_est);
      fscanf(IN, "%lf\n", &endPoint->newton_resid_d);
    }
    else // MPType == 1 || (*orig_prec >= 64 && MPType == 2)
    { // answer in multi precision
      endPoint->sol_prec = MPType == 2 ? *orig_prec : mpf_get_default_prec();
      endPoint->size_sol = size;
      endPoint->sol_d  = NULL;
      endPoint->sol_mp = (comp_mp *)bmalloc(size * sizeof(comp_mp));

      for (i = 0; i < size; i++)
      { // initialize to the correct precision
        init_mp2(endPoint->sol_mp[i], endPoint->sol_prec);
        mpf_inp_str(endPoint->sol_mp[i]->r, IN, 10);
        mpf_inp_str(endPoint->sol_mp[i]->i, IN, 10);
        fscanf(IN, "\n");
      }

      mpf_init2(endPoint->function_resid_mp, endPoint->sol_prec);
      mpf_inp_str(endPoint->function_resid_mp, IN, 10);
      fscanf(IN, "\n");

      fscanf(IN, "%lf\n", &endPoint->cond_est);

      mpf_init2(endPoint->newton_resid_mp, endPoint->sol_prec);
      mpf_inp_str(endPoint->newton_resid_mp, IN, 10);
      fscanf(IN, "\n");
    }

    fscanf(IN, "%lf\n", &endPoint->final_t);
    fscanf(IN, "%lf\n", &endPoint->accuracy_estimate);
    fscanf(IN, "%lf\n", &endPoint->first_increase);
    fscanf(IN, "%d\n",  &endPoint->cycle_num);
    fscanf(IN, "%d\n",  &endPoint->success);

    // initialize the other parts
    endPoint->sol_num = 0;
    endPoint->multiplicity = 1;
    endPoint->isFinite = 0;

    retVal = 0;
  }

  return retVal;
}

