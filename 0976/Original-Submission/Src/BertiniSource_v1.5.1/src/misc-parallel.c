// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"

// temporary mpf_t variables using in basic multiprecision arithmetic operations
// they provide a nice speed up since they do not have to be initialized & cleared all the time
// so they are pointers so that each thread can have its own set of them to use
mpf_t  *_tempMPF1, *_tempMPF2, *_tempMPF3;
int mp_need_init = 1;

/******************************************************/
void init_trackingStats(trackingStats *S)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes trackingStats                              *
\***************************************************************/
{
  S->numPoints = S->successes = S->failures = S->junkCount = S->nanCount = S->infCount = S->securityCount = 
      S->sizeCount = S->PSEGCount = S->precCount = S->cycleCount = S->stepCount = S->refineCount = S->otherCount = 0;

  return;
}

/******************************************************/
void add_trackingStats(trackingStats *Total, trackingStats *Input, int num_input)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: adds trackingStats together                            *
\***************************************************************/
{
  int i;

  // initialize Total to all 0
  init_trackingStats(Total);

  // add them up
  for (i = 0; i < num_input; i++)
  {
    Total->successes += Input[i].successes;
    Total->failures += Input[i].failures;
    Total->junkCount += Input[i].junkCount;
    Total->nanCount += Input[i].nanCount;
    Total->infCount += Input[i].infCount;
    Total->securityCount += Input[i].securityCount;
    Total->sizeCount += Input[i].sizeCount;
    Total->PSEGCount += Input[i].PSEGCount;
    Total->precCount += Input[i].precCount;
    Total->cycleCount += Input[i].cycleCount;
    Total->stepCount += Input[i].stepCount;
    Total->refineCount += Input[i].refineCount;
    Total->otherCount += Input[i].otherCount;
  }
  
  // the total number of points is success + failures
  Total->numPoints = Total->successes + Total->failures;

  return;
}

/******************************************************/
void reproduceInputFile(char *outputName, char *inName, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates the file needed to reproduce this run          *
\***************************************************************/
{
  char ch;
  FILE *OUT = fopen(outputName, "w");
  FILE *FUNC = fopen(inName, "r");
  if (FUNC == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n", inName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  fprintf(OUT, "\n\n*************** input file needed to reproduce this run ***************\n\n");

  // print the configurations
  printConfigValues(OUT, T, trackType, genType, randomSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // print the function
  fprintf(OUT, "\nINPUT\n\n");
  ch = fgetc(FUNC);
  while (ch != EOF)
  {
    fprintf(OUT, "%c", ch);
    ch = fgetc(FUNC);
  }
  fprintf(OUT, "\n");

  printVersion(OUT);

  // close the files
  fclose(OUT);
  fclose(FUNC);

  return;
}

/******************************************************/
double maxDiffMat_d(mat_d A, mat_d B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the maximum modulus of the entry in A - B      *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, m = A->rows, n = A->cols;
  double tempD, retVal = 0;
  comp_d tempC;

  // do error checking
  if (A->rows != B->rows || A->cols != B->cols)
  {
    printf("ERROR: Attempting to subtract matrices with dimensions that do not match!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      sub_d(tempC, &A->entry[i][j], &B->entry[i][j]);
      tempD = d_abs_d(tempC);
      if (tempD > retVal)
        retVal = tempD;
    }

  return retVal;
}

/******************************************************/
double infNormVec_d(vec_d X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the vector                *
* NOTES:                                                        *
\***************************************************************/
// (JDH - 6/21/06) rewritten to avoid evaluating d_abs_d evaluate twice
{
  int j; 
  double temp, max = 0;

  for (j = 0; j < X->size; j++)
  {
    temp = norm_sqr_d(&X->coord[j]);
    if (max < temp)
      max = temp;
  }
  return sqrt(max);
}

/******************************************************/
double infNormMatRow_d(mat_d A, int i)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the ith row of the matrix *
* NOTES:                                                        *
\***************************************************************/
// (JDH - 6/21/06) - rewritten to avoid evaluating d_abs_d twice
{
  int j; 
  double temp, max = 0;

  for (j = 0; j < A->cols; j++)
  {
    temp = d_abs_d(&A->entry[i][j]);
    if (max < temp)
      max = temp;
  }
  return max;
}

/******************************************************/
double infNormMat_d(mat_d M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the matrix                *
* NOTES:                                                        *
\***************************************************************/
// (JDH - 6/21/06) - rewritten to avoid repeatedly setting things to 0
{
  int i, j, numRows = M->rows, numCols = M->cols;
  double norm, max = 0.0;

  for (i = 0; i < numRows; i++)
  {
    norm = 0.0; 
    for (j = 0; j < numCols; j++)
      norm += norm_sqr_d(&M->entry[i][j]);

    if (max < norm)
      max = norm;
  }
  max = sqrt(max);

  return max;
}

double frobNormMat_d(mat_d M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the Frobenius norm of the matrix               *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, numRows = M->rows, numCols = M->cols;
  double norm = 0;

  for (i = 0; i < numRows; i++)
    for (j = 0; j < numCols; j++)
      norm += norm_sqr_d(&M->entry[i][j]);

  return sqrt(norm);
}

/******************************************************/
double conditionNumber_d(mat_d M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: condition number of matrix                     *
* NOTES:                                                        *
\***************************************************************/
{
  if ((M->rows == 1) && (M->cols == 1))
  {
    return 1.0;
  }
  else
  { // to avoid having 2 copies of the same function
    double tempD1, tempD2;
    return conditionNumber_and_norms_d(M, &tempD1, &tempD2);
  }
}

/******************************************************/
double conditionNumber_and_norms_d(mat_d M, double *norm_M, double *norm_M_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: condition number of matrix with norms of it &  *
* its inverse                                                   *
* NOTES: uses SVD                                               *
\***************************************************************/
{
  double condNum, tol = 1e-14;

  // find the infinity norm of the matrix
  *norm_M = infNormMat_d(M);

  // find the minimum singular value
  min_svd_d(norm_M_inv, M, tol);

  if (*norm_M_inv > 0)
  {
    *norm_M_inv = 1.0 / *norm_M_inv;
  }
  else
  { // report norm of the inverse as 1 / epsilon
    *norm_M_inv = 1e16;
  }

  // calculate CN
  condNum = *norm_M * (*norm_M_inv);

  if (condNum < 1)
    condNum = 1;

  return condNum;
}

/*****************************************************/
void mat_d_to_mp_rat(mat_mp MP, mpq_t ***RAT, mat_d D, int digits, int prec, int need_to_init_mp, int need_to_init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts a _d to _mp & mpq_t using digits & prec       *
\***************************************************************/
{
  int i, j;

  // check to see if MP needs initialized
  if (need_to_init_mp)
  { // this will initialize all of the matrix MP, while RAT is a specific size so we can initialize 
    // each entry of it in the loop
    init_mat_mp2(MP, D->rows, D->cols, prec);
    need_to_init_mp = 0; // since we have initialized all of MP here, we do not need to in the loop
  }
  else
  { // make sure MP is large enough
    increase_size_mat_mp(MP, D->rows, D->cols);
  }

  MP->rows = D->rows;
  MP->cols = D->cols;

  // find the rational numbers & convert to multiprecision
  for (i = 0; i < MP->rows; i++)
    for (j = 0; j < MP->cols; j++)
      comp_d_to_mp_rat(&MP->entry[i][j], RAT[i][j], &D->entry[i][j], digits, prec, need_to_init_mp, need_to_init_rat); // store the rational numbers out to 'digits' digits

  return;
}

/*****************************************************/
void comp_d_to_mp_rat(comp_mp MP, mpq_t *RAT, comp_d D, int digits, int prec, int need_to_init_mp, int need_to_init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts comp_d to mpq_t & comp_mp using digits & prec *
\***************************************************************/
{
  int mantissa_string_size;
  size_t size, numer_zeros, denom_zeros;
  long int exp;
  char *strFrac = NULL, *strOut = NULL;

  // error check on digits
  if (digits <= 0 || digits >= 16)
  { // set digits to 15 - that is, 15 digits after the decimal place for a total of 16 digits displayed
    digits = 15;
  }

  // find the size of the mantissa
  mantissa_string_size = 3 + digits; // sign + 1st digit + digits + '\0'
  // allocate for strOut
  strOut = (char *)bmalloc(mantissa_string_size * sizeof(char));

  // do real part

  // convert to mantissa and exp
  d_get_str(strOut, &exp, digits, D->r);

  // find the sizes
  size = outStr_to_frac_size(&numer_zeros, &denom_zeros, strOut, exp);
  // allocate memory
  strFrac = (char *)bmalloc(size * sizeof(char));
  // setup strFrac
  outStr_to_frac(strFrac, strOut, numer_zeros, exp, denom_zeros);

  // set RAT[0] using strFrac
  if (need_to_init_rat)
  { // initialize
    mpq_init(RAT[0]);
  }
  mpq_set_str(RAT[0], strFrac, 10); // in base 10
  mpq_canonicalize(RAT[0]);

  // free strFrac
  free(strFrac);

  // do imaginary part

  // convert to mantissa and exp
  d_get_str(strOut, &exp, digits, D->i);

  // find the sizes
  size = outStr_to_frac_size(&numer_zeros, &denom_zeros, strOut, exp);
  // allocate memory
  strFrac = (char *)bmalloc(size * sizeof(char));
  // setup strFrac
  outStr_to_frac(strFrac, strOut, numer_zeros, exp, denom_zeros);

  // set RAT[1] using strFrac
  if (need_to_init_rat)
  { // initialize
    mpq_init(RAT[1]);
  }
  mpq_set_str(RAT[1], strFrac, 10); // in base 10
  mpq_canonicalize(RAT[1]);

  if (need_to_init_mp)
  { // initialize MP to correct precision
    init_mp2(MP, prec);
  }
  else if (mpf_get_prec(MP->r) != prec)
  { // change precision
    change_prec_mp(MP, prec);
  }

  // set MP
  mpf_set_q(MP->r, RAT[0]);
  mpf_set_q(MP->i, RAT[1]);

  // adjust the _d so that they match
  D->r = mpf_get_d(MP->r);
  D->i = mpf_get_d(MP->i);

  return;
}

/*****************************************************/
void initMP(int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  mpf_set_default_prec(prec);

  if (mp_need_init)
  { // everything needs allocated & originally initialized
    int i, max = max_threads();

    _tempMPF1 = (mpf_t *)bmalloc(max * sizeof(mpf_t));
    _tempMPF2 = (mpf_t *)bmalloc(max * sizeof(mpf_t));
    _tempMPF3 = (mpf_t *)bmalloc(max * sizeof(mpf_t));
  
    for (i = 0; i < max; i++)
    {
      mpf_init(_tempMPF1[i]);
      mpf_init(_tempMPF2[i]);
      mpf_init(_tempMPF3[i]);
    }

    mp_need_init = 0;
  }
  else
  { // changing precision for the associated variables
    int oid = thread_num();
    mpf_set_prec(_tempMPF1[oid], prec); 
    mpf_set_prec(_tempMPF2[oid], prec);
    mpf_set_prec(_tempMPF3[oid], prec);
  }
}

/*****************************************************/
void clearMP()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, max = max_threads();

  for (i = 0; i < max; i++)
  { // clear
    mpf_clear(_tempMPF1[i]); 
    mpf_clear(_tempMPF2[i]);
    mpf_clear(_tempMPF3[i]); 
  }
  // if anything needs done after this clear, it needs initialized again
  mp_need_init = 1;
}

/******************************************************/
double maxDiffMat_mp(mat_mp A, mat_mp B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the maximum modulus of the entry in A - B      *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, m = A->rows, n = A->cols;
  double tempD, retVal = 0;
  comp_mp tempC;

  // do error checking
  if (A->rows != B->rows || A->cols != B->cols)
  {
    printf("ERROR: Attempting to subtract matrices with dimensions that do not match!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  init_mp(tempC);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      sub_mp(tempC, &A->entry[i][j], &B->entry[i][j]);
      tempD = d_abs_mp(tempC);
      if (tempD > retVal)
        retVal = tempD;
    }

  clear_mp(tempC);

  return retVal;
}

/******************************************************/
double infNormVec_mp(vec_mp X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the vector                *
* NOTES:                                                        *
\***************************************************************/
{
  int j; double temp, max = 0.0;

  for (j = 0; j < X->size; j++)
  {
    temp = d_abs_mp(&X->coord[j]);
    if (max < temp)
      max = temp;
  }
  return max;
}

void infNormVec_mp2(mpf_t norm, vec_mp X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the vector                *
* NOTES:                                                        *
\***************************************************************/
{
  int j; 
  comp_mp tempComp;
  init_mp(tempComp);

  mpf_set_ui(norm, 0);

  for (j = 0; j < X->size; j++)
  {
    mpf_mul(tempComp->r, X->coord[j].r, X->coord[j].r);
    mpf_mul(tempComp->i, X->coord[j].i, X->coord[j].i);
    mpf_add(tempComp->r, tempComp->r, tempComp->i);

    if (mpfr_less_p(norm, tempComp->r))
      mpf_set(norm, tempComp->r);
  }
  mpf_sqrt(norm, norm);

  clear_mp(tempComp);

  return;
}

/******************************************************/
double infNormMatRow_mp(mat_mp A, int i)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the ith row of the matrix *
* NOTES:                                                        *
\***************************************************************/
// (JDH - 6/21/06) - rewritten to avoid evaluating d_abs_mp twice
{
  int j; double temp, max = 0.0;

  for (j = 0; j < A->cols; j++)
  {
    temp = d_abs_mp(&A->entry[i][j]);
    if (max < temp)
      max = temp;
  }
  return max;
}

/******************************************************/
double infNormMat_mp(mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the matrix                *
* NOTES:                                                        *
\***************************************************************/
// (JDH - 6/21/06) - rewritten to avoid repeatedly setting things to 0
{
  int  i, j, numRows = M->rows, numCols = M->cols;
  double norm, max = 0;
  comp_mp tempComp;
  mpf_t tempNorm;

  mpf_init(tempNorm);
  init_mp(tempComp);

  for (i = 0; i < numRows; i++)
  {
    mpf_set_ui(tempNorm, 0);
    for (j = 0; j < numCols; j++)
    { 
      mpf_mul(tempComp->r, M->entry[i][j].r, M->entry[i][j].r);
      mpf_mul(tempComp->i, M->entry[i][j].i, M->entry[i][j].i);
      mpf_add(tempComp->r, tempComp->r, tempComp->i);
      mpf_add(tempNorm, tempNorm, tempComp->r);
    }
    norm = mpf_get_d(tempNorm);
    
    if (max < norm)
      max = norm;
  }
  max = sqrt(max);

  mpf_clear(tempNorm);
  clear_mp(tempComp);

  return max;
}

/******************************************************/
void infNormMat_mp2(mpf_t norm, mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the infinity norm of the matrix                *
* NOTES:                                                        *
\***************************************************************/
{
  int  i, j, numRows = M->rows, numCols = M->cols;
  mpf_t tempNorm;
  comp_mp tempComp;

  mpf_init(tempNorm);
  init_mp(tempComp);

  mpf_set_ui(norm, 0);

  for (i = 0; i < numRows; i++)
  {
    mpf_set_ui(tempNorm, 0);
    for (j = 0; j < numCols; j++)
    {
      mpf_mul(tempComp->r, M->entry[i][j].r, M->entry[i][j].r);
      mpf_mul(tempComp->i, M->entry[i][j].i, M->entry[i][j].i);
      mpf_add(tempComp->r, tempComp->r, tempComp->i);
      mpf_add(tempNorm, tempNorm, tempComp->r);
    }
    if (mpfr_less_p(norm, tempNorm))
      mpf_set(norm, tempNorm);
  }
  mpf_sqrt(norm, norm);

  mpf_clear(tempNorm);
  clear_mp(tempComp);

  return;
}

double frobNormMat_mp(mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: the Frobenius norm of the matrix               *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, numRows = M->rows, numCols = M->cols;
  double norm = 0;
  mpf_t norm_mp, tempNorm;

  mpf_init(norm_mp);
  mpf_init(tempNorm);
  
  for (i = 0; i < numRows; i++)
    for (j = 0; j < numCols; j++)
    {
      norm_sqr_mp(tempNorm, &M->entry[i][j]);
      mpf_add(norm_mp, norm_mp, tempNorm);
    }

  mpf_sqrt(norm_mp, norm_mp);
  norm = mpf_get_d(norm_mp);

  mpf_clear(norm_mp);
  mpf_clear(tempNorm);

  return norm;
}


/******************************************************/
double conditionNumber_mp(mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  if ((M->rows == 1) && (M->cols == 1))
  {
    return 1.0;
  }
  else
  { // to avoid having 2 copies of the same function
    double tempD1, tempD2;
    return conditionNumber_and_norms_mp(M, &tempD1, &tempD2);
  }
}

/******************************************************/
double conditionNumber_and_norms_mp(mat_mp M, double *norm_M, double *norm_M_inv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: condition number of matrix with norms of it &  *
* its inverse                                                   *
* NOTES: Uses SVD                                               *
\***************************************************************/
{
  int prec = (int) mpf_get_prec(M->entry[0][0].r);
  mpf_t min_sv;
  double condNum;

  mpf_init2(min_sv, prec);

  // find the infinity norm of the matrix
  *norm_M = infNormMat_mp(M);

  // find the minimum singular value
  min_svd_mp_prec(min_sv, M, prec);

  if (mpf_cmp_ui(min_sv, 0) > 0)
  {
    mpf_ui_div(min_sv, 1, min_sv);
    *norm_M_inv = mpf_get_d(min_sv);
  }
  else
  { // report norm of the inverse as 1 / epsilon
    prec = prec_to_digits(prec);
    if (prec > 300)
      *norm_M_inv = 1e300;
    else
      *norm_M_inv = pow(10, prec);
  }
    
  // calculate CN
  condNum = *norm_M * (*norm_M_inv);

  if (condNum < 1)
    condNum = 1;

  mpf_clear(min_sv);

  return condNum;
}

/********************************************************/
double d_abs_mp(comp_mp a)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: sqrt(a->r * a->r + a->i * a->i)                *
* NOTES:                                                        *
\***************************************************************/
{
  double retVal;

  mpf_t re, im;
  mpf_init(re); 
  mpf_init(im);

  mpf_mul(re, a->r, a->r);
  mpf_mul(im, a->i, a->i);
  mpf_add(re, re, im);
  mpf_sqrt(re, re);

  retVal = mpf_get_d(re);

  mpf_clear(re);
  mpf_clear(im);

  return retVal;
}

/********************************************************/
void mul_mat_vec_d(vec_d Res, mat_d M, vec_d X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Res = M*X                                              *
\***************************************************************/
// (JDH - 12/19/06) redone algorithm to make as efficient as possible
{
  int i, j, rows = M->rows, cols = M->cols;

  if (cols != X->size)
  {
    printf("WARNING:  Attempting to multiply a matrix (%d x %d) with a vector (%d) in which the dimensions do not match!\n", rows, cols, X->size);
  }

  if (Res != X) // do the multiplication in place
  { // setup Res
    change_size_vec_d(Res, rows);
 
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Res->coord[i], &M->entry[i][j], &X->coord[j]);
      }
    }
  }
  else // Res == X
  { // need to use a temporary vector
    vec_d v;
    // copy X to v
    init_vec_d(v, X->size);
    vec_cp_d(v, X);

    // setup Res
    change_size_vec_d(Res, rows);

    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Res->coord[i], &M->entry[i][j], &v->coord[j]);
      }
    }
    clear_vec_d(v);
  }
  // set the size
  Res->size = rows;

  return;
}

/********************************************************/
void mul_mat_vec_mp(vec_mp Res, mat_mp M, vec_mp X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Res = M*X                                              *
\***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;

  if (cols != X->size)
  {
    printf("WARNING:  Attempting to multiply a matrix (%d x %d) with a vector (%d) in which the dimensions do not match!\n", rows, cols, X->size);
  }

  if (Res != X) // do the multiplication in place
  { // setup Res
    change_size_vec_mp(Res, rows);

    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Res->coord[i], &M->entry[i][j], &X->coord[j]);
      }
    }
  }
  else // Res == X
  { // need to use a temporary vector
    vec_mp tempVec;
    // copy X to tempVec
    init_vec_mp(tempVec, X->size);
    vec_cp_mp(tempVec, X);

    // setup Res
    change_size_vec_mp(Res, rows);

    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Res->coord[i], &M->entry[i][j], &tempVec->coord[j]);
      }
    }
    clear_vec_mp(tempVec);
  }
  // set the size
  Res->size = rows;

  return;
}

/********************************************************/
void add_vec_d(vec_d Res, vec_d x, vec_d y)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = x->size;

  // increase size of Res
  increase_size_vec_d(Res, size);
  // do the addition
  for (i = 0; i < size; i++)
    add_d(&Res->coord[i], &x->coord[i], &y->coord[i]);
  // set the size
  Res->size = size;

  return;
}

/********************************************************/
void add_vec_mp(vec_mp Res, vec_mp x, vec_mp y)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = x->size;

  // increase size of Res
  increase_size_vec_mp(Res, size);
  // do the addition
  for (i = 0; i < size; i++)
    add_mp(&Res->coord[i], &x->coord[i], &y->coord[i]);
  // set the size
  Res->size = size;

  return;
}

/********************************************************/
double d_vec_abs_d(vec_d v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i; 
  double s = 0;

  for (i = 0; i < v->size; i++)
    s += d_abs_d(&v->coord[i]);

  return s;
}
    
/********************************************************/
double d_vec_abs_mp(vec_mp v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i; 
  double t;
  mpf_t s, temp;
  
  mpf_init_set_ui(s, 0);
  mpf_init(temp);

  for (i = 0; i < v->size; i++)
  {
    mpf_abs_mp(temp, &v->coord[i]); 
    mpf_add(s, s, temp);
  }

  t = mpf_get_d(s);

  mpf_clear(s);
  mpf_clear(temp);
  
  return t;
}

/********************************************************/
void print_d(FILE *fp, int digits, comp_d Z)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int size;
  char *fmt = NULL;

  // error checking
  if (digits <= 0 || digits >= 16)
  { // set digits to 15 - that is, 15 digits after the decimal place for a total of 16 digits displayed
    digits = 15;
  }

  // find the size needed
  size = 1 + snprintf(NULL, 0, "%%.%de %%.%de", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%de %%.%de", digits, digits);
  // print output
  fprintf(fp, fmt, Z->r, Z->i); 
  // release memory
  free(fmt);

  return;
}

void print_Matlab_d(FILE *fp, int digits, comp_d Z)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print Z using Matlab notation                          *
\***************************************************************/
{
  int size;
  char *fmt = NULL;

  // error checking
  if (digits <= 0 || digits >= 16)
  { // set digits to 15 - that is, 15 digits after the decimal place for a total of 16 digits displayed
    digits = 15;
  }

  // find the size needed
  size = 1 + snprintf(NULL, 0, "%%.%de", digits);
  // allocate size
  fmt = (char *)bmalloc(size * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%de", digits);

  // print real part
  fprintf(fp, fmt, Z->r);

  // print imag part
  if (Z->i >= 0)
    fprintf(fp, "+");
  fprintf(fp, fmt, Z->i);
  fprintf(fp, "*i");

  // release memory
  free(fmt);

  return;
}

/********************************************************/
void printMat_d(FILE *fp, int digits, mat_d M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i, j, rows = M->rows, cols = M->cols;

  for (i = 0; i < rows; i++) 
  {
    for (j = 0; j < cols; j++) 
    {
      fprintf(fp, "<"); 
      print_d(fp, digits, &M->entry[i][j]);
      fprintf(fp, "> ");
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  return;
}

void printMat_Matlab_d(FILE *fp, int digits, mat_d M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints M using Matlab format                           *
\***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;

  if (rows <= 0 || cols <= 0)
  {
    fprintf(fp, "\n");
  }
  else
  { // rows > 0 && cols > 0
    for (i = 0; i < rows; i++)
    {
      if (i == 0) fprintf(fp, "[");

      fprintf(fp, "[");
      for (j = 0; j < cols; j++)
      {
        print_Matlab_d(fp, digits, &M->entry[i][j]);
        if (j + 1 < cols)
          fprintf(fp, ", ");
        else
          fprintf(fp, "]");
      }

      if (i == rows - 1)
       fprintf(fp, "];\n\n");
      else
       fprintf(fp, ";\n");
    }
  }

  return;
}

/********************************************************/
void printVec_d(FILE *fp, int digits, vec_d v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i, size = v->size;
    
  for (i = 0; i < size; i++)
  {
    fprintf(fp, "<");
    print_d(fp, digits, &v->coord[i]);
    fprintf(fp, "> ");
  } 
  fprintf(fp, "\n");
}

void printVec_Matlab_d(FILE *fp, int digits, vec_d v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints v using Matlab format                           *
\***************************************************************/
{
  int i, size = v->size;

  if (size <= 0)
  {
    fprintf(fp, "\n");
  }
  else
  { // size > 0
    for (i = 0; i < size; i++)
    {
      if (i == 0) fprintf(fp, "[");

      fprintf(fp, "[");
      print_Matlab_d(fp, digits, &v->coord[i]);
      fprintf(fp, "]");

      if (i == size - 1)
       fprintf(fp, "];\n\n");
      else
       fprintf(fp, ";\n");
    }
  }

  return;
}

/********************************************************/
void printPoint_d(FILE *fp, int digits, point_d v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  printVec_d(fp, digits, v);

  return;
}

void printPoint_Matlab_d(FILE *fp, int digits, point_d v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints v using Matlab format                           *
\***************************************************************/
{
  printVec_Matlab_d(fp, digits, v);
  
  return;
}

/********************************************************/
void print_mp(FILE *fp, int digits, comp_mp Z)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  char *ch1 = NULL, *ch2 = NULL;
  long e1, e2;

  // error checking
  if (digits < 0)
  { // set digits to 0 so that the full precision is displayed
    digits = 0;
  }

  if (mpfr_number_p(Z->r) && mpfr_number_p(Z->i))
  {
    ch1 = mpf_get_str(NULL, &e1, 10, digits, Z->r);
    ch2 = mpf_get_str(NULL, &e2, 10, digits, Z->i);

    if (ch1[0] != '-')
      if (ch2[0] != '-')
        fprintf(fp, "0.%se%ld 0.%se%ld", ch1, e1, ch2, e2);
      else 
        fprintf(fp, "0.%se%ld -0.%se%ld", ch1, e1, &ch2[1], e2);
    else
      if (ch2[0] != '-')
        fprintf(fp, "-0.%se%ld 0.%se%ld", &ch1[1], e1, ch2, e2);
      else
        fprintf(fp, "-0.%se%ld -0.%se%ld", &ch1[1], e1, &ch2[1], e2);

    free(ch1);
    free(ch2);
  }
  else
    fprintf(fp, "NaN NaN");

  return;
}

void print_rat(FILE *fp, mpq_t *Z)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  mpq_out_str(fp, 10, Z[0]);
  fprintf(fp, " ");
  mpq_out_str(fp, 10, Z[1]);

  return;
}

void print_Matlab_mp(FILE *fp, int digits, comp_mp Z)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print Z using Matlab notation                          *
\***************************************************************/
{
  int base = 10;

  if (mpfr_number_p(Z->r) && mpfr_number_p(Z->i))
  {
    mpf_out_str(fp, base, digits, Z->r);
    if (mpfr_sgn(Z->i) >= 0)
      fprintf(fp, "+");
    mpf_out_str(fp, base, digits, Z->i);
    fprintf(fp, "*i");
  }
  else
    fprintf(fp, "NaN+NaN*i");

  return;
}

/********************************************************/
void printMat_mp(FILE *fp, int digits, mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i, j, rows = M->rows, cols = M->cols;

  for (i = 0; i < rows; i++) 
  {
    for (j = 0; j < cols; j++) 
    {
      fprintf(fp, "<");
      print_mp(fp, digits, &M->entry[i][j]);
      fprintf(fp, "> ");
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");

  return;
}

void printMat_Matlab_mp(FILE *fp, int digits, mat_mp M)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints M using Matlab format                           *
\***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;

  if (rows <= 0 || cols <= 0)
  {
    fprintf(fp, "\n");
  }
  else
  { // rows > 0 && cols > 0
    for (i = 0; i < rows; i++)
    {
      if (i == 0) fprintf(fp, "[");

      fprintf(fp, "[");
      for (j = 0; j < cols; j++)
      {
        print_Matlab_mp(fp, digits, &M->entry[i][j]);
        if (j + 1 < cols)
          fprintf(fp, ", ");
        else
          fprintf(fp, "]");
      }

      if (i == rows - 1)
       fprintf(fp, "];\n\n");
      else
       fprintf(fp, ";\n");
    }
  }

  return;
}

/********************************************************/
void printVec_mp(FILE *fp, int digits, vec_mp v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  int i, size = v->size;
    
  for (i = 0; i < size; i++)
  {  
    fprintf(fp, "<");
    print_mp(fp, digits, &v->coord[i]);
    fprintf(fp, "> ");
  }
  fprintf(fp, "\n");

  return;
}

void printVec_Matlab_mp(FILE *fp, int digits, vec_mp v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints v using Matlab format                           *
\***************************************************************/
{
  int i, size = v->size;

  if (size <= 0)
  {
    fprintf(fp, "\n");
  }
  else
  { // size > 0
    for (i = 0; i < size; i++)
    {
      if (i == 0) fprintf(fp, "[");

      fprintf(fp, "[");
      print_Matlab_mp(fp, digits, &v->coord[i]);
      fprintf(fp, "]");

      if (i == size - 1)
       fprintf(fp, "];\n\n");
      else
       fprintf(fp, ";\n");
    }
  }

  return;
}

/********************************************************/
void printPoint_mp(FILE *fp, int digits, point_mp v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  printVec_mp(fp, digits, v);

  return;
}

void printPoint_Matlab_mp(FILE *fp, int digits, point_mp v)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  printVec_Matlab_mp(fp, digits, v);

  return;
}

/************************************************************************/
void hermiteInterpCW_d(point_d Res, _comp_d *T, _point_d *Y, _point_d *dHdT, comp_d T0, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/************************************************************************/
/* See hermiteInterpCW_mp below for details - this is just a copy. 	*/
/************************************************************************/
{ int    i, j, c, d;
  _comp_d **F = NULL, *z = NULL, *dd = NULL;
  comp_d t1, t2, U;   

  // initialize data
  F = (_comp_d **)bmalloc(2 * n * sizeof(_comp_d *));
  z = (_comp_d *)bmalloc(2 * n * sizeof(_comp_d));
  dd = (_comp_d *)bmalloc(2 * n * sizeof(_comp_d));
  for (i = 0; i < 2*n; i++)
  {
    F[i] = (_comp_d *)bmalloc(2 * n * sizeof(_comp_d));
  }

  d = Y[0].size;

  // increase size on Res
  increase_size_point_d(Res, d);
  // set the size
  Res->size = d;

  for (c=0; c<d; c++) 
  {
    for (i=0; i<n ;i++) {
      set_d(&F[2*i][0], &Y[i].coord[c]);      /*  F[2*i][0]    = Y[i][c];    */
      set_d(&F[2*i+1][0], &Y[i].coord[c]);    /*  F[2*i+1][0]  = Y[i][c];    */
      set_d(&F[2*i+1][1], &dHdT[i].coord[c]); /*  F[2*i+1][1]  = dHdT[i][c]; */
      set_d(&z[2*i], &T[i]);                  /*  z[2*i]       = T[i];       */
      set_d(&z[2*i+1], &T[i]);                /*  z[2*i+1]     = T[i];       */
    }

    for (i=1; i<n; i++) {
      sub_d(t1, &F[2*i][0], &F[2*i-1][0]);
      sub_d(t2, &z[2*i], &z[2*i-1]);
      div_d(&F[2*i][1], t1, t2);
      /* F[2*i][1] = (F[2*i][0]-F[2*i-1][0])/(z[2*i]-z[2*i-1]); */
    }
    for (i=2; i<2*n; i++)
      for (j=2; j<=i; j++) {
        sub_d(t1, &F[i][j-1], &F[i-1][j-1]);
        sub_d(t2, &z[i], &z[i-j]);
        div_d(&F[i][j], t1, t2);
        /* F[i][j] = (F[i][j-1]-F[i-1][j-1])/(z[i]-z[i-j]); */
      }
    for (i=0; i<2*n; i++) {
      set_d(&dd[i], &F[i][i]);  /*  dd[i] = F[i][i]; */
    }
    set_d(U, &dd[2*n-1]); /* U = dd[2*n-1]; */
    i = n-1;

    while (i>=1) {
      sub_d(t1, T0, &T[i]);
      mul_d(t1, t1, U);
      add_d(t1, t1, &dd[2*i]);
      sub_d(t2, T0, &T[i-1]);
      mul_d(t1, t1, t2);
      add_d(U, t1, &dd[2*i-1]);
      /* U = ( U*(T0-T[i])+dd[2*i] )*(T0-T[i-1]) + dd[2*i-1]; */
      i--;
    }
    sub_d(t1, T0, &T[0]);
    mul_d(t1, t1, U);
    add_d(&Res->coord[c], t1, &dd[0]);
    /* Res[c] = U*(T0-T[0])+dd[0]; */
  }

  // clear
  for (i = 0; i < 2*n; i++)
  {
    free(F[i]);
  }
  free(F);
  free(z);
  free(dd);

  return;
}

/************************************************************************/
void hermiteInterpCW_mp(point_mp Res, _comp_mp *T, _point_mp *Y, _point_mp *dHdT, comp_mp T0, int n)
/************************************************************************/
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
/* Input: T values, T[0],...,T[n-1],
        points Y[i] with H(T[i], Y[i]) = 0,
        and the derivatives w.r.t. 'T', dHdT[i].
        A point Z.
 Output: Since H(t, y) is analytic in t in a neighborhood,
        it has a power series expansion in each coordinate, on
        that neighborhood: H(t,y) = [ f0(t), f1(t),..., fd(t) ].
        For j=0...d we will use the given T-values, the j-th coordinates 
        of the given points on H(t,y), and the j-th coordinates of the 
        derivatives to interpolate polynomials
        [ g0(t),..., gd(t)] ~ [ f0(t),..., fd(t) ]. 
        Finally, we will evaluate and return:
          Res = [ g0(Z[0]), g1(Z[1]), ..., gd(Z[d]) ].	
*/
/************************************************************************/
{ 
  int i, j, c, d;

  _comp_mp **__F = NULL, *__z = NULL, *__dd = NULL;
  comp_mp   __U, __t1, __t2;

  init_mp(__U); init_mp(__t1); init_mp(__t2);

  __F = (_comp_mp **)bmalloc(2 * n * sizeof(_comp_mp *));
  __z = (_comp_mp *)bmalloc(2 * n * sizeof(_comp_mp));
  __dd = (_comp_mp *)bmalloc(2 * n * sizeof(_comp_mp));
  for (i=0; i<2*n; i++) 
  {
    init_mp(&__z[i]); init_mp(&__dd[i]);
    __F[i] = (_comp_mp *)bmalloc(2 * n * sizeof(_comp_mp));
    for (j=0; j<2*n; j++)
      init_mp(&__F[i][j]);
  }

  d = Y[0].size;

  // increase size on Res
  increase_size_point_mp(Res, d);
  // set the size
  Res->size = d;

  for (c=0; c<d; c++) { 
    for (i=0; i<n ;i++) {
      set_mp(&__F[2*i][0], &Y[i].coord[c]);      /*  F[2*i][0]    = Y[i][c];    */
      set_mp(&__F[2*i+1][0], &Y[i].coord[c]);    /*  F[2*i+1][0]  = Y[i][c];    */
      set_mp(&__F[2*i+1][1], &dHdT[i].coord[c]); /*  F[2*i+1][1]  = dHdT[i][c]; */
      set_mp(&__z[2*i], &T[i]);                  /*  z[2*i]       = T[i];       */
      set_mp(&__z[2*i+1], &T[i]);                /*  z[2*i+1]     = T[i];       */
    }
  
    for (i=1; i<n; i++) {
      sub_mp(__t1, &__F[2*i][0], &__F[2*i-1][0]);
      sub_mp(__t2, &__z[2*i], &__z[2*i-1]);
      div_mp(&__F[2*i][1], __t1, __t2);
      /* F[2*i][1] = (F[2*i][0]-F[2*i-1][0])/(z[2*i]-z[2*i-1]); */
    }  
    for (i=2; i<2*n; i++)
      for (j=2; j<=i; j++) {
        sub_mp(__t1, &__F[i][j-1], &__F[i-1][j-1]);
	sub_mp(__t2, &__z[i], &__z[i-j]);
	div_mp(&__F[i][j], __t1, __t2);
        /* F[i][j] = (F[i][j-1]-F[i-1][j-1])/(z[i]-z[i-j]); */
      }
    for (i=0; i<2*n; i++) {
      set_mp(&__dd[i], &__F[i][i]);  /*  dd[i] = F[i][i]; */
    }
    set_mp(__U, &__dd[2*n-1]); /* U = dd[2*n-1]; */
    i = n-1;
 
    while (i>=1) {
      sub_mp(__t1, T0, &T[i]);
      mul_mp(__t1, __t1, __U);
      add_mp(__t1, __t1, &__dd[2*i]);
      sub_mp(__t2, T0, &T[i-1]);
      mul_mp(__t1, __t1, __t2);
      add_mp(__U, __t1, &__dd[2*i-1]);
      /* U = ( U*(T0-T[i])+dd[2*i] )*(T0-T[i-1]) + dd[2*i-1]; */
      i--;
    }
    sub_mp(__t1, T0, &T[0]);
    mul_mp(__t1, __t1, __U);
    add_mp(&Res->coord[c], __t1, &__dd[0]);
    /* Res[c] = U*(T0-T[0])+dd[0]; */
  }

  clear_mp(__U); clear_mp(__t1); clear_mp(__t2);
  for (i = 0; i < 2*n; i++)
  {
    clear_mp(&__z[i]); clear_mp(&__dd[i]);
    for (j=0; j<2*n; j++)
      clear_mp(&__F[i][j]);

    free(__F[i]);
  }
  free(__z);
  free(__dd);
  free(__F);

  return;  
}
	
void pow_rdouble_d(comp_d res, comp_d base, double e)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  double theta, rho;

  // find theta
  theta = atan2(base->i, base->r);
  // find rho
  rho = d_abs_d(base);
  // find rho^e
  rho = pow(rho, e);

  // res->r = rho^e * cos(e * theta)
  res->r = rho * cos(e * theta);
  // res->r = rho^e * sin(e * theta)
  res->i = rho * sin(e * theta);

  return;
}

void pow_rdouble_mp(comp_mp res, comp_mp base, double e)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{ 
  mpf_t theta, rho;

  mpf_init(theta);
  mpf_init(rho);

  // find rho
  mpf_abs_mp(rho, base);

  // find rho^e
  mpf_set_d(theta, e);
  mpfr_pow(rho, rho, theta, __gmp_default_rounding_mode);

  // find e *theta
  mpfr_atan2(res->r, base->i, base->r, __gmp_default_rounding_mode);
  mpf_mul(theta, theta, res->r);

  // res->r = (rho^e) * cos(e * theta)
  mpfr_cos(res->r, theta, __gmp_default_rounding_mode);
  mpf_mul(res->r, res->r, rho);

  // res->i = (rho^e) * sin(e * theta)
  mpfr_sin(res->i, theta, __gmp_default_rounding_mode);
  mpf_mul(res->i, res->i, rho);

  mpf_clear(theta);
  mpf_clear(rho);

  return;
}

void pow_rmpf_mp(comp_mp res, comp_mp base, mpf_t e)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  mpf_t theta, rho;

  mpf_init(theta);
  mpf_init(rho);

  // find rho
  mpf_abs_mp(rho, base); 

  // find rho^2
  mpfr_pow(rho, rho, e, __gmp_default_rounding_mode);

  // find e * theta
  mpfr_atan2(theta, base->i, base->r, __gmp_default_rounding_mode);
  mpf_mul(theta, theta, e);

  // res->r = (rho^e) * cos(e * theta)
  mpfr_cos(res->r, theta, __gmp_default_rounding_mode);
  mpf_mul(res->r, res->r, rho);

  // res->i = (rho^e) * sin(e * theta)
  mpfr_sin(res->i, theta, __gmp_default_rounding_mode);
  mpf_mul(res->i, res->i, rho);

  mpf_clear(theta);
  mpf_clear(rho);

  return;
}

void twoNormVec_d(vec_d V, double *retVal)  /* Computes the 2-norm of a vector which is stored as a vec_d. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = V->size;
  double sum = 0;

  for (i = 0; i < size; i++)  /* Square the real and imaginary parts of each element of the vector and add them to sum. */
    sum += norm_sqr_d(&V->coord[i]);

  *retVal = sqrt(sum);  /* Take the square root of sum and that's what you return. */

  return;
}

void outerProduct_d(mat_d Res, vec_d V, vec_d W)  /* Computes the outer product of V and W. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;
  comp_d tmp;

  // set Res to the correct size
  change_size_mat_d(Res, V->size, W->size);
  Res->rows = V->size;
  Res->cols = W->size;

  for (i = 0; i < Res->cols; i++)
    for (j = 0; j < Res->rows; j++)
    {
      conjugate_d(tmp, &V->coord[i]);  /* Notice the conjugation!!! */
      mul_d(&Res->entry[j][i], tmp, &W->coord[j]);
    }

  return;
}

void transpose_d(mat_d Res, mat_d M)  /* Stores CONJUGATE transpose of M in Res. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;

  if (Res != M)
  { // setup Res
    change_size_mat_d(Res, cols, rows);

    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_d(&Res->entry[i][j], &M->entry[j][i]); // Notice the conjugation!!
      }
  }
  else // Res = M
  { // need to use a temporary matrix
    mat_d tempMat;
    // copy M to tempMat
    init_mat_d(tempMat, rows, cols);
    mat_cp_d(tempMat, M);

    // setup Res
    change_size_mat_d(Res, cols, rows);
 
    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_d(&Res->entry[i][j], &tempMat->entry[j][i]); // Notice the conjugation!!
      }

    // clear tempMat
    clear_mat_d(tempMat);
  }
  // set the size
  Res->rows = cols;
  Res->cols = rows;
 
  return;
}

void find_opposite_phase_d(comp_d opposite_phase, comp_d number)  /* Given complex number "number", computes a complex number s.t. the product of the two is real. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  double tmp;  

  set_d(opposite_phase, number);

  if (number->i == 0.0)  /* If number is real, just return 1.0. */
  {
    if (number->r < 0.0)
      set_double_d(opposite_phase, -1.0, 0.0);
    if (number->r > 0.0)
      set_double_d(opposite_phase, 1.0, 0.0);
    return;
  }

  if (number->r != 0.0)
  {
    tmp = -1.0*(number->i/number->r);
    opposite_phase->r = cos(atan(tmp)); 
    opposite_phase->i = sin(atan(tmp));
  }

  return;
}

void make_vec_random_d(vec_d p, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  double norm, max = 0;

  // setup p
  change_size_vec_d(p, size);
  p->size = size;

  if (size > 0)
  {
    for (i = 0; i < size; i++)
    { // find ith random coord and its norm
      norm = get_comp_rand_d(&p->coord[i]);
      if (norm > max)
        max = norm;
    }

    // normalize to unit inf-norm
    norm = 1 / max;
    for (i = 0; i < size; i++)
    { 
      mul_rdouble_d(&p->coord[i], &p->coord[i], norm);
    }
  }
  
  return;
}

void make_vec_random_d2(vec_d p, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: random unit vector using 2-norm                        *
\***************************************************************/
{
  int i;
  double norm = 0;

  // setup p
  change_size_vec_d(p, size);
  p->size = size;

  if (size > 0)
  {
    for (i = 0; i < size; i++)
    { // find ith random coord and its norm
      get_comp_rand_d(&p->coord[i]);
      norm += norm_sqr_d(&p->coord[i]);
    }

    // normalize to unit 2-norm
    norm = 1 / sqrt(norm);
    for (i = 0; i < size; i++)
    {
      mul_rdouble_d(&p->coord[i], &p->coord[i], norm);
    }
  }

  return;
}

void make_matrix_ID_d(mat_d A, int rows, int cols)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;

  // make A the same size
  change_size_mat_d(A, rows, cols);
  A->rows = rows;
  A->cols = cols;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      set_double_d(&A->entry[i][j], i == j, 0);
    }

  return;
}

void make_matrix_random_d(mat_d A, int rows, int cols)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random conjugate-orthonormal matrix          *
* if rows <= cols, Q*QH = I_rows                                *
* if rows >= cols, QH*Q = I_cols                                *
\***************************************************************/
{
  if (rows > 0 && cols > 0)
  { // make sure that A has both rows & columns
    if (cols == rows)
    { // generate a cols x rows (square) conjugate-orthonormal matrix 
      make_square_matrix_orth_d(A, cols); // cols == rows
    }
    else if (cols < rows)
    { // generate a rows x rows conjugate-orthonormal matrix and simply ignore the rest of the matrix
      make_square_matrix_orth_d(A, rows);
      A->cols = cols;
    }
    else // cols > rows
    { // generate a cols x cols conjugate-orthonormal matrix and shrink to cols x rows and then conjugate-transpose to rows x cols
      make_square_matrix_orth_d(A, cols);
      A->cols = rows;
      transpose_d(A, A);
    }
  }
  else
  { // setup A
    change_size_mat_d(A, rows, cols);
    A->rows = rows;
    A->cols = cols;
  }

  return;
}

void make_matrix_random_real_d(mat_d A, int rows, int cols)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random real orthonormal matrix               *
* if rows <= cols, Q*QT = I_rows                                *
* if rows >= cols, QT*Q = I_cols                                *
\***************************************************************/
{
  if (rows > 0 && cols > 0)
  { // make sure that A has both rows & columns
    if (cols == rows)
    { // generate a cols x rows (square) orthonormal real matrix
      make_square_matrix_orth_real_d(A, cols); // cols == rows
    }
    else if (cols < rows)
    { // generate a rows x rows orthonormal real matrix and simply ignore the rest of the matrix
      make_square_matrix_orth_real_d(A, rows);
      A->cols = cols;
    }
    else // cols > rows
    { // generate a cols x cols orthonormal real matrix and shrink to cols x rows and then transpose to rows x cols
      make_square_matrix_orth_real_d(A, cols);
      A->cols = rows;
      transpose_d(A, A);
    }
  }
  else
  { // setup A
    change_size_mat_d(A, rows, cols);
    A->rows = rows;
    A->cols = cols;
  }

  return;
}

void make_square_matrix_orth_d(mat_d A, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random square conjugate-orthonormal matrix   *
* using Stewart's method: A is product of Householder matrices  *
\***************************************************************/
{
  int i, j, k;
  double norm, tempD;
  comp_d gamma;
  vec_d x, u;

  // initialize x & u
  init_vec_d(x, size);
  init_vec_d(u, size);

  // set to size x size 'gamma * Id'
  change_size_mat_d(A, size, size);
  A->rows = A->cols = size;
  // find gamma and normalize it
  norm = get_comp_rand_d(gamma);
  gamma->r /= norm;
  gamma->i /= norm;
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (i == j)
      {
        set_d(&A->entry[i][j], gamma);
      }
      else
      {
        set_zero_d(&A->entry[i][j]);
      }

  for (j = 0; j < size - 1; j++)
  { // generate a random vector u of unit length in 2-norm' (u[j:size])
    norm = 0;
    for (k = j; k < size; k++)
    {
      get_comp_rand_d(&u->coord[k]);
      norm += norm_sqr_d(&u->coord[k]);
    }
    // find normalization constant and normalize
    norm = 1 / sqrt(norm);
    for (k = j; k < size; k++)
    { // normalize u[k]
      mul_rdouble_d(&u->coord[k], &u->coord[k], norm);
    }

    // compute x = A*u, taking into account sizes 
    for (i = 0; i < size; i++)
    {
      set_zero_d(&x->coord[i]);
      for (k = j; k < size; k++)
      {
        sum_mul_d2(&x->coord[i], &A->entry[i][k], &u->coord[k], tempD);
      }
    }

    // then update A := A - 2*(A*u)*(uH)
    for (i = 0; i < size; i++)
      for (k = j; k < size; k++)
      {
        householder_mult_right_d2(&A->entry[i][k], &A->entry[i][k], &x->coord[i], &u->coord[k], tempD);
      }
  }

  // clear x & u
  clear_vec_d(x);
  clear_vec_d(u);

  return;
}

void make_square_matrix_orth_real_d(mat_d A, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random square orthonormal real matrix        *
* using Stewart's method: A is product of Householder matrices  *
\***************************************************************/
{
  int i, j, k;
  double norm, tempD;
  comp_d gamma;
  vec_d x, u;

  // initialize x & u
  init_vec_d(x, size);
  init_vec_d(u, size);

  // set to size x size 'gamma * Id'
  change_size_mat_d(A, size, size);
  A->rows = A->cols = size;
  // set gamma to be either 1 or -1 randomly
  set_one_d(gamma);
  if (rand() % 2) 
    neg_d(gamma, gamma);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (i == j)
      {
        set_d(&A->entry[i][j], gamma);
      }
      else
      {
        set_zero_d(&A->entry[i][j]);
      }

  for (j = 0; j < size - 1; j++)
  { // generate a random vector u of unit length in 2-norm' (u[j:size])
    norm = 0;
    for (k = j; k < size; k++)
    {
      get_comp_rand_real_d(&u->coord[k]);
      norm += norm_sqr_d(&u->coord[k]);
    }

    // find normalization constant and normalize
    norm = 1 / sqrt(norm);
    for (k = j; k < size; k++)
    { // normalize u[k]
      mul_rdouble_d(&u->coord[k], &u->coord[k], norm);
    }

    // compute x = A*u, taking into account sizes
    for (i = 0; i < size; i++)
    {
      set_zero_d(&x->coord[i]);
      for (k = j; k < size; k++)
      {
        sum_mul_d2(&x->coord[i], &A->entry[i][k], &u->coord[k], tempD);
      }
    }

    // then update A := A - 2*(A*u)*(uH)
    for (i = 0; i < size; i++)
      for (k = j; k < size; k++)
      {
        householder_mult_right_d2(&A->entry[i][k], &A->entry[i][k], &x->coord[i], &u->coord[k], tempD);
      }
  }

  // clear x & u
  clear_vec_d(x);
  clear_vec_d(u);

  return;
}

void make_elem_d(mat_d A, int num1, int num2)  /* Make A into an ID matrix of its present dimensions, but then swap rows num1 and num2 and cols num1 and num2 to create an elementary "row swap" or "column swap" matrix. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  make_matrix_ID_d(A, A->rows, A->cols);
  set_double_d(&A->entry[num1][num1], 0.0, 0.0);
  set_double_d(&A->entry[num2][num2], 0.0, 0.0);
  set_double_d(&A->entry[num1][num2], 1.0, 0.0);
  set_double_d(&A->entry[num2][num1], 1.0, 0.0);

  return;
}

void mat_mul_d(mat_d Res, mat_d A, mat_d B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, rows = A->rows, inner = A->cols, cols = B->cols;

  if (A->cols != B->rows)
    printf("WARNING:  Attempting to multiply matrices with dimensions that do not match!\n");

  if (Res != A)
  {
    if (Res != B) // do multiplication in place
    { // setup Res
      change_size_mat_d(Res, rows, cols);

      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_d(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_d(&Res->entry[i][j], &A->entry[i][k], &B->entry[k][j]); // Res_i,j += A_i,k * B_k,j
          }
        }
    }
    else // Res == B 
    { // copy B to M
      mat_d M;
      init_mat_d(M, B->rows, B->cols);
      mat_cp_d(M, B);

      // setup Res
      change_size_mat_d(Res, rows, cols);

      // multiply Res = A * M
      for (j = 0; j < cols; j++)
        for (i = 0; i < rows; i++)
        {
          set_zero_d(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_d(&Res->entry[i][j], &A->entry[i][k], &M->entry[k][j]); // Res_i,j += A_i,k * M_k,j
          }
        }

      // clear M
      clear_mat_d(M);
    }
  }
  else // Res == A 
  {
    if (Res != B)
    { // copy A to M
      mat_d M; 
      init_mat_d(M, A->rows, A->cols);
      mat_cp_d(M, A);

      // setup Res
      change_size_mat_d(Res, rows, cols);

      // multiply Res = M * B
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_d(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_d(&Res->entry[i][j], &M->entry[i][k], &B->entry[k][j]); // Res_i,j += M_i,k * B_k,j
          }
        }

      // clear M
      clear_mat_d(M);
    }
    else // Res == A == B
    { // copy A == B to M
      mat_d M;
      init_mat_d(M, A->rows, A->cols); 
      mat_cp_d(M, A); // M == A == B

      // setup Res
      change_size_mat_d(Res, rows, cols);

      // multiply Res = M * M
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_d(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_d(&Res->entry[i][j], &M->entry[i][k], &M->entry[k][j]); // Res_i,j += M_i,k * M_k,j
          }
        }
    }
  }
  // set the size
  Res->rows = rows;
  Res->cols = cols;

  return;
}

void Gram_Schmidt_d(mat_d A)  /* Perform Gram Schmidt orthogonalization on A, overwriting A with result. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int colnum, j, k, transposed = 0;
  vec_d tmpVec;
  mat_d Q, R;
  comp_d tmp;
  double norm;
  double TOL = 0.000000000000001;  // VERY BAD NEWS - SEE LATER COMMENTS!!!!!!!!!!!!!!!!!!!!!

  if (A->rows < A->cols)
  {
    transposed = 1;
    transpose_d(A, A);
  }

  // initialize
  init_vec_d(tmpVec, A->rows);
  init_mat_d(Q, A->rows, A->cols);
  init_mat_d(R, A->cols, A->cols);

  tmpVec->size = A->rows;
  Q->rows = A->rows;
  Q->cols = R->rows = R->cols = A->cols;

  for (colnum = 0; colnum < A->cols; colnum++)
  {
    for (j = 0; j < tmpVec->size; j++)
      set_d(&tmpVec->coord[j], &A->entry[j][colnum]);
    twoNormVec_d(tmpVec, &norm);
    if (norm == 0.0)
      set_double_d(&R->entry[colnum][colnum], 0.0000000001, 0.0);  // VERY BAD FIX!!!!!!  BETTER ALGORITHM IS NEEDED!!!!!!!!!!!!!!!!!!!!!
    if (norm != 0.0)
      set_double_d(&R->entry[colnum][colnum], norm, 0.0);   
    for (j = 0; j < A->rows; j++)
      div_d(&Q->entry[j][colnum], &A->entry[j][colnum], &R->entry[colnum][colnum]);
    for (j = colnum+1; j < A->cols; j++)
    {  
      set_double_d(&R->entry[colnum][j], 0.0, 0.0);
      for (k = 0; k < A->rows; k++)
      {
	conjugate_d(tmp, &Q->entry[k][colnum]);
	mul_d(tmp, tmp, &A->entry[k][j]);
	add_d(&R->entry[colnum][j], &R->entry[colnum][j], tmp);
      }
      for (k = 0; k < A->rows; k++)
      {
	mul_d(tmp, &R->entry[colnum][j], &Q->entry[k][colnum]);
	sub_d(&A->entry[k][j], &A->entry[k][j], tmp);
      }
    }
  }

  for (colnum = 0; colnum < A->cols; colnum++)  // This sets near-0 columns to 0.  Later, it should be replaced with a call to Gram_Schmidt_mp().!!!!!!
  {
    for (j = 0; j < tmpVec->size; j++)
      set_d(&tmpVec->coord[j], &A->entry[j][colnum]);
    twoNormVec_d(tmpVec, &norm);
    if (norm < TOL)  //THIS IS BIG TROUBLE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
      for (j = 0; j < tmpVec->size; j++)
        set_double_d(&A->entry[j][colnum], 0.0, 0.0);
    }
  }

  if (transposed)
  {
    transpose_d(A, A);
  }

  // clear
  clear_vec_d(tmpVec);
  clear_mat_d(Q);
  clear_mat_d(R);

  return;
}

void twoNormVec_mp(vec_mp V, comp_mp retVal)  /* Computes the 2-norm of a vector */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  set_zero_mp(retVal);
  for (i = 0; i < V->size; i++)  /* Square the real and imaginary parts of each element of the vector and add them */
  {
    norm_sqr_mp(retVal->i, &V->coord[i]);
    mpf_add(retVal->r, retVal->r, retVal->i);
  }
  mpf_sqrt(retVal->r, retVal->r);  /* Take the square root and that's what you return. */
  mpf_set_ui(retVal->i, 0);

  return;
}

void twoNormVec_mp2(vec_mp V, mpf_t retVal)  /* Computes the 2-norm of a vector */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  mpf_t tempMPF;
  mpf_init(tempMPF);

  mpf_set_ui(retVal, 0);
  for (i = 0; i < V->size; i++)  /* Square the real and imaginary parts of each element of the vector and add them */
  {
    norm_sqr_mp(tempMPF, &V->coord[i]);
    mpf_add(retVal, retVal, tempMPF);
  }
  mpf_sqrt(retVal, retVal);  /* Take the square root and that's what you return. */

  mpf_clear(tempMPF);

  return;
}

void outerProduct_mp(mat_mp Res, vec_mp V, vec_mp W)  /* Computes the outer product of V and W. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;
  comp_mp tmp;
  init_mp(tmp);

  // set Res to the correct size
  change_size_mat_mp(Res, V->size, W->size);
  Res->rows = V->size;
  Res->cols = W->size;

  for (i = 0; i < Res->cols; i++)
    for (j = 0; j < Res->rows; j++)
    {
      conjugate_mp(tmp, &V->coord[i]);  /* Notice the conjugation!!! */
      mul_mp(&Res->entry[j][i], tmp, &W->coord[j]); 
    }

  clear_mp(tmp);

  return;
}

void transpose_mp(mat_mp Res, mat_mp M)  /* Stores CONJUGATE transpose of M in Res. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, rows = M->rows, cols = M->cols;

  if (Res != M)
  { // setup Res
    change_size_mat_mp(Res, cols, rows);

    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_mp(&Res->entry[i][j], &M->entry[j][i]); // Notice the conjugation!!
      }
  }
  else // Res = M
  { // need to use a temporary matrix
    mat_mp tempMat;
    // copy M to tempMat
    init_mat_mp2(tempMat, rows, cols, M->curr_prec);
    mat_cp_mp(tempMat, M); 

    // setup Res
    change_size_mat_mp(Res, cols, rows);

    for (i = 0; i < cols; i++)
      for (j = 0; j < rows; j++)
      {
        conjugate_mp(&Res->entry[i][j], &tempMat->entry[j][i]); // Notice the conjugation!!
      }

    // clear tempMat 
    clear_mat_mp(tempMat);
  }
  // set the size
  Res->rows = cols;
  Res->cols = rows;

  return;
}

void find_opposite_phase_mp(comp_mp opposite_phase, comp_mp number)  /* Given complex number "number", computes the complex number s.t. the product of the two is real. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  comp_mp tmp, tmp2;

  init_mp(tmp);  
  init_mp(tmp2);

  set_mp(opposite_phase, number);

  if (mpf_get_d(number->i) == 0.0)  /* If number is real, return 1.0! */
  {
    if (mpf_get_d(number->r) < 0.0)
      set_double_mp2(opposite_phase, -1.0, 0.0);
    if (mpf_get_d(number->r) > 0.0)
      set_double_mp2(opposite_phase, 1.0, 0.0);
  }
  else
  {
    mpf_div(tmp->r, number->i, number->r);
    mul_rdouble_mp(tmp, tmp, -1.0);
    mpfr_atan(tmp2->r, tmp->r, GMP_RNDN);
    mpfr_cos(opposite_phase->r, tmp2->r, GMP_RNDN);
    mpfr_atan(tmp2->r, tmp->r, GMP_RNDN);
    mpfr_sin(opposite_phase->i, tmp2->r, GMP_RNDN);
  }

  clear_mp(tmp);
  clear_mp(tmp2);

  return;
}

void make_matrix_random_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int rows, int cols, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random conjugate-orthonormal matrix          *
* if rows <= cols, Q*QH = I_rows                                *
* if rows >= cols, QH*Q = I_cols                                *
\***************************************************************/
{
  int i, j;
  mat_mp tempMat;

  init_mat_mp2(tempMat, 0, 0, max_prec);

  // initialize - if needed
  if (need_to_init_mp)
    init_mat_mp2(A_mp, rows, cols, curr_prec);

  if (need_to_init_rat)
    init_mat_rat(A_rat, rows, cols);

  // set to the correct size
  change_size_mat_d(A_d, rows, cols);
  change_size_mat_mp(A_mp, rows, cols);

  if (rows > 0 && cols > 0)
  { // make sure that A has both rows & columns
    if (cols == rows)
    { // generate a cols x rows (square) conjugate-orthonormal matrix
      make_square_matrix_orth_mp(tempMat, cols, max_prec); // cols == rows
    }
    else if (cols < rows)
    { // generate a rows x rows conjugate-orthonormal matrix and simply ignore the rest of the matrix
      make_square_matrix_orth_mp(tempMat, rows, max_prec);
      tempMat->cols = cols;
    }
    else // cols > rows
    { // generate a cols x cols conjugate-orthonormal matrix and shrink to cols x rows and then conjugate-transpose to rows x cols
      make_square_matrix_orth_mp(tempMat, cols, max_prec);
      tempMat->cols = rows;
      transpose_mp(tempMat, tempMat);
    }

    // copy tempMat to A_d, A_mp, A_rat
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpf_t_to_rat(A_rat[i][j][0], tempMat->entry[i][j].r);
        mpf_t_to_rat(A_rat[i][j][1], tempMat->entry[i][j].i);
        mpf_set_q(A_mp->entry[i][j].r, A_rat[i][j][0]);
        mpf_set_q(A_mp->entry[i][j].i, A_rat[i][j][1]);
        A_d->entry[i][j].r = mpq_get_d(A_rat[i][j][0]);
        A_d->entry[i][j].i = mpq_get_d(A_rat[i][j][1]);
      }
  }
  // set the size
  A_d->rows = A_mp->rows = rows;
  A_d->cols = A_mp->cols = cols;

  clear_mat_mp(tempMat);

  return;
}

void make_vec_random_mp(vec_mp p, int size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  mpf_t norm, max;

  mpf_init(norm);
  mpf_init(max);
  mpf_set_ui(max, 0);

  // setup p
  change_size_vec_mp(p, size);
  p->size = size;

  if (size > 0)
  {
    for (i = 0; i < size; i++)
    { // find ith random coord and its norm
      get_comp_rand_mp2(norm, &p->coord[i]);
      if (mpf_cmp(norm, max) > 0)
        mpf_set(max, norm);
    }

    // normalize to unit inf-norm
    mpf_ui_div(norm, 1, max);
    for (i = 0; i < size; i++)
    {
      mul_rmpf_mp(&p->coord[i], &p->coord[i], norm);     
    }
  }

  mpf_clear(norm);
  mpf_clear(max);

  return;
}

void make_vec_random_rat(vec_d p_d, vec_mp p_mp, mpq_t **p_rat, int size, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  vec_mp tempVec;

  init_vec_mp2(tempVec, size, max_prec);

  // initialize - if needed
  if (need_to_init_mp)
    init_vec_mp2(p_mp, size, curr_prec);

  if (need_to_init_rat)
  {
    for (i = 0; i < size; i++)
    {
      mpq_init(p_rat[i][0]);
      mpq_init(p_rat[i][1]);
    }
  }

  // set to the correct size
  change_size_vec_d(p_d, size);
  change_size_vec_mp(p_mp, size);

  if (size > 0)
  { // setup tempVec
    tempVec->size = size;
    make_vec_random_mp(tempVec, size);

    // copy tempVec to p_d, p_mp, p_rat
    for (i = 0; i < size; i++)
    {
      mpf_t_to_rat(p_rat[i][0], tempVec->coord[i].r);
      mpf_t_to_rat(p_rat[i][1], tempVec->coord[i].i);
      mpf_set_q(p_mp->coord[i].r, p_rat[i][0]);
      mpf_set_q(p_mp->coord[i].i, p_rat[i][1]);
      p_d->coord[i].r = mpq_get_d(p_rat[i][0]);
      p_d->coord[i].i = mpq_get_d(p_rat[i][1]);
    }
  }
  // set the size
  p_d->size = p_mp->size = size;

  clear_vec_mp(tempVec);

  return;
}
   
void make_matrix_random_mp(mat_mp A, int rows, int cols, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random conjugate-orthonormal matrix          *
* if rows <= cols, Q*QH = I_rows                                *
* if rows >= cols, QH*Q = I_cols                                *
\***************************************************************/
{
  if (rows > 0 && cols > 0)
  { // make sure that A has both rows & columns
    if (cols == rows)
    { // generate a cols x rows (square) conjugate-orthonormal matrix
      make_square_matrix_orth_mp(A, cols, prec); // cols == rows
    }
    else if (cols < rows)
    { // generate a rows x rows conjugate-orthonormal matrix and simply ignore the rest of the matrix
      make_square_matrix_orth_mp(A, rows, prec);
      A->cols = cols;
    }
    else // cols > rows
    { // generate a cols x cols conjugate-orthonormal matrix and shrink to cols x rows and then conjugate-transpose to rows x cols
      make_square_matrix_orth_mp(A, cols, prec);
      A->cols = rows;
      transpose_mp(A, A);
    }
  }
  else
  { // setup A
    change_size_mat_mp(A, rows, cols);
    A->rows = rows;
    A->cols = cols;
  }

  return;
}

void make_matrix_random_real_mp(mat_mp A, int rows, int cols, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random orthonormal real matrix               *
* if rows <= cols, Q*QT = I_rows                                *
* if rows >= cols, QT*Q = I_cols                                *
\***************************************************************/
{
  if (rows > 0 && cols > 0)
  { // make sure that A has both rows & columns
    if (cols == rows)
    { // generate a cols x rows (square) conjugate-orthonormal matrix
      make_square_matrix_orth_real_mp(A, cols, prec); // cols == rows
    }
    else if (cols < rows)
    { // generate a rows x rows conjugate-orthonormal matrix and simply ignore the rest of the matrix
      make_square_matrix_orth_real_mp(A, rows, prec);
      A->cols = cols;
    }
    else // cols > rows
    { // generate a cols x cols conjugate-orthonormal matrix and shrink to cols x rows and then conjugate-transpose to rows x cols
      make_square_matrix_orth_real_mp(A, cols, prec);
      A->cols = rows;
      transpose_mp(A, A);
    }
  }
  else
  { // setup A
    change_size_mat_mp(A, rows, cols);
    A->rows = rows;
    A->cols = cols;
  }

  return;
}

void make_square_matrix_orth_mp(mat_mp A, int size, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random square conjugate-orthonormal matrix   *
* using Stewart's method: A is product of Householder matrices  *
\***************************************************************/
{
  int i, j, k, curr_prec = mpf_get_default_prec();
  mpf_t norm, temp;
  comp_mp gamma;
  vec_mp tempVec1, tempVec2;

  // increase precision
  initMP(prec);

  // intialize
  mpf_init(norm);
  mpf_init(temp);
  init_mp(gamma);
  init_vec_mp(tempVec1, size);
  init_vec_mp(tempVec2, size);

  // initialize to the size x size 'gamma * Id'
  increase_size_mat_mp(A, size, size);
  A->rows = A->cols = size;
  // find gamma and normalize it
  get_comp_rand_mp(gamma);
  mpf_abs_mp(norm, gamma);
  mpf_div(gamma->r, gamma->r, norm);
  mpf_div(gamma->i, gamma->i, norm);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (i == j)
      {
        set_mp(&A->entry[i][j], gamma);
      }
      else
      {
        set_zero_mp(&A->entry[i][j]);
      }

  for (j = 0; j < size - 1; j++)
  { // generate a random vector of unit length in 2-norm (u[j:size])
    mpf_set_ui(norm, 0);
    for (k = j; k < size; k++)
    {
      get_comp_rand_mp(&tempVec1->coord[k]);
      norm_sqr_mp(temp, &tempVec1->coord[k]);
      mpf_add(norm, norm, temp);
    }
    // compute the norm and normalize
    mpf_ui_div(norm, 1, norm);
    mpf_sqrt(norm, norm); // d = 1/sqrt(norm)
    for (k = j; k < size; k++)
      mul_rmpf_mp(&tempVec1->coord[k], &tempVec1->coord[k], norm);

    // compute x = A*u, taking into account sizes
    for (i = 0; i < size; i++)
    {
      set_zero_mp(&tempVec2->coord[i]);
      for (k = j; k < size; k++)
      {
        sum_mul_mp(&tempVec2->coord[i], &A->entry[i][k], &tempVec1->coord[k]);
      }
    }

    // then update A := A - 2*(A*u)*(uH)
    for (i = 0; i < size; i++)
      for (k = j; k < size; k++)
      {
        householder_mult_right_mp(&A->entry[i][k], &A->entry[i][k], &tempVec2->coord[i], &tempVec1->coord[k]);
      }
  }

  // restore precision
  initMP(curr_prec);

  mpf_clear(norm);
  mpf_clear(temp);
  clear_mp(gamma);
  clear_vec_mp(tempVec1);
  clear_vec_mp(tempVec2);

  return;
}

void make_square_matrix_orth_real_mp(mat_mp A, int size, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: creates a random square orthonormal real matrix        *
* using Stewart's method: A is product of Householder matrices  *
\***************************************************************/
{
  int i, j, k, curr_prec = mpf_get_default_prec();
  mpf_t norm, temp;
  comp_mp gamma;
  vec_mp tempVec1, tempVec2;

  // increase precision
  initMP(prec);

  // intialize
  mpf_init(norm);
  mpf_init(temp);
  init_mp(gamma);
  init_vec_mp(tempVec1, size);
  init_vec_mp(tempVec2, size);

  // initialize to the size x size 'gamma * Id'
  increase_size_mat_mp(A, size, size);
  A->rows = A->cols = size;
  // set gamma to either 1 or -1 randomly
  set_one_mp(gamma);
  if (rand() % 2)
    neg_mp(gamma, gamma);
  for (i = 0; i < size; i++)
    for (j = 0; j < size; j++)
      if (i == j)
      {
        set_mp(&A->entry[i][j], gamma);
      }
      else
      {
        set_zero_mp(&A->entry[i][j]);
      }

  for (j = 0; j < size - 1; j++)
  { // generate a random vector of unit length in 2-norm (u[j:size])
    mpf_set_ui(norm, 0);
    for (k = j; k < size; k++)
    {
      get_comp_rand_real_mp(&tempVec1->coord[k]);
      norm_sqr_mp(temp, &tempVec1->coord[k]);
      mpf_add(norm, norm, temp);
    }
    // compute the norm and normalize
    mpf_ui_div(norm, 1, norm);
    mpf_sqrt(norm, norm); // d = 1/sqrt(norm)
    for (k = j; k < size; k++)
      mul_rmpf_mp(&tempVec1->coord[k], &tempVec1->coord[k], norm);

    // compute x = A*u, taking into account sizes
    for (i = 0; i < size; i++)
    {
      set_zero_mp(&tempVec2->coord[i]);
      for (k = j; k < size; k++)
      {
        sum_mul_mp(&tempVec2->coord[i], &A->entry[i][k], &tempVec1->coord[k]);
      }
    }

    // then update A := A - 2*(A*u)*(uH)
    for (i = 0; i < size; i++)
      for (k = j; k < size; k++)
      {
        householder_mult_right_mp(&A->entry[i][k], &A->entry[i][k], &tempVec2->coord[i], &tempVec1->coord[k]);
      }
  }

  // restore precision
  initMP(curr_prec);

  mpf_clear(norm);
  mpf_clear(temp);
  clear_mp(gamma);
  clear_vec_mp(tempVec1);
  clear_vec_mp(tempVec2);

  return;
}

void make_matrix_ID_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int rows, int cols, int prec, int max_prec, int init_mp, int init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;

  if (init_mp)
    init_mat_mp2(A_mp, rows, cols, prec);

  if (init_rat)
  {
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_init(A_rat[i][j][0]);
        mpq_init(A_rat[i][j][1]);
      }
  }

  // setup size 
  change_size_mat_d(A_d, rows, cols);
  change_size_mat_mp(A_mp, rows, cols);
  A_d->rows = A_mp->rows = rows;
  A_d->cols = A_mp->cols = cols;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    { // imag = 0
      mpq_set_str(A_rat[i][j][1], "0", 10);

      if (i == j)
      {
        mpq_set_str(A_rat[i][j][0], "1", 10);
      }
      else
      {
        mpq_set_str(A_rat[i][j][0], "0", 10);
      }

      mpf_set_q(A_mp->entry[i][j].r, A_rat[i][j][0]);
      mpf_set_q(A_mp->entry[i][j].i, A_rat[i][j][1]);
      A_d->entry[i][j].r = mpq_get_d(A_rat[i][j][0]);
      A_d->entry[i][j].i = mpq_get_d(A_rat[i][j][1]);
    }

  return;
}

void make_matrix_ID_mp(mat_mp A, int rows, int cols)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j;
 
  // make A the same size
  change_size_mat_mp(A, rows, cols);
  A->rows = rows; 
  A->cols = cols;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      mpf_set_ui(A->entry[i][j].r, i == j);
      mpf_set_ui(A->entry[i][j].i, 0);
    }

  return;
}

void make_elem_mp(mat_mp A, int num1, int num2)  /* Make A into an ID matrix of its present dimensions, but then swap rows num1 and num2 and cols num1 and num2 to create an elementary "row swap" or "column swap" matrix. */
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  make_matrix_ID_mp(A, A->rows, A->cols);
  set_double_mp2(&A->entry[num1][num1], 0.0, 0.0);
  set_double_mp2(&A->entry[num2][num2], 0.0, 0.0);
  set_double_mp2(&A->entry[num1][num2], 1.0, 0.0);
  set_double_mp2(&A->entry[num2][num1], 1.0, 0.0);

  return;
}

void mat_mul_mp(mat_mp Res, mat_mp A, mat_mp B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, j, k, rows = A->rows, inner = A->cols, cols = B->cols;

  if (A->cols != B->rows)
    printf("WARNING:  Attempting to multiply matrices with dimensions that do not match!\n");

  if (Res != A)
  {
    if (Res != B) // do multiplication in place
    { // setup Res
      change_size_mat_mp(Res, rows, cols);

      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_mp(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_mp(&Res->entry[i][j], &A->entry[i][k], &B->entry[k][j]); // Res_i,j += A_i,k * B_k,j
          }
        }
    }
    else // Res == B 
    { // copy B to tempMat
      mat_mp tempMat;
      init_mat_mp(tempMat, B->rows, B->cols);
      mat_cp_mp(tempMat, B);

      // setup Res
      change_size_mat_mp(Res, rows, cols);

      // multiply Res = A * tempMat
      for (j = 0; j < cols; j++)
        for (i = 0; i < rows; i++)
        {
          set_zero_mp(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_mp(&Res->entry[i][j], &A->entry[i][k], &tempMat->entry[k][j]); // Res_i,j += A_i,k * tempMat_k,j
          }
        }

      // clear tempMat
      clear_mat_mp(tempMat);
    }
  }
  else // Res == A 
  {
    if (Res != B)
    { // copy A to tempMat
      mat_mp tempMat; 
      init_mat_mp(tempMat, A->rows, A->cols);
      mat_cp_mp(tempMat, A);

      // setup Res
      change_size_mat_mp(Res, rows, cols);

      // multiply Res = tempMat * B
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_mp(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_mp(&Res->entry[i][j], &tempMat->entry[i][k], &B->entry[k][j]); // Res_i,j += tempMat_i,k * B_k,j
          }
        }

      // clear tempMat
      clear_mat_mp(tempMat);
    }
    else // Res == A == B
    { // copy A == B to tempMat
      mat_mp tempMat;
      init_mat_mp(tempMat, A->rows, A->cols);
      mat_cp_mp(tempMat, A);  // A == B

      // setup Res
      change_size_mat_mp(Res, rows, cols);

      // multiply Res = tempMat * tempMat
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          set_zero_mp(&Res->entry[i][j]);
          for (k = 0; k < inner; k++)
          {
            sum_mul_mp(&Res->entry[i][j], &tempMat->entry[i][k], &tempMat->entry[k][j]); // Res_i,j += tempMat_i,k * tempMat_k,j
          }
        }

      // clear tempMat
      clear_mat_mp(tempMat);
    }
  }
  // set the size
  Res->rows = rows;
  Res->cols = cols;

  return;
}

void mat_mul_rat(mpq_t ***Res, mpq_t ***A, mpq_t ***B, int A_rows, int A_cols, int B_rows, int B_cols, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: multiplies two 'matrices' of rational numbers          *
\***************************************************************/
{
  int i, j, curr_prec = mpf_get_default_prec();
  mat_mp Res_mat, A_mat, B_mat;

  initMP(max_prec);

  // initialize
  init_mat_mp2(Res_mat, A_rows, B_cols, max_prec);
  init_mat_mp2(A_mat, A_rows, A_cols, max_prec);
  init_mat_mp2(B_mat, B_rows, B_cols, max_prec);
  // setup sizes
  A_mat->rows = A_rows;
  A_mat->cols = A_cols;
  B_mat->rows = B_rows;
  B_mat->cols = B_cols;

  // setup A_mat
  for (i = 0; i < A_rows; i++)
    for (j = 0; j < A_cols; j++)
    {
      mpf_set_q(A_mat->entry[i][j].r, A[i][j][0]);
      mpf_set_q(A_mat->entry[i][j].i, A[i][j][1]);
    }

  // setup B_mat
  for (i = 0; i < B_rows; i++)
    for (j = 0; j < B_cols; j++)
    {
      mpf_set_q(B_mat->entry[i][j].r, B[i][j][0]);
      mpf_set_q(B_mat->entry[i][j].i, B[i][j][1]);
    }

  // Res_mat = A_mat * B_mat
  mat_mul_mp(Res_mat, A_mat, B_mat);

  // setup Res
  for (i = 0; i < A_rows; i++)
    for (j = 0; j < B_cols; j++)
    {
      mpf_t_to_rat(Res[i][j][0], Res_mat->entry[i][j].r);
      mpf_t_to_rat(Res[i][j][1], Res_mat->entry[i][j].i);
    }

  // clear
  clear_mat_mp(Res_mat);
  clear_mat_mp(A_mat);
  clear_mat_mp(B_mat);

  // set precision back
  initMP(curr_prec);

  return;
}

void vec_mat_mul_d(vec_d Res, vec_d b, mat_d A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: vector * matrix                                        *
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;

  if (rows != b->size)
  {
    printf("WARNING:  Attempting to multiply a vector (%d) with a matrix (%d x %d) in which the dimensions do not match!\n", b->size, rows, cols);
  }

  if (Res != b) // do the multiplication in place
  { // setup Res
    change_size_vec_d(Res, cols);

    for (j = 0; j < cols; j++)
    {
      set_zero_d(&Res->coord[j]);
      for (i = 0; i < rows; i++)
      {
        sum_mul_d(&Res->coord[j], &b->coord[i], &A->entry[i][j]);
      }
    }
  }
  else // Res == b
  { // need to use a temporary vector
    vec_d v;
    // copy b to v
    init_vec_d(v, b->size);
    vec_cp_d(v, b);

    // setup Res
    change_size_vec_d(Res, cols);

    for (j = 0; j < cols; j++)
    {
      set_zero_d(&Res->coord[j]);
      for (i = 0; i < rows; i++)
      {
        sum_mul_d(&Res->coord[j], &v->coord[i], &A->entry[i][j]);
      }
    }
    clear_vec_d(v);
  }
  // set the size
  Res->size = cols;

  return;
}

void vec_mat_mul_mp(vec_mp Res, vec_mp b, mat_mp A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: vector * matrix                                        *
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;

  if (rows != b->size)
  {
    printf("WARNING:  Attempting to multiply a vector (%d) with a matrix (%d x %d) in which the dimensions do not match!\n", b->size, rows, cols);
  }

  if (Res != b) // do the multiplication in place
  { // setup Res
    change_size_vec_mp(Res, cols);

    for (j = 0; j < cols; j++)
    {
      set_zero_mp(&Res->coord[j]);
      for (i = 0; i < rows; i++)
      {
        sum_mul_mp(&Res->coord[j], &b->coord[i], &A->entry[i][j]);
      }
    }
  }
  else // Res == b
  { // need to use a temporary vector
    vec_mp v;
    // copy b to v
    init_vec_mp(v, b->size);
    vec_cp_mp(v, b);

    // setup Res
    change_size_vec_mp(Res, cols);
 
    for (j = 0; j < cols; j++)
    {
      set_zero_mp(&Res->coord[j]);
      for (i = 0; i < rows; i++)
      {
        sum_mul_mp(&Res->coord[j], &v->coord[i], &A->entry[i][j]);
      }
    }
    clear_vec_mp(v);
  }
  // set the size
  Res->size = cols;

  return;
}

void vec_mat_mul_rat(mpq_t **Res, mpq_t **b, mpq_t ***A, int b_size, int A_rows, int A_cols, int init_Res)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: multiplies 'vector' and a 'matrix' of rational numbers *
*   used on in setupPatch_d_to_mp                               *
\***************************************************************/
{
  int i, j;
  mpq_t *tempMPQ = (mpq_t *)bmalloc(2 * sizeof(mpq_t));

  for (i = 0; i < 2; i++)
    mpq_init(tempMPQ[i]);

  for (j = 0; j < A_cols; j++)
  {
    if (init_Res)
    {
      mpq_init(Res[j][0]);
      mpq_init(Res[j][1]);
    }
    set_zero_rat(Res[j]);
    for (i = 0; i < A_rows; i++)
    {
      mpq_mul(tempMPQ[0], b[i][0], A[i][j][0]);
      mpq_mul(tempMPQ[1], b[i][1], A[i][j][1]);
      mpq_sub(tempMPQ[0], tempMPQ[0], tempMPQ[1]);
      mpq_add(Res[j][0], Res[j][0], tempMPQ[0]);

      mpq_mul(tempMPQ[0], b[i][0], A[i][j][1]);
      mpq_mul(tempMPQ[1], b[i][1], A[i][j][0]);
      mpq_add(tempMPQ[0], tempMPQ[0], tempMPQ[1]);
      mpq_add(Res[j][1], Res[j][1], tempMPQ[0]);
    }
  }

  for (i = 1; i >= 0; i--)
    mpq_clear(tempMPQ[i]);

  return;
}

int get_rand_int(int *x)
/***************************************************************\
* USAGE:  Sets x to a randomly chosen integer                   
* ARGUMENTS:  x (the integer)                                       
* RETURN VALUES:  the number of digits in x                             
* NOTES:  Used primarily in parse_and_hom.y - get_rand() is a  
*         more common choice
\***************************************************************/
{
  while (1)  //We cycle through randoms until we get one with between 6 and 11 digits.
  {
    *x=random();
    if ((*x>1e5) && (*x<1e6))
      return 6;
    if ((*x>1e6) && (*x<1e7))
      return 7;
    if ((*x>1e7) && (*x<1e8))
      return 8;
    if ((*x>1e8) && (*x<1e9))
      return 9;
//    if ((*x>1e9) && (*x<1e10))
//      return 10;
  }
}

double get_comp_rand_d(comp_d x)
/***************************************************************\
* USAGE: obtain a complex random number in double precision     *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the complex random number           *
* NOTES:                                                        *
\***************************************************************/
{
  double tempMod;

  // random numbers in [-1, 1]
  set_double_d(x, 2 * (rand() / (RAND_MAX + 1.0) - 0.5), 2 * (rand() / (RAND_MAX + 1.0) - 0.5));

  // find the sqrt of the modulus - the modulus of the complex random number that is returned
  tempMod = sqrt(d_abs_d(x));

  // divide by the sqrt of the modulus to keep it sufficiently away from zero
  x->r /= tempMod;
  x->i /= tempMod;

  return tempMod;
}

double get_comp_rand_real_d(comp_d x)
/***************************************************************\
* USAGE: obtain a random real number in double precision        *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the random number                   *
* NOTES:                                                        *
\***************************************************************/
{
  double tempMod;

  // random numbers in [-1, 1]
  set_double_d(x, 2 * (rand() / (RAND_MAX + 1.0) - 0.5), 0);

  // find the sqrt of the modulus - the modulus of the random number that is returned
  tempMod = sqrt(fabs(x->r));

  // divide by the sqrt of the modulus to keep it sufficiently away from zero
  x->r /= tempMod;

  return tempMod;
}

void get_comp_rand_mp2(mpf_t mod, comp_mp x)
/***************************************************************\
* USAGE: obtain a complex random number in multi precision      *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the complex random number           *
* NOTES:                                                        *
\***************************************************************/
{
  int base = 10, prec = mpf_get_prec(x->r);
  int num_digits = (int) ceil(log10(2) * prec + 1);
  char *str = (char *)bmalloc((2 * num_digits + 4) * sizeof(char));
  mpq_t frac;

  mpq_init(frac);

  // setup real part
  create_random_number_str(str, num_digits);
  mpq_set_str(frac, str, base);
  mpq_canonicalize(frac);
  mpf_set_q(x->r, frac);

  // setup imag part
  create_random_number_str(str, num_digits);
  mpq_set_str(frac, str, base);
  mpq_canonicalize(frac);
  mpf_set_q(x->i, frac);

  // find the sqrt of the modulus - the modulus of the complex random number that is returned
  mpf_abs_mp(mod, x);
  mpf_sqrt(mod, mod);

  // divide by the sqrt of the modulus to keep it sufficiently away from zero
  mpf_div(x->r, x->r, mod);
  mpf_div(x->i, x->i, mod);

  // clear memory
  free(str);
  mpq_clear(frac);

  return;
}

void get_comp_rand_real_mp2(mpf_t mod, comp_mp x)
/***************************************************************\
* USAGE: obtain a random real number in multi precision         *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the random number                   *
* NOTES:                                                        *
\***************************************************************/
{
  int base = 10, prec = mpf_get_prec(x->r);
  int num_digits = (int) ceil(log10(2) * prec + 1);
  char *str = (char *)bmalloc((2 * num_digits + 4) * sizeof(char));
  mpq_t frac;

  mpq_init(frac);

  // setup real part
  create_random_number_str(str, num_digits);
  mpq_set_str(frac, str, base);
  mpq_canonicalize(frac);
  mpf_set_q(x->r, frac);

  // setup imag part
  mpf_set_ui(x->i, 0);

  // find the sqrt of the modulus - the modulus of the random number that is returned
  mpf_abs(mod, x->r);
  mpf_sqrt(mod, mod);

  // divide by the sqrt of the modulus to keep it sufficiently away from zero
  mpf_div(x->r, x->r, mod);

  // clear memory
  free(str);
  mpq_clear(frac);

  return;
}

double get_comp_rand_mp(comp_mp x)
/***************************************************************\
* USAGE: obtain a complex random number in multi precision      *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the complex random number           *
* NOTES:                                                        *
\***************************************************************/
{
  mpf_t norm_mp;
  double norm_d;

  mpf_init(norm_mp);

  get_comp_rand_mp2(norm_mp, x);

  norm_d = mpf_get_d(norm_mp);

  mpf_clear(norm_mp);

  return norm_d;
}

double get_comp_rand_real_mp(comp_mp x)
/***************************************************************\
* USAGE: obtain a random real number in multi precision         *
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the random number                   *
* NOTES:                                                        *
\***************************************************************/
{
  mpf_t norm_mp;
  double norm_d;

  mpf_init(norm_mp);

  get_comp_rand_real_mp2(norm_mp, x);

  norm_d = mpf_get_d(norm_mp);

  mpf_clear(norm_mp);

  return norm_d;
}

void create_random_number_str(char *str, int num_digits)
/***************************************************************\
* USAGE: creates a random number in [-1,1] using num_digits for *
*    its numerator                                              *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, tempInt, counter = 0;

  // random sign
  tempInt = rand() % 2;
  if (tempInt)
  {
    str[counter] = '-';
    counter++;
  }

  for (i = 0; i < num_digits; i++)
  {
    tempInt = rand() % 10;
    str[counter] = 48 + tempInt; // ASCII for the digits
    counter++;
  }
  str[counter] = '/';
  counter++;
  str[counter] = '1';
  counter++;
  for (i = 0; i < num_digits; i++)
  {
    str[counter] = '0';
    counter++;
  }
  str[counter] = '\0';

  return;
}

double get_comp_rand_rat(comp_d x_d, comp_mp x_mp, mpq_t *x_rat, int curr_prec, int max_prec, int need_to_init_mp, int need_to_init_rat)
/***************************************************************\
* USAGE: obtain a complex random number in double, MP & rational*
* ARGUMENTS:                                                    *
* RETURN VALUES: modulus of the complex random number           *
* NOTES:                                                        *
\***************************************************************/
{
  int i, num_digits = (int) ceil(log10(2) * max_prec + 1), base = 10;
  char *str = (char *)bmalloc((2 * num_digits + 4) * sizeof(char));
  double rV;
  mpf_t mod;
  comp_mp tempComp;

  mpf_init2(mod, max_prec);
  init_mp2(tempComp, max_prec);

  // do real & imag parts
  for (i = 0; i < 2; i++)
  { // create a random number string
    create_random_number_str(str, num_digits);

    // set to a rational number
    if (need_to_init_rat)
      mpq_init(x_rat[i]);
    mpq_set_str(x_rat[i], str, base);
    mpq_canonicalize(x_rat[i]);
  }

  // approximate to floating point precision in max_prec
  mpf_set_q(tempComp->r, x_rat[0]);
  mpf_set_q(tempComp->i, x_rat[1]);

  // find the sqrt of the modulus
  mpf_abs_mp(mod, tempComp);
  mpf_sqrt(mod, mod);

  // find double approximation of mod - this is returned
  rV = mpf_get_d(mod);

  // divide by the sqrt of the modulus to keep it sufficiently away from zero
  mpf_div(tempComp->r, tempComp->r, mod);
  mpf_div(tempComp->i, tempComp->i, mod);

  // convert back to rat
  mpf_t_to_rat(x_rat[0], tempComp->r);
  mpf_t_to_rat(x_rat[1], tempComp->i);

  // setup x_mp & x_d
  if (need_to_init_mp)
  { // initialize to curr_prec
    init_mp2(x_mp, curr_prec);
  }
  else
  { // set to curr_prec
    setprec_mp(x_mp, curr_prec);
  }
  mpf_set_q(x_mp->r, x_rat[0]);
  mpf_set_q(x_mp->i, x_rat[1]);
  x_d->r = mpq_get_d(x_rat[0]);
  x_d->i = mpq_get_d(x_rat[1]);

  // clear memory
  free(str);
  mpf_clear(mod);
  clear_mp(tempComp);

  return rV;
}

void convert_point_data_d_to_mp(point_data_mp *dataMP, point_data_d *data)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = data->point->size;

  dataMP->cycle_num = data->cycle_num;
  change_size_point_mp(dataMP->point, size);
  dataMP->point->size = size;

  // copy point 
  for (i = 0; i < size; i++)
    d_to_mp(&dataMP->point->coord[i], &data->point->coord[i]);

  // copy time
  d_to_mp2(dataMP->time, data->time); // using a rational conversion for the time - more accurate but slower

  return;
}

void convert_point_data_mp_to_d(point_data_d *data, point_data_mp *dataMP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i, size = dataMP->point->size;

  data->cycle_num = dataMP->cycle_num;
  change_size_point_d(data->point, size);
  data->point->size = size;

  // copy point
  for (i = 0; i < size; i++)
    mp_to_d(&data->point->coord[i], &dataMP->point->coord[i]);

  // copy time
  mp_to_d(data->time, dataMP->time);

  return;
}

void point_data_cp_d(point_data_d *dest, point_data_d *src)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  point_cp_d(dest->point, src->point);
  //point_cp_d(dest->param, src->param);  //Commented to be consistent with the mp version.
  set_d(dest->time, src->time); 
  dest->cycle_num = src->cycle_num; 
  return;
}

void point_data_cp_mp(point_data_mp *dest, point_data_mp *src)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  point_cp_mp(dest->point, src->point);
  //point_cp_mp(dest->param, src->param);  //Commented because this somehow causes the power series endgame to seg fault -> no big loss.
  set_mp(dest->time, src->time);
  dest->cycle_num = src->cycle_num;
  return;
}

char* mpf_to_str(mpf_t MPF, int base)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: returns a pointer to the str representing MPF          *
*   if MPF is NaN, return str is "0"                            *
\***************************************************************/
{
  int size;
  char *str = NULL, *tmpStr = NULL;
  mp_exp_t exp;

  if (mpfr_number_p(MPF))
  { // find the mantissa
    tmpStr = mpf_get_str(NULL, &exp, base, 0, MPF);

    // find the size needed
    size = 1 + snprintf(NULL, 0, "0.%se%ld", tmpStr, exp); // +1 for '\0'
    // allocate str
    str = (char *)bmalloc(size * sizeof(char));

    // setup str
    if (mpf_sgn(MPF) >= 0)
    { 
      sprintf(str, "0.%se%ld", tmpStr, exp);
    }
    else
    { 
      sprintf(str, "-0.%se%ld", &tmpStr[1], exp); 
    }

    free(tmpStr); 
  }
  else
  {
    str = (char *)bmalloc(2 * sizeof(char));
    strcpy(str, "0");
  }
  
 return str;
}
 
size_t outStr_to_frac_size(size_t *numer_zeros, size_t *denom_zeros, char *strMan, long int exp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: total size & number of zeros in num & denom    *
* NOTES: finds the size of a fraction string created by a       *
* number with strMan as the mantissa and exp as the exponent    *
*   strMan either starts with '-' or a number and the rest of   *
*    the characters must be a number                            *
\***************************************************************/
{
  size_t retVal, isNeg, size_numer, size_denom;

  // find the initial size of the numerator & denominator
  size_numer = strlen(strMan);
  size_denom = 1;  // for '1' in denom

  // make sure the strMan has length
  if (size_numer == 0)
  {
    printf("ERROR: The mantissa has no digits in outStr_to_frac_size!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // determine if negative
  if (strMan[0] == '-')
  {
    isNeg = 1;
    size_numer--; // one of the 'digits' is actually a negative sign
  }
  else
    isNeg = 0;

  // determine how the exponent fits into this
  if (exp < 0) // this is needed since size_t is unsigned!!!!!
  { // negative exponent creates digits in the denominator
    *numer_zeros = 0;
    *denom_zeros = size_numer - exp;
  }
  else if (size_numer > exp)
  { // the numerator has more digits than the exponent so just subtract them to get the number of zeros needed
    *numer_zeros = 0;
    *denom_zeros = size_numer - exp;
  }
  else
  { // the numerator needs to be padded with (exp - size_numer) zeros
    *numer_zeros = exp - size_numer;
    *denom_zeros = 0;
  }

  size_numer += (*numer_zeros);
  size_denom += (*denom_zeros);
  
  // find the total number = negative sign + numerator + '/' + denominator + '\0' 
  retVal = isNeg + size_numer + 1 + size_denom + 1;

  return retVal;
}

void outStr_to_frac(char *strFrac, char *strMan, size_t numer_zeros, long int exp, size_t denom_zeros)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sets strFrac to a fraction string created by a         *
* number with strMan as the mantissa and exp as the exponent    *
*    strFrac should be large enough to hold this fraction       *
*      use outStr_to_frac_size to find the correct sizes!!!!    *
\***************************************************************/
{
  size_t i;

  // copy strMan to strFrac
  strcpy(strFrac, strMan);

  // add the correct number of zeros in the numerator
  for (i = 0; i < numer_zeros; i++)
  {
    strcat(strFrac, "0");
  }

  // put in "/1"
  strcat(strFrac, "/1");

  // add the correct number of zeros in the denominator
  for (i = 0; i < denom_zeros; i++)
  {
    strcat(strFrac, "0");
  }

  return;
}

void mpf_t_to_rat(mpq_t RAT, mpf_t MP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts a comp_mp to mpq_t                            *
* sets RAT to 0 if MP is not a valid number                     *
\***************************************************************/
{
  char *strFrac = NULL;

  if (mpfr_number_p(MP))
  { 
    char *strOut; // string for holding the mantissa for MP
    long int exp; // hold exponent
    size_t size, numer_zeros, denom_zeros; 

    // find the number to its full precision
    strOut = mpf_get_str(NULL, &exp, 10, 0, MP);

    // find the sizes
    size = outStr_to_frac_size(&numer_zeros, &denom_zeros, strOut, exp);
    // allocate 'size' for strFrac
    strFrac = (char *)bmalloc(size * sizeof(char));
    // setup strFrac
    outStr_to_frac(strFrac, strOut, numer_zeros, exp, denom_zeros);

    // free strOut
    mpfr_free_str(strOut);  
  }
  else
  { // set to 0
    strFrac = (char *)bmalloc(2 * sizeof(char));
    strcpy(strFrac, "0");
  }

  // setup RAT
  mpq_set_str(RAT, strFrac, 10); // set using frac in base 10
  mpq_canonicalize(RAT);  // dump common factors

  // release memory
  free(strFrac);

  return;
}

void change_prec_mp(comp_mp x, int new_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  if (new_prec != mpf_get_prec(x->r))
  {
    comp_mp tmp;

    init_mp2(tmp, mpf_get_prec(x->r));

    set_mp(tmp, x);
    setprec_mp(x, new_prec);
    set_mp(x, tmp);

    clear_mp(tmp);
  }

  return;
}

void change_prec_mp2(comp_mp x, int new_prec)  
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: uses rational number so answer is more exact - slower  *
\***************************************************************/
{
  if (new_prec != mpf_get_prec(x->r))
  {
    mpq_t rr, ri;

    mpq_init(rr);
    mpq_init(ri);

    // convert x to rational expressions
    mpf_t_to_rat(rr, x->r);
    mpf_t_to_rat(ri, x->i);

    // change precision on x
    setprec_mp(x, new_prec);

    // set back to x
    mpf_set_q(x->r, rr);
    mpf_set_q(x->i, ri);

    // free the memory
    mpq_clear(rr); 
    mpq_clear(ri);
  }

  return;
}

void mypause()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  printf("Pausing...");
  scanf("%d", &i);
  printf("Proceeding.\n");

  return;
}

void *bmalloc(size_t size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does malloc with error checking                        *
\***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate
    return NULL;
  }
  else
  { // try to allocate memory
    void *x = malloc(size);
    if (x == NULL)
    {
      printf("ERROR: malloc was unable to allocate memory (%d)!\n", (int) size);
      bexit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
}

void *bcalloc(size_t num, size_t size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does calloc with error checking                        *
\***************************************************************/
{
  if (num <= 0 || size <= 0)
  { // nothing to allocate
    return NULL;
  }
  else
  { // try to allocate memory
    void *x = calloc(num, size);
    if (x == NULL)
    {
      printf("ERROR: calloc was unable to allocate memory!\n");
      bexit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
}

void *brealloc(void *ptr, size_t size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does realloc with error checking                       *
\***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate - free memory and return NULL
    free(ptr);
    ptr = NULL;
  }
  else
  { // try to reallocate memory
    ptr = realloc(ptr, size);
    if (ptr == NULL)
    {
      printf("ERROR: realloc was unable to allocate memory!\n");
      bexit(ERROR_MEMORY_ALLOCATION);
    }
  }
  return ptr;
}

void bexit(int errorCode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: exits Bertini - either standard or using MPI           *
\***************************************************************/
{
  if (errorCode == 0)
    errorCode = ERROR_OTHER;

  printf("%s\n", BERTINI_QUIT_MESSAGE);
#ifdef _HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, errorCode);
#else
  exit(errorCode);
#endif
}

double amp_criterion_A(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Calculates criterion A from AMP paper                  *
\***************************************************************/
{
  return safety_digits + log10(norm_J_inv * eps * (norm_J + Phi));
}

double amp_criterion_B(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi, double tol, double residual, int maxIts, int it_number)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Calculates criterion B from AMP paper                  *
\***************************************************************/
{
  if (maxIts - it_number == 1) // make sure the denominator will not be 0
    return 0;
  else
    return safety_digits + log10(norm_J_inv * ((2 + eps) * norm_J + eps * Phi) + 1) + (-log10(tol) + log10(residual)) / (maxIts - it_number - 1.0);
}

double amp_criterion_B2(int safety_digits, double norm_J, double norm_J_inv, double eps, double Phi, double tol, double proportion, int maxIts)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Calculates criterion B from AMP2 paper                 *
\***************************************************************/
{
  return safety_digits + log10(norm_J_inv * ((2 + eps) * norm_J + eps * Phi) + 1) + (-log10(tol) + log10(proportion)) / maxIts;
}

double amp_criterion_C(int safety_digits, double norm_J_inv, double Psi, double tol, double size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Calculates criterion C from AMP paper                  *
\***************************************************************/
{
  return safety_digits - log10(tol) + log10(norm_J_inv * Psi + size);
}

void d_to_mp2(comp_mp m, comp_d d)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts d into m using rational representation        *
\***************************************************************/
{
  int digits_after_decimal = 15; // number of digits in double precision printed after the decimal
  int mantissa_string_size = 3 + digits_after_decimal; // sign + 1st digit + digits_after_decimal + '\0'
  size_t size, numer_zeros, denom_zeros;
  long int exp;
  char *strFrac = NULL, *strOut = NULL;
  mpq_t rat;

  // initialize rat and allocate for strOut
  mpq_init(rat);
  strOut = (char *)bmalloc(mantissa_string_size * sizeof(char));

  // do the real part

  // convert to mantissa and exp
  d_get_str(strOut, &exp, digits_after_decimal, d->r);

  // find the sizes
  size = outStr_to_frac_size(&numer_zeros, &denom_zeros, strOut, exp);
  // allocate memory
  strFrac = (char *)bmalloc(size * sizeof(char));
  // setup strFrac
  outStr_to_frac(strFrac, strOut, numer_zeros, exp, denom_zeros);

  // set rat using strFrac
  mpq_set_str(rat, strFrac, 10); // in base 10
  mpq_canonicalize(rat);  

  // set real value of m using rat
  mpf_set_q(m->r, rat);

  // free strFrac
  free(strFrac);

  // do the imaginary part

  // convert to mantissa and exp
  d_get_str(strOut, &exp, digits_after_decimal, d->i);

  // find the sizes
  size = outStr_to_frac_size(&numer_zeros, &denom_zeros, strOut, exp);
  // allocate memory
  strFrac = (char *)bmalloc(size * sizeof(char));
  // setup strFrac
  outStr_to_frac(strFrac, strOut, numer_zeros, exp, denom_zeros);

  // set rat using strFrac
  mpq_set_str(rat, strFrac, 10); // in base 10
  mpq_canonicalize(rat);

  // set imaginary value of m using rat
  mpf_set_q(m->i, rat);

  // clear rat and free memory
  mpq_clear(rat);
  free(strFrac);
  free(strOut);

  return;
}

void d_get_str(char *strOut, long int *exp, int digits, double D)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts D to strOut and exp like mpf_get_str          *
*   strOut needs to be a max of 18: sign + 16 digits + '\0'     *
\***************************************************************/
{
  size_t size;
  char *formatStr = NULL, *numStr = NULL;

  // error checking
  if (digits <= 0 || digits >= 16)
  { // set digits to 15 - that is, 15 digits after the decimal place for a total of 16 digits displayed
    digits = 15;
  }

  // find the size of the format string
  size = 1 + snprintf(NULL, 0, "%%.%de", digits);
  // allocate size
  formatStr = (char *)bmalloc(size * sizeof(char));
  // setup formatStr
  sprintf(formatStr, "%%.%de", digits);

  // find the size of the number string
  size = 1 + snprintf(NULL, 0, formatStr, D);
  // allocate size
  numStr = (char *)bmalloc(size * sizeof(char));
  // setup numStr
  sprintf(numStr, formatStr, D);

  // numStr is of the form: [-] x.x..digits..x e exp

  // setup strOut

  // check for negative sign
  if (numStr[0] == '-')
  {
    strcpy(strOut, "-");
    size = 1;
  }
  else
  {
    strcpy(strOut, "");
    size = 0;
  }
  // copy over the one digit before the '.'
  strncat(strOut, &numStr[size], 1);

  // skip past that digit and '.'
  size += 2;
  
  // copy over the 'digits' digits
  strncat(strOut, &numStr[size], digits);

  // skip past those digits and 'e'
  size += digits + 1;

  // setup exp - move decimal place to the left and so the +1 is to offset this
  *exp = atol(&numStr[size]) + 1;

  // free the memory
  free(formatStr);
  free(numStr);

  return;
}

void printPatchCoeff_mat(FILE *OUT, int MPType, int printType, mat_d patch_d, mat_mp patch_mp, mpq_t ***patch_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the patch coefficients to OUT so they could be  *
* used again                                                    *
\***************************************************************/
{
  int i, j, rows, cols;

  if (printType == 0)
  { // print using _d
    rows = patch_d->rows;
    cols = patch_d->cols;

    // print the number of rows & cols to OUT
    fprintf(OUT, "%d %d\n", rows, cols);

    // print the coeff to OUT
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        fprintf(OUT, "%.15e %.15e\n", patch_d->entry[i][j].r, patch_d->entry[i][j].i);
      }
  }
  else if (printType == 1)
  { // print using _mp
    rows = patch_mp->rows;
    cols = patch_mp->cols;

    // print the number of rows & cols to OUT
    fprintf(OUT, "%d %d\n", rows, cols);

    // print the coeff to OUT
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpf_out_str(OUT, 10, 0, patch_mp->entry[i][j].r);
        fprintf(OUT, " "); 
        mpf_out_str(OUT, 10, 0, patch_mp->entry[i][j].i);
        fprintf(OUT, "\n");
      }
  }
  else
  { // print using _rat
    rows = patch_mp->rows;
    cols = patch_mp->cols;

    // print the number of rows & cols to OUT
    fprintf(OUT, "%d %d\n", rows, cols);

    // print the coeff to OUT
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        mpq_out_str(OUT, 10, patch_rat[i][j][0]);
        fprintf(OUT, " ");
        mpq_out_str(OUT, 10, patch_rat[i][j][1]);
        fprintf(OUT, "\n");
      }
  }
  fprintf(OUT, "\n");

  return;
}

void printPatchCoeff(FILE *OUT, int MPType, basic_eval_data_d *BED_d, basic_eval_data_mp *BED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the patch coefficients to OUT so they could be  *
* used again                                                    *
\***************************************************************/
{
  if (MPType == 0)
  { // print using double precision 
    printPatchCoeff_mat(OUT, MPType, 0, BED_d->patch.patchCoeff, NULL, NULL);
  }
  else if (MPType == 1 || MPType == 2)
  { // print using rational representation
    printPatchCoeff_mat(OUT, MPType, 2, NULL, BED_mp->patch.patchCoeff, BED_mp->patch.patchCoeff_rat);
  }

  return;
}

void findFunctionResidual_conditionNumber_d(double *func_residual, double *cond_num, point_data_d *PD, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the residual of the function & condition number  *
\***************************************************************/
{
  eval_struct_d e;
  init_eval_struct_d(e, 0, 0, 0);

  // evaluate the homotopy
  eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func);

  // find function residual
  *func_residual = infNormVec_d(e.funcVals);

  // find condition number
  *cond_num = conditionNumber_d(e.Jv);

  // clear e
  clear_eval_struct_d(e);

  return;
}

void findFunctionResidual_d(double *func_residual, point_data_d *PD, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the residual of the function                     *
\***************************************************************/
{
  eval_struct_d e;
  init_eval_struct_d(e, 0, 0, 0);

  // evaluate the homotopy
  eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func);

  // find function residual
  *func_residual = infNormVec_d(e.funcVals);

  // clear e
  clear_eval_struct_d(e);

  return;
}

void findFunctionResidual_conditionNumber_mp(mpf_t func_residual, double *cond_num, point_data_mp *PD, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the residual of the function & condition number  *
\***************************************************************/
{
  eval_struct_mp e;
  init_eval_struct_mp(e, 0, 0, 0);
 
  // evaluate the homotopy
  eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func);

  // find function residual
  infNormVec_mp2(func_residual, e.funcVals);

  // find condition number
  *cond_num = conditionNumber_mp(e.Jv);

  // clear
  clear_eval_struct_mp(e);

  return;
}

void findFunctionResidual_mp(mpf_t func_residual, point_data_mp *PD, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the residual of the function                     *
\***************************************************************/
{
  eval_struct_mp e;
  init_eval_struct_mp(e, 0, 0, 0);

  // evaluate the homotopy
  eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func);

  // find function residual
  infNormVec_mp2(func_residual, e.funcVals);

  // clear
  clear_eval_struct_mp(e);

  return;
}

void readInPatch(FILE *IN, int MPType, int old_MPType, patch_eval_data_d *PED_d, patch_eval_data_mp *PED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads the patch information from IN and stores to PED  *
\***************************************************************/
{
  int i, j, rows, cols;

  // read in the number of rows & cols
  fscanf(IN, "%d %d\n", &rows, &cols);

  if (rows > 0 && cols > 0)
  {
    if (old_MPType == 0)
    { // read in the entries using double precision
  
      mat_d patch_d;
      init_mat_d(patch_d, rows, cols);
      patch_d->rows = rows;
      patch_d->cols = cols;

      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          fscanf(IN, "%lf %lf\n", &patch_d->entry[i][j].r, &patch_d->entry[i][j].i);
        }

      if (MPType == 0)
      { // make sure that these are the correct sizes
        if (PED_d->num_patches != rows || PED_d->patchCoeff->rows != rows || PED_d->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_d to patchCoeff
        mat_cp_d(PED_d->patchCoeff, patch_d);
      }
      else if (MPType == 1) 
      { // make sure that these are the correct sizes
        if (PED_mp->num_patches != rows || PED_mp->patchCoeff->rows != rows || PED_mp->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_d to patchCoeff
        j = (int) mpf_get_prec(PED_mp->patchCoeff->entry[0][0].r);
        mat_d_to_mp_rat(PED_mp->patchCoeff, PED_mp->patchCoeff_rat, patch_d, 16, j, 0, 0);
      }
      else if (MPType == 2)
      { // make sure that these are the correct sizes
        if (PED_mp->num_patches != rows || PED_d->num_patches != rows || PED_mp->patchCoeff->rows != rows || PED_d->patchCoeff->rows != rows || PED_mp->patchCoeff->cols != cols || PED_d->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_d to patchCoeff_mp
        j = (int) mpf_get_prec(PED_mp->patchCoeff->entry[0][0].r);
        mat_d_to_mp_rat(PED_mp->patchCoeff, PED_mp->patchCoeff_rat, patch_d, 16, j, 0, 0);
        // copy patch_d to patchCoeff_d
        mat_cp_d(PED_d->patchCoeff, patch_d);
      }
    }
    else
    { // read in the entries using rational representation
      mpq_t ***patch_rat = NULL;
      init_mat_rat(patch_rat, rows, cols);

      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        { // setup the real part
          mpq_inp_str(patch_rat[i][j][0], IN, 10);
          mpq_canonicalize(patch_rat[i][j][0]);
          // setup the imaginary part
          mpq_inp_str(patch_rat[i][j][1], IN, 10);
          mpq_canonicalize(patch_rat[i][j][1]);
          // scan rest of line
          fscanf(IN, "\n");
        }

      if (MPType == 0)
      { // make sure that these are the correct sizes
        if (PED_d->num_patches != rows || PED_d->patchCoeff->rows != rows || PED_d->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_rat to patchCoeff - already setup
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            PED_d->patchCoeff->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
            PED_d->patchCoeff->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
          }
      }
      else if (MPType == 1)
      { // make sure that these are the correct sizes
        if (PED_mp->num_patches != rows || PED_mp->patchCoeff->rows != rows || PED_mp->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_rat to patchCoeff_rat & patchCoeff - already setup
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpq_set(PED_mp->patchCoeff_rat[i][j][0], patch_rat[i][j][0]);
            mpq_set(PED_mp->patchCoeff_rat[i][j][1], patch_rat[i][j][1]);

            mpf_set_q(PED_mp->patchCoeff->entry[i][j].r, patch_rat[i][j][0]);
            mpf_set_q(PED_mp->patchCoeff->entry[i][j].i, patch_rat[i][j][1]);
          }
      }
      else if (MPType == 2)
      { // make sure that these are the correct sizes
        if (PED_mp->num_patches != rows || PED_d->num_patches != rows || PED_mp->patchCoeff->rows != rows || PED_d->patchCoeff->rows != rows || PED_mp->patchCoeff->cols != cols || PED_d->patchCoeff->cols != cols)
        {
          printf("ERROR: The patch dimensions do not match!\n");
          bexit(ERROR_INVALID_SIZE);
        }

        // copy patch_rat to patchCoeff_rat & patchCoeff_mp & patchCoeff_d - already setup
        for (i = 0; i < rows; i++)
          for (j = 0; j < cols; j++)
          {
            mpq_set(PED_mp->patchCoeff_rat[i][j][0], patch_rat[i][j][0]);
            mpq_set(PED_mp->patchCoeff_rat[i][j][1], patch_rat[i][j][1]);

            mpf_set_q(PED_mp->patchCoeff->entry[i][j].r, patch_rat[i][j][0]);
            mpf_set_q(PED_mp->patchCoeff->entry[i][j].i, patch_rat[i][j][1]);

            PED_d->patchCoeff->entry[i][j].r = mpq_get_d(patch_rat[i][j][0]);
            PED_d->patchCoeff->entry[i][j].i = mpq_get_d(patch_rat[i][j][1]);
          }
      }

      // free the memory
      clear_mat_rat(patch_rat, rows, cols);
    }
  }

  return; 
}

void init_all_mat_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes all mat_mp in an efficient way!!!          *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;

    mat_mp **next = (mat_mp **)bmalloc(num * sizeof(mat_mp *));

    // initialize arg_addr
    va_start(arg_addr, num);

    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (mat_mp *) va_arg(arg_addr, void *);

      (*next[i])->rows = (*next[i])->cols = 0;
    }

    // clear arg_addr
    va_end(arg_addr);

    // initilize the matrices
    for (i = 0; i < num; i++)
      init_mat_mp(*next[i], 0, 0);

    // free the pointers
    free(next);
  }

  return;
}

void init_all_mat_mp2(int prec, int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes all mat_mp in an efficient way!!!          *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;
 
    mat_mp **next = (mat_mp **)bmalloc(num * sizeof(mat_mp *));
 
    // initialize arg_addr
    va_start(arg_addr, num);
 
    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (mat_mp *) va_arg(arg_addr, void *);
 
      (*next[i])->rows = (*next[i])->cols = 0;
    }
 
    // clear arg_addr
    va_end(arg_addr);
 
    // initialize the matrices
    for (i = 0; i < num; i++)
      init_mat_mp2(*next[i], 0, 0, prec);

    // free the pointers
    free(next);
  }
 
  return;
}

void init_all_vec_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes all vec_mp in an efficient way!!!          *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;

    vec_mp **next = (vec_mp **)bmalloc(num * sizeof(vec_mp *));

    // initialize arg_addr
    va_start(arg_addr, num);

    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (vec_mp *) va_arg(arg_addr, void *);

      (*next[i])->size = 0;
    }

    // clear arg_addr
    va_end(arg_addr);

    // initialize the vector
    for (i = 0; i < num; i++)
      init_vec_mp(*next[i], 0);

    // free the pointers
    free(next);
  }

  return;
}

void init_all_vec_mp2(int prec, int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes all vec_mp in an efficient way!!!          *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;

    vec_mp **next = (vec_mp **)bmalloc(num * sizeof(vec_mp *));

    // initialize arg_addr
    va_start(arg_addr, num);

    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (vec_mp *) va_arg(arg_addr, void *);

      (*next[i])->size = 0;
    }

    // clear arg_addr
    va_end(arg_addr);

    // initialize the vector
    for (i = 0; i < num; i++)
      init_vec_mp2(*next[i], 0, prec);

    // free the pointers
    free(next);
  }

  return;
}

void init_all_point_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initializes all point_mp in an efficient way!!!        *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;

    point_mp **next = (point_mp **)bmalloc(num * sizeof(point_mp *));

    // initialize arg_addr
    va_start(arg_addr, num);

    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (point_mp *) va_arg(arg_addr, void *);

      (*next[i])->size = 0;
    }

    // clear arg_addr
    va_end(arg_addr);

    // initialize the point
    for (i = 0; i < num; i++)
      init_vec_mp(*next[i], 0);

    // free the pointers
    free(next);
  }

  return;
}

void clear_all_mat_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears all mat_mp in an efficient way!!!               *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;

    mat_mp **next = (mat_mp **)bmalloc(num * sizeof(mat_mp *));

    // initialize arg_addr
    va_start(arg_addr, num);

    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (mat_mp *) va_arg(arg_addr, void *);

      (*next[i])->rows = (*next[i])->cols = 0;
    }

    // clear arg_addr
    va_end(arg_addr);
 
    // clear the matrices
    for (i = 0; i < num; i++)
      clear_mat_mp(*next[i]);

    // free the pointers
    free(next);
  }

  return;
}

void clear_all_vec_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears all vec_mp in an efficient way!!!               *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;
 
    vec_mp **next = (vec_mp **)bmalloc(num * sizeof(vec_mp *));
 
    // initialize arg_addr
    va_start(arg_addr, num);
 
    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (vec_mp *) va_arg(arg_addr, void *);
 
      (*next[i])->size = 0;
    }
 
    // clear arg_addr
    va_end(arg_addr);
 
    // clear the vector
    for (i = 0; i < num; i++)
      clear_vec_mp(*next[i]);
 
    // free the pointers
    free(next);
  }
 
  return;
}

void clear_all_point_mp(int num, ...)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears all point_mp in an efficient way!!!             *
\***************************************************************/
{
  if (num > 0)
  {
    int i;
    va_list arg_addr;
 
    point_mp **next = (point_mp **)bmalloc(num * sizeof(point_mp *));
 
    // initialize arg_addr
    va_start(arg_addr, num);
 
    // setup the pointers
    for (i = 0; i < num; i++)
    {
      next[i] = (point_mp *) va_arg(arg_addr, void *);
 
      (*next[i])->size = 0; 
    }
 
    // clear arg_addr
    va_end(arg_addr);
 
    // clear the point
    for (i = 0; i < num; i++) 
      clear_point_mp(*next[i]);
 
    // free the pointers
    free(next);
  }
 
  return;
}

void move_to_patch_vec_d(point_d outPt, point_d inPt, vec_d patch)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined by patch               *
\***************************************************************/
{
  int i, size = inPt->size; 
  comp_d patchVal, recip;

  set_zero_d(patchVal);
  for (i = 0; i < size; i++)
  {
    sum_mul_d(patchVal, &inPt->coord[i], &patch->coord[i]);
  }
  
  // setup recip
  if (patchVal->r == 0 && patchVal->i == 0)
  { // generate a random perturbation so that we can divide
    get_comp_rand_d(recip);
    mul_rdouble_d(recip, recip, 1e-16);
    recip_d(recip, recip);
  }
  else
  { // reciprocate
    recip_d(recip, patchVal);
  }

  // setup outPt
  increase_size_point_d(outPt, size);
  outPt->size = size;
  for (i = 0; i < size; i++)
  {
    mul_d(&outPt->coord[i], &inPt->coord[i], recip);
  }

  return;
}

void move_to_patch_vec_mp(point_mp outPt, point_mp inPt, vec_mp patch)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined by patch               *
\***************************************************************/
{
  int i, size = inPt->size, num_digits = prec_to_digits(inPt->curr_prec); 
  mpf_t epsilon;
  comp_mp patchVal, recip;

  mpf_init(epsilon);
  init_mp(patchVal);
  init_mp(recip);

  set_zero_mp(patchVal);
  for (i = 0; i < size; i++)
  {
    sum_mul_mp(patchVal, &inPt->coord[i], &patch->coord[i]);
  }

  // setup recip
  if (mpfr_zero_p(patchVal->r) && mpfr_zero_p(patchVal->i))
  { // generate a random perturbation so that we can divide
    get_comp_rand_mp(recip);
    mpfr_ui_pow_ui(epsilon, 10, num_digits, __gmp_default_rounding_mode);
    mpf_ui_div(epsilon, 1, epsilon);
    mul_rmpf_mp(recip, recip, epsilon);
    recip_mp(recip, recip);
  }
  else
  { // reciprocate
    recip_mp(recip, patchVal);
  }

  // setup outPt
  increase_size_point_mp(outPt, size);
  outPt->size = size;
  for (i = 0; i < size; i++)
  {
    mul_mp(&outPt->coord[i], &inPt->coord[i], recip);
  }

  // clear memory
  mpf_clear(epsilon);
  clear_mp(patchVal);
  clear_mp(recip);

  return;
}

void move_to_patch_mat_d(point_d outPt, point_d inPt, mat_d patch_d, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined in patch_d & PPD       *
\***************************************************************/
{
  int i, j, var, count, num_var_gps = PPD->num_var_gp, num_gps = PPD->num_var_gp + PPD->num_hom_var_gp;

  if (num_gps == 0)
  { // there are no patches!
    point_cp_d(outPt, inPt);
  }
  else
  { // calculate the current value of the patch functions
    vec_d patchVals;
    init_vec_d(patchVals, 0);
    mul_mat_vec_d(patchVals, patch_d, inPt);

    count = 0;
    var = num_var_gps;

    // setup outPt
    increase_size_point_d(outPt, inPt->size);
    outPt->size = inPt->size;

    // adjust based on patchVals
    for (i = 0; i < num_gps; i++)
    { // see what type of group it is
      if (PPD->type[i] == 0)
      { // already homogenized - divide out by the value of the patch
        for (j = 0; j < PPD->size[i]; j++)
        {
          div_d(&outPt->coord[var], &inPt->coord[var], &patchVals->coord[i]);
          var++;
        }
      }
      else
      { // need to divide out by the value of the patch
        div_d(&outPt->coord[count], &inPt->coord[count], &patchVals->coord[i]);
        for (j = 0; j < PPD->size[i]; j++)
        {
          div_d(&outPt->coord[var], &inPt->coord[var], &patchVals->coord[i]);
          var++;
        }
        count++;
      }
    }
    // clear patchVals
    clear_vec_d(patchVals);
  }

  return;
}

void move_to_patch_d(point_d outPt, point_d inPt, basic_eval_data_d *BED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined in BED                 *
\***************************************************************/
{
  move_to_patch_mat_d(outPt, inPt, BED->patch.patchCoeff, &BED->preProcData);

  return;
}

void move_to_patch_mat_mp(point_mp outPt, point_mp inPt, mat_mp patch_mp, preproc_data *PPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined in patch_mp & PPD      *
\***************************************************************/
{
  int i, j, var, count, num_var_gps = PPD->num_var_gp, num_gps = PPD->num_var_gp + PPD->num_hom_var_gp;

  if (num_gps == 0)
  { // there are no patches!
    point_cp_mp(outPt, inPt);
  }
  else
  { // calculate the current value of the patch functions
    vec_mp patchVals;
    init_vec_mp(patchVals, 0);
    mul_mat_vec_mp(patchVals, patch_mp, inPt);

    count = 0;
    var = num_var_gps;

    // setup outPt
    increase_size_point_mp(outPt, inPt->size);
    outPt->size = inPt->size;

    // adjust based on patchVals
    for (i = 0; i < num_gps; i++)
    { // see what type of group it is
      if (PPD->type[i] == 0)
      { // already homogenized - divide out by the value of the patch
        for (j = 0; j < PPD->size[i]; j++)
        {
          div_mp(&outPt->coord[var], &inPt->coord[var], &patchVals->coord[i]);
          var++;
        }
      }
      else
      { // need to divide out by the value of the patch
        div_mp(&outPt->coord[count], &inPt->coord[count], &patchVals->coord[i]);
        for (j = 0; j < PPD->size[i]; j++)
        {
          div_mp(&outPt->coord[var], &inPt->coord[var], &patchVals->coord[i]);
          var++;
        }
        count++;
      }
    }
    // clear patchVals
    clear_vec_mp(patchVals);
  }

  return;
}

void move_to_patch_mp(point_mp outPt, point_mp inPt, basic_eval_data_mp *BED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: moves inPt to the patch defined in BED                 *
\***************************************************************/
{
  move_to_patch_mat_mp(outPt, inPt, BED->patch.patchCoeff, &BED->preProcData);

  return;
}

int checkWritePrivilege()
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - have write privilege, 1- no write privilege*
* NOTES: determines if this run of Bertini has write privileges *
\***************************************************************/
{
  int retVal = 0;
  FILE *TEMP = NULL;

  TEMP = fopen("tempFile", "w");

  if (TEMP == NULL)
  { // cannot write
    retVal = 1;
  }
  else
  { // can write
    retVal = 0;
    fclose(TEMP);
    remove("tempFile");
  }

  return retVal;
}

int prec_to_digits(int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: converts precision in bits to digits                   *
\***************************************************************/
{
  if (prec == 52)
    return 16;
  else
    return (int) floor(prec * log10(2.0) - 0.5);
}

int digits_to_prec(int digits)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: returns the smallest precision capable to handle digits*
\***************************************************************/
{
  if (digits <= 16)
    return 52;
  else
  {
    int prec_digits = (int) ceil((digits + 0.5) / (32 * log10(2.0)));
    if (prec_digits < 2)
    { // set to use 64-bit
      prec_digits = 64;
    }
    else
      prec_digits *= 32;

    return prec_digits;
  }
}

void extrinsicToIntrinsic_d(point_d out, point_d in, mat_d B_transpose, vec_d p)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the intrinic point associated to 'in'*
\***************************************************************/
{ // out = B^T*(in - p)

  vec_sub_d(out, in, p);
  mul_mat_vec_d(out, B_transpose, out);

  return;
}

void intrinsicToExtrinsic_d(point_d out, point_d in, mat_d B, vec_d p)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic point associated to 'in'*
\***************************************************************/
{ // out = p + B*in

  mul_mat_vec_d(out, B, in);
  vec_add_d(out, out, p);

  return;
}

void extrinsicToIntrinsic_mp(point_mp out, point_mp in, mat_mp B_transpose, vec_mp p)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the intrinic point associated to 'in'*
\***************************************************************/
{ // out = B^T*(in - p)

  vec_sub_mp(out, in, p);
  mul_mat_vec_mp(out, B_transpose, out);

  return;
}

void intrinsicToExtrinsic_mp(point_mp out, point_mp in, mat_mp B, vec_mp p)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic point associated to 'in'*
\***************************************************************/
{ // out = p + B*in

  mul_mat_vec_mp(out, B, in);
  vec_add_mp(out, out, p);

  return;
}

int isSamePoint(point_d endPt1_d, point_mp endPt1_mp, int prec1, point_d endPt2_d, point_mp endPt2_mp, int prec2, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - |endPt1 - endPt2| < tol, 0 - otherwise     *
* NOTES: determine if endPt1 and endPt2 are the same point (tol)*
\***************************************************************/
{
  int i, min_prec, size, size1, size2, retVal = 1;
  double norm, tol_sqr = tol * tol;

  // find the minimum precision
  min_prec = MIN(prec1, prec2);

  // find size
  if (prec1 < 64)
    size1 = endPt1_d->size;
  else
    size1 = endPt1_mp->size;

  if (prec2 < 64)
    size2 = endPt2_d->size;
  else
    size2 = endPt2_mp->size;

  // find the size of the point
  size = MIN(size1, size2); // should be the same

  if (prec1 < 64)
  { // first one is _d
    comp_d diff;

    if (prec2 < 64)
    { // second one is _d
      for (i = 0; i < size; i++)
      { // subtract ith coordinates
        sub_d(diff, &endPt1_d->coord[i], &endPt2_d->coord[i]);
        norm = norm_sqr_d(diff);

        // check against tol
        if (norm > tol_sqr)
        {
          retVal = 0;
          i = size;
        }
      }
    }
    else
    { // second one is _mp
      for (i = 0; i < size; i++)
      { // convert to _d
        mp_to_d(diff, &endPt2_mp->coord[i]);
        // subtract ith coordinates
        sub_d(diff, &endPt1_d->coord[i], diff);
        norm = norm_sqr_d(diff);

        // check against tol
        if (norm > tol_sqr)
        {
          retVal = 0;
          i = size;
        }
      }
    }
  }
  else
  { // first one is _mp
    if (prec2 < 64)
    { // second one is _d
      comp_d diff;

      for (i = 0; i < size; i++)
      { // convert to _d
        mp_to_d(diff, &endPt1_mp->coord[i]);
        // subtract ith coordinates
        sub_d(diff, diff, &endPt2_d->coord[i]);
        norm = norm_sqr_d(diff);

        // check against tol
        if (norm > tol_sqr)
        {
          retVal = 0; 
          i = size;
        }
      }
    }
    else
    { // both are in MP
      comp_mp t1, t2;

      init_mp2(t1, min_prec);
      init_mp2(t2, min_prec); 

      for (i = 0; i < size && retVal; i++)
      { // copy to t1 & t2
        set_mp(t1, &endPt1_mp->coord[i]);
        set_mp(t2, &endPt2_mp->coord[i]);

        // subtract
        sub_mp(t1, t1, t2);

        // find ||t1||^2
        mpf_mul(t1->r, t1->r, t1->r);
        mpf_mul(t1->i, t1->i, t1->i);
        mpf_add(t1->r, t1->r, t1->i);

        // check against tol
        if (mpf_get_d(t1->r) > tol_sqr)
          retVal = 0;
      }
      clear_mp(t1);
      clear_mp(t2);
    }
  }

  return retVal;
}

void findDiff_point(mpf_t norm_diff, point_d endPt1_d, point_mp endPt1_mp, int prec1, point_d endPt2_d, point_mp endPt2_mp, int prec2)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the norm of the difference of endPt1 and endPt2  *
\***************************************************************/
{
  int i, min_prec, size, size1, size2;

  // find the minimum precision
  min_prec = MIN(prec1, prec2);

  // find size
  if (prec1 < 64)
    size1 = endPt1_d->size;
  else
    size1 = endPt1_mp->size;

  if (prec2 < 64)
    size2 = endPt2_d->size;
  else
    size2 = endPt2_mp->size;

  size = MIN(size1, size2);

  if (prec1 < 64)
  { // first one is _d
    comp_d diff;
    double norm_d = 0;

    if (prec2 < 64)
    { // second one is _d
      for (i = 0; i < size; i++)
      { // subtract ith coordinates
        sub_d(diff, &endPt1_d->coord[i], &endPt2_d->coord[i]);
        diff->r = diff->r * diff->r + diff->i * diff->i;

        // update norm_d
        if (diff->r > norm_d)
          norm_d = diff->r;
      }
    }
    else
    { // second one is _mp
      for (i = 0; i < size; i++)
      { // convert to _d
        mp_to_d(diff, &endPt2_mp->coord[i]);
        // subtract ith coordinates
        sub_d(diff, &endPt1_d->coord[i], diff);
        diff->r = diff->r * diff->r + diff->i * diff->i;

        // update norm_d
        if (diff->r > norm_d)
          norm_d = diff->r;
      }
    }
    norm_d = sqrt(norm_d);

    // setup norm_diff 
    mpf_set_prec(norm_diff, 64);
    mpf_set_d(norm_diff, norm_d);
  }
  else
  { // first one is _mp
    if (prec2 < 64)
    { // second one is _d
      comp_d diff;
      double norm_d = 0;
      for (i = 0; i < size; i++)
      { // convert to _d
        mp_to_d(diff, &endPt1_mp->coord[i]);
        // subtract ith coordinates
        sub_d(diff, diff, &endPt2_d->coord[i]);
        diff->r = diff->r * diff->r + diff->i * diff->i;

        // update norm_d
        if (diff->r > norm_d)
          norm_d = diff->r;
      }
      norm_d = sqrt(norm_d);

      // setup norm_diff
      mpf_set_prec(norm_diff, 64);
      mpf_set_d(norm_diff, norm_d);
    }
    else
    { // both are in MP
      comp_mp t1, t2;
      init_mp2(t1, min_prec);
      init_mp2(t2, min_prec);

      // initialize norm_diff
      mpf_set_prec(norm_diff, min_prec);
      mpf_set_ui(norm_diff, 0);

      for (i = 0; i < size; i++)
      { // copy to t1 & t2
        set_mp(t1, &endPt1_mp->coord[i]);
        set_mp(t2, &endPt2_mp->coord[i]);

        // subtract
        sub_mp(t1, t1, t2);
  
        // find ||t1||^2
        mpf_mul(t1->r, t1->r, t1->r);
        mpf_mul(t1->i, t1->i, t1->i);
        mpf_add(t1->r, t1->r, t1->i);

        // update norm_diff
        if (mpf_cmp(t1->r, norm_diff) > 0)
          mpf_set(norm_diff, t1->r);
      }
      mpf_sqrt(norm_diff, norm_diff);

      // clear
      clear_mp(t1);
      clear_mp(t2);
    }
  }

  return;
}

void intrinsicToExtrinsicSlices_d(mat_d B_out, vec_d p_out, mat_d B_in, vec_d p_in)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in & p_in
// the conj-transpose of these ON vectors are B_out while p_out in the orthogonal projection of p_in onto B_in scaled to make p_out * p_in = 1 & p_out * B_in = 0
{
  int i, j, k, rows = B_in->rows, cols = B_in->cols;
  mat_d L, tempMat;
  comp_d sum, norm, tempComp;

  // error checking
  if (B_in->rows != p_in->size)
  {
    printf("ERROR: The sizes are not equal when converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (rows < cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  init_mat_d(L, rows, rows);
  init_mat_d(tempMat, rows, rows - cols - 1);

  // generate random vectors
  make_matrix_random_d(tempMat, rows, rows - cols - 1);

  // setup L = [B_in, p_in, tempMat]
  L->rows = L->cols = rows;
  for (j = 0; j < rows; j++)
    if (j < cols)
    { // copy B_in
      for (i = 0; i < rows; i++)
      {
        set_d(&L->entry[i][j], &B_in->entry[i][j]);
      }
    }
    else if (j == cols)
    { // copy p_in
      for (i = 0; i < rows; i++)
      {
        set_d(&L->entry[i][j], &p_in->coord[i]);
      }
    }
    else
    { // copy tempMat
      for (i = 0; i < rows; i++)
      {
        set_d(&L->entry[i][j], &tempMat->entry[i][j - cols - 1]);
      }
    }

  // orthonormalize L
  for (i = 0; i < rows; i++)
  { // orthogonalize
    for (j = 0; j < i; j++)
    { // update L(:,i) = L(:,i) - (L(:,j)'*L(:,i) / L(:,j)'*L(:,j)) * L(:,j)
      set_zero_d(sum);
      set_zero_d(norm);
      for (k = 0; k < rows; k++)
      {
        conjugate_d(tempComp, &L->entry[k][j]);
        sum_mul_d(sum, tempComp, &L->entry[k][i]);
        sum_mul_d(norm, tempComp, &L->entry[k][j]);
      }
      neg_d(sum, sum);
      sum->r /= norm->r;
      sum->i /= norm->r;
      for (k = 0; k < rows; k++)
      {
        sum_mul_d(&L->entry[k][i], sum, &L->entry[k][j]);
      }
    }
    // normalize
    set_zero_d(sum);
    for (k = 0; k < rows; k++)
    {
      conjugate_d(tempComp, &L->entry[k][i]);
      sum_mul_d(sum, tempComp, &L->entry[k][i]);
    }
    sum->r = 1 / sqrt(sum->r);
    sum->i = 0;
    for (k = 0; k < rows; k++)
    {    
      mul_d(&L->entry[k][i], sum, &L->entry[k][i]);
    }
  }

  // setup B_out^T
  change_size_mat_d(B_out, rows, rows - cols - 1);
  B_out->rows = rows;
  B_out->cols = rows - cols - 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < B_out->cols; j++)
    {
      set_d(&B_out->entry[i][j], &L->entry[i][j + cols + 1]);
    }
  // setup B_out
  transpose_d(B_out, B_out);

  // setup p_out
  change_size_vec_d(p_out, rows);
  p_out->size = rows;
  for (i = 0; i < rows; i++)
  {
    conjugate_d(&p_out->coord[i], &L->entry[i][cols]);
  }
  // correct p_out
  set_zero_d(sum);
  for (i = 0; i < rows; i++)
  { 
    sum_mul_d(sum, &p_out->coord[i], &p_in->coord[i]);
  }
  recip_d(sum, sum);
  for (i = 0; i < rows; i++)
  {
    mul_d(&p_out->coord[i], &p_out->coord[i], sum);
  }

  // clear
  clear_mat_d(L);
  clear_mat_d(tempMat);

  return;
}

void intrinsicToExtrinsicMat_d(mat_d B_out, mat_d B_in)
/***************************************************************\
* USAGE: used for changing orig = B_in * new to B_out * orig = 0*
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in 
// the conj-transpose of these ON vectors are B_out 
{
  int i, j, k, rows = B_in->rows, cols = B_in->cols;
  mat_d L, tempMat;
  comp_d sum, norm, tempComp;

  // error checking
  if (rows < cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  init_mat_d(L, rows, rows);
  init_mat_d(tempMat, rows, rows - cols);

  // generate random vectors
  make_matrix_random_d(tempMat, rows, rows - cols);

  // setup L = [B_in, tempMat]
  L->rows = L->cols = rows;
  for (j = 0; j < rows; j++)
    if (j < cols)
    { // copy B_in
      for (i = 0; i < rows; i++)
      {
        set_d(&L->entry[i][j], &B_in->entry[i][j]);
      }
    }
    else
    { // copy tempMat
      for (i = 0; i < rows; i++)
      {
        set_d(&L->entry[i][j], &tempMat->entry[i][j - cols]);
      }
    }

  // orthonormalize L
  for (i = 0; i < rows; i++)
  { // orthogonalize
    for (j = 0; j < i; j++)
    { // update L(:,i) = L(:,i) - (L(:,j)'*L(:,i) / L(:,j)'*L(:,j)) * L(:,j)
      set_zero_d(sum);
      set_zero_d(norm);
      for (k = 0; k < rows; k++)
      {
        conjugate_d(tempComp, &L->entry[k][j]);
        sum_mul_d(sum, tempComp, &L->entry[k][i]);
        sum_mul_d(norm, tempComp, &L->entry[k][j]);
      }
      neg_d(sum, sum);
      sum->r /= norm->r;
      sum->i /= norm->r;
      for (k = 0; k < rows; k++)
      {
        sum_mul_d(&L->entry[k][i], sum, &L->entry[k][j]);
      }
    }
    // normalize
    set_zero_d(sum);
    for (k = 0; k < rows; k++)
    {
      conjugate_d(tempComp, &L->entry[k][i]);
      sum_mul_d(sum, tempComp, &L->entry[k][i]);
    }
    sum->r = 1 / sqrt(sum->r);
    sum->i = 0;
    for (k = 0; k < rows; k++)
    {
      mul_d(&L->entry[k][i], sum, &L->entry[k][i]);
    }
  }

  // setup B_out^T
  change_size_mat_d(B_out, rows, rows - cols);
  B_out->rows = rows;
  B_out->cols = rows - cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < B_out->cols; j++)
    {
      set_d(&B_out->entry[i][j], &L->entry[i][j + cols]);
    }
  // setup B_out
  transpose_d(B_out, B_out);

  // clear
  clear_mat_d(L);
  clear_mat_d(tempMat);

  return;
}

void intrinsicToExtrinsicSlices_mp(mat_mp B_out, vec_mp p_out, mat_mp B_in, vec_mp p_in)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in & p_in
// the conj-transpose of these ON vectors are B_out while p_out in the orthogonal projection of p_in onto B_in scaled to make p_out * p_in = 1 & p_out * B_in = 0
{
  int i, j, k, rows = B_in->rows, cols = B_in->cols, curr_prec;
  mat_mp L, tempMat;
  comp_mp sum, norm, tempComp;

  // error checking
  if (B_in->rows != p_in->size)
  {
    printf("ERROR: The sizes are not equal when converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (rows < cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  init_mp(sum); init_mp(norm); init_mp(tempComp);
  init_mat_mp(L, rows, rows);
  init_mat_mp(tempMat, rows, rows - cols - 1);

  // generate random vectors
  curr_prec = mpf_get_prec(sum->r);
  make_matrix_random_mp(tempMat, rows, rows - cols - 1, curr_prec);

  // setup L = [B_in, p_in, tempMat]
  L->rows = L->cols = rows;
  for (j = 0; j < rows; j++)
    if (j < cols)
    { // copy B_in
      for (i = 0; i < rows; i++)
      {
        set_mp(&L->entry[i][j], &B_in->entry[i][j]);
      }
    }
    else if (j == cols)
    { // copy p_in
      for (i = 0; i < rows; i++)
      {
        set_mp(&L->entry[i][j], &p_in->coord[i]);
      }
    }
    else
    { // copy tempMat
      for (i = 0; i < rows; i++)
      {
        set_mp(&L->entry[i][j], &tempMat->entry[i][j - cols - 1]);
      }
    }

  // orthonormalize L
  for (i = 0; i < rows; i++)
  { // orthogonalize
    for (j = 0; j < i; j++)
    { // update L(:,i) = L(:,i) - (L(:,j)'*L(:,i) / L(:,j)'*L(:,j)) * L(:,j)
      set_zero_mp(sum);
      set_zero_mp(norm);
      for (k = 0; k < rows; k++)
      {
        conjugate_mp(tempComp, &L->entry[k][j]);
        sum_mul_mp(sum, tempComp, &L->entry[k][i]);
        sum_mul_mp(norm, tempComp, &L->entry[k][j]);
      }
      neg_mp(sum, sum);
      mpf_div(sum->r, sum->r, norm->r);
      mpf_div(sum->i, sum->i, norm->r);
      for (k = 0; k < rows; k++)
      {
        sum_mul_mp(&L->entry[k][i], sum, &L->entry[k][j]);
      }
    }
    // normalize
    set_zero_mp(sum);
    for (k = 0; k < rows; k++)
    {
      conjugate_mp(tempComp, &L->entry[k][i]);
      sum_mul_mp(sum, tempComp, &L->entry[k][i]);
    }
    mpf_ui_div(sum->r, 1, sum->r);
    mpf_set_ui(sum->i, 0);
    for (k = 0; k < rows; k++)
    {
      mul_mp(&L->entry[k][i], sum, &L->entry[k][i]);
    }
  }

  // setup B_out^T
  change_size_mat_mp(B_out, rows, rows - cols - 1);
  B_out->rows = rows;
  B_out->cols = rows - cols - 1;
  for (i = 0; i < rows; i++)
    for (j = 0; j < B_out->cols; j++)
    {
      set_mp(&B_out->entry[i][j], &L->entry[i][j + cols + 1]);
    }
  // setup B_out
  transpose_mp(B_out, B_out);

  // setup p_out
  change_size_vec_mp(p_out, rows);
  p_out->size = rows;
  for (i = 0; i < rows; i++)
  {
    conjugate_mp(&p_out->coord[i], &L->entry[i][cols]);
  }
  // correct p_out
  set_zero_mp(sum);
  for (i = 0; i < rows; i++)
  {
    sum_mul_mp(sum, &p_out->coord[i], &p_in->coord[i]);
  }
  recip_mp(sum, sum);
  for (i = 0; i < rows; i++)
  {
    mul_mp(&p_out->coord[i], &p_out->coord[i], sum);
  }

  // clear
  clear_mat_mp(L); clear_mat_mp(tempMat);
  clear_mp(sum); clear_mp(norm); clear_mp(tempComp);
 
  return;
}

void intrinsicToExtrinsicMat_mp(mat_mp B_out, mat_mp B_in)
/***************************************************************\
* USAGE: used for changing orig = B_in * new to B_out * orig = 0*
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in
// the conj-transpose of these ON vectors are B_out
{
  int i, j, k, rows = B_in->rows, cols = B_in->cols, curr_prec;
  mat_mp L, tempMat;
  comp_mp sum, norm, tempComp;

  // error checking
  if (rows < cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize
  init_mp(sum); init_mp(norm); init_mp(tempComp);
  init_mat_mp(L, rows, rows);
  init_mat_mp(tempMat, rows, rows - cols);

  // generate random vectors
  curr_prec = mpf_get_prec(sum->r);
  make_matrix_random_mp(tempMat, rows, rows - cols, curr_prec);

  // setup L = [B_in, tempMat]
  L->rows = L->cols = rows;
  for (j = 0; j < rows; j++)
    if (j < cols)
    { // copy B_in
      for (i = 0; i < rows; i++)
      {
        set_mp(&L->entry[i][j], &B_in->entry[i][j]);
      }
    }
    else
    { // copy tempMat
      for (i = 0; i < rows; i++)
      {
        set_mp(&L->entry[i][j], &tempMat->entry[i][j - cols]);
      }
    }

  // orthonormalize L
  for (i = 0; i < rows; i++)
  { // orthogonalize     
    for (j = 0; j < i; j++)
    { // update L(:,i) = L(:,i) - (L(:,j)'*L(:,i) / L(:,j)'*L(:,j)) * L(:,j)
      set_zero_mp(sum);
      set_zero_mp(norm);
      for (k = 0; k < rows; k++)
      {
        conjugate_mp(tempComp, &L->entry[k][j]);
        sum_mul_mp(sum, tempComp, &L->entry[k][i]);
        sum_mul_mp(norm, tempComp, &L->entry[k][j]);
      }
      neg_mp(sum, sum);
      mpf_div(sum->r, sum->r, norm->r);
      mpf_div(sum->i, sum->i, norm->r);
      for (k = 0; k < rows; k++)
      {
        sum_mul_mp(&L->entry[k][i], sum, &L->entry[k][j]);
      }
    }
    // normalize
    set_zero_mp(sum);
    for (k = 0; k < rows; k++)
    {
      conjugate_mp(tempComp, &L->entry[k][i]);
      sum_mul_mp(sum, tempComp, &L->entry[k][i]);
    }
    mpf_ui_div(sum->r, 1, sum->r);
    mpf_set_ui(sum->i, 0);
    for (k = 0; k < rows; k++)
    {
      mul_mp(&L->entry[k][i], sum, &L->entry[k][i]);
    }
  }

  // setup B_out^T
  change_size_mat_mp(B_out, rows, rows - cols);
  B_out->rows = rows;
  B_out->cols = rows - cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < B_out->cols; j++)
    {
      set_mp(&B_out->entry[i][j], &L->entry[i][j + cols]);
    }
  // setup B_out
  transpose_mp(B_out, B_out);

  // clear
  clear_mat_mp(L); clear_mat_mp(tempMat);
  clear_mp(sum); clear_mp(norm); clear_mp(tempComp);

  return;
}

void intrinsicToExtrinsicSlices_rat(mpq_t ****B_out, int *B_out_rows, int *B_out_cols, mpq_t ***p_out, int *p_out_size, mpq_t ***B_in, mpq_t **p_in, int B_in_rows, int B_in_cols, int p_in_size, int curr_prec, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in & p_in
// the conj-transpose of these ON vectors are B_out while p_out in the orthogonal projection of p_in onto B_in scaled to make p_out * p_in = 1 & p_out * B_in = 0
{
  int i, j;
  mat_mp B_in_mp, B_out_mp;
  vec_mp p_in_mp, p_out_mp;

  // error checking
  if (B_in_rows != p_in_size)
  {
    printf("ERROR: The sizes are not equal when converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }
  else if (B_in_rows < B_in_cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize globally to the maximum precision
  initMP(max_prec);

  init_mat_mp(B_in_mp, B_in_rows, B_in_cols); 
  init_vec_mp(p_in_mp, p_in_size);
  init_mat_mp(B_out_mp, 0, 0);
  init_vec_mp(p_out_mp, 0);

  // setup B_in_mp & p_in_mp
  B_in_mp->rows = p_in_mp->size = B_in_rows;
  B_in_mp->cols = B_in_cols;
  for (i = 0; i < B_in_rows; i++)
  {
    for (j = 0; j < B_in_cols; j++)
    {
      mpf_set_q(B_in_mp->entry[i][j].r, B_in[i][j][0]);
      mpf_set_q(B_in_mp->entry[i][j].i, B_in[i][j][1]);
    }
    mpf_set_q(p_in_mp->coord[i].r, p_in[i][0]);
    mpf_set_q(p_in_mp->coord[i].i, p_in[i][1]);
  }

  // setup B_out_mp & p_out_mp
  intrinsicToExtrinsicSlices_mp(B_out_mp, p_out_mp, B_in_mp, p_in_mp);

  // setup B_out
  *B_out_rows = B_out_mp->rows;
  *B_out_cols = B_out_mp->cols;
  init_mat_rat(*B_out, *B_out_rows, *B_out_cols);
  for (i = 0; i < *B_out_rows; i++)
    for (j = 0; j < *B_out_cols; j++)
    {
      mpf_t_to_rat((*B_out)[i][j][0], B_out_mp->entry[i][j].r);
      mpf_t_to_rat((*B_out)[i][j][1], B_out_mp->entry[i][j].i);
    }

  // setup p_out
  *p_out_size = p_out_mp->size; 
  init_vec_rat(*p_out, *p_out_size);
  for (i = 0; i < *p_out_size; i++)
  {
    mpf_t_to_rat((*p_out)[i][0], p_out_mp->coord[i].r);
    mpf_t_to_rat((*p_out)[i][1], p_out_mp->coord[i].i);
  }

  // clear
  clear_mat_mp(B_in_mp); clear_mat_mp(B_out_mp);
  clear_vec_mp(p_in_mp); clear_vec_mp(p_out_mp);

  // set globally back to the current precision
  initMP(curr_prec);

  return;
}

void intrinsicToExtrinsicMat_rat(mpq_t ****B_out, int *B_out_rows, int *B_out_cols, mpq_t ***B_in, int B_in_rows, int B_in_cols, int curr_prec, int max_prec)
/***************************************************************\
* USAGE: used for changing orig = B_in * new to B_out * orig = 0*
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup 'out' to be the extrinic slices associated to 'in*
\***************************************************************/
// use a modified Gram-Schmidt method to produce orthonormal vectors to the columns of B_in
// the conj-transpose of these ON vectors are B_out
{
  int i, j;
  mat_mp B_in_mp, B_out_mp;

  // error checking
  if (B_in_rows < B_in_cols)
  {
    printf("ERROR: There are not enough rows for converting intrinsic slices to extrinsic slices!\n");
    bexit(ERROR_INVALID_SIZE);
  }

  // initialize globally to the maximum precision
  initMP(max_prec);

  init_mat_mp(B_in_mp, B_in_rows, B_in_cols); 
  init_mat_mp(B_out_mp, 0, 0);

  // setup B_in_mp
  B_in_mp->rows = B_in_rows;
  B_in_mp->cols = B_in_cols;
  for (i = 0; i < B_in_rows; i++)
    for (j = 0; j < B_in_cols; j++)
    {
      mpf_set_q(B_in_mp->entry[i][j].r, B_in[i][j][0]);
      mpf_set_q(B_in_mp->entry[i][j].i, B_in[i][j][1]);
    }

  // setup B_out_mp
  intrinsicToExtrinsicMat_mp(B_out_mp, B_in_mp);

  // setup B_out
  *B_out_rows = B_out_mp->rows;
  *B_out_cols = B_out_mp->cols;
  init_mat_rat(*B_out, *B_out_rows, *B_out_cols);
  for (i = 0; i < *B_out_rows; i++)
    for (j = 0; j < *B_out_cols; j++)
    {
      mpf_t_to_rat((*B_out)[i][j][0], B_out_mp->entry[i][j].r);
      mpf_t_to_rat((*B_out)[i][j][1], B_out_mp->entry[i][j].i);
    }

  // clear
  clear_mat_mp(B_in_mp); clear_mat_mp(B_out_mp);

  // set globally back to the current precision
  initMP(curr_prec);

  return;
}

void printRawMat(FILE *OUT, mat_d M_d, mat_mp M_mp, mpq_t ***M_rat, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the matrix to OUT so that they could be reused  *
\***************************************************************/
{
  int i, j, rows, cols;

  if (MPType == 0 || MPType == 2)
  { // find sizes using M_d
    rows = M_d->rows;
    cols = M_d->cols;
  }
  else
  { // find sizes using M_mp
    rows = M_mp->rows; 
    cols = M_mp->cols; 
  }

  // print the number of rows & cols to OUT
  fprintf(OUT, "%d %d\n", rows, cols);

  // print M to OUT
  if (MPType == 0)
  { // print M_d
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        printRawComp(OUT, &M_d->entry[i][j], NULL, NULL, MPType);
      }
  }
  else if (MPType == 1)
  { // print M_mp
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        printRawComp(OUT, NULL, &M_mp->entry[i][j], NULL, MPType);
      }
  }
  else
  { // print M_rat
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        printRawComp(OUT, &M_d->entry[i][j], &M_mp->entry[i][j], M_rat[i][j], MPType);
      } 
  }

  fprintf(OUT, "\n");

  return;
}

void printRawComp(FILE *OUT, comp_d C_d, comp_mp C_mp, mpq_t *C_rat, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the number to OUT so that they could be reused  *
\***************************************************************/
{
  int digits = 0, base = 10;

  if (MPType == 0)
  { // print C_d
    fprintf(OUT, "%.15e %.15e\n", C_d->r, C_d->i);
  }
  else if (MPType == 1)
  { // print C_mp
    mpf_out_str(OUT, base, digits, C_mp->r);
    fprintf(OUT, " ");
    mpf_out_str(OUT, base, digits, C_mp->i);
    fprintf(OUT, "\n");
  }
  else
  { // print C_rat
    mpq_out_str(OUT, base, C_rat[0]);
    fprintf(OUT, " ");
    mpq_out_str(OUT, base, C_rat[1]);
    fprintf(OUT, "\n");
  }

  return;
}

void printRawVec(FILE *OUT, vec_d V_d, vec_mp V_mp, mpq_t **V_rat, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the vector to OUT so that they could be reused  *
\***************************************************************/
{
  int i, size;

  if (MPType == 0 || MPType == 2)
  { // find size using V_d
    size = V_d->size;
  }
  else
  { // find size using V_mp
    size = V_mp->size;
  }

  // print the size to OUT
  fprintf(OUT, "%d\n", size);

  // print V to OUT
  if (MPType == 0)
  { // print V_d
    for (i = 0; i < size; i++)
    {
      printRawComp(OUT, &V_d->coord[i], NULL, NULL, MPType);
    }
  }
  else if (MPType == 1)
  { // print V_mp
    for (i = 0; i < size; i++)
    {
      printRawComp(OUT, NULL, &V_mp->coord[i], NULL, MPType);
    }
  }
  else
  { // print V_rat
    for (i = 0; i < size; i++)
    {
      printRawComp(OUT, &V_d->coord[i], &V_mp->coord[i], V_rat[i], MPType);
    }
  }

  fprintf(OUT, "\n");

  return;
}

void setupRawComp(FILE *IN, comp_d C_d, comp_mp C_mp, mpq_t *C_rat, int C_Type, int inputType, int mp_need_init, int rat_need_init)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads in the number from IN and sets up C              *
\***************************************************************/
{
  int digits = 16, base = 10;

  if (inputType == 0)
  { // read in using double precision
    comp_d in_d;
    digits = 16;

    fscanf(IN, "%lf %lf\n", &in_d->r, &in_d->i);

    // setup C
    if (C_Type == 0)
    { // copy in_d to C_d
      set_d(C_d, in_d);
    }
    else if (C_Type == 1)
    { // copy in_d to C_mp
      if (mp_need_init)
      {
        init_mp(C_mp);
      }
      d_to_mp(C_mp, in_d);
    }
    else
    { // copy in_d to C_d, C_mp & C_rat
      if (mp_need_init)
      {
        init_mp(C_mp);
      }

      set_d(C_d, in_d);
      comp_d_to_mp_rat(C_mp, C_rat, C_d, digits, (int) mpf_get_prec(C_mp->r), 0, rat_need_init);
    }
  }
  else if (inputType == 1)
  { // read in using MP
    comp_mp in_mp;
    init_mp(in_mp);

    mpf_inp_str(in_mp->r, IN, base);
    mpf_inp_str(in_mp->i, IN, base);
    // scan rest of line
    fscanf(IN, "\n");

    // setup C
    if (C_Type == 0)
    { // copy in_mp to C_d
      mp_to_d(C_d, in_mp);
    }
    else if (C_Type == 1)
    { // copy in_mp to C_mp
      if (mp_need_init)
      {
        init_mp(C_mp);
      }
      set_mp(C_mp, in_mp);
    }
    else
    { // copy in_mp to C_d, C_mp & C_rat
      if (mp_need_init)
      {
        init_mp(C_mp);
      }
      if (rat_need_init)
      {
        mpq_init(C_rat[0]);
        mpq_init(C_rat[1]);
      }

      mpf_t_to_rat(C_rat[0], in_mp->r);
      mpf_t_to_rat(C_rat[1], in_mp->i);

      mpf_set_q(C_mp->r, C_rat[0]);
      mpf_set_q(C_mp->i, C_rat[1]);
      
      C_d->r = mpq_get_d(C_rat[0]);
      C_d->i = mpq_get_d(C_rat[1]);
    }

    clear_mp(in_mp);
  }
  else
  { // read in using rational 
    mpq_t in_rat[2];
    mpq_init(in_rat[0]); mpq_init(in_rat[1]);

    mpq_inp_str(in_rat[0], IN, base);
    mpq_inp_str(in_rat[1], IN, base);
    // scan rest of line
    fscanf(IN, "\n");

    // setup C
    if (C_Type == 0)
    { // copy in_rat to C_d
      C_d->r = mpq_get_d(in_rat[0]);
      C_d->i = mpq_get_d(in_rat[1]);
    }
    else if (C_Type == 1)
    { // copy in_rat to C_mp
      if (mp_need_init)
      {
        init_mp(C_mp);
      }

      mpf_set_q(C_mp->r, in_rat[0]);
      mpf_set_q(C_mp->i, in_rat[1]);
    }
    else
    { // copy in_rat to C_d, C_mp & C_rat
      if (mp_need_init)
      {
        init_mp(C_mp);
      }
      if (rat_need_init)
      {
        mpq_init(C_rat[0]);
        mpq_init(C_rat[1]);
      }

      mpq_set(C_rat[0], in_rat[0]);
      mpq_set(C_rat[1], in_rat[1]);

      mpf_set_q(C_mp->r, C_rat[0]);
      mpf_set_q(C_mp->i, C_rat[1]);

      C_d->r = mpq_get_d(C_rat[0]);
      C_d->i = mpq_get_d(C_rat[1]);
    }

    mpq_clear(in_rat[0]); mpq_clear(in_rat[1]);
  }

  return;
}

void setupRawVec(FILE *IN, vec_d V_d, vec_mp V_mp, mpq_t ***V_rat, int V_Type, int inputType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads in the vector from IN and sets up V              *
\***************************************************************/
{
  int i, size;

  // read in the size
  fscanf(IN, "%d\n", &size);

  // setup V
  if (V_Type == 0)
  { // setup V_d
    init_vec_d(V_d, size);
    V_d->size = size;

    for (i = 0; i < size; i++)
    {
      setupRawComp(IN, &V_d->coord[i], NULL, NULL, V_Type, inputType, 0, 0);
    }
  }
  else if (V_Type == 1)
  { // seutp V_mp
    init_vec_mp(V_mp, size);
    V_mp->size = size;

    for (i = 0; i < size; i++)
    {
      setupRawComp(IN, NULL, &V_mp->coord[i], NULL, V_Type, inputType, 0, 0);
    }
  }
  else
  { // setup V_d, V_mp & V_rat
    init_vec_d(V_d, size);
    init_vec_mp(V_mp, size);
    init_vec_rat(*V_rat, size);
    V_d->size = V_mp->size = size;

    for (i = 0; i < size; i++)
    {
      setupRawComp(IN, &V_d->coord[i], &V_mp->coord[i], (*V_rat)[i], V_Type, inputType, 0, 0);
    }
  }

  return;
}

void setupRawMat(FILE *IN, mat_d M_d, mat_mp M_mp, mpq_t ****M_rat, int M_Type, int inputType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: reads in the vector from IN and sets up M              *
\***************************************************************/
{
  int i, j, rows, cols;

  // read in the sizes
  fscanf(IN, "%d %d\n", &rows, &cols);

  // setup M
  if (M_Type == 0)
  { // setup M_d
    init_mat_d(M_d, rows, cols);
    M_d->rows = rows;
    M_d->cols = cols;

    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        setupRawComp(IN, &M_d->entry[i][j], NULL, NULL, M_Type, inputType, 0, 0);
      }
  }
  else if (M_Type == 1)
  { // seutp M_mp
    init_mat_mp(M_mp, rows, cols);
    M_mp->rows = rows;
    M_mp->cols = cols;

    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        setupRawComp(IN, NULL, &M_mp->entry[i][j], NULL, M_Type, inputType, 0, 0);
      }
  }
  else
  { // setup M_d, M_mp & M_rat
    init_mat_d(M_d, rows, cols);
    init_mat_mp(M_mp, rows, cols);
    init_mat_rat(*M_rat, rows, cols);
    M_d->rows = M_mp->rows = rows;
    M_d->cols = M_mp->cols = cols;

    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
      {
        setupRawComp(IN, &M_d->entry[i][j], &M_mp->entry[i][j], (*M_rat)[i][j], M_Type, inputType, 0, 0);
      }
  }

  return;
}

void setup_omp_file(FILE ***Fptr, FILE *Orig, char *name, int max)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the files needed for tracking with OpenMP        *
\***************************************************************/
{
  int i, size;
  char *str = NULL;

  if (max <= 0)
  {
    printf("\n\nERROR: The number of threads (%d) needs to be positive when setting up for tracking!\n", max);
    bexit(ERROR_CONFIGURATION);
  }
  else if (max == 1)
  { // allocate space to hold pointers to the files
    *Fptr = (FILE **)bmalloc(max * sizeof(FILE *));
    // simply point to Orig
    (*Fptr)[0] = Orig;
  }
  else
  { // allocate space to hold pointers to the files
    *Fptr = (FILE **)bmalloc(max * sizeof(FILE *));

    // setup each of the files
    for (i = 0; i < max; i++)
    {
      size = 1 + snprintf(NULL, 0, "%s_%d", name, i);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "%s_%d", name, i);
      (*Fptr)[i] = fopen(str, "w+");
    }
  }

  // free str
  free(str);

  return;
}

void combine_omp_file(FILE *F, FILE ***Fptr, int max)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: combines the information in Fptr into F then closes    *
* Fptr and releases memory                                      *
\***************************************************************/
{
  // if max == 1, Fptr[0] just points to F
  if (max == 1)
  { // simply set to NULL
    (*Fptr)[0] = NULL;
  }
  else
  { // copy data over to MIDOUT
    int i;
    char ch;
    for (i = 0; i < max; i++)
    { // rewind to beginning for reading
      rewind((*Fptr)[i]);
      // copy over to F
      ch = fgetc((*Fptr)[i]);
      while (ch != EOF)
      {
        fprintf(F, "%c", ch);
        ch = fgetc((*Fptr)[i]);
      }
      // close file
      fclose((*Fptr)[i]);
    }
  }

  // clear memory
  free(*Fptr);

  return;
}

void delete_omp_file(int max_threads, char *fName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: deletes the temporary files created for using OpenMP   *
\***************************************************************/
{
  int i;
  size_t size;
  char *str = NULL;
  FILE *TEST = NULL;

  for (i = 0; i < max_threads; i++)
  {
    size = 1 + snprintf(NULL, 0, "%s_%d", fName, i);
    str = (char *)brealloc(str, size * sizeof(char));
    sprintf(str, "%s_%d", fName, i);
    // clear data from file
    TEST = fopen(str, "w");
    fclose(TEST);
    remove(str);
  }

  free(str);

  return;
}

int sort_order_d(const void *vp, const void *vq)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: provides the order function for qsort for sortStruct_d *
\***************************************************************/
{
  sortStruct_d *a = (sortStruct_d *)vp;
  sortStruct_d *b = (sortStruct_d *)vq;
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

int sort_order_mp(const void *vp, const void *vq)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: provides the order function for qsort for sortStruct_mp*
\***************************************************************/
{
  int rV = 0;
  sortStruct_mp *a = (sortStruct_mp *)vp;
  sortStruct_mp *b = (sortStruct_mp *)vq;
  mpf_t x;

  mpf_init(x);

  mpf_sub(x, a->norm, b->norm);

  // return the sign of x
  rV = mpfr_sgn(x);

  mpf_clear(x);

  return rV;
}

double factorial(int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: n!                                             *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  double retVal = 1;

  if (n > 1)
  { // update retVal
    for (i = 2; i <= n; i++)
      retVal *= i;
  }

  return retVal;
}

double factorial_array(int *array, int length)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: array! = array[0]! * ... * array[length - 1]!  *
* NOTES:                                                        *
\***************************************************************/
{
  int i; double retVal = 1;

  for (i = 0; i < length; i++)
    if (array[i] > 1)
    {
      retVal *= factorial(array[i]);
    }

  return retVal;
}

void factorial2(mpf_t factorial, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: n!                                             *
* NOTES:                                                        *
\***************************************************************/
{
  int i;

  // initialize to 1
  mpf_set_ui(factorial, 1);

  if (n > 1)
  {
    for (i = 2; i <= n; i++)
      mpf_mul_ui(factorial, factorial, i);
  }

  return;
}

void factorial_array2(mpf_t rV, int *array, int length)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: array! = array[0]! * ... * array[length - 1]!  *
* NOTES:                                                        *
\***************************************************************/
{
  int i;
  mpf_t factorial;
  mpf_init(factorial);

  // intialize to 1
  mpf_set_ui(rV, 1);

  for (i = 0; i < length; i++)
    if (array[i] > 1)
    {
      factorial2(factorial, array[i]);
      mpf_mul(rV, rV, factorial);
    }

  mpf_clear(factorial);

  return;
}

double combination(int d, int n)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds C(d,n) - combination of d choose n               *
\***************************************************************/
{
  int i;
  double retVal;

  // error checking
  if (d >= n && n >= 0)
  { // find the smallest of n & d - n since C(d,n) = C(d,d-n)
    retVal = 1;
    n = MIN(n, d - n);

    for (i = 0; i < n; i++)
    { // update retVal
      retVal /= i + 1;
      retVal *= d - i;
    }
  }
  else
  { // no combinations
    retVal = 0;
  }

  return retVal;
}

int residual_to_digits_d(double residual)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines the number of digits described by residual  *
\***************************************************************/
{
  int digits;

  if (residual == 0)
  { // no error - all digits are correct
    digits = prec_to_digits(52);
  }
  else
  { // do the calculation
    digits = floor(-log10(residual));
  }

  if (digits < 1) // fail safe checking
    digits = 1;

  return digits;
}

int residual_to_digits_mp(mpf_t residual, int curr_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: determines the number of digits described by residual  *
\***************************************************************/
{
  int digits;

  if (mpfr_zero_p(residual))
  { // no error - all digits are correct
    digits = prec_to_digits(curr_prec);
  }
  else
  { // do the calculation
    mpf_t tempMPF;
    mpf_init2(tempMPF, curr_prec);

    mpfr_log10(tempMPF, residual, __gmp_default_rounding_mode);
    mpf_neg(tempMPF, tempMPF);
    mpfr_floor(tempMPF, tempMPF);
    digits = (int) mpf_get_si(tempMPF);

    mpf_clear(tempMPF);
  }

  if (digits < 1) // fail safe checking
    digits = 1;

  return digits;
}

int scanRestOfLine(FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether we are at the EOF or not               *
* NOTES: scan in the rest of the line                           *
\***************************************************************/
{
  char ch = fgetc(IN);

  while (ch != '\n' && ch != EOF)
    ch = fgetc(IN);

  return (ch == EOF);
}

int printRestOfLine(FILE *OUT, FILE *IN)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether we are at the EOF or not               *
* NOTES: print the rest of the line to OUT                      *
\***************************************************************/
{
  int endLoop = 0;
  char ch;

  do
  { // read in the next character
    ch = fgetc(IN);

    // determine what to do
    if (ch == EOF)
    { // end the loop
      endLoop = 1;
    }
    else if (ch == '\n')
    { // we are at the end of the line
      fprintf(OUT, "\n");
      endLoop = 1;
    }
    else if (ch == '\r')
    { // we have a carriage return - ignore it
      endLoop = 0;
    }
    else
    { // print the character
      fprintf(OUT, "%c", ch);
      endLoop = 0; 
    }
  } while (!endLoop);

  return (ch == EOF);
}

void setupInput(char *outputName, char *inputName)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: removes extra '\r' and puts in '\n' after ';'          *
\***************************************************************/
{
  int cont = 1;
  char ch;
  FILE *OUT = fopen("tempBertiniFile", "w"), *IN = fopen(inputName, "r");

  // verify files
  if (IN == NULL)
  {
    printf("\n\nERROR: '%s' does not exist!!!\n\n", inputName);
    remove("tempBertiniFile");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // put in end-of-line characters after each ';'
  do
  { // read in the next character
    ch = fgetc(IN);

    if (ch == EOF)
    { // end loop
      cont = 0;
    }
    else if (ch == '%') // commented out
    { // ignore rest of line
      scanRestOfLine(IN);
    }
    else if (ch == ';')
    { // put in a ';\n'
      fprintf(OUT, ";\n");
    }
    else if (ch != '\r') // ignore extra carriage returns
    { // print the character
      fprintf(OUT, "%c", ch);
    }

  } while (cont);

  // close IN & OUT
  fclose(IN);
  fclose(OUT);

  // setup outputName
  IN = fopen("tempBertiniFile", "r");
  OUT = fopen(outputName, "w");

  if (OUT == NULL)
  {
    printf("\n\nERROR: '%s' is an invalid name!!!\n\n", outputName);
    remove("tempBertiniFile");
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // remove extra spaces at beginning of the each line
  cont = 1;
  do
  { // read in first character of the next line
    ch = fgetc(IN);

    if (ch == EOF)
    { // end loop
      cont = 0;
    }
    else if (ch != '\n') // ignore blank lines
    { // remove extra spaces
      while (ch == ' ' && ch != EOF)
        ch = fgetc(IN);

      // see if at EOF
      if (ch == EOF)
      { // end loop
        cont = 0;
      }
      else if (ch != '\n') // ignore blank lines
      { // print the character
        fprintf(OUT, "%c", ch);

        // print the rest of the line
        cont = !printRestOfLine(OUT, IN);
      }
    }

  } while (cont);

  // close files
  fclose(OUT);
  fclose(IN);

  // remove temporary file
  remove("tempBertiniFile");

  return;
}

void printVersion(FILE *OUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints version of Bertini/GMP/MPFR & authors           *
\***************************************************************/
{
  fprintf(OUT, "*************** version information ***************\n");
  fprintf(OUT, "Bertini(TM) v%s\nGMP v%d.%d.%d, MPFR v%s\n", BERTINI_VERSION_STRING, __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
  fprintf(OUT, "\nAuthors:\nD.J. Bates, J.D. Hauenstein,\nA.J. Sommese, C.W. Wampler\n\n");

  return;
}

void normalize_vec_d(vec_d x1, vec_d x)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: normalize the vector x in the two norm sense           *
\***************************************************************/
{
  int i, size = x->size;
  double norm = 0;

  for (i = 0; i < size; i++)
    norm += norm_sqr_d(&x->coord[i]);
  norm = sqrt(norm);

  if (norm > 0)
  { // normalize x
    increase_size_vec_d(x1, x->size);
    x1->size = x->size;

    norm = 1 / norm;
    for (i = 0; i < size; i++)
      mul_rdouble_d(&x1->coord[i], &x->coord[i], norm);
  }
  else
  { // leave x alone
    vec_cp_d(x1, x);
  }

  return;
}

void normalize_vec_mp(vec_mp x1, vec_mp x)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: normalize the vector x in the two norm sense           *
\***************************************************************/
{
  int i, size = x->size;
  mpf_t norm, temp;

  mpf_init(norm); mpf_init(temp);
  mpf_set_ui(norm, 0);

  for (i = 0; i < size; i++)
  {
    norm_sqr_mp(temp, &x->coord[i]);
    mpf_add(norm, norm, temp);
  }
  mpf_sqrt(norm, norm);

  if (mpf_cmp_ui(norm, 0) > 0)
  { // normalize x
    increase_size_vec_mp(x1, x->size);
    x1->size = x->size;

    mpf_ui_div(norm, 1, norm);
    for (i = 0; i < size; i++)
      mul_rmpf_mp(&x1->coord[i], &x->coord[i], norm);
  }
  else
  { // leave x alone
    vec_cp_mp(x1, x);
  }

  mpf_clear(norm); mpf_clear(temp);
 
  return;
}

double point_cp_d_norm(point_d outPt, point_d inPt)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm of the point                              *
* NOTES: copy the point and find its norm                       *
\***************************************************************/
{
  int i, size = inPt->size;
  double test, norm = 0;

  change_size_point_d(outPt, size);
  outPt->size = size;

  for (i = 0; i < size; i++)
  {
    set_d(&outPt->coord[i], &inPt->coord[i]);
    test = norm_sqr_d(&outPt->coord[i]);
    if (norm < test)
      norm = test;
  }
  norm = sqrt(norm);

  return norm;
}

void point_cp_mp_norm(mpf_t norm, point_mp outPt, point_mp inPt, int currPrec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm of the point                              *
* NOTES: copy the point and find its norm                       *
\***************************************************************/ 
{
  int i, size = inPt->size;
  mpf_t temp;

  mpf_init2(temp, currPrec);

  mpf_set_prec(norm, currPrec);
  mpf_set_ui(norm, 0);

  setprec_point_mp(outPt, currPrec);
  change_size_point_mp(outPt, size);
  outPt->size = size;

  for (i = 0; i < size; i++)
  {
    set_mp(&outPt->coord[i], &inPt->coord[i]);
    norm_sqr_mp(temp, &outPt->coord[i]);
    if (mpf_cmp(norm, temp) < 0)
      mpf_set(norm, temp);
  }
  mpf_sqrt(norm, norm);

  mpf_clear(temp);

  return;
}

void point_cp_mp_norm2(mpf_t norm, point_mp inOutPt, int newPrec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm of the point                              *
* NOTES: increase the precision of the point and find its norm  *
\***************************************************************/
{
  int i, size = inOutPt->size;
  mpf_t temp;

  mpf_init2(temp, newPrec);

  mpf_set_prec(norm, newPrec);
  mpf_set_ui(norm, 0);

  for (i = 0; i < size; i++)
  {
    change_prec_mp(&inOutPt->coord[i], newPrec);
    norm_sqr_mp(temp, &inOutPt->coord[i]);
    if (mpf_cmp(norm, temp) < 0)
      mpf_set(norm, temp);
  }
  mpf_sqrt(norm, norm);

  mpf_clear(temp);

  return;
}

void point_d_to_mp_norm(mpf_t norm, point_mp outPt, point_d inPt, int currPrec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm of the point                              *
* NOTES: copy the point and find its norm                       *
\***************************************************************/
{
  int i, size = inPt->size;
  mpf_t temp;

  mpf_init2(temp, currPrec);

  mpf_set_prec(norm, currPrec);
  mpf_set_ui(norm, 0);

  setprec_point_mp(outPt, currPrec);
  change_size_point_mp(outPt, size);
  outPt->size = size;

  for (i = 0; i < size; i++)
  {
    d_to_mp(&outPt->coord[i], &inPt->coord[i]);
    norm_sqr_mp(temp, &outPt->coord[i]);
    if (mpf_cmp(norm, temp) < 0)
      mpf_set(norm, temp);
  }
  mpf_sqrt(norm, norm);

  mpf_clear(temp);

  return;
}

void mat_perp_d(mat_d Cp, mat_d C)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Compute Cp so that Cp*C = ID (C->rows >= C->cols)      *
\***************************************************************/
{ // ASSUME Cp != C
  int i, j, rows = C->rows, cols = C->cols;
  int its = 100; 
  double rank_tol = 1e-14, tol_conv_Jacobi = 1e-14, tol_QR = 1e-14, tol_sign = 1e-20, largeChange = 1e12;
  mat_d U, E, V, Einv;

  if (rows < cols)
  {
    printf("ERROR: The number of rows cannot be less than the number of cols!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (Cp == C)
  {
    printf("ERROR: The output matrix cannot be the same as the input matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // intialize
  init_mat_d(U, 0, 0);
  init_mat_d(E, 0, 0);
  init_mat_d(V, 0, 0);
  init_mat_d(Einv, cols, rows);
  Einv->rows = cols;
  Einv->cols = rows;

  // compute the SVD
  svd_jacobi_d(U, E, V, C, its, rank_tol, tol_conv_Jacobi, tol_QR, tol_sign, largeChange);

  // setup Einv
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
      if (i == j)
      { // recip
        recip_d(&Einv->entry[i][j], &E->entry[j][i]);
      }
      else
      { // set to 0
        set_zero_d(&Einv->entry[i][j]);
      }

  // set E = V * Einv
  change_size_mat_d(E, cols, rows);
  E->rows = cols; 
  E->cols = rows;
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
      if (j >= cols)
      { // set to 0
        set_zero_d(&E->entry[i][j]);
      }
      else
      { // multiply
        mul_d(&E->entry[i][j], &V->entry[i][j], &Einv->entry[j][j]);
      }

  // tranpose U
  transpose_d(Einv, U);

  // compute Cp = V * Einv * U'
  mat_mul_d(Cp, E, Einv);

  clear_mat_d(U);
  clear_mat_d(E);
  clear_mat_d(V);
  clear_mat_d(Einv);

  return;
}

void mat_perp_mp(mat_mp Cp, mat_mp C, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Compute Cp so that Cp*C = ID (C->rows >= C->cols)      *
\***************************************************************/
{ // ASSUME Cp != C
  int i, j, rows = C->rows, cols = C->cols;
  int its = 100, digits = prec_to_digits(prec) - 4;
  size_t size;
  char *str = NULL;
  mpf_t rank_tol, tol_conv_Jacobi, tol_QR, tol_sign, largeChange;
  mat_mp U, E, V, Einv;

  if (rows < cols)
  {
    printf("ERROR: The number of rows cannot be less than the number of cols!\n");
    bexit(ERROR_CONFIGURATION);
  }
  else if (Cp == C)
  {
    printf("ERROR: The output matrix cannot be the same as the input matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // intialize
  init_mat_mp2(U, 0, 0, prec);
  init_mat_mp2(E, 0, 0, prec);
  init_mat_mp2(V, 0, 0, prec);
  init_mat_mp2(Einv, cols, rows, prec);
  Einv->rows = cols;
  Einv->cols = rows;

  // setup tolerances
  mpf_init2(rank_tol, prec);
  mpf_init2(tol_conv_Jacobi, prec);
  mpf_init2(tol_QR, prec);
  mpf_init2(tol_sign, prec);
  mpf_init2(largeChange, prec);

  size = 1 + snprintf(NULL, 0, "1e-%d", digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", digits);
  mpf_set_str(rank_tol, str, 10);
  mpf_set_str(tol_conv_Jacobi, str, 10);
  mpf_set_str(tol_QR, str, 10);

  size = 1 + snprintf(NULL, 0, "1e-%d", 2 * (digits - 2));
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e-%d", 2 * (digits - 2));
  mpf_set_str(tol_sign, str, 10);

  size = 1 + snprintf(NULL, 0, "1e%d", digits);
  str = (char *)brealloc(str, size * sizeof(char));
  sprintf(str, "1e%d", digits);
  mpf_set_str(largeChange, str, 10);

  // compute the SVD
  svd_jacobi_mp(U, E, V, C, its, rank_tol, tol_conv_Jacobi, tol_QR, tol_sign, largeChange);

  // setup Einv
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
      if (i == j)
      { // recip
        recip_mp(&Einv->entry[i][j], &E->entry[j][i]);
      }
      else
      { // set to 0
        set_zero_mp(&Einv->entry[i][j]);
      }

  // set E = V * Einv
  change_size_mat_mp(E, cols, rows);
  E->rows = cols;
  E->cols = rows;
  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
      if (j >= cols)
      { // set to 0
        set_zero_mp(&E->entry[i][j]);
      }
      else
      { // multiply
        mul_mp(&E->entry[i][j], &V->entry[i][j], &Einv->entry[j][j]);
      }

  // tranpose U
  transpose_mp(Einv, U);

  // compute Cp = V * Einv * U'
  mat_mul_mp(Cp, E, Einv);

  clear_mat_mp(U);
  clear_mat_mp(E);
  clear_mat_mp(V);
  clear_mat_mp(Einv);
  mpf_clear(rank_tol);
  mpf_clear(tol_conv_Jacobi);
  mpf_clear(tol_QR);
  mpf_clear(tol_sign);
  mpf_clear(largeChange);

  return;
}

void mat_perp_rat(mat_d Cp_d, mat_mp Cp_mp, mpq_t ***Cp_rat, mpq_t ***C_rat, int rows, int cols, int curr_prec, int max_prec, int init_rat)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Compute Cp so that Cp*C = ID (C->rows >= C->cols)      *
\***************************************************************/
{
  int i, j;
  mat_mp tempCp, tempC;

  // set to the max precision
  initMP(max_prec);
  init_mat_mp(tempC, rows, cols);
  init_mat_mp(tempCp, cols, rows);
  tempC->rows = tempCp->cols = rows;
  tempC->cols = tempCp->rows = cols;

  // setup tempC
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      rat_to_mp(&tempC->entry[i][j], C_rat[i][j]);

  // setup tempCP
  mat_perp_mp(tempCp, tempC, max_prec);

  // set back to the current precision
  initMP(curr_prec);

  // setup Cp
  change_size_mat_d(Cp_d, cols, rows);
  setprec_mat_mp(Cp_mp, curr_prec);
  change_size_mat_mp(Cp_mp, cols, rows);
  if (init_rat)
    init_mat_rat(Cp_rat, cols, rows);
  Cp_d->rows = Cp_mp->rows = cols;
  Cp_d->cols = Cp_mp->cols = rows;

  for (i = 0; i < cols; i++)
    for (j = 0; j < rows; j++)
    {
      mp_to_rat(Cp_rat[i][j], &tempCp->entry[i][j]);
      rat_to_mp(&Cp_mp->entry[i][j], Cp_rat[i][j]);
      rat_to_d(&Cp_d->entry[i][j], Cp_rat[i][j]);
    }

  clear_mat_mp(tempC);
  clear_mat_mp(tempCp);

  return;
}

void print_comp_out_d(FILE *FP, comp_d A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints A to FP so that it can be read back in easily   *
\***************************************************************/
{
  print_d(FP, 0, A);
  fprintf(FP, "\n");

  return;
}

void print_vec_out_d(FILE *FP, vec_d A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints A to FP so that it can be read back in easily   *
\***************************************************************/
{
  int i, size = A->size;
  
  fprintf(FP, "%d\n", size);
  for (i = 0; i < size; i++)
  { 
    print_d(FP, 0, &A->coord[i]);
    fprintf(FP, "\n");
  }

  return;
}

void print_mat_out_d(FILE *FP, mat_d A)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints A to FP so that it can be read back in easily   *
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;
  
  fprintf(FP, "%d %d\n", rows, cols);
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      print_d(FP, 0, &A->entry[i][j]);
      fprintf(FP, "\n");
    }
    
  return;
}

void setup_comp_in_d(comp_d A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  fscanf(FP, "%lf%lf", &A->r, &A->i);
  scanRestOfLine(FP);

  return;
}

void setup_vec_in_d(vec_d A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i, size;

  fscanf(FP, "%d\n", &size);
  change_size_vec_d(A, size);
  A->size = size;
  for (i = 0; i < size; i++)
  {
    fscanf(FP, "%lf%lf", &A->coord[i].r, &A->coord[i].i);
    scanRestOfLine(FP);
  }

  return;
}

void setup_mat_in_d(mat_d A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i, j, rows, cols;

  fscanf(FP, "%d%d\n", &rows, &cols);
  change_size_mat_d(A, rows, cols);
  A->rows = rows;
  A->cols = cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      fscanf(FP, "%lf%lf", &A->entry[i][j].r, &A->entry[i][j].i);
      scanRestOfLine(FP);
    }

  return;
}

void print_comp_out_mp(FILE *FP, comp_mp A)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  print_mp(FP, 0, A);
  fprintf(FP, "\n");

  return;
}

void print_vec_out_mp(FILE *FP, vec_mp A)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  int i, size = A->size;

  fprintf(FP, "%d\n", size);
  for (i = 0; i < size; i++)
  {
    print_mp(FP, 0, &A->coord[i]);
    fprintf(FP, "\n");
  }

  return;
}

void print_mat_out_mp(FILE *FP, mat_mp A)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;

  fprintf(FP, "%d %d\n", rows, cols);
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      print_mp(FP, 0, &A->entry[i][j]);
      fprintf(FP, "\n");
    }

  return;
}

void setup_comp_in_mp(comp_mp A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  mpf_inp_str(A->r, FP, 10);
  mpf_inp_str(A->i, FP, 10);
  scanRestOfLine(FP);

  return;
}

void setup_vec_in_mp(vec_mp A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i, size;

  fscanf(FP, "%d\n", &size);
  change_size_vec_mp(A, size);
  A->size = size;
  for (i = 0; i < size; i++)
  {
    mpf_inp_str(A->coord[i].r, FP, 10);
    mpf_inp_str(A->coord[i].i, FP, 10);
    scanRestOfLine(FP);
  }

  return;
}

void setup_mat_in_mp(mat_mp A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i, j, rows, cols;

  fscanf(FP, "%d%d\n", &rows, &cols);
  change_size_mat_mp(A, rows, cols);
  A->rows = rows;
  A->cols = cols;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      mpf_inp_str(A->entry[i][j].r, FP, 10);
      mpf_inp_str(A->entry[i][j].i, FP, 10);
      scanRestOfLine(FP);
    }

  return;
}

void print_comp_out_rat(FILE *FP, mpq_t *A)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  print_rat(FP, A);
  fprintf(FP, "\n");

  return;
}

void print_vec_out_rat(FILE *FP, mpq_t **A, int size)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  int i;

  fprintf(FP, "%d\n", size);
  for (i = 0; i < size; i++)
  {
    print_rat(FP, A[i]);
    fprintf(FP, "\n");
  }

  return;
}

void print_mat_out_rat(FILE *FP, mpq_t ***A, int rows, int cols)
/***************************************************************\ 
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES:                                                * 
* NOTES: prints A to FP so that it can be read back in easily   * 
\***************************************************************/
{
  int i, j;

  fprintf(FP, "%d %d\n", rows, cols);
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    {
      print_rat(FP, A[i][j]);
      fprintf(FP, "\n");
    }

  return;
}

void setup_comp_in_rat(mpq_t *A, FILE *FP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  init_rat(A);
  mpq_inp_str(A[0], FP, 10);
  mpq_canonicalize(A[0]);
  mpq_inp_str(A[1], FP, 10);
  mpq_canonicalize(A[1]);
  scanRestOfLine(FP);

  return;
}

void setup_vec_in_rat(mpq_t ***A, FILE *FP, int *size)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i;

  fscanf(FP, "%d\n", size);
  init_vec_rat(*A, *size);
  for (i = 0; i < *size; i++)
  {
    mpq_inp_str((*A)[i][0], FP, 10);
    mpq_canonicalize((*A)[i][0]);
    mpq_inp_str((*A)[i][1], FP, 10);
    mpq_canonicalize((*A)[i][1]);
    scanRestOfLine(FP);
  }

  return;
}

void setup_mat_in_rat(mpq_t ****A, FILE *FP, int *rows, int *cols)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup A from FP                                        *
\***************************************************************/
{
  int i, j;

  fscanf(FP, "%d%d\n", rows, cols);
  init_mat_rat(*A, *rows, *cols);
  for (i = 0; i < *rows; i++)
    for (j = 0; j < *cols; j++)
    {
      mpq_inp_str((*A)[i][j][0], FP, 10);
      mpq_canonicalize((*A)[i][j][0]);
      mpq_inp_str((*A)[i][j][1], FP, 10);
      mpq_canonicalize((*A)[i][j][1]);
      scanRestOfLine(FP);
    }

  return;
}


