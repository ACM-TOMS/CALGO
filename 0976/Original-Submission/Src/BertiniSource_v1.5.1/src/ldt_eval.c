// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "localdim.h"

void setupSysStruct_from_deriv(systemStruct *sys, prog_deriv_t *SLP);
void setupDerivs_from_sys(funcStruct *derivs, funcStruct *subDerivs, int *totalOpCount, int **memLoc, int ***derivAddr, systemStruct *sys);
void setupProg_order_derivs(prog_deriv_t *deriv);

int evalDeriv_d(point_d funcVals, point_d derivVals, point_d linVals, point_d linDerivVals, point_d vars, prog_deriv_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate function and all partial derivatives          *
\***************************************************************/
{
  int i, j, begin, end, count, skipUpdate = 1, oid = thread_num();

  // check to see if the memory is the correct size
  if (Prog->memSize_d != Prog->memSizeNeeded)
  { // need to do update steps
    skipUpdate = 0;

    // reallocate memory
    Prog->mem_d = (_comp_d *)brealloc(Prog->mem_d, Prog->memSizeNeeded * sizeof(_comp_d));
    Prog->memSize_d = Prog->memSizeNeeded; // size of memory
  }

  // run the program
  if (skipUpdate)
  { // we are able to skip the update steps since they have been done before
    j = Prog->numInstAtEndUpdate;
  }
  else
  { // setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSizeNeeded; i++)
    {
      Prog->mem_d[i].i = 0;

      if (i >= begin && i < end)
        Prog->mem_d[i].r = mpq_get_d(Prog->nums[i - begin].rat);
      else
        Prog->mem_d[i].r = 0;
    }

    // setup I
    set_double_d(&Prog->mem_d[Prog->IAddr], 0, 1);
    // setup Pi
    set_double_d(&Prog->mem_d[Prog->IAddr + 1], M_PI, 0);

    // need to do the update steps
    j = 0;
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_d(&Prog->mem_d[i], &vars->coord[i - begin]);
  }

  // do the evaluations
  evalInsts_d(Prog->mem_d, Prog->prog, j, Prog->size, oid);

  // Gather the output for the evaluation
  change_size_vec_d(funcVals, Prog->numFuncs);
  funcVals->size = Prog->numFuncs;
  count = 0;
  for (i = 0; i < Prog->order; i++)
    count += Prog->numDerivs[i];
  change_size_vec_d(derivVals, count);
  derivVals->size = count;

  for (i = 0; i < funcVals->size; i++)
  {
    set_d(&funcVals->coord[i], &Prog->mem_d[Prog->evalFuncs + i]);
  }
  count = 0;
  for (i = 0; i < Prog->order; i++)
    for (j = 0; j < Prog->numDerivs[i]; j++)
    {
      set_d(&derivVals->coord[count], &Prog->mem_d[Prog->evalJFuncs[i] + j]);
      count++;
    }

  // evaluate the linear values = coeff * vars - const, if needed
  if (Prog->numLinears > 0)
  { // need to setup linVals & linDerivVals
    count = Prog->linearCoeff_d->rows;
    increase_size_vec_d(linVals, count);
    linVals->size = count;
    count *= Prog->linearCoeff_d->cols;
    increase_size_vec_d(linDerivVals, count);
    linDerivVals->size = count;

    count = 0;
    for (i = 0; i < Prog->linearCoeff_d->rows; i++)
    {
      neg_d(&linVals->coord[i], &Prog->linearConst_d->coord[i]);
      for (j = 0; j < Prog->linearCoeff_d->cols; j++)
      {
        sum_mul_d(&linVals->coord[i], &Prog->linearCoeff_d->entry[i][j], &vars->coord[j]);
        set_d(&linDerivVals->coord[count], &Prog->linearCoeff_d->entry[i][j]);
        count++;
      }
    }
  }
  else
  { // no linears
    linVals->size = linDerivVals->size = 0;
  }

  return 0;
}

int evalDeriv_mp(point_mp funcVals, point_mp derivVals, point_mp linVals, point_mp linDerivVals, point_mp vars, prog_deriv_t *Prog)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate function and all partial derivatives          *
\***************************************************************/
{
  int i, j, begin, end, count, skipUpdate = 1, oid = thread_num(), currPrec = Prog->precision;

  // check to see if the memory is the correct size
  if (Prog->memSize_mp != Prog->memSizeNeeded)
  { // need to do update steps
    skipUpdate = 0;

    // reallocate memory
    Prog->mem_mp = (_comp_mp *)brealloc(Prog->mem_mp, Prog->memSizeNeeded * sizeof(_comp_mp));
    Prog->memSize_mp = Prog->memSizeNeeded; // size of memory

    // initialize memory
    for (i = 0; i < Prog->memSize_mp; i++)
    {
      init_mp2(&Prog->mem_mp[i], currPrec);
    }
  }
  else if (mpfr_get_prec(Prog->mem_mp[0].r) != currPrec) // check to see if the precision is correct
  { // need to do update steps
    skipUpdate = 0;

    // set precision & set to zero
    for (i = 0; i < Prog->memSize_mp; i++)
    {
      setprec_mp(&Prog->mem_mp[i], currPrec);
    }

    if (Prog->linearType == 2)
    { // change precision on linear coeffs & consts
      change_prec_mat_mp_rat(Prog->linearCoeff_mp, currPrec, Prog->linearCoeff_rat);
      change_prec_vec_mp_rat(Prog->linearConst_mp, currPrec, Prog->linearConst_rat);
    }
  }

  // run the program
  if (skipUpdate)
  { // we are able to skip the update steps since they have been done before
    j = Prog->numInstAtEndUpdate;
  }
  else
  { // setup the numbers
    begin = Prog->numAddr;
    end = begin + Prog->numNums;
    for (i = 0; i < Prog->memSizeNeeded; i++)
    {
      mpf_set_ui(Prog->mem_mp[i].i, 0);

      if (i >= begin && i < end)
        mpf_set_q(Prog->mem_mp[i].r, Prog->nums[i - begin].rat);
      else
        mpf_set_ui(Prog->mem_mp[i].r, 0);
    }

    // setup I
    mpf_set_ui(Prog->mem_mp[Prog->IAddr].r, 0);
    mpf_set_ui(Prog->mem_mp[Prog->IAddr].i, 1);
    // setup Pi
    mpfr_const_pi(Prog->mem_mp[Prog->IAddr + 1].r, __gmp_default_rounding_mode);
    mpf_set_ui(Prog->mem_mp[Prog->IAddr + 1].i, 0);

    // need to do the update steps
    j = 0;
  }

  // copy in the values of the variables
  begin = Prog->inpVars;
  end = begin + Prog->numVars;
  for (i = begin; i < end; i++)
  {
    set_mp(&Prog->mem_mp[i], &vars->coord[i - begin]);
  }

  // do the evaluations
  evalInsts_mp(Prog->mem_mp, Prog->prog, j, Prog->size, oid);

  // Gather the output for the evaluation
  setprec_vec_mp(funcVals, currPrec);
  change_size_vec_mp(funcVals, Prog->numFuncs);
  funcVals->size = Prog->numFuncs;
  count = 0;
  for (i = 0; i < Prog->order; i++)
    count += Prog->numDerivs[i];
  setprec_vec_mp(derivVals, currPrec);
  change_size_vec_mp(derivVals, count);
  derivVals->size = count;

  for (i = 0; i < funcVals->size; i++)
  {
    set_mp(&funcVals->coord[i], &Prog->mem_mp[Prog->evalFuncs + i]);
  }
  count = 0;
  for (i = 0; i < Prog->order; i++)
    for (j = 0; j < Prog->numDerivs[i]; j++)
    {
      set_mp(&derivVals->coord[count], &Prog->mem_mp[Prog->evalJFuncs[i] + j]);
      count++;
    }

  // evaluate the linear values = coeff * vars - const, if needed
  if (Prog->numLinears > 0)
  { // need to setup linVals & linDerivVals
    count = Prog->linearCoeff_mp->rows;
    increase_size_vec_mp(linVals, count);
    linVals->size = count;
    count *= Prog->linearCoeff_mp->cols;
    increase_size_vec_mp(linDerivVals, count);
    linDerivVals->size = count;

    count = 0;
    for (i = 0; i < Prog->linearCoeff_mp->rows; i++)
    {
      neg_mp(&linVals->coord[i], &Prog->linearConst_mp->coord[i]);
      for (j = 0; j < Prog->linearCoeff_mp->cols; j++)
      {
        sum_mul_mp(&linVals->coord[i], &Prog->linearCoeff_mp->entry[i][j], &vars->coord[j]);
        set_mp(&linDerivVals->coord[count], &Prog->linearCoeff_mp->entry[i][j]);
        count++;
      }
    }
  }
  else
  { // no linears
    linVals->size = linDerivVals->size = 0;
  }

  return 0;
}

void add_patches_to_deriv(prog_deriv_t *deriv, mat_d patch_d, mat_mp patch_mp, mpq_t ***patch_rat, int patchMP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: assumes deriv is already setup - just add patch eqns   *
\***************************************************************/
{
  int i, count;

  if (patchMP == 0)
  { // setup const_d
    vec_d const_d;

    // find out how many
    count = patch_d->rows;

    // setup const_d
    init_vec_d(const_d, count);
    for (i = 0; i < count; i++)
      set_one_d(&const_d->coord[i]);

    // setup deriv
    add_linears_to_deriv(deriv, patch_d, const_d, NULL, NULL, NULL, NULL, patchMP);

    // clear const_d
    clear_vec_d(const_d);
  }
  else if (patchMP == 1)
  { // setup const_mp
    vec_mp const_mp;

    // find out how many
    count = patch_mp->rows;

    // setup const_mp
    init_vec_mp(const_mp, count);
    for (i = 0; i < count; i++)
      set_one_mp(&const_mp->coord[i]);

    // setup deriv
    add_linears_to_deriv(deriv, NULL, NULL, patch_mp, const_mp, NULL, NULL, patchMP);

    // clear const_mp
    clear_vec_mp(const_mp);
  }
  else
  { // setup const
    vec_d const_d;
    vec_mp const_mp;
    mpq_t **const_rat;

    // find out how many
    count = patch_d->rows;

    // setup const
    init_vec_d(const_d, count);
    init_vec_mp(const_mp, count);
    init_vec_rat(const_rat, count);
    for (i = 0; i < count; i++)
    {
      set_one_d(&const_d->coord[i]);
      set_one_mp(&const_mp->coord[i]);
      set_one_rat(const_rat[i]);
    }

    // setup deriv
    add_linears_to_deriv(deriv, patch_d, const_d, patch_mp, const_mp, patch_rat, const_rat, patchMP);

    // clear const
    clear_vec(const_d, const_mp, const_rat, 2);
  }

  return;
}

void add_slices_patch_to_deriv(prog_deriv_t *deriv, mat_d slices_d, mat_mp slices_mp, mpq_t ***slices_rat, vec_d patch_d, vec_mp patch_mp, mpq_t **patch_rat, int MP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: assumes deriv is already setup - just add slice/patch  *
\***************************************************************/
{
  int i, j, count, cols;

  if (MP == 0)
  { // setup coeff_d & const_d
    mat_d coeff_d;
    vec_d const_d;

    // find out how many
    count = slices_d->rows + 1;
    cols = slices_d->cols;

    // setup coeff_d & const_d
    init_mat_d(coeff_d, count, cols);
    init_vec_d(const_d, count);
    coeff_d->rows = const_d->size = count;
    coeff_d->cols = cols;
    for (i = 0; i < count; i++)
      if (i < count - 1)
      { // const == 0
        set_zero_d(&const_d->coord[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
          set_d(&coeff_d->entry[i][j], &slices_d->entry[i][j]);
      }
      else
      { // const == 1
        set_one_d(&const_d->coord[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
          set_d(&coeff_d->entry[i][j], &patch_d->coord[j]);
      }

    // setup deriv
    add_linears_to_deriv(deriv, coeff_d, const_d, NULL, NULL, NULL, NULL, MP);

    // clear coeff_d & const_d
    clear_mat_d(coeff_d);
    clear_vec_d(const_d);
  }
  else if (MP == 1)
  { // setup coeff_mp & const_mp
    mat_mp coeff_mp;
    vec_mp const_mp;

    // find out how many
    count = slices_mp->rows + 1;
    cols = slices_mp->cols;

    // setup coeff_mp & const_mp
    init_mat_mp(coeff_mp, count, cols);
    init_vec_mp(const_mp, count);
    coeff_mp->rows = const_mp->size = count;
    coeff_mp->cols = cols;
    for (i = 0; i < count; i++)
      if (i < count - 1)
      { // const == 0
        set_zero_mp(&const_mp->coord[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
          set_mp(&coeff_mp->entry[i][j], &slices_mp->entry[i][j]);
      }
      else
      { // const == 1
        set_one_mp(&const_mp->coord[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
          set_mp(&coeff_mp->entry[i][j], &patch_mp->coord[j]);
      }

    // setup deriv
    add_linears_to_deriv(deriv, NULL, NULL, coeff_mp, const_mp, NULL, NULL, MP);

    // clear coeff_mp & const_mp
    clear_mat_mp(coeff_mp);
    clear_vec_mp(const_mp);
  }
  else
  { // setup const
    mat_d coeff_d;
    mat_mp coeff_mp;
    mpq_t ***coeff_rat;

    vec_d const_d;
    vec_mp const_mp;
    mpq_t **const_rat;

    // find out how many
    count = slices_d->rows + 1;
    cols = slices_d->cols;

    // setup coeff & const
    init_mat_d(coeff_d, count, cols); init_mat_mp(coeff_mp, count, cols); init_mat_rat(coeff_rat, count, cols);
    init_vec_d(const_d, count);       init_vec_mp(const_mp, count);       init_vec_rat(const_rat, count);
    coeff_d->rows = coeff_mp->rows = const_d->size = const_mp->size = count;
    coeff_d->cols = coeff_mp->cols = cols;
    for (i = 0; i < count; i++)
      if (i < count - 1)
      { // const == 0
        set_zero_d(&const_d->coord[i]);
        set_zero_mp(&const_mp->coord[i]);
        set_zero_rat(const_rat[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
        {
          set_d(&coeff_d->entry[i][j], &slices_d->entry[i][j]);
          set_mp(&coeff_mp->entry[i][j], &slices_mp->entry[i][j]);
          set_rat(coeff_rat[i][j], slices_rat[i][j]);
        }
      }
      else
      { // const == 1
        set_one_d(&const_d->coord[i]);
        set_one_mp(&const_mp->coord[i]);
        set_one_rat(const_rat[i]);

        // setup coeff
        for (j = 0; j < cols; j++)
        {
          set_d(&coeff_d->entry[i][j], &patch_d->coord[j]);
          set_mp(&coeff_mp->entry[i][j], &patch_mp->coord[j]);
          set_rat(coeff_rat[i][j], patch_rat[j]);
        }
      }

    // setup deriv
    add_linears_to_deriv(deriv, coeff_d, const_d, coeff_mp, const_mp, coeff_rat, const_rat, MP);

    // clear const
    clear_mat(coeff_d, coeff_mp, coeff_rat, 2);
    clear_vec(const_d, const_mp, const_rat, 2);
  }

  return;
}

void add_linears_to_deriv(prog_deriv_t *deriv, mat_d linears_d, vec_d const_d, mat_mp linears_mp, vec_mp const_mp, mpq_t ***linears_rat, mpq_t **const_rat, int linearsMP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: assumes deriv is already setup - just add linear eqns  *
\***************************************************************/
{
  int i, j, rows, cols;

  if (linearsMP == 0)
  { // setup _d
    deriv->numLinears = rows = linears_d->rows;
    cols = linears_d->cols;

    // setup linearType
    deriv->linearType = 0;
  
    // initialize to correct size
    init_mat_d(deriv->linearCoeff_d, rows, cols);
    deriv->linearCoeff_d->rows = rows;
    deriv->linearCoeff_d->cols = cols;
    init_vec_d(deriv->linearConst_d, rows);
    deriv->linearConst_d->size = rows;

    // setup coeff & const
    for (i = 0; i < rows; i++)
    { // setup const
      set_d(&deriv->linearConst_d->coord[i], &const_d->coord[i]);
      for (j = 0; j < cols; j++)
      { // setup coeff
        set_d(&deriv->linearCoeff_d->entry[i][j], &linears_d->entry[i][j]);
      }
    }
  }
  else if (linearsMP == 1)
  { // setup _mp
    deriv->numLinears = rows = linears_mp->rows;
    cols = linears_mp->cols;

    // setup linearType
    deriv->linearType = 1;

    // initialize to correct size
    init_mat_mp(deriv->linearCoeff_mp, rows, cols);
    deriv->linearCoeff_mp->rows = rows;
    deriv->linearCoeff_mp->cols = cols;
    init_vec_mp(deriv->linearConst_mp, rows);
    deriv->linearConst_mp->size = rows;

    // setup coeff & const
    for (i = 0; i < rows; i++)
    { // setup const
      set_mp(&deriv->linearConst_mp->coord[i], &const_mp->coord[i]);
      for (j = 0; j < cols; j++)
      { // setup coeff
        set_mp(&deriv->linearCoeff_mp->entry[i][j], &linears_mp->entry[i][j]);
      }
    }
  }
  else
  { // setup for AMP
    deriv->numLinears = rows = linears_d->rows;
    cols = linears_d->cols;

    // setup linearType
    deriv->linearType = 2;

    // initialize to correct size
    init_mat_d(deriv->linearCoeff_d, rows, cols);
    init_mat_mp2(deriv->linearCoeff_mp, rows, cols, deriv->precision);
    init_mat_rat(deriv->linearCoeff_rat, rows, cols);
    deriv->linearCoeff_d->rows = deriv->linearCoeff_mp->rows = rows;
    deriv->linearCoeff_d->cols = deriv->linearCoeff_mp->cols = cols;
    init_vec_d(deriv->linearConst_d, rows);
    init_vec_mp2(deriv->linearConst_mp, rows, deriv->precision);
    init_vec_rat(deriv->linearConst_rat, rows);
    deriv->linearConst_d->size = deriv->linearConst_mp->size = rows;

    // setup coeff & const
    for (i = 0; i < rows; i++)
    { // setup const
      set_rat(deriv->linearConst_rat[i], const_rat[i]);
      rat_to_d(&deriv->linearConst_d->coord[i], deriv->linearConst_rat[i]);
      rat_to_mp(&deriv->linearConst_mp->coord[i], deriv->linearConst_rat[i]);
      for (j = 0; j < cols; j++)
      { // setup coeff
        set_rat(deriv->linearCoeff_rat[i][j], linears_rat[i][j]);
        rat_to_d(&deriv->linearCoeff_d->entry[i][j], deriv->linearCoeff_rat[i][j]);
        rat_to_mp(&deriv->linearCoeff_mp->entry[i][j], deriv->linearCoeff_rat[i][j]);
      }
    }
  }

  return;
}

void setup_deriv_from_SLP(prog_deriv_t *deriv, prog_t *SLP)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup deriv from SLP                                   *
\***************************************************************/
{
  int i, currPrec = SLP->precision;

  // setup prog
  deriv->size = SLP->size;
  deriv->prog = (int *)bmalloc(deriv->size * sizeof(int));
  for (i = 0; i < deriv->size; i++)
    deriv->prog[i] = SLP->prog[i];
 
  // setup mem
  deriv->memSizeNeeded = SLP->memSize;
  deriv->memSize_d = deriv->memSize_mp = 0;
  deriv->mem_d = NULL;
  deriv->mem_mp = NULL;

  // seutp nums
  deriv->numNums = SLP->numNums;
  deriv->precision = SLP->precision;
  deriv->nums = (num_t *)bmalloc(deriv->numNums * sizeof(num_t));
  for (i = 0; i < deriv->numNums; i++)
  { // initialize & copy
    mpq_init(deriv->nums[i].rat);
    mpq_set(deriv->nums[i].rat, SLP->nums[i].rat);

    mpf_init2(deriv->nums[i].real, currPrec);
    mpf_set(deriv->nums[i].real, SLP->nums[i].real);
  }

  // assume no linears
  deriv->numLinears = 0;
  deriv->linearType = -1;

  // setup order
  deriv->order = 1; // the original SLP contains all partial derivatives

  // setup instruction numbers
  deriv->numInstAtEndUpdate = SLP->numInstAtEndUpdate;

  // setup other values
  deriv->numVars = SLP->numVars;
  deriv->numNums = SLP->numNums;
  deriv->numConsts = SLP->numConsts;
  deriv->numFuncs = SLP->numFuncs;
  deriv->numSubfuncs = SLP->numSubfuncs;
  deriv->numDerivs = (int *)bmalloc(deriv->order * sizeof(int));
  deriv->numDerivs[0] = SLP->numFuncs * SLP->numVars; 
  deriv->numSubDerivs = (int *)bmalloc(deriv->order * sizeof(int));
  deriv->numSubDerivs[0] = SLP->numSubfuncs * SLP->numVars;

  deriv->inpVars = SLP->inpVars;
  deriv->IAddr = SLP->IAddr;
  deriv->numAddr = SLP->numAddr;
  deriv->constAddr = SLP->constAddr;

  deriv->evalFuncs = SLP->evalFuncs;
  deriv->evalSubs = SLP->evalSubs;
  deriv->evalJFuncs = (int *)bmalloc(deriv->order * sizeof(int));
  deriv->evalJFuncs[0] = SLP->evalJVars;
  deriv->evalJSubs = (int *)bmalloc(deriv->order * sizeof(int));
  deriv->evalJSubs[0] = SLP->evalJSubsV;

  // setup sys from deriv
  setupSysStruct_from_deriv(&deriv->sys, deriv);

  // setup the derivative structures
  deriv->subFuncDerivs = (funcStruct **)bmalloc(deriv->order * sizeof(funcStruct *));
  deriv->subFuncDerivs[0] = (funcStruct *)bmalloc(deriv->numSubDerivs[0] * sizeof(funcStruct));
  deriv->funcDerivs = (funcStruct **)bmalloc(deriv->order * sizeof(funcStruct *));
  deriv->funcDerivs[0] = (funcStruct *)bmalloc(deriv->numDerivs[0] * sizeof(funcStruct));

  // finish setting everything else up
  setupDerivs_from_sys(deriv->funcDerivs[0], deriv->subFuncDerivs[0], &deriv->totalOpCount, &deriv->memLoc, &deriv->derivAddr, &deriv->sys);
  setupProg_order_derivs(deriv);

  return;
}

void clear_deriv(prog_deriv_t *deriv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear deriv                                            *
\***************************************************************/
{ 
  int i;

  // clear nums
  clearNums(&deriv->nums, deriv->numNums);

  // clear mem_d
  if (deriv->memSize_d > 0)
  { // release mem_d
    free(deriv->mem_d);
  }
  
  // clear mem_mp
  if (deriv->memSize_mp > 0)
  { // clear mem_mp
    for (i = deriv->memSize_mp - 1; i >= 0; i--)
    {
      clear_mp(&deriv->mem_mp[i]);
    }
    free(deriv->mem_mp);
  }

  // clear linears
  clear_mat(deriv->linearCoeff_d, deriv->linearCoeff_mp, deriv->linearCoeff_rat, deriv->linearType);
  clear_vec(deriv->linearConst_d, deriv->linearConst_mp, deriv->linearConst_rat, deriv->linearType); 

  // clear memLoc & derivAddr
  free(deriv->memLoc);
  for (i = 0; i < deriv->totalOpCount; i++)
    free(deriv->derivAddr[i]);
  free(deriv->derivAddr); 

  // clear funcDerivs & subFuncDerivs
  for (i = 0; i < deriv->order; i++)
  {
    clearFuncStruct(deriv->subFuncDerivs[i], deriv->numSubDerivs[i]);
    clearFuncStruct(deriv->funcDerivs[i], deriv->numDerivs[i]);
  }
  free(deriv->subFuncDerivs);
  free(deriv->funcDerivs);
  free(deriv->numSubDerivs);
  free(deriv->numDerivs);

  // clear sys
  clearSystemStructure(&deriv->sys);

  // clear other structures
  free(deriv->prog);
  free(deriv->evalJFuncs);
  free(deriv->evalJSubs);

  return;
}

void setupSysStruct_from_deriv(systemStruct *sys, prog_deriv_t *deriv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup sys from deriv                                   *
\***************************************************************/
{
  int i, j, inst, count, boolCont, isSubFunc, isFunc, subFuncNum, funcNum, total = deriv->numSubfuncs + deriv->numFuncs;

  // setup variables
  sys->numVars = deriv->numVars;
  sys->varsAddr = deriv->inpVars;

  // setup path variables - none
  sys->numPathVars = sys->pathVarsAddr = 0;

  // setup parameters - none
  sys->numParams = sys->paramAddr = 0;
  sys->params = NULL;

  // setup functions
  sys->numFuncs = deriv->numFuncs;
  sys->funcAddr = deriv->evalFuncs;
  sys->funcs = (funcStruct *)bmalloc(sys->numFuncs * sizeof(funcStruct));
  sys->funcOrder = (int *)bmalloc(sys->numFuncs * sizeof(int));
  for (i = 0; i < sys->numFuncs; i++)
  {
    sys->funcs[i].num_ops = 1;
    sys->funcs[i].ops = (func_ops *)bmalloc(sys->funcs[i].num_ops * sizeof(func_ops));

    sys->funcOrder[i] = i;
  }

  // setup subfunctions
  sys->numSubfuncs = deriv->numSubfuncs;
  sys->subFuncAddr = deriv->evalSubs;
  sys->subFuncs = (funcStruct *)bmalloc(sys->numSubfuncs * sizeof(funcStruct));
  sys->subFuncOrder = (int *)bmalloc(sys->numSubfuncs * sizeof(int));
  for (i = 0; i < sys->numSubfuncs; i++)
  {
    sys->subFuncs[i].num_ops = 1;
    sys->subFuncs[i].ops = (func_ops *)bmalloc(sys->subFuncs[i].num_ops * sizeof(func_ops));

    sys->subFuncOrder[i] = i;
  }

  // setup other information
  sys->firstFreeMemLoc = deriv->memSizeNeeded;
  sys->numConstants = deriv->numConsts;
  sys->constAddr = deriv->constAddr; 
  sys->numNumbers = deriv->numNums;
  sys->numAddr = deriv->numAddr;

  // setup update instructions, if needed
  if (deriv->numInstAtEndUpdate > 0)
  { // setup the update instructions
    sys->numUpdate = 1;
    sys->updateOps = (func_ops *)bmalloc(sys->numUpdate * sizeof(func_ops));

    inst = count = 0;
    do
    { // setup the update operation
      sys->updateOps[count].op = deriv->prog[inst];
      sys->updateOps[count].memLoc = deriv->prog[inst+1];

      if (isUnary(deriv->prog[inst]))
      { // unary operation
        sys->updateOps[count].in[0] = deriv->prog[inst+2];
        sys->updateOps[count].in[1] = -1;

        // increment inst
        inst += 3;
      }
      else
      { // binary operation
        sys->updateOps[count].in[0] = deriv->prog[inst+2];
        sys->updateOps[count].in[1] = deriv->prog[inst+3];

        // increment inst
        inst += 4;
      }

      // increment count
      count++;

      // see if there are more update operations
      if (inst < deriv->numInstAtEndUpdate)
      { // see if we need more room
        if (count >= sys->numUpdate)
        { // allocate more memory
          sys->numUpdate *=  2;
          sys->updateOps = (func_ops *)brealloc(sys->updateOps, sys->numUpdate * sizeof(func_ops));
        }
      }
    } while (inst < deriv->numInstAtEndUpdate);

    // set everything exactly
    sys->numUpdate = count;
    sys->updateOps = (func_ops *)brealloc(sys->updateOps, sys->numUpdate * sizeof(func_ops));
  }
  else
  { // no update instructions
    sys->numUpdate = 0;
    sys->updateOps = NULL;
  }

  // setup the functions & subfunctions
  inst = deriv->numInstAtEndUpdate;
  for (i = 0; i < total; i++)
  { // find the end of the next function/subfunction
    j = count = isSubFunc = isFunc = 0;
    boolCont = 1;
    subFuncNum = funcNum = -1;

    do
    { // determine if the operation calculated the end of a subfunction
      if (sys->subFuncAddr <= deriv->prog[inst+1+j] && deriv->prog[inst+1+j] < sys->subFuncAddr + sys->numSubfuncs)
      { // this is the last operation of a subfunction
        boolCont = 0;
        isSubFunc = 1;
        subFuncNum = deriv->prog[inst+1+j] - sys->subFuncAddr;
      }

      // determine if the operation calculated the end of a function
      if (boolCont && sys->funcAddr <= deriv->prog[inst+1+j] && deriv->prog[inst+1+j] < sys->funcAddr + sys->numFuncs)
      { // this is the last operation of a function
        boolCont = 0;
        isFunc = 1;
        funcNum = deriv->prog[inst+1+j] - sys->funcAddr;
      }

      // increment count - number of operations needed
      count++;

      // increment j - move past this operation
      if (isUnary(deriv->prog[inst+j]))
        j += 3; // unary
      else
        j += 4; // binary

    } while (boolCont); 

    // store the subfunction/function
    if (isSubFunc)
    { // store the size and allocte the ops for this subfunction
      sys->subFuncOrder[subFuncNum] = i;
      sys->subFuncs[subFuncNum].num_ops = count;
      sys->subFuncs[subFuncNum].ops = (func_ops *)bmalloc(count * sizeof(func_ops));

      // store the ops
      for (count = 0; count < sys->subFuncs[subFuncNum].num_ops; count++)
      {
        sys->subFuncs[subFuncNum].ops[count].op = deriv->prog[inst];
        sys->subFuncs[subFuncNum].ops[count].memLoc = deriv->prog[inst+1];
        sys->subFuncs[subFuncNum].ops[count].in[0] = deriv->prog[inst+2];

        if (isUnary(deriv->prog[inst]))
        {
          sys->subFuncs[subFuncNum].ops[count].in[1] = -1;
          inst += 3;
        }
        else
        {
          sys->subFuncs[subFuncNum].ops[count].in[1] = deriv->prog[inst+3];
          inst += 4;
        }

        // initialize other data
        sys->subFuncs[subFuncNum].ops[count].lastUsed = -1;
        sys->subFuncs[subFuncNum].ops[count].numDiffInst = 0;
      }
    }
    else if (isFunc)
    { // store the size and allocte the ops for this function
      sys->funcOrder[funcNum] = i;
      sys->funcs[funcNum].num_ops = count;
      sys->funcs[funcNum].ops = (func_ops *)bmalloc(count * sizeof(func_ops));

      // store the ops
      for (count = 0; count < sys->funcs[funcNum].num_ops; count++)
      {
        sys->funcs[funcNum].ops[count].op = deriv->prog[inst];
        sys->funcs[funcNum].ops[count].memLoc = deriv->prog[inst+1];
        sys->funcs[funcNum].ops[count].in[0] = deriv->prog[inst+2];

        if (isUnary(deriv->prog[inst]))
        {
          sys->funcs[funcNum].ops[count].in[1] = -1;
          inst += 3;
        }
        else
        {
          sys->funcs[funcNum].ops[count].in[1] = deriv->prog[inst+3];
          inst += 4;
        }

        // initialize other data
        sys->funcs[funcNum].ops[count].lastUsed = -1;
        sys->funcs[funcNum].ops[count].numDiffInst = 0;
      }
    }
  }

  // put the subfunctions & functions into order - subfunctions first and then functions
  // this will simplify the derivative calculations
  for (i = 0; i < total; i++)
    if (i < sys->numSubfuncs)
    { // set to i
      sys->subFuncOrder[i] = i;
    }
    else
    { // set to i
      sys->funcOrder[i - sys->numSubfuncs] = i;
    }

  // setup other values
  sys->varGpSizes = NULL;

  return;
}

void setupDerivs_from_sys(funcStruct *derivs, funcStruct *subDerivs, int *totalOpCount, int **memLoc, int ***derivAddr, systemStruct *sys)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the first order partial derivatives              *
\***************************************************************/
{
  int i, j, currLoc;

  // setup totalOpCount
  *totalOpCount = sys->numVars + sys->numConstants + sys->numNumbers;
  // total number of things we need the derivative of
  for (i = 0; i < sys->numSubfuncs; i++)
    *totalOpCount += sys->subFuncs[i].num_ops;
  for (i = 0; i < sys->numFuncs; i++)
    *totalOpCount += sys->funcs[i].num_ops;
 
  // allocate memLoc & derivAddr 
  *memLoc = (int *)bmalloc(*totalOpCount * sizeof(int));
  *derivAddr = (int **)bmalloc(*totalOpCount * sizeof(int *));
  for (i = 0; i < *totalOpCount; i++)
  { // variables
    (*derivAddr)[i] = (int *)bmalloc(sys->numVars * sizeof(int));
    for (j = 0; j < sys->numVars; j++)
      (*derivAddr)[i][j] = -2;
  }

  // find the first free memory location based on the numbers, update, subfunction & function instructions
  sys->firstFreeMemLoc = sys->numNumbers + sys->numAddr;
  for (i = 0; i < sys->numUpdate; i++)
    sys->firstFreeMemLoc = MAX(sys->updateOps[i].memLoc, sys->firstFreeMemLoc);
  for (i = 0; i < sys->numSubfuncs; i++)
    for (j = 0; j < sys->subFuncs[i].num_ops; j++)
      sys->firstFreeMemLoc = MAX(sys->subFuncs[i].ops[j].memLoc, sys->firstFreeMemLoc);
  for (i = 0; i < sys->numFuncs; i++)
    for (j = 0; j < sys->funcs[i].num_ops; j++)
      sys->firstFreeMemLoc = MAX(sys->funcs[i].ops[j].memLoc, sys->firstFreeMemLoc);
  sys->firstFreeMemLoc++;

  // setup the subfunction derivatives
  currLoc = 0;
  for (i = 0; i < sys->numSubfuncs; i++)
  { // setup the ith subfunction partial derivatives
    diff_sys_vars_subfuncs_old(&subDerivs[currLoc], sys, i, *memLoc, *derivAddr, *totalOpCount);
    currLoc += sys->numVars;
  }
  
  // setup the function derivatives
  currLoc = 0;
  for (i = 0; i < sys->numFuncs; i++)
  { // setup the ith function partial derivatives
    diff_sys_vars_funcs_old(&derivs[currLoc], sys, i, *memLoc, *derivAddr, *totalOpCount);
    currLoc += sys->numVars;
  }

  return;
}

void setupNext_derivs(prog_deriv_t *deriv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the next set of partial derivatives              *
\***************************************************************/
{
  int i, j, k, l, count, currLoc, currFunc, currVar, newOpCount = deriv->totalOpCount;
  int *array = (int *)bmalloc(deriv->numVars * sizeof(int));

  // setup totalOpCount - add on the operations for the last set of derivatives
  i = deriv->order - 1;
  for (j = 0; j < deriv->numSubDerivs[i]; j++)
    newOpCount += deriv->subFuncDerivs[i][j].num_ops;
  for (j = 0; j < deriv->numDerivs[i]; j++)
    newOpCount += deriv->funcDerivs[i][j].num_ops;

  // increase the size of memLoc & derivAddr
  deriv->memLoc = (int *)brealloc(deriv->memLoc, newOpCount * sizeof(int));
  deriv->derivAddr = (int **)brealloc(deriv->derivAddr, newOpCount * sizeof(int *));
  for (i = deriv->totalOpCount; i < newOpCount; i++)
  {
    deriv->derivAddr[i] = (int *)bmalloc(deriv->numVars * sizeof(int));
    for (j = 0; j < deriv->numVars; j++)
      deriv->derivAddr[i][j] = -2;
  }

  // setup memLoc for the new operations
  count = deriv->totalOpCount;
  i = deriv->order - 1;
  for (j = 0; j < deriv->numSubDerivs[i]; j++)
    for (k = 0; k < deriv->subFuncDerivs[i][j].num_ops; k++)
    {
      deriv->memLoc[count] = deriv->subFuncDerivs[i][j].ops[k].memLoc;
      count++;
    }
  for (j = 0; j < deriv->numDerivs[i]; j++)
    for (k = 0; k < deriv->funcDerivs[i][j].num_ops; k++)
    {
      deriv->memLoc[count] = deriv->funcDerivs[i][j].ops[k].memLoc;
      count++;
    }

  // reallocate memory
  deriv->numDerivs = (int *)brealloc(deriv->numDerivs, (deriv->order + 1) * sizeof(int));
  deriv->numSubDerivs = (int *)brealloc(deriv->numSubDerivs, (deriv->order + 1) * sizeof(int));
  deriv->evalJFuncs = (int *)brealloc(deriv->evalJFuncs, (deriv->order + 1) * sizeof(int));
  deriv->evalJSubs = (int *)brealloc(deriv->evalJSubs, (deriv->order + 1) * sizeof(int));
  deriv->funcDerivs = (funcStruct **)brealloc(deriv->funcDerivs, (deriv->order + 1) * sizeof(funcStruct *));
  deriv->subFuncDerivs = (funcStruct **)brealloc(deriv->subFuncDerivs, (deriv->order + 1) * sizeof(funcStruct *));

  // setup the number of derivatives needed
  deriv->numDerivs[deriv->order] = deriv->numDerivs[deriv->order - 1] * (deriv->numVars + deriv->order) / (deriv->order + 1);
  deriv->numSubDerivs[deriv->order] = deriv->numSubDerivs[deriv->order - 1] * (deriv->numVars + deriv->order) / (deriv->order + 1);

  // allocate the subFuncDerivs & funcDerivs
  deriv->funcDerivs[deriv->order] = (funcStruct *)bmalloc(deriv->numDerivs[deriv->order] * sizeof(funcStruct));
  deriv->subFuncDerivs[deriv->order] = (funcStruct *)bmalloc(deriv->numSubDerivs[deriv->order] * sizeof(funcStruct));

  // setup where the derivatives are stored
  deriv->evalJFuncs[deriv->order] = deriv->sys.firstFreeMemLoc;
  deriv->sys.firstFreeMemLoc += deriv->numDerivs[deriv->order];
  deriv->evalJSubs[deriv->order] = deriv->sys.firstFreeMemLoc;
  deriv->sys.firstFreeMemLoc += deriv->numSubDerivs[deriv->order];

  // setup count to be the end of the other operations
  count = deriv->totalOpCount;

  // compute the derivatives for the subfunctions
  if (deriv->numSubfuncs > 0)
  { // loop over the subfunction derivatives to create the next order
    currLoc = 0;
    k = deriv->numSubDerivs[deriv->order - 1] / deriv->numSubfuncs;
    for (currFunc = 0; currFunc < deriv->numSubfuncs; currFunc++)
    { // reset array & currVar
      start_array(array, deriv->numVars, deriv->order);
      currVar = 0;

      for (j = 0; j < k; j++)
      { // need to compute the derivative of this from currVar to numVars
 
        // find the end of this function
        count += deriv->subFuncDerivs[deriv->order - 1][currFunc * k + j].num_ops;

        for (i = currVar; i < deriv->numVars; i++)
        {
          diff_funcStruct_old(&deriv->subFuncDerivs[deriv->order - 1][currFunc * k + j], count, i, deriv->evalJSubs[deriv->order] + currLoc, deriv->sys.numAddr, deriv->sys.numAddr + 1, &deriv->sys.firstFreeMemLoc, deriv->memLoc, deriv->derivAddr, newOpCount, &deriv->subFuncDerivs[deriv->order][currLoc].num_ops, &deriv->subFuncDerivs[deriv->order][currLoc].ops);

          // increment the deriv location
          currLoc++;
        }

        // move to the next setup
        advance_array(array, deriv->numVars);

        // find currVar - last non-zero entry in array
        for (l = deriv->numVars - 1; l >= 0; l--)
          if (array[l] != 0)
          { // this is the last non-zero entry - store it and exit loop
            currVar = l;
            break;
          }
      }
    }
  }

  // compute the derivatives for the functions
  currLoc = 0;
  k = deriv->numDerivs[deriv->order - 1] / deriv->numFuncs;
  for (currFunc = 0; currFunc < deriv->numFuncs; currFunc++)
  { // reset array & currVar
    start_array(array, deriv->numVars, deriv->order);
    currVar = 0;

    for (j = 0; j < k; j++)
    { // need to compute the derivative of this from currVar to numVars

      // find the end of this function
      count += deriv->funcDerivs[deriv->order - 1][currFunc * k + j].num_ops;

      for (i = currVar; i < deriv->numVars; i++)
      {
        diff_funcStruct_old(&deriv->funcDerivs[deriv->order - 1][currFunc * k + j], count, i, deriv->evalJFuncs[deriv->order] + currLoc, deriv->sys.numAddr, deriv->sys.numAddr + 1, &deriv->sys.firstFreeMemLoc, deriv->memLoc, deriv->derivAddr, newOpCount, &deriv->funcDerivs[deriv->order][currLoc].num_ops, &deriv->funcDerivs[deriv->order][currLoc].ops);

        // increment the deriv location
        currLoc++;
      }

      // move to the next setup
      advance_array(array, deriv->numVars);

      // find currVar - last non-zero entry in array
      for (l = deriv->numVars - 1; l >= 0; l--)
        if (array[l] != 0) 
        { // this is the last non-zero entry - store it and exit loop
          currVar = l; 
          break; 
        }
    }
  }

  // update the number of operations
  deriv->totalOpCount = newOpCount;

  // increment the order of the derivatives
  deriv->order++;

  // setup the SLP for this order
  setupProg_order_derivs(deriv);
  free(array);

  return;
}

void checkExp_derivs(prog_deriv_t *deriv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: see if the exponent operations can be simplified       *
\***************************************************************/
{
  int i, j, memSize = deriv->sys.firstFreeMemLoc, oneAddr = deriv->sys.numAddr + 1; // '1' is the second number
  int allowNegExp = 0;
  comp_d *mem = (comp_d *)bmalloc(memSize * sizeof(comp_d));

  // set mem to zero
  for (i = 0; i < memSize; i++)
    set_zero_d(mem[i]);

  // setup the numbers
  for (i = 0; i < deriv->sys.numNumbers; i++)
    mem[i + deriv->sys.numAddr]->r = mpq_get_d(deriv->nums[i].rat);

  // setup I
  set_double_d(mem[deriv->sys.constAddr], 0, 1);
  // setup Pi
  set_double_d(mem[deriv->sys.constAddr + 1], M_PI, 0);

  // setup variabes (random)
  for (i = 0; i < deriv->sys.numVars; i++)
    get_comp_rand_d(mem[i + deriv->sys.varsAddr]);

  // go through the update instructions
  checkExp_ops(deriv->sys.updateOps, deriv->sys.numUpdate, mem, oneAddr, allowNegExp);
  
  // go through the subfunctions
  for (i = 0; i < deriv->numSubfuncs; i++)
    checkExp_ops(deriv->sys.subFuncs[i].ops, deriv->sys.subFuncs[i].num_ops, mem, oneAddr, allowNegExp);

  // go through the functions
  for (i = 0; i < deriv->numFuncs; i++)
    checkExp_ops(deriv->sys.funcs[i].ops, deriv->sys.funcs[i].num_ops, mem, oneAddr, allowNegExp);

  // go through the derivatives
  for (j = 0; j < deriv->order; j++)
  { // go through the subfunction derivatives of order j
    for (i = 0; i < deriv->numSubDerivs[j]; i++)
      checkExp_ops(deriv->subFuncDerivs[j][i].ops, deriv->subFuncDerivs[j][i].num_ops, mem, oneAddr, allowNegExp);

    // go through the function derivatives of order j
    for (i = 0; i < deriv->numDerivs[j]; i++)
      checkExp_ops(deriv->funcDerivs[j][i].ops, deriv->funcDerivs[j][i].num_ops, mem, oneAddr, allowNegExp);
  }

  // clear mem
  free(mem);

  return;
}

void setupProg_order_derivs(prog_deriv_t *deriv)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the SLP in order for the current 'order'         *
\***************************************************************/
{
  int i, j, k, count = 0;
  int *prog = NULL;

  // check the exponents
  checkExp_derivs(deriv);

  // count the number of locations needed
  for (i = 0; i < deriv->sys.numUpdate; i++)
    if (isUnary(deriv->sys.updateOps[i].op))
      count += 3;
    else
      count += 4;
  for (j = 0; j < deriv->sys.numSubfuncs; j++)
    for (i = 0; i < deriv->sys.subFuncs[j].num_ops; i++)
      if (isUnary(deriv->sys.subFuncs[j].ops[i].op))
        count += 3;
      else 
        count += 4;
  for (j = 0; j < deriv->sys.numFuncs; j++)
    for (i = 0; i < deriv->sys.funcs[j].num_ops; i++)
      if (isUnary(deriv->sys.funcs[j].ops[i].op))
        count += 3;
      else
        count += 4;
  for (j = 0; j < deriv->order; j++)
    for (k = 0; k < deriv->numSubDerivs[j]; k++)
      for (i = 0; i < deriv->subFuncDerivs[j][k].num_ops; i++)
        if (isUnary(deriv->subFuncDerivs[j][k].ops[i].op))
          count += 3;
        else
          count += 4;
  for (j = 0; j < deriv->order; j++)
    for (k = 0; k < deriv->numDerivs[j]; k++)
      for (i = 0; i < deriv->funcDerivs[j][k].num_ops; i++)
        if (isUnary(deriv->funcDerivs[j][k].ops[i].op))
          count += 3;
        else
          count += 4;

  // allocate prog
  prog = (int *)bmalloc(count * sizeof(int));

  // setup prog
  count = 0;
  // update instructions
  for (i = 0; i < deriv->sys.numUpdate; i++)
  {
    prog[count] = deriv->sys.updateOps[i].op;
    prog[count+1] = deriv->sys.updateOps[i].memLoc;
    prog[count+2] = deriv->sys.updateOps[i].in[0];
    if (isUnary(prog[count]))
      count += 3;
    else
    {
      prog[count+3] = deriv->sys.updateOps[i].in[1];
      count += 4;
    }
  }

  // subfunctions
  for (j = 0; j < deriv->sys.numSubfuncs; j++)
    for (i = 0; i < deriv->sys.subFuncs[j].num_ops; i++)
    {
      prog[count] = deriv->sys.subFuncs[j].ops[i].op;
      prog[count+1] = deriv->sys.subFuncs[j].ops[i].memLoc;
      prog[count+2] = deriv->sys.subFuncs[j].ops[i].in[0];
      if (isUnary(prog[count]))
        count += 3;
      else
      {
        prog[count+3] = deriv->sys.subFuncs[j].ops[i].in[1];
        count += 4;
      }
    }

  // functions
  for (j = 0; j < deriv->sys.numFuncs; j++)
    for (i = 0; i < deriv->sys.funcs[j].num_ops; i++)
    {
      prog[count] = deriv->sys.funcs[j].ops[i].op;
      prog[count+1] = deriv->sys.funcs[j].ops[i].memLoc;
      prog[count+2] = deriv->sys.funcs[j].ops[i].in[0];
      if (isUnary(prog[count]))
        count += 3;
      else
      {
        prog[count+3] = deriv->sys.funcs[j].ops[i].in[1];
        count += 4;
      }
    }

  for (j = 0; j < deriv->order; j++)
  { // derivatives of subfunctions of order j
    for (k = 0; k < deriv->numSubDerivs[j]; k++)
      for (i = 0; i < deriv->subFuncDerivs[j][k].num_ops; i++)
      {
        prog[count] = deriv->subFuncDerivs[j][k].ops[i].op;
        prog[count+1] = deriv->subFuncDerivs[j][k].ops[i].memLoc;
        prog[count+2] = deriv->subFuncDerivs[j][k].ops[i].in[0];
        if (isUnary(prog[count]))
          count += 3; 
        else 
        { 
          prog[count+3] = deriv->subFuncDerivs[j][k].ops[i].in[1];
          count += 4;
        }
      }

    // derivatives of functions of order j
    for (k = 0; k < deriv->numDerivs[j]; k++)
      for (i = 0; i < deriv->funcDerivs[j][k].num_ops; i++)
      {
        prog[count] = deriv->funcDerivs[j][k].ops[i].op;
        prog[count+1] = deriv->funcDerivs[j][k].ops[i].memLoc;
        prog[count+2] = deriv->funcDerivs[j][k].ops[i].in[0];
        if (isUnary(prog[count]))
          count += 3;  
        else  
        {  
          prog[count+3] = deriv->funcDerivs[j][k].ops[i].in[1];
          count += 4;
        } 
      }
  }

  // copy prog to deriv->prog
  free(deriv->prog);
  deriv->prog = prog;
  deriv->size = count;
  deriv->memSizeNeeded = deriv->sys.firstFreeMemLoc;
  prog = NULL;

  return;
}

int change_prec_prog_deriv(void const *ED, int prec)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: changes the precision on prog_deriv                    *
 \***************************************************************/
{
  int i, j;
  prog_deriv_t *Prog = (prog_deriv_t *)ED;
	
  // change precision if needed
  if (Prog->precision != prec)
  { // update precision -- this will automatically change mem_mp during the next evaluation
    Prog->precision = prec;

    if (Prog->linearType == 2) 
    { // change linearCoeff_mp
      for (i = 0; i < Prog->linearCoeff_mp->rows; i++)
        for (j = 0; j < Prog->linearCoeff_mp->cols; j++)
        {
          setprec_mp(&Prog->linearCoeff_mp->entry[i][j], prec);
          rat_to_mp(&Prog->linearCoeff_mp->entry[i][j], Prog->linearCoeff_rat[i][j]);
        }

        // change linearConst_mp
	for (i = 0; i < Prog->linearConst_mp->size; i++)
        {
          setprec_mp(&Prog->linearConst_mp->coord[i], prec);
          rat_to_mp(&Prog->linearConst_mp->coord[i], Prog->linearConst_rat[i]);
        }
    }
  }

  Prog = NULL;

  return 0;
}



