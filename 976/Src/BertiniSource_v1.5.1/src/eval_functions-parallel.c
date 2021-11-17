// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

//The main eval function.
int eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
{
  int retVal = 0;
  retVal = eval_func(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ED);
  return retVal;
}

int patch_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{ // evaluate the patch

  patch_eval_data_d *P = (patch_eval_data_d *)ED;
  int i, j, num_patches = P->num_patches, size = vars->size;

  // setup the sizes
  change_size_vec_d(funcVals, num_patches);
  change_size_mat_d(Jv, num_patches, size);
  change_size_mat_d(Jp, num_patches, 1);
  funcVals->size = Jv->rows = Jp->rows = num_patches;
  Jv->cols = size;
  Jp->cols = 1;

  // evaluate the patch = coeff*vars - 1
  // Jv = patchCoeff
  for (i = 0; i < num_patches; i++)
  { // initialize funcVals to -1
    set_double_d(&funcVals->coord[i], -1.0, 0.0);
    for (j = 0; j < size; j++)
    { // funcVals += coeff * vars
      sum_mul_d(&funcVals->coord[i], &P->patchCoeff->entry[i][j], &vars->coord[j]);
      // Jv = coeff 
      set_d(&Jv->entry[i][j], &P->patchCoeff->entry[i][j]);
    }
    // Jp == 0
    set_zero_d(&Jp->entry[i][0]);
  }

  return 0;
}

int start_system_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{ // evaluate the start system

  start_system_eval_data_d *SS = (start_system_eval_data_d *)ED;
  int i, j, k, l, beg;
  comp_d tempComp, tempComp2;

  // initialize tempComp
  set_zero_d(tempComp);
  set_zero_d(tempComp2);

  // setup the sizes
  change_size_vec_d(funcVals, SS->size_r);
  change_size_mat_d(Jv, SS->size_r, vars->size);
  change_size_mat_d(Jp, SS->size_r, 1);
  funcVals->size = Jv->rows = Jp->rows = SS->size_r;
  Jv->cols = vars->size;
  Jp->cols = 1;

  if (SS->startSystemType == 0) // total degree start system is needed
  {
    for (i = 0; i < SS->size_r; i++) // find y_(i+1) ^ d_i - y_0 ^ d_i, i = 0, .. SS->size_r, where SS->size_r == vars->size - 1
    {
      if (SS->degrees[i] > 0)
      {
        exp_d(tempComp, &vars->coord[0], SS->degrees[i] - 1); // y_0 ^ (d_i - 1)
        exp_d(tempComp2, &vars->coord[i+1], SS->degrees[i] - 1); // y_(i+1) ^ (d_i - 1)
      }
      else
      {
        set_zero_d(tempComp);
        set_zero_d(tempComp2);
      }

      // zero out row of Jv since we are only storing to a couple of locations in this row
      for (j = 0; j < Jv->cols; j++)
      {
        set_zero_d(&Jv->entry[i][j]);
      }

      mul_rdouble_d(&Jv->entry[i][0], tempComp, -SS->degrees[i]); // -d_i * y_0 ^ (d_i - 1)
      mul_rdouble_d(&Jv->entry[i][i+1], tempComp2, SS->degrees[i]); // d_i * y_(i+1) ^ (d_i - 1)

      mul_d(tempComp, tempComp, &vars->coord[0]); // y_0 ^ d_i
      mul_d(tempComp2, tempComp2, &vars->coord[i+1]); // y_(i+1) ^ d_i

      sub_d(&funcVals->coord[i], tempComp2, tempComp);  // y_(i+1) ^ d_i - y_0 ^ d_i

      // Jp == 0
      set_zero_d(&Jp->entry[i][0]);
    }
  }
  else if (SS->startSystemType == 1) // product of linears is needed
  { // allocate space based on the maximum degree to store the linears
    comp_d *tempLinears = (comp_d *)bmalloc(SS->max_degree * sizeof(comp_d));

    beg = 0; // the beginning of the linears for the next function
    for (i = 0; i < SS->size_r; i++)
    {
      for (j = 0; j < SS->degrees[i]; j++)
      { // calculate the linear = sum(coeff * vars) - the linears do not have constant terms
        mul_d(tempLinears[j], &vars->coord[0], SS->coeff[beg + j][0]);
        for (k = 1; k < vars->size; k++)
        {
          sum_mul_d(tempLinears[j], &vars->coord[k], SS->coeff[beg + j][k]);
        }
      }

      // multiply the linears to find the start system value for the ith spot
      set_d(&funcVals->coord[i], tempLinears[0]);
      for (j = 1; j < SS->degrees[i]; j++)
      {
        mul_d(&funcVals->coord[i], &funcVals->coord[i], tempLinears[j]);
      }

      // calcuate the ith row of Jv = sum(prod(L_j, j != i) * coeff(L_i), i)
      for (j = 0; j < vars->size; j++)
      {
        set_zero_d(&Jv->entry[i][j]);
        for (k = 0; k < SS->degrees[i]; k++)
        {
          set_d(tempComp, SS->coeff[beg + k][j]);
          for (l = 0; l < k; l++)
          {
            mul_d(tempComp, tempComp, tempLinears[l]);
          }
          for (l = k + 1; l < SS->degrees[i]; l++)
          {
            mul_d(tempComp, tempComp, tempLinears[l]);
          }
          add_d(&Jv->entry[i][j], &Jv->entry[i][j], tempComp);
        }
      }
      beg += SS->degrees[i];

      // Jp == 0
      set_zero_d(&Jp->entry[i][0]);
    }

    // free tempLinears
    free(tempLinears);
  }

  return 0;
}

int square_system_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{ // evaluates the squared system F = [I A.*W]*[P]*f(By)

  square_system_eval_data_d *Sq = (square_system_eval_data_d *)ED; // to avoid having to cast every time

  if (Sq->noChanges)
  { // if the square system provides no changes, we just need to evaluate Prog
    evalProg_d_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, Sq->Prog);
  }
  else
  { // the square system provides changes, so we need to do the squaring procedure
    int i, j;
    comp_d tempComp, sum;
    mat_d mat_IAW, Jv_f, tempMat;
    vec_d f_old, f_new, x;

    // initialize
    init_mat_d(mat_IAW, 0, 0); init_mat_d(Jv_f, 0, 0); init_mat_d(tempMat, 0, 0);
    init_vec_d(f_old, 0); init_vec_d(f_new, 0); init_vec_d(x, 0);
    set_zero_d(tempComp);

    // calculate x = B*y
    mul_mat_vec_d(x, Sq->B, vars);

    // evalute the original function f at x
    evalProg_d_std(f_old, parVals, parDer, Jv_f, Jp, x, pathVars, Sq->Prog);

    // do the row swaps to f_old and store as f_new
    change_size_vec_d(f_new, f_old->size);
    f_new->size = f_old->size;
    for (i = 0; i < f_old->size; i++)
    {
      set_d(&f_new->coord[i], &f_old->coord[Sq->P[i]]); // f_(i) goes to f_(P[i])
    }

    // setup mat_IAW based on W & A
    change_size_mat_d(mat_IAW, Sq->size_r, f_old->size);
    mat_IAW->rows = Sq->size_r;
    mat_IAW->cols = f_old->size; // number of original functions f

    for (i = 0; i < Sq->size_r; i++)
    { // first part of the row is the 'identity row'
      for (j = 0; j < Sq->size_r; j++)
        if (i != j)
        {
          set_zero_d(&mat_IAW->entry[i][j]);
        }
        else
        {
          set_double_d(&mat_IAW->entry[i][j], 1.0, 0.0);
        }

      // last part of the row is the is A .* 'W' part
      for (j = Sq->size_r; j < f_old->size; j++)
      {
        exp_d(tempComp, &vars->coord[0], Sq->W[i][j-Sq->size_r]); // y_0 ^ W[i][j]
        mul_d(&mat_IAW->entry[i][j], &Sq->A->entry[i][j-Sq->size_r], tempComp); // A[i][j] * y_0 ^ W[i][j]
      }
    }

    // calculate F = [I A.*W][P]f(By) - this will set the size of funcVals appropriately
    mul_mat_vec_d(funcVals, mat_IAW, f_new);  // could be made more efficient based on structure of [I A.*W]

    // Jv evaluation
 
    // df/dy = df/dx*B
    mat_mul_d(Jv, Jv_f, Sq->B); // this should be an n x (r+1) matrix!

    // do the row swaps to find P*df/dy
    change_size_mat_d(tempMat, Jv->rows, Jv->cols);
    tempMat->rows = Jv->rows;
    tempMat->cols = Jv->cols;
    for (i = 0; i < Jv->rows; i++)
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&tempMat->entry[i][j], &Jv->entry[Sq->P[i]][j]);  // f_(i) goes to f_(P[i])
      }

    // multiply these matrices
    mat_mul_d(Jv, mat_IAW, tempMat); // [I A.* W][P]*df/dx(By)*B - should be an r x (r+1) matrix!
  
    // the y_0 variables in W also affect Jv, but only in the first column - the column associated with y_0
    change_size_mat_d(Jp, Sq->size_r, 1);
    Jp->rows = Sq->size_r;
    Jp->cols = 1;
    for (i = 0; i < Sq->size_r; i++)
    { // calculate the (i,0) entry of [0 A.*dW/dy][P]*f(By) = sum(A[i][j]*(W[i][j] * y_0 ^ (W[i][j] - 1)) * f_new[j], j = r+1, .. , n)
      set_zero_d(sum);
      for (j = Sq->size_r; j < f_new->size; j++)
      {
        if (Sq->W[i][j - Sq->size_r] > 0)
        {
          exp_d(tempComp, &vars->coord[0], Sq->W[i][j-Sq->size_r] - 1); // y_0 ^ (W[i][j] - 1)
          mul_rdouble_d(tempComp, tempComp, Sq->W[i][j-Sq->size_r]);    // W[i][j] * y_0 ^ (W[i][j] - 1)
        }
        else
        {
          set_zero_d(tempComp);
        }

        mul_d(tempComp, tempComp, &Sq->A->entry[i][j-Sq->size_r]); // A[i][j] * W[i][j] * y_0 ^ (W[i][j] - 1)
        sum_mul_d(sum, tempComp, &f_new->coord[j]); // sum += f_j' * A[i][j] * W[i][j] * y_0 ^ (W[i][j] - 1)
      }
      // add on this value to Jv in the appropriate spot
      add_d(&Jv->entry[i][0], &Jv->entry[i][0], sum);

      // Jp == 0
      set_zero_d(&Jp->entry[i][0]);
    }

    clear_mat_d(mat_IAW); clear_mat_d(Jv_f); clear_mat_d(tempMat);
    clear_vec_d(f_old); clear_vec_d(f_new); clear_vec_d(x);
  }

  return 0;
}

int basic_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{ // evaluates the function as described in basic_eval_data_d: 'top' squareSystem*s + gamma*(1-s)*startSystem and 'bottom' patch  

  basic_eval_data_d *BED = (basic_eval_data_d *)ED; // to avoid having to cast every time
  int i, j, size, top;
  comp_d one_minus_s, gamma_s;
  vec_d F, startSys, patchValues;
  mat_d Jv_F, Jv_start, Jv_Patch; 
  // we assume that the only parameter is s = t and setup parVals & parDer accordingly. Also, F, startSys & patchValues do not depend on a parameter

  // initialize
  init_vec_d(F, 0); init_vec_d(startSys, 0); init_vec_d(patchValues, 0);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_start, 0, 0); init_mat_d(Jv_Patch, 0, 0); 
 
  set_one_d(one_minus_s);
  sub_d(one_minus_s, one_minus_s, pathVars);
  mul_d(gamma_s, BED->startSystem.gamma, pathVars);

  // evalute the 'square' system F
  square_system_eval_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, &BED->squareSystem); // Jp is ignored
  // evaluate the patch
  patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);  // Jp is ignored
  // evaluate the start system
  start_system_eval_d(startSys, parVals, parDer, Jv_start, Jp, vars, pathVars, &BED->startSystem); // Jp is ignored

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1

  // combine everything

  // find funcVals = 'top' F*(1-s) + gamma*s*startSys and 'bottom' patchValues
  // find Jv = 'top' Jv_F*(1-s) + gamma*s*Jv_start and 'bottom' Jv_Patch
  // find Jp = 'top' -F + gamma*startSys and 'bottom' all zeros
  size = F->size + BED->patch.num_patches;
  change_size_mat_d(Jv, size, vars->size);
  change_size_mat_d(Jp, size, 1);
  change_size_vec_d(funcVals, size); 
  funcVals->size = Jv->rows = Jp->rows = size;
  top = F->size;
  Jv->cols = vars->size;  // this should be square!!!
  Jp->cols = 1;
  for (i = 0; i < size; i++)
  {
    if (i < top)
    { // funcVals = F*(1-s) + gamma*s*startSys
      mul_d(&funcVals->coord[i], &F->coord[i], one_minus_s);
      sum_mul_d(&funcVals->coord[i], &startSys->coord[i], gamma_s);
      // Jp = -F + gamma*startSys
      mul_d(&Jp->entry[i][0], BED->startSystem.gamma, &startSys->coord[i]);
      sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &F->coord[i]);
      // Jv = Jv_F*(1-s) + gamma*s*Jv_start
      for (j = 0; j < vars->size; j++)
      {
        mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_s);
        sum_mul_d(&Jv->entry[i][j], &Jv_start->entry[i][j], gamma_s);
      }
    }
    else
    { // funcVals = patchValues
      set_d(&funcVals->coord[i], &patchValues->coord[i - top]);
      // Jp = 0
      set_zero_d(&Jp->entry[i][0]);
      // Jv = Jv_Patch
      for (j = 0; j < vars->size; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_Patch->entry[i - top][j]);
      }
    }
  }

  clear_mat_d(Jv_F); clear_mat_d(Jv_start); clear_mat_d(Jv_Patch);
  clear_vec_d(F); clear_vec_d(startSys); clear_vec_d(patchValues);

  return 0;
}

int userHom_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{
  int retVal = 0;
  retVal = evalProg_d_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ((basic_eval_data_d *)ED)->squareSystem.Prog);
  return retVal;
}

int paramHom_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
{ // SLP + patches
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;
  int i, j, numPatch, numFuncs, numVars, numParams;
  vec_d patchValues;
  mat_d Jv_Patch;

  init_vec_d(patchValues, 0);
  init_mat_d(Jv_Patch, 0, 0);

  // evaluate the patch
  patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);

  numPatch = patchValues->size;
  numVars = Jv_Patch->cols;

  // evaluate the SLP
  evalProg_d_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, BED->squareSystem.Prog);
  numFuncs = funcVals->size;
  numParams = Jp->cols;

  // setup funcVals, Jv & Jp correctly
  increase_size_vec_d(funcVals, numFuncs + numPatch);
  increase_size_mat_d(Jv, numFuncs + numPatch, numVars);
  increase_size_mat_d(Jp, numFuncs + numPatch, numParams);
  funcVals->size = Jv->rows = Jp->rows = numFuncs + numPatch;
  Jv->cols = numVars;
  Jp->cols = numParams;
  for (i = 0; i < numPatch; i++)
  { // setup funcVals
    set_d(&funcVals->coord[i + numFuncs], &patchValues->coord[i]);
    for (j = 0; j < numVars; j++)
    { // setup Jv
      set_d(&Jv->entry[i + numFuncs][j], &Jv_Patch->entry[i][j]);
    }
    for (j = 0; j < numParams; j++)
    { // setup Jp
      set_zero_d(&Jp->entry[i + numFuncs][j]);
    }
  }

  clear_vec_d(patchValues);
  clear_mat_d(Jv_Patch);

  return 0;
}

////// THE MP VERSIONS ///////

int eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
{
  int retVal = 0;
  retVal = eval_func(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ED);
  return retVal;
}

int patch_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{ // evaluate the patch

  patch_eval_data_mp *P = (patch_eval_data_mp *)ED;
  int i, j, num_patches = P->num_patches, size = vars->size;

  // setup the sizes
  change_size_vec_mp(funcVals, num_patches);
  change_size_mat_mp(Jv, num_patches, size);
  change_size_mat_mp(Jp, num_patches, 1);
  funcVals->size = Jv->rows = Jp->rows = num_patches;
  Jv->cols = size;
  Jp->cols = 1;

  // evaluate the patch = coeff*vars - 1
  // Jv = patchCoeff
  for (i = 0; i < num_patches; i++)
  { // initialize funcVals to -1
    set_one_mp(&funcVals->coord[i]);
    neg_mp(&funcVals->coord[i], &funcVals->coord[i]);
    for (j = 0; j < size; j++)
    { // funcVals += coeff * vars
      sum_mul_mp(&funcVals->coord[i], &P->patchCoeff->entry[i][j], &vars->coord[j]);
      // Jv = coeff
      set_mp(&Jv->entry[i][j], &P->patchCoeff->entry[i][j]);
    }
    // Jp == 0
    set_zero_mp(&Jp->entry[i][0]);
  }

  return 0;
}

int start_system_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{ // evaluate the start system
 
  start_system_eval_data_mp *SS = (start_system_eval_data_mp *)ED;
  int i, j, k, l, beg;
  comp_mp tempComp1, tempComp2;

  init_mp(tempComp1); init_mp(tempComp2);

  // setup the sizes
  change_size_vec_mp(funcVals, SS->size_r);
  change_size_mat_mp(Jv, SS->size_r, vars->size);
  change_size_mat_mp(Jp, SS->size_r, 1);
  funcVals->size = Jv->rows = Jp->rows = SS->size_r;
  Jv->cols = vars->size;
  Jp->cols = 1;
 
  if (SS->startSystemType == 0) // total degree start system is needed
  {
    for (i = 0; i < SS->size_r; i++) // find y_(i+1) ^ d_i - y_0 ^ d_i, i = 0, .. SS->size_r, where SS->size_r == vars->size - 1
    {
      if (SS->degrees[i] > 0)
      {
        exp_mp_int(tempComp1, &vars->coord[0], SS->degrees[i] - 1); // y_0 ^ (d_i - 1)
        exp_mp_int(tempComp2, &vars->coord[i+1], SS->degrees[i] - 1); // y_(i+1) ^ (d_i - 1)
      }
      else
      {
        set_zero_mp(tempComp1);
        set_zero_mp(tempComp2);
      }
 
      // zero out row of Jv since we are only storing to a couple of locations in this row
      for (j = 0; j < Jv->cols; j++)
      {
        set_zero_mp(&Jv->entry[i][j]);
      }

      mul_rdouble_mp(&Jv->entry[i][0], tempComp1, -SS->degrees[i]); // -d_i * y_0 ^ (d_i - 1)
      mul_rdouble_mp(&Jv->entry[i][i+1], tempComp2, SS->degrees[i]); // d_i * y_(i+1) ^ (d_i - 1)
 
      mul_mp(tempComp1, tempComp1, &vars->coord[0]); // y_0 ^ d_i
      mul_mp(tempComp2, tempComp2, &vars->coord[i+1]); // y_(i+1) ^ d_i
 
      sub_mp(&funcVals->coord[i], tempComp2, tempComp1);  // y_(i+1) ^ d_i - y_0 ^ d_i

      // Jp == 0
      set_zero_mp(&Jp->entry[i][0]);
    }
  }
  else if (SS->startSystemType == 1) // product of linears is needed
  { // allocate space based on the maximum degree to store the linears
    comp_mp *tempLinears = (comp_mp *)bmalloc(SS->max_degree * sizeof(comp_mp));

    for (i = 0; i < SS->max_degree; i++)
      init_mp(tempLinears[i]);

    beg = 0; // the beginning of the linears for the next function
    for (i = 0; i < SS->size_r; i++)
    {
      for (j = 0; j < SS->degrees[i]; j++)
      { // calculate the linear = sum(coeff * vars) - the linears do not have constant terms
        mul_mp(tempLinears[j], &vars->coord[0], SS->coeff[beg + j][0]);
        for (k = 1; k < vars->size; k++)
        {
          sum_mul_mp(tempLinears[j], &vars->coord[k], SS->coeff[beg + j][k]);
        }
      }
 
      // multiply the linears to find the start system value for the ith spot
      set_mp(&funcVals->coord[i], tempLinears[0]);
      for (j = 1; j < SS->degrees[i]; j++)
      {
        mul_mp(&funcVals->coord[i], &funcVals->coord[i], tempLinears[j]);
      }

      // calcuate the ith row of Jv = sum(prod(L_j, j != i) * coeff(L_i), i)
      for (j = 0; j < vars->size; j++)
      {
        set_zero_mp(&Jv->entry[i][j]);
        for (k = 0; k < SS->degrees[i]; k++)
        {
          if (d_oneNorm_mp(SS->coeff[beg+k][j]) > 0)
          { // only need to multiply if this coefficient is non-zero
            set_mp(tempComp1, SS->coeff[beg + k][j]);
            for (l = 0; l < k; l++)
            {
              mul_mp(tempComp1, tempComp1, tempLinears[l]);
            }
            for (l = k + 1; l < SS->degrees[i]; l++)
            {
              mul_mp(tempComp1, tempComp1, tempLinears[l]);
            }
            add_mp(&Jv->entry[i][j], &Jv->entry[i][j], tempComp1);
          }
        }
      }
      beg += SS->degrees[i];

      // Jp == 0
      set_zero_mp(&Jp->entry[i][0]);
    }

    // clear tempLinears
    for (i = SS->max_degree - 1; i >= 0; i--)
    {
      clear_mp(tempLinears[i]);
    }
    free(tempLinears);
  }

  clear_mp(tempComp1); clear_mp(tempComp2);

  return 0;
}

int square_system_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{ // evaluates the squared system F = [I A.*W]*[P]*f(By)

  square_system_eval_data_mp *Sq = (square_system_eval_data_mp *)ED; // to avoid having to cast every time

  if (Sq->noChanges)
  { // if the square system provides no changes, we just need to evaluate Prog
    evalProg_mp_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, Sq->Prog);
  }
  else
  { // the square system provides changes, so we need to do the squaring procedure
    int i, j;
    comp_mp tempComp1, tempComp2;
    vec_mp tempVec1, tempVec2;
    mat_mp tempMat1, tempMat2;

    init_mp(tempComp1); init_mp(tempComp2);
    init_vec_mp(tempVec1, 0); init_vec_mp(tempVec2, 0);
    init_mat_mp(tempMat1, 0, 0); init_mat_mp(tempMat2, 0, 0);

    // calculate x = B*y
    mul_mat_vec_mp(tempVec2, Sq->B, vars);

    // evalute the original function f at x
    evalProg_mp_std(tempVec1, parVals, parDer, tempMat2, Jp, tempVec2, pathVars, Sq->Prog);

    // do the row swaps to f_old and store as f_new
    change_size_vec_mp(tempVec2, tempVec1->size);
    tempVec2->size = tempVec1->size;
    for (i = 0; i < tempVec2->size; i++)
    {
      set_mp(&tempVec2->coord[i], &tempVec1->coord[Sq->P[i]]); // f_(i) goes to f_(P[i])
    }

    // setup mat_IAW based on W & A
    change_size_mat_mp(tempMat1, Sq->size_r, tempVec2->size);
    tempMat1->rows = Sq->size_r;
    tempMat1->cols = tempVec2->size; // number of original functions f

    for (i = 0; i < Sq->size_r; i++)
    { // first part of the row is the 'identity row'
      for (j = 0; j < Sq->size_r; j++)
        if (i != j)
        {
          set_zero_mp(&tempMat1->entry[i][j]);
        }
        else
        {
          set_double_mp(&tempMat1->entry[i][j], 1.0, 0.0);
        }

      // last part of the row is the is A .* 'W' part
      for (j = Sq->size_r; j < tempVec2->size; j++)
      {
        exp_mp_int(tempComp1, &vars->coord[0], Sq->W[i][j-Sq->size_r]); // y_0 ^ W[i][j]
        mul_mp(&tempMat1->entry[i][j], &Sq->A->entry[i][j-Sq->size_r], tempComp1); // A[i][j] * y_0 ^ W[i][j]
      }
    }

    // calculate F = [I A.*W][P]f(By) - this will set the size of funcVals appropriately
    mul_mat_vec_mp(funcVals, tempMat1, tempVec2);  // could be made more efficient based on structure of [I A.*W]
  
    // Jv evaluation
 
    // df/dy = df/dx*B
    mat_mul_mp(Jv, tempMat2, Sq->B); // this should be an n x (r+1) matrix!

    // do the row swaps to find P*df/dy
    change_size_mat_mp(tempMat2, Jv->rows, Jv->cols);
    tempMat2->rows = Jv->rows;
    tempMat2->cols = Jv->cols;
    for (i = 0; i < Jv->rows; i++)
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&tempMat2->entry[i][j], &Jv->entry[Sq->P[i]][j]);  // f_(i) goes to f_(P[i])
      }

    // multiply these matrices
    mat_mul_mp(Jv, tempMat1, tempMat2); // [I A.* W][P]*df/dx(By)*B - should be an r x (r+1) matrix!
 
    // the y_0 variables in W also affect Jv, but only in the first column - the column associated with y_0
    change_size_mat_mp(Jp, Sq->size_r, 1);
    Jp->rows = Sq->size_r;
    Jp->cols = 1;
    for (i = 0; i < Sq->size_r; i++)
    { // calculate the (i,0) entry of [0 A.*dW/dy][P]*f(By) = sum(A[i][j]*(W[i][j] * y_0 ^ (W[i][j] - 1)) * f_new[j], j = r+1, .. , n)
      set_zero_mp(tempComp2);
      for (j = Sq->size_r; j < tempVec2->size; j++)
      {
        if (Sq->W[i][j - Sq->size_r] > 0)
        {
          exp_mp_int(tempComp1, &vars->coord[0], Sq->W[i][j-Sq->size_r] - 1); // y_0 ^ (W[i][j] - 1)
          mul_rdouble_mp(tempComp1, tempComp1, Sq->W[i][j-Sq->size_r]);    // W[i][j] * y_0 ^ (W[i][j] - 1)
        }
        else
        {
          set_zero_mp(tempComp1);
        }
 
        mul_mp(tempComp1, tempComp1, &Sq->A->entry[i][j-Sq->size_r]); // A[i][j] * W[i][j] * y_0 ^ (W[i][j] - 1)
        sum_mul_mp(tempComp2, tempComp1, &tempVec2->coord[j]); // sum += f_j' * A[i][j] * W[i][j] * y_0 ^ (W[i][j] - 1)
      }
      // add on this value to Jv in the appropriate spot
      add_mp(&Jv->entry[i][0], &Jv->entry[i][0], tempComp2);

      // Jp == 0
      set_zero_mp(&Jp->entry[i][0]);
    }

    clear_mp(tempComp1); clear_mp(tempComp2);
    clear_vec_mp(tempVec1); clear_vec_mp(tempVec2);
    clear_mat_mp(tempMat1); clear_mat_mp(tempMat2);
  }

  return 0;
}

int basic_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{ // evaluates the function as described in basic_eval_data_mp: 'top' squareSystem*s + gamma*(1-s)*startSystem and 'bottom' patch

  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED; // to avoid having to cast every time
  int i, j, size, top;
  // we assume that the only parameter is s = t and setup parVals & parDer accordingly. Also, F, startSys & patchValues do not depend on a parameter

  comp_mp one_minus_s, gamma_s;
  vec_mp F, startSys, patchValues;
  mat_mp Jv_F, Jv_start, Jv_Patch;

  // initialize
  init_mp(one_minus_s); init_mp(gamma_s);
  init_vec_mp(F, 0); init_vec_mp(startSys, 0); init_vec_mp(patchValues, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_start, 0, 0); init_mat_mp(Jv_Patch, 0, 0);
  
  set_one_mp(one_minus_s);
  sub_mp(one_minus_s, one_minus_s, pathVars); // 1 - s
  mul_mp(gamma_s, BED->startSystem.gamma, pathVars); // gamma * s

  // evalute the 'square' system F
  square_system_eval_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, &BED->squareSystem); // Jp is ignored
  // evaluate the patch
  patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);  // Jp is ignored
  // evaluate the start system
  start_system_eval_mp(startSys, parVals, parDer, Jv_start, Jp, vars, pathVars, &BED->startSystem); // Jp is ignored

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_point_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);       // ds/dt = 1

  // combine everything 

  // find funcVals = 'top' F*(1-s) + gamma*s*startSys and 'bottom' patchValues
  // find Jv = 'top' Jv_F*(1-s) + gamma*s*Jv_start and 'bottom' Jv_Patch
  // find Jp = 'top' -F + gamma*startSys and 'bottom' all zeros
  size = F->size + BED->patch.num_patches;
  change_size_mat_mp(Jv, size, vars->size);
  change_size_mat_mp(Jp, size, 1);
  change_size_vec_mp(funcVals, size);
  funcVals->size = Jv->rows = Jp->rows = size;
  top = F->size;
  Jv->cols = vars->size;  // this should be square!!!
  Jp->cols = 1;
  for (i = 0; i < size; i++)
  { 
    if (i < top)
    { // funcVals = F*(1-s) + gamma*s*startSys
      mul_mp(&funcVals->coord[i], &F->coord[i], one_minus_s);
      sum_mul_mp(&funcVals->coord[i], &startSys->coord[i], gamma_s);
      // Jp = -F + gamma*startSys
      mul_mp(&Jp->entry[i][0], BED->startSystem.gamma, &startSys->coord[i]);
      sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], &F->coord[i]);
      // Jv = Jv_F*(1-s) + gamma*s*Jv_start
      for (j = 0; j < vars->size; j++)
      {
        mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_s);
        sum_mul_mp(&Jv->entry[i][j], &Jv_start->entry[i][j], gamma_s);
      }
    }
    else
    { // funcVals = patchValues
      set_mp(&funcVals->coord[i], &patchValues->coord[i - top]);
      // Jp = 0
      set_zero_mp(&Jp->entry[i][0]);
      // Jv = Jv_Patch  
      for (j = 0; j < vars->size; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_Patch->entry[i - top][j]);
      }
    }
  }

  clear_mp(one_minus_s); clear_mp(gamma_s);
  clear_vec_mp(F); clear_vec_mp(startSys); clear_vec_mp(patchValues);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_start); clear_mat_mp(Jv_Patch); 

  return 0;
}

int userHom_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{
  int retVal = 0;
  retVal = evalProg_mp_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ((basic_eval_data_mp *)ED)->squareSystem.Prog);
  return retVal;
}

int paramHom_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
{ // SLP + patches
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;
  int i, j, numPatch, numFuncs, numVars, numParams;
  vec_mp patchValues;
  mat_mp Jv_Patch;

  init_vec_mp(patchValues, 0);
  init_mat_mp(Jv_Patch, 0, 0);

  // evaluate the patch
  patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);

  numPatch = patchValues->size;
  numVars = Jv_Patch->cols;

  // evaluate the SLP
  evalProg_mp_std(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, BED->squareSystem.Prog);

  numFuncs = funcVals->size;
  numParams = Jp->cols;

  // setup funcVals, Jv & Jp correctly
  increase_size_vec_mp(funcVals, numFuncs + numPatch);
  increase_size_mat_mp(Jv, numFuncs + numPatch, numVars);
  increase_size_mat_mp(Jp, numFuncs + numPatch, numParams);
  funcVals->size = Jv->rows = Jp->rows = numFuncs + numPatch;
  Jv->cols = numVars;
  Jp->cols = numParams;
  for (i = 0; i < numPatch; i++)
  { // setup funcVals
    set_mp(&funcVals->coord[i + numFuncs], &patchValues->coord[i]);
    for (j = 0; j < numVars; j++)
    { // setup Jv
      set_mp(&Jv->entry[i + numFuncs][j], &Jv_Patch->entry[i][j]);
    }
    for (j = 0; j < numParams; j++)
    { // setup Jp
      set_zero_mp(&Jp->entry[i + numFuncs][j]);
    }
  }

  clear_vec_mp(patchValues);
  clear_mat_mp(Jv_Patch);

  return 0;
}


