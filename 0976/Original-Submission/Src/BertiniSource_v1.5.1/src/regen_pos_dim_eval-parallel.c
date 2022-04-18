// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "regen_pos_dim.h"

void topID_mul_mat_vec_d(vec_d Res, mat_d A, vec_d X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume the top of A is the Identity matrix             *
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;

  if (A->cols != X->size)
    printf("WARNING:  Attempting to multiply a matrix (%d x %d) with a vector (%d) in which the dimensions do not match!\n", rows, cols, X->size);

  // make sure we have enough room
  increase_size_vec_d(Res, rows);

  // compute the bottom of Res
  for (i = 0; i < rows; i++)
    if (i < cols)
    { // Res[i] = X[i]
      set_d(&Res->coord[i], &X->coord[i]);
    }
    else
    { // Res[i] = SUM(A[i][j] * X[j])
      set_zero_d(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      { 
        sum_mul_d(&Res->coord[i], &A->entry[i][j], &X->coord[j]); 
      }
    }

  // set the size correctly
  change_size_vec_d(Res, rows);
  Res->size = rows;

  return;
}

void topID_mat_mul_d(mat_d Res, mat_d A, mat_d B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume the top of B is the Identity matrix and Res != B*
\***************************************************************/
{
  int i, j, k, rows = A->rows, cols = B->cols, inner = A->cols, bottom = B->cols;

  if (A->cols != B->rows)
    printf("WARNING:  Attempting to multiply matrices with dimensions that do not match!\n");
  else if (Res == B)
  {
    printf("ERROR: The right-hand matrix cannot be the same as the output matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // make sure we have enough room
  increase_size_mat_d(Res, rows, cols);

  // do the randomizing upward
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    { // Res[i][j] = A[i][j] + ...
      set_d(&Res->entry[i][j], &A->entry[i][j]);
      for (k = bottom; k < inner; k++)
      {
        sum_mul_d(&Res->entry[i][j], &A->entry[i][k], &B->entry[k][j]); // Res_i,j += A_i,k * B_k,j
      }
    }

  // set the size correctly
  change_size_mat_d(Res, rows, cols);
  Res->rows = rows;
  Res->cols = cols;

  return;
}

void topID_mul_mat_vec_mp(vec_mp Res, mat_mp A, vec_mp X)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume the top of A is the Identity matrix             *
\***************************************************************/
{
  int i, j, rows = A->rows, cols = A->cols;

  if (A->cols != X->size)
    printf("WARNING:  Attempting to multiply a matrix (%d x %d) with a vector (%d) in which the dimensions do not match!\n", rows, cols, X->size);

  // make sure we have enough room
  increase_size_vec_mp(Res, rows);

  // compute the bottom of Res
  for (i = 0; i < rows; i++)
    if (i < cols)
    { // Res[i] = X[i]
      set_mp(&Res->coord[i], &X->coord[i]);
    }
    else
    { // Res[i] = SUM(A[i][j] * X[j])
      set_zero_mp(&Res->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Res->coord[i], &A->entry[i][j], &X->coord[j]);
      }
    }

  // set the size correctly
  change_size_vec_mp(Res, rows);
  Res->size = rows;

  return;
}

void topID_mat_mul_mp(mat_mp Res, mat_mp A, mat_mp B)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Assume the top of B is the Identity matrix and Res != B*
\***************************************************************/
{
  int i, j, k, rows = A->rows, cols = B->cols, inner = A->cols, bottom = B->cols;

  if (A->cols != B->rows)
    printf("WARNING:  Attempting to multiply matrices with dimensions that do not match!\n");
  else if (Res == B)
  {
    printf("ERROR: The right-hand matrix cannot be the same as the output matrix!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // make sure we have enough room
  increase_size_mat_mp(Res, rows, cols);

  // do the randomizing upward
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
    { // Res[i][j] = A[i][j] + ...
      set_mp(&Res->entry[i][j], &A->entry[i][j]);
      for (k = bottom; k < inner; k++)
      {
        sum_mul_mp(&Res->entry[i][j], &A->entry[i][k], &B->entry[k][j]); // Res_i,j += A_i,k * B_k,j
      }
    }

  // set the size correctly
  change_size_mat_mp(Res, rows, cols);
  Res->rows = rows;
  Res->cols = cols;

  return;
}

void regen_pos_dim_square_d(vec_d F_sqr, mat_d Jv_sqr, point_d vars, point_d F, mat_d Jv_F, int codim_index, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Turn the system into the square system for the codim   *
\***************************************************************/
{ // F_sqr = [A*.W] F
  int i, j, k, exp, F_size, cols, codim = RPD->codim[codim_index].codim, num_vars = vars->size;
  comp_d tempComp, hom_var;
  vec_d F_new;
  mat_d Jv_F_new, sqrMat, sqrMat_deriv;

  // initialize
  set_zero_d(tempComp); set_zero_d(hom_var);
  F_size = F->size;
  cols = Jv_F->cols;

  init_vec_d(F_new, F_size);
  F_new->size = F_size;

  init_mat_d(Jv_F_new, F_size, cols);
  Jv_F_new->rows = F_size;
  Jv_F_new->cols = cols;

  init_mat_d(sqrMat, codim, F_size);
  init_mat_d(sqrMat_deriv, codim, F_size);
  sqrMat->rows = sqrMat_deriv->rows = codim;
  sqrMat->cols = sqrMat_deriv->cols = F_size;

  // reorder F to F_new and Jv_F to Jv_F_new
  for (i = 0; i < F_size; i++)
  { // do the swap on F
    k = RPD->P[i];
    set_d(&F_new->coord[i], &F->coord[k]);
    for (j = 0; j < cols; j++)
    { // do the swap of Jv_F
      set_d(&Jv_F_new->entry[i][j], &Jv_F->entry[k][j]);
    }
  }

  // compute hom_var
  set_d(hom_var, RPD->homVarConst_d);
  for (i = 0; i < num_vars; i++)
  {
    sum_mul_d(hom_var, &vars->coord[i], &RPD->H_d->coord[i]);
  }

  // setup sqrMat  = [A * h^W] & sqrMat_deriv = [W * A * h^(W-1)]
  for (i = 0; i < codim; i++)
    for (j = 0; j < F_size; j++)
      if (i > j)
      { // set to 0
        set_zero_d(&sqrMat->entry[i][j]);
        set_zero_d(&sqrMat_deriv->entry[i][j]);
      }
      else if (i == j)
      { // set to A
        set_d(&sqrMat->entry[i][j], &RPD->A_d[codim_index]->entry[i][j]);
        // set deriv to 0 
        set_zero_d(&sqrMat_deriv->entry[i][j]);
      }
      else // i < j
      { // compute A * h^W & set deriv to W * A * h^(W-1) & 
        exp = RPD->W[codim_index][i][j - i - 1];
        if (exp > 0)
        { // compute h^(W-1)
          exp_d(tempComp, hom_var, exp - 1);
          // compute A * h^(W-1)
          mul_d(tempComp, tempComp, &RPD->A_d[codim_index]->entry[i][j]);
          // compute sqrMat = A * h^W
          mul_d(&sqrMat->entry[i][j], tempComp, hom_var);
          // compute deriv = W * A * h^(W-1)
          mul_rdouble_d(&sqrMat_deriv->entry[i][j], tempComp, exp);
        }
        else
        { // h^0 = 1 
          set_d(&sqrMat->entry[i][j], &RPD->A_d[codim_index]->entry[i][j]);
          set_zero_d(&sqrMat_deriv->entry[i][j]);
        }
      }

  // compute F_sqr = [A * h^W] * F_new & Jv_sqr = d([A * h^W] * F_new) = [A * h^W] * Jv_F_new + [W * A * h^(W-1)] * F_new * d(h)
  change_size_vec_d(F_sqr, codim);
  change_size_mat_d(Jv_sqr, codim, cols);
  F_sqr->size = Jv_sqr->rows = codim;
  Jv_sqr->cols = cols;
  for (i = 0; i < codim; i++)
  { // compute F_sqr[i]
    set_zero_d(&F_sqr->coord[i]);
    for (j = i; j < F_size; j++)
    { // += sqrMat * F_new
      sum_mul_d(&F_sqr->coord[i], &sqrMat->entry[i][j], &F_new->coord[j]);
    }

    for (k = 0; k < cols; k++)
    { // compute Jv_F_sqr[i][k]
      mul_d(&Jv_sqr->entry[i][k], &sqrMat->entry[i][i], &Jv_F_new->entry[i][k]);

      for (j = i + 1; j < F_size; j++)
      { // += sqrMat * Jv_F_new
        sum_mul_d(&Jv_sqr->entry[i][k], &sqrMat->entry[i][j], &Jv_F_new->entry[j][k]);

        // += d(h^W) * F_new
        mul_d(tempComp, &sqrMat_deriv->entry[i][j], &F_new->coord[j]);
        sum_mul_d(&Jv_sqr->entry[i][k], tempComp, &RPD->H_d->coord[k]);
      }
    }
  }

  clear_vec_d(F_new);
  clear_mat_d(Jv_F_new); clear_mat_d(sqrMat); clear_mat_d(sqrMat_deriv);

  return;
}

void regen_pos_dim_square_old_new_d(vec_d top_F, mat_d Jv_top_F, mat_d Jp_top_F, comp_d extra_F, vec_d Jv_extra_F, vec_d vars, comp_d pathVars, vec_d F, mat_d Jv_F, int codim_index, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computes top_F = [old_F]*gamma*t + one_minus_t*[new_F] *
*   and extra_F to be the next function being added             *
\***************************************************************/
{
  int i, j, topSize;

  if (RPD->sameA)
  { // we have the same A so that we can be efficient in the computation

    // compute the current one
    regen_pos_dim_square_d(top_F, Jv_top_F, vars, F, Jv_F, codim_index, RPD);

    // save the actual size
    topSize = top_F->size - 1;

    // setup extra_F
    set_d(extra_F, &top_F->coord[topSize]);

    // setup Jv_extra_F
    change_size_vec_d(Jv_extra_F, vars->size);
    Jv_extra_F->size = vars->size;
    for (i = 0; i < vars->size; i++)
    {
      set_d(&Jv_extra_F->coord[i], &Jv_top_F->entry[topSize][i]);
    }

    // setup top_F, Jv_top_F & Jp_top_F
    change_size_mat_d(Jp_top_F, topSize, 1);
    top_F->size = Jv_top_F->rows = Jp_top_F->rows = topSize;
    Jv_top_F->cols = vars->size;
    Jp_top_F->cols = 1;
 
    // Jp = 0
    for (i = 0; i < topSize; i++)
    { 
      set_zero_d(&Jp_top_F->entry[i][0]);
    }
  }
  else
  { // different A's
    vec_d old_F, new_F;
    mat_d Jv_old_F, Jv_new_F;
    init_vec_d(old_F, 0); init_vec_d(new_F, 0);
    init_mat_d(Jv_old_F, 0, 0); init_mat_d(Jv_new_F, 0, 0);

    // compute the current one
    regen_pos_dim_square_d(new_F, Jv_new_F, vars, F, Jv_F, codim_index, RPD);

    // see if the first codim
    if (codim_index == 0)
    { // set top_F, Jv_top_F & Jp_top_F to 'empty'
      change_size_vec_d(top_F, 0);
      change_size_mat_d(Jv_top_F, 0, vars->size);
      change_size_mat_d(Jp_top_F, 0, 1);
      top_F->size = Jv_top_F->rows = Jp_top_F->rows = 0;
      Jv_top_F->cols = vars->size;
      Jp_top_F->cols = 1;

      // setup extra_F
      set_d(extra_F, &new_F->coord[0]);
     
      // setup Jv_extra_F
      change_size_vec_d(Jv_extra_F, vars->size);
      Jv_extra_F->size = vars->size;
      for (i = 0; i < vars->size; i++)
      {
        set_d(&Jv_extra_F->coord[i], &Jv_new_F->entry[0][i]);
      }
    }
    else
    { // compute gamma * t & 1 - t
      comp_d gamma_t, one_minus_t;
      mul_d(gamma_t, RPD->gamma_d, pathVars);
      set_one_d(one_minus_t);
      sub_d(one_minus_t, one_minus_t, pathVars);

      // compute old one
      regen_pos_dim_square_d(old_F, Jv_old_F, vars, F, Jv_F, codim_index - 1, RPD);

      // setup top_F, Jv_top_F & Jp_top_F
      change_size_vec_d(top_F, old_F->size);
      change_size_mat_d(Jv_top_F, old_F->size, vars->size);
      change_size_mat_d(Jp_top_F, old_F->size, 1);
      top_F->size = Jv_top_F->rows = Jp_top_F->rows = old_F->size;
      Jv_top_F->cols = vars->size;
      Jp_top_F->cols = 1;

      // top_F = old_F * gamma * t + (1-t) * new_F
      // Jp = old_F * gamma - new_F
      // Jv = Jv_old_F * gamma * t + (1-t) * Jv_new_F
      for (i = 0; i < old_F->size; i++)
      { // top_F
        mul_d(&top_F->coord[i], &old_F->coord[i], gamma_t);
        sum_mul_d(&top_F->coord[i], &new_F->coord[i], one_minus_t);
        // Jp
        mul_d(&Jp_top_F->entry[i][0], &old_F->coord[i], RPD->gamma_d);
        sub_d(&Jp_top_F->entry[i][0], &Jp_top_F->entry[i][0], &new_F->coord[i]);
        // Jv
        for (j = 0; j < vars->size; j++)
        {
          mul_d(&Jv_top_F->entry[i][j], &Jv_old_F->entry[i][j], gamma_t);
          sum_mul_d(&Jv_top_F->entry[i][j], &Jv_new_F->entry[i][j], one_minus_t);
        }
      }

      // setup extra_F
      set_d(extra_F, &new_F->coord[old_F->size]);

      // setup Jv_extra_F
      change_size_vec_d(Jv_extra_F, vars->size);
      Jv_extra_F->size = vars->size;
      for (i = 0; i < vars->size; i++)
      {
        set_d(&Jv_extra_F->coord[i], &Jv_new_F->entry[old_F->size][i]);
      }
    }

    // clear
    clear_vec_d(old_F); clear_vec_d(new_F);
    clear_mat_d(Jv_old_F); clear_mat_d(Jv_new_F);
  }

  return;
}

void regen_pos_dim_patch_d(comp_d patch, vec_d Jv_patch, point_d vars, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate the patch                                     *
\***************************************************************/
{ 
  int i, num_vars = vars->size;

  // compute patch: p * vars - 1 OR p * vars - hom_var
  if (RPD->PPD.num_var_gp)
  { // evaluate p * vars - 1
    set_one_d(patch);
    neg_d(patch, patch);

    for (i = 0; i < num_vars; i++)
    { // update patch
      sum_mul_d(patch, &RPD->patchCoeff_d->coord[i], &vars->coord[i]);
      // setup Jv_patch
      set_d(&Jv_patch->coord[i], &RPD->patchCoeff_d->coord[i]);
    }
  }
  else
  { // compute hom_var
    comp_d hom_var;
    set_d(hom_var, RPD->homVarConst_d);
    for (i = 0; i < num_vars; i++)
    {
      sum_mul_d(hom_var, &vars->coord[i], &RPD->H_d->coord[i]);
    }

    // evaluate p * vars - hom_var
    neg_d(patch, hom_var);

    for (i = 0; i < num_vars; i++)
    { // update patch
      sum_mul_d(patch, &RPD->patchCoeff_d->coord[i], &vars->coord[i]);
      // setup Jv_patch
      sub_d(&Jv_patch->coord[i], &RPD->patchCoeff_d->coord[i], &RPD->H_d->coord[i]);
    }
  }

  return;
}

int regen_pos_dim_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for pos-dim regeneration  *
\***************************************************************/
{
  regen_pos_dim_t *RPD = (regen_pos_dim_t *)ED;
  int i, j, orig_vars_size, codim_index = RPD->curr_codim, num_vars = vars->size; 
  int codim = RPD->codim[codim_index].codim, useIntrinsicSlice = RPD->codim[codim_index].useIntrinsicSlice;
  int size_top_F = codim - 1, num_linears = RPD->new_degrees[codim_index];

  comp_d one_minus_t, gamma_t, linear_prod, extra_F;
  vec_d F, top_F, Jv_extra_F, linears, d_linear_prod, tempVec;
  mat_d Jv_F, Jv_top_F, Jp_top_F; 

  // initialize
  init_vec_d(F, 0); init_vec_d(top_F, size_top_F); init_vec_d(Jv_extra_F, RPD->new_variables);
  init_vec_d(linears, num_linears); init_vec_d(d_linear_prod, RPD->new_variables); init_vec_d(tempVec, 0);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_top_F, size_top_F, RPD->new_variables); init_mat_d(Jp_top_F, size_top_F, 1);

  top_F->size = size_top_F;
  Jv_extra_F->size = RPD->new_variables;
  linears->size = num_linears;
  d_linear_prod->size = RPD->new_variables;
  Jv_top_F->rows = Jp_top_F->rows = size_top_F; 
  Jv_top_F->cols = RPD->new_variables;
  Jp_top_F->cols = 1;

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_d(gamma_t, RPD->gamma_d, pathVars);

  if (useIntrinsicSlice)
  { // using intrinsic formulation
    vec_d orig_vars;
    init_vec_d(orig_vars, 0);

    // convert vars to original variables (p + B*vars)
    intrinsicToExtrinsic_d(orig_vars, vars, RPD->codim[codim_index].B_d, RPD->codim[codim_index].p_d);

    // store the size
    orig_vars_size = orig_vars->size;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_d prog_vars;
      init_vec_d(prog_vars, 0);
      // topID_mul_mat_vec_d(prog_vars, RPD->C_d, orig_vars);
      mul_mat_vec_d(prog_vars, RPD->C_d, orig_vars);

      // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to orig_vars
      // topID_mat_mul_d(Jv_F, Jv_F, RPD->C_d);
      mat_mul_d(Jv_F, Jv_F, RPD->C_d);

      clear_vec_d(prog_vars);
    }
    else
    { // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, RPD->Prog);
    }

    // compute the top_F & extra_F from F & Jv_F
    regen_pos_dim_square_old_new_d(top_F, Jv_top_F, Jp_top_F, extra_F, Jv_extra_F, orig_vars, pathVars, F, Jv_F, codim_index, RPD);

    // setup tempVec to all 1's
    increase_size_vec_d(tempVec, num_linears);
    tempVec->size = num_linears;
    for (i = 0; i < num_linears; i++)
      set_one_d(&tempVec->coord[i]);

    // compute each of the linears and their product, and product over j != i
    set_one_d(linear_prod);
    for (i = 0; i < num_linears; i++)
    { // compute linear[i] 
      set_zero_d(&linears->coord[i]); 
      for (j = 0; j < orig_vars_size; j++) 
      { 
        sum_mul_d(&linears->coord[i], RPD->coeff_d[codim_index][i][j], &orig_vars->coord[j]);
      } 
      mul_d(linear_prod, linear_prod, &linears->coord[i]);

      // update the product over j != i
      for (j = 0; j < num_linears; j++)
        if (i != j)
        {
          mul_d(&tempVec->coord[j], &tempVec->coord[j], &linears->coord[i]);
        }
    } 

    // find d(linear_prod) = sum { d(linear_j) * product k != j }
    for (i = 0; i < orig_vars_size; i++)
    { // find the deriv w.r.t. x_i
      set_zero_d(&d_linear_prod->coord[i]);
      // loop over the linears
      for (j = 0; j < num_linears; j++)
      { // += d(l_j) * prod(k!=j)
        sum_mul_d(&d_linear_prod->coord[i], RPD->coeff_d[codim_index][j][i], &tempVec->coord[j]);
      }
    }

    // setup funcVals = 'TOP' [top_F] 'BOTTOM' [linear_prod] * gamma * t + (1-t) * [extra_F]
    // setup Jv = 'TOP' [Jv_top_F] 'BOTTOM' [d_linear_prod] * gamma * t + (1-t) * [Jv_extra_F]
    // setup Jp = 'TOP' [Jp_top_F] 'BOTTOM' [linear_prod] * gamma - [new_F]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, orig_vars_size); 
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jp->rows = num_vars;
    Jv->cols = orig_vars_size;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_top_F)
      { // funcVals
        set_d(&funcVals->coord[i], &top_F->coord[i]);
        // Jp
        set_d(&Jp->entry[i][0], &Jp_top_F->entry[i][0]);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &Jv_top_F->entry[i][j]);
        }
      }
      else
      { // funcVals
        mul_d(&funcVals->coord[i], linear_prod, gamma_t);
        sum_mul_d(&funcVals->coord[i], one_minus_t, extra_F);
        // Jp
        mul_d(&Jp->entry[i][0], linear_prod, RPD->gamma_d);
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], extra_F);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          mul_d(&Jv->entry[i][j], &d_linear_prod->coord[j], gamma_t);
          sum_mul_d(&Jv->entry[i][j], one_minus_t, &Jv_extra_F->coord[j]);
        }
      }

    // convert Jv to the input vars
    mat_mul_d(Jv, Jv, RPD->codim[codim_index].B_d);

    // clear memory
    clear_vec_d(orig_vars);
  } 
  else
  { // using extrinsic coordinates
    int num_extra = num_vars - codim - 1;
    comp_d patchVal;
    vec_d extraLinears, Jv_patch;

    init_vec_d(extraLinears, num_extra); init_vec_d(Jv_patch, num_vars);
    extraLinears->size = num_extra;
    Jv_patch->size = num_vars;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_d prog_vars;
      init_vec_d(prog_vars, 0);
      // topID_mul_mat_vec_d(prog_vars, RPD->C_d, vars);
      mul_mat_vec_d(prog_vars, RPD->C_d, vars);

      // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to vars
      // topID_mat_mul_d(Jv_F, Jv_F, RPD->C_d);
      mat_mul_d(Jv_F, Jv_F, RPD->C_d);

      clear_vec_d(prog_vars);
    }
    else
    { // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, RPD->Prog);
    }

    // compute the top_F & extra_F from F & Jv_F
    regen_pos_dim_square_old_new_d(top_F, Jv_top_F, Jp_top_F, extra_F, Jv_extra_F, vars, pathVars, F, Jv_F, codim_index, RPD);

    // setup tempVec to all 1's
    increase_size_vec_d(tempVec, num_linears);
    tempVec->size = num_linears;
    for (i = 0; i < num_linears; i++)
      set_one_d(&tempVec->coord[i]);

    // compute each of the linears and their product, and product over j != i
    set_one_d(linear_prod);
    for (i = 0; i < num_linears; i++)
    { // compute linear[i]
      set_zero_d(&linears->coord[i]);
      for (j = 0; j < num_vars; j++)
      {
        sum_mul_d(&linears->coord[i], RPD->coeff_d[codim_index][i][j], &vars->coord[j]);
      }
      mul_d(linear_prod, linear_prod, &linears->coord[i]);

      // update the product over j != i
      for (j = 0; j < num_linears; j++)
        if (i != j)
        {
          mul_d(&tempVec->coord[j], &tempVec->coord[j], &linears->coord[i]);
        }
    }

    // find d(linear_prod) = sum { d(linear_j) * product k != j }
    for (i = 0; i < num_vars; i++)
    { // find the deriv w.r.t. x_i
      set_zero_d(&d_linear_prod->coord[i]);
      // loop over the linears
      for (j = 0; j < num_linears; j++)
      { // += d(l_j) * prod(k!=j)
        sum_mul_d(&d_linear_prod->coord[i], RPD->coeff_d[codim_index][j][i], &tempVec->coord[j]);
      }
    }

    // evaluate the extra linears
    for (i = 0; i < num_extra; i++)
    { // regular linears
      set_zero_d(&extraLinears->coord[i]);
      for (j = 0; j < num_vars; j++)
      { // += coeff * vars
        sum_mul_d(&extraLinears->coord[i], RPD->coeff_d[codim + i][0][j], &vars->coord[j]);
      }
    }

    // evaluate the patch
    regen_pos_dim_patch_d(patchVal, Jv_patch, vars, RPD);

    // setup funcVals = 'TOP' [top_F] 'MIDDLE' [linear_prod] * gamma * t + (1-t) * [extra_F] 'BOTTOM' [extraLinears, patchVal]
    // setup Jv = 'TOP' [Jv_top_F] 'MIDDLE' [d_linear_prod] * gamma * t + (1-t) * [Jv_extra_F] 'BOTTOM' [Jv_extraLinears, Jv_patch]
    // setup Jp = 'TOP' [Jp_top_F] 'MIDDLE' [linear_prod] * gamma - [extra_F] 'BOTTOM' [0, 0]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, num_vars);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_top_F)
      { // funcVals
        set_d(&funcVals->coord[i], &top_F->coord[i]);
        // Jp
        set_d(&Jp->entry[i][0], &Jp_top_F->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &Jv_top_F->entry[i][j]);
        }
      }
      else if (i < codim)
      { // funcVals
        mul_d(&funcVals->coord[i], linear_prod, gamma_t);
        sum_mul_d(&funcVals->coord[i], one_minus_t, extra_F);
        // Jp
        mul_d(&Jp->entry[i][0], linear_prod, RPD->gamma_d);
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], extra_F);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          mul_d(&Jv->entry[i][j], &d_linear_prod->coord[j], gamma_t);
          sum_mul_d(&Jv->entry[i][j], one_minus_t, &Jv_extra_F->coord[j]);
        }
      }
      else if (i < num_vars - 1)
      { // funcVals
        set_d(&funcVals->coord[i], &extraLinears->coord[i - codim]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], RPD->coeff_d[i][0][j]);
        }
      }
      else
      { // funcVals
        set_d(&funcVals->coord[i], patchVal);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], &Jv_patch->coord[j]);
        }
      }

    clear_vec_d(extraLinears); clear_vec_d(Jv_patch);
  }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  // clear memory
  clear_vec_d(F); clear_vec_d(top_F); clear_vec_d(Jv_extra_F);
  clear_vec_d(linears); clear_vec_d(d_linear_prod); clear_vec_d(tempVec);
  clear_mat_d(Jv_F); clear_mat_d(Jv_top_F); clear_mat_d(Jp_top_F);

  return 0;
}

int regen_pos_dim_moving_linear_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluator to move linears for regeneration             *
\***************************************************************/
{
  regen_pos_dim_t *RPD = (regen_pos_dim_t *)ED;
  int i, j, orig_vars_size, codim_index = RPD->curr_codim, next_codim_index = RPD->curr_codim + 1, num_vars = vars->size, curr_deg = 0, next_deg = RPD->moving_degree;
  int codim = RPD->codim[codim_index].codim, next_codim = RPD->codim[next_codim_index].codim, useIntrinsicSlice = RPD->codim[next_codim_index].useIntrinsicSlice;
  int size_F = codim;

  comp_d one_minus_t, gamma_t, curr_linear, next_linear;
  vec_d F, old_F, d_linear;
  mat_d Jv_F, Jv_old_F;

  init_vec_d(F, 0); init_vec_d(old_F, size_F); init_vec_d(d_linear, RPD->new_variables);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_old_F, size_F, num_vars);
  F->size = Jv_F->rows = Jv_F->cols = 0;
  old_F->size = Jv_old_F->rows = size_F;
  Jv_old_F->cols = num_vars;
  d_linear->size = RPD->new_variables;

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_d(gamma_t, RPD->gamma_d, pathVars);

  if (useIntrinsicSlice)
  { // use intrinsic formulation
    vec_d orig_vars; 
    init_vec_d(orig_vars, 0);

    // convert vars to original vars (p + B*vars)
    intrinsicToExtrinsic_d(orig_vars, vars, RPD->codim[next_codim_index].B_d, RPD->codim[next_codim_index].p_d);

    // store the size
    orig_vars_size = orig_vars->size;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_d prog_vars;
      init_vec_d(prog_vars, 0);
      // topID_mul_mat_vec_d(prog_vars, RPD->C_d, orig_vars);
      mul_mat_vec_d(prog_vars, RPD->C_d, orig_vars);

      // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to orig_vars
      // topID_mat_mul_d(Jv_F, Jv_F, RPD->C_d);
      mat_mul_d(Jv_F, Jv_F, RPD->C_d);

      clear_vec_d(prog_vars);
    }
    else
    { // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, RPD->Prog);
    }
    
    // compute the old F = [ [A[c].*W] * F]
    regen_pos_dim_square_d(old_F, Jv_old_F, orig_vars, F, Jv_F, codim_index, RPD);

    // compute curr_linear & next_linear as well as d_linear = d(curr_linear * gamma * t + (1-t) * next_linear)
    set_zero_d(curr_linear);
    set_zero_d(next_linear);
    for (i = 0; i < orig_vars_size; i++)
    { // update curr_linear
      sum_mul_d(curr_linear, RPD->coeff_d[next_codim_index][curr_deg][i], &orig_vars->coord[i]);
      // update next_linear
      sum_mul_d(next_linear, RPD->coeff_d[next_codim_index][next_deg][i], &orig_vars->coord[i]);

      // compute d_linear
      mul_d(&d_linear->coord[i], RPD->coeff_d[next_codim_index][curr_deg][i], gamma_t);
      sum_mul_d(&d_linear->coord[i], RPD->coeff_d[next_codim_index][next_deg][i], one_minus_t);
    }

    // setup funcVals = 'TOP' [old_F] 'BOTTOM' [curr_linear] * gamma * t + (1-t) * [next_linear]
    // setup Jv = 'TOP' [Jv_old_F] 'BOTTOM' [d_linear]
    // setup Jp = 'TOP' [0] 'BOTTOM' [curr_linear] * gamma - [next_linear]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, orig_vars_size);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jp->rows = num_vars;
    Jv->cols = orig_vars_size;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_F)
      { // funcVals
        set_d(&funcVals->coord[i], &old_F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &Jv_old_F->entry[i][j]);
        }
      }
      else
      { // funcVals
        mul_d(&funcVals->coord[i], curr_linear, gamma_t);
        sum_mul_d(&funcVals->coord[i], one_minus_t, next_linear);
        // Jp
        mul_d(&Jp->entry[i][0], curr_linear, RPD->gamma_d);
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], next_linear);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &d_linear->coord[j]);
        }
      }

    // convert Jv to the input vars
    mat_mul_d(Jv, Jv, RPD->codim[next_codim_index].B_d);

    // clear memory
    clear_vec_d(orig_vars);
  }
  else
  { // using extrinsic coordinates
    int num_extra = num_vars - next_codim - 1;
    comp_d patchVal;
    vec_d extraLinears, Jv_patch;

    init_vec_d(extraLinears, num_extra); 
    init_vec_d(Jv_patch, num_vars);
    extraLinears->size = num_extra;
    Jv_patch->size = num_vars;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_d prog_vars;
      init_vec_d(prog_vars, 0);
      // topID_mul_mat_vec_d(prog_vars, RPD->C_d, vars);
      mul_mat_vec_d(prog_vars, RPD->C_d, vars);

      // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to vars
      // topID_mat_mul_d(Jv_F, Jv_F, RPD->C_d);
      mat_mul_d(Jv_F, Jv_F, RPD->C_d);

      clear_vec_d(prog_vars);
    }
    else
    { // evaluate F
      evalProg_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, RPD->Prog);
    }

    // compute the old F = [ [A[c].*W] * F]
    regen_pos_dim_square_d(old_F, Jv_old_F, vars, F, Jv_F, codim_index, RPD);

    // compute curr_linear & next_linear as well as d_linear = d(curr_linear * gamma * t + (1-t) * next_linear)
    set_zero_d(curr_linear);
    set_zero_d(next_linear);
    for (i = 0; i < num_vars; i++)
    { // update curr_linear
      sum_mul_d(curr_linear, RPD->coeff_d[next_codim_index][curr_deg][i], &vars->coord[i]);
      // update next_linear
      sum_mul_d(next_linear, RPD->coeff_d[next_codim_index][next_deg][i], &vars->coord[i]);

      // compute d_linear
      mul_d(&d_linear->coord[i], RPD->coeff_d[next_codim_index][curr_deg][i], gamma_t);
      sum_mul_d(&d_linear->coord[i], RPD->coeff_d[next_codim_index][next_deg][i], one_minus_t);
    }

    // evaluate the extra linears
    for (i = 0; i < num_extra; i++)
    { // regular linears
      set_zero_d(&extraLinears->coord[i]);
      for (j = 0; j < num_vars; j++)
      { // += coeff * vars
        sum_mul_d(&extraLinears->coord[i], RPD->coeff_d[next_codim + i][0][j], &vars->coord[j]);
      }
    }

    // evaluate the patch
    regen_pos_dim_patch_d(patchVal, Jv_patch, vars, RPD);

    // setup funcVals = 'TOP' [old_F] 'MIDDLE' [curr_linear] * gamma * t + (1-t) * [next_linear] 'BOTTOM' [extraLinears, patchVal]
    // setup Jv = 'TOP' [Jv_old_F] 'MIDDLE' [d_linear] 'BOTTOM' 'BOTTOM' [Jv_extraLinears, Jv_patch]
    // setup Jp = 'TOP' [0] 'MIDDLE' [curr_linear] * gamma - [next_linear] 'BOTTOM' [0, 0]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, num_vars);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_F)
      { // funcVals
        set_d(&funcVals->coord[i], &old_F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &Jv_old_F->entry[i][j]);
        }
      }
      else if (i < next_codim)
      { // funcVals
        mul_d(&funcVals->coord[i], curr_linear, gamma_t);
        sum_mul_d(&funcVals->coord[i], one_minus_t, next_linear);
        // Jp
        mul_d(&Jp->entry[i][0], curr_linear, RPD->gamma_d);
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], next_linear);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_d(&Jv->entry[i][j], &d_linear->coord[j]);
        }
      }
      else if (i < num_vars - 1)
      { // funcVals
        set_d(&funcVals->coord[i], &extraLinears->coord[i - next_codim]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], RPD->coeff_d[i][0][j]);
        }
      }
      else
      { // funcVals
        set_d(&funcVals->coord[i], patchVal);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], &Jv_patch->coord[j]);
        }
      }

    clear_vec_d(extraLinears); clear_vec_d(Jv_patch);
  }
  
  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  // clear memory
  clear_vec_d(F); clear_vec_d(old_F); clear_vec_d(d_linear);
  clear_mat_d(Jv_F); clear_mat_d(Jv_old_F);

  return 0;
}

// MP EVALUATION FUNCTIONS //

void regen_pos_dim_square_mp(vec_mp F_sqr, mat_mp Jv_sqr, point_mp vars, point_mp F, mat_mp Jv_F, int codim_index, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Turn the system into the square system for the codim   *
\***************************************************************/
{ // F_sqr = [A*.W] F
  int i, j, k, exp, F_size, cols, codim = RPD->codim[codim_index].codim, num_vars = vars->size;
  comp_mp tempComp, hom_var;
  vec_mp F_new;
  mat_mp Jv_F_new, sqrMat, sqrMat_deriv;

  // initialize
  init_mp(tempComp); init_mp(hom_var);
  F_size = F->size;
  cols = Jv_F->cols;

  init_vec_mp(F_new, F_size);
  F_new->size = F_size;

  init_mat_mp(Jv_F_new, F_size, cols);
  Jv_F_new->rows = F_size;
  Jv_F_new->cols = cols;

  init_mat_mp(sqrMat, codim, F_size);
  init_mat_mp(sqrMat_deriv, codim, F_size);
  sqrMat->rows = sqrMat_deriv->rows = codim;
  sqrMat->cols = sqrMat_deriv->cols = F_size;

  // reorder F to F_new and Jv_F to Jv_F_new
  for (i = 0; i < F_size; i++)
  { // do the swap on F
    k = RPD->P[i];
    set_mp(&F_new->coord[i], &F->coord[k]);
    for (j = 0; j < cols; j++)
    { // do the swap of Jv_F
      set_mp(&Jv_F_new->entry[i][j], &Jv_F->entry[k][j]);
    }
  }

  // compute hom_var
  set_mp(hom_var, RPD->homVarConst_mp);
  for (i = 0; i < num_vars; i++)
  {
    sum_mul_mp(hom_var, &vars->coord[i], &RPD->H_mp->coord[i]);
  }

  // setup sqrMat  = [A * h^W] & sqrMat_deriv = [W * A * h^(W-1)]
  for (i = 0; i < codim; i++)
    for (j = 0; j < F_size; j++)
      if (i > j)
      { // set to 0
        set_zero_mp(&sqrMat->entry[i][j]);
        set_zero_mp(&sqrMat_deriv->entry[i][j]);
      }
      else if (i == j)
      { // set to A
        set_mp(&sqrMat->entry[i][j], &RPD->A_mp[codim_index]->entry[i][j]);
        // set deriv to 0
        set_zero_mp(&sqrMat_deriv->entry[i][j]);
      }
      else // i < j
      { // compute A * h^W & set deriv to W * A * h^(W-1) &
        exp = RPD->W[codim_index][i][j - i - 1];
        if (exp > 0)
        { // compute h^(W-1)
          exp_mp_int(tempComp, hom_var, exp - 1);
          // compute A * h^(W-1)
          mul_mp(tempComp, tempComp, &RPD->A_mp[codim_index]->entry[i][j]);
          // compute sqrMat = A * h^W
          mul_mp(&sqrMat->entry[i][j], tempComp, hom_var);
          // compute deriv = W * A * h^(W-1)
          mul_rdouble_mp(&sqrMat_deriv->entry[i][j], tempComp, exp);
        }
        else
        { // h^0 = 1
          set_mp(&sqrMat->entry[i][j], &RPD->A_mp[codim_index]->entry[i][j]);
          set_zero_mp(&sqrMat_deriv->entry[i][j]);
        }
      }

  // compute F_sqr = [A * h^W] * F_new & Jv_sqr = d([A * h^W] * F_new) = [A * h^W] * Jv_F_new + [W * A * h^(W-1)] * F_new * d(h)
  change_size_vec_mp(F_sqr, codim);
  change_size_mat_mp(Jv_sqr, codim, cols);
  F_sqr->size = Jv_sqr->rows = codim;
  Jv_sqr->cols = cols;
  for (i = 0; i < codim; i++)
  { // compute F_sqr[i]
    set_zero_mp(&F_sqr->coord[i]);
    for (j = i; j < F_size; j++)
    { // += sqrMat * F_new
      sum_mul_mp(&F_sqr->coord[i], &sqrMat->entry[i][j], &F_new->coord[j]);
    }

    for (k = 0; k < cols; k++)
    { // compute Jv_F_sqr[i][k]
      mul_mp(&Jv_sqr->entry[i][k], &sqrMat->entry[i][i], &Jv_F_new->entry[i][k]);

      for (j = i + 1; j < F_size; j++)
      { // += sqrMat * Jv_F_new
        sum_mul_mp(&Jv_sqr->entry[i][k], &sqrMat->entry[i][j], &Jv_F_new->entry[j][k]);

        // += d(h^W) * F_new
        mul_mp(tempComp, &sqrMat_deriv->entry[i][j], &F_new->coord[j]);
        sum_mul_mp(&Jv_sqr->entry[i][k], tempComp, &RPD->H_mp->coord[k]);
      }
    }
  }

  clear_mp(tempComp); clear_mp(hom_var);
  clear_vec_mp(F_new);
  clear_mat_mp(Jv_F_new); clear_mat_mp(sqrMat); clear_mat_mp(sqrMat_deriv);

  return;
}

void regen_pos_dim_square_old_new_mp(vec_mp top_F, mat_mp Jv_top_F, mat_mp Jp_top_F, comp_mp extra_F, vec_mp Jv_extra_F, vec_mp vars, comp_mp pathVars, vec_mp F, mat_mp Jv_F, int codim_index, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: computes top_F = [old_F]*gamma*t + one_minus_t*[new_F] *
*   and extra_F to be the next function being added             *
\***************************************************************/
{
  int i, j, topSize;

  if (RPD->sameA)
  { // we have the same A so that we can be efficient in the computation

    // compute the current one
    regen_pos_dim_square_mp(top_F, Jv_top_F, vars, F, Jv_F, codim_index, RPD);

    // save the actual size
    topSize = top_F->size - 1;

    // setup extra_F
    set_mp(extra_F, &top_F->coord[topSize]);

    // setup Jv_extra_F
    change_size_vec_mp(Jv_extra_F, vars->size);
    Jv_extra_F->size = vars->size;
    for (i = 0; i < vars->size; i++)
    {
      set_mp(&Jv_extra_F->coord[i], &Jv_top_F->entry[topSize][i]);
    }

    // setup top_F, Jv_top_F & Jp_top_F
    change_size_mat_mp(Jp_top_F, topSize, 1);
    top_F->size = Jv_top_F->rows = Jp_top_F->rows = topSize;
    Jv_top_F->cols = vars->size;
    Jp_top_F->cols = 1;

    // Jp = 0
    for (i = 0; i < topSize; i++)
    {
      set_zero_mp(&Jp_top_F->entry[i][0]);
    }
  }
  else
  { // different A's
    vec_mp old_F, new_F;
    mat_mp Jv_old_F, Jv_new_F;
    init_vec_mp(old_F, 0); init_vec_mp(new_F, 0);
    init_mat_mp(Jv_old_F, 0, 0); init_mat_mp(Jv_new_F, 0, 0);

    // compute the current one
    regen_pos_dim_square_mp(new_F, Jv_new_F, vars, F, Jv_F, codim_index, RPD);

    // see if the first codim
    if (codim_index == 0)
    { // set top_F, Jv_top_F & Jp_top_F to 'empty'
      change_size_vec_mp(top_F, 0);
      change_size_mat_mp(Jv_top_F, 0, vars->size);
      change_size_mat_mp(Jp_top_F, 0, 1);
      top_F->size = Jv_top_F->rows = Jp_top_F->rows = 0;
      Jv_top_F->cols = vars->size;
      Jp_top_F->cols = 1;

      // setup extra_F
      set_mp(extra_F, &new_F->coord[0]);

      // setup Jv_extra_F
      change_size_vec_mp(Jv_extra_F, vars->size);
      Jv_extra_F->size = vars->size;
      for (i = 0; i < vars->size; i++)
      {
        set_mp(&Jv_extra_F->coord[i], &Jv_new_F->entry[0][i]);
      }
    }
    else
    { // compute gamma * t & 1 - t
      comp_mp gamma_t, one_minus_t;
      init_mp(gamma_t); init_mp(one_minus_t);
      mul_mp(gamma_t, RPD->gamma_mp, pathVars);
      set_one_mp(one_minus_t);
      sub_mp(one_minus_t, one_minus_t, pathVars);

      // compute old one
      regen_pos_dim_square_mp(old_F, Jv_old_F, vars, F, Jv_F, codim_index - 1, RPD);

      // setup top_F, Jv_top_F & Jp_top_F
      change_size_vec_mp(top_F, old_F->size);
      change_size_mat_mp(Jv_top_F, old_F->size, vars->size);
      change_size_mat_mp(Jp_top_F, old_F->size, 1);
      top_F->size = Jv_top_F->rows = Jp_top_F->rows = old_F->size;
      Jv_top_F->cols = vars->size;
      Jp_top_F->cols = 1;

      // top_F = old_F * gamma * t + (1-t) * new_F
      // Jp = old_F * gamma - new_F
      // Jv = Jv_old_F * gamma * t + (1-t) * Jv_new_F
      for (i = 0; i < old_F->size; i++)
      { // top_F
        mul_mp(&top_F->coord[i], &old_F->coord[i], gamma_t);
        sum_mul_mp(&top_F->coord[i], &new_F->coord[i], one_minus_t);
        // Jp
        mul_mp(&Jp_top_F->entry[i][0], &old_F->coord[i], RPD->gamma_mp);
        sub_mp(&Jp_top_F->entry[i][0], &Jp_top_F->entry[i][0], &new_F->coord[i]);
        // Jv
        for (j = 0; j < vars->size; j++)
        {
          mul_mp(&Jv_top_F->entry[i][j], &Jv_old_F->entry[i][j], gamma_t);
          sum_mul_mp(&Jv_top_F->entry[i][j], &Jv_new_F->entry[i][j], one_minus_t);
        }
      }

      // setup extra_F
      set_mp(extra_F, &new_F->coord[old_F->size]);

      // setup Jv_extra_F
      change_size_vec_mp(Jv_extra_F, vars->size);
      Jv_extra_F->size = vars->size;
      for (i = 0; i < vars->size; i++)
      {
        set_mp(&Jv_extra_F->coord[i], &Jv_new_F->entry[old_F->size][i]);
      }

      clear_mp(gamma_t); clear_mp(one_minus_t);
    }

    // clear
    clear_vec_mp(old_F); clear_vec_mp(new_F);
    clear_mat_mp(Jv_old_F); clear_mat_mp(Jv_new_F);
  }

  return;
}

void regen_pos_dim_patch_mp(comp_mp patch, vec_mp Jv_patch, point_mp vars, regen_pos_dim_t *RPD)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluate the patch                                     *
\***************************************************************/
{
  int i, num_vars = vars->size;

  // compute patch: p * vars - 1 OR p * vars - hom_var
  if (RPD->PPD.num_var_gp)
  { // evaluate p * vars - 1
    set_one_mp(patch);
    neg_mp(patch, patch);

    for (i = 0; i < num_vars; i++)
    { // update patch
      sum_mul_mp(patch, &RPD->patchCoeff_mp->coord[i], &vars->coord[i]);
      // setup Jv_patch
      set_mp(&Jv_patch->coord[i], &RPD->patchCoeff_mp->coord[i]);
    }
  }
  else
  { // compute hom_var
    comp_mp hom_var;
    init_mp(hom_var);
    set_mp(hom_var, RPD->homVarConst_mp);
    for (i = 0; i < num_vars; i++)
    {
      sum_mul_mp(hom_var, &vars->coord[i], &RPD->H_mp->coord[i]);
    }

    // evaluate p * vars - hom_var
    neg_mp(patch, hom_var);

    for (i = 0; i < num_vars; i++)
    { // update patch
      sum_mul_mp(patch, &RPD->patchCoeff_mp->coord[i], &vars->coord[i]);
      // setup Jv_patch
      sub_mp(&Jv_patch->coord[i], &RPD->patchCoeff_mp->coord[i], &RPD->H_mp->coord[i]);
    }
  }

  return;
}

int regen_pos_dim_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for pos-dim regeneration  *
\***************************************************************/
{
  regen_pos_dim_t *RPD = (regen_pos_dim_t *)ED;
  int i, j, orig_vars_size, codim_index = RPD->curr_codim, num_vars = vars->size;
  int codim = RPD->codim[codim_index].codim, useIntrinsicSlice = RPD->codim[codim_index].useIntrinsicSlice;
  int size_top_F = codim - 1, num_linears = RPD->new_degrees[codim_index];

  comp_mp one_minus_t, gamma_t, linear_prod, extra_F;
  vec_mp F, top_F, Jv_extra_F, linears, d_linear_prod, tempVec;
  mat_mp Jv_F, Jv_top_F, Jp_top_F;

  // initialize
  init_mp(one_minus_t); init_mp(gamma_t); init_mp(linear_prod); init_mp(extra_F);
  init_vec_mp(F, 0); init_vec_mp(top_F, size_top_F); init_vec_mp(Jv_extra_F, RPD->new_variables);
  init_vec_mp(linears, num_linears); init_vec_mp(d_linear_prod, RPD->new_variables); init_vec_mp(tempVec, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_top_F, size_top_F, RPD->new_variables); init_mat_mp(Jp_top_F, size_top_F, 1);

  top_F->size = size_top_F;
  Jv_extra_F->size = RPD->new_variables;
  linears->size = num_linears;
  d_linear_prod->size = RPD->new_variables;
  Jv_top_F->rows = Jp_top_F->rows = size_top_F;
  Jv_top_F->cols = RPD->new_variables;
  Jp_top_F->cols = 1;

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_mp(gamma_t, RPD->gamma_mp, pathVars);

  if (useIntrinsicSlice)
  { // using intrinsic formulation
    vec_mp orig_vars;
    init_vec_mp(orig_vars, 0);

    // convert vars to original variables (p + B*vars)
    intrinsicToExtrinsic_mp(orig_vars, vars, RPD->codim[codim_index].B_mp, RPD->codim[codim_index].p_mp);

    // store the size
    orig_vars_size = orig_vars->size;
    
    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_mp prog_vars;
      init_vec_mp(prog_vars, 0);
      // topID_mul_mat_vec_mp(prog_vars, RPD->C_mp, orig_vars);
      mul_mat_vec_mp(prog_vars, RPD->C_mp, orig_vars);

      // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to orig_vars
      // topID_mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);
      mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);

      clear_vec_mp(prog_vars);
    }
    else
    { // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, RPD->Prog);
    }

    // compute the top_F & extra_F from F & Jv_F
    regen_pos_dim_square_old_new_mp(top_F, Jv_top_F, Jp_top_F, extra_F, Jv_extra_F, orig_vars, pathVars, F, Jv_F, codim_index, RPD);

    // setup tempVec to all 1's
    increase_size_vec_mp(tempVec, num_linears);
    tempVec->size = num_linears;
    for (i = 0; i < num_linears; i++)
      set_one_mp(&tempVec->coord[i]);

    // compute each of the linears and their product, and product over j != i
    set_one_mp(linear_prod);
    for (i = 0; i < num_linears; i++)
    { // compute linear[i]
      set_zero_mp(&linears->coord[i]);
      for (j = 0; j < orig_vars_size; j++)
      {
        sum_mul_mp(&linears->coord[i], RPD->coeff_mp[codim_index][i][j], &orig_vars->coord[j]);
      }
      mul_mp(linear_prod, linear_prod, &linears->coord[i]);

      // update the product over j != i
      for (j = 0; j < num_linears; j++)
        if (i != j)
        {
          mul_mp(&tempVec->coord[j], &tempVec->coord[j], &linears->coord[i]);
        }
    }

    // find d(linear_prod) = sum { d(linear_j) * product k != j }
    for (i = 0; i < orig_vars_size; i++)
    { // find the deriv w.r.t. x_i
      set_zero_mp(&d_linear_prod->coord[i]);
      // loop over the linears
      for (j = 0; j < num_linears; j++)
      { // += d(l_j) * prod(k!=j)
        sum_mul_mp(&d_linear_prod->coord[i], RPD->coeff_mp[codim_index][j][i], &tempVec->coord[j]);
      }
    }

    // setup funcVals = 'TOP' [top_F] 'BOTTOM' [linear_prod] * gamma * t + (1-t) * [extra_F]
    // setup Jv = 'TOP' [Jv_top_F] 'BOTTOM' [d_linear_prod] * gamma * t + (1-t) * [Jv_extra_F]
    // setup Jp = 'TOP' [Jp_top_F] 'BOTTOM' [linear_prod] * gamma - [new_F]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, orig_vars_size);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jp->rows = num_vars;
    Jv->cols = orig_vars_size;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_top_F)
      { // funcVals
        set_mp(&funcVals->coord[i], &top_F->coord[i]);
        // Jp
        set_mp(&Jp->entry[i][0], &Jp_top_F->entry[i][0]);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &Jv_top_F->entry[i][j]);
        }
      }
      else
      { // funcVals
        mul_mp(&funcVals->coord[i], linear_prod, gamma_t);
        sum_mul_mp(&funcVals->coord[i], one_minus_t, extra_F);
        // Jp
        mul_mp(&Jp->entry[i][0], linear_prod, RPD->gamma_mp);
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], extra_F);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          mul_mp(&Jv->entry[i][j], &d_linear_prod->coord[j], gamma_t);
          sum_mul_mp(&Jv->entry[i][j], one_minus_t, &Jv_extra_F->coord[j]);
        }
      }

    // convert Jv to the input vars
    mat_mul_mp(Jv, Jv, RPD->codim[codim_index].B_mp);

    // clear memory
    clear_vec_mp(orig_vars);
  }
  else
  { // using extrinsic coordinates
    int num_extra = num_vars - codim - 1;
    comp_mp patchVal;
    vec_mp extraLinears, Jv_patch;

    init_mp(patchVal);
    init_vec_mp(extraLinears, num_extra); init_vec_mp(Jv_patch, num_vars);
    extraLinears->size = num_extra;
    Jv_patch->size = num_vars;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_mp prog_vars;
      init_vec_mp(prog_vars, 0);
      // topID_mul_mat_vec_mp(prog_vars, RPD->C_mp, vars);
      mul_mat_vec_mp(prog_vars, RPD->C_mp, vars);

      // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to vars
      // topID_mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);
      mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);

      clear_vec_mp(prog_vars);
    }
    else
    { // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, RPD->Prog);
    }

    // compute the top_F & extra_F from F & Jv_F
    regen_pos_dim_square_old_new_mp(top_F, Jv_top_F, Jp_top_F, extra_F, Jv_extra_F, vars, pathVars, F, Jv_F, codim_index, RPD);

    // setup tempVec to all 1's
    increase_size_vec_mp(tempVec, num_linears);
    tempVec->size = num_linears;
    for (i = 0; i < num_linears; i++)
      set_one_mp(&tempVec->coord[i]);

    // compute each of the linears and their product, and product over j != i
    set_one_mp(linear_prod);
    for (i = 0; i < num_linears; i++)
    { // compute linear[i]
      set_zero_mp(&linears->coord[i]);
      for (j = 0; j < num_vars; j++)
      {
        sum_mul_mp(&linears->coord[i], RPD->coeff_mp[codim_index][i][j], &vars->coord[j]);
      }
      mul_mp(linear_prod, linear_prod, &linears->coord[i]);

      // update the product over j != i
      for (j = 0; j < num_linears; j++)
        if (i != j)
        {
          mul_mp(&tempVec->coord[j], &tempVec->coord[j], &linears->coord[i]);
        }
    }

    // find d(linear_prod) = sum { d(linear_j) * product k != j }
    for (i = 0; i < num_vars; i++)
    { // find the deriv w.r.t. x_i
      set_zero_mp(&d_linear_prod->coord[i]);
      // loop over the linears
      for (j = 0; j < num_linears; j++)
      { // += d(l_j) * prod(k!=j)
        sum_mul_mp(&d_linear_prod->coord[i], RPD->coeff_mp[codim_index][j][i], &tempVec->coord[j]);
      }
    }

    // evaluate the extra linears
    for (i = 0; i < num_extra; i++)
    { // regular linears
      set_zero_mp(&extraLinears->coord[i]);
      for (j = 0; j < num_vars; j++)
      { // += coeff * vars
        sum_mul_mp(&extraLinears->coord[i], RPD->coeff_mp[codim + i][0][j], &vars->coord[j]);
      }
    }

    // evaluate the patch
    regen_pos_dim_patch_mp(patchVal, Jv_patch, vars, RPD);

    // setup funcVals = 'TOP' [top_F] 'MIDDLE' [linear_prod] * gamma * t + (1-t) * [extra_F] 'BOTTOM' [extraLinears, patchVal]
    // setup Jv = 'TOP' [Jv_top_F] 'MIDDLE' [d_linear_prod] * gamma * t + (1-t) * [Jv_extra_F] 'BOTTOM' [Jv_extraLinears, Jv_patch]
    // setup Jp = 'TOP' [Jp_top_F] 'MIDDLE' [linear_prod] * gamma - [extra_F] 'BOTTOM' [0, 0]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, num_vars);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_top_F)
      { // funcVals
        set_mp(&funcVals->coord[i], &top_F->coord[i]);
        // Jp
        set_mp(&Jp->entry[i][0], &Jp_top_F->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &Jv_top_F->entry[i][j]);
        }
      }
      else if (i < codim)
      { // funcVals
        mul_mp(&funcVals->coord[i], linear_prod, gamma_t);
        sum_mul_mp(&funcVals->coord[i], one_minus_t, extra_F);
        // Jp
        mul_mp(&Jp->entry[i][0], linear_prod, RPD->gamma_mp);
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], extra_F);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          mul_mp(&Jv->entry[i][j], &d_linear_prod->coord[j], gamma_t);
          sum_mul_mp(&Jv->entry[i][j], one_minus_t, &Jv_extra_F->coord[j]);
        }
      }
      else if (i < num_vars - 1)
      { // funcVals
        set_mp(&funcVals->coord[i], &extraLinears->coord[i - codim]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], RPD->coeff_mp[i][0][j]);
        }
      }
      else
      { // funcVals
        set_mp(&funcVals->coord[i], patchVal);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], &Jv_patch->coord[j]);
        }
      }

    clear_mp(patchVal);
    clear_vec_mp(extraLinears); clear_vec_mp(Jv_patch);
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]); // ds/dt = 1

  // clear memory
  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(linear_prod);
  clear_vec_mp(F); clear_vec_mp(top_F); clear_vec_mp(Jv_extra_F);
  clear_vec_mp(linears); clear_vec_mp(d_linear_prod); clear_vec_mp(tempVec);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_top_F); clear_mat_mp(Jp_top_F);

  return 0;
}

int regen_pos_dim_moving_linear_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluator to move linears for regeneration             *
\***************************************************************/
{
  regen_pos_dim_t *RPD = (regen_pos_dim_t *)ED;
  int i, j, orig_vars_size, codim_index = RPD->curr_codim, next_codim_index = RPD->curr_codim + 1, num_vars = vars->size, curr_deg = 0, next_deg = RPD->moving_degree;
  int codim = RPD->codim[codim_index].codim, next_codim = RPD->codim[next_codim_index].codim, useIntrinsicSlice = RPD->codim[next_codim_index].useIntrinsicSlice;
  int size_F = codim;

  comp_mp one_minus_t, gamma_t, curr_linear, next_linear;
  vec_mp F, old_F, d_linear;
  mat_mp Jv_F, Jv_old_F;

  init_mp(one_minus_t); init_mp(gamma_t); init_mp(curr_linear); init_mp(next_linear);
  init_vec_mp(F, 0); init_vec_mp(old_F, size_F); init_vec_mp(d_linear, RPD->new_variables);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_old_F, size_F, num_vars);
  F->size = Jv_F->rows = Jv_F->cols = 0;
  old_F->size = Jv_old_F->rows = size_F;
  Jv_old_F->cols = num_vars;
  d_linear->size = RPD->new_variables;

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_mp(gamma_t, RPD->gamma_mp, pathVars);

  if (useIntrinsicSlice)
  { // use intrinsic formulation
    vec_mp orig_vars;
    init_vec_mp(orig_vars, 0);

    // convert vars to original vars (p + B*vars)
    intrinsicToExtrinsic_mp(orig_vars, vars, RPD->codim[next_codim_index].B_mp, RPD->codim[next_codim_index].p_mp);

    // store the size
    orig_vars_size = orig_vars->size;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_mp prog_vars;
      init_vec_mp(prog_vars, 0);
      // topID_mul_mat_vec_mp(prog_vars, RPD->C_mp, orig_vars);
      mul_mat_vec_mp(prog_vars, RPD->C_mp, orig_vars);

      // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to orig_vars
      // topID_mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);
      mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);

      clear_vec_mp(prog_vars);
    }
    else
    { // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, RPD->Prog);
    }

    // compute the old F = [ [A[c].*W] * F]
    regen_pos_dim_square_mp(old_F, Jv_old_F, orig_vars, F, Jv_F, codim_index, RPD);

    // compute curr_linear & next_linear as well as d_linear = d(curr_linear * gamma * t + (1-t) * next_linear)
    set_zero_mp(curr_linear);
    set_zero_mp(next_linear);
    for (i = 0; i < orig_vars_size; i++)
    { // update curr_linear
      sum_mul_mp(curr_linear, RPD->coeff_mp[next_codim_index][curr_deg][i], &orig_vars->coord[i]);
      // update next_linear
      sum_mul_mp(next_linear, RPD->coeff_mp[next_codim_index][next_deg][i], &orig_vars->coord[i]);

      // compute d_linear
      mul_mp(&d_linear->coord[i], RPD->coeff_mp[next_codim_index][curr_deg][i], gamma_t);
      sum_mul_mp(&d_linear->coord[i], RPD->coeff_mp[next_codim_index][next_deg][i], one_minus_t);
    }

    // setup funcVals = 'TOP' [old_F] 'BOTTOM' [curr_linear] * gamma * t + (1-t) * [next_linear]
    // setup Jv = 'TOP' [Jv_old_F] 'BOTTOM' [d_linear]
    // setup Jp = 'TOP' [0] 'BOTTOM' [curr_linear] * gamma - [next_linear]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, orig_vars_size);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jp->rows = num_vars;
    Jv->cols = orig_vars->size;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_F)
      { // funcVals
        set_mp(&funcVals->coord[i], &old_F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &Jv_old_F->entry[i][j]);
        }
      }
      else
      { // funcVals
        mul_mp(&funcVals->coord[i], curr_linear, gamma_t);
        sum_mul_mp(&funcVals->coord[i], one_minus_t, next_linear);
        // Jp
        mul_mp(&Jp->entry[i][0], curr_linear, RPD->gamma_mp);
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], next_linear);
        // Jv
        for (j = 0; j < orig_vars_size; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &d_linear->coord[j]);
        }
      }

    // convert Jv to the input vars
    mat_mul_mp(Jv, Jv, RPD->codim[next_codim_index].B_mp);

    // clear memory
    clear_vec_mp(orig_vars);
  }
  else
  { // using extrinsic coordinates
    int num_extra = num_vars - next_codim - 1;
    comp_mp patchVal;
    vec_mp extraLinears, Jv_patch;

    init_mp(patchVal);
    init_vec_mp(extraLinears, num_extra); init_vec_mp(Jv_patch, num_vars);
    extraLinears->size = num_extra;
    Jv_patch->size = num_vars;

    // evaluate the original functions
    if (RPD->new_variables != RPD->orig_variables)
    { // convert orig_vars to prog_vars
      vec_mp prog_vars;
      init_vec_mp(prog_vars, 0);
      // topID_mul_mat_vec_mp(prog_vars, RPD->C_mp, vars);
      mul_mat_vec_mp(prog_vars, RPD->C_mp, vars);

      // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, prog_vars, pathVars, RPD->Prog);

      // convert Jv_F to vars
      // topID_mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);
      mat_mul_mp(Jv_F, Jv_F, RPD->C_mp);

      clear_vec_mp(prog_vars);
    }
    else
    { // evaluate F
      evalProg_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, RPD->Prog);
    }

    // compute the old F = [ [A[c].*W] * F]
    regen_pos_dim_square_mp(old_F, Jv_old_F, vars, F, Jv_F, codim_index, RPD);

    // compute curr_linear & next_linear as well as d_linear = d(curr_linear * gamma * t + (1-t) * next_linear)
    set_zero_mp(curr_linear);
    set_zero_mp(next_linear);
    for (i = 0; i < num_vars; i++)
    { // update curr_linear
      sum_mul_mp(curr_linear, RPD->coeff_mp[next_codim_index][curr_deg][i], &vars->coord[i]);
      // update next_linear
      sum_mul_mp(next_linear, RPD->coeff_mp[next_codim_index][next_deg][i], &vars->coord[i]);

      // compute d_linear
      mul_mp(&d_linear->coord[i], RPD->coeff_mp[next_codim_index][curr_deg][i], gamma_t);
      sum_mul_mp(&d_linear->coord[i], RPD->coeff_mp[next_codim_index][next_deg][i], one_minus_t);
    }

    // evaluate the extra linears
    for (i = 0; i < num_extra; i++)
    { // regular linears
      set_zero_mp(&extraLinears->coord[i]);
      for (j = 0; j < num_vars; j++)
      { // += coeff * vars
        sum_mul_mp(&extraLinears->coord[i], RPD->coeff_mp[next_codim + i][0][j], &vars->coord[j]);
      }
    }

    // evaluate the patch
    regen_pos_dim_patch_mp(patchVal, Jv_patch, vars, RPD);

    // setup funcVals = 'TOP' [old_F] 'MIDDLE' [curr_linear] * gamma * t + (1-t) * [next_linear] 'BOTTOM' [extraLinears, patchVal]
    // setup Jv = 'TOP' [Jv_old_F] 'MIDDLE' [d_linear] 'BOTTOM' 'BOTTOM' [Jv_extraLinears, Jv_patch]
    // setup Jp = 'TOP' [0] 'MIDDLE' [curr_linear] * gamma - [next_linear] 'BOTTOM' [0, 0]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, num_vars);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < size_F)
      { // funcVals
        set_mp(&funcVals->coord[i], &old_F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &Jv_old_F->entry[i][j]);
        }
      }
      else if (i < next_codim)
      { // funcVals
        mul_mp(&funcVals->coord[i], curr_linear, gamma_t);
        sum_mul_mp(&funcVals->coord[i], one_minus_t, next_linear);
        // Jp
        mul_mp(&Jp->entry[i][0], curr_linear, RPD->gamma_mp);
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], next_linear);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // setup Jv[i][j]
          set_mp(&Jv->entry[i][j], &d_linear->coord[j]);
        }
      }
      else if (i < num_vars - 1)
      { // funcVals
        set_mp(&funcVals->coord[i], &extraLinears->coord[i - next_codim]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], RPD->coeff_mp[i][0][j]);
        }
      }
      else
      { // funcVals
        set_mp(&funcVals->coord[i], patchVal);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], &Jv_patch->coord[j]);
        }
      }

    clear_mp(patchVal);
    clear_vec_mp(extraLinears); clear_vec_mp(Jv_patch);
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]); // ds/dt = 1

  // clear memory
  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(curr_linear); clear_mp(next_linear);
  clear_vec_mp(F); clear_vec_mp(old_F); clear_vec_mp(d_linear);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_old_F);

  return 0;
}



