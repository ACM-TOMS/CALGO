// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"

int standard_witness_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for witness set           *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;
 
  mat_d Jv_F, Jv_linears;
  vec_d F, linears;
  comp_d tempComp, tempComp2, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2); set_zero_d(hom_var);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_linears, 0, 0);
  init_vec_d(F, 0); init_vec_d(linears, 0);

  // find the homogenizing coordinate
  set_d(hom_var, W->codim[codim_index].homVarConst_d);
  for (i = 0; i < W->new_variables; i++)
  {
    sum_mul_d(hom_var, &W->codim[codim_index].H_d->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (W->new_variables != W->orig_variables)
  { // convert vars to prog_vars
    point_d prog_vars;
    init_point_d(prog_vars, 0);

    mul_mat_vec_d(prog_vars, W->C_d, vars);
    // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, W->Prog);
    // convert Jv to vars
    mat_mul_d(Jv, Jv, W->C_d);

    clear_point_d(prog_vars);
  }
  else
  { // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_d(F, W->num_funcs);
  change_size_mat_d(Jv_F, W->num_funcs, W->new_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_d(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_d(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_d(tempComp2, tempComp2, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_d(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_d(tempComp2); // d/d(hom_var) = 0
        set_d(tempComp, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->new_variables; k++)
      {
        sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_d(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_d->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_d(linears, W->new_variables - codim);
  change_size_mat_d(Jv_linears, W->new_variables - codim, W->new_variables);
  linears->size = Jv_linears->rows = W->new_variables - codim;
  Jv_linears->cols = W->new_variables;
  for (i = 0; i < linears->size; i++)
    if (i + 1 < linears->size)
    { // the first W->new_variables - codim - 1 linears are of the form B * vars
      set_zero_d(&linears->coord[i]);
      for (j = 0; j < vars->size; j++)
      { // update linears[i]
        sum_mul_d(&linears->coord[i], &W->codim[codim_index].B_d->entry[i][j], &vars->coord[j]);
        // update Jv_linears[i][j]
        set_d(&Jv_linears->entry[i][j], &W->codim[codim_index].B_d->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var
      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_d(&linears->coord[i]);
      }
      else
      { // == hom_var
        set_d(&linears->coord[i], hom_var);
      }
      neg_d(&linears->coord[i], &linears->coord[i]);
      for (j = 0; j < vars->size; j++)
      { // update linears[i]
        sum_mul_d(&linears->coord[i], &W->codim[codim_index].p_d->coord[j], &vars->coord[j]);
        // update Jv_linears[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_d(&Jv_linears->entry[i][j], &W->codim[codim_index].p_d->coord[j]);
        }
        else
        { // need to worry about p and H
          sub_d(&Jv_linears->entry[i][j], &W->codim[codim_index].p_d->coord[j], &W->codim[codim_index].H_d->coord[j]);
        }
      }
    }
  // combine everything!
  // funcVals = [F, linears]
  // Jv = [Jv_F, Jv_linears]
  // Jp = [0, 0]
  change_size_point_d(funcVals, W->new_variables);
  change_size_mat_d(Jv, W->new_variables, W->new_variables);
  change_size_mat_d(Jp, W->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->new_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);
      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else
    { // funcVals
      set_d(&funcVals->coord[i], &linears->coord[i-codim]);
      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_linears->entry[i-codim][j]);
      }
    }
  
  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_d(Jv_F); clear_mat_d(Jv_linears);
  clear_vec_d(F); clear_vec_d(linears);

  return 0;
}

int membership_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE: slice moving when we can have intrinsic reduction in   *
*  the number of variables that are being used                  *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for membership test            *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;

  mat_d Jv_F, Jv_linears_start, Jv_linears_finish;
  vec_d F, linears_start, linears_finish;
  comp_d tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_linears_start, 0, 0); init_mat_d(Jv_linears_finish, 0, 0);
  init_vec_d(F, 0); init_vec_d(linears_start, 0); init_vec_d(linears_finish, 0);
    
  // setup one_minus_t & gamma_t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  mul_d(gamma_t, W->gamma_d, pathVars);

  // find the homogenizing coordinate
  set_d(hom_var, W->codim[codim_index].homVarConst_d);
  for (i = 0; i < W->new_variables; i++)
  {
    sum_mul_d(hom_var, &W->codim[codim_index].H_d->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (W->new_variables != W->orig_variables)
  { // convert vars to prog_vars
    point_d prog_vars;
    init_point_d(prog_vars, 0);

    mul_mat_vec_d(prog_vars, W->C_d, vars);
    // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, W->Prog);
    // convert Jv to vars
    mat_mul_d(Jv, Jv, W->C_d);

    clear_point_d(prog_vars);
  }
  else
  { // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_d(F, W->num_funcs);
  change_size_mat_d(Jv_F, W->num_funcs, W->new_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_d(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_d(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_d(tempComp2, tempComp2, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_d(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_d(tempComp2); // d/d(hom_var) = 0
        set_d(tempComp, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->new_variables; k++)
      {
        sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_d(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_d->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_d(linears_start, W->new_variables - codim);
  change_size_vec_d(linears_finish, W->new_variables - codim);
  change_size_mat_d(Jv_linears_start, W->new_variables - codim, W->new_variables);
  change_size_mat_d(Jv_linears_finish, W->new_variables - codim, W->new_variables);
  linears_start->size = linears_finish->size = Jv_linears_start->rows = Jv_linears_finish->rows = W->new_variables - codim;
  Jv_linears_start->cols = Jv_linears_finish->cols = W->new_variables;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first W->new_variables - codim - 1 linears_start and linears_finish are of the form B * vars and Mat * vars - Vec
      set_zero_d(&linears_start->coord[i]);
      neg_d(&linears_finish->coord[i], &W->targetSliceVec_d->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears_start[i]
        sum_mul_d(&linears_start->coord[i], &W->codim[codim_index].B_d->entry[i][j], &vars->coord[j]);
        // update linears_finish[i]
        sum_mul_d(&linears_finish->coord[i], &W->targetSliceMat_d->entry[i][j], &vars->coord[j])
        // update Jv_linears_start[i][j]
        set_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].B_d->entry[i][j]);
        // update Jv_linears_finsish[i][j]
        set_d(&Jv_linears_finish->entry[i][j], &W->targetSliceMat_d->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var for both

      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_d(&linears_start->coord[i]);
      }
      else
      { // == hom_var
        set_d(&linears_start->coord[i], hom_var);
      }
      neg_d(&linears_start->coord[i], &linears_start->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears_start[i]
        sum_mul_d(&linears_start->coord[i], &W->codim[codim_index].p_d->coord[j], &vars->coord[j]);
        // update Jv_linears_start[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_d->coord[j]);
        }
        else
        { // need to worry about p and H
         sub_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_d->coord[j], &W->codim[codim_index].H_d->coord[j]);
        }
      }
    }

  // combine everything!
  // funcVals = [F, linears_start * t * gamma + (1 - t) * linears_finish, 'patch']
  // Jv = [Jv_F, Jv_linears_start * t * gamma + (1 - t) * Jv_linears_finish, 'Jv_patch']
  // Jp = [0, linears_start * gamma - linears_finish, 0]
  change_size_point_d(funcVals, W->new_variables);
  change_size_mat_d(Jv, W->new_variables, W->new_variables);
  change_size_mat_d(Jp, W->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->new_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else if (i + 1 < funcVals->size)
    { // funcVals
      mul_d(&funcVals->coord[i], &linears_finish->coord[i-codim], one_minus_t);
      sum_mul_d(&funcVals->coord[i], gamma_t, &linears_start->coord[i-codim]);

      // Jp
      neg_d(&Jp->entry[i][0], &linears_finish->coord[i-codim]);
      sum_mul_d(&Jp->entry[i][0], W->gamma_d, &linears_start->coord[i-codim]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        mul_d(&Jv->entry[i][j], &Jv_linears_finish->entry[i-codim][j], one_minus_t);
        sum_mul_d(&Jv->entry[i][j], gamma_t, &Jv_linears_start->entry[i-codim][j]);
      }
    }
    else
    { // funcVals
      set_d(&funcVals->coord[i], &linears_start->coord[i-codim]);
 
      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_linears_start->entry[i-codim][j]);
      }
    }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_d(Jv_F); clear_mat_d(Jv_linears_start); clear_mat_d(Jv_linears_finish);
  clear_vec_d(F); clear_vec_d(linears_start); clear_vec_d(linears_finish);

  return 0;
}

int decomposition_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE: totally extrinsic slice moving                         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for pure-dimensional decomp    *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;

  mat_d Jv_F, Jv_linears_start, Jv_linears_finish;
  vec_d F, linears_start, linears_finish;
  comp_d tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_linears_start, 0, 0); init_mat_d(Jv_linears_finish, 0, 0);
  init_vec_d(F, 0); init_vec_d(linears_start, 0); init_vec_d(linears_finish, 0);

  // setup one_minus_t & gamma_t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  mul_d(gamma_t, W->gamma_d, pathVars);

  // find the homogenizing coordinate
  set_d(hom_var, W->codim[codim_index].homVarConst_d);
  for (i = 0; i < W->orig_variables; i++)
  {
    sum_mul_d(hom_var, &W->codim[codim_index].H_d->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);

  // setup F & Jv_F using P
  change_size_vec_d(F, W->num_funcs);
  change_size_mat_d(Jv_F, W->num_funcs, W->orig_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->orig_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_d(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->orig_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_d(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_d(tempComp2, tempComp2, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_d(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_d(tempComp2); // d/d(hom_var) = 0
        set_d(tempComp, &W->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->orig_variables; k++)
      {
        sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_d(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_d->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_d(linears_start, W->orig_variables - codim);
  change_size_vec_d(linears_finish, W->orig_variables - codim);
  change_size_mat_d(Jv_linears_start, W->orig_variables - codim, W->orig_variables);
  change_size_mat_d(Jv_linears_finish, W->orig_variables - codim, W->orig_variables);
  linears_start->size = linears_finish->size = Jv_linears_start->rows = Jv_linears_finish->rows = W->orig_variables - codim;
  Jv_linears_start->cols = Jv_linears_finish->cols = W->orig_variables;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first W->orig_variables - codim - 1 linears_start and linears_finish are of the form B * vars and Mat * vars - Vec
      set_zero_d(&linears_start->coord[i]);
      neg_d(&linears_finish->coord[i], &W->targetSliceVec_d->coord[i]);

      for (j = 0; j < W->orig_variables; j++)
      { // update linears_start[i]
        sum_mul_d(&linears_start->coord[i], &W->codim[codim_index].B_d->entry[i][j], &vars->coord[j]);
        // update linears_finish[i]
        sum_mul_d(&linears_finish->coord[i], &W->targetSliceMat_d->entry[i][j], &vars->coord[j])
        // update Jv_linears_start[i][j]
        set_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].B_d->entry[i][j]);
        // update Jv_linears_finsish[i][j]
        set_d(&Jv_linears_finish->entry[i][j], &W->targetSliceMat_d->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var for start

      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_d(&linears_start->coord[i]);
      }
      else
      { // == hom_var
        set_d(&linears_start->coord[i], hom_var);
      }
      neg_d(&linears_start->coord[i], &linears_start->coord[i]);

      for (j = 0; j < W->orig_variables; j++)
      { // update linears_start[i]
        sum_mul_d(&linears_start->coord[i], &W->codim[codim_index].p_d->coord[j], &vars->coord[j]);
        // update Jv_linears_start[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_d->coord[j]);
        }
        else
        { // need to worry about p and H
         sub_d(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_d->coord[j], &W->codim[codim_index].H_d->coord[j]);
        }
      }
    }

  // combine everything!
  // funcVals = [F, linears_start * t * gamma + (1 - t) * linears_finish, 'patch']
  // Jv = [Jv_F, Jv_linears_start * t * gamma + (1 - t) * Jv_linears_finish, 'Jv_patch']
  // Jp = [0, linears_start * gamma - linears_finish, 0]
  change_size_point_d(funcVals, W->orig_variables);
  change_size_mat_d(Jv, W->orig_variables, W->orig_variables);
  change_size_mat_d(Jp, W->orig_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->orig_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else if (i + 1 < funcVals->size)
    { 
      k = i - codim;

      // funcVals
      mul_d(&funcVals->coord[i], &linears_finish->coord[k], one_minus_t);
      sum_mul_d(&funcVals->coord[i], gamma_t, &linears_start->coord[k]);

      // Jp
      neg_d(&Jp->entry[i][0], &linears_finish->coord[k]);
      sum_mul_d(&Jp->entry[i][0], W->gamma_d, &linears_start->coord[k]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        mul_d(&Jv->entry[i][j], &Jv_linears_finish->entry[k][j], one_minus_t);
        sum_mul_d(&Jv->entry[i][j], gamma_t, &Jv_linears_start->entry[k][j]);
      }
    }
    else
    { 
      k = i - codim;

      // funcVals
      set_d(&funcVals->coord[i], &linears_start->coord[k]);

      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_linears_start->entry[k][j]);
      }
    }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_d(Jv_F); clear_mat_d(Jv_linears_start); clear_mat_d(Jv_linears_finish);
  clear_vec_d(F); clear_vec_d(linears_start); clear_vec_d(linears_finish);

  return 0;
}

int standard_witness_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for witness set           *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;
    
  mat_mp Jv_F, Jv_linears;
  vec_mp F, linears;
  comp_mp tempComp, tempComp2, hom_var;

  init_mp(tempComp); init_mp(tempComp2); init_mp(hom_var);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_linears, 0, 0);
  init_vec_mp(F, 0); init_vec_mp(linears, 0);

  // find the homogenizing coordinate
  set_mp(hom_var, W->codim[codim_index].homVarConst_mp);
  for (i = 0; i < W->new_variables; i++)
  {
    sum_mul_mp(hom_var, &W->codim[codim_index].H_mp->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (W->new_variables != W->orig_variables)
  { // convert vars to prog_vars
    point_mp prog_vars;
    init_point_mp(prog_vars, 0);

    mul_mat_vec_mp(prog_vars, W->C_mp, vars);
    // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, W->Prog);
    // convert Jv to vars
    mat_mul_mp(Jv, Jv, W->C_mp);

    clear_point_mp(prog_vars);
  }
  else
  { // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_mp(F, W->num_funcs);
  change_size_mat_mp(Jv_F, W->num_funcs, W->new_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_mp(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_mp_int(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_mp(tempComp2, tempComp2, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_mp(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_mp(tempComp2); // d/d(hom_var) = 0
        set_mp(tempComp, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->new_variables; k++)
      {
        sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_mp->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_mp(linears, W->new_variables - codim);
  change_size_mat_mp(Jv_linears, W->new_variables - codim, W->new_variables);
  linears->size = Jv_linears->rows = W->new_variables - codim;
  Jv_linears->cols = W->new_variables;
  for (i = 0; i < linears->size; i++)
    if (i + 1 < linears->size)
    { // the first W->new_variables - codim - 1 linears are of the form B * vars
      set_zero_mp(&linears->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears[i]
        sum_mul_mp(&linears->coord[i], &W->codim[codim_index].B_mp->entry[i][j], &vars->coord[j]);
        // update Jv_linears[i][j]
        set_mp(&Jv_linears->entry[i][j], &W->codim[codim_index].B_mp->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var

      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_mp(&linears->coord[i]);
      }
      else
      { // == hom_var
        set_mp(&linears->coord[i], hom_var);
      }
      neg_mp(&linears->coord[i], &linears->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears[i]
        sum_mul_mp(&linears->coord[i], &W->codim[codim_index].p_mp->coord[j], &vars->coord[j]);
        // update Jv_linears[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_mp(&Jv_linears->entry[i][j], &W->codim[codim_index].p_mp->coord[j]);
        }
        else
        { // need to worry about p and H
          sub_mp(&Jv_linears->entry[i][j], &W->codim[codim_index].p_mp->coord[j], &W->codim[codim_index].H_mp->coord[j]);
        }
      }
    }

  // combine everything!
  // funcVals = [F, linears]
  // Jv = [Jv_F, Jv_linears]
  // Jp = [0, 0]
  change_size_point_mp(funcVals, W->new_variables);
  change_size_mat_mp(Jv, W->new_variables, W->new_variables);
  change_size_mat_mp(Jp, W->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->new_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else
    { // funcVals
      set_mp(&funcVals->coord[i], &linears->coord[i-codim]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_linears->entry[i-codim][j]);
      }
    }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_point_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);        // ds/dt = 1

  clear_mp(tempComp); clear_mp(tempComp2); clear_mp(hom_var);
  clear_vec_mp(F); clear_vec_mp(linears);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_linears);

  return 0;
}

int membership_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for membership test            *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;

  mat_mp Jv_F, Jv_linears_start, Jv_linears_finish;
  vec_mp F, linears_start, linears_finish;
  comp_mp tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  init_mp(tempComp); init_mp(one_minus_t); init_mp(gamma_t); init_mp(tempComp2); init_mp(hom_var);
  init_vec_mp(F, 0); init_vec_mp(linears_start, 0); init_vec_mp(linears_finish, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_linears_start, 0, 0); init_mat_mp(Jv_linears_finish, 0, 0);

  // setup one_minus_t & gamma_t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  mul_mp(gamma_t, W->gamma_mp, pathVars);

  // find the homogenizing coordinate
  set_mp(hom_var, W->codim[codim_index].homVarConst_mp);
  for (i = 0; i < W->new_variables; i++)
  {
    sum_mul_mp(hom_var, &W->codim[codim_index].H_mp->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (W->new_variables != W->orig_variables)
  { // convert vars to prog_vars
    point_mp prog_vars;
    init_point_mp(prog_vars, 0);

    mul_mat_vec_mp(prog_vars, W->C_mp, vars);
    // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, W->Prog);
    // convert Jv to vars
    mat_mul_mp(Jv, Jv, W->C_mp);

    clear_point_mp(prog_vars);
  }
  else
  { // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_mp(F, W->num_funcs);
  change_size_mat_mp(Jv_F, W->num_funcs, W->new_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_mp(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->new_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_mp_int(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_mp(tempComp2, tempComp2, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_mp(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_mp(tempComp2); // d/d(hom_var) = 0
        set_mp(tempComp, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->new_variables; k++)
      {
        sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_mp->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_mp(linears_start, W->new_variables - codim);
  change_size_vec_mp(linears_finish, W->new_variables - codim);
  change_size_mat_mp(Jv_linears_start, W->new_variables - codim, W->new_variables);
  change_size_mat_mp(Jv_linears_finish, W->new_variables - codim, W->new_variables);
  linears_start->size = linears_finish->size = Jv_linears_start->rows = Jv_linears_finish->rows = W->new_variables - codim;
  Jv_linears_start->cols = Jv_linears_finish->cols = W->new_variables;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first W->new_variables - codim - 1 linears_start and linears_finish are of the form B * vars and Mat * vars - Vec
      set_zero_mp(&linears_start->coord[i]);
      neg_mp(&linears_finish->coord[i], &W->targetSliceVec_mp->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears_start[i]
        sum_mul_mp(&linears_start->coord[i], &W->codim[codim_index].B_mp->entry[i][j], &vars->coord[j]);
        // update linears_finish[i]
        sum_mul_mp(&linears_finish->coord[i], &W->targetSliceMat_mp->entry[i][j], &vars->coord[j])
        // update Jv_linears_start[i][j]
        set_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].B_mp->entry[i][j]);
        // update Jv_linears_finsish[i][j]
        set_mp(&Jv_linears_finish->entry[i][j], &W->targetSliceMat_mp->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var for both

      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_mp(&linears_start->coord[i]);
      }
      else
      { // == hom_var
        set_mp(&linears_start->coord[i], hom_var);
      }
      neg_mp(&linears_start->coord[i], &linears_start->coord[i]);

      for (j = 0; j < vars->size; j++)
      { // update linears_start[i]
        sum_mul_mp(&linears_start->coord[i], &W->codim[codim_index].p_mp->coord[j], &vars->coord[j]);
        // update Jv_linears_start[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_mp->coord[j]);
        }
        else
        { // need to worry about p and H
          sub_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_mp->coord[j], &W->codim[codim_index].H_mp->coord[j]);
        }
      }
    }

  // combine everything!
  // funcVals = [F, linears_start * t * gamma + (1 - t) * linears_finish, 'patch']
  // Jv = [Jv_F, Jv_linears_start * t * gamma + (1 - t) * Jv_linears_finish, 'Jv_patch']
  // Jp = [0, linears_start * gamma - linears_finish, 0]
  change_size_vec_mp(funcVals, W->new_variables);
  change_size_mat_mp(Jv, W->new_variables, W->new_variables);
  change_size_mat_mp(Jp, W->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->new_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else if (i + 1 < funcVals->size)
    { // funcVals
      mul_mp(&funcVals->coord[i], &linears_finish->coord[i-codim], one_minus_t);
      sum_mul_mp(&funcVals->coord[i], gamma_t, &linears_start->coord[i-codim]);

      // Jp
      neg_mp(&Jp->entry[i][0], &linears_finish->coord[i-codim]);
      sum_mul_mp(&Jp->entry[i][0], W->gamma_mp, &linears_start->coord[i-codim]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        mul_mp(&Jv->entry[i][j], &Jv_linears_finish->entry[i-codim][j], one_minus_t);
        sum_mul_mp(&Jv->entry[i][j], gamma_t, &Jv_linears_start->entry[i-codim][j]);
      }
    }
    else
    { // funcVals
      set_mp(&funcVals->coord[i], &linears_start->coord[i-codim]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_linears_start->entry[i-codim][j]);
      }
    }

  // set parVals & parDer correctly
  change_size_vec_mp(parVals, 1);
  change_size_point_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);         // ds/dt = 1

  clear_mp(tempComp); clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(tempComp2); clear_mp(hom_var);
  clear_vec_mp(F); clear_vec_mp(linears_start); clear_vec_mp(linears_finish);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_linears_start); clear_mat_mp(Jv_linears_finish);

  return 0;
}

int decomposition_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE: totally extrinsic slice moving                         *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for pure-dimensional decomp    *
\***************************************************************/
{
  witness_t *W = (witness_t *)ED;
  int i, j, k, codim_index = W->curr_codim_index;
  int codim = W->codim[codim_index].codim;

  mat_mp Jv_F, Jv_linears_start, Jv_linears_finish;
  vec_mp F, linears_start, linears_finish;
  comp_mp tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  init_mp(tempComp); init_mp(tempComp2); init_mp(hom_var); init_mp(one_minus_t); init_mp(gamma_t);
  init_vec_mp(F, 0); init_vec_mp(linears_start, 0); init_vec_mp(linears_finish, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_linears_start, 0, 0); init_mat_mp(Jv_linears_finish, 0, 0);

  // setup one_minus_t & gamma_t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  mul_mp(gamma_t, W->gamma_mp, pathVars);

  // find the homogenizing coordinate
  set_mp(hom_var, W->codim[codim_index].homVarConst_mp);
  for (i = 0; i < W->orig_variables; i++)
  {
    sum_mul_mp(hom_var, &W->codim[codim_index].H_mp->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, W->Prog);

  // setup F & Jv_F using P
  change_size_vec_mp(F, W->num_funcs);
  change_size_mat_mp(Jv_F, W->num_funcs, W->orig_variables);
  F->size = Jv_F->rows = W->num_funcs;
  Jv_F->cols = W->orig_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_mp(&F->coord[i], &funcVals->coord[W->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[W->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = W->orig_variables;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < W->num_funcs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (W->codim[codim_index].W[i][j-codim] > 0)
      {
        exp_mp_int(tempComp2, hom_var, W->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_mp(tempComp2, tempComp2, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_mp(tempComp2, tempComp2, W->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_mp(tempComp2); // d/d(hom_var) = 0
        set_mp(tempComp, &W->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < W->orig_variables; k++)
      {
        sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &W->codim[codim_index].H_mp->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_mp(linears_start, W->orig_variables - codim);
  change_size_vec_mp(linears_finish, W->orig_variables - codim);
  change_size_mat_mp(Jv_linears_start, W->orig_variables - codim, W->orig_variables);
  change_size_mat_mp(Jv_linears_finish, W->orig_variables - codim, W->orig_variables);
  linears_start->size = linears_finish->size = Jv_linears_start->rows = Jv_linears_finish->rows = W->orig_variables - codim;
  Jv_linears_start->cols = Jv_linears_finish->cols = W->orig_variables;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first W->orig_variables - codim - 1 linears_start and linears_finish are of the form B * vars and Mat * vars - Vec
      set_zero_mp(&linears_start->coord[i]);
      neg_mp(&linears_finish->coord[i], &W->targetSliceVec_mp->coord[i]);

      for (j = 0; j < W->orig_variables; j++)
      { // update linears_start[i]
        sum_mul_mp(&linears_start->coord[i], &W->codim[codim_index].B_mp->entry[i][j], &vars->coord[j]);
        // update linears_finish[i]
        sum_mul_mp(&linears_finish->coord[i], &W->targetSliceMat_mp->entry[i][j], &vars->coord[j])
        // update Jv_linears_start[i][j]
        set_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].B_mp->entry[i][j]);
        // update Jv_linears_finsish[i][j]
        set_mp(&Jv_linears_finish->entry[i][j], &W->targetSliceMat_mp->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var for both

      // setup constant term
      if (W->PPD.num_var_gp)
      { // == 1
        set_one_mp(&linears_start->coord[i]);
      }
      else
      { // == hom_var
        set_mp(&linears_start->coord[i], hom_var);
      }
      neg_mp(&linears_start->coord[i], &linears_start->coord[i]);

      for (j = 0; j < W->orig_variables; j++)
      { // update linears_start[i]
        sum_mul_mp(&linears_start->coord[i], &W->codim[codim_index].p_mp->coord[j], &vars->coord[j]);
        // update Jv_linears_start[i][j]
        if (W->PPD.num_var_gp)
        { // only worry about p
          set_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_mp->coord[j]);
        }
        else
        { // need to worry about p and H
         sub_mp(&Jv_linears_start->entry[i][j], &W->codim[codim_index].p_mp->coord[j], &W->codim[codim_index].H_mp->coord[j]);
        }
      }
    }

  // combine everything!
  // funcVals = [F, linears_start * t * gamma + (1 - t) * linears_finish, 'patch']
  // Jv = [Jv_F, Jv_linears_start * t * gamma + (1 - t) * Jv_linears_finish, 'Jv_patch']
  // Jp = [0, 0, linears_start * gamma - linears_finish, 0]
  change_size_vec_mp(funcVals, W->orig_variables);
  change_size_mat_mp(Jv, W->orig_variables, W->orig_variables);
  change_size_mat_mp(Jp, W->orig_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = W->orig_variables;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < codim)
    { // funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      } 	
    }
    else if (i + 1 < funcVals->size)
    {
      k = i - codim;

      // funcVals
      mul_mp(&funcVals->coord[i], &linears_finish->coord[k], one_minus_t);
      sum_mul_mp(&funcVals->coord[i], gamma_t, &linears_start->coord[k]);

      // Jp
      neg_mp(&Jp->entry[i][0], &linears_finish->coord[k]);
      sum_mul_mp(&Jp->entry[i][0], W->gamma_mp, &linears_start->coord[k]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        mul_mp(&Jv->entry[i][j], &Jv_linears_finish->entry[k][j], one_minus_t);
        sum_mul_mp(&Jv->entry[i][j], gamma_t, &Jv_linears_start->entry[k][j]);
      }
    }
    else
    {
      k = i - codim;

      // funcVals
      set_mp(&funcVals->coord[i], &linears_start->coord[k]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_linears_start->entry[k][j]);
      }
    }

  // set parVals & parDer correctly
  change_size_vec_mp(parVals, 1);
  change_size_point_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);        // ds/dt = 1

  clear_mp(tempComp); clear_mp(tempComp2); clear_mp(hom_var); clear_mp(one_minus_t); clear_mp(gamma_t);
  clear_vec_mp(F); clear_vec_mp(linears_start); clear_vec_mp(linears_finish);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_linears_start); clear_mat_mp(Jv_linears_finish);

  return 0;
}

int slice_mover_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE: slice moving using a parallel slice                    *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for membership test            *
\***************************************************************/
{
  membership_slice_moving_t *SM = (membership_slice_moving_t *)ED;
  int i, j, k;

  mat_d Jv_F;
  vec_d F, Bx;
  comp_d one_minus_t, gamma_t, tempComp, gamma_minus_one, patch;

  init_mat_d(Jv_F, 0, 0);
  init_vec_d(F, 0); init_vec_d(Bx, 0);

  // setup one_minus_t, gamma_t, tempComp = 1 - t + gamma*t, and gamma_minus_one
  set_one_d(one_minus_t);
  sub_d(gamma_minus_one, SM->gamma_d, one_minus_t); // gamma - 1
  sub_d(one_minus_t, one_minus_t, pathVars); // 1 - t
  mul_d(gamma_t, SM->gamma_d, pathVars); // gamma * t
  add_d(tempComp, one_minus_t, gamma_t); // 1 - t + gamma * t

  // evaluate the function
  evalProg_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, SM->Prog);

  // see if we need to randomized to the correct size
  if (SM->A_d->rows >= 0 && SM->A_d->cols > 0)
  { // find F = [I A]F and Jv_F = [I A]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = SM->A_d->rows;
    Jv_F->cols = vars->size;
    for (i = 0; i < F->size; i++)
    { // add on the randomized part (i.e. [0 A]F & [0 A]Jv)
      for (j = 0; j < SM->A_d->cols; j++)
      { // F[i] += A[i][j] * F[j]
        sum_mul_d(&F->coord[i], &SM->A_d->entry[i][j], &F->coord[j + F->size]);

        // Jv_F[i][k] += A[i][j] * Jv_F[j][k]
        for (k = 0; k < vars->size; k++)
        {
          sum_mul_d(&Jv_F->entry[i][k], &SM->A_d->entry[i][j], &Jv_F->entry[j + F->size][k]);
        }
      }
    }
  }

  // evaluate the linears
  change_size_vec_d(Bx, SM->B_d->rows);
  Bx->size = SM->B_d->rows;
  for (i = 0; i < Bx->size; i++)
  {
    set_zero_d(&Bx->coord[i]);
    for (j = 0; j < SM->B_d->cols; j++)
    {
      sum_mul_d(&Bx->coord[i], &SM->B_d->entry[i][j], &vars->coord[j]);
    }
  }

  // evaluate the patch
  set_one_d(patch);
  neg_d(patch, patch);
  for (j = 0; j < SM->p_d->size; j++)
  {
    sum_mul_d(patch, &SM->p_d->coord[j], &vars->coord[j]);
  }

  // combine everything
  // funcVals = [F, (1-t+gamma_t)*Bx - (gamma_t*V_start + (1-t)*V_target), patch]
  // Jv = [Jv_F, (1-t+gamma_t)*B, p]
  // Jp = [0, (-1 + gamma)*Bx - gamma*V_start + V_target, 0]
  change_size_vec_d(funcVals, vars->size);
  change_size_mat_d(Jv, vars->size, vars->size);
  change_size_mat_d(Jp, vars->size, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = vars->size;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < F->size)
    { // funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else if (i + 1 < funcVals->size)
    { // funcVals
      mul_d(&funcVals->coord[i], gamma_t, &SM->startSliceVec_d->coord[i - F->size]); 
      sum_mul_d(&funcVals->coord[i], one_minus_t, &SM->targetSliceVec_d->coord[i - F->size]); 
      neg_d(&funcVals->coord[i], &funcVals->coord[i]);
      sum_mul_d(&funcVals->coord[i], &Bx->coord[i - F->size], tempComp);

      // Jp
      mul_d(&Jp->entry[i][0], SM->gamma_d, &SM->startSliceVec_d->coord[i - F->size]);
      neg_d(&Jp->entry[i][0], &Jp->entry[i][0]);
      add_d(&Jp->entry[i][0], &Jp->entry[i][0], &SM->targetSliceVec_d->coord[i - F->size]);
      sum_mul_d(&Jp->entry[i][0], gamma_minus_one, &Bx->coord[i - F->size]);

      // Jv
      for (j = 0; j < Jv->cols; j++)  
        if (j < SM->B_d->cols)
        {
          mul_d(&Jv->entry[i][j], tempComp, &SM->B_d->entry[i - F->size][j]);
        }
        else
        {
          set_zero_d(&Jv->entry[i][j]);
        }
    }
    else
    { // funcVals
      set_d(&funcVals->coord[i], patch);

      // Jp
      set_zero_d(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
        if (j < SM->p_d->size)
        {
          set_d(&Jv->entry[i][j], &SM->p_d->coord[j]);
        }
        else
        {
          set_zero_d(&Jv->entry[i][j]);
        }
    }
  
  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_d(Jv_F);
  clear_vec_d(F); clear_vec_d(Bx);

  return 0;
}

int slice_mover_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE: slice moving using a parallel slice                    *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: slice moving evaluation for membership test            *
\***************************************************************/
{
  membership_slice_moving_t *SM = (membership_slice_moving_t *)ED;
  int i, j, k;

  mat_mp Jv_F;
  vec_mp F, Bx;
  comp_mp one_minus_t, gamma_t, tempComp, gamma_minus_one, patch;

  init_mat_mp(Jv_F, 0, 0);
  init_vec_mp(F, 0); init_vec_mp(Bx, 0);
  init_mp(one_minus_t); init_mp(gamma_t); init_mp(tempComp); init_mp(gamma_minus_one); init_mp(patch);

  // setup one_minus_t, gamma_t, tempComp = 1 - t + gamma*t, and gamma_minus_one
  set_one_mp(one_minus_t);
  sub_mp(gamma_minus_one, SM->gamma_mp, one_minus_t); // gamma - 1
  sub_mp(one_minus_t, one_minus_t, pathVars); // 1 - t
  mul_mp(gamma_t, SM->gamma_mp, pathVars); // gamma * t
  add_mp(tempComp, one_minus_t, gamma_t); // 1 - t + gamma * t

  // evaluate the function
  evalProg_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, SM->Prog);

  // see if we need to randomized to the correct size
  if (SM->A_mp->rows >= 0 && SM->A_mp->cols > 0)
  { // find F = [I A]F and Jv_F = [I A]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = SM->A_mp->rows;
    Jv_F->cols = vars->size;
    for (i = 0; i < F->size; i++)
    { // add on the randomized part (i.e. [0 A]F & [0 A]Jv)
      for (j = 0; j < SM->A_mp->cols; j++)
      { // F[i] += A[i][j] * F[j]
        sum_mul_mp(&F->coord[i], &SM->A_mp->entry[i][j], &F->coord[j + F->size]);

        // Jv_F[i][k] += A[i][j] * Jv_F[j][k]
        for (k = 0; k < vars->size; k++)
        {
          sum_mul_mp(&Jv_F->entry[i][k], &SM->A_mp->entry[i][j], &Jv_F->entry[j + F->size][k]);
        }
      }
    }
  }

  // evaluate the linears
  change_size_vec_mp(Bx, SM->B_mp->rows);
  Bx->size = SM->B_mp->rows;
  for (i = 0; i < Bx->size; i++)
  {
    set_zero_mp(&Bx->coord[i]);
    for (j = 0; j < SM->B_mp->cols; j++)
    {
      sum_mul_mp(&Bx->coord[i], &SM->B_mp->entry[i][j], &vars->coord[j]);
    }
  }

  // evaluate the patch
  set_one_mp(patch);
  neg_mp(patch, patch);
  for (j = 0; j < SM->p_mp->size; j++)
  {
    sum_mul_mp(patch, &SM->p_mp->coord[j], &vars->coord[j]);
  }

  // combine everything
  // funcVals = [F, (1-t+gamma_t)*Bx - (gamma_t*V_start + (1-t)*V_target), patch]
  // Jv = [Jv_F, (1-t+gamma_t)*B, p]
  // Jp = [0, (-1 + gamma)*Bx - gamma*V_start + V_target, 0]
  change_size_vec_mp(funcVals, vars->size);
  change_size_mat_mp(Jv, vars->size, vars->size);
  change_size_mat_mp(Jp, vars->size, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = vars->size;
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < F->size)
    { // funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else if (i + 1 < funcVals->size)
    { // funcVals
      mul_mp(&funcVals->coord[i], gamma_t, &SM->startSliceVec_mp->coord[i - F->size]);
      sum_mul_mp(&funcVals->coord[i], one_minus_t, &SM->targetSliceVec_mp->coord[i - F->size]);
      neg_mp(&funcVals->coord[i], &funcVals->coord[i]);
      sum_mul_mp(&funcVals->coord[i], &Bx->coord[i - F->size], tempComp);

      // Jp
      mul_mp(&Jp->entry[i][0], SM->gamma_mp, &SM->startSliceVec_mp->coord[i - F->size]);
      neg_mp(&Jp->entry[i][0], &Jp->entry[i][0]);
      add_mp(&Jp->entry[i][0], &Jp->entry[i][0], &SM->targetSliceVec_mp->coord[i - F->size]);
      sum_mul_mp(&Jp->entry[i][0], gamma_minus_one, &Bx->coord[i - F->size]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
        if (j < SM->B_mp->cols)
        {
          mul_mp(&Jv->entry[i][j], tempComp, &SM->B_mp->entry[i - F->size][j]);
        }
        else
        {
          set_zero_mp(&Jv->entry[i][j]);
        }
    }
    else
    { // funcVals
      set_mp(&funcVals->coord[i], patch);

      // Jp
      set_zero_mp(&Jp->entry[i][0]);

      // Jv
      for (j = 0; j < Jv->cols; j++)
        if (j < SM->p_mp->size)
        {
          set_mp(&Jv->entry[i][j], &SM->p_mp->coord[j]);
        }
        else
        {
          set_zero_mp(&Jv->entry[i][j]);
        }
    }

  // set parVals & parDer correctly
  change_size_vec_mp(parVals, 1);
  change_size_point_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_mp(Jv_F);
  clear_vec_mp(F); clear_vec_mp(Bx);
  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(tempComp); clear_mp(gamma_minus_one); clear_mp(patch);

  return 0;
}

