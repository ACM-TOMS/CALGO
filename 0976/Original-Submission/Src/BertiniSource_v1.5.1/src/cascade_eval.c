// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"

int standard_codim_cascade_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for cascade               *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)ED;
  int i, j, k, codim_index = CD->curr_codim_index;
  int codim = CD->codim[codim_index].codim,
      rank = CD->system_rank;
  int curr_loc = CD->system_rank - codim; // number of active t variables in T, with codim >= 2

  mat_d Jv_F, RTB;
  vec_d F, RTBx;
  comp_d tempComp, tempComp2, Bx_curr_loc, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2);
  init_mat_d(Jv_F, 0, 0); init_mat_d(RTB, 0, 0);
  init_vec_d(F, 0); init_vec_d(RTBx, 0);

  // find the homogenizing coordinate
  set_d(hom_var, CD->homVarConst_d);
  for (i = 0; i < CD->new_variables; i++)
  {
    sum_mul_d(hom_var, &CD->H_d->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (CD->new_variables != CD->orig_variables)
  { // convert vars to prog_vars
    point_d prog_vars;
    init_point_d(prog_vars, 0);
  
    mul_mat_vec_d(prog_vars, CD->C_d, vars);
    // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, CD->Prog);
    // convert Jv to vars
    mat_mul_d(Jv, Jv, CD->C_d);

    clear_point_d(prog_vars);
  }
  else
  { // evaluate F
    evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, CD->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_d(F, CD->num_funcs);
  change_size_mat_d(Jv_F, CD->num_funcs, CD->new_variables);
  F->size = Jv_F->rows = CD->num_funcs;
  Jv_F->cols = CD->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_d(&F->coord[i], &funcVals->coord[CD->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[CD->P[i]][j]);
    }
  }

  // see if we need to randomize F & Jv_F to the correct size
  if (CD->num_funcs != CD->system_rank)
  { // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = CD->system_rank;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < rank; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = rank; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
        if (CD->W[i][j-rank] > 0)
        {
          exp_d(tempComp2, hom_var, CD->W[i][j-rank] - 1); // hom_var ^ (W[i][j] - 1)
          mul_d(tempComp2, tempComp2, &CD->A_d->entry[i][j-rank]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mul_rdouble_d(tempComp2, tempComp2, CD->W[i][j-rank]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
        }
        else // W[i][j] == 0
        {
          set_zero_d(tempComp2); // d/d(hom_var) = 0
          set_d(tempComp, &CD->A_d->entry[i][j-rank]); // A[i][j] * hom_var ^ 0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j]*hom_var^W[i][j]
        sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

        // Jv_F[i][k] += A[i][j]*hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_d(&Jv_F->entry[i][k], tempComp2, &CD->H_d->coord[k]);
        }
      }
    }
  }

  // setup T_d
  for (i = 0; i < CD->T_d->size; i++)
    if (i < curr_loc)
    {
      set_one_d(&CD->T_d->coord[i]);
    }
    else if (i == curr_loc)
    {
      set_d(&CD->T_d->coord[i], pathVars);
    }
    else
    {
      set_zero_d(&CD->T_d->coord[i]);
    }

  // evaluate R*T*B & R*T*B*x
  change_size_vec_d(RTBx, CD->system_rank);
  change_size_mat_d(RTB, CD->system_rank, CD->new_variables);
  RTB->rows = RTBx->size = CD->system_rank;
  RTB->cols = CD->new_variables;
  for (i = 0; i < rank; i++)
  {
    set_zero_d(&RTBx->coord[i]);
    for (j = 0; j < CD->new_variables; j++)
    { // find RTB[i][j]
      set_zero_d(&RTB->entry[i][j]);
      for (k = 0; k <= curr_loc && k < CD->T_d->size; k++)
      { // only need to evaluate the top since the bottom of T is 0!
        mul_d(tempComp, &CD->R_d->entry[i][k], &CD->T_d->coord[k]);
        sum_mul_d(&RTB->entry[i][j], tempComp, &CD->B_d->entry[k][j]);
      }
      // update RTBx[i]
      sum_mul_d(&RTBx->coord[i], &RTB->entry[i][j], &vars->coord[j]);
    }
  }

  // evaluate Bx_curr_loc
  set_zero_d(Bx_curr_loc);
  if (curr_loc < CD->B_d->rows)
  { // only evaluate if t is really used in the evaluation
    for (j = 0; j < CD->new_variables; j++)
    {
      sum_mul_d(Bx_curr_loc, &CD->B_d->entry[curr_loc][j], &vars->coord[j]); 
    }
  }

  // combine everything
  // funcVals = F + W'*RTBx and 'bottom' is patch
  // Jv = Jv_F + W'*RTB + d(W')/dx*RTBx and 'bottom' is patch coefficients
  // Jp = W'*R*[0,..,1,..,0]*B*x and 'bottom' is 0
  change_size_vec_d(funcVals, CD->new_variables);
  change_size_mat_d(Jv, CD->new_variables, CD->new_variables);
  change_size_mat_d(Jp, CD->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = CD->new_variables; // == system_rank + 1
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < rank)
    { // find funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);
      if (CD->W_prime[i] > 0)
      { // find (RTBx)_i * W_prime[i] * hom_var ^ (W_prime[i] - 1) and hom_var ^ W_prime[i]
        exp_d(tempComp2, hom_var, CD->W_prime[i] - 1); // hom_var ^ (W_prime[i] - 1)
        mul_d(tempComp, tempComp2, hom_var); // hom_var ^ W_prime[i]

        mul_d(tempComp2, tempComp2, &RTBx->coord[i]); // hom_var^(W_prime - 1) * (RTBx)_i
        mul_rdouble_d(tempComp2, tempComp2, CD->W_prime[i]); // W_prime[i] * hom_var ^ (W_prime[i] - 1) * (RTBx)_i
      }
      else
      { // W_prime == 0
        set_zero_d(tempComp2); // d/d(hom_var) hom_var^0 = 0
        set_one_d(tempComp); // hom_var ^ 0
      }
      sum_mul_d(&funcVals->coord[i], tempComp, &RTBx->coord[i]); // F + W'*RTBx

      // setup Jv
      for (j = 0; j < CD->new_variables; j++)
      { // setup Jv[i][j]
        mul_d(&Jv->entry[i][j], tempComp, &RTB->entry[i][j]); // hom_var ^ W_prime[i] * RTB[i][j]
        add_d(&Jv->entry[i][j], &Jv->entry[i][j], &Jv_F->entry[i][j]); // Jv_F + W'*RTB

        sum_mul_d(&Jv->entry[i][j], tempComp2, &CD->H_d->coord[j]); // Jv += d(W')/dx_j * RTBx
      }
  
      // setup Jp
      mul_d(&Jp->entry[i][0], tempComp, &CD->R_d->entry[i][curr_loc]); // W'*R[i][curr_loc]
      mul_d(&Jp->entry[i][0], &Jp->entry[i][0], Bx_curr_loc);
    }
    else
    { // evaluate the patch and setup the bottom of Jv & Jp 
      // patch = p * vars - 1 or p * vars - hom_var

      // setup constant term
      if (CD->PPD.num_var_gp)
      { // == 1
        set_one_d(&funcVals->coord[i]);
      }
      else
      { // == hom_var
        set_d(&funcVals->coord[i], hom_var);
      }
      neg_d(&funcVals->coord[i], &funcVals->coord[i]);

      for (j = 0; j < CD->new_variables; j++)
      { // update funcVals
        sum_mul_d(&funcVals->coord[i], &CD->p_d->coord[j], &vars->coord[j]);
        // setup Jv
        if (CD->PPD.num_var_gp)
        { // only worry about p
          set_d(&Jv->entry[i][j], &CD->p_d->coord[j]);
        }
        else
        { // need to worry about p and H
         sub_d(&Jv->entry[i][j], &CD->p_d->coord[j], &CD->H_d->coord[j]);
        }
      }
      // setup Jp
      set_zero_d(&Jp->entry[i][0]);
    }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1


  clear_mat_d(Jv_F); clear_mat_d(RTB);
  clear_vec_d(F); clear_vec_d(RTBx);

  return 0;
}

int initial_codim_cascade_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initial function evaluation for cascade                *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)ED;
  int i, j, rank = CD->system_rank;

  vec_d SS;
  mat_d Jv_SS;
  comp_d tempComp, gamma_t, one_minus_t;

  set_zero_d(tempComp);
  init_vec_d(SS, 0);
  init_mat_d(Jv_SS, 0, 0);

  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars); // 1 - t
  mul_d(gamma_t, CD->gamma_d, pathVars); // gamma * t

  // evaluate the cascade function 'F + RTBx' where T = identity
  standard_codim_cascade_eval_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ED);

  // evaluate the start system
  change_size_vec_d(SS, CD->system_rank);
  change_size_mat_d(Jv_SS, CD->system_rank, CD->new_variables);
  SS->size = Jv_SS->rows = CD->system_rank;
  Jv_SS->cols = CD->new_variables;
  for (i = 0; i < rank; i++)
  {
    if (CD->new_degrees[i] > 0)
    {
      exp_d(&SS->coord[i], &vars->coord[0], CD->new_degrees[i] - 1);
      neg_d(&SS->coord[i], &SS->coord[i]); // - x_0 ^ (d_i - 1)
      exp_d(tempComp, &vars->coord[i+1], CD->new_degrees[i] - 1); // x_(i+1) ^ (d_i - 1)
    }
    else
    {
      set_zero_d(&SS->coord[i]);
      set_zero_d(tempComp);
    }

    // zero out row of Jv_SS since we are only storing to a couple of locations in this row
    for (j = 0; j < Jv_SS->cols; j++)
    {
      set_zero_d(&Jv_SS->entry[i][j]);
    }

    // d_i * -x_0 ^ (d_i - 1)
    mul_rdouble_d(&Jv_SS->entry[i][0], &SS->coord[i], CD->new_degrees[i]);
    // d_i * x_(i+1)^(d_i - 1)
    mul_rdouble_d(&Jv_SS->entry[i][i+1], tempComp, CD->new_degrees[i]);

    // -x_0 ^ d_i
    mul_d(&SS->coord[i], &SS->coord[i], &vars->coord[0]);
    // x_(i+1) ^ d_i - x_0 ^ d_i
    sum_mul_d(&SS->coord[i], tempComp, &vars->coord[i+1]);
  }

  // combine everything - sizes are already setup!
  // funcVals = funcVals * one_minus_t + gamma_t * SS and 'bottom' is patch - already setup
  // Jv = Jv * one_minus_t + gamma_t * Jv_SS and 'bottom' is patch coefficients - already setup
  // Jp = -funcVals + gamma * SS and 'bottom' is 0 - already setup
  for (i = 0; i < rank; i++)
  { // setup Jp
    mul_d(&Jp->entry[i][0], CD->gamma_d, &SS->coord[i]);
    sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &funcVals->coord[i]);

    // find funcVals
    mul_d(&funcVals->coord[i], &funcVals->coord[i], one_minus_t);
    sum_mul_d(&funcVals->coord[i], gamma_t, &SS->coord[i]);

    // setup Jv
    for (j = 0; j < CD->new_variables; j++)
    { // setup Jv[i][j]
      mul_d(&Jv->entry[i][j], &Jv->entry[i][j], one_minus_t);
      sum_mul_d(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
    }
  }

  // parVals & parDer are already setup correctly

  clear_vec_d(SS);
  clear_mat_d(Jv_SS);

  return 0;
}

/// MP EVALUATION FUNCTIONS ///

int standard_codim_cascade_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for cascade               *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)ED;
  int i, j, k, codim_index = CD->curr_codim_index;
  int codim = CD->codim[codim_index].codim,
      rank = CD->system_rank;
  int curr_loc = CD->system_rank - codim; // number of active t variables in T, with codim >= 2

  mat_mp Jv_F, RTB;
  vec_mp F, RTBx;
  comp_mp tempComp, tempComp2, Bx_curr_loc, hom_var;

  init_mat_mp(Jv_F, 0, 0); init_mat_mp(RTB, 0, 0);
  init_vec_mp(F, 0); init_vec_mp(RTBx, 0);
  init_mp(tempComp); init_mp(Bx_curr_loc); init_mp(tempComp2); init_mp(hom_var);
 
  // find the homogenizing coordinate
  set_mp(hom_var, CD->homVarConst_mp);
  for (i = 0; i < CD->new_variables; i++)
  {
    sum_mul_mp(hom_var, &CD->H_mp->coord[i], &vars->coord[i]);
  }

  // evaluate the original function
  if (CD->new_variables != CD->orig_variables)
  { // convert vars to prog_vars
    point_mp prog_vars;
    init_point_mp(prog_vars, 0);

    mul_mat_vec_mp(prog_vars, CD->C_mp, vars);
    // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, CD->Prog);
    // convert Jv to vars
    mat_mul_mp(Jv, Jv, CD->C_mp);

    clear_point_mp(prog_vars);
  }
  else
  { // evaluate F
    evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, CD->Prog);
  }

  // setup F & Jv_F using P
  change_size_vec_mp(F, CD->num_funcs);
  change_size_mat_mp(Jv_F, CD->num_funcs, CD->new_variables);
  F->size = Jv_F->rows = CD->num_funcs;
  Jv_F->cols = CD->new_variables;
  for (i = 0; i < F->size; i++)
  { // put F(i) = f(P[i])
    set_mp(&F->coord[i], &funcVals->coord[CD->P[i]]);
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[CD->P[i]][j]);
    }
  }

  // see if we need to randomize F & Jv_F to the correct size
  if (CD->num_funcs != CD->system_rank)
  { // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = CD->system_rank;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < rank; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = rank; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and A[i][j]*hom_var^W[i][j]
        if (CD->W[i][j-rank] > 0)
        {
          exp_mp_int(tempComp2, hom_var, CD->W[i][j-rank] - 1); // hom_var ^ (W[i][j] - 1)
          mul_mp(tempComp2, tempComp2, &CD->A_mp->entry[i][j-rank]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mpf_mul_ui(tempComp2->r, tempComp2->r, CD->W[i][j-rank]); // W[i][j]*A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mpf_mul_ui(tempComp2->i, tempComp2->i, CD->W[i][j-rank]);
        }
        else // W[i][j] == 0
        {
          set_zero_mp(tempComp2); // == 0
          set_mp(tempComp, &CD->A_mp->entry[i][j-rank]); // A[i][j] * hom_var^0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j]*hom_var^W[i][j]
        sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);

        // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + W[i][j]*A[i][j]*hom_var^(W[i][j] - 1)*F[j] * d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &CD->H_mp->coord[k]);
        }
      }
    }
  }

  // setup T_mp
  for (i = 0; i < CD->T_mp->size; i++)
    if (i < curr_loc)
    {
      set_one_mp(&CD->T_mp->coord[i]);
    }
    else if (i == curr_loc)
    {
      set_mp(&CD->T_mp->coord[i], pathVars);
    }
    else
    {
      set_zero_mp(&CD->T_mp->coord[i]);
    }

  // evaluate R*T*B & R*T*B*x
  change_size_vec_mp(RTBx, CD->system_rank);
  change_size_mat_mp(RTB, CD->system_rank, CD->new_variables);
  RTB->rows = RTBx->size = CD->system_rank;
  RTB->cols = CD->new_variables;
  for (i = 0; i < rank; i++)
  {
    set_zero_mp(&RTBx->coord[i]);
    for (j = 0; j < CD->new_variables; j++)
    { // find RTB[i][j]
      set_zero_mp(&RTB->entry[i][j]);
      for (k = 0; k <= curr_loc && k < CD->T_mp->size; k++)
      { // only need to evaluate the top since the bottom of T is 0!
        mul_mp(tempComp, &CD->R_mp->entry[i][k], &CD->T_mp->coord[k]);
        sum_mul_mp(&RTB->entry[i][j], tempComp, &CD->B_mp->entry[k][j]);
      }
      // update RTBx[i]
      sum_mul_mp(&RTBx->coord[i], &RTB->entry[i][j], &vars->coord[j]);
    }
  }

  // evaluate Bx_curr_loc
  set_zero_mp(Bx_curr_loc);
  if (curr_loc < CD->B_mp->rows)
  { // only evaluate if t is really used in the evaluation
    for (j = 0; j < CD->new_variables; j++)
    {
      sum_mul_mp(Bx_curr_loc, &CD->B_mp->entry[curr_loc][j], &vars->coord[j]);
    }
  }

  // combine everything
  // funcVals = F + W'*RTBx and 'bottom' is patch
  // Jv = Jv_F + W'*RTB + d(W')/dx*RTBx and 'bottom' is patch coefficients
  // Jp = W'*R*[0,..,1,..,0]*B*x and 'bottom' is 0
  change_size_vec_mp(funcVals, CD->new_variables);
  change_size_mat_mp(Jv, CD->new_variables, CD->new_variables);
  change_size_mat_mp(Jp, CD->new_variables, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = CD->new_variables; // == system_rank + 1
  Jp->cols = 1;
  for (i = 0; i < funcVals->size; i++)
    if (i < rank)
    { // find funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);
      if (CD->W_prime[i] > 0)
      { // find (RTBx)_i * W_prime[i] * hom_var ^ (W_prime[i] - 1) and hom_var ^ W_prime[i]
        exp_mp_int(tempComp2, hom_var, CD->W_prime[i] - 1); // hom_var ^ (W_prime[i] - 1)
        mul_mp(tempComp, tempComp2, hom_var); // hom_var ^ W_prime[i]

        mul_mp(tempComp2, tempComp2, &RTBx->coord[i]); // hom_var^(W_prime[i] - 1) * (RTBx)_i

        mpf_mul_ui(tempComp2->r, tempComp2->r, CD->W_prime[i]); // (RTBx)_i * W_prime[i] * hom_var ^ (W_prime[i] - 1)
        mpf_mul_ui(tempComp2->i, tempComp2->i, CD->W_prime[i]);
      }
      else
      { // W_prime[i] == 0
        set_zero_mp(tempComp2); // == 0
        set_one_mp(tempComp); // hom_var ^ 0 == 1
      }
      sum_mul_mp(&funcVals->coord[i], tempComp, &RTBx->coord[i]); // F + W'*RTBx

      // setup Jv
      for (j = 0; j < CD->new_variables; j++)
      { // setup Jv[i][j]
        mul_mp(&Jv->entry[i][j], tempComp, &RTB->entry[i][j]); // x_0 ^ W_prime[i] * RTB[i][j]
        add_mp(&Jv->entry[i][j], &Jv->entry[i][j], &Jv_F->entry[i][j]); // Jv_F + W'*RTB

        sum_mul_mp(&Jv->entry[i][j], tempComp2, &CD->H_mp->coord[j]); // Jv += d(W')/dx_j * RTBx
      }

      // setup Jp
      mul_mp(&Jp->entry[i][0], tempComp, &CD->R_mp->entry[i][curr_loc]); // W'*R[i][curr_loc]
      mul_mp(&Jp->entry[i][0], &Jp->entry[i][0], Bx_curr_loc);
    }
    else
    { // evaluate the patch and setup the bottom of Jv & Jp
      // patch = p * vars - 1 or p * vars - hom_var

      // setup constant term
      if (CD->PPD.num_var_gp)
      { // == 1
        set_one_mp(&funcVals->coord[i]);
      }
      else
      { // == hom_var
        set_mp(&funcVals->coord[i], hom_var);
      }
      neg_mp(&funcVals->coord[i], &funcVals->coord[i]);

      for (j = 0; j < CD->new_variables; j++)
      { // update funcVals
        sum_mul_mp(&funcVals->coord[i], &CD->p_mp->coord[j], &vars->coord[j]);
        // setup Jv
        if (CD->PPD.num_var_gp)
        { // only worry about p
          set_mp(&Jv->entry[i][j], &CD->p_mp->coord[j]);
        }
        else
        { // need to worry about p and H
         sub_mp(&Jv->entry[i][j], &CD->p_mp->coord[j], &CD->H_mp->coord[j]);
        }
      }

      // setup Jp
      set_zero_mp(&Jp->entry[i][0]);
    }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);        // ds/dt = 1

  clear_mp(tempComp); clear_mp(Bx_curr_loc); clear_mp(tempComp2); clear_mp(hom_var);
  clear_vec_mp(F); clear_vec_mp(RTBx);
  clear_mat_mp(Jv_F); clear_mat_mp(RTB);

  return 0;
}

int initial_codim_cascade_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initial function evaluation for cascade                *
\***************************************************************/
{
  cascade_t *CD = (cascade_t *)ED;

  int i, j, rank = CD->system_rank;

  vec_mp SS;
  mat_mp Jv_SS;
  comp_mp tempComp, gamma_t, one_minus_t;

  init_mp(tempComp); init_mp(gamma_t); init_mp(one_minus_t);
  init_vec_mp(SS, 0);
  init_mat_mp(Jv_SS, 0, 0);

  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars); // 1 - t
  mul_mp(gamma_t, CD->gamma_mp, pathVars); // gamma * t

  // evaluate the cascade function 'F + RTBx' where T = identity
  standard_codim_cascade_eval_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, ED);

  // evaluate the start system
  change_size_vec_mp(SS, CD->system_rank);
  change_size_mat_mp(Jv_SS, CD->system_rank, CD->new_variables);
  SS->size = Jv_SS->rows = CD->system_rank;
  Jv_SS->cols = CD->new_variables;
  for (i = 0; i < rank; i++)
  {
    if (CD->new_degrees[i] > 0)
    {
      exp_mp_int(&SS->coord[i], &vars->coord[0], CD->new_degrees[i] - 1);
      neg_mp(&SS->coord[i], &SS->coord[i]); // - x_0 ^ (d_i - 1)
      exp_mp_int(tempComp, &vars->coord[i+1], CD->new_degrees[i] - 1); // x_(i+1) ^ (d_i - 1)
    }
    else
    {
      set_zero_mp(&SS->coord[i]);
      set_zero_mp(tempComp);
    }

    // zero out row of Jv since we are only storing to a couple of locations in this row
    for (j = 0; j < Jv_SS->cols; j++)
    {
      set_zero_mp(&Jv_SS->entry[i][j]);
    }

    // d_i * -x_0 ^ (d_i - 1)
    mpf_mul_ui(Jv_SS->entry[i][0].r, SS->coord[i].r, CD->new_degrees[i]);
    mpf_mul_ui(Jv_SS->entry[i][0].i, SS->coord[i].i, CD->new_degrees[i]);
    // d_i * x_(i+1)^(d_i - 1)
    mpf_mul_ui(Jv_SS->entry[i][i+1].r, tempComp->r, CD->new_degrees[i]);
    mpf_mul_ui(Jv_SS->entry[i][i+1].i, tempComp->i, CD->new_degrees[i]);

    // -x_0 ^ d_i
    mul_mp(&SS->coord[i], &SS->coord[i], &vars->coord[0]);
    // x_(i+1) ^ d_i - x_0 ^ d_i
    sum_mul_mp(&SS->coord[i], tempComp, &vars->coord[i+1]);
  }

  // combine everything - sizes are already setup
  // funcVals = funcVals * one_minus_t + gamma_t * SS and 'bottom' is patch - already setup
  // Jv = Jv * one_minus_t + gamma_t * Jv_SS and 'bottom' is patch coefficients - already setup
  // Jp = -funcVals + gamma * SS and 'bottom' is 0 - already setup
  for (i = 0; i < rank; i++)
  { // setup Jp
    mul_mp(&Jp->entry[i][0], CD->gamma_mp, &SS->coord[i]);
    sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], &funcVals->coord[i]);

    // find funcVals
    mul_mp(&funcVals->coord[i], &funcVals->coord[i], one_minus_t);
    sum_mul_mp(&funcVals->coord[i], gamma_t, &SS->coord[i]);

    // setup Jv
    for (j = 0; j < CD->new_variables; j++)
    { // setup Jv[i][j]
      mul_mp(&Jv->entry[i][j], &Jv->entry[i][j], one_minus_t);
      sum_mul_mp(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
    }
  }

  // parVals & parDer are already setup correctly

  clear_mp(tempComp); clear_mp(gamma_t); clear_mp(one_minus_t);
  clear_vec_mp(SS);
  clear_mat_mp(Jv_SS);

  return 0;
}


