// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"
#include "dimbydim.h"

int standard_dimbydim_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for dim-by-dim            *
\***************************************************************/
{
  codim_t *CD = (codim_t *)ED;
  int i, j, k, codim_index = CD->curr_codim_index;
  int codim = CD->codim[codim_index].codim,
      useIntrinsic = CD->codim[codim_index].useIntrinsicSlice;

  mat_d Jv_F, Jv_SS;
  vec_d F, SS;
  comp_d tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_SS, 0, 0);
  init_vec_d(F, CD->num_funcs); init_vec_d(SS, 0);

  // setup one_minus_t & gamma_t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  mul_d(gamma_t, CD->gamma_d, pathVars);

  if (useIntrinsic)
  { // convert vars to the original variables
    point_d orig_vars;
    init_point_d(orig_vars, 0);

    intrinsicToExtrinsic_d(orig_vars, vars, CD->codim[codim_index].B_d, CD->codim[codim_index].p_d);  

    // find the homogenizing coordinate
    set_d(hom_var, CD->codim[codim_index].homVarConst_d);
    for (i = 0; i < CD->new_variables; i++)
    {
      sum_mul_d(hom_var, &CD->codim[codim_index].H_d->coord[i], &orig_vars->coord[i]);
    }

    // evaluate the original functions
    if (CD->new_variables != CD->orig_variables)
    { // convert orig_vars to prog_vars
      point_d prog_vars;
      init_point_d(prog_vars, 0);

      mul_mat_vec_d(prog_vars, CD->C_d, orig_vars);
      // evaluate F
      evalProg_d(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, CD->Prog);
      // convert Jv to orig_vars
      mat_mul_d(Jv, Jv, CD->C_d);

      clear_point_d(prog_vars);
    }
    else
    { // evaluate F
      evalProg_d(funcVals, parVals, parDer, Jv, Jp, orig_vars, pathVars, CD->Prog);
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

    // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = codim;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = codim; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
        if (CD->codim[codim_index].W[i][j-codim] > 0)
        {
          exp_d(tempComp2, hom_var, CD->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
          mul_d(tempComp2, tempComp2, &CD->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mul_rdouble_d(tempComp2, tempComp2, CD->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
        }
        else // W[i][j] == 0
        {
          set_zero_d(tempComp2); // d/d(hom_var) = 0
          set_d(tempComp, &CD->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
        sum_mul_d(&F->coord[i], &F->coord[j], tempComp);
        // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_d(&Jv_F->entry[i][k], tempComp2, &CD->codim[codim_index].H_d->coord[k]);
        }
      }
    }

    // evaluate the start system
    change_size_vec_d(SS, codim);
    change_size_mat_d(Jv_SS, codim, CD->new_variables);
    SS->size = Jv_SS->rows = codim;
    Jv_SS->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    {  
      if (CD->new_degrees[i] > 0)
      {
        exp_d(&SS->coord[i], &orig_vars->coord[0], CD->new_degrees[i] - 1);
        neg_d(&SS->coord[i], &SS->coord[i]); // - x_0 ^ (d_i - 1)
        exp_d(tempComp, &orig_vars->coord[i+1], CD->new_degrees[i] - 1); // x_(i+1) ^ (d_i - 1)
      }
      else
      {
        set_zero_d(&SS->coord[i]);
        set_zero_d(tempComp);
      }

      // zero out row of Jv since we are only storing to a couple of locations in this row
      for (j = 0; j < Jv_SS->cols; j++)
      {
        set_zero_d(&Jv_SS->entry[i][j]);
      }

      // d_i * -x_0 ^ (d_i - 1)
      mul_rdouble_d(&Jv_SS->entry[i][0], &SS->coord[i], CD->new_degrees[i]);
      // d_i * x_(i+1)^(d_i - 1)
      mul_rdouble_d(&Jv_SS->entry[i][i+1], tempComp, CD->new_degrees[i]);

      // -x_0 ^ d_i
      mul_d(&SS->coord[i], &SS->coord[i], &orig_vars->coord[0]);
      // x_(i+1) ^ d_i - x_0 ^ d_i
      sum_mul_d(&SS->coord[i], tempComp, &orig_vars->coord[i+1]);
    }

    // combine everything! 
    // funcVals = F * (1 - t) + gamma * t * SS
    // Jv = Jv_F * (1 - t) + gamma * t * Jv_SS
    // Jp = -F + gamma * SS
    change_size_vec_d(funcVals, codim);
    change_size_mat_d(Jv, codim, CD->new_variables);
    change_size_mat_d(Jp, codim, 1);
    funcVals->size = Jv->rows = Jp->rows = codim;
    Jv->cols = CD->new_variables;
    Jp->cols = 1;
    for (i = 0; i < codim; i++)
    { // funcVals
      mul_d(&funcVals->coord[i], &F->coord[i], one_minus_t);
      sum_mul_d(&funcVals->coord[i], gamma_t, &SS->coord[i]);

      // Jp
      mul_d(&Jp->entry[i][0], CD->gamma_d, &SS->coord[i]);
      sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &F->coord[i]);

      // Jv
      for (j = 0; j < CD->new_variables; j++)
      {
        mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
        sum_mul_d(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
      }
    }

    // convert Jv to the intrinsic variables - B is of the form [[I][B]]
    Jv->rows = Jv->cols = codim;
    for (i = 0; i < codim; i++)
      for (j = 0; j < codim; j++)
        for (k = codim; k < CD->new_variables; k++)
        {
          sum_mul_d(&Jv->entry[i][j], &Jv->entry[i][k], &CD->codim[codim_index].B_d->entry[k][j]);
        }

    clear_vec_d(orig_vars);
  }
  else
  { 
    vec_d linears;
    mat_d Jv_linears;

    init_vec_d(linears, 0);
    init_mat_d(Jv_linears, 0, 0);

    // find the homogenizing coordinate
    set_d(hom_var, CD->codim[codim_index].homVarConst_d);
    for (i = 0; i < CD->new_variables; i++)
    {
      sum_mul_d(hom_var, &CD->codim[codim_index].H_d->coord[i], &vars->coord[i]);
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

    // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = codim;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = codim; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
        if (CD->codim[codim_index].W[i][j-codim] > 0)
        {
          exp_d(tempComp2, hom_var, CD->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
          mul_d(tempComp2, tempComp2, &CD->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mul_rdouble_d(tempComp2, tempComp2, CD->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
        }
        else // W[i][j] == 0
        {
          set_zero_d(tempComp2); // d/d(hom_var) = 0
          set_d(tempComp, &CD->codim[codim_index].A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
        sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

        // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_d(&Jv_F->entry[i][k], tempComp2, &CD->codim[codim_index].H_d->coord[k]);
        }
      }
    }

    // evaluate the start system
    change_size_vec_d(SS, codim);
    change_size_mat_d(Jv_SS, codim, CD->new_variables);
    SS->size = Jv_SS->rows = codim;
    Jv_SS->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
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

      // zero out row of Jv since we are only storing to a couple of locations in this row
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

    // evaluate the linears
    change_size_vec_d(linears, CD->new_variables - codim);
    change_size_mat_d(Jv_linears, CD->new_variables - codim, CD->new_variables);
    linears->size = Jv_linears->rows = CD->new_variables - codim;
    Jv_linears->cols = CD->new_variables;
    for (i = 0; i < linears->size; i++)
      if (i + 1 < linears->size)
      { // the first CD->new_variables - codim - 1 linears are of the form B * vars
        set_zero_d(&linears->coord[i]);
        for (j = 0; j < vars->size; j++)
        { // update linears[i]
          sum_mul_d(&linears->coord[i], &CD->codim[codim_index].B_d->entry[i][j], &vars->coord[j]);
          // update Jv_linears[i][j]
          set_d(&Jv_linears->entry[i][j], &CD->codim[codim_index].B_d->entry[i][j]);
        }
      }
      else
      { // the last linear is the patch p * vars - 1 or p * vars - hom_var

        // setup constant term
        if (CD->PPD.num_var_gp)
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
          sum_mul_d(&linears->coord[i], &CD->codim[codim_index].p_d->coord[j], &vars->coord[j]);
          // update Jv_linears[i][j]
          if (CD->PPD.num_var_gp)
          { // only worry about p
            set_d(&Jv_linears->entry[i][j], &CD->codim[codim_index].p_d->coord[j]);
          }
          else
          { // need to worry about p and H
            sub_d(&Jv_linears->entry[i][j], &CD->codim[codim_index].p_d->coord[j], &CD->codim[codim_index].H_d->coord[j]);
          }
        }
      }

    // combine everything!
    // funcVals = [F * (1 - t) + gamma * t * SS, linears]
    // Jv = [Jv_F * (1 - t) + gamma * t * Jv_SS, Jv_linears]
    // Jp = [-F + gamma * SS, 0]
    change_size_vec_d(funcVals, CD->new_variables);
    change_size_mat_d(Jv, CD->new_variables, CD->new_variables);
    change_size_mat_d(Jp, CD->new_variables, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = CD->new_variables;
    Jp->cols = 1;
    for (i = 0; i < funcVals->size; i++)
      if (i < codim)
      { // funcVals
        mul_d(&funcVals->coord[i], &F->coord[i], one_minus_t);
        sum_mul_d(&funcVals->coord[i], gamma_t, &SS->coord[i]);

        // Jp
        neg_d(&Jp->entry[i][0], &F->coord[i]);
        sum_mul_d(&Jp->entry[i][0], CD->gamma_d, &SS->coord[i]);

        // Jv
        for (j = 0; j < Jv->cols; j++)
        {
          mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
          sum_mul_d(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
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

    clear_vec_d(linears);
    clear_mat_d(Jv_linears);
  }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  clear_mat_d(Jv_F); clear_mat_d(Jv_SS);
  clear_vec_d(F); clear_vec_d(SS);

  return 0;
}

int standard_dimbydim_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for dim-by-dim            *
\***************************************************************/
{
  codim_t *CD = (codim_t *)ED;
  int i, j, k, codim_index = CD->curr_codim_index;
  int codim = CD->codim[codim_index].codim,
      useIntrinsic = CD->codim[codim_index].useIntrinsicSlice;

  mat_mp Jv_F, Jv_SS;
  vec_mp F, SS;
  comp_mp tempComp, tempComp2, one_minus_t, gamma_t, hom_var;

  init_mp(tempComp); init_mp(one_minus_t); init_mp(gamma_t); init_mp(tempComp2); init_mp(hom_var);
  init_vec_mp(F, 0); init_vec_mp(SS, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_SS, 0, 0);

  // setup one_minus_t & gamma_t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  mul_mp(gamma_t, CD->gamma_mp, pathVars);

  if (useIntrinsic)
  { // convert vars to the original variables
    point_mp orig_vars;
    init_point_mp(orig_vars, 0);

    intrinsicToExtrinsic_mp(orig_vars, vars, CD->codim[codim_index].B_mp, CD->codim[codim_index].p_mp);

    // find the homogenizing coordinate
    set_mp(hom_var, CD->codim[codim_index].homVarConst_mp);
    for (i = 0; i < CD->new_variables; i++)
    {
      sum_mul_mp(hom_var, &CD->codim[codim_index].H_mp->coord[i], &orig_vars->coord[i]);
    }

    // evaluate the original functions
    if (CD->new_variables != CD->orig_variables)
    { // convert orig_vars to prog_vars
      point_mp prog_vars;
      init_point_mp(prog_vars, 0);

      mul_mat_vec_mp(prog_vars, CD->C_mp, orig_vars);
      // evaluate F
      evalProg_mp(funcVals, parVals, parDer, Jv, Jp, prog_vars, pathVars, CD->Prog);
      // convert Jv to orig_vars
      mat_mul_mp(Jv, Jv, CD->C_mp);

      clear_point_mp(prog_vars);
    }
    else
    { // evaluate F
      evalProg_mp(funcVals, parVals, parDer, Jv, Jp, orig_vars, pathVars, CD->Prog);
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

    // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = codim;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = codim; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
        if (CD->codim[codim_index].W[i][j-codim] > 0)
        {
          exp_mp_int(tempComp2, hom_var, CD->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
          mul_mp(tempComp2, tempComp2, &CD->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mpf_mul_ui(tempComp2->r, tempComp2->r, CD->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
          mpf_mul_ui(tempComp2->i, tempComp2->i, CD->codim[codim_index].W[i][j-codim]);
        }
        else // W[i][j] == 0
        {
          set_zero_mp(tempComp2); // d/d(hom_var) = 0
          set_mp(tempComp, &CD->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
        sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);
        // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &CD->codim[codim_index].H_mp->coord[k]);
        }
      }
    }

    // evaluate the start system
    change_size_vec_mp(SS, codim);
    change_size_mat_mp(Jv_SS, codim, CD->new_variables);
    SS->size = Jv_SS->rows = codim;
    Jv_SS->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    {
      if (CD->new_degrees[i] > 0)
      {
        exp_mp_int(&SS->coord[i], &orig_vars->coord[0], CD->new_degrees[i] - 1);
        neg_mp(&SS->coord[i], &SS->coord[i]); // - x_0 ^ (d_i - 1)
        exp_mp_int(tempComp, &orig_vars->coord[i+1], CD->new_degrees[i] - 1); // x_(i+1) ^ (d_i - 1)
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
      mul_mp(&SS->coord[i], &SS->coord[i], &orig_vars->coord[0]);
      // x_(i+1) ^ d_i - x_0 ^ d_i
      sum_mul_mp(&SS->coord[i], tempComp, &orig_vars->coord[i+1]);
    }

    // combine everything!
    // funcVals = F * (1 - t) + gamma * t * SS
    // Jv = Jv_F * (1 - t) + gamma * t * Jv_SS
    // Jp = -F + gamma * SS
    change_size_vec_mp(funcVals, codim);
    change_size_mat_mp(Jv, codim, CD->new_variables);
    change_size_mat_mp(Jp, codim, 1);
    funcVals->size = Jv->rows = Jp->rows = codim;
    Jv->cols = CD->new_variables;
    Jp->cols = 1;
    for (i = 0; i < codim; i++)
    { // funcVals
      mul_mp(&funcVals->coord[i], &F->coord[i], one_minus_t);
      sum_mul_mp(&funcVals->coord[i], gamma_t, &SS->coord[i]);

      // Jp
      neg_mp(&Jp->entry[i][0], &F->coord[i]);
      sum_mul_mp(&Jp->entry[i][0], CD->gamma_mp, &SS->coord[i]);

      // Jv
      for (j = 0; j < CD->new_variables; j++)
      {
        mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
        sum_mul_mp(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
      }
    }

    // convert Jv to the intrinsic variables - B is of the form [[I][B]]
    Jv->rows = Jv->cols = codim;
    for (i = 0; i < codim; i++)
      for (j = 0; j < codim; j++)
        for (k = codim; k < CD->new_variables; k++)
        {
          sum_mul_mp(&Jv->entry[i][j], &Jv->entry[i][k], &CD->codim[codim_index].B_mp->entry[k][j]);
        }

    clear_point_mp(orig_vars);
  }
  else
  {
    vec_mp linears;
    mat_mp Jv_linears;

    init_vec_mp(linears, 0);
    init_mat_mp(Jv_linears, 0, 0);

    // find the homogenizing coordinate
    set_mp(hom_var, CD->codim[codim_index].homVarConst_mp);
    for (i = 0; i < CD->new_variables; i++)
    {
      sum_mul_mp(hom_var, &CD->codim[codim_index].H_mp->coord[i], &vars->coord[i]);
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

    // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
    F->size = Jv_F->rows = codim;
    Jv_F->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
    { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
      for (j = codim; j < CD->num_funcs; j++)
      { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
        if (CD->codim[codim_index].W[i][j-codim] > 0)
        {
          exp_mp_int(tempComp2, hom_var, CD->codim[codim_index].W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
          mul_mp(tempComp2, tempComp2, &CD->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
          mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

          mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
          mpf_mul_ui(tempComp2->r, tempComp2->r, CD->codim[codim_index].W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
          mpf_mul_ui(tempComp2->i, tempComp2->i, CD->codim[codim_index].W[i][j-codim]);
        }
        else // W[i][j] == 0
        {
          set_zero_mp(tempComp2); // d/d(hom_var) = 0
          set_mp(tempComp, &CD->codim[codim_index].A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
        }

        // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
        sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);
        // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
        for (k = 0; k < CD->new_variables; k++)
        {
          sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
          sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &CD->codim[codim_index].H_mp->coord[k]);
        }
      }
    }

    // evaluate the start system
    change_size_vec_mp(SS, codim);
    change_size_mat_mp(Jv_SS, codim, CD->new_variables);
    SS->size = Jv_SS->rows = codim;
    Jv_SS->cols = CD->new_variables;
    for (i = 0; i < codim; i++)
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

    // evaluate the linears
    change_size_vec_mp(linears, CD->new_variables - codim);
    change_size_mat_mp(Jv_linears, CD->new_variables - codim, CD->new_variables);
    linears->size = Jv_linears->rows = CD->new_variables - codim;
    Jv_linears->cols = CD->new_variables;
    for (i = 0; i < linears->size; i++)
      if (i + 1 < linears->size)
      { // the first CD->new_variables - codim - 1 linears are of the form B * vars
        set_zero_mp(&linears->coord[i]);
        for (j = 0; j < vars->size; j++)
        { // update linears[i]
          sum_mul_mp(&linears->coord[i], &CD->codim[codim_index].B_mp->entry[i][j], &vars->coord[j]);
          // update Jv_linears[i][j]
          set_mp(&Jv_linears->entry[i][j], &CD->codim[codim_index].B_mp->entry[i][j]);
        }
      }
      else
      { // the last linear is the patch p * vars - 1 or p * vars - hom_var

        // setup constant term
        if (CD->PPD.num_var_gp)
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
          sum_mul_mp(&linears->coord[i], &CD->codim[codim_index].p_mp->coord[j], &vars->coord[j]);
          // update Jv_linears[i][j]
          if (CD->PPD.num_var_gp)
          { // only worry about p
            set_mp(&Jv_linears->entry[i][j], &CD->codim[codim_index].p_mp->coord[j]);
          }
          else
          { // need to worry about p and H
            sub_mp(&Jv_linears->entry[i][j], &CD->codim[codim_index].p_mp->coord[j], &CD->codim[codim_index].H_mp->coord[j]);
          }
        }
      }

    // combine everything!
    // funcVals = [F * (1 - t) + gamma * t * SS, linears]
    // Jv = [Jv_F * (1 - t) + gamma * t * Jv_SS, Jv_linears]
    // Jp = [-F + gamma * SS, 0]
    change_size_vec_mp(funcVals, CD->new_variables);
    change_size_mat_mp(Jv, CD->new_variables, CD->new_variables);
    change_size_mat_mp(Jp, CD->new_variables, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = CD->new_variables;
    Jp->cols = 1;
    for (i = 0; i < funcVals->size; i++)
      if (i < codim)
      { // funcVals
        mul_mp(&funcVals->coord[i], &F->coord[i], one_minus_t);
        sum_mul_mp(&funcVals->coord[i], gamma_t, &SS->coord[i]);

        // Jp
        neg_mp(&Jp->entry[i][0], &F->coord[i]);
        sum_mul_mp(&Jp->entry[i][0], CD->gamma_mp, &SS->coord[i]);

        // Jv
        for (j = 0; j < Jv->cols; j++)
        {
          mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
          sum_mul_mp(&Jv->entry[i][j], gamma_t, &Jv_SS->entry[i][j]);
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

    clear_vec_mp(linears);
    clear_mat_mp(Jv_linears);
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]); // ds/dt = 1

  clear_mp(tempComp); clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(tempComp2); clear_mp(hom_var);
  clear_vec_mp(F); clear_vec_mp(SS);
  clear_mat_mp(Jv_F); clear_mat_mp(Jv_SS);

  return 0;
}

