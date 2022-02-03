// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"

int general_slice_moving_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: eval for moving of general slices                      *
\***************************************************************/
{
  general_slice_moving_t *M = (general_slice_moving_t *)ED;
  int i, j, k, codim = M->curr_codim, numVars = M->orig_variables, numFuncs = M->num_funcs;
  mat_d Jv_F, Jv_linears_start, Jv_linears_end;
  vec_d F, linears_start, linears_end;
  comp_d tempComp, tempComp2, hom_var;

  set_zero_d(tempComp); set_zero_d(tempComp2); set_zero_d(hom_var);
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_linears_start, 0, 0); init_mat_d(Jv_linears_end, 0, 0);
  init_vec_d(F, 0); init_vec_d(linears_start, 0); init_vec_d(linears_end, 0);

  // find the homogenizing coordinate
  set_d(hom_var, M->homVarConst_d);
  for (i = 0; i < numVars; i++)
  {
    sum_mul_d(hom_var, &M->H_d->coord[i], &vars->coord[i]);
  }

  // evaluate the function
  evalProg_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, M->Prog);

  // setup F & Jv_F using P
  change_size_vec_d(F, numFuncs);
  change_size_mat_d(Jv_F, numFuncs, numVars);
  F->size = Jv_F->rows = numFuncs;
  Jv_F->cols = numVars;
  for (i = 0; i < numFuncs; i++)
  { // put F(i) = f(P[i])
    set_d(&F->coord[i], &funcVals->coord[M->P[i]]);
    for (j = 0; j < numVars; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[M->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = numVars;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < numFuncs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (M->W[i][j-codim] > 0)
      {
        exp_d(tempComp2, hom_var, M->W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_d(tempComp2, tempComp2, &M->A_d->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_d(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_d(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_d(tempComp2, tempComp2, M->W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_d(tempComp2); // d/d(hom_var) = 0
        set_d(tempComp, &M->A_d->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_d(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < numVars; k++)
      {
        sum_mul_d(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_d(&Jv_F->entry[i][k], tempComp2, &M->H_d->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_d(linears_start, numVars - codim);
  change_size_vec_d(linears_end, numVars - codim);
  change_size_mat_d(Jv_linears_start, numVars - codim, numVars);
  change_size_mat_d(Jv_linears_end, numVars - codim, numVars);
  linears_start->size = linears_end->size = Jv_linears_start->rows = Jv_linears_end->rows = numVars - codim;
  Jv_linears_start->cols = Jv_linears_end->cols = numVars;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first numVars - codim - 1 linears are of the form B * vars
      set_zero_d(&linears_start->coord[i]);
      set_zero_d(&linears_end->coord[i]);
      for (j = 0; j < numVars; j++)
      { // update linears_start[i] & linears_end[i]
        sum_mul_d(&linears_start->coord[i], &M->B_start_d->entry[i][j], &vars->coord[j]);
        sum_mul_d(&linears_end->coord[i], &M->B_end_d->entry[i][j], &vars->coord[j]);
        // update Jv_linears_start[i][j] & Jv_linears_end[i][j]
        set_d(&Jv_linears_start->entry[i][j], &M->B_start_d->entry[i][j]);
        set_d(&Jv_linears_end->entry[i][j], &M->B_end_d->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var
      // setup constant term
      if (M->PPD.num_var_gp)
      { // == 1
        set_one_d(&linears_start->coord[i]);
        set_one_d(&linears_end->coord[i]);
      }
      else
      { // == hom_var
        set_d(&linears_start->coord[i], hom_var);
        set_d(&linears_end->coord[i], hom_var);
      }
      neg_d(&linears_start->coord[i], &linears_start->coord[i]);
      neg_d(&linears_end->coord[i], &linears_end->coord[i]);
      for (j = 0; j < numVars; j++)
      { // update linears_start[i] & linears_end[i]
        sum_mul_d(&linears_start->coord[i], &M->patch_start_d->coord[j], &vars->coord[j]);
        sum_mul_d(&linears_end->coord[i], &M->patch_end_d->coord[j], &vars->coord[j]);
        // update Jv_linears[i][j]
        if (M->PPD.num_var_gp)
        { // only worry about p
          set_d(&Jv_linears_start->entry[i][j], &M->patch_start_d->coord[j]);
          set_d(&Jv_linears_end->entry[i][j], &M->patch_end_d->coord[j]);
        }
        else
        { // need to worry about p and H
          sub_d(&Jv_linears_start->entry[i][j], &M->patch_start_d->coord[j], &M->H_d->coord[j]);
          sub_d(&Jv_linears_end->entry[i][j], &M->patch_end_d->coord[j], &M->H_d->coord[j]);
        }
      }
    }
  // combine everything!
  // funcVals = [F, linears_start*gamma*t + linears_end*(1-t)]
  // Jv = [Jv_F, Jv_linears_end*gamma*t + Jv_linears_end*(1-t)]
  // Jp = [0, linears_start*gamma - linears_end]
  change_size_point_d(funcVals, numVars);
  change_size_mat_d(Jv, numVars, numVars);
  change_size_mat_d(Jp, numVars, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = numVars;
  Jp->cols = 1;
  mul_d(tempComp, M->gamma_d, pathVars); // gamma*t
  set_one_d(tempComp2); sub_d(tempComp2, tempComp2, pathVars) ; // 1 - t
  for (i = 0; i < numVars; i++)
    if (i < codim)
    { // funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);
      // Jp
      set_zero_d(&Jp->entry[i][0]);
      // Jv
      for (j = 0; j < numVars; j++)
      {
        set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else
    { // funcVals
      mul_d(&funcVals->coord[i], &linears_start->coord[i-codim], tempComp);
      sum_mul_d(&funcVals->coord[i], &linears_end->coord[i-codim], tempComp2);
      // Jp
      mul_d(&Jp->entry[i][0], &linears_start->coord[i-codim], M->gamma_d);
      sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &linears_end->coord[i-codim]);
      // Jv
      for (j = 0; j < numVars; j++)
      {
        mul_d(&Jv->entry[i][j], &Jv_linears_start->entry[i-codim][j], tempComp);
        sum_mul_d(&Jv->entry[i][j], &Jv_linears_end->entry[i-codim][j], tempComp2);
      }
    }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_d(Jv_F); clear_mat_d(Jv_linears_start); clear_mat_d(Jv_linears_end);
  clear_vec_d(F); clear_vec_d(linears_start); clear_vec_d(linears_end);
  
  return 0;
}

int general_slice_moving_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: eval for moving of general slices                      *
\***************************************************************/
{
  general_slice_moving_t *M = (general_slice_moving_t *)ED;
  int i, j, k, codim = M->curr_codim, numVars = M->orig_variables, numFuncs = M->num_funcs;
  mat_mp Jv_F, Jv_linears_start, Jv_linears_end;
  vec_mp F, linears_start, linears_end;
  comp_mp tempComp, tempComp2, hom_var;

  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_linears_start, 0, 0); init_mat_mp(Jv_linears_end, 0, 0);
  init_vec_mp(F, 0); init_vec_mp(linears_start, 0); init_vec_mp(linears_end, 0);
  init_mp(tempComp); init_mp(tempComp2); init_mp(hom_var);
  set_zero_mp(tempComp); set_zero_mp(tempComp2); set_zero_mp(hom_var);

  // find the homogenizing coordinate
  set_mp(hom_var, M->homVarConst_mp);
  for (i = 0; i < numVars; i++)
  {
    sum_mul_mp(hom_var, &M->H_mp->coord[i], &vars->coord[i]);
  }

  // evaluate the function
  evalProg_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, M->Prog);

  // setup F & Jv_F using P
  change_size_vec_mp(F, numFuncs);
  change_size_mat_mp(Jv_F, numFuncs, numVars);
  F->size = Jv_F->rows = numFuncs;
  Jv_F->cols = numVars;
  for (i = 0; i < numFuncs; i++)
  { // put F(i) = f(P[i])
    set_mp(&F->coord[i], &funcVals->coord[M->P[i]]);
    for (j = 0; j < numVars; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[M->P[i]][j]);
    }
  }

  // find F = [I A.*W]F and Jv_F = [I A.*W]Jv_F (i.e. randomize to the correct size)
  F->size = Jv_F->rows = codim;
  Jv_F->cols = numVars;
  for (i = 0; i < codim; i++)
  { // add on the randomized part (i.e [0 A.*W]F & [0 A.*W]Jv_F)
    for (j = codim; j < numFuncs; j++)
    { // find F[j]*W[i][j]*A[i][j]*hom_var^(W[i][j] - 1) and find A[i][j]*hom_var^W[i][j]
      if (M->W[i][j-codim] > 0)
      {
        exp_mp_int(tempComp2, hom_var, M->W[i][j-codim] - 1); // hom_var ^ (W[i][j] - 1)
        mul_mp(tempComp2, tempComp2, &M->A_mp->entry[i][j-codim]); // A[i][j] * hom_var^(W[i][j] - 1)
        mul_mp(tempComp, tempComp2, hom_var); // A[i][j] * hom_var ^ W[i][j]

        mul_mp(tempComp2, tempComp2, &F->coord[j]); // A[i][j]*hom_var^(W[i][j] - 1)*F[j]
        mul_rdouble_mp(tempComp2, tempComp2, M->W[i][j-codim]); // W[i][j]*A[i][j]*hom_var^([i][j]-1)*F[j]
      }
      else // W[i][j] == 0
      {
        set_zero_mp(tempComp2); // d/d(hom_var) = 0
        set_mp(tempComp, &M->A_mp->entry[i][j-codim]); // A[i][j] * hom_var ^ 0 = A[i][j]
      }

      // F[i] += F[j] * A[i][j] * hom_var^W[i][j]
      sum_mul_mp(&F->coord[i], &F->coord[j], tempComp);

      // Jv_F[i][k] += A[i][j] * hom_var^W[i][j] * Jv_F[j][k] + A[i][j]*W[i][j]*hom_var^(W[i][j] - 1)*F[j]* d(hom_var)/dx_k
      for (k = 0; k < numVars; k++)
      {
        sum_mul_mp(&Jv_F->entry[i][k], tempComp, &Jv_F->entry[j][k]);
        sum_mul_mp(&Jv_F->entry[i][k], tempComp2, &M->H_mp->coord[k]);
      }
    }
  }

  // evaluate the linears
  change_size_vec_mp(linears_start, numVars - codim);
  change_size_vec_mp(linears_end, numVars - codim);
  change_size_mat_mp(Jv_linears_start, numVars - codim, numVars);
  change_size_mat_mp(Jv_linears_end, numVars - codim, numVars);
  linears_start->size = linears_end->size = Jv_linears_start->rows = Jv_linears_end->rows = numVars - codim;
  Jv_linears_start->cols = Jv_linears_end->cols = numVars;
  for (i = 0; i < linears_start->size; i++)
    if (i + 1 < linears_start->size)
    { // the first numVars - codim - 1 linears are of the form B * vars
      set_zero_mp(&linears_start->coord[i]);
      set_zero_mp(&linears_end->coord[i]);
      for (j = 0; j < numVars; j++)
      { // update linears_start[i] & linears_end[i]
        sum_mul_mp(&linears_start->coord[i], &M->B_start_mp->entry[i][j], &vars->coord[j]);
        sum_mul_mp(&linears_end->coord[i], &M->B_end_mp->entry[i][j], &vars->coord[j]);
        // update Jv_linears_start[i][j] & Jv_linears_end[i][j]
        set_mp(&Jv_linears_start->entry[i][j], &M->B_start_mp->entry[i][j]);
        set_mp(&Jv_linears_end->entry[i][j], &M->B_end_mp->entry[i][j]);
      }
    }
    else
    { // the last linear is the patch p * vars - 1 or p * vars - hom_var
      // setup constant term
      if (M->PPD.num_var_gp)
      { // == 1
        set_one_mp(&linears_start->coord[i]);
        set_one_mp(&linears_end->coord[i]);
      }
      else
      { // == hom_var
        set_mp(&linears_start->coord[i], hom_var);
        set_mp(&linears_end->coord[i], hom_var);
      }
      neg_mp(&linears_start->coord[i], &linears_start->coord[i]);
      neg_mp(&linears_end->coord[i], &linears_end->coord[i]);
      for (j = 0; j < numVars; j++)
      { // update linears_start[i] & linears_end[i]
        sum_mul_mp(&linears_start->coord[i], &M->patch_start_mp->coord[j], &vars->coord[j]);
        sum_mul_mp(&linears_end->coord[i], &M->patch_end_mp->coord[j], &vars->coord[j]);
        // update Jv_linears[i][j]
        if (M->PPD.num_var_gp)
        { // only worry about p
          set_mp(&Jv_linears_start->entry[i][j], &M->patch_start_mp->coord[j]);
          set_mp(&Jv_linears_end->entry[i][j], &M->patch_end_mp->coord[j]);
        }
        else
        { // need to worry about p and H
          sub_mp(&Jv_linears_start->entry[i][j], &M->patch_start_mp->coord[j], &M->H_mp->coord[j]);
          sub_mp(&Jv_linears_end->entry[i][j], &M->patch_end_mp->coord[j], &M->H_mp->coord[j]);
        }
      }
    }
  // combine everything!
  // funcVals = [F, linears_start*gamma*t + linears_end*(1-t)]
  // Jv = [Jv_F, Jv_linears_end*gamma*t + Jv_linears_end*(1-t)]
  // Jp = [0, linears_start*gamma - linears_end]
  change_size_point_mp(funcVals, numVars);
  change_size_mat_mp(Jv, numVars, numVars);
  change_size_mat_mp(Jp, numVars, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = numVars;
  Jp->cols = 1;
  mul_mp(tempComp, M->gamma_mp, pathVars); // gamma*t
  set_one_mp(tempComp2); sub_mp(tempComp2, tempComp2, pathVars) ; // 1 - t
  for (i = 0; i < numVars; i++)
    if (i < codim)
    { // funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);
      // Jp
      set_zero_mp(&Jp->entry[i][0]);
      // Jv
      for (j = 0; j < numVars; j++)
      {
        set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
      }
    }
    else
    { // funcVals
      mul_mp(&funcVals->coord[i], &linears_start->coord[i-codim], tempComp);
      sum_mul_mp(&funcVals->coord[i], &linears_end->coord[i-codim], tempComp2);
      // Jp
      mul_mp(&Jp->entry[i][0], &linears_start->coord[i-codim], M->gamma_mp);
      sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], &linears_end->coord[i-codim]);
      // Jv
      for (j = 0; j < numVars; j++)
      {
        mul_mp(&Jv->entry[i][j], &Jv_linears_start->entry[i-codim][j], tempComp);
        sum_mul_mp(&Jv->entry[i][j], &Jv_linears_end->entry[i-codim][j], tempComp2);
      }
    }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);        // ds/dt = 1

  clear_mat_mp(Jv_F); clear_mat_mp(Jv_linears_start); clear_mat_mp(Jv_linears_end);
  clear_vec_mp(F); clear_vec_mp(linears_start); clear_vec_mp(linears_end);
  clear_mp(tempComp); clear_mp(tempComp2); clear_mp(hom_var);

  return 0;
}

void initialize_setup_general_slice(general_slice_moving_t *M, witness_t *W_start, int start_codim_index, witness_t *W_end, int end_codim_index, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: initialize and setup to move a general slice           *
\***************************************************************/
{
  int i, j, prec;

  // verify that the codimensions are the same
  if (W_start->codim[start_codim_index].codim != W_end->codim[end_codim_index].codim)
  {
    printf("\nERROR: The codimensions are not the same!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // copy SLP
  M->Prog = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(M->Prog, W_start->Prog);
  // copy PPD
  cp_preproc_data(&M->PPD, &W_start->PPD);
  // copy codim
  M->curr_codim = W_start->codim[start_codim_index].codim;
  // copy orig_variables
  M->orig_variables = W_start->orig_variables;
  // copy num_funcs
  M->num_funcs = W_start->num_funcs;
  // copy curr_precision
  M->curr_precision = W_start->curr_precision;
  // copy P
  M->P = (int *)bmalloc(M->num_funcs * sizeof(int));
  for (i = 0; i < M->num_funcs; i++)
    M->P[i] = W_start->P[i];

  // setup structures
  if (MPType == 0)
  { // setup gamma
    set_d(M->gamma_d, W_start->gamma_d);
    // setup A
    init_mat_d(M->A_d, 0, 0);
    mat_cp_d(M->A_d, W_start->codim[start_codim_index].A_d);
    M->A_rows = M->A_d->rows;
    M->A_cols = M->A_d->cols;
    // setup W
    M->W = (int **)bmalloc(M->A_rows * sizeof(int *));
    for (i = 0; i < M->A_rows; i++)
    {
      M->W[i] = (int *)bmalloc(M->A_cols * sizeof(int));
      for (j = 0; j < M->A_cols; j++)
        M->W[i][j] = W_start->codim[start_codim_index].W[i][j];
    }
    // setup H
    init_vec_d(M->H_d, 0);
    vec_cp_d(M->H_d, W_start->codim[start_codim_index].H_d);
    // setup homVarConst
    set_d(M->homVarConst_d, W_start->codim[start_codim_index].homVarConst_d);
    // setup B_start
    init_mat_d(M->B_start_d, 0, 0);
    mat_cp_d(M->B_start_d, W_start->codim[start_codim_index].B_d);
    // setup B_end
    init_mat_d(M->B_end_d, 0, 0);
    mat_cp_d(M->B_end_d, W_end->codim[end_codim_index].B_d);
    // setup patch_start
    init_vec_d(M->patch_start_d, 0);
    vec_cp_d(M->patch_start_d, W_start->codim[start_codim_index].p_d);
    // setup patch_end
    init_vec_d(M->patch_end_d, 0);
    vec_cp_d(M->patch_end_d, W_end->codim[end_codim_index].p_d);
  }
  else if (MPType == 1)
  { // setup gamma
    init_mp(M->gamma_mp);
    set_mp(M->gamma_mp, W_start->gamma_mp);
    // setup A
    init_mat_mp(M->A_mp, 0, 0);
    mat_cp_mp(M->A_mp, W_start->codim[start_codim_index].A_mp);
    M->A_rows = M->A_mp->rows;
    M->A_cols = M->A_mp->cols;
    // setup W
    M->W = (int **)bmalloc(M->A_rows * sizeof(int *));
    for (i = 0; i < M->A_rows; i++)
    {
      M->W[i] = (int *)bmalloc(M->A_cols * sizeof(int));
      for (j = 0; j < M->A_cols; j++)
        M->W[i][j] = W_start->codim[start_codim_index].W[i][j];
    }
    // setup H
    init_vec_mp(M->H_mp, 0);
    vec_cp_mp(M->H_mp, W_start->codim[start_codim_index].H_mp);
    // setup homVarConst
    init_mp(M->homVarConst_mp);
    set_mp(M->homVarConst_mp, W_start->codim[start_codim_index].homVarConst_mp);
    // setup B_start
    init_mat_mp(M->B_start_mp, 0, 0);
    mat_cp_mp(M->B_start_mp, W_start->codim[start_codim_index].B_mp);
    // setup B_end
    init_mat_mp(M->B_end_mp, 0, 0);
    mat_cp_mp(M->B_end_mp, W_end->codim[end_codim_index].B_mp);
    // setup patch_start
    init_vec_mp(M->patch_start_mp, 0);
    vec_cp_mp(M->patch_start_mp, W_start->codim[start_codim_index].p_mp);
    // setup patch_end
    init_vec_mp(M->patch_end_mp, 0);
    vec_cp_mp(M->patch_end_mp, W_end->codim[end_codim_index].p_mp);
  }
  else
  { // setup gamma
    init_mp2(M->gamma_mp, M->curr_precision);
    init_rat(M->gamma_rat);
    set_d(M->gamma_d, W_start->gamma_d);
    set_mp(M->gamma_mp, W_start->gamma_mp);
    set_rat(M->gamma_rat, W_start->gamma_rat);
    // setup A
    M->A_rows = W_start->codim[start_codim_index].A_rows;
    M->A_cols = W_start->codim[start_codim_index].A_cols;
    init_mat_d(M->A_d, M->A_rows, M->A_cols);
    init_mat_mp2(M->A_mp, M->A_rows, M->A_cols, M->curr_precision);
    init_mat_rat(M->A_rat, M->A_rows, M->A_cols);
    mat_cp_d(M->A_d, W_start->codim[start_codim_index].A_d);
    mat_cp_mp(M->A_mp, W_start->codim[start_codim_index].A_mp);
    mat_cp_rat(M->A_rat, W_start->codim[start_codim_index].A_rat, M->A_rows, M->A_cols);
    // setup W
    M->W = (int **)bmalloc(M->A_rows * sizeof(int *));
    for (i = 0; i < M->A_rows; i++)
    {
      M->W[i] = (int *)bmalloc(M->A_cols * sizeof(int));
      for (j = 0; j < M->A_cols; j++)
        M->W[i][j] = W_start->codim[start_codim_index].W[i][j];
    }
    // setup H
    init_vec_d(M->H_d, 0); 
    init_vec_mp2(M->H_mp, 0, M->curr_precision);
    init_vec_rat(M->H_rat, W_start->codim[start_codim_index].H_d->size);
    vec_cp_d(M->H_d, W_start->codim[start_codim_index].H_d);
    vec_cp_mp(M->H_mp, W_start->codim[start_codim_index].H_mp);
    vec_cp_rat(M->H_rat, W_start->codim[start_codim_index].H_rat, M->H_d->size);
    // setup homVarConst
    init_mp2(M->homVarConst_mp, M->curr_precision);
    init_rat(M->homVarConst_rat);
    set_d(M->homVarConst_d, W_start->codim[start_codim_index].homVarConst_d);
    set_mp(M->homVarConst_mp, W_start->codim[start_codim_index].homVarConst_mp);
    set_rat(M->homVarConst_rat, W_start->codim[start_codim_index].homVarConst_rat);
    // setup B_start
    init_mat_d(M->B_start_d, 0, 0);
    init_mat_mp2(M->B_start_mp, 0, 0, M->curr_precision);
    init_mat_rat(M->B_start_rat, W_start->codim[start_codim_index].B_d->rows, W_start->codim[start_codim_index].B_d->cols);
    mat_cp_d(M->B_start_d, W_start->codim[start_codim_index].B_d);
    mat_cp_mp(M->B_start_mp, W_start->codim[start_codim_index].B_mp);
    mat_cp_rat(M->B_start_rat, W_start->codim[start_codim_index].B_rat, M->B_start_d->rows, M->B_start_d->cols);
    // setup B_end
    init_mat_d(M->B_end_d, 0, 0);
    init_mat_mp2(M->B_end_mp, 0, 0, M->curr_precision);
    init_mat_rat(M->B_end_rat, W_end->codim[end_codim_index].B_d->rows, W_end->codim[end_codim_index].B_d->cols); 
    mat_cp_d(M->B_end_d, W_end->codim[end_codim_index].B_d);
    mat_cp_mp(M->B_end_mp, W_end->codim[end_codim_index].B_mp);
    mat_cp_rat(M->B_end_rat, W_end->codim[end_codim_index].B_rat, M->B_end_d->rows, M->B_end_d->cols);
    // setup patch_start
    init_vec_d(M->patch_start_d, 0);
    init_vec_mp2(M->patch_start_mp, 0, M->curr_precision);
    init_vec_rat(M->patch_start_rat, W_start->codim[start_codim_index].p_d->size);
    vec_cp_d(M->patch_start_d, W_start->codim[start_codim_index].p_d);
    vec_cp_mp(M->patch_start_mp, W_start->codim[start_codim_index].p_mp);
    vec_cp_rat(M->patch_start_rat, W_start->codim[start_codim_index].p_rat, M->patch_start_d->size);
    // setup patch_end
    init_vec_d(M->patch_end_d, 0);
    init_vec_mp2(M->patch_end_mp, 0, M->curr_precision);
    init_vec_rat(M->patch_end_rat, W_end->codim[end_codim_index].p_d->size);
    vec_cp_d(M->patch_end_d, W_end->codim[end_codim_index].p_d);
    vec_cp_mp(M->patch_end_mp, W_end->codim[end_codim_index].p_mp);
    vec_cp_rat(M->patch_end_rat, W_end->codim[end_codim_index].p_rat, M->patch_end_d->size);

    // make sure everything is in the proper precision
    prec = M->curr_precision;
    M->curr_precision--;
    change_general_slice_prec(M, prec);
  }

  return;
}

int change_general_slice_prec(void const *ED, int prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change precision for general slice mover               *
\***************************************************************/
{
  general_slice_moving_t *M = (general_slice_moving_t *)ED;

  // set the SLP to the correct precision
  M->Prog->precision = prec;

  if (M->curr_precision != prec)
  { // change precision
    M->curr_precision = prec;
    // gamma
    change_prec_mp_rat(M->gamma_mp, prec, M->gamma_rat);
    // A
    change_prec_mat_mp_rat(M->A_mp, prec, M->A_rat);
    // H
    change_prec_vec_mp_rat(M->H_mp, prec, M->H_rat);
    // homVarConst
    change_prec_mp_rat(M->homVarConst_mp, prec, M->homVarConst_rat);
    // B_start
    change_prec_mat_mp_rat(M->B_start_mp, prec, M->B_start_rat);
    // B_end
    change_prec_mat_mp_rat(M->B_end_mp, prec, M->B_end_rat);
    // patch_start
    change_prec_vec_mp_rat(M->patch_start_mp, prec, M->patch_start_rat);
    // patch_end
    change_prec_vec_mp_rat(M->patch_end_mp, prec, M->patch_end_rat);

    // set precision
    M->curr_precision = prec;
  }

  return 0;
}

void clear_general_slice(general_slice_moving_t *M, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear general slice mover                              *
\***************************************************************/
{
  int i;

  // clear Prog
  clearProg(M->Prog, MPType, 0);
  free(M->Prog);
  // clear PPD
  preproc_data_clear(&M->PPD);
  // clear P
  free(M->P);
  // clear W
  for (i = 0; i < M->A_rows; i++)
    free(M->W[i]);
  free(M->W);
  // clear gamma
  clear_d_mp_rat(M->gamma_d, M->gamma_mp, M->gamma_rat, MPType);
  // clear A
  clear_mat(M->A_d, M->A_mp, M->A_rat, MPType);
  // clear H
  clear_vec(M->H_d, M->H_mp, M->H_rat, MPType);
  // clear homVarConst
  clear_d_mp_rat(M->homVarConst_d, M->homVarConst_mp, M->homVarConst_rat, MPType);
  // clear B_start
  clear_mat(M->B_start_d, M->B_start_mp, M->B_start_rat, MPType);
  // clear B_end
  clear_mat(M->B_end_d, M->B_end_mp, M->B_end_rat, MPType);
  // clear patch_start
  clear_vec(M->patch_start_d, M->patch_start_mp, M->patch_start_rat, MPType);
  // clear patch_end
  clear_vec(M->patch_end_d, M->patch_end_mp, M->patch_end_rat, MPType);

  return; 
}

int general_slice_moving_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
 * NOTES: compute the dehom point                                *
\***************************************************************/
{
  general_slice_moving_t *M = NULL;
  *out_prec = in_prec;

  if (in_prec < 64)
  { // compute out_d
    M = (general_slice_moving_t *)ED_d;
    getDehomPoint_d(out_d, in_d, in_d->size, &M->PPD);
  }
  else
  { // compute out_mp
    M = (general_slice_moving_t *)ED_mp;
    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
    getDehomPoint_mp(out_mp, in_mp, in_mp->size, &M->PPD);
  }

  M = NULL;

  return 0;
}

int general_slice_moving_track(endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, endgame_data_t *endGame, general_slice_moving_t *M, point_d startPt_d, point_mp startPt_mp, int startPt_prec, int pathNum_startPt, int checkSharpen, tracker_config_t *T, FILE *OUT, FILE *MIDOUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: whether slice moving was a success or not      *
* NOTES: moves the slice                                        *
\***************************************************************/
{
  int retVal = 0, rankType = 1;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  int (*change_prec)(void const *, int) = NULL;
  int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

  ptr_to_eval_d = &general_slice_moving_eval_d;
  ptr_to_eval_mp = &general_slice_moving_eval_mp;
  change_prec = &change_general_slice_prec;
  find_dehom = &general_slice_moving_dehom;

  // setup for tracking
  T->first_step_of_path = 1;
  T->endgameOnly = 0;

  if (T->MPType == 0)
  { // track the path in double precision
    point_data_d startPt;
    init_point_data_d(&startPt, startPt_d->size);

    // setup for tracking the path
    point_cp_d(startPt.point, startPt_d);
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, M, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(pathNum_startPt, rankType, NULL, &endPt_d->corank, &endPt_d->smallest_nonzero_SV, &endPt_d->largest_zero_SV, endGame, &startPt, OUT, MIDOUT, T, M, M, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // see if we should look to sharpen this endPt
    if (checkSharpen && endGame->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endGame, T, OUT, M, M, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer for the path and find retVal
    retVal = printMembershipFooter_d(&endGame->PD_d, endGame->condition_number, endGame->function_residual_d, endGame->latest_newton_residual_d, endGame->t_val_at_latest_sample_point_d, endGame->error_at_latest_sample_point_d, OUT, endGame->retVal, T);

    // finishing setup endPt_d & endGame correctly
    point_cp_d(endPt_d->endPt, endGame->PD_d.point);
    point_cp_d(endPt_d->last_approx, endGame->last_approx_d);
    set_d(endPt_d->finalT, endGame->PD_d.time);
    endPt_d->cond_num = endGame->condition_number;
    endPt_d->retVal = endGame->retVal = retVal;

    clear_point_data_d(&startPt);
  }
  else if (T->MPType == 1)
  { // track the path in fixed multi precision
    point_data_mp startPt;
    init_point_data_mp(&startPt, startPt_mp->size);

    // setup for tracking the path
    point_cp_mp(startPt.point, startPt_mp);
    set_one_mp(startPt.time);

    // print the header for the path
    printPathHeader_mp(OUT, &startPt, T, pathNum_startPt, M, ptr_to_eval_mp);

    // track the path
    zero_dim_track_path_rank_mp(pathNum_startPt, rankType, NULL, &endPt_mp->corank, &endPt_mp->smallest_nonzero_SV, &endPt_mp->largest_zero_SV, endGame, &startPt, OUT, MIDOUT, T, M, ptr_to_eval_mp, find_dehom);

    // see if we should look to sharpen this endPt
    if (checkSharpen && endGame->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endGame, T, OUT, M, M, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer in multi precision and find retVal
    retVal = printMembershipFooter_mp(&endGame->PD_mp, endGame->condition_number, endGame->first_increase, endGame->function_residual_mp, endGame->latest_newton_residual_mp, endGame->t_val_at_latest_sample_point_mp, endGame->error_at_latest_sample_point_mp, OUT, endGame->retVal, T);

    // finishing setup endPt_mp & endGame correctly
    point_cp_mp(endPt_mp->endPt, endGame->PD_mp.point);
    point_cp_mp(endPt_mp->last_approx, endGame->last_approx_mp);
    set_mp(endPt_mp->finalT, endGame->PD_mp.time);
    endPt_mp->cond_num = endGame->condition_number;
    endPt_mp->retVal = endGame->retVal = retVal;

    clear_point_data_mp(&startPt);
  }
  else
  { // track the path using AMP
    point_data_d startPt;
    init_point_data_d(&startPt, 0);

    // setup for tracking the path
    if (startPt_prec < 64)
    {
      point_cp_d(startPt.point, startPt_d);
    }
    else
    {
      point_mp_to_d(startPt.point, startPt_mp);
    }
    set_one_d(startPt.time);

    // print the header for the path
    printPathHeader_d(OUT, &startPt, T, pathNum_startPt, M, ptr_to_eval_d);

    // track the path
    zero_dim_track_path_rank_d(pathNum_startPt, rankType, NULL, &endPt_amp->corank, &endPt_amp->smallest_nonzero_SV, &endPt_amp->largest_zero_SV, endGame, &startPt, OUT, MIDOUT, T, M, M, ptr_to_eval_d, ptr_to_eval_mp, change_prec, find_dehom);

    // see if we should look to sharpen this endPt
    if (checkSharpen && endGame->retVal == 0 && T->sharpenDigits > 0)
    { // use the sharpener for after an endgame
      sharpen_endpoint_endgame(endGame, T, OUT, M, M, ptr_to_eval_d, ptr_to_eval_mp, change_prec);
    }

    // print the footer for the path and find retVal
    if (endGame->prec < 64)
    { // print the footer in double precision and find retVal
      retVal = printMembershipFooter_d(&endGame->PD_d, endGame->condition_number, endGame->function_residual_d, endGame->latest_newton_residual_d, endGame->t_val_at_latest_sample_point_d, endGame->error_at_latest_sample_point_d, OUT, endGame->retVal, T);
    }
    else
    { // print the footer in multi precision and find retVal
      retVal = printMembershipFooter_mp(&endGame->PD_mp, endGame->condition_number, endGame->first_increase, endGame->function_residual_mp, endGame->latest_newton_residual_mp, endGame->t_val_at_latest_sample_point_mp, endGame->error_at_latest_sample_point_mp, OUT, endGame->retVal, T);
    }

    // finishing setup endPt_amp & endGame correctly
    endPt_amp->curr_prec = endGame->prec;
    if (endGame->prec < 64)
    { // copy _d
      point_cp_d(endPt_amp->endPt_d, endGame->PD_d.point);
      set_d(endPt_amp->finalT_d, endGame->PD_d.time);
    }
    else
    { // copy _mp
      setprec_point_mp(endPt_amp->endPt_mp, endPt_amp->curr_prec);
      setprec_mp(endPt_amp->finalT_mp, endPt_amp->curr_prec);
      point_cp_mp(endPt_amp->endPt_mp, endGame->PD_mp.point);
      set_mp(endPt_amp->finalT_mp, endGame->PD_mp.time);
    }
    endPt_amp->last_approx_prec = endGame->last_approx_prec;
    if (endGame->last_approx_prec < 64)
    { // copy _d
      point_cp_d(endPt_amp->last_approx_d, endGame->last_approx_d);
    }
    else
    { // copy _mp
      setprec_point_mp(endPt_amp->last_approx_mp, endPt_amp->last_approx_prec);
      point_cp_mp(endPt_amp->last_approx_mp, endGame->last_approx_mp);
    }
    endPt_amp->cond_num = endGame->condition_number;
    endPt_amp->retVal = endGame->retVal = retVal;

    clear_point_data_d(&startPt);
  }

  return retVal;
}


