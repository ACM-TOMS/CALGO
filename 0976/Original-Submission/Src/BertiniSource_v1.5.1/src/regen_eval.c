// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"
#include "regeneration.h"

int regen_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for regeneration          * 
\***************************************************************/
{
  regen_t *regen = (regen_t *) ED; // to avoid casting the void const * every time
  int i, j, k, l, m, level_num = regen->curr_level_num, num_funcs = regen->num_funcs, num_vars = vars->size;
  int level = regen->level[level_num].level, // the starting level for the regeneration
      depth = regen->level[level_num].depth, // the depth of this level of the regeneration
      useIntrinsic = regen->level[level_num].useIntrinsicSlice, // whether we are using intrinsic slicing for this level of the regeneration
      num_var_gps = regen->num_var_gps,      // total number of variable groups
      num_variables = regen->num_variables;  // total number of variables
  int num_orig_vars = num_variables - 1; // total number of orig variables
  int end = level + depth;

  double tempDouble;
  mat_d Jv_F;
  vec_d F;
  comp_d one_minus_t, gamma_t, prod_gamma, tempComp, sum, **linears = (comp_d **)bmalloc(num_var_gps * sizeof(comp_d *));
  // initialize
  init_mat_d(Jv_F, 0, 0);
  init_vec_d(F, 0); 
  for (i = 0; i < num_var_gps; i++)
    linears[i] = NULL;

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_d(gamma_t, regen->gamma_d, pathVars);

  if (useIntrinsic)   
  { // using intrinsic formulation 
    vec_d orig_vars;
    init_vec_d(orig_vars, 0);

    // convert vars to the original variables (p + B*vars)
    intrinsicToExtrinsic_d(orig_vars, vars, regen->level[level_num].B_d, regen->level[level_num].p_d); 

    // remove the extra hom variable
    orig_vars->size--;

    // evaluate the 'square' system
    regen_square_eval_d(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, regen);

    // set the size back
    orig_vars->size++;

    // calculate funcVals = [F, F * (1-t) + gamma * t * prod_of_linears]
    // calculate Jp = [0, -F + gamma * prod_of_linears]
    // calculate Jv = [Jv_F, Jv_F * (1-t) + gamma * t * d(prod_of_linears)]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, num_variables);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jp->rows = Jv->rows = num_vars;
    Jv->cols = orig_vars->size; // convert back to input vars after it is setup
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < level)
      { // funcVals
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_d(&Jv->entry[i][num_orig_vars]);
      }
      else
      { // adjust the values based on level & depth
        // to do this, we need to evaluate the approriate linears - go through all the linears needed for this function

        // initialize prod_gamma = gamma
        set_d(prod_gamma, regen->gamma_d);
        // loop over each variable group
        for (j = 0; j < num_var_gps; j++)
        { // allocate linears[j]
          linears[j] = (comp_d *)brealloc(linears[j], regen->degrees[i][j] * sizeof(comp_d));
          // loop over each degree associated to function i & variable group j
          for (k = 0; k < regen->degrees[i][j]; k++) 
          { // setup linears[j][k]
            set_zero_d(linears[j][k]);
            for (l = 0; l < num_variables; l++) // for each original variable
            { // multiply variable by coefficient and add on
              sum_mul_d2(linears[j][k], regen->coeff_d[i][j][k][l], &orig_vars->coord[l], tempDouble);
            }
            // gamma_prod *= linears
            mul_d(prod_gamma, prod_gamma, linears[j][k]); 
          }
        }

        // funcVals
        mul_d(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
        sum_mul_d2(&funcVals->coord[i], &F->coord[i], one_minus_t, tempDouble); // += F * (1-t)
        // Jp
        sub_d(&Jp->entry[i][0], prod_gamma, &F->coord[i]); // -F + gamma * prod_of_linears
        // Jv
        for (j = 0; j < num_variables; j++)
        { // find d(prod_of_linears, w.r.t. x_j)
          set_zero_d(&Jv->entry[i][j]);
          // loop over the variable groups
          for (k = 0; k < num_var_gps; k++) 
          {
            set_zero_d(sum);
            for (l = 0; l < regen->degrees[i][k]; l++)
            {
              set_d(tempComp, regen->coeff_d[i][k][l][j]);
              for (m = 0; m < regen->degrees[i][k]; m++)
                if (m != l)
                {
                  mul_d(tempComp, tempComp, linears[k][m]);
                }
              add_d(sum, sum, tempComp);
            }

            for (l = 0; l < num_var_gps; l++) // loop over variable groups not k
              if (l != k)
                for (m = 0; m < regen->degrees[i][l]; m++)
                {
                  mul_d(sum, sum, linears[l][m]);
                }
            add_d(&Jv->entry[i][j], &Jv->entry[i][j], sum);
          }
          // multiply by gamma*t
          mul_d(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);

          if (j < num_orig_vars)
          { // add to this Jv_F*(1-t)
            sum_mul_d2(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t, tempDouble);
          }
        }
      }

    // convert Jv to the input vars
    mat_mul_d(Jv, Jv, regen->level[level_num].B_d);

    // clear memory
    clear_vec_d(orig_vars);
  }
  else
  { // the input vars are the original vars - evaluate normally

    // remove the extra hom variable
    vars->size--;

    // evaluate the 'square' system
    regen_square_eval_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, regen);

    // set the size back
    vars->size++;

    // calculate funcVals = [F, F * (1-t) + gamma * t * prod_of_linears, extra_linears, patchValues]
    // calculate Jp = [0, -F + gamma * prod_of_linears, 0, 0]
    // calculate Jv = [Jv_F, Jv_F * (1-t) + gamma * t * d(prod_of_linears), d(extra_linears), Jv_Patch)]
    change_size_vec_d(funcVals, num_funcs + num_var_gps + 1);
    change_size_mat_d(Jv, num_funcs + num_var_gps + 1, num_vars);
    change_size_mat_d(Jp, num_funcs + num_var_gps + 1, 1);
    funcVals->size = Jp->rows = Jv->rows = num_funcs + num_var_gps + 1;
    Jv->cols = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_variables; i++)
      if (i < level)
      { // funcVals
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_d(&Jv->entry[i][num_orig_vars]); 
      }
      else if (i < end) // level <= i < end
      { // initialize prod_gamma = gamma
        set_d(prod_gamma, regen->gamma_d);
        // loop over each variable group
        for (j = 0; j < num_var_gps; j++) 
        { // allocate linears[j]
          linears[j] = (comp_d *)brealloc(linears[j], regen->degrees[i][j] * sizeof(comp_d));
          for (k = 0; k < regen->degrees[i][j]; k++) // for each degree associated to function i and variable group j
          { // setup linears[j][k]
            set_zero_d(linears[j][k]);
            for (l = 0; l < num_vars; l++) // for each variable
            { // multiply variable by coefficient and add on
              sum_mul_d2(linears[j][k], regen->coeff_d[i][j][k][l], &vars->coord[l], tempDouble);
            }
            // prod_gamma *= linears
            mul_d(prod_gamma, prod_gamma, linears[j][k]);
          }
        }

        // funcVals
        mul_d(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
        sum_mul_d2(&funcVals->coord[i], &F->coord[i], one_minus_t, tempDouble); // += F * (1-t)
        // Jp
        sub_d(&Jp->entry[i][0], prod_gamma, &F->coord[i]); // -F + gamma * prod_of_linears
        // Jv
        for (j = 0; j < num_vars; j++)
        { // find d(prod_of_linears, w.r.t. x_j)
          set_zero_d(&Jv->entry[i][j]);
          for (k = 0; k < num_var_gps; k++) // loop over the variable groups
          {
            set_zero_d(sum);
            for (l = 0; l < regen->degrees[i][k]; l++)
            {
              set_d(tempComp, regen->coeff_d[i][k][l][j]);
              for (m = 0; m < regen->degrees[i][k]; m++)
                if (m != l)
                {
                  mul_d(tempComp, tempComp, linears[k][m]);
                }
              add_d(sum, sum, tempComp);
            }

            for (l = 0; l < num_var_gps; l++) // loop over variable groups not k
              if (l != k)
                for (m = 0; m < regen->degrees[i][l]; m++)
                {
                  mul_d(sum, sum, linears[l][m]);
                }
            add_d(&Jv->entry[i][j], &Jv->entry[i][j], sum);
          }
          // multiply by gamma*t
          mul_d(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);
          if (j < num_orig_vars)
          { // add to this Jv_F*(1-t)
            sum_mul_d2(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t, tempDouble);
          }
        }
      }
      else if (i < num_funcs) // end <= i < num_funcs
      { // main slices

        // funcVals & Jv
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // mutlipy variable by coefficient and add on
          sum_mul_d2(&funcVals->coord[i], &regen->mainCoeff_d->entry[i][j], &vars->coord[j], tempDouble);
          set_d(&Jv->entry[i][j], &regen->mainCoeff_d->entry[i][j]);
        }

        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else if (i < num_funcs + num_var_gps) // num_funcs <= i < num_funcs + num_var_gps
      { // evaluate the patch
  
        k = i - num_funcs;
        // funcVals & Jv
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patchValues
          sum_mul_d2(&funcVals->coord[i], &regen->patchCoeff_d->entry[k][j], &vars->coord[j], tempDouble);
          // setup Jv_Patch
          set_d(&Jv->entry[i][j], &regen->patchCoeff_d->entry[k][j]);
        }

        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else
      { // evaluate the hom var patch
        // funcVals & Jv
        set_neg_one_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patch
          sum_mul_d2(&funcVals->coord[i], &regen->main_homVar_d->coord[j], &vars->coord[j], tempDouble);
          // setup Jv
          set_d(&Jv->entry[i][j], &regen->main_homVar_d->coord[j]);
        }

        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
  }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  // release the memory
  clear_mat_d(Jv_F);
  clear_vec_d(F); 
  for (i = num_var_gps - 1; i >= 0; i--)
    free(linears[i]);
  free(linears);

  return 0;
}

int regen_moving_linear_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluator to move linears for regeneration             *
\***************************************************************/
{
  regen_t *regen = (regen_t *) ED; // to avoid casting the void const * every time
  int i, j, var_gp, next_deg, level_num = regen->curr_level_num, num_funcs = regen->num_funcs, num_vars = vars->size;
  int level = regen->level[level_num].level, // the starting level for the regeneration
      depth = regen->level[level_num].depth, // the depth of this level of the regeneration
      useIntrinsic = regen->level[level_num].useIntrinsicSlice, // whether we are using intrinsic slicing for this level of the regeneration
      num_var_gps = regen->num_var_gps,      // total number of variable groups
      num_variables = regen->num_variables;  // total number of variables
  int num_orig_vars = num_variables - 1; // total number of orig variables
  int end = level + depth;

  mat_d Jv_F;
  vec_d F;
  comp_d one_minus_t, gamma_t, curr_linear, next_linear;
  double tempDouble;

  init_mat_d(Jv_F, 0, 0);
  init_vec_d(F, 0); 

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_d(gamma_t, regen->gamma_d, pathVars);

  if (useIntrinsic)
  { // use intrinsic formulation
    vec_d orig_vars;
    init_vec_d(orig_vars, 0);

    // convert vars to the original vars (p + B*vars)
    intrinsicToExtrinsic_d(orig_vars, vars, regen->level[level_num].B_d, regen->level[level_num].p_d);

    // remove the extra hom variable
    orig_vars->size--;

    // evaluate the 'square' system
    regen_square_eval_d(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, regen);

    // set the size back
    orig_vars->size++;

    // calculate funcVals = [F, (1-t) * next_linears + gamma * t * curr_linears]
    // calculate Jp = [0, -next_linears + gamma * curr_linears]
    // calculate Jv = [Jv_F, (1-t) * d(next_linears) + gamma * t * d(curr_linears)]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, num_variables);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jp->rows = Jv->rows = num_vars;
    Jv->cols = orig_vars->size; // convert back to input vars after it is setup
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < level)
      { // funcVals
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_d(&Jv->entry[i][num_orig_vars]);
      }
      else
      { // find out what variable group it is currently on and where it is moving to
        var_gp = regen->curr_linear[i];
        next_deg = regen->curr_linear_degree[i]; // where it is moving to

        // calculate curr_linear & next_linear and setup Jv
        set_zero_d(curr_linear);
        set_zero_d(next_linear);
        for (j = 0; j < num_variables; j++)
        { // multiply variable by coefficient and add on
          sum_mul_d2(curr_linear, &regen->mainCoeff_d->entry[i][j], &orig_vars->coord[j], tempDouble);
          sum_mul_d2(next_linear, regen->coeff_d[i][var_gp][next_deg][j], &orig_vars->coord[j], tempDouble);

          // setup Jv[i][j]
          mul_d(&Jv->entry[i][j], gamma_t, &regen->mainCoeff_d->entry[i][j]); // gamma * t * d(curr_linear)
          sum_mul_d2(&Jv->entry[i][j], one_minus_t, regen->coeff_d[i][var_gp][next_deg][j], tempDouble); // += (1-t) * d(next_linear)
        }
        // funcVals
        mul_d(&funcVals->coord[i], gamma_t, curr_linear); // gamma * t * curr_linear
        sum_mul_d2(&funcVals->coord[i], one_minus_t, next_linear, tempDouble); // += (1-t) * next_linear
        // Jp
        mul_d(&Jp->entry[i][0], regen->gamma_d, curr_linear); // gamma * curr_linear
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], next_linear); // -= next_linear
      }

    // convert Jv to the input vars
    mat_mul_d(Jv, Jv, regen->level[level_num].B_d);

    // clear memory
    clear_vec_d(orig_vars);
  }
  else
  { // the input vars are the original vars - evaluate normally

    // remove the extra hom variable
    vars->size--;

    // evaluate the 'square' system
    regen_square_eval_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, regen);

    // set the size back
    vars->size++;

    // calculate funcVals = [F, (1-t) * next_linears + gamma * t * curr_linears, extra_linears, patchValues]
    // calculate Jp = [0, -next_linears + gamma * curr_linears, 0, 0]
    // calculate Jv = [Jv_F, (1-t) * d(next_linears) + gamma * t * d(curr_linears), d(extra_linears), Jv_Patch]
    change_size_vec_d(funcVals, num_funcs + num_var_gps + 1);
    change_size_mat_d(Jv, num_funcs + num_var_gps + 1, num_vars);
    change_size_mat_d(Jp, num_funcs + num_var_gps + 1, 1);
    funcVals->size = Jp->rows = Jv->rows = num_funcs + num_var_gps + 1;
    Jv->cols = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_variables; i++)
      if (i < level)
      { // funcVals
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_d(&Jv->entry[i][num_orig_vars]);
      }
      else if (i < end) // level <= i < end
      { // find out what variable group it is currently on and where it is moving to
        var_gp = regen->curr_linear[i];
        next_deg = regen->curr_linear_degree[i]; // where it is moving to

        // calculate curr_linear & next_linear and setup Jv
        set_zero_d(curr_linear);
        set_zero_d(next_linear);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_d2(curr_linear, &regen->mainCoeff_d->entry[i][j], &vars->coord[j], tempDouble);
          sum_mul_d2(next_linear, regen->coeff_d[i][var_gp][next_deg][j], &vars->coord[j], tempDouble);

          // setup Jv[i][j]
          mul_d(&Jv->entry[i][j], gamma_t, &regen->mainCoeff_d->entry[i][j]); // gamma * t * d(curr_linear)
          sum_mul_d2(&Jv->entry[i][j], one_minus_t, regen->coeff_d[i][var_gp][next_deg][j], tempDouble); // += (1-t) * d(next_linear)
        }
        // funcVals
        mul_d(&funcVals->coord[i], gamma_t, curr_linear); // gamma * t * curr_linear
        sum_mul_d2(&funcVals->coord[i], one_minus_t, next_linear, tempDouble); // += (1-t) * next_linear
        // Jp
        mul_d(&Jp->entry[i][0], regen->gamma_d, curr_linear); // gamma * curr_linear
        sub_d(&Jp->entry[i][0], &Jp->entry[i][0], next_linear); // -= next_linear
      }
      else if (i < num_funcs) // end <= i < num_funcs
      { // main slices

        // funcVals & Jv
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_d2(&funcVals->coord[i], &regen->mainCoeff_d->entry[i][j], &vars->coord[j], tempDouble);
          set_d(&Jv->entry[i][j], &regen->mainCoeff_d->entry[i][j]);
        }
        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else if (i < num_funcs + num_var_gps) // num_funcs <= i < num_funcs + num_var_gps
      { // evaluate the patch
 
        var_gp = i - num_funcs;
        // funcVals & Jv
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_d2(&funcVals->coord[i], &regen->patchCoeff_d->entry[var_gp][j], &vars->coord[j], tempDouble);
          set_d(&Jv->entry[i][j], &regen->patchCoeff_d->entry[var_gp][j]);
        }
        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else
      { // evaluate the hom var patch
        // funcVals & Jv
        set_neg_one_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patch
          sum_mul_d2(&funcVals->coord[i], &regen->main_homVar_d->coord[j], &vars->coord[j], tempDouble);
          // setup Jv
          set_d(&Jv->entry[i][j], &regen->main_homVar_d->coord[j]);
        }
        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
  }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  clear_mat_d(Jv_F);
  clear_vec_d(F);

  return 0;
}

int regen_square_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: square function evaluation for regeneration            * 
\***************************************************************/
{
  regen_t *regen = (regen_t *)ED;
  int level_num = regen->curr_level_num;
  int end = regen->level[level_num].level + regen->level[level_num].depth;

  if (regen->noChanges)
  { // evaluate the SLP efficiently - from function 0 to end
    square_system_eval_data_d *SSED = (square_system_eval_data_d *)regen->square_d;

    // evaluate functions - a total of 'end' functions (0 to 'end - 1')
    evalProg_eff_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, SSED->Prog, 0, end, regen->startSub, regen->endSub, regen->startFunc, regen->endFunc, regen->startJvsub, regen->endJvsub, regen->startJv, regen->endJv, regen->subFuncsBelow); 
  }
  else 
  { // evaluate the square system using square_d
    eval_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, regen->square_d, regen->square_eval_d);

    // reorder the functions from bottom to top
    int i, j, numFuncs = funcVals->size, numVars = Jv->cols, numParams = Jp->cols;
    int N = (numFuncs % 2) ? (numFuncs - 1) / 2 : numFuncs / 2;
    comp_d tempD;

    for (i = 0; i < N; i++)
    { // swap ith and (numFuncs - 1 - i)th entry/row
      set_d(tempD, &funcVals->coord[i]);
      set_d(&funcVals->coord[i], &funcVals->coord[numFuncs - 1 - i]);
      set_d(&funcVals->coord[numFuncs - 1 - i], tempD);

      for (j = 0; j < numVars; j++)
      {
        set_d(tempD, &Jv->entry[i][j]);
        set_d(&Jv->entry[i][j], &Jv->entry[numFuncs - 1 - i][j]);
        set_d(&Jv->entry[numFuncs - 1 - i][j], tempD);
      }

      for (j = 0; j < numParams; j++)
      {
        set_d(tempD, &Jp->entry[i][j]);
        set_d(&Jp->entry[i][j], &Jp->entry[numFuncs - 1 - i][j]);
        set_d(&Jp->entry[numFuncs - 1 - i][j], tempD);
      }
    }
  }

  return 0;
}

int regen_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: standard function evaluation for regeneration          * 
\***************************************************************/
{
  regen_t *regen = (regen_t *) ED; // to avoid casting the void const * every time
  int i, j, k, l, m, level_num = regen->curr_level_num, num_funcs = regen->num_funcs, num_vars = vars->size;
  int level = regen->level[level_num].level, // the starting level for the regeneration
      depth = regen->level[level_num].depth, // the depth of this level of the regeneration
      useIntrinsic = regen->level[level_num].useIntrinsicSlice, // whether we are using intrinsic slicing for this level of the regeneration
      num_var_gps = regen->num_var_gps,      // total number of variable groups
      num_variables = regen->num_variables;  // total number of variables
  int num_orig_vars = num_variables - 1; // total number of orig variables
  int end = level + depth;

  mat_mp Jv_F;
  vec_mp F;
  comp_mp one_minus_t, gamma_t, prod_gamma, tempComp, sum, **linears = (comp_mp **)bmalloc(num_var_gps * sizeof(comp_mp *));

  // initialize MP
  init_mp(one_minus_t); init_mp(gamma_t); init_mp(tempComp); init_mp(sum); init_mp(prod_gamma);
  init_mat_mp(Jv_F, 0, 0);
  init_vec_mp(F, 0);
  for (i = 0; i < num_var_gps; i++)
    linears[i] = NULL;

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_mp(gamma_t, regen->gamma_mp, pathVars);

  if (useIntrinsic)
  { // using intrinsic formulation
    point_mp orig_vars;
    init_point_mp(orig_vars, 0);

    // convert vars to the original variables (p + B*vars)
    intrinsicToExtrinsic_mp(orig_vars, vars, regen->level[level_num].B_mp, regen->level[level_num].p_mp);

    // remove the extra hom variable
    orig_vars->size--;

    // evaluate the 'square' system
    regen_square_eval_mp(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, regen);

    // set the size back
    orig_vars->size++;

    // calculate funcVals = [F, F * (1-t) + gamma * t * prod_of_linears]
    // calculate Jp = [0, -F + gamma * prod_of_linears]
    // calculate Jv = [Jv_F, Jv_F * (1-t) + gamma * t * d(prod_of_linears)]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, num_variables);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jp->rows = Jv->rows = num_vars;
    Jv->cols = num_variables; // convert back to input vars after it is setup
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < level)
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_mp(&Jv->entry[i][num_orig_vars]);
      }
      else
      { // adjust the values based on level & depth
        // to do this, we need to evaluate the approriate linears - go through all the linears needed for this function

        // initialize prod_gamma = gamma
        set_mp(prod_gamma, regen->gamma_mp);
        // loop over each variable group
        for (j = 0; j < num_var_gps; j++)
        { // allocate linears[j]
          linears[j] = (comp_mp *)brealloc(linears[j], regen->degrees[i][j] * sizeof(comp_mp));
          // loop over each degree associated to function i & variable group j
          for (k = 0; k < regen->degrees[i][j]; k++)
          { // setup linears[j][k]
            init_mp(linears[j][k]);
            set_zero_mp(linears[j][k]);
            for (l = 0; l < num_variables; l++) // for each original variable
            { // multiply variable by coefficient and add on
              sum_mul_mp(linears[j][k], regen->coeff_mp[i][j][k][l], &orig_vars->coord[l]);
            }
            // gamma_prod *= linears
            mul_mp(prod_gamma, prod_gamma, linears[j][k]);
          }
        }

        // funcVals
        mul_mp(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
        sum_mul_mp(&funcVals->coord[i], &F->coord[i], one_minus_t); // += F * (1-t)
        // Jp
        sub_mp(&Jp->entry[i][0], prod_gamma, &F->coord[i]); // -F + gamma * prod_of_linears
        // Jv
        for (j = 0; j < num_variables; j++)
        { // find d(prod_of_linears, w.r.t. x_j)
          set_zero_mp(&Jv->entry[i][j]);
          // loop over the variable groups
          for (k = 0; k < num_var_gps; k++)
          {
            set_zero_mp(sum);
            for (l = 0; l < regen->degrees[i][k]; l++)
            {
              set_mp(tempComp, regen->coeff_mp[i][k][l][j]);
              for (m = 0; m < regen->degrees[i][k]; m++)
                if (m != l)
                {
                  mul_mp(tempComp, tempComp, linears[k][m]);
                }
              add_mp(sum, sum, tempComp);
            }

            for (l = 0; l < num_var_gps; l++) // loop over variable groups not k
              if (l != k)
                for (m = 0; m < regen->degrees[i][l]; m++)
                {
                  mul_mp(sum, sum, linears[l][m]);
                }
            add_mp(&Jv->entry[i][j], &Jv->entry[i][j], sum);
          }
          // multiply by gamma*t
          mul_mp(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);
          if (j < num_orig_vars)
          { // add to this Jv_F*(1-t)
            sum_mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
          }
        }

        // clear MP linears
        for (j = 0; j < num_var_gps; j++)
          for (k = 0; k < regen->degrees[i][j]; k++)
            clear_mp(linears[j][k]);
      }

    // convert Jv to the input vars
    mat_mul_mp(Jv, Jv, regen->level[level_num].B_mp);

    // clear memory
    clear_vec_mp(orig_vars);
  }
  else
  { // the input vars are the original vars - evaluate normally

    // remove the extra hom variable
    vars->size--;

    // evaluate the 'square' system
    regen_square_eval_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, regen);

    // set the size back
    vars->size++;

    // calculate funcVals = [F, F * (1-t) + gamma * t * prod_of_linears, extra_linears, patchValues]
    // calculate Jp = [0, -F + gamma * prod_of_linears, 0, 0]
    // calculate Jv = [Jv_F, Jv_F * (1-t) + gamma * t * d(prod_of_linears), d(extra_linears), Jv_Patch)]
    change_size_vec_mp(funcVals, num_funcs + num_var_gps + 1);
    change_size_mat_mp(Jv, num_funcs + num_var_gps + 1, num_vars);
    change_size_mat_mp(Jp, num_funcs + num_var_gps + 1, 1);
    funcVals->size = Jp->rows = Jv->rows = num_funcs + num_var_gps + 1;
    Jv->cols = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_variables; i++)
      if (i < level)
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_mp(&Jv->entry[i][num_orig_vars]);
      }
      else if (i < end) // level <= i < end
      { // initialize prod_gamma = gamma
        set_mp(prod_gamma, regen->gamma_mp);
        // loop over each variable group
        for (j = 0; j < num_var_gps; j++)
        { // allocate linears[j]
          linears[j] = (comp_mp *)brealloc(linears[j], regen->degrees[i][j] * sizeof(comp_mp));
          for (k = 0; k < regen->degrees[i][j]; k++) // for each degree associated to function i and variable group j
          { // setup linears[j][k]
            init_mp(linears[j][k]);
            set_zero_mp(linears[j][k]);
            for (l = 0; l < num_vars; l++) // for each variable
            { // multiply variable by coefficient and add on
              sum_mul_mp(linears[j][k], regen->coeff_mp[i][j][k][l], &vars->coord[l]);
            }
            // prod_gamma *= linears
            mul_mp(prod_gamma, prod_gamma, linears[j][k]);
          }
        }

        // funcVals
        mul_mp(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
        sum_mul_mp(&funcVals->coord[i], &F->coord[i], one_minus_t); // += F * (1-t)
        // Jp
        sub_mp(&Jp->entry[i][0], prod_gamma, &F->coord[i]); // -F + gamma * prod_of_linears
        // Jv
        for (j = 0; j < num_vars; j++)
        { // find d(prod_of_linears, w.r.t. x_j)
          set_zero_mp(&Jv->entry[i][j]);
          for (k = 0; k < num_var_gps; k++) // loop over the variable groups
          {
            set_zero_mp(sum);
            for (l = 0; l < regen->degrees[i][k]; l++)
            {
              set_mp(tempComp, regen->coeff_mp[i][k][l][j]);
              for (m = 0; m < regen->degrees[i][k]; m++)
                if (m != l)
                {
                  mul_mp(tempComp, tempComp, linears[k][m]);
                }
              add_mp(sum, sum, tempComp);
            }

            for (l = 0; l < num_var_gps; l++) // loop over variable groups not k
              if (l != k)
                for (m = 0; m < regen->degrees[i][l]; m++)
                {
                  mul_mp(sum, sum, linears[l][m]);
                }
            add_mp(&Jv->entry[i][j], &Jv->entry[i][j], sum);
          }
          // multiply by gamma*t
          mul_mp(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);
          if (j < num_orig_vars)
          { // add to this Jv_F*(1-t)
            sum_mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
          }
        }

        // clear MP linears
        for (j = 0; j < num_var_gps; j++)
          for (k = 0; k < regen->degrees[i][j]; k++)
            clear_mp(linears[j][k]);
      }
      else if (i < num_funcs) // end <= i < num_funcs
      { // main slices

        // funcVals & Jv
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // mutlipy variable by coefficient and add on
          sum_mul_mp(&funcVals->coord[i], &regen->mainCoeff_mp->entry[i][j], &vars->coord[j]);
          set_mp(&Jv->entry[i][j], &regen->mainCoeff_mp->entry[i][j]);
        }

        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else if (i < num_funcs + num_var_gps) // num_funcs <= i < num_funcs + num_var_gps
      { // evaluate the patch

        k = i - num_funcs;
        // funcVals & Jv
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patchValues
          sum_mul_mp(&funcVals->coord[i], &regen->patchCoeff_mp->entry[k][j], &vars->coord[j]);
          // setup Jv_Patch
          set_mp(&Jv->entry[i][j], &regen->patchCoeff_mp->entry[k][j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else
      { // evaluate the hom var patch
        // funcVals & Jv
        set_neg_one_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patch
          sum_mul_mp(&funcVals->coord[i], &regen->main_homVar_mp->coord[j], &vars->coord[j]);
          // setup Jv
          set_mp(&Jv->entry[i][j], &regen->main_homVar_mp->coord[j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]); // ds/dt = 1

  // release the memory
  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(tempComp); clear_mp(sum); clear_mp(prod_gamma);
  clear_mat_mp(Jv_F);
  clear_vec_mp(F);
  for (i = num_var_gps - 1; i >= 0; i--)
    free(linears[i]);
  free(linears);

  return 0;
}

int regen_moving_linear_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: evaluator to move linears for regeneration             *
\***************************************************************/
{
  regen_t *regen = (regen_t *) ED; // to avoid casting the void const * every time
  int i, j, var_gp, next_deg, level_num = regen->curr_level_num, num_funcs = regen->num_funcs, num_vars = vars->size;
  int level = regen->level[level_num].level, // the starting level for the regeneration
      depth = regen->level[level_num].depth, // the depth of this level of the regeneration
      useIntrinsic = regen->level[level_num].useIntrinsicSlice, // whether we are using intrinsic slicing for this level of the regeneration
      num_var_gps = regen->num_var_gps,      // total number of variable groups
      num_variables = regen->num_variables;  // total number of variables
  int num_orig_vars = num_variables - 1; // total number of orig variables
  int end = level + depth;

  mat_mp Jv_F;
  vec_mp F;
  comp_mp one_minus_t, gamma_t, curr_linear, next_linear;

  init_mp(one_minus_t); init_mp(gamma_t); init_mp(curr_linear); init_mp(next_linear);
  init_mat_mp(Jv_F, 0, 0);
  init_vec_mp(F, 0);

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma*t
  mul_mp(gamma_t, regen->gamma_mp, pathVars);

  if (useIntrinsic)
  { // use intrinsic formulation
    vec_mp orig_vars;
    init_vec_mp(orig_vars, 0);

    // convert vars to the original vars (p + B*vars)
    intrinsicToExtrinsic_mp(orig_vars, vars, regen->level[level_num].B_mp, regen->level[level_num].p_mp);

    // remove the extra hom variable
    orig_vars->size--;

    // evaluate the 'square' system
    regen_square_eval_mp(F, parVals, parDer, Jv_F, Jp, orig_vars, pathVars, regen);

    // set the size back
    orig_vars->size++;

    // calculate funcVals = [F, (1-t) * next_linears + gamma * t * curr_linears]
    // calculate Jp = [0, -next_linears + gamma * curr_linears]
    // calculate Jv = [Jv_F, (1-t) * d(next_linears) + gamma * t * d(curr_linears)]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, orig_vars->size);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jp->rows = Jv->rows = num_vars;
    Jv->cols = orig_vars->size; // convert back to input vars after it is setup
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < level)
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_mp(&Jv->entry[i][num_orig_vars]);
      }
      else
      { // find out what variable group it is currently on and where it is moving to
        var_gp = regen->curr_linear[i];
        next_deg = regen->curr_linear_degree[i]; // where it is moving to

        // calculate curr_linear & next_linear and setup Jv
        set_zero_mp(curr_linear);
        set_zero_mp(next_linear);
        for (j = 0; j < num_variables; j++)
        { // multiply variable by coefficient and add on
          sum_mul_mp(curr_linear, &regen->mainCoeff_mp->entry[i][j], &orig_vars->coord[j]);
          sum_mul_mp(next_linear, regen->coeff_mp[i][var_gp][next_deg][j], &orig_vars->coord[j]);

          // setup Jv[i][j]
          mul_mp(&Jv->entry[i][j], gamma_t, &regen->mainCoeff_mp->entry[i][j]); // gamma * t * d(curr_linear)
          sum_mul_mp(&Jv->entry[i][j], one_minus_t, regen->coeff_mp[i][var_gp][next_deg][j]); // += (1-t) * d(next_linear)
        }
        // funcVals
        mul_mp(&funcVals->coord[i], gamma_t, curr_linear); // gamma * t * curr_linear
        sum_mul_mp(&funcVals->coord[i], one_minus_t, next_linear); // += (1-t) * next_linear
        // Jp
        mul_mp(&Jp->entry[i][0], regen->gamma_mp, curr_linear); // gamma * curr_linear
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], next_linear); // -= next_linear
      }

    // convert Jv to the input vars
    mat_mul_mp(Jv, Jv, regen->level[level_num].B_mp);

    // clear memory
    clear_vec_mp(orig_vars);
  }
  else
  { // the input vars are the original vars - evaluate normally

    // remove the extra hom variable
    vars->size--;

    // evaluate the 'square' system
    regen_square_eval_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, regen);

    // set the size back
    vars->size++;

    // calculate funcVals = [F, (1-t) * next_linears + gamma * t * curr_linears, extra_linears, patchValues]
    // calculate Jp = [0, -next_linears + gamma * curr_linears, 0, 0]
    // calculate Jv = [Jv_F, (1-t) * d(next_linears) + gamma * t * d(curr_linears), d(extra_linears), Jv_Patch]
    change_size_vec_mp(funcVals, num_funcs + num_var_gps + 1);
    change_size_mat_mp(Jv, num_funcs + num_var_gps + 1, num_vars);
    change_size_mat_mp(Jp, num_funcs + num_var_gps + 1, 1);
    funcVals->size = Jp->rows = Jv->rows = num_funcs + num_var_gps + 1;
    Jv->cols = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_variables; i++)
      if (i < level)
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_orig_vars; j++)
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        set_zero_mp(&Jv->entry[i][num_orig_vars]);
      }
      else if (i < end) // level <= i < end
      { // find out what variable group it is currently on and where it is moving to
        var_gp = regen->curr_linear[i];
        next_deg = regen->curr_linear_degree[i]; // where it is moving to

        // calculate curr_linear & next_linear and setup Jv
        set_zero_mp(curr_linear);
        set_zero_mp(next_linear);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_mp(curr_linear, &regen->mainCoeff_mp->entry[i][j], &vars->coord[j]);
          sum_mul_mp(next_linear, regen->coeff_mp[i][var_gp][next_deg][j], &vars->coord[j]);

          // setup Jv[i][j]
          mul_mp(&Jv->entry[i][j], gamma_t, &regen->mainCoeff_mp->entry[i][j]); // gamma * t * d(curr_linear)
          sum_mul_mp(&Jv->entry[i][j], one_minus_t, regen->coeff_mp[i][var_gp][next_deg][j]); // += (1-t) * d(next_linear)
        }
        // funcVals
        mul_mp(&funcVals->coord[i], gamma_t, curr_linear); // gamma * t * curr_linear
        sum_mul_mp(&funcVals->coord[i], one_minus_t, next_linear); // += (1-t) * next_linear
        // Jp
        mul_mp(&Jp->entry[i][0], regen->gamma_mp, curr_linear); // gamma * curr_linear
        sub_mp(&Jp->entry[i][0], &Jp->entry[i][0], next_linear); // -= next_linear
      }
      else if (i < num_funcs) // end <= i < num_funcs
      { // main slices

        // funcVals & Jv
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_mp(&funcVals->coord[i], &regen->mainCoeff_mp->entry[i][j], &vars->coord[j]);
          set_mp(&Jv->entry[i][j], &regen->mainCoeff_mp->entry[i][j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else if (i < num_funcs + num_var_gps) // num_funcs <= i < num_funcs + num_var_gps
      { // evaluate the patch

        var_gp = i - num_funcs;
        // funcVals & Jv
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // multiply variable by coefficient and add on
          sum_mul_mp(&funcVals->coord[i], &regen->patchCoeff_mp->entry[var_gp][j], &vars->coord[j]);
          set_mp(&Jv->entry[i][j], &regen->patchCoeff_mp->entry[var_gp][j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else
      { // evaluate the hom var patch
        // funcVals & Jv
        set_neg_one_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update patch
          sum_mul_mp(&funcVals->coord[i], &regen->main_homVar_mp->coord[j], &vars->coord[j]);
          // setup Jv
          set_mp(&Jv->entry[i][j], &regen->main_homVar_mp->coord[j]); 
        } 
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]); // ds/dt = 1

  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(curr_linear); clear_mp(next_linear);
  clear_mat_mp(Jv_F);
  clear_vec_mp(F);

  return 0;
}

int regen_square_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: square function evaluation for regeneration            * 
\***************************************************************/
{
  regen_t *regen = (regen_t *)ED;
  int level_num = regen->curr_level_num;
  int end = regen->level[level_num].level + regen->level[level_num].depth;

  if (regen->noChanges)
  { // evaluate the SLP efficiently - from function 0 to end
    square_system_eval_data_mp *SSED = (square_system_eval_data_mp *)regen->square_mp;
  
    // evaluate functions - a total of 'end' functions (0 to 'end - 1')
    evalProg_eff_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, SSED->Prog, 0, end, regen->startSub, regen->endSub, regen->startFunc, regen->endFunc, regen->startJvsub, regen->endJvsub, regen->startJv, regen->endJv, regen->subFuncsBelow);
  }
  else 
  { // evaluate the square system using square_mp
    eval_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, regen->square_mp, regen->square_eval_mp);

    // reorder the functions from bottom to top
    int i, j, numFuncs = funcVals->size, numVars = Jv->cols, numParams = Jp->cols;
    int N = (numFuncs % 2) ? (numFuncs - 1) / 2 : numFuncs / 2;
    comp_mp tempMP;

    init_mp(tempMP);

    for (i = 0; i < N; i++)
    { // swap ith and (numFuncs - 1 - i)th entry/row
      set_mp(tempMP, &funcVals->coord[i]);
      set_mp(&funcVals->coord[i], &funcVals->coord[numFuncs - 1 - i]);
      set_mp(&funcVals->coord[numFuncs - 1 - i], tempMP);

      for (j = 0; j < numVars; j++)
      {
        set_mp(tempMP, &Jv->entry[i][j]);
        set_mp(&Jv->entry[i][j], &Jv->entry[numFuncs - 1 - i][j]);
        set_mp(&Jv->entry[numFuncs - 1 - i][j], tempMP);
      }

      for (j = 0; j < numParams; j++)
      {
        set_mp(tempMP, &Jp->entry[i][j]);
        set_mp(&Jp->entry[i][j], &Jp->entry[numFuncs - 1 - i][j]);
        set_mp(&Jp->entry[numFuncs - 1 - i][j], tempMP);
      }
    }

    clear_mp(tempMP);
  }

  return 0;
}


