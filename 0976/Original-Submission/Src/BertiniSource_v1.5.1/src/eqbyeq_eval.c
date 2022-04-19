// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"
#include "eqbyeq.h"

int witness_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to find the original witness sets  *
\***************************************************************/
{
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;
  int i, j, k, l, deg, count, num_sub = BED->EqD->curr_stage_num, num_vars = vars->size;
  int startFunc = BED->EqD->witnessData_d[num_sub].startFunction,// starting function
      depth = BED->EqD->witnessData_d[num_sub].depth;            // depth for this subsystem

  comp_d one_minus_t, gamma_t, prod_gamma, tempComp, *curr_linears = (comp_d *)bmalloc(BED->EqD->witnessData_d[num_sub].num_linears * sizeof(comp_d));

  vec_d F, orig_vars;
  mat_d Jv_F;
  init_vec_d(F, 0); init_vec_d(orig_vars, 0);
  init_mat_d(Jv_F, 0, 0);

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_d(gamma_t, BED->EqD->gamma_d, pathVars);

  // convert vars to the original variables (p + B*vars)
  intrinsicToExtrinsic_d(orig_vars, vars, BED->EqD->witnessData_d[num_sub].B, BED->EqD->witnessData_d[num_sub].p);

  // evaluate the 'square' system
  eqbyeq_square_eval_d(F, parVals, parDer, Jv, Jp, orig_vars, pathVars, BED, startFunc, startFunc + depth);

  // pick out the rows of Jv that are needed and store to Jv_F
  change_size_mat_d(Jv_F, depth, Jv->cols);
  Jv_F->rows = depth;
  Jv_F->cols = Jv->cols;
  for (i = 0; i < depth; i++)
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_d(&Jv_F->entry[i][j], &Jv->entry[startFunc + i][j]);
    }
  // convert Jv_F to the input vars
  mat_mul_d(Jv_F, Jv_F, BED->EqD->witnessData_d[num_sub].B);

  // evaluate the linears
  count = BED->EqD->witnessData_d[num_sub].num_linears;
  for (i = 0; i < count; i++)
  { // intialize to -1
    set_one_d(curr_linears[i]);
    neg_d(curr_linears[i], curr_linears[i]);
    for (j = 0; j < num_vars; j++)
    { // multiply variable by coefficient and add on
      sum_mul_d(curr_linears[i], BED->EqD->witnessData_d[num_sub].startSystemCoeff[i][j], &vars->coord[j]);
    }
  }

  // calculate funcVals = [F * (1-t) + gamma * t * prod_of_linears]
  // calculate Jp = [-F + gamma * prod_of_linears]
  // calculate Jv = [Jv_F * (1-t) + gamma * t * d(prod_of_linears)]
  change_size_vec_d(funcVals, depth);
  change_size_mat_d(Jv, depth, depth);
  change_size_mat_d(Jp, depth, 1);
  funcVals->size = Jv->rows = Jv->cols = Jp->rows = depth;
  Jp->cols = 1;
  count = 0;
  for (i = 0; i < depth; i++)
  { // initialize prod_gamma = gamma
    set_d(prod_gamma, BED->EqD->gamma_d);
    // loop over the linears
    deg = BED->EqD->degrees[startFunc + i][0];
    for (j = 0; j < deg; j++)
    { // multiply by linear
      mul_d(prod_gamma, prod_gamma, curr_linears[count]);
      count++;
    }

    // funcVals
    mul_d(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
    sum_mul_d(&funcVals->coord[i], &F->coord[startFunc + i], one_minus_t); // += F * (1-t)
    // Jp
    sub_d(&Jp->entry[i][0], prod_gamma, &F->coord[startFunc + i]); // -F + gamma * prod_of_linears
    // Jv
    count -= deg;
    for (j = 0; j < num_vars; j++)
    { // find d(prod_of_linears w.r.t. x_j)
      set_zero_d(&Jv->entry[i][j]);
      for (k = 0; k < deg; k++)
      {
        set_d(tempComp, BED->EqD->witnessData_d[num_sub].startSystemCoeff[count+k][j]);
        for (l = 0; l < deg; l++)
          if (l != k)
          {
            mul_d(tempComp, tempComp, curr_linears[count+l]);
          }
        add_d(&Jv->entry[i][j], &Jv->entry[i][j], tempComp);
      }
      // multiply by gamma_t
      mul_d(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);
      // add on Jv_F * (1-t)
      sum_mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
    }
    // update count
    count += deg;
  }

  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1
 
  // release the memory
  clear_vec_d(F); clear_vec_d(orig_vars);
  clear_mat_d(Jv_F);
  free(curr_linears);

  return 0;
}

int standard_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to find the next stage             *
\***************************************************************/
{
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;
  int i, j, k, offset, stage = BED->EqD->curr_stage_num, num_vars = BED->EqD->num_vars, num_funcs = BED->EqD->num_funcs;
  int depth_x = BED->EqD->stageData_d[stage].depth_x,  // depth for the x-variables
      depth_y = BED->EqD->stageData_d[stage].depth_y,  // depth for the y-variables
      useIntrinsic = BED->EqD->stageData_d[stage].useIntrinsicSlice; // whether we are using intrinsic slicing for this stage
  int depth_sum = depth_x + depth_y;

  comp_d one_minus_t, gamma_t;
  point_d orig_vars_x, orig_vars_y;
  vec_d F;
  mat_d Jv_F, Jp_temp;

  init_point_d(orig_vars_x, 0); init_point_d(orig_vars_y, 0);
  init_vec_d(F, 0); 
  init_mat_d(Jv_F, 0, 0); init_mat_d(Jp_temp, 0, 0);

  // find 1 - t
  set_one_d(one_minus_t);
  sub_d(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_d(gamma_t, BED->EqD->gamma_d, pathVars);

  if (useIntrinsic)
  { // setup B = B0 + t*(gamma*B1 - B0) & p = p0 + t*(gamma*p1 - p0), Jp_temp = (gamma*p1 - p0) + (gamma*B1 - B0)*vars
    vec_d p;
    mat_d B;

    init_vec_d(p, BED->EqD->stageData_d[stage].p1->size);
    init_mat_d(B, BED->EqD->stageData_d[stage].B1->rows, BED->EqD->stageData_d[stage].B1->cols);
    change_size_mat_d(Jp_temp, BED->EqD->stageData_d[stage].p1->size, 1);
    
    B->rows = BED->EqD->stageData_d[stage].B1->rows;
    B->cols = BED->EqD->stageData_d[stage].B1->cols; // == vars->size
    Jp_temp->rows = p->size = BED->EqD->stageData_d[stage].p1->size; // == B->rows
    Jp_temp->cols = 1;

    for (i = 0; i < B->rows; i++)
    {
      set_zero_d(&Jp_temp->entry[i][0]);
      for (j = 0; j < B->cols; j++)
      { // find (gamma*B1 - B0)[i][j]
        mul_d(&B->entry[i][j], BED->EqD->gamma_d, &BED->EqD->stageData_d[stage].B1->entry[i][j]);
        sub_d(&B->entry[i][j], &B->entry[i][j], &BED->EqD->stageData_d[stage].B0->entry[i][j]);

        // Jp[i] += (gamma*B1 - B0)[i][j] * vars[j]
        sum_mul_d(&Jp_temp->entry[i][0], &B->entry[i][j], &vars->coord[j]);

        // setup B[i][j]
        mul_d(&B->entry[i][j], &B->entry[i][j], pathVars);
        add_d(&B->entry[i][j], &B->entry[i][j], &BED->EqD->stageData_d[stage].B0->entry[i][j]);
      }
      // find (gamma*p1 - p0)[i]
      mul_d(&p->coord[i], BED->EqD->gamma_d, &BED->EqD->stageData_d[stage].p1->coord[i]);
      sub_d(&p->coord[i], &p->coord[i], &BED->EqD->stageData_d[stage].p0->coord[i]);

      // Jp[i] += (gamma*p1 - p0)[i]
      add_d(&Jp_temp->entry[i][0], &Jp_temp->entry[i][0], &p->coord[i]);

      // setup p[i]
      mul_d(&p->coord[i], &p->coord[i], pathVars);
      add_d(&p->coord[i], &p->coord[i], &BED->EqD->stageData_d[stage].p0->coord[i]);
    }

    // convert vars to [x,y]
    intrinsicToExtrinsic_d(orig_vars_x, vars, B, p);

    // seup x & y
    change_size_point_d(orig_vars_y, num_vars);
    orig_vars_x->size = orig_vars_y->size = num_vars;
    for (i = 0; i < num_vars; i++)
    {
      set_d(&orig_vars_y->coord[i], &orig_vars_x->coord[i + num_vars]);
    }

    // allocate space for funcVals & Jv_F
    change_size_vec_d(funcVals, depth_sum);
    change_size_mat_d(Jv_F, depth_sum, 2 * num_vars);
    funcVals->size = Jv_F->rows = depth_sum;
    Jv_F->cols = 2 * num_vars;

    // evaluate the function for the x variables
    eqbyeq_square_eval_d(F, parVals, parDer, Jv, Jp, orig_vars_x, pathVars, BED, 0, depth_x);

    // setup top of funcVals & Jv_F
    for (i = 0; i < depth_x; i++)
    { // setup funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);
      // setup Jv_F
      for (j = 0; j < num_vars; j++)
      {
        set_d(&Jv_F->entry[i][j], &Jv->entry[i][j]);
        set_zero_d(&Jv_F->entry[i][j + num_vars]); // only depends on x variables
      }
    }

    // evaluate the function for the y variables
    eqbyeq_square_eval_d(F, parVals, parDer, Jv, Jp, orig_vars_y, pathVars, BED, depth_x, depth_sum);

    // setup bottom of funcVals & Jv_F
    for (i = depth_x; i < depth_sum; i++)
    { // setup funcVals
      set_d(&funcVals->coord[i], &F->coord[i]);
      // setup Jv_F
      for (j = 0; j < num_vars; j++)
      {
        set_d(&Jv_F->entry[i][j + num_vars], &Jv->entry[i][j]);
        set_zero_d(&Jv_F->entry[i][j]); // only depends on y variables
      }
    }

    // we already have funcVals = [F_x,F_y]

    // use the structure of Jv_F to find Jp = Jv_F * Jp_temp and Jv = Jv_F * B
    change_size_mat_d(Jv, depth_sum, B->cols);
    change_size_mat_d(Jp, depth_sum, 1);
    Jv->rows = Jp->rows = depth_sum;
    Jv->cols = B->cols; // == vars->size
    Jp->cols = 1;
    for (i = 0; i < depth_sum; i++)
    { // find offset
      offset = (i < depth_x) ? 0 : num_vars;

      for (j = 0; j < B->cols; j++)
      { // find Jv[i][j]
        set_zero_d(&Jv->entry[i][j]);
        for (k = 0; k < num_vars; k++)
        {
          sum_mul_d(&Jv->entry[i][j], &Jv_F->entry[i][k + offset], &B->entry[k + offset][j]);
        }
      }
      // find Jp[i][0]
      set_zero_d(&Jp->entry[i][0]);
      for (k = 0; k < num_vars; k++)
      {
        sum_mul_d(&Jp->entry[i][0], &Jv_F->entry[i][k + offset], &Jp_temp->entry[k + offset][0]);
      }
    }

    // clear memory
    clear_vec_d(p);
    clear_mat_d(B);
  }
  else
  { // using extrinsic formulation
    comp_d tempComp;
    vec_d patchValues;
    mat_d Jv_Patch;

    init_vec_d(patchValues, 0);
    init_mat_d(Jv_Patch, 0, 0);

    // setup the sizes for funcVals, Jv & Jp
    change_size_vec_d(funcVals, 2 * num_vars);
    change_size_mat_d(Jv, 2 * num_vars, 2 * num_vars);
    change_size_mat_d(Jp, 2 * num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = 2 * num_vars;
    Jp->cols = 1;

    // split vars into x & y
    change_size_vec_d(orig_vars_x, num_vars);
    change_size_vec_d(orig_vars_y, num_vars);
    orig_vars_x->size = orig_vars_y->size = num_vars;
    for (i = 0; i < num_vars; i++)
    {
      set_d(&orig_vars_x->coord[i], &vars->coord[i]);
      set_d(&orig_vars_y->coord[i], &vars->coord[i+num_vars]);
    }

    // evaluate the function for the x variables
    eqbyeq_square_eval_d(F, parVals, parDer, Jv_F, Jp_temp, orig_vars_x, pathVars, BED, 0, depth_x);
    // evaluate the patch for the x variables
    patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp_temp, orig_vars_x, pathVars, &BED->patch);

    // setup the top of funcVals, Jv & Jp
    // calculate funcVals(top) = [F,L_x(top)*gamma*t + (1-t)(x-y),L_x(rest),patch_x]
    // calculate Jp(top) = [0,-(x-y) + gamma*L_x(top),0,0]
    // calculate Jv(top) = [Jv_F,Jv_linears(top)*gamma*t + (1-t)d(x-y),Jv_linears(rest),Jv_Patch]
    for (i = 0; i < num_vars; i++)
      if (i < depth_x) // functions already solved for
      { // funcVals
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
          // y-coordinates
          set_zero_d(&Jv->entry[i][j + num_vars]);
        }
      }
      else if (i < depth_sum) // linears moving
      { // initialize tempComp - to be linears
        set_zero_d(tempComp);
        // evaluate the ith linear in x-variables
        for (j = 0; j < num_vars; j++)
        { // add on to tempComp
          sum_mul_d(tempComp, BED->EqD->coeff_d[i][j], &orig_vars_x->coord[j]);
          // setup Jv - x-coordinates
          mul_d(&Jv->entry[i][j], BED->EqD->coeff_d[i][j], gamma_t); // coeff * gamma * t
          // setup Jv - y-coordinates
          set_zero_d(&Jv->entry[i][j + num_vars]);
        }
        // funcVals
        sub_d(&funcVals->coord[i], &orig_vars_x->coord[i], &orig_vars_y->coord[i]); // x_i - y_i
        mul_d(&funcVals->coord[i], &funcVals->coord[i], one_minus_t); // (x_i - y_i) * (1-t)
        sum_mul_d(&funcVals->coord[i], tempComp, gamma_t); // += linear * gamma * t
        // Jp
        sub_d(&Jp->entry[i][0], &orig_vars_y->coord[i], &orig_vars_x->coord[i]); // y_i - x_i
        sum_mul_d(&Jp->entry[i][0], tempComp, BED->EqD->gamma_d); // (y_i - x_i) + gamma * linear
        // adjust Jv for i
        add_d(&Jv->entry[i][i], &Jv->entry[i][i], one_minus_t);
        sub_d(&Jv->entry[i][i + num_vars], &Jv->entry[i][i + num_vars], one_minus_t);
      }
      else if (i < num_funcs) // linears not moving
      { // funcVals - evaluate linear for x-variables
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // add on to funcVals
          sum_mul_d(&funcVals->coord[i], BED->EqD->coeff_d[i][j], &orig_vars_x->coord[j]);
          // setup Jv - x-coordinates
          set_d(&Jv->entry[i][j], BED->EqD->coeff_d[i][j]);
          // setup Jv - y-coordinates
          set_zero_d(&Jv->entry[i][j + num_vars]);
        }
        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else // patch
      { // funcVals
        set_d(&funcVals->coord[i], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_d(&Jv->entry[i][j], &Jv_Patch->entry[i - num_funcs][j]);
          // y-coordinates
          set_zero_d(&Jv->entry[i][j + num_vars]);
        }
      }

    // evaluate the function for the y variables
    eqbyeq_square_eval_d(F, parVals, parDer, Jv_F, Jp_temp, orig_vars_y, pathVars, BED, depth_x, depth_sum);
    // evaluate the patch for the y variables
    patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp_temp, orig_vars_y, pathVars, &BED->patch);

    // setup the bottom of funcVals, Jv & Jp
    // calculate funcVals(bottom) = [F_y,L_y(bottom)*gamma*t + (1-t)(x-y),L_y(rest),patch_y]
    // calculate Jp(bottom) = [0,-(x-y) + gamma*L_y(bottom),0,0]
    // calculate Jv(bottom) = [Jv_F,Jv_linears(bottom)*gamma*t + (1-t)d(x-y),Jv_linears(rest),Jv_Patch]
    for (i = 0; i < num_vars; i++)
      if (i < depth_x) // linears moving
      { // initialize tempComp - to be linear
        set_zero_d(tempComp);
        // evaluate the ith linear in x-variables
        for (j = 0; j < num_vars; j++)
        { // add on to tempComp
          sum_mul_d(tempComp, BED->EqD->coeff_d[i][j], &orig_vars_y->coord[j]);
          // setup Jv - y-coordinates
          mul_d(&Jv->entry[i + num_vars][j + num_vars], BED->EqD->coeff_d[i][j], gamma_t); // coeff * gamma * t
          // setup Jv - x-coordinates
          set_zero_d(&Jv->entry[i + num_vars][j]);
        }
        // funcVals
        sub_d(&funcVals->coord[i + num_vars], &orig_vars_x->coord[i], &orig_vars_y->coord[i]); // x_i - y_i
        mul_d(&funcVals->coord[i + num_vars], &funcVals->coord[i + num_vars], one_minus_t); // (x_i - y_i) * (1-t)
        sum_mul_d(&funcVals->coord[i + num_vars], tempComp, gamma_t); // += linear * gamma * t
        // Jp
        sub_d(&Jp->entry[i + num_vars][0], &orig_vars_y->coord[i], &orig_vars_x->coord[i]); // y_i - x_i
        sum_mul_d(&Jp->entry[i + num_vars][0], tempComp, BED->EqD->gamma_d); // (y_i - x_i) + gamma * linear
        // adjust Jv for i
        add_d(&Jv->entry[i + num_vars][i], &Jv->entry[i + num_vars][i], one_minus_t);
        sub_d(&Jv->entry[i + num_vars][i + num_vars], &Jv->entry[i + num_vars][i + num_vars], one_minus_t);
      }
      else if (i < depth_sum) // functions already solved for
      { // funcVals
        set_d(&funcVals->coord[i + num_vars], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i + num_vars][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_zero_d(&Jv->entry[i + num_vars][j]);
          // y-coordinates
          set_d(&Jv->entry[i + num_vars][j + num_vars], &Jv_F->entry[i][j]);
        }
      }
      else if (i < num_funcs) // linears not moving
      { // funcVals
        set_zero_d(&funcVals->coord[i + num_vars]);
        for (j = 0; j < num_vars; j++)
        { // add on to funcVals
          sum_mul_d(&funcVals->coord[i + num_vars], BED->EqD->coeff_d[i][j], &orig_vars_y->coord[j]);
          // setup Jv - y-coordinates
          set_d(&Jv->entry[i + num_vars][j + num_vars], BED->EqD->coeff_d[i][j]);
          // setup Jv - x-coordinates
          set_zero_d(&Jv->entry[i + num_vars][j]);
        }
        // Jp
        set_zero_d(&Jp->entry[i + num_vars][0]);
      }
      else // patch
      { // funcVals
        set_d(&funcVals->coord[i + num_vars], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_d(&Jp->entry[i + num_vars][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_zero_d(&Jv->entry[i + num_vars][j]);
          // y-coordinates
          set_d(&Jv->entry[i + num_vars][j + num_vars], &Jv_Patch->entry[i - num_funcs][j]);
        }
      }

    // clear memory
    clear_vec_d(patchValues);
    clear_mat_d(Jv_Patch);
  }

  // setup parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  clear_point_d(orig_vars_x); clear_point_d(orig_vars_y);
  clear_vec_d(F);
  clear_mat_d(Jv_F); clear_mat_d(Jp_temp);

  return 0;
}

int stage_sort_eqbyeq_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to sort the stage                  *
*   we know x = y  and t = 0                                    *
\***************************************************************/
// since we need to do SVD, we use this evaluator so that the slices 'x - y' are not used
{
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;
  int i, j, stage = BED->EqD->curr_stage_num, num_vars = BED->EqD->num_vars, num_funcs = BED->EqD->num_funcs;
  int depth_x = BED->EqD->stageData_d[stage].depth_x,  // depth for the x-variables
      depth_y = BED->EqD->stageData_d[stage].depth_y,  // depth for the y-variables
      useIntrinsic = BED->EqD->stageData_d[stage].useIntrinsicSlice; // whether we are using intrinsic slicing for this stage
  int depth_sum = depth_x + depth_y;

  if (useIntrinsic)
  {
    point_d orig_vars;
    init_point_d(orig_vars, 0);

    // now convert to the original variables (p + B*vars)
    intrinsicToExtrinsic_d(orig_vars, vars, BED->EqD->stageData_d[stage].B, BED->EqD->stageData_d[stage].p);

    // evaluate the 'square' function
    eqbyeq_square_eval_d(funcVals, parVals, parDer, Jv, Jp, orig_vars, pathVars, BED, 0, depth_sum);

    // pick out the correct rows of funcVals and Jv
    funcVals->size = Jv->rows = depth_sum;

    // convert Jv to the input vars
    mat_mul_d(Jv, Jv, BED->EqD->stageData_d[stage].B);

    // setup Jp == 0
    change_size_mat_d(Jp, depth_sum, 1);
    Jp->rows = depth_sum;
    Jp->cols = 1;
    for (i = 0; i < depth_sum; i++)
    {
      set_zero_d(&Jp->entry[i][0]);
    }

    // clear memory
    clear_point_d(orig_vars);
  }
  else
  { // use extrinsic formulation
    vec_d F, patchValues;
    mat_d Jv_F, Jv_Patch;

    init_vec_d(F, 0); init_vec_d(patchValues, 0);
    init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_Patch, 0, 0);
  
    // evaluate the 'square' function
    eqbyeq_square_eval_d(F, parVals, parDer, Jv_F, Jp, vars, pathVars, BED, 0, depth_sum);
    // evaluate the patch 
    patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);

    // calculate funcVals = [F, L, patch]
    // calculate Jp = [0, 0, 0]
    // calculate Jv = [Jv_F, Jv_linears, Jv_Patch]
    change_size_vec_d(funcVals, num_vars);
    change_size_mat_d(Jv, num_vars, num_vars);
    change_size_mat_d(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < depth_sum) // F
      { // funcVals 
        set_d(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        }
      }
      else if (i < num_funcs) // linears
      { // funcVals - linears
        set_zero_d(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update funcVals
          sum_mul_d(&funcVals->coord[i], BED->EqD->coeff_d[i][j], &vars->coord[j]);
          // setup Jv
          set_d(&Jv->entry[i][j], BED->EqD->coeff_d[i][j]);
        }
        // Jp
        set_zero_d(&Jp->entry[i][0]);
      }
      else // patch
      { // funcVals
        set_d(&funcVals->coord[i], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_d(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_d(&Jv->entry[i][j], &Jv_Patch->entry[i - num_funcs][j]);
        }
      }

    // clear memory
    clear_vec_d(F); clear_vec_d(patchValues);
    clear_mat_d(Jv_F); clear_mat_d(Jv_Patch);
  }

  // setup parVals & parDer correctly
  change_size_vec_d(parVals, 1);
  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars);   // s = t
  set_one_d(&parDer->coord[0]); // ds/dt = 1

  return 0;
}

int eqbyeq_square_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED, int startFunc, int endFunc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: square function evaluation for equation-by-equation    * 
\***************************************************************/
{
  basic_eval_data_d *BED = (basic_eval_data_d *)ED;

  if (BED->EqD->noChanges)
  { // evaluate the SLP efficiently - from function startFunc to endFunc
    evalProg_eff_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, BED->squareSystem.Prog, startFunc, endFunc, BED->EqD->startSub, BED->EqD->endSub, BED->EqD->startFunc, BED->EqD->endFunc, BED->EqD->startJvsub, BED->EqD->endJvsub, BED->EqD->startJv, BED->EqD->endJv, BED->EqD->subFuncsBelow);
  }
  else
  { // evaluate the square system
    square_system_eval_d(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, &BED->squareSystem);
  }

  return 0;
}

int witness_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to find the original witness sets  *
\***************************************************************/
{
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;
  int i, j, k, l, deg, count, num_sub = BED->EqD->curr_stage_num, num_vars = vars->size;
  int startFunc = BED->EqD->witnessData_mp[num_sub].startFunction,// starting function
      depth = BED->EqD->witnessData_mp[num_sub].depth;            // depth for this subsystem

  comp_mp one_minus_t, gamma_t, prod_gamma, tempComp, *curr_linears = (comp_mp *)bmalloc(BED->EqD->witnessData_mp[num_sub].num_linears * sizeof(comp_mp));

  // initialize MP
  init_mp(one_minus_t); init_mp(gamma_t); init_mp(prod_gamma); init_mp(tempComp);
  for (i = 0; i < BED->EqD->witnessData_mp[num_sub].num_linears; i++)
    init_mp(curr_linears[i]);

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_mp(gamma_t, BED->EqD->gamma_mp, pathVars);

  vec_mp F, orig_vars;
  mat_mp Jv_F;
  init_vec_mp(F, 0); init_vec_mp(orig_vars, 0);
  init_mat_mp(Jv_F, 0, 0);

  // convert vars to the original variables (p + B*vars)
  intrinsicToExtrinsic_mp(orig_vars, vars, BED->EqD->witnessData_mp[num_sub].B, BED->EqD->witnessData_mp[num_sub].p);

  // evaluate the 'square' system
  eqbyeq_square_eval_mp(F, parVals, parDer, Jv, Jp, orig_vars, pathVars, BED, startFunc, startFunc + depth);

  // pick out the rows of Jv that are needed and store to Jv_F
  change_size_mat_mp(Jv_F, depth, Jv->cols);
  Jv_F->rows = depth;
  Jv_F->cols = Jv->cols;
  for (i = 0; i < depth; i++)
    for (j = 0; j < Jv_F->cols; j++)
    {
      set_mp(&Jv_F->entry[i][j], &Jv->entry[startFunc + i][j]);
    }
  // convert Jv_F to the input vars
  mat_mul_mp(Jv_F, Jv_F, BED->EqD->witnessData_mp[num_sub].B);

  // evaluate the linears
  count = BED->EqD->witnessData_mp[num_sub].num_linears;
  for (i = 0; i < count; i++)
  { // intialize to -1
    set_one_mp(curr_linears[i]);
    neg_mp(curr_linears[i], curr_linears[i]);
    for (j = 0; j < num_vars; j++)
    { // multiply variable by coefficient and add on
      sum_mul_mp(curr_linears[i], BED->EqD->witnessData_mp[num_sub].startSystemCoeff[i][j], &vars->coord[j]);
    }
  }

  // calculate funcVals = [F * (1-t) + gamma * t * prod_of_linears]
  // calculate Jp = [-F + gamma * prod_of_linears]
  // calculate Jv = [Jv_F * (1-t) + gamma * t * d(prod_of_linears)]
  change_size_vec_mp(funcVals, depth);
  change_size_mat_mp(Jv, depth, depth);
  change_size_mat_mp(Jp, depth, 1);
  funcVals->size = Jp->rows = Jv->rows = Jv->cols = depth;
  Jp->cols = 1;
  count = 0;
  for (i = 0; i < depth; i++)
  { // intialize prod_gamma = gamma
    set_mp(prod_gamma, BED->EqD->gamma_mp);
    // loop over the linears
    deg = BED->EqD->degrees[startFunc + i][0];
    for (j = 0; j < deg; j++)
    { // multiply by linear
      mul_mp(prod_gamma, prod_gamma, curr_linears[count]);
      count++;
    }

    // funcVals
    mul_mp(&funcVals->coord[i], prod_gamma, pathVars); // gamma * prod_of_linears * t
    sum_mul_mp(&funcVals->coord[i], &F->coord[startFunc + i], one_minus_t); // += F * (1-t)
    // Jp
    sub_mp(&Jp->entry[i][0], prod_gamma, &F->coord[startFunc + i]); // -F + gamma * prod_of_linears
    // Jv
    count -= deg;
    for (j = 0; j < num_vars; j++)
    { // find d(prod_of_linears w.r.t. x_j) 
      set_zero_mp(&Jv->entry[i][j]);
      for (k = 0; k < deg; k++)
      {
        set_mp(tempComp, BED->EqD->witnessData_mp[num_sub].startSystemCoeff[count+k][j]);
        for (l = 0; l < deg; l++)
          if (l != k)
          {
            mul_mp(tempComp, tempComp, curr_linears[count+l]);
          }
        add_mp(&Jv->entry[i][j], &Jv->entry[i][j], tempComp);
      }
      // multiply by gamma_t
      mul_mp(&Jv->entry[i][j], &Jv->entry[i][j], gamma_t);
      // add on Jv_F * (1-t)
      sum_mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_t);
    }
    // update count
    count += deg;
  }

  // set parVals & parDer correctly 
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]);      // ds/dt = 1

  // release the memory
  for (i = BED->EqD->witnessData_mp[num_sub].num_linears - 1; i >= 0; i--)
  {
    clear_mp(curr_linears[i]);
  }
  free(curr_linears);
  clear_mp(one_minus_t); clear_mp(gamma_t); clear_mp(prod_gamma); clear_mp(tempComp);
  clear_vec_mp(F); clear_vec_mp(orig_vars);
  clear_mat_mp(Jv_F);

  return 0;
}

int standard_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to find the next stage             *
\***************************************************************/
{
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;
  int i, j, k, offset, stage = BED->EqD->curr_stage_num, num_vars = BED->EqD->num_vars, num_funcs = BED->EqD->num_funcs;
  int depth_x = BED->EqD->stageData_mp[stage].depth_x,  // depth for the x-variables
      depth_y = BED->EqD->stageData_mp[stage].depth_y,  // depth for the y-variables
      useIntrinsic = BED->EqD->stageData_mp[stage].useIntrinsicSlice; // whether we are using intrinsic slicing for this stage
  int depth_sum = depth_x + depth_y;

  comp_mp one_minus_t, gamma_t;
  point_mp orig_vars_x, orig_vars_y;
  vec_mp F;
  mat_mp Jv_F, Jp_temp;

  // initialize MP
  init_mp(one_minus_t); init_mp(gamma_t);
  init_point_mp(orig_vars_x, 0); init_point_mp(orig_vars_y, 0);
  init_vec_mp(F, 0);
  init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jp_temp, 0, 0);

  // find 1 - t
  set_one_mp(one_minus_t);
  sub_mp(one_minus_t, one_minus_t, pathVars);
  // find gamma * t
  mul_mp(gamma_t, BED->EqD->gamma_mp, pathVars);

  if (useIntrinsic)
  { // setup B = B0 + t*(gamma*B1 - B0) & p = p0 + t*(gamma*p1 - p0), Jp_temp = (gamma*p1 - p0) + (gamma*B1 - B0)*vars
    vec_mp p;
    mat_mp B;

    init_vec_mp(p, BED->EqD->stageData_mp[stage].p1->size);
    init_mat_mp(B, BED->EqD->stageData_mp[stage].B1->rows, BED->EqD->stageData_mp[stage].B1->cols);
    change_size_mat_mp(Jp_temp, BED->EqD->stageData_mp[stage].p1->size, 1);

    B->rows = BED->EqD->stageData_mp[stage].B1->rows;
    B->cols = BED->EqD->stageData_mp[stage].B1->cols; // vars->size
    Jp_temp->rows = p->size = BED->EqD->stageData_mp[stage].p1->size; // == B->rows
    Jp_temp->cols = 1;

    for (i = 0; i < B->rows; i++)
    {
      set_zero_mp(&Jp_temp->entry[i][0]);
      for (j = 0; j < B->cols; j++)
      { // find (gamma*B1 - B0)[i][j]
        mul_mp(&B->entry[i][j], BED->EqD->gamma_mp, &BED->EqD->stageData_mp[stage].B1->entry[i][j]);
        sub_mp(&B->entry[i][j], &B->entry[i][j], &BED->EqD->stageData_mp[stage].B0->entry[i][j]);

        // Jp[i] += (gamma*B1 - B0)[i][j] * vars[j]
        sum_mul_mp(&Jp_temp->entry[i][0], &B->entry[i][j], &vars->coord[j]);

        // setup B[i][j]
        mul_mp(&B->entry[i][j], &B->entry[i][j], pathVars);
        add_mp(&B->entry[i][j], &B->entry[i][j], &BED->EqD->stageData_mp[stage].B0->entry[i][j]);
      }
      // find (gamma*p1 - p0)[i]
      mul_mp(&p->coord[i], BED->EqD->gamma_mp, &BED->EqD->stageData_mp[stage].p1->coord[i]);
      sub_mp(&p->coord[i], &p->coord[i], &BED->EqD->stageData_mp[stage].p0->coord[i]);

      // Jp[i] += (gamma*p1 - p0)[i]
      add_mp(&Jp_temp->entry[i][0], &Jp_temp->entry[i][0], &p->coord[i]);

      // setup p[i]
      mul_mp(&p->coord[i], &p->coord[i], pathVars);
      add_mp(&p->coord[i], &p->coord[i], &BED->EqD->stageData_mp[stage].p0->coord[i]);
    }

    // convert vars to [x,y]
    intrinsicToExtrinsic_mp(orig_vars_x, vars, B, p);

    // seup x & y
    change_size_point_mp(orig_vars_y, num_vars);
    orig_vars_x->size = orig_vars_y->size = num_vars;
    for (i = 0; i < num_vars; i++)
    {
      set_mp(&orig_vars_y->coord[i], &orig_vars_x->coord[i + num_vars]);
    }

    // allocate space for funcVals & Jv_F
    change_size_vec_mp(funcVals, depth_sum);
    change_size_mat_mp(Jv_F, depth_sum, 2 * num_vars);
    funcVals->size = Jv_F->rows = depth_sum;
    Jv_F->cols = 2 * num_vars;

    // evaluate the function for the x variables
    eqbyeq_square_eval_mp(F, parVals, parDer, Jv, Jp, orig_vars_x, pathVars, BED, 0, depth_x);

    // setup top of funcVals & Jv_F
    for (i = 0; i < depth_x; i++)
    { // setup funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);
      // setup Jv_F
      for (j = 0; j < num_vars; j++)
      {
        set_mp(&Jv_F->entry[i][j], &Jv->entry[i][j]);
        set_zero_mp(&Jv_F->entry[i][j + num_vars]); // only depends on x variables
      }
    }

    // evaluate the function for the y variables
    eqbyeq_square_eval_mp(F, parVals, parDer, Jv, Jp, orig_vars_y, pathVars, BED, depth_x, depth_sum);

    // setup bottom of funcVals & Jv_F
    for (i = depth_x; i < depth_sum; i++)
    { // setup funcVals
      set_mp(&funcVals->coord[i], &F->coord[i]);
      // setup Jv_F
      for (j = 0; j < num_vars; j++)
      {
        set_mp(&Jv_F->entry[i][j + num_vars], &Jv->entry[i][j]);
        set_zero_mp(&Jv_F->entry[i][j]); // only depends on y variables
      }
    }

    // we already have funcVals = [F_x,F_y]

    // use the structure of Jv_F to find Jp = Jv_F * Jp_temp and Jv = Jv_F * B
    change_size_mat_mp(Jv, depth_sum, B->cols);
    change_size_mat_mp(Jp, depth_sum, 1);
    Jv->rows = Jp->rows = depth_sum;
    Jv->cols = B->cols; // == vars->size
    Jp->cols = Jp_temp->cols; // == 1
    for (i = 0; i < depth_sum; i++)
    { // find offset
      offset = (i < depth_x) ? 0 : num_vars;

      for (j = 0; j < B->cols; j++)
      { // find Jv[i][j]
        set_zero_mp(&Jv->entry[i][j]);
        for (k = 0; k < num_vars; k++)
        {
          sum_mul_mp(&Jv->entry[i][j], &Jv_F->entry[i][k + offset], &B->entry[k + offset][j]);
        }
      }
      // find Jp[i][0]
      set_zero_mp(&Jp->entry[i][0]);
      for (k = 0; k < num_vars; k++)
      {
        sum_mul_mp(&Jp->entry[i][0], &Jv_F->entry[i][k + offset], &Jp_temp->entry[k + offset][0]);
      }
    }

    // clear memory
    clear_vec_mp(p);
    clear_mat_mp(B);
  }
  else
  { // using extrinsic formulation
    comp_mp tempComp;
    vec_mp patchValues;
    mat_mp Jv_Patch;

    init_mp(tempComp);
    init_vec_mp(patchValues, 0);
    init_mat_mp(Jv_Patch, 0, 0);

    // setup the sizes for funcVals, Jv & Jp
    change_size_vec_mp(funcVals, 2 * num_vars);
    change_size_mat_mp(Jv, 2 * num_vars, 2 * num_vars);
    change_size_mat_mp(Jp, 2 * num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = 2 * num_vars;
    Jp->cols = 1;

    // split vars into x & y
    change_size_vec_mp(orig_vars_x, num_vars);
    change_size_vec_mp(orig_vars_y, num_vars);
    orig_vars_x->size = orig_vars_y->size = num_vars;
    for (i = 0; i < num_vars; i++)
    {
      set_mp(&orig_vars_x->coord[i], &vars->coord[i]);
      set_mp(&orig_vars_y->coord[i], &vars->coord[i+num_vars]);
    }

    // evaluate the function for the x variables
    eqbyeq_square_eval_mp(F, parVals, parDer, Jv_F, Jp_temp, orig_vars_x, pathVars, BED, 0, depth_x);
    // evaluate the patch for the x variables
    patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp_temp, orig_vars_x, pathVars, &BED->patch);

    // setup the top of funcVals, Jv & Jp
    // calculate funcVals(top) = [F,L_x(top)*gamma*t + (1-t)(x-y),L_x(rest),patch_x]
    // calculate Jp(top) = [0,-(x-y) + gamma*L_x(top),0,0]
    // calculate Jv(top) = [Jv_F,Jv_linears(top)*gamma*t + (1-t)d(x-y),Jv_linears(rest),Jv_Patch]
    for (i = 0; i < num_vars; i++)
      if (i < depth_x) // functions already solved for
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
          // y-coordinates
          set_zero_mp(&Jv->entry[i][j + num_vars]);
        }
      }
      else if (i < depth_sum) // linears moving
      { // initialize tempComp - to be linear
        set_zero_mp(tempComp);
        // evaluate the ith linear in x-variables
        for (j = 0; j < num_vars; j++)
        { // add on to tempComp
          sum_mul_mp(tempComp, BED->EqD->coeff_mp[i][j], &orig_vars_x->coord[j]);
          // setup Jv - x-coordinates
          mul_mp(&Jv->entry[i][j], BED->EqD->coeff_mp[i][j], gamma_t); // coeff * gamma * t
          // setup Jv - y-coordinates
          set_zero_mp(&Jv->entry[i][j + num_vars]);
        }
        // funcVals
        sub_mp(&funcVals->coord[i], &orig_vars_x->coord[i], &orig_vars_y->coord[i]); // x_i - y_i
        mul_mp(&funcVals->coord[i], &funcVals->coord[i], one_minus_t); // (x_i - y_i) * (1-t)
        sum_mul_mp(&funcVals->coord[i], tempComp, gamma_t); // += linear * gamma * t
        // Jp
        sub_mp(&Jp->entry[i][0], &orig_vars_y->coord[i], &orig_vars_x->coord[i]); // y_i - x_i
        sum_mul_mp(&Jp->entry[i][0], tempComp, BED->EqD->gamma_mp); // (y_i - x_i) + gamma * linear
        // adjust Jv for i
        add_mp(&Jv->entry[i][i], &Jv->entry[i][i], one_minus_t);
        sub_mp(&Jv->entry[i][i + num_vars], &Jv->entry[i][i + num_vars], one_minus_t);
      }
      else if (i < num_funcs) // linears not moving
      { // funcVals - evaluate linear for x-variables
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // add on to funcVals
          sum_mul_mp(&funcVals->coord[i], BED->EqD->coeff_mp[i][j], &orig_vars_x->coord[j]);
          // setup Jv - x-coordinates
          set_mp(&Jv->entry[i][j], BED->EqD->coeff_mp[i][j]);
          // setup Jv - y-coordinates
          set_zero_mp(&Jv->entry[i][j + num_vars]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else // patch
      { // funcVals
        set_mp(&funcVals->coord[i], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_mp(&Jv->entry[i][j], &Jv_Patch->entry[i - num_funcs][j]);
          // y-coordinates
          set_zero_mp(&Jv->entry[i][j + num_vars]);
        }
      }

    // evaluate the function for the y variables
    eqbyeq_square_eval_mp(F, parVals, parDer, Jv_F, Jp_temp, orig_vars_y, pathVars, BED, depth_x, depth_sum);
    // evaluate the patch for the y variables
    patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp_temp, orig_vars_y, pathVars, &BED->patch);

    // setup the bottom of funcVals, Jv & Jp
    // calculate funcVals(bottom) = [F_y,L_y(bottom)*gamma*t + (1-t)(x-y),L_y(rest),patch_y]
    // calculate Jp(bottom) = [0,-(x-y) + gamma*L_y(bottom),0,0]
    // calculate Jv(bottom) = [Jv_F,Jv_linears(bottom)*gamma*t + (1-t)d(x-y),Jv_linears(rest),Jv_Patch]
    for (i = 0; i < num_vars; i++)
      if (i < depth_x) // linears moving
      { // initialize tempComp - to be linear
        set_zero_mp(tempComp);
        // evaluate the ith linear in x-variables
        for (j = 0; j < num_vars; j++)
        { // add on to tempComp
          sum_mul_mp(tempComp, BED->EqD->coeff_mp[i][j], &orig_vars_y->coord[j]);
          // setup Jv - y-coordinates
          mul_mp(&Jv->entry[i + num_vars][j + num_vars], BED->EqD->coeff_mp[i][j], gamma_t); // coeff * gamma * t
          // setup Jv - x-coordinates
          set_zero_mp(&Jv->entry[i + num_vars][j]);
        }
        // funcVals
        sub_mp(&funcVals->coord[i + num_vars], &orig_vars_x->coord[i], &orig_vars_y->coord[i]); // x_i - y_i
        mul_mp(&funcVals->coord[i + num_vars], &funcVals->coord[i + num_vars], one_minus_t); // (x_i - y_i) * (1-t)
        sum_mul_mp(&funcVals->coord[i + num_vars], tempComp, gamma_t); // += linear * gamma * t
        // Jp
        sub_mp(&Jp->entry[i + num_vars][0], &orig_vars_y->coord[i], &orig_vars_x->coord[i]); // y_i - x_i
        sum_mul_mp(&Jp->entry[i + num_vars][0], tempComp, BED->EqD->gamma_mp); // (y_i - x_i) + gamma * linear
        // adjust Jv for i
        add_mp(&Jv->entry[i + num_vars][i], &Jv->entry[i + num_vars][i], one_minus_t);
        sub_mp(&Jv->entry[i + num_vars][i + num_vars], &Jv->entry[i + num_vars][i + num_vars], one_minus_t);
      }
      else if (i < depth_sum) // functions already solved for
      { // funcVals
        set_mp(&funcVals->coord[i + num_vars], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i + num_vars][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_zero_mp(&Jv->entry[i + num_vars][j]);
          // y-coordinates
          set_mp(&Jv->entry[i + num_vars][j + num_vars], &Jv_F->entry[i][j]);
        }
      }
      else if (i < num_funcs) // linears not moving
      { // funcVals
        set_zero_mp(&funcVals->coord[i + num_vars]);
        for (j = 0; j < num_vars; j++)
        { // add on to funcVals
          sum_mul_mp(&funcVals->coord[i + num_vars], BED->EqD->coeff_mp[i][j], &orig_vars_y->coord[j]);
          // setup Jv - y-coordinates
          set_mp(&Jv->entry[i + num_vars][j + num_vars], BED->EqD->coeff_mp[i][j]);
          // setup Jv - x-coordinates
          set_zero_mp(&Jv->entry[i + num_vars][j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i + num_vars][0]);
      }
      else // patch
      { // funcVals
        set_mp(&funcVals->coord[i + num_vars], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_mp(&Jp->entry[i + num_vars][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        { // x-coordinates
          set_zero_mp(&Jv->entry[i + num_vars][j]);
          // y-coordinates
          set_mp(&Jv->entry[i + num_vars][j + num_vars], &Jv_Patch->entry[i - num_funcs][j]);
        }
      }

    // clear memory
    clear_mp(tempComp);
    clear_vec_mp(patchValues);
    clear_mat_mp(Jv_Patch);
  }

  // set parVals & parDer correctly
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars);   // s = t
  set_one_mp(&parDer->coord[0]);          // ds/dt = 1

  // clear MP
  clear_mp(one_minus_t); clear_mp(gamma_t);
  clear_point_mp(orig_vars_x); clear_point_mp(orig_vars_y);
  clear_vec_mp(F);
  clear_mat_mp(Jv_F); clear_mat_mp(Jp_temp);

  return 0;
}

int stage_sort_eqbyeq_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: function evaluation to sort the stage                  *
*   we know x = y  and t = 0                                    *
\***************************************************************/
// since we need to do SVD, we use this evaluator so that the slices 'x - y' are not used
{
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;
  int i, j, stage = BED->EqD->curr_stage_num, num_vars = BED->EqD->num_vars, num_funcs = BED->EqD->num_funcs;
  int depth_x = BED->EqD->stageData_mp[stage].depth_x,  // depth for the x-variables
      depth_y = BED->EqD->stageData_mp[stage].depth_y,  // depth for the y-variables
      useIntrinsic = BED->EqD->stageData_mp[stage].useIntrinsicSlice; // whether we are using intrinsic slicing for this stage
  int depth_sum = depth_x + depth_y;

  if (useIntrinsic)
  {
    point_mp orig_vars;
    init_point_mp(orig_vars, 0);

    // now convert to the original variables (p + B*vars)
    intrinsicToExtrinsic_mp(orig_vars, vars, BED->EqD->stageData_mp[stage].B, BED->EqD->stageData_mp[stage].p);

    // evaluate the function
    eqbyeq_square_eval_mp(funcVals, parVals, parDer, Jv, Jp, orig_vars, pathVars, BED, 0, depth_sum);

    // pick out the correct rows of funcVals and Jv
    funcVals->size = Jv->rows = depth_sum;

    // convert Jv to the input vars
    mat_mul_mp(Jv, Jv, BED->EqD->stageData_mp[stage].B);

    // setup Jp == 0
    change_size_mat_mp(Jp, depth_sum, 1);
    Jp->rows = depth_sum;
    Jp->cols = 1;
    for (i = 0; i < depth_sum; i++)
    {
      set_zero_mp(&Jp->entry[i][0]);
    }

    // clear MP
    clear_point_mp(orig_vars);
  }
  else
  { // use extrinsic formulation
    vec_mp F, patchValues;
    mat_mp Jv_F, Jv_Patch;

    init_vec_mp(F, 0); init_vec_mp(patchValues, 0);
    init_mat_mp(Jv_F, 0, 0); init_mat_mp(Jv_Patch, 0, 0);

    // evaluate the function
    eqbyeq_square_eval_mp(F, parVals, parDer, Jv_F, Jp, vars, pathVars, BED, 0, depth_sum);
    // evaluate the patch
    patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);

    // calculate funcVals = [F, L, patch]
    // calculate Jp(bottom) = [0, 0, 0]
    // calculate Jv(bottom) = [Jv_F, Jv_linears, Jv_Patch]
    change_size_vec_mp(funcVals, num_vars);
    change_size_mat_mp(Jv, num_vars, num_vars);
    change_size_mat_mp(Jp, num_vars, 1);
    funcVals->size = Jv->rows = Jv->cols = Jp->rows = num_vars;
    Jp->cols = 1;
    for (i = 0; i < num_vars; i++)
      if (i < depth_sum) // F
      { // funcVals
        set_mp(&funcVals->coord[i], &F->coord[i]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], &Jv_F->entry[i][j]);
        }
      }
      else if (i < num_funcs) // linears
      { // funcVals - linears
        set_zero_mp(&funcVals->coord[i]);
        for (j = 0; j < num_vars; j++)
        { // update funcVals
          sum_mul_mp(&funcVals->coord[i], BED->EqD->coeff_mp[i][j], &vars->coord[j]);
          // setup Jv
          set_mp(&Jv->entry[i][j], BED->EqD->coeff_mp[i][j]);
        }
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
      }
      else // patch
      { // funcVals
        set_mp(&funcVals->coord[i], &patchValues->coord[i - num_funcs]);
        // Jp
        set_zero_mp(&Jp->entry[i][0]);
        // Jv
        for (j = 0; j < num_vars; j++)
        {
          set_mp(&Jv->entry[i][j], &Jv_Patch->entry[i - num_funcs][j]);
        }
      }

    // clear memory
    clear_vec_mp(F); clear_vec_mp(patchValues);
    clear_mat_mp(Jv_F); clear_mat_mp(Jv_Patch);
  }

  // setup parVals & parDer correctly
  change_size_vec_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
  set_zero_mp(&parVals->coord[0]);   // s = t = 0
  set_one_mp(&parDer->coord[0]);     // ds/dt = 1

  return 0;
}

int eqbyeq_square_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp vars, comp_mp pathVars, void const *ED, int startFunc, int endFunc)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: square function evaluation for equation-by-equation    * 
\***************************************************************/
{
  basic_eval_data_mp *BED = (basic_eval_data_mp *)ED;

  if (BED->EqD->noChanges)
  { // evaluate the SLP efficiently - from function startFunc to endFunc
    evalProg_eff_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, BED->squareSystem.Prog, startFunc, endFunc, BED->EqD->startSub, BED->EqD->endSub, BED->EqD->startFunc, BED->EqD->endFunc, BED->EqD->startJvsub, BED->EqD->endJvsub, BED->EqD->startJv, BED->EqD->endJv, BED->EqD->subFuncsBelow);
  }
  else
  { // evaluate the square system
    square_system_eval_mp(funcVals, parVals, parDer, Jv, Jp, vars, pathVars, &BED->squareSystem);
  }
  
  return 0;
} 

