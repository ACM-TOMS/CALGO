// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

// Verify that the given point is a smooth point on a generically reduced irreducible component 

void setup_slice_points(membership_slice_moving_t *sliceMover, int numSlices, int numVars, int numPoints, point_d *Pts_d, point_mp *Pts_mp, int prec, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup slices through the points                        *
\***************************************************************/
{
  // initialize B
  if (T->MPType == 0)
  { // setup B_d
    init_mat_d(sliceMover->B_d, numSlices, numVars);
    sliceMover->B_d->rows = numSlices;
    sliceMover->B_d->cols = numVars;
  }
  else if (T->MPType == 1)
  { // setup B_mp
    init_mat_mp(sliceMover->B_mp, numSlices, numVars);
    sliceMover->B_mp->rows = numSlices;
    sliceMover->B_mp->cols = numVars;
  }
  else
  { // setup B_d, B_mp, B_rat
    init_mat_d(sliceMover->B_d, numSlices, numVars);
    init_mat_mp(sliceMover->B_mp, numSlices, numVars);
    init_mat_rat(sliceMover->B_rat, numSlices, numVars);
    sliceMover->B_d->rows = sliceMover->B_mp->rows = numSlices;
    sliceMover->B_d->cols = sliceMover->B_mp->cols = numVars;
  }

  if (numPoints == 1)
  { // setup completely random
    if (T->MPType == 0)
    { // setup B_d
      make_matrix_random_d(sliceMover->B_d, numSlices, numVars);
    }
    else if (T->MPType == 1)
    { // setup B_mp
      make_matrix_random_mp(sliceMover->B_mp, numSlices, numVars, T->Precision);
    }
    else
    { // setup B_d, B_mp, B_rat
      make_matrix_random_rat(sliceMover->B_d, sliceMover->B_mp, sliceMover->B_rat, numSlices, numVars, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);
    }
  }
  else
  { // setup based on the given points
    int i, j, corank, nullity;

    if (T->MPType == 0 || T->MPType == 2)
    { // setup using Pts_d
      int digits = 16, prec = T->Precision;
      comp_d rand;
      vec_d x, b;
      mat_d P, U, E, V;

      init_vec_d(x, 0);
      init_vec_d(b, numPoints);
      init_mat_d(P, numPoints, numVars);
      init_mat_d(U, 0, 0);
      init_mat_d(E, 0, 0);
      init_mat_d(V, 0, 0);
      b->size = P->rows = numPoints;
      P->cols = numVars;

      // setup the rand, b, P
      get_comp_rand_d(rand);
      for (i = 0; i < numPoints; i++)
      {
        set_d(&b->coord[i], rand);
        for (j = 0; j < numVars; j++)
          set_d(&P->entry[i][j], &Pts_d[i]->coord[j]);
      }

      // compute SVD & nullity
      corank = svd_jacobi_d_prec(U, E, V, P, MAX(T->final_tol_times_mult, T->sing_val_zero_tol));

      if (numPoints >= numVars)
        nullity = corank;
      else
        nullity = (numVars - numPoints) + corank;

      // verify that we have enough freedom
      if (numSlices > nullity)
      {
        printf("ERROR: The null space of the points is not large enough (%d > %d)!\n", numSlices, nullity); 
        bexit(ERROR_CONFIGURATION);
      }

      // solve P*x = b using SVD
      transpose_d(U, U);
      mul_mat_vec_d(x, U, b);
      j = E->cols - nullity;
      for (i = 0; i < j; i++)
        recip_d(&E->entry[i][i], &E->entry[i][i]);
      transpose_d(E, E);
      mul_mat_vec_d(x, E, x);
      mul_mat_vec_d(x, V, x);

      if (T->MPType == 0)
      { // setup B_d = random element in null space of P (spanned by last columns of V) + x
        for (i = 0; i < numSlices; i++)
        {
          make_vec_random_d(b, numVars);
          for (j = 0; j < numVars - nullity; j++)
            set_zero_d(&b->coord[j]);
          mul_mat_vec_d(b, V, b);
          for (j = 0; j < numVars; j++)
          {
            add_d(&sliceMover->B_d->entry[i][j], &x->coord[j], &b->coord[j]);
          }
        }
      }
      else
      { // setup B_d, B_mp, B_rat = random element in null space of P (spanned by last columns of V) + x
        for (i = 0; i < numSlices; i++)
        {
          make_vec_random_d(b, numVars);
          for (j = 0; j < numVars - nullity; j++)
            set_zero_d(&b->coord[j]);
          mul_mat_vec_d(b, V, b);
          for (j = 0; j < numVars; j++)
          {
            add_d(&sliceMover->B_d->entry[i][j], &x->coord[j], &b->coord[j]);
            comp_d_to_mp_rat(&sliceMover->B_mp->entry[i][j], sliceMover->B_rat[i][j], &sliceMover->B_d->entry[i][j], digits, prec, 0, 0);
          }
        }
      }

      // clear memory
      clear_vec_d(x);
      clear_vec_d(b);
      clear_mat_d(U);
      clear_mat_d(E);
      clear_mat_d(V);
      clear_mat_d(P);
    }
    else  
    { // setup using Pts_mp
      comp_mp rand;
      vec_mp x, b;
      mat_mp P, U, E, V;

      init_mp(rand);
      init_vec_mp(x, 0);
      init_vec_mp(b, numPoints);
      init_mat_mp(P, numPoints, numVars);
      init_mat_mp(U, 0, 0);
      init_mat_mp(E, 0, 0);
      init_mat_mp(V, 0, 0);
      b->size = P->rows = numPoints;
      P->cols = numVars;

      // setup the rand, b, P
      get_comp_rand_mp(rand);
      for (i = 0; i < numPoints; i++)
      {
        set_mp(&b->coord[i], rand);
        for (j = 0; j < numVars; j++)
          set_mp(&P->entry[i][j], &Pts_mp[i]->coord[j]);
      }

      // compute SVD & nullity
      corank = svd_jacobi_mp_prec(U, E, V, P, MAX(T->final_tol_times_mult, T->sing_val_zero_tol), T->Precision);
      if (numPoints >= numVars)
        nullity = corank;
      else
        nullity = (numVars - numPoints) + corank;

      // verify that we have enough freedom
      if (numSlices > nullity)
      {
        printf("ERROR: The null space of the points is not large enough (%d > %d)!\n", numSlices, nullity);
        bexit(ERROR_CONFIGURATION);
      }

      // solve P*x = b using SVD
      transpose_mp(U, U);
      mul_mat_vec_mp(x, U, b);
      j = E->cols - nullity;
      for (i = 0; i < j; i++)
        recip_mp(&E->entry[i][i], &E->entry[i][i]);
      transpose_mp(E, E);
      mul_mat_vec_mp(x, E, x);
      mul_mat_vec_mp(x, V, x);

      // setup B_mp
      for (i = 0; i < numSlices; i++)
      {
        corank = numVars - 1 - i;
        for (j = 0; j < numVars; j++)
        {
          add_mp(&sliceMover->B_mp->entry[i][j], &x->coord[j], &V->entry[j][corank]);
        }
      }

      // clear memory
      clear_mp(rand);
      clear_vec_mp(x);
      clear_vec_mp(b);
      clear_mat_mp(U);
      clear_mat_mp(E);
      clear_mat_mp(V);
      clear_mat_mp(P);
    }
  }

  return;
}

void setup_slice_moving(membership_slice_moving_t *sliceMover, int nullity, int numPoints, point_d *Pts_d, point_mp *Pts_mp, int prec, tracker_config_t *T, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup sliceMover and adjust points to vanish on patch  *
\***************************************************************/
{
  int i, j, dim = nullity - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;

  // setup sliceMover
  sliceMover->curr_codim = T->numVars - nullity;
  sliceMover->orig_variables = W->orig_variables;
  sliceMover->curr_precision = W->curr_precision;

  // initialize A even though it will not be setup here
  if (T->MPType == 0 || T->MPType == 2)
  {
    init_mat_d(sliceMover->A_d, 0, 0);
  }
  if (T->MPType == 1 || T->MPType == 2)
  {
    init_mat_mp2(sliceMover->A_mp, 0, 0, sliceMover->curr_precision);
  }
  sliceMover->A_rat = NULL;
  sliceMover->A_d->rows = sliceMover->A_d->cols = sliceMover->A_mp->rows = sliceMover->A_mp->cols = 0;
  sliceMover->startSliceVec_init = 0;
  sliceMover->targetSliceVec_init = 0;

  // setup gamma
  if (T->MPType == 0)
  { // setup gamma_d
    get_comp_rand_d(sliceMover->gamma_d);
  }
  else if (T->MPType == 1)
  { // setup gamma_mp
    init_mp(sliceMover->gamma_mp);
    get_comp_rand_mp(sliceMover->gamma_mp);
  }
  else
  { // setup gamma_d, gamma_mp, gamma_rat
    get_comp_rand_rat(sliceMover->gamma_d, sliceMover->gamma_mp, sliceMover->gamma_rat, sliceMover->curr_precision, T->AMP_max_prec, 1, 1);
  }

  // setup p - random vector
  if (T->MPType == 0)
  { // setup _d
    init_vec_d(sliceMover->p_d, T->numVars);
    make_vec_random_d(sliceMover->p_d, T->numVars);
  }
  else if (T->MPType == 1)
  { // setup _mp
    init_vec_mp(sliceMover->p_mp, T->numVars);
    make_vec_random_mp(sliceMover->p_mp, T->numVars);
  }
  else
  { // setup _rat
    init_vec_d(sliceMover->p_d, T->numVars);
    init_vec_mp(sliceMover->p_mp, T->numVars);
    init_vec_rat(sliceMover->p_rat, T->numVars);
    make_vec_random_rat(sliceMover->p_d, sliceMover->p_mp, sliceMover->p_rat, T->numVars, sliceMover->curr_precision, T->AMP_max_prec, 0, 0);
  }

  // adjust Pts so that Pts * p  = 1
  if (prec < 64) 
  { // input data is _d
    comp_d tempComp;

    for (j = 0; j < numPoints; j++)
    { // compute Pts_d * p_d
      set_zero_d(tempComp);
      for (i = 0; i < T->numVars; i++)
        sum_mul_d(tempComp, &Pts_d[j]->coord[i], &sliceMover->p_d->coord[i]); 
    
      // normalize Pt_d
      for (i = 0; i < T->numVars; i++)
        div_d(&Pts_d[j]->coord[i], &Pts_d[j]->coord[i], tempComp);
    }
  }
  else
  { // input data is _mp
    comp_mp tempComp;
    init_mp2(tempComp, prec);

    for (j = 0; j < numPoints; j++)
    { // compute Pts_mp * p_mp
      set_zero_mp(tempComp);
      for (i = 0; i < T->numVars; i++)
        sum_mul_mp(tempComp, &Pts_mp[j]->coord[i], &sliceMover->p_mp->coord[i]); 

      // normalize Pt_mp
      for (i = 0; i < T->numVars; i++)
        div_mp(&Pts_mp[j]->coord[i], &Pts_mp[j]->coord[i], tempComp);
    }

    clear_mp(tempComp);
  }

  // setup B - random slices through points
  setup_slice_points(sliceMover, dim, T->numVars, numPoints, Pts_d, Pts_mp, prec, T);

  // setup K_rat
  sliceMover->K_rat = NULL;
  sliceMover->K_rows = sliceMover->K_cols = 0;

  return;
}

int isosingularDimTest(membership_slice_moving_t *sliceMover, point_d *endPts_d, point_mp *endPts_mp, int *endPts_prec, int nullity, int numPoints, point_d *Pts_d, point_mp *Pts_mp, int prec, tracker_config_t *T, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - success, 1 - error                         *
* NOTES: perform the isosingular dimension test                 *
\***************************************************************/
{
  int i, j, its = 0, rV = 0, maxIts = 20;
  endgame_data_t endPt;
  FILE *OUT = fopen("output_isosingular", "w"), *MIDOUT = fopen("midout_isosingular", "w");

  // initialize endPt
  init_endgame_data(&endPt, T->Precision);

  // setup sliceMover
  setup_slice_moving(sliceMover, nullity, numPoints, Pts_d, Pts_mp, prec, T, W);

  // setup the slice using the first point
  if (prec < 64)
    setup_slice_moving_slice(sliceMover, Pts_d[0], NULL, prec, T->MPType, T->AMP_max_prec);
  else
    setup_slice_moving_slice(sliceMover, NULL, Pts_mp[0], prec, T->MPType, T->AMP_max_prec);

  // final setup - setup gamma
  final_setup_slice_moving(sliceMover, W->Prog, T->MPType, T->AMP_max_prec, 1);

  // set startSliceVec == targetSliceVec & targetSliceVec == 0
  if (T->MPType == 0)
  {
    for (i = 0; i < sliceMover->targetSliceVec_d->size; i++)
    {
      set_d(&sliceMover->startSliceVec_d->coord[i], &sliceMover->targetSliceVec_d->coord[i]);
      set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
    }
  }
  else if (T->MPType == 1)
  {
    for (i = 0; i < sliceMover->targetSliceVec_mp->size; i++)
    {
      set_mp(&sliceMover->startSliceVec_mp->coord[i], &sliceMover->targetSliceVec_mp->coord[i]);
      set_zero_mp(&sliceMover->targetSliceVec_mp->coord[i]);
    }
  }
  else
  {
    for (i = 0; i < sliceMover->targetSliceVec_d->size; i++)
    {
      set_d(&sliceMover->startSliceVec_d->coord[i], &sliceMover->targetSliceVec_d->coord[i]);
      set_zero_d(&sliceMover->targetSliceVec_d->coord[i]);
      set_mp(&sliceMover->startSliceVec_mp->coord[i], &sliceMover->targetSliceVec_mp->coord[i]);
      set_zero_mp(&sliceMover->targetSliceVec_mp->coord[i]);
      set_rat(sliceMover->startSliceVec_rat[i], sliceMover->targetSliceVec_rat[i]);
      set_zero_rat(sliceMover->targetSliceVec_rat[i]);
    }
  }
 
  for (j = 0; j < numPoints && rV == 0; j++)
  { // try the jth point
    rV = 1;
    its = 0;
    while (its < maxIts && rV)
    { // move the slice
      if (prec < 64)
        rV = slice_moving_track(&endPt, sliceMover, Pts_d[j], NULL, prec, 0, 0, T, OUT, MIDOUT);
      else
        rV = slice_moving_track(&endPt, sliceMover, NULL, Pts_mp[j], prec, 0, 0, T, OUT, MIDOUT);

      its++;
    }

    // setup endPts_d[j] or endPts_mp[j]
    endPts_prec[j] = endPt.prec;
    if (endPts_prec[j] < 64)
    { // setup endPts_d[j]
      init_point_d(endPts_d[j], endPt.PD_d.point->size);
      vec_cp_d(endPts_d[j], endPt.PD_d.point);
    }
    else
    { // setup endPts_mp[j]
      init_point_mp2(endPts_mp[j], endPt.PD_d.point->size, endPt.prec);
      vec_cp_mp(endPts_mp[j], endPt.PD_mp.point);
    }

    if (!rV)
    { // need to verify we still have a solution
      int maxPrec = MAX(endPt.prec, endPt.last_approx_prec);

      if (maxPrec < 64)
      { // use _d
        rV = !nonsolutions_check_d(1, 0, endPt.PD_d.point, endPt.last_approx_d, endPt.PD_d.time, T->funcResTol, T->ratioTol, sliceMover->Prog);
      }
      else
      { // use _mp
        if (endPt.prec < 64)
        { // move to _mp
          setprec_point_mp(endPt.PD_mp.point, maxPrec);
          point_d_to_mp(endPt.PD_mp.point, endPt.PD_d.point);
        }
        
        if (endPt.last_approx_prec < 64)
        { // move to _mp
          setprec_point_mp(endPt.last_approx_mp, maxPrec);
          point_d_to_mp(endPt.last_approx_mp, endPt.last_approx_d);
        }
        
        rV = !nonsolutions_check_mp(1, 0, endPt.PD_mp.point, endPt.last_approx_mp, endPt.PD_mp.time, T->funcResTol, T->ratioTol, sliceMover->Prog);
      }
    }
  }

  // close OUT & MIDOUT
  fclose(OUT);
  fclose(MIDOUT);

  // delete midout_isosingular
  remove("midout_isosingular");

  // clear endPt
  clear_endgame_data(&endPt);

  return rV;
}





