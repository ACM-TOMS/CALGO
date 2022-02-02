// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "dimbydim.h"

void dimbydimOutputChart(codim_t *CD, FILE *fp, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int codim_index, total_paths = 0, num_codim = CD->num_codim;

  fprintf(fp, "\n\n*********** Dimension-by-Dimension Summary ***********\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|codim|   paths   |witness superset| nonsingular | singular |nonsolutions| inf endpoints | other bad endpoints\n");
  fprintf(fp, "----------------------------------------------------------------------------------------------------------------\n");

  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { 
    // add to the total paths tracked
    total_paths += CD->codim[codim_index].num_paths;

    fprintf(fp, "| %-4d|   %-8d|   %-13d|  %-11d|  %-8d|  %-10d|   %-12d|  %d\n", CD->codim[codim_index].codim, CD->codim[codim_index].num_paths, CD->codim[codim_index].num_superset, CD->codim[codim_index].num_nonsing, CD->codim[codim_index].num_sing, CD->codim[codim_index].num_nonsolns, CD->codim[codim_index].num_inf, CD->codim[codim_index].num_bad);
  }
  fprintf(fp, "----------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|total|   %d\n\n", total_paths);

  fprintf(fp, "****************************************************\n\n");

  return;
}

void dimbydimFindOrigVarsDehom_d(point_d orig_vars_d, point_d dehom_d, point_d P_d, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the original coordinates and dehom  *
\***************************************************************/
{
  int i;

  if (CD->codim[codim_index].useIntrinsicSlice)
  { // first convert to extrinsic
    point_d ext_d;
    init_point_d(ext_d, 0);

    intrinsicToExtrinsic_d(ext_d, P_d, CD->codim[codim_index].B_d, CD->codim[codim_index].p_d);

    // change ext_d if we had an 'intrinsic' homogeneous coordinate - when the variable group was already homogenized
    if (CD->PPD.num_hom_var_gp)
    { // find the 'intrinsic' dehom coordinate
      comp_d dehomCoord;
      set_d(dehomCoord, CD->codim[codim_index].homVarConst_d);
      for (i = 0; i < ext_d->size; i++)
      {
        sum_mul_d(dehomCoord, &CD->codim[codim_index].H_d->coord[i], &ext_d->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (isnan(dehomCoord->r) || isnan(dehomCoord->i) || isinf(dehomCoord->r) || isinf(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < ext_d->size; i++)
        {
          set_double_d(&ext_d->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (d_abs_d(dehomCoord) == 0)
        { // generate a random perturbation so that we can divide
          get_comp_rand_d(dehomCoord);
          mul_rdouble_d(dehomCoord, dehomCoord, 1e-16);
          recip_d(dehomCoord, dehomCoord);
        }
        else
        { // reciprocate dehomCoord
          recip_d(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < ext_d->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_d(&ext_d->coord[i], &ext_d->coord[i], dehomCoord);
        }
      }
    }

    // find orig_vars_d from ext_d
    if (CD->new_variables != CD->orig_variables)
    {
      mul_mat_vec_d(orig_vars_d, CD->C_d, ext_d);
    }
    else
    {
      point_cp_d(orig_vars_d, ext_d);
    }

    clear_point_d(ext_d);
  }
  else
  { // find the actual 'new_variables'
    if (CD->PPD.num_var_gp)
    { // the variable group was un-homogenized - find orig_vars as expected
   
      // convert to original coordinates
      if (CD->new_variables != CD->orig_variables)
      {
        mul_mat_vec_d(orig_vars_d, CD->C_d, P_d);
      }
      else
      {
        point_cp_d(orig_vars_d, P_d);
      }
    }
    else
    { // the variable group was homogenized - remove the intrinsic dehom coordinate and then find orig_vars as expected
      comp_d dehomCoord;
      point_d tempPoint;

      init_point_d(tempPoint, P_d->size);
      tempPoint->size = P_d->size;
      // find the 'intrinsic' dehom coordinate
      set_d(dehomCoord, CD->codim[codim_index].homVarConst_d);
      for (i = 0; i < P_d->size; i++)
      {
        sum_mul_d(dehomCoord, &CD->codim[codim_index].H_d->coord[i], &P_d->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (isnan(dehomCoord->r) || isnan(dehomCoord->i) || isinf(dehomCoord->r) || isinf(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < tempPoint->size; i++)
        {
          set_double_d(&tempPoint->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (d_abs_d(dehomCoord) == 0)
        { // generate a random perturbation so that we can divide
          get_comp_rand_d(dehomCoord);
          mul_rdouble_d(dehomCoord, dehomCoord, 1e-16);
          recip_d(dehomCoord, dehomCoord);
        }
        else
        { // reciprocate dehomCoord
          recip_d(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < tempPoint->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_d(&tempPoint->coord[i], &P_d->coord[i], dehomCoord);
        }
      }

      // convert tempPoint to orig_vars_d 
      if (CD->new_variables != CD->orig_variables)
      {
        mul_mat_vec_d(orig_vars_d, CD->C_d, tempPoint);
      }
      else
      {
        point_cp_d(orig_vars_d, tempPoint);
      }

      clear_point_d(tempPoint);
    }
  }

  // find dehom_d using orig_vars_d
  getDehomPoint_d(dehom_d, orig_vars_d, orig_vars_d->size, &CD->PPD);

  return;
}

void dimbydimFindOrigVarsDehom_mp(point_mp orig_vars_mp, point_mp dehom_mp, point_mp P_mp, codim_t *CD, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the original coordinates and dehom  *
\***************************************************************/
{
  int i;

  if (CD->codim[codim_index].useIntrinsicSlice)
  { // first convert to extrinsic
    point_mp ext_mp;
    init_point_mp2(ext_mp, 0, CD->curr_precision);

    intrinsicToExtrinsic_mp(ext_mp, P_mp, CD->codim[codim_index].B_mp, CD->codim[codim_index].p_mp);

    // change ext_mp if we had an 'intrinsic' homogeneous coordinate - when the variable group was already homogenized
    if (CD->PPD.num_hom_var_gp)
    { // find the 'intrinsic' dehom coordinate
      comp_mp dehomCoord;
      init_mp2(dehomCoord, CD->curr_precision);
      set_mp(dehomCoord, CD->codim[codim_index].homVarConst_mp);
      for (i = 0; i < ext_mp->size; i++)
      {
        sum_mul_mp(dehomCoord, &CD->codim[codim_index].H_mp->coord[i], &ext_mp->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (!mpfr_number_p(dehomCoord->r) || !mpfr_number_p(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < ext_mp->size; i++)
        {
          set_double_mp(&ext_mp->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (mpfr_zero_p(dehomCoord->r) && mpfr_zero_p(dehomCoord->i))
        { // generate a random perturbation so that we can divide
          mpf_t epsilon;
          mpf_init2(epsilon, CD->curr_precision);

          get_comp_rand_mp(dehomCoord);
          mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(CD->curr_precision), __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(dehomCoord, dehomCoord, epsilon);
          recip_mp(dehomCoord, dehomCoord);

          mpf_clear(epsilon);
        }
        else
        { // reciprocate dehomCoord
          recip_mp(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < ext_mp->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_mp(&ext_mp->coord[i], &ext_mp->coord[i], dehomCoord);
        }
      }

      clear_mp(dehomCoord);
    }    

    // find orig_vars_mp from ext_mp
    if (CD->new_variables != CD->orig_variables)
    {
      mul_mat_vec_mp(orig_vars_mp, CD->C_mp, ext_mp);
    }
    else
    {
      point_cp_mp(orig_vars_mp, ext_mp);
    }

    clear_point_mp(ext_mp);
  }
  else
  { // find the actual 'new_variables'
    if (CD->PPD.num_var_gp)
    { // the variable group was un-homogenized - find orig_vars as expected

      // convert to original homogeneous coordinates
      if (CD->new_variables != CD->orig_variables)
      {
        mul_mat_vec_mp(orig_vars_mp, CD->C_mp, P_mp);
      }
      else
      {
        point_cp_mp(orig_vars_mp, P_mp);
      }
    }
    else
    { // the variable group was homogenized - remove the intrinisic dehom coordinate and then find orig_vars as expected
      comp_mp dehomCoord;
      point_mp tempPoint;

      init_mp2(dehomCoord, CD->curr_precision);
      init_point_mp2(tempPoint, P_mp->size, CD->curr_precision);
      tempPoint->size = P_mp->size;
      // find the 'intrinsic' dehom coordinate
      set_mp(dehomCoord, CD->codim[codim_index].homVarConst_mp);
      for (i = 0; i < P_mp->size; i++)
      {
        sum_mul_mp(dehomCoord, &CD->codim[codim_index].H_mp->coord[i], &P_mp->coord[i]);
      }

      // make sure dehomCoord is a valid number
      if (!mpfr_number_p(dehomCoord->r) || !mpfr_number_p(dehomCoord->i))
      { // not a valid number
        for (i = 0; i < tempPoint->size; i++)
        {
          set_double_mp(&tempPoint->coord[i], -1e199, -1e199);
        }
      }
      else
      { // we have a number - determine if it is 0
        if (mpfr_zero_p(dehomCoord->r) && mpfr_zero_p(dehomCoord->i)) 
        { // generate a random perturbation so that we can divide
          mpf_t epsilon;
          mpf_init2(epsilon, CD->curr_precision);

          get_comp_rand_mp(dehomCoord);
          mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(CD->curr_precision), __gmp_default_rounding_mode);
          mpf_ui_div(epsilon, 1, epsilon);
          mul_rmpf_mp(dehomCoord, dehomCoord, epsilon);
          recip_mp(dehomCoord, dehomCoord);
 
          mpf_clear(epsilon);
        }
        else
        { // reciprocate dehomCoord
          recip_mp(dehomCoord, dehomCoord);
        }

        // find the coordinates
        for (i = 0; i < tempPoint->size; i++)
        { // multiply by recip to produce de-hom coord.
          mul_mp(&tempPoint->coord[i], &P_mp->coord[i], dehomCoord);
        }
      }

      // convert tempPoint to orig_vars_mp
      if (CD->new_variables != CD->orig_variables)
      {
        mul_mat_vec_mp(orig_vars_mp, CD->C_mp, tempPoint);
      }
      else
      {
        point_cp_mp(orig_vars_mp, tempPoint);
      }

      clear_mp(dehomCoord);
      clear_point_mp(tempPoint);
    }
  }

  // find dehom_mp from orig_vars_mp
  getDehomPoint_mp(dehom_mp, orig_vars_mp, orig_vars_mp->size, &CD->PPD);

  return;
}

