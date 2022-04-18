// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

// this file implements the AMP tracker that is described in the adaptive precision papers
// of Bates, Hauenstein, Sommese & Wampler

int constant_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int euler_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int heun_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int runge_kutta_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));

int constant_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int euler_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int heun_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));
int runge_kutta_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int AMP2_success_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T);
int AMP2_criterion_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T);

int AMP2_success_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T, int *prec_decreases, int max_prec_decreases);
int AMP2_criterion_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T);


double minTime_d(int digits, comp_d t_d, comp_mp t_mp, int prec_in)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: minimum time based on the precision            *
* NOTES:                                                        *
\***************************************************************/
{
  double tempD = pow(10, 2 - digits) * (prec_in < 64 ? d_abs_d(t_d) : d_abs_mp(t_mp));

  // do not let double precision handle time smaller than 1e-150
  if (tempD < 1e-150)
    tempD = 1e-150;
 
  return tempD;
}

void minTime_mp(mpf_t minPrecTime, int digits, comp_d t_d, comp_mp t_mp, int prec_in)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: minimum time based on the precision            *
* NOTES:                                                        *
\***************************************************************/
{
  double tempD;
  mpf_t tempMPF;
  mpf_init(tempMPF);

  mpfr_ui_pow_ui(tempMPF, 10, digits - 2, __gmp_default_rounding_mode);
  mpf_ui_div(tempMPF, 1, tempMPF);

  if (prec_in < 64)
  {
    tempD = d_abs_d(t_d);
    mpf_set_d(minPrecTime, tempD);
  }
  else
  {
    mpf_abs_mp(minPrecTime, t_mp);
  }
  mpf_mul(minPrecTime, minPrecTime, tempMPF);

  mpf_clear(tempMPF);

  return;
}

int AMPtrack(point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, comp_d fT_d, comp_mp fT_mp, int time_prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: Tracks from Start->time to finalT                      *
* This will take care of all of the changes in precision and    *
* will always either reach finalT or the path tracking will have*
* completely failed                                             * 
\***************************************************************/
{
  if (T->odePredictor >= 3)
    return AMPtrack_error(Final_d, Final_mp, prec_out, time_first_increase, Start_d, Start_mp, prec_in, fT_d, fT_mp, time_prec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  int cont, tempInt, retVal = 0, curr_prec = prec_in, numSteps = 0, consecSuccess = 0, curr_digits = prec_to_digits(prec_in), need_to_refine = 0;
  int prec_decreases = 0, max_prec_decreases = 10;
  double absDistLeft_d, tempD, sizeProportion, norm_J = 0, norm_J_inv = 0, residual_old = 0;
  mpf_t absDistLeft_mp, tempMPF, minPrecTime_mp;
  comp_d nextTime_d, distLeft_d, finalT_d, dT_d;
  comp_mp nextTime_mp, distLeft_mp, finalT_mp, dT_mp, currStepSize_mp; // currStepSize here to extend precision easily
  vec_d dZ_d;
  vec_mp dZ_mp;
  point_data_d currPoint_d, tempPoint_d;
  point_data_mp currPoint_mp, tempPoint_mp;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  // initialize _d
  init_vec_d(dZ_d, 0);
  init_point_data_d(&currPoint_d, 0); init_point_data_d(&tempPoint_d, 0);
  init_eval_struct_d(e_d, 0, 0, 0);

  // make sure that everything is setup to the correct precision
  if (curr_prec < 64)
  { // we can start the MP at 64 bits
    T->Precision = 64;
  } 
  else
  { // start at curr_prec
    T->Precision = curr_prec;
  }
  // set everything to this precision
  initMP(T->Precision);
  change_prec(ED_mp, T->Precision);
  // initialize MP data structures to this precision
  init_point_data_mp(&currPoint_mp, 0); init_point_data_mp(&tempPoint_mp, 0);
  mpf_init(absDistLeft_mp); mpf_init(tempMPF); mpf_init(minPrecTime_mp);
  init_mp(nextTime_mp); init_mp(distLeft_mp); init_mp(finalT_mp); init_mp(dT_mp); init_mp(currStepSize_mp);
  init_vec_mp(dZ_mp, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  // setup tempPoint
  if (curr_prec < 64)
  { // start in double precision
    point_data_cp_d(&currPoint_d, Start_d);
    tempD = d_abs_d(currPoint_d.time);
  }
  else
  { // start in multi precision
    point_data_cp_mp(&currPoint_mp, Start_mp);
    tempD = d_abs_mp(currPoint_mp.time);
  }

  // setup finalT
  if (time_prec < 64)
  { // ending time is input in double precision
    set_d(finalT_d, fT_d);
    d_to_mp(finalT_mp, fT_d);
  }
  else
  { // end time is input in multi precision
    set_mp(finalT_mp, fT_mp);
    mp_to_d(finalT_d, fT_mp);
  }

  // initialize currStepSize_mp
  set_zero_mp(currStepSize_mp);
  mpf_set_d(currStepSize_mp->r, T->currentStepSize);

  // determine if precision needs to be raised based on the current newton tolerance, final tolerance (if in endgame) and current step size - ignore the criterion rules initially
  if (curr_prec < 64)
    tempD = -log10(T->currentStepSize);
  else
  {
    mpfr_log10(tempMPF, currStepSize_mp->r, __gmp_default_rounding_mode);
    tempD = -mpf_get_d(tempMPF); // -log10(currStepSize_mp)
  }
  AMP2_update(&cont, &tempInt, &T->currentStepSize, currStepSize_mp->r, curr_prec, 0, 0, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, tempD, tempD, currPoint_d.time, currPoint_mp.time);
  // NOTE: Since we are not given approximations to P0 and criterion C, this may tell us to decrease precision, but we will ignore it!! 
  // We only want to make sure that the precision is appropriate for the newton and final tolerance and the step size

  // check to see if precision needs to be increased
  if (tempInt > curr_prec)
  { // need to increase precision - in particular, have to use MP
    fprintf(OUT, "Changing Prec at beginning: %d\n", tempInt);
    T->Precision = tempInt;
    initMP(T->Precision);
    change_prec(ED_mp, T->Precision);
    setprec_point_data_mp(&tempPoint_mp, T->Precision);
    mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
    setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
    setprec_vec_mp(dZ_mp, T->Precision);
    setprec_eval_struct_mp(e_mp, T->Precision);

    // change precision on the ones with valuable information
    change_prec_mp2(currStepSize_mp, T->Precision);
    change_prec_mp2(finalT_mp, T->Precision);
    if (curr_prec < 64)
    { // convert _d to _mp
      setprec_point_data_mp(&currPoint_mp, T->Precision);
      convert_point_data_d_to_mp(&currPoint_mp, &currPoint_d);

      // store the time of this increase
      *time_first_increase = d_abs_d(currPoint_d.time);
    }
    else
    { // increase _mp
      change_prec_point_data_mp(&currPoint_mp, T->Precision);
    }

    // store the current precision now that everything is updated
    curr_prec = T->Precision;
    curr_digits = prec_to_digits(curr_prec);
  }

  // so, curr_prec is large enough based on the step size and the tolerance(s)

  // loop to initially refine using Newton iterations and check criterion inside of newton
  // this loop will terminate when currPoint has been refined or we run out of precision
  cont = 1;
  residual_old = 1e300;
  while (cont)
  {
    if (curr_prec < 64)
    { // perform the newton corrections in double precision - only updates after each successful iteration
      retVal = newton_d_amp(&norm_J, &norm_J_inv, curr_digits, currPoint_d.point, currPoint_d.time, T, OUT, &e_d, ED_d, eval_func_d);

      // check for errors
      if (retVal)
      { // need to increase precision (_d to _mp) - copy to 64-bit and continue
        fprintf(OUT, "Increase Prec (refine): %d to %d (%d)\n", 52, 64, numSteps);
        *time_first_increase = d_abs_d(currPoint_d.time);
        convert_point_data_d_to_mp(&currPoint_mp, &currPoint_d);
        curr_prec = T->Precision = 64;
        curr_digits = prec_to_digits(curr_prec);
      }
    }
    else
    { // perform the newton corrections in multi precision - only updates after each successful iteration
      retVal = newton_mp_amp(&norm_J, &norm_J_inv, curr_digits, currPoint_mp.point, currPoint_mp.time, T, OUT, &e_mp, ED_mp, eval_func_mp);
      // check for errors
      if (retVal)
      { // need to increase precision - if possible
        if (curr_prec + 32 <= T->AMP_max_prec)
        { // raise precision
          fprintf(OUT, "Increase Prec (refine): %d to %d (%d)\n", curr_prec, curr_prec + 32, numSteps);
          curr_prec = T->Precision = curr_prec + 32;
          curr_digits = prec_to_digits(curr_prec);

          initMP(T->Precision);
          change_prec(ED_mp, T->Precision);
          // raise precision on other MP datatypes
          setprec_point_data_mp(&tempPoint_mp, T->Precision);
          mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
          setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
          setprec_vec_mp(dZ_mp, T->Precision);
          setprec_eval_struct_mp(e_mp, T->Precision);

          // increase precision on the MP datatypes that contain important info
          change_prec_mp2(finalT_mp, T->Precision);
          change_prec_mp2(currStepSize_mp, T->Precision);
          change_prec_point_data_mp(&currPoint_mp, T->Precision);
        }
        else
        { // we are out of precision - exit loop
          retVal = retVal_max_prec_reached;
          cont = 0;
        }
      }
    }

    // check for completion
    if (!retVal)
    { // exit loop since we are refined
      cont = 0;
    }
  }

  // so we are ready to start the tracking if we have not run out of precision
  if (retVal == retVal_max_prec_reached)
    cont = 0; // do not enter loop
  else
    cont = 1; // want to start tracking

  // initialize numSteps & consecSuccess
  numSteps = consecSuccess = 0;

  // main tracking loop
  while (cont)
  { // initialize for loop
    need_to_refine = 0;

    // see whether we are to perform a step in double or multi precision
    if (curr_prec < 64) // do a step in double precision
    { // print out the current time
      if (T->outputLevel > 0)
      {
        fprintf(OUT, "Time = "); print_d(OUT, 16, currPoint_d.time); fprintf(OUT, "\n");
        if (T->screenOut)
        {
          fprintf(stdout, "Time = "); print_d(stdout, 16, currPoint_d.time); fprintf(stdout, "\n");
        }
      }

      // find the distance to the end = end time - current time
      sub_d(distLeft_d, finalT_d, currPoint_d.time);
      absDistLeft_d = d_abs_d(distLeft_d);

      // check to see how close we are to the target
      if (finalT_d->r == currPoint_d.time->r && finalT_d->i == currPoint_d.time->i)
      { // we are exactly at the end - exit loop
        retVal = cont = 0;
        break;
      }
      else if (absDistLeft_d < T->currentStepSize)
      { // we can take a smaller step than our current step size to hit the target time
        set_d(nextTime_d, finalT_d);
        sub_d(dT_d, nextTime_d, currPoint_d.time); // dT = nextTime - currTime
      }
      else
      { // the distance left is larger than the current step size, so we head in the correct direction
        tempD = T->currentStepSize / absDistLeft_d;
        mul_rdouble_d(dT_d, distLeft_d, tempD); // this will put dT = T->currentStepSize * sign(distanceleft)
        add_d(nextTime_d, currPoint_d.time, dT_d); // nextTime = current_time + dT
      }

      // prediction step
      if (T->odePredictor == -1)
      { // use heun method
        retVal = constant_d_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      } 
      else if (T->odePredictor == 1)
      { // use heun method
        retVal = heun_d_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 2)
      { // use standard Runge-Kutta predictor (RK4)
        retVal = runge_kutta_d_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else
      { // use euler predictor
        retVal = euler_d_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }

      // check for prediction errors
      if (retVal)
      { // prediction errors always cause an increase in precision - move from _d to _mp
        fprintf(OUT, "Increase Prec (prediction): %d to %d (%d)\n", 52, 64, numSteps);
        consecSuccess = 0;
        *time_first_increase = d_abs_d(currPoint_d.time);
        convert_point_data_d_to_mp(&currPoint_mp, &currPoint_d);

        // cut the step size since we had an error and then convert to MP
        T->currentStepSize *= T->step_fail_factor;
        mpf_set_d(currStepSize_mp->r, T->currentStepSize);
       
        // increase precision to 64-bits
        curr_prec = T->Precision = 64;
        curr_digits = prec_to_digits(curr_prec);

        // since we increased precision, we will need to refine
        need_to_refine = 1;
      }
      else
      { // check criterion A & C
        if (curr_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
        {
          fprintf(OUT, "NOTE: Criterion A has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }
        else 
        { // find norm in double precision
          tempD = infNormVec_d(currPoint_d.point);
          if (curr_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, tempD))
          {
            fprintf(OUT, "NOTE: Criterion C has been violated.\n");
            retVal = retVal_higher_prec_needed;
          }
        }

        // complete prediction and do newton correction if there are no criterion problems
        if (!retVal)
        { // find tempPoint = currPoint + dZ, tempTime = nextTime
          vec_add_d(tempPoint_d.point, currPoint_d.point, dZ_d);
          set_d(tempPoint_d.time, nextTime_d);
          // do newton correction on tempPoint
          retVal = newton_d_amp(&norm_J, &norm_J_inv, curr_digits, tempPoint_d.point, tempPoint_d.time, T, OUT, &e_d, ED_d, eval_func_d);
        }

        if (retVal == retVal_higher_prec_needed || retVal == retVal_Failed_to_converge)
        { // we have either a criterion problem or newton convergence failure - recover from this error
          if (retVal == retVal_Failed_to_converge)
          { // recover from newton iterations failing to converge properly
            retVal = AMP2_convergence_error_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, &T->currentStepSize, currPoint_d.time, T); 
          }
          else
          { // recover from a criterion problem
            retVal = AMP2_criterion_error_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, &T->currentStepSize, currPoint_d.time, T);
          }
          
          // determine if the precision needs to be increased
          if (retVal > 0)
          { // we need to go to MP
            *time_first_increase = d_abs_d(currPoint_d.time);
            // determine the precision
            retVal = digits_to_prec(retVal);
            if (retVal < 64)
              retVal = 64;

            if (retVal <= T->AMP_max_prec)
            { // setup MP at this precision
              fprintf(OUT, "Increase Prec (newton/criterion): %d to %d (%d)\n", 52, retVal, numSteps);
              curr_prec = T->Precision = retVal;
              curr_digits = prec_to_digits(curr_prec);

              initMP(T->Precision);
              change_prec(ED_mp, T->Precision);

              // raise precision on other MP datatypes
              setprec_point_data_mp(&tempPoint_mp, T->Precision);
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
              setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
              setprec_vec_mp(dZ_mp, T->Precision);
              setprec_eval_struct_mp(e_mp, T->Precision);

              // setup currPoint_mp
              setprec_point_data_mp(&currPoint_mp, T->Precision);
              convert_point_data_d_to_mp(&currPoint_mp, &currPoint_d);

              // setup currStepSize
              setprec_mp(currStepSize_mp, T->Precision);
              mpf_set_d(currStepSize_mp->r, T->currentStepSize);

              // increase precision on the MP datatypes that contain important info
              change_prec_mp2(finalT_mp, T->Precision);

              // since we increased precision, we will need to refine
              need_to_refine = 1;
            }
            else
            { // we are out of precision - exit loop
              retVal = retVal_max_prec_reached;
              cont = 0;
              break;
            }
          }
        }
        else
        { // update currPoint_d with tempPoint_d since this step was successful
          point_data_cp_d(&currPoint_d, &tempPoint_d);

          // check for other errors - e.g. going to infinity
          if (retVal)
          { // immediately exit
            cont = 0;
            break;
          }

          // we had successful step - consider adjusting the step size or precision
          consecSuccess++;
          retVal = AMP2_success_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, &T->currentStepSize, currPoint_d.time, T);

          // determine if precision needs to be increased - hopefully not since this was a successful step!
          if (retVal > 0)
          { // we need to go to MP
            *time_first_increase = d_abs_d(currPoint_d.time);
            // determine the precision
            retVal = digits_to_prec(retVal);
            if (retVal < 64)
              retVal = 64;

            if (retVal <= T->AMP_max_prec)
            { // setup MP at this precision
              fprintf(OUT, "Increase Prec (success): %d to %d (%d)\n", 52, retVal, numSteps);
              curr_prec = T->Precision = retVal;
              curr_digits = prec_to_digits(curr_prec);

              initMP(T->Precision);
              change_prec(ED_mp, T->Precision);
              // raise precision on other MP datatypes
              setprec_point_data_mp(&tempPoint_mp, T->Precision);
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
              setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
              setprec_vec_mp(dZ_mp, T->Precision);
              setprec_eval_struct_mp(e_mp, T->Precision);

              // setup currPoint_mp
              setprec_point_data_mp(&currPoint_mp, T->Precision);
              convert_point_data_d_to_mp(&currPoint_mp, &currPoint_d);

              // setup currStepSize
              setprec_mp(currStepSize_mp, T->Precision);
              mpf_set_d(currStepSize_mp->r, T->currentStepSize);

              // increase precision on the MP datatypes that contain important info
              change_prec_mp2(finalT_mp, T->Precision);

              // since we increased precision, we will need to refine
              need_to_refine = 1;
            }
            else
            { // we are out of precision - exit loop
              retVal = retVal_max_prec_reached;
              cont = 0;
              break;
            }
          }
        }
      }
    }
    else // do a step in multi precision 
    { // print out the current time
      if (T->outputLevel > 0)
      {
        fprintf(OUT, "Time = "); print_mp(OUT, 0, currPoint_mp.time); fprintf(OUT, "\n");
        if (T->screenOut)
        {
          fprintf(stdout, "Time = "); print_mp(stdout, 0, currPoint_mp.time); fprintf(stdout, "\n");
        }
      }

      // find the distance to the end = end time - current time
      sub_mp(distLeft_mp, finalT_mp, currPoint_mp.time);
      mpf_abs_mp(absDistLeft_mp, distLeft_mp);

      // check to see how close we are to the target
      if (mpfr_equal_p(finalT_mp->r, currPoint_mp.time->r) && mpfr_equal_p(finalT_mp->i, currPoint_mp.time->i))
      { // we are exactly at the end!
        retVal = cont = 0;
        break;
      }
      else if (mpfr_less_p(absDistLeft_mp, currStepSize_mp->r))
      { // we can take a smaller step than our current step size to hit the target time
        set_mp(nextTime_mp, finalT_mp);
        sub_mp(dT_mp, nextTime_mp, currPoint_mp.time); // dT = nextTime - currTime
      }
      else
      { // the distance left is larger than the current step size, so we head in the correct direction will a step size of the current step size
        mpf_div(tempMPF, currStepSize_mp->r, absDistLeft_mp);
        mul_rmpf_mp(dT_mp, distLeft_mp, tempMPF); // this will put dT = currentStepSize * sign(distanceleft)
        add_mp(nextTime_mp, currPoint_mp.time, dT_mp); // nextTime = current_time + dT
      }

      // prediction step
      if (T->odePredictor == -1)
      { // use heun method
        retVal = constant_mp_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 1)
      { // use heun method
        retVal = heun_mp_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 2)
      { // use standard Runge-Kutta predictor (RK4)
        retVal = runge_kutta_mp_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else
      { // use euler predictor
        retVal = euler_mp_amp(&norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }

      // check for prediction errors
      if (retVal)
      { // prediction errors always cause an increase in precision
        consecSuccess = 0;
        if (curr_prec + 32 <= T->AMP_max_prec)
        { // raise precision
          fprintf(OUT, "Increase Prec(prediction): %d to %d (%d)\n", curr_prec, curr_prec + 32, numSteps);
          curr_prec = T->Precision = curr_prec + 32;
          curr_digits = prec_to_digits(curr_prec);

          initMP(T->Precision);
          change_prec(ED_mp, T->Precision);
          // raise precision on other MP datatypes
          setprec_point_data_mp(&tempPoint_mp, T->Precision);
          mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
          setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
          setprec_vec_mp(dZ_mp, T->Precision);
          setprec_eval_struct_mp(e_mp, T->Precision);

          // increase precision on the MP datatypes that contain important info
          change_prec_mp2(currStepSize_mp, T->Precision);
          change_prec_mp2(finalT_mp, T->Precision);
          change_prec_point_data_mp2(&currPoint_mp, T->Precision);

          // cut the step size since we had an error
          mpf_set_d(currStepSize_mp->i, T->step_fail_factor);
          mpf_mul(currStepSize_mp->r, currStepSize_mp->r, currStepSize_mp->i);

          // since we increased precision, we will need to refine
          need_to_refine = 1;
        }
        else
        { // we are out of precision - exit loop
          retVal = retVal_max_prec_reached;
          cont = 0;
          break;
        }
      }
      else
      { // check criterion A & C
        if (curr_digits < amp_criterion_A(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi))
        {
          fprintf(OUT, "NOTE: Criterion A has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }
        else
        { // find norm in double precision
          tempD = infNormVec_mp(currPoint_mp.point);
          if (curr_digits < amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, tempD))
          {
            fprintf(OUT, "NOTE: Criterion C has been violated.\n");
            retVal = retVal_higher_prec_needed;
          }
        }

        // complete prediction and do newton correction if there are no criterion problems
        if (!retVal)
        { // find tempPoint = currPoint + dZ, tempTime = nextTime
          vec_add_mp(tempPoint_mp.point, currPoint_mp.point, dZ_mp);
          set_mp(tempPoint_mp.time, nextTime_mp);
          // do newton correction on tempPoint
          retVal = newton_mp_amp(&norm_J, &norm_J_inv, curr_digits, tempPoint_mp.point, tempPoint_mp.time, T, OUT, &e_mp, ED_mp, eval_func_mp);
        }

        if (retVal == retVal_higher_prec_needed || retVal == retVal_Failed_to_converge)
        { // we have either a criterion problem or newton convergence failure - recover from this error
          if (retVal == retVal_Failed_to_converge)
          { // recover from newton iterations failing to converge properly
            retVal = AMP2_convergence_error_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, currStepSize_mp->r, currPoint_mp.time, T);
          }
          else
          { // recover from a criterion problem
            retVal = AMP2_criterion_error_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, currStepSize_mp->r, currPoint_mp.time, T);
          }

          // check to see if AMP2 wants more precision
          if (retVal > curr_digits)
          { // find the new precision
            retVal = digits_to_prec(retVal);
            fprintf(OUT, "Increase Prec (newton/criterion): %d", curr_prec);
            if (retVal <= curr_prec)
            { // we just increase in a 32-bit packet - but this should not happen!
              curr_prec += 32;
            }
            else
            { // increase to retVal
              curr_prec = retVal;
            }
            fprintf(OUT, " to %d (%d)\n", curr_prec, numSteps);
            // verify the new precision is not too large
            if (curr_prec <= T->AMP_max_prec)
            { // raise precision
              T->Precision = curr_prec;
              curr_digits = prec_to_digits(curr_prec);

              initMP(T->Precision);
              change_prec(ED_mp, T->Precision);
              // raise precision on other MP datatypes
              setprec_point_data_mp(&tempPoint_mp, T->Precision);
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
              setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
              setprec_vec_mp(dZ_mp, T->Precision);
              setprec_eval_struct_mp(e_mp, T->Precision);

              // increase precision on the MP datatypes that contain important info
              change_prec_mp2(currStepSize_mp, T->Precision);
              change_prec_mp2(finalT_mp, T->Precision);
              change_prec_point_data_mp2(&currPoint_mp, T->Precision);

              // since we increased precision, we will need to refine
              need_to_refine = 1;
            }
            else
            { // we are out of precision - exit loop
              retVal = retVal_max_prec_reached;
              cont = 0;
              break;
            }
          }
        }
        else
        { // update currPoint_mp with tempPoint_mp since this step was successful
          point_data_cp_mp(&currPoint_mp, &tempPoint_mp);

          // check for other errors - e.g. going to infinity
          if (retVal)
          { // immediately exit
            cont = 0;
            break;
          }

          // we had successful step - consider adjusting the step size or precision
          consecSuccess++;
          retVal = AMP2_success_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, currStepSize_mp->r, currPoint_mp.time, T, &prec_decreases, max_prec_decreases);

          if (retVal > 0)
          { // we need to adjust the precision
            retVal = digits_to_prec(retVal);
            if (retVal < 64)
            { // we can move to double precision
              fprintf(OUT, "Decrease Prec (success): %d to %d (%d)\n", curr_prec, 52, numSteps);
              curr_prec = T->Precision = 52;
              curr_digits = prec_to_digits(curr_prec);

              // copy over current step size
              T->currentStepSize = mpf_get_d(currStepSize_mp->r);
              // copy over currPoint
              convert_point_data_mp_to_d(&currPoint_d, &currPoint_mp);

              // since we changed precision, we will need to refine
              need_to_refine = 1;
            }
            else if (retVal <= T->AMP_max_prec)
            { // setup MP at this precision
              fprintf(OUT, "Change Prec (success): %d to %d (%d)\n", curr_prec, retVal, numSteps);
              curr_prec = T->Precision = retVal;
              curr_digits = prec_to_digits(curr_prec);

              initMP(T->Precision);
              change_prec(ED_mp, T->Precision);
              // change precision on other MP datatypes
              setprec_point_data_mp(&tempPoint_mp, T->Precision);
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision);
              setprec_mp(nextTime_mp, T->Precision); setprec_mp(distLeft_mp, T->Precision); setprec_mp(dT_mp, T->Precision);
              setprec_vec_mp(dZ_mp, T->Precision);
              setprec_eval_struct_mp(e_mp, T->Precision);

              // change precision on the MP datatypes that contain important info
              change_prec_mp2(finalT_mp, T->Precision);
              change_prec_mp2(currStepSize_mp, T->Precision);
              change_prec_point_data_mp2(&currPoint_mp, T->Precision);

              // since we changed precision, we will need to refine
              need_to_refine = 1;
            }
            else
            { // we are out of precision - exit loop
              retVal = retVal_max_prec_reached;
              cont = 0;
              break;
            }
          }
        }
      }
    }

    // see if we need to refine
    if (need_to_refine)
    { // do a basic refine - i.e. do not check AMP conditions
      if (curr_prec == 52) 
      { // refine in double precision
        refine_d_basic(&currPoint_d, T, OUT, &e_d, ED_d, eval_func_d);
      }
      else
      { // refine in multi precision
        refine_mp_basic(&currPoint_mp, T, OUT, &e_mp, ED_mp, eval_func_mp);
      }
    }

    // increment the number of steps and verify we have not taken too many steps
    numSteps++;
    if (numSteps > T->maxNumSteps)
    { // too many steps - exit loop
      retVal = retVal_too_many_steps;
      cont = 0;
      break;
    }
  }

  // copy currPoint to Final and setup prec_out
  *prec_out = curr_prec;
  if (curr_prec < 64)
  { // copy to Final_d
    point_data_cp_d(Final_d, &currPoint_d);
  }
  else
  { // copy to Final_mp - make sure precision is correct
    setprec_point_data_mp(Final_mp, curr_prec);
    point_data_cp_mp(Final_mp, &currPoint_mp);

    // record the step size
    T->currentStepSize = mpf_get_d(currStepSize_mp->r);
  }

  // clear
  clear_vec_d(dZ_d);
  clear_point_data_d(&currPoint_d); clear_point_data_d(&tempPoint_d);
  clear_eval_struct_d(e_d);

  clear_point_data_mp(&currPoint_mp); clear_point_data_mp(&tempPoint_mp);
  mpf_clear(absDistLeft_mp); mpf_clear(tempMPF); mpf_clear(minPrecTime_mp);
  clear_mp(nextTime_mp); clear_mp(distLeft_mp); clear_mp(finalT_mp); clear_mp(dT_mp); clear_mp(currStepSize_mp);
  clear_vec_mp(dZ_mp);
  clear_eval_struct_mp(e_mp);

  // display how many steps it took
  if (T->outputLevel >= 0)
  {
    fprintf(OUT, "numSteps = %d\n", numSteps);
    if (T->screenOut)
      printf("numSteps = %d\n", numSteps);
  }

  return retVal;
}

int constant_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: constant prediction                                    *
\***************************************************************/
{
  int i, retVal = 0;

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // evaluate for output
    eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
    }
  }

  // set the proportional constant - called 'a' in AMP2
  *sizeProportion = 0;

  // set norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 1;

  // setup dZ 
  change_size_vec_d(dZ, P->point->size);
  for (i = 0; i < P->point->size; i++)
    set_zero_d(&dZ->coord[i]); 

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
    fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dZ);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_d(stdout, 10, P->point);
      printf("dX = ");  printVec_d(stdout, 10, dZ);
    }
  }

  return retVal;
}

int euler_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion, dZ         *
* NOTES: finds dZ for the euler step in t of size dT            *
\***************************************************************/
{
  int i, j, rows, cols, new_col, retVal = 0;
  double cond_num;
  comp_d tempComp;

  // Get an estimate via Euler's method.
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "H_of_X = ");  printPoint_d(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_d(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_d(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_d(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_d(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_d(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_d(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_d(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_d(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_d(stdout, 10, e->Jp);
    }
  }

  // find dH/dt = dH/ds * ds/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
  rows = e->Jp->rows; // number of functions
  cols = e->Jp->cols; // number of parameters
  new_col = e->Jv->cols; // number of variables

  // set the sizes of Jv & f
  increase_size_mat_d(e->Jv, rows + 1, new_col + 1);
  increase_size_vec_d(e->funcVals, rows + 1);
  e->funcVals->size = e->Jv->rows = rows + 1;
  e->Jv->cols = new_col + 1;
  for (i = 0; i <= rows; i++)
  {
    if (i < rows)
    { // setup b
      set_zero_d(&e->funcVals->coord[i]);
      // setup the dH/dt entry
      set_zero_d(&e->Jv->entry[i][new_col]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
    }
    else
    { // setup the last entry of b
      set_d(&e->funcVals->coord[i], dT);
      // setup the last row as [0 1]
      for (j = 0; j <= new_col; j++)
      {
        set_zero_d(&e->Jv->entry[i][j]);
      }
      e->Jv->entry[i][new_col].r = 1;
    }
  }

  // do matrixSolve & cond_num together to calculate dX = (Jv)^(-1)*b
  retVal = matrixSolve_cond_num_norms_d(dZ, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  T->steps_since_last_CN = 0;
  T->latest_cond_num_exp = ceil(log10(cond_num));

  // find dZ if matrixSolve was successful
  if (!retVal)
  { // find the proportional constant for the euler prediction - called 'a' in AMP2
    *sizeProportion = 0;
    rows = dZ->size;
    for (i = 0; i < rows; i++)
    { // find |dZ[i]/dT|
      div_d(tempComp, &dZ->coord[i], dT);
      cond_num = d_abs_d(tempComp);
      if (cond_num > *sizeProportion)
        *sizeProportion = cond_num;
    }

    // remove the dT term from bottom of dZ
    dZ->size = dZ->size - 1;

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_d(stdout, 10, P->point);
        printf("dX = ");  printVec_d(stdout, 10, dZ);
      }
    }
  }
  else
  {
    fprintf(OUT, "NOTE: matrixSolve has failed when doing euler prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing euler prediction.\n");
  }

  return retVal;
}

int newton_d_amp(double *norm_J, double *norm_J_inv, int num_digits, point_d P, comp_d t, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does newton corrections on P until convergence and     *
* and checks criterion B & C along the way                      *
* only updates P after all of the error checking is completed   *
\***************************************************************/
{
  int its = 0, retVal = 0;
  double sizeP = 0;

  while (!retVal && its < T->maxNewtonIts)
  { // do a newton iteration
    retVal = newton_residual_d(e->parDer, &T->latest_newton_residual_d, 1, norm_J, norm_J_inv, P, t, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    { // matrixSolve failed - exit with error
      fprintf(OUT, "NOTE: matrixSolve has failed when doing newton correction.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed when doing newton correction.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    { // print out information about newton residual
      if (T->outputLevel > 1)
      {
        fprintf(OUT, "residual = %e\n", T->latest_newton_residual_d);
        if (T->screenOut)
          printf("residual = %e\n", T->latest_newton_residual_d);
      }
      
      // check for convergence and criterion B & C if newton worked
      if (T->latest_newton_residual_d < T->currentNewtonTol)
      { // we have converged - update P and exit loop
        vec_sub_d(P, P, e->parDer);
        retVal = 0;
        its = T->maxNewtonIts;
      }
      else
      { // find the size of P
        sizeP = infNormVec_d(P);

        // check criterion B & C
        if ((T->maxNewtonIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, *norm_J, *norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, T->maxNewtonIts, its))
        {
          fprintf(OUT, "NOTE: Criterion B has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }
        else if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, *norm_J_inv, T->AMP_Psi, T->currentNewtonTol, sizeP))
        {
          fprintf(OUT, "NOTE: Criterion C has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }

        // check for infinite path
        if ((T->endgameSwitch) && (sizeP > T->goingToInfinity))
          retVal = retVal_going_to_infinity;

        // update P if no errors
        if (!retVal)
        {
          vec_sub_d(P, P, e->parDer);
        }

        // increment the number of iterations completed
        its++;
      }
    }
  }

  // check for convergence if there were no other errors
  if (!retVal && T->latest_newton_residual_d > T->currentNewtonTol)
  { // failed to converge, but each iteration was successful
    retVal = retVal_Failed_to_converge;
  }

  return retVal;  
}

int constant_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: constant prediction step                               *
\***************************************************************/
{
  int i, retVal = 0;

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  { // evaluate for output
    eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_mp(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 10, e->Jp);
    }
  }

  // set the proportional constant - called 'a' in AMP2
  *sizeProportion = 0;

  // set norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 1;

  // setup dZ
  change_size_vec_mp(dZ, P->point->size);
  for (i = 0; i < P->point->size; i++)
    set_zero_mp(&dZ->coord[i]);

  // Spit out most computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 10, P->point);
    fprintf(OUT, "dX = ");  printVec_mp(OUT, 10, dZ);
    if (T->screenOut)
    {
      printf("P->point = ");  printPoint_mp(stdout, 10, P->point);
      printf("dX = ");  printVec_mp(stdout, 10, dZ);
    }
  }

  return retVal;
}

int euler_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion,dZ(for AMP2)*
* NOTES: finds dZ for the euler step in t of size dT            *
\***************************************************************/
{
  int i, j, rows, cols, new_col, retVal = 0;
  double cond_num;
  comp_mp tempComp;

  init_mp(tempComp);

  // Get an estimate via Euler's method.
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  // Spit out most some of the computed info if outputLevel >= 3.
  if (T->outputLevel > 2)
  {
    fprintf(OUT, "H_of_X = ");  printPoint_mp(OUT, 10, e->funcVals);
    fprintf(OUT, "newParam = ");  printVec_mp(OUT, 10, e->parVals);
    fprintf(OUT, "parDer = ");  printVec_mp(OUT, 10, e->parDer);
    fprintf(OUT, "Jv = ");  printMat_mp(OUT, 10, e->Jv);
    fprintf(OUT, "Jp = ");  printMat_mp(OUT, 10, e->Jp);
    if (T->screenOut)
    {
      printf("H_of_X = ");  printPoint_mp(stdout, 10, e->funcVals);
      printf("newParam = ");  printVec_mp(stdout, 10, e->parVals);
      printf("parDer = ");  printVec_mp(stdout, 10, e->parDer);
      printf("Jv = ");  printMat_mp(stdout, 10, e->Jv);
      printf("Jp = ");  printMat_mp(stdout, 10, e->Jp);
    }
  }

  // find dH/dt = dH/ds * ds/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
  rows = e->Jp->rows; // number of functions
  cols = e->Jp->cols; // number of parameters
  new_col = e->Jv->cols; // number of variables

  // set the sizes of Jv & f
  increase_size_mat_mp(e->Jv, rows + 1, new_col + 1);
  increase_size_vec_mp(e->funcVals, rows + 1);
  e->funcVals->size = e->Jv->rows = rows + 1;
  e->Jv->cols = new_col + 1;
  for (i = 0; i <= rows; i++)
  {
    if (i < rows)
    { // setup b
      set_zero_mp(&e->funcVals->coord[i]);
      // setup the dH/dt entry
      set_zero_mp(&e->Jv->entry[i][new_col]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }  
    }
    else
    { // setup the last entry of b
      set_mp(&e->funcVals->coord[i], dT);
      // setup the last row as [0 1]
      for (j = 0; j <= new_col; j++)
      {
        set_zero_mp(&e->Jv->entry[i][j]);
      }
      mpf_set_ui(e->Jv->entry[i][new_col].r, 1);
    }
  }

  // do matrixSolve & cond_num together to calculate dX = (Jv)^(-1)*b
  retVal = matrixSolve_cond_num_norms_mp(dZ, e->Jv, e->funcVals, &cond_num, norm_J, norm_J_inv);
  T->steps_since_last_CN = 0;
  T->latest_cond_num_exp = ceil(log10(cond_num));

  // find dZ if matrixSolve was successful
  if (!retVal)
  { // find the proportional constant for the euler prediction - called 'a' in AMP2
    *sizeProportion = 0;
    rows = dZ->size;
    for (i = 0; i < rows; i++)
    { // find |dZ[i]/dT|
      div_mp(tempComp, &dZ->coord[i], dT);
      cond_num = d_abs_mp(tempComp);
      if (cond_num > *sizeProportion)
        *sizeProportion = cond_num;
    }

    // remove the dT term from bottom of dZ
    dZ->size = dZ->size - 1;

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_mp(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_mp(stdout, 10, P->point);
        printf("dX = ");  printVec_mp(stdout, 10, dZ);
      }
    }
  }
  else
  {
    fprintf(OUT, "NOTE: matrixSolve has failed when doing euler prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing euler prediction.\n");
  }

  clear_mp(tempComp);

  return retVal;
}

int newton_mp_amp(double *norm_J, double *norm_J_inv, int num_digits, point_mp P, comp_mp t, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does newton corrections on P until convergence and     *
* and checks criterion B & C along the way                      *
* only updates P after all of the error checking is completed   *
\***************************************************************/
{
  int its = 0, retVal = 0;
  double sizeP = 0;

  while (!retVal && its < T->maxNewtonIts)
  { // do a newton iteration
    retVal = newton_residual_mp(e->parDer, T->latest_newton_residual_mp, 1, norm_J, norm_J_inv, P, t, e, ED, eval_func);
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    { // matrixSolve failed - exit with error
      fprintf(OUT, "NOTE: matrixSolve has failed when doing newton correction.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve has failed when doing newton correction.\n");
      retVal = retVal_higher_prec_needed;
    }
    else
    { // print out information about newton residual
      if (T->outputLevel > 1)
      {
        fprintf(OUT, "residual = %e\n", T->latest_newton_residual_d);
        if (T->screenOut)
          printf("residual = %e\n", T->latest_newton_residual_d);
      }

      // check for convergence and criterion B & C if newton worked
      if (T->latest_newton_residual_d < T->currentNewtonTol)
      { // we have converged - update P and exit loop
        vec_sub_mp(P, P, e->parDer);
        retVal = 0;
        its = T->maxNewtonIts;
      }
      else
      { // find the size of P
        sizeP = infNormVec_mp(P);

        // check criterion B & C
        if ((T->maxNewtonIts != its + 1) && num_digits < amp_criterion_B(T->AMP_safety_digits_1, *norm_J, *norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, T->latest_newton_residual_d, T->maxNewtonIts, its))
        {
          fprintf(OUT, "NOTE: Criterion B has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }
        else if (num_digits < amp_criterion_C(T->AMP_safety_digits_2, *norm_J_inv, T->AMP_Psi, T->currentNewtonTol, sizeP))
        {
          fprintf(OUT, "NOTE: Criterion C has been violated.\n");
          retVal = retVal_higher_prec_needed;
        }

        // check for infinite path
        if ((T->endgameSwitch) && (sizeP > T->goingToInfinity))
          retVal = retVal_going_to_infinity;

        // update P if no errors
        if (!retVal)
        {
          vec_sub_mp(P, P, e->parDer);
        }

        // increment the number of iterations completed
        its++;
      }
    }
  }

  // check for convergence if there were no other errors
  if (!retVal && T->latest_newton_residual_d > T->currentNewtonTol)
  { // failed to converge, but each iteration was successful
    retVal = retVal_Failed_to_converge;
  }

  return retVal;
}

int AMP2_success_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - adjusting step size                       *
*               > 0 - number of digits that it should use       *
* NOTES: looks at the data and decide whether to increase       *
* precision or shrink the step size                             *
\***************************************************************/
// < cSecInc, the only option is to increase precision
// == cSecInc, the only options are to increase step length or increase precision
// > cSecInc - reset back to 0
{
  int retVal, bestDigits, bestPrec, digits_C, currPrec = 52, currDigits = prec_to_digits(52);
  double eta_minSS, eta_maxSS, minStepSize, maxStepSize, P0;
  mpf_t stepSize_mp;

  mpf_init(stepSize_mp);

  // Find the minimum step size
  minStepSize = *currStepSize * T->step_fail_factor;

  // find P0
  P0 = amp_criterion_B2(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, sizeProportion, T->maxNewtonIts);

  // approximate rule C - (take the norm of the current point to be 1) - just need an approximate value rule C
  digits_C = floor(amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, 1) + 0.5);

  if (*consecSuccess < T->cSecInc)
  { // leave the step size alone and only increase precision if necessary
    maxStepSize = *currStepSize;

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }
  else if (*consecSuccess == T->cSecInc)
  { // increase the maximum step size
    maxStepSize = *currStepSize * T->step_success_factor;
    if (maxStepSize > T->maxStepSize)
      maxStepSize = T->maxStepSize;

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }
  else // > T->cSecInc
  { // set back to 0
    *consecSuccess = 0;

    // leave the step size alone and only increase precision if necessary
    maxStepSize = *currStepSize;

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }

  // find eta associated with minSS & maxSS
  eta_minSS = -log10(minStepSize);
  eta_maxSS = -log10(maxStepSize);

  // run AMP2 to see what we should do
  AMP2_update(&bestDigits, &bestPrec, currStepSize, stepSize_mp, currPrec, P0, digits_C, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, currTime, NULL);

  if (bestPrec > currPrec)
  { // we need to update currStepSize since this will be used by calling function
    *currStepSize = mpf_get_d(stepSize_mp);
  }

  // clear MP
  mpf_clear(stepSize_mp);

  // setup retVal
  if (bestPrec != currPrec)
  { // we have decided to change precision - return the number of digits it should use
    retVal = bestDigits;
  }
  else
  { // the current precision is optimal - this only adjusted the step size
    retVal = -1;
  }

  return retVal;
}

int AMP2_convergence_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - cutting step size                         *
*               > 0 - number of digits that it should use       * 
* NOTES: Newton iterations simply failed to converge - decrease *
* step size, if this will be too small, increase precision      *
\***************************************************************/
{
  int retVal = 0, currDigits, currPrec = 52;
  double newStepSize;

  // find the number of digits
  currDigits = prec_to_digits(currPrec);

  // reset consecSuccess
  *consecSuccess = 0;

  newStepSize = *currStepSize * T->step_fail_factor;

  // verify double precision can handle this step size
  if (newStepSize > minTime_d(currDigits, currTime, NULL, currPrec))
  { // we can safely cut the step size
    *currStepSize = newStepSize;
    retVal = -1;
  }
  else
  { // the step size is too small for double precision - increase to more digits
    retVal = currDigits + 1;
    // shrink the step size a little
    *currStepSize = *currStepSize * (1 + T->step_fail_factor) / 2.0;
  }

  return retVal;
}

int AMP2_criterion_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, double *currStepSize, comp_d currTime, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - cutting step size                         *
*               > 0 - number of digits that it should use       *
* NOTES: looks at the data and decide whether to increase       *
* precision or shrink the step size                             *
\***************************************************************/
// There was a criterion problem - look to stay in double precision by cutting the 
// step size, but only if it is optimal
{
  int retVal, bestDigits, bestPrec, currDigits, digits_C, currPrec = 52;
  double eta_minSS, eta_maxSS, minStepSize, maxStepSize, P0;
  mpf_t stepSize_mp;

  mpf_init(stepSize_mp);

  // find the number of digits
  currDigits = prec_to_digits(currPrec);

  // reset consecSuccess
  *consecSuccess = 0;

  // Find the minimum step size possible for double precision
  minStepSize = minTime_d(currDigits, currTime, NULL, currPrec);

  // find the new step size
  maxStepSize = *currStepSize * T->step_fail_factor; // this guarantees that the step size will get small enough at some point to move us to MP

  if (maxStepSize < minStepSize)
  { // we need to increase precision
    bestPrec = 64;
    bestDigits = 18;
    mpf_set_d(stepSize_mp, *currStepSize * (1 + T->step_fail_factor) / 2.0);
  }
  else
  { // double precision still works

    // find eta associated withh minSS & maxSS
    eta_minSS = -log10(minStepSize);
    eta_maxSS = -log10(maxStepSize);

    // find P0
    P0 = amp_criterion_B2(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, sizeProportion, T->maxNewtonIts);

    // approximate rule C - (take the norm of the current point to be 1) - just need an approximate value rule C
    digits_C = floor(amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, 1) + 0.5);

    // run AMP2 to see what we should do
    AMP2_update(&bestDigits, &bestPrec, currStepSize, stepSize_mp, currPrec, P0, digits_C, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, currTime, NULL);
  }

  if (bestPrec > currPrec)
  { // we need to update currStepSize since this will be used by calling function
    *currStepSize = mpf_get_d(stepSize_mp);
  }

  // clear MP
  mpf_clear(stepSize_mp);

  // setup retVal
  if (bestPrec != currPrec)
  { // we have decided to change precision - return the number of digits it should use
    retVal = bestDigits;
  }
  else
  { // the current precision is optimal - this only adjusted the step size
    retVal = -1;
  }

  return retVal;
}

int AMP2_success_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T, int *prec_decreases, int max_prec_decreases)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - adjusting step size                       *
*               > 0 - number of digits that it should use       *
* NOTES: looks at the data and decide whether to adjust the     *
* precision and/or adjust the step size                         *
\***************************************************************/
// < cSecInc, the only option is to increase precision
// == cSecInc, the only options are to increase step length or increase precision
// between cSecInc and 2*cSecInc, the only option is to increase precision
// == 2*cSecInc, all options available - decrease/increase precision, increase step length
// > 2*cSecInc - reset back to 0 and have the only option to be to increase precision
{
  int retVal = 0, bestDigits, bestPrec, digits_C, currDigits = prec_to_digits(curr_prec);
  double eta_minSS, eta_maxSS, P0, stepSize_d;
  mpf_t tempMPF, minStepSize, maxStepSize;

  // initialize MP
  mpf_init(tempMPF); mpf_init(minStepSize); mpf_init(maxStepSize);

  // Find the minimum step size
  mpf_set_d(tempMPF, T->step_fail_factor);
  mpf_mul(minStepSize, currStepSize, tempMPF);

  // find P0
  P0 = amp_criterion_B2(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, sizeProportion, T->maxNewtonIts);

  // approximate rule C - (take the norm of the current point to be 1) - just need an approximate value rule C
  digits_C = floor(amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, 1) + 0.5);
 
  if (*consecSuccess < T->cSecInc)
  { // leave the step size alone and only increase precision if necessary
    mpf_set(maxStepSize, currStepSize);

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }
  else if (*consecSuccess == T->cSecInc)
  { // increase the maximum step size
    mpf_set_d(tempMPF, T->step_success_factor);
    mpf_mul(maxStepSize, currStepSize, tempMPF);
    if (mpf_get_d(maxStepSize) > T->maxStepSize)
      mpf_set_d(maxStepSize, T->maxStepSize);

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }
  else if (*consecSuccess < 2 * T->cSecInc)
  { // leave the step size alone and only increase precision if necessary
    mpf_set(maxStepSize, currStepSize);

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }
  else if (*consecSuccess == 2 * T->cSecInc)
  { // increase the maximum step size
    mpf_set_d(tempMPF, T->step_success_factor);
    mpf_mul(maxStepSize, currStepSize, tempMPF);
    if (mpf_get_d(maxStepSize) > T->maxStepSize)
      mpf_set_d(maxStepSize, T->maxStepSize);

    if (*prec_decreases >= max_prec_decreases)
    { // use digits_C to keep the precision atleast at the current level since we have already exhausted the number of times we can decrease the precision
      // NOTE: This is used to avoid jumping up and down in precision constantly without any progress - if we have this many increases and then decreases,
      //       we should just use the higher precision and track it correctly!
      digits_C = MAX(digits_C, currDigits);
    }
    // else leave digits_C alone - allow for the possibility of decreasing precision
  }
  else // > 2 * T->cSecInc
  { // set back to 0
    *consecSuccess = 0;

    // leave the step size alone and only increase precision if necessary
    mpf_set(maxStepSize, currStepSize);

    // use digits_C to keep the precision atleast at the current level
    digits_C = MAX(digits_C, currDigits);
  }

  // find eta associated withh minSS & maxSS
  mpfr_log10(tempMPF, minStepSize, __gmp_default_rounding_mode);
  eta_minSS = -mpf_get_d(tempMPF); // -log10(minStepSize)
  mpfr_log10(tempMPF, maxStepSize, __gmp_default_rounding_mode);
  eta_maxSS = -mpf_get_d(tempMPF); // -log10(maxStepSize)

  // run AMP2 to see what we should do
  AMP2_update(&bestDigits, &bestPrec, &stepSize_d, currStepSize, curr_prec, P0, digits_C, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, NULL, currTime);

  if (bestPrec < 64)
  { // we need to update currStepSize since this will be used by calling function
    mpf_set_d(currStepSize, stepSize_d);
  }

  if (bestPrec < curr_prec)
  { // update the number of decreases
    *prec_decreases = *prec_decreases + 1;
  }

  // clear MP
  mpf_clear(tempMPF); mpf_clear(minStepSize); mpf_clear(maxStepSize);

  // setup retVal
  if (bestPrec != curr_prec)
  { // we have decided to change precision - return the number of digits it should use
    retVal = bestDigits;
  }
  else
  { // the current precision is optimal - this only (possibly) adjusted the step size
    retVal = -1;
  }

  return retVal;
}

int AMP2_convergence_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - cutting step size                         *
*               > 0 - number of digits that it should use       *
* NOTES: Newton iterations simply failed to converge - decrease *
* step size, if this will be too small, increase precision      *
\***************************************************************/
{
  int retVal = -1, curr_digits = prec_to_digits(curr_prec);
  mpf_t tempMPF, newStepSize;

  mpf_init(tempMPF); mpf_init(newStepSize);

  // reset consecSuccess
  *consecSuccess = 0;

  // find the new step size if we are to stay in this precision
  mpf_set_d(tempMPF, T->step_fail_factor);
  mpf_mul(newStepSize, tempMPF, currStepSize); // currStepSize * step_fail_factor

  // find the minimum step size this precision can handle and make sure the current step size is larger than that
  minTime_mp(tempMPF, curr_digits, NULL, currTime, curr_prec);
  if (mpf_cmp(newStepSize, tempMPF) > 0)
  { // we can safely cut the step size
    mpf_set(currStepSize, newStepSize);
    retVal = -1;
  }
  else
  { // the step size is too small for the current precision - increase to more digits
    retVal = curr_digits + 5;
    // shrink the step size a little
    mpf_set_d(tempMPF, (1 + T->step_fail_factor) / 2.0);
    mpf_mul(currStepSize, tempMPF, currStepSize);
  }

  mpf_clear(tempMPF); mpf_clear(newStepSize);

  return retVal;
}

int AMP2_criterion_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: -1 - cutting step size                         *
*               > 0 - number of digits that it should use       *
* NOTES: looks at the data and decides whether to increase      *
* precision or shrink the step size                             *
\***************************************************************/
// if consecSuccess >= T->cSecInc, the step size was just increased, so decrease the step size if possible
// otherwise, since we are in multi precision, we want to try to increase precision unless the
// precision is significantly larger than criterion 'B2' in which case, we want to
// cut the step size unless the minimum size tells us we must increase precision
{
  int retVal = 0, extra_digits_before_cut_size, digits_C, currDigits = prec_to_digits(curr_prec);
  double eta, P0;
  mpf_t tempMPF, newStepSize;

  // initialize MP
  mpf_init(tempMPF); mpf_init(newStepSize);

  if (*consecSuccess >= T->cSecInc)
  { // the step size was just raised - try cutting the step size first
    extra_digits_before_cut_size = 0; // set favor towards cutting step size rather than increasing precision
  }
  else
  { // set favor towards increasing precision
    extra_digits_before_cut_size = 5;
  }
  // reset consecSuccess
  *consecSuccess = 0;

  // find the new step size if we are to stay in this precision
  mpf_set_d(tempMPF, T->step_fail_factor);
  mpf_mul(newStepSize, tempMPF, currStepSize); // currStepSize * step_fail_factor

  mpfr_log10(tempMPF, newStepSize, __gmp_default_rounding_mode);
  eta = -mpf_get_d(tempMPF); // -log10(newStepSize)
  P0 = amp_criterion_B2(T->AMP_safety_digits_1, norm_J, norm_J_inv, T->AMP_eps, T->AMP_Phi, T->currentNewtonTol, sizeProportion, T->maxNewtonIts);

  digits_C = floor(amp_criterion_C(T->AMP_safety_digits_2, norm_J_inv, T->AMP_Psi, T->currentNewtonTol, 1) + 0.5);

  // see if current precision still works
  if (currDigits + eta / T->maxNewtonIts > P0 + extra_digits_before_cut_size && currDigits > digits_C)
  { // we should cut the step size - verify that it would not be too small
    minTime_mp(tempMPF, currDigits, NULL, currTime, curr_prec);
    if (mpf_cmp(newStepSize, tempMPF) > 0)
    { // we can safely cut the step size
      mpf_set(currStepSize, newStepSize);
      retVal = -1;
    }
    else
    { // the step size is too small for the current precision - increase to more digits
      retVal = currDigits + 5;
    }
  }
  else
  { // the current precision does not suffice
    retVal = currDigits + 5;
  }

  // check to see if the current precision sufficed above
  if (retVal > 0)
  { // increase precision - shrink the step size a little and then determine the number of digits needed
    mpf_set_d(tempMPF, (1 + T->step_fail_factor) / 2.0);
    mpf_mul(currStepSize, tempMPF, currStepSize);

    mpfr_log10(tempMPF, currStepSize, __gmp_default_rounding_mode);
    eta = -mpf_get_d(tempMPF); // -log10(currStepSize)

    retVal = ceil(P0 - eta / T->maxNewtonIts);
    // verify that digits are actually being added
    if (retVal <= currDigits)
      retVal = currDigits + 5;
  }

  mpf_clear(tempMPF); mpf_clear(newStepSize);

  return retVal;
}

int minDigits_from_eta(double eta, comp_d currTime_d, comp_mp currTime_mp, int prec_in)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: minimum number of digits needed based on eta   *
* NOTES:                                                        *
\***************************************************************/
{
  int minDigits;
  double norm_time;

  // find log10(currentTime)
  if (prec_in < 64)
    norm_time = d_abs_d(currTime_d);
  else
    norm_time = d_abs_mp(currTime_mp);
  norm_time = log10(norm_time);

  // the minimum digits is eta + log10(currentTime) + 2 
  minDigits = ceil(eta) + ceil(norm_time) + 2;

  return minDigits;
}

double AMP2_cost(int prec, double eta)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: cost per unit for prec (bits)                  *
* NOTES:                                                        *
\***************************************************************/
{
  double cost, stepSize;

  // find the numerator
  if (prec < 64)
    cost = 1;
  else
    cost = 0.039 * prec + 10.35;

  // find the denominator
  stepSize = pow(10, -eta);

  // divide
  cost = cost / stepSize;

  return cost;
}

void AMP2_minimize_cost(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, double eta_max, double P0, int maxNewtonIts, int minPrecNeeded, int maxPrecNeeded)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: newDigits, newPrec, and updates the step size  *
* NOTES: looks at the data and decide whether to adjust         *
* precision and/or adjust the step size                         *
\***************************************************************/
{
  int i, tempInt, count = 0, *prec = NULL;
  double *eta = NULL, *cost = NULL;

  // error checking
  if (maxPrecNeeded < minPrecNeeded)
  { // set maxPrecNeeded = minPrecNeeded
    maxPrecNeeded = minPrecNeeded;
  }

  // so we have minPrecNeeded <= maxPrecNeeded  
  // count the number of precision values between minPrecNeeded & maxPrecNeeded
  if (maxPrecNeeded == minPrecNeeded)
  { // equal, so only 1 precision is needed
    count = 1;
  }
  else // maxPrecNeeded > minPrecNeeded
  { // set tempInt to the next precision past minPrecNeeded
    if (minPrecNeeded < 64)
      tempInt = 64;
    else
      tempInt = minPrecNeeded + 32;
    // loop until we have went past maxPrecNeeded
    count = 1;
    while (tempInt <= maxPrecNeeded)
    { // increment count
      count++;
      // update tempInt
      tempInt += 32;
    }
  }

  // allocate memory
  prec = (int *)bmalloc(count * sizeof(int));
  eta = (double *)bmalloc(count * sizeof(double));
  cost = (double *)bmalloc(count * sizeof(double));

  // setup the first location
  prec[0] = minPrecNeeded;
  eta[0] = (P0 - prec_to_digits(prec[0])) * maxNewtonIts;
  if (eta[0] < eta_max)
    eta[0] = eta_max;
  cost[0] = AMP2_cost(prec[0], eta[0]);

  // set tempInt to the next precision past minPrecNeeded
  if (minPrecNeeded < 64)
    tempInt = 64;
  else
    tempInt = minPrecNeeded + 32;
  // setup the rest of the structures using tempInt & count
  for (i = 1; i < count; i++)
  { // setup the next location
    prec[i] = tempInt;
    eta[i] = (P0 - prec_to_digits(prec[i])) * maxNewtonIts;
    if (eta[i] < eta_max)
      eta[i] = eta_max;
    cost[i] = AMP2_cost(prec[i], eta[i]);

    // update tempInt
    tempInt += 32;
  }

  // minimize the cost
  tempInt = 0;
  for (i = 1; i < count; i++)
    if (cost[i] < cost[tempInt])
    { // smaller cost - update tempInt
      tempInt = i;
    }

  // setup newDigits & newPrec
  *newPrec = prec[tempInt];
  *newDigits = prec_to_digits(*newPrec);

  // setup step size
  if (*newPrec < 64)
  { // setup currStepSize_d
    *currStepSize_d = pow(10, -eta[tempInt]);
  }
  else
  { // setup currStepSize_mp
    mpf_set_prec(currStepSize_mp, *newPrec);
    mpf_set_d(currStepSize_mp, -eta[tempInt]);
    mpfr_ui_pow(currStepSize_mp, 10, currStepSize_mp, __gmp_default_rounding_mode);
  }

  // free the memory
  free(prec);
  free(eta);
  free(cost);

  return;
}

void AMP2_update(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, int currPrec, double P0, int digits_C, double newtonTol, int maxNewtonIts, double finalTol, int checkFinalTol, double eta_minSS, double eta_maxSS, comp_d currTime_d, comp_mp currTime_mp)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: newDigits, newPrec, and updates the step size  *
* NOTES: looks at the data and decide whether to adjust         *
* precision and/or adjust the step size                         *
\***************************************************************/
{
  int minPrecNeeded, maxPrecNeeded, minDigitsNeeded, maxDigitsNeeded, digits_tau, digits_final, digits_B2, digits_stepSize; 
  int tempInt, extra_criterion_digits_before_decrease_precision, currDigits = prec_to_digits(currPrec);

  // find the number of extra digits we would need before signaling a decrease in precision
  if (currPrec < 64)
  { // already in double precision so this does not matter
    extra_criterion_digits_before_decrease_precision = 0;
  }
  else
  { // work out a formula based on the current precision
    extra_criterion_digits_before_decrease_precision = (currPrec / 32.0) * 2;
  }

  // find the minimum digits based on the newton tolerance
  digits_tau = ceil(-log10(newtonTol));

  // find the minimum digits based on the final tolerance
  if (checkFinalTol)
    digits_final = ceil(-log10(finalTol)) + 1;
  else
    digits_final = 0;

  // find the minimum digits based on rule B2
  digits_B2 = ceil(P0 - eta_minSS / maxNewtonIts);
  // check to see if B2 suggests a decrease in precision
  if (digits_to_prec(digits_B2) < currPrec)
  { // B2 suggets a cut - so we make sure that we have at least digits_B2+extra digits to actually suggest to decrease the precision
    // otherwise, leave the number of digits at the current precision
    tempInt = digits_B2 + extra_criterion_digits_before_decrease_precision;
    digits_B2 = MAX(digits_B2, tempInt);
    // do not go more than currDigits since digits_B2 originally says we can cut
    digits_B2 = MIN(digits_B2, currDigits);
  }

  // check to see if C suggests a decrease in precision
  if (digits_to_prec(digits_C) < currPrec)
  { // C suggets a cut - so we make sure that we have at least digits_C+extra digits to actually suggest to decrease the precision
    // otherwise, leave the number of digits at the current precision
    tempInt = digits_C + extra_criterion_digits_before_decrease_precision;
    digits_C = MAX(digits_C, tempInt);
    // do not go more than currDigits since digits_C originally says we can cut
    digits_C = MIN(digits_C, currDigits);
  }

  // find the minimum digits to handle the step sizes
  digits_stepSize = minDigits_from_eta(eta_minSS, currTime_d, currTime_mp, currPrec);
  tempInt = minDigits_from_eta(eta_maxSS, currTime_d, currTime_mp, currPrec);
  digits_stepSize = MAX(digits_stepSize, tempInt);

  // check to see if the step size digits suggest a decrease in precision
  if (digits_to_prec(digits_stepSize) < currPrec)
  { // step size digits suggest a cut - so we make sure that we have at least digits+extra digits to actually suggest to decrease precision
    // otherwise, leave the number of digits at the current precision
    tempInt = digits_stepSize + extra_criterion_digits_before_decrease_precision;
    digits_stepSize = MAX(digits_stepSize, tempInt);
    // do not go more than currDigits since digits originally says we can cut
    digits_stepSize = MIN(digits_stepSize, currDigits);
  }

  // find the minimum number of digits that satisfies all of the previous conditions
  minDigitsNeeded = MAX(digits_tau, digits_final);
  minDigitsNeeded = MAX(minDigitsNeeded, digits_B2);
  minDigitsNeeded = MAX(minDigitsNeeded, digits_C);
  minDigitsNeeded = MAX(minDigitsNeeded, digits_stepSize);

  // find the precision associated with minDigitsNeeded and then find the number of digits associated with that precision
  minPrecNeeded = digits_to_prec(minDigitsNeeded);
  minDigitsNeeded = prec_to_digits(minPrecNeeded);

  // find the number of digits that is needed for the largest step size
  tempInt = ceil(P0 - eta_maxSS / maxNewtonIts);

  // find the maximum number of digits that would possibly be needed
  maxDigitsNeeded = MAX(minDigitsNeeded, tempInt);

  // find the precision associated with maxDigitsNeeded and then find the number of digits associated with that precision
  maxPrecNeeded = digits_to_prec(maxDigitsNeeded);
  maxDigitsNeeded = prec_to_digits(maxPrecNeeded);

  // find the best number of digits, precision and the current step size to minimize cost per unit based on minPrecNeeded & maxPrecNeeded
  AMP2_minimize_cost(newDigits, newPrec, currStepSize_d, currStepSize_mp, eta_maxSS, P0, maxNewtonIts, minPrecNeeded, maxPrecNeeded);

  return;
}

int heun_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion, dZ         *
* NOTES: finds dZ for the heun step in t of size dT             *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0;
  double cond_num, norm, norm_inv;
  comp_d newTime;
  point_d dX[2];

  for (i = 0; i < 2; i++)
    init_point_d(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 2; k++)
  { // find dH/dt = dH/dp * dp/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters
    new_col = e->Jv->cols; // number of variables

    // set the sizes of Jv & f
    increase_size_mat_d(e->Jv, rows + 1, new_col + 1);
    increase_size_vec_d(e->funcVals, rows + 1);
    e->funcVals->size = e->Jv->rows = rows + 1;
    e->Jv->cols = new_col + 1;
    for (i = 0; i <= rows; i++)
    {
      if (i < rows)
      { // setup b
        set_zero_d(&e->funcVals->coord[i]);
        // setup the dH/dt entry
        set_zero_d(&e->Jv->entry[i][new_col]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_d(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
        }
      }
      else
      { // setup the last entry of b
        set_d(&e->funcVals->coord[i], dT);
        // setup the last row as [0 1]
        for (j = 0; j <= new_col; j++)
        {
          set_zero_d(&e->Jv->entry[i][j]);
        }
        e->Jv->entry[i][new_col].r = 1;
      }
    }

    // do matrixSolve & cond_num together to calculate dX = Jv^-1 * b
    retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, e->funcVals, &cond_num, &norm, &norm_inv);
    T->steps_since_last_CN = 0;
    T->latest_cond_num_exp = ceil(log10(cond_num));

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // update norm_J & norm_J_inv
      if (*norm_J < norm)
        *norm_J = norm;
      if (*norm_J_inv < norm_inv)
        *norm_J_inv = norm_inv;

      // determine if more function evaluations are needed
      if (k == 0)
      { // find the next point: x + dX[0]
        increase_size_vec_d(dX[1], P->point->size);
        rows = dX[1]->size = P->point->size;
        // find next point
        for (i = 0; i < rows; i++)
        {
          add_d(&dX[1]->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
        }
        // find new time = t + dT
        add_d(newTime, P->time, dT);

        // do the next function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[1], newTime, ED, eval_func);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/2 (dX[0] + dX[1]) and the proportional constant for the prediction - called 'a' in AMP2
    *sizeProportion = 0;
    increase_size_vec_d(dZ, P->point->size);
    rows = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // adjust a
      for (j = 0; j < 2; j++)
      {
        cond_num = dX[j]->coord[i].r * dX[j]->coord[i].r + dX[j]->coord[i].i * dX[j]->coord[i].i;
        if (cond_num > *sizeProportion)
          *sizeProportion = cond_num;
      }

      // find dZ[i]
      add_d(&dZ->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
      mul_rdouble_d(&dZ->coord[i], &dZ->coord[i], 0.5);
    }
    // take sqrt(a) and then divide by |dT|
    *sizeProportion = sqrt(*sizeProportion);
    *sizeProportion /= d_abs_d(dT);

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_d(stdout, 10, P->point);
        printf("dX = ");  printVec_d(stdout, 10, dZ);
      }
    }
  }
  else
  { // display error message
    fprintf(OUT, "NOTE: matrixSolve has failed when doing heun prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing heun prediction.\n");
  }

  for (i = 0; i < 2; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int runge_kutta_d_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion, dZ         *
* NOTES: finds dZ for the Runge-Kutta step in t of size dT      *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0;
  double cond_num, norm, norm_inv;
  comp_d newTime;
  point_d dX[4];

  for (i = 0; i < 4; i++)
    init_point_d(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 4; k++)
  { // find dH/dt = dH/dp * dp/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters 
    new_col = e->Jv->cols; // number of variables

    // set the sizes of Jv & f
    increase_size_mat_d(e->Jv, rows + 1, new_col + 1);
    increase_size_vec_d(e->funcVals, rows + 1);
    e->funcVals->size = e->Jv->rows = rows + 1;
    e->Jv->cols = new_col + 1;
    for (i = 0; i <= rows; i++)
    {
      if (i < rows)
      { // setup b
        set_zero_d(&e->funcVals->coord[i]);
        // setup the dH/dt entry
        set_zero_d(&e->Jv->entry[i][new_col]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_d(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
        }
      }
      else
      { // setup the last entry of b
        set_d(&e->funcVals->coord[i], dT);
        // setup the last row as [0 1]
        for (j = 0; j <= new_col; j++)
        {
          set_zero_d(&e->Jv->entry[i][j]);
        }
        e->Jv->entry[i][new_col].r = 1;
      }
    }

    // do matrixSolve & cond_num together to calculate dX = Jv^-1 * b
    retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, e->funcVals, &cond_num, &norm, &norm_inv);
    T->steps_since_last_CN = 0;
    T->latest_cond_num_exp = ceil(log10(cond_num));

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // update norm_J & norm_J_inv
      if (*norm_J < norm)
        *norm_J = norm;
      if (*norm_J_inv < norm_inv)
        *norm_J_inv = norm_inv;

      // determine if more function evaluations are needed
      if (k < 3)
      { // find the next point: x + a * dX[k], where a = 1/2 when k == 0 || k == 1 and a = 1 when k == 2
        increase_size_vec_d(dX[3], P->point->size);
        rows = dX[3]->size = P->point->size;
        if (k < 2)
        { // find next point
          for (i = 0; i < rows; i++)
          {
            mul_rdouble_d(&dX[3]->coord[i], &dX[k]->coord[i], 0.5);
            add_d(&dX[3]->coord[i], &dX[3]->coord[i], &P->point->coord[i]);
          }
          // find new time = t + a * dT
          mul_rdouble_d(newTime, dT, 0.5);
          add_d(newTime, newTime, P->time);
        }
        else // k == 2
        { // find next point
          for (i = 0; i < rows; i++)
          {
            add_d(&dX[3]->coord[i], &dX[k]->coord[i], &P->point->coord[i]);
          }
          // find new time = t + dT
          add_d(newTime, P->time, dT);
        }
        // do the next function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[3], newTime, ED, eval_func);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/6 (dX[0] + 2(dX[1] + dX[2]) + dX[3]) and the proportional constant for the prediction - called 'a' in AMP2
    *sizeProportion = 0;
    increase_size_vec_d(dZ, P->point->size);
    rows = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // adjust a
      for (j = 0; j < 4; j++)
      {
        cond_num = dX[j]->coord[i].r * dX[j]->coord[i].r + dX[j]->coord[i].i * dX[j]->coord[i].i;
        if (cond_num > *sizeProportion)
          *sizeProportion = cond_num;
      }

      // find dZ[i]
      add_d(&dZ->coord[i], &dX[0]->coord[i], &dX[3]->coord[i]);
      add_d(&dX[1]->coord[i], &dX[1]->coord[i], &dX[2]->coord[i]);
      mul_rdouble_d(&dX[1]->coord[i], &dX[1]->coord[i], 2);
      add_d(&dZ->coord[i], &dZ->coord[i], &dX[1]->coord[i]);
      mul_rdouble_d(&dZ->coord[i], &dZ->coord[i], 1.0 / 6.0);
    }
    // take sqrt(a) and then divide by |dT|
    *sizeProportion = sqrt(*sizeProportion);
    *sizeProportion /= d_abs_d(dT);

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_d(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_d(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_d(stdout, 10, P->point);
        printf("dX = ");  printVec_d(stdout, 10, dZ);
      }
    }
  }
  else
  { // display error message
    fprintf(OUT, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
  }

  for (i = 0; i < 4; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int heun_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion,dZ(for AMP2)*
* NOTES: finds dZ for the heun step in t of size dT             *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0;
  double cond_num, norm, norm_inv;
  mpf_t tempMPF;
  comp_mp newTime;
  point_mp dX[2];

  // initialize MP
  mpf_init(tempMPF);
  init_mp(newTime);
  for (i = 0; i < 2; i++)
    init_point_mp(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 2; k++)
  { // find dH/dt = dH/dp * dp/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters
    new_col = e->Jv->cols; // number of variables

    // set the sizes of Jv & f
    increase_size_mat_mp(e->Jv, rows + 1, new_col + 1);
    increase_size_vec_mp(e->funcVals, rows + 1);
    e->funcVals->size = e->Jv->rows = rows + 1;
    e->Jv->cols = new_col + 1;
    for (i = 0; i <= rows; i++)
    {
      if (i < rows)
      { // setup b
        set_zero_mp(&e->funcVals->coord[i]);
        // setup the dH/dt entry
        set_zero_mp(&e->Jv->entry[i][new_col]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_mp(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
        }
      }
      else
      { // setup the last entry of b
        set_mp(&e->funcVals->coord[i], dT);
        // setup the last row as [0 1]
        for (j = 0; j <= new_col; j++)
        {
          set_zero_mp(&e->Jv->entry[i][j]);
        }
        mpf_set_ui(e->Jv->entry[i][new_col].r, 1);
      }
    }

    // do matrixSolve & cond_num together to calculate dX = Jv^-1 * b
    retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, e->funcVals, &cond_num, &norm, &norm_inv);
    T->steps_since_last_CN = 0;
    T->latest_cond_num_exp = ceil(log10(cond_num));

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // update norm_J & norm_J_inv
      if (*norm_J < norm)
        *norm_J = norm;
      if (*norm_J_inv < norm_inv)
        *norm_J_inv = norm_inv;

      // determine if more function evaluations are needed
      if (k == 0)
      { // find the next point: x + dX[0]
        increase_size_vec_mp(dX[1], P->point->size);
        rows = dX[1]->size = P->point->size;
        // find next point
        for (i = 0; i < rows; i++)
        {
          add_mp(&dX[1]->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
        }
        // find new time = t + dT
        add_mp(newTime, P->time, dT);

        // do the next function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[1], newTime, ED, eval_func);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/2 (dX[0] + dX[1]) and the proportional constant for the prediction - called 'a' in AMP2
    mpf_set_ui(tempMPF, 0);
    increase_size_vec_mp(dZ, P->point->size);
    rows = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // adjust a
      for (j = 0; j < 2; j++)
      {
        mpf_mul(newTime->r, dX[j]->coord[i].r, dX[j]->coord[i].r);
        mpf_mul(newTime->i, dX[j]->coord[i].i, dX[j]->coord[i].i);
        mpf_add(newTime->r, newTime->r, newTime->i);
        if (mpf_cmp(newTime->r, tempMPF) > 0)
          mpf_set(tempMPF, newTime->r);
      }

      // find dZ[i]
      add_mp(&dZ->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
      mpf_div_ui(dZ->coord[i].r, dZ->coord[i].r, 2);
      mpf_div_ui(dZ->coord[i].i, dZ->coord[i].i, 2);
    }
    // take sqrt(a) and then divide by |dT|
    mpf_sqrt(tempMPF, tempMPF);
    mpf_abs_mp(newTime->r, dT);
    mpf_div(tempMPF, tempMPF, newTime->r);
    *sizeProportion = mpf_get_d(tempMPF);

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_mp(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_mp(stdout, 10, P->point);
        printf("dX = ");  printVec_mp(stdout, 10, dZ);
      }
    }
  }
  else
  { // display error message
    fprintf(OUT, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
  }

  // clear MP
  mpf_clear(tempMPF);
  clear_mp(newTime);
  for (i = 0; i < 2; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int runge_kutta_mp_amp(double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    * 
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion,dZ(for AMP2)*
* NOTES: finds dZ for the runge-kutta step in t of size dT      *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0;
  double cond_num, norm, norm_inv;
  mpf_t tempMPF;
  comp_mp newTime;
  point_mp dX[4];

  // initialize MP
  mpf_init(tempMPF);
  init_mp(newTime);
  for (i = 0; i < 4; i++)
    init_point_mp(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 4; k++)
  { // find dH/dt = dH/dp * dp/dt and setup Jv = [[Jv dH/dt][0 1]] and b = [[0][dT]]
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters
    new_col = e->Jv->cols; // number of variables

    // set the sizes of Jv & f
    increase_size_mat_mp(e->Jv, rows + 1, new_col + 1);
    increase_size_vec_mp(e->funcVals, rows + 1);
    e->funcVals->size = e->Jv->rows = rows + 1;
    e->Jv->cols = new_col + 1;
    for (i = 0; i <= rows; i++)
    {
      if (i < rows)
      { // setup b
        set_zero_mp(&e->funcVals->coord[i]);
        // setup the dH/dt entry
        set_zero_mp(&e->Jv->entry[i][new_col]);
        for (j = 0; j < cols; j++)
        {
          sum_mul_mp(&e->Jv->entry[i][new_col], &e->Jp->entry[i][j], &e->parDer->coord[j]);
        }
      }
      else
      { // setup the last entry of b
        set_mp(&e->funcVals->coord[i], dT);
        // setup the last row as [0 1]
        for (j = 0; j <= new_col; j++)
        {
          set_zero_mp(&e->Jv->entry[i][j]);
        }
        mpf_set_ui(e->Jv->entry[i][new_col].r, 1);
      }
    }

    // do matrixSolve & cond_num together to calculate dX = Jv^-1 * b
    retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, e->funcVals, &cond_num, &norm, &norm_inv);
    T->steps_since_last_CN = 0;
    T->latest_cond_num_exp = ceil(log10(cond_num));

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // update norm_J & norm_J_inv
      if (*norm_J < norm)
        *norm_J = norm;
      if (*norm_J_inv < norm_inv)
        *norm_J_inv = norm_inv;

      // determine if more function evaluations are needed
      if (k < 3)
      { // find the next point: x + a * dX[k], where a = 1/2 when k == 0 || k == 1 and a = 1 when k == 2
        increase_size_vec_mp(dX[3], P->point->size);
        rows = dX[3]->size = P->point->size;
        if (k < 2)
        { // find next point
          for (i = 0; i < rows; i++)
          {
            mpf_div_ui(dX[3]->coord[i].r, dX[k]->coord[i].r, 2);
            mpf_div_ui(dX[3]->coord[i].i, dX[k]->coord[i].i, 2);
            add_mp(&dX[3]->coord[i], &dX[3]->coord[i], &P->point->coord[i]);
          }
          // find new time = t + a * dT
          mpf_div_ui(newTime->r, dT->r, 2);
          mpf_div_ui(newTime->i, dT->i, 2);
          add_mp(newTime, newTime, P->time);
        }
        else // k == 2
        { // find next point
          for (i = 0; i < rows; i++)
          {
            add_mp(&dX[3]->coord[i], &dX[k]->coord[i], &P->point->coord[i]);
          }
          // find new time = t + dT
          add_mp(newTime, P->time, dT);
        }
        // do the next function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[3], newTime, ED, eval_func);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/6 (dX[0] + 2(dX[1] + dX[2]) + dX[3]) and the proportional constant for the prediction - called 'a' in AMP2
    mpf_set_ui(tempMPF, 0);
    increase_size_vec_mp(dZ, P->point->size);
    rows = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // adjust a
      for (j = 0; j < 4; j++)
      {
        mpf_mul(newTime->r, dX[j]->coord[i].r, dX[j]->coord[i].r);
        mpf_mul(newTime->i, dX[j]->coord[i].i, dX[j]->coord[i].i);
        mpf_add(newTime->r, newTime->r, newTime->i);
        if (mpf_cmp(newTime->r, tempMPF) > 0)
          mpf_set(tempMPF, newTime->r);
      }

      // find dZ[i]
      add_mp(&dZ->coord[i], &dX[0]->coord[i], &dX[3]->coord[i]);
      add_mp(&dX[1]->coord[i], &dX[1]->coord[i], &dX[2]->coord[i]);
      mpf_mul_ui(dX[1]->coord[i].r, dX[1]->coord[i].r, 2);
      mpf_mul_ui(dX[1]->coord[i].i, dX[1]->coord[i].i, 2);
      add_mp(&dZ->coord[i], &dZ->coord[i], &dX[1]->coord[i]);
      mpf_div_ui(dZ->coord[i].r, dZ->coord[i].r, 6);
      mpf_div_ui(dZ->coord[i].i, dZ->coord[i].i, 6);
    }
    // take sqrt(a) and then divide by |dT|
    mpf_sqrt(tempMPF, tempMPF);
    mpf_abs_mp(newTime->r, dT);
    mpf_div(tempMPF, tempMPF, newTime->r);
    *sizeProportion = mpf_get_d(tempMPF);

    // Spit out most computed info if outputLevel >= 3.
    if (T->outputLevel > 2)
    {
      fprintf(OUT, "P->point = ");  printPoint_mp(OUT, 10, P->point);
      fprintf(OUT, "dX = ");  printVec_mp(OUT, 10, dZ);
      if (T->screenOut)
      {
        printf("P->point = ");  printPoint_mp(stdout, 10, P->point);
        printf("dX = ");  printVec_mp(stdout, 10, dZ);
      }
    }
  }
  else
  { // display error message
    fprintf(OUT, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing runge-kutta prediction.\n");
  }

  // clear MP
  mpf_clear(tempMPF);
  clear_mp(newTime);
  for (i = 0; i < 4; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int refine_d_basic(point_data_d *P, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: basic refine in double precision - only update P if    *
* we are better than when we started                            *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, maxIts = 15;
  int correct_digits_orig = 0, correct_digits_old = 0, correct_digits = 0, correct_digits_stop = prec_to_digits(52) - 4;
  int consecDigits = 0, consecStop = 2, min_digits_for_consecStop = floor(0.75 * correct_digits_stop);
  int min_digits_for_update = floor(-log10(T->currentNewtonTol) + 0.5);

  point_data_d P_temp;
  point_d best_approx;
  int best_approx_digits = 0;

  init_point_d(best_approx, 0);

  // copy P to P_temp
  init_point_data_d(&P_temp, P->point->size);
  point_data_cp_d(&P_temp, P);

  while (cont)
  { // do a newton iteration
    retVal = newton_iteration_d(&T->latest_newton_residual_d, 0, NULL, NULL, P_temp.point, P_temp.time, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in refine.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve failed in refine.");

      clear_point_data_d(&P_temp);
      clear_point_d(best_approx);

      return retVal;
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "refining residual = %.7e\n", T->latest_newton_residual_d);
      if (T->screenOut)
        printf("refining residual = %.7e\n", T->latest_newton_residual_d);
    }

    // find the number of digits correct on this iteration
    if (T->latest_newton_residual_d == 0)
    { // set to the digits_stop
      correct_digits = correct_digits_stop;
    }
    else
    { // calculate the number of digits correct
      correct_digits = floor(-log10(T->latest_newton_residual_d) + 0.5);
    }

    if (its == 0)
    { // initialize
      best_approx_digits = correct_digits_orig = correct_digits_old = correct_digits;
      point_cp_d(best_approx, P_temp.point);
    }
    else
    { // compare against best_approx
      if (correct_digits >= best_approx_digits)
      { // update best_approx
        point_cp_d(best_approx, P_temp.point);
        best_approx_digits = correct_digits;
      }

      // check for leveling off
      if (correct_digits == correct_digits_old && correct_digits >= min_digits_for_consecStop)
      { // increment consecDigits
        consecDigits++;
      }
      else
      { // reset consecDigits
        consecDigits = 0;
      }
    }

    // check to see if we should stop
    if (correct_digits >= correct_digits_stop || consecDigits >= consecStop || correct_digits + 1 < correct_digits_old)
    { // we should stop
      cont = 0;
    }
    else
    { // update digits_old
      correct_digits_old = correct_digits;

      its++;
      // check to see if we have done too many iterations
      if (its > maxIts)
      {
        cont = 0;
      }
    }
  }

  // check to see if we should update P
  if (best_approx_digits >= correct_digits_stop || T->latest_newton_residual_d < T->final_tolerance || (best_approx_digits >= correct_digits_orig && best_approx_digits >= min_digits_for_update))
  { // update P since it is either converged or better than when we started
    point_cp_d(P->point, best_approx);
  }

  // check for convergence - as long as it is better than either the correction tolerance (based on precision) or the final tolerance, we can call refine a success
  if (best_approx_digits >= correct_digits_stop || T->latest_newton_residual_d < T->final_tolerance)
  { // success
    retVal = 0;
  }
  else
  { // failure
    retVal = retVal_refining_failed;
  }

  clear_point_data_d(&P_temp);
  clear_point_d(best_approx);

  return retVal;
}

int refine_mp_basic(point_data_mp *P, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: basic refine in double precision - only update P if    *
* we are better than when we started                            *
\***************************************************************/
{
  int its = 0, cont = 1, retVal = 0, maxIts = 15, best_approx_digits = 0;
  int correct_digits_orig = 0, correct_digits_old = 0, correct_digits = 0, correct_digits_stop = prec_to_digits(T->Precision) - 4;
  int consecDigits = 0, consecStop = 2, min_digits_for_consecStop = floor(0.75 * correct_digits_stop);
  int min_digits_for_update = floor(-log10(T->currentNewtonTol) + 0.5);
  point_data_mp P_temp;
  point_mp best_approx; 

  init_point_mp(best_approx, 0);

  // copy P to P_temp
  init_point_data_mp(&P_temp, P->point->size);
  point_data_cp_mp(&P_temp, P);

  while (cont)
  { // do a newton iteration
    retVal = newton_iteration_mp(T->latest_newton_residual_mp, 0, NULL, NULL, P_temp.point, P_temp.time, e, ED, eval_func);

    // check for newton success (i.e. if matrixSolve failed)
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in refine.\n");
      if (T->screenOut)
        fprintf(stderr, "NOTE: matrixSolve failed in refine.");

      // clear MP
      clear_point_data_mp(&P_temp);
      clear_point_mp(best_approx);

      return retVal;
    }

    if (T->outputLevel > 1)
    {
      fprintf(OUT, "refining residual = "); mpf_out_str(OUT, 10, 7, T->latest_newton_residual_mp); fprintf(OUT, "\n");
      if (T->screenOut)
      {
        printf("refining residual = "); mpf_out_str(stdout, 10, 7, T->latest_newton_residual_mp); printf("\n");
      }
    }

    // find double approximation
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);

    // find the number of digits correct on this iteration
    if (T->latest_newton_residual_d == 0)
    { // set to digits_stop
      correct_digits = correct_digits_stop;
    }
    else
    { // calculate the number of digits correct
      correct_digits = floor(-log10(T->latest_newton_residual_d) + 0.5);
    }

    if (its == 0)
    { // initialize
      best_approx_digits = correct_digits_orig = correct_digits_old = correct_digits;
      point_cp_mp(best_approx, P_temp.point);
    }
    else
    { // compare against best_approx
      if (correct_digits >= best_approx_digits)
      { // update best_approx
        point_cp_mp(best_approx, P_temp.point);
        best_approx_digits = correct_digits;
      }

      // check for leveling off
      if (correct_digits == correct_digits_old && correct_digits >= min_digits_for_consecStop)
      { // increment consecDigits 
        consecDigits++;
      }
      else
      { // reset consecDigits
        consecDigits = 0;
      }
    }

    // check to see if we should stop
    if (correct_digits >= correct_digits_stop || consecDigits >= consecStop || correct_digits + 1 < correct_digits_old)
    { // we should stop
      cont = 0;
    }
    else
    { // update digits_old
      correct_digits_old = correct_digits;

      its++;
      // check to see if we have done too many iterations
      if (its > maxIts)
      {
        cont = 0;
      }
    }
  }

  // check to see if we should update P
  if (best_approx_digits >= correct_digits_stop || T->latest_newton_residual_d < T->final_tolerance || (best_approx_digits >= correct_digits_orig && best_approx_digits >= min_digits_for_update))
  { // update P since it is either converged or better than when we started
    point_cp_mp(P->point, best_approx);
  }

  // check for convergence - as long as it is better than either the correction tolerance (based on precision) or the final tolerance, we can call refine a success
  if (best_approx_digits >= correct_digits_stop || T->latest_newton_residual_d < T->final_tolerance)
  { // success
    retVal = 0;
  }
  else
  { // failure
    retVal = retVal_refining_failed;
  }

  // clear MP
  clear_point_data_mp(&P_temp);
  clear_point_mp(best_approx);

  return retVal;
}

int refine_digits_amp(int outputLevel, int digitsCorrect, double *latest_newton_residual_d, mpf_t latest_newton_residual_mp, int digitsCorrect_in, point_data_d *out_d, point_data_mp *out_mp, int *prec_out, point_data_d *in_d, point_data_mp *in_mp, int prec_in, comp_d time, FILE *OUT, eval_struct_d *e_d, eval_struct_mp *e_mp, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: success: 0, otherwise: retVal_refining_failed  *
* NOTES: refines 'in' so that 'out' has 'digitsCorrect'         *
* uses time in double precision so that we have a consistent    *
* way to increase precision                                     *
\***************************************************************/
{
  int retVal = 0, its = 0, cont = 0, correct_digits = 0, prec_digits = 0, curr_prec = 0, maxIts = 2 * log10(digitsCorrect) / log10(2) + 20;
  point_data_mp tempPoint;
  mpf_t tempMPF;

  // find curr_prec (atleast 64-bits)
  curr_prec = MAX(prec_in, 64);

  // initialize tempMPF
  mpf_init2(tempMPF, curr_prec);

  // initialize tempPoint
  init_point_data_mp2(&tempPoint, 0, curr_prec);
  if (prec_in < 64)
  { // copy in_d to tempPoint
    convert_point_data_d_to_mp(&tempPoint, in_d);

    // setup time
    d_to_mp(tempPoint.time, time);
  }
  else
  { // copy in_mp to tempPoint
    point_data_cp_mp(&tempPoint, in_mp);

    // setup time
    d_to_mp(tempPoint.time, time);
  }

  // setup correct_digits
  correct_digits = MAX(digitsCorrect_in, 1); // just an idea of how many digits should be correct upon input

  // main loop
  do 
  { // find the precision needed for the next iteration - first calculate the '32-bit' multiplier
    prec_digits = (int) ceil((correct_digits + 1.5) / (32 * log10(2)));
    if (prec_digits < 1)
      prec_digits = 64; // need precision to be atleast 64-bits
    else
      prec_digits *= 4 * 32;
    // find the precision that we should use
    if (prec_digits > 3 * curr_prec + 32)
      prec_digits = 3 * curr_prec + 32;

    // make sure that we are atleast increasing the precision every third iteration - this makes sure that precision is not the reason for not refining properly
    if (its > 0 && !(its % 3) && curr_prec >= prec_digits)
      prec_digits = curr_prec + 32;

    // see if we have enough preciion
    if (curr_prec < prec_digits)
    { // precision needs increase
      curr_prec = prec_digits;

      // set everything to this precision
      initMP(curr_prec);
      change_prec(ED_mp, curr_prec);
      setprec_eval_struct_mp(*e_mp, curr_prec);
      mpf_set_prec(latest_newton_residual_mp, curr_prec);
      mpf_set_prec(tempMPF, curr_prec);

      // change tempPoint
      change_prec_point_mp(tempPoint.point, curr_prec);
      setprec_mp(tempPoint.time, curr_prec);
      d_to_mp(tempPoint.time, time);
    }

    // perform a newton iteration
    retVal = newton_iteration_mp(latest_newton_residual_mp, 0, NULL, NULL, tempPoint.point, tempPoint.time, e_mp, ED_mp, eval_func_mp);

    // check for newton success
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in refine, but Bertini can still continue.\n");
      retVal = retVal_refining_failed;
      cont = 0; // need to quit
    }
    else
    { // print residual
      if (outputLevel > 1)
      {
        fprintf(OUT, "refining residual = "); 
        mpf_out_str(OUT, 10, 7, latest_newton_residual_mp);
        fprintf(OUT, "\n");
      }

      // find the number of correct digits
      correct_digits = residual_to_digits_mp(latest_newton_residual_mp, curr_prec);

      // check for convergence
      if (correct_digits >= digitsCorrect)
      { // we have the correct number of digits!
        retVal = cont = 0;
      }
      else
      { // update the number of iteration and check to see if we should continue
        its++;
 
        if (its > maxIts)
        { // too many iterations - need to quit
          retVal = retVal_refining_failed;
          cont = 0;
        }
        else
        { // continue on
          cont = 1;
        }
      }
    }
  } while (cont);

  // calculate the exact precision required based on digitsCorrect
  if (prec_in < 64 && digitsCorrect < prec_to_digits(52))
  { // can use double precision
    prec_digits = 52;
  }
  else
  { // need to use MP
    prec_digits = (int) ceil((digitsCorrect + 0.5) / (32 * log10(2)));
    if (prec_digits < 2)
    { // use 64-bit precision
      prec_digits = 64;
    }
    else
      prec_digits *= 32;

    // make sure we have atleast prec_in
    prec_digits = MAX(prec_digits, prec_in);
  }

  // record this precision
  *prec_out = prec_digits;

  // see what precision to set everything to
  prec_digits = MAX(*prec_out, 64);

  // set back to this precision
  initMP(prec_digits);
  change_prec(ED_mp, prec_digits);
  setprec_eval_struct_mp(*e_mp, prec_digits);
  
  // set latest_newton_residual_mp to this precision
  mpf_set(tempMPF, latest_newton_residual_mp);
  mpf_set_prec(latest_newton_residual_mp, prec_digits);
  mpf_set(latest_newton_residual_mp, tempMPF);

  // setup the things that will be returned
  if (*prec_out < 64)
  { // setup _d
    convert_point_data_mp_to_d(out_d, &tempPoint);
    *latest_newton_residual_d = mpf_get_d(latest_newton_residual_mp);
  }
  else
  { // setup _mp
    setprec_point_mp(out_mp->point, *prec_out);
    setprec_mp(out_mp->time, *prec_out);

    point_data_cp_mp(out_mp, &tempPoint);
  }

  // clear MP
  clear_point_data_mp(&tempPoint);
  mpf_clear(tempMPF); 

  return retVal;
}

int refine_digits_d(int outputLevel, int digitsCorrect, double *latest_newton_residual_d, int digitsCorrect_in, point_data_d *out, point_data_d *in, comp_d time, FILE *OUT, eval_struct_d *e_d, void const *ED_d, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: success: 0, otherwise: retVal_refining_failed  *
* NOTES: refines 'in' so that 'out' has 'digitsCorrect'         *
\***************************************************************/
{
  int retVal = 0, its = 0, cont = 0, correct_digits = 0, maxIts = 2 * log10(digitsCorrect) / log10(2) + 20;
  point_data_d tempPoint;

  // setup tempPoint
  init_point_data_d(&tempPoint, in->point->size);
  point_data_cp_d(&tempPoint, in);

  // setup correct_digits
  correct_digits = MAX(digitsCorrect_in, 1); // just an idea of how many digits should be correct upon input

  do
  { // perform a newton iteration
    retVal = newton_iteration_d(latest_newton_residual_d, 0, NULL, NULL, tempPoint.point, time, e_d, ED_d, eval_func_d);

    // check for newton success
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in refine, but Bertini can still continue.\n");
      retVal = retVal_refining_failed;
      cont = 0; // need to quit
    }
    else
    { // print residual
      if (outputLevel > 1)
        fprintf(OUT, "refining residual = %.7e\n", *latest_newton_residual_d);

      // find the number of correct digits
      correct_digits = residual_to_digits_d(*latest_newton_residual_d);

      // check for convergence
      if (correct_digits >= digitsCorrect)
      { // we have the correct number of digits!
        retVal = cont = 0;
      }
      else
      { // update the number of iteration and check to see if we should continue
        its++;

        if (its > maxIts)
        { // too many iterations - need to quit
          retVal = retVal_refining_failed;
          cont = 0;
        }
        else
        { // continue on
          cont = 1;
        }
      }
    }
  } while (cont);

  // copy tempPoint to out
  point_data_cp_d(out, &tempPoint);
  clear_point_data_d(&tempPoint);

  return retVal;
}

int refine_digits_mp(int outputLevel, int digitsCorrect, mpf_t latest_newton_residual_mp, int digitsCorrect_in, point_data_mp *out, point_data_mp *in, int prec_in, comp_mp time, FILE *OUT, eval_struct_mp *e_mp, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: success: 0, otherwise: retVal_refining_failed  *
* NOTES: refines 'in' so that 'out' has 'digitsCorrect'         *
\***************************************************************/
{
  int retVal = 0, its = 0, cont = 0, correct_digits = 0, maxIts = 2 * log10(digitsCorrect) / log10(2) + 20;
  point_data_mp tempPoint;

  // setup tempPoint
  init_point_data_mp(&tempPoint, in->point->size);
  point_data_cp_mp(&tempPoint, in);

  // setup correct_digits
  correct_digits = MAX(digitsCorrect_in, 1); // just an idea of how many digits should be correct upon input

  do
  { // perform a newton iteration
    retVal = newton_iteration_mp(latest_newton_residual_mp, 0, NULL, NULL, tempPoint.point, time, e_mp, ED_mp, eval_func_mp);

    // check for newton success
    if (retVal)
    {
      fprintf(OUT, "NOTE: matrixSolve failed in refine, but Bertini can still continue.\n");
      retVal = retVal_refining_failed;
      cont = 0; // need to quit
    }
    else
    { // print residual
      if (outputLevel > 1)
      {
        fprintf(OUT, "refining residual = ");
        mpf_out_str(OUT, 10, 7, latest_newton_residual_mp);
        fprintf(OUT, "\n");
      }

      // find the number of correct digits
      correct_digits = residual_to_digits_mp(latest_newton_residual_mp, prec_in);

      // check for convergence
      if (correct_digits >= digitsCorrect)
      { // we have the correct number of digits!
        retVal = cont = 0;
      }
      else
      { // update the number of iteration and check to see if we should continue
        its++;

        if (its > maxIts)
        { // too many iterations - need to quit
          retVal = retVal_refining_failed;
          cont = 0;
        }
        else
        { // continue on
          cont = 1;
        }
      }
    }
  } while (cont);

  // copy tempPoint to out
  point_data_cp_mp(out, &tempPoint);
  clear_point_data_mp(&tempPoint);

  return retVal;
}


