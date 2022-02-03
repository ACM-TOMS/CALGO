// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

// this file implements the AMP tracker with error approximation that is described in the adaptive precision paper
// of Bates, Hauenstein & Sommese

int newton_error_d_amp(int endgameSwitch, int consecSuccess, double errorEstimate, double *norm_J, double *norm_J_inv, int num_digits, point_d P, comp_d t, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int newton_error_mp_amp(int endgameSwitch, int consecSuccess, mpf_t errorEstimate, double *norm_J, double *norm_J_inv, int num_digits, point_mp P, comp_mp t, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int heun_euler_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int heun_euler_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rkn34_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkn34_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rk34_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[5], double a_minus_b[5], double c[5], double d[5][4], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk34_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[5], mpf_t a_minus_b[5], mpf_t c[5], mpf_t d[5][4], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rkf45_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkf45_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rkck45_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkck45_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rk45_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[6], double a_minus_b[6], double c[6], double d[6][5], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk45_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[6], mpf_t a_minus_b[6], mpf_t c[6], mpf_t d[6][5], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rkdp56_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkdp56_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rk56_methods_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[8], double b[8], double c[8], double d[8][7], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk56_methods_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[8], mpf_t b[8], mpf_t c[8], mpf_t d[8][7], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rkv67_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rkv67_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int rk67_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[10], double a_minus_b[10], double c[10], double d[10][9], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *));
int rk67_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[10], mpf_t a_minus_b[10], mpf_t c[10], mpf_t d[10][9], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *));

int AMP3_success_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, int order, double *currStepSize, comp_d currTime, tracker_config_t *T);
int AMP3_success_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, int order, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T, int *prec_decreases, int max_prec_decreases);

int AMP3_criterion_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, int order, double *currStepSize, comp_d currTime, tracker_config_t *T);
int AMP3_criterion_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, int order, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T);

void AMP3_update(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, int currPrec, double P0, int digits_C, int order, double newtonTol, int maxNewtonIts, double finalTol, int checkFinalTol, double eta_minSS, double eta_maxSS, comp_d currTime_d, comp_mp currTime_mp);

int AMPtrack_error(point_data_d *Final_d, point_data_mp *Final_mp, int *prec_out, double *time_first_increase, point_data_d *Start_d, point_data_mp *Start_mp, int prec_in, comp_d fT_d, comp_mp fT_mp, int time_prec, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
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
  int cont, tempInt, retVal = 0, curr_prec = prec_in, numSteps = 0, consecSuccess = 0, curr_digits = prec_to_digits(prec_in), need_to_refine = 0, odeOrder = 0;
  int prec_decreases = 0, max_prec_decreases = 10;
  double absDistLeft_d, tempD, sizeProportion, norm_J = 0, norm_J_inv = 0, residual_old = 0, predictorError_d = 0;
  mpf_t absDistLeft_mp, tempMPF, minPrecTime_mp, predictorError_mp;
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

  // determine the order of the predictor that we are using for error control
  if (T->odePredictor == 3)
  { // use Euler-Heun method
    odeOrder = 1;
  }
  else if (T->odePredictor == 4)
  { // use Runge-Kutta-Norsett (RKN34)
    odeOrder = 3;
  }
  else if (T->odePredictor == 5)
  { // use Runge-Kuta-Fehlberg (RKF45)
    odeOrder = 4;
  }
  else if (T->odePredictor == 6)
  { // use Runge-Kuta-Cash-Karp (RKCK45)
    odeOrder = 4;
  }
  else if (T->odePredictor == 7)
  { // use Runge-Kutta-Dormand-Prince (RKDP56)
    odeOrder = 5;
  }
  else 
  { // use Runge-Kutta-Verner (RKV67)
    odeOrder = 6;
  }

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
  mpf_init(absDistLeft_mp); mpf_init(tempMPF); mpf_init(minPrecTime_mp); mpf_init(predictorError_mp);
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
    mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
          mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
      if (T->odePredictor == 3)
      { // use Euler-Heun method
        retVal = heun_euler_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 4)
      { // use Runge-Kutta-Norsett (RKN34)
        retVal = rkn34_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 5)
      { // use Runge-Kutta-Fehlberg (RKF45)
        retVal = rkf45_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 6)
      { // use Runge-Kutta-Cash-Karp (RKCK45)
        retVal = rkck45_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 7)
      { // use Runge-Kutta-Dormand-Prince (RKDP56)
        retVal = rkdp56_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
      }
      else if (T->odePredictor == 8)
      { // use Runge-Kutta-Verner (RKV67)
        retVal = rkv67_d_amp(&predictorError_d, &norm_J, &norm_J_inv, &sizeProportion, dZ_d, &currPoint_d, T, OUT, dT_d, &e_d, ED_d, eval_func_d);
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
          retVal = newton_error_d_amp(T->endgameSwitch, consecSuccess, predictorError_d, &norm_J, &norm_J_inv, curr_digits, tempPoint_d.point, tempPoint_d.time, T, OUT, &e_d, ED_d, eval_func_d);
        }

        if (retVal == retVal_higher_prec_needed || retVal == retVal_Failed_to_converge)
        { // we have either a criterion problem or newton convergence failure - recover from this error
          if (retVal == retVal_Failed_to_converge)
          { // recover from newton iterations failing to converge properly
            retVal = AMP2_convergence_error_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, &T->currentStepSize, currPoint_d.time, T);
          }
          else
          { // recover from a criterion problem
            retVal = AMP3_criterion_error_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, odeOrder, &T->currentStepSize, currPoint_d.time, T);
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
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
          retVal = AMP3_success_d(&consecSuccess, norm_J, norm_J_inv, sizeProportion, odeOrder, &T->currentStepSize, currPoint_d.time, T);

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
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
      if (T->odePredictor == 3)
      { // use Euler-Heun method
        retVal = heun_euler_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 4)
      { // use Runge-Kutta-Norsett (RKN34)
        retVal = rkn34_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 5)
      { // use Runge-Kutta-Fehlberg (RKF45)
        retVal = rkf45_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 6)
      { // use Runge-Kutta-Cash-Karp (RKCK45)
        retVal = rkck45_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else if (T->odePredictor == 7) 
      { // use Runge-Kutta-Dormand-Prince (RKDP56)
        retVal = rkdp56_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
      }
      else
      { // use Runge-Kutta-Verner (RKV67)
        retVal = rkv67_mp_amp(predictorError_mp, &norm_J, &norm_J_inv, &sizeProportion, dZ_mp, &currPoint_mp, T, OUT, dT_mp, &e_mp, ED_mp, eval_func_mp);
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
          mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
          retVal = newton_error_mp_amp(T->endgameSwitch, consecSuccess, predictorError_mp, &norm_J, &norm_J_inv, curr_digits, tempPoint_mp.point, tempPoint_mp.time, T, OUT, &e_mp, ED_mp, eval_func_mp);
        }

        if (retVal == retVal_higher_prec_needed || retVal == retVal_Failed_to_converge)
        { // we have either a criterion problem or newton convergence failure - recover from this error
          if (retVal == retVal_Failed_to_converge)
          { // recover from newton iterations failing to converge properly
            retVal = AMP2_convergence_error_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, currStepSize_mp->r, currPoint_mp.time, T);
          }
          else
          { // recover from a criterion problem
            retVal = AMP3_criterion_error_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, odeOrder, currStepSize_mp->r, currPoint_mp.time, T);
          }

          // check to see if AMP3 wants more precision
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
              mpf_set_prec(absDistLeft_mp, T->Precision); mpf_set_prec(tempMPF, T->Precision); mpf_set_prec(minPrecTime_mp, T->Precision); mpf_set_prec(predictorError_mp, T->Precision);
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
          retVal = AMP3_success_mp(&consecSuccess, curr_prec, norm_J, norm_J_inv, sizeProportion, odeOrder, currStepSize_mp->r, currPoint_mp.time, T, &prec_decreases, max_prec_decreases);

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
  mpf_clear(absDistLeft_mp); mpf_clear(tempMPF); mpf_clear(minPrecTime_mp); mpf_clear(predictorError_mp);
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

int newton_error_d_amp(int endgameSwitch, int consecSuccess, double errorEstimate, double *norm_J, double *norm_J_inv, int num_digits, point_d P, comp_d t, tracker_config_t *T, FILE *OUT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does newton corrections on P until convergence and     *
* and checks criterion B & C along the way                      *
* only updates P after all of the error checking is completed   *
\***************************************************************/
{
  int its = 0, retVal = 0, cont = 1;
  double sizeP = 0;

  // determine if newton iterations are needed
  if (endgameSwitch || errorEstimate >= T->currentNewtonTol || consecSuccess >= T->cSecInc - 1 || consecSuccess >= T->maxStepsBeforeNewton)
  { // run newton iterations
    cont = 1;
  } 
  else
  { // apporimation was good enough
    T->latest_newton_residual_d = errorEstimate;
    cont = 0;
  }

  while (cont && !retVal && its < T->maxNewtonIts)
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

int newton_error_mp_amp(int endgameSwitch, int consecSuccess, mpf_t errorEstimate, double *norm_J, double *norm_J_inv, int num_digits, point_mp P, comp_mp t, tracker_config_t *T, FILE *OUT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does newton corrections on P until convergence and     *
* and checks criterion B & C along the way                      *
* only updates P after all of the error checking is completed   *
\***************************************************************/
{
  int its = 0, retVal = 0, cont = 1;
  double sizeP = 0;

  // determine if newton iterations are needed
  if (endgameSwitch || mpf_get_d(errorEstimate) >= T->currentNewtonTol || consecSuccess >= T->cSecInc - 1 || consecSuccess >= T->maxStepsBeforeNewton)
  { // run newton iterations
    cont = 1;
  }
  else
  { // apporimation was good enough
    mpf_set(T->latest_newton_residual_mp, errorEstimate);
    T->latest_newton_residual_d = mpf_get_d(T->latest_newton_residual_mp);
    cont = 0;
  }

  while (cont && !retVal && its < T->maxNewtonIts)
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

int heun_euler_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion, dZ         *
* NOTES: finds dZ for the heun step in t of size dT             *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0, digits = prec_to_digits(52) - 3;
  double cond_num;
  comp_d newTime;
  point_d Y, dX[2];

  init_point_d(Y, 0);
  for (i = 0; i < 2; i++)
    init_point_d(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 2; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt 
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters
    new_col = e->Jv->cols; // number of variables
    increase_size_vec_d(Y, rows);
    Y->size = rows;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    { // do matrixSolve & cond_num together to calculate dX = Jv^-1 * Y
      retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // do matrixSolve to calculate dX = Jv^-1 * Y
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k == 0)
      { // compute dX[0] = dX[0] * dT
        // find the next point: x + dX[0]
        increase_size_vec_d(dX[1], P->point->size);
        rows = dX[1]->size = P->point->size;
        // find next point
        for (i = 0; i < rows; i++)
        {
          mul_d(&dX[0]->coord[i], &dX[0]->coord[i], dT);
          add_d(&dX[1]->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
        }
        // find new time = t + dT
        add_d(newTime, P->time, dT);

        // do the next function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[1], newTime, ED, eval_func);
      }
      else
      { // compute dX[1] = dX[1] * dT
        vec_mulcomp_d(dX[1], dX[1], dT);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/2 (dX[0] + dX[1]), error = || 1/2 (dX[0] - dX[1]) || = sizeProportion * || dT ||^2
    *sizeProportion = 0;
    increase_size_vec_d(dZ, P->point->size);
    increase_size_vec_d(Y, P->point->size);
    rows = dZ->size = Y->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // find Y[i] & dX[i]
      sub_d(&Y->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
      add_d(&dZ->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);

      mul_rdouble_d(&Y->coord[i], &Y->coord[i], 0.5);
      mul_rdouble_d(&dZ->coord[i], &dZ->coord[i], 0.5);
    }
    // compute the norm of Y
    *error = infNormVec_d(Y);

    // compute sizeProportion
    cond_num = norm_sqr_d(dT);
    if (*error == 0 || -log10(*error) >= digits || -log10(cond_num) >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional size
      *sizeProportion = 1;
    }
    else
    { // compute the proportional size
      *sizeProportion = pow(10, log10(*error) - log10(cond_num));
    }

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

  clear_point_d(Y);
  for (i = 0; i < 2; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int heun_euler_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: norm_J, norm_J_inv, sizeProportion, dZ         *
* NOTES: finds dZ for the heun step in t of size dT             *
\***************************************************************/
{
  int i, j, k, rows, cols, new_col, retVal = 0, digits = prec_to_digits(T->Precision) - 5;
  double cond_num, logError;
  comp_mp tempComp;
  point_mp Y, dX[2];

  init_mp(tempComp);
  init_point_mp(Y, 0);
  for (i = 0; i < 2; i++)
    init_point_mp(dX[i], 0);

  // initialize norm_J & norm_J_inv
  *norm_J = *norm_J_inv = 0;

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 2; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    rows = e->Jp->rows; // number of functions
    cols = e->Jp->cols; // number of parameters
    new_col = e->Jv->cols; // number of variables
    increase_size_vec_mp(Y, rows);
    Y->size = rows;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    { // do matrixSolve & cond_num together to calculate dX = Jv^-1 * Y
      retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // do matrixSolve to calculate dX = Jv^-1 * Y
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }

    // check for matrixSolve failures
    if (retVal)
    { // error occured - exit loop
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k == 0)
      { // compute dX[0] = dX[0] * dT
        // find the next point: x + dX[0]
        increase_size_vec_mp(dX[1], P->point->size);
        rows = dX[1]->size = P->point->size;
        // find next point
        for (i = 0; i < rows; i++)
        {
          mul_mp(&dX[0]->coord[i], &dX[0]->coord[i], dT);
          add_mp(&dX[1]->coord[i], &dX[0]->coord[i], &P->point->coord[i]);
        }
        // find new time = t + dT
        add_mp(tempComp, P->time, dT);

        // do the next function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, dX[1], tempComp, ED, eval_func);
      }
      else
      { // compute dX[1] = dX[1] * dT
        vec_mulcomp_mp(dX[1], dX[1], dT);
      }
    }
  }

  // make sure no errors occured
  if (!retVal)
  { // find dZ = 1/2 (dX[0] + dX[1]), error = || 1/2 (dX[0] - dX[1]) || = sizeProportion * || dT ||^2
    *sizeProportion = 0;
    increase_size_vec_mp(dZ, P->point->size);
    increase_size_vec_mp(Y, P->point->size);
    rows = dZ->size = Y->size = P->point->size;
    for (i = 0; i < rows; i++)
    { // find Y[i] & dX[i]
      sub_mp(&Y->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);
      add_mp(&dZ->coord[i], &dX[0]->coord[i], &dX[1]->coord[i]);

      div_ui_mp(&Y->coord[i], &Y->coord[i], 2);
      div_ui_mp(&dZ->coord[i], &dZ->coord[i], 2);
    }
    // compute the norm of Y
    infNormVec_mp2(error, Y);
    mpfr_log10(tempComp->r, error, __gmp_default_rounding_mode);
    logError = mpf_get_d(tempComp->r);

    // compute sizeProportion
    norm_sqr_mp(tempComp->r, dT);
    mpfr_log10(tempComp->i, tempComp->r, __gmp_default_rounding_mode);
    cond_num = mpf_get_d(tempComp->i);
    if (mpfr_zero_p(error) || -logError >= digits || -cond_num >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional error
      *sizeProportion = 1;
    }
    else
    {
      mpf_div(tempComp->i, error, tempComp->r);
      *sizeProportion = mpf_get_d(tempComp->i);
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing heun prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing heun prediction.\n");
  }

  clear_mp(tempComp);
  clear_point_mp(Y);
  for (i = 0; i < 2; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkn34_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKN34 method - Norsett method                 *
\***************************************************************/
{
  int rV;
  double a[5] = {25.0/162, 32.0/135, 256.0/567, 0, 11.0/70};
  double a_minus_b[5] = {-41.0/4050, -244.0/1755, 256.0/567, -448.0/975, 11.0/70};
  double c[5] = {0, 3.0/8, 9.0/16, 25.0/32, 1};
  double d[5][4] = {{0, 0, 0, 0},
                    {3.0/8, 0, 0, 0},
                    {0, 9.0/16, 0, 0},
                    {-125.0/672, 325.0/336, 0, 0},
                    {371.0/891, -200.0/297, 1120.0/891, 0}};

  // run the Runge-Kutta34 method
  rV = rk34_methods_d_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkn34_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKN34 method - Norsett method                 *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[5] = {25, 32, 256, 0, 11};
  int a_denom[5] = {162, 135, 567, 1, 70};
  int a_minus_b_numer[5] = {-41, -244, 256, -448, 11};
  int a_minus_b_denom[5] = {4050, 1755, 567, 975, 70};
  int c_numer[5] = {0, 3, 9, 25, 1};
  int c_denom[5] = {1, 8, 16, 32, 1};
  int d_numer[5][4] = {{0, 0, 0, 0},
                       {3, 0, 0, 0},
                       {0, 9, 0, 0},
                       {-125, 325, 0, 0},
                       {371, -200, 1120, 0}};
  int d_denom[5][4] = {{1, 1, 1, 1},
                       {8, 1, 1, 1},
                       {1, 16, 1, 1},
                       {672, 336, 1, 1},
                       {891, 297, 891, 1}};

  mpf_t a[5], a_minus_b[5], c[5], d[5][4];
  for (i = 0; i < 5; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 4; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta34 method
  rV = rk34_methods_mp_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 5; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 4; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk34_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[5], double a_minus_b[5], double c[5], double d[5][4], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK34 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(52) - 3;
  double cond_num;
  comp_d newTime;
  point_d Y, dX[5];

  init_point_d(Y, 0);
  for (i = 0; i < 5; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 5; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 4)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_d(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_d(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_d(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rdouble_d(newTime, dT, c[k+1]);
        add_d(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
        vec_mulcomp_d(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^4
    // find dZ = sum(k = 0..5, a_k * dX[k])
    *sizeProportion = 0;
    increase_size_vec_d(Y, P->point->size);
    increase_size_vec_d(dZ, P->point->size);
    rows = Y->size = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    {
      mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rdouble_d(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 5; k++)
      {
        mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_d(&Y->coord[i], &Y->coord[i], newTime);

        mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
        add_d(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    *error = infNormVec_d(Y);

    // compute sizeProportion
    cond_num = pow(norm_sqr_d(dT), 2);
    if (*error == 0 || -log10(*error) >= digits || -log10(cond_num) >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional size
      *sizeProportion = 1;
    }
    else
    { // compute the proportional size
      *sizeProportion = pow(10, log10(*error) - log10(cond_num));
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk34 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk34 prediction.\n");
  }

  clear_point_d(Y);
  for (i = 0; i < 5; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rk34_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[5], mpf_t a_minus_b[5], mpf_t c[5], mpf_t d[5][4], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK34 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(T->Precision) - 5;
  double cond_num, logError;
  comp_mp newTime;
  point_mp Y, dX[5];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 5; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 5; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 4)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_mp(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_mp(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_mp(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rmpf_mp(newTime, dT, c[k+1]);
        add_mp(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
       vec_mulcomp_mp(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^4
    // find dZ = sum(k = 0..5, a_k * dX[k])
    increase_size_vec_mp(Y, dX[0]->size);
    increase_size_vec_mp(dZ, dX[0]->size);
    rows = Y->size = dZ->size = dX[0]->size;
    for (i = 0; i < rows; i++)
    {
      mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rmpf_mp(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 5; k++)
      {
        mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_mp(&Y->coord[i], &Y->coord[i], newTime);

        mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
        add_mp(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    infNormVec_mp2(error, Y);
    mpfr_log10(newTime->r, error, __gmp_default_rounding_mode);
    logError = mpf_get_d(newTime->r);

    // compute sizeProportion
    abs_mp(newTime, dT);
    mpf_pow_ui(newTime->i, newTime->r, 4);
    mpfr_log10(newTime->r, newTime->i, __gmp_default_rounding_mode);
    cond_num = mpf_get_d(newTime->r);
    if (mpfr_zero_p(error) || -logError >= digits || -cond_num >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional error
      *sizeProportion = 1;
    }
    else
    {
      mpf_div(newTime->r, error, newTime->i);
      *sizeProportion = mpf_get_d(newTime->r);
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk34 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk34 prediction.\n");
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 5; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkf45_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKF45 method - Kincaid/Cheney pg 544-545      *
\***************************************************************/
{
  int rV;
  double a[6] = {16.0/135, 0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55};
  double a_minus_b[6] = {1.0/360, 0, -128.0/4275, -2197.0/75240, 1.0/50, 2.0/55};
  double c[6] = {0, 1.0/4, 3.0/8, 12.0/13, 1, 1.0/2};
  double d[6][5] = {{0, 0, 0, 0, 0},
                    {1.0/4, 0, 0, 0, 0},
                    {3.0/32, 9.0/32, 0, 0, 0},
                    {1932.0/2197, -7200.0/2197, 7296.0/2197, 0, 0},
                    {439.0/216, -8, 3680.0/513, -845.0/4104, 0},
                    {-8.0/27, 2, -3544.0/2565, 1859.0/4104, -11.0/40}};

  // run the Runge-Kutta45 method
  rV = rk45_methods_d_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkf45_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKF45 method - Kincaid/Cheney pg 544-545      *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[6] = {16, 0, 6656, 28561, -9, 2};
  int a_denom[6] = {135, 1, 12825, 56430, 50, 55};
  int a_minus_b_numer[6] = {1, 0, -128, -2197, 1, 2};
  int a_minus_b_denom[6] = {360, 1, 4275, 75240, 50, 55};
  int c_numer[6] = {0, 1, 3, 12, 1, 1};
  int c_denom[6] = {1, 4, 8, 13, 1, 2};
  int d_numer[6][5] = {{0, 0, 0, 0, 0},
                       {1, 0, 0, 0, 0},
                       {3, 9, 0, 0, 0},
                       {1932, -7200, 7296, 0, 0},
                       {439, -8, 3680, -845, 0},
                       {-8, 2, -3544, 1859, -11}};
  int d_denom[6][5] = {{1, 1, 1, 1, 1},
                       {4, 1, 1, 1, 1},
                       {32, 32, 1, 1, 1},
                       {2197, 2197, 2197, 1, 1},
                       {216, 1, 513, 4104, 1},
                       {27, 1, 2565, 4104, 40}};

  mpf_t a[6], a_minus_b[6], c[6], d[6][5];
  for (i = 0; i < 6; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 5; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta45 method
  rV = rk45_methods_mp_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 6; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 5; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rkck45_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKCK45 method - Cash/Karp method              *
\***************************************************************/
{
  int rV;
  double a[6] = {37.0/378, 0, 250.0/621, 125.0/594, 0, 512.0/1771};
  double a_minus_b[6] = {-277.0/64512, 0, 6925.0/370944, -6925.0/202752, -277.0/14336, 277.0/7084};
  double c[6] = {0, 1.0/5, 3.0/10, 3.0/5, 1, 7.0/8};
  double d[6][5] = {{0, 0, 0, 0, 0},
                    {1.0/5, 0, 0, 0, 0},
                    {3.0/40, 9.0/40, 0, 0, 0},
                    {3.0/10, -9.0/10, 6.0/5, 0, 0},
                    {-11.0/54, 5.0/2, -70.0/27, 35.0/27, 0},
                    {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096}};

  // run the Runge-Kutta45 method
  rV = rk45_methods_d_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkck45_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKCK45 method - Cash/Karp method              *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[6] = {37, 0, 250, 125, 0, 512};
  int a_denom[6] = {378, 1, 621, 594, 1, 1771};
  int a_minus_b_numer[6] = {-277, 0, 6925, -6925, -277, 277};
  int a_minus_b_denom[6] = {64512, 1, 370944, 202752, 14336, 7084};
  int c_numer[6] = {0, 1, 3, 3, 1, 7};
  int c_denom[6] = {1, 5, 10, 5, 1, 8};
  int d_numer[6][5] = {{0, 0, 0, 0, 0},
                       {1, 0, 0, 0, 0},
                       {3, 9, 0, 0, 0},
                       {3, -9, 6, 0, 0},
                       {-11, 5, -70, 35, 0},
                       {1631, 175, 575, 44275, 253}};
  int d_denom[6][5] = {{1, 1, 1, 1, 1},
                       {5, 1, 1, 1, 1},
                       {40, 40, 1, 1, 1},
                       {10, 10, 5, 1, 1},
                       {54, 2, 27, 27, 1},
                       {55296, 512, 13824, 110592, 4096}};

  mpf_t a[6], a_minus_b[6], c[6], d[6][5];
  for (i = 0; i < 6; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);

    for (j = 0; j < 5; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  }

  // run the Runge-Kutta45 method
  rV = rk45_methods_mp_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 6; i++)
  {
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 5; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk45_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[6], double a_minus_b[6], double c[6], double d[6][5], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(52) - 3;
  double cond_num;
  comp_d newTime;
  point_d Y, dX[6];
  
  init_point_d(Y, 0);
  for (i = 0; i < 6; i++)
    init_point_d(dX[i], 0);
    
  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)
  
  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);
  
  for (k = 0; k < 6; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }
    
    if (k == 0) 
    {
      retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 5)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_d(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);  
          set_d(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_d(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rdouble_d(newTime, dT, c[k+1]);
        add_d(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
        vec_mulcomp_d(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    *sizeProportion = 0;
    increase_size_vec_d(Y, P->point->size);
    increase_size_vec_d(dZ, P->point->size);
    rows = Y->size = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    {
      mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rdouble_d(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 6; k++)
      {
        mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_d(&Y->coord[i], &Y->coord[i], newTime);

        mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
        add_d(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    *error = infNormVec_d(Y);

    // compute sizeProportion
    cond_num = pow(norm_sqr_d(dT), 2.5);
    if (*error == 0 || -log10(*error) >= digits || -log10(cond_num) >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional size
      *sizeProportion = 1;
    }
    else
    { // compute the proportional size
      *sizeProportion = pow(10, log10(*error) - log10(cond_num));
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk45 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk45 prediction.\n");
  }

  clear_point_d(Y);
  for (i = 0; i < 6; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rk45_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[6], mpf_t a_minus_b[6], mpf_t c[6], mpf_t d[6][5], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK45 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(T->Precision) - 5;
  double cond_num, logError;
  comp_mp newTime;
  point_mp Y, dX[6];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 6; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 6; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 5)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_mp(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_mp(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_mp(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rmpf_mp(newTime, dT, c[k+1]);
        add_mp(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
       vec_mulcomp_mp(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    increase_size_vec_mp(Y, dX[0]->size);
    increase_size_vec_mp(dZ, dX[0]->size);
    rows = Y->size = dZ->size = dX[0]->size;
    for (i = 0; i < rows; i++)
    {
      mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rmpf_mp(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 6; k++)
      {
        mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_mp(&Y->coord[i], &Y->coord[i], newTime);

        mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
        add_mp(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    infNormVec_mp2(error, Y);
    mpfr_log10(newTime->r, error, __gmp_default_rounding_mode);
    logError = mpf_get_d(newTime->r);

    // compute sizeProportion
    abs_mp(newTime, dT);
    mpf_pow_ui(newTime->i, newTime->r, 5);
    mpfr_log10(newTime->r, newTime->i, __gmp_default_rounding_mode);
    cond_num = mpf_get_d(newTime->r);
    if (mpfr_zero_p(error) || -logError >= digits || -cond_num >= 2 * digits) 
    { // numerics are too bad to even try to calculate the proportional error
      *sizeProportion = 1;
    }
    else
    {
      mpf_div(newTime->r, error, newTime->i);
      *sizeProportion = mpf_get_d(newTime->r);
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk45 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk45 prediction.\n");
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 6; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkdp56_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKDP56 method - Dormand/Prince method         *
\***************************************************************/
{   
  int rV;
  double a[8] = {821.0/10800, 0, 19683.0/71825, 175273.0/912600, 395.0/3672, 785.0/2704, 3.0/50, 0};
  double a_minus_b[8] = {13.0/2400, 0, -19683.0/618800, 2401.0/31200, -65.0/816, 15.0/416, 521.0/5600, -1.0/10};
  double c[8] = {0, 1.0/10, 2.0/9, 3.0/7, 3.0/5, 4.0/5, 1, 1};
  double d[8][7] = {{0, 0, 0, 0, 0, 0, 0},
                    {1.0/10, 0, 0, 0, 0, 0, 0},
                    {-2.0/81, 20.0/81, 0, 0, 0, 0, 0},
                    {615.0/1372, -270.0/343, 1053.0/1372, 0, 0, 0, 0},
                    {3243.0/5500, -54.0/55, 50949.0/71500, 4998.0/17875, 0, 0, 0},
                    {-26492.0/37125, 72.0/55, 2808.0/23375, -24206.0/37125, 338.0/459, 0, 0},
                    {5561.0/2376, -35.0/11, -24117.0/31603, 899983.0/200772, -5225.0/1836, 3925.0/4056, 0},
                    {465467.0/266112, -2945.0/1232, -5610201.0/14158144, 10513573.0/3212352, -424325.0/205632, 376225.0/454272, 0}};

  // run the Runge-Kutta56 method
  rV = rk56_methods_d_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkdp56_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKDP56 method - Dormand/Prince method         *
\***************************************************************/
{
  int i, j, rV;

  int a_numer[8] = {821, 0, 19683, 175273, 395, 785, 3, 0};
  int a_denom[8] = {10800, 1, 71825, 912600, 3672, 2704, 50, 1};
  int a_minus_b_numer[8] = {13, 0, -19683, 2401, -65, 15, 521, -1};
  int a_minus_b_denom[8] = {2400, 1, 618800, 31200, 816, 416, 5600, 10};
  int c_numer[8] = {0, 1, 2, 3, 3, 4, 1, 1};
  int c_denom[8] = {1, 10, 9, 7, 5, 5, 1, 1};
  int d_numer[8][7] = {{0, 0, 0, 0, 0, 0, 0},
                          {1, 0, 0, 0, 0, 0, 0},
                          {-2, 20, 0, 0, 0, 0, 0},
                          {615, -270, 1053, 0, 0, 0, 0},
                          {3243, -54, 50949, 4998, 0, 0, 0},
                          {-26492, 72, 2808, -24206, 338, 0, 0},
                          {5561, -35, -24117, 899983, -5225, 3925, 0},
                          {465467, -2945, -5610201, 10513573, -424325, 376225, 0}};
  int d_denom[8][7] = {{1, 1, 1, 1, 1, 1, 1},
                       {10, 1, 1, 1, 1, 1, 1},
                       {81, 81, 1, 1, 1, 1, 1},
                       {1372, 343, 1372, 1, 1, 1, 1},
                       {5500, 55, 71500, 17875, 1, 1, 1},
                       {37125, 55, 23375, 37125, 459, 1, 1},
                       {2376, 11, 31603, 200772, 1836, 4056, 1},
                       {266112, 1232, 14158144, 3212352, 205632, 454272, 1}};

  mpf_t a[8], a_minus_b[8], c[8], d[8][7];
  for (i = 0; i < 8; i++)
  {
    mpf_init(a[i]); mpf_set_si(a[i], a_numer[i]); mpf_div_ui(a[i], a[i], a_denom[i]);
    mpf_init(a_minus_b[i]); mpf_set_si(a_minus_b[i], a_minus_b_numer[i]); mpf_div_ui(a_minus_b[i], a_minus_b[i], a_minus_b_denom[i]);
    mpf_init(c[i]); mpf_set_si(c[i], c_numer[i]); mpf_div_ui(c[i], c[i], c_denom[i]);
  
    for (j = 0; j < 7; j++)
    {
      mpf_init(d[i][j]); mpf_set_si(d[i][j], d_numer[i][j]); mpf_div_ui(d[i][j], d[i][j], d_denom[i][j]);
    }
  } 
   
  // run the Runge-Kutta56 method
  rV = rk56_methods_mp_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);
    
  for (i = 0; i < 8; i++)
  { 
    mpf_clear(a[i]);
    mpf_clear(a_minus_b[i]);
    mpf_clear(c[i]);
    for (j = 0; j < 7; j++) 
      mpf_clear(d[i][j]);
  }

  return rV;
} 

int rk56_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[8], double a_minus_b[8], double c[8], double d[8][7], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK56 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(52) - 3;;
  double cond_num;
  comp_d newTime;
  point_d Y, dX[8];

  init_point_d(Y, 0);
  for (i = 0; i < 8; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 8; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 7)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_d(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_d(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_d(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rdouble_d(newTime, dT, c[k+1]);
        add_d(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
        vec_mulcomp_d(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    *sizeProportion = 0;
    increase_size_vec_d(Y, P->point->size);
    increase_size_vec_d(dZ, P->point->size);
    rows = Y->size = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    {
      mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rdouble_d(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 8; k++)
      {
        mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_d(&Y->coord[i], &Y->coord[i], newTime);

        mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
        add_d(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    *error = infNormVec_d(Y);

    // compute sizeProportion
    cond_num = pow(norm_sqr_d(dT), 3);
    if (*error == 0 || -log10(*error) >= digits || -log10(cond_num) >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional size
      *sizeProportion = 1;
    }
    else
    { // compute the proportional size
      *sizeProportion = pow(10, log10(*error) - log10(cond_num));
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk56 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk56 prediction.\n");
  }

  clear_point_d(Y);
  for (i = 0; i < 8; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rk56_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[8], mpf_t a_minus_b[8], mpf_t c[8], mpf_t d[8][7], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK56 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(T->Precision) - 5;
  double cond_num, logError;
  comp_mp newTime;
  point_mp Y, dX[8];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 8; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 8; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 7)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_mp(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_mp(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_mp(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rmpf_mp(newTime, dT, c[k+1]);
        add_mp(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
       vec_mulcomp_mp(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    increase_size_vec_mp(Y, dX[0]->size);
    increase_size_vec_mp(dZ, dX[0]->size);
    rows = Y->size = dZ->size = dX[0]->size;
    for (i = 0; i < rows; i++)
    {
      mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rmpf_mp(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 8; k++)
      {
        mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_mp(&Y->coord[i], &Y->coord[i], newTime);

        mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
        add_mp(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    infNormVec_mp2(error, Y);
    mpfr_log10(newTime->r, error, __gmp_default_rounding_mode);
    logError = mpf_get_d(newTime->r);

    // compute sizeProportion
    abs_mp(newTime, dT);
    mpf_pow_ui(newTime->i, newTime->r, 6);
    mpfr_log10(newTime->r, newTime->i, __gmp_default_rounding_mode);
    cond_num = mpf_get_d(newTime->r);
    if (mpfr_zero_p(error) || -logError >= digits || -cond_num >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional error
      *sizeProportion = 1;
    }
    else
    {
      mpf_div(newTime->r, error, newTime->i);
      *sizeProportion = mpf_get_d(newTime->r);
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk56 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk56 prediction.\n");
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 8; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int rkv67_d_amp(double *err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKV67 method - Verner method                  *
\***************************************************************/
{
  int rV;
  double a[10], a_minus_b[10], c[10], d[10][9];

a[0] = 0.0471556184862722217043176510883817567956932673866498263717863225066525727846555841973346401995827631062036755590558062534375047415234260388676624358999573840784918991178510655224244139232889913788333823544199823935451154623648793450503723002665841976598046959664007442181911565685230727013689593392982;
a[1] = 0;
a[2] = 0;
a[3] = 0.257505642984341518959643610103768758098639571258883586561584933570028952897166612637862059798925621609396907594245725006327476849765387793842037461290921370430264152893947313918692476051185173305706391851963836018457313459645894963869257643822167873501457534217624142946250030154297352910582419606143;
a[4] = 0.262166539774126204771386309576452771112914046736971319203584006643876812614319887681344771902494347386463073904093142446268038453678784950161124772440805335094258253338958758076230419920570973073193926469534220996564683552221244076754949769609944998538014278529773548771150093620971914747932175414550;
a[5] = 0.152160926567385574032313319916511753552259518225131674367096139045695165543077407261967511436494354991858803879708884258567824761923036313221569840599729116798542290424111444095404797846598555249340492159421091497761496624438637702252254645961603915777481188523759942071760307205888241753387037725744;
a[6] = 0.493996917003248424690717589322787684429585308571357948791678035642674335879997785662081156612226972186240516630330279108098440721497194764675468796945045160954090374742547119955506272386697362858820688622832619000271491070870242751247344166824766097061388163627211199461983342148211784138952128223367;
a[7] = -0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817;
a[8] = 0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153;
a[9] = 0;

a_minus_b[0] = 0.00254701187993104541699947511358977898137393503651067684217948288343659637043927033143148723131125158014121649938900878979858451840352977282793140381518240008963362033067932018366864275478963166454694523116781829305324081676250402401051739684741551799350892643994927331531661450056821725901262420206562;
a_minus_b[1] = 0;
a_minus_b[2] = 0;
a_minus_b[3] = -0.00965839487279574909126661599061503187520108884519528045821452710507663251123547574163069610539919389963989845238921904972854937606654792208299092082065749981679263365462373361333668366236519725184736326271819113885999268066494138424342590942206813226759750600367884910177105022385487684166900829606156;
a_minus_b[4] = 0.0420647097563969027734147319113774614805948824559621257010925807491537576600651305254295976225991591885211868879798925804616780192287010286942291836640369052613351689178787207331605464157207419903756767065167872280544436640338531219590767268952538968560279938316779313158803412904221702386669877100151;
a_minus_b[5] = -0.0666822437469301090659987634347776289054953577328478620986708910886522673765780883582819579617699092650805636589852395875931279767916652676612036293647748766874420794927278867725689176383726462247671413895793673537650626399046027049371138790153665954418043148360367277640961073763921508851946773219422;
a_minus_b[6] = 0.265009746462128136352900200346432447893361385631363781535494557795019991145447029415600259733188324877703429783214090406452549289091431065680284500120422909432522977485881928123107977267794033274265297233044517672944951473291445755463803660526233224001811988404937384560301642803203641418032914674657;
a_minus_b[7] = -0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817;
a_minus_b[8] = 0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153;
a_minus_b[9] = -0.0202951846633562822276705479381043035855420443667990862261517658249536055689205887319585505702055732014823934917746960666918489454776188166903872302377514709236100830696726470857731850092255075866785330602598147948274804639773599730786794693464008286038012269720794348562965112960546349366261206596328;

c[0] = 0;
c[1] = 0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[2] = 0.108888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889;
c[3] = 0.163333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
c[4] = 0.455500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[5] = 0.609509448997838131708700442148602494963796759286605673642939302386787496187674982310265299374033791406090583319455261759148494730544940521238957202640385592141084318199590688285015525265373754263537303527537491030484005373807797476632985866792541212386248529601147679387596870467556824874147722141014;
c[6] = 0.884000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[7] = 0.925000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[8] = 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
c[9] = 1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;

d[0][0] = 0;
d[0][1] = 0;
d[0][2] = 0;
d[0][3] = 0;
d[0][4] = 0;
d[0][5] = 0;
d[0][6] = 0;
d[0][7] = 0;
d[0][8] = 0;
d[1][0] = 0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
d[1][1] = 0;
d[1][2] = 0;
d[1][3] = 0;
d[1][4] = 0;
d[1][5] = 0;
d[1][6] = 0;
d[1][7] = 0;
d[1][8] = 0;
d[2][0] = -1.07679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012346;
d[2][1] = 1.18567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901235;
d[2][2] = 0;
d[2][3] = 0;
d[2][4] = 0;
d[2][5] = 0;
d[2][6] = 0;
d[2][7] = 0;
d[2][8] = 0;
d[3][0] = 0.0408333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
d[3][1] = 0;
d[3][2] = 0.122500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000;
d[3][3] = 0;
d[3][4] = 0;
d[3][5] = 0;
d[3][6] = 0;
d[3][7] = 0;
d[3][8] = 0;
d[4][0] = 0.638913923625572678050812161599333610995418575593502707205331112036651395251978342357351103706788837984173261141191170345689296126613910870470637234485630987088713036234902124114952103290295710120783007080383173677634319033735943356934610578925447730112453144523115368596418159100374843815077051228655;
d[4][1] = 0;
d[4][2] = -2.45567263822365680966264056643065389421074552269887546855476884631403581840899625156184922948771345272802998750520616409829237817576009995835068721366097459391920033319450229071220324864639733444398167430237401082882132444814660558100791336942940441482715535193669304456476468138275718450645564348188;
d[4][3] = 2.27225871459808413161182840483132028321532694710537276134943773427738442315701790920449812578092461474385672636401499375260308204914618908788004997917534360683048729695960016659725114535610162432319866722199083715118700541441066222407330279050395668471470220741357767596834652228238234069137859225323;
d[4][4] = 0;
d[4][5] = 0;
d[4][6] = 0;
d[4][7] = 0;
d[4][8] = 0;
d[5][0] = -2.66157737501875713111925929786181811927886110621279491101469713930509612941385509594430505557023368554958603063493231787317799674230845336779694159618434634321299954680555583166102818777787674050555772208977921672782109407296127769462859055244305774353176625390759049478562797393162386458986585502767;
d[5][1] = 0;
d[5][2] = 10.8045138864561376956539665536553283848216695313861604686310006136290059392424850728726033575720412843056269529183992807963144686325598710039491665420757292749882852097214396218692807538243757762105076420989844983852423950323217516089733038440975804125387258140136771071745292560309970797058687391147;
d[5][3] = -8.35391465739619941196804854781929169154088257610338993266161781790756725509495082412344107605516870213559043057090330898930263408054993414287063406173197791833013415175411049756858203990308618552762172172286334472369559218815026430538538071155831209364354686322679692202799997090849267285704422656095;
d[5][4] = 0.820487594956656979142041734174383920961870910216630048688253645970444941453995829505408073427394894785640091606891607825314656920843457027957366318480980578695932807037817395645344999121960904086209105241195554096758296602597587867673653286696330637022835832721857989026695559276676282615189064614936;
d[5][5] = 0;
d[5][6] = 0;
d[5][7] = 0;
d[5][8] = 0;
d[6][0] = 6.06774143469677099271836018387727671467926640176887306359909941550772632909985821015131404165728812111281168328257663180016064215895470690198142759429883346899999316325970007167545285470047913584083460611017718063130166260596491115835032527423915877478243890066807867827599446375139672648735819080763;
d[6][1] = 0;
d[6][2] = -24.7112736359110857973420348529074600180341771629677540813665554596587326431239932816061210269509909363811041380357177346279533645394890272496798875375316795980503341386187230715331278943723231323700151408836595698414095201199019330839199397683748689498911238776565489210858226525290942123615012204069;
d[6][3] = 20.4275179307888939404577311174834661269660931649405558961215868278847940059501212187832475382874361032894160252654902848069119955759637823202167209261918250178828181639508851014288273511443537985849961840706061832183847204574232612353992610081640468028660624990850592258752345152578286375382830911574;
d[6][4] = -1.90615797881664715062409678435275701087897468130154449732332455905089888161653727987321175840553252820469326122182450914811549396345902609080269405892288783291419739936108774599029695502963878984287699984458255527075364785177970366374833241148008090451090875831710847373226588365269422264785656055595;
d[6][5] = 1.00617224924206801479004033589947418726779227755986961896919377531711118969055113254477120541179924018356969070947532716899622076802956411828443307596390894408172021076922564441914464355712898778706135054745876126247678490829346435391868589745174427675353123622051949066685955717256307098371649899780;
d[6][6] = 0;
d[6][7] = 0;
d[6][8] = 0;
d[7][0] = 12.0546700762532029950910945289277831164804252410675022456504679687895689139960146039443288283750336379411096035825954447760936334776440213949443362851020478648773614432833479739455614899658192055871179070263606189309615200393669693029537872090663785823585220176083264357375138384872002392464453156991;
d[7][1] = 0;
d[7][2] = -49.7547849504689893280725761533144475832181722583426930533724074758711662694085370746679713708359769175525414749792979166191877040702950632035050470020509416758662857711589354027529455771937047867243267859099031357708288506308006555337876940221477891866864320545092351815961637286145785810654620897866;
d[7][3] = 41.1428886386046766325969841671015735420861993168074196817767197622564349824524439840170095382691179158858647537995469065281803314037669775100891174117560375843405111217539168488777686205987779081576532479718384543237590330774192826000651008343724280693109569745879303050600351583118738080778029866863;
d[7][4] = -4.46176014997400418564191160348481537505077261504527943671205073698503173782183498817306004035577444077774338030344144153116606679374731727848251310165306184178863494208558304215340239111889056504039756992455997594732226828711730693128980332843828825010470766013799305694315017127922117864625331835752;
d[7][5] = 2.04233482223917495982171707770860854373769387362950704675226264453594041136972385119392178374708596422006064736887487517564055580246077638154888337104289952170304921640066775478538257502794816329229535477961389665728962587873947154164966547885933962957589901056567981878726157221164947050255341552053;
d[7][6] = -0.0983484366540610737953080169387022440353735581164564840949921627257463005878103763142287391994861597167501494682778683295607498198293948045947769641969814532660010681934141327023647172799499252723421539433498581938590600776077609795910561717120688444542382881147083210454966691169237581150863097617698;
d[7][7] = 0;
d[7][8] = 0;
d[8][0] = 10.1381465228818078764184514198168903076855566197446625878810302234504456896166594014779323382898853196615353532636712517123358396606673215965993494780266380948135104983713885563789845413826484979295440534003324926587827813571056844823251010555486301677111922636549074237942689192388927633686568909689;
d[8][1] = 0;
d[8][2] = -42.6411360317175021462284600673663573062458353006064366534390081892169164380334190644420973788613121461309320350305871333287662229028588588586066541960004850314875558898834391865978268773660940651739460069750441404794804377117425966233263555407232236441737361866796074494069051528927970522898826975093;
d[8][3] = 35.7638400399225700713502117802316005403429826863865820113565153536300404598612532898965354818454694648874725072682749267373157568233054310410187975820486383990095600494571284321105493325839468611507648538420067322279253828587432646374273647616387422018633170284509612151810792245462886190112774102403;
d[8][4] = -4.34802284039290765334037029690824594371028341639116518865331287370628888237874447914620105044959090476557742632948537833967528093187874711244268408364492971056602670317493656831158047324970779777483136577206179452089410068354309842483825320643538174054369045758176371157191719967939556428314730251243;
d[8][5] = 2.00986226837703589544194359301182755477072972921068058748766475060301022678114409555285903080395717143894554579171389259003194150520251106072972691640555512448806839231389179903663905816067082512966845731905978567106788706342683552320444934837670908237905142357004848267401763746008483191555662392919;
d[8][6] = 0.348749046033827240595382285305314587914041325248182841863580609691330339196441065124674708479092059422639244419745826821459351853291206642381263602952792439394331494334732664490055557979195367934536312102247333268705820911625949222101271883714392689492479451747299183798095407000894534464887545496646;
d[8][7] = -0.271439005104831284237158714091029740757191643592506186496469874451621395043334308463703130107500964514083189383333386192701386007728864369679799299788209315651887841418765697106821139490659689195736303916540408826107333795616038816893578302119868756728613523161845144468638835673968132187348470613286;
d[8][8] = 0;
d[9][0] = -45.0300720342986771243532240507376963515145444719937644824079824393198401727221975824532835519581554198813324085760134180570787821182548227893314627083718036426491592921518503051297688605921310892888780618959209224202919510149533695706542210627124350751677569019920750022536127060269294301531394376504;
d[9][1] = 0;
d[9][2] = 187.327243765458884075241820615420199738351561133306193730472042811504816158488076206432002582649614412920275002256672174701510739653751114685455028887061290838474955260652565481309225324569424103179332013504226474976360325989130909739828895637852868788385760452234355087997545517424616009154634705758;
d[9][3] = -154.028823693501869059672862103451040258175735549044835668465248022257814538167848348770086411169018627838612812103724152543810567982425588659759737860803806260371145723240966294468826535262918904290878915967941233375788211216428688383225716305408979344377644058147926817639142528931780694303713181806;
d[9][4] = 18.5646530634753623385949233295843913676513806454050558851549018133209261097165524415954345328693490705442847387683039834385800075030766667530814210089033384064285859281150173647602629309861028182523307187208812695652827444150045993819313815469829521981879407525662683342644692190497470756749789063511;
d[9][5] = -7.14180967929507885492542049682355119282094179375667055538129282740352869535208365846699842107061792731728890413284789902812520734947617347293091211988712457899696810183968746594422973631393050228738091580973160939356720771138507764901376187221457776450991563372884707638237400874442565117025303276579;
d[9][6] = 1.30880857816137862511476270600769669650828003608402109062757866415544113803750094166293126867882849157267438378760931148892381029332880348348566279309810523711373192846492121947333687661345357443547516144848602064800429953863162648113342205550017119748161538906822547401311450722877269079749204011337;
d[9][7] = 0;
d[9][8] = 0;

  // run the Runge-Kutta56 method
  rV = rk67_methods_d_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  return rV;
}

int rkv67_mp_amp(mpf_t err, double *norm_J, double *norm_J_inv, double *sizeProportion, vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the RKV67 method - Verner method                  *
\***************************************************************/
{
  int i, j, rV;

  mpf_t a[10], a_minus_b[10], c[10], d[10][9];
  for (i = 0; i < 10; i++)
  {
    mpf_init(a[i]); mpf_init(a_minus_b[i]); mpf_init(c[i]);
    for (j = 0; j < 9; j++)
      mpf_init(d[i][j]);
  }

mpf_set_str(a[0],"0.0471556184862722217043176510883817567956932673866498263717863225066525727846555841973346401995827631062036755590558062534375047415234260388676624358999573840784918991178510655224244139232889913788333823544199823935451154623648793450503723002665841976598046959664007442181911565685230727013689593392982",10);
mpf_set_str(a[1],"0",10);
mpf_set_str(a[2],"0",10);
mpf_set_str(a[3],"0.257505642984341518959643610103768758098639571258883586561584933570028952897166612637862059798925621609396907594245725006327476849765387793842037461290921370430264152893947313918692476051185173305706391851963836018457313459645894963869257643822167873501457534217624142946250030154297352910582419606143",10);
mpf_set_str(a[4],"0.262166539774126204771386309576452771112914046736971319203584006643876812614319887681344771902494347386463073904093142446268038453678784950161124772440805335094258253338958758076230419920570973073193926469534220996564683552221244076754949769609944998538014278529773548771150093620971914747932175414550",10);
mpf_set_str(a[5],"0.152160926567385574032313319916511753552259518225131674367096139045695165543077407261967511436494354991858803879708884258567824761923036313221569840599729116798542290424111444095404797846598555249340492159421091497761496624438637702252254645961603915777481188523759942071760307205888241753387037725744",10);
mpf_set_str(a[6],"0.493996917003248424690717589322787684429585308571357948791678035642674335879997785662081156612226972186240516630330279108098440721497194764675468796945045160954090374742547119955506272386697362858820688622832619000271491070870242751247344166824766097061388163627211199461983342148211784138952128223367",10);
mpf_set_str(a[7],"-0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817",10);
mpf_set_str(a[8],"0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153",10);
mpf_set_str(a[9],"0",10);

mpf_set_str(a_minus_b[0],"0.00254701187993104541699947511358977898137393503651067684217948288343659637043927033143148723131125158014121649938900878979858451840352977282793140381518240008963362033067932018366864275478963166454694523116781829305324081676250402401051739684741551799350892643994927331531661450056821725901262420206562",10);
mpf_set_str(a_minus_b[1],"0",10);
mpf_set_str(a_minus_b[2],"0",10);
mpf_set_str(a_minus_b[3],"-0.00965839487279574909126661599061503187520108884519528045821452710507663251123547574163069610539919389963989845238921904972854937606654792208299092082065749981679263365462373361333668366236519725184736326271819113885999268066494138424342590942206813226759750600367884910177105022385487684166900829606156",10);
mpf_set_str(a_minus_b[4],"0.0420647097563969027734147319113774614805948824559621257010925807491537576600651305254295976225991591885211868879798925804616780192287010286942291836640369052613351689178787207331605464157207419903756767065167872280544436640338531219590767268952538968560279938316779313158803412904221702386669877100151",10);
mpf_set_str(a_minus_b[5],"-0.0666822437469301090659987634347776289054953577328478620986708910886522673765780883582819579617699092650805636589852395875931279767916652676612036293647748766874420794927278867725689176383726462247671413895793673537650626399046027049371138790153665954418043148360367277640961073763921508851946773219422",10);
mpf_set_str(a_minus_b[6],"0.265009746462128136352900200346432447893361385631363781535494557795019991145447029415600259733188324877703429783214090406452549289091431065680284500120422909432522977485881928123107977267794033274265297233044517672944951473291445755463803660526233224001811988404937384560301642803203641418032914674657",10);
mpf_set_str(a_minus_b[7],"-0.294303117140325044155724474409270342913873557068949824155892157793208589321399942529392350490770742507468780043904195897239587492035877794842026930898306757155719243786459153030937692324381551585672032706362933754854709864697658247349641770296612440560976646488147527504393807577038069940675914010817",10);
mpf_set_str(a_minus_b[8],"0.0813174723249510999973459944013676189247818448899554688601627203842807496021826650888022105410466832273058024764703588245403019636480479340741636237218483898000722732690434514626793121960404957197771512481911838482546096951567594081754632438115453580228307856233779500350588778791457036884531937017153",10);
mpf_set_str(a_minus_b[9],"-0.0202951846633562822276705479381043035855420443667990862261517658249536055689205887319585505702055732014823934917746960666918489454776188166903872302377514709236100830696726470857731850092255075866785330602598147948274804639773599730786794693464008286038012269720794348562965112960546349366261206596328",10);

mpf_set_str(c[0],"0",10);
mpf_set_str(c[1],"0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[2],"0.108888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889",10);
mpf_set_str(c[3],"0.163333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333",10);
mpf_set_str(c[4],"0.455500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[5],"0.609509448997838131708700442148602494963796759286605673642939302386787496187674982310265299374033791406090583319455261759148494730544940521238957202640385592141084318199590688285015525265373754263537303527537491030484005373807797476632985866792541212386248529601147679387596870467556824874147722141014",10);
mpf_set_str(c[6],"0.884000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[7],"0.925000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[8],"1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(c[9],"1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);

mpf_set_str(d[0][0],"0",10);
mpf_set_str(d[0][1],"0",10);
mpf_set_str(d[0][2],"0",10);
mpf_set_str(d[0][3],"0",10);
mpf_set_str(d[0][4],"0",10);
mpf_set_str(d[0][5],"0",10);
mpf_set_str(d[0][6],"0",10);
mpf_set_str(d[0][7],"0",10);
mpf_set_str(d[0][8],"0",10);
mpf_set_str(d[1][0],"0.00500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(d[1][1],"0",10);
mpf_set_str(d[1][2],"0",10);
mpf_set_str(d[1][3],"0",10);
mpf_set_str(d[1][4],"0",10);
mpf_set_str(d[1][5],"0",10);
mpf_set_str(d[1][6],"0",10);
mpf_set_str(d[1][7],"0",10);
mpf_set_str(d[1][8],"0",10);
mpf_set_str(d[2][0],"-1.07679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012345679012346",10);
mpf_set_str(d[2][1],"1.18567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901234567901235",10);
mpf_set_str(d[2][2],"0",10);
mpf_set_str(d[2][3],"0",10);
mpf_set_str(d[2][4],"0",10);
mpf_set_str(d[2][5],"0",10);
mpf_set_str(d[2][6],"0",10);
mpf_set_str(d[2][7],"0",10);
mpf_set_str(d[2][8],"0",10);
mpf_set_str(d[3][0],"0.0408333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333",10);
mpf_set_str(d[3][1],"0",10);
mpf_set_str(d[3][2],"0.122500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
mpf_set_str(d[3][3],"0",10);
mpf_set_str(d[3][4],"0",10);
mpf_set_str(d[3][5],"0",10);
mpf_set_str(d[3][6],"0",10);
mpf_set_str(d[3][7],"0",10);
mpf_set_str(d[3][8],"0",10);
mpf_set_str(d[4][0],"0.638913923625572678050812161599333610995418575593502707205331112036651395251978342357351103706788837984173261141191170345689296126613910870470637234485630987088713036234902124114952103290295710120783007080383173677634319033735943356934610578925447730112453144523115368596418159100374843815077051228655",10);
mpf_set_str(d[4][1],"0",10);
mpf_set_str(d[4][2],"-2.45567263822365680966264056643065389421074552269887546855476884631403581840899625156184922948771345272802998750520616409829237817576009995835068721366097459391920033319450229071220324864639733444398167430237401082882132444814660558100791336942940441482715535193669304456476468138275718450645564348188",10);
mpf_set_str(d[4][3],"2.27225871459808413161182840483132028321532694710537276134943773427738442315701790920449812578092461474385672636401499375260308204914618908788004997917534360683048729695960016659725114535610162432319866722199083715118700541441066222407330279050395668471470220741357767596834652228238234069137859225323",10);
mpf_set_str(d[4][4],"0",10);
mpf_set_str(d[4][5],"0",10);
mpf_set_str(d[4][6],"0",10);
mpf_set_str(d[4][7],"0",10);
mpf_set_str(d[4][8],"0",10);
mpf_set_str(d[5][0],"-2.66157737501875713111925929786181811927886110621279491101469713930509612941385509594430505557023368554958603063493231787317799674230845336779694159618434634321299954680555583166102818777787674050555772208977921672782109407296127769462859055244305774353176625390759049478562797393162386458986585502767",10);
mpf_set_str(d[5][1],"0",10);
mpf_set_str(d[5][2],"10.8045138864561376956539665536553283848216695313861604686310006136290059392424850728726033575720412843056269529183992807963144686325598710039491665420757292749882852097214396218692807538243757762105076420989844983852423950323217516089733038440975804125387258140136771071745292560309970797058687391147",10);
mpf_set_str(d[5][3],"-8.35391465739619941196804854781929169154088257610338993266161781790756725509495082412344107605516870213559043057090330898930263408054993414287063406173197791833013415175411049756858203990308618552762172172286334472369559218815026430538538071155831209364354686322679692202799997090849267285704422656095",10);
mpf_set_str(d[5][4],"0.820487594956656979142041734174383920961870910216630048688253645970444941453995829505408073427394894785640091606891607825314656920843457027957366318480980578695932807037817395645344999121960904086209105241195554096758296602597587867673653286696330637022835832721857989026695559276676282615189064614936",10);
mpf_set_str(d[5][5],"0",10);
mpf_set_str(d[5][6],"0",10);
mpf_set_str(d[5][7],"0",10);
mpf_set_str(d[5][8],"0",10);
mpf_set_str(d[6][0],"6.06774143469677099271836018387727671467926640176887306359909941550772632909985821015131404165728812111281168328257663180016064215895470690198142759429883346899999316325970007167545285470047913584083460611017718063130166260596491115835032527423915877478243890066807867827599446375139672648735819080763",10);
mpf_set_str(d[6][1],"0",10);
mpf_set_str(d[6][2],"-24.7112736359110857973420348529074600180341771629677540813665554596587326431239932816061210269509909363811041380357177346279533645394890272496798875375316795980503341386187230715331278943723231323700151408836595698414095201199019330839199397683748689498911238776565489210858226525290942123615012204069",10);
mpf_set_str(d[6][3],"20.4275179307888939404577311174834661269660931649405558961215868278847940059501212187832475382874361032894160252654902848069119955759637823202167209261918250178828181639508851014288273511443537985849961840706061832183847204574232612353992610081640468028660624990850592258752345152578286375382830911574",10);
mpf_set_str(d[6][4],"-1.90615797881664715062409678435275701087897468130154449732332455905089888161653727987321175840553252820469326122182450914811549396345902609080269405892288783291419739936108774599029695502963878984287699984458255527075364785177970366374833241148008090451090875831710847373226588365269422264785656055595",10);
mpf_set_str(d[6][5],"1.00617224924206801479004033589947418726779227755986961896919377531711118969055113254477120541179924018356969070947532716899622076802956411828443307596390894408172021076922564441914464355712898778706135054745876126247678490829346435391868589745174427675353123622051949066685955717256307098371649899780",10);
mpf_set_str(d[6][6],"0",10);
mpf_set_str(d[6][7],"0",10);
mpf_set_str(d[6][8],"0",10);
mpf_set_str(d[7][0],"12.0546700762532029950910945289277831164804252410675022456504679687895689139960146039443288283750336379411096035825954447760936334776440213949443362851020478648773614432833479739455614899658192055871179070263606189309615200393669693029537872090663785823585220176083264357375138384872002392464453156991",10);
mpf_set_str(d[7][1],"0",10);
mpf_set_str(d[7][2],"-49.7547849504689893280725761533144475832181722583426930533724074758711662694085370746679713708359769175525414749792979166191877040702950632035050470020509416758662857711589354027529455771937047867243267859099031357708288506308006555337876940221477891866864320545092351815961637286145785810654620897866",10);
mpf_set_str(d[7][3],"41.1428886386046766325969841671015735420861993168074196817767197622564349824524439840170095382691179158858647537995469065281803314037669775100891174117560375843405111217539168488777686205987779081576532479718384543237590330774192826000651008343724280693109569745879303050600351583118738080778029866863",10);
mpf_set_str(d[7][4],"-4.46176014997400418564191160348481537505077261504527943671205073698503173782183498817306004035577444077774338030344144153116606679374731727848251310165306184178863494208558304215340239111889056504039756992455997594732226828711730693128980332843828825010470766013799305694315017127922117864625331835752",10);
mpf_set_str(d[7][5],"2.04233482223917495982171707770860854373769387362950704675226264453594041136972385119392178374708596422006064736887487517564055580246077638154888337104289952170304921640066775478538257502794816329229535477961389665728962587873947154164966547885933962957589901056567981878726157221164947050255341552053",10);
mpf_set_str(d[7][6],"-0.0983484366540610737953080169387022440353735581164564840949921627257463005878103763142287391994861597167501494682778683295607498198293948045947769641969814532660010681934141327023647172799499252723421539433498581938590600776077609795910561717120688444542382881147083210454966691169237581150863097617698",10);
mpf_set_str(d[7][7],"0",10);
mpf_set_str(d[7][8],"0",10);
mpf_set_str(d[8][0],"10.1381465228818078764184514198168903076855566197446625878810302234504456896166594014779323382898853196615353532636712517123358396606673215965993494780266380948135104983713885563789845413826484979295440534003324926587827813571056844823251010555486301677111922636549074237942689192388927633686568909689",10);
mpf_set_str(d[8][1],"0",10);
mpf_set_str(d[8][2],"-42.6411360317175021462284600673663573062458353006064366534390081892169164380334190644420973788613121461309320350305871333287662229028588588586066541960004850314875558898834391865978268773660940651739460069750441404794804377117425966233263555407232236441737361866796074494069051528927970522898826975093",10);
mpf_set_str(d[8][3],"35.7638400399225700713502117802316005403429826863865820113565153536300404598612532898965354818454694648874725072682749267373157568233054310410187975820486383990095600494571284321105493325839468611507648538420067322279253828587432646374273647616387422018633170284509612151810792245462886190112774102403",10);
mpf_set_str(d[8][4],"-4.34802284039290765334037029690824594371028341639116518865331287370628888237874447914620105044959090476557742632948537833967528093187874711244268408364492971056602670317493656831158047324970779777483136577206179452089410068354309842483825320643538174054369045758176371157191719967939556428314730251243",10);
mpf_set_str(d[8][5],"2.00986226837703589544194359301182755477072972921068058748766475060301022678114409555285903080395717143894554579171389259003194150520251106072972691640555512448806839231389179903663905816067082512966845731905978567106788706342683552320444934837670908237905142357004848267401763746008483191555662392919",10);
mpf_set_str(d[8][6],"0.348749046033827240595382285305314587914041325248182841863580609691330339196441065124674708479092059422639244419745826821459351853291206642381263602952792439394331494334732664490055557979195367934536312102247333268705820911625949222101271883714392689492479451747299183798095407000894534464887545496646",10);
mpf_set_str(d[8][7],"-0.271439005104831284237158714091029740757191643592506186496469874451621395043334308463703130107500964514083189383333386192701386007728864369679799299788209315651887841418765697106821139490659689195736303916540408826107333795616038816893578302119868756728613523161845144468638835673968132187348470613286",10);
mpf_set_str(d[8][8],"0",10);
mpf_set_str(d[9][0],"-45.0300720342986771243532240507376963515145444719937644824079824393198401727221975824532835519581554198813324085760134180570787821182548227893314627083718036426491592921518503051297688605921310892888780618959209224202919510149533695706542210627124350751677569019920750022536127060269294301531394376504",10);
mpf_set_str(d[9][1],"0",10);
mpf_set_str(d[9][2],"187.327243765458884075241820615420199738351561133306193730472042811504816158488076206432002582649614412920275002256672174701510739653751114685455028887061290838474955260652565481309225324569424103179332013504226474976360325989130909739828895637852868788385760452234355087997545517424616009154634705758",10);
mpf_set_str(d[9][3],"-154.028823693501869059672862103451040258175735549044835668465248022257814538167848348770086411169018627838612812103724152543810567982425588659759737860803806260371145723240966294468826535262918904290878915967941233375788211216428688383225716305408979344377644058147926817639142528931780694303713181806",10);
mpf_set_str(d[9][4],"18.5646530634753623385949233295843913676513806454050558851549018133209261097165524415954345328693490705442847387683039834385800075030766667530814210089033384064285859281150173647602629309861028182523307187208812695652827444150045993819313815469829521981879407525662683342644692190497470756749789063511",10);
mpf_set_str(d[9][5],"-7.14180967929507885492542049682355119282094179375667055538129282740352869535208365846699842107061792731728890413284789902812520734947617347293091211988712457899696810183968746594422973631393050228738091580973160939356720771138507764901376187221457776450991563372884707638237400874442565117025303276579",10);
mpf_set_str(d[9][6],"1.30880857816137862511476270600769669650828003608402109062757866415544113803750094166293126867882849157267438378760931148892381029332880348348566279309810523711373192846492121947333687661345357443547516144848602064800429953863162648113342205550017119748161538906822547401311450722877269079749204011337",10);
mpf_set_str(d[9][7],"0",10);
mpf_set_str(d[9][8],"0",10);

  // run the Runge-Kutta67 method
  rV = rk67_methods_mp_amp(err, norm_J, norm_J_inv, sizeProportion, a, a_minus_b, c, d, dZ, P, T, OUT, dT, e, ED, eval_func);

  for (i = 0; i < 10; i++)
  {
    mpf_clear(a[i]); mpf_clear(a_minus_b[i]); mpf_clear(c[i]);
    for (j = 0; j < 9; j++)
      mpf_clear(d[i][j]);
  }

  return rV;
}

int rk67_methods_d_amp(double *error, double *norm_J, double *norm_J_inv, double *sizeProportion, double a[10], double a_minus_b[10], double c[10], double d[10][9], vec_d dZ, point_data_d *P, tracker_config_t *T, FILE *OUT, comp_d dT, eval_struct_d *e, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK67 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(52) - 3;;
  double cond_num;
  comp_d newTime;
  point_d Y, dX[10];

  init_point_d(Y, 0);
  for (i = 0; i < 10; i++)
    init_point_d(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 10; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_d(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_d(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_d(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_d(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_d(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
    else
    { // just do matrixSolve
      retVal = matrixSolve_d(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 9)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_d(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_d(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_d(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rdouble_d(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_d(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rdouble_d(newTime, dT, c[k+1]);
        add_d(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_d(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
        vec_mulcomp_d(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    *sizeProportion = 0;
    increase_size_vec_d(Y, P->point->size);
    increase_size_vec_d(dZ, P->point->size);
    rows = Y->size = dZ->size = P->point->size;
    for (i = 0; i < rows; i++)
    {
      mul_rdouble_d(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rdouble_d(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 10; k++)
      {
        mul_rdouble_d(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_d(&Y->coord[i], &Y->coord[i], newTime);

        mul_rdouble_d(newTime, &dX[k]->coord[i], a[k]);
        add_d(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    *error = infNormVec_d(Y);

    // compute sizeProportion
    cond_num = pow(norm_sqr_d(dT), 3.5);
    if (*error == 0 || -log10(*error) >= digits || -log10(cond_num) >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional size
      *sizeProportion = 1;
    }
    else
    { // compute the proportional size
      *sizeProportion = pow(10, log10(*error) - log10(cond_num));
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk67 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk67 prediction.\n");
  }

  clear_point_d(Y);
  for (i = 0; i < 10; i++)
    clear_point_d(dX[i]);

  return retVal;
}

int rk67_methods_mp_amp(mpf_t error, double *norm_J, double *norm_J_inv, double *sizeProportion, mpf_t a[10], mpf_t a_minus_b[10], mpf_t c[10], mpf_t d[10][9], vec_mp dZ, point_data_mp *P, tracker_config_t *T, FILE *OUT, comp_mp dT, eval_struct_mp *e, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: runs the RK67 method using a, a-b, c, & d              *
\***************************************************************/
{
  int i, j, k, rows, cols, retVal = 0, digits = prec_to_digits(T->Precision) - 5;
  double cond_num, logError;
  comp_mp newTime;
  point_mp Y, dX[10];

  // initialize MP
  init_mp(newTime);
  init_point_mp(Y, 0);
  for (i = 0; i < 10; i++)
    init_point_mp(dX[i], 0);

  // find dX[k] = dT*f(x + sum(j<k, d[k][j] * dX[j]),t + c[k]*dT) = -(Jv^-1)*(Jp*dp/dt)

  // do the first function evaluation
  eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, P->point, P->time, ED, eval_func);

  for (k = 0; k < 10; k++)
  { // find Y = -dH/dt = -dH/dp * dp/dt
    increase_size_vec_mp(Y, e->Jp->rows);
    rows = Y->size = e->Jp->rows;
    cols = e->Jp->cols;
    for (i = 0; i < rows; i++)
    {
      set_zero_mp(&Y->coord[i]);
      for (j = 0; j < cols; j++)
      {
        sum_mul_mp(&Y->coord[i], &e->Jp->entry[i][j], &e->parDer->coord[j]);
      }
      neg_mp(&Y->coord[i], &Y->coord[i]);
    }

    if (k == 0)
    {
      retVal = matrixSolve_cond_num_norms_mp(dX[k], e->Jv, Y, &cond_num, norm_J, norm_J_inv);
      T->steps_since_last_CN = 0;
      T->latest_cond_num_exp = ceil(log10(cond_num));
    }
     else
    { // just do matrixSolve
      retVal = matrixSolve_mp(dX[k], e->Jv, Y);
    }
    // matrixSolve calculated dX[k] = (dH/dX)^(-1) * Y

    // check for matrixSolve failures
    if (retVal)
    {
      break;
    }
    else
    { // determine if more function evaluations are needed
      if (k < 9)
      { // compute dX[k] = dX[k] * dT
        // find Y = x + sum(j <= k, d[k+1][j] * dX[j])
        increase_size_vec_mp(Y, dX[k]->size);
        rows = Y->size = dX[k]->size;
        for (i = 0; i < rows; i++)
        {
          mul_mp(&dX[k]->coord[i], &dX[k]->coord[i], dT);
          set_mp(&Y->coord[i], &P->point->coord[i]);
          for (j = 0; j <= k; j++)
          {
            mul_rmpf_mp(newTime, &dX[j]->coord[i], d[k+1][j]);
            add_mp(&Y->coord[i], &Y->coord[i], newTime);
          }
        }

        // find newTime = t + c[k+1]*dT
        mul_rmpf_mp(newTime, dT, c[k+1]);
        add_mp(newTime, newTime, P->time);

        // do the next - for k+1 - function evaluation
        eval_mp(e->funcVals, e->parVals, e->parDer, e->Jv, e->Jp, Y, newTime, ED, eval_func);
      }
      else
      { // calculate dX[k] = dT * dX[k]
       vec_mulcomp_mp(dX[k], dX[k], dT);
      }
    }
  }

  if (!retVal)
  { // find error estimate: sum(k = 0..5, (a_k - b_k) * dX[k]) sizeProportion * || dT ||^5
    // find dZ = sum(k = 0..5, a_k * dX[k])
    increase_size_vec_mp(Y, dX[0]->size);
    increase_size_vec_mp(dZ, dX[0]->size);
    rows = Y->size = dZ->size = dX[0]->size;
    for (i = 0; i < rows; i++)
    {
      mul_rmpf_mp(&Y->coord[i], &dX[0]->coord[i], a_minus_b[0]);
      mul_rmpf_mp(&dZ->coord[i], &dX[0]->coord[i], a[0]);
      for (k = 1; k < 10; k++)
      {
        mul_rmpf_mp(newTime, &dX[k]->coord[i], a_minus_b[k]);
        add_mp(&Y->coord[i], &Y->coord[i], newTime);

        mul_rmpf_mp(newTime, &dX[k]->coord[i], a[k]);
        add_mp(&dZ->coord[i], &dZ->coord[i], newTime);
      }
    }
    // compute the norm of Y
    infNormVec_mp2(error, Y);
    mpfr_log10(newTime->r, error, __gmp_default_rounding_mode);
    logError = mpf_get_d(newTime->r);

    // compute sizeProportion
    abs_mp(newTime, dT);
    mpf_pow_ui(newTime->i, newTime->r, 7);
    mpfr_log10(newTime->r, newTime->i, __gmp_default_rounding_mode);
    cond_num = mpf_get_d(newTime->r);
    if (mpfr_zero_p(error) || -logError >= digits || -cond_num >= 2 * digits)
    { // numerics are too bad to even try to calculate the proportional error
      *sizeProportion = 1;
    }
    else
    {
      mpf_div(newTime->r, error, newTime->i);
      *sizeProportion = mpf_get_d(newTime->r);
    }

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
    fprintf(OUT, "NOTE: matrixSolve has failed when doing rk67 prediction.\n");
    if (T->screenOut)
      fprintf(stderr, "NOTE: matrixSolve has failed when doing rk67 prediction.\n");
  }

  // clear MP
  clear_mp(newTime);
  clear_point_mp(Y);
  for (i = 0; i < 10; i++)
    clear_point_mp(dX[i]);

  return retVal;
}

int AMP3_success_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, int order, double *currStepSize, comp_d currTime, tracker_config_t *T)
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

  // run AMP3 to see what we should do
  AMP3_update(&bestDigits, &bestPrec, currStepSize, stepSize_mp, currPrec, P0, digits_C, order, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, currTime, NULL);

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

int AMP3_success_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, int order, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T, int *prec_decreases, int max_prec_decreases)
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

  // run AMP3 to see what we should do
  AMP3_update(&bestDigits, &bestPrec, &stepSize_d, currStepSize, curr_prec, P0, digits_C, order, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, NULL, currTime);

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

void AMP3_minimize_cost(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, double eta_max, double P0, int order, int maxNewtonIts, int minPrecNeeded, int maxPrecNeeded)
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
  eta[0] = (P0 - prec_to_digits(prec[0])) * maxNewtonIts / (order + 1.0);
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
    eta[i] = (P0 - prec_to_digits(prec[i])) * maxNewtonIts / (order + 1.0);
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

void AMP3_update(int *newDigits, int *newPrec, double *currStepSize_d, mpf_t currStepSize_mp, int currPrec, double P0, int digits_C, int order, double newtonTol, int maxNewtonIts, double finalTol, int checkFinalTol, double eta_minSS, double eta_maxSS, comp_d currTime_d, comp_mp currTime_mp)
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
  digits_B2 = ceil(P0 - (order + 1) * eta_minSS / maxNewtonIts);
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
  tempInt = ceil(P0 - (order + 1) * eta_maxSS / maxNewtonIts);

  // find the maximum number of digits that would possibly be needed
  maxDigitsNeeded = MAX(minDigitsNeeded, tempInt);

  // find the precision associated with maxDigitsNeeded and then find the number of digits associated with that precision
  maxPrecNeeded = digits_to_prec(maxDigitsNeeded);
  maxDigitsNeeded = prec_to_digits(maxPrecNeeded);

  // find the best number of digits, precision and the current step size to minimize cost per unit based on minPrecNeeded & maxPrecNeeded
  AMP3_minimize_cost(newDigits, newPrec, currStepSize_d, currStepSize_mp, eta_maxSS, P0, order, maxNewtonIts, minPrecNeeded, maxPrecNeeded);

  return;
}

int AMP3_criterion_error_d(int *consecSuccess, double norm_J, double norm_J_inv, double sizeProportion, int order, double *currStepSize, comp_d currTime, tracker_config_t *T)
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

    // run AMP3 to see what we should do
    AMP3_update(&bestDigits, &bestPrec, currStepSize, stepSize_mp, currPrec, P0, digits_C, order, T->currentNewtonTol, T->maxNewtonIts, T->final_tolerance, T->endgameSwitch, eta_minSS, eta_maxSS, currTime, NULL);
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

int AMP3_criterion_error_mp(int *consecSuccess, int curr_prec, double norm_J, double norm_J_inv, double sizeProportion, int order, mpf_t currStepSize, comp_mp currTime, tracker_config_t *T)
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
  if (currDigits + (order + 1) * eta / T->maxNewtonIts > P0 + extra_digits_before_cut_size && currDigits > digits_C)
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

    retVal = ceil(P0 - (order + 1) * eta / T->maxNewtonIts);
    // verify that digits are actually being added
    if (retVal <= currDigits)
      retVal = currDigits + 5;
  }

  mpf_clear(tempMPF); mpf_clear(newStepSize);

  return retVal;
}


