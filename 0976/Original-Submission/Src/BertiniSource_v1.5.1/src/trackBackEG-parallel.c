// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

void find_minimum_separation(int *numPoints, double *minDist_d, double **startPoint_norm_d, point_d **startPts_d, mpf_t minDist_mp, mpf_t **startPoint_norm_mp, point_mp **startPts_mp, tracker_config_t *T, FILE *START)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: minimum distance & norms of the start points   *
* NOTES: find the minimum separation in the start points        *
\***************************************************************/
{
  int i, j, k, indexI, indexJ, numVars = T->numVars;
  double minDist_sqr_d;

  // find the number of points
  fscanf(START, "%d", numPoints);
  scanRestOfLine(START);

  if (T->MPType == 0 || T->MPType == 2)
  { // sort using double precision
    double tempD, currDist;
    comp_d tempComp;
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(*numPoints * sizeof(sortStruct_d));

    // setup startPoint_norm_d
    *startPoint_norm_d = (double *)bmalloc(*numPoints * sizeof(double));
    // setup startPts_d
    *startPts_d = (point_d *)bmalloc(*numPoints * sizeof(point_d));

    for (i = 0; i < *numPoints; i++)
    { // setup sortPts[i]
      sortPts[i].path_num = i;

      // read in the ith start point and find its norm
      init_point_d((*startPts_d)[i], numVars);
      (*startPts_d)[i]->size = numVars;
      sortPts[i].norm = 0;
      for (j = 0; j < numVars; j++)
      {
        fscanf(START, "%lf%lf", &(*startPts_d)[i]->coord[j].r, &(*startPts_d)[i]->coord[j].i);
        scanRestOfLine(START);

        currDist = norm_sqr_d(&(*startPts_d)[i]->coord[j]);
        if (sortPts[i].norm < currDist)
          sortPts[i].norm = currDist;
      }
      (*startPoint_norm_d)[i] = sortPts[i].norm = sqrt(sortPts[i].norm);
    }

    // sort the structure
    qsort(sortPts, *numPoints, sizeof(sortStruct_d), sort_order_d);

    // do the final comparison
    minDist_sqr_d = *minDist_d = 1e300;
    for (i = 0; i < *numPoints; i++)
    {
      indexI = sortPts[i].path_num;
      j = i + 1;
      while ((j < *numPoints) && (sortPts[j].norm - sortPts[i].norm < *minDist_d))
      { // subtract to see if closer
        currDist = 0;
        indexJ = sortPts[j].path_num;
        for (k = 0; k < numVars; k++)
        { // subtract kth coord & find its norm
          sub_d(tempComp, &(*startPts_d)[indexI]->coord[k], &(*startPts_d)[indexJ]->coord[k]);
          tempD = norm_sqr_d(tempComp);
          if (currDist < tempD)
          { // update currDist
            currDist = tempD;

            // see if we are larger than the minimum - if so, we can move past this one
            if (currDist > minDist_sqr_d)
              k = numVars;
          }
        }
        // update minDist_d
        if (currDist < minDist_sqr_d)
        {
          minDist_sqr_d = currDist;
          *minDist_d = sqrt(minDist_sqr_d);
        }

        // increment j
        j++;
      }
    }

    // release memory
    free(sortPts);
  }
  else
  { // sort using fixed precision
    int cont;
    mpf_t tempMPF, currDist;
    point_mp tempVec;
    sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(*numPoints * sizeof(sortStruct_mp));

    mpf_init(tempMPF);
    mpf_init(currDist);
    init_vec_mp(tempVec, numVars);
    tempVec->size = numVars;

    // setup startPoint_norm_mp
    *startPoint_norm_mp = (mpf_t *)bmalloc(*numPoints * sizeof(mpf_t));
    // setup startPts_mp
    *startPts_mp = (point_mp *)bmalloc(*numPoints * sizeof(point_mp));

    for (i = 0; i < *numPoints; i++)
    { // setup sortPts[i]
      sortPts[i].path_num = i;

      // read in the ith start point and find its norm
      init_point_mp((*startPts_mp)[i], numVars);
      (*startPts_mp)[i]->size = numVars;
      mpf_init(sortPts[i].norm);
      mpf_set_ui(sortPts[i].norm, 0);
      for (j = 0; j < numVars; j++)
      {
        mpf_inp_str((*startPts_mp)[i]->coord[j].r, START, 10);
        mpf_inp_str((*startPts_mp)[i]->coord[j].i, START, 10);
        scanRestOfLine(START);

        norm_sqr_mp(currDist, &(*startPts_mp)[i]->coord[j]);
        if (mpf_cmp(sortPts[i].norm, currDist) < 0)
          mpf_set(sortPts[i].norm, currDist);
      }
      mpf_sqrt(sortPts[i].norm, sortPts[i].norm);

      // setup startPoint_norm_mp
      mpf_init((*startPoint_norm_mp)[i]);
      mpf_set((*startPoint_norm_mp)[i], sortPts[i].norm);
    }

    // sort the structure
    qsort(sortPts, *numPoints, sizeof(sortStruct_mp), sort_order_mp);

    // do the final comparison
    mpf_set_d(minDist_mp, 1e300);
    for (i = 0; i < *numPoints; i++)
    {
      indexI = sortPts[i].path_num;
      cont = 1;
      j = i;
      do
      { // increment the counter - start at i + 1
        j++;

        if (j < *numPoints)
        { // find difference in norm
          mpf_sub(tempMPF, sortPts[j].norm, sortPts[i].norm);
          if (mpf_cmp(tempMPF, minDist_mp) > 0)
            cont = 0;
        }
        else
          cont = 0;

        // check to see if we can continue
        if (cont)
        { // subtract to see if closer
          mpf_set_ui(currDist, 0);
          indexJ = sortPts[j].path_num;
          for (k = 0; k < numVars; k++)
          { // subtract kth coord & find its norm
            sub_mp(&tempVec->coord[k], &(*startPts_mp)[indexI]->coord[k], &(*startPts_mp)[indexJ]->coord[k]);
            norm_sqr_mp(tempMPF, &tempVec->coord[k]);
            mpf_sqrt(tempMPF, tempMPF);
            if (mpf_cmp(currDist, tempMPF) < 0)
            { // update currDist
              mpf_set(currDist, tempMPF);

              // see if we are larger than the minimum - if so, we can move past this one
              if (mpf_cmp(currDist, minDist_mp) > 0)
                k = numVars;
            }
          }
          // update minDist_mp
          if (mpf_cmp(currDist, minDist_mp) < 0)
            mpf_set(minDist_mp, currDist);
        }
      } while (cont);

      // clear norm
      mpf_clear(sortPts[i].norm);
    }

    // release memory
    free(sortPts);
    clear_vec_mp(tempVec);
  }

  rewind(START);

  return;
}

void track_samples_trackBack(int minCycleTrackBack, double trackBack_final_tol, trackBack_samples_t *EGsample, point_data_d *endSamples_d, point_data_mp *endSamples_mp, int samplesPrec, point_d bdryPt_d, point_mp bdryPt_mp, int bdryPrec, point_d startPt_d, point_mp startPt_mp, int startPrec, int cycle, int samples_per_loop, tracker_config_t *T, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: tracks the sample points for 'parallel' endgame        *
\***************************************************************/
{
  int i, retVal, num = 0;
  double currStepSize = T->currentStepSize;
  eval_struct_d e_d;
  eval_struct_mp e_mp;

  init_eval_struct_d(e_d, 0, 0, 0);
  init_eval_struct_mp(e_mp, 0, 0, 0);

  if (T->MPType == 0)
  { // use _d
    int *notUsed = (int *)bmalloc(cycle * sizeof(int));
    double *norm = (double *)bmalloc(cycle * sizeof(double));
    point_d *pts = (point_d *)bmalloc(cycle * sizeof(point_d)), *midPts = (point_d *)bmalloc(cycle * sizeof(point_d));
        
    comp_d time;
    point_data_d bdryPD, startPD;

    init_point_data_d(&bdryPD, 0);
    init_point_data_d(&startPD, 0);
    
    // setup the other ones
    if (cycle >= minCycleTrackBack)
    {
      for (i = 1; i < cycle; i++)
      { // track back to the endgame boundary
        set_double_d(time, T->endgameBoundary, 0);
        T->currentStepSize = currStepSize;
        retVal = track_d(&bdryPD, &endSamples_d[i * samples_per_loop], time, T, OUT, ED_d, eval_func_d);

        // check for success
        if (retVal == 0)
        { // track to the start
          set_one_d(time);
          T->currentNewtonTol = T->basicNewtonTol;
          T->minStepSize = T->minStepSizeBeforeEndGame;
          retVal = track_d(&startPD, &bdryPD, time, T, OUT, ED_d, eval_func_d);

          // check for success
          if (retVal == 0)
          { // refine the value
            refine_d_basic(&startPD, T, OUT, &e_d, ED_d, eval_func_d);

            // store the values
            notUsed[num] = 1;
            init_point_d(pts[num], 0);
            init_point_d(midPts[num], 0);
            norm[num] = point_cp_d_norm(pts[num], startPD.point);
            point_cp_d(midPts[num], bdryPD.point);

            // increment num
            num++;
          }
        }
      }
    }

    // setup values
    EGsample->numSamples = num;
    EGsample->samplePts_prec = 52;
    EGsample->normSamples_d = (double *)brealloc(norm, num * sizeof(double));
    EGsample->normSamples_mp = NULL;
    EGsample->samplePts_d = (point_d *)brealloc(pts, num * sizeof(point_d));
    EGsample->samplePts_mp = NULL;
    EGsample->samplePts_notUsed = (int *)brealloc(notUsed, num * sizeof(int));
    EGsample->midPt_prec = 52;
    EGsample->midPt_d = (point_d *)brealloc(midPts, num * sizeof(point_d));
    EGsample->midPt_mp = NULL;

    // clear memory
    notUsed = NULL;
    norm = NULL;
    pts = NULL;
    midPts = NULL;
    clear_point_data_d(&bdryPD);
    clear_point_data_d(&startPD);
  }
  else if (T->MPType == 1)
  { // use _mp
    int *notUsed = (int *)bmalloc(cycle * sizeof(int));
    mpf_t *norm = (mpf_t *)bmalloc(cycle * sizeof(mpf_t));
    point_mp *pts = (point_mp *)bmalloc(cycle * sizeof(point_mp)), *midPts = (point_mp *)bmalloc(cycle * sizeof(point_mp));

    comp_mp time;
    point_data_mp bdryPD, startPD;

    init_mp(time);
    init_point_data_mp(&bdryPD, 0);
    init_point_data_mp(&startPD, 0);

    // setup the other ones
    if (cycle >= minCycleTrackBack)
    {
      for (i = 1; i < cycle; i++)
      { // track back to the endgame boundary
        set_double_mp(time, T->endgameBoundary, 0);
        T->currentStepSize = currStepSize;
        retVal = track_mp(&bdryPD, &endSamples_mp[i * samples_per_loop], time, T, OUT, ED_mp, eval_func_mp);

        // check for success
        if (retVal == 0)
        { // track to the start
          set_one_mp(time);
          T->currentNewtonTol = T->basicNewtonTol;
          T->minStepSize = T->minStepSizeBeforeEndGame;
          retVal = track_mp(&startPD, &bdryPD, time, T, OUT, ED_mp, eval_func_mp);

          // check for success
          if (retVal == 0)
          { // refine the value
            refine_mp_basic(&startPD, T, OUT, &e_mp, ED_mp, eval_func_mp);

            notUsed[num] = 1;
            mpf_init(norm[num]);
            init_point_mp(pts[num], 0);
            init_point_mp(midPts[num], 0);  
            point_cp_mp_norm(norm[num], pts[num], startPD.point, T->Precision);
            point_cp_mp(midPts[num], bdryPD.point);

            // increment num
            num++;
          }
        }
      }
    }

    // setup values
    EGsample->numSamples = num;
    EGsample->samplePts_prec = T->Precision;
    EGsample->normSamples_d = NULL;
    EGsample->normSamples_mp = (mpf_t *)brealloc(norm, num * sizeof(mpf_t));
    EGsample->samplePts_d = NULL;
    EGsample->samplePts_mp = (point_mp *)brealloc(pts, num * sizeof(point_mp));
    EGsample->samplePts_notUsed = (int *)brealloc(notUsed, num * sizeof(int));
    EGsample->midPt_prec = T->Precision;
    EGsample->midPt_d = NULL;
    EGsample->midPt_mp = (point_mp *)brealloc(midPts, num * sizeof(point_mp));

    // clear memory
    notUsed = NULL;
    norm = NULL;
    pts = NULL;
    midPts = NULL;
    clear_mp(time);
    clear_point_data_mp(&bdryPD);
    clear_point_data_mp(&startPD);
  }
  else
  { // use amp
    int new_startPrec, new_bdryPrec, digitsCorrect = ceil(-log10(trackBack_final_tol) + 3.5), currDigits = ceil(-log10(T->basicNewtonTol));
    int *pt_prec = (int *)bmalloc(cycle * sizeof(int)), *midPt_prec = (int *)bmalloc(cycle * sizeof(int)), *notUsed = (int *)bmalloc(cycle * sizeof(int));
    double first_increase, *norm_d = (double *)bmalloc(cycle * sizeof(double));
    mpf_t *norm_mp = (mpf_t *)bmalloc(cycle * sizeof(mpf_t));
    point_d *pts_d = (point_d *)bmalloc(cycle * sizeof(point_d)), *midPts_d = (point_d *)bmalloc(cycle * sizeof(point_d));
    point_mp *pts_mp = (point_mp *)bmalloc(cycle * sizeof(point_mp)), *midPts_mp = (point_mp *)bmalloc(cycle * sizeof(point_mp));

    comp_d time_d;
    comp_mp time_mp;
    point_data_d bdryPD_d, startPD_d;
    point_data_mp bdryPD_mp, startPD_mp;

    init_mp(time_mp);
    init_point_data_d(&bdryPD_d, 0);
    init_point_data_d(&startPD_d, 0);
    init_point_data_mp(&bdryPD_mp, 0)
    init_point_data_mp(&startPD_mp, 0);

    // setup the other ones
    if (cycle >= minCycleTrackBack)
    {   
      for (i = 1; i < cycle; i++)
      { // track back to the endgame boundary
        set_double_d(time_d, T->endgameBoundary, 0);
        setprec_mp(time_mp, MAX(64, samplesPrec));
        d_to_mp(time_mp, time_d);
        T->currentStepSize = currStepSize;

        if (samplesPrec < 64)
          retVal = AMPtrack(&bdryPD_d, &bdryPD_mp, &new_bdryPrec, &first_increase, &endSamples_d[i * samples_per_loop], NULL, samplesPrec, time_d, time_mp, samplesPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
        else
          retVal = AMPtrack(&bdryPD_d, &bdryPD_mp, &new_bdryPrec, &first_increase, NULL, &endSamples_mp[i * samples_per_loop], samplesPrec, time_d, time_mp, samplesPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

        // check for success
        if (retVal == 0)
        { // track to the start
          setprec_mp(time_mp, MAX(64, new_bdryPrec));
          set_one_d(time_d);
          set_one_mp(time_mp);

          T->currentNewtonTol = T->basicNewtonTol;
          T->minStepSize = T->minStepSizeBeforeEndGame;

          retVal = AMPtrack(&startPD_d, &startPD_mp, &new_startPrec, &first_increase, &bdryPD_d, &bdryPD_mp, new_bdryPrec, time_d, time_mp, new_bdryPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
  
          // check for success
          if (retVal == 0)
          { // refine the value to the correct tolerance
            refine_digits_amp(T->outputLevel, digitsCorrect, &T->latest_newton_residual_d, T->latest_newton_residual_mp, currDigits, &startPD_d, &startPD_mp, &new_startPrec, &startPD_d, &startPD_mp, new_startPrec, time_d, OUT, &e_d, &e_mp, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

            notUsed[num] = 1;
            midPt_prec[num] = new_bdryPrec;
            if (new_bdryPrec < 64)
            {
              init_point_d(midPts_d[num], 0);
              point_cp_d(midPts_d[num], bdryPD_d.point);
            }
            else
            {
              init_point_mp2(midPts_mp[num], 0, new_bdryPrec);
              point_cp_mp(midPts_mp[num], bdryPD_mp.point);
            }
            pt_prec[num] = new_startPrec;
            if (new_startPrec < 64)
            {
              init_point_d(pts_d[num], 0);
              norm_d[num] = point_cp_d_norm(pts_d[num], startPD_d.point);
            }
            else
            {
              mpf_init2(norm_mp[num], new_startPrec);
              init_point_mp2(pts_mp[num], 0, new_startPrec);
              point_cp_mp_norm(norm_mp[num], pts_mp[num], startPD_mp.point, new_startPrec);
            }

            // increment num
            num++;
          }
        }
      }
    }

    // setup the number of samples
    EGsample->numSamples = num;

    // setup midPt
    new_bdryPrec = 52;
    for (i = 0; i < num; i++)
      new_bdryPrec = MAX(new_bdryPrec, midPt_prec[i]);

    if (new_bdryPrec < 64)
    { // all of mid points in double precision
      EGsample->midPt_prec = new_bdryPrec;
      EGsample->midPt_d = (point_d *)brealloc(midPts_d, num * sizeof(point_d));
      EGsample->midPt_mp = NULL;

      // clear memory
      free(midPt_prec);
      free(midPts_mp);
      midPts_d = NULL;
    } 
    else
    { // set the correct precision
      for (i = 0; i < num; i++)
        if (midPt_prec[i] < 64)
        { // move from _d to _mp
          init_point_mp2(midPts_mp[i], 0, new_bdryPrec);
          point_d_to_mp(midPts_mp[i], midPts_d[i]);

          // clear
          clear_point_d(midPts_d[i]);
        }
        else if (midPt_prec[i] < new_bdryPrec)
        { // increase prec
          change_prec_point_mp(midPts_mp[i], new_bdryPrec);
        }

      // setup the values
      EGsample->midPt_prec = new_bdryPrec;
      EGsample->midPt_d = NULL;
      EGsample->midPt_mp = (point_mp *)brealloc(midPts_mp, num * sizeof(point_mp));

      // clear memory
      free(midPt_prec);
      free(midPts_d);
      midPts_mp = NULL;
    }

    // setup pts
    new_startPrec = 52;
    for (i = 0; i < num; i++)
      new_startPrec = MAX(new_startPrec, pt_prec[i]);

    if (new_startPrec < 64)
    { // all in double precision already
      EGsample->samplePts_prec = new_startPrec;
      EGsample->normSamples_d = (double *)brealloc(norm_d, num * sizeof(double));
      EGsample->normSamples_mp = NULL;
      EGsample->samplePts_d = (point_d *)brealloc(pts_d, num * sizeof(point_d));
      EGsample->samplePts_mp = NULL;
      EGsample->samplePts_notUsed = (int *)brealloc(notUsed, num * sizeof(int));

      // clear memory
      free(pt_prec);
      notUsed = NULL;
      norm_d = NULL;
      free(norm_mp);
      pts_d = NULL;
      free(pts_mp);
    }
    else
    { // set the correct precision
      for (i = 0; i < num; i++)
        if (pt_prec[i] < 64)
        { // move from _d to _mp
          mpf_init2(norm_mp[i], new_startPrec);
          init_point_mp2(pts_mp[i], 0, new_startPrec);
          point_d_to_mp_norm(norm_mp[i], pts_mp[i], pts_d[i], new_startPrec);

          // clear
          clear_point_d(pts_d[i]);
        }
        else if (pt_prec[i] < new_startPrec)
        { // increase prec
          mpf_set_prec(norm_mp[i], new_startPrec);
          point_cp_mp_norm2(norm_mp[i], pts_mp[i], new_startPrec);
        }

      // setup the values
      EGsample->samplePts_prec = new_startPrec;
      EGsample->normSamples_d = NULL;
      EGsample->normSamples_mp = (mpf_t *)brealloc(norm_mp, num * sizeof(mpf_t));
      EGsample->samplePts_d = NULL;
      EGsample->samplePts_mp = (point_mp *)brealloc(pts_mp, num * sizeof(point_mp));
      EGsample->samplePts_notUsed = (int *)brealloc(notUsed, num * sizeof(int));

      // clear memory
      free(pt_prec);
      notUsed = NULL;
      free(norm_d);
      norm_mp = NULL;
      free(pts_d);
      pts_mp = NULL;
    }

    // clear memory
    clear_mp(time_mp);
    clear_point_data_d(&bdryPD_d);
    clear_point_data_mp(&bdryPD_mp);
    clear_point_data_d(&startPD_d);
    clear_point_data_mp(&startPD_mp);
  }

  clear_eval_struct_d(e_d);
  clear_eval_struct_mp(e_mp);

  return;
}

void setup_samples_trackBack(int minCycleTrackBack, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d *endSamples_d, point_data_mp *endSamples_mp, int samplesPrec, point_d bdryPt_d, point_mp bdryPt_mp, int bdryPrec, point_d startPt_d, point_mp startPt_mp, int startPrec, int cycle, int samples_per_loop, tracker_config_t *T, double time_first_increase, FILE *OUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup EGsamples from the completed endgame             *
\***************************************************************/
{
  if (EGsamples->endPt.prec < 64)
  { // save using _d
    EGsamples->endPt.latest_newton_residual_d = T->latest_newton_residual_d;
    EGsamples->endPt.t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
    EGsamples->endPt.error_at_latest_sample_point_d = T->error_at_latest_sample_point;
    findFunctionResidual_conditionNumber_d(&EGsamples->endPt.function_residual_d, &EGsamples->endPt.condition_number, &EGsamples->endPt.PD_d, ED_d, eval_func_d);
  }
  else
  { // save using _mp
    mpf_clear(EGsamples->endPt.function_residual_mp);
    mpf_init2(EGsamples->endPt.function_residual_mp, EGsamples->endPt.prec);

    mpf_clear(EGsamples->endPt.latest_newton_residual_mp);
    mpf_init2(EGsamples->endPt.latest_newton_residual_mp, EGsamples->endPt.prec);

    mpf_clear(EGsamples->endPt.t_val_at_latest_sample_point_mp);
    mpf_init2(EGsamples->endPt.t_val_at_latest_sample_point_mp, EGsamples->endPt.prec);

    mpf_clear(EGsamples->endPt.error_at_latest_sample_point_mp);
    mpf_init2(EGsamples->endPt.error_at_latest_sample_point_mp, EGsamples->endPt.prec);

    mpf_set(EGsamples->endPt.latest_newton_residual_mp, T->latest_newton_residual_mp);
    mpf_set_d(EGsamples->endPt.t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
    mpf_set_d(EGsamples->endPt.error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
    findFunctionResidual_conditionNumber_mp(EGsamples->endPt.function_residual_mp, &EGsamples->endPt.condition_number, &EGsamples->endPt.PD_mp, ED_mp, eval_func_mp);
  }

  // setup the sample locations
  track_samples_trackBack(minCycleTrackBack, trackBack_final_tol, EGsamples, endSamples_d, endSamples_mp, samplesPrec, bdryPt_d, bdryPt_mp, bdryPrec, startPt_d, startPt_mp, startPrec, cycle, samples_per_loop, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  return;
}

int trackBackEG_d(int pathNum, int minCycleTrackBack, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d **endSamples, int *cycle, int *samples_per_loop, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the endpoint and then tracks the points back      *
\***************************************************************/
{
  int i, retVal = 0;
  comp_d endTime;
  point_data_d bdryPt;

  // initialize
  init_point_data_d(&bdryPt, 0);
  *cycle = *samples_per_loop = 0;
  *endSamples = NULL;

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_d(endTime, T->endgameBoundary, 0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0;

  // setup last_approx in case of failure
  EGsamples->endPt.last_approx_prec = 52;
  point_cp_d(EGsamples->endPt.last_approx_d, Start->point);

  // track to the endgame boundary
  retVal = track_d(&bdryPt, Start, endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < bdryPt.point->size; i++)
    fprintf(midOUT, "%.15e %.15e\n", bdryPt.point->coord[i].r, bdryPt.point->coord[i].i);
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;
  
  // check for success
  if (retVal)
  { // failure
    point_data_cp_d(&EGsamples->endPt.PD_d, &bdryPt);
    fprintf(OUT, "NOTE: The endgame never started!\n");
    T->t_val_at_latest_sample_point = EGsamples->endPt.PD_d.time->r;
    // NULL out the rest of the structures
    EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
    EGsamples->normSamples_d = NULL;
    EGsamples->normSamples_mp = NULL;
    EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
    EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
    EGsamples->samplePts_notUsed = NULL;
  }
  else
  { // success - so run the acutal endgame
    eval_struct_d e_d;
    init_eval_struct_d(e_d, 0, 0, 0);

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // setup endTime
    set_double_d(endTime, T->targetT, 0);

    // do the endgame in double precision
    EGsamples->endPt.last_approx_prec = 52;
    retVal = CauchyEG_main_d2(&EGsamples->endPt.PD_d, EGsamples->endPt.last_approx_d, endSamples, cycle, samples_per_loop, endTime, &bdryPt, T, OUT, ED, eval_func, find_dehom);

    if (retVal == 0)
    { // setup the appropriate information if successful
      setup_samples_trackBack(minCycleTrackBack, trackBack_final_tol, EGsamples, *endSamples, NULL, 52, bdryPt.point, NULL, 52, Start->point, NULL, 52, *cycle, *samples_per_loop, T, 0, OUT, ED, NULL, eval_func, NULL, NULL);
    }
    else
    { // NULL out the rest of the structures
      EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
      EGsamples->normSamples_d = NULL;
      EGsamples->normSamples_mp = NULL;
      EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
      EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
      EGsamples->samplePts_notUsed = NULL;
    }

    // clear
    clear_eval_struct_d(e_d);
  }

  clear_point_data_d(&bdryPt);

  return retVal;
}

//// Multi precision ////

int trackBackEG_mp(int pathNum, int minCycleTrackBack, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_mp **endSamples, int *cycle, int *samples_per_loop, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the endpoint and then tracks the points back      *
\***************************************************************/
{
  int i, retVal = 0;
  comp_mp endTime;
  point_data_mp bdryPt;

  init_point_data_mp(&bdryPt, 0);
  *cycle = *samples_per_loop = 0;
  *endSamples = NULL;

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  init_mp(endTime);
  set_double_mp(endTime, T->endgameBoundary, 0);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0;

  // setup last_approx in case of failure
  EGsamples->endPt.last_approx_prec = T->Precision;
  point_cp_mp(EGsamples->endPt.last_approx_mp, Start->point);

  // track to the endgame boundary
  retVal = track_mp(&bdryPt, Start, endTime, T, OUT, ED, eval_func);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  for (i = 0; i < bdryPt.point->size; i++)
  {
    print_mp(midOUT, 0, &bdryPt.point->coord[i]);
    fprintf(midOUT, "\n");
  }
  fprintf(midOUT, "\n");

  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check for success
  if (retVal)
  { // failure
    point_data_cp_mp(&EGsamples->endPt.PD_mp, &bdryPt);
    fprintf(OUT, "NOTE: The endgame never started!\n");
    T->t_val_at_latest_sample_point = mpf_get_d(EGsamples->endPt.PD_mp.time->r);
    // NULL out the rest of the structures
    EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
    EGsamples->normSamples_d = NULL;
    EGsamples->normSamples_mp = NULL;
    EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
    EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
    EGsamples->samplePts_notUsed = NULL;
  }
  else
  { // success - so run the acutal endgame
    eval_struct_mp e_mp;
    init_eval_struct_mp(e_mp, 0, 0, 0);

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // setup endTime
    set_double_mp(endTime, T->targetT, 0);

    // do the endgame in fixed precision
    EGsamples->endPt.last_approx_prec = T->Precision;
    retVal = CauchyEG_main_mp2(&EGsamples->endPt.PD_mp, EGsamples->endPt.last_approx_mp, endSamples, cycle, samples_per_loop, endTime, &bdryPt, T, OUT, ED, eval_func, find_dehom);

    if (retVal == 0)
    { // setup the appropriate information if successful
      setup_samples_trackBack(minCycleTrackBack, trackBack_final_tol, EGsamples, NULL, *endSamples, T->Precision, NULL, bdryPt.point, T->Precision, NULL, Start->point, T->Precision, *cycle, *samples_per_loop, T, 0, OUT, NULL, ED, NULL, eval_func, NULL);
    }
    else
    { // NULL out the rest of the structures
      EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
      EGsamples->normSamples_d = NULL;
      EGsamples->normSamples_mp = NULL;
      EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
      EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
      EGsamples->samplePts_notUsed = NULL;
    }

    // clear
    clear_eval_struct_mp(e_mp);
  }

  clear_mp(endTime);
  clear_point_data_mp(&bdryPt);

  return retVal;
}

/// Adaptive precision ///

int trackBackEG_amp(int pathNum, int minCycleTrackBack, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d **endSamples_d, point_data_mp **endSamples_mp, int *samples_prec, int *cycle, int *samples_per_loop, point_data_d *Start_d, point_data_mp *Start_mp, int startPrec, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find the endpoint and then tracks the points back      *
\***************************************************************/
{
  int i, bdryPrec, retVal = 0;
  comp_d endTime_d;
  comp_mp endTime_mp;
  point_data_d bdryPt_d;
  point_data_mp bdryPt_mp;

  // set the precision
  T->Precision = MAX(64, startPrec);
  initMP(T->Precision);
  change_prec(ED_mp, T->Precision);

  // initialize memory
  init_mp(endTime_mp);
  init_point_data_d(&bdryPt_d, 0);
  init_point_data_mp(&bdryPt_mp, 0);
  *cycle = *samples_per_loop = 0;
  *endSamples_d = NULL;
  *endSamples_mp = NULL;
  *samples_prec = 52;

  // first, we need to track to the endgameBoundary - it is assumed the T->endgameBoundary >= 0 == T->targetT
  set_double_d(endTime_d, T->endgameBoundary, 0);
  d_to_mp(endTime_mp, endTime_d);
  T->currentStepSize = T->maxStepSize;
  T->currentNewtonTol = T->basicNewtonTol;
  T->minStepSize = T->minStepSizeBeforeEndGame;
  T->endgameSwitch = 0;
  T->error_at_latest_sample_point = 0;

  // setup last_approx in case of failure
  EGsamples->endPt.last_approx_prec = startPrec;
  if (startPrec < 64)
  {
    point_cp_d(EGsamples->endPt.last_approx_d, Start_d->point);
  }
  else
  {
    setprec_point_mp(EGsamples->endPt.last_approx_mp, startPrec);
    point_cp_mp(EGsamples->endPt.last_approx_mp, Start_mp->point);
  }

  // track to the endgame boundary
  retVal = AMPtrack(&bdryPt_d, &bdryPt_mp, &bdryPrec, &EGsamples->endPt.first_increase, Start_d, Start_mp, startPrec, endTime_d, endTime_mp, startPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);

  // print the ending point to midOUT
  fprintf(midOUT, "%d\n", pathNum);
  if (bdryPrec == 52)
  { // print to midOUT using _d
    for (i = 0; i < bdryPt_d.point->size; i++)
      fprintf(midOUT, "%.15e %.15e\n", bdryPt_d.point->coord[i].r, bdryPt_d.point->coord[i].i);
    fprintf(midOUT, "\n");
  }
  else
  { // print to midOUT using _mp
    for (i = 0; i < bdryPt_mp.point->size; i++)
    {
      print_mp(midOUT, 0, &bdryPt_mp.point->coord[i]);
      fprintf(midOUT, "\n");
    }
    fprintf(midOUT, "\n");
  }
  // since we have printed to midOUT, this path is considered to be "in the endgame"
  T->endgameSwitch = 1;

  // check for success
  if (retVal)
  { // failure
    EGsamples->endPt.prec = bdryPrec;
    if (bdryPrec < 64)
    { // setup Final_d
      point_data_cp_d(&EGsamples->endPt.PD_d, &bdryPt_d);
      T->t_val_at_latest_sample_point = EGsamples->endPt.PD_d.time->r;
    }
    else
    { // setup Final_mp
      setprec_point_data_mp(&EGsamples->endPt.PD_mp, bdryPrec);
      point_data_cp_mp(&EGsamples->endPt.PD_mp, &bdryPt_mp);
      T->t_val_at_latest_sample_point = mpf_get_d(EGsamples->endPt.PD_mp.time->r);
    }
    // NULL out the rest of the structures
    EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
    EGsamples->normSamples_d = NULL;
    EGsamples->normSamples_mp = NULL;
    EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
    EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
    EGsamples->samplePts_notUsed = NULL; 
    fprintf(OUT, "NOTE: The endgame never started!\n");
  }
  else
  { // success - so run the acutal endgame
    eval_struct_d e_d;
    eval_struct_mp e_mp;

    init_eval_struct_d(e_d, 0, 0, 0);
    init_eval_struct_mp(e_mp, 0, 0, 0);

    // set the endgame tolerances
    T->currentNewtonTol = T->endgameNewtonTol;
    T->minStepSize = T->minStepSizeDuringEndGame;

    // setup endTime
    set_double_d(endTime_d, T->targetT, 0);
    d_to_mp(endTime_mp, endTime_d);

    // do the endgame
    retVal = CauchyEG_main_amp2(&EGsamples->endPt.prec, &EGsamples->endPt.PD_d, &EGsamples->endPt.PD_mp, EGsamples->endPt.last_approx_d, EGsamples->endPt.last_approx_mp, &EGsamples->endPt.last_approx_prec, endSamples_d, endSamples_mp, cycle, samples_per_loop, &EGsamples->endPt.first_increase, endTime_d, endTime_mp, &bdryPt_d, &bdryPt_mp, bdryPrec, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
    *samples_prec = EGsamples->endPt.prec;

    if (retVal == 0)
    { // setup the appropriate information if successful
      setup_samples_trackBack(minCycleTrackBack, trackBack_final_tol, EGsamples, *endSamples_d, *endSamples_mp, EGsamples->endPt.prec, bdryPt_d.point, bdryPt_mp.point, bdryPrec, Start_d->point, Start_mp->point, startPrec, *cycle, *samples_per_loop, T, 0, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }
    else
    { // NULL out the rest of the structures
      EGsamples->numSamples = EGsamples->samplePts_prec = EGsamples->midPt_prec = 0;
      EGsamples->normSamples_d = NULL;
      EGsamples->normSamples_mp = NULL;
      EGsamples->samplePts_d = EGsamples->midPt_d = NULL;
      EGsamples->samplePts_mp = EGsamples->midPt_mp = NULL;
      EGsamples->samplePts_notUsed = NULL;
    }

    // clear
    clear_eval_struct_d(e_d);
    clear_eval_struct_mp(e_mp);
  }

  clear_mp(endTime_mp);
  clear_point_data_d(&bdryPt_d);
  clear_point_data_mp(&bdryPt_mp);

  return retVal;
}

//// zero dim track back endgame ////

int check_point_trackBack(int *indexI, int *indexJ, double tol, trackBack_samples_t **EGsamples, int numSamples, double norm_d, mpf_t norm_mp, point_d testPt_d, point_mp testPt_mp, int testPt_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 0 - not found, 1 - found                       *
* NOTES: check to see if the start point has been found         *
\***************************************************************/
{
  int i, j, num, retVal = 0;
  trackBack_samples_t *testEG = NULL;

  if (testPt_prec < 64)
  { // compare using double precision
    for (i = 0; i < numSamples; i++)
    { // loop over the sample points
      testEG = &(*EGsamples)[i];
      if (testEG->samplePts_prec < 64)
      { // sample points in double precision
        num = testEG->numSamples;
        for (j = 0; j < num; j++)
        { // determine if this one has not been used and the difference in norm is small enough
          if (testEG->samplePts_notUsed[j] && fabs(testEG->normSamples_d[j] - norm_d) < tol)
          { // need to compare the values
            if (isSamePoint(testEG->samplePts_d[j], NULL, testEG->samplePts_prec, testPt_d, testPt_mp, testPt_prec, tol))
            { // we have a match!
              retVal = 1;
              *indexI = i;
              *indexJ = j;
              // break out of loops
              j = testEG->numSamples;
              i = numSamples;
            }
          }
        }
      }
      else
      { // sample points in fixed precision
        num = testEG->numSamples;
        for (j = 0; j < num; j++)
        { // determine if this one has not been used and the difference in norm is small enough
          if (testEG->samplePts_notUsed[j] && fabs(mpf_get_d(testEG->normSamples_mp[j]) - norm_d) < tol)
          { // need to compare the values
            if (isSamePoint(NULL, testEG->samplePts_mp[j], testEG->samplePts_prec, testPt_d, testPt_mp, testPt_prec, tol))
            { // we have a match!
              retVal = 1;
              *indexI = i;
              *indexJ = j;
              // break out of loops
              j = testEG->numSamples;
              i = numSamples;
            }
          }
        }
      }
    }
  }
  else
  { // compare using multi precision
    for (i = 0; i < numSamples; i++)
    { // loop over the sample points
      testEG = &(*EGsamples)[i];
      if (testEG->samplePts_prec < 64)
      { // sample points in double precision
        num = testEG->numSamples;
        for (j = 0; j < num; j++)
        { // determine if this one has not been used and the difference in norm is small enough
          if (testEG->samplePts_notUsed[j] && fabs(testEG->normSamples_d[j] - mpf_get_d(norm_mp)) < tol)
          { // need to compare the values
            if (isSamePoint(testEG->samplePts_d[j], NULL, testEG->samplePts_prec, testPt_d, testPt_mp, testPt_prec, tol))
            { // we have a match!
              retVal = 1;
              *indexI = i;
              *indexJ = j;
              // break out of loops
              j = testEG->numSamples;
              i = numSamples;
            }
          }
        }
      }
      else
      { // sample points in fixed precision
        mpf_t diff;
        mpf_init2(diff, testPt_prec);

        num = testEG->numSamples;
        for (j = 0; j < num; j++)
        { // determine if this one has not been used and the difference in norm is small enough
          if (testEG->samplePts_notUsed[j])
          { // compare the norms
            mpf_sub(diff, testEG->normSamples_mp[j], norm_mp);
            if (mpf_get_d(diff) < tol)          
            { // need to compare the values
              if (isSamePoint(NULL, testEG->samplePts_mp[j], testEG->samplePts_prec, testPt_d, testPt_mp, testPt_prec, tol))
              { // we have a match!
                retVal = 1;
                *indexI = i;
                *indexJ = j;
                // break out of loops
                j = testEG->numSamples;
                i = numSamples;
              }
            }
          }
        }

        mpf_clear(diff);
      }
    }
  }

  testEG = NULL;

  return retVal;
}

void zero_dim_trackBack_path_d(int pathNum, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame in double precision        *
\***************************************************************/
{
  int i, cycle = 0, samples_per_loop = 0;

  EGsamples->endPt.pathNum = pathNum;
  EGsamples->endPt.codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;
  if (T->MPType == 2)
  { // track using AMP
    int samples_prec = 52;
    point_data_d *endSamples_d = NULL;
    point_data_mp *endSamples_mp = NULL;

    EGsamples->endPt.retVal = trackBackEG_amp(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples_d, &endSamples_mp, &samples_prec, &cycle, &samples_per_loop, Pin, NULL, 52, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    if (samples_prec < 64)
    {
      if (endSamples_d != NULL)
      {
        for (i = samples_per_loop * cycle - 1; i >= 0; i--)
          clear_point_data_d(&endSamples_d[i]);
        free(endSamples_d);
      }
    }
    else
    {
      if (endSamples_mp != NULL)
      {
        for (i = samples_per_loop * cycle - 1; i >= 0; i--)
          clear_point_data_mp(&endSamples_mp[i]);
        free(endSamples_mp);
      }
    }
  }
  else
  { // track using double precision
    point_data_d *endSamples = NULL;

    EGsamples->endPt.prec = EGsamples->endPt.last_approx_prec = 52;
    EGsamples->endPt.first_increase = 0;

    EGsamples->endPt.retVal = trackBackEG_d(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples, &cycle, &samples_per_loop, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);

    // clear
    if (endSamples != NULL)
    {
      for (i = samples_per_loop * cycle - 1; i >= 0; i--)
        clear_point_data_d(&endSamples[i]);
      free(endSamples);
    }
  }

  return;
}

void zero_dim_trackBack_d(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame for the start points       *
\***************************************************************/
{
  int i, j, size, rV, indexI, indexJ, numPoints, same_pathNum, numSamples = 0;
  double trackBack_final_tol, minDist, *startPoint_norm = NULL;
  point_d *startPts = NULL;
  point_data_d startPD;
  trackBack_samples_t *EGsamples = NULL;
  FILE *NONSOLN = NULL;

  init_point_data_d(&startPD, 0);

  // setup NONSOLN
  if (!((basic_eval_data_d *)ED_d)->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // setup the start points
  find_minimum_separation(&numPoints, &minDist, &startPoint_norm, &startPts, NULL, NULL, NULL, T, START);
  trackCount->numPoints = numPoints;

  // setup the track back final tol based on the minimum distance
  trackBack_final_tol = 1e-2 * minDist;

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // track each of the start points
  for (i = 0; i < numPoints; i++)  
  { // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, numPoints);

    // setup startPD and clear startPts[i]
    point_cp_d(startPD.point, startPts[i]);
    set_one_d(startPD.time);
    clear_point_d(startPts[i]);

    // print the header of the path to OUT
    printPathHeader_d(OUT, &startPD, T, i, ED_d, eval_func_d);

    // determine if the start point is already known & find the indices if it is found
    rV = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, startPoint_norm[i], NULL, startPD.point, NULL, 52);
    if (rV)
    { // the point was found
      // print the mid point to MIDOUT
      fprintf(MIDOUT, "%d\n", i);
      if (EGsamples[indexI].midPt_prec < 64)
      { // print midPt_d
        size = EGsamples[indexI].midPt_d[indexJ]->size;
        for (j = 0; j < size; j++)
          fprintf(MIDOUT, "%.15e %.15e\n", EGsamples[indexI].midPt_d[indexJ]->coord[j].r, EGsamples[indexI].midPt_d[indexJ]->coord[j].i);
        fprintf(MIDOUT, "\n");
      }
      else
      { // print midPt_mp
        size = EGsamples[indexI].midPt_mp[indexJ]->size;
        for (j = 0; j < size; j++)
        {
          print_mp(MIDOUT, 0, &EGsamples[indexI].midPt_mp[indexJ]->coord[j]);
          fprintf(MIDOUT, "\n");
        }
        fprintf(MIDOUT, "\n");
      }

      // store the path number & update it
      same_pathNum = EGsamples[indexI].endPt.pathNum;
      EGsamples[indexI].endPt.pathNum = i;  

      // print the footer
      if (EGsamples[indexI].endPt.prec < 64)
      { // print footer in double precision
        printPathFooter_d(trackCount, &EGsamples[indexI].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
      }
      else
      { // print footer in multi precision
        printPathFooter_mp(trackCount, &EGsamples[indexI].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
      }

      // send the path number back
      EGsamples[indexI].endPt.pathNum = same_pathNum;

      // say that this sample has been used
      EGsamples[indexI].samplePts_notUsed[indexJ] = 0;
    }
    else
    { // increase the size of EGsamples
      EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));

      if (T->MPType == 0)
      { // initialize for double precision tracking
        init_endgame_data(&EGsamples[numSamples].endPt, 52);
      }
      else
      { // initialize for AMP tracking
        init_endgame_data(&EGsamples[numSamples].endPt, 64);
      }

      // track the path
      zero_dim_trackBack_path_d(i, trackBack_final_tol, &EGsamples[numSamples], &startPD, OUT, MIDOUT, T, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

      // check to see if it should be sharpened
      if (EGsamples[numSamples].endPt.retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&EGsamples[numSamples].endPt, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
      }

      // print the footer
      if (EGsamples[numSamples].endPt.prec < 64)
      { // print footer in double precision
        printPathFooter_d(trackCount, &EGsamples[numSamples].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_d);
      }
      else
      { // print footer in multi precision
        printPathFooter_mp(trackCount, &EGsamples[numSamples].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);
      }

      // increment the size
      numSamples++;
    }
  }

  if (!((basic_eval_data_d *)ED_d)->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }

  // clear memory
  for (i = numSamples - 1; i >= 0; i--)
  { // clear EGsamples[i]
    clear_trackBack_sample(&EGsamples[i]);
  }
  free(startPoint_norm);
  free(startPts);
  clear_point_data_d(&startPD);

  return;
}

void zero_dim_trackBack_path_mp(int pathNum, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame in fixed multi precision   *
\***************************************************************/
{
  int i, cycle = 0, samples_per_loop = 0;
  point_data_mp *endSamples = NULL;

  EGsamples->endPt.pathNum = pathNum;
  EGsamples->endPt.codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;

  EGsamples->endPt.prec = T->Precision;
  EGsamples->endPt.first_increase = 0;

  // track using fixed multi precision
  EGsamples->endPt.retVal = trackBackEG_mp(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples, &cycle, &samples_per_loop, Pin, T, OUT, MIDOUT, ED, eval_func_mp, find_dehom);

  // clear
  if (endSamples != NULL)
  {
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_mp(&endSamples[i]);
    free(endSamples);
  }

  return;
}

void zero_dim_trackBack_mp(trackingStats *trackCount, FILE *OUT, FILE *RAWOUT, FILE *MIDOUT, FILE *START, FILE *FAIL, int pathMod, tracker_config_t *T, void const *ED_mp, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame for the start points       *
\***************************************************************/
{
  int i, j, size, rV, indexI, indexJ, numPoints, same_pathNum, numSamples = 0;
  double trackBack_final_tol;
  mpf_t minDist, *startPoint_norm = NULL;
  point_mp *startPts = NULL;
  point_data_mp startPD;
  trackBack_samples_t *EGsamples = NULL;
  FILE *NONSOLN = NULL;

  mpf_init(minDist);
  init_point_data_mp(&startPD, 0);

  // setup NONSOLN
  if (!((basic_eval_data_mp *)ED_mp)->squareSystem.noChanges)
  { // setup NONSOLN
    NONSOLN = fopen("nonsolutions", "w");
    fprintf(NONSOLN, "                                    \n\n");
  }

  // setup the start points
  find_minimum_separation(&numPoints, NULL, NULL, NULL, minDist, &startPoint_norm, &startPts, T, START);
  trackCount->numPoints = numPoints;

  // setup the track back final tol based on the minimum distance
  trackBack_final_tol = 1e-2 * mpf_get_d(minDist);

  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);

  // track each of the start points
  for (i = 0; i < numPoints; i++)
  { // print the path number if needed
    if (pathMod > 0 && !(i % pathMod))
      printf("Tracking path %d of %d\n", i, numPoints);

    // setup startPD and clear startPts[i]
    point_cp_mp(startPD.point, startPts[i]);
    set_one_mp(startPD.time);
    clear_point_mp(startPts[i]);

    // print the header of the path to OUT
    printPathHeader_mp(OUT, &startPD, T, i, ED_mp, eval_func_mp);

    // determine if the start point is already known & find the indices if it is found
    rV = check_point_trackBack(&indexI, &indexJ, trackBack_final_tol, &EGsamples, numSamples, 0, startPoint_norm[i], NULL, startPD.point, T->Precision);
    if (rV)
    { // the point was found

      // print the mid point to MIDOUT
      size = EGsamples[indexI].midPt_mp[indexJ]->size;
      fprintf(MIDOUT, "%d\n", i);
      for (j = 0; j < size; j++)
      {
        print_mp(MIDOUT, 0, &EGsamples[indexI].midPt_mp[indexJ]->coord[j]);
        fprintf(MIDOUT, "\n");
      }
      fprintf(MIDOUT, "\n");

      // store the path number & update it
      same_pathNum = EGsamples[indexI].endPt.pathNum;
      EGsamples[indexI].endPt.pathNum = i;

      // print the footer
      printPathFooter_mp(trackCount, &EGsamples[indexI].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);

      // send the path number back
      EGsamples[indexI].endPt.pathNum = same_pathNum;

      // say that this sample has been used
      EGsamples[indexI].samplePts_notUsed[indexJ] = 0;
    }
    else
    { // increase the size of EGsamples
      EGsamples = (trackBack_samples_t *)brealloc(EGsamples, (numSamples + 1) * sizeof(trackBack_samples_t));

      init_endgame_data(&EGsamples[numSamples].endPt, T->Precision);

      // track the path
      zero_dim_trackBack_path_mp(i, trackBack_final_tol, &EGsamples[numSamples], &startPD, OUT, MIDOUT, T, ED_mp, eval_func_mp, find_dehom);

      // check to see if it should be sharpened
      if (EGsamples[numSamples].endPt.retVal == 0 && T->sharpenDigits > 0)
      { // use the sharpener for after an endgame
        sharpen_endpoint_endgame(&EGsamples[numSamples].endPt, T, OUT, NULL, ED_mp, NULL, eval_func_mp, NULL);
      }

      // print the footer
      printPathFooter_mp(trackCount, &EGsamples[numSamples].endPt, T, OUT, RAWOUT, FAIL, NONSOLN, ED_mp);

      // increment the size
      numSamples++;
    }
  }

  if (!((basic_eval_data_mp *)ED_mp)->squareSystem.noChanges)
  { // complete NONSOLN
    rewind(NONSOLN);
    fprintf(NONSOLN, "%d", trackCount->junkCount);
    fclose(NONSOLN);
  }
  // clear memory
  for (i = numSamples - 1; i >= 0; i--)
  { // clear EGsamples[i]
    clear_trackBack_sample(&EGsamples[i]);
  }
  free(startPoint_norm);
  free(startPts);
  mpf_clear(minDist);
  clear_point_data_mp(&startPD);

  return;
}

void zero_dim_trackBack_path_rank_d(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_d *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame in double precision        *
\***************************************************************/
{
  int i, cycle = 0, samples_per_loop = 0;

  EGsamples->endPt.pathNum = pathNum;
  EGsamples->endPt.codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;
  if (T->MPType == 2)
  { // track using AMP
    int samples_prec = 52;
    comp_d finalTime_d;
    comp_mp finalTime_mp;
    point_data_d *endSamples_d = NULL;
    point_data_mp *endSamples_mp = NULL;
  
    init_mp(finalTime_mp);
    set_double_d(finalTime_d, T->targetT, 0);
    d_to_mp(finalTime_mp, finalTime_d);

    // track using AMP
    EGsamples->endPt.retVal = trackBackEG_amp(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples_d, &endSamples_mp, &samples_prec, &cycle, &samples_per_loop, Pin, NULL, 52, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);

    if (EGsamples->endPt.retVal == 0 || EGsamples->endPt.retVal == retVal_EG_failed_to_converge)
    { // continue on with rank determination
      EGsamples->endPt.retVal = Cauchy_rank_main_amp(&EGsamples->endPt.condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, EGsamples->endPt.retVal, &EGsamples->endPt.prec, &EGsamples->endPt.first_increase, &EGsamples->endPt.PD_d, &EGsamples->endPt.PD_mp, EGsamples->endPt.last_approx_d, EGsamples->endPt.last_approx_mp, EGsamples->endPt.last_approx_prec, finalTime_d, finalTime_mp, &endSamples_d[0], &endSamples_mp[0], samples_prec, cycle, samples_per_loop, T, OUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec);
    }
    else
    { // failure
      EGsamples->endPt.condition_number = -1;
      *smallest_nonzero_SV = *largest_zero_SV = 0;
    } 

    if (samples_prec < 64)
    {
      if (endSamples_d != NULL)
      {
        for (i = samples_per_loop * cycle - 1; i >= 0; i--)
          clear_point_data_d(&endSamples_d[i]);
        free(endSamples_d);
      }
    }
    else
    {
      if (endSamples_mp != NULL)
      {
        for (i = samples_per_loop * cycle - 1; i >= 0; i--)
          clear_point_data_mp(&endSamples_mp[i]);
        free(endSamples_mp);
      }
    }
    clear_mp(finalTime_mp);
  }
  else
  { // set last_approx_d
    comp_d finalTime;
    point_data_d *endSamples = NULL;

    set_double_d(finalTime, T->targetT, 0);

    // track using double precision
    EGsamples->endPt.prec = EGsamples->endPt.last_approx_prec = 52;
    EGsamples->endPt.first_increase = 0;

    EGsamples->endPt.retVal = trackBackEG_d(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples, &cycle, &samples_per_loop, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);

    if (EGsamples->endPt.retVal == 0 || EGsamples->endPt.retVal == retVal_EG_failed_to_converge)
    { // continue on with rank determination
      EGsamples->endPt.retVal = Cauchy_rank_main_d(&EGsamples->endPt.condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, EGsamples->endPt.retVal, &EGsamples->endPt.PD_d, EGsamples->endPt.last_approx_d, finalTime, &endSamples[0], cycle, samples_per_loop, T, OUT, ED_d, eval_func_d);
    }

    // clear
    if (endSamples != NULL)
    {
      for (i = samples_per_loop * cycle - 1; i >= 0; i--)
        clear_point_data_d(&endSamples[i]);
      free(endSamples);
    }
  }

  return;
}

void zero_dim_trackBack_path_rank_mp(int pathNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, double trackBack_final_tol, trackBack_samples_t *EGsamples, point_data_mp *Pin, FILE *OUT, FILE *MIDOUT, tracker_config_t *T, void const *ED, int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does the 'parallel' endgame in fixed multi precision   *
\***************************************************************/
{
  int i, cycle = 0, samples_per_loop = 0;
  comp_mp finalTime;
  point_data_mp *endSamples = NULL;

  init_mp(finalTime);
  set_double_mp(finalTime, T->targetT, 0);

  EGsamples->endPt.pathNum = pathNum;
  EGsamples->endPt.codim = 0; // zero dimensional - this is ignored

  T->first_step_of_path = 1;

  EGsamples->endPt.prec = EGsamples->endPt.last_approx_prec = T->Precision;
  EGsamples->endPt.first_increase = 0;

  // track using fixed multi precision
  EGsamples->endPt.retVal = trackBackEG_mp(pathNum, T->minCycleTrackBack, trackBack_final_tol, EGsamples, &endSamples, &cycle, &samples_per_loop, Pin, T, OUT, MIDOUT, ED, eval_func_mp, find_dehom);

  if (EGsamples->endPt.retVal == 0 || EGsamples->endPt.retVal == retVal_EG_failed_to_converge)
  { // continue on with rank determination
    EGsamples->endPt.retVal = Cauchy_rank_main_mp(&EGsamples->endPt.condition_number, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, EGsamples->endPt.retVal, &EGsamples->endPt.PD_mp, EGsamples->endPt.last_approx_mp, finalTime, &endSamples[0], cycle, samples_per_loop, T, OUT, ED, eval_func_mp);
  }

  // clear
  if (endSamples != NULL)
  {
    for (i = samples_per_loop * cycle - 1; i >= 0; i--)
      clear_point_data_mp(&endSamples[i]);
    free(endSamples);
  }

  clear_mp(finalTime);

  return;
}

