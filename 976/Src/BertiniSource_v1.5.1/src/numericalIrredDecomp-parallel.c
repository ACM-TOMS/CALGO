// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"
#include "regen_pos_dim.h"

void numIrredDecompMainData(char *mainName, witness_t *W, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom);

void clear_sliceMover(membership_slice_moving_t *sliceMover, int MPType);

int junkRemoval_mem(witness_t *W, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode);

void numericalIrredDecomp(unsigned int currentSeed, int MPType, int genType, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does numerical irreducible decomposition using either  *
* cascade or dim-by-dim to generate the witness superset        *
\***************************************************************/
// 3 step process: 
// 1. Generate witness superset (genType == 0 - use cascade, otherwise, use dim-by-dim)
// 2. Do junk removal
// 3. Do the break-up of the pure dimensional witness sets into their irreducible components
{
  int i, userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0, topDimension = 0;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  tracker_config_t T;
  witness_t witnessSet;
  int **fullRankProgInfo = NULL; // 0 - just needs NULLed out, 1 - deflated properly and needs to be cleared, -1 - did not deflate properly but needs to be cleared
  prog_t ***fullRankProgs = NULL;
  endpoint_data_d **endPts_d = NULL;
  endpoint_data_mp **endPts_mp = NULL;
  endpoint_data_amp **endPts_amp = NULL;
  membership_slice_moving_t *sliceMover = NULL;

  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &constructWitnessSet, &supersetOnly, &paramHom, MPType);

  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision

#ifdef _OPENMP
  #pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }

  // if we are using MPI - we need to tell the workers that they are doing irreducible decomposition
  if (num_processes > 1)
  { // broadcast numerical irreducible decomposition 
#ifdef _HAVE_MPI
    worker_info sendType;
    sendType.dataType = IRREDDECOMP;
    bcast_worker_info(&sendType, my_id, headnode);
#endif
  }

  // use cascade/dim-by-dim/regeneration to generate the witness superset
  if (genType == 2 || useRegen)
  { // use regeneration
    topDimension = regen_pos_dim_main(&witnessSet, regenStartLevel, maxCodim, specificCodim, &T, pathMod, midpoint_tol, intrinsicCutoffMultiplier, my_id, num_processes, headnode);
  }
  else if (genType == 1)
  { // use dim-by-dim
    topDimension = 0; // unknown if we are solving for the top dimension
    dimbydim_main(&witnessSet, maxCodim, specificCodim, &T, pathMod, midpoint_tol, intrinsicCutoffMultiplier, my_id, num_processes, headnode);
  }
  else
  { // use cascade
    topDimension = cascade_main(&witnessSet, maxCodim, specificCodim, &T, pathMod, midpoint_tol, intrinsicCutoffMultiplier, my_id, num_processes, headnode);
  }

  // remove all of the extra multiple points
  for (i = 0; i < witnessSet.num_codim; i++)
  {
    multiplicity_witness(&witnessSet, i, T.MPType, T.final_tol_times_mult);
  }

  // setup to do everything extrinsic (use original variables)
  setupWitnessTotallyExtrinisic(&witnessSet, T.MPType, T.AMP_max_prec);

  // print a witness superset file
  witnessSupersetOutput(&witnessSet, T.MPType);

  // if we are using MPI - we need to tell the workers if they are doing decompositions
  if (num_processes > 1)
  { // broadcast numerical irreducible decomposition 
#ifdef _HAVE_MPI
    worker_info sendType;
    if (supersetOnly)
      sendType.dataType = STOPCODE;
    else
      sendType.dataType = IRREDDECOMP;
    bcast_worker_info(&sendType, my_id, headnode);
#endif
  }

  // determine if we should do the decomposition
  if (!supersetOnly)
  { // setup the collection of fullRankProgs and do junk removal (turn witness superset into witness set)
    junkRemoval(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, &T, pathMod, midpoint_tol, reducedOnly, specificCodim, topDimension, my_id, num_processes, headnode);

    // do the break-up into irreducible components
    pureDecomp(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, &witnessSet, &T, pathMod, my_id, num_processes, headnode);

    // display decomposition chart
    numIrredDecompChart(&witnessSet, stdout, T.MPType, reducedOnly);

    // create output files
    numIrredDecompOutput(&witnessSet, &T, 1, genType, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 1

    // display deflation errors
    displayDeflationSummary(fullRankProgInfo, &witnessSet);

    // clear the deflation information
    clear_sliceMover_fullRankProgs(&sliceMover, &fullRankProgs, &fullRankProgInfo, &endPts_d, &endPts_mp, &endPts_amp, &witnessSet, T.MPType);
  }

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();
 
  return;
}

void junkRemoval(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup fullRankProgs and endPts and then remove junk    *
* from W to make it a true witness set                          *
\***************************************************************/
{
  if (num_processes > 1)
  {
#ifdef _HAVE_MPI
    junkRemoval_par(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, pathMod, midpoint_tol, reducedOnly, specificCodim, topDimension, my_id, num_processes, headnode);
#endif
  }
  else
  {
    junkRemoval_seq(sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, W, T, pathMod, midpoint_tol, reducedOnly, specificCodim, topDimension, my_id, num_processes, headnode);
  }

  return;
}

void junkRemoval_seq(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, tracker_config_t *T, int pathMod, double midpoint_tol, int reducedOnly, int specificCodim, int topDimension, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup fullRankProgs and endPts and then remove junk    *
* from W to make it a true witness set                          *
\***************************************************************/
{
  int i, j, k, num_sing, num_codim = W->num_codim, *isJunk = NULL, *singEndpoints = NULL;
  FILE *OUT = NULL;
  char outName[] = "output_junkRemoval", midName[] = "midpath_data";

  // make sure that we should continue
  if (num_codim < 1)
  { // no codimensions for junk removal!
    return;
  }

  // seutp OUT
  OUT = fopen(outName, "w");

  // allocate for each codimension
  *sliceMover = (membership_slice_moving_t *)bmalloc(num_codim * sizeof(membership_slice_moving_t));
  *fullRankProgs = (prog_t ***)bmalloc(num_codim * sizeof(prog_t **));
  *fullRankProgInfo = (int **)bmalloc(num_codim * sizeof(int *));
  if (T->MPType == 0)
    *endPts_d = (endpoint_data_d **)bmalloc(num_codim * sizeof(endpoint_data_d *));
  else if (T->MPType == 1)
    *endPts_mp = (endpoint_data_mp **)bmalloc(num_codim * sizeof(endpoint_data_mp *));
  else
    *endPts_amp = (endpoint_data_amp **)bmalloc(num_codim * sizeof(endpoint_data_amp *));

  // setup sliceMover for each codimension
  for (i = 0; i < num_codim; i++)
    basic_setup_slice_moving(&(*sliceMover)[i], W, i, T->MPType, T->AMP_max_prec);

  if (specificCodim == 0 || topDimension)
  { // everything at the top dimension is as expected

    // remove the singular points from the top, if we are only looking for reduced components
    if (reducedOnly)
    { // remove all singular endpoints
      isJunk = (int *)brealloc(isJunk, W->codim[0].num_sing * sizeof(int));
      for (j = 0; j < W->codim[0].num_sing; j++)
        isJunk[j] = 1;

      // remove the singular endpoints
      remove_junk_points(W, 0, T->MPType, isJunk);
    }

    // everything at the top are witness points, so we can setup everything
    W->codim[0].deflations_needed = (int *)bmalloc(W->codim[0].num_set * sizeof(int));
    if (T->MPType == 0)
      deflate_for_junkRemoval(&(*fullRankProgs)[0], &(*fullRankProgInfo)[0], &(*endPts_d)[0], NULL, NULL, &(*sliceMover)[0], W, 0, -1, T, OUT);
    else if (T->MPType == 1)
      deflate_for_junkRemoval(&(*fullRankProgs)[0], &(*fullRankProgInfo)[0], NULL, &(*endPts_mp)[0], NULL, &(*sliceMover)[0], W, 0, -1, T, OUT);
    else
      deflate_for_junkRemoval(&(*fullRankProgs)[0], &(*fullRankProgInfo)[0], NULL, NULL, &(*endPts_amp)[0], &(*sliceMover)[0], W, 0, -1, T, OUT);

    // look to sharpen the witness points, if needed
    if (T->sharpenDigits > 0)
    { // sharpen the top dimensional witness points
      if (T->MPType == 0)
        sharpen_deflated(W, 0, (*fullRankProgs)[0], (*fullRankProgInfo)[0], &(*sliceMover)[0], (*endPts_d)[0], NULL, NULL, T, OUT); 
      else if (T->MPType == 1)
        sharpen_deflated(W, 0, (*fullRankProgs)[0], (*fullRankProgInfo)[0], &(*sliceMover)[0], NULL, (*endPts_mp)[0], NULL, T, OUT); 
      else
        sharpen_deflated(W, 0, (*fullRankProgs)[0], (*fullRankProgInfo)[0], &(*sliceMover)[0], NULL, NULL, (*endPts_amp)[0], T, OUT); 
    }
  }

  // loop through the other codimensions removing the junk points
  for (i = 0; i < num_codim; i++)
    if (i > 0 || (i == 0 && specificCodim > 0 && topDimension == 0))
    { // find the number of singular endpoints for this codimension
      num_sing = W->codim[i].num_sing;

      // setup isJunk - non-singular endpoints are always not junk so we only need to check the ones that are singular
      isJunk = (int *)brealloc(isJunk, num_sing * sizeof(int));
      for (j = 0; j < num_sing; j++)
        isJunk[j] = 0;

      // setup singEndpoints - indices of the singular endpoints
      singEndpoints = (int *)brealloc(singEndpoints, num_sing * sizeof(int));
      k = 0;
      for (j = 0; j < num_sing; j++)
      { // find the index of the jth singular endpoint
        while (W->codim[i].witnessPt_types[k] != SINGULAR)
          k++;

        singEndpoints[j] = k;

        // increment k
        k++;
      }

      // display messages
      printf("\nRemoving junk points from codimension %d: %d endpoints to check.\n", W->codim[i].codim, W->codim[i].num_sing);
      fprintf(OUT, "\n*****************************************************\n");
      fprintf(OUT, "Removing junk points from codimension %d.\n", W->codim[i].codim);
      fprintf(OUT, "*****************************************************\n");

      // loop through the singular endpoints and check if they are junk
      for (j = 0; j < num_sing; j++)
      { // loop through each of the smaller codim to see if this singular endpoint is a member of that codim
        isJunk[j] = 0;

        // print the path number if needed
        if (pathMod > 0 && !(j % pathMod))
          printf("Checking %d of %d\n", j, num_sing);

        // see if this singular points is junk
        if (reducedOnly)
        { // the singular point is junk
          isJunk[j] = 1;
        }
        else
        { // do the membership test
          if (T->junkRemovalTest == 0)
          { // use the standard membership test for junk removal
            isJunk[j] = junkRemoval_mem(W, singEndpoints[j], i, *sliceMover, *fullRankProgs, *fullRankProgInfo, *endPts_d, *endPts_mp, *endPts_amp, T, OUT, midName, my_id, num_processes, headnode);
          }
          else // T->junkRemovalTest == 1
          { // use the local dimension test for junk removal
            isJunk[j] = junkRemoval_ldt(W, singEndpoints[j], i, *sliceMover, *fullRankProgs, *fullRankProgInfo, *endPts_d, *endPts_mp, *endPts_amp, T, OUT, midName, specificCodim, my_id, num_processes, headnode);
          }
        }
      }

      if (specificCodim > 0)
      { // count the number that failed
        k = 0;
        for (j = 0; j < num_sing; j++)
          if (isJunk[j] < 0)
          { // increment and say it is junk
            k++; 
            isJunk[j] = 1;
          }

        if (k > 0)
        {
          printf("\nNOTE: The local dimension test was terminated before it could determine if %d point%s %s isolated.\n", k, k == 1 ? "" : "s", k == 1 ? "was" : "were");
          printf("      Consider either increasing MaxLDTDepth or using MaxCodimension instead of SpecificCodimension.\n");
        }
      }

    // remove the junk points
    remove_junk_points(W, i, T->MPType, isJunk);

    // everything at this codim are witness points, so we can setup everything
    W->codim[i].deflations_needed = (int *)bmalloc(W->codim[i].num_set * sizeof(int));
    if (T->MPType == 0)
      deflate_for_junkRemoval(&(*fullRankProgs)[i], &(*fullRankProgInfo)[i], &(*endPts_d)[i], NULL, NULL, &(*sliceMover)[i], W, i, -1, T, OUT);
    else if (T->MPType == 1)
      deflate_for_junkRemoval(&(*fullRankProgs)[i], &(*fullRankProgInfo)[i], NULL, &(*endPts_mp)[i], NULL, &(*sliceMover)[i], W, i, -1, T, OUT);
    else
      deflate_for_junkRemoval(&(*fullRankProgs)[i], &(*fullRankProgInfo)[i], NULL, NULL, &(*endPts_amp)[i], &(*sliceMover)[i], W, i, -1, T, OUT);

    // look to sharpen the witness points, if needed
    if (T->sharpenDigits > 0)
    { // sharpen the witness points
      if (T->MPType == 0)
        sharpen_deflated(W, i, (*fullRankProgs)[i], (*fullRankProgInfo)[i], &(*sliceMover)[i], (*endPts_d)[i], NULL, NULL, T, OUT); 
      else if (T->MPType == 1)
        sharpen_deflated(W, i, (*fullRankProgs)[i], (*fullRankProgInfo)[i], &(*sliceMover)[i], NULL, (*endPts_mp)[i], NULL, T, OUT); 
      else
        sharpen_deflated(W, i, (*fullRankProgs)[i], (*fullRankProgInfo)[i], &(*sliceMover)[i], NULL, NULL, (*endPts_amp)[i], T, OUT); 
    }
  }

  // print the witness set chart
  witnessSetOutputChart(W, stdout, T->MPType, reducedOnly);

  // free isJunk & singEndpoints
  free(isJunk);
  free(singEndpoints);

  // close OUT
  fclose(OUT);

  return;
}

int junkRemoval_mem(witness_t *W, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform junk removal using membership testing          *
* Loop through the codimensions below this one                  *
\***************************************************************/
{
  int k, isJunk = 0;

  for (k = 0; k < pathNum_codim_index && !isJunk; k++)
  { // perform membership test with the pure-dim witness set for codim_index k and the given point 
    if (T->MPType == 0)
      isJunk = junkRemoval_membershipTest(W, k, pathNum, pathNum_codim_index, &sliceMover[k], fullRankProgs[k], fullRankProgInfo[k], endPts_d[k], NULL, NULL, T, OUT, midName, my_id, num_processes, headnode);
    else if (T->MPType == 1)
      isJunk = junkRemoval_membershipTest(W, k, pathNum, pathNum_codim_index, &sliceMover[k], fullRankProgs[k], fullRankProgInfo[k], NULL, endPts_mp[k], NULL, T, OUT, midName, my_id, num_processes, headnode);
    else
      isJunk = junkRemoval_membershipTest(W, k, pathNum, pathNum_codim_index, &sliceMover[k], fullRankProgs[k], fullRankProgInfo[k], NULL, NULL, endPts_amp[k], T, OUT, midName, my_id, num_processes, headnode);
  }

  return isJunk;
}

void remove_junk_points(witness_t *W, int codim_index, int MPType, int *isJunk)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: removes junk from codim_index                          *
\***************************************************************/
{
  int i, k, singCount, num_junk = 0, num_sing = W->codim[codim_index].num_sing;
  int *temp_types = NULL, *temp_mult = NULL;
  endpoint_data_d *tempPts_d = NULL;
  endpoint_data_mp *tempPts_mp = NULL;
  endpoint_data_amp *tempPts_amp = NULL;

  // count the number of junk points for this codim
  for (i = 0; i < num_sing; i++)
    if (isJunk[i])
      num_junk++;

  if (num_junk > 0)
  { // remove the junk endpoints from this codim
    temp_types = W->codim[codim_index].witnessPt_types;
    W->codim[codim_index].witnessPt_types = (int *)bmalloc((W->codim[codim_index].num_set - num_junk) * sizeof(int));

    temp_mult = W->codim[codim_index].multiplicities;
    W->codim[codim_index].multiplicities = (int *)bmalloc((W->codim[codim_index].num_set - num_junk) * sizeof(int));

    if (MPType == 0)
    { // remove junk from witnessPts_d
      tempPts_d = W->codim[codim_index].witnessPts_d;
      W->codim[codim_index].witnessPts_d = (endpoint_data_d *)bmalloc((W->codim[codim_index].num_set - num_junk) * sizeof(endpoint_data_d));

      // copy back the good points
      k = singCount = 0;
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // see if it is good
        if (temp_types[i] == NON_SINGULAR)
        { // copy over the non-singular endpoints
          init_endpoint_data_d(&W->codim[codim_index].witnessPts_d[k]);
          endpoint_data_cp_d(&W->codim[codim_index].witnessPts_d[k], &tempPts_d[i]);
          W->codim[codim_index].witnessPt_types[k] = temp_types[i];
          W->codim[codim_index].multiplicities[k] = temp_mult[i];
          k++;
        }
        else // endpoint is singular
        { // check to see if it is junk
          if (!isJunk[singCount])
          { // copy over the singular endpoint that is not junk
            init_endpoint_data_d(&W->codim[codim_index].witnessPts_d[k]);
            endpoint_data_cp_d(&W->codim[codim_index].witnessPts_d[k], &tempPts_d[i]);
            W->codim[codim_index].witnessPt_types[k] = temp_types[i];
            W->codim[codim_index].multiplicities[k] = temp_mult[i];
            k++;
          }
          // increment singCount since we have checked out a singular endpoint
          singCount++;
        }
        // clear old points
        clear_endpoint_data_d(&tempPts_d[i]);
      }
      // clear tempPts_d
      free(tempPts_d);
    }
    else if (MPType == 1)
    { // remove junk from witnessPts_mp
      tempPts_mp = W->codim[codim_index].witnessPts_mp;
      W->codim[codim_index].witnessPts_mp = (endpoint_data_mp *)bmalloc((W->codim[codim_index].num_set - num_junk) * sizeof(endpoint_data_mp));

      // copy back the good points and clear the old points
      k = singCount = 0;
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // see if it is good
        if (temp_types[i] == NON_SINGULAR)
        { // copy over the non-singular endpoints
          init_endpoint_data_mp(&W->codim[codim_index].witnessPts_mp[k]);
          endpoint_data_cp_mp(&W->codim[codim_index].witnessPts_mp[k], &tempPts_mp[i]);
          W->codim[codim_index].witnessPt_types[k] = temp_types[i];
          W->codim[codim_index].multiplicities[k] = temp_mult[i];
          k++;
        }
        else // endpoint is singular
        { // check to see if it is junk
          if (!isJunk[singCount])
          { // copy over the singular endpoint that is not junk
            init_endpoint_data_mp(&W->codim[codim_index].witnessPts_mp[k]);
            endpoint_data_cp_mp(&W->codim[codim_index].witnessPts_mp[k], &tempPts_mp[i]);
            W->codim[codim_index].witnessPt_types[k] = temp_types[i];
            W->codim[codim_index].multiplicities[k] = temp_mult[i];
            k++;
          }
          // increment singCount since we have checked out a singular endpoint
          singCount++;
        }
        // clear old points
        clear_endpoint_data_mp(&tempPts_mp[i]);
      }
      // clear tempPts_mp
      free(tempPts_mp);
    }
    else
    { // remove junk from witnessPts_amp
      tempPts_amp = W->codim[codim_index].witnessPts_amp;
      W->codim[codim_index].witnessPts_amp = (endpoint_data_amp *)bmalloc((W->codim[codim_index].num_set - num_junk) * sizeof(endpoint_data_amp));

      // copy back the good points and clear the old points
      k = singCount = 0;
      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // see if it is good
        if (temp_types[i] == NON_SINGULAR)
        { // copy over the non-singular endpoints
          init_endpoint_data_amp(&W->codim[codim_index].witnessPts_amp[k], tempPts_amp[i].curr_prec, tempPts_amp[i].last_approx_prec);
          endpoint_data_cp_amp(&W->codim[codim_index].witnessPts_amp[k], &tempPts_amp[i]);
          W->codim[codim_index].witnessPt_types[k] = temp_types[i];
          W->codim[codim_index].multiplicities[k] = temp_mult[i];
          k++;
        }
        else // endpoint is singular
        { // check to see if it is junk
          if (!isJunk[singCount])
          { // copy over the singular endpoint that is not junk
            init_endpoint_data_amp(&W->codim[codim_index].witnessPts_amp[k], tempPts_amp[i].curr_prec, tempPts_amp[i].last_approx_prec);
            endpoint_data_cp_amp(&W->codim[codim_index].witnessPts_amp[k], &tempPts_amp[i]);
            W->codim[codim_index].witnessPt_types[k] = temp_types[i];
            W->codim[codim_index].multiplicities[k] = temp_mult[i];
            k++;
          }
          // increment singCount since we have checked out a singular endpoint
          singCount++;
        }
        // clear old points
        clear_endpoint_data_amp(&tempPts_amp[i]);
      }
      // clear tempPts_amp
      free(tempPts_amp);
    }

    // clear temp_types & temp_mult
    free(temp_types);
    free(temp_mult);

    // update the counts
    W->codim[codim_index].num_set -= num_junk;
    W->codim[codim_index].num_sing -= num_junk;
  }

  return;
}

void witnessSetOutputChart(witness_t *W, FILE *fp, int MPType, int reducedOnly)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int codim_index, num_codim = W->num_codim;

  if (reducedOnly)
    fprintf(fp, "\n\n*********** Multiplicity 1 Witness Set Summary ************\n\n");
  else
    fprintf(fp, "\n\n*************** Witness Set Summary ****************\n\n");

  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|codim| witness points | nonsingular | singular \n");
  fprintf(fp, "-------------------------------------------------\n");

  for (codim_index = 0; codim_index < num_codim; codim_index++)
  {
    fprintf(fp, "| %-4d|   %-13d|  %-11d|  %-8d\n", W->codim[codim_index].codim, W->codim[codim_index].num_set, W->codim[codim_index].num_nonsing, W->codim[codim_index].num_sing);
  }
  fprintf(fp, "-------------------------------------------------\n\n");

  fprintf(fp, "****************************************************\n\n");

  return;
}

void sharpen_deflated(witness_t *W, int codim_index, prog_t **fullRankProgs, int *fullRankProgInfo, membership_slice_moving_t *sliceMover, endpoint_data_d *endPts_d, endpoint_data_mp *endPts_mp, endpoint_data_amp *endPts_amp, tracker_config_t *T, FILE *OUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sharpens the deflated components                       *
\***************************************************************/
{
  int i, retVal, setupGamma = 1, num_paths = W->codim[codim_index].num_set;
  int prec_in, prec_out;
  point_data_d in_d, out_d;
  point_data_mp in_mp, out_mp;

  // initialize
  init_point_data_d(&in_d, 0); init_point_data_d(&out_d, 0);
  init_point_data_mp(&in_mp, 0); init_point_data_mp(&out_mp, 0);

  // make sure that sliceMover is setup properly - only need the start slice to be == 0
  setup_random_slice_moving(sliceMover, T->MPType, T->AMP_max_prec);

  // loop through the start points
  for (i = 0; i < num_paths; i++)
  { // verify that deflation worked
    if (fullRankProgInfo[i] != -1)
    { // finish the setup for sliceMover for this endpoint
      final_setup_slice_moving(sliceMover, fullRankProgs[i], T->MPType, T->AMP_max_prec, setupGamma);

      if (T->MPType == 0)
      { // use double precision sharpening

        // setup in_d
        prec_in = 52;
        point_cp_d(in_d.point, endPts_d[i].endPt);
        set_one_d(in_d.time);

        //  do the actual sharpening
        retVal = sharpen_zero_main(T->final_tolerance, NULL, 52, &out_d, NULL, &prec_out, &in_d, NULL, prec_in, T, OUT, sliceMover, sliceMover, slice_mover_eval_d, slice_mover_eval_mp, slice_mover_change_prec);

        // determine if successful
        if (!retVal)
        { // copy back to endPts
          point_cp_d(endPts_d[i].endPt, out_d.point);        

          // project onto the proper number of coordinates and copy to W
          out_d.point->size = W->codim[codim_index].witnessPts_d[i].endPt->size;
          point_cp_d(W->codim[codim_index].witnessPts_d[i].endPt, out_d.point);
        }
      }
      else if (T->MPType == 1)
      { // use fixed multi precision sharpening

        // setup in_mp
        prec_in = T->Precision;
        point_cp_mp(in_mp.point, endPts_mp[i].endPt);
        set_one_mp(in_mp.time);

        //  do the actual sharpening
        retVal = sharpen_zero_main(T->final_tolerance, NULL, 52, NULL, &out_mp, &prec_out, NULL, &in_mp, prec_in, T, OUT, sliceMover, sliceMover, slice_mover_eval_d, slice_mover_eval_mp, slice_mover_change_prec);

        // determine if successful
        if (!retVal)
        { // copy back to endPts
          point_cp_mp(endPts_mp[i].endPt, out_mp.point);        

          // project onto the proper number of coordinates and copy to W
          out_mp.point->size = W->codim[codim_index].witnessPts_mp[i].endPt->size;
          point_cp_mp(W->codim[codim_index].witnessPts_mp[i].endPt, out_mp.point);
        }
      }
      else
      { // use adaptive precision sharpening

        // setup in
        prec_in = endPts_amp[i].curr_prec;
        if (prec_in < 64)
        { // setup in_d
          point_cp_d(in_d.point, endPts_amp[i].endPt_d);
          set_one_d(in_d.time);
        }
        else
        { // setup in_mp
          setprec_point_data_mp(&in_mp, prec_in);
          point_cp_mp(in_mp.point, endPts_amp[i].endPt_mp);
          set_one_mp(in_mp.time);
        }

        //  do the actual sharpening
        retVal = sharpen_zero_main(T->final_tolerance, NULL, 52, &out_d, &out_mp, &prec_out, &in_d, &in_mp, prec_in, T, OUT, sliceMover, sliceMover, slice_mover_eval_d, slice_mover_eval_mp, slice_mover_change_prec);

        // determine if successful
        if (!retVal)
        { // copy back to endPts
          if (prec_out < 64)
          { // copy using out_d 
            endPts_amp[i].curr_prec = prec_out;
            point_cp_d(endPts_amp[i].endPt_d, out_d.point);
 
            // project onto the proper number of coordinates and copy to W
            if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              out_d.point->size = W->codim[codim_index].witnessPts_amp[i].endPt_d->size;
            else
              out_d.point->size = W->codim[codim_index].witnessPts_amp[i].endPt_mp->size;
            W->codim[codim_index].witnessPts_amp[i].curr_prec = prec_out;
            point_cp_d(W->codim[codim_index].witnessPts_amp[i].endPt_d, out_d.point);
          }
          else
          { // copy using out_mp
            endPts_amp[i].curr_prec = prec_out;
            setprec_point_mp(endPts_amp[i].endPt_mp, prec_out);
            point_cp_mp(endPts_amp[i].endPt_mp, out_mp.point);
 
            // project onto the proper number of coordinates and copy to W
            if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              out_mp.point->size = W->codim[codim_index].witnessPts_amp[i].endPt_d->size;
            else
              out_mp.point->size = W->codim[codim_index].witnessPts_amp[i].endPt_mp->size;
            W->codim[codim_index].witnessPts_amp[i].curr_prec = prec_out;
            setprec_point_mp(W->codim[codim_index].witnessPts_amp[i].endPt_mp, prec_out);
            point_cp_mp(W->codim[codim_index].witnessPts_amp[i].endPt_mp, out_mp.point);
          }
        }
      }
    } 
  }

  // clear 
  clear_point_data_d(&in_d); clear_point_data_d(&out_d);
  clear_point_data_mp(&in_mp); clear_point_data_mp(&out_mp);

  return;
}

//////// main worker parallel control function //////////////

void worker_numericalIrredDecomp(int my_id, int num_processes, int headnode, int dataType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: does numerical irreducible decomposition - 'worker'    *
\***************************************************************/
{
#ifdef _HAVE_MPI
  worker_info recvType;

  // recv the worker info
  bcast_worker_info(&recvType, my_id, headnode);

  // determine how the witness supersets are being generated
  if (recvType.dataType == DIMBYDIM)
  { // do dimension-by-dimension
    worker_dimbydim(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == CASCADE)
  { // do cascade 
    worker_cascade(my_id, num_processes, headnode, recvType.dataType);
  }
  else if (recvType.dataType == REGEN_POS_DIM)
  { // do regen pos dim
    worker_regen_pos_dim(my_id, num_processes, headnode, recvType.dataType);
  }
  else
  {
    printf("ERROR: The data type is incorrect!!\n");
    bexit(ERROR_CONFIGURATION);
  }

  // determine if decomposition is needed
  bcast_worker_info(&recvType, my_id, headnode);

  if (recvType.dataType != STOPCODE)
  { // now do the junk removal & break up into irreducible components
    worker_witness_superset_decomposition(my_id, num_processes, headnode);
  }
#endif

  return;
}

////////////// setup witness set to be totally extrinsic //////////////

void setupWitnessTotallyExtrinisic(witness_t *W, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: change the witness points to the original homogeneous  *
* coordinates as well as updating H, B, & p                     *
\***************************************************************/
{
  int i, j, k, l, rows, cols, C_prec = W->curr_precision;
  mat_d C_perp_d, N_d;
  mat_mp C_perp_mp, N_mp;
  mpq_t ***C_perp_rat = NULL, ***N_rat = NULL;

  // make sure that there is something to do!
  if (W->orig_variables != W->new_variables)
  { // convert witess points
    if (MPType == 0)
    { // convert _d
      for (i = 0; i < W->num_codim; i++)
        for (j = 0; j < W->codim[i].num_set; j++)
        { // convert witness points
          mul_mat_vec_d(W->codim[i].witnessPts_d[j].endPt, W->C_d, W->codim[i].witnessPts_d[j].endPt);
          mul_mat_vec_d(W->codim[i].witnessPts_d[j].last_approx, W->C_d, W->codim[i].witnessPts_d[j].last_approx);
        }
    }
    else if (MPType == 1)
    { // convert _mp
      for (i = 0; i < W->num_codim; i++)
        for (j = 0; j < W->codim[i].num_set; j++)
        { // convert witness points
          mul_mat_vec_mp(W->codim[i].witnessPts_mp[j].endPt, W->C_mp, W->codim[i].witnessPts_mp[j].endPt);
          mul_mat_vec_mp(W->codim[i].witnessPts_mp[j].last_approx, W->C_mp, W->codim[i].witnessPts_mp[j].last_approx);
        }
    }
    else
    { // convert _amp
      for (i = 0; i < W->num_codim; i++)
        for (j = 0; j < W->codim[i].num_set; j++)
        { // convert endPt
          if (W->codim[i].witnessPts_amp[j].curr_prec < 64)
          {
            mul_mat_vec_d(W->codim[i].witnessPts_amp[j].endPt_d, W->C_d, W->codim[i].witnessPts_amp[j].endPt_d);
          }
          else
          { // make sure C_mp is set to high enough precision
            if (W->codim[i].witnessPts_amp[j].curr_prec > C_prec)
            { // set the precision
              C_prec = W->codim[i].witnessPts_amp[j].curr_prec;
              for (k = 0; k < W->C_mp->rows; k++)
                for (l = 0; l < W->C_mp->cols; l++)
                {
                  setprec_mp(&W->C_mp->entry[k][l], C_prec);
                  mpf_set_q(W->C_mp->entry[k][l].r, W->C_rat[k][l][0]);
                  mpf_set_q(W->C_mp->entry[k][l].i, W->C_rat[k][l][1]);
                }
            }
            mul_mat_vec_mp(W->codim[i].witnessPts_amp[j].endPt_mp, W->C_mp, W->codim[i].witnessPts_amp[j].endPt_mp);
          }

          // convert last_approx
          if (W->codim[i].witnessPts_amp[j].last_approx_prec < 64)
          {
            mul_mat_vec_d(W->codim[i].witnessPts_amp[j].last_approx_d, W->C_d, W->codim[i].witnessPts_amp[j].last_approx_d);
          }
          else
          { // make sure C_mp is set to high enough precision
            if (W->codim[i].witnessPts_amp[j].curr_prec > C_prec)
            { // set the precision
              C_prec = W->codim[i].witnessPts_amp[j].curr_prec; 
              for (k = 0; k < W->C_mp->rows; k++)
                for (l = 0; l < W->C_mp->cols; l++)
                {
                  setprec_mp(&W->C_mp->entry[k][l], C_prec);
                  mpf_set_q(W->C_mp->entry[k][l].r, W->C_rat[k][l][0]);
                  mpf_set_q(W->C_mp->entry[k][l].i, W->C_rat[k][l][1]);
                }
            }
            mul_mat_vec_mp(W->codim[i].witnessPts_amp[j].last_approx_mp, W->C_mp, W->codim[i].witnessPts_amp[j].last_approx_mp);
          }
        }
    }

    // setup C_perp
    if (MPType == 0)
    { // setup C_perp_d
      init_mat_d(C_perp_d, W->C_d->cols, W->C_d->rows);
      mat_perp_d(C_perp_d, W->C_d);
    }
    else if (MPType == 1)
    { // setup C_perp_mp
      init_mat_mp2(C_perp_mp, W->C_mp->cols, W->C_mp->rows, W->curr_precision);
      mat_perp_mp(C_perp_mp, W->C_mp, W->curr_precision);
    }
    else
    { // setup C_perp_d, _mp, _rat
      init_mat_d(C_perp_d, W->C_d->cols, W->C_d->rows);
      init_mat_mp2(C_perp_mp, W->C_mp->cols, W->C_mp->rows, W->curr_precision);
      init_mat_rat(C_perp_rat, W->C_d->cols, W->C_d->rows);
      mat_perp_rat(C_perp_d, C_perp_mp, C_perp_rat, W->C_rat, W->C_d->rows, W->C_d->cols, W->curr_precision, max_prec, 0);
    }

    // setup N
    if (MPType == 0)
    { // setup N_d
      init_mat_d(N_d, 0, 0);
      intrinsicToExtrinsicMat_d(N_d, W->C_d);
    }
    else if (MPType == 1)
    { // setup N_mp
      init_mat_mp2(N_mp, 0, 0, W->curr_precision);
      intrinsicToExtrinsicMat_mp(N_mp, W->C_mp);
    }
    else
    { // setup N_d, N_mp & N_rat
      intrinsicToExtrinsicMat_rat(&N_rat, &rows, &cols, W->C_rat, W->C_d->rows, W->C_d->cols, W->curr_precision, max_prec);

      // setup N_d & N_mp
      init_mat_d(N_d, rows, cols);
      init_mat_mp2(N_mp, rows, cols, W->curr_precision);
      N_d->rows = N_mp->rows = rows;
      N_d->cols = N_mp->cols = cols;
      for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
        {
          mpf_set_q(N_mp->entry[i][j].r, N_rat[i][j][0]);
          mpf_set_q(N_mp->entry[i][j].i, N_rat[i][j][1]);
          N_d->entry[i][j].r = mpq_get_d(N_rat[i][j][0]);
          N_d->entry[i][j].i = mpq_get_d(N_rat[i][j][1]);
        }
    }

    // update H = H * C_perp, p = p * C_perp & B = [B * C_perp, N]
    if (MPType == 0)
    { // update _d values
      for (i = 0; i < W->num_codim; i++)
      { // update H_d, p_d & B_d
        vec_mat_mul_d(W->codim[i].H_d, W->codim[i].H_d, C_perp_d);
        vec_mat_mul_d(W->codim[i].p_d, W->codim[i].p_d, C_perp_d);
        mat_mul_d(W->codim[i].B_d, W->codim[i].B_d, C_perp_d);

        // put N on the bottom of B
        rows = W->codim[i].B_d->rows;
        cols = W->codim[i].B_d->cols;
        increase_size_mat_d(W->codim[i].B_d, rows + N_d->rows, cols);
        W->codim[i].B_d->rows = rows + N_d->rows;
        for (k = 0; k < N_d->rows; k++)
          for (j = 0; j < cols; j++)
          {
            set_d(&W->codim[i].B_d->entry[k + rows][j], &N_d->entry[k][j]);
          }
      }
    }
    else if (MPType == 1)
    { // update _mp values
      for (i = 0; i < W->num_codim; i++)
      { // update H_mp, p_mp & B_mp
        vec_mat_mul_mp(W->codim[i].H_mp, W->codim[i].H_mp, C_perp_mp);
        vec_mat_mul_mp(W->codim[i].p_mp, W->codim[i].p_mp, C_perp_mp);
        mat_mul_mp(W->codim[i].B_mp, W->codim[i].B_mp, C_perp_mp);

        // put N on the bottom of B
        rows = W->codim[i].B_mp->rows;
        cols = W->codim[i].B_mp->cols;
        increase_size_mat_mp(W->codim[i].B_mp, rows + N_mp->rows, cols);
        W->codim[i].B_mp->rows = rows + N_mp->rows;
        for (k = 0; k < N_mp->rows; k++)
          for (j = 0; j < cols; j++)
          {
            set_mp(&W->codim[i].B_mp->entry[k + rows][j], &N_mp->entry[k][j]);
          }
      }
    }
    else
    { // update amp values
      mpq_t **tempVec = NULL, ***tempMat = NULL;

      for (i = 0; i < W->num_codim; i++)
      { // update H
        rows = W->codim[i].H_d->size;
        tempVec = W->codim[i].H_rat;

        // setup H_rat to be the correct size
        init_vec_rat(W->codim[i].H_rat, W->orig_variables);

        // find H_rat = tempVec * C_perp_rat
        vec_mat_mul_rat(W->codim[i].H_rat, tempVec, C_perp_rat, rows, C_perp_d->rows, C_perp_d->cols, 0);

        // setup H_d & H_mp
        increase_size_vec_d(W->codim[i].H_d, W->orig_variables);
        increase_size_vec_mp(W->codim[i].H_mp, W->orig_variables);
        rows = W->codim[i].H_d->size = W->codim[i].H_mp->size = W->orig_variables;
        for (j = 0; j < rows; j++)
        {
          mpf_set_q(W->codim[i].H_mp->coord[j].r, W->codim[i].H_rat[j][0]);
          mpf_set_q(W->codim[i].H_mp->coord[j].i, W->codim[i].H_rat[j][1]);
          W->codim[i].H_d->coord[j].r = mpq_get_d(W->codim[i].H_rat[j][0]);
          W->codim[i].H_d->coord[j].i = mpq_get_d(W->codim[i].H_rat[j][1]);
        }

        // clear tempVec
        clear_vec_rat(tempVec, W->new_variables);

        // update p
        rows = W->codim[i].p_d->size;
        tempVec = W->codim[i].p_rat;

        // setup p_rat to be the correct size
        init_vec_rat(W->codim[i].p_rat, W->orig_variables);

        // find p_rat = tempVec * C_perp_rat
        vec_mat_mul_rat(W->codim[i].p_rat, tempVec, C_perp_rat, rows, C_perp_d->rows, C_perp_d->cols, 0);

        // setup p_d & p_mp
        increase_size_vec_d(W->codim[i].p_d, W->orig_variables);
        increase_size_vec_mp(W->codim[i].p_mp, W->orig_variables);
        rows = W->codim[i].p_d->size = W->codim[i].p_mp->size = W->orig_variables;
        for (j = 0; j < rows; j++)
        {
          mpf_set_q(W->codim[i].p_mp->coord[j].r, W->codim[i].p_rat[j][0]);
          mpf_set_q(W->codim[i].p_mp->coord[j].i, W->codim[i].p_rat[j][1]);
          W->codim[i].p_d->coord[j].r = mpq_get_d(W->codim[i].p_rat[j][0]);
          W->codim[i].p_d->coord[j].i = mpq_get_d(W->codim[i].p_rat[j][1]);
        }

        // clear tempVec
        clear_vec_rat(tempVec, W->new_variables);

        // update B
        rows = W->codim[i].B_d->rows;
        cols = W->codim[i].B_d->cols;
        tempMat = W->codim[i].B_rat;

        // setup B_rat to be the correct size
        init_mat_rat(W->codim[i].B_rat, rows + N_d->rows , W->orig_variables);

        // find B_rat = [tempMat * C_perp, N]
        mat_mul_rat(W->codim[i].B_rat, tempMat, C_perp_rat, rows, cols, C_perp_d->rows, C_perp_d->cols, max_prec);
        for (k = 0; k < N_d->rows; k++)
          for (j = 0; j < N_d->cols; j++)
          {
            mpq_set(W->codim[i].B_rat[k + rows][j][0], N_rat[k][j][0]);
            mpq_set(W->codim[i].B_rat[k + rows][j][1], N_rat[k][j][1]);
          }

        // setup B_d & B_mp
        increase_size_mat_d(W->codim[i].B_d, rows + N_d->rows, W->orig_variables);
        increase_size_mat_mp(W->codim[i].B_mp, rows + N_d->rows, W->orig_variables);
        W->codim[i].B_d->rows = W->codim[i].B_mp->rows = rows + N_d->rows;
        W->codim[i].B_d->cols = W->codim[i].B_mp->cols = W->orig_variables;
        for (j = 0; j < W->codim[i].B_d->rows; j++)
          for (k = 0; k < W->codim[i].B_d->cols; k++)
          {
            mpf_set_q(W->codim[i].B_mp->entry[j][k].r, W->codim[i].B_rat[j][k][0]);
            mpf_set_q(W->codim[i].B_mp->entry[j][k].i, W->codim[i].B_rat[j][k][1]);
            W->codim[i].B_d->entry[j][k].r = mpq_get_d(W->codim[i].B_rat[j][k][0]);
            W->codim[i].B_d->entry[j][k].i = mpq_get_d(W->codim[i].B_rat[j][k][1]);
          }

        // clear tempMat
        clear_mat_rat(tempMat, rows, cols);
      }
    }

    // set the number of variables
    W->new_variables = W->orig_variables;

    // clear the matrics
    clear_mat(C_perp_d, C_perp_mp, C_perp_rat, MPType);
    clear_mat(N_d, N_mp, N_rat, MPType);
  }

  // make sure all codimensions have the same patch
  rows = 1;
  cols = W->orig_variables;
  if (MPType == 0)
  { // setup the same patch in double precision
    mat_d new_patch;
    init_mat_d(new_patch, rows, cols); 
    make_matrix_random_d(new_patch, rows, cols);

    for (i = 0; i < W->num_codim; i++)
    { // copy the patch
      for (j = 0; j < cols; j++)
      {
        set_d(&W->codim[i].p_d->coord[j], &new_patch->entry[0][j]);
      }
      // move the points to this patch
      for (j = 0; j < W->codim[i].num_set; j++)
      { // move both end points and last approximations to the patch
        move_to_patch_mat_d(W->codim[i].witnessPts_d[j].endPt, W->codim[i].witnessPts_d[j].endPt, new_patch, &W->PPD);
        move_to_patch_mat_d(W->codim[i].witnessPts_d[j].last_approx, W->codim[i].witnessPts_d[j].last_approx, new_patch, &W->PPD);
      }
    }

    clear_mat_d(new_patch);
  }
  else if (MPType == 1)
  { // setup the same patch in fixed multi precision
    mat_mp new_patch;
    init_mat_mp(new_patch, rows, cols);
    make_matrix_random_mp(new_patch, rows, cols, W->curr_precision);

    for (i = 0; i < W->num_codim; i++)
    { // copy the patch
      for (j = 0; j < cols; j++)
      {
        set_mp(&W->codim[i].p_mp->coord[j], &new_patch->entry[0][j]);
      }
      // move the points to this patch
      for (j = 0; j < W->codim[i].num_set; j++)
      { // move both end points and last approximations to the patch
        move_to_patch_mat_mp(W->codim[i].witnessPts_mp[j].endPt, W->codim[i].witnessPts_mp[j].endPt, new_patch, &W->PPD);
        move_to_patch_mat_mp(W->codim[i].witnessPts_mp[j].last_approx, W->codim[i].witnessPts_mp[j].last_approx, new_patch, &W->PPD);
      }
    }

    clear_mat_mp(new_patch);
  }
  else
  { // setup the same patch in _d, _mp, _rat
    mat_d new_patch_d;
    mat_mp new_patch_mp;
    mpq_t ***new_patch_rat = NULL;

    C_prec = W->curr_precision;
    init_mat_d(new_patch_d, rows, cols);
    init_mat_mp2(new_patch_mp, rows, cols, C_prec);
    init_mat_rat(new_patch_rat, rows, cols);
    make_matrix_random_rat(new_patch_d, new_patch_mp, new_patch_rat, rows, cols, W->curr_precision, max_prec, 0, 0);

    for (i = 0; i < W->num_codim; i++)
    { // copy the patch
      for (j = 0; j < cols; j++)
      {
        mpq_set(W->codim[i].p_rat[j][0], new_patch_rat[0][j][0]);
        mpq_set(W->codim[i].p_rat[j][1], new_patch_rat[0][j][1]);
        set_mp(&W->codim[i].p_mp->coord[j], &new_patch_mp->entry[0][j]);
        set_d(&W->codim[i].p_d->coord[j], &new_patch_d->entry[0][j]);
      }

      // move the points to this patch
      for (j = 0; j < W->codim[i].num_set; j++)
      { // move end point to the patch
        if (W->codim[i].witnessPts_amp[j].curr_prec < 64)
        { // move endPt_d
          move_to_patch_mat_d(W->codim[i].witnessPts_amp[j].endPt_d, W->codim[i].witnessPts_amp[j].endPt_d, new_patch_d, &W->PPD);
        }
        else
        { // make sure new_patch_mp has high enough precision
          if (W->codim[i].witnessPts_amp[j].curr_prec > C_prec)
          { // increase precision
            C_prec = W->codim[i].witnessPts_amp[j].curr_prec;
            for (k = 0; k < cols; k++)
            {
              setprec_mp(&new_patch_mp->entry[0][k], C_prec);
              mpf_set_q(new_patch_mp->entry[0][k].r, new_patch_rat[0][k][0]);
              mpf_set_q(new_patch_mp->entry[0][k].i, new_patch_rat[0][k][1]);
            }
          }
          move_to_patch_mat_mp(W->codim[i].witnessPts_amp[j].endPt_mp, W->codim[i].witnessPts_amp[j].endPt_mp, new_patch_mp, &W->PPD);
        }

        // move last approx to the patch
        if (W->codim[i].witnessPts_amp[j].last_approx_prec < 64)
        { // move last_approx_d
          move_to_patch_mat_d(W->codim[i].witnessPts_amp[j].last_approx_d, W->codim[i].witnessPts_amp[j].last_approx_d, new_patch_d, &W->PPD);
        }
        else
        { // make sure new_patch_mp has high enough precision
          if (W->codim[i].witnessPts_amp[j].last_approx_prec > C_prec)
          { // increase precision
            C_prec = W->codim[i].witnessPts_amp[j].curr_prec;
            for (k = 0; k < cols; k++)
            {
              setprec_mp(&new_patch_mp->entry[0][k], C_prec);
              mpf_set_q(new_patch_mp->entry[0][k].r, new_patch_rat[0][k][0]);
              mpf_set_q(new_patch_mp->entry[0][k].i, new_patch_rat[0][k][1]);
            }
          }
          move_to_patch_mat_mp(W->codim[i].witnessPts_amp[j].last_approx_mp, W->codim[i].witnessPts_amp[j].last_approx_mp, new_patch_mp, &W->PPD);
        }
      }
    }
    clear_mat(new_patch_d, new_patch_mp, new_patch_rat, MPType);
  }

  return;
}

////////////////// main numerical irreducible decomposition output /////////////////////////////

void numIrredDecompChart(witness_t *W, FILE *fp, int MPType, int reducedOnly)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int i, j, count, codim_index, num_codim = W->num_codim, num_components, classified, unclassified;
  int *degrees = NULL;

  if (reducedOnly)
    fprintf(fp, "\n\n********* Multiplicity 1 Witness Set Decomposition *********\n\n");
  else
    fprintf(fp, "\n\n************* Witness Set Decomposition *************\n\n");

  fprintf(fp, "| dimension | components | classified | unclassified\n");
  fprintf(fp, "-----------------------------------------------------\n");

  for (codim_index = 0; codim_index < num_codim; codim_index++)
    if (W->codim[codim_index].num_set > 0)
    { // determine classified & unclassified
      num_components = W->codim[codim_index].num_components;
      classified = unclassified = 0;
      for (i = 0; i < W->codim[codim_index].num_set; i++)
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < num_components)
          classified++;
        else
          unclassified++;

      fprintf(fp, "|   %-8d|   %-9d|   %-9d|  %d\n", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, num_components, classified, unclassified);
    }

  fprintf(fp, "-----------------------------------------------------\n\n");

  fprintf(fp, "************** Decomposition by Degree **************\n\n");

  for (codim_index = 0; codim_index < num_codim; codim_index++)
  { // determine the degree of each component
    num_components = W->codim[codim_index].num_components;

    if (num_components > 0)
    {
      degrees = (int *)brealloc(degrees, num_components * sizeof(int));
      for (i = 0; i < num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      // find the minimum and maximum degree
      classified = 0;
      unclassified = W->codim[codim_index].num_set + 1;
      for (i = 0; i < num_components; i++)
      { // find maximum degree
        if (classified < degrees[i])
          classified = degrees[i]; 

        // find minimum degree
        if (unclassified > degrees[i])
          unclassified = degrees[i];
      }

      // display header
      fprintf(fp, "Dimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, num_components);
      if (num_components == 1)
        fprintf(fp, "\n");
      else
        fprintf(fp, "s\n");
      fprintf(fp, "-----------------------------------------------------\n");

      // display the summary
      for (i = unclassified; i <= classified; i++)
      { // count the number that have degree == i
        count = 0;
        for (j = 0; j < num_components; j++)
          if (degrees[j] == i)
            count++;

        if (count > 0)
        { // display the number
          fprintf(fp, "   degree %d: %d component", i, count);
          if (count == 1)
            fprintf(fp, "\n");
          else
            fprintf(fp, "s\n");
        }
      }

      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "*****************************************************\n\n");

  free(degrees);

  return;
}

void numIrredDecompOutput(witness_t *W, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: create the output files                                *
\***************************************************************/
{
  // generate main_data
  numIrredDecompMainData("main_data", W, T, trackType, genType, randomSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // generate witness_data
  numIrredDecompWitnessData("witness_data", W, T->MPType);

  return;
}

void numIrredDecompWitnessData(char *witnessName, witness_t *W, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: generates the file containing the witness data         *
\***************************************************************/
{
  int codim_index, count, num_codim = W->num_codim;
  FILE *OUT = fopen(witnessName, "w");

  // print the number of variables
  fprintf(OUT, "%d\n", W->orig_variables);

  // find the number of codimensions that have witness points and print this to OUT
  count = 0;
  for (codim_index = 0; codim_index < num_codim; codim_index++)
    if (W->codim[codim_index].num_set > 0)
      count++;
  fprintf(OUT, "%d\n", count);

  // loop through to print the witnesss points for each of the codimensions that have witness points
  for (codim_index = 0; codim_index < num_codim; codim_index++)
    if (W->codim[codim_index].num_set > 0)
    {
      printCodimWitnessSet(OUT, &W->codim[codim_index], MPType);
    }

  // print the end of the witness set and the MPType for the codim structures
  fprintf(OUT, "%d\n\n%d\n", -1, MPType);

  // loop through to print the structures for each of the codimensions that have witness points
  for (codim_index = 0; codim_index < num_codim; codim_index++)
    if (W->codim[codim_index].num_set > 0)
    {
      printCodimWitnessStructures(OUT, &W->codim[codim_index], MPType);
    }

  fclose(OUT);

  return;
}

void numIrredDecompMainData(char *mainName, witness_t *W, tracker_config_t *T, int trackType, int genType, unsigned int randomSeed, int pathMod, int userHom, int useRegen, int regenStartLevel, int maxCodim, int specificCodim, double intrinsicCutoffMultiplier, int reducedOnly, int constructWitnessSet, int supersetOnly, int paramHom)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the data to mainName (usually main_data)        *
\***************************************************************/
{
  int i, j, k, size, print_count, codim_index, codim, dim, num_set, num_codim = W->num_codim;
  char ch, **name_table = (char **)bmalloc(W->orig_variables * sizeof(char *));
  FILE *NAMES = NULL, *OUT = fopen(mainName, "w");

  // read in the names of the variables
  NAMES = fopen("names.out", "r");
  if (NAMES == NULL)
  {
    printf("ERROR: 'names.out' does not exist!\n");
    bexit(ERROR_FILE_NOT_EXIST);
  }
  for (i = 0; i < W->orig_variables; i++)
  { // initial allocation
    size = 1;
    name_table[i] = (char *)bmalloc(size * sizeof(char));
    // read in name
    while ((name_table[i][size - 1] = fgetc(NAMES)) != '\n')
    {
      size++;
      name_table[i] = (char *)brealloc(name_table[i], size * sizeof(char));
    }
    name_table[i][size - 1] = '\0';
  }
  fclose(NAMES);

  fprintf(OUT, "Number of variables: %d\n", W->orig_variables - W->PPD.num_var_gp);
  fprintf(OUT, "Variables: ");
  for (i = W->PPD.num_var_gp; i < W->orig_variables; i++)
    fprintf(OUT, " %s", name_table[i]);
  fprintf(OUT, "\nRank: %d\n", W->system_rank);

  if (T->MPType == 0)
  { // print answers in double precision
    point_d dehom;
    init_point_d(dehom, 0);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    {
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular grouped by component
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_d[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_d(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }

        // print extraneous non-singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_nonsing)
        {
          fprintf(OUT, "\nUNCLASSIFIED NONSINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_d[i].cond_num);

              // find the dehomogenized point
              witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_d(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }
        }
      }

      // print the singular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_sing > 0)
      { // print singular grouped by component
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_d[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_d(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]);
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }

        // print extraneous singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_sing)
        {
          fprintf(OUT, "\nUNCLASSIFIED SINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_d[i].cond_num);

              // find the dehomogenized point
              witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_d(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]);
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }
        }
      }
    }
    clear_point_d(dehom);
  }
  else if (T->MPType == 1)
  { // print answers in multi precision
    point_mp dehom;
    init_point_mp(dehom, 0);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    {
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular group by component
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_mp[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_mp(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }
      
        // print extraneous non-singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_nonsing)
        { 
          fprintf(OUT, "\nUNCLASSIFIED NONSINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_mp[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_mp(OUT, 0, &dehom->coord[k]);
               fprintf(OUT, "\n");
              } 
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }
        }
      }

      // print the singular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_sing > 0)
      { // print singular grouped by component
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_mp[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_mp(OUT, 0, &dehom->coord[k]);
                fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]);
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }

        // print extraneous non-singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_sing)
        {
          fprintf(OUT, "\nUNCLASSIFIED SINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_mp[i].cond_num);
              // find the dehomogenized point
              witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
              // print dehomogenized point
              for (k = 0; k < dehom->size; k++)
              {
                print_mp(OUT, 0, &dehom->coord[k]);
               fprintf(OUT, "\n");
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }
        }
      }
    }

    // clear dehom
    clear_point_mp(dehom);
  }
  else
  { // print answers using AMP
    point_d dehom_d;
    point_mp dehom_mp;
    init_point_d(dehom_d, 0);
    init_point_mp2(dehom_mp, 0, 64);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    {
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular group by component
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_amp[i].cond_num);

              if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              { // find the dehomogenized point using double precision
                witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
                // print dehomogenized point
                for (k = 0; k < dehom_d->size; k++)
                {
                  print_d(OUT, 0, &dehom_d->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              else
              { // find the dehomogenized point using multi precision

                // set the precision so that the calculations are correct
                initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
                setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // find the dehomogenized point
                witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // print dehomogenized point
                for (k = 0; k < dehom_mp->size; k++)
                {
                  print_mp(OUT, 0, &dehom_mp->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }

        // print extraneous non-singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_nonsing)
        {
          fprintf(OUT, "\nUNCLASSIFIED NONSINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_amp[i].cond_num);

              if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              { // find the dehomogenized point using double precision
                witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
                // print dehomogenized point
                for (k = 0; k < dehom_d->size; k++)
                {
                  print_d(OUT, 0, &dehom_d->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              else
              { // find the dehomogenized point using multi precision

                // set the precision so that the calculations are correct
                initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
                setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // find the dehomogenized point
                witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // print dehomogenized point
                for (k = 0; k < dehom_mp->size; k++)
                {
                  print_mp(OUT, 0, &dehom_mp->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]); // better be == 0!!!!
            }
        }
      }

      // print the singular solutions, if they exist
      print_count = 0;
      if (W->codim[codim_index].num_sing > 0)
      { // print singular grouped by component
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        for (j = 0; j < W->codim[codim_index].num_components; j++)
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] == j && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              print_count++;

              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Component number: %d\n", W->codim[codim_index].component_nums[i]);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_amp[i].cond_num);

              if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              { // find the dehomogenized point using double precision
                witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
                // print dehomogenized point
                for (k = 0; k < dehom_d->size; k++)
                {
                  print_d(OUT, 0, &dehom_d->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              else
              { // find the dehomogenized point using multi precision

                // set the precision so that the calculations are correct
                initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
                setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // find the dehomogenized point
                witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // print dehomogenized point
                for (k = 0; k < dehom_mp->size; k++)
                {
                  print_mp(OUT, 0, &dehom_mp->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]);
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }
          
        // print extraneous non-singular solutions, if they exist
        if (print_count != W->codim[codim_index].num_sing)
        {
          fprintf(OUT, "\nUNCLASSIFIED SINGULAR SOLUTIONS\n");
          for (i = 0; i < num_set; i++)
            if (W->codim[codim_index].component_nums[i] < 0 && W->codim[codim_index].witnessPt_types[i] == SINGULAR)
            {
              fprintf(OUT, "---------------\n");
              fprintf(OUT, "Path number: %d\n", i);
              fprintf(OUT, "Estimated condition number: %e\n", W->codim[codim_index].witnessPts_amp[i].cond_num);

              if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
              { // find the dehomogenized point using double precision
                witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
                // print dehomogenized point
                for (k = 0; k < dehom_d->size; k++)
                {
                  print_d(OUT, 0, &dehom_d->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              else
              { // find the dehomogenized point using multi precision

                // set the precision so that the calculations are correct
                initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
                setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // find the dehomogenized point
                witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

                // print dehomogenized point
                for (k = 0; k < dehom_mp->size; k++)
                {
                  print_mp(OUT, 0, &dehom_mp->coord[k]);
                  fprintf(OUT, "\n");
                }
              }
              fprintf(OUT, "Multiplicity: %d\n", W->codim[codim_index].multiplicities[i]);
              fprintf(OUT, "Deflations needed: %d\n", W->codim[codim_index].deflations_needed[i]);
            }
        }
      }
    }

    // clear dehom
    clear_point_d(dehom_d);
    clear_point_mp(dehom_mp);
  }

  // print the input onto the bottom
  NAMES = fopen("func_input", "r");

  fprintf(OUT, "\n\n*************** input file needed to reproduce this run ***************\n\n");

  // print the configurations
  printConfigValues(OUT, T, trackType, genType, randomSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom);

  // print the function
  fprintf(OUT, "\nINPUT\n\n");
  ch = fgetc(NAMES);
  while (ch != EOF)
  {
    fprintf(OUT, "%c", ch);
    ch = fgetc(NAMES);
  }
  fprintf(OUT, "\n");

  // print the version of Bertini
  printVersion(OUT);

  // close files
  fclose(OUT);
  fclose(NAMES);

  // clear memory
  for (i = W->orig_variables - 1; i >= 0; i--)
    free(name_table[i]);
  free(name_table);

  return;
}

void witnessSupersetOutput(witness_t *W, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints the witness supserset to 'witness_superset'     *
\***************************************************************/
{
  int i, k, codim_index, codim, dim, num_set, num_codim = W->num_codim;
  FILE *OUT = fopen("witness_superset", "w");

  if (MPType == 0)
  { // print in double precision
    point_d dehom;
    init_point_d(dehom, 0);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    { 
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular 
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
          { // find the dehomogenized point
            witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
            // print dehomogenized point
            for (k = 0; k < dehom->size; k++)
            {
              print_d(OUT, 0, &dehom->coord[k]);
              fprintf(OUT, "\n");
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
          }
      }

      // print the singular solutions, if they exist
      if (W->codim[codim_index].num_sing > 0)
      { // print singular 
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == SINGULAR)
          { // find the dehomogenized point
            witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
            // print dehomogenized point
            for (k = 0; k < dehom->size; k++)
            {
              print_d(OUT, 0, &dehom->coord[k]);
              fprintf(OUT, "\n");
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]);
          }
      }
    }
    clear_point_d(dehom);
  }
  else if (MPType == 1)
  { // print in multi precision
    point_mp dehom;
    init_point_mp(dehom, 0);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    {
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
          { // find the dehomogenized point
            witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
            // print dehomogenized point
            for (k = 0; k < dehom->size; k++)
            {
              print_mp(OUT, 0, &dehom->coord[k]);
              fprintf(OUT, "\n");
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
          }
      }

      // print the singular solutions, if they exist
      if (W->codim[codim_index].num_sing > 0)
      { // print singular
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == SINGULAR)
          { // find the dehomogenized point
            witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
            // print dehomogenized point
            for (k = 0; k < dehom->size; k++)
            {
              print_mp(OUT, 0, &dehom->coord[k]);
              fprintf(OUT, "\n");
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]);
          }
      }
    }
    clear_point_mp(dehom);
  }
  else
  { // print answers using AMP
    point_d dehom_d;
    point_mp dehom_mp;
    init_point_d(dehom_d, 0);
    init_point_mp2(dehom_mp, 0, 64);

    for (codim_index = 0; codim_index < num_codim; codim_index++)
    {
      codim = W->codim[codim_index].codim;
      dim = W->orig_variables - codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
      num_set = W->codim[codim_index].num_set;

      fprintf(OUT, "\n----------DIMENSION %d----------\n", dim);
      // print the nonsingular solutions, if they exist
      if (W->codim[codim_index].num_nonsing > 0)
      { // print non-singular
        fprintf(OUT, "\nNONSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == NON_SINGULAR)
          {
            if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
            { // find the dehomogenized point using double precision
              witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom_d->size; k++)
              {
                print_d(OUT, 0, &dehom_d->coord[k]);
                fprintf(OUT, "\n");
              }
            }
            else
            { // find the dehomogenized point using multi precision

              // set the precision so that the calculations are correct
              initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
              setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

              // find the dehomogenized point
              witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

              // print dehomogenized point
              for (k = 0; k < dehom_mp->size; k++)
              {
                print_mp(OUT, 0, &dehom_mp->coord[k]);
                fprintf(OUT, "\n");
              }
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]); // better be == 1!!!!
          }
      }

      // print the singular solutions, if they exist
      if (W->codim[codim_index].num_sing > 0)
      { // print singular
        fprintf(OUT, "\nSINGULAR SOLUTIONS\n");
        fprintf(OUT, "---------------\n");

        for (i = 0; i < num_set; i++)
          if (W->codim[codim_index].witnessPt_types[i] == SINGULAR)
          {
            if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
            { // find the dehomogenized point using double precision
              witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
              // print dehomogenized point
              for (k = 0; k < dehom_d->size; k++)
              {
                print_d(OUT, 0, &dehom_d->coord[k]);
                fprintf(OUT, "\n");
              }
            }
            else
            { // find the dehomogenized point using multi precision

              // set the precision so that the calculations are correct
              initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
              setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);

              // find the dehomogenized point
              witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);

              // print dehomogenized point
              for (k = 0; k < dehom_mp->size; k++)
              {
                print_mp(OUT, 0, &dehom_mp->coord[k]);
                fprintf(OUT, "\n");
              }
            }
            fprintf(OUT, "Multiplicity: %d\n\n", W->codim[codim_index].multiplicities[i]);
          }
      }
    }
    clear_point_d(dehom_d);
    clear_point_mp(dehom_mp);
  }

  fprintf(OUT, "\n");
  fclose(OUT);

  return;
}

void witnessFindDehom_d(point_d dehom_d, point_d P_d, witness_t *W, int codim_index)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the dehom coordinates (user input)  *
\***************************************************************/
{
  int i;
  point_d orig_vars;
  init_point_d(orig_vars, 0);

  if (W->PPD.num_var_gp)
  { // the variable group was un-homogenized - find orig_vars == P_d
    point_cp_d(orig_vars, P_d);
  }
  else
  { // the variable group was homogenized - remove the intrinic dehom coordinate and then find orig_vars
    comp_d dehomCoord;

    increase_size_point_d(orig_vars, P_d->size);
    orig_vars->size = P_d->size;
    // find the 'intrinsic' dehom coordinate
    set_d(dehomCoord, W->codim[codim_index].homVarConst_d);
    for (i = 0; i < P_d->size; i++)
    {
      sum_mul_d(dehomCoord, &W->codim[codim_index].H_d->coord[i], &P_d->coord[i]);
    }

    // make sure dehomCoord is a valid number
    if (isnan(dehomCoord->r) || isnan(dehomCoord->i) || isinf(dehomCoord->r) || isinf(dehomCoord->i))
    { // not a valid number
      for (i = 0; i < orig_vars->size; i++)
      {
        set_double_d(&orig_vars->coord[i], -1e199, -1e199);
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
      for (i = 0; i < orig_vars->size; i++)
      { // multiply by recip to produce de-hom coord.
        mul_d(&orig_vars->coord[i], &P_d->coord[i], dehomCoord);
      }
    }
  }

  // find dehom_d using orig_vars
  getDehomPoint_d(dehom_d, orig_vars, orig_vars->size, &W->PPD);

  // clear
  clear_point_d(orig_vars);

  return;
}

void witnessFindDehom_mp(point_mp dehom_mp, point_mp P_mp, witness_t *W, int codim_index, int P_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: finds the point in the dehom coordinates (user input)  *
\***************************************************************/
{
  int i;
  point_mp orig_vars;
  init_point_mp2(orig_vars, 0, P_prec);

  if (W->PPD.num_var_gp)
  { // the variable group was un-homogenized - find orig_vars == P_d
    point_cp_mp(orig_vars, P_mp);
  }
  else
  { // the variable group was homogenized - remove the intrinic dehom coordinate and then find orig_vars
    comp_mp dehomCoord;
    init_mp2(dehomCoord, P_prec);

    increase_size_point_mp(orig_vars, P_mp->size);
    orig_vars->size = P_mp->size;
    // find the 'intrinsic' dehom coordinate
    set_mp(dehomCoord, W->codim[codim_index].homVarConst_mp);
    for (i = 0; i < P_mp->size; i++)
    {
      sum_mul_mp(dehomCoord, &W->codim[codim_index].H_mp->coord[i], &P_mp->coord[i]);
    }

    // make sure dehomCoord is a valid number
    if (!mpfr_number_p(dehomCoord->r) || !mpfr_number_p(dehomCoord->i))
    { // not a valid number
      for (i = 0; i < orig_vars->size; i++)
      {
        set_double_mp(&orig_vars->coord[i], -1e199, -1e199);
      }
    }
    else
    { // we have a number - determine if it is 0
      if (mpfr_zero_p(dehomCoord->r) && mpfr_zero_p(dehomCoord->i))
      { // generate a random perturbation so that we can divide
        mpf_t epsilon;
        mpf_init2(epsilon, P_prec);

        get_comp_rand_mp(dehomCoord);
        mpfr_ui_pow_ui(epsilon, 10, prec_to_digits(P_prec), __gmp_default_rounding_mode);
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
      for (i = 0; i < orig_vars->size; i++)
      { // multiply by recip to produce de-hom coord.
        mul_mp(&orig_vars->coord[i], &P_mp->coord[i], dehomCoord);
      }
    }

    clear_mp(dehomCoord);
  }

  // find dehom_mp using orig_vars
  getDehomPoint_mp(dehom_mp, orig_vars, orig_vars->size, &W->PPD);

  clear_point_mp(orig_vars);

  return;
}

void deflate_for_junkRemoval(prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, membership_slice_moving_t *sliceMover, witness_t *W, int codim_index, int component_number, tracker_config_t *T, FILE *OUT)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup fullRankProgs and endPts to be reduced           *
*  Assume that codim_index has only true witness points         *
*  If component_number == -1, setup all, otherwise component    *
\***************************************************************/
{
  int i, j, k, rV, input_prec, output_prec, newRandomizedProg, randomizedProgUsed = 0, num_set = 0, num_sing = 0, validComponent = 0;
  point_data_d inputPD_d, outputPD_d;
  point_data_mp inputPD_mp, outputPD_mp;
  prog_t *randomizedProg = NULL;
  prog_t *tempProgs = NULL;

  // setup num_set & num_sing
  if (0 <= component_number && component_number < W->codim[codim_index].num_components)
  { // we have a valid component -- setup num_set & num_sing based on this component
    validComponent = 1;
    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if (W->codim[codim_index].component_nums[i] == component_number)
      {
        num_set++;
        if (W->codim[codim_index].multiplicities[i] > 1)
          num_sing++;
      }
  }
  else
  { // setup for whole codim
    num_set = W->codim[codim_index].num_set;
    num_sing = W->codim[codim_index].num_sing;
  }

  // allocate memory
  tempProgs = (prog_t *)bmalloc(num_sing * sizeof(prog_t));

  init_point_data_d(&inputPD_d, 0);
  init_point_data_d(&outputPD_d, 0);
  init_point_data_mp(&inputPD_mp, 0);
  init_point_data_mp(&outputPD_mp, 0);

  // setup the randomized SLP - does [I A]f
  if (W->codim[codim_index].A_rows > 0 && W->codim[codim_index].A_cols > 0)
  { // need to randomize
    randomizedProg = (prog_t *)bmalloc(1 * sizeof(prog_t));
    randomize_SLP(randomizedProg, W->Prog, W->codim[codim_index].A_rat, W->codim[codim_index].A_rows, W->codim[codim_index].A_cols);
    // a new randomized prog was made
    newRandomizedProg = 1;
  }
  else
  { // point to original SLP 
    randomizedProg = W->Prog;
    // a new randomized prog was not made
    newRandomizedProg = 0;
  }

  // everything are witness points - so we can setup fullRankProgs and endPts for each one
  *fullRankProgs = (prog_t **)bmalloc(num_set * sizeof(prog_t *));
  *fullRankProgInfo = (int *)bmalloc(num_set * sizeof(int));
  if (T->MPType == 0)
    *endPts_d = (endpoint_data_d *)bmalloc(num_set * sizeof(endpoint_data_d));
  else if (T->MPType == 1)
    *endPts_mp = (endpoint_data_mp *)bmalloc(num_set * sizeof(endpoint_data_mp));
  else
    *endPts_amp = (endpoint_data_amp *)bmalloc(num_set * sizeof(endpoint_data_amp));

  i = j = 0;
  for (k = 0; k < W->codim[codim_index].num_set; k++)
    if (validComponent == 0 || W->codim[codim_index].component_nums[k] == component_number)
    { // setup progs[j] & endPts[j]
      if (T->MPType == 0)
      { // copy all of the endPt data
        init_endpoint_data_d(&(*endPts_d)[j]);
        endpoint_data_cp_d(&(*endPts_d)[j], &W->codim[codim_index].witnessPts_d[k]);

        // see if we need to deflate
        if (W->codim[codim_index].witnessPt_types[k] == NON_SINGULAR)
        { // simply point to the randomized prog
          (*fullRankProgs)[j] = randomizedProg;
          // setup Info
          if (newRandomizedProg && !randomizedProgUsed) 
          { // we have a new SLP that has not been used before - so it will need cleared
            (*fullRankProgInfo)[j] = 1;
          }
          else
          { // this does not need cleared out
            (*fullRankProgInfo)[j] = 0;
          }
          // we have used the randomized SLP
          randomizedProgUsed = 1;

          // setup the number of deflations needed == 0
          W->codim[codim_index].deflations_needed[k] = 0;
        }
        else
        { // deflate

          // setup inputPD_d
          point_cp_d(inputPD_d.point, W->codim[codim_index].witnessPts_d[k].endPt);
          set_d(inputPD_d.time, W->codim[codim_index].witnessPts_d[k].finalT);
          input_prec = 52;

          // do the deflation
          rV = deflation(&W->codim[codim_index].deflations_needed[k], &tempProgs[i], &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, (*endPts_d)[j].corank, (*endPts_d)[j].smallest_nonzero_SV, (*endPts_d)[j].largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, W->codim[codim_index].witnessPts_d[k].last_approx, NULL, 52, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[k] - 1);

          // setup fullRankProgs & endPts
          (*fullRankProgs)[j] = &tempProgs[i];
          point_cp_d((*endPts_d)[j].endPt, outputPD_d.point);

          // check for deflation success
          if (rV == 0)
            (*fullRankProgInfo)[j] = 1; // deflated properly
          else
          {
            (*fullRankProgInfo)[j] = -1; // did not deflate properly
            printf("\nNOTE: A witness point did not deflate properly.\n");
          }

          // increment i
          i++;
        }
      }
      else if (T->MPType == 1)
      {
        init_endpoint_data_mp(&(*endPts_mp)[j]);
        endpoint_data_cp_mp(&(*endPts_mp)[j], &W->codim[codim_index].witnessPts_mp[k]);

        // see if we need to deflate
        if (W->codim[codim_index].witnessPt_types[k] == NON_SINGULAR)
        { // simply point to the randomized prog
          (*fullRankProgs)[j] = randomizedProg;
          // setup Info
          if (newRandomizedProg && !randomizedProgUsed)
          { // we have a new SLP that has not been used before - so it will need cleared
            (*fullRankProgInfo)[j] = 1;
          }
          else
          { // this does not need cleared out
            (*fullRankProgInfo)[j] = 0;
          }
          // we have used the randomized SLP
          randomizedProgUsed = 1; 

          // setup the number of deflations needed == 0
          W->codim[codim_index].deflations_needed[k] = 0;
        }
        else
        { // deflate

          // setup inputPD_mp
          point_cp_mp(inputPD_mp.point, W->codim[codim_index].witnessPts_mp[k].endPt);
          set_mp(inputPD_mp.time, W->codim[codim_index].witnessPts_mp[k].finalT);
          input_prec = T->Precision;

          // do the deflation
          rV = deflation(&W->codim[codim_index].deflations_needed[k], &tempProgs[i], &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, (*endPts_mp)[j].corank, (*endPts_mp)[j].smallest_nonzero_SV, (*endPts_mp)[j].largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, NULL, W->codim[codim_index].witnessPts_mp[k].last_approx, T->Precision, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[k] - 1);

          // setup fullRankProgs & endPts
          (*fullRankProgs)[j] = &tempProgs[i];
          point_cp_mp((*endPts_mp)[j].endPt, outputPD_mp.point);

          // check for deflation success
          if (rV == 0)
            (*fullRankProgInfo)[j] = 1; // deflated properly
          else
          {
            (*fullRankProgInfo)[j] = -1; // did not deflate properly
            printf("NOTE: A witness point did not deflate properly.\n");
            printf("      Consider tightening the final tolerance and/or use adaptive precision.\n");
          }

          // increment i
          i++;
        }
      }
      else
      {
        init_endpoint_data_amp(&(*endPts_amp)[j], W->codim[codim_index].witnessPts_amp[k].curr_prec, W->codim[codim_index].witnessPts_amp[k].last_approx_prec);
        endpoint_data_cp_amp(&(*endPts_amp)[j], &W->codim[codim_index].witnessPts_amp[k]);
 
        // see if we need to deflate
        if (W->codim[codim_index].witnessPt_types[k] == NON_SINGULAR)
        { // simply point to the randomized prog
          (*fullRankProgs)[j] = randomizedProg;
          // setup Info
          if (newRandomizedProg && !randomizedProgUsed)
          { // we have a new SLP that has not been used before - so it will need cleared
            (*fullRankProgInfo)[j] = 1;
          }
          else
          { // this does not need cleared out
            (*fullRankProgInfo)[j] = 0;
          }
          // we have used the randomized SLP
          randomizedProgUsed = 1; 

          // setup the number of deflations needed == 0
          W->codim[codim_index].deflations_needed[k] = 0;
        }
        else
        { // deflate

          // setup inputPD
          input_prec = W->codim[codim_index].witnessPts_amp[k].curr_prec;
          if (input_prec < 64)
          {
            point_cp_d(inputPD_d.point, W->codim[codim_index].witnessPts_amp[k].endPt_d);
            set_d(inputPD_d.time, W->codim[codim_index].witnessPts_amp[k].finalT_d);
          }
          else
          {
            setprec_point_mp(inputPD_mp.point, input_prec);
            point_cp_mp(inputPD_mp.point, W->codim[codim_index].witnessPts_amp[k].endPt_mp);

            setprec_mp(inputPD_mp.time, input_prec);
            set_mp(inputPD_mp.time, W->codim[codim_index].witnessPts_amp[k].finalT_mp);
          }

          // do the deflation
          rV = deflation(&W->codim[codim_index].deflations_needed[k], &tempProgs[i], &outputPD_d, &outputPD_mp, &output_prec, randomizedProg, (*endPts_amp)[j].corank, (*endPts_amp)[j].smallest_nonzero_SV, (*endPts_amp)[j].largest_zero_SV, &inputPD_d, &inputPD_mp, input_prec, W->codim[codim_index].witnessPts_amp[k].last_approx_d, W->codim[codim_index].witnessPts_amp[k].last_approx_mp, W->codim[codim_index].witnessPts_amp[k].last_approx_prec, sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols, T, OUT, W->codim[codim_index].multiplicities[k] - 1);

          // setup fullRankProgs & endPts
          (*fullRankProgs)[j] = &tempProgs[i];
          if (output_prec < 64)
          {
            point_cp_d((*endPts_amp)[j].endPt_d, outputPD_d.point);
          }
          else
          { // set precision correctly
            if ((*endPts_amp)[j].curr_prec != output_prec)
            {
              setprec_point_mp((*endPts_amp)[j].endPt_mp, output_prec);
            }
            point_cp_mp((*endPts_amp)[j].endPt_mp, outputPD_mp.point);
          }

          // check for deflation success
          if (rV == 0)
            (*fullRankProgInfo)[j] = 1; // deflated properly
          else
          {
            (*fullRankProgInfo)[j] = -1; // did not deflate properly
            printf("\nNOTE: A witness point did not deflate properly.\n");
            printf("      Consider tightening the final tolerance.\n");
          }

          // increment i
          i++;
        }
      }
      // increment j
      j++;
    }

  // clear data
  clear_point_data_d(&inputPD_d);
  clear_point_data_d(&outputPD_d);
  clear_point_data_mp(&inputPD_mp);
  clear_point_data_mp(&outputPD_mp);

  // see if we can clear randomizedProg
  if (newRandomizedProg && !randomizedProgUsed)
  { // clear out randomizedProg
    clearProg(randomizedProg, T->MPType, 0);
    free(randomizedProg);
  }
  else
  { // set to NULL
    randomizedProg = NULL;
  }

  tempProgs = NULL; // the memory is pointed to using fullRankProgs and will be cleared later

  return;
}

void sort_endpoint_data(int *mult, endpoint_data_d *Pts_d, endpoint_data_mp *Pts_mp, endpoint_data_amp *Pts_amp, int numPoints, int MPType, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: sort the endpoint data                                 *
\***************************************************************/
{
  int i, j, size, cont, indexI = 0, indexJ = 0;
  mpf_t norm_diff;

  mpf_init(norm_diff);

  // initialize mult to all zeros
  for (i = 0; i < numPoints; i++)
    mult[i] = 0;

  if (MPType == 0)
  { // look at _d
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(numPoints * sizeof(sortStruct_d));
    for (i = 0; i < numPoints; i++)
    { // initialize sortPts[i]
      sortPts[i].path_num = i;
      // find norm using double precision
      size = Pts_d[i].endPt->size;
      sortPts[i].norm = 0;
      for (j = 0; j < size; j++)
        sortPts[i].norm += norm_sqr_d(&Pts_d[i].endPt->coord[j]);
      sortPts[i].norm = sqrt(sortPts[i].norm);
    }

    // sort
    qsort(sortPts, numPoints, sizeof(sortStruct_d), sort_order_d);

    // loop through to check the points
    for (i = 0; i < numPoints; i++)
    { // find the index
      indexI = sortPts[i].path_num;
      // see if we need to classify this one
      if (!mult[indexI])
      { // setup mult
        mult[indexI] = 1;

        // loop through to find the multiplicity
        cont = 1;
        j = i;
        do
        { // increment the counter (start at i + 1)
          j++;

          if (j < numPoints)
          { // check norm diff
            indexJ = sortPts[j].path_num;
            if (sortPts[j].norm - sortPts[i].norm > tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we should compare
          if (cont && !mult[indexJ])
          { // find the difference
            findDiff_point(norm_diff, Pts_d[indexI].endPt, NULL, 52, Pts_d[indexJ].endPt, NULL, 52);
            if (mpf_get_d(norm_diff) < tol)
            { // i & j are equal
              mult[indexI]++;    // increase multiplicity
              mult[indexJ] = -1; // j is dead
            }
          }
        } while (cont);
      }
    }
    // clear memory
    free(sortPts);
  }
  else if (MPType == 1)
  { // look at _mp
    sortStruct_mp *sortPts = (sortStruct_mp *)bmalloc(numPoints * sizeof(sortStruct_mp));
    for (i = 0; i < numPoints; i++)
    { // initialize sortPts[i]
      sortPts[i].path_num = i;
      // find norm using fixed precision
      size = Pts_mp[i].endPt->size;
      mpf_init(sortPts[i].norm);
      mpf_set_ui(sortPts[i].norm, 0);
      for (j = 0; j < size; j++)
      {
        norm_sqr_mp(norm_diff , &Pts_mp[i].endPt->coord[j]);
        mpf_add(sortPts[i].norm, sortPts[i].norm, norm_diff);
      }
      mpf_sqrt(sortPts[i].norm, sortPts[i].norm);
    }

    // sort
    qsort(sortPts, numPoints, sizeof(sortStruct_mp), sort_order_mp);

    // loop through to check the points
    for (i = 0; i < numPoints; i++)
    { // find the index
      indexI = sortPts[i].path_num;
      // see if we need to classify this one
      if (!mult[indexI])
      { // setup mult
        mult[indexI] = 1;

        // loop through to find the multiplicity
        cont = 1;
        j = i;
        do
        { // increment the counter (start at i + 1)
          j++;

          if (j < numPoints)
          { // check norm diff
            indexJ = sortPts[j].path_num;
            mpf_sub(norm_diff, sortPts[j].norm, sortPts[i].norm);
            if (mpf_get_d(norm_diff) > tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we should compare
          if (cont && !mult[indexJ])
          { // find the difference
            findDiff_point(norm_diff, NULL, Pts_mp[indexI].endPt, mpf_get_default_prec(), NULL, Pts_mp[indexJ].endPt, mpf_get_default_prec());
            if (mpf_get_d(norm_diff) < tol)
            { // i & j are equal
              mult[indexI]++;    // increase multiplicity
              mult[indexJ] = -1; // j is dead
            }
          }
        } while (cont);
      }
    }
    // clear memory
    for (i = numPoints - 1; i >= 0; i--)
      mpf_clear(sortPts[i].norm);
    free(sortPts);
  }
  else
  { // look at _amp
    sortStruct_d *sortPts = (sortStruct_d *)bmalloc(numPoints * sizeof(sortStruct_d));
    mpf_t tempNorm;
    mpf_init(tempNorm);

    for (i = 0; i < numPoints; i++)
    { // initialize sortPts[i]
      sortPts[i].path_num = i;
      // find norm
      sortPts[i].norm = 0;
      if (Pts_amp[i].curr_prec < 64)
      { // use _d
        size = Pts_amp[i].endPt_d->size;
        for (j = 0; j < size; j++)
          sortPts[i].norm += norm_sqr_d(&Pts_amp[i].endPt_d->coord[j]);
        sortPts[i].norm = sqrt(sortPts[i].norm);
      }
      else
      { // use _mp
        mpf_set_prec(norm_diff, Pts_amp[i].curr_prec);
        mpf_set_prec(tempNorm, Pts_amp[i].curr_prec);
        mpf_set_ui(tempNorm, 0);

        size = Pts_amp[i].endPt_mp->size;
        for (j = 0; j < size; j++)
        {
          norm_sqr_mp(norm_diff , &Pts_amp[i].endPt_mp->coord[j]);
          mpf_add(tempNorm, tempNorm, norm_diff);
        }
        mpf_sqrt(tempNorm, tempNorm);
        sortPts[i].norm = mpf_get_d(tempNorm);
      }
    }

    // sort
    qsort(sortPts, numPoints, sizeof(sortStruct_d), sort_order_d);

    // loop through to check the points
    for (i = 0; i < numPoints; i++)
    { // find the index
      indexI = sortPts[i].path_num;
      // see if we need to classify this one
      if (!mult[indexI])
      { // setup mult
        mult[indexI] = 1;

        // loop through to find the multiplicity
        cont = 1;
        j = i;
        do
        { // increment the counter (start at i + 1)
          j++;

          if (j < numPoints)
          { // check norm diff
            indexJ = sortPts[j].path_num;
            if (sortPts[j].norm - sortPts[i].norm > tol)
              cont = 0;
          }
          else
            cont = 0;

          // check to see if we should compare
          if (cont && !mult[indexJ])
          { // find the difference
            findDiff_point(norm_diff, Pts_amp[indexI].endPt_d, Pts_amp[indexI].endPt_mp, Pts_amp[indexI].curr_prec, Pts_amp[indexJ].endPt_d, Pts_amp[indexJ].endPt_mp, Pts_amp[indexJ].curr_prec);
            if (mpf_get_d(norm_diff) < tol)
            { // i & j are equal
              mult[indexI]++;    // increase multiplicity
              mult[indexJ] = -1; // j is dead
            }
          }
        } while (cont);
      }
    }
    // free memory
    mpf_clear(tempNorm);
    free(sortPts);
  }

  mpf_clear(norm_diff);

  return;
}

void multiplicity_witness(witness_t *W, int codim_index, int MPType, double tol)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: find multiplicity for the witness set and remove the   *
* multiple points                                               *
\***************************************************************/
{
  int i, j, count, num_points = W->codim[codim_index].num_set;
  int *mult = (int *)bmalloc(num_points * sizeof(int)), *temp_types = NULL;
  endpoint_data_d *tempPts_d = NULL;
  endpoint_data_mp *tempPts_mp = NULL;
  endpoint_data_amp *tempPts_amp = NULL;

  // initialize mult to all zeros
  for (i = 0; i < num_points; i++)
    mult[i] = 0;

  // perform the sortin
  sort_endpoint_data(mult, W->codim[codim_index].witnessPts_d, W->codim[codim_index].witnessPts_mp, W->codim[codim_index].witnessPts_amp, W->codim[codim_index].num_set, MPType, tol);

  // count the number of points to keep
  count = 0;
  for (i = 0; i < num_points; i++)
    if (mult[i] > 0)
      count++;

  // remove the extra multiple points and setup multiplicities
  W->codim[codim_index].multiplicities = (int *)bmalloc(count * sizeof(int));

  if (count == num_points)
  { // leave all alone
    for (i = 0; i < num_points; i++)
      W->codim[codim_index].multiplicities[i] = mult[i];
  }
  else
  { // some need removed
    temp_types = W->codim[codim_index].witnessPt_types;
    W->codim[codim_index].witnessPt_types = (int *)bmalloc(count * sizeof(int));

    if (MPType == 0)
    { // remove extra from witnessPts_d 
      tempPts_d = W->codim[codim_index].witnessPts_d;
      W->codim[codim_index].witnessPts_d = (endpoint_data_d *)bmalloc(count * sizeof(endpoint_data_d));

      // copy back the good points
      j = 0;
      for (i = 0; i < num_points; i++)
      { // see if it needs copied over
        if (mult[i] > 0)
        { // copy this over
          init_endpoint_data_d(&W->codim[codim_index].witnessPts_d[j]);
          endpoint_data_cp_d(&W->codim[codim_index].witnessPts_d[j], &tempPts_d[i]);
          W->codim[codim_index].witnessPt_types[j] = temp_types[i];
          W->codim[codim_index].multiplicities[j] = mult[i];
          // increment j
          j++;
        }
        // clear old point
        clear_endpoint_data_d(&tempPts_d[i]);
      }
      // clear tempPts_d
      free(tempPts_d);
    }
    else if (MPType == 1)
    { // remove extra from witnessPts_mp
      tempPts_mp = W->codim[codim_index].witnessPts_mp;
      W->codim[codim_index].witnessPts_mp = (endpoint_data_mp *)bmalloc(count * sizeof(endpoint_data_mp));

      // copy back the good points
      j = 0;
      for (i = 0; i < num_points; i++)
      { // see if it needs copied over
        if (mult[i] > 0)
        { // copy this over
          init_endpoint_data_mp(&W->codim[codim_index].witnessPts_mp[j]);
          endpoint_data_cp_mp(&W->codim[codim_index].witnessPts_mp[j], &tempPts_mp[i]);
          W->codim[codim_index].witnessPt_types[j] = temp_types[i];
          W->codim[codim_index].multiplicities[j] = mult[i];
          // increment j
          j++;
        }
        // clear old point
        clear_endpoint_data_mp(&tempPts_mp[i]);
      }
      // clear tempPts_mp
      free(tempPts_mp);
    }
    else
    { // remove extra from witnessPts_amp
      tempPts_amp = W->codim[codim_index].witnessPts_amp;
      W->codim[codim_index].witnessPts_amp = (endpoint_data_amp *)bmalloc(count * sizeof(endpoint_data_amp));

      // copy back the good points
      j = 0;
      for (i = 0; i < num_points; i++)
      { // see if needs copied over
        if (mult[i] > 0)
        { // copy this over
          init_endpoint_data_amp(&W->codim[codim_index].witnessPts_amp[j], tempPts_amp[i].curr_prec, tempPts_amp[i].last_approx_prec);
          endpoint_data_cp_amp(&W->codim[codim_index].witnessPts_amp[j], &tempPts_amp[i]);
          W->codim[codim_index].witnessPt_types[j] = temp_types[i];
          W->codim[codim_index].multiplicities[j] = mult[i];
          // increment j
          j++;
        }
        // clear old point
        clear_endpoint_data_amp(&tempPts_amp[i]);
      }
      // clear tempPts_amp
      free(tempPts_amp);
    }

    // update counts
    j = 0;
    for (i = 0; i < count; i++)
      if (W->codim[codim_index].witnessPt_types[i] == SINGULAR)
        j++;

    W->codim[codim_index].num_set = count;
    W->codim[codim_index].num_sing = j;
    W->codim[codim_index].num_nonsing = count - j;
  }

  free(mult);

  return;
}

void clear_sliceMover_fullRankProgs(membership_slice_moving_t **sliceMover, prog_t ****fullRankProgs, int ***fullRankProgInfo, endpoint_data_d ***endPts_d, endpoint_data_mp ***endPts_mp, endpoint_data_amp ***endPts_amp, witness_t *W, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear fullRankProgs and endPts                         *
\***************************************************************/
{
  int i, num_codim = W->num_codim;

  // clear each codim
  for (i = num_codim - 1; i >= 0; i--)
  {
    if (MPType == 0)
    {
      clear_sliceMover_fullRankProgs_codim(&(*sliceMover)[i], &(*fullRankProgs)[i], &(*fullRankProgInfo)[i], &(*endPts_d)[i], NULL, NULL, W->codim[i].num_set, MPType);
      free((*endPts_d)[i]);
    }
    else if (MPType == 1)
    {
      clear_sliceMover_fullRankProgs_codim(&(*sliceMover)[i], &(*fullRankProgs)[i], &(*fullRankProgInfo)[i], NULL, &(*endPts_mp)[i], NULL, W->codim[i].num_set, MPType);

      free((*endPts_mp)[i]);
    }
    else
    {
      clear_sliceMover_fullRankProgs_codim(&(*sliceMover)[i], &(*fullRankProgs)[i], &(*fullRankProgInfo)[i], NULL, NULL, &(*endPts_amp)[i], W->codim[i].num_set, MPType);

      free((*endPts_amp)[i]);
    }
    free((*fullRankProgs)[i]);
    free((*fullRankProgInfo)[i]);
  }

  if (MPType == 0)
    free(*endPts_d);
  else if (MPType == 1)
    free(*endPts_mp);
  else
    free(*endPts_amp);

  free(*sliceMover);
  free(*fullRankProgs);
  free(*fullRankProgInfo);

  return;
}

void clear_sliceMover_fullRankProgs_codim(membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, int num_set, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear fullRankProgs and endPts                         *
\***************************************************************/
{
  int i;

  // clear for each point
  for (i = num_set - 1; i >= 0; i--)
  { // clear fullRankProgs[i] && endPts
    if (MPType == 0)
    { // clear endPts_d[i]
      clear_fullRankProg_endPt(&(*fullRankProgs)[i], (*fullRankProgInfo)[i], &(*endPts_d)[i], NULL, NULL, MPType);
    }
    else if (MPType == 1)
    { // clear endPts_mp[i]
      clear_fullRankProg_endPt(&(*fullRankProgs)[i], (*fullRankProgInfo)[i], NULL, &(*endPts_mp)[i], NULL, MPType);
    }
    else
    { // clear endPts_amp[i]
      clear_fullRankProg_endPt(&(*fullRankProgs)[i], (*fullRankProgInfo)[i], NULL, NULL, &(*endPts_amp)[i], MPType);
    }
  }

  // clear sliceMover
  clear_sliceMover(sliceMover, MPType);

  return;
}

void clear_fullRankProg_endPt(prog_t **fullRankProg, int fullRankProgInfo, endpoint_data_d *endPt_d, endpoint_data_mp *endPt_mp, endpoint_data_amp *endPt_amp, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clear fullRankProg and endPt                           *
\***************************************************************/
{
  if (MPType == 0)
  { // clear endPt_d
    clear_endpoint_data_d(endPt_d);
  }
  else if (MPType == 1)
  { // clear endPt_mp
    clear_endpoint_data_mp(endPt_mp);
  }
  else
  { // clear endPt_amp
    clear_endpoint_data_amp(endPt_amp);
  }

  // clear fullRankProg, if needed
  if (fullRankProgInfo && *fullRankProg != NULL)
  { // clear the SLP
    clearProg(*fullRankProg, MPType, 0);
  }
  else
  { // set to NULL
    *fullRankProg = NULL;
  }

  return;
}

void clear_sliceMover(membership_slice_moving_t *sliceMover, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: clears sliceMover                                      *
\***************************************************************/
{
  // set Prog to NULL
  sliceMover->Prog = NULL;

  // clear gamma
  clear_d_mp_rat(sliceMover->gamma_d, sliceMover->gamma_mp, sliceMover->gamma_rat, MPType);

  // clear A
  clear_mat(sliceMover->A_d, sliceMover->A_mp, sliceMover->A_rat, MPType);

  // clear B
  clear_mat(sliceMover->B_d, sliceMover->B_mp, sliceMover->B_rat, MPType);

  // clear p
  clear_vec(sliceMover->p_d, sliceMover->p_mp, sliceMover->p_rat, MPType);
  
  // clear startSliceVec
  if (sliceMover->startSliceVec_init)
  {
    clear_vec(sliceMover->startSliceVec_d, sliceMover->startSliceVec_mp, sliceMover->startSliceVec_rat, MPType);
  }

  // clear targetSliceVec
  if (sliceMover->targetSliceVec_init)
  {
    clear_vec(sliceMover->targetSliceVec_d, sliceMover->targetSliceVec_mp, sliceMover->targetSliceVec_rat, MPType);
  }

  // clear K_rat
  clear_mat_rat(sliceMover->K_rat, sliceMover->K_rows, sliceMover->K_cols);  

  return;
}

void displayDeflationSummary(int **fullRankProgInfo, witness_t *W)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: displays deflation summary                             *
\***************************************************************/
{
  int i, j, total_count = 0, total_error = 0;
  int *def_count = (int *)bmalloc(W->num_codim * sizeof(int));
  int *def_error = (int *)bmalloc(W->num_codim * sizeof(int));

  // count how many points were deflated and how many failed
  for (i = 0; i < W->num_codim; i++)
  { // initialize counts
    def_count[i] = def_error[i] = 0;

    for (j = 0; j < W->codim[i].num_set; j++)
      if (W->codim[i].deflations_needed[j] > 0)
      { // deflation was used
        def_count[i]++;

        if (fullRankProgInfo[i][j] == -1)
        { // deflation failed
          def_error[i]++;
        }
      }

    // update totals
    total_count += def_count[i];
    total_error += def_error[i];
  }

  if (total_count > 0)
    printf("Witness Points Deflated: %d\n", total_count);

  if (total_error > 0)
  { // print error summary
    printf(" Deflation failures: %d - try adjusting FinalTol and/or RatioTolerance.\n", total_error);
  }

  free(def_count);
  free(def_error);

  return;
}

