// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"
#include "localdim.h"

// This file contains the junk removal process using the local dimension to determine if the point is junk or not
int junkRemoval_is_isolated(witness_t *W, int pathNum, int pathNum_codim_index,  tracker_config_t *T);

int junkRemoval_ldt(witness_t *W, int pathNum, int pathNum_codim_index, membership_slice_moving_t *sliceMover, prog_t ***fullRankProgs, int **fullRankProgInfo, endpoint_data_d **endPts_d, endpoint_data_mp **endPts_mp, endpoint_data_amp **endPts_amp, tracker_config_t *T, FILE *OUT, char *midName, int specificCodim, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: perform junk removal using local dimension test        *
\***************************************************************/
{
  int isJunk = 0;

  // do the local dimension test
  isJunk = junkRemoval_is_isolated(W, pathNum, pathNum_codim_index, T);

  // see if we have a satisfactory answer
  if (isJunk < 0 && specificCodim == 0)
  { // we need to switch to a membership test
    isJunk = junkRemoval_mem(W, pathNum, pathNum_codim_index, sliceMover, fullRankProgs, fullRankProgInfo, endPts_d, endPts_mp, endPts_amp, T, OUT, midName, my_id, num_processes, headnode);
  }

  return isJunk;
}

int junkRemoval_is_isolated(witness_t *W, int pathNum, int pathNum_codim_index,  tracker_config_t *T)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES: 1 - is junk (not isolated), -1 - unknown       *
* NOTES: use local dimension test                               *
\***************************************************************/
{
  int rV = 0, isJunk = 0, mult = 0, corank = 0;

  // do the first set of testing based on the corank of the Jacobian (1st order LDT)
  // Junk = 1 if mult < 1 + corank
  mult = W->codim[pathNum_codim_index].multiplicities[pathNum];
  if (T->MPType == 0)
    corank = W->codim[pathNum_codim_index].witnessPts_d[pathNum].corank;
  else if (T->MPType == 1)
    corank = W->codim[pathNum_codim_index].witnessPts_mp[pathNum].corank;
  else
    corank = W->codim[pathNum_codim_index].witnessPts_amp[pathNum].corank;

  // see if we have already have enough to classify it as junk or we need to continue with higher order LDT
  if (mult < 1 + corank)
  { // already junk!
    isJunk = 1;
  }
  else
  { // run 'is_isolated' to perform the higher order LDT
    prog_deriv_t deriv;

    // setup deriv from the SLP
    setup_deriv_from_SLP(&deriv, W->Prog);

    // add slices and patches to deriv  
    add_slices_patch_to_deriv(&deriv, W->codim[pathNum_codim_index].B_d, W->codim[pathNum_codim_index].B_mp, W->codim[pathNum_codim_index].B_rat, W->codim[pathNum_codim_index].p_d, W->codim[pathNum_codim_index].p_mp, W->codim[pathNum_codim_index].p_rat, T->MPType);

    // perform the local dimension test
    if (T->MPType == 0)
      rV = is_isolated(&mult, &deriv, W->codim[pathNum_codim_index].witnessPts_d[pathNum].endPt, NULL, 52, W->codim[pathNum_codim_index].witnessPts_d[pathNum].last_approx, NULL, 52, W->codim[pathNum_codim_index].multiplicities[pathNum], T, 0);
    else if (T->MPType == 1)
      rV = is_isolated(&mult, &deriv, NULL, W->codim[pathNum_codim_index].witnessPts_mp[pathNum].endPt, T->Precision, NULL, W->codim[pathNum_codim_index].witnessPts_mp[pathNum].last_approx, T->Precision, W->codim[pathNum_codim_index].multiplicities[pathNum], T, 0);
    else
      rV = is_isolated(&mult, &deriv, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_d, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].endPt_mp, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].curr_prec, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].last_approx_d, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].last_approx_mp, W->codim[pathNum_codim_index].witnessPts_amp[pathNum].last_approx_prec, W->codim[pathNum_codim_index].multiplicities[pathNum], T, 0);

    if (rV < 0) // don't know if junk
      isJunk = rV;
    else // know an answer
      isJunk = !rV;

    // clear deriv
    clear_deriv(&deriv);
  }

  return isJunk;
}



