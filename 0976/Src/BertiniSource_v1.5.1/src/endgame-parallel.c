// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"

int endgame_d(int endgameNumber, int pathNum, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame  
	    retVal = PSEG_d(pathNum, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func_d, find_dehom);
            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_d(pathNum, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func_d, find_dehom);
            break;
    default:
	    printf("ERROR: Inappropriate endgame number chosen.\n");
            bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}

int endgame_rank_d(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_d *Final, point_d last_approx, point_data_d *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame
            retVal = PSEG_rank_d(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func_d, find_dehom);
            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_rank_d(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func_d, find_dehom);
            break;
    default:
            printf("ERROR: Inappropriate endgame number chosen.\n");
            bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}

//////////////////MP VERSIONS ///////////////////////////

int endgame_mp(int endgameNumber, int pathNum, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame
            retVal = PSEG_mp(pathNum, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func, find_dehom);
            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_mp(pathNum, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func, find_dehom);
            break;
    default:
             printf("ERROR: Inappropriate endgame number chosen.\n");
             bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}

int endgame_rank_mp(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, point_data_mp *Final, point_mp last_approx, point_data_mp *Start, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED, int (*eval_func)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame
            retVal = PSEG_rank_mp(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func, find_dehom);
            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_rank_mp(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, Final, last_approx, Start, T, OUT, midOUT, ED, eval_func, find_dehom);
            break;
    default:
            printf("ERROR: Inappropriate endgame number chosen.\n");
            bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}

//////// Adaptive Precision //////////

int endgame_amp(int endgameNumber, int pathNum, int *prec, double *first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame
//            #ifdef _REACTIVE
//              retVal = PSEG_reactive_amp(pathNum, prec, first_increase, Final_d, Final_mp, Start_d, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp);
//            #else
              retVal = PSEG_amp(pathNum, prec, first_increase, Final_d, Final_mp, last_approx_prec, last_approx_d, last_approx_mp, Start_d, NULL, 52, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
//            #endif

            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_amp(pathNum, prec, first_increase, Final_d, Final_mp, last_approx_prec, last_approx_d, last_approx_mp, Start_d, NULL, 52, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
            break;
    default:
            printf("ERROR: Inappropriate endgame number chosen.\n");
            bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}

int endgame_rank_amp(int endgameNumber, int pathNum, double *condNum, int rankType, int *rankDef, int *corank, double *smallest_nonzero_SV, double *largest_zero_SV, int *prec, double *first_increase, point_data_d *Final_d, point_data_mp *Final_mp, int *last_approx_prec, point_d last_approx_d, point_mp last_approx_mp, point_data_d *Start_d, tracker_config_t *T, FILE *OUT, FILE *midOUT, void const *ED_d, void const *ED_mp, int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *), int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *), int (*change_prec)(void const *, int), int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int retVal = 0;

  switch (endgameNumber)
  {
    case 1: // Power series endgame
            retVal = PSEG_rank_amp(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, prec, first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, Start_d, NULL, 52, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
            break;
    case 2: // Cauchy endgame
    case 3: // track-back endgame - somehow got through
            retVal = CauchyEG_rank_amp(pathNum, condNum, rankType, rankDef, corank, smallest_nonzero_SV, largest_zero_SV, prec, first_increase, Final_d, Final_mp, last_approx_d, last_approx_mp, last_approx_prec, Start_d, NULL, 52, T, OUT, midOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
            break;
    default:
            printf("ERROR: Inappropriate endgame number chosen.\n");
            bexit(ERROR_CONFIGURATION);
  }

  return retVal;
}


