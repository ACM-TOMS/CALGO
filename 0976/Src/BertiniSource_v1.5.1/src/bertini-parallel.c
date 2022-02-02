// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "parallel.h"
#include "ppParse.h"

int main(int argC, char *args[])
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  bclock_t time1, time2;
  int trackType, genType = 0, MPType, userHom, sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
  int my_id, num_processes, headnode = 0; // headnode is always 0
  double parse_time;
  unsigned int currentSeed;
  char *inputName = NULL, *startName = NULL;

  // initialize MPI and find the number of processes and their id numbers
  initialize_mpi(argC, args, &num_processes, &my_id);

  if (my_id == headnode)
  { // headnode controls the overall execution
    printf("\n   Bertini(TM) v%s", BERTINI_VERSION_STRING);
    printf("\n  (%s)\n\n", BERTINI_DATE_STRING);
    printf(" D.J. Bates, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
    printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());

    // check for write privilege
    if (checkWritePrivilege())
    {
      printf("ERROR: Bertini does not have write privileges!\n");
      bexit(ERROR_WRITE_PRIVILEGE);
    }

    if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--help") || !strcmp(args[1], "-help"))) // help
    { // print information about Bertini
      printf("\nThis is Bertini v%s, developed by Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_VERSION_STRING);
      printf("See http://bertini.nd.edu for details about Bertini.\n\n");
    }
    else if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--version") || !strcmp(args[1], "-version"))) // version
    { // version information already printed -- simply exit

    }
    else
    { // standard execution
      bclock(&time1);

      // setup inputName
      if (argC >= 2 && args[1] != NULL)
      { // inputName is args[1]
        inputName = (char *)bmalloc((strlen(args[1]) + 1) * sizeof(char));
        strcpy(inputName, args[1]);

        // setup startName
        if (argC >= 3 && args[2] != NULL)
        { // startName is args[2]
          startName = (char *)bmalloc((strlen(args[2]) + 1) * sizeof(char));
          strcpy(startName, args[2]);
        }
        else
        { // default to 'start'
          startName = (char *)bmalloc(6 * sizeof(char));
          strcpy(startName, "start");
        }
      }
      else
      { // default to 'input' & 'start'
        inputName = (char *)bmalloc(6 * sizeof(char));
        strcpy(inputName, "input");
        startName = (char *)bmalloc(6 * sizeof(char));
        strcpy(startName, "start");
      }

      // parse the input file and seed the random number generator
      parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);

      // remove the output files from possibly previous runs - delete only files not needed
      remove_output_files(trackType, sharpenOnly, 1);

      bclock(&time2);
      totalTime(&parse_time, time1, time2);

      // call the main functions
      if (sharpenOnly)
      { // sharpen the endpoints from a previous run
        sharpen_process_main(MPType, trackType, currentSeed, my_id, num_processes, headnode);
      }
      else
      { // do either function evaluation, zero dimensional or positive dimensional tracking
        if (trackType == 0)
        { // zero dimensional tracking
          zero_dim_main(MPType, parse_time, currentSeed, startName, my_id, num_processes, headnode);
        }
        else if (1 <= trackType && trackType <= 7)
        { // positive dimensional tracking
          pos_dim_main(trackType, genType, MPType, currentSeed, startName, my_id, num_processes, headnode);
        }
        else if (trackType == -4 || trackType == -3)
        { // function evaluation
          function_eval_main(trackType == -3, MPType, currentSeed, startName, my_id, num_processes, headnode);
        }
        else if (trackType == -2 || trackType == -1)
        { // newton evaluation
          newton_eval_main(trackType == -1, MPType, currentSeed, startName, my_id, num_processes, headnode);
        }
      }

      free(inputName);
      free(startName);

      if (remove_temp)  // remove temporary files
        remove_temp_files();
    }
  }
  else
  { // worker process
#ifdef _HAVE_MPI
    parallel_diff_worker(my_id, num_processes, headnode);

    worker_process_main(my_id, num_processes, headnode);
#endif
  }

  // finalize MPI and exit
  finalize_mpi();

  return 0;
}


