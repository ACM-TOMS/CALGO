// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "cascade.h"

void cascadeOutputChart(cascade_t *CD, FILE *fp, int MPType)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int codim_index, total_paths = 0, num_codim = CD->num_codim;

  fprintf(fp, "\n\n****************** Cascade Summary *******************\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|codim|   paths   |witness superset| nonsingular | singular |nonsolutions| inf endpoints | other endpoints | bad endpoints\n");
  fprintf(fp, "--------------------------------------------------------------------------------------------------------------------------\n");

  for (codim_index = 0; codim_index < num_codim; codim_index++)
  {
    // add to the total paths tracked
    total_paths += CD->codim[codim_index].num_paths;

    fprintf(fp, "| %-4d|   %-8d|   %-13d|  %-11d|  %-8d|  %-10d|   %-12d|  %-15d|  %d\n", CD->codim[codim_index].codim, CD->codim[codim_index].num_paths, CD->codim[codim_index].num_superset, CD->codim[codim_index].num_nonsing, CD->codim[codim_index].num_sing, CD->codim[codim_index].num_nonsolns, CD->codim[codim_index].num_inf, CD->codim[codim_index].num_other, CD->codim[codim_index].num_bad);
  }
  fprintf(fp, "--------------------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|total|   %d\n\n", total_paths);

  fprintf(fp, "****************************************************\n\n");

  return;
}

