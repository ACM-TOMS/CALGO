// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "eqbyeq.h"

void eqbyeqOutputChart_d(eqData_t *EqD, FILE *fp, int infRemoved)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int stage, discarded = 0, total_paths = 0, num_stages = EqD->num_subsystems;

  fprintf(fp, "\n\n********** Witness set generation summary **********\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|subsystem|  paths  | nonsing endpoints | total discarded | sing endpoints |");
  if (infRemoved)
    fprintf(fp, " inf endpoints |");
  fprintf(fp, " higher dim'l | other bad endpoints\n");
  fprintf(fp, "-------------------------------------------------------------------------------------------------------------------------------\n");

  for (stage = 0; stage < num_stages; stage++)
  { // add to the total paths tracked
    total_paths += EqD->witnessData_d[stage].num_paths;

    // find the number discarded at this level
    discarded = EqD->witnessData_d[stage].num_sing + EqD->witnessData_d[stage].num_bad;
    if (infRemoved)
      discarded += EqD->witnessData_d[stage].num_inf;

    fprintf(fp, "|   %-6d|  %-7d| %-18d| %-16d| %-15d|", stage, EqD->witnessData_d[stage].num_paths, EqD->witnessData_d[stage].num_nonsing, discarded, EqD->witnessData_d[stage].num_sing);

    if (infRemoved)
      fprintf(fp, " %-14d|", EqD->witnessData_d[stage].num_inf);
    fprintf(fp, "  %-12d| %d", EqD->witnessData_d[stage].num_higher_dim, EqD->witnessData_d[stage].num_bad);
    if (num_stages == 1 && EqD->witnessData_d[stage].num_bad > 0)
      fprintf(fp, " ** see failure summary");

    fprintf(fp, "\n");
  }
  fprintf(fp, "-------------------------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|  total  |  %d\n\n", total_paths);
  fprintf(fp, "****************************************************\n\n");

  if (num_stages > 1)
  {
    discarded = total_paths = 0;

    fprintf(fp, "\n\n*********** Equation-by-Equation summary ***********\n\n");
    fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
    fprintf(fp, "|stage|   paths   | nonsing endpoints | total discarded | sing endpoints |");
    if (infRemoved)
      fprintf(fp, " inf endpoints |");
    fprintf(fp, " higher dim'l | other bad endpoints\n");
    fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");

    for (stage = 1; stage < num_stages; stage++)
    { // add to the total paths tracked
      total_paths += EqD->stageData_d[stage].num_paths;

      // find the number discarded at this level
      discarded = EqD->stageData_d[stage].num_sing + EqD->stageData_d[stage].num_bad;
      if (infRemoved)
        discarded += EqD->stageData_d[stage].num_inf;

      fprintf(fp, "| %-4d|   %-8d| %-18d| %-16d| %-15d|", EqD->stageData_d[stage].depth_x + EqD->stageData_d[stage].depth_y - 1, EqD->stageData_d[stage].num_paths, EqD->stageData_d[stage].num_nonsing, discarded, EqD->stageData_d[stage].num_sing);

      if (infRemoved)
        fprintf(fp, " %-14d|", EqD->stageData_d[stage].num_inf);
      fprintf(fp, "  %-12d| %d", EqD->stageData_d[stage].num_higher_dim, EqD->stageData_d[stage].num_bad);
   
      if (stage + 1 == num_stages && EqD->stageData_d[stage].num_bad > 0)
        fprintf(fp, " ** see failure summary");
      fprintf(fp, "\n");
    }
    fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "|total|   %d\n\n", total_paths);
    fprintf(fp, "****************************************************\n\n");
  }

  return;
}

//////////////////////// MP FUNCTIONS /////////////////////////

void eqbyeqOutputChart_mp(eqData_t *EqD, FILE *fp, int infRemoved)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: prints output chart                                    *
\***************************************************************/
{
  int stage, discarded = 0, total_paths = 0, num_stages = EqD->num_subsystems;

  fprintf(fp, "\n\n********** Witness set generation summary **********\n\n");
  fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
  fprintf(fp, "|subsystem|  paths  | nonsing endpoints | total discarded | sing endpoints |");
  if (infRemoved)
    fprintf(fp, " inf endpoints |");
  fprintf(fp, " higher dim'l | other bad endpoints\n");
  fprintf(fp, "-------------------------------------------------------------------------------------------------------------------------------\n");  

  for (stage = 0; stage < num_stages; stage++)
  { // add to the total paths tracked
    total_paths += EqD->witnessData_mp[stage].num_paths;

    // find the number discarded at this level
    discarded = EqD->witnessData_mp[stage].num_sing + EqD->witnessData_mp[stage].num_bad;
    if (infRemoved)
      discarded += EqD->witnessData_mp[stage].num_inf;

    fprintf(fp, "|   %-6d|  %-7d| %-18d| %-16d| %-15d|", stage, EqD->witnessData_mp[stage].num_paths, EqD->witnessData_mp[stage].num_nonsing, discarded, EqD->witnessData_mp[stage].num_sing);

    if (infRemoved)
      fprintf(fp, " %-14d|", EqD->witnessData_mp[stage].num_inf);
    fprintf(fp, "  %-12d| %d", EqD->witnessData_mp[stage].num_higher_dim, EqD->witnessData_mp[stage].num_bad);
    if (num_stages == 1 && EqD->witnessData_mp[stage].num_bad > 0)
      fprintf(fp, " ** see failure summary");

    fprintf(fp, "\n");
  }
  fprintf(fp, "-------------------------------------------------------------------------------------------------------------------------------\n");
  fprintf(fp, "|  total  |  %d\n\n", total_paths);
  fprintf(fp, "****************************************************\n\n");

  if (num_stages > 1)
  {
    discarded = total_paths = 0;

    fprintf(fp, "\n\n*********** Equation-by-Equation summary ***********\n\n");
    fprintf(fp, "NOTE: nonsingular vs singular is based on rank deficiency and identical endpoints\n\n");
    fprintf(fp, "|stage|   paths   | nonsing endpoints | total discarded | sing endpoints |");
    if (infRemoved)
      fprintf(fp, " inf endpoints |");
    fprintf(fp, " higher dim'l | other bad endpoints\n");
    fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");

    for (stage = 1; stage < num_stages; stage++)
    { // add to the total paths tracked
      total_paths += EqD->stageData_mp[stage].num_paths;

      // find the number discarded at this level
      discarded = EqD->stageData_mp[stage].num_sing + EqD->stageData_mp[stage].num_bad;
      if (infRemoved)
        discarded += EqD->stageData_mp[stage].num_inf;

      fprintf(fp, "| %-4d|   %-8d| %-18d| %-16d| %-15d|", EqD->stageData_mp[stage].depth_x + EqD->stageData_mp[stage].depth_y - 1, EqD->stageData_mp[stage].num_paths, EqD->stageData_mp[stage].num_nonsing, discarded, EqD->stageData_mp[stage].num_sing);

      if (infRemoved)
        fprintf(fp, " %-14d|", EqD->stageData_mp[stage].num_inf);
      fprintf(fp, "  %-12d| %d", EqD->stageData_mp[stage].num_higher_dim, EqD->stageData_mp[stage].num_bad);

      if (stage + 1 == num_stages && EqD->stageData_mp[stage].num_bad > 0)
        fprintf(fp, " ** see failure summary");
      fprintf(fp, "\n");
    }
    fprintf(fp, "-----------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(fp, "|total|   %d\n\n", total_paths);
    fprintf(fp, "****************************************************\n\n");
  }

  return;
}

