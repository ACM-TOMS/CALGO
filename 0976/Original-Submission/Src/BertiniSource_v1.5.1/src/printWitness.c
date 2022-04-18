// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "cascade.h"
#include "pos_dim.h"
#include "parallel.h"

// Print witness set
void printWitnessMenu(witness_t *W, int MPType, int max_prec);

void printWitnessMain(unsigned int currentSeed, int MPType, int my_id, int num_processes, int headnode)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the linear system and witness points to files    *
\***************************************************************/
{
  int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, constructWitnessSet = 0, supersetOnly = 0, paramHom = 0;
  double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
  tracker_config_t T;
  witness_t witnessSet;

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

  // quit if we are using MPI
  if (num_processes > 1)
  { // exit since parallel sampling is not implemented
#ifdef _HAVE_MPI
    printf("ERROR: Parallel printing is not implemented. Please use sequential version!\n");
    bexit(ERROR_OTHER);
#endif
  }

  // setup witnessSet
  setupWitnessDataFromFile("witness_data", "witness_data_old", "preproc_data", "deg.out", &witnessSet, &T, 1);

  // error checking
  if (witnessSet.Prog->numPathVars > 0)
  { // path variable present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }
  if (witnessSet.Prog->numPars > 0)
  { // parameter present
    rename("witness_data_old", "witness_data");

    printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
    bexit(ERROR_INPUT_SYSTEM);
  }

  // setup the rest of T, if needed
  if (T.MPType == 1)
  { // initialize latest_newton_residual_mp
    mpf_init(T.latest_newton_residual_mp);
  }
  else if (T.MPType == 2)
  { // initialize latest_newton_residual_mp
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);

    // setup eps, Phi & Psi
    T.AMP_eps = (double) T.numVars * T.numVars;
    T.AMP_Phi = T.AMP_bound_on_degree * (T.AMP_bound_on_degree - 1) * T.AMP_bound_on_abs_vals_of_coeffs;
    T.AMP_Psi = T.AMP_bound_on_degree * T.AMP_bound_on_abs_vals_of_coeffs;
  }

  // create output files now that everything is setup
  numIrredDecompOutput(&witnessSet, &T, 4, 2, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, constructWitnessSet, supersetOnly, paramHom); // trackType == 4

  // go to the menu
  printWitnessMenu(&witnessSet, T.MPType, T.AMP_max_prec);

  // clear witnessSet
  witness_clear(&witnessSet, T.MPType);

  // clear T
  tracker_config_clear(&T);

  // clear MP
  clearMP();

  return;
}

void printWitnessMenu(witness_t *W, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: menu for printing witness to file                      *
\***************************************************************/
{
  int i, j, codim_index, min_deg, max_deg, count, rV, dim_number, component_number, selection_made, size_of_string = 255;
  int *degrees = NULL, *dim = (int *)bmalloc(W->num_codim * sizeof(int)), *codim_good = (int *)bmalloc(W->num_codim * sizeof(int));
  char ch, *tempStr = NULL, *pointsFile = NULL, *linearFile = NULL;

  // find the dimension and good dimensions, and make sure one exists
  rV = 0;
  for (codim_index = 0; codim_index < W->num_codim; codim_index++)
  { // see if there are classified components for this codim
    codim_good[codim_index] = (W->codim[codim_index].num_components > 0 ? 1 : 0);

    if (codim_good[codim_index])
      rV = 1;

    // determine the dimension of this codim
    dim[codim_index] = W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp;
  }

  if (!rV)
  { // there are no classified components!!
    printf("\nThere are no classified components to print!\n\n");
    free(dim);
    free(codim_good);
    return;
  }

  // so we have atleast one classified component
  do
  { // initialize
    selection_made = 0;

    // print title
    printf("\n\n*************** Components to Print ****************\n\n");

    // display a catalog of the available components in each codim
    for (codim_index = 0; codim_index < W->num_codim; codim_index++)
    { // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      // find the minimum and maximum degree
      max_deg = 0;
      min_deg = W->codim[codim_index].num_set + 1;
      for (i = 0; i < W->codim[codim_index].num_components; i++)
      { // find maximum degree
        if (max_deg < degrees[i])
          max_deg = degrees[i];

        // find minimum degree
        if (min_deg > degrees[i])
          min_deg = degrees[i];
      }

      if (W->codim[codim_index].num_components > 0)
      {
        printf("Dimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
        if (W->codim[codim_index].num_components == 1)
          printf("\n");
        else
          printf("s\n");
        printf("-----------------------------------------------------\n");

        // display the summary
        for (i = min_deg; i <= max_deg; i++)
        { // count the number that have degree == i
          count = 0;
          for (j = 0; j < W->codim[codim_index].num_components; j++)
            if (degrees[j] == i)
              count++;

          if (count > 0)
          { // display the number
            printf("   degree %d: %d component", i, count);
            if (count == 1)
              printf("\n");
            else
              printf("s\n");
          }
        }
        printf("\n");
      }
    }

    printf("\nPlease select a dimension to print (-1 to quit): ");
    rV = scanf("%d", &dim_number);

    if (rV < 0)
    { // at EOF - so we need to quit
      dim_number = -1;
      selection_made = 1;
    }
    else
    { // we are not at EOF - flush the buffer
      do
      {
        ch = getchar();
      } while (ch != EOF && ch != '\n');

      if (rV == 0)
      { // invalid input
        printf("\nThe input was not read in correctly!\n");
        selection_made = 0;
      }
      else if (dim_number == -1)
      { // dim_number is valid
        selection_made = 1;
      }
      else
      { // verify that the dimension is one of the good ones
        selection_made = 0;
        for (j = 0; j < W->num_codim; j++)
          if (codim_good[j] && dim[j] == dim_number)
          { // selection is valid
            selection_made = 1;
            break;
          }

        if (!selection_made)
        { // dim_number is not valid
          printf("\nThe dimension %d is not valid!\n", dim_number);
        }
      }
    }
  } while (!selection_made);

  // so, either dim_number == -1 OR dim_number is one of the ones that have a classified component
  if (dim_number != -1)
  { // find the codim_index for the selected dimension
    codim_index = 0;
    for (i = 0; i < W->num_codim; i++)
    { // see if the dimension agrees
      if (dim[i] == dim_number)
      { // store the index and exit loop
        codim_index = i;
        break;
      }
    }

    do
    { // initialize
      selection_made = 0;

      // determine the degree of each component
      degrees = (int *)brealloc(degrees, W->codim[codim_index].num_components * sizeof(int));
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        degrees[i] = 0;

      for (i = 0; i < W->codim[codim_index].num_set; i++)
      { // increment degree[component_nums[i]]
        if (0 <= W->codim[codim_index].component_nums[i] && W->codim[codim_index].component_nums[i] < W->codim[codim_index].num_components)
          degrees[W->codim[codim_index].component_nums[i]]++;
      }

      printf("\nDimension %d: %d classified component", W->orig_variables - W->codim[codim_index].codim - W->PPD.num_var_gp - W->PPD.num_hom_var_gp, W->codim[codim_index].num_components);
      if (W->codim[codim_index].num_components == 1)
        printf("\n");
      else
        printf("s\n");
      printf("-----------------------------------------------------\n");

      // display the summary of components
      for (i = 0; i < W->codim[codim_index].num_components; i++)
        printf("   component %d has degree %d\n", i, degrees[i]);
      printf("\n");

      printf("\nPlease select a component to print (-1 to quit, -2 to print all): ");
      rV = scanf("%d", &component_number);

      if (rV < 0)
      { // at EOF - so we need to quit
        component_number = -1;
        selection_made = 1;
      }
      else
      { // we are not at EOF - flush the buffer
        do
        {
          ch = getchar();
        } while (ch != EOF && ch != '\n');

        if (rV == 0)
        { // invalid input
          printf("\nThe input was not read in correctly!\n");
          selection_made = 0;
        }
        else if (component_number < -2 || component_number >= W->codim[codim_index].num_components)
        { // component_number is not valid
          printf("\nThe component %d is not valid!\n", component_number);
          selection_made = 0;
        }
        else
        { // component_number is valid
          selection_made = 1;
        }
      }
    } while (!selection_made);

    // so, either component_number == -2 OR component_number == -1 OR 0 <= component_number < num_components
    if (component_number == -2 || (0 <= component_number && component_number < W->codim[codim_index].num_components))
    { // find the name of the files to use : linear & points
      tempStr = (char *)bmalloc(((int) log10(size_of_string) + 10) * sizeof(char));
      snprintf(tempStr, size_of_string + 10, "%%%ds", size_of_string);
      pointsFile = (char *)bmalloc((size_of_string + 1) * sizeof(char));
      linearFile = (char *)bmalloc((size_of_string + 1) * sizeof(char));
      for (i = 0; i <= size_of_string; i++)
      {
        pointsFile[i] = '\0';
        linearFile[i] = '\0';
      }
      do
      { // initialize
        selection_made = 0;

        printf("Enter the name of the file where to write the witness points (max of %d characters): ", size_of_string);
        rV = scanf(tempStr, pointsFile);

        if (rV < 0)
        { // at EOF - setup to be 'witness_points'
          sprintf(pointsFile, "witness_points");
          selection_made = 1;
        }
        else
        { // we are not at EOF - flush the buffer
          do
          {
            ch = getchar();
          } while (ch != EOF && ch != '\n');

          // check to see if it is a valid file name
          FILE *TEMP = fopen(pointsFile, "w");
          if (TEMP == NULL)
          { // not valid
            printf("\nThe name \"%s\" is not valid!\n\n", pointsFile);
            selection_made = 0;
          }
          else
          { // valid - so close
            fclose(TEMP);
            selection_made = 1;
          }
        }
      } while (!selection_made);

      do
      { // initialize
        selection_made = 0;

        printf("Enter the name of the file where to write the linear system (max of %d characters): ", size_of_string);
        rV = scanf(tempStr, linearFile);

        if (rV < 0)
        { // at EOF - setup to be 'linear_system'
          sprintf(linearFile, "linear_system");
          selection_made = 1;
        }
        else
        { // we are not at EOF - flush the buffer
          do
          {
            ch = getchar();
          } while (ch != EOF && ch != '\n');

          // check to see if it is a valid file name
          FILE *TEMP = fopen(linearFile, "w");
          if (TEMP == NULL)
          { // not valid
            printf("\nThe name \"%s\" is not valid!\n\n", linearFile);
            selection_made = 0;
          }
          else
          { // valid - so close
            fclose(TEMP);
            selection_made = 1;
          }
        }
      } while (!selection_made);

      // create the files
      printf("\nWriting the witness points to '%s' and the linear system to '%s'.\n\n", pointsFile, linearFile);
      printLinearSystem(W, W->codim[codim_index].codim, component_number, pointsFile, linearFile, MPType, max_prec);
 
      free(pointsFile);
      free(linearFile);
      free(tempStr);
    }
  }

  free(degrees);
  free(dim);
  free(codim_good);

  return;
}


void printLinearSystem(witness_t *W, int codim, int component_num, char *pointsName, char *linearName, int MPType, int max_prec)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: print the linear system and witness points to files    *
*  component_num == -2 : print all points for codimension       *
\***************************************************************/
{
  int i, j, size, num_vars = W->orig_variables, codim_index, num_slices, num_points;
  FILE *NAMES, *PTS, *LINEAR;

  // find the codim index
  codim_index = -1;
  for (i = 0; i < W->num_codim; i++)
    if (W->codim[i].codim == codim)
    {
      codim_index = i;
      break;
    }
  if (codim_index == -1)
  {
    printf("\n\nERROR: Codimension %d does not exist!\n", codim);
    bexit(ERROR_CONFIGURATION);
  }

  // find the number of slices
  if (MPType == 0 || MPType == 2)
    num_slices = W->codim[codim_index].B_d->rows;
  else
    num_slices = W->codim[codim_index].B_mp->rows;
  
  // find the number of points to print
  if (component_num == -2)
    num_points = W->codim[codim_index].num_set;
  else
  { // find the number of points in the component'
    num_points = 0;
    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if (W->codim[codim_index].component_nums[i] == component_num)
        num_points++;
  }

  if (num_points == 0)
  {
    printf("\n\nERROR: Invalid component number (%d).\n", component_num);
    bexit(ERROR_CONFIGURATION);
  }

  // print the linear system (if needed)
  if (num_slices > 0)
  { // table for names of variables 
    char **name_table = (char **)bmalloc(num_vars * sizeof(char *));

    // open the files
    NAMES = fopen("names.out", "r");
    LINEAR = fopen(linearName, "w");
    if (NAMES == NULL)
    {
      printf("ERROR: 'names.out' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }
    else if (LINEAR == NULL)
    {
      printf("\n\nERROR: '%s' is an invalid name!!\n\n", linearName);
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // read in the name of the variables
    for (i = 0; i < num_vars; i++)
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

    // print the variables
    if (W->PPD.num_var_gp)
      fprintf(LINEAR, "variable_group ");
    else
      fprintf(LINEAR, "hom_variable_group ");
    for (i = W->PPD.num_var_gp; i < num_vars; i++)
      fprintf(LINEAR, "%s%s", name_table[i], i+1 == num_vars ? ";\n" : ",");
  
    // print the function names
    fprintf(LINEAR, "function ");
    for (i = 0; i < num_slices; i++)
      fprintf(LINEAR, "linear%d%s", i+1, i+1 == num_slices ? ";\n" : ",");

    // print the constant names
    for (i = 0; i < num_slices; i++)
    {
      fprintf(LINEAR, "constant ");
      for (j = 0; j < num_vars; j++)
        fprintf(LINEAR, "const%d_%d%s", i+1, j+1, j+1 == num_vars ? ";\n" : ",");
    }
    fprintf(LINEAR, "\n");

    // print the constants
    if (MPType == 0)
    { // print using _d
      for (i = 0; i < num_slices; i++)
        for (j = 0; j < num_vars; j++)
          if (W->codim[codim_index].B_d->entry[i][j].i >= 0)
            fprintf(LINEAR, "const%d_%d = %.15e+%.15e*I;\n", i+1, j+1, W->codim[codim_index].B_d->entry[i][j].r,  W->codim[codim_index].B_d->entry[i][j].i);
          else
            fprintf(LINEAR, "const%d_%d = %.15e%.15e*I;\n", i+1, j+1, W->codim[codim_index].B_d->entry[i][j].r,  W->codim[codim_index].B_d->entry[i][j].i);
    }
    else if (MPType == 1)
    { // print using _mp
      for (i = 0; i < num_slices; i++)
        for (j = 0; j < num_vars; j++)
        {
          fprintf(LINEAR, "const%d_%d = ", i+1, j+1);
          mpf_out_str(LINEAR, 10, 0, W->codim[codim_index].B_mp->entry[i][j].r);
          if (mpfr_sgn(W->codim[codim_index].B_mp->entry[i][j].i) >= 0)
            fprintf(LINEAR, "+");
          mpf_out_str(LINEAR, 10, 0, W->codim[codim_index].B_mp->entry[i][j].i);
          fprintf(LINEAR, "*I;\n");
        }
    }
    else
    { // print using _rat
      comp_mp temp;
      init_mp2(temp, max_prec);
      for (i = 0; i < num_slices; i++)
        for (j = 0; j < num_vars; j++)
        {
          fprintf(LINEAR, "const%d_%d = ", i+1, j+1);
          rat_to_mp(temp, W->codim[codim_index].B_rat[i][j]);
          mpf_out_str(LINEAR, 10, 0, temp->r);
          if (mpfr_sgn(temp->i) >= 0)
            fprintf(LINEAR, "+");
          mpf_out_str(LINEAR, 10, 0, temp->i);
          fprintf(LINEAR, "*I;\n");
        }
      clear_mp(temp);
    }

    // print the linears
    fprintf(LINEAR, "\n");
    for (i = 0; i < num_slices; i++)
    {
      fprintf(LINEAR, "linear%d = ", i+1);
      for (j = 0; j < num_vars; j++)
        if (j == 0 && W->PPD.num_var_gp)
          fprintf(LINEAR, "const%d_%d%s", i+1, j+1, " + ");
        else
          fprintf(LINEAR, "const%d_%d*%s%s", i+1, j+1, name_table[j], j+1 == num_vars ? ";\n" : " + ");
    }

    // print the bottom of file    
    fprintf(LINEAR, "\nEND;\n\n");

    // close the files
    fclose(NAMES);
    fclose(LINEAR);

    // clear memory
    for (i = 0; i < num_vars; i++)
      free(name_table[i]);
    free(name_table);
  }

  // print the points  
  PTS = fopen(pointsName, "w");
  if (PTS == NULL)
  {
    printf("\n\nERROR: '%s' is an invalid name!!\n\n", pointsName);
    bexit(ERROR_FILE_NOT_EXIST);
  }

  // print the number of points and then the points
  fprintf(PTS, "%d\n\n", num_points);
  if (MPType == 0)
  { // print using _d
    vec_d dehom;
    init_vec_d(dehom, 0);

    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if (component_num == -2 || W->codim[codim_index].component_nums[i] == component_num)
      { // find the dehomogenized point
        witnessFindDehom_d(dehom, W->codim[codim_index].witnessPts_d[i].endPt, W, codim_index);
        // print dehomogenized point
        for (j = 0; j < dehom->size; j++)
        {
          print_d(PTS, 0, &dehom->coord[j]);
          fprintf(PTS, "\n");
        }
        fprintf(PTS, "\n");
      }

    clear_vec_d(dehom);
  }
  else if (MPType == 1)
  { // print using _mp
    vec_mp dehom;
    init_vec_mp(dehom, 0);

    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if (component_num == -2 || W->codim[codim_index].component_nums[i] == component_num)
      { // find the dehomogenized point
        witnessFindDehom_mp(dehom, W->codim[codim_index].witnessPts_mp[i].endPt, W, codim_index, W->curr_precision);
        // print dehomogenized point
        for (j = 0; j < dehom->size; j++)
        {
          print_mp(PTS, 0, &dehom->coord[j]);
          fprintf(PTS, "\n");
        }
      }

    clear_vec_mp(dehom);
  }
  else
  { // print
    vec_d dehom_d;
    vec_mp dehom_mp;
    init_vec_d(dehom_d, 0);
    init_vec_mp(dehom_mp, 0);

    for (i = 0; i < W->codim[codim_index].num_set; i++)
      if (component_num == -2 || W->codim[codim_index].component_nums[i] == component_num)
      { // find the dehomogenized point
        if (W->codim[codim_index].witnessPts_amp[i].curr_prec < 64)
        { // find the dehomogenized point using double precision
          witnessFindDehom_d(dehom_d, W->codim[codim_index].witnessPts_amp[i].endPt_d, W, codim_index);
          // print dehomogenized point
          for (j = 0; j < dehom_d->size; j++)
          {
            print_d(PTS, 0, &dehom_d->coord[j]);
            fprintf(PTS, "\n");
          }
          fprintf(PTS, "\n");
        }
        else
        { // find the dehomogenized point using multi precision

          // set the precision so that the calculations are correct
          initMP(W->codim[codim_index].witnessPts_amp[i].curr_prec);
          setprec_point_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].curr_prec);
          // find the dehomogenized point
          witnessFindDehom_mp(dehom_mp, W->codim[codim_index].witnessPts_amp[i].endPt_mp, W, codim_index, W->codim[codim_index].witnessPts_amp[i].curr_prec);
          // print dehomogenized point
          for (j = 0; j < dehom_mp->size; j++)
          {
            print_mp(PTS, 0, &dehom_mp->coord[j]);
            fprintf(PTS, "\n");
          }
          fprintf(PTS, "\n");
        }
      }

    clear_vec_d(dehom_d);
    clear_vec_mp(dehom_mp);
  }

  fclose(PTS);

  return;
}

