// Copyright (C) 2013 Daniel J. Bates, Jonathan D. Hauenstein, Andrew J. Sommese, and Charles W. Wampler

#include "bertini.h"
#include "ppParse.h"

void setupParameterHomotopy(FILE *OUT, FILE* IN, int paramHom, int readInOld, int maxPrec, preProcessArray *ppArray, variablegroupArray *vargpArray)
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES: setup the parameter homotopy in OUT from IN            *
\***************************************************************/
{
  int i, lineNumber = 1, count = 0, numParams = ppArray->types[PARAMETERTYPE].numType, currPrec = mpf_get_default_prec();
  char ch;
  FILE *PARAM = NULL, *PARAM2 = NULL;

  // initialize MP
  initMP(maxPrec);

  // verify parameters exist
  if (numParams < 1)
  { // error!
    printf("ERROR: There must be at least one parameter!\n");
    bexit(ERROR_INPUT_SYNTAX);
  }

  if (paramHom == 1)
  { // ab initio parameter homotopy 
    if (!readInOld)
    { // setup parameter values 
      comp_d rand_d;
      comp_mp rand_mp;
      mpq_t rand_rat[2];

      // initialize memory
      init_mp(rand_mp);
      init_rat(rand_rat);

      // generate random parameter values - put in 'start_parameters'
      PARAM = fopen("start_parameters", "w+");
      fprintf(PARAM, "%d\n\n", numParams);
      // print values to PARAM
      for (i = 0; i < numParams; i++)
      { // generate a random number
        get_comp_rand_rat(rand_d, rand_mp, rand_rat, maxPrec, maxPrec, 0, 0);
        // print number to file
        print_mp(PARAM, 0, rand_mp);
        fprintf(PARAM, "\n");
      }
      fprintf(PARAM, "\n");
      rewind(PARAM);

      // clear memory
      clear_mp(rand_mp);
      clear_rat(rand_rat);
    }
    else
    { // open parameter values
      PARAM = fopen("start_parameters", "r");
      if (PARAM == NULL)
      { // error
        printf("\nERROR: 'start_parameters' does not exist!\n");
        bexit(ERROR_FILE_NOT_EXIST);
      }
    }

    // move past the number of parameters
    fscanf(PARAM, "%d\n\n", &i);
    if (i != numParams)
    { // error
      printf("ERROR: The number of parameters in 'start_parameters' is incorrect!\n");
      bexit(ERROR_INVALID_SIZE);
    }

    // setup input system
    while (count < numParams)
    { // determine if this line defines parameters
      if (ppArray->types[PARAMETERTYPE].lineNumber[count] == lineNumber)
      { // scan in 'parameter ' & replace with 'constant ' in OUT
        fscanf(IN, "parameter ");
        fprintf(OUT, "constant ");
        while ((ch = fgetc(IN)) != '\n')
          fprintf(OUT, "%c", ch); 
        fprintf(OUT, "\n");
 
        while (count < numParams && ppArray->types[PARAMETERTYPE].lineNumber[count] == lineNumber)
        { // print the constant value of the parameter
          fprintf(OUT, "%s=", ppArray->types[PARAMETERTYPE].name[count]);
          while ((ch = fgetc(PARAM)) != ' ')
            fprintf(OUT, "%c", ch);
          fprintf(OUT, "+I*");
          while ((ch = fgetc(PARAM)) != '\n')
            fprintf(OUT, "%c", ch);
          fprintf(OUT, ";\n");

          // increment count
          count++;
        }
      }
      else 
      { // copy line to OUT
        while ((ch = fgetc(IN)) != '\n')
          fprintf(OUT, "%c", ch);
        fprintf(OUT, "\n");
      }

      // increment lineNumber
      lineNumber++;
    }
    // print the rest to OUT
    while ((ch = fgetc(IN)) != EOF)
      fprintf(OUT, "%c", ch);

    // rewind, clear ppArray and reparse
    rewind(OUT);
    clear_preProcessArray(ppArray);
    preProcessParse(ppArray, OUT, 0, 1); // no parameter homotopy, but setup for reuse

    // verify information for standard tracking
    verify_standard_homotopy(ppArray, vargpArray, 0); // zero-dim
  }
  else
  { // moving from a general parameter value to a specific parameter value
    point_mp paramStart, paramFinal;
    init_point_mp(paramStart, numParams);
    init_point_mp(paramFinal, numParams);
    paramStart->size = paramFinal->size = numParams;

    // setup homotopy using parameter values in 'start_parameters' and 'final_parameters'
    PARAM = fopen("start_parameters", "r");
    PARAM2 = fopen("final_parameters", "r");

    // verify that the files exist
    if (PARAM == NULL)
    { // 'start_parameters' does not exist
      printf("ERROR: The file 'start_parameters' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }
    else if (PARAM2 == NULL)
    { // 'final_parameters' does not exist
      printf("ERROR: The file 'final_parameters' does not exist!\n");
      bexit(ERROR_FILE_NOT_EXIST);
    }

    // read in the number of parameters in each file and verify the counts
    fscanf(PARAM, "%d", &i);
    if (i != numParams)
    {
      printf("ERROR: The number of parameters in 'start_parameters' is not correct!\n");
      bexit(ERROR_INVALID_SIZE);
    }
    fscanf(PARAM2, "%d", &i);
    if (i != numParams)
    {
      printf("ERROR: The number of parameters in 'final_parameters' is not correct!\n");
      bexit(ERROR_INVALID_SIZE);
    }

    // read in the actual parameters
    for (i = 0; i < numParams; i++)
    {
      mpf_inp_str(paramStart->coord[i].r, PARAM, 10);
      mpf_inp_str(paramStart->coord[i].i, PARAM, 10);
      scanRestOfLine(PARAM);
      mpf_inp_str(paramFinal->coord[i].r, PARAM2, 10);
      mpf_inp_str(paramFinal->coord[i].i, PARAM2, 10);
      scanRestOfLine(PARAM2);
    } 

    // setup input system
    fprintf(OUT, "pathvariable pathVarParam1234;\n");
    while (count < numParams)
    { // determine if this line defines parameters
      if (ppArray->types[PARAMETERTYPE].lineNumber[count] == lineNumber)
      { // copy parameter line to OUT
        while ((ch = fgetc(IN)) != '\n')
          fprintf(OUT, "%c", ch);
        fprintf(OUT, "\n");

        while (count < numParams && ppArray->types[PARAMETERTYPE].lineNumber[count] == lineNumber)
        { // print the parameter homotopy
          fprintf(OUT, "%s=(", ppArray->types[PARAMETERTYPE].name[count]);
          mpf_out_str(OUT, 10, 0, paramFinal->coord[count].r);
          fprintf(OUT, "+I*");
          mpf_out_str(OUT, 10, 0, paramFinal->coord[count].i);
          fprintf(OUT, ")*(1-pathVarParam1234) + pathVarParam1234*(");
          mpf_out_str(OUT, 10, 0, paramStart->coord[count].r);
          fprintf(OUT, "+I*");
          mpf_out_str(OUT, 10, 0, paramStart->coord[count].i);
          fprintf(OUT, ");\n");

          // increment count
          count++;
        }
      }
      else
      { // copy line to OUT
        while ((ch = fgetc(IN)) != '\n')
          fprintf(OUT, "%c", ch);
        fprintf(OUT, "\n");
      }

      // increment lineNumber
      lineNumber++;
    }
    // print the rest to OUT
    while ((ch = fgetc(IN)) != EOF)
      fprintf(OUT, "%c", ch);

    // rewind, clear ppArray and reparse
    rewind(OUT);
    clear_preProcessArray(ppArray);
    preProcessParse(ppArray, OUT, 0, 1); // no parameter homotopy but setup for reusing

    // close files
    fclose(PARAM);
    fclose(PARAM2);

    // clear memory
    clear_point_mp(paramStart);
    clear_point_mp(paramFinal);
  }

  // set precision back 
  initMP(currPrec);

  return;
}

