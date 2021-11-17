#include "symbolics/isosingular.hpp"



int isosingular_deflation(int *num_deflations, int **deflation_sequence,
			  BertiniRealConfig & program_options,  
						  boost::filesystem::path inputFile,
						  boost::filesystem::path witness_point_filename,
						  boost::filesystem::path output_name,
						  int max_deflations)
{
	int ii, success = 0, nullSpace = 0, *declarations = NULL;
	char ch, *strStabilizationTest = NULL;
	FILE *IN = NULL, *OUT = NULL;
	
	
	// remove previous files.
	
	
	//open the input file.
	partition_parse(&declarations, inputFile, "func_input_real" , "config_real" ,0); // the 0 means not self conjugate mode
																					 //why is this necessary?
	
	
	// setup input file to test for stabilization
	stabilization_input_file("input_stabilization_test", "func_input_real", "config_real");
	
	
	std::cout << "Performing a stabilization test\n";
	
	remove("isosingular_summary");
	

	
	// you cannot call stabilization test in parallel
	int blabla;
	parse_input_file("input_stabilization_test", &blabla);
	// seed, mptype, startName, my_id,  num_processes, headnode
	witnessGeneration(0, blabla, const_cast<char*>(witness_point_filename.c_str()), 0,1,0);
	initMP(mpf_get_default_prec());
	
	// read in the file
	
	IN = safe_fopen_read("isosingular_summary");
	fscanf(IN, "%d%d", &nullSpace, &success);     // set the success flag
	fclose(IN);
	// setup the first entry in the deflation sequence
	*num_deflations = 0;
	*deflation_sequence = (int *)br_malloc((*num_deflations + 1) * sizeof(int));
	(*deflation_sequence)[*num_deflations] = nullSpace;
	
	// loop until successful or run out of iterations
	while (*num_deflations < max_deflations && !success)
	{ // create input file for deflation
		isosingular_deflation_iteration(declarations, "func_input_real", program_options, nullSpace, *num_deflations + 1);
		
		// setup input file to test for stabilization
		stabilization_input_file("input_stabilization_test", "func_input_real", "config_real");
		
		// you cannot call stabilization test in parallel
		int blabla;
		parse_input_file("input_stabilization_test", &blabla);
					// seed, mptype, startName, my_id,  num_processes, headnode
		witnessGeneration(2565279, blabla, const_cast<char*>(witness_point_filename.c_str()), 0,1,0);
		initMP(mpf_get_default_prec());
		
		

		
		
		// read in the file
		IN = safe_fopen_read("isosingular_summary");
		fscanf(IN, "%d%d", &nullSpace, &success); // get the success indicator
		
		// setup the next entry in the deflation sequence
		(*num_deflations)++;
		*deflation_sequence = (int *)br_realloc(*deflation_sequence, (*num_deflations + 1) * sizeof(int));
		(*deflation_sequence)[*num_deflations] = nullSpace;
		
		if ((*deflation_sequence)[*num_deflations] > (*deflation_sequence)[*num_deflations - 1])
		{
			printf("ERROR: The deflation sequence must be a nonincreasing sequence of nonnegative integers!\n\n");
			success = 0;
			break;
		}
	}
	
	if (success)
	{ // print deflation sequence
		printf("\nIsosingular deflation was successful!\n\n");
		printf("Number of deflations: %d\n", *num_deflations);
		printf("Deflation sequence: ");
		for (ii = 0; ii < *num_deflations; ii++)
			printf("%d, ", (*deflation_sequence)[ii]);
		printf("%d, ...\n\n", (*deflation_sequence)[*num_deflations]);
		
		
		// create deflated system
		
		OUT = safe_fopen_write(output_name);
		fprintf(OUT, "CONFIG\n");
		IN = safe_fopen_read("config_real");
		while ((ch = fgetc(IN)) != EOF)
			fprintf(OUT, "%c", ch);
		fclose(IN);
		fprintf(OUT, "END;\nINPUT\n");
		IN = safe_fopen_read("func_input_real");
		while ((ch = fgetc(IN)) != EOF)
			fprintf(OUT, "%c", ch);
		fclose(IN);
		fprintf(OUT, "END;\n\n");
		fclose(OUT);
		
		std::cout << "Deflated system printed to '" << output_name.string() << "'." << std::endl << std::endl;
		
	}
	else if (*num_deflations >= max_deflations)
	{
		printf("ERROR: The maximum number of deflations (%d) was exceeded without stabilization!\n\n", max_deflations);
	}
	// delete temporary files
	remove("func_input_real");
	remove("config_real");
	
	// clear memory
	free(strStabilizationTest);
	free(declarations);
	
	return success;
}










void isosingular_deflation_iteration(int *declarations,
									 boost::filesystem::path inputOutputName,
				      BertiniRealConfig & program_options, 
									 int nullSpaceDim,
									 int deflation_number)
{
	int ii, numVars = 0, numFuncs = 0, numConstants = 0, minorSize = 0;
	int *lineVars = NULL, *lineFuncs = NULL, *lineConstants = NULL, *degrees = NULL;
	char ch, *str = NULL, **vars = NULL, **funcs = NULL, **consts = NULL;
	FILE *IN = NULL, *OUT = NULL;
	
	// test for existence of input file
	IN = safe_fopen_read(inputOutputName.c_str());
	fclose(IN);
	
	// move the file & open it
	rename(inputOutputName.c_str(), "deflation_input_file");
	IN = safe_fopen_read("deflation_input_file");
	// setup variables
	if (declarations[0] > 0)
	{ // using variable_group
		parse_names(&numVars, &vars, &lineVars, IN,const_cast< char *>("variable_group"), declarations[0]);
	}
	else
	{ // using hom_variable_group
		parse_names(&numVars, &vars, &lineVars, IN,const_cast< char *>("hom_variable_group"), declarations[1]);
	}
	
	// setup constants
	rewind(IN);
	parse_names(&numConstants, &consts, &lineConstants, IN,const_cast< char *>("constant"), declarations[8]);
	
	// setup functions
	rewind(IN);
	parse_names(&numFuncs, &funcs, &lineFuncs, IN, const_cast< char *>("function"), declarations[9]);
	rewind(IN);


	// read in the degrees
	degrees = (int *)br_malloc(numFuncs * sizeof(int));
	OUT = safe_fopen_read("deg.out");
	for (ii = 0; ii < numFuncs; ii++)
		fscanf(OUT, "%d", &degrees[ii]);
	fclose(OUT);


	printf("\nPerforming an isosingular deflation\n");

	switch (program_options.symbolic_engine())
	{
		case SymEngine::Matlab:
		{
			// setup Matlab script
			OUT = safe_fopen_write("matlab_deflate.m");
			minorSize = numVars - declarations[1] - nullSpaceDim + 1;
			createMatlabDeflation(OUT, numVars, vars, lineVars, numConstants, consts, lineConstants, numFuncs, funcs, lineFuncs, IN, minorSize, degrees, deflation_number);

			fclose(OUT);

			// run Matlab script
			std::stringstream converter;
			converter << program_options.matlab_command() << "matlab_deflate";
			system(converter.str().c_str());
			converter.clear(); converter.str("");
			break;
		}

    	case SymEngine::Python:
		{
			// Create Python Script
			OUT = safe_fopen_write("python_deflate.py");	    
			minorSize = numVars - declarations[1] - nullSpaceDim + 1;
			createPythonDeflation(OUT, numVars, vars, lineVars, numConstants, consts, lineConstants, numFuncs, funcs, lineFuncs, IN, minorSize, degrees, deflation_number);
			fclose(OUT);

			//run Python script
			std::cout << "inkoving `python python_deflate.py`\n";
			system("python python_deflate.py");
			break;
		}
	} // switch
	
	WaitOnGeneratedFile("deflation_polynomials_declaration");
	WaitOnGeneratedFile("deflation_polynomials");
	
	
	
	// setup new file
	OUT = safe_fopen_write(inputOutputName.c_str());
	rewind(IN);
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	
	
	IN = safe_fopen_read("deflation_polynomials_declaration");
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	
	
	
	
	IN = safe_fopen_read("deflation_polynomials");
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	fclose(OUT);
	
	// increment the number of function declarations
	declarations[9]++;
	
	// clear memory
	free(str);
	free(degrees);
	for (ii = 0; ii < numVars; ii++)
		free(vars[ii]);
	free(vars);
	free(lineVars);
	for (ii = 0; ii < numConstants; ii++)
		free(consts[ii]);
	free(consts);
	free(lineConstants);
	for (ii = 0; ii < numFuncs; ii++)
		free(funcs[ii]);
	free(funcs);
	free(lineFuncs);
	
	return;
}

void createMatlabDeflation(FILE *OUT, int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs, FILE *IN, int minorSize, int *degrees, int deflation_number)
/***************************************************************\
 * USAGE: setup a Matlab script to perform the deflation         *
 * ARGUMENTS: input file, declaration name, and number of lines  *
 * RETURN VALUES: Matlab script                                  *
 * NOTES:                                                        *
 \***************************************************************/
{
	int ii, lineNumber = 1, cont = 1, declares = 0, strLength = 0, strSize = 1;
	char ch;
	char *str = (char *)br_malloc(strSize * sizeof(char));
	
	// setup Bertini constants in Matlab
	fprintf(OUT, "syms I Pi;\nI=1i;\n");
	
	// setup variables
	fprintf(OUT, "syms");
	for (ii = 0; ii < numVars; ii++)
		fprintf(OUT, " %s", vars[ii]);
	fprintf(OUT, ";\nX = [");
	for (ii = 0; ii < numVars; ii++)
		fprintf(OUT, " %s", vars[ii]);
	fprintf(OUT, "];\n");
	
	// setup constants
	if (numConstants > 0)
	{
		fprintf(OUT, "syms");
		for (ii = 0; ii < numConstants; ii++)
			fprintf(OUT, " %s", consts[ii]);
		fprintf(OUT, ";\n");
	}
	
	// setup degrees
	fprintf(OUT, "deg = [");
	for (ii = 0; ii < numFuncs; ii++)
		fprintf(OUT, " %d", degrees[ii]);
	fprintf(OUT, "];\n");
	
	// copy lines which do not declare items or define constants (keep these as symbolic objects)
	while (cont)
	{ // see if this line number declares items
		declares = 0;
		for (ii = 0; ii < numVars; ii++)
			if (lineNumber == lineVars[ii])
				declares = 1;
		for (ii = 0; ii < numConstants; ii++)
			if (lineNumber == lineConstants[ii])
				declares = 1;
		for (ii = 0; ii < numFuncs; ii++)
			if (lineNumber == lineFuncs[ii])
				declares = 1;
		
		if (declares)
		{ // move past this line
			do
			{ // read in character
				ch = fgetc(IN);
			} while (ch != '\n' && ch != EOF);
		}
		else
		{ // check to see if this defines a constant - line must be of the form '[NAME]=[EXPRESSION];' OR EOF
			ch = fgetc(IN);
			if (ch != EOF)
			{ // definition line
				strLength = 0;
				do
				{ // add to str
					if (strLength + 1 == strSize)
					{ // increase strSize
						strSize *= 2;
						str = (char *)br_realloc(str, strSize * sizeof(char));
					}
					str[strLength] = ch;
					strLength++;
				} while ((ch = fgetc(IN)) != '=');
				str[strLength] = '\0';
				strLength++;
				
				// compare against constants
				declares = 0;
				for (ii = 0; ii < numConstants; ii++)
					if (strcmp(str, consts[ii]) == 0)
						declares = 1;
				
				if (declares)
				{ // move past this line
					do
					{ // read in character
						ch = fgetc(IN);
					} while (ch != '\n' && ch != EOF);
				}
				else
				{ // print line
					fprintf(OUT, "%s=", str);
					do
					{ // read in character & print it
						ch = fgetc(IN);
						if (ch != EOF)
							fprintf(OUT, "%c", ch);
					} while (ch != '\n' && ch != EOF);
				}
			}
		}
		
		// increment lineNumber
		lineNumber++;
		
		// test for EOF
		if (ch == EOF)
			cont = 0;
	}
	
	// setup functions
	fprintf(OUT, "\nF = [");
	for (ii = 0; ii < numFuncs; ii++)
		fprintf(OUT, " %s", funcs[ii]);
	fprintf(OUT, "];\n");
	
	// compute the jacobian
	fprintf(OUT, "J = jacobian(F,X);\n");
	
	// find the combinations of the rows & columns
	fprintf(OUT, "R = nchoosek(1:%d,%d);\nC = nchoosek(1:%d,%d);\n", numFuncs, minorSize, numVars, minorSize); // changed to nchoosek april 16, 2013 DAB.
	fprintf(OUT, "r = size(R,1);\nc = size(C,1);\n");
	
	// loop over rows & columns printing the nonzero minors to a file
	fprintf(OUT, "count = 0;\n");
	fprintf(OUT, "OUT = fopen('deflation_polynomials','w');\n");
	fprintf(OUT, "for j = 1:r\n");
	fprintf(OUT, "  for k = 1:c\n");
	fprintf(OUT, "    A = simplify(det(J(R(j,:),C(k,:))));\n");
	fprintf(OUT, "    if A ~= 0\n");
	fprintf(OUT, "      count = count + 1;\n");
	fprintf(OUT, "		curr_eqn = char(2*A/prod(deg(R(j,:))));%%\n");
	fprintf(OUT, "		i_locations = regexp(curr_eqn,'[\\W\\s]i[\\W\\s]');\n");
	fprintf(OUT, "		curr_eqn(i_locations+1) = 'I';\n");
	fprintf(OUT, "      fprintf(OUT, 'f_%d_%cd = %cs;%cn', count, curr_eqn);\n", deflation_number, '%', '%', '\\');
	fprintf(OUT, "    end;\n");
	fprintf(OUT, "  end;\n");
	fprintf(OUT, "end;\n");
	fprintf(OUT, "fclose(OUT);\n");
	
	// print number of minors to a file
	fprintf(OUT, "OUT = fopen('deflation_polynomials_declaration','w');\n");
	fprintf(OUT, "fprintf(OUT, 'function ');\n");
	fprintf(OUT, "for j = 1:count\n");
	fprintf(OUT, "  fprintf(OUT, 'f_%d_%cd', j);\n", deflation_number, '%');
	fprintf(OUT, "  if j == count");
	fprintf(OUT, "    fprintf(OUT, ';%cn');\n", '\\');
	fprintf(OUT, "  else");
	fprintf(OUT, "    fprintf(OUT, ',');\n");
	fprintf(OUT, "  end;\n");
	fprintf(OUT, "end;\n");
	fprintf(OUT, "\nexit\n");
	// clear memory
	free(str);
	
	return;
}


void createPythonDeflation(FILE *OUT, int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs, FILE *IN, int minorSize, int *degrees, int deflation_number)
/***************************************************************
 * USAGE: setup a Python script to perform the deflation         *
 * ARGUMENTS: input file, declaration name, and number of lines  *
 * RETURN VALUES: Python script                                  *
 * NOTES:                                                        *
 ***************************************************************/
{
	int ii, lineNumber = 1, cont = 1, declares = 0, strLength = 0, strSize = 1;
	char ch;
	char *str = (char *)br_malloc(strSize * sizeof(char));

	//setup Python libraries
	fprintf(OUT, "import numpy as np\nfrom numpy import array\n");
	fprintf(OUT, "import math\nimport itertools\nfrom itertools import combinations\n");
	fprintf(OUT, "import sympy as sp\nfrom sympy import *\n");
	fprintf(OUT, "import scipy as sci\nimport pickle\n\n");

	// setup nchoosek function
	fprintf(OUT,"def my_nchoosek( r, n ):\n");
	fprintf(OUT,"\tcnt = 0\n");
	fprintf(OUT,"\tfor i in combinations(range(1,r+1),n):\n");
	fprintf(OUT,"\t\tif cnt == 0:\n");
	fprintf(OUT,"\t\t\tR = Matrix([i])\n");
	fprintf(OUT,"\t\telse:\n");
	fprintf(OUT,"\t\t\tR_i = Matrix([i])\n");
	fprintf(OUT,"\t\t\tR = R.row_insert(cnt,R_i)\n");
	fprintf(OUT,"\t\tcnt = cnt + 1\n");
	fprintf(OUT,"\treturn R;\n\n");

	// setup Bertini constants in Python
	fprintf(OUT, "from numpy import pi\nfrom sympy import sin, cos, Matrix, simplify, var\n");
	fprintf(OUT, "from sympy.abc import I, pi\nI=1j\n");

	// setup variables
	fprintf(OUT, "\n# variables\n");
	fprintf(OUT, "var('");
	for (ii = 0; ii < numVars; ii++)
	{
		fprintf(OUT, "%s", vars[ii]);
		if (ii < numVars-1)
			fprintf(OUT,", ");
	}
  
	fprintf(OUT, "')\nX = Matrix([");
	for (ii = 0; ii < numVars; ii++)
	{
		fprintf(OUT, "%s", vars[ii]);
		if (ii < numVars - 1)
			fprintf(OUT, ", ");
	}
	fprintf(OUT, "])\n");
	
	// setup constants
	if (numConstants > 0)
	{
		fprintf(OUT, "\n# Constants \n");
		fprintf(OUT, "var('");
		for (ii = 0; ii < numConstants; ii++)
		{
			fprintf(OUT, "%s", consts[ii]);
			if (ii < numConstants-1)
				fprintf(OUT, ", ");
		}
		fprintf(OUT, "')\n");
	}

	// setup degrees
	fprintf(OUT, "\n# Degrees of the equations\n");
	fprintf(OUT, "deg = np.array([");
	for (ii = 0; ii < numFuncs; ii++)
	{
		fprintf(OUT, "%d", degrees[ii]);
		if (ii < numFuncs - 1)
			fprintf(OUT, ", ");
	}
	fprintf(OUT, "])\n");

	// copy lines which do not declare items or define constants (keep these as symbolic objects)
	while (cont)
	{ // see if this line number declares items
		declares = 0;
		for (ii = 0; ii < numVars; ii++)
		{
			if (lineNumber == lineVars[ii])
				declares = 1;
		}

		for (ii = 0; ii < numConstants; ii++)
		{
			if (lineNumber == lineConstants[ii])
				declares = 1;
		}
		for (ii = 0; ii < numFuncs; ii++)
		{
			if (lineNumber == lineFuncs[ii])
				declares = 1;
		}

		if (declares)
		{ // move past this line
			do
			{ // read in character
				ch = fgetc(IN);
			} 
			while (ch != '\n' && ch != EOF);
		}
		else
		{ // check to see if this defines a constant - line must be of the form '[NAME]=[EXPRESSION];' OR EOF
			ch = fgetc(IN);
			if (ch != EOF)
			{ // definition line
				strLength = 0;
				do
				{ // add to str
					if (strLength + 1 == strSize)
					{ // increase strSize
					strSize *= 2;
					str = (char *)br_realloc(str, strSize * sizeof(char));
					}
					str[strLength] = ch;
					strLength++;
				} 
				while ((ch = fgetc(IN)) != '=');

				str[strLength] = '\0';
				strLength++;

				// compare against constants
				declares = 0;
				for (ii = 0; ii < numConstants; ii++)
				{
					if (strcmp(str, consts[ii]) == 0)
						declares = 1;
				}

				if (declares)
				{ // move past this line
					do
					{ // read in character
						ch = fgetc(IN);
					} 
					while (ch != '\n' && ch != EOF);
				}
				else
				{ // print line
					fprintf(OUT, "%s=Matrix([", str);
					do
					{ // read in character & print it
						ch = fgetc(IN);
						if (ch != EOF)
						{
						if (ch == '^')
						{
						fprintf(OUT, "**");
						}
						else if (ch == ';')
						{
						fprintf(OUT, "])");
						}
						else
						{
						fprintf(OUT, "%c", ch);
						}
						}
					} 
					while (ch != '\n' && ch != EOF);
				} // else
			} // if
		} // else

		// increment lineNumber
		lineNumber++;

		// test for EOF
		if (ch == EOF)
			cont = 0;
	}  // while

	// setup functions
	fprintf(OUT, "\n# Matrix of the functions\n");
	fprintf(OUT, "\nF = sp.Matrix(");
	for (ii = 0; ii < numFuncs; ii++)
	{
		if (ii == 0)
			fprintf(OUT, "%s)\n", funcs[ii]);
		else
			fprintf(OUT, "F = F.col_insert(%i,%s)\n", ii,funcs[ii]);
	}

	fprintf(OUT, "\n");
	
	// compute the jacobian
	fprintf(OUT, "J = F.jacobian(X)\n");
	
	// find the combinations of the rows & columns
	fprintf(OUT, "\n# getting combinations of the outputs\n");
	fprintf(OUT, "R = my_nchoosek(%d,%d);\nC = my_nchoosek(%d,%d)\n", numFuncs, minorSize, numVars, minorSize); // changed to nchoosek april 16, 2013 DAB.
	fprintf(OUT, "\n# getting the sizes of R and C\n");
	fprintf(OUT, "r_t = Matrix([R.shape])\nc_t = Matrix([C.shape])\n");
	fprintf(OUT, "\n# Turning the sizes into an array\n");
	fprintf(OUT, "rv = np.array([r_t])\ncv = np.array([c_t])\n");
	fprintf(OUT, "\n# Getting the first element in the array of sizes\n");
	fprintf(OUT, "r = rv.item(0)\nc = cv.item(0)\n");
	
	// loop over rows & columns printing the nonzero minors to a file
	fprintf(OUT, "count = 0;\n");
	fprintf(OUT, "\n# Opening the first output file\n");
	fprintf(OUT, "fo = open(\"deflation_polynomials\",\"wb\")\n");
	fprintf(OUT, "for j in range(1,r+1):\n");
		fprintf(OUT, "\tfor k in range(1,c+1):\n");
			fprintf(OUT, "\t\ttemp_R = np.array(R.row(j-1))\n");
			fprintf(OUT, "\t\ttemp_C = np.array(C.row(k-1))\n");
	fprintf(OUT, "\n# Calculating whether temp_R is a 1x1 Matrix or not\n");
			fprintf(OUT, "\t\ttR_s = np.array([temp_R.shape])\n");
	fprintf(OUT, "\n# Choosing the rows in the determinant\n");
			fprintf(OUT, "\t\tif tR_s.item(0)*tR_s.item(1) > 1:\n");
				fprintf(OUT, "\t\t\tmyJ_r = Matrix([J.row(temp_R.item(0)-1),J.row(temp_R.item(1)-1)])\n");
			fprintf(OUT, "\t\telse:\n");
				fprintf(OUT, "\t\t\tmyJ_r = Matrix([J.row(temp_R.item(0)-1)])\n");
	fprintf(OUT, "\n# Calculating whether temp_C is singular or not\n");
			fprintf(OUT, "\t\ttC_s = np.array([temp_C.shape])\n");
	fprintf(OUT, "\n# Choosing the columns of the determinant\n");
			fprintf(OUT, "\t\tif tC_s.item(0)*tC_s.item(1) > 1:\n");
				fprintf(OUT, "\t\t\tmyA_1a = np.array(myJ_r.col(temp_C.item(0)-1))\n");
				fprintf(OUT, "\t\t\tmyA_1b = np.array(myJ_r.col(temp_C.item(1)-1))\n");
				fprintf(OUT, "\t\t\tmyA_1 = np.array([[myA_1a.item(0),myA_1b.item(0)],[myA_1a.item(1),myA_1b.item(1)]])\n");
			fprintf(OUT, "\t\telse:\n");
				fprintf(OUT, "\t\t\tmyA_1 = np.array([myJ_r.col(temp_C.item(0)-1)])\n");
			fprintf(OUT, "\t\tmyA = Matrix(myA_1)\n");
			fprintf(OUT, "\t\tmydet = myA.det()\n");
			fprintf(OUT, "\t\tA = simplify(mydet)\n");
			fprintf(OUT, "\t\tif A != 0:\n");
				fprintf(OUT, "\t\t\tcount = count + 1\n");
				fprintf(OUT, "\t\t\ttx = temp_R.item(0)-1\n");
				fprintf(OUT, "\t\t\tif temp_R.size > 1:\n");
					fprintf(OUT, "\t\t\t\tty = temp_R.item(1)-1\n");
				fprintf(OUT, "\t\t\telse:\n");
					fprintf(OUT, "\t\t\t\tty = 0\n");
				fprintf(OUT, "\t\t\tmydeg = np.array([deg.item(tx),deg.item(ty)])\n");
				fprintf(OUT, "\t\t\tmyprod = np.multiply(mydeg.item(0),mydeg.item(1))\n");
				fprintf(OUT, "\t\t\tmyfunc = 2*A/myprod\n");
				fprintf(OUT, "\t\t\tmyfunc2 = simplify(myfunc)\n");
				fprintf(OUT, "\t\t\tcurr_eqn = str(myfunc)\n");
				fprintf(OUT, "\t\t\ttf_i = curr_eqn.find(\"1j\")\n");
				fprintf(OUT, "\t\t\tif tf_i != -1:\n");
					fprintf(OUT, "\t\t\t\tcurr_eqn = str.replace(\"1j\",\"I\")\n");
				fprintf(OUT, "\t\t\tfunc_ida = \"f_%d_%%i = \" %% count\n", deflation_number);
				fprintf(OUT, "\t\t\tfunc_idb = str(curr_eqn)\n");
				fprintf(OUT, "\t\t\tfunc_id = func_ida+func_idb+\";\\n\"\n");
				fprintf(OUT, "\t\t\tfo.write(func_id.replace(\"**\",\"^\"))\n");
	fprintf(OUT, "fo.close()\n\n");
	
	fprintf(OUT, "\n# Opening the second output file\n");
	fprintf(OUT, "fo = open(\"deflation_polynomials_declaration\",\"wb\")\n");
	fprintf(OUT, "mystr4 = \"\"\n");
	fprintf(OUT, "for j in range(1,count+1):\n");
		fprintf(OUT, "\tmystr2a = \"f_%d_%%i\" %% (j)\n", deflation_number);
		fprintf(OUT, "\tif j == count:\n");
			fprintf(OUT, "\t\tmystr2b = \";\\n\"\n");
		fprintf(OUT, "\telse:\n");
			fprintf(OUT, "\t\tmystr2b = \", \"\n");
		fprintf(OUT, "\tmystr3 = \"%%s%%s\" %% (mystr2a,mystr2b)\n");
		fprintf(OUT, "\tmystr4 = \"%%s%%s\" %% (mystr4,mystr3)\n");


	fprintf(OUT, "\n\nfo.write(\"function \")\n");
	fprintf(OUT, "fo.write(mystr4)\n");
	fprintf(OUT, "fo.write(\"\\n\")\n");
	fprintf(OUT, "fo.close()\n\n");

	// clear memory
	free(str);
	
	return;
}




void stabilization_input_file(boost::filesystem::path outputFile,
							  boost::filesystem::path funcInput,
							  boost::filesystem::path configInput)
{
	char ch;
	FILE *OUT = safe_fopen_write(outputFile.c_str());
	
	// setup configurations in OUT
	fprintf(OUT, "CONFIG\n");
	
	FILE *IN = safe_fopen_read(configInput.c_str());
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	fprintf(OUT, "TrackType: 6;\nMultOneOnly: 1;\nDeleteTempFiles: 0;\nTargetTolMultiplier: 1e3; %% this line added so that computed solutions fed in from BR will not falsely cause to fail -- DAB\nEND;\nINPUT\n");
	
	// setup system in OUT
	IN = safe_fopen_read(funcInput.c_str());
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	fprintf(OUT, "END;\n");
	fclose(OUT);
	
	return;
}





void check_declarations(int *declarations){
	
	// error check based on declarations
	if ( (declarations[0] + declarations[1]) > 1)
	{ // multiple variable groups
		printf("\n\n%d variable_groups, %d hom_variable_groups\nERROR: Please use either one 'variable_group' or one 'hom_variable_group'!\n\n",declarations[0],declarations[1]);
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[2] > 0)
	{ // 'variable'
		printf("\n\nERROR: The system needs to be defined using either 'variable_group' or 'hom_variable_group'.\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[3] > 0)
	{ // 'pathvariable'
		printf("\n\nERROR: The system should not depend on pathvariables!\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[4] > 0)
	{ // 'parameter'
		printf("\n\nERROR: The system should not depend on parameters!\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[5] > 0)
	{ // 'definedSubfunction'
		printf("\n\nERROR: The system should not depend on defined subfunctions!\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[6] > 0)
	{ // 'random'
		printf("\n\nERROR: The system should not involve 'random' numbers!\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	else if (declarations[7] > 0)
	{ // 'random_real'
		printf("\n\nERROR: The system should not involve 'random_real' numbers!\n\n");
		bexit(ERROR_CONFIGURATION);
	}
	
}

