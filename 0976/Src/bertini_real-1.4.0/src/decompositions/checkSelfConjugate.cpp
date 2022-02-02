#include "decompositions/checkSelfConjugate.hpp"






bool checkSelfConjugate(vec_mp test_point,
						BertiniRealConfig & program_options,
						boost::filesystem::path input_file)
{
	
	
	// setup input file
	int *declarations = NULL;
	partition_parse(&declarations, input_file, "func_input_real", "config_real",0); // the 0 means not self conjugate
	free(declarations);
	
	
	//check existence of the required witness_data file.
	FILE *IN = safe_fopen_read("witness_data");
	fclose(IN);
	
	
	
	
	//we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
	membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_sc(test_point);//,fmt
	
	
	int blabla;
	parse_input_file("input_membership_test", &blabla);
	membershipMain(13423, blabla, 0, 1, 0);
	initMP(mpf_get_default_prec());
	
	
	std::vector<int> component_numbers = read_incidence_matrix();
	
	
	
	
	if (component_numbers[0]==component_numbers[1]) {
//		printf("component IS self conjugate\n");
		return true;
	}
	else
	{
//		printf("component is NOT self conjugate\n");
		return false;
	}
	
	
}



//TODO: this function has an error, in that it will only return the last component index the points lie on.
std::vector<int> read_incidence_matrix()
{
	
	std::vector<int> component_numbers;
	

	
	
	//open incidence_matrix file
	FILE *IN = safe_fopen_read("incidence_matrix");
	
	
	
	
	int num_nonempty_codims;
	fscanf(IN, "%d",&num_nonempty_codims);      // number of nonempty codimensions
	
	int total_num_components = 0; // create and initialize counter
	for (int ii = 0; ii<num_nonempty_codims; ii++) {
		int codim, num_components;
		fscanf(IN, "%d",&codim);      // codimension  (iterated for each codimension)
		fscanf(IN, "%d",&num_components);  // number of components (is whatever)
		total_num_components = total_num_components+num_components;
	}
	
	int num_pts;
	fscanf(IN, "%d", &num_pts);    // number of points
								  //	printf("reading incidence for %d pts\n",num_pts);
	component_numbers.resize(num_pts);
	
	//and then a binary matrix indicating membership on which component
	//from the appendices of the book:
	//	% Binary matrix with one row per point, columns corresponding to components.
	//	%                0 if given point is not on given component;
	//	%                1 else .
	for (int jj=0; jj<num_pts; jj++) {
		int component_number=-10, component_indicator;;
		for (int ii=0; ii<total_num_components; ii++)  // i don't this is correct if there is more than one nonempty codimension.
												//TODO: check this iterator limit  !!!
		{
			fscanf(IN, "%d", &component_indicator);
			if (component_indicator==1) {  //then is on this component
				component_number = ii;
			}
		}
		
		if (component_number==-10) {
			printf("it appears the candidate point lies on NO COMPONENT.\n");
		}
		
		component_numbers[jj]=component_number;
	}
	
	fclose(IN);
	
	
	return component_numbers;
}






int write_member_points_sc(vec_mp point_to_write)
{
	FILE *OUT = NULL;
	int ii;
	
	
	remove("member_points");
	OUT = safe_fopen_write("member_points");
	
	
	vec_mp result;
	init_vec_mp(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"2\n\n");
	for(ii=0;ii<result->size;ii++)
	{
		print_mp(OUT,0,&result->coord[ii]);  fprintf(OUT,"\n");
	}
	
	
	comp_mp temp; init_mp(temp);
	fprintf(OUT,"\n");
	for(ii=0;ii<result->size;ii++)
	{
		conjugate_mp(temp, &result->coord[ii]);
		print_mp(OUT,0,temp);  fprintf(OUT,"\n");
	}
	fclose(OUT);
	
	clear_vec_mp(result);  clear_mp(temp);
	
	return 0;
}


int write_member_points_singlept(vec_mp point_to_write)
{
	FILE *OUT = NULL;
	
	
	
	OUT = safe_fopen_write("member_points");
	
	
	vec_mp result;
	init_vec_mp(result,0);
	dehomogenize(&result,point_to_write);
	
	fprintf(OUT,"1\n\n");
	for(int ii=0; ii<result->size; ii++){
		print_mp(OUT,0,&result->coord[ii]);
		fprintf(OUT,"\n");
	}
	
	
	fclose(OUT);
	
	
	clear_vec_mp(result);
	
	return 0;
}











void membership_test_input_file(boost::filesystem::path outputFile,
                                boost::filesystem::path funcInput,
                                boost::filesystem::path configInput,
                                int  tracktype)
{
	char ch;
	FILE *OUT = safe_fopen_write(outputFile), *IN = NULL;
	
	
	// setup configurations in OUT
	fprintf(OUT, "CONFIG\n");
	IN = safe_fopen_read(configInput);
	
	
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);

	fprintf(OUT, "TrackType: %d;\nDeleteTempFiles: 0;\nEND;\nINPUT\n",tracktype);

	
	// setup system in OUT
	IN = safe_fopen_read(funcInput);
	
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	fprintf(OUT, "END;\n");
	fclose(OUT);
	
	return;
}







int get_incidence_number(vec_mp test_point,
						 BertiniRealConfig & program_options,
						 boost::filesystem::path input_file)
{
	

	
	// setup input file
	int *declarations = NULL;
	partition_parse(&declarations, input_file, "func_input_real", "config_real",0); // the 0 means not self conjugate
	free(declarations);
	
	
	//check existence of the required witness_data file.
	FILE *IN = NULL;
	IN = safe_fopen_read("witness_data"); fclose(IN);
	
	
	
	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
	membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	
	//setup  member_points file, including both the first witness point, and its complex-conjugate
	write_member_points_singlept(test_point);
	
	
//	std::vector<std::string> command_line_options;
//	command_line_options.push_back("input_membership_test");
//	
	int blabla;
	parse_input_file("input_membership_test", &blabla);
	membershipMain(13423, blabla, 0, 1, 0);
	initMP(mpf_get_default_prec());
	
//	bertini_main_wrapper(command_line_options, 1, 0,0);
	
	
	std::vector<int> component_number = read_incidence_matrix();
	
	return component_number[0];
}







