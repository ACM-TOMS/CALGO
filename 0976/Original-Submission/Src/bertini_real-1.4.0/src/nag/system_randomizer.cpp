#include "nag/system_randomizer.hpp"







void SystemRandomizer::randomize(vec_d randomized_func_vals, mat_d randomized_jacobian,
								  vec_d func_vals, mat_d jacobian_vals,
								  comp_d hom_var)
{
	bool bail_out = false;
	
	if (!is_ready()) {
		std::cout << "trying to randomize_d, but is not set up!" << std::endl;
		bail_out = true;
	}
	
	if (func_vals->size != this->num_original_funcs) {
		std::cout << "mismatch in number of expected input functions (" << num_original_funcs << ") and actual inputted-number of functions (" << func_vals->size << ")." << std::endl;
		bail_out = true;
	}
	
	if (jacobian_vals->rows != this->num_original_funcs) {
		std::cout << "mismatch in size of expected input jacobian (" << num_original_funcs << ") and actual number of rows in jacobian (" << jacobian_vals->rows << ")." << std::endl;
		bail_out = true;
	}
	
	
	if (bail_out)
		br_exit(-9509);
	
	
	if (is_square()) {
		mat_cp_d(randomized_jacobian,jacobian_vals);
		vec_cp_d(randomized_func_vals,func_vals);
		return;
	}
	
	
	
	//ensure outputs are of correct size.
	increase_size_mat_d(randomized_jacobian, num_randomized_funcs, jacobian_vals->cols);
	randomized_jacobian->rows = num_randomized_funcs; randomized_jacobian->cols = jacobian_vals->cols;
	
	increase_size_vec_d(randomized_func_vals, num_randomized_funcs);
	randomized_func_vals->size = num_randomized_funcs;
	
	// do a little precomputation
	//0th entry is 1, having been set previously.
	for (int ii=1; ii<=max_degree_deficiency; ii++) {
		mul_d(&temp_homogenizer_d->coord[ii],&temp_homogenizer_d->coord[ii-1],hom_var);
	}
	
	increase_size_mat_d(temp_jac_d,num_original_funcs,jacobian_vals->cols);
	temp_jac_d->rows = num_original_funcs; temp_jac_d->cols = jacobian_vals->cols;
	
	// we do both function and jacobian randomization in the same loop for optimization.
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		// this must be done in a loop because each output (randomized function) has a different homogeneous structure.
		//
		//TODO: optimization could be done by looking at the previous function's structure and omitting previously done calculations (for the second and subsequent functions)
		
		
		// copy the current randomizer coefficients into a single row matrix
		for (int jj=0; jj<num_original_funcs; jj++) {
			set_d(&single_row_input_d->entry[0][jj], &randomizer_matrix_d->entry[ii][jj]); //TODO this is repetitive, and wasteful.  optimize away.
		} // at least is is just setting, not multiplying
		
		
		
		
		/////////////////
		//
		//  functions
		//
		//////////////
		increase_size_vec_d(temp_funcs_d, num_original_funcs);  temp_funcs_d->size = num_original_funcs;
		
		for (int jj=0; jj<num_original_funcs; jj++) {
			// structure_matrix[ii][jj] gives the degree deficiency
			if (structure_matrix[ii][jj]>0) { // if must homogenize at least one degree.
				mul_d(&temp_funcs_d->coord[jj],&func_vals->coord[jj], &temp_homogenizer_d->coord[ structure_matrix[ii][jj] ]);
			}
			else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
			{
				set_d(&temp_funcs_d->coord[jj], &func_vals->coord[jj]);
			}
			else // yeah...  maybe optimize away?
			{
				set_d(&temp_funcs_d->coord[jj], &func_vals->coord[jj]);
			}
			
		}
		
		mul_mat_vec_d(temp_vec_d, single_row_input_d, temp_funcs_d);
		set_d(&randomized_func_vals->coord[ii], &temp_vec_d->coord[0]);
		
		
		
		
		/////////////////
		//
		//  jacobian
		//
		//////////////
		
		
		for (int kk=0; kk<jacobian_vals->cols; kk++) { // kk indexes the variables in jacobian_vals (columns in jacobian)
			for (int jj=0; jj<num_original_funcs; jj++) { // jj indexes the original functions (rows in jacobian)
													  // structure_matrix[ii][jj] gives the degree deficiency
				if (structure_matrix[ii][jj]>0) { // must homogenize at least one degree.
					mul_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk], &temp_homogenizer_d->coord[ structure_matrix[ii][jj] ]);
				}
				else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
				{
					set_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				else
				{
					set_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				
				
			}
		}
		
		
		// actually randomize here
		mat_mul_d(temp_mat_d, single_row_input_d, temp_jac_d);
		
		
		//copy the output into the returned value.
		for (int kk=0; kk<jacobian_vals->cols; kk++) {
			set_d(&randomized_jacobian->entry[ii][kk], &temp_mat_d->entry[0][kk]);
		}
		
		
		// for the last step, the first variable is the hom_var, and it must use the product rule.
		for (int jj=0; jj<num_original_funcs; jj++) {
			if (structure_matrix[ii][jj]>0) { // only if we need to actually homogenize this function
											  // we abuse the temp_funcs here and use as a temp storage area.
				mul_d(&temp_funcs_d->coord[jj], //d•h^(d-1)
					  &integer_coeffs_d->coord[ structure_matrix[ii][jj] ], &temp_homogenizer_d->coord[ structure_matrix[ii][jj]-1 ]); // this could be optimized if degree deficiency == 1.
				mul_d(&temp_funcs_d->coord[jj], &temp_funcs_d->coord[jj], &func_vals->coord[jj]);
				//temp_funcs_ = d•h^(d-1)•f_{jj}
			}
			else
			{
				set_zero_d(&temp_funcs_d->coord[jj]);
			}
		}
		
		//randomize
		mul_mat_vec_d(temp_vec_d, single_row_input_d, temp_funcs_d); // this produces a single number as output.
																	 // M = R•(hommed_f)
																	 // combine for power rule
		add_d(&randomized_jacobian->entry[ii][0], &randomized_jacobian->entry[ii][0], &temp_vec_d->coord[0]);
		
	}
	
	
	
}


void SystemRandomizer::randomize(vec_mp randomized_func_vals, mat_mp randomized_jacobian,
								  vec_mp func_vals, mat_mp jacobian_vals,
								  comp_mp hom_var)
{
	bool bail_out = false;
	
	if (!is_ready()) {
		std::cout << "trying to randomize_d, but is not set up!" << std::endl;
		bail_out = true;
	}
	
	if (func_vals->size != this->num_original_funcs) {
		std::cout << "mismatch in number of expected input functions (" << num_original_funcs << ") and actual inputted-number of functions (" << func_vals->size << ")." << std::endl;
		bail_out = true;
	}
	
	if (jacobian_vals->rows != this->num_original_funcs) {
		std::cout << "mismatch in size of expected input jacobian (" << num_original_funcs << ") and actual number of rows in jacobian (" << jacobian_vals->rows << ")." << std::endl;
		bail_out = true;
	}
	
	
	if (bail_out)
		br_exit(-9510);
	
	
	
	if (is_square()) {
		mat_cp_mp(randomized_jacobian,jacobian_vals);
		vec_cp_mp(randomized_func_vals,func_vals);
		return;
	}
	
	//ensure outputs are of correct size.
	increase_size_mat_mp(randomized_jacobian, num_randomized_funcs, jacobian_vals->cols);
	randomized_jacobian->rows = num_randomized_funcs; randomized_jacobian->cols = jacobian_vals->cols;
	
	increase_size_vec_mp(randomized_func_vals, num_randomized_funcs);
	randomized_func_vals->size = num_randomized_funcs;
	
	// do a little precomputation
	//0th entry is 1, having been set previously.
	for (int ii=1; ii<=max_degree_deficiency; ii++) {
		mul_mp(&temp_homogenizer_mp->coord[ii],&temp_homogenizer_mp->coord[ii-1],hom_var);
	}
	
	increase_size_mat_mp(temp_jac_mp,num_original_funcs,jacobian_vals->cols);
	temp_jac_mp->rows = num_original_funcs;
	temp_jac_mp->cols = jacobian_vals->cols;
	
	// we do both function and jacobian randomization in the same loop for optimization.
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		// this must be done in a loop because each output (randomized function) has a different homogeneous structure.
		//
		// optimization could be done by looking at the previous function's structure and omitting previously done calculations (for the second and subsequent functions)
		
		
		// copy the current randomizer coefficients into a single row matrix
		for (int jj=0; jj<num_original_funcs; jj++) {
			set_mp(&single_row_input_mp->entry[0][jj], &randomizer_matrix_mp->entry[ii][jj]); // this is repetitive, and wasteful.  optimize away.
		} // at least is is just setting, not multiplying
		
		
		
		
		/////////////////
		//
		//  functions
		//
		//////////////
		
		increase_size_vec_mp(temp_funcs_mp, num_original_funcs);  temp_funcs_mp->size = num_original_funcs;
		
		for (int jj=0; jj<num_original_funcs; jj++) {
			// structure_matrix[ii][jj] gives the degree deficiency
			if (structure_matrix[ii][jj]>0) { // if must homogenize at least one degree.
				mul_mp(&temp_funcs_mp->coord[jj],&func_vals->coord[jj], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj] ]);
			}
			else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
			{
				set_mp(&temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
			}
			else // yeah...  maybe optimize away?
			{
				set_mp(&temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
			}
			
		}
		
		mul_mat_vec_mp(temp_vec_mp, single_row_input_mp, temp_funcs_mp);
		set_mp(&randomized_func_vals->coord[ii], &temp_vec_mp->coord[0]);
		
		
		
		
		/////////////////
		//
		//  jacobian
		//
		//////////////
		
		
		for (int kk=0; kk<jacobian_vals->cols; kk++) { // kk indexes the variables in jacobian_vals (columns in jacobian)
			for (int jj=0; jj<num_original_funcs; jj++) { // jj indexes the original functions (rows in jacobian)
													  // structure_matrix[ii][jj] gives the degree deficiency
				if (structure_matrix[ii][jj]>0) { // must homogenize at least one degree.
					mul_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj] ]);
				}
				else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
				{
					set_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				else
				{
					set_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				
				
			}
		}
		
		
		// actually randomize here
		mat_mul_mp(temp_mat_mp, single_row_input_mp, temp_jac_mp);
		
		
		//copy the output into the returned value.
		for (int kk=0; kk<jacobian_vals->cols; kk++) {
			set_mp(&randomized_jacobian->entry[ii][kk], &temp_mat_mp->entry[0][kk]);
		}
		
		
		// for the last step, the first variable is the hom_var, and it must use the product rule.
		for (int jj=0; jj<num_original_funcs; jj++) {
			if (structure_matrix[ii][jj]>0) { // only if we need to actually homogenize this function
											  // we abuse the temp_funcs here and use as a temp storage area.
				mul_mp(&temp_funcs_mp->coord[jj], //d•h^(d-1)
					   &integer_coeffs_mp->coord[ structure_matrix[ii][jj] ], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj]-1 ]); // this could be optimized if degree deficiency == 1.
				mul_mp(&temp_funcs_mp->coord[jj], &temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
				//temp_funcs_ = d•h^(d-1)•f_{jj}
			}
			else
			{
				set_zero_mp(&temp_funcs_mp->coord[jj]);
			}
		}
		
		//randomize
		mul_mat_vec_mp(temp_vec_mp, single_row_input_mp, temp_funcs_mp);
		// M = R•(hommed_f)
		// combine for power rule
		add_mp(&randomized_jacobian->entry[ii][0], &randomized_jacobian->entry[ii][0], &temp_vec_mp->coord[0]);
		
	}
}


void SystemRandomizer::change_prec(int new_prec)
{
	if (!is_ready()) {
		throw std::logic_error("trying to change precision when not set up!");
	}
	change_prec_mat_mp(randomizer_matrix_mp, new_prec);
	mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
	
	change_prec_vec_mp(temp_homogenizer_mp,new_prec);
	change_prec_vec_mp(temp_funcs_mp,new_prec);
	change_prec_mat_mp(temp_jac_mp,new_prec);
	
	
	change_prec_mat_mp(single_row_input_mp,new_prec);
	
	change_prec_vec_mp(integer_coeffs_mp,new_prec);
	change_prec_mat_mp(temp_mat_mp,new_prec);
	change_prec_vec_mp(temp_vec_mp,new_prec);
	
	
}

void SystemRandomizer::setup(int num_desired_rows, int num_funcs)
{
	if (num_desired_rows<=0) {
		throw std::logic_error("requested a randomizer for <= 0 output functions...");
	}
	if (num_funcs<=0) {
		throw std::logic_error("requested a randomizer for <= 0 input functions...");
	}
	
	if (num_desired_rows>num_funcs) {
		
		std::stringstream ss;
		ss << "requested randomizer which would UP-randomize." << std::endl;
		ss << num_desired_rows << " > " << num_funcs << ", which is not allowed" << std::endl;
		
		throw std::logic_error(ss.str());
	}
	
	setup_indicator = false; // reset;
	randomized_degrees.resize(0);
	original_degrees.resize(0);
	max_base_degree = 0;
	
	//get unique degrees
	int *unique_degrees = new int[num_funcs];
	
	
	FILE *IN = safe_fopen_read("deg.out"); //open the deg.out file for reading.
	int num_unique_degrees = 0;
	int occurrence_counter;
	int tempdegree;
	//TODO: this doesn't read correctly if there is more than one variable group.
	for (int ii=0; ii<num_funcs; ++ii) {
		fscanf(IN,"%d\n",&tempdegree); // read data
		original_degrees.push_back(tempdegree);
		
		occurrence_counter = 0; // set the counter for how many times the current degree has already been found.
		for (int jj=0; jj<ii; jj++) {
			if (original_degrees[jj]==tempdegree) { // if previously stored degree is same as current one
				occurrence_counter++; // increment counter
			}
		}
		
		if (occurrence_counter==0) { // if did not find already in list
			unique_degrees[num_unique_degrees] = tempdegree; // add to list of unique degrees.
			num_unique_degrees++; // have one more unique degree
		} // re: jj
		
		if (tempdegree>max_base_degree) {
			max_base_degree = tempdegree;
		}
		
	}// re: ii
	fclose(IN);
	
	
	
	num_randomized_funcs = num_desired_rows;
	num_original_funcs = num_funcs;
	
	

	square_indicator = false;
	
	//sort the unique degrees into decreasing order
	qsort(unique_degrees, num_unique_degrees, sizeof(int), compare_integers_decreasing);
	
	//count how many of each unique degree there are.
	int *num_of_each_degree = (int *) br_malloc(num_unique_degrees*sizeof(int));
	for (int ii=0; ii<num_unique_degrees; ii++) {
		num_of_each_degree[ii] = 0;
		for (int jj=0; jj<num_funcs; ++jj) {
			if (unique_degrees[ii]==original_degrees[jj]) {
				num_of_each_degree[ii]++;
			}
		}
	}
	
	
	
	
	//resize the matrix
	change_size_mat_mp(randomizer_matrix_full_prec,num_desired_rows,num_funcs);
	randomizer_matrix_full_prec->rows = num_desired_rows;
	randomizer_matrix_full_prec->cols = num_funcs;
	
	structure_matrix.resize(num_desired_rows);
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); ++iter) {
		iter->resize(num_funcs);
	}
	
	
	
	
	if (num_desired_rows==num_funcs) {
		for (int ii=0; ii<num_desired_rows; ii++) {
			randomized_degrees.push_back(original_degrees[ii]);
			
			for (int jj=0; jj<num_funcs; jj++) {
				if (ii==jj) {
					structure_matrix[ii][jj] = max_base_degree - original_degrees[ii];
				}
				else{
					structure_matrix[ii][jj] = 0;
				}
			}
		}
		make_matrix_ID_mp(randomizer_matrix_full_prec,num_funcs,num_funcs);
		square_indicator = true;
	}
	else{
		
		int counter = 0;
		int current_degree_index = 0; // start at the end
		
		for (int ii=0; ii<num_desired_rows; ii++) {
		
			counter++;
			if (counter>num_of_each_degree[current_degree_index]) {
				current_degree_index++;
				counter = 1;
			}
			
			int current_degree = unique_degrees[current_degree_index];
			randomized_degrees.push_back(current_degree);
			
			int encountered_current_degree = 0;
			for (int jj=0; jj<num_funcs; jj++) {
				if ( original_degrees[jj]<= current_degree ) {
					encountered_current_degree++;
					if (encountered_current_degree >= counter){
						get_comp_rand_real_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
						structure_matrix[ii][jj] = current_degree - original_degrees[jj];  // deficiency level
					}
					else{
						set_zero_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
						structure_matrix[ii][jj] = 0;
					}
				}
				else
				{
					set_zero_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
					structure_matrix[ii][jj] = 0;
				}
			}
			

		}
		
		square_indicator = false;
	}
	
	
	max_degree_deficiency = unique_degrees[0] - unique_degrees[num_unique_degrees-1];
	
	
	
	mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
	mat_mp_to_d(randomizer_matrix_d,randomizer_matrix_full_prec);
	
	
	free(num_of_each_degree);
	delete [] unique_degrees;
	
	
	
	
	
	
	setup_temps();
	
	
}


void SystemRandomizer::setup_temps()
{
	change_size_vec_d(integer_coeffs_d,max_degree_deficiency+1);  integer_coeffs_d->size = max_degree_deficiency+1;
	change_size_vec_mp(integer_coeffs_mp,max_degree_deficiency+1); integer_coeffs_mp->size = max_degree_deficiency+1;
	
	for (int ii=0; ii<=max_degree_deficiency; ii++) {
		integer_coeffs_d->coord[ii].r = ii; integer_coeffs_d->coord[ii].i = 0;
		set_zero_mp(&integer_coeffs_mp->coord[ii]); // initialize
		mpf_set_d(integer_coeffs_mp->coord[ii].r,ii); // set to integer value
	}
	
	
	change_size_mat_d(single_row_input_d,1,num_original_funcs);  single_row_input_d->rows = 1; single_row_input_d->cols = num_original_funcs;
	change_size_mat_mp(single_row_input_mp,1,num_original_funcs);  single_row_input_mp->rows = 1; single_row_input_mp->cols = num_original_funcs;
	
	change_size_vec_d(temp_homogenizer_d,max_degree_deficiency+1);  temp_homogenizer_d->size = max_degree_deficiency+1;
	change_size_vec_mp(temp_homogenizer_mp,max_degree_deficiency+1); temp_homogenizer_mp->size = max_degree_deficiency+1;
	set_one_d(&temp_homogenizer_d->coord[0]);
	set_one_mp(&temp_homogenizer_mp->coord[0]);
	
	
	
	
	setup_indicator = true;
}




void SystemRandomizer::send(int target, ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "SystemRandomizer::send" << std::endl;
#endif
	int sendme = setup_indicator;
	MPI_Send(&sendme,1,MPI_INT,target,UNUSED,mpi_config.comm());

	if (!setup_indicator) {
		std::cout << "bailing on sending upsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	buffer[0] = square_indicator;
	buffer[1] = num_randomized_funcs;
	buffer[2] = num_original_funcs;
	buffer[3] = max_base_degree;
	buffer[4] = max_degree_deficiency;
	MPI_Send(buffer,5,MPI_INT,target, SYSTEM_RANDOMIZER ,mpi_config.comm());
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_send = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_send];
	
	
	
	if (randomized_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = randomized_degrees.begin(); iter!=randomized_degrees.end(); iter++) {
			buffer2[cnt] = *iter;
			cnt++;
		}
		
	}
	
	int offset = num_randomized_funcs;
	if (original_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = original_degrees.begin(); iter!=original_degrees.end(); iter++) {
			buffer2[cnt+offset] = *iter;
			cnt++;
		}
	}
	
	
	offset += num_original_funcs;
	
	int cnt = 0;
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); iter++) {
		for (auto jter = iter->begin(); jter != iter->end(); jter++) {
			buffer2[cnt+offset] = *jter;
			cnt++;
		}
		
	}

	
	
	
	MPI_Send(buffer2, size_to_send, MPI_INT, target, SYSTEM_RANDOMIZER, mpi_config.comm());
	delete [] buffer2;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		send_mat_mp(randomizer_matrix_full_prec, target);
	}
	
	
	
}


void SystemRandomizer::receive(int source, ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "SystemRandomizer::receive" << std::endl;
#endif
	
	MPI_Status statty_mc_gatty;
	
	int recvme;
	MPI_Recv(&recvme,1,MPI_INT,source,UNUSED,mpi_config.comm(),&statty_mc_gatty);
	setup_indicator = recvme;
	if (!setup_indicator) {
		std::cout << "bailing from getting unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	MPI_Recv(buffer,5,MPI_INT,source,SYSTEM_RANDOMIZER,mpi_config.comm(),&statty_mc_gatty);
	
	square_indicator = buffer[0];
	num_randomized_funcs = buffer[1];
	num_original_funcs = buffer[2];
	max_base_degree = buffer[3];
	max_degree_deficiency = buffer[4];
	
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_receive = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_receive];
	
	
	
	
	MPI_Recv(buffer2, size_to_receive, MPI_INT, source, SYSTEM_RANDOMIZER, mpi_config.comm(),&statty_mc_gatty);
	
	int cnt = 0;
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		randomized_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	for (int ii=0; ii<num_original_funcs; ii++) {
		original_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	structure_matrix.resize(num_randomized_funcs);
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		for (int jj=0; jj<num_original_funcs; jj++) {
			structure_matrix[ii].push_back(buffer2[cnt]);
			cnt++;
		}
	}
	
	
	delete [] buffer2;
	
	change_size_mat_mp(randomizer_matrix_full_prec,num_randomized_funcs,num_original_funcs);
	randomizer_matrix_full_prec->rows = num_randomized_funcs;
	randomizer_matrix_full_prec->cols = num_original_funcs;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		receive_mat_mp(randomizer_matrix_full_prec, source);
	}
	
	
	setup_temps();
	// temps and things set up after you get the matrix and a few parameters.

	
}



void SystemRandomizer::bcast_send(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "SystemRandomizer::bcast_send" << std::endl;
#endif
	
	
	int sendme = setup_indicator;
	MPI_Bcast(&sendme,1,MPI_INT,mpi_config.head(),mpi_config.comm());
	
	if (!setup_indicator) {
		std::cout << "bailing on sending unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	buffer[0] = square_indicator;
	buffer[1] = num_randomized_funcs;
	buffer[2] = num_original_funcs;
	buffer[3] = max_base_degree;
	buffer[4] = max_degree_deficiency;
	MPI_Bcast(buffer,5,MPI_INT,mpi_config.head(),mpi_config.comm());
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_send = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_send];
	
	
	
	if (randomized_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = randomized_degrees.begin(); iter!=randomized_degrees.end(); iter++) {
			buffer2[cnt] = *iter;
			cnt++;
		}
		
	}
	
	int offset = num_randomized_funcs;
	if (original_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = original_degrees.begin(); iter!=original_degrees.end(); iter++) {
			buffer2[cnt+offset] = *iter;
			cnt++;
		}
	}
	
	
	offset += num_original_funcs;
	
	int cnt = 0;
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); iter++) {
		for (auto jter = iter->begin(); jter != iter->end(); jter++) {
			buffer2[cnt+offset] = *jter;
			cnt++;
		}
		
	}
	
	

	MPI_Bcast(buffer2, size_to_send, MPI_INT, mpi_config.head(), mpi_config.comm());
	delete [] buffer2;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		bcast_mat_mp(randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
	}
	
	

	
	
	
	
}


void SystemRandomizer::bcast_receive(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "SystemRandomizer::bcast_receive" << std::endl;
#endif
	
	
	
	
	
	int recvme;
	MPI_Bcast(&recvme,1,MPI_INT,mpi_config.head(),mpi_config.comm());
	setup_indicator = recvme;
	if (!setup_indicator) {
		std::cout << "bailing from getting unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	MPI_Bcast(buffer,5,MPI_INT,mpi_config.head(),mpi_config.comm());
	
	square_indicator = buffer[0];
	num_randomized_funcs = buffer[1];
	num_original_funcs = buffer[2];
	max_base_degree = buffer[3];
	max_degree_deficiency = buffer[4];
	
	
	delete[] buffer;

	
	
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_receive = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_receive];
	MPI_Bcast(buffer2, size_to_receive, MPI_INT, mpi_config.head(), mpi_config.comm());
	
	int cnt = 0;
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		randomized_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	for (int ii=0; ii<num_original_funcs; ii++) {
		original_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	structure_matrix.resize(num_randomized_funcs);
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		for (int jj=0; jj<num_original_funcs; jj++) {
			structure_matrix[ii].push_back(buffer2[cnt]);
			cnt++;
		}
	}
	
	
	delete [] buffer2;
	
	change_size_mat_mp(randomizer_matrix_full_prec,num_randomized_funcs,num_original_funcs);
	randomizer_matrix_full_prec->rows = num_randomized_funcs;
	randomizer_matrix_full_prec->cols = num_original_funcs;
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		bcast_mat_mp(randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
		mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
		mat_mp_to_d(randomizer_matrix_d,randomizer_matrix_full_prec);
	}
	
	
	setup_temps();
	// temps and things set up after you get the matrix and a few parameters.

}






