#include "nag/nid.hpp"





void NumericalIrreducibleDecomposition::populate(tracker_config_t * T)
{
	FILE *IN;
	

	
	IN = safe_fopen_read("witness_data");
	
	
	int num_nonempty_codims;

	
	fscanf(IN,"%d %d", &num_variables_, &num_nonempty_codims);
	
	
	
	vec_mp prev_approx;
	init_vec_mp(prev_approx,num_variables_); prev_approx->size = num_variables_;
	
	vec_mp temp_vec;
	init_vec_mp(temp_vec,num_variables_); temp_vec->size = num_variables_;
	
	
	//create counters and initialize to zeros
	std::vector<int> component_counter(num_nonempty_codims,0),codim_indicator(num_nonempty_codims,0);
	
	
	for (int ii=0; ii<num_nonempty_codims; ii++) {
		int current_codimension,num_points_this_dim;
		fscanf(IN,"%d %d",&current_codimension,&num_points_this_dim);
		
		
		//std::cout << "getting codimension " << current_codimension << " with " << num_points_this_dim << " points " << std::endl;
		
		codim_indicator[ii] = current_codimension;
		
		int current_dimension = num_variables_-1-current_codimension;
		
		nonempty_dimensions.push_back(current_dimension);
		
		for (int jj=0; jj<num_points_this_dim; jj++) {
			
			
			int precision;
			// last approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables_; kk++) {
				mpf_inp_str(temp_vec->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(temp_vec->coord[kk].i, IN, 10);
			}
			
			
			 
			// previous approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables_; kk++) {
				mpf_inp_str(prev_approx->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(prev_approx->coord[kk].i, IN, 10);
			}
			
			
			WitnessPointMetadata meta(current_dimension);
			
			meta.set_from_file(IN);
			

			int index = add_solution(temp_vec, meta);
			
			int num_existing_pts_this_comp = map_lookup_with_default( dimension_component_counter[current_dimension], meta.component_number(), 0 ); // the right hand 0 is the default if not found
			
			dimension_component_counter[current_dimension][meta.component_number()] = num_existing_pts_this_comp+1;
			index_tracker[current_dimension][meta.component_number()].push_back(index);
			
			
			if (meta.component_number()+1 > component_counter[ii]) {
				component_counter[ii]  = meta.component_number()+1;  // keeps track of the number of components per dimension
			}
			
			
			
		}
	}
	
	
	
	int should_be_minus_one;
	fscanf(IN,"%d",&should_be_minus_one);
	
	if (should_be_minus_one!=(-1)) {
		std::cerr << "did not parse top of file correctly.  got " << should_be_minus_one << " instead of -1\n";
		br_exit(-1);
	}
	
	int numbertype;
	fscanf(IN,"%d",&numbertype);
	
	if (numbertype==2 || numbertype==1) {
		change_prec_vec_mp(temp_vec,T->AMP_max_prec); 
	}
	else {
		change_prec_vec_mp(temp_vec,64);
	}
	
	
	mpq_t *temp_rat = (mpq_t *) br_malloc(2*sizeof(mpq_t)); // create and allocate two rationals.
	
	
	for (int mm=0; mm<num_nonempty_codims; mm++) {
		//BEGIN REPEATING BLOCK, one per nonempty codimension
		int current_codimension;
		current_codimension = codim_indicator[mm];
		int current_dimension = num_variables_-1-current_codimension;
		
		int num_rows_randomization, num_cols_randomization;
		fscanf(IN,"%d %d",&num_rows_randomization,&num_cols_randomization);
		
		mat_mp randomization_matrix;  init_mat_mp(randomization_matrix,num_rows_randomization,num_cols_randomization);
		for (int ii = 0; ii < num_rows_randomization; ii++) {
			for (int jj=0; jj<num_cols_randomization; jj++) {
				if (numbertype==2) {
					setup_comp_in_rat(temp_rat,IN);
					rat_to_mp(&randomization_matrix->entry[ii][jj],temp_rat); clear_rat(temp_rat);
				}
				else{
					mpf_inp_str(randomization_matrix->entry[ii][jj].r, IN, 10);
					mpf_inp_str(randomization_matrix->entry[ii][jj].i, IN, 10);
				}
			}
		}
		clear_mat_mp(randomization_matrix);
		
		
		//MATRIX W FOR HOMOGENIZATION
		//  same length as the randomization matrix
		homogenization_matrix_.resize(num_rows_randomization);
		for (int ii=0; ii<num_rows_randomization; ii++){
			homogenization_matrix_[ii].resize(num_cols_randomization);
			for (int jj=0; jj<num_cols_randomization; jj++){
				fscanf(IN,"%d",&homogenization_matrix_[ii][jj]);}
		
		}
		
		
		int num_entries_hom_patch_eqn;
		fscanf(IN,"%d",&num_entries_hom_patch_eqn);
		vec_mp homogenization_patch_eqn;  init_vec_mp(homogenization_patch_eqn,num_entries_hom_patch_eqn); homogenization_patch_eqn->size = num_entries_hom_patch_eqn;
		//		//VECTOR H FOR HOMOGENIZATION
		for (int ii = 0; ii<num_entries_hom_patch_eqn; ii++) {
			if (numbertype==2) {
				setup_comp_in_rat(temp_rat,IN);
				rat_to_mp(&temp_vec->coord[ii],temp_rat);  clear_rat(temp_rat);
			}
			else{
				mpf_inp_str(temp_vec->coord[ii].r, IN, 10);
				mpf_inp_str(temp_vec->coord[ii].i, IN, 10);
			}
		}
		clear_vec_mp(homogenization_patch_eqn);
		
		
		
		
		//   HOMVARCONST
		comp_mp hom_variable_constant;  init_mp(hom_variable_constant);
		setup_comp_in_rat(temp_rat,IN);
		rat_to_mp(hom_variable_constant,temp_rat); clear_rat(temp_rat);
		clear_mp(hom_variable_constant);
		
		
		
		//MATRIX B FOR LINEAR SLICE COEFFICIENTS
		int num_linears, num_lin_entries;
		fscanf(IN,"%d %d",&num_linears,&num_lin_entries);
		
		change_size_vec_mp(temp_vec,num_lin_entries);
		temp_vec->size = num_lin_entries;
		
		for (int ii = 0; ii < num_linears; ii++) {
			for (int jj=0; jj< num_lin_entries; jj++) {
				if (numbertype==2) {
					setup_comp_in_rat(temp_rat,IN);
					rat_to_mp(&temp_vec->coord[jj],temp_rat); clear_rat(temp_rat);
				}
				else{
					mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
					mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
				}
			}
			
			add_linear_w_meta(temp_vec, WitnessLinearMetadata(current_dimension));
			
		}
		
		
		
		
		int num_patch_coefficients;
		fscanf(IN,"%d",&num_patch_coefficients);
		
		
		
		change_size_vec_mp(temp_vec,num_patch_coefficients);
		temp_vec->size = num_patch_coefficients;
		
		//		// PATCH COEFFICIENTS
		for (int ii = 0; ii<num_patch_coefficients; ii++) {
			if (numbertype==2) {
				setup_comp_in_rat(temp_rat,IN);
				rat_to_mp(&temp_vec->coord[ii],temp_rat);  clear_rat(temp_rat);
			}
			else{
				mpf_inp_str(temp_vec->coord[ii].r, IN, 10);
				mpf_inp_str(temp_vec->coord[ii].i, IN, 10);
			}
			
		}
		
		add_patch_w_meta(temp_vec, WitnessPatchMetadata(current_dimension));
		
		//END REPEATED BLOCK (ONE FOR EACH NONEMPTY CODIM).
	}
	
	fclose(IN);
	
	free(temp_rat);
	
	
	
	
	clear_vec_mp(temp_vec);
	clear_vec_mp(prev_approx);
	
	
	
	
	
	
	
	
	
	
	
	
	return;
}






WitnessSet NumericalIrreducibleDecomposition::choose(BertiniRealConfig & options)
{
#ifdef functionentry_output
	std::cout << "NumericalIrreducibleDecomposition::choose" << std::endl;
#endif

	int target_dimension = options.target_dimension();
	int target_component = options.target_component();
	
	
	
	if (target_dimension == -1) {
		//try to get it by inference.  this is default behaviour.
		
		if (nonempty_dimensions.size()==1) {
			
			
			target_dimension = nonempty_dimensions[0];
			if (target_component==-1) { // want all components
				return best_possible_automatic_set(options);
			}
			else if(target_component==-2) // this is default
			{
				if (dimension_component_counter[target_dimension].size()==1) {
					return form_specific_witness_set(target_dimension,0);
				}
				else
				{
					return choose_set_interactive(options);
				}
			}
			else{ // want a specific witness set
				
//				std::cout << map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0) << " " << dimension_component_counter[target_dimension][target_component] << std::endl;
				
				
				if (map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0)==0) {
					std::cout << "you asked for component " << target_component << " (of dimension " << nonempty_dimensions[0] << ") which does not exist" << std::endl;
					WitnessSet W(num_variables());
					W.set_dimension(nonempty_dimensions[0]);
					W.set_component_number(-1);
					return W;
				}
				else{
					return form_specific_witness_set(nonempty_dimensions[0],target_component);
				}
			}
		}
		else{
			return choose_set_interactive(options);
		}
	}
	else{
		//want a specific dimension
		
		
		if (dimension_component_counter.find(target_dimension)==dimension_component_counter.end()) {
			std::cout << "there are no components of dimension " << options.target_dimension() << std::endl;
			
			WitnessSet W(num_variables());
			W.set_dimension(options.target_dimension());
			W.set_component_number(-1);
			return W;
		}
		else{
			if (target_component==-2) {
				if (dimension_component_counter[target_dimension].size()==1) { // this line is fragile because if something is looked up and doesn't exist, and entry is added.  this is why there is the function map_lookup_with_default.
					return form_specific_witness_set(target_dimension,0);
				}
				else{
					return choose_set_interactive(options);
				}
			}
			else if (target_component==-1) {
				return best_possible_automatic_set(options); // this may eventually call the interactive chooser
			}
			else{
				// want both a specific dimension and component number
				if (map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0)==0) {
					std::cout << "you asked for a component (" << target_component << ") which does not exist" << std::endl;
					WitnessSet W(num_variables());
					W.set_dimension(target_dimension);
					W.set_component_number(-1);
					return W;
				}
				else{
					return form_specific_witness_set(target_dimension,target_component);
				}
			}
		}
		
		
		
	}
	
	
}













WitnessSet NumericalIrreducibleDecomposition::best_possible_automatic_set(BertiniRealConfig & options)
{
#ifdef functionentry_output
	std::cout << "NumericalIrreducibleDecomposition::best_possible_automatic_set" << std::endl;
#endif
	int target_dimension = options.target_dimension();
	
	WitnessSet W(num_variables()); // create blank witness set
	W.set_dimension(target_dimension);
	W.set_component_number(-1);
	
	
	if (index_tracker.find(target_dimension)==index_tracker.end()) {
		std::cout << "must return empty set from auto constructor.  have no components of dimension " << target_dimension << std::endl;
		return W;
	}
	
	
	
	std::vector<int> components_with_no_deflations_needed;
	for (auto iter=index_tracker[target_dimension].begin(); iter!=index_tracker[target_dimension].end(); ++iter) {
		// iterate over components for target dimension
		if (iter->second.size()==0) {
			std::cout << "detected a witness set with no points.  this should be impossible...  by definition a witness set has at least one point in it." << std::endl;
			mypause();
		}
		else{
			if (point_metadata[iter->second[0]].num_deflations_needed()==0) {
				components_with_no_deflations_needed.push_back(iter->first);
			}
		}
		
	}
	
	if (components_with_no_deflations_needed.size()==0) { // everything needs deflation
		return choose_set_interactive(options);
	}
	
	std::cout << "checking for self-conjugate components" << std::endl;
	
	int sc_counter = 0;
	// need only copy the points, as all witness sets of a dimension have the same linears and patches.
	for (auto iter=components_with_no_deflations_needed.begin(); iter!=components_with_no_deflations_needed.end(); ++iter) {
		// iterate over deflation-free components for target dimension
		
		int current_index = index_tracker[target_dimension][*iter][0]; // guaranteed to exist, b/c nonempty.  already checked.
		
		if (checkSelfConjugate(point(current_index), options, options.input_filename())==true) {
			std::cout << "dim " << target_dimension << ", comp " << *iter << " is self-conjugate" << std::endl;
			for (int ii=0; ii<dimension_component_counter[target_dimension][*iter]; ++ii) {
				W.add_point( point(index_tracker[target_dimension][*iter][ii]) );
			}
			sc_counter++;
		}
	}
	
	
	
	
	
	if (sc_counter==0) { // found 0 self-conjugate components.
		std::cout << color::green() << "found 0 self-conjugate components, which do not need deflation." << color::console_default() << std::endl;
		return choose_set_interactive(options);
	}
	else{
		
		for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
			if (linear_metadata[ii].dimension() == target_dimension) {
				W.add_linear(linear(ii));
			}
		}
		
		for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
			if (patch_metadata[ii].dimension() == target_dimension) {
				W.add_patch(patch(ii));
			}
		}
		
		
		return W;
	}
		

	
}







// lets the user choose a set, and returns a copy of it.
WitnessSet NumericalIrreducibleDecomposition::choose_set_interactive(BertiniRealConfig & options)
{
#ifdef functionentry_output
	std::cout << "NumericalIrreducibleDecomposition::choose_set_interactive" << std::endl;
#endif
	int target_dimension = options.target_dimension();
	
	std::cout << "the nonempty dimensions:" << std::endl;
	for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
		std::cout << *iter << " ";
	}
	std::cout << std::endl;
	
	
		
	if (target_dimension==-1) { // display all dimensions, which is default
		
		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
				
				
				int first_index = *(jter->second.begin());
				std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity() << ", deflations needed: " << point_metadata[first_index].num_deflations_needed() << std::endl;
			}
			std::cout << std::endl;
		}
		
		std::set<int> valid_choices;
		for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
			valid_choices.insert(*iter);
		}
		
		target_dimension = get_int_choice("\n\nchoose a single dimension:\n",valid_choices);
		options.set_target_dimension(target_dimension);
	}
		
	
	std::cout << "dimension " << target_dimension << std::endl;
	for (auto jter = index_tracker[target_dimension].begin(); jter!= index_tracker[target_dimension].end(); ++jter) {
		int first_index = *(jter->second.begin());
		std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity() << ", deflations needed: " << point_metadata[first_index].num_deflations_needed() << ", degree " << dimension_component_counter[target_dimension][jter->first] << std::endl;
		
	}
	std::cout << std::endl;
	
	
	
	int num_components = get_int_choice("how many components would you like to decompose? (0 for all mult-one components)\n",0,dimension_component_counter[target_dimension].size());
	
	
	std::set<int> chosen_few;
	
	if (num_components==0) {
		for (unsigned int ii=0; ii<dimension_component_counter[target_dimension].size(); ++ii) {
			if (point_metadata[index_tracker[target_dimension][ii][0]].multiplicity()==1) {
				chosen_few.insert(ii);
			}
			
		}
	}
	else{
		std::set<int> valid_choices;
		for (unsigned int ii=0; ii<dimension_component_counter[target_dimension].size(); ++ii) {
			valid_choices.insert(ii);
		}
		
		for (int ii=0; ii<num_components; ii++) {
			int new_choice = get_int_choice("choose a component:\n",valid_choices);
			chosen_few.insert(new_choice);
			valid_choices.erase(new_choice);
		}
		
	}
	
	
	
	
	//3. return the set via a copy.
	WitnessSet W(num_variables());
	
	W.set_dimension(target_dimension);
	if (chosen_few.size()==1) {
		W.set_component_number(*chosen_few.begin());
	}
	else{
		W.set_component_number(-1);
	}
	
	for (auto iter=chosen_few.begin(); iter!=chosen_few.end(); ++iter) {
		
		//iterate over each point in the component
		for (auto jter=index_tracker[target_dimension][*iter].begin(); jter!=index_tracker[target_dimension][*iter].end(); ++jter) {
			W.add_point( point(*jter) );
		}
		
	}
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension() == target_dimension) {
			W.add_linear(linear(ii));
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension() == target_dimension) {
			W.add_patch(patch(ii));
		}
	}
	
	
	return W;
}





WitnessSet NumericalIrreducibleDecomposition::form_specific_witness_set(int dim, int comp)
{
	WitnessSet W(num_variables());
	
	W.set_dimension(dim);
	W.set_component_number(comp);
	
	
	// check to make sure actually have the set we need.
	if (dimension_component_counter.find(dim)==dimension_component_counter.end()) {
		std::cout << "trying to construct witness set of dimension " << dim << " but there are no components of that dimension." << std::endl;
		return W;
	}
	
	
	if (dimension_component_counter[dim].find(comp) == dimension_component_counter[dim].end()) {
		std::cout << "trying to constuct witness set for dimension " << dim << ", component " << comp << ", but that component does not exist" << std::endl;
		return W;
	}
	
	// ok, it exists.  lets form it.
	
	for (auto iter = index_tracker[dim][comp].begin(); iter!= index_tracker[dim][comp].end(); ++iter) {
		//iter points to an index into the vertices stored in the Vertex set.
		W.add_point( point(*iter) );
	}
	
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension() == dim) {
			W.add_linear(linear(ii));
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension() == dim) {
			W.add_patch(patch(ii));
		}
	}
	
	return W;
}



