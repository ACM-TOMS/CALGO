#include "nag/witness_set.hpp"








int WitnessSet::Parse(const boost::filesystem::path witness_set_file, const int num_vars)
{
	
	
	
	
	FILE *IN = safe_fopen_read(witness_set_file);
	
	
	int temp_num_patches, patch_size, temp_num_linears, temp_num_points, num_vars_in_linears;
	
	
	fscanf(IN, "%d %d %d", &temp_num_points, &dim_, &comp_num_); scanRestOfLine(IN);
	
	
	
	this->num_vars_ = num_vars;
	this->num_natty_vars_ = num_vars;
	
	
	vec_mp temp_vec;  init_vec_mp2(temp_vec, num_vars,1024); temp_vec->size = num_vars;
	
	for (int ii=0; ii < temp_num_points; ii++) {
		
		//read the witness points into memory
		for (int jj=0; jj < num_vars; ++jj) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10); // 10 is the base
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			
			scanRestOfLine(IN);
		}
		
		add_point(temp_vec);
	}
	
	
	
	
	
	fscanf(IN, "%d %d", &temp_num_linears, &num_vars_in_linears);  scanRestOfLine(IN);
	
	
	for (int ii=0; ii < temp_num_linears; ii++) {
		change_size_vec_mp(temp_vec,num_vars_in_linears);
		temp_vec->size = num_vars_in_linears;
		
		//read the witness linears into memory
		for (int jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		add_linear(temp_vec);
	}
	
	
	
	
	fscanf(IN, "%d %d", &temp_num_patches, &patch_size); scanRestOfLine(IN);
	
	if (temp_num_patches>1) {
		std::cerr << temp_num_patches << " patches detected.  this probably indicates a problem." << std::endl;
		std::cerr << "the file being read: " << witness_set_file << std::endl;
		std::cerr << "trying to read " << num_vars << " variables." << std::endl;
		mypause();
	}
	
	
	for (int ii=0; ii < temp_num_patches; ii++) {
		
		change_size_vec_mp(temp_vec,patch_size);
		temp_vec->size = patch_size;
		
		//read the patch into memory
		for (int jj=0; jj < patch_size; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		add_patch(temp_vec);
	}
	
	fclose(IN);
	
	
	clear_vec_mp(temp_vec);
	
	return 0;
}










void WitnessSet::write_homogeneous_coordinates(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_points()); // print the header line
	
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		vec_mp &curr_point = point(ii);
		for ( int jj=0; jj< curr_point->size; jj++) {
			print_mp(OUT,0,&( curr_point->coord[jj]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void WitnessSet::write_dehomogenized_coordinates(boost::filesystem::path filename) const
{

	
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%zu\n\n",num_points()); // print the header line
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (this->num_synth_vars()>0) {
			dehomogenize(&result,point(ii), num_natty_vars_);
		}
		else{
			dehomogenize(&result,point(ii));
		}
		
		for (int jj=0; jj<num_natty_vars_-1; jj++) {
			print_mp(OUT, 0, &result->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	clear_vec_mp(result);
	return;
}



void WitnessSet::write_dehomogenized_coordinates(boost::filesystem::path filename,std::set<unsigned int> indices) const
{
	
	for (auto ii=indices.begin(); ii!=indices.end(); ++ii) {
		if (*ii >= this->num_points()) {
			std::cout << "requested to print out-of-range point index " << *ii << " to a dehomogenized file." << std::endl;
			std::cout << "[this WitnessSet contains " << num_points() << " points.]" << std::endl;
			br_exit(66190);
		}
	}
	
	
	
	
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%lu\n\n",indices.size()); // print the header line
	for (auto ii=indices.begin(); ii!=indices.end(); ++ii) {
		if (this->num_synth_vars()>0) {
			dehomogenize(&result,this->point(*ii), num_natty_vars_);
		}
		else{
			dehomogenize(&result,this->point(*ii));
		}
		
		for (int jj=0; jj<num_natty_vars_-1; jj++) {
			print_mp(OUT, 0, &result->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	clear_vec_mp(result);
	return;
}



void WitnessSet::write_linears(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_linears()); // print the header line
	
	for (unsigned int ii=0; ii<num_linears(); ++ii) {
		for (int jj=0; jj<linear(ii)->size; jj++) {
			print_mp(OUT, 0, &linear(ii)->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void WitnessSet::print_patches(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_patches()); // print the header line
	
	for (unsigned int ii=0; ii<num_patches(); ++ii) {
		fprintf(OUT,"%d\n",patch(ii)->size);
		for (int jj=0; jj<patch(ii)->size; jj++) {
			print_mp(OUT, 0, &patch(ii)->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void WitnessSet::read_patches_from_file(boost::filesystem::path filename)
{
	
	FILE *IN = safe_fopen_read(filename);
	
	WitnessSet::reset_patches();
	int curr_num_patches;
	fscanf(IN,"%d\n",&curr_num_patches);
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		int curr_size;
		fscanf(IN,"%d\n",&curr_size);
		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			mpf_inp_str(temp_patch->coord[jj].r, IN, 10);
			mpf_inp_str(temp_patch->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		WitnessSet::add_patch(temp_patch);
	}
	
	clear_vec_mp(temp_patch);
	fclose(IN);
	
	return;
}



void WitnessSet::print_to_screen() const
{
	
	vec_mp dehom;  init_vec_mp(dehom,1); dehom->size = 1;
	
	std::stringstream varname;
	
	std::cout << "witness set has " << num_vars_ << " total variables, " << num_natty_vars_ << " natural variables." << std::endl;
	
	
	std::cout << "dim " << dim_ << ", comp " << comp_num_ << std::endl;
	std::cout << "input file name " << input_filename_ << std::endl;

	
	printf("******\n%zu points\n******\n",num_points());
	std::cout << color::green();
	for (unsigned ii=0; ii<num_points(); ii++) {
		
		dehomogenize(&dehom, point(ii), num_natty_vars_);
		
		varname << "point_" << ii;
		
		print_point_to_screen_matlab(dehom,varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::blue();
	printf("******\n%zu linears\n******\n",num_linears());
	
	for (unsigned ii=0; ii<num_linears(); ii++) {
		varname << "linear_" << ii;
		print_point_to_screen_matlab(linear(ii),varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::cyan();
	printf("******\n%zu patches\n******\n",num_patches());
	
	for (unsigned ii=0; ii<num_patches(); ii++) {
		varname << "patch_" << ii;
		print_point_to_screen_matlab(patch(ii),varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << "variable names:\n";
	for (unsigned ii=0; ii< num_var_names(); ii++) {
		std::cout << name(ii) << "\n";
	}
	printf("\n\n");
	
	
	clear_vec_mp(dehom);
}


void WitnessSet::print_to_file(boost::filesystem::path filename) const
{
	// print back into the same format we parse from.
	
	
	FILE *OUT = safe_fopen_write(filename);
	
	
	fprintf(OUT, "%zu %d %d\n\n", num_points(), dim_, comp_num_);
	


	for (unsigned int ii=0; ii < num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		//print the witness points into file
		for (int jj=0; jj < curr_point->size; ++jj) {
			mpf_out_str(OUT,10,0,curr_point->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,curr_point->coord[jj].i);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	
	fprintf(OUT, "%zu %d\n", num_linears(), num_vars_);
	
	
	for (unsigned int ii=0; ii < num_linears(); ii++) {

		for (int jj=0; jj < linear(ii)->size; jj++) {
			mpf_out_str(OUT,10,0,linear(ii)->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,linear(ii)->coord[jj].i);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	fprintf(OUT, "%zu %d\n", num_patches(), num_vars_);
	
	for (unsigned int ii=0; ii < num_patches(); ii++) {
		
		vec_mp & curr_patch = patch(ii);
		
		//read the patch into memory
		for (int jj=0; jj < curr_patch->size; jj++) {
			mpf_out_str(OUT,10,0,curr_patch->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,curr_patch->coord[jj].i);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	fclose(OUT);
	
	
	
	
	
	return;
}


void WitnessSet::only_natural_vars()
{
	WitnessSet::only_first_vars(num_natty_vars_);
}


void WitnessSet::only_first_vars(int num_vars)
{
	
	vec_mp tempvec;  init_vec_mp2(tempvec, num_vars, 1024);
	tempvec->size = num_vars;
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		
		for (int jj=0; jj<num_vars; jj++) {
			set_mp(&tempvec->coord[jj], &curr_point->coord[jj]);
		}
		
		change_size_vec_mp(curr_point, num_vars);  curr_point->size = num_vars;
		vec_cp_mp(curr_point, tempvec);
	}
	
	
	this->num_vars_ = num_vars;
	
	int patch_size_counter = 0, trim_from_here =  0;
	for (unsigned int ii=0; ii<num_patches(); ii++) {
		patch_size_counter += patch(ii)->size;
		if (patch_size_counter == num_vars)
		{
			trim_from_here = ii+1;
		}
	}
	
	if (trim_from_here==0) {
		std::cerr << "problem: the sum of the patch sizes never equalled the number of variables to trim to...\nhence, the trimming operation could not complete." << std::endl;
		this->print_to_screen();
		deliberate_segfault();
	}
	
	for (unsigned int ii=trim_from_here; ii<num_patches(); ii++) {
		clear_vec_mp(patch(ii));
	}
	
	patch_mp_ = (vec_mp *) br_realloc(patch_mp_, trim_from_here* sizeof(vec_mp));
	num_patches_ = trim_from_here;
	
	for (unsigned int ii=0; ii<num_linears(); ii++) {
		linear(ii)->size = num_vars;
	}
	
	clear_vec_mp(tempvec);
	return;
}



void WitnessSet::sort_for_real(tracker_config_t * T)
{
	

	
	
	int *real_indicator = new int[num_points()];
	int counter = 0;
	
	vec_mp result; init_vec_mp(result,num_natty_vars_-1);
	result->size = num_natty_vars_-1;
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		for (int jj=1; jj<num_natty_vars_; jj++) {
			div_mp(&result->coord[jj-1], &curr_point->coord[jj], &curr_point->coord[0]);
		}
		real_indicator[ii] = checkForReal_mp(result, T->real_threshold);
		
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	vec_mp *tempvec = (vec_mp *)br_malloc(counter * sizeof(vec_mp));
	
	counter = 0;  // reset
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		if (real_indicator[ii]==1) {
			
			init_vec_mp2(tempvec[counter],this->num_vars_,1024); tempvec[counter]->size = this->num_vars_;
			vec_cp_mp(tempvec[counter],curr_point);
			counter++;
		}
		else{
			
		}
	}
	
	clear_vec_mp(result);
	
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		clear_vec_mp(point(ii));
	}
	free(pts_mp_);
	
	pts_mp_ = tempvec;
	num_pts_ = counter;
	
	delete[] real_indicator;
	return;
}






// T is necessary for the tolerances.
void WitnessSet::sort_for_unique(tracker_config_t * T)
{
	
	if (num_vars_==0) {
		throw std::logic_error("sorting witness set with 0 variables for uniqueness");
	}
	
	int curr_uniqueness;
	int num_good_pts = 0;
	std::vector<int> is_unique;
	
	for (unsigned int ii = 0; ii<num_points(); ++ii) {
		vec_mp &curr_point = point(ii);
		
		curr_uniqueness = 1;
		
		int prev_size_1 = curr_point->size;  curr_point->size = num_natty_vars_; // cache and change to natural number
		
		for (unsigned int jj=ii+1; jj<num_points(); ++jj) {
			vec_mp & inner_point = point(jj);
			int prev_size_2 = inner_point->size; inner_point->size = num_natty_vars_; // cache and change to natural number
			if ( isSamePoint_homogeneous_input(curr_point,inner_point,T->final_tol_times_mult) ){
				curr_uniqueness = 0;
			}
			inner_point->size = prev_size_2; // restore
		}
		curr_point->size = prev_size_1; // restore
		
		if (curr_uniqueness==1) {
			is_unique.push_back(1);
			num_good_pts++;
		}
		else {
			is_unique.push_back(0);
		}
		
	}
	
	
	
	vec_mp *transferme = (vec_mp *)br_malloc(num_good_pts*sizeof(vec_mp));
	int counter = 0;
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (is_unique[ii]==1) {
			
			init_vec_mp2(transferme[counter],num_vars_,1024);  transferme[counter]->size = num_vars_;
			vec_cp_mp(transferme[counter], point(ii));
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		std::logic_error("counter mismatch");
	}
	
	if (num_points()>0) {
		for (unsigned int ii=0; ii<num_points(); ii++) {
			clear_vec_mp(point(ii));
		}
		free(pts_mp_);
	}
	
	
	num_pts_ = num_good_pts;
	pts_mp_ = transferme;
	
	return;
}





void WitnessSet::sort_for_inside_sphere(comp_mp radius, vec_mp center)
{
	
	
	
	
	int num_good_pts = 0;
	
	
	std::vector<int> is_ok;
	
	vec_mp temp_vec; init_vec_mp(temp_vec,0);
	comp_mp temp; init_mp(temp);
	
	for (unsigned int ii = 0; ii<num_points(); ++ii) {
		
		
		
		dehomogenize(&temp_vec, point(ii),num_natty_vars_);
		temp_vec->size = center->size;
		
		norm_of_difference(temp->r, temp_vec, center);
		if ( mpf_cmp(temp->r, radius->r) < 0   ){
			is_ok.push_back(1);
			num_good_pts++;
		}
		else
		{
			is_ok.push_back(0);
		}
		
		
		
	}
	
	
	
	vec_mp *transferme = (vec_mp *)br_malloc(num_good_pts*sizeof(vec_mp));
	int counter = 0;
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (is_ok[ii]==1) {
			init_vec_mp2(transferme[counter],0,1024);  transferme[counter]->size = 0;
			vec_cp_mp(transferme[counter], point(ii));
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		br_exit(271);
	}
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		clear_vec_mp(point(ii));
	}
	free(pts_mp_);
	
	
	num_pts_ = num_good_pts;
	pts_mp_ = transferme;
	
	clear_vec_mp(temp_vec);
	clear_mp(temp);
	return;
}




void WitnessSet::merge(const WitnessSet & W_in, tracker_config_t * T)
{
	
	//error checking first
	if ( (num_vars_==0) && (W_in.num_vars_!=0) && (num_points()==0)) {
		num_vars_ = W_in.num_vars_;
		num_natty_vars_ = W_in.num_natty_vars_;
	}

	
	
	
	if (W_in.num_natty_vars_ != this->num_natty_vars_) {
		std::stringstream ss;
		ss << "merging two witness sets with differing numbers of natural variables. "<< W_in.num_natural_variables() <<" merging set, "<< this->num_natural_variables() << " existing\n",
		throw std::logic_error(ss.str());
	}
	
	//just mindlessly add the linears.  up to user to ensure linears get merged correctly.  no way to know what they want...
	for (unsigned int ii = 0; ii<W_in.num_linears(); ii++) {
		WitnessSet::add_linear(W_in.linear(ii));
	}
	
	
	for (unsigned int ii = 0; ii<W_in.num_points(); ii++) {
		int is_new = 1;
		vec_mp & in_point = W_in.point(ii);
		for (unsigned int jj = 0; jj<num_points(); jj++){
			vec_mp & curr_point = this->point(jj);
			
			//cache the sizes
			int in_size = in_point->size;
			int curr_size = curr_point->size;
			
			
			in_point->size = this->num_natural_variables();
			curr_point->size = this->num_natural_variables();
			
			if (isSamePoint_homogeneous_input(curr_point, in_point, T->final_tol_times_mult)) {
				is_new = 0;
				
				in_point->size = in_size;
				curr_point->size = curr_size;
				break;
			}
			
			// restore the sizes
			in_point->size = in_size;
			curr_point->size = curr_size;
		}
		
		if (is_new==1)
			WitnessSet::add_point( (in_point) );
	}
	
	
	for (unsigned int ii = 0; ii<W_in.num_patches(); ii++) {
		int is_new = 1;
		
		for (unsigned int jj = 0; jj<this->num_patches(); jj++){
			if ( this->patch(jj)->size ==  W_in.patch(ii)->size) {
				if (isSamePoint_inhomogeneous_input(patch(jj), W_in.patch(ii), T->final_tol_times_mult)) {
					is_new = 0;
					break;
				}
			}
		}
		
		if (is_new==1)
			WitnessSet::add_patch(W_in.patch(ii));
	}
	
	
	
	return;
}//re: merge_witness_sets
















void WitnessSet::send(ParallelismConfig & mpi_config, int target) const
{
    //send all the data to the other party
    

//
//	
//	std::vector< std::string > variable_names;
//	
//	boost::filesystem::path input_filename;
//	function input_file;
//	// end data members
    
    int *buffer = (int *) br_malloc(8*sizeof(int));
    buffer[0] = dimension();
    buffer[1] = component_number();
    buffer[2] = incidence_number();
    buffer[3] = num_variables();
    buffer[4] = num_natural_variables();
    buffer[5] = num_points();
    buffer[6] = num_linears();
    buffer[7] = num_patches();
    
    MPI_Send(buffer, 8, MPI_INT, WITNESS_SET, target, mpi_config.comm());
    
    free(buffer);
    
    for (unsigned int ii=0; ii<num_linears(); ii++) {
        send_vec_mp( linear(ii),target);
    }
    for (unsigned int ii=0; ii<num_patches(); ii++) {
        send_vec_mp( patch(ii),target);
    }
    for (unsigned int ii=0; ii<num_points(); ii++) {
        send_vec_mp( point(ii) ,target);
    }
    
    char * namebuffer = (char *) br_malloc(1024*sizeof(char));
    
	
    free(namebuffer);
    return;
}

void WitnessSet::receive(int source, ParallelismConfig & mpi_config)
{
    MPI_Status statty_mc_gatty;
    
	int *buffer = new int[8];
    
    
    MPI_Recv(buffer, 8, MPI_INT, WITNESS_SET, source, mpi_config.comm(), &statty_mc_gatty);
    
    set_dimension(buffer[0]);
    comp_num_ = buffer[1];
    incid_num_ = buffer[2];
    num_vars_ = buffer[3];
    num_natty_vars_ = buffer[4];
    unsigned int temp_num_pts = buffer[5];
    unsigned int temp_num_linears = buffer[6];
    unsigned int temp_num_patches = buffer[7];
    
	delete [] buffer;
	
    vec_mp tempvec; init_vec_mp2(tempvec,0,1024);
    
    for (unsigned int ii=0; ii<temp_num_linears; ii++) {
        receive_vec_mp(tempvec,source);
        add_linear(tempvec);
    }
    
    for (unsigned int ii=0; ii<temp_num_patches; ii++) {
        receive_vec_mp(tempvec,source);
        add_patch(tempvec);
    }
    
	for (unsigned int ii=0; ii<temp_num_pts; ii++) {
        receive_vec_mp(tempvec,source);
        add_point(tempvec);
    }
    
    clear_vec_mp(tempvec);
    return;
}







