#include "decompositions/decomposition.hpp"








int Decomposition::add_witness_set(const WitnessSet & W, VertexType add_type, VertexSet & V)
{
#ifdef functionentry_output
	std::cout << "Decomposition::add_witness_set" << std::endl;
#endif
	
	
    V.set_curr_input(W.input_filename());
    
    Vertex temp_vertex;
    temp_vertex.set_type(add_type);
    
    for (unsigned int ii=0; ii<W.num_points(); ii++) {
        vec_cp_mp(temp_vertex.point(), W.point(ii));
        this->index_in_vertices_with_add(V, temp_vertex);
    }
    
    return 0;
}


int Decomposition::add_vertex(VertexSet & V, Vertex source_vertex)
{
#ifdef functionentry_output
	std::cout << "Decomposition::add_vertex" << std::endl;
#endif
	
	int current_index = V.add_vertex(source_vertex);
	
	
	return current_index;
}






int Decomposition::index_in_vertices(VertexSet & V,
									 vec_mp testpoint) const
{
#ifdef functionentry_output
	std::cout << "Decomposition::index_in_vertices" << std::endl;
#endif
	int index = -1;
	
	
	// first we search the non-removed points.
	index = V.search_for_point(testpoint);
    
    
	return index;
}






int Decomposition::index_in_vertices_with_add(VertexSet &V,
											  Vertex vert)
{
	int index = Decomposition::index_in_vertices(V, vert.point());
	
	if (index==-1) {
		index = Decomposition::add_vertex(V, vert);
	}
	
	return index;
	
}











int Decomposition::setup(boost::filesystem::path INfile)
{
	boost::filesystem::path directoryName = INfile.parent_path();

	std::stringstream converter;
	std::string tempstr;
	std::ifstream fin(INfile.c_str());
	
	
	getline(fin, tempstr);
	set_input_filename(directoryName / tempstr);
	
	getline(fin, tempstr);
	converter << tempstr;
	converter >> this->num_variables_ >> this->dim_;
	converter.clear(); converter.str("");
	

	comp_mp temp; init_mp2(temp, 1024);
	vec_mp tempvec; init_vec_mp2(tempvec, this->num_variables_,1024);
	tempvec->size = this->num_variables_;
	

	for (int ii=0; ii<dimension(); ii++) {
		int temp_size;
		fin >> temp_size;
		
		change_size_vec_mp(tempvec,temp_size);  tempvec->size = temp_size;
		for (int jj=0;jj<temp_size;jj++)
		{
			std::string re, im;
			fin >> re >> im ;
			mpf_set_str(tempvec->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(tempvec->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		Decomposition::add_projection(tempvec);
	}
	
	
	
	
	int curr_num_patches = -1;
	fin >> curr_num_patches;
	
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		

		int curr_size = -1;
		fin >> curr_size;

		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			std::string re, im;
			fin >> re >> im; // this line is correct
		
			mpf_set_str(temp_patch->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(temp_patch->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		Decomposition::add_patch(temp_patch);
	}
	

	fin >> tempstr;
	mpf_set_str(temp->r, const_cast<char *>(tempstr.c_str()), 10);

	fin >> tempstr;
	mpf_set_str(temp->i, const_cast<char *>(tempstr.c_str()), 10);
	set_sphere_radius(temp);

	int curr_size;
	fin >> curr_size;
	if (curr_size!=num_variables_-1)
		std::cout << "sphere has incorrect number of variables (" << curr_size << ") when reading decomposition into memory\n";

	change_size_vec_mp(tempvec,curr_size); tempvec->size = curr_size;
	for (int jj=0; jj<curr_size; jj++) {
		std::string re, im;
		fin >> re >> im; // this line is correct
	
		mpf_set_str(tempvec->coord[jj].r, const_cast<char *>(re.c_str()), 10);
		mpf_set_str(tempvec->coord[jj].i, const_cast<char *>(im.c_str()), 10);
	}
	set_sphere_center(tempvec);


	fin >> curr_size;
	change_size_vec_mp(tempvec,curr_size); tempvec->size = curr_size;
	for (int jj=0; jj<curr_size; jj++) {
		std::string re, im;
		fin >> re >> im; // this line is correct
	
		mpf_set_str(tempvec->coord[jj].r, const_cast<char *>(re.c_str()), 10);
		mpf_set_str(tempvec->coord[jj].i, const_cast<char *>(im.c_str()), 10);
	}
	SetCritSliceValues(tempvec);



	if (fin.eof())
		std::cout << "premature end of file when reading decomposition from file";

	fin.close();
	

	clear_vec_mp(tempvec);
	clear_vec_mp(temp_patch);
	clear_mp(temp);
	
	
	return 0;
}






void Decomposition::print(boost::filesystem::path base) const
{
	
#ifdef functionentry_output
	std::cout << "Decomposition::print" << std::endl;
#endif
	
	if (dimension() != num_curr_projections()) {
//		std::cout << "Decomposition was short projections\nneeded	" << this->dimension << " but had " << num_curr_projections << std::endl;;
	}
	
	boost::filesystem::create_directory(base);
	
	FILE *OUT = safe_fopen_write(base / "decomp");
	
	fprintf(OUT,"%s\n",input_filename().filename().c_str());
	
	fprintf(OUT,"%d %d\n\n",num_variables(), dimension());
	
	
	for (int ii=0; ii<num_curr_projections(); ii++) {
		fprintf(OUT,"%d\n",pi_[ii]->size);
		for(int jj=0;jj<pi_[ii]->size;jj++)
		{
			print_mp(OUT, 0, &pi_[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fprintf(OUT,"%zu\n\n",this->num_patches()); // print the header line
	
	for (unsigned ii=0; ii<this->num_patches(); ++ii) {
		vec_mp & curr_patch = patch(ii);
		fprintf(OUT,"%d\n",curr_patch->size);
		for (int jj=0; jj<curr_patch->size; jj++) {
			print_mp(OUT, 0, &curr_patch->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
    
    fprintf(OUT,"\n\n");
    
    print_mp(OUT, 0, const_cast<_comp_mp*>(this->sphere_radius_));
    
    fprintf(OUT, "\n%d\n",this->sphere_center_->size); 
    for (int jj=0; jj<this->sphere_center_->size; jj++) {
        print_mp(OUT, 0, &this->sphere_center_->coord[jj]);
        fprintf(OUT, "\n");
    }
    fprintf(OUT,"\n");
    
    fprintf(OUT, "\n%d\n",this->crit_slice_values->size);
    for (int jj=0; jj<this->crit_slice_values->size; jj++) {
        print_mp(OUT, 0, &this->crit_slice_values->coord[jj]);
        fprintf(OUT, "\n");
    }
	fclose(OUT);
	
	
}




int Decomposition::read_sphere(const boost::filesystem::path & bounding_sphere_filename)
{
#ifdef functionentry_output
	std::cout << "Decomposition::read_sphere" << std::endl;
#endif
	
	if (num_variables_<2) {
		std::cout << "during read of sphere, Decomposition of dimension	" << dimension() << " has " << num_variables() << " variables!" << std::endl;
		mypause();
	}

	change_size_vec_mp(this->sphere_center_, num_variables()-1); //destructive resize
	sphere_center_->size = num_variables()-1;
	
	
	FILE *IN = safe_fopen_read(bounding_sphere_filename);
	
	mpf_inp_str(sphere_radius_->r, IN, 10);
	mpf_set_str(sphere_radius_->i,"0",10);
	
	
	for (int jj=1; jj<num_variables(); jj++) {
		mpf_inp_str(sphere_center_->coord[jj-1].r, IN, 10);
		mpf_set_str(sphere_center_->coord[jj-1].i,"0",10);
	}
	
	
	fclose(IN);
	
	have_sphere_ = true;
	
	
	return SUCCESSFUL;
}



void Decomposition::compute_sphere_bounds(const WitnessSet & W_crit)
{
	
#ifdef functionentry_output
	std::cout << "Decomposition::compute_sphere_bounds" << std::endl;
#endif

	
	int num_vars = W_crit.num_natural_variables()-1;
	
	change_size_vec_mp(this->sphere_center_, num_vars); //destructive resize
	sphere_center_->size = num_vars;
	
	if (W_crit.num_points() == 0) {
		set_zero_mp(sphere_radius_);
		mpf_set_str(sphere_radius_->r, "3.0",10);
		for (int ii=0; ii<num_vars; ii++) {
			set_zero_mp(&sphere_center_->coord[ii]);
		}
		this->have_sphere_ = true;
		return;
	}
	
	vec_mp(temp_vec); init_vec_mp2(temp_vec,0,1024);
	if (W_crit.num_points()==1)
	{
		set_zero_mp(sphere_radius_);
		mpf_set_str(sphere_radius_->r, "3.0",10);
		dehomogenize(&temp_vec, W_crit.point(0));
		
		for (int ii=0; ii<num_vars; ii++) {
			set_mp(&sphere_center_->coord[ii], &temp_vec->coord[ii]);
		}

		real_threshold(sphere_center_, 1e-13);

		clear_vec_mp(temp_vec);
		this->have_sphere_ = true;
		return;
	}
	
	
	//	W_crit.print_to_screen();
	
	
	
	
	comp_mp temp_rad; init_mp2(temp_rad,1024);
	set_zero_mp(temp_rad);
	
	set_one_mp(sphere_radius_);
	neg_mp(sphere_radius_,sphere_radius_); // set to impossibly low value.
	
	
	comp_mp temp_mp;  init_mp2(temp_mp,1024);
	
	
	vec_mp(cumulative_sum); init_vec_mp2(cumulative_sum,num_vars,1024);
	cumulative_sum->size = num_vars;
	
	for (int ii=0; ii<num_vars; ii++) {
		set_zero_mp(&cumulative_sum->coord[ii]);
	}
	
	
	for (unsigned int ii=0; ii<W_crit.num_points(); ii++) {
		dehomogenize(&temp_vec, W_crit.point(ii), num_vars+1);
		temp_vec->size = num_vars;
		add_vec_mp(cumulative_sum, cumulative_sum, temp_vec);
	}
	
	set_zero_mp(temp_mp);
	mpf_set_d(temp_mp->r, double(W_crit.num_points()));
	
	
	
	for (int ii=0; ii<num_vars; ii++) {
		div_mp(&sphere_center_->coord[ii], &cumulative_sum->coord[ii], temp_mp);
	}
	
	
	
	for (unsigned int ii=0; ii<W_crit.num_points(); ii++) {
		dehomogenize(&temp_vec, W_crit.point(ii), num_vars+1);
		temp_vec->size = num_vars;
		vec_sub_mp(temp_vec, temp_vec, sphere_center_);
		
		
		twoNormVec_mp(temp_vec, temp_mp);
		mpf_abs_mp(temp_rad->r, temp_mp);
		
		
		if (mpf_cmp(sphere_radius_->r, temp_rad->r) < 0){
			set_mp(sphere_radius_, temp_rad);
		}
	}
	
	
	
	mpf_set_str(temp_rad->r,"2.0",10);
	mpf_set_str(temp_rad->i,"0.0",10);
	mul_mp(sphere_radius_,temp_rad,sphere_radius_);  // double the radius to be safe.
	
	
	clear_mp(temp_mp); clear_mp(temp_rad);
	clear_vec_mp(temp_vec); clear_vec_mp(cumulative_sum);
	
	
	
	this->have_sphere_ = true;

	real_threshold(sphere_center_,1e-13);


	return;
}



void Decomposition::copy_data_from_witness_set(const WitnessSet & W)
{
	// set some member information.
	set_input_filename(W.input_filename());
	set_num_variables(W.num_variables());
	set_component_number(W.component_number());
	
	this->W_ = W;
}

void Decomposition::output_main(const boost::filesystem::path base) const
{
#ifdef functionentry_output
	std::cout << "Decomposition::output_main" << std::endl;
#endif

	FILE *OUT;
	boost::filesystem::path backupdir = base;
	backupdir += "_bak";
	if (boost::filesystem::exists(base)) {
		
		if (boost::filesystem::exists(backupdir)) {
			boost::filesystem::remove_all(backupdir);
		}
		boost::filesystem::rename(base, backupdir);
	}
	boost::filesystem::create_directory(base);
	
	
	copyfile("witness_data",base / "witness_data"); // this is wrong for nested decompositions.
	
	W_.print_to_file(base / "witness_set");
	
	
// TODO:  this should be a write call, not a copy!
	if (input_filename().filename().string().compare("unset")) {
		copyfile(input_filename(), base / input_filename().filename());
	}
	else{
//		std::cout << "not copying inputfile because name was unset -- " << input_filename << std::endl;
	}
	
	
	this->print(base); // using polymorphism and virtualism here!
	
	PrintPointTypeMapping(base / "vertex_types");


	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%s\n",base.c_str());
	fprintf(OUT,"%d\n",2);//remove this
	fprintf(OUT,"%d\n",dimension());
	fclose(OUT);
	
	
//	if (boost::filesystem::exists(backupdir)) {
//		boost::filesystem::remove_all(backupdir);
//	}
}



void Decomposition::send(int target, ParallelismConfig & mpi_config) const
{
#ifdef functionentry_output
	std::cout << "Decomposition::send" << std::endl;
#endif

	
	int * buffer2;
	
	
	//pack and send numbers of things.
	buffer2 = new int[7];
	buffer2[0] = num_variables();
	buffer2[1] = dimension();
	buffer2[2] = component_number();
	buffer2[3] = num_curr_projections();

	buffer2[4] = num_patches();
	buffer2[5] = int(have_sphere_);
	int strleng = input_filename().string().size() + 1;
	buffer2[6] = strleng;

	MPI_Send(buffer2, 7, MPI_INT, target, 6, mpi_config.comm());
	delete [] buffer2;
	
	
	
	
	
	if (num_curr_projections()>0) {
		for (int ii=0; ii<num_curr_projections(); ii++) {
			send_vec_mp(pi_[ii],target);
		}
	}
	
	
	
	
	if ( num_patches()>0) {

		for (unsigned int ii=0; ii<num_patches(); ii++) {
			send_vec_mp(patch(ii),target);
		}
	}
	
	
	
	
	if (have_sphere_) {

		send_vec_mp(sphere_center_,target);
		send_comp_mp(sphere_radius_,target);
	}
	
	
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		memcpy(buffer, input_filename().c_str(), strleng);
		MPI_Send(buffer, strleng, MPI_CHAR, target, 7, mpi_config.comm());
		delete [] buffer;
	}
	
	randomizer_->send(target,mpi_config);

	send_vec_mp(crit_slice_values, target);
	
	return;
}



void Decomposition::receive(int source, ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "Decomposition::receive" << std::endl;
#endif

	
	
	MPI_Status statty_mc_gatty;
	
	
	
	vec_mp tempvec;  init_vec_mp2(tempvec, 0, 1024);
	
	int * buffer2;
	
	
	
	
	
	
	buffer2 = new int[7];
	
	MPI_Recv(buffer2, 7, MPI_INT, source, 6, mpi_config.comm(), &statty_mc_gatty);
	set_num_variables(buffer2[0]);
	dim_ = buffer2[1];
	comp_num_ = buffer2[2];
	int temp_num_projections = buffer2[3];
	int temp_num_patches = buffer2[4];
	have_sphere_ = bool(buffer2[5]);
	int strleng = buffer2[6];

	delete [] buffer2;
	
	
	
	
	if (temp_num_projections>0) {

		for (int ii=0; ii<temp_num_projections; ii++) {
			receive_vec_mp(tempvec,source);
			add_projection(tempvec);
		}
	}
	
	
	
	

	
	
	if (temp_num_patches>0) {
		for (int ii=0; ii<temp_num_patches; ii++) {
			receive_vec_mp(tempvec,source);
			add_patch(tempvec);
		}
	}
	
	
	
	if (have_sphere_) {
		receive_vec_mp(sphere_center_,source);
		receive_comp_mp(sphere_radius_,source);
	}
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		MPI_Recv(buffer, strleng, MPI_CHAR, source, 7, mpi_config.comm(), &statty_mc_gatty);
		
		set_input_filename(buffer);
		delete [] buffer;
	}
	
	

	

	clear_vec_mp(tempvec);
	
	randomizer_->receive(source,mpi_config);
	
	receive_vec_mp(crit_slice_values, source);

	return;
}










