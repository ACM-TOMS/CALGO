#include "containers/vertex_set.hpp"





int VertexSet::search_for_point(vec_mp testpoint)
{
    int index = -1;
	
    index = search_for_active_point(testpoint);
    
    if (index==-1) {
        index = search_for_removed_point(testpoint);
    }
    
    return index;
}


int VertexSet::search_for_active_point(vec_mp testpoint)
{

    // dehomogenize the testpoint into the internal temp container.
    for (int jj=1; jj<num_natural_variables_; jj++) {
		div_mp(&checker_1_->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
	//	WTB: a faster comparison search.
	int current_index = 0;
	for (auto curr_vert=vertices_.begin(); curr_vert!=vertices_.end(); curr_vert++) {
		
		
		
		if (!curr_vert->is_removed()) {
			vec_mp& current_point = curr_vert->point();
			
			// dehomogenize the current point under investigation
			for (int jj=1; jj<num_natural_variables_; jj++) {
				div_mp(&checker_2_->coord[jj-1], &(current_point)->coord[jj], &(current_point)->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1_, checker_2_, same_point_tolerance_)){
				return current_index;
			}
		}
		current_index++;
	}
    
    return -1;
}



int VertexSet::search_for_removed_point(vec_mp testpoint)
{
    
    // dehomogenize the testpoint into the internal temp container.
    
    for (int jj=1; jj<num_natural_variables_; jj++) {
		div_mp(&checker_1_->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
    //	WTB: a faster comparison search.
	int current_index = 0;
	for (auto curr_vert = vertices_.begin(); curr_vert!=vertices_.end(); curr_vert++) {
		
		
		if (vertices_[current_index].is_removed()) {
			
			vec_mp & current_point = vertices_[current_index].point();
			
			for (int jj=1; jj<num_natural_variables_; jj++) {
				div_mp(&checker_2_->coord[jj-1], &(current_point)->coord[jj], &(current_point)->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1_, checker_2_, same_point_tolerance_)){
				return current_index;
			}
			
		}
		current_index++;
	}
    
    return -1;
}





int VertexSet::compute_downstairs_crit_midpts(const WitnessSet & W,
                                               vec_mp crit_downstairs,
                                               vec_mp midpoints_downstairs,
                                               std::vector< int > & index_tracker,
                                               vec_mp pi,
											  tracker_config_t * T)
{
    
    
	
	int retVal = SUCCESSFUL;
	
    int proj_index = this->get_proj_index(pi);
    
	vec_mp projection_values; init_vec_mp2(projection_values, int(W.num_points()) ,1024);
	projection_values->size = W.num_points();
	
	for (unsigned int ii=0; ii<W.num_points(); ii++){
		vec_mp & curr_point = W.point(ii);
		
		for (int jj=0; jj<curr_point->size; jj++) {
			
			if (!(mpfr_number_p(curr_point->coord[jj].r) && mpfr_number_p(curr_point->coord[jj].i))) {
				std::cout << color::red();
				std::cout << "there was NAN in a coordinate for a point in the projections to sort :(" << std::endl;
				print_point_to_screen_matlab(curr_point, "bad_point");
				
				print_point_to_screen_matlab(pi,"pi");
				std::cout << color::console_default();
				return CRITICAL_FAILURE;
			}
			
		}
		
        
        int curr_index = search_for_point(curr_point);
        
        if (curr_index < 0) {
            std::cout << color::red() << "trying to retrieve projection value from a non-stored point" << color::console_default() << std::endl;
			print_point_to_screen_matlab(curr_point,"curr_point");
			
			print_to_screen();
            mypause();
        }
        
		
        
        set_mp(&projection_values->coord[ii], & (vertices_[curr_index].projection_values())->coord[proj_index]);
        
		
		
		if (!(mpfr_number_p(projection_values->coord[ii].r) && mpfr_number_p(projection_values->coord[ii].i)))
		{
			print_comp_matlab(&projection_values->coord[ii],"not_a_number");
			print_point_to_screen_matlab(pi,"pi");
		}
		
	}
	

	real_threshold(projection_values, T->real_threshold);
	
	
	change_size_vec_mp(crit_downstairs,1); // destructive resize
	crit_downstairs->size = 1;
	
	retVal = sort_increasing_by_real(crit_downstairs, index_tracker, projection_values, 1e-30);
	
	clear_vec_mp(projection_values); // done with this data.  clear it.
	
	int num_midpoints = crit_downstairs->size - 1;
	
	if (num_midpoints<1) {
        // work is done, simply return
		return retVal;
	}
	
	
	comp_mp half; init_mp2(half,1024);
	mpf_set_str( half->r, "0.5", 10);
	mpf_set_str( half->i, "0.0", 10);

	change_size_vec_mp(midpoints_downstairs, num_midpoints);
	midpoints_downstairs->size = num_midpoints;
	
	comp_mp temp; init_mp2(temp,1024);
	for (int ii=0; ii<num_midpoints; ii++){
		add_mp(temp, &crit_downstairs->coord[ii], &crit_downstairs->coord[ii+1]);
		mul_mp(&midpoints_downstairs->coord[ii], temp, half);
	}
	
	real_threshold(midpoints_downstairs, T->real_threshold);
	
	clear_mp(temp);
	clear_mp(half);
	
	
	return retVal;
}




std::vector<int> VertexSet::assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value)
{
    if (this->curr_projection_<0) {
        std::cout << color::red() << "trying to assert projection value (current index) without having index set" << color::console_default() << std::endl;
        br_exit(-91621);
    }
    return assert_projection_value(relevant_indices,new_value,this->curr_projection_);
}

std::vector<int> VertexSet::assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index)
{
	
    
    if ( proj_index > num_projections_ ) {
		throw std::out_of_range("trying to assert projection value, but index of projection is larger than possible");
    }
	
	std::vector<int> bad_indices;
	
	comp_mp temp; init_mp(temp);
	
	
    for (std::set<int>::iterator ii=relevant_indices.begin(); ii!=relevant_indices.end(); ii++) {
        //*ii
        
        sub_mp(temp, &(vertices_[*ii].projection_values())->coord[proj_index], new_value);
        if (fabs(mpf_get_d(temp->r))>0.0001) {
            std::cout << "trying to assert projection value of " << mpf_get_d(new_value->r)
				      << " but original value is " << mpf_get_d((vertices_[*ii].projection_values())->coord[proj_index].r) << std::endl;
            std::cout << "point index is " << *ii << std::endl;
			bad_indices.push_back(*ii);
        }
        
        set_mp(&(vertices_[*ii].projection_values())->coord[proj_index], new_value);
    }
    
    
    clear_mp(temp);
    return bad_indices;
}





int VertexSet::add_vertex(const Vertex & source_vertex)
{
	
	
//	if (this->num_projections <= 0) {
//		std::cout << color::red() << "Vertex set has no projections!!!" << color::console_default() << std::endl;
//		br_exit(-9644);
//	}
//	
//	if (curr_input_index<0 && source_vertex.input_filename_index == -1) {
//		std::cout << color::red() << "adding points to Vertex set, but input file index unset." << color::console_default() << std::endl;
//		br_exit(6711);
//	}
	
	
	vertices_.push_back(source_vertex);
	
	
	if ((vertices_[num_vertices_].projection_values())->size < num_projections_) {
		increase_size_vec_mp((vertices_[num_vertices_].projection_values()), num_projections_ );
		(vertices_[num_vertices_].projection_values())->size = num_projections_;
	}
	
	for (int ii=0; ii<this->num_projections_; ii++){
		bool compute_proj_val = false;
		
		if (! (mpfr_number_p( (vertices_[num_vertices_].projection_values())->coord[ii].r) && mpfr_number_p( (vertices_[num_vertices_].projection_values())->coord[ii].i)  ) )
		{
			compute_proj_val = true;
		}
		//		else if ( false )//yeah, i dunno what else right yet.
		//		{
		//			compute_proj_val = true;
		//		}
		
		
		if (compute_proj_val==true) {
			projection_value_homogeneous_input(&(vertices_[num_vertices_].projection_values())->coord[ii],
											   vertices_[num_vertices_].point(),
											   projections_[ii]);
            // real_threshold(&(vertices_[num_vertices_].projection_values())->coord[ii],1e-13);
		}
		
	}
	
    
	
	if (vertices_[num_vertices_].input_filename_index() == -1)
	{
		vertices_[num_vertices_].set_input_filename_index(curr_input_index_);
	}
		
	
	this->num_vertices_++;
	return this->num_vertices_-1;
}


void VertexSet::print_to_screen() const
{
	
	vec_mp temp; init_vec_mp(temp,0); temp->size = 0;
	printf("Vertex set has %zu vertices:\n\n",num_vertices_);
	for (unsigned int ii=0; ii<this->num_vertices_; ++ii) {
		std::stringstream ss;  ss << "vertex_" << ii;
		dehomogenize(&temp, vertices_[ii].point(), num_natural_variables());
		temp->size = num_natural_variables()-1;
		
		print_point_to_screen_matlab(temp,ss.str());
		print_point_to_screen_matlab(vertices_[ii].get_projection_values(),"projection_values");
		printf("type: %d\n", vertices_[ii].type());
	}
}


int VertexSet::setup_vertices(boost::filesystem::path INfile)
//setup the Vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	unsigned int temp_num_vertices;
	int num_vars;
	int tmp_num_projections;
	int tmp_num_filenames;
	fscanf(IN, "%u %d %d %d\n\n", &temp_num_vertices, &tmp_num_projections, &num_natural_variables_, &tmp_num_filenames);
	
	
	vec_mp temp_vec; init_vec_mp2(temp_vec,num_natural_variables_,1024);
	for (int ii=0; ii<tmp_num_projections; ii++) {
		for (int jj=0; jj<num_natural_variables_; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
		}
		add_projection(temp_vec);
	}
	clear_vec_mp(temp_vec);
	
	scanRestOfLine(IN);
	scanRestOfLine(IN);
	
	
	for (int ii=0; ii<tmp_num_filenames; ii++) {
		int tmp_size;
		fscanf(IN,"%d\n",&tmp_size);
		
		char * buffer = new char[tmp_size];
		fgets(buffer, tmp_size, IN);
		boost::filesystem::path temppath = buffer;
		this->filenames_.push_back(temppath);
		
		delete [] buffer;
	}
	
	
	
	Vertex temp_vertex;
	
	for (unsigned int ii=0; ii<temp_num_vertices; ii++)
	{
		fscanf(IN, "%d\n", &num_vars);
		if ((temp_vertex.point())->size != num_vars) {
			change_size_vec_mp(temp_vertex.point(),num_vars); (temp_vertex.point())->size = num_vars;
		}
		
		for (int jj=0; jj<num_vars; jj++)
		{
			mpf_inp_str((temp_vertex.point())->coord[jj].r, IN, 10);
			mpf_inp_str((temp_vertex.point())->coord[jj].i, IN, 10);
		}
		
		int temp_num;
		fscanf(IN,"%d\n",&temp_num);
		increase_size_vec_mp(temp_vertex.projection_values(),temp_num);
		(temp_vertex.projection_values())->size = temp_num;
		for (int jj=0; jj<temp_num; jj++) {
			mpf_inp_str((temp_vertex.projection_values())->coord[jj].r, IN, 10);
			mpf_inp_str((temp_vertex.projection_values())->coord[jj].i, IN, 10);
		}
		
		int temp_int;
		fscanf(IN,"%d\n",&temp_int);
		temp_vertex.set_input_filename_index(temp_int);
		
		fscanf(IN,"%d\n",&temp_int);
	   temp_vertex.set_type(static_cast<VertexType>(temp_int));
		
		VertexSet::add_vertex(temp_vertex);
	}
	
	
	
	fclose(IN);
	
	if (this->num_vertices_!=temp_num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %zu\tactual: %u\n",num_vertices_,temp_num_vertices); // this is totally impossible.
		br_exit(25943);
	}
	
	return num_vertices_;
}







void VertexSet::print(boost::filesystem::path const& outputfile) const
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%zu %d %d %lu\n\n",num_vertices_,num_projections_, num_natural_variables_, filenames_.size());
	
	
	
	for (int ii=0; ii<num_projections_; ii++) {
		for (int jj=0; jj<num_natural_variables_; jj++) {
			print_mp(OUT, 0, &projections_[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	
	for (unsigned int ii=0; ii!=filenames_.size(); ii++) {
		int strleng = filenames_[ii].string().size() + 1; // +1 for the null character
		char * buffer = new char[strleng];
		memcpy(buffer, filenames_[ii].c_str(), strleng);
		fprintf(OUT,"%d\n",strleng);
		fprintf(OUT,"%s\n",buffer);
		delete [] buffer;
	}
	
	for (unsigned int ii = 0; ii < num_vertices_; ii++)
	{ // output points
		fprintf(OUT,"%d\n", (vertices_[ii].get_point())->size);
		for(int jj=0;jj<(vertices_[ii].get_point())->size;jj++) {
			print_mp(OUT, 0, &(vertices_[ii].get_point())->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",(vertices_[ii].get_projection_values())->size);
		for(int jj=0;jj<(vertices_[ii].get_projection_values())->size;jj++) {
			print_mp(OUT, 0, &(vertices_[ii].get_projection_values())->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",vertices_[ii].input_filename_index());
		
		fprintf(OUT,"\n");
		fprintf(OUT,"%d\n\n",vertices_[ii].type());
	}
	
	
	
	fclose(OUT);
	
}



int VertexSet::set_curr_input(boost::filesystem::path const& el_nom)
{
	
	int nom_index = -1;
    
	if ( el_nom.string().compare("unset_filename")==0 )
	{
		throw std::logic_error("trying to set curr_input from unset_filename");
	}
	
	
	int counter = 0;
	for (std::vector<boost::filesystem::path>::iterator ii = filenames_.begin(); ii!= filenames_.end(); ++ii) {
		
		if (*ii == el_nom) {
			nom_index = counter;
			break;
		}
		counter++;
	}
	
	if (nom_index==-1) {
		filenames_.push_back(el_nom);
		nom_index = counter;
	}
	
	curr_input_index_ = nom_index;
	return nom_index;
}

int VertexSet::get_proj_index(vec_mp proj) const
{
    int init_size = proj->size;
    proj->size = this->num_natural_variables_;
    
    
    int proj_index = -1;
    
    for (int ii=0; ii<num_projections_; ii++) {
		if (isSamePoint_inhomogeneous_input(projections_[ii],proj,same_point_tolerance_)) {
			proj_index = ii;
			break;
		}
	}
    
    proj->size = init_size;
    
    return proj_index;
}

int VertexSet::set_curr_projection(vec_mp new_proj){
    
    int proj_index = get_proj_index(new_proj);
    
    if (proj_index==-1) {
        int init_size = new_proj->size;
        new_proj->size = this->num_natural_variables_;
        
		proj_index = add_projection(new_proj);
        
        new_proj->size = init_size;
	}
	
    curr_projection_ = proj_index;
    
//        std::cout << "curr_projection is now: " << curr_projection << std::endl;
	return proj_index;
}



int VertexSet::add_projection(vec_mp proj){
	
	if (this->num_projections_==0) {
		projections_ = (vec_mp *) br_malloc(sizeof(vec_mp));
	}
	else{
		this->projections_ = (vec_mp *)br_realloc(this->projections_, (this->num_projections_+1) * sizeof(vec_mp));
	}
	
	
	init_vec_mp2(this->projections_[num_projections_],num_natural_variables_,T_->AMP_max_prec);
	this->projections_[num_projections_]->size = num_natural_variables_;
	
	
	if (proj->size != num_natural_variables_) {
		vec_mp tempvec;
		init_vec_mp2(tempvec,num_natural_variables_,T_->AMP_max_prec); tempvec->size = num_natural_variables_;
		for (int kk=0; kk<num_natural_variables_; kk++) {
			set_mp(&tempvec->coord[kk], &proj->coord[kk]);
		}
		vec_cp_mp(projections_[num_projections_], tempvec);
		clear_vec_mp(tempvec);
	}
	else
	{
		vec_cp_mp(projections_[num_projections_], proj);
	}
	

	num_projections_++;
	
	return num_projections_;
}


void VertexSet::reset()
{
	for (int ii=0; ii<num_projections_; ii++) {
		clear_vec_mp(projections_[ii]);
	}
	num_projections_ = 0;
	
	filenames_.resize(0);
	
	vertices_.resize(0);
	num_vertices_ = 0;
	
	clear();
	init();
}


void VertexSet::init()
{
	T_ = NULL;
	
	same_point_tolerance_ = 1e-8;
	num_projections_ = 0;
	projections_ = NULL;
	curr_projection_ = -1;
	
	curr_input_index_ = -2;
	
	this->num_vertices_ = 0;
	this->num_natural_variables_ = 0;
	
	init_vec_mp(checker_1_,0);
	init_vec_mp(checker_2_,0);
	
	
	
	init_mp(this->diff_);

	mpf_init(abs_);
	mpf_init(zerothresh_);
	mpf_set_d(zerothresh_, 1e-8);
}



void VertexSet::copy(const VertexSet &other)
{
	set_tracker_config(other.T());
	
	
	this->curr_projection_ = other.curr_projection_;
	if (this->num_projections_==0 && other.num_projections_>0) {
		projections_ = (vec_mp *) br_malloc(other.num_projections_*sizeof(vec_mp));
	}
	else if(this->num_projections_>0 && other.num_projections_>0) {
		projections_ = (vec_mp *) br_realloc(projections_,other.num_projections_*sizeof(vec_mp));
	}
	else if (this->num_projections_>0 && other.num_projections_==0){
		for (int ii=0; ii<this->num_projections_; ii++) {
			clear_vec_mp(projections_[ii]);
		}
		free(projections_);
	}
	
	for (int ii=0; ii<other.num_projections_; ii++) {
		if (ii>=this->num_projections_){
			init_vec_mp2(projections_[ii],1,1024); projections_[ii]->size = 1;
		}
		vec_cp_mp(projections_[ii],other.projections_[ii]);
	}
	
	this->num_projections_ = other.num_projections_;
	
	
	curr_input_index_ = other.curr_input_index_;
	filenames_ = other.filenames_;
	
	
	this->num_vertices_ = other.num_vertices_;
	this->num_natural_variables_ = other.num_natural_variables_;
	
	this->vertices_ = other.vertices_;
	
	vec_cp_mp(this->checker_1_,other.checker_1_);
	vec_cp_mp(this->checker_2_,other.checker_2_);
	
}


void VertexSet::clear()
{
	clear_vec_mp(checker_1_);
	clear_vec_mp(checker_2_);
	
	clear_mp(diff_);
	
	mpf_clear(abs_);
	mpf_clear(zerothresh_);
	
	for (int ii=0; ii<num_projections_; ii++) {
		clear_vec_mp(projections_[ii]);
	}
	free(projections_);
}
void VertexSet::send(int target, ParallelismConfig & mpi_config) const
{
	

	int num_filenames = filenames_.size();
	
	
	int * buffer2 = new int[6];
	buffer2[0] = num_natural_variables_;
	buffer2[1] = num_projections_;
	buffer2[2] = curr_projection_;
	buffer2[3] = num_filenames;
	buffer2[4] = curr_input_index_;
	buffer2[5] = num_vertices_;
	
	MPI_Send(buffer2, 6, MPI_INT, target, VERTEX_SET,  mpi_config.comm());
	
	delete [] buffer2;
	
	auto send_me = same_point_tolerance_;
	MPI_Send(&send_me,1,MPI_DOUBLE,target,VERTEX_SET, mpi_config.comm());
//	std::cout << "sending " << num_projections << " projections" << std::endl;
	
	buffer2 = new int[num_projections_];
	for (int ii=0; ii<num_projections_; ii++) {
		buffer2[ii] = projections_[ii]->size;
	}
	MPI_Send(buffer2, num_projections_, MPI_INT, target, VERTEX_SET,  mpi_config.comm());
	delete [] buffer2;
	
	for (int ii=0; ii<num_projections_; ii++) {
		send_vec_mp(projections_[ii],target);
	}

//	std::cout << "sending " << num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<num_filenames; ii++) {
		char * buffer;
		
		
		int strleng = filenames_[ii].string().size()+1;
		buffer = new char[strleng];
		memcpy(buffer, filenames_[ii].string().c_str(), strleng-1); // this sucks
		buffer[strleng-1] = '\0';
		
		MPI_Send(&strleng, 1, MPI_INT, target, VERTEX_SET,  mpi_config.comm());
//		std::cout << "sending filename length " << strleng << " " << filenames[ii].string() << std::endl;
		MPI_Send(&buffer[0], strleng, MPI_CHAR, target, VERTEX_SET,  mpi_config.comm());
		
		delete [] buffer;
		
	}
	
	
	for (unsigned int ii=0; ii<num_vertices_; ii++) {
		GetVertex(ii).send(target, mpi_config);
	}
	
	
	

	
	
	
	
	return;
}


void VertexSet::receive(int source, ParallelismConfig & mpi_config)
{
	MPI_Status statty_mc_gatty;
	
	int * buffer2 = new int[6];
	MPI_Recv(buffer2, 6, MPI_INT, source, VERTEX_SET,  mpi_config.comm(), &statty_mc_gatty);
	
	
	
	int temp_num_natural_variables = buffer2[0];
	int temp_num_projections = buffer2[1];
	curr_projection_ = buffer2[2];
	int temp_num_filenames = buffer2[3];
	curr_input_index_ = buffer2[4];
	unsigned int temp_num_vertices = buffer2[5];
	
	delete [] buffer2;
	
	MPI_Recv(&same_point_tolerance_,1,MPI_DOUBLE,source,VERTEX_SET, mpi_config.comm(), &statty_mc_gatty);
	
	set_num_vars(temp_num_natural_variables);
	
	
//	std::cout << "receiving " << temp_num_projections << " projections" << std::endl;
	
	buffer2 = new int[temp_num_projections];
	
	MPI_Recv(buffer2, temp_num_projections, MPI_INT, source, VERTEX_SET,  mpi_config.comm(), &statty_mc_gatty);

	
	
	vec_mp tempvec; init_vec_mp2(tempvec, 0, 1024);
	for (int ii=0; ii<temp_num_projections; ii++) {
//		std::cout << "recving " << ii << "th proj" << std::endl;
		change_size_vec_mp(tempvec,buffer2[ii]); tempvec->size = buffer2[ii];
		receive_vec_mp(tempvec,source);
		add_projection(tempvec);
//		print_point_to_screen_matlab(tempvec,"tempvec_recvd_proj");
		
	}
	clear_vec_mp(tempvec);
	delete [] buffer2;
	
	if (num_projections_!=temp_num_projections) {
		std::cout << "num_projections doesn't match!" << std::endl;
	}
	
//	std::cout << "receiving " << temp_num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<temp_num_filenames; ii++) {
		char * buffer; int strleng;
		
		MPI_Recv(&strleng, 1, MPI_INT, source, VERTEX_SET,  mpi_config.comm(), &statty_mc_gatty);
		
		buffer = new char[strleng];
//		std::cout << "recving filename length " << strleng << std::endl;
		MPI_Recv(&buffer[0], strleng, MPI_CHAR, source, VERTEX_SET,  mpi_config.comm(), &statty_mc_gatty);
		filenames_.push_back(boost::filesystem::path(std::string(buffer)));
		
		delete [] buffer;
		
	}
	
	
	
	

	
	for (unsigned int ii=0; ii<temp_num_vertices; ii++) {
		Vertex tempvert;
		tempvert.receive(source, mpi_config);
		add_vertex(tempvert);
	}
	
	if (num_vertices_ != temp_num_vertices) {
		std::cout << "logical inconsistency.  do not have correct num vertices." << std::endl;
	}
	
	
	

	
	return;
}




















