#include "nag/solvers/solver.hpp"


extern _comp_d  **mem_d;
extern _comp_mp **mem_mp;
extern int *size_d;  // size of mem_d
extern int *size_mp;  // size of mem_mp
extern int *mem_needs_init_d; // determine if mem_d has been initialized
extern int *mem_needs_init_mp; // determine if mem_mp has been initialized





void StraightLineProgramGlobalPointers::capture_globals()
{
	
	local_mem_d = mem_d;
	local_mem_mp = mem_mp;
	local_size_d = size_d;  // size of mem_d
	local_size_mp = size_mp;  // size of mem_mp
	local_mem_needs_init_d = mem_needs_init_d; // determine if mem_d has been initialized
	local_mem_needs_init_mp = mem_needs_init_mp; // determine if mem_mp has been initialized
	
//	if (local_mem_d==NULL) {
//		std::cout << "local_mem_d is NULL" << std::endl;
//	}
//	if (local_mem_mp==NULL) {
//		std::cout << "local_mem_mp is NULL" << std::endl;
//	}
//	
//	if (local_size_d==NULL) {
//		std::cout << "local_size_d is NULL" << std::endl;
//	}
//	
//	if (local_size_mp==NULL) {
//		std::cout << "local_size_mp is NULL" << std::endl;
//	}
//	
//	if (local_mem_needs_init_d==NULL) {
//		std::cout << "local_mem_needs_init_d is NULL" << std::endl;
//	}
//	if (local_mem_needs_init_mp==NULL) {
//		std::cout << "local_mem_needs_init_mp is NULL" << std::endl;
//	}
	
}

void StraightLineProgramGlobalPointers::set_globals_to_this()
{
	mem_d							= local_mem_d;
	mem_mp						= local_mem_mp;
	size_d						= local_size_d;  // size of mem_d
	size_mp						= local_size_mp;  // size of mem_mp
	mem_needs_init_d	= local_mem_needs_init_d; // determine if mem_d has been initialized
	mem_needs_init_mp = local_mem_needs_init_mp; // determine if mem_mp has been initialized
}









void adjust_tracker_AMP(tracker_config_t * T, int num_variables)
{
    
	
    T->AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T->AMP_Phi = T->AMP_bound_on_degree*(T->AMP_bound_on_degree-1.0)*T->AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
    T->AMP_Psi = T->AMP_bound_on_degree*T->AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
																			 // initialize latest_newton_residual_mp to the maximum precision
    
}




void SolverConfiguration::init()
{
	total_num_paths_tracked = 0;
	
	robust = false;
	
	PPD.num_funcs = 0;
	PPD.num_hom_var_gp = 0;
	PPD.num_var_gp = 0;
	PPD.size = NULL;
	PPD.type = NULL;
	

	path_number_modulus = 20;
	
	verbose_level_ = 0;
	
	
	use_midpoint_checker = 0;
	
	midpoint_tol = 1e-6;
	
	use_gamma_trick = 0;
	
	
}


std::vector<int> SolverOutput::get_cyclenums_noninfinite_w_mult()
{
	std::vector<int> cycle_nums;

	for (auto index = ordering.begin(); index != ordering.end(); ++index) {
		//index->second is the input index.  index->first is the index in vertices.  sorted by input index.
		if (metadata[index->first].is_finite)
			cycle_nums.push_back(metadata[index->first].CycleNumber());
	}
	return cycle_nums;
}


void SolverOutput::get_noninfinite_w_mult_full(WitnessSet & W_transfer)
{
	get_noninfinite_w_mult(W_transfer);
	
	get_patches_linears(W_transfer);
	
	set_witness_set_nvars(W_transfer);
	
}


void SolverOutput::get_noninfinite_w_mult(WitnessSet & W_transfer)
{
	for (auto index = ordering.begin(); index != ordering.end(); ++index) {
		//index->second is the input index.  index->first is the index in vertices.  sorted by input index.
		if (metadata[index->first].is_finite) {
			W_transfer.add_point(vertices_[index->first].point());
		}
	}
	
	set_witness_set_nvars(W_transfer);
}


void SolverOutput::get_nonsing_finite_multone(WitnessSet & W_transfer)
{
	for (auto index = ordering.begin(); index != ordering.end(); ++index) {
		//index->second is the input index.  index->first is the index in vertices.  sorted by input index.
		if ( (metadata[index->first].is_finite) && (!metadata[index->first].is_singular) && (metadata[index->first].multiplicity==1) ) {
			W_transfer.add_point(vertices_[index->first].point());
		}
	}
	
	set_witness_set_nvars(W_transfer);
}

void SolverOutput::get_multpos(std::map<int, WitnessSet> & W_transfer)
{
	
	//for each multiplicity, construct the witness_sets
	for (auto mult_ind = occuring_multiplicities.begin(); mult_ind!=occuring_multiplicities.end(); ++mult_ind) {
		
		int num_added_points = 0;
		for (auto index = metadata.begin(); index != metadata.end(); ++index) {
			if ((index->multiplicity== *mult_ind) && (index->is_finite))  {
				W_transfer[*mult_ind].add_point(vertices_[index->output_index].point());
				num_added_points++;
			}
		}
		
		if (num_added_points>0) { // avoids adding an empty witness set when all the points were infinite.
			set_witness_set_nvars(W_transfer[*mult_ind]);
		}
		
	}
	
	
}

void SolverOutput::get_multpos_full(std::map<int, WitnessSet> & W_transfer)
{
	
	get_multpos(W_transfer);
	
	for (auto iter = W_transfer.begin(); iter!=W_transfer.end(); ++iter) {
		get_patches_linears(iter->second);
	}
	
	
}



void SolverOutput::get_sing(WitnessSet & W_transfer)
{
	for (auto index = ordering.begin(); index != ordering.end(); ++index) {
		//index->second is the input index.  index->first is the index in vertices.  sorted by input index.
		if ( (metadata[index->first].is_singular) ) {
			W_transfer.add_point(vertices_[index->first].point());
		}
	}
	set_witness_set_nvars(W_transfer);
}


void SolverOutput::get_sing_finite(WitnessSet & W_transfer)
{
	for (auto index = ordering.begin(); index != ordering.end(); ++index) {
		//index->second is the input index.  index->first is the index in vertices.  sorted by input index.
		if ( (metadata[index->first].is_singular && metadata[index->first].is_finite) ) {
			W_transfer.add_point(vertices_[index->first].point());
		}
	}
	set_witness_set_nvars(W_transfer);
}













int Solver::send(ParallelismConfig & mpi_config)
{
	
	
	int *buffer = new int[4];
	
	buffer[0] = this->num_variables;
	buffer[1] = this->num_steps;
	buffer[2] = this->verbose_level_;
	buffer[3] = this->MPType;
	
	MPI_Bcast(buffer, 4, MPI_INT, 0, mpi_config.comm());
	
	send_preproc_data(&this->preProcData);
	
	
	delete[] buffer;
	
	
	return SUCCESSFUL;
}

int Solver::receive(ParallelismConfig & mpi_config)
{
	
    
	int *buffer = new int[4];
	
    
	MPI_Bcast(buffer, 4, MPI_INT, 0, mpi_config.comm());
	
	this->num_variables = buffer[0];
	this->num_steps = buffer[1];
	this->verbose_level(buffer[2]);
	this->MPType = buffer[3];
	
	receive_preproc_data(&this->preProcData);
	have_PPD = true;
	
	
	
	
	
	delete[] buffer;
	return SUCCESSFUL;
}













//////////////
//
//	SOLVER MP
//
////////////


void SolverMultiplePrecision::clear()
{

	patch_eval_data_clear_mp(&this->patch);
	
	
	
	clear_mp(gamma);
	
	if (MPType==2) {
		clear_rat(gamma_rat);
		free(gamma_rat);
	}
	
	
	if (have_SLP && received_mpi) { // other wise don't have it, or someone else is responsible for clearing it.
		clearProg(this->SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete[] SLP;
	}
	
	if (received_mpi) {
		;
	}
}



int SolverMultiplePrecision::send(ParallelismConfig & mpi_config)
{
	Solver::send(mpi_config);
	
    
    int num_SLP;
    if (this->have_SLP)
        num_SLP = 1;
    else
        num_SLP = 0;
    
    MPI_Bcast(&num_SLP, 1, MPI_INT, 0, mpi_config.comm());
    
    for (int ii=0; ii<num_SLP; ii++) {
        bcast_prog_t(this->SLP, MPType, 0, 0);
    }
    
    
    
	send_patch_mp(&this->patch);
	
	
	bcast_comp_mp(this->gamma, 0,0);
	if (this->MPType==2) {
        bcast_comp_rat(gamma_rat, 0, 0);
	}
	else{
		
	}
	
	
	int * buffer;
	buffer = new int[1];
	buffer[0] = this->curr_prec;
	MPI_Bcast(buffer, 1, MPI_INT, 0, mpi_config.comm());
	delete[] buffer;
	
	
	
	randomizer_->bcast_send(mpi_config);
	
	return SUCCESSFUL;
}


int SolverMultiplePrecision::receive(ParallelismConfig & mpi_config)
{
	
	Solver::receive(mpi_config);
	
	
    if (this->have_SLP) {
        clearProg(this->SLP, this->MPType, 1);
        this->have_SLP = false;
    }
    
    received_mpi = true;
    int num_SLP;
    MPI_Bcast(&num_SLP, 1, MPI_INT, 0, mpi_config.comm()); // get the number of SLP's to receieve
    
    if (num_SLP>0) {
        prog_t * _SLP = new prog_t[num_SLP];//(prog_t *) br_malloc(num_SLP*sizeof(prog_t));
        for (int ii=0; ii<num_SLP; ii++) {
            //			std::cout << "worker bcasting the SLP, MPType" << this->MPType << std::endl;
            bcast_prog_t(&_SLP[ii], this->MPType, 1, 0); // last two arguments are: myid, headnode
														 //			std::cout << "worker copying the SLP" << std::endl;
            this->SLP = &_SLP[ii];
            //			cp_prog_t(this->SLP, &_SLP[ii]);
            //			std::cout << "worker copied the SLP" << std::endl;
        }
        
        this->have_SLP = true;
        initEvalProg(this->MPType);
		SLP_memory.capture_globals();
    }
    
	receive_patch_mp(&this->patch);
	
	bcast_comp_mp(this->gamma, 1,0);
	
	if (this->MPType==2) {
        bcast_comp_rat(gamma_rat, 1, 0);
	}
	else{

	}
	

	
	
	int *buffer = new int[1];
	
	MPI_Bcast(buffer, 1, MPI_INT, 0, mpi_config.comm());
	this->curr_prec = buffer[0];
	
	delete[] buffer;
	
	
	randomizer_ = std::make_shared<SystemRandomizer> (*(new SystemRandomizer));
	randomizer_->bcast_receive(mpi_config);
	
	return SUCCESSFUL;
}







//////////////
//
//	SOLVER D
//
////////////



void SolverDoublePrecision::clear(){
	
	if (MPType==0) {
	}

	//	std::cout << "clearing double patch" << std::endl;
	patch_eval_data_clear_d(& this->patch);
	
	clear_d(gamma);
	
	
	if (have_SLP && received_mpi) {
		clearProg(this->SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete[] SLP;
	}
}





int SolverDoublePrecision::send(ParallelismConfig & mpi_config)
{
	Solver::send(mpi_config);
	
    if (this->MPType == 0) {
        int num_SLP;
        if (this->have_SLP)
            num_SLP = 1;
        else
            num_SLP = 0;
        
        MPI_Bcast(&num_SLP, 1, MPI_INT, 0, mpi_config.comm());
        
        for (int ii=0; ii<num_SLP; ii++) {
            //		std::cout << "master bcasting the SLP, MPType" << this->MPType << std::endl;
            bcast_prog_t(this->SLP, MPType, 0, 0);
        }
    }
    
    
    
    
	send_patch_d(&(this->patch));
    
	bcast_comp_d(this->gamma, 0,0);
		
	
	if (MPType==0) {
		randomizer_->bcast_send(mpi_config);
	}

	
	
	return SUCCESSFUL;
}


int SolverDoublePrecision::receive(ParallelismConfig & mpi_config)
{
    
	
	Solver::receive(mpi_config);
    
    if (this->MPType == 0) {
        received_mpi = true;
        if (this->have_SLP) {
            clearProg(this->SLP, this->MPType, 1);
            this->have_SLP = false;
        }
        
        
        int num_SLP;
        MPI_Bcast(&num_SLP, 1, MPI_INT, 0, mpi_config.comm()); // get the number of SLP's to receieve
        
        if (num_SLP>0) {
            prog_t * _SLP = new prog_t[num_SLP];
            for (int ii=0; ii<num_SLP; ii++) {
                bcast_prog_t(&_SLP[ii], this->MPType, 1, 0); // last two arguments are: myid, headnode
                this->SLP = &_SLP[ii];
				
            }
            
            this->have_SLP = true;
            initEvalProg(this->MPType);
			SLP_memory.capture_globals();
        }
    }
    else if (MPType==2) {
        this->setup(BED_mp->SLP, BED_mp->randomizer());
		if (have_SLP) {
			SLP_memory.capture_globals();
		}
		
    }
    
	
    
	receive_patch_d(&this->patch); // the receiving part of the broadcast
	
	bcast_comp_d(this->gamma, 1, 0);
	
	
	if (MPType==0) {
		randomizer_ = std::make_shared<SystemRandomizer>( *( new SystemRandomizer));
		randomizer_->bcast_receive(mpi_config);
	}
	else
	{
		randomizer_ = BED_mp->randomizer(); // i believe this was already accomplished
	}

	
	return SUCCESSFUL;
}








//////





void get_tracker_config(SolverConfiguration & solve_options,int MPType)
{
	
	//necessary for the setupConfig call
	double intrinsicCutoffMultiplier;
	int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, supersetOnly = 0, paramHom = 0;
	//end necessaries for the setupConfig call.
	int constructWitnessSet = 0;
	
    setupConfig(&solve_options.T,
                &solve_options.midpoint_tol,
                &userHom,
                &useRegen,
                &regenStartLevel,
                &maxCodim,
                &specificCodim,
                &solve_options.path_number_modulus,
                &intrinsicCutoffMultiplier,
                &reducedOnly,
                &constructWitnessSet,
                &supersetOnly,
                &paramHom,
                MPType);
	
	mpf_init2(solve_options.T.latest_newton_residual_mp, solve_options.T.AMP_max_prec);
	
	
	
	return;
}










void master_solver(SolverOutput & solve_out, const WitnessSet & W,
                   SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
                   SolverConfiguration & solve_options)
{
	
	
	initMP(mpf_get_default_prec());
	
    int num_crossings = 0;
	
	solve_out.num_variables = W.num_variables();
	solve_out.num_natural_vars = W.num_natural_variables();
	
	solve_out.copy_patches(W); // copy the patches over from the original witness set
	solve_out.copy_names(W);
	
	
	if (solve_options.use_parallel()) {
		
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
		
		
		int *settings_buffer = (int *) br_malloc(2*sizeof(int));
		settings_buffer[0] = solve_options.robust;
		settings_buffer[1] = solve_options.use_gamma_trick;
		
		MPI_Bcast(settings_buffer,2,MPI_INT, 0, solve_options.comm());
		free(settings_buffer);
		
		switch (solve_options.T.MPType) {
			case 1:
				
				ED_mp->send(solve_options);
                //				std::cout << "master done sending mp type" << std::endl;
				break;
				
			default:
				
				ED_d->send(solve_options);
                //				std::cout << "master done sending double type " << solve_options.T.MPType <<  std::endl;
				break;
		}
	}
	
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	post_process_t *endPoints = (post_process_t *)br_malloc(W.num_points() * sizeof(post_process_t)); //overallocate, expecting full
	
	
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "output",
						&midOUT, "midpath_data");
	
	if (solve_options.use_parallel()) {
		
		
		master_tracker_loop(&trackCount, OUT, midOUT,
                            W,
                            endPoints,
                            ED_d, ED_mp,
                            solve_options);
	}
	else{
		serial_tracker_loop(&trackCount, OUT, midOUT,
                             W,
                             endPoints,
                             ED_d, ED_mp,
                             solve_options);
	}
	
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	
	
	// check for path crossings
	if (solve_options.use_midpoint_checker==1) {
		midpoint_checker(trackCount.numPoints, solve_options.T.numVars,solve_options.midpoint_tol, &num_crossings);
	}
	
	
	// post process
	switch (solve_options.T.MPType) {
		case 0:
			solve_out.post_process(endPoints, trackCount.successes, &ED_d->preProcData, &solve_options.T, solve_options);
			break;
			
		default:
			solve_out.post_process(endPoints, trackCount.successes, &ED_mp->preProcData, &solve_options.T, solve_options);
			break;
	}
	
    //clear the endpoints here
	for (int ii=0; ii<trackCount.successes; ii++) {
		clear_post_process_t(&endPoints[ii],W.num_variables());
	}
	free(endPoints);
}






/**
 sets the start_pts structure to hold all points in W
 
 \param startPts the value being set.  should be NULL or uninitialized input.
 \param W the WitnessSet input
 
 */
void generic_set_start_pts(point_data_d ** startPts,
                           const WitnessSet & W)
{
	
	*startPts = (point_data_d *)br_malloc(W.num_points() * sizeof(point_data_d));
	
	for (unsigned int ii = 0; ii < W.num_points(); ii++)
	{ // setup startPts[ii]
		init_point_data_d(&(*startPts)[ii], W.num_variables()); // also performs initialization on the point inside startPts
		change_size_vec_d((*startPts)[ii].point,W.num_variables());
		(*startPts)[ii].point->size = W.num_variables();
		
		//1 set the coordinates
		vec_mp_to_d((*startPts)[ii].point, W.point(ii));
		
		//2 set the start time to 1.
		set_one_d((*startPts)[ii].time);
	}
}



void generic_set_start_pts(point_data_mp ** startPts,
                           const WitnessSet & W)
{
	
	(*startPts) = (point_data_mp *)br_malloc(W.num_points() * sizeof(point_data_mp));
	
	for (unsigned int ii = 0; ii < W.num_points(); ii++)
	{ // setup startPts[ii]
		init_point_data_mp(&(*startPts)[ii], W.num_variables()); // also performs initialization on the point inside startPts
		change_size_vec_mp((*startPts)[ii].point,W.num_variables());
		(*startPts)[ii].point->size = W.num_variables();
		
		//1 set the coordinates
		vec_cp_mp((*startPts)[ii].point, W.point(ii));
		
		//2 set the start time to 1.
		set_one_mp((*startPts)[ii].time);
	}
}






void serial_tracker_loop(trackingStats *trackCount,
                          FILE * OUT, FILE * MIDOUT,
                          const WitnessSet & W,  // was the startpts file pointer.
                          post_process_t *endPoints,
                          SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
                          SolverConfiguration & solve_options)
{
	
	
    
	
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
    int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	int (*change_prec)(void const *, int) = NULL;
	int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	bool randomized;
	switch (solve_options.T.MPType) {
		case 0:
			curr_eval_d = ED_d->evaluator_function_d;
			curr_eval_mp = ED_d->evaluator_function_mp;
			change_prec = ED_d->precision_changer;
			find_dehom = ED_d->dehomogenizer;
			randomized = ED_d->IsRandomized();
			break;
            
		default:
			curr_eval_d = ED_mp->evaluator_function_d;
			curr_eval_mp = ED_mp->evaluator_function_mp;
			change_prec = ED_mp->precision_changer;
			find_dehom = ED_mp->dehomogenizer;
			randomized = ED_mp->IsRandomized();
			break;
			
	}
    
	
	
	
	point_data_d *startPts_d = NULL;
	generic_set_start_pts(&startPts_d, W);
	
	point_data_mp *startPts_mp = NULL;
	generic_set_start_pts(&startPts_mp, W);
	
	
	
    // setup the rest of the structures
	endgame_data_t EG; //this will hold the temp solution data produced for each individual track
	init_endgame_data(&EG, solve_options.T.Precision);
	
	
	
	
	trackCount->numPoints = W.num_points();
	int solution_counter = 0;
	
	
	
	// track each of the start points
	
	for (unsigned int ii = 0; ii < W.num_points(); ii++)
	{
		if ((solve_options.verbose_level()>=0) && (solve_options.path_number_modulus!=0) )
		{
			if ((ii%solve_options.path_number_modulus)==0 && ( unsigned(solve_options.path_number_modulus)<W.num_points() )) {
				std::cout << color::gray();
				std::cout << "tracking path " << ii << " of " << W.num_points() << std::endl;
				std::cout << color::console_default();
			}
			
		}
		
		solve_options.increment_num_paths_tracked();
		
		
		
		//        print_point_to_screen_matlab(startPts_d[ii].point,"start");
        
		if (solve_options.robust==true) {
			robust_track_path(ii, &EG,
                              &startPts_d[ii], &startPts_mp[ii],
                              OUT, MIDOUT,
                              solve_options, ED_d, ED_mp,
                              curr_eval_d, curr_eval_mp, change_prec, find_dehom);
		}
		else{
//            boost::timer::auto_cpu_timer t;
            // track the path
			generic_track_path(ii, &EG,
                               &startPts_d[ii], &startPts_mp[ii],
                               OUT, MIDOUT,
                               &solve_options.T, ED_d, ED_mp,
                               curr_eval_d, curr_eval_mp, change_prec, find_dehom);
            
		}
		
		// check to see if it should be sharpened
		if (EG.retVal == 0 && solve_options.T.sharpenDigits > 0)
		{ // use the sharpener for after an endgame
			sharpen_endpoint_endgame(&EG, &solve_options.T, OUT, ED_d, ED_mp, curr_eval_d, curr_eval_mp, change_prec);
		}
		
		
		
		int issoln;
		
		if (randomized){
			switch (solve_options.T.MPType) {
				case 0:
	                issoln = ED_d->is_solution_checker_d(&EG,  &solve_options.T, ED_d);
	                
					break;
					
				default:
					
					if (EG.prec<64){
						issoln = ED_mp->is_solution_checker_d(&EG,  &solve_options.T, ED_d); }
					else {
						issoln = ED_mp->is_solution_checker_mp(&EG, &solve_options.T, ED_mp); }
					break;
			}
		}
		else
		{
			issoln = true;
		}
		
		
		//get the terminal time in double form
		comp_d time_to_compare;
		if (EG.prec < 64) {
			set_d(time_to_compare,EG.PD_d.time);}
		else {
			mp_to_d(time_to_compare, EG.PD_mp.time); }
		
		
		if ((EG.retVal != 0 && time_to_compare->r > solve_options.T.minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
			
			trackCount->failures++;
			
			if (solve_options.verbose_level()>=3) {
                
				printf("\nthere was a path failure tracking witness point %d\nretVal = %d; issoln = %d\n",ii,EG.retVal, issoln);
				
				print_path_retVal_message(EG.retVal);
				
				if (solve_options.verbose_level() >= 5) {
					if (EG.prec < 64)
						print_point_to_screen_matlab(EG.PD_d.point,"bad_terminal_point");
					else
						print_point_to_screen_matlab(EG.PD_mp.point,"bad_terminal_point");
				}
			}
		}
		else
		{
			//otherwise converged, but may have still had non-zero retval due to other reasons.
			endgamedata_to_endpoint(&endPoints[solution_counter], &EG);
			trackCount->successes++;
			solution_counter++; // probably this could be eliminated
		}
		
	}// re: for (ii=0; ii<W.num_points ;ii++)
	clear_endgame_data(&EG);
	
	
	//clear the data structures.
    for (unsigned int ii = 0; ii < W.num_points(); ii++)
    { // clear startPts[ii]
        clear_point_data_d(&startPts_d[ii]);
		clear_point_data_mp(&startPts_mp[ii]);
    }
    free(startPts_d);
	free(startPts_mp);
    
	
	
}



void master_tracker_loop(trackingStats *trackCount,
                         FILE * OUT, FILE * MIDOUT,
                         const WitnessSet & W,  // was the startpts file pointer.
                         post_process_t *endPoints,
                         SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
                         SolverConfiguration & solve_options)
{
	
	
	point_data_d *startPts_d = NULL;
	point_data_mp *startPts_mp = NULL;
	
	
	
	
	switch (solve_options.T.MPType) {
		case 1:
			generic_set_start_pts(&startPts_mp, W);
			break;
			
		default:
			generic_set_start_pts(&startPts_d, W);
			break;
	}
	
	
	
	
	


	
	
	trackCount->numPoints = W.num_points();
	int solution_counter = 0;
	
	int total_number_points = W.num_points();
	MPI_Bcast(&total_number_points, 1, MPI_INT, solve_options.head(), solve_options.comm());
	
	
	int max_incoming = get_num_at_a_time(solve_options.num_procs()-1,total_number_points);
	// setup the rest of the structures
	endgame_data_t * EG_receives = (endgame_data_t *) br_malloc(max_incoming*sizeof(endgame_data_t)); //this will hold the temp solution data produced for each individual track
	for (int ii=0; ii<max_incoming; ii++) {
		init_endgame_data(&EG_receives[ii], 64);
	}
	
	
	
	// seed the workers
	int next_index = 0;
	
	for (unsigned int ii=1; int(ii)<solve_options.num_procs() && next_index< int(W.num_points()); ii++) {
		int next_worker = solve_options.activate_next_worker();
        
		int num_packets = get_num_at_a_time(solve_options.num_procs()-1,total_number_points-next_index);
		send_start_points(next_worker, num_packets,
                          startPts_d,
                          startPts_mp,
                          next_index,// gets muted here
                          solve_options);
		
	}
	
	
    
	while (next_index<total_number_points)
	{
        
        
		receive_endpoints(trackCount,
                          &EG_receives, max_incoming,
						  solution_counter,
						  endPoints,
						  ED_d, ED_mp,
						  solve_options);
		
		
		
        int next_worker = solve_options.activate_next_worker();
        
		
        send_start_points(next_worker, get_num_at_a_time(solve_options.num_procs()-1,total_number_points-next_index),
                          startPts_d,
                          startPts_mp,
                          next_index,
                          solve_options);
		
	}// re: for (ii=0; ii<W.num_points ;ii++)
	
	while (solve_options.have_active()) {
//		std::cout << "waiting to receive from active worker" << std::endl;
		receive_endpoints(trackCount,
						  &EG_receives, max_incoming,
						  solution_counter,
						  endPoints,
						  ED_d, ED_mp,
						  solve_options);
	} // re: while have active
	
	
	solve_options.send_all_available(0);
	
	//clear the data structures.
	switch (solve_options.T.MPType) {
		case 1:
			for (unsigned int ii = 0; ii < W.num_points(); ii++)
				clear_point_data_mp(&startPts_mp[ii]);
			free(startPts_mp);
			break;
			
		default:
			for (unsigned int ii = 0; ii < W.num_points(); ii++)
				clear_point_data_d(&startPts_d[ii]);
			free(startPts_d);
			break;
	}
    
}


void worker_tracker_loop(trackingStats *trackCount,
                         FILE * OUT, FILE * MIDOUT,
                         SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
                         SolverConfiguration & solve_options)
{
	
	
	int total_number_points;
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
    int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	int (*change_prec)(void const *, int) = NULL;
	int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;

	switch (solve_options.T.MPType) {
		case 0:
			curr_eval_d = ED_d->evaluator_function_d;
			curr_eval_mp = ED_d->evaluator_function_mp;
			change_prec = ED_d->precision_changer;
			find_dehom = ED_d->dehomogenizer;
			break;
			
		default:
			curr_eval_d = ED_mp->evaluator_function_d;
			curr_eval_mp = ED_mp->evaluator_function_mp;
			change_prec = ED_mp->precision_changer;
			find_dehom = ED_mp->dehomogenizer;
			break;
			
	}
	
	
	
	
	point_data_d *startPts_d;
	point_data_mp *startPts_mp;
	
	switch (solve_options.T.MPType) {
		case 1:
			startPts_mp = (point_data_mp *) br_malloc(1*sizeof(point_data_mp));
			init_point_data_mp(&startPts_mp[0], ED_mp->num_variables);
			break;
			
		default:
			startPts_d = (point_data_d *) br_malloc(1*sizeof(point_data_d));
			init_point_data_d(&startPts_d[0], ED_d->num_variables);
			break;
	}
	
	
    // setup the rest of the structures
	endgame_data_t * EG = (endgame_data_t *) br_malloc(1*sizeof(endgame_data_t)); //this will hold the temp solution data produced for each individual track
	init_endgame_data(&EG[0], solve_options.T.Precision);
	
	int *indices_incoming = (int *) br_malloc(1*sizeof(int));
	
	MPI_Status statty_mc_gatty;
	int max_num_allocated = 1;
	
	int numStartPts = 1;
	
	MPI_Bcast(&total_number_points, 1, MPI_INT, solve_options.head(), solve_options.comm());
	while (1)
	{
        //		std::cout << "worker" << solve_options.id() << " receiving work" << std::endl;
		MPI_Recv(&numStartPts, 1, MPI_INT, solve_options.head(), NUMPACKETS, solve_options.comm(), &statty_mc_gatty);
		// recv next set of start points
		
		if (numStartPts==0) {
			break;
		}
		
		if (numStartPts>max_num_allocated) {
			switch (solve_options.T.MPType) {
					
				case 1:
					for (int zz=0; zz<max_num_allocated; zz++) {
						clear_point_data_mp(&startPts_mp[zz]);
					}
					startPts_mp = (point_data_mp *) br_realloc(startPts_mp, numStartPts*sizeof(point_data_mp));
					for (int zz=0; zz<numStartPts; zz++) {
						init_point_data_mp(&startPts_mp[zz], ED_mp->num_variables);
					}
					break;
					
				default:
					for (int zz=0; zz<max_num_allocated; zz++) {
						clear_point_data_d(&startPts_d[zz]);
					}
					startPts_d = (point_data_d *) br_realloc(startPts_d, numStartPts*sizeof(point_data_d));
					for (int zz=0; zz<numStartPts; zz++) {
						init_point_data_d(&startPts_d[zz], ED_d->num_variables);
					}
					break;
			}
			
			indices_incoming = (int *) br_realloc(indices_incoming, numStartPts*sizeof(int));
			
			for (int zz=0; zz<max_num_allocated; zz++) {
				clear_endgame_data(&EG[zz]);
			}
			EG = (endgame_data_t *) br_realloc(EG,numStartPts*sizeof(endgame_data_t));
			for (int zz=0; zz<numStartPts; zz++) {
				init_endgame_data(&EG[zz], solve_options.T.Precision);
			}
			max_num_allocated = numStartPts;
		}
		
		MPI_Recv(indices_incoming, numStartPts, MPI_INT, solve_options.head(), INDICES, solve_options.comm(), &statty_mc_gatty);
		
		for (int ii=0; ii<numStartPts; ii++) {
			switch (solve_options.T.MPType) {
				case 1:
					receive_vec_mp(startPts_mp[ii].point, solve_options.head());
					set_one_mp(startPts_mp[ii].time);
					break;
					
				default:
					receive_vec_d(startPts_d[ii].point, solve_options.head());
					set_one_d(startPts_d[ii].time);
					break;
			}
		}
        
		
		// track each of the start points
		for (int ii = 0; ii < numStartPts; ii++)
		{
			int current_index =  indices_incoming[ii];
			
			
			if ((solve_options.verbose_level()>=0) && (solve_options.path_number_modulus!=0) && (solve_options.path_number_modulus < total_number_points) )
			{
				if ((current_index%solve_options.path_number_modulus)==0) {
					std::cout << color::gray();
					printf("tracking path %d of %d\n",current_index, total_number_points);
					std::cout << color::console_default();
				}
				
			}
			
            
            
			
            if (solve_options.robust==true) {
				//                boost::timer::auto_cpu_timer t;
                robust_track_path(indices_incoming[ii], &EG[ii],
                                  &startPts_d[ii], &startPts_mp[ii],
                                  OUT, MIDOUT,
                                  solve_options, ED_d, ED_mp,
                                  curr_eval_d, curr_eval_mp, change_prec, find_dehom);
                
            }
            else{
                // track the path
                generic_track_path(indices_incoming[ii], &EG[ii],
                                   &startPts_d[ii], &startPts_mp[ii],
                                   OUT, MIDOUT,
                                   &solve_options.T, ED_d, ED_mp,
                                   curr_eval_d, curr_eval_mp, change_prec, find_dehom);
            }
			
            
			
			// check to see if it should be sharpened
			if (EG[ii].retVal == 0 && solve_options.T.sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[ii], &solve_options.T, OUT, ED_d, ED_mp, curr_eval_d, curr_eval_mp, change_prec);
			}
			
			
            
		}// re: for (ii=0; ii<W.num_points ;ii++)
		
		
		MPI_Send(&numStartPts, 1, MPI_INT, solve_options.head(), NUMPACKETS, solve_options.comm());
		
		send_recv_endgame_data_t(&EG, &numStartPts, solve_options.T.MPType, solve_options.head(), 1); //1 is 'issending'
		
		
	}
    
	
	switch (solve_options.T.MPType) {
		case 1:
			for (int ii=0; ii<max_num_allocated; ii++) {
				clear_point_data_mp(&startPts_mp[ii]);
				clear_endgame_data(&EG[ii]);
			}
			free(startPts_mp);
			break;
			
		default:
			for (int ii=0; ii<max_num_allocated; ii++) {
				clear_point_data_d(&startPts_d[ii]);
				clear_endgame_data(&EG[ii]);
			}
			free(startPts_d);
			break;
	}
	free(EG);
	free(indices_incoming);
}


int get_num_at_a_time(int num_workers, int num_points)
{
	int num_packets = 1 + ((num_points - 1) / num_workers);
	
	num_packets = 1 + ((num_packets - 1) / 25);
	
	
	return num_packets;
}

void send_start_points(int next_worker, int num_packets,
                       point_data_d *startPts_d,
                       point_data_mp *startPts_mp,
                       int & next_index,
                       SolverConfiguration & solve_options)
{
	MPI_Send(&num_packets, 1, MPI_INT, next_worker, NUMPACKETS, solve_options.comm());
	
	int *indices_outgoing = new int[num_packets];
	
	for (int ii=0; ii<num_packets; ii++) {
		indices_outgoing[ii] = next_index;
		next_index++;
	}
	
    
	MPI_Send(indices_outgoing, num_packets, MPI_INT, next_worker, INDICES, solve_options.comm());
	
	
	
	for (int ii=indices_outgoing[0]; ii<=indices_outgoing[num_packets-1]; ii++) {
		if (solve_options.T.MPType==1) {
			send_vec_mp( startPts_mp[ii].point, next_worker);
			
		}
		else
		{
			send_vec_d( startPts_d[ii].point, next_worker);
			
		}
	}
	
	delete[] indices_outgoing;
	return;
}


int receive_endpoints(trackingStats *trackCount,
                      endgame_data_t **EG_receives, int & max_incoming,
                      int & solution_counter,
                      post_process_t *endPoints,
                      SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
                      SolverConfiguration & solve_options)
{
	
	
	//now to receive data
	int num_incoming;
	MPI_Status statty_mc_gatty;

	MPI_Recv(&num_incoming, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, solve_options.comm(), &statty_mc_gatty);
	
	if (statty_mc_gatty.MPI_SOURCE>= int(solve_options.num_procs()) || statty_mc_gatty.MPI_SOURCE<0) {
		std::cout << statty_mc_gatty.MPI_SOURCE << ", " << solve_options.num_procs() << std::endl;
		mypause();
	}
	
	bool randomized;
	switch (solve_options.T.MPType) {
		case 0:
			randomized = ED_d->IsRandomized();
			break;
		default:
			randomized = ED_mp->IsRandomized();
			break;
	}
				
	
	if (num_incoming > max_incoming) {
		std::cout << "the impossible happened -- want to receive more endpoints than max" << std::endl;
	}
	

	int incoming_id = send_recv_endgame_data_t(EG_receives, &num_incoming, solve_options.T.MPType, statty_mc_gatty.MPI_SOURCE, 0); // the trailing 0 indicates receiving

	
	solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
	

	
	for (int ii=0; ii<num_incoming; ii++) {
		int issoln;

		if (randomized)
		{
			switch (solve_options.T.MPType) {
				case 0:
					issoln = ED_d->is_solution_checker_d( &(*EG_receives)[ii],  &solve_options.T, ED_d);
					
					break;
					
				default:
					
					if ((*EG_receives)[ii].prec<64){
						issoln = ED_mp->is_solution_checker_d( &((*EG_receives)[ii]), &solve_options.T, ED_d); } // this function call is a reference!
					else {
						issoln = ED_mp->is_solution_checker_mp( &((*EG_receives)[ii]), &solve_options.T, ED_mp); } // this function call is a reference!
					break;
			}
		}
		else
			issoln = true;
		
		
		//get the terminal time in double form
		comp_d time_to_compare;
		if ((*EG_receives)[ii].prec < 64) {
			set_d(time_to_compare,(*EG_receives)[ii].PD_d.time);}
		else {
			mp_to_d(time_to_compare, (*EG_receives)[ii].PD_mp.time); }
		
		solve_options.increment_num_paths_tracked();
		if (((*EG_receives)[ii].retVal != 0 && time_to_compare->r > solve_options.T.minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
			
			trackCount->failures++;
			
			if (solve_options.verbose_level()>=1) {
				printf("\nthere was a path failure tracking witness point %d\nretVal = %d; issoln = %d\n",
					   (*EG_receives)[ii].pathNum, (*EG_receives)[ii].retVal, issoln);
				
				print_path_retVal_message((*EG_receives)[ii].retVal);
				
				//                (point_d out_d, point_mp out_mp,
				//                 int *out_prec,
				//                 point_d in_d, point_mp in_mp,
				//                 int in_prec,
				//                 void const *ED_d, void const *ED_mp)
				//
				//
				if (solve_options.verbose_level() >= 4) {
					if ((*EG_receives)[ii].prec < 64){
                        int out_prec;
                        vec_d temp; init_vec_d(temp,0);
                        ED_d->dehomogenizer(temp, NULL, &out_prec, (*EG_receives)[ii].PD_d.point,NULL,52,ED_d, NULL);
						print_point_to_screen_matlab(temp,"bad_terminal_point");
                        print_comp_matlab((*EG_receives)[ii].PD_d.time,"time");
                        clear_vec_d(temp);
                    }
					else{
                        int out_prec;
                        vec_mp temp; init_vec_mp(temp,0);
                        ED_mp->dehomogenizer(NULL,temp, &out_prec, NULL,(*EG_receives)[ii].PD_mp.point,72,NULL,ED_mp);
						print_point_to_screen_matlab(temp,"bad_terminal_point");
                        print_comp_matlab((*EG_receives)[ii].PD_mp.time,"time");
                        clear_vec_mp(temp);
						
                    }
				}
			}
            
		}
		else
		{
			//otherwise converged, but may have still had non-zero retval due to other reasons.
			
			
			//this conversion of type is total crap. i mean, really, what is the point?  all it does is put the data into a more useless format and waste time.
			endgamedata_to_endpoint(&endPoints[solution_counter], &((*EG_receives)[ii]));
			
			trackCount->successes++;
			solution_counter++; // probably this could be eliminated
		}
	}
	
	
	
	
	
	return incoming_id;
}







void generic_track_path(int pathNum, endgame_data_t *EG_out,
                        point_data_d *Pin, point_data_mp *Pin_mp,
                        FILE *OUT, FILE *MIDOUT,
                        tracker_config_t *T,
                        void const *ED_d, void const *ED_mp,
                        int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
                        int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
                        int (*change_prec)(void const *, int),
                        int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
//	print_tracker(T);
	
	EG_out->pathNum = pathNum;
	EG_out->codim = 0; // this is ignored
	
    T->first_step_of_path = 1;
    T->endgameSwitch = 0;
	
    if (T->MPType == 2)
    { // track using AMP
        
        change_prec(ED_mp,64);
        T->Precision = 64;
        EG_out->prec = EG_out->last_approx_prec = 52;
        
        EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
        
        if (EG_out->prec == 52)
        { // copy over values in double precision
            EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
            EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
            EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
            findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
        }
        else
        { // make sure that the other MP things are set to the correct precision
            mpf_clear(EG_out->function_residual_mp);
            mpf_init2(EG_out->function_residual_mp, EG_out->prec);
            
            mpf_clear(EG_out->latest_newton_residual_mp);
            mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);
            
            mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
            mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);
            
            mpf_clear(EG_out->error_at_latest_sample_point_mp);
            mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);
            
            // copy over the values
            mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
            mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
            mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
            findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
        }
    }
    else if (T->MPType == 0)
    { // track using double precision
        EG_out->prec = EG_out->last_approx_prec = 52;
        
        EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);  // WHERE THE ACTUAL TRACKING HAPPENS
        EG_out->first_increase = 0;
        // copy over values in double precision
        EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
        EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
        EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
        findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
    }
    else if (T->MPType == 1)
    {

        // track using MP
        EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin_mp, T, OUT, MIDOUT, ED_mp, eval_func_mp, find_dehom);
        
        
        EG_out->prec = EG_out->last_approx_prec = T->Precision;
        EG_out->first_increase = 0;
        
        // copy over the values
        mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
        mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
        mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
        findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
        
    }
    
	
//	if (EG_out->retVal !=0) {
//		
//		comp_d time_to_compare;
//		if (EG_out->prec < 64) {
//			set_d(time_to_compare,EG_out->PD_d.time);}
//		else {
//			mp_to_d(time_to_compare, EG_out->PD_mp.time); }
//		
//		
//		vec_d solution_as_double; init_vec_d(solution_as_double,0);
//		if (EG_out->prec < 64){
//			int out_prec;
//			SolverDoublePrecision * evalll = (SolverDoublePrecision *)ED_d;
//			evalll->dehomogenizer(solution_as_double, NULL, &out_prec, EG_out->PD_d.point,NULL,52,evalll, NULL);
//			
//			
//		}
//		else{
//			SolverMultiplePrecision * evalll = (SolverMultiplePrecision *)ED_mp;
//			
//			int out_prec;
//			vec_mp temp2; init_vec_mp(temp2,0);
//			evalll->dehomogenizer(NULL,temp2, &out_prec, NULL,EG_out->PD_mp.point,EG_out->prec,NULL,evalll);
//			
//			vec_mp_to_d(solution_as_double,temp2);
//			clear_vec_mp(temp2);
//			
//		}
//		
//		if (1) {
//			print_point_to_screen_matlab(solution_as_double,"failed_candidate_solution");
//			print_comp_matlab(time_to_compare,"corresponding_time");
//		}
//		
//		clear_vec_d(solution_as_double);
//	}
	
	
	return;
}



void robust_track_path(int pathNum, endgame_data_t *EG_out,
                       point_data_d *Pin, point_data_mp *Pin_mp,
                       FILE *OUT, FILE *MIDOUT,
                       SolverConfiguration & solve_options,
                       SolverDoublePrecision *ED_d, SolverMultiplePrecision *ED_mp,
                       int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
                       int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
                       int (*change_prec)(void const *, int),
                       int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
	
    //	std::cout << "using robust tracker" << std::endl;
	EG_out->pathNum = pathNum;
	EG_out->codim = 0; // this is ignored
	
	
	tracker_config_t * T = &solve_options.T;
	
//	print_tracker(T);
	
	int iterations=0, max_iterations = 3;
	
	solve_options.backup_tracker_config("robust_init");
	
	
	
	std::map <int,int> setting_increments;
    
	solve_options.increment_num_paths_tracked();
	
	
	EG_out->retVal = -876; // set to bad return value
	while ((iterations<max_iterations) && (EG_out->retVal!=0)) {
		
		if (solve_options.verbose_level()>=4) {
			std::cout << color::gray() << "\t\tpath " << pathNum << ", pass " << iterations << color::console_default() << std::endl;
		}
		// reset a few things here
		
		EG_out->retVal = 0;
		T->first_step_of_path = 1;
		T->endgameSwitch = 0;
		
		if (T->MPType == 2)
		{ // track using AMP
			
			change_prec(ED_mp,64);
			T->Precision = 64;
			
			
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase,
                                         &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec,
                                         EG_out->last_approx_d, EG_out->last_approx_mp,
                                         Pin,
                                         T,
                                         OUT, MIDOUT,
                                         ED_d, ED_mp,
                                         eval_func_d, eval_func_mp, change_prec, find_dehom);
			
			if (EG_out->prec == 52)
			{ // copy over values in double precision
				EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
				EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
				EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
				findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
			}
			else
			{ // make sure that the other MP things are set to the correct precision
				mpf_clear(EG_out->function_residual_mp);
				mpf_init2(EG_out->function_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->latest_newton_residual_mp);
				mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
				mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);
				
				mpf_clear(EG_out->error_at_latest_sample_point_mp);
				mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);
				
				// copy over the values
				mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
				mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
				mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
				findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			}
		}
		else if (T->MPType == 0)
		{ // track using double precision
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);  // WHERE THE ACTUAL TRACKING HAPPENS
			EG_out->first_increase = 0;
			// copy over values in double precision
			EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
			EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
			EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
			findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
		}
		else if (T->MPType == 1)
		{
			EG_out->pathNum = pathNum;
			EG_out->codim = 0; // zero dimensional - this is ignored
			
			T->first_step_of_path = 1;
			
			// track using MP
			EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin_mp, T, OUT, MIDOUT, ED_mp, eval_func_mp, find_dehom);
			
			
			EG_out->prec = EG_out->last_approx_prec = T->Precision;
			EG_out->first_increase = 0;
			
			// copy over the values
			mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
			mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
			mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
			findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			
		}
		
		
		
		
		
		
		
		
		
		
		// get how many times we have changed settings due to this type of failure.
		int current_retval_counter = map_lookup_with_default( setting_increments, EG_out->retVal, 0 ); // how many times have we encountered this retval?
		
		if ( EG_out->retVal!=0 ) {  // ||   EG_out->retVal==-50
			
			vec_d solution_as_double; init_vec_d(solution_as_double,0);
			if (EG_out->prec < 64){
				int out_prec;
				
				ED_d->dehomogenizer(solution_as_double, NULL, &out_prec, EG_out->PD_d.point,NULL,52,ED_d, NULL);
				
				
			}
			else{
				int out_prec;
				vec_mp temp2; init_vec_mp(temp2,0);
				ED_mp->dehomogenizer(NULL,temp2, &out_prec, NULL,EG_out->PD_mp.point,EG_out->prec,NULL,ED_mp);
				
				vec_mp_to_d(solution_as_double,temp2);
				clear_vec_mp(temp2);
				
			}
			
			
			
			comp_d time_to_compare;
			if (EG_out->prec < 64) {
				set_d(time_to_compare,EG_out->PD_d.time);}
			else {
				mp_to_d(time_to_compare, EG_out->PD_mp.time); }
			
			
			if (solve_options.verbose_level()>=4) {
				print_point_to_screen_matlab(solution_as_double,"failed_candidate_solution");
				print_comp_matlab(time_to_compare,"time");
			}
			
			
			
			// if
			if ( (time_to_compare->r < std::max(1e-3,1e-2*solve_options.T.endgameBoundary)) && (infNormVec_d(solution_as_double) > solve_options.T.finiteThreshold)) {
				if (solve_options.verbose_level()>=2) {
					print_point_to_screen_matlab(solution_as_double,"big_solution");
					print_comp_matlab(time_to_compare,"at_time");
					std::cout << "discarding non-finite solution.\n\n" << std::endl;
				}
				EG_out->retVal=0;
				break;
			}
			
			clear_vec_d(solution_as_double);
			
			
			if (solve_options.verbose_level()>=3){
				print_path_retVal_message(EG_out->retVal);
				std::cout << color::red() << "solution had non-zero retVal " << EG_out->retVal << " (" << current_retval_counter << ")th occurrence on iteration " << iterations << "." << color::console_default() << std::endl;
			}
			
			switch (EG_out->retVal) {
					
				case -20: // // refining failed
				case -50: // some other failure
					
					
				case -100:   //this is higher precision needed.
							 // break deliberately omitted
					
				case 100:
					
                    solve_options.T.endgameNumber = 2;
                    
					
					if (iterations>0) { // changing predictor
						solve_options.T.odePredictor  = (solve_options.T.odePredictor+1) %9;
					}
					if (iterations>=2) {
						solve_options.T.maxStepSize *= 0.1;
					}
					
					break;
					
				case -10:
					
					solve_options.T.maxNumSteps *=2; // factor of 2 each time
					break;
					
					
				case -3: // minstepsize
					
					solve_options.T.minStepSizeBeforeEndGame *= 1e-1;
					solve_options.T.minStepSizeDuringEndGame *= 1e-1;
					solve_options.T.minStepSize *=  1e-1;
					
					if (iterations<1) {
						solve_options.T.endgameNumber = 2;
						solve_options.T.odePredictor  = MIN(8,solve_options.T.odePredictor+1);
					}
					
					break;
					
				case -4: // securitymax
					
					if (iterations<2) {
						solve_options.T.securityMaxNorm *= 10;  // exponential increase by 10's
																//						std::cout << "increasing securityMaxNorm to " << solve_options.T.securityMaxNorm << std::endl;
					}
					else
					{
						// on the third try, go to security level 1.
						solve_options.T.securityLevel = 1; // just turn on security level 1
														   //						std::cout << "setting securityLevel to 1" << std::endl;
					}
                    
					
					break;
				case -2:
					if (iterations<2) {
						solve_options.T.goingToInfinity *= 10;  // exponential increase by 10's
					}
					else
					{
						solve_options.T.goingToInfinity *= 10;  // exponential increase by 10's
																// on the manyth try, go to security level 1.
						solve_options.T.securityLevel = 1;
					}
					
					break;
					
				case -200: // cycle number too high
					solve_options.T.endgameNumber = 2;
					solve_options.T.cycle_num_max +=3 ;
					break;
				default:
					
					std::cout << color::red() << "retVal was of unprogrammed robust changing: " << EG_out->retVal << color::console_default() << std::endl;
					print_path_retVal_message(EG_out->retVal);
 					break;
			}
		}
		
		// increment the counter for how many times we have changed this type of setting.
		setting_increments[EG_out->retVal] = current_retval_counter + 1;
		
		
		iterations++;
	} // re: while
	
	if (solve_options.verbose_level()>=3) {
		if (iterations==0 && EG_out->retVal==0) {
			std::cout << "success path " << pathNum << std::endl;
		}
		if (iterations>1 && EG_out->retVal==0) {
			std::cout << "resolution of " << pathNum << " was successful" << std::endl;
		}
		
		if (iterations>1 && EG_out->retVal!=0) {
			std::cout << "resolution of path " << pathNum << " failed, terminal retVal " << EG_out->retVal << std::endl;
			//		print_tracker(T);
			//		mypause();
		}
	}

	solve_options.restore_tracker_config("robust_init");
	
	return;
} // re: robust_track_path










void generic_setup_patch(patch_eval_data_d *P, const WitnessSet & W)
{
//	std::cout << "setting up double patch " << std::endl;
	
	
	if (W.num_patches()==0) {
		std::cerr << "the number of patches in input W is 0.  this is not allowed, the number must be positive.\n" << std::endl;
		br_exit(1800);
	}
	
	int total_num_vars_in_patches = 0;
	for (unsigned int ii=0; ii<W.num_patches(); ++ii) {
		total_num_vars_in_patches += W.patch(ii)->size;
	}
	
	if (total_num_vars_in_patches > W.num_variables()) {
		std::cout << "parity mismatch in patches ("<< total_num_vars_in_patches <<") and number of variables (" << W.num_variables() <<")." << std::endl;
		for (unsigned int ii=0; ii<W.num_patches(); ++ii) {
			std::cout << W.patch(ii)->size << " ";
		}
		std::cout << std::endl;
		br_exit(4012);
	}
	
	
	
	P->num_patches = W.num_patches();
	init_mat_d(P->patchCoeff, int(W.num_patches()), W.num_variables());
	P->patchCoeff->rows = W.num_patches(); P->patchCoeff->cols = W.num_variables();
	
	int varcounter = 0;
	for (unsigned int jj=0; jj<W.num_patches(); jj++) {
		vec_mp & curr_patch = W.patch(jj);
		
		for (int ii=0; ii<varcounter; ii++) {
			set_zero_d(&P->patchCoeff->entry[jj][ii]);
		}
		
		int offset = varcounter;
		for (int ii = 0; ii < curr_patch->size ; ii++){
			mp_to_d(&P->patchCoeff->entry[jj][ii+offset],&curr_patch->coord[ii]);
			varcounter++;
		}
		
		for (int ii=varcounter; ii<W.num_variables(); ii++) {
			set_zero_d(&P->patchCoeff->entry[jj][ii]);
		}
	}
}


void generic_setup_patch(patch_eval_data_mp *P, const WitnessSet & W)
{
	
	
	if (W.num_patches()==0) {
		std::cerr << "the number of patches in input W is 0.  this is not allowed, the number must be positive.\n" << std::endl;
		br_exit(1801);
	}
	
	
    int total_num_vars_in_patches = 0;
	for (unsigned int ii=0; ii<W.num_patches(); ++ii) {
		total_num_vars_in_patches += W.patch(ii)->size;
	}
	
	if (total_num_vars_in_patches > W.num_variables()) {
		std::cout << "parity mismatch in patches ("<< total_num_vars_in_patches <<") and number of variables (" << W.num_variables() <<")." << std::endl;
		for (unsigned int ii=0; ii<W.num_patches(); ++ii) {
			std::cout << W.patch(ii)->size << " ";
		}
		std::cout << std::endl;
		br_exit(4013);
	}
	
	
	
//	if (total_num_vars_in_patches < W.num_variables) {
//		std::cout << "parity mismatch in patches ("<< total_num_vars_in_patches <<") and number of variables ("<< W.num_variables <<")." << std::endl;\
//	}
	

	
	
	init_mat_rat(P->patchCoeff_rat, int(W.num_patches()), W.num_variables());
	init_mat_mp2(P->patchCoeff, int(W.num_patches()), W.num_variables(), mpf_get_default_prec());
	
	P->curr_prec = mpf_get_default_prec();
	P->num_patches = W.num_patches();
	P->patchCoeff->rows = W.num_patches();
	P->patchCoeff->cols = W.num_variables();
	
	
	
	
	
	
	
	int varcounter = 0;
	for (unsigned int jj=0; jj<W.num_patches(); jj++) {
		vec_mp & curr_patch = W.patch(jj);
		for (int ii=0; ii<varcounter; ii++) {
			set_zero_mp(&P->patchCoeff->entry[jj][ii]);
			set_zero_rat(P->patchCoeff_rat[jj][ii]);
		}
		
		int offset = varcounter;
		for (int ii = 0; ii < curr_patch->size; ii++){
			set_mp(&P->patchCoeff->entry[jj][ii+offset],&curr_patch->coord[ii]);
			mp_to_rat(P->patchCoeff_rat[jj][ii+offset], &curr_patch->coord[ii]);
			varcounter++;
		}
		
		for (int ii=varcounter; ii<W.num_variables(); ii++) {
			set_zero_mp(&P->patchCoeff->entry[jj][ii]);
			set_zero_rat(P->patchCoeff_rat[jj][ii]);
		}
	}
}



int generic_setup_files(FILE ** OUT, boost::filesystem::path outname,
                        FILE ** MIDOUT, boost::filesystem::path midname)
{
	
	*OUT = safe_fopen_write(outname);  // open the main output files.
    *MIDOUT = safe_fopen_write(midname);
	
	return SUCCESSFUL;
}




















