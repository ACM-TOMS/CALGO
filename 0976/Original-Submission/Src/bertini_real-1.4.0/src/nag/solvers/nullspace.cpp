#include "nag/solvers/nullspace.hpp"







void NullspaceConfiguration::clear()
{
	
	
	
	
	
	
	if (num_jac_equations>0) {
		if (side_ == nullspace_handedness::LEFT) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					clear_vec_mp(starting_linears[ii][jj]);
				}
				free(starting_linears[ii]);
			}
		}
		else
		{
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					clear_vec_mp(starting_linears[ii][jj]);
				}
				free(starting_linears[ii]);
			}
		}
		
		free(starting_linears);
	}
	
	
	num_jac_equations = 0;
	max_degree = 0;
	
	if (num_additional_linears >0 ) {
		for (int ii=0; ii<num_additional_linears; ii++){
			clear_vec_mp(additional_linears_starting[ii]);
			clear_vec_mp(additional_linears_terminal[ii]);
		}
		free(additional_linears_terminal);
		free(additional_linears_starting);
		
		num_additional_linears = 0;
	}
	
	
	if (num_v_linears>0) {
		for (int ii=0; ii<num_v_linears; ii++)
			clear_vec_mp(v_linears[ii]);
		free(v_linears);
		clear_vec_mp(v_patch);
		
		num_v_linears = 0;
	}
	
	
	if (num_projections>0) {
		for (int ii=0; ii<num_projections; ii++)
			clear_vec_mp(target_projection[ii]);
		free(target_projection);
		
		num_projections = 0;
	}
	
}




void NullspaceConfiguration::print()
{
	std::stringstream converter;
	
	std::cout << "*******************\n\tns_config:\n\n";
	
	std::cout << "#nat: " << num_natural_vars << " " << " #synth: " << num_synth_vars << " #v: " << num_v_vars << std::endl;
	
	
	
	std::cout << "# additional linears: " << this->num_additional_linears << std::endl;
	std::cout << "# jacobian equations: " << this->num_jac_equations << std::endl;
	std::cout << "# randomized equations: " << this->randomizer_->num_rand_funcs() << std::endl;
	std::cout << "max_degree:\t" << max_degree << std::endl;
	
	std::cout << *(this->randomizer_) << std::endl;
	
	
	for (int ii=0; ii<num_additional_linears; ii++) {
		converter << "additional_linear_starting_" << ii+1;
		print_point_to_screen_matlab(additional_linears_starting[ii],converter.str());
		converter.str("");
	}
	
	
	for (int ii=0; ii<num_additional_linears; ii++) {
		converter << "additional_linear_terminal_" << ii+1;
		print_point_to_screen_matlab(additional_linears_terminal[ii],converter.str());
		converter.str("");
	}
	
	
	if (side_==nullspace_handedness::LEFT) {
		for (int ii=0; ii<num_jac_equations; ii++) {
			for (int jj=0; jj<this->max_degree; jj++) {
				converter << "starting_linear_" << ii+1 << "_" << jj+1;
				print_point_to_screen_matlab(starting_linears[ii][jj],converter.str());
				converter.str("");
			}
		}
	}
	else
	{
		for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
			for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
				converter << "starting_linear_" << ii << "_" << jj;
				print_point_to_screen_matlab(starting_linears[ii][jj],converter.str());
				converter.str("");
			}
		}
	}
	
	
	for (int ii=0; ii<num_jac_equations; ii++) {
		converter << "v_linear_" << ii+1;
		print_point_to_screen_matlab(v_linears[ii],converter.str());
		converter.str("");
	}
	
	for (int ii=0; ii<num_projections; ii++) {
		converter << "pi" << ii;
		print_point_to_screen_matlab(target_projection[ii],converter.str());
		converter.str("");
	}
	std::cout << "\n********************\n" << std::endl;
	
}


///////////////
//
//   end the NullspaceConfiguration
//
/////////////






///////////////
//
//   begin nullspace_eval_data_mp
//
/////////////

void nullspacejac_eval_data_mp::set_sidedness(int which_side)
{
	this->side_ = which_side;
	set_function_handles(which_side);
}

void nullspacejac_eval_data_mp::set_function_handles(int which_side)
{
	
	if (which_side==nullspace_handedness::LEFT) {
		this->evaluator_function_d = &nullspacejac_left_eval_d;
		this->evaluator_function_mp = &nullspacejac_left_eval_mp;
	}
	else if(which_side==nullspace_handedness::RIGHT)
	{
		this->evaluator_function_d = &nullspacejac_right_eval_d;
		this->evaluator_function_mp = &nullspacejac_right_eval_mp;
	}
	else
	{
		throw std::logic_error("sidedness incorrect in set_function_handles");
	}
}

void nullspacejac_eval_data_mp::init()
{
	this->is_solution_checker_d = &check_issoln_nullspacejac_d;
	this->is_solution_checker_mp = &check_issoln_nullspacejac_mp;
	
	this->evaluator_function_d = NULL;
	this->evaluator_function_mp = NULL;
	
	this->precision_changer = &change_nullspacejac_eval_prec;
	this->dehomogenizer = &nullspacejac_dehom;
	
	side_ = -1;
	
	
	additional_linears_terminal = NULL;
	additional_linears_terminal_full_prec = NULL;
	
	additional_linears_starting = NULL;
	additional_linears_starting_full_prec = NULL;
	
	
	starting_linears = NULL; // outer layer should have as many as there are randomized equations
							 // inside layer has number corresponding to randomized_degrees
	starting_linears_full_prec = NULL; // outer layer should have as many as there are randomized equations
									   // inside layer has number corresponding to randomized_degrees
	
	v_linears = NULL;         // should be as many in here as there are randomized equations
	v_linears_full_prec = NULL;         // should be as many in here as there are randomized equations
	
	init_vec_mp(v_patch,0);
	
	SLP_derivative = new prog_deriv_t;
	
	init_mat_mp(jac_with_proj,0,0);
	

	
	
	
	target_projection = NULL; //
	target_projection_full_prec = NULL; //
	
	
	if (this->MPType==2) {
		init_vec_mp2(v_patch_full_prec,0,1024);
		init_mat_mp2(jac_with_proj_full_prec,0,0,1024);
	}
	
	
}


int nullspacejac_eval_data_mp::send(ParallelismConfig & mpi_config)
{
	
	int solver_choice = NULLSPACE;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	SolverMultiplePrecision::send(mpi_config);
	
	
	int *buffer = new int[12];
	
	buffer[0] = num_additional_linears;
	buffer[1] = num_jac_equations;
	buffer[2] = max_degree;
	buffer[3] = num_v_linears;
	buffer[4] = num_projections;
	
	buffer[5] = num_natural_vars;
	buffer[6] = num_v_vars;
	
	buffer[7] = target_dim;
	buffer[8] = ambient_dim;
	buffer[9] = target_crit_codim;
	
	buffer[10] = num_synth_vars;
	buffer[11] = side_;
	// now can actually send the data.
	
	
	
	MPI_Bcast(buffer,12,MPI_INT, mpi_config.id(), mpi_config.comm());
	
	delete[] buffer;
	
	
	if (this->MPType==2){
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				// send the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal_full_prec[ii], mpi_config.id(), mpi_config.head());
				
				// send the full precision starting additional linear
				bcast_vec_mp(this->additional_linears_starting_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else {} // num_additional_linears == 0
		
		
		
		if (num_jac_equations>0) {
			
			if (side_==nullspace_handedness::LEFT) {
				for (int ii=0; ii<num_jac_equations; ii++) {
					for (int jj=0; jj<max_degree; jj++) {
						bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
					}
				}
			}
			else{
				for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
					for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
						bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
					}
				}
			}
			
		}
		else{}
		
		
		
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				bcast_vec_mp(v_linears_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}
		
		bcast_vec_mp(v_patch_full_prec, mpi_config.id(), mpi_config.head());
		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				bcast_vec_mp(target_projection_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}
	}
	else {
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
				
				// receive the full precision starting additional linear
				bcast_vec_mp(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else {} // num_additional_linears == 0
		
		
		// recieve the post-randomizer-matrix
		
		
		if (side_==nullspace_handedness::LEFT) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		else{
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		
		
		for (int ii=0; ii<num_v_linears; ii++) {
			bcast_vec_mp(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
		
		
		bcast_vec_mp(v_patch, mpi_config.id(), mpi_config.head());
		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				bcast_vec_mp(target_projection[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}

	}
	
	
	return SUCCESSFUL;
}




int nullspacejac_eval_data_mp::receive(ParallelismConfig & mpi_config)
{
	int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, mpi_config.comm());
	
	if (buffer[0] != NULLSPACE) {
		std::cout << "worker failed to confirm it is receiving the NULLSPACE type eval data" << std::endl;
		mpi_config.abort(777);
	}
	
	SolverMultiplePrecision::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	MPI_Bcast(buffer,12,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	num_additional_linears = buffer[0];
	num_jac_equations = buffer[1];
	max_degree = buffer[2];
	num_v_linears = buffer[3];
	num_projections = buffer[4];
	
	num_natural_vars = buffer[5];
	num_v_vars = buffer[6];
	
	target_dim = buffer[7];
	ambient_dim = buffer[8];
	target_crit_codim = buffer[9];
	
	num_synth_vars = buffer[10];
	
	set_sidedness(buffer[11]);
	
	delete[] buffer;
	
	
	
	
	if (this->MPType==2) {
		if (num_additional_linears>0) {
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_terminal_full_prec = (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting_full_prec = (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				init_vec_mp2(this->additional_linears_terminal_full_prec[ii],0,1024);
				init_vec_mp(this->additional_linears_terminal[ii],0);
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->additional_linears_terminal[ii],this->additional_linears_terminal_full_prec[ii]);
				
				// receive the full precision starting additional linear
				init_vec_mp2(this->additional_linears_starting_full_prec[ii],0,1024);
				init_vec_mp(this->additional_linears_starting[ii],0);
				bcast_vec_mp(this->additional_linears_starting_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->additional_linears_starting[ii],this->additional_linears_starting_full_prec[ii]);
				
			}
		}
		else {} // num_additional_linears == 0
		
		
		
		
		if (num_jac_equations>0) {
			
			if (side_==nullspace_handedness::LEFT) {
				this->starting_linears_full_prec = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
				this->starting_linears = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
				
				for (int ii=0; ii<num_jac_equations; ii++) {
					
					this->starting_linears_full_prec[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
					this->starting_linears[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
					for (int jj=0; jj<max_degree; jj++) {
						init_vec_mp(this->starting_linears[ii][jj],0);
						init_vec_mp2(this->starting_linears_full_prec[ii][jj],0,1024);
						
						// recieve the starting linearsin full prec and convert
						bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
						vec_cp_mp(this->starting_linears[ii][jj],starting_linears_full_prec[ii][jj]);
					}
				}
			}
			else{
				this->starting_linears_full_prec = (vec_mp **) br_malloc(randomizer()->num_rand_funcs()*sizeof(vec_mp *));
				this->starting_linears = (vec_mp **) br_malloc(randomizer()->num_rand_funcs()*sizeof(vec_mp *));
				
				for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
					
					int curr_degree = randomizer()->randomized_degree(ii)-1;
					
					this->starting_linears_full_prec[ii] = (vec_mp *) br_malloc(curr_degree*sizeof(vec_mp ));
					this->starting_linears[ii] = (vec_mp *) br_malloc(curr_degree*sizeof(vec_mp ));
					for (int jj=0; jj<curr_degree; jj++) {
						init_vec_mp(this->starting_linears[ii][jj],0);
						init_vec_mp2(this->starting_linears_full_prec[ii][jj],0,1024);
						
						// recieve the starting linearsin full prec and convert
						bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
						vec_cp_mp(this->starting_linears[ii][jj],starting_linears_full_prec[ii][jj]);
					}
				}
			}
			
			
		}
		else{}
		
		if (num_v_linears>0) {
			this->v_linears = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
			this->v_linears_full_prec = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_mp(this->v_linears[ii],0);
				init_vec_mp2(this->v_linears_full_prec[ii],0,1024);
				// receive the full precision v_linears
				bcast_vec_mp(v_linears_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->v_linears[ii],v_linears_full_prec[ii]);
			}
		}
		else{}
		
		init_vec_mp2(v_patch_full_prec,0,1024);
		init_vec_mp(v_patch,0);
		bcast_vec_mp(v_patch_full_prec, mpi_config.id(), mpi_config.head());
		vec_cp_mp(this->v_patch, v_patch_full_prec);
		
		init_mat_mp2(jac_with_proj_full_prec, num_natural_vars-1, randomizer()->num_rand_funcs()+num_projections,1024);
		
		
		init_mat_mp(jac_with_proj, num_natural_vars-1, randomizer()->num_rand_funcs()+num_projections);
		jac_with_proj->rows = jac_with_proj_full_prec->rows = num_natural_vars-1;
		jac_with_proj->cols = jac_with_proj_full_prec->cols = randomizer()->num_rand_funcs()+num_projections;
		

		
		
		if (num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			this->target_projection_full_prec = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				init_vec_mp2(this->target_projection[ii],0,1024);
				
				bcast_vec_mp(target_projection_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->target_projection[ii],target_projection_full_prec[ii]);
			}
		}
		else{}
	}
	else{ // MPType == 1
		if (num_additional_linears>0) {
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				init_vec_mp2(this->additional_linears_terminal[ii],0,1024);
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
				
				// receive the full precision starting additional linear
				init_vec_mp2(this->additional_linears_starting[ii],0,1024);
				bcast_vec_mp(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
				
			}
		}
		else {} // num_additional_linears == 0
		
		
		
		if (side_==nullspace_handedness::LEFT) {
			this->starting_linears = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
			
			for (int ii=0; ii<num_jac_equations; ii++) {
				this->starting_linears[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
				for (int jj=0; jj<max_degree; jj++) {
					init_vec_mp(this->starting_linears[ii][jj],num_natural_vars);  starting_linears[ii][jj]->size = num_natural_vars;
					// recieve the starting linearsin full prec and convert
					bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		else{
			this->starting_linears = (vec_mp **) br_malloc(randomizer()->num_rand_funcs()*sizeof(vec_mp *));
			
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				
				int curr_degree = randomizer()->randomized_degree(ii)-1;
				this->starting_linears[ii] = (vec_mp *) br_malloc(curr_degree*sizeof(vec_mp ));
				for (int jj=0; jj<curr_degree; jj++) {
					init_vec_mp(this->starting_linears[ii][jj],num_natural_vars);  starting_linears[ii][jj]->size = num_natural_vars;
					// recieve the starting linearsin full prec and convert
					bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
				
			}
		}
		
		
		this->v_linears = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
		for (int ii=0; ii<num_v_linears; ii++) {
			init_vec_mp(this->v_linears[ii],0);
			// receive the full precision v_linears
			bcast_vec_mp(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
		
		init_vec_mp(v_patch,0);
		bcast_vec_mp(v_patch, mpi_config.id(), mpi_config.head());
		//		std::cout << "a";
		
		init_mat_mp(jac_with_proj,num_natural_vars-1,randomizer()->num_rand_funcs()+num_projections);
		jac_with_proj->rows = num_natural_vars-1;
		jac_with_proj->cols = randomizer()->num_rand_funcs()+num_projections;
		
		
		if (num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				bcast_vec_mp(target_projection[ii], mpi_config.id(), mpi_config.head());
			}
		}
	}
	
	setup_deriv_from_SLP(SLP_derivative, this->SLP);
	while (SLP_derivative->order < 2) {
		setupNext_derivs(SLP_derivative);
	}
	
	
	
	
	return SUCCESSFUL;
}

int nullspacejac_eval_data_mp::setup(prog_t * _SLP,
                                     NullspaceConfiguration *ns_config,
                                     WitnessSet & W,
                                     SolverConfiguration & solve_options)
{
	
	verbose_level(solve_options.verbose_level());
	
	SolverMultiplePrecision::setup(_SLP, ns_config->randomizer());
	
	generic_setup_patch(&patch,W);
	

	
	
	setup_deriv_from_SLP(SLP_derivative, SLP);
	while (SLP_derivative->order < 2) {
		setupNext_derivs(SLP_derivative);
	}
	
	
	comp_d temp;
	if (this->MPType==2) {
		if (solve_options.use_gamma_trick==1){
			get_comp_rand_rat(temp, this->gamma, this->gamma_rat, 64, solve_options.T.AMP_max_prec, 0, 0);
		}
		else{
			set_one_mp(this->gamma);
			set_one_rat(this->gamma_rat);
		}
	}
	else{
        if (solve_options.use_gamma_trick==1)
            get_comp_rand_mp(this->gamma); // set gamma to be random complex value
        else{
            set_one_mp(this->gamma);
        }
    }
    
    set_sidedness(ns_config->side());
	
	
	num_jac_equations = ns_config->num_jac_equations;
	target_dim = ns_config->target_dim;
	ambient_dim = ns_config->ambient_dim;
	target_crit_codim = ns_config->target_crit_codim;
	
	num_v_vars = ns_config->num_v_vars;
	num_natural_vars = ns_config->num_natural_vars;
	num_synth_vars = ns_config->num_synth_vars;
	num_variables = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	
	num_v_linears = ns_config->num_v_linears;   //
	
	num_projections = ns_config->num_projections;
	
	num_additional_linears = ns_config->num_additional_linears;
	
	max_degree = ns_config->max_degree;

	this->set_randomizer(ns_config->randomizer());
	
	
	
	
	// set up the vectors to hold the linears.
	
	target_projection = (vec_mp *)br_malloc(ns_config->num_projections * sizeof(vec_mp));
	
	init_mat_mp(jac_with_proj, ns_config->num_natural_vars-1,ns_config->num_v_vars);
	jac_with_proj->rows = ns_config->num_natural_vars-1;
	jac_with_proj->cols = ns_config->num_v_vars;
	
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_mp(target_projection[ii],ns_config->num_natural_vars+ ns_config->num_synth_vars);
		target_projection[ii]->size =  ns_config->num_natural_vars+ ns_config->num_synth_vars;
		vec_cp_mp(target_projection[ii], ns_config->target_projection[ii]);
	}
	
	
	vec_cp_mp(v_patch, ns_config->v_patch);
	
	
	
	
	this->v_linears = (vec_mp *) br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp(v_linears[ii],ns_config->num_v_vars);
		v_linears[ii]->size = ns_config->num_v_vars;
		vec_cp_mp(v_linears[ii], ns_config->v_linears[ii]);
	}
	
	
	additional_linears_terminal = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
	additional_linears_starting = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_mp(additional_linears_terminal[ii], ns_config->num_natural_vars+ ns_config->num_synth_vars);
		additional_linears_terminal[ii]->size = ns_config->num_natural_vars+ ns_config->num_synth_vars;
		vec_cp_mp(additional_linears_terminal[ii],ns_config->additional_linears_terminal[ii]);
		
		init_vec_mp(additional_linears_starting[ii], ns_config->num_natural_vars+ ns_config->num_synth_vars);
		additional_linears_starting[ii]->size = ns_config->num_natural_vars+ ns_config->num_synth_vars;
		vec_cp_mp(additional_linears_starting[ii],ns_config->additional_linears_starting[ii]);
	}
	
	
	if (side_ == nullspace_handedness::LEFT) {
		starting_linears = (vec_mp **)br_malloc(ns_config->num_jac_equations*sizeof(vec_mp *));
		
		for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
			starting_linears[ii] = (vec_mp *)br_malloc(ns_config->max_degree*sizeof(vec_mp));
			for (int jj=0; jj<ns_config->max_degree; jj++) {
				init_vec_mp(starting_linears[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
				starting_linears[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
				
				vec_cp_mp(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
			}
		}
	}
	else{ // right hand nullspace
		starting_linears = (vec_mp **)br_malloc(ns_config->randomizer()->num_rand_funcs()*sizeof(vec_mp *));
		
		for (int ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ++ii) {
			int curr_degree = ns_config->randomizer()->randomized_degree(ii)-1;
			
			starting_linears[ii] = (vec_mp *)br_malloc(curr_degree*sizeof(vec_mp));
			for (int jj=0; jj<curr_degree; jj++) {
				init_vec_mp(starting_linears[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
				starting_linears[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
				
				vec_cp_mp(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
			}
		}
	}
	
	
	
	if (this->MPType==2) {

		
		
		// set up the vectors to hold the linears.
		
		target_projection_full_prec = (vec_mp *)br_malloc(ns_config->num_projections * sizeof(vec_mp));
		
		init_mat_mp2(jac_with_proj_full_prec, ns_config->num_natural_vars-1,ns_config->num_v_vars,solve_options.T.AMP_max_prec);
		jac_with_proj_full_prec->rows = ns_config->num_natural_vars-1;
		jac_with_proj_full_prec->cols = ns_config->num_v_vars;
		
		for (int ii=0; ii<ns_config->num_projections; ii++) {
			init_vec_mp2(target_projection_full_prec[ii],ns_config->num_natural_vars+ ns_config->num_synth_vars,solve_options.T.AMP_max_prec);
			target_projection_full_prec[ii]->size =  ns_config->num_natural_vars+ ns_config->num_synth_vars;
			vec_cp_mp(target_projection_full_prec[ii], ns_config->target_projection[ii]);

		}
		
		
		vec_cp_mp(v_patch_full_prec, ns_config->v_patch);
		
		
		
		
		this->v_linears_full_prec = (vec_mp *) br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
		for (int ii=0; ii<ns_config->num_v_linears; ii++) {
			init_vec_mp2(v_linears_full_prec[ii],ns_config->num_v_vars,solve_options.T.AMP_max_prec);
			v_linears_full_prec[ii]->size = ns_config->num_v_vars;
			vec_cp_mp(v_linears_full_prec[ii], ns_config->v_linears[ii]);
		}
		
		
		additional_linears_terminal_full_prec = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
		additional_linears_starting_full_prec = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
		
		for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
			init_vec_mp2(additional_linears_terminal_full_prec[ii], ns_config->num_natural_vars+ ns_config->num_synth_vars,solve_options.T.AMP_max_prec);
			additional_linears_terminal_full_prec[ii]->size = ns_config->num_natural_vars+ ns_config->num_synth_vars;
			vec_cp_mp(additional_linears_terminal_full_prec[ii],ns_config->additional_linears_terminal[ii]);
			
			init_vec_mp2(additional_linears_starting_full_prec[ii], ns_config->num_natural_vars+ ns_config->num_synth_vars,solve_options.T.AMP_max_prec);
			additional_linears_starting_full_prec[ii]->size = ns_config->num_natural_vars+ ns_config->num_synth_vars;
			vec_cp_mp(additional_linears_starting_full_prec[ii],ns_config->additional_linears_starting[ii]);
		}
		
		
		if (side_ == nullspace_handedness::LEFT) {
			starting_linears_full_prec = (vec_mp **)br_malloc(ns_config->num_jac_equations*sizeof(vec_mp *));
			
			for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
				starting_linears_full_prec[ii] = (vec_mp *)br_malloc(ns_config->max_degree*sizeof(vec_mp));
				for (int jj=0; jj<ns_config->max_degree; jj++) {
					init_vec_mp(starting_linears_full_prec[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
					starting_linears_full_prec[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
					
					vec_cp_mp(starting_linears_full_prec[ii][jj], ns_config->starting_linears[ii][jj]);
				}
			}
		}
		else{ // right hand nullspace
			starting_linears_full_prec = (vec_mp **)br_malloc(ns_config->randomizer()->num_rand_funcs()*sizeof(vec_mp *));
			
			for (int ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ++ii) {
				int curr_degree = ns_config->randomizer()->randomized_degree(ii)-1;
				
				starting_linears_full_prec[ii] = (vec_mp *)br_malloc(curr_degree*sizeof(vec_mp));
				for (int jj=0; jj<curr_degree; jj++) {
					init_vec_mp(starting_linears_full_prec[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
					starting_linears_full_prec[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
					
					vec_cp_mp(starting_linears_full_prec[ii][jj], ns_config->starting_linears[ii][jj]);
				}
			}
		}
	}
	
	
	return SUCCESSFUL;
}




void nullspacejac_eval_data_mp::print()
{
    std::cout << "num_jac_equations " << num_jac_equations << std::endl;
    std::cout << "target_dim " << target_dim << std::endl;
    
    std::cout << "ambient_dim " << ambient_dim << std::endl;
    std::cout << "target_crit_codim " << target_crit_codim << std::endl;
    
    std::cout << "num_natural_vars " << num_natural_vars << std::endl;
    std::cout << "num_synth_vars " << num_synth_vars << std::endl;
    std::cout << "num_v_vars " << num_v_vars << std::endl;
    std::cout << "num_randomized_eqns " << randomizer()->num_rand_funcs() << std::endl;
    std::cout << "max_degree " << max_degree << std::endl;
    
    

	if (MPType==2) {
		if (side_==nullspace_handedness::LEFT) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					std::stringstream name;
					name << "starting" << ii << "_" << jj;
					print_point_to_screen_matlab(starting_linears_full_prec[ii][jj],name.str());
				}
			}
		}
		else{
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					std::stringstream name;
					name << "starting" << ii << "_" << jj;
					print_point_to_screen_matlab(starting_linears_full_prec[ii][jj],name.str());
				}
			}
		}
	}
	else{
		if (side_==nullspace_handedness::LEFT) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					std::stringstream name;
					name << "starting" << ii << "_" << jj;
					print_point_to_screen_matlab(starting_linears[ii][jj],name.str());
				}
			}
		}
		else{
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					std::stringstream name;
					name << "starting" << ii << "_" << jj;
					print_point_to_screen_matlab(starting_linears[ii][jj],name.str());
				}
			}
		}
	}
	

    
    for (int ii=0; ii< num_additional_linears; ii++) {
		if (MPType==2) {
			print_point_to_screen_matlab(additional_linears_terminal_full_prec[ii],"add_term");
			print_point_to_screen_matlab(additional_linears_starting_full_prec[ii],"add_start");
		}
		else{
			print_point_to_screen_matlab(additional_linears_terminal[ii],"add_term");
			print_point_to_screen_matlab(additional_linears_starting[ii],"add_start");
		}
    }
    
    
    for (int ii=0; ii<num_v_linears; ii++) {
		if (MPType==2) {
			print_point_to_screen_matlab(v_linears_full_prec[ii],"v_linears");
		}
		else
		{
			print_point_to_screen_matlab(v_linears[ii],"v_linears");
		}
        
	}
    
	if (MPType==2) {
		print_point_to_screen_matlab(v_patch_full_prec,"v_patch");
		print_matrix_to_screen_matlab(jac_with_proj_full_prec,"jac_with_proj");
	}
	else
	{
		print_point_to_screen_matlab(v_patch,"v_patch");
		print_matrix_to_screen_matlab(jac_with_proj,"jac_with_proj");
	}
    
    
	std::cout << *randomizer() << std::endl;
    
    
    print_comp_matlab(gamma,"gamma");
    
}



///////////////
//
//   end nullspace_eval_data_mp
//
/////////////
















///////////////
//
//   begin nullspace_eval_data_d
//
/////////////


void nullspacejac_eval_data_d::set_sidedness(int which_side)
{
	this->side_ = which_side;
	set_function_handles(which_side);
}

void nullspacejac_eval_data_d::set_function_handles(int which_side)
{
	
	if (which_side==nullspace_handedness::LEFT) {
		this->evaluator_function_d = &nullspacejac_left_eval_d;
		this->evaluator_function_mp = &nullspacejac_left_eval_mp;
	}
	else if(which_side==nullspace_handedness::RIGHT)
	{
		this->evaluator_function_d = &nullspacejac_right_eval_d;
		this->evaluator_function_mp = &nullspacejac_right_eval_mp;
	}
	else
	{
		throw std::logic_error("sidedness incorrect in set_function_handles");
	}
}




void nullspacejac_eval_data_d::init()
{
	
	if (this->MPType==2){
		this->BED_mp = new nullspacejac_eval_data_mp(2);
        SolverDoublePrecision::BED_mp = this->BED_mp;
    }
    
	else{this->BED_mp = NULL; SolverDoublePrecision::BED_mp = NULL;}
	
	
	side_ = -1;
	
	this->is_solution_checker_d = &check_issoln_nullspacejac_d;
	this->is_solution_checker_mp = &check_issoln_nullspacejac_mp;
	this->evaluator_function_d = NULL;
	this->evaluator_function_mp = NULL;
	this->precision_changer = &change_nullspacejac_eval_prec;
	this->dehomogenizer = &nullspacejac_dehom;
	
	additional_linears_terminal = NULL;
	
	additional_linears_starting = NULL;
	
	
	
	starting_linears = NULL; // outer layer should have as many as there are randomized equations
							 // inside layer has number corresponding to randomized_degrees
	
	
	
	v_linears = NULL;         // should be as many in here as there are randomized equations
	
	init_vec_d(v_patch,0);
	
	init_mat_d(jac_with_proj,0,0);
	

	
	target_projection = NULL; //
	
}


int nullspacejac_eval_data_d::send(ParallelismConfig & mpi_config)
{
    
    int solver_choice = NULLSPACE;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
    
    if (this->MPType==2) {
		this->BED_mp->send(mpi_config);
	}
    
    
	
	
	//send the base class stuff.
	SolverDoublePrecision::send(mpi_config);
	
	
	int *buffer = new int[12];
	
	buffer[0] = num_additional_linears;
	buffer[1] = num_jac_equations;
	buffer[2] = max_degree;
	buffer[3] = num_v_linears;
	buffer[4] = num_projections;
	
	buffer[5] = num_natural_vars;
	buffer[6] = num_v_vars;
	
	buffer[7] = target_dim;
	buffer[8] = ambient_dim;
	buffer[9] = target_crit_codim;
	
	buffer[10] = num_synth_vars;
	buffer[11] = side_;
	
	// now can actually send the data.
	
	MPI_Bcast(buffer,12,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	delete[] buffer;
	
	
	
	if (num_additional_linears>0) {
		for (int ii=0; ii<num_additional_linears; ii++) {
			bcast_vec_d(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
			bcast_vec_d(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else {} // num_additional_linears == 0
	
	
	
	
	if (num_jac_equations>0) {
		
		if (side_ == nullspace_handedness::LEFT) {
			
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		else{
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		
		
		
	}
	else{}
	
	
	
	
	if (num_v_linears>0) {
		for (int ii=0; ii<num_v_linears; ii++) {
			bcast_vec_d(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	bcast_vec_d(v_patch, mpi_config.id(), mpi_config.head());
	
	
	if (num_projections>0) {
		for (int ii=0; ii<num_projections; ii++) {
			bcast_vec_d(target_projection[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	
	return SUCCESSFUL;
}

int nullspacejac_eval_data_d::receive(ParallelismConfig & mpi_config)
{
    int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, mpi_config.comm());
	
	if (buffer[0] != NULLSPACE){
		std::cout << "worker failed to confirm it is receiving the nullspace type eval data" << std::endl;
		mpi_config.abort(777);
	}
    
    
    
	if (this->MPType==2) {
		this->BED_mp->receive(mpi_config);
	}
    
    
	
	
	SolverDoublePrecision::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	
	
	
	MPI_Bcast(buffer,12,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	
	num_additional_linears = buffer[0];
	num_jac_equations = buffer[1];
	max_degree = buffer[2];
	num_v_linears = buffer[3];
	num_projections = buffer[4];
	
	num_natural_vars = buffer[5];
	num_v_vars = buffer[6];
	
	target_dim = buffer[7];
	ambient_dim = buffer[8];
	target_crit_codim = buffer[9];
	
	num_synth_vars = buffer[10];
	set_sidedness(buffer[11]);
	
	delete[] buffer;
	
	
	
	
	if (num_additional_linears>0) {
		
		this->additional_linears_terminal						= (vec_d *) br_malloc(num_additional_linears*sizeof(vec_d));
		this->additional_linears_starting						= (vec_d *) br_malloc(num_additional_linears*sizeof(vec_d));
		
		for (int ii=0; ii<num_additional_linears; ii++) {
			
			// receive the full precision terminal additional linear
			bcast_vec_d(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
			
			// receive the full precision starting additional linear
			bcast_vec_d(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
			
		}
	}
	else {} // num_additional_linears == 0
	
	
	
	
	if (num_jac_equations>0) {
		if (side_==nullspace_handedness::LEFT) {
		
			this->starting_linears = (vec_d **) br_malloc(num_jac_equations*sizeof(vec_d *));
			
			for (int ii=0; ii<num_jac_equations; ii++) {
				
				this->starting_linears[ii] = (vec_d *) br_malloc(max_degree*sizeof(vec_d ));
				for (int jj=0; jj<max_degree; jj++) {
					init_vec_d(this->starting_linears[ii][jj],0);
					
					// recieve the starting linearsin full prec and convert
					bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		else{
			
			this->starting_linears = (vec_d **) br_malloc(randomizer()->num_rand_funcs()*sizeof(vec_d *));
			
			for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
				
				this->starting_linears[ii] = (vec_d *) br_malloc( (randomizer()->randomized_degree(ii)-1)*sizeof(vec_d ));
				for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
					init_vec_d(this->starting_linears[ii][jj],0);
					
					// recieve the starting linearsin full prec and convert
					bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
	}
	else{}
	
	if (num_v_linears>0) {
		this->v_linears = (vec_d *) br_malloc(num_v_linears*sizeof(vec_d));
		for (int ii=0; ii<num_v_linears; ii++) {
			init_vec_d(this->v_linears[ii],0);
			// receive the full precision v_linears
			bcast_vec_d(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	bcast_vec_d(v_patch, mpi_config.id(), mpi_config.head());
	
	
	
	
	
	if (num_projections>0) {
		this->target_projection = (vec_d *) br_malloc(num_projections*sizeof(vec_d));
		
		for (int ii=0; ii<num_projections; ii++) {
			init_vec_d(this->target_projection[ii],0);
			
			bcast_vec_d(target_projection[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	
	
	
	if (this->MPType==2)
	{
		this->SLP_derivative = this->BED_mp->SLP_derivative;
	}
	else{
		SLP_derivative = new prog_deriv_t;
		
		setup_deriv_from_SLP(SLP_derivative, this->SLP);
		while (SLP_derivative->order < 2) {
			setupNext_derivs(SLP_derivative);
		}
	}
	
	
	
	return SUCCESSFUL;
}




int nullspacejac_eval_data_d::setup(prog_t * _SLP,
									NullspaceConfiguration *ns_config,
									WitnessSet & W,
									SolverConfiguration & solve_options)
{
	
	
	
	verbose_level(solve_options.verbose_level());
	
	generic_setup_patch(&patch,W);
	
	
	
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_d(this->gamma); // set gamma to be random complex value
	else
		set_one_d(this->gamma);
	
	
	set_sidedness(ns_config->side());
	
	
	num_jac_equations = ns_config->num_jac_equations;
	target_dim = ns_config->target_dim;
	ambient_dim = ns_config->ambient_dim;
	target_crit_codim = ns_config->target_crit_codim;
	
	num_v_vars = ns_config->num_v_vars;
	num_natural_vars = ns_config->num_natural_vars;
	num_synth_vars = ns_config->num_synth_vars;
	num_variables = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	
	
	max_degree = ns_config->max_degree;


	
	num_v_linears = ns_config->num_v_linears;   //
	
	

	
	
	
	// set up the vectors to hold the linears.
	num_projections = ns_config->num_projections;
	target_projection = (vec_d *)br_malloc(ns_config->num_projections * sizeof(vec_d));
	
	init_mat_d(jac_with_proj,ns_config->num_natural_vars-1,ns_config->num_v_vars);
	jac_with_proj->rows = ns_config->num_natural_vars-1;
	jac_with_proj->cols = ns_config->num_v_vars;
	
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_d(target_projection[ii],ns_config->num_natural_vars+ns_config->num_synth_vars);
		target_projection[ii]->size =  ns_config->num_natural_vars+ns_config->num_synth_vars;
		vec_mp_to_d(target_projection[ii], ns_config->target_projection[ii]); // copy in the projection
	}
	
	
	vec_mp_to_d(v_patch, ns_config->v_patch);
	
	
	
	
	this->v_linears = (vec_d *) br_malloc(ns_config->num_v_linears*sizeof(vec_d));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_d(v_linears[ii],ns_config->num_v_vars);
		v_linears[ii]->size = ns_config->num_v_vars;
		vec_mp_to_d(v_linears[ii], ns_config->v_linears[ii]);
	}
	
	num_additional_linears = ns_config->num_additional_linears;
	additional_linears_terminal = (vec_d *) br_malloc(ns_config->num_additional_linears*sizeof(vec_d));
	additional_linears_starting = (vec_d *) br_malloc(ns_config->num_additional_linears*sizeof(vec_d));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_d(additional_linears_terminal[ii], ns_config->num_natural_vars+ns_config->num_synth_vars);
		additional_linears_terminal[ii]->size = ns_config->num_natural_vars+ns_config->num_synth_vars;
		vec_mp_to_d(additional_linears_terminal[ii],ns_config->additional_linears_terminal[ii]);
		
		init_vec_d(additional_linears_starting[ii], ns_config->num_natural_vars+ns_config->num_synth_vars);
		additional_linears_starting[ii]->size = ns_config->num_natural_vars+ns_config->num_synth_vars;
		vec_mp_to_d(additional_linears_starting[ii],ns_config->additional_linears_starting[ii]);
	}
	
	
	
	
	if (side_ == nullspace_handedness::LEFT) {
		starting_linears = (vec_d **)br_malloc(ns_config->num_jac_equations*sizeof(vec_d *));
		
		for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
			starting_linears[ii] = (vec_d *)br_malloc(ns_config->max_degree*sizeof(vec_d));
			for (int jj=0; jj<ns_config->max_degree; jj++) {
				init_vec_d(starting_linears[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
				starting_linears[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
				
				vec_mp_to_d(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
			}
		}
	}
	else{ // right hand nullspace
		starting_linears = (vec_d **)br_malloc(ns_config->randomizer()->num_rand_funcs()*sizeof(vec_d *));
		
		for (int ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ++ii) {
			int curr_degree = ns_config->randomizer()->randomized_degree(ii)-1;
			
			starting_linears[ii] = (vec_d *)br_malloc(curr_degree*sizeof(vec_d));
			for (int jj=0; jj<curr_degree; jj++) {
				init_vec_d(starting_linears[ii][jj],ns_config->num_synth_vars+ns_config->num_natural_vars);
				starting_linears[ii][jj]->size = ns_config->num_synth_vars+ns_config->num_natural_vars;
				
				vec_mp_to_d(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
			}
		}
	}
	
	
	
	SolverDoublePrecision::setup(_SLP, ns_config->randomizer());
	
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(_SLP,ns_config, W, solve_options);
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
		
		this->SLP_derivative = this->BED_mp->SLP_derivative;
	}
	else{
		SLP_derivative = new prog_deriv_t;
		
		setup_deriv_from_SLP(SLP_derivative, this->SLP);
		while (SLP_derivative->order < 2) {
			setupNext_derivs(SLP_derivative);
		}
	}
	
	
	
	
	return SUCCESSFUL;
}




void nullspacejac_eval_data_d::print()
{
    
    std::cout << "num_jac_equations " << num_jac_equations << std::endl;
    std::cout << "target_dim " << target_dim << std::endl;
    
    std::cout << "ambient_dim " << ambient_dim << std::endl;
    std::cout << "target_crit_codim " << target_crit_codim << std::endl;
    
    std::cout << "num_natural_vars " << num_natural_vars << std::endl;
    std::cout << "num_synth_vars " << num_synth_vars << std::endl;
    std::cout << "num_v_vars " << num_v_vars << std::endl;
	
    std::cout << "num_randomized_eqns " << randomizer()->num_rand_funcs() << std::endl;
    std::cout << "max_degree " << max_degree << std::endl;
	

	
	for (int ii=0; ii<num_jac_equations; ii++) {
		print_point_to_screen_matlab(v_linears[ii],"v_lin");
	}
	
	if (side_==nullspace_handedness::LEFT) {
		for (int ii=0; ii<num_jac_equations; ii++) {
			for (int jj=0; jj<max_degree; jj++) {
				print_point_to_screen_matlab(starting_linears[ii][jj],"starting_lin");
			}
		}
	}
	else{
		for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
			for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
				print_point_to_screen_matlab(starting_linears[ii][jj],"starting_lin");
			}
		}
	}
	
	
	std::cout << randomizer()->num_rand_funcs() << " " << num_additional_linears << " " << num_jac_equations << " " << patch.num_patches << std::endl;
	
	
	std::cout << *randomizer() << std::endl;
}


///////////////
//
//   end nullspace_eval_data_d
//
/////////////






















int nullspacejac_solver_master_entry_point(int							MPType,
										   WitnessSet					&W, // carries with it the start points, and the linears.
										   SolverOutput				& solve_out, // new data goes in here
										   NullspaceConfiguration				*ns_config,
										   SolverConfiguration			& solve_options)
{
	
    
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(NULLSPACE);
	}
	
	
	
	
	
	prog_t SLP;
	//	// setup a straight-line program, using the file(s) created by the parser.  the input file must already be parsed
	setupProg(&SLP, solve_options.T.Precision, solve_options.T.MPType);
	
	
	solve_options.T.numVars = W.num_variables();
	
	nullspacejac_eval_data_d *ED_d = NULL;
	nullspacejac_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new nullspacejac_eval_data_d(0);
			
			ED_d->setup(&SLP,
						ns_config,
						W,
						solve_options);
			break;
			
		case 1:
			ED_mp = new nullspacejac_eval_data_mp(1);
			
			ED_mp->setup(&SLP,
						 ns_config,
						 W,
						 solve_options);
			break;
		case 2:
			ED_d = new nullspacejac_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(&SLP,
						ns_config,
						W,
						solve_options);
			
			
			
			adjust_tracker_AMP(& (solve_options.T), W.num_variables()); // & initialize latest_newton_residual_mp
			
			break;
		default:
			break;
	}
	
	
	
	master_solver(solve_out, W,
                  ED_d, ED_mp,
                  solve_options);
	
	
	
	switch (solve_options.T.MPType) {
		case 0:
			delete ED_d;
			break;
			
		case 1:
			delete ED_mp;
			break;
		case 2:
			delete ED_d;
			break;
		default:
			br_exit(396);
			break;
	}
	
	
	clearProg(&SLP, solve_options.T.MPType, 1); // 1 means call freeprogeval()
	return SUCCESSFUL;
	
}






void nullspace_slave_entry_point(SolverConfiguration & solve_options)
{
	
	
	// already received the flag which indicated that this worker is going to be performing the nullspace calculation.
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	int *settings_buffer = (int *) br_malloc(2*sizeof(int));
	MPI_Bcast(settings_buffer,2,MPI_INT, 0,solve_options.comm());
	solve_options.robust = settings_buffer[0];
	solve_options.use_gamma_trick = settings_buffer[1];
	free(settings_buffer);
	
	
	nullspacejac_eval_data_d *ED_d = NULL;
	nullspacejac_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new nullspacejac_eval_data_d(0);
			ED_d->receive(solve_options);
			break;
			
		case 1:
			ED_mp = new nullspacejac_eval_data_mp(1);
			ED_mp->receive(solve_options);
			
			break;
		case 2:
			ED_d = new nullspacejac_eval_data_d(2);
			ED_mp = ED_d->BED_mp;
			ED_d->receive(solve_options);
			break;
		default:
			break;
	}
	
    
	
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "nullspace_left_output",
                        &midOUT, "nullspace_left_midpath_data");
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	

	worker_tracker_loop(&trackCount, OUT, midOUT,
						ED_d, ED_mp,
						solve_options);
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	switch (solve_options.T.MPType) {
		case 0:
			delete ED_d;
			break;
			
		case 1:
			delete ED_mp;
			break;
            
		case 2:
			delete ED_d;
			
		default:
			break;
	}
	
	//clear data
}






int nullspacejac_right_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	
	BED->temp_vars.ensure_have_scalars(6);
	BED->temp_vars.ensure_have_vectors(13);
	BED->temp_vars.ensure_have_matrices(10);
	
	
	int offset;
	comp_d *one_minus_s = &BED->temp_vars.scalars[0], *gamma_s = &BED->temp_vars.scalars[1];
	
	set_one_d(*one_minus_s);
	sub_d(*one_minus_s, *one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_d(*gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	
	comp_d *running_prod = &BED->temp_vars.scalars[2];
	comp_d *temp = &BED->temp_vars.scalars[3], *temp2 = &BED->temp_vars.scalars[4], *temp3 = &BED->temp_vars.scalars[5];
	
	
	
	
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	vec_d *curr_x_vars = &BED->temp_vars.vectors[0];
	increase_size_vec_d(*curr_x_vars, BED->num_natural_vars);
	(*curr_x_vars)->size = BED->num_natural_vars;
	for (int ii=0; ii<BED->num_natural_vars; ii++)
		set_d(&(*curr_x_vars)->coord[ii], &current_variable_values->coord[ii]);
	
	vec_d *curr_v_vars = &BED->temp_vars.vectors[1];
	increase_size_vec_d(*curr_v_vars, BED->num_v_vars);
	(*curr_v_vars)->size = BED->num_v_vars;
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_d(&(*curr_v_vars)->coord[ii], &current_variable_values->coord[ii+BED->num_natural_vars]);
	
	
	
	
	vec_d *patchValues = &BED->temp_vars.vectors[2];
	vec_d *temp_function_values = &BED->temp_vars.vectors[3];
	
	
	vec_d *AtimesF = &BED->temp_vars.vectors[4];
	vec_d *linprod_x = &BED->temp_vars.vectors[5];  change_size_vec_d(*linprod_x, BED->num_jac_equations);
	(*linprod_x)->size = BED->num_jac_equations;
	
	vec_d *linprod_times_gamma_s = &BED->temp_vars.vectors[6]; change_size_vec_d(*linprod_times_gamma_s,BED->num_jac_equations);
	(*linprod_times_gamma_s)->size = BED->num_jac_equations;
	vec_d *tempvec = &BED->temp_vars.vectors[7];
	vec_d *tempvec2 = &BED->temp_vars.vectors[8];
	
	vec_d *target_function_values = &BED->temp_vars.vectors[9];
	vec_d *target_function_values_times_oneminus_s = &BED->temp_vars.vectors[10];
	
	vec_d *start_function_values = &BED->temp_vars.vectors[11];
	increase_size_vec_d(*start_function_values,BED->num_jac_equations); (*start_function_values)->size = BED->num_jac_equations;
	
	
	
	vec_d *deriv_vals = &BED->temp_vars.vectors[12];
	
	
	
	
	
	
	mat_d *Jv_Patch = &BED->temp_vars.matrices[0];
	
	mat_d *tempmat5 = &BED->temp_vars.matrices[1];
	
	
	
	mat_d *lin_func_vals = &BED->temp_vars.matrices[2]; change_size_mat_d(*lin_func_vals,BED->randomizer()->num_rand_funcs(), BED->max_degree);
	(*lin_func_vals)->rows = BED->randomizer()->num_rand_funcs(); (*lin_func_vals)->cols = BED->max_degree;
	//this is an overallocation, in that not each has the same number of entries in it, but it's easier this way.
	
	
	mat_d *AtimesJ = &BED->temp_vars.matrices[3];
	
	
	mat_d *linprod_derivative_wrt_x = &BED->temp_vars.matrices[4];
	increase_size_mat_d(*linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_natural_vars);
	(*linprod_derivative_wrt_x)->rows = BED->num_jac_equations; (*linprod_derivative_wrt_x)->cols = BED->num_natural_vars;
	
	
	
	mat_d* tempmat4 = &BED->temp_vars.matrices[5];
	
	// create temp matrices
	mat_d *tempmat1 = &BED->temp_vars.matrices[6];
	mat_d *tempmat2 = &BED->temp_vars.matrices[7];
	mat_d *tempmat3 = &BED->temp_vars.matrices[8];
	
	// these resizes should be moved to just before the matrices are used.
	increase_size_mat_d(*tempmat1,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat1)->rows = BED->num_jac_equations; (*tempmat1)->cols = BED->num_v_vars;
	
	increase_size_mat_d(*tempmat2,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat2)->rows = BED->num_jac_equations; (*tempmat2)->cols = BED->num_v_vars;
	
	
	mat_d *temp_jacobian_functions = &BED->temp_vars.matrices[9];
	
	
	
	
	
	
	// the main evaluations for $x$
	evalDeriv_d(*temp_function_values, *deriv_vals, *tempvec, *tempvec2, *curr_x_vars, BED->SLP_derivative);
	//tempvec, and tempvec2, are throwaway variables used to capture output from this which we don't care about -- the linears added to the system.  perhaps we could optimize this?
	
	
	// evaluate the patch
	patch_eval_d(    *patchValues, parVals, parDer, *Jv_Patch, Jp, *curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	
	
	
	
	
	
	
	
	//resize output variables to correct size
	increase_size_vec_d(funcVals,BED->num_variables);
	increase_size_mat_d(Jv, BED->num_variables, BED->num_variables);
	increase_size_mat_d(Jp, BED->num_variables, 1);
	funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
	Jv->cols = BED->num_variables;  //  <-- this must be square
	Jp->cols = 1; // this is a column vector, or nx1 matrix, for this evaluator
	
	
	
	//////
	// initialize stuff to all 0's
	///////
	
	for (int ii=0; ii<Jv->rows; ii++){
		for (int jj=0; jj<Jv->cols; jj++){
			set_zero_d(&Jv->entry[ii][jj]);
			// initialize entire matrix to 0 // this is kinda wasteful, but safe.  wasteful because most entries are set in this function, but safe because we only have to set entries which we know aren't 0, and rest are 0 by default.
		}
	}
	
	for (int ii = 0; ii<BED->num_variables; ii++)
		set_zero_d(&Jp->entry[ii][0]);  // initialize entire matrix to 0.  same comment regarding wastefulness.
	
	
	
	
	
	///////////
	// orig eqns
	/////////////
	
	
	// resize temp_jacobians to be Nxn, to hold the first derivatives.
	change_size_mat_d(*temp_jacobian_functions, (*temp_function_values)->size, (BED->num_natural_vars+BED->num_synth_vars));
	(*temp_jacobian_functions)->rows = (*temp_function_values)->size;
	(*temp_jacobian_functions)->cols = (BED->num_natural_vars+BED->num_synth_vars);
	
	
	// unpack the first derivatives
	int derivative_offset = 0;
	for (int ii = 0; ii<(*temp_function_values)->size; ii++) {
		for (int jj = 0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			set_d(&(*temp_jacobian_functions)->entry[ii][jj],&(*deriv_vals)->coord[derivative_offset]);
			derivative_offset++; // the end of this loop will leave this offset for use later.
		}
	} // the remainer of deriv_vals contains the second derivatives, which will be used later in this function
	
	
	
	
	
	// randomize
	BED->randomizer()->randomize(*AtimesF,*AtimesJ,*temp_function_values,*temp_jacobian_functions, &(*curr_x_vars)->coord[0]);
	
	for (int ii=0; ii<(*AtimesF)->size; ii++)  // for each function, after (real) randomization
		set_d(&funcVals->coord[ii], &(*AtimesF)->coord[ii]);
	
	
	
	
	
	//////////////
	//  orig function jacobian values
	////////////////
	
	
	
	
	// set the jacobian equations for orig into Jv
	
	// copy the jacobian into the return value for the evaluator
	for (int ii=0; ii< (*AtimesJ)->rows; ii++) // for every function
		for (int jj=0; jj< (*AtimesJ)->cols; jj++)
			set_d(&Jv->entry[ii][jj],&(*AtimesJ)->entry[ii][jj]);
	
	
	
	
	
	
	
	//////////
	// the additional linears.  there are $r-\ell$ of them.
	///////////
	
	
	
	offset = BED->randomizer()->num_rand_funcs();
	for (int ii=0; ii< BED->num_additional_linears; ii++) {
		
		dot_product_d(*temp, BED->additional_linears_terminal[ii], *curr_x_vars); //temp = terminal(x)
		neg_d(&Jp->entry[offset+ii][0], *temp);  // Jp = -terminal(x)
		mul_d(*temp3, *temp, *one_minus_s);        // temp3 = (1-s)*terminal(x)
		
		dot_product_d(*temp,  BED->additional_linears_starting[ii], *curr_x_vars); // temp = starting(x)
		mul_d(*temp2, BED->gamma, *temp);         // temp2 = gamma*starting(x)
		
		add_d(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], *temp2);   // Jp = -terminal + gamma*start
		
		mul_d(*temp2, *gamma_s, *temp);            // temp2 = gamma*s*starting(x)
		
		add_d(&funcVals->coord[offset+ii],*temp2, *temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (int jj=0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			mul_d(*temp, *gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]);
			mul_d(*temp2,*one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_d(&Jv->entry[ii+offset][jj], *temp, *temp2);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	

	
	
	
	
	
	
	
	
	
	
	// done to here.  work below.
	
	
	
	/////////////////////////////
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	////////////////////////////////
	

	
	increase_size_mat_d(BED->jac_with_proj, BED->randomizer()->num_rand_funcs()+BED->num_projections, BED->num_natural_vars-1);
	BED->jac_with_proj->rows = BED->randomizer()->num_rand_funcs()+BED->num_projections;
	BED->jac_with_proj->cols = BED->num_natural_vars-1;
	
	
	
	// copy into beginning of a temp matrix stored in BED
	for (int ii=0; ii< BED->randomizer()->num_rand_funcs(); ii++)// for every function
		for (int jj=1; jj<BED->num_natural_vars; jj++) // for only the natural variables, omitting the hom var and synth_vars
			set_d(&BED->jac_with_proj->entry[ii][jj-1], &(*AtimesJ)->entry[ii][jj]);
	
	
	// concatenate in the values of the projections, homogenized, to get [Jv | pi]
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=1; jj<BED->num_natural_vars; jj++) {
			set_d(&BED->jac_with_proj->entry[ii+BED->randomizer()->num_rand_funcs()][jj-1],&BED->target_projection[ii]->coord[jj]);
			// set, to  below of the jacobian matrix
		}
	}
	
	
	
	
	//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jac_hom_w_proj");
	mul_mat_vec_d(*target_function_values, BED->jac_with_proj, (*curr_v_vars));
	
	
	
	//	print_point_to_screen_matlab(target_function_values,"target_f");
	vec_mulcomp_d(*target_function_values_times_oneminus_s, *target_function_values, *one_minus_s);
	
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	
	// the product of the linears
	
	for (int jj=0; jj<BED->randomizer()->num_rand_funcs(); jj++) { // for the entries down below which depend on x.
		
		//perform the $x$ evaluation
		set_one_d(&(*linprod_x)->coord[jj]); // initialize to 1 for multiplication
		for (int ii=0; ii<BED->randomizer()->randomized_degree(jj)-1; ++ii) {
			dot_product_d(&(*lin_func_vals)->entry[jj][ii],BED->starting_linears[jj][ii],*curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_d(&(*linprod_x)->coord[jj], &(*linprod_x)->coord[jj], &(*lin_func_vals)->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_d(*temp, BED->v_linears[jj], (*curr_v_vars));
		
		//now set the combined $x,v$ value
		mul_d(&(*start_function_values)->coord[jj], &(*linprod_x)->coord[jj], *temp); // start_function_values = linprod_x(x)*linear(v)
	}
	
	int asdf = BED->randomizer()->num_rand_funcs();
	for (int jj=0; jj<BED->num_projections; jj++){
		dot_product_d(&(*start_function_values)->coord[asdf+jj], BED->v_linears[jj+asdf], *curr_v_vars);
	}
	//done creating the start function values for jacobian equations
	
	vec_mulcomp_d(*linprod_times_gamma_s, *start_function_values, *gamma_s); // sets the value linprod_x*gamma*s
	
	
	
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int jj=0; jj<BED->num_jac_equations; jj++) {
		
		neg_d( &Jp->entry[jj + offset][0], &(*target_function_values)->coord[jj]);  // Jp = -target
		mul_d( *temp, BED->gamma, &(*start_function_values)->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_d( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], *temp); // Jp = Jp + gamma*start = -target + gamma*start
	}
	// mix the start and target into final function value to return.
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		add_d(&funcVals->coord[ii+offset], &(*target_function_values_times_oneminus_s)->coord[ii], &(*linprod_times_gamma_s)->coord[ii]);
	}
	
	
	// done with non-patch function values and Jp
	
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE start and target WRT V
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->randomizer()->num_rand_funcs(); ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_d(*temp, &BED->v_linears[ii]->coord[jj], &(*linprod_x)->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_d(*temp2, *temp, *gamma_s);                                         // temp2 = gamma*s*M_ij * linprod_x(x)
			
			mul_d(*temp, *one_minus_s, &BED->jac_with_proj->entry[ii][jj]);        // temp = (1-s)*([Jf; Jpi])_ij
			
			add_d(&Jv->entry[ii+offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2);       // Jv = temp + temp2
		}
	}
	
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_d(*temp, &BED->target_projection[ii]->coord[jj+1],*one_minus_s); // moving to the target projection
			mul_d(*temp2, &BED->v_linears[ii+BED->randomizer()->num_rand_funcs()]->coord[jj], *gamma_s); // moving away from teh starting linear
			add_d(&Jv->entry[ii+BED->randomizer()->num_rand_funcs() + offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2); //(1-s)*pi_i(jj) + gamma*s*M_ij
		}
	}
	
	
	
	
	//////////////////
	// now the x derivatives corresponding to the linprod start system.  these will be used after the next block
	//////////////////////
	
	
	
	// an implementation of the product rule
	for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
		
		for (int kk=0; kk<(BED->num_natural_vars+BED->num_synth_vars); kk++) { // for each variable
			set_zero_d(&(*linprod_derivative_wrt_x)->entry[mm][kk]); // initialize to 0 for the sum
			
			int curr_degree = BED->randomizer()->randomized_degree(mm)-1;
			for (int ii=0; ii<curr_degree; ++ii) { //  for each linear
				
				set_d(*running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (int jj=0; jj<curr_degree; jj++) {
					if (jj!=ii) {
						mul_d(*running_prod,*running_prod,&(*lin_func_vals)->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_d(&(*linprod_derivative_wrt_x)->entry[mm][kk],&(*linprod_derivative_wrt_x)->entry[mm][kk],*running_prod);
				
			}// re:ii
			
			dot_product_d(*temp, BED->v_linears[mm], (*curr_v_vars));  // these two lines multiply by  (v_linear •	v)
			mul_d(&(*linprod_derivative_wrt_x)->entry[mm][kk], &(*linprod_derivative_wrt_x)->entry[mm][kk], *temp);
		} // re: kk
		
	} // re: mm

	
	
	
	
	
	
	
	
	
	
	
	
	
	///////////////
	// DIFFERENTIATE THE derivative of the target jacobian system wrt $h,x_i$.
	/////////////////
	
	
	//adjust the offset to remain vertically correct
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	
	
	//only need to do this size change once.
	increase_size_mat_d(*tempmat1,BED->randomizer()->num_base_funcs(),BED->num_natural_vars+BED->num_synth_vars-1);
	(*tempmat1)->rows = BED->randomizer()->num_base_funcs();
	(*tempmat1)->cols = BED->num_natural_vars+BED->num_synth_vars-1;
	
	
	// we'll pack the second derivatives into tempmat1, multiply by the randomization matrix, and then multiply by the v variables
	
	
	
	// you only need the following homogenizing randomization matrices if the system is not square
	if (!BED->randomizer()->is_square()) { // homogenize and randomize
		mat_cp_d(*tempmat2, *(BED->randomizer()->get_mat_d()) );
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			for (int nn=0; nn<BED->randomizer()->num_base_funcs(); nn++) {
				for (int pp=0; pp<BED->randomizer()->deficiency(mm,nn); pp++) {
					// homogenize
					mul_d(&( (*tempmat2)->entry[mm][nn]), &( (*curr_x_vars)->coord[0]),&( (*tempmat2)->entry[mm][nn]));
				}
			}
		}
		
		
		
		mat_cp_d(*tempmat4, *(BED->randomizer()->get_mat_d()) );
		
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			for (int nn=0; nn<BED->randomizer()->num_base_funcs(); nn++) {
				
				int degree_deficiency = BED->randomizer()->deficiency(mm,nn);
				
				for (int pp=0; pp<degree_deficiency-1; pp++) {  // -1 here because this is for the product rule
					// homogenize
					mul_d(&( (*tempmat4)->entry[mm][nn]), &( (*curr_x_vars)->coord[0]),&( (*tempmat4)->entry[mm][nn]));
				}
				
				set_zero_d(*temp);
				(*temp)->r = degree_deficiency;  //  n * h^{n-1}
				
				mul_d(&( (*tempmat4)->entry[mm][nn]), *temp ,&( (*tempmat4)->entry[mm][nn]));
				
				
			}
		}
		
		for (int mm=0; mm<BED->randomizer()->num_base_funcs(); mm++) {
			for (int nn=1; nn<BED->num_natural_vars; nn++) {
				set_d(&(*temp_jacobian_functions)->entry[mm][nn-1],&(*temp_jacobian_functions)->entry[mm][nn]);
			}
		}
		(*temp_jacobian_functions)->cols = BED->num_natural_vars-1;
		
	}
	
	
	
	
	
	
	
	// the hom variable
	int current_variable_index = 0;
	int local_offset = derivative_offset-1; // initialize.  we count up from here
	for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each base function
		int entry_counter = 0;
		
		for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
			for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
				local_offset++;//  MUST always increment.
				
				if (jj==0) {
					continue;
				}
				
				if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
					continue;
				}

				set_d(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);

				entry_counter++;
			}
		}
	} // re: making the interior matrix of second derivatives
	
//	print_matrix_to_screen_matlab(*tempmat1,"second_derivatives");
	//need to homogenize properly, with the product rule
	

	// ok, now tempmat1 has in it all the (non-homogenized) derivatives of the derivatives of the system f.

	if (!BED->randomizer()->is_square()) { // homogenize and randomize
		mat_mul_d(*tempmat3, *tempmat2, *tempmat1); // randomize
		
		
		// not the other part for the product rule with the homogenizing variable
		mat_mul_d(*tempmat5, *tempmat4, *temp_jacobian_functions);  // n*h^{n-1} * f_{bla}
		
		
		
		
		mul_mat_vec_d(*tempvec, *tempmat3, *curr_v_vars);
		mul_mat_vec_d(*tempvec2, *tempmat5, *curr_v_vars);
		
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			add_d(&(*tempvec2)->coord[mm], &(*tempvec)->coord[mm], &(*tempvec2)->coord[mm]);
		}

	}
	else
	{
		mul_mat_vec_d(*tempvec2, *tempmat1, (*curr_v_vars));
	}
	
	
	
	vec_mulcomp_d(*tempvec, *tempvec2, *one_minus_s);// now tempvec has the derivatives wrt $x$. (for all jac eqns)
	
	// combine this derivative and the linprod derivative to get the values for the $x$ portion of the jacobian matrix to return
	for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
		mul_d(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
		add_d(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
	}
	
	
	// the remainder are all 0, cuz just have v in them.  this is probably unnecessary because the matrix was zero'd earlier.  (wasn't it?)
	for (int mm=0; mm<BED->num_projections; mm++) {
		set_zero_d(&Jv->entry[mm+offset+BED->randomizer()->num_rand_funcs()][current_variable_index]);
	}
	//////  end the calculation for the hom var
	
	
	
	
	
	
	
	
	
	
	
	// THE X VARS (non-homogenizing vars)
	
	// start this loop at 1, because we computed the homogenizing variable differently (product rule)
	for (current_variable_index=1; current_variable_index<BED->num_natural_vars+BED->num_synth_vars; current_variable_index++) {
//		std::cout << "variable " << current_variable_index << std::endl;
		
		
		int local_offset = derivative_offset-1; // initialize this counter.  we count up from here.  must be -1 cause ended -1 when unpacking first derivatives earlier
		for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each base function
			
			int entry_counter = 0;
			
			
			for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
				for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
					local_offset++;//  MUST always increment.
					
					
					if (ii==0) { // skip the mixed derivatives involving the hom var // they all come first
						continue;
					}
					
					if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
						continue;
					}
					
					// made it here, so this entry in deriv_vals is indeed one of the second derivatives we need
					set_d(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
					
					entry_counter++;
				}
			}
		} // re: making the interior matrix of second derivatives
		
		
		
		
		
		// ok, now tempmat1 has in it all the (non-homogenized) derivatives of the derivatives of the system f.

		if (!BED->randomizer()->is_square()) { // homogenize and randomize
			mat_mul_d(*tempmat3, *tempmat2, *tempmat1); // randomize
			
			mul_mat_vec_d(*tempvec2, *tempmat3, (*curr_v_vars)); // recall v has been truncated to avoid unnecessary 0*number multiplication
		}
		else
		{
			mul_mat_vec_d(*tempvec2, *tempmat1, (*curr_v_vars)); // recall v has been truncated to avoid unnecessary 0*number multiplication
		}
		
		
		vec_mulcomp_d(*tempvec, *tempvec2, *one_minus_s);// now tempvec has the derivatives wrt $x$. (for all jac eqns)
		
		// combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			mul_d(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
			add_d(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
			//  (1-s)*R ∂/∂x Jf + gamma s ∂/∂x linprod
		}
	
		
		
		// the remainder are all 0, cuz just have v in them.  this is probably unnecessary because the matrix was zero'd earlier.  (wasn't it?)
		for (int mm=0; mm<BED->num_projections; mm++) {
			set_zero_d(&Jv->entry[mm+offset+BED->randomizer()->num_rand_funcs()][current_variable_index]);
		}
	}//re: current_variable_index for diff
	
	// END DIFF wrt x
	
	

	
	//set the X PATCH values
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations;
	
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_d(&funcVals->coord[ii+offset], &(*patchValues)->coord[ii]);
		
		// Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &(*Jv_Patch)->entry[ii][jj]);
	}
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 2151;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		throw std::logic_error("logic error thrown");
	}
	
	// V patch
	set_one_d(*temp2);
	dot_product_d(*temp, BED->v_patch, (*curr_v_vars));
	sub_d(&funcVals->coord[BED->num_variables-1], *temp, *temp2);  // f = patch*v-1
	
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_d(&Jv->entry[BED->num_variables-1][BED->num_natural_vars+ii], &BED->v_patch->coord[ii]);
	
	
	
	// finally, set parVals & parDer correctly
	
	change_size_point_d(parVals, 1);  change_size_vec_d(parDer, 1);
	parVals->size = parDer->size = 1;
	
	set_d(&parVals->coord[0], pathVars); // s = t
	set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	if (BED->verbose_level()==10 || BED->verbose_level()==-10) {
		printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		
		//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"R");
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
		
		print_point_to_screen_matlab(*curr_x_vars,"x1");
		print_point_to_screen_matlab(*curr_v_vars,"V1");
		print_point_to_screen_matlab(funcVals,"F");
		print_matrix_to_screen_matlab(Jv,"Jv");
		print_matrix_to_screen_matlab(Jp,"Jp");
		//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
		//			//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
		//		print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");
		//	std::cout << "\n\n**************\n\n";
		if (BED->verbose_level()<0) {
			mypause();
		}
		
	}
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_x_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	//	printf("exiting eval\n");
	return 0;
}




int nullspacejac_right_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	BED->temp_vars.ensure_have_scalars(6);
	BED->temp_vars.ensure_have_vectors(13);
	BED->temp_vars.ensure_have_matrices(10);
	
	
	int offset;
	comp_mp *one_minus_s = &BED->temp_vars.scalars[0], *gamma_s = &BED->temp_vars.scalars[1];
	
	set_one_mp(*one_minus_s);
	sub_mp(*one_minus_s, *one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_mp(*gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	
	comp_mp *running_prod = &BED->temp_vars.scalars[2];
	comp_mp *temp = &BED->temp_vars.scalars[3], *temp2 = &BED->temp_vars.scalars[4], *temp3 = &BED->temp_vars.scalars[5];
	
	
	
	
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	vec_mp *curr_x_vars = &BED->temp_vars.vectors[0];
	increase_size_vec_mp(*curr_x_vars, BED->num_natural_vars);
	(*curr_x_vars)->size = BED->num_natural_vars;
	for (int ii=0; ii<BED->num_natural_vars; ii++)
		set_mp(&(*curr_x_vars)->coord[ii], &current_variable_values->coord[ii]);
	
	vec_mp *curr_v_vars = &BED->temp_vars.vectors[1];
	increase_size_vec_mp(*curr_v_vars, BED->num_v_vars);
	(*curr_v_vars)->size = BED->num_v_vars;
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&(*curr_v_vars)->coord[ii], &current_variable_values->coord[ii+BED->num_natural_vars]);
	
	
	
	
	vec_mp *patchValues = &BED->temp_vars.vectors[2];
	vec_mp *temp_function_values = &BED->temp_vars.vectors[3];
	
	
	vec_mp *AtimesF = &BED->temp_vars.vectors[4];
	vec_mp *linprod_x = &BED->temp_vars.vectors[5];  change_size_vec_mp(*linprod_x, BED->num_jac_equations);
	(*linprod_x)->size = BED->num_jac_equations;
	
	vec_mp *linprod_times_gamma_s = &BED->temp_vars.vectors[6]; change_size_vec_mp(*linprod_times_gamma_s,BED->num_jac_equations);
	(*linprod_times_gamma_s)->size = BED->num_jac_equations;
	vec_mp *tempvec = &BED->temp_vars.vectors[7];
	vec_mp *tempvec2 = &BED->temp_vars.vectors[8];
	
	vec_mp *target_function_values = &BED->temp_vars.vectors[9];
	vec_mp *target_function_values_times_oneminus_s = &BED->temp_vars.vectors[10];
	
	vec_mp *start_function_values = &BED->temp_vars.vectors[11];
	increase_size_vec_mp(*start_function_values,BED->num_jac_equations); (*start_function_values)->size = BED->num_jac_equations;
	
	
	
	vec_mp *deriv_vals = &BED->temp_vars.vectors[12];
	
	
	
	
	
	
	mat_mp *Jv_Patch = &BED->temp_vars.matrices[0];
	mat_mp *tempmat5 = &BED->temp_vars.matrices[1];
	
	mat_mp *lin_func_vals = &BED->temp_vars.matrices[2]; change_size_mat_mp(*lin_func_vals,BED->randomizer()->num_rand_funcs(), BED->max_degree);
	(*lin_func_vals)->rows = BED->randomizer()->num_rand_funcs(); (*lin_func_vals)->cols = BED->max_degree;
	//this is an overallocation, in that not each has the same number of entries in it, but it's easier this way.
	
	
	mat_mp *AtimesJ = &BED->temp_vars.matrices[3];
	
	
	mat_mp *linprod_derivative_wrt_x = &BED->temp_vars.matrices[4];
	increase_size_mat_mp(*linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_natural_vars);
	(*linprod_derivative_wrt_x)->rows = BED->num_jac_equations; (*linprod_derivative_wrt_x)->cols = BED->num_natural_vars;
	
	
	
	
	mat_mp* tempmat4 = &BED->temp_vars.matrices[5];
	
	
	// create temp matrices
	mat_mp *tempmat1 = &BED->temp_vars.matrices[6];
	mat_mp *tempmat2 = &BED->temp_vars.matrices[7];
	mat_mp *tempmat3 = &BED->temp_vars.matrices[8];
	
	increase_size_mat_mp(*tempmat1,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat1)->rows = BED->num_jac_equations; (*tempmat1)->cols = BED->num_v_vars;
	
	increase_size_mat_mp(*tempmat2,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat2)->rows = BED->num_jac_equations; (*tempmat2)->cols = BED->num_v_vars;
	
	
	mat_mp *temp_jacobian_functions = &BED->temp_vars.matrices[9];
	
	
	
	
	
	
	// the main evaluations for $x$
	evalDeriv_mp(*temp_function_values, *deriv_vals, *tempvec, *tempvec2, *curr_x_vars, BED->SLP_derivative);
	//tempvec, and tempvec2, are throwaway variables used to capture output from this which we don't care about -- the linears added to the system.  perhaps we could optimize this?
	
	
	// evaluate the patch
	patch_eval_mp(    *patchValues, parVals, parDer, *Jv_Patch, Jp, *curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	
	
	
	
	
	
	
	
	//resize output variables to correct size
	increase_size_vec_mp(funcVals,BED->num_variables);
	increase_size_mat_mp(Jv, BED->num_variables, BED->num_variables);
	increase_size_mat_mp(Jp, BED->num_variables, 1);
	funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
	Jv->cols = BED->num_variables;  //  <-- this must be square
	Jp->cols = 1; // this is a column vector, or nx1 matrix, for this evaluator
	
	
	
	//////
	// initialize stuff to all 0's
	///////
	
	for (int ii=0; ii<Jv->rows; ii++){
		for (int jj=0; jj<Jv->cols; jj++){
			set_zero_mp(&Jv->entry[ii][jj]);
			// initialize entire matrix to 0 // this is kinda wasteful, but safe.  wasteful because most entries are set in this function, but safe because we only have to set entries which we know aren't 0, and rest are 0 by default.
		}
	}
	
	for (int ii = 0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);  // initialize entire matrix to 0.  same comment regarding wastefulness.
	
	
	
	
	
	///////////
	// orig eqns
	/////////////
	
	
	// resize temp_jacobians to be Nxn, to hold the first derivatives.
	change_size_mat_mp(*temp_jacobian_functions, (*temp_function_values)->size, (BED->num_natural_vars+BED->num_synth_vars));
	(*temp_jacobian_functions)->rows = (*temp_function_values)->size;
	(*temp_jacobian_functions)->cols = (BED->num_natural_vars+BED->num_synth_vars);
	
	
	// unpack the first derivatives
	int derivative_offset = 0;
	for (int ii = 0; ii<(*temp_function_values)->size; ii++) {
		for (int jj = 0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			set_mp(&(*temp_jacobian_functions)->entry[ii][jj],&(*deriv_vals)->coord[derivative_offset]);
			derivative_offset++; // the end of this loop will leave this offset for use later.
		}
	} // the remainer of deriv_vals contains the second derivatives, which will be used later in this function
	
	
	
	
	
	// randomize
	BED->randomizer()->randomize(*AtimesF,*AtimesJ,*temp_function_values,*temp_jacobian_functions, &(*curr_x_vars)->coord[0]);
	
	for (int ii=0; ii<(*AtimesF)->size; ii++)  // for each function, after (real) randomization
		set_mp(&funcVals->coord[ii], &(*AtimesF)->coord[ii]);
	
	
	
	
	
	//////////////
	//  orig function jacobian values
	////////////////
	
	
	
	
	// set the jacobian equations for orig into Jv
	
	// copy the jacobian into the return value for the evaluator
	for (int ii=0; ii< (*AtimesJ)->rows; ii++) // for every function
		for (int jj=0; jj< (*AtimesJ)->cols; jj++)
			set_mp(&Jv->entry[ii][jj],&(*AtimesJ)->entry[ii][jj]);
	
	
	
	
	
	
	
	//////////
	// the additional linears.  there are $r-\ell$ of them.
	///////////
	
	
	
	offset = BED->randomizer()->num_rand_funcs();
	for (int ii=0; ii< BED->num_additional_linears; ii++) {
		
		dot_product_mp(*temp, BED->additional_linears_terminal[ii], *curr_x_vars); //temp = terminal(x)
		neg_mp(&Jp->entry[offset+ii][0], *temp);  // Jp = -terminal(x)
		mul_mp(*temp3, *temp, *one_minus_s);        // temp3 = (1-s)*terminal(x)
		
		dot_product_mp(*temp,  BED->additional_linears_starting[ii], *curr_x_vars); // temp = starting(x)
		mul_mp(*temp2, BED->gamma, *temp);         // temp2 = gamma*starting(x)
		
		add_mp(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], *temp2);   // Jp = -terminal + gamma*start
		
		mul_mp(*temp2, *gamma_s, *temp);            // temp2 = gamma*s*starting(x)
		
		add_mp(&funcVals->coord[offset+ii],*temp2, *temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (int jj=0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			mul_mp(*temp, *gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]);
			mul_mp(*temp2,*one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_mp(&Jv->entry[ii+offset][jj], *temp, *temp2);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// done to here.  work below.
	
	
	
	/////////////////////////////
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	////////////////////////////////
	
	
	
	increase_size_mat_mp(BED->jac_with_proj, BED->randomizer()->num_rand_funcs()+BED->num_projections, BED->num_natural_vars-1);
	BED->jac_with_proj->rows = BED->randomizer()->num_rand_funcs()+BED->num_projections;
	BED->jac_with_proj->cols = BED->num_natural_vars-1;
	
	
	
	// copy into beginning super duper temp matrix stored in BED
	for (int ii=0; ii< BED->randomizer()->num_rand_funcs(); ii++)// for every function
		for (int jj=1; jj<BED->num_natural_vars; jj++) // for only the natural variables, omitting the hom var and synth_vars
			set_mp(&BED->jac_with_proj->entry[ii][jj-1], &(*AtimesJ)->entry[ii][jj]);
	
	
	// concatenate in the values of the projections, homogenized, to get [Jv | pi]
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=1; jj<BED->num_natural_vars; jj++) {
			set_mp(&BED->jac_with_proj->entry[ii+BED->randomizer()->num_rand_funcs()][jj-1],&BED->target_projection[ii]->coord[jj]);
			// set, to  below of the jacobian matrix
		}
	}
	
	
	
	
	//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jac_hom_w_proj");
	mul_mat_vec_mp(*target_function_values, BED->jac_with_proj, (*curr_v_vars));
	
	
	
	//	print_point_to_screen_matlab(target_function_values,"target_f");
	vec_mulcomp_mp(*target_function_values_times_oneminus_s, *target_function_values, *one_minus_s);
	
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	
	// the product of the linears
	
	for (int jj=0; jj<BED->randomizer()->num_rand_funcs(); jj++) { // for the entries down below which depend on x.
		
		//perform the $x$ evaluation
		set_one_mp(&(*linprod_x)->coord[jj]); // initialize to 1 for multiplication
		for (int ii=0; ii<BED->randomizer()->randomized_degree(jj)-1; ++ii) {
			dot_product_mp(&(*lin_func_vals)->entry[jj][ii],BED->starting_linears[jj][ii],*curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_mp(&(*linprod_x)->coord[jj], &(*linprod_x)->coord[jj], &(*lin_func_vals)->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_mp(*temp, BED->v_linears[jj], (*curr_v_vars));
		
		//now set the combined $x,v$ value
		mul_mp(&(*start_function_values)->coord[jj], &(*linprod_x)->coord[jj], *temp); // start_function_values = linprod_x(x)*linear(v)
	}
	
	int asdf = BED->randomizer()->num_rand_funcs();
	for (int jj=0; jj<BED->num_projections; jj++){
		dot_product_mp(&(*start_function_values)->coord[asdf+jj], BED->v_linears[jj+asdf], *curr_v_vars);
	}
	//done creating the start function values for jacobian equations
	
	vec_mulcomp_mp(*linprod_times_gamma_s, *start_function_values, *gamma_s); // sets the value linprod_x*gamma*s
	
	
	
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int jj=0; jj<BED->num_jac_equations; jj++) {
		
		neg_mp( &Jp->entry[jj + offset][0], &(*target_function_values)->coord[jj]);  // Jp = -target
		mul_mp( *temp, BED->gamma, &(*start_function_values)->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_mp( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], *temp); // Jp = Jp + gamma*start = -target + gamma*start
	}
	// mix the start and target into final function value to return.
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		add_mp(&funcVals->coord[ii+offset], &(*target_function_values_times_oneminus_s)->coord[ii], &(*linprod_times_gamma_s)->coord[ii]);
	}
	
	
	// done with non-patch function values and Jp
	
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE start and target WRT V
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->randomizer()->num_rand_funcs(); ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_mp(*temp, &BED->v_linears[ii]->coord[jj], &(*linprod_x)->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_mp(*temp2, *temp, *gamma_s);                                         // temp2 = gamma*s*M_ij * linprod_x(x)
			
			mul_mp(*temp, *one_minus_s, &BED->jac_with_proj->entry[ii][jj]);        // temp = (1-s)*([Jf; Jpi])_ij
			
			add_mp(&Jv->entry[ii+offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2);       // Jv = temp + temp2
		}
	}
	
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_mp(*temp, &BED->target_projection[ii]->coord[jj+1],*one_minus_s); // moving to the target projection
			mul_mp(*temp2, &BED->v_linears[ii+BED->randomizer()->num_rand_funcs()]->coord[jj], *gamma_s); // moving away from teh starting linear
			add_mp(&Jv->entry[ii+BED->randomizer()->num_rand_funcs() + offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2); //(1-s)*pi_i(jj) + gamma*s*M_ij
		}
	}
	
	
	
	
	//////////////////
	// now the x derivatives corresponding to the linprod start system.  these will be used after the next block
	//////////////////////
	
	
	
	// an implementation of the product rule
	for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
		
		for (int kk=0; kk<(BED->num_natural_vars+BED->num_synth_vars); kk++) { // for each variable
			set_zero_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk]); // initialize to 0 for the sum
			
			int curr_degree = BED->randomizer()->randomized_degree(mm)-1;
			for (int ii=0; ii<curr_degree; ++ii) { //  for each linear
				
				set_mp(*running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (int jj=0; jj<curr_degree; jj++) {
					if (jj!=ii) {
						mul_mp(*running_prod,*running_prod,&(*lin_func_vals)->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk],&(*linprod_derivative_wrt_x)->entry[mm][kk],*running_prod);
				
			}// re:ii
			
			dot_product_mp(*temp, BED->v_linears[mm], (*curr_v_vars));  // these two lines multiply by  (v_linear •	v)
			mul_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk], &(*linprod_derivative_wrt_x)->entry[mm][kk], *temp);
		} // re: kk
		
	} // re: mm
	
	
	
	
	
	
	///////////////
	// DIFFERENTIATE THE derivative of the target jacobian system wrt $h,x_i$.
	/////////////////
	
	
	//adjust the offset to remain vertically correct
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	
	
	//only need to do this size change once.
	increase_size_mat_mp(*tempmat1,BED->randomizer()->num_base_funcs(),BED->num_natural_vars+BED->num_synth_vars-1);
	(*tempmat1)->rows = BED->randomizer()->num_base_funcs();
	(*tempmat1)->cols = BED->num_natural_vars+BED->num_synth_vars-1;
	
	
	// we'll pack the second derivatives into tempmat1, multiply by the randomization matrix, and then multiply by the v variables
	
	
	
	// you only need the following homogenizing randomization matrices if the system is not square
	if (!BED->randomizer()->is_square()) { // homogenize and randomize
		mat_cp_mp(*tempmat2, *(BED->randomizer()->get_mat_mp()) );
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			for (int nn=0; nn<BED->randomizer()->num_base_funcs(); nn++) {
				for (int pp=0; pp<BED->randomizer()->deficiency(mm,nn); pp++) {
					// homogenize
					mul_mp(&( (*tempmat2)->entry[mm][nn]), &( (*curr_x_vars)->coord[0]),&( (*tempmat2)->entry[mm][nn]));
				}
			}
		}
		
		
		
		mat_cp_mp(*tempmat4, *(BED->randomizer()->get_mat_mp()) );
		
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			for (int nn=0; nn<BED->randomizer()->num_base_funcs(); nn++) {
				
				int degree_deficiency = BED->randomizer()->deficiency(mm,nn);
				
				for (int pp=0; pp<degree_deficiency-1; pp++) {  // -1 here because this is for the product rule
																// homogenize
					mul_mp(&( (*tempmat4)->entry[mm][nn]), &( (*curr_x_vars)->coord[0]),&( (*tempmat4)->entry[mm][nn]));
				}
				
				set_zero_mp(*temp);
				mpf_set_d((*temp)->r ,degree_deficiency);
				mul_mp(&( (*tempmat4)->entry[mm][nn]), *temp ,&( (*tempmat4)->entry[mm][nn]));
				
				
			}
		}
		
		
		for (int mm=0; mm<BED->randomizer()->num_base_funcs(); mm++) {
			for (int nn=1; nn<BED->num_natural_vars; nn++) {
				set_mp(&(*temp_jacobian_functions)->entry[mm][nn-1],&(*temp_jacobian_functions)->entry[mm][nn]);
			}
		}
		(*temp_jacobian_functions)->cols = BED->num_natural_vars-1;
		
		
	}
	
	
	
	// the hom variable
	int current_variable_index = 0;
	int local_offset = derivative_offset-1; // initialize.  we count up from here
	for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each base function
		int entry_counter = 0;
		
		for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
			for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
				local_offset++;//  MUST always increment.
				
				if (jj==0) {
					continue; // skip the ∂/∂h^2 entry
				}
				
				if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
					continue;
				}
				
				set_mp(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
				
				entry_counter++;
			}
		}
	} // re: making the interior matrix of second derivatives
	
	//need to homogenize properly, with the product rule
	
	
	// ok, now tempmat1 has in it all the (non-homogenized) derivatives of the derivatives of the system f.
	//  we are trying to get d/dx( [(R*Jf)^T | h^d•pi^T] v) = ( R*d/dx(Jf) )^T v_truncate.
	
	if (!BED->randomizer()->is_square()) { // homogenize and randomize
		mat_mul_mp(*tempmat3, *tempmat2, *tempmat1); // randomize
		
		
		// not the other part for the product rule with the homogenizing variable
		mat_mul_mp(*tempmat5, *tempmat4, *temp_jacobian_functions);  // n*h^{n-1} * f_{bla}
		
		
		
		mul_mat_vec_mp(*tempvec, *tempmat3, *curr_v_vars);
		mul_mat_vec_mp(*tempvec2, *tempmat5, *curr_v_vars);
		
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			add_mp(&(*tempvec2)->coord[mm], &(*tempvec)->coord[mm], &(*tempvec2)->coord[mm]);
		}
	}
	else
	{
		mul_mat_vec_mp(*tempvec2, *tempmat1, (*curr_v_vars));
	}
	
	
	
	vec_mulcomp_mp(*tempvec, *tempvec2, *one_minus_s);// now tempvec has the derivatives wrt $x$. (for all jac eqns)
	
	// combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
	for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
		mul_mp(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
		add_mp(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
	}
	
	
	// the remainder are all 0, cuz just have v in them.  this is probably unnecessary because the matrix was zero'd earlier.  (wasn't it?)
	for (int mm=0; mm<BED->num_projections; mm++) {
		set_zero_mp(&Jv->entry[mm+offset+BED->randomizer()->num_rand_funcs()][current_variable_index]);
	}
	//////  end the calculation for the hom var
	
	
	
	
	
	
	
	
	
	
	
	// THE X VARS (non-homogenizing vars)
	
	// start this loop at 1, because we computed the homogenizing variable differently
	for (current_variable_index=1; current_variable_index<BED->num_natural_vars+BED->num_synth_vars; current_variable_index++) {
		//		std::cout << "variable " << current_variable_index << std::endl;
		
		
		int local_offset = derivative_offset-1; // initialize this counter.  we count up from here.  must be -1 cause ended -1 when unpacking first derivatives earlier
		for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each base function
			
			int entry_counter = 0;
			
			
			for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
				for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
					local_offset++;//  MUST always increment.
					
					
					if (ii==0) { // skip the mixed derivatives involving the hom var
						continue;
					}
					
					if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
						continue;
					}
					
					// made it here, so this entry in deriv_vals is indeed one of the second derivatives we need
					set_mp(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
					
					entry_counter++;
					//
					
				}
			}
		} // re: making the interior matrix of second derivatives
		
		
		
		
		
		// ok, now tempmat1 has in it all the (non-homogenized) derivatives of the derivatives of the system f.
		//  we are trying to get d/dx( [(R*Jf)^T | h^d•pi^T] v) = ( R*d/dx(Jf) )^T v_truncate.
		
		if (!BED->randomizer()->is_square()) { // homogenize and randomize
			mat_mul_mp(*tempmat3, *tempmat2, *tempmat1); // randomize
		}
		else
		{
			mat_mp * supertemporary = tempmat1; // create a super temporary mat_mp*  for the swap
			tempmat1 = tempmat3;
			tempmat3 = supertemporary;
		}
		
		
		mul_mat_vec_mp(*tempvec, *tempmat3, (*curr_v_vars)); // recall v has been truncated to avoid unnecessary 0*number multiplication
		
		if (BED->randomizer()->is_square()) { // swap back.  this is because tempmat1 is already the correct size
			mat_mp * supertemporary = tempmat1;
			tempmat1 = tempmat3;
			tempmat3 = supertemporary;
		}
		
		vec_mulcomp_mp(*tempvec, *tempvec, *one_minus_s);// now tempvec has the derivatives wrt $x$. (for all jac eqns)
		
		// combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		for (int mm=0; mm<BED->randomizer()->num_rand_funcs(); mm++) {
			mul_mp(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
			add_mp(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
		}
		
		
		
		// the remainder are all 0, cuz just have v in them.  this is probably unnecessary because the matrix was zero'd earlier.  (wasn't it?)
		for (int mm=0; mm<BED->num_projections; mm++) {
			set_zero_mp(&Jv->entry[mm+offset+BED->randomizer()->num_rand_funcs()][current_variable_index]);
		}
	}//re: current_variable_index for diff
	
	// END DIFF wrt x
	
	
	//set the X PATCH values
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations;
	
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_mp(&funcVals->coord[ii+offset], &(*patchValues)->coord[ii]);
		
		// Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &(*Jv_Patch)->entry[ii][jj]);
	}
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 2151;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		throw std::logic_error("logic error thrown");
	}
	
	// V patch
	set_one_mp(*temp2);
	dot_product_mp(*temp, BED->v_patch, (*curr_v_vars));
	sub_mp(&funcVals->coord[BED->num_variables-1], *temp, *temp2);  // f = patch*v-1
	
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&Jv->entry[BED->num_variables-1][BED->num_natural_vars+ii], &BED->v_patch->coord[ii]);
	
	
	
	// finally, set parVals & parDer correctly
	
	change_size_point_mp(parVals, 1);  change_size_vec_mp(parDer, 1);
	parVals->size = parDer->size = 1;
	
	set_mp(&parVals->coord[0], pathVars); // s = t
	set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
	
	if (BED->verbose_level()==10 || BED->verbose_level()==-10) {
		print_comp_matlab(pathVars,"t");
		//	print_matrix_to_screen_matlab(jac_homogenizing_matrix,"jac_hom_1044");
		//	print_matrix_to_screen_matlab(BED->post_randomizer_matrix,"S");
		//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"R");
		//
		//
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
		print_point_to_screen_matlab(*curr_x_vars,"x1");
		print_point_to_screen_matlab(*curr_v_vars,"V1");
		print_point_to_screen_matlab(funcVals,"F");
		print_matrix_to_screen_matlab(Jv,"Jv");
		print_matrix_to_screen_matlab(Jp,"Jp");
		//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
		//			//these values are set in this function:  point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp
		//		print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");
		//	std::cout << "\n\n**************\n\n";
		if (BED->verbose_level()<0) {
			mypause();
		}
		
	}
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_x_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_mp(dehommed);
#endif
	
	//	printf("exiting eval\n");
	return 0;
}











int nullspacejac_left_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	
	BED->temp_vars.ensure_have_scalars(6);
	BED->temp_vars.ensure_have_vectors(13);
	BED->temp_vars.ensure_have_matrices(10);
	
	
	int offset;
	comp_d *one_minus_s = &BED->temp_vars.scalars[0], *gamma_s = &BED->temp_vars.scalars[1];
	
	set_one_d(*one_minus_s);
	sub_d(*one_minus_s, *one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_d(*gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	
	comp_d *running_prod = &BED->temp_vars.scalars[2];
	comp_d *temp = &BED->temp_vars.scalars[3], *temp2 = &BED->temp_vars.scalars[4], *temp3 = &BED->temp_vars.scalars[5];
	
	
	
	
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	vec_d *curr_x_vars = &BED->temp_vars.vectors[0];
	increase_size_vec_d(*curr_x_vars, BED->num_natural_vars);
	(*curr_x_vars)->size = BED->num_natural_vars;
	for (int ii=0; ii<BED->num_natural_vars; ii++)
		set_d(&(*curr_x_vars)->coord[ii], &current_variable_values->coord[ii]);
	
	vec_d *curr_v_vars = &BED->temp_vars.vectors[1];
	increase_size_vec_d(*curr_v_vars, BED->num_v_vars);
	(*curr_v_vars)->size = BED->num_v_vars;
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_d(&(*curr_v_vars)->coord[ii], &current_variable_values->coord[ii+BED->num_natural_vars]);
	
	
	
	
	vec_d *patchValues = &BED->temp_vars.vectors[2];
	vec_d *temp_function_values = &BED->temp_vars.vectors[3];
	
	
	vec_d *AtimesF = &BED->temp_vars.vectors[4];
	vec_d *linprod_x = &BED->temp_vars.vectors[5];  change_size_vec_d(*linprod_x, BED->num_jac_equations);
	(*linprod_x)->size = BED->num_jac_equations;
	
	vec_d *linprod_times_gamma_s = &BED->temp_vars.vectors[6]; change_size_vec_d(*linprod_times_gamma_s,BED->num_jac_equations);
	(*linprod_times_gamma_s)->size = BED->num_jac_equations;
	vec_d *tempvec = &BED->temp_vars.vectors[7];
	vec_d *tempvec2 = &BED->temp_vars.vectors[8];
	
	vec_d *target_function_values = &BED->temp_vars.vectors[9];
	vec_d *target_function_values_times_oneminus_s = &BED->temp_vars.vectors[10];
	
	vec_d *start_function_values = &BED->temp_vars.vectors[11];
	increase_size_vec_d(*start_function_values,BED->num_jac_equations); (*start_function_values)->size = BED->num_jac_equations;
	
	
	
	vec_d *deriv_vals = &BED->temp_vars.vectors[12];
	
	
	
	
	
	
	mat_d *Jv_Patch = &BED->temp_vars.matrices[0];
	mat_d *tempmat = &BED->temp_vars.matrices[1]; change_size_mat_d(*tempmat,BED->num_variables-1,BED->num_variables-1);
	(*tempmat)->rows = (*tempmat)->cols = BED->num_variables-1; // change the size indicators
	
	mat_d *lin_func_vals = &BED->temp_vars.matrices[2]; change_size_mat_d(*lin_func_vals,BED->num_jac_equations, BED->max_degree);
	(*lin_func_vals)->rows = BED->num_jac_equations; (*lin_func_vals)->cols = BED->max_degree;
	
	
	
	mat_d *AtimesJ = &BED->temp_vars.matrices[3];
	
	
	mat_d *linprod_derivative_wrt_x = &BED->temp_vars.matrices[4];
	increase_size_mat_d(*linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_natural_vars);
	(*linprod_derivative_wrt_x)->rows = BED->num_jac_equations; (*linprod_derivative_wrt_x)->cols = BED->num_natural_vars;
	
	
	
	
	mat_d* Jf_pi_homogenized = &BED->temp_vars.matrices[5];  increase_size_mat_d(*Jf_pi_homogenized, BED->num_natural_vars-1, BED->num_v_vars); // set up temp matrix
	(*Jf_pi_homogenized)->rows = BED->num_natural_vars-1; (*Jf_pi_homogenized)->cols = BED->num_v_vars;
	
	
	// create temp matrices
	mat_d *tempmat1 = &BED->temp_vars.matrices[6];
	mat_d *tempmat2 = &BED->temp_vars.matrices[7];
	mat_d *tempmat3 = &BED->temp_vars.matrices[8];
	
	increase_size_mat_d(*tempmat1,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat1)->rows = BED->num_jac_equations; (*tempmat1)->cols = BED->num_v_vars;
	
	increase_size_mat_d(*tempmat2,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat2)->rows = BED->num_jac_equations; (*tempmat2)->cols = BED->num_v_vars;
	
	
	mat_d *temp_jacobian_functions = &BED->temp_vars.matrices[9];
	
	
	
	
	
	
	// the main evaluations for $x$
	evalDeriv_d(*temp_function_values, *deriv_vals, *tempvec, *tempvec2, *curr_x_vars, BED->SLP_derivative);
	//tempvec, and tempvec2, are throwaway variables used to capture output from this which we don't care about -- the linears added to the system.  perhaps we could optimize this?
	
	
	// evaluate the patch
	patch_eval_d(    *patchValues, parVals, parDer, *Jv_Patch, Jp, *curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	
	
	
	
	
	
	
	
	//resize output variables to correct size
	increase_size_vec_d(funcVals,BED->num_variables);
	increase_size_mat_d(Jv, BED->num_variables, BED->num_variables);
	increase_size_mat_d(Jp, BED->num_variables, 1);
	funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
	Jv->cols = BED->num_variables;  //  <-- this must be square
	Jp->cols = 1; // this is a column vector, or nx1 matrix, for this evaluator
	
	
	
	//////
	// initialize stuff to all 0's
	///////
	
	for (int ii=0; ii<Jv->rows; ii++){
		for (int jj=0; jj<Jv->cols; jj++){
			set_zero_d(&Jv->entry[ii][jj]);
			// initialize entire matrix to 0 // this is kinda wasteful, but safe.  wasteful because most entries are set in this function, but safe because we only have to set entries which we know aren't 0, and rest are 0 by default.
		}
	}
	
	for (int ii = 0; ii<BED->num_variables; ii++)
		set_zero_d(&Jp->entry[ii][0]);  // initialize entire matrix to 0.  same comment regarding wastefulness.
	
	
	
	
	
	///////////
	// orig eqns
	/////////////
	
	
	// resize temp_jacobians to be Nxn, to hold the first derivatives.
	change_size_mat_d(*temp_jacobian_functions, (*temp_function_values)->size, (BED->num_natural_vars+BED->num_synth_vars));
	(*temp_jacobian_functions)->rows = (*temp_function_values)->size;
	(*temp_jacobian_functions)->cols = (BED->num_natural_vars+BED->num_synth_vars);
	
	
	// unpack the first derivatives
	int derivative_offset = 0;
	for (int ii = 0; ii<(*temp_function_values)->size; ii++) {
		for (int jj = 0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			set_d(&(*temp_jacobian_functions)->entry[ii][jj],&(*deriv_vals)->coord[derivative_offset]);
			derivative_offset++; // the end of this loop will leave this offset for use later.
		}
	} // the remainer of deriv_vals contains the second derivatives, which will be used later in this function
	
	
	
	
	
	// randomize
	BED->randomizer()->randomize(*AtimesF,*AtimesJ,*temp_function_values,*temp_jacobian_functions, &(*curr_x_vars)->coord[0]);
	
	for (int ii=0; ii<(*AtimesF)->size; ii++)  // for each function, after (real) randomization
		set_d(&funcVals->coord[ii], &(*AtimesF)->coord[ii]);
	
	
	
	
	
	//////////////
	//  orig function jacobian values
	////////////////
	
	
	
	
	// set the jacobian equations for orig into Jv
	
	// copy the jacobian into the return value for the evaluator
	for (int ii=0; ii< (*AtimesJ)->rows; ii++) // for every function
		for (int jj=0; jj< (*AtimesJ)->cols; jj++)
			set_d(&Jv->entry[ii][jj],&(*AtimesJ)->entry[ii][jj]);
	
	
	
	
	
	
	
	//////////
	// the additional linears.  there are $r-\ell$ of them.
	///////////
	
	
	
	offset = BED->randomizer()->num_rand_funcs();
	for (int ii=0; ii< BED->num_additional_linears; ii++) {
		
		dot_product_d(*temp, BED->additional_linears_terminal[ii], *curr_x_vars); //temp = terminal(x)
		neg_d(&Jp->entry[offset+ii][0], *temp);  // Jp = -terminal(x)
		mul_d(*temp3, *temp, *one_minus_s);        // temp3 = (1-s)*terminal(x)
		
		dot_product_d(*temp,  BED->additional_linears_starting[ii], *curr_x_vars); // temp = starting(x)
		mul_d(*temp2, BED->gamma, *temp);         // temp2 = gamma*starting(x)
		
		add_d(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], *temp2);   // Jp = -terminal + gamma*start
		
		mul_d(*temp2, *gamma_s, *temp);            // temp2 = gamma*s*starting(x)
		
		add_d(&funcVals->coord[offset+ii],*temp2, *temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (int jj=0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			mul_d(*temp, *gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]);
			mul_d(*temp2,*one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_d(&Jv->entry[ii+offset][jj], *temp, *temp2);
		}
	}
	
	
	
	
	
	
	/////////////////////////////
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	////////////////////////////////
	
	
	
	// HOMOGENIZE the jacobian matrix for the supplemental equations.
	
	for (int ii=0; ii<(*temp_jacobian_functions)->rows; ++ii) {
		int degree_deficiency = BED->max_degree+1-BED->randomizer()->base_degree(ii);
		
		if (degree_deficiency==0) { // jacobian value already set with correct homogenization
			continue;
		}
		else{
			set_d(*temp,&(*curr_x_vars)->coord[0]);
			for (int jj=1; jj<degree_deficiency; ++jj) {
				mul_d(*temp,*temp,&(*curr_x_vars)->coord[0])
			}//TODO: optimize this away.  we use this all the time, and this calculation is repeated many times.
			
			for (int jj=0; jj<(*temp_jacobian_functions)->cols; ++jj) {
				mul_d(&(*temp_jacobian_functions)->entry[ii][jj],&(*temp_jacobian_functions)->entry[ii][jj],*temp);
			}
		}
	}
	
	
	
	if (!BED->randomizer()->is_square()) {
		mat_mul_d(*tempmat1,*(BED->randomizer()->get_mat_d()),*temp_jacobian_functions);
	}
	else
	{
		mat_d *supertemporary = tempmat1;
		tempmat1 = temp_jacobian_functions;
		temp_jacobian_functions = supertemporary;
	}
	
	change_size_mat_d(BED->jac_with_proj,BED->num_natural_vars-1,BED->randomizer()->num_rand_funcs()+BED->num_projections);
	BED->jac_with_proj->rows = BED->num_natural_vars-1;
	BED->jac_with_proj->cols = BED->randomizer()->num_rand_funcs()+BED->num_projections;
	
	
	
	// transpose into beginning columns jac_with_proj
	for (int ii=0; ii< BED->randomizer()->num_rand_funcs(); ii++)// for every function
		for (int jj=1; jj<BED->num_natural_vars; jj++) // for only the natural variables, omitting the hom var and synth_vars
			set_d(&BED->jac_with_proj->entry[jj-1][ii], &(*tempmat1)->entry[ii][jj]);
	
	
	
	
	// concatenate in the values of the projections, homogenized, to get [Jv | pi]
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=1; jj<BED->num_natural_vars; jj++) {
			set_d(*temp, &BED->target_projection[ii]->coord[jj]); // seed the homogenization loop with the projection coefficient
			for (int mm=0; mm<BED->max_degree; mm++) { // max_degree is highest degree of a derivative.
				mul_d(*temp,*temp,&(*curr_x_vars)->coord[0]); // homogenize
															   //TODO: this could be optimized
			}
			set_d(&BED->jac_with_proj->entry[jj-1][ii+BED->randomizer()->num_rand_funcs()],*temp); // set into the transpose, to the right of the jacobian matrix
		}
	} // this verified correct for paraboloid and torus base systems using matlab
	
	
	
	
	//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jac_hom_w_proj");
	
	mul_mat_vec_d(*target_function_values, BED->jac_with_proj, (*curr_v_vars));
	
	//	print_point_to_screen_matlab(target_function_values,"target_f");
	vec_mulcomp_d(*target_function_values_times_oneminus_s, *target_function_values, *one_minus_s);
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	// the product of the linears
	for (int jj=0; jj<BED->num_jac_equations; jj++) {
		
		//perform the $x$ evaluation
		set_one_d(&(*linprod_x)->coord[jj]); // initialize to 1 for multiplication
		for (int ii=0; ii<BED->max_degree; ++ii) {
			dot_product_d(&(*lin_func_vals)->entry[jj][ii],BED->starting_linears[jj][ii],*curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_d(&(*linprod_x)->coord[jj], &(*linprod_x)->coord[jj], &(*lin_func_vals)->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_d(*temp, BED->v_linears[jj], (*curr_v_vars));
		
		//now set the combined $x,v$ value
		mul_d(&(*start_function_values)->coord[jj], &(*linprod_x)->coord[jj], *temp); // start_function_values = linprod_x(x)*linear(v)
		mul_d(&(*linprod_times_gamma_s)->coord[jj], *gamma_s, &(*start_function_values)->coord[jj]); // sets the value linprod_x*gamma*s
		
		
		neg_d( &Jp->entry[jj + offset][0], &(*target_function_values)->coord[jj]);  // Jp = -target
		mul_d( *temp, BED->gamma, &(*start_function_values)->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_d( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], *temp); // Jp = Jp + gamma*start = -target + gamma*start
	}
	
	
	// mix the start and target into final function value to return.
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		add_d(&funcVals->coord[ii+offset], &(*target_function_values_times_oneminus_s)->coord[ii], &(*linprod_times_gamma_s)->coord[ii]);
	}
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE start and target WRT V
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_d(*temp, &BED->v_linears[ii]->coord[jj], &(*linprod_x)->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_d(*temp2, *temp, *gamma_s);                                         // temp2 = gamma*s*M_ij * linprod_x(x)
			
			mul_d(*temp, *one_minus_s, &BED->jac_with_proj->entry[ii][jj]);        // temp = (1-s)*([Jf^T pi^T])_ij
			add_d(&Jv->entry[ii+offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2);       // Jv = temp + temp2
		}
	}
	
	
	//////////////////
	// now the x derivatives corresponding to the linprod start system
	//////////////////////
	
	
	
	// an implementation of the product rule
	for (int mm=0; mm< BED->num_jac_equations; mm++) {
		for (int kk=0; kk<(BED->num_natural_vars+BED->num_synth_vars); kk++) { // for each variable
			set_zero_d(&(*linprod_derivative_wrt_x)->entry[mm][kk]); // initialize to 0 for the sum
			
			for (int ii=0; ii<BED->max_degree; ++ii) { //  for each linear
				
				set_d(*running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (int jj=0; jj<BED->max_degree; jj++) {
					if (jj!=ii) {
						mul_d(*running_prod,*running_prod,&(*lin_func_vals)->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_d(&(*linprod_derivative_wrt_x)->entry[mm][kk],&(*linprod_derivative_wrt_x)->entry[mm][kk],*running_prod);
				
			}// re:ii
			
			dot_product_d(*temp, BED->v_linears[mm], (*curr_v_vars));  // these two lines multiply by  (v_linear •	v)
			mul_d(&(*linprod_derivative_wrt_x)->entry[mm][kk], &(*linprod_derivative_wrt_x)->entry[mm][kk], *temp);
		} // re: kk
	} // re: mm
	
	
	
	///////////////
	// DIFFERENTIATE THE derivative of the target jacobian system wrt $h,x_i$.
	/////////////////
	
	//  compute ∂/∂h
	
	// tempmat1 will hold the second derivatives, omitting the homvar column (hence the minus 1), and omitting all synth variables, too
	increase_size_mat_d(*tempmat1,(*temp_function_values)->size,BED->num_natural_vars-1);
	(*tempmat1)->rows = (*temp_function_values)->size; (*tempmat1)->cols = BED->num_natural_vars-1;
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	
	
	// this one is different from the x derivatives because h is the homogenizing variable for the problem,
	// and there are terms which contain only constants (pi) and the homvar.
	
	int current_variable_index = 0; // later we will loop from 1:num_natural_vars
	int local_offset = derivative_offset-1; // local_offset indexes into the vector containing second derivatives.
											//  the -1 is because we increment this at the beginning of the nested for loop
	
	for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each function
		
		int degree_deficiency = BED->max_degree+1 - BED->randomizer()->base_degree(base_func_index);// +1 is because maxdegree is deg of derivatives,
																							 //and base_degrees is deg of natural function before differentiation.
		int entry_counter = 0; // indexes into the column we will set.  easier this way.
		for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
			for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) { // it's essentially an upper triangular matrix passed through reshape(-,[],1)
				local_offset++;//  MUST always increment.   cannot 'continue' before this.
				
				if (jj==0) {// for this particular nullspace implementation calculation, we omit the hom var.  one would not do this if they needed ALL second derivatives.
					continue; //
				}
				
				if (ii!=current_variable_index && jj!=current_variable_index) {
					continue;
				}
				
				// g = h^n * function = h^n df_{curr_func}/dx_{curr_var}
				// dg/dh = n*h^(n-1)function + h^n dfunction/dh
				if (degree_deficiency==0) { // no need to further homogenize, n=0 in above formula
					set_d(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
				}
				else{ // must homogenize, n>0 in above formula
					
					//					std::cout << "homogenizing deg def " << degree_deficiency << std::endl;
					//
					//					print_comp_matlab(&deriv_vals->coord[local_offset],"second_der");
					// first, homogenize against the derivative number of times according to deficiency.
					//   h^n (∂^2)f/(∂x∂x)
					set_d(*temp,&(*deriv_vals)->coord[local_offset]); // seed, grab the first derivative
					for (int zz = 0; zz<degree_deficiency; zz++) {
						mul_d(*temp,*temp,&(*curr_x_vars)->coord[0]); // loop to multiply
					}
					// temp = h^n*∂f/∂x
					
					
					// next, homogenize against the non-derivative part (function)
					//   n h^(n-1) f
					set_d(*temp3, &(*deriv_vals)->coord[base_func_index*(BED->num_natural_vars+BED->num_synth_vars)+jj]); // seed
					for (int zz=0; zz<degree_deficiency-1; zz++) {  // i emphasize the -1 in upper limit for this loop (degree_deficiency-1)
						mul_d(*temp3, *temp3, &(*curr_x_vars)->coord[0]); // loop multiply
					}
					
					set_zero_d(*temp2);
					(*temp2)->r = degree_deficiency;
					//					temp2->r = degree_deficiency; temp2->i = 0; // temp2 = n
					mul_d(*temp3,*temp2,*temp3);  // multiply by the coefficient according to the degree.
												   // temp3 = n*h^{n-1} f
					
					// put 'em together
					add_d(&(*tempmat1)->entry[base_func_index][entry_counter], *temp3, *temp);
					
				}
				entry_counter++;
			}
		}
	}
	// ok, now tempmat1 has in it all the derivatives of the derivatives of the system f.
	//  we are trying to get d/dx( [(R*Jf)^T | pi^T] v).
	
	if (!BED->randomizer()->is_square()) {
		mat_mul_d(*tempmat2, *(BED->randomizer()->get_mat_d()), *tempmat1); // randomize
		
		nonconj_transpose(*tempmat3,*tempmat2);
		
		
	}
	else
	{
		nonconj_transpose(*tempmat3,*tempmat1);
		
	}
	
	
	
	
	
	
	//	print_matrix_to_screen_matlab(tempmat3,"hom_deriv_mat");
	increase_size_mat_d(*tempmat3,(*tempmat3)->rows, BED->num_v_vars); // nondestructive resize
	(*tempmat3)->cols = BED->num_v_vars;
	
	
	// homogenize the projection terms for the concatenated jacobian
	set_one_d(*temp2);
	(*temp2)->r = BED->max_degree;
	for (int mm=0; mm<BED->max_degree-1; mm++)  // remember, max_degree is the highest degree occurring in the DERIVATIVES
		mul_d(*temp2,*temp2,&(*curr_x_vars)->coord[0]); // we will use this in a line or two to homogenize
														 // temp2 = n*h^{n-1}
	
	
	for (int mm=0; mm<BED->num_projections; mm++) {
		for (int nn=1; nn<BED->num_natural_vars; nn++) {
			mul_d(&(*tempmat3)->entry[nn-1][BED->randomizer()->num_rand_funcs()+mm], &BED->target_projection[mm]->coord[nn],*temp2);
			// entry = pi_mm[nn]*n*h^{n-1}
		}
	}
	
	
	
	
	mul_mat_vec_d(*tempvec, *tempmat3, (*curr_v_vars)); // multiply against the $v$ variables
														 //	print_point_to_screen_matlab(tempvec,"diff_h_times_v");
	vec_mulcomp_d(*tempvec, *tempvec, *one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable 0. (for all jac eqns), (with current time taken into account)
	
	// now, combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
	
	for (int mm=0; mm<BED->num_jac_equations; mm++) {
		mul_d(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
		add_d(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
	}
	
	// end diff WRT hom_var
	
	
	
	// THE X VARS
	
	change_size_mat_d(*tempmat1,BED->randomizer()->num_base_funcs(),BED->num_natural_vars+BED->num_synth_vars);
	
	
	(*curr_v_vars)->size -= BED->num_projections; // truncate.  this is for optimization, not theoretical reasons
												  //
	for (current_variable_index=1; current_variable_index<BED->num_natural_vars+BED->num_synth_vars; current_variable_index++) {
		
		
		local_offset = derivative_offset-1; // initialize.  we count up from here
		for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each function
			
			int degree_deficiency = BED->max_degree+1 - BED->randomizer()->base_degree(base_func_index);
			int entry_counter = 0;
			for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
				for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
					local_offset++;//  MUST always increment.
					
					if (ii==0) {
						continue;
					}
					
					if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
						continue;
					}
					
					
					set_d(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
					
					
					// g = h^n * f
					// dg/dx = h^n df/dx
					if (degree_deficiency>0){
						for (int zz = 0; zz<degree_deficiency; zz++) {
							mul_d(&(*tempmat1)->entry[base_func_index][entry_counter],&(*tempmat1)->entry[base_func_index][entry_counter],&(*curr_x_vars)->coord[0]); // loop to multiply
						}
					}
					
					entry_counter++;
				}
			}
		}
		
		
		
		
		
		// ok, now tempmat1 has in it all the (homogenized) derivatives of the derivatives of the system f.
		//  we are trying to get d/dx( [(R*Jf)^T | h^d•pi^T] v) = ( R*d/dx(Jf) )^T v_truncate.
		
		if (!BED->randomizer()->is_square()) {
			mat_mul_d(*tempmat2, *(BED->randomizer()->get_mat_d()), *tempmat1); // randomize
			nonconj_transpose(*tempmat3,*tempmat2);
		}
		else
		{
			nonconj_transpose(*tempmat3,*tempmat1);
		}
		
		
		
		// no need to copy in the homogenized pi coefficients, because the derivatives wrt x are all 0.
		
		//		print_matrix_to_screen_matlab(tempmat3,"xxx");
		mul_mat_vec_d(*tempvec, *tempmat3, (*curr_v_vars)); // recall v has been truncated to avoid unnecessary 0*number multiplication
															 //		print_point_to_screen_matlab(tempvec,"diff_x_times_v");
		vec_mulcomp_d(*tempvec, *tempvec, *one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable 0. (for all jac eqns)
		
		// combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		for (int mm=0; mm<BED->num_jac_equations; mm++) {
			mul_d(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
			add_d(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
		}
		
	}//re: current_variable_index for numerical diff
	
	// END DIFF wrt x
	
	
	
	
	
	
	(*curr_v_vars)->size += BED->num_projections; // un-truncate.  the values of the trailing variables remain unaffected by previous truncation.
	
	
	
	//set the X PATCH values
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations;
	
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_d(&funcVals->coord[ii+offset], &(*patchValues)->coord[ii]);
		
		// Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &(*Jv_Patch)->entry[ii][jj]);
	}
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 2151;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		throw std::logic_error("logic error thrown");
	}
	
	// V patch
	set_one_d(*temp2);
	dot_product_d(*temp, BED->v_patch, (*curr_v_vars));
	sub_d(&funcVals->coord[BED->num_variables-1], *temp, *temp2);  // f = patch*v-1
	
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_d(&Jv->entry[BED->num_variables-1][BED->num_natural_vars+ii], &BED->v_patch->coord[ii]);
	
	
	
	// finally, set parVals & parDer correctly
	
	change_size_point_d(parVals, 1);  change_size_vec_d(parDer, 1);
	parVals->size = parDer->size = 1;
	
	set_d(&parVals->coord[0], pathVars); // s = t
	set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	if (BED->verbose_level()==10 || BED->verbose_level()==-10) {
		printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		//	print_matrix_to_screen_matlab(jac_homogenizing_matrix,"jac_hom_1044");
		//	print_matrix_to_screen_matlab(BED->post_randomizer_matrix,"S");
		//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"R");
		//
		//
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
		print_point_to_screen_matlab(*curr_x_vars,"x1");
		print_point_to_screen_matlab(*curr_v_vars,"V1");
		print_point_to_screen_matlab(funcVals,"F");
		print_matrix_to_screen_matlab(Jv,"Jv");
		print_matrix_to_screen_matlab(Jp,"Jp");
		//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
		//			//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
//		print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");
		//	std::cout << "\n\n**************\n\n";
		if (BED->verbose_level()<0) {
			mypause();
		}
		
	}
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_x_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	//	printf("exiting eval\n");
	return 0;
}




int nullspacejac_left_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	BED->temp_vars.ensure_have_scalars(6);
	BED->temp_vars.ensure_have_vectors(13);
	BED->temp_vars.ensure_have_matrices(10);
	
	
	int offset;
	comp_mp *one_minus_s = &BED->temp_vars.scalars[0], *gamma_s = &BED->temp_vars.scalars[1];
	
	set_one_mp(*one_minus_s);
	sub_mp(*one_minus_s, *one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_mp(*gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	
	comp_mp *running_prod = &BED->temp_vars.scalars[2];
	comp_mp *temp = &BED->temp_vars.scalars[3], *temp2 = &BED->temp_vars.scalars[4], *temp3 = &BED->temp_vars.scalars[5];
	
	
	
	
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	vec_mp *curr_x_vars = &BED->temp_vars.vectors[0];
	increase_size_vec_mp(*curr_x_vars, BED->num_natural_vars);
	(*curr_x_vars)->size = BED->num_natural_vars;
	for (int ii=0; ii<BED->num_natural_vars; ii++)
		set_mp(&(*curr_x_vars)->coord[ii], &current_variable_values->coord[ii]);
	
	vec_mp *curr_v_vars = &BED->temp_vars.vectors[1];
	increase_size_vec_mp(*curr_v_vars, BED->num_v_vars);
	(*curr_v_vars)->size = BED->num_v_vars;
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&(*curr_v_vars)->coord[ii], &current_variable_values->coord[ii+BED->num_natural_vars]);
	
	
	
	
	vec_mp *patchValues = &BED->temp_vars.vectors[2];
	vec_mp *temp_function_values = &BED->temp_vars.vectors[3];
	
	
	vec_mp *AtimesF = &BED->temp_vars.vectors[4];
	vec_mp *linprod_x = &BED->temp_vars.vectors[5];  change_size_vec_mp(*linprod_x, BED->num_jac_equations);
	(*linprod_x)->size = BED->num_jac_equations;
	
	vec_mp *linprod_times_gamma_s = &BED->temp_vars.vectors[6]; change_size_vec_mp(*linprod_times_gamma_s,BED->num_jac_equations);
	(*linprod_times_gamma_s)->size = BED->num_jac_equations;
	vec_mp *tempvec = &BED->temp_vars.vectors[7];
	vec_mp *tempvec2 = &BED->temp_vars.vectors[8];
	
	vec_mp *target_function_values = &BED->temp_vars.vectors[9];
	vec_mp *target_function_values_times_oneminus_s = &BED->temp_vars.vectors[10];
	
	vec_mp *start_function_values = &BED->temp_vars.vectors[11];
	increase_size_vec_mp(*start_function_values,BED->num_jac_equations); (*start_function_values)->size = BED->num_jac_equations;
	
	
	
	vec_mp *deriv_vals = &BED->temp_vars.vectors[12];
	
	
	
	
	
	
	mat_mp *Jv_Patch = &BED->temp_vars.matrices[0];
	mat_mp *tempmat = &BED->temp_vars.matrices[1]; change_size_mat_mp(*tempmat,BED->num_variables-1,BED->num_variables-1);
	(*tempmat)->rows = (*tempmat)->cols = BED->num_variables-1; // change the size indicators
	
	mat_mp *lin_func_vals = &BED->temp_vars.matrices[2]; change_size_mat_mp(*lin_func_vals,BED->num_jac_equations, BED->max_degree);
	(*lin_func_vals)->rows = BED->num_jac_equations; (*lin_func_vals)->cols = BED->max_degree;
	
	
	
	mat_mp *AtimesJ = &BED->temp_vars.matrices[3];
	
	
	mat_mp *linprod_derivative_wrt_x = &BED->temp_vars.matrices[4];
	increase_size_mat_mp(*linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_natural_vars);
	(*linprod_derivative_wrt_x)->rows = BED->num_jac_equations; (*linprod_derivative_wrt_x)->cols = BED->num_natural_vars;
	
	
	
	
	mat_mp* Jf_pi_homogenized = &BED->temp_vars.matrices[5];  increase_size_mat_mp(*Jf_pi_homogenized, BED->num_natural_vars-1, BED->num_v_vars); // set up temp matrix
	(*Jf_pi_homogenized)->rows = BED->num_natural_vars-1; (*Jf_pi_homogenized)->cols = BED->num_v_vars;
	
	
	// create temp matrices
	mat_mp *tempmat1 = &BED->temp_vars.matrices[6];
	mat_mp *tempmat2 = &BED->temp_vars.matrices[7];
	mat_mp *tempmat3 = &BED->temp_vars.matrices[8];
	
	increase_size_mat_mp(*tempmat1,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat1)->rows = BED->num_jac_equations; (*tempmat1)->cols = BED->num_v_vars;
	
	increase_size_mat_mp(*tempmat2,BED->num_jac_equations,BED->num_v_vars);
	(*tempmat2)->rows = BED->num_jac_equations; (*tempmat2)->cols = BED->num_v_vars;
	
	
	mat_mp *temp_jacobian_functions = &BED->temp_vars.matrices[9];
	

	
	
	
	
	// the main evaluations for $x$
	evalDeriv_mp(*temp_function_values, *deriv_vals, *tempvec, *tempvec2, *curr_x_vars, BED->SLP_derivative);
	//tempvec, and tempvec2, are throwaway variables used to capture output from this which we don't care about -- the linears added to the system.  perhaps we could optimize this?
	
	
	// evaluate the patch
	patch_eval_mp(    *patchValues, parVals, parDer, *Jv_Patch, Jp, *curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	
	
	
	
	
	
	
	
	//resize output variables to correct size
	increase_size_vec_mp(funcVals,BED->num_variables);
	increase_size_mat_mp(Jv, BED->num_variables, BED->num_variables);
	increase_size_mat_mp(Jp, BED->num_variables, 1);
	funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
	Jv->cols = BED->num_variables;  //  <-- this must be square
	Jp->cols = 1; // this is a column vector, or nx1 matrix, for this evaluator
	
	
	
	//////
	// initialize stuff to all 0's
	///////
	
	for (int ii=0; ii<Jv->rows; ii++){
		for (int jj=0; jj<Jv->cols; jj++){
			set_zero_mp(&Jv->entry[ii][jj]);
			// initialize entire matrix to 0 // this is kinda wasteful, but safe.  wasteful because most entries are set in this function, but safe because we only have to set entries which we know aren't 0, and rest are 0 by default.
		}
	}
	
	for (int ii = 0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);  // initialize entire matrix to 0.  same comment regarding wastefulness.
	
	
	
	
	
	///////////
	// orig eqns
	/////////////
	
	
	// resize temp_jacobians to be Nxn, to hold the first derivatives.
	change_size_mat_mp(*temp_jacobian_functions, (*temp_function_values)->size, (BED->num_natural_vars+BED->num_synth_vars));
	(*temp_jacobian_functions)->rows = (*temp_function_values)->size;
	(*temp_jacobian_functions)->cols = (BED->num_natural_vars+BED->num_synth_vars);
	
	
	// unpack the first derivatives
	int derivative_offset = 0;
	for (int ii = 0; ii<(*temp_function_values)->size; ii++) {
		for (int jj = 0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			set_mp(&(*temp_jacobian_functions)->entry[ii][jj],&(*deriv_vals)->coord[derivative_offset]);
			derivative_offset++; // the end of this loop will leave this offset for use later.
		}
	} // the remainer of deriv_vals contains the second derivatives, which will be used later in this function
	
	
	
	
	
	// randomize
	BED->randomizer()->randomize(*AtimesF,*AtimesJ,*temp_function_values,*temp_jacobian_functions, &(*curr_x_vars)->coord[0]);
	
	for (int ii=0; ii<(*AtimesF)->size; ii++)  // for each function, after (real) randomization
		set_mp(&funcVals->coord[ii], &(*AtimesF)->coord[ii]);
	
	
	
	
	
	//////////////
	//  orig function jacobian values
	////////////////
	
	
	
	
	// set the jacobian equations for orig into Jv
	
	// copy the jacobian into the return value for the evaluator
	for (int ii=0; ii< (*AtimesJ)->rows; ii++) // for every function
		for (int jj=0; jj< (*AtimesJ)->cols; jj++)
			set_mp(&Jv->entry[ii][jj],&(*AtimesJ)->entry[ii][jj]);
	
	
	
	
	
	
	
	//////////
	// the additional linears.  there are $r-\ell$ of them.
	///////////
	
	
	
	offset = BED->randomizer()->num_rand_funcs();
	for (int ii=0; ii< BED->num_additional_linears; ii++) {
		
		dot_product_mp(*temp, BED->additional_linears_terminal[ii], *curr_x_vars); //temp = terminal(x)
		neg_mp(&Jp->entry[offset+ii][0], *temp);  // Jp = -terminal(x)
		mul_mp(*temp3, *temp, *one_minus_s);        // temp3 = (1-s)*terminal(x)
		
		dot_product_mp(*temp,  BED->additional_linears_starting[ii], *curr_x_vars); // temp = starting(x)
		mul_mp(*temp2, BED->gamma, *temp);         // temp2 = gamma*starting(x)
		
		add_mp(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], *temp2);   // Jp = -terminal + gamma*start
		
		mul_mp(*temp2, *gamma_s, *temp);            // temp2 = gamma*s*starting(x)
		
		add_mp(&funcVals->coord[offset+ii],*temp2, *temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (int jj=0; jj<(BED->num_natural_vars+BED->num_synth_vars); jj++) {
			mul_mp(*temp, *gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]);
			mul_mp(*temp2,*one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_mp(&Jv->entry[ii+offset][jj], *temp, *temp2);
		}
	}
	
	
	
	
	
	
	/////////////////////////////
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	////////////////////////////////
	
	
	
	// HOMOGENIZE the jacobian matrix for the supplemental equations.
	
	for (int ii=0; ii<(*temp_jacobian_functions)->rows; ++ii) {
		int degree_deficiency = BED->max_degree+1-BED->randomizer()->base_degree(ii);
		
		if (degree_deficiency==0) { // jacobian value already set with correct homogenization
			continue;
		}
		else{
			set_mp(*temp,&(*curr_x_vars)->coord[0]);
			for (int jj=1; jj<degree_deficiency; ++jj) {
				mul_mp(*temp,*temp,&(*curr_x_vars)->coord[0])
			}//TODO: optimize this away.  we use this all the time, and this calculation is repeated many times.
			
			for (int jj=0; jj<(*temp_jacobian_functions)->cols; ++jj) {
				mul_mp(&(*temp_jacobian_functions)->entry[ii][jj],&(*temp_jacobian_functions)->entry[ii][jj],*temp);
			}
		}
	}
	
	
	
	if (!BED->randomizer()->is_square()) {
		mat_mul_mp(*tempmat1,*(BED->randomizer()->get_mat_mp()),*temp_jacobian_functions);
	}
	else
	{
		mat_mp * supertemporary = temp_jacobian_functions;
		temp_jacobian_functions = tempmat1;
		tempmat1 = supertemporary;
//			mat_cp_mp(*tempmat1,*temp_jacobian_functions);
	}
	
	increase_size_mat_mp(BED->jac_with_proj,BED->num_natural_vars-1,BED->randomizer()->num_rand_funcs()+BED->num_projections);
	BED->jac_with_proj->rows = BED->num_natural_vars-1;
	BED->jac_with_proj->cols = BED->randomizer()->num_rand_funcs()+BED->num_projections;
	
	// transpose into beginning columns jac_with_proj
	for (int ii=0; ii< (*tempmat1)->rows; ii++)// for every function
		for (int jj=1; jj<BED->num_natural_vars; jj++) // for only the natural variables, omitting the hom var and synth_vars
			set_mp(&BED->jac_with_proj->entry[jj-1][ii], &(*tempmat1)->entry[ii][jj]);
	
	
	
	
	// concatenate in the values of the projections, homogenized, to get [Jv | pi]
	for (int ii=0; ii<BED->num_projections; ii++) {
		for (int jj=1; jj<BED->num_natural_vars; jj++) {
			set_mp(*temp, &BED->target_projection[ii]->coord[jj]); // seed the homogenization loop with the projection coefficient
			for (int mm=0; mm<BED->max_degree; mm++) { // max_degree is highest degree of a derivative.
				mul_mp(*temp,*temp,&(*curr_x_vars)->coord[0]); // homogenize
														 // this could be optimized
			}
			set_mp(&BED->jac_with_proj->entry[jj-1][ii+BED->randomizer()->num_rand_funcs()],*temp); // set into the transpose, to the right of the jacobian matrix
		}
	} // this verified correct for paraboloid and torus base systems using matlab
	
	
	
	
	//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jac_hom_w_proj");
	
	mul_mat_vec_mp(*target_function_values, BED->jac_with_proj, (*curr_v_vars));
	
	//	print_point_to_screen_matlab(target_function_values,"target_f");
	vec_mulcomp_mp(*target_function_values_times_oneminus_s, *target_function_values, *one_minus_s);
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	// the product of the linears
	for (int jj=0; jj<BED->num_jac_equations; jj++) {
		
		//perform the $x$ evaluation
		set_one_mp(&(*linprod_x)->coord[jj]); // initialize to 1 for multiplication
		for (int ii=0; ii<BED->max_degree; ++ii) {
			dot_product_mp(&(*lin_func_vals)->entry[jj][ii],BED->starting_linears[jj][ii],*curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_mp(&(*linprod_x)->coord[jj], &(*linprod_x)->coord[jj], &(*lin_func_vals)->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_mp(*temp, BED->v_linears[jj], (*curr_v_vars));
		
		//now set the combined $x,v$ value
		mul_mp(&(*start_function_values)->coord[jj], &(*linprod_x)->coord[jj], *temp); // start_function_values = linprod_x(x)*linear(v)
		mul_mp(&(*linprod_times_gamma_s)->coord[jj], *gamma_s, &(*start_function_values)->coord[jj]); // sets the value linprod_x*gamma*s
		
		
		neg_mp( &Jp->entry[jj + offset][0], &(*target_function_values)->coord[jj]);  // Jp = -target
		mul_mp( *temp, BED->gamma, &(*start_function_values)->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_mp( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], *temp); // Jp = Jp + gamma*start = -target + gamma*start
	}
	
	
	// mix the start and target into final function value to return.
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		add_mp(&funcVals->coord[ii+offset], &(*target_function_values_times_oneminus_s)->coord[ii], &(*linprod_times_gamma_s)->coord[ii]);
	}
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE start and target WRT V
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	for (int ii=0; ii<BED->num_jac_equations; ii++) {
		for (int jj=0; jj<BED->num_v_vars; jj++) {
			mul_mp(*temp, &BED->v_linears[ii]->coord[jj], &(*linprod_x)->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_mp(*temp2, *temp, *gamma_s);                                         // temp2 = gamma*s*M_ij * linprod_x(x)
			
			mul_mp(*temp, *one_minus_s, &BED->jac_with_proj->entry[ii][jj]);        // temp = (1-s)*([Jf^T pi^T])_ij
			add_mp(&Jv->entry[ii+offset][jj+(BED->num_natural_vars+BED->num_synth_vars)], *temp, *temp2);       // Jv = temp + temp2
		}
	}
	
	
	//////////////////
	// now the x derivatives corresponding to the linprod start system
	//////////////////////
	
	
	
	// an implementation of the product rule
	for (int mm=0; mm< BED->num_jac_equations; mm++) {
		for (int kk=0; kk<(BED->num_natural_vars+BED->num_synth_vars); kk++) { // for each variable
			set_zero_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk]); // initialize to 0 for the sum
			
			for (int ii=0; ii<BED->max_degree; ++ii) { //  for each linear
				
				set_mp(*running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (int jj=0; jj<BED->max_degree; jj++) {
					if (jj!=ii) {
						mul_mp(*running_prod,*running_prod,&(*lin_func_vals)->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk],&(*linprod_derivative_wrt_x)->entry[mm][kk],*running_prod);
				
			}// re:ii
			
			dot_product_mp(*temp, BED->v_linears[mm], (*curr_v_vars));  // these two lines multiply by  (v_linear •	v)
			mul_mp(&(*linprod_derivative_wrt_x)->entry[mm][kk], &(*linprod_derivative_wrt_x)->entry[mm][kk], *temp);
		} // re: kk
	} // re: mm
	
	
	
	///////////////
	// DIFFERENTIATE THE derivative of the target jacobian system wrt $h,x_i$.
	/////////////////
	
	//  compute ∂/∂h
	
	// tempmat1 will hold the second derivatives, omitting the homvar column (hence the minus 1), and omitting all synth variables, too
	increase_size_mat_mp(*tempmat1,(*temp_function_values)->size,BED->num_natural_vars-1);
	(*tempmat1)->rows = (*temp_function_values)->size; (*tempmat1)->cols = BED->num_natural_vars-1;
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears;
	
	
	// this one is different from the x derivatives because h is the homogenizing variable for the problem,
	// and there are terms which contain only constants (pi) and the homvar.
	
	int current_variable_index = 0; // later we will loop from 1:num_natural_vars
	int local_offset = derivative_offset-1; // local_offset indexes into the vector containing second derivatives.
											//  the -1 is because we increment this at the beginning of the nested for loop
	
	for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each function
		
		int degree_deficiency = BED->max_degree+1 - BED->randomizer()->base_degree(base_func_index);// +1 is because maxdegree is deg of derivatives,
																							 //and base_degrees is deg of natural function before differentiation.
		int entry_counter = 0; // indexes into the column we will set.  easier this way.
		for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
			for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) { // it's essentially an upper triangular matrix passed through reshape(-,[],1)
				local_offset++;//  MUST always increment.   cannot 'continue' before this.
				
				if (jj==0) {// for this particular nullspace implementation calculation, we omit the hom var.  one would not do this if they needed ALL second derivatives.
					continue; //
				}
				
				if (ii!=current_variable_index && jj!=current_variable_index) {
					continue;
				}
				
				// g = h^n * function = h^n df_{curr_func}/dx_{curr_var}
				// dg/dh = n*h^(n-1)function + h^n dfunction/dh
				if (degree_deficiency==0) { // no need to further homogenize, n=0 in above formula
					set_mp(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
				}
				else{ // must homogenize, n>0 in above formula
					
					//					std::cout << "homogenizing deg def " << degree_deficiency << std::endl;
					//
					//					print_comp_matlab(&deriv_vals->coord[local_offset],"second_der");
					// first, homogenize against the derivative number of times according to deficiency.
					//   h^n (∂^2)f/(∂x∂x)
					set_mp(*temp,&(*deriv_vals)->coord[local_offset]); // seed, grab the first derivative
					for (int zz = 0; zz<degree_deficiency; zz++) {
						mul_mp(*temp,*temp,&(*curr_x_vars)->coord[0]); // loop to multiply
					}
					// temp = h^n*∂f/∂x
					
					
					// next, homogenize against the non-derivative part (function)
					//   n h^(n-1) f =
					//					print_comp_matlab(&deriv_vals->coord[base_func_index*BED->num_natural_vars+jj],"first_der");
					set_mp(*temp3, &(*deriv_vals)->coord[base_func_index*(BED->num_natural_vars+BED->num_synth_vars)+jj]); // seed
					for (int zz=0; zz<degree_deficiency-1; zz++) {  // i emphasize the -1 in upper limit for this loop (degree_deficiency-1)
						mul_mp(*temp3, *temp3, &(*curr_x_vars)->coord[0]); // loop multiply
					}
					
					set_zero_mp(*temp2);
					mpf_set_d((*temp2)->r ,degree_deficiency);
//					temp2->r = degree_deficiency; temp2->i = 0; // temp2 = n
					mul_mp(*temp3,*temp2,*temp3);  // multiply by the coefficient according to the degree.
											   // temp3 = n*h^{n-1} f
					
					// put 'em together
					add_mp(&(*tempmat1)->entry[base_func_index][entry_counter], *temp3, *temp);
					
				}
				entry_counter++;
			}
		}
	}
	// ok, now tempmat1 has in it all the derivatives of the derivatives of the system f.
	//  we are trying to get d/dx( [(R*Jf)^T | pi^T] v).
	
	if (!BED->randomizer()->is_square()) {
		mat_mul_mp(*tempmat2, *(BED->randomizer()->get_mat_mp()), *tempmat1); // randomize
		nonconj_transpose(*tempmat3,*tempmat2);
	}
	else
	{
		nonconj_transpose(*tempmat3,*tempmat1);
	}
	
	
	
	
	
	
	//	print_matrix_to_screen_matlab(tempmat3,"hom_deriv_mat");
	increase_size_mat_mp(*tempmat3,(*tempmat3)->rows, BED->num_v_vars); // nondestructive resize
	(*tempmat3)->cols = BED->num_v_vars;
	
	
	// homogenize the projection terms for the concatenated jacobian
	set_one_mp(*temp2);
	mpf_set_d((*temp2)->r, BED->max_degree);
//	temp2->r = BED->max_degree;
	for (int mm=0; mm<BED->max_degree-1; mm++)  // remember, max_degree is the highest degree occurring in the DERIVATIVES
		mul_mp(*temp2,*temp2,&(*curr_x_vars)->coord[0]); // we will use this in a line or two to homogenize
												   // temp2 = n*h^{n-1}
	
	
	for (int mm=0; mm<BED->num_projections; mm++) {
		for (int nn=1; nn<BED->num_natural_vars; nn++) {
			mul_mp(&(*tempmat3)->entry[nn-1][BED->randomizer()->num_rand_funcs()+mm], &BED->target_projection[mm]->coord[nn],*temp2);
			// entry = pi_mm[nn]*n*h^{n-1}
		}
	}
	

	
	
	mul_mat_vec_mp(*tempvec, *tempmat3, (*curr_v_vars)); // multiply against the $v$ variables
												   //	print_point_to_screen_matlab(tempvec,"diff_h_times_v");
	vec_mulcomp_mp(*tempvec, *tempvec, *one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable 0. (for all jac eqns), (with current time taken into account)
	
	// now, combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
	
	for (int mm=0; mm<BED->num_jac_equations; mm++) {
		mul_mp(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
		add_mp(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
	}
	
	// end diff WRT hom_var
	
	
	
	change_size_mat_mp(*tempmat1,BED->randomizer()->num_base_funcs(),BED->num_natural_vars+BED->num_synth_vars);
	
	// THE X VARS
	(*curr_v_vars)->size -= BED->num_projections; // truncate.  this is for optimization, not theoretical reasons
											   //
	for (current_variable_index=1; current_variable_index<BED->num_natural_vars+BED->num_synth_vars; current_variable_index++) {
		
		
		local_offset = derivative_offset-1; // initialize.  we count up from here
		for (int base_func_index=0; base_func_index<BED->randomizer()->num_base_funcs(); base_func_index++) { // iterate over each function
			
			int degree_deficiency = BED->max_degree+1 - BED->randomizer()->base_degree(base_func_index);
			int entry_counter = 0;
			for (int ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++){
				for (int jj=ii; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) {
					local_offset++;//  MUST always increment.
					
					if (ii==0) {
						continue;
					}
					
					if ( (ii!=current_variable_index) && (jj!=current_variable_index) ) {
						continue;
					}
					
					
					set_mp(&(*tempmat1)->entry[base_func_index][entry_counter], &(*deriv_vals)->coord[local_offset]);
					
					
					// g = h^n * f
					// dg/dx = h^n df/dx
					if (degree_deficiency>0){
						for (int zz = 0; zz<degree_deficiency; zz++) {
							mul_mp(&(*tempmat1)->entry[base_func_index][entry_counter],&(*tempmat1)->entry[base_func_index][entry_counter],&(*curr_x_vars)->coord[0]); // loop to multiply
						}
					}
					
					entry_counter++;
				}
			}
		}
		
		
		
		
		
		// ok, now tempmat1 has in it all the (homogenized) derivatives of the derivatives of the system f.
		//  we are trying to get d/dx( [(R*Jf)^T | h^d•pi^T] v) = ( R*d/dx(Jf) )^T v_truncate.
		
		if (!BED->randomizer()->is_square()) {
			mat_mul_mp(*tempmat2, *(BED->randomizer()->get_mat_mp()), *tempmat1); // randomize
			nonconj_transpose(*tempmat3,*tempmat2);
		}
		else
		{
			nonconj_transpose(*tempmat3,*tempmat1);
		}
		
		
		
		// no need to copy in the homogenized pi coefficients, because the derivatives wrt x are all 0.
		
		//		print_matrix_to_screen_matlab(tempmat3,"xxx");
		mul_mat_vec_mp(*tempvec, *tempmat3, (*curr_v_vars)); // recall v has been truncated to avoid unnecessary 0*number multiplication
													   //		print_point_to_screen_matlab(tempvec,"diff_x_times_v");
		vec_mulcomp_mp(*tempvec, *tempvec, *one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable 0. (for all jac eqns)
		
		// combine this derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		for (int mm=0; mm<BED->num_jac_equations; mm++) {
			mul_mp(*temp, *gamma_s, &(*linprod_derivative_wrt_x)->entry[mm][current_variable_index]);
			add_mp(&Jv->entry[mm+offset][current_variable_index], &(*tempvec)->coord[mm], *temp);
		}
		
	}//re: current_variable_index for numerical diff
	
	// END DIFF wrt x
	
	
	
	
	
	
	(*curr_v_vars)->size += BED->num_projections; // un-truncate.  the values of the trailing variables remain unaffected by previous truncation.
	
	
	
	//set the X PATCH values
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations;
	
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_mp(&funcVals->coord[ii+offset], &(*patchValues)->coord[ii]);
		
		// Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_natural_vars+BED->num_synth_vars; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &(*Jv_Patch)->entry[ii][jj]);
	}
	
	
	offset = BED->randomizer()->num_rand_funcs() + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 2701;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		deliberate_segfault();
	}
	
	// V patch
	set_one_mp(*temp2);
	dot_product_mp(*temp, BED->v_patch, (*curr_v_vars));
	sub_mp(&funcVals->coord[BED->num_variables-1], *temp, *temp2);  // f = patch*v-1
	
	for (int ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&Jv->entry[BED->num_variables-1][BED->num_natural_vars+ii], &BED->v_patch->coord[ii]);
	
	
	// finally, set parVals & parDer correctly
	
	change_size_point_mp(parVals, 1);  change_size_vec_mp(parDer, 1);
	parVals->size = parDer->size = 1;
	
	set_mp(&parVals->coord[0], pathVars); // s = t
	set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	
	if (BED->verbose_level()==10) {
        print_comp_matlab(pathVars,"t");
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
//        print_point_to_screen_matlab(curr_x_vars,"currxvars");
//        print_point_to_screen_matlab(curr_v_vars,"currvvars");
		print_point_to_screen_matlab(funcVals,"F_mp");
		//	print_point_to_screen_matlab(parVals,"parVals");
		//	print_point_to_screen_matlab(parDer,"parDer");
//		print_matrix_to_screen_matlab(Jv,"Jv_mp");
		//	print_matrix_to_screen_matlab(Jp,"Jp");
		
//        print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
		//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
//		print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");
		
	}
	
	
	
	
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize_mp(&dehommed,curr_x_vars);
	mpf_out_str (BED->FOUT, 10, 15, pathVars->r);
	fprintf(BED->FOUT," ");
	mpf_out_str (BED->FOUT, 10, 15, pathVars->i);
	fprintf(BED->FOUT," ");
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		mpf_out_str (BED->FOUT, 10, 15, dehommed->coord[ii].r);
		fprintf(BED->FOUT," ");
		mpf_out_str (BED->FOUT, 10, 15, dehommed->coord[ii].i);
		fprintf(BED->FOUT," ");
		
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_mp(dehommed);
#endif
	
	
	return 0;
}














int nullspacejac_dehom(point_d out_d, point_mp out_mp,
					   int *out_prec,
					   point_d in_d, point_mp in_mp,
					   int in_prec,
					   void const *ED_d, void const *ED_mp)
{
	
	
	
	*out_prec = in_prec;
	
	
	
	if (in_prec < 64)
	{ // compute out_d
		nullspacejac_eval_data_d *BED_d = (nullspacejac_eval_data_d *)ED_d;
		
		change_size_vec_d(out_d,BED_d->num_natural_vars-1);
		out_d->size = BED_d->num_natural_vars-1;
		
		
		for (int ii=0; ii<BED_d->num_natural_vars-1; ++ii) {
			div_d(&out_d->coord[ii],&in_d->coord[ii+1],&in_d->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		BED_d = NULL;
		
	}
	else
	{ // compute out_mp
		
		nullspacejac_eval_data_mp *BED_mp = (nullspacejac_eval_data_mp *)ED_mp;
		// set prec on out_mp
		setprec_point_mp(out_mp, *out_prec);
		
		change_size_vec_mp(out_mp,BED_mp->num_natural_vars-1);
		out_mp->size = BED_mp->num_natural_vars-1;
		
		for (int ii=0; ii<BED_mp->num_natural_vars-1; ++ii) {
			div_mp(&out_mp->coord[ii],&in_mp->coord[ii+1],&in_mp->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		BED_mp = NULL;
		
	}
	
	
	
	
	
	return 0;
}









int change_nullspacejac_eval_prec(void const *ED, int new_prec)
{
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii, jj;
	
	
	if (new_prec != BED->curr_prec){
		// change the precision for the patch
		changePatchPrec_mp(new_prec, &BED->patch);
		
		if (BED->verbose_level() >=8)
			printf("prec  %lu\t-->\t%d\n",BED->curr_prec, new_prec);
		
		BED->SLP->precision = new_prec;
		
		change_prec_prog_deriv(BED->SLP_derivative,new_prec);
		
		BED->curr_prec = new_prec;
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		BED->randomizer()->change_prec(new_prec);

		
		BED->temp_vars.change_prec(new_prec);
		
		for (ii=0; ii<BED->num_additional_linears; ii++) {
			change_prec_point_mp(BED->additional_linears_terminal[ii],new_prec);
			vec_cp_mp(BED->additional_linears_terminal[ii],BED->additional_linears_terminal_full_prec[ii]);
			
			
			change_prec_point_mp(BED->additional_linears_starting[ii],new_prec);
			vec_cp_mp(BED->additional_linears_starting[ii],BED->additional_linears_starting_full_prec[ii]);
		}
		
		
		if (BED->side()==nullspace_handedness::LEFT) {
			for (ii=0; ii<BED->num_jac_equations; ++ii) {
				for (jj=0; jj<BED->max_degree; jj++) {
					change_prec_point_mp(BED->starting_linears[ii][jj],new_prec);
					vec_cp_mp(BED->starting_linears[ii][jj], BED->starting_linears_full_prec[ii][jj]);
				}
			}
		}
		else{
			for (ii=0; ii<BED->randomizer()->num_rand_funcs(); ++ii) {
				for (jj=0; jj<BED->randomizer()->randomized_degree(ii)-1; jj++) {
					change_prec_point_mp(BED->starting_linears[ii][jj],new_prec);
					vec_cp_mp(BED->starting_linears[ii][jj], BED->starting_linears_full_prec[ii][jj]);
				}
			}
		}
		
		
		
		for (ii=0; ii<BED->num_v_linears; ii++) {
			change_prec_point_mp(BED->v_linears[ii],new_prec);
			vec_cp_mp(BED->v_linears[ii],BED->v_linears_full_prec[ii]);
		}
		
		
		change_prec_point_mp(BED->v_patch,new_prec);
		vec_cp_mp(BED->v_patch,BED->v_patch_full_prec);
		
		change_prec_mat_mp(BED->jac_with_proj,new_prec);
		mat_cp_mp(BED->jac_with_proj,BED->jac_with_proj_full_prec);
	}
	
	
	return 0;
}






int check_issoln_nullspacejac_d(endgame_data_t *EG,
								tracker_config_t *T,
								void const *ED)
{
	nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	

	
	double n1, n2, max_rat;
	point_d f;
	eval_struct_d e;
	//
	//	mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
	init_point_d(f, 1);
	init_eval_struct_d(e,0, 0, 0);
	
	max_rat = T->ratioTol;
	
	// setup threshold based on given threshold and precision
	//	if (num_digits > 300)
	//		num_digits = 300;
	//	num_digits -= 2;
	double tol = MAX(T->funcResTol, 1e-10);
	
	
	if (EG->prec>=64){
		vec_d terminal_pt;  init_vec_d(terminal_pt,1);
		vec_mp_to_d(terminal_pt,EG->PD_mp.point);
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, terminal_pt, EG->PD_d.time, BED->SLP);
		clear_vec_d(terminal_pt);
	}
	else{
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, BED->SLP);
		
		
	}
	
	
	if (EG->last_approx_prec>=64) {
		vec_d prev_pt;  init_vec_d(prev_pt,1);
		vec_mp_to_d(prev_pt,EG->last_approx_mp);
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, prev_pt, EG->PD_d.time, BED->SLP);
		clear_vec_d(prev_pt);}
	else{
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_d, EG->PD_d.time, BED->SLP);
	}
	
	
	
	
	
	
	
	
	//	print_point_to_screen_matlab(EG->PD_d.point,"soln");
	//	print_point_to_screen_matlab(e.funcVals,"howfaroff");	// compare the function values
	int isSoln = 1;
	for (int ii = 0; (ii < BED->SLP->numFuncs) && isSoln; ii++)
	{
		n1 = d_abs_d( &e.funcVals->coord[ii]); // corresponds to final point
		n2 = d_abs_d( &f->coord[ii]); // corresponds to the previous point
		
		
		if (tol <= n1 && n1 <= n2)
		{ // compare ratio
			if (n1 > max_rat * n2){ // seriously what is the point of this
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (d) 1 coord %d\n",ii);
			}
		}
		else if (tol <= n2 && n2 <= n1)
		{ // compare ratio
			if (n2 > max_rat * n1){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (d) 2 coord %d\n",ii);
			}
		}
	}
	
	
	
	
	if (!isSoln) {
		
		print_point_to_screen_matlab(e.funcVals,"terminal_func_vals");
		print_point_to_screen_matlab(f,"prev");
		
		printf("tol was %le\nmax_rat was %le\n",tol,max_rat);
	}
	
	
	

	
	clear_eval_struct_d(e);
	clear_vec_d(f);
	
	
	return isSoln;
	
}


int check_issoln_nullspacejac_mp(endgame_data_t *EG,
								 tracker_config_t *T,
								 void const *ED)
{
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii;
	
	vec_mp curr_v_vars; init_vec_mp(curr_v_vars, BED->num_v_vars);
	curr_v_vars->size = BED->num_v_vars;
	
	vec_mp curr_x_vars; init_vec_mp(curr_x_vars, BED->num_natural_vars);
	curr_x_vars->size = BED->num_natural_vars;
	
	
	
	for (ii = 0; ii < T->numVars; ii++)
	{
		if (!(mpfr_number_p(EG->PD_mp.point->coord[ii].r) && mpfr_number_p(EG->PD_mp.point->coord[ii].i)))
		{
			printf("got not a number\n");
			print_point_to_screen_matlab(EG->PD_mp.point,"bad solution");
			return 0;
		}
	}
	
	
	
	mpf_t n1, n2, zero_thresh, max_rat;
	mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
	
	point_mp f; init_point_mp(f, 1);f->size = 1;
	eval_struct_mp e; init_eval_struct_mp(e, 0, 0, 0);
	
	mpf_set_d(max_rat, T->ratioTol);
	
	
	int num_digits = prec_to_digits((int) mpf_get_default_prec());
	// setup threshold based on given threshold and precision
	if (num_digits > 300)
		num_digits = 300;
	num_digits -= 4;
	double tol = MAX(T->funcResTol, pow(10,-num_digits));
	mpf_set_d(zero_thresh, tol);
	
	
	for (ii=0; ii<BED->num_natural_vars; ii++)
		set_mp(&curr_x_vars->coord[ii], &EG->PD_mp.point->coord[ii]);
	
	for (ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&curr_v_vars->coord[ii], &EG->PD_mp.point->coord[ii+BED->num_natural_vars]);
	
	
	
	//this one guaranteed by entry condition
	//	lin_to_lin_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, ED);
	evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, BED->SLP);
	
	//	print_point_to_screen_matlab(EG->PD_mp.point,"soln");
	//	print_point_to_screen_matlab(e.funcVals,"howfaroff");
	
	if (EG->last_approx_prec < 64) { // copy to _mp
		point_d_to_mp(EG->last_approx_mp, EG->last_approx_d);
	}
	
	evalProg_mp(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_mp, EG->PD_mp.time, BED->SLP);
	//	lin_to_lin_eval_mp(f,          e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_mp, EG->PD_mp.time, ED);
	// compare the function values
	int isSoln = 1;
	for (ii = 0; ii < BED->SLP->numFuncs && isSoln; ii++)
	{
		mpf_abs_mp(n1, &e.funcVals->coord[ii]);
		mpf_abs_mp(n2, &f->coord[ii]);
		
		//		mpf_out_str(NULL,10,9,n1);
		
		if ( (mpf_cmp(zero_thresh, n1) <= 0) &&  (mpf_cmp(n1, n2) <= 0) )
		{ // compare ratio
			mpf_mul(n2, max_rat, n2);
			if (mpf_cmp(n1, n2) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 1\n");
				mpf_out_str(NULL,10,0,max_rat);
			}
		}
		else if ( (mpf_cmp(zero_thresh, n2) <= 0) &&  (mpf_cmp(n2, n1) <= 0) )
		{ // compare ratio
			mpf_mul(n1, max_rat, n1);
			if (mpf_cmp(n2, n1) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 2\n");
			}
		}
	}
	
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	clear_vec_mp(f);
	clear_vec_mp(curr_x_vars);
	clear_vec_mp(curr_v_vars);
	
	
	
	return isSoln;
	
}







int check_isstart_nullspacejac_d(vec_d testpoint,
								 tracker_config_t *T,
								 void const *ED)
{
	nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED;
	eval_struct_d e;
	init_eval_struct_d(e,0, 0, 0);
	
	comp_d time;
	set_one_d(time);
	
	
	double tol = T->funcResTol;
	
	BED->evaluator_function_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, testpoint, time, ED);
	
	int isSoln = 1;
	
	for (int ii = 0; (ii < e.funcVals->size) && isSoln; ii++) // function by function
	{
		if (tol <= d_abs_d( &e.funcVals->coord[ii])){ // compare
			isSoln = 0;
			print_point_to_screen_matlab(testpoint,"invalid_startpoint");
			print_point_to_screen_matlab(e.funcVals,"start_residual");
			std::cout << tol << " is tol" << std::endl;
			break;
		}
		
	}
	
	
	clear_eval_struct_d(e);
	
	return isSoln;
	
}


int check_isstart_nullspacejac_mp(vec_mp testpoint,
								 tracker_config_t *T,
								 void const *ED)
{
	
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED;
	
	
	eval_struct_mp e;
	init_eval_struct_mp(e,0, 0, 0);
	
	comp_mp time; init_mp(time);
	set_one_mp(time);
	
	
	double tol = T->funcResTol;
	
	BED->evaluator_function_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, testpoint, time, ED);
	
	int isSoln = 1;
	
	for (int ii = 0; (ii < e.funcVals->size) && isSoln; ii++) // function by function
	{
		if (tol <= d_abs_mp( &e.funcVals->coord[ii])){ // compare
			isSoln = 0;
			print_point_to_screen_matlab(testpoint,"invalid_startpoint");
			print_point_to_screen_matlab(e.funcVals,"start_residual");
			break;
		}
		
	}
	
	clear_mp(time);
	
	clear_eval_struct_mp(e);
	
	return isSoln;
	
}



void check_nullspace_evaluator(point_mp current_values,
							   void const *ED)
{
	int ii;
	printf("checking homogeneousness of double evaluator\n");
	nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
																	//initialize
	eval_struct_d e_d; init_eval_struct_d(e_d, 0, 0, 0);
	eval_struct_d e_d2; init_eval_struct_d(e_d2, 0, 0, 0);
	
	
	
	
	comp_d time;
	
	
	
	
	point_d tempvec;  init_point_d(tempvec,0);
	point_d tempvec2;  init_point_d(tempvec2,0);
	make_vec_random_d(tempvec, current_values->size);
	
	
	comp_d lambda;
	lambda->r = 2; lambda->i = 0;
//	get_comp_rand_d(lambda);
	
	
	vec_cp_d(tempvec2, tempvec);
	for (ii=0; ii<BED->num_natural_vars+BED->num_synth_vars; ii++) {
		mul_d(&tempvec2->coord[ii],&tempvec->coord[ii],lambda);
	}
	
	print_point_to_screen_matlab(tempvec,"input");
	printf("lambda = %lf+1i*%lf\n",lambda->r, lambda->i);
	
	set_zero_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"f1");
	print_point_to_screen_matlab(e_d2.funcVals,"f2");
	print_matrix_to_screen_matlab(e_d.Jv,"j1");
	print_matrix_to_screen_matlab(e_d2.Jv,"j2");
	
	set_one_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"g1");
	print_point_to_screen_matlab(e_d2.funcVals,"g2");
	print_matrix_to_screen_matlab(e_d.Jv,"k1");
	print_matrix_to_screen_matlab(e_d2.Jv,"k2");
	
	
	
	
	
	
	
	

	vec_cp_d(tempvec2, tempvec);
	for (ii=BED->num_natural_vars+BED->num_synth_vars; ii<BED->num_variables; ii++) {
		mul_d(&tempvec2->coord[ii],&tempvec->coord[ii],lambda);
	}
	
	print_point_to_screen_matlab(tempvec,"input");
	printf("lambda = %lf+1i*%lf\n",lambda->r, lambda->i);
	
	set_zero_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"f1");
	print_point_to_screen_matlab(e_d2.funcVals,"f2");
	print_matrix_to_screen_matlab(e_d.Jv,"j1");
	print_matrix_to_screen_matlab(e_d2.Jv,"j2");
	
	set_one_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"g1");
	print_point_to_screen_matlab(e_d2.funcVals,"g2");
	print_matrix_to_screen_matlab(e_d.Jv,"k1");
	print_matrix_to_screen_matlab(e_d2.Jv,"k2");
	
	
	
	
	vec_cp_d(tempvec2, tempvec);
	for (ii=0; ii<BED->num_variables; ii++) {
		mul_d(&tempvec2->coord[ii],&tempvec->coord[ii],lambda);
	}
	
	print_point_to_screen_matlab(tempvec,"input");
	printf("lambda = %lf+1i*%lf\n",lambda->r, lambda->i);
	
	set_zero_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"f1");
	print_point_to_screen_matlab(e_d2.funcVals,"f2");
	print_matrix_to_screen_matlab(e_d.Jv,"j1");
	print_matrix_to_screen_matlab(e_d2.Jv,"j2");
	
	set_one_d(time);
	BED->evaluator_function_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, time, ED);
	BED->evaluator_function_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec2, time, ED);
	
	print_point_to_screen_matlab(e_d.funcVals,"g1");
	print_point_to_screen_matlab(e_d2.funcVals,"g2");
	print_matrix_to_screen_matlab(e_d.Jv,"k1");
	print_matrix_to_screen_matlab(e_d2.Jv,"k2");
	
	
	
	
	
	
	
	clear_eval_struct_d(e_d);
	clear_eval_struct_d(e_d2);
	clear_vec_d(tempvec);
	
	return;
	
}
































