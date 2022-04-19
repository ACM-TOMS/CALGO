#include "nag/solvers/sphere_intersection.hpp"


void SphereConfiguration::set_memory(SolverConfiguration & solve_options)
{
	
	//TODO: should i assume here that the input file is already parsed??  (as i already do)
	this->MPType = solve_options.T.MPType;
	solve_options.T.numVars = setupProg(SLP, solve_options.T.Precision, solve_options.T.MPType);
	//make randomizer matrix here
	SLP_memory.capture_globals();
	SLP_memory.set_globals_null();
	have_mem = true;
	
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
}












void sphere_eval_data_mp::init()
{
	this->is_solution_checker_d = &check_issoln_sphere_d;
	this->is_solution_checker_mp = &check_issoln_sphere_mp;
	this->evaluator_function_d = &sphere_eval_d;
	this->evaluator_function_mp = &sphere_eval_mp;
	this->precision_changer = &change_sphere_eval_prec;
	this->dehomogenizer = &sphere_dehom;
	
	
	
	this->num_static_linears = 0; // we will copy these out of the witness set.
	static_linear = static_linear_full_prec = NULL;
	
	
	starting_linear = (vec_mp *) br_malloc(2*sizeof(vec_mp));;
	for (int ii=0; ii<2; ii++) {
		init_vec_mp(starting_linear[ii],0);
	}
	
	init_vec_mp(center, 0);
	init_mp(radius);
	
	if (MPType==2) {
		init_vec_mp2(center_full_prec,0,1024);
		init_mp2(radius_full_prec,1024);
		
		starting_linear_full_prec = (vec_mp *) br_malloc(2*sizeof(vec_mp));
		for (int ii=0; ii<2; ii++) {
			init_vec_mp2(starting_linear_full_prec[ii],0,1024);
		}
		
		init_mp2(two_full_prec,1024);
		set_zero_mp(two_full_prec);
		mpf_set_str(two_full_prec->r, "2.0", 10);
	}
	
	init_mp(two);
	set_zero_mp(two);
	mpf_set_str(two->r, "2.0", 10);
	
    num_natural_vars = 0;
}


int sphere_eval_data_mp::send(ParallelismConfig & mpi_config)
{
	
	int solver_choice = SPHERE_SOLVER;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	SolverMultiplePrecision::send(mpi_config);
	
	
	// now can actually send the data.
	
	
	int *buffer = new int[2];
	
	buffer[0] = num_natural_vars;
	buffer[1] = num_static_linears;
	
	MPI_Bcast(buffer,2,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	delete[] buffer;
	
	
	if (this->MPType==2){
		for (int ii=0; ii<num_static_linears; ii++) {
			bcast_vec_mp(static_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
		}
		
		for (int ii=0; ii<2; ii++) {
			bcast_vec_mp(starting_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
		}
		
		bcast_vec_mp(center_full_prec, mpi_config.id(), mpi_config.head());
		
		bcast_comp_mp(radius_full_prec, mpi_config.id(), mpi_config.head());
		
		
	}
	else {
		for (int ii=0; ii<num_static_linears; ii++) {
			bcast_vec_mp(static_linear[ii], mpi_config.id(), mpi_config.head());
		}
		
		for (int ii=0; ii<2; ii++) {
			bcast_vec_mp(starting_linear[ii], mpi_config.id(), mpi_config.head());
		}
		
		bcast_vec_mp(center, mpi_config.id(), mpi_config.head());
		
		bcast_comp_mp(radius, mpi_config.id(), mpi_config.head());
		
	}
	
	
	return SUCCESSFUL;
}

int sphere_eval_data_mp::receive(ParallelismConfig & mpi_config)
{
	
	int *buffer = new int[2];
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	
	if (buffer[0] != SPHERE_SOLVER) {
		std::cout << "worker failed to confirm it is receiving the SPHERE_SOLVER type eval data" << std::endl;
		mpi_config.abort(777);
	}
	
	SolverMultiplePrecision::receive(mpi_config);
	
	
	// now can actually receive the data from whoever.
	MPI_Bcast(buffer, 2, MPI_INT, mpi_config.head(), mpi_config.comm());
	
	num_natural_vars = buffer[0];
	num_static_linears = buffer[1];
	
	delete[] buffer;
	
	//starting linears already created and initted
	static_linear = (vec_mp *) br_malloc(num_static_linears*sizeof(vec_mp));
	
	if (this->MPType==2) {
		static_linear_full_prec = (vec_mp *) br_malloc(num_static_linears*sizeof(vec_mp));
		
		for (int ii=0; ii<num_static_linears; ii++) {
			
			init_vec_mp(static_linear[ii],1);
			init_vec_mp2(static_linear_full_prec[ii],1,1024);
			
			bcast_vec_mp(static_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
			
			vec_cp_mp(static_linear[ii],static_linear_full_prec[ii]);
			
		}
		
		for (int ii=0; ii<2; ii++) {
			
			bcast_vec_mp(starting_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
			
			vec_cp_mp(starting_linear[ii],starting_linear_full_prec[ii]);
		}
		
		bcast_vec_mp(center_full_prec, mpi_config.id(), mpi_config.head());
		
		bcast_comp_mp(radius_full_prec, mpi_config.id(), mpi_config.head());
		
		vec_cp_mp(center, center_full_prec);
		set_mp(radius, radius_full_prec);
	}
	else{ // MPType == 1
		for (int ii=0; ii<num_static_linears; ii++) {
			init_vec_mp(static_linear[ii],1);
			bcast_vec_mp(static_linear[ii], mpi_config.id(), mpi_config.head());
		}
		
		for (int ii=0; ii<2; ii++) {
			bcast_vec_mp(starting_linear[ii], mpi_config.id(), mpi_config.head());
		}
		
		bcast_vec_mp(center, mpi_config.id(), mpi_config.head());
		bcast_comp_mp(radius, mpi_config.id(), mpi_config.head());
	}
	
	
	return SUCCESSFUL;
	
}

int sphere_eval_data_mp::setup(SphereConfiguration & config,
							   const WitnessSet & W,
							   SolverConfiguration & solve_options)
{
	
	
	if (config.randomizer().use_count()==0) {
		std::cout << "don't have randomizer set up!" << std::endl;
		br_exit(-97621);
	}
	
	if (!config.have_mem) {
		std::cout << "don't have memory!" << std::endl;
		br_exit(-3231);
	}
	
	this->SLP_memory = config.SLP_memory;
	
	num_natural_vars = W.num_natural_variables();
	num_variables = W.num_variables();
	
	
	// set up the vectors to hold the linears.
	if (this->num_static_linears==0) {
		static_linear = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
	}
	else
	{
		static_linear = (vec_mp *) br_realloc(static_linear, W.num_linears()*sizeof(vec_mp));
	}
	
	for (unsigned int ii=0; ii<W.num_linears(); ii++) {
		init_vec_mp(static_linear[ii],0);
		vec_cp_mp(static_linear[ii],W.linear(ii));
	}
	
	
	
	set_mp(this->radius, config.radius);
	vec_cp_mp(this->center, config.center);
	if (this->center->size < W.num_variables()) {
        int old_size = this->center->size;
        increase_size_vec_mp(this->center, W.num_variables());
        this->center->size = W.num_variables();
        
        for (int ii=old_size; ii<W.num_variables(); ii++) {
            set_zero_mp(&this->center->coord[ii]);
        }
    }
    
    
	if (this->MPType==2) {
		
		set_mp(this->radius_full_prec, config.radius);
		vec_cp_mp(this->center_full_prec, config.center);
        if (this->center_full_prec->size < W.num_variables()) {
            int old_size = this->center_full_prec->size;
            increase_size_vec_mp(this->center_full_prec, W.num_variables());
            this->center_full_prec->size = W.num_variables();
            
            for (int ii=old_size; ii<W.num_variables(); ii++) {
                set_zero_mp(&this->center_full_prec->coord[ii]);
            }
        }
        
        
		
		if (this->num_static_linears==0) {
			static_linear_full_prec = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
		}
		else
		{
			static_linear_full_prec = (vec_mp *) br_realloc(static_linear_full_prec, W.num_linears()*sizeof(vec_mp));
		}
		
		for (unsigned int ii=0; ii<W.num_linears(); ii++) {
			init_vec_mp2(static_linear_full_prec[ii],0,1024);
			
			vec_cp_mp(static_linear_full_prec[ii], W.linear(ii));
		}
	}
	
	this->num_static_linears = W.num_linears();
	
	
	
	
	
	for (int ii=0; ii<2; ii++) {
		vec_cp_mp( this->starting_linear[ii], config.starting_linear[ii]);
		
		if (MPType==2) {
			vec_cp_mp(this->starting_linear_full_prec[ii], config.starting_linear[ii]);
		}
	}
	
	
	
	
	
	// the usual
	
	verbose_level(solve_options.verbose_level());
	
	SolverMultiplePrecision::setup(config.SLP, config.randomizer());
	
	generic_setup_patch(&patch,W);
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_mp(this->gamma); // set gamma to be random complex value
	else{
		set_one_mp(this->gamma);
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
	
	

	
	return 0;
}























void sphere_eval_data_d::init()
{
	
	if (this->MPType==2)
	{
		this->BED_mp = new sphere_eval_data_mp(2);
		SolverDoublePrecision::BED_mp = this->BED_mp;                   //   <---------  you gotta do this cuz of masking problems.
	}
	else{
		this->BED_mp = NULL;
	}
	
	
	this->is_solution_checker_d = &check_issoln_sphere_d;
	this->is_solution_checker_mp = &check_issoln_sphere_mp;
	this->evaluator_function_d = &sphere_eval_d;
	this->evaluator_function_mp = &sphere_eval_mp;
	this->precision_changer = &change_sphere_eval_prec;
	this->dehomogenizer = &sphere_dehom;
	
	starting_linear = (vec_d *) br_malloc(2*sizeof(vec_d));
	init_vec_d(starting_linear[0],0);
	init_vec_d(starting_linear[1],0);
	num_static_linears = 0;
	static_linear = NULL;
	
	two->r = double(2);
	two->i = double(0);
	
	init_vec_d(center,0);
	
    num_natural_vars = 0;
}




int sphere_eval_data_d::send(ParallelismConfig & mpi_config)
{
	
	int solver_choice = SPHERE_SOLVER;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
    
    if (this->MPType==2) {
		this->BED_mp->send(mpi_config);
	}
    
    
	//send the base class stuff.
	SolverDoublePrecision::send(mpi_config);
	
	
	
	int *buffer = new int[2];
	
	buffer[0] = num_natural_vars;
	buffer[1] = num_static_linears;
	
	MPI_Bcast(buffer,2,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	delete[] buffer;
	
	for (int ii=0; ii<num_static_linears; ii++) {
		bcast_vec_d(static_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	for (int ii=0; ii<2; ii++) {
		bcast_vec_d(starting_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	bcast_vec_d(center, mpi_config.id(), mpi_config.head());
	
	bcast_comp_d(radius, mpi_config.id(), mpi_config.head());
	
	
	
	return SUCCESSFUL;
}

int sphere_eval_data_d::receive(ParallelismConfig & mpi_config)
{
	
    int *buffer = new int[1];
	
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	if (buffer[0] != SPHERE_SOLVER){
		std::cout << "worker failed to confirm it is receiving the double SPHERE_SOLVER type eval data" << std::endl;
		mpi_config.abort(777);
	}
    
    
    
	if (this->MPType==2) {
		this->BED_mp->receive(mpi_config);
	}
    
    
	
	SolverDoublePrecision::receive(mpi_config);
	
	
	
	// now can actually receive the data from whoever.
	MPI_Bcast(buffer, 2, MPI_INT, mpi_config.head(), mpi_config.comm());
	
	num_natural_vars = buffer[0];
	num_static_linears = buffer[1];
	
	delete[] buffer;
	
	
	static_linear = (vec_d *) br_malloc(num_static_linears*sizeof(vec_d));
	
	for (int ii=0; ii<num_static_linears; ii++) {
		init_vec_d(static_linear[ii],1);
		bcast_vec_d(static_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	for (int ii=0; ii<2; ii++) {
		bcast_vec_d(starting_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	bcast_vec_d(center, mpi_config.id(), mpi_config.head());
	bcast_comp_d(radius, mpi_config.id(), mpi_config.head());
	
	return SUCCESSFUL;
}


int sphere_eval_data_d::setup(SphereConfiguration & config,
							  const WitnessSet & W,
							  SolverConfiguration & solve_options)
{
	
	SLP_memory = config.SLP_memory;
	
	// set up the vectors to hold the two linears.
	
	
	mp_to_d(radius, config.radius);
	vec_mp_to_d(center, config.center);
    
    if (this->center->size < W.num_variables()) {
        int old_size = this->center->size;
        increase_size_vec_d(this->center, W.num_variables());
        this->center->size = W.num_variables();
        
        for (int ii=old_size; ii<W.num_variables(); ii++) {
            set_zero_d(&this->center->coord[ii]);
        }
    }
    
    
	if (this->num_static_linears==0) {
		static_linear = (vec_d *) br_malloc(W.num_linears()*sizeof(vec_d));
	}
	else
	{
		static_linear = (vec_d *) br_realloc(static_linear, W.num_linears()*sizeof(vec_d));
	}
	
	for (unsigned int ii=0; ii<W.num_linears(); ii++) {
		init_vec_d(static_linear[ii],0);
		vec_mp_to_d(static_linear[ii],W.linear(ii));
	}
	num_static_linears = W.num_linears();
	
    num_natural_vars = W.num_natural_variables();
	num_variables = W.num_variables();
	
	
	
	for (int ii=0; ii<2; ii++) {
		vec_mp_to_d(starting_linear[ii],config.starting_linear[ii]);
	}
	
	verbose_level(solve_options.verbose_level());
	
	SolverDoublePrecision::setup(config.SLP, config.randomizer());
	
	generic_setup_patch(&patch,W);
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_d(this->gamma); // set gamma to be random complex value
	else
		set_one_d(this->gamma);
	
	
	
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(config, W, solve_options);
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
	}
	
	return 0;
}



















int sphere_solver_master_entry_point(const WitnessSet						&W, // carries with it the start points, and the linears.
									 SolverOutput							&solve_out, // new data goes in here
									 SphereConfiguration &		config,
									 SolverConfiguration		& solve_options)
{
	
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(SPHERE_SOLVER);
	}
	
	sphere_eval_data_d *ED_d = NULL;
	sphere_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new sphere_eval_data_d(0);
			
			ED_d->setup(config,
						W,
						solve_options);
			break;
			
		case 1:
			ED_mp = new sphere_eval_data_mp(1);
			
			ED_mp->setup(config,
						 W,
						 solve_options);
//			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //    <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new sphere_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(config,
						W,
						solve_options);
			
			
			
			adjust_tracker_AMP(& (solve_options.T), W.num_variables());
			
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
		{
			std::stringstream printme;
			printme << "trying to use an MPType which is not 0, 1, or 2.  yours is " << solve_options.T.MPType << std::endl;
			throw std::logic_error(printme.str());
			break;
		}
	}
	
	
	for (unsigned int jj=0; jj<W.num_linears(); jj++)
	{
		solve_out.add_linear(W.linear(jj));
	}
	
	return SUCCESSFUL;
	
}





int sphere_slave_entry_point(SolverConfiguration & solve_options)
{
	
	
	// already received the flag which indicated that this worker is going to be performing the nullspace calculation.
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	int *settings_buffer = (int *) br_malloc(2*sizeof(int));
	MPI_Bcast(settings_buffer,2,MPI_INT, 0,solve_options.comm());
	solve_options.robust = settings_buffer[0];
	solve_options.use_gamma_trick = settings_buffer[1];
	free(settings_buffer);
	
	sphere_eval_data_d *ED_d = NULL;
	sphere_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new sphere_eval_data_d(0);
			ED_d->receive(solve_options);
			break;
			
		case 1:
			ED_mp = new sphere_eval_data_mp(1);
			ED_mp->receive(solve_options);
			
			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //   <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new sphere_eval_data_d(2);
			ED_mp = ED_d->BED_mp;
			ED_d->receive(solve_options);
			
			
			
			
			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //   <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		default:
			break;
	}
	
    
	
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "output",
                        &midOUT, "midpath_data");
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	worker_tracker_loop(&trackCount, OUT, midOUT,
						ED_d, ED_mp,
						solve_options);
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	//clear data
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
	
	
	return SUCCESSFUL;
}







int sphere_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp,
				  point_d current_variable_values, comp_d pathVars,
				  void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	sphere_eval_data_d *BED = (sphere_eval_data_d *)ED; // to avoid having to cast every time
	
	
	BED->SLP_memory.set_globals_to_this();
	
	int ii, jj, mm; // counters
	int offset;
	comp_d one_minus_s, gamma_s;
	comp_d temp, temp2;
	comp_d func_val_sphere, func_val_start;
	
	set_one_d(one_minus_s);
	sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	vec_d patchValues; init_vec_d(patchValues, 0);
	vec_d temp_function_values; init_vec_d(temp_function_values,0);
	vec_d AtimesF; init_vec_d(AtimesF,BED->randomizer()->num_rand_funcs()); AtimesF->size = BED->randomizer()->num_rand_funcs();// declare  // initialize
	
	
	
	mat_d temp_jacobian_functions; init_mat_d(temp_jacobian_functions,BED->randomizer()->num_base_funcs(),BED->num_variables);
	temp_jacobian_functions->rows = BED->randomizer()->num_base_funcs(); temp_jacobian_functions->cols = BED->num_variables;
	mat_d temp_jacobian_parameters; init_mat_d(temp_jacobian_parameters,0,0);
	mat_d Jv_Patch; init_mat_d(Jv_Patch, 0, 0);
	mat_d AtimesJ; init_mat_d(AtimesJ,BED->randomizer()->num_rand_funcs(),BED->num_variables);
	AtimesJ->rows = BED->randomizer()->num_rand_funcs(); AtimesJ->cols = BED->num_variables;
	
	
	//set the sizes
	change_size_vec_d(funcVals,BED->num_variables); funcVals->size = BED->num_variables;
	change_size_mat_d(Jv, BED->num_variables, BED->num_variables); Jv->rows = Jv->cols = BED->num_variables; //  -> this should be square!!!
	
	for (ii=0; ii<BED->num_variables; ii++)
		for (jj=0; jj<BED->num_variables; jj++)
			set_zero_d(&Jv->entry[ii][jj]);
	
	
	
	// evaluate the SLP to get the system's whatnot.
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	
	// evaluate the patch
	patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	// note that you can only really do this AFTER you are done calling other evaluators.
	// set parVals & parDer correctly
	
	// i.e. these must remain here, or below.  \/
	change_size_point_d(parVals, 1);
	change_size_vec_d(parDer, 1);
	change_size_mat_d(Jp, BED->num_variables, 1); Jp->rows = BED->num_variables; Jp->cols = 1;
	for (ii=0; ii<BED->num_variables; ii++)
		set_zero_d(&Jp->entry[ii][0]);
	
	
	parVals->size = parDer->size = 1;
	set_d(&parVals->coord[0], pathVars); // s = t
	set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	///////////////////////////
	//
	// the original (randomized) functions.
	//
	///////////////////////////////////
	
	BED->randomizer()->randomize(AtimesF,AtimesJ,temp_function_values,temp_jacobian_functions,&current_variable_values->coord[0]);

	
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real orthogonal) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	for (ii = 0; ii < BED->randomizer()->num_rand_funcs(); ii++)
		for (jj = 0; jj < BED->num_variables; jj++)
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	//Jp is 0 for the equations.
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	///////////////////
	//
	//  the sphere equation.
	//
	//////////////////////////
	
	offset = BED->randomizer()->num_rand_funcs();
	
	mul_d(func_val_sphere, BED->radius, BED->radius);
	neg_d(func_val_sphere, func_val_sphere);
	mul_d(func_val_sphere, func_val_sphere, &current_variable_values->coord[0]);
	mul_d(func_val_sphere, func_val_sphere, &current_variable_values->coord[0]);
	//f_sph = -r^2*h^2
	
	
	
	for (int ii=1; ii<BED->num_natural_vars; ii++) {
		mul_d(temp2, &BED->center->coord[ii-1], &current_variable_values->coord[0]); // temp2 = c_{i-1}*h
		
		sub_d(temp, &current_variable_values->coord[ii], temp2);  // temp = x_i - h*c_{i-1}
		mul_d(temp2, temp, temp);                                 // temp2 = (x_i - h*c_{i-1})^2
		add_d(func_val_sphere, func_val_sphere, temp2);           // f_sph += (x_i - h*c_{i-1})^2
	}
	
	
	
	set_one_d(func_val_start);
	for (mm=0; mm<2; ++mm) {
		dot_product_d(temp, BED->starting_linear[mm], current_variable_values);
		mul_d(func_val_start, func_val_start, temp);
		//f_start *= L_i (x)
	}
	
	
	// combine the function values
	mul_d(temp, one_minus_s, func_val_sphere);
	mul_d(temp2, gamma_s, func_val_start);
	add_d(&funcVals->coord[offset], temp, temp2);
	// f = (1-t) f_sph + gamma t f_start
	
	
	
	
	//// / / / / / /    now the derivatives wrt x
	
	//  first we store the derivatives of the target function, the sphere.  then we will add the part for the linear product start.
	
	
	
	//ddx for sphere
	
	
	
	for (int ii=1; ii<BED->num_natural_vars; ii++) {
		mul_d(temp2, &BED->center->coord[ii-1], &current_variable_values->coord[0]); // temp2 = c_{i-1}*h
		sub_d(temp, &current_variable_values->coord[ii], temp2) // temp = x_i - c_{i-1}*h
		mul_d(&Jv->entry[offset][ii], BED->two, temp); // Jv = 2*(x_i - c_{i-1}*h)
		mul_d(&Jv->entry[offset][ii], &Jv->entry[offset][ii], one_minus_s); // Jv = (1-t)*2*(x_i - c_{i-1}*h)
																			// multiply these entries by (1-t)
		
		mul_d(temp2, &BED->center->coord[ii-1], temp);  // temp2 = c_{i-1} * ( x_i - c_{i-1} * h )
		add_d(&Jv->entry[offset][0], &Jv->entry[offset][0], temp2); // Jv[0] += c_{i-1} * ( x_i - c_{i-1} * h )
	}
	
	
	
	// the homogenizing var deriv
	mul_d(temp, &current_variable_values->coord[0], BED->radius);
	mul_d(temp, temp, BED->radius);  // temp = r^2 h
	
	add_d(&Jv->entry[offset][0], &Jv->entry[offset][0], temp); // Jv[0] = \sum_{i=1}^n {c_{i-1} * ( x_i - c_{i-1} * h )} + r^2 h
	neg_d(&Jv->entry[offset][0], &Jv->entry[offset][0]); // Jv[0] = -Jv[0]
	mul_d(&Jv->entry[offset][0], &Jv->entry[offset][0], BED->two);  // Jv[0] *= 2
	mul_d(&Jv->entry[offset][0], &Jv->entry[offset][0], one_minus_s);  // Jv[0] *= (1-t)
	
	// f = \sum{ ( x_i - c_{i-1} * h )^2 } - r^2 h^2
	//Jv[0] = -2(1-t) (  \sum_{i=1}^n {  c_{i-1} * ( x_i - c_{i-1} * h )  } + r^2 h )
	
	
	
	
	
	// a hardcoded product rule for the two linears.
	for (int ii=0; ii<BED->num_variables; ii++) {
		
		dot_product_d(temp, BED->starting_linear[0], current_variable_values);
		mul_d(temp, temp, &BED->starting_linear[1]->coord[ii]);
		
		dot_product_d(temp2, BED->starting_linear[1], current_variable_values);
		mul_d(temp2, temp2, &BED->starting_linear[0]->coord[ii]);
		
		add_d(temp, temp, temp2);
		mul_d(temp2, temp, gamma_s);
		
		//temp2 = gamma s * (L_1(x) * L_0[ii] + L_0(x) * L_1[ii])
		
		//temp2 now has the value of the derivative of the start system wrt x_i
		
		add_d(&Jv->entry[offset][ii], &Jv->entry[offset][ii], temp2);
	}
	
	
	
	
	//// // / / /// // //     finally, the Jp entry for sphere equation's homotopy.
	//Jp = -f_sph + gamma f_start
	neg_d(&Jp->entry[offset][0], func_val_sphere);
	mul_d(temp, BED->gamma, func_val_start);
	add_d(&Jp->entry[offset][0], &Jp->entry[offset][0], temp);
	
	
	
	
	
	
	
	
	
	
	//////////////
	//
	// function values for the static linears
	//
	////////////////////
	
	offset++;
	for (mm=0; mm<BED->num_static_linears; ++mm) {
		dot_product_d(&funcVals->coord[mm+offset], BED->static_linear[mm], current_variable_values);
	}
	
	for (mm=0; mm<BED->num_static_linears; ++mm) {
		for (ii=0; ii<BED->num_variables; ii++) {
			set_d(&Jv->entry[mm+offset][ii], &BED->static_linear[mm]->coord[ii]);
		}
	}
	
	//Jp is 0 for the static linears
	
	
	
	//////////////
	//
	// the entries for the patch equations.
	//
	////////////////////
	if (offset+BED->num_static_linears != BED->num_variables-BED->patch.num_patches) {
		std::cout << color::red() << "mismatch in offset!\nleft: " <<
		offset+BED->num_static_linears << " right " << BED->num_variables-BED->patch.num_patches << color::console_default() << std::endl;
		mypause();
	}
	
	offset = BED->num_variables-BED->patch.num_patches;
	for (ii=0; ii<BED->patch.num_patches; ii++)
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	
	
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{  // Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	//Jp is 0 for the patch.
	
	
	
	
	
	// done!  yay!
	
	if (BED->verbose_level()==16 || BED->verbose_level()==-16) {
		//uncomment to see screen output of important variables at each solve step.
		print_point_to_screen_matlab(BED->center,"center");
		print_comp_matlab(BED->radius,"radius");
		
		printf("gamma = %lf+1i*%lf;\n", BED->gamma->r, BED->gamma->i);
		printf("time = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		print_point_to_screen_matlab(current_variable_values,"currvars");
		print_point_to_screen_matlab(funcVals,"F");
		print_matrix_to_screen_matlab(Jv,"Jv");
		print_matrix_to_screen_matlab(Jp,"Jp");
		
	}
	
	
	BED->SLP_memory.set_globals_null();
	
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);
	clear_vec_d(AtimesF);
	
	
	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	clear_mat_d(Jv_Patch);
	clear_mat_d(AtimesJ);
	
#ifdef printpathsphere
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,current_variable_values);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	
	return 0;
}




//this derived from basic_eval_d
int sphere_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	//	print_comp_mp_matlab(pathVars,"pathvars");
	
	
	sphere_eval_data_mp *BED = (sphere_eval_data_mp *)ED; // to avoid having to cast every time
	
	BED->SLP_memory.set_globals_to_this();
	
	int ii, jj, mm; // counters
	int offset;
	comp_mp one_minus_s, gamma_s;
	comp_mp temp, temp2;
	comp_mp func_val_sphere, func_val_start;
	init_mp(one_minus_s);  init_mp(gamma_s);
	init_mp(temp);  init_mp(temp2);
	init_mp(func_val_start);
	init_mp(func_val_sphere);
	
	set_one_mp(one_minus_s);
	sub_mp(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
	mul_mp(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	vec_mp patchValues; init_vec_mp(patchValues, 0);
	vec_mp temp_function_values; init_vec_mp(temp_function_values,0);
	vec_mp AtimesF; init_vec_mp(AtimesF,BED->randomizer()->num_rand_funcs()); AtimesF->size = BED->randomizer()->num_rand_funcs();// declare  // initialize
	
	
	
	mat_mp temp_jacobian_functions; init_mat_mp(temp_jacobian_functions,BED->randomizer()->num_base_funcs(),BED->num_variables);
	temp_jacobian_functions->rows = BED->randomizer()->num_base_funcs(); temp_jacobian_functions->cols = BED->num_variables;
	mat_mp temp_jacobian_parameters; init_mat_mp(temp_jacobian_parameters,0,0);
	mat_mp Jv_Patch; init_mat_mp(Jv_Patch, 0, 0);
	mat_mp AtimesJ; init_mat_mp(AtimesJ,BED->randomizer()->num_rand_funcs(),BED->num_variables);
	AtimesJ->rows = BED->randomizer()->num_rand_funcs(); AtimesJ->cols = BED->num_variables;
	
	
	//set the sizes
	change_size_vec_mp(funcVals,BED->num_variables); funcVals->size = BED->num_variables;
	change_size_mat_mp(Jv, BED->num_variables, BED->num_variables); Jv->rows = Jv->cols = BED->num_variables; //  -> this should be square!!!
	
	for (ii=0; ii<BED->num_variables; ii++)
		for (jj=0; jj<BED->num_variables; jj++)
			set_zero_mp(&Jv->entry[ii][jj]);
	
	
	
	// evaluate the SLP to get the system's whatnot.
	evalProg_mp(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	
	// evaluate the patch
	patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	// note that you can only really do this AFTER you are done calling other evaluators.
	// set parVals & parDer correctly
	
	// i.e. these must remain here, or below.  \/
	change_size_point_mp(parVals, 1);
	change_size_vec_mp(parDer, 1);
	change_size_mat_mp(Jp, BED->num_variables, 1); Jp->rows = BED->num_variables; Jp->cols = 1;
	for (ii=0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);
	
	
	parVals->size = parDer->size = 1;
	set_mp(&parVals->coord[0], pathVars); // s = t
	set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	///////////////////////////
	//
	// the original (randomized) functions.
	//
	///////////////////////////////////
	
	BED->randomizer()->randomize(AtimesF,AtimesJ,temp_function_values,temp_jacobian_functions,&current_variable_values->coord[0]);
	
	
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real orthogonal) randomization
		set_mp(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	for (ii = 0; ii < BED->randomizer()->num_rand_funcs(); ii++)
		for (jj = 0; jj < BED->num_variables; jj++)
			set_mp(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	//Jp is 0 for the equations.
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	///////////////////
	//
	//  the sphere equation.
	//
	//////////////////////////
	
	offset = BED->randomizer()->num_rand_funcs();
	
	mul_mp(func_val_sphere, BED->radius, BED->radius);
	neg_mp(func_val_sphere, func_val_sphere);
	mul_mp(func_val_sphere, func_val_sphere, &current_variable_values->coord[0]);
	mul_mp(func_val_sphere, func_val_sphere, &current_variable_values->coord[0]);
	//f_sph = -r^2*h^2
	
	
	
	for (int ii=1; ii<BED->num_natural_vars; ii++) {
		mul_mp(temp2, &BED->center->coord[ii-1], &current_variable_values->coord[0]); // temp2 = c_{i-1}*h
		
		sub_mp(temp, &current_variable_values->coord[ii], temp2);  // temp = x_i - h*c_{i-1}
		mul_mp(temp2, temp, temp);                                 // temp2 = (x_i - h*c_{i-1})^2
		add_mp(func_val_sphere, func_val_sphere, temp2);           // f_sph += (x_i - h*c_{i-1})^2
	}
	
	
	
	set_one_mp(func_val_start);
	for (mm=0; mm<2; ++mm) {
		dot_product_mp(temp, BED->starting_linear[mm], current_variable_values);
		mul_mp(func_val_start, func_val_start, temp);
		//f_start *= L_i (x)
	}
	
	
	// combine the function values
	mul_mp(temp, one_minus_s, func_val_sphere);
	mul_mp(temp2, gamma_s, func_val_start);
	add_mp(&funcVals->coord[offset], temp, temp2);
	// f = (1-t) f_sph + gamma t f_start
	
	
	
	
	//// / / / / / /    now the derivatives wrt x
	
	//  first we store the derivatives of the target function, the sphere.  then we will add the part for the linear product start.
	
	
	
	//ddx for sphere
	
	
	
	for (int ii=1; ii<BED->num_natural_vars; ii++) {
		mul_mp(temp2, &BED->center->coord[ii-1], &current_variable_values->coord[0]); // temp2 = c_{i-1}*h
		sub_mp(temp, &current_variable_values->coord[ii], temp2) // temp = x_i - c_{i-1}*h
		mul_mp(&Jv->entry[offset][ii], BED->two, temp); // Jv = 2*(x_i - c_{i-1}*h)
		mul_mp(&Jv->entry[offset][ii], &Jv->entry[offset][ii], one_minus_s); // Jv = (1-t)*2*(x_i - c_{i-1}*h)
		
		
		mul_mp(temp2, &BED->center->coord[ii-1], temp);  // temp2 = c_{i-1} * ( x_i - c_{i-1} * h )
		add_mp(&Jv->entry[offset][0], &Jv->entry[offset][0], temp2); // Jv[0] += c_{i-1} * ( x_i - c_{i-1} * h )
	}
	// multiply these entries by (1-t)
	
	
	// the homogenizing var deriv
	mul_mp(temp, &current_variable_values->coord[0], BED->radius);
	mul_mp(temp, temp, BED->radius);  // temp = r^2 h
	
	add_mp(&Jv->entry[offset][0], &Jv->entry[offset][0], temp); // Jv[0] = \sum_{i=1}^n {c_{i-1} * ( x_i - c_{i-1} * h )} + r^2 h
	neg_mp(&Jv->entry[offset][0], &Jv->entry[offset][0]); // Jv[0] = -Jv[0]
	mul_mp(&Jv->entry[offset][0], &Jv->entry[offset][0], BED->two);  // Jv[0] *= 2
	mul_mp(&Jv->entry[offset][0], &Jv->entry[offset][0], one_minus_s);  // Jv[0] *= (1-t)
	
	// f = \sum{ ( x_i - c_{i-1} * h )^2 } - r^2 h^2
	//Jv = -2(1-t) (  \sum_{i=1}^n {  c_{i-1} * ( x_i - c_{i-1} * h )  } + r^2 h )
	
	
	// a hardcoded product rule for the two linears.
	for (int ii=0; ii<BED->num_variables; ii++) {
		
		dot_product_mp(temp, BED->starting_linear[0], current_variable_values);
		mul_mp(temp, temp, &BED->starting_linear[1]->coord[ii]);
		
		dot_product_mp(temp2, BED->starting_linear[1], current_variable_values);
		mul_mp(temp2, temp2, &BED->starting_linear[0]->coord[ii]);
		
		add_mp(temp, temp, temp2);
		mul_mp(temp2, temp, gamma_s);
		
		//temp2 = gamma s * (L_1(x) * L_0[ii] + L_0(x) * L_1[ii])
		
		//temp2 now has the value of the derivative of the start system wrt x_i
		
		add_mp(&Jv->entry[offset][ii], &Jv->entry[offset][ii], temp2);
	}
	
	
	
	// finally, the Jp entry for sphere equation's homotopy.
	//Jp = -f_sph + gamma f_start
	neg_mp(&Jp->entry[offset][0], func_val_sphere);
	mul_mp(temp, BED->gamma, func_val_start);
	add_mp(&Jp->entry[offset][0], &Jp->entry[offset][0], temp);
	
	
	
	
	
	
	
	//////////////
	//
	// function values for the static linears
	//
	////////////////////
	
	offset++;
	for (mm=0; mm<BED->num_static_linears; ++mm) {
		dot_product_mp(&funcVals->coord[mm+offset], BED->static_linear[mm], current_variable_values);
	}
	
	for (mm=0; mm<BED->num_static_linears; ++mm) {
		for (ii=0; ii<BED->num_variables; ii++) {
			set_mp(&Jv->entry[offset+mm][ii], &BED->static_linear[mm]->coord[ii]);
		}
	}
	
	//Jp is 0 for the static linears
	
	
	
	//////////////
	//
	// the entries for the patch equations.
	//
	////////////////////
	if (offset+BED->num_static_linears != BED->num_variables-BED->patch.num_patches) {
		std::cout << color::red() << "mismatch in offset!\nleft: " <<
		offset+BED->num_static_linears << " right " << BED->num_variables-BED->patch.num_patches << color::console_default() << std::endl;
		mypause();
	}
	
	offset = BED->num_variables-BED->patch.num_patches;
	for (ii=0; ii<BED->patch.num_patches; ii++)
		set_mp(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	
	
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{  // Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	//Jp is 0 for the patch.
	
	
	
	
	
	// done!  yay!
	
	if (BED->verbose_level()==16 || BED->verbose_level()==-16) {
		//uncomment to see screen output of important variables at each solve step.
		
		print_comp_matlab(pathVars, "t_mp");
		print_comp_matlab(BED->gamma, "gamma_mp");
		print_point_to_screen_matlab(current_variable_values,"currvars_mp");
		print_point_to_screen_matlab(funcVals,"F_mp");
		print_matrix_to_screen_matlab(Jv,"Jv_mp");
		print_matrix_to_screen_matlab(Jp,"Jp_mp");
		
		
	}
	
	
	BED->SLP_memory.set_globals_null();
	
	clear_mp(temp);
	clear_mp(temp2);
	clear_mp(gamma_s);
	clear_mp(one_minus_s);
	clear_mp(func_val_sphere);
	clear_mp(func_val_start);
	
	clear_vec_mp(patchValues);
	clear_vec_mp(temp_function_values);
	clear_vec_mp(AtimesF);
	
	
	clear_mat_mp(temp_jacobian_functions);
	clear_mat_mp(temp_jacobian_parameters);
	clear_mat_mp(Jv_Patch);
	clear_mat_mp(AtimesJ);
	
	
	
	
	return 0;
}





int sphere_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
{
	sphere_eval_data_d *BED_d = NULL;
	sphere_eval_data_mp *BED_mp = NULL;
	
	*out_prec = in_prec;
	
	
	
	if (in_prec < 64)
	{ // compute out_d
		sphere_eval_data_d *BED_d = (sphere_eval_data_d *)ED_d;
		
		comp_d denom;
		change_size_vec_d(out_d,in_d->size-1);
		out_d->size = in_d->size-1;
		
		set_d(denom, &in_d->coord[0]);
		
		for (int ii=0; ii<BED_d->num_variables-1; ++ii) {
			set_d(&out_d->coord[ii],&in_d->coord[ii+1]);
			div_d(&out_d->coord[ii],&out_d->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		
		
		//		print_point_to_screen_matlab(in_d,"in");
		//		print_point_to_screen_matlab(out_d,"out");
		
		
	}
	else
	{ // compute out_mp
		sphere_eval_data_mp *BED_mp = (sphere_eval_data_mp *)ED_mp;
		
		setprec_point_mp(out_mp, *out_prec);
		
		comp_mp denom; init_mp(denom);
		change_size_vec_mp(out_mp,in_mp->size-1);
		out_mp->size = in_mp->size-1;
		
		set_mp(denom, &in_mp->coord[0]);
		
		for (int ii=0; ii<BED_mp->num_variables-1; ++ii) {
			set_mp(&out_mp->coord[ii],&in_mp->coord[ii+1]);
			div_mp(&out_mp->coord[ii],&out_mp->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		
		clear_mp(denom);
		
		
		// set prec on out_mp
		
		
		//		print_point_to_screen_matlab(in_mp,"in");
		//		print_point_to_screen_matlab(out_mp,"out");
		
	}
	
	
	BED_d = NULL;
	BED_mp = NULL;
	
	
	
	
	return 0;
}




int change_sphere_eval_prec(void const *ED, int new_prec)
{
	
	sphere_eval_data_mp *BED = (sphere_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	
	BED->SLP->precision = new_prec;
	// change the precision for the patch
	changePatchPrec_mp(new_prec, &BED->patch);
	
	
	if (new_prec != BED->curr_prec){
		
		if (BED->verbose_level() >=8){
			std::cout << color::brown();
			printf("prec  %ld\t-->\t%d\n",BED->curr_prec, new_prec);
			std::cout << color::console_default();
		}
		
		BED->curr_prec = new_prec;
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		for (int ii=0; ii<1; ++ii){
			change_prec_point_mp(BED->starting_linear[ii],new_prec);
			vec_cp_mp(BED->starting_linear[ii], BED->starting_linear_full_prec[ii]);
		}
		
		for (int ii=0; ii<BED->num_static_linears; ++ii){
			change_prec_point_mp(BED->static_linear[ii],new_prec);
			vec_cp_mp(BED->static_linear[ii], BED->static_linear_full_prec[ii]);
		}
		
		change_prec_point_mp(BED->center, new_prec);
		vec_cp_mp(BED->center, BED->center_full_prec);
		
		setprec_mp(BED->radius, new_prec);
		set_mp(BED->radius, BED->radius_full_prec);
		
		setprec_mp(BED->two, new_prec);
		set_mp(BED->two, BED->two_full_prec);
		
		BED->randomizer()->change_prec(new_prec);
		
	}
	
	
	return 0;
}





int check_issoln_sphere_d(endgame_data_t *EG,
						  tracker_config_t *T,
						  void const *ED)
{
	sphere_eval_data_d *BED = (sphere_eval_data_d *)ED; // to avoid having to cast every time
	
	BED->SLP_memory.set_globals_to_this();
	int ii;
	
	
	
	
	double n1, n2, max_rat;
	point_d f;
	eval_struct_d e;
	
	
	init_point_d(f, 1);
	init_eval_struct_d(e,0, 0, 0);
	
	max_rat = MAX(T->ratioTol,0.99999);
	
	// setup threshold based on given threshold and precision
	//	if (num_digits > 300)
	//		num_digits = 300;
	//	num_digits -= 2;
	double tol = MAX(T->funcResTol, 1e-8);
	
	
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
		vec_mp_to_d(prev_pt,EG->PD_mp.point);
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, prev_pt, EG->PD_d.time, BED->SLP);
		clear_vec_d(prev_pt);}
	else{
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_d, EG->PD_d.time, BED->SLP);
	}
	
	// compare the function values
	int isSoln = 1;
	for (ii = 0; (ii < BED->SLP->numFuncs) && isSoln; ii++)
	{
		n1 = d_abs_d( &e.funcVals->coord[ii]); // corresponds to final point
		n2 = d_abs_d( &f->coord[ii]); // corresponds to the previous point
		
		if (tol <= n1 && n1 <= n2)
		{ // compare ratio
			if (n1 > max_rat * n2){
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
	
	BED->SLP_memory.set_globals_null();
	
	
	return isSoln;
	
}


int check_issoln_sphere_mp(endgame_data_t *EG,
						   tracker_config_t *T,
						   void const *ED)
{
	sphere_eval_data_mp *BED = (sphere_eval_data_mp *)ED; // to avoid having to cast every time
	BED->SLP_memory.set_globals_to_this();
	int ii;
	
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
	
	mpf_set_d(max_rat, MAX(T->ratioTol,0.99999));
	
	point_mp f; init_point_mp(f, 1);f->size = 1;
	eval_struct_mp e; init_eval_struct_mp(e, 0, 0, 0);
	
	
	int num_digits = prec_to_digits((int) mpf_get_default_prec());
	// setup threshold based on given threshold and precision
	if (num_digits > 300)
		num_digits = 300;
	num_digits -= 4;
	double my_guarantor = MAX(T->funcResTol, 1e-8);
	double tol = MAX(my_guarantor, pow(10,-num_digits));
	mpf_set_d(zero_thresh, tol);
	
	//this one guaranteed by entry condition
	//	sphere_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, ED);
	evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, BED->SLP);
	//	print_point_to_screen_matlab(e.funcVals,"howfaroff");
	
	if (EG->last_approx_prec < 64)
	{ // copy to _mp
		vec_mp temp_vec;  init_vec_mp(temp_vec,0);
		point_d_to_mp(temp_vec, EG->last_approx_d);
		evalProg_mp(f, e.parVals, e.parDer, e.Jv, e.Jp, temp_vec, EG->PD_mp.time, BED->SLP);
		clear_vec_mp(temp_vec);
	}
	else
	{
		evalProg_mp(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_mp, EG->PD_mp.time, BED->SLP);
	}
	
	// compare the function values
	int isSoln = 1;
	for (ii = 0; ii < BED->SLP->numFuncs && isSoln; ii++)
	{
		mpf_abs_mp(n1, &e.funcVals->coord[ii]);
		mpf_abs_mp(n2, &f->coord[ii]);
		
		if ( (mpf_cmp(zero_thresh, n1) <= 0) &&  (mpf_cmp(n1, n2) <= 0) )
		{ // compare ratio
			mpf_mul(n2, max_rat, n2);
			if (mpf_cmp(n1, n2) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 1\n");
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
	
	
	
	if (!isSoln) {
		
		print_point_to_screen_matlab(e.funcVals,"terminal_func_vals");
		printf("tol was %le\nmax_rat was ",tol);
		mpf_out_str(NULL,10,15,max_rat);
		std::cout << "\n";
	}
	
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	clear_vec_mp(f);
	
	BED->SLP_memory.set_globals_null();
	return isSoln;
	
}



















