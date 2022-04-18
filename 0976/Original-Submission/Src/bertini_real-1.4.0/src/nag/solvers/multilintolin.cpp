#include "nag/solvers/multilintolin.hpp"




void multilintolin_eval_data_mp::init()
{
	this->is_solution_checker_d = &check_issoln_multilintolin_d;
	this->is_solution_checker_mp = &check_issoln_multilintolin_mp;
	this->evaluator_function_d = &multilin_to_lin_eval_d;
	this->evaluator_function_mp = &multilin_to_lin_eval_mp;
	this->precision_changer = &change_multilintolin_eval_prec;
	this->dehomogenizer = &multilintolin_dehom;
	
	this->num_linears = 0;
	old_linear_full_prec = old_linear = NULL;
	current_linear_full_prec = current_linear = NULL;
	
}





int multilintolin_eval_data_mp::setup(MultilinConfiguration & config,
									  const WitnessSet & W,
									  vec_mp * target_linears,
									  SolverConfiguration & solve_options)
{
	
	
	if (!config.randomizer()->is_ready()) {
		std::cout << "don't have multilin randomizer!" << std::endl;
		deliberate_segfault();
	}
	
	if (!config.have_mem) {
		std::cout << "don't have multilin SLP memory!" << std::endl;
		deliberate_segfault();
	}
	
	this->SLP_memory = config.SLP_memory;
	
	
	num_variables = W.num_variables();
	
	
	// set up the vectors to hold the two linears.
	if (this->num_linears==0) {
		current_linear = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
		old_linear = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
	}
	else
	{
		current_linear = (vec_mp *) br_realloc(current_linear, W.num_linears()*sizeof(vec_mp));
		old_linear = (vec_mp *) br_realloc(old_linear, W.num_linears()*sizeof(vec_mp));
	}
	
	
	//actually do the transfer
	for (unsigned int ii=0; ii<W.num_linears(); ii++) {
		init_vec_mp(current_linear[ii],0);
		init_vec_mp(old_linear[ii],0);
		vec_cp_mp(current_linear[ii],target_linears[ii]);
		vec_cp_mp(old_linear[ii],W.linear(ii));
	}
	
	
	if (this->MPType==2) {
		if (this->num_linears==0) {
			current_linear_full_prec = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
			old_linear_full_prec = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp));
		}
		else
		{
			current_linear_full_prec = (vec_mp *) br_realloc(current_linear_full_prec, W.num_linears()*sizeof(vec_mp));
			old_linear_full_prec = (vec_mp *) br_realloc(old_linear_full_prec, W.num_linears()*sizeof(vec_mp));
		}
		
		for (unsigned int ii=0; ii<W.num_linears(); ii++) {
			init_vec_mp2(current_linear_full_prec[ii],0,1024);
			init_vec_mp2(old_linear_full_prec[ii],0,1024);
			
			vec_cp_mp(current_linear_full_prec[ii],target_linears[ii]);
			vec_cp_mp(old_linear_full_prec[ii],W.linear(ii));
		}
	}
	
	this->num_linears= W.num_linears();
	
	
	
	
	
	
	
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









int multilintolin_eval_data_mp::send(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "multilintolin_eval_data_mp::send()" << std::endl;
#endif
	
	int solver_choice = MULTILIN;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	SolverMultiplePrecision::send(mpi_config);
	
	
	int *buffer = new int[1];
	
	
	buffer[0] = num_linears;
	
	// now can actually send the data.
	
	MPI_Bcast(buffer,1,MPI_INT, mpi_config.head(), mpi_config.comm());
	
	delete[] buffer;
	
	if (this->MPType==2){
		for (int ii=0; ii<num_linears; ii++) {
			bcast_vec_mp(old_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
			bcast_vec_mp(current_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else {
		for (int ii=0; ii<num_linears; ii++) {
			bcast_vec_mp(old_linear[ii], mpi_config.id(), mpi_config.head());
			bcast_vec_mp(current_linear[ii], mpi_config.id(), mpi_config.head());
		}
		
	}
	
	
	return SUCCESSFUL;
}




int multilintolin_eval_data_mp::receive(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "multilintolin_eval_data_mp::receive()" << std::endl;
#endif
	
	
	int *buffer = new int[1];
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	
	if (buffer[0] != MULTILIN) {
		std::cout << "worker failed to confirm it is receiving the multilin type eval data" << std::endl;
		mpi_config.abort(777);
	}
	
	SolverMultiplePrecision::receive(mpi_config);
	
	
	// now can actually receive the data from whoever.
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	num_linears = buffer[0];
	delete[] buffer;
	

	old_linear     = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
	current_linear = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
	
	if (this->MPType==2) {
		old_linear_full_prec = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
		current_linear_full_prec = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
		
		for (int ii=0; ii<num_linears; ii++) {
			init_vec_mp(old_linear[ii],1);
			init_vec_mp(current_linear[ii],1);
			
			init_vec_mp2(old_linear_full_prec[ii],1,1024);
			init_vec_mp2(current_linear_full_prec[ii],1,1024);
			
			bcast_vec_mp(old_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
			bcast_vec_mp(current_linear_full_prec[ii], mpi_config.id(), mpi_config.head());
			
			vec_cp_mp(old_linear[ii],old_linear_full_prec[ii]);
			vec_cp_mp(current_linear[ii],current_linear_full_prec[ii]);
		}
		
		
	}
	else{ // MPType == 1
		for (int ii=0; ii<num_linears; ii++) {
			init_vec_mp(old_linear[ii],1);
			init_vec_mp(current_linear[ii],1);
			bcast_vec_mp(old_linear[ii], mpi_config.id(), mpi_config.head());
			bcast_vec_mp(current_linear[ii], mpi_config.id(), mpi_config.head());
		}
	}
	

	return SUCCESSFUL;
}


















void multilintolin_eval_data_d::init()
{
	
	if (this->MPType==2){
		this->BED_mp = new multilintolin_eval_data_mp(2);
		SolverDoublePrecision::BED_mp = this->BED_mp;                   //   <---------  you gotta do this cuz of masking problems.
	}
	else{
		this->BED_mp = NULL;
	}
	
	
	
	this->is_solution_checker_d = &check_issoln_multilintolin_d;
	this->is_solution_checker_mp = &check_issoln_multilintolin_mp;
	this->evaluator_function_d = &multilin_to_lin_eval_d;
	this->evaluator_function_mp = &multilin_to_lin_eval_mp;
	this->precision_changer = &change_multilintolin_eval_prec;
	this->dehomogenizer = &multilintolin_dehom;
	
	this->num_linears = 0;
	old_linear = NULL;
	current_linear = NULL;
	
}







int multilintolin_eval_data_d::setup(MultilinConfiguration & config,
									 const WitnessSet & W,
									 vec_mp * target_linears,
									 SolverConfiguration & solve_options)
{
	
	SLP_memory = config.SLP_memory;
	// set up the vectors to hold the two linears.
	if (this->num_linears==0) {
		current_linear = (vec_d *) br_malloc(W.num_linears()*sizeof(vec_d));
		old_linear = (vec_d *) br_malloc(W.num_linears()*sizeof(vec_d));
	}
	else
	{
		current_linear = (vec_d *) br_realloc(current_linear, W.num_linears()*sizeof(vec_d));
		old_linear = (vec_d *) br_realloc(old_linear, W.num_linears()*sizeof(vec_d));
	}
	
	for (unsigned int ii=0; ii<W.num_linears(); ii++) {
		init_vec_d(current_linear[ii],0);
		init_vec_d(old_linear[ii],0);
		vec_mp_to_d(current_linear[ii],target_linears[ii]);
		vec_mp_to_d(old_linear[ii],W.linear(ii));
	}
	
	num_linears = W.num_linears();
	
	num_variables = W.num_variables();
	
	
	verbose_level(solve_options.verbose_level());
	
	
	
	generic_setup_patch(&patch,W);
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_d(this->gamma); // set gamma to be random complex value
	else
		set_one_d(this->gamma);
	
	
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(config, W, target_linears, solve_options);
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
	}
	
	
	SolverDoublePrecision::setup(config.SLP, config.randomizer());
	return 0;
}













int multilintolin_eval_data_d::send(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
    std::cout << "multilintolin_eval_data_d::send()" << std::endl;
#endif
    int solver_choice = MULTILIN;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	// send the confirmation integer, to ensure that we are sending the correct type.
    
    if (this->MPType==2) {
		this->BED_mp->send(mpi_config);
	}
    
    
	//send the base class stuff.
	SolverDoublePrecision::send(mpi_config);
	
	// now can actually send the data.
	int *buffer = new int[1];
	buffer[0] = num_linears;
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head() , mpi_config.comm());
	delete[] buffer;
	
	for (int ii=0; ii<num_linears; ii++) {
		bcast_vec_d(old_linear[ii], mpi_config.id(), mpi_config.head());
		bcast_vec_d(current_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	return SUCCESSFUL;
}

int multilintolin_eval_data_d::receive(ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "multilintolin_eval_data_d::receive()" << std::endl;
#endif
	
    int *buffer = new int[1];
	
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	if (buffer[0] != MULTILIN){
		std::cout << "worker failed to confirm it is receiving the multilin type eval data" << std::endl;
		mpi_config.abort(777);
	}
    
    
    
	if (this->MPType==2) {
		this->BED_mp->receive(mpi_config);
	}
    
    

	SolverDoublePrecision::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	
	
	
	MPI_Bcast(buffer, 1, MPI_INT, mpi_config.head(), mpi_config.comm());
	num_linears = buffer[0];
	delete[] buffer;
	

	
	old_linear     = (vec_d *) br_malloc(num_linears*sizeof(vec_d));
	current_linear = (vec_d *) br_malloc(num_linears*sizeof(vec_d));
	
	for (int ii=0; ii<num_linears; ii++) {
		init_vec_d(old_linear[ii],0);
		init_vec_d(current_linear[ii],0);
		bcast_vec_d(old_linear[ii], mpi_config.id(), mpi_config.head());
		bcast_vec_d(current_linear[ii], mpi_config.id(), mpi_config.head());
	}
	
	return SUCCESSFUL;
}









int multilin_solver_master_entry_point(const WitnessSet & W, // carries with it the start points, and the linears.
									   SolverOutput & solve_out, // new data goes in here
									   vec_mp * target_linears,
									   MultilinConfiguration &		config,
									   SolverConfiguration		& solve_options)
{
	
	bool prev_parallel_state = solve_options.force_no_parallel();
	
	if ( int(W.num_points()) < solve_options.num_procs()-1) {
		solve_options.force_no_parallel(true);
	}
	
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(MULTILIN);
	}
	
	multilintolin_eval_data_d *ED_d = NULL;
	multilintolin_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new multilintolin_eval_data_d(0);
			
			ED_d->setup(config,
						W,
						target_linears,
						solve_options);
			break;
			
		case 1:
			ED_mp = new multilintolin_eval_data_mp(1);
			
			ED_mp->setup(config,
						 W,
						 target_linears,
						 solve_options);
			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //    <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new multilintolin_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(config,
						W,
						target_linears,
						solve_options);
			
			
			
			adjust_tracker_AMP(& (solve_options.T), W.num_variables());
			// initialize latest_newton_residual_mp
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
			std::cout << "MPTYPE is not 0, 1, 2, cannot clear... " << std::endl;
			br_exit(398);
			break;
	}
	
	
	for (unsigned int jj=0; jj<W.num_linears(); jj++)
	{
		solve_out.add_linear(target_linears[jj]);
	}
	
	solve_options.force_no_parallel(prev_parallel_state);

	
	return SUCCESSFUL;
	
}








int multilin_slave_entry_point(SolverConfiguration & solve_options)
{
	
	
	// already received the flag which indicated that this worker is going to be performing the nullspace calculation.
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	int *settings_buffer = (int *) br_malloc(2*sizeof(int));
	MPI_Bcast(settings_buffer,2,MPI_INT, 0,solve_options.comm());
	solve_options.robust = settings_buffer[0];
	solve_options.use_gamma_trick = settings_buffer[1];
	free(settings_buffer);
	
	multilintolin_eval_data_d *ED_d = NULL;
	multilintolin_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new multilintolin_eval_data_d(0);
			ED_d->receive(solve_options);
			break;
			
		case 1:
			ED_mp = new multilintolin_eval_data_mp(1);
			ED_mp->receive(solve_options);
			
			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //    <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new multilintolin_eval_data_d(2);
			ED_mp = ED_d->BED_mp;
			ED_d->receive(solve_options);
			
			
			
			
			// initialize latest_newton_residual_mp
//			mpf_init(solve_options.T.latest_newton_residual_mp);   //    <------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		default:
			break;
	}
	
    
	
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "multilin_output",
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








int multilin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, build for bertini_real
	
	// uncomment to see the time at each step.
	//	printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
	
	multilintolin_eval_data_d *BED = (multilintolin_eval_data_d *)ED; // to avoid having to cast every time
	
	
	
	int offset;
	comp_d one_minus_s, gamma_s;
	comp_d temp, temp2;
	
	
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
	
	for (int ii=0; ii<BED->num_variables; ii++) {
		for (int jj=0; jj<BED->num_variables; jj++) {
			set_zero_d(&Jv->entry[ii][jj]);
		}
	}
	
	BED->SLP_memory.set_globals_to_this();
	// evaluate the SLP to get the system's whatnot.
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	BED->SLP_memory.set_globals_null();
	
	// evaluate the patch
	patch_eval_d(patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	// note that you can only really do this AFTER you are done calling other evaluators.
	// set parVals & parDer correctly
	
	// i.e. these MUST remain here, or below.  \/
	change_size_point_d(parVals, 1);
	change_size_vec_d(parDer, 1);
	change_size_mat_d(Jp, BED->num_variables, 1); Jp->rows = BED->num_variables; Jp->cols = 1;
	for (int ii=0; ii<BED->num_variables; ii++)
		set_zero_d(&Jp->entry[ii][0]);
	
	
	parVals->size = parDer->size = 1;
	set_d(&parVals->coord[0], pathVars); // s = t
	set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	///////// / / / /  /   /
	// combine everything
	///////// / / / /  /   /
	
	// randomize
	BED->randomizer()->randomize(AtimesF,AtimesJ,temp_function_values,temp_jacobian_functions,&current_variable_values->coord[0]);
	
	for (int ii=0; ii<AtimesF->size; ii++)  // for each function, after (real orthogonal) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	
	
	
	
	
	
	//////////////
	//
	// function values and Jp for the linears
	//
	////////////////////
	
	
	offset = BED->num_variables-BED->patch.num_patches-BED->num_linears;
	for (int ii=0; ii<BED->num_linears; ++ii) {
		dot_product_d(temp,BED->current_linear[ii],current_variable_values);
		
		mul_d(&funcVals->coord[ii+offset],one_minus_s,temp);
		neg_d(&Jp->entry[offset+ii][0],temp);
		
		dot_product_d(temp2,BED->old_linear[ii],current_variable_values);
		
		mul_d(temp,temp2,BED->gamma);
		add_d(&Jp->entry[offset+ii][0],&Jp->entry[offset+ii][0],temp); // Jp = -curr_lin_val + gamma*old_lin_val
		
		mul_d(temp2,temp2,gamma_s);
		add_d(&funcVals->coord[ii+offset],&funcVals->coord[ii+offset],temp2); // f = (1-s)*curr_lin_val + gamma*s*old_lin_val
	}
	
	
	//////////////
	//
	//set the function PATCH values
	//
	////////////////////
	
	
	
	offset = BED->num_variables-BED->patch.num_patches;
	for (int ii=0; ii<BED->patch.num_patches; ii++)
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//////////////
	//
	// set the JACOBIAN values.
	//
	////////////////////
	
	//first, the entries related to the functions
	
	for (int ii = 0; ii < BED->randomizer()->num_rand_funcs(); ii++)
		for (int jj = 0; jj < BED->num_variables; jj++)
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	
	
	
	//////////////
	//
	// Entries for the linears we are homotoping
	//
	////////////////////
	
	
	offset = BED->num_variables-BED->patch.num_patches-BED->num_linears;
	for (int ii=0; ii<BED->num_linears; ++ii) {
		for (int jj=0; jj<BED->num_variables; jj++) {
			mul_d(temp,&BED->old_linear[ii]->coord[jj],gamma_s);
			mul_d(temp2,&BED->current_linear[ii]->coord[jj],one_minus_s);
			add_d(&Jv->entry[offset+ii][jj], temp,temp2);
		}
	}
	
	
	

	
	
	//////////////
	//
	// the entries in the jacobian for the patch equations.
	//
	////////////////////
	
	offset = BED->num_variables - BED->patch.num_patches;
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{  // Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	// done!  yay!
	
	if (BED->verbose_level()==12) {
		//uncomment to see screen output of important variables at each solve step.
		printf("gamma = %lf+1i*%lf;\n", BED->gamma->r, BED->gamma->i);
		printf("time = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		print_matrix_to_screen_matlab( AtimesJ,"R*jac");
		print_point_to_screen_matlab(current_variable_values,"currvars");
		print_point_to_screen_matlab(BED->current_linear[0],"new");
		print_point_to_screen_matlab(BED->old_linear[0],"old");
		print_point_to_screen_matlab(funcVals,"F");
		print_matrix_to_screen_matlab(Jv,"Jv");
		print_matrix_to_screen_matlab(Jp,"Jp");
		std::cout << BED->randomizer() << std::endl;//
		
	}
	
	
	
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);
	clear_vec_d(AtimesF);
	
	
	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	clear_mat_d(Jv_Patch);
	clear_mat_d(AtimesJ);
	
#ifdef printpathmultilintolin
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
int multilin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	//	print_comp_mp_matlab(pathVars,"pathvars");
	
	
	multilintolin_eval_data_mp *BED = (multilintolin_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	
	comp_mp one_minus_s, gamma_s; init_mp(one_minus_s); init_mp(gamma_s);
	comp_mp temp, temp2; init_mp(temp); init_mp(temp2);
	
	int offset;
	
	
	
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
	
	for (int ii=0; ii<BED->num_variables; ii++)
		for (int jj=0; jj<BED->num_variables; jj++)
			set_zero_mp(&Jv->entry[ii][jj]);
	
	
	
	// evaluate the SLP to get the system's whatnot.
	BED->SLP_memory.set_globals_to_this();
	evalProg_mp(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	BED->SLP_memory.set_globals_null();
	
	// evaluate the patch
	patch_eval_mp(    patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	// note that you can only really do this AFTER you are done calling other evaluators.
	// set parVals & parDer correctly
	
	// i.e. these must remain here, or below.  \/
	change_size_point_mp(parVals, 1);
	change_size_vec_mp(parDer, 1);
	change_size_mat_mp(Jp, BED->num_variables, 1); Jp->rows = BED->num_variables; Jp->cols = 1;
	
	for (int ii=0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);
	
	
	parVals->size = parDer->size = 1;
	set_mp(&parVals->coord[0], pathVars); // s = t
	set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	///////// / / / /  /   /
	// combine everything
	///////// / / / /  /   /
	
	BED->randomizer()->randomize(AtimesF,AtimesJ,temp_function_values,temp_jacobian_functions,&current_variable_values->coord[0]);

	//perform the randomization multiplications
//	mat_mul_mp(AtimesJ,BED->randomizer_matrix,temp_jacobian_functions);
//	mul_mat_vec_mp(AtimesF,BED->randomizer_matrix, temp_function_values ); // set values of AtimesF (A is randomization matrix)
	
	for (int ii=0; ii<AtimesF->size; ii++)  // for each function, after (real orthogonal) randomization
		set_mp(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	
	
	
	
	
	
	//////////////
	//
	// function values and Jp for the linears
	//
	////////////////////
	
	
	offset = BED->num_variables-BED->patch.num_patches-BED->num_linears;
	for (int ii=0; ii<BED->num_linears; ++ii) {
		dot_product_mp(temp,BED->current_linear[ii],current_variable_values);
		
		mul_mp(&funcVals->coord[ii+offset],one_minus_s,temp);
		neg_mp(&Jp->entry[offset+ii][0],temp);
		
		dot_product_mp(temp2,BED->old_linear[ii],current_variable_values);
		
		mul_mp(temp,temp2,BED->gamma);
		add_mp(&Jp->entry[offset+ii][0],&Jp->entry[offset+ii][0],temp); // Jp = -curr_lin_val + gamma*old_lin_val
		
		mul_mp(temp2,temp2,gamma_s);
		add_mp(&funcVals->coord[ii+offset],&funcVals->coord[ii+offset],temp2); // f = (1-s)*curr_lin_val + gamma*s*old_lin_val
	}
	
	
	//////////////
	//
	//set the function PATCH values
	//
	////////////////////
	
	
	
	offset = BED->num_variables-BED->patch.num_patches;
	for (int ii=0; ii<BED->patch.num_patches; ii++)
		set_mp(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//////////////
	//
	// set the JACOBIAN values.
	//
	////////////////////
	
	//first, the entries related to the functions
	
	for (int ii = 0; ii < BED->randomizer()->num_rand_funcs(); ii++)
		for (int jj = 0; jj < BED->num_variables; jj++)
			set_mp(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	
	
	
	//////////////
	//
	// Entries for the linears we are homotoping
	//
	////////////////////
	
	
	offset = BED->num_variables-BED->patch.num_patches-BED->num_linears;
	for (int ii=0; ii<BED->num_linears; ++ii) {
		for (int jj=0; jj<BED->num_variables; jj++) {
			mul_mp(temp,&BED->old_linear[ii]->coord[jj],gamma_s);
			mul_mp(temp2,&BED->current_linear[ii]->coord[jj],one_minus_s);
			add_mp(&Jv->entry[offset+ii][jj], temp,temp2);
		}
	}
	
	
	
	
	
	
	//////////////
	//
	// the entries in the jacobian for the patch equations.
	//
	////////////////////
	
	offset = BED->num_variables - BED->patch.num_patches;
	for (int ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{  // Jv = Jv_Patch
		for (int jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	// done!  yay!
	
	if (BED->verbose_level()==12) {
		//uncomment to see screen output of important variables at each solve step.
		print_comp_matlab(pathVars, "t");
		printf("BED->num_linears = %d\n",BED->num_linears);
		print_point_to_screen_matlab(parVals,"parVals_mp");
		print_matrix_to_screen_matlab( AtimesJ,"R*jac_mp");
		print_point_to_screen_matlab(current_variable_values,"currvars_mp");
		for (int ii=0; ii<BED->num_linears; ii++) {
			print_point_to_screen_matlab(BED->current_linear[ii],"new_mp");
			print_point_to_screen_matlab(BED->old_linear[ii],"old_mp");
		}
		print_point_to_screen_matlab(funcVals,"F_mp");
		print_matrix_to_screen_matlab(Jv,"Jv_mp");
		print_matrix_to_screen_matlab(Jp,"Jp_mp");
		std::cout << BED->randomizer() << std::endl;
	}
	
	
	
	
	
	
	
	
	clear_mp(one_minus_s);
	clear_mp(gamma_s);
	clear_mp(temp);
	clear_mp(temp2);
	
	
	clear_vec_mp(patchValues);
	clear_vec_mp(temp_function_values);
	clear_vec_mp(AtimesF);
	
	
	clear_mat_mp(temp_jacobian_functions);
	clear_mat_mp(temp_jacobian_parameters);
	clear_mat_mp(Jv_Patch);
	clear_mat_mp(AtimesJ);
	
#ifdef printpathmultilintolin
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize_mp(&dehommed,current_variable_values);
	mpf_out_str (BED->FOUT, 10, 15, pathVars->r);
	fprintf(BED->FOUT," ");
	mpf_out_str (BED->FOUT, 10, 15, pathVars->i);
	fprintf(BED->FOUT," ");
	for (int ii=0; ii<BED->num_variables-1; ++ii) {
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





int multilintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
{
	multilintolin_eval_data_d *BED_d = NULL;
	multilintolin_eval_data_mp *BED_mp = NULL;
	
	*out_prec = in_prec;
	
	
	
	if (in_prec < 64)
	{ // compute out_d
		multilintolin_eval_data_d *BED_d = (multilintolin_eval_data_d *)ED_d;
		
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
		multilintolin_eval_data_mp *BED_mp = (multilintolin_eval_data_mp *)ED_mp;
		
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




int change_multilintolin_eval_prec(void const *ED, int new_prec)
{
	
	multilintolin_eval_data_mp *BED = (multilintolin_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	
	BED->SLP->precision = new_prec;
	// change the precision for the patch
	changePatchPrec_mp(new_prec, &BED->patch);
	
	
	if (new_prec != BED->curr_prec){
		
		if (BED->verbose_level() >=8){
			std::cout << color::brown();
			printf("prec  %lu\t-->\t%d\n",BED->curr_prec, new_prec);
			std::cout << color::console_default();
		}
		
		BED->curr_prec = new_prec;
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		int ii;
		for(ii=0; ii<BED->num_linears; ++ii){
			change_prec_point_mp(BED->current_linear[ii],new_prec);
			vec_cp_mp(BED->current_linear[ii], BED->current_linear_full_prec[ii]);
			
			change_prec_point_mp(BED->old_linear[ii],new_prec);
			vec_cp_mp(BED->old_linear[ii], BED->old_linear_full_prec[ii]);
		}
		BED->randomizer()->change_prec(new_prec);
		
	}
	
	return 0;
}





int check_issoln_multilintolin_d(endgame_data_t *EG,
								 tracker_config_t *T,
								 void const *ED)
{
	
	
	
	
	if (EG->last_approx_prec<64) {
		if (EG->last_approx_d->size==0) {
			return 0;
		}
	}
	else{
		if (EG->last_approx_mp->size==0) {
			return 0;
		}
	}
	
	
	multilintolin_eval_data_d *BED = (multilintolin_eval_data_d *)ED; // to avoid having to cast every time
	
	BED->SLP_memory.set_globals_to_this();
	
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
	for (int ii = 0; (ii < BED->SLP->numFuncs) && isSoln; ii++)
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


int check_issoln_multilintolin_mp(endgame_data_t *EG,
								  tracker_config_t *T,
								  void const *ED)
{
	
	
	multilintolin_eval_data_mp *BED = (multilintolin_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	if (EG->last_approx_prec<64) {
		if (EG->last_approx_d->size==0) {
			return 0;
		}
	}
	else{
		if (EG->last_approx_mp->size==0) {
			return 0;
		}
	}
	
	
	
	BED->SLP_memory.set_globals_to_this();

	
	for (int ii = 0; ii < T->numVars; ii++)
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
	//	multilin_to_lin_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, ED);
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
	for (int ii = 0; ii < BED->SLP->numFuncs && isSoln; ii++)
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
		mpf_out_str (NULL, 10, 8, max_rat);
		printf("\n");
	}
	
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	clear_vec_mp(f);
	
	BED->SLP_memory.set_globals_null();
	return isSoln;
	
}















