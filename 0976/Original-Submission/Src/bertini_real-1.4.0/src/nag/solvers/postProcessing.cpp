#include "nag/solvers/postProcessing.hpp"




void SolverOutput::post_process(post_process_t *endPoints, int num_pts_to_check,
								 preproc_data *preProcData, tracker_config_t *T,
								 const SolverConfiguration & solve_options)
{
	
	int num_nat_vars = num_natural_vars;
	
	
	// sets the multiplicity and solution number in the endPoints data
	//direct from the bertini library:
	findMultSol(endPoints, num_pts_to_check, num_nat_vars, preProcData, T->final_tol_times_mult);
	
	
	int num_singular_solns = 0, num_real_solns = 0, num_finite_solns = 0;
	
		
	if (0) { /// this code represents an attempt to get bertini's own postprocessing to do the work for me.  since in bertini_real, we only care about the first group (and in face some variables aren't even listed), the bertini PP fails to correctly set the isFinite bit.  so this remains in an impossible-to-get-to if-else.
		
		std::vector<int> origErrorIsInf;
		origErrorIsInf.resize(num_pts_to_check);
		
		std::vector<double> origErrorEst;
		origErrorEst.resize(num_pts_to_check);
		
		point_d *dehomPoints_d = T->MPType == 1 ? NULL : (point_d *)bmalloc(num_pts_to_check * sizeof(point_d));
		point_mp *dehomPoints_mp = T->MPType == 0 ? NULL : (point_mp *)bmalloc(num_pts_to_check * sizeof(point_mp));
		
		
		// setup the dehomogenized points and error estimates
		for (int ii = 0; ii < num_pts_to_check; ii++)
		{ // setup dehom
			if (endPoints[ii].sol_prec < 64)
			{ // setup _d
				init_point_d(dehomPoints_d[ii], 0);
				getDehomPoint_comp_d(dehomPoints_d[ii], &origErrorIsInf[ii], &origErrorEst[ii], endPoints[ii].sol_d, num_nat_vars+1, preProcData, endPoints[ii].accuracy_estimate);
			}
			else
			{ // setup _mp
				initMP(endPoints[ii].sol_prec);
				init_point_mp(dehomPoints_mp[ii], 0);
				getDehomPoint_comp_mp(dehomPoints_mp[ii], &origErrorIsInf[ii], &origErrorEst[ii], endPoints[ii].sol_mp, num_nat_vars+1, preProcData, endPoints[ii].accuracy_estimate);
			}
		}
		if (preProcData->num_var_gp > 0)
		{ // find the finite solutions
			findFiniteSol(endPoints, dehomPoints_d, dehomPoints_mp, num_pts_to_check, 1000, preProcData, T->finiteThreshold);
		}
		else
		{ // set finite as -1
			for (int ii = 0; ii < num_pts_to_check; ii++)
				endPoints[ii].isFinite = -1;
		}
		
		// determine which ones are real
		findRealSol(endPoints, dehomPoints_d, dehomPoints_mp, num_pts_to_check, 1000, preProcData, T->real_threshold);
		
		// determine which ones are singular
		findSingSol(endPoints, dehomPoints_d, dehomPoints_mp, num_pts_to_check, 1000, preProcData, T->cond_num_threshold, T->final_tol_times_mult, 0);
		
		
		for (int ii=0; ii<num_pts_to_check; ii++) {
			if (endPoints[ii].isFinite==1) {
				num_finite_solns++;
			}
			
			if (endPoints[ii].isSing==1) {
				num_singular_solns++;
			}
			
			if (endPoints[ii].isReal==1) {
				num_real_solns++;
			}
			
		}
		
		
		for (int ii = 0; ii < num_pts_to_check; ii++)
		{
			if (endPoints[ii].sol_prec < 64)
			{
				clear_point_d(dehomPoints_d[ii]);
			}
			else
			{
				clear_point_mp(dehomPoints_mp[ii]);
			}
		}
		if (T->MPType != 1)
			free(dehomPoints_d);
		if (T->MPType != 0)
			free(dehomPoints_mp);
		
		
	}
	else{
		//sets the singularity flag in endPoints.
		//custom, derived from bertini's analagous call.
		num_singular_solns = BRfindSingularSolns(endPoints, num_pts_to_check, num_nat_vars, T);
		
		//sets the finite flag in endPoints.
		//custom, derived from bertini's analagous call.
		num_finite_solns = BRfindFiniteSolns(endPoints, num_pts_to_check, num_nat_vars, T);
		
		
		num_real_solns   = BRfindRealSolns(endPoints,num_pts_to_check,num_nat_vars,T);
	}
			
				
	
	
	
	
	if (solve_options.verbose_level()>=3)
		printf("%d finite solutions, %d singular solutions, %d real solutions\n",num_finite_solns, num_singular_solns, num_real_solns);
	
	
	if (solve_options.verbose_level()>=3) {
		for (int ii=0; ii<num_pts_to_check; ++ii) {
			//		int success;      // success flag
			//		int multiplicity; // multiplicity
			//		int isReal;       // real flag:  0 - not real, 1 - real
			//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
			//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
			printf("solution %d, success %d, multi %d, isFinite %d, isSing %d, isReal %d, cycle_num %d\n",ii,endPoints[ii].success,endPoints[ii].multiplicity,endPoints[ii].isFinite,endPoints[ii].isSing,endPoints[ii].isReal,endPoints[ii].cycle_num);
			
		}
	}
	

	
	std::vector<  std::pair< long long,long long > > soln_indices;
	for (int ii=0; ii<num_pts_to_check; ii++) {
		soln_indices.push_back(std::pair<long long,long long>(ii,endPoints[ii].path_num));
	}
	
	
	// sort the indices into the post_process_t array, based on the path numbers.  this restores the order permuted by multiple processors.
	std::sort(soln_indices.begin(), soln_indices.end(),
			  boost::bind(&std::pair<long long, long long>::second, _1) <
			  boost::bind(&std::pair<long long, long long>::second, _2));
	
	
	
	Vertex temp_vertex;
	change_size_point_mp(temp_vertex.point(),num_variables);
	(temp_vertex.point())->size = num_variables;
	
	
	//first, lets take care of the multiplicity 1 solutions
	for (int ii=0; ii<num_pts_to_check; ii++) {
		int curr_ind = soln_indices[ii].first;
		
		if (endPoints[curr_ind].multiplicity!=1) {
			continue;
		}
		
		
		endpoint_to_vec_mp(temp_vertex.point(), &endPoints[curr_ind]);
		
		SolutionMetadata meta;

		meta.set_real(endPoints[curr_ind].isReal);
		meta.set_finite(endPoints[curr_ind].isFinite);
		meta.set_singular(endPoints[curr_ind].isSing);
		meta.set_multiplicity(endPoints[curr_ind].multiplicity);
		meta.set_successful(endPoints[curr_ind].success);
		meta.SetCycleNumber(endPoints[curr_ind].cycle_num);
		meta.set_output_index(this->num_vertices_);
		meta.add_input_index(endPoints[curr_ind].path_num);
		
		add_solution(temp_vertex, meta);
	}
	
	
	
	
	
	//now we deal with the multiplicity > 1 solutions
	for (int ii=0; ii<num_pts_to_check; ii++) {
		int curr_ind = soln_indices[ii].first;
		
		if (endPoints[curr_ind].multiplicity==1) {
			continue;
		}
		
		if (endPoints[curr_ind].multiplicity<1) {
			endpoint_to_vec_mp(temp_vertex.point(), &endPoints[curr_ind]);
			continue;
		}
		
		
		
		if ( find(occuring_multiplicities.begin(),occuring_multiplicities.end(),endPoints[curr_ind].multiplicity)==occuring_multiplicities.end()) {
			occuring_multiplicities.push_back(endPoints[curr_ind].multiplicity);
		}
		
		endpoint_to_vec_mp(temp_vertex.point(), &endPoints[curr_ind]);
		
		SolutionMetadata meta;
		meta.set_finite(endPoints[curr_ind].isFinite);
		meta.set_singular(endPoints[curr_ind].isSing);
		meta.set_multiplicity(endPoints[curr_ind].multiplicity);
		meta.set_successful(endPoints[curr_ind].success);
		meta.SetCycleNumber(endPoints[curr_ind].cycle_num);

		for (int jj=0; jj<num_pts_to_check; jj++) {
			int inner_ind = soln_indices[jj].first;
			
			if (endPoints[inner_ind].sol_num==endPoints[curr_ind].sol_num) {
				meta.add_input_index(endPoints[inner_ind].path_num);
			}
		}
		meta.set_output_index(this->num_vertices());
		
		add_solution(temp_vertex, meta);
	}
	
	
	std::sort(occuring_multiplicities.begin(), occuring_multiplicities.end());
	for (unsigned int ii=0; ii<num_vertices(); ii++) {
		for (auto jj = metadata[ii].input_index.begin(); jj!=metadata[ii].input_index.end(); ++jj) {
			ordering.push_back(std::pair<long long, long long>(ii,*jj));
		}
	}
	
	// sort the ordering.  now properly respects the vertex set WRT multiplicity.
	std::sort(ordering.begin(), ordering.end(),
			  boost::bind(&std::pair<long long, long long>::second, _1) <
			  boost::bind(&std::pair<long long, long long>::second, _2));

	return;
}






int BRfindSingularSolns(post_process_t *endPoints, int num_sols, int num_vars,
												tracker_config_t *T )
{
	
	if (num_vars<=1) {
		std::cout << "requesting to find singular solutions, but using " << num_vars << " variables." << std::endl;
	}
	
	
	int sing_count=0;
	
	for (int ii = 0; ii < num_sols; ii++){
		if ( (endPoints[ii].cond_est >  T->cond_num_threshold) || (endPoints[ii].cond_est < 0.0) || (endPoints[ii].multiplicity!=1))
			endPoints[ii].isSing = 1;
		else
			endPoints[ii].isSing = 0;
		
		if (endPoints[ii].isSing)
		{
			sing_count++;
		}
	}
	
	return sing_count;
}


int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
											tracker_config_t *T )
{
	
	if (num_vars<=1) {
		std::cout << "requesting to find finite solutions, but using " << num_vars << " variables." << std::endl;
	}
	
	int finite_count=0;
	
	
	//initialize temp stuffs
	comp_d dehom_coord_recip_d;
	comp_mp dehom_coord_recip_mp; init_mp(dehom_coord_recip_mp);
	vec_d dehom_d;   init_vec_d(dehom_d,num_vars-1);   dehom_d->size = num_vars-1;
	vec_mp dehom_mp; init_vec_mp(dehom_mp,num_vars-1); dehom_mp->size = num_vars-1;
	
	
	
	for (int ii = 0; ii < num_sols; ii++){
		if (endPoints[ii].sol_prec<64) {
			set_d(dehom_coord_recip_d,endPoints[ii].sol_d[0]);
			recip_d(dehom_coord_recip_d,dehom_coord_recip_d);
			for (int jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_d(&dehom_d->coord[jj],dehom_coord_recip_d,endPoints[ii].sol_d[jj+1])
			}
			
			if (infNormVec_d(dehom_d) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
//			print_point_to_screen_matlab(dehom_d,"soln");
			
		}
		else // high precision, do mp
		{
			change_prec_point_mp(dehom_mp,endPoints[ii].sol_prec);
			setprec_mp(dehom_coord_recip_mp,endPoints[ii].sol_prec);
			recip_mp(dehom_coord_recip_mp,endPoints[ii].sol_mp[0]);
			for (int jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_mp(&dehom_mp->coord[jj],dehom_coord_recip_mp,endPoints[ii].sol_mp[jj+1])
			}
			
			if (infNormVec_mp(dehom_mp) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
//			print_point_to_screen_matlab(dehom_mp,"soln");
		}
	}
	
	clear_vec_d(dehom_d);
	clear_vec_mp(dehom_mp);
	clear_mp(dehom_coord_recip_mp);
	
	return finite_count;
}


int BRfindRealSolns(post_process_t *endPoints, int num_sols, int num_vars,
					  tracker_config_t *T )
{
	if (num_vars<=1) {
		std::cout << "requesting to find real solutions, but using " << num_vars << " variables." << std::endl;
	}
	
	
	int real_count=0;
	
	
	//initialize temp stuffs
	comp_d dehom_coord_recip_d;
	comp_mp dehom_coord_recip_mp; init_mp(dehom_coord_recip_mp);
	vec_d dehom_d;   init_vec_d(dehom_d,num_vars-1);   dehom_d->size = num_vars-1;
	vec_mp dehom_mp; init_vec_mp(dehom_mp,num_vars-1); dehom_mp->size = num_vars-1;
	
	
	
	for (int ii = 0; ii < num_sols; ii++){
		if (endPoints[ii].sol_prec<64) {
			set_d(dehom_coord_recip_d,endPoints[ii].sol_d[0]);
			recip_d(dehom_coord_recip_d,dehom_coord_recip_d);
			for (int jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_d(&dehom_d->coord[jj],dehom_coord_recip_d,endPoints[ii].sol_d[jj+1])
			}
			
			endPoints[ii].isReal = checkForReal_d(dehom_d, T->real_threshold);
			
		}
		else // high precision, do mp
		{
			change_prec_point_mp(dehom_mp,endPoints[ii].sol_prec);
			setprec_mp(dehom_coord_recip_mp,endPoints[ii].sol_prec);
			recip_mp(dehom_coord_recip_mp,endPoints[ii].sol_mp[0]);
			for (int jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_mp(&dehom_mp->coord[jj],dehom_coord_recip_mp,endPoints[ii].sol_mp[jj+1])
			}
			
			endPoints[ii].isReal = checkForReal_mp(dehom_mp, T->real_threshold);
			
		}
		
		if (endPoints[ii].isReal) {
			real_count++;
		}
	}
	
	clear_vec_d(dehom_d);
	clear_vec_mp(dehom_mp);
	clear_mp(dehom_coord_recip_mp);
	
	return real_count;
}






void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG)
{
	
	
	//	printf("endgame2endpoint\n");
	int num_vars;
	endPoint->path_num = EG->pathNum;
	
	endPoint->sol_prec = EG->prec;
	endPoint->cond_est = EG->condition_number;
	endPoint->final_t = EG->t_val_at_latest_sample_point_d;//???
	endPoint->first_increase = EG->first_increase;
	
	
	if (EG->prec==52) {
		num_vars = EG->PD_d.point->size;
		endPoint->sol_d  = (comp_d *)br_malloc(num_vars * sizeof(comp_d));
		endPoint->sol_mp = NULL;
		
		for (int ii=0; ii<num_vars; ii++) {
			endPoint->sol_d[ii]->r = EG->PD_d.point->coord[ii].r;
			endPoint->sol_d[ii]->i = EG->PD_d.point->coord[ii].i;
		}
		
		endPoint->size_sol = num_vars;
		endPoint->function_resid_d = EG->function_residual_d;  // the function residual
		endPoint->newton_resid_d = EG->latest_newton_residual_d;
		endPoint->cycle_num = EG->PD_d.cycle_num;
		endPoint->accuracy_estimate = EG->error_at_latest_sample_point_d;//
		
	}
	else
	{
		num_vars = EG->PD_mp.point->size;
		
		endPoint->sol_d  = NULL;
		endPoint->sol_mp = (comp_mp *)br_malloc(num_vars * sizeof(comp_mp));
		for (int ii=0; ii<num_vars; ii++) {
			init_mp2(endPoint->sol_mp[ii],EG->prec);
			mpf_set(endPoint->sol_mp[ii]->r,EG->PD_mp.point->coord[ii].r);
			mpf_set(endPoint->sol_mp[ii]->i,EG->PD_mp.point->coord[ii].i);
		}
		endPoint->size_sol = num_vars;
		mpf_init2(endPoint->function_resid_mp, EG->prec); mpf_init2(endPoint->newton_resid_mp, EG->prec);
		
		mpf_set(endPoint->function_resid_mp,EG->function_residual_mp); //this is undoubtedly incorrect
		mpf_set(endPoint->newton_resid_mp,EG->latest_newton_residual_mp);
		endPoint->cycle_num = EG->PD_mp.cycle_num;
		endPoint->accuracy_estimate = mpf_get_d(EG->error_at_latest_sample_point_mp);
	}
	
	
	if (EG->retVal==0 || EG->retVal==-22 || EG->retVal==-50 || EG->retVal==-21) {//EG->retVal==-50
		endPoint->success = 1;
	}
	else if (EG->retVal == retVal_sharpening_failed){
		endPoint->success = retVal_sharpening_failed;
	}
	else if (EG->retVal == retVal_sharpening_singular_endpoint){
		endPoint->success = retVal_sharpening_singular_endpoint;
	}
	else{
		std::cout << "setting endPoint->success = -1 because retVal==" << EG->retVal << std::endl;
		endPoint->success = -1;
	}
	
	
	
	endPoint->sol_num = -1; // set up for post-processing
	endPoint->multiplicity = 1;
	endPoint->isFinite = 0;
	//	// the post_process_t structure is used in post-processing //
	//	typedef struct
	//	{
	//		int path_num;     // path number of the solution
	//		int sol_num;      // solution number
	//		comp_d  *sol_d;   // solution
	//		comp_mp *sol_mp;
	//		int sol_prec;     // precision of the solution
	//		int size_sol;     // the number of entries in sol
	//		double function_resid_d;  // the function residual
	//		mpf_t  function_resid_mp;
	//		double cond_est;  // the estimate of the condition number
	//		double newton_resid_d;    // the newton residual
	//		mpf_t  newton_resid_mp;
	//		double final_t;   // the final value of time
	//		double accuracy_estimate; // accuracy estimate between extrapolations
	//		double first_increase;    // time value of the first increase in precision
	//		int cycle_num;    // cycle number used in extrapolations
	//		int success;      // success flag
	//		int multiplicity; // multiplicity
	//		int isReal;       // real flag:  0 - not real, 1 - real
	//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
	//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
	//	} post_process_t;
	
	
	
	//	typedef struct
	//	{
	//		int prec;
	//		point_data_d PD_d;
	//		point_data_mp PD_mp;
	//
	//		int last_approx_prec;       // precision of the last approximation
	//		point_d last_approx_d;      // last approximation to the end point
	//		point_mp last_approx_mp;    // last approximation to the end point
	//
	//		int retVal;
	//		int pathNum;
	//		int codim;
	//		double first_increase;
	//		double condition_number;
	//		double function_residual_d;
	//		mpf_t  function_residual_mp;
	//		double latest_newton_residual_d;
	//		mpf_t  latest_newton_residual_mp;
	//		double t_val_at_latest_sample_point_d;
	//		mpf_t  t_val_at_latest_sample_point_mp;
	//		double error_at_latest_sample_point_d;
	//		mpf_t  error_at_latest_sample_point_mp;
	//	} endgame_data_t;
	
	
	
	
}




void endpoint_to_vec_mp(vec_mp veccie, post_process_t *endPoint)
{
	change_size_point_mp(veccie, endPoint->size_sol);  veccie->size = endPoint->size_sol;
	change_prec_point_mp(veccie,MAX(endPoint->sol_prec,64));
	
	if (endPoint->sol_prec<64) {
		//copy out of the double structure.
		for (int jj=0; jj<endPoint->size_sol; jj++)
			d_to_mp(&veccie->coord[jj],endPoint->sol_d[jj]);
		
	}
	else{
		//copy out of the mp structure.
		for (int jj=0; jj<endPoint->size_sol; jj++)
			set_mp(&veccie->coord[jj],endPoint->sol_mp[jj]);
		
	}
}
