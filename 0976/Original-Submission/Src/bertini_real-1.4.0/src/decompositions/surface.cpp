#include "decompositions/surface.hpp"












void Surface::main(VertexSet & V,
                                 const WitnessSet & W_surf,
                                 vec_mp *pi,
                                 BertiniRealConfig & program_options,
                                 SolverConfiguration & solve_options)
{
    
	
#ifdef functionentry_output
	std::cout << "surface_main" << std::endl;
#endif
	
	Decomposition::copy_data_from_witness_set(W_surf);
	
	
    
	add_projection(pi[0]); // add to *this
	add_projection(pi[1]); // add to *this
	
	
	
	
	
	
    
	//copy the patch over into this object
	this->copy_patches(W_surf);
	

	
	
	beginning_stuff( W_surf, program_options, solve_options);
	
	
    
    // get the witness points for the critical curve.
    WitnessSet W_critcurve;
	std::map< int, WitnessSet> higher_multiplicity_witness_sets;
	
	
    compute_critcurve_witness_set(W_critcurve,
								  higher_multiplicity_witness_sets,
                                  W_surf,
                                  0,
                                  program_options,
                                  solve_options);
    
	
	
	
	
	
	
	
	
	
	// get the critical points and the sphere intersection points for the critical curve
    WitnessSet W_critcurve_crit;
    compute_critcurve_critpts(W_critcurve_crit, // the computed value
                              W_critcurve,
                              0,
                              program_options,
                              solve_options);
    
    
    this->crit_curve_.add_witness_set(W_critcurve_crit,Critical,V);
    
    if (0)
    {
	    WitnessSet W_critcurve1;
		std::map< int, WitnessSet> higher_multiplicity_witness_sets1;
		
		
	    compute_critcurve_witness_set(W_critcurve1,
									  higher_multiplicity_witness_sets1,
	                                  W_surf,
	                                  1,
	                                  program_options,
	                                  solve_options);
		
		
	    compute_critcurve_critpts(W_critcurve_crit, // the computed value
	                              W_critcurve1,
	                              1,
	                              program_options,
	                              solve_options);
    }
	
	
	
	///////////////////////////////
	
	WitnessSet points_which_needed_no_deflation;
	std::map< SingularObjectMetadata, WitnessSet > split_sets;
	
	deflate_and_split(split_sets,
					  higher_multiplicity_witness_sets,
					  points_which_needed_no_deflation,
					  program_options,
					  solve_options);
	
	
	
	WitnessSet W_singular_crit;
	compute_singular_crit(W_singular_crit,
						  split_sets,
						  V,
						  program_options,
						  solve_options);
	///////////////////////////////
	
	
	
	// merge together the critical points from both the critical curve and the sphere intersection curve.
    WitnessSet W_total_crit;
    
	W_total_crit.merge(W_critcurve_crit,&solve_options.T);
	W_total_crit.merge(W_singular_crit,&solve_options.T);
	
				
				
	if (have_sphere()) {
		std::cout << "sorting for inside sphere" << std::endl;
		W_total_crit.sort_for_inside_sphere(sphere_radius(), sphere_center());
	}
	else
	{
		this->compute_sphere_bounds(W_total_crit); // sets the radius and center in this Decomposition.  Must propagate to the constituent decompositions as well.   fortunately, i have a method for that!!!
	}
	crit_curve_.copy_sphere_bounds(*this); // copy the bounds into the critcurve.
				
				
				
	
	std::cout << color::bold('m') << "intersecting critical curve with sphere" << color::console_default() << std::endl;
	
	W_critcurve_crit.set_input_filename("input_critical_curve");
	
	
	WitnessSet W_sphere_intersection;
	W_sphere_intersection.set_input_filename("input_critical_curve");
	// now get the sphere intersection critical points and ends of the interval
	crit_curve_.get_sphere_intersection_pts(&W_sphere_intersection,  // the returned value
									  W_critcurve,       // all else here is input
									  program_options,
									  solve_options);
	
	
	
	W_sphere_intersection.sort_for_real(&solve_options.T);
	W_sphere_intersection.sort_for_unique(&solve_options.T);
	
	
	crit_curve_.add_witness_set(W_sphere_intersection,Critical,V);
	
	W_total_crit.merge(W_sphere_intersection,&solve_options.T);
				
				
	
	
	// now we get the critical points for the sphere intersection curve.
	
    
    
    
    // make the input file
	create_sphere_system(W_surf.input_filename(),
                         "input_surf_sphere",
                         sphere_radius(),
                         sphere_center(),
                         W_surf);
    
    // get witness points
	WitnessSet W_intersection_sphere;
	compute_sphere_witness_set(W_surf,
                               W_intersection_sphere, // output
                               program_options,
                               solve_options);
	
    
    // compute critical points
    WitnessSet W_sphere_crit;
	compute_sphere_crit(W_intersection_sphere,
                        W_sphere_crit, // output
                        program_options,
                        solve_options);
    
	
	
    this->sphere_curve_.add_witness_set(W_sphere_crit,Critical,V);
    
	
	
	W_total_crit.merge(W_sphere_crit,&solve_options.T);
	
	W_total_crit.set_input_filename("W_total_crit_nonexistant_filename");
	W_total_crit.sort_for_unique(&solve_options.T);
    

    
    
	
	
	
	
    
	
	compute_critical_curve(W_critcurve, // all input.
                           W_total_crit,
                           V,
                           program_options,
                           solve_options);
	
		
	
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	


	
    // actually perform the interslice on the bounding sphere curve.
	compute_bounding_sphere(W_intersection_sphere, // the witness points we will track from
                            W_total_crit,          // the critical points for the both sphere and critical curve.
                            V,                     // vertex set.  it goes almost everywhere.
                            program_options, solve_options); // configuration
    
	
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	
	
	compute_singular_curves(W_total_crit,
							split_sets,
							V,
							program_options, solve_options);
	
	
	
    
    
    
    
    
    
	
	
	
	
	
	
	
	// compute the downstairs crit and midpoints for slicing
	vec_mp midpoints_downstairs, crit_downstairs;  std::vector< int > index_tracker;
	init_vec_mp(midpoints_downstairs,0); init_vec_mp(crit_downstairs,0);
    
    V.compute_downstairs_crit_midpts(W_total_crit, crit_downstairs, midpoints_downstairs, index_tracker, pi[0],&solve_options.T);
	
	
	
	
	SetCritSliceValues(crit_downstairs);
	
	if (program_options.verbose_level()>=0) {
        std::cout << color::green() << "the pi[0] critical projection values at which we will be slicing:\n\n" << color::console_default();
		print_point_to_screen_matlab(crit_downstairs,"crit_down");
        
	}
	
    if (program_options.verbose_level()>=2) {
		W_total_crit.print_to_screen();
	}
	
	
	
	
	
	
	// get the midpoint slices
	program_options.merge_edges(true);
	compute_slices(W_surf, V,
                   midpoints_downstairs, this->mid_slices_,
                   program_options, solve_options, "mid");
	
	
	
	//incremental output
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	// get the critical slices
	
	program_options.merge_edges(true);
	compute_slices(W_surf, V,
                   crit_downstairs, this->crit_slices_,
                   program_options, solve_options, "crit");
	
    
	//incremental output
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	//connect the dots - the final routine
	connect_the_dots(V, program_options, solve_options);
	
	
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	clear_vec_mp(crit_downstairs); clear_vec_mp(midpoints_downstairs);
    
	
	if (program_options.verbose_level()>=0) {
		std::cout << "decomposed surface has " << num_faces() << " faces" << std::endl;
	}
	
	
	
	
}








void Surface::beginning_stuff(const WitnessSet & W_surf,
                                            BertiniRealConfig & program_options,
                                            SolverConfiguration & solve_options)
{
	
#ifdef functionentry_output
	std::cout << "surface::beginning_stuff" << std::endl;
#endif
	
	if (1) {
		// perform an isosingular deflation
		// note: you do not need witness_data to perform isosingular deflation
		if (program_options.verbose_level()>=2)
			printf("performing isosingular deflation\n");
		
		

		boost::filesystem::path temp_path = program_options.input_filename();
		
		std::stringstream converter;
		converter << "_dim_" << W_surf.dimension() << "_comp_" << W_surf.component_number() << "_deflated";
		temp_path += converter.str();
		
		
		
		
		std::set<unsigned int> zeroonly;
		zeroonly.insert(0);
		W_surf.write_dehomogenized_coordinates("witness_points_dehomogenized",zeroonly); // write the points to file
		
		
		int num_deflations, *deflation_sequence = NULL;
		
		isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
							  program_options.input_filename(),
							  "witness_points_dehomogenized",
							  temp_path, // the desired output name
							  program_options.max_deflations());
		free(deflation_sequence);
		
		
		program_options.set_input_deflated_filename(temp_path);
		converter.clear(); converter.str("");
	}
	else {
		program_options.set_input_deflated_filename(program_options.input_filename());
	}
	
	
	
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename());
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	if (1) {
		if (program_options.verbose_level()>=2)
			printf("checking if component is self-conjugate\n");
		checkSelfConjugate( W_surf.point(0),program_options, W_surf.input_filename());
		
		//regenerate the various files, since we ran bertini since then and many files were deleted.
		parse_input_file(program_options.input_deflated_filename());
	}
	
	
	
	if (program_options.user_sphere()) {
		read_sphere(program_options.bounding_sphere_filename());
	}
	
	
	randomizer()->setup(W_surf.num_variables()-1-this->dimension(),solve_options.PPD.num_funcs);
	
	

	
	
	
	
	if (!verify_projection_ok(W_surf,
                              this->randomizer(),
                              &(this->pi()),
                              solve_options))
	{
		std::cout << "the projections being used appear to suffer rank deficiency with Jacobian matrix..." << std::endl;
		mypause();
	}
    
	
}














void Surface::compute_critcurve_witness_set(WitnessSet & W_critcurve,
														  std::map<int, WitnessSet> & higher_multiplicity_witness_sets,
                                                          const WitnessSet & W_surf,
                                                          int pi_ind,
                                                          BertiniRealConfig & program_options,
                                                          SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_critcurve_witness_set" << std::endl;
#endif
	
	
	
	std::cout << color::bold('m') << "computing witness points for the critical curve" << color::console_default() << std::endl;
    
	
	
	// find witness points on the critical curve.
	
	
	
	solve_options.use_gamma_trick = 0;
	
	
	NullspaceConfiguration ns_config;
    
	SolverOutput solve_out;
	
	compute_crit_nullspace(solve_out,	// the returned value
                           W_surf,            // input the original witness set
                           this->randomizer(),
                           &(this->pi()),
                           2,  // dimension of ambient complex object
                           2,   //  target dimension to find
                           1,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	
	
	solve_out.get_nonsing_finite_multone(W_critcurve);
	solve_out.get_patches_linears(W_critcurve);
	solve_out.copy_names(W_critcurve);
	
	
	solve_out.get_multpos_full(higher_multiplicity_witness_sets);
	//now we have a map of multiplicities and witness sets.  for each point of each multiplicity, we need to perform iso. defl.  this will enable us to get ahold of the singular curves.
	

	if (program_options.verbose_level()>=2) {
		for (auto iter = higher_multiplicity_witness_sets.begin(); iter!=higher_multiplicity_witness_sets.end(); iter++) {
			std::cout << "found " << iter->second.num_points() << " points of multiplicity " << iter->first << std::endl;
			
			iter->second.print_to_screen();
		}
	}
	

	
	
	
	
	
	
	//this uses both pi[0] and pi[1] to compute witness points
	
	W_critcurve.write_dehomogenized_coordinates("W_curve"); // write the points to file
	
	W_critcurve.set_input_filename("input_critical_curve");
	
	// this system describes the system for the critical curve
	create_nullspace_system("input_critical_curve",
                            boost::filesystem::path(program_options.called_dir()) / program_options.input_deflated_filename(),
                            program_options, &ns_config);
    
    
    
    
    // adjust the size of the linears to match the number of variables.  this should be a method in the witness set data type.
    
    for (unsigned int ii=0; ii<W_critcurve.num_linears(); ii++) {
		vec_mp & curr_linear = W_critcurve.linear(ii);
		
		increase_size_vec_mp(curr_linear, W_critcurve.num_variables()); curr_linear->size = W_critcurve.num_variables();
		
		for (int jj=W_surf.num_variables(); jj<W_critcurve.num_variables(); jj++) {
			set_zero_mp(&(curr_linear)->coord[jj]);
		}
	}
	
	
	W_critcurve.only_natural_vars();// trim off the synthetic variables for the input witness set for computing critical points
	int blabla;
	parse_input_file("input_critical_curve", &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	crit_curve_.randomizer()->setup(W_critcurve.num_variables() - W_critcurve.num_patches() - 1,solve_options.PPD.num_funcs);
	
	
	
    return;
    
}


void Surface::compute_critcurve_critpts(WitnessSet & W_critcurve_crit,  // the computed value
                                      WitnessSet & W_critcurve, // input witness set
                                      int pi_ind,
                                      BertiniRealConfig & program_options, //as usual, goes everywhere.
                                      SolverConfiguration & solve_options) // wtb: a way to make this global
{
#ifdef functionentry_output
	std::cout << "surface::compute_critcurve_critpts" << std::endl;
#endif
	
	
	
	
	
	
	
	
	
	solve_options.use_gamma_trick = 0;
	
	NullspaceConfiguration ns_config; // this is set up in the nullspace call.
	
	SolverOutput solve_out;


	//get crit points of the surface.
	int blabla;
	
	parse_input_file(W_critcurve.input_filename(),&blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	int temp_degree = solve_options.T.AMP_bound_on_degree;
	std::ifstream fin("deg.out");
	int max_degree = 0;
	for (int ii=0; ii<solve_options.PPD.num_funcs; ii++) {
		fin >> temp_degree;
		if (temp_degree>max_degree) {
			max_degree = temp_degree;
		}
	}
	fin.close();
	

	
	std::cout << color::bold('m') << "\ncomputing critical points of the critical curve WRT pi_" << pi_ind <<  color::console_default() << std::endl;
	solve_options.T.AMP_bound_on_degree = max_degree;
	compute_crit_nullspace(solve_out, // the returned value
						   W_critcurve,            // input the original witness set
						   crit_curve_.randomizer(),
						   &(this->pi())+pi_ind,
						   1,  // dimension of ambient complex object
						   1,   //  target dimension to find
						   1,   // COdimension of the critical set to find.
						   program_options,
						   solve_options,
						   &ns_config);
	
	
	solve_options.T.AMP_bound_on_degree = temp_degree; //reset
	
	WitnessSet W_temp;
	solve_out.get_noninfinite_w_mult_full(W_temp);
	ns_config.clear();
	solve_out.reset();
	W_critcurve_crit.merge(W_temp,&solve_options.T);//(*&!@U#H*DB(F*&^@#*&$^(*#&YFNSD




	W_critcurve_crit.set_input_filename("input_critical_curve");
	W_critcurve_crit.sort_for_real(&solve_options.T);
	W_critcurve_crit.sort_for_unique(&solve_options.T);
	
	
	if (have_sphere()) {
		W_critcurve_crit.sort_for_inside_sphere(sphere_radius(), sphere_center());
	}
	
    
    return;
}




void Surface::compute_critical_curve(const WitnessSet & W_critcurve,
                                                   const WitnessSet & W_critpts,
                                                   VertexSet & V,
                                                   BertiniRealConfig & program_options,
                                                   SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_critical_curve" << std::endl;
#endif
	
	
    std::cout << color::bold('m') << "interslicing critical curve" << color::console_default() << std::endl;
    
	
    if (program_options.verbose_level() >= 3)
		W_critpts.print_to_screen();
    
	program_options.merge_edges(false);
	
    
	//fluff up the projection to have 0 entries for all the synthetic variables.
	vec_mp temp_proj; init_vec_mp2(temp_proj,0,1024);
	vec_cp_mp(temp_proj, this->pi(0));
	std::cout << temp_proj->size << std::endl;
	increase_size_vec_mp(temp_proj, W_critcurve.num_variables()); // nondestructive increase in size
    temp_proj->size = W_critcurve.num_variables();
    
	for ( int ii=this->num_variables(); ii<W_critcurve.num_variables(); ii++) {
		set_zero_mp(&temp_proj->coord[ii]);
	}
	
	
    
	
	
	crit_curve_.interslice(W_critcurve,
                          W_critpts,
                          &temp_proj,
                          program_options,
                          solve_options,
                          V);
	
	clear_vec_mp(temp_proj);
	
	crit_curve_.add_projection(this->pi(0));
	
	crit_curve_.set_num_variables(W_critcurve.num_variables());
	
    std::cout << color::magenta() << "done decomposing critical curve" << color::console_default() << std::endl;
	return;
}













void Surface::deflate_and_split(std::map< SingularObjectMetadata, WitnessSet > & split_sets,
											  std::map<int, WitnessSet > & higher_multiplicity_witness_sets,
											WitnessSet & points_which_needed_no_deflation,
											  BertiniRealConfig & program_options,
											  SolverConfiguration & solve_options)
{
	

	if (higher_multiplicity_witness_sets.size()==0) {
		return;
	}

	points_which_needed_no_deflation.copy_skeleton(higher_multiplicity_witness_sets.begin()->second);
	
	for (auto iter = higher_multiplicity_witness_sets.begin(); iter!=higher_multiplicity_witness_sets.end(); ++iter) {
		
		if (program_options.verbose_level()>=0) {
			std::cout << std::endl << color::magenta() << "splitting points for multiplicity " << iter->first << " singular curve" << color::console_default() << std::endl;
		}
		
		
		int num_this_multiplicity = 0;
		
		
		
		
		
		WitnessSet active_set = iter->second; // seed the loop/  this sucks because it duplicates data
		active_set.only_first_vars(this->num_variables());
		
		WitnessSet W_only_one_witness_point;
		W_only_one_witness_point.copy_skeleton(active_set);
		
		
		while (active_set.has_points()) {
			
			WitnessSet W_reject; // these will be populated in the find matching call.
			SingularObjectMetadata curr_index(iter->first,num_this_multiplicity);
			std::stringstream converter;
			
			converter << "_singcurve_mult_" << iter->first << "_" << num_this_multiplicity;
			
			boost::filesystem::path singcurve_filename = program_options.input_filename();
			singcurve_filename += converter.str(); converter.clear(); converter.str("");
			
			
			W_only_one_witness_point.reset_points();
			W_only_one_witness_point.add_point( active_set.point(0) ); // exists by entrance condition
			W_only_one_witness_point.real_threshold_points(solve_options.T.real_threshold);
			W_only_one_witness_point.write_dehomogenized_coordinates("singular_witness_points_dehomogenized"); // write the points to file
			
			
			
			int num_deflations, *deflation_sequence = NULL;
			isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
								  program_options.input_filename(), // start from the beginning.
								  "singular_witness_points_dehomogenized",
								  singcurve_filename,
								  program_options.max_deflations());
			
			
			if (num_deflations==0) {
				points_which_needed_no_deflation.add_point(active_set.point(0));
				std::cout << "found a point which did not need deflation!!!" << std::endl;
				print_point_to_screen_matlab(active_set.point(0),"anomaly");
			}
			free(deflation_sequence);
			
			active_set.set_input_filename(singcurve_filename);
			active_set.set_dimension(1);
			
			int blabla;
			parse_input_file(singcurve_filename,&blabla);
			preproc_data_clear(&solve_options.PPD); // ugh this sucks
			parse_preproc_data("preproc_data", &solve_options.PPD);
			
			std::cout << "testing points for deflation validity" << std::endl;
			
			
			find_matching_singular_witness_points(split_sets[curr_index],
												  W_reject, //W_reject contains the points of this multiplicity which DO NOT satisfy this deflation.  they must be deflated again.
												  active_set,//input witness set
												  solve_options);
			
			split_sets[curr_index].set_input_filename(singcurve_filename);
			
			
			if (W_reject.has_points()) {
				std::cout << "found that current singular witness set had " << W_reject.num_points() << " non-deflated points" << std::endl;
			}
			
			//TODO: this is an ideal place for a swap operator.
			active_set = W_reject;
			num_this_multiplicity++;
			
		}
		
		
		
		
		
		
	}
	
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		std::cout << "multiplicity " << iter->first.multiplicity() << " " << iter->first.index() << std::endl;
	}
	
}



void Surface::compute_singular_crit(WitnessSet & W_singular_crit,
												  const std::map<SingularObjectMetadata, WitnessSet> & split_sets,
												  VertexSet & V,
												  BertiniRealConfig & program_options,
												  SolverConfiguration & solve_options)
{
	
	
	
	W_singular_crit.set_num_variables(this->num_variables());
	W_singular_crit.set_num_natural_variables(this->num_variables());
	W_singular_crit.copy_patches(*this);
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		
		std::cout << std::endl << color::magenta() << "getting critical points for multiplicity " << iter->first.multiplicity() << " " << iter->first.index() << " singular curve" << color::console_default() << std::endl;
		
		
		
		
		int blabla;
		parse_input_file(iter->second.input_filename(),&blabla);
		preproc_data_clear(&solve_options.PPD); // ugh this sucks
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		std::ifstream fin("deg.out");
		int max_degree = 0;
		int temp_degree;
		for (int ii=0; ii<solve_options.PPD.num_funcs; ii++) {
			fin >> temp_degree;
			if (temp_degree>max_degree) {
				max_degree = temp_degree;
			}
		}
		fin.close();
		
		temp_degree = solve_options.T.AMP_bound_on_degree;
		solve_options.T.AMP_bound_on_degree = max_degree;
		
		
		
		
		singular_curves_[iter->first].randomizer()->setup(iter->second.num_variables()-iter->second.num_patches()-1, solve_options.PPD.num_funcs);

		
		
		NullspaceConfiguration ns_config;
		SolverOutput solve_out;
		
		
		compute_crit_nullspace(solve_out,                   // the returned value
							   iter->second,            // input the witness set with linears
							   singular_curves_[iter->first].randomizer(),
							   &(this->pi()),
							   1,                                // dimension of ambient complex object
							   1,                                //  target dimension to find
							   1,                                // COdimension of the critical set to find.
							   program_options,
							   solve_options,
							   &ns_config);
		
		
		solve_options.T.AMP_bound_on_degree = temp_degree;
		
		WitnessSet W_this_round;
		solve_out.get_noninfinite_w_mult_full(W_this_round);
		
		ns_config.clear();
		
		W_this_round.sort_for_unique(&solve_options.T); // this could be made to be unnecessary, after rewriting a bit of solverout
		W_this_round.sort_for_real(&solve_options.T);
		W_this_round.set_input_filename(iter->second.input_filename());
		
		if (have_sphere()) {
			W_this_round.sort_for_inside_sphere(sphere_radius(), sphere_center());
		}
		singular_curves_[iter->first].add_witness_set(W_this_round,Critical,V); // creates the curve Decomposition for this multiplicity
		
		W_singular_crit.merge(W_this_round,&solve_options.T);
		
	}
	
}





void Surface::compute_singular_curves(const WitnessSet & W_total_crit,
													const std::map< SingularObjectMetadata, WitnessSet> & split_sets,
													VertexSet & V,
													BertiniRealConfig & program_options,
													SolverConfiguration & solve_options)
{
	program_options.merge_edges(false);
	
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		
		singular_curves_[iter->first].set_input_filename(iter->second.input_filename());
		singular_curves_[iter->first].copy_sphere_bounds(*this);
		
		
		std::cout << "getting sphere intersection with singular curve " << iter->first.multiplicity() << " " << iter->first.index() << std::endl;
		WitnessSet W_sphere_intersection;
		// now get the sphere intersection critical points and ends of the interval
		singular_curves_[iter->first].get_sphere_intersection_pts(&W_sphere_intersection,  // the returned value
															iter->second,       // all else here is input
															program_options,
															solve_options);
		
		W_sphere_intersection.set_input_filename(iter->second.input_filename());
		
		
//		W_sphere_intersection.only_first_vars(this->num_variables); // throw out the extra variables.
		W_sphere_intersection.sort_for_real(&solve_options.T);
		W_sphere_intersection.sort_for_unique(&solve_options.T);
		
		singular_curves_[iter->first].add_witness_set(W_sphere_intersection,Critical,V); // creates the curve Decomposition for this multiplicity
		
		W_sphere_intersection.merge(W_total_crit,&solve_options.T);
		W_sphere_intersection.set_input_filename("should_have_already_been_added_elsewhere");
		
		
		
		
		
		
		std::cout << "interslicing singular curve " << std::endl;
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		singular_curves_[iter->first].interslice(iter->second,
												W_sphere_intersection,
												&(this->pi()),
												program_options,
												solve_options,
												V);
		singular_curves_[iter->first].add_projection(this->pi(0));
		
		num_singular_curves_++;
		
		
		
	}
	
}



// will compute a randomizer matrix since you don't provide one. must have current PPD in solve_options for this to work correctly
int find_matching_singular_witness_points(WitnessSet & W_match,
										  WitnessSet & W_reject,
										  const WitnessSet & W,
										  SolverConfiguration & solve_options)
{
	
	if (W.has_no_points()) {
		std::cout << color::red() << "input witness set for find_matching_  has NO points, but hypothetically it does..." << color::console_default() << std::endl;
		return TOLERABLE_FAILURE;
	}
	
	
	// assumes input file for W is already parsed.
	
	
	
	prog_t SLP;
	setupProg(&SLP, solve_options.T.Precision, 2);
	
	
	comp_mp zerotime; init_mp(zerotime);
	set_zero_mp(zerotime);
	
	
	
	eval_struct_mp ED; init_eval_struct_mp(ED, 0, 0, 0);
	
	tracker_config_t * T = &solve_options.T;
	double tol = MAX(T->final_tol_times_mult, T->sing_val_zero_tol);
	
	mat_mp U, E, V; init_mat_mp(U, 0, 0); init_mat_mp(E, 0, 0); init_mat_mp(V, 0, 0);
	
	
	
	
	evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, W.point(0), zerotime, &SLP);
	int hypothesis_corank = svd_jacobi_mp_prec(U, E, V, ED.Jv, tol, T->Precision); // this wraps around svd_jacobi_mp.
	
	
	
	std::vector< bool > validity_flag;
	validity_flag.push_back(true);
	for (unsigned int zz = 1;zz<W.num_points();++zz) // by hypothesis, the first (0th) point satisfies the deflation
	{
		
		vec_mp & curr_point = W.point(zz);
		
		evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, curr_point, zerotime, &SLP);
		
		// first, check that the point satifies the system.
		if (d_vec_abs_mp(ED.funcVals)>T->final_tol_times_mult) { // is this the correct measure of the vector to compare?
			validity_flag.push_back(false);
			continue;
		}
		
		
		// now, check the rank and make sure is same as first (0th) point.
		int corank = svd_jacobi_mp_prec(U, E, V, ED.Jv, tol, T->Precision); // this wraps around svd_jacobi_mp.
		
		if (corank != hypothesis_corank) {
			
			validity_flag.push_back(false);
			continue;
		}
		else{
			validity_flag.push_back(true);
		}
	}
	
	
	for (unsigned int zz=0; zz<W.num_points(); zz++) {
		if (validity_flag[zz]==true) { // trivially true for first point -- it generated the deflation!
			W_match.add_point(W.point(zz));
		}
		else
		{
			std::cout << "adding reject point" << std::endl;
			W_reject.add_point(W.point(zz));
		}
	}
	
	
	
	W_match.copy_skeleton(W);
	W_reject.copy_skeleton(W);
	
	W_match.copy_linears(W);
	W_reject.copy_linears(W);
	
	W_match.copy_patches(W);
	W_reject.copy_patches(W);
	
	
	
	clear_mat_mp(U); clear_mat_mp(E); clear_mat_mp(V);
	clear_mp(zerotime);
	clear_eval_struct_mp(ED); clearProg(&SLP, solve_options.T.MPType, 1);
	
	
	return SUCCESSFUL;
}







void Surface::compute_sphere_witness_set(const WitnessSet & W_surf,
                                                       WitnessSet & W_intersection_sphere,
                                                       BertiniRealConfig & program_options,
                                                       SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_sphere_witness_set" << std::endl;
#endif
    
	
	if (program_options.verbose_level()>=0) {
		std::cout << color::magenta() << "getting sphere witness set" << color::console_default() << std::endl;
	}
	
	//build up the start system
	solve_options.robust = true;
	solve_options.use_gamma_trick = 0;
	
	int blabla;
	
	
	
	parse_input_file(W_surf.input_filename(), &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
    
    
	if (!this->randomizer()->is_ready()) {
		throw std::logic_error("randomizer not ready");
	}
	
	
	sphere_curve_.randomizer()->setup(W_surf.num_variables()-W_surf.num_patches()-W_surf.num_linears(), solve_options.PPD.num_funcs);
	
	
	
	MultilinConfiguration ml_config(solve_options,this->randomizer()); // copies in the randomizer matrix and sets up the SLP & globals.
	
    
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
    
	init_vec_mp2(multilin_linears[0],W_surf.num_variables(),solve_options.T.AMP_max_prec);
	multilin_linears[0]->size = W_surf.num_variables();
    
	init_vec_mp2(multilin_linears[1],0,solve_options.T.AMP_max_prec);
	vec_cp_mp(multilin_linears[1],W_surf.linear(1));
	
    
    
	WitnessSet W_sphere;
	SphereConfiguration sp_config(sphere_curve_.randomizer());
	
	
	for (int ii=0; ii<2; ii++) {
		
		for (int jj=0; jj<W_surf.num_variables(); jj++) {
			get_comp_rand_mp(&multilin_linears[0]->coord[jj]);
		}
		
		vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[0]); // copy in the first multilinear to the new witness set we are computing.
		
		WitnessSet W_temp;
		
		
		
		
		
		SolverOutput fillme;
		multilin_solver_master_entry_point(W_surf,         // WitnessSet
                                           fillme, // the new data is put here!
                                           multilin_linears,
                                           ml_config,
                                           solve_options);
		
		//get stuff from fillme into W_temp, or whatever
		fillme.get_noninfinite_w_mult(W_temp);
		
		W_sphere.merge(W_temp,&solve_options.T);
		
	}
	
	clear_vec_mp(multilin_linears[0]);
	clear_vec_mp(multilin_linears[1]);
	free(multilin_linears);
	
	
	
	W_sphere.add_linear(W_surf.linear(1));
	W_sphere.copy_patches(W_surf); //.patch_mp[0]
	
	
	
	// need to actually move to the sphere system now.
	
	
	sp_config.set_memory(solve_options);
	sp_config.set_center(this->sphere_center());
	sp_config.set_radius(this->sphere_radius());
	
	
	
	
	
	
	
	
	SolverOutput fillme;
	sphere_solver_master_entry_point(W_sphere,
									 fillme,
                                     sp_config,
                                     solve_options);
	
	fillme.get_noninfinite_w_mult_full(W_intersection_sphere);
	
	W_intersection_sphere.set_input_filename("input_surf_sphere");
	
}




void Surface::compute_sphere_crit(const WitnessSet & W_intersection_sphere,
                                                WitnessSet & W_sphere_crit,
                                                BertiniRealConfig & program_options,
                                                SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_sphere_crit" << std::endl;
#endif
	
	
	int blabla;
	parse_input_file("input_surf_sphere", &blabla); // having already been written to disk
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
    
    
    
	sphere_curve_.randomizer()->setup(W_intersection_sphere.num_variables()-W_intersection_sphere.num_patches()-1, solve_options.PPD.num_funcs);
    
	
	
	solve_options.use_gamma_trick = 0;
	
	NullspaceConfiguration ns_config;
	
	std::cout << color::magenta() << "computing critical points of sphere curve" << color::console_default() << std::endl;
	
	SolverOutput solve_out;
	
	compute_crit_nullspace(solve_out,                   // the returned value
                           W_intersection_sphere,            // input the witness set with linears
                           sphere_curve_.randomizer(),
                           &(this->pi(0)),
                           1,                                // dimension of ambient complex object
                           1,                                //  target dimension to find
                           1,                                // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	
	solve_out.get_noninfinite_w_mult_full(W_sphere_crit);
	
	ns_config.clear();
	
	W_sphere_crit.sort_for_unique(&solve_options.T);
	W_sphere_crit.sort_for_real(&solve_options.T);
    
	W_sphere_crit.set_input_filename("input_surf_sphere");
	
	return;
}



void Surface::compute_bounding_sphere(const WitnessSet & W_intersection_sphere,
                                                    const WitnessSet & W_crit,
                                                    VertexSet & V,
                                                    BertiniRealConfig & program_options,
                                                    SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_bounding_sphere" << std::endl;
#endif
    
	program_options.merge_edges(false);
	
	comp_mp temp; init_mp(temp);
	set_zero_mp(temp);
	mpf_set_d(temp->r, 1.1);
	mul_mp(temp, temp, this->sphere_radius())
	this->sphere_curve_.set_sphere_radius(temp);
	this->sphere_curve_.set_sphere_center(this->sphere_center());
	
	std::cout << color::magenta() << "entering interslice for sphere" << color::console_default() << std::endl;
	// then feed it to the interslice algorithm
	this->sphere_curve_.interslice(W_intersection_sphere, // the witness set with a single linear and some patches.
                                  W_crit, // the critical points
                                  &(pi(0)),
                                  program_options,
                                  solve_options,
                                  V);
	
	sphere_curve_.set_input_filename(W_intersection_sphere.input_filename());
	this->sphere_curve_.add_projection(pi(0));
	this->sphere_curve_.set_num_variables(this->num_variables());
	
	
	std::cout << color::magenta() << "done decomposing sphere curve" << color::console_default() << std::endl;
	return;
}















void Surface::compute_slices(const WitnessSet W_surf,
                                           VertexSet & V,
                                           vec_mp projection_values_downstairs,
                                           std::vector< Curve > & slices,
                                           BertiniRealConfig & program_options,
                                           SolverConfiguration & solve_options,
                                           std::string kindofslice)
{
#ifdef functionentry_output
	std::cout << "surface::compute_slices" << std::endl;
#endif
	
	
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
	
	
	init_vec_mp2(multilin_linears[0], W_surf.num_variables(),1024);
	init_vec_mp2(multilin_linears[1], W_surf.num_variables(),1024);
	vec_cp_mp(multilin_linears[0],pi(0));
	vec_cp_mp(multilin_linears[1],W_surf.linear(1));
	
	
	
	
	
	slices.resize(projection_values_downstairs->size);
	
	int blabla;
	parse_input_file(W_surf.input_filename(), & blabla);
	
	
	
	MultilinConfiguration ml_config(solve_options); // copies in the randomizer matrix and sets up the SLP & globals.
	
	for (int ii=0; ii<projection_values_downstairs->size; ii++){
		
		if (program_options.verbose_level()>=0) {
			std::cout << color::magenta() << "decomposing " << kindofslice << " slice " << ii << " of " << projection_values_downstairs->size << color::console_default() << std::endl;
			print_comp_matlab(&projection_values_downstairs->coord[ii], "target_proj");
		}
		
		
		solve_options.backup_tracker_config("surface_slice");
		
		
		WitnessSet slice_witness_set; // deliberately scoped variable
		
		
		
		std::stringstream converter;
		converter << ii;
		
		int blabla;
		parse_input_file(W_surf.input_filename(), &blabla);
		
		preproc_data_clear(&solve_options.PPD); // ugh this sucks
		parse_preproc_data("preproc_data", &solve_options.PPD);
		

		ml_config.set_randomizer(this->randomizer());
		neg_mp(&multilin_linears[0]->coord[0], &projection_values_downstairs->coord[ii]);
		
		
		if (program_options.verbose_level()>=1) {
			std::cout << color::green() << "getting slice witness points and linear" << color::console_default() << std::endl;
		}
		

		
		if (program_options.quick_run()<=1)
			solve_options.robust = true;
		else
			solve_options.robust = false;
		
		
		
		
		SolverOutput fillme;
		multilin_solver_master_entry_point(W_surf,         // WitnessSet
										   fillme, // the new data is put here!
										   multilin_linears,
										   ml_config,
										   solve_options);
		

		
		fillme.get_noninfinite_w_mult(slice_witness_set);
		fillme.reset();
		
		boost::filesystem::path slicename = W_surf.input_filename();
		slicename += "_"; slicename += kindofslice; slicename += "slice_"; slicename += converter.str();
		create_sliced_system(W_surf.input_filename(), slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		

		
		
		slice_witness_set.set_num_variables(W_surf.num_variables());
		slice_witness_set.set_input_filename(slicename);
		slice_witness_set.add_linear(multilin_linears[1]);
		slice_witness_set.add_patch(W_surf.patch(0));
		slice_witness_set.copy_names(W_surf);
		slice_witness_set.set_dimension(1);
		
		
		
		slices[ii].set_input_filename(slicename);
		slices[ii].copy_sphere_bounds(*this);
		
		
		
		
		// the memory for the multilin system will get erased in this call...
		

		
		
		if (program_options.verbose_level()>=1) {
			std::cout << color::green() << "computing slice" << color::console_default() << std::endl;
		}
		
		
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		slices[ii].computeCurveSelfConj(slice_witness_set,
									   &(pi(1)),
									   V,
									   program_options, solve_options);

		

		
		slice_witness_set.reset();
		solve_options.restore_tracker_config("surface_slice");
		
		slices[ii].add_projection(pi(1));
		slices[ii].set_num_variables(this->num_variables());
		
		
		
        // does it matter speedwise whether i do this before or after the copy immediately above?  i think the answer is no.
        V.assert_projection_value(slices[ii].all_edge_indices(), &projection_values_downstairs->coord[ii], 0); // the 0 is an index into the number of projections.
		
		this->output_main(program_options.output_dir());
		V.print(program_options.output_dir()/ "V.vertex");
		
		
		if (program_options.verbose_level()>=0) {
			std::cout << color::magenta() << "DONE decomposing " << kindofslice << "slice " << ii << color::console_default() << std::endl;
		}
		
		
	} // re: for loop
	
	
	clear_vec_mp(multilin_linears[0]);
	clear_vec_mp(multilin_linears[1]);
	free(multilin_linears);
	return;
}




















void Surface::connect_the_dots(VertexSet & V,
                                             BertiniRealConfig & program_options,
                                             SolverConfiguration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::connect_the_dots" << std::endl;
#endif
	
	std::cout << color::bold('m') << "***\n\nCONNECT THE DOTS\n\n***" << color::console_default() << std::endl;
    
	
	MidpointConfiguration md_config;
	md_config.setup(*this, solve_options); // yep, pass 'this' object into another call. brilliant.
	
	
	if (solve_options.use_parallel())
		master_connect(V, md_config, solve_options, program_options);
	else
		serial_connect(V, md_config, solve_options, program_options);
	
    
	
	
	return;
}


void Surface::serial_connect(VertexSet & V, MidpointConfiguration & md_config, SolverConfiguration & solve_options, BertiniRealConfig & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::serial_connect" << std::endl;
#endif
	
    
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	
	for (unsigned int ii=0; ii!=mid_slices_.size(); ii++) { // each edge of each midslice will become a Face.  degenerate edge => degenerate Face.
		
		
		for (unsigned int jj=0; jj<mid_slices_[ii].num_edges(); jj++) {
			
			
			
			//make Face
			
			Face F = make_face(ii, jj, V, md_config, solve_options, program_options);
			

			
			
			if (!F.is_degenerate())
			{
				add_face(F);
				
//					boost::filesystem::path::rename(program_options.output_dir() / "F.faces");
//					this->print_faces(program_options.output_dir() + "_bak" / "F.faces");
				this->print_faces(program_options.output_dir() / "F.faces");
//					this->output_main(program_options.output_dir());
//					V.print(program_options.output_dir()/ "V.vertex");
			}
			

		}
	}
	
	return;
}



void Surface::master_connect(VertexSet & V, MidpointConfiguration & md_config, SolverConfiguration & solve_options, BertiniRealConfig & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::master_connect" << std::endl;
#endif
	
	
    
	MPI_Status statty_mc_gatty;
	
	solve_options.call_for_help(MIDPOINT_SOLVER); // sets available workers, too
	
	MPI_Barrier(solve_options.comm());
	
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	MPI_Barrier(solve_options.comm());
	
	
	md_config.bcast_send(solve_options);
	
	
	
	if (V.num_vertices() > 1e5) {
		std::cout << color::red() << "attempting to transmit over 1e5 points to all workers..." << color::console_default() << std::endl;
	}
	
	
	MPI_Barrier(solve_options.comm());
	
	
	//seed the workers
	for (int ii=1; ii<solve_options.num_procs(); ii++) {
		this->send(ii, solve_options);
		V.send(ii, solve_options);
	}
	
	
	this->output_main(program_options.output_dir());
	V.print(program_options.output_dir()/ "V.vertex");
	
	// this loop is semi-self-seeding
	for (unsigned int ii=0; ii!=mid_slices_.size(); ii++) { // each edge of each midslice will become a Face.  degenerate edge => degenerate Face.
		
		
		for (unsigned int jj=0; jj<mid_slices_[ii].num_edges(); jj++) {
			
			if (mid_slices_[ii].get_edge(jj).is_degenerate()) {
				continue;
			}
			
			int next_worker = solve_options.activate_next_worker();
			
			int send_num_faces = 1;// num_faces doubles as the keep_going signal.  if 0, the worker halts.
			MPI_Send(&send_num_faces, 1, MPI_INT, next_worker, NUMPACKETS, solve_options.comm());
			
			
			master_face_requester(ii,jj, next_worker, solve_options);
			
			
			
			if (solve_options.have_available()) {
				
				continue;
			}
			else
			{
				//perform blocking receive of the Face.
				int recv_num_faces;
				MPI_Recv(&recv_num_faces, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, solve_options.comm(), &statty_mc_gatty);
				bool added_face = false;
				for (int qq = 0; qq<recv_num_faces; qq++) {
					Face F;
					F.receive(statty_mc_gatty.MPI_SOURCE, solve_options);
					if (!F.is_degenerate()) {
						add_face(F);
						added_face = true;
					}
					
				}
				
				solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
				
				
				
				//TODO: this needs to be inside of a better protector, in the sense that we shouldn't do it after EVERY Face, but rather at thoughtful times.
				if (added_face) { // only the faces file really needs to be updated.  the rest stay the same.  this is very wasteful.
//					boost::filesystem::path::rename(program_options.output_dir() / "F.faces");
//					this->print_faces(program_options.output_dir() + "_bak" / "F.faces");
					this->print_faces(program_options.output_dir() / "F.faces");
//					this->output_main(program_options.output_dir());
//					V.print(program_options.output_dir()/ "V.vertex");
				}
				
				
			}
		}
	}
	
	
	//cleanup
	
	
	while (solve_options.have_active()) {
		int recv_num_faces;
		MPI_Recv(&recv_num_faces, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, solve_options.comm(), &statty_mc_gatty);
		for (int ii=0; ii<recv_num_faces; ii++) {
			Face F;
			F.receive(statty_mc_gatty.MPI_SOURCE, solve_options);
			add_face(F);
			solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
		}
	}
	
	solve_options.send_all_available(0);
	
	
	return;
}





void Surface::worker_connect(SolverConfiguration & solve_options, BertiniRealConfig & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::worker_connect" << std::endl;
#endif
	
	MPI_Barrier(solve_options.comm());
	
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	MPI_Barrier(solve_options.comm());
	solve_options.robust = true;
	

	MidpointConfiguration md_config;
	md_config.bcast_receive(solve_options);
	//receive the md_config from the master.  it holds the three SLP's, as well as everything needed to run the system except:
	//	• system types
	//	• patches
	//	• starting point
	//which are all updated from another call later.
	
	
	MPI_Barrier(solve_options.comm());
	
	this->receive(solve_options.head(), solve_options);
	
	
	
	VertexSet V;
	V.set_tracker_config(&solve_options.T);
	
	V.receive(solve_options.head(), solve_options);
	
	
	
	
	
	
	while (1) {
		
		// get the continue or discontinue signal.
		int keep_going;
		MPI_Status statty_mc_gatty;
		MPI_Recv(&keep_going, 1, MPI_INT, solve_options.head(), NUMPACKETS, solve_options.comm(), &statty_mc_gatty);
		if (keep_going==0) {
			break;
		}
		
		
		
		int ii, jj;
		// get the indices of the Face to make.
		worker_face_requester(ii,jj,solve_options);
		
		
		//make the Face
		Face F = make_face(ii, jj, V, md_config, solve_options, program_options);
		
		
		//send the Face back to master.
		int send_num_faces = 1;
		MPI_Send(&send_num_faces, 1, MPI_INT, solve_options.head(), NUMPACKETS, solve_options.comm());
		
		F.send(solve_options.head(), solve_options);
	}
	
	
	
	return;
}




void Surface::master_face_requester(int ii, int jj, int next_worker, ParallelismConfig & mpi_config) const
{
#ifdef functionentry_output
	std::cout << "surface::master_face_requester" << std::endl;
#endif
	
	
	
	
	int * buffer = new int[2];
	buffer[0] = ii;
	buffer[1] = jj;
	MPI_Send(&buffer[0], 2, MPI_INT, next_worker, DATA_TRANSMISSION, mpi_config.comm());
	delete [] buffer;
}


void Surface::worker_face_requester(int & ii, int & jj, ParallelismConfig & mpi_config) const
{
#ifdef functionentry_output
	std::cout << "surface::worker_face_requester" << std::endl;
#endif
	
	
	
	
	int * buffer = new int[2];
	MPI_Status statty_mc_gatty;
	
	MPI_Recv(&buffer[0], 2, MPI_INT, mpi_config.head(), DATA_TRANSMISSION, mpi_config.comm(), &statty_mc_gatty);
	ii = buffer[0];
	jj = buffer[1];
	
	delete [] buffer;
	return;
}

Face Surface::make_face(int ii, int jj, VertexSet & V,
									  MidpointConfiguration & md_config,
									  SolverConfiguration & solve_options, BertiniRealConfig & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::make_face" << std::endl;
#endif
	
	
	
	Curve & current_midslice = mid_slices_[ii];
	Curve & left_critslice = crit_slices_[ii];
	Curve & right_critslice = crit_slices_[ii+1];
	
	// assert some solver options
	solve_options.use_gamma_trick = 0;
	
	if (program_options.quick_run()<=1)
		solve_options.robust = true;
	else
		solve_options.robust = false;
    
	
	//create the Face
	Face F;
	

	
	if (current_midslice.get_edge(jj).is_degenerate()) {
		return F;
	}
	
	
	
	
	
	
	
	
	
	
	
	F.crit_slice_index(ii); // the index of which midslice this Face came from.
	F.midpt( current_midslice.get_edge(jj).midpt() ); // index the point
	
	
	std::cout << color::magenta() << "\n\n\n*****************************\nmidslice " << ii << " / " << this->mid_slices_.size() <<  ", edge " << jj << " / " << current_midslice.num_edges() << color::console_default() << "\n";
	
	
	std::cout << color::brown() << "current midpoint: " <<  current_midslice.get_edge(jj).midpt()  << " " << color::console_default() << "\n";
	
	
	
	// perform a search to find the top and bottom edges in the appropriate curve.
	
	
	// get the type of system for the top and bottom edges.  this is determined by reading the system name for the midpoints.
	// info on the files from which the points came from.
	int bottom_input_index = V[current_midslice.get_edge(jj).left()].input_filename_index();
	int top_input_index = V[current_midslice.get_edge(jj).right()].input_filename_index();
	
	bool bail_out = false;
	if (md_config.systems.find(V.filename(bottom_input_index).string())==md_config.systems.end()) {
		std::cout << "bottom system is " << V.filename(bottom_input_index) << ", which is not in md_config" << std::endl;
		bail_out = true;
	}
	if (md_config.systems.find(V.filename(top_input_index).string())==md_config.systems.end()) {
		std::cout << "top system is " << V.filename(top_input_index) << ", which is not in md_config" << std::endl;
		bail_out = true;
	}
	
	
	md_config.system_name_bottom = V.filename(bottom_input_index).filename().string();
	md_config.system_name_top = V.filename(top_input_index).filename().string();
	md_config.system_name_mid = this->input_filename().filename().string();
	
	
	if (program_options.verbose_level()>=0) {
		std::cout << md_config.system_name_top << " " << md_config.system_name_bottom << std::endl;
	}
	
	
	
	const Curve * top_curve = curve_with_name(md_config.system_name_top);
	const Curve * bottom_curve = curve_with_name(md_config.system_name_bottom);
	
	//man, i hate checking for null...
	
	if (top_curve==NULL) {
		std::cout << "did not find matching top curve: " << md_config.system_name_top << std::endl;
		bail_out = true;
	}
	
	if (bottom_curve==NULL) {
		std::cout << "did not find matching bottom curve: " << md_config.system_name_bottom << std::endl;
		bail_out = true;
	}
	
	
	
	if (left_critslice.num_edges() == 0) {
		std::cout << "critslice " << ii << " has no edges!" << std::endl;
		bail_out = true;
	}
	if (right_critslice.num_edges() == 0) {
		std::cout << "critslice " << ii+1 << " has no edges!" << std::endl;
		bail_out = true;
	}
	
	
	// get the bottom and top edges for this Face.
	
	int num_bottom_vars = md_config.num_bottom_vars();
	int num_top_vars = md_config.num_top_vars();
	
	F.set_bottom_edge( bottom_curve->edge_w_midpt(current_midslice.get_edge(jj).left()) ); // index the *edge*
	F.system_name_bottom( md_config.system_name_bottom );
	
	F.set_top_edge( top_curve->edge_w_midpt(current_midslice.get_edge(jj).right()) ); // index the *edge*
	F.system_name_top( md_config.system_name_top );
	
	if (F.bottom_edge() < 0) { // this would happen because the bottom point was of type PROBLEMATIC
		std::cout << "bottom edge is set to negative" << std::endl;
		bail_out = true;
	}
	
	if (F.top_edge() < 0) { // this would happen because the bottom point was of type PROBLEMATIC
		std::cout << "top edge is set to negative" << std::endl;
		bail_out = true;
	}
		
	
	if (num_bottom_vars==0) {
		std::cout << "0 bottom variables" << std::endl;
		bail_out = true;
	}
	
	if (num_top_vars==0) {
		std::cout << "0 top variables" << std::endl;
		bail_out = true;
	}
	
	
	
	if (bail_out) {
		std::cout << color::red() << "bailing out " << ii << " " << jj << "." << std::endl;
		
		std::cout << "tracking from these point indices:" << std::endl;
		std::cout <<  current_midslice.get_edge(jj).left()  << " " << current_midslice.get_edge(jj).midpt()  << " "  << current_midslice.get_edge(jj).right() << color::console_default() << std::endl;
		
		return F;
	}
	
	
	
	
	
	comp_mp temp, temp2, temp3; init_mp2(temp,1024); init_mp2(temp2,1024); init_mp2(temp3,1024);
	comp_mp numer, denom; init_mp2(numer,1024); init_mp2(denom,1024);
	comp_mp proj_top, proj_bottom, proj_mid; init_mp2(proj_mid,1024); init_mp2(proj_bottom,1024); init_mp2(proj_top,1024);
	vec_mp found_point; init_vec_mp2(found_point, this->num_variables(),1024); found_point->size = this->num_variables();
	
	
	//copy in the start point as three points concatenated.
	WitnessSet W_midtrack;
	vec_mp blank_point;  init_vec_mp2(blank_point, 0,1024);
	W_midtrack.add_point(blank_point);
	clear_vec_mp(blank_point);
	
	W_midtrack.set_num_variables(this->num_variables() + num_bottom_vars + num_top_vars);
	W_midtrack.set_num_natural_variables(this->num_variables());
	change_size_vec_mp( W_midtrack.point(0), W_midtrack.num_variables());
	W_midtrack.point(0)->size = W_midtrack.num_variables(); // destructive resize
	
	
	// mid
	int var_counter = 0;
	for (int kk=0; kk<this->num_variables(); kk++) {
		set_mp(&(W_midtrack.point(0)->coord[kk]), &(V[current_midslice.get_edge(jj).midpt()].point())->coord[kk]);
		var_counter++;
	}
	
	// bottom
	int offset = var_counter;
	for (int kk=0; kk<num_bottom_vars; kk++) {
		set_mp(&(W_midtrack.point(0)->coord[kk+offset]), &(V[current_midslice.get_edge(jj).left()].point())->coord[kk]); // y0
		var_counter++;
	}
	
	// top
	offset = var_counter;
	for (int kk=0; kk<num_top_vars; kk++) {
		set_mp(&(W_midtrack.point(0)->coord[kk+offset]), &(V[current_midslice.get_edge(jj).right()].point())->coord[kk]); // y2
		var_counter++;
	}
	
	
	
	
	
	
	
	//copy in the patches appropriate for the systems we will be tracking on.  this could be improved.
	W_midtrack.reset_patches();
	
	for (unsigned int qq = 0; qq<this->num_patches(); qq++) {
		W_midtrack.add_patch(this->patch(qq));
	}
	
	for (unsigned int qq = 0; qq<bottom_curve->num_patches(); qq++) {
		W_midtrack.add_patch(bottom_curve->patch(qq));
	}
	
	for (unsigned int qq = 0; qq<top_curve->num_patches(); qq++) {
		W_midtrack.add_patch(top_curve->patch(qq));
	}
	
	

	
	
	
	
	
	
	// make u, v target values.
	
	

	set_mp(md_config.crit_val_left,   &(V[ left_critslice.get_edge(0).midpt() ].projection_values())->coord[0]);
	set_mp(md_config.crit_val_right,  &(V[ right_critslice.get_edge(0).midpt() ].projection_values())->coord[0]);
	
	
	// the u direction corresponds to pi[0].
	for (int zz=0; zz<2; zz++) { // go left (zz=0) and right (zz=1)
		if (zz==0) {
			std::cout << "\n      <<=========   going left" << std::endl;
			std::cout << "tracking from these point indices:" << std::endl;
			std::cout <<  current_midslice.get_edge(jj).left()  << " " << current_midslice.get_edge(jj).midpt()  << " "  << current_midslice.get_edge(jj).right() << std::endl;
		}
		else{
			std::cout << "\n\n           going right   =======>> " << std::endl;
			std::cout << "tracking from these point indices:" << std::endl;
			std::cout <<  current_midslice.get_edge(jj).left()  << " " << current_midslice.get_edge(jj).midpt()  << " "  << current_midslice.get_edge(jj).right() << std::endl;
		}
		
		
		
		//track
		int final_top_ind = -1, final_bottom_ind = -2; // indexes in V of the bottom and top points of the left or right edge.
		
		
		if (zz==0) {
			
			set_zero_mp(md_config.u_target);
			final_bottom_ind = bottom_curve->get_edge(F.bottom_edge()).left();
			final_top_ind = top_curve->get_edge(F.top_edge()).left();
			
			
		}
		else{ // zz==1, and going right
			
			set_one_mp(md_config.u_target);
			final_bottom_ind = bottom_curve->get_edge(F.bottom_edge()).right();
			final_top_ind = top_curve->get_edge(F.top_edge()).right();
			
		}
		
		Curve & target_critslice = crit_slices_[ii+zz];
		
		
		// check for degeneracy
		if (final_bottom_ind==final_top_ind) {
			// can simply set the top or bottom edge to be this one.  know it goes there.
			std::cout << "current edge is degenerate, " << final_bottom_ind <<"=="<<final_top_ind << std::endl;
			
			
			int current_edge = target_critslice.edge_w_midpt(final_bottom_ind);
			
			if (current_edge<0) {
				std::cout << "unable to find a degenerate edge in crit_slices[" << ii+zz << "] with midpoint " << final_bottom_ind << std::endl;
//				std::cout << "making new degenerate edge" << std::endl;
//				edge E(final_bottom_ind,final_bottom_ind,final_bottom_ind);
//				current_edge = target_critslice.add_edge(E);
				continue;
			}
			
			
			if (zz==0){
				F.add_left_edge(current_edge);
			}
			else {
				F.add_right_edge(current_edge);
			}
			continue; // go to next zz value, or next midslice edge, or whatever
		}
		
		std::cout << "final top: " << final_top_ind << ", final bottom:	" << final_bottom_ind << std::endl;
		
		// get the projection values of top and bottom final points.
		set_mp(proj_top, &(V[ final_top_ind ].projection_values())->coord[1]);
		set_mp(proj_bottom, &(V[ final_bottom_ind ].projection_values())->coord[1]);
		
		
		
		// i think the projection values have already been thresholded
		real_threshold(proj_top,solve_options.T.real_threshold);
		real_threshold(proj_bottom,solve_options.T.real_threshold);
		
		
		//initialize current index trackers.
		int current_bottom_ind = final_bottom_ind;
		int current_top_ind = -12131; // initialize to impossible value;
		
		
		
		// candidates
		std::set< int > found_edges;
		std::set< int > possible_edges;
		std::set< int > remaining_possible_edges;
		
		
		for (unsigned int rr = 0; rr< target_critslice.num_edges(); rr++)
		{
			possible_edges.insert(rr);
			remaining_possible_edges.insert(rr);
		}
		
		
		while ((current_top_ind != final_top_ind) && (remaining_possible_edges.size()>0)) // int yy=0; yy<target_critslice.num_edges; yy++
		{
			
			std::cout << "target bottom: " << final_bottom_ind << " current bottom: " << current_bottom_ind << " current top: " << current_top_ind << " final top: " << final_top_ind << std::endl;
			
			
			
			std::vector< int > candidates; // indices of candidates for next one.
			
			

			
			int candidate_counter = 0;
			for (std::set<int>::iterator poss_iter=remaining_possible_edges.begin(); poss_iter != remaining_possible_edges.end(); poss_iter++) {
				
				int qq = *poss_iter;
				
				bool matches_end = ((target_critslice.get_edge(qq).left() == current_bottom_ind) || (target_critslice.get_edge(qq).right() == current_bottom_ind));
				bool already_found = (found_edges.find(qq)!=found_edges.end());
				
				
				// we gotta be moving from lower to higher...  so temp > temp2 is required
				bool correct_interval = false;
				if (matches_end) {
					set_mp(temp , &(V[ target_critslice.get_edge(qq).midpt()].projection_values())->coord[1]);
					set_mp(temp2, &(V[ final_bottom_ind].projection_values())->coord[1]);
					set_mp(temp3, &(V[ final_top_ind].projection_values())->coord[1]);
					correct_interval =  ( mpf_get_d(temp3->r) > mpf_get_d(temp->r)) && (mpf_get_d(temp->r) > mpf_get_d(temp2->r)) ;
				}
				
				// we gotta be moving from lower to higher...  so temp > temp2 is required
				bool within_bounds = false;
				if (matches_end) {
					set_mp(temp , &(V[ target_critslice.get_edge(qq).right()].projection_values())->coord[1]);
					set_mp(temp2, &(V[ final_bottom_ind].projection_values())->coord[1]);
					set_mp(temp3, &(V[ final_top_ind].projection_values())->coord[1]);
					within_bounds =  ( mpf_get_d(temp3->r) >= mpf_get_d(temp->r)) && (mpf_get_d(temp->r) > mpf_get_d(temp2->r)) ;
				}

				bool degenerate = target_critslice.get_edge(qq).is_degenerate();
				
				if ( (!already_found) && matches_end && correct_interval && (!degenerate) && within_bounds ) { //
					candidates.push_back(qq);
					
					if (program_options.verbose_level()>=1) {
						std::cout << "candidate [" << candidate_counter << "] = " <<
						target_critslice.get_edge(qq).left() << " " << target_critslice.get_edge(qq).midpt() << " " << target_critslice.get_edge(qq).right() <<  std::endl;
					}
					
					
					candidate_counter++;
				}
				else
				{
					//                            if (!correct_interval) {
					//                                print_comp_matlab(temp3,"final_top_proj_1");
					//                                print_comp_matlab(temp ,"critslices[].proj_val_1");
					//                                print_comp_matlab(temp2,"final_bottom_proj_1");
					//                            }
					if (program_options.verbose_level()>=4) {
						std::cout << "edge " << qq << " excluded: " << correct_interval << " correct_interval, " << matches_end << " matches_end, " << already_found << "  already_found" << std::endl;
					}
					
					
				}
				
			}

			
			
			
			
			if (candidate_counter==0) {
				std::cout << color::red() << "found 0 candidates for bottom index " << current_bottom_ind << color::console_default() << std::endl;
				break; // out of the while loop
			}
			
			for (int qq=0; qq<candidate_counter; qq++)
			{
				int current_edge = candidates[qq];
				
				if (program_options.verbose_level()>=0) {
					//					std::cout << "Face #: " << this->num_faces << ", zz: " << zz << ", current_edge: " << current_edge << std::endl;
					std::cout << "tracking to these indices: " << final_bottom_ind << " " << target_critslice.get_edge(current_edge).midpt() << " " << final_top_ind << std::endl;
				}
				
				
				
				
				
				set_mp(proj_mid, &(V[ target_critslice.get_edge(current_edge).midpt() ].projection_values())->coord[1]);
				
				
				

					real_threshold(proj_mid,solve_options.T.real_threshold);

				
				
				
				sub_mp(denom, proj_top, proj_bottom); // p2(w2) - p2(w0);
				sub_mp(numer, proj_mid, proj_bottom); // p2(e.w) - p2(w0);
				div_mp(md_config.v_target, numer, denom); // [p2(e.w) - p2(w0)] / [p2(w2) - p2(w0)]
				

				real_threshold(md_config.v_target,solve_options.T.real_threshold);

				
				
				
				if (solve_options.verbose_level() >= 0)
				{

					
					print_comp_matlab(md_config.u_target,"u_target");
					print_comp_matlab(md_config.v_target,"v_target");
					

				}
				
				
				if (solve_options.verbose_level() >= 1){
					print_comp_matlab(proj_top,"upper");
					print_comp_matlab(proj_bottom,"lower");
					print_comp_matlab(proj_mid,"mid");
					
					print_comp_matlab(numer,"numer");
					print_comp_matlab(denom,"denom");
					
					print_comp_matlab(md_config.crit_val_left,"proj_val_left");
					print_comp_matlab(md_config.crit_val_right,"proj_val_right");
					
					
					set_one_mp(temp);
					sub_mp(temp2, temp, md_config.v_target);
					mul_mp(temp, temp2, proj_bottom);
					
					mul_mp(temp2, md_config.v_target, proj_top);
					
					add_mp(temp3, temp, temp2);
					
					print_comp_matlab(temp3,"proj_1_target_mid");
					
				}
				
				
				SolverOutput fillme;
				WitnessSet W_new;
				midpoint_solver_master_entry_point(W_midtrack, // carries with it the start points, and the linears.
												   fillme, // new data goes in here
												   md_config,
												   solve_options);
				
				fillme.get_noninfinite_w_mult_full(W_new);
				
				
				
				//print some display to screen
				if (solve_options.verbose_level() >= 3)
				{
					print_point_to_screen_matlab((V[current_midslice.get_edge(jj).right()].point()),"top_start");
					print_point_to_screen_matlab((V[current_midslice.get_edge(jj).left()].point()),"bottom_start");
					print_point_to_screen_matlab((V[target_critslice.get_edge(current_edge).midpt() ].point()),"midpoint_target");
				}
				
				
				
				
				
				// should get a single point back from this solver.
				
				if (W_new.num_points()==0) {
					std::cout << color::red() << "midpoint tracker did not return any points :(" << color::console_default() << std::endl;
					remaining_possible_edges.erase(current_edge);
					continue;
				}
				
				
				// get only the midpoint coordinates out of the returned point
				for (int tt = 0; tt<this->num_variables(); tt++) {
					set_mp(&found_point->coord[tt], & W_new.point(0)->coord[tt]);
				}
				
				//need to look the found point up in vertex set V
				int found_index = index_in_vertices(V, found_point);

				
				
				
				
				
				
				
				
				if (solve_options.verbose_level()>=0) {
					vec_mp top_found, bottom_found;  init_vec_mp(top_found,md_config.num_top_vars());
					init_vec_mp(bottom_found,md_config.num_bottom_vars());

					int offset = md_config.num_mid_vars();
					// get only the bottom coordinates out of the returned point
					for (int tt = 0; tt<md_config.num_bottom_vars(); tt++) {
						set_mp(&bottom_found->coord[tt], & W_new.point(0)->coord[offset+tt]);
					}
					
					offset += md_config.num_bottom_vars();
					// get only the top coordinates out of the returned point
					for (int tt = 0; tt<md_config.num_top_vars(); tt++) {
						set_mp(&top_found->coord[tt], & W_new.point(0)->coord[offset+tt]);
					}
					print_point_to_screen_matlab(found_point,"found_point");
					std::cout << index_in_vertices(V, bottom_found) << " bottom_found" << std::endl;
					std::cout << index_in_vertices(V, top_found) << " top_found" << std::endl;
					std::cout << "found_index of point: " << found_index << std::endl;
					clear_vec_mp(bottom_found);
					clear_vec_mp(top_found);
				}
				
				if (solve_options.verbose_level()>=1) {
					
					projection_value_homogeneous_input(temp, found_point, pi(1));
					print_comp_matlab(temp, "found_point_proj_val");
					vec_mp tempvec; init_vec_mp(tempvec,0);
					dehomogenize(&tempvec, found_point);
					print_point_to_screen_matlab(tempvec, "found_point");
					clear_vec_mp(tempvec);
				}
				
				
				
				//search among the current edge possibilities for the one containing the found (mid) point
				int index_in_set = -1;
				for (std::set<int>::iterator possibility_iter=remaining_possible_edges.begin(); possibility_iter!=remaining_possible_edges.end(); possibility_iter++) {
					if (found_index == target_critslice.get_edge(*possibility_iter).midpt()) {
						index_in_set = *possibility_iter;
						break;
					}
				}
				
				
				// search edges for the found point as a removed point.
				if (index_in_set < 0) {
					index_in_set = target_critslice.edge_w_removed(found_index);
					
					if (index_in_set>=0 && solve_options.verbose_level()>=0) {
						std::cout << color::green() << "found point as removed point from edge " << index_in_set << color::console_default() << std::endl;
					}
				}
				

				//search among the overall set of possibilities.  this includes ones we have already tracked to.
				if (index_in_set < 0) {
					for (std::set<int>::iterator possibility_iter=possible_edges.begin(); possibility_iter!=possible_edges.end(); possibility_iter++) {
						if (found_index == target_critslice.get_edge(*possibility_iter).midpt()) {
							index_in_set = *possibility_iter;
							break;
						}
					}
				}
				
				if (index_in_set < 0 && solve_options.verbose_level()>=0) {
					std::cout << color::blue() << "did not find the indexed point as the midpoint of any current possibility." << color::console_default() << std::endl;
				}
				
				//perhaps check for the point as left or right point of an edge?
				
				if (index_in_set>=0 && found_edges.find(index_in_set)==found_edges.end())
				{
					int next_edge = index_in_set; // index the *edge*
					
					if (program_options.verbose_level()>=0) {
						std::cout << "added_edge " << next_edge << ", l m r: " << target_critslice.get_edge(next_edge).left() << " " << target_critslice.get_edge(next_edge).midpt() << " " << target_critslice.get_edge(next_edge).right() << "\n\n";
					}
					
					
					
					if ( (next_edge<0) || !(  (target_critslice.get_edge(next_edge).left()!=current_bottom_ind) ||  target_critslice.get_edge(next_edge).right()!=current_bottom_ind))  {
						continue;
					}
					
					if (target_critslice.get_edge(next_edge).left()==current_bottom_ind) {
						current_bottom_ind = current_top_ind = target_critslice.get_edge(next_edge).right(); // the upper value
					}
					else {
						current_bottom_ind = current_top_ind = target_critslice.get_edge(next_edge).left(); // the lower value
					}
					
					// keep track of those edges we found.
					found_edges.insert(next_edge);
					
					// erase both the currently testing edge, and the found one, from the list of possibilities.
					remaining_possible_edges.erase(current_edge);
					remaining_possible_edges.erase(next_edge);
					
					
					// add the next edge to the set we can connect together.
					if (zz==0) {
						F.add_left_edge(next_edge);
					}
					else
					{
						F.add_right_edge(next_edge);
					}
					
					break;
				}
				else
				{
					//didn't find, so simply remove from the list of possibilities.
					remaining_possible_edges.erase(current_edge);
				}
				
			}//re: for qq
			
		} // re: while...
		
		
		
		
	}
	
	clear_vec_mp(found_point);
	
	clear_mp(temp); clear_mp(temp2); clear_mp(temp3);
	clear_mp(numer); clear_mp(denom);
	
	clear_mp(proj_top); clear_mp(proj_bottom); clear_mp(proj_mid);
	
	return F;
}





void Surface::send(int target, ParallelismConfig & mpi_config) const
{
#ifdef functionentry_output
	std::cout << "surface::send" << std::endl;
#endif
	
	
	
	
	Decomposition::send(target, mpi_config);
	int * buffer = new int[4];
	
	buffer[0] = num_faces_;
	
	buffer[1] = mid_slices_.size();
	buffer[2] = crit_slices_.size();
	buffer[3] = num_singular_curves_;
	
	MPI_Send(buffer, 4, MPI_INT, target, SURFACE, mpi_config.comm());
	delete [] buffer;
	
	
	

	
	for (auto f=faces_.begin(); f!=faces_.end(); f++) {
		f->send(target, mpi_config);
	}
	
	
	for (auto s=mid_slices_.begin(); s!=mid_slices_.end(); s++) {
		s->send(target, mpi_config);
	}
	
	for (auto s=crit_slices_.begin(); s!=crit_slices_.end(); s++) {
		s->send(target, mpi_config);
	}
	
	crit_curve_.send(target, mpi_config);
	sphere_curve_.send(target, mpi_config);
	
	
	buffer = new int[2*num_singular_curves_];
	int counter = 0;
	for (auto iter = singular_curves_.begin(); iter!= singular_curves_.end(); ++iter) {
		buffer[counter] = iter->first.multiplicity();
		buffer[counter+1] = iter->first.index();
		counter+=2;
	}
	MPI_Send(buffer, 2*num_singular_curves_, MPI_INT, target, SURFACE, mpi_config.comm());
	delete [] buffer;
	
	for (auto iter = singular_curves_.begin(); iter!= singular_curves_.end(); ++iter) {
		iter->second.send(target, mpi_config);
	}
	
	
	return;
}


void Surface::receive(int source, ParallelismConfig & mpi_config)
{
	
#ifdef functionentry_output
	std::cout << "surface::receive" << std::endl;
#endif
	
	
	
	
	MPI_Status statty_mc_gatty;
	
	Decomposition::receive(source, mpi_config);
	
	
	int * buffer = new int[4];
	MPI_Recv(buffer, 4, MPI_INT, source, SURFACE, mpi_config.comm(), &statty_mc_gatty);
	int b, c, d;
	b = buffer[0];
	c = buffer[1];
	d = buffer[2];
	num_singular_curves_ = buffer[3];
	delete [] buffer;
	
	
	
	
	for (int ii=0; ii<b; ii++) {
		Face F;
		F.receive(source, mpi_config);
		add_face(F);
	}
	
	//TODO: rewrite this to prevent unnecessary copy operation.
	mid_slices_.resize(c);
	for (auto ii=mid_slices_.begin(); ii!=mid_slices_.end(); ii++) {
		ii->receive(source,mpi_config);
	}
	
	
	crit_slices_.resize(d);
	for (auto ii=crit_slices_.begin(); ii!=crit_slices_.end(); ii++) {
		ii->receive(source, mpi_config);
	}
	
	crit_curve_.receive(source, mpi_config);
	sphere_curve_.receive(source, mpi_config);
	
	
	buffer = new int[2*num_singular_curves_];//
	MPI_Recv(buffer, 2*num_singular_curves_, MPI_INT, source, SURFACE, mpi_config.comm(), &statty_mc_gatty);
	
	// receive and unpack the buffer at same time
	int counter = 0;
	for (unsigned int ii=0; ii<num_singular_curves_; ii++) {
		singular_curves_[SingularObjectMetadata(buffer[counter],buffer[counter+1])].receive(source, mpi_config);
		counter+=2;
	}
	
	delete [] buffer;
	return;
}




















void Surface::print(boost::filesystem::path base) const
{
	
#ifdef functionentry_output
	std::cout << "surface::print" << std::endl;
#endif
	
	
	
	
	Decomposition::print(base);
	
	
	
	boost::filesystem::path summaryname = base;
	summaryname /= "S.surf";
	FILE *OUT = safe_fopen_write(summaryname);
	fprintf(OUT,"%d 0 %ld %ld\n\n", num_faces(), mid_slices_.size(), crit_slices_.size());
	fprintf(OUT,"%ld\n",singular_curves_.size());
	for (auto iter = singular_curves_.begin(); iter!=singular_curves_.end(); ++iter) {
		fprintf(OUT,"%d %d ",iter->first.multiplicity(),iter->first.index());
	}
	fprintf(OUT,"\n");
	// what more to print here?
	fclose(OUT);
	
	summaryname = base;
	summaryname /= "F.faces";
	print_faces(summaryname);
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	for (unsigned int ii=0; ii!=mid_slices_.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices_[ii].output_main(specific_loc);
	}
	
	for (unsigned int ii=0; ii!=crit_slices_.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices_[ii].output_main(specific_loc);
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve_.output_main(specific_loc);
	
	specific_loc = curve_location;
	specific_loc += "_sphere";
	sphere_curve_.output_main(specific_loc);
	
	for (auto iter = singular_curves_.begin(); iter!=singular_curves_.end(); ++iter) {
		
		std::stringstream converter;
		converter << curve_location.string() << "_singular_mult_" << iter->first.multiplicity() << "_" << iter->first.index();
		
		specific_loc = converter.str();
		converter.clear(); converter.str("");
		
		iter->second.output_main(specific_loc);
	}
	
	
}







void Surface::print_faces(boost::filesystem::path outputfile) const
{
#ifdef functionentry_output
	std::cout << "surface::print_faces" << std::endl;
#endif
	
	
	
	
    //	std::cout << "printing faces to file " << outputfile << std::endl;
	
	
	std::ofstream fout(outputfile.c_str());
	fout << faces_.size() << std::endl << std::endl;
	for (auto iter = faces_.begin(); iter!=faces_.end(); ++iter) {
		fout << *iter << std::endl;
	}
	fout.close();
	
}






void Surface::read_faces(boost::filesystem::path load_from_me)
{
#ifdef functionentry_output
	std::cout << "surface::read_faces" << std::endl;
#endif
	
	
	
	std::ifstream fin(load_from_me.c_str());
	int temp_num_faces;
	fin >> temp_num_faces;
	for (int ii=0; ii<temp_num_faces; ii++) {
		Face F;
		
		fin >> F;
		add_face(F);
	}
	
	fin.close();

	
	
	return;
}






void Surface::setup(boost::filesystem::path base)
{
	Decomposition::setup(base / "decomp");
	
	std::vector<SingularObjectMetadata > singular_multiplicities;
	int temp_num_crit, temp_num_mid;
	
	read_summary(singular_multiplicities,temp_num_mid, temp_num_crit, base / "S.surf");
	
	read_faces(base / "F.faces");
	
	mid_slices_.resize(temp_num_mid);
	crit_slices_.resize(temp_num_crit);
	
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	std::stringstream converter;
	
	for (int ii=0; ii<temp_num_mid; ii++) {
		
		
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices_[ii].setup(specific_loc);
	}
	
	for (int ii=0; ii<temp_num_crit; ii++) {
		
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices_[ii].setup(specific_loc);
	}
	
	
	for (auto iter = singular_multiplicities.begin(); iter!=singular_multiplicities.end(); ++iter) {
		
		converter << iter->multiplicity() << "_" << iter->index();
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_singular_mult_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		singular_curves_[*iter].setup(specific_loc);
		num_singular_curves_++;
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve_.setup(specific_loc);
	
	specific_loc = curve_location;
	specific_loc += "_sphere";
	sphere_curve_.setup(specific_loc);
	
	
	
	return;
	
}





void read_summary(std::vector<SingularObjectMetadata > & singular_multiplicities, int & temp_num_mid, int & temp_num_crit, boost::filesystem::path INfile)
{
	FILE *IN = safe_fopen_read(INfile);
	int temp_num_faces, temp_num_edges;
	
	fscanf(IN,"%d %d %d %d", &temp_num_faces, &temp_num_edges, &temp_num_mid, &temp_num_crit);
	
	int temp_num_sing;
	fscanf(IN,"%d",&temp_num_sing);
	int temp_mult,temp_index;
	for (int ii=0; ii<temp_num_sing; ii++) {
		fscanf(IN,"%d %d",&temp_mult, &temp_index);
		singular_multiplicities.push_back(SingularObjectMetadata(temp_mult,temp_index));
	}
	fclose(IN);
	
	
	return;
}



















