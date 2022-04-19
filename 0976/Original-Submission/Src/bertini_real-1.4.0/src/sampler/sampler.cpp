#include "sampler.hpp"








void sampler_configuration::SetDefaults()
{
	no_duplicates = true;
	use_distance_condition = false;
	
	target_num_samples = 10;
			
	stifle_membership_screen = 1;
	stifle_text = " > /dev/null ";
	
	max_num_ribs = 20;
	min_num_ribs = 3;

	minimum_num_iterations = 2;
	maximum_num_iterations = 10;
	
	mpf_init(TOL);
	mpf_set_d(TOL, 1e-1); // this should be made adaptive to the span of the projection values or the endpoints
	
	use_gamma_trick = 0;

	mode = Mode::AdaptivePredMovement;
}



void sampler_configuration::splash_screen()
{
	printf("\n Sampler module for Bertini_real(TM) v%s\n\n", VERSION);
	printf(" D.A. Brake, \n with D.J. Bates, W. Hao, \n J.D. Hauenstein, A.J. Sommese, and C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	printf("See the website at www.bertinireal.com\n\n");
	printf("Send email to %s for assistance.\n\n",PACKAGE_BUGREPORT);
	
	
	
	bertini_splash_screen();
	
	
}


void sampler_configuration::print_usage()
{
	std::cout << "Bertini_real has the following options:\n";
	std::cout << "option name(s)\t\t\targument\n\n";
	std::cout << "-ns -nostifle\t\t\t   --\n";
	std::cout << "-v -version\t\t\t   -- \n";
	std::cout << "-h -help\t\t\t   --\n";
	std::cout << "-t -tol -tolerance \t\tdouble > 0\n";
	std::cout << "-verb\t\t\t\tint\n";
	std::cout << "-minits \t\t\tint minimum number of passes for adaptive curve or surface refining\n";
	std::cout << "-maxits \t\t\tint maximum number of passes for adaptive curve or surface refining\n";
	std::cout << "-maxribs \t\t\tint maximum number of ribs for adaptive surface refining\n";
	std::cout << "-minribs \t\t\tint minimum number of ribs for adaptive surface refining\n";
	std::cout << "-gammatrick -g \t\t\tbool\n";
	std::cout << "-numsamples \t\t\tint number samples per edge\n";
	std::cout << "-mode -m \t\t\tchar sampling mode.  ['a'] adaptive by movement, 'd' adaptive by distance, 'f' fixed, \n";
	std::cout << "\n\n\n";
	std::cout.flush();
	return;
}

int  sampler_configuration::parse_commandline(int argc, char **argv)
{
	// this code created based on gnu.org's description of getopt_long
	int choice;
	
	while (1)
	{
		static struct option long_options[] =
		{
			/* These options set a flag. */
			{"nostifle", no_argument,       0, 's'},
			{"ns", no_argument,       0, 's'},
			{"help",		no_argument,			 0, 'h'},
			{"h",		no_argument,			 0, 'h'},
			{"version",		no_argument,			 0, 'v'},
			{"v",		no_argument,			 0, 'v'},
			{"verb",		required_argument,			 0, 'V'},
			{"tolerance",		required_argument,			 0, 't'},
			{"tol",		required_argument,			 0, 't'},
			{"t",		required_argument,			 0, 't'},
			{"minits",		required_argument,			 0, 'l'},
			{"maxits",		required_argument,			 0, 'm'},
			{"maxribs",		required_argument,			 0, 'R'},
			{"minribs",		required_argument,			 0, 'r'},
			{"gammatrick",		required_argument,			 0, 'g'},
			{"g",		required_argument,			 0, 'g'},
			{"numsamples",		required_argument,			 0, 'n'},
			{"nd", no_argument,0,'d'},
			{"m",		required_argument,			 0, 'M'},
			{"mode",		required_argument,			 0, 'M'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		choice = getopt_long_only (argc, argv, "bdf:svt:g:V:l:m:R:r:hM:", // colon requires option, two is optional
															 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (choice == -1)
			break;
		
		switch (choice)
		{				
			case 'd':
				no_duplicates = false;
				break;
				
				
			case 'n':
				target_num_samples = atoi(optarg);
				
				if (target_num_samples <= 3) {
					std::cout << "The number of desired samples must be larger than 3, but you provided " << target_num_samples << std::endl;
					exit(0);
				}
				break;
				
			case 's':
				this->stifle_text = "\0";
				break;
				
			case 'v':
				printf("\n Sampler module for Bertini_real(TM) version %s\n\n", VERSION);
				std::cout << "for help, use option '-h'\n\n";
				exit(0);
				break;
				
			case 't':
				
				mpf_set_str(this->TOL,optarg,10);
				break;
				
			case 'g':
				this->use_gamma_trick = atoi(optarg);
				if (! (this->use_gamma_trick==0 || this->use_gamma_trick==1) ) {
					printf("value for 'gammatrick' or 'g' must be 1 or 0\n");
					exit(0);
				}
				break;
				
			case 'V':
				this->verbose_level(atoi(optarg));
				break;

			case 'l':
				this->minimum_num_iterations = atoi(optarg);
				break;

			case 'm':
				this->maximum_num_iterations = atoi(optarg);
				break;
			
			case 'M':
			{
				std::string curr_opt{optarg};
				if (curr_opt.size()!=1){
					std::cout << "option to 'mode' must be a single character.  valid modes are 'a' and 'f'\n";
					exit(0);
				}
				
				switch (curr_opt[0])
				{
					case 'd':
						mode = Mode::AdaptiveConsecDistance;
						break;

					case 'a':
						mode = Mode::AdaptivePredMovement;
						break;

					case 'f':
						mode = Mode::Fixed;
						break;
				}
				break;
			}

			case 'R':
				this->max_num_ribs = atoi(optarg);
				break;

			case 'r':
				this->min_num_ribs = atoi(optarg);
				break;

			case 'h':
				
				sampler_configuration::print_usage();
				exit(0);
				break;
				
			case '?':
				/* getopt_long already printed an error message. */
				break;
				
			default:
				sampler_configuration::print_usage();
				exit(0);
		}
	}
	
	/* Instead of reporting ‘--verbose’
	 and ‘--brief’ as they are encountered,
	 we report the final status resulting from them. */
	
	
	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
			printf ("%s ", argv[optind++]);
		putchar ('\n');
	}
	
	return 0;
}







int main(int argC, char *args[])
{
	
	boost::timer::auto_cpu_timer t;
	
	
	MPI_Init(&argC,&args);
	
	
	boost::filesystem::path inputName, witnessSetName, samplingNamenew;
	
	
	
	
	
	
	
	sampler_configuration sampler_options;
	
	sampler_options.splash_screen();
	sampler_options.parse_commandline(argC, args);
    
	
	
	int MPType, dimension;
	
	boost::filesystem::path directoryName;
	
	get_dir_mptype_dimen( directoryName, MPType, dimension);
	
	
	
	
	
    
	
	witnessSetName = directoryName / "WitnessSet";
	samplingNamenew = directoryName;
	
	
	SolverConfiguration solve_options;
	
	
	//only one will be used.  i don't know how to avoid this duplicity.
	Curve C;
	Surface surf_input;
	
	Decomposition * decom_pointy; // this feels unnecessary
	switch (dimension) {
		case 1:
		{
			C.setup(directoryName);
			decom_pointy = &C;
			
		}
			break;
			
		case 2:
		{
			surf_input.setup(directoryName);
			decom_pointy = &surf_input;
			
		}
			break;
		case -1:
		{
			std::cout << "sampler unable to proceed\n";
			return 12461;
		}
		default:
		{
			std::cout << "sampler not capable of sampling anything but dimension 1 or 2.  this is of dim " << dimension << std::endl;
			return 12462;
		}
			break;
	}
	
	
	
	
	common_sampler_startup(*decom_pointy,
						   sampler_options,
						   solve_options);
	
	
	
	VertexSet V(decom_pointy->num_variables());
	
	V.set_tracker_config(&solve_options.T);
	V.setup_vertices(directoryName / "V.vertex"); //setup V structure from V.vertex
	V.set_same_point_tolerance(1e1*solve_options.T.real_threshold);
	
	
	/////////
	////////
	//////
	////
	//
	//  Generate new sampling data
	//
    

	switch (dimension) {
		case 1:
		{
			switch (sampler_options.mode){
				case sampler_configuration::Mode::Fixed:
					C.fixed_sampler(V,
									sampler_options,
									solve_options,
									sampler_options.target_num_samples);
					break;
				case sampler_configuration::Mode::AdaptiveConsecDistance:

					C.adaptive_sampler_distance(V,
												sampler_options,
												solve_options);
					break;

				case sampler_configuration::Mode::AdaptivePredMovement:
					C.adaptive_sampler_movement(V,
												sampler_options,
												solve_options);
					break;
			} // switch
			C.output_sampling_data(directoryName);
			V.print(directoryName / "V_samp.vertex");
			
			break;
		}
			
			
		case 2:
		{
			switch (sampler_options.mode){
				case sampler_configuration::Mode::Fixed:
				{
					surf_input.fixed_sampler(V,
											 sampler_options,
											 solve_options);
					
					break;
				}
				case sampler_configuration::Mode::AdaptivePredMovement:
					std::cout << "adaptive by movement not implemented for surfaces, using adaptive by distance\n\n";
				case sampler_configuration::Mode::AdaptiveConsecDistance:
				{
					surf_input.AdaptiveSampler(V,
											 sampler_options,
											 solve_options);
					break;
				}

			} // switch

			surf_input.output_sampling_data(directoryName);
			V.print(directoryName / "V_samp.vertex");

			break;
		}	
		default:
			break;
	}
    
	
	
	//
	//   Done with the main call
	////
	/////
	///////
	////////

	
	clearMP();
	MPI_Finalize();
	
	return 0;
}



void common_sampler_startup(const Decomposition & D,
							sampler_configuration & sampler_options,
							SolverConfiguration & solve_options)
{
	
	parse_input_file(D.input_filename()); // restores all the temp files generated by the parser, to this folder.  i think this can be removed?  but the PPD and tracker config both depend on it...  so maybe not.
	
    
	
	// set up the solver configuration
	get_tracker_config(solve_options, solve_options.T.MPType);
	
	initMP(solve_options.T.Precision);
	
	
	parse_preproc_data("preproc_data", &solve_options.PPD);
    
	
	
	
	
	solve_options.verbose_level(sampler_options.verbose_level());

	
	
	solve_options.T.ratioTol = 1; // manually assert to be more permissive.  i don't really like this.
	
	
	
	solve_options.use_midpoint_checker = 0;
	solve_options.use_gamma_trick = 0;
    solve_options.robust = true;
	
}








//dehomogenizes, takes the average, computes the projection.
//takes in the full projection \pi, including the homogenizing coordinate.
void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi){
	int ii;
	
	if (left->size < pi->size) {
		printf("left point too short in estimate new projection value\n");
		deliberate_segfault();
	}
	
	if (right->size < pi->size) {
		printf("left point too short in estimate new projection value\n");
		deliberate_segfault();
	}
    //	print_point_to_screen_matlab(left,"left");
    //	print_point_to_screen_matlab(right,"right");
    //	print_point_to_screen_matlab(pi,"pi");
	
	vec_mp dehom_left, dehom_right;
	init_vec_mp(dehom_left,pi->size-1);   dehom_left->size  = pi->size-1;
	init_vec_mp(dehom_right,pi->size-1);  dehom_right->size = pi->size-1;
	
	dehomogenize(&dehom_left,left,pi->size);
	dehomogenize(&dehom_right,right,pi->size);
	
	
	comp_mp temp, temp2, half; init_mp(temp); init_mp(temp2);  init_mp(half);
	
	mpf_set_d(half->r, 0.5); mpf_set_d(half->i, 0.0);
	
    
	set_zero_mp(result);                                           // result = 0;  initialize
	
	for (ii = 0; ii<pi->size-1; ii++) {
		add_mp(temp,&dehom_left->coord[ii],&dehom_right->coord[ii]); //  a = (x+y)
		mul_mp(temp2, temp, half);                                   //  b = a/2
		mul_mp(temp,&pi->coord[ii+1],temp2);                          //  a = b.pi
		set_mp(temp2,result);                                        //  b = result
		add_mp(result, temp, temp2);                                  //  result = a+b
	}
	// in other words, result += (x+y)/2 \cdot pi


    real_threshold(result,1e-9);

	
	clear_mp(temp); clear_mp(temp2); clear_mp(half);
	clear_vec_mp(dehom_right);clear_vec_mp(dehom_left);
	
	return;
}



//dehomogenizes, takes the average, computes the projection.
//takes in the full projection \pi, including the homogenizing coordinate.
// this version returns the estimated point as well.
void estimate_new_projection_value(comp_mp result, vec_mp estimated_point, vec_mp left, vec_mp right, vec_mp pi){

	
	if (left->size != right->size) {
		printf("attempting to estimate_new_projection_value on vectors of different size\n");
		br_exit(98128);
	}
	
	comp_mp temp;  init_mp(temp);
	
	comp_mp half; init_mp(half);
	mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);
	
	
	
	vec_mp dehom_left, dehom_right;
	init_vec_mp(dehom_left,left->size-1);   dehom_left->size  = left->size-1;
	init_vec_mp(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left,pi->size);
	dehomogenize(&dehom_right,right,pi->size);
	
	
	
	vec_add_mp(estimated_point, dehom_left, dehom_right);
	vec_mulcomp_mp(estimated_point, estimated_point, half);
	
	set_zero_mp(result);                                           // result = 0;  initialize
	
	for (int ii = 0; ii<pi->size-1; ii++) {                                 //  b = a/2
		mul_mp(temp,&pi->coord[ii+1],&estimated_point->coord[ii]);                          //  a = b.pi
		add_mp(result, result, temp);                                  //  result = a+b
	}
    
	// in other words, result += (x+y)/2 \cdot pi
    
	
	//  i think this thresholding should be moved to outside this call.
	mpf_t zerothresh; mpf_init(zerothresh);
	mpf_set_d(zerothresh, 1e-9);
	if (mpf_cmp(result->i, zerothresh) < 0){
		mpf_set_str(result->i, "0.0", 10);
	}
	
	mpf_clear(zerothresh);
	
	clear_mp(temp); clear_mp(half);
	
	clear_vec_mp(dehom_right);clear_vec_mp(dehom_left);
	
	return;
}















void triangulate_two_ribs_by_angle_optimization(const std::vector< int > & rib1, const std::vector< int > & rib2,
											  VertexSet & V, double real_thresh,
											  std::vector< Triangle> & current_samples)
{
#ifdef functionentry_output
	std::cout << "triangulate_by_distance_binning" << std::endl;
#endif
	
	
	bool bail_out = false;
	
	if (rib1.size()==0) {
		std::cout << "rib1 had 0 size!" << std::endl;
		bail_out = true;
	}
	if (rib2.size()==0) {
		std::cout << "rib2 had 0 size!" << std::endl;
		bail_out = true;
	}
	
	if (rib1.size()==1 && rib2.size()==1) {
		std::cout << "both ribs have size 1!" << std::endl;
		bail_out = true;
	}
	
	if (bail_out) {
		return;
	}
	
	
	int num_vars = V.num_natural_variables();
	
	unsigned int num_vectors_needed = 9;
	vec_mp * bulk_vectors = (vec_mp *) br_malloc(num_vectors_needed* sizeof(vec_mp));
	
	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		init_vec_mp(bulk_vectors[ii],num_vars);  bulk_vectors[ii]->size = num_vars;
		
	}
	
	vec_mp *A = &bulk_vectors[0], *B = &bulk_vectors[1], *C = &bulk_vectors[2], *D = &bulk_vectors[3];


	vec_mp *AB = &bulk_vectors[4];
	vec_mp *AC = &bulk_vectors[5];
	vec_mp *BC = &bulk_vectors[6];
	vec_mp *DA = &bulk_vectors[7];
	vec_mp *DC = &bulk_vectors[8];

	
	
	// seed the advancing loop.
	dehomogenize(B, V[rib1[0]].point(), num_vars);
	(*B)->size = num_vars-1;
	real_threshold(*B,real_thresh);
	dehomogenize(D, V[rib2[0]].point(), num_vars);
	(*D)->size = num_vars-1;
	real_threshold(*D,real_thresh);
	
	
	comp_mp cos_angle_CAB, cos_angle_BCA, cos_angle_ABC;
	init_mp(cos_angle_CAB); init_mp(cos_angle_BCA); init_mp(cos_angle_ABC);
	
	comp_mp cos_angle_DCA, cos_angle_ADC, cos_angle_CAD;
	init_mp(cos_angle_DCA); init_mp(cos_angle_ADC); init_mp(cos_angle_CAD);
	
	
	comp_mp length_AB, length_AC, length_BC, length_DA, length_DC;
	init_mp(length_AB); init_mp(length_AC); init_mp(length_BC); init_mp(length_DA); init_mp(length_DC);
	
	comp_mp dot_CAB, dot_BCA, dot_ABC;
	init_mp(dot_CAB); init_mp(dot_BCA); init_mp(dot_ABC);
	
	comp_mp dot_DCA, dot_ADC, dot_CAD;
	init_mp(dot_DCA); init_mp(dot_ADC); init_mp(dot_CAD);
	
	comp_mp temp;  init_mp(temp);
	comp_d temp_d;
	
	unsigned int curr_index_rib1 = 0, curr_index_rib2 = 0;
	bool moved_1 = true, moved_2 = true;  //this is an intial condition to get them both set properly.  all subsequent iterations have only one as moved==true, and the other is always false.
	while (curr_index_rib1 < rib1.size()-1 // neither rib size is 0, so this -1 is ok, won't underflow
		   &&
		   curr_index_rib2 < rib2.size()-1)
	{

#ifdef debug_compile
		int a =rib1[curr_index_rib1];
		int b =rib1[curr_index_rib1+1];
		int c =rib2[curr_index_rib2];
		int d =rib2[curr_index_rib2+1];
		
		std::cout << rib1[curr_index_rib1] << " " << rib1[curr_index_rib1+1] << std::endl;
		std::cout << rib2[curr_index_rib2] << " " << rib2[curr_index_rib2+1] << std::endl;
#endif
		
		if (moved_1) {
			vec_mp * temp_vec = A; // swap
			A = B;
			B = temp_vec;
			dehomogenize(B, V[rib1[curr_index_rib1+1]].point(), num_vars);
			(*B)->size = num_vars-1;
			real_threshold(*B,real_thresh);

			
			vec_sub_mp(*AB, *B,*A);
			
			
			
			twoNormVec_mp(*AB, length_AB);
		}
							   
		if (moved_2){ // moved_smaller
			vec_mp * temp_vec = C; // swap
			C = D;
			D = temp_vec;
			dehomogenize(D, V[rib2[curr_index_rib2+1]].point(), num_vars);
			(*D)->size = num_vars-1;
			
			real_threshold(*D,real_thresh);
			
			vec_sub_mp(*DC, *C,*D);
			twoNormVec_mp(*DC, length_DC);
		}
		
		
		
//		A           --->           B
//		  ***********************
//		  *- <--              *
//		  * -  --           *
//		  *  -   --       *    --
//		| *   -         *   ---
//		| *    -      *  <--
//		| *     -   *
//		\/*      -*
//		  *     * -
//		  *   *    -
//		  * *       -
//		  *-----------
//		C    <--       D
		
		
		
//        AC, BC, DA,
		
		vec_sub_mp(*AC, *C,*A);
		vec_sub_mp(*BC, *C,*B);
		vec_sub_mp(*DA, *A,*D);
		
		
#ifdef debug_compile
		print_point_to_screen_matlab(*A,"A");
		print_point_to_screen_matlab(*B,"B");
		print_point_to_screen_matlab(*C,"C");
		print_point_to_screen_matlab(*D,"D");
		
		print_point_to_screen_matlab(*AB,"AB");
		print_point_to_screen_matlab(*AC,"AC");
		print_point_to_screen_matlab(*BC,"BC");
		print_point_to_screen_matlab(*DA,"DA");
		print_point_to_screen_matlab(*DC,"DC");
#endif
		
		
		// now have the 5 vectors for the test computed.  (5 because the two triangles share a common leg)
		
		dot_product_mp(dot_CAB, *AC, *AB);
		
		dot_product_mp(dot_BCA, *AC, *BC);
		
		dot_product_mp(dot_ABC, *AB, *BC);
		neg_mp(dot_ABC,dot_ABC);
		
		dot_product_mp(dot_DCA, *DC, *AC);
		
		dot_product_mp(dot_ADC, *DA, *DC);
		
		
		dot_product_mp(dot_CAD, *DA, *AC);
		neg_mp(dot_CAD,dot_CAD);
		
		twoNormVec_mp(*AC, length_AC);
		twoNormVec_mp(*BC, length_BC);
		twoNormVec_mp(*DA, length_DA);
		
		
		double thresh = 3;
		int advance = 0;
		bool aspect_ok_rib1 = true, aspect_ok_rib2 = true;
		
		div_mp(temp,length_BC,length_AB);
		mp_to_d(temp_d, temp);
//		double aspect11 = temp_d->r;
		if ( temp_d->r  > thresh) {
			aspect_ok_rib1 = false;
		}
		div_mp(temp,length_BC,length_AC);
		mp_to_d(temp_d, temp);
//		double aspect12 = temp_d->r;
		if ( temp_d->r  > thresh) {
			aspect_ok_rib1 = false;
		}
//
		
		
		div_mp(temp,length_DA,length_DC);
		mp_to_d(temp_d, temp);
//		double aspect21 = temp_d->r;
		if (temp_d->r > thresh) {
			aspect_ok_rib2 = false;
		}
		div_mp(temp,length_DA,length_AC);
		mp_to_d(temp_d, temp);
//		double aspect22 = temp_d->r;
		if (temp_d->r > thresh) {
			aspect_ok_rib2 = false;
		}

		
		
		
		
		
		// the below lines computes the sum of the squares of the differences of the absolute values of the cosines of two of the angles in one of the triangles.
		
		//CAB
		double total_error_rib1 = compute_square_of_difference_from_sixtydegrees(temp, length_AC, length_AB, dot_CAB);
#ifdef debug_compile
		print_comp_matlab(dot_CAB,"CAB");
		std::cout << c << " " << a << " " << b << " " << total_error_rib1 << std::endl;
#endif
		
		
		//BCA
		double angle_BCA = compute_square_of_difference_from_sixtydegrees(temp, length_BC, length_AC, dot_BCA);
#ifdef debug_compile
		print_comp_matlab(dot_BCA,"BCA");
		std::cout << b << " " << c << " " << a << " " << angle_BCA << std::endl;
#endif
		total_error_rib1 += angle_BCA;
		
		//ABC
		double angle_ABC = compute_square_of_difference_from_sixtydegrees(temp, length_BC, length_AB, dot_ABC);
		total_error_rib1 += angle_ABC;
#ifdef debug_compile
		print_comp_matlab(dot_ABC,"ABC");
		std::cout << a << " " << b << " " << c << " " << angle_ABC << std::endl;
#endif
		
		
		
		
		
		
		//DCA
		double total_error_rib2 = compute_square_of_difference_from_sixtydegrees(temp, length_DC, length_AC, dot_DCA);
#ifdef debug_compile
		print_comp_matlab(dot_DCA,"DCA");
		std::cout << d << " " << c << " " << a << " " << total_error_rib2 << std::endl;
#endif
		
		//ADC
		
		double angle_ADC = compute_square_of_difference_from_sixtydegrees(temp, length_DA, length_DC, dot_ADC);
		
		total_error_rib2+= angle_ADC;
#ifdef debug_compile
		print_comp_matlab(dot_ADC,"ADC");
		std::cout << a << " " << d << " " << c << " " << angle_ADC << std::endl;
#endif
		//CAD
		double angle_CAD = compute_square_of_difference_from_sixtydegrees(temp, length_AC, length_DA, dot_CAD);
		total_error_rib2 += angle_CAD;
		
#ifdef debug_compile
		print_comp_matlab(dot_CAD,"CAD");
		std::cout << c << " " << a << " " << d << " " << angle_CAD << std::endl;
#endif
		
		
		
		
		
		if ((total_error_rib1>0.25) && (total_error_rib2>0.25)){
			// both triangles more or less equilateral, or both are bad
			if ( mpf_cmp(length_DA->r,length_BC->r)<0)
				advance = 2;
			else
				advance = 1;
		}
		else{
			if (total_error_rib1 < total_error_rib2) {
				advance = 1;
			}
			else
				advance = 2;
		}
		
		
//		if(!aspect_ok_rib1 && aspect_ok_rib2) {
//			advance = 2;
//		}
//		else if(!aspect_ok_rib2 && aspect_ok_rib1) {
//			advance = 1;
//		}
//		else if ( (!aspect_ok_rib2 && !aspect_ok_rib1) )
//		{
//			
//			
//			
//				
//			
//		}
//		else
//		{
//			
//			
//			if (total_error_rib1 < total_error_rib2) {
//				advance = 1;
//			}
//			else{
//				advance = 2;
//			}
//		}
		
		
		
		
		if (advance==1) { // if the 1 Triangle is more equilateral than the 2 Triangle.
			current_samples.push_back(
									  Triangle(// Triangle A B C
											   rib1[curr_index_rib1], //A
											   rib1[curr_index_rib1+1], //B
											   rib2[curr_index_rib2]) //C
									  );
			moved_1 = true;  curr_index_rib1++;
			moved_2 = false;
			
		}
		else
		{
			current_samples.push_back(
									  Triangle(// Triangle C A D
											   rib2[curr_index_rib2], //C
											   rib1[curr_index_rib1], //A
											   rib2[curr_index_rib2+1]) //D
									  );
			moved_1 = false;
			moved_2 = true; curr_index_rib2++;
		}
		
	} // re: while loop.
	
	
	
	// now down here, we have triangulated until one of the ribs is on its last point, so there is no more testing that can be done.  you simply have to connect the rest into triangles.
	
	const std::vector< int > *exhausted_rib, *rib_still_going;
	unsigned int index_still_going, terminal_index;
	bool flip;
	
	if (curr_index_rib1==rib1.size()-1) {
		exhausted_rib = &rib1;
		rib_still_going = &rib2;
		index_still_going = curr_index_rib2;
		terminal_index = curr_index_rib1;
		flip = false;
	}
	else
	{
		exhausted_rib = &rib2;
		rib_still_going = &rib1;
		index_still_going = curr_index_rib1;
		terminal_index = curr_index_rib2;
		flip = true;
	}
	
	
	for (; index_still_going < rib_still_going->size()-1; index_still_going++) { // initializer deliberately empty
		
		long long v1, v2, v3;
		if (flip) {
			v1 = rib_still_going->at(index_still_going);
			v2 = rib_still_going->at(index_still_going+1);
		}
		else
		{
			v1 = rib_still_going->at(index_still_going+1);
			v2 = rib_still_going->at(index_still_going);
		}
		v3 = exhausted_rib->at(terminal_index);
		
		current_samples.push_back( Triangle(v1,v2,v3) );
	}
	
	
	
	
	clear_mp(cos_angle_CAB); clear_mp(cos_angle_BCA); clear_mp(cos_angle_ABC);
	clear_mp(cos_angle_DCA); clear_mp(cos_angle_ADC); clear_mp(cos_angle_CAD);
	
	
	clear_mp(length_AB); clear_mp(length_AC); clear_mp(length_BC); clear_mp(length_DA); clear_mp(length_DC);
	clear_mp(dot_CAB); clear_mp(dot_BCA); clear_mp(dot_ABC);
	clear_mp(dot_DCA); clear_mp(dot_ADC); clear_mp(dot_CAD);
	
	clear_mp(temp);
	// clean up at the end.  i wish scope was deletion!
	
	for (unsigned int ii=0; ii<num_vectors_needed; ii++) {
		clear_vec_mp(bulk_vectors[ii]);
	}
	free(bulk_vectors);
	
	
	
	return;
}





double compute_square_of_difference_from_sixtydegrees(comp_mp temp, comp_mp length1, comp_mp length2, comp_mp dot_prod)
{
	comp_d cos;
	
	mul_mp(temp, length1, length2);
	div_mp(temp, dot_prod, temp);
	mp_to_d(cos, temp);
	double error = fabs(cos->r - 0.5);
	
	
//	double angle = acos(cos->r);
//	double pi_over_three = acos(-1)/3;
//	double error = fabs(angle-pi_over_three);

	
//	div_mp(temp, length1,length2)
//	mp_to_d(cos,temp);
//	double error2 = fabs(cos->r - 1);
//	
//	div_mp(temp, length2,length1)
//	mp_to_d(cos,temp);
//	double error3 = fabs(cos->r - 1);
//	
//	double error = error1 + error2 + error3;
	
	return error*error;//*error*error
}





// if x> 0.5
// 	cycle_num = c2;
// 	pi_out = 1;
// 	pi_mid = 0.5;
// 	p = (x-0.5)*2;
// 	r = pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
// else
// 	cycle_num = c1;
// 	pi_out = 0;
// 	pi_mid = 0.5;
// 	p = x*2;
// 	r = pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
// end
void ScaleByCycleNum(comp_mp result, comp_mp input, int cycle_num_l, int cycle_num_r)
{
	comp_mp one;  init_mp(one); set_one_mp(one);
	comp_mp two;  init_mp(two); mpf_set_str(two->r, "2.0", 10); mpf_set_str(two->i, "0.0", 10);
	comp_mp half;  init_mp(half); mpf_set_str(half->r, "0.5", 10); mpf_set_str(half->i, "0.0", 10);

	comp_mp temp1, temp2; init_mp(temp1); init_mp(temp2);
	comp_mp p; init_mp(p);

	// pi_out + (pi_mid - pi_out) * (1-p)^cycle_num;
	if ( mpf_get_d(input->r)>=0.5 )
	{
		// p = (x-0.5)*2;
		sub_mp(temp1, input, half);
		mul_mp(p, temp1 ,two);

		sub_mp(temp1, one, p);
		exp_mp_int(temp2, temp1, cycle_num_r);
		mul_mp(temp1, half, temp2);
		// 1 - (0.5) * (1-p)^cycle_num_r;


		sub_mp(result, one, temp1);
	}
	else
	{
		// p = x*2
		sub_mp(temp1, half, input);
		mul_mp(p, temp1, two);
		sub_mp(temp1, one, p);
		exp_mp_int(temp2, temp1, cycle_num_l);
		mul_mp(result, half, temp2);
		// r = 0.5 * (1-p)^cycle_num
	}

	clear_mp(one); clear_mp(half); clear_mp(temp1); clear_mp(temp2); clear_mp(p); clear_mp(two);
}







void set_witness_set_mp(WitnessSet & W, vec_mp new_linear, vec_mp new_point)
{
	W.reset_points();
	W.add_point(new_point);
	
	W.reset_linears();
	W.add_linear(new_linear);	
}








int get_dir_mptype_dimen(boost::filesystem::path & Dir_Name, int & MPType, int & dimension){
    
	std::string tempstr;
	std::ifstream fin("Dir_Name");

	if (!fin.is_open())
	{
		std::cout << color::red() << "did not find a decomposition in this directory.  currently uses file `Dir_Name` to record the name of the directory of the decomposition.  please ensure you have a completed decomposition in this directory, and `Dir_Name` is intact, specifying the correct location.\n" << color::console_default();
		dimension = -1;
		MPType = -1;
		Dir_Name = "missing";
		return MPType;
	}

	fin >> tempstr;
	fin >> MPType;
	fin >> dimension;
	
	Dir_Name = tempstr;
	Dir_Name = Dir_Name.filename();
	return MPType;
}




















