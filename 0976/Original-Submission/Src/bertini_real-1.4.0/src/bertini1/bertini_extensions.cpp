
#include "bertini1/bertini_extensions.hpp"








std::string enum_lookup(int flag, int hint)
{
	if (hint)
	{
		std::stringstream s;
	for (const auto& t : VertexTypes)
		if (flag & t)
			s << t << " ";
	return s.str();
	}
	else
	{
	switch (flag) {
		case SUCCESSFUL:
			return "SUCCESSFUL";
			break;
			
		case CRITICAL_FAILURE:
			return "CRITICAL_FAILURE";
			break;
			
		case TOLERABLE_FAILURE:
			return "TOLERABLE_FAILURE";
			break;
			
		case NULLSPACE:
			return "NULLSPACE";
			break;
			
		case LINPRODTODETJAC:
			return "LINPRODTODETJAC";
			break;
			
		case DETJACTODETJAC:
			return "DETJACTODETJAC";
			break;
			
		case LINTOLIN:
			return "LINTOLIN";
			break;
			
		case MULTILIN:
			return "MULTILIN";
			break;
			
		case MIDPOINT_SOLVER:
			return "MIDPOINT_SOLVER";
			break;
			
		case SPHERE_SOLVER:
			return "SPHERE_SOLVER";
			break;
			
		case BERTINI_MAIN:
			return "BERTINI_MAIN";
			break;
		case TERMINATE:
			return "TERMINATE";
			break;
			
		case INITIAL_STATE:
			return "INITIAL_STATE";
			break;
			
			
		case PARSING:
			return "PARSING";
			break;
			
		case TYPE_CONFIRMATION:
			return "TYPE_CONFIRMATION";
			break;
			
		case DATA_TRANSMISSION:
			return "DATA_TRANSMISSION";
			break;
			
		case NUMPACKETS:
			return "NUMPACKETS";
			break;
			
		case INACTIVE:
			return "INACTIVE";
			break;
			
		case VEC_MP:
			return "VEC_MP";
			break;
			
		case VEC_D:
			return "VEC_D";
			break;
			
		case MAT_MP:
			return "MAT_MP";
			break;
			
		case MAT_D:
			return "MAT_D";
			break;
			
		case COMP_MP:
			return "COMP_MP";
			break;
			
		case COMP_D:
			return "COMP_D";
			break;
			
		case INDICES:
			return "INDICES";
			break;
			
			
		default:
			break;
	}
	}
	return "unknown...  check out data_type.cpp";
}



int bertini_main_wrapper(const std::vector<std::string> & options, int num_processes, int my_id, int headnode)
{
	if (my_id == headnode)
	{ // headnode controls the overall execution
		
		
		std::cout << "calling bertini's main with options:" << std::endl;
		// convert the vector to char array
		std::vector< char const * > args;  args.push_back(NULL);
		for (auto iter=options.begin(); iter!=options.end(); iter++) {
			std::cout << *iter << std::endl;
			
			args.push_back(iter->c_str());
		}
		
		int argC = options.size()+1;
		
		bclock_t time1, time2;
		int trackType, genType = 0, MPType, userHom, sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
		double parse_time;
		unsigned int currentSeed;
		char *inputName = NULL, *startName = NULL;
		
		
		
		// standard execution
		bclock(&time1);
		
		// setup inputName
		if (argC >= 2 && args[1] != NULL)
		{ // inputName is args[1]
			inputName = (char *)bmalloc((strlen(args[1]) + 1) * sizeof(char));
			strcpy(inputName, args[1]);
			
			// setup startName
			if (argC >= 3 && args[2] != NULL)
			{ // startName is args[2]
				startName = (char *)bmalloc((strlen(args[2]) + 1) * sizeof(char));
				strcpy(startName, args[2]);
			}
			else
			{ // default to 'start'
				startName = (char *)bmalloc(6 * sizeof(char));
				strcpy(startName, "start");
			}
		}
		else
		{ // default to 'input' & 'start'
			inputName = (char *)bmalloc(6 * sizeof(char));
			strcpy(inputName, "input");
			startName = (char *)bmalloc(6 * sizeof(char));
			strcpy(startName, "start");
		}
		
		// parse the input file and seed the random number generator
		parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
		
		// remove the output files from possibly previous runs - delete only files not needed
		remove_output_files(trackType, sharpenOnly, 1);
		
		bclock(&time2);
		totalTime(&parse_time, time1, time2);
		
		// call the main functions
		if (sharpenOnly)
		{ // sharpen the endpoints from a previous run
			sharpen_process_main(MPType, trackType, currentSeed, my_id, num_processes, headnode);
		}
		else
		{ // do either function evaluation, zero dimensional or positive dimensional tracking
			if (trackType == 0)
			{ // zero dimensional tracking
				zero_dim_main(MPType, parse_time, currentSeed, startName, my_id, num_processes, headnode);
			}
			else if (1 <= trackType && trackType <= 7)
			{ // positive dimensional tracking
				pos_dim_main(trackType, genType, MPType, currentSeed, startName, my_id, num_processes, headnode);
			}
			else if (trackType == -4 || trackType == -3)
			{ // function evaluation
				function_eval_main(trackType == -3, MPType, currentSeed, startName, my_id, num_processes, headnode);
			}
			else if (trackType == -2 || trackType == -1)
			{ // newton evaluation
				newton_eval_main(trackType == -1, MPType, currentSeed, startName, my_id, num_processes, headnode);
			}
		}
		
		free(inputName);
		free(startName);
		
		if (remove_temp)  // remove temporary files
			remove_temp_files();
		
	}
	else
	{ // worker process
		parallel_diff_worker(my_id, num_processes, headnode);
		
		worker_process_main(my_id, num_processes, headnode);
	}
	
	initMP(mpf_get_default_prec());
	
	return 0;
	
	
}







void bertini_splash_screen()
{
	printf("\n  Library-linked Bertini(TM) v%s", BERTINI_VERSION_STRING);
	printf("\n   (%s)\n\n", BERTINI_DATE_STRING);
	printf(" D.J. Bates, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n", __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
}


void * br_malloc(size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does malloc with error checking                        *
 \***************************************************************/
{
	if (size <= 0)
	{ // nothing to allocate
		return NULL;
	}
	else
	{ // try to allocate memory
		void *x = malloc(size);
		if (x == NULL)
		{
			//			raise(SIGINT);
			printf("ERROR: bertini_real's malloc was unable to allocate memory (%d)!\n", (int) size);
			br_exit(ERROR_MEMORY_ALLOCATION);
		}
		return x;
	}
}

void *br_realloc(void *ptr, size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does realloc with error checking                       *
 \***************************************************************/
{
	if (size <= 0)
	{ // nothing to allocate - free memory and return NULL
		free(ptr);
		ptr = NULL;
	}
	else
	{ // try to reallocate memory
		ptr = realloc(ptr, size);
		if (ptr == NULL)
		{
			printf("ERROR: bertini_real's realloc was unable to re-allocate memory!\n");
			br_exit(ERROR_MEMORY_ALLOCATION);
		}
	}
	return ptr;
}





bool is_identity(const mat_d M)
{
	
	
	if (M->rows!=M->cols) {
		return false;
	}
	
	comp_d one;
	set_one_d(one);
	
	comp_d temp;
	for (int ii=0; ii<M->rows; ii++) {
		for (int jj=0; jj<M->cols; jj++) {
			if (ii==jj) {
				sub_d(temp, &M->entry[ii][jj], one);
				if (d_oneNorm_d(temp)>0) {
					return false;
				}
			}
			else{
				if (d_oneNorm_d(&M->entry[ii][jj])>0) {
					return false;
				}
			}
			
		}
	}
	
	return true;
}



bool is_identity(mat_mp M)
{
	
	if (M->rows!=M->cols) {
		return false;
	}
	
	
	comp_mp one; init_mp(one); set_one_mp(one);
	comp_mp temp;  init_mp(temp);
	
	
	for (int ii=0; ii<M->rows; ii++) {
		for (int jj=0; jj<M->cols; jj++) {
			if (ii==jj) {
				sub_mp(temp, &M->entry[ii][jj], one);
				if (d_oneNorm_mp(temp)>0) {
					clear_mp(one); clear_mp(temp);
					return false;
				}
			}
			else{
				if (d_oneNorm_mp(&M->entry[ii][jj])>0) {
					clear_mp(one); clear_mp(temp);
					return false;
				}
			}
			
		}
	}
	
	clear_mp(one); clear_mp(temp);
	
	return true;
}






void norm_of_difference(mpf_t result, const vec_mp left, const vec_mp right)
{
	if (left->size!=right->size || left->size == 0) {
		printf("attempting to take difference of two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		
		print_point_to_screen_matlab(left,"left");
		print_point_to_screen_matlab(right,"right");
		deliberate_segfault();
	}
	
	
	vec_mp difference;  init_vec_mp2(difference, left->size,MIN(left->curr_prec,right->curr_prec));difference->size = left->size;
	comp_mp temp; init_mp2(temp,MIN(left->curr_prec,right->curr_prec));
	
	for (int ii = 0;  ii< left->size; ++ii) {
		sub_mp(&difference->coord[ii], &left->coord[ii], &right->coord[ii]);
	}
	
	twoNormVec_mp(difference, temp);
	
	mpf_abs_mp(result, temp);
	clear_vec_mp(difference);
	clear_mp(temp);
	return;
}


void norm_of_difference_mindim(mpf_t result, const vec_mp left, const vec_mp right)
{
	using std::min;
	auto ls = left->size;
	auto rs = right->size;
	auto m = min(ls,rs);

	vec_mp a, b;
	init_vec_mp(a,m);
	vec_cp_mp(a, left);
	a->size = m;

	init_vec_mp(b,m);
	vec_cp_mp(b, right);
	b->size = m;

	norm_of_difference(result, a, b);
	return;
}


void dehomogenize(point_d *result, const point_d dehom_me)
{
	comp_d denom;
	change_size_vec_d(*result,dehom_me->size-1);
	(*result)->size = dehom_me->size-1;
	set_d(denom, &dehom_me->coord[0]);
	
	for (int ii=0; ii<dehom_me->size-1; ++ii) {
		div_d(&(*result)->coord[ii],&(dehom_me)->coord[ii+1],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	return;
}

void dehomogenize(point_mp *result, const point_mp dehom_me)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		br_exit(977);
	}
	
	comp_mp denom; init_mp2(denom,dehom_me->curr_prec);
	change_size_vec_mp((*result),dehom_me->size-1);
	
	(*result)->size = dehom_me->size-1;
	
	set_mp(denom, &dehom_me->coord[0]);
	for (int ii=0; ii<dehom_me->size-1; ++ii) {
		div_mp(&(*result)->coord[ii],&(dehom_me)->coord[ii+1],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	clear_mp(denom);
	return;
}

void dehomogenize(point_mp *result, const point_mp dehom_me, int num_variables)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		br_exit(977);
	}
	
	comp_mp denom; init_mp2(denom,dehom_me->curr_prec);
	change_size_vec_mp((*result),dehom_me->size-1);
	
	(*result)->size = dehom_me->size-1;
	
	set_mp(denom, &dehom_me->coord[0]);
	for (int ii=0; ii<num_variables-1; ++ii) {
		div_mp(&(*result)->coord[ii],&(dehom_me)->coord[ii+1],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	for (int ii=num_variables-1; ii<dehom_me->size-1; ++ii) {
		set_mp( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
	}
	
	clear_mp(denom);
	return;
}




void dehomogenize(point_d *result, const point_d dehom_me, int num_variables)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		br_exit(977);
	}
	
	comp_d denom;
	change_size_point_d((*result),dehom_me->size-1);
	(*result)->size = dehom_me->size-1;
	
	set_d(denom, &dehom_me->coord[0]);
	for (int ii=0; ii<num_variables-1; ++ii) {
		div_d( &(*result)->coord[ii],&(dehom_me)->coord[ii+1],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	for (int ii=num_variables-1; ii<dehom_me->size-1; ++ii) {
		set_d( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
	}
	
	return;
}


void nonconj_transpose(mat_d Res, const mat_d M)  /* Stores NON-CONJUGATE transpose of M in Res. */
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	int i, j, rows = M->rows, cols = M->cols;
	
	if (Res != M)
	{ // setup Res
		change_size_mat_d(Res, cols, rows);
		
		for (i = 0; i < cols; i++)
			for (j = 0; j < rows; j++)
			{
				set_d(&Res->entry[i][j], &M->entry[j][i]);
			}
	}
	else // Res = M
	{ // need to use a temporary matrix
		mat_d tempMat;
		// copy M to tempMat
		init_mat_d(tempMat, rows, cols);
		mat_cp_d(tempMat, M);
		
		// setup Res
		change_size_mat_d(Res, cols, rows);
		
		for (i = 0; i < cols; i++)
			for (j = 0; j < rows; j++)
			{
				set_d(&Res->entry[i][j], &tempMat->entry[j][i]);
			}
		
		// clear tempMat
		clear_mat_d(tempMat);
	}
	// set the size
	Res->rows = cols;
	Res->cols = rows;
	
	return;
}

void nonconj_transpose(mat_mp Res, const mat_mp M)  /* Stores NON-CONJUGATE transpose of M in Res. */
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	int i, j, rows = M->rows, cols = M->cols;
	
	if (Res != M)
	{ // setup Res
		change_size_mat_mp(Res, cols, rows);
		
		for (i = 0; i < cols; i++)
			for (j = 0; j < rows; j++)
			{
				set_mp(&Res->entry[i][j], &M->entry[j][i]);
			}
	}
	else // Res = M
	{ // need to use a temporary matrix
		mat_mp tempMat;
		// copy M to tempMat
		init_mat_mp2(tempMat, rows, cols, M->curr_prec);
		mat_cp_mp(tempMat, M);
		
		// setup Res
		change_size_mat_mp(Res, cols, rows);
		
		for (i = 0; i < cols; i++)
			for (j = 0; j < rows; j++)
			{
				set_mp(&Res->entry[i][j], &tempMat->entry[j][i]);
			}
		
		// clear tempMat
		clear_mat_mp(tempMat);
	}
	// set the size
	Res->rows = cols;
	Res->cols = rows;
	
	return;
}



void dot_product_d(comp_d result, const vec_d left, const vec_d right)
{
	if (left->size!=right->size) {
		printf("attempting to dot_d two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		br_exit(5901);
	}
	
	set_zero_d(result);
	
	comp_d temp;
	for (int ii=0; ii<left->size; ++ii) {
		mul_d(temp,&left->coord[ii],&right->coord[ii]);
		add_d(result,result,temp);
	}
}

void dot_product_mp(comp_mp result, const vec_mp left, const vec_mp right)
{
	if (left->size!=right->size) {
		printf("attempting to dot_mp two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		br_exit(5902);
	}
	
	set_zero_mp(result);
	comp_mp temp; init_mp2(temp,MIN(left->curr_prec,right->curr_prec));
	for (int ii=0; ii<left->size; ++ii) {
		mul_mp(temp,&left->coord[ii],&right->coord[ii]);
		add_mp(result,result,temp);
	}
	clear_mp(temp);
}




void dot_product_mindim(comp_d result, const vec_d left, const vec_d right)
{
	
	set_zero_d(result);
	comp_d temp;
	for (int ii=0; ii<MIN(left->size,right->size); ++ii) {
		mul_d(temp,&left->coord[ii],&right->coord[ii]);
		add_d(result,result,temp);
	}
	
}


void dot_product_mindim(comp_mp result, const vec_mp left, const vec_mp right)
{
	
	set_zero_mp(result);
	comp_mp temp; init_mp2(temp,MIN(left->curr_prec,right->curr_prec));
	for (int ii=0; ii<MIN(left->size,right->size); ++ii) {
		mul_mp(temp,&left->coord[ii],&right->coord[ii]);
		add_mp(result,result,temp);
	}
	clear_mp(temp);
}





int take_determinant_d(comp_d determinant, const mat_d source_matrix)
{
	
	
	
	if (source_matrix->cols!=source_matrix->rows) {
		printf("source matrix is not square! (%d rows, %d columns)\n",source_matrix->rows,source_matrix->cols);
		exit(-108);
	}
	if (source_matrix->cols==0) {
		printf("source matrix has 0 entries!");
		exit(-109);
	}
	
	int num_variables = source_matrix->cols;
	int ii;
	
	
	mat_d intermediate; init_mat_d(intermediate,0,0);
	
	int *rwnm = NULL;
	vec_d garbage; init_vec_d(garbage,0);
	vec_d zerovec; init_vec_d(zerovec,0); change_size_vec_d(zerovec,num_variables); zerovec->size = num_variables;
	
	
	for (ii=0; ii<num_variables; ii++) {
		set_zero_d(&zerovec->coord[ii]);
	}
	
	
	
	double tol = TOL_DOUBLE_PRECISION; //  these should be for realsies
	double largeChange = LARGECHANGE_DOUBLEPRECISION;
	
	// returns x, intermediate.
	
	//	print_matrix_to_screen_matlab(source_matrix,"detme");
	
	int sign;
	
	int retval = LU_matrixSolve_d(garbage, intermediate, &rwnm, &sign, const_cast<_mat_d*>(source_matrix), zerovec,tol,largeChange);
	//the solution is in intermediate.
	//error check.  solution failed if retval!=0
	if (retval!=0) {
		//		printf("LU decomposition failed (d)\n");
		//		print_matrix_to_screen_matlab(intermediate,"failed_result");
		//		print_matrix_to_screen_matlab(source_matrix,"source_matrix");
		//		deliberate_segfault();
		set_zero_d(determinant);
	}
	else{
		//compute the determinant
		set_one_d(determinant); // initialize
		for (ii=0; ii<num_variables; ii++) {
			mul_d(determinant,determinant,&intermediate->entry[rwnm[ii]][ii]);
		}
		determinant->r = determinant->r*sign;
		determinant->i = determinant->i*sign;
	}
	// this verified correct via 20 samples in matlab.  dab.
	
	free(rwnm);
	clear_vec_d(garbage);
	clear_vec_d(zerovec);
	clear_mat_d(intermediate);
	
	return 0;
}

int take_determinant_mp(comp_mp determinant, const mat_mp source_matrix)
{
	
	
	
	if (source_matrix->cols!=source_matrix->rows) {
		printf("source matrix is not square! (%d rows, %d columns)\n",source_matrix->rows,source_matrix->cols);
		exit(-108);
	}
	if (source_matrix->cols==0) { // no need to check rows -- know they are equal already
		printf("source matrix has 0 entries!");
		exit(-109);
	}
	
	int num_variables = source_matrix->cols;
	
	mat_mp intermediate; init_mat_mp(intermediate,num_variables,num_variables);
	intermediate->rows = intermediate->cols = num_variables;
	vec_mp zerovec; init_vec_mp(zerovec,source_matrix->cols); change_size_vec_mp(zerovec,num_variables); zerovec->size = num_variables;
	vec_mp garbage; init_vec_mp(garbage,source_matrix->cols); garbage->size = source_matrix->cols;
	
	int sign;
	int *rwnm = NULL;
	
	mpf_t tol;  mpfr_init(tol); // = 1e-14
	mpf_t largeChange; mpfr_init(largeChange); //  = 1e11
	
	mpf_set_d(tol, TOL_MP); // tol is the minimum acceptable 1-norm for each row during decomposition
	mpf_set_d(largeChange, LARGECHANGE_MP);
	
	
	
	
	for (int ii=0; ii<num_variables; ii++) {
		set_zero_mp(&zerovec->coord[ii]);
	}
	
	//  these should be for realsies
	
	// returns x, intermediate.
	
	//	print_matrix_to_screen_matlab(tempmat,"tempmat");
	
	
	int retval = LU_matrixSolve_mp(garbage, intermediate, &rwnm, &sign, const_cast<_mat_mp*>(source_matrix), zerovec,tol,largeChange);
	//the solution is in intermediate.
	//error check.  solution failed if retval!=0
	if (retval!=0) {
		//		printf("LU decomposition failed (mp)\n");
		//		print_matrix_to_screen_matlab_mp(intermediate,"failed_result");
		//		print_matrix_to_screen_matlab_mp(source_matrix,"source_matrix");
		//		deliberate_segfault();
		set_zero_mp(determinant);
	}
	else{
		//compute the determinant
		set_one_mp(determinant); // initialize
		for (int ii=0; ii<num_variables; ii++) {
			mul_mp(determinant,determinant,&intermediate->entry[rwnm[ii]][ii]);
		}
		if (sign==-1)
		{
			neg_mp(determinant,determinant);
		}
	}
	//	print_matrix_to_screen_matlab(source_matrix,"detme");
	//	printf("candidate=%lf+1i*%lf;det(detme)\n",determinant->r,determinant->i);
	//	mypause();
	// this verified correct via 20 samples in matlab.  dab.
	
	free(rwnm);
	mpf_clear(tol);
	mpf_clear(largeChange);
	clear_mat_mp(intermediate);
	clear_vec_mp(zerovec);
	clear_vec_mp(garbage);
	return 0;
}



/**
 computes the projection value given a homogeneous input.
 
 double type
 */
void projection_value_homogeneous_input(comp_d result, const vec_d input, const vec_d projection)
{
	set_zero_d(result);
	comp_d temp;
	for (int ii=0; ii<projection->size; ii++) {
		mul_d(temp, &input->coord[ii], &projection->coord[ii]);
		add_d(result, result, temp);
	}
	set_d(temp, result);
	div_d(result, temp, &input->coord[0]);
	
	return;
}

/**
 computes the projection value given a homogeneous input.
 
 mp type
 */
void projection_value_homogeneous_input(comp_mp result, const vec_mp input, const vec_mp projection)
{
	
	
	set_zero_mp(result);
	comp_mp temp; init_mp2(temp,MIN(input->curr_prec,projection->curr_prec));
	for (int ii=0; ii<projection->size; ii++) {
		mul_mp(temp, &input->coord[ii], &projection->coord[ii]);
		add_mp(result, result, temp);
	}
	set_mp(temp, result);
	div_mp(result, temp, &input->coord[0]);
	clear_mp(temp);
	
}



int isSamePoint_inhomogeneous_input(const point_d left, const point_d right, double tolerance){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_inhom_d with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
		//		exit(-287);
	}
	
	
	int indicator = isSamePoint(const_cast<_point_d*>(left),NULL,52,const_cast<_point_d*>(right),NULL,52,tolerance);
	
	
	return indicator;
}


int isSamePoint_inhomogeneous_input(const point_mp left, const point_mp right, double tolerance){
	
	if (left->size!=right->size) {
		std::stringstream ss;
		ss << "attempting to isSamePoint_inhom_mp with disparate sized points.\n";
		ss << "left: " << left->size << "\t right: " << right->size;
		throw std::logic_error(ss.str());
	}
	
	comp_mp temp1, temp2;  init_mp(temp1); init_mp(temp2);
	comp_mp one;  init_mp(one); set_one_mp(one);
	//double infNormVec_mp(vec_mp X)
	for (int ii=0; ii<left->size; ii++) {
		abs_mp(temp1, &left->coord[ii]);
		abs_mp(temp2, &right->coord[ii]); // take the absolute values of the two coordinates.
		
		// compare them.
		if ( mpf_cmp(temp1->r, temp2->r)<0 ) {  // if left smaller than right,
			
			if (mpf_cmp(temp2->r,one->r)>0) { // if the max abs value bigger than 1.
											//scale by temp2
				div_mp(temp1, &left->coord[ii], temp2);
				div_mp(temp2, &right->coord[ii], temp2);
				
				sub_mp(temp1, temp1, temp2);  // temp1 = temp1-temp2
				abs_mp(temp2, temp1);  // temp2 = |temp1|
				
				if (mpf_get_d(temp2->r) > tolerance) {
					clear_mp(temp1); clear_mp(temp2); clear_mp(one);
					return 0;
				}
				
			}
			else{ // no need to scale, because bigger coord was smaller than 1.
				sub_mp(temp1, &left->coord[ii], &right->coord[ii]); // take difference
				abs_mp(temp2, temp1); // absolute value that difference
				if (mpf_get_d(temp2->r)> tolerance) { // if bigger than allowed tolerance, reject as not the same point.
					clear_mp(temp1); clear_mp(temp2); clear_mp(one);
					return 0;
				}
			}
		}
		else {  // if right smaller than left,
			
			if (mpf_cmp(temp1->r,one->r)>0) { // if the larger absolute value is bigger than 1.
											//scale by temp1, cuz it was bigger than temp2 and > 1
				
				div_mp(temp2, &right->coord[ii], temp1);
				div_mp(temp1, &left->coord[ii], temp1);
				
				sub_mp(temp1, temp1, temp2);
				
				abs_mp(temp2, temp1);  // temp2 = |temp1|
				
				if (mpf_get_d(temp2->r) > tolerance) {
					clear_mp(temp1); clear_mp(temp2); clear_mp(one);
					return 0;
				}
			}
			else{
				sub_mp(temp1, &left->coord[ii], &right->coord[ii]); // take difference
				abs_mp(temp2, temp1); // absolute value that difference
				if (mpf_get_d(temp2->r)> tolerance) { // if bigger than allowed tolerance, reject as not the same point.
					clear_mp(temp1); clear_mp(temp2); clear_mp(one);
					return 0;
				}
			}
		}
		
		
	}
	
	clear_mp(temp1); clear_mp(temp2); clear_mp(one);
	return 1; //made it all the way down here, must be the same!!!
	
}



int isSamePoint_homogeneous_input(const point_d left, const point_d right, double tolerance){
	
	if (left->size!=right->size) {
		std::stringstream ss;
		ss << "attempting to isSamePoint_hom_d with disparate sized points.\n";
		ss << "left: " << left->size << "\t right: " << right->size;
		throw std::logic_error(ss.str());
	}
	
	vec_d dehom_left;  init_vec_d(dehom_left,left->size-1);  dehom_left->size = left->size-1;
	vec_d dehom_right; init_vec_d(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left);
	dehomogenize(&dehom_right,right);
	
	int indicator = isSamePoint(dehom_left,NULL,52,dehom_right,NULL,52,tolerance);
	
	clear_vec_d(dehom_left); clear_vec_d(dehom_right);
	
	return indicator;
}


int isSamePoint_homogeneous_input(const point_mp left, const point_mp right, double tolerance){
	
	if (left->size!=right->size) {
		std::stringstream ss;
		ss << "attempting to isSamePoint_hom_mp with disparate sized points.\n";
		ss << "left: " << left->size << "\t right: " << right->size;
		throw std::logic_error(ss.str());
	}
	
	vec_mp dehom_left;  init_vec_mp2(dehom_left,left->size-1,left->curr_prec);  dehom_left->size = left->size-1;
	vec_mp dehom_right; init_vec_mp2(dehom_right,right->size-1,right->curr_prec); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left);
	dehomogenize(&dehom_right,right);
	

	int indicator = isSamePoint_inhomogeneous_input(dehom_left, dehom_right, tolerance);
	
	
	clear_vec_mp(dehom_left); clear_vec_mp(dehom_right);
	
	return indicator;
}



void real_threshold(comp_mp blabla, double threshold)
{
	
	comp_d temp;
	mp_to_d(temp, blabla);
	
   if (fabs(temp->r) < threshold) {
		mpf_set_str( blabla->r, "0.0", 10);
	}
	
	if (fabs(temp->i) < threshold) {
		mpf_set_str( blabla->i, "0.0", 10);
	}
	
	return;
}



void real_threshold(vec_mp blabla, double threshold)
{
	if (blabla->size == 0) {
		return;
	}
	
	comp_d temp;
	for (int ii=0; ii<blabla->size; ii++) {
		mp_to_d(temp, &blabla->coord[ii]);
        
       if (fabs(temp->r) < threshold) {
			mpf_set_str( blabla->coord[ii].r, "0.0", 10);
		}
		
		if (fabs(temp->i) < threshold) {
			mpf_set_str( blabla->coord[ii].i, "0.0", 10);
		}
	}
	return;
}


void real_threshold(mat_mp blabla, double threshold)
{
	if ( (blabla->rows == 0) || (blabla->cols == 0) ) {
		return;
	}
	
	comp_d temp;
	for (int jj=0; jj<blabla->cols; jj++) {
		for (int ii=0; ii<blabla->rows; ii++) {
			mp_to_d(temp, &blabla->entry[ii][jj]);
//            if ( fabs(temp->r) < threshold) {
//				mpf_set_str( blabla->entry[ii][jj].r, "0.0", 10);
//			}
			if ( fabs(temp->i) < threshold) {
				mpf_set_str( blabla->entry[ii][jj].i, "0.0", 10);
			}
		}
	}
	return;
}




void print_point_to_screen_matlab(const vec_d M, std::string name)
{
	
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		printf(" %.8le+1i*%.8le;\n",M->coord[kk].r,M->coord[kk].i);
	}
	printf("];\n\n");
}

void print_point_to_screen_matlab(const vec_mp M, std::string name)
{
	
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		mpf_out_str (NULL, 10, 8, M->coord[kk].r);
		printf("+1i*");
		mpf_out_str (NULL, 10, 8, M->coord[kk].i);
		printf(";\n");
	}
	printf("];\n\n");
}



void print_matrix_to_screen_matlab(const mat_d M, std::string name)
{

	
	printf("%%matrix '%s' has dimensions %dx%d\n", name.c_str(), M->rows,M->cols);
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (int jj = 0; jj < M->cols; jj++)
		{
			printf("%.4le+1i*%.4le ",M->entry[kk][jj].r,M->entry[kk][jj].i );
		}
		if (kk!= M->rows-1) {
			printf(";...\n");
		}
		
	}
	printf("];\n\n");
}
void print_matrix_to_screen_matlab(const mat_mp M, std::string name)
{

	
	printf("%%matrix '%s' has dimensions %dx%d\n",name.c_str(), M->rows,M->cols);
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (int jj = 0; jj < M->cols; jj++)
		{
			
			mpf_out_str (NULL, 10, 8, M->entry[kk][jj].r);
			printf("+1i*");
			mpf_out_str (NULL, 10, 8, M->entry[kk][jj].i); // base 10 , 7 digits
			printf("\t");
		}
		printf(";\n");
	}
	printf("];\n\n");
}


void print_comp_matlab(const comp_mp M, std::string name){
	printf("%s=",name.c_str());
	mpf_out_str (NULL, 10, 8, M->r);
	printf("+1i*");
	mpf_out_str (NULL, 10, 8, M->i); // base 10, 6 digits
	printf("\n");
	return;
}

void print_comp_matlab(const comp_d M, std::string name){
	printf("%s=%.5le+1i*%.5le\n",name.c_str(),M->r,M->i);
	return;
}

void print_path_retVal_message(int retVal){
	
	
	if (retVal==100) {
		printf("max_prec_reached\n");
	}
	if (retVal==-50) {
		printf("reached_minTrackT\nrelevant setting name is 'NbhdRadius'\n");
	}
	else if (retVal==-200){
		printf("cycle_num_too_high\n");
	}
	else if (retVal==-15){
		printf("PSEG_failed\n");
	}
	else if (retVal==-2){
		printf("going_to_infinity\n");
	}
	else if (retVal==-4){
		printf("security_max\n");
	}
	else if (retVal==-3){
		printf("step_size_too_small\n");
	}
	else if (retVal==-10){
		printf("too_many_steps\n");
	}
	else if (retVal==-20){
		printf("refining_failed\n");
	}
	else if (retVal==-100){
		printf("higher_prec_needed\n");
	}
	else if (retVal==-99){
		printf("retVal_NAN\n");
	}
	else if (retVal==-98){
		printf("retVal_Bertini_Junk\n");
	}
	else if (retVal==-97){
		printf("Failed_to_converge (used in newton iterations)\n");
	}
	else if (retVal==-22){
		printf("sharpening_singular_endpoint\nthis is used when the sharpening sees that the endpoint is singular and cannot sharpen it\n");
	}
	else if (retVal==-21){
		printf("sharpening_failed\nthis is used when the sharpening of an endpoint does not reach the desired tolerance\n");
	}
	else if (retVal==-22){
		printf("higher_dim\nthis is used in regeneration when an endpoint lies on a higher dimensional component\n");
	}
	
	
	return;
}



int get_num_vars_PPD(const preproc_data PPD){
	int num_vars = 0; // initialize
	

	
	//run through each variable group
	for (int ii=0; ii<(PPD.num_var_gp+PPD.num_hom_var_gp); ii++) {
		num_vars+= PPD.size[ii];
		num_vars+= PPD.type[ii];
	}
	
	return num_vars;
	
}



void cp_patch_mp(patch_eval_data_mp *PED, const patch_eval_data_mp PED_input)
{
	PED->num_patches = PED_input.num_patches;
	
	// set the current precision
	PED->curr_prec = PED_input.curr_prec;
	
	// initialize patchCoeff to this preicision
	init_mat_mp2(PED->patchCoeff, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols, PED->curr_prec);
	init_mat_rat(PED->patchCoeff_rat, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols);
	
	// setup patchCoeff
	mat_cp_mp_rat(PED->patchCoeff, PED->patchCoeff_rat, PED_input.patchCoeff, PED_input.patchCoeff_rat);
	
	return;
}


void cp_patch_d(patch_eval_data_d *PED, const patch_eval_data_d PED_input)
{
	PED->num_patches = PED_input.num_patches;
	
	
	// initialize patchCoeff to this preicision
	init_mat_d(PED->patchCoeff, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols);
	
	// setup patchCoeff
	mat_cp_d(PED->patchCoeff, PED_input.patchCoeff);
	
	return;
}

void cp_preproc_data(preproc_data *PPD, const preproc_data & PPD_input)
{
	
	PPD->num_funcs = PPD_input.num_funcs;
	PPD->num_hom_var_gp = PPD_input.num_hom_var_gp;
	PPD->num_var_gp = PPD_input.num_var_gp;
	
	int total_gp = PPD->num_hom_var_gp + PPD->num_var_gp;
	
	
	
	if (total_gp==0) {
		return;
	}
	else{
		PPD->size = (int *)br_malloc(total_gp * sizeof(int));
		PPD->type = (int *)br_malloc(total_gp * sizeof(int));
		
		for (int i = 0; i < total_gp; i++)
		{
			PPD->size[i] = PPD_input.size[i];
			PPD->type[i] = PPD_input.type[i];
		}
	}
	
	
	
	return;
}


void clear_post_process_t(post_process_t * endPoint, int num_vars)
{
	
	if (endPoint->sol_prec >1)
	{
		if (endPoint->sol_prec >= 64)
		{ // clear _mp
			mpf_clear(endPoint->function_resid_mp);
			mpf_clear(endPoint->newton_resid_mp);
			for (int j = 0; j < num_vars; j++)
			{
				clear_mp(endPoint->sol_mp[j]);
			}
		}
		free(endPoint->sol_mp);
		free(endPoint->sol_d);
	}
}


void print_tracker(const tracker_config_t * T)
{
    std::cout << "BEGIN TRACKER CONFIG:\n\n" << std::endl;
	std::cout << "numVars: " << T->numVars << std::endl;
	std::cout << "numPathVars: " << T->numPathVars << std::endl;
	std::cout << "numParams: " << T->numParams << std::endl;
	std::cout << "numFuncs: " << T->numFuncs << std::endl;
	
	std::cout << "maxStepSize: " << T->maxStepSize << std::endl;
	std::cout << "minStepSizeBeforeEndGame: " << T->minStepSizeBeforeEndGame << std::endl;
	std::cout << "minStepSizeDuringEndGame: " << T->minStepSizeDuringEndGame << std::endl;
	std::cout << "minStepSize: " << T->minStepSize << std::endl;
	std::cout << "currentStepSize: " << T->currentStepSize << std::endl;
	std::cout << "first_step_of_path: " << T->first_step_of_path << std::endl;
	std::cout << "minTrackT: " << T->minTrackT << std::endl;
	
	std::cout << "basicNewtonTol: " << T->basicNewtonTol << std::endl;
	std::cout << "endgameNewtonTol: " << T->endgameNewtonTol << std::endl;
	std::cout << "final_tolerance: " << T->final_tolerance << std::endl;
	
	std::cout << "cSecInc: " << T->cSecInc << std::endl;
	std::cout << "maxNewtonIts: " << T->maxNewtonIts << std::endl;
	std::cout << "MPType: " << T->MPType << std::endl;
	std::cout << "Precision: " << T->Precision << std::endl;
	std::cout << "outputLevel: " << T->outputLevel << std::endl;
	std::cout << "screenOut: " << T->screenOut << std::endl;
	
	std::cout << "targetT: " << T->targetT << std::endl;
	std::cout << "endgameBoundary: " << T->endgameBoundary << std::endl;
	std::cout << "endgameSwitch: " << T->endgameSwitch << std::endl;
	
	std::cout << "goingToInfinity: " << T->goingToInfinity << std::endl;
	std::cout << "maxNumSteps: " << T->maxNumSteps << std::endl;
	std::cout << "endgameNumber: " << T->endgameNumber << std::endl;
	
	std::cout << "latest_cond_num_exp: " << T->latest_cond_num_exp << std::endl;
	std::cout << "steps_since_last_CN: " << T->steps_since_last_CN << std::endl;
	
	std::cout << "power_series_sample_factor: " << T->power_series_sample_factor << std::endl;
	std::cout << "cycle_num_max: " << T->cycle_num_max << std::endl;
	std::cout << "num_PSEG_sample_points: " << T->num_PSEG_sample_points << std::endl;
	
	std::cout << "latest_newton_residual_d: " << T->latest_newton_residual_d << std::endl;
	
	
	std::cout << "t_val_at_latest_sample_point: " << T->t_val_at_latest_sample_point << std::endl;
	std::cout << "error_at_latest_sample_point: " << T->error_at_latest_sample_point << std::endl;
	std::cout << "final_tolerance: " << T->final_tolerance << std::endl;
	
	std::cout << "real_threshold: " << T->real_threshold << std::endl;
	std::cout << "endgameOnly: " << T->endgameOnly << std::endl;
	
	std::cout << "AMP_bound_on_abs_vals_of_coeffs: " << T->AMP_bound_on_abs_vals_of_coeffs << std::endl;
	std::cout << "AMP_bound_on_degree: " << T->AMP_bound_on_degree << std::endl;
	std::cout << "AMP_eps: " << T->AMP_eps << std::endl;
	std::cout << "AMP_Phi: " << T->AMP_Phi << std::endl;
	std::cout << "AMP_Psi: " << T->AMP_Psi << std::endl;
	
	std::cout << "AMP_safety_digits_1: " << T->AMP_safety_digits_1 << std::endl;
	std::cout << "AMP_safety_digits_2: " << T->AMP_safety_digits_2 << std::endl;
	std::cout << "AMP_max_prec: " << T->AMP_max_prec << std::endl;
	
	std::cout << "sing_val_zero_tol: " << T->sing_val_zero_tol << std::endl;
	std::cout << "cond_num_threshold: " << T->cond_num_threshold << std::endl;
	
	std::cout << "step_fail_factor: " << T->step_fail_factor << std::endl;
	std::cout << "step_success_factor: " << T->step_success_factor << std::endl;
	
	std::cout << "max_num_pts_for_trace: " << T->max_num_pts_for_trace << std::endl;
	std::cout << "max_num_mon_linears: " << T->max_num_mon_linears << std::endl;
	std::cout << "max_num_bad_loops_in_mon: " << T->max_num_bad_loops_in_mon << std::endl;
	
	std::cout << "final_tol_multiplier: " << T->final_tol_multiplier << std::endl;
	std::cout << "final_tol_times_mult: " << T->final_tol_times_mult << std::endl;
	
	std::cout << "sharpenDigits: " << T->sharpenDigits << std::endl;
	std::cout << "sharpenOnly: " << T->sharpenOnly << std::endl;
	
	std::cout << "regen_remove_inf: " << T->regen_remove_inf << std::endl;
	std::cout << "regen_higher_dim_check: " << T->regen_higher_dim_check << std::endl;
	std::cout << "sliceBasicNewtonTol: " << T->sliceBasicNewtonTol << std::endl;
	std::cout << "sliceEndgameNewtonTol: " << T->sliceEndgameNewtonTol << std::endl;
	std::cout << "sliceFinalTol: " << T->sliceFinalTol << std::endl;
	
	std::cout << "minCycleTrackBack: " << T->minCycleTrackBack << std::endl;
	std::cout << "junkRemovalTest: " << T->junkRemovalTest << std::endl;
	std::cout << "maxDepthLDT: " << T->maxDepthLDT << std::endl;
	std::cout << "odePredictor: " << T->odePredictor << std::endl;
	
	std::cout << "securityLevel: " << T->securityLevel << std::endl;
	std::cout << "securityMaxNorm: " << T->securityMaxNorm << std::endl;
	
	std::cout << "cutoffCycleTime: " << T->cutoffCycleTime << std::endl;
	std::cout << "cutoffRatioTime: " << T->cutoffRatioTime << std::endl;
	std::cout << "finiteThreshold: " << T->finiteThreshold << std::endl;
	std::cout << "funcResTol: " << T->funcResTol << std::endl;
	std::cout << "ratioTol: " << T->ratioTol << std::endl;
	std::cout << "maxStepsBeforeNewton: " << T->maxStepsBeforeNewton << std::endl;
    
    std::cout << "END TRACKER CONFIG" << std::endl << std::endl << std::endl;
	
	return;
}




//TODO this sort should be optimized.  it is sloppy and wasteful right now.
int sort_increasing_by_real(vec_mp projections_sorted, 
							std::vector< int > & index_tracker, 
							const vec_mp projections_input,
							double distinct_thresh){
	
	
	if (projections_input->size == 0) {
		change_size_vec_mp(projections_sorted,1);
		projections_sorted->size = 0;
		return -1;
	}
	
	
	
	for (int ii=0; ii<projections_input->size; ii++) {
		if (!(mpfr_number_p(projections_input->coord[ii].r) && mpfr_number_p(projections_input->coord[ii].i))) {
			std::cout << "there was NAN in the projections to sort :(" << std::endl;
			print_point_to_screen_matlab(projections_input, "projections_input");
			
			return -51;
		}
	}
	
	
	
	
	
	std::vector< int > index_tracker_non_unique;
	std::vector< double > projvals_as_doubles;
	
	
	vec_mp projections_sorted_non_unique;
	init_vec_mp2(projections_sorted_non_unique,projections_input->size,1024);
	projections_sorted_non_unique->size = projections_input->size;
	
	
	
	std::set<int> unsorted_indices;
	for (int ii=0; ii<projections_input->size; ii++) {
		unsorted_indices.insert(ii);
	}
	
	
	
	
	//sort by size
	for (int ii=0; ii<projections_input->size; ii++) { // for each of the projection values input
		double min = 1e20; // reset this bogus value
		
		int indicator = -1;
		
		// this loop finds the minimum projection value
		for (std::set<int>::iterator set_iter = unsorted_indices.begin(); set_iter!=unsorted_indices.end(); set_iter++) {
			
			double curr = mpf_get_d(projections_input->coord[*set_iter].r); // convert projection value to a double for comparison
			if ( curr < min) { // compare
				indicator = *set_iter;
				min = curr;
			}
		}
		
		if (indicator==-1) { // if min value was larger than a huge number
			printf("min projection value was *insanely* large\n");
		}
		
		unsorted_indices.erase(indicator);
		
		projvals_as_doubles.push_back(min);
		index_tracker_non_unique.push_back(indicator);
		set_mp( &projections_sorted_non_unique->coord[ii],&projections_input->coord[indicator]);
	}
	
	
	
	
	// filter for uniqueness
	
	change_size_vec_mp(projections_sorted,1); projections_sorted->size = 1;
	
	index_tracker.push_back(index_tracker_non_unique[0]);
	set_mp(&projections_sorted->coord[0],&projections_sorted_non_unique->coord[0])
	int unique_counter = 1;
	for (int ii=1; ii<projections_input->size; ii++) {
		
		double compare_me_lower, compare_me_upper;
		if (fabs(projvals_as_doubles[ii-1])>1) {
			compare_me_lower = 1;
			compare_me_upper = projvals_as_doubles[ii]/projvals_as_doubles[ii-1];
		}
		else{
			compare_me_lower = projvals_as_doubles[ii-1];
			compare_me_upper = projvals_as_doubles[ii];
		}
		
		
		if ( fabs(compare_me_upper-compare_me_lower) < distinct_thresh) { //fabs( projvals_as_doubles[ii-1]-projvals_as_doubles[ii])
			continue;
		}
		else
		{
			increase_size_vec_mp(projections_sorted,unique_counter+1); projections_sorted->size = unique_counter+1;
			set_mp(&projections_sorted->coord[unique_counter],&projections_sorted_non_unique->coord[ii]);
			unique_counter++;
			
			index_tracker.push_back(index_tracker_non_unique[ii]);
		}
	}
	
	
	clear_vec_mp(projections_sorted_non_unique);
	
	return 0;
}



void send_patch_d(const patch_eval_data_d * patch)
{
	comp_d *patch_coeff = NULL;
	patch_eval_data_d_int PED_int;
	MPI_Datatype mpi_comp_d, mpi_patch_d_int;
	
	// setup mpi_comp_d & mpi_patch_d_int
	create_comp_d(&mpi_comp_d);
	create_patch_eval_data_d_int(&mpi_patch_d_int);
	// setup PED_int
	cp_patch_d_int(&PED_int, const_cast<patch_eval_data_d *>(patch), &patch_coeff, 0);
	
	// broadcast patch structures
	MPI_Bcast(&PED_int, 1, mpi_patch_d_int, 0, MPI_COMM_WORLD);
	MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, 0, MPI_COMM_WORLD);
	
	// free memory
	MPI_Type_free(&mpi_comp_d);
	MPI_Type_free(&mpi_patch_d_int);
	free(patch_coeff);
}


void receive_patch_d(patch_eval_data_d * patch)
{
	comp_d *patch_coeff = NULL;
	patch_eval_data_d_int PED_int;
	MPI_Datatype mpi_comp_d, mpi_patch_d_int;
	
	// setup mpi_comp_d & mpi_patch_d_int
	create_comp_d(&mpi_comp_d);
	create_patch_eval_data_d_int(&mpi_patch_d_int);
	
	// recv patch structures
	MPI_Bcast(&PED_int, 1, mpi_patch_d_int, 0, MPI_COMM_WORLD);
	// setup patch_coeff
	patch_coeff = (comp_d *)br_malloc(PED_int.patchCoeff_rows * PED_int.patchCoeff_cols * sizeof(comp_d));
	MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, 0, MPI_COMM_WORLD);
	
	
	// setup patch
	cp_patch_d_int(patch, &PED_int, &patch_coeff, 1);  // patch_coeff is freed in here
	
	// free mpi_comp_d & mpi_patch_d_int
	MPI_Type_free(&mpi_comp_d);
	MPI_Type_free(&mpi_patch_d_int);
}


void send_patch_mp(const patch_eval_data_mp * patch)
{
	char *patchStr = NULL;
	patch_eval_data_mp_int PED_int;
	MPI_Datatype mpi_patch_int;
	
	// setup mpi_patch_int
	create_patch_eval_data_mp_int(&mpi_patch_int);
	// setup PED_int
	cp_patch_mp_int(&PED_int, const_cast<patch_eval_data_mp *>(patch), &patchStr, 0, 0);
	
	// send PED_int
	MPI_Bcast(&PED_int, 1, mpi_patch_int, 0, MPI_COMM_WORLD);
	// send patchStr
	MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// clear memory
	free(patchStr);
	MPI_Type_free(&mpi_patch_int);
}

void receive_patch_mp(patch_eval_data_mp * patch)
{
	char *patchStr = NULL;
	patch_eval_data_mp_int PED_int;
	MPI_Datatype mpi_patch_int;
	
	// setup mpi_patch_int
	create_patch_eval_data_mp_int(&mpi_patch_int);
	// recv PED_int
	MPI_Bcast(&PED_int, 1, mpi_patch_int, 0, MPI_COMM_WORLD);
	
	// setup patchStr
	patchStr = (char *)br_malloc(PED_int.totalLength * sizeof(char));
	// recv patchStr
	MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// setup _mp patch
	cp_patch_mp_int(patch, &PED_int, &patchStr, 1, 1);
	// last 3 arguments:   ,    ,   intype -- 0 if sender, !0, else
	
	
	// free mpi_patch_int (patchStr is freed in cp_patch_mp_int)
	MPI_Type_free(&mpi_patch_int);
}



void send_preproc_data(const preproc_data *PPD){
	
	
	int *buffer = new int[3];
	
	buffer[0] = PPD->num_funcs;
	buffer[1] = PPD->num_hom_var_gp;
	buffer[2] = PPD->num_var_gp;
	
	MPI_Bcast(buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);
	
	int size = PPD->num_hom_var_gp + PPD->num_var_gp;
	MPI_Bcast(PPD->type, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(PPD->size, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	delete [] buffer;
}

void receive_preproc_data(preproc_data *PPD){
	
	
	int *buffer = new int[3];
	
	MPI_Bcast(buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);
	
	PPD->num_funcs = buffer[0];
	PPD->num_hom_var_gp = buffer[1];
	PPD->num_var_gp = buffer[2];
	
	
	
	int num_groups = PPD->num_hom_var_gp + PPD->num_var_gp;
	
	PPD->type = (int *) br_malloc(num_groups * sizeof(int));
	PPD->size = (int *) br_malloc(num_groups * sizeof(int));
	MPI_Bcast(PPD->type, num_groups, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(PPD->size, num_groups, MPI_INT, 0, MPI_COMM_WORLD);
	
	delete [] buffer;
}











void send_mat_d(const mat_d A, int target)
{
    int num_entries;
    MPI_Datatype mpi_mat_d_int, mpi_comp_d;
    mat_d_int A_int;
    comp_d *entries = NULL;
    
    // create the datatypes mpi_mat_d_int & mpi_comp_d
    create_mat_d_int(&mpi_mat_d_int);
    create_comp_d(&mpi_comp_d);
    
    // setup A_int and entries
    cp_mat_d_int(&A_int, const_cast<_mat_d*>(A), &entries, 0);
    num_entries = A->rows * A->cols;
    
    // send A_int
    MPI_Send(&A_int, 1, mpi_mat_d_int, target, MAT_D, MPI_COMM_WORLD);
    // send entries
    MPI_Send(entries, num_entries, mpi_comp_d, target, MAT_D, MPI_COMM_WORLD);
    
    // clear entries
    free(entries);
    
    
    // clear mpi_mat_d_int & mpi_comp_d
    MPI_Type_free(&mpi_mat_d_int);
    MPI_Type_free(&mpi_comp_d);
    
    return;
}
void receive_mat_d(mat_d A, int source)
{
    MPI_Status statty_mc_gatty;
    int num_entries;
    MPI_Datatype mpi_mat_d_int, mpi_comp_d;
    mat_d_int A_int;
    comp_d *entries = NULL;
    
    // create the datatypes mpi_mat_d_int & mpi_comp_d
    create_mat_d_int(&mpi_mat_d_int);
    create_comp_d(&mpi_comp_d);
    
    // recv A_int
    MPI_Recv(&A_int, 1, mpi_mat_d_int, source, MAT_D, MPI_COMM_WORLD, &statty_mc_gatty);
    // setup A and entries
    init_mat_d(A, A_int.rows, A_int.cols);
    
    num_entries = A_int.rows * A_int.cols;
    entries = (comp_d *)bmalloc(num_entries * sizeof(comp_d));
    // recv entries
    MPI_Recv(entries, num_entries, mpi_comp_d, source, MAT_D, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup A
    cp_mat_d_int(A, &A_int, &entries, 1);
    
    
    // clear mpi_mat_d_int & mpi_comp_d
    MPI_Type_free(&mpi_mat_d_int);
    MPI_Type_free(&mpi_comp_d);
    
    return;
}



void send_mat_mp(const mat_mp A, int target)
{
    MPI_Datatype mpi_mat_mp_int;
    mat_mp_int A_int;
    char *Astr = NULL;
    
    // create the datatypes mpi_mat_mp_int
    create_mat_mp_int(&mpi_mat_mp_int);
    
    
    cp_mat_mp_int(&A_int, const_cast<_mat_mp*>(A), &Astr, 1, 0);
    
    // send A_int and Astr
    MPI_Send(&A_int, 1, mpi_mat_mp_int, target, MAT_MP, MPI_COMM_WORLD);
    MPI_Send(Astr, A_int.totalLength, MPI_CHAR, target, MAT_MP, MPI_COMM_WORLD);
    
    // clear Astr
    free(Astr);
    
    
    // clear mpi_mat_mp_int
    MPI_Type_free(&mpi_mat_mp_int);
    
    return;
}
void receive_mat_mp(mat_mp A, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_mat_mp_int;
    mat_mp_int A_int;
    char *Astr = NULL;
    
    // create the datatypes mpi_mat_mp_int
    create_mat_mp_int(&mpi_mat_mp_int);
    
    // recv A_int and Astr
    MPI_Recv(&A_int, 1, mpi_mat_mp_int, source, MAT_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    Astr = (char *)bmalloc(A_int.totalLength * sizeof(char));
    MPI_Recv(Astr, A_int.totalLength, MPI_CHAR, source, MAT_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup A and clear Astr
    cp_mat_mp_int(A, &A_int, &Astr, 1, 1);
    
    
    // clear mpi_mat_mp_int
    MPI_Type_free(&mpi_mat_mp_int);
    
    return;
}





void send_mat_rat(const mat_d A_d, const mat_mp A_mp, const mpq_t ***A_rat, int target)
{
    MPI_Datatype mpi_mat_rat;
    mat_rat_int A_int;
    int rows, cols;
    char *ratStr = NULL;
    
    // create the datatype mpi_mat_rat
    create_mat_rat_int(&mpi_mat_rat);
    
    // setup A_int & ratStr
    rows = A_d->rows;
    cols = A_d->cols;
    cp_mat_rat_int(&A_int, A_rat, &ratStr, rows, cols, 1, 0);
    
    // send A_int & ratStr
    MPI_Send(&A_int, 1, mpi_mat_rat, target, MAT_RAT, MPI_COMM_WORLD);
    MPI_Send(ratStr, A_int.totalLength, MPI_CHAR, target, MAT_RAT, MPI_COMM_WORLD);
    
    // clear ratStr
    free(ratStr);
    
    
    // clear mpi_mat_rat
    MPI_Type_free(&mpi_mat_rat);
    
    return;
}
void receive_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_mat_rat;
    mat_rat_int A_int;
    int rows, cols;
    char *ratStr = NULL;
    
    // create the datatype mpi_mat_rat
    create_mat_rat_int(&mpi_mat_rat);
    
    // recv A_int & ratStr
    MPI_Recv(&A_int, 1, mpi_mat_rat, source, MAT_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    ratStr = (char *)bmalloc(A_int.totalLength * sizeof(char));
    MPI_Recv(ratStr, A_int.totalLength, MPI_CHAR, source, MAT_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup A_rat and clear ratStr
    cp_mat_rat_int(A_rat, &A_int, &ratStr, A_int.rows, A_int.cols, 1, 1);
    
    // setup A_d & A_mp
    for (rows = 0; rows < A_int.rows; rows++)
        for (cols = 0; cols < A_int.cols; cols++)
        {
            mpf_set_q(A_mp->entry[rows][cols].r, A_rat[rows][cols][0]);
            mpf_set_q(A_mp->entry[rows][cols].i, A_rat[rows][cols][1]);
            A_d->entry[rows][cols].r = mpq_get_d(A_rat[rows][cols][0]);
            A_d->entry[rows][cols].i = mpq_get_d(A_rat[rows][cols][1]);
        }
    
    
    // clear mpi_mat_rat
    MPI_Type_free(&mpi_mat_rat);
    
    return;
}




void send_vec_d(const vec_d b, int target)
{
    MPI_Datatype mpi_point_d_int, mpi_comp_d;
    point_d_int b_int;
    comp_d *entries = NULL;
    
    // create the datatype mpi_point_d_int & mpi_comp_d
    create_point_d_int(&mpi_point_d_int);
    create_comp_d(&mpi_comp_d);
    
    // setup b_int and entries
    cp_point_d_int(&b_int, const_cast<_point_d*>(b), &entries, 0, 0, 0);
    
    // send b_int
    MPI_Send(&b_int, 1, mpi_point_d_int, target, VEC_D, MPI_COMM_WORLD);
    // send entries
    MPI_Send(entries, b_int.size, mpi_comp_d, target, VEC_D, MPI_COMM_WORLD);
    
    // clear entries
    free(entries);
    
    
    // clear mpi_point_d_int & mpi_comp_d
    MPI_Type_free(&mpi_point_d_int);
    MPI_Type_free(&mpi_comp_d);
    
    return;
}
void receive_vec_d(vec_d b, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_point_d_int, mpi_comp_d;
    point_d_int b_int;
    comp_d *entries = NULL;
    
    // create the datatype mpi_point_d_int & mpi_comp_d
    create_point_d_int(&mpi_point_d_int);
    create_comp_d(&mpi_comp_d);
    
    // recv b_int
    MPI_Recv(&b_int, 1, mpi_point_d_int, source, VEC_D, MPI_COMM_WORLD, &statty_mc_gatty);
    
    entries = (comp_d *)bmalloc(b_int.size * sizeof(comp_d));
    // recv entries
    MPI_Recv(entries, b_int.size, mpi_comp_d, source, VEC_D, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup b
    cp_point_d_int(b, &b_int, &entries, 1, 0, 1);
    
    
    // clear mpi_point_d_int & mpi_comp_d
    MPI_Type_free(&mpi_point_d_int);
    MPI_Type_free(&mpi_comp_d);
    
    return;
}




void send_vec_mp(const vec_mp b, int target)
{
    MPI_Datatype mpi_vec_mp_int;
    point_mp_int b_int;
    char *bstr = NULL;
    
    // create the datatypes mpi_vec_mp_int
    create_point_mp_int(&mpi_vec_mp_int);
    
    // setup b_int and bstr
    cp_point_mp_int(&b_int, const_cast<_point_mp*>(b), &bstr, 0, 0, 0);
    
    // send b_int and bstr
    MPI_Send(&b_int, 1, mpi_vec_mp_int, target, VEC_MP, MPI_COMM_WORLD);
    MPI_Send(bstr, b_int.totalLength, MPI_CHAR, target, VEC_MP, MPI_COMM_WORLD);
    
    // clear bstr
    free(bstr);
    
    
    // clear mpi_vec_mp_int
    MPI_Type_free(&mpi_vec_mp_int);
    
    return;
}
void receive_vec_mp(vec_mp b, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_vec_mp_int;
    point_mp_int b_int;
    char *bstr = NULL;
    
    // create the datatypes mpi_vec_mp_int
    create_point_mp_int(&mpi_vec_mp_int);
    
    // recv b_int and bstr
    MPI_Recv(&b_int, 1, mpi_vec_mp_int, source, VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    bstr = (char *)bmalloc(b_int.totalLength * sizeof(char));
    MPI_Recv(bstr, b_int.totalLength, MPI_CHAR, source, VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup b and clear bstr
    cp_point_mp_int(b, &b_int, &bstr, 1, 0, 1); // first 1 indicates to free bstr.  second is init_point.  third is intype.  see doc in copy_functions.c
    
    
    // clear mpi_vec_mp_int
    MPI_Type_free(&mpi_vec_mp_int);
    
    return;
}





void send_vec_rat(const mpq_t ***b, int size, int target)
{
    MPI_Datatype mpi_point_rat;
    point_rat_int b_int;
    char *ratStr = NULL;
    
    // create the datatype mpi_point_rat
    create_point_rat_int(&mpi_point_rat);
    
    // setup b_int & ratStr
    b_int.size = size;
    cp_vec_rat_char(&ratStr, b, &b_int.totalLength, size, 0, 0);
    
    // send b_int
    MPI_Send(&b_int, 1, mpi_point_rat, target, VEC_RAT, MPI_COMM_WORLD);
    
    // send ratStr
    MPI_Send(ratStr, b_int.totalLength, MPI_CHAR, target, VEC_RAT, MPI_COMM_WORLD);
    
    // clear ratStr
    free(ratStr);
    
    
    // clear mpi_point_rat
    MPI_Type_free(&mpi_point_rat);
    
    return;
}
void receive_vec_rat(mpq_t ***b, int size, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_point_rat;
    point_rat_int b_int;
    char *ratStr = NULL;
    
    // create the datatype mpi_point_rat
    create_point_rat_int(&mpi_point_rat);
    
    // recv b_int
    MPI_Recv(&b_int, 1, mpi_point_rat, source, VEC_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup & recv ratStr
    ratStr = (char *)bmalloc(b_int.totalLength * sizeof(char));
    MPI_Recv(ratStr, b_int.totalLength, MPI_CHAR, source, VEC_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup b - clears all structures
    cp_vec_rat_char(b, &ratStr, &b_int.totalLength, b_int.size, 1, 1);
    
    
    // clear mpi_point_rat
    MPI_Type_free(&mpi_point_rat);
    
    return;
}











void send_comp_d(const comp_d c, int target)
{
    MPI_Datatype mpi_comp_d;
    create_comp_d(&mpi_comp_d);
    
    MPI_Send(c, 1, mpi_comp_d, target, COMP_D, MPI_COMM_WORLD);
    
    MPI_Type_free(&mpi_comp_d);
    
    return;
}
void receive_comp_d(comp_d c, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_comp_d;
    create_comp_d(&mpi_comp_d);
    
    MPI_Recv(c, 1, mpi_comp_d, source, COMP_D, MPI_COMM_WORLD, &statty_mc_gatty);
    
    MPI_Type_free(&mpi_comp_d);
    
    return;
}




void send_comp_num_d(const comp_d *c, int num, int target)
{
    MPI_Datatype mpi_comp_d;
    create_comp_d(&mpi_comp_d);
    
    MPI_Send(c, num, mpi_comp_d, target, COMP_D, MPI_COMM_WORLD);
    
    MPI_Type_free(&mpi_comp_d);
    
    return;
}
void receive_comp_num_d(comp_d *c, int num, int source)
{
    MPI_Status statty_mc_gatty;
    MPI_Datatype mpi_comp_d;
    create_comp_d(&mpi_comp_d);
    
    MPI_Recv(c, num, mpi_comp_d, source, COMP_D, MPI_COMM_WORLD, &statty_mc_gatty);
    
    MPI_Type_free(&mpi_comp_d);
    
    return;
}




void send_comp_mp(const comp_mp c, int target)
{
    char *str = NULL;
    comp_mp_int c_int;
    MPI_Datatype mpi_comp_mp_int;
    create_comp_mp_int(&mpi_comp_mp_int);
    
    // send data
    cp_comp_mp_int(&c_int, const_cast<_comp_mp*>(c), &str, 0, 0);
    // send c_int
    MPI_Send(&c_int, 1, mpi_comp_mp_int, target, COMP_MP, MPI_COMM_WORLD);
    // send str
    MPI_Send(str, c_int.totalLength, MPI_CHAR, target, COMP_MP, MPI_COMM_WORLD);
    
    // clear str
    free(str);
    
    
    MPI_Type_free(&mpi_comp_mp_int);
    
    return;
}
void receive_comp_mp(comp_mp c, int source)
{
    MPI_Status statty_mc_gatty;
    char *str = NULL;
    comp_mp_int c_int;
    MPI_Datatype mpi_comp_mp_int;
    create_comp_mp_int(&mpi_comp_mp_int);
    
    // recv data
    MPI_Recv(&c_int, 1, mpi_comp_mp_int, source, COMP_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    // setup & recv str
    str = (char *)bmalloc(c_int.totalLength * sizeof(char));
    MPI_Recv(str, c_int.totalLength, MPI_CHAR, source, COMP_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup c
    cp_comp_mp_int(c, &c_int, &str, 1, 1);
    
    
    MPI_Type_free(&mpi_comp_mp_int);
    
    return;
}





void send_comp_num_mp(const comp_mp *c, int num, int target)
{
    int i, j, total = 0, currLoc = 0;
    comp_mp_int *c_int = (comp_mp_int *)bmalloc(num * sizeof(comp_mp_int));
    char *str = NULL, *tempStr = NULL;
    MPI_Datatype mpi_comp_mp_int;
    create_comp_mp_int(&mpi_comp_mp_int);
    
    // send data
    for (i = 0; i < num; i++)
    { // setup c_int[i]
        cp_comp_mp_int(&c_int[i], const_cast<_comp_mp*>(c[i]), &tempStr, 0, 0);
        // update
        total += c_int[i].totalLength;
        str = (char *)brealloc(str, total * sizeof(char));
        for (j = 0; j < c_int[i].totalLength; j++)
        {
            str[currLoc] = tempStr[j];
            currLoc++;
        }
        free(tempStr);
    }
    // send c_int
    MPI_Send(c_int, num, mpi_comp_mp_int, target, COMP_MP, MPI_COMM_WORLD);
    // send str
    MPI_Send(str, total, MPI_CHAR, target, COMP_MP, MPI_COMM_WORLD);
    
    // clear data
    free(str);
    tempStr = NULL;
    
    
    // free data
    free(c_int);
    MPI_Type_free(&mpi_comp_mp_int);
    
    return;
}
void receive_comp_num_mp(comp_mp *c, int num, int source)
{
    MPI_Status statty_mc_gatty;
    int i, total = 0, currLoc = 0;
    comp_mp_int *c_int = (comp_mp_int *)bmalloc(num * sizeof(comp_mp_int));
    char *str = NULL, *tempStr = NULL;
    MPI_Datatype mpi_comp_mp_int;
    create_comp_mp_int(&mpi_comp_mp_int);
    
    // recv data
    MPI_Recv(c_int, num, mpi_comp_mp_int, source, COMP_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    // setup & recv str
    for (i = 0; i < num; i++)
        total += c_int[i].totalLength;
    str = (char *)bmalloc(total * sizeof(char));
    MPI_Recv(str, total, MPI_CHAR, source, COMP_MP, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup c
    for (i = 0; i < num; i++)
    { // setup c[i]
        tempStr = &str[currLoc];
        cp_comp_mp_int(&c[i], &c_int[i], &tempStr, 0, 1);
        // update currLoc
        currLoc += c_int[i].totalLength;
    }
    
    // clear data
    free(str);
    tempStr = NULL;
    
    
    // free data
    free(c_int);
    MPI_Type_free(&mpi_comp_mp_int);
    
    return;
}







void send_comp_rat(const mpq_t c[2], int target)
{
    comp_rat_int c_int;
    char *str = NULL;
    MPI_Datatype mpi_comp_rat_int;
    create_comp_rat_int(&mpi_comp_rat_int);
    
    // send data
    cp_comp_rat_int(&c_int, const_cast<mpq_t*>(c), &str, 0, 0);
    // send c_int
    MPI_Send(&c_int, 1, mpi_comp_rat_int, target, COMP_RAT, MPI_COMM_WORLD);
    // send str
    MPI_Send(str, c_int.length[0] + c_int.length[1], MPI_CHAR, target, COMP_RAT, MPI_COMM_WORLD);
    
    // clear str
    free(str);
    
    
    MPI_Type_free(&mpi_comp_rat_int);
    
    return;}
void receive_comp_rat(mpq_t c[2], int source)
{
    MPI_Status statty_mc_gatty;
    comp_rat_int c_int;
    char *str = NULL;
    MPI_Datatype mpi_comp_rat_int;
    create_comp_rat_int(&mpi_comp_rat_int);
    
    // recv data
    MPI_Recv(&c_int, 1, mpi_comp_rat_int, source, COMP_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    // setup & recv str
    str = (char *)bmalloc((c_int.length[0] + c_int.length[1]) * sizeof(char));
    MPI_Recv(str, c_int.length[0] + c_int.length[1], MPI_CHAR, source, COMP_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup c
    cp_comp_rat_int(c, &c_int, &str, 1, 1);
    
    
    MPI_Type_free(&mpi_comp_rat_int);
    
    return;
}









void send_comp_num_rat(const mpq_t c[][2], int num, int target)
{
    int i, j, total = 0, currLoc = 0;
    comp_rat_int *c_int = (comp_rat_int *)bmalloc(num * sizeof(comp_rat_int));
    char *str = NULL, *tempStr = NULL;
    MPI_Datatype mpi_comp_rat_int;
    create_comp_rat_int(&mpi_comp_rat_int);
    
    // send data
    for (i = 0; i < num; i++)
    { // setup c_int[i]
        cp_comp_rat_int(&c_int[i], const_cast<mpq_t*>(c[i]), &tempStr, 0, 0);
        // update
        total += c_int[i].length[0] + c_int[i].length[1];
        str = (char *)brealloc(str, total * sizeof(char));
        for (j = 0; j < c_int[i].length[0] + c_int[i].length[1]; j++)
        {
            str[currLoc] = tempStr[j];
            currLoc++;
        }
        free(tempStr);
    }
    // send c_int
    MPI_Send(c_int, num, mpi_comp_rat_int, target, COMP_RAT, MPI_COMM_WORLD);
    // send str
    MPI_Send(str, total, MPI_CHAR, target, COMP_RAT, MPI_COMM_WORLD);
    
    // clear str
    free(str);
    
    MPI_Type_free(&mpi_comp_rat_int);
    
    return;
}
void receive_comp_num_rat(mpq_t c[][2], int num, int source)
{
    MPI_Status statty_mc_gatty;
    int i, total = 0, currLoc = 0;
    comp_rat_int *c_int = (comp_rat_int *)bmalloc(num * sizeof(comp_rat_int));
    char *str = NULL, *tempStr = NULL;
    MPI_Datatype mpi_comp_rat_int;
    create_comp_rat_int(&mpi_comp_rat_int);
    
    // recv data
    MPI_Recv(c_int, num, mpi_comp_rat_int, source, COMP_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    // setup & recv str
    for (i = 0; i < num; i++)
        total += c_int[i].length[0] + c_int[i].length[1];
    str = (char *)bmalloc(total * sizeof(char));
    MPI_Recv(str, total, MPI_CHAR, source, COMP_RAT, MPI_COMM_WORLD, &statty_mc_gatty);
    
    // setup c
    for (i = 0; i < num; i++)
    { // setup c[i]
        tempStr = &str[currLoc];
        cp_comp_rat_int(c[i], &c_int[i], &tempStr, 0, 1);
        // update currLoc
        currLoc += c_int[i].length[0] + c_int[i].length[1];
    }
    
    
    MPI_Type_free(&mpi_comp_rat_int);
    
    return;
}





int compare_integers_decreasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left<right) {
		return 1;
	}
	else if (left>right){
		return -1;
	}
	else{
		return 0;
	}
	
}

int compare_integers_increasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left>right) {
		return 1;
	}
	else if(left < right){
		return -1;
	}
	else{
		return 0;
	}
	
}




