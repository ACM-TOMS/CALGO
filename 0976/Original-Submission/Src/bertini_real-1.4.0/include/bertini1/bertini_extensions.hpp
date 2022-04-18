#ifndef BERTINI_EXTENSIONS_H
#define BERTINI_EXTENSIONS_H

/** 
 \file bertini_extensions.hpp 
 
 \brief Direct extensions to Bertini, including single-target/source MPI sending and receiving
 */

#include "bertini1/bertini_headers.hpp"
#include "io/fileops.hpp"

#define TOL_DOUBLE_PRECISION 1e-13
#define LARGECHANGE_DOUBLEPRECISION 1e14


#define TOL_MP 1e-40
#define LARGECHANGE_MP 1e50


enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};


enum {SYSTEM_CRIT = -1600, SYSTEM_SPHERE};

using VertexType = int;

constexpr VertexType Unset = 0;
constexpr VertexType Critical = 1;
constexpr VertexType Semicritical = 2;
constexpr VertexType Midpoint = 4;
constexpr VertexType Isolated = 8;
constexpr VertexType New = 16;
constexpr VertexType Curve_sample_point = 32;
constexpr VertexType Surface_sample_point = 64;
constexpr VertexType Removed = 128;
constexpr VertexType Problematic = 256;

// //The following lets us use words instead of numbers to indicate vertex type.
// enum VertexType {Unset= 100, Critical, Semicritical, Midpoint, Isolated, New, Curve_sample_point, Surface_sample_point, Removed, Problematic};

constexpr VertexType VertexTypes[]{Unset, Critical, Semicritical, Midpoint, Isolated, New, Curve_sample_point, Surface_sample_point, Removed, Problematic};



// enum for worker mode choice
enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN, MIDPOINT_SOLVER, SPHERE_SOLVER, BERTINI_MAIN};

enum {TERMINATE = 2000, INITIAL_STATE};



enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, VEC_RAT, MAT_RAT, COMP_RAT, INDICES, DECOMPOSITION, CURVE, SURFACE, EDGE, CELL, FACE, UNUSED, VERTEX_SET, WITNESS_SET, VERTEX, SYSTEM_RANDOMIZER};



enum {INACTIVE = 500, ACTIVE};
enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};



extern "C"{
	void parallel_diff_worker(int my_id, int num_processes, int headnode);
}



/**
 look up an integer, for what it means in the world of enumerations.
 */
std::string enum_lookup(int flag, int hint=0);





int bertini_main_wrapper(const std::vector<std::string> & options,  int num_processes, int my_id, int headnode);



/**
 Print the Bertini splash screen to the screen.
 */
void bertini_splash_screen();



/**
 \brief bertini_real's malloc.
 
 will call br_exit if it fails to malloc correctly.
 
 \return a pointer to the memory allocated.
 \param size the amount of memory to allocate.
 */
void *br_malloc(size_t size);

/**
 \brief bertini_real's realloc.
 
 will call br_exit if it fails to realloc correctly.
 
 \return a pointer to memory reallocated
 \param ptr the pointer to be reallocated
 \param size the new size of the memory to be reallocated.
 */
void *br_realloc(void *ptr, size_t size);


/**
 \brief test whether a matrix is the identity
 
 \return boolean indicating whether the input matrix is the identity.
 \param M test matrix
 */
bool is_identity(const mat_d M);

/**
 \brief test whether a matrix is the identity
 
 \return boolean indicating whether the input matrix is the identity.
 \param M test matrix
 */
bool is_identity(const mat_mp M);






//function prototypes for bertini_real data clearing etc.


/**
 \brief compute the (multiple-precision) 2-norm of the difference of two vectors of the same length
 
 \param result the returned computed norm
 \param left First input vector
 \param right Second input vector
 */
void norm_of_difference(mpf_t result, const vec_mp left, const vec_mp right);


/**
 \brief compute the (multiple-precision) 2-norm of the difference of two vectors of the differing length, using the minimal dimension.
 
 \param result the returned computed norm
 \param left First input vector
 \param right Second input vector
 */
void norm_of_difference_mindim(mpf_t result, const vec_mp left, const vec_mp right);


/**
 \brief Dehomogenize a vector assuming there is a single variable group and the leading coordinate is the homogenizing variable.
 
 \param result The computed value.
 \param dehom_me The vector to dehomogenize
 */
void dehomogenize(point_d *result, const point_d dehom_me);
/**
 \brief Dehomogenize only the first n variables of a vector assuming there is a single variable group and the leading coordinate is the homogenizing variable.
 
 \param result The computed value.
 \param dehom_me The vector to dehomogenize
 \param num_variables the total number of variables, including the homogenizing variable.
 */
void dehomogenize(point_d *result, const point_d dehom_me, int num_variables);



/**
 \brief Dehomogenize a vector assuming there is a single variable group and the leading coordinate is the homogenizing variable.
 
 \param result The computed value.
 \param dehom_me The vector to dehomogenize
 */
void dehomogenize(point_mp *result, const point_mp dehom_me);
/**
 \brief Dehomogenize only the first n variables of a vector assuming there is a single variable group and the leading coordinate is the homogenizing variable.
 
 \param result The computed value.
 \param dehom_me The vector to dehomogenize
 \param num_variables the total number of variables, including the homogenizing variable.
 */
void dehomogenize(point_mp *result, const point_mp dehom_me, int num_variables);



/**
 \brief Compute the non-conjugate transpose of a complex matrix M.
 
 \param Res the result of the operation
 \param M the input matrix
 */
void nonconj_transpose(mat_d Res, const mat_d M);
/**
 \brief Compute the non-conjugate transpose of a complex matrix M.
 
 \param Res the result of the operation
 \param M the input matrix
 */
void nonconj_transpose(mat_mp Res, const mat_mp M);



/**
 \brief Compute the non-conjugate dot product of two vectors of the same size only
 
 \param result The computed value.
 \param one The first input vector.
 \param two The second input vector.
 */
void dot_product_d(comp_d result, const vec_d one, const vec_d two);
/**
 \brief Compute the non-conjugate dot product of two vectors of the same size only
 
 \param result The computed value.
 \param one The first input vector.
 \param two The second input vector.
 */
void dot_product_mp(comp_mp result, const vec_mp one, const vec_mp two);


/**
 \brief Compute the non-conjugate dot product of two complex vectors of possibly differing length, only for the first computable number of variables.
 
 \param result The computed result
 \param left The first input vector.
 \param right The second input vector.
 */
void dot_product_mindim(comp_d result, const vec_d left, const vec_d right);
/**
 \brief Compute the non-conjugate dot product of two complex vectors of possibly differing length, only for the first computable number of variables.
 
 \param result The computed result
 \param left The first input vector.
 \param right The second input vector.
 */
void dot_product_mindim(comp_mp result, const vec_mp left, const vec_mp right);



/**
 \brief Compute the determinant of a matrix
 
 \return The integer 0.  Seems dumb.
 \param determinant The computed value.
 \param source_matrix The matrix to compute the determinant of.
 */
int take_determinant_d(comp_d determinant, const mat_d source_matrix);
/**
 \brief Compute the determinant of a matrix
 
 \return The integer 0.  Seems dumb.
 \param determinant The computed value.
 \param source_matrix The matrix to compute the determinant of.
 */
int take_determinant_mp(comp_mp determinant, const mat_mp source_matrix);



/**
 \brief Compute the linear projection value of a point, given homogenenous coordinates for the point.
 
 Dehomogenizes the point before inner-producting it with the projection vector.
 \param result The computed value.
 \param input The input point.
 \param projection The vector representing the linear projection.
 */
void projection_value_homogeneous_input(comp_d result, const vec_d input, const vec_d projection);
/**
 \brief Compute the linear projection value of a point, given homogenenous coordinates for the point.
 
 Dehomogenizes the point before inner-producting it with the projection vector.
 \param result The computed value.
 \param input The input point.
 \param projection The vector representing the linear projection.
 */
void projection_value_homogeneous_input(comp_mp result, const vec_mp input, const vec_mp projection);


/**
 \brief tests whether two point given in already-dehomogenized form are the same, using a defined threshold.
 
 no need to dehomogenize - is inhomogeneous by assumption
 
 \return boolean integer indicating whether ||left-right||_2<tol
 \param left The first input
 \param right The second input
 \param tolerance The L_2 separation distance for the two points to be considered the same
 */
int isSamePoint_inhomogeneous_input(const point_d left, const point_d right, double tolerance);
/**
 \brief tests whether two point given in already-dehomogenized form are the same, using a defined threshold.
 
  no need to dehomogenize - is inhomogeneous by assumption
 
 \return boolean integer indicating whether ||left-right||_2<tol
 \param left The first input
 \param right The second input
 \param tolerance The L_2 separation distance for the two points to be considered the same
 */
int isSamePoint_inhomogeneous_input(const point_mp left, const point_mp right, double tolerance);



/**
 \brief tests whether two point given in homogenized form are the same, using a defined threshold.
 
 dehomogenizes the points before testing.
 
 \return boolean integer indicating whether ||left-right||_2<tol
 \param left The first input
 \param right The second input
 \param tolerance The L_2 separation tolerance for the points
 */
int isSamePoint_homogeneous_input(const point_d left, const point_d right, double tolerance);
/**
 \brief tests whether two point given in homogenized form are the same, using a defined threshold.
 
  dehomogenizes the points before testing.
 
 \return boolean integer indicating whether ||left-right||_2<tol
 \param left The first input
 \param right The second input
 \param tolerance The L_2 separation tolerance for the points
 */
int isSamePoint_homogeneous_input(const point_mp left, const point_mp right, double tolerance);


/**
 \brief thresholds a number \f$x\f$ so that if \f$|Im{x}|<\eps \f$, we set \f$Im(x) = 0\f$.
 
 \param blabla the input AND output value.  it changes the input directly.
 \param threshold the value of epsilon, so that if the imaginary part is large (in abolute value), we set it to 0.
 */
void real_threshold(comp_mp blabla, double threshold);
/**
 \brief thresholds a number \f$x\f$ so that if \f$|Im{x}|<\eps \f$, we set \f$Im(x) = 0\f$.
 
 \param blabla the input AND output value.  it changes the input directly.
 \param threshold the value of epsilon, so that if the imaginary part is large (in abolute value), we set it to 0.
 */
void real_threshold(vec_mp blabla, double threshold);
/**
 \brief thresholds a number \f$x\f$ so that if \f$|Im{x}|<\eps \f$, we set \f$Im(x) = 0\f$.
 
 \param blabla the input AND output value.  it changes the input directly.
 \param threshold the value of epsilon, so that if the imaginary part is large (in abolute value), we set it to 0.
 */
void real_threshold(mat_mp blabla, double threshold);



/**
 \brief prints a message to screen based on a Bertini retVal.
 
 \param retVal The returned integer value from a Bertini tracker loop. 0 means total success, is something else else.
 */
void print_path_retVal_message(int retVal);

/**
 retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
 */
/**
 \brief Computes the total number of variables in a Bertini setup, from a preproc_data.
 
 takes the sum of the numbers in all variable groups including hom and non-hom, plus the sum of the number of variable_groups (the number of homogenizing coordinates).
 
 \return The computed number of variables
 \param PPD the preproc_data from which to compute.  This must be populated from the parsed input file elsewhere.
 */
int get_num_vars_PPD(const preproc_data PPD);


/**
 \brief copy a bertini patch structure
 
 \param PED The patch into which to copy.
 \param PED_input The patch from which to copy.
 */
void cp_patch_mp(patch_eval_data_mp *PED, const patch_eval_data_mp PED_input);
/**
 \brief copy a bertini patch structure
 
 \param PED The patch into which to copy.
 \param PED_input The patch from which to copy.
 */
void cp_patch_d(patch_eval_data_d *PED, const patch_eval_data_d PED_input);



/**
 \brief Copy a Bertini preproc_data structure
 
 \param PPD The Output result preproc_data
 \param PPD_input the input preproc_data, from which we copy.
 */
void cp_preproc_data(preproc_data *PPD, const preproc_data & PPD_input);


/**
 \brief Clear a post_process_t struct from Bertini.
 
 \param endPoint The struct to clear.
 \param num_vars The number of variables appearing in it.  Don't worry, I hate it, too.
 */
void clear_post_process_t(post_process_t * endPoint, int num_vars);


/**
 \brief Prints all fields from a tracker_config_t structure to the screen.
 
 \param T The pointer to the tracker_config_t.
 */
void print_tracker(const tracker_config_t * T);


/**
 \brief Sort a vector of real-valued comp_mp's by their real values.
 
 If two values are closer than some smaller than 1e-10, they are called the same thing...
 \todo remove the hard-coded tolerance for declaring two values the same, and make it user-controllable.
 
 \return The integer number 0.  Seems dumb.
 \param projections_sorted The output value.
 \param index_tracker The order the inputs get sorted into.  It's a permutation vector.
 \param projections_input The projection values you want to sort.
 \param distinct_thresh The separation value for determining if two values are distinct or not.  You should be able to set this pretty small, but it should also be a function of the level of sharpening /  accuracy requested from the tracker.
 */
int sort_increasing_by_real(vec_mp projections_sorted, std::vector< int > & index_tracker, const vec_mp projections_input, double distinct_thresh);


/**
 \brief A comparitor for integers, to sort them using qsort into decreasing order
 
 \return integer indicating which is larger.  if left < right, return 1.  if right < left return -1. return 0 if equal.
 \param left_in pointer to input
 \param right_in pointer to input
 */
int compare_integers_decreasing(const void * left_in, const void * right_in);
/**
 \brief A comparitor for integers, to sort them using qsort into increasing order
 
 \return integer indicating which is larger.  if left > right, return 1.  if right > left return -1. return 0 if equal.
 \param left_in pointer to input
 \param right_in pointer to input
 */
int compare_integers_increasing(const void * left_in, const void * right_in);




/**
 \brief Broadcast send a patch to everyone
 
 \param patch a pointer to a patch to send
 */
void send_patch_mp   (const patch_eval_data_mp * patch);
/**
 \brief Broadcast receive a patch from 0
 
 \param patch a pointer to a patch to receive into
 */
void receive_patch_mp(patch_eval_data_mp * patch);



/**
 \brief Broadcast send a patch to everyone
 
 \param patch a pointer to a patch to send
 */
void send_patch_d   (const patch_eval_data_d * patch);
/**
 \brief Broadcast receive a patch from 0
 
 \param patch a pointer to a patch to receive into
 */
void receive_patch_d(patch_eval_data_d * patch);




/**
 \brief Broadcast send a preproc_data to everyone
 
 \param PPD a pointer to a preproc_data to send
 */
void send_preproc_data(const preproc_data *PPD);
/**
 \brief Broadcast receive a preproc_data from 0
 
 \param PPD a pointer to a preproc_data to receive into
 */
void receive_preproc_data(preproc_data *PPD);




/**
 \brief send a matrix to a single target.
 
 \param A matrix to send
 \param target Where to send it relative to MPI_COMM_WORLD
 */
void send_mat_d(const mat_d A, int target);
/**
 \brief receive a matrix from a single source.
 
 \param A matrix to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_mat_d(mat_d A, int source);


/**
 \brief send a matrix to a single target.
 
 \param A matrix to send
 \param target Where to send it relative to MPI_COMM_WORLD
 */
void send_mat_mp(const mat_mp A, int target);
/**
 \brief receive a matrix from a single source.
 
 \param A matrix to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_mat_mp(mat_mp A, int source);


/**
 \brief send a matrix to a single target.
 
 \param A_d The double matrix to send
 \param A_mp The mp matrix to send
 \param A_rat The rational matrix to send
 \param target Where to send it relative to MPI_COMM_WORLD
 */
void send_mat_rat(const mat_d A_d, const mat_mp A_mp, const mpq_t ***A_rat, int target);
/**
 \brief Simultaneously receive a double, mp, and rational matrix from a single source.
 
 \param A_d The double matrix to receive into
 \param A_mp The mp matrix to receive into
 \param A_rat The rational matrix to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_mat_rat(mat_d A_d, mat_mp A_mp, mpq_t ***A_rat, int source);



/**
 \brief send a vector to a single target.
 
 \param b matrix to send
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_vec_d(const vec_d b, int target);
/**
 \brief receive a vector from a single source.
 
 \param b matrix to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_vec_d(vec_d b, int source);

/**
 \brief send a vector to a single target.
 
 \param b matrix to send
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_vec_mp(const vec_mp b, int target);
/**
 \brief receive a vector from a single source.
 
 \param b matrix to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_vec_mp(vec_mp b, int source);


/**
 \brief send a vector to a single target.
 
 \param b matrix to send
 \param size the number of entries
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_vec_rat(const mpq_t ***b, int size, int target);
/**
 \brief receive a vector from a single source.
 
 \param b matrix to receive into
 \param size the number of entries
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_vec_rat(mpq_t ***b, int size, int source);





/**
 \brief send a complex number to a single target.
 
 \param c Number to send
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_comp_d(const comp_d c, int target);
/**
 \brief Receive a complex number from a single source.
 
 \param c Number to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_comp_d(comp_d c, int source);










/**
 \brief send a complex number to a single target.
 
 \param c Number to send
 \param num How many there are
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_comp_num_d(const comp_d *c, int num, int target);
/**
 \brief Receive an array of complex numbers from a single source.
 
 \param c Number to receive into
 \param num How many there are
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_comp_num_d(comp_d *c, int num, int source);








/**
 \brief send a complex number to a single target.
 
 \param c Number to send
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_comp_mp(const comp_mp c, int target);
/**
 \brief Receive a complex number from a single source.
 
 \param c Number to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_comp_mp(comp_mp c, int source);


/**
 \brief send an array of complex numbers to a single target.
 
 \param c array of numbers to send
 \param num How many of them there are.
 \param target Where to send them, relative to MPI_COMM_WORLD
 */
void send_comp_num_mp(const comp_mp *c, int num, int target);
/**
 \brief Receive an array of complex numbers from a single target.
 
 \param c array of numbers to receive into
 \param num How many of them there are.
 \param source Where to receive them from, relative to MPI_COMM_WORLD
 */
void receive_comp_num_mp(comp_mp *c, int num, int source);


/**
 \brief send a single complex number to a single target.
 
 \param c number to send
 \param target Where to send it, relative to MPI_COMM_WORLD
 */
void send_comp_rat(const mpq_t c[2], int target);
/**
 \brief Receive a complex number from a single source.
 
 \param c Number to receive into
 \param source Where to receive it from, relative to MPI_COMM_WORLD
 */
void receive_comp_rat(mpq_t c[2], int source);


/**
 \brief send an array of complex numbers to a single target.
 
 \param c array of numbers to send
 \param num How many of them there are.
 \param target Where to send them, relative to MPI_COMM_WORLD
 */
void send_comp_num_rat(const mpq_t c[][2], int num, int target);
/**
 \brief Receive an array of complex numbers from a single target.
 
 \param c array of numbers to receive into
 \param num How many of them there are.
 \param source Where to receive them from, relative to MPI_COMM_WORLD
 */
void receive_comp_num_rat(const mpq_t c[][2], int num, int source);





/**
 \brief Print a vector to the screen for copypasta into Matlab's command window or a file.
 
 \param M the vector to print
 \param name the string of the name to give it.
 */
void print_point_to_screen_matlab(const vec_d M, std::string name);
/**
 \brief Print a vector to the screen for copypasta into Matlab's command window or a file.
 
 \param M the vector to print
 \param name the string of the name to give it.
 */
void print_point_to_screen_matlab(const vec_mp M, std::string name);
/**
 \brief Print a matrix to the screen for copypasta into Matlab's command window or a file.
 
 \param M the matrix to print
 \param name the string of the name to give it.
 */
void print_matrix_to_screen_matlab(const mat_d M, std::string name);
/**
 \brief Print a matrix to the screen for copypasta into Matlab's command window or a file.
 
 \param M the matrix to print
 \param name the string of the name to give it.
 */
void print_matrix_to_screen_matlab(const mat_mp M, std::string name);



/**
 \brief Print a complex number to the screen for copypasta into Matlab's command window or a file.
 
 \param M the number to print
 \param name the string of the name to give it.
 */
void print_comp_matlab(const comp_mp M,std::string name);
/**
 \brief Print a complex number to the screen for copypasta into Matlab's command window or a file.
 
 \param M the number to print
 \param name the string of the name to give it.
 */
void print_comp_matlab(const comp_d M,std::string name);





#endif


