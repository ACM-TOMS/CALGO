#pragma once

/** \file data_type.hpp */

#include "bertini1/bertini_extensions.hpp"
#include "programConfiguration.hpp"
/**
 \brief Comprehensive system randomization, based on deg.out.
 
 The SystemRandomizer is created based on the desired size to randomize down to, and the degrees of the functions, which are contained in the 'deg.out' file, which must pre-exist.
 
 This class does not keep track of the desired mp mode, and always populates all three randomizer matrices.
 
 This class is capable of randomizing for systems with a single hom_variable_group (set hom_var = 1), and a single variable_group.
 
 */
class SystemRandomizer
{
private:
	bool square_indicator; ///< a boolean indicating whether the system is square.
	bool setup_indicator; ///< a boolean indicating whether the randomizer is ready to use.
	
	mat_mp randomizer_matrix_full_prec; ///< holds the full precision randomizer matrix, from which we downsample when changing precision in MP mode.
	mat_mp randomizer_matrix_mp; ///< holds the randomizer matrix to current precision.
	mat_d randomizer_matrix_d; ///< holds the randomization matrix in double precision.
	
	int num_randomized_funcs; ///< the number of functions to which we randomize.
	int num_original_funcs; ///< the number of function from which we randomize.
	
	int max_base_degree; ///< the highest degree of any function occurring in the system.
	int max_degree_deficiency; ///< the greatest occuring deficiency in function degree relative to the max.
	std::vector<int> randomized_degrees; ///< a vector of integers keeping track of the degrees of the output functions.
	std::vector<int> original_degrees; ///< a vector of integers keeping track of the degrees of the input functions.
	
	std::vector<std::vector<int> > structure_matrix; ///< matrix of integers indicating the degree deficiency of each input function relative to the degree of the output functions.
	
	
	vec_mp integer_coeffs_mp;
	mat_mp single_row_input_mp;
	vec_mp temp_homogenizer_mp;
	vec_mp temp_funcs_mp;
	mat_mp temp_jac_mp;
	mat_mp temp_mat_mp;
	vec_mp temp_vec_mp;
	
	
	vec_d integer_coeffs_d;
	mat_d single_row_input_d;
	vec_d temp_homogenizer_d;
	vec_d temp_funcs_d;
	mat_d temp_jac_d;
	mat_d temp_mat_d;
	vec_d temp_vec_d;
	
	
public:
	
	
	/**
	 not all functions are of the same degree, and when randomizing, if there is a deficiency of a function with respect to another, you must homogenize.  this helps give that info.
	 
	 \param randomized_func The index of the output function.
	 \param base_index The index of the input function.
	 \return the level of deficiency of randomized_func with respect to base_index.
	 */
	int deficiency(unsigned int randomized_func, unsigned int base_index) const
	{
		return structure_matrix[randomized_func][base_index];
	}
	
	
	/**
	 get a reference to the degrees of the output functions for this randomizer.
	 
	 \return the degrees in vector form.
	 */
	const std::vector<int>& rand_degrees() const
	{
		return randomized_degrees;
	}
	
	
	/**
	 get the degree of a particular randomized function.
	 
	 \throws out_of_range if there isn't any such function.
	 \param index the index of the output function.
	 \return the degree of the output function
	 */
	int randomized_degree(unsigned int index) const
	{
		if (index>= randomized_degrees.size()) {
			throw std::out_of_range("trying to access a randomized degree out of range");
		}
		
		return randomized_degrees[index];
	}
	
	

	
	/**
	 output to a stream.  only really usable with std::cout, as it calls print_matrix_to_screen_matlab
	 
	 \param os the stream to put this text on.
	 \param s the system randomizer to write.
	 */
	friend std::ostream & operator<<(std::ostream &os, SystemRandomizer & s)
	{
		os << "square: " << s.square_indicator << ", is_ready: " << s.setup_indicator << std::endl;
		os << "num_rand: " <<  s.num_randomized_funcs << ", num_orig " << s.num_original_funcs << std::endl;
		os << "max_base_deg: " << s.max_base_degree << ", max_deficiency " << s.max_degree_deficiency << std::endl;
		
		print_matrix_to_screen_matlab(s.randomizer_matrix_d,"rand_mat");
		
		for (auto iter = s.structure_matrix.begin(); iter!=s.structure_matrix.end(); ++iter) {
			for (auto jter = iter->begin(); jter!=iter->end(); ++jter) {
				os << *jter << " ";
			}
			os << std::endl;
		}
		os << std::endl << std::endl;
		
		print_point_to_screen_matlab(s.integer_coeffs_d,"int_coeffs");
		
		std::cout << "randomized" << std::endl;
		for (auto iter = s.randomized_degrees.begin(); iter!=s.randomized_degrees.end(); ++iter) {
			std::cout << *iter << " ";
		}
		os << std::endl << std::endl;
		
		std::cout << "original" << std::endl;
		for (auto iter = s.original_degrees.begin(); iter!=s.original_degrees.end(); ++iter) {
			std::cout << *iter << " ";
		}
		os << std::endl;
		
		//TODO: add mode output here
		return os;
	}
	
	
	/**
	 constructor
	 */
	SystemRandomizer()
	{
		init();
	}
	
	
	/**
	 assignment operator
	 \param other another SystemRandomizer to copy from.
	 */
	SystemRandomizer & operator=( const SystemRandomizer & other)
	{
		copy(other);
		return *this;
	}
	
	/**
	 copy constructor
	 */
	SystemRandomizer(const SystemRandomizer & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	/**
	 destructor
	 
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
	~SystemRandomizer()
	{
		clear_mat_mp(randomizer_matrix_full_prec);
		clear_mat_mp(randomizer_matrix_mp);
		clear_mat_d(randomizer_matrix_d);
		
		clear_vec_mp(temp_homogenizer_mp);
		clear_vec_d(temp_homogenizer_d);
		
		
		
		
		clear_vec_mp(integer_coeffs_mp);
		clear_vec_d(integer_coeffs_d);
		
		
		
		clear_mat_d(single_row_input_d);
		clear_mat_mp(single_row_input_mp);
		
		
		
		clear_vec_d(temp_funcs_d);
		clear_vec_mp(temp_funcs_mp);
		
		clear_mat_d(temp_jac_d);
		clear_mat_mp(temp_jac_mp);
		
		
		clear_mat_d(temp_mat_d);
		clear_mat_mp(temp_mat_mp);
		
		clear_vec_d(temp_vec_d);
		clear_vec_mp(temp_vec_mp);
		
		
	}
	
	/** 
	 \brief returns the number of original functions
	 \return int num_original_funcs
	 */
	int num_base_funcs() const
	{
		return num_original_funcs;
	}
	
	
	/**
	 \brief returns the number of randomized functions
	 \return int num_randomized_funcs
	 */
	int num_rand_funcs() const
	{
		return num_randomized_funcs;
	}
	
	
	/**
	 \brief returns the highest degree in the system
	 \return int max_base_degree
	 */
	int max_degree() const
	{
		return max_base_degree;
	}
	
	/**
	 \brief returns the greatest occurring deficiency
	 \return int max_degree_deficiency
	 */
	int max_deficiency() const
	{
		return max_degree_deficiency;
	}
	
	
	/**
	 \brief indicates whether the randomizer is square
	 \return bool square_indicator
	 */
	bool is_square() const
	{
		return square_indicator;
	}
	
	
	/**
	 \brief indicates whether the SystemRandomizer is ready to go.
	 \return bool setup_indicator
	 */
	bool is_ready() const
	{
		return setup_indicator;
	}
	
	/**
	 \brief gets the degree of function with index (loc).
	 \return int degree of original function with input index.
	 \param loc the index of the base function to get degree of.
	 */
	int base_degree(unsigned int loc) const
	{
		
		
		
		if (loc>=original_degrees.size()) {
			throw std::out_of_range("trying to access an original degree out of range");
		}
		
		return original_degrees[loc];
	}
	
	/**
	 \brief gets a pointer to the full precision randomizer matrix.
	 \return a pointer to the full precision randomizer matrix.
	 */
	mat_mp * get_mat_full_prec()
	{
		return &randomizer_matrix_full_prec;
	}
	
	/**
	 \brief gets a pointer to the double randomizer matrix.
	 \return a pointer to the double randomizer matrix.
	 */
	mat_d * get_mat_d()
	{
		return &randomizer_matrix_d;
	}
	
	
	/**
	 \brief gets a pointer to the current precision MP randomizer matrix.
	 \return a pointer to the current precision MP randomizer matrix.
	 */
	mat_mp * get_mat_mp()
	{
		return &randomizer_matrix_mp;
	}
	
	
	/**
	 \brief changes the precision of the SystemRandomizer
	 \param new_prec new precision.
	 */
	void change_prec(int new_prec);
	
	
	
	/**
	 \brief randomizes, by taking the input function values and jacobian, and multiplying the randomizer matrix.
	 \param randomized_func_vals output argument, set to R*f.
	 \param randomized_jacobian output jacobian, set to R*J.
	 \param func_vals input function values, probably produced by evaluating an SLP.
	 \param jacobian_vals input jacobian, probably produced by evaluating an SLP.
	 \param hom_var the coordinate of the single homogenizing variable.
	 */
	void randomize(vec_d randomized_func_vals, mat_d randomized_jacobian,
				   vec_d func_vals, mat_d jacobian_vals,
				   comp_d hom_var);
	
	
	/**
	 \brief randomizes, by taking the input function values and jacobian, and multiplying the randomizer matrix.
	 \param randomized_func_vals output argument, set to R*f.
	 \param randomized_jacobian output jacobian, set to R*J.
	 \param func_vals input function values, probably produced by evaluating an SLP.
	 \param jacobian_vals input jacobian, probably produced by evaluating an SLP.
	 \param hom_var the coordinate of the single homogenizing variable.
	 */
	void randomize(vec_mp randomized_func_vals, mat_mp randomized_jacobian,
				   vec_mp func_vals, mat_mp jacobian_vals,
				   comp_mp hom_var);
	
	
	/**
	 \brief sets up randomizer using 'deg.out'.
	 
	 This setup function parses "deg.out" for the degrees of the functions, and sets up the internals of this class object, for randomizing a system.
	 
	 \param num_desired_rows The number of output functions.  Must be bigger than the input num_funcs.
	 \param num_funcs The original number of functions in the system to be randomized.
	 */
	void setup(int num_desired_rows, int num_funcs);
	
	/** 
	 \brief sets up the temporaries.
	 set up the internal temporary variables.
	 */
	void setup_temps();
	
	
	/**
	 \brief single-target send
	 
	 send SystemRandomizer to a single target.
	 
	 \param target the id of the target.  
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void send(int target, ParallelismConfig & mpi_config);
	
	
	
	/**
	 \brief single-source receive
	 
	 receive SystemRandomizer from a single source.
	 
	 \param source the id of the source.
	 \param mpi_config container holding the mpi_config for the caller.
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	/**
	 \brief collective MPI_COMM_WORLD broadcast send
	 
	 \see VertexSet::bcast_receive
	 
	 Send the SystemRandomizer to everyone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief Collective MPI_COMM_WORLD broadcast receive
	 
	 Receive the SystemRandomizer from someone in MPI_COMM_WORLD.
	 
	 \param mpi_config The current state of MPI, as represented in Bertini_real
	 */
	void bcast_receive(ParallelismConfig & mpi_config);
	
	
protected:
	
	
	/**
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
	void init()
	{
		num_randomized_funcs = -1223;
		num_original_funcs = -976;
		
		max_degree_deficiency = -1232;
		max_base_degree = -1321;
		
		setup_indicator = false;
		
		init_mat_mp2(randomizer_matrix_full_prec,0,0,1024);
		init_mat_mp(randomizer_matrix_mp,0,0);
		init_mat_d(randomizer_matrix_d,0,0);
		
		
		init_vec_mp(integer_coeffs_mp,0);
		init_vec_d(integer_coeffs_d,0);
		
		
		
		init_mat_d(single_row_input_d,0,0);
		init_mat_mp(single_row_input_mp,0,0);
		
		
		
		init_vec_d(temp_homogenizer_d,0);
		init_vec_mp(temp_homogenizer_mp,0);
		
		init_vec_d(temp_funcs_d,0);
		init_vec_mp(temp_funcs_mp,0);
		
		init_mat_d(temp_jac_d,0,0);
		init_mat_mp(temp_jac_mp,0,0);
		
		
		init_mat_d(temp_mat_d,0,0);
		init_mat_mp(temp_mat_mp,0,0);
		
		init_vec_d(temp_vec_d,0);
		init_vec_mp(temp_vec_mp,0);
		
		
	}
	
	
	
	/**
	 because the SystemRandomizer uses Bertini data types, there are lots of pointers, and this must be done explicitly.
	 */
	void copy(const SystemRandomizer & other)
	{

		
		
		if (other.is_ready()) {
			square_indicator = other.square_indicator;
			
			num_original_funcs = other.num_original_funcs;
			num_randomized_funcs = other.num_randomized_funcs;
			
			max_base_degree = other.max_base_degree;
			max_degree_deficiency = other.max_degree_deficiency;
			
			mat_cp_d(randomizer_matrix_d, other.randomizer_matrix_d);
			
			if (other.randomizer_matrix_mp->rows>0 && other. randomizer_matrix_mp->cols>0) {
				change_prec_mat_mp(randomizer_matrix_mp, mpf_get_prec(other.randomizer_matrix_mp->entry[0][0].r));
			}
			mat_cp_mp(randomizer_matrix_mp, other.randomizer_matrix_mp);
			
			mat_cp_mp(randomizer_matrix_full_prec, other.randomizer_matrix_full_prec);
			
			// these vectors are guaranteed to be nonempty by virtue of other being already setup.
			
			original_degrees.resize(0);
			for (auto iter=other.original_degrees.begin(); iter!=other.original_degrees.end(); ++iter) {
				original_degrees.push_back(*iter);
			}
			
			randomized_degrees.resize(0);
			for (auto iter=other.randomized_degrees.begin(); iter!=other.randomized_degrees.end(); ++iter) {
				randomized_degrees.push_back(*iter);
			}
			

			
			structure_matrix.resize(other.structure_matrix.size());
			for (unsigned int ii=0; ii<other.structure_matrix.size(); ++ii) {
				structure_matrix[ii].resize(0);
				for (auto jter = other.structure_matrix[ii].begin(); jter!= other.structure_matrix[ii].end(); ++jter) {
					structure_matrix[ii].push_back(*jter);
				}
			}
			
			setup_temps();
		}
		else
		{
			setup_indicator = false;
			
			randomized_degrees.resize(0);
			original_degrees.resize(0);
			structure_matrix.resize(0);
			
			num_original_funcs = num_randomized_funcs = max_base_degree = max_degree_deficiency = -1;
			
			change_size_mat_mp(randomizer_matrix_full_prec,0,0);
			change_size_mat_mp(randomizer_matrix_mp,0,0);
			change_size_mat_d(randomizer_matrix_d,0,0);
			randomizer_matrix_d->rows = randomizer_matrix_d->cols = randomizer_matrix_mp->rows = randomizer_matrix_mp->cols = randomizer_matrix_full_prec->rows = randomizer_matrix_full_prec->cols = 0;
			
			
			
		}
		
	}
	
	

};//re: randomizer class


