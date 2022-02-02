#ifndef NULLSPACE_H_
#define NULLSPACE_H_

/** 
 \file nullspace.hpp
 
 \brief Methods for computing points such that a vector is in the left or right nullspace of a jacobian matrix, coupled with the system which generated it, and several linear projections.
 */

#include "io/fileops.hpp"
#include "nag/solvers/multilintolin.hpp"
#include "nag/solvers/nullspace.hpp"
#include "double_odometer.hpp"
#include "symbolics/derivative_systems.hpp"

/**
 \brief the main function for computing critical sets.
 
 \return SUCCESSFUL
 \param solve_out The output class for all solvers.
 \param W								input WitnessSet.
 \param randomizer	randomizes the system down to the correct number of equations.
 \param pi							the set of projections to use.
 \param ambient_dim The dimension of the object containing this critical set.
 \param target_dim The dimension of the critical object
 \param target_crit_dim The dimension of the critical object.
 \param program_options				holds the configuration for the main program.  is a pointer so that it is mutable.
 \param solve_options					holds the configuration for any solvers called.  is a pointer so that it is mutable.
 \param ns_config	NullspaceConfiguration object.  this is populated in this method.  must be empty as input.
 */
int compute_crit_nullspace(SolverOutput & solve_out, // the returned value
								const WitnessSet & W,
								std::shared_ptr<SystemRandomizer> randomizer,
								vec_mp *pi,
								int ambient_dim,
								int target_dim, // this should also be the number of vectors in the *pi entry
								int target_crit_dim,
								BertiniRealConfig & program_options,
								SolverConfiguration & solve_options,
								NullspaceConfiguration *ns_config);






/**
 \brief the main function for computing critical sets.
 
 \return SUCCESSFUL
 \param solve_out The output class for all solvers.
 \param W								input WitnessSet.
 \param randomizer	randomizes the system down to the correct number of equations.
 \param pi							the set of projections to use.
 \param ambient_dim The dimension of the object containing this critical set.
 \param target_dim The dimension of the critical object
 \param target_crit_dim The dimension of the critical object.
 \param program_options				holds the configuration for the main program.  is a pointer so that it is mutable.
 \param solve_options					holds the configuration for any solvers called.  is a pointer so that it is mutable.
 \param ns_config	NullspaceConfiguration object.  this is populated in this method.  must be empty as input.
 */
int compute_crit_nullspace_left(SolverOutput & solve_out, // the returned value
						   const WitnessSet & W,
						   std::shared_ptr<SystemRandomizer> randomizer,
						   vec_mp *pi,
						   int ambient_dim,
						   int target_dim, // this should also be the number of vectors in the *pi entry
						   int target_crit_dim,
						   BertiniRealConfig & program_options,
						   SolverConfiguration & solve_options,
						   NullspaceConfiguration *ns_config);












/**
 \brief performs the setup for the NullspaceConfiguration which is used in the compute_crit_nullspace method, and is passed into the solverNullspace.
 
 \param ns_config						the data structure we are setting up.
 \param pi the projections to use.  there could be more than used.
 \param ambient_dim the dimension of the containing object.
 \param target_dim					the dimension of the object we are detecting.
 \param target_crit_codim COdimension of the critical set.
 \param max_degree  computed -- the highest degree of any derivative of the system passed in.
 \param randomizer		how the system is randomized to the correct number of equations.
 \param W										the input WitnessSet
 \param solve_options The current solver setup.
 */
void nullspace_config_setup_left(NullspaceConfiguration *ns_config,
							vec_mp *pi, // an array of projections, the number of which is the target dimensions
							int ambient_dim,
							int target_dim,
							int target_crit_codim,
							int *max_degree, // a pointer to the value
							std::shared_ptr<SystemRandomizer> randomizer,
							const WitnessSet & W,
							SolverConfiguration & solve_options);














/**
 \brief the main function for computing critical sets via a single right nullspace vector.
 
 \return SUCCESSFUL
 \param solve_out The output class for all solvers.
 \param W								input WitnessSet.
 \param randomizer	randomizes the system down to the correct number of equations.
 \param pi							the set of projections to use.
 \param ambient_dim The dimension of the object containing this critical set.
 \param program_options				holds the configuration for the main program.  is a pointer so that it is mutable.
 \param solve_options					holds the configuration for any solvers called.  is a pointer so that it is mutable.
 \param ns_config	NullspaceConfiguration object.  this is populated in this method.  must be empty as input.
 */
int compute_crit_nullspace_right(SolverOutput & solve_out, // the returned value
						   const WitnessSet & W,
						   std::shared_ptr<SystemRandomizer> randomizer,
						   vec_mp *pi,
						   int ambient_dim,
						   BertiniRealConfig & program_options,
						   SolverConfiguration & solve_options,
						   NullspaceConfiguration *ns_config);






/**
 \brief performs the setup for the NullspaceConfiguration which is used in the compute_crit_nullspace method, and is passed into the solverNullspace.
 
 \param ns_config						the data structure we are setting up.
 \param pi the projections to use.  there could be more than used.
 \param ambient_dim the dimension of the containing object.
 \param max_degree  computed -- the highest degree of any derivative of the system passed in.
 \param randomizer		how the system is randomized to the correct number of equations.
 \param W										the input WitnessSet
 \param solve_options The current solver setup.
 */
void nullspace_config_setup_right(NullspaceConfiguration *ns_config,
								  vec_mp *pi, // an array of projections, the number of which is the target dimensions
								  int ambient_dim,
								  int *max_degree, // a pointer to the value
								  std::shared_ptr<SystemRandomizer> randomizer,
								  const WitnessSet & W,
								  SolverConfiguration & solve_options);










/**
 \brief Put a few finishing details on the output from compute_crit_nullspace
 
 \param solve_out Returning object from solver.
 \param W input witness set
 \param ns_config The nullspace solver config object.
 */
void ns_concluding_modifications(SolverOutput & solve_out,
								 const WitnessSet & W,
								 NullspaceConfiguration * ns_config);








/**
 \brief Create a Bertini input file with the nullspace system.
 
 \param output_name The desired name of the output file.
 \param input_name The name of the file from which to construct the new file.
 \param program_options The current state of Bertini_real
 \param ns_config nullspace config object.
 */
void create_nullspace_system(boost::filesystem::path output_name,
							 boost::filesystem::path input_name,
							 BertiniRealConfig & program_options,
							 NullspaceConfiguration *ns_config);





/**
 \brief Create a matlab file which will take the determinant of the jacobian matrix, and write it to a text file.
 
 
 \param output_name the desired output file's name
 \param input_name bertini input file out of which to create the new file.
 \param ns_config the nullspace configuration.
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables appear
 \param numConstants the number of constants
 \param consts The names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions
 \param funcs The names of the functions
 \param lineFuncs the lines on which the functions appear
 */
void create_matlab_determinantal_system(boost::filesystem::path output_name,
										boost::filesystem::path input_name,
										NullspaceConfiguration *ns_config,
										int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs);


/**
 \brief Create a python file which will take the determinant of the jacobian matrix, and write it to a text file.
 
 
 \param output_name the desired output file's name
 \param input_name bertini input file out of which to create the new file.
 \param ns_config the nullspace configuration.
 \param numVars the number of variables
 \param vars The names of the variables
 \param lineVars the lines on which the variables appear
 \param numConstants the number of constants
 \param consts The names of the constants
 \param lineConstants the lines on which the constants appear
 \param numFuncs the number of functions
 \param funcs The names of the functions
 \param lineFuncs the lines on which the functions appear
 */
void create_python_determinantal_system( FILE *OUT, FILE *IN,
					 NullspaceConfiguration *ns_config,
					 int numVars, char **vars, int *lineVars,
					 int numConstants, char **consts, int *lineConstants,
					 int numFuncs, char **funcs, int *lineFuncs);


#endif


