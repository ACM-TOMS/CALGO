#ifndef SOLVER_NULLSPACE_H
#define SOLVER_NULLSPACE_H


/** \file solver_nullspace.hpp */


#include "bertini1/bertini_headers.hpp"



#include "nag/solvers/solver.hpp"
#include "io/fileops.hpp"
#include "programConfiguration.hpp"
#include "nag/solvers/postProcessing.hpp"


namespace nullspace_handedness{
	enum {LEFT=0, RIGHT};
}



/**
 \brief config class for left nullspace solver
 
 \todo Rewrite nullspace code to make more of it members of this class.
 */
class NullspaceConfiguration
{
	
protected:
	std::shared_ptr<SystemRandomizer> randomizer_; ///< randomizer for the underlying supplied system

	int side_;
	
public:
	
	void set_side(int new_side){ side_ = new_side;};
	int side(){ return side_;};
	
	std::shared_ptr<SystemRandomizer> randomizer()
	{
		return randomizer_;
	}
	
	
	void set_randomizer(std::shared_ptr<SystemRandomizer> randy)
	{
		randomizer_ = randy;
	}
	
	int target_dim;   ///< r			the dimension of the real set we are looking for
	int ambient_dim;  ///< k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    ///< the COdimension of the target crit space.  must be at least one (1), and at most target_dim (r).
	int num_jac_equations;    ///< the number of equations (after the randomization)
	int num_projections;   ///< the number of projections stored
	
	int num_v_vars;  ///< N, the number of variables in original problem statement (including homogenizing variables)
	int num_natural_vars;  ///< the number of natural variables, including the homogenizing variable.
	int num_synth_vars;   ///< how many synthetic variable the input system to the overall method has. these would be the result of a previous nullspace computation, e.g.
	
	
	int max_degree;    ///< the max degree of differentiated (randomized) functions
	
	
	vec_mp **starting_linears;	///< outer layer should have as many as there are randomized equations (N-k).  inside layer has number corresponding to max of randomized_degrees
				
	
	int num_additional_linears; ///<  these linears are for slicing, when finding witness points for a positive dimensional critical set.
	vec_mp *additional_linears_terminal; ///<these linears are for slicing, when finding witness points for a positive dimensional critical set.
	vec_mp *additional_linears_starting; ///<these linears are for slicing, when finding witness points for a positive dimensional critical set.
	
	
	int num_v_linears;   ///< # is  (N), cause there are N equations in the subsystem.
	vec_mp *v_linears;    ///< the actual v-linears, which are part of the start system.
	
	vec_mp v_patch; ///<  patch equation for the v-variables, to avoid the stupid all-zero degenerate solution. length of this should be N-k+ell
	
	vec_mp *target_projection; ///< # of these should be ell

	
	void clear();
	
	NullspaceConfiguration(){
		init();
	}
	
	void print();
	
	~NullspaceConfiguration()
	{
		clear();
	}
	
	
private:
	
	void init()
	{
		target_dim = ambient_dim = target_crit_codim = num_jac_equations = -1;
		num_natural_vars = num_v_vars = num_synth_vars = -1;
		max_degree = -1;
		
		starting_linears = NULL;
		
		additional_linears_starting = additional_linears_terminal = NULL;
		
		num_v_linears = -1;
		v_linears = NULL;
		
		target_projection = NULL;
	}
};






// the mp version
// this must be defined before the double version, because double has mp.

/**
 \brief the mp version of the eval data for left nullspace method.
 */
class nullspacejac_eval_data_mp : public SolverMultiplePrecision
{
	int side_;
	
public:
	
	int side(){ return side_;};
	
	
	prog_deriv_t * SLP_derivative; ///< pointer to a SLP-like class for evaluating a system with higher-order derivatives.
	
	
	int num_jac_equations; ///< how many jacobian equations there are.
	int target_dim;   ///< r			the dimension of the real set we are looking for
	int ambient_dim;  ///< k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    ///< ell.  must be at least one (1), and at most target_dim (r).
	
	int num_natural_vars;  ///< N   number of variables in original problem statement (including homogenizing variables)
	int num_v_vars;  ///< (N-k) + (k-ell+1)
	int num_synth_vars; ///< how many pre-existing synthetic variables there are.
	
	
	int max_degree;						///< the max degree of differentiated (randomized) functions

	
	int num_additional_linears; ///< the number of slicing linears, for finding witness sets of positive-dimensional critical sets.
	vec_mp *additional_linears_terminal; ///< slicing linears
	vec_mp *additional_linears_terminal_full_prec;///< slicing linears, in full precision for AMP
	
	vec_mp *additional_linears_starting; ///< starting slicing linears.
	vec_mp *additional_linears_starting_full_prec; ///< starting slicing linears in full precision for AMP.
	
	
	
	
	
	
	vec_mp **starting_linears; ///< starting linear-product linears. outer layer should have as many as there are randomized equations. inside layer has number corresponding to randomized_degrees.
	vec_mp **starting_linears_full_prec; ///< starting linear-product linears. outer layer should have as many as there are randomized equations. inside layer has number corresponding to randomized_degrees.  in full precision for AMP.
	
	
	int num_v_linears; ///< number of starting v-linears for the linear-product move.
	vec_mp *v_linears;         ///< starting v-linears for the linear-product move.  should be as many in here as there are randomized equations
	vec_mp *v_linears_full_prec;         ///< should be as many in here as there are randomized equations
	
	vec_mp v_patch; ///< the patch equation for the new v-variables.
	vec_mp v_patch_full_prec; ///< the v-patch in full precision for AMP.
	

	
	mat_mp jac_with_proj; ///< matrix holding the jacobian with the projection vectors.  essentially a stored temporary.
	mat_mp jac_with_proj_full_prec; ///< matrix holding the jacobian with the projection vectors.  essentially a stored temporary.
	
	
	
	int num_projections; ///< how many projections there are.
	vec_mp *target_projection; ///< linear projection vectors for the nullity condition on the jacobian matrix.
	vec_mp *target_projection_full_prec; ///< linear projection vectors for the nullity condition on the jacobian matrix.
	
	
	// default initializer
	nullspacejac_eval_data_mp() : SolverMultiplePrecision(){
		
		reset_counters();
		
		init();
	}
	
	nullspacejac_eval_data_mp(int mp) : SolverMultiplePrecision(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~nullspacejac_eval_data_mp(){
		
		nullspacejac_eval_data_mp::clear();
		// no need to reset the counters.
	}
	
	void print();
    
    
    
    
	void reset_counters()
	{
		num_jac_equations = 0;
		target_dim = 0;   // r			the dimension of the real set we are looking for
		ambient_dim = 0;  // k			the dimension of the complex component we are looking IN.
		target_crit_codim = 0;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		num_v_vars = -771;  // N   number of variables in original problem statement (including homogenizing variables)
		num_natural_vars = -791;  // N-k+\ell
		num_synth_vars = -1341;
		
		
		max_degree = 0;						// the max degree of differentiated (randomized) functions
		
		
		num_projections = 0;
		num_v_linears = 0;
		num_additional_linears = 0;
	} // re: reset_counters
	
	
	nullspacejac_eval_data_mp & operator=(const nullspacejac_eval_data_mp & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	nullspacejac_eval_data_mp(const nullspacejac_eval_data_mp & other)
	{
		
		
		SolverMultiplePrecision();
		nullspacejac_eval_data_mp();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	/**
	 \brief MPI broadcast send to all workers in communicator, for the nullspace solver.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI.
	 */
	int send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief MPI broadcast recive from the head of the communicator, for the nullspace solver.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI.
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	
	void set_sidedness(int which_side);
	
	void set_function_handles(int which_side);
	
	/**
	 \brief populate the eval_data form the config object.
	 
	 \return SUCCESSFUL
	 \param _SLP a pointer to the underlying SLP.
	 \param ns_config a pointer to the NullspaceConfiguration from which to set up.
	 \param W input witness set containing the patches, etc.
	 \param solve_options The current state of the solver.
	 */
	int setup(prog_t * _SLP,
			  NullspaceConfiguration *ns_config,
			  WitnessSet & W,
			  SolverConfiguration & solve_options);
	
	
	
private:
	
	
	void clear()
	{
		
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_mp(additional_linears_terminal[ii]);
			}
			free(additional_linears_terminal);
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_mp(additional_linears_starting[ii]);
			}
			free(additional_linears_starting);
		}
		
		
	
		
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
		
		
		
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				clear_vec_mp(v_linears[ii]);
			}
			free(v_linears);
		}
		
		
		clear_vec_mp(v_patch);
		
		clear_mat_mp(jac_with_proj);
		

		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(target_projection[ii]);
			}
			free(target_projection);
		}
		
		
		
		if (this->MPType==2) {
			
			if (num_additional_linears>0) {
				for (int ii=0; ii<num_additional_linears; ii++) {
					clear_vec_mp(additional_linears_terminal_full_prec[ii]);
				}
				free(additional_linears_terminal_full_prec);
				
				for (int ii=0; ii<num_additional_linears; ii++) {
					clear_vec_mp(additional_linears_starting_full_prec[ii]);
				}
				free(additional_linears_starting_full_prec);
			}
			
			
			
			
			if (num_jac_equations>0) {
				if (side_ == nullspace_handedness::LEFT) {
					for (int ii=0; ii<num_jac_equations; ii++) {
						for (int jj=0; jj<max_degree; jj++) {
							clear_vec_mp(starting_linears_full_prec[ii][jj]);
						}
						free(starting_linears_full_prec[ii]);
					}
				}
				else
				{
					for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
						for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
							clear_vec_mp(starting_linears_full_prec[ii][jj]);
						}
						free(starting_linears_full_prec[ii]);
					}
				}
				
				free(starting_linears_full_prec);
			}
			
			
			
			
			for (int ii=0; ii<num_v_linears; ii++)
				clear_vec_mp(v_linears_full_prec[ii]);
			
			free(v_linears_full_prec);
			
			clear_vec_mp(v_patch_full_prec);
			clear_mat_mp(jac_with_proj_full_prec);

			
			
			if (num_projections>0) {
				for (int ii=0; ii<num_projections; ii++) {
					clear_vec_mp(target_projection_full_prec[ii]);
				}
				free(target_projection_full_prec);
			}
		}
		
		clear_deriv(SLP_derivative);
		delete this->SLP_derivative;
		
	} // re: clear
	
	void init();
	
	void copy(const nullspacejac_eval_data_mp & other)
	{
		
		clear();
		reset_counters();
		
		
		
		
		if (other.num_additional_linears>0) {
			
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_terminal_full_prec = (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting_full_prec = (vec_mp *) br_malloc(other.num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<other.num_additional_linears; ii++) {
				vec_cp_mp(this->additional_linears_terminal[ii],						other.additional_linears_terminal[ii]);
				vec_cp_mp(this->additional_linears_terminal_full_prec[ii],	other.additional_linears_terminal_full_prec[ii]);
				
				vec_cp_mp(this->additional_linears_starting[ii],						other.additional_linears_starting[ii]);
				vec_cp_mp(this->additional_linears_starting_full_prec[ii],	other.additional_linears_starting_full_prec[ii]);
			}
		}
		else {} // other.num_additional_linears == 0
		
		this->num_additional_linears = other.num_additional_linears;
		
		if (other.num_jac_equations) {
			this->starting_linears_full_prec = (vec_mp **) br_malloc(other.num_jac_equations*sizeof(vec_mp *));
			this->starting_linears = (vec_mp **) br_malloc(other.num_jac_equations*sizeof(vec_mp *));
			
			for (int ii=0; ii<other.num_jac_equations; ii++) {
				this->starting_linears_full_prec[ii] = (vec_mp *) br_malloc(other.max_degree*sizeof(vec_mp ));
				this->starting_linears[ii] = (vec_mp *) br_malloc(other.max_degree*sizeof(vec_mp ));
				for (int jj=0; jj<other.max_degree; jj++) {
					init_vec_mp(this->starting_linears[ii][jj],0);
					init_vec_mp2(this->starting_linears_full_prec[ii][jj],0,1024);
					vec_cp_mp(this->starting_linears[ii][jj],other.starting_linears[ii][jj]);
					vec_cp_mp(this->starting_linears_full_prec[ii][jj],other.starting_linears_full_prec[ii][jj]);
				}
			}
		}
		else{}
		
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->max_degree = other.max_degree;
		
		
		if (other.num_v_linears) {
			this->v_linears = (vec_mp *) br_malloc(other.num_v_linears*sizeof(vec_mp));
			this->v_linears_full_prec = (vec_mp *) br_malloc(other.num_v_linears*sizeof(vec_mp));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_mp(this->v_linears[ii],0);
				init_vec_mp2(this->v_linears_full_prec[ii],0,1024);
				vec_cp_mp(this->v_linears[ii],other.v_linears[ii]);
				vec_cp_mp(this->v_linears_full_prec[ii],other.v_linears_full_prec[ii]);
			}
		}
		else{}
		
		vec_cp_mp(this->v_patch, other.v_patch);
		vec_cp_mp(this->v_patch_full_prec, other.v_patch_full_prec);
		
		mat_cp_mp(this->jac_with_proj, other.jac_with_proj);
		mat_cp_mp(this->jac_with_proj_full_prec, other.jac_with_proj_full_prec);
		

		
		
		
		
		if (other.num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
			this->target_projection_full_prec = (vec_mp *) br_malloc(other.num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				init_vec_mp2(this->target_projection[ii],0,1024);
				vec_cp_mp(this->target_projection[ii],other.target_projection[ii]);
				vec_cp_mp(this->target_projection_full_prec[ii],other.target_projection_full_prec[ii]);
			}
		}
		else{}
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->target_dim = other.target_dim;   // r			the dimension of the real set we are looking for
		this->ambient_dim = other.ambient_dim;  // k			the dimension of the complex component we are looking IN.
		this->target_crit_codim = other.target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		this->num_v_vars = other.num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
		this->num_natural_vars = other.num_natural_vars;  // N-k+\ell
		this->num_synth_vars = other.num_synth_vars;
		
		this->max_degree = other.max_degree;						// the max degree of differentiated (randomized) functions
		
		
		this->num_projections = other.num_projections;
		this->num_v_linears = other.num_v_linears;
		this->num_additional_linears = other.num_additional_linears;
		
		
		
	} // re: copy
	
	
	
	
}; // re: class nullspace_eval_data_mp










// the double version
// this must be defined after the mp version, because double has mp.

/**
 \brief double-format of the left-nullspace evaluator data.
 
 \see nullspacejac_eval_data_mp
 */
class nullspacejac_eval_data_d : public SolverDoublePrecision
{
	int side_;
	
public:
	
	int side()
	{
		return side_;
	}
	
	nullspacejac_eval_data_mp *BED_mp; ///< pointer to the MP eval data. used only for AMP
	
	
	prog_deriv_t * SLP_derivative;
	
	
	int num_jac_equations;
	int target_dim;   // r			the dimension of the real set we are looking for
	int ambient_dim;  // k			the dimension of the complex component we are looking IN.
	int target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
	
	int num_natural_vars;  // N   number of variables in original problem statement (including homogenizing variables)
	int num_v_vars;  // (N-k) + (k-\ell+1)
	int num_synth_vars;
	
	int max_degree;						// the max degree of differentiated (randomized) functions
	
	int num_additional_linears;
	vec_d *additional_linears_terminal;
	
	vec_d *additional_linears_starting;
	
	
	
	vec_d **starting_linears; // outer layer should have as many as there are randomized equations
							  // inside layer has number corresponding to randomized_degrees
	
	
	int num_v_linears;
	vec_d *v_linears;         // should be as many in here as there are randomized equations
	
	vec_d v_patch;
	

	mat_d jac_with_proj;
	
	
	
	int num_projections;
	vec_d *target_projection; // # of these should be target_dim (for now)
	
	
	// default initializer
	nullspacejac_eval_data_d() : SolverDoublePrecision(){
		
		reset_counters();
		
		init();
	}
	
	nullspacejac_eval_data_d(int mp) : SolverDoublePrecision(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~nullspacejac_eval_data_d(){
		
		nullspacejac_eval_data_d::clear();
		// no need to reset the counters.
	}
	
	void print();
	
	
	void reset_counters()
	{
		num_jac_equations = 0;
		target_dim = 0;   // r			the dimension of the real set we are looking for
		ambient_dim = 0;  // k			the dimension of the complex component we are looking IN.
		target_crit_codim = 0;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		num_synth_vars = -1631;
		num_v_vars = -941;  // N   number of variables in original problem statement (including homogenizing variables)
		num_natural_vars = -746;  // N-k+\ell
		
		max_degree = -112;						// the max degree of differentiated (randomized) functions
		
		
		num_projections = 0;
		num_v_linears = 0;
		num_additional_linears = 0;
	} // re: reset_counters
	
	
	nullspacejac_eval_data_d & operator=(const nullspacejac_eval_data_d & other)
	{
		clear(); // this is wasteful, but safe for now
		reset_counters();// this is wasteful, but safe for now
		
		copy(other);
		return *this;
	}
	
	nullspacejac_eval_data_d(const nullspacejac_eval_data_d & other)
	{
		init();
		SolverDoublePrecision();
		nullspacejac_eval_data_d();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	/**
	 \brief MPI broadcast send to all workers in communicator, for the nullspace solver.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI.
	 */
	int send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief MPI broadcast recive from the head of the communicator, for the nullspace solver.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI.
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	
	
	void set_sidedness(int which_side);
	
	void set_function_handles(int which_side);
	
	
	
	/**
	 \brief populate the eval_data form the config object.
	 
	 \return SUCCESSFUL
	 \param _SLP a pointer to the underlying SLP.
	 \param ns_config a pointer to the NullspaceConfiguration from which to set up.
	 \param W input witness set containing the patches, etc.
	 \param solve_options The current state of the solver.
	 */
	int setup(prog_t * _SLP,
			  NullspaceConfiguration *ns_config,
			  WitnessSet & W,
			  SolverConfiguration & solve_options);
	
	
private:
	
	
	void clear()
	{
		
		
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_d(additional_linears_terminal[ii]);
			}
			free(additional_linears_terminal);
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				clear_vec_d(additional_linears_starting[ii]);
			}
			free(additional_linears_starting);
		}
		
		
		
		
		if (num_jac_equations>0) {
			if (side_ == nullspace_handedness::LEFT) {
				for (int ii=0; ii<num_jac_equations; ii++) {
					for (int jj=0; jj<max_degree; jj++) {
						clear_vec_d(starting_linears[ii][jj]);
					}
					free(starting_linears[ii]);
				}
			}
			else
			{
				for (int ii=0; ii<randomizer()->num_rand_funcs(); ii++) {
					for (int jj=0; jj<randomizer()->randomized_degree(ii)-1; jj++) {
						clear_vec_d(starting_linears[ii][jj]);
					}
					free(starting_linears[ii]);
				}
			}
			
			free(starting_linears);
		}
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				clear_vec_d(v_linears[ii]);
			}
			free(v_linears);
		}
		
		
		clear_vec_d(v_patch);
		clear_mat_d(jac_with_proj);
		

		
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_d(target_projection[ii]);
			}
			free(target_projection);
		}
		
		
		
		if (this->MPType==2) {
			delete this->BED_mp;
		}
		else{
			clear_deriv(SLP_derivative);
			delete this->SLP_derivative;
		}
	} // re: clear
	
	void init();
	
	void copy(const nullspacejac_eval_data_d & other)
	{
		
		clear();
		reset_counters();
		
		
		
		
		if (other.num_additional_linears>0) {
			
			
			this->additional_linears_terminal						= (vec_d *) br_malloc(other.num_additional_linears*sizeof(vec_d));
			this->additional_linears_starting						= (vec_d *) br_malloc(other.num_additional_linears*sizeof(vec_d));
			
			for (int ii=0; ii<other.num_additional_linears; ii++) {
				vec_cp_d(this->additional_linears_terminal[ii],						other.additional_linears_terminal[ii]);
				
				vec_cp_d(this->additional_linears_starting[ii],						other.additional_linears_starting[ii]);
			}
		}
		else {} // other.num_additional_linears == 0
		
		this->num_additional_linears = other.num_additional_linears;
		
		
		
		
		if (other.num_jac_equations) {
			this->starting_linears = (vec_d **) br_malloc(other.num_jac_equations*sizeof(vec_d *));
			
			for (int ii=0; ii<other.num_jac_equations; ii++) {
				this->starting_linears[ii] = (vec_d *) br_malloc(other.max_degree*sizeof(vec_d ));
				for (int jj=0; jj<other.max_degree; jj++) {
					init_vec_d(this->starting_linears[ii][jj],0);
					vec_cp_d(this->starting_linears[ii][jj],other.starting_linears[ii][jj]);
				}
			}
		}
		else{}
		
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->max_degree = other.max_degree;
		
		
		if (other.num_v_linears) {
			this->v_linears = (vec_d *) br_malloc(other.num_v_linears*sizeof(vec_d));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_d(this->v_linears[ii],0);
				vec_cp_d(this->v_linears[ii],other.v_linears[ii]);
			}
		}
		else{}
		
		vec_cp_d(this->v_patch, other.v_patch);
		
		mat_cp_d(this->jac_with_proj, other.jac_with_proj);
		
		
		
		
		if (other.num_projections>0) {
			this->target_projection = (vec_d *) br_malloc(other.num_projections*sizeof(vec_d));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_d(this->target_projection[ii],0);
				vec_cp_d(this->target_projection[ii],other.target_projection[ii]);
			}
		}
		else{}
		
		
		this->num_jac_equations = other.num_jac_equations;
		this->target_dim = other.target_dim;   // r			the dimension of the real set we are looking for
		this->ambient_dim = other.ambient_dim;  // k			the dimension of the complex component we are looking IN.
		this->target_crit_codim = other.target_crit_codim;    // \ell.  must be at least one (1), and at most target_dim (r).
		
		this->num_v_vars = other.num_v_vars;  // N   number of variables in original problem statement (including homogenizing variables)
		this->num_natural_vars = other.num_natural_vars;  // N-k+\ell
		this->num_synth_vars = other.num_synth_vars;
		
		this->max_degree = other.max_degree;						// the max degree of differentiated (randomized) functions
		
		
		this->num_projections = other.num_projections;
		this->num_v_linears = other.num_v_linears;
		this->num_additional_linears = other.num_additional_linears;
		
		
		
	} // re: copy
	
};








/** 
 \brief the main function for finding critical conditions WRT a projection
 
 \todo remove MPType as an input.
 
 \return SUCCESSFUL
 \param MPType The operating MPtype.
 \param W input witness set, carrying with it the start points.
 \param solve_out Computed solutions and metadata go here.
 \param ns_config The nullspace config setup, previously computed.
 \param solve_options The current state of the solver.
 */
int nullspacejac_solver_master_entry_point(int MPType,
										   WitnessSet & W,
										   SolverOutput & solve_out,
										   NullspaceConfiguration				*ns_config,
										   SolverConfiguration		& solve_options);









/**
 \brief Evaluator function for the left nullspace solver.
 
 \todo explain with diagram how this works
 
 this function makes use of the TemporariesMultiplePrecision class for persistence of temporaries.
 
 \return the number 0.
 \param funcVals the computed function values.
 \param parVals the computed parameter values.
 \param parDer the computed derivatives with respect to parameters.
 \param Jv the computed Jacobian with respect to the variables.
 \param Jp the computed Jacobian with respect to the parameters (time).
 \param current_variable_values The input variable values.
 \param pathVars the current time
 \param ED a pointer from which we type-cast, into the correct type.
 */
int nullspacejac_left_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);


/**
 \brief Evaluator function for the left nullspace solver.
 
 \see nullspacejac_eval_d
 
 this function makes use of the TemporariesMultiplePrecision class for persistence of temporaries.
 
 \return the number 0.
 \param funcVals the computed function values.
 \param parVals the computed parameter values.
 \param parDer the computed derivatives with respect to parameters.
 \param Jv the computed Jacobian with respect to the variables.
 \param Jp the computed Jacobian with respect to the parameters (time).
 \param current_variable_values The input variable values.
 \param pathVars the current time
 \param ED a pointer from which we type-cast, into the correct type.
 */
int nullspacejac_left_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);













/**
 \brief Evaluator function for the right nullspace solver.
 
 \todo explain with diagram how this works
 
 this function makes use of the TemporariesMultiplePrecision class for persistence of temporaries.
 
 \return the number 0.
 \param funcVals the computed function values.
 \param parVals the computed parameter values.
 \param parDer the computed derivatives with respect to parameters.
 \param Jv the computed Jacobian with respect to the variables.
 \param Jp the computed Jacobian with respect to the parameters (time).
 \param current_variable_values The input variable values.
 \param pathVars the current time
 \param ED a pointer from which we type-cast, into the correct type.
 */
int nullspacejac_right_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);





/**
 \brief Evaluator function for the right nullspace solver.
 
 \see nullspacejac_right_eval_d
 
 this function makes use of the TemporariesMultiplePrecision class for persistence of temporaries.
 
 \return the number 0.
 \param funcVals the computed function values.
 \param parVals the computed parameter values.
 \param parDer the computed derivatives with respect to parameters.
 \param Jv the computed Jacobian with respect to the variables.
 \param Jp the computed Jacobian with respect to the parameters (time).
 \param current_variable_values The input variable values.
 \param pathVars the current time
 \param ED a pointer from which we type-cast, into the correct type.
 */
int nullspacejac_right_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);
















/**
 \brief dehomogenization method for the left nullspace jacobian solver.
 
 returned (filled) type inferred by in_prec -- if in_prec<64, populate the double, else the mp.
 
 This function fits the format for all Bertini dehomogenizers.
 
 \return the number 0.
 \param out_d returned double values
 \param out_mp returned mp values, after dehomogenization
 \param out_prec the precision of the output, and is set to = in_prec.
 \param in_d input in double format, should only be populated if in_prec<64.
 \param in_mp input point in mp format, populated if in_rec >= 64.
 \param in_prec precision of the input point.
 \param ED_d input evaluator, needed to get some other parameters.
 \param ED_mp input evaluator, needed to get some other parameters.
 */
int nullspacejac_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);



/**
 \brief change the precision of an MP nullspace evaluator.
 
 This function fits the format for all Bertini precision changers.
 
 \return the number 0.
 \param ED pointer to the evaluator data to change.
 \param new_prec the precision to change to.
 */
int change_nullspacejac_eval_prec(void const *ED, int new_prec);



/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_nullspacejac_d(endgame_data_t *EG,
								tracker_config_t *T,
								void const *ED);

/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \see check_issoln_nullspacejac_d
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_nullspacejac_mp(endgame_data_t *EG,
								 tracker_config_t *T,
								 void const *ED);



/**
 \brief check whether an input point is a valid start point, by residuals of the functions.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param testpoint the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the nullspace eval data in double form.
 */
int check_isstart_nullspacejac_d(vec_d testpoint,
								 tracker_config_t *T,
								 void const *ED);

/**
 \brief check whether an input point is a valid start point, by residuals of the functions.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param testpoint the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the nullspace eval data in double form.
 */
int check_isstart_nullspacejac_mp(vec_mp testpoint,
								  tracker_config_t *T,
								  void const *ED);



/**
 \brief a small function for verifying the homogeneousness of the evaluator (for correctness), etc.
 
 This function has no real purpose except for development.
 
 \param current_values the input
 \param ED the midtrack eval data.
 */
void check_nullspace_evaluator(point_mp current_values,
							   void const *ED);



/**
 \brief the no-output entry point for a worker to help with a nullspace solve.  
 
 this function is called *after* the worker has received the call-for-help broadcast from the master.
 \param solve_options The current state of the solver.
 */
void nullspace_slave_entry_point(SolverConfiguration & solve_options);







#endif


