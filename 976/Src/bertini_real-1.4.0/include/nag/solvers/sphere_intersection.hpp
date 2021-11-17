#ifndef _SPHERESOLVER_H
#define _SPHERESOLVER_H


/** \file solver_sphere.hpp */


#include "bertini1/bertini_headers.hpp"

#include "nag/solvers/solver.hpp"
#include "io/fileops.hpp"
#include "programConfiguration.hpp"
#include "nag/solvers/postProcessing.hpp"






/**
 \brief The config class for the sphere intersection solver.
 
 This config class is for intersection of a sphere with another system. The config holds the system, randomizer, and parameters such as the center and radius.  It is passed in to the main() method, with a witness set containing the start points, whereafter the master calls for help and the paths are tracked.
 
 \see nullspace_config
 \see midpoint_config
 \see multilin_config
 */
class SphereConfiguration
{
	std::shared_ptr<SystemRandomizer> randomizer_; ///< randomizer for the main system being intersected.

public:
	
	
	/**
	 \brief get a shared pointer to the randomizer
	 
	 \return a shared pointer to the randomizer
	 */
	std::shared_ptr<SystemRandomizer> randomizer()
	{
		return randomizer_;
	}
	
	
	
	StraightLineProgramGlobalPointers SLP_memory; ///< the memory containing the temps for the evaluation of the SLP
	prog_t * SLP; ///< a pointer to the SLP in memory
	bool have_mem; ///< whether have the memory set up for SLP
	
	
	int MPType; ///< current operating MP mode
	
	
	
	
	vec_mp *starting_linear;  ///< the set of starting linears for the linear product.  There should damn well be two of them.
	
	
	
	vec_mp center; ///< the center of the sphere in complex coordinates.
	comp_mp radius;///< the radius of the sphere.
	
	
	
	SphereConfiguration(SolverConfiguration & solve_options,
				  const WitnessSet & W)
	{
		init();
		
		set_memory(solve_options);
		
		make_randomizer(solve_options, W);
	}
	
	
	
	SphereConfiguration(SolverConfiguration & solve_options,
				  std::shared_ptr<SystemRandomizer> _random)
	{
		init();
		
		set_memory(solve_options);
		
		set_randomizer(_random);
	}
	
	
	SphereConfiguration(SolverConfiguration & solve_options)
	{
		init();
		set_memory(solve_options);
	}
	
	
	SphereConfiguration(std::shared_ptr<SystemRandomizer> _random)
	{
		init();
		set_randomizer(_random);
	}
	
	/**
	 \brief have the config make its own randomizer for the supplied system.
	 
	 \param solve_options The current state of the solver.
	 \param W the input witness set with the patches, and numbers of variables and linears.
	 */
	void make_randomizer(const SolverConfiguration & solve_options, const WitnessSet & W)
	{
		randomizer_ = std::make_shared<SystemRandomizer>(*(new SystemRandomizer));
		randomizer_->setup(W.num_variables()-W.num_linears()-W.num_patches(), solve_options.PPD.num_funcs);
	}
	
	/**
	 \brief Get the SLP, and capture the global memory pointers.  Get the PPD.
	 
	 the input file must already be parsed out in the current working folder.
	 
	 \param solve_options The current state of the solver.
	 */
	void set_memory(SolverConfiguration & solve_options);
	
	
	/**
	 \brief set this randomizer to point to the one passed in.
	 
	 \param _random Pointer to the randomizer to use.
	 */
	void set_randomizer(std::shared_ptr<SystemRandomizer> _random)
	{
		randomizer_ = _random;
	}
	
	/**
	 \brief Set the center of the sphere.
	 
	 \param new_center The center of the sphere.
	 */
	void set_center(vec_mp new_center)
	{
		vec_cp_mp(center, new_center);
	}
	
	
	/**
	 \brief set the radius of the sphere.
	 
	 \param new_radius The new radius of the sphere.
	 */
	void set_radius(comp_mp new_radius)
	{
		set_mp(this->radius, new_radius);
	}
	
	
	
	
	
	SphereConfiguration()
	{
		init();
	}
	
	SphereConfiguration(const SphereConfiguration & other)
	{
		init();
		copy(other);
	}
	
	SphereConfiguration & operator=(const SphereConfiguration & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	
	~SphereConfiguration(){
		clear();
	}
	
protected:
	
	
	void copy(const SphereConfiguration & other)
	{
		this->randomizer_ = other.randomizer_;
		//TODO: write this function.
		//ah yes, this problem.
		//		this->SLP = other.SLP;// this needs to be a deep copy
	}
	
	void init()
	{
		
		SLP = new prog_t;
		have_mem = false;

		
		starting_linear = (vec_mp *) br_malloc(2*sizeof(vec_mp));
		init_vec_mp2(starting_linear[0],0,1024);
		init_vec_mp2(starting_linear[1],0,1024);
		
		init_vec_mp2(center,0,1024);
		init_mp2(radius,1024);
	}
	
	
	void clear()
	{
		clear_vec_mp(center);
		clear_mp(radius);
		
		
		clear_vec_mp(starting_linear[0]);
		clear_vec_mp(starting_linear[1]);
		free(starting_linear);
		
		
		SLP_memory.set_globals_to_this();
		clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete SLP;
	}
};



















// the mp version
// this must be defined before the double version, because double has mp.

/**
 \brief The multiple-precision evaluator data for the sphere intersection solver.
 */
class sphere_eval_data_mp : public SolverMultiplePrecision
{
public:
	
    int num_natural_vars; ///< the number of non-synthetic variables, including the homogenizing variable
    
	comp_mp two, two_full_prec; ///< the number two.
	
	//there had better be two starting linears.
	vec_mp *starting_linear;					///< the two starting linears for start system
	vec_mp *starting_linear_full_prec;			///< the two starting linears for start system, in full precision for AMP.
	
	
	
	vec_mp *static_linear; ///< additional linears to NOT move.  These are linear slices.
	vec_mp *static_linear_full_prec;///< additional linears to NOT move.  These are linear slices.  In full precision for AMP
	int num_static_linears; ///< how many static linears there are.
	
	
	vec_mp center, center_full_prec;  ///< the center of the sphere
	comp_mp radius, radius_full_prec; ///< the radius of the sphere
	
	
	// default initializer
	sphere_eval_data_mp() : SolverMultiplePrecision(){
		init();
	}
	
	sphere_eval_data_mp(int mp) : SolverMultiplePrecision(mp){
		this->MPType = mp;
		init();
	}
	
	
	sphere_eval_data_mp(const sphere_eval_data_mp & other) : SolverMultiplePrecision(other)
	{
		init();
		copy(other);
	}
	
	sphere_eval_data_mp & operator=(const sphere_eval_data_mp & other)
	{
		
		init();
		copy(other);
		return *this;
	}
	
	
	
	virtual ~sphere_eval_data_mp(){
		sphere_eval_data_mp::clear();
	}
	
	
	
	virtual void print()
	{
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	/**
	 \brief MPI-broadcast send for the eval data.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	/**
	 \brief MPI-broadcast receive for the eval data.
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief set up eval data from a config object, a witness set, and the current solver state.
	 
	 \return the number 0.
	 \param config The sphere config object, with the system, radius, center, etc.
	 \param W witness set containing the static linears.
	 \param solve_options The current state of the solver.
	 */
	int setup(SphereConfiguration & config,
			  const WitnessSet & W,
			  SolverConfiguration & solve_options);
	
	
	
protected:
	
	
	void clear()
	{
		clear_mp(two);
		clear_vec_mp(center);
		clear_mp(radius);
		
		for (int ii=0; ii<2; ii++) {
			clear_vec_mp(starting_linear[ii]);
		}
		free(starting_linear);
		
		if (num_static_linears>0) {
			for (int ii=0; ii<num_static_linears; ii++) {
				clear_vec_mp(static_linear[ii]);
			}
			
			free(static_linear);
		}
		
		
		
		
		if (this->MPType == 2) {
			clear_vec_mp(center_full_prec);
			clear_mp(radius_full_prec);
			clear_mp(two_full_prec);
			
			for (int ii=0; ii<2; ii++) {
				clear_vec_mp(starting_linear_full_prec[ii]);
			}
			free(starting_linear_full_prec);
			
			if (num_static_linears>0) {
				for (int ii=0; ii<num_static_linears; ii++) {
					clear_vec_mp(static_linear_full_prec[ii]);
				}
				
				free(static_linear_full_prec);
			}
			
			
		}
		
		
		
	} // re: clear
	
	void init();
	
	
	void copy(const sphere_eval_data_mp & other)
	{
		SolverMultiplePrecision::copy(other);
		
		
		if (this->num_static_linears==0) {
			static_linear = (vec_mp *) br_malloc(other.num_static_linears*sizeof(vec_mp));
		}
		else
		{
			static_linear = (vec_mp *) br_realloc(static_linear, other.num_static_linears*sizeof(vec_mp));
		}
		
		for (int ii=0; ii<other.num_static_linears; ii++) {
			init_vec_mp(static_linear[ii],0);
			vec_cp_mp(static_linear[ii],other.static_linear[ii]);
		}
		
		for (int ii=0; ii<2; ii++) {
			vec_cp_mp(starting_linear[ii], other.starting_linear[ii]);
		}
		
		set_mp(this->radius, other.radius);
		vec_cp_mp(this->center, other.center);
		
		if (this->MPType==2) {
			
			set_mp(this->radius_full_prec, other.radius_full_prec);
			vec_cp_mp(this->center_full_prec, other.center_full_prec);
			
			if (this->num_static_linears==0) {
				static_linear_full_prec = (vec_mp *) br_malloc(other.num_static_linears*sizeof(vec_mp));
			}
			else
			{
				static_linear_full_prec = (vec_mp *) br_realloc(static_linear_full_prec, other.num_static_linears*sizeof(vec_mp));
			}
			
			for (int ii=0; ii<other.num_static_linears; ii++) {
				init_vec_mp2(static_linear_full_prec[ii],0,1024);
				vec_cp_mp(static_linear_full_prec[ii],other.static_linear_full_prec[ii]);
			}
			
			for (int ii=0; ii<2; ii++) {
				vec_cp_mp(starting_linear_full_prec[ii], other.starting_linear_full_prec[ii]);
			}
		}
		
		this->num_static_linears = other.num_static_linears;
	} // re: copy
	
	
	
}; // re: class nullspace_eval_data_mp














// the mp version
// this must be defined before the double version, because double has mp.
class sphere_eval_data_d : public SolverDoublePrecision
{
public:
	
    int num_natural_vars; ///< the number of non-synthetic variables, including the homogenizing variable
    
	sphere_eval_data_mp * BED_mp; ///< a pointer to the MP eval data, for AMP mode.
	
	vec_d *starting_linear;	///< the two starting linears for start system
	
	comp_d two; ///< the number two
	
	vec_d *static_linear; ///< unmoving linear slices.
	int num_static_linears; ///< how many unmoving linear slices there are.
	
	
	vec_d center; ///< the center of the sphere.
	comp_d radius; ///< the radius of the sphere.
	
	
	
	// default initializer
	sphere_eval_data_d() : SolverDoublePrecision(){
		init();
	}
	
	sphere_eval_data_d(int mp) : SolverDoublePrecision(mp)
	{
		this->MPType = mp;
		init();
	}
	
	
	sphere_eval_data_d(const sphere_eval_data_d & other) : SolverDoublePrecision(other)
	{
		init();
		copy(other);
	}
	
	sphere_eval_data_d & operator=(const sphere_eval_data_d & other)
	{
		copy(other);
		return *this;
	}
	
	
	
	virtual ~sphere_eval_data_d(){
		sphere_eval_data_d::clear();
	}
	
	
	
	virtual void print()
	{
		SolverDoublePrecision::print();
		
		std::cout << "sphere evaluator data (double):" << std::endl;
		for (int ii=0; ii<num_static_linears; ii++) {
			
			std::stringstream name;
			name << "static_linear_" << ii;
			print_point_to_screen_matlab(static_linear[ii],name.str());
		}
		
		for (int ii=0; ii<2; ii++) {
			
			std::stringstream name;
			name << "starting_linear_" << ii;
			print_point_to_screen_matlab(starting_linear[ii],name.str());
		}
		
		print_point_to_screen_matlab(center,"center");
		
		print_comp_matlab(radius,"radius");
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	/**
	 \brief MPI-broadcast send for the eval data.
	 
	 \param mpi_config The current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief MPI-broadcast receive for the eval data.
	 
	 \param mpi_config The current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief set up eval data from a config object, a witness set, and the current solver state.
	 
	 \return the number 0.
	 \param config The sphere config object, with the system, radius, center, etc.
	 \param W witness set containing the static linears.
	 \param solve_options The current state of the solver.
	 */
	int setup(SphereConfiguration & config,
			  const WitnessSet & W,
			  SolverConfiguration & solve_options);
	
	
	
protected:
	void init();
	
	void clear()
	{
		
		clear_vec_d(center);
		
		for (int ii=0; ii<2; ii++) {
			clear_vec_d(starting_linear[ii]);
		}
		free(starting_linear);
		
		if (num_static_linears>0) {
			for (int ii=0; ii<num_static_linears; ii++) {
				clear_vec_d(static_linear[ii]);
			}
			
			free(static_linear);
		}
		
		
		
		if (this->MPType==2)
		{
			delete this->BED_mp;
		}
	}
	
	
	
	void copy(const sphere_eval_data_d & other)
	{
		SolverDoublePrecision::copy(other);
		
		
		if (this->num_static_linears==0) {
			static_linear = (vec_d *) br_malloc(other.num_static_linears*sizeof(vec_d));
		}
		else
		{
			static_linear = (vec_d *) br_realloc(static_linear, other.num_static_linears*sizeof(vec_d));
		}
		
		for (int ii=0; ii<other.num_static_linears; ii++) {
			init_vec_d(static_linear[ii],0);
			vec_cp_d(static_linear[ii],other.static_linear[ii]);
		}
		
		for (int ii=0; ii<2; ii++) {
			vec_cp_d(starting_linear[ii], other.starting_linear[ii]);
		}
		
		set_d(this->radius, other.radius);
		vec_cp_d(this->center, other.center);
	}
	
};


/**
 \brief Main entry point for intersecting a sphere with another system.
 
 To use this function, you must already have two ready-to-move generic linear slices for a positive dimensional component.
 
 \return SUCCESSFUL
 \param W input witness set with start points and linears.
 \param solve_out The returned constructed data and metadata goes here.
 \param config The SphereConfiguration object with the SLP, etc.
 \param solve_options The current state of the solver.
 */
int sphere_solver_master_entry_point(const WitnessSet & W,
									 SolverOutput & solve_out,
									 SphereConfiguration & config,
									 SolverConfiguration & solve_options);


/**
 \brief Have workers help track paths in a sphere intersection.
 
 This function has no returned data, and is called after the workers have received the call for help from the head.
 
 \return SUCCESSFUL
 \param solve_options The current state of MPI
 */
int sphere_slave_entry_point(SolverConfiguration & solve_options);




/**
 \brief Evaluator function for the left nullspace solver.
 
 \todo explain with diagram how this works
 
 this function makes use of the temps_d class for persistence of temporaries.
 
 \return the number 0.
 \param funcVals the computed function values.
 \param parVals the computed parameter values.
 \param parDer the computed derivatives with respect to parameters.
 \param Jv the computed Jacobian with respect to the variables.
 \param Jp the computed Jacobian with respect to the parameters (time).
 \param vars The input variable values.
 \param pathVars the current time
 \param ED a pointer from which we type-cast, into the correct type.
 */
int sphere_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);



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
int sphere_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);


/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_sphere_d(endgame_data_t *EG,
						  tracker_config_t *T,
						  void const *ED);

/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_sphere_mp(endgame_data_t *EG,
						   tracker_config_t *T,
						   void const *ED);


/**
 \brief dehomogenization method for the sphere intersection solver.
 
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
int sphere_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);



/**
 \brief change the precision of an MP sphere intersection evaluator.
 
 This function fits the format for all Bertini precision changers.
 
 \return the number 0.
 \param ED pointer to the evaluator data to change.
 \param new_prec the precision to change to.
 */
int change_sphere_eval_prec(void const *ED, int new_prec);










#endif


