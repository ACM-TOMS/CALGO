
#ifndef _MULTILINSOLVER_H
#define _MULTILINSOLVER_H

/** \file solver_multilintolin.hpp */

#include "bertini1/bertini_headers.hpp"



#include "nag/solvers/solver.hpp"
#include "io/fileops.hpp"
#include "programConfiguration.hpp"
#include "nag/solvers/postProcessing.hpp"






/**
 \brief config class for the multilin solver, so that SLP can be set up independently from the eval data.
 */
class MultilinConfiguration
{
	
	std::shared_ptr<SystemRandomizer> randomizer_; ///< pointer to randomizer
	
public:
	
	/**
	 \brief get a shared pointer to the randomizer
	 
	 \return a shared pointer to the randomizer
	 */
	std::shared_ptr<SystemRandomizer> randomizer()
	{
		return randomizer_;
	}
	
	
	/**
	 \brief set up this object's randomizer to point to the same place as the input randomizer.
	 
	 \param _random a pointer to the randomizer this should point to.
	 */
	void set_randomizer(std::shared_ptr<SystemRandomizer> _random)
	{
		randomizer_ = _random;
	}
	
	
	
	
	StraightLineProgramGlobalPointers SLP_memory; ///< the memory for the SLP
	prog_t * SLP; ///< a pointer to the SLP for the solver
	
	int MPType; ///< M.O.
	
	bool have_mem; ///< do we have the memory setup?
	
	
	MultilinConfiguration(SolverConfiguration & solve_options,
					const WitnessSet & W)
	{
		init();
		
		set_memory(solve_options);
		
		make_randomizer(solve_options, W);
	}
	
	
	
	MultilinConfiguration(SolverConfiguration & solve_options,
					std::shared_ptr<SystemRandomizer> _random)
	{
		init();
		
		set_memory(solve_options);
		
		set_randomizer(_random);
	}
	
	
	MultilinConfiguration(SolverConfiguration & solve_options)
	{
		init();
		set_memory(solve_options);
	}
	
	
	MultilinConfiguration(std::shared_ptr<SystemRandomizer> _random)
	{
		init();
		set_randomizer(_random);
	}
	
	
	/**
	 \brief in the setup of this object, make a new system randomizer.
	 
	 \param solve_options the current state of the solver.
	 \param W input witness set, witn numbers of linears, patches, and variables.
	 */
	void make_randomizer(const SolverConfiguration & solve_options, const WitnessSet & W)
	{
		
		randomizer_ = std::make_shared<SystemRandomizer>(*(new SystemRandomizer));
		randomizer_->setup(W.num_variables()-W.num_linears()-W.num_patches(), solve_options.PPD.num_funcs);
		
	}
	
	
	/**
	 \brief set up the program, and capture the memory, followed by setting up the PPD.
	 \param solve_options The current state of the solver.
	 */
	void set_memory(SolverConfiguration & solve_options)
	{
		
		//TODO: should i assume here that the input file is already parsed??
		this->MPType = solve_options.T.MPType;
		solve_options.T.numVars = setupProg(SLP, solve_options.T.Precision, solve_options.T.MPType);
		//make randomizer matrix here
		SLP_memory.capture_globals();
		SLP_memory.set_globals_null();
		have_mem = true;
		
		
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
	}
	
	
	
	
	
	
	/**
	 \brief a mediocre function, which purports to copy from one MultilinConfiguration to another.  this is broken.
	 \param other The MultilinConfiguration from which to copy.
	 */
	void copy(const MultilinConfiguration & other)
	{
		this->randomizer_ = other.randomizer_;
		//TODO: write this function
		//ah yes, this problem.
		//		this->SLP = other.SLP;// this needs to be a deep copy
	}
	
	void init()
	{

		SLP = new prog_t;
		have_mem = false;
		
	}
	
	
	void clear()
	{
		
		
		SLP_memory.set_globals_to_this();
		clearProg(SLP, this->MPType, 1); // 1 means call freeprogeval()
		delete SLP;
	}
	
	
	
	MultilinConfiguration()
	{
		init();
	}
	
	MultilinConfiguration(const MultilinConfiguration & other)
	{
		init();
		copy(other);
	}
	
	MultilinConfiguration & operator=(const MultilinConfiguration & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	
	~MultilinConfiguration(){
		clear();
	}
};

















/**
\brief the evaluator data type for the multilin solver.
 */
// this must be defined before the double version, because double has mp.
class multilintolin_eval_data_mp : public SolverMultiplePrecision
{
public:
	
	
	int num_linears; ///< how many linears it currently has.
	
	
	
	vec_mp *current_linear;		///< the linear we are moving to, in current precision
	vec_mp *current_linear_full_prec; ///< the linear we are moving to, in full precision
	vec_mp *old_linear;				///< the linear we are moving away from, in current precision
	vec_mp *old_linear_full_prec;			///< the linear we are moving away from, in full precision
	
	
	
	// default initializer
	multilintolin_eval_data_mp() : SolverMultiplePrecision(){
		std::cout << "instantiating multilin eval d without declaring mp type"	 << std::endl;
		mypause();
		init();
	}
	
	multilintolin_eval_data_mp(int mp) : SolverMultiplePrecision(mp){
		this->MPType = mp;
		init();
	}
	
	
	multilintolin_eval_data_mp(const multilintolin_eval_data_mp & other) : SolverMultiplePrecision(other)
	{
		this->MPType = other.MPType;
		init();
		copy(other);
	}
	
	multilintolin_eval_data_mp & operator=(const multilintolin_eval_data_mp & other)
	{
		this->MPType = other.MPType;
		copy(other);
		return *this;
	}
	
	
	
	virtual ~multilintolin_eval_data_mp(){
		multilintolin_eval_data_mp::clear();
	}
	
	
	
	virtual void print()
	{
		
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	
	/**
	 \brief bcast send for the multilin solver
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	/**
	 \brief bcast receive for the multilin solver
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief set up in memory the system for solving, from the config object passed in, and the witness set.
	 
	 \return The number 0.
	 \param config The config object for this evaluator, containing the SLP etc.
	 \param W The input witness set containing start points
	 \param target_linears The linears we move to.
	 \param solve_options The current state of the solver.
	 */
	int setup(MultilinConfiguration & config,
			  const WitnessSet & W,
			  vec_mp * target_linears,
			  SolverConfiguration & solve_options);
	
	
	
protected:
	
	
	void clear()
	{
		
		
		for (int ii=0; ii<num_linears; ii++) {
			clear_vec_mp(current_linear[ii]);
			clear_vec_mp(old_linear[ii]);
		}
		free(current_linear);
		free(old_linear);
		
		
		if (this->MPType == 2) {
			for (int ii=0; ii<num_linears; ii++) {
				clear_vec_mp(current_linear_full_prec[ii]);
				clear_vec_mp(old_linear_full_prec[ii]);
			}
			free(current_linear_full_prec);
			free(old_linear_full_prec);
		}
		
		
		
	} // re: clear
	
	void init();
	
	
	void copy(const multilintolin_eval_data_mp & other)
	{
		SolverMultiplePrecision::copy(other);
		
		
		if (this->num_linears==0) {
			current_linear = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
			old_linear = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
		}
		else
		{
			current_linear = (vec_mp *) br_realloc(current_linear, other.num_linears*sizeof(vec_mp));
			old_linear = (vec_mp *) br_realloc(old_linear, other.num_linears*sizeof(vec_mp));
		}
		
		for (int ii=0; ii<other.num_linears; ii++) {
			init_vec_mp(current_linear[ii],0);
			init_vec_mp(old_linear[ii],0);
			vec_cp_mp(current_linear[ii],other.current_linear[ii]);
			vec_cp_mp(old_linear[ii],other.old_linear[ii]);
		}
		
		
		
		if (this->MPType==2) {
			if (this->num_linears==0) {
				current_linear_full_prec = (vec_mp *) br_malloc(other.num_linears*sizeof(vec_mp));
				old_linear_full_prec = (vec_mp *) br_malloc(num_linears*sizeof(vec_mp));
			}
			else
			{
				current_linear_full_prec = (vec_mp *) br_realloc(current_linear_full_prec, other.num_linears*sizeof(vec_mp));
				old_linear_full_prec = (vec_mp *) br_realloc(old_linear_full_prec, other.num_linears*sizeof(vec_mp));
			}
			
			for (int ii=0; ii<other.num_linears; ii++) {
				init_vec_mp2(current_linear_full_prec[ii],0,1024);
				init_vec_mp2(old_linear_full_prec[ii],0,1024);
				vec_cp_mp(current_linear_full_prec[ii],other.current_linear_full_prec[ii]);
				vec_cp_mp(old_linear_full_prec[ii],other.old_linear_full_prec[ii]);
			}
		}
		
		this->num_linears= other.num_linears;
	} // re: copy
	
	
	
}; // re: class nullspace_eval_data_mp














// the mp version
// this must be defined before the double version, because double has mp.
class multilintolin_eval_data_d : public SolverDoublePrecision
{
public:
	
	multilintolin_eval_data_mp * BED_mp; ///< a pointer to the MP eval data, for AMP mode
	int num_linears; ///< the number of linears.
	
	
	
	vec_d *current_linear;///< the linears we move TO
	vec_d *old_linear; ///< the linears we move away FROM.
	
	
	
	// default initializer
	multilintolin_eval_data_d() : SolverDoublePrecision(){
		std::cout << "instantiating multilin eval d without declaring mp type"	 << std::endl;
		mypause();
		init();
	}
	
	multilintolin_eval_data_d(int mp) : SolverDoublePrecision(mp)
	{
		this->MPType = mp;
		
		init();
	}
	
	
	multilintolin_eval_data_d(const multilintolin_eval_data_d & other) : SolverDoublePrecision(other)
	{
		this->MPType = other.MPType;
		init();
		copy(other);
	}
	
	multilintolin_eval_data_d & operator=(const multilintolin_eval_data_d & other)
	{
		this->MPType = other.MPType;
		copy(other);
		return *this;
	}
	
	
	
	virtual ~multilintolin_eval_data_d(){
		clear();
	}
	
	
	/**
	 \brief print useful information to the screen
	 */
	virtual void print()
	{
		SolverDoublePrecision::print();
		
		std::cout << "multilintolin evaluator data (double):" << std::endl;
		for (int ii=0; ii<num_linears; ii++) {
			
			std::stringstream name;
			name << "old_linear_" << ii;
			print_point_to_screen_matlab(old_linear[ii],name.str());
		}
		
		for (int ii=0; ii<num_linears; ii++) {
			
			std::stringstream name;
			name << "current_linear_" << ii;
			print_point_to_screen_matlab(current_linear[ii],name.str());
		}
	}
	
	
	
	
	
	// MPI SENDS AND RECEIVES
	/**
	 \brief bcast send for the multilin solver
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief bcast receive for the multilin solver
	 
	 \return SUCCESSFUL
	 \param mpi_config The current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	

	/**
	 \brief set up in memory the system for solving, from the config object passed in, and the witness set.
	 
	 \return The number 0.
	 \param config The config object for this evaluator, containing the SLP etc.
	 \param W The input witness set containing start points
	 \param target_linears The linears we move to.
	 \param solve_options The current state of the solver.
	 */
	int setup(MultilinConfiguration & config,
			  const WitnessSet & W,
			  vec_mp * target_linears,
			  SolverConfiguration & solve_options);
	
	
	
protected:
	void init();
	
	void clear()
	{
		
		for (int ii=0; ii<num_linears; ii++) {
			clear_vec_d(current_linear[ii]);
			clear_vec_d(old_linear[ii]);
		}
		free(current_linear);
		free(old_linear);
		
		if (this->MPType==2)
		{
			delete this->BED_mp;
		}
	}
};



/**
 \brief The main way to get either the head or a serial process into moving linears around.
 
 \return SUCCESSFUL
 \param W the input witness set including start points, starting linears, and patches.
 \param solve_out computed results and metadata go here.
 \param new_linears the linears to which we wish to move.
 \param config the multilin config object, which is passed in as argument to allow more efficient management of data.
 \param solve_options the current state of the solver.
 */
int multilin_solver_master_entry_point(const WitnessSet & W, // carries with it the start points, and the linears.
									   SolverOutput & solve_out, // new data goes in here
									   vec_mp * new_linears,
									   MultilinConfiguration &		config,
									   SolverConfiguration		& solve_options);



/**
 \brief how to get a worker to cooperate to move linears around.
 
 a slave comes to this function empty-handed and leaves empty-handed, but helps do the work to move from one set of linears to another, by tracking some of the paths and sending the results to the head.
 
 \return SUCCESSFUL
 \param solve_options the current state of the solver config.
 */
int multilin_slave_entry_point(SolverConfiguration & solve_options);





/**
 \brief Evaluator function for the multilin solver.
 
 \todo explain with diagram how this works
 
 this function makes use of the TemporariesMultiplePrecision class for persistence of temporaries.
 
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
int multilin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);



/**
 \brief Evaluator function for the multilin solver.
 
 \see multilin_to_lin_eval_d
 
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
int multilin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);



/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_multilintolin_d(endgame_data_t *EG,
								 tracker_config_t *T,
								 void const *ED);

/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_multilintolin_mp(endgame_data_t *EG,
								  tracker_config_t *T,
								  void const *ED);




/**
 \brief dehomogenization method for this solver.
 

 assumes have only a single non-homogeneous variable group with a leading homogenizing coordinate, and dehomogenizes.
 
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
int multilintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);




/**
 \brief change the precision of an MP multilin evaluator.
 
 This function fits the format for all Bertini precision changers.
 
 \return the number 0.
 \param ED pointer to the evaluator data to change.
 \param new_prec the precision to change to.
 */
int change_multilintolin_eval_prec(void const *ED, int new_prec);









#endif


