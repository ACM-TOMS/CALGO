#ifndef SOLVER_MIDPOINT_H
#define SOLVER_MIDPOINT_H

/** 
 \file solver_midpoint_tracker.hpp 
 
 contains the classes and method for the triple-system tracker, midpoint_tracker.
 */
 
#include "bertini1/bertini_headers.hpp"



#include "nag/solvers/solver.hpp"
#include "io/fileops.hpp"
#include "programConfiguration.hpp"
#include "nag/solvers/postProcessing.hpp"


//forward declarations:
class MidpointConfiguration;
class midpoint_eval_data_mp;
class midpoint_eval_data_d;










class Surface; // forward declaration



/**
 \brief configuration class for the midpoint solver.
 
 */
class MidpointConfiguration
{
public:
	friend class CompleteSystem;
	
	std::map<std::string, CompleteSystem> systems; ///< the systems available for tracking.
	
	
    
    int MPType; ///< the operating MP type
    
    
	
	comp_mp v_target; ///< set during the loop in connect the dots, the target value to move to for v
	comp_mp u_target; ///< set during the loop in connect the dots, the target value to move to for u
	
	comp_mp crit_val_left; ///< set during the loop in connect the dots, the projection value of the right critslice
	comp_mp crit_val_right; ///< set during the loop in connect the dots, the projection value of the right critslice
	
	std::string system_name_mid; ///< the current name of the mid system -- acts as index into systems
	std::string system_name_bottom;///< the current name of the  bottom system -- acts as index into systems
	std::string system_name_top; ///< the current name of the top system -- acts as index into systems
	
	
	vec_mp *pi; ///< the projections used in the decomposition
	int num_projections; ///< the number of projections.  should match the dimension.
    
    
	
	MidpointConfiguration(){
		init();
	}
	
	~MidpointConfiguration(){
		clear();
	}
	
	
	MidpointConfiguration & operator=(const MidpointConfiguration & other)
	{
		init();
		
		copy(other);
		return *this;
	}
	
	MidpointConfiguration(const MidpointConfiguration & other)
	{
		init();
		
		copy(other);
		
	}
	
	

	
	/**
	 \brief display core information to the screen
	 
	 */
    void print()
	{

		print_comp_matlab(u_target,"u_target");
		print_comp_matlab(v_target,"v_target");
		
		print_comp_matlab(crit_val_right,"crit_val_right");
		print_comp_matlab(crit_val_left,"crit_val_left");
		

	}
	
	/**
	 \brief set up the data fields from the input surface decomposition.
	 
	 \param surf the surface we will connect dots on, to build faces.
	 \param solve_options the current state of the solver.
	 */
	void setup(const Surface & surf,
               SolverConfiguration & solve_options);
    
	
	/**
	 \brief broadcast mpi send this config to everyone in the communicator.
	 
	 \param mpi_config the current state of mpi
	 */
    void bcast_send(ParallelismConfig & mpi_config);
	
	/**
	 \brief broadcast receive a config from headnode in mpi communicator.
	 
	 \param mpi_config the current state of mpi.
	 
	 */
    void bcast_receive(ParallelismConfig & mpi_config);

	
	/**
	 \brief add a projection to this config.
	 
	 \param proj the projection vector to add.
	 */
    void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else {
			this->pi = (vec_mp *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_mp));
		}
		
		init_vec_mp2(pi[num_projections],0,1024);
		vec_cp_mp(pi[num_projections], proj);
		
		num_projections++;
		
	}
	
	/**
	 \brief get the number of variables in the top system
	 
	 \return the number of variables in the top system.
	 */
	int num_top_vars()
	{
		return systems[system_name_top].num_variables();
	}
	
	
	/**
	 \brief get the number of variables in the bottom edge's system.
	 
	 \return the number of variables in the bottom system
	 */
	int num_bottom_vars()
	{
		return systems[system_name_bottom].num_variables();
	}
	
	
	/**
	 \brief get the number of variables for the middle of the face.
	 
	 \return the number of variables for the middle of the face.
	 */
	int num_mid_vars()
	{
		return systems[system_name_mid].num_variables();
	}
	
	
	
	
	
	
private:
    
	void copy(const MidpointConfiguration & other){
		
		this->MPType = other.MPType;
		
		this->systems = other.systems;
		

		
		system_name_mid = other.system_name_mid;
		system_name_top = other.system_name_top;
		system_name_bottom = other.system_name_bottom;
		

		
		set_mp(v_target,other.v_target);
		set_mp(u_target,other.u_target);
		set_mp(crit_val_left,other.crit_val_left);
		set_mp(crit_val_right,other.crit_val_right);
		

        
        for (int ii=0; ii<other.num_projections; ii++) {
            add_projection(other.pi[ii]);
        }
        
	}
	
	
	
	
	void init();
	
	void clear()
	{
		
		clear_mp(v_target);
		clear_mp(u_target);
		clear_mp(crit_val_left);
		clear_mp(crit_val_right);
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(pi[ii]);
			}
			free(pi);
		}
		
	}
	
    
	
    
	
};





/**
 \brief the mp version of the derived solver type for the midtrack solver.
 
 */

// this must be defined before the double version, because double has mp.
class midpoint_eval_data_mp : public SolverMultiplePrecision
{
public:
	
	int num_mid_vars; ///< the number of variables for the middle of the face.
	int num_bottom_vars; ///< the number of variables for the bottom edge.
	int num_top_vars; ///< the number of variables for the top edge.
	
    
	StraightLineProgramGlobalPointers top_memory; ///< the memory for the top system
	StraightLineProgramGlobalPointers mid_memory; ///< the memory for the middle system
	StraightLineProgramGlobalPointers bottom_memory; ///< the memory for the bottom system
	
	// these are all merely pointers, and should only be assigned to memory set by the SLP creation routine inside of setupProg(), by the MidpointConfiguration::setup() call
	
	prog_t *SLP_top; ///< pointer to SLP for top edge.
	prog_t *SLP_bottom; ///< pointer to SLP for bottom edge.
	prog_t *SLP_mid; ///< pointer to SLP for middle of face
		// these are to be freed by the midpoint_setup object.
	
	
	std::shared_ptr<SystemRandomizer> randomizer_bottom; ///< pointer to randomizer for bottom edge.
	std::shared_ptr<SystemRandomizer> randomizer_top; ///< pointer to randomizer for top edge.
	
	vec_mp *pi; ///< projection vectors
	int num_projections; ///< projection vectors.  this should be 2 when filled.
	
	//patch already lives in the base class.
	
	comp_mp v_target; ///< the target value for v.
	comp_mp u_target; ///< the target value for u.
	
    
	
	comp_mp crit_val_left; ///< the projection value of the left critslice.
	comp_mp crit_val_right; ///< the projection value of the right critslice.
	
	
	comp_mp u_start; ///< initial value for u
	comp_mp v_start;///< initial value for v
	
	comp_mp one; /// < the number 1
	comp_mp zero; ///< the number 0.
	comp_mp half;///< the number 1/2
	
	
	vec_mp *pi_full_prec; ///< full precision version of pi
	comp_mp v_target_full_prec;///< full precision version of v_target
	comp_mp u_target_full_prec;///< full precision version of u_target
	
	
	comp_mp crit_val_left_full_prec;///< full precision version of crit_val_left
	comp_mp crit_val_right_full_prec;///< full precision version of crit_val_right
	
	comp_mp half_full_prec;///< full precision version of half
	comp_mp u_start_full_prec;///< full precision version of u_start
	comp_mp v_start_full_prec;///< full precision version of v_start
	comp_mp one_full_prec;///< full precision version of one
	comp_mp zero_full_prec;///< full precision version of zero
	
	
	
	
	// default initializer
	midpoint_eval_data_mp() : SolverMultiplePrecision(){
		
		reset_counters();
		
		init();
	}
	
	midpoint_eval_data_mp(int mp) : SolverMultiplePrecision(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~midpoint_eval_data_mp(){
		
		clear();
		// no need to reset the counters.
	}
	
	
	
	/**
	 \brief print some relevant data to the screen
	 */
	void print()
	{
        
	};
	
	/**
	 \brief reset the number of projections
	 
	 */
	void reset_counters()
	{
		this->num_projections = 0;
	} // re: reset_counters
	
	
	midpoint_eval_data_mp & operator=(const midpoint_eval_data_mp & other)
	{
		copy(other);
		return *this;
	}
	
	midpoint_eval_data_mp(const midpoint_eval_data_mp & other)
	{
		init();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	
	/**
	 \brief a non-usable broadcast send method.  midpoint_solver is intended for use in what appears to be a serial mode, even though it is used by many processors.  
	 
	 \see Surface::make_face
	 
	 \todo convert the name to bcast_send
	 
	 \return SUCCESSFUL
	 \param mpi_config the current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief a non-usable broadcast receive method.  midpoint_solver is intended for use in what appears to be a serial mode, even though it is used by many processors.
	 
	 \see Surface::make_face
	 
	 \todo convert the name to bcast_receive
	 
	 \return SUCCESSFUL
	 \param mpi_config the current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief main setup call, taking in MidpointConfiguration and a witness set, and producing a ready-to-go solver.
	 
	 \return SUCCESSFUL
	 \param md_config Configuration for the midpoint tracker
	 \param W input witness set
	 \param solve_options the current state of the solver
	 */
	int setup(MidpointConfiguration & md_config,
              const WitnessSet & W,
              SolverConfiguration & solve_options);
	
	
	
	/**
	 \brief add a projection to the solver.
	 \param proj the new projection to add.
	 */
	void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else {
			this->pi = (vec_mp *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_mp));
		}
		
		init_vec_mp(pi[num_projections],0);
		vec_cp_mp(pi[num_projections], proj);
		
		if (this->MPType==2){
			if (this->num_projections==0) {
				this->pi_full_prec = (vec_mp *) br_malloc(sizeof(vec_mp));
			}
			else {
				this->pi_full_prec = (vec_mp *) br_realloc(pi_full_prec,(this->num_projections+1) *sizeof(vec_mp));
			}
			
			init_vec_mp2(pi_full_prec[num_projections],0,1024);
			vec_cp_mp(pi_full_prec[num_projections], proj);
		}
		num_projections++;
		
	}
	
	
	
	
	/**
	 \brief reset to empty state.
	 */
	void clear()
	{
        
		if (this->num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				clear_vec_mp(pi[ii]);
			}
			free(pi);
		}
		
		clear_mp(v_target);
		clear_mp(u_target);

		clear_mp(crit_val_left);
		clear_mp(crit_val_right);
		
		clear_mp(u_start);
		clear_mp(v_start);
		
		clear_mp(half);
		clear_mp(one);
		clear_mp(zero);
		
		
		if (this->MPType==2){
            
			if (this->num_projections>0) {
				for (int ii=0; ii<num_projections; ii++) {
					clear_vec_mp(pi_full_prec[ii]);
				}
				free(pi_full_prec);
			}
			
			clear_mp(v_target_full_prec);
			clear_mp(u_target_full_prec);

			clear_mp(crit_val_left_full_prec);
			clear_mp(crit_val_right_full_prec);
			
			clear_mp(u_start_full_prec);
			clear_mp(v_start_full_prec);
			
			clear_mp(half_full_prec);
			clear_mp(one_full_prec);
			clear_mp(zero_full_prec);
			
		}
        
		
		
	} // re: clear
	
	
	/**
	 \brief initialization of the solver.  declares some function pointers, etc.
	 */
	void init(); // in the .cpp file.  must be there, not here.
    
	
	void copy(const midpoint_eval_data_mp & other)
	{
		SolverMultiplePrecision::copy(other);
		
		this->num_mid_vars = other.num_mid_vars;
		this->num_top_vars = other.num_top_vars;
		this->num_bottom_vars = other.num_bottom_vars;
		
		this->bottom_memory = other.bottom_memory;
		this->top_memory = other.top_memory;
		this->mid_memory = other.mid_memory;
		
		this->SLP_mid = other.SLP_mid;
		this->SLP_top = other.SLP_top;
		this->SLP_bottom = other.SLP_bottom;
		
		for (int ii=0; ii<other.num_projections; ii++) {
			if (other.MPType==2) {
				add_projection(other.pi_full_prec[ii]);
			}
			else{
				add_projection(other.pi[ii]);
			}
		}
		
		//patch already lives in the base class.
		
		set_mp(v_target,other.v_target);
		set_mp(u_target,other.u_target);
		
		randomizer_bottom = other.randomizer_bottom;
		randomizer_top = other.randomizer_top;
		
		
		set_mp(crit_val_left,other.crit_val_left);
		set_mp(crit_val_right,other.crit_val_right);
		
		set_mp(u_start, other.u_start);
		set_mp(v_start, other.v_start);
		
		
		if (this->MPType==2) {
			set_mp(v_target_full_prec,other.v_target_full_prec);
			set_mp(u_target_full_prec,other.u_target_full_prec);
			
			
			set_mp(crit_val_left_full_prec,other.crit_val_left_full_prec);
			set_mp(crit_val_right_full_prec,other.crit_val_right_full_prec);
			
			set_mp(u_start_full_prec, other.u_start_full_prec);
			set_mp(v_start_full_prec, other.v_start_full_prec);
		}
        
	} // re: copy
	
	
	
	
}; // re: class nullspace_eval_data_mp










// the double version
// this must be defined after the mp version, because double has mp.
class midpoint_eval_data_d : public SolverDoublePrecision
{
public:
	
	midpoint_eval_data_mp *BED_mp; ///< pointer to an MP type midpoint eval data. used only for AMP
	
	int num_mid_vars; ///< the number of variables in the mid system
	int num_top_vars;///< the number of variables in the top system
	int num_bottom_vars;///< the number of variables in the bottom system
	
	StraightLineProgramGlobalPointers bottom_memory;///< the memory for the bottom system
	StraightLineProgramGlobalPointers top_memory;///< the memory for the top system
	StraightLineProgramGlobalPointers mid_memory;///< the memory for the mid system
	
	prog_t *SLP_bottom; ///< a pointer to the SLP for the bottom system
	prog_t *SLP_top;///< a pointer to the SLP for the top system
	prog_t *SLP_mid;///< a pointer to the SLP for the mid system
	
	vec_d *pi; ///< the projection being used.  should be exactly 2 when fully populated.
	int num_projections; ///< the number of projections.
	
	std::shared_ptr<SystemRandomizer> randomizer_bottom; ///< a pointer to the randomizer for the bottom system
	std::shared_ptr<SystemRandomizer> randomizer_top;///< a pointer to the randomizer for the top system
	
	
	//patch already lives in the base class.
	
	comp_d v_target; ///< the target value for u, should be between 0 and 1.
	comp_d u_target;///< the target value for v, should be between 0 and 1.
	
	
	
	comp_d crit_val_left; ///< the projection value of the left critslice.
	comp_d crit_val_right; ///< the projection value of the right critslice.
	
	comp_d half; ///< the number 1/2
	comp_d u_start; ///< the initial value of u, probably 1/2
	comp_d v_start;///< the initial value of v, probably 1/2
	comp_d one; ///< the number 1
	comp_d zero;///< the number 0.
	
	
	// default initializer
	midpoint_eval_data_d() : SolverDoublePrecision(){
		
		reset_counters();
		
		init();
	}
	
	midpoint_eval_data_d(int mp) : SolverDoublePrecision(mp){
		
		this->MPType = mp;
		reset_counters();
		
		init();
	}
	
	virtual ~midpoint_eval_data_d(){
		
		midpoint_eval_data_d::clear();
		// no need to reset the counters.
	}
	
	/**
	 \brief print some relevant data to the screen
	 
	 */
	void print()
	{
		std::cout << "num_variables " << num_variables << std::endl;
		std::cout << "num_mid_vars " << num_mid_vars << std::endl;
		std::cout << "SLP_top "      << SLP_top->size << std::endl;
		std::cout << "SLP_bottom " << SLP_bottom->size << std::endl;
		std::cout << "SLP_mid " << SLP_mid->size << std::endl;
        
		print_comp_matlab(u_target,"u_target");
		print_comp_matlab(v_target,"v_target");
		
		print_comp_matlab(crit_val_left,"crit_val_left");
		print_comp_matlab(crit_val_right,"crit_val_right");
		
		
		std::cout << num_projections << " projections: " << std::endl;
		for (int ii=0; ii<2; ii++) {
			print_point_to_screen_matlab(this->pi[ii],"pi");
		}
		
	};
	
	
	
	void reset_counters()
	{
		num_projections=0;
	} // re: reset_counters
	
	
	midpoint_eval_data_d & operator=(const midpoint_eval_data_d & other)
	{
		this->MPType = other.MPType;
		copy(other);
		return *this;
	}
	
	midpoint_eval_data_d(const midpoint_eval_data_d & other)
	{
		init();
		SolverDoublePrecision();
		midpoint_eval_data_d();
		
		//no need to clear or reset counters, as this is a new object.
		copy(other);
	}
	
	
	// MPI SENDS AND RECEIVES
	/**
	 \brief a non-usable broadcast send.  this type is not intended to be sent about. 
	 
	 \see surface_decompositon::make_face
	 
	 \return SUCCESSFUL
	 \param mpi_config the current state of MPI
	 */
	int send(ParallelismConfig & mpi_config);
	
	/**
	 \brief a non-usable broadcast receive.  this type is not intended to be sent about.
	 
	 \see surface_decompositon::make_face
	 
	 \return SUCCESSFUL
	 \param mpi_config the current state of MPI
	 */
	int receive(ParallelismConfig & mpi_config);
	
	
	/**
	 \brief main setup call, taking in MidpointConfiguration and a witness set, and producing a ready-to-go solver.
	 
	 \return SUCCESSFUL
	 \param md_config Configuration for the midpoint tracker
	 \param W input witness set
	 \param solve_options the current state of the solver
	 */
	int setup(MidpointConfiguration & md_config,
              const WitnessSet & W,
              SolverConfiguration & solve_options);
	
	/**
	 \brief add a projection to the solver.
	 \param proj the new projection to add.
	 */
	void add_projection(vec_mp proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_d *) br_malloc(sizeof(vec_d));
		}
		else {
			this->pi = (vec_d *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_d));
		}
		
		init_vec_d(pi[num_projections],0);
		vec_mp_to_d(pi[num_projections], proj);
		
		num_projections++;
	}
	
	/**
	 \brief add a projection to the solver.
	 \param proj the new projection to add.
	 */
	void add_projection(vec_d proj)
	{
		if (this->num_projections==0) {
			this->pi = (vec_d *) br_malloc(sizeof(vec_d));
		}
		else {
			this->pi = (vec_d *) br_realloc(pi,(this->num_projections+1) *sizeof(vec_d));
		}
		
		init_vec_d(pi[num_projections],0);
		vec_cp_d(pi[num_projections], proj);
		
		num_projections++;
	}
private:
	
	
	void clear()
	{
        
		//yeah, there's some stuff to clear here.
		if (this->MPType==2) {
			delete this->BED_mp;
		}
		
		
		for (int ii=0; ii<num_projections; ii++) {
			clear_vec_d(pi[ii]);
		}
		free(pi);
		

	} // re: clear
	
	void init();
	
	void copy(const midpoint_eval_data_d & other)
	{
		
		
		SolverDoublePrecision::copy(other);
		
		this->num_mid_vars = other.num_mid_vars;
		this->num_top_vars = other.num_top_vars;
		this->num_bottom_vars = other.num_bottom_vars;
		
		this->top_memory = other.top_memory;
		this->bottom_memory = other.bottom_memory;
		this->mid_memory = other.mid_memory;
		
		this->SLP_bottom = other.SLP_bottom;
		this->SLP_mid = other.SLP_mid;
		this->SLP_top = other.SLP_top;
        
		
		for (int ii=0; ii<other.num_projections; ii++) {
			add_projection(other.pi[ii]);
		}
		
		//patch already lives in the base class.
		
		set_d(v_target,other.v_target);
		set_d(u_target,other.u_target);
		
		randomizer_bottom = other.randomizer_bottom;
		randomizer_top = other.randomizer_top;
		
		set_d(crit_val_left,other.crit_val_left);
		set_d(crit_val_right,other.crit_val_right);
		
		set_d(u_start, other.u_start);
		set_d(v_start, other.v_start);
	}
};







/** 
 \brief the main function for connecting a midpoint to the right or left, or at least finding the point to which it SHOULD connect.
 
 \return SUCCESSFUL
 \param W the input witness set with the start points and patches
 \param solve_out The output method, carrying solutions and metadata.
 \param md_config the configuration for the track about to happen
 \param solve_options the current state of the solver
 */

int midpoint_solver_master_entry_point(const WitnessSet & W, // carries with it the start points, and the linears.
                                       SolverOutput & solve_out, // new data goes in here
                                       MidpointConfiguration & md_config,
                                       SolverConfiguration		& solve_options);





/**
 \brief Evaluator function for the midpoint solver.
 
 \todo: explain with diagram how this works
 
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
int midpoint_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED);


/**
 \brief Evaluator function for the midpoint solver.
 
 \todo: explain with diagram how this works
 
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
int midpoint_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED);




/**
 \brief dehomogenization method for this solver.
 

 dehomogenizes only the mid variables, which come first.
 
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
int midpoint_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp);


/**
 \brief change the precision of an MP midpoint evaluator.
 
 This function fits the format for all Bertini precision changers.
 
 \return the number 0.
 \param ED pointer to the evaluator data to change.
 \param new_prec the precision to change to.
 */
int change_midpoint_eval_prec(void const *ED, int new_prec);


/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_midpoint_d(endgame_data_t *EG,
                            tracker_config_t *T,
                            void const *ED);

/**
 \brief check whether an input point is a solution, by residuals and ratio tolerances.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param EG the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_issoln_midpoint_mp(endgame_data_t *EG,
                             tracker_config_t *T,
                             void const *ED);




/**
 \brief check whether an input point is a valid start point, by residuals of the functions.
 
 \return a boolean indicating whether it is in fact a solution to the system.
 \param testpoint the solution to check.
 \param T the current tracker config, has tolerances.
 \param ED the midtrack eval data.
 */
int check_isstart_midpoint_d(point_d testpoint,
                             tracker_config_t *T,
                             void const *ED);

/**
 \brief a small function for verifying the homogeneousness of the evaluator (for correctness), etc.
 
 This function has no real purpose except for development.
 
 \param current_values the input
 \param ED the midtrack eval data.
 */
void check_midpoint_evaluator(point_mp current_values,
                              void const *ED);





/**
 \brief slave entry point, though it should not be used currently, as the solver is NOT parallel functional.
 
 The mid tracker is called in a loop, indeed, but is not called in parallel, because you only track one point for one set of parameter values.  Nontetheless, this function is here.
 
 \param solve_options The current state of the solver.
 */
void midpoint_slave_entry_point(SolverConfiguration & solve_options);










#endif


