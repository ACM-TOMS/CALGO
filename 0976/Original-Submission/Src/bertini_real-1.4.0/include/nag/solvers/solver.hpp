#ifndef SOLVER_MAIN_HEADER_H
#define SOLVER_MAIN_HEADER_H

/** \file solver.hpp */



#include "bertini1/bertini_headers.hpp"

#include "io/fileops.hpp"
#include "programConfiguration.hpp"
#include "nag/solvers/postProcessing.hpp"
#include "nag/system_randomizer.hpp"
#include "decompositions/decomposition.hpp"
#include "limbo.hpp"

///////////
//
//    SOLVER CONFIGURATION
//
//////////



extern _comp_d  **mem_d;
extern _comp_mp **mem_mp;
extern int *size_d;  // size of mem_d
extern int *size_mp;  // size of mem_mp
extern int *mem_needs_init_d; // determine if mem_d has been initialized
extern int *mem_needs_init_mp; // determine if mem_mp has been initialized






/**
 \brief A class that enables use of > 1 SLP for evaluation.
 
 Because Bertini uses some global pointers to perform evaluation of an SLP, some trickery had to be developed in order to hold multiple SLP's simultaneously.  This class is that trickery.  It is capable of setting these Bertini-globals to the appropriate locations in memory for a particular SLP, and then making them go NULL so that the memory is essentially protected.

 this class is ignorant of MPtype.
 
 */
class StraightLineProgramGlobalPointers
{
public:
	
	/**
	 \brief Find the values of the globals, and copy them to the local, internal pointers in this class.
	 */
	void capture_globals();
	
	/**
	 \brief  Set the global pointers to be those in this object.
	 */
	void set_globals_to_this();
	
	/**
	 \brief Override the global values to NULL.  it prevents accidental erasure.  This does NOT clear the memory, so beware memory leaks!!!
	 */
	void set_globals_null()
	{
		mem_d				= NULL;
		mem_mp				= NULL;
		size_d				= NULL;  // size of mem_d
		size_mp				= NULL;  // size of mem_mp
		mem_needs_init_d	= NULL; // determine if mem_d has been initialized
		mem_needs_init_mp	= NULL; // determine if mem_mp has been initialized
	}
	
	
	
	StraightLineProgramGlobalPointers()
	{
		init();
	}
	
	
	
	StraightLineProgramGlobalPointers & operator=( const StraightLineProgramGlobalPointers & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	StraightLineProgramGlobalPointers(const StraightLineProgramGlobalPointers & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	
	
	
	
protected:
	
	void init()
	{
		this->local_mem_d = NULL;
		this->local_mem_mp = NULL;
		this->local_size_d = NULL;  // size of mem_d
		this->local_size_mp = NULL;  // size of mem_mp
		this->local_mem_needs_init_d = NULL; // determine if mem_d has been initialized
		this->local_mem_needs_init_mp = NULL; // determine if mem_mp has been initialized
	}
	
	void copy(const StraightLineProgramGlobalPointers & other)
	{
		if (other.local_mem_d!=NULL)
			this->local_mem_d = other.local_mem_d;
		
		if (other.local_mem_mp!=NULL)
			this->local_mem_mp = other.local_mem_mp;
		
		if (other.local_size_d!=NULL)
			this->local_size_d = other.local_size_d;  // size of mem_d
		
		if (other.local_size_mp!=NULL)
			this->local_size_mp = other.local_size_mp;  // size of mem_mp
		
		if (other.local_mem_needs_init_d!=NULL)
			this->local_mem_needs_init_d = other.local_mem_needs_init_d; // determine if mem_d has been initialized
		
		if (other.local_mem_needs_init_mp!=NULL)
			this->local_mem_needs_init_mp = other.local_mem_needs_init_mp; // determine if mem_mp has been initialized
	}
	
	
private:
	
	_comp_d  **local_mem_d;
	_comp_mp **local_mem_mp;
	int *local_size_d;  // size of mem_d
	int *local_size_mp;  // size of mem_mp
	int *local_mem_needs_init_d; // determine if mem_d has been initialized
	int *local_mem_needs_init_mp; // determine if mem_mp has been initialized
};











/**
 \brief holds an SLP, SystemRandomizer, and StraightLineProgramGlobalPointers, for complete encapsulation of the system into one unit.
 
 \ingroup mpienabled
 
 complete with send/receive methods.
 
 */
class CompleteSystem
{
	
	boost::filesystem::path input_filename_; ///< the name of the input file generating this system
	StraightLineProgramGlobalPointers memory_; ///< the memory for the SLP
	prog_t * SLP_; ///< the actual SLP.  a pointer to it.  comes with have_SLP, which indicates whether this object owns the SLP.
	bool have_SLP_; ///< indicator of whether this object owns the SLP
	
	
	bool have_randomizer_; ///< indicator whether this object owns the randomizer
	std::shared_ptr<SystemRandomizer> randomizer_; ///< pointer to a randomizer.  comes with have_randomizer to indicate ownership.
	
	int num_variables_; ///< the number of variables in the system.
	
	
	
	int MPType; ///< the MP type.
	
public:
	
	
	/**
	 get a the pointer to the SLP for this system
	 \return A pointer to the SLP.
	 */
	prog_t * SLP_pointer()
	{
		return SLP_;
	}
	
	/**
	 get the memory for the system
	 \return the SLP memory object stored internally to this class.
	 */
	StraightLineProgramGlobalPointers memory()
	{
		return memory_;
	}
	
	/**
	 get the input file name
	 \return the name of the input file which was parsed into this CompleteSystem
	 */
	boost::filesystem::path input_filename() const
	{
		return input_filename_;
	}
	
	/**
	 set the input file name
	 \param new_name the name of the input file which was parsed into this CompleteSystem
	 */
	void set_input_filename(const boost::filesystem::path & new_name)
	{
		input_filename_ = new_name;
	}
	
	/**
	 \brief get the number of variables in the system
	 
	 \return the number of variables in the system
	 */
	int num_variables() const
	{
		return num_variables_;
	}
	
	int num_variables()
	{
		return num_variables_;
	}
	/**
	 \brief get a shared pointer to the randomizer
	 \return a shared pointer to the randomizer object stored in this CompleteSystem.
	 */
	std::shared_ptr<SystemRandomizer> randomizer()
	{
		return randomizer_;
	}
	
	
	/**
	 \brief send the complete system to everyone in the communicator
	 
	 
	 \param mpi_config the current state of MPI
	 */
	void bcast_send(ParallelismConfig & mpi_config)
	{
		
		
		int * buffer = new int[2];
		buffer[0] = MPType;
		buffer[1] = num_variables();
		MPI_Bcast(buffer, 2, MPI_INT, mpi_config.head(), mpi_config.comm());
		
		delete [] buffer;
		
		//need something here  to send randomizer to everyone.
		randomizer_->bcast_send(mpi_config);
		
		memory_.set_globals_to_this();
		bcast_prog_t(SLP_, this->MPType, 0, 0); // last two arguments are: myid, headnode
		memory_.set_globals_null();
	}
	
	
	/**
	 \brief receive from the head node
	 
	 
	 \param mpi_config the current state of MPI
	 */
	void bcast_receive(ParallelismConfig & mpi_config)
	{
		
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
		}
		else
		{
			SLP_ = new prog_t;
		}
		
		int * buffer = new int[2];
		
		MPI_Bcast(buffer, 2, MPI_INT, mpi_config.head(), mpi_config.comm());
		
		
		MPType = buffer[0];
		num_variables_ = buffer[1];
		
		delete [] buffer;
		
		
		// here, need something to receive from head.
		randomizer_ = std::make_shared<SystemRandomizer> (*(new SystemRandomizer));
		
		randomizer_->bcast_receive(mpi_config);
		
		
		
		bcast_prog_t(SLP_, MPType, 1, 0); // last two arguments are: myid, headnode
		initEvalProg(MPType);
		memory_.capture_globals();
		memory_.set_globals_null();
		
		have_SLP_ = true;
	}
	
	
	
	CompleteSystem()
	{
		init();
	}
	
	CompleteSystem(const Decomposition & D, tracker_config_t * T)
	{
		init();
		get_system(D,T);
	}
	
	
	~CompleteSystem()
	{
		clear();
	}
	
	
	CompleteSystem & operator=( const CompleteSystem & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	CompleteSystem(const CompleteSystem & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	
	void init()
	{
		
		have_randomizer_ = false;
		
		num_variables_ = 0;
		input_filename_ = "unset_input_filename_complete_system";
		MPType = -1;
		
		have_SLP_ = false;
	}
	
	
	void copy(const CompleteSystem & other)
	{
		
		this->num_variables_ = other.num_variables();
		
		this->MPType = other.MPType;
		
		set_input_filename(other.input_filename());
		
		if (this->have_SLP_) {
			this->memory_.set_globals_to_this();
			clearProg(this->SLP_, this->MPType, 1); // 1 means call freeprogeval()
			
			if (!other.have_SLP_) {
				delete this->SLP_;
				this->have_SLP_ = false;
			}
			
		}
		
		if (other.have_SLP_) {
			
			if (!this->have_SLP_)
			{
				this->SLP_ = new prog_t;
			}
			
			
			
			cp_prog_t(this->SLP_, other.SLP_);
			this->memory_.set_globals_null();
			initEvalProg(this->MPType);
			this->memory_.capture_globals();
			
			this->have_SLP_ = true;
		}
		
		
		randomizer_ = other.randomizer_; // copy the shared pointer value.
		
		
		
	}
	
	
	void clear()
	{

		
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
			delete SLP_;
			have_SLP_ = false;
		}
		
		
	}
	
	
	/**
	 \brief set up a system, based on a Decomposition.
	 
	 gets the randomizer as a pointer to a randomizer, based on the pointer in D.  sets up the SLP as local to this thing, and captures the memory to local as well.
	 
	 \param D the Decomposition for which to set up a complete system
	 \param T the current state of the tracker
	 */
	void get_system(const Decomposition & D, tracker_config_t * T)
	{
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
		}
		else{
			SLP_ = new prog_t;
		}
		
		
		this->randomizer_ = D.randomizer();
		have_randomizer_ = false;
		int blabla;  // i would like to move this.
		parse_input_file(D.input_filename(), &blabla);
		num_variables_ = D.num_variables();
		set_input_filename(D.input_filename());
		this->MPType = T->MPType;
		
		
		
		
		//	// setup a straight-line program, using the file(s) created by the parser
		int numVars = setupProg(SLP_, T->Precision, T->MPType);
		if (num_variables() != numVars) {
			std::cout << "numvars is incorrect...  setupprog gives " << numVars << ", but should be " << num_variables() << "." << std::endl;
			mypause();
			br_exit(-57189); // this should be a throw or something
		}
		
		preproc_data PPD;
		parse_preproc_data("preproc_data", &PPD);
		
		
		memory_.capture_globals();
		memory_.set_globals_null();
		
		
		
		preproc_data_clear(&PPD);
		
		have_SLP_ = true;
		
	}
	
	
	
	/**
	 \brief set up a system, from a filename and a tracker_config_t.
	 
	 doesn't do any randomization.  sets up the SLP as local to this thing, and captures the memory to local as well.
	 
	 \param input_filename the name of the file for which to set up a complete system
	 \param T the current state of the tracker
	 */
	void get_system(const boost::filesystem::path & input_filename, tracker_config_t * T)
	{
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
		}
		else{
			SLP_ = new prog_t;
		}
		
		
		int blabla;  // i would like to move this.
		parse_input_file(input_filename, &blabla);
		
		
		set_input_filename(input_filename);
		this->MPType = T->MPType;
		
		
		
		
		//	// setup a straight-line program, using the file(s) created by the parser
		num_variables_ = setupProg(SLP_, T->Precision, T->MPType);
		
		preproc_data PPD;
		parse_preproc_data("preproc_data", &PPD);
		
		
		memory_.capture_globals();
		memory_.set_globals_null();
		
		
		
		preproc_data_clear(&PPD);
		
		have_SLP_ = true;
		
	}
	
	/**
	 \brief set up a system, based on a Decomposition.
	 
	 gets the randomizer as a pointer to a randomizer, based on the pointer in D.  sets up the SLP as local to this thing, and captures the memory to local as well.
	 
	 \param new_input_name The name of the file to parse into this CompleteSystem.
	 \param randy A pointer to a randomizer to use.
	 \param T The current state of the tracker.
	 */
	void get_system(const boost::filesystem::path & new_input_name, std::shared_ptr<SystemRandomizer> randy, tracker_config_t * T)
	{
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
		}
		else{
			SLP_ = new prog_t;
		}
		
		
		
		this->randomizer_ = randy; // this could be eliminated via the use of smart pointers.  i'm doing that right now!!

		int blabla;  // i would like to move this.
		parse_input_file(new_input_name, &blabla);
		set_input_filename(new_input_name);
		this->MPType = T->MPType;
		
		
		
		
		//	// setup a straight-line program, using the file(s) created by the parser
		num_variables_ = setupProg(SLP_, T->Precision, T->MPType);
		have_SLP_ = true;
		
		preproc_data PPD;
		parse_preproc_data("preproc_data", &PPD);
		
		
		memory_.capture_globals();
		memory_.set_globals_null();
		
		
		
		preproc_data_clear(&PPD);
	}
	
	
	
	
	
	/**
	 \brief set up a system, based on a Decomposition.
	 
	 gets the randomizer as a pointer to a randomizer, based on the pointer in D.  sets up the SLP as local to this thing, and captures the memory to local as well.
	 
	 \param new_input_name The name of the file to parse into this CompleteSystem.
	 \param num_func_in The number of functions naturally appearing in the system.
	 \param num_func_out The number of functions for output.  If less than number in, will result in randomizer.  If greater, will result in logical inconsistency.
	 \param T The current state of the tracker (pointer to).
	 */
	void get_system(const boost::filesystem::path & new_input_name, int num_func_in, int num_func_out, tracker_config_t * T)
	{
		if (have_SLP_) {
			memory_.set_globals_to_this();
			clearProg(SLP_, this->MPType, 1); // 1 means call freeprogeval()
		}
		else{
			SLP_ = new prog_t;
		}
		

		
		
		
		
		int blabla;  // i would like to move this.
		parse_input_file(new_input_name, &blabla);
		set_input_filename(new_input_name);
		this->MPType = T->MPType;
		
		
		//	// setup a straight-line program, using the file(s) created by the parser
		num_variables_ = setupProg(SLP_, T->Precision, T->MPType);
		have_SLP_ = true;
		
		
		randomizer_ = std::make_shared<SystemRandomizer>(*(new SystemRandomizer));
		randomizer_->setup(num_func_out, num_func_in);
		
		
		preproc_data PPD;
		parse_preproc_data("preproc_data", &PPD);
		
		
		memory_.capture_globals();
		memory_.set_globals_null();
		
		
		
		preproc_data_clear(&PPD);
	}
	
	
};







/**
 \brief A flexible container for holding comp_mp, vec_mp, and mat_mp, allocated on the fly, but not cleared.  
 
 
 \see TemporariesDoublePrecision
 
 When an object of this type is stored over iterations of a function by being clever, you can spare many many init/clear calls from happening, which tend to be expensive for MP type
 */
class TemporariesMultiplePrecision
{
public:
	comp_mp * scalars; ///< the temporary comp_mp values.
	vec_mp * vectors;///< the temporary vec_mp values.
	mat_mp * matrices;///< the temporary mat_mp values.
	
	int num_scalars; ///< how many scalars have been allocated
	int num_vectors;///< how many vectors have been allocated
	int num_matrices;///< how many matrices have been allocated
	
	int curr_prec; ///< the current precision of all temps.  they must be all the same precision.
	
	
	
	
	
	TemporariesMultiplePrecision(){init();}
	~TemporariesMultiplePrecision(){clear();}
	
	TemporariesMultiplePrecision & operator=(const TemporariesMultiplePrecision & other){
		copy(other);
		
		return *this;
	}
	
	TemporariesMultiplePrecision(const TemporariesMultiplePrecision & other){
		init();
		copy(other);
	}
	
	
	
	
	/**
	 \brief purge this temp object of all scalars
	 */
	inline void clear_scalars()
	{
		if (num_scalars>0) {
			for (int ii=0; ii<num_scalars; ii++) {
				clear_mp(scalars[ii]);
			}
			free(scalars);
		}
		num_scalars=0;
	}
	
	/**
	 \brief purge this temp object of all matrices
	 */
	void clear_matrices()
	{
		if (num_matrices>0) {
			for (int ii=0; ii<num_matrices; ii++) {
				clear_mat_mp(matrices[ii]);
			}
			free(matrices);
		}
		num_matrices = 0;
	}
	
	/**
	 \brief purge this temp object of all vectors
	 */
	void clear_vectors()
	{
		if (num_vectors>0) {
			for (int ii=0; ii<num_vectors; ii++) {
				clear_vec_mp(vectors[ii]);
			}
			free(vectors);
		}
		num_vectors = 0;
	}
	
	
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_scalars(int num_to_require)
	{
		if (num_scalars<num_to_require) {
			for (int ii=num_scalars; ii<num_to_require; ii++) {
				add_scalar();
			}
		}
	}
	
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_vectors(int num_to_require)
	{
		if (num_vectors<num_to_require) {
			for (int ii=num_vectors; ii<num_to_require; ii++) {
				add_vector();
			}
		}
	}
	
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_matrices(int num_to_require)
	{
		if (num_matrices<num_to_require) {
			for (int ii=num_matrices; ii<num_to_require; ii++) {
				add_matrix();
			}
		}
	}
	
	
	
	/**
	 \brief add a scalar to the lineup, regardless of how many there are.
	 
	 \return the index of the added scalar
	 */
	int add_scalar()
	{
		int new_index = num_scalars;
		
		if (num_scalars==0) {
			scalars = (comp_mp *) br_malloc(sizeof(comp_mp));
		}
		else{
			scalars = (comp_mp *) br_realloc(scalars,(num_scalars+1)*sizeof(comp_mp));
		}
		init_mp(scalars[new_index]);
		
		num_scalars++;
		return new_index;
	}
	
	/**
	 \brief add a vector to the lineup, regardless of how many there are.
	 
	 \return the index of the added vector
	 */
	int add_vector()
	{
		int new_index = num_vectors;
		if (num_vectors==0) {
			vectors = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else {
			vectors = (vec_mp *) br_realloc(vectors,(num_vectors+1)*sizeof(vec_mp));
		}
		init_vec_mp(vectors[new_index],1);
		vectors[new_index]->size = 1;
		
		num_vectors++;
		return new_index;
		
		
	}
	
	/**
	 \brief add a matrix to the lineup, regardless of how many there are.
	 
	 \return the index of the added matrix
	 */
	int add_matrix()
	{
		int new_index = num_matrices;
		if (num_matrices==0) {
			matrices = (mat_mp *) br_malloc(sizeof(mat_mp));
		}
		else {
			matrices = (mat_mp *) br_realloc(matrices,(num_matrices+1)*sizeof(mat_mp));
		}
		init_mat_mp(matrices[new_index],1,1);
		matrices[new_index]->rows = matrices[new_index]->cols = 1;
		
		num_matrices++;
		return new_index;
	}
	
	
	/**
	 \brief change the precision of all member scalars, vectors, and matrices.
	 
	 \param new_prec The precision to change to.
	 */
	void change_prec(int new_prec)
	{
		if (new_prec==curr_prec) {
			return;
		}
		
		for (int ii=0; ii<num_scalars; ii++) {
			change_prec_mp(scalars[ii],new_prec);
		}
		for (int ii=0; ii<num_vectors; ii++) {
			change_prec_vec_mp(vectors[ii],new_prec);
		}
		for (int ii=0; ii<num_matrices; ii++) {
			change_prec_mat_mp(matrices[ii],new_prec);
		}
		
		curr_prec = new_prec;
	}
	
	
	
private:
	
	void init()
	{
		scalars = NULL;
		vectors = NULL;
		matrices = NULL;
		
		num_scalars = 0;
		num_vectors = 0;
		num_matrices = 0;
		
		curr_prec = mpfr_get_default_prec();
	}
	
	void copy(const TemporariesMultiplePrecision &other)
	{
		for (int ii=0; ii<other.num_scalars; ii++) {
			int new_index = this->add_scalar();
			set_mp(this->scalars[new_index],other.scalars[ii]);
		}
		
		for (int ii=0; ii<other.num_vectors; ii++) {
			int new_index = this->add_vector();
			vec_cp_mp(this->vectors[new_index],other.vectors[ii]);
		}
		
		for (int ii=0; ii<other.num_matrices; ii++) {
			int new_index = this->add_matrix();
			mat_cp_mp(this->matrices[new_index],other.matrices[ii]);
		}
		
		
		this->change_prec(other.curr_prec);
	}
	
	
	void clear()
	{
		clear_matrices();
		clear_vectors();
		clear_scalars();
	}
};


/**
 \brief A flexible container for holding comp_d, vec_d, and mat_d, allocated on the fly, but not cleared.
 
 
 \see TemporariesMultiplePrecision
 
 When an object of this type is stored over iterations of a function by being clever, you can spare many many init/clear calls from happening, which tend to be expensive for MP type, and not so bad for doubles.  nonetheless, here it is.
 */
class TemporariesDoublePrecision
{
public:
	comp_d * scalars;///< the temporary comp_mp values.
	vec_d * vectors;///< the temporary vec_mp values.
	mat_d * matrices;///< the temporary mat_mp values.
	
	int num_scalars;///< how many scalars have been allocated
	int num_vectors;///< how many vectors have been allocated
	int num_matrices;///< how many matrices have been allocated
	
	
	TemporariesDoublePrecision(){init();}
	~TemporariesDoublePrecision(){clear();}
	
	TemporariesDoublePrecision & operator=(const TemporariesDoublePrecision & other){
		copy(other);
		
		return *this;
	}
	
	TemporariesDoublePrecision(const TemporariesDoublePrecision & other){
		init();
		copy(other);
	}
	
	
	
	
	/**
	 \brief purge this temp object of all scalars
	 */
	inline void clear_scalars()
	{
		if (num_scalars>0) {
			for (int ii=0; ii<num_scalars; ii++) {
				clear_d(scalars[ii]);
			}
			free(scalars);
		}
		num_scalars=0;
	}
	
	/**
	 \brief purge this temp object of all matrices
	 */
	void clear_matrices()
	{
		if (num_matrices>0) {
			for (int ii=0; ii<num_matrices; ii++) {
				clear_mat_d(matrices[ii]);
			}
			free(matrices);
		}
		num_matrices = 0;
	}
	
	/**
	 \brief purge this temp object of all vectors
	 */
	void clear_vectors()
	{
		if (num_vectors>0) {
			for (int ii=0; ii<num_vectors; ii++) {
				clear_vec_d(vectors[ii]);
			}
			free(vectors);
		}
		num_vectors = 0;
	}
	
	
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_scalars(int num_to_require)
	{
		if (num_scalars<num_to_require) {
			for (int ii=num_scalars; ii<num_to_require; ii++) {
				add_scalar();
			}
		}
	}
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_vectors(int num_to_require)
	{
		if (num_vectors<num_to_require) {
			for (int ii=num_vectors; ii<num_to_require; ii++) {
				add_vector();
			}
		}
	}
	/**
	 \brief ensure that this temp object has AT LEAST as many as the input requirement.  Does not ensure equality of quantity.
	 
	 \param num_to_require The number the user wants the temp object to have.
	 */
	void ensure_have_matrices(int num_to_require)
	{
		if (num_matrices<num_to_require) {
			for (int ii=num_matrices; ii<num_to_require; ii++) {
				add_matrix();
			}
		}
	}
	
	
	
	/**
	 \brief add a scalar to the lineup, regardless of how many there are.
	 
	 \return the index of the added scalar
	 */
	int add_scalar()
	{
		int new_index = num_scalars;
		
		if (num_scalars==0) {
			scalars = (comp_d *) br_malloc(sizeof(comp_d));
		}
		else{
			scalars = (comp_d *) br_realloc(scalars,(num_scalars+1)*sizeof(comp_d));
		}
		
		num_scalars++;
		return new_index;
	}
	
	/**
	 \brief add a vector to the lineup, regardless of how many there are.
	 
	 \return the index of the added vector
	 */
	int add_vector()
	{
		int new_index = num_vectors;
		if (num_vectors==0) {
			vectors = (vec_d *) br_malloc(sizeof(vec_d));
		}
		else {
			vectors = (vec_d *) br_realloc(vectors,(num_vectors+1)*sizeof(vec_d));
		}
		init_vec_d(vectors[new_index],1);
		vectors[new_index]->size = 1;
		
		num_vectors++;
		return new_index;
		
		
	}
	
	/**
	 \brief add a matrix to the lineup, regardless of how many there are.
	 
	 \return the index of the added matrix
	 */
	int add_matrix()
	{
		int new_index = num_matrices;
		if (num_matrices==0) {
			matrices = (mat_d *) br_malloc(sizeof(mat_d));
		}
		else {
			matrices = (mat_d *) br_realloc(matrices,(num_matrices+1)*sizeof(mat_d));
		}
		init_mat_d(matrices[new_index],1,1);
		matrices[new_index]->rows = matrices[new_index]->cols = 1;
		
		num_matrices++;
		return new_index;
	}
	
	
	
	
	
private:
	
	void init()
	{
		scalars = NULL;
		vectors = NULL;
		matrices = NULL;
		
		num_scalars = 0;
		num_vectors = 0;
		num_matrices = 0;
	}
	
	void copy(const TemporariesDoublePrecision &other)
	{
		for (int ii=0; ii<other.num_scalars; ii++) {
			int new_index = this->add_scalar();
			set_d(this->scalars[new_index],other.scalars[ii]);
		}
		
		for (int ii=0; ii<other.num_vectors; ii++) {
			int new_index = this->add_vector();
			vec_cp_d(this->vectors[new_index],other.vectors[ii]);
		}
		
		for (int ii=0; ii<other.num_matrices; ii++) {
			int new_index = this->add_matrix();
			mat_cp_d(this->matrices[new_index],other.matrices[ii]);
		}
		
	}
	
	
	void clear()
	{
		clear_matrices();
		clear_vectors();
		clear_scalars();
	}
};






/**
 
 \brief Holds the current state of the solver, including the all-important tracker_config_t struct.
 
 This class offers setups for the solver, which are system-independent, including methods for storing a configuration for restoration later.
 
 */
class SolverConfiguration : public ParallelismConfig
{
	
	int verbose_level_; ///< controls how much info is printed to screen
	
	
	std::map< std::string, std::pair<bool, tracker_config_t> > tracker_config_backups_;
	
public:
	
	
	
	bool robust; ///< whether to use robust mode
	tracker_config_t T; ///< the ubiquitous Bertini tracker configuration

	
	preproc_data PPD; ///< the structure of the current variable groups
	

	
	int path_number_modulus; ///< for display of path number to screen
	
	
	
	int use_midpoint_checker; ///< whether to use the midpoint checker.  use of the checker is currently broken in bertini_real.
	double midpoint_tol; ///< how far apart midpoints must be to be considered distinct.
	
	int use_gamma_trick;///< whether to use the gamma trick for start systems.
	
	
	
	/**
	 \brief get the level of verbosity
	 
	 \return the level of verbosity
	 */
	inline int verbose_level() const
	{
		return verbose_level_;
	}
	
	/**
	 \brief set the level of verbosity
	 
	 \param new_level the new level of verbosity
	 \return The verbosity level, which you just put in.
	 */
	int verbose_level(int new_level)
	{
		return verbose_level_ = new_level;
	}
	
	
	
	
	
	/**
	 \brief copy the stored tracker config to the active one.
	 \param lookup_key A string by which to reference the tracker config.  Tracker configs must be restored from before they are backed up to again.
	 */
	void restore_tracker_config(std::string lookup_key)
	{
		if (tracker_config_backups_.find(lookup_key)==tracker_config_backups_.end()){
			throw std::logic_error("trying to restore from a backup tracker config which does not exist: " + lookup_key);
		}
		else{
			tracker_config_clear(&T);
			cp_tracker_config_t(&T, &tracker_config_backups_[lookup_key].second);
			tracker_config_backups_[lookup_key].first = true;
		}
	}
	
	/**
	 \brief copy the active tracker config to the backup one.
	 
	 the tracker config is backed up, and locked for overwriting.  It *must* be restored before it is backed up to again.
	 
	 \param lookup_key A string by which to refer to the tracker config.
	 */
	void backup_tracker_config(std::string lookup_key)
	{
		
		if (tracker_config_backups_.find(lookup_key)==tracker_config_backups_.end()) {
			tracker_config_backups_[lookup_key].first = false;
			mpf_init2(tracker_config_backups_[lookup_key].second.latest_newton_residual_mp, T.AMP_max_prec);
			cp_tracker_config_t(&(tracker_config_backups_[lookup_key].second), &T);
		}
		else{ // already have one of this name.  check and see if available.
			if (tracker_config_backups_[lookup_key].first == false) {
				std::stringstream ss;
				ss << "trying to back up to currently locked tracker backup: " << lookup_key;
				throw std::logic_error(ss.str());
			}
			else{
				tracker_config_clear(&(tracker_config_backups_[lookup_key].second));
				cp_tracker_config_t(&(tracker_config_backups_[lookup_key].second), &T);
				tracker_config_backups_[lookup_key].first = false;
			}
		}
	}
	
	
	/**
	 \brief tick the number of paths tracked by this process, under this solver_config.  if a multiple of 500, print to screen.
	 
	 \return The total number of paths tracked by this process so far.
	 */
	int increment_num_paths_tracked()
	{
		total_num_paths_tracked++;
		if (verbose_level()>=0 && (total_num_paths_tracked%500)==0 && ParallelismConfig::is_head()) {
			std::cout << "\t\t\t\t\ttracked " << total_num_paths_tracked << " paths total." << std::endl;
		}
		return total_num_paths_tracked;
	}
	
	SolverConfiguration(){
		init();
	}
	~SolverConfiguration(){
		tracker_config_clear(&this->T);
		preproc_data_clear(&this->PPD);
		
		for (auto iter = tracker_config_backups_.begin(); iter!=tracker_config_backups_.end(); iter++) {
			tracker_config_clear( &iter->second.second);
		}
	}
	 
	
	SolverConfiguration & operator=( const SolverConfiguration & other)
	{
		init();
		copy(other);
		return *this;
	}
	
	SolverConfiguration(const SolverConfiguration & other)
	{
		init();
		copy(other);
	} // re: copy
	
	
	void copy(const SolverConfiguration & other)
	{
		cp_tracker_config_t(&this->T, &other.T);
		
		this->total_num_paths_tracked = other.total_num_paths_tracked;
		this->robust = other.robust;
		
		cp_preproc_data(&(this->PPD), other.PPD);
		

		
		this->path_number_modulus = other.path_number_modulus;
		
		this->verbose_level(other.verbose_level());
		
		this->use_midpoint_checker = other.use_midpoint_checker;
		this->midpoint_tol = other.midpoint_tol;
		
		this->use_gamma_trick = other.use_gamma_trick;
		

	}
	
	
	
	void init();
	
	/**
	 \brief read in the preproc_data from file "preproc_data", and store in this object's PPD field.
	 */
	void get_PPD()
	{
		parse_preproc_data("preproc_data", &this->PPD);
	}
	
	
	
	
private:
	long long total_num_paths_tracked;
	
};






class SolverOutput; //  that'd be a forward declaration, jim.



/**
 \brief Metadata for solutions produced by the tracker, and stored in SolverOutput.
 */
class SolutionMetadata
{

friend SolverOutput;
	
	std::vector<long long> input_index; ///< the indices of the start points which ended here.
	long long output_index; ///< the output index of the endpoint solution
	int multiplicity; ///< how many paths ended here.
	bool is_finite; ///< flag for whether has been declared finite.
	bool is_singular;///< flag for whether has been declared singular.
	bool is_successful;///< flag for whether tracker was successful, or gave up for some reason -- wish this stored the reason too.
	bool is_real;///< flag for whether has been declared real.
	int cycle_number_ = -1;  ///< the cycle number of the path taken to the solution.  This number can regularize the path, if used correctly.
	
public:
	
	void SetCycleNumber(int c)
	{
		cycle_number_ = c;
	}

	int CycleNumber() const
	{
		return cycle_number_;
	}
	
	/**
	 \brief set the finiteness state
	 \param state set is_finite to the state.
	 */
	void set_finite(bool state){
		is_finite = state;
	}
	
	/**
	  \brief set the singularness state
	  \param state set is_singular to the state.
	  */
	void set_singular(bool state){
		is_singular = state;
	}
	
	/**
	 \brief set the successfulness state
	 \param state set is_successful to the state.
	 */
	void set_successful(bool state){
		is_successful = state;
	}
	
	/**
	 \brief set the multiplicity
	 \param state set multiplicity to the state.
	 */
	void set_multiplicity(int state){
		multiplicity = state;
	}
	
	/**
	 \brief set the realness state
	 \param state set is_real to the state.
	 */
	void set_real(int state) {
		is_real = state;
	}
	
	/**
	 \brief add another input start index which tracked to this solution
	 \param new_ind the index to add.
	 */
	void add_input_index(long long new_ind){
		input_index.push_back(new_ind);
	}
	
	/**
	 \brief set the output index for the found solution.  multiplicity>1 solutions are stored only once...
	 \param new_ind the index to assert
	 */
	void set_output_index(long long new_ind){
		output_index = new_ind;
	}
	

	friend std::ostream & operator<<(std::ostream &os, const SolutionMetadata & t)
	{
		
		
		os << t.output_index << std::endl;
		for (auto iter=t.input_index.begin(); iter!=t.input_index.end(); ++iter) {
			os << *iter << ' ';
		}
		os << '\n';
		
		os << t.multiplicity << ' ' << t.is_finite << ' ' << t.is_singular << ' ' << t.is_successful << ' ' << t.CycleNumber();
		
		return os;
	}
};


/**
 \brief Class for turning the output from a solver into a more useable form, namely into witness sets ultimately.
 
 This class acts as a VertexSet, holding vertices and metadata.  The main call is post_process()
 */
class SolverOutput : public PatchHolder, public LinearHolder, public NameHolder, public VertexSet
{

	
	
private:
	
	std::vector< SolutionMetadata > metadata;
	
	std::vector< std::pair<long long, long long> > ordering; /// created in the post-processing.

	std::vector< int > occuring_multiplicities;
	
public:
	
	int num_variables; ///< How many total variables there are.
	int num_natural_vars; ///< how many main variables there are, including the homogenizing variable.
	
	
	int num_higher_multiplicities()
	{
		return occuring_multiplicities.size();
	}
	
	void reset()
	{
		VertexSet::reset();
		reset_names();
		reset_linears();
		reset_patches();
		
		ordering.resize(0);
		metadata.resize(0);
		occuring_multiplicities.resize(0);
		
		
	}
	
	
	
	/**
	 \brief Set the num_variables and num_natural_vars to be the number in this SolverOutput object.
	 
	 \param W_transfer The input mutable witness set.
	 */
	void set_witness_set_nvars(WitnessSet & W_transfer)
	{
		W_transfer.set_num_variables(this->num_variables);
		W_transfer.set_num_natural_variables(this->num_natural_vars);
	}
	
	
	/**
	 \brief Add a solution and metadata to this solver output object.
	 
	 \param temp_vert The input vertex to pass in.
	 \param meta The solution metadata, containing input indices, output indices, etc.
	 */
	void add_solution(const Vertex & temp_vert, const SolutionMetadata & meta)
	{

		add_vertex(temp_vert);
		metadata.push_back(meta);

	}
	
	/**
	 \brief Get the finite solutions, including those which are singular or multiple, and put them in a witness set.
	 
	 This 'full' version also gets the linears and patches.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_noninfinite_w_mult_full(WitnessSet & W_transfer);
	
	std::vector<int> get_cyclenums_noninfinite_w_mult();
	
	/**
	 \brief Get the finite solutions, including those which are singular or multiple, and put them in a witness set.
	 	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_noninfinite_w_mult(WitnessSet & W_transfer);
	
	/**
	 \brief Get the nonsingular, finite, multiplicity one solutions, and put them in a witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_nonsing_finite_multone(WitnessSet & W_transfer);
	
	/**
	 \brief Assemble a map of multiplicities and points, for those solutions with multiplicity>1.
	 
	 \see get_multpos_full
	 
	 \param W_transfer The input mutable witness sets to populate.
	 */
	void get_multpos(std::map<int, WitnessSet> & W_transfer);
		
	
	/**
	 \brief Assemble a map of multiplicities and points, for those solutions with multiplicity>1.
	 
	 This 'full' version gets the linears and patches as well.
	 
	 \see get_multpos
	 
	 \param W_transfer The input mutable witness sets to populate.
	 */
	void get_multpos_full(std::map<int, WitnessSet> & W_transfer);

	
	/**
	 \brief Get the singular points from a solution set, and put them into a witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_sing(WitnessSet & W_transfer);
	
	
	
	
	/**
	 \brief Get the singular and finite points from a solution set, and put them into a witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_sing_finite(WitnessSet & W_transfer);
	
	
	
	
	/**
	 Put the patches and linears from this into the input witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_patches_linears(WitnessSet & W_transfer)
	{
		get_patches(W_transfer);
		get_linears(W_transfer);
	}
	
	/**
	 Put the patches only from this into the input witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_patches(WitnessSet & W_transfer)
	{
		W_transfer.copy_patches(*this);
	}
	
	/**
	 Put the linears only from this into the input witness set.
	 
	 \param W_transfer The input mutable witness set to populate.
	 */
	void get_linears(WitnessSet & W_transfer)
	{
		W_transfer.copy_linears(*this);
	}
	
	
	
	/**
	 \brief Bertini_real's version of post-processing.  options are set via the SolverConfiguration.
	 
	 \param endPoints The input for the method, having been converted into this format previously.
	 \param num_pts_to_check The number of endPoints.
	 \param preProcData structure containing the variable groups.
	 \param T The current tracker configuration.
	 \param solve_options The current state of the solver.
	 */
	void post_process(post_process_t *endPoints, int num_pts_to_check,
					  preproc_data *preProcData, tracker_config_t *T,
					  const SolverConfiguration & solve_options);
	
	
};



/**
 \brief Base class from whence all solver classes are derived.
 */
class Solver
{
	
protected:
	
	std::shared_ptr<SystemRandomizer> randomizer_; ///< Pointer to a randomizer.
	
	
	int verbose_level_;  ///< how verbose to be
	
	
public:
	
	bool IsRandomized() const
	{
		return !randomizer_->is_square();
	}
	
	
	/**
	 \brief get the level of verbosity
	 
	 \return the level of verbosity
	 */
	inline int verbose_level() const
	{
		return verbose_level_;
	}
	
	/**
	 \brief set the level of verbosity
	 
	 \param new_level the new level of verbosity
	 \return The verbosity level, which you just put in.
	 */
	int verbose_level(int new_level)
	{
		return verbose_level_ = new_level;
	}
	
	
	
	
	
	
	/**
	 \brief get a shared pointer to the randomizer
	 
	 \return a shared pointer to the randomizer
	 */
	std::shared_ptr<SystemRandomizer> randomizer()
	{
		return randomizer_;
	}
	
	
	void set_randomizer(std::shared_ptr<SystemRandomizer> randy)
	{
		randomizer_ = randy;
	}
	
	boost::filesystem::path preproc_file;
	boost::filesystem::path function_file;
	
	int num_steps; ///< the number of evaluations made using this evaluator
	
	
	
	
	
    bool received_mpi; ///< Whether this Solver has received its contents via MPI
	
	int MPType; ///< the multiple precision type for solve
	preproc_data preProcData; ///< information related to the SLP for system
	StraightLineProgramGlobalPointers SLP_memory;
	prog_t *SLP; ///< the SLP
	bool have_SLP, have_PPD;
	
	int num_variables; ///< the grand total number of variables to be computed on in the problem at hand
	
	// these function handles are required for all solvers.
	int (*evaluator_function_d) (point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *); // function handle to evaluator to use
	int (*evaluator_function_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *);
	
	int (*precision_changer)(void const *ED, int new_prec);
	
	int (*dehomogenizer)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *);
	
	int (*is_solution_checker_d)(endgame_data_t *, tracker_config_t *, void const *);
	int (*is_solution_checker_mp)(endgame_data_t *, tracker_config_t *, void const *);
	
	Solver(){
		init();
	}
	
	
	Solver(int new_mp_type){
		init();
		this->MPType = new_mp_type;
	}
	
	virtual ~Solver()
	{
		Solver::clear();
	}
	
	virtual int send(ParallelismConfig & mpi_config);
	
	virtual int receive(ParallelismConfig & mpi_config);
	
	virtual void print()
	{};
	
	
	void setup(prog_t * _SLP, std::shared_ptr<SystemRandomizer> randy)
	{
		setupPreProcData(const_cast<char *>(preproc_file.c_str()), &this->preProcData);
		have_PPD = true;
		
		//TODO: if there is a way, remove the above setup call, as it refers to parsing the preprocdata file, which is dumb.
		
		
		this->SLP = _SLP; // assign the pointer
		this->have_SLP = true;
		
		randomizer_ = randy;
	}
	
	void setup()
	{
		setupPreProcData(const_cast<char *>(preproc_file.c_str()), &this->preProcData);
	}
	
protected:
	
	void clear()
	{
		if (have_PPD) {
			preproc_data_clear(&this->preProcData);
			have_PPD = false;
		}
		
	}
	
	void init()
	{

		
        received_mpi = false;
        
		
		this->preproc_file = "preproc_data";
		this->function_file = "func_input";
		
		SLP = NULL;
		have_SLP = false;
		have_PPD = false;
		num_variables = 0;
		
		this->num_steps = 0;
		this->verbose_level_ = 0;
		
		// initialize the function handles.
		is_solution_checker_d = NULL;
		is_solution_checker_mp = NULL;
		evaluator_function_d = NULL;
		evaluator_function_mp = NULL;
		precision_changer = NULL;
		dehomogenizer = NULL;
	}
	// these referenced functions will need to be programmed into the derived classes.
	
	
	void copy(const Solver & other){
		cp_preproc_data(&this->preProcData, other.preProcData);
		
		this->SLP = other.SLP;
		
		randomizer_ = other.randomizer_; // merely copy the pointer
		
		this->is_solution_checker_d = other.is_solution_checker_d;
		this->is_solution_checker_mp = other.is_solution_checker_mp;
		
		this->evaluator_function_mp = other.evaluator_function_mp;
		this->evaluator_function_d = other.evaluator_function_d;
		this->precision_changer = other.precision_changer;
		this->dehomogenizer = other.dehomogenizer;
		
		this->num_variables = other.num_variables;
		this->MPType = other.MPType;
		this->verbose_level_ = other.verbose_level_;
		this->num_steps = other.num_steps;
	}
};  // re: generic Solver BASE class









class SolverMultiplePrecision : public Solver
{
	
public:
	
	patch_eval_data_mp patch; ///< patch in x
	
	comp_mp gamma;    ///< randomizer
	mpq_t *gamma_rat; ///< randomizer
	
	mpfr_prec_t curr_prec;
	
	TemporariesMultiplePrecision temp_vars;
	
	SolverMultiplePrecision() : Solver(){
		this->MPType = 2;
		init();
	}; // re: default constructor
	
	
	SolverMultiplePrecision(int new_mp_type) : Solver(new_mp_type){
		this->MPType = new_mp_type;
		init();
	}; // re: default constructor
	
	
	void init(){
		init_mp(gamma);
		
		if (this->MPType == 2 ) {
			gamma_rat = (mpq_t *)br_malloc(2 * sizeof(mpq_t));
			init_rat(gamma_rat);
		}
		
		curr_prec = mpfr_get_default_prec();
	}
	
	
	
	virtual ~SolverMultiplePrecision(){
		SolverMultiplePrecision::clear();
	}// re: default destructor
	
	
	
	
	SolverMultiplePrecision & operator=( const SolverMultiplePrecision & other)
	{
		Solver::operator= (other);
		copy(other);
		return *this;
	}  // re: assigment
	
	
	SolverMultiplePrecision(const SolverMultiplePrecision & other) : Solver(other){
		copy(other);
	} // re: copy
	
	virtual int send(ParallelismConfig & mpi_config);
	
	
	virtual int receive(ParallelismConfig & mpi_config);
	
	virtual void print()
	{};
	
	
	void setup(prog_t * _SLP, std::shared_ptr<SystemRandomizer> randy)
	{
		Solver::setup(_SLP, randy);
	}
	
	void setup()
	{
		Solver::setup();
	}
protected:
	
	void copy(const SolverMultiplePrecision & other){
		cp_patch_mp(&this->patch, other.patch);
	
		set_mp(this->gamma, other.gamma);
		
		if (other.MPType==2) {
			set_rat(this->gamma_rat, other.gamma_rat);
		}
		
		this->curr_prec = other.curr_prec;
	}
	
	void clear();
};








class SolverDoublePrecision : public Solver
{
	
public:
	
	patch_eval_data_d patch; ///< patch in x
	
	comp_d gamma;    ///< randomizer
	
	SolverMultiplePrecision *BED_mp; // why even have this?
    
	TemporariesDoublePrecision temp_vars;
    
	SolverDoublePrecision() : Solver(){
		init();
	} // re: default constructor
	
	
	SolverDoublePrecision(int new_mp_type) : Solver(new_mp_type){
		init();
	} // re: default constructor
	
	
	
	virtual ~SolverDoublePrecision(){
		SolverDoublePrecision::clear();
	}// re: default destructor
	
	
	SolverDoublePrecision & operator=( const SolverDoublePrecision & other)
	{
		Solver::operator= (other);
		
		copy(other);
		return *this;
	}  // re: assigment
	
	
	SolverDoublePrecision(const SolverDoublePrecision & other) : Solver(other){
		copy(other);
	} // re: copy
	
	virtual int send(ParallelismConfig & mpi_config);
	
	
	virtual int receive(ParallelismConfig & mpi_config);
	
	virtual void print()
	{};
	
	
	void setup(prog_t * _SLP, std::shared_ptr<SystemRandomizer> randy)
	{
		Solver::setup(_SLP, randy);
	}
	
	void setup()
	{
		Solver::setup();
		
		if (MPType==2) {
			
		}
	}
protected:
	
	
	
	
	
	void init()
	{
		Solver::init();
		
		init_d(gamma);
		
		if (MPType == 0) {
		}
	}
	
	void copy(const SolverDoublePrecision & other){
		
		cp_patch_d(&this->patch, other.patch);
		set_d(this->gamma, other.gamma);
	}
	
	void clear();
	
};








/**
 \brief just before calling a Solver, call this to adjust the tracker_config.
 
 This function ensures that the correct things inside T are intialized to the correct values, including at least one MP variable, and several other settings.
 
 \param T the tracker_config_t to adjust
 \param num_variables The number of variables in the current system
 */
void adjust_tracker_AMP(tracker_config_t * T, int num_variables);




void generic_set_start_pts(point_data_d ** startPts,
						   const WitnessSet & W);

void generic_set_start_pts(point_data_mp ** startPts,
						   const WitnessSet & W);

void generic_setup_patch(patch_eval_data_d *P, const WitnessSet & W); // for mp type 0
void generic_setup_patch(patch_eval_data_mp *P, const WitnessSet & W);// for my type 2
void generic_setup_patch(patch_eval_data_mp *P, const WitnessSet & W, int prec); // for mp type 1

int generic_setup_files(FILE ** OUT, boost::filesystem::path outname,
						FILE ** MIDOUT, boost::filesystem::path midname);


/** 
 \brief reads the tracker_config_t from file, calling Bertini's setupConfig(), which reads the "config" file from disk.
 
 \param solve_options the current state of the Solver
 \param MPType the current operating MP mode.
 
 */
void get_tracker_config(SolverConfiguration &solve_options,int MPType);



/**
 \brief The main function call for a master process to call in order to track paths; ED must already be set up. 
 
 automatically switches between parallel and serial mode depending on solve_options.
 
 \param solve_out The way data is transferred out of this function
 \param W const input witness set, with points and patches
 \param ED_d a pointer to an already-set-up evaluator in double.
 \param ED_mp a pointer to an already-set-up evaluator in mp.
 \param solve_options the current state of the Solver.
 */
void master_solver(SolverOutput & solve_out, const WitnessSet & W,
				   SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
				   SolverConfiguration & solve_options);


/**
 \brief the mail loop for solving in serial.
 
 \param trackCount collects statistics regarding paths
 \param OUT open output file
 \param midOUT open file for midpath data
 \param endPoints the output from this function, the solutions
 \param W the input witness set, including start points.
 \param ED_d already populated abstract evaluator data in double.
 \param ED_mp already populated abstract evaluator data in mp.
 \param solve_options The current state of the Solver config.
 */
void serial_tracker_loop(trackingStats *trackCount,
						  FILE * OUT, FILE * midOUT,
						  const WitnessSet & W,  // was the startpts file pointer.
						  post_process_t *endPoints,
						  SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
						  SolverConfiguration & solve_options);

/**
 \brief the main loop for solving in parallel.
 
 \param trackCount collects statistics regarding paths
 \param OUT open output file
 \param MIDOUT open file for midpath data
 \param W the input witness set, including start points.
 \param endPoints the output from this function, the solutions
 \param ED_d already populated abstract evaluator data in double.
 \param ED_mp already populated abstract evaluator data in mp.
 \param solve_options The current state of the Solver config.
 */
void master_tracker_loop(trackingStats *trackCount,
						 FILE * OUT, FILE * MIDOUT,
						 const WitnessSet & W,  // was the startpts file pointer.
						 post_process_t *endPoints,
						 SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
						 SolverConfiguration & solve_options);


/**
 \brief figure out how many paths to pass out this time
 
 this number is calculated as n = int( 1 + (( int(1 + ((num_points - 1) / num_workers)) - 1) / 10))
 and uses rounding.
 
 \return the number of packets of data to be sent, using the method.
 \param num_workers how many workers there are in this ring.
 \param num_points how many paths there are left.
 */
int get_num_at_a_time(int num_workers, int num_points);



/**
 \brief send start points to a worker
 
 
 \param next_worker the id of the worker to whom to send the work.
 \param num_packets The number of start points to send
 \param startPts_d the pointers to the data to send, in double
 \param startPts_mp the pointers to the data to send, in MP
 \param next_index the bottom index.  points sent will be consecutive, always.
 \param solve_options the current state of the Solver config.
 */
void send_start_points(int next_worker, int num_packets,
					   point_data_d *startPts_d,
					   point_data_mp *startPts_mp,
					   int & next_index,
					   SolverConfiguration & solve_options);


/**
 \brief as a mater, receive solution points from a worker
 
 
 also checks for isSoln at the same time, after has received the end points.
 
 \return the MPI ID of the worker received from.
 \param trackCount collects statistics
 \param EG_receives a pointer to a set of pointers into which to receive the data
 \param max_incoming the largest number of points received so far.
 \param solution_counter count the number of successful tracks
 \param endPoints the container into which to receive the data.
 \param ED_d pointer to the double evaluator_data
 \param ED_mp pointer to the mp evaluator_data
 \param solve_options The current state of the Solver config.
 */
int receive_endpoints(trackingStats *trackCount,
					  endgame_data_t **EG_receives, int & max_incoming,
					  int & solution_counter,
					  post_process_t *endPoints,
					  SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
					  SolverConfiguration & solve_options);



/**
 \brief as a worker, receive starts, track, and send the solutions
 
 \param trackCount keeps track of statistics
 \param OUT open file into which we can write
 \param MIDOUT open file into which to print midpath data.
 \param ED_d double format evaluator data.
 \param ED_mp MP format evaluator data.
 \param solve_options the current state of the Solver config.
 */
void worker_tracker_loop(trackingStats *trackCount,
						 FILE * OUT, FILE * MIDOUT,
						 SolverDoublePrecision * ED_d, SolverMultiplePrecision * ED_mp,
						 SolverConfiguration & solve_options);



/**
 \brief a fairly low level functon, which calls trackpath from bertini
 
 \param pathNum the ID of the path
 \param EG_out the output data structure
 \param Pin the input point data, double format.  may be NULL
 \param Pin_mp the input point data, mp format, may be NULL depending on MPtype
 \param OUT open file into which to print data
 \param MIDOUT open file into which to print midpath data.
 \param T the current tracker config
 \param ED_d double format evaluator data
 \param ED_mp mp format evaluator data
 \param eval_func_d pointer to the double evaluator function.
 \param eval_func_mp pointer to the mp evaluator function.
 \param change_prec pointer to the precision changing function
 \param find_dehom pointer to the dehomogenizing function
 */
void generic_track_path(int pathNum, endgame_data_t *EG_out,
						point_data_d *Pin, point_data_mp *Pin_mp,
						FILE *OUT, FILE *MIDOUT,
						tracker_config_t *T,
						void const *ED_d, void const *ED_mp,
						int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
						int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
						int (*change_prec)(void const *, int),
						int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));

/**
 \brief a fairly low level functon, which calls trackpath from bertini inside of a loop which checks retVal, and adjusts T with reruns as necessary, to get the path to declare success
 
 \param pathNum the ID of the path
 \param EG_out the output data structure
 \param Pin the input point data, double format.  may be NULL
 \param Pin_mp the input point data, mp format, may be NULL depending on MPtype
 \param OUT open file into which to print data
 \param MIDOUT open file into which to print midpath data.
 \param solve_options the current state of the Solver
 \param ED_d double format Solver derived type
 \param ED_mp mp format Solver derived type
 \param eval_func_d pointer to the double evaluator function.
 \param eval_func_mp pointer to the mp evaluator function.
 \param change_prec pointer to the precision changing function
 \param find_dehom pointer to the dehomogenizing function
 */

void robust_track_path(int pathNum, endgame_data_t *EG_out,
					   point_data_d *Pin, point_data_mp *Pin_mp,
					   FILE *OUT, FILE *MIDOUT,
					   SolverConfiguration & solve_options,
					   SolverDoublePrecision *ED_d, SolverMultiplePrecision *ED_mp,
					   int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
					   int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
					   int (*change_prec)(void const *, int),
					   int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *));














#endif

