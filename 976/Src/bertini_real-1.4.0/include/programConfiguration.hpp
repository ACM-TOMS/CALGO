#ifndef PROGRAM_CONFIGURATION
#define PROGRAM_CONFIGURATION

/** \file programConfiguration.hpp */
#include "config.h"

#include <map>
#include <getopt.h> 
#include <queue>

#include <algorithm>
#include <boost/timer/timer.hpp>


#include "bertini1/bertini_extensions.hpp"

#include "io/fileops.hpp"


enum {BERTINIREAL=-9000,CRIT=-8999};

enum class SymEngine{Matlab, Python}; // this is to indicate which symbolic engine will run the decompositions; currently default is Matlab

///////////
//
//    PROGRAM CONFIGURATION
//
//////////


/**
 \brief Carries the configuration of MPI, including ID number, communicator, id of headnode, etc.
 
 
 */
class ParallelismConfig
{
	
protected:
	
	bool force_no_parallel_ = false;
	int headnode_ = 0;
	int my_id_, my_id_global_;
	int numprocs_;
	MPI_Comm   my_communicator_;
	
	int worker_level_; // higher worker level means more tedious work, in a sense.  worker_level 0 is uber-master.  worker_level 1 will be the next level down in management, etc.  the exact usage of this is relative to the process being run.
	
	
	
	std::map< int, int> worker_status_;
	
	std::queue< int > available_workers_;
	
	
	
public:
	
	/**
	 \brief get the state of whether forcing no parallel at the current time
	 \return indicator of whether parallelism currently turned off
	 */
	bool force_no_parallel()
	{
		return force_no_parallel_;
	}
	
	
	/**
	 \brief set the value of force_no_parallel, to turn on/off parallelism
	 
	 \param new_val the new value to set it to
	 */
	void force_no_parallel(bool new_val)
	{
		force_no_parallel_ = new_val;
	}
	
	
	
	/**
	 default constructor
	 */
	ParallelismConfig(){
		init();
	}
	
	
	/**
	 \brief determine whether the process is the head.
	 
	 \return Indicator of whether the process thinks it is currently the head.
	 */
	inline bool is_head() const 
	{
		if (my_id_ == headnode_)
			return true;
		else
			return false;
	}
	
	
	/**
	 
	 \brief Determine whether the process thinks it should be in parallel mode.
	 
	 \return Indicator of whether to use parallel mode.
	 */
	inline bool use_parallel() const {
		if (numprocs_>1 && force_no_parallel_!=true)
			return true;
		else
			return false;
	}
	
	
	
	
	/**
	 \brief Call MPI_Abort using the internally stored communicator.
	 
	 \param why An integer to feed to MPI_Abort.
	 */
	void abort(int why){
		MPI_Abort(my_communicator_,why);
	}
	
	/**
	 Get the current communicator
	 \return the currently stored communicator
	 */
	inline MPI_Comm comm() const {return my_communicator_;}
	
	/**
	 \brief Get the ID of the supervisor
	 
	 \return the ID of the head node, supervisor, or whatever you want to call it.
	 */
	inline int head() const {return headnode_;}
	
	/**
	 \brief Get the ID of this process
	 
	 \return the ID
	 */
	inline int id() const { return my_id_;}
	
	/**
	 \brief  Get the worker level
	 
	 \return the worker level
	 */
	inline int level() const {return worker_level_;}
	
	/**
	 \brief Get how many workers there are in the current communicator
	 
	 \return The number of workers.
	 */
	inline int num_procs() const { return numprocs_;}
	
	
	
	
	/**
	 \brief Get how many workers there are in the current communicator
	 
	 \return The number of workers.
	 */
	inline int size() const { return numprocs_;}
	
	
	/**
	 \brief Set up the vector of available workers, based on how many processors there are.
	 \throws logic_error if the set of active workers is not empty.
	 
	 It also sets up the workers to be listed as inactive.
	 */
	void init_active_workers();
	
	/**
	 \brief Get the next available worker, relist it as active, and return its id.
	 
	 Available workers are stored as a queue of integers, and work is assigned to the front of the vector.  This method pops the front entry of the queue, relists it as active, and returns its ID.
	 \return the ID of the next worker.
	 */
	int activate_next_worker();
	
	/**
	 \brief Take an active worker, relist it as inactive, and make it available.  If the worker is not active, it calls MPI_Abort
	 
	 \param worker_id The id of the worker to deactivate.
	 */
	void deactivate(int worker_id);
	
	
	/**
	 \brief Individually send a number to all available workers.
	 
	 pops available workers off the queue, and deliveres the them numtosend.
	 
	 \param numtosend The integer number to send.
	 */
	void send_all_available(int numtosend);
	
	
	/**
	 \brief MPI_Bcast a number to everyone in the ring, calling for help.
	 
	 \param solver_type The case-index of the type of help the head wants.
	 */
	void call_for_help(int solver_type);
	
	/**
	 \brief check if there are available workers.
	 
	 \return a boolean indicating if there are available workers.
	 */
	bool have_available();	
	/**
	 \brief check if there are workers working.
	 
	 \return A boolean indicating whether there are workers with the status 'ACTIVE'.
	 */
	bool have_active();
	
	/**
	 \brief Get the number of workers listed as active.
	 
	 \return the number of workers listed as Active.
	 */
	int num_active();
	
	
private:
	
	/**
	\brief initialize the parallelism config

	called in the constructor.
	*/
	void init();
	
};



/**
 \brief Base class for program configuations.
 
 Both sampler_configuration and BertiniRealConfig inherit from this.
 */
class ProgramConfigBase : public ParallelismConfig
{
	
private:
	
	
	int verbose_level_ = 0;
	
	boost::filesystem::path called_dir_;
	boost::filesystem::path working_dir_;
	boost::filesystem::path output_dir_;
	
	boost::timer::cpu_timer timer_;
protected:
	
	/**
	 \brief set the name of the called directory
	 \param new_name the new name of the called directory.
	 */
	void set_called_dir(boost::filesystem::path new_name)
	{
		called_dir_ = new_name;
	}
	
public:
	
	/**
	 get the directory from which the program was called
	 
	 \return the directory in which the user was, when the ProgramConfigBase was created.
	 */
	boost::filesystem::path called_dir() const
	{
		return called_dir_;
	}
	
	
	/**
	 get the working directory
	 
	 \return the working directory
	 */
	boost::filesystem::path working_dir() const
	{
		return working_dir_;
	}
	
	
	/**
	 \brief set the name of the working directory
	 
	 \param new_name the new name of the working directory
	 */
	void working_dir(boost::filesystem::path new_name)
	{
		working_dir_ = new_name;
	}
	
	
	
	/**
	 get the output directory
	 
	 \return the output directory
	 */
	boost::filesystem::path output_dir() const
	{
		return output_dir_;
	}
	
	
	/**
	 \brief set the name of the output directory
	 
	 \param new_name the new name of the output directory
	 */
	void output_dir(boost::filesystem::path new_name)
	{
		output_dir_ = new_name;
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
	
	
	void move_to_temp();
	
	void move_to_called();
	
	/**
	\brief Write some structured meta data to a file, including things like runtime, number of processors used, and version of program used.
	*/
	void PrintMetadata(boost::filesystem::path const& filename) const;

};


/**
\brief Write a file contanining the indices of vertex types.
*/
void PrintPointTypeMapping(boost::filesystem::path const& filename);

/**
 \brief holds the current state of configuration for Bertini_real.
 
 */
class BertiniRealConfig : public ProgramConfigBase
{
	bool orthogonal_projection_;
	bool compute_cycle_numbers_ = false; ///< whether we should compute cycle numbers for edges and faces
	bool debugwait_; ///< flag for whether to wait 30 seconds before starting, and print the master process ID to screen.
	int max_deflations_; ///< the maximum allowable number of deflation iterations before it gives up.
	
	bool stifle_membership_screen_; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text_; ///<
	
	int quick_run_ = false;  ///< indicator of whether to use the robust solver wherever possible
	
	
	bool user_sphere_; ///< flag for whether to read the sphere from a file, rather than compute it.
	
	bool user_projection_; ///< indicator for whether to read the projection from a file, rather than randomly choose it.
	
	bool merge_edges_; ///< a mode switch, indicates whether should be merging.
	
	
	
	int primary_mode_; ///< mode of operation -- bertini_real is default, but there is also crit method for computing critical points.
	
	
	boost::filesystem::path bounding_sphere_filename_; ///< name of file to read if user_sphere==true
	boost::filesystem::path projection_filename_; ///name of file to read if user_projection==true
	boost::filesystem::path input_filename_; ///< name of the input file to read -- by default it's "input"
	
	boost::filesystem::path input_deflated_filename_; ///< the name of the file post-deflation
	
	
	
	int target_dimension_;  ///< the dimension to shoot for
	int target_component_;  ///< the integer index of the component to decompose.  by default, it's -2, which indicates 'ask me'.
	
	
	std::string matlab_command_; ///< the string for how to call matlab.
	
	bool use_gamma_trick_; ///< indicator for whether to use the gamma trick in a particular solver.
	
        SymEngine engine_; ///< the symbolic class variable that indicates which symbolic engine the user desires. Default is currently Matlab	
	
	
	
public:
	
	/**
	 get the mode for the program.  by default, it's bertini_real
	 */
	int primary_mode() const
	{
		return primary_mode_;
	}
	
	/**
	 get the path to the input_deflated file.
	 
	 \return the path to the input_deflated file.
	 */
	boost::filesystem::path input_deflated_filename() const{
		return input_deflated_filename_;
	}
	
	
	
	/**
	 set the path to the input_deflated file.
	 
	 \param new_name the path to the input_deflated file.
	 */
	void set_input_deflated_filename(const boost::filesystem::path & new_name)
	{
		input_deflated_filename_ = new_name;
	}
	
	
	/**
	 get the path to the projection file.
	 
	 \return the path to the projection file.
	 */
	boost::filesystem::path projection_filename() const{
		return projection_filename_;
	}
	
	/**
	 get the path to the Bertini input file.
	 
	 \return the path to the input file.
	 */
	boost::filesystem::path input_filename() const{
		return input_filename_;
	}
	
	
	/**
	 get the path to the sphere file.
	 
	 \return the path to the sphere file.
	 */
	boost::filesystem::path bounding_sphere_filename() const{
		return bounding_sphere_filename_;
	}
	
	
	/**
	 get the command for calling Matlab via system(), which I hate doing anyway...  ugh.
	 
	 \return the string for calling Matlab.
	 */
	std::string matlab_command() const
	{
		return matlab_command_;
	}
	
	
	
	
	
	/**
	 get the target dimension.  by default, this is -1, which is 'ask the user'
	 \return the dimension.
	 */
	int target_dimension() const{
		return target_dimension_;
	}
	
	
	
	/**
	 set the target dimension.
	 \param new_dim the dimension.
	 */
	void set_target_dimension(int new_dim)
	{
		target_dimension_ = new_dim;
	}
	
	
	/**
	 get the number of the target component to decompose.
	 \return the target component.  default is -2, which is 'ask the user'
	 */
	int target_component() const
	{
		return target_component_;
	}
	
	/**
	 get whether are supposed to use the gamma trick.  default is no.
	 
	 \return whether using gamma trick
	 */
	bool use_gamma_trick() const
	{
		return use_gamma_trick_;
	}
	
	
	/** 
	 \brief query whether should merge edges
	 \return whether we should merge edges or not
	 */
	bool merge_edges()
	{
		return merge_edges_;
	}
	
	
	/**
	 \brief set whether should merge edges
	 \param new_val set whether we should merge edges or not
	 */
	void merge_edges(bool new_val)
	{
		merge_edges_ = new_val;
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 \brief whether to read projection from file.
	 \return whether to read the projection from a user-created and specified file.
	 */
	bool user_projection()
	{
		return user_projection_;
	}
	
	/**
	 \brief set whether to read projection from file.
	 \param new_val whether to read the projection from a user-created and specified file.
	 */
	void user_projection(bool new_val)
	{
		user_projection_ = new_val;
	}
	
	
	/**
	 \brief whether to read sphere parameters from file.
	 \return whether to read the sphere parameters from a user-created and specified file.
	 */
	bool user_sphere()
	{
		return user_sphere_;
	}
	
	
	
	/**
	 \brief set whether to read sphere parameters from file.
	 \param new_val whether to read the sphere parameters from a user-created and specified file.
	 */
	void user_sphere(bool new_val)
	{
		user_sphere_ = new_val;
	}
	
	/**
	 \brief get the level of quick.  higher == faster (less robust)
	 \return the level of quickness
	 */
	int quick_run()
	{
		return quick_run_;
	}
	
	/**
	 \brief set the level of quickness
	 \param new_val the new level
	 */
	void quick_run(int new_val)
	{
		quick_run_ = new_val;
	}
	
	
	
	
	/**
	 \brief the stifling text for system commands
	 \return string to append to system commands to stifle screen output.
	 */
	std::string stifle_text() const
	{
		return stifle_text_;
	}
	
	/**
	 \brief set the stifling text
	 \param new_val the new stifling text
	 */
	void stifle_text(std::string new_val)
	{
		stifle_text_ = new_val;
	}
	
	
	
	/**
	 \brief should we wait for 30 seconds before running program?
	 \return whether we should.  true==yes
	 */
	bool debugwait() const
	{
		return debugwait_;
	}
	
	/**
	 \brief set value for question -- should we wait for 30 seconds before running program?
	 \param new_val the new value for the debugwait parameter
	 */
	void debugwait(bool new_val)
	{
		debugwait_ = new_val;
	}
	
	
	/**
	 how many times it is ok to deflate
	 \return how many times to deflate, at maximum.  default is 10.
	 */
	int max_deflations() const
	{
		return max_deflations_;
	}
	
	/**
	 \brief set the maximum number of deflations
	 \param new_val the new value for the maximum number of deflations
	 */
	void max_deflations(int new_val)
	{
		max_deflations_ = new_val;
	}
	
	
	/**
	 \brief query whether we should stifle some screen output
	 \return whether we should.
	 */
	bool stifle_membership_screen() const
	{
		return stifle_membership_screen_;
	}
	
	/**
	 \brief set whether should stifle the membership testing screen output
	 \param new_val The new value for setting, whether to stifle screen output.
	 */
	void stifle_membership_screen(bool new_val)
	{
		stifle_membership_screen_ = new_val;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 get whether should use orthogonal projection.  default is yes
	 \return whether we are using an orthogonal projection.  default internally is yes.
	 */
	bool orthogonal_projection()
	{
		return orthogonal_projection_;
	}

    /**
    \brief Get which symbolic engine is being used. 

    The default (currently) is Matlab
    
    \return whether we are using Matlab or Python. default is Matlab
    */

    SymEngine symbolic_engine() const
    {
		return engine_;
	}


   /** 
   \brief set which symbolic engine is being used. 
   */
  	void symbolic_engine(SymEngine new_engine)
	{
		engine_ = new_engine;
	}
	
	/** 
	 \brief  display options to user. */
	void print_usage();
	
	
	/** 
	 \brief get the BertiniRealConfig from the command line. 
	 
	 \return the number 0.
	 \param argC the command count, from main()
	 \param args the input command string from main()
	 */
	int parse_commandline(int argC, char *args[]);
	
	
	/** 
	 \brief check to make sure files are in place, etc.  
	 
	 checks for write priveledges, and for the existence of bounding_sphere_filename and projection_filename.
	 
	 
	 \return the number 0.
	 */
	int startup();
	
	
	/** 
	 \brief displays the bertini_real splash screen */
	void splash_screen();
	
	/** 
	 \brief prints the current configuration to the screen, and pauses. */
	void display_current_options();
	
	
	BertiniRealConfig() : ProgramConfigBase()
	{
		
        init();
	};
	
	void init();


}; //re: BertiniRealConfig




class sampler_configuration : public ProgramConfigBase
{
public:
	enum class Mode{
		Fixed,
		AdaptiveConsecDistance,
		AdaptivePredMovement
	};


	int stifle_membership_screen; ///< boolean controlling whether stifle_text is empty or " > /dev/null"
	std::string stifle_text; ///< the text to append to system() commands to stifle screen output
	
	int minimum_num_iterations; ///< the minimum number of passes for iterative adaptive sampling
	int maximum_num_iterations; ///< the maximum number of passes for iterative adaptive sampling
	
	int use_gamma_trick; ///< indicator for whether to use the gamma trick.
	mpf_t TOL; ///< the distance-tolerance for spatial-adaptive sampling
	
	bool no_duplicates; ///< a flag for whether to never duplicate points in the vertex_set as it is constructed.
	
	bool use_distance_condition; ///< switch for adaptive modes, between distance or movement breaking of while loop.
	Mode mode; ///< mode switch between adaptive and fixed-number.
	int target_num_samples; ///< the number of samples per cell, more or less.
	
	int max_num_ribs;
	int min_num_ribs;
	
	/** 
	 \brief get the sampler_configuration from the command line. */
	int  parse_commandline(int argc, char **argv);
	
	/**
	 \brief print a splash opening message, including the version number and author list
	 */
	void splash_screen();
	
	/**
	 \brief print a message to screen about how to use the program
	 */
	void print_usage();
	
	/**
	 \brief parse the command line for options, using getopt_long_*
	 */
	int  parse_options(int argc, char **argv);
	
	
	/**
	 \brief default constructor, contains some default settings.
	 */

	void SetDefaults();


	sampler_configuration()
	{
		SetDefaults();
	};
	
	~sampler_configuration()
	{
		mpf_clear(TOL);
	}

	
	
	
};






/**
 \brief splits the bertini input file into several files for later use.
 
 calls MPI_Bcast(&PARSING, 1, MPI_INT, 0, MPI_COMM_WORLD); to be able to let workers carry through.  sadly, the parse_input() method in Bertini calls an MPI_Bcast, which will trip up the workers if it is not caught.
 
 basically, the files made are
 • num.out
 • arr.out
 • deg.out
 • names.out
 • func_input
 • preproc_data
 and maybe others.
 
 \param filename the name of the file to parse.
 */
void parse_input_file(boost::filesystem::path filename);

/**
 \brief splits the bertini input file into several files for later use.
 
 calls MPI_Bcast(&PARSING, 1, MPI_INT, 0, MPI_COMM_WORLD); to be able to let workers carry through.  sadly, the parse_input() method in Bertini calls an MPI_Bcast, which will trip up the workers if it is not caught.
 
 basically, the files made are
 • num.out
 • arr.out
 • deg.out
 • names.out
 • func_input
 • preproc_data
 and maybe others.
 
 \param filename the name of the file to parse.
 \param MPType a set-integer by pointer, this function splits the file and gets the MPType
 */

void parse_input_file(boost::filesystem::path filename, int * MPType);


/**
 \brief a wrapper around setupPreProcData(), and populates a preproc_data
 
 \param filename the name of the file to parse.
 \param PPD the preproc_data to populate
 */
void parse_preproc_data(boost::filesystem::path filename, preproc_data *PPD);









/**
 \brief reads in projection from file if user specified, creates one otherwise.
 
 currently defaults to create a random real projection with homogeneous value 0;
 
 \param pi the projection vectors to fill.  must be initted already, but not necessarily the correct size.
 \param program_options The current state of Bertini_real.
 \param num_vars how many variables to set up, including the homogenizing variable.
 \param num_projections how many proj vectors to set up.  again, these must already be allocated outside this call.
 */
void get_projection(vec_mp *pi,
					BertiniRealConfig program_options,
					int num_vars,
					int num_projections);




















#endif


