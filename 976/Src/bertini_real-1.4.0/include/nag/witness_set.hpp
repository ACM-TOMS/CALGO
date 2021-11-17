#pragma once

#include "containers/holders.hpp"
#include "programConfiguration.hpp"
#include "io/color.hpp"

/** 
 \brief a woefully incomplete class to contain systems which bertini will parse.
 
 This class is intended to hold what Bertini would need to produce a straight-line program for evaluation.
 */
class Function
{
	std::string func;  ///< symbolic representation of function (straight from input file).
										 // this class is woefully incomplete.
};




/**
 \brief witness set holds points, patches, and linears, with the names of the variables.
 
 The witness set class collects points, patches, and linears into one object.  It offers methods for sorting for real points only, for sorting to contain only unique points.
 
 A witness set gets two numbers of variables, one is the total number appearing in it, and the other [more importantly] is the numebr of natural variables contained therein.  The witness set is assumed to be in correspondence to a variable group with a single leading homogenizing variable, and automatically dehomogenizes points for uniqueness and reality testing.
 */
class WitnessSet : public PatchHolder, public LinearHolder, public PointHolder, public NameHolder
{
	
protected:
	
	//begin data members
	

	int dim_;
	int comp_num_;
	int incid_num_;
	
	int num_vars_;
	int num_natty_vars_;
	


	boost::filesystem::path input_filename_;
	Function input_file_;
	
	
	// end data members
	
	
public:
	
	/**
	 \brief get the name of the bertini input file for this witness set
	 
	 \return the path of the file
	 */
	inline boost::filesystem::path input_filename() const
	{
		return input_filename_;
	}
	
	/**
	 \brief set the name of the bertini input file
	 
	 \param new_input_filename The new name of the file
	 */
	void set_input_filename(boost::filesystem::path new_input_filename)
	{
		input_filename_ = new_input_filename;
	}
	
	/**
	 \brief get the dimension of the set represented by the witness set.
	 
	 \return the integer dimension of the component.
	 */
	inline int dimension() const
	{
		return dim_;
	}
	
	/**
	 \brief set the dimension of the witness set
	 
	 \param new_dim the dimension of the set
	 */
	void set_dimension(int new_dim)
	{
		dim_ = new_dim;
	}
	
	
	
	
	
	
	/**
	 \brief get the component number of the set represented by the witness set.
	 
	 \return the index the component.
	 */
	inline int component_number() const
	{
		return comp_num_;
	}
	
	
	/**
	 \brief sets the component number
	 
	 \param new_comp_num The new component number to set.
	 */
	void set_component_number(int new_comp_num)
	{
		comp_num_ = new_comp_num;
	}
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of variables in the set
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_vars_;
	}
	
	
	
	/**
	 \brief sets the total number of variables for the set
	 
	 \param new_num_vars the new total number of variables
	 */
	void set_num_variables(int new_num_vars)
	{
		num_vars_ = new_num_vars;
	}
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of natural variables in the set (those in the first group, including the homogenizing variable if present).
	 
	 \return the number of natural variables (those in the first variable group.
	 */
	inline int num_natural_variables() const
	{
		return num_natty_vars_;
	}
	
	/**
	 \brief set the number of natural variables.  
	 
	 \param new_num_nat_vars The new number of natural variables
	 */
	void set_num_natural_variables(int new_num_nat_vars)
	{
		num_natty_vars_ = new_num_nat_vars;
	}
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the incidence number for the witness set.  by default, is negative value.
	 
	 \return the incidence number, which is used for reading membership from bertini membership testing.
	 */
	inline int incidence_number() const
	{
		return incid_num_;
	}
	
	
	
	/**
	 \brief set the incidence number, probably after having determined it somehow.
	 
	 \param new_incidence the new incidence number to assign to the witness set.
	 */
	void set_incidence_number(int new_incidence)
	{
		incid_num_ = new_incidence;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// overloaded operators
	
	
	// default constructor
	
	WitnessSet(int nvar)
	{
		init();
		num_vars_ = num_natty_vars_ = nvar;
	};
	
	WitnessSet(){
		init();
	}
	
	
	~WitnessSet(){ // the destructor
		
		clear();
		
	}
	
	
	// assignment
	/**
	 \brief custom assignment call.
	 
	 \return the new witness set
	 \param other the other witness set to copy from.
	 */
	WitnessSet& operator=( const WitnessSet& other)
	{
		reset();
		copy(other);
    return *this;
  }
	
	
	/**
	 \brief copy operator.  
	 
	 must be explicitly declared because the underlying c structures use pointers.
	 
	 \param other the other witness set to copy from.
	 */
	WitnessSet(const WitnessSet & other)
	{
		init();
		copy(other);
	}
	
	
	/**
	 \brief get her ready for action.
	 */
	void init()
	{
		input_filename_ = "unset_filename";
		
		num_vars_ = 0;
		num_natty_vars_ = 0;
		
		
		incid_num_ = -1;
		comp_num_ = dim_ = -1;
		
		
		
	}
	
	
	/** 
	 \brief perform a total deep copy of the witness set
	 
	 \param other the other witness set to copy from.
	 */
	void copy(const WitnessSet & other)
	{
		
		copy_skeleton(other);

		copy_names(other);
		copy_points(other);
		copy_patches(other);
		copy_linears(other);
	}
	
	
	
	/**
	 \brief copy only the WitnessSet data members, but not any inherited members.
	 
	 \param other the other witness set to copy from.
	 */
    void copy_skeleton(const WitnessSet & other)
	{
		this->input_filename_ = other.input_filename_;
		
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		this->incid_num_ = other.incid_num_;
		
		this->num_vars_ = other.num_vars_;
		this->num_natty_vars_ = other.num_natty_vars_;
	}

    

    

	/**
	 get the number of synthetic variables
	 
	 \return the number of synthetic variables held.
	 */
    int num_synth_vars() const
    {
        return num_vars_ - num_natty_vars_;
    }
    
	
	/**
	 \brief get rid of any variables which are synthetic
	 */
	void only_natural_vars();
	
	/**
	 \brief trim off all but a number of variables.  trims linears, patches, and points.
	 
	 \param num_vars the number of variables to keep
	 */
	void only_first_vars(int num_vars);
	
	/**
	 \brief keep only the real points (upon dehomogenization)
	 
	 \param T the current state of the tracker, for the real tolerance.
	 */
	void sort_for_real(tracker_config_t * T);
	
	/**
	 \brief keep only the uniqie points (upon dehomogenization)
	 
	 \param T the current state of the tracker, for the unique tolerance.
	 */
	void sort_for_unique(tracker_config_t * T);
	
	
	
	/**
	 \brief keep only the points which are inside the sphere (upon dehomogenization)
	 
	 \param radius The radius of the sphere
	 \param center The center of the sphere.
	 */
	void sort_for_inside_sphere(comp_mp radius, vec_mp center);
	
	
	/**
	 \brief read in the witness set from a file, which MUST be formatted correctly.  for details on the format, see this code, or \see print_to_file()
	 
	 \return 0, for no good reason.
	 \param witness_set_file the path of the file to parse into this object
	 \param num_vars the number of variables in the witness set.  sadly, you have to set this manually at the time, as the header does not contain the information.  the is due to Bertini reasons.
	 */
	int  Parse(const boost::filesystem::path witness_set_file, const int num_vars);
	
	
	/**
	 \brief empty the set.  calls clear()
	 */
	void reset()
	{
		clear();
	}
	
	
	
	
	/**
	 \brief clear the entire contents of the witness set
	 */
	void clear()
	{
		
		reset_names();
		
		reset_points();
		
		reset_linears();
		
		reset_patches();

		init();
	}
	
	

	
	

	

	
	/**
	 \brief merges another witness set into this, not checking for uniqueness at all.
	 
	 Should you want to straight-up merge the contents of two witness sets with the same number of natural variables, you may, using this function.  All the points, linears, and patches will be copied from the input into the existing one on which you call this function.
	 
	 \param W_in The witness set containing data you want to copy.
	 \param T Bertini's tracker_config_t object, with settings for telling whether two points are the same.
	 */
	void merge(const WitnessSet & W_in, tracker_config_t * T);
	
	

	
	/**
	 \brief prints some information about the witness set to the screen
	 
	 Print variable information, linears, and patches to screen  
	 This is potentially a very large amount of data depending on the set, and should be done sparingly
	 */
	void print_to_screen() const;
	
	/**
	 \brief print the witness set into a file, which can be read back in to the same format.
	 
	 Print the witness set to a file of filename
	 
	 the format is:
	 
	 [
	 num_points dim comp_num
	 
	 point 1 - in standard bertini format
	 
	 point 2 - in standard bertini format
	 
	 ...
	 
	 last point - in standard bertini format
	 
	 num_linears num_variables
	 
	 linear 1 - in standard bertini point/vec format
	 
	 ...
	 
	 linear last - in standard bertini point/vec format.
	 
	 num_patches num_variables \todo this is incorrect.
	 
	 patch 1 - in standard bertini point/vec format
	 
	 ...
	 
	 patch last - in standard bertini point/vec format.
	 ]
	 
	 
	 \param filename the name of the file to which to write.
	 */
	void print_to_file(boost::filesystem::path filename) const;
	
	/**
	 writes the linears in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void write_linears(boost::filesystem::path filename) const;
	
	/**
	 writes the patches in point form to file filename
	 
	 \param filename the name of the file to be written.
	 */
	void print_patches(boost::filesystem::path filename) const;
	
	/**
	 read patches from a suitably set up file.
	 
	 \param filename the name of the file to read.
	 */
	void read_patches_from_file(boost::filesystem::path filename);
	
	/**
	 write all the points to a file, without dehomogenizing.
	 
	 \param filename the name of the file to which to write.
	 */
	void write_homogeneous_coordinates(boost::filesystem::path filename) const;
	
	/**
	 write all the points to a file, first dehomogenizing.
	 
	 \param filename the name of the file to which to write.
	 */
	void write_dehomogenized_coordinates(boost::filesystem::path filename) const;
	
	/**
	 write all the points to a file, first dehomogenizing.
	 
	 \param filename The name of the file to which to write.
	 \param indices Set of indices of points to write to file.
	 */
	void write_dehomogenized_coordinates(boost::filesystem::path filename,std::set<unsigned int> indices) const;

	
	
	
	/**
	 individual send, relative to MPI_COMM_WORLD
	 
	 \param target the ID target of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void send(ParallelismConfig & mpi_config, int target) const;
	
	/**
	 individual receive, relative to MPI_COMM_WORLD
	 
	 \param source the ID source of the communication
	 \param mpi_config the current MPI state, as implemented in bertini_real
	 */
    void receive(int source, ParallelismConfig & mpi_config);
};










