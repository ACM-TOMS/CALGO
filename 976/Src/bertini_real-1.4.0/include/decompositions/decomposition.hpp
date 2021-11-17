#pragma once

#include "containers/holders.hpp"
#include "nag/system_randomizer.hpp"
#include "nag/witness_set.hpp"
#include "containers/vertex_set.hpp"

/**
 \brief base Decomposition class.  curves and surfaces inherit from this.
 
 The Decomposition class holds the basic information for any dimensional decomposition -- a witness set which generated it, the number of variables, the dimension, which component it represents, the projections, the randomizer, the sphere, the input file name.
 
 */
class Decomposition : public PatchHolder
{


	
	
	
public:
	
	
    
/**
 \brief add a projection vector to the decomposition
 
 Decompositions in Bertini_real are computed with respect to linear projections, which are stored as vectors. They are repeated in each decomposition.  
 
 \param proj the projection to add.
 */
	void add_projection(vec_mp proj){
		if (this->num_curr_projections_==0) {
			pi_ = (vec_mp *) br_malloc(sizeof(vec_mp));
		}
		else{
			this->pi_ = (vec_mp *)br_realloc(this->pi_, (this->num_curr_projections_+1) * sizeof(vec_mp));
		}
		
		init_vec_mp2(this->pi_[num_curr_projections_],proj->size,proj->curr_prec);
		this->pi_[num_curr_projections_]->size = proj->size;
		
		vec_cp_mp(pi_[num_curr_projections_], proj);
		num_curr_projections_++;
		
	}
	
	
	/**
	 \brief commit a set of points to the Vertex set, associating them with the input file for this decomposition.
	 
	 \todo rename this function to something more accurately descriptive
	 
	 \return the number 0.  this seems pointless.
	 \param W the witness set containing the points to add.
	 \param add_type the type the points will inherit.  {\em e.g.} CRITICAL
	 \param V the Vertex set to add the points to.
	 */
    int add_witness_set(const WitnessSet & W, VertexType add_type, VertexSet & V);
    
	
	/**
	 \brief Find the index of a testpoint.
	 
	 Search the Vertex set passed in for testpoint.  This happens relative to this decomposition for no reason whatsoever.
	 \todo change this to be a property of the Vertex set
	 
	 \return the index of the point, or -1 if not found.
	 \param V the Vertex set in which to search
	 \param testpoint the point for which to search
	 */
	int index_in_vertices(VertexSet &V,
                          vec_mp testpoint) const;
	
	
	/**
	 \brief search for a testvertex, and add it to Vertex set, and this decomposition, if not found.
	 
	 Search the passed in Vertex set for the testvertex -- and add it to the Vertex set, and its index to the decomposition if not found.
	 
	 \return the index of the point, or -1 if not found.
	 \param V the Vertex set in which to search
	 \param vert a Vertex with point for which to search.
	 */
	int index_in_vertices_with_add(VertexSet &V,
                                   Vertex vert);
	
	/** 
	 set up the base decomposition class from a file
	 
	 \return the number 0. seems useless.
	 \param INfile the name of the file to read from
	 */
	int setup(boost::filesystem::path INfile);
	
	
	/**
	 base method for printing decomposition to file.
	 
	 \param outputfile the name of the file to which to write.
	 */
	virtual void print(boost::filesystem::path outputfile) const;
	
	
	
	/**
	 read the bounding sphere from a file containing the radius and center coordinates.
	 
	 format is:
	 
	 [radius
	 
	 x_0
	 x_1
	 ...
	 x_n]
	 
	 
	 \return SUCCESSFUL value.
	 \param bounding_sphere_filename the name of the file.
	 */
	int read_sphere(const boost::filesystem::path & bounding_sphere_filename);
	
	
	/**
	 \brief set the sphere for this decomposition according to the input witness set.
	 
	 Because we need to capture parts of the component which go to inifinity, we intersect the component with a sphere containing all the critical points, which are the input to this method.  
	 
	 This method computes the centroid of the set -- which becomes the center of the sphere -- and the distance from the centroid to the outermost critical point -- 3 times which becomes the radius of the sphere.
	 
	 \param W_crit the witness set containing the critical points to capture inside the sphere
	 */
	void compute_sphere_bounds(const WitnessSet & W_crit);
	
	
	/**
	 \brief copy the bounds from another Decomposition
	 
	 Sub-decompositions will often want to inherit the bounding sphere of another.  This method lets you copy from one to another.
	 
	 \throws invalid_argument if the input Decomposition has no sphere yet
	 
	 \param other the Decomposition which already holds sphere bounds.
	 */
	void copy_sphere_bounds(const Decomposition & other)
	{
		if (!other.have_sphere_) {
			throw std::invalid_argument("trying to copy sphere bounds from a Decomposition which does not have them set!");
		}
		
		set_mp(this->sphere_radius_, other.sphere_radius_);
		vec_cp_mp(this->sphere_center_, other.sphere_center_);
		this->have_sphere_ = true;
	}
	
	
	/**
	 \brief the main way to print a Decomposition to a file.
	 
	 This method backs up the existing folder to one siffixed with "_bak", and creates a new folder with the correct name, to which it prints the Decomposition in text file format.
	 
	 \param base the base folder name to print the Decomposition.
	 */
	void output_main(const boost::filesystem::path base) const;
	
	
	/**
	 \brief copy the component number, filename, number of variables, and witness set itself into the Decomposition.
	 
	 \param W the witness set with the data.
	 */
	void copy_data_from_witness_set(const WitnessSet & W);
	
	
	/**
	 reset Decomposition to empty.
	 */
	void reset()
	{
		clear();
		init();
	}
	
	
	
	/**
	 default constructor.
	 */
	Decomposition(){
		init();
	}
	
	
	/**
	 default destructor.
	 */
	virtual ~Decomposition()
	{
		this->clear();
	}
	
	
	/**
	 assignment
	 \return The assigned decomposition.
	 \param other The input decomposition, from which to assign.
	 */
	Decomposition & operator=(const Decomposition& other){
		
		this->init();
		
		this->copy(other);
		
		return *this;
	}
	
	
	/**
	 copy-constructor
	 \param other Another decomposition from which to copy-construct.
	 */
	Decomposition(const Decomposition & other){
		
		this->init();
		
		this->copy(other);
	}
	
	
	/**
	 \brief single target MPI send.
	 
	 Send base Decomposition to another process.
	 
	 \param target the ID of the worker to which to send.
	 \param mpi_config the current configuration of MPI
	 */
	void send(int target, ParallelismConfig & mpi_config) const;
	
	/**
	 \brief single source MPI receive.
	 
	 Receive base Decomposition from another process.
	 
	 \param source the ID of the process from which to receive
	 \param mpi_config the current configuration of MPI
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	
	
	
	
	
	
	
	
	
	

	/**
	 \brief get a shared pointer to the randomizer within
	 
	 \return the shared pointer to the system randomizer
	 */
	std::shared_ptr<SystemRandomizer> randomizer() const
	{
		return randomizer_;
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 \brief get the number of variables in the Decomposition
	 
	 \return the number of variables
	 */
	inline int num_variables() const
	{
		return num_variables_;
	}
	
	/**
	 \brief set the number of variables
	 
	 \param new_num_variables the new number of variables
	 \return the number of variables
	 */
	int set_num_variables(int new_num_variables)
	{
		return num_variables_ = new_num_variables;
	}
	
	
	
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
	 \brief get the witness set associated with the Decomposition
	 
	 \return the witness set which generated the Decomposition
	 */
	const WitnessSet& get_W() const
	{
		return W_;
	}
	
	
	/**
	 \brief get the number of current projections
	 \return the number of currently held projections
	 */
	inline int num_curr_projections() const
	{
		return num_curr_projections_;
	}
	
	
	/**
	 \brief get a pointer to the beginning of the array of projections.  
	 \throws out of range if set of projections is empty when this is requested
	 \return pointer to the 0th projection.  will be NULL if have no projections.
	 */
	inline vec_mp& pi() const
	{
		if (num_curr_projections_>0) {
			return pi_[0];
		}
		else
		{
			throw std::out_of_range("trying get pointer for projections, but have no projections");
		}
	}
	
	
	/**
	 \brief get a pointer to the ith projections.
	 \throws out of range if trying to get out of range projection
	 \param index The index of the projection you want to get.
	 \return pointer to the ith projection.  will throw if out of range
	 */
	inline vec_mp& pi(int index) const
	{
		if (index >= num_curr_projections_) {
			throw std::out_of_range("trying to access an out of range projection in Decomposition");
		}
		else
		{
			return pi_[index];
		}
	}
	
	
	
	/**
	 \brief test for whether the sphere is set
	 \return Whether have the sphere radius and center set, whether from file or computed from a set of points.
	 */
	inline bool have_sphere() const
	{
		return have_sphere_;
	}
	
	/**
	 get pointer to the sphere radius
	 
	 \return pointer to the sphere radius
	 */
	comp_mp& sphere_radius()
	{
		return sphere_radius_;
	}
	
	/**
	 \brief get pointer to the sphere's center
	 
	 \return pointer to the sphere's center
	 */
	vec_mp& sphere_center()
	{
		return sphere_center_;
	}
	
	
	/**
	 \brief set the radius of the sphere
	 \param new_radius the radius of the sphere
	 */
	void set_sphere_radius(comp_mp new_radius)
	{
		set_mp(sphere_radius_, new_radius);
	}
	
	
	/**
	 \brief set the center of the sphere
	 \param new_center the center of the sphere
	 */
	void set_sphere_center(vec_mp new_center)
	{
		vec_cp_mp(sphere_center_,new_center);
	}
	
	
	
	/**
	 make a deep copy of another Decomposition
	 
	 \param other The other Decomposition, from which to clone into this.
	 */
	void clone(const Decomposition & other)
	{
		PatchHolder::copy(other);
		
		
		this->randomizer_ = std::make_shared<SystemRandomizer>( *other.randomizer_ );
		
		this->W_ = other.W_;
		
		this->input_filename_ = other.input_filename_;
		
		this->num_variables_ = other.num_variables_;
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		
		
		if (this->num_curr_projections_==0) {
			this->pi_ = (vec_mp *) br_malloc(other.num_curr_projections_ * sizeof(vec_mp));
		}
		else{
			for (int ii=0; ii<num_curr_projections_; ii++) {
				clear_vec_mp(pi_[ii]);
			}
			this->pi_ = (vec_mp *) br_realloc(this->pi_,other.num_curr_projections_ * sizeof(vec_mp));
		}
		
		this->num_curr_projections_ = other.num_curr_projections_;
		for (int ii = 0; ii<other.num_curr_projections_; ii++) {
			init_vec_mp2(this->pi_[ii],other.pi_[ii]->size,other.pi_[ii]->curr_prec);
			this->pi_[ii]->size = other.pi_[ii]->size;
			vec_cp_mp(this->pi_[ii], other.pi_[ii])
		}
		
		
		
		
		
		copy_sphere_bounds(other);

	}
	
	
	void SetCritSliceValues(vec_mp & p)
	{	
		vec_cp_mp(crit_slice_values, p);
	}

	vec_mp& CritSliceValues() 
	{
		return crit_slice_values;
	}

protected:
	
	vec_mp crit_slice_values;
	vec_mp sphere_center_; ///< the center of the sphere.
	comp_mp sphere_radius_; ///< the radius of the sphere.
	bool have_sphere_; ///< indicates whether the Decomposition has the radius set, or needs one still.
	
	
	
	int num_curr_projections_; ///< the number of projections stored in the Decomposition.  should match the dimension when complete.
	vec_mp	*pi_; ///< the projections used to decompose.  first ones are used to decompose nested objects.
	
	
	

	
	
	
	/**
	 \brief set the witness set.
	 
	 \param new_w the new witness set to set.
	 */
	void set_W(const WitnessSet & new_w)
	{
		W_ = new_w;
	}
	
	
	
	WitnessSet W_; ///< generating witness set
	
	boost::filesystem::path input_filename_; ///< the name of the text file in which the system resides.
	
	
	int dim_; ///< the dimension of the Decomposition
	int comp_num_; ///< the component number.
	
	
	//	function input_file;
	int num_variables_; ///< the number of variables in the Decomposition
	
	
	std::shared_ptr<SystemRandomizer> randomizer_; ///< the randomizer for the Decomposition.
	

	
	
	int add_vertex(VertexSet &V, Vertex source_vertex);
	
	
	
	void init(){
		
		randomizer_ = std::make_shared<SystemRandomizer> (*(new SystemRandomizer()));

		input_filename_ = "unset";
		pi_ = NULL;
		
		
		num_curr_projections_ = 0;
		num_variables_ = 0;
		dim_ = -1;
		comp_num_ = -1;
		
		init_mp2(sphere_radius_,1024);
		init_vec_mp2(sphere_center_,0,1024);
		sphere_center_->size = 0;
		have_sphere_ = false;
		
		set_one_mp(sphere_radius_);
		neg_mp(sphere_radius_,sphere_radius_);

		init_vec_mp(crit_slice_values,0);
		crit_slice_values->size = 0;
	}
	
	
	
	
	void copy(const Decomposition & other)
	{
		PatchHolder::copy(other);
		
		
		this->randomizer_ = other.randomizer_;
		
		this->W_ = other.W_;
		
		this->input_filename_ = other.input_filename_;
		
		this->num_variables_ = other.num_variables_;
		this->dim_ = other.dim_;
		this->comp_num_ = other.comp_num_;
		
		
		if (this->num_curr_projections_==0) {
			this->pi_ = (vec_mp *) br_malloc(other.num_curr_projections_ * sizeof(vec_mp));
		}
		else{
			for (int ii=0; ii<num_curr_projections_; ii++) {
				clear_vec_mp(pi_[ii]);
			}
			this->pi_ = (vec_mp *) br_realloc(this->pi_,other.num_curr_projections_ * sizeof(vec_mp));
		}
		
		this->num_curr_projections_ = other.num_curr_projections_;
		for (int ii = 0; ii<other.num_curr_projections_; ii++) {
			init_vec_mp2(this->pi_[ii],other.pi_[ii]->size,1024);
			this->pi_[ii]->size = other.pi_[ii]->size;
			vec_cp_mp(this->pi_[ii], other.pi_[ii])
		}
		
		
		if (other.have_sphere_) {
			copy_sphere_bounds(other);
		}
		else{
			this->have_sphere_ = false;
		}
		
		
		return;
	}
	
	void clear()
	{
		if (num_curr_projections_>0){
			for (int ii=0; ii<num_curr_projections_; ii++)
				clear_vec_mp(pi_[ii]);
			free(pi_);
		}
		num_curr_projections_ = 0;
		
		clear_mp(sphere_radius_);
		clear_vec_mp(sphere_center_);
		
		clear_vec_mp(crit_slice_values);
	}
	
	
	
	

	
}; // end Decomposition




