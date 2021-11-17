#pragma once


#include "bertini1/bertini_extensions.hpp"
#include "programConfiguration.hpp"


/**
 \brief 0-cell 
 
 a bertini_real Vertex, a 0-cell.  contains a point, its projection values, its type, whether its been removed, and an index into a set of filenames contained in the VertexSet.
 
 \todo remove the metadata from this, and instead track it in the Vertex set, much like the solver_output
 */
class Vertex
{

private:
	point_mp pt_mp_; ///< the main data for this class -- the point.
	
	
	vec_mp  projection_values_; ///< a vector containing the projection values.
	
	VertexType type_;  ///< See enum.
	int input_filename_index_; ///< index into the vertex_set's vector of filenames.
	
public:
	
	
	/**
	 \brief get the index of the originating file name
	 
	 \return the integer index
	 */
	inline int input_filename_index() const
	{
		return input_filename_index_;
	}
	
	
	/**
	 \brief set the index
	 
	 \param new_index the new index to set in the Vertex
	 */
	void set_input_filename_index(int new_index)
	{
		input_filename_index_ = new_index;
	}
	
	
	/**
	 \brief get the projection values.
	 \return the projection values in vec_mp form
	 */
	inline vec_mp& projection_values()
	{
		return projection_values_;
	}
	
	
	const vec_mp& get_projection_values() const
	{
		return projection_values_;
	}
	
	
	
	/**
	 \brief set the type of the Vertex, completely overwriting the old one.
	 
	 \param new_type the new type for the Vertex
	 */
	void set_type(VertexType const& new_type)
	{
		type_ = new_type;
	}
	
	/**
	 \brief set the type of the Vertex
	 
	 \param new_type the new type for the Vertex
	 */
	void add_type(VertexType const& new_type)
	{
		type_ |= new_type;
	}

	/**
	 \brief get the type of the Vertex
	 
	 \return the type, in integer form
	 */
	VertexType type() const
	{
		return type_;
	}
	
	/**
	 \brief set the type of the Vertex
	 
	 \param new_type the new type for the Vertex
	 */
	void remove_type(VertexType const& new_type)
	{
		type_ ^= new_type;
	}
	
	/**
	 \brief set the Vertex to be 'removed'
	 
	 \param new_val the new value for the flag
	 */
	void set_removed(bool new_val)
	{
		if (new_val)
			type_ |= Removed;
		else
			type_ ^= Removed;
	}
	
	
	/**
	 \brief query whether the Vertex has been set to 'removed'
	 
	 \return whether it has been removed.
	 */
	bool is_removed() const
	{
		return is_type(Removed);
	}
	
	/**
	 \brief query whether the Vertex has been set to 'removed'
	 
	 \return whether it has been removed.
	 */
	bool is_type(VertexType test_type) const
	{
		return type_ & test_type;
	}


	Vertex()
	{
		init();
	}
	
	~Vertex()
	{
		clear();
		
	}
	
	Vertex & operator=(const Vertex & other)
	{
		copy(other);
		return *this;
	}
	
	Vertex(const Vertex& other)
	{
		init();
		copy(other);
	}
	
	/**
	 /brief prints the Vertex to the screen
	 Prints the Vertex to the screen
	 */
	void print() const
	{
		print_point_to_screen_matlab(pt_mp_,"point");
		print_point_to_screen_matlab(projection_values_,"projection_values");
		std::cout << "type: " << type_ << std::endl;
	}
	
	
	/**
	 \brief sets the point.
	 
	 Set the vertex's point to the input.
	 \param new_point the input point to be set.
	 */
	void set_point(const vec_mp new_point);
	
	
	const vec_mp & get_point() const
	{
		return pt_mp_;
	}
	
	
	/**
	 \brief get the point in the Vertex
	 
	 \return the point in vec_mp form
	 */
	vec_mp& point()
	{
		return pt_mp_;
	}
	
	
	/**
	 \brief get the point in the Vertex
	 
	 \return the point in vec_mp form
	 */
	const vec_mp& point() const
	{
		return pt_mp_;
	}


	/**
	 \brief single target mpi send.
	 
	 Send the Vertex to a single target.
	 
	 \see Vertex::receive
	 
	 \param target The MPI ID of the target for this send
	 \param mpi_config current mpi settings
	 */
	void send(int target, ParallelismConfig & mpi_config) const;
	
	
	/**
	 \brief single source receive.
	 
	 Receive a Vertex from a single source.
	 
	 \see Vertex::send
	 
	 \param source The MPI ID of hte source.
	 \param mpi_config current mpi settings
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
private:
	
	void clear()
	{
		clear_vec_mp(this->pt_mp_);
		clear_vec_mp(this->projection_values_);
	}
	
	void copy(const Vertex & other)
	{
		set_point(other.pt_mp_);
		
		vec_cp_mp(this->projection_values_, other.projection_values_);
		this->type_ = other.type_;
		
		this->input_filename_index_ = other.input_filename_index_;
	}
	
	
	void init()
	{
		init_point_mp2(this->projection_values_,0,1024);
		init_point_mp2(this->pt_mp_,1,64);
		this->pt_mp_->size = 1;
		set_zero_mp(&pt_mp_->coord[0]);
		this->projection_values_->size = 0;
		this->type_ = Unset;
		
		this->input_filename_index_ = -1;
	}
};



