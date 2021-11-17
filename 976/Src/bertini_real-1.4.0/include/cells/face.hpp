#pragma once

#include "cell.hpp"


/**
 \brief the 2-Cell
 
 The Face data type.  Extends the functionality of the basic Cell type, adding left and right indices of connected edges, top and bottom indices and system names, left and right projection values, and an index for which midslice it came from.  The midpoint member is inherited from the base class Cell.
 */
class Face : public Cell
{
	
	std::vector<int>	left_edges_;  ///< index into vertices
	std::vector<int>	right_edges_; ///< index into vertices
	
	int top_edge_index_;			///<  index of the top edge in the appropriate Decomposition, as indicated by system_name_top.
	int bottom_edge_index_;			///< index of the bottom edge in the appropriate Decomposition, as indicated by system_name_bottom.
	
	
	std::string system_name_bottom_; ///< the plain text name of the bottom edge's system.  e.g. crit_curve
	std::string system_name_top_; ///< the plain text name of the top edge's system.  e.g. crit_curve
	
	int crit_slice_index_; ///< which midpoint slice this Face came from.
	
	
	
	comp_mp left_crit_val_; ///< the pi_0 projection value of the left crit slice
	comp_mp right_crit_val_; ///< the pi_0 projection value of the right crit slice
	
public:
	
	/**
	 \brief get a pointer to the comp_mp of the left critical value.
	 \return pointer to a comp_mp of the critical value of the left edge.
	 */
	const comp_mp * left_crit_val() const
	{
		return & left_crit_val_;
	}
	
	/**
	 \brief set the projection value of the left critical edge
	 \param new_left_crit_val the new value of the projection for left edge
	 */
	void set_left_crit_val(comp_mp new_left_crit_val)
	{
		set_mp(left_crit_val_,new_left_crit_val);
	}
	
	/**
	 \brief get a pointer to the comp_mp of the right critical value.
	 \return pointer to a comp_mp of the critical value of the right edge.
	 */
	const comp_mp * right_crit_val() const
	{
		return & right_crit_val_;
	}
	
	
	/**
	 \brief set the projection value of the right critical edge
	 \param new_right_crit_val the new value of the projection for right edge
	 */
	void set_right_crit_val(comp_mp new_right_crit_val)
	{
		set_mp(right_crit_val_,new_right_crit_val);
	}
	
	
	/**
	 \brief add an edge to set of found edges.
	 \param index the index of the edge
	 */
	void add_left_edge(int index)
	{
		left_edges_.push_back(index);
	}
	
	
	
	/**
	 \brief add an edge to set of found edges.
	 \param index the index of the edge
	 */
	void add_right_edge(int index)
	{
		right_edges_.push_back(index);
	}
	
	
	
	
	/**
	 \brief get the bottom edge index
	 \return the index of the bottom edge
	 */
	int bottom_edge() const
	{
		return bottom_edge_index_;
	}
	
	
	/**
	 \brief set the bottom edge index
	 \param index the new index to set
	 \return the index of the bottom edge
	 */
	int set_bottom_edge(int index)
	{
		return bottom_edge_index_ = index;
	}
	
	
	/**
	 \brief get the top edge index
	 \return the index of the top edge
	 */
	int top_edge() const
	{
		return top_edge_index_;
	}
	
	
	/**
	 \brief get the top edge index
	 \param index the new index
	 \return the index of the top edge
	 */
	int set_top_edge(int index)
	{
		return top_edge_index_ = index;
	}
	
	
	
	
	
	
	
	/**
	 \brief set the name of the top system
	 \param new_name the new name of the top system
	 */
	void system_name_top(std::string new_name)
	{
		system_name_top_ = new_name;
	}
	
	/**
	 \brief get the name of the top system
	 \return the string name of the top system
	 */
	std::string system_name_top() const
	{
		return system_name_top_;
	}
	
	
	
	/**
	 \brief set the name of the top system
	 \param new_name the new name of the bottom system
	 */
	void system_name_bottom(std::string new_name)
	{
		system_name_bottom_ = new_name;
	}
	
	/**
	 \brief get the name of the bottom system
	 \return the string name of the bottom system
	 */
	std::string system_name_bottom() const
	{
		return system_name_bottom_;
	}
	
	
	
	/**
	 \brief get the index of the crit slice to which this corresponds
	 \return the index of the critical slice to which this Cell corresponds
	 */
	int crit_slice_index() const
	{
		return crit_slice_index_;
	}
	
	
	/**
	 \brief set the index of the crit slice to which this corresponds
	 \param new_index the new index
	 \return the index of the critical slice to which this Cell corresponds
	 */
	int crit_slice_index(int new_index)
	{
		return crit_slice_index_ = new_index;
	}
	
	
	/**
	 \brief get the index of the left edge
	 \return the index
	 \param index the index to look up.
	 */
	int left_edge(unsigned int index) const
	{
		if (index>=left_edges_.size()) {
			throw std::out_of_range("trying to access left edge, index out of range");
			return -1;
		}
		else
		{
			return left_edges_[index];
		}
	}

	
	/**
	 \brief get the index of the right edge
	 \return the index
	 \param index the index to look up.
	 \throws out_of_range, if the index exceeds the size of the right edges
	 */
	int right_edge(unsigned int index) const
	{
		if (index>=right_edges_.size()) {
			throw std::out_of_range("trying to access right edge, index out of range");
			return -1;
		}
		else
		{
			return right_edges_[index];
		}
	}
	
	
	
	
	
	/**
	 \brief get the number of left edges
	 \return the number of left edges
	 */
	unsigned int num_left() const
	{ return left_edges_.size();	 ///< the number of left mapped edges.
	}
	
	/**
	 \brief get the number of right edges
	 \return the number of right edges
	 */
	unsigned int num_right() const
	{ return right_edges_.size();	 ///< the number of right mapped edges.
	}
	
	
	

	
	

	

	
	friend std::ostream & operator<<(std::ostream &os, const Face & f)
	{
		os << f.midpt() << std::endl;
		os << f.crit_slice_index_ << std::endl << f.top_edge_index_ << " " << f.bottom_edge_index_ << std::endl;
		os << f.system_name_top_ << " " << f.system_name_bottom_ << std::endl;
		
		os << f.num_left() << std::endl;
		for (auto jj=f.left_edges_.begin(); jj!=f.left_edges_.end(); jj++) {
			os << *jj << " ";
		}
		os << std::endl;
		
		os << f.num_right() << std::endl;
		for (auto jj=f.right_edges_.begin(); jj!=f.right_edges_.end(); jj++) {
			os << *jj << " ";
		}
		os << std::endl << std::endl;
		
		return os;
	}
	
	
	friend std::istream & operator>>(std::istream &os, Face & f)
	{
		
		f.read_from_stream(os);
		return os;
	}
	
	/**
	 \brief read a Face from an input stream.
	 
	 \param os the istream input stream to pass into this Face.
	 */
	virtual void read_from_stream( std::istream &os )
	{
		
		
		int tmp;
		os >> tmp;  midpt(tmp);
		os >> crit_slice_index_ >> top_edge_index_ >> bottom_edge_index_;
		os >> system_name_top_ >> system_name_bottom_;
		
		unsigned int temp_num_left;
		os >> temp_num_left;
		left_edges_.resize(temp_num_left);
		for (unsigned int jj=0; jj<temp_num_left; jj++) {
			os >> left_edges_[jj];
		}
		

		unsigned temp_num_right;
		os >> temp_num_right;
		right_edges_.resize(temp_num_right);
		for (unsigned int jj=0; jj<temp_num_right; jj++) {
			os >> right_edges_[jj];
		}
	}
	
	
	Face() : Cell()
	{
		init();
	}
	
	~Face(){
		clear_mp(left_crit_val_);
		clear_mp(right_crit_val_);
	}
	
	Face(const Face & other){
		init();
		copy(other);
	}
	
	Face& operator=(const Face & other){

		init();
		
		copy(other);
		return *this;
	}
	
	
	/**
	 \brief single-target MPI send of a Face.
	 
	 
	 \param target The ID of the MPI target to send the Face to.
	 \param mpi_config The current state of MPI
	 */
	void send(int target, ParallelismConfig & mpi_config) const;
	
	
	
	/**
	 \brief single-source MPI receive of a Face.
	 
	 \param source The ID of the MPI source from whom to receive a Face.
	 \param mpi_config The current state of MPI
	 */
	void receive(int source, ParallelismConfig & mpi_config);
	
	/**
	 \brief Test a Face for degeneracy, based on which crit slice the Face is from.
	 
	 \return Boolean, true if degenerate (crit_slice_index < 0), false if not.
	 */
	bool is_degenerate() const
	{
		if (crit_slice_index_<0) {
			return true;
		}
		else{
			return false;
		}
	}
	
	bool is_malformed() const
	{
		return (left_edges_.size()==0 || right_edges_.size()==0);
	}
	
private:
	
	void init()
	{
		system_name_top_ = "UNSET_TOP";
		system_name_bottom_ = "UNSET_BOTTOM";
		
		
		init_mp2(left_crit_val_,1024);
		init_mp2(right_crit_val_,1024);
		
		left_edges_.resize(0);
		right_edges_.resize(0);
		top_edge_index_ = bottom_edge_index_ = -1;
		crit_slice_index_ = -1;

		
	}
	
	
	
	void copy(const Face & other)
	{

				
		Cell::copy(other);
		
		this->system_name_bottom_ = other.system_name_bottom_;
		this->system_name_top_ = other.system_name_top_;
		
		
		this->crit_slice_index_ = other.crit_slice_index_;
		
		this->left_edges_.clear();
		for (auto ii = other.left_edges_.begin(); ii!=other.left_edges_.end(); ii++) {
			this->left_edges_.push_back(*ii);
		}
		
		this->right_edges_.clear();
		for (auto ii = other.right_edges_.begin(); ii!=other.right_edges_.end(); ii++) {
			this->right_edges_.push_back(*ii);
		}
		
		this->top_edge_index_ = other.top_edge_index_;
		this->bottom_edge_index_ = other.bottom_edge_index_;
		
		
		set_mp(this->left_crit_val_, other.left_crit_val_);
		set_mp(this->right_crit_val_, other.right_crit_val_);
	}
	
};





