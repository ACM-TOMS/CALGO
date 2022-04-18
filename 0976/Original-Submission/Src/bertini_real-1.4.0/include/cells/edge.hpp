#pragma once

/** \file edge.hpp */

#include "cell.hpp"

/**
 \brief 1-cell.
 
 the edge data type.  has three indices: left, right, midpt.
 */
class Edge : public Cell
{
	int left_;  ///< index into vertices
	int right_; ///< index into vertices
	
	std::vector< int > removed_points_;
	
	
public:
	
	
	typedef std::vector< int >::iterator removed_iterator;
	typedef std::vector< int >::const_iterator removed_const_iterator;
	
	/**
	 \return get an iterator to the beginning of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_iterator removed_begin() {return removed_points_.begin();}
	
	/**
	 \return get an iterator to the beginning of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_const_iterator removed_begin() const {return removed_points_.begin();}
	
	/**
	 \return get an iterator to the end of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_iterator removed_end() {return removed_points_.end();}
	
	/**
	 \return get an iterator to the end of the set of removed points, which were merged away in the merge operation.
	 */
	inline removed_const_iterator removed_end() const {return removed_points_.end();}
	
	
	
	
	/**
	 \brief adds a point as a removed point.  tacks on to the end of the vector
	 
	 \param new_removed_point the index of the point to add
	 \return the index of the point
	 */
	int add_removed_point(int new_removed_point)
	{
		removed_points_.push_back(new_removed_point);
		return new_removed_point;
	}
	
	
	/**
	 \brief get the right point
	 
	 \return the index of the right point
	 */
	inline int right() const
	{
		return right_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_right the new index of the left point
	 \return the index of the left point
	 */
	int right(int new_right)
	{
		return right_ = new_right;
	}
	
	
	
	
	
	/**
	 \brief get the left point
	 
	 \return the index of the left point
	 */
	inline int left() const
	{
		return left_;
	}
	
	
	/**
	 \brief set the left point
	 \param new_left the new index of the left point
	 \return the index of the left point
	 */
	int left(int new_left)
	{
		return left_ = new_left;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 default constructor
	 */
	Edge() : Cell()
	{
		left_ = right_ = -1;
	}
	
	
	/**
	 \brief construct edge from left mid and right indices
	 \param new_left the new left index for constructed edge
	 \param new_midpt the new mid index for constructed edge
	 \param new_right the new right index for constructed edge
	 */
	Edge(int new_left, int new_midpt, int new_right)
	{
		left(new_left);
		right(new_right);
		midpt(new_midpt);
	}
	
	
	
	
	/**
	 check whether the edge is degenerate
	 \return true if the edge is degenerate, false if not.
	 */
	inline bool is_degenerate() const
	{
		if ((left() == right()) && (left()==midpt()) && (right()==midpt()))
			return true;
		else
			return false;
	}
	
	
	
	/**
	 \brief send to a single target
	 \param target to whom to send this edge
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, ParallelismConfig & mpi_config) const
	{
		int * buffer = new int[4];
		
		buffer[0] = left();
		buffer[1] = midpt();
		buffer[2] = right();
		buffer[3] = removed_points_.size();
		
		MPI_Send(buffer, 4, MPI_INT, target, EDGE, mpi_config.comm());
		
		delete [] buffer;
		
		buffer = new int[removed_points_.size()];
		for (unsigned int ii=0; ii!=removed_points_.size(); ii++) {
			buffer[ii] = removed_points_[ii];
		}
		MPI_Send(buffer, removed_points_.size(), MPI_INT, target, EDGE, mpi_config.comm());
		delete [] buffer;
		
		
	}
	
	
	/**
	 \brief receive from a single source
	 \param source from whom to receive this edge
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, ParallelismConfig & mpi_config) 
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[4];
		MPI_Recv(buffer, 4, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		
		left(buffer[0]);
		midpt(buffer[1]);
		right(buffer[2]);
		int temp_num_removed = buffer[3];
		
		
		delete [] buffer;
		
		buffer = new int[temp_num_removed];
		MPI_Recv(buffer, temp_num_removed, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		for (int ii=0; ii<temp_num_removed; ii++) {
			removed_points_.push_back(buffer[ii]);
		}
		
		delete [] buffer;
		
	}
	
	
	/**
	 \brief get edge from input stream.  this function is defunct, and needs implementation apparently.
	 \param is the stream from whom to read
	 */
	virtual void read_from_stream( std::istream &is )
	{
		
		is >> left_ >> midpt_ >> right_;
	}
	
	friend std::ostream & operator<<(std::ostream & out, const Edge & E){
		out << E.left_ << " " << E.midpt_ << " " << E.right_;
		return out;
	}
	
};



/**
\brief Holds metadata for an edge.  Namely, the cycle numbers.

Edges produce cycle numbers as you track from the generic midpoint out to either boundary point.  This metadata class holds this information.
*/
class EdgeCycleNumbers
{
private:
	int cycle_number_left_ = -1;
	int cycle_number_right_ = -1;
public:
	EdgeCycleNumbers(int ell, int arr) : cycle_number_left_(ell), cycle_number_right_(arr)
	{	}

	EdgeCycleNumbers() = default;

	void CycleNumLeft(int c)
	{
		cycle_number_left_ = c;
	}

	void CycleNumRight(int c)
	{
		cycle_number_right_ = c;
	}

	int CycleNumLeft() const
	{
		return cycle_number_left_;
	}

	int CycleNumRight() const
	{
		return cycle_number_right_;
	}



	/**
	 \brief send to a single target
	 \param target to whom to send this edge
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, ParallelismConfig & mpi_config) const
	{
		int * buffer = new int[2];
		
		buffer[0] = CycleNumLeft();
		buffer[1] = CycleNumRight();
		
		MPI_Send(buffer, 2, MPI_INT, target, EDGE, mpi_config.comm());
		
		delete [] buffer;
	}
	
	
	/**
	 \brief receive from a single source
	 \param source from whom to receive this edge
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, ParallelismConfig & mpi_config) 
	{
		MPI_Status statty_mc_gatty;
		int * buffer = new int[2];
		MPI_Recv(buffer, 2, MPI_INT, source, EDGE, mpi_config.comm(), &statty_mc_gatty);
		
		CycleNumLeft(buffer[0]);
		CycleNumRight(buffer[1]);

		delete [] buffer;		
	}

	friend std::ostream& operator<<(std::ostream& out, EdgeCycleNumbers const& e)
	{
		out << e.CycleNumLeft() << " " << e.CycleNumRight();
		return out;
	}
};




