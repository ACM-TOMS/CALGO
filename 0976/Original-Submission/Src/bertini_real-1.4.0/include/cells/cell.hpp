#pragma once

/** \file cell.hpp */

#include "programConfiguration.hpp"

/**
 \brief The base cell class.  All n-cells should inherit from this.
 
 Contains a midpoint, methods to send and receive, etc.
 */
class Cell
{
	
protected:

	int midpt_; ///< index into Vertex set
	
	
public:
	
	/**
	 \brief get the midpoint
	 
	 \return the index of the midpoint
	 */
	inline int midpt() const
	{
		return midpt_;
	}
	
	/**
	 \brief set the midpoint
	 
	 \param new_mid the new index of the midpoint
	 \return the new index of the midpoint
	 */
	int midpt(int new_mid)
	{
		return midpt_ = new_mid;
	}
	
	
	
	
	/**
	 \brief get a cell from an input stream
	 \return reference to the input stream, so you can chain inputs together.
	 \param is the input stream to get from.
	 \param c Cell to get from stream.
	 */
	friend std::istream & operator>>(std::istream &is,  Cell & c)
	{
		is >> c.midpt_;
		return is;
	}
	
	/**
	 \brief copy to another cell
	 \param other the other cell to copy to
	 */
	inline void copy(const Cell & other){
		midpt(other.midpt());
	}
	
	/**
	 \brief send cell to target in communicator
	 \param target the integer id of the target of the send
	 \param mpi_config the current state of parallelism
	 */
	void send(int target, ParallelismConfig & mpi_config) const
	{
		int buffer = midpt();
		MPI_Send(&buffer, 1, MPI_INT, target, CELL, mpi_config.comm());
	}
	
	
	/**
	 \brief receive cell from source in communicator
	 \param source the integer id of the source of the send
	 \param mpi_config the current state of parallelism
	 */
	void receive(int source, ParallelismConfig & mpi_config)
	{
		MPI_Status statty_mc_gatty;
		int buffer;
		MPI_Recv(&buffer, 1, MPI_INT, source, CELL, mpi_config.comm(), &statty_mc_gatty);
		midpt(buffer);
	}
	
	/**
	 \brief virtual read_from_stream function which enforces all cells to have this functtion.
	 \param is input stream from whom to read
	 */
	virtual void read_from_stream( std::istream &is ) = 0;
	
};






