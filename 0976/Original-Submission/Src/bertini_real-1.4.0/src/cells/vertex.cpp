#include "cells/vertex.hpp"

void Vertex::set_point(const vec_mp new_point)
{
	change_prec_vec_mp(this->pt_mp_, new_point->curr_prec);
	vec_cp_mp(this->pt_mp_, new_point);
}




void Vertex::send(int target, ParallelismConfig & mpi_config) const
{
	
	send_vec_mp(pt_mp_, target);
	
	send_vec_mp(projection_values_, target);
	
	int * buffer = (int *) br_malloc(2*sizeof(int));
	buffer[0] = type_;
	buffer[1] = input_filename_index_;
	
	MPI_Send(buffer, 2, MPI_INT, target, VERTEX, mpi_config.comm());
	free(buffer);
	
}


void Vertex::receive(int source, ParallelismConfig & mpi_config)
{
	MPI_Status statty_mc_gatty;
	int * buffer = (int *) br_malloc(2*sizeof(int));
	
	
	receive_vec_mp(pt_mp_, source);
	receive_vec_mp(projection_values_, source);
	
	MPI_Recv(buffer, 2, MPI_INT, source, VERTEX, mpi_config.comm(), &statty_mc_gatty);
	
	type_ = static_cast<VertexType>(buffer[0]);
	input_filename_index_ = buffer[1];

	free(buffer);
}

