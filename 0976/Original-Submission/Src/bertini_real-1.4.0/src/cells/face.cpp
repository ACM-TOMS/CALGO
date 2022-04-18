#include "cells/face.hpp"




void Face::send(int target, ParallelismConfig & mpi_config) const
{
#ifdef functionentry_output
	std::cout << "Face::send" << std::endl;
#endif
	
	
	
	
	Cell::send(target,mpi_config);
	
	int * buffer = new int[7];
	
	buffer[0] = num_left();
	buffer[1] = num_right();
	buffer[2] = top_edge_index_;
	buffer[3] = bottom_edge_index_;
	buffer[4] = system_name_bottom_.size();
	buffer[5] = system_name_top_.size();
	buffer[6] = crit_slice_index_;
	
	MPI_Send(buffer, 7, MPI_INT, target, FACE, mpi_config.comm());
	delete[] buffer;
	
	
	if (num_left()>0) {
		buffer = new int[num_left()];
		for (unsigned int ii=0; ii<num_left(); ii++) {
			buffer[ii] = left_edges_[ii];
		}
		MPI_Send(buffer, num_left(), MPI_INT, target, FACE, mpi_config.comm());
		delete[] buffer;
	}
	
	
	
	if (num_right()>0) {
		buffer = new int[num_right()];
		
		for (unsigned int ii=0; ii<num_right(); ii++) {
			buffer[ii] = right_edges_[ii];
		}
		MPI_Send(buffer, num_right(), MPI_INT, target, FACE, mpi_config.comm());
		delete[] buffer;
	}
	
	
	
	
	
	
	
	std::string sendme = system_name_bottom_;
	sendme.append(system_name_top_);
	
	int num_to_send = sendme.size()+1;
	char * charbuff = new char[num_to_send];
	strcpy(charbuff, sendme.c_str());
	charbuff[num_to_send-1] = '\0';
	MPI_Send(&charbuff[0], num_to_send, MPI_CHAR, target, FACE, mpi_config.comm());
	delete [] charbuff;
	
	
	send_comp_mp(left_crit_val_, target);
	send_comp_mp(right_crit_val_, target);
	
	
	return;
	
}

void Face::receive(int source, ParallelismConfig & mpi_config)
{
#ifdef functionentry_output
	std::cout << "Face::receive" << std::endl;
#endif
	
	

	Cell::receive(source,mpi_config);
	
	MPI_Status statty_mc_gatty;
	int * buffer= new int[7];
	
	MPI_Recv(buffer, 7, MPI_INT, source, FACE, mpi_config.comm(), &statty_mc_gatty);
	
	int tmp_size_left = buffer[0];
	int tmp_size_right = buffer[1];
	top_edge_index_ = buffer[2];
	bottom_edge_index_ = buffer[3];
	int nchars_name_bottom = buffer[4];
	int nchars_name_top = buffer[5];
	crit_slice_index_ = buffer[6];
	
	delete[] buffer;
	
	
	if (tmp_size_left>0) {
		int * buffer2 = new int[tmp_size_left];
		MPI_Recv(buffer2, tmp_size_left, MPI_INT, source, FACE, mpi_config.comm(), &statty_mc_gatty);
		for (int ii=0; ii<tmp_size_left; ii++) {
			add_left_edge(buffer2[ii]);
		}
		delete[] buffer2;
	}
	
	
	
	if (tmp_size_right>0) {
		int * buffer3 = new int[tmp_size_right];
		MPI_Recv(buffer3, tmp_size_right, MPI_INT, source, FACE, mpi_config.comm(), &statty_mc_gatty);
		for (int ii=0; ii<tmp_size_right; ii++) {
			add_right_edge(buffer3[ii]);
		}
		delete[] buffer3;
	}
	
	
	
	
	
	
	
	
	
	char * charbuff = new char[nchars_name_bottom+nchars_name_top+1];
	
	MPI_Recv(charbuff, nchars_name_bottom+nchars_name_top+1, MPI_CHAR, source, FACE, mpi_config.comm(), &statty_mc_gatty);
	
	std::stringstream converter;
	for (int jj=0; jj<nchars_name_bottom; ++jj) {
		converter << charbuff[jj];
	}
	system_name_bottom_ = converter.str();
	converter.clear();
	converter.str("");
	
	int offset = nchars_name_bottom;
	for (int jj=0; jj<nchars_name_top; ++jj) {
		converter << charbuff[jj+offset];
	}
	system_name_top_ = converter.str();
	
	delete [] charbuff;
	
	
	receive_comp_mp(left_crit_val_,source);
	receive_comp_mp(right_crit_val_,source);
	
	return;
	
}






