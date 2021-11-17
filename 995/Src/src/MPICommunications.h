#ifndef FEA_MESHER2D_MPICOMMUNICATIONS_H
#define FEA_MESHER2D_MPICOMMUNICATIONS_H

#include <mpi.h>
#include <vector>

namespace MPICommunications {
    //Wrapper for MPI_Init
    void initialize();

    //Wrapper for MPI_Finalize
    void finalize();

    int myRank();
    int numberOfProcesses();

    //Helper functions for determining the MPI_Datatype for the templated functions
    MPI_Datatype getType(int value);
    MPI_Datatype getType(size_t value);
    MPI_Datatype getType(long value);
    MPI_Datatype getType(double value);

    //The following function are wrappers for their corresponding MPI call
    //For example, Gatherv performs an MPI_Gatherv operation

    template<typename T>
    void Send(std::vector<T>& send_buffer, int size, int destination);

    template <typename T>
    void Send(MPI_Datatype datatype, T& value, int destination);

    template<typename T>
    void Recv(std::vector<T>& recv_buffer, int size, int source);

    template <typename T>
    void Recv(MPI_Datatype datatype, T& value, int source);

    template<typename T>
    void Scatter(std::vector<T>& send_buffer, T& recv_value, int root);

    template<typename T>
    void Scatterv(std::vector<T>& send_buffer, std::vector<T>& recv_buffer, int root);

    template<typename T>
    void Gather(T value, std::vector<T>& recv_buffer, int root);

    template<typename T>
    void Gatherv(MPI_Datatype datatype, const std::vector<T>& send_buffer, std::vector<T>& recv_buffer, int root);

    template<typename T>
    void Gatherv(const std::vector<T>& send_buffer, std::vector<T>& recv_buffer, int root);

    template<typename T>
    void Broadcast(T& value, int root);

    template<typename T>
    void Broadcast(std::vector<T>& buffer, int buffer_size, int root);

    template<typename T>
    void Broadcast(std::vector<T>& buffer, int root);
} //MPICommunications namespace

#include "MPICommunications.hpp"


#endif //FEA_MESHER2D_MPICOMMUNICATIONS_H
