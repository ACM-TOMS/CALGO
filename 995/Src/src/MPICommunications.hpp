#include <mpi.h>
#include <stdexcept>

inline void MPICommunications::initialize() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if(not initialized)
        MPI_Init(NULL, NULL);
}

inline void MPICommunications::finalize() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if(not finalized)
        MPI_Finalize();
}

inline int MPICommunications::myRank() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

inline int MPICommunications::numberOfProcesses() {
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
}

inline MPI_Datatype MPICommunications::getType(int value) {
    return MPI_INT;
}

inline MPI_Datatype MPICommunications::getType(size_t value) {
    return MPI_LONG;
}

inline MPI_Datatype MPICommunications::getType(long value) {
    return MPI_LONG;
}

inline MPI_Datatype MPICommunications::getType(double value) {
    return MPI_DOUBLE;
}

template<typename T>
void MPICommunications::Send(std::vector<T>& send_buffer, int size, int destination) {
    int error = MPI_Send(send_buffer.data(), size*sizeof(T), MPI_CHAR, destination, 0, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Send ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template <typename T>
void MPICommunications::Send(MPI_Datatype datatype, T& value, int destination) {
    int error = MPI_Send(&value, 1, datatype, destination, 0, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Send ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Recv(std::vector<T>& recv_buffer, int size, int source) {
    recv_buffer.clear();
    recv_buffer.resize(size);
    int error = MPI_Recv(recv_buffer.data(), size*sizeof(T), MPI_CHAR, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Recv ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template <typename T>
void MPICommunications::Recv(MPI_Datatype datatype, T& value, int source) {
    int error = MPI_Recv(&value, 1, datatype, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Recv ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Gather(T value, std::vector<T>& recv_buffer, int root) {
    recv_buffer.resize(numberOfProcesses());
    int error = MPI_Gather(&value, 1, getType(value), recv_buffer.data(), 1, getType(value), root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Gather ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Gatherv(MPI_Datatype datatype, const std::vector<T>& send_buffer, std::vector<T>& recv_buffer,
                                int root) {
    int send_count = static_cast<int>(send_buffer.size());
    int processes = numberOfProcesses();
    std::vector<int> recv_counts(processes, 0);
    Gather(send_count, recv_counts, root);

    std::vector<int> map;
    if(myRank() == root) {
        map.assign(processes+1, 0);
        for(int i=1; i<processes+1; ++i)
            map[i] = map[i-1] + recv_counts[i-1];
        recv_buffer.resize(map.back());
    }

    int error = MPI_Gatherv(send_buffer.data(), send_count, datatype, recv_buffer.data(), recv_counts.data(),
                            map.data(), datatype, root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Gatherv ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Gatherv(const std::vector<T>& send_buffer, std::vector<T>& recv_buffer, int root) {
    int send_count = static_cast<int>(send_buffer.size()*sizeof(T));
    int processes = numberOfProcesses();
    std::vector<int> recv_counts(processes, 0);
    Gather(send_count, recv_counts, root);

    std::vector<int> map;
    if(myRank() == root) {
        map.assign(processes+1, 0);
        for(int i=1; i<processes+1; ++i)
            map[i] = map[i-1] + recv_counts[i-1];
        recv_buffer.resize(map.back()/sizeof(T));
    }

    int error = MPI_Gatherv(send_buffer.data(), send_count, MPI_CHAR, recv_buffer.data(), recv_counts.data(),
                            map.data(), MPI_CHAR, root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Gatherv ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Scatter(std::vector<T>& send_buffer, T& recv_value, int root) {
    int error = MPI_Scatter(send_buffer.data(), 1, getType(recv_value), &recv_value, 1, getType(recv_value), root,
                            MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Scatter ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Scatterv(std::vector<T>& send_buffer, std::vector<T>& recv_buffer, int root) {
    std::vector<int> send_counts;
    std::vector<int> displacements;
    int local_size = 0;

    if(myRank() == root) {
        int total_size = static_cast<int>(send_buffer.size());
        int processes = numberOfProcesses();
        send_counts.assign(processes, total_size / processes);
        for (int i=0; i<processes; ++i)
            if(i < (total_size % processes))
                ++send_counts[i];
        displacements.assign(processes, 0);
        for (int i=1; i<processes; ++i)
            displacements[i] = displacements[i-1] + send_counts[i-1];
    }

    Scatter(send_counts, local_size, root);
    recv_buffer.resize(local_size);
    for(int& s : send_counts)
        s *= sizeof(T);
    for(int& d : displacements)
        d *= sizeof(T);
    local_size *= sizeof(T);

    int error = MPI_Scatterv(send_buffer.data(), send_counts.data(), displacements.data(), MPI_CHAR, recv_buffer.data(),
                             local_size, MPI_CHAR, root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Scatterv ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Broadcast(T& value, int root) {
    int error = MPI_Bcast(&value, 1, getType(value), root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Broadcast ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Broadcast(std::vector<T>& buffer, int buffer_size, int root) {
    if(static_cast<int>(buffer.size()) != buffer_size)
        buffer.resize(buffer_size);

    int error = MPI_Bcast(buffer.data(), buffer_size, getType(T()), root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Broadcast ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}

template<typename T>
void MPICommunications::Broadcast(std::vector<T>& buffer, int root) {
    int size = 0;
    if(myRank() == root)
        size = static_cast<int>(buffer.size());
    Broadcast(size, root);

    if(myRank() != root) {
        buffer.clear();
        buffer.resize(size);
    }

    int error = MPI_Bcast(buffer.data(), size*sizeof(T), MPI_CHAR, root, MPI_COMM_WORLD);
    if(MPI_SUCCESS != error)
        throw std::logic_error("MPI Broadcast ERROR " + std::to_string(error) + " on Rank " + std::to_string(myRank()));
}
