#ifndef PARAPERM_IMPL_PARAPERM_IMPL_H
#define PARAPERM_IMPL_PARAPERM_IMPL_H

#include <mpi.h>

#include <algorithm>
#include <utility>
#include <vector>

#include <boost/mpi/datatype.hpp>

#include "check_mpi_return_value.h"
#include "random_number_generator.h"

namespace paraperm
{
    template <typename T>
    struct Paraperm<T>::Impl
    {
        typedef std::vector<T> vector_type;

        Impl();
       ~Impl();

        MPI_Comm comm;
        int N, r;

        MPI_Datatype mdt;

        vector_type perm;

        T n, m, pos, count;

        vector_type temp;

        void reset(MPI_Comm comm, T n);

        void phase_1();
        void phase_2();
        void phase_3();

        void finalize();

        void free_communicator();
    };

    template <typename T>
    Paraperm<T>::Impl::Impl()
        : count(0), pos(0), comm(MPI_COMM_NULL)
    {
        mdt = boost::mpi::get_mpi_datatype<T>();
    }

    template <typename T>
    Paraperm<T>::Impl::~Impl()
    {
        free_communicator();
    }

    template <typename T>
    void Paraperm<T>::Impl::reset(MPI_Comm comm, T n)
    {
        int mpi_ret; 

        free_communicator();

        mpi_ret = MPI_Comm_dup(comm, &(this->comm));
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error duplicating MPI communicator!");

        mpi_ret = MPI_Comm_size(this->comm, &N);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error getting the number of MPI processes!");

        mpi_ret = MPI_Comm_rank(this->comm, &r);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error getting MPI process rank!");

        vector_type().swap(perm);
        vector_type().swap(temp);
        this->n = n;

        m = n / N;
        if (n % N)
            m++;

        pos = r * m;
        count = m;
        if ((r + 1) * m > n)
            count = n - pos;
        if (pos >= n)
            count = 0;
    }

    template <typename T>
    void Paraperm<T>::Impl::phase_1()
    {
        int mpi_ret;

        utils::random_number_generator<T> rng(0, N - 1, r);

        std::vector<std::pair<int, T> > pairs(count + 1);

        for (T k = 0; k < count; ++k)
            pairs[k] = std::make_pair(rng(), pos + k);

        // terminator
        pairs[count] = std::make_pair(N, 0);

        std::sort(pairs.begin(), pairs.end());

        // send buffer
        std::vector<T> sendbuf(count, 0);
        for (T k = 0; k < count; ++k)
            sendbuf[k] = pairs[k].second;

        // send counts
        std::vector<int> sendcnts(N, 0);
        T k = 0;
        for (T r_ = 0; r_ < N; r_++) {
            while (r_ == pairs[k].first) {
                sendcnts[r_]++;
                k++;
            }
        }
           
        // no longer needed
        std::vector<std::pair<int, T> >().swap(pairs);

        // send displacements
        std::vector<int> sdispls(N, 0);
        for (T r_ = 1; r_ < N; r_++) 
            sdispls[r_] = sdispls[r_ - 1] + sendcnts[r_ - 1];
           
        // receive counts
        std::vector<int> recvcnts(N);
        mpi_ret = MPI_Alltoall(&(sendcnts[0]), 1, MPI_INT, &(recvcnts[0]), 1, MPI_INT, comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error executing all-to-all MPI operation!");

        // receive displacements
        std::vector<int> rdispls(N, 0);
        for (T r_ = 1; r_ < N; r_++) 
            rdispls[r_] = rdispls[r_ - 1] + recvcnts[r_ - 1];
           
        // total receive count
        T total = 0;
        for (T r_ = 0; r_ < N; r_++) 
            total += recvcnts[r_];

        temp.resize(total);
        mpi_ret = MPI_Alltoallv(
                &(sendbuf[0]), &(sendcnts[0]), &(sdispls[0]), mdt,
                &(temp[0]),    &(recvcnts[0]), &(rdispls[0]), mdt, comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error executing all-to-all MPI operation!");

        mpi_ret = MPI_Barrier(comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error executing MPI barrier!");
    }

    template <typename T>
    void Paraperm<T>::Impl::phase_2()
    {
        utils::random_number_generator<T> rng(r);

        if (temp.size() > 1) {
            for (T k = temp.size() - 1; k > 0; --k) {
                T l = rng() % (k + 1);
                std::swap(temp[k], temp[l]);
            }
        }

        int mpi_ret = MPI_Barrier(comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error executing MPI barrier!");
    }

    template <typename T>
    void Paraperm<T>::Impl::phase_3()
    {
        int mpi_ret; 

        T size = temp.size();
        T first;

        mpi_ret = MPI_Scan(&size, &first, 1, mdt, MPI_SUM, comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error counting parallel prefix sum via MPI!");

        first -= size;

        T last = first + size - 1;
        int r_ = first / m;
        T first_ = first;
        T remains = count;

        perm.resize(m);

        std::vector<MPI_Request> requests;

        enum Tag { ELEMENTS = 1, POSITION = 2 };

        while (true) {
            T last_ = (r_ + 1) * m - 1;
            last_ = std::min(last_, last);
            T count_ = last_ - first_ + 1;

            if (r_ == r) {
                for (T k = first_; k <= last_; ++k)
                    perm[k - pos] = temp[k - first];
                remains -= count_;
            }
            else {
                MPI_Request request;

                T buf[2] = { first_, count_ };
                mpi_ret = MPI_Isend(buf, 2, mdt, r_, POSITION, comm, &request);
                CHECK_MPI_RETURN_VALUE(mpi_ret, "Error sending data via MPI!");
                requests.push_back(request);

                mpi_ret = MPI_Isend(&(temp[first_ - first]), count_, mdt, r_, ELEMENTS, comm, &request);
                CHECK_MPI_RETURN_VALUE(mpi_ret, "Error sending data via MPI!");
                requests.push_back(request);
            }

            r_++;
            first_ += count_;

            if (first_ > last)
                break;
        }

        while (remains > 0) {
            T buf[2];
            MPI_Status status;
            mpi_ret = MPI_Recv(buf, 2, mdt, MPI_ANY_SOURCE, POSITION, comm, &status);
            CHECK_MPI_RETURN_VALUE(mpi_ret, "Error receiving data via MPI!");

            T first_ = buf[0];
            T count_ = buf[1];

            mpi_ret = MPI_Recv(&(perm[first_ - pos]), count_, mdt, status.MPI_SOURCE, ELEMENTS, comm, &status);
            CHECK_MPI_RETURN_VALUE(mpi_ret, "Error receiving data via MPI!");

            remains -= count_;
        }

        for (int k = 0; k < requests.size(); ++k) {
            MPI_Status status;
            mpi_ret = MPI_Wait(&(requests[k]), &status);
            CHECK_MPI_RETURN_VALUE(mpi_ret, "Error waiting for MPI request comletion!");
        }

        vector_type().swap(temp);

        mpi_ret = MPI_Barrier(comm);
        CHECK_MPI_RETURN_VALUE(mpi_ret, "Error executing MPI barrier!");
    }

    template <typename T>
    void Paraperm<T>::Impl::finalize()
    {
        if (perm.size() > count)
            perm.resize(count);

        free_communicator();
    }

    template <typename T>
    void Paraperm<T>::Impl::free_communicator()
    {
        if (MPI_COMM_NULL != comm) {
            int mpi_ret = MPI_Comm_free(&comm);
            CHECK_MPI_RETURN_VALUE(mpi_ret, "Error freeing MPI communicator!");
            comm = MPI_COMM_NULL;
        }
    }
}

#endif
