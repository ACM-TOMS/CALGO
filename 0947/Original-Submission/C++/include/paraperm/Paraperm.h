#ifndef PARAPERM_PARAPERM_H
#define PARAPERM_PARAPERM_H

#include <mpi.h>
#include <cstdint>
#include <vector>
#include <boost/noncopyable.hpp>

namespace paraperm
{
    template <typename T = uintmax_t>
    class Paraperm : boost::noncopyable
    {
        public:
            typedef T value_type;
            typedef std::vector<T> vector_type;

            Paraperm();
           ~Paraperm();

            void generate(MPI_Comm comm, T n);

            const vector_type& perm() const;
            T pos() const;
            T count() const;

        private:
            struct Impl;
            Impl* pimpl_;
    };
}

#include "impl/Paraperm.h"

#endif
