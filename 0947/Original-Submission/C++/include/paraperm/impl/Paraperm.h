#ifndef PARAPERM_IMPL_PARAPERM_H
#define PARAPERM_IMPL_PARAPERM_H

#include <mpi.h>

#include "Paraperm_Impl.h"

namespace paraperm
{
    template <typename T>
    Paraperm<T>::Paraperm() 
        : pimpl_(new Impl)
    {
    }

    template <typename T>
    Paraperm<T>::~Paraperm()
    {
        delete pimpl_;
    }

    template <typename T>
    void Paraperm<T>::generate(MPI_Comm comm, T n)
    {
        int mpi_ret;

        pimpl_->reset(comm, n);
        pimpl_->phase_1();
        pimpl_->phase_2();
        pimpl_->phase_3();
        pimpl_->finalize();
    }

    template <typename T>
    const typename Paraperm<T>::vector_type& Paraperm<T>::perm() const
    {
        return pimpl_->perm;
    }

    template <typename T>
    T Paraperm<T>::pos() const
    {
        return pimpl_->pos;
    }

    template <typename T>
    T Paraperm<T>::count() const
    {
        return pimpl_->count;
    }

}

#endif
