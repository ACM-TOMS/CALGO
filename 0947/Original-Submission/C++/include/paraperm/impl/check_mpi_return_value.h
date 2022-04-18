#ifndef PARAPERM_IMPL_CHECK_MPI_RETURN_VALUE_H
#define PARAPERM_IMPL_CHECK_MPI_RETURN_VALUE_H

#define CHECK_MPI_RETURN_VALUE(return_value, error_message)       \
    do {                                             \
        if (MPI_SUCCESS != return_value)             \
            throw std::runtime_error(error_message); \
    } while (0)

#endif
