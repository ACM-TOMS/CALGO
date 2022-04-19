#include <mpi.h>
#include <cstdint>
#include <paraperm/Paraperm.h>

typedef paraperm::Paraperm<> Paraperm;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);

    Paraperm paraperm;

    const Paraperm::value_type n = (1UL << 24) * N; 
    paraperm.generate(MPI_COMM_WORLD, n);

    const Paraperm::vector_type& perm = paraperm.perm();
    const Paraperm::value_type pos = paraperm.pos();
    const Paraperm::value_type count = paraperm.count();

    // do whatever with the generated permutation

    MPI_Finalize();
    return 0;
}
