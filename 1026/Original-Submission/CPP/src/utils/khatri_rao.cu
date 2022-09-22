#ifndef CP_CALS_KHATRI_RAO_CU
#define CP_CALS_KHATRI_RAO_CU

// #include "stdio.h"

namespace cals
{
  namespace mttkrp
  {
    __global__ void khatri_rao_cuda_kernel(double const *const A,
                                           double const *const B,
                                           unsigned int IA,
                                           unsigned int IB,
                                           unsigned int JK,
                                           double *const K)
    {
      const unsigned int n_threads = blockDim.x * gridDim.x;
      const unsigned int IAB = IA * IB;
      const unsigned int n_elements = IAB * JK;

      //const unsigned int work_per_thread = ceil((float)n_elements / (float)n_threads);
      //const unsigned int thread_idx = blockIdx.x * blockDim.x * work_per_thread + threadIdx.x * work_per_thread;
      const unsigned int thread_idx = blockIdx.x * blockDim.x + threadIdx.x;

      //for (auto i = thread_idx; i < thread_idx + work_per_thread; i ++)
      for (auto i = thread_idx; i < n_elements; i += n_threads)
      {
        //if (i < n_elements)
        //{
        const unsigned int j = i / IAB;
        const unsigned modulo = i % IAB;
        const unsigned int iA = modulo / IB;
        const unsigned int iB = modulo % IB;

        K[i] = A[iA + j * IA] * B[iB + j * IB];
        //}
      }
    }

    void khatri_rao_cuda(double *A,
                         double *B,
                         unsigned int IA,
                         unsigned int IB,
                         unsigned int JK,
                         double *K)
    {
      dim3 block_dim;
      dim3 grid_dim;
      cudaOccupancyMaxPotentialBlockSize(reinterpret_cast<int *>(&grid_dim),
                                         reinterpret_cast<int *>(&block_dim),
                                         khatri_rao_cuda_kernel);
      //printf("Blocks: %d %d %d \n", block_dim.x, block_dim.y, block_dim.z);
      //printf("Grids : %d %d %d \n", grid_dim.x, grid_dim.y, grid_dim.z);
      khatri_rao_cuda_kernel<<<grid_dim, block_dim>>>(A, B, IA, IB, JK, K);
    }
  }
}
#endif //CP_CALS_KHATRI_RAO_CU
