#ifndef CALS_CUDA_UTILS_H
#define CALS_CUDA_UTILS_H

#if CUDA_ENABLED

#include <cuda.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

int const cuda_n_streams = 2;

namespace cuda
{
  void identify_devices();

  cudaStream_t create_stream();

  void initialize_stream(cudaStream_t &stream);

  void destroy_stream(cudaStream_t &stream);

  void allocate(double* &d_A, long int size);

  void deallocate(double* &d_A);

  void allocate_async(double* &host_ptr, long int size);

  void deallocate_async(double* &host_ptr);

  void send_to_device(double* h_A, double* d_A, int size);

  void send_to_device_async(double* h_A, double* d_A, int size, cudaStream_t &stream);

  void receive_from_device(double* d_A, double* h_A, int size);

  void receive_from_device_async(double* d_A, double* h_A, int size, cudaStream_t &stream);
}

namespace cuda::cublas
{
  cublasHandle_t create_handle();

  void destroy_handle(cublasHandle_t &handle);
}

namespace cuda
{
  struct CudaParams
  {
    cublasHandle_t handle;
    cudaStream_t streams[cuda_n_streams]{};
    double const one{1.0};
    double const zero{0.0};

    CudaParams()
    {
      handle = cuda::cublas::create_handle();
      for (auto & stream : streams)
        initialize_stream(stream);
    }

    ~CudaParams()
    {
      cuda::cublas::destroy_handle(handle);
      for (auto & stream : streams)
        destroy_stream(stream);
    }
  };
}
#endif

#endif //CALS_CUDA_UTILS_H
