#if CUDA_ENABLED

#include "cuda_utils.h"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>

#include "tensor.h"

namespace cuda {
void identify_devices() {
  int n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  std::cout << "Number of devices: " << n_devices << std::endl;

  for (auto id = 0; id < n_devices; id++) {
    cudaDeviceProp prop{};
    cudaGetDeviceProperties(&prop, id);
    printf("Device Number: %d\n", id);
    printf("  Device name: %s\n", prop.name);
    printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
    printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
    printf("  Peak Memory Bandwidth (GB/s): %f\n\n", 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
  }
}

void allocate(double *&d_A, long int size) {
  if (cudaMalloc(reinterpret_cast<void **>(&d_A), size * sizeof(double)) != cudaSuccess)
    std::cerr << "Device memory allocation error" << std::endl;
}

void deallocate(double *&d_A) {
  auto error = cudaFree(d_A);
  if (error != cudaSuccess)
    std::cerr << "Memory free error: " << error << std::endl;
}

void allocate_async(double *&host_ptr, long int size) {
  auto error = cudaHostAlloc((void **)&host_ptr, size * sizeof(double), cudaHostAllocDefault);

  if (error != cudaSuccess) {
    std::cerr << "Cannot allocate memory on host: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}

void deallocate_async(double *&host_ptr) {
  auto error = cudaFreeHost(host_ptr);

  if (error != cudaSuccess) {
    std::cerr << "Cannot allocate memory on host: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}

void send_to_device(double *h_A, double *d_A, int size) {
  auto status = cublasSetVector(size, sizeof(double), h_A, 1, d_A, 1);

  if (status != CUBLAS_STATUS_SUCCESS) {
    std::cerr << "Data transfer error to Device. Return value: " << status << std::endl;
    exit(EXIT_FAILURE);
  }
}

void send_to_device_async(double *h_A, double *d_A, int size, cudaStream_t &stream) {
  auto error = cudaMemcpyAsync(d_A, h_A, size * sizeof(double), cudaMemcpyHostToDevice, stream);
  if (error != cudaSuccess) {
    std::cerr << "Data transfer error to Device: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}

void receive_from_device(double *d_A, double *h_A, int size) {
  auto status = cublasGetVector(size, sizeof(double), d_A, 1, h_A, 1);

  if (status != CUBLAS_STATUS_SUCCESS) {
    std::cerr << "Data transfer error from Device: " << status << std::endl;
    exit(EXIT_FAILURE);
  }
}

void receive_from_device_async(double *d_A, double *h_A, int size, cudaStream_t &stream) {
  auto error = cudaMemcpyAsync(h_A, d_A, size * sizeof(double), cudaMemcpyDeviceToHost, stream);

  if (error != cudaSuccess) {
    std::cerr << "Data transfer error from Device: " << error << std::endl;
    exit(EXIT_FAILURE);
  }
}

void initialize_stream(cudaStream_t &stream) {
  auto status = cudaStreamCreate(&stream);

  if (status != cudaSuccess)
    std::cerr << "Stream creation failed." << std::endl;
}

cudaStream_t create_stream() {
  cudaStream_t stream;
  auto status = cudaStreamCreate(&stream);

  if (status != cudaSuccess)
    std::cerr << "Stream creation failed." << std::endl;

  return stream;
}

void destroy_stream(cudaStream_t &stream) {
  auto status = cudaStreamDestroy(stream);

  if (status != cudaSuccess)
    std::cerr << "Stream creation failed." << std::endl;
}
} // namespace cuda

namespace cuda::cublas {
cublasHandle_t create_handle() {
  cublasHandle_t handle;
  auto status = cublasCreate_v2(&handle);

  if (status != CUBLAS_STATUS_SUCCESS)
    std::cerr << "CUBLAS initialization error" << std::endl;

  return handle;
}

void destroy_handle(cublasHandle_t &handle) {
  auto status = cublasDestroy(handle);

  if (status != CUBLAS_STATUS_SUCCESS)
    std::cerr << "CUBLAS shutdown error";
}
} // namespace cuda::cublas

namespace cuda {}

int test_func() {
  double *a = nullptr;

  auto handle = cuda::cublas::create_handle();

  auto t = cals::Tensor(5, {10, 10, 10});
  t.randomize();
  t.print();

  std::cout << "T n elements: " << t.get_n_elements() << std::endl;

  cuda::allocate(a, t.get_n_elements());

  cuda::send_to_device(t.get_data(), a, t.get_n_elements());

  double zero = 0.0;
  cublasDscal(handle, t.get_n_elements(), &zero, a, 1);

  cuda::receive_from_device(a, t.get_data(), t.get_n_elements());

  t.print();

  cuda::deallocate(a);

  cuda::cublas::destroy_handle(handle);

  std::cout << "Alles gut" << std::endl;
  return EXIT_SUCCESS;
}
#endif
