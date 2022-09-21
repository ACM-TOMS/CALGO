#include "utils/mttkrp.h"

#include <stack>
#include <fstream>
#include <iostream>

#include "utils/utils.h"

using std::vector;
using std::string;
using cals::Matrix;

using MTTKRP_TIMERS = cals::MttkrpTimers::TIMERS;

namespace cals::mttkrp
{

  void khatri_rao_cuda(double *A,
                       double *B,
                       unsigned int IA,
                       unsigned int IB,
                       unsigned int JK,
                       double *K);

  MttkrpLut read_lookup_table(vector<dim_t> const &modes, int threads, bool gpu, bool suppress_warning)
  {
    auto lut = MttkrpLut();
    lut.lut_v = LUT_v(modes.size());

    std::string dir{};
    if (gpu)
    {
      dir = string(SOURCE_DIR) +
            "/data/GPU" +
            "/lookup_tables/" +
            cals::utils::mode_string(modes) + "/";
    } else
    {
      std::string blas_name = CALS_BACKEND;
      dir = string(SOURCE_DIR) +
            "/data/" +
            blas_name +
            "/lookup_tables/" +
            cals::utils::mode_string(modes) +
            "/" + std::to_string(threads) + "/";
    }

    for (auto m = 0lu; m < modes.size(); m++)
    {
      auto file = std::ifstream(dir + std::to_string(m), std::ios::in);
      if (file.is_open())
      {
        int rank = 0, mttkrp_method = 0;

        while (file >> rank >> mttkrp_method)
          lut.lut_v[m].insert(std::pair<int, int>(rank, mttkrp_method));
      } else
      {
        if (!suppress_warning)
          std::cout << "Lookup table missing for tensor: "
                    << utils::mode_string(modes)
                    << ", Threads: " << threads
                    << std::endl;
        return MttkrpLut();
      }
    }
    if (!lut.lut_v.empty())
    {
      lut.keys_v.reserve(lut.lut_v[0].size());
      for (auto &[key, val] : lut.lut_v[0])
        lut.keys_v.emplace_back(key);
    }
    return lut;
  }

  struct Mttkrp2StepGemmParams
  {
    // Intermediate Matrix
    dim_t inter_rows{};

    // DGEMM parameters
    CBLAS_TRANSPOSE trans_A{};
    dim_t n_blocks{};
    dim_t block_offset{};
    dim_t block_rows{};
    dim_t next_block{};
    dim_t stride{};
    dim_t B{};
  };

  struct Mttkrp2StepGemvParams
  {
    // DGEMV parameters
    CBLAS_TRANSPOSE trans_A{};
    dim_t A_rows{};
    dim_t A_cols{};
    dim_t stride{};
    dim_t x{};
    dim_t y{};
  };

  Matrix &khatri_rao(const Matrix &A, const Matrix &B, Matrix &K, KrpParams &params)
  {

#if true
    {
      const dim_t cols = K.get_cols();
      const dim_t rowsA = A.get_rows();
      const dim_t rowsB = B.get_rows();

#pragma omp parallel for // NOLINT(openmp-use-default-none)
      for (dim_t col = 0; col < cols; col++)
      {
#pragma omp parallel for // NOLINT(openmp-use-default-none)
        for (dim_t rowa = 0; rowa < rowsA; rowa++)
        {
          auto *const K_pointer = K.get_data() + col * K.get_col_stride() + rowa * rowsB;
          auto const factor = A(rowa, col);

// #pragma clang loop interleave_count(8)  // This needs testing.
          for (dim_t rowb = 0; rowb < rowsB; rowb++)
            K_pointer[rowb] = factor * B[rowb + rowsB * col];
        }
      }

      params.flops += 1llu * K.get_cols() * K.get_rows();
      params.memops += 1llu * K.get_cols() * K.get_rows() +
                       1llu * A.get_rows() * A.get_cols() +
                       1llu * B.get_rows() * B.get_cols();
      return K;
    }
#else
    //#if CUDA_ENABLED
    //    {
    //      const int cols = K.get_cols();
    //      const int rowsA = A.get_rows();
    //      const int rowsB = B.get_rows();
    //
    //      auto handle = cuda::cublas::create_handle();
    //      for (int col = 0; col < cols; col++)
    //        for (int row = 0; row < rowsA; row++)
    //        {
    //          auto *K_pointer = K.get_cudata() + col * K.get_col_stride() + row * rowsB;
    //          cublasDcopy(handle, rowsB, B.get_cudata() + col * B.get_col_stride(), 1, K_pointer, 1);
    //          const double val = A.at(row,col);
    //          cublasDscal(handle, rowsB, &val, K_pointer, 1);
    //        }
    //      cuda::cublas::destroy_handle(handle);
    //    }
    //    return K;
    //#else
    //    {
    //      const int cols = K.get_cols();
    //      const int rowsA = A.get_rows();
    //      const int rowsB = B.get_rows();
    //
    //#pragma omp parallel for // NOLINT(openmp-use-default-none)
    //      for (int col = 0; col < cols; col++)
    //        for (int row = 0; row < rowsA; row++)
    //        {
    //          auto *K_pointer = K.get_data() + col * K.get_col_stride() + row * rowsB;
    //          cblas_dcopy(rowsB, B.get_data() + col * B.get_col_stride(), 1, K_pointer, 1);
    //          cblas_dscal(rowsB, A.at(row, col), K_pointer, 1);
    //        }
    //    }
    //#if CUDA_ENABLED
    //    K.send_to_device();
    //#endif
    //    return K;
    //#endif
#endif
  };

  // Internal recursive implementation of the Khatri-Rao product.
  Matrix &
  khatri_rao_rec(std::stack<Matrix *> &remaining_targets, vector<Matrix> &workspace, int w_index, KrpParams &params)
  {
    if (!remaining_targets.empty())
    {
      auto *B = remaining_targets.top();
      remaining_targets.pop();
      auto &A = workspace[w_index - 1];
      auto &K = workspace[w_index];

      K.resize(A.get_rows() * B->get_rows(), A.get_cols());
      if (params.cuda)
      {
#if CUDA_ENABLED
        auto &handle = params.cuda_params.handle;
        auto &stream = params.cuda_params.streams[0];
        cublasSetStream(handle, stream);
        khatri_rao_cuda(A.get_cudata(), B->get_cudata(),
                        A.get_rows(), B->get_rows(), K.get_cols(),
                        K.get_cudata());
        cudaStreamSynchronize(stream);  // Make sure KRP is done
#else
        std::cerr << "Not compiled with CUDA support" << std::endl;
        exit(EXIT_FAILURE);
#endif
      } else
        khatri_rao(A, *B, K, params);

      params.flops += K.get_n_elements();
      params.memops += 1llu * A.get_n_elements() + B->get_n_elements() + K.get_n_elements();

      return khatri_rao_rec(remaining_targets, workspace, ++w_index, params);
    } else return workspace[w_index - 1];
  }

  // Compute ⨀ (n != mode) u.factor[n]
  Matrix &khatri_rao(Ktensor &u, vector<Matrix> &workspace, dim_t mode, KrpParams &params)
  {
    params.flops = 0;
    params.memops = 0;

    // TODO see if you can preallocate size of stack
    std::stack<Matrix *> remaining_targets;
    for (dim_t i = 0; i < u.get_n_modes(); i++)
      if (i != mode)
        remaining_targets.push(&(u.get_factor(i)));

    int w_index = 0;
    auto *u0 = remaining_targets.top();
    remaining_targets.pop();
    auto *u1 = remaining_targets.top();
    remaining_targets.pop();

    workspace[w_index].resize(u0->get_rows() * u1->get_rows(), u0->get_cols());

    if (params.cuda)
    {
#if CUDA_ENABLED
      auto &handle = params.cuda_params.handle;
      auto &stream = params.cuda_params.streams[0];
      cublasSetStream(handle, stream);
      khatri_rao_cuda((*u0).get_cudata(), (*u1).get_cudata(),
                      (*u0).get_rows(), (*u1).get_rows(), workspace[w_index].get_cols(),
                      workspace[w_index].get_cudata());
      cudaStreamSynchronize(stream);  // Make sure KRP is done
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    } else
      khatri_rao(*u0, *u1, workspace[w_index], params);

    params.flops += workspace[w_index].get_n_elements();
    params.memops += 1llu * (*u0).get_n_elements() + (*u1).get_n_elements() + workspace[w_index].get_n_elements();

    return khatri_rao_rec(remaining_targets, workspace, ++w_index, params);
  }

  Matrix &mttkrp_impl(const Tensor &X,
                      Ktensor &u,
                      vector<Matrix> &workspace,
                      int mode,
                      MttkrpParams &params)
  {
    // Explicitly generate the Khatri-Rao product
    params.mttkrp_timers[MTTKRP_TIMERS::MT_KRP].start();

    params.krp_params.cuda = params.cuda;
    auto &krp = khatri_rao(u, workspace, mode, params.krp_params);

    params.mttkrp_timers[MTTKRP_TIMERS::MT_KRP].stop();

    // Compute the matrix product
    params.mttkrp_timers[MTTKRP_TIMERS::MT_GEMM].start();

    // Figure out how to implicitly unfold the tensor as one or more matrix blocks.
    Unfolding unfolding = X.implicit_unfold(mode);

    if (params.cuda)
    {
#if CUDA_ENABLED
      auto &handle = params.cuda_params.handle;
      auto const &one = params.cuda_params.one;
      auto const &zero = params.cuda_params.zero;
      auto &stream = params.cuda_params.streams[0];

      cublasSetStream(handle, stream);
      if (mode != 0)
      {
        auto &factor = u.get_factor(mode);
        cudaMemsetAsync(factor.get_cudata(), 0, factor.get_n_elements() * sizeof(double), stream);
      }
      for (auto block_idx = 0; block_idx < unfolding.n_blocks; block_idx++)
      {
        // Locate block of the mode-th unfolding of X.
        double const *X_mode_blk_ptr = X.get_cudata() + block_idx * unfolding.block_offset;
        int ldX_mode = unfolding.stride;

        // Locate block of krp.
        double const *krp_blk_ptr = krp.get_cudata() + block_idx * unfolding.cols;

        auto &G = u.get_factor(mode);

        if (mode == 0)
          cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, G.get_rows(), G.get_cols(), unfolding.cols,
                      &one, X_mode_blk_ptr, ldX_mode, krp_blk_ptr, krp.get_col_stride(),
                      &zero, G.get_cudata(), G.get_col_stride());
        else
          cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, G.get_rows(), G.get_cols(), unfolding.cols,
                      &one, X_mode_blk_ptr, ldX_mode, krp_blk_ptr, krp.get_col_stride(),
                      &one, G.get_cudata(), G.get_col_stride());
      }
      u.get_factor(mode).receive_from_device_async(stream);
      cudaStreamSynchronize(stream);  // Make sure KRP is done
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    } else
    {
      if (mode != 0)
        u.get_factor(mode).zero();
      for (dim_t block_idx = 0; block_idx < unfolding.n_blocks; block_idx++)
      {
        /**
         MODE 0
         ====================
         Compute
           G := G + X * K

         where
           - G and K are column-major
           - X is column-major

         No transpositions needed.

         MODES 1, 2, ...
         ====================
         Compute
           G := G + X * K  or rather  G := G + (X')' * K

         where
           - X is row-major
           - G, K are col-major

         We transpose_copy X, since we are using col-major BLAS
         **/

        // Locate block of the mode-th unfolding of X.
        double const *X_mode_blk_ptr = X.get_data() + block_idx * unfolding.block_offset;
        dim_t ldX_mode = unfolding.stride;

        // Locate block of krp.
        double const *krp_blk_ptr = krp.get_data() + block_idx * unfolding.cols;

        auto &G = u.get_factor(mode);

        if (mode == 0)
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      G.get_rows(), G.get_cols(), unfolding.cols,
                      1.0, X_mode_blk_ptr, ldX_mode, krp_blk_ptr, krp.get_col_stride(),
                      0.0, G.get_data(), G.get_col_stride());
        else
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                      G.get_rows(), G.get_cols(), unfolding.cols,
                      1.0, X_mode_blk_ptr, ldX_mode, krp_blk_ptr, krp.get_col_stride(),
                      1.0, G.get_data(), G.get_col_stride());
      }
    }

    params.mttkrp_timers[MTTKRP_TIMERS::MT_GEMM].stop();

    auto &G = u.get_factor(mode);
    // GEMM  --  blocks * (2 * m * n * k + 2 * m * n)
    params.flops = 1llu * unfolding.n_blocks *
                   (2llu * G.get_rows() * G.get_cols() * unfolding.cols +
                    2llu * G.get_rows() * G.get_cols());
    params.memops = 1llu * unfolding.n_blocks *
                    (1llu * G.get_rows() * unfolding.cols +
                     1llu * G.get_cols() * unfolding.cols +
                     1llu * G.get_rows() * G.get_cols());

    params.flops += params.krp_params.flops;
    params.memops += params.krp_params.memops;

    return u.get_factor(mode);
  }

  Matrix &mttkrp_twostep_impl(const Tensor &X,
                              Ktensor &ktensor,
                              Matrix &workspace,
                              Mttkrp2StepGemmParams &gemm_params,
                              Mttkrp2StepGemvParams &gemv_params,
                              MttkrpParams &params)
  {
    params.mttkrp_timers[MTTKRP_TIMERS::TS_GEMM].start();
    auto &B = ktensor.get_factor(gemm_params.B);
    // Perform TTM
    workspace.resize(gemm_params.inter_rows, ktensor.get_components());

    if (params.cuda)
    {
#if CUDA_ENABLED
      auto const &handle = params.cuda_params.handle;
      auto const &one = params.cuda_params.one;
      auto const &zero = params.cuda_params.zero;

      cublasOperation_t opA = (gemm_params.trans_A == CblasNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
//      gemm_params.trans_A == CblasNoTrans ? opA = CUBLAS_OP_N : opA = CUBLAS_OP_T;

      for (auto blk = 0; blk < gemm_params.n_blocks; blk++)
      {
        cublasSetStream(handle, params.cuda_params.streams[blk % cuda_n_streams]);
        cublasDgemm(handle, opA, CUBLAS_OP_N,
                    gemm_params.block_rows, workspace.get_cols(), B.get_rows(), &one,
                    X.get_cudata() + blk * gemm_params.block_offset,
                    gemm_params.stride, // careful with the stride. Is it always this?
                    B.get_cudata(), B.get_col_stride(),
                    &zero, workspace.get_cudata() + blk * gemm_params.block_rows, workspace.get_col_stride());
      }
      for (auto &stream : params.cuda_params.streams)
        cudaStreamSynchronize(stream);  // Make sure all DGEMMS are done
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    } else
    {
#pragma omp parallel for // NOLINT(openmp-use-default-none)
      for (dim_t blk = 0; blk < gemm_params.n_blocks; blk++)
        cblas_dgemm(CblasColMajor, gemm_params.trans_A, CblasNoTrans,
                    gemm_params.block_rows, workspace.get_cols(), B.get_rows(), 1.0,
                    X.get_data() + blk * gemm_params.block_offset,
                    gemm_params.stride, // careful with the stride. Is it always this?
                    B.get_data(), B.get_col_stride(),
                    0.0, workspace.get_data() + blk * gemm_params.block_rows, workspace.get_col_stride());
    }

    params.mttkrp_timers[MTTKRP_TIMERS::TS_GEMM].stop();

    // Perform series of TTVs
    params.mttkrp_timers[MTTKRP_TIMERS::TS_GEMV].start();

    auto &x = ktensor.get_factor(gemv_params.x);
    auto &y = ktensor.get_factor(gemv_params.y);

    if (params.cuda)
    {
#if CUDA_ENABLED
      auto const &handle = params.cuda_params.handle;
      auto const &one = params.cuda_params.one;
      auto const &zero = params.cuda_params.zero;

      cublasOperation_t opA = (gemv_params.trans_A == CblasNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
      //gemv_params.trans_A == CblasNoTrans ? opA = CUBLAS_OP_N : opA = CUBLAS_OP_T;

      bool flag = false;
      double *temp;
      if (opA == CUBLAS_OP_N)
      {
        flag = true;
        opA = CUBLAS_OP_T;
        std::swap(gemv_params.A_rows, gemv_params.A_cols);
        cuda::allocate(temp, workspace.get_n_elements());
      }

      for (auto col = 0; col < workspace.get_cols(); col++)
      {
        cublasSetStream(handle, params.cuda_params.streams[col % cuda_n_streams]);
        if (flag)
        {
          cublasDgeam(handle, CUBLAS_OP_T, CUBLAS_OP_N,
                      gemv_params.A_rows, gemv_params.A_cols, &one,
                      workspace.get_cudata()+ col * workspace.get_col_stride(), gemv_params.A_cols,
                      &zero, temp + col * workspace.get_col_stride(), gemv_params.A_rows,
                      temp + col * workspace.get_col_stride(), gemv_params.A_rows);
          cublasDgemv(handle, opA, gemv_params.A_rows, gemv_params.A_cols, &one,
                      temp + col * workspace.get_col_stride(), gemv_params.A_rows,
                      x.get_cudata() + col * x.get_col_stride(), 1, &zero,
                      y.get_cudata() + col * y.get_col_stride(), 1);
        }
        else
          cublasDgemv(handle, opA, gemv_params.A_rows, gemv_params.A_cols, &one,
                      workspace.get_cudata() + col * workspace.get_col_stride(), gemv_params.stride,
                      x.get_cudata() + col * x.get_col_stride(), 1, &zero,
                      y.get_cudata() + col * y.get_col_stride(), 1);
      }

      for (auto &stream : params.cuda_params.streams)
        cudaStreamSynchronize(stream);  // Make sure all DGEMVs are done
      if (flag)
        cuda::deallocate(temp);
      y.receive_from_device();
#else
      std::cerr << "Not compiled with CUDA support" << std::endl;
      exit(EXIT_FAILURE);
#endif
    } else
    {
#pragma omp parallel for // NOLINT(openmp-use-default-none)
      for (dim_t col = 0; col < workspace.get_cols(); col++)
        cblas_dgemv(CblasColMajor, gemv_params.trans_A, gemv_params.A_rows, gemv_params.A_cols, 1.0,
                    workspace.get_data() + col * workspace.get_col_stride(), gemv_params.stride,
                    x.get_data() + col * x.get_col_stride(), 1, 0.0,
                    y.get_data() + col * y.get_col_stride(), 1);
    }

    params.mttkrp_timers[MTTKRP_TIMERS::TS_GEMV].stop();

    // GEMM  --  blocks * (2 * m * n * k + 2 * m * n)
    params.flops = 1llu * gemm_params.n_blocks *
                   (2llu * gemm_params.block_rows * workspace.get_cols() * B.get_rows() +
                    2llu * gemm_params.block_rows * workspace.get_cols());
    params.memops = 1llu * gemm_params.n_blocks *
                    (1llu * gemm_params.block_rows * B.get_rows() +
                     1llu * B.get_rows() * workspace.get_cols() +
                     1llu * gemm_params.block_rows * workspace.get_cols());

    // GEMVs --  cols * (2 * m * n)
    params.flops += workspace.get_cols() * 2llu * gemv_params.A_rows * gemv_params.A_cols;
    params.memops +=
        workspace.get_cols() * (1llu * gemv_params.A_rows * gemv_params.A_cols + 2llu * gemv_params.A_cols);

    return y;
  }

  Matrix &mttkrp_twostep(const Tensor &X,
                         Ktensor &u,
                         Matrix &workspace,
                         int mode,
                         MttkrpParams &params)
  {
    assert(u.get_n_modes() == 3);  // Twostep method only implemented for 3D tensors

    auto &U_0 = u.get_factor(0);
    auto &U_1 = u.get_factor(1);
    auto &U_2 = u.get_factor(2);

    Mttkrp2StepGemmParams gemm_params;
    Mttkrp2StepGemvParams gemv_params;

    if (mode == 2)
    {
      if (params.method == MTTKRP_METHOD::TWOSTEP0) // Option 0 (JK x I) * (I * R)
      {
        gemm_params.inter_rows = U_1.get_rows() * U_2.get_rows();
        gemm_params.trans_A = CblasTrans;
        gemm_params.stride = U_0.get_rows();
        gemm_params.block_rows = gemm_params.inter_rows;
        gemm_params.n_blocks = 1;
        gemm_params.B = 0;
        gemv_params.trans_A = CblasTrans;
        gemv_params.A_rows = U_1.get_rows();
        gemv_params.A_cols = U_2.get_rows();
        gemv_params.stride = U_1.get_rows();
        gemv_params.x = 1;
        gemv_params.y = 2;
      } else  // Option 1 (IK x J) * (J * R)
      {
        gemm_params.inter_rows = U_0.get_rows() * U_2.get_rows();
        gemm_params.trans_A = CblasNoTrans;
        gemm_params.block_offset = U_0.get_rows() * U_1.get_rows();
        gemm_params.stride = U_0.get_rows();
        gemm_params.block_rows = U_0.get_rows();
        gemm_params.n_blocks = U_2.get_rows();
        gemm_params.B = 1;

        gemv_params.trans_A = CblasTrans;
        gemv_params.A_rows = U_0.get_rows();
        gemv_params.A_cols = U_2.get_rows();
        gemv_params.stride = U_0.get_rows();
        gemv_params.x = 0;
        gemv_params.y = 2;
      }
    } else if (mode == 1)
    {
      if (params.method == MTTKRP_METHOD::TWOSTEP0)  // Option 0 (IJ x K) * (K * R)
      {
        gemm_params.inter_rows = U_0.get_rows() * U_1.get_rows();
        gemm_params.trans_A = CblasNoTrans;
        gemm_params.stride = U_0.get_rows() * U_1.get_rows();
        gemm_params.block_rows = gemm_params.inter_rows;
        gemm_params.n_blocks = 1;
        gemm_params.B = 2;
        gemv_params.trans_A = CblasTrans;
        gemv_params.A_rows = U_0.get_rows();
        gemv_params.A_cols = U_1.get_rows();
        gemv_params.stride = U_0.get_rows();
        gemv_params.x = 0;
        gemv_params.y = 1;
      } else  // Option 1 (JK x I) * (I * R)
      {
        gemm_params.inter_rows = U_1.get_rows() * U_2.get_rows();
        gemm_params.trans_A = CblasTrans;
        gemm_params.stride = U_0.get_rows();
        gemm_params.block_rows = gemm_params.inter_rows;
        gemm_params.n_blocks = 1;
        gemm_params.B = 0;
        gemv_params.trans_A = CblasNoTrans;
        gemv_params.A_rows = U_1.get_rows();
        gemv_params.A_cols = U_2.get_rows();
        gemv_params.stride = U_1.get_rows();
        gemv_params.x = 2;
        gemv_params.y = 1;
      }
    } else if (mode == 0)
    {
      if (params.method == MTTKRP_METHOD::TWOSTEP0) // Option 0 (IJ x K) * (K * R)
      {
        gemm_params.inter_rows = U_0.get_rows() * U_1.get_rows();
        gemm_params.trans_A = CblasNoTrans;
        gemm_params.stride = U_0.get_rows() * U_1.get_rows();
        gemm_params.block_rows = gemm_params.inter_rows;
        gemm_params.n_blocks = 1;
        gemm_params.B = 2;
        gemv_params.trans_A = CblasNoTrans;
        gemv_params.A_rows = U_0.get_rows();
        gemv_params.A_cols = U_1.get_rows();
        gemv_params.stride = U_0.get_rows();
        gemv_params.x = 1;
        gemv_params.y = 0;
      } else  // Option 1 (IK x J) * (J * R)
      {
        gemm_params.inter_rows = U_0.get_rows() * U_2.get_rows();
        gemm_params.trans_A = CblasNoTrans;
        gemm_params.block_offset = U_0.get_rows() * U_1.get_rows();
        gemm_params.stride = U_0.get_rows();
        gemm_params.block_rows = U_0.get_rows();
        gemm_params.n_blocks = U_2.get_rows();
        gemm_params.B = 1;

        gemv_params.trans_A = CblasNoTrans;
        gemv_params.A_rows = U_0.get_rows();
        gemv_params.A_cols = U_2.get_rows();
        gemv_params.stride = U_0.get_rows();
        gemv_params.x = 2;
        gemv_params.y = 0;
      }
    } else
    {
      std::cerr << "Illegal mode given (" << mode << ") for a " << X.get_n_modes() << "D tensor." << std::endl;
      abort();
    }

    return mttkrp_twostep_impl(X, u, workspace, gemm_params, gemv_params, params);
  }

  Matrix &mttkrp(const Tensor &X,
                 Ktensor &u,
                 vector<Matrix> &workspace,
                 int mode,
                 MttkrpParams &params)
  {
    // Make sure all timers are zero, bc some of the timers are not going to be set later, depending on chosen branch.
    for (auto i = 0; i < params.mttkrp_timers.LENGTH; i++)
      params.mttkrp_timers.timers[i].reset();

    if (X.get_n_modes() != 3)
      mttkrp_impl(X, u, workspace, mode, params);
    else
    {
      if (params.method == MTTKRP_METHOD::MTTKRP)
        mttkrp_impl(X, u, workspace, mode, params);
      else if ((params.method == MTTKRP_METHOD::TWOSTEP0)
               || (params.method == MTTKRP_METHOD::TWOSTEP1))
        mttkrp_twostep(X, u, workspace[0], mode, params);  // workspace[0] is size maxmode1*maxmode2 x rank
      else if (params.method == MTTKRP_METHOD::AUTO)
      {
          if (!params.lut.lut_v.empty() && !params.lut.keys_v.empty())  // If LUT exists, read it
          {
            auto key = std::lower_bound(params.lut.keys_v.cbegin(), params.lut.keys_v.cend(), u.get_components());
            if (key == params.lut.keys_v.end())
              key = params.lut.keys_v.end() - 1;
            auto val = params.lut.lut_v[mode].at(*key);
            if (val == MTTKRP_METHOD::MTTKRP)
              mttkrp_impl(X, u, workspace, mode, params);
            else if (val == MTTKRP_METHOD::TWOSTEP0 || val == MTTKRP_METHOD::TWOSTEP1)
            {
              params.method = static_cast<MTTKRP_METHOD>(val);
              mttkrp_twostep(X, u, workspace[0], mode, params);
              params.method = MTTKRP_METHOD::AUTO;
            }
          } else  // If LUT does not exist, revert to common sense defaults
          {
            if (params.cuda)  // If CUDA, MTTKRP is the best alternative
              mttkrp_impl(X, u, workspace, mode, params);
            else
              if (get_threads() != 1) // If multi-threaded, never Twostep0
              {
                if (mode == 1)
                  mttkrp_impl(X, u, workspace, mode, params);
                else
                {
                  params.method = MTTKRP_METHOD::TWOSTEP1;
                  mttkrp_twostep(X, u, workspace[0], mode, params);  // workspace[0] is (maxmode1*maxmode2 x rank)
                  params.method = MTTKRP_METHOD::AUTO;
                }
              } else  // If single threaded, always Twostep0
              {
                params.method = MTTKRP_METHOD::TWOSTEP0;
                mttkrp_twostep(X, u, workspace[0], mode, params);  // workspace[0] is (maxmode1*maxmode2 x rank)
                params.method = MTTKRP_METHOD::AUTO;
              }
          }
      } else
      {
        std::cerr << "Invalid MTTKRP method selected." << std::endl;
        abort();
      }
    }
    return u.get_factor(mode);
  }

}
