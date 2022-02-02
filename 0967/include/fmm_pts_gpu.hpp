#ifndef _CUDA_FUNC_HPP_
#define _CUDA_FUNC_HPP_

#ifdef __cplusplus
extern "C" {
#endif
  void  in_perm_gpu_f(char* precomp_data, float *  input_data, char* buff_in , size_t*  input_perm, size_t vec_cnt, size_t M_dim0, cudaStream_t* stream);
  void  in_perm_gpu_d(char* precomp_data, double*  input_data, char* buff_in , size_t*  input_perm, size_t vec_cnt, size_t M_dim0, cudaStream_t* stream);

  void out_perm_gpu_f(char* precomp_data, float * output_data, char* buff_out, size_t* output_perm, size_t vec_cnt, size_t M_dim1, cudaStream_t* stream);
  void out_perm_gpu_d(char* precomp_data, double* output_data, char* buff_out, size_t* output_perm, size_t vec_cnt, size_t M_dim1, cudaStream_t* stream);
#ifdef __cplusplus
}
#endif

template <class Real_t>
void  in_perm_gpu(char* precomp_data, Real_t*  input_data, char* buff_in , size_t*  input_perm, size_t vec_cnt, size_t M_dim0, cudaStream_t* stream);

template <class Real_t>
void out_perm_gpu(char* precomp_data, Real_t* output_data, char* buff_out, size_t* output_perm, size_t vec_cnt, size_t M_dim1, cudaStream_t* stream);

template<> inline void  in_perm_gpu<float >(char* precomp_data, float *  input_data, char* buff_in , size_t*  input_perm, size_t vec_cnt, size_t M_dim0, cudaStream_t* stream){
  in_perm_gpu_f (precomp_data,  input_data, buff_in ,  input_perm, vec_cnt, M_dim0, stream);
}
template<> inline void  in_perm_gpu<double>(char* precomp_data, double*  input_data, char* buff_in , size_t*  input_perm, size_t vec_cnt, size_t M_dim0, cudaStream_t* stream){
  in_perm_gpu_d (precomp_data,  input_data, buff_in ,  input_perm, vec_cnt, M_dim0, stream);
}

template<> inline void out_perm_gpu<float >(char* precomp_data, float * output_data, char* buff_out, size_t* output_perm, size_t vec_cnt, size_t M_dim1, cudaStream_t* stream){
  out_perm_gpu_f(precomp_data, output_data, buff_out, output_perm, vec_cnt, M_dim1, stream);
}
template<> inline void out_perm_gpu<double>(char* precomp_data, double* output_data, char* buff_out, size_t* output_perm, size_t vec_cnt, size_t M_dim1, cudaStream_t* stream){
  out_perm_gpu_d(precomp_data, output_data, buff_out, output_perm, vec_cnt, M_dim1, stream);
}

#endif //_CUDA_FUNC_HPP_
