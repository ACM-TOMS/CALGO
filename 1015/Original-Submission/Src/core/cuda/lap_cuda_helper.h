#pragma once

#include <iostream>

#ifndef checkCudaErrors
#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

// Error Code string definitions here
typedef struct
{
	char const *error_string;
	int  error_id;
} s_CudaErrorStr;

inline void __checkCudaErrors(cudaError_t err, const char *file, const int line)
{
	if ((cudaError_t)0 != err)
	{
		std::cerr << "checkCudaErrors() error = " << err << " \"" << cudaGetErrorString(err) << "\" from file <" << file << ">, line " << line << "." << std::endl;
		exit(EXIT_FAILURE);
	}
}
#endif
