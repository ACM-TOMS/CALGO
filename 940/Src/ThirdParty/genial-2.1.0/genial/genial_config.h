//GENIAL - GENeric Image & Array Library
//Copyright (C) 2005  IENT - RWTH Aachen
//
//This program is free software; you can redistribute it and/or
//modify it under the terms of the GNU General Public License
//as published by the Free Software Foundation; either version 2
//of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program; if not, write to the Free Software
//Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#ifndef GENIAL_CONFIG_H
#define GENIAL_CONFIG_H


#if (defined(_MSC_VER) && defined(__ICL))
#pragma warning(disable : 1572)
#ifndef NDEBUG
#pragma warning(disable : 963 964)
#endif
#endif

#if (defined(_MSC_VER) && !defined(__ICL))
#pragma warning(disable : 4200 4244 4521 4793 4996 )
#ifndef NDEBUG
#pragma warning(disable : 4799 )
#endif
#endif


// Uncomment to simulate the library optimizations, even in debug mode.
// But assertions won't work anymore.
#ifndef NDEBUG
//#define NDEBUG
#endif


// Uncomment to try to force inlining
#if defined(__ICL) || defined(_MSC_VER)
  #define inline __forceinline
#elif defined(__ICC)
  //#define inline __forceinline
  //#define inline __inline
  //#define inline __attribute__((__always_inline__)) // Have obviously no effect on ICC, bad...
#elif defined(__GNUC__)
  //#define inline __inline
  //#define inline __attribute__((__always_inline__)) // GCC4.2 seems to respect this attribute, unlike former version...
#endif


// Static/Shared library setting
#if defined(__ICL) || defined(_MSC_VER)
  #ifdef GENIAL_SHARED // should only be eventually defined in the library project/makefile
    #ifdef GENIAL_EXPORT // should only be eventually defined in the library project/makefile
    #define GENIAL_API __declspec(dllexport)
    #else
    #define GENIAL_API __declspec(dllimport)
    #endif
  #else
    #define GENIAL_API
  #endif
#elif defined(__GNUC__)
  #define GENIAL_API
#else
  #define GENIAL_API
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Uncomment for SIMD calculation
#define HAS_MMX
#define HAS_SEE
#define HAS_SSE2
//#define HAS_SSE3


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the maximum number of threads that can be simultaneously executed by built-in multi-threading implementation.
// Multi-threading causes a small overhead involved in synchhronizing threads.
//   0 : Automatically retreaves the number of physical cores (only works with Intel & AMD processors)
//   1 : No muti-threading
//   2,3,4,...: Maximum number of simultaneous threads
// Default value if not defined: 0
#ifndef MAX_THREADS
#define MAX_THREADS 1
#endif

// Uncomment to use the pthread library instead of the native libraries of Windows and MacOS.
// For Windows, the pthreads-w32 library has to be installed.
#if !defined(PTHREADS)
  //#define PTHREADS
#endif

// Uncomment to use OpenMP instead of the built-in multi-threading implementation.
// Do not forget to set the corresponding option for your compiler
// (ICL=/Qopenmp, VC2005=/openmp, GCC4.2=-openmp,...)
#if !defined(OPENMP)
//#define OPENMP
#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define it to precompile the FFT algorithms
// and let it defined to link your programme with the compiled library
//#define FFT_PRECOMPILE

// Set the complexity level of the fft algorithm. Possible values: 2,3,4,5,8,16,32,64.
// Default value if not defined: 8.
#ifndef FFT_LEVEL
//#define FFT_LEVEL 64
#endif

// Define it to use multi-threaded algorithms for 2D-FFT
// Multi-threading causes an overhead involved in synchronization, but increases
// execution speed for big matrices on systems with several processors/cores.
#ifndef FFT_THREADING
//#define FFT_THREADING
#endif

// Uncomment if the FFT has to be used on signals which lengths are only powers of 2.
//  => quicker compilation and smaller executable.
#ifndef FFT_ONLY_POW2
//#define FFT_ONLY_POW2
#endif

// Set the number of twiddles tables that can be kept in memory.
// Increase this number to maximize execution speed, if you
// foresee to alternately calculate the FFT on different array sizes.
// Default value if not defined: 1
#ifndef FFT_TWIDDLES
#define FFT_TWIDDLES 5
#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define it to precompile the BLAS routines
// and let it defined to link your programme with the compiled library
//#define BLAS_PRECOMPILE

// Define it to use multi-threaded algorithms for matrix multiplications
// Multi-threading causes an overhead involved in synchronization, but increases
// execution speed for big matrices on systems with several processors/cores.
#ifndef BLAS_THREADING
//#define BLAS_THREADING
#endif

// BLAS_M, BLAS_N, BLAS_K and BLAS_K2 are dimensions of block matrices used to optimize the cache memory.
// Try various combinations to optimize speed on your system.
// - BLAS_M and BLAS_N must be multiple of 12.
// - BLAS_K and BLAS_K2 must be smaller or equal to 25.
// - BLAS_K2 can be small if you do not care about execution speed for small matrices.
#if !defined(BLAS_M) && !defined(BLAS_N) && !defined(BLAS_K) && !defined(BLAS_K2)
  #if defined(HAS_SSE) || defined(HAS_SSE2) || defined(HAS_SSE3) // seems best for Pentium4 and seems ok for Core2Duo
    #define BLAS_M   48
    #define BLAS_N   48
    #define BLAS_K   15
    #define BLAS_K2  15
  #else
    #define BLAS_M   48 // probably not the best, but you obviously do not care about speed if you do not use SIMD instructions.
    #define BLAS_N   48
    #define BLAS_K   15
    #define BLAS_K2  15
  #endif
#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define it to precompile the Motion Estimation routines during the installation
// and let it defined to link your programme with the compiled library
//#define MOTION_ESTIMATION_PRECOMPILE

// Set the maximum height of search areas. Possible values: 8,16,32,64.
// Default value if not defined: 32
// Set it to a low value if you won't use bigger search areas
// => quicker compilation and smaller executable.
#ifndef MOTION_ESTIMATION_LEVEL
//#define MOTION_ESTIMATION_LEVEL 64
#endif



#endif
