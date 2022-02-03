#ifndef _COMMON_H
#define _COMMON_H

#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))

#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)

/* Reduce the number of registers in use to avoid register spilling
 * and extra load/store when hardware counters are being measured */
#ifdef HWC

# define c00  0.5

# define c14z 0.5
# define c13z 0.5
# define c12z 0.5
# define c11z 0.5
# define c10z 0.5
# define c9z  0.5
# define c8z  0.5
# define c7z  0.5
# define c6z  0.5
# define c5z  0.5
# define c4z  0.5
# define c3z  0.5
# define c2z  0.5
# define c1z  0.5

# define c14x 0.5
# define c13x 0.5
# define c12x 0.5
# define c11x 0.5
# define c10x 0.5
# define c9x  0.5
# define c8x  0.5
# define c7x  0.5
# define c6x  0.5
# define c5x  0.5
# define c4x  0.5
# define c3x  0.5
# define c2x  0.5
# define c1x  0.5

# define c14y 0.5
# define c13y 0.5
# define c12y 0.5
# define c11y 0.5
# define c10y 0.5
# define c9y  0.5
# define c8y  0.5
# define c7y  0.5
# define c6y  0.5
# define c5y  0.5
# define c4y  0.5
# define c3y  0.5
# define c2y  0.5
# define c1y  0.5

#else

#if defined(SSE) // TODO

# define VECSIZE 2

#elif defined(AVX)

# define VECSIZE 4

# define c00v _mm256_set1_pd(0.5)

# define c4zv _mm256_set1_pd(1.4)
# define c3zv _mm256_set1_pd(1.3)
# define c2zv _mm256_set1_pd(1.2)
# define c1zv _mm256_set1_pd(1.1)

# define c4xv _mm256_set1_pd(2.4)
# define c3xv _mm256_set1_pd(2.3)
# define c2xv _mm256_set1_pd(2.2)
# define c1xv _mm256_set1_pd(2.1)

# define c4yv _mm256_set1_pd(3.4)
# define c3yv _mm256_set1_pd(3.3)
# define c2yv _mm256_set1_pd(3.2)
# define c1yv _mm256_set1_pd(3.1)

#elif defined(MIC)

# define VECSIZE 8

# define c00v _mm512_set1_pd(0.5)

# define c4zv _mm512_set1_pd(1.4)
# define c3zv _mm512_set1_pd(1.3)
# define c2zv _mm512_set1_pd(1.2)
# define c1zv _mm512_set1_pd(1.1)

# define c4xv _mm512_set1_pd(2.4)
# define c3xv _mm512_set1_pd(2.3)
# define c2xv _mm512_set1_pd(2.2)
# define c1xv _mm512_set1_pd(2.1)

# define c4yv _mm512_set1_pd(3.4)
# define c3yv _mm512_set1_pd(3.3)
# define c2yv _mm512_set1_pd(3.2)
# define c1yv _mm512_set1_pd(3.1)

#else

# define c00 0.5 //(0.5/(fac*fac))

# define c14z 1.14
# define c13z 1.13
# define c12z 1.12
# define c11z 1.11
# define c10z 1.10
# define c9z  1.9
# define c8z  1.8
# define c7z  1.7
# define c6z  1.6
# define c5z  1.5
# define c4z  1.4
# define c3z  1.3
# define c2z  1.2
# define c1z  1.1

# define c14x 2.14
# define c13x 2.13
# define c12x 2.12
# define c11x 2.11
# define c10x 2.10
# define c9x  2.9
# define c8x  2.8
# define c7x  2.7
# define c6x  2.6
# define c5x  2.5
# define c4x  2.4
# define c3x  2.3
# define c2x  2.2
# define c1x  2.1

# define c14y 3.14
# define c13y 3.13
# define c12y 3.12
# define c11y 3.11
# define c10y 3.10
# define c9y  3.9
# define c8y  3.8
# define c7y  3.7
# define c6y  3.6
# define c5y  3.5
# define c4y  3.4
# define c3y  3.3
# define c2y  3.2
# define c1y  3.1

#endif // VEC

#endif // HWC

#ifdef VECSIZE
#define PAD(_LENGTH_) (_LENGTH_ % VECSIZE)
#else
#define PAD(_LENGTH_) 0
#endif

#define VERSION "0.4.9-" __DATE__

#endif // _COMMON_H

