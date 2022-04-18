/* Stencil naive codes for stencilprobe test
 * Raúl de la Cruz (delacruz@bsc.es)
 * Mauricio Araya-Polo (mauricio.araya@bsc.es)
 * Barcelona Supercomputing Center (BSC)
 *
 * These stencils are based on Chombo's heattut code.
 * (Only the spatial differential operator is computed)
 * Five different lengths of stencils are implemented:
 *
 *   * 7pt stencil (1pt / Original Chombo's heattut)
 *   * 13pt stencil (2pt length per direction)
 *   * 25pt stencil (4pt length per direction)
 *   * 43pt stencil (7pt length per direction)
 *   * 85pt stencil (14pt length per direction)
 *
 * All the points in the stencil are weighted with
 * specific coefficients (cl[z|x|y] coeffs in common.h)
 * (l is the order of the stencil and z, x, y the axies).
 * The operation applied to the domain is the following:
 *
 *                        Y
 *             clz *        * cly
 *                 ·      ·
 *             c2z *     * c2y       Anext[i,j,k] = c00 *  A0[i  ,j  ,k  ] +
 *                 |   /                            c1z * (A0[i+1,j  ,k  ] + A0[i-1,j  ,k  ]) +
 *             c1z *  * c1y                         ···
 *             c00 |/  c1x c2x clx                  clz * (A0[i+l,j  ,k  ] + A0[i-l,j  ,k  ]) +
 *     *···*---*---O---*---*···*                    c1x * (A0[i  ,j+1,k  ] + A0[i  ,j-1,k  ]) +
 *   clx c2x c1x  /|            X                   ···
 *          c1y *  * c1z                            clx * (A0[i  ,j+l,k  ] + A0[i  ,j-l,k  ]) +
 *             /   |                                c1y * (A0[i  ,j  ,k+1] + A0[i  ,j  ,k-1]) +
 *       c2y *     * c2z                            ···
 *          ·      ·                                cly * (A0[i  ,j  ,k+l] + A0[i  ,j  ,k-l])
 *    cly *        * clz                             
 *               Z
 *
 * Notice that seismic nomenclature is used for these codes:
 * (Z unit-stride dimension, X and Y largest stride dimensions)
 *
 * Three versions for each stencil are implemented:
 *  * Fussion: 1 loop  (by default)
 *  * Fission: 2 loops (OPTS=-DFISSION_2LOOPS)
 *  * Fission: 3 loops (OPTS=-DFISSION_3LOOPS)
 *
 * Version date: November 2009
 *
 * Changelog
 * - September 2011: Added support for factored and expanded
 *   expressions of the stencil (under research)
 *     Use: -DFACTORED
 * - January 2013: Intel pragmas for autovectorization
 * - June 2013: Initial SSE/AVX support for Sandy Bridge (only STENCIL=4)
 *     Use: -DSSE and -DAVX respectively
 * - July 2013: Temporal/non-temporal stores for SSE/AVX
 *     Use: -DSTREAM to active non-temporal stores (cache-bypass)
 * - August 2013: Intel MIC implementation with mm512 intrinsics
 *     Use: -DMIC
 *
 */


#ifndef _STENCIL_H_
#define _STENCIL_H_


#ifdef SSE
# include <xmmintrin.h>
#endif

#if defined(AVX) || defined(MIC)
# include <immintrin.h>
#endif


#if !defined(FACTORED)
# if (!defined(SSE) && !defined(AVX) && !defined(MIC))

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_ ## length( Anext, A0 )

/* Stencil for fission in 2 loops */
#   define STENCIL1_2LOOPS( Anext, A0, length ) STENCIL1_2LOOPS_ ## length( Anext, A0 )
#   define STENCIL2_2LOOPS( Anext, A0, length ) STENCIL2_2LOOPS_ ## length( Anext, A0 )

/* Stencil for fission in 3 loops */
#   define STENCIL1_3LOOPS( Anext, A0, length ) STENCIL1_3LOOPS_ ## length( Anext, A0 )
#   define STENCIL2_3LOOPS( Anext, A0, length ) STENCIL2_3LOOPS_ ## length( Anext, A0 )
#   define STENCIL3_3LOOPS( Anext, A0, length ) STENCIL3_3LOOPS_ ## length( Anext, A0 )

#   define IINC 1

# elif defined(SSE)

#   error "SSE version not yet implemented!"

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_SSE_ ## length( Anext, A0 )

#   define IINC 2

# elif defined(AVX)

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_AVX_ ## length( Anext, A0 )

/* Stencil for fission in 2 loops */
#   define STENCIL1_2LOOPS( Anext, A0, length ) STENCIL1_2LOOPS_AVX_ ## length( Anext, A0 )
#   define STENCIL2_2LOOPS( Anext, A0, length ) STENCIL2_2LOOPS_AVX_ ## length( Anext, A0 )

/* Stencil for fission in 3 loops */
#   define STENCIL1_3LOOPS( Anext, A0, length ) STENCIL1_3LOOPS_AVX_ ## length( Anext, A0 )
#   define STENCIL2_3LOOPS( Anext, A0, length ) STENCIL2_3LOOPS_AVX_ ## length( Anext, A0 )
#   define STENCIL3_3LOOPS( Anext, A0, length ) STENCIL3_3LOOPS_AVX_ ## length( Anext, A0 )

#   define IINC 4

#   ifndef STREAM
#     define _mm256_stream_pd _mm256_store_pd
#   endif

# elif defined(MIC)

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_MIC_ ## length( Anext, A0 )

/* Stencil for fission in 2 loops */
#   define STENCIL1_2LOOPS( Anext, A0, length ) STENCIL1_2LOOPS_MIC_ ## length( Anext, A0 )
#   define STENCIL2_2LOOPS( Anext, A0, length ) STENCIL2_2LOOPS_MIC_ ## length( Anext, A0 )

/* Stencil for fission in 3 loops */
#   define STENCIL1_3LOOPS( Anext, A0, length ) STENCIL1_3LOOPS_MIC_ ## length( Anext, A0 )
#   define STENCIL2_3LOOPS( Anext, A0, length ) STENCIL2_3LOOPS_MIC_ ## length( Anext, A0 )
#   define STENCIL3_3LOOPS( Anext, A0, length ) STENCIL3_3LOOPS_MIC_ ## length( Anext, A0 )

#   define IINC 8

#   ifndef STREAM
#     define _mm512_storenr_pd _mm512_store_pd
#   endif

# endif
#else
# if (!defined(SSE) && !defined(AVX) && !defined(MIC))

/* Stencil factored case */
/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_FACTORED_ ## length( Anext, A0 )

/* Stencil for fission in 2 loops */
#   define STENCIL1_2LOOPS( Anext, A0, length ) STENCIL1_2LOOPS_FACTORED_ ## length( Anext, A0 )
#   define STENCIL2_2LOOPS( Anext, A0, length ) STENCIL2_2LOOPS_FACTORED_ ## length( Anext, A0 )

/* Stencil for fission in 3 loops */
#   define STENCIL1_3LOOPS( Anext, A0, length ) STENCIL1_3LOOPS_FACTORED_ ## length( Anext, A0 )
#   define STENCIL2_3LOOPS( Anext, A0, length ) STENCIL2_3LOOPS_FACTORED_ ## length( Anext, A0 )
#   define STENCIL3_3LOOPS( Anext, A0, length ) STENCIL3_3LOOPS_FACTORED_ ## length( Anext, A0 )

#   define IINC 1

# elif defined(SSE)

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_FACTORED_SSE_ ## length( Anext, A0 )

#   define IINC 2

# elif defined(AVX)

/* Stencil for 1 loop */
#   define STENCIL( Anext, A0, length )        STENCIL_FACTORED_AVX_ ## length( Anext, A0 )

#   define IINC 4

#   ifndef STREAM
#     define _mm256_stream_pd _mm256_store_pd
#   endif

# endif
#endif



/* 7pt stencil */
#define STENCIL_1( Anext, A0 )                   \
          Anext[Index3D (nx, ny, i, j, k)] =     \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_2LOOPS_1( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)];

#define STENCIL2_2LOOPS_1( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_3LOOPS_1( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)];

#define STENCIL2_3LOOPS_1( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)];

#define STENCIL3_3LOOPS_1( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

/* 7pt stencil (FACTORED) */
#define STENCIL_FACTORED_1( Anext, A0 )          \
          Anext[Index3D (nx, ny, i, j, k)] =     \
             c1y * (A0[Index3D (nx, ny, i, j, k + 1)] + A0[Index3D (nx, ny, i, j, k - 1)]) + \
             c1x * (A0[Index3D (nx, ny, i, j + 1, k)] + A0[Index3D (nx, ny, i, j - 1, k)]) + \
             c1z * (A0[Index3D (nx, ny, i + 1, j, k)] + A0[Index3D (nx, ny, i - 1, j, k)]) + \
             c00 * A0[Index3D (nx, ny, i, j, k)];


/* 13pt stencil */
#define STENCIL_2( Anext, A0 )                         \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_2LOOPS_2( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)];

#define STENCIL2_2LOOPS_2( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_3LOOPS_2( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)];

#define STENCIL2_3LOOPS_2( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)];

#define STENCIL3_3LOOPS_2( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

/* 13pt stencil (FACTORED) */
#define STENCIL_FACTORED_2( Anext, A0 )                         \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * (A0[Index3D (nx, ny, i, j, k + 1)] + A0[Index3D (nx, ny, i, j, k - 1)]) + \
             c2y * (A0[Index3D (nx, ny, i, j, k + 2)] + A0[Index3D (nx, ny, i, j, k - 2)]) + \
             c1x * (A0[Index3D (nx, ny, i, j + 1, k)] + A0[Index3D (nx, ny, i, j - 1, k)]) + \
             c2x * (A0[Index3D (nx, ny, i, j + 2, k)] + A0[Index3D (nx, ny, i, j - 2, k)]) + \
             c1z * (A0[Index3D (nx, ny, i + 1, j, k)] + A0[Index3D (nx, ny, i - 1, j, k)]) + \
             c2z * (A0[Index3D (nx, ny, i + 2, j, k)] + A0[Index3D (nx, ny, i - 2, j, k)]) + \
             c00 * A0[Index3D (nx, ny, i, j, k)];


/* 25pt stencil */
#define STENCIL_4( Anext, A0 )                         \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_2LOOPS_4( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)];

#define STENCIL2_2LOOPS_4( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_3LOOPS_4( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)];

#define STENCIL2_3LOOPS_4( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)];

#define STENCIL3_3LOOPS_4( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

/* Vectorized version - SSE complaint */
#define STENCIL_SSE_4( Anext, A0 )             \
          __m128d A0zm4 = _mm_load_pd( &A0[Index3D(nx, ny, i-2, j, k)] ); \
          __m128d A0zce = _mm_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m128d A0zp4 = _mm_load_pd( &A0[Index3D(nx, ny, i+2, j, k)] ); \
\
          __m128d shufv = _mm_shuffle_pd( A0zm4, A0zce, 0x0F ); /* 00001111 */ \
          __m128d A0zm3 = _mm_shuffle_pd( A0zm4, shufv, 0x99 ); /* 10011001 */ \
          __m128d A0zm2 = _mm_shuffle_pd( A0zm4, A0zce, 0x4E ); /* 01001110 */ \
          __m128d A0zm1 = _mm_shuffle_pd( shufv, A0zce, 0x99 ); /* 10011001 */ \
          __m128d shufv = _mm_shuffle_pd( A0zce, A0zp4, 0x0F ); /* 00001111 */ \
          __m128d A0zp1 = _mm_shuffle_pd( A0zce, shufv, 0x99 ); /* 10011001 */ \
          __m128d A0zp2 = _mm_shuffle_pd( A0zce, A0zp4, 0x4E ); /* 01001110 */ \
          __m128d A0zp3 = _mm_shuffle_pd( shufv, A0zp4, 0x99 ); /* 10011001 */ \
  \
          __m128d tempz = _mm_add_pd( _mm_mul_pd( A0zce, c00v ), \
                                      _mm_add_pd( _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0zm1, c1zv ), \
                                                                          _mm_mul_pd( A0zm2, c2zv )), \
                                                              _mm_add_pd( _mm_mul_pd( A0zm3, c3zv ), \
                                                                          _mm_mul_pd( A0zm4, c4zv ))), \
                                                  _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0zp1, c1zv ), \
                                                                          _mm_mul_pd( A0zp2, c2zv )), \
                                                              _mm_add_pd( _mm_mul_pd( A0zp3, c3zv ), \
                                                                          _mm_mul_pd( A0zp4, c4zv ))))); \
\
          __m128d A0xm4 = _mm_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m128d A0xm3 = _mm_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m128d A0xm2 = _mm_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m128d A0xm1 = _mm_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m128d A0xp1 = _mm_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m128d A0xp2 = _mm_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m128d A0xp3 = _mm_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m128d A0xp4 = _mm_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m128d tempx = _mm_add_pd( _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0xm1, c1xv ), \
                                                              _mm_mul_pd( A0xm2, c2xv )), \
                                                  _mm_add_pd( _mm_mul_pd( A0xm3, c3xv ), \
                                                              _mm_mul_pd( A0xm4, c4xv ))), \
                                      _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0xp1, c1xv ), \
                                                              _mm_mul_pd( A0xp2, c2xv )), \
                                                  _mm_add_pd( _mm_mul_pd( A0xp3, c3xv ), \
                                                              _mm_mul_pd( A0xp4, c4xv )))); \
\
          __m128d A0ym4 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m128d A0ym3 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m128d A0ym2 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m128d A0ym1 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m128d A0yp1 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m128d A0yp2 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m128d A0yp3 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m128d A0yp4 = _mm_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m128d tempy = _mm_add_pd( _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0ym1, c1yv ), \
                                                              _mm_mul_pd( A0ym2, c2yv )), \
                                                  _mm_add_pd( _mm_mul_pd( A0ym3, c3yv ), \
                                                              _mm_mul_pd( A0ym4, c4yv ))), \
                                      _mm_add_pd( _mm_add_pd( _mm_mul_pd( A0yp1, c1yv ), \
                                                              _mm_mul_pd( A0yp2, c2yv )), \
                                                  _mm_add_pd( _mm_mul_pd( A0yp3, c3yv ), \
                                                              _mm_mul_pd( A0yp4, c4yv )))); \
\
          __m128d Anextv = _m_add_pd( _m_add_pd( tempz, tempx ), tempy ); \
          _mm_store_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv ); \


/* 25pt stencil (FACTORED) - Vectorized version SSE */
#define STENCIL_FACTORED_SSE_4( Anext, A0 )                \
          __m128 A0zm4 = _mm_load_ps( &A0[Index3D(nx, ny, i-4, j, k)] ); \
          __m128 A0zce = _mm_load_ps( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m128 A0zp4 = _mm_load_ps( &A0[Index3D(nx, ny, i+4, j, k)] ); \
\
          __m128 shufv = _mm_shuffle_ps( A0zm4, A0zce, 0x0F ); /* 00001111 */ \
          __m128 A0zm3 = _mm_shuffle_ps( A0zm4, shufv, 0x99 ); /* 10011001 */ \
          __m128 A0zm2 = _mm_shuffle_ps( A0zm4, A0zce, 0x4E ); /* 01001110 */ \
          __m128 A0zm1 = _mm_shuffle_ps( shufv, A0zce, 0x99 ); /* 10011001 */ \
          __m128 shufv = _mm_shuffle_ps( A0zce, A0zp4, 0x0F ); /* 00001111 */ \
          __m128 A0zp1 = _mm_shuffle_ps( A0zce, shufv, 0x99 ); /* 10011001 */ \
          __m128 A0zp2 = _mm_shuffle_ps( A0zce, A0zp4, 0x4E ); /* 01001110 */ \
          __m128 A0zp3 = _mm_shuffle_ps( shufv, A0zp4, 0x99 ); /* 10011001 */ \
  \
          __m128 tempz = _mm_add_ps( _mm_mul_ps( A0zce, c00v ),                                                  \
                                     _mm_add_ps( _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0zm1, A0zp1 ), c1zv ),     \
                                                             _mm_mul_ps( _mm_add_ps( A0zm2, A0zp2 ), c2zv )),    \
                                                 _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0zm3, A0zp3 ), c3zv ),     \
                                                             _mm_mul_ps( _mm_add_ps( A0zm4, A0zp4 ), c4zv ) ))); \
\
          __m128 A0xm4 = _mm_load_ps( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m128 A0xm3 = _mm_load_ps( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m128 A0xm2 = _mm_load_ps( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m128 A0xm1 = _mm_load_ps( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m128 A0xp1 = _mm_load_ps( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m128 A0xp2 = _mm_load_ps( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m128 A0xp3 = _mm_load_ps( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m128 A0xp4 = _mm_load_ps( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m128 tempx = _mm_add_ps( _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0xm1, A0xp1 ), c1xv ),   \
                                                 _mm_mul_ps( _mm_add_ps( A0xm2, A0xp2 ), c2xv )),  \
                                     _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0xm3, A0xp3 ), c3xv ),   \
                                                 _mm_mul_ps( _mm_add_ps( A0xm4, A0xp4 ), c4xv ))); \
\
          __m128 A0ym4 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m128 A0ym3 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m128 A0ym2 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m128 A0ym1 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m128 A0yp1 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m128 A0yp2 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m128 A0yp3 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m128 A0yp4 = _mm_load_ps( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m128 tempy = _mm_add_ps( _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0ym1, A0yp1 ), c1yv ),   \
                                                 _mm_mul_ps( _mm_add_ps( A0ym2, A0yp2 ), c2yv )),  \
                                     _mm_add_ps( _mm_mul_ps( _mm_add_ps( A0ym3, A0yp3 ), c3yv ),   \
                                                 _mm_mul_ps( _mm_add_ps( A0ym4, A0yp4 ), c4yv ))); \
\
          __m128 Anextv = _m_add_ps( _m_add_ps( tempz, tempx ), tempy ); \
          _mm_store_ps( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


/* Vectorized version - AVX complaint */
#define STENCIL_AVX_4( Anext, A0 )             \
          __m256d A0zm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i-4, j, k)] ); \
          __m256d A0zce = _mm256_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m256d A0zp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i+4, j, k)] ); \
\
          __m256d A0zm2 = _mm256_permute2f128_pd( A0zm4, A0zce, 0x21 ); /* 100001 */ \
          __m256d A0zm3 = _mm256_shuffle_pd( A0zm4, A0zm2, 0x5 );       /* 0101   */ \
          __m256d A0zm1 = _mm256_shuffle_pd( A0zm2, A0zce, 0x5 );       /* 0101   */ \
          __m256d A0zp2 = _mm256_permute2f128_pd( A0zce, A0zp4, 0x21 ); /* 100001 */ \
          __m256d A0zp1 = _mm256_shuffle_pd( A0zce, A0zp2, 0x5 );       /* 0101   */ \
          __m256d A0zp3 = _mm256_shuffle_pd( A0zp2, A0zp4, 0x5 );       /* 0101   */ \
\
          __m256d tempz = _mm256_add_pd( _mm256_mul_pd( A0zce, c00v ), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zm1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zm2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zm3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zm4, c4zv ))),   \
                                                        _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zp1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zp2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zp3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zp4, c4zv ))))); \
\
          __m256d A0xm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m256d A0xm3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m256d A0xm2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m256d A0xm1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m256d A0xp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m256d A0xp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m256d A0xp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m256d A0xp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m256d tempx = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xm1, c1xv ),    \
                                                                       _mm256_mul_pd( A0xm2, c2xv )),   \
                                                        _mm256_add_pd( _mm256_mul_pd( A0xm3, c3xv ),    \
                                                                       _mm256_mul_pd( A0xm4, c4xv ))),  \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xp1, c1xv ),    \
                                                                       _mm256_mul_pd( A0xp2, c2xv )),   \
                                                        _mm256_add_pd( _mm256_mul_pd( A0xp3, c3xv ),    \
                                                                       _mm256_mul_pd( A0xp4, c4xv )))); \
\
          __m256d A0ym4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m256d A0ym3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m256d A0ym2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m256d A0ym1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m256d A0yp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m256d A0yp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m256d A0yp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m256d A0yp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m256d tempy = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0ym1, c1yv ), \
                                                                       _mm256_mul_pd( A0ym2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0ym3, c3yv ), \
                                                                       _mm256_mul_pd( A0ym4, c4yv ))), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0yp1, c1yv ), \
                                                                       _mm256_mul_pd( A0yp2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0yp3, c3yv ), \
                                                                       _mm256_mul_pd( A0yp4, c4yv )))); \
\
          __m256d Anextv = _mm256_add_pd( _mm256_add_pd( tempz, tempx ), tempy ); \
          _mm256_stream_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv ); \


#define STENCIL_FACTORED_AVX_4( Anext, A0 )             \
          __m256d A0zm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i-4, j, k)] ); \
          __m256d A0zce = _mm256_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m256d A0zp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i+4, j, k)] ); \
\
          __m256d shufv = _mm256_shuffle_pd( A0zm4, A0zce, 0x0F ); /* 00001111 */ \
          __m256d A0zm3 = _mm256_shuffle_pd( A0zm4, shufv, 0x99 ); /* 10011001 */ \
          __m256d A0zm2 = _mm256_shuffle_pd( A0zm4, A0zce, 0x4E ); /* 01001110 */ \
          __m256d A0zm1 = _mm256_shuffle_pd( shufv, A0zce, 0x99 ); /* 10011001 */ \
                  shufv = _mm256_shuffle_pd( A0zce, A0zp4, 0x0F ); /* 00001111 */ \
          __m256d A0zp1 = _mm256_shuffle_pd( A0zce, shufv, 0x99 ); /* 10011001 */ \
          __m256d A0zp2 = _mm256_shuffle_pd( A0zce, A0zp4, 0x4E ); /* 01001110 */ \
          __m256d A0zp3 = _mm256_shuffle_pd( shufv, A0zp4, 0x99 ); /* 10011001 */ \
\
          __m256d tempz = _mm256_add_pd( _mm256_mul_pd( A0zce, c00v ),                                                  \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0zm1, A0zp1 ), c1zv ),     \
                                                                       _mm256_mul_pd( _mm256_add_pd( A0zm2, A0zp2 ), c2zv )),    \
                                                        _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0zm3, A0zp3 ), c3zv ),     \
                                                                       _mm256_mul_pd( _mm256_add_pd( A0zm4, A0zp4 ), c4zv ) ))); \
\
          __m256d A0xm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m256d A0xm3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m256d A0xm2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m256d A0xm1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m256d A0xp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m256d A0xp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m256d A0xp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m256d A0xp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m256d tempx = _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0xm1, A0xp1 ), c1xv ),   \
                                                        _mm256_mul_pd( _mm256_add_pd( A0xm2, A0xp2 ), c2xv )),  \
                                         _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0xm3, A0xp3 ), c3xv ),   \
                                                        _mm256_mul_pd( _mm256_add_pd( A0xm4, A0xp4 ), c4xv ))); \
\
          __m256d A0ym4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m256d A0ym3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m256d A0ym2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m256d A0ym1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m256d A0yp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m256d A0yp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m256d A0yp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m256d A0yp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m256d tempy = _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0ym1, A0yp1 ), c1yv ),   \
                                                        _mm256_mul_pd( _mm256_add_pd( A0ym2, A0yp2 ), c2yv )),  \
                                         _mm256_add_pd( _mm256_mul_pd( _mm256_add_pd( A0ym3, A0yp3 ), c3yv ),   \
                                                        _mm256_mul_pd( _mm256_add_pd( A0ym4, A0yp4 ), c4yv ))); \
\
          __m256d Anextv = _mm256_add_pd( _mm256_add_pd( tempz, tempx ), tempy ); \
          _mm256_stream_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL1_2LOOPS_AVX_4( Anext, A0 )                 \
          __m256d A0xp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m256d A0xp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m256d A0xp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m256d A0xp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m256d tempx = _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xp1, c1xv ),    \
                                                        _mm256_mul_pd( A0xp2, c2xv )),   \
                                         _mm256_add_pd( _mm256_mul_pd( A0xp3, c3xv ),    \
                                                        _mm256_mul_pd( A0xp4, c4xv ))); \
\
          __m256d A0ym4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m256d A0ym3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m256d A0ym2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m256d A0ym1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m256d A0yp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m256d A0yp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m256d A0yp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m256d A0yp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m256d tempy = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0ym1, c1yv ), \
                                                                       _mm256_mul_pd( A0ym2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0ym3, c3yv ), \
                                                                       _mm256_mul_pd( A0ym4, c4yv ))), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0yp1, c1yv ), \
                                                                       _mm256_mul_pd( A0yp2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0yp3, c3yv ), \
                                                                       _mm256_mul_pd( A0yp4, c4yv )))); \
\
          __m256d Anextv = _mm256_add_pd( tempx, tempy ); \
          _mm256_store_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL2_2LOOPS_AVX_4( Anext, A0 )                 \
          __m256d A0zm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i-4, j, k)] ); \
          __m256d A0zce = _mm256_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m256d A0zp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i+4, j, k)] ); \
\
          __m256d A0zm2 = _mm256_permute2f128_pd( A0zm4, A0zce, 0x21 ); /* 100001 */ \
          __m256d A0zm3 = _mm256_shuffle_pd( A0zm4, A0zm2, 0x5 );       /* 0101   */ \
          __m256d A0zm1 = _mm256_shuffle_pd( A0zm2, A0zce, 0x5 );       /* 0101   */ \
          __m256d A0zp2 = _mm256_permute2f128_pd( A0zce, A0zp4, 0x21 ); /* 100001 */ \
          __m256d A0zp1 = _mm256_shuffle_pd( A0zce, A0zp2, 0x5 );       /* 0101   */ \
          __m256d A0zp3 = _mm256_shuffle_pd( A0zp2, A0zp4, 0x5 );       /* 0101   */ \
\
          __m256d Anextv = _mm256_load_pd( &Anext[Index3D(nx, ny, i, j, k)] ); \
\
          __m256d tempz = _mm256_add_pd( _mm256_mul_pd( A0zce, c00v ), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zm1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zm2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zm3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zm4, c4zv ))),   \
                                                        _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zp1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zp2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zp3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zp4, c4zv ))))); \
\
          __m256d A0xm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m256d A0xm3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m256d A0xm2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m256d A0xm1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
\
          __m256d tempx = _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xm1, c1xv ),    \
                                                        _mm256_mul_pd( A0xm2, c2xv )),   \
                                         _mm256_add_pd( _mm256_mul_pd( A0xm3, c3xv ),    \
                                                        _mm256_mul_pd( A0xm4, c4xv ))); \
\
          Anextv = _mm256_add_pd( tempz, _mm256_add_pd( tempx, Anextv ) ); \
          _mm256_stream_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL1_3LOOPS_AVX_4( Anext, A0 )                 \
          __m256d A0ym4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m256d A0ym3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m256d A0ym2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m256d A0ym1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m256d A0yp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m256d A0yp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m256d A0yp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m256d A0yp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
          __m256d tempy = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0ym1, c1yv ), \
                                                                       _mm256_mul_pd( A0ym2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0ym3, c3yv ), \
                                                                       _mm256_mul_pd( A0ym4, c4yv ))), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0yp1, c1yv ), \
                                                                       _mm256_mul_pd( A0yp2, c2yv )), \
                                                        _mm256_add_pd( _mm256_mul_pd( A0yp3, c3yv ), \
                                                                       _mm256_mul_pd( A0yp4, c4yv )))); \
\
          _mm256_store_pd( &Anext[Index3D(nx, ny, i, j, k)], tempy );


#define STENCIL2_3LOOPS_AVX_4( Anext, A0 )                 \
          __m256d A0xm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m256d A0xm3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m256d A0xm2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m256d A0xm1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m256d A0xp1 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m256d A0xp2 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m256d A0xp3 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m256d A0xp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m256d Anextv = _mm256_load_pd( &Anext[Index3D(nx, ny, i, j, k)] ); \
\
          __m256d tempx = _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xp1, c1xv ),    \
                                                                       _mm256_mul_pd( A0xp2, c2xv )),   \
                                                        _mm256_add_pd( _mm256_mul_pd( A0xp3, c3xv ),    \
                                                                       _mm256_mul_pd( A0xp4, c4xv ))), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0xm1, c1xv ),    \
                                                                       _mm256_mul_pd( A0xm2, c2xv )),   \
                                                        _mm256_add_pd( _mm256_mul_pd( A0xm3, c3xv ),    \
                                                                       _mm256_mul_pd( A0xm4, c4xv )))); \
\
          Anextv = _mm256_add_pd( tempx, Anextv ); \
          _mm256_store_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL3_3LOOPS_AVX_4( Anext, A0 )                 \
          __m256d A0zm4 = _mm256_load_pd( &A0[Index3D(nx, ny, i-4, j, k)] ); \
          __m256d A0zce = _mm256_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m256d A0zp4 = _mm256_load_pd( &A0[Index3D(nx, ny, i+4, j, k)] ); \
\
          __m256d A0zm2 = _mm256_permute2f128_pd( A0zm4, A0zce, 0x21 ); /* 100001 */ \
          __m256d A0zm3 = _mm256_shuffle_pd( A0zm4, A0zm2, 0x5 );       /* 0101   */ \
          __m256d A0zm1 = _mm256_shuffle_pd( A0zm2, A0zce, 0x5 );       /* 0101   */ \
          __m256d A0zp2 = _mm256_permute2f128_pd( A0zce, A0zp4, 0x21 ); /* 100001 */ \
          __m256d A0zp1 = _mm256_shuffle_pd( A0zce, A0zp2, 0x5 );       /* 0101   */ \
          __m256d A0zp3 = _mm256_shuffle_pd( A0zp2, A0zp4, 0x5 );       /* 0101   */ \
\
          __m256d Anextv = _mm256_load_pd( &Anext[Index3D(nx, ny, i, j, k)] ); \
\
          __m256d tempz = _mm256_add_pd( _mm256_mul_pd( A0zce, c00v ), \
                                         _mm256_add_pd( _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zm1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zm2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zm3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zm4, c4zv ))),   \
                                                        _mm256_add_pd( _mm256_add_pd( _mm256_mul_pd( A0zp1, c1zv ),     \
                                                                                      _mm256_mul_pd( A0zp2, c2zv )),    \
                                                                       _mm256_add_pd( _mm256_mul_pd( A0zp3, c3zv ),     \
                                                                                      _mm256_mul_pd( A0zp4, c4zv ))))); \
\
          Anextv = _mm256_add_pd( tempz, Anextv ); \
          _mm256_stream_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


/* Vectorized version - MIC complaint */
#define STENCIL_MIC_4( Anext, A0 )             \
          __m512d A0zm8 = _mm512_load_pd( &A0[Index3D(nx, ny, i-8, j, k)] ); \
          __m512d A0zce = _mm512_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m512d A0zp8 = _mm512_load_pd( &A0[Index3D(nx, ny, i+8, j, k)] ); \
\
          __m512d A0zm4 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 8  ); \
          __m512d A0zm3 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 10 ); \
          __m512d A0zm2 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 12 ); \
          __m512d A0zm1 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 14 ); \
          __m512d A0zp1 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 2  ); \
          __m512d A0zp2 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 4  ); \
          __m512d A0zp3 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 6  ); \
          __m512d A0zp4 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 8  ); \
\
          __m512d temp1 = _mm512_mul_pd( A0zce, c00v ); \
          __m512d temp2 = _mm512_mul_pd( A0zm1, c1zv ); \
                  temp1 = _mm512_fmadd_pd( A0zm2, c2zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zm3, c3zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zm4, c4zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zp1, c1zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zp2, c2zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zp3, c3zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zp4, c4zv, temp1 ); \
\
          __m512d A0xm4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m512d A0xm3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m512d A0xm2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m512d A0xm1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
          __m512d A0xp1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m512d A0xp2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m512d A0xp3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m512d A0xp4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
                  temp1 = _mm512_fmadd_pd( A0xm1, c1xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xm2, c2xv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0xm3, c3xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xm4, c4xv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0xp1, c1xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xp2, c2xv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0xp3, c3xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xp4, c4xv, temp2 ); \
\
          __m512d A0ym4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m512d A0ym3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m512d A0ym2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m512d A0ym1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m512d A0yp1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m512d A0yp2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m512d A0yp3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m512d A0yp4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
                  temp1 = _mm512_fmadd_pd( A0ym1, c1yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0ym2, c2yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0ym3, c3yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0ym4, c4yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0yp1, c1yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0yp2, c2yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0yp3, c3yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0yp4, c4yv, temp2 ); \
\
          __m512d Anextv = _mm512_add_pd( temp1, temp2 ); \
          _mm512_storenr_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL1_2LOOPS_MIC_4( Anext, A0 )                 \
          __m512d A0xp1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+1, k)] ); \
          __m512d A0xp2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+2, k)] ); \
          __m512d A0xp3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+3, k)] ); \
          __m512d A0xp4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j+4, k)] ); \
\
          __m512d temp1 = _mm512_mul_pd( A0xp1, c1xv ); \
          __m512d temp2 = _mm512_mul_pd( A0xp2, c2xv ); \
                  temp1 = _mm512_fmadd_pd( A0xp3, c3xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xp4, c4xv, temp2 ); \
\
          __m512d A0ym4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-4)] ); \
          __m512d A0ym3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-3)] ); \
          __m512d A0ym2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-2)] ); \
          __m512d A0ym1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k-1)] ); \
          __m512d A0yp1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+1)] ); \
          __m512d A0yp2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+2)] ); \
          __m512d A0yp3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+3)] ); \
          __m512d A0yp4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j, k+4)] ); \
\
                  temp1 = _mm512_fmadd_pd( A0ym1, c1yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0ym2, c2yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0ym3, c3yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0ym4, c4yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0yp1, c1yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0yp2, c2yv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0yp3, c3yv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0yp4, c4yv, temp2 ); \
\
          __m512d Anextv = _mm512_add_pd( temp1, temp2 ); \
          _mm512_store_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


#define STENCIL2_2LOOPS_MIC_4( Anext, A0 )                 \
          __m512d A0zm8 = _mm512_load_pd( &A0[Index3D(nx, ny, i-8, j, k)] ); \
          __m512d A0zce = _mm512_load_pd( &A0[Index3D(nx, ny,   i, j, k)] ); \
          __m512d A0zp8 = _mm512_load_pd( &A0[Index3D(nx, ny, i+8, j, k)] ); \
\
          __m512d A0zm4 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 8  ); \
          __m512d A0zm3 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 10 ); \
          __m512d A0zm2 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 12 ); \
          __m512d A0zm1 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zce, (__m512i)A0zm8, 14 ); \
          __m512d A0zp1 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 2  ); \
          __m512d A0zp2 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 4  ); \
          __m512d A0zp3 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 6  ); \
          __m512d A0zp4 = (__m512d)_mm512_alignr_epi32( (__m512i)A0zp8, (__m512i)A0zce, 8  ); \
\
          __m512d Anextv = _mm512_load_pd( &Anext[Index3D(nx, ny, i, j, k)] ); \
\
          __m512d temp1 = _mm512_mul_pd( A0zce, c00v ); \
          __m512d temp2 = _mm512_mul_pd( A0zm1, c1zv ); \
                  temp1 = _mm512_fmadd_pd( A0zm2, c2zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zm3, c3zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zm4, c4zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zp1, c1zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zp2, c2zv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0zp3, c3zv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0zp4, c4zv, temp1 ); \
\
          __m512d A0xm4 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-4, k)] ); \
          __m512d A0xm3 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-3, k)] ); \
          __m512d A0xm2 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-2, k)] ); \
          __m512d A0xm1 = _mm512_load_pd( &A0[Index3D(nx, ny, i, j-1, k)] ); \
\
                  temp1 = _mm512_fmadd_pd( A0xm1, c1xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xm2, c2xv, temp2 ); \
                  temp1 = _mm512_fmadd_pd( A0xm3, c3xv, temp1 ); \
                  temp2 = _mm512_fmadd_pd( A0xm4, c4xv, temp2 ); \
\
          Anextv = _mm512_add_pd( temp1, _mm512_add_pd( temp2, Anextv ) ); \
          _mm512_storenr_pd( &Anext[Index3D(nx, ny, i, j, k)], Anextv );


/* 25pt stencil (FACTORED) */
#define STENCIL_FACTORED_4( Anext, A0 )                \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * (A0[Index3D (nx, ny, i, j, k + 1)] + A0[Index3D (nx, ny, i, j, k - 1)]) + \
             c2y * (A0[Index3D (nx, ny, i, j, k + 2)] + A0[Index3D (nx, ny, i, j, k - 2)]) + \
             c3y * (A0[Index3D (nx, ny, i, j, k + 3)] + A0[Index3D (nx, ny, i, j, k - 3)]) + \
             c4y * (A0[Index3D (nx, ny, i, j, k + 4)] + A0[Index3D (nx, ny, i, j, k - 4)]) + \
             c1x * (A0[Index3D (nx, ny, i, j + 1, k)] + A0[Index3D (nx, ny, i, j - 1, k)]) + \
             c2x * (A0[Index3D (nx, ny, i, j + 2, k)] + A0[Index3D (nx, ny, i, j - 2, k)]) + \
             c3x * (A0[Index3D (nx, ny, i, j + 3, k)] + A0[Index3D (nx, ny, i, j - 3, k)]) + \
             c4x * (A0[Index3D (nx, ny, i, j + 4, k)] + A0[Index3D (nx, ny, i, j - 4, k)]) + \
             c1z * (A0[Index3D (nx, ny, i + 1, j, k)] + A0[Index3D (nx, ny, i - 1, j, k)]) + \
             c2z * (A0[Index3D (nx, ny, i + 2, j, k)] + A0[Index3D (nx, ny, i - 2, j, k)]) + \
             c3z * (A0[Index3D (nx, ny, i + 3, j, k)] + A0[Index3D (nx, ny, i - 3, j, k)]) + \
             c4z * (A0[Index3D (nx, ny, i + 4, j, k)] + A0[Index3D (nx, ny, i - 4, j, k)]) + \
             c00 * A0[Index3D (nx, ny, i, j, k)];


/* 43pt stencil */
#define STENCIL_7( Anext, A0 )                         \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k + 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k + 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k + 7)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k - 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k - 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k - 7)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j + 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j + 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j + 7, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j - 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j - 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j - 7, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i + 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i + 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i + 7, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i - 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i - 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i - 7, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_2LOOPS_7( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k + 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k + 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k + 7)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k - 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k - 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k - 7)] + \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j + 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j + 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j + 7, k)];

#define STENCIL2_2LOOPS_7( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j - 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j - 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j - 7, k)] + \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i + 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i + 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i + 7, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i - 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i - 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i - 7, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_3LOOPS_7( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * A0[Index3D (nx, ny, i, j, k + 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k + 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k + 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k + 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k + 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k + 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k + 7)] + \
             c1y * A0[Index3D (nx, ny, i, j, k - 1)] + \
             c2y * A0[Index3D (nx, ny, i, j, k - 2)] + \
             c3y * A0[Index3D (nx, ny, i, j, k - 3)] + \
             c4y * A0[Index3D (nx, ny, i, j, k - 4)] + \
             c5y * A0[Index3D (nx, ny, i, j, k - 5)] + \
             c6y * A0[Index3D (nx, ny, i, j, k - 6)] + \
             c7y * A0[Index3D (nx, ny, i, j, k - 7)];

#define STENCIL2_3LOOPS_7( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1x * A0[Index3D (nx, ny, i, j + 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j + 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j + 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j + 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j + 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j + 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j + 7, k)] + \
             c1x * A0[Index3D (nx, ny, i, j - 1, k)] + \
             c2x * A0[Index3D (nx, ny, i, j - 2, k)] + \
             c3x * A0[Index3D (nx, ny, i, j - 3, k)] + \
             c4x * A0[Index3D (nx, ny, i, j - 4, k)] + \
             c5x * A0[Index3D (nx, ny, i, j - 5, k)] + \
             c6x * A0[Index3D (nx, ny, i, j - 6, k)] + \
             c7x * A0[Index3D (nx, ny, i, j - 7, k)];

#define STENCIL3_3LOOPS_7( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] +=          \
             c1z * A0[Index3D (nx, ny, i + 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i + 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i + 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i + 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i + 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i + 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i + 7, j, k)] + \
             c1z * A0[Index3D (nx, ny, i - 1, j, k)] + \
             c2z * A0[Index3D (nx, ny, i - 2, j, k)] + \
             c3z * A0[Index3D (nx, ny, i - 3, j, k)] + \
             c4z * A0[Index3D (nx, ny, i - 4, j, k)] + \
             c5z * A0[Index3D (nx, ny, i - 5, j, k)] + \
             c6z * A0[Index3D (nx, ny, i - 6, j, k)] + \
             c7z * A0[Index3D (nx, ny, i - 7, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

/* 43pt stencil (FACTORED) */
#define STENCIL_FACTORED_7( Anext, A0 )                \
          Anext[Index3D (nx, ny, i, j, k)] =           \
             c1y * (A0[Index3D (nx, ny, i, j, k + 1)] + A0[Index3D (nx, ny, i, j, k - 1)]) + \
             c2y * (A0[Index3D (nx, ny, i, j, k + 2)] + A0[Index3D (nx, ny, i, j, k - 2)]) + \
             c3y * (A0[Index3D (nx, ny, i, j, k + 3)] + A0[Index3D (nx, ny, i, j, k - 3)]) + \
             c4y * (A0[Index3D (nx, ny, i, j, k + 4)] + A0[Index3D (nx, ny, i, j, k - 4)]) + \
             c5y * (A0[Index3D (nx, ny, i, j, k + 5)] + A0[Index3D (nx, ny, i, j, k - 5)]) + \
             c6y * (A0[Index3D (nx, ny, i, j, k + 6)] + A0[Index3D (nx, ny, i, j, k - 6)]) + \
             c7y * (A0[Index3D (nx, ny, i, j, k + 7)] + A0[Index3D (nx, ny, i, j, k - 7)]) + \
             c1x * (A0[Index3D (nx, ny, i, j + 1, k)] + A0[Index3D (nx, ny, i, j - 1, k)]) + \
             c2x * (A0[Index3D (nx, ny, i, j + 2, k)] + A0[Index3D (nx, ny, i, j - 2, k)]) + \
             c3x * (A0[Index3D (nx, ny, i, j + 3, k)] + A0[Index3D (nx, ny, i, j - 3, k)]) + \
             c4x * (A0[Index3D (nx, ny, i, j + 4, k)] + A0[Index3D (nx, ny, i, j - 4, k)]) + \
             c5x * (A0[Index3D (nx, ny, i, j + 5, k)] + A0[Index3D (nx, ny, i, j - 5, k)]) + \
             c6x * (A0[Index3D (nx, ny, i, j + 6, k)] + A0[Index3D (nx, ny, i, j - 6, k)]) + \
             c7x * (A0[Index3D (nx, ny, i, j + 7, k)] + A0[Index3D (nx, ny, i, j - 7, k)]) + \
             c1z * (A0[Index3D (nx, ny, i + 1, j, k)] + A0[Index3D (nx, ny, i - 1, j, k)]) + \
             c2z * (A0[Index3D (nx, ny, i + 2, j, k)] + A0[Index3D (nx, ny, i - 2, j, k)]) + \
             c3z * (A0[Index3D (nx, ny, i + 3, j, k)] + A0[Index3D (nx, ny, i - 3, j, k)]) + \
             c4z * (A0[Index3D (nx, ny, i + 4, j, k)] + A0[Index3D (nx, ny, i - 4, j, k)]) + \
             c5z * (A0[Index3D (nx, ny, i + 5, j, k)] + A0[Index3D (nx, ny, i - 5, j, k)]) + \
             c6z * (A0[Index3D (nx, ny, i + 6, j, k)] + A0[Index3D (nx, ny, i - 6, j, k)]) + \
             c7z * (A0[Index3D (nx, ny, i + 7, j, k)] + A0[Index3D (nx, ny, i - 7, j, k)]) + \
             c00 * A0[Index3D (nx, ny, i, j, k)];


/* 85pt stencil */
#define STENCIL_14( Anext, A0 )                          \
          Anext[Index3D (nx, ny, i, j, k)] =             \
             c1y  * A0[Index3D (nx, ny, i, j, k + 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k + 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k + 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k + 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k + 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k + 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k + 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k + 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k + 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k + 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k + 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k + 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k + 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k + 14)] + \
             c1y  * A0[Index3D (nx, ny, i, j, k - 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k - 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k - 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k - 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k - 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k - 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k - 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k - 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k - 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k - 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k - 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k - 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k - 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k - 14)] + \
             c1x  * A0[Index3D (nx, ny, i, j + 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j + 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j + 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j + 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j + 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j + 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j + 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j + 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j + 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j + 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j + 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j + 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j + 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j + 14, k)] + \
             c1x  * A0[Index3D (nx, ny, i, j - 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j - 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j - 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j - 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j - 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j - 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j - 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j - 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j - 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j - 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j - 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j - 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j - 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j - 14, k)] + \
             c1z  * A0[Index3D (nx, ny, i + 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i + 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i + 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i + 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i + 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i + 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i + 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i + 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i + 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i + 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i + 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i + 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i + 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i + 14, j, k)] + \
             c1z  * A0[Index3D (nx, ny, i - 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i - 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i - 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i - 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i - 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i - 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i - 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i - 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i - 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i - 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i - 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i - 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i - 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i - 14, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_2LOOPS_14( Anext, A0 )                  \
          Anext[Index3D (nx, ny, i, j, k)] =             \
             c1y  * A0[Index3D (nx, ny, i, j, k + 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k + 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k + 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k + 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k + 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k + 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k + 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k + 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k + 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k + 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k + 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k + 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k + 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k + 14)] + \
             c1y  * A0[Index3D (nx, ny, i, j, k - 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k - 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k - 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k - 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k - 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k - 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k - 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k - 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k - 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k - 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k - 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k - 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k - 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k - 14)] + \
             c1x  * A0[Index3D (nx, ny, i, j + 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j + 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j + 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j + 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j + 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j + 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j + 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j + 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j + 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j + 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j + 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j + 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j + 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j + 14, k)];

#define STENCIL2_2LOOPS_14( Anext, A0 )                  \
          Anext[Index3D (nx, ny, i, j, k)] +=            \
             c1x  * A0[Index3D (nx, ny, i, j - 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j - 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j - 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j - 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j - 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j - 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j - 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j - 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j - 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j - 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j - 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j - 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j - 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j - 14, k)] + \
             c1z  * A0[Index3D (nx, ny, i + 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i + 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i + 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i + 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i + 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i + 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i + 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i + 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i + 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i + 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i + 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i + 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i + 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i + 14, j, k)] + \
             c1z  * A0[Index3D (nx, ny, i - 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i - 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i - 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i - 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i - 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i - 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i - 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i - 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i - 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i - 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i - 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i - 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i - 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i - 14, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

#define STENCIL1_3LOOPS_14( Anext, A0 )                  \
          Anext[Index3D (nx, ny, i, j, k)] =             \
             c1y  * A0[Index3D (nx, ny, i, j, k + 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k + 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k + 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k + 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k + 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k + 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k + 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k + 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k + 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k + 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k + 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k + 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k + 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k + 14)] + \
             c1y  * A0[Index3D (nx, ny, i, j, k - 1)] +  \
             c2y  * A0[Index3D (nx, ny, i, j, k - 2)] +  \
             c3y  * A0[Index3D (nx, ny, i, j, k - 3)] +  \
             c4y  * A0[Index3D (nx, ny, i, j, k - 4)] +  \
             c5y  * A0[Index3D (nx, ny, i, j, k - 5)] +  \
             c6y  * A0[Index3D (nx, ny, i, j, k - 6)] +  \
             c7y  * A0[Index3D (nx, ny, i, j, k - 7)] +  \
             c8y  * A0[Index3D (nx, ny, i, j, k - 8)] +  \
             c9y  * A0[Index3D (nx, ny, i, j, k - 9)] +  \
             c10y * A0[Index3D (nx, ny, i, j, k - 10)] + \
             c11y * A0[Index3D (nx, ny, i, j, k - 11)] + \
             c12y * A0[Index3D (nx, ny, i, j, k - 12)] + \
             c13y * A0[Index3D (nx, ny, i, j, k - 13)] + \
             c14y * A0[Index3D (nx, ny, i, j, k - 14)];

#define STENCIL2_3LOOPS_14( Anext, A0 )                  \
          Anext[Index3D (nx, ny, i, j, k)] +=            \
             c1x  * A0[Index3D (nx, ny, i, j + 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j + 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j + 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j + 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j + 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j + 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j + 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j + 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j + 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j + 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j + 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j + 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j + 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j + 14, k)] + \
             c1x  * A0[Index3D (nx, ny, i, j - 1, k)] +  \
             c2x  * A0[Index3D (nx, ny, i, j - 2, k)] +  \
             c3x  * A0[Index3D (nx, ny, i, j - 3, k)] +  \
             c4x  * A0[Index3D (nx, ny, i, j - 4, k)] +  \
             c5x  * A0[Index3D (nx, ny, i, j - 5, k)] +  \
             c6x  * A0[Index3D (nx, ny, i, j - 6, k)] +  \
             c7x  * A0[Index3D (nx, ny, i, j - 7, k)] +  \
             c8x  * A0[Index3D (nx, ny, i, j - 8, k)] +  \
             c9x  * A0[Index3D (nx, ny, i, j - 9, k)] +  \
             c10x * A0[Index3D (nx, ny, i, j - 10, k)] + \
             c11x * A0[Index3D (nx, ny, i, j - 11, k)] + \
             c12x * A0[Index3D (nx, ny, i, j - 12, k)] + \
             c13x * A0[Index3D (nx, ny, i, j - 13, k)] + \
             c14x * A0[Index3D (nx, ny, i, j - 14, k)];

#define STENCIL3_3LOOPS_14( Anext, A0 )                  \
          Anext[Index3D (nx, ny, i, j, k)] +=            \
             c1z  * A0[Index3D (nx, ny, i + 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i + 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i + 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i + 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i + 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i + 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i + 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i + 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i + 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i + 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i + 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i + 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i + 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i + 14, j, k)] + \
             c1z  * A0[Index3D (nx, ny, i - 1, j, k)] +  \
             c2z  * A0[Index3D (nx, ny, i - 2, j, k)] +  \
             c3z  * A0[Index3D (nx, ny, i - 3, j, k)] +  \
             c4z  * A0[Index3D (nx, ny, i - 4, j, k)] +  \
             c5z  * A0[Index3D (nx, ny, i - 5, j, k)] +  \
             c6z  * A0[Index3D (nx, ny, i - 6, j, k)] +  \
             c7z  * A0[Index3D (nx, ny, i - 7, j, k)] +  \
             c8z  * A0[Index3D (nx, ny, i - 8, j, k)] +  \
             c9z  * A0[Index3D (nx, ny, i - 9, j, k)] +  \
             c10z * A0[Index3D (nx, ny, i - 10, j, k)] + \
             c11z * A0[Index3D (nx, ny, i - 11, j, k)] + \
             c12z * A0[Index3D (nx, ny, i - 12, j, k)] + \
             c13z * A0[Index3D (nx, ny, i - 13, j, k)] + \
             c14z * A0[Index3D (nx, ny, i - 14, j, k)] + \
             c00 * A0[Index3D (nx, ny, i, j, k)];

/* 85pt stencil (FACTORED) */
#define STENCIL_FACTORED_14( Anext, A0 )                 \
          Anext[Index3D (nx, ny, i, j, k)] =             \
             c1y  * (A0[Index3D (nx, ny, i, j, k + 1)] +  A0[Index3D (nx, ny, i, j, k - 1)]) +  \
             c2y  * (A0[Index3D (nx, ny, i, j, k + 2)] +  A0[Index3D (nx, ny, i, j, k - 2)]) +  \
             c3y  * (A0[Index3D (nx, ny, i, j, k + 3)] +  A0[Index3D (nx, ny, i, j, k - 3)]) +  \
             c4y  * (A0[Index3D (nx, ny, i, j, k + 4)] +  A0[Index3D (nx, ny, i, j, k - 4)]) +  \
             c5y  * (A0[Index3D (nx, ny, i, j, k + 5)] +  A0[Index3D (nx, ny, i, j, k - 5)]) +  \
             c6y  * (A0[Index3D (nx, ny, i, j, k + 6)] +  A0[Index3D (nx, ny, i, j, k - 6)]) +  \
             c7y  * (A0[Index3D (nx, ny, i, j, k + 7)] +  A0[Index3D (nx, ny, i, j, k - 7)]) +  \
             c8y  * (A0[Index3D (nx, ny, i, j, k + 8)] +  A0[Index3D (nx, ny, i, j, k - 8)]) +  \
             c9y  * (A0[Index3D (nx, ny, i, j, k + 9)] +  A0[Index3D (nx, ny, i, j, k - 9)]) +  \
             c10y * (A0[Index3D (nx, ny, i, j, k + 10)] + A0[Index3D (nx, ny, i, j, k - 10)]) + \
             c11y * (A0[Index3D (nx, ny, i, j, k + 11)] + A0[Index3D (nx, ny, i, j, k - 11)]) + \
             c12y * (A0[Index3D (nx, ny, i, j, k + 12)] + A0[Index3D (nx, ny, i, j, k - 12)]) + \
             c13y * (A0[Index3D (nx, ny, i, j, k + 13)] + A0[Index3D (nx, ny, i, j, k - 13)]) + \
             c14y * (A0[Index3D (nx, ny, i, j, k + 14)] + A0[Index3D (nx, ny, i, j, k - 14)]) + \
             c1x  * (A0[Index3D (nx, ny, i, j + 1, k)] +  A0[Index3D (nx, ny, i, j - 1, k)]) +  \
             c2x  * (A0[Index3D (nx, ny, i, j + 2, k)] +  A0[Index3D (nx, ny, i, j - 2, k)]) +  \
             c3x  * (A0[Index3D (nx, ny, i, j + 3, k)] +  A0[Index3D (nx, ny, i, j - 3, k)]) +  \
             c4x  * (A0[Index3D (nx, ny, i, j + 4, k)] +  A0[Index3D (nx, ny, i, j - 4, k)]) +  \
             c5x  * (A0[Index3D (nx, ny, i, j + 5, k)] +  A0[Index3D (nx, ny, i, j - 5, k)]) +  \
             c6x  * (A0[Index3D (nx, ny, i, j + 6, k)] +  A0[Index3D (nx, ny, i, j - 6, k)]) +  \
             c7x  * (A0[Index3D (nx, ny, i, j + 7, k)] +  A0[Index3D (nx, ny, i, j - 7, k)]) +  \
             c8x  * (A0[Index3D (nx, ny, i, j + 8, k)] +  A0[Index3D (nx, ny, i, j - 8, k)]) +  \
             c9x  * (A0[Index3D (nx, ny, i, j + 9, k)] +  A0[Index3D (nx, ny, i, j - 9, k)]) +  \
             c10x * (A0[Index3D (nx, ny, i, j + 10, k)] + A0[Index3D (nx, ny, i, j - 10, k)]) + \
             c11x * (A0[Index3D (nx, ny, i, j + 11, k)] + A0[Index3D (nx, ny, i, j - 11, k)]) + \
             c12x * (A0[Index3D (nx, ny, i, j + 12, k)] + A0[Index3D (nx, ny, i, j - 12, k)]) + \
             c13x * (A0[Index3D (nx, ny, i, j + 13, k)] + A0[Index3D (nx, ny, i, j - 13, k)]) + \
             c14x * (A0[Index3D (nx, ny, i, j + 14, k)] + A0[Index3D (nx, ny, i, j - 14, k)]) + \
             c1z  * (A0[Index3D (nx, ny, i + 1, j, k)] +  A0[Index3D (nx, ny, i - 1, j, k)]) +  \
             c2z  * (A0[Index3D (nx, ny, i + 2, j, k)] +  A0[Index3D (nx, ny, i - 2, j, k)]) +  \
             c3z  * (A0[Index3D (nx, ny, i + 3, j, k)] +  A0[Index3D (nx, ny, i - 3, j, k)]) +  \
             c4z  * (A0[Index3D (nx, ny, i + 4, j, k)] +  A0[Index3D (nx, ny, i - 4, j, k)]) +  \
             c5z  * (A0[Index3D (nx, ny, i + 5, j, k)] +  A0[Index3D (nx, ny, i - 5, j, k)]) +  \
             c6z  * (A0[Index3D (nx, ny, i + 6, j, k)] +  A0[Index3D (nx, ny, i - 6, j, k)]) +  \
             c7z  * (A0[Index3D (nx, ny, i + 7, j, k)] +  A0[Index3D (nx, ny, i - 7, j, k)]) +  \
             c8z  * (A0[Index3D (nx, ny, i + 8, j, k)] +  A0[Index3D (nx, ny, i - 8, j, k)]) +  \
             c9z  * (A0[Index3D (nx, ny, i + 9, j, k)] +  A0[Index3D (nx, ny, i - 9, j, k)]) +  \
             c10z * (A0[Index3D (nx, ny, i + 10, j, k)] + A0[Index3D (nx, ny, i - 10, j, k)]) + \
             c11z * (A0[Index3D (nx, ny, i + 11, j, k)] + A0[Index3D (nx, ny, i - 11, j, k)]) + \
             c12z * (A0[Index3D (nx, ny, i + 12, j, k)] + A0[Index3D (nx, ny, i - 12, j, k)]) + \
             c13z * (A0[Index3D (nx, ny, i + 13, j, k)] + A0[Index3D (nx, ny, i - 13, j, k)]) + \
             c14z * (A0[Index3D (nx, ny, i + 14, j, k)] + A0[Index3D (nx, ny, i - 14, j, k)]) + \
             c00 * A0[Index3D (nx, ny, i, j, k)];


/* To enable vector aligned versions of
 * stencil codes they must be vectorized
 * manually. For nontemporal stores, writes
 * must be aligned and therefore is not possible
 * for an unaligned vectorized version of the code
 *
 * A future handmade vectorized version of the
 * stencilprobe might solve all these issues
 */
//_Pragma("vector aligned") \
//_Pragma("vector nontemporal (Anext)") \
//_Pragma("simd") \

/* Stencil body for naive with 1 loop */
#define NAIVE_1LOOP( Anext, A0, LENGTH,            \
                     xs, xe, ys, ye, zs, ze )      \
    for (k = zs; k < ze; k++) {                    \
      for (j = ys; j < ye; j++) {                  \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL( Anext, A0, LENGTH )             \
        }                                          \
      }                                            \
    }

/* Stencil body for naive fission with 2 loops */
#define NAIVE_FISSION_2LOOPS( Anext, A0, LENGTH,   \
                         xs, xe, ys, ye, zs, ze )  \
    for (k = zs; k < ze; k++) {                    \
      for (j = ys; j < ye; j++) {                  \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL1_2LOOPS( Anext, A0, LENGTH )     \
        }                                          \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL2_2LOOPS( Anext, A0, LENGTH )     \
        }                                          \
      }                                            \
    }

/* Stencil body for naive fission with 2 loops */
#define NAIVE_FISSION_3LOOPS( Anext, A0, LENGTH,   \
                         xs, xe, ys, ye, zs, ze )  \
    for (k = zs; k < ze; k++) {                    \
      for (j = ys; j < ye; j++) {                  \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL1_3LOOPS( Anext, A0, LENGTH )     \
        }                                          \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL2_3LOOPS( Anext, A0, LENGTH )     \
        }                                          \
_Pragma("ivdep") \
_Pragma("vector always") \
        for (i = xs; i < xe; i+=IINC) {            \
          STENCIL3_3LOOPS( Anext, A0, LENGTH )     \
        }                                          \
      }                                            \
    }


// Add as many entries as versions
#if defined(FISSION_2LOOPS)
# define NAIVE_( ... ) NAIVE_FISSION_2LOOPS( __VA_ARGS__ )
#elif defined(FISSION_3LOOPS)
# define NAIVE_( ... ) NAIVE_FISSION_3LOOPS( __VA_ARGS__ )
#else
# define NAIVE_( ... ) NAIVE_1LOOP( __VA_ARGS__ )
#endif


#endif /* _PROBE_H_ */

