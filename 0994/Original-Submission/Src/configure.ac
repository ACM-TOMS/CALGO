#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([mindistance], [0.59], [gquintan@uji.es])
AM_INIT_AUTOMAKE(-Wall -Werror foreign)
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h stdint.h stdlib.h])

# Check and add appropriate OpenMP flags.
AC_OPENMP
AC_SUBST(AM_CFLAGS,"$OPENMP_CFLAGS")

AC_ARG_ENABLE(vectorization,
	     AS_HELP_STRING([--enable-vectorization[=type] Specify the vectorization type:
                      yes (default): Automatically detect vectorization support,
                      sse: Add support for SSE,
                      avx2: Add support for AVX2,
                      avx512: Add support for AVX512]),
 [ if   test "x$enableval" = "xyes" ||  test "x$enableval" = "xauto"  ; then
      vectorization=auto
      
   elif test "x$enableval" = "xsse" ; then
      vectorization=sse2
      AC_DEFINE(USE_VECTORIZATION, 1, "Vectorization was requested by user")
      AC_DEFINE(USE_SSE, 1, "SSE vectorization requested by user")
   elif test "x$enableval" = "xavx2" ; then
      vectorization=avx2
      AC_DEFINE(USE_VECTORIZATION, 1, "Vectorization was requested by user")
      AC_DEFINE(USE_AVX2, 1, "AVX2 vectorization requested by user")
   elif test "x$enableval" = "xavx512" ; then
      vectorization=avx512
      AC_DEFINE(USE_VECTORIZATION, 1, "Vectorization was requested by user")
      AC_DEFINE(USE_AVX512, 1, "AVX512 vectorization requested by user")
   else
      AC_MSG_ERROR(["Unknown vectorization type"])
   fi
 ])

#AS_IF([test "x$enable_vectorization" != "xyes" || test "x$enable_vectorization" == "xauto"  ], [
AS_IF([ test "x$enable_vectorization" == "xyes" ], [
  dnl Do the stuff needed for enabling the feature

  # Check and add appropriate vectorization flags.
  #AX_EXT
  AX_CHECK_X86_FEATURES

  AC_CHECK_HEADERS([x86intrin.h])

  #AC_DEFINE(USE_VECTORIZATION, 1, "Automatic vectorization support was requested by user")

  # Get best between AVX512, AVX2, SSE.
      AX_GCC_X86_CPU_SUPPORTS(avx512f,
        [
  		AC_DEFINE(USE_VECTORIZATION, 1, "Automatic vectorization support was requested by user")
  		AC_DEFINE(USE_AVX512, 1, "Automatic vectorization support was automatically detected and set to AVX512")
        ],
	# If not, check AVX2
        [
      		AX_GCC_X86_CPU_SUPPORTS(avx2,
        		[
  				AC_DEFINE(USE_VECTORIZATION, 1, "Automatic vectorization support was requested by user")
  				AC_DEFINE(USE_AVX2, 1, "Automatic vectorization support was automatically detected and set to AVX2")
        		],
			# If not, check SSE
			[
      				AX_GCC_X86_CPU_SUPPORTS(sse2,
        				[
  							AC_DEFINE(USE_VECTORIZATION, 1, "Automatic vectorization support was requested by user")
  							AC_DEFINE(USE_SSE, 1, "Automatic vectorization support was automatically detected and set to SSE2") 
        				] )
			] )
	] )
],

[
AX_CHECK_X86_FEATURES
AC_CHECK_HEADERS([x86intrin.h])
]

)

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC

AC_CONFIG_FILES([Makefile src/Makefile matrices/Makefile test/Makefile])
AC_OUTPUT
