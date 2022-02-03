#if 0
    $Id: netcdf-cpp.in,v 1.5 1997/05/14 15:17:24 steve Exp $ 

    NB: Do not use comments in this file because it has to be processed by
    both FORTRAN and C compilers.
#endif


#ifndef UD_NETCDF_CPP_INC
#define UD_NETCDF_CPP_INC


#if 0
    FORTRAN data type of <data> or <value> argument of corresponding netCDF
    version 3 functions:

    If a macro is not defined, then the corresponding type does not exist
    (e.g. if <NF_BYTE_T> is not defined, then the FORTRAN language does not 
    have a signed, 8-bit data type).
#endif

#if 0
    The FORTRAN data type corresponding to a C "signed char" (e.g. 
    INTEGER*1).
#endif
#define	NF_BYTE_T byte

#if 0
    The FORTRAN data type corresponding to a C "short" (e.g. INTEGER*2).
#endif
#define	NF_SHORT_T integer*2

#if 0
    The FORTRAN data type corresponding to a C "int" (e.g. INTEGER).
#endif
#define	NF_INT_T integer

#if 0
    The FORTRAN data type corresponding to a C "float" (e.g. REAL).
#endif
#define	NF_FLOAT_T real*4

#if 0
    The FORTRAN data type corresponding to a C "double" (e.g. REAL).
#endif
#define	NF_DOUBLE_T real*8

#if 0
    Whether the system uses something besides the IEEE floating-point
    format to represent floating-point values.
#endif
/* #undef NO_IEEE_FLOAT */

#if 0
    FORTRAN data types corresponding to netCDF version 2 data types,
    e.g. INTEGER*1, INTEGER*2, INTEGER, REAL, DOUBLEPRECISION.  If your
    system does not have a corresponding type, then leave it undefined.
#endif
#ifndef NO_NETCDF_2
#   define	NCBYTE_T byte
#   define	NCSHORT_T integer*2
#   define	NCLONG_T integer*4
#   define	NCFLOAT_T real*4
#   define	NCDOUBLE_T real*8
#endif	/* NO_NETCDF_2 not defined */


#if 0
	END OF CUSTOMIZATION
#endif


#ifdef	NF_BYTE_T
#   define NF_INT1_T	NF_BYTE_T
#else
#   define FORTRAN_HAS_NO_BYTE
#endif

#ifdef	NF_SHORT_T
#   define NF_INT2_T	NF_SHORT_T
#else
#   define FORTRAN_HAS_NO_SHORT
#endif

#ifdef	NF_FLOAT_T
#   define NF_REAL_T	NF_FLOAT_T
#else
#   define FORTRAN_HAS_NO_FLOAT
#endif


#endif	/* header-file lockout */
