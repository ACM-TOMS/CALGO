#if !defined(GAMPSI_H)		/* file ignored if included more than once */
#define GAMPSI_H

/***********************************************************************

This header file for programs written in ANSI/ISO Standard C (1989) or
Standard C++ (1998) defines an interface from those languages to the
Fortran routines for Gamma(x) and psi(x), as described in the article

    James S. Ball and Nelson H. F. Beebe,
    ``Algorithm xxx: Quadruple-Precision $\Gamma(x)$ and $\psi(x)$
    Functions for Real Arguments'',
    ACM Transactions on Mathematical Software 26(??), ??--?? (2000).

On UNIX and POSIX-compliant systems, user code should be compiled with
-I/usr/local/include, and linked with -L/usr/local/lib -lgampsi,
assuming the GNU standard installation locations.  In most systems,
additional -l switches will be needed to supply the Fortran library
names.

[01-Aug-2000]
***********************************************************************/

#if defined(__cplusplus)
extern "C"
{
#endif

/* Handle mapping of Fortran names to their equivalents in C */
#if defined(ardent) || defined(fortran_names_uppercase)
    /* Stardent (now defunct) uppercased Fortran names */
#define aeps   	AEPS
#define afpmax 	AFPMAX
#define ainf   	AINF
#define airan  	AIRAN
#define algam	ALGAM
#define alog2  	ALOG2
#define anan   	ANAN
#define astore 	ASTORE
#define astorf 	ASTORF
#define deps   	DEPS
#define derbit 	DERBIT
#define dfloat 	DFLOAT
#define dfpmax 	DFPMAX
#define dgamma 	DGAMMA
#define dinf   	DINF
#define diran  	DIRAN
#define dlgam  	DLGAM
#define dlog2  	DLOG2
#define dnan   	DNAN
#define dpsi   	DPSI
#define dpsiln 	DPSILN
#define dpsum  	DPSUM
#define dran   	DRAN
#define dstore 	DSTORE
#define dstorf 	DSTORF
#define ffalse	FFALSE
#define ftrue	FTRUE
#define gamma  	GAMMA
#define iceil 	ICEIL
#define idceil 	IDCEIL
#define iqceil 	IQCEIL
#define isainf 	ISAINF
#define isanan 	ISANAN
#define isdinf 	ISDINF
#define isdnan 	ISDNAN
#define isqinf 	ISQINF
#define isqnan 	ISQNAN
#define psi    	PSI
#define psiln  	PSILN
#define pythag 	PYTHAG
#define qeps   	QEPS
#define qepsln 	QEPSLN
#define qfloat 	QFLOAT
#define qfpmax 	QFPMAX
#define qgamma 	QGAMMA
#define qinf   	QINF
#define qiran  	QIRAN
#define qlgam  	QLGAM
#define qlog2  	QLOG2
#define qnan   	QNAN
#define qpsi   	QPSI
#define qpsiln 	QPSILN
#define qstore 	QSTORE
#define qstorf 	QSTORF
#define ran    	RAN
#define trapit 	TRAPIT
#elif defined(_AIX) || defined(__hpux) || defined(fortran_names_same)
    /* IBM RS/6000 AIX and HP HP-UX use identical names in C and Fortran */
#else
    /* Everyone else adds a trailing underscore to Fortran names */
#define aeps   	aeps_
#define afpmax 	afpmax_
#define ainf   	ainf_
#define airan  	airan_
#define algam	algam_
#define alog2  	alog2_
#define anan   	anan_
#define astore 	astore_
#define astorf 	astorf_
#define deps   	deps_
#define derbit 	derbit_
#define dfloat 	dfloat_
#define dfpmax 	dfpmax_
#define dgamma 	dgamma_
#define dinf   	dinf_
#define diran  	diran_
#define dlgam  	dlgam_
#define dlog2  	dlog2_
#define dnan   	dnan_
#define dpsi   	dpsi_
#define dpsiln 	dpsiln_
#define dpsum  	dpsum_
#define dran   	dran_
#define dstore 	dstore_
#define dstorf 	dstorf_
#define ffalse	ffalse_
#define ftrue	ftrue_
#define gamma  	gamma_
#define iceil 	iceil_
#define idceil 	idceil_
#define iqceil 	iqceil_
#define isainf 	isainf_
#define isanan 	isanan_
#define isdinf 	isdinf_
#define isdnan 	isdnan_
#define isqinf 	isqinf_
#define isqnan 	isqnan_
#define psi    	psi_
#define psiln	psiln_
#define pythag 	pythag_
#define qeps   	qeps_
#define qepsln 	qepsln_
#define qfloat 	qfloat_
#define qfpmax 	qfpmax_
#define qgamma 	qgamma_
#define qinf   	qinf_
#define qiran  	qiran_
#define qlgam  	qlgam_
#define qlog2  	qlog2_
#define qnan   	qnan_
#define qpsi   	qpsi_
#define qpsiln	qpsiln_
#define qstore 	qstore_
#define qstorf 	qstorf_
#define ran    	ran_
#define trapit 	trapit_
#endif

/* Handle mappings of Fortran datatypes to C/C++ equivalents.  We
   define these with macros, rather than typedefs, in order to be able
   to test for user-defined overriding values. */

#if !defined(fortran_double_precision)
#define fortran_double_precision double
#endif

#if !defined(fortran_integer)
#define fortran_integer int
#endif

#if !defined(fortran_logical)
    /*******************************************************************
       NB: Fortran implementations differ (sometimes even on the same
       O/S!) about the representation of .TRUE. and .FALSE..

       Most UNIX Fortran compilers produce 0 for .FALSE. and 1 for
       .TRUE., but a few (Compaq/DEC f77, f90, f95; PGI pgf77, pgf90,
       pghpf) produce 0 for .FALSE. and -1 for .TRUE..

       In C and C++, these conform to the zero-for-false,
       nonzero-for-true, language convention.

       Because it may sometimes be necessary to pass a LOGICAL
       constant from a C function to a Fortran routine, the
       representation matters.  To avoid system dependence, use the
       functions ftrue() and ffalse() to represent .TRUE. and
       .FALSE. in such calls.
    *******************************************************************/
#define fortran_logical int
#endif

#if !defined(fortran_quadruple_precision)
#define fortran_quadruple_precision long double
#endif

#if !defined(fortran_real)
#define fortran_real float
#endif

    extern fortran_double_precision deps(const fortran_double_precision * x_);
    extern fortran_double_precision derbit(const fortran_double_precision * relerr_,
					   const fortran_double_precision * ulp_);
    extern fortran_double_precision dfloat(const fortran_integer * n_);
    extern fortran_double_precision dfpmax(void);
    extern fortran_double_precision dgamma(const fortran_double_precision * x_);
    extern fortran_double_precision dgamma(const fortran_double_precision * x_);
    extern fortran_double_precision dinf(void);
    extern fortran_double_precision diran(const fortran_double_precision * x_,
					  const fortran_double_precision * y_);
    extern fortran_double_precision dlgam(const fortran_double_precision * x_);
    extern fortran_double_precision dlog2(const fortran_double_precision * x_);
    extern fortran_double_precision dnan(void);
    extern fortran_double_precision dpsi(const fortran_double_precision * x_);
    extern fortran_double_precision dpsiln(const fortran_double_precision * x_);
    extern fortran_double_precision dran(void);
    extern fortran_double_precision dstorf(const fortran_double_precision * x_);
    extern fortran_double_precision pythag(const fortran_double_precision * a_,
					   const fortran_double_precision * b_);

    extern fortran_integer iceil(const fortran_real * x_);
    extern fortran_integer idceil(const fortran_double_precision * x_);
    extern fortran_integer iqceil(const fortran_quadruple_precision * x_);
    extern fortran_integer trapit(const fortran_integer * sig_,
				  const fortran_integer * code_,
				  const fortran_integer * sigcon_);

    extern fortran_logical ffalse(void);
    extern fortran_logical ftrue(void);
    extern fortran_logical isainf(const fortran_real * x_);
    extern fortran_logical isanan(const fortran_real * x_);
    extern fortran_logical isdinf(const fortran_double_precision * x_);
    extern fortran_logical isdnan(const fortran_double_precision * x_);
    extern fortran_logical isqinf(const fortran_quadruple_precision * x_);
    extern fortran_logical isqnan(const fortran_quadruple_precision * x_);

    extern fortran_quadruple_precision qeps(const fortran_quadruple_precision * x_);
    extern fortran_quadruple_precision qepsln(const fortran_quadruple_precision *
					      x_);
    extern fortran_quadruple_precision qfloat(const fortran_integer * n_);
    extern fortran_quadruple_precision qfpmax(void);
    extern fortran_quadruple_precision qgamma(const fortran_quadruple_precision *
					      x_);
    extern fortran_quadruple_precision qinf(void);
    extern fortran_quadruple_precision qiran(const fortran_quadruple_precision * x_,
					     const fortran_quadruple_precision * y_);
    extern fortran_quadruple_precision qlgam(const fortran_quadruple_precision * x_);
    extern fortran_quadruple_precision qlog2(const fortran_quadruple_precision * x_);
    extern fortran_quadruple_precision qnan(void);
    extern fortran_quadruple_precision qpsi(const fortran_quadruple_precision * x_);
    extern fortran_quadruple_precision qpsiln(const fortran_quadruple_precision * x_);
    extern fortran_quadruple_precision qstorf(const fortran_quadruple_precision *
					      x_);

    extern fortran_real aeps(const fortran_real * x_);
    extern fortran_real afpmax(void);
    extern fortran_real ainf(void);
    extern fortran_real airan(const fortran_real * x_, const fortran_real * y_);
    extern fortran_real algam(const fortran_real * x_);
    extern fortran_real alog2(const fortran_real * x_);
    extern fortran_real anan(void);
    extern fortran_real astorf(const fortran_real * x_);
    extern fortran_real gamma(const fortran_real * x_);
    extern fortran_real psi(const fortran_real * x_);
    extern fortran_real psiln(const fortran_real * x_);
    extern fortran_real ran(void);

    extern void astore(const fortran_real * x_);
    extern void dpsum(fortran_double_precision * sumneg_,
		      fortran_double_precision * sumpos_,
		      const fortran_double_precision * term_);
    extern void dstore(const fortran_double_precision * x_);
    extern void qstore(const fortran_quadruple_precision * x_);

#if defined(__cplusplus)
}
#endif

#endif				/* !defined(GAMPSI_H) */
