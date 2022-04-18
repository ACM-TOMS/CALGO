#if !defined(GJLLOG_H)		/* file ignored if included more than once */
#define GJLLOG_H

/***********************************************************************

This header file for programs written in ANSI/ISO Standard C (1989) or
Standard C++ (1998) defines an interface from those languages to the
Fortran routines for Gauss-Jacobi and Gauss-Laguerre quadrature
routines with logarithmic weights, as described in the articles

    James S. Ball and Nelson H. F. Beebe,
    ``Efficient Hig-Accuracy Quadrature for Two Classes of Logarithmic Weight
    Functions'',
    ACM Transactions on Mathematical Software xx(??), ??--?? (20xx).

    James S. Ball and Nelson H. F. Beebe,
    ``Algorithm xxx: GAUSSLOG --- A Package of Routines for
    Generating Gaussian Quadrature for Two Classes of Logarithmic
    Weight Functions and the Associated Orthogonal Polynomials'',
    ACM Transactions on Mathematical Software 26(??), ??--?? (2000).

On UNIX and POSIX-compliant systems, user code should be compiled with
-I/usr/local/include, and linked with -L/usr/local/lib -lgjl,
assuming the GNU standard installation locations.  In most systems,
additional -l switches will be needed to supply the Fortran library
names.

[04-Nov-2003] -- minor updates
[06-Apr-2000]
***********************************************************************/

#if defined(__cplusplus)
extern "C"
{
#endif

/* Handle mapping of Fortran names to their equivalents in C */
#if defined(ardent) || defined(fortran_names_uppercase)
    /* Stardent (now defunct) uppercased Fortran names */
#define acopy	ACOPY
#define aeps	AEPS
#define agamma	AGAMMA
#define apsi	APSI
#define avsum	AVSUM

#define dcopy	DCOPY
#define deps	DEPS
#define dgamma	DGAMMA
#define dpsi	DPSI
#define dvsum	DVSUM
#define getnan	GETNAN

#define qcopy	QCOPY
#define qeps	QEPS
#define qgamma	QGAMMA
#define qpsi	QPSI
#define qvsum	QVSUM

#define gjqf	GJQF
#define gjqfd	GJQFD
#define gjqrc	GJQRC
#define glqf	GLQF
#define glqfd	GLQFD
#define glqrc	GLQRC
#define tql1	TQL1

#define qgjqf	QGJQF
#define qgjqfd	QGJQFD
#define qgjqrc	QGJQRC
#define qglqf	QGLQF
#define qglqfd	QGLQFD
#define qglqrc	QGLQRC
#define qtql1	QTQL1
#elif defined(_AIX) || defined(__hpux) || defined(fortran_names_same)
    /* IBM RS/6000 AIX and HP HP-UX use identical names in C and Fortran */
#else
    /* Everyone else adds a trailing underscore to Fortran names */
#define acopy	acopy_
#define aeps	aeps_
#define agamma	agamma_
#define apsi	apsi_
#define avsum	avsum_

#define dcopy	dcopy_
#define deps	deps_
#define dgamma	dgamma_
#define dpsi	dpsi_
#define dvsum	dvsum_
#define getnan	getnan_

#define qcopy	qcopy_
#define qeps	qeps_
#define qgamma	qgamma_
#define qpsi	qpsi_
#define qvsum	qvsum_

#define gjqf	gjqf_
#define gjqfd	gjqfd_
#define gjqrc	gjqrc_
#define glqf	glqf_
#define glqfd	glqfd_
#define glqrc	glqrc_
#define tql1	tql1_

#define qgjqf	qgjqf_
#define qgjqfd	qgjqfd_
#define qgjqrc	qgjqrc_
#define qglqf	qglqf_
#define qglqfd	qglqfd_
#define qglqrc	qglqrc_
#define qtql1	qtql1_
#endif

/* Handle mappings of Fortran datatypes to C/C++ equivalents.  We
   define these with macros, rather than typedefs, in order to be able
   to test for user-defined overriding values. */

#if !defined(fortran_single_precision)
#define fortran_single_precision float
#endif

#if !defined(fortran_double_precision)
#define fortran_double_precision double
#endif

#if !defined(fortran_quadruple_precision)
#define fortran_quadruple_precision long double
#endif

#if !defined(fortran_integer)
#define fortran_integer int
#endif

    void acopy(const fortran_integer * n_,
	       fortran_single_precision dx_[],
	       const fortran_integer * incx_,
	       fortran_single_precision dy_[],
	       const fortran_integer * incy_);

    fortran_single_precision aeps(const fortran_single_precision * x_);

    fortran_single_precision agamma(const fortran_single_precision * x_);

    fortran_single_precision apsi(const fortran_single_precision * x_);

    fortran_single_precision avsum(fortran_single_precision * x_,
	       const fortran_integer * n_);

    void dcopy(const fortran_integer * n_,
	       fortran_double_precision dx_[],
	       const fortran_integer * incx_,
	       fortran_double_precision dy_[],
	       const fortran_integer * incy_);

    fortran_double_precision deps(const fortran_double_precision * x_);

    fortran_double_precision dgamma(const fortran_double_precision * x_);

    fortran_double_precision dpsi(const fortran_double_precision * x_);

    fortran_double_precision dvsum(fortran_double_precision * x_,
				   const fortran_integer * n_);

    fortran_double_precision getnan(void);

    void gjqf(fortran_double_precision x_[],
	      fortran_double_precision w_[],
	      fortran_double_precision y_[],
	      fortran_double_precision z_[],
	      const fortran_double_precision * alpha_,
	      const fortran_double_precision * beta_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void gjqfd(fortran_double_precision x_[],
	       fortran_double_precision w_[],
	       fortran_double_precision deltaw_[],
	       fortran_double_precision deltax_[],
	       const fortran_double_precision * alpha_,
	       const fortran_double_precision * beta_,
	       const fortran_integer * nquad_,
	       fortran_integer * ierr_);

    void gjqrc(fortran_double_precision a_[],
	      fortran_double_precision b_[],
	      fortran_double_precision s_[],
	      fortran_double_precision t_[],
	      const fortran_double_precision * alpha_,
	      const fortran_double_precision * beta_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void glqf(fortran_double_precision x_[],
	      fortran_double_precision w_[],
	      fortran_double_precision wxm1_[],
	      fortran_double_precision y_[],
	      fortran_double_precision z_[],
	      const fortran_double_precision * alpha_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void glqfd(fortran_double_precision x_[],
	       fortran_double_precision w_[],
	       fortran_double_precision deltaw_[],
	       fortran_double_precision deltax_[],
	       const fortran_double_precision * alpha_,
	       const fortran_integer * nquad_,
	       fortran_integer * ierr_);

    void glqrc(fortran_double_precision a_[],
	      fortran_double_precision b_[],
	      fortran_double_precision s_[],
	      fortran_double_precision t_[],
	      const fortran_double_precision * alpha_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void tql1(const fortran_integer * n_,
	      fortran_double_precision d_[],
	      fortran_double_precision e_[],
	      fortran_integer * ierr_);

#if defined(HAVE_LONG_DOUBLE)

    void qcopy(const fortran_integer * n_,
	       fortran_quadruple_precision dx_[],
	       const fortran_integer * incx_,
	       fortran_quadruple_precision dy_[],
	       const fortran_integer * incy_);

    fortran_quadruple_precision qeps(const fortran_quadruple_precision * x_);

    fortran_quadruple_precision qgamma(const fortran_quadruple_precision * x_);

    fortran_quadruple_precision qpsi(const fortran_quadruple_precision * x_);

    fortran_quadruple_precision qvsum(fortran_quadruple_precision * x_,
	       const fortran_integer * n_);

    void qgjqf(fortran_quadruple_precision x_[],
	      fortran_quadruple_precision w_[],
	      fortran_quadruple_precision y_[],
	      fortran_quadruple_precision z_[],
	      const fortran_quadruple_precision * alpha_,
	      const fortran_quadruple_precision * beta_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void qgjqfd(fortran_quadruple_precision x_[],
	       fortran_quadruple_precision w_[],
	       fortran_quadruple_precision deltaw_[],
	       fortran_quadruple_precision deltax_[],
	       const fortran_quadruple_precision * alpha_,
	       const fortran_quadruple_precision * beta_,
	       const fortran_integer * nquad_,
	       fortran_integer * ierr_);

    void qgjqrc(fortran_quadruple_precision a_[],
	      fortran_quadruple_precision b_[],
	      fortran_quadruple_precision s_[],
	      fortran_quadruple_precision t_[],
	      const fortran_quadruple_precision * alpha_,
	      const fortran_quadruple_precision * beta_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void qglqf(fortran_quadruple_precision x_[],
	      fortran_quadruple_precision w_[],
	      fortran_quadruple_precision wxm1_[],
	      fortran_quadruple_precision y_[],
	      fortran_quadruple_precision z_[],
	      const fortran_quadruple_precision * alpha_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void qglqfd(fortran_quadruple_precision x_[],
	       fortran_quadruple_precision w_[],
	       fortran_quadruple_precision deltaw_[],
	       fortran_quadruple_precision deltax_[],
	       const fortran_quadruple_precision * alpha_,
	       const fortran_integer * nquad_,
	       fortran_integer * ierr_);

    void qglqrc(fortran_quadruple_precision a_[],
	      fortran_quadruple_precision b_[],
	      fortran_quadruple_precision s_[],
	      fortran_quadruple_precision t_[],
	      const fortran_quadruple_precision * alpha_,
	      const fortran_integer * nquad_,
	      fortran_integer * ierr_);

    void qtql1(const fortran_integer * n_,
	      fortran_quadruple_precision d_[],
	      fortran_quadruple_precision e_[],
	      fortran_integer * ierr_);
#endif /* defined(HAVE_LONG_DOUBLE) */

#if defined(__cplusplus)
}
#endif

#endif				/* !defined(GJLLOG_H) */
