/**********************    OMP_Talbot_pack_DE.c    ************************
 *                                                                        *
 *                                                                        *
 *                  parallel version of TALBOT SUITE DE                   *
 *                                                                        *
 *                APPLICATION OF MODIFIED TALBOT'S METHODS                *
 *                     TO SOLVE DIFFERENTIAL PROBLEMS                     *
 *                                                                        *
 *           THE LAPLACE TRANSFORM IS GIVEN BY NUMERICAL SAMPLES          *
 *                             AND NOT AS A FUNCTION                      *
 *                                                                        *
 *                                                                        *
 *           OpenMP-based parallel version of TALBOT SUITE DE             *
 *                                                                        *
 *                            path: code/SRC/FUN_DE/                      *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>       VERSION 4.0     May 25th, 2016       <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 *                                                                        *
 *                       AUTHOR: Mariarosaria Rizzardi                    *
 *                                                                        *
 *                  mariarosaria.rizzardi@uniparthenope.it                *
 *                                                                        *
 *                  DiST - Dept. of Science and Technology                *
 *                  "Parthenope" University, Naples (Italy)               *
 *                                                                        *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 *  Algorithm's steps are:                                                *
 *                                                                        *
 *      1) compute Talbot's parameters at (Tmin + Tmax)/2                 *
 *         by means of COM_TalbotPAR function;                            *
 *                                                                        *
 *      2) compute Laplace Transform samples U(x,s) for s on              *
 *         Talbot's contour;                                              *
 *                                                                        *
 *      3) for all x,t do                                                 *
 *              approximate the Inverse Laplace Transform u(x,t)          *
 *              by means of OMP_TalbotSUM11_DE, OMP_TalbotSUM12_DE or     *
 *              OMP_TalbotSUM13_DE functions.                             *
 *                                                                        *
 *              COM_TalbotPAR and OMP_TalbotSUM11_DE, OMP_TalbotSUM12_DE  *
 *              or OMP_TalbotSUM13_DE are skill-level functions.          *
 *              The first is from Talbot Suite (file COM_Talbot_pack.c).  *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "A modification of Talbot's method for the simultaneous   *
 *               approximation of several values of the Inverse Laplace   *
 *               Transform".                                              *
 *               ACM Trans. Math. Soft., vol. 21, no. 4, Dec. 1995,       *
 *               pp. 347-371.                                             *
 *                                                                        *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *               ACM TRANS. MATH. SOFTWARE, VOL. xx, NO. x, month year,   *
 *               pp ##.                                                   *
 *                                                                        *
 **************************************************************************/

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "../COM/main_include.h"
#include "../COM/COM_Talbot_pack.h"
#include "../COM_DE/COM_Talbot_pack_DE.h"
#include "OMP_Talbot_pack_DE.h"


int OMP_Talbot11_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS)
/**************************************************************************

    OMP_Talbot11_DE   DRIVER FUNCTION (user level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version

    PURPOSE
    =======
    This function provides a numerical approximation to the Inverse Laplace
    Transform u(x,t) computed at each value of the  Xval,Tval  arrays.

    A coarse-grained parallelism is implemented; parallelization strategy is
    data distribution.


    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_Talbot11_DE (LTsamples, sigma0, NXval, Xval, NTval, Tval,
                                     tol, NUMft, IFAIL, Nsings, SINGS, MULT,
                                     Tmin, Tmax, THREADS);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTsamples - (double complex function pointer) user-defined function according
                to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

                The function returns the LT samples, computed in parallel, by
                solving the problem obtained by the application of the Laplace
                transform method to the original differential problem.

    sigma0    - (double) abscissa of (absolute) convergence of the Laplace Transform
                function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    Xval      - (double array) values for x where the Inverse Laplace Transform
                u(x,t) is approximated.
                It must be dimensioned at least "NXval".

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. Its components must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |u(x,t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays SINGS and MULT.

    SINGS     - (double complex array) singularities of the Laplace
                Transform function. Only singularities with non-negative
                imaginary parts are required; their complex conjugates are
                unnecessary.
                It must be dimensioned at least "Nsings".

    MULT      - (unsigned integer array) multiplicities of those singularities
                (in SINGS) which are poles, zero otherwise.
                It must be dimensioned as SINGS.

    Tmin,Tmax - (double) endpoints of the interval enclosing the t values.
                Method's parameters are computed at (Tmin + Tmax)/2.

    THREADS   - (integer) number of parallel OpenMP threads to be used in
                parallel regions.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                      u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t(k))


    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : compute the Talbot parameters (in COM_Talbot_pack.c).

    OMP_TalbotSUM11_DE: approximate the Inverse Laplace Transform.

    LTsamples     : (user-defined function) Laplace Transform function
                    according to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU;
    unsigned int NOPTS;
    int IFAIL_tot;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;


    /* 1) Compute Talbot's parameters at the middle point of the interval [Tmin, Tmax] */
    COM_TalbotPAR (sigma0,(Tmin + Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);

    /* 2) Compute the correction to NOPTS */
    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);


    /* 3) Compute, in parallel, sample points on Talbot's contour */
    double complex *S = (double complex*)malloc(NOPTS*sizeof(double complex));
    if ( S == NULL )
    {   fprintf(stderr, "\n***   ERROR in OMP_Talbot1_DE: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");
        exit(1);
    }

    double thetak,    piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
    unsigned int k;
    #pragma omp parallel for    default   (shared)   \
                                private   (thetak,k) \
                                num_threads (THREADS)
    for (k=1; k<NOPTS; k++)
    {   thetak = k*piN;
        S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
    }

    S[0] = CONSIG + CONLAM;


    /* 4) Compute the matrix FS of LT samples on Talbot's contour **/
    double complex *FS = (*LTsamples)(NXval,Xval,NOPTS,S,tol,THREADS);
    free(S);


    /* 5) Compute, in parallel and with the same LT samples:
          for all x and for all t
                    u(x,t),
                    local error flags
                    and a global error flag
     */
    IFAIL_tot = OMP_TalbotSUM11_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,Tval,NUMft,IFAIL,THREADS);
    free(FS);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int OMP_Talbot12_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS)
/**************************************************************************

    OMP_Talbot12_DE   DRIVER FUNCTION (user level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version

    PURPOSE
    =======
    This function provides a numerical approximation to the Inverse Laplace
    Transform u(x,t) computed at each value of the  Xval,Tval  arrays.

    A fine-grained parallelism is implemented; parallelization strategy is
    task distribution, i.e. the summation process has been parallelized.


    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_Talbot12_DE (LTsamples, sigma0, NXval, Xval, NTval, Tval,
                                     tol, NUMft, IFAIL, Nsings, SINGS, MULT,
                                     Tmin, Tmax, THREADS);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTsamples - (double complex function pointer) user-defined function according
                to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

                The function returns the LT samples, computed in parallel, by
                solving the problem obtained by the application of the Laplace
                transform method to the original differential problem.

    sigma0    - (double) abscissa of (absolute) convergence of the Laplace Transform
                function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    Xval      - (double array) values for x where the Inverse Laplace Transform
                u(x,t) is approximated.
                It must be dimensioned at least "NXval".

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. Its components must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |u(x,t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays SINGS and MULT.

    SINGS     - (double complex array) singularities of the Laplace
                Transform function. Only singularities with non-negative
                imaginary parts are required; their complex conjugates are
                unnecessary.
                It must be dimensioned at least "Nsings".

    MULT      - (unsigned integer array) multiplicities of those singularities
                (in SINGS) which are poles, zero otherwise.
                It must be dimensioned as SINGS.

    Tmin,Tmax - (double) endpoints of the interval enclosing the t values.
                Method's parameters are computed at (Tmin + Tmax)/2.

    THREADS   - (integer) number of parallel OpenMP threads to be used in
                parallel regions.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                      u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t(k))


    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : compute the Talbot parameters (in COM_Talbot_pack.c).

    OMP_TalbotSUM12_DE: approximate the Inverse Laplace Transform.

    LTsamples     : (user-defined function) Laplace Transform function
                    according to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU;
    unsigned int NOPTS;
    int IFAIL_tot;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;


    /* 1) Compute Talbot's parameters at the middle point of the interval [Tmin, Tmax] */
    COM_TalbotPAR (sigma0,(Tmin + Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);

    /* 2) Compute the correction to NOPTS */
    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);


    /* 3) Compute, in parallel, sample points on Talbot's contour */
    double complex *S = (double complex*)malloc(NOPTS*sizeof(double complex));
    if ( S == NULL )
    {   fprintf(stderr, "\n***   ERROR in OMP_Talbot1_DE: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");
        exit(1);
    }

    double thetak,    piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
    unsigned int k;
    #pragma omp parallel for    default   (shared)   \
                                private   (thetak,k) \
                                num_threads (THREADS)
    for (k=1; k<NOPTS; k++)
    {   thetak = k*piN;
        S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
    }

    S[0] = CONSIG + CONLAM;


    /* 4) Compute the matrix FS of LT samples on Talbot's contour **/
    double complex *FS = (*LTsamples)(NXval,Xval,NOPTS,S,tol,THREADS);
    free(S);


    /* 5) Compute, in parallel and with the same LT samples:
          for all x and for all t
                    u(x,t),
                    local error flags
                    and a global error flag
     */
    IFAIL_tot = OMP_TalbotSUM12_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,Tval,NUMft,IFAIL,THREADS);
    free(FS);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int OMP_Talbot13_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS1, int THREADS2)
/**************************************************************************

    OMP_Talbot13_DE   DRIVER FUNCTION (user level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version

    PURPOSE
    =======
    This function provides a numerical approximation to the Inverse Laplace
    Transform u(x,t) computed at each value of the  Xval,Tval  arrays.

    A hybrid parallelism is implemented by means of OpenMP nested parallelism,
    that must be enabled. Outer parallelization strategy is data distribution,
    inner parallelization strategy is task distribution.


    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_Talbot13_DE (LTsamples, sigma0, NXval, Xval, NTval, Tval,
                                     tol, NUMft, IFAIL, Nsings, SINGS, MULT,
                                     Tmin, Tmax, THREADS1, THREADS2);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTsamples - (double complex function pointer) user-defined function according
                to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

                The function returns the LT samples, computed in parallel, by
                solving the problem obtained by the application of the Laplace
                transform method to the original differential problem.

    sigma0    - (double) abscissa of (absolute) convergence of the Laplace Transform
                function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    Xval      - (double array) values for x where the Inverse Laplace Transform
                u(x,t) is approximated.
                It must be dimensioned at least "NXval".

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. Its components must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |u(x,t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays SINGS and MULT.

    SINGS     - (double complex array) singularities of the Laplace
                Transform function. Only singularities with non-negative
                imaginary parts are required; their complex conjugates are
                unnecessary.
                It must be dimensioned at least "Nsings".

    MULT      - (unsigned integer array) multiplicities of those singularities
                (in SINGS) which are poles, zero otherwise.
                It must be dimensioned as SINGS.

    Tmin,Tmax - (double) endpoints of the interval enclosing the t values.
                Method's parameters are computed at (Tmin + Tmax)/2.

    THREADS1,THREADS2 - (integer) number of parallel OpenMP threads to be used in
                nested parallel regions of the summation step, respectively for
                outer and inner region.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                         u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t)


    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : compute the Talbot parameters (in COM_Talbot_pack.c).

    OMP_TalbotSUM13_DE: approximate the Inverse Laplace Transform.

    LTsamples     : (user-defined function) Laplace Transform function
                    according to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol, int THREADS)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU;
    unsigned int NOPTS;
    int IFAIL_tot;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;


    /* 1) Compute Talbot's parameters at the middle point of the interval [Tmin, Tmax] */
    COM_TalbotPAR (sigma0,(Tmin + Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);

    /* 2) Compute the correction to NOPTS */
    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);


    /* 3) Compute, in parallel, sample points on Talbot's contour */
    double complex *S = (double complex*)malloc(NOPTS*sizeof(double complex));
    if ( S == NULL )
    {   fprintf(stderr, "\n***   ERROR in OMP_Talbot1_DE: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");
        exit(1);
    }

    double thetak,    piN = 4.0*atan(1.0)/NOPTS; /* pi over NOPTS (step) */
    unsigned int k;
    #pragma omp parallel for    default   (shared)   \
                                private   (thetak,k) \
                                num_threads (THREADS1*THREADS2)
    for (k=1; k<NOPTS; k++)
    {   thetak = k*piN;
        S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
    }

    S[0] = CONSIG + CONLAM;


    /* 4) Compute the matrix FS of LT samples on Talbot's contour **/
    double complex *FS = (*LTsamples)(NXval,Xval,NOPTS,S,tol,THREADS1*THREADS2);
    free(S);


    /* 5) Compute, in parallel and with the same LT samples:
          for all x and for all t
                    u(x,t),
                    local error flags
                    and a global error flag
     */
    IFAIL_tot = OMP_TalbotSUM13_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,Tval,NUMft,IFAIL,THREADS1,THREADS2);
    free(FS);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int OMP_TalbotSUM11_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS)
/**************************************************************************

    OMP_TalbotSUM11_DE   SUMMATION FUNCTION (skill level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version - coarse-grained parallelism

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace Transform
    u(x,t) evaluated at each value of the  Xval,Tval  arrays. This is accomplished
    according to the modified Talbot method.

    A coarse-grained parallelism is implemented; parallelization strategy is
    data distribution.

    The composite Trapezoidal rule, approximating the contour integral for u(x,t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_TalbotSUM11_DE (CONLAM, CONSIG, CONNU, NOPTS, NXval, FF,
                                        NTval, Tval, NUMft, IFAIL, THREADS);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of terms in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    FF        - (double complex array) row-wise matrix, of size (NXval,NOPTS),
                containing the Laplace Transform samples on the Talbot contour
                to be used in summation to invert the Laplace Transform.
                The row-wise matrix FF is stored in a mono-dimensional array F
                as
                        FF(jX,jS) = F(j) = F( jX*NOPTS + jS )
                where
                        jX is the integer quotient  j/NOPTS
                        jS is the integer remainder j%NOPTS

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. It must be dimensioned at least "NTval".

    THREADS   - (integer) number of parallel OpenMP threads to be used in
                parallel regions.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                         u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t)

    NUMft and IFAIL are row-wise matrices MM, of size (NXval,NTval), stored in
    a mono-dimensional array M as:
                M(jX,jT) = M(j) = M( jX*NTval + jT )
        so that
                jX is the integer quotient  j/NTval
                jT is the integer remainder j%Ntval


    REQUIRED FUNCTIONS
    ==================
    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.
    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    double pi = 4.*atan(1.0);
    double PIOVN = pi/(double)NOPTS;
    double TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN, ft_loc;
    int          K,   NM1 = NOPTS - 1;
    int          IFAIL_tot = 0;
    int          mod, myid, IFAIL_tot_loc, IFAIL_loc;
    unsigned int j, jX, jT,   NTOT = NXval*NTval;
    unsigned int NTOTloc, STARTloc, ENDloc;


    /* PARALLEL SECTION */
    #pragma omp parallel default   (shared)                                                     \
                         private   (TAU,PSI,C,U,BR,BI,DBR,DBI,K,TETA,ALFA,BETA,G,H,EAT,SIGN,OF) \
                         private   (NTOTloc,j,jT,jX)                                            \
                         private   (mod,myid,STARTloc,ENDloc,IFAIL_tot_loc,ft_loc,IFAIL_loc)    \
                         reduction (|| : IFAIL_tot)                                             \
                         num_threads (THREADS)
    {
        /* index local distribution in output data */
        NTOTloc = NTOT/THREADS;    mod = NTOT%THREADS;
        #ifdef _OPENMP
            myid = omp_get_thread_num();
        #else
            myid = 0;
        #endif
        if (myid < mod)
        {
            NTOTloc  = NTOTloc+1;
            STARTloc = myid*NTOTloc;
        }
        else
            STARTloc = myid*NTOTloc + mod;
        ENDloc = STARTloc + NTOTloc - 1;

        IFAIL_tot_loc = 0;               /* local global error flag */
        for (j=STARTloc; j<=ENDloc; j++) /* loop on local values */
        {
            ft_loc = 0.0;    IFAIL_loc = 0; /* initialize output local values */

            jX = j/NTval; /* integer quotient:  jX is the index on x-values */
            jT = j%NTval; /* integer remainder: jT is the index on t-values */

            TAU = CONLAM*Tval[jT];

            if (NM1 > 0)
            {
                PSI=PIOVN*TAU*CONNU; C=cos(PSI);
                BR=0.; BI=0.; DBR=0.; DBI=0.;

                /* begin  Goertzel-Reinsch algorithm   for   Clenshaw sums  */
                if (C <= 0.)
                {
                    U = +4.*pow(cos(PSI/2),2);
                    for (K=NM1; K>=1; K--)
                    {
                        TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                                      BETA=TETA+ALFA*(ALFA-1.)/TETA;
                        G=creal(FF[jX*NOPTS+K]); H=cimag(FF[jX*NOPTS+K]); /* FF is a row-wise matrix of size (NXval,NOPTS) */
                        EAT = exp(ALFA*TAU);
                        BR = DBR - BR;
                        BI = DBI - BI;
                        DBR = U*BR - DBR + EAT*(G*CONNU - H*BETA);
                        DBI = U*BI - DBI + EAT*(H*CONNU + G*BETA);
                    }
                    BR = DBR - BR;
                    BI = DBI - BI;
                    DBR = U*BR - DBR;
                }
                else /* C=cos(PSI) > 0 */
                {
                    U = -4.*pow(sin(PSI/2),2);
                    for (K=NM1; K>=1; K--)
                    {
                        TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                                      BETA=TETA+ALFA*(ALFA-1.)/TETA;
                        G=creal(FF[jX*NOPTS+K]); H=cimag(FF[jX*NOPTS+K]); /* FF is a row-wise matrix of size (NXval,NOPTS) */
                        EAT = exp(ALFA*TAU);
                        BR = DBR + BR;
                        BI = DBI + BI;
                        DBR = U*BR + DBR + EAT*(G*CONNU - H*BETA);
                        DBI = U*BI + DBI + EAT*(H*CONNU + G*BETA);
                    }
                    BR = DBR + BR;
                    BI = DBI + BI;
                    DBR = U*BR + DBR;
                }
                /* end   Goertzel-Reinsch algorithm   for   Clenshaw sums  */

                ft_loc = DBR - BR*U/2 - BI*sin(PSI);
            } /* end-if (NM1 > 0) */

            /* compute   T0 (the first term in the summation) */
            ft_loc = ft_loc + CONNU*exp(TAU)*creal(FF[jX*NOPTS])/2.0;
            if ( ft_loc != 0.0 )
            {
                if ( ft_loc > 0.0 )
                    SIGN = +1.0;
                else
                    SIGN = -1.0;

                /* overflow   problem   handling */
                OF = CONSIG*Tval[jT] + log(CONLAM*fabs( ft_loc )/NOPTS);
                if (OF > log(DBL_MAX))
                    IFAIL_loc = 1;
                else
                    ft_loc = exp(OF)*SIGN;
            }

            NUMft[j] = ft_loc;   IFAIL[j] = IFAIL_loc;
            IFAIL_tot_loc = IFAIL_tot_loc || IFAIL_loc;

        } /* end-loop on local values for (j=STARTloc; j<=ENDloc; j++) */

        IFAIL_tot = IFAIL_tot || IFAIL_tot_loc; /* reduction for the global "or" */

    } /* end omp parallel section */

    return IFAIL_tot;
/**************************************************************************/
}


int OMP_TalbotSUM12_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS)
/**************************************************************************

    OMP_TalbotSUM12_DE   SUMMATION FUNCTION (skill level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version - fine-grained parallelism

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace Transform
    u(x,t) evaluated at each value of the  Xval,Tval  arrays. This is accomplished
    according to the modified Talbot method.

    A fine-grained parallelism is implemented; parallelization strategy is
    task distribution, i.e. the summation process has been parallelized.

    The composite Trapezoidal rule, approximating the contour integral for u(x,t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_TalbotSUM12_DE (CONLAM, CONSIG, CONNU, NOPTS, NXval, FF,
                                        NTval, Tval, NUMft, IFAIL, THREADS);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of terms in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    FF        - (double complex array) row-wise matrix, of size (NXval,NOPTS),
                containing the Laplace Transform samples on the Talbot contour
                to be used in summation to invert the Laplace Transform.
                The row-wise matrix FF is stored in a mono-dimensional array F
                as
                        FF(jX,jS) = F(j) = F( jX*NOPTS + jS )
                where
                        jX is the integer quotient  j/NOPTS
                        jS is the integer remainder j%NOPTS

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. It must be dimensioned at least "NTval".

    THREADS   - (integer) number of parallel OpenMP threads to be used in
                parallel regions.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                         u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t)

    NUMft and IFAIL are row-wise matrices MM, of size (NXval,NTval), stored in
    a mono-dimensional array M as:
                M(jX,jT) = M(j) = M( jX*NTval + jT )
        so that
                jX is the integer quotient  j/NTval
                jT is the integer remainder j%Ntval


    REQUIRED FUNCTIONS
    ==================
    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.
    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    double pi = 4.*atan(1.0);
    double PIOVN = pi/(double)NOPTS;
    double TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN, sum;

    unsigned int j, jX, jT,   NTOT = NXval*NTval;
    int IFAIL_tot=0;


    for (j=0; j<NTOT; j++) /* loop on u(x,t) */
        {
            jX = j/NTval; /* integer quotient:  jX is the index on x-values */
            jT = j%NTval; /* integer remainder: jT is the index on t-values */

            IFAIL[j] = 0;
            NUMft[j] = 0.0;
            sum = 0.0;
            TAU=CONLAM*Tval[jT];

            if (NOPTS > 1)
            {
                int K, STARTloc, ENDloc;
                unsigned int Nloc, mod, myid; double sum_loc;
                PSI=PIOVN*TAU*CONNU;   C=cos(PSI);

                /* PARALLEL SECTION */
                #pragma omp parallel default   (shared)                                   \
                                     private   (BR,BI,DBR,DBI,U,K,TETA,ALFA,BETA,G,H,EAT) \
                                     private   (Nloc,mod,myid,STARTloc,ENDloc,sum_loc)    \
                                     reduction (+ : sum)                                  \
                                     num_threads (THREADS)
                {
                    /* local summation index distribution (jS) */
                    Nloc=NOPTS/THREADS; mod=NOPTS%THREADS;

                    #ifdef _OPENMP
                        myid = omp_get_thread_num();
                    #else
                        myid = 0;
                    #endif

                    if (myid < mod)
                    {
                        Nloc = Nloc+1;
                        STARTloc = myid*Nloc;
                    }
                    else
                        STARTloc = myid*Nloc + mod;

                    ENDloc = STARTloc + Nloc - 1;

                    /* T0 is computed outside the parallel section */
                    if (myid == 0)
                        STARTloc = STARTloc + 1;

                    /* begin   Goertzel-Reinsch algorithm   for   Clenshaw sums  */
                    BR=0.; BI=0.; DBR=0.; DBI=0.;
                    if (C <= 0.) /* C=cos(PSI) <= 0 */
                    {
                        U = +4.*pow(cos(PSI/2),2);
                        for (K=ENDloc; K>=STARTloc; K--)
                        {
                            TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                                          BETA=TETA+ALFA*(ALFA-1.)/TETA;
                            /* FF(jX,jS) = FF[jX*NOPTS + jS] */
                            G=creal(FF[jX*NOPTS + K]); H=cimag(FF[jX*NOPTS + K]);
                            EAT = exp(ALFA*TAU);
                            BR = DBR - BR;
                            BI = DBI - BI;
                            DBR = U*BR - DBR + EAT*(G*CONNU - H*BETA);
                            DBI = U*BI - DBI + EAT*(H*CONNU + G*BETA);
                        }
                    }
                    else         /* C=cos(PSI) > 0 */
                    {
                        U = -4.*pow(sin(PSI/2),2);
                        for (K=ENDloc; K>=STARTloc; K--)
                        {
                            TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                            BETA=TETA+ALFA*(ALFA-1.)/TETA;
                            /* FF(jX,jS) = FF[jX*NOPTS + jS] */
                            G=creal(FF[jX*NOPTS + K]); H=cimag(FF[jX*NOPTS + K]);
                            EAT = exp(ALFA*TAU);
                            BR = DBR + BR;
                            BI = DBI + BI;
                            DBR = U*BR + DBR + EAT*(G*CONNU - H*BETA);
                            DBI = U*BI + DBI + EAT*(H*CONNU + G*BETA);
                        }
                    }

                    DBR = DBR - BR*U/2;
                    DBI = DBI - BI*U/2;
                    BR = BR*sin(PSI);
                    BI = BI*sin(PSI);

                    sum_loc = DBR*cos(STARTloc*PSI) - BR*sin(STARTloc*PSI) - BI*cos(STARTloc*PSI) - DBI*sin(STARTloc*PSI);
                    /* end   Goertzel-Reinsch algorithm   for   Clenshaw sums  */

                    sum += sum_loc; /* reduction for the global sum */

                } /* end parallel section */

            } /* end if (NOPTS > 1) */

            /* Master thread computes T0 (the first term in the summation) */
            NUMft[j] = sum + CONNU*exp(TAU)*creal(FF[jX*NOPTS])/2.0; /* FF(jX,0) */
            if (NUMft[j] != 0.0)
            {
                if (NUMft[j] > 0.0)
                    SIGN = +1.0;
                else
                    SIGN = -1.0;

                /* overflow   problem   handling */
                OF = CONSIG*Tval[jT] + log(CONLAM*fabs(NUMft[j])/NOPTS);
                if (OF > log(DBL_MAX))
                {
                    IFAIL[j] = 1;
                    IFAIL_tot = 1; /* logical or of IFAIL(jX,jT) */
                }
                else
                    NUMft[j] = exp(OF)*SIGN;
            }
        } /* end for-loop j = jX*NTval+jT  to compute u(jX,jT) */

    return IFAIL_tot;
/**************************************************************************/
}


int OMP_TalbotSUM13_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS1, int THREADS2)
/**************************************************************************

    OMP_TalbotSUM13_DE   SUMMATION FUNCTION (skill level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION
                          OpenMP-based version - nested parallelism

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace Transform
    u(x,t) evaluated at each value of the  Xval,Tval  arrays. This is accomplished
    according to the modified Talbot method.

    A hybrid parallelism is implemented by means of OpenMP nested parallelism,
    that must be enabled. Outer parallelization strategy is data distribution,
    inner parallelization strategy is task distribution.

    The composite Trapezoidal rule, approximating the contour integral for u(x,t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_TalbotSUM13_DE (CONLAM, CONSIG, CONNU, NOPTS, NXval, FF,
                                        NTval, Tval, NUMft, IFAIL, THREADS1, THREADS2);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of terms in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    FF        - (double complex array) row-wise matrix, of size (NXval,NOPTS),
                containing the Laplace Transform samples on the Talbot contour
                to be used in summation to invert the Laplace Transform.
                The row-wise matrix FF is stored in a mono-dimensional array F
                as
                        FF(jX,jS) = F(j) = F( jX*NOPTS + jS )
                where
                        jX is the integer quotient  j/NOPTS
                        jS is the integer remainder j%NOPTS

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. It must be dimensioned at least "NTval".

    THREADS1, THREADS2  - (integer) number of parallel OpenMP threads to be used in
                parallel regions. For nested parallelism, the former refers to outer
                parallelism and the latter to inner parallelism.


    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)).

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)):

                                / 0 no error
                  IFAIL(h,k) = |
                                \ 1 an overflow occurs in u(x(h),t(k)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                         u(x(h),t(k)) = u(x(h),t(k))/exp(sigma0*t)

    NUMft and IFAIL are row-wise matrices MM, of size (NXval,NTval), stored in
    a mono-dimensional array M as:
                M(jX,jT) = M(j) = M( jX*NTval + jT )
        so that
                jX is the integer quotient  j/NTval
                jT is the integer remainder j%Ntval


    REQUIRED FUNCTIONS
    ==================
    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.
    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    double pi = 4.*atan(1.0);
    double PIOVN = pi/(double)NOPTS;
    double TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN, sum, sum_loc;

    unsigned int j, jX, jT,   NTOT = NXval*NTval;
    unsigned int Nloc, mod, myid;
    int K, STARTloc, ENDloc,   IFAIL_tot = 0;


    /* outer loop on u(x,t) */
    #pragma omp parallel for num_threads(THREADS1)                              \
                             private (j,jX,jT,sum,TAU,PSI,C)                          \
                             private (BR,BI,DBR,DBI,U,K,TETA,ALFA,BETA,G,H,EAT) \
                             private (Nloc,mod,myid,STARTloc,ENDloc,sum_loc)
        for (j=0; j<NTOT; j++)
        {
            IFAIL[j] = 0;
            NUMft[j] = 0.0;
            sum = 0.0;
            jX = j/NTval; /* integer quotient:  jX is the index on x-values */
            jT = j%NTval; /* integer remainder: jT is the index on t-values */
            TAU=CONLAM*Tval[jT];

            if (NOPTS > 1)
            {
                PSI=PIOVN*TAU*CONNU;   C=cos(PSI);

                #pragma omp parallel num_threads (THREADS2)                               \
                                     private   (BR,BI,DBR,DBI,U,K,TETA,ALFA,BETA,G,H,EAT) \
                                     private   (Nloc,mod,myid,STARTloc,ENDloc,sum_loc)    \
                                     reduction (+ : sum)
                {
                    /* local summation index distribution (jS) */
                    Nloc=NOPTS/THREADS2; mod=NOPTS%THREADS2;

                    #ifdef _OPENMP
                        myid = omp_get_thread_num();
                    #else
                        myid = 0;
                    #endif
                    if (myid < mod)
                    {   Nloc = Nloc+1;
                        STARTloc = myid*Nloc;
                    }
                    else
                        STARTloc = myid*Nloc + mod;
                    ENDloc = STARTloc + Nloc - 1;

                    /* T0 is computed outside the parallel section */
                    if (myid == 0)
                        STARTloc = STARTloc + 1;

                    /* begin   Goertzel-Reinsch algorithm   for   Clenshaw sums  */
                    BR=0.; BI=0.; DBR=0.; DBI=0.;
                    if (C <= 0.) /* C=cos(PSI) <= 0 */
                    {
                        U = +4.*pow(cos(PSI/2),2);
                        for (K=ENDloc; K>=STARTloc; K--)
                        {
                            TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                                          BETA=TETA+ALFA*(ALFA-1.)/TETA;
                            /* FF(jX,jS) = FF[jX*NOPTS + jS] */
                            G=creal(FF[jX*NOPTS + K]); H=cimag(FF[jX*NOPTS + K]);
                            EAT = exp(ALFA*TAU);
                            BR = DBR - BR;
                            BI = DBI - BI;
                            DBR = U*BR - DBR + EAT*(G*CONNU - H*BETA);
                            DBI = U*BI - DBI + EAT*(H*CONNU + G*BETA);
                        }
                    }
                    else         /* C=cos(PSI) > 0 */
                    {
                        U = -4.*pow(sin(PSI/2),2);
                        for (K=ENDloc; K>=STARTloc; K--)
                        {
                            TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                            BETA=TETA+ALFA*(ALFA-1.)/TETA;
                            /* FF(jX,jS) = FF[jX*NOPTS + jS] */
                            G=creal(FF[jX*NOPTS + K]); H=cimag(FF[jX*NOPTS + K]);
                            EAT = exp(ALFA*TAU);
                            BR = DBR + BR;
                            BI = DBI + BI;
                            DBR = U*BR + DBR + EAT*(G*CONNU - H*BETA);
                            DBI = U*BI + DBI + EAT*(H*CONNU + G*BETA);
                        }
                    }

                    DBR = DBR - BR*U/2;
                    DBI = DBI - BI*U/2;
                    BR = BR*sin(PSI);
                    BI = BI*sin(PSI);

                    sum_loc = DBR*cos(STARTloc*PSI) - BR*sin(STARTloc*PSI) - BI*cos(STARTloc*PSI) - DBI*sin(STARTloc*PSI);
                    /* end   Goertzel-Reinsch algorithm   for   Clenshaw sums  */

                    sum += sum_loc; /* reduction for the global sum */

                } /* end parallel section */

            } /* end if (NOPTS > 1) */

            /* Master thread computes T0 (the first term in the summation) */
            NUMft[j] = sum + CONNU*exp(TAU)*creal(FF[jX*NOPTS])/2.0; /* FF(jX,0) */
            if (NUMft[j] != 0.0)
            {
                if (NUMft[j] > 0.0)
                    SIGN = +1.0;
                else
                    SIGN = -1.0;

                /* overflow   problem   handling */
                OF = CONSIG*Tval[jT] + log(CONLAM*fabs(NUMft[j])/NOPTS);
                if (OF > log(DBL_MAX))
                {
                    IFAIL[j] = 1;
                    IFAIL_tot = 1; /* logical or of IFAIL(jX,jT) */
                }
                else
                    NUMft[j] = exp(OF)*SIGN;
            }
        } /* end for-loop j = jX*NTval+jT  to compute u(jX,jT) */

    return IFAIL_tot;
/**************************************************************************/
}



