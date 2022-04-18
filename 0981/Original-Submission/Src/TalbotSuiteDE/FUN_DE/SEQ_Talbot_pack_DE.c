/**********************    SEQ_Talbot_pack_DE.c    ************************
 *                                                                        *
 *                sequential version of TALBOT SUITE DE                   *
 *                                                                        *
 *     APPLICATION OF TALBOT'S METHODS TO SOLVE DIFFERENTIAL PROBLEMS     *
 *                                                                        *
 *           THE LAPLACE TRANSFORM IS GIVEN BY NUMERICAL SAMPLES          *
 *                             AND NOT AS A FUNCTION                      *
 *                                                                        *
 *                                                                        *
 *                          SEQUENTIAL IMPLEMENTATION                     *
 *                                                                        *
 *                            path: code/SRC/FUN_DE/                      *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>       VERSION 4.0     May 25th, 2016       <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *                       AUTHOR: Mariarosaria Rizzardi                    *
 *                                                                        *
 *                  mariarosaria.rizzardi@uniparthenope.it                *
 *                                                                        *
 *                  DiST - Dept. of Science and Technology                *
 *                  "Parthenope" University, Naples (Italy)               *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/

#include "../COM/main_include.h"
#include "../COM/COM_Talbot_pack.h"
#include "../COM_DE/COM_Talbot_pack_DE.h"
#include "SEQ_Talbot_pack_DE.h"


int SEQ_Talbot1_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol),
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax)
/**************************************************************************

    SEQ_Talbot1_DE   DRIVER FUNCTION (user level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides a numerical approximation to the Inverse Laplace
    Transform u(x,t) computed at each value of the  Xval,Tval  arrays.
    This is accomplished according to the modified Talbot method, described in:

    Rizzardi M. - "A modification of Talbot's method for the simultaneous
                   approximation of several values of the Inverse Laplace Transform".
                   ACM Trans. Math. Soft., vol. 21, no. 4, Dec. 1995, pp. 347-371.

    Algorithm's steps are:
        1) compute Talbot's parameters at (Tmin + Tmax)/2
           by means of COM_TalbotPAR function;

        2) compute Laplace Transform samples U(x,s) for s on
           Talbot's contour;

        3) for all x,t do
                approximate the Inverse Laplace Transform u(x,t)
                by means of SEQ_TalbotSUM1_DE function.

    COM_TalbotPAR() and SEQ_TalbotSUM1_DE() are skill-level functions:
    the former is from Talbot Suite (file COM_Talbot_pack.c).


    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_Talbot1_DE (LTsamples, sigma0, NXval, Xval, NTval, Tval,
                                    tol, NUMft, IFAIL, Nsings, SINGS, MULT,
                                    Tmin, Tmax);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTsamples - (double complex function pointer) user-defined function according
                to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol)
                The function returns the LT samples by solving the problem obtained
                by the application of the Laplace transform method to the original
                differential problem.

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

    SEQ_TalbotSUM1_DE: approximate the Inverse Laplace Transform.

    LTsamples     : (user-defined function) Laplace Transform function
                    according to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU;
    unsigned int NOPTS;
    int IFAIL_tot;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;

    /* 1) Compute Talbot's parameters at the middle point of the interval [Tmin, Tmax] */
    COM_TalbotPAR (sigma0,(Tmin+Tmax)/2,tol,Nsings,SINGS,MULT, &CONLAM,&CONSIG,&CONNU,&NOPTS);

    /* 2) Compute the correction to NOPTS */
    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4);

    /* 3) Compute sample points on Talbot's contour */
    double complex *S = (double complex*)malloc(NOPTS*sizeof(double complex));
    if ( S == NULL )
    {   fprintf(stderr, "\n***   ERROR in SEQ_Talbot1_DE: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");
        exit(1);
    }
    double thetak,   pi = 4.0*atan(1.0);
    unsigned int k;
    for (k=1; k<NOPTS; k++)
    {   thetak = (pi/NOPTS)*k;
        S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
    }
    S[0] = CONSIG + CONLAM;

    /* 4) Compute the matrix of Laplace Transform samples FS(i,j) on Talbot's contour **/
    double complex *FS = (*LTsamples)(NXval,Xval,NOPTS,S,tol);
    free(S);

    /* 5) Compute, for all x and for all t, u(x,t), local error flags
                   and a global error flag with the same LT samples */
    IFAIL_tot = SEQ_TalbotSUM1_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,Tval,NUMft,IFAIL);
    free(FS);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int SEQ_Talbot2_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol),
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[])
/**************************************************************************

    SEQ_Talbot2_DE   DRIVER FUNCTION (user level)

        IMPLEMENTATION OF CLASSICAL TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides numerical approximations to the Inverse Laplace
    Transform u(x,t) computed at each value of the  Xval,Tval  arrays.
    This is accomplished according to the classical Talbot method, described in:

        Talbot A. - "The accurate numerical inversion of Laplace Transforms".
                     J. Inst. Maths. Applics. (1979), n.23, pp.97-120.

        Murli A., Rizzardi M. - "Algorithm 682: Talbot's method for the
                                 Laplace inversion problem".
                                 ACM Trans. Math. Soft., vol. 16,
                                 no. 2, June 1990, pp. 158-168.

    Algorithm's sketch:

        for each t in Tval do

            1) compute Talbot's parameters, at t, by means of
               COM_TalbotPAR function;

            2) compute Laplace Transform samples U(x,s) for s on
               Talbot's contour;

            3) approximate the Inverse Laplace Transform u(x,t) for all x
               by means of SEQ_TalbotSUM2_DE function.


    COM_TalbotPAR() and SEQ_TalbotSUM2_DE() are skill-level functions:
    the former is from Talbot Suite (file code/SRC/COM/COM_Talbot_pack.c).

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_Talbot2_DE  (LTsamples, sigma0, NXval, Xval, NTval, Tval,
                                    tol, NUMft, IFAIL, Nsings, SINGS, MULT);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTsamples - (double complex function pointer) user-defined function according
                to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol)
                The function returns the LT samples by solving the problem obtained
                by the application of the Laplace transform method to the original
                differential problem.


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
                f(t) is approximated. Its components must be positive numbers.
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

    SEQ_TalbotSUM2_DE: approximate the Inverse Laplace Transform.

    LTsamples     : (user-defined function) Laplace Transform function
                    according to the following prototype:
                    double complex* (*LTsamples) (unsigned int NXval, double Xval[],
                                                  unsigned int NOPTS, double complex S[],
                                                  double tol)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU, TVALUE;
    int IFAIL_tot = 0;
    unsigned int NOPTS, jT;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;

    /* For each t */
    for (jT=0; jT<NTval; jT++)
    {
        TVALUE = *(Tval+jT);

        /* 1) Compute Talbot's parameters at each TVALUE */
        COM_TalbotPAR (sigma0, TVALUE, tol, Nsings, SINGS, MULT, &CONLAM, &CONSIG, &CONNU, &NOPTS);

        /* 2) Compute sample points on Talbot's contour at each TVALUE */
        double complex *S = (double complex*)malloc(NOPTS*sizeof(double complex));
        if ( S == NULL )
        {   fprintf(stderr, "\n***   ERROR in SEQ_Talbot2_DE: DYNAMIC ALLOCATION OF S IS FAILED. ***\n");
            exit(1);
        }
        double thetak,   pi = 4.0*atan(1.0);
        unsigned int k;
        for (k=1; k<NOPTS; k++)
        {   thetak = (pi/NOPTS)*k;
            S[k] = CONSIG + CONLAM*thetak/tan(thetak) + I*CONLAM*CONNU*thetak; /* in math.h cot(x) doesn't exist */
        }
        S[0] = CONSIG + CONLAM;

        /* 3) Compute the matrix of Laplace Transform samples FS(i,j) on Talbot's contour **/
        double complex *FS = (*LTsamples)(NXval,Xval,NOPTS,S,tol);
        free(S);

        /* 4) Compute, for all x and for TVALUE, u(x,TVALUE), local error flags
              and a global error flag with the LT samples */
        *(IFAIL+jT) = SEQ_TalbotSUM2_DE (CONLAM,CONSIG,CONNU,NOPTS,NXval,FS,NTval,TVALUE,jT,NUMft,IFAIL);
        free(FS);
        IFAIL_tot = IFAIL_tot || *(IFAIL+jT);
    }
    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int SEQ_TalbotSUM1_DE (double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[])
/**************************************************************************

    SEQ_TalbotSUM1_DE   SUMMATION FUNCTION (skill level)

        IMPLEMENTATION OF MODIFIED TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace Transform
    u(x,t) evaluated at each value of the  Xval,Tval  arrays.
    This is accomplished according to the modified Talbot method, described in:

        Rizzardi M. - "A modification of Talbot's method for the simultaneous
                       approximation of several values of the Inverse Laplace
                       Transform". ACM Trans. Math. Soft., vol. 21, no. 4,
                       December 1995, pp. 347-371.

    The composite Trapezoidal rule, approximating the contour integral for u(x,t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_TalbotSUM1_DE (CONLAM, CONSIG, CONNU, NOPTS, NXval, FF, NTval, Tval,
                                    NUMft, IFAIL);

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
                containing the Laplace Transform samples on the Talbot contour.

    NTval     - (unsigned integer) number of t values where the Inverse
                Laplace Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                u(x,t) is approximated. It must be dimensioned at least "NTval".

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
    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.
    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int NM1, K, IFAIL_tot;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN;

    IFAIL_tot=0; PIOVN=pi/(double)NOPTS; NM1=NOPTS - 1;

    unsigned int jT, jX; /* for loop indices */
    for (jX=0; jX<NXval; jX++)
        for (jT=0; jT<NTval; jT++)
        {
            NUMft[jX*NTval+jT]=0.; IFAIL[jX*NTval+jT]=0; /* NUMft and IFAIL are row-wise matrices of size (NXval,NTval) */
            TAU=CONLAM*Tval[jT];

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

                NUMft[jX*NTval+jT] = DBR - BR*U/2 - BI*sin(PSI); /* NUMft and IFAIL are row-wise matrices of size (NXval,NTval) */
            }

            /* compute   T0 (the first term in the summation) */
            NUMft[jX*NTval+jT] = NUMft[jX*NTval+jT] + CONNU*exp(TAU)*creal(FF[jX*NOPTS])/2.0;

            if ( NUMft[jX*NTval+jT] != 0.0 )
            {
                if ( NUMft[jX*NTval+jT] > 0.0 )
                    SIGN = +1.0;
                else
                    SIGN = -1.0;

                /* overflow   problem   handling */
                OF = CONSIG*Tval[jT] + log(CONLAM*fabs( NUMft[jX*NTval+jT] )/NOPTS);
                if (OF > log(DBL_MAX))
                    IFAIL[jX*NTval+jT] = 1;
                else
                    NUMft[jX*NTval+jT] = exp(OF)*SIGN;
            }

            IFAIL_tot = IFAIL_tot || IFAIL[jX*NTval+jT];

        } /* END for (jT=0; jT<NTval; jT++)
             END for (jX=0; jX<NXval; jX++) */

    return IFAIL_tot;
/**************************************************************************/
}


int SEQ_TalbotSUM2_DE (double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double TVALUE, unsigned int jT, double NUMft[], int IFAIL[])
/**************************************************************************

    SEQ_TalbotSUM2_DE   SUMMATION FUNCTION (skill level)

        IMPLEMENTATION OF CLASSICAL TALBOT'S METHOD FOR DIFFERENTIAL EQUATIONS
                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes the numerical approximation to the Inverse Laplace Transform
    u(x,t) evaluated at a single value of t (TVALUE) and for all the values x.
    This is accomplished according to the classical Talbot method, described in:

        Talbot A. - "The accurate numerical inversion of Laplace Transforms".
                     J. Inst. Maths. Applics. (1979), n.23, pp.97-120.

        Murli A., Rizzardi M. - "Algorithm 682: Talbot's method for the
                     Laplace Inversion problem". ACM Trans. Math. Soft.,
                     vol. 16, no. 2, June 1990, pp.158-168.

    The composite Trapezoidal rule, approximating the contour integral for u(x,t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_TalbotSUM2_DE (CONLAM, CONSIG, CONNU, NOPTS, NXval, FF, NTval,
                                       TVALUE, jT, NUMft, IFAIL);

    where IFAIL_tot is an error indicator computed as a logical "or" among all
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of addends in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    NXval     - (unsigned integer) number of x values where the Inverse Laplace
                Transform function is approximated.

    FF        - (double complex array) row-wise matrix, of size (NXval,NOPTS),
                containing the Laplace Transform samples on the Talbot contour.

    NTval     - (unsigned integer) number of t values where the Inverse Laplace
                Transform function is approximated.

    TVALUE    - (double) value for t where the Inverse Laplace Transform
                u(x,t) is going to be approximated.

    jT        - (unsigned integer) index corresponding to the current value of t
                (TVALUE). It locates a column in the output matrices, NUMft and
                IFAIL.

    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) row-wise matrix of size (NXval,NTval) containing
                the approximations to the values u(x(h),t(k)). Only the column
                related to jT is returned.

    IFAIL     - (integer array) row-wise matrix of size (NXval,NTval) containing
                the error flags at each u(x(h),t(k)). Only the column related to
                jT is returned:

                                / 0 no error
                  IFAIL(h,jT) = |
                                \ 1 an overflow occurs in u(x(h),t(jT)) so that,
                                    to avoid Inf as result, the returned value
                                    is scaled as
                                         u(x(h),t(jT)) = u(x(h),t(jT))/exp(sigma0*t(jT))

    REQUIRED FUNCTIONS
    ==================
    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.
    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int NM1, K, IFAIL_tot = 0;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN;


    PIOVN=pi/(double)NOPTS; TAU=CONLAM*TVALUE;
    NM1=NOPTS - 1;

    unsigned int jX;
    for (jX=0; jX<NXval; jX++)
    {
        IFAIL[jX*NTval+jT]=0; NUMft[jX*NTval+jT]=0.;
        if (NM1 > 0)
        {
            PSI=PIOVN*TAU*CONNU; C=cos(PSI);
            BR=0.; BI=0.; DBR=0.; DBI=0.;

            /* begin   Goertzel-Reinsch algorithm   for   Clenshaw sums  */
            if (C <= 0.) /* C=cos(PSI) <= 0 */
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
            else         /* C=cos(PSI) > 0 */
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
            /* end  Goertzel-Reinsch algorithm   for   Clenshaw sums  */
            NUMft[jX*NTval+jT] = DBR - BR*U/2 - BI*sin(PSI); /* NUMft and IFAIL are row-wise matrices of size (NXval,NTval) */
        }

        /* Computes   T0 (the first term in the summation) */
        NUMft[jX*NTval+jT] = NUMft[jX*NTval+jT] + CONNU*exp(TAU)*creal(FF[jX*NOPTS])/2.0;
        if (NUMft[jX*NTval+jT] != 0.0)
        {
            if (NUMft[jX*NTval+jT] > 0.0)
                SIGN = +1.0;
            else
                SIGN = -1.0;

            /* overflow   problem   handling */
            OF = CONSIG*TVALUE + log(CONLAM*fabs(NUMft[jX*NTval+jT])/NOPTS);
            if (OF > log(DBL_MAX))
                IFAIL[jX*NTval+jT] = 1;
            else
                NUMft[jX*NTval+jT] = exp(OF)*SIGN;
        }
        IFAIL_tot = IFAIL_tot || IFAIL[jX*NTval+jT];
    } /* END for (jX=0; jX<NXval; jX++) */

    return IFAIL_tot;
/**************************************************************************/
}

