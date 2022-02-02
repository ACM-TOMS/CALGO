/************************    SEQ_Talbot_pack.c    *************************
 *                                                                        *
 *                  sequential version of TALBOT SUITE                    *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *             THE NUMERICAL INVERSION OF LAPLACE TRANSFORMS              *                                                                        *
 *                                                                        *
 *                           SEQUENTIAL PACKAGE                           *
 *                                  FOR                                   *
 *                     TALBOT'S METHOD IMPLEMENTATION                     *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>        VERSION 2.0    May 25, 2016          <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *           AUTHOR: Mariarosaria Rizzardi                                *
 *                                                                        *
 *                   DiST - Department of Science and Technology          *
 *                          "Parthenope" University, Naples (Italy)       *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * L. ANTONELLI, S. CORSARO, Z. MARINO,                                   *
 * M. RIZZARDI: "Algorithm 944: TALBOT SUITE: A PARALLEL IMPLEMENTATIONS  *
 *                              OF TALBOT'S METHOD FOR THE NUMERICAL      *
 *                              INVERSION OF LAPLACE TRANSFORMS".         *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. 40,    *
 *                                 NO. 4, June 2014, (29:2--29:18).       *
 *                                                                        *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP).               *
 **************************************************************************/


#include "../../../TalbotSuiteDE/COM/main_include.h"
#include "../../../TalbotSuiteDE/COM/COM_Talbot_pack.h"
#include "SEQ_Talbot_pack.h"



int SEQ_Talbot1 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                 double Tmin, double Tmax)
/**************************************************************************

    SEQUENTIAL PACKAGE: SEQ_Talbot1 DRIVER FUNCTION (user level)

                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides a numerical approximation to the Inverse Laplace
    Transform f(t) computed at each value of the  Tval  array.
    This is accomplished according to the modified Talbot method, described in:

    Rizzardi M. - "A modification of Talbot's method for the simultaneous
                   approximation of several values of the Inverse Laplace Transform".
                   ACM Trans. Math. Soft., vol. 21, no. 4, Dec. 1995, pp. 347-371.

    Algorithm's steps are:
        1) compute Talbot's parameters at (Tmin + Tmax)/2
           by means of COM_TalbotPAR function;

        2) for all t in Tval do
                approximate the Inverse Laplace Transform f(t)
                by means of SEQ_TalbotSUM1 function.

    COM_TalbotPAR() and SEQ_TalbotSUM1() are skill-level functions
    implementing the modified Talbot method.

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_Talbot1 (LTpt, sigma0, NTval, Tval, tol, NUMft, IFAIL,
                                 Nsings, SINGS, MULT);

    where IFAIL_tot is an error indicator computed as a logical "or" among
    the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTpt      - (double complex function pointer) pointer to the Laplace
                Transform function to be inverted.
                It is a user defined function.

    sigma0    - (double) abscissa of (absolute) convergence of the
                Laplace Transform function.

    NTval     - (unsigned integer) number of t values where the Inverse
                Laplace Transform function is approximated.

    Tval      - (double array) values for t where the Inverse
                Laplace Transform f(t) is approximated.
                Its components must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |f(t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays SINGS and MULT.

    SINGS     - (double complex array) singularities of the Laplace
                Transform function. Only singularities with non-negative
                imaginary parts are required; their complex conjugates are
                unnecessary.
                It must be dimensioned at least "Nsings".

    MULT      - (unsigned integer array) multiplicities of those
                singularities (in SINGS) which are poles, zero otherwise.
                It must be dimensioned as SINGS.

    Tmin,Tmax - (double) endpoints of the interval enclosing the t values.
                Method's parameters are computed at (Tmin + Tmax)/2.

    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) approximations to the values f(t) for t in Tval.
                It must be dimensioned at least "NTval".

    IFAIL     - (integer array) error flags at each t in Tval:

                                / 0 no error
                    IFAIL[k] = |
                                \ 1 an overflow occurs in f(t) so that, to
                                    avoid Inf as result, the returned value
                                    is scaled as
                                        NUMft = f(t)/exp(sigma0*t)

                It must be dimensioned at least "NTval".

    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : to compute the Talbot parameters (in COM_Talbot_pack.c).

    SEQ_TalbotSUM1: to approximate the Inverse Laplace Transform.

    LTpt          : (user defined function) Laplace Transform function
                    according to the following prototype:
                        double complex LTpt (double complex s)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU;
    int IFAIL_tot;
    unsigned int NOPTS;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;

    /* Talbot's parameters at the middle point of [Tmin, Tmax] */
    COM_TalbotPAR (sigma0, (Tmin + Tmax)/2, tol, Nsings, SINGS, MULT, &CONLAM, &CONSIG, &CONNU, &NOPTS);

    /* Does not compute the correction to NOPTS */
/*    NOPTS = COM_TalbotNcorr(Tmin,Tmax,sigma0,CONLAM,CONSIG,CONNU,NOPTS,1e-4); */

    /* compute, for all t, f(t), local error flags and a global error flag */
    IFAIL_tot = SEQ_TalbotSUM1 (LTpt, CONLAM, CONSIG, CONNU, NOPTS, NTval, Tval, NUMft, IFAIL);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int SEQ_Talbot2 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[])
/**************************************************************************

    SEQUENTIAL PACKAGE: SEQ_Talbot2 DRIVER FUNCTION (user level)

                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides numerical approximations to the Inverse Laplace
    Transform f(t) computed at each value of the  Tval  array.
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

            2) approximate the Inverse Laplace Transform f(t)
               by means of SEQ_TalbotSUM2 function.

    COM_TalbotPAR() and SEQ_TalbotSUM2() are skill-level functions
    implementing the classical Talbot method.

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_Talbot2 (LTpt, sigma0, NTval, Tval, tol, NUMft, IFAIL,
                                 Nsings, SINGS, MULT);

    where IFAIL_tot is an error indicator computed as a logical or
    among the values of IFAIL.

    INPUT PARAMETERS
    ================
    LTpt      - (double complex function pointer) pointer to the Laplace
                Transform function to be inverted.
                It is a user defined function.

    sigma0    - (double) abscissa of (absolute) convergence of the
                Laplace Transform function.

    NTval     - (unsigned integer) number of t values where the Inverse
                Laplace Transform function is approximated.

    Tval      - (double array) values for t where the Inverse
                Laplace Transform f(t) is approximated.
                Its components must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |f(t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays SINGS and MULT.

    SINGS     - (double complex array) singularities of the Laplace
                Transform function. Only singularities with non-negative
                imaginary parts are required; their complex conjugates are
                unnecessary.
                It must be dimensioned at least "Nsings".

    MULT      - (unsigned integer array) multiplicities of those
                singularities (in SINGS) which are poles, zero otherwise.
                It must be dimensioned as SINGS.

    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) approximations to the values f(t) for t in Tval.
                It must be dimensioned at least "NTval".

    IFAIL     - (integer array) error flags at each t in Tval:

                                / 0  no error
                    IFAIL[k] = |
                                \ 1  an overflow occurs in f(t) so that, to
                                     avoid Inf as result, the returned value
                                     is scaled as
                                        NUMft = f(t)/exp(sigma0*t)

                It must be dimensioned at least "NTval".

    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : to compute the Talbot parameters (in COM_Talbot_pack.c).

    SEQ_TalbotSUM2: to approximate the Inverse Laplace Transform.

    LTpt          : (user defined function) Laplace Transform function
                    according to the following prototype:
                        double complex LTpt (double complex s)

 **************************************************************************/
{
    double CONLAM, CONSIG, CONNU, TVALUE;
    int IFAIL_tot = 0;
    unsigned int NOPTS, k;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;

    /* compute, for each t, Talbot's parameters, f(t) and a local error flag */
    for (k=0; k<NTval; k++)
    {
        TVALUE = *(Tval+k);
        COM_TalbotPAR (sigma0, TVALUE, tol, Nsings, SINGS, MULT, &CONLAM, &CONSIG, &CONNU, &NOPTS);

        *(IFAIL+k) = SEQ_TalbotSUM2 (LTpt, CONLAM, CONSIG, CONNU, NOPTS, TVALUE, NUMft+k);
        IFAIL_tot = IFAIL_tot || *(IFAIL+k);
    }
    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int SEQ_TalbotSUM1 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    unsigned int NTval, double *Tval, double *NUMft, int *IFAIL)
/**************************************************************************

    SEQUENTIAL PACKAGE: SEQ_TalbotSUM1 SUMMATION FUNCTION (skill level)

                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace
    Transform f(t) evaluated at each value of the  Tval  array.
    This is accomplished according to the modified Talbot method, described in:

        Rizzardi M. - "A modification of Talbot's method for the simultaneous
                       approximation of several values of the Inverse Laplace
                       Transform". ACM Trans. Math. Soft., vol. 21, no. 4,
                       December 1995, pp. 347-371.

    The composite Trapezoidal rule, approximating the contour integral for f(t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented. This is accomplished for all t in the
    Tval array.

    CALLING SEQUENCE
    ================
        IFAIL_tot = SEQ_TalbotSUM1 (LTpt, CONLAM, CONSIG, CONNU, NOPTS,
                                    NTval, Tval, NUMft, IFAIL);

    INPUT PARAMETERS
    ================
    LTpt      - (double complex function pointer) pointer to the Laplace
                Transform function to be inverted.
                It is a user defined function.

    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of terms in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    NTval     - (unsigned integer) number of t values where the Inverse
                Laplace Transform function is approximated.

    Tval      - (double array) values for t where the Inverse Laplace Transform
                f(t) is approximated. It must be dimensioned at least "NTval".

    OUTPUT PARAMETERS
    =================
    NUMft     - (double array) approximations to the values f(t) for t in Tval.
                It must be dimensioned at least "NTval".

    IFAIL     - (integer array) error flags at each t in Tval:

                                / 0  no error
                    IFAIL[k] = |
                                \ 1  an overflow occurs in f(t) so that, to
                                     avoid Inf as result, the returned value
                                     is scaled as
                                         NUMft[k] = f(t)/exp(sigma0*t)

                It must be dimensioned at least "NTval".

    REQUIRED FUNCTIONS
    ==================
    LTpt        : (user defined function) Laplace Transform function
                  according to the following prototype:
                        double complex LTpt (double complex s)

    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.

    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int NM1, K, IFAIL_tot;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN;
    double complex S, FF;

    IFAIL_tot=0; PIOVN=pi/(double)NOPTS; NM1=NOPTS - 1;

    unsigned int jT;
    for (jT=0; jT<NTval; jT++)
    {
        NUMft[jT]=0.; IFAIL[jT]=0;
        TAU=CONLAM*Tval[jT];

        if (NM1 > 0)
        {
            PSI=PIOVN*TAU*CONNU; C=cos(PSI);
            BR=0.; BI=0.; DBR=0.; DBI=0.;

            /* BEGIN  GOERTZEL-REINSCH ALGORITHM   FOR   CLENSHAW SUMS  */
            if (C <= 0.)
            {
                U = +4.*pow(cos(PSI/2),2);
                for (K=NM1; K>=1; K--)
                {
                    TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                    BETA=TETA+ALFA*(ALFA-1.)/TETA;
                    S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                    FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
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
            else /* C=COS(PSI) > 0 */
            {
                U = -4.*pow(sin(PSI/2),2);
                for (K=NM1; K>=1; K--)
                {
                    TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                    BETA=TETA+ALFA*(ALFA-1.)/TETA;
                    S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                    FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
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
            /* END  GOERTZEL-REINSCH ALGORITHM   FOR   CLENSHAW SUMS  */
            *(NUMft+jT)=DBR - BR*U/2 - BI*sin(PSI);
        }

        /* COMPUTES   T0 (the first term in the summation) */
        S = CONLAM + CONSIG + I*0.0;   FF = (*LTpt)( S );
        *(NUMft+jT) = *(NUMft+jT) + CONNU*exp(TAU)*creal(FF)/2.0;

        if ( *(NUMft+jT) != 0.0 )
        {
            if ( *(NUMft+jT) > 0.0 )
                SIGN = +1.0;
            else
                SIGN = -1.0;

            /* OVERFLOW   PROBLEM   HANDLING */
            OF = CONSIG*Tval[jT] + log(CONLAM*fabs( *(NUMft+jT) )/NOPTS);
            if (OF > log(DBL_MAX))
                *(IFAIL+jT) = 1;
            else
                *(NUMft+jT) = exp(OF)*SIGN;
        }

        IFAIL_tot = IFAIL_tot || IFAIL[jT];

    } /* END for (jT=0; jT<NTval; jT++) */

    return IFAIL_tot;
/**************************************************************************/
}


int SEQ_TalbotSUM2 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    double TVALUE, double *NUMft)
/**************************************************************************

    SEQUENTIAL PACKAGE: SEQ_TalbotSUM2 SUMMATION FUNCTION (skill level)

                        DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes the numerical approximation to the Inverse Laplace
    Transform f(t) evaluated at a single value of t (TVALUE).
    This is accomplished according to the classical Talbot method, described in:

        Talbot A. - "The accurate numerical inversion of Laplace Transforms".
                     J. Inst. Maths. Applics. (1979), n.23, pp.97-120.

        Murli A., Rizzardi M. - "Algorithm 682: Talbot's method for the
                     Laplace Inversion problem". ACM Trans. Math. Soft.,
                     vol. 16, no. 2, June 1990, pp.158-168.

    The composite Trapezoidal rule, approximating the contour integral for f(t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented.

    CALLING SEQUENCE
    ================
        IFAIL = SEQ_TalbotSUM2 (LTpt, CONLAM, CONSIG, CONNU, NOPTS,
                                TVALUE, &NUMft);

    where IFAIL is an error indicator in f(t)

                              / 0  no error
                    IFAIL  = |
                              \ 1  an overflow occurs in f(t) so that, to
                                   avoid Inf as result, the returned value
                                   is scaled as
                                       NUMft = f(t)/exp(sigma0*t)

    INPUT PARAMETERS
    ================
    LTpt      - (double complex function pointer) pointer to the Laplace
                Transform function to be inverted.
                It is a user defined function.

    CONLAM, CONSIG, CONNU - (double) geometrical parameters for the Talbot
                integration contour (respectively lambda, sigma and nu in
                Talbot's original paper). Their values may be computed by
                means of the COM_TalbotPAR function.

    NOPTS     - (unsigned integer) number of points required by the quadrature
                rule, i.e. the number of addends in the Clenshaw sum. Its
                value may be computed by means of the COM_TalbotPAR function.

    TVALUE    - (double) value for t where the Inverse Laplace Transform
                f(t) is approximated.

    OUTPUT PARAMETERS
    =================
    NUMft     - (pointer to double) approximation to the Inverse Laplace
                Transform f(t) computed at t=TVALUE.

    REQUIRED FUNCTIONS
    ==================
    LTpt        : (user defined function) Laplace Transform function
                  according to the following prototype:
                        double complex LTpt (double complex s)

    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.

    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int NM1, K, IFAIL;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN;
    double complex S, FF;

    IFAIL=0; *NUMft=0.;
    PIOVN=pi/(double)NOPTS; TAU=CONLAM*TVALUE;
    NM1=NOPTS - 1;

    if (NM1 > 0)
    {
        PSI=PIOVN*TAU*CONNU; C=cos(PSI);
        BR=0.; BI=0.; DBR=0.; DBI=0.;

        /* BEGIN  GOERTZEL-REINSCH ALGORITHM   FOR   CLENSHAW SUMS  */
        if (C <= 0.)
        {
            U = +4.*pow(cos(PSI/2),2);
            for (K=NM1; K>=1; K--)
            {
                TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                BETA=TETA+ALFA*(ALFA-1.)/TETA;
                S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
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
        else  /* C=COS(PSI) > 0 */
        {
            U = -4.*pow(sin(PSI/2),2);
            for (K=NM1; K>=1; K--)
            {
                TETA=K*PIOVN; ALFA=TETA*cos(TETA)/sin(TETA);
                BETA=TETA+ALFA*(ALFA-1.)/TETA;
                S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
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
        /* END  GOERTZEL-REINSCH ALGORITHM   FOR   CLENSHAW SUMS  */
        *NUMft=DBR - BR*U/2 - BI*sin(PSI);
    }

    /* COMPUTES   T0 (the first term in the summation) */
    S = CONLAM + CONSIG + I*0.0;
    FF = (*LTpt)( S );
    *NUMft = *NUMft + CONNU*exp(TAU)*creal(FF)/2.0;
    if (*NUMft != 0.0)
    {
        if (*NUMft > 0.0)
            SIGN = +1.0;
        else
            SIGN = -1.0;

        /* OVERFLOW   PROBLEM   HANDLING */
        OF = CONSIG*TVALUE + log(CONLAM*fabs(*NUMft)/NOPTS);
        if (OF > log(DBL_MAX))
            IFAIL = 1;
        else
            *NUMft = exp(OF)*SIGN;
    }

    return IFAIL;
/**************************************************************************/
}

