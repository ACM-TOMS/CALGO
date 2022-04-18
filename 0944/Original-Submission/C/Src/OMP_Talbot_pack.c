/************************    OMP_Talbot_pack.c    *************************
 *                                                                        *
 *                              TALBOT SUITE                              *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *   SEQUENTIAL AND PARALLEL NUMERICAL INVERSION OF LAPLACE TRANSFORMS    *
 *                                                                        *
 *                                                                        *
 *              OMP-BASED IMPLEMENTATION OF TALBOT'S METHOD               *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>        VERSION 1.0    Sept 13th, 2012       <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *           AUTHORS: Laura Antonelli (1), Stefania Corsaro (2-1),        *
 *                    Zelda Marino (2), Mariarosaria Rizzardi (3)         *
 *                                                                        *
 *                   (1) ICAR  - National Research Council of Italy       *
 *                                                                        *
 *                   (2) DSMRE - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 *                   (3) DSA   - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * Antonelli L., Corsaro S.,                                              *
 * Marino Z., Rizzardi M. - "Talbot Suite: a parallel software collection *
 *                    for the numerical inversion of Laplace Transforms". *
 *                    ACM Trans. Math. Softw., vol.##, no.#, month year,  *
 *                    pp. ##-##.                                          *
 *                                                                        *
 **************************************************************************/


#ifdef _OPENMP
    #include <omp.h>
#endif

#include "main_include.h"
#include "COM_Talbot_pack.h"
#include "OMP_Talbot_pack.h"


int OMP_Talbot1 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                 double Tmin, double Tmax)
/**************************************************************************

    OMP-BASED IMPLEMENTATION: OMP_Talbot1 DRIVER FUNCTION (user level)

                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides numerical approximations to the Inverse Laplace Transform f(t)
    computed at each value of the  Tval  array.
    This is accomplished according to the modified Talbot method, described in:

        Rizzardi M. - "A modification of Talbot's method for the simultaneous
                       approximation of several values of the Inverse Laplace Transform".
                       ACM Trans. Math. Soft., vol. 21, no. 4, Dec. 1995, pp. 347-371.

    A coarse grain OpenMP-based parallelism is implemented; parallelization
    strategy is data distribution.

    Algorithm's sketch:

        1) compute Talbot's parameters at (Tmin + Tmax)/2
           by means of COM_TalbotPAR function;

        2) for all t in Tval do
                approximate the Inverse Laplace Transform f(t)
                by means of OMP_TalbotSUM1 function.

    COM_TalbotPAR() and OMP_TalbotSUM1() are skill-level functions
    implementing the modified Talbot method.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_Talbot1 (LTpt, sigma0, NTval, Tval, tol, NUMft, IFAIL,
                                 Nsings, SINGS, MULT, Tmin, Tmax);

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

    OMP_TalbotSUM1: to approximate the Inverse Laplace Transform.

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

    /* compute, for all t, f(t), local error flags and a global error flag */
    IFAIL_tot = OMP_TalbotSUM1 (LTpt, CONLAM, CONSIG, CONNU, NOPTS, NTval, Tval, NUMft, IFAIL);

    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int OMP_Talbot2 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[])
/**************************************************************************

    OMP-BASED IMPLEMENTATION: OMP_Talbot2 DRIVER FUNCTION (user level)

                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides numerical approximations to the Inverse Laplace Transform f(t)
    computed at each value of the  Tval  array.
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
               by means of OMP_TalbotSUM2 function.

    COM_TalbotPAR() and OMP_TalbotSUM2() are skill-level functions
    implementing the classical Talbot method.
    A fine grain OMP-based parallelism is implemented in OMP_TalbotSUM2();
    the summation process has been parallelized.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_Talbot2 (LTpt, sigma0, NTval, Tval, tol, NUMft, IFAIL,
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

    OMP_TalbotSUM2: to approximate the Inverse Laplace Transform.

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

        *(IFAIL+k) = OMP_TalbotSUM2 (LTpt, CONLAM, CONSIG, CONNU, NOPTS, TVALUE, NUMft+k);
        IFAIL_tot = IFAIL_tot || *(IFAIL+k);
    }
    return IFAIL_tot; /* global error flag */
/**************************************************************************/
}


int OMP_TalbotSUM1 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    unsigned int NTval, double *Tval, double *NUMft, int *IFAIL)
/**************************************************************************

    OMP-BASED IMPLEMENTATION: OMP_TalbotSUM1 SUMMATION FUNCTION (skill level)

                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes numerical approximations to the Inverse Laplace Transform f(t)
    evaluated at each value of the  Tval  array.
    This is accomplished according to the modified Talbot method, described in:

        Rizzardi M. - "A modification of Talbot's method for the simultaneous
                       approximation of several values of the Inverse Laplace
                       Transform". ACM Trans. Math. Soft., vol. 21, no. 4,
                       December 1995, pp. 347-371.

    A coarse grain parallelism (data distribution) is implemented.

    The composite Trapezoidal rule, approximating the contour integral for f(t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been implemented. This is accomplished for all t in the
    Tval array.

    CALLING SEQUENCE
    ================
        IFAIL_tot = OMP_TalbotSUM1 (LTpt, CONLAM, CONSIG, CONNU, NOPTS,
                                    NTval, Tval, NUMft, IFAIL);
    where IFAIL_tot is an error indicator computed as a logical or
    among the values of IFAIL.

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
    int THREADS, NM1, K, IFAIL_tot;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN;
    double complex S, FF;

    #ifdef _OPENMP
        THREADS = omp_get_max_threads(); /* activated threads */
    #else
        THREADS = 1; /* set sequential */
    #endif

    IFAIL_tot=0; PIOVN=pi/(double)NOPTS; NM1=NOPTS - 1;

    unsigned int jT, NTloc, mod, STARTloc, ENDloc; int myid, IFAIL_loc, IFAIL_t_loc;
    double ft_loc;

    #ifdef _OPENMP
        #pragma omp parallel default   (shared)                                                           \
                             private   (TAU,PSI,C,U,BR,BI,DBR,DBI,K,TETA,ALFA,BETA,S,FF,G,H,EAT,SIGN,OF)  \
                             private   (NTloc,mod,myid,STARTloc,ENDloc,IFAIL_loc,jT,ft_loc,IFAIL_t_loc)   \
                             reduction (|| : IFAIL_tot)                                                   \
                             num_threads (THREADS)
    #endif
    {
        /* local t values index distribution */
        NTloc=NTval/THREADS; mod=NTval%THREADS;

        #ifdef _OPENMP
            myid = omp_get_thread_num();
        #else
            myid = 0;
        #endif

        if (myid < mod)
        {
            NTloc = NTloc+1;
            STARTloc = myid*NTloc;
        }
        else
            STARTloc = myid*NTloc + mod;

        ENDloc = STARTloc + NTloc - 1;

        IFAIL_loc=0;

        /* loop on local values of t */
        for (jT=STARTloc; jT<=ENDloc; jT++)
        {
            ft_loc=0.0; IFAIL_t_loc=0;
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
                else /* C=cos(PSI) > 0 */
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
                /* end   Goertzel-Reinsch algorithm   for   Clenshaw sums  */

                ft_loc=DBR - BR*U/2 - BI*sin(PSI);
            }

            /* compute   T0 (the first term in the summation) */
            S = CONLAM + CONSIG + I*0.0;   FF = (*LTpt)( S );
            ft_loc = ft_loc + CONNU*exp(TAU)*creal(FF)/2.0;
            if ( ft_loc != 0.0 )
            {
                if ( ft_loc > 0.0 )
                    SIGN = +1.0;
                else
                    SIGN = -1.0;

                /* overflow   problem   handling */
                OF = CONSIG*Tval[jT] + log(CONLAM*fabs( ft_loc )/NOPTS);
                if (OF > log(DBL_MAX))
                    IFAIL_t_loc = 1;
                else
                    ft_loc = exp(OF)*SIGN;
            }

            NUMft[jT]=ft_loc; IFAIL[jT]=IFAIL_t_loc;
            IFAIL_loc = IFAIL_loc || IFAIL_t_loc;

        } /* END for (jT=STARTloc; jT<=ENDloc; jT++) */

        #ifdef _OPENMP
            #pragma omp barrier
        #endif
        IFAIL_tot = IFAIL_tot || IFAIL_loc; /* reduction for the global "or" */

    }

    return IFAIL_tot;
/**************************************************************************/
}


int OMP_TalbotSUM2 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    double TVALUE, double *NUMft)
/**************************************************************************

    OMP-BASED IMPLEMENTATION: OMP_TalbotSUM2 SUMMATION FUNCTION (skill level)

                              DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function computes the numerical approximation to the Inverse Laplace Transform
    f(t) evaluated at a single value of t (TVALUE).
    This is accomplished according to the classical Talbot method, described in:

        Talbot A. - "The accurate numerical inversion of Laplace Transforms".
                     J. Inst. Maths. Applics. (1979), n.23, pp.97-120.

        Murli A., Rizzardi M. - "Algorithm 682: Talbot's method for the
                     Laplace Inversion problem". ACM Trans. Math. Soft.,
                     vol. 16, no. 2, June 1990, pp.158-168.

    The composite Trapezoidal rule, approximating the contour integral for f(t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been parallelized with a fine grain parallel strategy,
    according to the algorithm described in:

        De Rosa, M.A. Giunta G., Rizzardi M. - "Parallel Talbot's algorithm
                     for distributed memory machines", Parallel Computing,
                     (1995), vol. 21, pp.783-801.

    CALLING SEQUENCE
    ================
        IFAIL = OMP_TalbotSUM2 (LTpt, CONLAM, CONSIG, CONNU, NOPTS,
                                TVALUE, NUMft);

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

    omp_get_max_threads, omp_get_thread_num: OpenMP runtime library functions.

    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.

    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int THREADS, K, IFAIL;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN, sum;
    double complex S, FF;

    #ifdef _OPENMP
        THREADS = omp_get_max_threads(); /* activated threads */
    #else
        THREADS = 1; /* set sequential */
    #endif

    IFAIL=0; *NUMft=0.;
    PIOVN=pi/(double)NOPTS; TAU=CONLAM*TVALUE; sum=0.0;

    if (NOPTS > 1)
    {
        unsigned int Nloc, mod, myid, STARTloc, ENDloc; double sum_loc;
        PSI=PIOVN*TAU*CONNU; C=cos(PSI);

        #ifdef _OPENMP
            #pragma omp parallel default   (shared)                                        \
                                 private   (BR,BI,DBR,DBI,U,K,TETA,ALFA,BETA,S,FF,G,H,EAT) \
                                 private   (Nloc,mod,myid,STARTloc,ENDloc,sum_loc)         \
                                 reduction (+ : sum)                                       \
                                 num_threads (THREADS)
        #endif
        {
            /* local summation index distribution */
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
                    S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                    FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
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
                    S = CONLAM*ALFA+CONSIG + I*CONLAM*CONNU*TETA;
                    FF = (*LTpt)( S ); G=creal(FF); H=cimag(FF);
                    EAT = exp(ALFA*TAU);
                    BR = DBR + BR;
                    BI = DBI + BI;
                    DBR = U*BR + DBR + EAT*(G*CONNU - H*BETA);
                    DBI = U*BI + DBI + EAT*(H*CONNU + G*BETA);
                }
            }

            DBR = DBR - BR*U/2;
            DBI = DBI - BI*U/2;
            BR = BR*sin(PSI);;
            BI = BI*sin(PSI);

            sum_loc = DBR*cos(STARTloc*PSI) - BR*sin(STARTloc*PSI) - BI*cos(STARTloc*PSI) - DBI*sin(STARTloc*PSI);
            /* end   Goertzel-Reinsch algorithm   for   Clenshaw sums  */

            #ifdef _OPENMP
                #pragma omp barrier
            #endif
            sum += sum_loc; /* reduction for the global sum */
        }
    }

    /* Master thread computes T0 (the first term in the summation) */
    S = CONLAM + CONSIG + I*0.0;
    FF = (*LTpt)( S );
    *NUMft = sum + CONNU*exp(TAU)*creal(FF)/2.0;
    if (*NUMft != 0.0)
    {
        if (*NUMft > 0.0)
            SIGN = +1.0;
        else
            SIGN = -1.0;

        /* overflow   problem   handling */
        OF = CONSIG*TVALUE + log(CONLAM*fabs(*NUMft)/NOPTS);
        if (OF > log(DBL_MAX))
            IFAIL = 1;
        else
            *NUMft = exp(OF)*SIGN;
    }

    return IFAIL;
/**************************************************************************/
}

