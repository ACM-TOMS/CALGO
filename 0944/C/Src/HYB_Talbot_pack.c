/************************    HYB_Talbot_pack.c    *************************
 *                                                                        *
 *                              TALBOT SUITE                              *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *   SEQUENTIAL AND PARALLEL NUMERICAL INVERSION OF LAPLACE TRANSFORMS    *
 *                                                                        *
 *                                                                        *
 *           MPI/OMP-BASED HYB IMPLEMENTATION OF TALBOT'S METHOD          *
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

#include "mpi.h"

#include "main_include.h"
#include "HYB_Talbot_pack.h"
#include "COM_Talbot_pack.h"


int HYB_Talbot3 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval,  double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                 double Tmin, double Tmax)
/**************************************************************************

    MPI/OMP-BASED HYB IMPLEMENTATION: HYB_Talbot3 DRIVER FUNCTION (user level)

                          DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides numerical approximations to the Inverse Laplace Transform f(t)
    computed at each value of the  Tval  array.
    This is accomplished according to the modified Talbot method, described in:

        Rizzardi M. - "A modification of Talbot's method for the simultaneous
                       approximation of several values of the Inverse Laplace
                       Transform". ACM Trans. Math. Soft., vol. 21, no. 4,
                       Dec. 1995, pp. 347-371.

    Algorithm's sketch:

        1) compute Talbot's parameters at (Tmin + Tmax)/2
           by means of COM_TalbotPAR function;

        2) for all t in Tval do
                approximate the Inverse Laplace Transform f(t)
                by means of HYB_TalbotSUM3 function.

    COM_TalbotPAR() and HYB_TalbotSUM3() are skill-level functions
    implementing the modified Talbot method.

    A hybrid MPI/OpenMP parallel algorithm is implemented by using MPI
    for the coarse grain parallelism (data distribution) and by using
    OMP for the fine grain parallelism (summation parallelization).

    HYB_Talbot3() is a node function, called by each process activated in the
    MPI communicator; all the variables are local so that no communication is required.
    Data distribution has to be performed before calling this function.

    CALLING SEQUENCE
    ================
        IFAIL_tot = HYB_Talbot3 (LTpt, sigma0, NTval, Tval, tol, NUMft, IFAIL,
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
                They must be positive numbers.
                It must be dimensioned at least "NTval".

    tol       - (double) tolerance to the error in the result,
                in terms of absolute or relative error as follows:
                       absolute error <= tol   if   |f(t)| <= 1
                or
                       relative error <= tol   otherwise.

    Nsings    - (unsigned integer) size of the arrays: SINGS and MULT.

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

                                / 0  no error
                    IFAIL[k] = |
                                \ 1  an overflow occurs in f(t) so that, to
                                     avoid Inf as result, the returned value
                                     is scaled as
                                        NUMft = f(t)/exp(sigma0*t)

                It must be dimensioned at least "NTval".

    REQUIRED FUNCTIONS
    ==================
    COM_TalbotPAR : to compute the talbot parameters (in COM_Talbot_pack.c).

    HYB_TalbotSUM3: to approximate the Inverse Laplace Transform.

    LTpt          : (user defined function) Laplace Transform function
                    according to the following prototype:
                        double complex LTpt (double complex s)

 **************************************************************************/
{
    double CONSIG, CONLAM, CONNU;
    int IFAIL_tot = 0;
    unsigned int NOPTS;

    /* default value if tol is wrong */
    if (tol <= 0)   tol=1e-8;

    /* Talbot's parameters at the middle point of [Tmin, Tmax] */
    COM_TalbotPAR (sigma0, (Tmin + Tmax)/2, tol, Nsings, SINGS, MULT, &CONLAM, &CONSIG, &CONNU, &NOPTS);

    /* compute, for all t: f(t), local error flags and a global error flag */
    IFAIL_tot = HYB_TalbotSUM3 (LTpt, CONLAM, CONSIG, CONNU, NOPTS, NTval, Tval, NUMft, IFAIL);

    return IFAIL_tot;
}


int HYB_TalbotSUM3 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    unsigned int NTval, double *Tval, double *NUMft, int *IFAIL)
/**************************************************************************

    MPI/OMP-BASED HYB IMPLEMENTATION: HYB_TalbotSUM3 SUMMATION FUNCTION (skill level)

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

    This is a node function, called by each process activated in the MPI communicator.

    The composite Trapezoidal rule, approximating the contour integral for f(t),
    leads to the real part of a complex Clenshaw sum. In order to compute it
    the Goertzel algorithm, in the Reinsch stable version and in double precision
    real arithmetic, has been parallelized with a fine grain OMP-based parallel
    strategy, according to the algorithm described in:

        De Rosa, M.A. Giunta G., Rizzardi M. - "Parallel Talbot's algorithm
                     for distributed memory machines", Parallel Computing,
                     (1995), vol. 21, pp.783-801.

    The fine grain parallel algorithm is repeated for each value of t in Tval.

    CALLING SEQUENCE
    ================
        IFAIL_tot = HYB_TalbotSUM3 (LTpt, CONLAM, CONSIG, CONNU, NOPTS,
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

    omp_get_max_threads, omp_get_thread_num: OpenMP runtime library functions.

    abs, atan, cos, exp, fabs, log, pow, sin: math intrinsic functions.

    cimag, creal: complex intrinsic functions.

 **************************************************************************/
{
    int THREADS, K, jT, IFAIL_tot;
    double pi = 4.*atan(1.0);
    double PIOVN, TAU, PSI, C, BR, BI, DBR, DBI, U, TETA, ALFA, BETA;
    double G, H, EAT, OF, SIGN, sum;
    double complex S, FF;

    #ifdef _OPENMP
        THREADS = omp_get_max_threads(); /* activated threads */
    #else
        THREADS = 1; /* set sequential */
    #endif

    IFAIL_tot=0;  PIOVN=pi/(double)NOPTS;

    /* compute, for all t: f(t), local error flags and a global error flag */
    for (jT=0; jT<NTval; jT++)
    {
        NUMft[jT]=0.; IFAIL[jT]=0.;
        TAU=CONLAM*Tval[jT]; sum=0.0;

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
        *(NUMft+jT) = sum + CONNU*exp(TAU)*creal(FF)/2.0;

        if ( *(NUMft+jT) != 0.0)
        {
            if ( *(NUMft+jT) > 0.0)
                SIGN = +1.0;
            else
                SIGN = -1.0;

            /* overflow   problem   handling */
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
