/************************    COM_Talbot_pack.c    *************************
 *                                                                        *
 *                              TALBOT SUITE                              *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *   SEQUENTIAL AND PARALLEL NUMERICAL INVERSION OF LAPLACE TRANSFORMS    *
 *                                                                        *
 *                                                                        *
 *                        SHARED UTILITY FUNCTIONS                        *
 *                                   FOR                                  *
 *                    TALBOT'S METHOD IMPLEMENTATIONS                     *
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


#include "main_include.h"
#include "COM_Talbot_pack.h"


void COM_TalbotPAR (double sigma0, double TVALUE, double tol,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                    double *CONLAM, double *CONSIG, double *CONNU, unsigned int *NOPTS)
/**************************************************************************

    SHARED UTILITY PACKAGE: COM_TalbotPAR FUNCTION (skill level)

                            DOUBLE PRECISION VERSION

    PURPOSE
    =======
    This function provides values of the contour parameters
                (lambda, sigma, nu in Talbot's original paper)
    and of the accuracy parameter (N) according to Talbot's method for
    the numerical inversion of Laplace Transforms.
    These values can be used in any summation module of Talbot Suite
    in order to estimate, by a contour integration, the Inverse Laplace
    Transform computed at TVALUE.

    CALLING SEQUENCE
    ================
        COM_TalbotPAR (sigma0, TVALUE, tol, Nsings, SINGS, MULT,
                       CONLAM, CONSIG, CONNU, NOPTS);

    INPUT PARAMETERS
    ================
    sigma0    - (double) abscissa of (absolute) convergence of the
                Laplace Transform function to be inverted.

    TVALUE    - (double) value for t where the Inverse Laplace Transform
                f(t) is approximated. it must be a positive number.

    tol       - (double) tolerance to the error in the result, in terms
                of absolute or relative error as follows:
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
    CONLAM, CONSIG, CONNU - (pointer to double) geometrical parameters for the
                Talbot integration contour (respectively lambda, sigma and nu
                in Talbot's original paper).

    NOPTS     - (pointer to unsigned integer) number of points required by
                the quadrature rule, i.e. the number of addends in the
                Clenshaw sum.

    REQUIRED FUNCTIONS
    ==================
    COM_TalbotINV : auxiliary function to approximate the "principal
                inverse" applying a suitable real Newton process.

    atan, atan2, ceil, log10, pow, tan: math intrinsic functions.

    cimag, creal: complex intrinsic functions.

    min, max    : macros (defined in COM_Talbot_pack.h).

    REMARKS
    =======
    1      Numbers in brackets, in the comments, refer to formulas in the
           Talbot's original paper:
              Talbot A. - "The accurate numerical inversion of Laplace
                           Transforms".
                           J. Inst. Maths. Applics. (1979), n.23, pp.97-120.

    2      The abscissa of convergence sigma0 is a general parameter of
           all the numerical inversion methods and it is usually required
           among the input parameters; for this reason it has been
           inserted here too, although it could be computed starting from
           the singularities (SINGS) of the Laplace Transform function.

    3      The Talbot machine precision constant (IC), in CALGO #682, is
           defined as 3/4 of the equivalent decimal precision for mantissas
           of the floating-point arithmetic system.
           According to the IEEE Standard 754 for binary floating-point
           arithmetic, it may be computed as
                    IC = (unsigned int)(p+1)*log10(2);
           where  p  is the number of significant bits for mantissas.
           For portability it may be also computed by means of the machine
           epsilon (DBL_EPSILON in float.h) as
                    IC = floor( -log10(DBL_EPSILON/2)*3./4. );
           Here, for simplicity, the same value is obtained by means of
           DBL_DIG (in float.h) simply by
                    IC = (unsigned int)DBL_DIG*0.75;

 **************************************************************************/
{
    unsigned char CASE1;
    int DD, KMAX, J, L, ID, MD=0, N0=0, N1=0, N2=0, ICM1, N0L;
    double OMEGA=0.0, PMAX=0.0, RD=0.0, SID=0.0, TETD=0.0;
    double ETAJ=0.0, TETJ=0.0, RJ=0.0, SRD=0.0, V=0.0, PK=0.0, FI=0.0, PMU=0.0;
    double E, UNRO, P, Q, U, GAMM, ETA, Y;

    double pi = 4.0*atan(1.0);
    unsigned int IC, DECDIG;

    IC = (unsigned int)DBL_DIG*0.75;
    CASE1=1; /* true */
    OMEGA=0.4*(IC+1);

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       compute
       "PMAX" = maximum real parts of singularities                    (58)
       "KMAX" = array-index of the singularity having the largest
                multiplicity
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    PMAX=creal(SINGS[0]); KMAX=0;
    for (J=1; J<Nsings; J++)
    {
        if (creal(SINGS[J]) > PMAX)  PMAX=creal(SINGS[J]);
        if (MULT[J] > MULT[KMAX]) KMAX=J;
    }

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       compute the "dominant singularity"

              SD = SRD + I*SID         ("I" is the imaginary unit)

       such that   im(SD)/arg(SD) is maximum over all the complex
       singularities (shifted by "sigma0") with im(SD) > 0.              (63)

       otherwise   im(SD) is set to zero and arg(SD) is set to "pi"
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
    RD=0.; ID=0; SID=0.; TETD=pi;
    for (J=0; J<Nsings; J++)
        if (cimag(SINGS[J]) > 0.0)
        {
            CASE1=0; /* false */
            ETAJ=atan2(-creal(SINGS[J]) + sigma0,  cimag(SINGS[J]));
            TETJ=ETAJ + pi/2.;
            RJ=cimag(SINGS[J])/TETJ;
            if (RD < RJ)
            {
                RD=RJ;
                TETD=TETJ;
                ID=J;
            }
        }

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       for an appropriate choice of the geometrical parameters:
       in the case of real singularities, always "Case 1" holds.
       otherwise it computes the shifted "dominant singularity"
       and the quantities "V" and "OMEGA".                    (65) and (66)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if ( CASE1 )
    {
        *CONLAM = OMEGA/TVALUE;
        *CONSIG = sigma0;
        *CONNU  = 1.0;
    }
    else /* CASE1 is false */
    {
        SRD = creal(SINGS[ID]) - sigma0;
        SID = cimag(SINGS[ID]);
        MD  = MULT[ID];
        V   = SID*TVALUE;
        OMEGA = min(OMEGA + V/2.,  2.*(IC+1)/3.);

        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            it checks for "Case 1" or "Case 2"
            when singularities are complex
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
        if (1.8*V <= OMEGA*TETD)
        {
            /* geometrical parameters for "case 1"    (68) */
            CASE1=1;
            *CONLAM = OMEGA/TVALUE;
            *CONSIG = sigma0;
            *CONNU  = 1.0;
        }
        else
        {
            /* geometrical parameters for "Case 2"    (74) */
            CASE1 = 0;
            PK  = 1.6 + 12./(V+25.);
            FI  = 1.05 + 1050./(max(553.0, 800.-V));
            PMU = (OMEGA/TVALUE + sigma0 - PMAX)/(PK/FI - 1./tan(FI));
            *CONLAM = PK*PMU/FI;
            *CONSIG = PMAX-PMU/tan(FI);
            *CONNU  = SID/PMU;
        }
    }

    /* required accurate decimal digits in the result */
    DECDIG = ceil(-log10(tol));

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        compute   NOPTS = max (N0, N1, N2)     (78)
        at first "N1" is computed              (87)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
    E = (2.3*DECDIG+OMEGA)/(*CONLAM*TVALUE);

    if (E <= 4.4)
        UNRO = (16. + 4.3 * E)/(24.8 - 2.5 * E);
    else if (E <= 10)
        UNRO=(50. + 3. * E)/(129. / E - 4.0);
    else
        UNRO=(44. + 19. * E)/(256. / E + 0.4);

    N1 = *CONLAM * TVALUE * (UNRO + (*CONNU-1.)/2.);
    N1++;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        compute "N0" and "N2"
        for complex singularities               (81) and (89)
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
    if (SID != 0.0)
    {
        DD = DECDIG + 2*(min(MD-1, 1)) + MD/4;
        if (CASE1 && (MD != 0))
        {
            P = SRD/(*CONLAM);
            Q = SID/(*CONLAM);
            U = COM_TalbotINV(P,Q,TETD,pi);
            if (U > 0.0)
            {
                N0=(2.3*DD+SRD*TVALUE)/U;
                N0++;
            }
        }
        GAMM = (*CONSIG - sigma0)/(*CONLAM);
        Y = V/1000.;
        ETA = min(1.78, 1.236+0.0064*pow(1.78,DD));
        ETA = ETA*(1.09-Y*(0.92-0.8*Y));
        N2  = ETA*(*CONNU)*(2.3*DD+OMEGA)/(3.+4.*GAMM+1./exp(GAMM));
        N2++;
    }
    else
    {
        /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            compute "N0" and "N2"
            for real singularities              (82) and (89)
        >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/
        DD = DECDIG+2*(min(MULT[KMAX]-1, 1)) + MULT[KMAX]/4;
        N2 = (2.3*DD + OMEGA)/2.;
        N2++;
        ICM1 = IC-1;
        if (DD >= ICM1)
        {
            Q = 0.0;
            for (L=0; L<Nsings; L++)
                if (MULT[L] != 0)
                {
                    DD  = DECDIG+2*(min(MULT[L]-1, 1)) + MULT[L]/4;
                    SRD = creal(SINGS[L]) - sigma0;
                    if ( DD >= ICM1 && SRD < 0.0)
                    {
                        P = SRD/(*CONLAM);
                        U = COM_TalbotINV(P,Q,TETD,pi);
                        if (U > 0.0)
                        {
                            N0L = (2.3*DD + TVALUE*SRD)/U;
                            if (N0L > N0) N0=N0L;
                        }
                    }
                }
            N0++;
        }
    }
    *NOPTS = max(N0, N1);
    *NOPTS = max(N2, *NOPTS);
}



double COM_TalbotINV (double P, double Q, double TETA, double pi)
/**************************************************************************

    SHARED UTILITY PACKAGE: COM_TalbotINV FUNCTION (internal utility)


    PURPOSE
    =======
    Given a complex number s = P+i*Q, the equation
                s = z / (1 - cexp(-z))
    is solved with respect to z applying Newton's method for real roots.
    The returned value is the opposite of the real part of z.

    CALLING SEQUENCE
    ================
           U = COM_TalbotINV (P, Q, TETA, pi);

    where U is an approximation to -real(z), where z is the
                    principal inverse of s.


    INPUT PARAMETERS
    ================
    P, Q, TETA    - (float) real part, imaginary part and argument of s
                    respectively. It is supposed that
                                  P <= 0.0
                                  Q >= 0.0
                                  pi/2 <= TETA <= pi

    pi            - (double) the value of pi, already computed as 4*atan(1)
                    in the calling function.

   OUTPUT PARAMETERS
   =================
   -

 **************************************************************************/
{
    double ANG, ARG, B, C, E, G, R, S, X, Y, DY, U;
    int K;

    R=P*P+Q*Q;
    R=sqrt(R);

    /* set the initial guess */
    X = 13.0/(5.0-2*P-Q-0.45*exp(P));
    Y = 2*pi-X;
    K = 1;

    /* Newton's method iterations */
    do
    {
        ANG = Y - TETA;
        S = sin(ANG);
        C = cos(ANG);
        B = Q-Y;
        E = -B/S;
        ARG = E/R;
        if (ARG <= 0.0)
        {
            U = -1.0;
            break;
        }
        U = log(ARG);
        G = P+B*C/S+U;
        DY = B*G/(1.0+E*(E-2*C));

        K = K+1;
        Y = Y+DY;
    } while ( fabs(DY) >= 0.0001 && K <= 100);
    /* Test for convergence (tolerance is 1.e-4) and
            for number of iterations (maximum is 100) */

    return U;
}

