/**********************    COM_Talbot_pack_DE.c    ************************
 *                                                                        *
 *                           TALBOT SUITE DE                              *
 *                                                                        *
 *     APPLICATION OF TALBOT'S METHODS TO SOLVE DIFFERENTIAL PROBLEMS     *
 *                                                                        *
 *           THE LAPLACE TRANSFORM IS GIVEN BY NUMERICAL SAMPLES          *
 *                             AND NOT AS A FUNCTION                      *
 *                                                                        *
 *                                                                        *
 *                            COMMON SHARED FUNCTION                      *
 *                                                                        *
 *                            path: code/SRC/COM_DE/                      *
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


#include <math.h>

/* >>>   MACRO   <<< */
#ifndef max
    #define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif


unsigned int COM_TalbotNcorr(double Tmin, double Tmax, double sigma0, double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, double tol)
/**************************************************************************

    SHARED UTILITY PACKAGE: COM_TalbotNcorr FUNCTION


    PURPOSE
    =======
    Compute the correction to the accuracy parameter (NOPTS) for the
    modified Talbot method.

 **************************************************************************/
{
    double H = (Tmax-Tmin)/2;
    double r;
    unsigned int N1max, N2max;

    /* N1max: 1st CORRECTION TO NOPTS */
    r=18; /* or: r=-log(tol); */
    N1max = floor( NOPTS + H*(CONLAM*(CONNU+1)/2+CONSIG/r) + 0.5 );   /*(4.4b)*/


    /* N2max: 2nd CORRECTION TO NOPTS */
    double c, v, y, yp, expm1;
    c = (CONNU-1.0)/2;
    v = (1.0 - (CONSIG-sigma0)/CONLAM)/2;

    /* initialization */
    r  = (1.0-v)/(0.5+c);
    expm1=exp(r)-1.0;
    y  = r*( 1.0/expm1 - c );
    yp = ( exp(r)*(1.0-r) - 1.0 )/expm1/expm1 - c;
    int niter = 0; /* number of iterations */

    /* iterations (for safe < 100) */
    while ( fabs(y-v) > tol  &&  niter < 100 )
    {   niter++;
        r  = r + (v-y)/yp;
        y  = r*( 1.0/expm1 - c );
        yp = ( exp(r)*(1.0-r) - 1.0 )/expm1/expm1 - c;
    }

    /* r is the root for N2max */
    N2max = floor( NOPTS + H*((CONSIG+CONLAM)/r-CONLAM*(CONNU-1)/2) + 0.5 ); /*(4.6b)*/

    /* max between N1max, N2max */
    return max(N1max,N2max);
}

