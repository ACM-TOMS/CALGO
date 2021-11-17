/************************    COM_Talbot_pack.h    *************************
 *                                                                        *
 *                              TALBOT SUITE                              *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *   SEQUENTIAL AND PARALLEL NUMERICAL INVERSION OF LAPLACE TRANSFORMS    *
 *                                                                        *
 *                                                                        *
 *                              HEADER FILE                               *
 *                    FOR THE SHARED UTILITY FUNCTIONS                    *
 *                             OF TALBOT SUITE                            *
 *                                                                        *
 *                                                                        *
 *  >>>>>>>>>>>>>        VERSION 1.0    Sept 13th, 2012       <<<<<<<<<<  *
 *                                                                        *
 *                                                                        *
 *                                                                        *
 *           AUTHORS: Laura Antonelli (1), Stefania Corsaro (2-1),        *
 *                    Zelda Marino (2), Mariarosaria Rizzardi (3)         *
 *                                                                        *
 *                   (1) ICAR Naples Branch -                             *
 *                                     National Research Council of Italy *
 *                   (2) DSMRE - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 *                   (3) DSA   - "Parthenope" University, Naples (Italy)  *
 *                                                                        *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * Antonelli L., Corsaro S.,                                              *
 * Marino Z., Rizzardi M. - "Talbot Suite: a parallel software collection *
 *                           for the numerical inversion of Laplace       *
 *                           Transforms".                                 *
 *                           ACM Trans. Math. Softw., vol. ##,            *
 *                           no. #, month year, pp. ##-##.                *
 *                                                                        *
 **************************************************************************/


/* >>>   MACROS   <<< */
#ifndef max
    #define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif


#ifndef min
    #define min(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif


/* >>>   TALBOT'S METHOD SHARED MODULES PROTOTYPES   <<< */

/* Talbot's parameters */
void COM_TalbotPAR (double sigma0, double TVALUE, double tol,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                    double *CONLAM, double *CONSIG, double *CONNU, unsigned int *NOPTS);



/* Opposite of the real part of the principal inverse Z
   of  S = Z / (1 - cexp(-Z)) where S=P+I*Q */
double COM_TalbotINV (double P, double Q, double TETA, double PI);
