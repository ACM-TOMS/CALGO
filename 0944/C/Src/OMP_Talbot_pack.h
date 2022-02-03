/************************    OMP_Talbot_pack.h    *************************
 *                                                                        *
 *                              TALBOT SUITE                              *
 *                                                                        *
 *                                   FOR                                  *
 *                                                                        *
 *   SEQUENTIAL AND PARALLEL NUMERICAL INVERSION OF LAPLACE TRANSFORMS    *
 *                                                                        *
 *                                                                        *
 *                              HEADER FILE                               *
 *               FOR THE OMP-BASED PARALLEL IMPLEMENTATION                *
 *                                   OF                                   *
 *                             TALBOT'S METHOD                            *
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


/* >>>   OMP-BASED PARALLEL PACKAGE PROTOTYPES   <<< */

/* Classical Talbot's method */
int OMP_Talbot2 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[]);


/* Modified Talbot's method */
int OMP_Talbot1 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                 double Tmin, double Tmax);


/* Classical Talbot's method summation function [fine-grain parallelism] */
int OMP_TalbotSUM2 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    double TVALUE, double *NUMft);


/* Modified Talbot's method summation function [coarse-grain parallelism] */
int OMP_TalbotSUM1 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    unsigned int NTval, double *TVALUE, double *NUMft, int *IFAIL);
