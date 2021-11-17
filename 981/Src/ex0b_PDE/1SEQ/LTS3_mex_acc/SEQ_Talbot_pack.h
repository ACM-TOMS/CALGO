/************************    SEQ_Talbot_pack.h    *************************
 *                                                                        *
 *                   HEADER FILE FOR SEQ_Talbot_pack.c                    *
 *                                                                        *
 **************************************************************************/


/* >>>   TALBOT'S METHOD SEQUENTIAL PACKAGE PROTOTYPES   <<< */

/* Classical Talbot's method (user-level function) */
int SEQ_Talbot2 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[]);


/* Modified Talbot's method (user-level function) */
int SEQ_Talbot1 (double complex (*LTpt)(double complex s), double sigma0,
                 unsigned int NTval, double Tval[], double tol, double *NUMft, int *IFAIL,
                 unsigned int Nsings, double complex SINGS[], unsigned int MULT[],
                 double Tmin, double Tmax);


/* Classical Talbot's method summation function (skill-level) */
int SEQ_TalbotSUM2 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    double TVALUE, double *NUMft);


/* Modified Talbot's method summation function (skill-level) */
int SEQ_TalbotSUM1 (double complex (*LTpt)(double complex s),
                    double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS,
                    unsigned int NTval, double *Tval, double *NUMft, int *IFAIL);
