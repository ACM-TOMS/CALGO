/************************    SEQ_Talbot_pack_DE.h    *************************
 *                                                                           *
 *                   HEADER FILE FOR SEQ_Talbot_pack_DE.c                    *
 *                                                                           *
 *****************************************************************************/


/* >>>   TALBOT'S METHOD FOR DIFFERENTIAL PROBLEMS - SEQUENTIAL PACKAGE PROTOTYPES   <<< */


/* Modified Talbot's method - user level function */
int SEQ_Talbot1_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol),
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax);

/* Classical Talbot's method - user level function */
int SEQ_Talbot2_DE (double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol),
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol,
                    double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[]);

/* Modified Talbot's method summation function - skill level function */
int SEQ_TalbotSUM1_DE (double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval,
                       double complex FF[], unsigned int NTval, double Tval[], double NUMft[], int IFAIL[]);

/* Classical Talbot's method summation function - skill level function */
int SEQ_TalbotSUM2_DE (double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double TVALUE, unsigned int jT, double NUMft[], int IFAIL[]);


