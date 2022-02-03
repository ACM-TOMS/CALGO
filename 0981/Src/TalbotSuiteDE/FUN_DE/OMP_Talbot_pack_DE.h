/**********************    OMP_Talbot_pack_DE.h    ***************************
 *                                                                           *
 *                   HEADER FILE FOR OMP_Talbot_pack_DE.c                    *
 *                                                                           *
 *****************************************************************************/


/** 1a) Modified Talbot's method: user-level function [coarse-grain parallelism] **/
int OMP_Talbot11_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS);

/*+ 1b) Modified Talbot's method: skill-level (summation) function [coarse-grain parallelism] **/
int OMP_TalbotSUM11_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS);



/** 2a) Modified Talbot's method: user-level function [fine-grain parallelism] **/
int OMP_Talbot12_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS);

/** 2b) Modified Talbot's method: skill-level (summation) function [fine-grain parallelism] **/
int OMP_TalbotSUM12_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS);



/** 3a) Modified Talbot's method: user-level function [hybrid parallelism] **/
int OMP_Talbot13_DE(double complex* (*LTsamples)(unsigned int NXval, double Xval[], unsigned int NOPTS, double complex S[], double tol, int THREADS),
                    /*    1st par:  user-defined function for LT samples  */
                    double sigma0, unsigned int NXval, double *Xval, unsigned int NTval, double *Tval, double tol, double *NUMft, int *IFAIL,
                    unsigned int Nsings, double complex SINGS[], unsigned int MULT[], double Tmin, double Tmax, int THREADS1, int THREADS2);

/** 3b) Modified Talbot's method: skill-level (summation) function [hybrid parallelism] **/
int OMP_TalbotSUM13_DE(double CONLAM, double CONSIG, double CONNU, unsigned int NOPTS, unsigned int NXval, double complex FF[],
                       unsigned int NTval, double *Tval, double NUMft[], int IFAIL[], int THREADS1, int THREADS2);
