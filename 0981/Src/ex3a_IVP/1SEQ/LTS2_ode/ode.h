/** ode.h **/

/** ODEfun_PT  : global type for a pointer to a function f(x,U,s) such that U'=f(x,U,s) where s is a parameter **/
typedef void (*ODEfun_PT) (double x, double U[], double Up[], double complex S);

/** INITcond_PT: global type for a pointer to a function for the intial conditions of the ODE problem **/
typedef double complex (*INITcond_PT) (double x0, double complex S, double U0[]);



int my_ode( ODEfun_PT f, double complex S, unsigned int NXval, double X[], double complex U[], double relerr, double abserr);


void ode ( ODEfun_PT f, double complex S, int neqn, double y[], double *t, double tout, double relerr, double abserr, int *iflag, double work[], int iwork[] );


void de ( ODEfun_PT f, double complex S, int neqn, double y[], double *t, double tout, double relerr, double abserr, int *iflag, double yy[],
          double wt[], double p[], double yp[], double ypout[], double phi[],
          double alpha[], double beta[], double sig[], double v[], double w[],
          double g[], int *phase1, double psi[], double *x, double *h, double *hold,
          int *start, double *told, double *delsgn, int *ns, int *nornd, int *k, int *kold, int *isnold );


void step ( double *x, double y[], ODEfun_PT f, double complex S, int neqn, double *h, double *eps, double wt[], int *start, double *hold,
            int *k, int *kold, int *crash, double phi[], double p[], double yp[], double psi[], double alpha[], double beta[],
            double sig[], double v[], double w[], double g[], int *phase1, int *ns, int *nornd );


void intrp ( double x, double y[], double xout, double yout[], double ypout[], int neqn, int kold, double phi[], double psi[] );


int i4_sign ( int i );


double r8_abs ( double x );


double r8_add ( double x, double y );


double r8_epsilon ( void );


double r8_max ( double x, double y );


double r8_min ( double x, double y );


double r8_sign ( double x );

