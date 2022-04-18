extern void oqs_quartic_solver(double coeff[5], complex double roots[4]);      
extern void oqs_quartic_solver_dl(long double coeff[5], complex long double roots[4]); 
extern void csolve_quartic_abramovitz_cmplx(double *coeff, complex double sol[4]);
extern void fast_quartic_solver(double coeff[5], complex double solqua[4]);
extern void CquarticRoots (double cc[5], int *nReal, complex double root[4]);
extern void solve_numrec (double coeff[5], complex double csol[4], int *ok);
extern void solve_numrecl(long double *coeff, complex long double *csol, int m, int *ok);
extern void CLDLT_quartic(double coeff[5], complex double roots[4]);      
extern void csolve_quartic_shmakov(double *coeff, complex double sol[4]);
