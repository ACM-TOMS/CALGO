/*	==================================================================
	RELIADIFF Laplace Inverses database
	==================================================================
 
	Here 89 Laplace Inverse Transforms are defined, to compare results
	obteined by RELIADIFF software on the 89 Laplace Transforms 
	contained in the file dbLaplace.c
	
	Each Inverse function is of the kind
	double gzXX(double)
	where XX is the number of the function in the database.

	===========================================================================
	ARGUMENTS
	===========================================================================
	
	each function requires in input: 
	t: 		double precision: the evaluation point of the Inverse function
	
	===========================================================================
	RETURN VALUE
	===========================================================================
	
	each function returns:
	double precision: the evaluation of the Inverse function in t
	
	==================================================================
	AUTHORS
	==================================================================

		Luisa D'Amore - University of Naples, Federico II
		Rosanna Campagna - University of Naples, Federico II
		Valeria Mele - University of Naples, Federico II
		Almerico Murli - CMCC and SPACI

	==================================================================	
 
*/

#include "dbL.h"

/*==================================================================*/
/* G(t) = (sqrt(pi)*((t/(b-a))^(n-(1/2)))*Besseli(n-(1/2),(1/2)*(b-a)*t))/( gamma(n)*exp((1/2)*(b+a)*t) with a=3/5, b=5/7, n=5 */
double gz1(double t){
	if(t==0) return 1;
	double a1=3./5.;
	double b1=5./7.;
	double bma=b1-a1;
	double bpa=b1+a1;
	int n=5;
	double alpha=1./2.;
	double nu=n-alpha;
	double x=(1./2.)*(bma)*t;
	double outbess=gsl_sf_bessel_Inu (nu, x);
	return (sqrt(M_PI)*(pow(t/bma,nu))*outbess)/(gsl_sf_gamma ((double)n)*exp((1./2.)*(bpa)*t));
}

/*==================================================================*/
/* G(t) = (sqrt(pi)*(t/(2a))^(n-1/2) besselj(n-1/2,a*t))/gamma(n) with a=3/5, n=5 */
double gz2(double t){
	double nd=5.;
	double a=3./5.;
	double alpha=1./2.;
	double arg=a*t;
	double arg2=t/(2*a);
	double pot=nd-alpha;
	double outbess=gsl_sf_bessel_Jnu (pot, arg);
	return (sqrt(M_PI)*(pow(arg2,pot))*outbess)/(gsl_sf_gamma (nd));
}

/*==================================================================*/
/* G(t) = (sqrt(pi)*((t/(2a))^(n-(1/2)))*BesselI(n-(1/2),a*t))/gamma(n) with a=3/5, n=5 */
double gz3(double t){
	double nd=5.;
	double a=3./5.;
	double alpha=1./2.;
	double pot=nd-alpha;
	double arg2=t/(2*a);
	double arg=a*t;
	double outbess=gsl_sf_bessel_Inu (pot, arg);
	return (sqrt(M_PI)*(pow(arg2,pot))*outbess)/(gsl_sf_gamma (nd));
}

/*==================================================================*/
/* G(t) =(a^nu)* besselI(nu,a*t) with a=3/5, nu=3 */
double gz4(double t){
	int nu=3;
	double a=3./5.;
	double arg=a*t;
	double outbess=gsl_sf_bessel_In(nu, arg);
	return (gsl_sf_pow_int(a,nu))*outbess;
}

/*==================================================================*/
/* G(t) = t^(n-1/2)/gamma(n+1/2)  with n=5 */
double gz5(double t){
	int n=5;
	double x=n+(1./2.);
	return pow(t,n-(1./2.))/(gsl_sf_gamma (x));
}

/*==================================================================*/
/* G(t) = t^(n-1)/((n-1)!exp(a*t)) with a=3/5, n=5 */
double gz6(double t){
	double a1=3./5.;
	return gsl_sf_pow_int(t,4)/(gsl_sf_fact(4)*exp(a1*t));
}

/*==================================================================*/
/* G(t) = t^(n-1)/(n-1)! with n=5 */
double gz7(double t){
	return gsl_sf_pow_int(t,4)/gsl_sf_fact(4);
}

/*==================================================================*/
/* G(t) =((n*a^n) besselj(n,a*t))/t with a=3/5    */		
double gz8(double t){ 
	if (t==0) return 1;
	int n=5;
	double a=3./5.;
	double arg=a*t;
	double outbess;
	outbess=gsl_sf_bessel_Jn (n, arg);
	return ((n*gsl_sf_pow_int(a,n))*outbess)/t;
}

/*==================================================================*/
/* G(t) = (at-sin(at))/a^3 with a=3/5 */
double gz9(double t){
	double a1=3./5.;
	return (a1*t-sin(a1*t))/gsl_sf_pow_int(a1,3);
}

/*==================================================================*/
/* G(t) = (sin(at)-atcos(at))/(2*a^3) with a=3/5 */
double gz10(double t){
	double a1=3./5.;
	return (sin(a1*t)-a1*t*cos(a1*t))/(2*gsl_sf_pow_int(a1,3));
}

/*==================================================================*/
/* G(t) = (1+a^2*t^2)*sin(at)-atcos(at) with a=3/5 */
double gz11(double t){
	double a1=3./5.;
  	return (1+gsl_sf_pow_int(a1,2)*gsl_sf_pow_int(t,2))*sin(a1*t)-a1*t*cos(a1*t);
}

/*==================================================================*/
/* G(t) =  (1/a)*sinh(a*t)-t with a=0.5 */
double gz12(double t ) {
    return 2*sinh(0.5*t)-t;
}

/*==================================================================*/
/* G(t) = sinh(at)-sin(at)/2*a^3 with a=3/5 */
double gz13(double t){
	double a1=3./5.;
	return (sinh(a1*t)-sin(a1*t))/(2*gsl_sf_pow_int(a1,3));
}

/*==================================================================*/
/* G(t) =(n*besselI(n,(1/2)*(a-b)*t))/(t*exp((1/2)*(b+a)*t)) with a=3/5, b=5/7, n=5 */
double gz14(double t){													
	if(t==0) return 1;	
	double n=5;
	double a=3./5.;
	double b=5./7.;
	double amb=a-b;
	double apb=b+a;
	double alpha=1./2.;
	double x=alpha*(apb)*t; 
	double arg=alpha*(amb)*t;
	double outbess;
	if(arg<0) arg=-arg;
	outbess=gsl_sf_bessel_In (n, arg);
	return (n*outbess)/(t*exp(x));
}

/*==================================================================*/
/* G(t) =  besselI(nu,a*t/2)/((a^nu)*exp(at/2)) with a=3/5, nu=3*/
double gz15(double t){
	double a=3./5.;
	double y=(a*t)/2;	
	int nu=3;
	double outbess=gsl_sf_bessel_In (nu, y);
	return (outbess)/(gsl_sf_pow_int(a,nu)*exp(y));
}

/*======================================================================*/
/* G(t) = (1-cos(at))/a^2 with a=3/5 */
double gz16(double t){
	double a1=3./5.;
	return (1-cos(a1*t))/gsl_sf_pow_int(a1,2);
}

/*==================================================================*/
/* G(t) = (cos(at)-cos(bt))/(b^2-a^2) with a=3/5, b=5/7 */
double gz17(double t){
	double a1=3./5.;
	double b1=5./7.;
	return (cos(a1*t)-cos(b1*t))/(gsl_sf_pow_int(b1,2)-gsl_sf_pow_int(a1,2));
}

/*==================================================================*/
/* G(t) = exp(-at)-exp(a*t/2)*(cos(1/2*sqrt(3)*a*t)-sqrt(3)sin(1/2*sqrt(3)*a*t)) with a=3/5 */
/*double gz18(double t){
	double a1=3./5.;
	return exp(-a1*t)-exp((a1*t)/2)*(cos(1./2.*sqrt(3.)*a1*t)-sqrt(3.)*sin(1./2.*sqrt(3.)*a1*t));
}
*/

/*==================================================================*/
/* G(t) = cosh(at)-cos(at)/2*a^2  with a=3/5 */
double gz18(double t){
	double a1=3./5.;
   	return (cosh(a1*t)-cos(a1*t))/(2*gsl_sf_pow_int(a1,2));
}

/*==================================================================*/  	
/* G(t) =-((b-c)/(exp(a*t))+(c-a)/(exp(b*t))+(a-b)/(exp(c*t)))/((a-b)*(b-c)*(c-a))  with a=3/5, b=5/7, c=-9/7 */	
double gz19(double t){	
	double a=3./5.;
	double b=5./7.;
	double c=-9./7.;
	return -((b-c)/exp(a*t)+(c-a)/exp(b*t)+(a-b)/exp(c*t))/((a-b)*(b-c)*(c-a));
}

/*======================================================================*/
/* G(t) =  t */  
double gz20(double t ) {
	return t; 
}  

/*==================================================================*/
/* G(t) = (exp(-at)-exp(-bt))/(b-a) with a=3/5, b=5/7 */
double gz21(double t){
	double a1=3./5.;
	double b1=5./7.;
	double num=(exp(-a1*t)-exp(-b1*t));
	double den=b1-a1;
	return num/den;
}

/*==================================================================*/
/* G(t) = t*cos(t) */
double gz22(double t ){
	return t*cos(t);
} 

/*==================================================================*/
/* G(t) = sin(t) */
double gz23(double t ) {
      return sin(t);
}
/*==================================================================*/
/* G(t) = 2/sqrt(3)*[exp(-x/2)*sin(x*sqrt(3)/2)] */  
double gz24(double t ) {
	return 2/sqrt(3.)*(exp(-t/2)*sin(t*sqrt(3.)/2)); 
} 

/*==================================================================*/
/* G(t) = (sin(at)+atcos(at))/(2*a) with a=3/5 */
double gz25(double t){
	double a1=3./5.;
	return (sin(a1*t)+a1*t*cos(a1*t))/(2*a1);
}

/*==================================================================*/ 
/* G(t) = sinh(at)/a  with a=3/5 */
double gz26(double t){
	double a1=3./5.;
	return sinh(a1*t)/(a1);
}

/*==================================================================*/
/* G(t) = t*sin(at)/(2a)  with a=3/5 */
double gz27(double t){
    double a1=3./5.;
	return (t*sin(a1*t))/(2.*a1);
}

/*==================================================================*/
/* G(t) = t*cos(at)  with a=3/5 */
double gz28(double t){
    double a1=3./5.;
	return t*cos(a1*t);
}

/*==================================================================*/
/* G(t) =  sin(bt)/(b*exp(at)) with a=3/5, b=5/7 */
double gz29(double t){
    double a1=3./5.;
    double b1=5./7.;
	return sin(b1*t)/(b1*exp(a1*t));
}

/*==================================================================*/
/* G(t) =t/exp(a*t) with a=3/5 */
double gz30(double t){
	double a=3./5.;
	return t/exp(a*t);
}

/*==================================================================*/

/* G(t) = exp(a^2 t)*(((b*Erf(a*sqrt(t)))/a)-1) + exp(b^2 t)*Erfc(b*sqrt(t)) with a=3/5, b=5/7 */
double gz31(double t){
    double a1=3./5.;
    double b1=5./7.;
	double a2=gsl_sf_pow_int(a1,2);
	double b2=gsl_sf_pow_int(b1,2);
    double x=exp(a2*t);
    double y=exp(b2*t);
    double w=a1*sqrt(t);
    double wb=b1*sqrt(t);
	double erfw=gsl_sf_erf (w);
	double erfcwb=gsl_sf_erfc (wb);
	return x*(((b1*(erfw))/a1)-1.)+y*erfcwb;
}

/*==================================================================*/
/* G(t) = (1./(exp((1/2)*(b+a)*t) )*(t*(besselI(0,(1/2)*(a-b)*t)) + besselI(1,(1/2)*(a-b)*t))  with a=3/5, b=5/7 */
double gz32(double t){
	double a=3./5.;
	double b=5./7.;
	double amb=a-b;
	double apb=b+a;
	double alpha=1./2.;
	double x=alpha*(apb)*t;
	double arg=alpha*(amb)*t;
	return (1./exp(x))*(t*(gsl_sf_bessel_I0 (arg)+gsl_sf_bessel_I1 (arg)));
}
/*==================================================================*/
/* G(t) = exp(-0.5*t)+t+exp(-0.2*t)*sin(t)*/  
double gz33(double t ) {
	return exp(-0.5*t)+t+exp(-0.2*t)*sin(t); 
}  

/*======================================================================*/
/* G(t) = 2*sqrt(t)/sqrt(pi) */
double gz34(double t){
	return 2.*sqrt(t)/sqrt(M_PI);
}

/*==================================================================*/
/* G(t) = exp(a^2 t)*(b-a*Erf(a*sqrt(t))-b*exp(b^2 t)*Erfc(b*sqrt(t)) with a=3/5, b=5/7 */
double gz35(double t){
	double a1=3./5.;
	double b1=5./7.;
	double xa=a1*sqrt(t);
	double xb=b1*sqrt(t);
	return exp(gsl_sf_pow_int(a1,2)*t)*(b1-a1*gsl_sf_erf (xa))-b1*exp(gsl_sf_pow_int(b1,2)*t)*gsl_sf_erfc (xb);
   }

/*==================================================================*/
/* G(t) = (exp(a^2 t)*Erf(a*sqrt(t)))/a with a=3/5 */
double gz36(double t){
	double a1=3./5.;
	double x=a1*sqrt(t);
	double y=gsl_sf_pow_int(a1,2)*t;
	return (exp(y)*gsl_sf_erf (x))/a1;
}

/*==================================================================*/
/* G(t) = Erf(sqrt(a-b)*sqrt(t))/sqrt(a-b)*exp(b*t) with a=3/5, b=5/7 */ 
/* G(t) = 1/2*sqrt(35)*exp(-5*t/7)*Erfi(2*sqrt(t)/sqrt(35))	 with a=3/5, b=5/7 */
/* G(t) = exp(-b*t)*Erfi(sqrt(t*(b-a)))/sqrt(b-a) with a=3/5, b=5/7 */	
double gz37(double t){
	double a=3./5.;
	double b=5./7.;
	double w=sqrt((b-a)*t);
	double dawson_w=gsl_sf_dawson(w);
	double erfi_w=exp(gsl_sf_pow_int(w,2))*2*dawson_w/sqrt(M_PI);
	return exp(-b*t)*erfi_w/sqrt(b-a);
}

/*==================================================================*/
/* G(t) = (n!Hermite(2n+1,sqrt(t)))/((2n+1)!*sqrt(M_PI)) with n=5 */
double gz38(double t){
	double f=gsl_sf_fact(5);
	double H=Hermite(11,sqrt(t));
	double num=f*H;
	double f1=gsl_sf_fact(11);
	double den=f1*sqrt(M_PI);
	return num/den;
}
		
/*======================================================================*/
/* G(t) = 1 */ 
double gz39(double t){
	return 1;
}

/*==================================================================*/
/* G(t) =  0.5*exp(-0.5*t) */ 
double gz40(double t ) { 
	return 0.5*exp(-0.5*t); 
}

/*==================================================================*/
/* G(t) =  (1+t)*exp(t) */
double gz41(double t ) { 
	return (1+t)*exp(t); 
}

/*==================================================================*/
/* G(t) =  exp(2t) */ 
double gz42(double t) {
	return exp(2*t); 
}

/*==================================================================*/
/* G(t) = ((a/exp(at))-(b/exp(bt)))/(a-b) with a=3/5, b=5/7 */
double gz43(double t){
	double a1=3./5.;
	double b1=5./7.;
	return ((a1/exp(a1*t))-(b1/exp(b1*t)))/(a1-b1);
}

/*==================================================================*/
/* G(t) = cos(at) with a=3/5 */
double gz44(double t){
	double a1=3./5.;
	return cos(a1*t);
}

/*==================================================================*/
//21
/* G(t) =  cos(bt)/exp(at) with a=3/5, b=5/7 */
double gz45(double t){
	double a1=3./5.;
	double b1=5./7.;
	return cos(b1*t)*exp(-a1*t);
}

/*==================================================================*/	
/* G(t) = cosh(at) with a=3/5 */
double gz46(double t){
	double a1=3./5.;
	return cosh(a1*t);
}

/*==================================================================*/
/* G(t) = BesselJ(0,a*t)  with a=3/5 */
double gz47(double t){
	double a1=3./5.;
	double a1t=a1*t;
	return gsl_sf_bessel_J0 (a1t); 
}

/*==================================================================*/
/* G(t) = sin(k*t)/t with k=9/11*/
double gz48(double t){
	if (t==0) return 1;
	double k=9./11.;
	return sin(k*t)/t;
}

/*==================================================================*/
/* G(t) = (exp(a^2 t)*Erfc(a*sqrt(t))) with a=3/5 */
double gz49(double t){
	double a1=3./5.;
	double x=a1*sqrt(t);
	double y=gsl_sf_pow_int(a1,2)*t;
	return (exp(y)*gsl_sf_erfc (x));
}

/*==================================================================*/
/* G(t) =  besselI(1,a*t)/(t*exp(at)) with a=3/5 */
double gz50(double t){
	if (t==0) return 1;
	double a=3./5.;
	double x=a*t;
	if (t==0) return 1;
	return gsl_sf_bessel_I1 (x)/(t*exp(x));
}

/*==================================================================*/
/* G(t) = Laguerre(n,t) with n=5 */
double gz51(double t){
	int n=5;
	return Laguerre(n,t);
}
		
/*======================================================================*/
/* G(t) = t^(k-1)/gamma(k) with k=9/11*/
double gz52(double t){
	if (t==0) return 1;
    double k=9./11.;
	return  pow(t,k-1)/(gsl_sf_gamma (k));

}

/*==================================================================*/
/* G(t) = t^(k-1)/(gamma(k)exp(a*t)) with a=3/5, k=9/11*/
double gz53(double t){ 
	if (t==0) return 1;
	double k=9./11.;
	double a1=3./5.;
	return pow(t,k-1)/((gsl_sf_gamma (k))*exp(a1*t));
}

/*======================================================================*/
/* G(t) = 1/sqrt(M_PI*t) */
double gz54(double t){
	if (t==0) return 1;
	return 1./(sqrt(M_PI*t));
}

/*==================================================================*/
/* G(t) = (1-2at)/(sqrt(t)*sqrt(pi)*exp(at))  with a=3/5 */
double gz55(double t){
	if (t==0) return 1;
	double a1=3./5.;
	return (1.-2.*a1*t)/(sqrt(t*M_PI)*exp(a1*t));
}


/*==================================================================*/
/* G(t) = (exp(-bt)-exp(-at))/(2(sqrt(pi)*t^(3/2))) with a=3/5, b=5/7 */
double gz56(double t){
	if (t==0) return 1;
    double a1=3./5.;
    double b1=5./7.;
	return (exp(-b1*t)-exp(-a1*t))/(2.*(sqrt(M_PI)*sqrt(gsl_sf_pow_int(t,3))));	
}

/*==================================================================*/
/* G(t) = (1/(sqrt(t)*sqrt(pi)))-a*exp(a^2 t)*Erfc(asqrt(t))  with a=3/5 */
double gz57(double t){
	if (t==0) return 1;
	double a1=3/5;
	double x=a1*sqrt(t);
	return (1./sqrt(t*M_PI))-a1*exp(gsl_sf_pow_int(a1,2)*t)*(gsl_sf_erfc (x));	
}

/*==================================================================*/
/* G(t) = (n!Hermite(2n,sqrt(t)))/((2n)!*sqrt(M_PI)*sqrt(t)) with n=5 */
double gz58(double t){
	if (t==0) return 1;
	int n=5;
	double f=gsl_sf_fact(n);
	double H=Hermite(2*n,sqrt(t));
	double num=f*H;
	double f1=gsl_sf_fact(2*n);
	double den=f1*sqrt(M_PI*t);
	return num/den;
}
			
/*==================================================================*/		
/*G(t) = G(t) =  (1/(sqrt(PI*t)))*exp(-a*t)*/
double gz59(double t) {
 	double a=1;
    double g;
	g=1./(sqrt(M_PI*t));
	return g*exp(-a*t);
}


/*======================================================================*/
/* G(t) = 1-ExpIntegralE(1,t)*/
double gz60(double t ){   
	return 1.-gsl_sf_expint_E1(t);
}	

/*==================================================================*/
/* G(t) = -EulerGamma-log(t) */
double gz61(double t){
	if (t==0) return 1;
	return -M_EULER-log(t);
}
		
/*==================================================================*/		
/* G(t) = ExpIntegralE(1,t), k=1*/ 
double gz62(double t ) { 
	return gsl_sf_expint_E1(t);
}

/*==================================================================*/
/* G(t) = (exp(-b*t)-exp(-a*t))/t with a=3/5, b=5/7 */
double gz63(double t){
	if(t==0) return 1;
	double a1=3./5.;
	double b1=5./7.;	
	return (exp(-b1*t)-exp(-a1*t))/t;
}

/*==================================================================*/
/* G(t) = (2*(1-cos(a*t)))/t  with a=3/5 */
double gz64(double t){
	if(t==0) return 1;
	double a1=3./5.;
	return (2*(1-cos(a1*t)))/t;
}

/*==================================================================*/
/* G(t) = (2*(1-cosh(a*t)))/t   with a=3/5 */
double gz65(double t){
	if(t==0) return 1;
	double a1=3./5.;	
	return (2*(1-cosh(a1*t)))/t;
}

/*======================================================================*/
/* G(t) = a*(BesselI(1,a*t)+BesselI(0,a*t))/exp(a*t)  with a=3/5 */
double gz66(double t){
      double a1=3./5.;
      double a1t=a1*t;
      return a1*(gsl_sf_bessel_I1 (a1t)+gsl_sf_bessel_I0 (a1t))/exp(a1t);
}

/*======================================================================*/
/* G(t) = ((t/k)^((mu-1)/2))*BesselI(mu-1,2*sqrt(kt)) with k=9/11, mu=4 */
double gz67(double t){
	int mu=4;
	double k=9./11.;
	double esp=(mu-1)/2.;
	double y=t*k;
	double x=t/k;
	double arg=2.*sqrt(y);
	double outbess=gsl_sf_bessel_In(mu-1, arg);
	return (pow(x,esp)*outbess);
}

/*==================================================================*/
/* G(t) = sinh(2*sqrt(k*t))/(sqrt(pi)*sqrt(k)) with k=9/11*/
double gz68(double t){
	double k=9./11.;
	double a=2.;
	double x=a*sqrt(k*t);
	return  sinh(x)/sqrt(M_PI*k);
}

/*==================================================================*/
/* G(t) = Cosh(2*sqrt(k*t))/(sqrt(pi)*sqrt(t)) with k=9/11*/
double gz69(double t){
	if (t==0) return 1;
	double k=9./11.;
	double a=2.;
	double x=a*sqrt(k*t);
	return cosh(x)/sqrt(M_PI*t);		
}

/*======================================================================*/
/* G(t) = (k*exp(-k^2/(4.*t)))/(2*(sqrt(pi)*t^(3/2))) with k=9/11*/
double gz70(double t){
	if (t==0) return 1;
	double k=9./11.;
	double tmp=2.;
	double tmp2=4.;
	double x;
	x=-gsl_sf_pow_int(k,2)/(tmp2*t);
	return (k*exp(x))/(tmp*(sqrt(M_PI)*pow(t,(3./2.))));
}

/*==================================================================*/
/* G(t) = integral_0_t( (1./sqrt(M_PI*u))* exp(-k^2/(4*u))* ((t-u)^(n-1))/(n-1)! du ) , n=5*/     
double f1(double u, void * params) {
	double t= *(double *) params;
	double k=9./11.;
	double a1=1./sqrt(M_PI*u);
	double a2=exp(-gsl_sf_pow_int(k,2)/(4*u));
	double a3=gsl_sf_pow_int(t-u,4)/gsl_sf_fact(4);
	return a1*a2*a3;
}

double gz71(double t){
	if(t==0) return 1;
	double result, error;
	double epsabs=0.0;
	double epsrel=1.0E-04;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(8000);    
	gsl_function F;
    F.function = &f1;
    F.params = &t;
	gsl_integration_qags(&F, 0, t, epsabs, epsrel, 8000, w, &result, &error);
	return result;
}
		

/*==================================================================*/
/* G(t) =2*sqrt(t)*exp(-k^2/(4.*t))/sqrt(pi)-k*Erfc(k/2*sqrt(t)) with k=9/11*/
double gz72(double t){
	double k=9./11.;
	double tmp=2.;
	double tmp2=4.;
	double x=-gsl_sf_pow_int(k,2)/(tmp2*t);
	double y=k/(tmp*sqrt(t));
	double w=tmp*sqrt(t)/sqrt(M_PI);
	return w*exp(x)-k*gsl_sf_erfc (y);
} 

/*==================================================================*/
/* G(t) =-exp(a*k)*exp(a^2*t)*Erfc(a*sqrt(t)+k/(2*sqrt(t)))+Erfc(k/(2*sqrt(t))) with a=3/5, k=9/11*/
double gz73(double t){
	double k=9./11.;
	double a=3./5.;
	double y=k/(2.*sqrt(t));
	double ay=a*sqrt(t)+y;
	return -exp(a*k)*exp(gsl_sf_pow_int(a,2)*t)*gsl_sf_erfc (ay)+gsl_sf_erfc (y);
}

/*==================================================================*/
/* G(t) = Erfc(k/(2*sqrt(t))) with k=9/11*/
double gz74(double t){
	if (t==0) return 1;
	double k=9./11.;
	double x;
	x=k/(2.*(sqrt(t)));
	return gsl_sf_erfc (x);
}

/*==================================================================*/
/* G(t) = exp(-(1+(t/2)))*Besseli(0,2*sqrt(t))*/ 
double gz75(double t){
    double arg=2.*sqrt(t);
	return exp(-(1.+(t/2.)))*gsl_sf_bessel_I0 (arg);
}

/*==================================================================*/
/* G(t) = exp(a*k)*exp(a^2*t)*Erfc(a*sqrt(t)+k/(2*sqrt(t))) with a=3/5, k=9/11 */
double gz76(double t){
        double k=9./11.;
        double a=3./5.;
        double y=k/(2.*sqrt(t));
        double ay=a*sqrt(t)+y;
	return exp(a*k)*exp(gsl_sf_pow_int(a,2)*t)*gsl_sf_erfc (ay);
}

/*==================================================================*/
/* G(t) = BesselJ(0,a*sqrt(t^2+2kt)) with a=3/5, k=9/11 */
double gz77(double t){
        double k=9./11.;
        double a=3./5.;
        double x=2.*k*t;
        double y=gsl_sf_pow_int(t,2);
		double w=a*sqrt(y+x);
        return gsl_sf_bessel_J0 (w); 
}
		
/*==================================================================*/	
/* G(t) =(exp(-k^2/(4.*t)))/(sqrt(pi)*sqrt(t)) with k=9/11*/
double gz78(double t){
	if (t==0) return 1;
	double k=9./11.;
	double tmp2=4.;
	double x=-gsl_sf_pow_int(k,2)/(tmp2*t);
	return exp(x)/sqrt(M_PI*t);
} 

/*==================================================================*/
/* G(t) =exp(-k^2/(4.*t))/(sqrt(pi)*sqrt(t))-a*exp(a*k)*exp(a^2*t)*Erfc(a*sqrt(t)+k/(2*sqrt(t))) with a=3/5, k=9/11*/
double gz79(double t){
	if(t==0) return 1;
	double k=9./11.;
	double a=3./5.;
	double x=-gsl_sf_pow_int(k,2)/(4.*t);
	double y=k/(2.*sqrt(t));
	double w=sqrt(t)*sqrt(M_PI);
	double ay=a*sqrt(t)+y;
	return (exp(x)/w)-a*exp(a*k)*exp(gsl_sf_pow_int(a,2)*t)*gsl_sf_erfc (ay);
}

/*==================================================================*/
/* G(t) =  1/sqrt(M_PI*t)*exp(-(1+t))*cosh(2*sqrt(t))*/
double gz80(double t ) { 
	return 1./sqrt(M_PI*t)*exp(-(1.+t))*cosh(2.*sqrt(t)); 
}

/*==================================================================*/
/* G(t) = ((t/k)^((mu-1)/2))*BesselJ(mu-1,2*sqrt(kt)) with k=9/11, mu=4 */
double gz81(double t){
	int mu=4;
	double k=9./11.;
	double esp=(mu-1)/2.;
	double y=t*k;
	double arg=2.*sqrt(y);
	double x=t/k;
	double outbess=gsl_sf_bessel_Jn(mu-1, arg);
	return pow(x,esp)*outbess;
}

/*==================================================================*/
/* G(t) = sin(2*sqrt(k*t))/(sqrt(pi)*sqrt(k)) with k=9/11*/
double gz82(double t){
	double k=9./11.;
	double a=2.;
	double x=a*sqrt(k*t);
	return sin(x)/sqrt(M_PI*k);
}

/*==================================================================*/
/* G(t) = Besselj0(2*sqrt(k*t)) with k=9/11*/
double gz83(double t){
    double k=9./11.;
	double x=2.*sqrt(k*t);
	return gsl_sf_bessel_J0 (x); 
}

/*==================================================================*/
/* G(t) = Cos(2*sqrt(k*t))/(sqrt(pi)*sqrt(t)) with k=9/11*/
double gz84(double t){
	if (t==0) return 1;
	double k=9./11.;
	double a=2.;
	double x=a*sqrt(k*t);
	return cos(x)/sqrt(M_PI*t);		
}


/*======================================================================*/
/*
NEXT 5 are the Inverse Laplace Transform of
function TEST 1, 2, 3, 4, 5 of the paper
"On the Numerical Inversion of Laplace Transforms: Comparison of Three
New Methods on Characteristic Problems from Applications".
DEAN G. DUFFY
*/
/*======================================================================*/

double gz85(double t){
	double PI2=gsl_sf_pow_int(M_PI,2);
	double A=(0.5-(exp(2.0)/(exp(2.0)-1.0)))*exp(-t) + 0.5;
	double n=0;
	double sum=0;
	double B;
	double FEX;
	do{
		n=n+1.0;
		B=sin((n*M_PI*t)-atan(n*M_PI));
		B=B/(n*(sqrt((gsl_sf_pow_int(n,2)*PI2)+1.0)));
		sum=sum+B;
	}
	while (n<=1.0E+7); 
	sum=sum/M_PI;
    FEX=A-sum;
    return FEX;
}

/*==================================================================*/

double ff(double x, void * params) {
	return x*tan(x)-1;
}

double bcalc(int i){
	int status;
    int iter = 0, max_iter = 1000000;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    double r = 0; //radix   	
	double x_lo, x_hi; //region
	gsl_function F;     
	double p=0.;

    F.function = &ff;
    F.params = &p;
	
	if(i==0){
		x_lo = 0;
		x_hi = M_PI/3.; //region
	}	
	else{
		x_lo = i*M_PI;
		x_hi = (i+1)*M_PI+1./2.; //region
	}
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);     
    do{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
	}while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);     
    return r;
}

double gz86(double t) { 
	double b,b2,a1,a21,a22;
	double sum=0.;
	double add=0.;
	int i=0;

	while (fabs(add)>=DBL_EPSILON*fabs(sum)){		
		sum=sum+add;
		b=bcalc(i);
		b2=gsl_pow_int(b,2);
		a1=100.+1./b2;      
		a21=2.*b*sin(b/2.)*exp(-b2*t); 
		a22=(2.+b2)*cos(b);		
		add=a1*a21/a22;
		i++;
	}
	return -1./2.+sum; 
}


/*==================================================================*/

double FST3(double u, void * params) {
	double m,theta, arg1, arg2, r, c;
	double t= *(double *) params;
	double FST3val;
	r=0.5;
	c=0.4;
	m=(1.+gsl_sf_pow_int(u,2))/(1.+(gsl_sf_pow_int(c,2)*gsl_sf_pow_int(u,2)));
	m=pow(m,0.25);
	theta=atan(u)-atan(u*c);
	theta=theta/2;
	arg1=-r*m*sqrt(u/2)*(cos(theta)-sin(theta));
	arg2=t*u-(r*m*sqrt(u/2)*(cos(theta)+sin(theta)));
	FST3val=exp(arg1)*sin(arg2)/u;
	return FST3val;
}

double gz87(double t) { 
	double a=0;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(8000);
    double result, error;
    double alpha = t;
	double epsabs=0.0;
	double epsrel=1.0E-4;
	double FEX;
	gsl_function F;
    F.function = &FST3;
    F.params = &alpha;
    gsl_integration_qagiu(&F, a, epsabs, epsrel, 8000,w, &result, &error); 
    gsl_integration_workspace_free (w);
	FEX=(result/M_PI)+0.5;
	return FEX; 
}

/*==================================================================*/

double FST4(double u, void * params) { 
	double FST4val;
	double t= *(double *) params;
	double k;
	double u1=4.*(2.0-sqrt(3.0));
	double u2=4.*(2.0+sqrt(3.0));      
	double arg=(u1-gsl_sf_pow_int(u,2))*(u2-gsl_sf_pow_int(u,2));
	if (arg>= 0.0){ 
		k=acos(.25*sqrt(arg));
		FST4val=(sin(u*t+2.*k)-sin(u*t-2*k))/u;
	}
	return FST4val;
}

double gz88(double t) { 
	double FEX;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (8000);       
    double alpha = t;
	double epsabs=0.0;
	double epsrel=1.0E-4;
	double u1=2.*sqrt(2.0-sqrt(3.0));
	double u2=2.*sqrt(2.0+sqrt(3.0));
	double a=0.0;
	double b=u1;
	double result, abserr, result1, result2;     
    gsl_function F;
    F.function = &FST4;
    F.params = &alpha;    
    gsl_integration_qag(&F, a, b, epsabs, epsrel, 8000, 6, w, &result, &abserr);
	result1=result;
	a=u2;
	b=4.0;
    gsl_integration_qag(&F, a, b, epsabs, epsrel, 8000, 6, w, &result, &abserr); 
    gsl_integration_workspace_free (w);
	result2=result;
	FEX=((-result1+result2)/M_PI)+1.0;
	return FEX; 
}

/*==================================================================*/ 

  double FST51(double u, void * params) {                  
	double t= *(double *) params;
	double FST51val;	  
	double N=0.5;
	double c=(1.-N)/N;
	double R=sqrt(gsl_sf_pow_int(u,2)+(gsl_sf_pow_int(N,2)*(gsl_sf_pow_int(c,2)-gsl_sf_pow_int(u,2)))); 
	double c2mu2=gsl_sf_pow_int(c,2)-gsl_sf_pow_int(u,2);
	FST51val=cosh(t*u)*(u*(sqrt((R+u)/2))+(sqrt(c2mu2)*sqrt((R-u)/2)))/(R*sqrt(c2mu2)*sqrt(u));   
	return FST51val;
}

double FST52(double u, void * params) {
	double t= *(double *) params;
	double FST52val;
	double N=0.5;
	double c=(1.-N)/N;
	double c2pu2=gsl_sf_pow_int(c,2)+gsl_sf_pow_int(u,2);
	FST52val=cos(t*u)*(u-sqrt(c2pu2))/(sqrt(u)*sqrt(c2pu2)*sqrt(N*sqrt(c2pu2)-u));   
	return FST52val;
}

double gz89(double t) { 
	double FEX;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(8000);     
    double result, error;
	double alpha = t;
	double epsabs=0.0;
	double epsrel=1.0E-4;
	double N=0.5;
	double b=sqrt((1.0-N)/(1.0+N));
	double c=(1.0-N)/N;
	double result1, result2;
    gsl_function F;
    F.function = &FST51;
    F.params = &alpha;
    gsl_integration_qags(&F, 0, c, epsabs, epsrel, 8000,w, &result, &error);
	result1=result;
    gsl_function G;
    G.function = &FST52;
    G.params = &alpha;
    gsl_integration_qag(&G, 0, b, epsabs, epsrel, 8000, 6, w, &result, &error);
	gsl_integration_workspace_free(w);
	result2=result;
	FEX=(result1+result2)*(2.0/M_PI);
	return FEX; 
}












	
	
	