/*	==================================================================
	RELIADIFF Laplace Transforms database
	==================================================================
 
	Here 89 Laplace Transforms are defined, to test RELIADIFF software 
	and compare its results with the 89 Laplace Inverse Transforms 
	contained in the file dbInvLaplace.c
	
	Each Transform function is of the kind
	T<double> fzXX(T<double> z)
	where XX is the number of the function in the database


	===========================================================================
	ARGUMENTS
	===========================================================================
	
	each function requires in input: 
	z: 		(TADIFF) double precision: the evaluation point of the Transform
	
	===========================================================================
	RETURN VALUE
	===========================================================================
	
	each function returns:
	(TADIFF) double precision: the evaluation of the Transform in z
	
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
          

/*======================================================================*/
/* F(z)= 1./(((z+a)^n)*((z+b)^n)) with a=3/5, b=5/7, n=5 */	
T<double> fz1(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	int n=5;
	T<double> a2=pow_int_T(z+a,n);
	T<double> b2=pow_int_T(z+b,n);
	return 1./(a2*b2);
}

/*==================================================================*/
/* F(z)= 1/((z^2+a^2)^n) with a=3/5, n=5  */	
T<double> fz2(T<double> z) {
	int n=5;
	double a=3./5.;
	T<double> den=pow_int_T(z,2)+gsl_sf_pow_int(a,2);
	return 1./(pow_int_T(den,n));
}

/*==================================================================*/
/* F(z)= 1/((z^2-a^2)^n) with a=3/5, n=5   */
T<double> fz3(T<double> z) {
	int n=5;
	double a=3./5.;
	T<double> den=pow_int_T(z,2)-gsl_sf_pow_int(a,2);
	return 1./(pow_int_T(den,n));
}
   
/*======================================================================*/
/* F(z)= ((z-sqrt(z^2-a^2))^nu)/(sqrt(z^2-a^2)) with a=3/5, nu=3 */	
T<double> fz4(T<double> z) {
	int nu=3;
	double a=3./5.;
	T<double> x=sqrt(pow_int_T(z,2)-gsl_sf_pow_int(a,2));
	T<double> num=pow_int_T(z-x,nu);
	return num/x;
}

/*======================================================================*/
/* F(z)=1/(z^(n+(1/2))) with n=5   */
T<double> fz5(T<double> z) {
	int n=5;
	double pot=n+(1./2.);
	return 1./pow(z,pot);
}

/*======================================================================*/
/* F(z)=1/((z+a)^5)  with a=3/5 */
T<double> fz6(T<double> z) {
	int n=5;      
	return 1./pow_int_T((z+3./5.),n);            
}

/*==================================================================*/
/* F(z)=1/z^5   */
T<double>  fz7(T<double> z) {
    return 1./pow(z,5);
}

/*==================================================================*/
/* F(z)= (sqrt(z^2+a^2)-z)^n, a=3/5 */		
T<double> fz8(T<double> z) {
   int n=5;
   double a=3./5.;
   T<double> num=sqrt(pow_int_T(z,2)+gsl_sf_pow_int(a,2));
   return pow_int_T(num-z,n);
}

/*======================================================================*/
/* F(z)=1/(z^2*(z^2+a^2)  with a=3/5   */
T<double> fz9(T<double> z){
	return 1./(pow_int_T(z,2)*(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2)));  
}

/*==================================================================*/
/* F(z)=1/(z^2+a^2)^2  with a=3/5   */
T<double> fz10(T<double> z) {
	return 1./pow_int_T(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2),2);      
}

/*==================================================================*/
/* F(z)=(8*a^3*z^2)/(z^2+a^2)^3  with a=3/5   */
T<double> fz11(T<double> z) {
	return (8*gsl_sf_pow_int(3./5.,3)*pow_int_T(z,2))/pow_int_T(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2),3);       
}

/*==================================================================*/
/* F(z)=(a^2)/(z^2*(z^2 - a^2))  with a=0.5   */
T<double> fz12(T<double> z) {
	T<double> tmp=pow_int_T(z,2)-0.25;
	return 0.25/(pow_int_T(z,2)*tmp);
}

/*==================================================================*/
/* F(z)=1/(z^4-a^4)  with a=3/5   */
T<double> fz13(T<double> z) {
	return 1./(pow_int_T(z,4)-gsl_sf_pow_int(3./5.,4));                    
}

/*==================================================================*/
/* F(z)= (a-b)^n/((sqrt(z+a)+sqrt(z+b))^2n) with a=3/5, b=5/7, n=5 */	
T<double> fz14(T<double> z) {									
	double a=3./5.;
	double b=5./7.;
	int n=5;
	int pot=2*n;
	T<double> sza=sqrt(z+a);
	T<double> szb=sqrt(z+b);
	return gsl_sf_pow_int(a-b,n)/(pow_int_T(sza+szb,pot));
}

/*==================================================================*/
/* F(z)=1/(((sqrt(z+a)+sqrt(z))^2nu)*(sqrt(z+a)*sqrt(z))) with a=3/5, nu=3 */	
T<double> fz15(T<double> z) {
	int nu=3;
	double a=3./5.;
	T<double> sza=sqrt(z+a);
	T<double> sz=sqrt(z);
	return 1./(pow_int_T(sza+sz,2*nu)*(sza*sz));
}

/*======================================================================*/
/* F(z)=1/(z*(z^2+a^2))  with a=3/5 */
T<double> fz16(T<double> z) {
	return 1./(z*(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2)));      
}

/*==================================================================*/
/* F(z)=z/((z^2+a^2)(z^2+b^2))  with a=3/5, b=5/7 */
T<double> fz17(T<double> z) {
	return z/((pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2))*(pow_int_T(z,2)+gsl_sf_pow_int(5./7.,2)));  
}

/*==================================================================*/
/* F(z)=(3*a^2)/(z^3+a^3)  with a=3/5   */
/*T<double> fz18(T<double> z) {
	return (3*gsl_sf_pow_int(3./5.,2))/(pow_int_T(z,3)+gsl_sf_pow_int(3./5.,3));         
}
*/

/*==================================================================*/
/* F(z)=z/(z^4-a^4)  with a=3/5   */
T<double> fz18(T<double> z) {
	return z/(pow_int_T(z,4)-gsl_sf_pow_int(3./5.,4));                    
}

/*==================================================================*/
/* F(z)=1/((z+a)(z+b)(z+c)) with a=3/5, b=5/7, c=-9/7   */
T<double>  fz19(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	double c=-9./7.;
	return 1./((z+a)*(z+b)*(z+c));
}

/*======================================================================*/
/* F(z)=1/(z^2)  */
T<double> fz20(T<double> z) { 
	return 1/pow_int_T(z,2);
}

/*==================================================================*/
/* F(z)=1/((z+a)(z+b))  with a=3/5, b=5/7 */
T<double> fz21(T<double> z) {
	return 1./((z+3./5.)*(z+5./7.));
}

/*==================================================================*/
/* F(z)=(z^2-1)/(z^2+1)^2 */
T<double> fz22(T<double> z){     
	T<double> tmp,pot;  
	pot=pow_int_T(z,2); 
	tmp=pow_int_T((pot+1),2); 
	return (pot-1)/tmp; 
}

/*==================================================================*/
/* F(z)=1/(z^2+1) */ 
T<double> fz23(T<double> z) { 
	return 1/(pow_int_T(z,2)+1);
} 

/*==================================================================*/
/* F(z)=1/(z^2+z+1) */  
T<double> fz24(T<double> z) { 
	return 1/(pow_int_T(z,2)+z+1); 
}

/*==================================================================*/
/* F(z)=z^2/(z^2+a^2)^2  with a=3/5 */
T<double> fz25(T<double> z) {
	return pow_int_T(z,2)/pow_int_T(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2),2);     
}

/*==================================================================*/
/* F(z)=1/(z^2-a^2)  with a=3/5   */
T<double> fz26(T<double> z) {
	return 1./(pow_int_T(z,2)-gsl_sf_pow_int(3./5.,2));     
}

/*==================================================================*/
/* F(z)= z/((z^2+a^2)^2)  with a=3/5 */
T<double> fz27(T<double> z) {
	return z/pow_int_T((pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2)),2);
}

/*==================================================================*/
/* F(z)= (z^2-a^2)/((z^2+a^2)^2)  with a=3/5   */
T<double> fz28(T<double> z) {
	return (pow_int_T(z,2)-gsl_sf_pow_int(3./5.,2))/pow_int_T((pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2)),2);
}

/*==================================================================*/
/* F(z)=1/((z+a)^2+b^2)  with a=3/5, b=5/7 */
T<double> fz29(T<double> z) {
	return 1./(pow_int_T(z+3./5.,2)+gsl_sf_pow_int(5./7.,2));
}

/*==================================================================*/
/* F(z)=1/(z+a)^2 with a=3/5 */
T<double>  fz30(T<double> z) {
	double a=3./5.;
	return 1./pow(z+a,2);
}

/*==================================================================*/
/* F(z)= (b^2-a^2)/(sqrt(z)*(z-a^2)*(sqrt(z)+b)) with a=3/5, b=5/7   */		
T<double> fz31(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	double a2=gsl_sf_pow_int(a,2);
	double b2=gsl_sf_pow_int(b,2);
	T<double> sz=sqrt(z);
	return (b2-a2)/(sz*(z-a2)*(sz+b));
}

/*==================================================================*/
/* F(z)= 1./(sqrt(z+a)*((z+b)^(3/2))) with a=3/5, b=5/7 */	
T<double> fz32(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	double pot=3./2.;
	T<double> b2=pow(z+b,pot);
	T<double> sza=sqrt(z+a);
	return 1./(sza*b2);
}

/*==================================================================*/
/* F(z)=1/(z+0.5)+1/(z^2)+1/(1+(z+0.2)^2)   */ 
T<double> fz33(T<double> z) { 
	if (z[0]==0) return 1;
	return 1./(z+0.5)+1./pow_int_T(z,2)+1./(1.+pow_int_T((z+0.2),2));
}

/*======================================================================*/
/* F(z)= 1./z^(3/2)  */			
T<double> fz34(T<double> z) {
	return 1./pow(z,3./2.);
}

/*==================================================================*/
/* F(z)= (b^2-a^2)/((z-a^2)*(sqrt(z)+b))  with a=3/5, b=5/7   */
T<double> fz35(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	return (gsl_sf_pow_int(b,2)-gsl_sf_pow_int(a,2))/((z-gsl_sf_pow_int(a,2))*(sqrt(z)+b));
}

/*==================================================================*/
/* F(z)= 1./(sqrt(z)*(z-a^2)) with a=3/5   */
T<double> fz36(T<double> z) {
	double a=3./5.;
	return 1./(sqrt(z)*(z-gsl_sf_pow_int(a,2)));
}

/*==================================================================*/
/* F(z)= 1./((z+b)*(sqrt(z+a))) with a=3/5, b=5/7 */	        
T<double> fz37(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	return 1./((z+b)*sqrt(z+a));
}

/*==================================================================*/
/*F(z)=((1-z)^n)/(z^(n+3/2)) with n=5  */
T<double> fz38(T<double> z){
	int n=5;
	T<double> tmp=pow_int_T(1-z,n);
	return tmp/pow(z,n+3./2.);
}
	
/*======================================================================*/
/* F(z)=1/z */
T<double> fz39(T<double> z) {
	return 1/z;
}

/*==================================================================*/
/* F(z)=1/(1+2*z) */  
T<double> fz40(T<double> z) {   
	return 1/(1+2*z); 
}

/*==================================================================*/
/* F(z)=z/(z-1)^2   */ 
T<double> fz41(T<double> z){    
	return z/pow_int_T((z-1),2); 
}

/*==================================================================*/
/* F(z)=1/(z-2)  */
T<double> fz42(T<double> z) {
	return 1./(z-2);
}

/*==================================================================*/
/* F(z)=z/((z+a)(z+b))  with a=3/5, b=5/7 */
T<double> fz43(T<double> z) {
	return z/((z+3./5.)*(z+5./7.));
}

/*==================================================================*/
/* F(z)=z/(z^2+a^2)  with a=3/5 */
T<double> fz44(T<double> z) {
	return z/(pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2));  
}

/*==================================================================*/
/* F(z)=(z+a)/((z+a)^2+b^2)  with a=3/5, b=5/7 */
T<double> fz45(T<double> z) {
	return (z+3./5.)/(pow_int_T(z+3./5.,2)+gsl_sf_pow_int(5./7.,2));            
}

/*==================================================================*/
/* F(z)=z/(z^2-a^2)  with a=3/5  */
T<double> fz46(T<double> z) {
	return z/(pow_int_T(z,2)-gsl_sf_pow_int(3./5.,2));          
}

/*==================================================================*/
/* F(z)=1./(sqrt(z^2+a^2))  with a=3/5 */
T<double> fz47(T<double> z) {
	double a=3./5.;
	return (1./(sqrt(pow_int_T(z,2)+gsl_sf_pow_int(a,2))));
}

/*==================================================================*/
/* F(z)=atan(k/z) with k=9/11 */
T<double>  fz48(T<double> z) {
    double k=9./11.;
    return atan(k/z);
}

/*==================================================================*/
/* F(z)= 1./(sqrt(z)*(sqrt(z)+a)) with a=3/5 */
T<double> fz49(T<double> z) {
	double a=3./5.;
	return 1./(sqrt(z)*(sqrt(z)+a));
}

/*==================================================================*/
/* F(z)= (sqrt(z+2a)-sqrt(z))/(sqrt(z+2a)+sqrt(z)) with a=3/5 */	
T<double> fz50(T<double> z) {
	double a=3./5.;
	T<double> sz2a=sqrt(z+2*a);
	T<double> sz=sqrt(z);
	return (sz2a-sz)/(sz2a+sz);
}

/*==================================================================*/
/*F(z)=(((z-1)/z)^n)/z with n=5   */
T<double> fz51(T<double> z){
	int n=5;
	T<double> tmp=pow_int_T((z-1)/z,n);
	return tmp/z;
}
		
/*======================================================================*/
/* F(z)=1/(z^k)  with k=9/11   */
T<double> fz52(T<double> z) {
	double k=9./11.;
	return 1./(pow(z,k));
}

/*==================================================================*/
/* F(z)=1/((z+a)^k)  with a=3/5, k=9/11 */
T<double>  fz53(T<double> z) {
	double k=9./11.;
	double a=3./5.;
   	return 1./pow((z+a),k);
}

/*======================================================================*/
/* F(z)= 1./sqrt(z)  */
T<double> fz54(T<double> z) {
	return 1./sqrt(z);
}

/*==================================================================*/
/* F(z)= z./(z+a)^(3/2)  with a=3/5 */
T<double> fz55(T<double> z) {
	return z/pow(z+3./5.,3./2.);
}

/*==================================================================*/
/* F(z)=sqrt(z+a)-sqrt(z+b)  with a=3/5, b=5/7 */
T<double> fz56(T<double> z) {
	return sqrt(z+3./5.)-sqrt(z+5./7.);
}

/*==================================================================*/
/* F(z)= 1./(sqrt(z)+a)  with a=3/5  */
T<double> fz57(T<double> z) {
	return 1./(sqrt(z)+3./5.);
}

/*==================================================================*/
/* F(z)= sqrt(z)/(z-a^2)  with a=3/5   */
/*T<double> fz59(T<double> z) {
	return sqrt(z)/(z-gsl_sf_pow_int(3./5.,2));
}
*/

/*==================================================================*/
/* F(z)=(1./(sqrt(z+a)*sqrt(z+b)))-1  with a=3/5, b=5/7   */
/*T<double> fz60(T<double> z) {
	double a=3./5.;
	double b=5./7.;
	return (1./(sqrt(z+a)*sqrt(z+b)))-1.;
}
*/

/*==================================================================*/
/*F(z)=((1-z)^n)/(z^(n+1/2)) with n=5   */
T<double> fz58(T<double> z){
	int n=5;
	if (z[0]==0) return (T<double>) 1;
	T<double> tmp=pow_int_T(1-z,n);
	return tmp/pow(z,n+1./2.);
}

/*==================================================================*/
/*F(z)= 1/sqrt(z+1)*/
T<double> fz59(T<double> z) {
 	double a=1;
 	return 1./sqrt(z+a);
}

/*======================================================================*/
/* F(z)=1/z*(1-log(1+z)/z)*/ 
T<double> fz60(T<double> z) { 
	if (z[0]==0) return 1;
	return 1./z*(1.-log(1.+z)/z);
}

/*==================================================================*/
/*F(z)=log(z)/z */       
T<double> fz61(T<double> z){
	if (z[0]==0) return (T<double>) 1;
	else return log(z)/z;
}

/*==================================================================*/
/* F(z)=log(1+z)/z, k=1*/ 
T<double> fz62(T<double> z) { 
	if (z[0]==0) return 1;
	return log(1.+z)/z;
}

/*==================================================================*/
/* F(z)= log((z+a)/(z+b))  with a=3/5, b=5/7 */
T<double> fz63(T<double> z) {
	return log((z+3./5.)/(z+5./7.)) ;
}

/*==================================================================*/
/* F(z)= log((z^2+a^2)/(z^2))  with a=3/5 */
T<double> fz64(T<double> z) {
	return log((pow_int_T(z,2)+gsl_sf_pow_int(3./5.,2))/pow_int_T(z,2)) ;           
}

/*==================================================================*/
/* F(z)= log((z^2-a^2)/(z^2))  with a=3/5   */
T<double> fz65(T<double> z) {
	return log((pow_int_T(z,2)-gsl_sf_pow_int(3./5.,2))/pow_int_T(z,2)) ;            
}

/*======================================================================*/
/* F(z)=(sqrt(z+2a)/sqrt(z))-1  with a=3/5 */
T<double> fz66(T<double> z) {
	double a=3./5.;
	return (sqrt(z+2*a)/sqrt(z))-1.;
}

/*======================================================================*/
/* F(z)= exp(k/z)/(z^mu) with mu=4, k=9/11  */
T<double> fz67(T<double> z) {
	int mu=4;
	double k=9./11.;
	return exp(k/z)/(pow_int_T(z,mu));
}

/*==================================================================*/
/* F(z)=exp(k/z)/z^(3/2) with k=9/11  */
T<double>  fz68(T<double> z) {
	double k=9./11.;
	double aa=3./2.;
	return exp(k/z)/pow(z,aa);
}

/*==================================================================*/
/* F(z)=exp(k/z)/sqrt(z) with k=9/11  */
T<double>  fz69(T<double> z) {
	double k=9./11.;
	return exp(k/z)/sqrt(z);
}

/*==================================================================*/
/* F(z)=exp(-k*sqrt(z)) with k=9/11 */
T<double>  fz70(T<double> z) {
	double k=9./11.;
	return exp(-k*sqrt(z));
}

/*==================================================================*/
/*F(z)=(z^(-1./2.-n))/exp(k*sqrt(z))  n=5, k=9/11  */
T<double> fz71(T<double> z){
	double n=5;
	double k=9./11.;
	double esp=-1./2.-n;
	return pow(z,esp)/exp(k*sqrt(z));
}

/*==================================================================*/
/* F(z)=1./((z^(3/2))*exp(k*sqrt(z))) with k=9/11   */
T<double>  fz72(T<double> z) {
    double k=9./11.;
	double a=3./2.;
    return 1./(pow(z,a)*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=a/((z*(sqrt(z)+a))*exp(k*sqrt(z))) with a=3/5, k=9/11 */			
T<double>  fz73(T<double> z) {
	double k=9./11.;
	double a=3./5.;
	return a/((z*(sqrt(z)+a))*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=1./(z*exp(k*sqrt(z))) with k=9/11 */				
T<double>  fz74(T<double> z) {
	double k=9./11.;
	return 1./(z*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=1/(z+0.5)*exp(-(2*z-1)/(2*z+1)) */ 
T<double> fz75(T<double> z) { 
	if (z[0]==0) return 1;
	return 1./(z+0.5)*exp(-(2.*z-1.)/(2.*z+1.));
}

/*==================================================================*/
/* F(z)=1/((sqrt(z)*(sqrt(z)+a))*exp(k*sqrt(z))) with a=3/5, k=9/11*/		
T<double>  fz76(T<double> z) {
	double k=9./11.;
	double a=3./5.;
	return 1./((sqrt(z)*(sqrt(z)+a))*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=1/((sqrt(z^2+a^2))*exp(k*(sqrt(z^2+a^2)-z))) with a=3/5, k=9/11*/
T<double> fz77(T<double> z) {
	double k=9./11.;
	double a=3./5.;
	T<double> x=sqrt(pow_int_T(z,2)+gsl_sf_pow_int(a,2));
	T<double> den=exp(k*(x-z));
	return 1./(x*den);
}

/*==================================================================*/
/* F(z)=1./(sqrt(z)*exp(k*sqrt(z))) with k=9/11 */
T<double>  fz78(T<double> z) {
	double k=9./11.;
	return 1./(sqrt(z)*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=1./((sqrt(z)+a)*exp(k*sqrt(z))) with a=3/5, k=9/11 */
T<double>  fz79(T<double> z) {
	double k=9./11.;
	double a=3./5.;
	return 1./((sqrt(z)+a)*exp(k*sqrt(z)));
}

/*==================================================================*/
/* F(z)=1/sqrt(1+z)*exp(-z/(1+z))*/ 
T<double> fz80(T<double> z) { 
	return 1./sqrt(1.+z)*exp(-z/(1.+z));
}

/*==================================================================*/
/* F(z)= 1/((z^mu)*exp(k/z)) with mu=4, k=9/11   */
T<double> fz81(T<double> z) {
	int mu=4;
	double k=9./11.;
	return 1./((pow_int_T(z,mu))*exp(k/z));
}

/*==================================================================*/
/* F(z)=1./(exp(k/z)*(z^(3/2))) with k=9/11   */
T<double>  fz82(T<double> z) {
	double k=9./11.;
	return 1./(exp(k/z)*pow(z,3./2.));
}

/*==================================================================*/
/* F(z)=1/(z*exp(k/z)) with k=9/11 */ 
T<double>  fz83(T<double> z) { 
	double k=9./11.;
	return 1./(z*exp(k/z));
}

/*==================================================================*/
/* F(z)=1/(sqrt(z)*exp(k/z)) with k=9/11   */
T<double>  fz84(T<double> z) {
	double k=9./11.;
	return 1./(sqrt(z)*exp(k/z));
}


/*======================================================================*/
/*
NEXT 5 are the Laplace Transform of
function TEST 1, 2, 3, 4, 5 of the paper
"On the Numerical Inversion of Laplace Transforms: Comparison of Three
New Methods on Characteristic Problems from Applications".
DEAN G. DUFFY
*/
/*======================================================================*/
/*F(z)=(1/(z*z+z))*[(1.0/(2*z))-exp(-2*z)/(1-exp(-2*z))]*/
T<double> fz85(T<double> z) {
	T<double>  B = (1.0/(2.*z))-(exp(-2.*z)/(1.-exp(-2.*z)));
	return (1./(z*z+z))*B;
} 

 /*==================================================================*/

T<double> fz86(T<double> z) { 
	double a=100; 
    T<double> sqrt1=sqrt(z);
    T<double> sqrtd2=sqrt(z)/2; 
    return ((a*z-1.)*(exp(sqrtd2)-exp(-sqrtd2))/2)/(z*(z*(exp(sqrt1)-exp(-sqrt1))/2+sqrt1*(exp(sqrt1)+exp(-sqrt1))/2));
} 

/*==================================================================*/
 /*F(z)= (1/z)*exp(-r*sqrt((z*(1+z))/(1+c*z))),  r=0.5, c=0.4*/
T<double> fz87(T<double> z) { 
	double r=0.5;
	double c=0.4;
	T<double> B =-r*sqrt((z*(1+z))/(1+c*z));
	return (1./z)*exp(B);
} 

/*==================================================================*/
/*F(z)= 1/(z*(sqrt(1+z^2+z^4/16)+(z/4)*sqrt(16+z^2))^2)*/
T<double> fz88(T<double> z) {
	T<double> Z2=z*z;
	T<double> Z4=Z2*Z2;
	T<double> ARG=1.+Z2+Z4/16.;
	T<double> A1=sqrt(ARG);
	T<double> A2=(z/4.)*sqrt(16.+Z2);
	return 1./(z*pow_int_T((A1+A2),2));
}

/*==================================================================*/
/*F(z)= (z-sqrt(z^2-1))/(sqrt(z)*sqrt(z^2-1)*sqrt(z-0.5*sqrt(z^2-1)))*/
T<double> fz89(T<double> z) {  
	T<double> B=z-sqrt(z*z-1.);
	T<double> BB=sqrt(z)*sqrt(z*z-1.)*sqrt(z-0.5*(sqrt(z*z-1.)));
	return B/BB;
}
