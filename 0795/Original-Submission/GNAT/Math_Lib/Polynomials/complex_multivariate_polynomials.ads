with Multivariate_Polynomials;
with Complex_Instantiation_Parameters;  use Complex_Instantiation_Parameters;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Natural_Graded_Lexicographical_Ordening; 
 use Natural_Graded_Lexicographical_Ordening;

package Complex_Multivariate_Polynomials is 
   new Multivariate_Polynomials
     (double_complex,Vector,CMPLX(0.0),clear,copy,equal,"+","-","-","*",
      convert,Plus_Cmplx,Min_Cmplx,Min_Cmplx,Mult_Cmplx,"<",">");
