with Multivariate_Laurent_Polynomials;
with Complex_Instantiation_Parameters;  use Complex_Instantiation_Parameters;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Integer_Graded_Lexicographical_Ordening; 
 use Integer_Graded_Lexicographical_Ordening;

package Complex_Multivariate_Laurent_Polynomials is 
   new Multivariate_Laurent_Polynomials
     (double_complex,Vector,CMPLX(0.0),clear,copy,equal,"+","-","-","*","/",
      convert,Plus_Cmplx,Min_Cmplx,Min_Cmplx,Mult_Cmplx,Div_Cmplx,"<",">");
