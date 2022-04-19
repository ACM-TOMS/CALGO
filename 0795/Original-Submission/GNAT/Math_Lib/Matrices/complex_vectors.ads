with Vectors,Complex_Numbers;           use Complex_Numbers;
with Complex_Instantiation_Parameters;  use Complex_Instantiation_Parameters;

package Complex_Vectors is new Vectors
     (double_complex,CMPLX(0.0),clear,copy,equal,"+","-","-","*",
      Plus_Cmplx,Min_Cmplx,Min_Cmplx,Mult_Cmplx);
