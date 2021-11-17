with Vectors;
with Floating_Point_Numbers;          use Floating_Point_Numbers;
with Float_Instantiation_Parameters;  use Float_Instantiation_Parameters;

package Float_Vectors is new Vectors
     (double_float,0.0,clear,copy,equal,"+","-","-","*",
      Plus_Float,Min_Float,Min_Float,Mult_Float);
