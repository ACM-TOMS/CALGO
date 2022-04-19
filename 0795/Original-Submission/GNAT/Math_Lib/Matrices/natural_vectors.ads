with Vectors;
with Natural_Instantiation_Parameters;  use Natural_Instantiation_Parameters;

package Natural_Vectors is new Vectors
     (natural,0,clear,copy,equal,"+","-","-","*",
      Plus_Nat,Min_Nat,Min_Nat,Mult_Nat);
