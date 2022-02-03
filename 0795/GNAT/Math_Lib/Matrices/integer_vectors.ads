with Vectors;
with Integer_Instantiation_Parameters;  use Integer_Instantiation_Parameters;

package Integer_Vectors is new Vectors
     (integer,0,clear,copy,equal,"+","-","-","*",
      Plus_Int,Min_Int,Min_Int,Mult_Int);
