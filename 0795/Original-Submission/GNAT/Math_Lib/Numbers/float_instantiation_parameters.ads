with Floating_Point_Numbers;     use Floating_Point_Numbers;

package Float_Instantiation_Parameters is

-- DESCRIPTION :
--   This package contains the routines necessary for
--   instantiating the packages Vectors and Polynomials.

  procedure clear   ( a : in out single_float );
  procedure clear   ( a : in out double_float );
  procedure copy    ( a : in single_float; b : in out single_float );
  procedure copy    ( a : in double_float; b : in out double_float );
  function  equal   ( a,b : single_float ) return boolean;
  function  equal   ( a,b : double_float ) return boolean;
  function  convert ( n : integer ) return single_float;
  function  convert ( n : integer ) return double_float;

  procedure Plus_Float ( a : in out single_float; b : in single_float );
  procedure Plus_Float ( a : in out double_float; b : in double_float ); 
                                                               -- a := a+b;
  procedure Min_Float  ( a : in out single_float; b : in single_float );
  procedure Min_Float  ( a : in out double_float; b : in double_float );
                                                               -- a := a-b;
  procedure Min_Float  ( a : in out single_float );             
  procedure Min_Float  ( a : in out double_float );
                                                               -- a := -a;
  procedure Mult_Float ( a : in out single_float; b : in single_float );
  procedure Mult_Float ( a : in out double_float; b : in double_float );
                                                               -- a := a*b;
  procedure Div_Float  ( a : in out single_float; b : in single_float ); 
  procedure Div_Float  ( a : in out double_float; b : in double_float );
                                                               -- a := a/b;

end Float_Instantiation_Parameters;
