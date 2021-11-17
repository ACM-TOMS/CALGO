with Complex_Numbers;   use Complex_Numbers;

package Complex_Instantiation_Parameters is

-- DESCRIPTION :
--   This package contains the routines necessary for
--   instantiating the packages Vectors and Polynomials.

  procedure clear   ( a : in out double_complex );
  procedure copy    ( a : in double_complex; b : in out double_complex );
  function  equal   ( a,b : double_complex ) return boolean;
  function  convert ( n : integer )   return double_complex;

  procedure Plus_Cmplx ( a : in out double_complex; b : in double_complex );
                                                                 -- a := a+b;
  procedure Min_Cmplx  ( a : in out double_complex; b : in double_complex );
                                                                 -- a := a-b;
  procedure Min_Cmplx  ( a : in out double_complex ); 
                                                                 -- a := -a;
  procedure Mult_Cmplx ( a : in out double_complex; b : in double_complex );
                                                                 -- a := a*b;
  procedure Div_Cmplx  ( a : in out double_complex; b : in double_complex );
                                                                 -- a := a/b;

end Complex_Instantiation_Parameters;
