with Floating_Point_Numbers;    use Floating_Point_Numbers;

package body Complex_Instantiation_Parameters is

  procedure clear ( a : in out double_complex ) is
  begin
    null;
  end clear;

  procedure copy ( a : in double_complex; b : in out double_complex ) is
  begin
    b := a;
  end copy;

  function  equal   ( a,b : double_complex ) return boolean is
  begin
    return a=b;
  end equal;

  function  convert ( n : integer )   return double_complex is
  begin
    return CMPLX(double_float(n));
  end convert;

  procedure Plus_Cmplx ( a : in out double_complex; b : in double_complex ) is
  begin
    a := a + b;
  end Plus_Cmplx;

  procedure Min_Cmplx  ( a : in out double_complex; b : in double_complex ) is
  begin
    a := a - b;
  end Min_Cmplx;

  procedure Min_Cmplx  ( a : in out double_complex ) is
  begin
    a := -a;
  end Min_Cmplx;

  procedure Mult_Cmplx ( a : in out double_complex; b : in double_complex ) is
  begin
    a := a*b;
  end Mult_Cmplx;

  procedure Div_Cmplx ( a : in out double_complex; b : in double_complex ) is
  begin
    a := a/b;
  end Div_Cmplx;

end Complex_Instantiation_Parameters;
