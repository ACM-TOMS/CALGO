package body Float_Instantiation_Parameters is

  procedure clear ( a : in out single_float ) is
  begin
    null;
  end clear;

  procedure clear ( a : in out double_float ) is
  begin
    null;
  end clear;

  procedure copy ( a : in single_float; b : in out single_float ) is
  begin
    b := a;
  end copy;

  procedure copy ( a : in double_float; b : in out double_float ) is
  begin
    b := a;
  end copy;

  function  equal   ( a,b : single_float ) return boolean is
  begin
    return a=b;
  end equal;

  function  equal   ( a,b : double_float ) return boolean is
  begin
    return a=b;
  end equal;

  function  convert ( n : integer )   return single_float is
  begin
    return single_float(n);
  end convert;

  function  convert ( n : integer )   return double_float is
  begin
    return double_float(n);
  end convert;

  procedure Plus_Float ( a : in out single_float; b : in single_float ) is
  begin
    a := a + b;
  end Plus_Float;

  procedure Plus_Float ( a : in out double_float; b : in double_float ) is
  begin
    a := a + b;
  end Plus_Float;

  procedure Min_Float  ( a : in out single_float; b : in single_float ) is
  begin
    a := a - b;
  end Min_Float;

  procedure Min_Float  ( a : in out double_float; b : in double_float ) is
  begin
    a := a - b;
  end Min_Float;

  procedure Min_Float  ( a : in out single_float ) is
  begin
    a := -a;
  end Min_Float;

  procedure Min_Float  ( a : in out double_float ) is
  begin
    a := -a;
  end Min_Float;

  procedure Mult_Float ( a : in out single_float; b : in single_float ) is
  begin
    a := a*b;
  end Mult_Float;

  procedure Mult_Float ( a : in out double_float; b : in double_float ) is
  begin
    a := a*b;
  end Mult_Float;

  procedure Div_Float ( a : in out single_float; b : in single_float ) is
  begin
    a := a/b;
  end Div_Float;

  procedure Div_Float ( a : in out double_float; b : in double_float ) is
  begin
    a := a/b;
  end Div_Float;

end Float_Instantiation_Parameters;
