package body Integer_Instantiation_Parameters is

  procedure clear   ( a : in out integer ) is
  begin
    null;
  end clear;

  procedure copy    ( a : in integer; b : in out integer ) is
  begin
    b := a;
  end copy;

  function  equal   ( a,b : integer ) return boolean is
  begin
    return a = b;
  end equal;

  function  convert ( n : natural )   return integer is
  begin
    return integer(n);
  end convert;

  procedure Plus_Int ( a : in out integer; b : in integer ) is
  begin
    a := a+b;
  end Plus_Int;

  procedure Min_Int  ( a : in out integer; b : in integer ) is
  begin
    a := a-b;
  end Min_Int;

  procedure Min_Int  ( a : in out integer ) is
  begin
    a := -a;
  end Min_Int;

  procedure Mult_Int ( a : in out integer; b : in integer ) is
  begin
    a := a*b;
  end Mult_Int;

end Integer_Instantiation_Parameters;
