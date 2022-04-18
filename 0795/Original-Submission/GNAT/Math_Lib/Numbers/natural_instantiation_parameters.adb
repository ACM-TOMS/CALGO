package body Natural_Instantiation_Parameters is

  procedure clear   ( a : in out natural ) is
  begin
    null;
  end clear;

  procedure copy    ( a : in natural; b : in out natural ) is
  begin
    b := a;
  end copy;

  function  equal   ( a,b : natural ) return boolean is
  begin
    return a = b;
  end equal;

  procedure Plus_Nat ( a : in out natural; b : in natural ) is
  begin
    a := a+b;
  end Plus_Nat;

  procedure Min_Nat  ( a : in out natural; b : in natural ) is
  begin
    a := a-b;
  end Min_Nat;

  procedure Min_Nat  ( a : in out natural ) is
  begin
    a := -a;
  end Min_Nat;

  procedure Mult_Nat ( a : in out natural; b : in natural ) is
  begin
    a := a*b;
  end Mult_Nat;

end Natural_Instantiation_Parameters;
