package Natural_Instantiation_Parameters is

-- DESCRIPTION :
--   This package contains the routines necessary for
--   instantiating the package Vectors.

  procedure clear   ( a : in out natural );
  procedure copy    ( a : in natural; b : in out natural );
  function  equal   ( a,b : natural ) return boolean;

  procedure Plus_Nat ( a : in out natural; b : in natural );  -- a := a+b;
  procedure Min_Nat  ( a : in out natural; b : in natural );  -- a := a-b;
  procedure Min_Nat  ( a : in out natural );                  -- a := -a;
  procedure Mult_Nat ( a : in out natural; b : in natural );  -- a := a*b;

end Natural_Instantiation_Parameters;
