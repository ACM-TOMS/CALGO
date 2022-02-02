package Integer_Instantiation_Parameters is

-- DESCRIPTION :
--   This package contains the routines necessary for
--   instantiating the packages Vectors and Polynomials.

  procedure clear   ( a : in out integer );
  procedure copy    ( a : in integer; b : in out integer );
  function  equal   ( a,b : integer ) return boolean;
  function  convert ( n : natural )   return integer;

  procedure Plus_Int ( a : in out integer; b : in integer );  -- a := a+b;
  procedure Min_Int  ( a : in out integer; b : in integer );  -- a := a-b;
  procedure Min_Int  ( a : in out integer );                  -- a := -a;
  procedure Mult_Int ( a : in out integer; b : in integer );  -- a := a*b;

end Integer_Instantiation_Parameters;
