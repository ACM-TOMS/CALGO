package Greatest_Common_Divisor is

 -- DESCRIPTION :
 --   This package contains two routines for the computation
 --   of the greatest common divisor of two integer numbers;
 --   in addition, a routine for the computation of the least
 --   common multiple is supplied.

  function gcd ( a,b : integer ) return integer;

   -- DESCRIPTION :
   --   returns the greatest common divisor of a and b.

  function lcm ( a,b : integer ) return integer;

   -- DESCRIPTION :
   --   returns the least common multiple of a and b.

  procedure gcd ( a,b : in integer; k,l,d : out integer );

   -- DESCRIPTION :
   --   Computes the greatest common divisor d of a and b;
   --   After gcd(a,b,k,l,d), there holds: k*a + l*b = d.

end Greatest_Common_Divisor;
