with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package Homogenization is

-- DESCRIPTION :
--   This package provides routines for constructing additional
--   equations to a system for projective transformations.
--   There is also a routine that isolates the homogeneous part
--   of a given polynomial system.

  function Homogeneous_Part ( p : Poly ) return Poly;
  function Homogeneous_Part ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   These functions isolate all terms having a degree equal to
  --   the degree of the polynomial.

  function Add_Equations ( s1 : Poly_Sys; s2 : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   The resulting polynomial system is the concatenation of s1 and s2.

  function Add_Equation ( s : Poly_Sys; p : Poly ) return Poly_Sys;

  -- DESCRIPTION :
  --   the resulting polynomial system is the concatenation
  --   of the system s and the polynomial p

  function Add_Random_Hyperplanes
               ( s : Poly_Sys; m : natural; re : boolean ) return Poly_Sys;

  -- DESCRIPTION :
  --   To the polynomial system s, m hyperplanes are added with
  --   randomly choosen coefficients;
  --   if re = true 
  --    then the coefficients will be floating point numbers;
  --    else the coefficients will be complex numbers.

  function Add_Standard_Hyperplanes
               ( s : Poly_Sys; m : natural ) return Poly_Sys;

  -- DESCRIPTION :
  --   If n = Number_Of_Unknowns(s(i)), for i in s'range,
  --    then m hyperplanes of the form
  --           x_(j+n) - 1 = 0 will be added, for j in 1..m,
  --         to the system s.

end Homogenization;
