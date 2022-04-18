with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;

package Polynomial_to_Laurent_Converters is

-- DESCRIPTION :
--   This package contains routines for converting ordinary polynomials
--   into Laurent polynomials.

  function Polynomial_to_Laurent_Polynomial
             ( p : Complex_Multivariate_Polynomials.Poly )
             return Complex_Multivariate_Laurent_Polynomials.Poly;

  -- DESCRIPTION :
  --   Transforms a polynomial into a Laurent polynomial.

  function Polynomial_to_Laurent_System ( p : Poly_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Transforms a polynomial system into a Laurent polynomial system.

end Polynomial_to_Laurent_Converters;
