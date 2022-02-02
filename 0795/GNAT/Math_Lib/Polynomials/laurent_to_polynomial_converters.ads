with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;

package Laurent_to_Polynomial_Converters is

-- DESCRIPTION :
--   This package contains routines for converting Laurent polynomials to
--   polynomials with natural exponent vectors by shifting when necessary.

  function Laurent_Polynomial_to_Polynomial
              ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
              return Complex_Multivariate_Polynomials.Poly;

  procedure Laurent_Polynomial_to_Polynomial
             ( l : in Complex_Multivariate_Laurent_Polynomials.Poly;
               t : out Complex_Multivariate_Laurent_Polynomials.Term;
               p : out Complex_Multivariate_Polynomials.Poly );

  function Laurent_Polynomial_to_Polynomial
	     ( l : Complex_Multivariate_Laurent_Polynomials.Poly;
               t : Complex_Multivariate_Laurent_Polynomials.Term )
             return Complex_Multivariate_Polynomials.Poly;

  -- DESCRIPTION :
  --   Transforms a Laurent polynomial into an ordinary polynomial
  --   by multiplying by an appropriate monomial t.

  function Laurent_to_Polynomial_System ( p : Laur_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Transforms a Laurent polynomial system into a polynomial system.

end Laurent_to_Polynomial_Converters;
