with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;

package Polynomial_Randomizers is

-- DESCRIPTION :
--   This package offers routines for randomizing and perturbing
--   the coefficients of polynomials.
--   Except for the last three functions, the monomial structure 
--   remains the same, only random (real or complex) coefficients 
--   will replace the existing ones.

  function Real_Randomize ( p : Poly ) return Poly;
  function Real_Randomize ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Generates random real coefficients in [-1.0,1.0].

  function Complex_Randomize ( p : Poly ) return Poly;
  function Complex_Randomize ( p : Poly_Sys ) return Poly_Sys;
 
  -- DESCRIPTION :
  --   The real and imaginary part of the randomly generated
  --   coefficients are in [-1.0,1.0]

  function Complex_Randomize1 ( p : Poly ) return Poly;
  function Complex_Randomize1 ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   Generates random complex coefficients with modulus one.

  function Complex_Perturb ( p : Poly ) return Poly;
  function Complex_Perturb ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   The coefficients of the polynomial system will be
  --   perturbed by adding random complex coefficients with modulus one.

  function Randomize_with_Real_Matrix     ( p : Poly_Sys ) return Poly_Sys;
  function Randomize_with_Complex_Matrix  ( p : Poly_Sys ) return Poly_Sys;
  function Randomize_with_Complex_Matrix1 ( p : Poly_Sys ) return Poly_Sys;

  -- DESCRIPTION :
  --   First a dense matrix with randomly generated (real, complex or 
  --   complex with modulus 1) coefficients will be constructed.
  --   This matrix will then be used for multiplication with p.
  --   Note that the solution set of the resulting polynomial system
  --   remains unchanged.

end Polynomial_Randomizers;
