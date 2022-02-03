with Complex_Multivariate_Laurent_Polynomials;
 use Complex_Multivariate_Laurent_Polynomials;
with Complex_Laurent_Polynomial_Systems; 
 use Complex_Laurent_Polynomial_Systems;

package Laurent_Polynomial_Randomizers is

-- DESCRIPTION :
--   This package offers routines for `randomizing' polynomials:
--   the monomial structure remains the same,
--   but random (real or complex) coefficients will be generated.

  function Real_Randomize ( p : Poly ) return Poly;
  function Real_Randomize ( p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The random real coefficients lie in [-1.0,1.0].

  function Complex_Randomize ( p : Poly ) return Poly;
  function Complex_Randomize ( p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The real and imaginary parts of the complex random coefficients 
  --   lie in [-1.0,1.0].

  function Complex_Randomize1 ( p : Poly ) return Poly;
  function Complex_Randomize1 ( p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The generated random complex coefficients have modulus 1.

  function Complex_Perturb ( p : Poly ) return Poly;
  function Complex_Perturb ( p : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The coefficients of the polynomial system will be
  --   perturbed by adding random complex coefficients with
  --   modulus one.

end Laurent_Polynomial_Randomizers;
