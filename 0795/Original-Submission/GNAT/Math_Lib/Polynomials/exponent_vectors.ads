with Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;
with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;

package Exponent_Vectors is

-- DESCRIPTION :
--   This package facilitates the management of exponent vectors.

-- DATA STRUCTURE : array of exponent vectors

  type Exponent_Vectors_Array is
    array ( integer range <> ) of Integer_Vectors_of_Vectors.Link_to_Vector;

-- CREATORS :

  function Create ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
                  return Integer_Vectors_of_Vectors.Vector;
  function Create ( p : Complex_Multivariate_Polynomials.Poly )
                  return Integer_Vectors_of_Vectors.Vector;

  -- DESCRIPTION :
  --   The range of the vector on return is 1..Number_of_Terms(p).
  --   This vector contains copies of all exponents of p, as ordered in p.

  function Create ( p : Poly_Sys ) return Exponent_Vectors_Array;
  function Create ( p : Laur_Sys ) return Exponent_Vectors_Array;

-- SELECTOR :

  function Position ( ev : Integer_Vectors_of_Vectors.Vector;
                      v : Integer_Vectors.Vector ) return integer;

  -- DESCRIPTION :
  --   Returns the position of v in the vector ev.
  --   If v does not occur in ev, then ev'last+1 will be returned.

-- DESTRUCTORS :

  procedure Clear ( v : in out Exponent_Vectors_Array );

  -- DESCRIPTION :
  --   Clears the allocated memory.

end Exponent_Vectors;
