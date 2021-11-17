with Permutations;                        use Permutations;

with Natural_Vectors,Integer_Vectors;
with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Float_Vectors,Complex_Vectors;  

with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;

with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;

package Permute_Operations is

-- DESCRIPTION :
--   This package provides permute operations on vectors,
--   on polynomials and on systems of polynomials.

  function "*" ( p : Permutation; v : Natural_Vectors.Vector )
	       return Natural_Vectors.Vector;

  function "*" ( p : Permutation; v : Integer_Vectors.Vector )
	       return Integer_Vectors.Vector;

  function "*" ( p : Permutation; v : Float_Vectors.Vector )
	       return Float_Vectors.Vector;

  function "*" ( p : Permutation; v : Complex_Vectors.Vector )
	       return Complex_Vectors.Vector;

  -- DESCRIPTION :
  --   returns the result of the permutation of p on the vector v.
  -- REQUIRED :
  --   p'range = v'range

  function Permutable ( v1,v2 : Natural_Vectors.Vector ) return boolean;
  function Permutable ( v1,v2 : Integer_Vectors.Vector ) return boolean;
  function Permutable ( v1,v2 :   Float_Vectors.Vector ) return boolean;
  function Permutable ( v1,v2 : Complex_Vectors.Vector ) return boolean;
  function Permutable ( v1,v2 :   Float_Vectors.Vector; tol : double_float )
                      return boolean;
  function Permutable ( v1,v2 : Complex_Vectors.Vector; tol : double_float )
                      return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists a permutation between the two vectors.
  --   If provided, tol is the tolerance for comparing two numeric values.

  function Sign_Permutable ( v1,v2 : Natural_Vectors.Vector ) return boolean;
  function Sign_Permutable ( v1,v2 : Integer_Vectors.Vector ) return boolean;
  function Sign_Permutable ( v1,v2 :   Float_Vectors.Vector ) return boolean;
  function Sign_Permutable ( v1,v2 : Complex_Vectors.Vector ) return boolean;
  function Sign_Permutable ( v1,v2 :   Float_Vectors.Vector;
                                         tol : double_float ) return boolean;
  function Sign_Permutable ( v1,v2 : Complex_Vectors.Vector; 
                                         tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Also permutations where the sign of one of the components can
  --   be changed, are checked.

  function "*" ( p : Permutation; t : Complex_Multivariate_Polynomials.Term )
	       return Complex_Multivariate_Polynomials.Term;

  function "*" ( p : Permutation; s : Complex_Multivariate_Polynomials.Poly )
	       return Complex_Multivariate_Polynomials.Poly;

  function "*" ( p : Permutation;
                 t : Complex_Multivariate_Laurent_Polynomials.Term )
	       return Complex_Multivariate_Laurent_Polynomials.Term;

  function "*" ( p : Permutation;
                 s : Complex_Multivariate_Laurent_Polynomials.Poly )
               return Complex_Multivariate_Laurent_Polynomials.Poly;

  -- DESCRIPTION :
  --   permutes the unknowns in the term t or the polynonomial s,
  --   according to the permuation p.

  function "*" ( s : Poly_Sys; p : Permutation ) return Poly_Sys;
  function "*" ( s : Laur_Sys; p : Permutation ) return Laur_Sys;

  function "*" ( p : Permutation; s : Poly_Sys ) return Poly_Sys;
  function "*" ( p : Permutation; s : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   s*p permutes the unknowns in the individual polynomials.
  --   p*s permutes the equations in the system.
  --   Watch out for sharing by this second type of operation!

end Permute_Operations;
