with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Solutions;                         use Solutions;

package Projective_Transformations is

-- DESCRIPTION :
--   This package provides routines for performing a projective
--   transformation on polynomials and polynomial systems, on
--   solutions and on lists of solutions.

  function Projective_Transformation ( p : Poly ) return Poly;
  function Projective_Transformation ( p : Poly_Sys ) return Poly_Sys;

  procedure Projective_Transformation ( p : in out Poly );
  procedure Projective_Transformation ( p : in out Poly_Sys );

  -- DESCRIPTION :
  --   An additional unknown is added so that each term has the same degree.

  function  Projective_Transformation ( s : Solution) return Solution;
  procedure Projective_Transformation ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   An additional component, equal to 1, is added to the solution vector.

  function  Affine_Transformation ( s : Solution ) return Solution;
  procedure Affine_Transformation ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   All components of the solution vector will be divided by the last
  --   component, which is afterwards cut off.

end Projective_Transformations;
