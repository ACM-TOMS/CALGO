with Integer_Matrices;                  use Integer_Matrices;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;

package Degrees_in_Sets_of_Unknowns is

-- DESCRIPTION :
--   This procedure provides routines for computing the degree
--   of a given complex polynomial w.r.t. to a given set.

  function Degree ( t : Term; s : Set ) return integer;
  function Degree ( p : Poly; s : Set ) return integer;

  -- DESCRIPTION :
  --   Returns the degree in the given set.

  function Degree_Table ( p : Poly_Sys; z : Partition ) return matrix;

  -- DESCRIPTION :
  --   The element (i,j) of the returned matrix contains Degree(p(i),z(j)).

end Degrees_in_Sets_of_Unknowns;
