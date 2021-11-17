with Complex_Vectors,Solutions;         use Complex_Vectors,Solutions;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;

package Random_Product_System is

-- DESCRIPTION :
--   This package enables the construction and the solution
--   of a polynomial system of which each polynomial is 
--   represented as a product of hyperplanes.

-- CONSTRUCTORS :

  procedure Init ( n : in natural );

  -- DESCRIPTION :
  --   The internal data for this package are initialised
  --   with n, this is the number of equations of the system
  -- NOTE :
  --   This operation must be the first one executed when
  --   using this package.

  procedure Add_Hyperplane ( i : in natural; h : in Vector );

  -- DESCRIPTION :
  --   The hyperplane h is added to the i-the equation of the
  --   random product system                   
  -- REQUIRED :                               __n_
  --   i <= n                                 \
  --   h : Vector(0..n) representing  h(0) +   >    h(j) x
  --                                          /___        j
  --                                           j=1

  function Dimension return natural;

  -- DESCRIPTION :
  --   returns the number of equations in the product system.

  function Number_of_Hyperplanes (i : natural) return natural;

  -- DESCRIPTION :
  --   returns the number of added hyperplanes for the i-th equation

  function Get_Hyperplane ( i,j : in natural ) return Vector;
  function Get_Hyperplane ( i,j : in natural ) return Link_to_Vector;

  -- DESCRIPTION :
  --   returns the j-th hyperplane for the i-th equation

  procedure Change_Hyperplane ( i,j : in natural; h : in Vector );

  -- DESCRIPTION :
  --   the (i,j)-th hyperplane will be changed into h

  procedure Solve ( sols : in out Solution_List; nl : out natural );

  -- DESCRIPTION :
  --   The random product system is solved and the solutions are
  --   put in the list sols.
  --   nl is the number of matrices that are factored.
  -- NOTE :
  --   All possible linear systems are factorized using Gaussian
  --   elimination, together with the estimation of the condition
  --   of the matrices.
  --   Systems with a bad condition are not solved.

  procedure Solve ( sols : in out Solution_List; 
                    nl : out natural; l : in List );

  -- DESCRIPTION :
  --   Cf. Solve, but only those linear systems are factorized,
  --   for which the linear equations are as indicated in the list
  --   of positions l.

  function Polynomial ( h : Vector ) return Poly;

  -- DESCRIPTION :
  --   returns the linear polynomial defined by the coefficients in h.

  function Polynomial_System return Poly_Sys;

  -- DESCRIPTION :
  --   A polynomial system is constructed by multiplying all
  --   the hyperplanes from the equations of the random product system.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   This procedure frees all memory space used by the random
  --   product system.

end Random_Product_System;
