with Complex_Numbers,Complex_Vectors;     use Complex_Numbers,Complex_Vectors;
with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;
with Complex_Matrices;                    use Complex_Matrices;

package Homotopy is

-- DESCRIPTION :
--   This package exports operations on a polynomial homotopy.

-- CONSTRUCTORS :

  procedure Create ( p,q : in Poly_Sys; k : in positive; 
                     a : in double_complex );

  -- DESCRIPTION :
  --   the following artificial homotopy is build :
  --     H(x,t) = a * ((1 - t)^k) * q + (t^k) * p.

  procedure Create ( p : in Poly_Sys; k : in integer );

  -- DESCRIPTION :
  --   Given a polynomial system p with dimension n*(n+1) and
  --   k an index, then t = x_k as continuation parameter.

-- SELECTOR :

  function Homotopy_System return Poly_Sys;

  -- DESCRIPTION :
  --   Returns the homotopy system in the unknowns (x,t).

-- SYMBOLIC ROUTINES :

  function Eval ( t : double_complex ) return Poly_Sys;

  -- DESCRIPTION :
  --   The homotopy is evaluated in t and a polynomial system is returned

  function Diff ( t : double_complex ) return Poly_Sys;

  -- DESCRIPTION :
  --   The homotopy is symbolically differentiated w.r.t. t.

-- NUMERIC ROUTINES :

  function Eval ( x : Vector; t : double_complex ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is evaluated in x and t and a vector is returned.

  function Diff ( x : Vector; t : double_complex ) return Vector;

  -- DESCRIPTION :
  --   The homotopy is differentiated wr.t. t and is evaluated in (x,t).

  function Diff ( x : Vector; t : double_complex ) return matrix;

  -- DESCRIPTION :
  --   The homotopy is differentiated to x and the Jacobi matrix
  --   of H(x,t) is returned.

  function Diff ( x : Vector; t : double_complex; k : natural ) return Vector;

  -- DESCRIPTION :
  --   The returning vector contains all derivatives from the homotopy
  --   to the unknown x_k; note that t = x_n+1.

  function Diff ( x : Vector; t : double_complex; k : natural ) return matrix;

  -- DESCRIPTION :
  --   The Jacobi matrix of the homotopy is returned where the kth
  --   column has been deleted; note that Diff(x,t,n+1) = Diff(x,t).

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   The homotopy is cleared.

 -- pragma inline (Eval,Diff);

end Homotopy;
