with Integer_Vectors,Transformations;    use Integer_Vectors,Transformations;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Complex_Multivariate_Laurent_Polynomials;
 use Complex_Multivariate_Laurent_Polynomials;

package Transforming_Laurent_Systems is

-- DESCRIPTION :
--   This package offers some routines for transforming
--   Laurent polynomials.

  procedure Shift ( p : in out Poly );
  function  Shift ( p : Poly ) return Poly;

  procedure Shift ( l : in out Laur_Sys );
  function  Shift ( l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Shifts the support of the polynomial so that the constant term
  --   belongs to p.
  --   This Shift does not change the term order in p!

  procedure Transform ( t : in Transfo; p : in out Poly );
  function  Transform ( t : Transfo; p : Poly ) return Poly;

  procedure Transform ( t : in Transfo; l : in out Laur_Sys );
  function  Transform ( t : Transfo; l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION : Application of the transformation t.

  function Maximal_Support ( p : Poly; v : Vector ) return integer;
  function Maximal_Support ( p : Poly; v : Link_to_Vector ) return integer;

  -- DESCRIPTION :
  --   Computes the value of the supporting function of the polynomial p,
  --   for the direction v.

  procedure Face ( i,m : in integer; p : in out Poly );
  function  Face ( i,m : integer; p : Poly ) return Poly;

  procedure Face ( i,m : in integer; l : in out Laur_Sys );
  function  Face ( i,m : integer; l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   returns only the terms t for which deg(t,xi) = m.

  procedure Face ( v : in Vector; m : in integer; p : in out Poly );
  function  Face ( v : Vector; m : integer; p : Poly ) return Poly;

  procedure Face ( v,m : in Vector; l : in out Laur_Sys );
  function  Face ( v,m : Vector; l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   Only the terms for which for the degrees d the following holds
  --    < d , v > = m, are left.

  procedure Reduce ( i : in integer; p : in out Poly );
  function  Reduce ( i : integer; p : Poly ) return Poly;

  procedure Reduce ( i : in integer; l : in out Laur_Sys );
  function  Reduce ( i : integer; l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The i-th component out of every monomial will be removed,
  --   so that the polynomials will have an unknown less.

  procedure Insert ( i,d : in integer; p : in out Poly );
  function  Insert ( i,d : integer; p : Poly ) return Poly;

  procedure Insert ( i,d : in integer; l : in out Laur_Sys );
  function  Insert ( i,d : integer; l : Laur_Sys ) return Laur_Sys;

  -- DESCRIPTION :
  --   The i-th component of each monomial will be inserted,
  --   using the value d.

end Transforming_Laurent_Systems;
