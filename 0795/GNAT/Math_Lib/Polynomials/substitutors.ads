with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package Substitutors is

-- DESCRIPTION :
--   This package contains routines for substituting 
--   equations into polynomials and polynomial systems.

  function  Substitute ( k : integer; c : double_complex; t : Term )
                       return Term;
  procedure Substitute ( k : in integer; c : in double_complex; 
                         t : in out Term );

  function  Substitute ( k : integer; c : double_complex; p : Poly )
                       return Poly;
  procedure Substitute ( k : in integer; c : in double_complex;
                         p : in out Poly );

  function  Substitute ( k : integer; c : double_complex; p : Poly_Sys )
                       return Poly_Sys;
  procedure Substitute ( k : in integer; c : in double_complex;
                         p : in out Poly_Sys );

   -- DESCRIPTION :
   --   These routines substitute the kth unknown of the term t or
   --   polynomial (system) p by a complex constant c.

   -- ON ENTRY :
   --   k         an unknown in the polynomial p;
   --   c         a complex constant;
   --   t         a term;
   --   p         a polynomial (system).

   -- ON RETURN :
   --   t         a term where the kth unknonw has been replaced by the
   --             complex constant c;
   --   p         a polynomial (system) where the kth unknown has been
   --             replaced by the complex constant c.

  function  Substitute ( k : integer; h : Vector; p : Poly ) return Poly;
  procedure Substitute ( k : in integer; h : in Vector; p : in out Poly );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomial p
  --   by a linear equation defined by h.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   h          a vector h(0..n), n = Number_of_Unknowns(p),
  --              defines h(0) + h(1)*x1 + ... + h(n)*xn;
  --   p          a polynomial.

  -- REQUIRED : h(k) /= 0.

  -- ON RETURN :
  --   p          a polynomial where the kth unknown has been replaced
  --              by a linear equation defined by h.

  function  Substitute ( k : integer; s,p : Poly ) return Poly;
  procedure Substitute ( k : in integer; s : in Poly; p : in out Poly );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomial p
  --   by a polynomial s.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   s          a polynomial so that xk = s;
  --   p          a polynomial.

  -- ON RETURN :
  --   p          a polynomial where the kth unknown has been replaced
  --              by the polynomial s.

  function  Substitute ( k : integer; h : Vector; p : Poly_Sys )
		       return Poly_Sys;
  procedure Substitute ( k : in integer; h : in Vector; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomials in the
  --   system p by a linear equation defined by h.

  -- ON ENTRY :
  --   k          an unknown in the polynomials in the system p;
  --   h          a vector h(0..n), n = Number_of_Unknowns(p(i)),
  --              defines h(0) + h(1)*x1 + ... + h(n)*xn;
  --   p          a polynomial system.

  -- REQUIRED : h(k) /= 0

  -- ON RETURN :
  --   p          a polynomial system where the kth unknown has been replaced
  --              by a linear equation defined by h.

  function  Substitute ( k : integer; s : Poly; p : Poly_Sys ) return Poly_Sys;
  procedure Substitute ( k : in integer; s : in Poly; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   These routines substitute the kth unknown of the polynomials in p
  --   by a polynomial s.

  -- ON ENTRY :
  --   k          an unknown in the polynomial p;
  --   s          a polynomial so that xk = s;
  --   p          a polynomial.

  -- ON RETURN :
  --   p          a polynomial system where the kth unknown has been replaced
  --              by the polynomial s.

end Substitutors;
