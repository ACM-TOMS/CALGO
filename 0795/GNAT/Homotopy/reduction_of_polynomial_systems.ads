with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Vectors;                   use Complex_Vectors;

package Reduction_of_Polynomial_Systems is

-- DESCRIPTION :
--   This package offers some routines for lowering the Bezout
--   number of a given - probably deficient - polynomial system.

  function Total_Degree ( p : Poly_Sys ) return natural;

  -- DESCRIPTION :
  --   Returns the total degree of the polynomial system,
  --   i.e. the product of the degrees of the polynomials.
 
  procedure Reduce ( p : in out Poly_Sys;
                     diagonal,inconsistent,infinite : in out boolean );

  -- DESCRIPTION :
  --   This procedure tries to lower the total degree of p by means
  --   of linear reduction.

  -- ON ENTRY : 
  --   p             a polynomial system.

  -- ON RETURN :
  --   p             a polynomial system with a possible lower total degree;
  --   diagonal      true if all leading terms in p are different;
  --   inconsistent  is true if the reduced system has equations `4=0';
  --   infinite      is true if some equations of the original system
  --                 disappeared during the reduction process.

  procedure Sparse_Reduce ( p : in out Poly_Sys;
                            inconsistent,infinite : in out boolean );

  -- DESCRIPTION :
  --   This procedure makes the coefficient matrix of p as sparse as
  --   possible.

  procedure Reduce ( p : in Poly_Sys; res : in out Poly_Sys;
                     cnt_eq : in out natural; max_eq : in natural;
                     cnt_sp : in out natural; max_sp : in natural;
                     cnt_rp : in out natural; max_rp : in natural );

  -- DESCRIPTION : 
  --   This procedure tries to lower the total degree of the system p
  --   by means of nonlinear reduction. 

  -- REQUIRED : the counters must equal 0, on entry.

  -- ON ENTRY :
  --   p              a polynomial system;
  --   cnt_eq         counts the number of equal degree substitutions;
  --   max_eq         limit on the number of equal degree substitutions;
  --   cnt_sp         counts the number of S-polynomial computations;
  --   max_sp         limit on the number of S-polynomial computations.
  --   cnt_rp         counts the number of R-polynomial computations;
  --   max_rp         limit on the number of R-polynomial computations.

  -- ON RETURN :
  --   res            the reduced system;
  --   cnt_eq         the number of equal degree substitutions;
  --   cnt_sp         the number of computed S-polynomials;
  --   cnt_rp         the number of computed R-polynomials.

  procedure Sparse_Reduce ( p : in Poly_Sys; res : in out Poly_Sys;
                            cnt_eq : in out natural; max_eq : in natural );

  -- DESCRIPTION :
  --   the polynomial system is reduced by computing S-polynomials.
  --   After each replacement, the coefficient matrix is made as sparse
  --   as possible.

end Reduction_of_Polynomial_Systems;
