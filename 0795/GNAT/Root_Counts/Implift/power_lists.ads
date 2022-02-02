with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Complex_Multivariate_Polynomials; 
with Complex_Multivariate_Laurent_Polynomials;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;

package Power_Lists is

-- DESCRIPTION :
--   This package offers routines dealing with the supports of polynomials.

  function Construct_Power_List
                ( p : Complex_Multivariate_Polynomials.Poly ) return List;
  function Construct_Power_List
                ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
                return List;

  -- DESCRIPTION : Returns the support of p.

  function Select_Terms
                ( p : Complex_Multivariate_Polynomials.Poly; l : List )
                return Complex_Multivariate_Polynomials.Poly;
  function Select_Terms
                ( p : Complex_Multivariate_Laurent_Polynomials.Poly; l : List )
                return Complex_Multivariate_Laurent_Polynomials.Poly;

  -- DESCRIPTION :
  --   Returns those terms in p whose vector of powers occurs in the list l.

  function Construct_Power_Lists 
                ( p : Poly_Sys ) return Array_of_Lists;
  function Construct_Power_Lists 
                ( p : Laur_Sys ) return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns the supports of the (Laurent) polynomial systems.

  function Select_Terms 
                ( p : Poly_Sys; al : Array_of_Lists ) return Poly_Sys;
  function Select_Terms 
                ( p : Laur_Sys; al : Array_of_Lists ) return Laur_Sys;

  -- DESCRIPTION :
  --   Returns those terms in each polynomial p(i) whose vector of powers
  --   occurs in the list al(i), for i in p'range = al'range.

end Power_Lists;
