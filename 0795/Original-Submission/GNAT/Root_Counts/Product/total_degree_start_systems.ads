with Complex_Polynomial_Systems;   use Complex_Polynomial_Systems;
with Complex_Vectors,Solutions;    use Complex_Vectors,Solutions;

package Total_Degree_Start_Systems is

-- DESCRIPTION :
--   This package constructs the start system based on the total degree.

  procedure Total_Degree_Info;

  -- DESCRIPTION :
  --   Displays information about the total degree on screen.

  procedure Start_System 
               ( p : in Poly_Sys; q : in out Poly_Sys; c : in Vector;
                 qsols : in out Solution_List);

  -- ON ENTRY :
  --   p         the polynomial system that has to be solved;
  --   c         a vector with constants.

  -- ON RETURN :
  --   q         the start system, with for i in q'range :
  --             q(i) = x_i^Degree(p(i)) - c(i);
  --   qsols     the solutions of the start system.

  procedure Start_System 
               ( p : in Poly_Sys; q : in out Poly_Sys;
                 qsols : in out Solution_List ); 

  -- DESCRIPTION :
  --   The meaning of the parameters is here the same, except that instead
  --   of the vector c, a random number generator will be used for choosing
  --   the c(i)'s on the unit circle.

end Total_Degree_Start_Systems;
