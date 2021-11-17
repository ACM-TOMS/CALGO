with Complex_Polynomial_Systems;   use Complex_Polynomial_Systems;
with Solutions;                    use Solutions;

package Construct_Random_Product_Start_System is

-- DESCRIPTION :
--   This package constructs a random product start
--   product system for a given polynomial system.

  procedure Build_Set_Structure ( p : in Poly_Sys );

  -- DESCRIPTION :
  --   This is a heuristic procedure for constructing a supporting
  --   set structure of the system p.

  procedure Build_Random_Product_System ( n : in natural );

  -- DESCRIPTION :
  --   Based on the set structure, a random linear-product system
  --   will be constructed.  The result is stored in the internal
  --   data manage by the package Random_Product_System.

  -- REQUIRED :
  --   The set structure may not be empty.

  procedure Construct ( p : in Poly_Sys; q : in out Poly_Sys;
		        sols : in out Solution_List );

  -- DESCRIPTION :
  --   Constructs a start system q, with more or less the same 
  --   structure as p.  A heuristic procedure will be used for
  --   constructing a supporting set structure.

end Construct_Random_Product_Start_System;
