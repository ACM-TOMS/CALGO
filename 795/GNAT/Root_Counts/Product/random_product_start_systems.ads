with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;
with Solutions;                       use Solutions;

package Random_Product_Start_Systems is

-- DESCRIPTION :
--   This package contains a routine for the construction
--   of random product start systems.

  procedure RPQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List; nl : out natural );

  -- DESCRIPTION :
  --   This routine constructs a start system q with the same structure
  --   as the system p.
  --   Solving q happens by factoring all possible linear systems.

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   q          a suitable start system for p;
  --   sols       the solutions of the start system;
  --   nl         the number of matrices that are factored.

  procedure GBQ ( p : in Poly_Sys; q : out Poly_Sys;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   This routine constructs a start system q with the same structure
  --   as the system p.
  --   Solving q happens by factoring only these linear systems that
  --   correspond to admissible products.

  -- ON ENTRY :
  --   p          a polynomial system.

  -- ON RETURN :
  --   q          a suitable start system for p;
  --   sols       the solutions of the start system.

end Random_Product_Start_Systems;
