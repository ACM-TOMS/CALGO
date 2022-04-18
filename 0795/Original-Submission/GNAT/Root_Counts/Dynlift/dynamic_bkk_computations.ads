with text_io,Solutions;              use text_io,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;

package Dynamic_BKK_Bound_Computations is

-- DESCRIPTION :
--   This package exports some routines for computing the BKK bound
--   and solving a random coefficient system by polyhedral continuation.
--   These function are black box routines: the user does not have to
--   worry about intermediate data structures.

  function BKK_by_Dynamic_Lifting ( p : Poly_Sys ) return natural;
  function BKK_by_Dynamic_Lifting ( file : file_type; p : Poly_Sys )
                                  return natural;
  -- DESCRIPTION :
  --   If a file is specified, then the mixed subdivision will be
  --   written on that file.

  function Solve_by_Dynamic_Lifting ( p : Poly_Sys ) return Solution_List;
  function Solve_by_Dynamic_Lifting ( file : file_type; p : Poly_Sys )
                                    return Solution_List;
  -- DESCRIPTION :
  --   If a file is specified, then intermediate results will be
  --   write on that file.

end Dynamic_BKK_Bound_Computations;
