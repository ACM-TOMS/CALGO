with text_io,Solutions;              use text_io,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;
with Symmetry_Group;                 use Symmetry_Group;

package Symmetric_BKK_Bound_Solvers is

-- DESCRIPTION :
--   This package offers functions for the black box computation of the
--   BKK bound of a given polynomial system.

  function Symmetric_BKK_Solve
              ( p : Poly_Sys; sign : boolean ) return Solution_List;

  function Symmetric_BKK_Solve
              ( p : Poly_Sys; grp : List_of_Permutations; sign : boolean )
              return Solution_List;

  function Symmetric_BKK_Solve 
              ( file : file_type; p : Poly_Sys; sign : boolean ) 
              return Solution_List;

  function Symmetric_BKK_Solve
              ( file : file_type; p : Poly_Sys;
                grp : List_of_Permutations; sign : boolean  )
              return Solution_List;

  -- DESCRIPTION :
  --   This is a black box computation of all generating solutions,
  --   based on the computation of the symmetric mixed subdivision.
  --   If a file is specified, then intermediate results will be
  --   write on that file.
  --   If no group is given, then the group of all permutations is
  --   assumed.  Sign = true means that there is also a sign symmetry
  --   to take into account.

end Symmetric_BKK_Bound_Solvers;
