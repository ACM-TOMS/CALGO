with text_io,Solutions;                  use text_io,Solutions;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Trees_of_Vectors;                   use Trees_of_Vectors;

package Mixed_Homotopy_Continuation is

-- DESCRIPTION :
--   This package provides solvers based on the principle
--   of mixed continuation.

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    bkk : out natural; sols : in out Solution_List );

  -- ON ENTRY :
  --   file         a file for writing intermediate results;
  --   p            a Laurent polynomial system.

  -- ON RETURN :
  --   bkk          the BKK bound of p;
  --   sols         the computed solutions.

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    tv : in Tree_of_Vectors; bkk : out natural;
                    sols : in out Solution_List );

  -- ON ENTRY :
  --   file         a file for writing intermediate results;
  --   p            a Laurent polynomial system;
  --   tv           the tree of vectors containing the useful directions.


  -- ON RETURN :
  --   bkk          the BKK bound of p;
  --   sols         the computed solutions.

end Mixed_Homotopy_Continuation;
