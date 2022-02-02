with text_io,Solutions;           use text_io,Solutions;
with Complex_Polynomial_Systems;  use Complex_Polynomial_Systems;
with Lists_of_Integer_Vectors;    use Lists_of_Integer_Vectors;

package Drivers_for_Set_Structures is

  procedure Set_Structure_Info;

  -- DESCRIPTION :
  --   Displays information on set structures on screen.

  procedure Driver_for_Set_Structure
               ( file : in file_type; p : in Poly_Sys;
                 b : in out natural; lpos : in out List;
                 q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Allows the interactive computation of a generalized Bezout number,
  --   with an optional construction of a start system.

  -- ON ENTRY :
  --   file      output file;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the set structure;
  --   lpos      a list of positions indicating the acceptable classes;
  --   q         a random product start system;
  --   qsols     the solutions of q.

end Drivers_for_Set_Structures;
