with text_io,Solutions;            use text_io,Solutions;
with Complex_Polynomial_Systems;   use Complex_Polynomial_Systems;

package Drivers_for_Multi_Homogenization is

  procedure Multi_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about multi-homogenization on screen.

  procedure Driver_for_Multi_Homogenization
               ( file : in file_type; p : in Poly_Sys; b : in out natural;
                 q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   This is an interactive driver for multi-homogenization.

  -- ON ENTRY :
  --   file      to write diagnostics on;
  --   p         a polynomial system.

  -- ON RETURN :
  --   b         a bound based on the degree structure;
  --   q         a random product start system;
  --   qsols     the solutions of q.

end Drivers_for_Multi_Homogenization;
