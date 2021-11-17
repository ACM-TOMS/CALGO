with text_io,Solutions;           use text_io,Solutions;
with Complex_Polynomial_Systems;  use Complex_Polynomial_Systems;

package Drivers_for_m_Homogenization is

  procedure m_Homogenization_Info;

  -- DESCRIPTION :
  --   Displays information about m-homogenization on screen.

  procedure Driver_for_m_Homogenization
                ( file : in file_type; p : in Poly_Sys; b : in out natural;
                  q : out Poly_Sys; qsols : out Solution_List );

  -- DESCRIPTION :
  --   Computation of an m-homogeneous Bezout number with the option
  --   of constructing an m-homogeneous start system.

  -- ON ENTRY :
  --   file       output file;
  --   p          a polynomial system.

  -- ON RETURN :
  --   b          m-homogeneous Bezout number;
  --   q          m-homogeneous start system;
  --   qsols      solutions of q.

end Drivers_for_m_Homogenization;
