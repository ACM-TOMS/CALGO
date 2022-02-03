with text_io,Solutions;                use text_io,Solutions;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;

package Drivers_for_Symmetric_Lifting is

  procedure Symmetric_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on symmetric lifting on screen.

  procedure Driver_for_Symmetric_Mixed_Volume_Computation 
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  mv : out natural );

  -- DESCRIPTION :
  --   This procedure presents an interactive driver for the computation
  --   of the mixed volume of a symmetric polynomial system, based on the
  --   construction of a symmetric mixed subdivision.

  -- ON ENTRY :
  --   file       a file to put useful statistics on;
  --   p          a polynomial system;
  --   byebye     if true, then a bye-bye message will appear on screen,
  --              if false, then no bye-bye.

  -- ON RETURN :
  --   q          a start system with randomly choosen coefficients,
  --              which can be used in a coefficient homotopy;
  --   qsols      the solutions of q;
  --   mv         the BKK bound of p and the number of solutions of q.

end Drivers_for_Symmetric_Lifting;
