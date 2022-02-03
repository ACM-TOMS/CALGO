with text_io,Solutions;                use text_io,Solutions;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;

package Drivers_for_Static_Lifting is

  procedure Static_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on static lifting on screen.

  procedure Driver_for_Mixed_Volume_Computation 
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  mv : out natural );

  -- DESCRIPTION :
  --   Interactive driver for the computation of the mixed volume.

  -- ON ENTRY :
  --   file       output file, must be opened for output;
  --   p          a polynomial system;
  --   byebye     if true, then a bye-bye message will appear on screen,
  --              if false, then no bye-bye.

  -- ON RETURN :
  --   q          a start system with randomly choosen coefficients,
  --              which can be used in a coefficient homotopy;
  --   qsols      the solutions of q;
  --   mv         mixed volume of p and the number of solutions of q.

end Drivers_for_Static_Lifting;
