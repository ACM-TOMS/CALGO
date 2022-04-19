with text_io,Solutions;               use text_io,Solutions;
with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;

package Drivers_for_Implicit_Lifting is

  procedure Implicit_Lifting_Info;

  -- DESCRIPTION :
  --   Displays information on the BKK Bound on screen.

  procedure Driver_for_Mixture_Bezout_BKK
                ( file : in file_type; p : in Poly_Sys; byebye : in boolean;
                  q : out Poly_Sys; qsols : out Solution_List;
                  b : out natural );

  -- DESCRIPTION :
  --   This driver allows to compute a mixture between a generalized Bezout
  --   number and the BKK bound.

  -- ON ENTRY :
  --   file       to write diagnostics on;
  --   p          a polynomial system;
  --   byebye     if true, then a bye-bye message will appear on screen,
  --              if false, then no bye-bye.

  -- ON RETURN :
  --   q          a start system;
  --   qsols      the solutions of q;
  --   b          Bezout-BKK bound for the number of finite solutions of p.

end Drivers_for_Implicit_Lifting;
