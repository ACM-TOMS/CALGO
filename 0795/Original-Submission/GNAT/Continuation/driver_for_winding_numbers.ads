with text_io;                            use text_io;
with Complex_Numbers,Solutions;          use Complex_Numbers,Solutions;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;

procedure Driver_for_Winding_Numbers
             ( file : in file_type; p : in Poly_Sys;
               sols : in out Solution_List );

-- DESCRIPTION :
--   Interactive driver for the computation of winding numbers.
--   The user will be asked to define a homotopy.

-- ON ENTRY :
--   file       file to write intermediate results on;
--   p          a polynomial system;
--   sols       solution list, with t < target value for continuation.

-- ON RETURN :
--   sols       refined solution list with appropriate winding number.
