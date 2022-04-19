with text_io,Solutions;              use text_io,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;

procedure Driver_for_Own_Start_System
             ( file : in file_type; p : in Poly_Sys;
               q : out Poly_Sys; qsols : in out Solution_List );

-- DESCRIPTION :
--   This procedure implements an interactive driver for reading
--   a start system delivered by user.

-- ON ENTRY :
--   file      to write diagnostics on;
--   p         a polynomial system.

-- ON RETURN :
--   q         a start system based on the chosen root count;
--   qsols     the solutions of q.
