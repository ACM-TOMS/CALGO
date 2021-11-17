with text_io,Solutions;              use text_io,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;

procedure Driver_for_Root_Counts
             ( file : in file_type; p,q : in out Poly_Sys;
               own : in boolean;
               qsols : in out Solution_List; roco : out natural );

-- DESCRIPTION :
--   This procedure implements an interactive driver for several
--   root counting methods.

-- ON ENTRY :
--   file      to write diagnostics on;
--   p         a polynomial system;
--   own       if true, then the user has the possibility to give
--             an own start system.

-- ON RETURN :
--   p         has eventually been made homogeneous;
--   q         a start system based on the chosen root count;
--   qsols     the solutions of q;
--   roco      the root count.
