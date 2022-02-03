with text_io,Solutions;                use text_io,Solutions;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;
with Partitions_of_Sets_of_Unknowns;   use Partitions_of_Sets_of_Unknowns;

procedure Driver_for_Interpolating_Homotopies
              ( file : in file_type; p : in Poly_Sys; z : in Partition;
                b : in out natural; q : out Poly_Sys;
                qsols : in out Solution_List );

-- DESCRIPTION :
--   This is an interactive driver for the construction of an interpolating
--   homotopy based on an m-homogeneous Bezout number.

-- ON ENTRY :
--   file       to write diagnostics on;
--   p          the polynomial system;
--   z          partition of the set of unknowns of p;
--   b          an m-homogeneous Bezout number.

-- ON RETURN :
--   b          number of interpolating roots;
--   q          an m-homogeneous start system;
--   qsols      solutions of q.
