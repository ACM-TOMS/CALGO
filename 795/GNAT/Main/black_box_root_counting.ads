with text_io;                           use text_io;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Solutions;                         use Solutions;

procedure Black_Box_Root_Counting 
               ( file : in out file_type;
                 p : in Poly_Sys; rc : out natural;
                 q : out Poly_Sys; qsols : out Solution_List;
                 rocotime,hocotime : out duration );

-- DESCRIPTION :
--   Calculates four different root counts: total degree, m-homogeneous
--   Bezout number, generalized Bezout number based on set structure,
--   and mixed volume.  Heuristics are used for the Bezout numbers.
--   Returns the start system with lowest root count and least amount
--   of work, which means that linear-product start systems are prefered,
--   when Bezout numbers equal the mixed volume.

-- ON ENTRY :
--   file        must be opened for output;
--   p           a polynomial system.

-- ON RETURN :
--   rc          root count, Bezout number or mixed volume;
--   q           start system;
--   qsols       solutions of q, Length_Of(qsols) = rc;
--   rocotime    elapsed user cpu time for computation of the root counts;
--   hocotime    elapsed user cpu time for construction of start system.
