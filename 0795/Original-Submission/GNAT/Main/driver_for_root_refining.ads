with text_io,Complex_Vectors;              use text_io,Complex_Vectors;
with Complex_Polynomial_Systems;           use Complex_Polynomial_Systems;
with Solutions;                            use Solutions;

procedure Driver_for_Root_Refining
             ( file : in file_type; scalp,p : in Poly_Sys; basis : in natural;
               scalvec : in Link_to_Vector; sols : in out Solution_List);

-- DESCRIPTION :
--   This is the driver for root refining after the continuation.

-- ON ENTRY :
--   file        to write diagnostics on;
--   scalp       the scaled system;
--   p           the original polynomial system;
--   basis       used for scaling;
--   scalvec     vector of coefficients used for scaling;
--   sols        a list of solutions of scalp.

-- ON RETURN :
--   sols        the corrected solutions.
