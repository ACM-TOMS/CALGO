with text_io;                          use text_io;
with Complex_Numbers,Solutions;        use Complex_Numbers,Solutions;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;

procedure Driver_for_Homotopy_Construction
             ( file : in file_type; p,q : in out Poly_Sys;
               qsols : in out Solution_List; target : out double_complex );

-- DESCRIPTION :
--   This is an interactive driver for the construction of an artificial
--   parameter homotopy.  The user ask to homogenize the polynomial system
--   and can determine the homotopy parameters and the target value for the
--   continuation parameter.  This target value is returned as output.
--   Note that the starting value for the continuation parameter is stored
--   with the solutions.
