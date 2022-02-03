with text_io;                        use text_io;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;

package Drivers_for_Reduction is

-- DESCRIPTION:
--   This package collects driver routines for reducing a polynomial
--   system w.r.t. its total degree.

  procedure Display_Info;

  -- DESCRIPTION :
  --   Displays information about reduction on screen.

  procedure Driver_for_Linear_Reduction
               ( file : in file_type; p : in out Poly_Sys; d : out natural );

  -- DESCRIPTION :
  --   The coefficient matrix of the system is triangulated.

  procedure Driver_for_Sparse_Linear_Reduction
               ( file : in file_type; p : in out Poly_Sys; d : out natural );

  -- DESCRIPTION :
  --   The coefficient matrix of the system is diagonalized. 

  procedure Driver_for_Nonlinear_Reduction
               ( file : in file_type; p : in out Poly_Sys; d : out natural );

  -- DESCRIPTION :
  --   Combinations of S-polynomials are used to lower the total degree.

  procedure Driver_for_Overconstrained_Reduction
               ( file : in file_type; p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Random combinations are added.

  procedure Driver_for_Reduction 
               ( file : in file_type; p : in out Poly_Sys; d : out natural;
                 exit_option : in boolean );

  -- DESCRIPTION :
  --   This is an interactive driver for the reduction procedures.

  -- ON ENTRY :
  --   file         a file to write intermediate results and diagnostics on;
  --   p            a polynomial system;
  --   exit_option  if true, then the leave-option will be shown.

  -- ON RETURN :
  --   p            the system in a reduced form;
  --   d            total degree of the new system.

end Drivers_for_Reduction;
