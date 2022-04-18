with text_io;                            use text_io;
with Float_Vectors;
with Float_Vectors_of_Vectors;
with Complex_Numbers,Solutions;          use Complex_Numbers,Solutions;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;

package Drivers_for_Polynomial_Continuation is

-- DESCRIPTION :
--   This package contains three drivers for two types of homotopies:
--   artificial and natural parameter.

  procedure Driver_for_Process_io ( file : in file_type; oc : out natural );

  -- DESCRIPTION :
  --   Choice of kind of output information during continuation.

  -- ON ENTRY :
  --   file     must be opened for output.

  -- ON RETURN :
  --   oc       number between 0 and 8 indicating the output code:
  --              0 : no intermediate output information during continuation;
  --              1 : only the final solutions at the end of the paths;
  --              2 : intermediate solutions at each step along the paths;
  --              3 : information of the predictor: t and step length;
  --              4 : information of the corrector: corrections and residuals;
  --              5 : intermediate solutions and information of the predictor;
  --              6 : intermediate solutions and information of the corrector;
  --              7 : information of predictor and corrector;
  --              8 : intermediate solutions, info of predictor and corrector.

  procedure Driver_for_Continuation_Parameters ( file : in file_type );

  -- DESCRIPTION :
  --   This procedure allows the user to determine all relevant parameters
  --   for the continuation.

  procedure Check_Continuation_Parameter ( sols : in out Solution_List );

  -- DESCRIPTION ;
  --   Reads the value of the continuation parameter for the first
  --   solution.  If different from zero, the user is given the
  --   opportunity to change it.

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; p : in Poly_Sys;
                  sols : out Solution_list; target : out double_complex );

  -- DESCRIPTION :
  --   This is a driver for the polynomial continuation routine
  --   with an artificial parameter homotopy.
  --   It reads the start system and start solutions and enables the
  --   user to determine all relevant parameters.
  
  -- ON ENTRY :
  --   file       to write diagnostics and results on;
  --   p          a polynomial system.

  -- ON RETURN :
  --   sols       the computed solutions.

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; p : in Poly_Sys; k : in natural;
                  target : in double_complex; sols : out Solution_list );

  -- DESCRIPTION :
  --   This is a driver for the polynomial continuation routine
  --   with a natural parameter homotopy.
  --   The start solutions will be read from file.
  --   A gentle interface makes it possible for the user to determine
  --   all relevant parameters.

  -- ON ENTRY :
  --   file       to write diagnostics and results on;
  --   p          a polynomial system, with n equations and n+1 unknowns;
  --   k          index of t = xk;
  --   target     target value for the continuation parameter.

  -- ON RETURN :
  --   sols       the computed solutions.

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; sols : in out Solution_List;
                  proj : in boolean;
                  target : double_complex := CMPLX(1.0) );

  -- DESCRIPTION :
  --   Given a homotopy, contained in the package Homotopy,
  --   the continuation procedure will be be carried out.
  --   The user may tune all continuation paramters.
   
  -- ON ENTRY :
  --   file       to write intermediate results and diagnostics on;
  --   sols       start solutions for the continuation;
  --   proj       true when a projective-perpendicular corrector will be used;
  --   target     target value for the continuation parameter.

  -- ON RETURN :
  --   sols       the computed solutions.

end Drivers_for_Polynomial_Continuation;
