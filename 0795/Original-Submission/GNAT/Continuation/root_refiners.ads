with text_io;                        use text_io;
with Floating_Point_Numbers;         use Floating_Point_Numbers;
with Complex_Vectors,Solutions;      use Complex_Vectors,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;
with Jacobi_Matrices;                use Jacobi_Matrices;

package Root_Refiners is

-- DESCRIPTION:
--   The aim of this package is to provide some routines to
--   perform root refining, to validate the approximate solutions.
--   It can be used as a postprocessor to check the computed roots,
--   or, as preprocessor, to find some suitable starting values for
--   the continuation procedure.

--   The basic root refiner is Newton's method.  There are two versions
--   provided : a silent and a reporting version.  The reporting version puts
--   intermediate results on file during the iterations, while the silent
--   version simply returns the refined solutions.

-- NEWTON'S METHOD :

  procedure Silent_Newton 
               ( p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jacobi;
                 zero : in out Solution; epsxa,epsfa : in double_float; 
                 numit : in out natural; max : in natural; fail : out boolean );

  procedure Reporting_Newton
               ( file : in file_type;
                 p_eval : in Eval_Poly_Sys; j_eval : in Eval_Jacobi;
                 zero : in out Solution; epsxa,epsfa : in double_float;
                 numit : in out natural; max : in natural; fail : out boolean );

  -- DESCRIPTION :
  --   Newton's method is applied to refine the approximation of a root.
  --   The stopping criteria are:
  --     * numit > max   (maximum number of iterations is exceeded);
  --     * zero.err < epsxa (accuracy for x is reached);
  --     * zero.res < epsfa (tolerance for residual is reached).
  --   When one of these conditions is fulfilled, the procedure stops.

  -- ON ENTRY :
  --   file      to write intermediate results on;
  --   p_eval    evaluable form of the polynomial system;
  --   j_eval    evaluable form of the Jacobian matrix;
  --   zero      starting value;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value of the residue;
  --   numit     number of iterations, to be initialized with zero;
  --   max       maximum number of iterations.

  -- ON RETURN :
  --   zero      refined root;
  --   numit     number of iterations performed;
  --   fail      is true when the desired precision is not reached.

-- APPLICATION of Newton's Method on a list of solutions.
--   The silent versions simply perform the calculations.  
--   The reporting root refiners allow the output of intermediate results and
--   produce a summary of the calculations.
--   With each version, an additional output parameter can be supplied to
--   contain only those solutions that satisfy the accuracy requirements.
  
  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural );

  procedure Silent_Root_Refiner
               ( p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural );

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural; wout : in boolean );

  procedure Reporting_Root_Refiner
               ( file : in file_type;
                 p : in Poly_Sys; sols,refsols : in out Solution_List;
                 epsxa,epsfa,tolsing : in double_float;
                 numit : in out natural; max : in natural; wout : in boolean );

  -- DESCRIPTION :
  --   The list of solutions sols is refined w.r.t. the system p.
  --   The multiplicity of each solution in sols is determined as follows:
  --     m = 0 : if the solution is singular and probably non isolated
  --             or if the solution lies at infinity ( in fact no solution );
  --     m = 1 : if the solution is regular;
  --     m > 1 : a multiple solution with multiplicity m.

  -- ON ENTRY :
  --   file      file for writing diagnostics on;
  --   p         a polynomial system;
  --   sols      the start solutions;
  --   epsxa     maximum absolute error on the zero;
  --   epsfa     maximum absolute value for the residue;
  --   tolsing   tolerance on inverse condition number for singular solution;
  --   numit     the number of iterations, to be initialized with zero;
  --   max       maximum number of iterations per zero;
  --   wout      has to be true when intermediate output is wanted.

  -- ON RETURN :
  --   sols      a list of computed solutions;
  --   refsols   only those solutions which satisfy the given accuracy;
  --   numit     the number of iterations.

end Root_Refiners;
