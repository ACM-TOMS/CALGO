with text_io,Solutions;                use text_io,Solutions;
with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;  use Complex_Numbers,Complex_Vectors;

package Process_io is

-- DESCRIPTION :
--   This package determines the output operations during the continuation.
--   The purpose is to concentrate all trivial output operations which could
--   possibly overload the coding of the continuation process.
--   Moreover, an uniform output format is achieved by this package.

  type output_code is ( nil,s,p,c,sp,sc,pc,spc );

  -- Explanation of the output_code :
  --   nil  :  nothing will be written during the continuation process
  --   s    :  all intermediate solutions are written
  --   p    :  predictor information is written
  --   c    :  corrector information is written
  --   sp, sc, pc and spc are combinations of s, p and c

  procedure Set_Output_Code ( u : in output_code );
 
  -- DESCRIPTION :
  --   Sets the status code for output during continuation.

  procedure Write_path ( n : in positive );
  procedure Write_path ( ft : in file_type; n : in positive ); 

  -- DESCRIPTION :
  --   The number of the paths is written on file or on standard output.

  procedure Write_block ( n : in positive );
  procedure Write_block ( ft : in file_type; n : in positive );

  -- DESCRIPTION :
  --   The block number is written on the output device

  procedure sWrite ( sol : in Solution);
  procedure sWrite ( ft : in file_type; sol : in Solution );

  -- DESCRIPTION :
  --   The solution is written on file or on standard output.

  procedure pWrite ( step : in double_float; t : in double_complex );
  procedure pWrite ( ft : in file_type;
                     step : in double_float; t : in double_complex );
  procedure pWrite ( step : in double_float; t : in double_complex;
                     sol : in Solution );
  procedure pWrite ( ft : in file_type; step : in double_float;
                     t : in double_complex; sol : in Solution );
  -- DESCRIPTION :
  --   The predictor information is written on file or on standard output.

  procedure cWrite ( normax,normrx,normaf,normrf : in double_float );
  procedure cWrite ( ft : in file_type;
                     normax,normrx,normaf,normrf : in double_float );

  -- DESCRIPTION :
  --   The norm of the correction on x and residual is written.

  -- ON ENTRY :
  --   ft       file type, must be created or opened for output,
  --            if not specified, then standard output will be taken;
  --   normax   absolute norm of the correction dx on the solution x: ||dx||;
  --   normrx   relative norm of the correction dx: ||dx||/||x||;
  --   normaf   absolute norm of the residual: ||f(x)||;
  --   normrf   relative norm of the residual: ||f(x)||/||x||.

  procedure cWrite ( rcond : in double_float; m : in natural );
  procedure cWrite ( ft : in file_type;
                     rcond : in double_float; m : in natural );
  -- DESCRIPTION :
  --   The estimate for the inverse condition number of the Jacobi matrix
  --   is written, jointly with the (estimated) multiplicity of the solution.

  procedure Write_Statistics ( nstep,nfail,niter,nsyst : in natural );
  procedure Write_Statistics ( ft : in file_type;
                               nstep,nfail,niter,nsyst : in natural );
  -- DESCRIPTION :
  --   This procedure writes statistical information after the
  --   computation of parts of the results.

  -- ON ENTRY :
  --   nstep     the number of predictor steps;
  --   nfail     the number of failures;
  --   niter     the number of corrector iterations;
  --   nsyst     the number of linear systems solved.

  procedure Write_Total_Statistics ( tnstep,tnfail,tniter,tnsyst : in natural );
  procedure Write_Total_Statistics ( ft : in file_type;
                                     tnstep,tnfail,tniter,tnsyst : in natural );
  -- DESCRIPTION
  --   This procedure writes statistical information after the 
  --   solution of the problem.

  -- ON ENTRY :
  --   tnstep     the total number of predictor steps;
  --   tnfail     the total number of failures;
  --   tniter     the total number of corrector iterations;
  --   tnsyst     the total number of linear systems solved.

  procedure Write_convergence_factor ( factor : in double_float );
  procedure Write_convergence_factor 
                ( ft : in file_type; factor : in double_float );
  -- DESCRIPTION :
  --   writes the convergence factor of the correction process

  procedure sWrite_Solutions ( sols : in Solution_List );
  procedure sWrite_Solutions ( ft : in file_type; sols : in Solution_List );

  -- DESCRIPTION :
  --   Writes down the computed solutions on standard output or on file.

end Process_io;
