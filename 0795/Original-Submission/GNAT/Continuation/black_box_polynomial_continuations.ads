with text_io;                       use text_io;
with Complex_Polynomial_Systems;    use Complex_Polynomial_Systems;
with Solutions;                     use Solutions;

package Black_Box_Polynomial_Continuations is

-- DESCRIPTION :
--   This package provides two procedure for performing polynomial
--   continuation in batch processing mode.  They mainly differ by the fact
--   that the homotopy might be already provided in the input parameter.

  procedure Black_Box_Polynomial_Continuation
                  ( infile,outfile : in file_type; pocotime : out duration );

  procedure Black_Box_Polynomial_Continuation
                  ( infile,outfile : in file_type; 
                    p,q : in Poly_Sys; sols : in out Solution_List;
                    pocotime : out duration );

  -- DESCRIPTION :
  --   Performs polynomial continuation with default settings.

  -- REQUIRED :
  --   Both files must be opened respectively for input and output.

  -- ON ENTRY :
  --   infile       the input file;
  --   outfile      file to write the results on;
  --   p            target polynomial system;
  --   q            a start system for solving p;
  --   sols         the solutions of q.

  -- ON RETURN :
  --   sols         the solutions of p;
  --   pocotime     elapsed user cpu time for polyhedral continuation.

end Black_Box_Polynomial_Continuations;
