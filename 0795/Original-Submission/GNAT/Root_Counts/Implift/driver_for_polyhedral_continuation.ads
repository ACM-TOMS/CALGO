with text_io;                         use text_io;
with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;

procedure Driver_for_Polyhedral_Continuation
               ( file : in file_type; p : in Poly_Sys; k : in natural;
                 byebye : in boolean;
                 q : out Poly_Sys; qfile,solsfile : in out file_type;
                 tosolve,ranstart,contrep : out boolean );

-- DESCRIPTION :
--   Prepares the settings for polyhedral continuation.

-- ON ENTRY :
--   file        output file, must be opened for output;
--   p           a polynomial system;
--   k           if k > 0, then q will not be written on file;
--   byebye      if true, then a bye-bye message will appear on screen,
--               if false, then no bye-bye.

-- ON RETURN :
--   q           system to solve, if tosolve=true;
--   qfile       output file for q, is opened for output;
--   solsfile    output file for solutions of q, is opened for output;
--   tosolve     true if polyhedral will be applied, false otherwise;
--   ranstart    true if q is random coefficient start system, false otherwise;
--   contrep     true if output information is wanted during continuation.
