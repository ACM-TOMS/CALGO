with text_io,Complex_Matrices;  use text_io,Complex_Matrices;

package Complex_Matrices_io is

-- DESCRIPTION :
--   This package contains routines for input and output
--   of matrices with complex entries.

  procedure get ( m : out Matrix );
  procedure get ( file : in file_type; m : out Matrix );

  -- DESCRIPTION :
  --   reads a matrix from standard input or on file.

  procedure put ( m : in Matrix );
  procedure put ( file : in file_type; m : in Matrix );
  procedure put ( m : in Matrix; fore,aft,exp : in natural );
  procedure put ( file : in file_type; m : in Matrix;
                  fore,aft,exp : in natural );

  -- DESCRIPTION :
  --   writes a matrix on standard output or on file.

end Complex_Matrices_io;
