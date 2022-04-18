with text_io,Integer_Matrices;  use text_io,Integer_Matrices;

package Integer_Matrices_io is

-- DESCRIPTION :
--   This package contains routines for input and output
--   of integer matrices.

  procedure get ( m : out Matrix );
  procedure get ( file : in file_type; m : out Matrix );

  procedure get ( m : out Matrix; rw1,rw2 : in integer );
  procedure get ( file : in file_type; m : out Matrix; rw1,rw2 : in integer );

  -- DESCRIPTION :
  --   Reads an integer matrix m or m(rw1..rw2,m'range(2))
  --   from standard input or on file.

  procedure put ( m : in Matrix );
  procedure put ( file : in file_type; m : in Matrix );
  procedure put ( m : in Matrix; width : in natural );
  procedure put ( file : in file_type; m : in Matrix; width : in natural );

  -- DESCRIPTION :
  --   Writes an integer matrix m on file or on standard output.
  --   The width of the columns can be given as parameter. 
  --   By default, only one additional spacing is used to separate
  --   the columns.

  procedure put ( m : in Matrix; rw1,rw2 : in integer );
  procedure put ( file : in file_type; m : in Matrix; rw1,rw2 : in integer );

  -- DESCRIPTION :
  --   Writes an integer matrix m or m(rw1..rw2,m'range(2))
  --   on standard output or on file.

end Integer_Matrices_io;
