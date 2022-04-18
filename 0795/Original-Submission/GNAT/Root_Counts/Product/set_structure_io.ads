with text_io;       use text_io;

package Set_Structure_io is

-- DESCRIPTION :
--   This package contains i/o operations of a set structure.

  procedure get;
  procedure get ( file : in file_type );

  -- DESCRIPTION :
  --   The set structure will be read from standard input or from file.
  --   Example: { x1 } { x1 x2 } { x2 }
  --    --> separation of unknowns and delimiters by spaces. 

  -- REQUIRED :
  --   The Symbol_Table must be initialized.

  procedure put;
  procedure put ( file : in file_type );

  -- DESCRIPTION :
  --   The set structure is written on standard output or on file.

end Set_Structure_io;
