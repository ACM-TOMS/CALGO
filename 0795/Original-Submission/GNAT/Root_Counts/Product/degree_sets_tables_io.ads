with text_io;              use text_io;
with Degree_Sets_Tables;   use Degree_Sets_Tables;

package Degree_Sets_Tables_io is

-- DESCRIPTION :
--   This package provides some output facilities for degree set tables.

  procedure put ( ase : in Array_of_Sets );
  procedure put ( file : in file_type; ase : in Array_of_Sets );

  -- DESCRIPTION :
  --   Writes the array of sets on file or on standard output,

  procedure put ( dst : in Degree_Sets_Table );
  procedure put ( file : in file_type; dst : in Degree_Sets_Table );

  -- DESCRIPTION :
  --   Writes the degree set table on file or on standard output,
  --   if the file is omitted.

end Degree_Sets_Tables_io;
