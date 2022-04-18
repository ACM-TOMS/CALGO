with text_io;       use text_io;
with Simplices;     use Simplices;

package Simplices_io is

-- DESCRIPTION :
--   This package provides an interface to the simplices.

  procedure get ( s : in out Simplex );
  procedure get ( n : in natural; s : in out Simplex );
  procedure get ( file : in file_type; s : in out Simplex );
  procedure get ( file : in file_type; n : in natural; s : in out Simplex );

   -- DESCRIPTION :
   --   Reads the dimension n if not specified as parameter,
   --   and then n integer vectors of length n, from standard input
   --   or from file.

  procedure put ( s : in Simplex );
  procedure put ( file : in file_type; s : in Simplex );

   -- DESCRIPTION :
   --   Writes the n vectors that span the simplex on standard output
   --   or on file.

end Simplices_io;
