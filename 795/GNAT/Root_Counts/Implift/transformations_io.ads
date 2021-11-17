with text_io;                use text_io;
with Transformations;        use Transformations;

package Transformations_io is

-- DESCRIPTION :
--   This package contains routines for input and output
--   of transformations.

  procedure get ( n : in natural; t : out Transfo );
  procedure get ( file : in file_type; n : in natural; t : out Transfo );

  -- DESCRIPTION :
  --   Reads n vectors from standard input or from file.
  --   These vectors are considered as the images of
  --   the basis vectors under the transformation t.

  procedure put ( t : in Transfo );
  procedure put ( file : in file_type; t : in Transfo );

  -- DESCRIPTION :
  --   Writes the images of the basis vectors under t.

end Transformations_io;
