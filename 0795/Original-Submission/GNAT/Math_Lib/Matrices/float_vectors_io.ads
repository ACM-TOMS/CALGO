with text_io,Float_Vectors;    use text_io,Float_Vectors;

package Float_Vectors_io is

-- DESCRIPTION :
--   This package contains routines for input and output
--   of vectors with float entries.

  procedure get ( v : in out Vector );
  procedure get ( file : in file_type; v : in out Vector );

  -- DESCRIPTION :
  --   float number will be read from standard input or from file,
  --   until all entries of v are filled.

  procedure get ( n : in natural; v : out Link_to_Vector );
  procedure get ( file : in file_type; n : in natural;
                  v : out Link_to_Vector );

  -- DESCRIPTION :
  --   n float numbers are read from standard input or from file.

  procedure put ( v : in Vector );
  procedure put ( file : in file_type; v : in Vector );
  procedure put ( v : in Link_to_Vector );
  procedure put ( file : in file_type; v : in Link_to_Vector );

  procedure put ( v : in Vector; fore,aft,exp : in natural );
  procedure put ( file : in file_type; v : in Vector;
                  fore,aft,exp : in natural );
  procedure put ( v : in Link_to_Vector; fore,aft,exp : in natural );
  procedure put ( file : in file_type; v : in Link_to_Vector;
                  fore,aft,exp : in natural );

  -- DESCRIPTION :
  --   the vector v is written on standard output or on file.

end Float_Vectors_io;
