with text_io;                    use text_io;
with Lists_of_Float_Vectors;     use Lists_of_Float_Vectors;

package Lists_of_Float_Vectors_io is

-- DESCRIPTION :
--   This package offers routines for Input/Output of
--   lists of integer vectors.

  procedure get ( n,m : in natural; l : out List );
  procedure get ( file : in file_type; n,m : in natural; l : out List );

  -- DESCRIPTION :
  --   Reads m integer vectors of length n from standard output or from file.

  procedure put ( l : in List );
  procedure put ( file : in file_type; l : in List );
  procedure put ( l : in List; fore,aft,exp : in natural );
  procedure put ( file : in file_type;
                  l : in List; fore,aft,exp : in natural );

  -- DESCRIPTION :
  --   Writes the vectors in l on standard output or on file.
  --   The fore, aft and exp determine the output format of the floats.

end Lists_of_Float_Vectors_io;
