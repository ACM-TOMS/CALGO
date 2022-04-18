with text_io,Integer_Vectors;             use text_io,Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;      use Arrays_of_Integer_Vector_Lists;

package Arrays_of_Integer_Vector_Lists_io is

-- DESCRIPTION :
--   This package offers routines for Input/Output of
--   arrays of lists of integer vectors.

  procedure get ( al : in out Link_to_Array_of_Lists );
  procedure get ( n : in natural; m : in Vector; al : out Array_of_Lists );
  procedure get ( file : in file_type; al : in out Link_to_Array_of_Lists );
  procedure get ( file : in file_type;
                  n : in natural; m : in Vector; al : out Array_of_Lists );

  -- DESCRIPTION :
  --   Reads a number of lists from standard input or from file;
  --   the ith list must contain m(i) integer vectors of length n.
  --   If n and m are not provided, then they are first read.

  procedure put ( al : in Array_of_Lists );
  procedure put ( file : in file_type; al : in Array_of_Lists );

  -- DESCRIPTION :
  --   Writes the lists in al on standard output or on file.

end Arrays_of_Integer_Vector_Lists_io;
