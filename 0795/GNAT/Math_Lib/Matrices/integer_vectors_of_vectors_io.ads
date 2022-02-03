with text_io;                          use text_io;
with Integer_Vectors_of_Vectors;       use Integer_Vectors_of_Vectors;

package Integer_Vectors_of_Vectors_io is

 -- DESCRIPTION :
 --   This package contains routines for input and output
 --   of vectors of vectors with integer entries.

  procedure get ( m : in natural; v : in out Vector );
  procedure get ( file : in file_type; m : in natural; v : in out Vector );

   -- DESCRIPTION :
   --   Integer numbers will be read from standard input or from file,
   --   until all entries of v are filled.
   --   The entries of v will be vectors of dimension m.

  procedure get ( n,m : in natural; v : out Link_to_Vector );
  procedure get ( file : in file_type; n,m : in natural;
                  v : out Link_to_Vector );

   -- DESCRIPTION :
   --   n integer vectors of dimension m are read from standard input
   --   or from file.

  procedure put ( v : in Vector );
  procedure put ( file : in file_type; v : in Vector );
  procedure put ( v : in Link_to_Vector );
  procedure put ( file : in file_type; v : in Link_to_Vector );

   -- DESCRIPTION :
   --   The vector v is written on standard output or on file.
   --   For each new element of v, a new line begins.

end Integer_Vectors_of_Vectors_io;
