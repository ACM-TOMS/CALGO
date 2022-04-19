with text_io;         use text_io;
with Solutions;       use Solutions;

package Solutions_io is

-- DESCRIPTION :
--   This routines provides routines for input and output of solutions.

-- FOR SOLUTION VECTORS ONLY :

  procedure get_vector ( s : in out Solution );
  procedure get_vector ( file : in file_type; s : in out Solution );

  -- DESCRIPTION :
  --   The input must contain the solution vector.

  procedure put_vector ( s : in Solution );
  procedure put_vector ( file : in file_type; s : in Solution );

  -- DESCRIPTION :
  --   On the output the solution vector will be written.

-- FOR SOLUTIONS :

  procedure get ( s : in out Solution );
  procedure get ( file : in file_type; s : in out Solution );

  -- DESCRIPTION :
  --   The input must contain the following : s.t, s.m and s.v(i), 
  --   a vector of s.n complex numbers

  procedure put ( s : in Solution );
  procedure put ( file : in file_type; s : in Solution );

  -- DESCRIPTION :
  --   On the output the following will be written :
  --   s.t, s.m and s.v, a vector of s.n complex numbers

-- FOR LISTS OF SOLUTIONS :

  procedure get ( sols : in out Solution_List );
  procedure get ( sols,sols_last : in out Solution_List );
  procedure get ( len,n : in natural; sols : in out Solution_List );
  procedure get ( len,n : in natural; sols,sols_last : in out Solution_List );
  procedure get ( file : in file_type; sols : in out Solution_List );
  procedure get ( file : in file_type; sols,sols_last : in out Solution_List );
  procedure get ( file : in file_type; len,n : in natural;
                  sols : in out Solution_List );
  procedure get ( file : in file_type; len,n : in natural;
                  sols,sols_last : in out Solution_List );

  -- DESCRIPTION :
  --   A solution list will be read.  If the length len and dimension n
  --   of the list is not supplied, then they will be read first.
  --   If the parameter sols_last is supplied, then this parameter contains
  --   the pointer to the last element of the list on return.
  --   The solutions should be in the appropriate format.

  procedure put ( sols : in Solution_List );
  procedure put ( len,n : in natural; sols : in Solution_List );
  procedure put ( file : in file_type; sols : in Solution_List );
  procedure put ( file : in file_type; len,n : in natural;
                  sols : in Solution_List );

  -- DESCRIPTION :
  --   The solutions are written on standard output or on file.
  --   First the length of the list and the dimension of the solutions
  --   will be put on file if they are supplied as parameter.

-- USER-FRIENDLY ROUTINES :

  procedure Display_Format;

  -- DESCRIPTION :
  --   Displays on screen the formatting rules as on-line help facility.

  procedure Read ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   Reads the solution list from file, displays the formatting information
  --   in case of exception and let the user try again.

end Solutions_io;
