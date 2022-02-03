with Floating_Point_Numbers;   use Floating_Point_Numbers;

package Numbers_io is

-- DESCRIPTION :
--   This package provides user friendly input routines.

  procedure Read_Positive ( p : out positive );
  procedure Read_Natural ( n : out natural );
  procedure Read_Integer ( i : out integer );
  procedure Read_Single_Float ( f : out single_float );
  procedure Read_Double_Float ( f : out double_float );

  -- DESCRIPTION :
  --   Reads a number from standard input.
  --   As long as the value obtained is not of the right type,
  --   the user will be asked to try again.

end Numbers_io;
