with text_io,Complex_Numbers;   use text_io,Complex_Numbers;

package Complex_Numbers_io is

-- DESCRIPTION :
--   This package provides io-routines for complex numbers.

  procedure get ( c : out double_complex );
  procedure get ( f : in file_type; c : out double_complex );

  procedure put ( c : in double_complex );
  procedure put ( f : in file_type; c : in double_complex );

  procedure put ( c : in double_complex; fore,aft,exp : in natural );
  procedure put ( f : in file_type;
                  c : in double_complex; fore,aft,exp : in natural );

end Complex_Numbers_io;
