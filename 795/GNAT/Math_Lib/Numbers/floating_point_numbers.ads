with text_io;

package Floating_Point_Numbers is

-- DESCRIPTION :
--   This package sets floating-point types to be independent
--   of the compiler predefined floating-point declarations.

  type single_float is digits 7;                  -- single precision
  package single_float_io is new text_io.float_io(single_float);

  type double_float is digits 15;                 -- double precision
  package double_float_io is new text_io.float_io(double_float);

end Floating_Point_Numbers;
