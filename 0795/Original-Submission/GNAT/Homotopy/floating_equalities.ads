with Floating_Point_Numbers;         use Floating_Point_Numbers;
with Complex_Numbers;                use Complex_Numbers;
with Float_Vectors,Complex_Vectors;

package Floating_Equalities is

-- DESCRIPTION :
--   This package collects functions to test equality of floating-point
--   (complex) numbers and vectors upon a given tolerance.

  function Is_Equal ( x,y,tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(x-y) <= tol, otherwise false is returned.

  function Is_Equal ( x,y : double_complex; tol : double_float )
                    return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(REAL_PART(x)-REAL_PART(y)) <= tol and
  --   abs(IMAG_PART(x)-IMAG_PART(y)) <= tol, otherwise false is returned.

  function Is_Equal ( x,y : Float_Vectors.Vector; tol : double_float )
                    return boolean;

  -- DESCRIPTION :
  --   Returns true if Is_Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

  function Is_Equal ( x,y : Complex_Vectors.Vector; tol : double_float )
                    return boolean;

  -- DESCRIPTION :
  --   Returns true if Is_Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

end Floating_Equalities;
