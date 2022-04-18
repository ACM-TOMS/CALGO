with Floating_Point_Numbers;    use Floating_Point_Numbers;

package Mathematical_Functions is

-- DESCRIPTION :
--   This package provides some special mathematical functions to ensure
--   a more portable version of the software.

-- CONSTANT :

  PI : constant :=  3.14159_26535_89793_23846_26433_83279_50288;

-- EXPONENTIAL AND LOGARITHMIC FUNCTIONS :

  function "**" ( x,y : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns x**y.

  function LOG10 ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns the usual logarithm.

  function SQRT ( x : double_float ) return double_float;

  -- DSECRIPTION :
  --   Returns the square root of x.

-- TRIGONIOMETRIC FUNCTIONS :

  function SIN ( x : double_float ) return double_float;
  function COS ( x : double_float ) return double_float;
  function TAN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns sine, cosine and tangens of x.

  function ARCSIN ( x : double_float ) return double_float;
  function ARCCOS ( x : double_float ) return double_float;
  function ARCTAN ( x : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns arcsin, arccos and argtan of x.

  function Radius ( x,y : double_float ) return double_float;
  function Angle  ( x,y : double_float ) return double_float;

end Mathematical_Functions;
