with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Integer_Vectors,Float_Vectors;
with Integer_Matrices,Float_Matrices;

package Float_to_Integer_Rounding is

-- DESCRIPTION :
--   This package provides operations to round the floating point solutions
--   of linear inequality systems to vectors with integer entries.

-- BASIC CONVERSION OPERATIONS :

  function Round ( v : Float_Vectors.Vector ) return Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the integer vector with rounded entries.

  function Convert_to_Float ( v : Integer_Vectors.Vector )
                            return Float_Vectors.Vector;

  -- DESCRIPTION :
  --   Converts the given integer vector to a floating point vector.

  function Convert_to_Float ( m : Integer_Matrices.matrix )
                            return Float_Matrices.matrix;

  -- DESCRIPTION :
  --   Converts an integer matrix to a floating point matrix.

-- CONVERSIONS BY CONTINUED FRACTIONS :

  procedure Fractions ( x,tol : in double_float; a,b : out integer );

  -- DESCRIPTION :
  --   On return, we have that abs(x - a/b) < tol.

  function Scale_to_Integer ( v : Float_Vectors.Vector; tol : in double_float )
                            return Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Scales the given rational vector to a vector of integer numbers by
  --   multiplication by the common denominator.

  function Scale_to_Integer
               ( ine : Float_Matrices.Matrix; v : Float_Vectors.Vector;
                 tol : in double_float ) return Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the first integer vector that satisfies the inequality system,
  --   derived by continued fractions.

end Float_to_Integer_Rounding;
