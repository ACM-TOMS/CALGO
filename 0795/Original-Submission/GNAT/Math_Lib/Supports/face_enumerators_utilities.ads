with Integer_Vectors;
with Integer_Vectors_of_Vectors;

package Face_Enumerators_Utilities is

-- DESCRIPTION :
--   This package contains utilities for the face enumerators.

  function Is_Zero ( v : Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the given vector equals the zero vector.

  function gcd ( v : Integer_Vectors.Vector ) return integer;

  -- DESCRIPTION :
  --   Returns the greatest common divisor gcd(v(v'first),..,v(v'last)).

  procedure Scale ( v : in out Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides each element in v by gcd(v), when gcd(v) /= 0 of course.

  function Is_In ( x : integer; v : Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists an entry in v, say v(k),
  --   such that v(k) = x.

  function In_Edge ( x,a,b : Integer_Vectors.Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector x lies between a and b.

  function In_Face ( k : in natural; face,x : Integer_Vectors.Vector;
                     pts : Integer_Vectors_of_Vectors.Vector )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if x lies in the given k-face, which contains entries
  --   to its elements in pts.

end Face_Enumerators_Utilities;
