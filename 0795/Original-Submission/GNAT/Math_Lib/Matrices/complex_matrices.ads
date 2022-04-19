with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;

package Complex_Matrices is

-- DESCRIPTION :
--   This package offers a data abstraction and operations for
--   working with matrices of floating point numbers.

  type Matrix is array ( integer range <>, integer range <> ) of double_complex;

-- MATRIX-MATRIX OPERATIONS :

  function "+" ( a,b : Matrix ) return Matrix;           -- return a+b
  function "-" ( a,b : Matrix ) return Matrix;           -- return a-b
  function "-" ( a   : Matrix ) return Matrix;           -- return -a
  function "*" ( a,b : Matrix ) return Matrix;           -- return a*b

  procedure Mult1 ( a : in out Matrix; b : in Matrix );  -- a := a*b
  procedure Mult2 ( a : in Matrix; b : in out Matrix );  -- b := a*b

-- MATRIX-VECTOR OPERATIONS :

  function "*" ( a : Matrix; v : Vector ) return Vector; -- return a*v
  procedure Mult ( a : in Matrix; v : in out Vector );   -- v := a*v

end Complex_Matrices;
