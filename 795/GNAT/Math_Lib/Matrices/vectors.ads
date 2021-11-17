generic

  type coefftp is private;     -- type of the items in the vectors

  zero : coefftp;              -- the neutral element for "+"
                             
  with procedure clear    ( a : in out coefftp );
  with procedure copy     ( a : in coefftp; b : in out coefftp );
  with function  equal    ( a,b : coefftp ) return boolean;

  with function "+"       ( a,b : coefftp ) return coefftp;  -- return a+b;
  with function "-"       ( a,b : coefftp ) return coefftp;  -- return a-b;
  with function "-"       ( a : coefftp )   return coefftp;  -- return -a;
  with function "*"       ( a,b : coefftp ) return coefftp;  -- return a*b;

  with procedure Plus_Coeff ( a : in out coefftp; b : in coefftp ); -- a := a+b;
  with procedure Min_Coeff  ( a : in out coefftp; b : in coefftp ); -- a := a-b;
  with procedure Min_Coeff  ( a : in out coefftp );                 -- a := -a;
  with procedure Mult_Coeff ( a : in out coefftp; b : in coefftp ); -- a := a*b;

package Vectors is

-- DESCRIPTION :
--   This package provides a data abstraction for vectors with coefficients
--   over any ring.
--   The same functionality is provided for pointers to vectors.

  type Vector is array ( integer range <> ) of coefftp;
  type Link_to_Vector is access Vector;

  Range_Error : exception; -- occurs when two vectors have incompatible ranges

  procedure Clear ( v : in out Vector );
  procedure Clear ( v : in out Link_to_Vector );

  function  Equal ( v1,v2 : Vector ) return boolean;
  function  Equal ( v1,v2 : Link_to_Vector ) return boolean;
  procedure Copy  ( v1: in Vector; v2 : in out Vector );
  procedure Copy  ( v1: in Link_to_Vector; v2 : in out Link_to_Vector );

  function "+" ( v1,v2 : Vector ) return Vector;            -- return v1+v2;
  function "+" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
  function "-" ( v : Vector ) return Vector;                -- return -v;
  function "-" ( v : Link_to_Vector ) return Link_to_Vector;
  function "-" ( v1,v2 : Vector ) return Vector;            -- return v1-v2;
  function "-" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
  function "*" ( v : Vector; a : coefftp ) return Vector;   -- return v*a;
  function "*" ( v : Link_to_Vector; a : coefftp ) return Link_to_Vector; 
  function "*" ( a : coefftp; v : Vector ) return Vector;   -- return a*v;
  function "*" ( a : coefftp; v : Link_to_Vector ) return Link_to_Vector;

  function "*" ( v1,v2 : Vector ) return coefftp;           
  function "*" ( v1,v2 : Link_to_Vector ) return coefftp;
     -- returns the inner product of the vectors v1 and v2

  function "*" ( v1,v2 : Vector ) return Vector;
  function "*" ( v1,v2 : Link_to_Vector ) return Link_to_Vector;
     -- returns the vector v, with v(k) = v1(k)*v2(k);

  function Sum ( v : Vector ) return coefftp;
  function Sum ( v : Link_to_Vector ) return coefftp;
     -- returns the sum of all components of v;

  procedure Plus_Vector ( v1 : in out Vector; v2 : in Vector ); -- v1 := v1+v2;
  procedure Plus_Vector ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector ); 
  procedure Min_Vector  ( v  : in out Vector);                  -- v  := -v;
  procedure Min_Vector  ( v  : in out Link_to_Vector );
  procedure Min_Vector  ( v1 : in out Vector; v2 : in Vector ); -- v1 := v1-v2;
  procedure Min_Vector  ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector );
  procedure Mult_Coeff  ( v  : in out Vector; a : in coefftp ); -- v  := v*a;
  procedure Mult_Coeff  ( v  : in out Link_to_Vector; a : in coefftp );
  procedure Mult_Vector ( v1 : in out Vector; v2 : in Vector ); -- v1 := v1*v2;
  procedure Mult_Vector ( v1 : in out Link_to_Vector; v2 : in Link_to_Vector );

end Vectors;
