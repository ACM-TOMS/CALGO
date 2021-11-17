with Integer_Vectors;   use Integer_Vectors;
with Transformations;   use Transformations;

package Integer_Vectors_Utilities is

-- DESCRIPTION :
--   This package offers utilities for working with
--   integer vectors.

  function Pivot ( v : Vector ) return integer;
  function Pivot ( v : Link_to_Vector ) return integer;

  -- DESCRIPTION :
  --   Returns the first nonzero entry out of v.
  --   If all entries of v are zero, then v'last+1 is returned.

  function gcd ( v : Vector ) return integer;
  function gcd ( v : Link_to_Vector ) return integer;

  -- DESCRIPTION :
  --   Returns the gcd of all entries in v.

  procedure Normalize ( v : in out Vector );
  procedure Normalize ( v : in out Link_to_Vector );

  -- DESCRIPTION :
  --   All elements in d are divided by the gcd of the array.

  function  Reduce ( v : Vector; i : integer ) return Vector;
  procedure Reduce ( v : in out Link_to_Vector; i : in integer );
  function  Reduce ( v : Link_to_Vector; i : integer ) return Link_to_Vector;

  -- DESCRIPTION :
  --   The i-th component will be deleted out of the vector.

  function  Insert ( v : Vector; i,a : integer ) return Vector;
  procedure Insert ( v : in out Link_to_Vector; i,a : in integer );
  function  Insert ( v : Link_to_Vector; i,a : integer ) return Link_to_Vector;

  -- DESCRIPTION :
  --   The i-th component will be inserted, using the value a.

  function  Insert_and_Transform
             ( v : Vector; i,a : integer; t : Transfo ) return Vector;
  procedure Insert_and_Transform
             ( v : in out Link_to_Vector; i,a : in integer; t : in Transfo );
  function  Insert_and_Transform
             ( v : Link_to_Vector; i,a : integer; t : Transfo )
             return Link_to_Vector;

  -- DESCRIPTION :
  --   Inserts the i-th component in the vector v,
  --   using the value a, and transforms the vector, applying t.

end Integer_Vectors_Utilities;
