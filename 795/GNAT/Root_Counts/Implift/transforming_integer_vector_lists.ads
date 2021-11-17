with Integer_Vectors;            use Integer_Vectors;
with Lists_of_Integer_Vectors;   use Lists_of_Integer_Vectors;
with Transformations;            use Transformations;

package Transforming_Integer_Vector_Lists is

-- DESCRIPTION :
--   This package provides unimodular transformations of
--   lists of integer vectors.

  procedure Shift ( l : in out List; v : in Vector );
  procedure Shift ( l : in out List; v : in Link_to_Vector );

  function  Shift ( l : List; v : Vector ) return List;
  function  Shift ( l : List; v : Link_to_Vector ) return List;

  -- DESCRIPTION :
  --   The list will be shifted: Shift(l,v) = { y-v | Is_In(l,y) }

  function "*"( l : List; t : Transfo ) return List;
  function "*"( t : Transfo; l : List ) return List;

  -- DESCRIPTION :
  --   Returns the transformed list of points.

  procedure Apply ( l : in out List; t : in Transfo );

  -- DESCRIPTION :
  --   Applies the transformation t to the list l.

  function  Reduce ( l : List; i : integer ) return List;
  procedure Reduce ( l : in out List; i : in integer );

  -- DESCRIPTION :
  --   Returns a list of vectors where the i-th component has been deleted.

  function  Insert ( l : List; i,a : integer ) return List;
  procedure Insert ( l : in out List; i,a : in integer );

  -- DESCRIPTION :
  --   Returns a list of vectors where the i-th component has been inserted,
  --   for all d in l: d(i) = a.

  function  Transform_and_Reduce ( t : Transfo; i : integer; l : List )
                                 return List;
  procedure Transform_and_Reduce ( t : in Transfo; i : in integer;
                                   l : in out List );

  -- DESCRIPTION :
  --   Transforms the list l and deletes the i-th component
  --   of every element in the transformed list.

  function  Insert_and_Transform
             ( l : List; i,a : integer; t : Transfo ) return List;
  procedure Insert_and_Transform
             ( l : in out List; i,a : in integer; t : in Transfo );

  -- DESCRIPTION :
  --   Inserts the i-th component of every element in the list l,
  --   using the value a, and transforms the list, applying t.

end Transforming_Integer_Vector_Lists;
