with Integer_Vectors;            use Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;   use Lists_of_Integer_Vectors;

package Lists_of_Vectors_Utilities is

-- DESCRIPTION :
--   This package offers some utilities for working with
--   lists of integer vectors.

  procedure Compute_Normal ( v : in Integer_Vectors_of_Vectors.Vector; 
                             n : out Link_to_Vector; deg : out natural );
  function  Compute_Normal ( v : Integer_Vectors_of_Vectors.Vector )
                           return Link_to_Vector;
  -- DESCRIPTION :
  --   Returns the normal vector to the space generated by the
  --   points in v.  If deg = 0, then more than one solution is possible.

  function Pointer_to_Last ( l : List ) return List;

  -- DESCRIPTION :
  --   Returns a pointer to the last element of the list.

  procedure Move_to_Front ( l : in out List; v : in Vector );

  -- DESCRIPTION :
  --   Searches the vector v in the list l.  When found, then this
  --   vector v is swapped with the first element of the list.

  function Difference ( l1,l2 : List ) return List;

  -- DESCRIPTION :
  --   Returns the list of points in l1 that do not belong to l2.

  function Different_Points ( l : List ) return List;

  -- DESCRIPTION :
  --   Returns a lists of all different points out of l.

  procedure Remove_Duplicates ( l : in out List );

  -- DESCRIPTION :
  --   Removes duplicate points out of the list l.

end Lists_of_Vectors_Utilities;