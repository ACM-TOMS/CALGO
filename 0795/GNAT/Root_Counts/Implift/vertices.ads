with Integer_Vectors;            use Integer_Vectors;
with Lists_of_Integer_Vectors;   use Lists_of_Integer_Vectors;
with Float_Vectors;

package Vertices is

-- DESCRIPTION :
--   This function offers some functions for computing the 
--   vertices that span the polytope, given its set of points.
--   Based on the functions in this package it is possible to
--   determine the dimension of conv(l), without the actual computation
--   of the convex hull.

  function Is_In_Hull ( point : Vector; l : List ) return boolean;
  function Is_In_Hull ( point : Float_Vectors.Vector; l : List )
                      return boolean;

  -- DESCRIPTION :
  --   This function determines whether the point belongs to convex
  --   hull of the points in the given list.

  function Vertex_Points ( l : List ) return List;

  -- DESCRIPTION :
  --   The vertex points of conv(l) are returned.

  function Extremal_Points ( l : List; v : Link_to_Vector ) return List;

  -- DESCRIPTION :
  --   Returns the vertex points of the list l w.r.t. the direction v and -v.

  function Extremal_Points ( k,n : natural; exl,l : List ) return List;
  function Max_Extremal_Points ( k,n : natural; exl,l : List ) return List;

  -- DESCRIPTION :
  --   There are already k linearly independent vertex points found,
  --   given in exl.  This routine tries to take one more linearly
  --   independent vertex point out of the list l.
  --   The first function stops when a degeneracy is discovered,
  --   while the second one still proceeds and returns one point
  --   more, when possible.

  function Extremal_Points ( n : natural; l : List ) return List;

  -- DESCRIPTION :
  --   Searches in the list l for linearly independent vertex points.
  --   If dim(conv(l)) = n then n+1 points will be returned,
  --   otherwise the number of points in the resulting list will be
  --   less but not necessarily equal to dim(conv(l))+1.
  --   The advantage of this routine lies in the fact that a degeneracy,
  --   i.e., dim(conv(l)) < n, is detected as soon as possible.

  function Max_Extremal_Points ( n : natural; l : list ) return List;

  -- DESCRIPTION :
  --   Does the same as the function above, except for the fact that
  --   the number of points in the resulting lists equals dim(conv(l))+1.
  --   The points in the resulting list can be regarded as a basis,
  --   i.e., origin and as many linearly independent points as dim(conv(l)),
  --   for the points in l.

end Vertices;
