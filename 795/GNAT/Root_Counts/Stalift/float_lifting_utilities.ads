with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Integer_Vectors;
with Float_Vectors;                     use Float_Vectors;
with Lists_of_Float_Vectors;            use Lists_of_Float_Vectors;
with Float_Mixed_Subdivisions;          use Float_Mixed_Subdivisions;
with Arrays_of_Float_Vector_Lists;      use Arrays_of_Float_Vector_Lists;

package Float_Lifting_Utilities is

-- DESCRIPTION :
--   This package provides some utilities for dealing with lifting functions.

  function Adaptive_Lifting ( l : Array_of_Lists ) return Vector;

  -- DESCRIPTION :
  --   Returns upper bounds for a random lifting, depending on the lengths
  --   of the lists in l.

  procedure Search_Lifting ( l : in List; pt : in Vector;
                             found : out boolean; lif : out double_float );

  -- DESCRIPTION :
  --   Searches the lifting of the point in the lifted list l.
  --   If found, then lif equals the lifting, otherwise lif is meaningless.

  function Search_and_Lift ( l : List; pt : Vector ) return Vector;

  -- DESCRIPTION :
  --   Given a lifted list of points and a unlifted vector, the function
  --   either returns the corresponding lifted vector from the list, or
  --   the same point, when there is no lifted point in l whose projection
  --   equals the given point pt.

  function Induced_Lifting ( mixsub : Mixed_Subdivision; k : natural;
                             pt : Vector ) return Vector;
  function Induced_Lifting
               ( n : natural; mix : Integer_Vectors.Vector;
                 points : Array_of_Lists; mixsub : Mixed_Subdivision )
               return Array_of_Lists;

  -- DESCRIPTION :
  --   Given a mixed subdivision for a tuple of supports,
  --   then the lifted points will be returned as induced by the
  --   subdivision.   When points do not occur in the mixed subdivision,
  --   they will be lifted conservatively.

  function Conservative_Lifting 
               ( mic : Mixed_Cell; k : natural; point : Vector )
               return double_float;
  function Conservative_Lifting 
               ( mixsub : Mixed_Subdivision; k : natural; point : Vector )
               return double_float;
  
  -- DESCRIPTION :
  --   Returns the value of the conservative lifting function of the point
  --   to be considered w.r.t. the kth polytope.

  -- REQUIRED :
  --   The given point must already be in the lifted space and its last
  --   coordinate must contain already a lower bound for the lifting value.

end Float_Lifting_Utilities;
