with Integer_Vectors;            use Integer_Vectors;
with Lists_of_Integer_Vectors;   use Lists_of_Integer_Vectors;

package Integer_Support_Functions is

-- DESCRIPTION :
--   This package provides support functions for polytopes spanned by
--   lists of integer vectors.

  function Maximal_Support ( l : List; v : Vector ) return integer;
  function Minimal_Support ( l : List; v : Vector ) return integer;

  -- DESCRIPTION :
  --   Returns respectively max <d,v> or min <d,v>, for all d in l.

  procedure Min_Max ( l : in List; k : in integer; min,max : in out integer );

  -- DESCRIPTION :
  --   Computes the minimum and maximum of the k-th entry of all
  --   points belonging to l.

  function Graded_Max ( l : List ) return Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the maximal element in the list w.r.t. the graded
  --   lexicographical ordening.

  -- REQUIRED : List is not empty.

  function Face ( l : List; v : Vector; m : integer ) return List;

  -- DESCRIPTION :
  --   Returns a list of vectors d for which <d,v> = m.
  --   For m = Maximal_Support(l,v), we get the face w.r.t. the outer normal v;
  --   for m = Minimal_Support(l,v), we get the face w.r.t. the inner normal v.

  function Inner_Face ( l : List; v : Vector ) return List;
  function Outer_Face ( l : List; v : Vector ) return List;

  -- DESCRIPTION :
  --   Returns the face of the list l, where v is inner or outer normal.

end Integer_Support_Functions;
