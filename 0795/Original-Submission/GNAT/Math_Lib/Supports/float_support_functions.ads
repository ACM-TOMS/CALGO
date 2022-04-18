with Floating_Point_Numbers;      use Floating_Point_Numbers;
with Float_Vectors;               use Float_Vectors;
with Lists_of_Float_Vectors;      use Lists_of_Float_Vectors;

package Float_Support_Functions is

-- DESCRIPTION :
--   This package provides support functions for polytopes spanned by
--   lists of floating-point vectors.

  function Minimal_Support ( l : List; v : Vector ) return double_float;
  function Maximal_Support ( l : List; v : Vector ) return double_float;

  -- DESCRIPTION : 
  --   Returns respectively min <d,v> or max <d,v>, for all d in l.

  function Face ( l : List; v : Vector; m,tol : double_float ) return List;

  -- DESCRIPTION :
  --   Returns a list of vectors d for which abs(<d,v> - m) < tol,
  --   with tol the tolerance on precision.

end Float_Support_Functions;
