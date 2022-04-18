with Integer_Vectors,Float_Vectors;
with Lists_of_Integer_Vectors;
with Lists_of_Float_Vectors;
with Arrays_of_Integer_Vector_Lists;
with Arrays_of_Float_Vector_Lists;
with Integer_Mixed_Subdivisions;
with Float_Mixed_Subdivisions;

package Float_Integer_Convertors is

-- DESCRIPTION :
--   This package provides routines to convert lists of integer and floating-
--   point vectors into lists of floating-point and integer vectors.
--   The conversion from float to integer is done by merely rounding.

  function Convert ( v : Integer_Vectors.Vector ) return Float_Vectors.Vector;
  function Convert ( v : Float_Vectors.Vector ) return Integer_Vectors.Vector;

  function Convert ( l : Lists_of_Integer_Vectors.List )
                   return Lists_of_Float_Vectors.List;
  function Convert ( l : Lists_of_Float_Vectors.List )
                   return Lists_of_Integer_Vectors.List;

  function Convert ( l : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Float_Vector_Lists.Array_of_Lists;
  function Convert ( l : Arrays_of_Float_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Integer_Vector_Lists.Array_of_Lists;

  function Convert ( m : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return Float_Mixed_Subdivisions.Mixed_Cell;
  function Convert ( m : Float_Mixed_Subdivisions.Mixed_Cell )
                   return Integer_Mixed_Subdivisions.Mixed_Cell;

  function Convert ( s : Integer_Mixed_Subdivisions.Mixed_Subdivision )
                   return Float_Mixed_Subdivisions.Mixed_Subdivision;
  function Convert ( s : Float_Mixed_Subdivisions.Mixed_Subdivision )
                   return Integer_Mixed_Subdivisions.Mixed_Subdivision;

end Float_Integer_Convertors;
