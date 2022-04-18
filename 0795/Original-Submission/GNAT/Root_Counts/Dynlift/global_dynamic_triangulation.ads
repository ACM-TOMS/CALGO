with Integer_Vectors,Simplices;    use Integer_Vectors,Simplices;
with Triangulations;               use Triangulations;
with Lists_of_Integer_Vectors;     use Lists_of_Integer_Vectors;

package Global_Dynamic_Triangulation is

-- DESCRIPTION :
--   This package contains some global aspects of the dynamic lifting
--   algorithm applied to the unmixed case.

  procedure Initial_Simplex ( pts : in List; order : in boolean; 
                              s : out Simplex; rest : in out List );

  -- DESCRIPTION :
  --   Searches for an initial simplex of the triangulation of the
  --   list of points pts.

  -- ON ENTRY :
  --   pts       a list of integer vectors;
  --   order     if true, then the first n+1 points of the list will be
  --             considered first for the construction of the initial cell.

  -- ON RETURN :
  --   s         an initial simplex, if Volume(s) = 0, then Volume(pts) = 0,
  --             the lifting values of the initial cell are all zero;
  --   rest      the rest of the points: pts - s.

  function Max_Extreme ( l : List; k : natural ) return Link_to_Vector;

  -- DESCRIPTION :
  --   Returns the vector in the list l with the maximal kth component.

  function Max_Extreme ( l : List; weights : vector ) return Link_to_Vector;

  -- DESCRIPTION :
  --   Return the element in l for which the weighted sum is maximal.

  function Max_Extreme ( l : List; n : natural; low,upp : integer )
                       return Link_to_Vector;

  -- DESCRIPTION :
  --   Generates first a vector with n random weights, in between low and upp,
  --   and gives the maximal extreme w.r.t. this weighted sum.

end Global_Dynamic_Triangulation;
