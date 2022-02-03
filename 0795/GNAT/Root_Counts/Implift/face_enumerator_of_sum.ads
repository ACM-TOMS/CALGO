with Integer_Vectors;                     use Integer_Vectors;
with Integer_Vectors_of_Vectors;

package Face_Enumerator_of_Sum is

-- DESCRIPTION :
--   This package contains a procedure which allows the efficient
--   enumeration of all facets of a certain type of the sum of
--   a tuple of polytopes.

  generic

    with procedure Process ( faces : in Integer_Vectors_of_Vectors.Vector;
                             cont : out boolean );

    -- DESCRIPTION :
    --   The parameter faces has the following meaning:
    --   faces(i) contains the entries of points in pts which span
    --   the k-face of the ith polytope, with k = typ(i).
    --   When cont is set to false, the enumeration stops.

  procedure Enumerate_Faces_of_Sum
               ( ind,typ : in vector; k : in natural;
                 pts : in Integer_Vectors_of_Vectors.Vector );

  -- DESCRIPTION :
  --   Enumerates the points which span faces of a certain type
  --   of the polytope in lexicographic order.

  -- ON ENTRY :
  --   ind       ind(i) indicates where the points of the ith polytope
  --             in the vector pts begin;
  --   typ       typ(i) contains the dimension k of the k-face of the ith
  --             polytope which spans the facet of the sum;
  --   k         sum of the entries in typ, must be less than
  --             or equal to the dimension of the space;
  --   pts       contains all points of the polytopes.

  -- ON RETURN :
  --   Each time a face of the sum has been found, Process is invoked.

end Face_Enumerator_of_Sum;
