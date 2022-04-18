with Integer_Vectors;                 use Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Integer_Matrices;                use Integer_Matrices;
with Lists_of_Integer_Vectors;        use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;  use Arrays_of_Integer_Vector_Lists;

package Frequency_Graph is

-- DESCRIPTION :
--   This package allows to compute the frequency of a vector w.r.t.
--   a tuple of point sets.

-- CREATORS :

  function Occurrences ( i : natural; l : List ) return natural;

  -- DESCRIPTION :
  --   Returns the number of times the ith component of vectors in l
  --   differs from zero.

  function Graph ( n : natural; supports : Array_of_Lists ) return matrix;

  -- DESCRIPTION :
  --   Computes the graph wich represents the occurencies of each variable
  --   in each support.  In particular, if Graph(i,j) = 0, then the ith
  --   unknown does not occur in the jth support, else the number of
  --   occurrences in different vectors is given.

-- MODIFIER :

  procedure Ignore ( m : in out matrix; point : in vector );

  -- DESCRIPTION :
  --   For each component i of the point, for which point(i) /= 0,
  --   the element m(i,j) := 1, for all j in m'range(2).
  --   This is to ignore the occurrences of the ith unknown and to give
  --   priority to those unknowns that have not been chosen yet.

-- SELECTORS :

  function Occurrence ( i : natural; m : matrix ) return natural;
  function Occurrence ( i : natural; m : matrix; col : natural; perm : vector )
                     return natural;

  -- DESCRIPTION :
  --   Returns the number of nonzero columns of ith row of m.
  --   The second function does not consider a number of columns
  --   equal to col and applies Occurrence to the columns
  --   perm(j) for j in col+1..m'last(2).

  function Occurrence ( v : vector; m : matrix ) return natural;
  function Occurrence ( v : vector; m : matrix; col : natural; perm : vector )
                      return natural;

  -- DESCRIPTION :
  --   Returns the minimal Occurrence(i,m), for all i, for which v(i) /= 0.
  --   The second function does not consider a number of columns
  --   equal to col and applies Occurrence to the columns
  --   perm(j) for j in col+1..m'last(2).

  function Lowest_Occurrence
               ( vec : Integer_Vectors_of_Vectors.Vector;
                 start : natural; m : matrix ) return natural;
  function Lowest_Occurrence
               ( vec : Integer_Vectors_of_Vectors.Vector;
                 start : natural; m : matrix;
                 col : natural; perm : vector  ) return natural;

  -- DESCRIPTION :
  --   Returns the index of the vector with lowest occurrence in
  --   vec(start..vec'last).
  --   The second function does not consider a number of columns
  --   equal to col and applies Lowest_Occurrence to the columns
  --   perm(j) for j in col+1..m'last(2).

-- CONSTRUCTORS :

  function  Sort ( l : List; m : matrix ) return List;
  procedure Sort ( l : in out List; m : in matrix );

  -- DESCRIPTION :
  --   Sorts the points in the list l, so that on return,
  --   the points are ordered, from low to high occurrence.

  function  Sort ( l : List; m : matrix;
                   col : natural; perm : vector ) return List;
  procedure Sort ( l : in out List; m : in matrix;
                   col : in natural; perm : in vector );

  -- DESCRIPTION :
  --   Sorts the points in the list l, so that on return,
  --   the points are ordered, from low to high occurrence.
  --   Only the columns perm(j), for j in col+1..m'last(2) will be
  --   considered when computing the occurrences.

end Frequency_Graph;
