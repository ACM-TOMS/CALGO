with Integer_Vectors,Integer_Matrices;  use Integer_Vectors,Integer_Matrices;

package Integer_Linear_Inequality_Solvers is

-- DESCRIPTION :
--   This package provides for incrementally solving systems of
--   linear systems w.r.t. linear inequalities.

  procedure Triangulate ( l : in natural; m : in matrix;
                          first,last : in natural; ineq : in out matrix );
  procedure Triangulate ( l : in natural; m : in matrix;
                          index : in natural; ineq : in out matrix );

  -- DESCRIPTION :
  --   Updates a matrix of inequalities, after elimination of the lth unknown.

  -- ON ENTRY :
  --   l         current unknown to be eliminated;
  --   m         m(1..l,m'range(2)) is upper triangular;
  --   first     indicates start in range of ineq to be updated;
  --   last      indicates end in range of ineq to be updated;
  --   index     indicates the inequality to be updated;
  --   ineq      in ineq(first..last) or in ineq(index),
  --             the first l-1 unknowns are already eliminated.

  -- ON RETURN :
  --   ineq      the updated inequalities.

  procedure Triangulate ( index,start : in natural; ineq : in out matrix );

  -- DESCRIPTION :
  --   Updates the inequality ineq(index), i.e. tries to eliminate as
  --   many unknowns in ineq(index) as possible, by making positive
  --   combinations with inequalities in ineq(ineq'first..index-1).

  -- ON ENTRY :
  --   index     current row in ineq that has to be updated;
  --   start     indicates in which column the first nonzero elements
  --             have to be found;
  --   ineq      matrix of inequalities of type <.,.> >= 0.

  -- ON RETURN :
  --   ineq      the updated matrix of inequalities.

end Integer_Linear_Inequality_Solvers;
