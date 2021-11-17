with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Integer_Vectors;
with Float_Vectors,Float_Matrices;     use Float_Vectors,Float_Matrices;

package Farkas_Lemma is

-- DESCRIPTION :
--   This procedure contains some routines for verifying the Farkas lemma.

  procedure Complementary_Slackness
                 ( tableau : in out matrix; lastcol : in integer;
                   rhs : in out vector; tol : in double_float;
                   solution : out vector; columns : out Integer_Vectors.Vector;
                   feasible : out boolean );

  procedure Complementary_Slackness
                 ( tableau : in out matrix; rhs : in out vector;
                   tol : in double_float; solution : out vector;
                   columns : out Integer_Vectors.Vector;
                   feasible : out boolean );

  -- DESCRIPTION :
  --   Solves the complementary slackness problem: determines
  --   whether there exists a positive combination of the columns
  --   such that the right hand side is satisfied.

  -- REQUIRED :
  --   rhs'range = solution'range = columns'range = tableau'range(1)

  -- ON ENTRY :
  --   tableau     inequalities as columns;
  --   lastcol     indicates the last significant column in the tableau,
  --                if not given, then lastcol = tableau'last(2);
  --   tol         tolerance to decide whether a number equals zero.
  --   rhs         right hand side vector;

  -- ON RETURN :
  --   tableau     modified tableau of inequalities;
  --   rhs         modified right hand side;
  --   solution    the computed solution;
  --   columns     indicates which columns has been used;
  --   feasible    if true then the solution is feasible.

end Farkas_Lemma;
