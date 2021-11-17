with Integer_Vectors_of_Vectors;   use Integer_Vectors_of_Vectors;
with Transformations;              use Transformations;
with Complex_Vectors; 
with Solutions;                    use Solutions;

package Binomials is

-- DESCRIPTION :
--   This package contains a routine for the solution of a binomial system,
--   where the i-th equation is defined as follows:
--
--         X^u_i - c_i = 0, for i = 1,2,..,n,
--   with
--        [ u ]_i a vector of degrees and [ c ]_i a vector of complex numbers,
--
--   using a multi_index notation:  X^u_i = x1^u_i1 x2^u_i2 . . . xn^u_in.

  procedure Factorize ( v : in out Vector; n : in natural;
                        t : in out Transfo );

  -- DESCRIPTION :
  --   This routines factorizes the binomial system defined by the
  --   degrees in the vector v.

  -- ON ENTRY :
  --   v          defines the degrees of binomial system;
  --   n          the number of unknowns to be eliminated.

  -- ON RETURN :
  --   v          the factorized binomial system;
  --   t          the transformations used to factorize.

  procedure Solve ( v : in Vector; cv : in Complex_Vectors.Vector;
                    n : in natural; sols : in out Solution_List );

  -- DESCRIPTION :
  --   This routine solves the binomial system defined by the degrees
  --   in the vector v and by the constants in the vector cv.

  -- REQUIRED :
  --   The vector cv contains no zero components!

  -- ON ENTRY :
  --   v          defines the degrees of binomial system;
  --   cv         is the vector of constants;
  --   n          the dimension of the system.

  -- ON RETURN :
  --   sols       is the solution list of the binomial system.

  procedure Residuals ( v : in Vector; cv : in Complex_Vectors.Vector;
			n : in natural; sols : in Solution_List;
			res : out Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   This routine computes the residuals of the solutions of
  --   the binomial system defined by v and cv.

  -- REQUIRED :
  --   The dimension of res equals the number of solutions in sols.

  -- ON ENTRY :
  --   v          defines the degrees of the binomial system;
  --   cv         the vector of constants;
  --   n          the dimension of the system;
  --   sols       is the solution list of the binomial system.

  -- ON RETURN :
  --   res        a vector of residuals.

end Binomials;
