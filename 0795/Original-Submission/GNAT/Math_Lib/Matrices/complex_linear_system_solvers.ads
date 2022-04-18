with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Integer_Vectors,Complex_Vectors;
with Complex_Numbers,Complex_Matrices;  use Complex_Numbers,Complex_Matrices;

package Complex_Linear_System_Solvers is

-- DESCRIPTION :
--   This package offers a few routines to solve linear systems of equations.
--   The code for lufac, lufco and lusolve is a literal translation from the
--   f77-linpack code.

  procedure Scale ( a : in out Matrix; b : in out Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Divides the ith equation in the system a*x = b by the largest
  --   element on the ith row of a, for i in a'range(1).

  -- REQUIRED : a'range(1) = b'range(1).

  procedure lufac ( a : in out Matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector; info : out integer );

  -- DESCRIPTION :
  --   lufac factors a complex matrix by gaussian elimination

  --   lufac is usually called by lufco, but it can be called
  --   directly with a saving of time if rcond is not needed.
  --   (time for lufco) = (1 + 9/n)*(time for lufac).

  -- ON ENTRY :
  --   a       complex matrix(1..n,1..n) to be factored
  --   n       the dimension of the matrix a

  -- ON RETURN :
  --   a       an upper triangular matrix and the multipliers
  --           which were used to obtain it.
  --           The factorization can be written a = l*u where
  --           l is a product of permutation and unit lower
  --           triangular matrices and u is upper triangular.
  --   ipvt    an integer vector of pivot indices
  --   info    = 0  normal value
  --           = k  if u(k,k) = 0.0.
  --                This is not an error for this routine,
  --                but it does indicate that lusolve will
  --                divide by zero if called.  Use rcond in
  --                lufco for a reliable indication of singularity.

  procedure lufco ( a : in out Matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector;
                    rcond : out double_float );

  -- DESCRIPTION :
  --   lufco factors a complex matrix by gaussian elimination
  --   and estimates the condition of the matrix.

  -- If rcond is not needed, lufac is slightly faster.
  -- To solve a*x = b, follow lufco by lusolve.

  -- ON ENTRY :
  --   a       complex matrix(1..n,1..n) to be factored
  --   n       the dimension of the matrix a

  -- ON RETURN :
  --   a       an upper triangular matrix and the multipliers 
  --           which are used to obtain it.
  --           The factorization can be written a = l*u, where
  --           l is a product of permutation and unit lower triangular
  --           matrices and u is upper triangular.
  --   ipvt    an integer vector of pivot indices
  --   rcond   an estimate of the reciprocal condition of a.
  --           For the system a*x = b, relative perturbations
  --           in a and b of size epsilon may cause relative
  --           perturbations in x of size epsilon/rcond.
  --           If rcond is so small that the logical expression
  --                  1.0 + rcond = 1.0 
  --           is true, than a may be singular to working precision.
  --           In particular, rcond is zero if exact singularity is
  --           detected or the estimate underflows.

  procedure lusolve ( a : in Matrix; n : in integer;
                      ipvt : in Integer_Vectors.Vector;
                      b : in out Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   lusolve solves the complex system a*x = b using the factors
  --   computed by lufac or lufco

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..n), the output from
  --           lufac or lufco
  --   n       the dimension of the matrix a
  --   ipvt    the pivot vector from lufac or lufco
  --   b       the right hand side vector

  -- ON RETURN :
  --   b       the solution vector x

  procedure Triangulate ( a : in out Matrix; n,m : in integer );

  -- DESCRIPTION :
  --   triangulate makes the n*m complex matrix a triangular using
  --   Gaussian elimination.

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..m)
  --   n       the number of rows of a
  --   m       the number of columns of a

  -- ON RETURN :
  --   a       the triangulated matrix

  procedure Diagonalize ( a : in out Matrix; n,m : in integer );

  -- DESCRIPTION :
  --   diagonalize makes the n*m complex matrix a diagonal using
  --   Gauss-Jordan.

  -- ON ENTRY :
  --   a       a complex matrix(1..n,1..m)
  --   n       the number of rows of a
  --   m       the number of columns of a

  -- ON RETURN :
  --   a       the diagonalized matrix

end Complex_Linear_System_Solvers;
